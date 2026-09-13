#!/usr/bin/env python
"""Runs the hicfilecpp equivalence harness: every case through hicstraw
(oracle.py, in its own environment) and through hicfilecpp (the
hicfilecpp-harness driver), compares the outputs, measures peak RSS and CPU
time of both processes, applies the gates and writes a report.

    python harness/run.py --driver BUILD/harness/hicfilecpp-harness \\
        --oracle-python ENV/bin/python \\
        --hicx-data ~/src/HiCExplorer-v4/hicexplorer/test/test_data \\
        --out REPORT_DIR [--filter REGEX] [--cases FILE ...]

Path tokens in case files: @hicx/ (HiCExplorer test data) and @data/ (this
repository's tests/data). Every file they name must be listed with its SHA-256
in harness/data_manifest.json.

A case file holds explicit "cases" and "generate" entries. A generate entry
names a file; the runner expands it, with harness/hicindex.py, into a header
case, one records case per unit, resolution and matrix type covering every
chromosome pair and every normalization the file holds, one vectors case per
unit and resolution, and a records case for the whole-genome matrix.

Gates: hicfilecpp may not use more CPU time than hicstraw, and its peak RSS
must stay within hicstraw's.
"""

import argparse
import glob
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import compare  # noqa: E402
import hicindex  # noqa: E402

PROCESS_TIMEOUT_S = 3600
MEASURE = None


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def run_measured(cmd, log_path):
    report = log_path + ".rusage.json"
    if os.path.exists(report):
        os.remove(report)
    with open(log_path, "wb") as log:
        proc = subprocess.Popen(
            [MEASURE, report, str(PROCESS_TIMEOUT_S), "--"] + cmd,
            stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
        )
        try:
            proc.wait()
        except BaseException:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            raise
    if proc.returncode != 0 or not os.path.exists(report):
        return {"exit": -1, "cpu_s": float("nan"), "rss_mb": float("nan"), "wall_s": float("nan")}
    with open(report) as handle:
        return json.load(handle)


class Resolver:
    def __init__(self, hicx, manifest):
        self.roots = {"@hicx/": hicx, "@data/": os.path.join(REPO, "tests", "data")}
        self.manifest = manifest

    def path(self, token):
        for prefix, root in self.roots.items():
            if token.startswith(prefix):
                return os.path.join(root, token[len(prefix):])
        return token

    def substitute(self, value):
        if isinstance(value, dict):
            return {k: self.substitute(v) for k, v in value.items()}
        if isinstance(value, list):
            return [self.substitute(v) for v in value]
        if isinstance(value, str) and value.startswith(("@hicx/", "@data/")):
            return self.path(value)
        return value

    def verify(self, token):
        entry = self.manifest.get(token)
        path = self.path(token)
        if entry is None:
            return f"{token} is not listed in data_manifest.json"
        if not os.path.exists(path):
            return f"{token} is missing ({path})"
        if os.path.getsize(path) != entry["size"] or sha256(path) != entry["sha256"]:
            return f"{token} does not match its manifest checksum"
        return None


def expand(entry, resolver):
    """Cases for one generate entry."""
    token = entry["file"]
    prefix = entry["prefix"]
    index = hicindex.read_index(resolver.path(token))
    chromosomes = [c for c in index["chromosomes"] if c[1] > 0]
    norms_at = {}
    for norm, chrom, unit, resolution in index["norm_vectors"]:
        norms_at.setdefault((unit, resolution), {}).setdefault(norm, set()).add(chrom)
    cases = [{"name": f"H.{prefix}", "op": "header", "file": token}]

    def extent(chrom, unit, resolution):
        if unit == "BP":
            return chrom[2]
        return index["frag_sites"][chrom[1]] + 1 + resolution

    units = [("BP", r) for r in index["bp"]] + [("FRAG", r) for r in index["frag"]]
    for unit, resolution in units:
        norms = ["NONE"] + sorted(norms_at.get((unit, resolution), {}))
        for matrix_type in ("observed", "oe"):
            queries = []
            for norm in norms:
                with_vectors = norms_at.get((unit, resolution), {}).get(norm, set())
                for i, a in enumerate(chromosomes):
                    for b in chromosomes[i:]:
                        key = f"{a[1]}_{b[1]}"
                        if key in index["matrices"] and (unit, resolution) not in index["zooms"][key]:
                            continue
                        if norm != "NONE" and not (a[1] in with_vectors and b[1] in with_vectors):
                            continue
                        queries.append([norm, a[0], b[0], 0, extent(a, unit, resolution),
                                        0, extent(b, unit, resolution)])
            cases.append({"name": f"R.{prefix}.{unit}.{resolution}.{matrix_type}", "op": "records",
                          "file": token, "matrix_type": matrix_type, "unit": unit,
                          "resolution": resolution, "queries": queries})
        vector_queries = []
        for norm in norms:
            with_vectors = norms_at.get((unit, resolution), {}).get(norm, set())
            for chrom in chromosomes:
                key = f"{chrom[1]}_{chrom[1]}"
                if key not in index["matrices"] or (unit, resolution) not in index["zooms"][key]:
                    continue
                if norm != "NONE" and chrom[1] not in with_vectors:
                    continue
                vector_queries.append([norm, chrom[0], chrom[1]])
        cases.append({"name": f"V.{prefix}.{unit}.{resolution}", "op": "vectors", "file": token,
                      "matrix_type": "oe", "unit": unit, "resolution": resolution,
                      "queries": vector_queries})
    whole = index["chromosomes"][0]
    for unit, bin_size in index["zooms"].get("0_0", []):
        cases.append({"name": f"R.{prefix}.All.{unit}.{bin_size}", "op": "records", "file": token,
                      "matrix_type": "observed", "unit": unit, "resolution": bin_size,
                      "queries": [["NONE", whole[0], whole[0], 0, whole[2], 0, whole[2]]]})
    return cases


def run_case(case, args, resolver, out_root):
    name = case["name"]
    work = os.path.join(out_root, "cases", re.sub(r"[^A-Za-z0-9_.-]", "_", name))
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    record = {"name": name, "op": case["op"], "problems": [], "notes": []}
    problem = resolver.verify(case["file"])
    if problem:
        record["problems"].append(problem)
        record["verdict"] = "FAIL"
        return record
    resolved = resolver.substitute(case)
    sides = {}
    for side in ("py", "cpp"):
        side_dir = os.path.join(work, side)
        out_dir = os.path.join(side_dir, "out")
        os.makedirs(out_dir)
        case_path = os.path.join(side_dir, "case.json")
        with open(case_path, "w") as handle:
            json.dump(resolved, handle)
        if side == "py":
            cmd = [args.oracle_python, os.path.join(HERE, "oracle.py"), "run", case_path, out_dir]
        else:
            cmd = [args.driver, case_path, out_dir]
        measure = run_measured(cmd, os.path.join(side_dir, "log.txt"))
        sides[side] = {"dir": out_dir, "measure": measure}
        if measure["exit"] != 0 or not os.path.exists(os.path.join(out_dir, "result.json")):
            record["problems"].append(f"{side} process failed with exit status {measure['exit']}")
    record["py"] = sides["py"]["measure"]
    record["cpp"] = sides["cpp"]["measure"]
    if record["problems"]:
        record["verdict"] = "FAIL"
        return record
    with open(os.path.join(sides["py"]["dir"], "result.json")) as handle:
        py_doc = json.load(handle)
    with open(os.path.join(sides["cpp"]["dir"], "result.json")) as handle:
        cpp_doc = json.load(handle)
    record["py"]["op_s"] = py_doc["op_seconds"]
    record["cpp"]["op_s"] = cpp_doc["op_seconds"]
    cmp = compare.compare_results(py_doc["result"], cpp_doc["result"], sides["py"]["dir"], sides["cpp"]["dir"])
    if py_doc["result"]["kind"] == "error":
        record["notes"].append("both report an error")
    if case["op"] == "records":
        record["notes"].append(f"{len(case['queries'])} queries, {sum(py_doc['result'].get('n', []))} records")
    record["problems"].extend(cmp.problems)
    record["float_class"] = cmp.float_class
    record["cpu_gate"] = record["cpp"]["cpu_s"] <= record["py"]["cpu_s"]
    record["rss_gate"] = record["cpp"]["rss_mb"] <= record["py"]["rss_mb"]
    if not record["cpu_gate"]:
        record["problems"].append(
            f"CPU gate: hicfilecpp {record['cpp']['cpu_s']:.3f} s > hicstraw {record['py']['cpu_s']:.3f} s")
    if not record["rss_gate"]:
        record["problems"].append(
            f"RSS gate: hicfilecpp {record['cpp']['rss_mb']:.1f} MB > hicstraw {record['py']['rss_mb']:.1f} MB")
    record["verdict"] = "PASS" if not record["problems"] else "FAIL"
    return record


def write_report(records, out_root):
    passed = sum(1 for r in records if r["verdict"].lower() == "pass")
    failed = sum(1 for r in records if r["verdict"].lower() == "fail")
    lines = [
        "# hicfilecpp equivalence report", "",
        f"Date: {time.strftime('%Y-%m-%d')}",
        f"Cases: {len(records)}, pass: {passed}, fail: {failed}", "",
        "| case | verdict | floats | notes | hicstraw CPU s | C++ CPU s | hicstraw RSS MB | C++ RSS MB |",
        "|---|---|---|---|---:|---:|---:|---:|",
    ]
    for r in records:
        py, cpp = r.get("py", {}), r.get("cpp", {})
        lines.append("| {} | {} | {} | {} | {:.3f} | {:.3f} | {:.1f} | {:.1f} |".format(
            r["name"], r["verdict"], r.get("float_class") or "-", "; ".join(r["notes"]) or "-",
            py.get("cpu_s", float("nan")), cpp.get("cpu_s", float("nan")),
            py.get("rss_mb", float("nan")), cpp.get("rss_mb", float("nan"))))
    problems = [r for r in records if r["problems"]]
    if problems:
        lines += ["", "## Problems", ""]
        for r in problems:
            for p in r["problems"][:20]:
                lines.append(f"- {r['name']}: {p}")
    with open(os.path.join(out_root, "report.md"), "w") as handle:
        handle.write("\n".join(lines) + "\n")
    with open(os.path.join(out_root, "report.json"), "w") as handle:
        json.dump(records, handle, indent=1)
    return passed, failed


def load_cases(paths, resolver):
    cases = []
    for path in paths:
        with open(path) as handle:
            document = json.load(handle)
        for entry in document.get("generate", []):
            cases.extend(expand(entry, resolver))
        cases.extend(document.get("cases", []))
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--driver", required=True)
    parser.add_argument("--oracle-python", required=True, help="python of the environment with hicstraw")
    parser.add_argument("--measure", default=None)
    parser.add_argument("--hicx-data", required=True)
    parser.add_argument("--out", default=None)
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--filter", default=None)
    parser.add_argument("--update-manifest", action="store_true")
    args = parser.parse_args()
    global MEASURE
    MEASURE = args.measure or os.path.join(os.path.dirname(os.path.abspath(args.driver)), "hicfilecpp-measure")

    manifest_path = os.path.join(HERE, "data_manifest.json")
    manifest = {}
    if os.path.exists(manifest_path):
        with open(manifest_path) as handle:
            manifest = json.load(handle)
    resolver = Resolver(os.path.abspath(os.path.expanduser(args.hicx_data)), manifest)
    case_files = args.cases or sorted(glob.glob(os.path.join(HERE, "cases", "*.json")))

    if args.update_manifest:
        tokens = set()
        for path in case_files:
            with open(path) as handle:
                document = json.load(handle)
            tokens.update(e["file"] for e in document.get("generate", []))
            tokens.update(c["file"] for c in document.get("cases", []))
        for token in sorted(tokens):
            path = resolver.path(token)
            manifest[token] = {"sha256": sha256(path), "size": os.path.getsize(path)}
        with open(manifest_path, "w") as handle:
            json.dump(dict(sorted(manifest.items())), handle, indent=1)
            handle.write("\n")
        print(f"manifest: {len(manifest)} files")
        return 0

    if not args.out:
        sys.exit("--out is required")
    if not os.access(MEASURE, os.X_OK):
        sys.exit(f"hicfilecpp-measure not found at {MEASURE}")
    cases = load_cases(case_files, resolver)
    names = [c["name"] for c in cases]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        sys.exit(f"duplicate case names: {duplicates}")
    if args.filter:
        cases = [c for c in cases if re.search(args.filter, c["name"])]
    os.makedirs(args.out, exist_ok=True)
    records = []
    for case in cases:
        record = run_case(case, args, resolver, args.out)
        records.append(record)
        py, cpp = record.get("py", {}), record.get("cpp", {})
        print(f"{record['verdict']} {record['name']} floats={record.get('float_class')} "
              f"cpu py={py.get('cpu_s', 0):.3f}s cpp={cpp.get('cpu_s', 0):.3f}s "
              f"rss py={py.get('rss_mb', 0):.0f}MB cpp={cpp.get('rss_mb', 0):.0f}MB", flush=True)
        for problem in record["problems"][:5]:
            print(f"    {problem}", flush=True)
        write_report(records, args.out)
    if not records:
        write_report(records, args.out)
        print("no cases selected")
        return 2
    passed, failed = write_report(records, args.out)
    print(f"cases {len(records)}: pass {passed}, fail {failed}; report {os.path.join(args.out, 'report.md')}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))
    sys.exit(main())
