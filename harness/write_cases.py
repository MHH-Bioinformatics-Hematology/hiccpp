"""Writer cases of the hiccpp harness, run from run.py.

A write case names a source (.hic or fixed-bin .cool), its resolution, a
.hic version, resolutions and normalizations. The runner

1. turns the source into the same contacts for both writers (writer_io.py prep);
2. runs Juicer tools pre on them (1.22.01 for version 8, 2.20.00 for
   version 9) with the case's resolutions and normalizations: the reference;
3. runs hiccpp-harness write on them (with via_addnorm, writing without
   normalizations and then running addNorm), measured against the reference:
   hiccpp may not use more CPU time;
4. writes the file twice more, on 1 thread and again on the case's threads,
   and requires identical bytes;
5. reads both files back through hicstraw (writer_io.py check): observed
   records at every resolution against Juicer's and, at the source resolution,
   against the source; normalization vectors and expected values to the ED gate;
6. reads the hiccpp file back through hic2cool 1.0.1 (version 8) and
   through Juicer tools dump, and requires the source pixels exactly.
"""

import json
import os
import re
import shutil
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(path):
    with open(path) as handle:
        return json.load(handle)


def _dump_to_rows(text_path, chr1, chr2, resolution):
    rows = []
    with open(text_path) as handle:
        for line in handle:
            parts = line.split()
            if len(parts) != 3:
                continue
            x, y, value = int(parts[0]), int(parts[1]), float(parts[2])
            if chr1 == chr2 and x > y:
                x, y = y, x
            rows.append([chr1, chr2, x // resolution, y // resolution, np.float32(value)])
    return rows


def run_write_case(case, args, resolver, out_root, run_measured):
    name = case["name"]
    work = os.path.join(out_root, "cases", re.sub(r"[^A-Za-z0-9_.-]", "_", name))
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    record = {"name": name, "op": "write", "problems": [], "notes": []}

    def fail(text):
        record["problems"].append(text)
        record["verdict"] = "FAIL"
        return record

    problem = resolver.verify(case["file"])
    if problem:
        return fail(problem)
    version = case["version"]
    jar = args.juicer8 if version == 8 else args.juicer9
    if not args.java or not jar:
        return fail(f"write cases need --java and --juicer{version}")
    source = resolver.path(case["file"])
    resolution = case["source_resolution"]
    inputs = os.path.join(work, "inputs")
    os.makedirs(inputs)
    prep = run_measured([args.oracle_python, os.path.join(HERE, "writer_io.py"), "prep", source, str(resolution),
                         inputs], os.path.join(work, "prep.log"))
    if prep["exit"] != 0:
        return fail("input preparation failed")
    source_info = _load(os.path.join(inputs, "source.json"))
    record["notes"].append(f"{source_info['pixels']} source pixels, "
                           f"{'integer' if source_info['integer_counts'] else 'float'} counts")

    resolutions = sorted(case["resolutions"], reverse=True)
    norms = case["normalizations"]
    juicer_dir = os.path.join(work, "juicer")
    os.makedirs(juicer_dir)
    reference = os.path.join(juicer_dir, "out.hic")
    cmd = [args.java, f"-Xmx{args.java_heap}", "-jar", jar, "pre", "-j", "1",
           "-r", ",".join(str(r) for r in resolutions)]
    cmd += ["-k", ",".join(norms)] if norms else ["-n"]
    cmd += [os.path.join(inputs, "contacts.txt"), reference, os.path.join(inputs, "genome.chrom.sizes")]
    record["py"] = run_measured(cmd, os.path.join(juicer_dir, "log.txt"))
    if record["py"]["exit"] != 0 or not os.path.exists(reference):
        return fail(f"Juicer tools pre failed with exit status {record['py']['exit']}")

    def write(threads, output, label):
        side = os.path.join(work, label)
        os.makedirs(side, exist_ok=True)
        spec = {"op": "write", "inputs": inputs, "output": output, "version": version, "resolutions": resolutions,
                "source_resolution": resolution, "normalizations": norms, "threads": threads,
                "genome": "genome.chrom.sizes", "via_addnorm": bool(case.get("via_addnorm", False))}
        case_path = os.path.join(side, "case.json")
        with open(case_path, "w") as handle:
            json.dump(spec, handle)
        return run_measured([args.driver, case_path, side], os.path.join(side, "log.txt"))

    threads = case.get("threads", 8)
    candidate = os.path.join(work, "cpp", "out.hic")
    record["cpp"] = write(threads, candidate, "cpp")
    if record["cpp"]["exit"] != 0 or not os.path.exists(candidate):
        return fail(f"hiccpp write failed with exit status {record['cpp']['exit']}")
    record["py"]["op_s"] = record["py"]["wall_s"]
    record["cpp"]["op_s"] = _load(os.path.join(work, "cpp", "result.json"))["op_seconds"]

    one = os.path.join(work, "cpp_t1", "out.hic")
    again = os.path.join(work, "cpp_again", "out.hic")
    write(1, one, "cpp_t1")
    write(threads, again, "cpp_again")
    with open(candidate, "rb") as a, open(one, "rb") as b, open(again, "rb") as c:
        first = a.read()
        if first != b.read():
            record["problems"].append(f"determinism: output on 1 thread differs from output on {threads}")
        if first != c.read():
            record["problems"].append("determinism: a repeated run gives different bytes")
    record["notes"].append(f"{len(first)} bytes, identical on 1 and {threads} threads and when repeated"
                           if not any(p.startswith("determinism") for p in record["problems"]) else "not deterministic")

    check_spec = os.path.join(work, "check.json")
    with open(check_spec, "w") as handle:
        json.dump({"inputs": inputs, "candidate": candidate, "reference": reference,
                   "source_resolution": resolution, "normalizations": norms}, handle)
    checked = os.path.join(work, "check.result.json")
    run = run_measured([args.oracle_python, os.path.join(HERE, "writer_io.py"), "check", check_spec, checked],
                       os.path.join(work, "check.log"))
    if run["exit"] != 0 or not os.path.exists(checked):
        return fail("hicstraw read-back failed")
    result = _load(checked)
    record["problems"].extend(f"hicstraw: {p}" for p in result["problems"][:20])
    record["float_class"] = result["float_class"]
    counts = result["counts"]
    record["notes"].append(f"hicstraw: {counts['records']} records in {counts['record_sets']} pair matrices, "
                           f"{counts['vectors']} vectors, {counts['expected']} expected vectors")
    record["worst_by_kind"] = result["worst_by_kind"]

    chroms = _load(os.path.join(inputs, "chroms.json"))
    readbacks = []
    if version == 8 and source_info["integer_counts"] and args.hic2cool_python:
        rows = os.path.join(work, "hic2cool.npy")
        run = run_measured([args.hic2cool_python, os.path.join(HERE, "hic2cool_readback.py"), candidate,
                            str(resolution), os.path.join(inputs, "chroms.json"), rows,
                            os.path.join(work, "hic2cool.cool")], os.path.join(work, "hic2cool.log"))
        readbacks.append(("hic2cool", run, rows, "all"))
    elif version == 8 and args.hic2cool_python:
        record["notes"].append("hic2cool skipped: it stores counts as int32")
    elif version > 8:
        record["notes"].append("hic2cool skipped: it reads version 8 only")
    dump_pairs = case.get("dump_pairs", [[0, 0]])
    rows = []
    dump_ok = {"exit": 0}
    for chr1, chr2 in dump_pairs:
        text = os.path.join(work, f"dump_{chr1}_{chr2}.txt")
        run = run_measured([args.java, f"-Xmx{args.java_heap}", "-jar", jar, "dump", "observed", "NONE", candidate,
                            chroms[chr1][0], chroms[chr2][0], "BP", str(resolution), text],
                           os.path.join(work, f"dump_{chr1}_{chr2}.log"))
        if run["exit"] != 0 or not os.path.exists(text):
            dump_ok = run
            break
        rows.extend(_dump_to_rows(text, chr1, chr2, resolution))
    dump_rows = os.path.join(work, "dump.npy")
    np.save(dump_rows, np.asarray(rows, dtype=np.float64).reshape(-1, 5))
    readbacks.append(("Juicer tools dump", dump_ok, dump_rows, ",".join(f"{a}_{b}" for a, b in dump_pairs)))
    for label, run, rows_path, pairs in readbacks:
        if run["exit"] != 0 or not os.path.exists(rows_path):
            record["problems"].append(f"{label} read-back failed")
            continue
        out = rows_path + ".result.json"
        compared = run_measured([args.oracle_python, os.path.join(HERE, "writer_io.py"), "compare-pixels", inputs,
                                 rows_path, pairs, out], rows_path + ".compare.log")
        if compared["exit"] != 0:
            record["problems"].append(f"{label} comparison failed")
            continue
        result = _load(out)
        record["problems"].extend(f"{label}: {p}" for p in result["problems"])
        record["notes"].append(f"{label}: {result['pixels']} pixels of {result['pairs']} pairs")

    record["cpu_gate"] = record["cpp"]["cpu_s"] <= record["py"]["cpu_s"]
    if not record["cpu_gate"]:
        record["problems"].append(
            f"CPU gate: hiccpp {record['cpp']['cpu_s']:.3f} s > Juicer tools {record['py']['cpu_s']:.3f} s")
    record["rss_gate"] = True
    record["verdict"] = "PASS" if not record["problems"] else "FAIL"
    return record
