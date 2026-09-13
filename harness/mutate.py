#!/usr/bin/env python
"""Proves that harness cases can fail: applies each mutation of
harness/mutations.json to a scratch copy of the hicstraw source, builds that
copy into a scratch directory, runs the cases the mutation targets with the
mutated hicstraw first on PYTHONPATH, and requires at least one of them to
fail. The installed hicstraw is never touched.

    python harness/mutate.py --driver BUILD/harness/hicfilecpp-harness \\
        --oracle-python ENV/bin/python --hicx-data HICX_TEST_DATA \\
        --straw-source STRAW/pybind11_python --scratch DIR \\
        [--include-dir DIR] [--library-dir DIR]

--include-dir and --library-dir locate libcurl, which hicstraw's source
includes; pass them when curl's headers are not on the compiler's path.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--driver", required=True)
    parser.add_argument("--oracle-python", required=True)
    parser.add_argument("--hicx-data", required=True)
    parser.add_argument("--straw-source", required=True)
    parser.add_argument("--scratch", required=True)
    parser.add_argument("--include-dir", default=None)
    parser.add_argument("--library-dir", default=None)
    parser.add_argument("--only", default=None, help="run only the mutations whose name matches this regex")
    args = parser.parse_args()

    with open(os.path.join(HERE, "mutations.json")) as handle:
        mutations = json.load(handle)["mutations"]
    if args.only:
        mutations = [m for m in mutations if re.search(args.only, m["name"])]
    env_build = dict(os.environ)
    if args.include_dir:
        env_build["CFLAGS"] = f"-I{args.include_dir}"
        env_build["CPPFLAGS"] = f"-I{args.include_dir}"
    if args.library_dir:
        env_build["LDFLAGS"] = f"-L{args.library_dir} -Wl,-rpath,{args.library_dir}"

    results = []
    for mutation in mutations:
        root = os.path.join(args.scratch, mutation["name"])
        shutil.rmtree(root, ignore_errors=True)
        source = os.path.join(root, "source")
        shutil.copytree(args.straw_source, source, ignore=shutil.ignore_patterns("build", "*.egg-info", "__pycache__"))
        target = os.path.join(source, "src", "straw.cpp")
        with open(target) as handle:
            text = handle.read()
        if text.count(mutation["old"]) != 1:
            results.append((mutation["name"], "INVALID (pattern not found exactly once)"))
            print(f"{mutation['name']}: {results[-1][1]}", flush=True)
            continue
        with open(target, "w") as handle:
            handle.write(text.replace(mutation["old"], mutation["new"]))
        site = os.path.join(root, "site")
        build = subprocess.run(
            [args.oracle_python, "-m", "pip", "install", "--quiet", "--no-build-isolation", "--no-deps",
             "--target", site, source],
            env=env_build, capture_output=True, text=True,
        )
        if build.returncode != 0:
            results.append((mutation["name"], "INVALID (mutated hicstraw does not build)"))
            print(f"{mutation['name']}: {results[-1][1]}\n{build.stderr[-2000:]}", flush=True)
            continue
        env = dict(os.environ)
        env["PYTHONPATH"] = site + os.pathsep + env.get("PYTHONPATH", "")
        proc = subprocess.run(
            [sys.executable, os.path.join(HERE, "run.py"), "--driver", args.driver,
             "--oracle-python", args.oracle_python, "--hicx-data", args.hicx_data,
             "--out", os.path.join(root, "report"), "--filter", mutation["cases"]],
            env=env, capture_output=True, text=True,
        )
        lines = proc.stdout.splitlines()
        verdicts = [line.split()[0] for line in lines if re.match(r"^(pass|fail)\b", line, re.I)]
        failed = sum(1 for v in verdicts if v.lower() == "fail")
        if not verdicts:
            status = "INVALID (no case ran)"
        elif failed > 0 and proc.returncode == 1:
            status = f"caught ({failed} of {len(verdicts)} cases fail)"
        else:
            status = f"MISSED ({len(verdicts)} cases pass)"
        results.append((mutation["name"], status))
        print(f"{mutation['name']}: {status}", flush=True)

    missed = [name for name, status in results if not status.startswith("caught")]
    print(f"mutations {len(results)}: caught {len(results) - len(missed)}, not caught {len(missed)}")
    return 0 if not missed else 1


if __name__ == "__main__":
    sys.exit(main())
