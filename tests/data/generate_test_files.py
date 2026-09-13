#!/usr/bin/env python
"""Regenerates the Juicer-written test files in this directory from the
5 kb contact records of HiCExplorer's SRR1791297_30.hic (sacCer3).

    python tests/data/generate_test_files.py --srr SRR1791297_30.hic \\
        --java JAVA --juicer8 juicer_tools_1.22.01.jar \\
        --juicer9 juicer_tools.2.20.00.jar --work DIR

Needs hicstraw (to read the source) and a Java 8 or newer runtime. Every pixel
becomes one contact of Juicer's "short with score" format at the start of its
bins, so Juicer bins it back into the same pixel at 5 kb and into the right
pixel at every resolution that is a multiple of 5 kb. Fragment numbers come
from a synthetic restriction site every 2 kb (sites.txt); the fragment maps
exist only to give the reader FRAG resolutions to read.

Files written into this directory:
  SRR1791297_30.juicer_tools_1.22.01.v8.hic       version 8, BP resolutions
  SRR1791297_30.juicer_tools_1.22.01.frag.v8.hic  version 8, BP and FRAG
  SRR1791297_30.juicer_tools_2.20.00.v9.hic       version 9, BP resolutions
"""

import argparse
import os
import shutil
import subprocess

import hicstraw

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE_RESOLUTION = 5000
SITE_SPACING = 2000


def write_inputs(srr, work):
    handle = hicstraw.HiCFile(srr)
    chromosomes = [c for c in handle.getChromosomes() if c.name.lower() != "all"]
    sizes = os.path.join(work, "sacCer3.chrom.sizes")
    with open(sizes, "w") as out:
        for c in chromosomes:
            out.write(f"{c.name}\t{c.length}\n")
    with open(os.path.join(work, "sites.txt"), "w") as out:
        for c in chromosomes:
            out.write(" ".join([c.name] + [str(p) for p in range(SITE_SPACING, c.length, SITE_SPACING)]) + "\n")
    contacts = os.path.join(work, "contacts.txt")
    with open(contacts, "w") as out:
        for i, a in enumerate(chromosomes):
            for b in chromosomes[i:]:
                mzd = handle.getMatrixZoomData(a.name, b.name, "observed", "NONE", "BP", SOURCE_RESOLUTION)
                for r in sorted(mzd.getRecords(0, a.length, 0, b.length), key=lambda r: (r.binX, r.binY)):
                    out.write(f"0 {a.name} {r.binX} {r.binX // SITE_SPACING} "
                              f"0 {b.name} {r.binY} {r.binY // SITE_SPACING} {r.counts:.9g}\n")
    return "sacCer3.chrom.sizes", contacts


def run(cmd, work, log):
    with open(os.path.join(work, log), "w") as out:
        subprocess.run(cmd, cwd=work, stdout=out, stderr=subprocess.STDOUT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--srr", required=True)
    parser.add_argument("--java", required=True)
    parser.add_argument("--juicer8", required=True)
    parser.add_argument("--juicer9", required=True)
    parser.add_argument("--work", required=True)
    args = parser.parse_args()
    os.makedirs(args.work, exist_ok=True)
    sizes, contacts = write_inputs(args.srr, args.work)
    java = [args.java, "-Xmx4g", "-jar"]
    outputs = {
        "SRR1791297_30.juicer_tools_1.22.01.v8.hic":
            java + [args.juicer8, "pre", "-j", "1", "-r", "1000000,250000,50000,10000", contacts, "v8.hic", sizes],
        "SRR1791297_30.juicer_tools_1.22.01.frag.v8.hic":
            java + [args.juicer8, "pre", "-j", "1", "-f", "sites.txt", "-r", "500000,50000,100f,20f",
                    contacts, "frag.v8.hic", sizes],
        "SRR1791297_30.juicer_tools_2.20.00.v9.hic":
            java + [args.juicer9, "pre", "-j", "1", "-k", "VC,VC_SQRT,KR,SCALE,INTER_SCALE,GW_SCALE",
                    "-r", "1000000,250000,50000,10000", contacts, "v9.hic", sizes],
    }
    for name, cmd in outputs.items():
        run(cmd, args.work, name + ".log")
        shutil.copyfile(os.path.join(args.work, cmd[-2]), os.path.join(HERE, name))
        print(name, os.path.getsize(os.path.join(HERE, name)))


if __name__ == "__main__":
    main()
