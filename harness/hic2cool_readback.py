#!/usr/bin/env python
"""Reads a version 8 .hic file back through hic2cool 1.0.1 and writes its
pixels as [chr1, chr2, bin1, bin2, count] rows, chromosomes numbered in the
order of chroms.json (the .hic file's order without "All").

    python harness/hic2cool_readback.py FILE.hic RESOLUTION CHROMS.json OUT.npy WORK.cool

Runs in an environment with hic2cool, h5py and numpy.
"""

import json
import sys

import h5py
import numpy as np
from hic2cool import hic2cool_convert


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    hic, resolution, chroms_path, out, cool = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
    hic2cool_convert(hic, cool, resolution, nproc=1, silent=True)
    with open(chroms_path) as handle:
        index = {name: i for i, (name, _) in enumerate(json.load(handle))}
    with h5py.File(cool, "r") as h5:
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in h5["chroms/name"][:]]
        chrom = np.array([index[names[c]] for c in h5["bins/chrom"][:]], dtype=np.int64)
        start = h5["bins/start"][:].astype(np.int64)
        bin1 = h5["pixels/bin1_id"][:].astype(np.int64)
        bin2 = h5["pixels/bin2_id"][:].astype(np.int64)
        count = h5["pixels/count"][:].astype(np.float64)
    rows = np.stack([chrom[bin1], chrom[bin2], start[bin1] // resolution, start[bin2] // resolution, count], axis=1)
    np.save(out, rows.astype(np.float64))


if __name__ == "__main__":
    main()
