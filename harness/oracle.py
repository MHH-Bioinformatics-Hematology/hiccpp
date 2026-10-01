#!/usr/bin/env python
"""The hicstraw side of the hiccpp harness: runs one case through
hicstraw 1.3.1 and writes the result in the layout of hiccpp-harness.

    python harness/oracle.py run CASE.json OUT_DIR

Runs in its own environment with hicstraw installed (built from aidenlab/straw,
pybind11_python); see harness/run.py.
"""

import json
import os
import sys
import time

import numpy as np

import hicstraw


def op_header(case, out):
    handle = hicstraw.HiCFile(case["file"])
    return {
        "kind": "header",
        "genome": handle.getGenomeID(),
        "resolutions": [int(r) for r in handle.getResolutions()],
        "chromosomes": [[c.name, int(c.index), int(c.length)] for c in handle.getChromosomes()],
    }


def op_records(case, out):
    handle = hicstraw.HiCFile(case["file"])
    bin_x, bin_y, counts, sizes = [], [], [], []
    for norm, chr1, chr2, gx0, gx1, gy0, gy1 in case["queries"]:
        mzd = handle.getMatrixZoomData(chr1, chr2, case["matrix_type"], norm, case["unit"], case["resolution"])
        records = mzd.getRecords(gx0, gx1, gy0, gy1)
        sizes.append(len(records))
        bin_x.extend(r.binX for r in records)
        bin_y.extend(r.binY for r in records)
        counts.extend(r.counts for r in records)
    np.save(os.path.join(out, "binX.npy"), np.asarray(bin_x, dtype=np.int32))
    np.save(os.path.join(out, "binY.npy"), np.asarray(bin_y, dtype=np.int32))
    np.save(os.path.join(out, "counts.npy"), np.asarray(counts, dtype=np.float32))
    return {"kind": "records", "n": sizes}


def op_matrices(case, out):
    handle = hicstraw.HiCFile(case["file"])
    shapes = []
    for k, (chr1, chr2, matrix_type, norm, unit, resolution, gx0, gx1, gy0, gy1) in enumerate(case["queries"]):
        mzd = handle.getMatrixZoomData(chr1, chr2, matrix_type, norm, unit, resolution)
        matrix = np.ascontiguousarray(np.asarray(mzd.getRecordsAsMatrix(gx0, gx1, gy0, gy1), dtype=np.float32))
        if matrix.ndim != 2:
            matrix = matrix.reshape(1, 1)
        np.save(os.path.join(out, f"matrix_{k}.npy"), matrix)
        shapes.append(list(matrix.shape))
    return {"kind": "matrices", "shapes": shapes}


def op_vectors(case, out):
    handle = hicstraw.HiCFile(case["file"])
    lengths = []
    for k, (norm, name, index) in enumerate(case["queries"]):
        mzd = handle.getMatrixZoomData(name, name, case["matrix_type"], norm, case["unit"], case["resolution"])
        expected = np.asarray(mzd.getExpectedValues(), dtype=np.float64).ravel()
        vector = np.asarray(mzd.getNormVector(index), dtype=np.float64).ravel()
        np.save(os.path.join(out, f"expected_{k}.npy"), expected)
        np.save(os.path.join(out, f"norm_{k}.npy"), vector)
        lengths.append([int(expected.size), int(vector.size)])
    return {"kind": "vectors", "lengths": lengths}


OPS = {"header": op_header, "records": op_records, "matrices": op_matrices, "vectors": op_vectors}


def main():
    if len(sys.argv) != 4 or sys.argv[1] != "run":
        sys.exit(__doc__)
    with open(sys.argv[2]) as handle:
        case = json.load(handle)
    out = sys.argv[3]
    start = time.perf_counter()
    try:
        result = OPS[case["op"]](case, out)
    except Exception as error:  # noqa: BLE001 - reported to the comparator
        result = {"kind": "error", "type": type(error).__name__, "message": str(error)}
    seconds = time.perf_counter() - start
    with open(os.path.join(out, "result.json"), "w") as handle:
        json.dump({"op_seconds": seconds, "result": result}, handle)


if __name__ == "__main__":
    main()
