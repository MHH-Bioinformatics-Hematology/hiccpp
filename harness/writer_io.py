#!/usr/bin/env python
"""Inputs and read-back checks of the hicfilecpp writer cases. Runs in the
oracle environment (hicstraw, h5py, numpy).

    python harness/writer_io.py prep SOURCE SOURCE_RESOLUTION OUT_DIR
    python harness/writer_io.py check CHECK.json OUT.json

prep turns a .hic file (read with hicstraw) or a fixed-bin .cool file (read
with h5py) into the same contacts for both writers:
  chroms.json   [[name, length], ...] in file order
  pixels.bin    records (chr1, chr2, bin1, bin2: int32; count: float32),
                grouped by chromosome pair in pair order, chr1 <= chr2
  pairs.json    [[chr1, chr2, first record, end record], ...]
  genome.chrom.sizes, contacts.txt
                Juicer tools pre input, one "short with score" line per pixel
                at the start of its bins

check reads a hicfilecpp-written file and the Juicer-written file of the same
contacts back through hicstraw and compares, at every resolution and for every
chromosome pair, the observed records with each other and, at the source
resolution, with the source pixels; for every normalization and chromosome,
the normalization vectors and the expected values; and the raw expected
values. Records must be identical (bit for bit) when the source counts are
integers, and agree to the ED gate otherwise; vectors must reach the ED gate.
"""

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import compare  # noqa: E402
import hicindex  # noqa: E402

PIXEL = np.dtype([("chr1", "<i4"), ("chr2", "<i4"), ("bin1", "<i4"), ("bin2", "<i4"), ("count", "<f4")])


def _from_hic(path, resolution):
    import hicstraw

    handle = hicstraw.HiCFile(path)
    chroms = [(c.name, int(c.length)) for c in handle.getChromosomes() if c.index > 0]
    parts = []
    for i, (name1, length1) in enumerate(chroms):
        for j in range(i, len(chroms)):
            name2, length2 = chroms[j]
            mzd = handle.getMatrixZoomData(name1, name2, "observed", "NONE", "BP", resolution)
            records = mzd.getRecords(0, length1, 0, length2)
            if not records:
                continue
            part = np.zeros(len(records), dtype=PIXEL)
            part["chr1"], part["chr2"] = i, j
            part["bin1"] = [r.binX // resolution for r in records]
            part["bin2"] = [r.binY // resolution for r in records]
            part["count"] = [r.counts for r in records]
            part.sort(order=["bin1", "bin2"])
            parts.append(part)
    return chroms, np.concatenate(parts)


def _from_cool(path, resolution):
    import h5py

    with h5py.File(path, "r") as h5:
        group = h5
        if "resolutions" in h5:
            group = h5["resolutions"][str(resolution)]
        if int(group.attrs["bin-size"]) != resolution:
            raise SystemExit(f"{path} has bin size {group.attrs['bin-size']}, not {resolution}")
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in group["chroms/name"][:]]
        lengths = [int(x) for x in group["chroms/length"][:]]
        chrom_of_bin = group["bins/chrom"][:].astype(np.int64)
        start_of_bin = group["bins/start"][:].astype(np.int64)
        bin1 = group["pixels/bin1_id"][:].astype(np.int64)
        bin2 = group["pixels/bin2_id"][:].astype(np.int64)
        count = group["pixels/count"][:]
    pixels = np.zeros(len(bin1), dtype=PIXEL)
    c1, c2 = chrom_of_bin[bin1], chrom_of_bin[bin2]
    swap = c1 > c2
    pixels["chr1"] = np.where(swap, c2, c1)
    pixels["chr2"] = np.where(swap, c1, c2)
    local1 = start_of_bin[bin1] // resolution
    local2 = start_of_bin[bin2] // resolution
    pixels["bin1"] = np.where(swap, local2, local1)
    pixels["bin2"] = np.where(swap, local1, local2)
    if not np.all(np.isfinite(count)):
        raise SystemExit(f"{path} has non-finite counts")
    pixels["count"] = count.astype(np.float32)
    order = np.lexsort((pixels["bin2"], pixels["bin1"], pixels["chr2"], pixels["chr1"]))
    return list(zip(names, lengths)), pixels[order]


def prep(source, resolution, out):
    os.makedirs(out, exist_ok=True)
    if source.endswith(".hic"):
        chroms, pixels = _from_hic(source, resolution)
    else:
        chroms, pixels = _from_cool(source, resolution)
    pixels.tofile(os.path.join(out, "pixels.bin"))
    keys = pixels["chr1"].astype(np.int64) * 1000000 + pixels["chr2"]
    starts = np.flatnonzero(np.r_[True, keys[1:] != keys[:-1]])
    ends = np.r_[starts[1:], len(pixels)]
    pairs = [[int(pixels["chr1"][s]), int(pixels["chr2"][s]), int(s), int(e)] for s, e in zip(starts, ends)]
    with open(os.path.join(out, "pairs.json"), "w") as handle:
        json.dump(pairs, handle)
    with open(os.path.join(out, "chroms.json"), "w") as handle:
        json.dump([[n, l] for n, l in chroms], handle)
    with open(os.path.join(out, "genome.chrom.sizes"), "w") as handle:
        for name, length in chroms:
            handle.write(f"{name}\t{length}\n")
    names = [n for n, _ in chroms]
    with open(os.path.join(out, "contacts.txt"), "w") as handle:
        for chr1, chr2, first, end in pairs:
            part = pixels[first:end]
            a, b = names[chr1], names[chr2]
            lines = [f"0 {a} {x * resolution} 0 0 {b} {y * resolution} 1 {float(c):.9g}\n"
                     for x, y, c in zip(part["bin1"], part["bin2"], part["count"])]
            handle.writelines(lines)
    integer = bool(np.all(np.floor(pixels["count"]) == pixels["count"]))
    with open(os.path.join(out, "source.json"), "w") as handle:
        json.dump({"pixels": int(len(pixels)), "integer_counts": integer, "pairs": len(pairs)}, handle)
    print(json.dumps({"pixels": int(len(pixels)), "integer_counts": integer, "pairs": len(pairs)}))


def _records(handle, chroms, i, j, resolution, norm="NONE"):
    (a, la), (b, lb) = chroms[i], chroms[j]
    records = handle.getMatrixZoomData(a, b, "observed", norm, "BP", resolution).getRecords(0, la, 0, lb)
    out = np.zeros(len(records), dtype=[("x", "<i4"), ("y", "<i4"), ("c", "<f4")])
    out["x"] = [r.binX for r in records]
    out["y"] = [r.binY for r in records]
    out["c"] = [r.counts for r in records]
    out.sort(order=["x", "y"])
    return out


def check(spec, out_path):
    import hicstraw

    cmp = compare.Comparison()
    counts = {"record_sets": 0, "records": 0, "vectors": 0, "expected": 0}
    inputs = spec["inputs"]
    with open(os.path.join(inputs, "chroms.json")) as handle:
        chroms = [tuple(c) for c in json.load(handle)]
    with open(os.path.join(inputs, "source.json")) as handle:
        integer = json.load(handle)["integer_counts"]
    with open(os.path.join(inputs, "pairs.json")) as handle:
        pairs = json.load(handle)
    pixels = np.fromfile(os.path.join(inputs, "pixels.bin"), dtype=PIXEL)
    candidate = hicstraw.HiCFile(spec["candidate"])
    reference = hicstraw.HiCFile(spec["reference"])
    source_resolution = spec["source_resolution"]

    header_c = [(c.name, int(c.length)) for c in candidate.getChromosomes()]
    header_r = [(c.name, int(c.length)) for c in reference.getChromosomes()]
    if header_c[1:] != header_r[1:] or header_c[1:] != chroms:
        cmp.fail("chromosomes", f"{header_c[:4]} vs {header_r[:4]}")
    if list(candidate.getResolutions()) != list(reference.getResolutions()):
        cmp.fail("resolutions", f"{candidate.getResolutions()} vs {reference.getResolutions()}")

    offsets = {}
    for chr1, chr2, first, end in pairs:
        offsets[(chr1, chr2)] = (first, end)
    for resolution in reference.getResolutions():
        for i in range(len(chroms)):
            for j in range(i, len(chroms)):
                if (i, j) not in offsets:
                    continue
                mine = _records(candidate, chroms, i, j, resolution)
                theirs = _records(reference, chroms, i, j, resolution)
                counts["record_sets"] += 1
                counts["records"] += len(mine)
                where = f"records {chroms[i][0]} {chroms[j][0]} {resolution}"
                if len(mine) != len(theirs) or not (np.array_equal(mine["x"], theirs["x"])
                                                   and np.array_equal(mine["y"], theirs["y"])):
                    cmp.fail(where, f"{len(mine)} records vs Juicer's {len(theirs)}, or different pixels")
                    continue
                if integer:
                    if mine["c"].tobytes() != theirs["c"].tobytes():
                        cmp.fail(where, "counts differ from Juicer's")
                    else:
                        cmp.floats("bit-exact")
                else:
                    klass = compare.float_agreement(mine["c"], theirs["c"])
                    if klass is None:
                        cmp.fail(where, "counts differ from Juicer's beyond the ED gate")
                    else:
                        cmp.floats(klass)
                if resolution == source_resolution:
                    first, end = offsets[(i, j)]
                    part = pixels[first:end]
                    source = np.zeros(len(part), dtype=mine.dtype)
                    source["x"] = np.minimum(part["bin1"], part["bin2"]) * resolution if i == j else part["bin1"] * resolution
                    source["y"] = np.maximum(part["bin1"], part["bin2"]) * resolution if i == j else part["bin2"] * resolution
                    source["c"] = part["count"]
                    source.sort(order=["x", "y"])
                    if source.tobytes() != mine.tobytes():
                        cmp.fail(where, "records differ from the source pixels")

    index_c = hicindex.read_index(spec["candidate"])
    index_r = hicindex.read_index(spec["reference"])
    wanted = set(spec["normalizations"])
    vectors_c = {v for v in index_c["norm_vectors"] if v[0] in wanted}
    vectors_r = {v for v in index_r["norm_vectors"] if v[0] in wanted}
    if vectors_c != vectors_r:
        cmp.fail("normalization vectors", f"{len(vectors_c)} entries vs Juicer's {len(vectors_r)}: "
                                          f"{sorted(vectors_c ^ vectors_r)[:6]}")
    worst = {}
    for norm, chrom, unit, resolution in sorted(vectors_c & vectors_r):
        name = chroms[chrom - 1][0]
        for handle_name, handle in (("candidate", candidate), ("reference", reference)):
            pass
        mzd_c = candidate.getMatrixZoomData(name, name, "oe", norm, unit, resolution)
        mzd_r = reference.getMatrixZoomData(name, name, "oe", norm, unit, resolution)
        for kind, a, b in (("vector", mzd_c.getNormVector(chrom), mzd_r.getNormVector(chrom)),
                           ("expected", mzd_c.getExpectedValues(), mzd_r.getExpectedValues())):
            a = np.asarray(a, dtype=np.float64)
            b = np.asarray(b, dtype=np.float64)
            counts["vectors" if kind == "vector" else "expected"] += 1
            klass = compare.float_agreement(a, b) if a.shape == b.shape else None
            where = f"{kind} {norm} {name} {resolution}"
            if klass is None:
                cmp.fail(where, f"differs beyond the ED gate (lengths {a.size} and {b.size})")
            else:
                cmp.floats(klass)
                worst[f"{kind} {norm}"] = max(worst.get(f"{kind} {norm}", "bit-exact"), klass,
                                              key=lambda k: compare.CLASS_RANK[k])
    for resolution in reference.getResolutions():
        for i, (name, _) in enumerate(chroms):
            if (i, i) not in offsets:
                continue
            a = np.asarray(candidate.getMatrixZoomData(name, name, "oe", "NONE", "BP", resolution).getExpectedValues())
            b = np.asarray(reference.getMatrixZoomData(name, name, "oe", "NONE", "BP", resolution).getExpectedValues())
            counts["expected"] += 1
            klass = compare.float_agreement(a, b) if a.shape == b.shape else None
            if klass is None:
                cmp.fail(f"expected NONE {name} {resolution}", f"differs (lengths {a.size} and {b.size})")
            else:
                cmp.floats(klass)
                worst["expected NONE"] = max(worst.get("expected NONE", "bit-exact"), klass,
                                             key=lambda k: compare.CLASS_RANK[k])
    with open(out_path, "w") as handle:
        json.dump({"problems": cmp.problems, "float_class": cmp.float_class, "counts": counts,
                   "worst_by_kind": worst}, handle, indent=1)


def compare_pixels(inputs, other_path, pairs, out_path):
    """Compares pixels read back by another tool, [chr1, chr2, bin1, bin2,
    count] rows in a .npy file, with the source pixels of the chromosome pairs
    named in `pairs` ("all", or "0_0,0_3")."""
    pixels = np.fromfile(os.path.join(inputs, "pixels.bin"), dtype=PIXEL)
    other = np.load(other_path).reshape(-1, 5)
    if pairs == "all":
        wanted = {(int(a), int(b)) for a, b in zip(pixels["chr1"], pixels["chr2"])}
    else:
        wanted = {tuple(int(v) for v in p.split("_")) for p in pairs.split(",")}
    problems = []
    intra = pixels["chr1"] == pixels["chr2"]
    a = np.stack([pixels["chr1"], pixels["chr2"],
                  np.where(intra, np.minimum(pixels["bin1"], pixels["bin2"]), pixels["bin1"]),
                  np.where(intra, np.maximum(pixels["bin1"], pixels["bin2"]), pixels["bin2"])],
                 axis=1).astype(np.int64)
    mask = np.array([(int(x), int(y)) in wanted for x, y in a[:, :2]]) if len(a) else np.zeros(0, bool)
    a = a[mask]
    a_counts = pixels["count"][mask]
    order = np.lexsort((a[:, 3], a[:, 2], a[:, 1], a[:, 0]))
    b = other[:, :4].astype(np.int64)
    b_counts = other[:, 4].astype(np.float32)
    order_b = np.lexsort((b[:, 3], b[:, 2], b[:, 1], b[:, 0]))
    if len(a) != len(b):
        problems.append(f"{len(b)} pixels read back, source has {len(a)}")
    elif not np.array_equal(a[order], b[order_b]):
        problems.append("pixel coordinates differ from the source")
    elif a_counts[order].tobytes() != b_counts[order_b].tobytes():
        problems.append("pixel counts differ from the source")
    with open(out_path, "w") as handle:
        json.dump({"problems": problems, "pixels": int(len(b)), "pairs": len(wanted)}, handle)


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == "prep" and len(sys.argv) == 5:
        prep(sys.argv[2], int(sys.argv[3]), sys.argv[4])
    elif len(sys.argv) == 4 and sys.argv[1] == "check":
        with open(sys.argv[2]) as handle:
            check(json.load(handle), sys.argv[3])
    elif len(sys.argv) == 6 and sys.argv[1] == "compare-pixels":
        compare_pixels(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
