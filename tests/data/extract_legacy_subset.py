#!/usr/bin/env python
"""Extracts the version 6 and version 7 test files in this directory from a
real version 7 .hic file, GSE63525_GM12878_insitu_primary+replicate_combined_30.hic
(Rao et al. 2014, GEO GSE63525, written by the Juicebox tools of 2014).

    python tests/data/extract_legacy_subset.py SOURCE.hic OUT_DIR \\
        [--chromosomes 21 22] [--resolutions 2500000 1000000 500000 250000]

No Juicer tools release that writes version 6 or 7 can be obtained any more
(docs/PROVENANCE.md), so the test files are cut out of a deposited file
instead of being written from pixels:

  GM12878_combined_30.chr21_chr22.v7.hic
      A byte-level subset of the source. The header keeps the genome, the
      attributes and the chromosome dictionary restricted to All and the
      chosen chromosomes (renumbered 0, 1, 2, ...), the chosen BP resolutions
      and no fragment resolutions. Every chromosome pair of the chosen
      chromosomes keeps its matrix at the chosen resolutions: the zoom
      metadata as stored, and every block's compressed bytes copied verbatim.
      The footer keeps the raw and normalized expected value vectors of the
      chosen resolutions as stored, with their per-chromosome normalization
      factors restricted to the kept chromosomes, and every normalization
      vector of a kept chromosome at a kept resolution, copied verbatim. Only
      file positions, sizes and chromosome indices are rewritten.

  GM12878_combined_30.chr21_chr22.v6.hic
      The same file with version 6 in the header and every block re-encoded
      in the version 6 record layout (int32 nRecords, then int32 binX,
      int32 binY, float32 counts per record; readBlock in hicstraw and
      Juicebox for version < 7): each version 7 block is inflated, decoded in
      hicstraw's order and deflated again with zlib's default level. No real
      version 6 deposit is at hand; this file carries the real pixels of the
      version 7 file in the version 6 block layout, which is the only place
      where the readers treat 6 differently from 7.

The layout is the one hicstraw 1.3.1 and Juicebox's DatasetReaderV2 read for
versions below 8: int32 sizes in the footer, 8 byte expected and vector
values, 'nBytes' spanning the master index and the raw expected values.
"""

import argparse
import os
import struct
import zlib


class Reader:
    def __init__(self, handle, offset):
        self.handle = handle
        handle.seek(offset)

    def get(self, fmt):
        size = struct.calcsize("<" + fmt)
        data = self.handle.read(size)
        if len(data) != size:
            raise EOFError
        return struct.unpack("<" + fmt, data)[0]

    def raw(self, n):
        data = self.handle.read(n)
        if len(data) != n:
            raise EOFError
        return data

    def cstr(self):
        out = bytearray()
        while True:
            byte = self.handle.read(1)
            if not byte:
                raise EOFError
            if byte == b"\0":
                return out.decode("utf-8")
            out += byte

    def tell(self):
        return self.handle.tell()


class Writer:
    def __init__(self):
        self.data = bytearray()

    def put(self, fmt, value):
        self.data += struct.pack("<" + fmt, value)

    def cstr(self, text):
        self.data += text.encode("utf-8") + b"\0"

    def raw(self, data):
        self.data += data

    def tell(self):
        return len(self.data)


def read_source(path):
    src = {}
    handle = open(path, "rb")
    r = Reader(handle, 0)
    if r.cstr() != "HIC":
        raise SystemExit("not a .hic file")
    src["version"] = r.get("i")
    if src["version"] != 7:
        raise SystemExit(f"expected a version 7 source, got {src['version']}")
    src["master"] = r.get("q")
    src["genome"] = r.cstr()
    src["attributes"] = [(r.cstr(), r.cstr()) for _ in range(r.get("i"))]
    src["chromosomes"] = [(r.cstr(), r.get("i")) for _ in range(r.get("i"))]
    src["bp"] = [r.get("i") for _ in range(r.get("i"))]

    f = Reader(handle, src["master"])
    f.get("i")  # nBytes
    src["matrices"] = []
    for _ in range(f.get("i")):
        key = f.cstr()
        src["matrices"].append((key, f.get("q"), f.get("i")))

    def expected(normalized):
        out = []
        for _ in range(f.get("i")):
            norm = f.cstr() if normalized else None
            unit = f.cstr()
            bin_size = f.get("i")
            n = f.get("i")
            values = f.raw(8 * n)
            factors = [(f.get("i"), f.get("d")) for _ in range(f.get("i"))]
            out.append((norm, unit, bin_size, n, values, factors))
        return out

    src["expected"] = expected(False)
    src["normalized_expected"] = expected(True)
    src["norm_index"] = []
    for _ in range(f.get("i")):
        norm = f.cstr()
        chrom = f.get("i")
        unit = f.cstr()
        resolution = f.get("i")
        src["norm_index"].append((norm, chrom, unit, resolution, f.get("q"), f.get("i")))
    return handle, src


def read_zooms(handle, position):
    m = Reader(handle, position)
    c1, c2 = m.get("i"), m.get("i")
    zooms = []
    for _ in range(m.get("i")):
        unit = m.cstr()
        stats = m.raw(20)  # zoomIndex, sumCounts, occupiedCellCount, stdDev, percent95
        bin_size = m.get("i")
        block_bin_count = m.get("i")
        block_column_count = m.get("i")
        blocks = [(m.get("i"), m.get("q"), m.get("i")) for _ in range(m.get("i"))]
        zooms.append((unit, stats, bin_size, block_bin_count, block_column_count, blocks))
    return c1, c2, zooms


def decode_v7_block(data):
    """The records of a version 7 or 8 block in hicstraw's readBlock order."""
    records = []
    n_records, x_offset, y_offset = struct.unpack_from("<iii", data, 0)
    use_short = data[12] == 0
    block_type = data[13]
    pos = 14
    if block_type == 1:
        (row_count,) = struct.unpack_from("<h", data, pos)
        pos += 2
        for _ in range(row_count):
            y, col_count = struct.unpack_from("<hh", data, pos)
            pos += 4
            for _ in range(col_count):
                (x,) = struct.unpack_from("<h", data, pos)
                pos += 2
                if use_short:
                    (c,) = struct.unpack_from("<h", data, pos)
                    pos += 2
                    counts = float(c)
                else:
                    (counts,) = struct.unpack_from("<f", data, pos)
                    pos += 4
                records.append((x_offset + x, y_offset + y, counts))
    elif block_type == 2:
        n_pts, w = struct.unpack_from("<ih", data, pos)
        pos += 6
        for i in range(n_pts):
            row = i // w
            col = i - row * w
            if use_short:
                (c,) = struct.unpack_from("<h", data, pos)
                pos += 2
                if c != -32768:
                    records.append((x_offset + col, y_offset + row, float(c)))
            else:
                (counts,) = struct.unpack_from("<f", data, pos)
                pos += 4
                if counts == counts:
                    records.append((x_offset + col, y_offset + row, counts))
    else:
        raise SystemExit(f"unknown block type {block_type}")
    return records


def encode_v6_block(records):
    raw = bytearray(struct.pack("<i", len(records)))
    for x, y, counts in records:
        raw += struct.pack("<iif", x, y, counts)
    return zlib.compress(bytes(raw))


def write_subset(handle, src, chromosomes, resolutions, version, path):
    old_index = {name: i for i, (name, _) in enumerate(src["chromosomes"])}
    missing = [c for c in chromosomes if c not in old_index]
    if missing:
        raise SystemExit(f"chromosomes not in the source: {missing}")
    kept = [0] + [old_index[c] for c in chromosomes]
    renumber = {old: new for new, old in enumerate(kept)}
    resolutions = [r for r in src["bp"] if r in set(resolutions)]

    out = Writer()
    out.cstr("HIC")
    out.put("i", version)
    master_at = out.tell()
    out.put("q", 0)
    out.cstr(src["genome"])
    out.put("i", len(src["attributes"]))
    for key, value in src["attributes"]:
        out.cstr(key)
        out.cstr(value)
    out.put("i", len(kept))
    for old in kept:
        name, length = src["chromosomes"][old]
        out.cstr(name)
        out.put("i", length)
    out.put("i", len(resolutions))
    for r in resolutions:
        out.put("i", r)
    out.put("i", 0)  # no fragment resolutions

    entries = []
    for key, position, _ in src["matrices"]:
        a, b = (int(x) for x in key.split("_"))
        if a == 0 or a not in renumber or b not in renumber:
            continue
        c1, c2, zooms = read_zooms(handle, position)
        kept_zooms = []
        for unit, stats, bin_size, bbc, bcc, blocks in zooms:
            if unit != "BP" or bin_size not in resolutions:
                continue
            new_blocks = []
            for number, block_position, size in blocks:
                handle.seek(block_position)
                data = handle.read(size)
                if version < 7:
                    data = encode_v6_block(decode_v7_block(zlib.decompress(data)))
                new_blocks.append((number, out.tell(), len(data)))
                out.raw(data)
            kept_zooms.append((unit, stats, bin_size, bbc, bcc, new_blocks))
        metadata = out.tell()
        out.put("i", renumber[c1])
        out.put("i", renumber[c2])
        out.put("i", len(kept_zooms))
        for unit, stats, bin_size, bbc, bcc, blocks in kept_zooms:
            out.cstr(unit)
            out.raw(stats)
            out.put("i", bin_size)
            out.put("i", bbc)
            out.put("i", bcc)
            out.put("i", len(blocks))
            for number, block_position, size in blocks:
                out.put("i", number)
                out.put("q", block_position)
                out.put("i", size)
        entries.append((f"{renumber[c1]}_{renumber[c2]}", metadata, out.tell() - metadata))

    vectors = []
    for norm, chrom, unit, resolution, position, size in src["norm_index"]:
        if chrom not in renumber or unit != "BP" or resolution not in resolutions:
            continue
        handle.seek(position)
        data = handle.read(size)
        vectors.append((norm, renumber[chrom], unit, resolution, out.tell(), len(data)))
        out.raw(data)

    def put_expected(section, normalized):
        chosen = [e for e in section if e[1] == "BP" and e[2] in resolutions]
        body.put("i", len(chosen))
        for norm, unit, bin_size, n, values, factors in chosen:
            if normalized:
                body.cstr(norm)
            body.cstr(unit)
            body.put("i", bin_size)
            body.put("i", n)
            body.raw(values)
            kept_factors = [(renumber[c], v) for c, v in factors if c in renumber]
            body.put("i", len(kept_factors))
            for chrom, value in kept_factors:
                body.put("i", chrom)
                body.put("d", value)

    master = out.tell()
    body = Writer()
    body.put("i", len(entries))
    for key, position, size in entries:
        body.cstr(key)
        body.put("q", position)
        body.put("i", size)
    put_expected(src["expected"], False)
    n_bytes = body.tell()
    put_expected(src["normalized_expected"], True)
    body.put("i", len(vectors))
    for norm, chrom, unit, resolution, position, size in vectors:
        body.cstr(norm)
        body.put("i", chrom)
        body.cstr(unit)
        body.put("i", resolution)
        body.put("q", position)
        body.put("i", size)
    out.put("i", n_bytes)
    out.raw(bytes(body.data))
    struct.pack_into("<q", out.data, master_at, master)
    with open(path, "wb") as handle_out:
        handle_out.write(out.data)
    return len(entries), len(vectors), len(out.data)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("source")
    parser.add_argument("out_dir")
    parser.add_argument("--chromosomes", nargs="+", default=["21", "22"])
    parser.add_argument("--resolutions", nargs="+", type=int, default=[2500000, 1000000, 500000, 250000])
    args = parser.parse_args()
    handle, src = read_source(args.source)
    stem = "GM12878_combined_30." + "_".join("chr" + c for c in args.chromosomes)
    for version in (7, 6):
        path = os.path.join(args.out_dir, f"{stem}.v{version}.hic")
        matrices, vectors, size = write_subset(handle, src, args.chromosomes, args.resolutions, version, path)
        print(f"{path}: {matrices} matrices, {vectors} normalization vectors, {size} bytes")


if __name__ == "__main__":
    main()
