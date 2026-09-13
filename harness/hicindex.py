"""A small, independent .hic index reader used only to expand harness cases:
which chromosomes, resolutions, zoom levels, normalization vectors and
normalized expected values a file holds. Its output never enters a
comparison; it only selects the queries both sides answer.
"""

import struct


class _Reader:
    def __init__(self, handle, offset):
        self.handle = handle
        self.handle.seek(offset)

    def get(self, fmt):
        size = struct.calcsize("<" + fmt)
        data = self.handle.read(size)
        if len(data) < size:
            raise EOFError
        return struct.unpack("<" + fmt, data)[0]

    def cstr(self):
        out = bytearray()
        while True:
            byte = self.handle.read(1)
            if not byte:
                raise EOFError
            if byte == b"\0":
                return out.decode("utf-8")
            out += byte

    def skip(self, n):
        self.handle.seek(n, 1)

    def at_end(self):
        position = self.handle.tell()
        more = self.handle.read(1)
        self.handle.seek(position)
        return not more


def read_index(path):
    index = {"chromosomes": [], "bp": [], "frag": [], "frag_sites": [], "matrices": {},
             "zooms": {}, "norm_vectors": set(), "normalized_expected": set()}
    with open(path, "rb") as handle:
        r = _Reader(handle, 0)
        r.cstr()
        version = r.get("i")
        index["version"] = version
        master = r.get("q")
        index["genome"] = r.cstr()
        if version > 8:
            r.get("q")
            r.get("q")
        for _ in range(r.get("i")):
            r.cstr()
            r.cstr()
        for i in range(r.get("i")):
            name = r.cstr()
            length = r.get("q") if version > 8 else r.get("i")
            index["chromosomes"].append((name, i, length))
        index["bp"] = [r.get("i") for _ in range(r.get("i"))]
        index["frag"] = [r.get("i") for _ in range(r.get("i"))]
        if index["frag"]:
            for _ in index["chromosomes"]:
                n = r.get("i")
                index["frag_sites"].append(n)
                r.skip(4 * n)

        f = _Reader(handle, master)
        f.get("q" if version > 8 else "i")
        for _ in range(f.get("i")):
            key = f.cstr()
            position = f.get("q")
            f.get("i")
            index["matrices"][key] = position
        value = 4 if version > 8 else 8
        try:
            for normalized in (False, True):
                if f.at_end():
                    break
                for _ in range(f.get("i")):
                    norm = f.cstr() if normalized else "NONE"
                    unit = f.cstr()
                    bin_size = f.get("i")
                    n = f.get("q") if version > 8 else f.get("i")
                    f.skip(n * value)
                    factors = f.get("i")
                    f.skip(factors * (4 + value))
                    if normalized:
                        index["normalized_expected"].add((norm, unit, bin_size))
            if not f.at_end():
                for _ in range(f.get("i")):
                    norm = f.cstr()
                    chrom = f.get("i")
                    unit = f.cstr()
                    resolution = f.get("i")
                    f.get("q")
                    f.get("q" if version > 8 else "i")
                    index["norm_vectors"].add((norm, chrom, unit, resolution))
        except EOFError:
            pass

        for key, position in index["matrices"].items():
            m = _Reader(handle, position)
            m.get("i")
            m.get("i")
            zooms = []
            for _ in range(m.get("i")):
                unit = m.cstr()
                m.skip(4 + 16)
                bin_size = m.get("i")
                m.skip(8)
                blocks = m.get("i")
                m.skip(16 * blocks)
                zooms.append((unit, bin_size))
            index["zooms"][key] = zooms
    return index
