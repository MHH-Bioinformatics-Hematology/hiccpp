"""Comparators of the hicfilecpp harness.

compare_results(py, cpp, py_dir, cpp_dir) compares the result documents and
arrays that oracle.py and hicfilecpp-harness write for the same case. Contact
records and matrices must be identical, bit for bit. Vectors are classified by
their worst relative difference and must reach the ED gate (three significant
digits, rel <= 1e-3).
"""

import math
import os

import numpy as np

FLOAT_CLASSES = [
    ("bit-exact", 0.0),
    ("rel<=1e-15", 1e-15),
    ("rel<=1e-12", 1e-12),
    ("rel<=1e-9", 1e-9),
    ("rel<=1e-6", 1e-6),
    ("ED", 1e-3),
]
CLASS_RANK = {name: i for i, (name, _) in enumerate(FLOAT_CLASSES)}


class Comparison:
    def __init__(self):
        self.problems = []
        self.float_class = None

    def fail(self, where, text):
        self.problems.append(f"{where}: {text}")

    def floats(self, klass):
        if self.float_class is None or CLASS_RANK[klass] > CLASS_RANK[self.float_class]:
            self.float_class = klass

    @property
    def ok(self):
        return not self.problems


def float_agreement(a, b):
    """The strictest class two float arrays of equal shape reach, or None.
    An exact zero on one side must be an exact zero on the other."""
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    if a.shape != b.shape:
        return None
    nan_a, nan_b = np.isnan(a), np.isnan(b)
    if not np.array_equal(nan_a, nan_b):
        return None
    keep = ~nan_a
    a, b = a[keep], b[keep]
    if a.tobytes() == b.tobytes():
        return "bit-exact"
    if not np.array_equal(a == 0, b == 0):
        return None
    inf = np.isinf(a) | np.isinf(b)
    if not np.array_equal(a[inf], b[inf]):
        return None
    a, b = a[~inf], b[~inf]
    diff = np.abs(a - b)
    scale = np.abs(b)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(diff == 0, 0.0, diff / np.where(scale == 0, 1.0, scale))
    worst = float(rel.max()) if rel.size else 0.0
    for name, bound in FLOAT_CLASSES[1:]:
        if worst <= bound:
            return name
    return None


def _load(directory, name):
    return np.load(os.path.join(directory, name), allow_pickle=False)


def _exact(cmp, where, a, b):
    if a.dtype != b.dtype or a.shape != b.shape:
        cmp.fail(where, f"dtype/shape {a.dtype}{a.shape} vs {b.dtype}{b.shape}")
        return
    if a.tobytes() != b.tobytes():
        differing = int(np.count_nonzero(a != b)) if a.shape == b.shape else -1
        cmp.fail(where, f"{differing} of {a.size} values differ")


def compare_results(py, cpp, py_dir, cpp_dir, cmp=None):
    cmp = cmp or Comparison()
    if py.get("kind") != cpp.get("kind"):
        cmp.fail("kind", f"hicstraw {py.get('kind')} ({py.get('message', '')}) vs "
                         f"hicfilecpp {cpp.get('kind')} ({cpp.get('message', '')})")
        return cmp
    kind = py["kind"]
    if kind == "error":
        return cmp
    if kind == "header":
        for key in ("genome", "resolutions", "chromosomes"):
            if py[key] != cpp[key]:
                cmp.fail(key, f"{py[key]!r} vs {cpp[key]!r}")
        return cmp
    if kind == "records":
        if py["n"] != cpp["n"]:
            differing = [i for i, (a, b) in enumerate(zip(py["n"], cpp["n"])) if a != b]
            cmp.fail("records per query", f"differ at queries {differing[:10]}")
            return cmp
        for name in ("binX.npy", "binY.npy", "counts.npy"):
            _exact(cmp, name, _load(py_dir, name), _load(cpp_dir, name))
        cmp.floats("bit-exact")
        return cmp
    if kind == "matrices":
        if py["shapes"] != cpp["shapes"]:
            cmp.fail("matrix shapes", f"{py['shapes']} vs {cpp['shapes']}")
            return cmp
        for k in range(len(py["shapes"])):
            name = f"matrix_{k}.npy"
            _exact(cmp, name, _load(py_dir, name), _load(cpp_dir, name))
        cmp.floats("bit-exact")
        return cmp
    if kind == "vectors":
        if py["lengths"] != cpp["lengths"]:
            cmp.fail("vector lengths", f"{py['lengths']} vs {cpp['lengths']}")
            return cmp
        for k in range(len(py["lengths"])):
            for stem in ("expected", "norm"):
                name = f"{stem}_{k}.npy"
                klass = float_agreement(_load(cpp_dir, name), _load(py_dir, name))
                if klass is None:
                    cmp.fail(name, "values differ beyond the ED gate")
                else:
                    cmp.floats(klass)
        return cmp
    cmp.fail("kind", f"unknown result kind {kind}")
    return cmp


def finite_or_nan(value):
    return value if isinstance(value, (int, float)) and math.isfinite(value) else float("nan")
