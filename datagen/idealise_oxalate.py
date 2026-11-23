#!/usr/bin/env python3
"""Replace the oxalate of M(ox)(NH3)4 start structures by an ideal chelate ring.

GFN2-xTB (and the molSimplify input) stretched the oxalate C-C bond to 1.7-1.9 A with C-O 1.18 A, i.e. towards
2 CO2; the r2SCAN-3c optimisations then broke the C-C bond. This script rebuilds the planar M-O-C-C-O ring in
the plane of the metal and the two coordinating O atoms (C-C 1.55, C-O 1.28 / 1.23 A, O-C-C 115 / 118 deg,
M-O 2.0 A) and leaves all other atoms in place. The files are marked with NO_XTB so preopt_xtb.py keeps them.

    python datagen/idealise_oxalate.py structures/xyz/fe3_nh3_4ox.xyz ...
"""
import math
import sys
from pathlib import Path

from geom import read_xyz

NO_XTB = "no xTB preopt: oxalate idealised by idealise_oxalate.py"
CC, CO_IN, CO_OUT, MO = 1.55, 1.28, 1.23, 2.0
A_IN, A_OUT = math.radians(115), math.radians(118)  # O(coord)-C-C and O(outer)-C-C


def ring_2d():
    """Ideal ring in 2D: C-C along y, centred at the origin, metal on the -x axis. Returns {label: (x, y)}."""
    c1 = (0.0, CC / 2)
    o1 = (-CO_IN * math.sin(A_IN), CC / 2 - CO_IN * math.cos(A_IN))
    oo = (CO_OUT * math.sin(A_OUT), CC / 2 - CO_OUT * math.cos(A_OUT))
    xm = o1[0] - math.sqrt(MO**2 - o1[1] ** 2)
    pts = {"C_a": c1, "O_a": o1, "Oo_a": oo, "M": (xm, 0.0)}
    pts.update({k.replace("_a", "_b"): (x, -y) for k, (x, y) in list(pts.items()) if k != "M"})
    return pts


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def unit(v):
    n = math.sqrt(sum(x * x for x in v))
    return tuple(x / n for x in v)


def idealise(el, xyz):
    """New coordinates with an ideal oxalate ring; the metal is atom 0."""
    m = xyz[0]
    cs = [i for i, e in enumerate(el) if e == "C"]
    assert len(cs) == 2, "expected exactly one oxalate and no other carbon"
    ox_o = [i for i, e in enumerate(el) if e == "O" and min(math.dist(xyz[i], xyz[c]) for c in cs) < 1.6]
    o_a, o_b = sorted(ox_o, key=lambda i: math.dist(xyz[i], m))[:2]  # coordinating O
    c_a = min(cs, key=lambda c: math.dist(xyz[c], xyz[o_a]))
    c_b = next(c for c in cs if c != c_a)
    outer = [i for i in ox_o if i not in (o_a, o_b)]
    oo_a = min(outer, key=lambda i: math.dist(xyz[i], xyz[c_a]))
    oo_b = next(i for i in outer if i != oo_a)

    u = unit(sub(tuple((p + q) / 2 for p, q in zip(xyz[o_a], xyz[o_b])), m))  # metal -> ring
    d = sub(xyz[o_a], xyz[o_b])
    v = unit(sub(d, tuple(sum(x * y for x, y in zip(d, u)) * x for x in u)))
    pts = ring_2d()
    origin = tuple(mi - pts["M"][0] * ui for mi, ui in zip(m, u))
    new = list(xyz)
    for label, idx in {"C_a": c_a, "C_b": c_b, "O_a": o_a, "O_b": o_b, "Oo_a": oo_a, "Oo_b": oo_b}.items():
        x, y = pts[label]
        new[idx] = tuple(o + x * ui + y * vi for o, ui, vi in zip(origin, u, v))
    ox = {c_a, c_b, o_a, o_b, oo_a, oo_b}
    clash = min(math.dist(new[i], new[j]) for i in ox for j in range(1, len(el)) if j not in ox)
    return new, clash


def main(paths):
    for path in map(Path, paths):
        el, xyz, _ = read_xyz(path)
        new, clash = idealise(el, xyz)
        body = [f"{e} {x:12.6f} {y:12.6f} {z:12.6f}" for e, (x, y, z) in zip(el, new)]
        path.write_text("\n".join([str(len(el)), NO_XTB, *body]) + "\n")
        print(f"{path.name}: oxalate idealised, closest contact to other ligands {clash:.2f} A")


if __name__ == "__main__":
    main(sys.argv[1:])
