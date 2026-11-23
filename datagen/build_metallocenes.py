#!/usr/bin/env python3
"""Start structures of the metallocenes M(Cp)2, M = Cr, Mn, Fe, Co, Ni (build_structures.HANDBUILT).

Staggered D5d sandwich: C-C 1.42 A, C-H 1.08 A (H in the ring plane), metal-ring centroid distance from gas-phase
electron diffraction (Haaland, Acc. Chem. Res. 1979, 12, 415). Manganocene gets a high-spin start (<id>.xyz) and a
low-spin start for the doublet (<id>_m2.xyz). The files are marked 'no xTB preopt', so preopt_xtb.py keeps them.

    python datagen/build_metallocenes.py structures
"""
import math
import sys
from pathlib import Path

CC, CH = 1.42, 1.08
CENTROID = {"cr2_cp2": 1.80, "mn2_cp2": 2.05, "fe2_cp2": 1.67, "co2_cp2": 1.74, "ni2_cp2": 1.83}  # high-spin Mn
EXTRA_STARTS = {("mn2_cp2", 2): 1.73}  # low-spin manganocene, from (Me5Cp)2Mn (Mn-C ~2.11 A)
NO_XTB = "no xTB preopt: metallocene built by build_metallocenes.py"


def metallocene(symbol, h):
    """(elements, coordinates) of a staggered D5d M(Cp)2 with the metal at the origin and rings at z = +-h."""
    r = CC / (2 * math.sin(math.pi / 5))
    el, xyz = [symbol], [(0.0, 0.0, 0.0)]
    for z, phase in ((h, 0.0), (-h, math.pi / 5)):
        for k in range(5):
            a = 2 * math.pi * k / 5 + phase
            el += ["C", "H"]
            xyz += [(r * math.cos(a), r * math.sin(a), z), ((r + CH) * math.cos(a), (r + CH) * math.sin(a), z)]
    return el, xyz


def write(path, el, xyz):
    body = [f"{e} {x:12.6f} {y:12.6f} {z:12.6f}" for e, (x, y, z) in zip(el, xyz)]
    Path(path).write_text("\n".join([str(len(el)), NO_XTB, *body]) + "\n")


def main(out):
    xyz_dir = Path(out) / "xyz"
    for cid, h in CENTROID.items():
        symbol = cid[:2].capitalize()
        write(xyz_dir / f"{cid}.xyz", *metallocene(symbol, h))
        print(f"{cid}.xyz: M-centroid {h:.2f} A")
    for (cid, mult), h in EXTRA_STARTS.items():
        write(xyz_dir / f"{cid}_m{mult}.xyz", *metallocene(cid[:2].capitalize(), h))
        print(f"{cid}_m{mult}.xyz: M-centroid {h:.2f} A")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "structures")
