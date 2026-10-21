#!/usr/bin/env python3
"""Pre-optimise the start structures of a manifest with GFN2-xTB (see docs/data_plan.md).

Each xyz file is optimised for the state it serves: <id>_m<mult>.xyz for that multiplicity, <id>.xyz for
the highest multiplicity without its own file. A result is accepted only if xtb terminated normally and
the bond graph (including metal-donor bonds, see graph) is unchanged; otherwise the input structure is kept and the
failure is logged in preopt_log.csv. Optimised files are marked in the comment line and skipped on re-runs.
"""
import argparse
import csv
import math
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from make_inputs import read_manifest

TAG = "GFN2-xTB preopt"
# Pyykko single-bond covalent radii (Angstrom)
RCOV = {"H": 0.32, "C": 0.75, "N": 0.71, "O": 0.63, "F": 0.64, "P": 1.11, "S": 1.03, "Cl": 0.99, "Br": 1.14,
        "Cr": 1.22, "Mn": 1.19, "Fe": 1.16, "Co": 1.11, "Ni": 1.10, "Ru": 1.25, "Rh": 1.25}
METALS = {"Cr", "Mn", "Fe", "Co", "Ni", "Ru", "Rh"}


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    atoms = [line.split() for line in lines[2:2 + n]]
    return [a[0] for a in atoms], [tuple(map(float, a[1:4])) for a in atoms], lines[1]


def graph(el, xyz):
    """Bonded index pairs: d < 1.3 (r_A + r_B); 1.35 (r_A + r_B) for metal bonds, which keeps long high-spin/JT
    bonds but not the second-shell atoms of five-membered chelate rings (M...C ~2.7 A)."""
    return {(i, j) for i in range(len(el)) for j in range(i)
            if math.dist(xyz[i], xyz[j]) < (1.35 if {el[i], el[j]} & METALS else 1.3) * (RCOV[el[i]] + RCOV[el[j]])}


def preopt(path, charge, mult, xtb):
    el, xyz, comment = read_xyz(path)
    if TAG in comment:
        return "done earlier"
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(path, Path(tmp) / "in.xyz")
        cmd = [xtb, "in.xyz", "--opt", "--gfn", "2", "--chrg", str(charge), "--uhf", str(mult - 1)]
        try:
            p = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, timeout=3600,
                               env={**os.environ, "OMP_NUM_THREADS": "1"})
        except subprocess.TimeoutExpired:
            return "failed: timeout"
        opt = Path(tmp) / "xtbopt.xyz"
        if "normal termination of xtb" not in p.stdout + p.stderr or not opt.is_file():
            return "failed: xtb"
        el2, xyz2, _ = read_xyz(opt)
        if graph(el, xyz) != graph(el2, xyz2):
            return "failed: bonding changed"
        body = [f"{e} {x:12.6f} {y:12.6f} {z:12.6f}" for e, (x, y, z) in zip(el2, xyz2)]
        Path(path).write_text("\n".join([str(len(el2)), f"{TAG} charge={charge} mult={mult}", *body]) + "\n")
    return "ok"


def tasks(row):
    """(xyz file, multiplicity) pairs for one manifest row."""
    xyz, mults = row["xyz"], sorted(row["mults"])
    out = [(xyz.with_name(f"{xyz.stem}_m{m}.xyz"), m) for m in mults if xyz.with_name(f"{xyz.stem}_m{m}.xyz").is_file()]
    rest = [m for m in mults if m not in {m for _, m in out}]
    return out + ([(xyz, max(rest))] if rest else [])


def main(argv=None):
    ap = argparse.ArgumentParser(description="Pre-optimise manifest structures with GFN2-xTB.")
    ap.add_argument("manifest", help="complexes.csv")
    ap.add_argument("--xtb", default="xtb", help="xtb executable")
    ap.add_argument("--nprocs", type=int, default=4, help="parallel xtb runs (one core each)")
    a = ap.parse_args(argv)

    todo = [(row["id"], f, row["charge"], m) for row in read_manifest(a.manifest) for f, m in tasks(row)]
    with ThreadPoolExecutor(a.nprocs) as pool:
        status = list(pool.map(lambda t: preopt(t[1], t[2], t[3], a.xtb), todo))
    with open(Path(a.manifest).parent / "preopt_log.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "file", "mult", "status"])
        w.writerows((cid, path.name, m, s) for (cid, path, _, m), s in zip(todo, status))
    failed = [f"{t[1].name}: {s}" for t, s in zip(todo, status) if s.startswith("failed")]
    print(f"{len(todo) - len(failed)}/{len(todo)} structures ok", *failed, sep="\n")


if __name__ == "__main__":
    main()
