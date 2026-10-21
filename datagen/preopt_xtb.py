#!/usr/bin/env python3
"""Pre-optimise the start structures of a manifest with GFN2-xTB (see docs/data_plan.md).

Each xyz file is optimised for the state it serves: <id>_m<mult>.xyz for that multiplicity, <id>.xyz for
the highest multiplicity without its own file. A result is accepted only if xtb terminated normally and the
bond graph (geom.graph, including metal-donor bonds) is unchanged; otherwise the input structure is kept and
the failure is logged in preopt_log.csv. Optimised files are marked in the comment line and skipped on re-runs.
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

from geom import graph, read_xyz
from make_inputs import read_manifest

TAG = "GFN2-xTB preopt"


def run_xtb(path, charge, mult, xtb):
    """GFN2-xTB optimisation of an xyz file: (elements, coordinates) or an error string."""
    el, xyz, _ = read_xyz(path)
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(path, Path(tmp) / "in.xyz")
        cmd = [xtb, "in.xyz", "--opt", "--gfn", "2", "--chrg", str(charge), "--uhf", str(mult - 1)]
        try:
            p = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, timeout=3600,
                               env={**os.environ, "OMP_NUM_THREADS": "1"})
        except subprocess.TimeoutExpired:
            return "timeout"
        opt = Path(tmp) / "xtbopt.xyz"
        if "normal termination of xtb" not in p.stdout + p.stderr or not opt.is_file():
            return "xtb error"
        el2, xyz2, _ = read_xyz(opt)
    changed = graph(el, xyz) ^ graph(el2, xyz2)
    if changed:
        return "bonding changed " + " ".join(
            f"{el[i]}{i + 1}-{el[j]}{j + 1}:{math.dist(xyz[i], xyz[j]):.2f}->{math.dist(xyz2[i], xyz2[j]):.2f}"
            for i, j in sorted(changed))
    return el2, xyz2


def preopt(path, charge, mults, xtb):
    """Optimise an xyz file in place, trying the multiplicities it serves from high to low (GFN2 often
    over-stretches high-spin bonds); return a status string."""
    comment = read_xyz(path)[2]
    if TAG in comment:
        return f"ok ({comment.removeprefix(TAG).strip()}, earlier run)"
    errors = []
    for mult in mults:
        res = run_xtb(path, charge, mult, xtb)
        if isinstance(res, str):
            errors.append(f"mult {mult}: {res}")
            continue
        body = [f"{e} {x:12.6f} {y:12.6f} {z:12.6f}" for e, (x, y, z) in zip(*res)]
        Path(path).write_text("\n".join([str(len(body)), f"{TAG} charge={charge} mult={mult}", *body]) + "\n")
        return f"ok (mult={mult})"
    return "failed: " + "; ".join(errors)


def tasks(row):
    """(xyz file, multiplicities it serves, highest first) for one manifest row."""
    xyz, mults = row["xyz"], sorted(row["mults"])
    own = {m: xyz.with_name(f"{xyz.stem}_m{m}.xyz") for m in mults}
    out = [(f, [m]) for m, f in own.items() if f.is_file()]
    rest = sorted((m for m in mults if not own[m].is_file()), reverse=True)
    return out + ([(xyz, rest)] if rest else [])


def main(argv=None):
    ap = argparse.ArgumentParser(description="Pre-optimise manifest structures with GFN2-xTB.")
    ap.add_argument("manifest", help="complexes.csv")
    ap.add_argument("--xtb", default="xtb", help="xtb executable")
    ap.add_argument("--nprocs", type=int, default=4, help="parallel xtb runs (one core each)")
    a = ap.parse_args(argv)

    todo = [(row["id"], f, row["charge"], ms) for row in read_manifest(a.manifest) for f, ms in tasks(row)]
    with ThreadPoolExecutor(a.nprocs) as pool:
        status = list(pool.map(lambda t: preopt(t[1], t[2], t[3], a.xtb), todo))
    with open(Path(a.manifest).parent / "preopt_log.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "file", "status"])
        w.writerows((cid, path.name, s) for (cid, path, _, _), s in zip(todo, status))
    failed = [f"{t[1].name}: {s}" for t, s in zip(todo, status) if s.startswith("failed")]
    print(f"{len(todo) - len(failed)}/{len(todo)} structures ok", *failed, sep="\n")


if __name__ == "__main__":
    main()
