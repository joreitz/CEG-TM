#!/usr/bin/env python3
"""Write ORCA 6 inputs for the FOD-TM reference data set (see docs/data_plan.md).

All inputs go into one folder as <name>_<charge>_<mult>.inp, e.g. fe2-h2o6-tpssh_p2_5.inp:
  name    complex id and job, joined by '-' (no '_' inside the name)
  charge  p<q> for positive, n<|q|> for zero and negative charges (p2, n0, n1)
  mult    spin multiplicity (1 singlet, 2 doublet, ...)
Coordinates are written into the input, so every .inp is self-contained. Outputs are expected next to the
inputs as <same stem>.out. Finished jobs are skipped; jobs whose prerequisites are not finished are reported
and skipped. The DLPNO jobs read the TPSSh orbitals from <tpssh stem>.gbw in the same folder (absolute path).
"""
import argparse
import csv
import sys
from pathlib import Path

from geom import graph, read_xyz

NORMAL = "ORCA TERMINATED NORMALLY"
CONVERGED = "OPTIMIZATION HAS CONVERGED"
COLUMNS = {"id", "metal", "ox", "charge", "dn", "row", "mults", "geometry", "ligand_class", "flags", "xyz"}


def cc(basis, tcut):
    """DLPNO-CCSD(T1) header; guess orbitals are taken from the TPSSh single point."""
    return (f"! {{hf}} DLPNO-CCSD(T1) TightPNO {basis} {basis}/C def2/J RIJCOSX TightSCF MORead\n"
            f'%moinp "{{gbw}}"\n%scf StabPerform true end\n%mdci TCutPNO {tcut} end')


# job -> (stage, ORCA header); {ks}/{hf} become RKS/RHF for singlets and UKS/UHF otherwise.
# Keywords to confirm in the pilot: TPSS0, r2SCAN0, D4 parameters for TPSS/TPSS0, DLPNO-CCSD(T1).
JOBS = {
    "opt": ("opt", "! {ks} r2SCAN-3c Opt TightSCF"),
    "tpssh": ("sp", "! {ks} TPSSh D4 def2-TZVP def2/J RIJCOSX TightSCF\n%scf StabPerform true end"),
    "tpss": ("sp", "! {ks} TPSS D4 def2-TZVP def2/J TightSCF"),
    "tpss0": ("sp", "! {ks} TPSS0 D4 def2-TZVP def2/J RIJCOSX TightSCF"),
    "r2scan0": ("sp", "! {ks} r2SCAN0 D4 def2-TZVP def2/J RIJCOSX TightSCF"),
    "fod": ("sp", "! {ks} TPSS def2-TZVP def2/J TightSCF FOD\n%scf SmearTemp 5000 end"),
    "cc_tz6": ("cc", cc("def2-TZVPP", "1e-6")),
    "cc_tz7": ("cc", cc("def2-TZVPP", "1e-7")),
    "cc_qz6": ("cc", cc("def2-QZVPP", "1e-6")),
}


def read_manifest(path):
    """Manifest rows with 'mults' as list of int, 'charge' as int and 'xyz' relative to the manifest."""
    path = Path(path)
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        missing = COLUMNS - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"{path}: missing columns {sorted(missing)}")
        rows = list(reader)
    ids = [r["id"] for r in rows]
    if len(ids) != len(set(ids)):
        sys.exit(f"{path}: duplicate ids")
    for r in rows:
        r["mults"] = [int(m) for m in r["mults"].split()]
        r["charge"] = int(r["charge"])
        r["xyz"] = path.parent / r["xyz"]
    return rows


def stem(cid, job, charge, mult):
    """File stem <name>_<charge>_<mult>, e.g. stem('fe2_h2o6', 'cc_tz6', 2, 5) = 'fe2-h2o6-cc-tz6_p2_5'."""
    return f"{cid}-{job}".replace("_", "-") + f"_{f'p{charge}' if charge > 0 else f'n{-charge}'}_{mult}"


def start_xyz(row, mult):
    """Start structure of a state: <id>_m<mult>.xyz if present, else <id>.xyz."""
    own = row["xyz"].with_name(f"{row['xyz'].stem}_m{mult}.xyz")
    return own if own.is_file() else row["xyz"]


def done(out):
    return out.is_file() and NORMAL in out.read_text(errors="ignore")


def last_geometry(text):
    """Elements and coordinates of the last 'CARTESIAN COORDINATES (ANGSTROEM)' block of an ORCA output."""
    i = text.rfind("CARTESIAN COORDINATES (ANGSTROEM)")
    el, xyz = [], []
    for line in text[i:].splitlines()[2:] if i >= 0 else []:
        p = line.split()
        if len(p) != 4:
            break
        el.append(p[0])
        xyz.append(tuple(map(float, p[1:])))
    return el, xyz


def geometry(stage, row, mult, folder):
    """(elements, coordinates) for a job, or a reason why the job is not ready yet."""
    el0, xyz0, _ = read_xyz(start_xyz(row, mult))
    if stage == "opt":
        return el0, xyz0
    out = folder / f"{stem(row['id'], 'opt', row['charge'], mult)}.out"
    text = out.read_text(errors="ignore") if out.is_file() else ""
    if NORMAL not in text or CONVERGED not in text:
        return "opt not converged"
    el, xyz = last_geometry(text)
    if el != el0:
        return "no final geometry in opt output"
    if graph(el0, xyz0) != graph(el, xyz):
        return "bonding changed during opt (e.g. ligand dissociated), inspect the opt output"
    if stage == "cc":
        tpssh = folder / stem(row["id"], "tpssh", row["charge"], mult)
        if not (done(tpssh.with_suffix(".out")) and tpssh.with_suffix(".gbw").is_file()):
            return "tpssh .out/.gbw not available"
    return el, xyz


def main(argv=None):
    ap = argparse.ArgumentParser(description="Write ORCA inputs for one stage of the FOD-TM data set.")
    ap.add_argument("manifest", help="complexes.csv")
    ap.add_argument("folder", help="folder for all inputs and outputs")
    ap.add_argument("--stage", required=True, choices=["opt", "sp", "cc"])
    ap.add_argument("--jobs", nargs="+", help="restrict to these jobs of the stage (default: all)")
    ap.add_argument("--nprocs", type=int, default=8)
    ap.add_argument("--maxcore", type=int, default=4000, help="MB per core")
    a = ap.parse_args(argv)

    jobs = [j for j, (stage, _) in JOBS.items() if stage == a.stage and (not a.jobs or j in a.jobs)]
    folder, written, skipped = Path(a.folder).resolve(), [], 0
    folder.mkdir(parents=True, exist_ok=True)
    for row in read_manifest(a.manifest):
        for mult in row["mults"]:
            for job in jobs:
                name = stem(row["id"], job, row["charge"], mult)
                if done(folder / f"{name}.out"):
                    skipped += 1
                    continue
                geo = geometry(a.stage, row, mult, folder)
                if isinstance(geo, str):
                    print(f"not ready: {name} ({geo})", file=sys.stderr)
                    continue
                ks, hf = ("RKS", "RHF") if mult == 1 else ("UKS", "UHF")
                gbw = folder / f"{stem(row['id'], 'tpssh', row['charge'], mult)}.gbw"
                coords = [f"  {e:2s} {x:14.8f} {y:14.8f} {z:14.8f}" for e, (x, y, z) in zip(*geo)]
                (folder / f"{name}.inp").write_text("\n".join([
                    JOBS[job][1].format(ks=ks, hf=hf, gbw=gbw), f"%pal nprocs {a.nprocs} end",
                    f"%maxcore {a.maxcore}", f"* xyz {row['charge']} {mult}", *coords, "*"]) + "\n")
                written.append(f"{name}.inp")
    if written:
        (folder / f"new_{a.stage}.txt").write_text("\n".join(written) + "\n")
    print(f"{a.stage}: {len(written)} inputs written to {folder} (list: new_{a.stage}.txt), "
          f"{skipped} finished jobs skipped")


if __name__ == "__main__":
    main()
