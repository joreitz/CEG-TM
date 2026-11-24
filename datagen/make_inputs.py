#!/usr/bin/env python3
"""Write ORCA 6 inputs for the FOD-TM reference data set (see docs/data_plan.md).

All inputs go into one folder as <name>_<charge>_<mult>.inp, e.g. fe2-h2o6-tpssh_p2_5.inp:
  name    complex id and job, joined by '-' (no '_' inside the name)
  charge  p<q> for positive, n<|q|> for zero and negative charges (p2, n0, n1)
  mult    spin multiplicity (1 singlet, 2 doublet, ...)
Coordinates are written into the input, so every .inp is self-contained. Outputs are expected next to the
inputs as <same stem>.out. Finished jobs are skipped; jobs whose prerequisites are not finished are reported
and skipped. Optimisations that finished without converging are reported; --restart continues them from
their last geometry. Per-state overrides (exclude a state, or accept a bonding change) are read from
state_overrides.csv next to the manifest. The DLPNO jobs read the TPSSh orbitals via %moinp "<tpssh stem>.gbw" (file name only): a submission
script that runs ORCA in a scratch directory must copy that .gbw there (suborca does this for %moinp files).
"""
import argparse
import csv
import sys
from pathlib import Path

from geom import graph, read_xyz

NORMAL = "ORCA TERMINATED NORMALLY"
CONVERGED = "OPTIMIZATION HAS CONVERGED"
OPT_MAXITER = 300  # ORCA's default (~3 N_atoms, at least 50) was too small for 27 of 649 states
OVERRIDES = "state_overrides.csv"  # id,mult,action,reason; action 'exclude' or 'keep' (accept bonding change)
COLUMNS = {"id", "metal", "ox", "charge", "dn", "row", "mults", "geometry", "ligand_class", "flags", "xyz"}


def cc(basis, tcut):
    """DLPNO-CCSD(T1) header; guess orbitals are taken from the TPSSh single point."""
    return (f"! {{hf}} DLPNO-CCSD(T1) TightPNO {basis} {basis}/C def2/J RIJCOSX TightSCF MORead\n"
            f'%moinp "{{gbw}}"\n{{stab}}\n%mdci TCutPNO {tcut} end')


# Stability analysis: open shells may relax to a lower UKS/UHF solution; singlets stay restricted (decision
# 2026-11-24): an RKS/RHF instability is only reported (collect.py flags it as multireference), because the
# broken-symmetry UKS/UHF solution is spin-contaminated and is not the singlet state.
# NoPop: ORCA 6.1.1 aborts in the population analysis after an unstable RKS without UKS restart (missing .scfr).
STAB = {True: "! NoPop\n%scf StabPerform true StabRestartUHFifUnstable false end", False: "%scf StabPerform true end"}

# job -> (stage, ORCA header); {ks}/{hf} become RKS/RHF for singlets and UKS/UHF otherwise; {stab} see STAB.
# Keywords to confirm in the pilot: TPSS0, r2SCAN0, D4 parameters for TPSS/TPSS0, DLPNO-CCSD(T1).
JOBS = {
    "opt": ("opt", "! {ks} r2SCAN-3c Opt TightSCF NOTRAH KDIIS SOSCF\n"  # helps SCF convergence of TM states
                   f"%geom MaxIter {OPT_MAXITER} end"),
    "tpssh": ("sp", "! {ks} TPSSh D4 def2-TZVP def2/J RIJCOSX TightSCF\n{stab}"),
    "tpss": ("sp", "! {ks} TPSS D4 def2-TZVP def2/J TightSCF"),
    "tpss0": ("sp", "! {ks} TPSS0 D4 def2-TZVP def2/J RIJCOSX TightSCF"),
    "r2scan0": ("sp", "! {ks} r2SCAN0 D4 def2-TZVP def2/J RIJCOSX TightSCF"),
    "fod": ("sp", "! {ks} TPSS def2-TZVP def2/J TightSCF FOD\n%scf SmearTemp 5000 end"),
    "cc_tz6": ("cc", cc("def2-TZVPP", "1e-6")),
    "cc_tz7": ("cc", cc("def2-TZVPP", "1e-7")),
    "cc_qz6": ("cc", cc("def2-QZVPP", "1e-6")),
}


def read_overrides(path):
    """{(id, mult): action} from the overrides file, {} if there is none."""
    if not Path(path).is_file():
        return {}
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    bad = [r for r in rows if r["action"] not in ("exclude", "keep")]
    if bad:
        sys.exit(f"{path}: unknown action in {bad[0]}")
    return {(r["id"], int(r["mult"])): r["action"] for r in rows}


def read_manifest(path):
    """Manifest rows with 'mults' as list of int, 'charge' as int and 'xyz' relative to the manifest.

    States marked 'exclude' in state_overrides.csv are removed from 'mults'; a complex with fewer than two
    states left has no spin gap and is dropped. 'keep' holds the states whose bonding change is accepted.
    """
    path = Path(path)
    overrides = read_overrides(path.parent / OVERRIDES)
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        missing = COLUMNS - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"{path}: missing columns {sorted(missing)}")
        rows = list(reader)
    ids = [r["id"] for r in rows]
    if len(ids) != len(set(ids)):
        sys.exit(f"{path}: duplicate ids")
    kept = []
    for r in rows:
        r["mults"] = [int(m) for m in r["mults"].split() if overrides.get((r["id"], int(m))) != "exclude"]
        r["keep"] = {m for m in r["mults"] if overrides.get((r["id"], m)) == "keep"}
        r["charge"] = int(r["charge"])
        r["xyz"] = path.parent / r["xyz"]
        if len(r["mults"]) >= 2:
            kept.append(r)
    return kept


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
    if graph(el0, xyz0) != graph(el, xyz) and mult not in row["keep"]:
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
    ap.add_argument("--restart", action="store_true",
                    help="opt stage: continue optimisations that finished unconverged from their last geometry")
    ap.add_argument("--nprocs", type=int, default=8)
    ap.add_argument("--maxcore", type=int, default=4000, help="MB per core")
    a = ap.parse_args(argv)

    jobs = [j for j, (stage, _) in JOBS.items() if stage == a.stage and (not a.jobs or j in a.jobs)]
    folder, written, skipped, unconverged = Path(a.folder).resolve(), [], 0, []
    folder.mkdir(parents=True, exist_ok=True)
    for row in read_manifest(a.manifest):
        for mult in row["mults"]:
            for job in jobs:
                name = stem(row["id"], job, row["charge"], mult)
                out = folder / f"{name}.out"
                text = out.read_text(errors="ignore") if out.is_file() else ""
                if NORMAL in text and (job != "opt" or CONVERGED in text):
                    skipped += 1
                    continue
                if job == "opt" and NORMAL in text:  # finished, but not converged
                    if not a.restart:
                        unconverged.append(name)
                        continue
                    geo = last_geometry(text)
                    if not geo[0]:
                        print(f"not ready: {name} (no geometry in the unconverged opt output)", file=sys.stderr)
                        continue
                else:
                    geo = geometry(a.stage, row, mult, folder)
                if isinstance(geo, str):
                    print(f"not ready: {name} ({geo})", file=sys.stderr)
                    continue
                ks, hf = ("RKS", "RHF") if mult == 1 else ("UKS", "UHF")
                gbw = f"{stem(row['id'], 'tpssh', row['charge'], mult)}.gbw"
                coords = [f"  {e:2s} {x:14.8f} {y:14.8f} {z:14.8f}" for e, (x, y, z) in zip(*geo)]
                (folder / f"{name}.inp").write_text("\n".join([
                    JOBS[job][1].format(ks=ks, hf=hf, gbw=gbw, stab=STAB[mult == 1]), f"%pal nprocs {a.nprocs} end",
                    f"%maxcore {a.maxcore}", f"* xyz {row['charge']} {mult}", *coords, "*"]) + "\n")
                written.append(f"{name}.inp")
    if written:
        (folder / f"new_{a.stage}.txt").write_text("\n".join(written) + "\n")
    if unconverged:
        print(f"{len(unconverged)} optimisations finished unconverged (use --restart): {' '.join(unconverged)}",
              file=sys.stderr)
    print(f"{a.stage}: {len(written)} inputs written to {folder} (list: new_{a.stage}.txt), "
          f"{skipped} finished jobs skipped")


if __name__ == "__main__":
    main()
