#!/usr/bin/env python3
"""Write ORCA 6 inputs for the FOD-TM reference data set (see docs/data_plan.md).

Layout: <root>/<id>/m<mult>/<job>/<job>.inp. The geometry is copied in as geom.xyz and guess
orbitals as guess.gbw, so every job directory is self-contained. Run each job in its own
directory as `orca <job>.inp > <job>.out`. Finished jobs are skipped; jobs whose prerequisites
are not finished yet are reported and skipped.
"""
import argparse
import csv
import shutil
import sys
from pathlib import Path

from geom import same_bonding

NORMAL = "ORCA TERMINATED NORMALLY"
CONVERGED = "OPTIMIZATION HAS CONVERGED"
COLUMNS = {"id", "metal", "ox", "charge", "dn", "row", "mults", "geometry", "ligand_class", "flags", "xyz"}


def cc(basis, tcut):
    """DLPNO-CCSD(T1) header; guess orbitals are taken from the TPSSh single point."""
    return (f"! {{hf}} DLPNO-CCSD(T1) TightPNO {basis} {basis}/C def2/J RIJCOSX TightSCF MORead\n"
            f'%moinp "guess.gbw"\n%scf StabPerform true end\n%mdci TCutPNO {tcut} end')


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


def done(out):
    return out.is_file() and NORMAL in out.read_text(errors="ignore")


def sources(stage, sdir, xyz, mult):
    """Files to copy into a job directory as {name: source}, or a reason why the job is not ready."""
    if stage == "opt":
        per_state = xyz.with_name(f"{xyz.stem}_m{mult}.xyz")  # optional state-specific start geometry
        files = {"geom.xyz": per_state if per_state.is_file() else xyz}
    else:
        opt = sdir / "opt" / "opt.out"
        if not (done(opt) and CONVERGED in opt.read_text(errors="ignore")):
            return "opt not converged"
        start, final = sdir / "opt" / "geom.xyz", sdir / "opt" / "opt.xyz"
        if start.is_file() and final.is_file() and not same_bonding(start, final):
            return "bonding changed during opt (e.g. ligand dissociated), inspect opt.xyz"
        files = {"geom.xyz": final}
        if stage == "cc":
            if not done(sdir / "tpssh" / "tpssh.out"):
                return "tpssh not finished"
            files["guess.gbw"] = sdir / "tpssh" / "tpssh.gbw"
    missing = [str(p) for p in files.values() if not p.is_file()]
    return f"missing {', '.join(missing)}" if missing else files


def main(argv=None):
    ap = argparse.ArgumentParser(description="Write ORCA inputs for one stage of the FOD-TM data set.")
    ap.add_argument("manifest", help="complexes.csv")
    ap.add_argument("root", help="run directory")
    ap.add_argument("--stage", required=True, choices=["opt", "sp", "cc"])
    ap.add_argument("--jobs", nargs="+", help="restrict to these jobs of the stage (default: all)")
    ap.add_argument("--nprocs", type=int, default=8)
    ap.add_argument("--maxcore", type=int, default=4000, help="MB per core")
    a = ap.parse_args(argv)

    jobs = [j for j, (stage, _) in JOBS.items() if stage == a.stage and (not a.jobs or j in a.jobs)]
    root, written, skipped = Path(a.root), [], 0
    for row in read_manifest(a.manifest):
        for mult in row["mults"]:
            sdir = root / row["id"] / f"m{mult}"
            for job in jobs:
                jdir = sdir / job
                if done(jdir / f"{job}.out"):
                    skipped += 1
                    continue
                files = sources(a.stage, sdir, row["xyz"], mult)
                if isinstance(files, str):
                    print(f"not ready: {jdir} ({files})", file=sys.stderr)
                    continue
                jdir.mkdir(parents=True, exist_ok=True)
                for name, src in files.items():
                    shutil.copy(src, jdir / name)
                ks, hf = ("RKS", "RHF") if mult == 1 else ("UKS", "UHF")
                (jdir / f"{job}.inp").write_text(
                    JOBS[job][1].format(ks=ks, hf=hf)
                    + f"\n%pal nprocs {a.nprocs} end\n%maxcore {a.maxcore}\n"
                    + f"* xyzfile {row['charge']} {mult} geom.xyz\n")
                written.append(str(jdir))
    if written:
        root.mkdir(parents=True, exist_ok=True)
        (root / f"jobs_{a.stage}.txt").write_text("\n".join(written) + "\n")
    print(f"{a.stage}: {len(written)} inputs written, {skipped} finished jobs skipped")


if __name__ == "__main__":
    main()
