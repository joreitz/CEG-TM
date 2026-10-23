#!/usr/bin/env python3
"""Build start structures of the FOD-TM candidate complexes with molSimplify (see docs/data_plan.md).

Writes <out>/complexes.csv, <out>/pilot.csv and <out>/xyz/<id>.xyz. CN4 complexes that may switch
between square-planar and tetrahedral also get <id>_m<lowest mult>.xyz (square planar), which
make_inputs.py uses for that state. Needs the molSimplify legacy CLI; the structures are versioned,
so this only has to be re-run when the candidate list changes.
"""
import argparse
import csv
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from geom import METALS as METAL_SYMBOLS, graph, read_xyz

MAX_ATOMS = 52
QMIN, QMAX = -1, 3  # gas-phase charge window (multiply charged anions are often unbound)
ROMAN = {1: "I", 2: "II", 3: "III"}
COORD = {"oct": 6, "spy": 5, "tbp": 5, "sqp": 4, "thd": 4}

# tag: (molSimplify core, oxidation state, d^n, row, multiplicities)
METALS = {
    "cr2": ("cr", 2, 4, "3d", "3 5"), "mn3": ("mn", 3, 4, "3d", "3 5"),
    "mn2": ("mn", 2, 5, "3d", "2 4 6"), "fe3": ("fe", 3, 5, "3d", "2 4 6"),
    "fe2": ("fe", 2, 6, "3d", "1 3 5"), "co3": ("co", 3, 6, "3d", "1 3 5"),
    "co2": ("co", 2, 7, "3d", "2 4"), "ni3": ("ni", 3, 7, "3d", "2 4"),
    "ni2": ("ni", 2, 8, "3d", "1 3"), "co1": ("co", 1, 8, "3d", "1 3"),
    "ru2": ("ru", 2, 6, "4d", "1 3 5"), "ru3": ("ru", 3, 5, "4d", "2 4 6"), "rh3": ("rh", 3, 6, "4d", "1 3 5"),
}

# SMILES ligands: name: (SMILES, 1-based connecting atoms, formal charge)
SMILES = {"pdi": ("CC(=N)c1cccc(n1)C(C)=N", [3, 9, 12], 0), "dtl": ("[S-]C=C[S-]", [1, 4], -2),
          "sme": ("C[S-]", [2], -1)}
CHARGE = {"water": 0, "ammonia": 0, "acetonitrile": 0, "en": 0, "pyridine": 0, "imidazole": 0, "bipy": 0,
          "terpy": 0, "ph3": 0, "pme3": 0, "dmpe": 0, "h2s": 0, "carbonyl": 0, "nitroso": 0, "chloride": -1,
          "fluoride": -1, "bromide": -1, "ncs": -1, "cyanide": -1, "acac": -1, "ox": -2, "cat": -2,
          "porphyrin": -2, **{k: v[2] for k, v in SMILES.items()}}

# tag: (geometry, ligands, occupancies, ligand classes); "sqp;thd" = state-dependent CN4 start geometry.
# A repeated ligand name fixes the placement order (molSimplify -ligloc): cis-M(bpy)2X2, trans-MX2(dmpe)2,
# trans-MX2L2 (square planar).
SETS = {
    "h2o6": ("oct", ["water"], [6], "O-donor"),
    "nh3_6": ("oct", ["ammonia"], [6], "N-sp3"),
    "mecn6": ("oct", ["acetonitrile"], [6], "N-sp"),
    "en3": ("oct", ["en"], [3], "N-sp3"),
    "h2o4cl2": ("oct", ["water", "chloride"], [4, 2], "O-donor;halide"),
    "h2o4f2": ("oct", ["water", "fluoride"], [4, 2], "O-donor;halide"),
    "nh3_4cl2": ("oct", ["ammonia", "chloride"], [4, 2], "N-sp3;halide"),
    "en2cl2": ("oct", ["en", "chloride"], [2, 2], "N-sp3;halide"),
    "py4cl2": ("oct", ["pyridine", "chloride"], [4, 2], "N-sp2;halide"),
    "py4ncs2": ("oct", ["pyridine", "ncs"], [4, 2], "N-sp2;N-sp"),
    "imd4cl2": ("oct", ["imidazole", "chloride"], [4, 2], "N-sp2;halide"),
    "bpy2cl2": ("oct", ["bipy", "chloride", "bipy", "chloride"], [1, 1, 1, 1], "N-sp2;halide"),
    "bpy2ncs2": ("oct", ["bipy", "ncs", "bipy", "ncs"], [1, 1, 1, 1], "N-sp2;N-sp"),
    "acac3": ("oct", ["acac"], [3], "O-donor"),
    "nh3_4ox": ("oct", ["ox", "ammonia"], [1, 4], "N-sp3;O-donor"),  # chelate first, else molSimplify fails
    "nh3_4cn2": ("oct", ["ammonia", "cyanide"], [4, 2], "N-sp3;C-donor"),
    "nh3_5co": ("oct", ["ammonia", "carbonyl"], [5, 1], "N-sp3;C-donor"),
    "nh3_4h2s2": ("oct", ["ammonia", "h2s"], [4, 2], "N-sp3;S-donor"),
    "nh3_4sme2": ("oct", ["ammonia", "sme"], [4, 2], "N-sp3;S-donor"),
    "ph3_4cl2": ("oct", ["ph3", "chloride"], [4, 2], "P-donor;halide"),
    "dmpe2cl2": ("oct", ["dmpe", "dmpe", "chloride", "chloride"], [1, 1, 1, 1], "P-donor;halide"),
    "nh3_4cat": ("oct", ["cat", "ammonia"], [1, 4], "N-sp3;non-innocent"),
    "cl2ph3_2": ("sqp;thd", ["chloride", "ph3", "chloride", "ph3"], [1, 1, 1, 1], "halide;P-donor"),
    "cl2pme3_2": ("sqp;thd", ["chloride", "pme3", "chloride", "pme3"], [1, 1, 1, 1], "halide;P-donor"),
    "br2ph3_2": ("sqp;thd", ["bromide", "ph3", "bromide", "ph3"], [1, 1, 1, 1], "halide;P-donor"),
    "cl2py2": ("sqp;thd", ["chloride", "pyridine", "chloride", "pyridine"], [1, 1, 1, 1], "halide;N-sp2"),
    "ncs2py2": ("sqp;thd", ["ncs", "pyridine", "ncs", "pyridine"], [1, 1, 1, 1], "N-sp;N-sp2"),
    "nh3_4": ("sqp;thd", ["ammonia"], [4], "N-sp3"),
    "en2": ("sqp", ["en"], [2], "N-sp3"),  # molSimplify misplaces chelates in tetrahedral starts
    "acac2": ("sqp", ["acac"], [2], "O-donor"),
    "bpycl2": ("sqp", ["bipy", "chloride"], [1, 2], "N-sp2;halide"),
    "nh3_2cn2": ("sqp;thd", ["ammonia", "cyanide", "ammonia", "cyanide"], [1, 1, 1, 1], "N-sp3;C-donor"),
    "co2ph3_2": ("sqp;thd", ["carbonyl", "ph3", "carbonyl", "ph3"], [1, 1, 1, 1], "C-donor;P-donor"),
    "co2pme3_2": ("sqp;thd", ["carbonyl", "pme3", "carbonyl", "pme3"], [1, 1, 1, 1], "C-donor;P-donor"),
    "cl4": ("thd", ["chloride"], [4], "halide"),
    "terpycl2": ("spy", ["terpy", "chloride"], [1, 2], "N-sp2;halide"),
    "pdicl2": ("spy", ["pdi", "chloride"], [1, 2], "non-innocent;halide"),
    "nh3_4cl": ("spy", ["ammonia", "chloride"], [4, 1], "N-sp3;halide"),
    "h2o4cl": ("spy", ["water", "chloride"], [4, 1], "O-donor;halide"),
    "nh3_3cl2": ("tbp", ["ammonia", "chloride"], [3, 2], "N-sp3;halide"),
    "pme3_3cl2": ("tbp", ["pme3", "chloride"], [3, 2], "P-donor;halide"),
    "co3pme3_2": ("tbp", ["carbonyl", "pme3"], [3, 2], "C-donor;P-donor"),
    "co3ph3_2": ("tbp", ["carbonyl", "ph3"], [3, 2], "C-donor;P-donor"),
    "por": ("sqp", ["porphyrin"], [1], "non-innocent"),
    "porcl": ("spy", ["porphyrin", "chloride"], [1, 1], "non-innocent;halide"),
    "porimd": ("spy", ["porphyrin", "imidazole"], [1, 1], "non-innocent;N-sp2"),
    "porco": ("spy", ["porphyrin", "carbonyl"], [1, 1], "non-innocent;C-donor"),
    "pornh3_2": ("oct", ["porphyrin", "ammonia"], [1, 2], "non-innocent;N-sp3"),
    "porno": ("spy", ["porphyrin", "nitroso"], [1, 1], "non-innocent"),
    "h2o5no": ("oct", ["water", "nitroso"], [5, 1], "O-donor;non-innocent"),
    "nh3_5no": ("oct", ["ammonia", "nitroso"], [5, 1], "N-sp3;non-innocent"),
    "dtl2": ("sqp", ["dtl"], [2], "non-innocent"),
    "bpycat2": ("oct", ["bipy", "cat"], [1, 2], "N-sp2;non-innocent"),
}

OCT_A = ["h2o6", "nh3_6", "mecn6", "en3", "h2o4cl2", "nh3_4cl2", "en2cl2", "py4cl2", "bpy2cl2", "acac3",
         "nh3_4cn2", "dmpe2cl2"]
OCT_B = OCT_A + ["py4ncs2", "imd4cl2", "bpy2ncs2"]
GRID = {
    "cr2": OCT_A + ["nh3_4", "acac2", "cl2py2", "en2", "bpycl2", "nh3_4cl", "terpycl2", "pme3_3cl2"],
    "mn3": [s for s in OCT_A if s != "mecn6"] + ["h2o4f2", "nh3_4ox", "nh3_4cat", "h2o4cl", "nh3_4cl", "terpycl2",
                                                 "porcl"],
    "mn2": OCT_B + ["nh3_4h2s2", "nh3_4sme2", "ph3_4cl2", "cl2ph3_2", "cl2pme3_2", "br2ph3_2", "cl2py2",
                    "ncs2py2", "acac2", "terpycl2", "pme3_3cl2", "nh3_4cl", "nh3_3cl2", "pdicl2"],
    "fe3": [s for s in OCT_B if s != "mecn6"] + ["h2o4f2", "nh3_4ox", "nh3_4cat", "cl4", "terpycl2", "nh3_4cl",
                                                 "h2o4cl", "nh3_3cl2", "nh3_4h2s2", "nh3_4sme2", "porcl",
                                                 "pornh3_2"],
    "fe2": OCT_B + ["h2o4f2", "nh3_4ox", "nh3_5co", "nh3_4h2s2", "nh3_4sme2", "ph3_4cl2", "cl2ph3_2",
                    "cl2pme3_2", "br2ph3_2", "cl2py2", "ncs2py2", "acac2", "terpycl2", "nh3_4cl", "nh3_3cl2",
                    "pme3_3cl2", "pdicl2", "por", "porimd", "porco", "pornh3_2"],
    "co3": [s for s in OCT_B if s != "py4ncs2"] + ["h2o4f2", "nh3_4ox", "nh3_4cat", "nh3_4h2s2", "nh3_4sme2",
                                                   "porcl", "pornh3_2"],
    "co2": OCT_B + ["nh3_4h2s2", "nh3_4sme2", "cl2ph3_2", "cl2pme3_2", "br2ph3_2", "cl2py2", "ncs2py2", "nh3_4",
                    "en2", "acac2", "bpycl2", "nh3_2cn2", "terpycl2", "pme3_3cl2", "nh3_4cl", "nh3_3cl2", "pdicl2",
                    "por"],
    "ni3": ["nh3_6", "en3", "h2o4cl2", "nh3_4cl2", "en2cl2", "bpy2cl2", "dmpe2cl2", "terpycl2", "pme3_3cl2"],
    "ni2": ["h2o6", "nh3_6", "en3", "h2o4cl2", "nh3_4cl2", "nh3_4h2s2", "nh3_4sme2", "en2cl2", "py4cl2", "bpy2cl2",
            "cl2ph3_2", "cl2pme3_2",
            "br2ph3_2", "cl2py2", "ncs2py2", "nh3_4", "en2", "acac2", "bpycl2", "nh3_2cn2", "terpycl2",
            "pme3_3cl2", "nh3_4cl", "nh3_3cl2", "por"],
    "co1": ["co2ph3_2", "co2pme3_2", "co3pme3_2", "co3ph3_2"],
    "ru2": ["h2o6", "nh3_6", "mecn6", "en3", "nh3_4cl2", "en2cl2", "bpy2cl2", "bpy2ncs2", "nh3_4cn2", "nh3_5co",
            "dmpe2cl2"],
    "ru3": ["h2o6", "nh3_6", "en3", "nh3_4cl2", "en2cl2", "bpy2cl2", "acac3", "dmpe2cl2"],
    "rh3": ["h2o6", "nh3_6", "en3", "nh3_4cl2", "en2cl2", "bpy2cl2", "acac3"],
}
# Ligand radicals / NO: charge, oxidation state, d^n and multiplicities do not follow the grid rule.
# NO is counted as neutral, so the formal oxidation state is that of the metal without NO.
EXTRA = {  # id: (metal tag, set, charge, ox, d^n, multiplicities, flags)
    "fe2_porno": ("fe2", "porno", 0, 2, 6, "2 4", "noninnocent"),
    "fe2_h2o5no": ("fe2", "h2o5no", 2, 2, 6, "2 4", "noninnocent"),
    "fe2_nh3_5no": ("fe2", "nh3_5no", 2, 2, 6, "2 4", "noninnocent"),
    "co2_nh3_5no": ("co2", "nh3_5no", 2, 2, 7, "1 3", "noninnocent"),
    "ru3_nh3_5no": ("ru3", "nh3_5no", 3, 3, 5, "1 3", "noninnocent"),
    "ni2_dtl2": ("ni2", "dtl2", 0, 2, 8, "1 3", "noninnocent;mr"),  # bis(dithiolene radical)
    "ni2_dtl2m": ("ni2", "dtl2", -1, 2, 8, "2 4", "noninnocent;mr"),
    "co3_bpycat2": ("co3", "bpycat2", 0, 3, 6, "2 4 6", "noninnocent;mr"),  # Co(III)(cat)(sq) <-> Co(II)(sq)2
}
PILOT = ["fe2_h2o6", "fe2_mecn6", "mn2_h2o6", "ni2_cl2ph3_2", "fe2_por"]


def candidates():
    """(metal tag, set tag, manifest row without xyz) for the grid and the extras, charge window applied."""
    rows = []
    for mtag, sets in GRID.items():
        core, ox, dn, row, mults = METALS[mtag]
        for stag in sets:
            geo, ligs, occ, lclass = SETS[stag]
            q = ox + sum(CHARGE[lig] * n for lig, n in zip(ligs, occ))
            flags = ["noninnocent"] if "non-innocent" in lclass else []
            if geo == "oct" and len(ligs) == 1 and dn in (4, 7):
                flags.append("jt")  # e_g^1 in one of the states
            if mtag == "ni2" and geo == "oct":
                flags.append("mr")  # octahedral d8 singlet is open shell
            if QMIN <= q <= QMAX:
                rows.append((mtag, stag, dict(id=f"{mtag}_{stag}", metal=core.capitalize(), ox=ox, charge=q, dn=dn,
                                              row=row, mults=mults, geometry=geo, ligand_class=lclass,
                                              flags=";".join(flags))))
    for cid, (mtag, stag, q, ox, dn, mults, flags) in EXTRA.items():
        core, row = METALS[mtag][0], METALS[mtag][3]
        rows.append((mtag, stag, dict(id=cid, metal=core.capitalize(), ox=ox, charge=q, dn=dn, row=row,
                                      mults=mults, geometry=SETS[stag][0], ligand_class=SETS[stag][3], flags=flags)))
    return rows


def build(molsimplify, mtag, stag, row, geo, mult, dest):
    """Run molSimplify for one geometry; return an error string or None."""
    core, ox = METALS[mtag][0], row["ox"]
    _, ligs, occ, _ = SETS[stag]
    names = [SMILES[lig][0] if lig in SMILES else lig for lig in ligs]
    smicat = [SMILES[lig][1] for lig in ligs if lig in SMILES]
    with tempfile.TemporaryDirectory() as tmp:
        cmd = [molsimplify, "legacy", "-core", core, "-oxstate", ROMAN[ox], "-geometry", geo,
               "-coord", str(COORD[geo]), "-lig", *names, "-ligocc", *map(str, occ), "-spin", str(mult),
               "-skipANN", "True", "-rundir", tmp]
        if smicat:
            cmd += ["-smicat", str(smicat).replace(" ", "")]
        if len(set(ligs)) < len(ligs):
            cmd += ["-ligloc", "True"]
        try:
            subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=tmp)
        except subprocess.TimeoutExpired:
            return "timeout"
        xyz = sorted(Path(tmp).rglob("*.xyz"))
        if not xyz:
            return "no structure"
        el, coords, _ = read_xyz(xyz[0])
        if len(el) > MAX_ATOMS:
            return f"{len(el)} atoms"
        metal = next(i for i, e in enumerate(el) if e in METAL_SYMBOLS)
        cn = sum(metal in bond for bond in graph(el, coords))
        if cn != COORD[geo]:
            return f"{geo}: {cn} metal contacts instead of {COORD[geo]}"
        tc = sorted(Path(tmp).rglob("terachem_input"))
        m = re.search(r"^charge\s+(-?\d+)", tc[0].read_text(), re.M) if tc else None
        if row["id"] not in EXTRA and (m is None or int(m[1]) != row["charge"]):
            return f"charge mismatch (molSimplify {m[1] if m else '?'}, manifest {row['charge']})"
        shutil.copy(xyz[0], dest)
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description="Build FOD-TM candidate structures with molSimplify.")
    ap.add_argument("out", help="output directory, e.g. structures")
    ap.add_argument("--molsimplify", default="molsimplify", help="molSimplify executable")
    ap.add_argument("--only", nargs="+", help="build only these ids")
    ap.add_argument("--nprocs", type=int, default=4, help="parallel molSimplify runs")
    a = ap.parse_args(argv)

    out = Path(a.out)
    (out / "xyz").mkdir(parents=True, exist_ok=True)

    def build_all(item):
        mtag, stag, row = item
        for old in (out / "xyz").glob(f"{row['id']}.xyz"):
            old.unlink()
        for old in (out / "xyz").glob(f"{row['id']}_m*.xyz"):
            old.unlink()
        mults = [int(m) for m in row["mults"].split()]
        jobs = [(row["geometry"], max(mults), out / "xyz" / f"{row['id']}.xyz")]
        if row["geometry"] == "sqp;thd":  # d4: square planar; otherwise square planar only for the lowest mult
            jobs = [("sqp", max(mults), jobs[0][2])] if row["dn"] == 4 else \
                [("thd", max(mults), jobs[0][2]), ("sqp", min(mults), out / "xyz" / f"{row['id']}_m{min(mults)}.xyz")]
        errors = [e for geo, mult, dest in jobs if (e := build(a.molsimplify, mtag, stag, row, geo, mult, dest))]
        print(f"{row['id']}: {'; '.join(errors) or 'ok'}", flush=True)
        return row, errors

    todo = [c for c in candidates() if not a.only or c[2]["id"] in a.only]
    with ThreadPoolExecutor(a.nprocs) as pool:
        results = list(pool.map(build_all, todo))
    failed = {row["id"]: errors for row, errors in results if errors}
    for cid in failed:  # never leave a partial set of files behind
        for old in (out / "xyz").glob(f"{cid}*.xyz"):
            if old.stem == cid or old.stem.startswith(f"{cid}_m"):
                old.unlink()
    # the manifest lists every candidate with a structure on disk, so partial rebuilds (--only) keep it complete
    rows = [{**row, "xyz": f"xyz/{row['id']}.xyz"} for _, _, row in candidates()
            if (out / "xyz" / f"{row['id']}.xyz").is_file()]
    for name, subset in (("complexes.csv", rows), ("pilot.csv", [r for r in rows if r["id"] in PILOT])):
        if subset:
            with open(out / name, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(subset[0]))
                w.writeheader()
                w.writerows(subset)
    print(f"{len(rows)} in manifest, {len(failed)} failed in this run",
          *(f"{cid}: {'; '.join(e)}" for cid, e in failed.items()), sep="\n")


if __name__ == "__main__":
    main()
