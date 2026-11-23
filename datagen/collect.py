#!/usr/bin/env python3
"""Collect ORCA results of the FOD-TM data set into states.csv and gaps.csv (see docs/data_plan.md).

states.csv: one row per complex and spin state, energies in Eh.
gaps.csv:   one row per spin gap relative to the lowest multiplicity of the complex, in kcal/mol.
Missing or suspicious values become NaN plus an entry in 'warn'; nothing is silently defaulted.
"""
import argparse
import csv
import math
import re
from pathlib import Path

from make_inputs import JOBS, NORMAL, read_manifest, stem

HARTREE = 627.509474  # kcal/mol per Eh
F_PNO = 1.5  # CPS extrapolation TCutPNO 1e-6 -> 1e-7 (Altun, Neese, Bistoni, JCTC 2020)
BETA = 2.97  # def2 TZ/QZ correlation-energy exponent (Neese, Valeev, JCTC 2011); check in pilot
SIGMA_CC, SIGMA_CC_MR = 1.5, 3.0  # kcal/mol: CCSD(T) MAE on SSE17, doubled if T1 > T1_MAX
T1_MAX = 0.05  # 3d-metal threshold (Jiang, DeYonker, Wilson, JCTC 2012)
S2_TOL = 0.1
FOD_TOL = 0.01  # covers the 4-decimal rounding of the printed occupation numbers
REQUIRED = ("tpssh", "fod", "cc_tz6", "cc_tz7", "cc_qz6")  # a state is 'ok' only if these are clean
AX = {"tpss": 0.0, "tpssh": 0.10, "tpss0": 0.25}  # Fock-exchange fractions of the a_x scan

PAT = {
    "E": r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)",
    "scf": r"Total Energy\s*:\s*(-?\d+\.\d+)",
    "s2": r"Expectation value of <S\*\*2>\s*:\s*(-?\d+\.\d+)",
    "nfod": r"N_FOD\s*=\s*(-?\d+\.\d+)",
    "t1": r"T1 diagnostic\s*\.*\s*(-?\d+\.\d+)",
}
ORB = re.compile(r"^\s*\d+\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+-?\d+\.\d+\s*$")  # NO OCC E(Eh) E(eV)


def last(key, text):
    """Last match of PAT[key] as float, NaN if absent."""
    m = re.findall(PAT[key], text)
    return float(m[-1]) if m else math.nan


def occupations(text):
    """Occupation numbers per spin channel, sorted by orbital energy, from the last ORBITAL ENERGIES block."""
    i = text.rfind("ORBITAL ENERGIES")
    chans, cur = [], None
    for line in text[i:].splitlines()[1:] if i >= 0 else []:
        m = ORB.match(line)
        if "SPIN" in line and "ORBITALS" in line:
            cur = []
            chans.append(cur)
        elif m:
            if cur is None:
                cur = []
                chans.append(cur)
            cur.append((float(m[2]), float(m[1])))
        elif cur and line.strip() and line.split()[0] != "NO":
            break
    return [[f for _, f in sorted(c)] for c in chans]


def nfod_from_occ(chans):
    """N_FOD = sum_i |delta_i - f_i| (Grimme/Hansen), delta_i = aufbau occupation at T = 0.

    Occupied orbitals contribute 1 - f_i, virtual ones f_i. One channel means restricted (f <= 2).
    Returns (N_FOD, electrons per channel).
    """
    if not chans:
        return math.nan, []
    fmax = 2.0 if len(chans) == 1 else 1.0
    total, nel = 0.0, []
    for occ in chans:
        n = sum(occ)
        k = round(n / fmax)
        total += sum(abs((fmax if i < k else 0.0) - f) for i, f in enumerate(occ))
        nel.append(n)
    return total, nel


def e_ref(r, cps=True, cbs=True):
    """E_SCF(QZ) + E_c(TZ, CPS-extrapolated) + [E_c,CBS(TZ/QZ) - E_c(TZ)] in Eh.

    cps=False / cbs=False give the non-extrapolated variants used for sigma_PNO / sigma_CBS.
    """
    ec6, ec7, q6 = r["Ec_tz6"], r["Ec_tz7"], r["Ec_qz6"]
    ec = ec7 + F_PNO * (ec7 - ec6) if cps else ec7
    ec_basis = (3**BETA * ec6 - 4**BETA * q6) / (3**BETA - 4**BETA) if cbs else q6
    return r["E_scf_qz"] + ec + ec_basis - ec6


def parse_state(folder, cid, charge, mult):
    """States.csv fields (without id/mult) from the outputs of one spin state in the job folder."""
    warn, txt = [], {}
    for job in JOBS:
        if job == "opt":
            continue
        out = folder / f"{stem(cid, job, charge, mult)}.out"
        t = out.read_text(errors="ignore") if out.is_file() else None
        if t is None or NORMAL not in t:
            warn.append(f"{job}:{'missing' if t is None else 'failed'}")
        txt[job] = t if t and NORMAL in t else ""

    def get(job, key):
        v = last(key, txt[job])
        if txt[job] and math.isnan(v):
            warn.append(f"{job}:no_{key}")
        return v

    r = {f"E_{j}": get(j, "E") for j in ("tpssh", "tpss", "tpss0", "r2scan0")}
    s = (mult - 1) / 2
    r["s2"] = get("tpssh", "s2") if mult > 1 else math.nan
    if abs(r["s2"] - s * (s + 1)) > S2_TOL:
        warn.append("tpssh:s2")

    r["nfod"] = get("fod", "nfod")
    chans = occupations(txt["fod"])
    r["nfod_check"], nel = nfod_from_occ(chans)
    if txt["fod"] and not chans:
        warn.append("fod:no_orbitals")
    if abs(r["nfod"] - r["nfod_check"]) > FOD_TOL:
        warn.append("fod:nfod_mismatch")
    if len(nel) == 2 and abs(nel[0] - nel[1] - (mult - 1)) > FOD_TOL:
        warn.append("fod:spin")  # smearing did not conserve N_alpha - N_beta

    for tag in ("tz6", "tz7", "qz6"):
        r[f"Ec_{tag}"] = get(f"cc_{tag}", "E") - get(f"cc_{tag}", "scf")
    r["E_scf_qz"] = last("scf", txt["cc_qz6"])
    t1 = [last("t1", txt[f"cc_{t}"]) for t in ("tz6", "tz7", "qz6")]
    r["t1"] = max((v for v in t1 if not math.isnan(v)), default=math.nan)
    if txt["cc_tz7"] and math.isnan(r["t1"]):
        warn.append("t1:missing")  # informative only, does not affect 'ok'
    r["E_ref"] = e_ref(r)
    r["ok"] = not any(w.split(":")[0] in REQUIRED for w in warn)
    r["warn"] = ";".join(warn)
    return r


def slope(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


def gaps(states, nevpt2=None):
    """Gap rows (kcal/mol) of every state relative to the lowest multiplicity of the same complex."""
    nevpt2, by_id, rows = nevpt2 or {}, {}, []
    for st in states:
        by_id.setdefault(st["id"], []).append(st)
    for cid, sts in by_id.items():
        lo = min(sts, key=lambda st: st["mult"])
        for hi in sts:
            if hi is lo:
                continue

            def d(f):
                return (f(hi) - f(lo)) * HARTREE

            ref = d(e_ref)
            s_pno = abs(ref - d(lambda st: e_ref(st, cps=False)))
            s_cbs = abs(ref - d(lambda st: e_ref(st, cbs=False)))
            mr = any(st["t1"] > T1_MAX for st in (lo, hi))
            s_mr = SIGMA_CC_MR if mr else SIGMA_CC
            keys = [(cid, lo["mult"]), (cid, hi["mult"])]
            if all(k in nevpt2 for k in keys):
                s_mr = max(s_mr, abs(ref - (nevpt2[keys[1]] - nevpt2[keys[0]]) * HARTREE) / 2)
            rows.append({
                "id": cid, "mult_lo": lo["mult"], "mult_hi": hi["mult"],
                "dE_ref": ref, "sigma": math.sqrt(s_pno**2 + s_cbs**2 + s_mr**2),
                "sigma_pno": s_pno, "sigma_cbs": s_cbs, "sigma_mr": s_mr,
                "dE_tpssh": d(lambda st: st["E_tpssh"]), "dE_r2scan0": d(lambda st: st["E_r2scan0"]),
                "dE_dax": slope(list(AX.values()), [d(lambda st, j=j: st[f"E_{j}"]) for j in AX]),
                "dnfod": hi["nfod"] - lo["nfod"], "mr_flag": mr, "ok": lo["ok"] and hi["ok"],
            })
    return rows


def write(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Collect ORCA results into states.csv and gaps.csv.")
    ap.add_argument("manifest", help="complexes.csv")
    ap.add_argument("folder", help="folder with all inputs and outputs")
    ap.add_argument("--nevpt2", help="optional CSV with columns id,mult,E (total energy in Eh)")
    ap.add_argument("--out", default=".", help="output directory")
    a = ap.parse_args(argv)

    states = [{"id": row["id"], "mult": m, **parse_state(Path(a.folder), row["id"], row["charge"], m)}
              for row in read_manifest(a.manifest) for m in row["mults"]]
    nevpt2 = {}
    if a.nevpt2:
        with open(a.nevpt2, newline="") as f:
            nevpt2 = {(r["id"], int(r["mult"])): float(r["E"]) for r in csv.DictReader(f)}
    rows = gaps(states, nevpt2)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    write(out / "states.csv", states)
    write(out / "gaps.csv", rows)
    print(f"{sum(s['ok'] for s in states)}/{len(states)} states ok, "
          f"{sum(g['ok'] for g in rows)}/{len(rows)} gaps ok -> {out}")


if __name__ == "__main__":
    main()
