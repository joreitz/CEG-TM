"""Tests for datagen/. ORCA output snippets use the assumed output format; confirm it in the pilot."""
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "datagen"))
import build_structures  # noqa: E402
import collect  # noqa: E402
import make_inputs  # noqa: E402
import preopt_xtb  # noqa: E402

UKS_ORBITALS = """
----------------
ORBITAL ENERGIES
----------------
                 SPIN UP ORBITALS
  NO   OCC          E(Eh)            E(eV)
   0   1.0000      -1.000000       -27.2114
   1   0.9000      -0.300000        -8.1634
   2   0.1000      -0.290000        -7.8913
   3   0.0000       0.100000         2.7211

                 SPIN DOWN ORBITALS
  NO   OCC          E(Eh)            E(eV)
   0   1.0000      -0.900000       -24.4903
   1   0.0000       0.200000         5.4423

********************************
"""


def test_occupations_and_nfod_uks():
    chans = collect.occupations(UKS_ORBITALS)
    assert chans == [[1.0, 0.9, 0.1, 0.0], [1.0, 0.0]]
    nfod, nel = collect.nfod_from_occ(chans)
    assert nfod == pytest.approx(0.2) and nel == pytest.approx([2.0, 1.0])


def test_nfod_restricted_and_integer():
    assert collect.nfod_from_occ([[2.0, 2.0, 1.8, 0.2, 0.0]])[0] == pytest.approx(0.4)
    assert collect.nfod_from_occ([[1.0, 1.0, 0.0], [1.0, 0.0]])[0] == 0.0
    assert math.isnan(collect.nfod_from_occ([])[0])


def test_patterns():
    text = ("Total Energy       :        -1263.12345678 Eh          -34371.1 eV\n"
            "Expectation value of <S**2>     :     6.0123\n"
            "N_FOD =     0.938222\n"
            "T1 diagnostic                              ...      0.0234\n"
            "FINAL SINGLE POINT ENERGY     -1264.876543210\n")
    assert collect.last("scf", text) == -1263.12345678
    assert collect.last("s2", text) == 6.0123
    assert collect.last("nfod", text) == 0.938222
    assert collect.last("t1", text) == 0.0234
    assert collect.last("E", text) == -1264.87654321
    assert math.isnan(collect.last("E", ""))


def state(mult, scf, ec, nfod=0.0, t1=0.01, ok=True):
    """Synthetic state: E_c(X) = ec + 0.5 X^-beta; PNO errors 0.001 / 0.0006 Eh at TCutPNO 1e-6 / 1e-7,
    which F = 1.5 extrapolates away exactly (e7 = F / (1 + F) * e6)."""
    tz, qz = ec + 0.5 * 3**-collect.BETA, ec + 0.5 * 4**-collect.BETA
    return {"id": "x", "mult": mult, "E_scf_qz": scf, "Ec_tz6": tz + 0.001, "Ec_tz7": tz + 0.0006,
            "Ec_qz6": qz + 0.001, "E_tpssh": scf, "E_tpss": scf, "E_tpss0": scf, "E_r2scan0": scf,
            "nfod": nfod, "t1": t1, "ok": ok}


def test_e_ref_recovers_cps_and_cbs_limit():
    assert collect.e_ref(state(1, -100.0, -2.0)) == pytest.approx(-102.0, abs=1e-10)


def test_gaps():
    lo, hi = state(1, -100.0, -2.0), state(5, -100.01, -1.99, nfod=0.3, t1=0.06)
    hi.update(E_tpss=-100.02, E_tpssh=-100.01, E_tpss0=-99.995)
    g = collect.gaps([lo, hi])[0]
    assert (g["mult_lo"], g["mult_hi"]) == (1, 5)
    assert g["dE_ref"] == pytest.approx(0.0, abs=1e-8)
    assert g["mr_flag"] and g["sigma_mr"] == collect.SIGMA_CC_MR
    assert g["dnfod"] == pytest.approx(0.3)
    assert g["dE_dax"] == pytest.approx(0.1 * collect.HARTREE)  # gap rises 0.001 Eh per 0.01 a_x
    nev = {("x", 1): -1.0, ("x", 5): -1.0 + 10 / collect.HARTREE}  # NEVPT2 gap 10 kcal/mol
    assert collect.gaps([lo, hi], nev)[0]["sigma_mr"] == pytest.approx(5.0)


def test_collect_end_to_end(tmp_path):
    (tmp_path / "c.csv").write_text("id,metal,ox,charge,dn,row,mults,geometry,ligand_class,flags,xyz\n"
                                    "fe,Fe,2,2,6,3d,1 3,oct,test,,fe.xyz\n")
    rks = "N_FOD =     0.200000\nORBITAL ENERGIES\n  NO OCC\n   0   2.0000  -1.0  -27.2\n" \
          "   1   1.9000  -0.3  -8.1\n   2   0.1000  -0.2  -5.4\n\n*****\n"
    for m, extra in ((1, rks), (3, "N_FOD =     0.200000\n" + UKS_ORBITALS)):
        for job in make_inputs.JOBS:
            (tmp_path / f"runs/fe/m{m}/{job}").mkdir(parents=True)
            (tmp_path / f"runs/fe/m{m}/{job}/{job}.out").write_text(
                f"Total Energy  :  -100.0 Eh\nFINAL SINGLE POINT ENERGY  -102.0\n{extra if job == 'fod' else ''}"
                "****ORCA TERMINATED NORMALLY****\n")
    collect.main([str(tmp_path / "c.csv"), str(tmp_path / "runs"), "--out", str(tmp_path)])
    rows = list(collect.csv.DictReader(open(tmp_path / "states.csv")))
    assert rows[0]["ok"] == "True" and float(rows[0]["nfod_check"]) == pytest.approx(0.2)
    assert rows[1]["ok"] == "False" and "fod:spin" in rows[1]["warn"]  # N_a - N_b = 1, not 2
    assert "tpssh:no_s2" in rows[1]["warn"]


def test_make_inputs_stages(tmp_path):
    (tmp_path / "fe.xyz").write_text("1\n\nFe 0 0 0\n")
    (tmp_path / "c.csv").write_text("id,metal,ox,charge,dn,row,mults,geometry,ligand_class,flags,xyz\n"
                                    "fe,Fe,2,2,6,3d,1 5,oct,test,,fe.xyz\n")
    run = tmp_path / "runs"
    make_inputs.main([str(tmp_path / "c.csv"), str(run), "--stage", "opt"])
    assert "! RKS r2SCAN-3c" in (run / "fe/m1/opt/opt.inp").read_text()
    assert "* xyzfile 2 5 geom.xyz" in (run / "fe/m5/opt/opt.inp").read_text()

    make_inputs.main([str(tmp_path / "c.csv"), str(run), "--stage", "sp"])  # opt not run yet
    assert not (run / "fe/m5/fod").exists()

    opt = run / "fe/m5/opt"
    (opt / "opt.out").write_text("THE OPTIMIZATION HAS CONVERGED\n****ORCA TERMINATED NORMALLY****\n")
    (opt / "opt.xyz").write_text("1\n\nFe 0 0 0\n")
    make_inputs.main([str(tmp_path / "c.csv"), str(run), "--stage", "sp"])
    fod = (run / "fe/m5/fod/fod.inp").read_text()
    assert "! UKS TPSS def2-TZVP" in fod and "SmearTemp 5000" in fod
    assert not (run / "fe/m1/fod").exists()


def test_candidates_consistent():
    rows = [r for _, _, r in build_structures.candidates()]
    assert len({r["id"] for r in rows}) == len(rows)
    assert all(build_structures.QMIN <= r["charge"] <= build_structures.QMAX for r in rows)
    for r in rows:  # electron-count parity: d^n plus ligand radicals/NO must match every multiplicity
        assert len({int(m) % 2 for m in r["mults"].split()}) == 1, r["id"]
    for mtag, stag, r in build_structures.candidates():
        if r["id"] not in build_structures.EXTRA:  # grid: d^n parity fixes the multiplicity parity
            assert all((int(m) - 1) % 2 == r["dn"] % 2 for m in r["mults"].split()), r["id"]
    assert set(build_structures.PILOT) <= {r["id"] for r in rows}


def test_graph_and_tasks(tmp_path):
    el = ["Fe", "Cl", "O", "H"]
    xyz = [(0, 0, 0), (2.3, 0, 0), (0, 2.2, 0), (0, 3.15, 0)]
    assert preopt_xtb.graph(el, xyz) == {(1, 0), (2, 0), (3, 2)}
    xyz[1] = (3.5, 0, 0)  # Fe-Cl dissociated
    assert (1, 0) not in preopt_xtb.graph(el, xyz)
    (tmp_path / "a.xyz").touch()
    (tmp_path / "a_m1.xyz").touch()
    row = {"xyz": tmp_path / "a.xyz", "mults": [1, 3, 5]}
    assert preopt_xtb.tasks(row) == [(tmp_path / "a_m1.xyz", 1), (tmp_path / "a.xyz", 5)]


Z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15, "S": 16, "Cl": 17, "Br": 35, "Cr": 24, "Mn": 25, "Fe": 26,
     "Co": 27, "Ni": 28, "Ru": 44, "Rh": 45}


@pytest.mark.skipif(not (ROOT := Path(__file__).resolve().parents[1] / "structures/complexes.csv").is_file(),
                    reason="structures not built")
def test_structures_electron_parity():
    for row in make_inputs.read_manifest(ROOT):
        for xyz in [row["xyz"], *row["xyz"].parent.glob(f"{row['xyz'].stem}_m*.xyz")]:
            el = preopt_xtb.read_xyz(xyz)[0]
            nel = sum(Z[e] for e in el) - row["charge"]
            assert all((nel - m + 1) % 2 == 0 for m in row["mults"]), xyz.name
