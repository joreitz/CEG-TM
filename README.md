# CEG-TM
Additive correction term for correlation estimation from geometry-based ligand field splitting for transition metals

- `docs/data_plan.md` – model (FOD-TM), reference data plan and pilot checklist
- `structures/` – candidate complexes: `complexes.csv` (manifest), `pilot.csv` (5-complex pilot), `xyz/` (GFN2-xTB pre-optimised start structures)
- `datagen/build_structures.py` – builds the candidate structures with molSimplify
- `datagen/preopt_xtb.py` – GFN2-xTB pre-optimisation with a bond-graph check
- `datagen/build_metallocenes.py` – start structures of the metallocenes M(Cp)2 (M = Cr–Ni)
- `datagen/idealise_oxalate.py` – rebuilds the oxalate ring of the M(ox)(NH3)4 start structures (xTB breaks it)
- `structures/state_overrides.csv` – states excluded from the data set or with an accepted bonding change, each with its reason
- `datagen/make_inputs.py` – ORCA inputs per complex, spin state and stage (`opt`, `sp`, `cc`), all in one folder
- `datagen/collect.py` – parses ORCA outputs into `states.csv` and `gaps.csv`
- `tests/` – `python -m pytest tests/`

Quick start (pilot):

```bash
python datagen/make_inputs.py structures/pilot.csv jobs --stage opt
```

All inputs land in one folder as `<name>_<charge>_<mult>.inp` (e.g. `fe2-h2o6-tpssh_p2_5.inp`); see `docs/data_plan.md`, section 8.
