# CEG-TM
Additive correction term for correlation estimation from geometry-based ligand field splitting for transition metals

- `docs/data_plan.md` – model (FOD-TM), reference data plan and pilot checklist
- `datagen/make_inputs.py` – ORCA inputs per complex, spin state and stage (`opt`, `sp`, `cc`)
- `datagen/collect.py` – parses ORCA outputs into `states.csv` and `gaps.csv`
- `tests/` – `python -m pytest tests/`
