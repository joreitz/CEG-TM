# FOD-TM: reference data plan (spin gaps, TPSSh)

Status: plan for data generation. The fit itself (`fodtm` package) comes later.

## 1. Model

Final energy: $E = E_\mathrm{DFT} + \Delta E$, evaluated separately for every electronic state.

$$\Delta E = \sum_M \Big[\underbrace{s_1(f)\,g(d^n,\mathrm{CN}_M)\,N_{\rm FOD}^{M}}_{\text{(I) static correlation at the metal}} + \underbrace{s_2(f)\,S_M(S_M+1)\,h(d^n)}_{\text{(II) spin / exchange bias}}\Big] + \underbrace{s_3(f)\sum_{L\in\text{1st sphere}} N_{\rm FOD}^{L}}_{\text{(III) ligand MR character}}$$

- $s_1, s_2, s_3$ are global parameters, one set per functional $f$.
- $N_{\rm FOD}^{M}$ and $N_{\rm FOD}^{L}$ are the FOD integrated over the Hirshfeld region of the metal and of the first-sphere donor atoms.
- $g(d^n,\mathrm{CN})$ is either a lookup table over $d^n$ with a linear CN term, or $g = 1 + a(\mathrm{CN}-6) + b(d^n-5)^2$.
  - CN is the D3 coordination number: $\mathrm{CN}_A=\sum_B [1+e^{-k_1(k_2(R_{A,\rm cov}+R_{B,\rm cov})/R_{AB}-1)}]^{-1}$, with $k_1=16$ and $k_2=4/3$.
- $h(d^n)$ is a lookup table over $d^n$, kept separate for 3d and 4d.
- $S_M$ is the nominal total spin for mononuclear complexes.

Expected signs:
- $s_1<0$.
- $s_2<0$ for (meta-)GGAs and $s_2>0$ for hybrids with a large Fock-exchange fraction.

**FOD definition (Grimme/Hansen):**

$$\rho_{\rm FOD}(\mathbf r)=\sum_i(\delta_1-\delta_2 f_i)\,|\phi_i(\mathbf r)|^2$$

- Occupied orbitals at $T=0$ contribute $(1-f_i)$, virtual orbitals contribute $f_i$, so $N_{\rm FOD}=\sum_i|\delta_i-f_i|$.
- The earlier form $\sum_i(\delta_i-f_i)$ with $\delta_i=0$ for virtuals integrates to exactly zero, so it is not used.

Open-shell states use UKS, so singly occupied orbitals of the majority spin do not count as FOD.

## 2. Decisions

| Item | Choice |
|---|---|
| FOD descriptor | TPSS/def2-TZVP, $T_\mathrm{el}$ = 5000 K, for **every** functional. This replaces the $a_x$-dependent temperature rule, so descriptors are functional independent. |
| Functionals to fit | TPSSh first, r2SCAN0 second |
| Target | adiabatic spin gaps (relative to the lowest multiplicity of each complex) |
| Size | ~250 mononuclear complexes, ~600 states |
| Software | ORCA 6 (all automated jobs); ORCA or OpenMolcas/autoCAS for the NEVPT2 subset; xtb optional for pre-optimisation |

## 3. Corrections to the original specification

1. **FOD sign:** see section 1.
2. **Only relative energies are fitted.** Absolute $E_\mathrm{ref}-E_\mathrm{DFT}$ contains method and basis offsets of the order of hartrees.
   - Each fit row is a gap within one complex: $y=\Delta E_\mathrm{ref}-\Delta E_\mathrm{DFT}$.
   - The descriptors enter as state differences.
3. **Identifiability.** With only two states of a $d^n$ class, $\Delta[S(S+1)]$ is constant, and $s_2h(d^n)$ becomes a per-class offset.
   - $s_1g$ is then determined only through the variation of $\Delta N_\mathrm{FOD}$ within the class.
   - Every class therefore needs ligands along the spectrochemical series.
   - Where possible, a class also needs intermediate-spin states, placed in low symmetry so that they are not orbitally degenerate.
4. **$d^n$ classes.**
   - $d^1$, $d^9$ and $d^{10}$ have one spin state only and give no gaps.
   - $d^3$ is left out for now: low-spin Cr(III) (²E) is a spin-flip state with multireference character, not a competing ground state.
   - Classes used: $d^4$–$d^8$ (3d) and $d^5$/$d^6$ (4d).
5. **Charge.** Multiply charged anions are often unbound in the gas phase (e.g. [Fe(CN)6]⁴⁻).
   - Rule: $q \ge -1$ and $|q| \le 3$.
   - Balance charge with heteroleptic ligand sets.
6. **AFQMC/DMC** is not available with the installed software. Hard cases get a larger σ or are dropped.
7. **DLPNO for spin states:** use iterative triples (T1) and extrapolate to the complete PNO space (Drosou, Mitsopoulou, Pantazis 2022). Semicanonical (T0) is a known error source.

## 4. Data set

Candidate set as built (`structures/complexes.csv`): 255 complexes, 649 states.

| $d^n$ | metals (complexes) | states (2S+1) | complexes | states |
|---|---|---|---|---|
| d4 | Cr(II) 20, Mn(III) 18 | 3, 5 | 38 | 76 |
| d5 | Mn(II) 29, Fe(III) 26 | 2, 4, 6 | 55 | 165 |
| d6 | Fe(II) 39, Co(III) 22 | 1, 3, 5 | 61 | 180 |
| d7 | Co(II) 34, Ni(III) 9 | 2, 4 | 43 | 86 |
| d8 | Ni(II) 27, Co(I) 4 | 1, 3 | 31 | 62 |
| 4d d5 | Ru(III) 9 | 2, 4, 6 | 9 | 26 |
| 4d d6 | Ru(II) 11, Rh(III) 7 | 1, 3, 5 | 18 | 54 |

- **Non-innocent ligands (25, flag `noninnocent`):** porphine (Fe, Mn, Co, Ni; with axial Cl/imidazole/CO/NH3/NO), {FeNO}⁷, {CoNO}⁸, {RuNO}, bis(imino)pyridine (Mn, Fe, Co), catecholate (Mn, Fe, Co), Ni bis(dithiolene) (0/−1), Co(bpy)(cat)(sq).
  - NO complexes deviate from the grid's state rule: Fe(II)+NO uses 2S+1 = 2, 4; Co/Ru+NO use 1, 3.
- **Geometry mix:** 171 CN6, 39 CN5 (square pyramidal/TBP), 45 CN4 (square planar/tetrahedral).
- **Ligand classes** (used for leave-ligand-class-out CV, number of complexes): halide 129, N-sp3 103, N-sp2 57, O-donor 45, P-donor 34, non-innocent 25, N-sp 20, C-donor (CN⁻, CO) 17, S-donor (H2S, MeS⁻) 12.
- **Size:** at most 52 atoms, because of the QZ DLPNO cost. Recalibrate after the pilot.
- **Sources:** combinatorial structures (molSimplify, see "Candidate structures" below). Real complexes are added by hand.
  - tmQM does not fit here: it contains only closed-shell complexes with an even electron count, so no d5/d7, and it has no oxidation states.
- **Start geometries per state:** if `xyz/<id>_m<mult>.xyz` exists, the opt stage uses it. Example: square-planar singlet vs. tetrahedral triplet for NiCl2(PH3)2.
- **Flags** (`flags` column, `;`-separated): `noninnocent`, `jt` (homoleptic octahedral d4/d7, e_g¹ in one state), `mr` (expected multireference: octahedral Ni(II) singlet, dithiolene, catecholate/semiquinone). Flagged states are kept and get priority for NEVPT2.
- **Negative controls** (N_FOD ≈ 0 in all states) are identified after the FOD runs, not in advance.
- **External validation (never trained on):** SSE17 (Radoń et al. 2024).
  - Its references are electronic spin-state energetics back-corrected from experiment.
  - Vertical entries need $E$(state j at the geometry of state i).
  - Add a scalar-relativistic (DKH2) DFT increment for this comparison only. In training, DFT and CC are both non-relativistic (def2-ECP for 4d), so the effect cancels.

## 5. Protocol per complex and spin state

RKS/RHF for singlets, UKS/UHF otherwise. The exact ORCA headers live in `datagen/make_inputs.py` (`JOBS`).

| stage | job | method | purpose |
|---|---|---|---|
| opt | `opt` | r2SCAN-3c | geometry of each state |
| sp | `tpssh` | TPSSh-D4/def2-TZVP, StabPerform | $E_\mathrm{DFT}$ (fit target), ⟨S²⟩ |
| sp | `tpss`, `tpss0` | TPSS-D4, TPSS0-D4 /def2-TZVP | $a_x$ = 0 / 0.10 / 0.25 → $d\Delta E/da_x$ |
| sp | `r2scan0` | r2SCAN0-D4/def2-TZVP | second functional |
| sp | `fod` | TPSS/def2-TZVP, SmearTemp 5000 K, `FOD` | $N_\mathrm{FOD}$; keep `fod.gbw` and `fod.scfp_fod` for the Hirshfeld partition |
| cc | `cc_tz6`, `cc_tz7`, `cc_qz6` | DLPNO-CCSD(T1) TightPNO; TCutPNO 1e-6/1e-7 (def2-TZVPP) and 1e-6 (def2-QZVPP); TPSSh orbitals as guess | reference |
| manual | NEVPT2 | CASSCF/SC-NEVPT2, double-shell active space | consensus for the MR subset |

**NEVPT2 subset (~15–20 %):** states with T1 > 0.05 (3d threshold of Jiang, DeYonker, Wilson 2012), high $N_\mathrm{FOD}$, or the flags `jt`/`mr`. Results go into `nevpt2.csv` (`id,mult,E` in Eh), independent of the program used.

The $a_x$ scan checks the fit: if the fitted term (I) correlates strongly with $d\Delta E/da_x$, term (II) is partly redundant.

## 6. Reference energy and uncertainty

The correlation energy of each DLPNO job is $E_c = E_\mathrm{final} - E_\mathrm{SCF}$.

- $E_c^\mathrm{CPS} = E_c(10^{-7}) + F\,[E_c(10^{-7}) - E_c(10^{-6})]$, with $F = 1.5$ (Altun, Neese, Bistoni 2020)
- $E_{c,\mathrm{CBS}} = (3^\beta E_c^{TZ} - 4^\beta E_c^{QZ}) / (3^\beta - 4^\beta)$, with $\beta = 2.97$ (def2 3/4), both at TCutPNO $10^{-6}$
- $E_\mathrm{ref} = E_\mathrm{SCF}(QZ) + E_c^\mathrm{CPS}(TZ) + [E_{c,\mathrm{CBS}} - E_c^{TZ}]$

Per gap, in kcal/mol:

- $\sigma_\mathrm{PNO} = |\Delta E_\mathrm{ref} - \Delta E_\mathrm{ref}(\text{no CPS})|$
- $\sigma_\mathrm{CBS} = |\Delta E_\mathrm{ref} - \Delta E_\mathrm{ref}(\text{QZ instead of CBS})|$
- $\sigma_\mathrm{MR}$:
  - 1.5 (the CCSD(T) MAE on SSE17), or 3.0 if T1 > 0.05 in either state;
  - if NEVPT2 is available: $\max(\sigma_\mathrm{MR}, |\Delta_\mathrm{CC}-\Delta_\mathrm{NEVPT2}|/2)$.
- $\sigma = \sqrt{\sigma_\mathrm{PNO}^2 + \sigma_\mathrm{CBS}^2 + \sigma_\mathrm{MR}^2}$
- Fit weight $w = 1/\sigma^2$.

$F$ and $\beta$ are set in `datagen/collect.py`. Check them in the pilot against ORCA's own `Extrapolate(3/4,def2)` and the `extrapolate_PNO` compound script.

## 7. Files

`complexes.csv` (input, one row per complex):
`id, metal, ox, charge, dn, row, mults, geometry, ligand_class, flags, xyz`
(`mults` is space-separated, e.g. `1 3 5`).

`states.csv` (output, one row per state, energies in Eh):

| column | meaning |
|---|---|
| `E_tpssh`, `E_tpss`, `E_tpss0`, `E_r2scan0` | DFT energies |
| `s2` | ⟨S²⟩ (TPSSh, UKS only) |
| `nfod` | $N_\mathrm{FOD}$ printed by ORCA |
| `nfod_check` | $N_\mathrm{FOD}$ recomputed from the printed occupations |
| `Ec_tz6`, `Ec_tz7`, `Ec_qz6`, `E_scf_qz` | DLPNO components |
| `E_ref` | reference energy |
| `t1` | max T1 diagnostic |
| `ok` | true if all required jobs (`tpssh`, `fod`, `cc_*`) are clean |
| `warn` | all warnings |

`gaps.csv` (output, one row per gap, kcal/mol): `id, mult_lo, mult_hi, dE_ref, sigma, sigma_pno, sigma_cbs, sigma_mr, dE_tpssh, dE_r2scan0, dE_dax, dnfod, mr_flag, ok`.

Atom-resolved descriptors ($N_\mathrm{FOD}^M$, $N_\mathrm{FOD}^L$, CN) will be produced later by `fodtm descriptors` from `fod.gbw`/`fod.scfp_fod` and the geometries.

**Open decision:** the implementation of the Hirshfeld reference densities.

## 8. Workflow

The start structures are already in `structures/` (section 8a), so ORCA work starts directly:

```bash
python datagen/make_inputs.py structures/pilot.csv runs --stage opt   # then run the jobs
python datagen/make_inputs.py structures/pilot.csv runs --stage sp    # needs converged opt
python datagen/make_inputs.py structures/pilot.csv runs --stage cc    # needs opt + tpssh
python datagen/collect.py structures/pilot.csv runs --out data [--nevpt2 nevpt2.csv]
```

Use `structures/complexes.csv` instead of `structures/pilot.csv` for the full set.

- Every stage writes `runs/jobs_<stage>.txt` listing the new job directories.
- Run each job inside its directory as `orca <job>.inp > <job>.out`, using the full ORCA path for parallel runs.
- Re-running a stage skips finished jobs. Jobs that are not ready yet are listed on stderr.
- `--jobs` restricts a stage, e.g. `--stage sp --jobs tpssh fod` to postpone r2SCAN0.

**Warnings written by `collect.py`:**

| warning | meaning |
|---|---|
| `<job>:missing` / `<job>:failed` | no output / not terminated normally |
| `<job>:no_<value>` | expected value not found in a finished output |
| `tpssh:s2` | \|⟨S²⟩ − S(S+1)\| > 0.1 |
| `fod:nfod_mismatch` | \|N_FOD(ORCA) − N_FOD(recomputed)\| > 0.01 |
| `fod:spin` | smearing changed $N_\alpha - N_\beta$, i.e. the multiplicity is not conserved |
| `t1:missing` | informative only; does not affect `ok` |

## 8a. Candidate structures

`structures/complexes.csv` with `structures/xyz/` is the candidate list. It is produced in two steps:

```bash
python datagen/build_structures.py structures --molsimplify <molsimplify>   # molSimplify 2.0, legacy CLI
python datagen/preopt_xtb.py structures/complexes.csv --xtb <xtb>           # GFN2-xTB 6.7.1
```

**`build_structures.py`** builds a grid of metal/oxidation state × ligand set (`METALS`, `SETS`, `GRID`), plus explicit non-innocent cases (`EXTRA`). The rules are:
- Charge window: −1 ≤ q ≤ 3.
- At most 52 atoms.
- The metal must have exactly the target coordination number (bond criterion as below), otherwise the structure is rejected.
- molSimplify's own charge must equal the manifest charge.
- Isomers:
  - A repeated ligand name fixes the placement order: cis-M(bpy)2X2, trans-MX2(dmpe)2, trans square-planar MX2L2.
  - Other octahedral MX2L4 come out trans by default.
- CN4 complexes with geometry `sqp;thd`:
  - d4: square planar for all states.
  - Otherwise: square planar for the lowest multiplicity (`<id>_m<mult>.xyz`) and tetrahedral for all other states.
- Ligands with radical or NO character (`EXTRA`):
  - NO is counted as neutral. The formal oxidation state and dⁿ are those of the metal without NO.
  - Ni dithiolene and Co catecholate/semiquinone carry explicit charges.
  - All of them are flagged `noninnocent`. Dithiolene and catecholate are also flagged `mr`.
- The complexes are hypothetical model complexes chosen to span the ligand field. Real complexes (e.g. literature spin-crossover systems) should be added by hand with the same manifest columns.

**`preopt_xtb.py`** optimises each file with GFN2-xTB for the state it serves:
- `<id>_m<mult>.xyz` for that multiplicity.
- `<id>.xyz` for the highest remaining multiplicity.

A result is accepted only if the bond graph is unchanged. The criterion is d < 1.3 (r_A + r_B), or 1.4 (r_A + r_B) for metal bonds, with Pyykkö covalent radii. Otherwise the molSimplify structure is kept. `structures/preopt_log.csv` lists the status of every file. GFN2-xTB only provides start geometries; the r2SCAN-3c optimisation decides the final structure of each state.

## 9. Pilot (5 complexes, `structures/pilot.csv`)

The pilot complexes are:
- `fe2_h2o6`: [Fe(H2O)6]²⁺
- `fe2_mecn6`: [Fe(NCMe)6]²⁺, the nitrile model; HCN is not in the molSimplify library
- `mn2_h2o6`: [Mn(H2O)6]²⁺ (negative control)
- `ni2_cl2ph3_2`: NiCl2(PH3)2 (square planar vs. tetrahedral)
- `fe2_por`: Fe(II) porphine (intermediate spin, MR test)

Check before the full run:
1. ORCA accepts TPSS0, r2SCAN0, D4 for TPSS/TPSS0, and DLPNO-CCSD(T1). Fallback for missing D4 parameters: D3BJ for all three $a_x$-scan functionals.
2. The FOD job really uses 5000 K and TPSS/def2-TZVP.
3. `collect.py` reports no `no_*` warnings, i.e. the parser patterns match the real output.
4. `nfod` and `nfod_check` agree, and there is no `fod:spin`. If `fod:spin` appears, smearing in ORCA does not keep separate Fermi levels per spin. Then the UKS-FOD protocol has to be revised before any production run.
5. Our CPS/CBS numbers reproduce ORCA's built-in extrapolations for one state.
6. Timings decide whether QZ is feasible for all complexes or only for a calibration subset.

## 10. Checks after the full run

- Condition number of the design matrix per $d^n$ class.
- At least 20 negative controls ($N_\mathrm{FOD}$ < 0.1 in all states).
- Correlation of $\Delta N_\mathrm{FOD}$ with $d\Delta E/da_x$.

## References

- S. Grimme, A. Hansen, *Angew. Chem. Int. Ed.* **2015**, 54, 12308 (FOD).
- A. Altun, F. Neese, G. Bistoni, *J. Chem. Theory Comput.* **2020**, 16, 6142 (PNO extrapolation).
- M. Drosou, C. A. Mitsopoulou, D. A. Pantazis, *J. Chem. Theory Comput.* **2022** (DLPNO for spin-state energetics), https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9202354/
- M. Radoń et al., *Chem. Sci.* **2024** (SSE17), https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11577268/
- W. Jiang, N. J. DeYonker, A. K. Wilson, *J. Chem. Theory Comput.* **2012**, 8, 460 (MR diagnostics for 3d metals).
- F. Neese, E. F. Valeev, *J. Chem. Theory Comput.* **2011**, 7, 33 (def2 extrapolation; verify β in the pilot).
- S. Grimme, A. Hansen, S. Ehlert, J.-M. Mewes, *J. Chem. Phys.* **2021**, 154, 064103 (r2SCAN-3c).
- ORCA manual, FOD: https://www.faccts.de/docs/orca/6.0/tutorials/prop/fod.html
