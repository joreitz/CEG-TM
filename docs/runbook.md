# Runbook: unattended operation (agreed 2026-11-23, user away for one week)

## Reporting
- Every evening: short summary in the chat (what finished, what failed, what was submitted, open questions).
- Push notification to the phone only for things that are completely unexpected or need a decision.

## Autonomy
Allowed without asking:
- Restart / resubmit optimisations and single points according to the convergence protocol below.
- Start the sp stage (ORCA 6.1.1) once the opts of a state are clean, then the cc stage.
- Write evaluations and reports.
- Add `exclude` / `keep` entries to `structures/state_overrides.csv` under the agreed rules
  (ligand dissociation, 4+2 elongation; spin contamination is handled by collect.py thresholds 0.2 / 0.5).

Ask first (push notification, then wait):
- Changing the method (functional, basis, thresholds, protocol in data_plan.md).
- Deleting data other than outputs that are superseded by a rerun.
- Merging pull requests.
- Anything expensive (NEVPT2 subset, full cc stage of the whole set before the pilot checks pass).

## Convergence / multireference protocol (user's instructions)
1. Learn first how experts handle convergence problems for transition-metal (TM) and TM–TM complexes
   (ORCA manual, literature, best practice). Focus (user, 2026-11-23): the correlated methods, i.e. DLPNO-CCSD(T)
   (reference determinant: UHF/QRO/UKS orbitals, stability, CC iterations), CASSCF (active space, orbital
   convergence) and NEVPT2; DFT SCF problems only as far as they feed into these.
2. Pick one difficult species and run test calculations varying convergence criteria and settings to find the
   most reliable settings.
3. Apply those settings to a test group: one complex from each category. Check that everything converges.
4. Only then apply them to all other complexes.

## Cluster
- marvin, workdir `/lustre/scratch/data/s6jjreit_hpc-agostic/test/CEG-TM`, submit only with `suborca`
  (`export HOSTNAME=$(hostname)`, `module load Python/3.13.5-GCCcore-14.3.0` first).
- Opts: `-c 12 -m 5 -t 8h -v 6.1.0`; energies (sp, cc): ORCA 6.1.1 (suborca default); tpssh `-o`, fod `-f`.
- At most 300 queued jobs per user; keep below ~280.

## Expert practice (research 2026-11-23; sources: ORCA 6.1 manual CASSCF/NEVPT2/MDCI, ORCA input library,
## Radoń et al. Chem. Sci. 2024 (SSE17), Drabik/Radoń JCTC 2022, Feldt EJIC 2022)
DLPNO-CCSD(T1), open shell:
- ORCA converts the UHF (or UKS) determinant to QROs automatically; (T1) is required for spin states.
- The reference must be a converged, stable determinant: TightSCF, stability analysis, check <S2>.
  UHF of TM complexes often lands in an excited / strongly contaminated solution -> try several guesses
  (PModel, MORead of the TPSSh orbitals, TRAH) and keep the lowest stable one.
- HF vs KS orbitals in the reference: no consistent improvement on SSE17 (CCSD(T) MAE 1.5 kcal/mol with HF);
  KS (BP86) orbitals helped for Co(III). Test both, keep HF unless KS is clearly more robust.
- CC iterations: default MaxIter 50, MaxDIIS 7, level shift 0.2 (%mdci MaxIter / MaxDIIS / LShift).
  Non-convergence usually points to an exotic electronic structure: inspect, do not just raise MaxIter.
CASSCF:
- Never start from the default guess. Start orbitals: AVAS(Valence-D) / AVAS(Double-D), or natural orbitals
  (UHF/QRO, MP2) via MORead; build large spaces from converged smaller ones.
- Active-space check: natural occupations between 0.02 and 1.98; orbitals at ~0 or ~2 cause convergence
  problems and should be removed. Metal 3d (+ sigma-bonding e_g partners, + double d-shell for late metals).
- Convergers: default SuperCI_PT (no level shift); if it oscillates, `! TRAH` (no AutoTRAH for CASSCF);
  alternative SuperCI + DIIS with a large level shift (ShiftUp/ShiftDn 2.0).
- RI: TrafoStep RI needs a /C or /JK aux basis (not /J only).
- Hand-picked active spaces (user, 2026-11-23): depending on convergence and on which orbitals take part, pick
  the active orbitals by hand. Identify them from the Loewdin MO compositions (metal d / donor p / ligand pi,pi*
  / 4d shell; print with %output Print[P_OrbPopMO_L] 1 end) and from UHF-natural / CASSCF natural occupations;
  move them into the active space with %scf rotate {i, j, 90} end on MORead orbitals; re-test convergence.
  Document every choice (orbital index, composition, occupation) in the evening summary; push-notify the user
  in doubtful cases (non-innocent ligands, porphyrins, NO, strongly covalent bonds).
NEVPT2:
- SC-NEVPT2 is intruder-free; FIC/DLPNO-NEVPT2 need RI (/C or /JK); RIJK more accurate than RIJCOSX.
- False intruders (negative denominators): tighten D3TPre/D4TPre, avoid cumulant approximations,
  last resort `imaginary 0.2` (FIC). Run NEVPT2 on a converged CASSCF (MORead NoIter).

## Convergence test plan (step 2): fe2_bpy2ncs2 (singlet / triplet / quintet, r2SCAN-3c geometries)
A  TPSSh/def2-TZVP single points (production sp job, gives the guess orbitals).
B  Reference SCF def2-TZVPP (RHF for the singlet, UHF otherwise), TightSCF + StabPerform:
   B1 PModel guess, default convergers;  B2 MORead TPSSh (production);  B3 MORead TPSSh + TRAH.
   Record: converged?, iterations, E, <S2>, stability, which solution is lowest.
C  DLPNO-CCSD(T1)/def2-TZVPP TightPNO TCutPNO 1e-6: C1 on the lowest stable UHF/RHF; C2 on TPSSh orbitals (NoIter).
   Record: CC iterations, T1, gaps C1 vs C2.
D  CASSCF/def2-TZVP (TrafoStep RI, def2/JK): D1 CAS(6,5) from AVAS(Valence-D); D2 CAS(6,5) from UHF natural
   orbitals; D3 AVAS(Double-D). SuperCI_PT, fallback TRAH. Record convergence and occupations.
E  SC-NEVPT2 and DLPNO-NEVPT2 on the converged CASSCF of each state; compare gaps with C.
Then step 3: one complex per category with the best settings; step 4: all.

## Log
- 2026-11-23 fe2_terpycl2 triplet (r2SCAN-3c, start geometry): SCF has several close UKS solutions.
  SOSCF + MaxIter 500 converges from PModel/PAtom/Hueckel to -2926.5110..-2926.5114 Eh with <S2> ~3.0 (ideal 2),
  spin Fe 2.2 / terpy -0.18 (no ligand radical). TRAH and BS(3,1) (Fe(III)/terpy radical start) both stall at a
  flat point at -2926.5095 Eh (gradient 5e-5..1e-4 for >300 macro iterations, 1 h limit). Quintet is ~28 kcal/mol
  lower. Opt resubmitted by hand with %scf MaxIter 500 (only change); the state will very likely fail the
  <S2> > 0.5 rule at the TPSSh stage.
- 2026-11-23 fe2_bpy2ncs2 round 1 (TPSSh/def2-TZVP; HF/def2-TZVPP; CASSCF/def2-TZVP):
  TPSSh stable for all states, <S2> 2.06 / 6.13; gaps vs singlet: triplet +16.7, quintet +14.3 kcal/mol.
  B1 PModel HF: singlet RHF UNSTABLE -> ORCA restarted as UHF, <S2> 2.97 (broken-symmetry singlet; unusable as
  CC reference). Triplet UHF stable but <S2> 4.70 (ideal 2), quintet 8.34 (ideal 6): UHF heavily contaminated.
  -> Production cc header (StabPerform for RHF singlets) would turn the singlet reference into a BS-UHF: fix
  needed before the cc stage (no UHF restart for singlets). Tested in round 2 (C1 = exact production job).
  D1 CAS(6,5) AVAS(Valence-D): converges in 12-15 macro iterations, no problems. Occupations singlet
  1.98/1.98/1.97/0.035/0.035, triplet 1.97/1.94/1.04/0.98/0.06. SC-NEVPT2 on CAS(6,5): quintet lowest,
  singlet +4.0, triplet +19.1 kcal/mol -> known high-spin bias of a 3d-only active space (exp.: LS ground state).
- Round 2 submitted: C1 production cc_tz6, C2a TPSSh/def2-TZVPP (KS reference), D3 AVAS(Double-D) CAS(6,10),
  comp-d1 Loewdin MO compositions of the CAS(6,5) orbitals (to hand-pick the sigma e_g partners for CAS(10,12)).
- comp-d1 (Loewdin MO compositions of the converged CAS(6,5) orbitals, triplet/quintet; singlet RHF job printed no
  per-MO block): active MOs 121-125 are 90-98 % Fe 3d. No clean sigma(e_g) partners among the inactive orbitals:
  the Fe-N sigma character is spread over ~6 doubly occupied MOs with only 3-6 % Fe d (25-48 % donor N); C2
  symmetry, canonical CASSCF orbitals. Plan: take the AVAS(Double-D) result (D3), rotate the two inactive MOs with
  the largest Fe-d x donor-N weight into the active space -> CAS(10,12) (3d + 3d' + 2 sigma), check natural
  occupations (sigma ~1.98, 4d ~0.02).
- 2026-11-23 evening: all opts finished (check_opts round 2, 655 states). Restarts converged; fe2_terpycl2 m3 now
  converged (<S2> 2.24, the MaxIter 500 SCF fix worked). New overrides: exclude mn2_nh3_4sme2 m4, fe2_en2cl2 m3,
  fe2_nh3_4sme2 m3 (dissociation), fe2_nh3_4ox m3, ru2_h2o6 m3 (4+2), co3_nh3_4ox m5 (oxalate C-C breaks again
  from the idealised start, 2.04 A); keep fe2_h2o4f2 m3 (Fe-O 2.48 A). Octahedral Ni(II) singlets: all 10 go
  square planar also with TPSSh-D4/def2-SVP -> singlets excluded, the 10 octahedral Ni(II) complexes dropped.
  Metallocenes: all 12 states converged, FeCp2 singlet Fe-C 2.05 A (exp. 2.06), <S2> clean.
  Resubmitted: fe3_nh3_4ox m4 (SCF stopped at 125 it. with dE 2e-8 -> %scf MaxIter 500), fe3_terpycl2 m4
  (4 h time limit after 8 cycles, ~50 s per SCF iteration -> restart from last geometry, -t 24h).
  Manifest now 248 complexes / 629 states. sp stage started: tpssh (-o) + fod (-f), 1251 inputs, throttled loop
  ~/ceg_sp_submit.sh on login01 (log ~/ceg_sp_submit.log). tpss/tpss0/r2scan0 postponed until the keywords are
  checked on one complex.
- 2026-11-23 ~21:40 message from the cloud session (with the user's consent) on the cc header: (1) singlets RHF
  without UHF restart, record the instability + T1 as MR flag (-> NEVPT2 subset); open shells UHF + stability;
  (2) check the TPSSh singlets for an RKS->UKS switch (tpssh header also has StabPerform; collect.py checks <S2>
  only for mult > 1); consistent choice = restricted singlets, rerun affected states without stability analysis;
  (3) tighter cc: VeryTightSCF, RIJK def2/JK instead of RIJCOSX, %mdci STol 1e-7, check frozen core, show
  runtime on the pilot before the full cc stage. Needs the user's confirmation before the cc stage.
  Findings: C1 singlet: RHF unstable (lowest stability eigenvalue -0.039 Eh), ORCA restarted as UHF and needed
  237 SCF cycles. ORCA 6.1.1 knows %scf StabRestartUHFifUnstable (default evidently true here). Frozen core:
  chemical core, Fe [Ne] frozen only -> 3s3p are correlated. CC tolerances printed: energy 1e-6, max residuum
  1e-5 (STol). Triplet UHF reference: DIIS did not converge in 50 it., TRAH took over.
  C3 submitted (fe2_bpy2ncs2, all 3 states, same resources as C1): RHF/UHF DLPNO-CCSD(T1) TightPNO def2-TZVPP
  def2-TZVPP/C def2/JK RIJK VeryTightSCF MORead TPSSh; singlet StabPerform + StabRestartUHFifUnstable false;
  %mdci TCutPNO 1e-6 STol 1e-7. Compare runtime and gaps with C1.
  TPSSh singlets: so far only fe2_bpy2ncs2 finished (stable). Check all finished tpssh _1.out for
  "wavefunction is unstable" as they come in.
- 2026-11-23 20:30 CET check: sp 46 ok / 1 failed. cr2_nh3_4 fod m3 (square-planar Cr(II) triplet): SCF with
  smearing oscillates (TRAH and SOSCF are disabled with fractional occupations), not converged after 125 it.
  -> resubmitted with SlowConv + %scf MaxIter 500 DIISMaxEq 15 (convergers only). fe3_nh3_4ox m4 opt converged
  with %scf MaxIter 500; its tpssh/fod inputs still to be written (after the sp loop has finished).
- 2026-11-23 ~20:40 CET second message from the cloud session: user wants Dunning instead of def2 for the DLPNO
  reference, non-relativistic (3d: cc-pwCVnZ; Ru/Rh: cc-pwCVnZ-PP + matching ECP (ECP28MDF, ORCA "SK-MCDHF-RSC"),
  check in output; ligands cc-pVnZ, no aug; AutoAux; otherwise as C3); collect.py BETA 3.0, E_SCF from QZ, F 1.5;
  data_plan.md to update. Needs the user's confirmation before the cc stage.
  C4 submitted (fe2_bpy2ncs2, 3 states, TZ -t 24h and QZ -t 72h, 12 cores, 8 GB/core):
  ! RHF|UHF DLPNO-CCSD(T1) TightPNO cc-pVnZ AutoAux RIJK VeryTightSCF MORead; %basis NewGTO Fe "cc-pwCVnZ";
  MORead TPSSh/def2-TZVP orbitals (basis projection); singlet StabRestartUHFifUnstable false; STol 1e-7.
  Open point for the proposal: S with plain cc-pVnZ (cc-pV(n+d)Z is the usual choice for 2nd-row atoms).
- 2026-11-23 23:35 CET check: sp 685 ok / 7 failed, loop at 919/1251. kwcheck fe2_h2o6 TPSS/TPSS0/r2SCAN0: all 9 ok
  (keywords fine). TPSSh singlets: 17 of 36 finished singlets were unstable -> ORCA re-converged to broken-symmetry
  UKS (<S2> 0.21-0.99, E lower by ~9-54 kcal/mol vs the first RKS energy), mostly Fe(II); collect.py would not see
  it (no <S2> check for mult 1). Push sent to the user (decision: restricted singlets without UHF restart, as
  proposed by the cloud session; opts were RKS for singlets, so restricted is consistent). Open-shell
  instabilities (e.g. cr2_acac2 m3 -> <S2> 2.97, fe3_cl4) follow the agreed rule (lowest stable UKS + S2 flags).
  C1 singlet (production cc header) crashed in MDCI after the RHF->UHF switch ("UHF-GUESS: PNO flag is
  detected") -> production header unusable for unstable singlets.
  Basis checks: Ru cc-pwCVTZ-PP + NewECP "SK-MCDHF-RSC" works, only O 1s frozen (Ru 4s4p correlated).
  cc-pV(T+d)Z / (Q+d)Z are not in the ORCA 6.1.1 library (all spellings rejected); built-in cc-pVTZ/QZ for Cl are
  the plain sets [5s4p2d1f] / [6s5p3d2f1g] -> (n+d) would need a basis file from the Basis Set Exchange (GTOName).
  Resubmitted: 5 fod nearly converged (dE ~1e-7 oscillation) with MaxIter 500; cr2_nh3_4 fod m3 (still oscillating
  after 500 it. with SlowConv) with MORead TPSSh gbw + level shift 0.3; fe2_nh3_4cl2 tpssh m3 (unstable in all
  5 stability tries) with TRAH.
- 2026-11-24 user decisions: (1) singlets restricted; RKS/RHF instability stored as MC/MR flag ("als solche
  abspeichern und fortfahren") -> make_inputs.STAB, collect.py mr_rks / mr_rhf (branch restricted-singlets,
  deployed to marvin, not yet committed). (2) P/S/Cl: plain cc-pVnZ (no (n+d) basis files), agreed.
- 2026-11-24 ~00:30: restricted-singlet test: StabRestartUHFifUnstable false alone aborts in the population
  analysis (missing .scfr, ORCA 6.1.1 bug); with ! NoPop it terminates normally, stays RKS and prints "wavefunction
  is unstable" (fe2_h2o6 RKS -1722.01130 vs BS-UKS -1722.01799 Eh = 4.2 kcal/mol; the earlier "9-54 kcal/mol"
  numbers were wrong: they compared SCF and final (with D4) energies). STAB[True] now "! NoPop" + flag. 27
  unstable singlet tpssh outputs removed (superseded) and 28 restricted reruns submitted. tpss/tpss0/r2scan0
  stage started (1875 inputs, loop ~/ceg_sp2_submit.sh, log ~/ceg_sp2_submit.log).
  Preliminary N_FOD statistics for the user (scratchpad nfod_prelim.png / nfod_stats.py): CN alone no
  correlation (rho -0.02); spin class dominates (IS 1.69 > LS 1.21 > HS 0.98); N_FOD vs <S2> deviation rho +0.52;
  RKS-unstable singlets median N_FOD 1.45 vs 0.61 (p < 1e-6).
