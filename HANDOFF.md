# HANDOFF — resume LyaLenser on BNL RACF (or anywhere)

Written 2026-09-11 on the workstation where the project started. Read this first, then `PROGRESS.md` (chronological
log and decisions), `MEMORY.md` (machine caveats), `PLAN.md` (measurement plan), `IMPLEMENTATION.md` (pipeline spec).
The repo is `github.com/slosar/LyaLenser` (private); everything below is committed there.

## What this project is
Detect weak lensing of the DESI DR1 Lyman-alpha forest by cross-correlating a quadratic (pair-based) lensing
estimator built from forest pixel pairs with CMB lensing maps (ACT DR6, Planck PR4), after deprojecting the forest's
response to long-wavelength density modes with a kernel-matched quasar template. Theory and forecasts are in
`report/main.pdf` (23 pages, Secs. 1-5 + appendices); the headline forecast with DR1's measured coverage and noise is
S/N ~ 1 for DR1 and ~4 for a complete DESI, so DR1 is a pipeline/upper-limit exercise.

## Where things stand (state of the tree)
- **Theory/report**: done through three adversarial reviews (`report/reviews/codex_review_{1,2,3}.md`), all
  incorporated. Mathematica checks in `mathematica/*.wls` all pass.
- **Data**: all downloaded and verified (see `MEMORY.md`, `code/fetch/fetch_data.sh`): DR1 Lya deltas (1028 files,
  428,403 forests), DR1 quasar catalogues + LSS randoms, ACT DR6 lensing release incl. 400 baseline sims, Planck PR4
  2018-like maps. On the workstation they live in `/data/LyaLenser/raw`; the code finds them through the
  `LYALENSER_DATA` environment variable (`code/paths.py`).
- **Pipeline (Stage A = mock validation of the estimator)**: `code/pipeline/` after three implementation rounds
  (two by gpt-5.6-sol, one by gpt-6-astra). 36 unit tests pass. The frozen acceptance run of round 3 gives
  **20/30 gates PASS; Stage A is NOT accepted and Stage B (real data) is blocked.** See `report/mock_validation.md`,
  `GATES.md`, `code/pipeline/NOTES.md` (all three iterations), `report/reviews/codex_review_4_code.md` (astra's code
  review of round 2), `report/reviews/impl_round3_summary.md`.
- **Diagnosis of the remaining failures** (by the Claude session, on the saved round-3 products in
  `/data/LyaLenser/mocks/iteration3/`): see `code/pipeline/ITERATION4.md` "Findings" — (F1) matched quasar template
  recovers only 0.87 of the same-range slab convergence; (F2) validation used the shared sightline/template selection;
  (F3) the normalisation was tested with A_true up to 10 where the remapping is non-linear, and the A=1 recovery is
  only known to ~40%; (F4) measured vs analytic xi' differ by 19% in the response integral; (F5) covariance ratio
  0.696 vs gate 0.7.

## What to do next (in order)
1. **Round 4 of Stage A** per `code/pipeline/ITERATION4.md` (7 items incl. the map-level template gate, disjoint
   selection baseline, A_true in {0, 0.5, 1, 2} with dense noiseless mocks, the xi' discrepancy, and the cluster entry
   points `--phase seed/collect`). The user's standing instruction: implementation rounds now use **gpt-6-astra**
   (via the Codex plugin, `codex-companion.mjs task --write --model gpt-6-astra --effort high`; on RACF without Codex,
   do it directly). Keep the GATES.md discipline: freeze choices on development seeds 100-104, declare gates before
   the fresh ensemble, absolute statistics only.
2. Run the 60-seed ensemble as condor jobs (`condor/README_RACF.md`, `condor/seed.sub`, `condor/run_seed.sh`;
   32 GB, 8 CPUs per seed) once the per-seed entry point exists.
3. When all gates pass: an astra review of code + validation, then **Stage B** (`IMPLEMENTATION.md` Sec. 3):
   DR1 deltas -> measured xi -> pair catalogue -> templates (ACT, Planck, matched from DR1 quasars) -> amplitudes,
   bands, curl, injections, 400 ACT random templates, jackknife; single slab first, then tomographic sub-slabs with
   disjoint selection. Expect ~1 sigma; the response term (~2 sigma) is the positive control.

## How to resume on RACF (done 2026-09-11; the paths below are the live ones)
```bash
cd /gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser
export LYALENSER_DATA=/gpfs/mnt/gpfs02/astro/workarea/anze/Data/LyaLenser   # raw/ and mocks/ rsynced from the workstation
export PATH=/gpfs/mnt/gpfs02/astro/workarea/anze/envs/lyalenser/bin:$PATH   # env on gpfs (home is at quota), numpy<2 pinned
cd code/pipeline && NUMBA_NUM_THREADS=4 python -m pytest -q tests       # 36 passed (2.7 min on the 4-core login node)
python run_mock_validation.py --phase rebuild --mock-root $LYALENSER_DATA/mocks/iteration3   # reproduces report/mock_validation.md from disk, no simulation
```
Machine caveats (quota, condor, Codex sandbox) are in `MEMORY.md`; the round-4 execution split (Codex implements at
smoke scale, this session runs the scale-1 campaign on HTCondor) is in `code/pipeline/ITERATION4.md`.
Not available on RACF: the Mathematica MCP (derivation checks; all scripts already verified), the Codex plugin
(reviews/implementation rounds). If you continue without Codex, do the round-4 work yourself and record in
`PROGRESS.md` that the adversarial-review step was skipped or done differently.

## Conventions and traps (short list; details in MEMORY.md and IMPLEMENTATION.md Sec. 0)
- Pair separation is theta_a - theta_b (points from b to a); the response sign is unit-tested (Mathematica check 11).
- Lensing: observed field = true field at theta_obs + alpha; mocks sample at +alpha, injections shift positions by
  -alpha. kappa = -laplacian(phi)/2, phi_lm = 2 kappa_lm / (l(l+1)).
- The public DR1 deltas are named DELTA_BLIND with BLINDING=desi_y1, but picca's desi_y1 blinding only alters the
  fiducial distance conversion inside picca_cf; the arrays are the measured fluctuations (verified in picca source).
- DESI data server returns intermittent 503s under parallel downloads; `wget --retry-on-http-error=503`.
- Forecast code (`code/three_tracer.py`) uses the released ACT N_L and the measured DR1 coverage/noise
  (`report/data_delta_summary.json`, produced by `code/data_checks/delta_summary.py`).
- Never treat a "PASS with SEM larger than the quantity" as a pass; the acceptance rows must show achieved bounds.

## Key files
`report/main.tex|pdf`, `report/reviews/`, `report/mock_validation.{md,json}`, `GATES.md`, `IMPLEMENTATION.md`,
`PLAN.md`, `code/pipeline/{ITERATION2,ITERATION3,ITERATION4}.md`, `code/pipeline/NOTES.md`, `code/paths.py`,
`code/fetch/fetch_data.sh`, `condor/`, `MEMORY.md`, `PROGRESS.md`.
