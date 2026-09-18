# RESUME HERE (written 2026-09-18, after iteration 11)

Everything is committed; the working tree is clean unless the attribution run below has finished and its JSON
still needs copying. Perlmutter is down until 2026-09-23; the whole Stage B chain runs on the workstation in the
NaMaster-capable env `/data/LyaLenser/envs/lyalenser/bin/python` (MEMORY.md).

**State of the measurement.** DR1 low-z forest lensing, single slab **1.96 < z < 3.0** (every DR1 pixel below
z = 3; z_eff = 2.348 = weighted mean pixel redshift, the templates' source plane), iteration 11:
**A = 0.346 +- 0.474** (Fisher 0.379), `report/lowz/lowz.pdf` (19 pages), `report/stageb/dr1_lowz_v7.{json,md}`;
iteration 10 (2.1 < z < 3.0, six slices) gave 0.271 +- 0.519. NOTES "Iteration 11" has the details.
Iteration 10 (user's four requests, 2026-09-18; `code/pipeline/NOTES.md` "Iteration 10" is the technical record):
NaMaster spectra (`code/stageb/nmt_spectra.py`), the redshift-evolving correlation table (`code/pipeline/xi_zevol.py`,
layered `XiTable`; b_F ∝ (1+z)^3.5; A x B / A x A ratio 0.84 -> 0.92 at fixed z), BGS + BOSS tracers in six slices
with per-coverage-class Wiener weights and the shared-object noise term (`code/stageb/lowz_catalogues.py`,
templates in `/data/LyaLenser/lowz_v2`), and the deflection-template validation against ACT DR6 and Planck PR4
(`code/stageb/deflection_cmb_check.py`, `template_prediction.py`, `cmb_map_checks.py`; `report/stageb/deflection_cmb_check.json`,
`cmb_map_checks*.json`): **A_L = 0.99 +- 0.04 (ACT), 0.92 +- 0.03 (Planck)** relative to the exact prediction for
the template (class fractions inside each overlap; ACT map masked on input), B/E <= 6 %, maps consistent on common
sky and by region. The first-pass 0.83 / 0.71 were two errors of mine (NOTES "Correction to the deflection
validation"). Nothing applied to the result.
No mock re-validation was run (user decision).

To reproduce the current result (~78 min, 34 GB on the workstation; `slurm/dr1_lowz_v7.sbatch` on Perlmutter,
repo clone `/global/cfs/cdirs/m4895/users/anze/LyaLenser_iter11` made with `gh repo clone`, data in
`LyaLenser_data/{lowz_v3,stageb/basis_dr1_ab_z196.h5}` already copied there):
```bash
cd code/stageb && NUMBA_NUM_THREADS=22 /data/LyaLenser/envs/lyalenser/bin/python run_dr1_lowz.py \
    --out $LYALENSER_DATA/stageb/dr1_lowz_v7 --lowz $LYALENSER_DATA/lowz_v3 --basis $LYALENSER_DATA/stageb/basis_dr1_ab_z196.h5 \
    --regions lya lyb --xi-correction spline --randoms 40 --nside-jk 8 --zmin 1.96 --zmax 3.0 --zeff 2.3476
```
Templates: `python lowz_catalogues.py --out $LYALENSER_DATA/lowz_v3 --zref 2.3476 --tracer-zmax 1.6` (11 tracers,
~6 min); basis: `python build_basis_dr1.py --out .../basis_dr1_ab_z196.h5 --regions lya lyb --zmin 1.96` (15 min);
checks: `cmb_bias_check.py`, `cmb_map_checks.py`, `deflection_cmb_check.py` with `--lowz .../lowz_v3` (~15 min).
Evolution diagnostic on a saved run: `python xi_zevol_dr1.py --run .../dr1_lowz_v4` (`report/stageb/xi_zevol_dr1.json`).

**In flight / to pick up:**
1. Attribution run done (`report/stageb/dr1_lowz_v6_flat.json`): new templates + flat table A = 0.274 +- 0.537,
   so the templates carry the whole v4 -> v6 change; the evolving table trims the error by 3 %.
2. The user's sanity checks (to be specified). Candidates already flagged: the 8 % Planck deficit and the mild
   rise of A_L with L (clustering amplitude below the fiducial model, or low-ell systematic power in the tracer
   auto-spectra; ELG 0.8-1.1 cross/auto 0.65), the 8 % residual A x B / A x A difference at fixed z, the degeneracy of the base
   amplitude with the spline correction in the evolving fit (quote the base-only fit).
3. Mock validation of region B, the five bands, the evolving table, the enlarged tracer set and the 1.96-2.1
   extension (the mock box starts at z = 2.1); the iteration-7 mocks live on Perlmutter (`mocks/iteration7/`).
   The paper (`Paper/`, its own git repo, Overleaf remote) still carries the iteration-10 numbers: update it to v7.
4. Older items: debias the response matrix for the kernel-fit attenuation, the joint slice fit, the r_perp tilt
   of the corrected projection, tomographic sub-slabs, GATES v8.

---

# HANDOFF — resume LyaLenser (written 2026-09-12 on BNL RACF; valid for RACF, NERSC Perlmutter or elsewhere)

Read this first. Then: `PROGRESS.md` (chronological log and decisions), `code/pipeline/NOTES.md` section
"Iteration 4" (the technical diagnosis of the round-3 failures and what was changed), `GATES.md` (the frozen
acceptance protocol the running campaign is judged by), `condor/CAMPAIGN.md` (the campaign's job list),
`MEMORY.md` (machine caveats). The repo is `github.com/slosar/LyaLenser` (private); everything is committed.

## What this project is
Detect weak lensing of the DESI DR1 Lyman-alpha forest by cross-correlating a quadratic (pair-based) lensing
estimator built from forest pixel pairs with CMB lensing maps (ACT DR6, Planck PR4), after deprojecting the forest's
response to long-wavelength density modes with a kernel-matched quasar template. Theory, forecasts and the estimator
are in `report/main.pdf` (done, three adversarial reviews). Forecast: S/N ~ 1 for DR1, ~4 for complete DESI, so DR1
is a pipeline / upper-limit exercise. **Stage A** = validation of the estimator on self-lensed mocks against the
gates in `GATES.md`; **Stage B** = the real DR1 x ACT/Planck measurement, allowed only after Stage A passes.

## State on 2026-09-15
- Stage A rounds 1-3 (Codex gpt-5.6-sol x2, gpt-6-astra x1) ended at 20/30 gates. Round 4 was started by astra and
  finished by the Claude session after the OpenAI spend cap killed the Codex task (user decision). Round-4 code:
  `code/pipeline/campaign4.py` (idempotent, provenance-hashed phases `dev-seed / freeze / seed / control / collect`),
  `random_streams.py`, `template_audit.py`, `tests/test_iteration4.py` (51 tests pass), `GATES.md` v4.
- The two round-3 failure mechanisms were found at smoke scale (0.25) and fixed before the campaign (NOTES.md):
  (1) the 0.87 template deficit = the mock's lognormal quasar model on raw 2 Mpc/h cells (b sigma = 2.8) coupling to
  the 40-bin radial normalisation -> lognormal input now smoothed radially by 8 Mpc/h, audit gated in 40<=L<=300;
  (2) the 8 % xi' excess = sightline quasars drawn inside their own forest range (self-proximity term) -> sightline
  quasars now drawn behind the slab. Verified: template band coefficient 0.99-1.01 on scale-1 dev seeds, xi'
  coefficient 1.08 -> 0.98 at smoke.
- **The scale-1 acceptance campaign is COMPLETE** (2026-09-12: RACF for dev seeds, freeze and sparse 0-6; NERSC
  Perlmutter preempt QOS for the rest, ~1.5 h wall). **39/48 gates pass; Stage B remains blocked by protocol.**
  `report/mock_validation.{md,json}` (+ `report/figures/mock_iteration4_*.pdf`) is the acceptance table; the
  reading of it is in NOTES.md "Scale-1 campaign". In short: template prerequisites and the xi' model are solved,
  the deprojected estimator is unbiased (null -0.63 +- 0.67 A, recovery 0.70 +- 0.64 A, covariance ratio 0.95);
  what remains is a 4.2 +- 1.6 % multiplicative normalisation bias (dense slope 1.042, linear in A) traced to the
  smoothed-derivative design of the measured xi (the frozen smoothing width was chosen on 5 sparse seeds with +-4 %
  precision), plus gates whose frozen tolerances are unattainable at this ensemble size (0.3 A bounds need ~750
  seeds) or ill-posed (one-SEM sampled coefficient, 0.03 endpoint means with SEM 0.16) and one with an inadequate
  reference (analytic baseline). Products: `$LYALENSER_DATA/mocks/iteration4/` on RACF and on NERSC (identical).
- **Iteration 5 done (2026-09-14, Perlmutter)**: model-shaped xi table (`xi_fit.py`), kernel cut r_perp >= 3, GATES
  v5 with N-matched tolerances. **34/38 required gates pass**; normalisation 0.980 +- 0.014 (table-independent);
  deprojection unbiased; the four failing rows are two mis-set 1 A precision bounds (protocol), the injection
  control interacting with the r_perp cut (test design), and the margin-0 sampled template at 2.2 SEM. Details:
  NOTES.md "Iteration 5", `report/mock_validation.md`. Signal profile: `report/signal_profile.json` (50 % of the
  information inside r_perp = 15 Mpc/h, < 1 % below 3; 80 % at r_par < 5).
- **Review 5 done (2026-09-14, gpt-6-astra, `report/reviews/codex_review_5.md`): Stage A not accepted.** Its
  blockers were real: the sampled template audit was double-masked (the corrected margin-0 coefficient is
  0.994 +- 0.013, my "Poisson edge" story was wrong), the independent response prediction ignored the r_perp >= 3
  cut, and the protocol was undecidable (1 A bound vs t x SEM = 1.2 A at N = 40) with the validation seeds
  examined twice. NOTES.md "Iteration 6" has the reading.
- **Iteration 6 implemented (commit 7416f4b; GATES.md v6; ITERATION6.md)**: every review-5 fix (single mask,
  cut in the prediction, generator patch side `mock.patch_side_rad`, zq, provenance pinning of `campaign_config`
  + basis marker, manifest-only `condor/reprovenance.py`), fixed ensemble with **fresh seeds** (sparse 1000-1399,
  N = 400, roles `full` 1000-1039 / `core`; dense 2000-2019; dev 3000-3004; diagnostics seed 1000), deprojected
  bound 0.5 A chosen with N (user decision: N = 400, "a publishable method"), `required` flag on every row,
  noise-free injection expectation (`pairs.accumulate(true_positions=...)`) as the injection gate. Tests 64/64 in the campaign control.
- **Iteration-6 scale-1 campaign COMPLETE (2026-09-14, Perlmutter preempt, 49 jobs, ~4 h wall): 28/36 required
  rows pass, Stage B blocked.** Products on Perlmutter only (`mocks/iteration6/`, 275 GB; not copied to RACF);
  `report/mock_validation.{md,json}` and the figures are the acceptance table. The review-5 fixes hold (sampled
  template 1.008 +- 0.006, dense normalisation 1.001 +- 0.008, injection expectation 0.979). What fails is one
  thing seen from several sides (NOTES.md "Scale-1 campaign", iteration 6): the deprojected null is -0.41 +-
  0.18 A **with the response off** and the seed's own templates (fixed templates: zero); the linear response is
  deprojected correctly (paired +0.02 +- 0.08); the bias is the quasar template's correlation with the
  response-free forest pair products (matched 10.3 +- 1.6, CMB 0.33 +- 0.14, combined through the deprojection
  coefficient 0.066 -> -0.34), i.e. the non-squeezed forest-forest-template bispectrum, <m_p m_q T> with m the
  forest's long-mode component (exact b_q^2 xi_ms^2 for the lognormal template; Gaussian four-point terms give the
  15-20 % excess of the response over P(DCD-C)P^T, which IS deprojected). Shape test p = 0.002 with the same
  signature; the contaminant sits in the same (r_perp, r_par) bins as the signal. Real-universe analogue: quasar
  b_2, F2, forest non-linearity — generic, not modelled by the squeezed-limit deprojection of the report.
  Secondary: CMB-template normalisation 0.887 +- 0.048 vs truth 0.959 +- 0.019 on the same seeds.
- **Long-mode test done (96 seeds; NOTES.md "Long-mode test result"; `report/longmode_iteration6_summary.json`)**:
  the m-only field (m = c P[delta_L], the forest's long modes, squared against the seed's own templates; no
  lensing, no response) scores matched +11.5 +- 2.3 and deprojected -1.12 +- 0.29 A: <m m T> != 0 is the
  mechanism. But subtracting m from the forest removes 45 % of the lensing signal and 70 % of the response
  (the 10 Mpc/h long modes carry 19-75 % of xi at r_perp 3-20 Mpc/h): projecting the long modes out is not a
  fix; the non-squeezed term needs a model or its own nuisance template.
- **Review 6 done** (`report/reviews/codex_review_6.md`; NOTES.md "Review 6"): Stage A (CMB path) not passed.
- **PIVOT (user decision 2026-09-15): low-redshift tracers.** The CMB cross-correlation cannot separate lensing
  from the second-order forest-forest-density term with the modelling at hand; the (A_lens, A_2) plane is parked
  in `CMB_FUTURE_WORK.md`. Iteration 7 = cross-correlation of the same pair-template estimator with DESI LRG /
  ELG / QSO in five redshift slices (z = 0.4-1.75), bias per slice from the angular auto-correlation on
  k < 0.2 h/Mpc, per-slice A_L combined optimally, inverse-variance sightline weights; gate = unbiased detection
  against a lognormal redshift tracer (GATES.md v7). Code: `code/pipeline/lowz.py`, `run_lowz_validation.py`,
  `campaign7.py` (own provenance, seeds 4000-4399 / dev 5000-5004, no dense mocks), `slurm/campaign7_perlmutter.sh`,
  `code/stageb/lowz_catalogues.py` (DESI DR1 LSS catalogues -> per-slice HEALPix templates and biases; LRG/ELG
  catalogues and randoms are in `raw/desi/lss_v1.5/`), `tests/test_iteration7.py`. Findings and fixes in NOTES.md
  "Iteration 7" (lognormal smoothing; half-pixel binning bug; pixel window; spherical shot noise).
- **Iteration-7 scale-1 campaign COMPLETE (2026-09-15 13:00 UTC)**: 400 seeds (RACF Condor seeds 4000-4007 and the
  controls; Perlmutter seeds 4008-4399; collect on NERSC; products complete in `mocks/iteration7/` on NERSC, partial
  on RACF). `report/lowz_validation.{md,json}`: **24/30 required rows pass; every lensing gate passes** (combined
  null -0.05 +- 0.08 A, recovery 0.91 +- 0.08, paired response 0.966 +- 0.025, per-slice 0.94-1.02, fixed template
  0.01 +- 0.08, curl, covariance 0.96, slopes, spectra, injection 0.980). The failing rows are the six tracer-bias
  precision rows (1-3 % residuals against a 2-SEM = 0.4 % rule at N = 400) and QSO fallbacks in 2-5 % of seeds;
  reading in NOTES.md "Iteration 7" / "Scale-1 campaign". Per-seed scatter 1.56 A -> sigma(A) ~ 0.3-0.5 on DR1.
- **Stage B readiness checked on DR1 (2026-09-15; NOTES.md "Stage B readiness")**: data complete on RACF; the
  real-data templates for all five slices are built (`$LYALENSER_DATA/lowz_split/`, biases in
  `report/stageb/dr1_tracer_biases.json`), the delta reader and the spherical band templates exist, and the
  estimator runs end to end on real forests (`code/stageb/dry_run_lowz.py`; `report/stageb/dry_run_disc190_30_12.json`).
  DR1 precision ~1 A (limit, not detection).
- **First DR1 measurement done (2026-09-15, single slab; `report/stageb/dr1_lowz.{md,json}`; NOTES.md "DR1 result")**:
  A = 0.01 +- 0.67 (A < 1.1 at 95 %), curl and random-template nulls pass, injection bookkeeping 1.014. Run with
  `code/stageb/run_dr1_lowz.py` (24 min on a Perlmutter node; products in `LyaLenser_data/stageb/dr1_lowz/` on
  NERSC; the deltas, `lowz_split` templates and the DESI-pixel basis are on NERSC too). Tracer biases from the
  split-catalogue cross-spectra, cross-checked against ACT kappa (`report/stageb/cmb_bias_check.json`).
- **Dedicated report** `report/lowz/lowz.tex` (compile with `pdflatex; bibtex lowz; pdflatex x2` in `report/lowz/`;
  no latexmk on RACF); figures by `code/lowz_report_figures.py` (needs `$LYALENSER_DATA/lowz_split/summary.json`
  and `stageb/dr1_lowz/xi.h5`, both on RACF). Missing bibtex entries for 2503.14745, 2306.06312, 2405.16593:
  run `pip install adstex; ADS_API_TOKEN=... adstex lowz.tex -o lowz.bib` (skill `adstex-references`).
- **Iteration 8 (2026-09-16), the response-kernel fit** (`code/pipeline/NOTES.md` "Iteration 8"): the
  two-parameter Kaiser table is rejected by the DR1 counts and biases the amplitude at first order. The
  continuum projection is now the pair-weight average over real forest pairs (`xi_fit.project_fine_sample`,
  `code/stageb/build_basis_dr1.py` -> `$LYALENSER_DATA/stageb/basis_dr1_v2.h5`), and the table carries a spline
  correction plus a same-wavelength term that enters the mean field but not the kernel (`xi_spline.py`,
  `cfg.xi_correction='spline'`). chi^2 1880 -> 200 over 810 cells; the kernel change raises the amplitude by
  1.077 (the adopted bicubic correction, 22 parameters; more parameters gain nothing on the data and cost
  response on the mocks). **DR1 with it: A = -0.101 +- 0.754**, `$LYALENSER_DATA/stageb/dr1_lowz_v3`,
  `report/stageb/dr1_lowz_v3.{json,md}`, `slurm/dr1_lowz_v2.sbatch`. The injection expectation reads 1.038
  instead of 1.014 only because the injection displaces the non-lensable same-wavelength term (1.014 with it
  removed, `slurm/injection_sw.sbatch`). Mock cost of the correction: `report/xi_correction_mocks_variants.json`.
- **Iteration 9 (2026-09-16)**: region B (the Lyb window) used as an extension of every sightline (A x A and
  A x B pairs, B x B dropped), science bands extended to L = 500, slice cross-talk measured. **DR1:
  A = 0.182 +- 0.663**, `$LYALENSER_DATA/stageb/dr1_lowz_v4`, basis `stageb/basis_dr1_ab.h5`
  (`build_basis_dr1.py --regions lya lyb`, which samples geometry by quasar via `--modulus`). Details and the
  attribution between the two changes: `code/pipeline/NOTES.md` "Iteration 9".
- **Where things run now**: Perlmutter went down 2026-09-16 for a week of maintenance, so the workstation holds
  everything needed for Stage B: `raw/desi/lya-deltas/delta-{lya,lyb}-0-0`, `lowz_split/` (pulled from RACF),
  `stageb/`. A full DR1 run is ~40 min and 23.5 GB with `NUMBA_NUM_THREADS=22`. The iteration-7 mocks are only
  on Perlmutter (284 GB), so mock validation waits for it to come back.
- **NEXT**: mock-validate the five-band basis and, if a Lyb window is added to the generator, region B; debias
  the response matrix for the kernel-fit attenuation before the full sample; the joint slice fit; understand the
  r_perp tilt of the corrected projection (NOTES iteration 8, item 5); tomographic sub-slabs, the Planck
  cross-check, bias uncertainties into the error, GATES v8 for the record.
- Per-seed scatter 3.7 A on the 400 deg^2 mock is sigma(A) ~0.8 on DR1: S/N ~1 for A = 1, as forecast; DR1 gives
  a limit and the response measurement, detection needs full DESI (report Table 2).

## Absolute rule while the campaign runs
Every campaign product carries a fingerprint of all `code/**/*.py` (tests included), `GATES.md`, `Config`, the A
grid, the stream names, `report/numbers3.json`, the ACT mask and the package versions. **Changing any of these
invalidates every product on every site** (`completion()` raises "provenance mismatch"). Do not edit them until
`collect` has run; documentation files (`*.md`, `report/` outputs) are safe to edit. Products made on different
machines with the same fingerprint merge by copying directories.

## Where things live
- RACF: repo `/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser`; `LYALENSER_DATA=/gpfs/mnt/gpfs02/astro/workarea/anze/Data/LyaLenser`
  (`raw/` = all DR1 deltas, catalogues, ACT, Planck, 98 GB complete; `mocks/iteration3/` = round-3 products;
  `mocks/iteration4/` = the running campaign; `condor_logs/iteration4/campaign.out` = driver log);
  python `/gpfs/mnt/gpfs02/astro/workarea/anze/envs/lyalenser/bin/python` (numpy 1.26.4 pinned, see MEMORY.md).
- Perlmutter: nothing yet. Inputs tarball for it: `$LYALENSER_DATA/lyalenser_iter4_inputs.tgz` on RACF (1.4 GB:
  `mocks/iteration4/{dev,freeze}`, the ACT mask and N_L). Recipe: `slurm/README_PERLMUTTER.md`.
- Workstation where the project started: `/data/LyaLenser` (raw + mocks), Mathematica, Codex plugin; see MEMORY.md.

## If you are the session on Perlmutter
1. `git clone git@github.com:slosar/LyaLenser.git $SCRATCH/LyaLenser` (commit >= baeebbb), then follow
   `slurm/README_PERLMUTTER.md` sections 1-2 exactly: conda env from `slurm/requirements-perlmutter.txt` (the
   seven fingerprinted package versions must match RACF), unpack the inputs tarball into `$LYALENSER_DATA`
   (get it from RACF: `scp anze@astrosub02.sdcc.bnl.gov:/gpfs/mnt/gpfs02/astro/workarea/anze/Data/LyaLenser/lyalenser_iter4_inputs.tgz $SCRATCH/`
   or Globus; ask the user if you cannot reach RACF), check the mask sha256, run the 51 tests, and run the
   no-compute provenance check (`--phase freeze` must print `COMPLETE freeze`). If it raises, a version differs:
   fix the environment, never the code.
2. Decide the split with RACF: if you can list `mocks/iteration4/{sparse,dense,controls}/*/complete.json` on
   RACF, skip those seeds (`SPARSE_SEEDS=... DENSE_SEEDS=...`); if not, run everything — duplicates are harmless
   (same provenance, identical content), only wasted core-hours. The dense seeds (10 x ~40 core-h) are the part
   RACF will take longest to serve, so at minimum run those.
3. `export NERSC_ACCOUNT=<allocation> LYALENSER_DATA=... LYALENSER_REPO=... LYALENSER_PYTHON=...` and
   `cd slurm && ./campaign4_perlmutter.sh` (four dependent jobs; per-phase logs in `$LYALENSER_DATA/slurm_logs/iteration4/`).
   One CPU node each for sparse (~5 h) and dense (~4 h), then controls and collect (~1 h).
4. When `collect` has run (on whichever site holds the complete set; copy directories to complete it, README
   section 4), read `report/mock_validation.md`: the acceptance table with PASS/FAIL per gate and
   `Stage_B_allowed` in the JSON. Write the outcome into NOTES.md ("Scale-1 campaign" subsection) and
   PROGRESS.md, commit `report/mock_validation.{md,json}` and `report/figures/mock_iteration4_*.pdf`, push.
5. If gates fail: diagnose before changing anything (NOTES.md documents the tools — the exact pixel-position
   expectation for xi, the template audit variants, random-sightline controls); a code change means a new
   GATES.md, a new freeze and a fresh campaign (never tune on validation seeds). If all gates pass: astra review
   (when Codex works), then Stage B per `IMPLEMENTATION.md` Sec. 3 — DR1 deltas -> measured xi -> pair catalogue
   -> templates (ACT, Planck, matched from DR1 quasars) -> amplitudes, bands, curl, injections, 400 ACT random
   templates, jackknife; single slab first, then tomographic sub-slabs with disjoint selection. Expect ~1 sigma;
   the response term (~2 sigma) is the positive control. Stage B needs `raw/` (on RACF and the workstation).

## If you are the session on RACF
- The driver survives the session: `pgrep -f '[c]ampaign4.sh'`; resume with
  `cd condor && export LYALENSER_DATA=... && nohup ./campaign4.sh > $LYALENSER_DATA/condor_logs/iteration4/campaign.out 2>&1 &`
  (skips completed phases). `condor_q -nobatch`, `condor_q -analyze <id>` ("would match if drained" = waiting for
  memory), `condor_q -af HoldReason` if a job shows `H`. Requests: 32 GB sparse/dev, 36 GB dense/controls; do not
  raise them (48 GB never matches), do not lower below the measured peaks (23.5 GB dev, 26.9 GB sparse, ~28 GB dense).
- Products arriving from Perlmutter: rsync into `$LYALENSER_DATA/mocks/iteration4/`; the driver's own `submit`
  skips by condor log, so it may resubmit a seed that Perlmutter finished — that job is a no-op.

## Conventions and traps (details in MEMORY.md, IMPLEMENTATION.md Sec. 0, NOTES.md)
- Pair separation is theta_a - theta_b; observed field = true field at theta_obs + alpha; mocks sample at +alpha,
  injections shift by -alpha; kappa = -laplacian(phi)/2.
- Absolute statistics only (F^-1(q - mf)); a "PASS" whose SEM exceeds the quantity is not a pass; the smoothing
  width and every numerical choice are frozen on development seeds 100-104 and never changed on validation seeds.
- numpy must be < 2 (NEP 50 promotion breaks the float32 pair tests; `np.trapz` removed) — `requirements.txt`.
- HTCondor on RACF: executable and logs on gpfs (never /tmp), see MEMORY.md; `pkill -f` of a pattern in your own
  command line kills your shell.
- Codex on RACF needs `network_access = true` in `~/.codex/config.toml` (no network namespaces); the Codex sandbox
  cannot submit condor jobs; the spend cap can kill a task silently — check the job log for "spend cap".
- The DR1 deltas are `DELTA_BLIND` with `BLINDING=desi_y1`, but that blinding only affects picca's fiducial
  distance conversion; the arrays are the measured fluctuations.

## Key files
`report/main.tex|pdf`, `report/reviews/`, `report/iteration4_smoke2/` (smoke report, not acceptance),
`report/mock_validation.{md,json}` (round 3 until the campaign's collect overwrites it), `GATES.md`,
`IMPLEMENTATION.md`, `PLAN.md`, `code/pipeline/{campaign4,template_audit,random_streams,mock,run_mock_validation}.py`,
`code/pipeline/{ITERATION4,NOTES}.md`, `condor/{CAMPAIGN.md,campaign4.sh,phase.sub,run_phase.sh}`,
`slurm/{README_PERLMUTTER.md,campaign4_perlmutter.sh,phase_list.sbatch,requirements-perlmutter.txt}`,
`code/paths.py`, `code/fetch/fetch_data.sh`, `MEMORY.md`, `PROGRESS.md`.
