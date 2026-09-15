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
- **NEXT (user decision)**: accept the lensing gate as met and write GATES v8 with an equivalence bound on the
  bias (5 %) for the record, or rerun with a changed rule (no code change needed); then Stage B on DR1: sightline
  set from the deltas on the sphere, `code/stageb/lowz_catalogues.py` templates (LRG done for 0.4-0.6: b = 1.68
  +- 0.03; run `--tracers LRG ELG QSO` for all slices), deflections via `templates.alpha_at`, per-slice and combined
  A_L with jackknife, injection expectation on the real geometry, bias cross-check against kappa_CMB x tracer.
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
