# Stage A, iteration 4 brief (draft, 2026-09-11): targeted fixes after the round-3 acceptance run

Round 3 (gpt-6-astra) produced an honest 20/30 with frozen gates (report/mock_validation.md, GATES.md,
report/reviews/impl_round3_summary.md). The reviewer's own diagnostics on the saved products
(/data/LyaLenser/mocks/iteration3/) localise the failures:

## Findings from the saved products
F1. **Matched template vs truth (map level).** Regressing the matched map on the saved `maps/kappa_slab`
    (exact 691->128 rebinning, inside the common mask) gives coefficients 0.454 / 0.664 / 0.866 for template margins
    0 / 150 / 300 Mpc/h (seed 0; 20-seed mean 0.69 +- 0.02 at margin 150). `maps/kappa_slab` is integrated over the
    full box (forest slab +- 300), so most of the deficit is definitional, but at equal range (margin 300) the
    template still recovers only 0.87 of the slab convergence. A 13% normalisation deficit of the template leaves 13%
    of the response contamination in the deprojected statistic. Origin unknown (candidates: effective bias of the
    lognormal sample vs the b_q used in the weights, RSD moving objects across the range edges, completeness/mean
    normalisation, mask edge handling).
F2. **Shared sightlines.** The validation used the shared selection (sightline quasars are in the template). All
    matched/deprojected failures have the signature of template-geometry coupling: matched raw mean field 480 +- 30
    (corrected 25.6 +- 8.9, 2.9 SEM), matched response-only excess 19 +- 9 over the prediction, deprojected null
    -2.4 +- 0.7, combined recovery -1.6 +- 0.7. The `disjoint_selection` flag exists in `generate_mock` but was not
    the validation baseline.
F3. **Normalisation at A_true = 1.** Sparse-ensemble truth-template absolute means: A0 -0.39 +- 0.31, A1 +0.42 +- 0.35,
    A5 4.42 +- 0.28, A10 9.54 +- 0.39 (N = 40). The slope gate (0.999) is carried by A = 10, where remapping is not
    linear (dense mock: A(1) = 0.877, A(5)/5 = 1.02, A(10)/10 = 1.07). The A_true in {5, 10} rule of ITERATION3.md was
    a mistake: displacements at A = 10 are comparable to pair separations. The A = 1 increment is 0.80 +- 0.46
    (sparse) and 0.88 (dense, single realisation): a 10-20% low response at the physical amplitude is not excluded.
F4. **Measured vs analytic xi derivative**: slope ratio 1.19 (dense). Either the analytic table or the measured
    smoothed derivative is off by ~20% in the response integral; on real data only the measured one exists.
F5. Covariance ratio 0.696 (gate 0.7-1.3): borderline; check after F2 (the matched/deprojected variance is inflated
    by the coupling).

## Required changes
1. **Map-level template gate (new).** For each seed, compute kappa_slab over exactly the template's chi range from
   the stored delta_m (or store it in the mock) and regress the matched map on it inside the mask: coefficient
   1 +- 0.03 with a Poisson-noise-free template (use the mock's continuous quasar intensity field) and 1 within SEM
   with the sampled quasars. Report the coefficient vs margin. This gate must pass before any deprojection gate is
   interpreted; find and fix the origin of the 0.87.
2. **Disjoint selection as the validation baseline.** Template quasars exclude sightline quasars (and their
   randoms are unaffected); the shared selection is run as a diagnostic on the same seeds. Report the shift.
3. **Normalisation gates at physical amplitude.** A_true in {0, 0.5, 1, 2} only. Dense noiseless variant
   (n_los = 100, P_N = 0) with >= 10 seeds at scale 1 for the percent-level normalisation gate (slope 1 +- 0.03 and
   absolute A(1) - A(0) = 1 +- 0.03); sparse noisy ensemble for consistency (2 SEM) and for the residual-bias gates.
4. **Resolve F4**: predict xi and its derivative on the discrete grid with the pixel window (the round-3
   `grid_covariance.py` machinery) and compare with the measured table on dense noiseless mocks; the residual
   between measured and predicted derivative over 5 < r_perp < 30 must be < 3% in the response integral, or the
   discrepancy must be explained and the measured-table procedure corrected (smoothing bias, binning, edge effects).
5. Re-run only what these changes require (development seeds for any new numerical choice, then the frozen
   ensemble); keep GATES.md discipline: update GATES.md before the run with the new gates and A grid.
6. Persist per-seed kappa_slab-over-template-range maps and the template coefficient so item 1 is reconstructible.

Everything else from ITERATION3.md stays as implemented. Do not modify IMPLEMENTATION.md, PLAN.md, report/main.tex,
mathematica/. Do not commit.
7. **Cluster entry points (for BNL RACF / HTCondor).** Add `--phase seed --seed N` (idempotent: skips if the seed's
   diagnostics JSON exists with matching provenance; writes only that seed's products) and `--phase collect`
   (assembles diagnostics, runs the stopping rule, extras and rebuild). Read the thread count from
   `NUMBA_NUM_THREADS` instead of hard-coding 24. Use `code/paths.py` / `LYALENSER_DATA` for every data path (no
   absolute `/data/LyaLenser` or `/home/anze` in code). See `condor/README_RACF.md`.
8. **Independent random streams.** The round-3 audit found that the density field and the CMB-noise generators
   are seeded with the same integer and therefore reuse the same normal-stream prefix. Derive every generator in a
   seed from one `numpy.random.SeedSequence(seed).spawn(k)` (field, forest noise, CMB noise, quasar sampling,
   random catalogues, ...), add a unit test that the first draws of any two streams differ, and record the stream
   layout in NOTES.md. This is a new numerical choice: it is frozen with the others before the ensemble.

## Execution model on BNL RACF (added 2026-09-11, the round is run from here)

The repository now lives at `/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser`; data at
`LYALENSER_DATA=/gpfs/mnt/gpfs02/astro/workarea/anze/Data/LyaLenser` (`raw/` and `mocks/` as before;
`mocks/iteration3/` holds the round-3 products). Python: `/gpfs/mnt/gpfs02/astro/workarea/anze/envs/lyalenser/bin/python`
(3.11; numpy 2.4, scipy 1.17, numba 0.67 — `np.trapz` was already replaced by `scipy.integrate.trapezoid`; do not
reintroduce numpy-1.x-only calls). Tests: `cd code/pipeline && python -m pytest -q tests` with `NUMBA_NUM_THREADS=4`.

**Hard constraints of the machine the implementation runs on**: the login node has 4 cores and 31 GB RAM, and the
Codex sandbox cannot submit HTCondor jobs. A scale-1 seed needs ~28 GB peak and 24 threads on the workstation, so
**no scale-1 mock may be generated in the implementation session**. Do all development at `scale <= 0.25`
(smoke level) with `NUMBA_NUM_THREADS=4`, and make every heavy phase a separately launchable, idempotent,
seed-parallel job that the Claude session submits to HTCondor (each job: 8 CPUs, 32 GB; 48 GB for dense variants;
`condor/seed.sub`, `condor/run_seed.sh`). Concretely the runner must expose:

- `--phase dev-seed --seed N` (N in 100-104) and `--phase freeze` (reads the five dev-seed products, applies the
  frozen choice rules, writes `development.json`/`frozen.json` and the hashes; refuses to run if any dev seed
  is missing);
- `--phase seed --seed N [--variant sparse|dense]` (sparse noisy: seeds 0-59 as in round 3; dense noiseless:
  n_los = 100, P_N = 0, A_true in {0, 0.5, 1, 2}, >= 10 seeds, its own seed range) and the map-level template gate
  products of item 1/6 written per seed;
- `--phase collect` (stopping rule, extras, controls that fit in one 8-CPU/48-GB job, `rebuild`, report,
  figures) — anything in `collect` that itself needs a scale-1 mock must instead be a `--phase control --name X`
  job so it can also go through condor.

Every phase reads `NUMBA_NUM_THREADS`; every product goes under `--mock-root` (default
`$LYALENSER_DATA/mocks/iteration4`); provenance (source/GATES/config hashes) is checked at `freeze`, `seed` and
`collect`, and a mismatch is an error, not a warning.

**Deliverable of this implementation session** (stop and report when done; do not attempt the campaign):
1. the code changes for items 1-8 with unit tests (36 existing tests still pass, new tests for the template gate,
   the disjoint baseline, the A grid, the RNG streams, the xi'-prediction comparison at smoke scale, and the
   idempotent seed phase);
2. `GATES.md` updated before any scale-1 run (new gates, new A grid, dense variant, unchanged tolerances
   unless justified in the file);
3. `condor/CAMPAIGN.md`: the exact ordered list of job commands (dev-seeds, freeze, seeds, dense seeds, controls,
   collect) with their resource requests, so the Claude session can submit them verbatim;
4. a smoke-scale end-to-end run (`scale 0.25`, 2 sparse + 2 dense seeds, dev seeds, freeze, collect) in
   `$LYALENSER_DATA/mocks/iteration4_smoke/` proving the phases chain, plus its `report/mock_validation.md`
   rendered from that smoke root (labelled as smoke, not acceptance);
5. NOTES.md iteration-4 section: what changed, the RNG layout, the origin of the 0.87 if found at smoke scale,
   and what remains to be checked at scale 1.

After the condor campaign the session is resumed for `collect`, the acceptance table and the notes.
