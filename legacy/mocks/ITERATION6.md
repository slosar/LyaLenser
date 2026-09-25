# Stage A, iteration 6 brief (2026-09-14; implemented by the Claude session after review 5, to be reviewed by gpt-6-astra)

Starting point: iteration 5 at scale 1 (`report/mock_validation.md` @ f040dd1) and the gpt-6-astra review
`report/reviews/codex_review_5.md`: Stage A not accepted. Blockers: the sampled template audit was double-masked
(corrected margin-0 coefficient 0.994 +- 0.013, so the "Poisson edge effect" diagnosis was wrong); the
independent response prediction ignored the r_perp >= 3 cut; the 1 A bound was undecidable at N = 40 and the
validation seeds had been examined twice. Significant: generator vs template pixel scale (nx rounded up), fitted
minus generator slope -0.0079 +- 0.0024 (a resolved, small table dependence), provenance holes, jackknife nside
bookkeeping, `zq`, mock-realism caveats for DR1.

User decisions (2026-09-14): bite the bullet, N = 400 fresh seeds ("sub-20 % error, a publishable method");
tolerance chosen with N before the run; then the astra review.

## Changes
1. Bug fixes: `template_audit.audit_mock` (unmasked sampled map into `band_regression`; `make_bundles` exports
   `matched_map_unmasked` and `side_rad`), `response._predict` (production selection r_perp_min/max, r_par_max),
   `mock.patch_side_rad` used by every map-level operator (bundles, dense templates, audits, spectra, random
   templates), `SightlineSet.zq` = observed quasar redshift.
2. Protocol (`campaign4.py`, GATES v6): ITERATION 6; seeds dev 3000-3004, sparse 1000-1399 (roles `full`
   1000-1039 with the A grid, `core` otherwise; every seed A in {0, 1} with response on/off), dense 2000-2019,
   diagnostics seed 1000; iteration 4-5 seeds refused; no extension/stopping; `required` flag on every row
   (`run_mock_validation.record`), `Stage_B_allowed` from required rows; deprojected bound 0.5 A
   (`DEPROJECTED_BOUND`); new paired-response row; covariance row reports jackknife nside/regions.
3. Injection control: `pairs.accumulate(true_positions=...)` expectation mode (delta delta -> xi at the true
   separation), `inject.injection_test(expectation=True)` and per-amplitude odd slopes; amplitudes
   +-0.25, +-0.5, +-1, +-2; required gate on the expectation at |A| = 0.25; measured slopes report-only.
4. Provenance: fingerprint records `campaign_config`, recursive `code/**/*.py`, iteration, seed ranges, bound;
   basis marker pinned in the freeze and in every downstream provenance; `condor/reprovenance.py` needs a manifest
   and verifies artifacts first.
5. Drivers: `slurm/campaign4_perlmutter.sh` and `condor/campaign4.sh` default to iteration 6 (seed ranges,
   DEV_SEEDS, EXTRAS_SEED, job names ly6-*).
6. Tests: `tests/test_iteration6.py` (prediction selection, expectation mode, paired slopes, single masking,
   patch side, required flag, fingerprint contents); `test_iteration4.py` updated for the seed ranges.

## Verification before the campaign
Unit tests 57/57 (RACF); smoke chain at scale 0.25 on Perlmutter (`CAMPAIGN_NAME=iteration6_smoke`, sparse
1000-1001, dense 2000-2001) end to end; then the scale-1 campaign `iteration6`.

## Deliverables
GATES v6 frozen, scale-1 campaign (400 + 20 seeds), `report/mock_validation.md`, NOTES.md iteration-6 section,
then the astra review of rounds 6.
