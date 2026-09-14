# Stage A iteration 6: prospective acceptance protocol

Written before any iteration-6 mock is generated at scale 1. The runner records this file's SHA256, the source
SHA256s (every `code/**/*.py`, recursively), the campaign configuration (`campaign_config`, i.e. the frozen choices
actually used), the ACT mask, the seed ranges, the package versions and the basis/freeze markers in every product;
a changed input invalidates a product. Iteration 5 (`report/mock_validation.md` at commit f040dd1, 34/38 required
rows by the v5 count, 31/35 by the literal required/report-only classification) and its adversarial review
(`report/reviews/codex_review_5.md`) are what this iteration answers; the diagnoses are in `code/pipeline/NOTES.md`.

## What changed since iteration 5 (all from review 5; no numerical choice of the estimator changed)
1. **Bug fixes.** (a) The sampled map-level template coefficient was computed on an already masked map and masked
   again inside the band regression (`template_audit.audit_mock`): the unmasked matched map now enters and the
   common mask is applied exactly once. (b) The independent response prediction (`response._predict`) accepted
   every pixel pair with r_perp <= 30 while the production kernel applies 3 <= r_perp <= 30: the prediction now
   applies the production selection. (c) Every map-level operator (band templates, matched-template pixelisation,
   audits, spectra, random templates) uses the generator's patch side nx dx / chi_ref instead of the nominal
   20 x scale degrees (nx is rounded up, 0.09 % at scale 1), so the deflection decomposition and the generator's
   deflection share one angular scale. (d) `SightlineSet.zq` is the observed redshift of each sightline quasar.
2. **Protocol.** One fixed ensemble, no data-dependent extension or stopping. Every acceptance row carries a
   `required` flag; `Stage_B_allowed` is the conjunction of the required rows only. The deprojected bound is set
   together with the ensemble size (below). The realisations of iterations 4-5 (seeds 0-59, 100-104, 200-209) are
   development material and are refused by the seed phase.
3. **Injection control.** The coordinate injection is judged on its noise-free expectation (delta_p delta_q
   replaced by the fitted xi at the true separation, production selection, exact per-sightline deflection against
   the 7-component band basis; `pairs.accumulate(true_positions=...)`), which must have odd slope 1 +- 0.05 at the
   smallest injected amplitude |A| = 0.25; its convergence with amplitude and the measured (single-realisation,
   noisy) slopes are reported. The measured slope alone is not a gate: its realisation scatter is ~0.3 (per-seed
   scatter of the physical truth slope in iteration 5).
4. **Provenance.** `campaign_config` is what the fingerprint records; the basis marker digest is pinned into the
   freeze and into every downstream product; `condor/reprovenance.py` requires an explicit migration manifest,
   verifies every artifact before rewriting and keeps the original provenance.
5. **Jackknife bookkeeping.** The midpoint HEALPix jackknife runs at nside 8 and refines to nside 16 when fewer
   than 30 regions are populated (the 20-degree patch: ~27 regions at nside 16); nside and the populated region
   count are recorded per fit and reported with the covariance row. DR1 region sizing is a Stage B design item.

## Ensemble and freeze
Development seeds 3000-3004 (chain check only; no numerical choice is made from them). Sparse acceptance seeds
1000-1399 (N = 400): seeds 1000-1039 carry the physical A grid {0, 0.5, 1, 2} (slope rows), every seed carries
A_true in {0, 1} with the response on and off; the deprojected null (A = 0) and recovery (A = 1) rows therefore
use the same 400 realisations and are correlated (declared here; the paired difference row is the
sample-variance-free normalisation). Seed 1000 also carries the margin, 100-random-template, injection, flag and
benchmark diagnostics. Dense noiseless seeds 2000-2019 (n_los = 100, P_N = 0, g on and off). Frozen numerical
settings as in iteration 5: fitted Kaiser-model table (no smoothing), kernel r_perp in [3, 30], fit range
3 <= r_perp < 30, r_par < 30, A grid {0, 0.5, 1, 2}.

Power: the per-seed scatter of the deprojected amplitude is ~3.7 A at scale 1 (iteration 5), so N = 400 gives
SEM ~0.19 A and t x SEM ~0.37 A; a 0.5 A bound leaves ~0.13 A for the residual. Scaled to the DR1 footprint the
per-seed scatter corresponds to sigma(A) ~0.8, so 0.5 A is ~0.6 of the DR1 statistical error. Rationale
(user decision 2026-09-14): a method that is publishable needs the deprojection bias resolved well below the
DR1 error bar; N = 400 costs ~40 node-hours on Perlmutter.

## Gates
| Gate | Requirement | Type |
|---|---|---|
| Continuous template vs same-range truth, 40 <= L <= 300, margins 0/150/300 | mean 1 +- 0.03 (N = 400) | required, prerequisite of the deprojected gates |
| Sampled template vs same-range truth, same band and margins (single mask) | mean within 2 SEM of 1 | required, prerequisite |
| Truth / CMB / matched normalisation slopes (sparse A-grid seeds, response off) | within 2 SEM of 1 | required |
| Fixed vs refitted baseline (A-grid seeds) | slope difference within 0.03 | required |
| Mean fields, varying and fixed templates, truth/CMB/matched | within 2 SEM of 0 | required |
| Absolute stochastic recovery (CMB, A = 1) | within 2 SEM of 1 | required |
| Response-only CMB and matched vs independent prediction (production selection) | difference within 2 SEM | required |
| Response-only deprojected null (N = 400) | difference within 2 SEM; mean within 2 SEM of 0; 95 % bound <= 0.5 A | required |
| Combined deprojected recovery (N = 400, same seeds) | within 2 SEM of 1; 95 % bound <= 0.5 A | required |
| Paired deprojected response A(1) - A(0), same realisation | within 2 SEM of 1 | required |
| Six-bin Hotelling shape | p > 0.01 | required |
| Covariance: seed scatter / RMS jackknife (nside and region count reported) | in [0.7, 1.3] | required |
| Spectra klkl, klkc, kckc | each within 10 % | required |
| Dense normalisation slope, g on and g off (N = 20) | 1 +- 0.05 | required |
| Dense absolute endpoints | A0 within 2 SEM of 0, A1 within 2 SEM of 1 | required |
| Fitted-model vs generator-model normalisation (dense, same catalogues) | slope ratio within 5 % | required |
| First moments | omitted/correct slope ratio within 0.05 of the catalogue prediction | required |
| Injection bookkeeping, noise-free expectation with the production selection | odd slope at |A| = 0.25 within 1 +- 0.05 | required |
| Analytic derivative convergence (Hankel provider, nk doubling); discrete-grid quadrature (128 -> 256 angles) | < 1 % | required |
| Unit and regression tests | all pass | required |
| Injection convergence with amplitude; measured injection slopes and their difference from the expectation; curl | report | report only |
| Measured vs fitted table (coarse, r_perp^3 weight); fitted b_F^2 and beta_F vs the generator (sparse and dense) | report | report only |
| Template full-resolution coefficients; flags and margins; 100 random templates; benchmark and memory | report | report only |

A single failed required gate blocks Stage B. Failures are reported with their diagnosis; no seed, tolerance or
numerical choice changes in response to a validation result. Report-only rows never block acceptance.

## Declared limits of this validation (Stage B design items, not Stage A gates)
The mocks are flat-sky FFT patches with identical full-slab forests, uniform weights and a Cartesian light cone;
DR1 has varying forest lengths, noise-dependent weights, a curved footprint and spherical operators. The
model-shaped kernel's projection is a geometry-averaged approximation (report Sec. 5). The mock normalisation
(dense slope) is not transferred to DR1 as a scalar; Stage B carries its own response measurement and injection
tests on the real geometry.
