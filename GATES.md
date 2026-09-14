# Stage A iteration 5: prospective acceptance protocol

Written before any iteration-5 mock is generated at scale 1. The runner records this file's SHA256, the source
SHA256s (every `code/**/*.py`), the configuration, the ACT mask and the package versions in every product; a
changed input invalidates a product. Iteration 4 (`report/mock_validation.md` at commit 261e07b, 39/48 gates) is
the reference this iteration is judged against; its diagnosis is in `code/pipeline/NOTES.md`.

## What changed since iteration 4
1. **The pair kernel's xi table is a fitted model, not a smoothed measurement.** The forest spectrum
   P_F = b_F^2 (1 + beta_F mu^2)^2 P_lin(k) F_NL(k, mu) is expanded in the three basis spectra P_lin F_NL mu^{2i};
   each is transformed to (xi_i, dxi_i/dr_perp) exactly as the pixels see it — for the mocks by the discrete-grid
   covariance with the trilinear pixel window and its aliases (the generator's own operators; `F_NL` is the
   generator's exp(-(k/k_p)^2)), for data by a Hankel transform with the line-of-sight pixel window — and pushed
   through the per-forest mean+slope continuum projection on the campaign's radial pixel grid, which commutes with
   d/dr_perp. Per sample, (b_F^2, beta_F) are fitted by weighted least squares of the 1 Mpc/h-binned model to the
   sample's measured 1 Mpc/h table over 3 <= r_perp < 30, r_par < 30 (weights = accumulated pair weights). The
   fitted table (xi and its analytic derivative) is the kernel table and the mean-field table. There is no
   smoothing width and no development-seed tuning: the development seeds only check that the chain runs.
   Phase `basis` computes the basis once per campaign; `freeze` records the fitted development parameters.
2. **Pair selection 3 <= r_perp <= 30 Mpc/h** (was 0-30). On the iteration-4 scale-1 products, r_perp < 3 carries
   < 1 % of the Fisher information of the estimator (`report/signal_profile.json`: 50 % within 15 Mpc/h, 75 %
   within 21; 80 % at r_par < 5) and it is where the table's small-scale shape matters most.
3. **Tolerances matched to the ensemble size** (user decision 2026-09-14; "5 % normalisation is acceptable if
   understood"): see the table. Report-only rows never block acceptance.
4. Everything else — mocks (radially smoothed lognormal quasars, sightline quasars behind the slab, disjoint
   selection baseline, RNG streams), templates, band basis, absolute statistics, jackknife, controls, A grid
   {0, 0.5, 1, 2}, seed ranges and the stopping rule — is as in iteration 4.

## Ensemble and freeze
Development seeds 100-104; recovery 0-19, null 20-39; if after the initial 40 either the deprojected null or the
combined recovery has a 95 % residual bound above 1 A, seeds 40-59 are added once to both ensembles; no other
extension. Dense noiseless seeds 200-209 (n_los = 100, P_N = 0, g on and off). Frozen numerical settings as in
iteration 4 except: no xi smoothing; kernel r_perp in [3, 30]; fit range 3 <= r_perp < 30, r_par < 30; forest model
`kaiser` (generator) for the mocks. Provenance and idempotency rules unchanged (`campaign4.py`).

## Gates
| Gate | Requirement | Type |
|---|---|---|
| Continuous template vs same-range truth, 40 <= L <= 300, margins 0/150/300 | mean 1 ± 0.03 (N = 60) | required, prerequisite of the deprojected gates |
| Sampled template vs same-range truth, same band and margins | mean within 2 SEM of 1 | required, prerequisite |
| Truth / CMB / matched normalisation slopes (sparse, response off) | within 2 SEM of 1 | required |
| Fixed vs refitted baseline | slope difference within 0.03 | required |
| Mean fields, varying and fixed templates, truth/CMB/matched | within 2 SEM of 0 | required |
| Absolute stochastic recovery (CMB, A = 1) | within 2 SEM of 1 | required |
| Response-only CMB and matched vs independent prediction | difference within 2 SEM | required |
| Response-only deprojected null | difference within 2 SEM; mean within 2 SEM of 0; 95 % bound <= 1 A | required |
| Combined deprojected recovery | within 2 SEM of 1; 95 % bound <= 1 A | required |
| Six-bin Hotelling shape | p > 0.01 | required |
| Covariance: seed scatter / RMS jackknife | in [0.7, 1.3] | required |
| Spectra klkl, klkc, kckc | each within 10 % | required |
| Dense normalisation slope, g on and g off (N = 10) | 1 ± 0.05 | required |
| Dense absolute endpoints | A0 within 2 SEM of 0, A1 within 2 SEM of 1 | required |
| Fitted-model vs generator-model normalisation (dense, same catalogues, kernel from the fitted vs the known table) | slope ratio within 5 % | required |
| First moments | omitted/correct slope ratio within 0.05 of the catalogue prediction | required |
| Injection bookkeeping | odd slope in [0.95, 1.05]; curl within jackknife | required |
| Analytic derivative convergence (Hankel provider, nk doubling); discrete-grid quadrature (128 -> 256 angles) | < 1 % | required |
| Unit and regression tests | all pass | required |
| Measured vs fitted table (coarse, r_perp^3 weight); fitted b_F^2 and beta_F vs the generator (sparse and dense) | report | report only |
| Template full-resolution coefficients; flags and margins; 100 random templates; benchmark and memory | report | report only |

A single failed required gate blocks Stage B. Failures are reported with their diagnosis; no seed, tolerance or
numerical choice changes in response to a validation result.
