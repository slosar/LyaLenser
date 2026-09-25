# Stage A iteration 7: prospective acceptance protocol (low-redshift tracers)

Written before any iteration-7 mock is generated at scale 1. The runner (`code/pipeline/campaign7.py`) records this
file's SHA256, every `code/**/*.py`, the campaign configuration, the tracer table, the ACT mask, the seed ranges
and the package versions in every product; a changed input invalidates a product. GATES v6 (the CMB-lensing
protocol, iteration 6) is archived in git (commit 222f3a8) and in `mocks/iteration6/freeze/GATES.md`; its
outcome, 28/36 required rows with an additive deprojection bias from the second-order forest-forest-tracer term,
is why this iteration exists (`CMB_FUTURE_WORK.md`, `code/pipeline/NOTES.md` iteration 6).

## What is validated
The forest pair-template estimator cross-correlated with density tracers at z < 1.8 (DESI LRG, ELG and QSO in
redshift slices), which share no density modes with the forest at z = 2.1-3, so the correlation with the pair
products is lensing only. Per slice: a kernel-weighted tracer map (forest-source lensing kernel W(chi; chi_ref),
per-object weight 1/(b nbar)), the bias b of every tracer from its own masked angular auto-spectrum on large
scales (40 <= ell <= 0.2 chi(z_mid), shot noise subtracted, mask coupling and pixel window on the theory), the
Wiener combination of the slice's tracers with the model covariance (theory signal + measured shot noise), and
the sum over slices as the combined estimate of kappa_lya. Sightline pixels are inverse-variance weighted
(1/(sigma_N^2 + sigma_F^2); the pixel diagonal of C^-1). The user's gate (2026-09-15): an unbiased detection
against a lognormal redshift tracer.

## Ensemble and freeze
Development seeds 5000-5004 (chain check only; no numerical choice is made from them; seed 5000 also supplies the
other-realisation "fixed" templates). Sparse acceptance seeds 4000-4399 (N = 400): seeds 4000-4039 carry the A
grid {0, 0.5, 1, 2} and the response-off pair (A = 0, 1); every seed carries A_true in {0, 1} with the response
on. The null (A = 0) and recovery (A = 1) rows use the same 400 realisations and are correlated (declared; the
paired-difference row is the sample-variance-free normalisation). Seed 4000 supplies the injection and benchmark
controls. No dense mocks (the truth-template normalisation 1.001 +- 0.008 was established in iteration 6 on
noiseless mocks with the same forest code). Frozen numerical settings: as iteration 5-6 (fitted Kaiser-model
table, kernel r_perp in [3, 30], fit range 3 <= r_perp < 30, r_par < 30, A grid {0, 0.5, 1, 2}), plus the tracer
table of `lowz.TRACERS` (8 tracers in 5 slices, z = 0.4-1.75, DR1-like biases and densities, 2-D lognormal of
the projected density smoothed transversely by 3 Mpc/h), bias band k <= 0.2 h/Mpc, 128-pixel templates on the
generator's patch side, nearest-pixel binning about the pixel centres (the half-pixel offset of the flat-sky
catalogue binning was fixed on 2026-09-15 before this freeze), science bands 40-100-200-300.

Power: the per-seed scatter of the combined amplitude is not known before the dev seeds; the 0.5 A bound of v6
is retained with N = 400 (t x SEM ~ 0.37 A if the scatter is 3.7 A as for the deprojected estimator; a smaller
scatter is expected because no template is subtracted).

## Gates
| Gate | Requirement | Type |
|---|---|---|
| Bias from the angular auto-spectrum, every tracer (N = 400) | b_used / b_true within 2 SEM of 1; no fallback to the generator bias at scale 1 | required |
| Null (A_true = 0), own templates, every slice | within 2 SEM of 0 | required |
| Recovery (A_true = 1), own templates, every slice | within 2 SEM of 1 | required |
| Combined null (A_true = 0), own templates (N = 400) | within 2 SEM of 0; 95 % bound <= 0.5 A | required |
| Combined recovery (A_true = 1), own templates (N = 400) | within 2 SEM of 1; 95 % residual bound <= 0.5 A | required |
| Paired combined response A(1) - A(0), same realisation | within 2 SEM of 1 | required |
| Fixed (other-realisation) combined template, A_true = 1 | within 2 SEM of 0 | required |
| Curl null, combined template | within 2 SEM of 0 | required |
| Covariance: seed scatter / RMS jackknife, combined | in [0.7, 1.3] | required |
| Combined and truth normalisation slopes (A-grid seeds, response on) | within 2 SEM of 1 | required |
| Response on minus off, combined template, A_true = 0 (A-grid seeds) | within 2 SEM of 0 | required |
| Foreground and total spectra klkl, klkc, kckc | each within 10 % | required |
| Injection bookkeeping, noise-free expectation with the production selection (truth template) | odd slope at abs(A) = 0.25 within 1 +- 0.05 | required |
| Unit and regression tests; benchmark and memory | all pass; positive throughput, peak < 40 GB | required |
| Jackknife-covariance combination of the slice amplitudes vs the combined-template amplitude | report | report only |
| Paired truth-template response; injection convergence and measured slopes; tracer weights, shot noise, measured/model auto-spectra | report | report only |

A single failed required gate blocks Stage B on the low-redshift path. Failures are reported with their
diagnosis; no seed, tolerance or numerical choice changes in response to a validation result.

## Declared limits (Stage B design items, not gates)
The low-redshift tracers are 2-D lognormal Poisson samples of the projected slice density (independent slices,
Limber), without redshift errors, magnification bias of the tracers, or catalogue systematics; the forest mock
is as in iteration 6 (identical full-slab forests, uniform weights, flat patch, Cartesian light cone). DR1 has
varying forest lengths and weights, a curved footprint and spherical operators, the real n(z) and masks of the
LSS catalogues, and photometric systematics shared between the tracer catalogues and the quasar sample; Stage B
carries its own bias measurement (auto-spectrum on k <= 0.2, cross-checked against the kappa_CMB cross-spectrum),
injection tests and random-template nulls on the real geometry. The magnification of the sightline quasars by
the foreground is in the mock (it modulates the sightline density, absorbed by the mean field).
