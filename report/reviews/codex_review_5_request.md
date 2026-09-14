# Review request 5 (2026-09-14): adversarial review of Stage A rounds 4 and 5

Model: gpt-6-astra, effort high, read-only. Output saved verbatim as `report/reviews/codex_review_5.md`.

## Context to read first
`HANDOFF.md`; `PROGRESS.md` (2026-09-11 onwards); `code/pipeline/NOTES.md` sections "Iteration 4" and
"Iteration 5" (the diagnoses: template deficit, xi' excess, provenance fix, normalisation, signal profile);
`GATES.md` (v5, the frozen protocol); `code/pipeline/ITERATION4.md`, `ITERATION5.md`; `report/mock_validation.md`
(iteration-5 acceptance table, 34/38 required gates) and `report/mock_validation.json`; `report/signal_profile.json`;
`IMPLEMENTATION.md` Sec. 0-1 (estimator conventions) and `report/main.pdf` Sec. 5 (pair-template estimator).

## Code under review (all written since the round-3 review `codex_review_4_code.md`)
- `code/pipeline/xi_fit.py` (new, round 5): basis spectra, pixel operators, continuum projection, (b_F^2, beta_F)
  fit, the kernel table. `grid_covariance.py` (`xi_from_mock_grid(power=...)`, `sample_covariance`).
- `code/pipeline/mock.py` (round 4-5 changes): radially smoothed lognormal quasars (`quasar_radial_smoothing`),
  sightline quasars drawn behind the slab (`sightline_proximity`, 1.6x oversampling), node-centred cell sampling,
  `grid_geometry`, RNG streams (`random_streams.py`).
- `code/pipeline/template_audit.py` (map-level template gate, band-limited `band_regression`).
- `code/pipeline/campaign4.py` (phases basis/dev-seed/freeze/seed/control/collect, provenance fingerprint,
  re-provenance in `condor/reprovenance.py`), `run_mock_validation.py` (`table_for`, `process_seed`, `build_rows`).
- `code/pipeline/pairs.py` (`r_perp_min`), `config.py`, `xi_model.py` (even boundary; `xi_from_data` now only
  supplies counts), `validation_stats.py`, `signal_profile.py`, tests `tests/test_iteration4.py`, `test_iteration5.py`.

## Questions (answer each with a verdict and the evidence in the code or the products)
1. Is the model-shaped table correct as the estimator's response kernel: the basis expansion in mu^{2i}, the pixel
   operators (grid B3 window with aliases for mocks; Hankel + LOS window for data), the claim that the per-forest
   mean+slope projection commutes with d/dr_perp and is applied consistently to xi (mean field) and xi' (response),
   the weighted least squares (weights = accumulated pair weights) and its fit range, and the use of the same fitted
   table for the mean field mf? Any hidden dependence that biases A?
2. The dense normalisation is 0.980 +- 0.014, linear in A, and identical (0.988) with the generator's own table.
   Is the remaining -2 % explained by something in the code (flat-patch FFT band templates, chi cos(theta) depths,
   cos(dec_mid) pair metric, trilinear resampling in the lensed mock, junk/curl basis), and does it matter for
   Stage B on DR1 where the footprint is not a 20-degree patch?
3. The fitted (b_F^2, beta_F) are offset +10 %/-9 % from the generator on every sample (SEM 1 %) while the
   cell-by-cell data/model ratio is 1.00 +- 0.01 at r_par < 3 and within 1-2 % averaged over r_par at every r_perp;
   NOTES attributes the r_par pattern to the discrete 0.55 Mpc/h pixel lags in the 1 Mpc/h bins. Is that right, is
   the 2-parameter fit degenerate as claimed, and is the kernel really insensitive (0.8 %)? Would a 3-parameter or
   lag-aware fit be needed for the data (Arinyo terms, LOS resolution)?
4. Mocks: are the round-4 changes physically justified and free of new artefacts — radial 8 Mpc/h smoothing of the
   lognormal input (b_q sigma = 1.3), sightline quasars behind the slab (forest in front of the quasar), disjoint
   selection with the shared diagnostic, node-centred sampling? Any remaining coupling between the sightline
   positions and the forest field that biases the measured xi or the template?
5. Template gate: continuous 1.004 +- 0.002 in 40 <= L <= 300 but the sampled margin-0 coefficient is 0.972 +- 0.013 in
   both iterations. Is the band-limited masked-map coefficient the right statistic, and what produces the margin-0
   deficit (Poisson edge effect with RSD, random-catalogue normalisation, completeness)?
6. Protocol: the deprojected rows fail only on a 1 A 95 % bound that t x SEM alone exceeds at N = 40, and the
   injection control shifts positions across the r_perp >= 3 cut. Judge the proposed corrections for iteration 6
   (bound 1.5 A or t x SEM + 0.3 A; injection control with r_perp_min = 0) and whether they, or anything else,
   would constitute tuning on validation results. Is the stopping rule / extension handling sound?
7. Provenance: fingerprint contents, the documented re-provenance of iteration-4 products, cross-site merging,
   idempotent phases. Any hole that could let a changed input go unnoticed?
8. Errors: midpoint jackknife (nside 8), covariance ratio 0.96, sigma_F reported only as a scale. Adequate for
   Stage B on DR1 (jackknife regions vs the DR1 footprint; the 400 ACT simulations; the 100 random templates)?
9. Anything that should block Stage B beyond the four failing rows, and the shortest path to an acceptable Stage A.

Format: numbered findings with severity (blocker / significant / minor), file and line references, and a final
verdict on Stage A readiness. Do not modify any file.
