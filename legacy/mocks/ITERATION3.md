# Stage A, iteration 3 (gpt-6-astra): make the acceptance run defensible

Authoritative input: report/reviews/codex_review_4_code.md (17 findings). The pair-compression core, the sign
conventions, the 11 accumulators, the jackknife and the outside-box spectra were verified there and must not change
in behaviour (keep the 17 unit tests green and add the new ones below). Everything else in that review is to be fixed.
Spec: IMPLEMENTATION.md (Sec. 0-2). Do not modify IMPLEMENTATION.md, PLAN.md, report/main.tex, mathematica/.
Write large products under /data/LyaLenser/mocks/ (writable). Do not commit.

## Required changes (numbers refer to the review findings)
1. **Absolute acceptance statistics (finding 1).** Every acceptance row must report the estimator as specified
   (score minus its own mean-field subtraction, divided by the response matrix), not a difference between two
   realisations. Paired differences may be reported in addition, clearly labelled "diagnostic".
2. **Complete band basis (finding 2).** Science bands [40,100], [100,200], [200,300] with tapers, plus a junk
   template, must satisfy sum_b w_b(L) + w_junk(L) = 1 for every L up to the grid Nyquist (junk = 1 - sum of science
   windows, including the taper wings). Add a unit test: a unit unfiltered signal at several L inside and across
   tapers, plus outside-band nuisance modes, is recovered as 1 in the corresponding band(s) with the junk at 1.
   `amplitude()` must reject fits missing any required component or with a singular response matrix.
3. **Baseline and calibration independence (finding 3).** Freeze the xi_from_data smoothing width and all numerical
   choices on a *development* seed set (seeds 100-104), document them, then validate on fresh seeds. Measure
   xi_from_data on each validation sample actually used (lensed, noisy, sparse), not on one noiseless dense unlensed
   patch, and report the fixed-vs-refitted baseline difference in A. Derive the 0.903 discrete-grid/pixel-window
   factor from the grid and pixel windows (state the formula) instead of matching it empirically, or drop the factor
   and let xi_from_data carry the normalisation. The normalisation patch must contain the 40-100 band (use the full
   scale-1 patch for the dense noiseless variant, with fewer seeds if needed).
4. **Numerical convergence at the amplitude level (finding 4).** Use the specified xi step of 0.25 Mpc/h (not 0.5)
   and demonstrate < 0.5% change in A when halving it; set the analytic provider's nk so that the derivative norm
   over 10 <= r_perp <= 30 changes by < 1% when doubled; test both in pytest on a small case.
5. **Continuum projection (finding 5):** joint weighted fit of intercept and slope per skewer (or slope coordinate
   centred at the weighted mean chi); unit test that both weighted residual moments vanish and that the projection is
   idempotent.
6. **Periodic interpolation (finding 11):** pad periodic endpoints (or implement periodic trilinear); unit test on
   [0,1,2,3] at 3.5 -> 1.5.
7. **Coordinate consistency (finding 10):** one flat-patch light-cone model for forest sampling (x = chi theta),
   quasar positions, convergence projections and pair separations, with cos(dec_mid) used consistently; test cross-
   field responses at several patch positions.
8. **RSD velocities (finding 12):** solve from the full-resolution Fourier modes (or low-pass before downsampling);
   compare the displacement spectrum with the full solution.
9. **Matched template (finding 7):** shot noise N_L^{ss} = (1+r) int dchi W^2/(b_q^2 nbar_3D chi^2) (reduces to
   (1+r) D^2 <W^2>/(b_q^2 n_2D) for uniform n; verify numerically); completeness normalised over the footprint, not
   over occupied random pixels; flat and spherical paths must agree on a test; matched mask imposed as a common map
   operator with the CMB mask.
10. **Wiener operator on the maps actually used (finding 6):** C_XX from the masked, noisy, pixelised map itself;
    S_L with the mock's actual template density, magnification and selection; independent random templates masked
    identically.
11. **Independent response prediction (finding 8):** pixel-level modulation delta_F -> (1 + (R_delta/2) delta_L(pixel))
    propagated through the continuum operator and the pair response (build the prediction from the stored delta_L
    map at every pixel, not one value per sightline), with the actual smoothing and map transfer; the deprojected
    prediction is computed, not set to zero.
12. **Mean-field protocol (finding 9):** fixed validation ensemble declared in advance (seeds 0-19 for recovery,
    seeds 20-39 for nulls); no post-hoc seeds; raw and mean-field-corrected results on identical seeds; a fixed-
    template conditional mean-field test (same template on all seeds) in addition to the varying-template one.
13. **Sub-slab selection (finding 15):** pixel slab labels honoured in both pair kernels, pixel pairs assigned by
    midpoint slab, one catalogue group per slab; test boundary-crossing pairs and a disjoint template/sightline
    selection.
14. **Poisson intensity (finding 13):** total counts drawn from the integrated realisation-dependent intensity.
15. **Margin test (finding 14):** vary the template margin within one larger correlated realisation (generate the
    density box with the largest margin once; build templates from sub-ranges), so that the test changes only the
    selection.
16. **Tests (finding 16):** FFT deflection tested through the implementation, band completeness, convergence, end-
    to-end matched template, continuum projection, RSD, and the acceptance statistics on a synthetic case;
    covariance-aware shape statistic (Hotelling) in the runner; injection fitted with the full seven-component basis.
17. **Persistence (finding 17):** save per-region q, F, mf partial sums, template operators/potentials, measured-xi
    tables, per-seed diagnostics and configuration provenance so every acceptance row can be rebuilt from disk.

## Precision of the acceptance gates
Per-mock signal-to-noise on a 400 deg^2 patch with the noisy template is small, so "1 within 2 SEM" with SEM ~ 1 is
not informative. For the normalisation gates use the linearity of the estimator: run the physical remapping at
A_true in {0, 1, 5, 10} on the fiducial sparse noisy mocks with the truth template and require the fitted slope to be
1 +- 0.03 (the error on the slope scales as 1/A_true); do the same with the filtered kappa_CMB and matched templates
(A_true scaling kappa_lya only). For residual-bias gates (response-only null, deprojected recovery) report the
absolute value with its SEM over the fixed ensemble AND the achieved bound, e.g. "|residual| < 0.4 A at 95%"; if the
bound is weaker than 0.3 A, increase the ensemble (seeds up to 59) rather than declaring a pass. State every gate
and tolerance before running the fresh ensemble, in a `GATES.md` file committed to the working tree before the run.

## Deliverables
Updated modules and tests; `GATES.md`; `report/mock_validation.md/json` rebuilt from the fresh ensemble with absolute
statistics; figures; iteration-3 section of NOTES.md with wall times, peak memory, every remaining failure and its
diagnosis. Print a final summary with the acceptance table.
