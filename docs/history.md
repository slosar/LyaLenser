# How the analysis got here

A short chronology of the decisions behind the fiducial analysis (September 2026). The technical log of every
iteration, with all intermediate numbers, is archived verbatim in `legacy/mocks/NOTES.md`; the products of the
intermediate versions are in `results/` under their tags.

1. **CMB-lensing cross-correlation (iterations 1-6).** The project started as a cross-correlation of the forest
   pair estimator with the ACT/Planck convergence maps. The mock campaigns showed that a CMB template picks up
   the intrinsic forest-forest-density three-point function of the forest's own redshift range on top of
   lensing, and that this term cannot be filtered out. Parked; see `legacy/README.md`.
2. **Low-redshift tracers (iteration 7).** Templates from tracers at z < 1.8 share no density modes with the
   forest, so their correlation with the pair products is lensing only. The Stage-A mock campaign passed every
   lensing row (`legacy/mocks/report/lowz_validation.md`).
3. **Response-kernel fit (iteration 8).** The two-parameter Kaiser table under-fits the measured correlation by
   5-20 per cent over the kernel range; a bicubic spline correction and a same-wavelength term were added
   (`lyalenser.xi_spline`), fitted on the measured 1 Mpc/h cells with the basis projected through real forests.
4. **Region B, wider bands, cross-talk (iteration 9).** The Lyb-window deltas were added as an extension of the
   forest (matched on LOS_ID); the science bands went from L = 300 to 500; the slice-to-slice response leakage
   was measured (`slice_crosstalk.py`) and later removed by the joint fit.
5. **NaMaster, redshift evolution, BGS and BOSS, template validation (iteration 10).** Pseudo-C_ell with
   analytic mode coupling replaced the f_sky scaling; the correlation table became redshift-evolving (layered
   in chi, `lyalenser.xi_zevol`); BGS and BOSS joined the tracer set (0.1 < z < 0.4 and the BOSS footprint);
   the deflection templates were validated against the ACT and Planck maps with the exact prediction for the
   Wiener-filtered template (`deflection_cmb_check.py`).
6. **The whole forest below z = 3, source plane at z_eff (iteration 11).** 1.96 <= z <= 3.0, z_eff = 2.348 the
   weight-averaged pixel redshift; tracers restricted to z < 1.6 to keep a 400 Mpc/h buffer.
7. **Quasar-forest cross-correlation (iteration 12).** The same estimator on quasar-sightline pairs with the
   evolving cross model (`lyalenser.xi_cross`), nearly uncorrelated with the auto statistic; combined with the
   joint jackknife covariance.
8. **Bands to L = 1000 (iteration 13; fiducial by decision).** Five bands 40-200 ... 800-1000, deflections at
   nside 2048; L = 1300 kept as a convergence test (no precision gain).
9. **Measured-covariance Wiener weights, one joint response matrix, robustness (iteration 14).** The Wiener
   weights use the measured tracer auto/cross spectra per coverage class (`lowz_v4`); all 55 slice components
   are fitted at once per statistic (`joint_response_fit.py`, `results/joint_fit_v4.json`); injection tests on
   the fiducial products; robustness rows (bands, r_perp cut, knot sets) on NERSC and the workstation.
   Result: A = 0.61 +- 0.33; templates x CMB 0.99 +- 0.04 (ACT), 0.93 +- 0.03 (Planck).
10. **Derivative maps for the source-distance dependence (iteration 15, 2026-09-28).** Iterations 11-14 expanded
   the lensing efficiency to first order about the source plane with ONE coefficient for every template, an
   effective lens distance weighted by the forest x CMB kernel product (a remnant of the CMB thread; the true
   coefficient differs by an order of magnitude between the lowest and the highest slice). Every template now
   carries its own derivative map, the same tracers weighted by dW/dchi_s and Wiener-combined with the same
   weights (`lowz_catalogues.py` -> `lowz_v5`), and the fit contracts the pair accumulators with alpha and dalpha
   (`amplitude._partials`). The fiducial products were refitted from the saved catalogues (`dr1_lowz_v8`,
   `dr1_qso_v2`, `joint_fit_v5.json`), the robustness rows on NERSC (`slurm/refit.sbatch`), and a row without
   the source-distance term was added (`joint_fit_v5_noderiv.json`).

Decisions taken along the way that a reader might question:
- **Jackknife errors only, no mocks or Fisher errors in the paper.** The validation is the template x CMB
  comparison and the injection tests; the Stage-A mocks (legacy) validated the estimator's bookkeeping.
- **Curl templates stay in the fit** although dropping them changes nothing (`no_curl_*` entries of the joint fit).
- **A 2 Mpc/h inner cut is not worth it**: the information density (`scale_sensitivity.py`) shows about 2 per cent
  of the error to gain.
- **Cross-fade (partition-of-unity) band windows** would remove the notches between bands that the junk
  component now absorbs; offered, not done (a few per cent of S/N).
- **Per-slice source-plane rescaling** (pixels at z_f are lensed only by lenses in front of them): implemented in
  iteration 15 as the derivative maps (first order in the pixel's distance offset, exact lens distribution).
