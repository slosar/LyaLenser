# Codex review 4 (gpt-6-astra, high): Stage A code + mock validation, after implementation round 2 — 2026-09-11

Verbatim output follows.

---

**Stage A has not genuinely passed. A third implementation round is required before Stage B.** The final JSON’s PASS flags reproduce its implemented predicates, but several predicates test a different statistic from the specification. Independent checks also expose numerical and mock-construction defects.

I read all eight modules, all tests, IMPLEMENTATION.md Sections 0–2, the relevant report sections, both implementation notes, and reviews 2 and 3. No files were modified; `git status --short` remains empty.

All **17 unit tests pass**. I disabled file-writing caches and used the existing Matplotlib cache for this read-only run.

The central algebra is substantially correct:

- The pair direction is \(\theta_a-\theta_b\), with the correct convergence-response sign.
- All 11 accumulators and their contractions into \(q\), \(F\), and \(mf\) match the explicit formulas in Section 1.2, including the symmetric cross-template term. The omitted \(g_1^2\,dm\,dc\) term is a specified approximation, not an implementation omission.
- The spherical deflection uses the already-divided healpy azimuthal derivative correctly. The flat-sky FFT sign and amplitude also pass an independent cosine test.
- The positive displacement in fixed-geometry mock sampling and negative coordinate-injection displacement have the correct signs.
- Jackknife regions really are assigned by spherical pair midpoint. The nside-16 fallback follows ITERATION2.md.
- The outside-box Gaussian pair implements the three restricted Limber spectra correctly for the fiducial slab.
- An independent analytic Gaussian-power-spectrum check verifies the correlation-transform normalization to \(8\times10^{-8}\) relative accuracy.

The findings follow, ordered roughly by importance.

1. **Blocker — “Combined recovery” removes the contamination it is supposed to test.**  
   [run_mock_validation.py:231](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:231), [run_mock_validation.py:240](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:240).

   The reported combined statistic is
   \[
   \hat A(A_{\rm true}=1,R_{\delta}=2)-\hat A(A_{\rm true}=0,R_{\delta}=2),
   \]
   using the same realization. This cancels additive response leakage and baseline errors, even when deprojection fails. Stochastic normalization, response-only tests, shape tests, and covariance also use paired differences. These are useful response diagnostics, but they do not replace the specified absolute-amplitude tests.

   I rebuilt the template bundles and fitted all 20 saved combined catalogues. The **absolute deprojected result is \(-0.602\pm3.067\) SEM**, rather than the published increment \(1.380\pm0.789\). Its scatter/RMS-jackknife ratio is 1.002; the published 0.865 describes the differenced statistic.

   **Fix:** report and gate absolute estimates with only the estimator’s prescribed mean-field subtraction. Keep paired differences as separately labelled diagnostics.

2. **Blocker — Tapered bands do not span the physical signal; a unit signal can return exactly 2.**  
   [templates.py:62](/home/anze/Dropbox/work/LyaLenser/code/pipeline/templates.py:62), [templates.py:83](/home/anze/Dropbox/work/LyaLenser/code/pipeline/templates.py:83).

   Science windows taper down inside their nominal ranges, while junk covers only \(L<40\) or \(L>300\). Consequently,
   \[
   \sum_b w_b(L)+w_{\rm junk}(L)\ne1.
   \]
   The missing taper contribution remains outside the fitted response model. A full \(F_{bc}\) cannot correct a missing basis component.

   My independent noiseless synthetic test used a unit unfiltered \(L=95\) signal plus an outside-band nuisance mode. The fitted science amplitude was **2.000000**, because its science template has half amplitude; junk fitted 1 correctly.

   **Fix:** construct complementary windows spanning the full conditional-mean field, and test unit recovery through taper transitions. Reject missing or degenerate required components; checking only `kind=="junk"` is insufficient.

3. **Blocker — The production baseline procedure is not validated, and the calibration evidence is circular.**  
   [run_mock_validation.py:175](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:175), [run_mock_validation.py:181](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:181), [run_mock_validation.py:223](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:223), [NOTES.md:22](/home/anze/Dropbox/work/LyaLenser/code/pipeline/NOTES.md:22).

   One measured-\(\xi\) table from the **unlensed seed-0, scale-0.15, noiseless dense mock** is reused for every amplitude and every full-size seed. Thus the “baseline absorption” test never measures the response absorbed by fitting \(\xi\) from the actual lensed sample.

   NOTES explicitly says the smoothing width was selected using normalization. The forest factor 0.903 is also an empirical matching factor: it changes covariance by \(0.903^2=0.8154\). Agreement with the analytic model is therefore not independent validation.

   Moreover, the normalization patch has **no modes in the 40–100 science band**. Its global slope 1.00695 hides successive interval slopes **0.8851, 0.8661, and 1.1326**.

   **Fix:** freeze choices using a separate development set; validate on fresh seeds with populated science bands. Refit production \(\xi\) on each relevant sample, explicitly compare fixed versus refitted baselines, and derive the discrete-grid/window prediction instead of matching it empirically.

4. **Significant — The required amplitude grid-convergence test actually fails.**  
   [xi_model.py:72](/home/anze/Dropbox/work/LyaLenser/code/pipeline/xi_model.py:72), [run_mock_validation.py:170](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:170), [test_xi_model.py:19](/home/anze/Dropbox/work/LyaLenser/code/pipeline/tests/test_xi_model.py:19).

   Reusing identical measured coarse counts and the saved normalization mock, I obtained:

   | \(\xi\) step, Mpc/h | Absolute \(A\), truth \(A_{\rm true}=1\) | \(A(1)-A(0)\) |
   |---:|---:|---:|
   | 0.5 | 1.680345 | 0.875588 |
   | 0.25 | 1.646720 | 0.856412 |
   | 0.125 | 1.637350 | 0.851929 |

   Halving the **specified** 0.25 step changes absolute \(A\) by **0.569%**, exceeding 0.5%. The runner instead uses 0.5.

   The analytic provider also retains an inadequate default `nk=420`: versus 1600 samples, its derivative changes by 331% in relative array norm over \(10\le r_\perp\le30,\ r_\parallel\le30\). Even 1600 versus 3200 changes that norm by 4.85% at fixed `kmax=35`.

   **Fix:** establish convergence of amplitudes and relevant response contractions. Replace the coincident-table-sample test, which cannot detect interpolation errors.

5. **Significant — Continuum subtraction is not a mean-and-slope projection.**  
   [mock.py:372](/home/anze/Dropbox/work/LyaLenser/code/pipeline/mock.py:372).

   The code subtracts a mean, then regresses against `chi-chi_ref`, whose mean is nonzero. The slope subtraction reintroduces a mean and incompletely removes the slope.

   For the saved radial grid, a synthetic pure linear trend retains **39.7% of its RMS** and **15.8% of its slope**. Saved unlensed skewers have maximum absolute mean approximately 0.199.

   **Fix:** jointly fit the weighted intercept and slope, or center the slope coordinate at the weighted mean distance. Verify both weighted residual moments vanish and projection is idempotent.

6. **Significant — The Wiener ratio is correct, but its spectra do not describe the maps used.**  
   [run_mock_validation.py:99](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:99), [run_mock_validation.py:120](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:120), [templates.py:187](/home/anze/Dropbox/work/LyaLenser/code/pipeline/templates.py:187).

   The runner uses theoretical unmasked \(C_{XX}\), not the required spectrum of the actual masked/noisy map. Its template-density argument remains the theory default 25 deg\(^{-2}\), whereas the saved fiducial mock records approximately 35.96 deg\(^{-2}\). Finite-random noise and actual selection/transfer effects are absent.

   CMB and matched maps also have different pixelizations and masks. The matched mask is returned but not imposed as a common map operator. Magnification changes the catalogues without adding its corresponding terms to the matched \(S_L\).

   **Fix:** apply a common map transfer and mask; use actual density, selection, noise, and magnification in \(S_L\) and \(C_{XX}\). Validate the masked conditional response. Mask independent random-template realizations identically as well.

7. **Significant — Matched-template shot noise is wrong, and the spherical path lacks validation.**  
   [templates.py:212](/home/anze/Dropbox/work/LyaLenser/code/pipeline/templates.py:212), [templates.py:234](/home/anze/Dropbox/work/LyaLenser/code/pipeline/templates.py:234).

   The object weights contain the correct \(W/(b\,\bar n_\chi\,\Omega_{\rm pix})\), random subtraction, and completeness division. But the spherical noise expression effectively uses
   \[
   (1+r)D\langle W^2\rangle/n_{\rm 2D}
   \]
   instead of \((1+r)D^2\langle W^2\rangle/(b^2n_{\rm 2D})\) even in the uniform case. It lacks the general radial-density/bias integral.

   For the fiducial extended slab, the implemented uniform-case result is approximately **0.0121 times the correct noise**.

   Completeness is normalized over *occupied random pixels*, not the footprint, biasing its normalization when random pixels are empty. The flat path returns no noise calculation and does not validate this spherical implementation.

   **Fix:** integrate the specified noise using the same radial selection and bias as the weights; define the footprint independently of random occupancy; test radial selection, completeness, finite-random noise, and flat/spherical agreement.

8. **Significant — The “independent response prediction” uses the wrong modulation field and transfer.**  
   [mock.py:357](/home/anze/Dropbox/work/LyaLenser/code/pipeline/mock.py:357), [amplitude.py:122](/home/anze/Dropbox/work/LyaLenser/code/pipeline/amplitude.py:122), [run_mock_validation.py:143](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:143).

   Forest response depends on the long density at each sampled forest pixel. The prediction instead assigns one density value to each sightline, taken at its selected quasar’s radial position, and multiplies the entire pair’s stored \(\xi G\) sum by those two values.

   Its density template is a smoothed projected field, while the spectral ratio uses the theory’s unsmoothed projected-density spectra. A mode-count-weighted scalar ratio per band is not the actual masked catalogue response. Deprojected prediction is simply set to zero.

   **Fix:** propagate the pixel-level modulation covariance through the continuum operators and pair response, or validate an independently computed equivalent operator. Include the actual smoothing, radial weighting, and map transfer.

9. **Significant — The original 20-seed mean-field gate failed; the amended PASS lacks sufficient protocol provenance.**  
   [run_mock_validation.py:249](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:249), [NOTES.md:180](/home/anze/Dropbox/work/LyaLenser/code/pipeline/NOTES.md:180).

   Combining the final mean/SEM with the saved extra-seed result reconstructs the original matched null:
   \[
   45.3545\pm22.1905,\qquad z=2.04387.
   \]
   That fails the predeclared \(2\,{\rm SEM}\) criterion. Adding seed 20 changes it to \(40.8972\pm21.5728\), or 1.89577 SEM.

   NOTES calls the extra seed predeclared, but the supplied artifacts do not establish prospective declaration before inspecting the failure. The final JSON also compares 21-seed corrected means against 20-seed raw means.

   **Fix:** preserve the original failure and document any prospective stopping rule. Use a fresh fixed validation ensemble; compare raw and corrected results on identical seeds. Separately test fixed-template conditional mean fields, since current templates and geometry vary across seeds.

10. **Significant — Forest, quasar, convergence, and pair geometries are inconsistent across the patch.**  
    [mock.py:311](/home/anze/Dropbox/work/LyaLenser/code/pipeline/mock.py:311), [mock.py:364](/home/anze/Dropbox/work/LyaLenser/code/pipeline/mock.py:364), [mock.py:393](/home/anze/Dropbox/work/LyaLenser/code/pipeline/mock.py:393).

    Forest sampling follows \(x=\chi\theta\). Quasar angular positions use \(x/\chi_{\rm ref}\), while convergence projections integrate fixed transverse grid columns. Those fields therefore do not follow the same radial rays.

    Angular grid conversion also uses fixed \(\cos30^\circ\), whereas pair separations use \(\cos\delta_{\rm mid}\). This error depends on distance across the patch, not merely the small separation of an individual pair.

    **Fix:** use one consistent light-cone/flat-patch coordinate model for all fields, projections, templates, and pair separations. Test local covariance and cross-field responses at several patch positions and distances.

11. **Significant — Periodic mock interpolation extrapolates across the boundary.**  
    [mock.py:284](/home/anze/Dropbox/work/LyaLenser/code/pipeline/mock.py:284), [mock.py:394](/home/anze/Dropbox/work/LyaLenser/code/pipeline/mock.py:394).

    Coordinates are reduced modulo the grid length, but `RegularGridInterpolator(..., fill_value=None)` extrapolates between the final grid sample and coordinates below the grid length. It does not interpolate to sample zero.

    Independent example: for periodic samples `[0,1,2,3]`, interpolation at 3.5 returns **3.5**, instead of **1.5**.

    **Fix:** implement periodic trilinear interpolation or pad periodic endpoint samples. Test boundary-crossing remappings.

12. **Significant — The RSD approximation aliases unresolved density into velocities.**  
    [mock.py:83](/home/anze/Dropbox/work/LyaLenser/code/pipeline/mock.py:83), [mock.py:338](/home/anze/Dropbox/work/LyaLenser/code/pipeline/mock.py:338).

    The Fourier velocity sign is correct. However, taking every fourth cell without an antialiasing filter folds high-frequency density into low-frequency modes before applying \(ik_\parallel/k^2\). Nearest-coarse-cell assignment adds another approximation. “Velocities are large-scale dominated” does not justify aliasing the input density.

    **Fix:** solve velocities from the original Fourier modes, or low-pass before downsampling and interpolate the result. Compare displacement spectra and quasar RSD transfer with the full-resolution solution.

13. **Minor — Quasar counts are Poisson, but not with the specified realization-dependent intensity.**  
    [mock.py:123](/home/anze/Dropbox/work/LyaLenser/code/pipeline/mock.py:123).

    Spatial probabilities are normalized by the realized total lognormal weight, while total object count is drawn from the fixed nominal mean. Thus a common multiplicative intensity change—including mean completeness—does not change expected total counts. The lognormal variance subtraction cancels from the normalized sampling probabilities.

    **Fix:** draw the total from the integrated specified intensity, or explicitly specify and validate a count-conditioned mock. The existing low-variance clustering test does not test this distinction.

14. **Significant — The margin test cannot validate the report’s beyond-Limber claim.**  
    [mock.py:314](/home/anze/Dropbox/work/LyaLenser/code/pipeline/mock.py:314), [run_mock_validation.py:281](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:281).

    Outside-box convergence is independent of box density by construction. The cross-boundary density correlation that the 150 Mpc/h margin is meant to suppress is therefore absent. Changing the margin also changes the generated volume and realization, rather than only changing template selection.

    This is partly a limitation inherited from the mock specification: aggregate convergence-spectrum agreement cannot validate that report claim.

    **Fix:** test changing template margins within one larger correlated density realization, or add an independent non-Limber calculation covering the actual windows and response kernels.

15. **Significant — Sub-slab selection is absent from the compression path.**  
    [pairs.py:107](/home/anze/Dropbox/work/LyaLenser/code/pipeline/pairs.py:107), [xi_model.py:97](/home/anze/Dropbox/work/LyaLenser/code/pipeline/xi_model.py:97).

    Neither kernel consumes slab labels or assigns pixel pairs by midpoint slab. Setting every pixel’s `slab=-1` still produces accepted pairs; I verified this directly. `Config.slabs` does not make the required tomographic catalogue behavior operational.

    **Fix:** implement valid-pixel selection and midpoint assignment to requested sub-slabs in both measurement paths. Test boundary-crossing pairs and disjoint template/sightline selections before the proposed DR1 tomography.

16. **Significant — Several tests verify implementation identities, not the required physical behavior.**  
    [test_templates.py:38](/home/anze/Dropbox/work/LyaLenser/code/pipeline/tests/test_templates.py:38), [test_xi_model.py:19](/home/anze/Dropbox/work/LyaLenser/code/pipeline/tests/test_xi_model.py:19), [run_mock_validation.py:291](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:291).

    The “flat-sky” unit test differentiates an analytic cosine without calling the FFT implementation. The band test checks disjointness but not completeness. Grid convergence compares identical transform samples. There are no end-to-end matched-template, continuum-projection, RSD, or acceptance-statistic tests.

    The shape p-value ignores covariance between bins. From the saved values, correlations reach approximately 0.57 in magnitude. A full-covariance five-contrast calculation gives \(\chi^2=5.088\), nominal \(p=0.405\); accounting for covariance estimation with a Gaussian Hotelling test gives \(p\approx0.565\). It still passes, but **0.783 is not the correctly calibrated result**.

    Injection bookkeeping also fits only unfiltered injection/curl/junk templates, not the seven-component science-band basis.

    **Fix:** add tests targeting these failure mechanisms and use covariance-aware shape statistics.

17. **Minor — Required reviewable outputs are not persisted.**  
    [amplitude.py:13](/home/anze/Dropbox/work/LyaLenser/code/pipeline/amplitude.py:13), [templates.py:133](/home/anze/Dropbox/work/LyaLenser/code/pipeline/templates.py:133), [run_mock_validation.py:322](/home/anze/Dropbox/work/LyaLenser/code/pipeline/run_mock_validation.py:322).

    Per-region \(q,F,mf\) exist only in memory. Constructed flat-sky templates omit the harmonic potential, and JSON omits per-seed mean-field, stochastic, and response predictions. Those omissions prevent several acceptance rows from being independently reconstructed from JSON alone.

    **Fix:** save the complete fit summaries, region partial sums, template operators/potentials, measured-\(\xi\) tables, and per-seed diagnostics with configuration provenance.

The acceptance arithmetic supports the following interpretation:

| Item | Re-derived result | What it establishes |
|---|---|---|
| 1. Tables/benchmark | Spectrum ratios 0.9801, 0.9819, 1.0279 | Useful aggregate map-spectrum and performance check; not numerical convergence or correct joint response. |
| 2. Baseline absorption | Slope ratio 1.00736 | Numerical agreement on the tuned, fixed-baseline calibration sample; absorption untested. |
| 3. Physical normalization | Global slopes 1.00695 and 0.98991 | Limited response check on a patch missing the lowest science band; local slopes and grid checks expose problems. |
| 4. Mean field | Matched 1.896 SEM after extra seed; originally 2.044 SEM | Final predicate passes; original 20-seed gate fails. Residual is poorly bounded. |
| 5. Injection | 1.04724; curl slope 0.0526 versus error 0.7667 | Useful but near-boundary bookkeeping check; no physical calibration. |
| 6. Stochastic normalization | \(-0.5468\pm1.6917\) | Two-SEM interval **[-3.93, 2.84]** admits zero and reversed response. Also uses an increment. |
| 7. Response-only | Matched observed \(42.62\pm8.03\), predicted 459.03; difference \(-416.42\pm410.34\) | Response detection is informative at 5.31 SEM. Prediction agreement is not: a factor-10.8 discrepancy passes. |
| 8. Combined recovery | Increment \(1.3803\pm0.7894\) | Two-SEM interval **[-0.20, 2.96]** includes zero; absolute estimator recovery is not tested. |
| 9. Shape | Published \(p=0.783\); covariance-aware test also passes | No detected inconsistency at this precision, for the differenced statistic. |
| 10. Random templates | Scatter/JK 0.921; scatter/\(\sigma_F\) 2.770 | Useful conditional error-scale diagnostic; explicitly no acceptance tolerance. |
| 11. Covariance | 0.86517 | Correct arithmetic for paired increments, not the requested absolute estimator. |

In particular, **410.34 is the SEM of the observed-minus-predicted response difference**, not a separately reported prediction uncertainty. The response row passes because this difference is noisy. Similarly, the deprojected response-only result \(0.726\pm1.910\) permits substantial residual contamination.

For the requested saved-mock truth check, seed 0 with response emulation off and the **unfiltered truth template**, jointly fitted with curl and junk, gives
\[
\hat A=1.0351,\qquad \sigma_{\rm JK}=0.9883,\qquad \sigma_F=0.3589.
\]
That is consistent with \(A=1\), but its precision cannot establish percent-level normalization.

**Verdict: require a third implementation round.** It must:

- Restore the specified absolute acceptance statistics; retain paired differences separately.
- Complete the tapered response basis and enforce valid science/curl/junk fits.
- Correct continuum projection, periodic interpolation, coordinate consistency, and RSD sampling.
- Fix matched-template noise/completeness and validate the actual masked Wiener operator.
- Implement and test sub-slab selection and persist all required fit diagnostics.
- Establish amplitude-level numerical convergence and independently validate pixel-level response predictions.
- Freeze numerical choices, tolerances, seed counts, and stopping rules before running a fresh validation set. Demonstrate meaningful normalization and residual-bias bounds; merely failing to reject unity or zero with enormous errors is insufficient.
- Validate margin and selection effects with controlled correlated realizations.

The pair-compression foundation is usable. The present **11/11 verdict is not a defensible authorization to start Stage B**.
