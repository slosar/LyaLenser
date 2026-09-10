# Codex adversarial review 3 (gpt-6-astra, high), confirmation pass on the revised estimator spec — 2026-09-10

Scope: report Sec. 5 and 4.2 (revised), PLAN.md Sec. 0, IMPLEMENTATION.md v2, mathematica/pair_template_checks.wls checks 9-10. Verbatim output follows.

---

**IMPLEMENTATION.md is not ready for handoff as an exact specification.** The revisions resolve several central objections, but remaining sign, normalization, mock-construction, and validation gaps can still produce biased amplitudes.

No files were modified. I read the requested revisions and ran the Mathematica checks: check 9 returns -kappa r xi xi', and check 10 returns zero, as intended. I also independently checked the sign, first-moment normalization, and variance counterexamples below.

| Review-2 item | Status | Confirmation |
|---|---|---|
| 1. Amplitude definition/filter | Partially resolves | h_L=S_L/C_L^{XX} is correct under the stated conditional model. Masking, estimated spectra, and data-dependent kernels remain unspecified. |
| 2. Band response | Partially resolves | Full F_bc is specified; omitted-band response and tapered-band windows are not. |
| 3. Mean-field oddness | Resolves | Report Sec. 5.2 and check 9 correctly withdraw cancellation by oddness. The acceptance test still needs correction. |
| 4. Injection calibration | Partially resolves | Coordinate shifts are correctly demoted to bookkeeping tests. Claims about isolating boundary effects and the required injection slope remain too strong. |
| 5. Projected covariance response | Partially resolves | Geometry averaging is acknowledged as an approximation; the proposed mocks do not adequately validate the production approximation. |
| 6. Compression | Resolves | The numerator identity is correct, and six explicit separation bins support the newly specified real-space diagnostic. |
| 7. Covariance | Partially resolves | Shared-pair correlations are acknowledged. The replacement "lower bound" claim is also false; covariance validation remains insufficient. |
| 8. Curl leakage | Resolves | Joint gradient-curl fitting correctly handles leakage within the fitted response model. |
| 9. Simulation ensembles | Partially resolves | Random-template and joint ensembles are distinguished. Conditional mean-field testing and covariance precision remain inadequate. |
| 10. Radial matching | Partially resolves | Relevant residuals are named, but the mock cannot validate the claimed margin or all density/tidal/RSD transfers. |
| 11. Map normalization/noise | Partially resolves | Angular-density and pixel-area factors are correct. Completeness, common transfer, and finite-random noise are missing. |
| 12. Magnification | Partially resolves | Intrinsic-response leakage is now acknowledged; the mock omits the correlation needed to measure it. |
| 13. Stochasticity | Resolves | The inappropriate universal 1-r^2 residual is replaced by the relevant cross-correlation assumptions. |
| 14. Shared sampling | Partially resolves | Early tomography is adopted, but executable disjoint-selection rules and representative mocks are absent. |
| 15. Higher-order response | Partially resolves | The second-order interpretation is corrected; no quantitative higher-order acceptance test is required. |
| 16. Costs/storage | Partially resolves | Operation costs and precision are specified. The report's approximately 1 GB estimate is low: the 36 float64 accumulators alone occupy 2.88 GB for 1e7 pairs. |

The remaining findings, with fixes:

1. **Blocker — The executable pair direction reverses the response sign.** IMPLEMENTATION.md Sec. 0, 1.4, 2.7 step 5. The specified `dx = ra_b-ra_a`, `dy = dec_b-dec_a` points from a to b, while the text says b to a and the contraction uses alpha_a - alpha_b. Direct check: for xi(r)=e^{-r^2/2}, r=1, alpha(x)=-kappa x, remapping gives dxi/dkappa=+0.60653; the specified contraction gives -0.60653. Fix before implementation: define separation as theta_a - theta_b or reverse the deflection difference consistently, and add an independent finite-difference sign test. Sampling at theta_obs+alpha in Sec. 2.7 is itself correct.

2. **Blocker — First-moment normalization is wrong; midpoint storage omits a first-order response.** IMPLEMENTATION.md Sec. 1.2-1.4, 2.5; report Sec. 5.2. m_eff should be m+2g1 m1+g1^2 m2, not the linear prescription given (which returns A=1.2 for a unit signal with constant g=1.2). Also g_p alpha_a - g_q alpha_b has a first-order cross term in (chi_p - chi_q) that is not stored.

3. **Significant — Wiener filter valid only for a more restricted operator than implemented.** Sec. 2.4-2.5, 3.3. h_L=S_L/C_L^{XX} is correct for homogeneous fixed spectra, but the actual masked/filtered map needs the full C_{kx} C_{xx}^{-1} operator, and `three_tracer.spectra()` hard-codes the original source plane/slab rather than the actual constructed map's spectrum.

4. **Significant — Baseline fitting subtracts response, not just adds noise.** Sec. 2.2, 2.8 steps 1-2. Same pair products determine the fitted baseline and the score, so dE[Q]/dA = R^T W (I-H) R, not R^T W R: response can be partially absorbed by the fit.

5. **Significant — Random subtraction alone does not produce a kernel-matched map.** Sec. 2.4 `matched_template`. Completeness c(theta,chi) leaves a multiplicative transfer uncorrected; finite randoms add noise (~1/20 of data Poisson term at 20x density).

6. **Significant — Joint mock lacks the correlations it's meant to validate.** Sec. 2.7 steps 1-4, 2.8 step 5. Independently generated kappa_Lya and kappa_rest force zero forest-lensing/slab cross-spectrum, removing the intrinsic magnification leakage the mock is supposed to test; smoothed delta_L (quasars) vs. unsmoothed delta_m (kappa_slab) kernels also mismatch.

7. **Significant — Acceptance tests (Sec. 2.8) can pass wrong normalization, wrong mean field, or failed deprojection.** Needs distinct tests: fixed-template mean-field subtraction, fixed-geometry physical remapping at several amplitudes, stochastic normalization, response-only null test, and combined recovery with independently predicted response-control expectation, all with predeclared tolerances.

8. **Significant — Band matrix doesn't cover omitted modes or the claimed Fourier shape test.** Sec. 2.4-2.5, 3.4. E[A]=A+F^{-1} f_out with omitted modes; six real-space bins don't support the advertised D/T vs. k diagnostic.

9. **Significant — Uncertainty prescription (1/F) is not a general lower bound**, shown by explicit counterexample (F=2, Var=1/4 < 1/F=1/2); jackknife pair assignment to sightline a only doesn't perform true spatial deletion.

10. **Significant — Tomographic selection/systematic gates contradict the adopted design.** Sec. 3.4, PLAN.md Sec. 0: sub-slab restriction conflicts with extended template slab; z<3.0 read vs. z=3.2 sub-slab; L-range and bias treatment inconsistent with promised inference procedure.

11. **Significant — Several numerical conventions need independent tests**: k_par integral normalization (1/pi vs 1/(2pi) depending on integration limits), `alm2map_der1` already divides by sin theta (risk of double division), angular-to-Cartesian mapping and grid convergence for the mock.

**Verdict:** the report is materially improved, but IMPLEMENTATION.md should not yet be handed over as an exact, acceptance-ready specification. Correct the pair sign, response moments, conditional/map operators, baseline-fitting response, mock covariance, and executable validation/selection rules first, with explicit acceptance gates before real-data inference.
