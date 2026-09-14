# Review request 6 (2026-09-14): adversarial review of Stage A round 6 and of the deprojection failure

Model: gpt-6-astra, effort high, read-only. Output saved verbatim as `report/reviews/codex_review_6.md`.

## Context to read first
`HANDOFF.md`; `PROGRESS.md` (2026-09-14 entries); `code/pipeline/NOTES.md` section "Iteration 6" and, in it,
"Scale-1 campaign" (the reading of the campaign and the diagnosis under review) and "Long-mode test";
`GATES.md` (v6, frozen before the run); `code/pipeline/ITERATION6.md`; `report/reviews/codex_review_5.md` (your
previous review; every blocker was implemented); `report/mock_validation.md` and `.json` (the iteration-6
acceptance table: 54 rows, 36 required, 28 pass; the JSON holds every per-seed fit under `sparse[*].fits`,
`predictions`, `fixed_null`, `shared`, `shape`); `report/longmode_iteration6_summary.json` (the long-mode test);
`report/main.pdf` Sec. 4 (three-tracer separation, squeezed limit) and Sec. 5 (pair-template estimator);
`IMPLEMENTATION.md` Sec. 0-1.

## Code under review (changes since review 5, commit 7416f4b onwards)
- `code/pipeline/template_audit.py` (single mask), `response.py` (production selection in `_predict`),
  `mock.py` (`patch_side_rad`, `zq`), `run_mock_validation.py` (`make_bundles` with the generator patch side,
  `process_seed` roles full/core, A in {0, 1} with response on/off for every seed, `build_rows` with `required`,
  the paired deprojected row, jackknife nside bookkeeping), `pairs.py` (`accumulate(true_positions=...)`
  expectation mode), `inject.py` (`injection_test(expectation=True)`), `campaign4.py` (iteration 6 seed ranges,
  `campaign_config` in the fingerprint, basis marker pinned, `DEPROJECTED_BOUND`), `condor/reprovenance.py`
  (manifest), `tests/test_iteration6.py`.
- `code/pipeline/longmode_diagnosis.py` (written after the campaign; reads the products, bypasses provenance).

## The result
The campaign ran on fresh seeds (400 sparse + 20 dense, scale 1) with the protocol frozen in GATES v6. Every row
that review 5 asked to fix passes (single-mask sampled template 1.008 +- 0.006; dense normalisation 1.001 +-
0.008 g on and off; injection expectation 0.979; fitted-vs-generator; first moments; spectra; covariance). The
deprojected rows fail, and the ensemble shows why: with the response OFF and the realisation's own templates the
deprojected null is already -0.41 +- 0.18 A (matched-template null +10.3 +- 1.6, CMB +0.33 +- 0.14; the
fixed-template and random-template rows are zero); the linear response is deprojected correctly (paired
A0_R1 - A0_R0 = +0.02 +- 0.08 vs 0.10 predicted) although the response itself exceeds P(DCD-C)P^T by 15-20 %
on both templates; the lensing response of the deprojected estimator is 0.88-0.89 +- 0.06. NOTES attributes the
response-off term to the non-squeezed part of the forest-forest-template bispectrum, <m_p m_q T> with m the
forest's long-wavelength component (exact b_q^2 xi_ms(p) xi_ms(q) for the lognormal template), and the response
excess to the Gaussian four-point terms <delta_L m_p><m_q T>; the six-bin shape failure (p = 0.002) has the
matching signature. The long-mode test (96 seeds) subtracts c P[delta_L] per seed with the kernel fixed: the
m-only field (m_p m_q against its own xi; no lensing, no response, the seed's own templates) scores matched
+11.5 +- 2.3 and deprojected -1.12 +- 0.29 A, which establishes <m m T> != 0 with the sign and size of the bias;
but the subtraction also removes 45 % of the lensing signal and 70 % of the response, because c P[delta_L]
(10 Mpc/h smoothed, 0.9 % of the pixel variance) carries 19-75 % of xi at r_perp = 3-20 Mpc/h: at the kernel's
separations nothing is squeezed relative to the template's modes.

## Questions (answer each with a verdict and the evidence in the code or the products)
1. Is the diagnosis right? Verify from the per-seed JSON that (a) the deprojected A0_R0 mean is -0.41 +- 0.18,
   (b) the paired response and lensing differences are as stated, (c) the regression of the deprojected amplitude
   on the cmb and matched amplitudes gives (1.04, -0.066) and reproduces the response-off null from the two
   mean-field rows, (d) the long-mode test removes the own-template null and the response excess to the extent
   NOTES claims. Is there an alternative explanation you can support (a bookkeeping error in `fit_save`'s mean
   field, a template-construction correlation, the completeness pattern, the ACT mask, the lognormal Poisson
   sampling, the disjoint/shared selection) that the fixed-template rows would NOT have exposed?
2. Derive, or check against a derivation, the ensemble mean of the estimator's score for (i) a Gaussian forest
   with a lognormal template, response off, (ii) the same with the modulation (1 + R delta_L / 2) on, and (iii) the
   real universe at tree level (forest-forest-quasar bispectrum with b_1, b_2, F2 and the forest's own
   nonlinearity). Which terms does the report's deprojection (one kernel-matched template, one coefficient) remove
   and which does it not? Is the claim that the pairs (3-30 Mpc/h) against template modes (30-250 Mpc/h) are not
   in the squeezed regime correct, and how large is the non-squeezed part for DESI DR1 quasars (b_q ~ 3.5,
   b_2 from the peak-background split) relative to the mock's lognormal (b_2 = b_q^2)?
3. What is the right fix for the estimator? Candidates: (a) a second nuisance template with the shape of the
   <m m T> term (quadratic in the quasar density), fitted jointly; (b) subtracting the conditional mean of each
   forest pixel given the quasar field (the forest's long modes predicted from the template) before forming
   pair products, and its residual (1 - r^2) — note the test above: removing the long modes removes half of the
   lensing signal, so judge whether any version of (b) can work; (c) restricting the pair kernel in r_par or k_par; (d) a modified
   response prediction that conditions on the full long-mode field; (e) something else. For each: what does it
   cost in S/N (report Table 2 numbers), what does it require from the data, and how would Stage A test it?
   Also judge whether the mock is even the right tool (linear Gaussian forest, lognormal template with b_2 =
   12): should the next iteration use a mock with a realistic quasar b_2 and a non-linear forest?
4. The CMB-template lensing normalisation is 0.887 +- 0.048 (paired, 400 seeds) against 0.959 +- 0.019 for the
   truth template on the same seeds and 1.001 +- 0.008 for the noiseless dense mocks, while the spectra row is
   within 2 %. Fluctuation or a real 4 % (truth) and 11 % (CMB) deficit of the noisy sparse sample, and where in
   the code would it come from (Wiener transfer S/C^XX measured on the noisy map; noise x kernel; the r_perp >=
   3 cut; the fitted table on the noisy sample)?
5. Round-6 code: is every review-5 fix implemented correctly (single mask, production selection in the
   prediction, patch side, zq, provenance, jackknife, expectation mode, required flag, paired row)? Is the
   noise-free injection expectation a valid replacement for the measured injection gate?
6. Protocol: with the deprojected bias now resolved at 2-3 SEM by the fixed ensemble, is anything in the way the
   result was read (all 400 seeds for null and recovery, correlated rows declared) a tuning or a multiple-look
   problem? What should GATES v7 require so that the next campaign is decisive?
7. Stage B: given this, what can be measured on DR1 x ACT/Planck now (response term, upper limit) and what must
   wait for the estimator fix? Anything else that should block Stage B?

Format: numbered findings with severity (blocker / significant / minor), file and line references, and a final
verdict: is the diagnosis established, what is the fix, is Stage A passable in one more iteration. Do not modify
any file.
