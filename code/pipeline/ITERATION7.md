# Stage A, iteration 7 brief (2026-09-14, draft for the user's decision; nothing implemented yet)

## Where we are
- Iteration 6 (fresh N = 400 + 20 seeds, GATES v6 frozen, Perlmutter): 28/36 required rows pass. Every fix asked
  for by review 5 holds. The deprojected estimator has an additive bias of about -0.4 A (2-3 SEM) that is present
  with the response switched OFF and only with the realisation's own quasar template; the linear response is
  deprojected correctly; the deprojected lensing normalisation is 0.88-0.89 +- 0.06.
- Mechanism (NOTES.md "Scale-1 campaign" and "Long-mode test result"; review 6 agrees it is established as a
  contributor, not yet closed quantitatively): the pair product delta_p delta_q contains the product of the
  forest's long-wavelength components, m_p m_q, and its correlation with a non-Gaussian quasar template,
  <m_p m_q T> = b_q^2 xi_ms(p) xi_ms(q) for the mock's lognormal tracer, is not removed by the one-coefficient
  deprojection (which removes the squeezed, linear-response term). With the modulation on, the Gaussian
  contractions <delta_L m_p><m_q T> are the 15-21 % excess of the response over P(DCD-C)P^T (review 6 gives the
  full ensemble-mean formulae and the exact conditional form: modulation acts on M(l) = Sigma + mu mu^T, not on
  the unconditional C). In the real universe the analogue is the quasar b_2 (4.6-5.0 for b_q = 3.5 against the
  lognormal's 12.25) plus tidal bias, F2 and the forest's own non-linearity; the mock's -0.4 A does not rescale
  to a DR1 number.
- The 96-seed long-mode test: the m-only field scores -1.12 +- 0.29 A deprojected (mechanism), and subtracting
  c P[delta_L] from the forest halves the raw lensing response but leaves the calibrated precision (scatter /
  response 3.72 -> 3.63) and the calibrated bias (-0.92 -> -0.90) unchanged: that particular subtraction is
  not a fix, but review 6 is right that it does not rule out conditional methods on S/N grounds (my earlier
  "45 % of the signal is lost" was the raw response, not S/N).
- Review 6 verdict: Stage A not passed; Stage B blocked; no trustworthy DR1 upper limit until the estimator is
  calibrated or the contaminant has an explicit likelihood; "Stage A could pass after one more substantive
  development-and-validation iteration".

## Decision needed: which estimator change (review 6, question 3)
| Option | What it is | Reviewer's verdict | Cost / risk |
|---|---|---|---|
| (e) conditional pair-moment subtraction (preferred by review 6) | Subtract M_0(Q) = E[d d^T | quasar field, A = 0] from the pair products, including its covariance part, with the lensing response of the subtracted statistic made explicit; staged: oracle slab density -> continuous intensity -> sampled quasars -> data-usable construction | Preferred development path; retains information that field subtraction throws away | Needs a model of the forest-quasar conditional moment (b_F, b_q, b_2, shot noise, RSD); its uncertainty becomes nuisance parameters; cost measured, not known |
| (a) second nuisance template with the <m m T> shape | Add a template quadratic in the quasar density (pair-covariance shape over separation, L, redshift) and fit its amplitude jointly with lensing | Promising as the first component of (e); a squared 2D quasar map is NOT automatically the right shape; needs tracer bias / cross-covariance, radial selection, shot-noise treatment (split catalogues) | S/N retention sqrt(1 - rho^2) with rho the nuisance-lensing overlap in the covariance metric; not computable from present products |
| (d) response prediction conditioned on the long field | Predict with M(l) instead of C | Required for a credible prediction; a calibration improvement, not a fix by itself | None on S/N; leaves the additive bias in place |
| (b) subtract a predicted forest mean | The 96-seed test with a scalar c P[delta_L] fails its null; a 3-D conditional mean with recomputed operators is not ruled out | Not ruled out in principle; the scalar oracle version is not it | Data version has quasar shot noise and reconstruction error |
| (c) cuts in r_par or k_par | Diagnostic only | Not a fix: the contaminant sits in the bins with the information | r_par < 5 keeps ~45 % of the S/N and keeps the contaminant |

My recommendation: (e) with (a) as its first concrete component and (d) as the prediction, developed on the
present mock (kept as the hard regression test) plus a second validation tier with independently variable
quasar b_2 / tidal bias and a non-linear redshift-space forest, and with GATES v7 written before the campaign.

## What GATES v7 must contain (review 6, question 6)
Frozen estimator + nuisance model + transfer estimator + cuts + development-derived calibration, then genuinely
fresh acceptance seeds; separate absolute own-template null, paired response-minus-complete-prediction, paired
lensing normalisation for truth / CMB / deprojected, absolute combined recovery; equivalence bounds with
uncertainty rather than 2-SEM consistency (the matched slope "passes" at -0.77 +- 1.47); mechanism controls that
vary tracer non-linearity independently of lensing and response (oracle density, continuous intensity, sampled
template, selection); shape tests against the full predicted lensing-plus-nuisance vector with joint covariance;
known-covariance tests of the production selection, heterogeneous projection, radial support, noisy-map
transfer; a power calculation for the whole correlated gate set (N = 400 gives only ~54 % power for a 0.5 A
bound at zero bias with per-seed scatter 3.7; 576 / 712 seeds for 80 / 90 %, before nuisance costs); immutable
production provenance with per-seed outputs.

## Secondary items from review 6 to fold into iteration 7
- CMB-template vs truth-template normalisation difference is 1.57 SEM (-0.072 +- 0.046): not established;
  suspects listed (estimated Wiener transfer S / C^XX from the noisy map; fitted kernel on the same noisy
  sample; projection / binning approximations). Test with fixed geometry and known covariance.
- `condor/reprovenance.py`: 1e-9 relative float tolerance accepted without a manifest entry; basis not in the
  migration inventory. Do not reuse as is.
- `response._predict` has no selected-slab argument (`cfg.slab_index`): fine for the full slab, wrong for
  tomography.
- Generator still uses the nominal 20 x scale side for the outside-box fields, CMB-noise pixel area, mask edge
  distance and random-catalogue support (`mock.py` 413-419, 470, 543).
- Test count is 64, not 57 (HANDOFF).
- The 96-seed subset is not independent confirmation; report the m-only six-bin vector next time.

## Questions for the user
1. Go with (e)+(a)+(d) as above, or a different option?
2. Second mock tier (realistic quasar b_2, non-linear forest): build it in iteration 7 or after the estimator
   change is validated on the present mock?
3. Ensemble size for the next campaign given the power calculation (~600 seeds ~ 60 node-h on Perlmutter).
4. Who implements: Claude (with astra review after) or astra (Codex now works on RACF in workspace-write mode)?
