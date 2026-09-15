# CMB-lensing cross-correlation: future work (parked 2026-09-15)

State when parked: Stage A iteration 6 (`report/mock_validation.md`, `code/pipeline/NOTES.md` "Iteration 6",
`report/reviews/codex_review_6.md`). The pair-template estimator cross-correlated with kappa_CMB and deprojected
with the kernel-matched quasar template is unbiased for the linear (squeezed) density response but carries an
additive bias from the second-order forest-forest-tracer term <m_p m_q T> (m = the forest's long-wavelength
component, T = the realisation's non-Gaussian quasar template): -0.41 +- 0.18 A on the mock with a lognormal
tracer (b_2 = b_q^2 = 12.25), with no lensing and no response. The linear response and the lensing signal act
through the same long modes at the kernel's separations, so the term cannot be filtered out; it has to be
modelled and marginalised. With the modelling available here, the CMB cross-correlation measures a combination
of lensing and the intrinsic forest-forest-density three-point function that a single template cannot separate.
The project therefore pivoted to low-redshift tracers (ITERATION7.md), where the term vanishes because the
template shares no density modes with the forest.

## Direction to pursue later: a constraint in the (A_lens, A_2) plane
Instead of claiming a lensing measurement from kappa_CMB alone, fit two amplitudes to the same pair products:

    E[delta_p delta_q | maps] = xi_pq + A_lens * (lensing dilation term, template = Wiener kappa_CMB)
                              + (kernel-matched response term, deprojected as now)
                              + A_2 * Q_pq

with Q_pq the quadratic template, the product of two forest-quasar cross-correlations integrated over the
template map, Q_pq = sum_x K(x) xi_Fq(p, x) xi_Fq(q, x). Its shape is fixed by xi_Fq(r), which DESI DR1 measures
at high significance (the Lya-QSO cross-correlation of the BAO analysis), so A_2 has no free shape. A_2 has a
physical reading: the effective second-order bias of the quasar field relative to the tree-level expectation
(b_2 ~ 4.6-5.0 at b_1 = 3.5 from the halo fit / peak-background split; F_2 and the forest's own non-linearity
enter the same shape at leading order; tidal bias does not and would need a second shape).

What the present products already say about it:
- The six-bin shape test (p = 0.002) shows the contaminant changes sign at r_par > 10 Mpc/h where the lensing
  signal does not: the two shapes are distinguishable, the degeneracy will be partial.
- The matched-template null was 6.5 SEM on a 400 deg^2 mock; DR1 has 25 times the area, so A_2 will be a
  detection of a new three-point statistic (complementary to the P1D x kappa bispectrum of Karacayli et al.),
  even with a real b_2 of ~5 rather than the mock's 12. A_lens stays an upper limit at DR1 depth (S/N ~ 1).
- The mock validates the joint fit cleanly because its A_2 is known exactly (b_q^2 for the lognormal tracer);
  it validates only the lognormal shape. Real-universe terms the single shape does not absorb would bias A_lens.

Steps when resumed:
1. Build Q_pq in `pairs.accumulate` (a second kernel accumulator, xi_Fq from the fitted forest-quasar table or
   from the generator in the mock), a template Q with its own F entries, joint (A_lens, A_2, curl, junk) fit.
2. Validate on the iteration-6 products (seeds 1000-1399 are frozen; a new campaign with GATES v7-CMB): A_2
   recovered at b_q^2, A_lens unbiased, degeneracy coefficient rho reported, S/N retention sqrt(1 - rho^2).
3. Reviewer's list to fold in (review 6): conditional response prediction with M(l) = Sigma + mu mu^T; a
   second mock tier with independent quasar b_2 / tidal bias and a non-linear forest; equivalence bounds and a
   power calculation in the gates (N = 400 gives ~54 % power for a 0.5 A bound); reprovenance float tolerance
   and basis inventory; `response._predict` without `slab_index`; nominal side in four generator places.
4. Stage B for the CMB thread would then report the (A_lens, A_2) ellipse for ACT and Planck, A_2 as the
   measurement, A_lens as the limit, alongside the low-redshift result.

## Alternative kept on record
Review 6's preferred development path: subtract the conditional pair second moment M_0(Q) = E[d d^T | quasar
field, A = 0] including its covariance part, staged oracle density -> continuous intensity -> sampled quasars ->
data-usable construction; potentially retains information that field subtraction throws away; cost to be
measured. The scalar long-mode subtraction test (`longmode_diagnosis.py`, 96 seeds) fails its null and
changes neither the calibrated precision nor the calibrated bias.
