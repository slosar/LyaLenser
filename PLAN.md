# PLAN v2: measuring Lyman-alpha forest lensing x CMB lensing in DESI DR1

Written 2026-09-10; revised the same day after codex review 1 (report/reviews/codex_review_1.md) and the
switch from a cell-based to a pair-template estimator (report Sec. 5). Status: awaiting codex review 2 of the
estimator; then the implementation plan (IMPLEMENTATION.md) is written and handed to codex.

## 0. Design decisions (revised after codex review 2)

1. **Pair-template estimator, no cells, no kappa map (report Sec. 5).** For every pair of pixels on distinct
   sightlines the response to a common deflection field is exactly grad xi . [alpha(theta_a) - alpha(theta_b)] at
   first order, for any L. The amplitude A of C_L^{kappa_lya X} = A S_L is estimated as a *conditional template
   coefficient*: template alpha_T = grad phi_T with phi_T(L) = 2 h_L X(L)/L^2, h_L = S_L / C_L^{XX} (Wiener filter of
   the map X for kappa_lya), Ahat = [sum w w d d R - b]/F. Band amplitudes use the full band response matrix F_bc.
   The mean field b is even under pair reversal and is computed exactly. 1/F is a lower bound on the variance.
2. **Compression.** Per sightline pair store scalars v, m, beta (plus first (chi_mid - chi_ref) moments for the
   redshift dependence of the deflection, and a few (r_perp, r_par) bins for the shape test). Numerator compresses
   further to one vector per sightline; F needs the pair level. Any change of xi model, weights, masks or bins needs
   a new pixel pass.
3. **Normalisation** is model dependent through grad xi, measured from the same data (distorted, geometry-averaged
   correlation function at r < 40 Mpc/h). Injections by shifting sightline coordinates test the geometric bookkeeping
   (pair re-association, signs, bands, mean field) but do NOT calibrate the physical response (review-2 item 4); the
   physical normalisation is validated on joint mocks that apply the same measured-xi procedure.
4. **Deprojection with a kernel-matched quasar template** kappa_s-hat = sum_i W_CMB(chi_i)/(b_q(z_i) nbar_3D chi_i^2)
   minus the random-catalogue mean, over a template slab that extends ~150 Mpc/h beyond the forest slab on each side
   (non-Limber inside/outside correlations); b_q(z) imposed from the quasar auto-correlation; magnification m(z) and
   quasar RSD modelled; the clean multipole range starts where the RSD correction is a few per cent (L >~ 100).
5. **Templates and controls:** signal X = kappa_CMB - kappa_s-hat (ACT DR6, Planck PR4); response control
   X = kappa_s-hat; curl template fitted jointly with its geometric leakage F_RT; random CMB-simulation templates
   as a conditional diagnostic; joint mocks for the deprojected covariance and the mean field of the real template.
6. **Tomographic sub-slabs with disjoint template/sightline selections are part of the first measurement**, not of a
   later phase; the single 2.1-3.0 slab is used for forecasts and for the first end-to-end run only.
7. **Validation:** joint mocks (sparse sightlines with the real geometry, continuum projection, DR1 noise
   distribution, forest response to long modes correlated with mock quasars and with the slab part of kappa_CMB,
   quasar RSD and magnification) generated with fixed observed geometry; sightline-shift lensing is used for
   pipeline tests only. Paired +-A injections with zero-injection subtraction.

## 1. Data (all public)

| Data | Source | Size (est.) |
|---|---|---|
| DESI DR1 Lya deltas (picca, per-HEALPix `Delta-*.fits.gz`, Lya/Lyb/CIII regions) | https://data.desi.lbl.gov/public/dr1/vac/dr1/lya-deltas/ | 10-40 GB (verify) |
| DESI DR1 QSO catalogue + LSS randoms | DESI DR1 LSS catalogues | few GB |
| DESI DR1 DLA catalogue | https://data.desi.lbl.gov/doc/releases/dr1/vac/dla-cnn-gp/ | small |
| ACT DR6 lensing: kappa alm, masks, N_L, 400 sims (baseline and tSZ-deprojected); follow the LAMBDA mask/transfer instructions | LAMBDA | ~30 GB with sims |
| Planck PR4 lensing: klm, mask, N_L, sims | PLA / LAMBDA (Carron+2022) | ~10-20 GB with sims |
| (phase 2) low-z galaxy sample for 5s-2 (DESI Legacy LRGs) | public | few GB |

Store under `/data/LyaLenser/`; record paths in MEMORY.md.

## 2. Pipeline modules (`code/pipeline/`)

| # | Module | What it does | Cost here |
|---|---|---|---|
| 0 | `fetch.py` | download + checksum | I/O |
| 1 | `prep.py` | read deltas -> (RA, Dec, z_q, chi, delta, w); DLA/metal masks; sub-slab assignment; KD-tree of sightlines; per-pixel noise distribution (for the forecast update) | minutes |
| 2 | `xi_model.py` | mean xi(r_perp, r_par) and grad xi: from the DESI DR1 fitted model with distortion, cross-checked against a direct pair estimate; tabulated on a fine grid for the kernel | minutes |
| 3 | `pairs.py` (numba, parallel) | one pass over pixel pairs on distinct sightlines (r_perp<30, abs(r_par)<30 Mpc/h): accumulate V_ab, M_ab, B_ab per sightline pair and per sub-slab; option to run with shifted sightline positions (injection) | to be benchmarked; arithmetic bound minutes, realistic 0.5-3 h |
| 4 | `templates.py` | kappa_CMB alm -> filtered kappa_T -> alpha_T at quasar positions (healpy alm2map_der1 or flat-sky FFT per patch); kernel-matched kappa_s-hat from slab quasars and randoms; rotated (curl) templates; simulation templates | minutes per template set |
| 5 | `amplitude.py` | Ahat, F, b for each template and L band from the pair catalogue; jackknife over sky regions; null distributions from sim templates | seconds per template |
| 6 | `inject.py` | generate alpha_inj realisations, shift sightlines by -+A alpha_inj, rerun module 3, recover A; tests bookkeeping, signs, band matrix and mean field (not the physical normalisation) | one pair pass per injection (5-10 injections) |
| 7 | `mocks.py` | joint mocks (forest with long-mode response, biased quasars with RSD and magnification, kappa_lya, kappa_CMB with its slab part from the box), fixed observed sightline geometry, DR1 noise distribution, continuum projection; run modules 2-5 with the measured-xi procedure; recovery of A, deprojection residual, mean field, covariance | one pair pass per mock (>= 20 mocks) |
| 8 | `nulls.py` | curl template, sim templates, splits (magnitude, redshift, NGC/SGC), sightline-density x kappa mean-field check, response amplitude vs Karacayli+2024 | minutes |

Dependencies: numba (present), healpy, astropy/fitsio, h5py, pymaster (only for auxiliary spectra: b_q from C^{qq}, q x g).

## 3. Resource assessment

- **This machine (24 threads, 62 GB, ~880 GB free, RTX 3050 6 GB) is sufficient for phases 1-2.** The pair pass is
  the only heavy step; the review is right that 1.4e11 pixel pairs at 1e8/s/thread is ~2 min, so the real cost is set
  by the kernel (2D interpolation of grad xi, scattered access); module 3 is benchmarked on one HEALPix pixel first
  and the count extrapolated. Injections and mocks multiply the pair pass by 10-30, which fits in days of wall time.
- **24 GB GPU: not required.** Would make the pair pass and a phase-3 CG solve for full C^-1 weighting cheap.
- **NERSC: not required for the measurement.** Needed for the official DESI Lya mocks and for regenerating deltas.

## 4. Phases

1. **Phase 1 (estimator validation on mocks, then first measurement):** modules 2-7 on joint mocks until the
   acceptance tests of IMPLEMENTATION.md pass; then modules 0-6 on DR1 x ACT DR6 and Planck PR4, first single slab
   end-to-end, then tomographic sub-slabs with disjoint selections; injections; jackknife; random-template
   diagnostic; joint-mock covariance. Deliverable: Ahat for signal, response and curl templates with mock-validated
   errors and a bound on the residual response.
2. **Phase 2 (systematics):** magnification slope from q x LRG; DLA/metal/continuum tests; tSZ-deprojected ACT map;
   released N_L and mask transfer functions; b_q(z) from C^{qq}; non-Limber RSD/tidal transfer tests.
3. **Phase 3:** full C^-1 weighting (CG), DESI mocks at NERSC, paper.

## 5. Open items from review 1 carried forward

- Higher-order lensing corrections to the cross-response (Gaussian third-order terms; lensing bispectrum).
- Joint likelihood forecast with b_q(z), m and the response amplitude free.
- Validation of the sparse-sampling noise model (independent modes with P_1D/nbar; l_max cutoff) against the discrete
  estimator on mocks.
- Sightline-density / quasar-template correlation (shared quasars) modelled in the joint mocks.
- Radial transfer mismatch between density, tidal and quasar-RSD responses within a sub-slab.

## 6. Codex review 2 (estimator only): done 2026-09-10, report/reviews/codex_review_2.md; all 16 items folded into Sec. 0 above, report Sec. 5 and IMPLEMENTATION.md. Original checklist kept for reference:

Review report Sec. 5 and this file's Sec. 0: (1) the exactness claim for Eq. (Rij) and the first-order Gaussian
estimator (Ahat), its normalisation F and mean field b with diagonal weights; (2) the compression to per-pair
(V, M, B) and the statement that no per-sightline compression exists; (3) the Wiener-filter template choice as the
Fisher-optimal amplitude estimator in the diagonal-weight approximation; (4) the injection-calibration argument
(geometric response, continuum projection along the LOS); (5) the curl-template null and the conditional nature of the
sim-template null ensemble; (6) the kernel-matched template (Eq. matched) incl. shot noise and magnification; (7) any
bias that survives: same-sightline terms, sightline-density correlations, second-order deflection terms.
