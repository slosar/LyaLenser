# PLAN v2: measuring Lyman-alpha forest lensing x CMB lensing in DESI DR1

Written 2026-09-10; revised the same day after codex review 1 (report/reviews/codex_review_1.md) and the
switch from a cell-based to a pair-template estimator (report Sec. 5). Status: awaiting codex review 2 of the
estimator; then the implementation plan (IMPLEMENTATION.md) is written and handed to codex.

## 0. Design decisions

1. **Pair-template estimator, no cells, no kappa map (report Sec. 5).** For every pair of pixels on distinct
   sightlines the response to a deflection field is exactly grad xi . [alpha(theta_a) - alpha(theta_b)] at first
   order, for any L. The amplitude estimator is Ahat = [sum w w d d R - b]/F with R = grad xi . dalpha_T, F = sum w w
   R^2, b = sum w w xi R. All pixel sums are done once per sightline pair, giving a catalogue of pair vectors
   V_ab (2), M_ab (3), B_ab (2); every template afterwards costs ~1e7 operations. This removes the cell window and the
   squeezed-limit sinc(L theta_ab/2) error of the cell estimator (review item 13), and the E/B window problems (15).
2. **Deprojection with a kernel-matched quasar template (report Eq. matched).** kappa_s-hat = sum over slab quasars of
   W_CMB(chi)/(b_q(chi) nbar_3D chi^2) minus mean. beta = 1 by construction; the only inputs are geometry and b_q(z).
   A single slab-averaged coefficient does not cancel arbitrary radial response weights (review item 4).
3. **Templates:** signal = Wiener-filtered (kappa_CMB - kappa_s-hat) for ACT DR6 and Planck PR4; response control =
   kappa_s-hat alone; curl null = rotated deflection; CMB null ensemble = ACT/Planck simulation maps as templates.
4. **Normalisation by injection:** shift the observed sightline coordinates of the data by -A alpha_inj, recompute the
   pair catalogue, recover A. This calibrates F through continuum distortion, weights and forest-length variations
   (review item 14). The analytic grad xi in the weights then affects only optimality.
5. **Tomographic sub-slabs** (e.g. 2.1-2.45, 2.45-2.8, 2.8-3.2) with disjoint template and sightline quasar
   selections where possible (review item 7); the single 2.1-3.0 slab is used only for forecasts.
6. **Validation:** (a) injection on data; (b) self-made joint mocks (lognormal forest + biased Poisson quasars +
   kappa fields from the same LSS, DR1 noise distribution, sparse sightlines) lensed by shifting sightlines, for
   recovery, deprojection and covariance; (c) later, DESI DR1 Lya mocks at NERSC with fixed-observed-geometry
   lensing (review item 16). Joint-covariance validation is in phase 1, not phase 2 (review item 17).

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
| 6 | `inject.py` | generate alpha_inj realisations, shift sightlines, rerun module 3, recover A; produces the response calibration and its scale dependence | one pair pass per injection (5-10 injections) |
| 7 | `mocks.py` | joint lognormal mocks (forest, quasars, kappa_lya, kappa_CMB), sparse sightlines with DR1 noise, lensed by sightline shifts; run modules 3-5; recovery, deprojection residual, covariance | one pair pass per mock (10-20 mocks) |
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

1. **Phase 1 (measurement + validation of the estimator):** modules 0-6 on DR1 x ACT DR6 and Planck PR4 with the
   single-slab configuration; injection calibration; jackknife + sim-template nulls; joint mocks for the deprojected
   amplitude covariance. Deliverable: Ahat for signal, response and curl templates with calibrated errors.
2. **Phase 2 (systematics):** tomographic sub-slabs with disjoint selections; magnification slope from q x LRG;
   DLA/metal/continuum tests; tSZ-deprojected ACT map; released N_L and mask transfer functions; b_q(z) from C^{qq}.
3. **Phase 3:** full C^-1 weighting (CG), DESI mocks at NERSC, paper.

## 5. Open items from review 1 carried forward

- Higher-order lensing corrections to the cross-response (Gaussian third-order terms; lensing bispectrum).
- Joint likelihood forecast with b_q(z), m and the response amplitude free.
- Validation of the sparse-sampling noise model (independent modes with P_1D/nbar; l_max cutoff) against the discrete
  estimator on mocks.
- Sightline-density / quasar-template correlation (shared quasars) modelled in the joint mocks.
- Radial transfer mismatch between density, tidal and quasar-RSD responses within a sub-slab.

## 6. For codex review 2 (estimator only)

Review report Sec. 5 and this file's Sec. 0: (1) the exactness claim for Eq. (Rij) and the first-order Gaussian
estimator (Ahat), its normalisation F and mean field b with diagonal weights; (2) the compression to per-pair
(V, M, B) and the statement that no per-sightline compression exists; (3) the Wiener-filter template choice as the
Fisher-optimal amplitude estimator in the diagonal-weight approximation; (4) the injection-calibration argument
(geometric response, continuum projection along the LOS); (5) the curl-template null and the conditional nature of the
sim-template null ensemble; (6) the kernel-matched template (Eq. matched) incl. shot noise and magnification; (7) any
bias that survives: same-sightline terms, sightline-density correlations, second-order deflection terms.
