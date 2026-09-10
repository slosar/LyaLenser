# PLAN: measuring Lyman-alpha forest lensing x CMB lensing in DESI DR1

Written 2026-09-10 after steps 1-3 (see report/main.pdf). Status: proposal, awaiting codex review.

## 0. Design decisions (and why)

1. **Real-space, cell-based estimator; no kappa map.** In each sky cell (HEALPix nside 64-128, 0.5-0.9 deg) measure the
   local cross-sightline correlation function of the forest, xi_cell(r_perp, r_par, phi), keeping the transverse pair
   angle phi. Its m=0 part gives a scalar "dilation" field, its m=2 part (cos 2phi, sin 2phi) a spin-2 "forest shear"
   field. Lensing predicts d xi = -(kappa + gamma cos 2phi) r_perp d xi_bar/d r_perp, i.e. the SAME radial template
   T(r) = -r_perp d xi_bar/d r_perp for both harmonics, built from the measured mean xi_bar. This is exact at linear
   order for sparse sightlines (no continuum-limit / aliasing approximations), weights are the picca pixel weights
   (C^-1-like), and the spin-2 field gives an E/B split: lensing (and the tidal response) are pure E, so B x kappa is a
   null test.
2. **Response deprojection with the quasars in the slab** (report Sec. 4): D_L = C_L^{X kappa} - beta C_L^{X q} for
   X = dilation, shear-E, with beta = C^{q kappa}/(C^{qq} - 1/nbar) cross-checked against the geometric value
   W_CMB Delta_chi / b_q. No forest response model is needed.
3. **CMB kappa from both ACT DR6 and Planck PR4**; ACT for the higher S/N per mode over its footprint, Planck for the full
   DESI footprint and as a systematic cross-check; use their public simulation suites for the null distribution.
4. **Validation by lensing the mocks ourselves**: lensing only remaps transverse sightline positions, so a mock is
   lensed by shifting each sightline by -alpha(theta) from a Gaussian kappa realisation correlated with a synthetic
   kappa_CMB. No new hydro or N-body runs are needed.

## 1. Data (all public)

| Data | Source | Size (est.) |
|---|---|---|
| DESI DR1 Lya deltas (picca, per-HEALPix `Delta-*.fits.gz`, regions Lya/Lyb/CIII) | https://data.desi.lbl.gov/public/dr1/vac/dr1/lya-deltas/ (doc: https://data.desi.lbl.gov/doc/releases/dr1/vac/lya-deltas/) | 10-40 GB (verify) |
| DESI DR1 QSO catalogue + LSS randoms (angular selection for q) | DESI DR1 LSS catalogues | few GB |
| DESI DR1 DLA catalogue (masking) | https://data.desi.lbl.gov/doc/releases/dr1/vac/dla-cnn-gp/ | small |
| ACT DR6 lensing: kappa alm, mask, N_L, 400 sims (baseline and tSZ-deprojected) | LAMBDA, ACT DR6 lensing products | ~30 GB with sims |
| Planck PR4 lensing: klm, mask, N_L, sims | PLA / LAMBDA (Carron, Mirmelstein & Lewis 2022) | ~10-20 GB with sims |
| (phase 2) low-z galaxy sample for quasar magnification slope, e.g. DESI Legacy LRGs | public | few GB |

Store under `/data/LyaLenser/` (not in Dropbox); record paths in MEMORY.md.

## 2. Pipeline modules (`code/pipeline/`)

| # | Module | What it does | Output | Cost here |
|---|---|---|---|---|
| 0 | `fetch.py` | download + checksum the data above | `/data/LyaLenser/raw` | I/O only |
| 1 | `prep.py` | read deltas: (RA, Dec, z_q, lambda, delta, weight); lambda -> chi; z-slab cuts (start: 2.1<z<3.0, then two slabs); DLA/metal masks; per-cell sightline lists (KD-tree) | HDF5, ~2 GB (1.7e8 pixels, float32) | minutes |
| 2 | `pairs.py` (numba, parallel) | for every pair of pixels on distinct sightlines with r_perp<30, abs(r_par)<30 Mpc/h: accumulate per cell, per (r_perp, r_par) bin, per slab: S0=sum w_i w_j d_i d_j, S2c=sum ... cos2phi, S2s=sum ... sin2phi, W, W2c, W2s (geometry of the sightline pattern); also global xi_bar(r_perp, r_par) and its r_perp derivative | arrays cells x bins x slabs, <1 GB | ~3e11 pair ops: 1-3 h on 24 threads |
| 3 | `fields.py` | per-cell fields: d xi_cell = xi_cell - xi_bar (m=0), (xi_2c, xi_2s) (m=2, geometry-corrected); template projection with the jackknife bin covariance -> scalar map kappa_d(theta) and spin-2 map gamma(theta); keep unprojected bins for the shape test | HEALPix maps | seconds |
| 4 | `longmodes.py` | kappa_CMB maps at nside 256 from alm, masks apodised; q overdensity map from slab quasars and randoms; beta_L from C^{q kappa}/(C^{qq}-shot) and from geometry (b_q from C^{qq}) | maps | minutes |
| 5 | `spectra.py` (pymaster) | pseudo-C_L: kappa_d x kappa, gamma_E x kappa, gamma_B x kappa (null), same x q, q x kappa, q x q; L bins of ~30 from 30 to ~300 (set by cell size); deprojection D_L; amplitude fit A_lens against the template across bins | tables | minutes |
| 6 | `covariance.py` | (a) jackknife over ~150 nside-8 regions from stored per-cell sums (no recomputation); (b) ACT/Planck kappa sims through step 5 for the null distribution; (c) rotations of kappa relative to the footprint | covariance, null distributions | ~1 h |
| 7 | `nulls.py` | B x kappa; sims; splits by quasar magnitude, redshift, region, NGC/SGC; D_L/T shape consistency across bins; response term B_FFq compared with Karacayli+2024 and Chiang+2017; sightline-density x kappa (mean field) | plots | minutes |
| 8 | `mocks.py` | self-made mock: lognormal delta_F on a 3D grid + Poisson quasars biased on the same field + kappa_lya, kappa_CMB from the same LSS (Gaussian, correct spectra); sparse Poisson sightlines with DR1 noise distribution; lens by shifting sightlines by -alpha; run modules 2-7; check A_lens recovery and deprojection. Later: DESI DR1 Lya mocks (NERSC) lensed the same way | | 1 mock ~ 1 run of the pipeline (1-3 h); 10-20 mocks over a few days |

Dependencies to add: `pymaster` (NaMaster; easiest via conda-forge `namaster`, else a `.venv`), `fitsio`/`astropy`, `healpy`, `numba` (present), `h5py`.

## 3. Resource assessment

- **This machine (Ryzen 9 7900, 24 threads, 62 GB RAM, ~880 GB free disk, RTX 3050 6 GB): sufficient for phases 1-2.**
  The dominant cost is the pair accumulation (~3e11 pair operations; numba at ~1e8 s^-1 per thread gives 1-3 h wall
  time), everything downstream is minutes. Memory: pixel arrays ~2 GB, per-cell sums <1 GB, maps negligible. Disk:
  <100 GB including simulation suites.
- **24 GB RTX GPU: not required.** It would accelerate the pair counting by 10-50x (custom CUDA kernel) and would be
  the natural platform for a phase-3 FFT-based optimal quadratic estimator on 3D grids (1e8-cell grids per 10x10 deg
  patch), but neither is on the critical path for a first measurement.
- **NERSC (A100 or CPU): not required for the measurement.** Needed only for (a) the official DESI DR1 Lya mocks
  (LyaCoLoRe/Saclay, ~TB, at NERSC) if we want validation beyond the self-made mocks, (b) regenerating deltas with
  different continuum-fitting choices via picca on the ~700k spectra (CPU nodes, hours), (c) DESI-internal catalogues
  if the public randoms prove insufficient. GPUs are not needed at any stage of this plan.

## 4. Phases

1. **Phase 1 (measurement, ~2 weeks of work):** modules 0-6 on DR1 x ACT DR6 and Planck PR4; deliverable: D_L for
   dilation and shear-E, B-mode null, A_lens with jackknife and sim errors, response term B_FFq.
2. **Phase 2 (validation and systematics):** module 8 self-made mocks; magnification slope from q x low-z galaxies;
   DLA/metal/continuum tests; tSZ-deprojected ACT map; two redshift slabs.
3. **Phase 3 (if a signal is present):** optimal FFT/QE variant for higher L, DESI mocks at NERSC, paper.

## 5. Known risks (from steps 1-3)

- Continuum fitting distorts xi at small r_par and low k_par (picca distortion matrix); the template T is built from the
  distorted xi_bar, which is self-consistent to first order, but the mock test must confirm it.
- Sightline density correlates with q (template quasars also provide pixels in the lower slab) and with foreground
  structure through quasar magnification: use template quasars from the slab and pixels from quasars behind it; monitor
  the sightline-density x kappa cross-spectrum.
- Quasar magnification enters beta at the ~10-20% level: needs the slope 5s-2 (phase 2).
- Cell size trades L_max against pair statistics; the expected signal is at L < 300 anyway.
- Expected significance is ~2 sigma (report Table 2); the response term should appear at ~7 sigma and is a useful
  positive control.

## 6. For the codex review (after restart)

Review targets, in order: (1) report Sec. 2.2-2.6 derivations and factors of 2 in N_kappa; (2) Sec. 4 deprojection
identity and its assumptions, especially the claim that the tidal response needs no separate template for k_par=0 modes
and the treatment of quasar magnification; (3) the noise model N=(P_N+P_1D)/n_eff with unchanged response and the
harmonic-mean weighting; (4) the real-space template claim in Sec. 0 above (same radial template for m=0 and m=2, exact
for sparse sightlines); (5) the resource estimate (pair-count scaling) and the mock-lensing-by-sightline-shift idea.
