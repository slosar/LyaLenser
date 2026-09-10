# IMPLEMENTATION: pair-template estimator pipeline

Exact specification for the implementer (revised after codex review 2, `report/reviews/codex_review_2.md`).
Symbols and equation labels refer to `report/main.tex` (Sec. 5 = estimator, Sec. 4.2 = matched template).
Code goes in `code/pipeline/`, tests in `code/pipeline/tests/`, data under `/data/LyaLenser/` (never in the repo).
Python 3.11 at `/home/anze/anaconda3/bin/python3` with numpy, scipy, numba, healpy, astropy, fitsio, h5py; if a
package is missing, create `code/pipeline/.venv` with `uv` and record it in `CLAUDE.md`. Reuse `code/cosmo.py`,
`code/forest_power.py`, `code/cross_spectrum.py` (`kernel`), `code/three_tracer.py` (`spectra`, `limber`); do not
reimplement them. Machine: 24 threads, 62 GB RAM; keep peak memory under 40 GB.

**Stage A (mocks only, no downloads) must be complete and its acceptance tests passing before Stage B (real data).**

---------------------------------------------------------------------------------------------------------------------
## 0. Conventions (fixed)

- Units: comoving Mpc/h; angles in radians internally (degrees in catalogues); `cosmo.chi(z)`, `cosmo.z_of_chi`.
- Positions: RA, Dec on the sphere. For a pair (a, b) the separation vector in the tangent plane at a:
  `dx = (ra_b - ra_a) cos(dec_mid)`, `dy = dec_b - dec_a` (radians; small-angle errors O(theta^2) < 1e-4 are
  acceptable for theta < 0.6 deg). `that_ab = (dx, dy)/theta` points from b to a. Deflections `alpha = (east, north)`
  in radians at each sightline position.
- Lensing convention: observed `dt(theta_obs) = delta(theta_obs + alpha(theta_obs))`, `alpha = grad phi`,
  `kappa = -nabla^2 phi/2`, on the sphere `phi_lm = 2 kappa_lm/(l(l+1))` (l >= 2; l = 0, 1 set to zero).
  A field known at source positions is observed at `theta_obs = theta_src - alpha(theta_src) + O(alpha^2)`, so
  **mock lensing and injections shift sightline positions by `-alpha`**.
- Pixel-pair geometry: `r_par = |chi_p - chi_q|`, `chi_mid = (chi_p + chi_q)/2`, `r_perp = chi_mid * theta_ab`.
- Forest model: `xi(r_perp, r_par)` and `xi_rp = d xi/d r_perp` on a regular grid `[0, 40] x [0, 40]` Mpc/h, step
  0.25, bilinear interpolation. The angular gradient of the pair covariance is `chi_mid * xi_rp * that_ab`
  (report Eq. Rij). `xi` is the *observed* (continuum-projected, weighted) mean correlation function; in Stage A it
  is measured from the mock skewers with the same code path used for data (Sec. 2.2).
- Weights: per pixel, inverse variance including intrinsic forest variance (picca-like). The estimator is diagonal
  in the weights.
- Reference plane: `chi_ref = chi(2.4)`; kernel ratio `g(chi) = W_lens(chi)/W_lens(chi_ref)` where
  `W_lens(chi_s) = <(chi_s - chi_l)/chi_s>` averaged over a fiducial lens distribution (use the mean lens distance of
  the CMB-lensing x forest-lensing kernel product, computed once from `three_tracer.spectra` ingredients); to first
  order `g(chi) ~ 1 + g1 (chi - chi_ref)`; store `g1`.

---------------------------------------------------------------------------------------------------------------------
## 1. Data structures (HDF5 unless stated)

### 1.1 `SightlineSet`
`qid int64[Nq]`, `ra, dec float64[Nq]` (deg), `zq float32[Nq]`, `pix_start int64[Nq+1]` (CSR; pixels sorted by chi),
`chi float32[Npix]`, `delta float32[Npix]`, `w float32[Npix]`, `slab int8[Npix]` (-1 unused). Attrs: `zmin, zmax,
description, chi_ref`.

### 1.2 `PairCatalogue` (one group per sub-slab: `slab0`, ..., or `all`)
For each unordered sightline pair a<b within `theta_max` with >= 1 accepted pixel pair:
- `a, b int32[Np]`, `thx, thy float32[Np]` (`that_ab`), `theta float32[Np]` (radians)
- `v float64[Np, NB]`, `m float64[Np, NB]`, `beta float64[Np, NB]`: the sums of report Eq. Vab,
  `v = sum w_p w_q d_p d_q chi_mid xi_rp`, `m = sum w_p w_q (chi_mid xi_rp)^2`, `beta = sum w_p w_q xi chi_mid xi_rp`,
  split into `NB = 6` shape bins: `r_perp` in {[0,10), [10,20), [20,30]} x `r_par` in {[0,10), [10,30]} Mpc/h.
- `v1, m1, beta1 float64[Np, NB]`: the same sums weighted by `(chi_mid - chi_ref)` (first moment) so that the
  redshift dependence of the deflection can be applied afterwards: effective `v_eff = v + g1 v1`, etc.
- `npair int32[Np]`.
Store float64 accumulators; the file may be written as float32 for `v, beta` (state the precision in attrs).
Same-sightline pixel pairs are never included. Pixel pairs are accepted if `r_par <= r_par_max` (30) and
`r_perp <= r_perp_max` (30). With sub-slabs a pixel pair contributes to the slab of `chi_mid`.

### 1.3 `Template`
`alpha float32[Nq, 2]` at the sightline positions; attrs: `name`, `kind` in {signal, response, curl, random,
injection, truth}, `Lmin, Lmax`, `filter` description, `source`. Also store the harmonic `phi_lm` used.

### 1.4 `AmplitudeResult`
For a set of templates {T_b}: scores `q_b = sum_pairs v_eff d_b`, response matrix `F_bc = sum_pairs m_eff d_b d_c`,
mean field `mf_b = sum_pairs beta_eff d_b`, with `d_b = thx (alpha_bx[a] - alpha_bx[b]) + thy (alpha_by[a] -
alpha_by[b])`; `A = F^{-1} (q - mf)`; `sigma_F = sqrt(diag(F^{-1}))` (lower bound); jackknife estimates from
per-region partial sums (region = HEALPix pixel of sightline a at `nside_jk`); all partial sums saved.

---------------------------------------------------------------------------------------------------------------------
## 2. Stage A modules

### 2.1 `pipeline/config.py`
Dataclass `Config`: `r_perp_max=30.0, r_par_max=30.0, xi_step=0.25, xi_max=40.0, shape_bins` (as in 1.2),
`nside_jk=8, nside_alpha=1024, lmax_alpha=2000, chi_ref=chi(2.4), slabs=[(2.1, 3.0)]`, paths, seeds.

### 2.2 `pipeline/xi_model.py`
Two providers with the same interface `XiTable(r_perp, r_par, xi, xi_rp)`:
- `xi_from_model(pf: ForestPower, cfg)`: analytic. `xi(r_perp, r_par) = int d^3k/(2pi)^3 P_F(k_par, k_perp)
  e^{i k.r}` = `int dk_par/(2pi) cos(k_par r_par) int k_perp dk_perp/(2pi) J_0(k_perp r_perp) P_F`; `xi_rp` with
  `-k_perp J_1`. Use a log grid in k with scipy Bessel functions; verify `xi(0,0)` against the direct variance
  integral and `xi_rp` against finite differences (1%).
- `xi_from_data(sl: SightlineSet, cfg)`: measure the mean observed correlation function on the
  `(r_perp, r_par)` grid from the sightlines themselves (standard pair count with weights, 1 Mpc/h bins, then smooth
  with a 2D Savitzky-Golay or spline and differentiate). This is the production path (report Sec. 5.2, normalisation
  through the measured xi). It needs its own numba pair kernel (share code with 2.3 via a mode flag).
- Both return numba-callable bilinear interpolators (`@njit` functions taking flat arrays).

### 2.3 `pipeline/pairs.py`
- `find_pairs(sl, theta_max) -> (a, b, thx, thy, theta)`: `scipy.spatial.cKDTree` on unit vectors, radius
  `2 sin(theta_max/2)`, `theta_max = r_perp_max / chi_min(slab)`; keep a<b.
- `accumulate(sl, pairs, xi_table, cfg, shifted_positions=None) -> PairCatalogue`: `@njit(parallel=True)` over
  pairs. Per pair: two-pointer sweep over the two chi-sorted pixel arrays with `|chi_p - chi_q| <= r_par_max`; per
  pixel pair compute `chi_mid`, `r_perp`, skip if `r_perp > r_perp_max`, interpolate `xi_rp` and `xi`, find the shape
  bin, accumulate the six quantities (0th and 1st chi moments). Flat arrays only inside the kernel. Results must be
  independent of the thread count.
- `benchmark(sl, fraction)`: pixel pairs per second per thread; extrapolated wall time; print and save.
- Tests: brute-force numpy on 20 sightlines x 50 pixels agrees to 1e-10 relative; 1 thread == 24 threads bitwise.

### 2.4 `pipeline/templates.py`
- `phi_from_kappa(kappa_alm, lmax)`: `phi_lm = 2 kappa_lm/(l(l+1))`, zero for l < 2.
- `alpha_at(positions, phi_lm, nside)`: `healpy.alm2map_der1` -> (map, d/dtheta, d/dphi / sin theta);
  `alpha_east = dphi/(sin theta) component`, `alpha_north = -dtheta component`; `healpy.get_interp_val` at the
  positions. Test: (i) for a low-l `phi_lm`, finite-difference divergence of the alpha map equals `-2 kappa` to a few
  per cent; (ii) flat-sky patch with `phi = phi0 cos(L x)`: `alpha_x = -phi0 L sin(L x)`, sign and amplitude.
- `wiener_filter(X_alm, S_L, C_XX_L)`: `h_L = S_L / C_XX_L`; band-limited versions with cosine-tapered top hats over
  `[Lmin, Lmax]`; model spectra from `three_tracer.spectra()` (signal `S_L = C^{kl kc} - C^{kl s}`; `C_XX` = total
  power of the map actually used including its noise).
- `curl(alpha)`: `(east, north) -> (-north, east)`.
- `matched_template(quasars, randoms, b_q_of_z, cfg) -> kappa_s_alm` (report Eq. matched): per-object weight
  `u_i = W_CMB(chi_i)/(b_q(z_i) nbar(chi_i) Omega_pix)` with `nbar(chi) = dN/(dchi dOmega)` from the smoothed
  data redshift distribution over the *template slab* (forest slab extended by 150 Mpc/h on each side); map =
  `sum_data u - (N_data/N_rand) sum_rand u` per pixel; pixels with fewer than `nmin_rand` randoms set to zero and
  recorded in a mask; `map2alm` with `lmax_alpha`. Shot noise `N_L^{ss} = int dchi W^2/(b_q^2 nbar_3D chi^2)` from
  the same `nbar`. Optional RSD/magnification corrections are *not* implemented in Stage A; the mock quantifies
  their size instead (2.7).
- `gaussian_realisations(Cls, nside, lmax, seed, n)`: correlated Gaussian fields from a covariance of spectra.

### 2.5 `pipeline/amplitude.py`
- `amplitude(cat, templates: list, g1, regions) -> AmplitudeResult` as in 1.4 (vectorised numpy; per-region partial
  sums of `q, F, mf` for jackknife; joint fits of gradient + curl templates by including the curl template in the
  list).
- `random_ensemble(cat, list_of_alpha)`: distribution of `A` for random templates, each with its own `mf` and `F`.
- `shape_test(cat, template)`: `A` per shape bin (uses `v[:, bin]`, `m[:, bin]`, `beta[:, bin]`) with jackknife
  errors; the six values must be mutually consistent for a lensing signal.

### 2.6 `pipeline/inject.py`
- `shift_positions(sl, alpha, A) -> SightlineSet` (`ra' = ra - A alpha_east/cos dec`, `dec' = dec - A alpha_north`).
- `injection_test(sl, xi_table, alpha_inj, A_list, cfg)`: for A in `A_list` (include +-A pairs and 0), shift, rebuild
  pairs and catalogue, fit the amplitude against `alpha_inj` and against its curl; report the linear slope from the
  paired +-A differences. **Purpose: bookkeeping, sign, band matrix and mean-field checks (report Sec. 5.5). It is
  not a physical calibration; do not describe it as one in outputs.**

### 2.7 `pipeline/mock.py` (Stage A data source; fixed observed geometry)
Flat-sky patch, `Lx = Ly = 20 deg` at (RA, Dec) = (180, 30), depth `chi(2.1) .. chi(3.0)` plus a 150 Mpc/h margin on
each side for the template slab; transverse grid 2 Mpc/h at `chi_ref`, radial grid 0.5 Mpc/h, float32; a `scale`
parameter shrinks the patch for tests.
1. `delta_m`: Gaussian, linear P(k, z=2.4) via FFT (single growth factor over the box is acceptable; document it).
2. `delta_F = b_F (1 + beta_F mu^2) delta_m` in Fourier space with the Kaiser+cutoff `ForestPower(model="kaiser")`
   shape, then the **response emulation** `delta_F -> delta_F (1 + (R_delta/2) delta_L)`, `delta_L` = `delta_m`
   smoothed with a 10 Mpc/h Gaussian, `R_delta = 2` (cfg).
3. Quasars: Poisson sample of `n_q (1 + b_q delta_L + f mu^2 term implemented as a redshift-space shift of the
   quasar chi by the LOS velocity from delta_m)`, clipped at zero, `b_q = 3.5`, over the template slab; a subset in
   the forest slab provides the sightlines with density `n_los` per deg^2 (default 22); each sightline covers the
   forest slab (simplification; `zq` = z at the far edge). Randoms: uniform, 20x. Magnification: optional flag that
   modulates the quasar density by `(1 + m kappa_q)` with `kappa_q` = the mock `kappa_lya` (same source plane),
   `m = 0.5`.
4. Convergences: `kappa_slab = int W_CMB(chi) delta_m dchi` from the box; `(kappa_lya, kappa_rest)` a correlated
   Gaussian pair with spectra `(C^{kl kl}, C^{kl kc} - C^{kl slab}, C^{kc kc} - C^{slab slab})` from
   `three_tracer.spectra`; `kappa_CMB = kappa_slab + kappa_rest`; CMB noise added as white noise at the level in
   `report/numbers3.json` when a noisy template is wanted. `alpha_lya` from `kappa_lya` (flat-sky FFT).
5. Lensing with fixed observed geometry: the observed sightline positions are the quasar positions; the skewer at
   observed `theta_obs` samples `delta_F` at `theta_obs + alpha_lya(theta_obs)` (trilinear interpolation), which is
   the exact first-order statement (single source plane in the mock; a `g(chi)` scaling option applies
   `alpha_lya * g(chi)` per pixel to test the first-moment correction).
6. Pixels of 0.55 Mpc/h; Gaussian pixel noise with per-skewer sigma log-normal (median such that
   `P_N = sigma^2 dchi = 0.33`, ln-scatter in `P_N` of 2); `w = 1/(sigma^2 + var(delta_F))`; continuum emulation:
   subtract each skewer's weighted mean and slope.
7. Save everything, including truths (`kappa_lya`, `alpha_lya` at sightlines, `kappa_slab`, `kappa_CMB`,
   `delta_L` at sightlines) and the quasar/random catalogues.
Per-mock cost target: field generation minutes; pair pass seconds to minutes (8800 sightlines).

### 2.8 `pipeline/run_mock_validation.py` (Stage A acceptance; commit `report/mock_validation.md` + json + figures)
1. Build `xi_from_data` on the mock (production path) and `xi_from_model`; compare; use `xi_from_data` below.
2. **Normalisation:** amplitude against the truth template `alpha_lya` (kind `truth`, filter = identity) on a
   noiseless high-density variant (`P_N = 0`, `n_los = 100`) and on the fiducial mock: `A = 1` within 10% and within
   the jackknife error respectively. Report the difference between `xi_from_data` and `xi_from_model`
   normalisations (this is the measured-xi systematic).
3. **Injection bookkeeping:** paired +-A injections on the fiducial mock: slope 1 within 5% relative to the
   truth-template normalisation of step 2 (i.e. the injection must reproduce the model response, not the physical
   one; state both numbers).
4. **Curl:** joint gradient+curl fit for the truth and signal templates; curl amplitude consistent with zero after
   accounting for `F_RT`.
5. **Deprojection:** amplitudes against (a) filtered `kappa_CMB`, (b) filtered `kappa_CMB - kappa_s_hat` with the
   matched template built from the mock quasars and randoms, (c) `kappa_s_hat` alone. Expect (a) biased high by the
   response term (predict it from the mock's `R_delta` and the catalogue), (b) `A = 1` within errors, (c) consistent
   with the injected response. Repeat (b) with the magnification flag on and with the template slab equal to the
   forest slab (no margin) to show the size of those effects.
6. **Random-template ensemble:** 100 Gaussian `kappa` realisations with the `kappa_CMB` spectrum; compare the std of
   `A` with the jackknife error and with `1/sqrt(F)`.
7. **Covariance and mean field:** >= 20 seeds; scatter of the step-5(b) amplitude vs per-mock errors; mean field
   check on unlensed mocks (`alpha_lya = 0`) with the real template: the mean of `A` must be zero within errors.
8. Benchmark line from 2.3.

---------------------------------------------------------------------------------------------------------------------
## 3. Stage B modules (real data)

### 3.1 `pipeline/fetch.py`
Download with resume and checksums into `/data/LyaLenser/raw/` and write `MANIFEST.json` (URL, size, md5):
- DESI DR1 Lya deltas: https://data.desi.lbl.gov/public/dr1/vac/dr1/lya-deltas/ (`delta-lya-0-0/Delta/Delta-*.fits.gz`
  and `delta_attributes`). Inspect one file first and document the HDU layout before writing the reader.
- DESI DR1 QSO catalogue and LSS randoms (public DR1 LSS release; record the exact URLs in MEMORY.md).
- DESI DR1 DLA catalogue.
- ACT DR6 lensing (LAMBDA): kappa alm (baseline and tSZ-deprojected), masks, N_L, simulations; follow the LAMBDA
  instructions on mask usage for cross-correlations.
- Planck PR4 lensing (Carron, Mirmelstein & Lewis 2022): klm, mask, N_L, sims.

### 3.2 `pipeline/desi_io.py`
`read_deltas(paths, zmin, zmax, cfg) -> SightlineSet` (wavelength -> z -> chi; keep `zmin <= z <= zmax`; DLA masks
+- dv; drop forests with < 50 pixels; per-pixel weights; sub-slab labels); `quasar_catalogue`, `randoms` -> (ra,
dec, z); also write `report/noise_distribution.json` (histogram of per-skewer `P_N`) for the forecast update.

### 3.3 `pipeline/cmb_io.py`
Readers for ACT DR6 and Planck PR4 alm, masks, N_L; simulation iterators; filtered `kappa_alm` given a filter and
the recommended mask treatment; total-power estimate `C_XX` of the map actually used (from the map itself).

### 3.4 `pipeline/run_dr1.py`
1. Single slab 2.1 < z < 3.0 end to end: deltas -> `xi_from_data` -> pairs -> catalogue (save; benchmark);
   templates: ACT signal, Planck signal, matched `kappa_s_hat` from DR1 quasars over the extended template slab with
   `b_q(z) = 0.278((1+z)^2 - 6.565) + 2.393`, curl of each, 400 ACT random templates; amplitudes with joint
   gradient+curl fits; bands [40,100], [100,200], [200,300] with the response matrix; shape test; injections (+-A,
   5 realisations); jackknife at nside 8.
2. Tomographic sub-slabs (2.1-2.45, 2.45-2.8, 2.8-3.2) with template quasars restricted to each sub-slab and
   sightlines from quasars whose forests cover it; the same outputs per sub-slab; the sightline-density x kappa
   diagnostic (cross-spectrum of the sightline count map with the ACT map).
Outputs: `report/dr1_run.json`, figures, and a markdown summary listing every number with its error and the tests
passed or failed.

---------------------------------------------------------------------------------------------------------------------
## 4. Acceptance
- `pytest code/pipeline/tests` passes (kernel vs brute force, thread determinism, xi tables, deflection sign and
  derivative, amplitude on a synthetic catalogue with a known `d`, per-sightline compression identity).
- Stage A: Sec. 2.8 items 2-7 within tolerance, summary committed.
- Every script runnable standalone from `code/pipeline/`; wall times logged; nothing written inside the repo except
  small json/markdown/figures.

## 5. Out of scope for the first implementation
Full C^-1 weighting; released N_L inside the Wiener filter (band-pass filters and map-estimated C_XX suffice);
non-Limber RSD/tidal transfer corrections (measured on the mock only); DESI mocks at NERSC; magnification slope
measurement from q x LRG.
