# IMPLEMENTATION: pair-template estimator pipeline

Exact specification for the implementer (revised after codex reviews 2 and 3, `report/reviews/codex_review_{2,3}.md`).
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
- Positions: RA, Dec on the sphere. For a pair (a, b) the separation vector is **theta_ab = theta_a - theta_b**,
  i.e. in the local tangent plane `dx = (ra_a - ra_b) cos(dec_mid)`, `dy = dec_a - dec_b` (radians; small-angle
  errors O(theta^2) < 1e-4 are acceptable for theta < 0.6 deg). `that_ab = (dx, dy)/theta` **points from b to a**.
  Deflections `alpha = (east, north)` in radians at each sightline position. The response of the pair is
  `xi_rp * that_ab . (alpha_a - alpha_b)`; for a constant convergence `alpha = -kappa theta` this gives
  `-kappa r xi_rp > 0` for decreasing xi, equal to `d/dkappa xi((1-kappa) r)` (Mathematica check 11). **A unit test
  must verify this sign by finite differences on a synthetic pair; reversing the direction flips the sign of A.**
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
  `W_lens(chi_s) = <(chi_s - chi_l)/chi_s>` averaged over a fiducial lens distribution (the mean lens distance of the
  CMB-lensing x forest-lensing kernel product, computed once from `three_tracer` ingredients); to first order
  `g(chi) ~ 1 + g1 (chi - chi_ref)`; store `g1`. A pixel p at chi_p is deflected by `g(chi_p) alpha_a`, so for a
  pixel pair (Mathematica check 12)
  `dalpha_pq = (alpha_a - alpha_b)(1 + g1 dm) + g1 (chi_p - chi_q)/2 (alpha_a + alpha_b)`, `dm = chi_mid - chi_ref`.
  With `G = chi_mid xi_rp`, `d = that.(alpha_a - alpha_b)`, `s = that.(alpha_a + alpha_b)`:
  `R = G [ d (1 + g1 dm) + g1 (chi_p - chi_q)/2 s ]` and
  `R^2 = G^2 [ d^2 (1 + 2 g1 dm + g1^2 dm^2) + g1 d s (chi_p - chi_q)(1 + g1 dm) + g1^2 s^2 (chi_p - chi_q)^2/4 ]`.

---------------------------------------------------------------------------------------------------------------------
## 1. Data structures (HDF5 unless stated)

### 1.1 `SightlineSet`
`qid int64[Nq]`, `ra, dec float64[Nq]` (deg), `zq float32[Nq]`, `pix_start int64[Nq+1]` (CSR; pixels sorted by chi),
`chi float32[Npix]`, `delta float32[Npix]`, `w float32[Npix]`, `slab int8[Npix]` (-1 unused). Attrs: `zmin, zmax,
description, chi_ref`.

### 1.2 `PairCatalogue` (one group per sub-slab: `slab0`, ..., or `all`)
For each unordered sightline pair a<b within `theta_max` with >= 1 accepted pixel pair:
- `a, b int32[Np]`, `thx, thy float32[Np]` (`that_ab`), `theta float32[Np]` (radians)
- Accumulators, each `float64[Np, NB]` with `NB = 6` shape bins (`r_perp` in {[0,10), [10,20), [20,30]} x
  `r_par` in {[0,10), [10,30]} Mpc/h), with `G = chi_mid xi_rp`, `dm = chi_mid - chi_ref`, `dc = chi_p - chi_q`:
  `v0 = sum w w d_p d_q G`, `v1 = sum ... G dm`, `vc = sum ... G dc/2`;
  `m0 = sum w w G^2`, `m1 = sum ... G^2 dm`, `m2 = sum ... G^2 dm^2`, `mc = sum ... G^2 dc`, `mcc = sum ... G^2 dc^2/4`;
  `beta0 = sum w w xi G`, `beta1 = sum ... xi G dm`, `betac = sum ... xi G dc/2`;
  (11 accumulators; the `mc dm` and `s^2` second-order cross terms are dropped). Then for a template with pair
  scalars `d, s`: score `q = sum [ (v0 + g1 v1) d + g1 vc s ]`, response `F = sum [ (m0 + 2 g1 m1 + g1^2 m2) d^2
  + g1 mc d s + g1^2 mcc s^2 ]`, mean field `mf = sum [ (beta0 + g1 beta1) d + g1 betac s ]`. For two templates
  b, c the response matrix uses `d_b d_c`, `(d_b s_c + d_c s_b)/2`, `s_b s_c` in the obvious way.
- `npair int32[Np]`.
Accumulate in float64; write float32 (`11 x 6 x 4 B x 1e7 pairs = 2.6 GB`; state the precision in attrs).
Same-sightline pixel pairs are never included. Pixel pairs are accepted if `r_par <= r_par_max` (30) and
`r_perp <= r_perp_max` (30). With sub-slabs a pixel pair contributes to the slab of `chi_mid`.

### 1.3 `Template`
`alpha float32[Nq, 2]` at the sightline positions; attrs: `name`, `kind` in {signal, response, curl, random,
injection, truth}, `Lmin, Lmax`, `filter` description, `source`. Also store the harmonic `phi_lm` used.

### 1.4 `AmplitudeResult`
For a set of templates {T_b}: pair scalars `d_b = that . (alpha_b[a] - alpha_b[b])`, `s_b = that . (alpha_b[a] +
alpha_b[b])`; scores `q_b`, response matrix `F_bc` and mean field `mf_b` from the formulas in 1.2;
`A = F^{-1} (q - mf)`; `sigma_F = sqrt(diag(F^{-1}))` is the Gaussian independent-pair approximation (neither an
upper nor a lower bound; report it only as a scale); jackknife estimates from per-region partial sums (region =
HEALPix pixel at `nside_jk` of the pair **midpoint**); all partial sums saved. The template list for any fit always
includes a "junk" band covering every multipole of the map outside the science bands (e.g. `[2,40)` and
`(300, lmax]`) so that the response of the science bands to omitted modes is marginalised rather than ignored, and
the curl template of each science band.

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
- Tests: (i) the `k_par` integral normalisation: `xi(0, 0)` from the table equals `int d^3k/(2pi)^3 P_F` computed
  independently (watch the factor 2 from integrating `k_par` over [0, inf) only); (ii) `xi_rp` vs finite differences;
  (iii) grid convergence (halve `xi_step`, change < 0.5% in `A` on a mock).
- Note on baseline fitting: `xi_from_data` fits a function of `(r_perp, r_par)` only, so its response to a
  band-limited deflection template is proportional to the footprint mean of the template's convergence, which is
  ~0 for L >= 40. The absorption of signal into the baseline is therefore expected to be negligible; it is measured
  (2.8 step 2, `xi_from_data` vs `xi_from_model` normalisations) rather than assumed.

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
- `alpha_at(positions, phi_lm, nside)`: `healpy.alm2map_der1(phi_lm, nside)` returns `(map, dtheta, dphi_over_sin)`
  where the third array is **already** `(1/sin theta) d/dphi`; do not divide by `sin theta` again.
  `alpha_east = dphi_over_sin`, `alpha_north = -dtheta` (theta is colatitude); `healpy.get_interp_val` at the
  positions. Test: (i) for a low-l `phi_lm`, finite-difference divergence of the alpha map equals `-2 kappa` to a few
  per cent; (ii) flat-sky patch with `phi = phi0 cos(L x)`: `alpha_x = -phi0 L sin(L x)`, sign and amplitude.
- `wiener_filter(X_alm, S_L, C_XX_L)`: `h_L = S_L / C_XX_L`; band-limited versions with cosine-tapered top hats over
  `[Lmin, Lmax]` (the taper is part of the band definition and enters `F_bc`). `S_L` from a parametrised version of
  `three_tracer.spectra(z1, z2, zs, ...)` evaluated for the **actual** forest slab, template slab and source plane
  (extend that function; do not use its hard-coded defaults). `C_XX_L` = the pseudo-`C_L` of the map actually used
  (masked, apodised, filtered), divided by the mask's `f_sky` factor, including its noise. On a masked sky
  `h_L X` is an approximation to the conditional mean `C_{kX} C_{XX}^{-1} X`; sightlines closer than 2 degrees to a
  mask edge are excluded from the amplitude (flag in the catalogue), and the residual normalisation error is
  measured on the mock with the real mask applied (2.7).
- `curl(alpha)`: `(east, north) -> (-north, east)`.
- `matched_template(quasars, randoms, b_q_of_z, cfg) -> kappa_s_alm` (report Eq. matched): per-object weight
  `u_i = W_CMB(chi_i)/(b_q(z_i) nbar(chi_i) Omega_pix)` with `nbar(chi) = dN/(dchi dOmega)` the mean radial density
  over the *template slab* (forest slab extended by 150 Mpc/h on each side); per pixel
  `kappa_s = [ sum_data u - (N_data/N_rand) sum_rand u ] / c_pix`, where the completeness `c_pix` is the random
  density in the pixel divided by its footprint mean, smoothed on 1 degree; pixels with `c_pix < 0.5` or fewer than
  `nmin_rand` randoms are set to zero and recorded in a mask; `map2alm` with `lmax_alpha`. Shot noise
  `N_L^{ss} = (1 + N_data/N_rand) int dchi W^2/(b_q^2 nbar_3D chi^2)`. RSD and magnification are not corrected in
  Stage A; the mock measures their size (2.7-2.8).
- `gaussian_realisations(Cls, nside, lmax, seed, n)`: correlated Gaussian fields from a covariance of spectra.

### 2.5 `pipeline/amplitude.py`
- `amplitude(cat, templates: list, g1, regions) -> AmplitudeResult` as in 1.4 (vectorised numpy; per-region partial
  sums of `q, F, mf` for jackknife; the template list always contains the science bands, their curl partners and
  the junk band(s), fitted jointly).
- `random_ensemble(cat, list_of_alpha)`: distribution of `A` for random templates, each with its own `mf` and `F`.
- `shape_test(cat, template)`: `A` per shape bin with jackknife errors; the six real-space values must be mutually
  consistent for a lensing signal (this replaces the Fourier-space D/T-vs-k test of the report, which the catalogue
  does not support).

### 2.6 `pipeline/inject.py`
- `shift_positions(sl, alpha, A) -> SightlineSet` (`ra' = ra - A alpha_east/cos dec`, `dec' = dec - A alpha_north`).
- `injection_test(sl, xi_table, alpha_inj, A_list, cfg)`: for A in `A_list` (include +-A pairs and 0), shift, rebuild
  pairs and catalogue, fit the amplitude against `alpha_inj` and against its curl; report the linear slope from the
  paired +-A differences. **Purpose: bookkeeping, sign, band matrix and mean-field checks (report Sec. 5.5). It is
  not a physical calibration; do not describe it as one in outputs.** The expected slope is the response of the
  *model* covariance (it equals the physical one only if `xi_from_data` equals the truth); compare it with the
  fixed-geometry remapping test of 2.8 step 3, which is the physical one.

### 2.7 `pipeline/mock.py` (Stage A data source; fixed observed geometry)
Flat-sky patch, `Lx = Ly = 20 deg` at (RA, Dec) = (180, 30), depth `chi(2.1) .. chi(3.0)` plus a 150 Mpc/h margin on
each side for the template slab; transverse grid 2 Mpc/h at `chi_ref`, radial grid 0.5 Mpc/h, float32; a `scale`
parameter shrinks the patch for tests.
1. `delta_m`: Gaussian, linear P(k, z=2.4) via FFT (single growth factor over the box is acceptable; document it).
2. `delta_F = b_F (1 + beta_F mu^2) delta_m` in Fourier space with the Kaiser+cutoff `ForestPower(model="kaiser")`
   shape, then the **response emulation** `delta_F -> delta_F (1 + (R_delta/2) delta_L)`, `delta_L` = `delta_m`
   smoothed with a 10 Mpc/h Gaussian, `R_delta = 2` (cfg). (The smoothing only defines which modes the forest
   responds to; the deprojection must remove them whatever the smoothing, because the template traces the same
   `delta_m`.)
3. Quasars: Poisson sample of a **lognormal** density `n_q exp(b_q delta_g - b_q^2 var(delta_g)/2)` with `delta_g` =
   `delta_m` on the 2 Mpc/h grid (unsmoothed, so that the quasar template traces the same field as `kappa_slab`),
   `b_q = 3.5`, over the template slab, with an RSD shift of each quasar's chi by the LOS velocity from `delta_m`
   (linear theory, `f = 0.97`); the sightlines are the quasars in the forest slab, thinned to `n_los` per deg^2
   (default 22); each sightline covers the forest slab (simplification; `zq` = z at the far edge). Randoms:
   uniform, 20x, with an optional completeness pattern `c(theta)` applied to both data and randoms to test the
   completeness division of 2.4. Magnification: flag that modulates the quasar density by `(1 + m kappa_q)` with
   `kappa_q` = the mock `kappa_lya`, `m = 0.5`.
4. Convergences: from the box, `kappa_slab = int W_CMB(chi) delta_m dchi` over the template slab and
   `kappa_lya_box = int W_lya(chi) delta_m dchi` over `chi < chi_ref` inside the box (lenses inside the slab in
   front of the source plane; this is what creates the intrinsic slab overlap and the magnification leakage).
   Outside the box, `(kappa_lya_rest, kappa_rest)` is a correlated Gaussian pair whose spectra are Limber integrals
   with the kernels restricted to `chi` outside the box (use `three_tracer.limber` with modified limits):
   `C^{kl,rest kl,rest}`, `C^{kl,rest kc,rest}`, `C^{kc,rest kc,rest}`. Then `kappa_lya = kappa_lya_box +
   kappa_lya_rest`, `kappa_CMB = kappa_slab + kappa_rest`; white CMB noise at the level in `report/numbers3.json`
   when a noisy template is wanted. `alpha_lya` from `kappa_lya` (flat-sky FFT). Optionally apply the real ACT
   footprint mask (rotated onto the patch) to `kappa_CMB` before filtering, to measure the mask-induced
   normalisation error of 2.4.
5. Lensing with fixed observed geometry: the observed sightline positions are the quasar positions; the skewer at
   observed `theta_obs` samples `delta_F` at `theta_obs + A_true * alpha_lya(theta_obs) * g(chi_pixel)` (trilinear
   interpolation), the exact first-order statement with the sign convention of Sec. 0; `A_true` (default 1) and the
   `g(chi)` option (default on, `g1` from Sec. 0) are parameters so that the physical response and the first-moment
   correction can be tested at several amplitudes.
6. Pixels of 0.55 Mpc/h; Gaussian pixel noise with per-skewer sigma log-normal (median such that
   `P_N = sigma^2 dchi = 0.33`, ln-scatter in `P_N` of 2); `w = 1/(sigma^2 + var(delta_F))`; continuum emulation:
   subtract each skewer's weighted mean and slope.
7. Save everything, including truths (`kappa_lya`, `alpha_lya` at sightlines, `kappa_slab`, `kappa_CMB`,
   `delta_L` at sightlines) and the quasar/random catalogues.
Per-mock cost target: field generation minutes; pair pass seconds to minutes (8800 sightlines).

### 2.8 `pipeline/run_mock_validation.py` (Stage A acceptance; commit `report/mock_validation.md` + json + figures)
Predeclared tolerances are in brackets; a failure blocks Stage B.
1. **Tables and benchmark.** Build `xi_from_data` (production) and `xi_from_model` on the fiducial mock; report the
   benchmark of the pair pass.
2. **Baseline absorption.** Normalisation of step 3 with `xi_from_data` vs `xi_from_model` [agree to 3%].
3. **Physical normalisation (fixed geometry).** Mocks lensed with `A_true` in {0, 0.5, 1, 2} (same seed), truth
   template `alpha_lya` (kind `truth`, unfiltered), joint fit with curl and junk bands: recovered slope
   `dA_hat/dA_true` [1 +- 0.05 on the noiseless high-density variant `P_N = 0, n_los = 100`; 1 within the jackknife
   error on the fiducial mock]. Repeat with the `g(chi)` option off and the first-moment terms disabled, and on with
   them enabled [both 1 +- 0.05; with `g` on and moments disabled the slope must deviate by the predicted amount].
4. **Mean field.** Unlensed mocks (`A_true = 0`), fixed templates (truth, filtered `kappa_CMB`, `kappa_s_hat`):
   mean of `A_hat` over >= 20 seeds [0 within 2 sigma of the mean's error]; compare with `mf/F` to show the
   subtraction is doing the work.
5. **Injection bookkeeping.** Paired +-A coordinate shifts on the fiducial mock; slope relative to the model response
   [1 +- 0.05]; curl amplitude [0 within errors].
6. **Stochastic normalisation.** Signal template from the noisy `kappa_CMB` (Wiener filter, science bands + junk +
   curl) on >= 20 seeds with the response emulation OFF: mean `A_hat` [1 within 2 sigma of the mean's error].
7. **Response-only null.** Response emulation ON, lensing OFF (`A_true = 0`): (a) filtered `kappa_CMB` template gives
   the predicted response bias (prediction from `R_delta`, the mock spectra and the catalogue response, computed
   independently before the fit) [agree within errors]; (b) `kappa_CMB - kappa_s_hat` template [0 within errors];
   (c) `kappa_s_hat` template [consistent with the predicted response amplitude].
8. **Combined recovery.** Response ON, lensing ON, magnification flag ON, completeness pattern ON, real mask ON:
   deprojected `A_hat` [1 within errors over >= 20 seeds]; report the shifts when each flag is switched off, and
   with the template slab equal to the forest slab (no margin).
9. **Shape test** on step 8 [six bins consistent, chi^2 p-value > 0.01].
10. **Random-template ensemble.** 100 Gaussian `kappa` realisations with the `kappa_CMB` spectrum on one mock;
    compare std(`A_hat`) with the jackknife error and `sigma_F` (report the ratios; no tolerance).
11. **Covariance.** Scatter of the step-8 amplitude across seeds vs the per-mock jackknife error [ratio within
    1 +- 0.3, given ~20 seeds].

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
1. Single slab 2.1 < z < 3.0 end to end (first run only; the tomographic run is the science configuration):
   deltas -> `xi_from_data` -> pairs -> catalogue (save; benchmark);
   templates: ACT signal, Planck signal, matched `kappa_s_hat` from DR1 quasars over the extended template slab with
   `b_q(z) = 0.278((1+z)^2 - 6.565) + 2.393`, curl of each, 400 ACT random templates; amplitudes with joint
   gradient+curl fits; bands [40,100], [100,200], [200,300] with the response matrix; shape test; injections (+-A,
   5 realisations); jackknife at nside 8.
2. Tomographic sub-slabs (2.1-2.45, 2.45-2.8, 2.8-3.2; read deltas over 2.1 < z < 3.2 for this run): forest pixels
   of a sub-slab from all quasars whose forests cover it; the matched template for a sub-slab from quasars in the
   sub-slab extended by 150 Mpc/h on each side, **excluding the quasars that provide sightlines to that sub-slab**
   (disjoint selection; the excluded fraction is reported); the same outputs per sub-slab; science bands from
   L = 100 for the clean range (RSD) with the 40-100 band reported separately; the sightline-density x kappa
   diagnostic (cross-spectrum of the sightline count map with the ACT map).
Outputs: `report/dr1_run.json`, figures, and a markdown summary listing every number with its error and the tests
passed or failed.

---------------------------------------------------------------------------------------------------------------------
## 4. Acceptance
- `pytest code/pipeline/tests` passes (kernel vs brute force, thread determinism, xi tables incl. the k_par
  normalisation, the pair-direction sign test by finite differences, `alm2map_der1` convention test, amplitude on a
  synthetic catalogue with a known `d` and `s`, per-sightline compression identity, first-moment identity).
- Stage A: Sec. 2.8 items 2-11 within the bracketed tolerances, summary committed.
- Every script runnable standalone from `code/pipeline/`; wall times logged; nothing written inside the repo except
  small json/markdown/figures.

## 5. Out of scope for the first implementation
Full C^-1 weighting; released N_L inside the Wiener filter (band-pass filters and map-estimated C_XX suffice);
non-Limber RSD/tidal transfer corrections (measured on the mock only); DESI mocks at NERSC; magnification slope
measurement from q x LRG.
