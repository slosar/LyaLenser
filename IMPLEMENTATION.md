# IMPLEMENTATION: pair-template estimator pipeline (DRAFT, pending codex review 2)

This is the exact specification handed to the implementer. Conventions, symbols and equation numbers refer to
`report/main.tex` (Sec. 5 = pair-template estimator, Sec. 4.2 = kernel-matched template). Code goes in
`code/pipeline/`, tests in `code/pipeline/tests/`, data under `/data/LyaLenser/` (never in the repo). Python 3.11,
numpy, scipy, numba, healpy, astropy, fitsio, h5py; use `/home/anze/anaconda3/bin/python3`; if a package is missing,
create `code/pipeline/.venv` with `uv` and record it in `CLAUDE.md`.

Two stages. **Stage A is self-contained and must be finished and tested first** (mock data only, no downloads).
Stage B adds the real-data readers and templates.

---------------------------------------------------------------------------------------------------------------------
## 0. Conventions (fixed; do not change)

- Units: comoving distances in Mpc/h; angles in radians internally (degrees only in catalogues); redshift z.
- Cosmology: `code/cosmo.py` (CAMB, Planck-2018-like). Reuse `chi(z)`, `z_of_chi`, `linear_pk_interp`.
- Sky geometry: full-sky positions (RA, Dec). Pair separation vector in the local tangent plane at sightline a:
  `theta_ab = (dx, dy)` with `dx = (ra_b - ra_a) * cos(dec_mid)`, `dy = dec_b - dec_a` (radians, small-angle;
  errors O(theta^2) ~ 1e-5 are negligible for theta < 0.5 deg). The pair direction is `that_ab = theta_ab/|theta_ab|`,
  pointing from b to a. Deflection vectors `alpha = (alpha_x, alpha_y)` in the same local (east, north) basis at each
  quasar.
- Lensing convention (report Sec. 2.1): observed field `dt(theta_obs) = delta(theta_obs + alpha(theta_obs))`,
  `alpha = grad phi`, `kappa = -nabla^2 phi / 2`; on the sphere `phi_lm = 2 kappa_lm / (l(l+1))`. A field known at
  source positions `theta_src` is observed at `theta_obs = theta_src - alpha(theta_src) + O(alpha^2)`.
  **Injection and mock lensing therefore shift sightline positions by `-alpha`.**
- Pixel pair geometry: for pixels p (on a, at chi_p) and q (on b, at chi_q): `r_par = |chi_p - chi_q|`,
  `chi_mid = (chi_p + chi_q)/2`, `r_perp = chi_mid * |theta_ab|`.
- Forest model: mean correlation function `xi(r_perp, r_par)` (isotropic in the transverse plane) and its derivative
  `xi_rp(r_perp, r_par) = d xi / d r_perp`, tabulated on a regular grid `r_perp, r_par in [0, 40]` Mpc/h, step
  0.25 Mpc/h, bilinear interpolation. The transverse angular gradient of the pair covariance is
  `grad_theta xi = chi_mid * xi_rp(r_perp, r_par) * that_ab` (report Eq. Rij).
- Weights `w` are per pixel (inverse variance including intrinsic forest variance, i.e. picca weights); the estimator
  uses diagonal weights only.

---------------------------------------------------------------------------------------------------------------------
## 1. Data structures

### 1.1 Sightline set (`SightlineSet`, HDF5 file)
- `qid  int64 [Nq]`  unique id
- `ra, dec  float64 [Nq]` degrees
- `zq  float32 [Nq]`
- `pix_start int64 [Nq+1]` CSR offsets into the pixel arrays (pixels of each sightline sorted by increasing chi)
- `chi  float32 [Npix]`, `delta float32 [Npix]`, `w float32 [Npix]`
- optional `slab int8 [Npix]` sub-slab index of each pixel (-1 = unused)
- attrs: `zmin, zmax, description`

### 1.2 Pair catalogue (`PairCatalogue`, HDF5 file), one group per sub-slab (`slab0`, `slab1`, ... or `all`)
For every unordered sightline pair (a<b) with `|theta_ab| * chi_ref <= r_perp_max` and at least one contributing pixel
pair:
- `a, b  int32 [Np]`
- `thx, thy float32 [Np]`  components of `that_ab` (unit vector, from b to a)
- `v  float64 [Np]`  = sum_{p,q} w_p w_q delta_p delta_q chi_mid xi_rp(r_perp, r_par)
- `m  float64 [Np]`  = sum_{p,q} w_p w_q (chi_mid xi_rp)^2
- `bm float64 [Np]`  = sum_{p,q} w_p w_q xi(r_perp, r_par) chi_mid xi_rp
- `npair int32 [Np]` number of pixel pairs used
Because grad xi is along `that_ab`, the vectors of report Eq. Vab are `V_ab = v that_ab`, `M_ab = m that that^T`,
`B_ab = bm that_ab`; storing scalars is exact.

Pixel pairs are included if `r_par <= r_par_max` (default 30 Mpc/h) and `r_perp <= r_perp_max` (default 30 Mpc/h);
same-sightline pairs are never included. If a sub-slab assignment is used, a pixel pair contributes to the slab of
`chi_mid`.

### 1.3 Template (`Template`)
- `alpha  float32 [Nq, 2]` deflection at each sightline position (east, north), radians
- attrs: `name, kind in {signal, response, curl, sim, injection}, Lmin, Lmax, source map description`

### 1.4 Amplitude result
For a template: `num = sum_pairs v * d`, `F = sum_pairs m * d^2`, `mf = sum_pairs bm * d`, with
`d = thx*(alpha_x[a]-alpha_x[b]) + thy*(alpha_y[a]-alpha_y[b])`; `A = (num - mf)/F`; `sigma_F = 1/sqrt(F)` (Gaussian
diagonal approximation); jackknife mean and error; per-region sums stored for reuse.

---------------------------------------------------------------------------------------------------------------------
## 2. Modules (Stage A)

### 2.1 `pipeline/config.py`
Dataclass `Config` with all tunables and defaults: `r_perp_max=30.0`, `r_par_max=30.0`, `xi_grid_step=0.25`,
`xi_grid_max=40.0`, `nside_jk=8`, `slabs=[(2.1,3.0)]`, paths. Loaded from a small YAML/TOML or python file.

### 2.2 `pipeline/xi_model.py`
`build_xi_table(cfg) -> XiTable` with arrays `r_perp, r_par, xi, xi_rp` on the grid. For Stage A build it from the
same model used to generate the mock (linear Kaiser P_F from `code/forest_power.py`, `model="kaiser"` or `"arinyo"`):
xi(r_perp, r_par) by 2D Hankel/FFT of P_F(k_par, k_perp) on a fine grid, then `xi_rp` by finite differences on a finer
grid and downsampling, or analytically from the k-space integral with a J_1 Bessel kernel. Provide
`XiTable.interp(r_perp, r_par)` and `XiTable.interp_rp(...)` as pure-numpy bilinear interpolators AND numba-jitted
scalar versions for the kernel. Test: xi at r=0 equals the variance integral of P_F; xi_rp matches finite differences
of xi to 1%.

### 2.3 `pipeline/pairs.py`
- `find_pairs(sl: SightlineSet, r_perp_max, chi_ref_min) -> (a, b, thx, thy, theta)`: cKDTree on unit vectors
  with search radius `2 sin(theta_max/2)`, `theta_max = r_perp_max/chi_ref_min` (chi at the lowest usable redshift);
  keep a<b.
- `accumulate(sl, pairs, xi: XiTable, cfg) -> PairCatalogue`: numba `@njit(parallel=True)` over pairs (`prange`).
  For each pair: two-pointer sweep over the two sorted chi arrays keeping `|chi_p - chi_q| <= r_par_max`; for each
  pixel pair compute `r_perp = chi_mid * theta`, skip if `r_perp > r_perp_max`, bilinear-interpolate `xi_rp` and
  `xi`, accumulate `v, m, bm, npair` (float64 accumulators). No python objects inside the kernel; pass flat arrays.
- `benchmark(sl_subset)`: report pixel pairs per second per thread and the extrapolated wall time for the full set.
- Determinism: results must not depend on thread count (accumulators are per pair, so this is automatic).
- Test: brute-force numpy double loop on ~20 sightlines with ~50 pixels each must agree with the kernel to 1e-10
  relative; a run with 1 thread and with 24 threads must be identical.

### 2.4 `pipeline/templates.py`
- `alpha_from_kappa_alm(kappa_alm, lmax, nside, positions) -> alpha [Nq,2]`: `phi_lm = 2 kappa_lm / (l(l+1))`
  (set l=0,1 to zero), `healpy.alm2map_der1(phi_lm, nside)` gives (phi, dphi/dtheta, dphi/dphi_az / sin theta);
  convert to (east, north): `alpha_east = d/dphi_az/sin(theta)`, `alpha_north = -d/dtheta`; interpolate to the
  sightline positions with `healpy.get_interp_val`. Use nside >= 1024 for lmax <= 2000. Test: for a pure
  `kappa_lm` monopole-free low-l field, check `div alpha = -2 kappa` numerically on a map (finite differences) to a
  few per cent, and check the sign convention against an explicit `phi = phi0 cos(L x)` flat-sky patch.
- `filter_alm(kappa_alm, f_L)`: multiply by an L-dependent filter (band-pass `[Lmin, Lmax]` top-hat with cosine
  edges, or the Wiener filter from model spectra). Model spectra come from `code/three_tracer.spectra()`
  (reuse, do not reimplement).
- `curl_template(alpha)`: rotate every vector by +90 degrees (east, north) -> (-north, east).
- `matched_template(quasar_cat, randoms, b_q_of_z, nside_map, cfg) -> kappa_s alm`: per-object weight
  `u_i = W_CMB(chi_i) / (b_q(z_i) nbar(chi_i) Omega_pix)`, with `W_CMB` from `code/cross_spectrum.kernel` and
  `nbar(chi) = dN/(dchi dOmega)` from the smoothed data distribution; map = sum over data of u minus
  (N_data/N_rand) sum over randoms of u, pixels with too few randoms masked to zero; then `map2alm`. Document the mask
  handling. (Stage A: the mock provides quasar positions and a uniform random catalogue.)
- `gaussian_kappa_realisation(Cl, nside, seed)` for injections and for the mock (Stage A).

### 2.5 `pipeline/amplitude.py`
- `amplitude(cat: PairCatalogue, alpha, regions=None) -> AmplitudeResult` implementing Sec. 1.4 with numpy
  (vectorised; 1e7 pairs is fine). `regions`: HEALPix region index of sightline a at `nside_jk`; compute per-region
  partial sums of `num, F, mf` so that leave-one-out jackknife amplitudes and their covariance across a list of
  templates are obtained without recomputation.
- `band_amplitudes(cat, kappa_alm, bands, ...)`: loop over L bands, build filtered templates, return A_b, F_b,
  jackknife errors.
- `null_ensemble(cat, list_of_alpha)`: amplitudes for many sim templates; returns the array and its std.

### 2.6 `pipeline/inject.py`
- `shift_positions(sl, alpha) -> SightlineSet` with `ra' = ra - alpha_east/cos(dec)`, `dec' = dec - alpha_north`
  (radians -> degrees). Everything else unchanged.
- `injection_test(sl, xi, alpha_inj, A_list, cfg)`: for each A in `A_list` (e.g. [0, 2, 5, 10] times a realistic
  amplitude), shift by `-A alpha_inj`, rebuild pairs and catalogue, measure the amplitude against `alpha_inj`; return
  the recovered `A_hat(A)` and the linear calibration slope. Also measure the curl-template amplitude on the same
  injected sets (must be consistent with zero).

### 2.7 `pipeline/mock.py` (Stage A data source)
A joint flat-sky mock on a patch, deliberately simple but containing every effect the estimator must handle:
- Patch: `Lx = Ly = 20 deg` centred at (RA, Dec) = (180, 30) degrees, depth `chi in [chi(2.1), chi(3.0)]`;
  transverse grid 2 Mpc/h at `chi_mid`, radial grid 0.5 Mpc/h (float32 arrays, ~3 GB; allow a `scale` parameter to
  shrink for tests).
- Matter field: Gaussian `delta_m` with linear P(k, z=2.4) (`cosmo.linear_pk_interp`), generated by FFT.
- Forest field: `delta_F = b_F (1 + beta_F mu^2) delta_m` in Fourier space with a Gaussian small-scale cutoff
  (use `forest_power.ForestPower(model="kaiser")` parameters), then **response emulation**
  `delta_F -> delta_F (1 + (R_delta/2) delta_L)` with `delta_L` = `delta_m` smoothed with a 10 Mpc/h Gaussian and
  `R_delta = 2` (cfg). This inserts the long-mode response that the deprojection must remove.
- Quasars: Poisson sample with density `n_q(1 + b_q delta_L)` clipped at 0, `b_q = 3.5`, mean 2D density such that
  the sightline density is `n_los` per deg^2 (default 22) and every sightline spans the full depth (simplification;
  a `zq` is assigned = z at the far edge). Randoms: uniform, 20x the quasars.
- Convergences (2D Gaussian fields on the patch, generated jointly with the correct cross-spectra from
  `code/three_tracer.spectra()`): `kappa_lya`, `kappa_rest` (the part of `kappa_CMB` from outside the slab), and
  `kappa_slab = int dchi W_CMB(chi) delta_m` computed from the box; `kappa_CMB = kappa_slab + kappa_rest`;
  `kappa_lya` correlated with `kappa_rest` with coefficient set by the model spectra (the slab does not lens itself:
  `kappa_lya` is independent of the box). Provide the deflection `alpha_lya` from `kappa_lya`.
- Lensing: sightline positions are shifted by `-alpha_lya(theta)` (source -> observed).
- Skewers: sample `delta_F` along each sightline at 0.55 Mpc/h pixels with trilinear interpolation; add Gaussian
  pixel noise with per-skewer `sigma` drawn log-normally (median `sigma_delta = 0.78 * sqrt(0.55/dchi_pix)`
  scaled so that `P_N = sigma^2 dchi = 0.33` Mpc/h at the median, ln-scatter `sigma_ln/2` in sigma); weights
  `w = 1/(sigma^2 + sigma_LSS^2)` with `sigma_LSS^2 = var(delta_F)`; continuum-fitting emulation: subtract from
  each skewer its weighted mean and slope in chi (projection P).
- Outputs: `SightlineSet`, quasar catalogue + randoms, `kappa_CMB`, `kappa_lya`, `alpha_lya` (on the grid and at
  the sightlines), all truths saved. A CMB noise realisation is added to `kappa_CMB` when a noisy template is wanted
  (white N_L from the calibrated values in `report/numbers3.json`).
- Provide `make_mock(seed, cfg) -> paths`; 10-20 seeds are needed for covariance checks.

### 2.8 `pipeline/run_mock_validation.py` (Stage A acceptance test; also the first result of the project)
1. Build xi table from the mock's own P_F; build pairs + catalogue for one mock; print the benchmark.
2. Amplitude against the **true** `alpha_lya` template: `A` must be 1 within the jackknife error and within 10% for
   the noiseless mock (`P_N = 0`, high `n_los` option). This is the normalisation test.
3. Injection test (2.6) on the mock: slope 1 within 5% (noiseless) and consistent with 1 (noisy).
4. Curl-template amplitude consistent with zero.
5. Deprojection: amplitude against the filtered `kappa_CMB` template, then against `kappa_CMB - kappa_s_hat` with the
   matched template built from the mock quasars. The first must be biased high by the predicted response term
   (`(R_delta/2) R_ka C^{dkc}/C^{klkc}`, of order 1.6 for the fiducial numbers), the second unbiased. Also the
   amplitude against `kappa_s_hat` alone (response control) must be consistent with the injected `R_delta/2`.
6. Null ensemble: 100 Gaussian `kappa` realisations with the `kappa_CMB` spectrum as templates; the std of `A` must
   agree with the jackknife error and with `1/sqrt(F)` to within the expected factors (report and discuss).
7. Covariance: repeat 2-5 over >= 10 seeds; compare scatter with the per-mock errors.
Write all numbers to `report/mock_validation.json` and a short markdown summary; make figures into `report/figures/`.

---------------------------------------------------------------------------------------------------------------------
## 3. Modules (Stage B: real data)

### 3.1 `pipeline/fetch.py`
Download with resume and checksums into `/data/LyaLenser/raw/`:
- DESI DR1 Lya deltas: https://data.desi.lbl.gov/public/dr1/vac/dr1/lya-deltas/ (directory `delta-lya-0-0/Delta/`,
  per-HEALPix `Delta-*.fits.gz`; also `delta_attributes`). Inspect the FITS structure (picca "ImageHDU" format:
  METADATA table with RA, DEC, Z, TARGETID plus LAMBDA/DELTA/WEIGHT/CONT arrays) before writing the reader.
- DESI DR1 QSO catalogue and LSS randoms (public DR1 LSS release; locate the exact files and record URLs in
  MEMORY.md).
- ACT DR6 lensing maps, masks, N_L and simulations from LAMBDA (baseline and tSZ-deprojected); follow the LAMBDA
  instructions for cross-correlation use.
- Planck PR4 lensing (Carron, Mirmelstein & Lewis 2022): klm, mask, N_L, sims.
Record every file's URL, size, md5 in `/data/LyaLenser/raw/MANIFEST.json`.

### 3.2 `pipeline/desi_io.py`
`read_deltas(paths, zmin, zmax, cfg) -> SightlineSet`: wavelength -> z -> chi (`cosmo.chi`), keep pixels with
`zmin <= z <= zmax`, apply DLA masks if a DLA catalogue is given (mask +-`dv` around each DLA), drop forests with
fewer than 50 pixels, store per-pixel weights. `quasar_catalogue(...)`, `randoms(...)` -> arrays (ra, dec, z).
Also write the per-pixel noise distribution to `report/noise_distribution.json` (histogram of `P_N,a`), for the
forecast update.

### 3.3 `pipeline/cmb_io.py`
Readers for ACT DR6 and Planck PR4 alm + masks + N_L; simulation iterators; a function producing the filtered
`kappa_alm` given a filter and applying the recommended mask treatment.

### 3.4 `pipeline/run_dr1.py`
Single-slab first run: read deltas (2.1<z<3.0) -> pairs -> catalogue (save); xi table from the DESI DR1 fitted
model (implement the same Kaiser+distortion-free model first; the distortion is absorbed by the injection
calibration); templates: ACT signal, Planck signal, matched `kappa_s_hat` from DR1 quasars in the slab with a
fiducial `b_q(z) = 0.278((1+z)^2 - 6.565) + 2.393`, curl, 400 ACT sims; amplitudes, band amplitudes for bands
[40,100], [100,200], [200,300]; injection calibration with 5 injections; jackknife at nside 8. Outputs to
`report/dr1_run.json` and figures.

---------------------------------------------------------------------------------------------------------------------
## 4. Testing and acceptance
- `pytest code/pipeline/tests` must pass: kernel vs brute force, thread determinism, xi table checks, deflection
  sign/derivative checks, amplitude on a synthetic catalogue with a known injected `d`.
- Stage A acceptance = Sec. 2.8 items 2-6 within the stated tolerances, with the mock validation summary committed.
- Every script runnable standalone from `code/pipeline/` and writing to `report/figures/` or `/data/LyaLenser/`.
- Log wall times; the pair kernel benchmark is reported in the summary.

## 5. Out of scope for the first implementation
Full C^-1 weighting; tomographic sub-slabs (keep the `slab` field so it can be added); DESI mocks at NERSC;
magnification slope measurement; released ACT/Planck N_L in the Wiener filter (band-pass filters are enough first).
