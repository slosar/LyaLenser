# The pipeline

Everything runs from `scripts/` with the conda environment of `docs/computing.md` and `LYALENSER_DATA`
pointing at the data tree of `docs/data.md`. `scripts/run_fiducial_chain.sh` chains stages 2 to 7 with the
fiducial options; the stage-1 basis tables are built once. Product directories are named with a version tag
(`lowz_v4`, `dr1_lowz_v7d`, `dr1_qso_v1d` are the fiducial ones); the summaries copied to `results/` carry
the same tags.

## Fiducial choices

| Choice | Value | Where it is set |
|---|---|---|
| Forest sample | every DR1 pixel with 1.96 <= z <= 3.0 in both delta regions (Lya and Lyb windows), picca weights | `run_dr1_lowz.py --zmin 1.96 --zmax 3.0 --regions lya lyb` |
| Source plane | z_eff = 2.3476, the weight-averaged pixel redshift; the templates are built for this plane | `lowz_catalogues.py --zref`, `run_dr1_lowz.py --zeff` |
| Pair selection | 3 <= r_perp <= 30 Mpc/h, abs(r_par) < 30 Mpc/h | `lyalenser.config.production_config`, `--rperp-max` |
| Correlation model | Kaiser basis through the DESI pixel window and the per-forest continuum projection, redshift-evolving, plus a bicubic spline correction (16 + 6 coefficients) and a same-wavelength term | `--xi-correction spline --xi-knots bicubic`, `lyalenser.tables` |
| Tracers | DESI DR1 LSS v1.5 BGS, LRG, ELG, QSO and BOSS DR12 CMASSLOWZTOT in five slices: 0.1-0.4 (BGS, BOSS), 0.4-0.6 (LRG, BOSS), 0.6-0.8 (LRG, BOSS), 0.8-1.1 (LRG, ELG, QSO), 1.1-1.6 (ELG, QSO) | `lyalenser.lowz.TRACERS`, `lowz_catalogues.py --tracer-zmax 1.6` |
| Template combination | per-slice Wiener weights from the measured NaMaster auto/cross spectra of the tracer maps, per coverage class | `lowz_catalogues.py --wiener measured` |
| Science bands | 40-200, 200-400, 400-600, 600-800, 800-1000 (cosine-tapered, 10-multipole ramps); each with a curl partner and one junk component per slice: 11 components per slice, 55 in the joint fit | `--bands 40 200 400 600 800 1000`, `lyalenser.templates` |
| Deflection evaluation | alpha = grad phi from the alm at nside 2048 | `--nside-alpha 2048` |
| Errors | jackknife over HEALPix nside-8 regions of the pair midpoints (302 regions) | `--nside-jk 8` |

## Stages

### 1. Basis tables (once per forest sample)
`build_basis_dr1.py --out $D/stageb/basis_dr1_ab_z196.h5 --regions lya lyb --zmin 1.96` (15 min) computes the
three Kaiser basis correlations P_lin(k) mu^{2i} through the DESI line-of-sight pixel and resolution windows
(`lyalenser.xi_fit.basis_tables`) and pushes them through the continuum projection averaged over real DR1
forest pairs. `build_basis_cross_dr1.py --basis ... --out $D/stageb/basis_cross_dr1_z196.h5` does the same
with the one-sided projection for quasar-forest pairs (`lyalenser.xi_cross`).

### 2. Templates from the tracers
`lowz_catalogues.py --out $D/lowz_v4 --zref 2.3476 --tracer-zmax 1.6 --wiener measured` (about 10 min).
For every (tracer, slice): the kernel-weighted HEALPix map at nside 512 (`lyalenser.templates.matched_template`),
its footprint and completeness from the randoms, the linear bias from the NaMaster cross-spectrum of two
random halves against the Limber prediction on 40 <= ell <= 0.2 chi (`lyalenser.lowz`, `lyalenser.nmt_spectra`),
the unit-bias map and its shot noise. The slice maps are combined with per-multipole Wiener weights inside each
coverage class (the set of tracers covering a pixel) and summed into the combined convergence template.
Products: `kappa_slice_*_alm.fits`, `kappa_combined_alm.fits`, `mask_*.fits`, `unitbias_*.fits`, `summary.json`
(biases, spectra, weights, slices).

### 3. Template validation against CMB lensing
`deflection_cmb_check.py --lowz $D/lowz_v4 --bands 40 200 400 600 800 1000 --out results/deflection_cmb_check_v4.json`
(about 10 min). Builds exactly the deflection field the estimator consumes on the sphere, cross-correlates its
E-mode with the ACT DR6 and Planck PR4 convergence maps with NaMaster (spin-1 x spin-0), and fits the amplitude
A_L of the measurement relative to the prediction for the Wiener-filtered template
(`lyalenser.template_prediction`, class fractions inside each mask overlap). Covariance from the Gaussian
estimate with model spectra. `cmb_bias_check.py` and `cmb_map_checks.py` are the older per-tracer bias and
map-consistency checks; `plot_template_cmb_cross.py` draws the paper figure from the JSON.

### 4. Forest auto-correlation measurement
`run_dr1_lowz.py --out $D/stageb/dr1_lowz_v7d --lowz $D/lowz_v4 --basis $D/stageb/basis_dr1_ab_z196.h5 --regions lya lyb --xi-correction spline --randoms 20 --nside-jk 8 --zmin 1.96 --zmax 3.0 --zeff 2.3476 --bands 40 200 400 600 800 1000 --nside-alpha 2048`
(75 min, 34 GB with `NUMBA_NUM_THREADS=22`). Steps:
1. `lyalenser.desi_io.read_deltas`: all forests of both regions into one `SightlineSet` (428k forests, 2.8e8 pixels), saved as `sightlines.h5`.
2. `lyalenser.xi_model.xi_from_data`: 1 Mpc/h cells of the pair counts in redshift bins; `lyalenser.xi_zevol.fit_evolving_table` fits the redshift-evolving model with the spline correction (`xi.h5`, parameters in the JSON under `xi_fit`).
3. `lyalenser.pairs.find_pairs` / `accumulate`: the pair catalogue (6.5e6 sightline pairs) with the accumulators of the estimator, G = chi dxi/dr_perp from the layered table at each pair's mean distance (`catalogue.h5`).
4. `lyalenser.templates.sphere_band_templates`: the band deflections of every slice and of the combined template at the sightlines.
5. `lyalenser.amplitude.amplitude` per slice and for the combined template with the jackknife (`fits.h5`), the curl null, the injection expectation with the template's own deflection, and `--randoms` fits against Gaussian random templates.
Summary: `dr1_lowz.json` and `dr1_lowz.md`, copied to `results/dr1_lowz_<tag>.{json,md}`.

### 5. Quasar-forest cross-correlation measurement
`run_dr1_qso_lowz.py --out $D/stageb/dr1_qso_v1d --auto-run $D/stageb/dr1_lowz_v7d --lowz $D/lowz_v4 --basis $D/stageb/basis_cross_dr1_z196.h5 --randoms 20 --sub-slabs 1.96 2.25 2.55 3.0 --spline-fixed-base --bands ... --nside-alpha 2048`
(16 min, 26 GB). The DR1 quasar catalogue (`lyalenser.qso_io`, one row per TARGETID, 1.9 < z < 3.1) against the
sightline set of the auto run: the signed (r_perp, r_par) cells, the evolving cross model with the forest side
fixed from the auto fit and the quasar bias, redshift offset and smoothing free (`lyalenser.xi_cross`), the
quasar-sightline pair catalogue (9.5e6 pairs) and the same amplitude fits. Summary `dr1_qso.json`.

### 6. Splits and combination
- `auto_subslabs.py --run $D/stageb/dr1_lowz_v7d --lowz $D/lowz_v4 --bands ...`: the auto amplitude in three ranges of the pair mean redshift (`auto_subslabs.json`).
- `joint_response_fit.py --auto ... --cross ... --lowz ... --bands ... --tag v4`: **the fiducial fit**. One 55 x 55 response matrix per statistic over all slice components (`lyalenser.joint_fit`), the collapsed amplitude with curl and junk marginalised, per-slice and per-band amplitudes, the block-diagonal and no-curl variants, and the auto x cross combination with the joint jackknife covariance. Writes `results/joint_fit_<tag>.json` and the response-matrix figure.
- `combine_auto_cross.py --auto ... --cross ... --out results/auto_cross_combination_<tag>.json`: the combination of the per-slice fits, overall, per band and per sub-slab (used for the redshift split and the robustness rows).

### 7. Validation and diagnostics on the products
- `injection_dr1.py --auto ... --cross ... --lowz ... --bands ... --tag v4_corrected`: shifts every position by -A alpha for A = +-0.25, +-0.5, preserves A/B labels, rebuilds the production-selected pairs and refits the combined template. Both the noise-free expectation and the data injection have paired jackknife slope errors. `--remove-same-wavelength` checks the bookkeeping of the same-wavelength term. These are 11-component diagnostics, not a physical calibration of the 55-component slice fit. See `docs/audit_20260925.md` for the correction to earlier injection runs.
- `random_template_null.py`: reuses cached pairs to run 100 independent Gaussian combined-template nulls per statistic; `summarise_audit.py` compares the original 20 and additional 80, ensemble scatter, and jackknife errors (`results/random_template_null_100_v4.json`).
- `scale_sensitivity.py --tag v4_corrected`: the unmarginalised response density on the (r_perp, r_par) plane, including the production source-distance terms, with an independent total-response check (`results/scale_sensitivity_v4_corrected.json`, `--plot-only` redraws).
- `source_distance_expansion.py --run $D/stageb/dr1_lowz_v7d --summary results/dr1_lowz_v7d.json`: the size of the linear source-distance expansion of the lensing efficiency (common g1 and effective lens distance, per-slice values, weighted pixel-distance moments, fractional deflection change and the estimated amplitude bias of using one g1 for all slices; `results/source_distance_expansion.json`, quoted in the paper's estimator section).
- `slice_crosstalk.py`, `xi_zevol_dr1.py`, `qso_xi_qa.py`, `injection_same_wavelength.py`, `dry_run_lowz.py`: older diagnostics kept for reference (slice-to-slice response leakage, the redshift evolution of the forest correlation, QA plots of the cross fit, the same-wavelength term, a one-disc dry run).

### 8. Robustness rows
`scripts/slurm/robustness.sbatch` (NERSC; `TAG`, `LOWZ`, `BANDS`, `EXTRA`, `STEPS` variables) reruns stages 4 and 5
with one choice changed: bands to 500 or 1300 (`lowz_catalogues.py --lmax 1300` first), `--rperp-max 20` or `40`,
`--xi-knots medium` or `fine`. Each pair of products is combined with `combine_auto_cross.py` into
`results/auto_cross_combination_<tag>.json`, and `paper_figures.py --robustness` draws the rows.

## Figures of the paper and the report

| Script | Figures |
|---|---|
| `report_figures.py` | `kernel_slices.pdf`, `tracer_spectra.pdf`, `biases.pdf`, `dr1_xi.pdf`, `xi_correction.pdf`, and the report-only figures |
| `plot_templates.py --combined-only` | `templates_map_combined.pdf` |
| `plot_template_cmb_cross.py --tag v4` | `template_cmb_cross.pdf` |
| `paper_xi_figures.py` | `xi_ff_fit.pdf`, `xi_qf_fit.pdf` |
| `scale_sensitivity.py --plot-only` | `scale_sensitivity.pdf` |
| `joint_response_fit.py --plot-only v4 auto` | `response_matrix_auto_v4.pdf` (the paper's `response_matrix.pdf`) |
| `paper_figures.py --split / --slices / --bands / --robustness` | `redshift_split.pdf`, `slice_split.pdf`, `band_split.pdf`, `robustness.pdf` |

Every figure environment in `Paper/main.tex` and `report/lowz.tex` carries a `%% To reproduce:` comment with
the exact command. The scripts write to `report/figures/` and, when the directory exists, to `Paper/figures/`.

## Package map

| Module | Contents |
|---|---|
| `paths`, `cosmo`, `lensing` | data locations (`LYALENSER_DATA`), the CAMB background and power spectra, the lensing kernel and Limber integrals |
| `config` | `Config` (all numerical choices, `production_config()` for the DR1 values) and `SightlineSet` |
| `desi_io`, `qso_io`, `cmb_maps` | DR1 deltas to `SightlineSet` (both regions, `load_sightlines`), the DR1 quasar catalogue, the ACT/Planck convergence maps and masks |
| `forest_power`, `xi_model`, `xi_fit`, `xi_spline`, `xi_zevol`, `xi_cross`, `tables` | the forest power spectrum model, the correlation table `XiTable` (layered in chi), the pair counts, the Kaiser basis with pixel windows and continuum projection, the spline correction, the redshift-evolving fits for forest-forest and quasar-forest, and the production wrappers `table_for` / `correction_for` / `read_xi` |
| `pairs` | numba pair finding and accumulation (`PairCatalogue`), jackknife regions |
| `templates` | phi from kappa, `alpha_at`, cosine bands, `sphere_band_templates` (bands, curls, junk), the kernel-weighted tracer map |
| `amplitude`, `joint_fit`, `inject` | the amplitude fit with partial sums per jackknife region (`AmplitudeResult`, `common_science`), the joint fit over slices, the injection test |
| `lowz`, `nmt_spectra`, `template_prediction` | tracer slices and Limber spectra, NaMaster bandpowers and covariances, the prediction of the template x CMB cross-spectrum |

## Things that bit us (keep in mind)

- The DR1 delta VAC has two forest regions chunked independently (`delta-lya-0-0`, `delta-lyb-0-0`); they match on LOS_ID only, never on file number.
- The DR1 per-cell scatter of the correlation is about half of 1/sqrt(sum w_p w_q): the picca weights are fine for relative weighting, not for absolute chi-square.
- NaMaster on a fragmented binary mask overshoots white (non band-limited) levels by 13-16 per cent; shot-noise levels are measured as pseudo-C_ell / mean(mask^2), signal spectra with NaMaster. Fields must be built with the bins' lmax and the same `n_iter` as the alm.
- The ACT mask goes squared on the kappa side for cross-spectra (release README); the ACT map is masked on input.
- The decoupled EE auto-spectrum of the science-window deflection is unusable above L about 600 (low-L power leaks through the fragmented mask), so the validation covariance uses model spectra.
- The band windows are disjoint cosine ramps, so the junk component absorbs the notches between bands (no bias, a few per cent of S/N).
- With the larger spline knot sets the base parameters (beta_F, gamma_beta) go degenerate with the correction while the kernel and the amplitude stay put; quote the amplitude, not the base parameters.
- `pkill -f <pattern>` from a shell whose command line contains the pattern kills that shell.
