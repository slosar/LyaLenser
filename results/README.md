# Committed results

JSON summaries of the pipeline products (`docs/pipeline.md`). Every number in the paper and in `report/lowz.tex`
is read from one of these files by the figure scripts or quoted from it. Tags: `v4` / `v7d` / `v1d` mark the
fiducial analysis (templates `lowz_v4`, forest run `dr1_lowz_v7d`, quasar run `dr1_qso_v1d`).

## Fiducial

| File | Content |
|---|---|
| `joint_fit_v4.json` | **The headline numbers.** Joint 55-component response fit of both statistics: `statistics/{auto,cross}` (component names, global, per-slice and per-band amplitudes with jackknife errors) and `combination/{global, per_slice, per_band, block_diagonal_global, no_curl_*}` (A = 0.61 +- 0.33 combined; 0.42 +- 0.43 auto; 0.86 +- 0.47 cross). `joint_response_fit.py`. |
| `dr1_lowz_v7d.json`, `.md` | The forest auto-correlation run: sample counts, the fitted correlation (`xi_fit`, `xi_fit_uncorrected`), per-slice and combined amplitude fits, curl nulls, injection expectation, random-template fits, wall time and memory. `run_dr1_lowz.py`. |
| `dr1_qso_v1d.json`, `.md` | The quasar-forest run, same structure plus the quasar sample and the sub-slab fits. `run_dr1_qso_lowz.py`. |
| `auto_subslabs_v7d.json` | Auto amplitude in three ranges of the pair mean redshift. `auto_subslabs.py`. |
| `auto_cross_combination_v4.json` | Combination of the per-slice fits: overall, per band, per slice and per sub-slab (redshift-split figure). `combine_auto_cross.py`. |
| `deflection_cmb_check_v4.json` | Templates against ACT DR6 and Planck PR4: per band and over the science window, convergence (spin-0) and deflection (spin-1 E and B), predictions, amplitudes A_L and chi-square. `deflection_cmb_check.py`. |
| `injection_v4.json`, `injection_v4_nosw.json` | Injection tests on the fiducial products (expectation and data, both statistics); `_nosw` with the same-wavelength term removed from the injected correlation (slope 0.991). `injection_dr1.py`. |
| `scale_sensitivity_v4.json` | Fisher-information density on the (r_perp, r_par) plane for both statistics. `scale_sensitivity.py`. |
| `dr1_tracer_biases.json`, `cmb_bias_check.json`, `cmb_bias_check_planck.json` | Tracer biases from the auto-spectra and from the cross-spectra with the CMB maps (bias figure). |

## Robustness rows (`paper_figures.py --robustness`)

`auto_cross_combination_{l500v4,l1300v4,rp20,rp40,kmed,kfine}.json`: bands to L = 500 and 1300, r_perp cut at
20 and 40 Mpc/h, spline knot sets medium (30 + 42 coefficients) and fine (42 + 63), each from a pair of
`dr1_lowz_<tag>` / `dr1_qso_<tag>` products of `scripts/slurm/robustness.sbatch`.

## Superseded versions (kept for the record and for the report)

`dr1_lowz{,_v3,_v4,_v5_lya_only,_v6,_v6_flat,_v7,_v7b,_v7c}`, `dr1_qso_v1{,b,c,_free}`, `auto_subslabs_v7{,b,c}`,
`auto_cross_combination{,_wide,_l1300}`, `deflection_cmb_check{,_v2,_wide,_l1300}`, `cmb_map_checks*`,
`cmb_bias_check*_v2`, `joint_fit_v3b`, `lowz_v{2,3,4}_summary`, `injection_same_wavelength_v3`,
`slice_crosstalk`, `xi_zevol_dr1`, `dry_run_disc190_30_12`: the intermediate iterations described in
`docs/history.md` (tags: v3 single region, v4 both regions, v5 Lya only, v6 NaMaster/BGS/BOSS, v7 whole
forest below z = 3, v7b bands to 1000, v7c lmax 1300, v7d measured-covariance Wiener; the quasar runs b/c/d
follow the same sequence).
