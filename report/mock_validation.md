# Stage A iteration 3 validation

**20/30 gates PASS. Stage B remains blocked.**

Recovery and null amplitudes are absolute `F^-1(q-mf)`. Each mean-field subtraction belongs to the sample being fitted. Paired changes are labelled diagnostics.

The scalar common-science estimate constrains the three science-band amplitudes to one value while fitting the three curls and junk as separate nuisance components. The full seven-component response and amplitude vector are also saved.

The reported bound is `abs(mean − target) + t(0.975, N−1) × SEM`, with sample SEM. Amplitude bounds have A units; slope bounds are dimensionless. Normalization decisions use the predeclared ±0.03 mean-slope tolerance, with uncertainty also reported. A two-SEM consistency alone does not establish precise normalization.

Recovery seeds (40): [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59].

Null seeds (40): [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59].

| Gate | Measurement | Frozen tolerance | Result |
|---|---|---|---|
| truth normalization slope | 0.999183 ± 0.0273732 SEM; 95% residual bound 0.0561848 (slope units); N=40 | &#124;slope-1&#124; <= 0.03 | PASS |
| cmb normalization slope | 1.04582 ± 0.066089 SEM; 95% residual bound 0.179499 (slope units); N=40 | &#124;slope-1&#124; <= 0.03 | FAIL |
| matched normalization slope | 1.4748 ± 0.862855 SEM; 95% residual bound 2.22009 (slope units); N=40 | &#124;slope-1&#124; <= 0.03 | FAIL |
| fixed versus refitted baseline | Fixed slope 1.0118 ± 0.0539595 SEM; 95% residual bound 0.120947 (slope units); N=40; refitted slope 0.999183 ± 0.0273732 SEM; 95% residual bound 0.0561848 (slope units); N=40; diagnostic refitted − fixed slope -0.0126212 ± 0.055523 SEM; 95% residual bound 0.124927 (slope units); N=40; diagnostic A_true=1 difference 0.153295 ± 0.14608 SEM; 95% residual bound 0.44877 A; N=40 | &#124;slope difference&#124; <= 0.03 | PASS |
| varying template mean field: truth | -0.476432 ± 0.463564 SEM; 95% residual bound 1.41408 A; N=40; raw -4.00673 ± 2.50156 SEM; 95% residual bound 9.06662 A; N=40 | &#124;mean&#124; <= 2 SEM | PASS |
| fixed template mean field: truth | 0.226692 ± 0.455374 SEM; 95% residual bound 1.14777 A; N=40; raw -2.19538 ± 2.65932 SEM; 95% residual bound 7.57436 A; N=40 | &#124;mean&#124; <= 2 SEM | PASS |
| varying template mean field: cmb | -0.264206 ± 1.33898 SEM; 95% residual bound 2.97254 A; N=40; raw -1.39784 ± 8.79584 SEM; 95% residual bound 19.1891 A; N=40 | &#124;mean&#124; <= 2 SEM | PASS |
| fixed template mean field: cmb | -0.0292683 ± 0.89486 SEM; 95% residual bound 1.83929 A; N=40; raw -1.13027 ± 6.09536 SEM; 95% residual bound 13.4593 A; N=40 | &#124;mean&#124; <= 2 SEM | PASS |
| varying template mean field: matched | 25.612 ± 8.87865 SEM; 95% residual bound 43.5707 A; N=40; raw 480.184 ± 30.0622 SEM; 95% residual bound 540.991 A; N=40 | &#124;mean&#124; <= 2 SEM | FAIL |
| fixed template mean field: matched | -8.65689 ± 9.65828 SEM; 95% residual bound 28.1926 A; N=40; raw -22.5998 ± 34.2376 SEM; 95% residual bound 91.8519 A; N=40 | &#124;mean&#124; <= 2 SEM | PASS |
| absolute stochastic recovery | -0.437805 ± 1.25867 SEM; 95% residual bound 3.9837 A; N=40 | &#124;mean-1&#124; <= 2 SEM (precision established separately by slope) | PASS |
| absolute response-only: cmb | 1.36026 ± 1.34815 SEM; 95% residual bound 4.08715 A; N=40; predicted 1.18778 ± 0.417906 SEM; 95% residual bound 2.03308 A; N=40; observed − predicted 0.172482 ± 1.07035 SEM; 95% residual bound 2.33747 A; N=40 | &#124;observed-predicted&#124; <= 2 SEM | PASS |
| absolute response-only: matched | 52.1544 ± 9.43084 SEM; 95% residual bound 71.23 A; N=40; predicted 33.1558 ± 2.09324 SEM; 95% residual bound 37.3898 A; N=40; observed − predicted 18.9985 ± 8.76079 SEM; 95% residual bound 36.7189 A; N=40 | &#124;observed-predicted&#124; <= 2 SEM | FAIL |
| absolute response-only: deprojected | -2.96601 ± 0.95635 SEM; 95% residual bound 4.90041 A; N=40; predicted -1.41655 ± 0.266485 SEM; 95% residual bound 1.95557 A; N=40; observed − predicted -1.54946 ± 0.826059 SEM; 95% residual bound 3.22032 A; N=40 | &#124;observed-predicted&#124; <= 2 SEM; absolute null bound95 <= 0.3 A | FAIL |
| absolute combined deprojected recovery | -1.63392 ± 0.71317 SEM; 95% residual bound 4.07644 A; N=40 | &#124;mean-1&#124; <= 2 SEM and residual bound95 <= 0.3 A | FAIL |
| six-bin Hotelling shape | Hotelling T²=8.22936; F=1.47706; df=[5, 35]; p=0.222214 | p > 0.01 | PASS |
| absolute covariance | {"scatter": 4.51048116872763, "rms_jk": 6.4812099853175456, "ratio": 0.6959319600731374} | 0.7 <= scatter/RMS jackknife <= 1.3 | FAIL |
| outside-box and total spectra | klkl: 0.987398 ± 0.00836918 SEM; klkc: 0.988494 ± 0.00956604 SEM; kckc: 1.0187 ± 0.00612389 SEM | each total ratio within 10% | PASS |
| analytic derivative convergence | {"nk": 6400, "relative_norm_change": 2.7652501512515034e-05} | doubling nk changes derivative norm by <1% | PASS |
| discrete-grid covariance quadrature | {"nangle": 128, "doubled": 256, "derivative_norm_change": 7.871991201188842e-07} | doubling angular samples changes derivative norm by <1% | PASS |
| unit and regression tests | 36 passed, 0 failed; 50.9036 s | all tests pass | PASS |
| xi interpolation convergence | {"A_025": 0.877295842008689, "A_0125": 0.8772376302938488, "fractional_change": 6.635797739401357e-05} | relative absolute A change < 0.005 | PASS |
| dense physical g=True | Slope 1.07026; descriptive regression error 0.0231019; ensemble SEM and 95% bound unavailable (one realization); absolute A=[0.0746338365052487, 0.877295842008689, 5.114264713521431, 10.695393731177832] | slope within 0.05 of 1 | FAIL |
| data versus analytic baseline | {"data_slope": 1.0702613585248606, "analytic_slope": 0.8966547732000394, "ratio": 1.1936158603218525, "analytic_A": [-0.5384569157335565, 0.06261713307471273, 3.9889586115249345, 8.272943510527902]} | slope ratio within 3% | FAIL |
| first moments | {"observed_ratio": 1.0208597374579793, "predicted_ratio": 1.0181112008246955} | omitted/correct slope ratio within 5% of catalogue prediction | PASS |
| dense physical g=False | Slope 1.07158; descriptive regression error 0.0229881; ensemble SEM and 95% bound unavailable (one realization); absolute A=[0.021643859068786128, 0.8134072563817436, 5.075177674028912, 10.648222912984219] | slope within 0.05 of 1 | FAIL |
| injection bookkeeping | Science odd slope 0.991957; curl slope -0.0934685; frozen RMS individual-band JK error scale 3.02845; joint slope-error audit below | odd slope within 0.05 of 1; curl within jackknife error | PASS |
| flag and common-realization margin diagnostics | magnification off: absolute A=-17.2997, diagnostic shift=-11.4388; completeness off: absolute A=-17.6984, diagnostic shift=-11.8375; real_mask off: absolute A=-1.18693, diagnostic shift=4.67395; nested margins persisted | all requested diagnostics persisted | PASS |
| 100 random templates | -0.599917 ± 0.391589 SEM; 95% residual bound 1.37691 A; N=100; scatter/RMS JK=0.957035; scatter/RMS σF=5.12854 | finite absolute statistics and ratios (report only) | PASS |
| benchmark and memory | 8.12452e+08 pixel pairs/s; 24 threads; peak 25.8568 GiB (27.7636 GB) | 24 Numba threads; positive throughput; peak <40 GB | PASS |

## Prospective stopping decision

After the initial 40 seeds, recovery: -2.40071 ± 0.963322 SEM; 95% residual bound 5.41697 A; N=20; response-only null: -4.56773 ± 1.53847 SEM; 95% residual bound 7.78778 A; N=20.

Extension triggered: True. Added seeds: [40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59]. The decision was made once; no seeds beyond 59 were used.

## Absolute amplitude points and shape

All per-seed amplitude-grid points, raw/corrected fits, jackknife errors, predicted responses, and shape estimates are retained in [the JSON report](mock_validation.json). Figures show ensemble SEM, except individual combined-recovery points, which show midpoint-jackknife errors.

- [Physical normalization](figures/mock_physical_normalisation.pdf)
- [Absolute combined recovery](figures/mock_combined_ensemble.pdf)
- [Six-bin shape](figures/mock_shape.pdf)
- [Spectra](figures/mock_spectra.pdf)

Shape bins 0–5 are, in order, transverse separations 0–10, 10–20, and 20–30 Mpc/h, each split into radial separations 0–10 and 10–30 Mpc/h.

## Same-realization margin diagnostic

All margins below use seed 0 with the same larger density realization, forest, sightlines, and CMB map. Differences from the 150 Mpc/h selection are diagnostics. Flag and margin checks each use one realization, so an ensemble SEM or 95% bound is unavailable.

| Margin (Mpc/h) | Absolute A | Diagnostic change from 150 |
|---|---|---|
| 0 | -5.32034 | 0.540537 |
| 150 | -5.86088 | 0 |
| 300 | -9.3153 | -3.45443 |

## Coordinate injection bookkeeping

The values below are absolute fits on one realization. On this symmetric injection grid, the reported odd slope equals the least-squares slope with a free intercept. This tests coordinate bookkeeping; physical normalization is assessed by the fresh remapping ensemble. An ensemble SEM or 95% bound is unavailable for this single-realization control.

| Injected A | Absolute fitted A | Absolute mean curl amplitude |
|---|---|---|
| -2 | -1.65178 | -1.10318 |
| -1 | -0.74718 | -1.427 |
| 0 | 0.38581 | -1.51676 |
| 1 | 1.31042 | -1.54363 |
| 2 | 2.27921 | -1.5122 |

Joint curl-slope jackknife audit: slope -0.0934685, slope error 0.167859; within one sigma: PASS. The frozen RMS individual-band error scale is 3.02845. The audit uses the common union of saved midpoint-HEALPix regions and applies the original jackknife variance formula to the leave-region slopes. It is reported separately from the frozen acceptance decision.

## Wall time and memory

```json
{
  "ensemble_s": 14711.106358294375,
  "extras_s": 5040.3598235887475,
  "this_invocation_wall_s": 19823.410996453837,
  "initial_development_wall_s": 297.108978160657,
  "expanded_development_scan_wall_s": 2423.861297130119,
  "development_note": "initial and expanded development jobs overlapped; their wall times are not additive"
}
```

Successful campaign peak RSS: 25.8568 GiB (27.7636 GB).

The initial development process was interrupted with SIGKILL; its last observed RSS is not its known termination peak. The allocation repair preceded validation. Overlapping development-job durations must not be summed as elapsed campaign time.

## Provenance and reconstruction

Numerical choices were frozen on development seeds 100–104. [GATES.md](../GATES.md) predates validation; its hash, source hashes, configuration hash and library versions are recorded below. The reporting formatter does not recompute fits or change PASS/FAIL decisions.

```json
{
  "UTC_start": "2026-09-11T03:37:16.619495+00:00",
  "recovery_seeds": [
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19
  ],
  "null_seeds": [
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    32,
    33,
    34,
    35,
    36,
    37,
    38,
    39
  ],
  "extension_rule": "all seeds 40-59 once iff either required residual bound95 >0.3",
  "GATES_sha256": "d3849680a1f088083beac26f9ab6fc21390322413d1d5a20e2ed89a86375c7a8",
  "source_sha256": {
    "__init__.py": "da209807a22513917f33dcb62bb1e2841fe131ed679e8d96ab1e5e5c16012e09",
    "amplitude.py": "c3853e662ae37f3f4972b7c44c715883f1a72aefe0b3f45d792155ae4be155bf",
    "config.py": "9a5c5c2a71f174cbf26ac2381a522104e49c6ec51f81b09a88de90b44c8edddb",
    "grid_covariance.py": "8b653ff312f06a22b66555216cbf6aaf72194d88b273456e0b1ed8dd309b32a4",
    "inject.py": "f78282ed6656ee824c4413351e6cab9584abd0a7e44ae6332639b822e266b0c4",
    "mock.py": "7a5aa9d15abd7dcdc42aea85a1751ffdb98a0b826e511041bda84f1079ef7e9e",
    "pairs.py": "f3bc8186d39990f887b3c0ff3754931e31d75d45132abea113f5162bdda27c03",
    "response.py": "41f6bf47bbadd97ff913a882e695802f2c8e92cda6780140a25c56ba0f07c025",
    "run_mock_validation.py": "eb89a7c69ba6e88683005c71b659358d28d2ef643c442158547bf7d56aaf7f64",
    "templates.py": "06336e8f94bc6b23bcce39e92db7765c9a7169d3088c70c4a102a5384d57413f",
    "validation_stats.py": "df609e53f00c89d74a1d2ef01a5f582f8a1bc0d76ccc0411264010d38ac2595d",
    "xi_model.py": "e921ef9d97956f4899caf94e222aa2d084368579f02c2dc6a52dff7adbe51565"
  },
  "versions": {
    "numpy": "1.26.4",
    "scipy": "1.16.1",
    "numba": "0.59.0",
    "healpy": "1.18.1",
    "h5py": "3.9.0",
    "camb": "2.0.4",
    "fitsio": "1.2.1"
  },
  "python": "3.11.7 (main, Dec 15 2023, 18:12:31) [GCC 11.2.0]",
  "config_sha256": "a22cf484f5831c75bb5bf3ce79c7b3c3827fa842fa00462ce5a82809c8bf7e37",
  "development_seeds": [
    100,
    101,
    102,
    103,
    104
  ],
  "xi_smoothing": 2.0,
  "xi_step": 0.25,
  "analytic_nk": 6400,
  "threads": 24
}
```

From `code/pipeline/`, run `python run_mock_validation.py --phase rebuild`; then, from the repository root, run `python report/audit_injection_jackknife.py` and `python report/render_iteration3.py`. This uses the persisted ensemble, controls and diagnostics without generating new mocks.

## Remaining failures

- **cmb normalization slope**: Slope of absolute amplitudes at A_true={0,1,5,10}; per-seed intercept fitted freely.
- **matched normalization slope**: Slope of absolute amplitudes at A_true={0,1,5,10}; per-seed intercept fitted freely.
- **varying template mean field: matched**: The frozen tolerance was not met.
- **absolute response-only: matched**: Prediction uses stored delta_L at every pixel and P_a (D_a C_ab D_b - C_ab) P_b^T. Discrete-grid, phase-averaged trilinear covariance truncated at saved table support.
- **absolute response-only: deprojected**: Prediction uses stored delta_L at every pixel and P_a (D_a C_ab D_b - C_ab) P_b^T. Discrete-grid, phase-averaged trilinear covariance truncated at saved table support.
- **absolute combined deprojected recovery**: The frozen tolerance was not met.
- **absolute covariance**: The frozen tolerance was not met.
- **dense physical g=True**: One full scale-1 noiseless dense realization. Error is regression error, not an ensemble SEM.
- **data versus analytic baseline**: The empirical 0.903 factor is removed. A continuum- and grid-mismatched analytic baseline can fail this gate.
- **dense physical g=False**: One full scale-1 noiseless dense realization. Error is regression error, not an ensemble SEM.

### Additional unresolved independence issue

A source audit after the freeze found that the field and CMB-noise generators initialize separate random generators with the same seed. Their first normal draws reuse the same stream prefix. The saved small-case reproduction confirms exact prefix reuse. The effect after Fourier filtering and array reshaping has not been isolated; it is not claimed to explain the measured offsets.

This issue is outside the 17 listed findings and is not a post-hoc change to the acceptance table. The scientific source and ensemble were kept frozen. A new development/freeze/validation cycle with independent random streams is required before the stochastic results can certify independent instrumental noise. The audit is included in the JSON report.

Numerical diagnoses and implementation limitations are recorded in the iteration-3 section of [NOTES.md](../code/pipeline/NOTES.md). No failed gate was repaired by changing validation seeds or tolerances.
