# Stage A iteration 7 (low-redshift tracers): SMOKE — not acceptance

Scale 0.25; 2 seeds; 19/30 required gates pass. **Stage B remains blocked.**

Absolute statistics F^-1(q-mf); every seed carries A_true = 0 and 1 with the response on (null and recovery rows share realisations). Templates: Wiener-combined lognormal tracers of the same realisation per redshift slice, summed over slices. Smoke precision cannot certify scale-1 gates.

| Gate | Measurement | Frozen tolerance | Type | Result |
|---|---|---|---|---|
| bias from the angular auto-spectrum: ELG_0.8_1.1 | 1 ± 0 SEM; bound95 0; N=2 | &#124;mean-1&#124; <= 2 SEM | required | PASS |
| bias from the angular auto-spectrum: ELG_1.1_1.6 | 1 ± 0 SEM; bound95 0; N=2 | &#124;mean-1&#124; <= 2 SEM | required | PASS |
| bias from the angular auto-spectrum: LRG_0.4_0.6 | 1 ± 0 SEM; bound95 0; N=2 | &#124;mean-1&#124; <= 2 SEM | required | PASS |
| bias from the angular auto-spectrum: LRG_0.6_0.8 | 1.0729 ± 0.008645 SEM; bound95 0.1827; N=2 | &#124;mean-1&#124; <= 2 SEM | required | FAIL |
| bias from the angular auto-spectrum: LRG_0.8_1.1 | 1 ± 0 SEM; bound95 0; N=2 | &#124;mean-1&#124; <= 2 SEM | required | PASS |
| bias from the angular auto-spectrum: QSO_0.8_1.1 | 1 ± 0 SEM; bound95 0; N=2 | &#124;mean-1&#124; <= 2 SEM | required | PASS |
| bias from the angular auto-spectrum: QSO_1.1_1.6 | 1 ± 0 SEM; bound95 0; N=2 | &#124;mean-1&#124; <= 2 SEM | required | PASS |
| bias from the angular auto-spectrum: QSO_1.6_1.75 | 1 ± 0 SEM; bound95 0; N=2 | &#124;mean-1&#124; <= 2 SEM | required | PASS |
| null (A_true = 0), own template, slice0 | -0.389486 ± 4.601 SEM; bound95 58.85; N=2 | &#124;mean&#124; <= 2 SEM | required | PASS |
| recovery (A_true = 1), own template, slice0 | 2.42273 ± 3.862 SEM; bound95 50.49; N=2 | &#124;mean-1&#124; <= 2 SEM | required | PASS |
| null (A_true = 0), own template, slice1 | 5.79632 ± 2.758 SEM; bound95 40.85; N=2 | &#124;mean&#124; <= 2 SEM | required | FAIL |
| recovery (A_true = 1), own template, slice1 | 7.63763 ± 2.519 SEM; bound95 38.65; N=2 | &#124;mean-1&#124; <= 2 SEM | required | FAIL |
| null (A_true = 0), own template, slice2 | 20.2534 ± 10.77 SEM; bound95 157.1; N=2 | &#124;mean&#124; <= 2 SEM | required | PASS |
| recovery (A_true = 1), own template, slice2 | 23.0541 ± 7.485 SEM; bound95 117.2; N=2 | &#124;mean-1&#124; <= 2 SEM | required | FAIL |
| null (A_true = 0), own template, slice3 | -6.37202 ± 16.55 SEM; bound95 216.6; N=2 | &#124;mean&#124; <= 2 SEM | required | PASS |
| recovery (A_true = 1), own template, slice3 | -9.35727 ± 13.9 SEM; bound95 186.9; N=2 | &#124;mean-1&#124; <= 2 SEM | required | PASS |
| null (A_true = 0), own template, slice4 | -161.148 ± 13.23 SEM; bound95 329.3; N=2 | &#124;mean&#124; <= 2 SEM | required | FAIL |
| recovery (A_true = 1), own template, slice4 | -229.76 ± 12.26 SEM; bound95 386.5; N=2 | &#124;mean-1&#124; <= 2 SEM | required | FAIL |
| combined null (A_true = 0), own templates | 2.61686 ± 0.4843 SEM; bound95 8.77; N=2 | &#124;mean&#124; <= 2 SEM and bound95 <= 0.5 A | required | FAIL |
| combined recovery (A_true = 1), own templates | 5.19972 ± 0.6566 SEM; bound95 12.54; N=2 | &#124;mean-1&#124; <= 2 SEM and residual bound95 <= 0.5 A | required | FAIL |
| paired combined response A(1) - A(0), same realisation | 2.58287 ± 1.141 SEM; bound95 16.08; N=2 | &#124;mean-1&#124; <= 2 SEM | required | PASS |
| jackknife-covariance combination of the slice amplitudes minus the combined-template amplitude (A_true = 1) | nan ± nan SEM; bound95 nan; N=0 | report only | report only | FAIL |
| paired truth-template response A(1) - A(0) | 0.575707 ± 0.5059 SEM; bound95 6.853; N=2 | report only | report only | PASS |
| fixed (other-realisation) combined template, A_true = 1 | 3.59414 ± 1.215 SEM; bound95 19.03; N=2 | &#124;mean&#124; <= 2 SEM | required | FAIL |
| curl null, combined template (A_true = 1) | 1.80681 ± 4.777 SEM; bound95 62.51; N=2 | &#124;mean&#124; <= 2 SEM | required | PASS |
| absolute covariance, combined template | {"scatter": 0.9285448631815918, "rms_jk": 9.852723884278175, "ratio": 0.09424245255296915, "jackknife": {"nside": [16], "regions_mean": 4.0}} | 0.7 <= scatter/RMS jackknife <= 1.3 | required | FAIL |
| combined normalisation slope (A grid, response on) | 2.42376 ± 1.337 SEM; bound95 18.41; N=2 | &#124;slope-1&#124; <= 2 SEM | required | PASS |
| truth normalisation slope (A grid, response on) | 0.356735 ± 1.083 SEM; bound95 14.41; N=2 | &#124;slope-1&#124; <= 2 SEM | required | PASS |
| response on minus off, combined template, A_true = 0 (A-grid seeds) | 0.190077 ± 0.9004 SEM; bound95 11.63; N=2 | &#124;mean&#124; <= 2 SEM | required | PASS |
| foreground and total spectra klkl, klkc, kckc | {"ratios": [1.1317599260465823, 1.2187957252382065, 1.065671772047043], "sem": [0.2503610286798069, 0.43070640387691694, 0.3584969477526798]} | each total ratio within 10% | required | FAIL |
| unit and regression tests | {"exit_code": 0, "wall_s": 363.28954350599815, "output": ["73 passed in 357.07s (0:05:57)"]} | all pass | required | PASS |
| injection bookkeeping (noise-free expectation, production selection) | {"expected_slope_small_amplitude": 0.9801337428437474, "expected_slopes_by_amplitude": {"0.25": 0.9801337428437474, "0.5": 0.9796325106717235, "1.0": 0.9778331736761744, "2.0": 0.9707162443214317}, "amplitude": 0.25} | odd slope at &#124;A&#124; = 0.25 within 1 ± 0.05 | required | PASS |
| injection convergence and measured slopes | {"expected_by_amplitude": {"0.25": 0.9801337428437474, "0.5": 0.9796325106717235, "1.0": 0.9778331736761744, "2.0": 0.9707162443214317}, "measured_by_amplitude": {"0.25": 1.022019475664937, "0.5": 1.453192287762021, "1.0": 1.2137162190790902, "2.0": 1.0284871615431455}, "measured_combined": 1.083263 | report: expected slope vs &#124;A&#124;; measured odd slopes; curl | report only | PASS |
| benchmark and memory | {"benchmark": {"sightline_pairs": 271, "pixel_pairs": 32168622, "threads": 8, "wall_s": 0.15362361900042742, "pixel_pairs_per_s": 209398933.63604784, "pixel_pairs_per_s_per_thread": 26174866.70450598, "extrapolated_full_wall_s": 1.5362361900042742}, "peak_gib": 0.3007774353027344} | report: positive throughput; requested threads; peak <40 GB | report only | PASS |

## Per-seed combined amplitudes

| Seed | A(0) combined | A(1) combined | A(1) truth | jackknife combination A(1) ± err | slices A(1) |
|---|---|---|---|---|---|
| 4000 | 2.133 | 5.856 | 5.553 | nan ± nan | 6.28, 10.16, 15.57, -23.25, -217.50 |
| 4001 | 3.101 | 4.543 | 10.014 | nan ± nan | -1.44, 5.12, 30.54, 4.54, -242.02 |

## Tracer bias fits (first seed)

| Tracer | b true | b fit ± sigma | objects |
|---|---|---|---|
| ELG_0.8_1.1 | 1.30 | 0.895 ± 0.274 | 8124 |
| ELG_1.1_1.6 | 1.40 | 0.764 ± 0.334 | 11705 |
| LRG_0.4_0.6 | 1.90 | 1.377 ± 0.313 | 4482 |
| LRG_0.6_0.8 | 2.10 | 2.271 ± 0.331 | 5666 |
| LRG_0.8_1.1 | 2.30 | 1.411 ± 0.393 | 3481 |
| QSO_0.8_1.1 | 1.70 | 0.000 ± 921450.161 | 989 |
| QSO_1.1_1.6 | 2.10 | 0.000 ± 1712547.505 | 1367 |
| QSO_1.6_1.75 | 2.50 | 0.000 ± 2585441.276 | 537 |
