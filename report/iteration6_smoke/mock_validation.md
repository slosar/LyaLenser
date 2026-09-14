# Stage A iteration 6: SMOKE — not acceptance

Scale 0.25; 2 sparse and 2 dense seeds; 15/36 required gates pass. **Stage B remains blocked.**

All recovery/null statistics are absolute F^-1(q-mf). SEMs and 95% residual bounds use the seed ensemble (every seed carries A_true = 0 and 1, so the null and recovery rows share realisations). Shared/disjoint shifts and cell-offset comparisons are diagnostics. Smoke precision cannot certify scale-1 gates.

| Gate | Measurement | Frozen tolerance | Type | Result |
|---|---|---|---|---|
| truth normalization slope | 1.51428 ± 1.064 SEM; bound95 14.03; N=2 | &#124;slope-1&#124; <= 2 SEM | required | PASS |
| cmb normalization slope | 2.12711 ± 0.8479 SEM; bound95 11.9; N=2 | &#124;slope-1&#124; <= 2 SEM | required | PASS |
| matched normalization slope | 5.27703 ± 32.08 SEM; bound95 411.9; N=2 | &#124;slope-1&#124; <= 2 SEM | required | PASS |
| fixed versus refitted baseline | {"fixed_absolute_slope": {"n": 2, "mean": 1.621803826555416, "sem": 1.038622044194336, "target": 1.0, "residual": 0.621803826555416, "bound95": 13.818748163860272, "slopes": [0.5831817823610801, 2.660425870749752], "intercepts": [3.1970946488486365, 7.754595606512538]}, "refitted_absolute_slope": {"n": 2, "mean": 1.514284858229069, "sem": 1.0637363414001029, "target": 1.0, "residual": 0.5142848582290691, "bound95": 14.030336597642005, "slopes": [0.4505485168289663, 2.578021199629172], "intercepts": [3.1891387118176975, 7.81967659729646]}, "diagnostic_difference": {"n": 2, "mean": -0.10751896832634686, "sem": 0.025114297205766935, "target": 0.0, "residual": -0.10751896832634686, "bound95": 0.4266263704344261, "absolute_A1_refit_minus_fixed": {"n": 2, "mean": -0.0509021626764814, "sem": 0.10163907724863684, "target": 0.0, "residual": -0.0509021626764814, "bound95": 1.3423490874196984}}} | &#124;slope difference&#124; <= 0.03 | required | FAIL |
| varying template mean field: truth | 5.28653 ± 2.051 SEM; bound95 31.35; N=2 | &#124;mean&#124; <= 2 SEM | required | FAIL |
| fixed template mean field: truth | 5.26804 ± 0.1665 SEM; bound95 7.383; N=2 | &#124;mean&#124; <= 2 SEM | required | FAIL |
| varying template mean field: cmb | 47.2098 ± 10.97 SEM; bound95 186.6; N=2 | &#124;mean&#124; <= 2 SEM | required | FAIL |
| fixed template mean field: cmb | 6.40086 ± 10.56 SEM; bound95 140.5; N=2 | &#124;mean&#124; <= 2 SEM | required | PASS |
| varying template mean field: matched | 2.3862 ± 126.3 SEM; bound95 1607; N=2 | &#124;mean&#124; <= 2 SEM | required | PASS |
| fixed template mean field: matched | 2.68845 ± 103.9 SEM; bound95 1322; N=2 | &#124;mean&#124; <= 2 SEM | required | PASS |
| absolute stochastic recovery | 50.9753 ± 6.213 SEM; bound95 128.9; N=2 | &#124;mean-1&#124; <= 2 SEM (precision established separately by slope) | required | FAIL |
| absolute response-only: cmb | 59.3131 ± 9.753 SEM; bound95 183.2; N=2 | &#124;observed-predicted&#124; <= 2 SEM | required | FAIL |
| absolute response-only: matched | 23.5714 ± 125.3 SEM; bound95 1615; N=2 | &#124;observed-predicted&#124; <= 2 SEM | required | PASS |
| absolute response-only: deprojected | 75.0352 ± 13.72 SEM; bound95 249.3; N=2 | &#124;observed-predicted&#124; <= 2 SEM; &#124;mean&#124; <= 2 SEM and absolute null bound95 <= 0.5 A | required | FAIL |
| absolute combined deprojected recovery | 78.4598 ± 12.14 SEM; bound95 231.8; N=2 | &#124;mean-1&#124; <= 2 SEM and residual bound95 <= 0.5 A | required | FAIL |
| paired deprojected response A(1) - A(0), same realisation | 3.42457 ± 1.572 SEM; bound95 22.39; N=2 | &#124;mean-1&#124; <= 2 SEM | required | FAIL |
| six-bin Hotelling shape | {"p_value": 0.0, "T2": NaN, "reason": "insufficient nonsingular ensemble covariance"} | p > 0.01 | required | FAIL |
| absolute covariance | {"scatter": 17.17400373808258, "rms_jk": 30.059331747979403, "ratio": 0.5713368441478085, "jackknife": {"nside": [16], "regions_mean": 4.0}} | 0.7 <= scatter/RMS jackknife <= 1.3 | required | FAIL |
| outside-box and total spectra | {"ratios": [0.9765094706630597, 1.0852849901929744, 1.3082043538207417], "sem": [0.2584788456010089, 0.41503930096458447, 0.5267515024287468]} | each total ratio within 10% | required | FAIL |
| analytic derivative convergence | {"relative_norm_change": 2.7652501514417962e-05} | <1% | required | PASS |
| discrete-grid covariance quadrature | {"relative_norm_change": 7.85163532944286e-07} | <1% | required | PASS |
| unit and regression tests | {"exit_code": 0, "wall_s": 201.42836777993944, "output": ["64 passed in 195.60s (0:03:15)"]} | all pass | required | PASS |
| injection bookkeeping (noise-free expectation, production selection) | {"expected_slope_small_amplitude": 0.9773373643509963, "expected_slopes_by_amplitude": {"0.25": 0.9773373643509963, "0.5": 0.9768636613473187, "1.0": 0.9749791492391623, "2.0": 0.9671665363401791}, "amplitude": 0.25} | odd slope at &#124;A&#124; = 0.25 within 1 ± 0.05 | required | PASS |
| injection convergence and measured slopes | {"expected_by_amplitude": {"0.25": 0.9773373643509963, "0.5": 0.9768636613473187, "1.0": 0.9749791492391623, "2.0": 0.9671665363401791}, "measured_by_amplitude": {"0.25": 0.5892671989611733, "0.5": 0.6917104340582618, "1.0": 0.7461585396451738, "2.0": 0.7184874369743788}, "measured_combined": 0.7209157827750264, "expected_combined": 0.9692131379216274, "measured_minus_expected_by_amplitude": {"0.25": -0.388070165389823, "0.5": -0.2851532272890569, "1.0": -0.22882060959398853, "2.0": -0.24867909936580035}, "curl_slope": 0.09669997784097803, "curl_error": 7.807402477592741} | report: expected slope vs &#124;A&#124;; measured odd slopes; measured minus expected; curl | report only | PASS |
| flag and same-realization margin diagnostics | {"flags": {"magnification": {"absolute_A": -12.757187984304805, "baseline_A": 90.6036356739502, "diagnostic_shift": -103.360823658255}, "completeness": {"absolute_A": -30.017552048623315, "baseline_A": 90.6036356739502, "diagnostic_shift": -120.62118772257351}, "real_mask": {"absolute_A": 41.18794485713939, "baseline_A": 90.6036356739502, "diagnostic_shift": -49.4156908168108}}, "margins": {"0.0": 85.97735042942995, "150.0": 90.6036356739502, "300.0": 39.7422659373436}} | report: persist all (diagnostic) | report only | PASS |
| 100 random-template diagnostic | 3.64593 ± 2.762 SEM; bound95 9.126; N=100 | report: finite statistics; no numerical tolerance | report only | PASS |
| benchmark and memory | {"benchmark": {"sightline_pairs": 402, "pixel_pairs": 46749087, "threads": 24, "wall_s": 0.09258543502073735, "pixel_pairs_per_s": 504929171.5217313, "pixel_pairs_per_s_per_thread": 21038715.480072137, "extrapolated_full_wall_s": 0.9258543502073735}, "peak_gib": 0.2708740234375} | report: positive throughput; requested threads; peak <40 GB | report only | PASS |
| template coefficient continuous (40<=L<=300), margin 0 | 1.03709 ± 0.01129 SEM; bound95 0.1805; N=2 | 1 ± .03 | required | FAIL |
| template coefficient continuous (all pixel scales, diagnostic), margin 0 | 0.992258 ± 0.0003481 SEM; bound95 0.01216; N=2 | report only | report only | PASS |
| template coefficient sampled (40<=L<=300), margin 0 | 1.0434 ± 0.2481 SEM; bound95 3.196; N=2 | 1 within 2 SEM | required | PASS |
| template coefficient sampled (all pixel scales, diagnostic), margin 0 | 0.904479 ± 0.1148 SEM; bound95 1.554; N=2 | report only | report only | PASS |
| template coefficient continuous (40<=L<=300), margin 150 | 1.10539 ± 0.06091 SEM; bound95 0.8793; N=2 | 1 ± .03 | required | FAIL |
| template coefficient continuous (all pixel scales, diagnostic), margin 150 | 1.02693 ± 0.03333 SEM; bound95 0.4504; N=2 | report only | report only | PASS |
| template coefficient sampled (40<=L<=300), margin 150 | 1.11584 ± 0.2195 SEM; bound95 2.905; N=2 | 1 within 2 SEM | required | PASS |
| template coefficient sampled (all pixel scales, diagnostic), margin 150 | 1.08117 ± 0.02377 SEM; bound95 0.3833; N=2 | report only | report only | PASS |
| template coefficient continuous (40<=L<=300), margin 300 | 1.05791 ± 0.03011 SEM; bound95 0.4405; N=2 | 1 ± .03 | required | FAIL |
| template coefficient continuous (all pixel scales, diagnostic), margin 300 | 1.00897 ± 0.01637 SEM; bound95 0.217; N=2 | report only | report only | PASS |
| template coefficient sampled (40<=L<=300), margin 300 | 1.21558 ± 0.2227 SEM; bound95 3.046; N=2 | 1 within 2 SEM | required | PASS |
| template coefficient sampled (all pixel scales, diagnostic), margin 300 | 1.12178 ± 0.08783 SEM; bound95 1.238; N=2 | report only | report only | PASS |
| fitted b_F2 (sparse A=0 samples) versus generator 0.0169 | 0.0161923 ± 0.001301 SEM; bound95 0.01724; N=2 | report only | report only | PASS |
| fitted beta_F (sparse A=0 samples) versus generator 1.6 | 1.6964 ± 0.08836 SEM; bound95 1.219; N=2 | report only | report only | PASS |
| dense physical g=True | 0.883868 ± 0.06224 SEM; bound95 0.9069; N=2 | slope 1 ± .05; >=20 scale-1 seeds | required | FAIL |
| dense absolute endpoints g=True | {"A0_absolute": {"n": 2, "mean": -0.0978794289764386, "sem": 0.3237853887892703, "target": 0.0, "residual": -0.0978794289764386, "bound95": 4.211962869598172}, "A1_absolute": {"n": 2, "mean": 0.8185593903759613, "sem": 0.33114642415772094, "target": 1.0, "residual": -0.18144060962403874, "bound95": 4.3890548727094245}, "difference_of_absolute_means": 0.9164388193523999, "sem_from_absolute_ensembles": 0.4631359759572389} | A0 within 2 SEM of 0 and A1 within 2 SEM of 1 | required | FAIL |
| measured versus fitted table g=True | 1.04353 ± 0.02165 SEM; bound95 0.3187; N=2 | report only (coarse, r_perp^3 weight, fit range) | report only | PASS |
| dense fitted b_F2 versus generator 0.0169, g=True | 0.0191725 ± 0.0001321 SEM; bound95 0.003951; N=2 | report only | report only | PASS |
| dense fitted beta_F versus generator 1.6, g=True | 1.47219 ± 0.01003 SEM; bound95 0.2553; N=2 | report only | report only | PASS |
| fitted-model versus generator-model normalisation g=True | {"fitted": {"n": 2, "mean": 0.883868096735922, "sem": 0.06223636572386231, "target": 1.0, "residual": -0.11613190326407796, "bound95": 0.9069199082029374, "slopes": [0.8216317310120598, 0.9461044624597844], "intercepts": [-0.38036264442157774, 0.2168065182372786]}, "generator": {"n": 2, "mean": 0.9710987567122072, "sem": 0.11204757212123051, "target": 1.0, "residual": -0.028901243287792755, "bound95": 1.4526006348802887, "slopes": [0.8590511845909767, 1.0831463288334378], "intercepts": [-0.442648736881247, 0.5188558944442929]}, "ratio": 0.9101732348297746} | slope ratio within 5% | required | FAIL |
| first moments | {"observed_ratio": 1.0099855274436431, "predicted_ratio": 1.0181514233500892} | omitted/correct slope ratio agrees with prediction within .05 | required | PASS |
| dense physical g=False | 0.866358 ± 0.05396 SEM; bound95 0.8193; N=2 | slope 1 ± .05; >=20 scale-1 seeds | required | FAIL |
| dense absolute endpoints g=False | {"A0_absolute": {"n": 2, "mean": 0.016505808000640013, "sem": 0.22443656570816836, "target": 0.0, "residual": 0.016505808000640013, "bound95": 2.8682427622303224}, "A1_absolute": {"n": 2, "mean": 0.9169536404409683, "sem": 0.21548842391481146, "target": 1.0, "residual": -0.08304635955903172, "bound95": 2.821086392151696}, "difference_of_absolute_means": 0.9004478324403282, "sem_from_absolute_ensembles": 0.31113828576400954} | A0 within 2 SEM of 0 and A1 within 2 SEM of 1 | required | FAIL |
| measured versus fitted table g=False | 1.04353 ± 0.02165 SEM; bound95 0.3187; N=2 | report only (coarse, r_perp^3 weight, fit range) | report only | PASS |
| dense fitted b_F2 versus generator 0.0169, g=False | 0.0191725 ± 0.0001321 SEM; bound95 0.003951; N=2 | report only | report only | PASS |
| dense fitted beta_F versus generator 1.6, g=False | 1.47219 ± 0.01003 SEM; bound95 0.2553; N=2 | report only | report only | PASS |
| fitted-model versus generator-model normalisation g=False | {"fitted": {"n": 2, "mean": 0.8663580874634409, "sem": 0.053962380242265906, "target": 1.0, "residual": -0.1336419125365591, "bound95": 0.8192989639599879, "slopes": [0.8123957072211749, 0.9203204677057067], "intercepts": [-0.1629835729459287, 0.22921718381930178]}, "generator": {"n": 2, "mean": 0.9524444713197556, "sem": 0.10248648638360623, "target": 1.0, "residual": -0.047555528680244374, "bound95": 1.3497698073879052, "slopes": [0.8499579849361494, 1.0549309577033619], "intercepts": [-0.2221833067867608, 0.5410532412229446]}, "ratio": 0.9096153251464318} | slope ratio within 5% | required | FAIL |

## Template coefficients versus margin

Band columns are the gated 40<=L<=300 cross/auto coefficients; the other columns are full-resolution pixel regressions (diagnostic).

| Seed | Margin | Continuous band | Sampled band | Continuous | Sampled | Nominal (no radial normalisation) | Old half-cell offset |
|---|---|---|---|---|---|---|---|
| 1000 | 0 | 1.04838 | 0.795259 | 0.992606 | 1.0193 | 0.975031 | 0.935732 |
| 1000 | 150 | 1.1663 | 0.896358 | 1.06025 | 1.10494 | 1.06504 | 1.00476 |
| 1000 | 300 | 1.08802 | 0.992828 | 1.02534 | 1.20961 | 1.03093 | 0.966823 |
| 1001 | 0 | 1.0258 | 1.29155 | 0.99191 | 0.789658 | 0.998677 | 0.937273 |
| 1001 | 150 | 1.04447 | 1.33532 | 0.993603 | 1.0574 | 1.00575 | 0.948131 |
| 1001 | 300 | 1.0278 | 1.43832 | 0.992591 | 1.03394 | 1.00291 | 0.936593 |

## Shared-selection diagnostic

| Seed | Absolute disjoint A(1), response on | Absolute shared A(1), response on | Shared minus disjoint |
|---|---|---|---|
| 1000 | 90.6036 | 72.3649 | -18.2387 |
| 1001 | 66.3159 | 72.3013 | 5.98536 |

## Reconstruction and limits

Per-seed HDF5 retains maps over each template range, continuum intensity, masks, fitted templates, xi tables, catalogues and q/F/mf partial sums. Completion JSON files contain source/GATES/config provenance and artifact hashes. `collect` validates these before rebuilding this report; it performs no simulation.

See [iteration-4 notes](../code/pipeline/NOTES.md) for the grid-cell diagnosis, derivative comparison approximations and remaining scale-1 checks. Full absolute points, SEMs, predictions and stopping decision are in [JSON](mock_validation.json).
