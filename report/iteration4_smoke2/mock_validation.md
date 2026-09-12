# Stage A iteration 4: SMOKE — not acceptance

Scale 0.25; 2 sparse and 2 dense seeds. **Stage B remains blocked.**

All recovery/null statistics are absolute F^-1(q-mf). SEMs and 95% residual bounds use independent seed ensembles. Shared/disjoint shifts and cell-offset comparisons are diagnostics. Smoke precision cannot certify scale-1 gates.

| Gate | Measurement | Frozen tolerance | Result |
|---|---|---|---|
| truth normalization slope | -0.0269536 ± 0.1371 SEM; bound95 2.769; N=2 | &#124;slope-1&#124; <= 2 SEM | FAIL |
| cmb normalization slope | -0.0687912 ± 0.2656 SEM; bound95 4.443; N=2 | &#124;slope-1&#124; <= 2 SEM | FAIL |
| matched normalization slope | -8.89756 ± 2.891 SEM; bound95 46.64; N=2 | &#124;slope-1&#124; <= 2 SEM | FAIL |
| fixed versus refitted baseline | {"fixed_absolute_slope": {"n": 2, "mean": -0.06152112924558406, "sem": 0.10407014663131582, "target": 1.0, "residual": -1.061521129245584, "bound95": 2.3838577192935917, "slopes": [0.04254901738573178, -0.1655912758768999], "intercepts": [0.1530822647299607, -0.10300342097283946]}, "refitted_absolute_slope": {"n": 2, "mean": -0.026953593381722818, "sem": 0.13711187228655153, "target": 1.0, "residual": -1.0269535933817229, "bound95": 2.7691251144501763, "slopes": [0.11015827890482871, -0.16406546566827435], "intercepts": [0.16129791401142524, -0.08554815127080076]}, "diagnostic_difference": {"n": 2, "mean": 0.03456753586386124, "sem": 0.03304172565523569, "target": 0.0, "residual": 0.03456753586386124, "bound95": 0.4544024668843069, "absolute_A1_refit_minus_fixed": {"n": 2, "mean": 0.06277959100320436, "sem": 0.029084871462036876, "target": 0.0, "residual": 0.06277959100320436, "bound95": 0.43233792253265596}}} | &#124;slope difference&#124; <= 0.03 | FAIL |
| varying template mean field: truth | 0.0030569 ± 0.1735 SEM; bound95 2.207; N=2 | &#124;mean&#124; <= 2 SEM | PASS |
| fixed template mean field: truth | -0.551506 ± 0.08198 SEM; bound95 1.593; N=2 | &#124;mean&#124; <= 2 SEM | FAIL |
| varying template mean field: cmb | -0.549286 ± 1.483 SEM; bound95 19.39; N=2 | &#124;mean&#124; <= 2 SEM | PASS |
| fixed template mean field: cmb | 4.48471 ± 2.792 SEM; bound95 39.96; N=2 | &#124;mean&#124; <= 2 SEM | PASS |
| varying template mean field: matched | -4.90732 ± 13.68 SEM; bound95 178.8; N=2 | &#124;mean&#124; <= 2 SEM | PASS |
| fixed template mean field: matched | -0.748539 ± 11.29 SEM; bound95 144.2; N=2 | &#124;mean&#124; <= 2 SEM | PASS |
| absolute stochastic recovery | 0.302158 ± 1.396 SEM; bound95 18.43; N=2 | &#124;mean-1&#124; <= 2 SEM (precision established separately by slope) | PASS |
| absolute response-only: cmb | -0.0828685 ± 0.2877 SEM; bound95 3.738; N=2 | &#124;observed-predicted&#124; <= 2 SEM | PASS |
| absolute response-only: matched | -5.13196 ± 19.62 SEM; bound95 254.4; N=2 | &#124;observed-predicted&#124; <= 2 SEM | PASS |
| absolute response-only: deprojected | 0.0350273 ± 0.8459 SEM; bound95 10.78; N=2 | &#124;observed-predicted&#124; <= 2 SEM; absolute null bound95 <= 0.3 A | FAIL |
| absolute combined deprojected recovery | 1.27899 ± 1.156 SEM; bound95 14.97; N=2 | &#124;mean-1&#124; <= 2 SEM and residual bound95 <= 0.3 A | FAIL |
| six-bin Hotelling shape | {"p_value": 0.0, "T2": NaN, "reason": "insufficient nonsingular ensemble covariance"} | p > 0.01 | FAIL |
| absolute covariance | {"scatter": 1.635110330376695, "rms_jk": 8.143742232113077, "ratio": 0.20078119908179226} | 0.7 <= scatter/RMS jackknife <= 1.3 | FAIL |
| outside-box and total spectra | {"ratios": [0.8677684480163544, 1.0891624815996857, 1.5467345131537378], "sem": [0.07565282801504092, 0.12834183012374084, 0.14680334903595238]} | each total ratio within 10% | FAIL |
| analytic derivative convergence | {"relative_norm_change": 2.7652501512514845e-05} | <1% | PASS |
| discrete-grid covariance quadrature | {"relative_norm_change": 7.851635329442861e-07} | <1% | PASS |
| unit and regression tests | {"exit_code": 0, "wall_s": 306.42261356301606, "output": ["51 passed in 297.52s (0:04:57)"]} | all pass | PASS |
| injection bookkeeping | {"slope": 1.5255066529343249, "curl_slope": -0.912853958856183, "curl_error": 1.6787032217539062} | science slope within .05; curl within JK | FAIL |
| flag and same-realization margin diagnostics | {"flags": {"magnification": {"absolute_A": -0.9478483458453829, "baseline_A": 2.435187041157394, "diagnostic_shift": -3.383035387002777}, "completeness": {"absolute_A": 1.123696909060568, "baseline_A": 2.435187041157394, "diagnostic_shift": -1.311490132096826}, "real_mask": {"absolute_A": 1.0114518567024735, "baseline_A": 2.435187041157394, "diagnostic_shift": -1.4237351844549204}}, "margins": {"0.0": 1.3686626402504427, "150.0": 2.435187041157394, "300.0": 1.08867440395811}} | persist all (diagnostic) | PASS |
| 100 random-template diagnostic | 0.348002 ± 0.3774 SEM; bound95 1.097; N=100 | finite statistics; no numerical tolerance | PASS |
| benchmark and memory | {"benchmark": {"sightline_pairs": 311, "pixel_pairs": 37451832, "threads": 8, "wall_s": 0.2787547241896391, "pixel_pairs_per_s": 134354070.98077095, "pixel_pairs_per_s_per_thread": 16794258.87259637, "extrapolated_full_wall_s": 2.787547241896391}, "peak_gib": 0.24979400634765625} | positive throughput; requested threads; peak <40 GB | PASS |
| template coefficient continuous (40<=L<=300), margin 0 | 0.963633 ± 0.03346 SEM; bound95 0.4616; N=2 | 1 ± .03 | FAIL |
| template coefficient continuous (all pixel scales, diagnostic), margin 0 | 0.971676 ± 0.0039 SEM; bound95 0.07787; N=2 | report only | PASS |
| template coefficient sampled (40<=L<=300), margin 0 | 1.04602 ± 0.6063 SEM; bound95 7.75; N=2 | 1 within SEM | PASS |
| template coefficient sampled (all pixel scales, diagnostic), margin 0 | 0.94781 ± 0.107 SEM; bound95 1.412; N=2 | report only | PASS |
| template coefficient continuous (40<=L<=300), margin 150 | 0.929191 ± 0.002637 SEM; bound95 0.1043; N=2 | 1 ± .03 | FAIL |
| template coefficient continuous (all pixel scales, diagnostic), margin 150 | 0.958896 ± 0.004382 SEM; bound95 0.09678; N=2 | report only | PASS |
| template coefficient sampled (40<=L<=300), margin 150 | 1.10634 ± 0.3103 SEM; bound95 4.049; N=2 | 1 within SEM | PASS |
| template coefficient sampled (all pixel scales, diagnostic), margin 150 | 0.910745 ± 0.1053 SEM; bound95 1.428; N=2 | report only | PASS |
| template coefficient continuous (40<=L<=300), margin 300 | 0.935956 ± 0.006436 SEM; bound95 0.1458; N=2 | 1 ± .03 | FAIL |
| template coefficient continuous (all pixel scales, diagnostic), margin 300 | 0.978177 ± 0.002092 SEM; bound95 0.04841; N=2 | report only | PASS |
| template coefficient sampled (40<=L<=300), margin 300 | 1.00537 ± 0.02976 SEM; bound95 0.3835; N=2 | 1 within SEM | PASS |
| template coefficient sampled (all pixel scales, diagnostic), margin 300 | 1.07193 ± 0.04549 SEM; bound95 0.65; N=2 | report only | PASS |
| dense physical g=True | 0.94504 ± 0.01664 SEM; bound95 0.2664; N=2 | slope 1 ± .03; >=10 scale-1 seeds | FAIL |
| dense absolute endpoints g=True | {"A0_absolute": {"n": 2, "mean": 0.8095664697038232, "sem": 1.6084709396016879, "target": 0.0, "residual": 0.8095664697038232, "bound95": 21.247127540884172}, "A1_absolute": {"n": 2, "mean": 1.7129678832223012, "sem": 1.6290841764651869, "target": 1.0, "residual": 0.7129678832223012, "bound95": 21.41244496227084}, "difference_of_absolute_means": 0.903401413518478, "sem_from_absolute_ensembles": 2.289343577873883} | A0=0 and A1=1 within .03; difference of absolute means 1 ± .03 | FAIL |
| measured versus grid-predicted derivative g=True | 0.978383 ± 0.004563 SEM; bound95 0.07959; N=2 | response integral residual <3% | PASS |
| xi interpolation g=True | {"changes": [0.0013218777285716893, 0.0013830674242980001]} | all absolute amplitude changes <.005 | PASS |
| data versus analytic baseline g=True | {"data": {"n": 2, "mean": 0.945039962743588, "sem": 0.016642755615126514, "target": 1.0, "residual": -0.054960037256411964, "bound95": 0.2664262974806143, "slopes": [0.9283972071284615, 0.9616827183587145], "intercepts": [-0.8149846631247594, 2.40561593655632]}, "analytic": {"n": 2, "mean": 0.9004092503775462, "sem": 0.012329776994417507, "target": 1.0, "residual": -0.09959074962245384, "bound95": 0.25625542046807304, "slopes": [0.8880794733831286, 0.9127390273719637], "intercepts": [-1.0451328910109845, 2.346634714817503]}, "ratio": 1.0495671411053673} | slope ratio within 3% | FAIL |
| first moments | {"observed_ratio": 1.0260034372310163, "predicted_ratio": 1.0181480053914966} | omitted/correct slope ratio agrees with prediction within .05 | PASS |
| dense physical g=False | 0.957188 ± 0.02538 SEM; bound95 0.3654; N=2 | slope 1 ± .03; >=10 scale-1 seeds | FAIL |
| dense absolute endpoints g=False | {"A0_absolute": {"n": 2, "mean": 0.805081395702399, "sem": 1.6243228623689594, "target": 0.0, "residual": 0.805081395702399, "bound95": 21.44406024302981}, "A1_absolute": {"n": 2, "mean": 1.7257027119649986, "sem": 1.6502629778630116, "target": 1.0, "residual": 0.7257027119649986, "bound95": 21.69428197764653}, "difference_of_absolute_means": 0.9206213162625996, "sem_from_absolute_ensembles": 2.315554503206475} | A0=0 and A1=1 within .03; difference of absolute means 1 ± .03 | FAIL |
| measured versus grid-predicted derivative g=False | 0.978383 ± 0.004563 SEM; bound95 0.07959; N=2 | response integral residual <3% | PASS |
| xi interpolation g=False | {"changes": [0.0015519597901489113, 0.0013991230794441034]} | all absolute amplitude changes <.005 | PASS |
| data versus analytic baseline g=False | {"data": {"n": 2, "mean": 0.9571884608302574, "sem": 0.025384526759708592, "target": 1.0, "residual": -0.04281153916974256, "bound95": 0.3653525333160392, "slopes": [0.9318039340705488, 0.982572987589966], "intercepts": [-0.8321224152275171, 2.4171027281371664]}, "analytic": {"n": 2, "mean": 0.9099820486249333, "sem": 0.02189285168568322, "target": 1.0, "residual": -0.09001795137506674, "bound95": 0.36819300715770026, "slopes": [0.88808919693925, 0.9318749003106165], "intercepts": [-1.0680060298410692, 2.3646809698075906]}, "ratio": 1.0518762015983254} | slope ratio within 3% | FAIL |

## Template coefficients versus margin

Band columns are the gated 40<=L<=300 cross/auto coefficients; the other columns are full-resolution pixel regressions (diagnostic).

| Seed | Margin | Continuous band | Sampled band | Continuous | Sampled | Nominal (no radial normalisation) | Old half-cell offset |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 0.930168 | 0.439729 | 0.975576 | 0.840769 | 1.00138 | 0.922896 |
| 0 | 150 | 0.931829 | 1.41665 | 0.963277 | 1.01609 | 0.977131 | 0.905184 |
| 0 | 300 | 0.92952 | 1.03513 | 0.98027 | 1.02644 | 0.984244 | 0.92261 |
| 1 | 0 | 0.997097 | 1.65232 | 0.967777 | 1.05485 | 0.942277 | 0.906882 |
| 1 | 150 | 0.926554 | 0.796037 | 0.954514 | 0.805395 | 0.940062 | 0.899031 |
| 1 | 300 | 0.942392 | 0.975619 | 0.976085 | 1.11742 | 0.964406 | 0.915175 |

## Shared-selection diagnostic

| Seed | Absolute disjoint A(1), response on | Absolute shared A(1), response on | Shared minus disjoint |
|---|---|---|---|
| 0 | 2.43519 | 1.75991 | -0.675278 |
| 1 | 0.122792 | 0.751534 | 0.628742 |

## Reconstruction and limits

Per-seed HDF5 retains maps over each template range, continuum intensity, masks, fitted templates, xi tables, catalogues and q/F/mf partial sums. Completion JSON files contain source/GATES/config provenance and artifact hashes. `collect` validates these before rebuilding this report; it performs no simulation.

See [iteration-4 notes](../code/pipeline/NOTES.md) for the grid-cell diagnosis, derivative comparison approximations and remaining scale-1 checks. Full absolute points, SEMs, predictions and stopping decision are in [JSON](mock_validation.json).
