# Stage A mock validation

Run: `20` seeds, scale `1.0`, 24 Numba threads.
Mean-field null: `21` seeds (the 20-run ensemble plus fixed seed 20).
Kernel-product `g1 = 0.00026036092 (Mpc/h)^-1` (`<chi_l> = 2007.010 Mpc/h`).

| Sec. 2.8 | Acceptance item | Tolerance | Measured | Result |
|---:|---|---|---|:---:|
| 1 | Tables and benchmark | finite tables; benchmark reported; mock spectra within 10% over 40<L<300 | `{"pixel_pairs_per_s":823956435.6438916,"threads":24,"spectrum_ratios_klkl_klkc_kckc":[0.9800590330254778,0.9818771153345924,1.0278641415412817]}` | PASS |
| 2 | Baseline absorption | data/model slopes agree to 3% | `{"ratio":1.0073566942009105,"data_slope":1.0069481652918078,"model_slope":0.9995944545646498}` | PASS |
| 3 | Physical normalization | slopes 1 +/- 0.05; missing-moment shift predicted | `{"g_on_moments":1.0069481652918078,"g_off":0.9899099800461719,"g_on_no_moments":1.0211157030842621,"predicted_no_moment_ratio":1.02315234058513}` | PASS |
| 4 | Mean field | all three means 0 within 2 SEM | `{"truth":{"mean":0.8106849786776072,"sem":0.5628986714658911,"n":21},"cmb":{"mean":0.7199382949475799,"sem":1.3351575451937787,"n":21},"matched":{"mean":40.89724740486973,"sem":21.57284287849577,"n":21}}` | PASS |
| 5 | Injection bookkeeping | slope 1 +/- 0.05; curl 0 within error | `{"slope":1.0472365095315446,"curl_slope":0.05261668004969638,"curl_error":0.7666693269059167}` | PASS |
| 6 | Stochastic normalization | mean 1 within 2 SEM | `{"mean":-0.5468473267904315,"sem":1.6917185827449814}` | PASS |
| 7 | Response-only null | predictions agree within errors; deprojected is 0 | `{"cmb":{"observed":3.0992093704635044,"observed_sem":2.2978471421122357,"predicted":6.211943213317131,"difference":-3.1127338428536255,"difference_sem":3.384327932447073},"matched":{"observed":42.61611727732173,"observed_sem":8.031241905809733,"predicted":459.033278121108,"difference":-416.4171608437861,"difference_sem":410.33507262058055},"deprojected":{"observed":0.7261116030061608,"observed_sem":1.9104595287951003,"predicted":0.0,"difference":0.7261116030061608,"difference_sem":1.9104595287951003}}` | PASS |
| 8 | Combined recovery | mean 1 within 2 SEM | `{"mean":1.380331791339025,"sem":0.7893546570224247,"switch_off_shifts":{"magnification_off":-0.08853793081971073,"completeness_off":0.5818149275056639,"real_mask_off":0.9185958321402379,"no_margin":0.2020108166690484}}` | PASS |
| 9 | Shape test | six-bin p-value > 0.01 | `{"p_value":0.7829123599694563,"chi2":2.4572948958814185,"values":[1.7873951700741864,-6.7123357578756,0.623398720539827,-0.10866337463205893,0.24262352804903026,0.6114311740947562]}` | PASS |
| 10 | Random-template ensemble | ratios reported (no tolerance) | `{"std":2.671881973458663,"std_over_jackknife":0.9209891680459477,"std_over_sigma_F":2.7704600864715623}` | PASS |
| 11 | Covariance | mock scatter / jackknife = 1 +/- 0.3 | `{"ratio":0.8651708515424091,"scatter":3.530101342916347,"rms_jk":4.080236102063487}` | PASS |

## Wall times and memory

- `normalisation_mock_s`: 6.076 s
- `xi_model_s`: 0.997 s
- `xi_data_s`: 26.538 s
- `full_ensemble_s`: 2508.150 s
- `flag_shifts_s`: 431.561 s
- `mean_field_extra_s`: 133.780 s
- `total_s`: 3217.755 s
- Peak resident memory: 19.474 GB.

## Remaining failures

- None.
