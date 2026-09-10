# Stage A mock validation

Run: `20` seeds, mock scale `0.15`, 24 Numba threads.

| Sec. 2.8 | Acceptance test | Tolerance | Measured | Result |
|---:|---|---|---:|:---:|
| 1 | tables and benchmark | report only | 8.02678e+08 | — |
| 2 | baseline absorption: data/model normalisation | |ratio-1| <= 0.03 | 0.996443 | PASS |
| 3 | physical slope, g on + moments | 1 +/- 0.05 | 1.00212 | PASS |
| 3 | physical slope on fiducial mocks | 1 within jackknife error | 0.467695 +/- 0.097774 | FAIL |
| 3 | physical slope, g off + moments disabled | 1 +/- 0.05 | 0.945098 | FAIL |
| 3 | g on, moments disabled deviation | deviates by predicted amount | 1.00491 | PASS |
| 5 | paired coordinate injection slope | 1 +/- 0.05 | 1.12616 | FAIL |
| 5 | coordinate injection curl | 0 within errors | -0.0742059 | PASS |
| 4 | unlensed truth-template mean | 0 within 2 sigma | -20.1387 | FAIL |
| 4 | unlensed filtered kappa_CMB mean | 0 within 2 sigma | -1.17139 | FAIL |
| 4 | unlensed matched slab mean | 0 within 2 sigma | -1.16305 | FAIL |
| 6 | stochastic Wiener-template normalisation | 1 within 2 sigma | 1.41823 | PASS |
| 7 | response-only predicted amplitude | prediction within errors | 1.84859 | FAIL |
| 7 | deprojected response null | 0 within errors | 19.1751 | PASS |
| 7 | matched-slab response | consistent with prediction | 1.84859 | FAIL |
| 8 | combined deprojected recovery | 1 within errors | 0.756516 | FAIL |
| 8 | flag/slab switch-off shifts | report only | 0.28882 | — |
| 8 | magnification off shift | report only | not implemented | FAIL |
| 8 | completeness off shift | report only | not implemented | FAIL |
| 8 | real-mask off shift | report only | not implemented | FAIL |
| 8 | no template-slab margin shift | report only | not implemented | FAIL |
| 9 | six-bin shape consistency | p > 0.01 | 3.30975e-08 | FAIL |
| 10 | random-template std / jackknife | report only | 2.40012 | — |
| 10 | random-template std / sigma_F | report only | 9.45019 | — |
| 11 | mock scatter / jackknife | 1 +/- 0.3 | 12.278 | FAIL |

## Wall times

- `fiducial_mock_s`: 0.832 s
- `xi_model_s`: 0.567 s
- `xi_data_s`: 5.559 s
- `benchmark`: {'sightline_pairs': 3134, 'pixel_pairs': 416799713, 'threads': 24, 'wall_s': 0.5192614132538438, 'pixel_pairs_per_s': 802678000.6398149, 'pixel_pairs_per_s_per_thread': 33444916.69332562, 'extrapolated_full_wall_s': 0.5192614132538438}
- `ensemble_s`: 221.982 s
- `total_s`: 243.791 s

## Notes

This run uses the reduced patch controlled by `scale`; `/data` was read-only, so large products were written under `/tmp/LyaLenser/mocks`. See `code/pipeline/NOTES.md` for deviations and diagnoses.
