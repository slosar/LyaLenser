# DR1 Lya forest lensing x low-redshift tracers (single slab 2.1 < z < 3.0)

373756 forests, 3899521 sightline pairs, 11013 deg^2; xi fit b_F^2 = 0.0210, beta_F = 1.428; jackknife nside 8 (301 regions).

| Template | A | jackknife error | sigma_F | curl |
|---|---|---|---|---|
| combined | 0.010 | 0.670 | 0.601 | -0.178 +- 1.792 |
| slice_0.4_0.6 | -0.617 | 1.094 | 1.127 | -2.262 +- 3.406 |
| slice_0.6_0.8 | -0.366 | 1.367 | 1.080 | 1.867 +- 3.136 |
| slice_0.8_1.1 | 1.299 | 1.188 | 1.105 | 2.415 +- 3.433 |
| slice_1.1_1.6 | -1.071 | 1.779 | 1.654 | -5.524 +- 4.745 |
| slice_1.6_1.75 | 5.722 | 6.107 | 6.049 | -9.364 +- 18.370 |

Jackknife-covariance combination of the slices: A = 0.002 +- 0.638.
Random-template null (40 Gaussian realisations of the combined map's spectrum): mean -0.125 +- 0.102, scatter 0.644 against the RMS jackknife error 0.731.
Injection expectation (bookkeeping) odd slopes: {0.25: 1.0140220963253734, 0.5: 1.0139434356669481}.
Wall 24 min, peak 15.4 GB.
