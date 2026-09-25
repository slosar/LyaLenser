# DR1 Lya forest lensing x low-redshift tracers (single slab 2.1 < z < 3.0)

384014 forests, 5077968 sightline pairs, 11024 deg^2; xi fit b_F^2 = 0.0733, beta_F = 0.473, gamma_b = 3.49, gamma_beta = -1.37, gamma_S = 5.03 (z_ref 2.4); jackknife nside 8 (301 regions).

| Template | A | jackknife error | sigma_F | curl |
|---|---|---|---|---|
| combined | 0.271 | 0.519 | 0.398 | -0.271 +- 1.700 |
| slice_0.1_0.4 | 1.210 | 1.110 | 0.868 | 2.012 +- 3.579 |
| slice_0.4_0.6 | -0.400 | 0.781 | 0.776 | -1.224 +- 3.264 |
| slice_0.6_0.8 | -0.612 | 1.099 | 0.839 | -0.897 +- 3.581 |
| slice_0.8_1.1 | 1.509 | 0.987 | 0.867 | 2.077 +- 3.867 |
| slice_1.1_1.6 | -0.538 | 1.494 | 1.339 | -6.542 +- 6.103 |
| slice_1.6_1.75 | 0.736 | 5.632 | 5.230 | -6.085 +- 24.836 |

Jackknife-covariance combination of the slices: A = 0.220 +- 0.512.
Random-template null (40 Gaussian realisations of the combined map's spectrum): mean -0.094 +- 0.071, scatter 0.449 against the RMS jackknife error 0.514.
Injection expectation (bookkeeping) odd slopes: {0.25: 1.021741194299939, 0.5: 1.0215377034330586}.
Wall 57 min, peak 23.6 GB.
