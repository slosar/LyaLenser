# DR1 Lya forest lensing x low-redshift tracers (single slab 2.1 < z < 3.0)

427888 forests, 6538923 sightline pairs, 11043 deg^2; xi fit b_F^2 = 0.0551, beta_F = 0.588, gamma_b = 3.83, gamma_beta = -2.07, gamma_S = 6.18 (z_ref 2.4); jackknife nside 8 (301 regions).

| Template | A | jackknife error | sigma_F | curl |
|---|---|---|---|---|
| combined | 0.395 | 0.428 | 0.355 | 0.057 +- 2.120 |
| slice_0.1_0.4 | 0.826 | 0.911 | 0.763 | 0.877 +- 4.052 |
| slice_0.4_0.6 | -0.194 | 0.729 | 0.686 | -1.576 +- 3.820 |
| slice_0.6_0.8 | -0.266 | 0.950 | 0.762 | -1.174 +- 5.175 |
| slice_0.8_1.1 | 1.166 | 0.824 | 0.770 | 2.902 +- 5.247 |
| slice_1.1_1.6 | 1.090 | 1.249 | 1.163 | 2.263 +- 8.739 |

Jackknife-covariance combination of the slices: A = 0.398 +- 0.431.
Random-template null (20 Gaussian realisations of the combined map's spectrum): mean -0.202 +- 0.088, scatter 0.394 against the RMS jackknife error 0.446.
Injection expectation (bookkeeping) odd slopes: {0.25: 1.024694885764013, 0.5: 1.0245342254847003}.
Wall 75 min, peak 34.4 GB.
