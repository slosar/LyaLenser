# DR1 Lya forest lensing x low-redshift tracers (single slab 2.1 < z < 3.0)

427888 forests, 6538923 sightline pairs, 11043 deg^2; xi fit b_F^2 = 0.0551, beta_F = 0.588, gamma_b = 3.83, gamma_beta = -2.07, gamma_S = 6.18 (z_ref 2.4); jackknife nside 8 (301 regions).

| Template | A | jackknife error | sigma_F | curl |
|---|---|---|---|---|
| combined | 0.346 | 0.474 | 0.379 | -0.582 +- 1.603 |
| slice_0.1_0.4 | 0.827 | 1.034 | 0.820 | 1.459 +- 3.414 |
| slice_0.4_0.6 | -0.509 | 0.725 | 0.735 | -1.790 +- 3.106 |
| slice_0.6_0.8 | 0.113 | 1.039 | 0.797 | -1.054 +- 3.431 |
| slice_0.8_1.1 | 1.388 | 0.947 | 0.827 | 1.905 +- 3.679 |
| slice_1.1_1.6 | 0.024 | 1.431 | 1.295 | -6.874 +- 5.775 |

Jackknife-covariance combination of the slices: A = 0.208 +- 0.461.
Random-template null (40 Gaussian realisations of the combined map's spectrum): mean -0.063 +- 0.071, scatter 0.452 against the RMS jackknife error 0.485.
Injection expectation (bookkeeping) odd slopes: {0.25: 1.0257645318400148, 0.5: 1.0256102965605147}.
Wall 78 min, peak 34.3 GB.
