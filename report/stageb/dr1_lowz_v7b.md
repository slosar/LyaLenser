# DR1 Lya forest lensing x low-redshift tracers (single slab 2.1 < z < 3.0)

427888 forests, 6538923 sightline pairs, 11043 deg^2; xi fit b_F^2 = 0.0551, beta_F = 0.588, gamma_b = 3.83, gamma_beta = -2.07, gamma_S = 6.18 (z_ref 2.4); jackknife nside 8 (301 regions).

| Template | A | jackknife error | sigma_F | curl |
|---|---|---|---|---|
| combined | 0.344 | 0.401 | 0.336 | 0.025 +- 1.888 |
| slice_0.1_0.4 | 0.750 | 0.860 | 0.716 | 0.839 +- 3.632 |
| slice_0.4_0.6 | -0.315 | 0.680 | 0.647 | -1.323 +- 3.314 |
| slice_0.6_0.8 | -0.203 | 0.875 | 0.710 | -1.122 +- 4.468 |
| slice_0.8_1.1 | 1.137 | 0.823 | 0.742 | 2.823 +- 4.827 |
| slice_1.1_1.6 | 1.177 | 1.255 | 1.171 | 2.396 +- 8.777 |

Jackknife-covariance combination of the slices: A = 0.318 +- 0.403.
Random-template null (20 Gaussian realisations of the combined map's spectrum): mean -0.193 +- 0.083, scatter 0.373 against the RMS jackknife error 0.419.
Injection expectation (bookkeeping) odd slopes: {0.25: 1.0244938493334965, 0.5: 1.0243162778997044}.
Wall 73 min, peak 34.4 GB.
