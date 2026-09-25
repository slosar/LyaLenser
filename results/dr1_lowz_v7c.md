# DR1 Lya forest lensing x low-redshift tracers (single slab 2.1 < z < 3.0)

427888 forests, 6538923 sightline pairs, 11043 deg^2; xi fit b_F^2 = 0.0551, beta_F = 0.588, gamma_b = 3.83, gamma_beta = -2.07, gamma_S = 6.18 (z_ref 2.4); jackknife nside 8 (301 regions).

| Template | A | jackknife error | sigma_F | curl |
|---|---|---|---|---|
| combined | 0.404 | 0.399 | 0.333 | 0.274 +- 2.091 |
| slice_0.1_0.4 | 0.926 | 0.835 | 0.706 | 0.346 +- 3.984 |
| slice_0.4_0.6 | -0.230 | 0.679 | 0.640 | -0.765 +- 3.783 |
| slice_0.6_0.8 | -0.161 | 0.873 | 0.705 | 0.587 +- 5.198 |
| slice_0.8_1.1 | 1.065 | 0.821 | 0.738 | 2.676 +- 5.758 |
| slice_1.1_1.6 | 1.231 | 1.248 | 1.165 | -0.835 +- 10.946 |

Jackknife-covariance combination of the slices: A = 0.389 +- 0.401.
Random-template null (20 Gaussian realisations of the combined map's spectrum): mean -0.132 +- 0.101, scatter 0.453 against the RMS jackknife error 0.411.
Injection expectation (bookkeeping) odd slopes: {0.25: 1.024142924676148, 0.5: 1.0239632526692768}.
Wall 80 min, peak 37.1 GB.
