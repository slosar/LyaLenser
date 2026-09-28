# DR1 Lya forest lensing x low-redshift tracers (single slab 2.1 < z < 3.0)

Source-distance treatment: derivative maps of the templates (/data/LyaLenser/lowz_v5); effective ratios combined 2.19e-04, slice_0.1_0.4 6.35e-05, slice_0.4_0.6 1.30e-04, slice_0.6_0.8 2.08e-04, slice_0.8_1.1 3.34e-04, slice_1.1_1.6 6.50e-04 per Mpc/h.
427888 forests, 6538923 sightline pairs, 11043 deg^2; xi fit b_F^2 = 0.0551, beta_F = 0.588, gamma_b = 3.83, gamma_beta = -2.07, gamma_S = 6.18 (z_ref 2.4); jackknife nside 8 (301 regions).

| Template | A | jackknife error | sigma_F | curl |
|---|---|---|---|---|
| combined | 0.396 | 0.427 | 0.355 | 0.048 +- 2.126 |
| slice_0.1_0.4 | 0.832 | 0.917 | 0.768 | 0.858 +- 4.072 |
| slice_0.4_0.6 | -0.192 | 0.728 | 0.689 | -1.603 +- 3.832 |
| slice_0.6_0.8 | -0.263 | 0.951 | 0.764 | -1.162 +- 5.185 |
| slice_0.8_1.1 | 1.136 | 0.822 | 0.768 | 2.872 +- 5.235 |
| slice_1.1_1.6 | 1.108 | 1.232 | 1.142 | 2.074 +- 8.638 |

Jackknife-covariance combination of the slices: A = 0.395 +- 0.429.
Random-template null (20 Gaussian realisations of the combined map's spectrum): mean -0.202 +- 0.088, scatter 0.395 against the RMS jackknife error 0.446.
Injection expectation (bookkeeping) odd slopes: {0.25: 1.0249800293294333, 0.5: 1.024819132873986}.
Wall 62 min, peak 25.0 GB.
