# DR1 Lya forest lensing x low-redshift tracers (single slab 2.1 < z < 3.0)

373756 forests, 3899521 sightline pairs, 11013 deg^2; xi fit b_F^2 = 0.0300, beta_F = 1.024; jackknife nside 8 (301 regions).

| Template | A | jackknife error | sigma_F | curl |
|---|---|---|---|---|
| combined | -0.101 | 0.754 | 0.646 | -0.201 +- 1.961 |
| slice_0.4_0.6 | -0.841 | 1.226 | 1.211 | -2.566 +- 3.732 |
| slice_0.6_0.8 | -0.592 | 1.526 | 1.161 | 1.834 +- 3.459 |
| slice_0.8_1.1 | 1.431 | 1.312 | 1.187 | 2.931 +- 3.747 |
| slice_1.1_1.6 | -1.254 | 1.950 | 1.776 | -6.104 +- 5.188 |
| slice_1.6_1.75 | 6.531 | 6.835 | 6.499 | -9.356 +- 19.918 |

Jackknife-covariance combination of the slices: A = -0.082 +- 0.715.
Random-template null (40 Gaussian realisations of the combined map's spectrum): mean -0.119 +- 0.115, scatter 0.729 against the RMS jackknife error 0.811.
Injection expectation (bookkeeping) odd slopes: {0.25: 1.0381806753160467, 0.5: 1.038126182124997}.
Wall 25 min, peak 15.7 GB.
