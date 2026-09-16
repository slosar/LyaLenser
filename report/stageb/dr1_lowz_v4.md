# DR1 Lya forest lensing x low-redshift tracers (single slab 2.1 < z < 3.0)

384014 forests, 5077968 sightline pairs, 11024 deg^2; xi fit b_F^2 = 0.0296, beta_F = 1.037; jackknife nside 8 (301 regions).

| Template | A | jackknife error | sigma_F | curl |
|---|---|---|---|---|
| combined | 0.182 | 0.663 | 0.514 | -0.776 +- 2.169 |
| slice_0.4_0.6 | 0.188 | 1.031 | 0.972 | -1.126 +- 4.304 |
| slice_0.6_0.8 | -1.015 | 1.199 | 0.912 | -0.686 +- 3.693 |
| slice_0.8_1.1 | 2.130 | 1.085 | 0.939 | 2.495 +- 3.985 |
| slice_1.1_1.6 | -1.520 | 1.630 | 1.422 | -7.021 +- 6.329 |
| slice_1.6_1.75 | 1.095 | 5.677 | 5.178 | -8.276 +- 22.723 |

Jackknife-covariance combination of the slices: A = 0.353 +- 0.632.
Random-template null (40 Gaussian realisations of the combined map's spectrum): mean -0.125 +- 0.087, scatter 0.548 against the RMS jackknife error 0.635.
Injection expectation (bookkeeping) odd slopes: {0.25: 1.0419199862945878, 0.5: 1.0418549248103226}.
Wall 41 min, peak 23.5 GB.
