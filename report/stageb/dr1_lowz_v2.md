# DR1 Lya forest lensing x low-redshift tracers (single slab 2.1 < z < 3.0)

373756 forests, 3899521 sightline pairs, 11013 deg^2; xi fit b_F^2 = 0.0295, beta_F = 1.049; jackknife nside 8 (301 regions).

| Template | A | jackknife error | sigma_F | curl |
|---|---|---|---|---|
| combined | -0.086 | 0.746 | 0.651 | -0.197 +- 1.948 |
| slice_0.4_0.6 | -0.642 | 1.227 | 1.219 | -2.254 +- 3.704 |
| slice_0.6_0.8 | -0.564 | 1.533 | 1.169 | 1.731 +- 3.450 |
| slice_0.8_1.1 | 1.259 | 1.318 | 1.196 | 2.746 +- 3.788 |
| slice_1.1_1.6 | -1.183 | 1.981 | 1.788 | -5.922 +- 5.150 |
| slice_1.6_1.75 | 5.705 | 6.865 | 6.541 | -10.884 +- 20.185 |

Jackknife-covariance combination of the slices: A = -0.039 +- 0.708.
Random-template null (40 Gaussian realisations of the combined map's spectrum): mean -0.151 +- 0.113, scatter 0.714 against the RMS jackknife error 0.812.
Injection expectation (bookkeeping) odd slopes: {0.25: 1.0407581337695542, 0.5: 1.0406743981396591}.
Wall 25 min, peak 15.6 GB.
