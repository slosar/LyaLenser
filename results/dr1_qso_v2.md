# DR1 quasar x Lya forest lensing x low-redshift tracers

Source-distance treatment: derivative maps of the templates (/data/LyaLenser/lowz_v5); effective ratios combined 2.18e-04, slice_0.1_0.4 6.35e-05, slice_0.4_0.6 1.30e-04, slice_0.6_0.8 2.08e-04, slice_0.8_1.1 3.34e-04, slice_1.1_1.6 6.50e-04 per Mpc/h.
427888 forests, 583148 quasars, 9497606 quasar-sightline pairs; cross fit b_q = 3.204 (gamma_q 1.68), dr_par = 0.44, sigma_par = 4.92 Mpc/h, chi2 5729/11340.

| Template | A | jackknife error | sigma_F | curl |
|---|---|---|---|---|
| combined | 0.877 | 0.470 | 0.424 | -0.432 +- 2.993 |
| slice_0.1_0.4 | 0.954 | 1.110 | 0.912 | 2.212 +- 5.695 |
| slice_0.4_0.6 | 1.001 | 0.965 | 0.821 | -4.706 +- 5.573 |
| slice_0.6_0.8 | 1.449 | 1.163 | 0.911 | 1.813 +- 6.636 |
| slice_0.8_1.1 | 1.032 | 1.016 | 0.917 | 1.414 +- 6.651 |
| slice_1.1_1.6 | -1.523 | 1.656 | 1.371 | -4.317 +- 11.687 |
| sub-slab 1.96-2.25 | 0.972 | 0.889 | 0.706 | 0.815 +- 4.747 |
| sub-slab 2.25-2.55 | 1.390 | 0.861 | 0.677 | 0.517 +- 4.497 |
| sub-slab 2.55-3 | -0.051 | 1.018 | 0.858 | -3.685 +- 6.078 |

Random-template null (20): mean 0.181 +- 0.153, scatter 0.684 vs RMS jackknife 0.579.
Wall 32 min, peak 25.8 GB.
