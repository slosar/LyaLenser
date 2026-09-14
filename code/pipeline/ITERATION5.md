# Stage A, iteration 5 brief (2026-09-14; implemented by the Claude session, to be reviewed by gpt-6-astra)

Starting point: iteration 4 at scale 1, 39/48 gates (report/mock_validation.md @ 261e07b; NOTES.md "Scale-1
campaign"). Solved there: template prerequisites (continuous 1.004 +- 0.002), xi' model (0.9955 +- 0.002),
deprojection unbiased (null -0.63 +- 0.67, recovery 0.70 +- 0.64, covariance ratio 0.95). Open: a 4.2 +- 1.6 %
multiplicative normalisation bias, exactly linear in A, traced to the smoothed numerical derivative of the measured
xi (development slopes 0.52-0.91 across the smoothing widths; freeze at the grid edge); plus gates whose tolerances
were unattainable at the ensemble size (0.3 A bounds), ill-posed (one-SEM sampled coefficient, 0.03 on means with
SEM 0.16) or with an inadequate reference (continuum "analytic baseline", ratio 1.117).

User decisions (2026-09-14): Claude implements, astra reviews afterwards; tolerances matched to N; the general
Arinyo-capable model fit is implemented now (mocks use the generator's Kaiser + cutoff model); a ~5 % normalisation
systematic is acceptable if understood — do not fight it ad infinitum; report where the signal comes from in
separation.

## Changes
1. `xi_fit.py`: basis spectra (`BasisPower`: P_lin F_NL mu^{2i} with optional line-of-sight pixel/resolution
   window), `basis_tables` (providers `grid` = mock discrete-grid covariance with trilinear window and aliases,
   `hankel` = analytic transform with J0/J1), `project_fine` (per-forest mean+slope continuum projection of a
   stationary table on the actual radial pixel grid, applied identically to xi and dxi/dr_perp), `coarse_bin`,
   `fit_model_table` (weighted LSQ of (b_F^2, beta_F); a free 3-coefficient linear fit is reported as a
   model-adequacy diagnostic). `grid_covariance.xi_from_mock_grid(power=...)`. `mock.grid_geometry` shared by the
   generator and the basis phase.
2. `Config`: `forest_model`, `fit_rperp_min` (3), `r_perp_min` (3, enforced in the pair kernel), `los_pixel`,
   `los_resolution`.
3. `campaign4.py` (iteration 5): phase `basis`; `dev_seed` fits per A and records the parameters (no width scan);
   `freeze` records them; every sample's table comes from `run_mock_validation.table_for`; dense seeds compare the
   fitted-table kernel with the generator-model kernel on the same catalogues; collect rows per GATES.md v5.
4. `signal_profile.py` + `report/signal_profile.json|figures/signal_profile.pdf`: Fisher-information density.
5. GATES.md v5; drivers gain the `basis` phase.

## Verification before the campaign
Smoke chain at scale 0.25 (RACF HTCondor) end to end; unit tests `tests/test_iteration5.py` (analytic derivative
vs finite difference, basis recombination, projection commutation and constant removal, fit recovery, kernel cut,
grid geometry); the fitted table on the iteration-4 dense smoke mock: b_F^2 0.0179 (truth 0.0169), beta_F 1.54
(1.60), per-band shape residuals <= 3 %.

## Deliverables
GATES v5 frozen, scale-1 campaign on Perlmutter (`slurm/campaign4_perlmutter.sh`), `report/mock_validation.md`,
NOTES.md iteration-5 section with the acceptance reading and the understood normalisation, then the astra review of
rounds 4+5.
