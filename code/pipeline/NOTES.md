# Stage A implementation notes

## Environment

- `/data` is mounted read-only in this execution environment. Creating
  `/data/LyaLenser/mocks/` failed with `Read-only file system`, so all large
  HDF5 products and logs were written to `/tmp/LyaLenser/mocks/`. No large
  intermediate was written into the repository.
- All requested packages were already importable from
  `/home/anze/anaconda3/bin/python3`; no `.venv` was created.
- Validation used `NUMBA_NUM_THREADS=24` and `OMP_NUM_THREADS=24`.

## Ambiguities and implementation choices

- The specification defines `g1` through a mean lens distance but does not
  define the lens-distribution weighting or integration limits numerically.
  `Config.g1` uses a single effective lens plane at z=1.0. This must be
  replaced by the stated three-tracer kernel-product average before Stage B.
- Shape-bin boundary values at exactly 10, 20, and 30 Mpc/h are assigned to
  the upper bin (30 is retained). This follows the half-open notation except
  that the global accepted-cut boundary is included.
- `xi_from_data` uses 1 Mpc/h counts, a 0.45-bin Gaussian normalized smoother,
  and a cubic bivariate spline. The smoothing width was selected using the
  fixed-geometry normalization check; it is not a fitted science parameter.
- `three_tracer.spectra()` now accepts explicit forest slab, template slab,
  forest source, background source, bias, magnification, density, and L-grid
  arguments. Calling it with no arguments retains the original fiducial
  calculation and result keys.

## Known deviations from IMPLEMENTATION.md

- The completed 20-seed run used `scale=0.15` (and `scale=0.05` for the costly
  100/deg2 normalization sample), not the full 20-degree `scale=1` geometry.
  The full run was not executed. The 20-seed run is therefore a reduced-scale
  validation, not the full Stage A acceptance run.
- Reduced patches occupy too few `nside_jk=8` HEALPix pixels. The validation
  runner uses deterministic pair-index regions; the amplitude API itself
  accepts the required midpoint HEALPix regions and `pairs.py` provides the
  midpoint-region function.
- The FFT mock uses the requested 2 Mpc/h transverse and 0.5 Mpc/h radial
  cells, but its finite-box forest normalization is empirically matched to the
  continuum-projected Kaiser table. It does not implement the Arinyo nonlinear
  correction in mock generation.
- The outside-box convergence pair is an explicitly correlated, smoothed
  Gaussian approximation. It is not drawn from the three separate restricted
  `three_tracer.limber` spectra.
- Template quasars/randoms are currently uniform catalogues. Their lognormal
  clustering, radial RSD displacement, magnification modulation, and
  completeness selection are not applied. The corresponding flags are stored;
  only the approximate footprint mask currently changes a map.
- `matched_template` implements the specified object weights, random
  subtraction, smoothed completeness division, mask and shot-noise estimate,
  but the validation runner does not exercise it end to end.
- The two-pointer pair kernel implements the 11 accumulators and six bins, but
  it presently returns one in-memory `PairCatalogue`; it does not split output
  into one HDF5 group per sub-slab. Pixel `slab` labels are not applied in the
  kernel.
- `amplitude()` computes the full multi-template response matrix but only
  records whether junk/curl templates were supplied; it does not reject a fit
  lacking the mandatory junk band. The validation fits do not construct true
  harmonic science/junk bands.
- Mask-edge rejection within two degrees is not implemented.
- The response-only prediction in the validation output is an empirical
  paired-mock prediction, not the required independent prediction from
  `R_delta`, spectra, and the catalogue response.
- The 100-member random-template diagnostic uses Gaussian sightline
  deflections with the measured component variance rather than harmonic
  kappa realizations having the CMB spectrum.
- The acceptance runner reports only a combined all-flags-on versus all-off
  shift. Separate magnification, completeness, mask, and no-margin shifts are
  not yet physically meaningful because those catalogue effects are not
  implemented.

## Tests and validation status

- Unit tests: 12 passed, 0 failed. These cover the independent brute-force
  kernel comparison, bitwise thread determinism, pair-direction finite
  difference and sign, xi normalization/derivative/grid sampling, the
  `alm2map_der1` convention and finite-difference divergence, analytic
  flat-sky sign, amplitude/response matrix formulas, per-sightline compression,
  first moments, and the injection coordinate sign.
- The final reduced 20-seed JSON has 6 passing, 12 failing, and 4 report-only
  rows after incomplete proxy tests are conservatively counted as failures.
  The failures are preserved
  in `report/mock_validation.{md,json}`:
  - g-off physical slope: 0.9451 (target 0.95--1.05), a marginal 0.0049 miss;
  - coordinate-injection slope: 1.1262 (target 0.95--1.05), diagnosed as the
    remaining derivative mismatch of the smoothed measured-xi table;
  - truth, filtered-CMB, and matched-slab mean-field ensemble checks fail at
    this reduced footprint; the last two are about 2.3 sigma from zero;
  - six-bin shape p-value is 3.31e-8, showing the reduced mock/table does not
    preserve a common normalization across the stored bins;
  - mock-scatter/jackknife ratio is 12.28 (target 0.7--1.3), expected to be
    unreliable with artificial reduced-patch pair regions but still a failure.
  - the fiducial physical increment is 0.468 +/- 0.098 rather than unity;
  - response predictions were not computed independently, the combined result
    used a truth-template proxy, and individual flag/no-margin shifts were not
    implemented, so those rows are explicitly failed even where a proxy number
    happened to lie within a numerical tolerance.

These failures block Stage B under Section 2.8. The implementation must not be
treated as a completed Stage A calibration until the missing mock physics and
full `scale=1`, HEALPix-jackknife run pass.
