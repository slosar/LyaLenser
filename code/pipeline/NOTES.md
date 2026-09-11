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

## Iteration 2 (2026-09-10)

### Reviewer defects resolved

- Confirmed write access at the start of the iteration and wrote all final
  mocks/catalogues to `/data/LyaLenser/mocks/` (20 recovery mocks and pair
  catalogues, plus seed 20 for the mean-field null). The final products occupy
  about 3.8 GB.
- Quasars are now sampled in 3D from
  `exp(b_q delta_g - b_q^2 var(delta_g)/2)` using the unnormalised physical
  2 Mpc/h matter density, then shifted in chi by the linear LOS velocity.
  Sightlines are a thinning of forest-slab quasars; the matched catalogue uses
  the same clustered template-slab quasars, including sightline quasars, and
  uniform 20x randoms.
- The outside-box convergence pair is generated from the three restricted
  `three_tracer.limber` spectra. Across 20 scale-1 patches, the measured/theory
  ratios over `40 < L < 300` are 0.9801 (`kl-kl`), 0.9819 (`kl-kCMB`), and
  1.0279 (`kCMB-kCMB`).
- Magnification and completeness modify the data catalogue as specified;
  completeness also modifies randoms. The real-mask path reads the released
  ACT DR6 nside-4096 map, rotates an actual cutout onto the mock patch, masks
  CMB convergence, and rejects sightlines within 2 degrees of an edge at
  scale 1.
- `g1 = 2.60360917e-4 (Mpc/h)^-1` is computed from the `kl x kCMB`
  kernel-product lens distribution, whose mean lens distance is
  2007.010 Mpc/h. It is no longer a z=1 single-plane value.
- The fitted basis is three disjoint flat-sky science bands (40-100, 100-200,
  200-300), each science-band curl, and the outside-band junk template.
  `amplitude()` rejects any fit without a junk template. The response matrix
  fits all seven templates jointly.
- Response predictions are formed before the observed score fit from
  `R_delta`, the mock spectra, the stored `delta_L` map, and the catalogue's
  `xi*G` response to multiplicative forest modulation. The response-only CMB,
  matched, and deprojected checks all pass.
- Jackknife regions use pair-midpoint HEALPix pixels. Nside 8 produced fewer
  than 30 regions, so the full run used nside 16 as requested; individual
  mocks contain 26-27 occupied regions. No pair-index regions remain.

### Numerical and memory choices

- The 1 Mpc/h measured-correlation counts use a 0.76-bin normalized Gaussian
  smoother. The previous 0.45-bin value followed empty-cell derivative noise.
  With 0.76 bins, the data/model normalization ratio is 1.0074, the physical
  slope is 1.0069, and the independent coordinate-injection slope is 1.0472.
- The analytic Hankel transform uses 1600 log-k samples in validation. The
  iteration-1 360-sample call rang strongly beyond 10 Mpc/h and gave an
  invalid model normalization.
- The matter, forest, and lensing fields retain the fixed 2 x 2 x 0.5 Mpc/h
  grid. To remain below the host memory ceiling, the linear velocity solve is
  performed on every fourth cell of the same matter realization; velocities
  are large-scale dominated. FFTs use four workers while all Numba pair
  kernels use the required 24 threads. Sequential FFT products and blockwise
  response multiplication reduced peak resident memory from a killed 30.6 GB
  attempt to 19.474 GB.
- The fixed forest field factor 0.903 is the discrete-grid/pixel-window
  normalization relative to the continuum Kaiser transform; the matter field
  and quasar lognormal density are never normalized per realization.

### Runs and acceptance

- First scale-0.15 smoke: 146.084 s, 1.326 GB peak; it exposed the positional
  band-argument bug. Corrected scale-0.15 smoke: 130.410 s, 1.339 GB peak;
  all deterministic normalization/injection checks passed. Ensemble rows on
  the two-seed 3-degree patch were treated only as smoke diagnostics.
- Near-full scale-0.999, two-seed pilot: 361.505 s and 18.752 GB peak. It
  verified the full band/mask geometry before the acceptance launch.
- Scale-1 acceptance, seeds 0-19: normalization mock 6.076 s, analytic xi
  0.997 s, measured xi 26.538 s, 20-seed ensemble 2508.150 s, flag/no-margin
  shifts 431.561 s. The fixed extra seed-20 mean-field null took 133.780 s.
  Total internal acceptance campaign time was 3217.755 s; the main scale-1
  process alone was 3083.975 s. Peak resident memory was 19.474 GB.
- All 11 Section 2.8 rows pass. Combined recovery is 1.380 +/- 0.789, shape
  p=0.783, and scatter/RMS-jackknife is 0.865. The switch-off shifts relative
  to all flags on are -0.089 (magnification), +0.582 (completeness), +0.919
  (real mask), and +0.202 (no template margin).
- The mean-field matched-template result with 20 seeds was a marginal 2.044
  SEM fluctuation. A predeclared additional null seed gave the 21-seed result
  40.897 +/- 21.573 (1.896 SEM), while truth and CMB are also within 2 SEM.
  The original covariance summary used the median of strongly heteroscedastic
  jackknife errors; comparing variances requires RMS error. The corrected
  scatter/RMS-jackknife ratio is 0.865 and passes.

### Remaining failures and diagnoses

- None. Stage A acceptance is 11/11 PASS. The large matched-template null mean
  and its large uncertainty are retained in the report rather than hidden;
  the subtraction is demonstrated by comparison with the raw 20-seed mean.

## Iteration 3 (2026-09-10/11; supersedes iteration-2 acceptance claims)

The iteration-2 “11/11 PASS” claim is invalid. Its paired differences removed
additive errors from recovery, its taper basis was incomplete, and its baseline
was reused from a different sample. The historical text above is retained as
provenance; it is not evidence of Stage-A acceptance.

### Implementation and prospective protocol

- The public amplitude API now requires all three science bands, all three
  curl partners, and complementary junk. It rejects singular full response
  matrices. The old small algebra fixtures use the internal contraction fit;
  their q/F/mf identities remain tested. The old zero-junk fixture was
  singular, and the old “junk orthogonal to taper wings” assertion contradicted
  completeness. Those two fixtures were corrected, with explicit production
  rejection and unit-unfiltered-signal tests added.
- Pair sign, eleven accumulator formulas, moment contractions, spherical
  derivative convention, positive physical remapping, negative coordinate
  injection, midpoint HEALPix jackknife, and restricted outside-box Gaussian
  spectra are retained. Pixel validity and midpoint slab selection surround
  the existing accumulator arithmetic; `accumulate_slabs` persists one group
  per slab. Both xi and compression kernels consume pixel labels.
- The continuum operator is a joint weighted intercept/slope projection;
  its weighted coordinate is centered before fitting. Periodic trilinear
  interpolation wraps the last cell to the first. RSD uses the original
  Fourier modes, with no density decimation, and is interpolated at objects.
  Counts are Poisson draws from integrated realization-dependent intensity.
- All forest and projected-density/convergence samples use the same physical
  Cartesian light-cone rays. Local angular displacements use cos(dec) and
  preserve the pair kernel's cos(dec_mid) metric. The spherical-ray construction
  is needed because a globally Euclidean chart cannot have the spherical
  cos(dec) metric throughout a twenty-degree patch. The flat Fourier map
  potential remains an approximation whose amplitude error must be measured;
  the local metric test does not prove global statistical normalization.
- The empirical 0.903 forest factor is removed. Production xi is remeasured
  for each lensed/noisy/response variant. Development seeds 100–104 choose the
  smoothing width before validation, and fixed-versus-refitted baselines are
  compared on each validation realization. Analytic integration uses nk=6400,
  checked against 12800 over the response range. Tables now save their coarse
  numerator/denominator counts as well as the interpolated xi and derivative.
- Matched shot noise is integrated as (1+r) integral W²/(b² nbar_chi) dchi,
  where nbar_chi=nbar_3D chi². Completeness is normalized over the independently
  specified footprint, including empty random pixels. The spherical interface
  requires that footprint rather than inferring it from occupied pixels.
  CMB/matched maps share pixels and masks; CXX comes from the actual noisy map,
  S includes the actual selection and magnification and common mask transfer.
- The new independent response code forms P_a[D_a C_ab D_b−C_ab]P_b^T from
  every stored pixel's smoothed long-density value. Its scores do not use
  observed forest products, and the deprojected prediction is computed. The
  raw covariance is derived from the mock discrete power spectrum with the
  exact phase-averaged trilinear window (product of cubic B3 functions); the
  corresponding covariance cube and a 128-versus-256-angle convergence check
  are saved. Transverse-angle averaging and finite radial support remain
  approximations, to be assessed by the absolute prediction gate. The old sightline-level spectral helper
  remains solely for its historical algebra regression test and is not used
  by acceptance.
- All realizations use the largest 300 Mpc/h box. The validation mock retains
  the largest quasar/random catalogue; matched margins 0/150/300 select nested
  subsets while keeping the forest, sightlines, and CMB realization fixed.
- `GATES.md` predates the validation run. Provenance includes its hash, source
  hashes, frozen config, UTC launch, development results, and the prospective
  seed/extension rule. No git commit is made, as requested. Initial recovery
  seeds are 0–19 and null seeds 20–39; all 40–59 are added only if the declared
  residual precision rule requires it. Fixed-template nulls use development
  seed-100 maps/operators on every null realization. Shape uses the full
  covariance Gaussian Hotelling test. Injection uses the full seven-component
  basis. All recovery/null amplitudes are absolute, with SEM and achieved
  residual bounds; paired differences are only labelled diagnostics.
- HDF5 saves each variant's pair catalogue, xi/coarse counts, per-region
  q/F/mf, full fit matrices, templates/potentials, map operators, random maps,
  and response prediction catalogue. Per-seed JSON contains absolute and raw
  fits, predictions, shape values, timing, memory, and configuration. The
  `--phase rebuild` command reconstructs the report and figures from those
  saved statistics without generating any new mock.

### Development observations and remaining limitations

The first full development seed gave slopes below unity for widths 0.76,
0.9, and 1.1. The development scan was therefore broadened to include 1.3,
1.5, 1.7, and 2.0 before validation; its declared criterion remains the closest
mean slope to unity on development seeds only. This is not validation-seed
calibration. A weak-model-style choice based on fresh recovery seeds is
explicitly forbidden.

The A_true={0,1,5,10} requirement is followed. Physical remapping at large A
need not remain linear when displacements are comparable to correlation
lengths; a failed slope can therefore reflect finite remapping effects as
well as a derivative/baseline error. The requested linear regression is
reported regardless, with the individual absolute amplitude points retained
so that this distinction can be investigated. No claim that the precision
rule by itself proves linearity is made.

The absolute residual bound requirement may be unattainable at the fixed
maximum ensemble size. Such a precision failure must remain FAIL and block
Stage B, even when a noisy mean happens to be within two SEM of its target.

### Reproduction

Use `/home/anze/anaconda3/bin/python3` with `NUMBA_NUM_THREADS=24`,
`OPENBLAS_NUM_THREADS=4`, `OMP_NUM_THREADS=4`, and a writable
`MPLCONFIGDIR=/tmp/mpl-iter3`. From `code/pipeline/`:

```
python run_mock_validation.py --phase development
python run_mock_validation.py --phase freeze
python run_mock_validation.py --phase validation
python run_mock_validation.py --phase rebuild
```

The default large-product directory is `/data/LyaLenser/mocks/iteration3/`.
`validation` runs a pytest preflight, writes frozen provenance, and refuses to
resume after source, gate, or configuration changes. Completed seed JSONs
are resumable only under that identical provenance. `rebuild` performs no
new simulation. Logs for this campaign are stored alongside the products.

An intermediate pre-freeze regression run passed 34 tests in 130.29 seconds
while sharing the machine with development jobs; the final preflight passed
36 tests, as recorded below. The last observed initial-development peak was
28.85 GiB (30.98 GB); its terminal peak is unknown after the interruption.
Final campaign times, gate failures, and peak memory are recorded below.

### Development interruption and allocation fix

The initial development process exited with SIGKILL (137) before completing
seed 104. Seeds 100–103 and their diagnostics were intact; validation had
not started. The last instrumented peak was 28.85 GiB, but the process's final
peak is unavailable, so that number is not asserted to be its termination
peak. Memory pressure is plausible, not proven by the available evidence.
The interruption is recorded in `development_interruption.json`.

The resumed generator removes redundant float32 copies after float32 FFT
inverses and computes the lognormal density variance in radial blocks instead
of allocating a full float64 volume. These changes preserve the numerical
operation up to summation roundoff and reduce avoidable memory. The radial
projection batching described above is also used by the resumed seed and by
validation. Completed development mocks are retained, and only the unfinished
seed is regenerated. The two development jobs overlap in time; their recorded
wall durations must not be added to claim an elapsed campaign time.

### Frozen development result and preflight

The final width is 2.0 coarse (1 Mpc/h) bins, chosen from development seeds
100–104 only. Their mean physical slope is 0.97477250, SEM 0.06109119,
with a 95% residual bound of 0.19484384 in slope units. Individual slopes
are 0.98208106, 1.13693047, 1.04798743, 0.93554898, and 0.77131455.
The absolute-A fractional changes on halving xi spacing to 0.125 Mpc/h
are 0.00010496, 0.00026960, 0.00204271, 0.00259461, and 0.00072276.
The discrete-grid covariance derivative changes by 7.872e-7 in norm on
doubling angular samples from 128 to 256. These are development/numerical
checks, not substitutes for the fresh normalization and residual gates.

The authoritative final preflight passed **36 tests, 0 failures, in 50.904 s**
before the validation launch. This supersedes the intermediate 34-test timing
above. The original baseline had 17 passing tests. Nineteen added regression
tests cover weighted continuum moments/idempotency, periodic endpoints, actual
FFT deflections, complete-band unit recovery, rays and cross-field projections,
full-mode RSD, integrated Poisson counts, matched shot noise and flat/spherical
nonzero response, both slab kernels and disjoint catalogue selection,
amplitude/nk convergence, independent projected pixel response, absolute
statistics, full-covariance Hotelling calibration, and the exact trilinear
covariance window. The existing q/F/mf algebra assertions remain intact.

Provenance was written at 2026-09-11T03:37:16.619495+00:00, before seed 0.
The gate SHA256 is
`d3849680a1f088083beac26f9ab6fc21390322413d1d5a20e2ed89a86375c7a8`;
the configuration SHA256 is
`a22cf484f5831c75bb5bf3ce79c7b3c3827fa842fa00462ce5a82809c8bf7e37`.
The successful resumed development generator peaked at 23.953 GiB.
The development benchmark measured 8.031e8 pixel pairs/s on 24 threads
(3.346e7 per thread), extrapolating to 4.650 s for the sparse full catalogue.
The resumed initial scan took 297.109 s and the expanded scan 2423.861 s;
their overlap and the earlier interrupted attempt are retained in the
development timing records rather than hidden in a single additive total.

A disk-only audit reconstructs all 140 seed-0 stored full fits from their
per-region q/F/mf sums and verifies `solve(F, q-mf)` against the stored
amplitude vector. Its results are in `persistence_audit.json` alongside the
mocks. `report/render_iteration3.py` is a reporting-only second step after
`--phase rebuild`: it escapes table delimiters, labels slope bounds correctly,
and displays the saved raw/null and prediction statistics. It does not change
the frozen scientific source, numerical results, seeds, or gate decisions.

### Additional issue discovered during the frozen run

A read-only audit found that `_fft_fields` and `generate_mock` each initialize
`default_rng(seed)`. With CMB noise enabled, its first normal draw repeats the
prefix of the white-field normal stream. A small exact-prefix reproduction is
saved in `random_stream_audit.json`. Fourier filtering and the different array
shapes complicate the resulting map covariance, but they do not establish
statistical independence. This is an additional generator limitation beyond
the review's 17 findings. Its contribution to the observed amplitude offsets
has not been isolated, so it is not asserted to explain those offsets.

This was discovered after the freeze and after recovery seeds were inspected.
The validated scientific source and random draws are retained exactly; no
post-unblinding stream change is mixed into the ensemble. Independent child
random streams require a new development/freeze/validation cycle. The current
stochastic results cannot certify independent instrumental noise, regardless
of individual gate decisions. This issue remains open for the next iteration.

### Completed ensemble: absolute results and diagnoses

All seeds 0–59 completed under the original source/configuration/gate freeze.
After the first forty seeds, the combined residual bound was 5.416966 A and
the deprojected response-only null bound was 7.787777 A. Both exceed 0.3 A,
so the declared rule required all seeds 40–59. These twenty extension seeds
participate in both ensembles: recovery has 0–19 and 40–59; nulls have 20–59.
There are sixty distinct realizations, forty entries in each final ensemble,
and no seed above 59. The ensemble itself took **14711.106 s (4.086 h)**.

The disk audit reconstructs **2203 full fits** from per-region q/F/mf sums,
and independently reconstructs each observed raw and corrected common-science
scalar from its saved matrices. Every check passes. Scientific source, gate,
and preflight test-source hashes remain unchanged. This is recorded in
`final_persistence_audit.json`. It establishes numerical persistence, not
physical acceptance.

Seven of the eighteen ensemble gates fail:

1. **Filtered CMB normalization:** slope 1.045821 ± 0.066089 SEM; achieved
   95% slope residual bound 0.179499. The mean lies outside the frozen 0.97–1.03
   interval. The result is statistically consistent with unity, so it does
   not establish a 4.6% physical bias; it fails the requested tight calibration
   check. Finite-patch map transfer, masking/pixelization, baseline modelling,
   and the random-stream issue remain possible contributors that this run
   does not separate.
2. **Matched normalization:** slope 1.474798 ± 0.862855 SEM; slope bound
   2.220087. The mean fails the same interval and the calibration is very
   imprecise. Matched-template fluctuations are large even with A_true=10.
   No large multiplicative bias is claimed as statistically established;
   normalization at the required tolerance remains unverified.
3. **Varying matched mean field:** absolute null 25.611971 ± 8.878652 SEM,
   bound 43.570741 A, or 2.885 SEM from zero. The fixed-template matched null
   is −8.656895 ± 9.658276 SEM and passes. This pattern supports a residual
   connected to the template/data/selection correlation, rather than a
   demonstrated universal fixed-map offset. A single geometry-averaged xi
   subtraction does not automatically describe covariance conditional on
   clustered sightline selection and the matched density map. This is a
   diagnosis to investigate, not an isolated causal proof.
4. **Matched response prediction:** observed absolute A=52.154354 ± 9.430845
   SEM; prediction 33.155836 ± 2.093237 SEM; observed minus prediction
   18.998518 ± 8.760793 SEM, bound 36.718895 A. The discrepancy exceeds two
   SEM. Much of its scale is already present in the response-off matched
   null. The additional paired response-change diagnostic minus prediction
   is −6.613454 ± 3.412972 SEM; that diagnostic cannot replace the failed
   absolute gate.
5. **Deprojected response-only null:** absolute A=−2.966013 ± 0.956350 SEM,
   bound 4.900413 A. Observed minus predicted is −1.549463 ± 0.826059 SEM and
   meets the two-SEM prediction-agreement subtest, but the absolute bias and
   precision requirement fails. The independent prediction is −1.416550 ±
   0.266485 SEM; setting it to zero would conceal predicted leakage. The
   response-off deprojected null is already −2.347243 ± 0.844004 SEM. Both
   background subtraction and the response/template operator require further
   work before a sub-0.3 A residual can be claimed.
6. **Combined deprojected recovery:** absolute A=−1.633923 ± 0.713170 SEM;
   residual from unity −2.633923, bound 4.076445 A. Both the two-SEM consistency
   and precision subtests fail. Its paired response-on minus response-off
   diagnostic is −0.016886 ± 0.269019 SEM (bound 0.561028 A). This demonstrates
   why a small paired change cannot establish absolute recovery: the offset
   is retained in the actual estimator.
7. **Covariance:** absolute scatter 4.510481, RMS midpoint-jackknife error
   6.481210, ratio **0.6959319601**, below 0.7. This is a marginal numerical
   threshold failure, retained exactly. On this finite ensemble the jackknife
   errors are larger than the observed scatter. Sampling uncertainty,
   heterogeneous weights and finite jackknife regions have not been separated;
   the verified jackknife grouping is unchanged and is not retuned to pass.

The truth slope passes the frozen mean-slope gate at 0.999183 ± 0.027373 SEM,
with residual bound 0.056185 in slope units. A PASS on the predeclared point
tolerance is not a claim of 3% calibration at 95% confidence. The brief's
slope tolerance is applied to the ensemble mean, with its SEM and achieved
bound explicitly reported. The fixed-minus-refitted slope comparison passes:
refitted minus fixed is −0.012621 ± 0.055523 SEM. The diagnostic A_true=1
refitted-minus-fixed shift is +0.153295 ± 0.146080 SEM, bound 0.448770 A.
The Hotelling shape test passes (p=0.222214); total spectrum ratios are
0.987398, 0.988494, and 1.018703 and pass. These checks do not remove the
absolute failures above.

“Independent response prediction” means the covariance perturbation is built
from the grid power and stored pixel modulation, without using the observed
response score to calibrate it. Its estimator weights still use the measured
xi table, as does the observed fit. It is not statistically independent of
every datum and is not a full conditional simulation of baseline refitting.
The projected covariance uses transverse-angle averaging, finite 40 Mpc/h
support and an unconditional raw covariance; long/short-field conditional
correlations and changes in the empirically refitted baseline can matter.
For example, the deprojected paired response-change minus prediction is
+0.797780 ± 0.335977 SEM. Since the paired fits have separately refitted
baselines, this discrepancy does not isolate the covariance approximation
alone. These limitations require a fresh development cycle, not changes
selected using the present validation results.

The fresh full-scale dense interpolation check passes: A(0.25)=0.8772958420
and A(0.125)=0.8772376303, giving a fractional change of 6.6358e-5
(0.006636%, below 0.5%). These are absolute estimates on the same dense
A_true=1 sample. The dense A_true=0 catalogue has 1,397,463 sightline pairs;
its compression took 192.758 s. This explains the longer dense-control wall
times relative to the masked sparse ensemble. The saved matrices reproduce
the convergence comparison in `dense_convergence_audit.json`.

The completed dense g-on control adds two failures. Its absolute amplitude
points are [0.07463384, 0.87729584, 5.11426471, 10.69539373] at true amplitudes
[0,1,5,10], giving slope **1.07026136** with descriptive regression error
0.02310189. This exceeds the frozen 1±0.05 interval. One dense realization
does not supply an ensemble SEM or confidence bound. The points also show
that a single linear slope does not describe all local responses equally:
the diagnostic 0-to-1 increment is 0.802662, while the diagnostic 5-to-10
slope is 1.116226. Large physical remapping and independently refitted
baselines can change the operator as A_true changes; these data do not
isolate their individual contributions. This supports the earlier caveat
about treating A_true=10 physical remapping as exactly linear.

The dense analytic-baseline slope is **0.89665477** (descriptive regression
error 0.02410645), giving measured/analytic ratio **1.19361586**, outside the
3% tolerance. The analytic provider describes the continuum power model;
the simulated data also undergo finite-grid interpolation and continuum
projection. Those operations are included explicitly in the separate raw
response-covariance construction, but not in this analytic baseline provider.
Dropping the empirical factor does not make these providers equivalent. The
tiny interpolation-spacing change above rules out that spacing error as an
explanation of a 19.4% ratio discrepancy. This normalization comparison
remains FAIL and needs an analytic provider with the actual observation
operators or a justified alternative baseline comparison.

The **first-moment gate passes**: omitted/correct slope ratio 1.02085974,
catalogue prediction 1.01811120, difference 0.00274854. The scalar prediction
is reconstructed from its persisted full fit. Thus the observed relative
moment correction agrees with the verified contractions even though the
absolute dense normalization fails. The g-on audit is saved in
`dense_g_on_audit.json`.

Before the injection control ran, an additional reporting audit was specified
for its curl slope. The frozen implementation compares the slope against an
RMS of individual band-amplitude jackknife errors; that is an error scale,
not the joint uncertainty of a fitted slope. The audit fits the absolute
mean-curl values at [-2,-1,0,1,2] for each common leave-region sample, using
the union of the saved midpoint-HEALPix region IDs (an empty region leaves
that amplitude fit unchanged). It applies the unchanged jackknife variance
formula to those slopes and reports the resulting one-sigma check alongside
the frozen RMS check. It neither changes the injected samples nor retunes a
tolerance. The code is `report/audit_injection_jackknife.py`; its saved result
is `injection_joint_jackknife_audit.json` with the large products.

The **dense g-off/moments-off normalization also fails**. Its absolute points
are [0.02164386, 0.81340726, 5.07517767, 10.64822291], giving slope
1.07157610 with descriptive regression error 0.02298813. It differs from the
g-on slope by only 0.001315. The normalization discrepancy therefore persists
when both the radial lensing variation and first-moment contraction are
removed; it is not explained solely by those moment terms. Baseline and
physical-remapping nonlinearity remain the relevant unresolved comparisons.
The result is saved in `dense_g_off_audit.json`.

Coordinate injection **passes** with an odd/OLS slope of 0.99195702. Its
absolute fitted amplitudes at [-2,-1,0,1,2] are
[-1.65177963, -0.74718032, 0.38581018, 1.31041766, 2.27920646].
The curl slope is −0.09346845. The frozen RMS individual-band jackknife
scale is 3.02844638; the properly propagated joint slope jackknife error is
0.16785879. The curl check passes with both, so the weaker frozen error proxy
does not change this control's verdict. Both error definitions are reported;
the RMS scale must not be described as the slope's uncertainty. The joint
audit was registered at 2026-09-11T08:36:21.257695+00:00 before injection
products existed, with its code SHA256 recorded. Run
`python report/audit_injection_jackknife.py` before the report formatter to
reproduce the audit from the saved leave-region fits.

### Final campaign totals and diagnostic limits

The completed frozen table is **20/30 PASS, 10/30 FAIL**. The ten failed
gates are the seven ensemble failures above, dense g-on normalization,
dense g-off normalization, and data-versus-analytic normalization. Every
failure remains a failure; Stage B is blocked. Numerical derivative
convergence passes with relative norm change 2.76525e-5 on doubling nk from
6400 to 12800. The independent grid-covariance angular check also passes.

| Timing or resource | Measured value |
|---|---:|
| Final unit/regression preflight | 36 passed, 0 failed; 50.904 s |
| Sixty-seed ensemble | 14711.106 s (4 h 5 min 11 s) |
| Dense controls and remaining diagnostics | 5040.360 s (1 h 24 min 0 s) |
| Complete validation invocation, including preflight/setup | 19823.411 s (5 h 30 min 23 s) |
| Successful validation peak RSS | 25.856827 GiB = 27.763556 GB |
| Numba threads | 24 |
| Measured throughput | 8.12452e8 pixel pairs/s |
| Throughput per thread | 3.38522e7 pixel pairs/s |
| Sparse full-catalogue extrapolation | 4.789 s |

The peak above includes the dense and flag controls and is below 40 decimal
GB. It does not supply the unavailable terminal peak of the earlier killed
development process. The previously recorded 297.109 s resumed development
scan and 2423.861 s expanded scan overlap; neither their sum nor those values
alone describe the entire interrupted development elapsed time.

The flag diagnostics use seed 0. Relative to the all-flags-on absolute
A=−5.860876, magnification off gives −17.299665 (shift −11.438789),
completeness off gives −17.698383 (shift −11.837507), and ACT mask off gives
−1.186928 (shift +4.673947). Each flag can change catalogue selection and
downstream sampling. One realization does not determine an ensemble effect
size or its SEM. The nested margin selections give absolute A values
−5.320339, −5.860876, and −9.315305 for 0, 150, and 300 Mpc/h. Their shifts
from 150 are +0.540537, 0, and −3.454429. The larger underlying realization,
forest, sightlines, and CMB draw are shared throughout the margin test.

The 100 identically masked independent Gaussian templates give absolute
mean −0.599917 ± 0.391589 SEM, with achieved bound 1.376914 A. Their
scatter/RMS-jackknife ratio is 0.957035; scatter/RMS-sigma_F is 5.128539.
The naive inverse-response error is therefore much smaller than the observed
random-template scatter. The near-unit random-template jackknife ratio does
not erase the failed 0.695932 ratio for the physical combined ensemble:
the conditional template experiment and the realization-varying physical
ensemble have different covariance conditions. No error rescaling is fitted
to either result.

All **26 control fits** also reconstruct from their saved region sums, in
addition to the 2203 ensemble fits (2229 total). This audit is saved in
`control_persistence_audit.json`. The reports include both persistence audits,
the random-stream audit, and the independently registered joint curl audit.
The source/configuration/gate freeze was maintained through the complete run.

The review's numerical/core findings are accepted. The explicit departures
or interpretation limits are documented above: complementary junk requires
correcting the incompatible old orthogonality/zero-junk fixtures; spherical
rays preserve the requested varying cos(dec) metric more faithfully than one
globally Euclidean chart; large physical A need not be exactly linear; and
the frozen mean-slope tolerance does not itself imply a 95% precision claim.
The user explicitly prohibited commits, so GATES.md was written in the
working tree and hashed before validation without a git commit. The additional
random-stream independence issue and the loose original curl-error proxy
were disclosed; the latter passes the stronger joint-error audit as well.

### Final delivery verification

The explicit disk-only rebuild reproduced all 30 acceptance rows exactly.
The final JSON passes strict parsing without NaN or Infinity, and the primary
Markdown acceptance table is structurally valid. All four rebuilt PDF figures
were rendered and visually inspected: no clipping, overlap, or layout defects
were found. `final_delivery_audit.json` in the large-product directory records
these checks, PDF hashes, and unchanged frozen pipeline, tests, and GATES
hashes. The protected specification, plan, main TeX report, and Mathematica
files have no diff. No commit was made.
