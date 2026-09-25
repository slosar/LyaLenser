# Stage A, iteration 2: defects to fix before the acceptance run (reviewer notes, 2026-09-10)

The iteration-1 framework (modules, 12 unit tests) is kept. The mock and the validation runner do not yet implement
IMPLEMENTATION.md Sec. 2.7-2.8, and the acceptance failures are dominated by that. `/data/LyaLenser` is now writable
from the Codex sandbox (`~/.codex/config.toml` writable_roots); use `/data/LyaLenser/mocks/` as in the spec.

## Blocking defects (mock.py)
1. **Sightline sampling is wrong.** `prob = exp(3.5*proj/np.std(proj) - 0.5*3.5**2)` exponentiates a unit-variance
   field times the bias, i.e. a lognormal with ln-scatter 3.5, which concentrates the sightlines in a handful of
   cells. The spec (2.7 item 3) is a lognormal in the *physical* density: `n exp(b_q delta_g - b_q^2 var(delta_g)/2)`
   with `delta_g` the (unsmoothed, 2 Mpc/h) `delta_m` in the shell where the quasar sits, `b_q = 3.5`, no
   normalisation by the standard deviation. Sample quasars in 3D (position and chi) over the *template slab*
   (forest slab +- 150 Mpc/h) from that density, apply the linear-theory RSD shift of chi from the LOS velocity, then
   take the sightlines as the quasars in the forest slab thinned to `n_los` per deg^2. This is very likely the cause
   of the 0.47 +- 0.10 fiducial slope, the failed mean-field ensembles and the shape-test failure; re-measure them
   after the fix and report.
2. **Template quasars and randoms.** The matched template must be built from the clustered quasar sample of item 1
   (all quasars in the template slab, including the sightline quasars unless the disjoint-selection flag is on) and
   uniform randoms (20x), with the optional completeness pattern applied to both. Currently `qcat` is uniform, so the
   deprojection tests measure nothing.
3. **Outside-box convergences.** Replace the ad hoc `rest` field by the correlated Gaussian pair with the three
   restricted Limber spectra (`three_tracer.limber` with chi limits outside the box) as specified in 2.7 item 4:
   `C^{kl,rest kl,rest}`, `C^{kl,rest kc,rest}`, `C^{kc,rest kc,rest}`, generated on the patch by FFT; then
   `kappa_lya = kappa_lya_box + kappa_lya_rest`, `kappa_CMB = kappa_slab + kappa_rest`. Verify that the mock's
   `C^{kl kc}`, `C^{kl kl}`, `C^{kc kc}` measured on the patch agree with `three_tracer.spectra()` at the 10% level
   over 40 < L < 300 (add this as a check in the runner).
4. **Magnification flag** must modulate the quasar density by `(1 + m kappa_q)` with `kappa_q = kappa_lya` at the
   quasar positions, `m = 0.5`; **completeness flag** multiplies data and randoms by a smooth pattern `c(theta)`
   (e.g. 1 - 0.3 sin^2 structure); **real-mask flag** should use the actual ACT DR6 baseline mask
   (/data/LyaLenser/raw/act/baseline/mask_act_dr6_lensing_v1_healpix_nside_4096_baseline.fits; a flat-sky cut-out of
   a region of comparable size is sufficient), applied to `kappa_CMB` before filtering.
5. `g1` must be computed from the kernel-product average defined in IMPLEMENTATION.md Sec. 0 (not a single lens
   plane at z = 1); document the value.

## Blocking defects (runner and amplitude)
6. Science bands, junk band and curl partners must be real templates: band-pass filters on the patch (flat-sky FFT)
   for [40,100], [100,200], [200,300], junk = everything else up to the grid Nyquist, curl = 90-degree rotation of
   each; all fitted jointly through the response matrix. `amplitude()` must refuse a fit without the junk band.
7. Independent response prediction for tests 7(a) and 7(c): from `R_delta`, the mock's own `C^{delta_s kappa_CMB}`
   and `C^{delta_s delta_s}` (measured on the patch or from the spectra) and the catalogue's response to an amplitude
   modulation (response of the estimator to `delta_F -> (1+a) delta_F`, computed from the catalogue with a template
   built from the mock's `delta_L` map), computed before the fit.
8. Jackknife regions by pair **midpoint** in HEALPix pixels at `nside_jk = 8` on the full 20-degree patch (this gives
   ~10 regions; use `nside_jk = 16` if fewer than 30 regions result and say so); no artificial pair-index regions.
9. Mask-edge rejection of sightlines within 2 degrees of a mask edge when the real-mask flag is on.
10. Run the full configuration: `scale = 1`, >= 20 seeds, all acceptance items of Sec. 2.8 with their tolerances,
    outputs to `report/mock_validation.{md,json}` and figures. Keep the 15%-scale run as a smoke test only.

## Reporting
Update `NOTES.md` (keep the iteration-1 section, add an iteration-2 section), state the wall time of the full run and
the peak memory, and list any remaining failure with a diagnosis. Do not modify IMPLEMENTATION.md, PLAN.md,
report/main.tex or mathematica/. Do not commit.
