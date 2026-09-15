# Stage A, iteration 7 brief: the low-redshift pivot (2026-09-15, implemented by the Claude session; no Codex review by user decision)

## Decision (user, 2026-09-15)
After review 6 (`report/reviews/codex_review_6.md`): with the modelling available the CMB cross-correlation
measures a combination of lensing and the intrinsic second-order forest-forest-density term that one template
cannot separate. The two-parameter (A_lens, A_2) constraint is recorded for later in `CMB_FUTURE_WORK.md`. This
project pivots to cross-correlations with low-redshift tracers, where that term vanishes; no detection is
expected at DR1 depth, upper limits are the product. Requirements:
- tracers ELG, QSO and LRG, in redshift slices so that bias evolution is handled per slice;
- the bias of each tracer from its angular auto-correlation on large scales, k < 0.2 h/Mpc, converted to an
  ell_max per slice (ell_max = 0.2 chi(z_mid): 263, 349, 442, 566, 648);
- A_L measured in each slice and combined optimally for the total detection;
- sightlines C^-1 weighted;
- the gate: an unbiased detection against a lognormal redshift tracer (mock).

## What was built (details in NOTES.md "Iteration 7")
- `lowz.py`: tracer table (8 tracers in 5 slices, DR1-like), Limber slice covariances, Gaussian slice maps,
  2-D lognormal Poisson tracers of the 3 Mpc/h-smoothed projected density, kernel-weighted templates with the
  forest-source lensing kernel, bias fits on 40 <= ell <= 0.2 chi(z_mid), model-covariance Wiener combination
  of the tracers of a slice, sum over slices, joint jackknife covariance and optimal combination.
- `mock.generate_mock(lowz=True)`: the foreground of kappa_lya (and kappa_CMB) is the sum of the slice
  contributions plus the uncovered rest; the low-z catalogue and randoms are saved with the mock.
- `run_lowz_validation.py`, `campaign7.py`, `slurm/campaign7_perlmutter.sh`: the campaign (fresh seeds
  4000-4399, dev 5000-5004, no dense mocks), GATES v7 rows, controls (tests, injection expectation, benchmark).
- `code/stageb/lowz_catalogues.py`: the same construction on the sphere from the DESI DR1 LSS catalogues (LRG,
  ELG_LOPnotqso, QSO z < 1.75; weights; randoms; footprint from the randoms; pseudo-C_ell bias fit; per-slice
  and combined alm for Stage B).
- Two bugs found by the template-chain tests and fixed before any campaign product: the half-pixel offset of
  the flat-sky catalogue binning (affected the iteration-6 matched quasar map by ~3 % at the top of the band)
  and the pixel window squared in the Wiener numerator. The 2-D lognormal needed smoothing (like the quasars).

## Verification path
1. Template chain at scale 1 without the forest: b recovered to 1-3 %, combined template normalisation
   0.989 +- 0.011 (three realisations).
2. Unit tests (`tests/test_iteration7.py` + the existing suite).
3. Smoke chain at scale 0.25 (scale 0.1 is unusable for any band template: the 2-degree patch's fundamental
   multipole is 180, so the 40-100 band is empty and the response matrix singular). RACF Condor
   (`condor/campaign7.sh`, pool idle after the upgrade) or Perlmutter (`slurm/campaign7_perlmutter.sh`).
4. Scale-1 campaign: RACF Condor (`condor/campaign7.sh`, 400 seed jobs of 8 cores / 32 GB) or Perlmutter preempt
   (~30 node-hours; NERSC access needs the user's OTP to renew the sshproxy certificate).
5. `report/lowz_validation.{md,json}`: the acceptance table; the gate is the combined recovery row.

## Stage B on DR1 after the gate
Build the sightline set from the DR1 deltas (spherical geometry; the pair machinery works on ra/dec already),
the per-slice alm from `lowz_catalogues.py`, deflection templates via `templates.alpha_at`, jackknife on the
DR1 footprint, injection expectation on the real geometry, random-template nulls; per-slice A_L and the combined
value with its error; bias cross-checked against the kappa_CMB x tracer cross-spectrum (ACT/Planck).
