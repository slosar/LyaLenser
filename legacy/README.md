# Archived work

Two earlier threads of the project, kept for the record. **Nothing here is maintained or runnable as is**: the
code imports the old `code/pipeline` layout, which was refactored into the `lyalenser` package on 2026-09-25.
To run any of it, check out the last commit of the old layout (`git checkout 0a74cc1`) and follow its
`HANDOFF.md`. The mock products (32 GB) are under `$LYALENSER_DATA/mocks/` on the workstation and on NERSC.

| Directory | Content |
|---|---|
| `cmb_cross/` | The CMB-lensing cross-correlation thread: the first report (`report/main.tex`: literature review, the quadratic estimator in the sample-variance limit, the three-tracer separation, forecasts), its code (`three_tracer.py`, `cross_spectrum.py`, `recon_noise.py`, `make_plots.py`), the Mathematica derivation checks, the codex reviews, and `FUTURE_WORK.md` with the plan for resuming. |
| `mocks/` | Stage A: the lognormal mock generator (`pipeline/mock.py`), the provenance-hashed campaign machinery (`campaign4.py`, `campaign7.py`, `run_mock_validation.py`, `run_lowz_validation.py`), the acceptance protocols (`GATES.md`), the design documents (`PLAN.md`, `IMPLEMENTATION.md`), the per-iteration briefs and the full technical log (`NOTES.md`), the campaign reports (`report/mock_validation.md`, `report/lowz_validation.md`) and the HTCondor (BNL RACF) and Slurm (Perlmutter) drivers. |
| `notes/INSTRUCTIONS.md` | The original task statement. |

## What the CMB cross-correlation thread taught us

The aim was to detect lensing of the forest by cross-correlating the pair estimator (`lyalenser.pairs`,
`lyalenser.amplitude`, unchanged since) with the ACT DR6 and Planck PR4 convergence maps. Six mock iterations
(campaigns of up to 400 lognormal realisations of a 400 deg^2 patch, gated acceptance protocols, adversarial
reviews by codex/gpt-6-astra) led to the following conclusions.

1. **The estimator works.** The pair-template estimator recovers an injected lensing amplitude without bias,
   its jackknife errors are calibrated, and the linear response of the forest to long density modes (the
   dilation term that mimics lensing) is deprojected by a kernel-matched quasar template at the campaign's
   precision (`mocks/report/mock_validation.md`, iteration 6).
2. **A CMB template does not measure lensing alone.** On a mock with no lensing and no linear response the
   CMB cross-correlation still returns A = -0.41 +- 0.18. The term is the second-order forest-forest-tracer
   correlation <m_p m_q T> of the forest's own redshift range (m the long-wavelength forest component, T the
   realisation's non-Gaussian tracer template): it acts through the same long modes and separations as the
   lensing signal, so no filter removes it, and a scalar long-mode subtraction (`longmode_diagnosis.py`, 96 seeds)
   changes neither the precision nor the bias. With the modelling available, the CMB cross-correlation measures
   a combination of lensing and the intrinsic three-point function that one template cannot separate
   (review 6, `cmb_cross/report/reviews/codex_review_6.md`).
3. **Hence the pivot to low-redshift tracers.** A template built from objects at z < 1.6 shares no density mode
   with the forest at z > 1.96, so the contaminant vanishes by construction. That became the paper.
4. **If the CMB thread is resumed**, the way forward is a two-amplitude fit (A_lens, A_2) with a quadratic
   template Q_pq = sum_x K(x) xi_Fq(p, x) xi_Fq(q, x) whose shape is fixed by the measured forest-quasar
   cross-correlation; A_2 would be a detection of a new three-point statistic, A_lens an upper limit at DR1
   depth. The six-bin shape test (p = 0.002) shows the two shapes are distinguishable. Steps, the reviewer's list
   and an alternative (conditional pair-moment subtraction) are in `cmb_cross/FUTURE_WORK.md`.

Practical lessons that carried over to the low-redshift analysis:
- Injection tests check signs, cuts and the band bookkeeping; they do not test additive biases, because an
  injected displacement and its template are correlated by construction. Validation needs an independent
  handle (for the paper: the template x CMB comparison).
- Freeze the numerical choices before looking at the result and hash them into the products (the campaign
  provenance fingerprint); the paper's robustness rows are the data-side version of this.
- numpy < 2 for bit-identical numba pair kernels; NaMaster white levels overshoot on fragmented masks; the
  smoke-scale mocks need a patch large enough for the lowest band (L_fund below 100).
- From the first report: a quadratic estimator hardened against the forest's density response pays a factor
  of about 8 in noise (`cmb_cross/report/main.tex`, "Caveats"), which is why the density is deprojected with a
  template instead of hardened away.

## Reading order for the archive

`cmb_cross/report/main.pdf` (theory), `mocks/PLAN.md` and `mocks/IMPLEMENTATION.md` (design),
`mocks/GATES.md` (what was tested), `mocks/ITERATION7.md` and `mocks/NOTES.md` (what happened, in order),
`cmb_cross/FUTURE_WORK.md` (how to resume).
