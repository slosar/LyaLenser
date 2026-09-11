# Progress log

## Status summary
- [x] Clarifying questions answered, INSTRUCTIONS.md updated (2026-09-10)
- [x] Repo scaffolding (CLAUDE.md, PROGRESS.md, MEMORY.md, dirs, git, GitHub)
- [x] Step 1: literature review (report Sec. 1, report/literature_notes.md, main.bib via adstex, 32 refs)
- [x] Step 2: sample-variance-limited estimator formalism (report Sec. 2; Mathematica checks A, B, C pass)
- [x] Step 2 numerics: N_kappa(L), bias-hardened N, S/N (code/make_plots.py -> report/figures, numbers.json)
- [x] Report compiles cleanly (latexmk, 14 pages)
- [ ] Codex adversarial reviews run (see pending list)
- [x] Step 3: three-tracer (quasars, forest, CMB kappa) separation of lensing from response at the power/bispectrum level, plus realistic noise (report Sec. 4, code/three_tracer.py, mathematica/three_tracer_checks.wls)
- [ ] Follow-ups listed under "Open items"

## Decisions (from clarifying questions, 2026-09-10)
- Codex MCP down at start: proceed, queue reviews here.
- Estimator: flat sky, single source plane, 3D Gaussian delta_F with redshift-space P_F.
- Derive for discrete sightlines/pixels first, then dense-sightline continuum limit.
- Numerics: toy N_kappa(L) and S/N with fiducial P_F and DESI-like footprint.
- Report: single report/main.tex, adstex bibliography.
- Python: anaconda unless a package is missing.
- Git: local + private GitHub repo slosar/LyaLenser.
- Data: none needed yet.

## Next session
User intends to restart and have codex (gpt-6-astra high) review steps 1-3 and PLAN.md. Review checklist is in PLAN.md Sec. 6 and below. After the review: start PLAN.md phase 1, module 0 (data fetch to /data/LyaLenser).

## Codex reviews
- Review 1 (2026-09-10, gpt-6-astra high): report/reviews/codex_review_1.md. 18 items; continuum normalisation confirmed; blockers on broad-slab deprojection (fixed: kernel-matched template, report Eq. matched) and on the real-space template exactness (fixed: pair-template estimator, report Sec. 5). All items addressed or carried forward in PLAN.md Sec. 5 and report App. B.
- Review 2 (2026-09-10, gpt-6-astra high, estimator only): report/reviews/codex_review_2.md. 16 items; pair response and compression confirmed; blockers on the amplitude normalisation (fixed: conditional template coefficient with h_L = S_L/C^XX, band response matrix) and on injection calibration (reclassified as bookkeeping test; physical normalisation via measured xi validated on joint mocks). All items folded into report Sec. 5, PLAN.md Sec. 0, IMPLEMENTATION.md.
- Review 3 (2026-09-10, gpt-6-astra high, confirmation pass): report/reviews/codex_review_3.md. Two blockers in IMPLEMENTATION.md (pair-direction sign reversed relative to the report; first-moment normalisation linear instead of quadratic, missing cross term) plus 9 significant items; all fixed in IMPLEMENTATION.md v3, report Sec. 5, PLAN.md (Mathematica checks 11-12 added). Verdict before fixes: not ready; after fixes the Stage A acceptance tests (11 items with tolerances) are the gate.
- Implementation round 1 (2026-09-10, codex gpt-5.6-sol): framework + 12 passing unit tests delivered; acceptance run invalid (sandbox could not write /data; 15%-scale proxy; mock physics incomplete; sightline sampling bug). Defect list: code/pipeline/ITERATION2.md.
- Implementation round 2 (launched 2026-09-10, gpt-5.6-sol, /data now writable): fix ITERATION2.md items, full scale=1 x 20-seed acceptance run.
- **Decision (user, 2026-09-10): if the code is still defective after round 2, round 3 is done with gpt-6-astra.** Afterwards: gpt-6-astra review of code + mock validation before Stage B.

## Pending codex reviews (historical list, all covered by review 1)
1. Derivation of the mode coupling, squeezed limit and amplitude degeneracy (report Sec. 2.2-2.3, mathematica/*.wls).
2. Discrete optimal quadratic estimator and cubic cross-correlation estimator, list of biases (Sec. 2.4-2.5).
3. Continuum-limit N_kappa(L) formula incl. factors of 2 and mode counting (Sec. 2.6, code/recon_noise.py).
4. Bias-hardening formulae and the claim that hardened noise is set by transverse resolution (Sec. 2.7, Sec. 3).
5. Toy numbers: forest power model, CMB white-noise calibration, S/N (code/make_plots.py, Table 1).
6. Step 3: deprojection identity and its assumptions (Sec. 4.1-4.2), noise model N=(P_N+P_1D)/n_eff with unchanged response (Sec. 4.3), and the realistic-noise S/N table (Table 2).

## Open items (candidates for the "subsequent steps")
- Replace white CMB noise by public ACT DR6 / Planck PR4 N_L^{kk} curves.
- Refine Arinyo-i-Prats parameters to published z=2.4 values (currently approximate).
- Add pixel noise to P in N_kappa; add realistic forest n(z); drop single-source-plane approximation.
- Compute the intrinsic clustering bias (response of P_F to long modes x C_L^{delta kappa_cmb}) and the N1-type term.
- Design hardening against the tidal response as well as the density response.
- Discrete-sightline (aliasing) treatment via the pixel-covariance estimator on a mock.

## Log
- 2026-09-10: session start; asked clarifying questions; scaffolded repo; GitHub repo slosar/LyaLenser created.
- 2026-09-10: literature search (subagent, ADS-verified) -> report/literature_notes.md. Key finding: only Croft/Metcalf/Shaw group (5 papers 2018-2025); no CMB-lensing cross-correlation of forest lensing exists in the literature.
- 2026-09-10: derivation written (report Sec. 2), Mathematica checks pass, numerics run. Main result: for DESI-like 50 sightlines/deg^2, sample-variance N_kappa(100) ~ 2.5e-8 ~ C_L^{kk}; amplitude-hardened noise ~8x larger and independent of k_par,max. Toy S/N for cross-correlation with ACT-like kappa: ~37 (naive), ~15 (hardened); these are upper limits (no pixel noise, no aliasing, Gaussian).
- 2026-09-10: user question after step 2: does the response contamination force a switch to low-z galaxies (Shaw+2025)? Answer given: not yet; contamination can be removed with a high-z template. Step 3 defined: three-tracer separation at the spectrum/bispectrum level + realistic noise.
- 2026-09-10: step 3 done. Results (DR1-like: n_eff=25/deg^2, P_N=0.33, slab 2.1<z<3.0): response contamination = 3x lensing signal at L=100 (R_delta=2); deprojection with quasar template (beta = C^{d kc}/C^{d q} = 0.084, geometric) removes it at ~15% S/N cost; sampling (aliasing) noise dominates N_kappa; deprojected S/N ~2.1 (ACT or Planck), 2.8 with cleaned spectra, 5.1 for n_eff=50 + clean spectra. Consistent with Shaw+2025.
- 2026-09-10: user asked about C^-1 weighting per quasar. Added harmonic-mean noise model (Eq. harmonic in report): with a log-normal per-quasar noise distribution (median P_N=0.33, sigma_ln=1 or 2) the deprojected S/N for DR1-like goes 2.1 -> 2.2 / 2.3 (ACT); at n_eff=50: 3.8 -> 4.1 / 4.4. Gain limited by the aliasing floor P_1D/n_eff. Action item: use the actual DR1 per-pixel noise distribution.
- 2026-09-10: wrote PLAN.md (real-space cell-based estimator with E/B forest shear, quasar deprojection, ACT+Planck, self-lensed mocks). Resource assessment: runs on this machine; no GPU or NERSC needed for phases 1-2.
- 2026-09-10 (new session): codex plugin available (codex-cli 0.154.0, ChatGPT login). Launched adversarial review (gpt-6-astra, high) of report Secs. 2-4, code/, mathematica/, PLAN.md Sec. 6 checklist. Awaiting result.
- 2026-09-10: review 1 incorporated. Code: Limber to recombination (C_cc up 15% at L=100, CMB white noise recalibrated), response contamination factor 1/2 fixed (now 1.6x signal at L=100), beta with magnification, response sample variance, n_los=22, kernel-matched template, S/N over 40<=L<=300. Corrected DR1 forecast: matched-template S/N 1.4 (40-300) / 1.7 (all L) ACT; 1.6 with sigma_ln=2 weighting; response term 2.7 sigma; 50/deg^2 clean spectra 3.0-3.1 (40-300). Report: Sec. 2.5 fixes, Sec. 4.2/4.3/4.4 rewritten, new Sec. 5 (pair-template estimator), App. B (response to review). PLAN.md v2. Codex review 2 (estimator only) launched.
- 2026-09-10: review 2 incorporated (report Sec. 5 rewritten, Sec. 4.2 extended, App. B; Mathematica checks 9-10; PLAN.md Sec. 0 revised; IMPLEMENTATION.md rewritten with Stage A mock validation as gate). Review 3 launched.
- 2026-09-10: review 3 incorporated; IMPLEMENTATION.md v3 (sign, moments, junk band, completeness division, mock correlations, 11 acceptance tests). Stage A implementation handed to codex gpt-5.6.sol.
- 2026-09-10: data downloaded to /data/LyaLenser/raw (DESI DR1 deltas 1028 files verified by size, QSO catalogues + LSS randoms, ACT DR6 release + baseline products, Planck PR4 2018-like maps; ACT 400 sims still downloading). code/data_checks/delta_summary.py: 428,403 forests; coverage 22/deg^2 at z=2.1 -> 4.6 at z=3.0 (slab mean 12.5); median P_N 0.54 Mpc/h, ln-scatter 2.2. Forecast rerun with empirical n(z) and P_N,a ("DR1 empirical" row, report Table 2): matched-template S/N 0.9 (40<=L<=300), 1.1 (all L); response term 1.8 sigma. Report conclusion updated: DR1 gives pipeline + response + limits; detection needs full DESI.
- 2026-09-11: implementation round 2 finished (316 min): 17/17 unit tests; full scale=1 x 20-seed acceptance 11/11 reported PASS (normalisation slopes 1.007/0.990, injection 1.047, shape p=0.78, jackknife ratio 0.87, spectra within 3%). Weak passes: response prediction 459+-410, combined recovery 1.38+-0.79, stochastic normalisation -0.55+-1.69, matched mean field 41+-22; xi smoothing width tuned on the normalisation test (circular). Launched gpt-6-astra review of code + validation (round-3 decision depends on it; per user, round 3 would be astra).
- 2026-09-11: gpt-6-astra code review of round 2 (report/reviews/codex_review_4_code.md): core algebra verified; NOT passed. Blockers: paired-difference acceptance statistics (absolute deprojected recovery -0.6+-3.1), incomplete tapered band basis (unit signal -> 2), tuned baseline/0.903 factor; plus continuum projection, periodic interpolation, RSD aliasing, matched shot noise (x0.012), sub-slab selection, convergence step. Round 3 brief: code/pipeline/ITERATION3.md. Round 3 launched with gpt-6-astra per user decision.
- 2026-09-11: round 3 (gpt-6-astra, 6.8 h) finished: 36 tests; 20/30 frozen gates PASS; Stage B blocked (report/mock_validation.md, GATES.md, report/reviews/impl_round3_summary.md). Reviewer diagnostics on saved products: matched template recovers 0.87 of the same-range slab convergence (13% deficit, origin unknown); validation used shared sightline/template selection (all deprojection failures have the coupling signature); A_true=5,10 slope rule was a mistake (non-linear remapping), A=1 recovery 0.80+-0.46 sparse / 0.88 dense; measured-vs-analytic xi' ratio 1.19. Round-4 brief drafted: code/pipeline/ITERATION4.md (awaiting go-ahead; ~6 h of astra time).
