# Progress log

## Status summary
- [x] Clarifying questions answered, INSTRUCTIONS.md updated (2026-09-10)
- [x] Repo scaffolding (CLAUDE.md, PROGRESS.md, MEMORY.md, dirs, git, GitHub)
- [x] Step 1: literature review (report Sec. 1, report/literature_notes.md, main.bib via adstex, 32 refs)
- [x] Step 2: sample-variance-limited estimator formalism (report Sec. 2; Mathematica checks A, B, C pass)
- [x] Step 2 numerics: N_kappa(L), bias-hardened N, S/N (code/make_plots.py -> report/figures, numbers.json)
- [x] Report compiles cleanly (latexmk, 14 pages)
- [ ] Codex adversarial reviews run (see pending list)
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

## Pending codex reviews
1. Derivation of the mode coupling, squeezed limit and amplitude degeneracy (report Sec. 2.2-2.3, mathematica/*.wls).
2. Discrete optimal quadratic estimator and cubic cross-correlation estimator, list of biases (Sec. 2.4-2.5).
3. Continuum-limit N_kappa(L) formula incl. factors of 2 and mode counting (Sec. 2.6, code/recon_noise.py).
4. Bias-hardening formulae and the claim that hardened noise is set by transverse resolution (Sec. 2.7, Sec. 3).
5. Toy numbers: forest power model, CMB white-noise calibration, S/N (code/make_plots.py, Table 1).

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
