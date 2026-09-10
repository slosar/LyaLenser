# LyaLenser

Detecting weak lensing of the Lyman-alpha forest by cross-correlating DESI DR1 forest data with CMB lensing (ACT / Planck).

## Layout
- `INSTRUCTIONS.md` — the task contract from the user; read it first.
- `PROGRESS.md` — running log of what is done / in flight / pending (including pending codex reviews). Update it at every milestone.
- `MEMORY.md` — repo-level caveats and machine-specific notes needed to move the work elsewhere (distinct from Claude's personal memory).
- `report/` — single LaTeX report `main.tex` (one section per step), bibliography `main.bib` generated with adstex, figures in `report/figures/`.
- `code/` — python scripts that produce numbers and plots for the report. Each script is runnable standalone and writes figures to `report/figures/`.
- `mathematica/` — Wolfram Language scripts used to derive or verify results. Every script ends with an explicit `Print[...]` because the MCP only shows printed output.

## Tooling
- Python: system anaconda (`/home/anze/anaconda3/bin/python3`, 3.11) with numpy, scipy, matplotlib, healpy, astropy. Create a `.venv` with `uv` only if something is missing, and record it here.
- Mathematica: MCP server `mathematica` (`execute_mathematica`, `verify_derivation`); `wolframscript` also exists at `/usr/bin/wolframscript`.
- Codex: MCP server `codex`, model gpt-6-astra high, for adversarial review at crucial points. If it is down, log the checkpoint under "Pending codex reviews" in PROGRESS.md and continue.
- LaTeX: `latexmk -pdf main.tex` in `report/`. References via the `adstex-references` skill.
- Git: local repo plus private GitHub remote `slosar/LyaLenser`. Commit at milestones.

## Conventions
- Flat sky, Fourier convention: f(x) = ∫ d^2l/(2π)^2 f(l) e^{i l·x}; <f(l) g*(l')> = (2π)^2 δ(l−l') C_l.
- Convergence κ, lensing potential φ with κ = −∇²φ/2 (so κ(L) = L²φ(L)/2), deflection α = ∇φ.
- Forest field δ_F with anisotropic power spectrum P_F(k_∥, k_⊥) in the flat-sky, distant-observer limit; comoving transverse distance χ at the effective forest redshift.
- Notation and derivations in the report are the source of truth; code and Mathematica scripts reference the equation labels they implement.
