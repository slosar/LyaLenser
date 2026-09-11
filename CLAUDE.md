# LyaLenser

Detecting weak lensing of the Lyman-alpha forest by cross-correlating DESI DR1 forest data with CMB lensing (ACT / Planck).

## Layout
- `HANDOFF.md` — start here when resuming (state, next steps, how to run on RACF).
- `INSTRUCTIONS.md` — the task contract from the user.
- `PLAN.md` — the proposed measurement pipeline for DESI DR1 x ACT/Planck and the resource assessment (awaiting codex review).
- `PROGRESS.md` — running log of what is done / in flight / pending (including pending codex reviews). Update it at every milestone.
- `MEMORY.md` — repo-level caveats and machine-specific notes needed to move the work elsewhere (distinct from Claude's personal memory).
- `report/` — single LaTeX report `main.tex` (one section per step), bibliography `main.bib` generated with adstex, figures in `report/figures/`.
- `code/` — python scripts that produce numbers and plots for the report. Each script is runnable standalone and writes figures to `report/figures/`.
- `mathematica/` — Wolfram Language scripts used to derive or verify results. Every script ends with an explicit `Print[...]` because the MCP only shows printed output.

## Tooling
- Python: 3.11 with the packages in `requirements.txt` (workstation: system anaconda at `/home/anze/anaconda3/bin/python3`). Data location via `LYALENSER_DATA` (default `/data/LyaLenser`), see `code/paths.py`.
- Mathematica: MCP server `mathematica` (`execute_mathematica`, `verify_derivation`); `wolframscript` also exists at `/usr/bin/wolframscript`.
- Codex: Claude Code plugin (`codex-companion.mjs task ...`), gpt-6-astra high for adversarial reviews and (by user decision) for implementation rounds; job logs under `~/.claude/plugins/data/codex-openai-codex/state/`. Reviews are saved verbatim in `report/reviews/`. Not available on RACF.
- LaTeX: `latexmk -pdf main.tex` in `report/`. References via the `adstex-references` skill.
- Git: local repo plus private GitHub remote `slosar/LyaLenser`. Commit at milestones.

## Conventions
- Flat sky, Fourier convention: f(x) = ∫ d^2l/(2π)^2 f(l) e^{i l·x}; <f(l) g*(l')> = (2π)^2 δ(l−l') C_l.
- Convergence κ, lensing potential φ with κ = −∇²φ/2 (so κ(L) = L²φ(L)/2), deflection α = ∇φ.
- Forest field δ_F with anisotropic power spectrum P_F(k_∥, k_⊥) in the flat-sky, distant-observer limit; comoving transverse distance χ at the effective forest redshift.
- Notation and derivations in the report are the source of truth; code and Mathematica scripts reference the equation labels they implement.
