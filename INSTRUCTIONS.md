# Context

You are in a fresh folder for a new repo. This work will attempt to detect lensing of the Lyman-alpha forest data by cross-correlation of DESI DR1 and ACT or Planck data.

Note: `/home/anze/work` is a symlink to `/home/anze/Dropbox/work`, so `~/work/LyaLenser` and `~/Dropbox/work/LyaLenser` are the same directory.

# Contract

Use CLAUDE.md with the usual repo context. Use PROGRESS.md to track progress in case work gets interrupted (for example, you run out of tokens or machine reboots). Use MEMORY.md (repo level file, separate to your own memory) to store caveats, etc in case we need to move this work to a different machines.

You have access to Mathematica via MCP. Use mathematica to do or confirm any mathematical derivations. Note: The MCP only shows printed output, so every script must end in an explicit print.

You have access to codex via MCP. Use codex on gpt-6-astra high to adversarially review your results at crucial points.
If the codex MCP server is unavailable, proceed anyway and record each pending review checkpoint in PROGRESS.md (section "Pending codex reviews") so the reviews can be run once it is back. Do not stop work waiting for it.

Mathematical write-up should be in latex report in the report/ directory containing derivations and plots.
Conventions: a single evolving document `report/main.tex`, one section per step, literature review as the first section. Bibliography built with adstex from ADS (`report/main.bib`), never hand-written bibtex.

Any code used to make plots and show results should go into code/ directory. Use the system anaconda python (numpy, scipy, matplotlib, healpy, astropy are available). Only if a package is missing, create a `.venv` (with uv) and record it in CLAUDE.md.

Mathematica scripts should go into mathematica/ directory.

Version control: local git plus a private GitHub repository under the `slosar` account (created with `gh`). Commit at milestones.

Data: none is needed for the initial steps. No DESI DR1 deltas or CMB lensing maps are on this machine yet; data access (likely NERSC) will be specified when the subsequent steps are defined.

## Initial Steps

### Review of the existing literature.

Please search and review the existing literature on weak lensing of the Lyman-alpha forest. I don't think there should be much, but you should check.

### Write the lensing cross-correlation estimator in the limit of zero noise (sample variance limit -- totally not applicable but necessary to build intuition).

The effect of weak lensing is to make a pair of LyA forests appear at a different separation than they really are. So a pair of deltaF pixels will, on average, give a different correlation than if no lensing was present.
So an estimator for the kappa_lya map can be generated that is quadratic in delta_Fs.

So an estimator for C_(kappa_lya x kappa_cmb)(ell) is a cubic function W_ijk delta_F_i delta_F_j kappa_CMB_k.

Work out the maths and write down the formalism for the C_ell cross-correlation estimate.

Agreed scope for this step:
 - Flat sky. delta_F is a 3D Gaussian field with a redshift-space (anisotropic) power spectrum P_F(k_parallel, k_perp).
 - Single source plane: kappa is evaluated at one effective source redshift for the whole forest (no variation of the lensing kernel along the forest).
 - Derive the estimator for arbitrary discrete sightlines and pixels (general quadratic weights, hence general W_ijk), then take the dense-sightline continuum limit to obtain closed-form Hu-Okamoto-style expressions for the estimator, its normalisation and its reconstruction noise N_kappa(L).
 - Numerics: evaluate N_kappa(L) and a rough signal-to-noise for C_ell^{kappa_lya x kappa_cmb} with a fiducial forest power spectrum and a DESI-like footprint. No pixel noise or finite sightline density in the noise budget at this stage (sample-variance limit), but note where they would enter.

### Subsequent steps:
 - will write based on results from first two steps.
