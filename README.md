# LyaLenser

Weak lensing of the DESI DR1 Lyman-alpha forest, measured with a pair-based quadratic estimator against
deflection templates built from low-redshift galaxy and quasar catalogues (DESI DR1 LSS and BOSS DR12).
The paper lives in `Paper/` (a symlink to the Overleaf repository); this repository holds the code that
produced every number and figure in it.

**Result (fiducial, `results/joint_fit_v4.json`):** lensing amplitude relative to Planck-2018 LCDM
A = 0.61 +- 0.33 from the combination of the forest auto-correlation (0.42 +- 0.43) and the quasar-forest
cross-correlation (0.86 +- 0.47), jackknife errors. The deflection templates reproduce the ACT DR6 and
Planck PR4 convergence maps to 0.99 +- 0.04 and 0.93 +- 0.03 (`results/deflection_cmb_check_v4.json`).

## Layout

| Path | What it is |
|---|---|
| `lyalenser/` | The Python package: estimator, correlation-table fits, templates, data readers. |
| `scripts/` | Command-line drivers, one per pipeline stage, plus the paper/report figure scripts. `run_fiducial_chain.sh` runs the whole measurement. `slurm/` has the NERSC job scripts. |
| `tests/` | Unit and regression tests (`pytest tests`, about 2 minutes). |
| `results/` | The committed JSON products every quoted number comes from (see `results/README.md`). |
| `report/` | The internal technical report (`lowz.tex`), more detailed than the paper. |
| `docs/` | Documentation: [pipeline](docs/pipeline.md), [data](docs/data.md), [computing](docs/computing.md), [history](docs/history.md). |
| `legacy/` | Archived earlier threads: the CMB-lensing cross-correlation theory and the Stage-A mock campaigns. Not maintained; see `legacy/README.md` for what was learned there. |
| `Paper/` | Symlink to the paper's own git repository (Overleaf). |

## Quick start

```bash
# environment (python 3.11; NaMaster comes from conda-forge as `namaster`)
conda create -p ./env python=3.11 namaster -c conda-forge && conda activate ./env
pip install -r requirements.txt            # or: pip install -e .   (optional; the scripts find the package without it)
export LYALENSER_DATA=/data/LyaLenser       # raw data and products live outside the repository, see docs/data.md

pytest tests                                # 2 minutes
python scripts/paper_figures.py --robustness   # any figure script regenerates from results/ and $LYALENSER_DATA
```

The full measurement (`scripts/run_fiducial_chain.sh`) takes about two hours on a 24-thread workstation
with 40 GB of memory; `docs/pipeline.md` describes every stage, its inputs, outputs and cost.

## How the measurement works, in one paragraph

Lensing displaces the sightlines, so the observed correlation of two forest pixels at transverse separation
r_perp is the unlensed one evaluated at r_perp + delta r_perp. For every pixel pair the estimator forms the
product of the two deltas minus the fitted correlation, weighted by the derivative of that correlation
(the response kernel), and projects it on a template of the expected deflection difference. The template is
the gradient of the lensing potential estimated from foreground tracers at 0.1 < z < 1.6, Wiener-filtered
slice by slice and split into five multipole bands, each with a curl partner and a junk component that are
fitted and marginalised. The amplitude of the science bands relative to the template's own prediction is A.
Errors come from a HEALPix jackknife on the pair midpoints; the same estimator applied to quasar-forest pixel
pairs gives an almost independent second measurement. Details: `docs/pipeline.md`, the paper, and `report/lowz.tex`.
