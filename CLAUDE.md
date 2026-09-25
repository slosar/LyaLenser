# LyaLenser

Weak lensing of the DESI DR1 Lyman-alpha forest with deflection templates from low-redshift tracers.
Read `README.md` first; `docs/pipeline.md` explains every stage, `docs/data.md` the data tree,
`docs/computing.md` the environments, `docs/history.md` why the analysis looks the way it does.

## Layout
- `lyalenser/` package, `scripts/` drivers and figure scripts, `tests/` (`pytest tests`, 2 min), `results/` committed JSON
  products, `report/` the internal report (`lowz.tex`), `docs/`, `legacy/` archived CMB-lensing and mock threads.
- `Paper/` is a symlink to the paper's own git repository (Overleaf, branch `main`, REVTeX, `pdflatex` + `bibtex`).

## Rules
- Python: `/data/LyaLenser/envs/lyalenser/bin/python` (conda, NaMaster); data under `LYALENSER_DATA` (default `/data/LyaLenser`), never in the repository.
- Never quote a number in the paper or the report that is not in a committed file under `results/`; every figure environment carries a `%% To reproduce:` comment with the command.
- Bibliography with `adstex` (arXiv-ID cite keys), never hand-written BibTeX. Report: `cd report && latexmk -pdf lowz.tex`.
- Commit at milestones; the remote is the private `slosar/LyaLenser`.
- Conventions: flat sky, kappa = -nabla^2 phi / 2, alpha = grad phi; distances in Mpc/h; forest correlation from the fitted table (`lyalenser.tables`), never from the mocks.
