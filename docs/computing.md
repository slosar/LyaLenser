# Computing

## Environment

Python 3.11 with the packages of `requirements.txt`. NaMaster (`pymaster`) is easiest from conda-forge
(`namaster`); everything else installs with pip. numpy is pinned below 2: NEP 50 promotion changes the float32
arithmetic of the numba pair kernels and `np.trapz` was removed (we use `scipy.integrate.trapezoid`).

- Workstation: `/data/LyaLenser/envs/lyalenser/bin/python` (conda; numpy 1.26.4, numba 0.61.2, scipy 1.16,
  healpy 1.18.1, astropy 6.1, camb 2.0.4, fitsio, h5py, pymaster 2.7). Use this one; the base anaconda has no pymaster.
- NERSC Perlmutter: `/global/cfs/cdirs/m4895/users/anze/envs/lyalenser/bin/python`, repository clone
  `/global/cfs/cdirs/m4895/users/anze/LyaLenser_iter11` (keep it on git HEAD), allocation `m4895`.

The scripts add the repository root to `sys.path`, so no installation is needed; `pip install -e .` works too.

## Threads and memory

The pair accumulation and the correlation cells are numba kernels controlled by `NUMBA_NUM_THREADS`; NaMaster,
healpy and CAMB use OpenMP (`OMP_NUM_THREADS`). Typical settings and costs on the 24-thread, 62 GB workstation:

| Stage | Threads | Wall | Peak memory |
|---|---|---|---|
| `lowz_catalogues.py` | OMP 12 | 10 min | few GB |
| `deflection_cmb_check.py` | OMP 8 | 2 min | few GB |
| `run_dr1_lowz.py` (fiducial) | NUMBA 22, OMP 4 | 75 min | 34 GB |
| `run_dr1_qso_lowz.py` (fiducial) | NUMBA 16, OMP 4 | 16 min | 26 GB |
| `joint_response_fit.py` | NUMBA 12 | 2 min | 10 GB |
| `injection_dr1.py` | NUMBA 16 | 30 min | 30 GB |
| `pytest tests` | NUMBA 8 | 1-2 min | small |

Run the two measurements sequentially: their memory peaks do not fit together. On Perlmutter one CPU node
(`NUMBA_NUM_THREADS=64`) does the auto run in about 25 min; set `NUMBA_CACHE_DIR` under `$SCRATCH`.

## NERSC job scripts

`scripts/slurm/dr1_lowz.sbatch` runs the auto measurement; `scripts/slurm/robustness.sbatch` runs one
robustness variant (`sbatch --export=ALL,TAG=rp20,LOWZ=lowz_v4,BANDS="40 200 400 600 800 1000",EXTRA="--rperp-max 20" robustness.sbatch`;
`STEPS=cross` or `STEPS=auto` for one half). Logs go to `$LYALENSER_DATA/slurm_logs/stageb/`. During a maintenance
reservation the scheduler refuses jobs longer than the drain window (`ReqNodeNotAvail`); shorten `-t`.
The sshproxy certificate expires after 24 h and needs the OTP.

## Gotchas

- `pkill -f <pattern>` from a shell whose own command line contains the pattern kills that shell (exit 144), and
  the rest of a compound command is skipped silently. Use `pgrep -f '[p]attern'` and run `pkill` on its own.
- Session scratch under `/tmp` does not survive a reboot; anything worth keeping goes to `$LYALENSER_DATA`.
- Numba caches (`cache=True` kernels in `pairs`, `xi_model`, `xi_cross`) are invalidated when the files move;
  the first run after a checkout recompiles (a minute).
- `adstex` may write a literal Unicode alpha into a bibliography title; pdflatex then fails. Replace it with
  `{\ensuremath{\alpha}}` in the `.bib`.
- Build the report with `cd report && latexmk -pdf lowz.tex`; the paper with `pdflatex`/`bibtex` in `Paper/`
  (REVTeX; bibliography `references.bib` maintained with `adstex`).
