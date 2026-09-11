# Running LyaLenser on BNL RACF (HTCondor)

## 1. Move the code and data
```bash
# on RACF
git clone git@github.com:slosar/LyaLenser.git ~/LyaLenser          # code (do not copy the Dropbox folder)
export LYALENSER_DATA=/path/to/scratch/LyaLenser                   # any writable location; set it in ~/.bashrc
mkdir -p $LYALENSER_DATA
# on the workstation: data (raw ~80 GB incl. 60 GB of ACT sims, mocks ~10 GB)
rsync -avP --exclude 'mocks/val_scale0.5' /data/LyaLenser/raw  racf:$LYALENSER_DATA/
rsync -avP /data/LyaLenser/mocks racf:$LYALENSER_DATA/             # optional: existing mock products
```
Everything in the code reads `LYALENSER_DATA` (see `code/paths.py`; default `/data/LyaLenser`). If the ACT
simulations are not needed at first, skip `raw/act/baseline/simulations` (60 GB).

## 2. Python environment
```bash
conda create -n lyalenser python=3.11 -y && conda activate lyalenser
pip install -r ~/LyaLenser/requirements.txt          # numpy scipy numba healpy astropy camb h5py matplotlib pytest fitsio
cd ~/LyaLenser/code/pipeline && python -m pytest -q tests   # 36 tests, ~1 min
```
Mathematica and Codex are not needed on RACF (derivations and reviews stay on the workstation).

## 3. Job structure
The Stage A validation is embarrassingly parallel over mock seeds: one seed = generate mock (peak ~28 GB RSS at
scale 1, ~2 min on 24 threads) + pair catalogue + all fits (~3 min). The iteration-4 runner (`campaign4.py`, reached
through `run_mock_validation.py`) exposes idempotent, provenance-hashed phases `dev-seed`, `freeze`, `seed`
(`--variant sparse|dense`), `control --name X`, `collect`, `rebuild`; each honours `NUMBA_NUM_THREADS`. The exact
ordered job list with resource requests is `CAMPAIGN.md`. Submit every phase with the generic template:
```bash
export LYALENSER_DATA=/gpfs/mnt/gpfs02/astro/workarea/anze/Data/LyaLenser
mkdir -p $LYALENSER_DATA/condor_logs/iteration4
cd condor
condor_submit phase.sub -a 'tag=dev100' -a "args=--phase dev-seed --seed 100 --scale 1 --mock-root $LYALENSER_DATA/mocks/iteration4"
condor_submit phase.sub -a 'tag=dense200' -a 'mem=48 GB' -a 'disk=40 GB' -a "args=--phase seed --seed 200 --variant dense --scale 1 --mock-root $LYALENSER_DATA/mocks/iteration4"
```
(`phase.sub` + `run_phase.sh`; variables `campaign`, `tag`, `cpus`, `mem`, `disk`, `args`; logs land in
`$LYALENSER_DATA/condor_logs/<campaign>/<tag>.{out,err,log}` — they must be on gpfs, a `/tmp` path puts the job
on hold.) `seed.sub`/`run_seed.sh` are the older one-job-per-line variant for sparse seeds only.
Request 32 GB memory and 8 CPUs per sparse/dev job; dense noiseless variants (100 sightlines/deg^2) need ~1.5x
the pair-catalogue memory, request 48 GB for those and for the controls. Jobs typically idle ~30 min before the
pool serves group_astro.

## 4. Stage B (real DR1 data) on RACF
The DR1 pair pass over 428k forests is a single-node job (~1e11 pixel pairs, ~10 min at 8e8 pairs/s) with a
~3 GB catalogue; templates and amplitude fits are minutes. Injections and per-sub-slab runs parallelise over condor
the same way (one job per injection realisation / sub-slab).
