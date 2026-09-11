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
scale 1, ~2 min on 24 threads) + pair catalogue + all fits (~3 min). The round-3 runner exposes this only through
`--phase validation` (serial loop, 24 threads hard-coded); round 4 adds `--phase seed --seed N` (idempotent,
resumable, honours `NUMBA_NUM_THREADS`) and `--phase collect`. The DAG is then:
```
prepare  : python run_mock_validation.py --phase development ; --phase freeze      (once, one node)
seed_N   : python run_mock_validation.py --phase seed --seed N --mock-root $LYALENSER_DATA/mocks/iteration4   (N = 0..59, parallel)
collect  : python run_mock_validation.py --phase collect ; --phase rebuild        (one node)
```
Use the templates in this directory: `seed.sub` (one job per seed via `queue seed from seeds.txt`), `run_seed.sh`.
Request 32 GB memory and 8 CPUs per seed job (numba scales well to 8; 24 is not necessary when 60 jobs run at once).
Dense noiseless variants (100 sightlines/deg^2) need ~1.5x the pair-catalogue memory; request 48 GB for those.

## 4. Stage B (real DR1 data) on RACF
The DR1 pair pass over 428k forests is a single-node job (~1e11 pixel pairs, ~10 min at 8e8 pairs/s) with a
~3 GB catalogue; templates and amplitude fits are minutes. Injections and per-sub-slab runs parallelise over condor
the same way (one job per injection realisation / sub-slab).
