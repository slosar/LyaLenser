# Running the iteration-4 campaign on NERSC Perlmutter (CPU nodes)

The campaign is ~1000 core-hours of numba/CPU work in 27-28 GB jobs (60 sparse seeds of ~9 core-h, 10 dense
seeds of ~40 core-h, controls, collect). One Perlmutter CPU node (128 cores, 512 GB) runs 14 sparse seeds at a
time: sparse ~5 h, dense ~4 h on a second node, controls + collect ~1 h. Products made here merge with the RACF
products because every phase is provenance-hashed (code, GATES.md, config, ACT mask, package versions), so the
two sites can work on the same campaign concurrently. No DR1 data is needed for Stage A.

## 1. Inputs to transfer (~3.5 GB)
On the machine that holds the campaign root (RACF: `LYALENSER_DATA=/gpfs/mnt/gpfs02/astro/workarea/anze/Data/LyaLenser`):
```bash
cd $LYALENSER_DATA
tar czf lyalenser_iter4_inputs.tgz mocks/iteration4/dev mocks/iteration4/freeze \
    raw/act/baseline/mask_act_dr6_lensing_v1_healpix_nside_4096_baseline.fits raw/act/baseline/N_L_kk_act_dr6_lensing_v1_baseline.txt
```
`dev/` (five completed development seeds, 1.8 GB) and `freeze/` are required: every seed phase re-verifies their
hashes and reads `dev/100/mock.h5` and `freeze/raw.h5`. Copy the tarball to `$SCRATCH` on Perlmutter (scp or Globus).

## 2. Set up on Perlmutter
```bash
git clone git@github.com:slosar/LyaLenser.git $SCRATCH/LyaLenser      # commit >= 2aba3dd
module load python
conda create -p $SCRATCH/envs/lyalenser python=3.11 -y && conda activate $SCRATCH/envs/lyalenser
pip install -r $SCRATCH/LyaLenser/slurm/requirements-perlmutter.txt   # exact pins: the versions are part of the provenance
export LYALENSER_DATA=$SCRATCH/LyaLenser_data && mkdir -p $LYALENSER_DATA && tar xzf $SCRATCH/lyalenser_iter4_inputs.tgz -C $LYALENSER_DATA
sha256sum $LYALENSER_DATA/raw/act/baseline/mask_act_dr6_lensing_v1_healpix_nside_4096_baseline.fits
#   must be 2bca9299511222433b5db2f3bae29bbe042ba2e9a23578d603c26ad992f6f18a
cd $SCRATCH/LyaLenser/code/pipeline && NUMBA_NUM_THREADS=8 python -m pytest -q tests        # 51 passed (~5 min)
python run_mock_validation.py --phase freeze --scale 1 --mock-root $LYALENSER_DATA/mocks/iteration4
#   prints "COMPLETE freeze" without computing anything iff the provenance (code + versions + mask) matches RACF;
#   a RuntimeError here means a version or file differs: fix that before submitting anything.
```

## 3. Launch
```bash
export NERSC_ACCOUNT=<your allocation>            # e.g. desi
export LYALENSER_DATA=$SCRATCH/LyaLenser_data LYALENSER_REPO=$SCRATCH/LyaLenser LYALENSER_PYTHON=$SCRATCH/envs/lyalenser/bin/python
cd $SCRATCH/LyaLenser/slurm && ./campaign4_perlmutter.sh
```
This submits four jobs: sparse seeds (1 node, 14 concurrent x 9 threads), dense seeds (1 node, 10 x 12 threads),
controls (after sparse), collect (after all). Per-phase logs: `$LYALENSER_DATA/slurm_logs/iteration4/<tag>.out`;
a phase that already has `complete.json` with matching provenance is a no-op, so the driver can be resubmitted
after any failure. To run only a subset (e.g. dense seeds while RACF does the sparse ones):
`SPARSE_SEEDS="" ./campaign4_perlmutter.sh` or `DENSE_SEEDS="" ...`.

## 4. Merge back and collect
Copy finished products into the RACF root (directories are disjoint per seed, so rsync is safe; a seed finished
on both sites is identical in content up to floating-point noise, keep either):
```bash
rsync -a $LYALENSER_DATA/mocks/iteration4/{sparse,dense,controls} racf:/gpfs/mnt/gpfs02/astro/workarea/anze/Data/LyaLenser/mocks/iteration4/
```
Then on RACF `condor/campaign4.sh` skips whatever is complete and runs the rest, ending with `collect` (or run
`python run_mock_validation.py --phase collect --scale 1 --mock-root $LYALENSER_DATA/mocks/iteration4 --output ../../report`
directly on the login node: collect only reads products). Alternatively run collect on Perlmutter and copy
`report/mock_validation.{md,json}` and `report/figures/mock_iteration4_*.pdf` back.
