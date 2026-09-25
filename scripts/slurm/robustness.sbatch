#!/bin/bash
#SBATCH -A m4895
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH -J lya-robust
#SBATCH -o /global/cfs/cdirs/m4895/users/anze/LyaLenser_data/slurm_logs/stageb/robust_%x_%j.out
# Robustness rows of the paper (iteration 14): one variant per job, cross then auto. Variables:
#   TAG (product suffix), LOWZ (template dir), EXTRA (extra driver options, e.g. "--rperp-max 20"), BANDS,
#   STEPS ("cross auto" by default; set STEPS=cross or STEPS=auto to run one half).
export LYALENSER_DATA=/global/cfs/cdirs/m4895/users/anze/LyaLenser_data
export LYALENSER_REPO=/global/cfs/cdirs/m4895/users/anze/LyaLenser_iter11
export NUMBA_NUM_THREADS=64 OMP_NUM_THREADS=8 NUMBA_CACHE_DIR=$SCRATCH/numba_cache_$TAG
PY=/global/cfs/cdirs/m4895/users/anze/envs/lyalenser/bin/python
cd $LYALENSER_REPO/scripts
STEPS=${STEPS:-"cross auto"}
if [[ " $STEPS " == *" cross "* ]]; then
echo "== cross $TAG $(date)"
$PY run_dr1_qso_lowz.py --out $LYALENSER_DATA/stageb/dr1_qso_$TAG --auto-run $LYALENSER_DATA/stageb/dr1_lowz_v7d --lowz $LYALENSER_DATA/$LOWZ --basis $LYALENSER_DATA/stageb/basis_cross_dr1_z196.h5 \
    --randoms 0 --sub-slabs --spline-fixed-base --bands $BANDS --nside-alpha 2048 $EXTRA
fi
if [[ " $STEPS " == *" auto "* ]]; then
echo "== auto $TAG $(date)"
$PY run_dr1_lowz.py --out $LYALENSER_DATA/stageb/dr1_lowz_$TAG --lowz $LYALENSER_DATA/$LOWZ --basis $LYALENSER_DATA/stageb/basis_dr1_ab_z196.h5 \
    --regions lya lyb --xi-correction spline --randoms 0 --nside-jk 8 --zmin 1.96 --zmax 3.0 --zeff 2.3476 --bands $BANDS --nside-alpha 2048 $EXTRA
fi
echo "== done $TAG $(date)"
