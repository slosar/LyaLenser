#!/bin/bash
#SBATCH -A m4895
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 02:00:00
#SBATCH -J lya-refit
#SBATCH -o /global/cfs/cdirs/m4895/users/anze/LyaLenser_data/slurm_logs/stageb/refit_%x_%j.out
# Iteration 15: refit one robustness variant with the derivative-map templates from its saved catalogues.
# Variables: TAG (existing product suffix, e.g. rp20), NEW (new suffix, e.g. rp20_v5), LOWZ (template dir with
# dkappa maps), BANDS, EXTRA (the variant's driver options, e.g. "--rperp-max 20"), STEPS ("cross auto" default),
# AUTO_REF (auto run whose sightlines the cross driver reads; default dr1_lowz_v7d).
export LYALENSER_DATA=/global/cfs/cdirs/m4895/users/anze/LyaLenser_data
export LYALENSER_REPO=/global/cfs/cdirs/m4895/users/anze/LyaLenser_iter11
export NUMBA_NUM_THREADS=64 OMP_NUM_THREADS=8 NUMBA_CACHE_DIR=$SCRATCH/numba_cache_$NEW
PY=/global/cfs/cdirs/m4895/users/anze/envs/lyalenser/bin/python
cd $LYALENSER_REPO/scripts
STEPS=${STEPS:-"cross auto"}; AUTO_REF=${AUTO_REF:-dr1_lowz_v7d}
if [[ " $STEPS " == *" cross "* ]]; then
echo "== cross $TAG -> $NEW $(date)"
$PY run_dr1_qso_lowz.py --out $LYALENSER_DATA/stageb/dr1_qso_$NEW --refit-from $LYALENSER_DATA/stageb/dr1_qso_$TAG --auto-run $LYALENSER_DATA/stageb/$AUTO_REF --lowz $LYALENSER_DATA/$LOWZ --basis $LYALENSER_DATA/stageb/basis_cross_dr1_z196.h5 \
    --randoms 0 --sub-slabs --spline-fixed-base --bands $BANDS --nside-alpha 2048 $EXTRA
fi
if [[ " $STEPS " == *" auto "* ]]; then
echo "== auto $TAG -> $NEW $(date)"
$PY run_dr1_lowz.py --out $LYALENSER_DATA/stageb/dr1_lowz_$NEW --refit-from $LYALENSER_DATA/stageb/dr1_lowz_$TAG --lowz $LYALENSER_DATA/$LOWZ --basis $LYALENSER_DATA/stageb/basis_dr1_ab_z196.h5 \
    --regions lya lyb --xi-correction spline --randoms 0 --nside-jk 8 --zmin 1.96 --zmax 3.0 --zeff 2.3476 --bands $BANDS --nside-alpha 2048 $EXTRA
fi
echo "== done $NEW $(date)"
