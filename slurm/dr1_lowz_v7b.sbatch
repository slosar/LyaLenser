#!/bin/bash
#SBATCH -A m4895
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH -J dr1-lowz-v7b
#SBATCH -o /global/cfs/cdirs/m4895/users/anze/LyaLenser_data/slurm_logs/stageb/dr1_lowz_v7b_%j.out
# Iteration 13 (alternative band set 40-200 ... 800-1000, deflection at nside 2048), otherwise iteration 11: forest 1.96 <= z <= 3.0 (both delta regions), source plane at the weighted mean pixel redshift
# z_eff = 2.3476, tracers to z = 1.6 (templates lowz_v3), basis for the wider slab, evolving correlation table.
export LYALENSER_DATA=/global/cfs/cdirs/m4895/users/anze/LyaLenser_data
export LYALENSER_REPO=/global/cfs/cdirs/m4895/users/anze/LyaLenser_iter11
export NUMBA_NUM_THREADS=64 OMP_NUM_THREADS=8 NUMBA_CACHE_DIR=$SCRATCH/numba_cache
cd $LYALENSER_REPO/code/stageb
/global/cfs/cdirs/m4895/users/anze/envs/lyalenser/bin/python run_dr1_lowz.py \
    --out $LYALENSER_DATA/stageb/dr1_lowz_v7b --lowz $LYALENSER_DATA/lowz_v3 --basis $LYALENSER_DATA/stageb/basis_dr1_ab_z196.h5 \
    --regions lya lyb --xi-correction spline --randoms 40 --nside-jk 8 --zmin 1.96 --zmax 3.0 --zeff 2.3476 --bands 40 200 400 600 800 1000 --nside-alpha 2048 --randoms 20
