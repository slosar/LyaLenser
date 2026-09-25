#!/bin/bash
#SBATCH -A m4895
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH -t 02:00:00
#SBATCH -J lya-audit
#SBATCH -o /global/cfs/cdirs/m4895/users/anze/LyaLenser_data/audit_20260925/%x_%A_%a.out
# Array 0-4: 20 independent auto nulls each; 5-9: the corresponding cross nulls.
# Task 10: auto expectation, 11: lensable-only expectation, 12: scale density,
# 13: cross injection, 14: auto real-data injection (split for short queues).
set -euo pipefail
export LYALENSER_DATA=/global/cfs/cdirs/m4895/users/anze/LyaLenser_data
export NUMBA_NUM_THREADS=64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export NUMBA_CACHE_DIR=${SCRATCH}/lyalenser_audit_cache_${SLURM_ARRAY_TASK_ID}
export MPLCONFIGDIR=${SCRATCH}/lyalenser_audit_mpl_${SLURM_ARRAY_TASK_ID}
REPO=/global/cfs/cdirs/m4895/users/anze/LyaLenser_audit_20260925
PY=/global/cfs/cdirs/m4895/users/anze/envs/lyalenser/bin/python
cd "$REPO"
task=${SLURM_ARRAY_TASK_ID}
if (( task < 10 )); then
    stat=auto
    if (( task >= 5 )); then stat=cross; fi
    first=$(( (task % 5) * 20 ))
    "$PY" scripts/random_template_null.py --statistic "$stat" --start "$first" --stop "$((first+20))" \
        --out "$LYALENSER_DATA/audit_20260925/null_${stat}_${first}.json"
elif (( task == 10 )); then
    "$PY" scripts/injection_dr1.py --auto "$LYALENSER_DATA/stageb/dr1_lowz_v7d" \
        --cross "$LYALENSER_DATA/stageb/dr1_qso_v1d" --lowz "$LYALENSER_DATA/lowz_v4" \
        --bands 40 200 400 600 800 1000 --statistics auto --modes expectation --tag v4_corrected_auto_expectation
elif (( task == 11 )); then
    "$PY" scripts/injection_dr1.py --auto "$LYALENSER_DATA/stageb/dr1_lowz_v7d" \
        --lowz "$LYALENSER_DATA/lowz_v4" --bands 40 200 400 600 800 1000 \
        --tag v4_corrected_nosw --remove-same-wavelength --modes expectation --amplitudes -0.25 0.25
elif (( task == 12 )); then
    "$PY" scripts/scale_sensitivity.py --bands 40 200 400 600 800 1000 --tag v4_corrected
elif (( task == 13 )); then
    "$PY" scripts/injection_dr1.py --auto "$LYALENSER_DATA/stageb/dr1_lowz_v7d" \
        --cross "$LYALENSER_DATA/stageb/dr1_qso_v1d" --lowz "$LYALENSER_DATA/lowz_v4" \
        --bands 40 200 400 600 800 1000 --statistics cross --tag v4_corrected_cross
elif (( task == 14 )); then
    "$PY" scripts/injection_dr1.py --auto "$LYALENSER_DATA/stageb/dr1_lowz_v7d" \
        --lowz "$LYALENSER_DATA/lowz_v4" --bands 40 200 400 600 800 1000 \
        --statistics auto --modes real --tag v4_corrected_auto_real
fi
