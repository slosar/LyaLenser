#!/bin/bash
# HTCondor wrapper: one Stage A mock seed. Args: <seed> <mock-root> [threads]
set -e
SEED=$1; ROOT=$2; THREADS=${3:-8}
REPO=${LYALENSER_REPO:-/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser}
PYTHON=${LYALENSER_PYTHON:-/gpfs/mnt/gpfs02/astro/workarea/anze/envs/lyalenser/bin/python}
export NUMBA_NUM_THREADS=$THREADS OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MPLCONFIGDIR=/tmp/mpl-$USER
cd $REPO/code/pipeline
$PYTHON run_mock_validation.py --phase seed --seed $SEED --mock-root $ROOT
