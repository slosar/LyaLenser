#!/bin/bash
# Args: <seed> <mock-root> [threads] [sparse|dense] [scale]
set -euo pipefail
SEED=$1; ROOT=$2; THREADS=${3:-8}; VARIANT=${4:-sparse}; SCALE=${5:-1}
REPO=${LYALENSER_REPO:-/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser}
PYTHON=${LYALENSER_PYTHON:-/gpfs/mnt/gpfs02/astro/workarea/anze/envs/lyalenser/bin/python}
export NUMBA_NUM_THREADS=$THREADS OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MPLCONFIGDIR=/tmp/mpl-$USER
cd "$REPO/code/pipeline"
exec "$PYTHON" run_mock_validation.py --phase seed --seed "$SEED" --variant "$VARIANT" --scale "$SCALE" --mock-root "$ROOT"
