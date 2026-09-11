#!/bin/bash
# HTCondor wrapper for any iteration-4 phase. Args are passed verbatim to run_mock_validation.py,
# e.g.: run_phase.sh --phase seed --seed 7 --variant sparse --scale 1 --mock-root $LYALENSER_DATA/mocks/iteration4
set -euo pipefail
REPO=${LYALENSER_REPO:-/gpfs/mnt/gpfs02/astro/workarea/anze/work/LyaLenser}
PYTHON=${LYALENSER_PYTHON:-/gpfs/mnt/gpfs02/astro/workarea/anze/envs/lyalenser/bin/python}
export LYALENSER_DATA=${LYALENSER_DATA:-/gpfs/mnt/gpfs02/astro/workarea/anze/Data/LyaLenser}
export NUMBA_NUM_THREADS=${NUMBA_NUM_THREADS:-8} OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MPLCONFIGDIR=/tmp/lyalenser-mpl-$USER
echo "host=$(hostname) start=$(date -u +%FT%TZ) threads=$NUMBA_NUM_THREADS args=$*"
cd "$REPO/code/pipeline"
exec "$PYTHON" run_mock_validation.py "$@"
