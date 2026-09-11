#!/bin/bash
# HTCondor wrapper: one Stage A mock seed. Args: <seed> <mock-root> [threads]
set -e
SEED=$1; ROOT=$2; THREADS=${3:-8}
source ~/.bashrc; conda activate lyalenser
export NUMBA_NUM_THREADS=$THREADS OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MPLCONFIGDIR=/tmp/mpl-$USER
cd ~/LyaLenser/code/pipeline
python run_mock_validation.py --phase seed --seed $SEED --mock-root $ROOT
