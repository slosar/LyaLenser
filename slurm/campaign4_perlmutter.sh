#!/bin/bash
# Iteration-4 scale-1 campaign on Perlmutter: same phases and order as condor/campaign4.sh, as four Slurm jobs with
# dependencies. Development seeds and the freeze must already be in $LYALENSER_DATA/mocks/iteration4 (see README).
# Environment: NERSC_ACCOUNT, LYALENSER_DATA, LYALENSER_REPO, LYALENSER_PYTHON; optional SPARSE_SEEDS / DENSE_SEEDS
# (space-separated lists; empty string skips that job), CAMPAIGN_SCALE (1).
set -euo pipefail
: "${NERSC_ACCOUNT:?export NERSC_ACCOUNT=<allocation>}" "${LYALENSER_DATA:?}" "${LYALENSER_REPO:?}" "${LYALENSER_PYTHON:?}"
SCALE=${CAMPAIGN_SCALE:-1}
SPARSE=${SPARSE_SEEDS-$(seq 0 59)}
DENSE=${DENSE_SEEDS-$(seq 200 209)}
ROOT=$LYALENSER_DATA/mocks/iteration4
LISTS=$LYALENSER_DATA/slurm_logs/iteration4/lists; mkdir -p "$LISTS"
cd "$(dirname "$0")"
common="--scale $SCALE --mock-root $ROOT"

: > "$LISTS/sparse.list"; for s in $SPARSE; do echo "sparse$s|--phase seed --seed $s --variant sparse $common" >> "$LISTS/sparse.list"; done
: > "$LISTS/dense.list";  for s in $DENSE;  do echo "dense$s|--phase seed --seed $s --variant dense $common" >> "$LISTS/dense.list"; done
: > "$LISTS/controls.list"; for n in numerical injection flags random benchmark; do echo "ctl_$n|--phase control --name $n $common" >> "$LISTS/controls.list"; done
echo "collect|--phase collect $common --output $LYALENSER_REPO/report" > "$LISTS/collect.list"

deps=""
if [ -s "$LISTS/sparse.list" ]; then
  A=$(sbatch --parsable -A "$NERSC_ACCOUNT" -t 08:00:00 -J ly4-sparse phase_list.sbatch "$LISTS/sparse.list" 14 9)
  echo "sparse seeds: job $A ($(wc -l < "$LISTS/sparse.list") phases, 14 x 9 threads)"; deps="$deps:$A"
fi
if [ -s "$LISTS/dense.list" ]; then
  B=$(sbatch --parsable -A "$NERSC_ACCOUNT" -t 10:00:00 -J ly4-dense phase_list.sbatch "$LISTS/dense.list" 10 12)
  echo "dense seeds: job $B ($(wc -l < "$LISTS/dense.list") phases, 10 x 12 threads)"; deps="$deps:$B"
fi
# Controls read sparse seed 0; they follow the sparse job (or run immediately if sparse 0 is already complete).
cdep=""; [ -n "${A:-}" ] && cdep="--dependency=afterok:$A"
C=$(sbatch --parsable -A "$NERSC_ACCOUNT" -t 03:00:00 -J ly4-controls $cdep phase_list.sbatch "$LISTS/controls.list" 5 24)
echo "controls: job $C"; deps="$deps:$C"
D=$(sbatch --parsable -A "$NERSC_ACCOUNT" -t 01:00:00 -J ly4-collect --dependency=afterok${deps} phase_list.sbatch "$LISTS/collect.list" 1 32)
echo "collect: job $D (after${deps})"
echo "logs: $LYALENSER_DATA/slurm_logs/iteration4/<tag>.out ; squeue -u $USER"
