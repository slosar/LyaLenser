#!/bin/bash
# Iteration-4 scale-1 campaign on Perlmutter: the phases and order of condor/campaign4.sh as Slurm jobs with
# dependencies. Development seeds and the freeze must already be in $LYALENSER_DATA/mocks/iteration4 (README).
# Jobs are batched so that each finishes well inside 2 h: on the `preempt` QOS (default; it starts within
# minutes to hours while `regular` waits for days) a job becomes preemptible only after 2 h of runtime, so these
# jobs are never preempted, and every phase is idempotent anyway.
# Environment: NERSC_ACCOUNT, LYALENSER_DATA, LYALENSER_REPO, LYALENSER_PYTHON; optional NERSC_QOS (preempt),
# SPARSE_SEEDS / DENSE_SEEDS (space-separated; empty string skips), SPARSE_PER_JOB (12), DENSE_PER_JOB (2),
# CAMPAIGN_SCALE (1). Measured at scale 1: sparse seed 68 min at 8 threads, 26.9 GB; dev seed 62-114 min.
set -euo pipefail
: "${NERSC_ACCOUNT:?export NERSC_ACCOUNT=<allocation>}" "${LYALENSER_DATA:?}" "${LYALENSER_REPO:?}" "${LYALENSER_PYTHON:?}"
QOS=${NERSC_QOS:-preempt}
SCALE=${CAMPAIGN_SCALE:-1}
SPARSE=${SPARSE_SEEDS-$(seq 0 59)}
DENSE=${DENSE_SEEDS-$(seq 200 209)}
SPJ=${SPARSE_PER_JOB:-12}; DPJ=${DENSE_PER_JOB:-2}
ROOT=$LYALENSER_DATA/mocks/iteration4
LISTS=$LYALENSER_DATA/slurm_logs/iteration4/lists; mkdir -p "$LISTS"
cd "$(dirname "$0")"
common="--scale $SCALE --mock-root $ROOT"
sub() { sbatch --parsable -A "$NERSC_ACCOUNT" -q "$QOS" "$@"; }
deps=""

# Sparse seeds: SPJ concurrent per node, 128/SPJ threads each (~1 h per batch).
i=0; batch=(); sparse_jobs=""
flush_sparse() {
  [ ${#batch[@]} -gt 0 ] || return 0
  f=$LISTS/sparse_$i.list; printf '%s\n' "${batch[@]}" > "$f"
  j=$(sub -t 01:55:00 -J ly4-sparse$i phase_list.sbatch "$f" ${#batch[@]} $((128/SPJ)))
  echo "sparse batch $i (${#batch[@]} seeds, $((128/SPJ)) threads each): job $j"; deps="$deps:$j"; sparse_jobs="$sparse_jobs $j"
  i=$((i+1)); batch=()
}
for s in $SPARSE; do batch+=("sparse$s|--phase seed --seed $s --variant sparse $common"); [ ${#batch[@]} -ge "$SPJ" ] && flush_sparse; done
flush_sparse

# Dense seeds: DPJ concurrent per node, 128/DPJ threads each (pair accumulation dominates and scales with threads).
i=0; batch=()
flush_dense() {
  [ ${#batch[@]} -gt 0 ] || return 0
  f=$LISTS/dense_$i.list; printf '%s\n' "${batch[@]}" > "$f"
  j=$(sub -t 02:30:00 -J ly4-dense$i phase_list.sbatch "$f" ${#batch[@]} $((128/DPJ)))
  echo "dense batch $i (${#batch[@]} seeds, $((128/DPJ)) threads each): job $j"; deps="$deps:$j"
  i=$((i+1)); batch=()
}
for s in $DENSE; do batch+=("dense$s|--phase seed --seed $s --variant dense $common"); [ ${#batch[@]} -ge "$DPJ" ] && flush_dense; done
flush_dense

# Controls need sparse seed 0: run now if it is complete, otherwise after the batch that contains it.
: > "$LISTS/controls.list"; for n in numerical injection flags random benchmark; do echo "ctl_$n|--phase control --name $n $common" >> "$LISTS/controls.list"; done
cdep=""
if [ ! -f "$ROOT/sparse/0/complete.json" ]; then
  j0=$(grep -l '^sparse0|' "$LISTS"/sparse_*.list 2>/dev/null | head -1)
  [ -n "$j0" ] && { idx=${j0##*sparse_}; idx=${idx%.list}; jid=$(echo $sparse_jobs | awk -v n=$((idx+1)) '{print $n}'); cdep="--dependency=afterok:$jid"; }
fi
C=$(sub -t 01:55:00 -J ly4-controls $cdep phase_list.sbatch "$LISTS/controls.list" 5 24)
echo "controls: job $C ${cdep:+($cdep)}"; deps="$deps:$C"

echo "collect|--phase collect $common --output $LYALENSER_REPO/report" > "$LISTS/collect.list"
D=$(sub -t 00:40:00 -J ly4-collect --dependency=afterok${deps} phase_list.sbatch "$LISTS/collect.list" 1 32)
echo "collect: job $D (after${deps})"
echo "QOS $QOS; logs: $LYALENSER_DATA/slurm_logs/iteration4/<tag>.out ; squeue -u $USER"
