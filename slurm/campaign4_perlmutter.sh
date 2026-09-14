#!/bin/bash
# Stage A campaign on Perlmutter (iteration 5 default): basis -> dev seeds -> freeze -> sparse + dense seeds ->
# controls -> collect, as Slurm jobs with dependencies. Jobs are batched to finish inside 2 h (the `preempt` QOS
# requires -t >= 2 h and preempts only after 2 h of runtime); every phase is idempotent and provenance-checked, so
# the driver can simply be rerun after any failure (phases already complete are no-ops; their jobs are skipped).
# Environment: NERSC_ACCOUNT, LYALENSER_DATA, LYALENSER_REPO, LYALENSER_PYTHON; optional NERSC_QOS (preempt),
# CAMPAIGN_NAME (iteration5; mock root and log directory), CAMPAIGN_SCALE (1), SPARSE_SEEDS / DENSE_SEEDS
# (space-separated; empty string skips), SPARSE_PER_JOB (12), DENSE_PER_JOB (2).
# Measured at scale 1 (iteration 4): sparse seed 68 min at 8 threads, 26.9 GB; dev seed 62-114 min; dense 84 min at 64 threads.
set -euo pipefail
: "${NERSC_ACCOUNT:?export NERSC_ACCOUNT=<allocation>}" "${LYALENSER_DATA:?}" "${LYALENSER_REPO:?}" "${LYALENSER_PYTHON:?}"
QOS=${NERSC_QOS:-preempt}
NAME=${CAMPAIGN_NAME:-iteration5}
SCALE=${CAMPAIGN_SCALE:-1}
SPARSE=${SPARSE_SEEDS-$(seq 0 59)}
DENSE=${DENSE_SEEDS-$(seq 200 209)}
SPJ=${SPARSE_PER_JOB:-12}; DPJ=${DENSE_PER_JOB:-2}
ROOT=$LYALENSER_DATA/mocks/$NAME
LOGS=$LYALENSER_DATA/slurm_logs/$NAME; LISTS=$LOGS/lists; mkdir -p "$LISTS"
REPORT=$LYALENSER_REPO/report; [ "$NAME" = iteration5 ] || REPORT=$REPORT/$NAME
cd "$(dirname "$0")"
common="--scale $SCALE --mock-root $ROOT"
export CAMPAIGN_LOGS=$LOGS
sub() { sbatch --parsable -A "$NERSC_ACCOUNT" -q "$QOS" "$@"; }
done_phase() { [ -f "$ROOT/$1/complete.json" ]; }
dep() { [ -n "$1" ] && echo "--dependency=afterok:$1" || true; }

# Front: basis, development seeds, freeze (skipped when complete).
front=""
if ! done_phase basis; then
  echo "basis|--phase basis $common" > "$LISTS/basis.list"
  front=$(sub -t 02:00:00 -J ly5-basis phase_list.sbatch "$LISTS/basis.list" 1 32); echo "basis: job $front"
fi
devlist=""; for s in 100 101 102 103 104; do done_phase dev/$s || devlist="$devlist$s "; done
if [ -n "$devlist" ]; then
  : > "$LISTS/dev.list"; for s in $devlist; do echo "dev$s|--phase dev-seed --seed $s $common" >> "$LISTS/dev.list"; done
  front=$(sub -t 02:00:00 -J ly5-dev $(dep "$front") phase_list.sbatch "$LISTS/dev.list" $(wc -w <<< "$devlist") 25); echo "dev seeds ($devlist): job $front"
fi
if ! done_phase freeze; then
  echo "freeze|--phase freeze $common" > "$LISTS/freeze.list"
  front=$(sub -t 02:00:00 -J ly5-freeze $(dep "$front") phase_list.sbatch "$LISTS/freeze.list" 1 32); echo "freeze: job $front"
fi
deps=""

# Sparse seeds: SPJ concurrent per node, 128/SPJ threads each (~1 h per batch).
i=0; batch=(); sparse_jobs=""
flush_sparse() {
  [ ${#batch[@]} -gt 0 ] || return 0
  f=$LISTS/sparse_$i.list; printf '%s\n' "${batch[@]}" > "$f"
  j=$(sub -t 02:00:00 -J ly5-sparse$i $(dep "$front") phase_list.sbatch "$f" ${#batch[@]} $((128/SPJ)))
  echo "sparse batch $i (${#batch[@]} seeds, $((128/SPJ)) threads each): job $j"; deps="$deps:$j"; sparse_jobs="$sparse_jobs $j"
  i=$((i+1)); batch=()
}
for s in $SPARSE; do done_phase sparse/$s && continue; batch+=("sparse$s|--phase seed --seed $s --variant sparse $common"); [ ${#batch[@]} -ge "$SPJ" ] && flush_sparse; done
flush_sparse

# Dense seeds: DPJ concurrent per node, 128/DPJ threads each.
i=0; batch=()
flush_dense() {
  [ ${#batch[@]} -gt 0 ] || return 0
  f=$LISTS/dense_$i.list; printf '%s\n' "${batch[@]}" > "$f"
  j=$(sub -t 02:00:00 -J ly5-dense$i $(dep "$front") phase_list.sbatch "$f" ${#batch[@]} $((128/DPJ)))
  echo "dense batch $i (${#batch[@]} seeds, $((128/DPJ)) threads each): job $j"; deps="$deps:$j"
  i=$((i+1)); batch=()
}
for s in $DENSE; do done_phase dense/$s && continue; batch+=("dense$s|--phase seed --seed $s --variant dense $common"); [ ${#batch[@]} -ge "$DPJ" ] && flush_dense; done
flush_dense

# Controls need sparse seed 0: after the batch that contains it (or after the front if it is already complete).
: > "$LISTS/controls.list"; for n in numerical injection flags random benchmark; do echo "ctl_$n|--phase control --name $n $common" >> "$LISTS/controls.list"; done
cdep=$front
j0=$(grep -l '^sparse0|' "$LISTS"/sparse_*.list 2>/dev/null | head -1 || true)
if [ -n "$j0" ]; then idx=${j0##*sparse_}; idx=${idx%.list}; cdep=$(echo $sparse_jobs | awk -v n=$((idx+1)) '{print $n}'); fi
C=$(sub -t 02:00:00 -J ly5-controls $(dep "$cdep") phase_list.sbatch "$LISTS/controls.list" 5 24)
echo "controls: job $C"; deps="$deps:$C"

echo "collect|--phase collect $common --output $REPORT" > "$LISTS/collect.list"
D=$(sub -t 02:00:00 -J ly5-collect --dependency=afterok${deps} phase_list.sbatch "$LISTS/collect.list" 1 32)
echo "collect: job $D (after${deps})"
echo "QOS $QOS; logs: $LOGS/<tag>.out ; squeue -u $USER"
