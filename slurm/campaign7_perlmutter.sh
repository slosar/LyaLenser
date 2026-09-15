#!/bin/bash
# Stage A iteration 7 (low-redshift tracers) on Perlmutter: basis -> dev seeds -> freeze -> sparse seeds -> controls
# -> collect, as Slurm jobs with dependencies (preempt QOS, every job sized to finish inside 2 h). Every phase is
# idempotent and provenance-checked (campaign7.py), so the driver can be rerun after any failure.
# Environment: NERSC_ACCOUNT, LYALENSER_DATA, LYALENSER_REPO, LYALENSER_PYTHON; optional NERSC_QOS (preempt),
# CAMPAIGN_NAME (iteration7), CAMPAIGN_SCALE (1), SPARSE_SEEDS (default 4000-4399), DEV_SEEDS (5000-5004),
# SPARSE_PER_JOB (12); SKIP_CONTROLS=1 submits no controls and no collect (seeds only, e.g. when the extras seed
# and the controls run on RACF and collect runs later with the merged products).
set -euo pipefail
: "${NERSC_ACCOUNT:?export NERSC_ACCOUNT=<allocation>}" "${LYALENSER_DATA:?}" "${LYALENSER_REPO:?}" "${LYALENSER_PYTHON:?}"
QOS=${NERSC_QOS:-preempt}
NAME=${CAMPAIGN_NAME:-iteration7}
SCALE=${CAMPAIGN_SCALE:-1}
SPARSE=${SPARSE_SEEDS-$(seq 4000 4399)}
DEV=${DEV_SEEDS-$(seq 5000 5004)}
FIRST=${EXTRAS_SEED:-4000}
SPJ=${SPARSE_PER_JOB:-12}
ROOT=$LYALENSER_DATA/mocks/$NAME
LOGS=$LYALENSER_DATA/slurm_logs/$NAME; LISTS=$LOGS/lists; mkdir -p "$LISTS"
REPORT=$LYALENSER_REPO/report; [ "$NAME" = iteration7 ] || REPORT=$REPORT/$NAME
cd "$(dirname "$0")"
common="--scale $SCALE --mock-root $ROOT"
export CAMPAIGN_LOGS=$LOGS PHASE_SCRIPT=campaign7.py
sub() { sbatch --parsable -A "$NERSC_ACCOUNT" -q "$QOS" "$@"; }
done_phase() { [ -f "$ROOT/$1/complete.json" ]; }
dep() { [ -n "$1" ] && echo "--dependency=afterok:$1" || true; }

front=""
if ! done_phase basis; then
  echo "basis|--phase basis $common" > "$LISTS/basis.list"
  front=$(sub -t 02:00:00 -J ly7-basis phase_list.sbatch "$LISTS/basis.list" 1 32); echo "basis: job $front"
fi
devlist=""; for s in $DEV; do done_phase dev/$s || devlist="$devlist$s "; done
if [ -n "$devlist" ]; then
  : > "$LISTS/dev.list"; for s in $devlist; do echo "dev$s|--phase dev-seed --seed $s $common" >> "$LISTS/dev.list"; done
  front=$(sub -t 02:00:00 -J ly7-dev $(dep "$front") phase_list.sbatch "$LISTS/dev.list" $(wc -w <<< "$devlist") 25); echo "dev seeds ($devlist): job $front"
fi
if ! done_phase freeze; then
  echo "freeze|--phase freeze $common" > "$LISTS/freeze.list"
  front=$(sub -t 02:00:00 -J ly7-freeze $(dep "$front") phase_list.sbatch "$LISTS/freeze.list" 1 32); echo "freeze: job $front"
fi
deps=""
i=0; batch=(); sparse_jobs=""
flush_sparse() {
  [ ${#batch[@]} -gt 0 ] || return 0
  f=$LISTS/sparse_$i.list; printf '%s\n' "${batch[@]}" > "$f"
  j=$(sub -t 02:00:00 -J ly7-sparse$i $(dep "$front") phase_list.sbatch "$f" ${#batch[@]} $((128/SPJ)))
  echo "sparse batch $i (${#batch[@]} seeds, $((128/SPJ)) threads each): job $j"; deps="$deps:$j"; sparse_jobs="$sparse_jobs $j"
  i=$((i+1)); batch=()
}
for s in $SPARSE; do done_phase sparse/$s && continue; batch+=("sparse$s|--phase seed --seed $s $common"); [ ${#batch[@]} -ge "$SPJ" ] && flush_sparse; done
flush_sparse
if [ -n "${SKIP_CONTROLS:-}" ]; then echo "seeds only (SKIP_CONTROLS); QOS $QOS; logs: $LOGS/<tag>.out"; exit 0; fi
: > "$LISTS/controls.list"; for n in numerical injection benchmark; do echo "ctl_$n|--phase control --name $n $common" >> "$LISTS/controls.list"; done
cdep=$front
j0=$(grep -l "^sparse$FIRST|" "$LISTS"/sparse_*.list 2>/dev/null | head -1 || true)
if [ -n "$j0" ]; then idx=${j0##*sparse_}; idx=${idx%.list}; cdep=$(echo $sparse_jobs | awk -v n=$((idx+1)) '{print $n}'); fi
C=$(sub -t 02:00:00 -J ly7-controls $(dep "$cdep") phase_list.sbatch "$LISTS/controls.list" 3 32)
echo "controls: job $C"; deps="$deps:$C"
echo "collect|--phase collect $common --output $REPORT" > "$LISTS/collect.list"
D=$(sub -t 02:00:00 -J ly7-collect --dependency=afterok${deps} phase_list.sbatch "$LISTS/collect.list" 1 32)
echo "collect: job $D (after${deps})"
echo "QOS $QOS; logs: $LOGS/<tag>.out ; squeue -u $USER"
