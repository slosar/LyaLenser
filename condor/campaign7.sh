#!/bin/bash
# Stage A iteration 7 (low-redshift tracers) on HTCondor: basis -> dev seeds -> freeze -> sparse seeds -> controls
# -> collect (campaign7.py phases; idempotent, provenance-checked; rerun after any failure). Run from condor/ on
# the submit node: `nohup ./campaign7.sh > $LYALENSER_DATA/condor_logs/iteration7/campaign.out 2>&1 &`.
# Environment: LYALENSER_DATA (required); CAMPAIGN_NAME (iteration7), CAMPAIGN_SCALE (1), SPARSE_SEEDS
# ("4000 .. 4399"), DEV_SEEDS ("5000 .. 5004"), EXTRAS_SEED (4000; the controls wait for it), STOP_AFTER_DEV.
# Smoke: CAMPAIGN_NAME=iteration7_smoke CAMPAIGN_SCALE=0.25 SPARSE_SEEDS="4000 4001" ./campaign7.sh
set -euo pipefail
: "${LYALENSER_DATA:?export LYALENSER_DATA first}"
NAME=${CAMPAIGN_NAME:-iteration7}
SCALE=${CAMPAIGN_SCALE:-1}
SPARSE=${SPARSE_SEEDS:-$(seq 4000 4399)}
DEV=${DEV_SEEDS:-$(seq 5000 5004)}
FIRST=${EXTRAS_SEED:-4000}
ROOT=$LYALENSER_DATA/mocks/$NAME
LOGS=$LYALENSER_DATA/condor_logs/$NAME
REPORT=$(cd "$(dirname "$0")/.." && pwd)/report
[ "$NAME" = iteration7 ] || REPORT=$REPORT/$NAME
if [ "$SCALE" = 1 ]; then MEM="32 GB"; else MEM="16 GB"; fi
mkdir -p "$LOGS" "$ROOT"
cd "$(dirname "$0")"

submit() {  # submit <tag> <cpus> <mem> <disk> <args...>
  local tag=$1 cpus=$2 mem=$3 disk=$4; shift 4
  if [ -f "$LOGS/$tag.log" ] && grep -q 'return value 0' "$LOGS/$tag.log"; then echo "skip $tag (done)"; return; fi
  rm -f "$LOGS/$tag.log"
  condor_submit phase.sub -a "campaign=$NAME" -a "script=campaign7.py" -a "tag=$tag" -a "cpus=$cpus" -a "mem=$mem" -a "disk=$disk" \
    -a "args=$* --scale $SCALE --mock-root $ROOT" | tail -1 | sed "s/^/$tag: /"
}
wait_for() {
  local t rc=0
  for t in "$@"; do
    until [ -f "$LOGS/$t.log" ]; do sleep 30; done
    condor_wait -wait 172800 "$LOGS/$t.log" > /dev/null 2>&1 || true
    if grep -q 'return value 0' "$LOGS/$t.log"; then echo "$(date -u +%FT%TZ) $t ok"; else echo "$(date -u +%FT%TZ) $t FAILED"; tail -3 "$LOGS/$t.err"; rc=1; fi
  done
  return $rc
}

echo "$(date -u +%FT%TZ) $NAME scale $SCALE: basis"
submit basis 8 "$MEM" "20 GB" --phase basis
wait_for basis
echo "$(date -u +%FT%TZ) development seeds"
for s in $DEV; do submit dev$s 8 "$MEM" "20 GB" --phase dev-seed --seed $s; done
wait_for $(for s in $DEV; do echo dev$s; done)
[ -z "${STOP_AFTER_DEV:-}" ] || { echo "$(date -u +%FT%TZ) stopping after development seeds"; exit 0; }
echo "$(date -u +%FT%TZ) freeze"
submit freeze 8 "$MEM" "20 GB" --phase freeze
wait_for freeze
echo "$(date -u +%FT%TZ) seeds"
for s in $SPARSE; do submit sparse$s 8 "$MEM" "20 GB" --phase seed --seed $s; done
wait_for sparse$FIRST
for n in numerical injection benchmark; do submit ctl_$n 8 "$MEM" "40 GB" --phase control --name $n; done
wait_for $(for s in $SPARSE; do echo sparse$s; done) ctl_numerical ctl_injection ctl_benchmark
echo "$(date -u +%FT%TZ) collect"
submit collect 8 "$MEM" "40 GB" --phase collect --output "$REPORT"
wait_for collect
echo "$(date -u +%FT%TZ) CAMPAIGN $NAME DONE"
