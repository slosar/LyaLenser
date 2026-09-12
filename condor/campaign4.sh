#!/bin/bash
# Iteration-4 campaign driver (CAMPAIGN.md order) for HTCondor. Run from the repo's condor/ directory on the
# submit node, e.g. `nohup ./campaign4.sh > $LYALENSER_DATA/condor_logs/iteration4/campaign.out 2>&1 &`.
# Every phase is idempotent and provenance-checked, so the script can simply be rerun after a failure.
# Environment: LYALENSER_DATA (required); CAMPAIGN_NAME (default iteration4; mock root and log directory),
# CAMPAIGN_SCALE (1), SPARSE_SEEDS ("0 .. 59"), DENSE_SEEDS ("200 .. 209"), STOP_AFTER_DEV (unset; 1 = exit
# after the development seeds so that the freeze can be gated on a smoke run).
# Smoke chain: CAMPAIGN_NAME=iteration4_smoke2 CAMPAIGN_SCALE=0.25 SPARSE_SEEDS="0 1" DENSE_SEEDS="200 201" ./campaign4.sh
set -euo pipefail
: "${LYALENSER_DATA:?export LYALENSER_DATA first}"
NAME=${CAMPAIGN_NAME:-iteration4}
SCALE=${CAMPAIGN_SCALE:-1}
SPARSE=${SPARSE_SEEDS:-$(seq 0 59)}
DENSE=${DENSE_SEEDS:-$(seq 200 209)}
ROOT=$LYALENSER_DATA/mocks/$NAME
LOGS=$LYALENSER_DATA/condor_logs/$NAME
REPORT=$(cd "$(dirname "$0")/.." && pwd)/report
[ "$NAME" = iteration4 ] || REPORT=$REPORT/$NAME
# Memory requests: the pool is packed, slots with >= 48 GB free are rare (jobs asking for 48 GB idled for
# 9 h while 32 GB ones matched in minutes). Scale-1 peaks: dev/sparse 23.5 GB, dense ~28 GB (round 3).
if [ "$SCALE" = 1 ]; then MEM_SPARSE="32 GB"; MEM_DENSE=${MEM_DENSE:-"36 GB"}; CPUS_DENSE=${CPUS_DENSE:-8}
else MEM_SPARSE="16 GB"; MEM_DENSE="16 GB"; CPUS_DENSE=8; fi
mkdir -p "$LOGS" "$ROOT"
cd "$(dirname "$0")"

submit() {  # submit <tag> <cpus> <mem> <disk> <args...>
  local tag=$1 cpus=$2 mem=$3 disk=$4; shift 4
  if [ -f "$LOGS/$tag.log" ] && grep -q 'return value 0' "$LOGS/$tag.log"; then echo "skip $tag (done)"; return; fi
  rm -f "$LOGS/$tag.log"
  condor_submit phase.sub -a "campaign=$NAME" -a "tag=$tag" -a "cpus=$cpus" -a "mem=$mem" -a "disk=$disk" \
    -a "args=$* --scale $SCALE --mock-root $ROOT" | tail -1 | sed "s/^/$tag: /"
}
wait_for() {  # wait_for <tag>...; fails if any job did not return 0
  local t rc=0
  for t in "$@"; do
    until [ -f "$LOGS/$t.log" ]; do sleep 30; done
    condor_wait -wait 172800 "$LOGS/$t.log" > /dev/null 2>&1 || true
    if grep -q 'return value 0' "$LOGS/$t.log"; then echo "$(date -u +%FT%TZ) $t ok"; else echo "$(date -u +%FT%TZ) $t FAILED"; tail -3 "$LOGS/$t.err"; rc=1; fi
  done
  return $rc
}

echo "$(date -u +%FT%TZ) $NAME scale $SCALE: development seeds"
for s in 100 101 102 103 104; do submit dev$s 8 "$MEM_SPARSE" "20 GB" --phase dev-seed --seed $s; done
wait_for dev100 dev101 dev102 dev103 dev104
[ -z "${STOP_AFTER_DEV:-}" ] || { echo "$(date -u +%FT%TZ) stopping after development seeds (STOP_AFTER_DEV)"; exit 0; }

echo "$(date -u +%FT%TZ) freeze"
submit freeze 8 "$MEM_SPARSE" "20 GB" --phase freeze
wait_for freeze

echo "$(date -u +%FT%TZ) seeds"
for s in $SPARSE; do submit sparse$s 8 "$MEM_SPARSE" "20 GB" --phase seed --seed $s --variant sparse; done
for s in $DENSE; do submit dense$s $CPUS_DENSE "$MEM_DENSE" "40 GB" --phase seed --seed $s --variant dense; done
wait_for sparse0
for n in numerical injection flags random benchmark; do submit ctl_$n 8 "$MEM_DENSE" "40 GB" --phase control --name $n; done
wait_for $(for s in $SPARSE; do echo sparse$s; done) $(for s in $DENSE; do echo dense$s; done) \
         ctl_numerical ctl_injection ctl_flags ctl_random ctl_benchmark

echo "$(date -u +%FT%TZ) collect"
submit collect 8 "$MEM_DENSE" "40 GB" --phase collect --output "$REPORT"
wait_for collect
echo "$(date -u +%FT%TZ) CAMPAIGN $NAME DONE"
