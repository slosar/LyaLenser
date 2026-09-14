#!/bin/bash
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --job-name=lyalenser4
# Runs the phase commands listed in a file ("tag|args" per line) concurrently on one node.
# Usage: sbatch -A <account> -t <hh:mm:ss> phase_list.sbatch <listfile> <concurrency> <threads-per-phase>
# Requires LYALENSER_DATA, LYALENSER_REPO, LYALENSER_PYTHON in the submitting environment (sbatch exports it).
set -uo pipefail
LIST=$1; JOBS=$2; THREADS=$3
: "${LYALENSER_DATA:?}" "${LYALENSER_REPO:?}" "${LYALENSER_PYTHON:?}"
LOGS=${CAMPAIGN_LOGS:-$LYALENSER_DATA/slurm_logs/iteration5}; mkdir -p "$LOGS"
export NUMBA_NUM_THREADS=$THREADS OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MPLCONFIGDIR=$SCRATCH/mpl-$USER
cd "$LYALENSER_REPO/code/pipeline"
echo "$(date -u +%FT%TZ) node $(hostname): $(wc -l < "$LIST") phases, $JOBS concurrent, $THREADS threads each"
run_one() {
  local tag=${1%%|*} args=${1#*|}
  echo "host=$(hostname) start=$(date -u +%FT%TZ) threads=$NUMBA_NUM_THREADS args=$args" > "$LOGS/$tag.out"
  if "$LYALENSER_PYTHON" run_mock_validation.py $args >> "$LOGS/$tag.out" 2>&1; then echo "$(date -u +%FT%TZ) $tag ok"
  else echo "$(date -u +%FT%TZ) $tag FAILED"; tail -3 "$LOGS/$tag.out"; return 1; fi
}
export -f run_one; export LOGS LYALENSER_PYTHON NUMBA_NUM_THREADS
xargs -a "$LIST" -d '\n' -P "$JOBS" -I{} bash -c 'run_one "$1"' _ {}
rc=$?
echo "$(date -u +%FT%TZ) list $LIST finished, xargs exit $rc"
exit $rc
