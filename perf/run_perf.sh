#!/bin/bash -e
# Performance test runner for SKALA/GauXC (periodic water cells).
# Usage: OMP_NUM_THREADS=2 NP=2 ./run_perf.sh [logdir]
# Notes:
#  - OMP_NUM_THREADS must be set (Torch intra-op is tied to it); without it the
#    run oversubscribes and appears to hang.
#  - For GPU runs, ensure the offload/MPI rank-to-device mapping is set up.

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source /workspace/install/cp2k_env
export CP2K_DATA_DIR=/workspace/data
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
NP="${NP:-2}"
OUT="${1:-$here/logs}"
mkdir -p "$OUT"

cd "$here"
week=$(date +%s)
for inp in *.inp; do
  [ -e "$inp" ] || continue
  base="${inp%.inp}"
  start=$(date +%s)
  set +e
  mpiexec -n "$NP" --bind-to none /workspace/install/bin/cp2k.psmp \
    -i "$inp" -o "$OUT/$base.out" >/dev/null 2>&1
  ec=$?
  set -e
  end=$(date +%s)
  conv=$(grep -a 'SCF run converged in' "$OUT/$base.out" | tail -1 | sed 's/^ *//')
  abort=$(grep -a -m1 'ABORT' "$OUT/$base.out" | sed 's/^ *//')
  printf '%-28s exit=%-3s wall=%-6ss %s %s\n' "$base" "$ec" "$((end-start))" \
    "${conv:-[no convergence]}" "${abort:-}"
done
echo "total: $(( $(date +%s) - week ))s  (logs in $OUT)"
