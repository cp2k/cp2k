#!/usr/bin/env bash
# Run one NNP input twice with the SAME binary and the SAME settings and compare
# the energy trajectories step by step.
#
# The double-precision GPU path is reproducible to within the ordering of its
# floating-point reductions, not bit-identically: the descriptor and force
# accumulation uses atomicAdd on double, and the order in which blocks reach a
# given accumulator depends on how the device schedules them, which varies
# between runs. The expected signature is therefore agreement to a couple of
# units in the last place on some steps and exact agreement on the rest. A
# deviation far above that is a bug, so this reports the worst relative
# difference and gates on it.
#
# Usage:  check_reproducibility.sh <input.inp> [n_ranks]
# Env:    CP2K_BIN   path to the cp2k binary (required)
#         MPIEXEC    launcher (default: mpirun)
#         TOL        relative tolerance (default: 1e-14, a few ULP of a
#                    double; well below any physically meaningful change)
set -euo pipefail

input="${1:?usage: check_reproducibility.sh <input.inp> [n_ranks]}"
n_ranks="${2:-1}"
CP2K_BIN="${CP2K_BIN:?set CP2K_BIN to the cp2k binary}"
MPIEXEC="${MPIEXEC:-mpirun}"
TOL="${TOL:-1e-14}"

workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT

for run in a b; do
  "$MPIEXEC" -np "$n_ranks" "$CP2K_BIN" -i "$input" -o "$workdir/$run.out" > /dev/null 2>&1 || {
    echo "ERROR: run $run failed" >&2
    exit 1
  }
  grep 'ENERGY| Total FORCE_EVAL' "$workdir/$run.out" | awk '{print $NF}' > "$workdir/$run.ener"
done

echo "== reproducibility: $(basename "$input") (np=$n_ranks) =="

if [[ "$(wc -l < "$workdir/a.ener")" -ne "$(wc -l < "$workdir/b.ener")" ]]; then
  echo "  the two runs printed a different number of energies: FAIL"
  exit 1
fi

if diff -q "$workdir/a.ener" "$workdir/b.ener" > /dev/null; then
  echo "  identical energy trajectory across two runs (bit-for-bit): OK"
  exit 0
fi

paste "$workdir/a.ener" "$workdir/b.ener" |
  awk -v tol="$TOL" '
    { d = $1 - $2; if (d < 0) d = -d
      r = ($1 != 0) ? d / $1 : d; if (r < 0) r = -r
      if (r > worst) { worst = r; step = NR - 1 }
      n++ }
    END {
      printf "  %d steps compared, worst relative difference %.2e at step %d\n", n, worst, step
      if (worst > tol) {
        printf "  above the %s tolerance: FAIL\n", tol; exit 1
      }
      printf "  within the %s tolerance (reduction-order noise): OK\n", tol
    }'
