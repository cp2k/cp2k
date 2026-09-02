#!/usr/bin/env bash
# Run one NNP input twice, once on the CPU path and once on the GPU path, and
# compare the step-zero and final energies and the conserved quantity. The
# double-precision GPU path must agree with the CPU path to a small relative
# tolerance; a larger deviation is a bug, not floating-point noise. This is an
# agreement check, not a bit-exactness check: the device accumulates forces with
# atomic addition, so the last bits depend on the order the blocks arrive in and
# neither CPU-versus-GPU nor run-versus-run agreement is exact.
#
# The input must select the NNP method and run a few MD steps. USE_GPU is set on
# the command line via a CP2K preprocessor variable, so the same input file is
# used for both runs; add "USE_GPU ${GPU-OFF}" to the &NNP section of the input
# (CP2K's default separator is a plain '-', not the shell's ':-') and this
# script sets GPU=OFF then GPU=ON. A ready-made input is inputs/H2O-64_ab.inp
# (run with CP2K_DATA_DIR=<repo>/data so the NNP/ weight paths resolve).
#
# Usage:  verify_cpu_gpu_agreement.sh <input.inp> [n_ranks]
# Env:    CP2K_BIN   path to the cp2k binary (required)
#         MPIEXEC    launcher (default: mpirun)
#         COORDS     coordinate file for inputs that read ${COORDS-...}, needed
#                    when the run happens somewhere other than the repository
#                    root (the default path is relative to it)
#         TOL        relative tolerance (default: 1e-9)
#         PREC       GPU_PRECISION for the GPU run (default: DOUBLE). With
#                    MIXED the deviation is a physics gate, not agreement, so
#                    raise TOL to match (see README.md step 8).
set -euo pipefail

input="${1:?usage: verify_cpu_gpu_agreement.sh <input.inp> [n_ranks]}"
n_ranks="${2:-1}"
CP2K_BIN="${CP2K_BIN:?set CP2K_BIN to the cp2k binary}"
MPIEXEC="${MPIEXEC:-mpirun}"
TOL="${TOL:-1e-9}"
PREC="${PREC:-DOUBLE}"

workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT

run() {
  local mode="$1" out="$workdir/$1.out"
  local extra=()
  [ -n "${COORDS:-}" ] && extra=(-E "COORDS=$COORDS")
  "$MPIEXEC" -np "$n_ranks" "$CP2K_BIN" -E "GPU=$mode" -E "PREC=$PREC" "${extra[@]}" -i "$input" -o "$out" > /dev/null 2>&1 || {
    echo "ERROR: $mode run failed; see $out" >&2
    cat "$out" >&2 || true
    exit 1
  }
  # Step-zero and final total energy, plus the final conserved quantity.
  grep 'ENERGY| Total FORCE_EVAL' "$out" | head -1 | awk '{print $NF}'
  grep 'ENERGY| Total FORCE_EVAL' "$out" | tail -1 | awk '{print $NF}'
  grep 'Conserved quantity' "$out" | tail -1 | awk '{print $NF}' || echo "NA"
}

mapfile -t cpu < <(run OFF)
mapfile -t gpu < <(run ON)

# run() executes in a process substitution, so its exit does not stop this
# shell; an aborted run shows up as missing lines instead.
if [ "${#cpu[@]}" -lt 3 ] || [ "${#gpu[@]}" -lt 3 ]; then
  echo "ERROR: a run produced no parsable energies (see the messages above): FAIL" >&2
  exit 1
fi

echo "== CPU/GPU agreement: $(basename "$input") (np=$n_ranks, $PREC) =="
labels=("step-0 energy" "final energy" "conserved qty")
status=0
for i in 0 1 2; do
  a="${cpu[$i]}"
  b="${gpu[$i]}"
  [[ "$a" == "NA" || "$b" == "NA" ]] && {
    printf '  %-16s  (not present)\n' "${labels[$i]}"
    continue
  }
  rel="$(awk -v x="$a" -v y="$b" 'BEGIN{d=(x-y); if(x!=0)d/=x; if(d<0)d=-d; printf "%.3e", d}')"
  ok="$(awk -v r="$rel" -v t="$TOL" 'BEGIN{print (r<=t)?"OK":"FAIL"}')"
  [[ "$ok" == "FAIL" ]] && status=1
  printf '  %-16s cpu=%s gpu=%s  rel=%s  %s\n' "${labels[$i]}" "$a" "$b" "$rel" "$ok"
done
exit "$status"
