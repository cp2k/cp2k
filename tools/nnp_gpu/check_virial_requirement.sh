#!/usr/bin/env bash
# Assert that a MIXED force environment really does require STRESS_TENSOR on the
# NNP sub-force-eval, and that omitting it fails in the two different ways the
# input files claim it does.
#
# pv_availability is derived per force environment and is NOT inherited from the
# MIXED level. Dropping STRESS_TENSOR from the NNP section therefore leaves the
# NNP virial at zero while the MIXED level still prints a complete-looking stress
# tensor holding only the FIST term. That is the dangerous case: under a
# constant-volume ensemble nothing aborts and the pressure is quietly wrong.
# Under constant pressure the omission is caught and the run aborts.
#
# Both halves are checked here, because only checking one leaves the other free
# to regress:
#
#   NVT  minus the NNP STRESS_TENSOR -> completes, pressure moves a long way
#   NPT  minus the NNP STRESS_TENSOR -> aborts, "Virial evaluation not requested"
#
# A control NPT run WITH the stress tensor is included so a failure to start for
# some unrelated reason cannot masquerade as the expected abort.
#
# This gate is about the mixed-environment stress plumbing, not the device, so it
# runs on the CPU path (the reference input carries USE_GPU OFF) and needs no
# GPU. It is in tools/nnp_gpu/ because it is the virial half of the NNP GPU
# verification story; run it wherever the CPU binary is.
#
# Usage:  check_virial_requirement.sh [input]
#         default: tests/NNP/regtest-1/H2O-64_MIXED_C-NNP-FIST_MD-NpT.inc
#         The default is the shared input BODY, not the regtest stub: the stub
#         holds only an @INCLUDE, while the transforms below need the
#         &FORCE_EVAL, ENSEMBLE and BAROSTAT lines themselves.
# Env:    CP2K_BIN       path to the cp2k binary (required)
#         CP2K_DATA_DIR  data directory holding NNP/bulkH2O-jcp2020-cnnp
#                        (default: <repo>/data)
#         MPIEXEC        launcher (default: mpirun)
#         NRANKS         MPI ranks (default: 1)
#         MIN_SHIFT_BAR  smallest NVT pressure shift that counts as "moved"
#                        (default: 1e3; measurements on the committee H2O-64
#                        box range from 7.5e3 to 3.4e4 bar depending on the
#                        trajectory, so this gates on the omission being
#                        obvious rather than on a fixed number)
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
input="${1:-$repo_root/tests/NNP/regtest-1/H2O-64_MIXED_C-NNP-FIST_MD-NpT.inc}"
[ -f "$input" ] || {
  echo "ERROR: no such input: $input" >&2
  exit 2
}
export CP2K_DATA_DIR="${CP2K_DATA_DIR:-$repo_root/data}"

CP2K_BIN="${CP2K_BIN:?set CP2K_BIN to the cp2k binary}"
MPIEXEC="${MPIEXEC:-mpirun}"
NRANKS="${NRANKS:-1}"
MIN_SHIFT_BAR="${MIN_SHIFT_BAR:-1e3}"

workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT

input_dir="$(cd "$(dirname "$input")" && pwd)"

# The reference input reaches its coordinates through a path relative to its own
# directory. We run from a scratch directory, so pin it to an absolute path
# rather than making the caller export COORDS.
abs_coords() {
  awk -v dir="$input_dir" '
    /COORD_FILE_NAME/ {
      p = $2
      if (p !~ /^\//) p = dir "/" p
      printf "  COORD_FILE_NAME %s\n", p
      next
    }
    { print }'
}

# Drop STRESS_TENSOR from the NNP force environment ONLY. Each &FORCE_EVAL block
# is buffered so the decision can be made after METHOD has been seen; the MIXED
# and FIST environments keep theirs.
strip_nnp_stress() {
  awk '
    /^[[:space:]]*&FORCE_EVAL/ { inblk = 1; n = 0; isnnp = 0 }
    inblk {
      buf[++n] = $0
      if ($0 ~ /^[[:space:]]*METHOD[[:space:]]+NNP[[:space:]]*$/) isnnp = 1
      if ($0 ~ /^[[:space:]]*&END[[:space:]]+FORCE_EVAL/) {
        for (i = 1; i <= n; i++)
          if (!(isnnp && buf[i] ~ /^[[:space:]]*STRESS_TENSOR/)) print buf[i]
        inblk = 0
      }
      next
    }
    { print }'
}

# NPT_F -> NVT, and remove the barostat that only belongs to a constant-pressure
# run. The MIXED level keeps its own STRESS_TENSOR, so the pressure is still
# printed, which is exactly what makes the silent-wrong-answer case visible.
# Reached as make_variant's transform argument, never by name, which shellcheck
# cannot see.
# shellcheck disable=SC2329
to_nvt() {
  awk '
    /^[[:space:]]*&BAROSTAT/ { inbaro = 1; next }
    inbaro { if ($0 ~ /^[[:space:]]*&END[[:space:]]+BAROSTAT/) inbaro = 0; next }
    /^[[:space:]]*ENSEMBLE[[:space:]]+NPT_F/ { sub(/NPT_F/, "NVT") }
    { print }'
}

make_variant() { # name, transform-pipeline...
  local name="$1"
  shift
  local f="$workdir/$name.inp"
  if [ "$#" -eq 0 ]; then
    abs_coords < "$input" > "$f"
  else
    abs_coords < "$input" | "$@" > "$f"
  fi
  # Keep each variant's project name distinct so their outputs never collide.
  # An explicit backup suffix keeps sed -i portable across GNU and BSD.
  sed -i.orig "s/^\([[:space:]]*PROJECT[[:space:]][[:space:]]*\).*/\1$name/" "$f"
  rm -f "$f.orig"
  echo "$f"
}

run_variant() { # inp -> writes .out, echoes the exit status
  local inp="$1" name
  name="$(basename "$inp" .inp)"
  local rc=0
  (cd "$workdir" && "$MPIEXEC" -np "$NRANKS" "$CP2K_BIN" \
    -i "$inp" -o "$workdir/$name.out") > /dev/null 2>&1 || rc=$?
  echo "$rc"
}

first_pressure() { # out -> first printed pressure in bar, or NA
  awk '/MD\| Pressure \[bar\]/ { print $(NF); exit }' "$1" 2> /dev/null || true
}

echo "== virial requirement on the MIXED force environment =="
echo "   input:  $(basename "$input")"
echo "   binary: $CP2K_BIN  (np=$NRANKS)"
echo

status=0

# --- leg 1: NVT, both with and without the NNP stress tensor ---------------
nvt_ref="$(make_variant nvt_ref to_nvt)"
nvt_bad="$(make_variant nvt_nostress to_nvt)"
strip_nnp_stress < "$nvt_bad" > "$nvt_bad.tmp" && mv "$nvt_bad.tmp" "$nvt_bad"

rc_ref="$(run_variant "$nvt_ref")"
rc_bad="$(run_variant "$nvt_bad")"

if [ "$rc_ref" -ne 0 ]; then
  echo "  NVT reference run failed (rc=$rc_ref), cannot judge the rest: FAIL"
  sed -n '$p;/ABORT/p' "$workdir/nvt_ref.out" 2> /dev/null | head -5
  exit 1
fi

p_ref="$(first_pressure "$workdir/nvt_ref.out")"
p_bad="$(first_pressure "$workdir/nvt_nostress.out")"

if [ "$rc_bad" -ne 0 ]; then
  echo "  NVT without the NNP stress tensor aborted (rc=$rc_bad)"
  echo "    expected it to complete silently with a wrong pressure: FAIL"
  status=1
elif [ -z "$p_ref" ] || [ -z "$p_bad" ]; then
  echo "  no 'MD| Pressure [bar]' line in one of the NVT runs: FAIL"
  status=1
else
  shift_bar="$(awk -v a="$p_ref" -v b="$p_bad" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%.4e", d}')"
  ok="$(awk -v d="$shift_bar" -v m="$MIN_SHIFT_BAR" 'BEGIN{print (d>=m)?"OK":"FAIL"}')"
  [ "$ok" = FAIL ] && status=1
  printf '  NVT  with stress: %s bar\n' "$p_ref"
  printf '  NVT  without    : %s bar\n' "$p_bad"
  printf '  shift           : %s bar (>= %s required, no abort)  %s\n' \
    "$shift_bar" "$MIN_SHIFT_BAR" "$ok"
fi

# --- leg 2: NPT control, then NPT with the omission ------------------------
npt_ref="$(make_variant npt_ref)"
npt_bad="$(make_variant npt_nostress)"
strip_nnp_stress < "$npt_bad" > "$npt_bad.tmp" && mv "$npt_bad.tmp" "$npt_bad"

rc_npt_ref="$(run_variant "$npt_ref")"
if [ "$rc_npt_ref" -ne 0 ]; then
  echo "  NPT control run WITH the stress tensor failed (rc=$rc_npt_ref): FAIL"
  echo "    the abort below would not be attributable to the omission"
  status=1
else
  echo "  NPT  with stress: completes  OK"
fi

rc_npt_bad="$(run_variant "$npt_bad")"
npt_out="$workdir/npt_nostress.out"
if [ "$rc_npt_bad" -eq 0 ]; then
  echo "  NPT without the NNP stress tensor completed: FAIL"
  echo "    a constant-pressure run must not accept a partial virial"
  status=1
elif grep -qi 'Virial evaluation not requested' "$npt_out" 2> /dev/null; then
  echo "  NPT  without    : aborts with 'Virial evaluation not requested'  OK"
else
  echo "  NPT without the NNP stress tensor aborted (rc=$rc_npt_bad) but not with"
  echo "    the expected message: FAIL"
  grep -i -m3 -E 'ABORT|ERROR' "$npt_out" 2> /dev/null | sed 's/^/      /'
  status=1
fi

echo
[ "$status" -eq 0 ] && echo "PASS: the NNP sub-force-eval's STRESS_TENSOR is load-bearing in both ensembles." ||
  echo "FAIL: see above."
exit "$status"
