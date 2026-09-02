# NNP GPU verification runbook

The NNP GPU backend can only be built and exercised on a machine with a GPU and the matching
toolkit. This runbook is the ordered check to run there before relying on the backend on a new
architecture. Every step records the binary hash and the git commit, and every rebuild is a clean
rebuild: an incremental build can keep a stale object across a struct-layout or PARAMETER change.

The scripts alongside this file take their inputs from environment variables so they carry no
machine-specific paths. The `inputs/` directory holds the one CP2K input they drive: it lives here
rather than under `tests/` because it is parameterised through preprocessor variables the regression
harness cannot set, and only these scripts run it.

## 1. Build

```
GPU_ARCH=80 tools/nnp_gpu/build_gpu.sh       # A100 / A30; 120 for RTX 5090
```

The build is a stock CMake configure (`-DCP2K_USE_ACCEL=CUDA`, plus `-DCP2K_USE_MPI=ON` so step 3
can compare the code paths at several rank counts), which is the point: the backend must work from
the configure a user would run, with no extra preprocessor definitions.

Confirm the printed `cp2kflags` line contains `nnp_gpu` and note the binary `sha256`. Export the
binary for the remaining steps:

```
export CP2K_BIN=<path printed by build_gpu.sh>
```

## 2. Initialization check

Run any NNP input with `USE_GPU ON` for zero MD steps. The banner must read
`NNP| GPU offload: enabled (DOUBLE precision)`; a broken device setup aborts the run with a
descriptive message instead. The real kernel validation is step 3, which compares the CPU and GPU
paths on the same input.

## 3. CPU vs GPU agreement (double precision)

A ready-made input is `tools/nnp_gpu/inputs/H2O-64_ab.inp` (the H2O-64 committee regtest with
`USE_GPU ${GPU-OFF}` already in the `&NNP` section, where CP2K's preprocessor default separator is a
plain `-`, not the shell's `:-`). Every input body under `tests/NNP/regtest-1` carries the same
`${GPU-OFF}` variable, so other sizes (one water, H2O-256) can be driven the same way from that
directory. Run from the repository root, at 1, 2 and 4 MPI ranks, with `CP2K_DATA_DIR=<repo>/data`
(the input's coordinate path is relative to the working directory):

```
for np in 1 2 4; do
  tools/nnp_gpu/verify_cpu_gpu_agreement.sh tools/nnp_gpu/inputs/H2O-64_ab.inp "$np"
done
```

The double-precision GPU path must agree with the CPU path within the tolerance (`TOL`, default
`1e-9`) on the step-zero energy, the final energy and the conserved quantity. This is the primary
correctness gate.

## 4. Committee and biased committee

Repeat step 3 with a committee input and a biased-committee input. The biased run is the sensitive
check: the bias energy depends on the committee spread, so a per-member error that averages out of
the mean still shows up here.

## 5. Restart

Run a short MD, restart it, and confirm the restarted trajectory continues the uninterrupted one.
The device state must rebuild cleanly from the restart.

## 6. Stress and NPT

The force-scatter kernels accumulate the analytic virial on the device, so a stress-requesting
ensemble runs entirely on the GPU. Three checks:

- `NPT_F` with an analytic virial: compare all nine stress components, per step, against a CPU run
  of the same input.
- `NPT_I` with a numeric virial: the finite-difference displacement evaluations are energy-only, so
  this exercises the GPU energy path under cell perturbations.
- Analytic against finite-difference stress at step zero on the GPU path, which catches a virial
  that is self-consistent but wrong.

If the virial does not come out right on a new architecture, the way to keep going is `USE_GPU OFF`
for the affected runs, not a build switch: the device virial is not separately gated.

When the NNP runs inside `METHOD MIXED`, also run

```
tools/nnp_gpu/check_virial_requirement.sh
```

It asserts that the `STRESS_TENSOR` keyword of the NNP sub-force-eval is load-bearing: dropping it
shifts an NVT pressure by tens of kilobar with no diagnostic while NPT aborts. The script is
CPU-only and needs no GPU.

## 7. Reproducibility

Run the GPU twin from inside `tests/NNP/regtest-gpu`, where its `@INCLUDE` body and the relative
coordinate path resolve (the regtest-1 stub would exercise the deterministic CPU path and prove
nothing about device scheduling):

```
(cd tests/NNP/regtest-gpu && <repo>/tools/nnp_gpu/check_reproducibility.sh H2O-64_C-NNP_MD-gpu.inp 1)
```

Two runs of the same binary must agree to within the ordering of the floating-point reductions: the
accumulation uses `atomicAdd` on `double`, and the order in which blocks reach a given accumulator
follows the device's scheduling, which is not fixed between runs. Expect exact agreement on most
steps and one or two units in the last place on the rest, which is what the default `1e-14`
tolerance allows. A single force evaluation is normally exactly reproducible; the last-place
differences appear once a trajectory accumulates steps. Anything larger is a bug, not noise.

## 8. Mixed precision

Repeat steps 3, 5 and 6 with `GPU_PRECISION MIXED`. The single-precision angular symmetry functions
perturb the result well above the double-precision tolerance, so this step checks physics gates
rather than agreement: energy `1e-5`, forces `1e-4`, stress `1e-4`. Record the observed deviations,
which pin the mixed-precision regression tolerance and the accuracy statement in the documentation.
Step 7 does not apply: reproducibility is asserted for the double-precision path.

## 9. Performance sanity

Time a step on H2O-256 and a larger system (a few thousand atoms) on the CPU and GPU paths and
confirm the GPU is at least as fast as the CPU on the large system. Record the crossover size below
which the CPU path wins.

## 10. Evidence

Collect the logs, binary hashes and the summary table into a dated report and paste the summary into
the pull-request description.
