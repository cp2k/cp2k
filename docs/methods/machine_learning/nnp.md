# Neural Network Potentials

CP2K supports Behler-Parrinello high-dimensional neural network potentials (HDNNPs) with
atom-centred symmetry functions (ACSF) as descriptors. NNP can drive any CP2K run that goes through
`FORCE_EVAL`, from single-point energies through geometry optimisation, molecular dynamics, biased
MD via free-energy methods, and committee-of-models extrapolation diagnostics. The same NNP path
also backs the helium-solvent interaction in path-integral runs.

The implementation reads networks trained with the [n2p2](https://github.com/CompPhysVienna/n2p2)
format (`input.nn`, `scaling.data`, `weights.<element>.data`).

## Input

The NNP method is selected with `METHOD NNP` inside `FORCE_EVAL`. The network files and one or more
model definitions are configured under `&NNP`. The complete option list is generated from the source
and is linked from the input reference manual.

A minimal committee-NNP MD input looks like

```
&FORCE_EVAL
  METHOD NNP
  &NNP
    NNP_INPUT_FILE_NAME nnp-1/input.nn
    SCALE_FILE_NAME     nnp-1/scaling.data
    &MODEL
      WEIGHTS nnp-1/weights
    &END MODEL
    &MODEL
      WEIGHTS nnp-2/weights
    &END MODEL
    ! ... up to 8 typical for committee error bars ...
  &END NNP
  &SUBSYS
    &CELL
      ABC [angstrom] 12.42 12.42 12.42
    &END CELL
    &COORD
      ! ...
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
```

Worked examples covering NVT, NPT, biased MD, and restarts live under
[`tests/NNP/regtest-1/`](https://github.com/cp2k/cp2k/tree/master/tests/NNP/regtest-1). The
path-integral helium-solute coupling is exercised by
[`tests/Pimd/regtest-2/water_in_helium_nnp.inp`](https://github.com/cp2k/cp2k/tree/master/tests/Pimd/regtest-2).

## Tuning

The default ACSF spline grid (`RAD_SPLINE_N 8192`) is calibrated for radial cutoffs around 12 bohr.
Models with substantially larger cutoffs or stricter accuracy targets may want a larger value. The
cubic-Hermite value residual scales as O(1/n^4) and the force (derivative) residual as O(1/n^3); at
the default grid both stay near machine precision for a 12 bohr cutoff (value error ~1e-14, force
error ~1e-10), with the force term the binding constraint as the cutoff grows.

The cell-list neighbour search uses a Verlet skin so the chain is only rebuilt when an atom drifts
more than `skin/2` between force evaluations. By default the skin auto-selects
`MIN(0.5 bohr, 0.1 * cutoff)`. Long stable trajectories can raise it to reduce the rebuild rate at
the cost of a larger per-atom neighbour list:

```
&NNP
  ! ...
  VERLET_SKIN [bohr] 1.0
&END NNP
```

This is the analogue of LAMMPS's `neighbor <skin> bin` command. A negative value (the default)
selects the automatic heuristic; useful upper bounds sit near half the smallest perpendicular cell
width. See the input reference manual for the full list of NNP tuning keywords.

## Parallelism

NNP supports both MPI and OpenMP parallelism. MPI distributes the atoms across ranks, and OpenMP
parallelises the per-atom descriptor and force loops inside each rank.

## GPU acceleration

On an accelerated build the NNP force evaluation can run on the GPU: the neighbour walk,
symmetry-function descriptors, network evaluation and force assembly all execute on the device, with
one host-to-device transfer of positions and one device-to-host transfer of forces per step. It is
controlled per force environment by `USE_GPU`:

```
&NNP
  ! ...
  USE_GPU AUTO         ! AUTO (default) | ON | OFF
  GPU_PRECISION DOUBLE ! DOUBLE (default) | MIXED
&END NNP
```

`AUTO` uses a GPU when the build supports it, every rank sees a device, the model fits the
compiled-in device-table capacities, and the system is large enough for the device to pay off
(currently a thousand atoms per rank); otherwise it runs on the CPU, so the same input runs
unchanged on both. `ON` forces the GPU path and aborts if no device is available. `OFF` forces the
CPU path. Ranks use the device CP2K selected for the process at startup, so one rank per GPU is
recommended. The backend can be disabled at build time with `-DCP2K_ENABLE_NNP_GPU=OFF`.

The backend is written against CP2K's offload layer, so it builds for both CUDA
(`-DCP2K_USE_ACCEL=CUDA`) and HIP (`-DCP2K_USE_ACCEL=HIP`). It handles energies, forces and (for
stress-requesting ensembles such as NPT and variable cell) the analytic virial, which the
force-scatter kernels accumulate on the device alongside the forces. The path-integral
helium-solvent coupling and energy-only evaluations always run on the CPU.

`GPU_PRECISION DOUBLE` (the default) is designed to reproduce the CPU trajectory: the regression
tests hold the GPU and CPU energies (and, on the NPT tests, the pressure) to the documented
tolerances, and the manual checks in `tools/nnp_gpu` extend that to trajectories at several MPI-rank
counts and to the per-component stress. `MIXED` evaluates the angular symmetry functions in single
precision, keeping the radial functions, the network and every accumulator that sums over neighbours
in double precision. It is faster on GPUs whose single-precision throughput exceeds their
double-precision throughput, at the cost of a small perturbation of the results; see the testing
status below for how to measure that perturbation for your own model.

### Combining the NNP with a classical baseline

An HDNNP is often summed with a long-range electrostatic term, which `&FORCE_EVAL METHOD MIXED` does
by pairing an `&NNP` sub-force-eval with a `FIST` one. Two things are worth knowing.

**Every sub-force-eval needs its own `STRESS_TENSOR`.** Stress availability is decided per force
environment from that environment's own `STRESS_TENSOR` keyword and is not inherited from the
`MIXED` level. If the `MIXED` section requests an analytic stress but the `&NNP` sub-force-eval does
not, the NNP contributes a zero virial and CP2K still prints a complete-looking stress tensor built
from the classical term alone. On a 64-water box that omission moves the pressure by around ten
kilobar (measured `7.5e3` to `3.4e4` bar, trajectory dependent). Constant-pressure ensembles abort
on it, but NVT and NVE do not, so a stress-autocorrelation viscosity calculation will not notice.
`tests/NNP/regtest-gpu` and `tests/NNP/regtest-1` both carry a `MIXED` NNP-plus-FIST case that pins
this.

**Only one NNP environment per process may use the GPU.** The device state is per process and its
once-per-run uploads are keyed by network shape, so a second GPU-enabled `&NNP` section would share
the first one's weights while overwriting its descriptor tables. Under `USE_GPU AUTO` the second and
any further sections quietly fall back to the CPU path and say so in the output; under an explicit
`USE_GPU ON` the run aborts instead of silently ignoring what was asked for.

### Limitations

- Energy-only evaluations fall back to the CPU path (the device path is entered only when forces are
  requested).
- The GPU path uses double-precision atomic accumulation, so a result may differ from the CPU by a
  small amount that depends on the accumulation order; the regression tests allow a relative `1e-9`.
  The same applies between two runs of one binary: the order in which blocks reach an accumulator
  follows the device's scheduling, so a trajectory can differ by a unit or two in the last place.
  Runs are not bit-for-bit repeatable, and a simulation that needs exact repeatability should use
  the CPU path.

### Requirements and testing status

The device kernels use double-precision atomic addition (NVIDIA compute capability 6.0 / Pascal or
newer) and, on CUDA, opt into the large dynamic shared-memory carveout on compute capability 7.0 /
Volta or newer. The symmetry-function kernels fall back to a two-pass scheme without the carveout,
so they run on older devices; the network kernel does not, so a network whose per-block working set
exceeds the device's shared-memory limit aborts at the first force evaluation with the size it
needed on the error stream. Models larger than the compiled-in device-table capacities are caught
earlier, at initialization, where `AUTO` falls back to the CPU path and `ON` aborts naming the
limit. The GPU architecture is selected the usual way, with `-DCMAKE_CUDA_ARCHITECTURES=<number>`
(or the legacy `-DCP2K_WITH_GPU=<arch>`). A recent CUDA toolkit (11 or newer) or ROCm is expected.

The backend has been exercised on NVIDIA Ampere (A30, compute capability 8.0) from a stock
`-DCP2K_USE_ACCEL=CUDA` build with CUDA 12.9. On that machine, for a 64-water box with an
eight-member committee:

- the double-precision path reproduces the CPU step-zero energy exactly and the ten-step trajectory
  to a relative `1.2e-16`, at 1, 2 and 4 MPI ranks;
- the analytic virial agrees with the CPU virial to better than a relative `5e-11` over all 231
  stress components of a ten-step NPT run (three measurements of that worst component span `3.1e-11`
  to `4.7e-11`; the virial sums more device-side atomic contributions than any other output, so it
  is the figure most exposed to accumulation order);
- a 1000-step NVE run drifts by `8.1e-07` Ha/atom, the same as the CPU path to three digits;
- `GPU_PRECISION MIXED` shifts the step-zero energy by a relative `4.5e-11`, four orders inside the
  `1e-6` its regression tolerance allows;
- `NNP/regtest-gpu` and `NNP/regtest-1` both pass in full.

The mixed-precision deviation depends on the model and the system, so measure it for yours with
`tools/nnp_gpu/verify_cpu_gpu_agreement.sh` before trusting `MIXED` for production.

It has not been run on any other architecture. The AMD/HIP path compiles (verified with ROCm 6.2
hipcc for gfx90a through the CMake `-DCP2K_USE_ACCEL=HIP` configure) but has never run on AMD
hardware, and no CUDA device older than Ampere has been tried. Two things are worth knowing before
trying AMD: the kernels are written around a 32-thread block, so half the lanes of a 64-wide
wavefront sit idle throughout, and the large shared-memory opt-in is CUDA-only, so the angular
kernels stay at the 48 KB budget there. Neither costs accuracy, but both cost speed, and the block
geometry is the first thing to revisit when tuning for CDNA.

Before relying on the backend on a new architecture, build it there and work through
[`tools/nnp_gpu/README.md`](https://github.com/cp2k/cp2k/tree/master/tools/nnp_gpu), which compares
the GPU and CPU paths across MPI-rank counts, checks the analytic virial against finite differences,
bounds the run-to-run variation, and measures the speed-up.

## References

For background on the method and the existing implementation, the following links might be helpful

- <https://www.cp2k.org/tools:aml>
- <https://doi.org/10.1063/5.0160326>
- [](#Behler2007)
- [](#Behler2011)
- [](#Schran2020)
- [](#Schran2020b)
