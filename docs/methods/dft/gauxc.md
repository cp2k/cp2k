# GauXC

[GauXC](https://github.com/wavefunction91/GauXC) provides an external exchange-correlation (XC)
integrator for Quickstep. It can evaluate selected conventional functionals and GauXC OneDFT models,
including SKALA models, through the [GAUXC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.GAUXC)
section.

The CP2K interface currently provides two distinct paths:

- The default molecular-quadrature path evaluates the XC contribution through GauXC's atom-centered
  molecular grid. It is intended primarily for isolated calculations.
- The experimental native-grid path (`NATIVE_GRID T`) evaluates SKALA TorchScript models from the
  CP2K GPW real-space grid. It has a different support scope and should be treated separately from
  the molecular-quadrature path.

## Building CP2K with GauXC

Build CP2K with `-DCP2K_USE_GAUXC=ON`. The CP2K toolchain can install GauXC with
`--with-gauxc=install`; this also installs the LibTorch and model resources used by OneDFT/SKALA. An
MPI-enabled CP2K binary requires a GauXC installation with MPI support. See
[](../../technologies/libraries) for the build dependencies and LibTorch compatibility notes.

## Basic Input

GauXC is selected as the only XC functional in the
[XC_FUNCTIONAL](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL) section. The following fragment selects
a conventional functional evaluated by GauXC:

```text
&XC
  &XC_FUNCTIONAL
    &GAUXC
      FUNCTIONAL PBE
    &END GAUXC
  &END XC_FUNCTIONAL
&END XC
```

A non-`NONE` [MODEL](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.GAUXC.MODEL) selects a GauXC OneDFT
model. The model can be a `.fun` file or an installed model name. The underlying
[FUNCTIONAL](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.GAUXC.FUNCTIONAL) is optional in this case
and defaults to `PBE`:

```text
&XC
  &XC_FUNCTIONAL
    &GAUXC
      MODEL path/to/model.fun
    &END GAUXC
  &END XC_FUNCTIONAL
&END XC
```

`MODEL SKALA` obtains the model path from the `GAUXC_SKALA_MODEL` environment variable. The
molecular-quadrature OneDFT/SKALA path defaults to `GRID SUPERFINE` and `PRUNING_SCHEME UNPRUNED`
unless these settings are provided explicitly. These settings are recommended for force checks;
coarser grids are accuracy settings and should be converged for the target calculation.

## Molecular-Quadrature Path

The default path uses GauXC's molecular quadrature. It is the established interface for isolated
Quickstep calculations, including energy, XC potential, and nuclear-gradient calculations within its
supported scope.

### Periodic Reference Calculations

Periodic inputs are not a compact periodic GauXC implementation. To use the molecular path in a
periodic input as an isolated-cell reference calculation, explicitly set
[PERIODIC_REFERENCE](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.GAUXC.PERIODIC_REFERENCE) to `T`:

```text
&GAUXC
  PERIODIC_REFERENCE T
&END GAUXC
```

This reference path is restricted to all of the following:

- `PERIODIC XYZ`;
- Gamma-point calculations with one AO image;
- `METHOD GPW` with GTH pseudopotentials.

It uses a molecular quadrature and must not be used to validate compact periodic materials. Periodic
neighbor-cell AO blocks, k-points, compact-cell quadrature, and periodic stress tensors require a
dedicated periodic GauXC interface.

### OneDFT and SKALA Runtime Controls

[ONEDFT_ATOM_CHUNK_SIZE](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.GAUXC.ONEDFT_ATOM_CHUNK_SIZE)
controls atom-blocked Torch inference. A positive value selects the number of atoms per block, zero
disables chunking, and the default lets GauXC or the `GAUXC_ONEDFT_ATOM_CHUNK_SIZE` environment
variable choose the policy.

For MPI calculations,
[SKALA_RUNTIME](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.GAUXC.SKALA_RUNTIME) controls the
communicator used for OneDFT/SKALA energy and potential evaluation. `AUTO` uses the force-evaluation
communicator for closed-shell calculations and a replicated rank-local runtime for open-shell
calculations. The corresponding
[ONEDFT_GRADIENT_RUNTIME](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.GAUXC.ONEDFT_GRADIENT_RUNTIME)
setting defaults to a conservative replicated runtime for nuclear gradients. Select `MPI` only with
a GauXC installation that supports distributed OneDFT gradients.

### GAPW and Other Restrictions

- Conventional GauXC with `METHOD GAPW` requires all-electron potentials. With pseudopotentials,
  `METHOD GAPW` is available only for OneDFT/SKALA-style models and evaluates the molecular
  AO/valence-density XC term directly.
- OneDFT/SKALA with `METHOD GAPW` and GTH/ECP pseudopotentials currently supports energies only;
  forces and molecular virials are unavailable. CP2K's GAPW one-center XC correction is not used on
  this path.
- `METHOD GAPW_XC`, NLCC pseudopotentials with OneDFT/SKALA, and non-local `VDW_POTENTIAL`
  corrections are not supported by the molecular GauXC path.
- Higher-XC-derivative response and kernel properties are not available through GauXC. Real-time
  propagation is also unsupported.

`MOLECULAR_VIRIAL` is a finite-system diagnostic constructed from GauXC nuclear gradients; it is not
a periodic stress tensor.

## Experimental Native-Grid SKALA Path

`NATIVE_GRID T` bypasses GauXC's molecular quadrature and evaluates a SKALA TorchScript model from
CP2K's GPW real-space grid. It is an experimental path for energy, XC potential, and
nuclear-gradient/stress calculations with one GAUXC functional. Unlike the molecular path, it can
cover selected isolated and periodic GPW/GAPW calculations, including k-point density matrices.

A minimal native-grid input is:

```text
&XC
  &XC_FUNCTIONAL
    &GAUXC
      MODEL SKALA
      NATIVE_GRID T
      NATIVE_GRID_DIAGNOSTICS T
    &END GAUXC
  &END XC_FUNCTIONAL
&END XC
```

`NATIVE_GRID_DIAGNOSTICS T` prints the electron count, spin moment, and summed grid weights of the
feature block supplied to Torch. This is useful when validating a model or a periodic setup.

### CUDA and MPI

The CPU TorchScript path is limited to non-k-point calculations. K-point calculations require a
CUDA-capable LibTorch installation and:

```text
&GAUXC
  NATIVE_GRID T
  NATIVE_GRID_USE_CUDA T
  NATIVE_GRID_CUDA_DEVICE -1
&END GAUXC
```

A negative `NATIVE_GRID_CUDA_DEVICE` assigns the MPI-local rank to a visible CUDA device. A
non-negative value selects that visible device explicitly. CUDA TorchScript exports containing
hard-coded device constants may still require a rank-local `CUDA_VISIBLE_DEVICES` setting.

`NATIVE_GRID_ATOM_CHUNKS T` distributes the model evaluation in atom blocks for MPI calculations and
can reduce peak CUDA memory. `NATIVE_GRID_ATOM_CHUNK_MAX_ROWS` further limits the number of
atom-grid rows handled by one Torch call when needed.

### Current Scope

The native-grid path supports regular-grid GPW/GTH calculations and selected GAPW calculations. For
GAPW with GTH/ECP pseudopotentials, only regular-grid feature kinds are supported; PAW and
one-center GAPW pseudopotential contributions are not implemented. `METHOD GAPW_XC`, ROKS, ADMM, and
non-k-point multiple-image calculations are also outside the current scope.

Because this is an experimental interface, validate energies, forces, and stresses for the chosen
model and system before using it for production calculations.

## Troubleshooting

- `CP2K_GAUXC_STATUS_STDERR=1` mirrors GauXC status messages to standard error. This can be useful
  when a launcher or CI system does not retain the CP2K output file after an external-library
  failure.
- TorchScript models require a LibTorch installation compatible with CP2K's BLAS and OpenMP runtime.
  Pre-built LibTorch bundles can conflict with a CP2K build using oneMKL; use a compatible BLAS
  stack or build LibTorch against the selected dynamic stack when this occurs.
- [OUTPUT_PATH](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.GAUXC.OUTPUT_PATH) writes GauXC molecule
  and basis-set diagnostics to an existing directory. It requires GauXC to have been built with HDF5
  support.

## See Also

- [](../../technologies/libraries)
- [](gpw)
- [](gapw)
