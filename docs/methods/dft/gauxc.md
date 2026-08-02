# GauXC

[GauXC](https://github.com/wavefunction91/GauXC) provides an external exchange-correlation (XC)
integrator for Quickstep. It can evaluate selected conventional functionals and GauXC Skala models
through the [GAUXC](#CP2K_INPUT.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.GAUXC) section.

The CP2K interface currently provides two distinct paths:

- The default molecular-quadrature path evaluates the XC contribution through GauXC's atom-centered
  molecular grid. It is intended primarily for isolated calculations.
- The experimental native-grid path (`NATIVE_GRID T`) evaluates SKALA TorchScript models from the
  CP2K GPW real-space grid. It has a different support scope and should be treated separately from
  the molecular-quadrature path.

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

A non-`NONE` model selects a GauXC Skala model. The model can be a `.fun` file or an installed model
name. The underlying functional is optional in this case and defaults to `PBE`:

```text
&XC
  &XC_FUNCTIONAL
    &GAUXC
      MODEL path/to/model.fun
    &END GAUXC
  &END XC_FUNCTIONAL
&END XC
```

The `.fun` format and available model checkpoints are defined by GauXC rather than CP2K. `.fun`
files use TorchScript serialization. The official Skala 1.1 checkpoint can be downloaded together
with GauXC in toolchain, or from Hugging Face with the `hf` command provided by the
`huggingface_hub` package:

```bash
hf download microsoft/skala-1.1 skala-1.1.fun --local-dir .
```

Select the downloaded checkpoint with `MODEL ./skala-1.1.fun`. Alternatively, set
`GAUXC_SKALA_MODEL` to its path and use `MODEL SKALA`. See the
[GauXC/Skala model documentation](https://microsoft.github.io/skala/gauxc/c-library.html#download-checkpoint-from-huggingface)
for the model interface and available checkpoints. Obtain `.fun` files only from trusted sources.

The molecular-quadrature Skala path defaults to `GRID SUPERFINE` and `PRUNING_SCHEME UNPRUNED`
unless these settings are provided explicitly. These settings are recommended for force checks;
coarser grids are accuracy settings and should be converged for the target calculation.

## Molecular-Quadrature Path

The default path uses GauXC's molecular quadrature. It is the established interface for isolated
Quickstep calculations, including energy, XC potential, and nuclear-gradient calculations within its
supported scope.

### Periodic Reference Calculations

Periodic inputs are not a compact periodic GauXC implementation. To use the molecular path in a
periodic input as an isolated-cell reference calculation, explicitly set `PERIODIC_REFERENCE` to
`T`:

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

### Skala Runtime Controls

`MODEL_ATOM_CHUNK_SIZE` controls atom-blocked Torch inference. A positive value selects the number
of atoms per block, zero disables chunking, and the default lets GauXC or the
`GAUXC_ONEDFT_ATOM_CHUNK_SIZE` environment variable choose the policy.

For MPI calculations, `SKALA_RUNTIME` controls the communicator used for Skala energy and potential
evaluation. `AUTO` uses the force-evaluation communicator for closed-shell calculations and a
replicated rank-local runtime for open-shell calculations. The corresponding
`MODEL_GRADIENT_RUNTIME` setting defaults to a conservative replicated runtime for nuclear
gradients. Select `MPI` only with a GauXC installation that supports distributed Skala gradients.

### GAPW Density Representations

Conventional GauXC with `METHOD GAPW` requires all-electron potentials. Skala models additionally
support pseudopotential GAPW calculations and distinguish the following density representations:

- All-electron `METHOD GAPW` passes the GAPW AO density matrix to GauXC and retains CP2K's standard
  GAPW one-center XC contribution.
- A GTH/ECP kind with `GPW_TYPE` passes the molecular AO valence-density matrix directly to GauXC.
  This is the regular-grid, GPW-like pseudopotential route; no GAPW one-center XC correction is
  added.
- A GTH/ECP kind without `GPW_TYPE` uses a PAW-like GAPW representation. GauXC evaluates the smooth
  molecular AO term, while CP2K evaluates the Skala model separately for the hard and soft atomic
  densities and adds the one-center hard-minus-soft correction.
- `METHOD GAPW_XC` passes `rho_xc` rather than `rho` to GauXC and combines it with the corresponding
  one-center correction.

The shared one-center implementation reads `NATIVE_GRID_GAPW_DENSITY_PARTITION` even for the
molecular path. Its default, `HARD_MINUS_SOFT`, follows the CP2K GAPW XC construction; `HARD_ONLY`,
`SOFT_ONLY`, and `NONE` are diagnostic choices. The similarly named `NATIVE_GRID_ATOM_PARTITION`
does not affect GauXC molecular quadrature.

For a semilocal functional, the hard-minus-soft term is the usual GAPW one-center construction. A
Skala model also contains non-local descriptors, so evaluating the smooth, hard, and soft terms in
separate model calls is not a mathematical identity for a model evaluated on their combined density:
cross terms can be system dependent. `GPW_TYPE` is therefore the direct pseudopotential route when a
PAW-like one-center correction is not required; PAW-like results should be validated for the target
system.

The one-center term inherits the usual GAPW quadrature controls. Its convergence should be checked
with the kind-dependent `RADIAL_GRID`, `LEBEDEV_GRID`, and `HARD_EXP_RADIUS` settings and with
`GAPW_ACCURATE_XCINT`; established tighter settings should not be reduced for a Skala calculation.

Molecular Skala forces are available for these GAPW and GAPW_XC cases. CP2K currently evaluates the
GauXC molecular XC nuclear gradient for every GAPW method with a conservative central
finite-difference fallback and combines it with the CP2K one-center contribution where applicable.
`MOLECULAR_VIRIAL` is a finite-system diagnostic constructed from these nuclear gradients; it is not
a periodic stress tensor.

NLCC pseudopotentials with Skala and non-local `VDW_POTENTIAL` corrections are not supported by the
molecular GauXC path. Molecular NLCC would require the frozen-core density to enter the Skala
feature definition and its derivatives consistently. Higher-XC-derivative response and kernel
properties are not available through GauXC, and real-time propagation is also unsupported.

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

### CPU, CUDA, and MPI

The native-grid implementation supports both CPU and CUDA TorchScript evaluation, including k-point
calculations. CUDA evaluation is selected explicitly with:

```text
&GAUXC
  NATIVE_GRID T
  NATIVE_GRID_USE_CUDA T
  NATIVE_GRID_CUDA_DEVICE -1
&END GAUXC
```

A negative `NATIVE_GRID_CUDA_DEVICE` assigns the MPI-local rank to a visible CUDA device. A
non-negative value selects that visible device explicitly. CPU k-point calculations require a
compatible LibTorch/BLAS runtime; see Troubleshooting below.

`NATIVE_GRID_ATOM_CHUNKS T` distributes the model evaluation in atom blocks for MPI calculations and
can reduce peak CUDA memory. `NATIVE_GRID_ATOM_CHUNK_MAX_ROWS` further limits the number of
atom-grid rows handled by one Torch call when needed.

### Atom and GAPW Density Partitions

`NATIVE_GRID_ATOM_PARTITION` assigns native-grid rows to atomic feature blocks. `SMOOTH`, the
default, uses a differentiable Becke-like partition. `HARD` assigns each point to its nearest
periodic atom and is intended for legacy energy/VXC checks. Force and stress calculations use the
smooth partition internally because derivatives of the partition weights contribute to the Skala
response.

`NATIVE_GRID_GAPW_DENSITY_PARTITION` is independent of that spatial atom partition. It selects the
one-center density term for PAW-like `METHOD GAPW` and `METHOD GAPW_XC` calculations:
`HARD_MINUS_SOFT` is the default, while `HARD_ONLY`, `SOFT_ONLY`, and `NONE` are diagnostic
variants. Kinds marked with `GPW_TYPE` use only the regular-grid valence-density route and therefore
do not add this one-center correction.

### Current Scope

The native-grid path provides energy, VXC, nuclear forces, and analytical stress for regular-grid
GPW with GTH/ECP pseudopotentials, all-electron GAPW, pseudopotential GAPW with either `GPW_TYPE` or
the PAW-like one-center correction, and `METHOD GAPW_XC`. NLCC is supported for the native
regular-grid GPW representation.

For native-grid NLCC, CP2K adds the frozen-core density to each spin density in both real and
reciprocal space before constructing the Skala features. The model is evaluated once on the combined
valence-plus-core density, retaining density/core cross terms. **NLCC augments $\rho$ and
$\nabla\rho$, while $\tau$ remains valence-only.** This follows CP2K's meta-GGA-like NLCC
convention, but it is not identical to an all-electron frozen-core representation because no core
kinetic-energy density is supplied. Consequently, the effect of NLCC on agreement with an
all-electron reference must be validated for the target chemistry rather than assumed to be
systematically beneficial.

K-point density matrices use CP2K's standard weights and symmetry reduction. The tested scope
includes inversion-only reduction, full K290 reduction, and SPGLIB reduction for GPW/GTH,
all-electron GAPW, and PAW-like GAPW with GTH/ECP pseudopotentials.

ROKS, ADMM, and non-k-point multiple-image calculations remain outside the current scope.

Because this is an experimental interface, validate energies, forces, and stresses for the chosen
model and system before using it for production calculations.

## Troubleshooting

- `CP2K_GAUXC_STATUS_STDERR=1` mirrors GauXC status messages to standard error. This can be useful
  when a launcher or CI system does not retain the CP2K output file after an external-library
  failure.
- TorchScript models require a LibTorch installation compatible with CP2K's BLAS, ScaLAPACK, and
  OpenMP runtimes. In particular, an OpenBLAS-linked CP2K can fail inside LibTorch's batched matrix
  multiplication. Changing the load order does not repair an incompatible batched-BLAS interface.
  Use one consistent stack for CP2K, ScaLAPACK, and LibTorch instead.
- Do not preload oneMKL into an OpenBLAS-linked CP2K as a workaround: interposed complex BLAS
  symbols can break the ScaLAPACK k-point path. Rebuild CP2K and its numerical dependencies
  consistently against oneMKL, or use a mutually compatible OpenBLAS/LibTorch combination.
- `OUTPUT_PATH` writes GauXC molecule and basis-set diagnostics to an existing directory. It
  requires GauXC to have been built with HDF5 support.

## See Also

- [](../../technologies/libraries)
- [](gpw)
- [](gapw)
