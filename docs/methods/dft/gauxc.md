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
files use TorchScript serialization. The official Skala 1.1 Rev1 checkpoints can be downloaded
together with GauXC in the toolchain, or from Hugging Face with the `hf` command provided by the
`huggingface_hub` package:

```bash
hf download microsoft/skala-1.1 skala-1.1-rev1.fun --local-dir .
hf download microsoft/skala-1.1 skala-1.1-rev1-cuda.fun --local-dir .
```

Select a downloaded checkpoint explicitly with `MODEL ./skala-1.1-rev1.fun`. Alternatively, set
`GAUXC_SKALA_MODEL` to the CPU checkpoint and use `MODEL SKALA` for host execution. For
`INT_EXECUTION_SPACE DEVICE`, set `GAUXC_SKALA_CUDA_MODEL` to the CUDA checkpoint. CP2K falls back
to `GAUXC_SKALA_MODEL` for compatibility with existing configurations where that variable already
points to a device-compatible checkpoint; the CPU checkpoint itself is not device-compatible. For
the model interface and available checkpoints, consult
[GauXC/Skala model documentation](https://microsoft.github.io/skala/gauxc/c-library.html#download-checkpoint-from-huggingface)
and obtain `.fun` files only from trusted sources.

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
support pseudopotential GAPW calculations. `PSEUDOPOTENTIAL_GAPW_REPRESENTATION` selects their
density representation explicitly:

- `DIRECT_VALENCE` is the default. It evaluates the model on the direct valence density and is the
  GPW-like route for GTH/ECP kinds, irrespective of the kind's `GPW_TYPE` setting.
- `PAW_ONE_CENTER` reconstructs density, density gradient, and kinetic-energy density as smooth plus
  hard minus soft before Skala evaluation. Forming the nonlinear features after this sum retains all
  15 pairwise cross terms among the six smooth, hard, and soft spin-gradient components. Fourteen
  are additional to the smooth-density baseline. The same construction also retains the non-local
  couplings.
- `PAW_ONE_CENTER_SPLIT` keeps the legacy separately evaluated smooth plus hard-minus-soft energy as
  an explicit diagnostic. This expression is exact for semilocal GAPW XC but is not a mathematical
  identity for the non-local Skala model.
- `CP2K_DEFAULT` recovers the representation implied by the pre-existing `GPW_TYPE`, basis, and
  `FORCE_PAW` kind settings.

The selector applies only to pseudopotential GAPW kinds. All-electron `METHOD GAPW` retains its
all-electron AO density representation, and `METHOD GAPW_XC` selects CP2K's `rho_xc` density before
the corresponding one-center reconstruction.

`NATIVE_GRID_GAPW_DENSITY_PARTITION` controls the legacy one-center diagnostic; its default,
`HARD_MINUS_SOFT`, follows the CP2K GAPW XC construction. `HARD_ONLY`, `SOFT_ONLY`, and `NONE` are
diagnostic choices. The similarly named `NATIVE_GRID_ATOM_PARTITION` controls a native-grid atom
partition and does not affect GauXC molecular quadrature.

The one-center representation inherits the kind-dependent `RADIAL_GRID`, `LEBEDEV_GRID`, and
`HARD_EXP_RADIUS` controls, and established tighter settings should not be reduced for a Skala
calculation. `GAPW_ACCURATE_XCINT` keeps its normal role for classical GAPW XC, `CP2K_DEFAULT`, and
the split diagnostic. The combined `PAW_ONE_CENTER` Skala term instead uses its common reconstructed
quadrature and is independent of the legacy accurate-XCINT hard/soft weights and their derivatives.

The reconstructed expression conserves the electron number at the density-representation level. Its
numerical integral on the finite radial/Lebedev quadrature retains the usual molecular-grid error
and converges with the kind-dependent grid settings. CP2K does not rescale the density to force an
exact numerical integral: a density-dependent rescaling would modify the Skala functional and would
require additional VXC, force, and virial derivatives. `NATIVE_GRID_DIAGNOSTICS T` prints the
atom-composite electron integral for convergence checks.

Molecular Skala forces are available for these GAPW and GAPW_XC cases. The direct molecular GauXC
route evaluates its XC nuclear gradient through the configured GauXC gradient path. The
`PAW_ONE_CENTER` representation instead propagates the Skala feature adjoint analytically through
the CP2K smooth-field interpolation, one-center reconstruction, atom partition, and NLCC center
coordinates. `MOLECULAR_VIRIAL` is a finite-system diagnostic constructed from the nuclear
gradients; it is not a periodic stress tensor.

Direct molecular GauXC evaluation with NLCC pseudopotentials is not supported because GauXC does not
receive the frozen-core density and its derivatives. Molecular pseudopotential GAPW with
`PAW_ONE_CENTER` is a separate, supported route: CP2K evaluates the NLCC density and gradient on the
same atom-centered composite quadrature, adds them to the reconstructed primitive fields before
Skala constructs its features, and differentiates the core and grid-center coordinates analytically.
The kinetic-energy density remains valence-only. Non-local `VDW_POTENTIAL` corrections,
higher-XC-derivative response and kernel properties, and real-time propagation remain unsupported
through GauXC.

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

For the molecular `PAW_ONE_CENTER` representation, each rank assembles only its owned atoms and the
radial hard/soft fields are replicated once per kind for cross-atom overlap. Every active rank
evaluates its complete local atom blocks and the rank-local Skala energies are summed. Feature
adjoints remain on the owning rank, while the smooth plane-wave interpolation adjoint is summed
globally. Skala 1.1 retains an independent atom dimension through its non-local layers, so this
decomposition is exact. With CUDA and automatic device selection, MPI-local ranks execute their atom
blocks on distinct visible GPUs.

### Atom and GAPW Density Partitions

`NATIVE_GRID_ATOM_PARTITION` assigns native-grid rows to atomic feature blocks. `SMOOTH`, the
default, uses a differentiable Becke-like partition. `HARD` assigns each point to its nearest
periodic atom and is intended for legacy energy/VXC checks. Force and stress calculations use the
smooth partition internally because derivatives of the partition weights contribute to the Skala
response.

`NATIVE_GRID_GAPW_DENSITY_PARTITION` is independent of that spatial atom partition. It selects the
one-center density term for PAW-like `METHOD GAPW` and `METHOD GAPW_XC` calculations:
`HARD_MINUS_SOFT` is the default, while `HARD_ONLY`, `SOFT_ONLY`, and `NONE` are diagnostic
variants. `DIRECT_VALENCE` uses the regular-grid valence-density route irrespective of `GPW_TYPE`;
the legacy `CP2K_DEFAULT` representation may infer that route from the kind settings.

For `PAW_ONE_CENTER`, `NATIVE_GRID_GAPW_COMPOSITE_GRID` selects the grid on which the combined
primitive fields are formed. `ATOM_COMPOSITE`, the default, interpolates the smooth field to GAPW
radial/Lebedev grids and adds the hard-minus-soft fields there before constructing the nonlinear
Skala features. In periodic cells, a Becke-like partition over all atom images assigns the outer
energy quadrature. A separate partition over only the target atom's self-images defines its complete
periodic descriptor domain, avoiding a truncation of non-local descriptors at boundaries between
different atoms. `COMMON_GRID` reconstructs the same fields on the regular grid and is retained as a
cutoff-sensitive diagnostic reference. At matched cutoff, complete atom blocks are both
substantially closer to the corresponding non-periodic atom-composite limit and much less expensive
than resolving all hard-minus-soft detail on one global periodic grid. This selector does not change
`DIRECT_VALENCE`, molecular GauXC quadrature, or the all-electron AO density representation. For
pseudopotentials, molecular GauXC `DIRECT_VALENCE` and native `PAW_ONE_CENTER` remain distinct
density representations; agreement between periodic and non-periodic atom-composite calculations
does not imply equality with the direct AO-valence result.

### Current Scope

The native-grid path provides energy, VXC, nuclear forces, and analytical stress for regular-grid
GPW with GTH/ECP pseudopotentials, all-electron GAPW, pseudopotential GAPW with either `GPW_TYPE` or
the PAW-like one-center correction, and `METHOD GAPW_XC`. NLCC is supported for both the native
regular-grid representation and pseudopotential GAPW `PAW_ONE_CENTER`. Periodic `PAW_ONE_CENTER`
uses the atom-centered composite backend by default; its smooth lattice-image partitions,
native-grid interpolation, periodic images, and one-center fields are differentiated consistently
for nuclear forces and strain.

For native-grid NLCC, CP2K adds the frozen-core density to each spin density in both real and
reciprocal space before constructing the Skala features. In the molecular atom-composite route the
same core field is evaluated analytically on the partitioned radial/Lebedev rows. The model is
evaluated once on the combined valence-plus-core primitive fields, retaining density/core and
gradient cross terms. **NLCC augments $\rho$ and $\nabla\rho$, while $\tau$ remains valence-only.**
This follows CP2K's meta-GGA-like NLCC convention, but it is not identical to an all-electron
frozen-core representation because no core kinetic-energy density is supplied. Consequently, the
effect of NLCC on agreement with an all-electron reference must be validated for the target
chemistry rather than assumed to be systematically beneficial.

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
  OpenMP runtimes. Pre-built LibTorch distributions can bundle oneMKL symbols whose grouped
  SGEMM/DGEMM interface is incompatible with the same-named OpenBLAS entry points. For LP64 OpenBLAS
  builds, CP2K expands these grouped operations through the standard CBLAS interface. Other mixed
  BLAS interfaces still require a consistently built numerical stack.
- Do not preload oneMKL into an OpenBLAS-linked CP2K as a workaround: interposed complex BLAS
  symbols can break the ScaLAPACK k-point path.
- `OUTPUT_PATH` writes GauXC molecule and basis-set diagnostics to an existing directory. It
  requires GauXC to have been built with HDF5 support.

## See Also

- [](../../technologies/libraries)
- [](gpw)
- [](gapw)
