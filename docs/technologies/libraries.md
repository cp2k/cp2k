# Libraries

## BLAS and LAPACK (required, base functionality)

BLAS and LAPACK should be installed. Using vendor-provided libraries can make a very significant
difference (up to 100%, e.g., ACML, MKL, ESSL), not all optimized libraries are bug free. Use the
latest versions available, use the interfaces matching your compiler, and download all patches!

- The canonical BLAS and LAPACK can be obtained from the Netlib repository:
  - <http://www.netlib.org/blas/>
  - <http://www.netlib.org/lapack/> and see also
  - <http://www.netlib.org/lapack-dev/>
- Open fast alternatives, include:
  - <http://www.openblas.net>
  - <http://math-atlas.sourceforge.net>
  - <https://www.tacc.utexas.edu/research-development/tacc-software/gotoblas2>

Please note that the BLAS/LAPACK implementation used by CP2K needs to be thread-safe (OpenMP).
Examples are the sequential variant of the Intel MKL, the Cray libsci, the OpenBLAS OpenMP variant
and the reference BLAS/LAPACK packages. Usually the CMake step of CP2K will auto-detect the type of
BLAS and SCALAPACK and then use the right configuration to ensure the code is thread-safe; however,
if you encounter problems when compiling with MKL, try passing
`-DCP2K_BLAS_VENDOR=MKL -DCP2K_SCALAPACK_VENDOR=MKL` to CMake. MKL with multiple OpenMP threads in
CP2K requires that CP2K was compiled with the Intel compiler.

On the Mac, BLAS and LAPACK may be provided by Apple's Accelerate framework. If using this
framework, `-DCP2K_BLAS_VENDOR=Apple` must be passed to CMake to account for some interface
incompatibilities between Accelerate and reference BLAS/LAPACK.

## MPI and ScaLAPACK (required for MPI parallel builds)

MPI (version 3 or later) and SCALAPACK are needed for parallel code. (Use the latest versions
available and download all patches!).

:warning: Note that your MPI installation must match the used Fortran compiler. If your computing
platform does not provide MPI, there are several freely available alternatives:

- MPICH MPI: <https://www.mpich.org/> (may require `-fallow-argument-mismatch` when building with
  GCC 10)
- OpenMPI MPI: <http://www.open-mpi.org/>

For more information of ScaLAPACK, see <http://www.netlib.org/scalapack/>. ScaLAPACK can be part of
ACML (AMD) or cluster MKL (Intel); these libraries are recommended on the corresponding machines if
available.

CP2K assumes that the MPI library implements MPI version 3. Older versions of MPI (e.g., MPI 2.0)
are not supported. CP2K can make use of the `mpi_f08` module. If its use is requested, pass
`-DCP2K_USE_MPI_F08=ON` to CMake.

## FFTW (improved performance of FFTs)

FFTW can be used to improve FFT speed on a wide range of architectures. It is strongly recommended
to install and use FFTW3. The current version of CP2K works with FFTW 3.X (pass
`-DCP2K_USE_FFTW3=ON` to CMake). It can be downloaded from <http://www.fftw.org>.

FFTW is also provided by MKL. If you have MKL but still want to use standalone FFTW3, pass
`-DCP2K_USE_FFTW3_WITH_MKL=ON` to CMake.

:warning: Note that FFTW must know the Fortran compiler you will use in order to install properly
(e.g., `export F77=gfortran` before configure if you intend to use gfortran).

Since CP2K is OpenMP parallelized, the FFTW3 threading library libfftw3_omp (or libfftw3_threads) is
required. Pass `-DCP2K_ENABLE_FFTW3_OPENMP_SUPPORT=ON` or `-DCP2K_ENABLE_FFTW3_THREADS_SUPPORT=ON`
respectivly to CMake.

## LIBINT (enables methods including HF exchange)

- Hartree-Fock exchange requires the LIBINT package to calculate ERI.
- Pass `-DCP2K_USE_LIBINT2=ON` to CMake to enable LIBINT.
- It's always suggested to build LIBINT with toolchain or Spack. If you want to build it yourself,
  download a CP2K-configured LIBINT library from our
  [CP2K server](https://www.cp2k.org/static/downloads/), then build and install LIBINT by following
  the instructions provided [here](https://github.com/cp2k/cp2k/blob/master/tools/libint/README.md).
  Note that using a library configured for higher maximum angular momentum will increase build time
  and binary size of CP2K binary (assuming static linking).
- CP2K is not hardwired to these provided libraries and any other LIBINT library (version >= 2.5.0)
  should be compatible as long as it was compiled with `--enable-eri=1` and default ordering.
- Avoid debugging information (`-g` or `-g2` flag) for compiling LIBINT since this will increase
  library size by a large factor.

## LIBXSMM (improved performance for matrix multiplication)

- A library for matrix operations and deep learning primitives:
  <https://github.com/libxsmm/libxsmm/>.
- Pass `-DCP2K_USE_LIBXSMM=ON` to CMake to enable it.
- LIBXSMM can be used with both CUDA and HIP backends (see [](./accelerators/index.md)).

## LIBXC (wider choice of xc functionals)

- The latest version of LIBXC can be downloaded from <https://gitlab.com/libxc/libxc/-/releases>
- CP2K does not make use of fourth derivates such that LIBXC may be configured with
  `./configure --disable-lxc <other LIBXC configuration flags>`.
- During the installation, the directories `$(LIBXC_DIR)/lib` and `$(LIBXC_DIR)/include` are
  created.
- Pass `-DCP2K_USE_LIBXC=ON` to CMake.

## GauXC (xc integration library)

GauXC can be used to evaluate selected exchange-correlation functionals through an external
integrator.

- Pass `--with-gauxc=install` to the toolchain installer. The toolchain build enables GauXC
  OneDFT/SKALA support and therefore also installs libtorch.
- Pass `-DCP2K_USE_GAUXC=ON` to CMake.
- GauXC in CP2K is primarily an energy, potential, and nuclear-gradient path for isolated QS
  systems. A limited isolated-cell reference path can be enabled with `PERIODIC_REFERENCE T` for
  Gamma-only, single-image `METHOD GPW` calculations with GTH pseudopotentials in periodic
  `PERIODIC XYZ` inputs. This opt-in path keeps molecular validation inputs usable when CP2K is run
  with periodic boundary conditions, but it still uses GauXC's molecular quadrature. It must not be
  used as validation for compact periodic materials. k-points, periodic neighbour-cell AO blocks,
  compact-cell quadrature, GAPW/GAPW_XC, and periodic stress tensors require a dedicated periodic
  GauXC design. OneDFT/SKALA gradients under MPI are evaluated with a replicated single-rank GauXC
  runtime on each CP2K rank because GauXC does not yet provide distributed OneDFT gradients.
- OneDFT/SKALA is selected in the `&GAUXC` subsection with a non-`NONE` `MODEL`, for example a
  `.fun` model file or a GauXC-installed model name. The `FUNCTIONAL` keyword is optional for
  OneDFT/SKALA inputs and defaults to `PBE`; `MODEL SKALA` inputs do not need an explicit
  `FUNCTIONAL PBE` line.
- `ONEDFT_ATOM_CHUNK_SIZE` can be used to control the GauXC OneDFT/SKALA Torch atom blocking from
  CP2K. A positive value requests atom-by-atom chunks of that size, zero disables atom chunking, and
  the default leaves GauXC's model policy or `GAUXC_ONEDFT_ATOM_CHUNK_SIZE` environment setting in
  control.
- `GAUXC_GRADIENT_MPI_RUNTIME=1` is an experimental runtime override for testing GauXC versions that
  provide distributed OneDFT nuclear gradients. The default remains the conservative replicated
  single-rank gradient runtime for MPI calculations.
- `CP2K_GAUXC_STATUS_STDERR=1` mirrors GauXC status messages to standard error. This is useful when
  launcher or CI logs hide the CP2K output file after an external-library failure.
- Some OpenBLAS/libtorch combinations can be sensitive to BLAS symbol resolution for TorchScript
  models using batched matrix products. If a SKALA run crashes in `cblas_sgemm_batch`, use a
  compatible BLAS setup or ensure `libtorch_cpu.so` is loaded before `libopenblas.so`.
- `METHOD GAPW` with OneDFT/SKALA is a molecular GauXC matrix path. GauXC evaluates the full XC term
  directly on its molecular quadrature from the AO density. For pseudopotential inputs this is the
  smooth valence density, so CP2K's local/semi-local GAPW one-center XC correction is not used for
  OneDFT/SKALA. GAPW pseudopotential inputs currently support energies only in this path; nuclear
  gradients and molecular virials require a dedicated derivative of the molecular AO/valence-density
  XC path. NLCC pseudopotentials remain unsupported because the frozen core density would need a
  SKALA-consistent feature definition.
- `METHOD GAPW_XC` with GauXC remains disabled pending a dedicated design for the smooth-density and
  one-center XC terms. It must not be used for non-local OneDFT/SKALA models.
- A true compact-cell periodic GauXC path needs a new GauXC interface rather than only a CP2K input
  change. The missing pieces are an explicit cell/lattice descriptor, a unit-cell quadrature or
  CP2K-grid task interface, periodic AO images or an equivalent periodic density representation,
  return of the periodic XC matrix to CP2K's image/k-point layout, and separate force/stress
  derivatives. PBE-through-GauXC on compact Gamma-only `METHOD GPW` test cells should be the first
  quantitative validation target for that interface.
- A CP2K-native SKALA grid path is a separate research direction from periodic GauXC. The SKALA
  TorchScript protocol requires the meta-GGA grid features `density`, `grad`, and `kin`, the
  integration features `grid_coords` and `grid_weights`, and the per-atom packing metadata
  `atomic_grid_weights`, `atomic_grid_sizes`, `atomic_grid_size_bound_shape`, and
  `coarse_0_atomic_coords`. For Gamma-only periodic `METHOD GPW` this maps naturally to CP2K's
  regular real-space grid for the density, gradient, kinetic-energy density, grid coordinates, and
  weights. The non-local SKALA part still requires an explicit, validated atom-to-grid partition
  before it can be considered a compact periodic implementation. MPI runs gather the local GPW grid
  features into the same atom-partitioned feature block on every rank and evaluate the Torch model
  redundantly; this is a correctness path for distributed CP2K grids, not yet a scalable distributed
  SKALA integration algorithm.
- For such a native SKALA path the first useful implementation target is energy, VXC, and nuclear
  gradients on Gamma-only `METHOD GPW` with GTH pseudopotentials, `FUNCTIONAL PBE`, one `GAUXC`
  functional, and a single image. The force path combines CP2K's regular GPW density-response
  integration with the explicit Torch autograd derivative of the SKALA coarse atomic coordinates.
  The current atom partition is a fixed nearest-atom partition between assignment changes; periodic
  stress, k-points, ROKS, ADMM, and NLCC pseudopotentials should remain separate validation steps.
  `METHOD GAPW` is not enabled in this native grid path yet: a valid GAPW extension has to
  reconstruct not only an all-electron density on the CP2K grid, but also the corresponding density
  gradients and kinetic-energy density entering the SKALA features. GAPW pseudopotential
  augmentation corrections remain a separate methodological problem.
- The first native-grid GAPW target should be all-electron GAPW only. The required design is to
  build a SKALA feature block from the reconstructed all-electron density, gradient, and
  kinetic-energy density on a well-defined grid and validate it against the molecular all-electron
  GauXC matrix path before adding pseudopotential augmentation corrections, `GAPW_XC`, periodic
  stresses, or k-points.
- `NATIVE_GRID_DIAGNOSTICS T` prints the electron count, spin moment, and summed grid weights of the
  CP2K-native SKALA feature block passed to Torch. This is intended for compact RKS/UKS,
  wrapped-cell atom-partition, and MPI/OpenMP validation runs.
- `NATIVE_GRID_USE_CUDA T` evaluates the experimental CP2K-native SKALA Torch model on CUDA if the
  linked libtorch provides CUDA support. This is an explicit opt-in because CUDA-exported SKALA
  models can create internal ragged-index tensors on CUDA.
- `NATIVE_GRID_CUDA_DEVICE` controls the visible CUDA device used by the native-grid Torch path. The
  default value `0` preserves the previous behavior of using logical `cuda:0`. A negative value maps
  the MPI-local rank to a visible CUDA device, so multi-rank jobs can use multiple logical CUDA
  devices. The Torch wrapper loads TorchScript models through CPU and then moves them to the
  selected CUDA device, and guards Torch calls with that device so rank-local model internals and
  input tensors stay on the same GPU. Some CUDA-exported TorchScript models still embed `cuda:0`
  graph constants; those exports need per-rank `CUDA_VISIBLE_DEVICES` set before CP2K starts. When
  CUDA is enabled the native-grid diagnostic output also prints the selected logical and visible
  device for each real-space-grid rank.
- Molecular CDFT and mixed CDFT-CI energy calculations can be used with the GauXC matrix path. SKALA
  CDFT coverage is currently limited to smoke tests of the energy and constraint-potential path.
- Response/kernel properties requiring higher XC derivatives are not supported by the GauXC path and
  abort explicitly.
- Supported OneDFT/SKALA force checks use `GRID SUPERFINE` and `PRUNING_SCHEME UNPRUNED` by default.
  Coarser explicit GauXC grids are allowed, but should be treated as accuracy settings.
- `MOLECULAR_VIRIAL` is a finite-system force diagnostic from GauXC nuclear gradients, not a
  periodic stress tensor.
- SKALA regression tests are technical smoke and force-consistency checks. They do not constitute
  scientific validation of the SKALA model.

## PEXSI (low scaling SCF method)

The Pole EXpansion and Selected Inversion (PEXSI) method requires the PEXSI library and two
dependencies (ParMETIS or PT-Scotch and SuperLU_DIST).

- PEXSI is only available via a Spack build of CP2K.
- Pass `-DCP2K_USE_PEXSI=ON` to CMake.

Below are some additional hints that may help in the compilation process:

- For building PT-Scotch, the flag `-DSCOTCH_METIS_PREFIX` in `Makefile.inc` must not be set and the
  flag `-DSCOTCH_PTHREAD` must be removed.
- For building SuperLU_DIST with PT-Scotch, you must set the following in `make.inc`:

```shell
METISLIB = -lscotchmetis -lscotch -lscotcherr
PARMETISLIB = -lptscotchparmetis -lptscotch -lptscotcherr
```

## PLUMED (enables various enhanced sampling methods)

CP2K can be compiled with PLUMED 2.x by passing `-DCP2K_USE_PLUMED=ON` to CMake.

See <https://cp2k.org/howto:install_with_plumed> for full instructions.

## spglib (crystal symmetries tools)

A library for finding and handling crystal symmetries

- The spglib can be downloaded from <https://github.com/atztogo/spglib>
- For building CP2K with the spglib pass `-DCP2K_USE_SPGLIB=ON` to CMake.

## SIRIUS (plane wave calculations)

SIRIUS is a domain specific library for electronic structure calculations with plane wave method.

- The code is available at <https://github.com/electronic-structure/SIRIUS>
- For building CP2K with SIRIUS pass `-DCP2K_USE_SIRIUS=ON` to CMake.
- Pass `-DCP2K_USE_LIBVDWXC=ON` if support is activated in SIRIUS.
- Pass `-DCP2K_USE_SIRIUS_DFTD3=ON` when sirius is compiled with dftd3 support.
- Pass `-DCP2K_USE_SIRIUS_DFTD4=ON` when sirius is compiled with dftd4 support.
- Pass `-DCP2K_USE_SIRIUS_NLCG=ON` when sirius is compiled with nlcg support.
- Pass `-DCP2K_USE_SIRIUS_VCSQNM=ON` when sirius is compiled with variable cell relaxation support.
- See <https://electronic-structure.github.io/SIRIUS-doc/> for more information.

## COSMA (Distributed Communication-Optimal Matrix-Matrix Multiplication Algorithm)

- COSMA is an alternative for the pdgemm routine included in ScaLAPACK. The library supports both
  CPU and GPUs.
- Pass `-DCP2K_USE_COSMA=ON` to CMake to enable support for COSMA.
- See <https://github.com/eth-cscs/COSMA> for more information.

## LibVori (Voronoi Integration for Electrostatic Properties from Electron Density)

- LibVori is a library which enables the calculation of electrostatic properties (charge, dipole
  vector, quadrupole tensor, etc.) via integration of the total electron density in the Voronoi cell
  of each atom.
- Pass `-DCP2K_USE_VORI=ON` to CMake to enable support for LibVori.
- See <https://brehm-research.de/libvori> for more information.
- LibVori also enables support for the BQB file format for compressed trajectories, please see
  <https://brehm-research.de/bqb> for more information as well as the `bqbtool` to inspect BQB
  files.

## Torch (Machine Learning Framework needed for NequIP)

- The C++ API of PyTorch can be downloaded from https://pytorch.org/get-started/locally/.
- Pass `-DCP2K_USE_LIBTORCH=ON` to CMake to enable support for libtorch.

```{caution}
Note that currently pre-built libtorch bundle (up to 2.12.1) is not compatible with CP2K's external
oneMKL linking stack. If you build CP2K with MKL and want to enable libtorch, you may need to build
it by yourself.
```

## SPLA (Matrix-matrix multiplication offloading on GPU)

The SPLA library is a hard dependency of SIRIUS but can also be used as a standalone library. It
provides a generic interface to the blas gemm family with offloading on GPU. Offloading supports
both CUDA and ROCm (HIP), making the functionality available on both NVIDIA and AMD GPUs.

To make the functionality available, pass `-DCP2K_USE_SPLA_GEMM_OFFLOADING=ON` to CMake and compile
SPLA with Fortran interface and GPU support. Please note that only the functions replacing the dgemm
calls with `offload_dgemm` will eventually be offloaded to the GPU. The SPLA library has internal
criteria to decide if it is worth to do the operation on GPU or not. Calls to `offload_dgemm` also
accept pointers on GPU or a combination of them.

## DeePMD-kit (wider range of interaction potentials)

DeePMD-kit - Deep Potential Molecular Dynamics. Support for DeePMD-kit can be enabled by passing
`-DCP2K_USE_DEEPMD=ON` to CMake.

- DeePMD-kit C interface can be downloaded from
  <https://docs.deepmodeling.com/projects/deepmd/en/master/install/install-from-c-library.html>
- For more information see <https://github.com/deepmodeling/deepmd-kit.git>.

## ACE (atomic cluster expansion ML potentials)

Atomic cluster expansion for accurate and transferable interatomic potentials support can be enabled
by passing `-DCP2K_USE_ACE=ON` to CMake.

- the library files can be downloaded from <https://github.com/ICAMS/lammps-user-pace>
- use cmake/make to compile. There is no install, just ensure that the cp2k build process links in
  all three libraries (libpace, libyaml-cpp-pace and libcnpy). Access to ML-PACE/ace
  ML-PACE/ace-evaluator and yaml-cpp/include from the library is also needed (see toolchain for
  example).

## DFTD4 (dispersion correction)

- dftd4 - Generally Applicable Atomic-Charge Dependent London Dispersion Correction.
- Please always use the CMake-built dftd4 package rather than the Meson-built one for CP2K.
- For more information see <https://github.com/dftd4/dftd4>
- Pass `-DCP2K_USE_DFTD4=ON` to CMake.

## libsmeagol (electron transport calculation with current-induced forces)

- libsmeagol is an external library to compute electron transport properties using Non-Equilibrium
  Green Functions (NEGF) method. The library can be downloaded from
  <https://github.com/StefanoSanvitoGroup/libsmeagol>.
- libsmeagol depends on an MPI library and can only be linked with MPI parallel CP2K binaries.
- During the installation, the directories `$(LIBSMEAGOL_DIR)/lib` and `$(LIBGRPP_DIR)/obj` are
  created.
- Pass `-DCP2K_USE_LIBSMEAGOL=ON` to CMake.

## TREXIO (unified computational chemistry format)

TREXIO - Open-source file format and library. Support for TREXIO can be enabled by passing
`-DCP2K_USE_TREXIO=ON` to CMake.

- TREXIO library can be downloaded from <https://github.com/trex-coe/trexio>
- For more information see <https://trex-coe.github.io/trexio/index.html>.

## LibFCI (full-CI active-space solver)

LibFCI is an external library providing a full-CI solver for CP2K active-space calculations. Support
for LibFCI can be enabled by passing `-DCP2K_USE_LIBFCI=ON` to CMake.

- LibFCI can be downloaded from <https://github.com/DCM-Uni-Paderborn/libfci>

## GREENX (basically functionality for GreenX methods (RPA, GW, Laplace-MP2 etc.)

greenX - Open-source file format and library. Support for greenX can be enabled by passing
`-DCP2K_USE_GREENX=ON` to CMake.

- GREENX library can be downloaded from <https://github.com/nomad-coe/greenX>
- For more information see <https://nomad-coe.github.io/greenX/>.

## TBLITE (semiempirical method)

- tblite - Light-weight tight-binding framework
- With tblite you can calculate using GFN2-xTB method.
- Please always use the CMake-built tblite package rather than the Meson-built one for CP2K.
- For more information see <https://github.com/tblite/tblite>
- Pass `-DCP2K_USE_TBLITE=ON` to CMake.

## openPMD (structured output)

openPMD - Open-source data standard and library. Support for openPMD can be enabled in CMake via
`-DCP2K_USE_OPENPMD=ON`. CMake is the only supported way of enabling openPMD, use of `-D__OPENPMD`
as part of DFLAGS may or may not work.

- openPMD-api may be downloaded from <https://github.com/openPMD/openPMD-api/>, a equal to or
  greater than 0.16.1 is required.
- For more information see <https://openpmd-api.readthedocs.io>.
- The version of openPMD-api, determined by OPENPMDAPI_VERSION_GE, must be 0.16.1 or greater.
- openPMD-api must be built against MPI, determined by openPMD_HAVE_MPI.

## HDF5

- Pass `-DCP2K_USE_HDF5=ON` to CMake to enable HDF5 support.
- HDF5 is a hard dependency for SIRIUS and TREXIO, but can also be used by itself to allow
  read/write functionalities of QCSchema files in the active space module.

## MIMIC (multiscale simulations)

MiMiC - Multiscale simulation framework

- Interface realized through MCL library, which can be downloaded from
  <https://https://mimic-project.org>
- For more information about the framework and supported programs see <https://mimic-project.org>
- Pass `-DCP2K_USE_MIMIC=ON` to CMake

## libGint

libGint - A library for the calculation of the Hartree Fock exchange on GPUs

- Compared to a regular XF calculation, the changes needed in the input file are :
- &FORCE_EVAL &DFT &XC &HF HFX_LIBRARY libGint
- &FORCE_EVAL &DFT &XC &HF &MEMORY MAX_MEMORY X
- pass -DCP2K_USE_LIBGINT=ON\` to CMake.
