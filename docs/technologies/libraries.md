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
and the reference BLAS/LAPACK packages. If compiling with MKL, users must pass
`-DCP2K_BLAS_VENDOR=MKL -DCP2K_SCALAPACK_VENDOR=MKL` to CMake to ensure the code is thread-safe. MKL
with multiple OpenMP threads in CP2K requires that CP2K was compiled with the Intel compiler.

On the Mac, BLAS and LAPACK may be provided by Apple's Accelerate framework. If using this
framework, `-DCP2K_BLAS_VENDOR=Apple` must be passed to CMake to account for some interface
incompatibilities between Accelerate and reference BLAS/LAPACK.

## MPI and ScaLAPACK (required for MPI parallel builds)

MPI (version 3) and SCALAPACK are needed for parallel code. (Use the latest versions available and
download all patches!).

:warning: Note that your MPI installation must match the used Fortran compiler. If your computing
platform does not provide MPI, there are several freely available alternatives:

- MPICH2 MPI: <http://www-unix.mcs.anl.gov/mpi/mpich/> (may require `-fallow-argument-mismatch` when
  building with GCC 10)
- OpenMPI MPI: <http://www.open-mpi.org/>
- ScaLAPACK:
  - <http://www.netlib.org/scalapack/>
  - <http://www.netlib.org/lapack-dev/>
  - ScaLAPACK can be part of ACML or cluster MKL. These libraries are recommended if available.
  - Recently a [ScaLAPACK installer](http://www.netlib.org/scalapack/scalapack_installer.tgz) has
    been added that simplifies the installation.

CP2K assumes that the MPI library implements MPI version 3. Older versions of MPI (e.g., MPI 2.0)
are not supported. CP2K can make use of the mpi_f08 module. If its use is requested, pass
`-DCP2K_USE_MPI_F08=ON` to CMake.

## FFTW (improved performance of FFTs)

FFTW can be used to improve FFT speed on a wide range of architectures. It is strongly recommended
to install and use FFTW3. The current version of CP2K works with FFTW 3.X (pass
`-DCP2K_USE_FFTW3=ON` to CMake). It can be downloaded from <http://www.fftw.org>

FFTW is also provided by MKL. Pass `-DCP2K_USE_FFTW3_WITH_MKL=ON` to CMake.

:warning: Note that FFTW must know the Fortran compiler you will use in order to install properly
(e.g., `export F77=gfortran` before configure if you intend to use gfortran).

Since CP2K is OpenMP parallelized, the FFTW3 threading library libfftw3_threads (or libfftw3_omp) is
required. Pass `-DCP2K_ENABLE_FFTW3_THREADS_SUPPORT=ON` or `-DCP2K_ENABLE_FFTW3_OPENMP_SUPPORT=ON`
respectivly to CMake.

## LIBINT (enables methods including HF exchange)

- Hartree-Fock exchange requires the LIBINT package to be installed.
- Pass `-DCP2K_USE_LIBINT2=ON` to CMake to enable LIBINT.
- Recommended way to build LIBINT: Download a CP2K-configured LIBINT library from
  [libint-cp2k](https://github.com/cp2k/libint-cp2k). Build and install LIBINT by following the
  instructions provided there. Note that using a library configured for higher maximum angular
  momentum will increase build time and binary size of CP2K binary (assuming static linking).
- CP2K is not hardwired to these provided libraries and any other LIBINT library (version >= 2.5.0)
  should be compatible as long as it was compiled with `--enable-eri=1` and default ordering.
- Avoid debugging information (`-g` flag) for compiling LIBINT since this will increase library size
  by a large factor.

## LIBXSMM (improved performance for matrix multiplication)

- A library for matrix operations and deep learning primitives:
  <https://github.com/libxsmm/libxsmm/>.
- Pass `-DCP2K_USE_LIBXSMM=ON` to CMake to enable it.
- LIBXSMM can be used with both CUDA and HIP backends (see [](./accelerators/index.md)).

## LIBXC (wider choice of xc functionals)

- The version 5.1.0 (or later) of LIBXC can be downloaded from
  <https://www.tddft.org/programs/libxc>
- CP2K does not make use of fourth derivates such that LIBXC may be configured with
  `./configure --disable-lxc <other LIBXC configuration flags>`.
- During the installation, the directories `$(LIBXC_DIR)/lib` and `$(LIBXC_DIR)/include` are
  created.
- Pass `-DCP2K_USE_LIBXC=ON` to CMake.

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

SIRIUS is a domain specific library for electronic structure calculations.

- The code is available at <https://github.com/electronic-structure/SIRIUS>
- For building CP2K with SIRIUS pass `-DCP2K_USE_SIRIUS=ON` to CMake.
- Pass `-DCP2K_USE_LIBVDWXC=ON` if support is activated in SIRIUS.
- Pass `-DCP2K_USE_SIRIUS_DFTD4=ON` when sirius is compiled with dftd3 and dftd4 support.
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

## SpFFT (Sparse 3D FFT)

- SpFFT is a 3D FFT library for sparse frequency domain data written in C++ with support for MPI,
  OpenMP, CUDA and ROCm.
- Pass `-DCP2K_USE_SpFFT=ON` to CMake to enable support for SpFFT.
- See <https://github.com/eth-cscs/SpFFT> for more information.

## Torch (Machine Learning Framework needed for NequIP)

- The C++ API of PyTorch can be downloaded from https://pytorch.org/get-started/locally/.
- Pass `-DCP2K_USE_LIBTORCH=ON` to CMake to enable support for libtorch.

## matrix-matrix multiplication offloading on GPU using SPLA

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

## GREENX (basically functionality for GreenX methods (RPA, GW, Laplace-MP2 etc.)

greenX - Open-source file format and library. Support for greenX can be enabled by passing
`-DCP2K_USE_GREENX=ON` to CMake.

- GREENX library can be downloaded from <https://github.com/nomad-coe/greenX>
- For more information see <https://nomad-coe.github.io/greenX/>.

## TBLITE (semiempirical method)

- tblite - Light-weight tight-binding framework
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
