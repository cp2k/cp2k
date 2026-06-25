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
Examples are the sequential or thread variant of the Intel MKL, the Cray libsci, the OpenBLAS OpenMP
variant and the reference BLAS/LAPACK packages. Usually the CMake step of CP2K will auto-detect the
type of BLAS and SCALAPACK and then use the right configuration to ensure the code is thread-safe;
however, if detection is ambiguous, `-DCP2K_BLAS_VENDOR=MKL` and `-DCP2K_SCALAPACK_VENDOR=MKL` can
be used for a oneMKL installation

On the Mac, BLAS and LAPACK can be provided by either OpenBLAS or Apple's Accelerate framework.

## DBCSR (required, block-sparse matrix operations)

CP2K requires [DBCSR](https://github.com/cp2k/dbcsr/) for block-sparse matrix operations. It is
found automatically by CMake. For MPI builds, DBCSR must also have been built with MPI support. The
CP2K toolchain, Spack, and a compatible external DBCSR installation are supported ways to provide
it.

## MPI and ScaLAPACK (required for MPI parallel builds)

MPI (version 3 or later) and SCALAPACK are needed for parallel code. (Use the latest versions
available and download all patches!).

```{warning}
Note that the MPI installation must match the used Fortran compiler.
```

If your computing platform does not provide MPI, there are several freely available implementations:

- MPICH: <https://www.mpich.org/> (may require `-fallow-argument-mismatch` when building with GCC
  10\)
- OpenMPI: <http://www.open-mpi.org/>

For more information of ScaLAPACK, see <http://www.netlib.org/scalapack/>. ScaLAPACK can be part of
AOCL (AMD) or oneMKL (Intel); these libraries are recommended on the corresponding machines if
available.

CP2K assumes that the MPI library implements MPI version 3. Older versions of MPI (e.g., MPI 2.0)
are not supported. CP2K can make use of the `mpi_f08` module; pass `-DCP2K_USE_MPI_F08=ON` to CMake
to enable it.

## FFTW (improved performance of FFTs)

FFTW can be used to improve FFT speed on a wide range of architectures. It is strongly recommended
to install and use FFTW3. The current version of CP2K works with FFTW 3.X (pass
`-DCP2K_USE_FFTW3=ON` to CMake). It can be downloaded from <http://www.fftw.org>.

FFTW is also provided by MKL. If you have MKL but still want to use standalone FFTW3, pass
`-DCP2K_USE_FFTW3_WITH_MKL=ON` to CMake.

```{warning}
Note that FFTW must know the Fortran compiler you will use in order to install properly
(e.g., `export F77=gfortran` before configure if you intend to use gfortran).
```

Since CP2K is OpenMP parallelized, CP2K enables the FFTW3 OpenMP interface by default
(`-DCP2K_ENABLE_FFTW3_OPENMP_SUPPORT=ON`); the FFTW installation must therefore provide
`libfftw3_omp`. The alternative threads interface can be selected with
`-DCP2K_ENABLE_FFTW3_THREADS_SUPPORT=ON`, which requires `libfftw3_threads`.

```{important}
Support for FFTW is required for some features, especially systems with very large block sizes/grid
sizes. A future release of CP2K may make FFTW a hard dependency. Please consider CP2K to be compiled
with support for FFTW.
```

## LIBINT (ERI calculation for HFX)

[Libint2](https://github.com/evaleev/libint) provides the electron-repulsion integrals required for
Hartree--Fock exchange and related methods.

- Pass `-DCP2K_USE_LIBINT2=ON` to CMake to enable Libint2.
- The CP2K toolchain and Spack are the recommended ways to obtain a compatible Libint2 build. For a
  manual build, CP2K-configured Libint source packages are available from the
  [CP2K download server](https://www.cp2k.org/static/downloads/); see also
  [the CP2K Libint instructions](https://github.com/cp2k/cp2k/blob/master/tools/libint/README.md). A
  library configured for a higher maximum angular momentum increases CP2K compilation time and,
  particularly for static builds, the binary size.
- CP2K is not restricted to these source packages. A manually built installation must provide
  electron-repulsion integrals (`--enable-eri=1`) with Libint's default ordering, export a CMake
  package, and include the Fortran interface (`libint_f.mod`).
- Avoid compiling Libint with extensive debug information unless it is specifically required, since
  this can increase the library size substantially.

## LIBXS (improved performance for matrix multiplication)

- A library for matrix operations and deep learning primitives: <https://github.com/hfp/libxs/>.
- [LIBXS](https://github.com/hfp/libxs/) provides optimized matrix operations used by CP2K and is
  required when using CP2K's OpenCL backend.
- Pass `-DCP2K_USE_LIBXS=ON` to CMake to enable it.

## LIBXSTREAM (OpenCL offload runtime)

- [LIBXSTREAM](https://github.com/hfp/libxstream) provides the stream and memory-management layer
  used by CP2K's OpenCL offload backend.
- It is required automatically when configuring OpenCL acceleration with `-DCP2K_USE_ACCEL=OPENCL`;
  there is no separate `CP2K_USE_LIBXSTREAM` option.
- OpenCL builds also require LIBXS. For a manual build, make both LIBXSTREAM and the OpenCL
  development files discoverable by CMake, for example through `CMAKE_PREFIX_PATH`.
- See [](./accelerators/opencl.md) for OpenCL runtime and backend-specific requirements.

## LIBXSMM (JIT-kernel provider of libXS)

- A library that provide a JIT-kernel for LibXS: <https://github.com/libxsmm/libxsmm/>.
- Pass `-DCP2K_USE_LIBXSMM=ON` to CMake to enable it; this is only valid when `-DCP2K_USE_LIBXS=ON`
  is passed.
- The integration of LIBXS and LIBXSMM is done through `libxs_jit.F` that is provided by LIBXS but
  compiled by DBCSR and CP2K.
- LIBXSMM can be used with both CUDA and HIP backends (see [](./accelerators/index.md)).

## LIBXC (wider choice of xc functionals)

LIBXC is a library that provides wider choice of XC functionals.

- The latest version of LIBXC can be downloaded from <https://gitlab.com/libxc/libxc/-/releases>
- CP2K makes use of third derivates but does not use fourth derivates, so LIBXC may be configured
  with `cmake .. -DDISABLE_KXC=OFF <other LIBXC configuration flags>`.
- Pass `-DCP2K_USE_LIBXC=ON` to CMake.
- [LIBXSMM](https://github.com/libxsmm/libxsmm/) provides just-in-time kernels for LIBXS.
- Pass `-DCP2K_USE_LIBXSMM=ON` to CMake to enable it; this option is valid only when
  `-DCP2K_USE_LIBXS=ON` is also passed.
- The integration of LIBXS and LIBXSMM uses `libxs_jit.F`, which is provided by LIBXS and compiled
  by DBCSR and CP2K.
- LIBXSMM can be used with both CUDA and HIP backends; see [](./accelerators/index.md).

## GauXC (xc integration library)

GauXC can be used to evaluate selected exchange-correlation functionals through an external
integrator.

- Libtorch is required for OneDFT/SKALA support.
- Pass `-DCP2K_USE_GAUXC=ON` to CMake to enable GauXC. An MPI-enabled CP2K build requires a GauXC
  installation built with MPI support.
- GauXC support is primarily intended for isolated QS calculations. Coverage of compact periodic
  systems, k-points, periodic stress tensors, and GAPW-related paths is limited; consult the Input
  Reference for the current `&GAUXC` support and input restrictions.
- TorchScript-based GauXC models require a libtorch installation compatible with CP2K's BLAS and
  OpenMP runtime. Pre-built libtorch bundles can conflict with a CP2K build using oneMKL; use a
  compatible generic BLAS stack or rebuild libtorch against the selected dynamic stack when this
  occurs.
- `CP2K_GAUXC_STATUS_STDERR=1` mirrors GauXC status messages to standard error, which can be useful
  when launcher or CI logs do not preserve the CP2K output file after an external-library failure.

## PEXSI (low scaling SCF method)

The Pole EXpansion and Selected Inversion (PEXSI) method requires an MPI build and a compatible
PEXSI CMake package. PEXSI itself depends on a sparse-direct-solver and graph-partitioning stack,
typically SuperLU_DIST together with ParMETIS or PT-Scotch.

- Pass `-DCP2K_USE_PEXSI=ON` to CMake to enable PEXSI.
- Spack is the most convenient supported route for provisioning the complete PEXSI dependency stack.
  Manual builds are also possible when a compatible PEXSI installation and its dependencies are
  available to CMake. It's not supported to install PEXSI through toolchain

## PLUMED (enables various enhanced sampling methods)

CP2K can be compiled with PLUMED 2.x by passing `-DCP2K_USE_PLUMED=ON` to CMake.

See <https://cp2k.org/howto:install_with_plumed> for full instructions.

## spglib (crystal symmetries tools)

Spglib is a library for finding and handling crystal symmetries.

- The library can be downloaded from <https://github.com/atztogo/spglib>
- For building CP2K with the spglib pass `-DCP2K_USE_SPGLIB=ON` to CMake.

## SIRIUS (plane wave calculations)

SIRIUS is a domain specific library for electronic structure calculations with plane wave method.

- The code is available at <https://github.com/electronic-structure/SIRIUS>.
- SIRIUS support requires an MPI build. Pass `-DCP2K_USE_SIRIUS=ON` to CMake to enable it.
- SIRIUS has its own dependency stack, commonly including HDF5, SpFFT, SPLA, and eigensolver
  libraries. It's recommended to build SIRIUS through Spack to get all features enabled.
- Pass `-DCP2K_USE_LIBVDWXC=ON` when the selected SIRIUS build provides libvdwxc support.
- Pass `-DCP2K_USE_SIRIUS_DFTD3=ON` when SIRIUS was built with DFT-D3 support.
- Pass `-DCP2K_USE_SIRIUS_DFTD4=ON` when SIRIUS was built with DFT-D4 support.
- Pass `-DCP2K_USE_SIRIUS_NLCG=ON` when SIRIUS was built with NLCG support.
- Pass `-DCP2K_USE_SIRIUS_VCSQNM=ON` when SIRIUS was built with variable-cell-relaxation support.
- See <https://electronic-structure.github.io/SIRIUS-doc/> for build options and supported features.

## COSMA (Distributed Communication-Optimal Matrix-Matrix Multiplication Algorithm)

COSMA is an alternative for the pdgemm routine included in ScaLAPACK. The library supports both CPU
and GPUs.

- Pass `-DCP2K_USE_COSMA=ON` to CMake to enable support for COSMA.
- See <https://github.com/eth-cscs/COSMA> for more information.

## LibVori (Voronoi Integration for Electrostatic Properties from Electron Density)

LibVori is a library which enables the calculation of electrostatic properties (charge, dipole
vector, quadrupole tensor, etc.) via integration of the total electron density in the Voronoi cell
of each atom.

- Pass `-DCP2K_USE_VORI=ON` to CMake to enable support for LibVori.
- See <https://brehm-research.de/libvori> for more information.
- LibVori also enables support for the BQB file format for compressed trajectories, please see
  <https://brehm-research.de/bqb> for more information as well as the `bqbtool` to inspect BQB
  files.

## Torch (Machine Learning Framework needed for NequIP)

LibTorch is used by TorchScript-based functionality, including GauXC OneDFT/SKALA support.

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

SPLA support requires an MPI build and is enabled with `-DCP2K_USE_SPLA=ON`. To offload eligible
`dgemm` operations, additionally pass `-DCP2K_USE_SPLA_GEMM_OFFLOADING=ON` and enable CUDA or HIP.
SPLA must be built with its Fortran interface and a GPU backend. SPLA decides at runtime whether an
individual operation is suitable for offloading.

## DeePMD-kit (wider range of interaction potentials)

DeePMD-kit provides Deep Potential models. Support for its C interface can be enabled by passing
`-DCP2K_USE_DEEPMD=ON` to CMake.

- DeePMD-kit C interface can be downloaded from
  <https://docs.deepmodeling.com/projects/deepmd/en/master/install/install-from-c-library.html>
- For more information see <https://github.com/deepmodeling/deepmd-kit.git>.

## ACE (atomic cluster expansion ML potentials)

Atomic cluster expansion potentials from ML-PACE for accurate and transferable interatomic
potentials support can be enabled by passing `-DCP2K_USE_ACE=ON` to CMake.

- the library files can be downloaded from <https://github.com/ICAMS/lammps-user-pace>
- use cmake/make to compile. There is no install, just ensure that the cp2k build process links in
  all three libraries (libpace, libyaml-cpp-pace and libcnpy). Access to ML-PACE/ace
  ML-PACE/ace-evaluator and yaml-cpp/include from the library is also needed (see toolchain for
  example).

## DFTD4 (dispersion correction)

DFTD4 provides the Generally Applicable Atomic-Charge Dependent London Dispersion Correction.

- Please use the CMake-built dftd4 package rather than the Meson-built one for CP2K.
- For more information, see <https://github.com/dftd4/dftd4>
- Pass `-DCP2K_USE_DFTD4=ON` to CMake.

## libsmeagol (electron transport calculation with current-induced forces)

libsmeagol is an external library to compute electron transport properties using Non-Equilibrium
Green Functions (NEGF) method. The library can be downloaded from
<https://github.com/StefanoSanvitoGroup/libsmeagol>.

- libsmeagol depends on an MPI library and can only be linked with MPI parallel CP2K binaries.
- During the installation, the directories `$(LIBSMEAGOL_DIR)/lib` and `$(LIBGRPP_DIR)/obj` are
  created.
- Pass `-DCP2K_USE_LIBSMEAGOL=ON` to CMake.

## TREXIO (unified computational chemistry format)

TREXIO is an open-source file format and library. Support can be enabled by passing
`-DCP2K_USE_TREXIO=ON` to CMake. HDF5 is required.

- TREXIO library can be downloaded from <https://github.com/trex-coe/trexio>
- For more information see <https://trex-coe.github.io/trexio/index.html>.

## LibFCI (full-CI active-space solver)

LibFCI is an external library providing a full-CI solver for CP2K active-space calculations. Support
for LibFCI can be enabled by passing `-DCP2K_USE_LIBFCI=ON` to CMake.

- LibFCI can be downloaded from <https://github.com/DCM-Uni-Paderborn/libfci>

## GREENX (GreenX methods such as RPA, GW, and Laplace-MP2)

GreenX provides functionality for GreenX methods such as RPA, GW, and Laplace-MP2. Support can be
enabled by passing `-DCP2K_USE_GREENX=ON` to CMake.

- GREENX library can be downloaded from <https://github.com/nomad-coe/greenX>
- For more information see <https://nomad-coe.github.io/greenX/>.

## TBLITE (semiempirical method)

TBLITE is a lightweight tight-binding framework that provides the GFN2-xTB method.

- Please always use the CMake-built tblite package rather than the Meson-built one for CP2K.
- A CMake build of tblite from source also installs DFT-D4 and s-dftd3. Therefore, no separate DFTD4
  installation is needed when tblite is enabled; s-dftd3 also provides parameters for additional XC
  functionals.
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

- Its interface is realized through the MCL library, which can be downloaded from
  <https://https://mimic-project.org>
- For more information about the framework and supported programs see <https://mimic-project.org>
- Pass `-DCP2K_USE_MIMIC=ON` to CMake
