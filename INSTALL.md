# How to compile the CP2K code

## 1. Acquire the code

For users, the preferred method is to [download a release](https://github.com/cp2k/cp2k/releases/)
(use the versioned tarballs, `cp2k-X.Y.tar.bz2`). For developers, the preferred method is to
[download from Git](./README.md#downloading-cp2k-source-code).

For more details on downloading CP2K, see <https://www.cp2k.org/download>.

## 2. Install prerequisites

The easiest way to build CP2K with all its dependencies is as a
[Docker container](./tools/docker/README.md).

Alternatively, the [toolchain script](./tools/toolchain/install_cp2k_toolchain.sh) can also be run
directly.

For a complete introduction to the toolchain script, see the [README](./tools/toolchain/README.md).

The basic steps are:

- Read toolchain installation options:

```shell
cd tools/toolchain/
./install_cp2k_toolchain.sh --help
```

- Launch toolchain script (example option choice)

```shell
./install_cp2k_toolchain.sh --with-libxsmm=install --with-openblas=system \
     --with-fftw=system --with-reflapack=no  --enable-cuda
```

- Once the script has completed successfully, follow the instructions given at the end of its
  output. Note that the pre-built arch files provided by the toolchain are for the GNU compiler,
  users must adapt them for other compilers. It is possible to use the provided [arch files](./arch)
  as guidance.

There are [arch files](./arch) for a few specific platforms (e.g.
[Linux-gnu-x86_64](./arch/Linux-gnu-x86_64.psmp),
[Linux-intel-x86_64](./arch/Linux-intel-x86_64.psmp)) which include a toolchain build. Sourcing such
an arch file in the cp2k folder launches a toolchain build, e.g.

```
source ./arch/Linux-gnu-x86_64.psmp
```

After a successful toolchain build, run one of the suggested `make` commands

```
make -j ARCH=Linux-gnu-x86_64 VERSION=psmp
```

Check also the corresponding [HowTos](https://www.cp2k.org/howto/) for
[Apple M1 (macOS)](https://www.cp2k.org/howto:compile_on_macos/) and
[Cray XC40/50 (Piz Daint, CSCS)](https://www.cp2k.org/howto:compile_on_cray_cscs/).

Sub-points here discuss prerequisites needed to build CP2K. Copies of the recommended versions of
3rd party software can be downloaded from <https://www.cp2k.org/static/downloads/>.

Generally, CP2K supports only one version for each of its dependencies. These are defined by the
[toolchain scripts](./tools/toolchain/scripts/). Other versions might work too, but we don't test
them. So, your mileage may vary.

### 2a. GNU make (required, build system)

GNU make should be on your system (gmake or make on linux) and used for the build, go to
<https://www.gnu.org/software/make/make.html> download from <https://ftp.gnu.org/pub/gnu/make/>.

### 2b. Python (required, build system)

Python 3.5+ is needed to run the dependency generator. On most system Python is already installed.
For more information visit: <https://www.python.org>

### 2c. Fortran and C Compiler (required, build system)

A Fortran 2008 compiler and matching C99 compiler should be installed on your system. We have good
experience with gcc/gfortran (gcc >=4.6 works, later version recommended). Be aware that some
compilers have bugs that might cause them to fail (internal compiler errors, segfaults) or, worse,
yield a mis-compiled CP2K. Report bugs to compiler vendors; they (and we) have an interest in fixing
them. A list of tested compiler can be found [here](https://www.cp2k.org/dev:compiler_support).
Always run a `make -j test` (See point 5.) after compilation to identify these problems.

### 2d. BLAS and LAPACK (required, base functionality)

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
and the reference BLAS/LAPACK packages. If compiling with MKL, users must define `-D__MKL` to ensure
the code is thread-safe. MKL with multiple OpenMP threads in CP2K requires that CP2K was compiled
with the Intel compiler. If the `cpp` precompiler is used in a separate precompilation step in
combination with the Intel Fortran compiler, `-D__INTEL_LLVM_COMPILER` (`-D__INTEL_COMPILER`) must
be added explicitly (the Intel compiler sets `D__INTEL_LLVM_COMPILER` otherwise automatically).

On the Mac, BLAS and LAPACK may be provided by Apple's Accelerate framework. If using this
framework, `-D__ACCELERATE` must be defined to account for some interface incompatibilities between
Accelerate and reference BLAS/LAPACK.

When building on/for Windows using the Minimalist GNU for Windows (MinGW) environment, you must set
`-D__MINGW`, `-D__NO_STATM_ACCESS` and `-D__NO_SOCKETS` to avoid undefined references during
linking, respectively errors while printing the statistics.

### 2e. MPI and ScaLAPACK (optional, required for MPI parallel builds)

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
are not supported and the old flag `-D__MPI_VERSION` in the arch file will be ignored. CP2K can make
use of the mpi_f08 module. If its use is requested, set the flag `-D__MPI_F08`.

### 2f. FFTW (optional, improved performance of FFTs)

FFTW can be used to improve FFT speed on a wide range of architectures. It is strongly recommended
to install and use FFTW3. The current version of CP2K works with FFTW 3.X (use `-D__FFTW3`). It can
be downloaded from <http://www.fftw.org>

FFTW is also provided by MKL. Use `-D__FFTW3_MKL` to use the correct import path.

:warning: Note that FFTW must know the Fortran compiler you will use in order to install properly
(e.g., `export F77=gfortran` before configure if you intend to use gfortran).

:warning: Note that on machines and compilers which support SSE you can configure FFTW3 with
`--enable-sse2`. Compilers/systems that do not align memory (NAG f95, Intel IA32/gfortran) should
either not use `--enable-sse2` or otherwise set the define `-D__FFTW3_UNALIGNED` in the arch file.
Since CP2K is OpenMP parallelized, the FFTW3 threading library libfftw3_threads (or libfftw3_omp) is
required.

### 2g. LIBINT (optional, enables methods including HF exchange)

- Hartree-Fock exchange (optional, use `-D__LIBINT`) requires the LIBINT package to be installed.
- Recommended way to build LIBINT: Download a CP2K-configured LIBINT library from
  [libint-cp2k](https://github.com/cp2k/libint-cp2k). Build and install LIBINT by following the
  instructions provided there. Note that using a library configured for higher maximum angular
  momentum will increase build time and binary size of CP2K binary (assuming static linking).
- CP2K is not hardwired to these provided libraries and any other LIBINT library (version >= 2.5.0)
  should be compatible as long as it was compiled with `--enable-eri=1` and default ordering.
- Avoid debugging information (`-g` flag) for compiling LIBINT since this will increase library size
  by a large factor.
- In the arch file of CP2K: add `-D__LIBINT` to the `DFLAGS`. Add
  `-L$(LIBINT_DIR)/lib -lint2 -lstdc++` to `LIBS` and `-I$(LIBINT_DIR)/include` to `FCFLAGS`.
  `lstdc++` is needed if you use the GNU C++ compiler.
- Libint 1 is no longer supported and the previously needed flags `-D__LIBINT_MAX_AM` and
  `-D__LIBDERIV_MAX_AM1` are ignored.
- `-D__MAX_CONTR=4` (default=2) can be used to compile efficient contraction kernels up to l=4, but
  the build time will increase accordingly.

### 2h. LIBXSMM (optional, improved performance for matrix multiplication)

- A library for matrix operations and deep learning primitives:
  <https://github.com/libxsmm/libxsmm/>.
- Add `-D__LIBXSMM` to enable it, with suitable include and library paths, e.g.,
  `FCFLAGS += -I${LIBXSMM_DIR}/include -D__LIBXSMM` and
  `LIBS += -L${LIBXSMM_DIR}/lib -lxsmmf -lxsmmext -lxsmm -ldl`
- LIBSMM is not used if LIBXSMM is enabled.

### 2i. CUDA (optional, improved performance on GPU systems)

- Specify OFFLOAD_CC (e.g., `OFFLOAD_CC = nvcc`) and OFFLOAD_FLAGS (e.g.,
  `OFFLOAD_FLAGS = -O3 -g -w --std=c++11`) variables. Remember to include the support for the C++11
  standard.
- Use `-D__OFFLOAD_CUDA` to generally enable support for Nvidia GPUs
- Use the `-D__DBCSR_ACC` and `OFFLOAD_TARGET = cuda` to enable accelerator support for matrix
  multiplications.
- Add `-lstdc++ -lcudart -lnvrtc -lcuda -lcublas` to LIBS.
- Specify the GPU type (e.g., `GPUVER = P100`), possible values are K20X, K40, K80, P100, V100,
  A100, H100, A40.
- Specify the C++ compiler (e.g., `CXX = g++`) and the CXXFLAGS to support the C++11 standard.
- CUFFT 7.0 has a known bug and is therefore disabled by default. NVIDIA's webpage list a patch (an
  upgraded version cufft i.e. >= 7.0.35) - use this together with `-D__HAS_PATCHED_CUFFT_70`.
- Use `-D__OFFLOAD_PROFILING` to turn on Nvidia Tools Extensions. It requires to link
  `-lnvToolsExt`.
- Link to a BLAS/ScaLAPACK library that accelerates large DGEMMs (e.g., libsci_acc)
- Use `-D__NO_OFFLOAD_GRID` to disable the GPU backend of the grid library.
- Use `-D__NO_OFFLOAD_DBM` to disable the GPU backend of the sparse tensor library.
- Use `-D__NO_OFFLOAD_PW` to disable the GPU backend of FFTs and associated gather/scatter
  operations.

### 2j. LIBXC (optional, wider choice of xc functionals)

- The version 5.1.0 (or later) of LIBXC can be downloaded from
  <https://www.tddft.org/programs/libxc>
- CP2K does not make use of fourth derivates such that LIBXC may be configured with './configure
  --disable-lxc \<other LIBXC configuration flags>'.
- During the installation, the directories `$(LIBXC_DIR)/lib` and `$(LIBXC_DIR)/include` are
  created.
- Add `-D__LIBXC` to DFLAGS, `-I$(LIBXC_DIR)/include` to FCFLAGS and
  `-L$(LIBXC_DIR)/lib -lxcf03 -lxc` to LIBS.
- :warning: Note that the deprecated flags `-D__LIBXC2` and `-D__LIBXC3` are ignored.

### 2k. ELPA (optional, improved performance for diagonalization)

Library ELPA for the solution of the eigenvalue problem

- ELPA replaces the ScaLAPACK `SYEVD` to improve the performance of the diagonalization
- A version of ELPA can be downloaded from <http://elpa.rzg.mpg.de/software>.
- During the installation the `libelpa_openmp.a` is created.
- Minimal supported version of ELPA is 2018.05.001.
- Add `-D__ELPA` to `DFLAGS`
- Add `-D__ELPA_NVIDIA_GPU`, `-D__ELPA_AMD_GPU`, or `-D__ELPA_INTEL_GPU` to `DFLAGS` to enable GPU
  support for the respective vendor.
- Add `-I$(ELPA_INCLUDE_DIR)/modules` to `FCFLAGS`
- Add `-I$(ELPA_INCLUDE_DIR)/elpa` to `FCFLAGS`
- Add `-L$(ELPA_DIR)` to `LDFLAGS`
- Add `-lelpa` to `LIBS`
- For specific architectures it can be better to install specifically optimized kernels (see BG)
  and/or employ a higher optimization level to compile it.

### 2l. cuSOLVERMp (experimental, improved performance for diagonalization on Nvidia GPUs)

NVIDIA cuSOLVERMp is a high-performance, distributed-memory, GPU-accelerated library that provides
tools for the solution of dense linear systems and eigenvalue problems.

- cuSOLVERMp replaces the ScaLAPACK `SYEVD` to improve the performance of the diagonalization
- A version of cuSOLVERMp can be downloaded from <https://docs.nvidia.com/hpc-sdk/cusolvermp>.
- Add `-D__CUSOLVERMP` to `DFLAGS`
- Add `-lcusolverMp -lcusolver -lcal -lnvidia-ml` to `LIBS`

### 2m. DLA-Future (optional, improved performance for diagonalization on Nvidia and AMD GPUs)

[DLA-Future](https://github.com/eth-cscs/DLA-Future) is a high-performance, distributed-memory,
GPU-accelerated library that provides tools for the solution of eigenvalue problems, based on the
[pika](https://pikacpp.org/) runtime.
[DLA-Future-Fortran](https://github.com/eth-cscs/DLA-Future-Fortran) provides a Fortran interface to
DLA-Future.

- DLA-Future-Fortran replaces the ScaLAPACK `SYEVD`, `HEEVD`, and `HEGVD` to improve performance of
  the diagonalization
- DLA-Future is available at <https://github.com/eth-cscs/DLA-Future>
- DLA-Future-Fortran is available at <https://github.com/eth-cscs/DLA-Future-Fortran>
- DLA-Future is available via the [Spack](https://packages.spack.io/package.html?name=dla-future)
  package manager
- DLA-Future-Fortran is available via the
  [Spack](https://packages.spack.io/package.html?name=dla-future-fortran) package manager
- `-D__DLAF` is defined by CMake when `-DCP2K_USE_DLAF=ON`

### 2n. PEXSI (optional, low scaling SCF method)

The Pole EXpansion and Selected Inversion (PEXSI) method requires the PEXSI library and two
dependencies (ParMETIS or PT-Scotch and SuperLU_DIST).

- Download PEXSI (www.pexsi.org) and install it and its dependencies by following its README.md.
- PEXSI versions 0.10.x have been tested with CP2K. Older versions are not supported.
- PEXSI needs to be built with `make finstall`.

In the arch file of CP2K:

- Add `-lpexsi_${SUFFIX} -llapack -lblas -lsuperlu_dist_3.3 -lparmetis -lmetis`, and their paths
  (with `-L$(LIB_DIR)`) to LIBS.
- It is important that a copy of LAPACK and BLAS is placed before and after these libraries (replace
  `-llapack` and `-lblas` with the optimized versions as needed).
- In order to link in PT-Scotch instead of ParMETIS replace `-lparmetis -lmetis` with:
  `-lptscotchparmetis -lptscotch -lptscotcherr -lscotchmetis -lscotch -lscotcherr`
- Add `-I$(PEXSI_DIR)/fortran/` to FCFLAGS.
- Add `-D__LIBPEXSI` to DFLAGS.

Below are some additional hints that may help in the compilation process:

- For building PT-Scotch, the flag `-DSCOTCH_METIS_PREFIX` in `Makefile.inc` must not be set and the
  flag `-DSCOTCH_PTHREAD` must be removed.
- For building SuperLU_DIST with PT-Scotch, you must set the following in `make.inc`:

```shell
METISLIB = -lscotchmetis -lscotch -lscotcherr
PARMETISLIB = -lptscotchparmetis -lptscotch -lptscotcherr
```

### 2o. QUIP (optional, wider range of interaction potentials)

QUIP - QUantum mechanics and Interatomic Potentials Support for QUIP can be enabled via the flag
`-D__QUIP`.

For more information see <http://www.libatoms.org>.

### 2p. PLUMED (optional, enables various enhanced sampling methods)

CP2K can be compiled with PLUMED 2.x (`-D__PLUMED2`).

See <https://cp2k.org/howto:install_with_plumed> for full instructions.

### 2q. spglib (optional, crystal symmetries tools)

A library for finding and handling crystal symmetries

- The spglib can be downloaded from <https://github.com/atztogo/spglib>
- For building CP2K with the spglib add `-D__SPGLIB` to DFLAGS

### 2r. SIRIUS (optional, plane wave calculations)

SIRIUS is a domain specific library for electronic structure calculations.

- The code is available at <https://github.com/electronic-structure/SIRIUS>
- For building CP2K with SIRIUS add `-D__SIRIUS` to DFLAGS.
- Add `-D__LIBVDWXC` if support is activated in SIRIUS.
- See <https://electronic-structure.github.io/SIRIUS-doc/> for more information.

### 2s. FPGA (optional, plane wave FFT calculations)

- Use `-D__PW_FPGA` to enable FPGA support for PW (fft) calculations. Currently tested only for
  Intel Stratix 10 and Arria 10 GX1150 FPGAs.
- Supports single precision and double precision fft calculations with the use of dedicated APIs.
- Double precision is the default API chosen when set using the `-D__PW_FPGA` flag.
- Single precision can be set using an additional `-D__PW_FPGA_SP` flag along with the `-D__PW_FPGA`
  flag.
- Kernel code must be synthesized separately and copied to a specific location.
- See <https://github.com/pc2/fft3d-fpga/> for the kernel code and instructions for synthesis.
- Read `src/pw/fpga/README.md` for information on the specific location to copy the binaries to.
- Currently supported FFT3d sizes - 16^3, 32^3, 64^3.
- Include aocl compile flags and `-D__PW_FPGA -D__PW_FPGA_SP` to `CFLAGS`, aocl linker flags to
  `LDFLAGS` and aocl libs to `LIBS`.
- When building FPGA and OFFLOAD together then `-D__NO_OFFLOAD_PW` must be used.

### 2t. COSMA (Distributed Communication-Optimal Matrix-Matrix Multiplication Algorithm)

- COSMA is an alternative for the pdgemm routine included in ScaLAPACK. The library supports both
  CPU and GPUs.
- Add `-D__COSMA` to the DFLAGS to enable support for COSMA.
- See <https://github.com/eth-cscs/COSMA> for more information.

### 2u. LibVori (Voronoi Integration for Electrostatic Properties from Electron Density)

- LibVori is a library which enables the calculation of electrostatic properties (charge, dipole
  vector, quadrupole tensor, etc.) via integration of the total electron density in the Voronoi cell
  of each atom.
- Add `-D__LIBVORI` to the DFLAGS to enable support for LibVori.
- See <https://brehm-research.de/libvori> for more information.
- LibVori also enables support for the BQB file format for compressed trajectories, please see
  <https://brehm-research.de/bqb> for more information as well as the `bqbtool` to inspect BQB
  files.

### 2v. Torch (Machine Learning Framework needed for NequIP)

- The C++ API of PyTorch can be downloaded from https://pytorch.org/get-started/locally/.
- Add `-D__LIBTORCH` to the DFLAGS to enable support for libtorch.

### 2w. ROCM/HIP (Support for AMD GPU)

The code for the HIP based grid backend was developed and tested on Mi100 but should work out of the
box on Nvidia hardware as well.

- Use `-D__OFFLOAD_HIP` to generally enable support for AMD GPUs
- Use `-D__NO_OFFLOAD_GRID` to disable the GPU backend of the grid library.
- Use `-D__NO_OFFLOAD_DBM` to disable the GPU backend of the sparse tensor library.
- Use `-D__NO_OFFLOAD_PW` to disable the GPU backend of FFTs and associated gather/scatter
  operations.
- Add `-D__OFFLOAD_UNIFIED_MEMORY` to enable unified memory support (experimental and only supports
  Mi250X and above)
- Add `GPUVER=Mi50, Mi60, Mi100, Mi250`
- Add `OFFLOAD_CC = hipcc`
- Add `-lamdhip64` to the `LIBS` variable
- Add
  `OFFLOAD_FLAGS = '-munsafe-fp-atomics -fopenmp -m64 -pthread -fPIC -D__GRID_HIP -O2 --offload-arch=gfx908 --rocm-path=$(ROCM_PATH)'`
  where `ROCM_PATH` is the path where the rocm sdk resides. Architectures Mi300(A,X) (gfx1103),
  Mi250 (gfx90a), Mi100 (gfx908), Mi50 (gfx906) the hip backend for the grid library supports nvidia
  hardware as well. It uses the same code and can be used to validate the backend in case of access
  to Nvidia hardware only. To get the compilation working, follow the steps above and set the
  `OFFLOAD_FLAGS` with right `nvcc` parameters (see the cuda section of this document). The
  environment variable `HIP_PLATFORM` should be set to `HIP_PLATFORM=nvidia` to indicate to hipcc to
  use the nvcc compiler instead.
- Specify the C++ compiler (e.g., `CXX = g++`). Remember to set the CXXFLAGS flags to support C++11
  standard and OpenMP.
- When the HIP backend is enabled for DBCSR using `-D__DBCSR_ACC`, then add `-D__HIP_PLATFORM_AMD__`
  to `CXXFLAGS` and set `OFFLOAD_TARGET = hip`.
- Use `-D__OFFLOAD_PROFILING` to turn on the AMD ROC TX and Tracer libray. It requires to link
  `-lroctx64 -lroctracer64`.

### 2x. OpenCL Devices

OpenCL devices are currently supported for DBCSR and DBM/DBT, and can cover GPUs and other devices.
Kernels can be automatically tuned.

Note: the OpenCL backend uses some functionality from LIBXSMM (dependency). CP2K's offload-library
serving DBM/DBT and other libraries depends on DBCSR's OpenCL backend.

- Installing OpenCL and preparing the runtime environment
  - Installing an OpenCL runtime depends on the operating system and the device vendor. Debian for
    instance brings two packages called `opencl-headers` and `ocl-icd-opencl-dev` which can be
    present in addition to a vendor-specific installation. The OpenCL header files are only
    necessary if CP2K/DBCSR is compiled from source. Please note, some implementations ship with
    outdated OpenCL headers which can prevent using latest features (if an application discovers
    such features only at compile-time). When building from source, for instance `libOpenCL.so` is
    sufficient at link-time (ICD loader). However, an Installable Client Driver (ICD) is finally
    necessary at runtime.
  - Nvidia CUDA, AMD HIP, and Intel OneAPI are fully equipped with an OpenCL runtime (if
    `opencl-headers` package is not installed, CPATH can be needed to point into the former
    installation, similarly `LIBRARY_PATH` for finding `libOpenCL.so` at link-time). Installing a
    minimal or stand-alone OpenCL is also possible, e.g., following the instructions for Debian (or
    Ubuntu) as given for every [release](https://github.com/intel/compute-runtime/releases) of the
    [Intel Compute Runtime](https://github.com/intel/compute-runtime).
  - The environment variable `ACC_OPENCL_VERBOSE` prints information at runtime of CP2K about
    kernels generated (`ACC_OPENCL_VERBOSE=2`) or executed (`ACC_OPENCL_VERBOSE=3`) which can be
    used to check an installation.
- Building CP2K with OpenCL-based DBCSR
  - CP2K's toolchain supports `--enable-opencl` to select DBCSR's OpenCL backend. This can be
    combined with `--enable-cuda` (`--gpu-ver` is then imposed) to use a GPU for CP2K's GRID and PW
    components (no OpenCL support yet) with DBM's CUDA implementation to be preferred.
  - For manually writing an ARCH-file, add `-D__OPENCL` and `-D__DBCSR_ACC` to `CFLAGS` and add
    `-lOpenCL` to the `LIBS` variable, i.e., `OFFLOAD_CC` and `OFFLOAD_FLAGS` can duplicate `CC` and
    `CFLAGS` (no special offload compiler needed). Please also set `OFFLOAD_TARGET = opencl` to
    enable the OpenCL backend in DBCSR. For OpenCL, it is not necessary to specify a GPU version
    (e.g., `GPUVER = V100` would map/limit to
    `exts/dbcsr/src/acc/opencl/smm/params/tune_multiply_V100.csv`). In fact, `GPUVER` limits tuned
    parameters to the specified GPU, whereas by default all tuned parameters are embedded
    (`exts/dbcsr/src/acc/opencl/smm/params/*.csv`) and applied at runtime. If auto-tuned parameters
    are not available for DBCSR, well-chosen defaults will be used to populate kernels at runtime.
  - Auto-tuned parameters are embedded into the binary, i.e., CP2K does not rely on a hard-coded
    location. Setting `OPENCL_LIBSMM_SMM_PARAMS=/path/to/csv-file` environment variable can supply
    parameters for an already built application, or `OPENCL_LIBSMM_SMM_PARAMS=0` can disable using
    tuned parameters. Refer to <https://cp2k.github.io/dbcsr/> on how to tune kernels (parameters).
- Building CP2K with OpenCL-based DBM library
  - For manually writing an ARCH-file, add `-D__OFFLOAD_OPENCL` to `CFLAGS` in addition to following
    above instructions for "Building CP2K with OpenCL-based DBCSR". An additional Makefile rule can
    be necessary to transform OpenCL code into a ressource header file.

### 2y. matrix-matrix multiplication offloading on GPU using SPLA

The SPLA library is a hard dependency of SIRIUS but can also be used as a standalone library. It
provides a generic interface to the blas gemm family with offloading on GPU. Offloading supports
both CUDA and ROCM.

To make the functionality available, add the flag `-D__SPLA -D__OFFLOAD_GEMM` to the `DFLAGS`
variable and compile SPLA with Fortran interface and GPU support. Please note that only the
functions replacing the dgemm calls with `offload_dgemm` will eventually be offloaded to the GPU.
The SPLA library has internal criteria to decide if it is worth to do the operation on GPU or not.
Calls to `offload_dgemm` also accept pointers on GPU or a combination of them.

### 2y. libgrpp (optional, enables calculations with ECPs)

- libgrpp is a library for the calculation of integrals with GTOs and ECPs
- The libgrpp library can be found under <https://github.com/aoleynichenko/libgrpp>
- During the installation, the directories `$(LIBGRPP_DIR)/lib` and `$(LIBGRPP_DIR)/include` are
  created.
- Add `-D__LIBGRPP` to DFLAGS, `-I$(LIBGRPP_DIR)/include` to FCFLAGS and
  `-L$(LIBGRPP_DIR)/lib -llibgrpp` to LIBS

<!---
### 2y. LibMaxwell (External Maxwell Solver)

- LibMaxwell is a library to solve the time-dependent Maxwell equations
  and use the resulting electric field in MD runs or real-time propagation.
- Add `-D__LIBMAXWELL` to DFLAGS to enable support for LibMaxwell.
- See <https://brehm-research.de> for more information.
-->

### 2y. DeePMD-kit (optional, wider range of interaction potentials)

DeePMD-kit - Deep Potential Molecular Dynamics. Support for DeePMD-kit can be enabled via the flag
`-D__DEEPMD`.

- DeePMD-kit C interface can be downloaded from
  <https://docs.deepmodeling.com/projects/deepmd/en/master/install/install-from-c-library.html>
- For more information see <https://github.com/deepmodeling/deepmd-kit.git>.

### 2z. DFTD4 (optional, dispersion correction)

- dftd4 - Generally Applicable Atomic-Charge Dependent London Dispersion Correction.
- For more information see <https://github.com/dftd4/dftd4>
- Add `-D__DFTD4` to DFLAGS, `-ldftd4 -lmstore -lmulticharge -lmctc-lib` to LIBS and
  `-I'${DFTD4_DFTD4}/../..' -I'${DFTD4_DFTD4}' -I'${DFTD4_MCTC}'` to CFLAGS

### 2y. libsmeagol (optional, electron transport calculation with current-induced forces)

- libsmeagol is an external library to compute electron transport properties using Non-Equilibrium
  Green Functions (NEGF) method. The library can be downloaded from
  <https://github.com/StefanoSanvitoGroup/libsmeagol>.
- libsmeagol depends on an MPI library and can only be linked with MPI parallel CP2K binaries.
- During the installation, the directories `$(LIBSMEAGOL_DIR)/lib` and `$(LIBGRPP_DIR)/obj` are
  created.
- Add `-D__SMEAGOL` to DFLAGS, `-I$(LIBSMEAGOL_DIR)/obj` to FCFLAGS and
  `-L$(LIBSMEAGOL_DIR)/lib -lsmeagol` to LIBS

### 2z. TREXIO (optional, unified computational chemistry format)

TREXIO - Open-source file format and library. Support for TREXIO can be enabled via the flag
`-D__TREXIO`.

- TREXIO library can be downloaded from <https://github.com/trex-coe/trexio>
- For more information see <https://trex-coe.github.io/trexio/index.html>.

## 3. Compile

### 3a. ARCH files

The location of compiler and libraries needs to be specified. Examples for several common
architectures can be found in [arch folder](./arch/). The names of these files match
`architecture.version` e.g., [Linux-gnu-x86_64.psmp](./arch/Linux-gnu-x86_64.psmp). Alternatively,
<https://dashboard.cp2k.org> provides sample arch files as part of the testing reports (click on the
status field, search for 'ARCH-file').

Conventionally, there are six versions:

| VERSION | Meaning                          |
| ------- | -------------------------------- |
| sdbg    | OpenMP + debug settings          |
| sopt    | OpenMP + OMP_NUM_THREADS=1       |
| ssmp    | OpenMP                           |
| pdbg    | MPI + OpenMP + debug settings    |
| popt    | MPI + OpenMP + OMP_NUM_THREADS=1 |
| psmp    | MPI + OpenMP                     |

You'll need to modify one of these files to match your system's settings.

Some architecture files like the file [Linux-gnu-x86_64.psmp](./arch/Linux-gnu-x86_64.psmp) are
sourceable (see above), i.e.

```shell
source arch/Linux-gnu-x86_64.psmp     
```

will launch a build of the CP2K toolchain which will build all dependencies needed for compiling
CP2K. Building a `psmp` version will also create a `popt` CP2K binary and vice versa. The same is
true for the `ssmp` and `sopt` versions of CP2K.

You can now build CP2K using these settings (where -j N allows for a parallel build using N
processes):

```shell
make -j N ARCH=architecture VERSION=version
```

e.g.

```shell
make -j N ARCH=Linux-gnu-x86_64 VERSION=psmp
```

A CP2K binary should appear in the `./exe/ARCH/` folder.

All compiled files, libraries, binaries, etc. of all architectures and versions can be removed with

```shell
make distclean
```

To remove only `*.o` and `*.mod` files (i.e. keep CP2K binaries in exe) for a given ARCH/VERSION use

```shell
make ARCH=Linux-gnu-x86_64 VERSION=psmp clean
```

To remove everything for a given ARCH/VERSION use

```shell
make ARCH=Linux-gnu-x86_64 VERSION=psmp realclean
```

### 3b. Compilation Flags

The following flags should be present (or not) in the arch file, partially depending on installed
libraries (see 2.)

- `-D__parallel` builds an MPI parallel CP2K binary (implies the use and thus the availabiltity of
  the ScaLAPACK/BLACS libraries)
- `-D__LIBINT` use LIBINT (needed for HF exchange)
- `-D__LIBXC` use LIBXC
- `-D__LIBGRPP` use libgrpp (for calculations with ECPs)
- `-D__ELPA` use ELPA in place of SYEVD to solve the eigenvalue problem
- `-D__FFTW3` FFTW version 3 is recommended
- `-D__MKL` link the MKL library for linear algebra and/or FFT
- `-D__GRID_CORE=X` (with X=1..6) specific optimized core routines can be selected. Reasonable
  defaults are [provided](./src/grid/collocate_fast.f90) but trial-and-error might yield (a small
  ~10%) speedup.
- `-D__PILAENV_BLOCKSIZE`: can be used to specify the blocksize (e.g.,
  `-D__PILAENV_BLOCKSIZE=1024`), which is a hack to overwrite (if the linker allows this) the
  PILAENV function provided by ScaLAPACK. This can lead to much improved PDGEMM performance. The
  optimal value depends on hardware (GPU?) and precise problem. Alternatively, Cray provides an
  environment variable to this effect (e.g., `export LIBSCI_ACC_PILAENV=4000`)
- `-D__STATM_RESIDENT` or `-D__STATM_TOTAL` toggles memory usage reporting between resident memory
  and total memory
- `-D__CRAY_PM_ACCEL_ENERGY` or `-D__CRAY_PM_ENERGY` switch on energy profiling on Cray systems
- `-D__NO_ABORT` to avoid calling abort, but STOP instead (useful for coverage testing, and to avoid
  core dumps on some systems)
- `-D__HDF5` enables hdf5 support. This is a hard dependency for SIRIUS and TREXIO, but can also be
  used by itself to allow read/write functionalities of QCSchema files in the active space module
- `-D__TREXIO` enables TREXIO I/O support

Features useful to deal with legacy systems

- `-D__NO_MPI_THREAD_SUPPORT_CHECK` - Workaround for MPI libraries that do not declare they are
  thread safe (serialized).
- `-D__NO_SOCKETS` disables the socket interface in case of troubles compiling on systems that do
  not support POSIX sockets.
- `-D__HAS_IEEE_EXCEPTIONS` disables trapping temporarily for libraries like ScaLAPACK.
- The Makefile automatically compiles in the path to the data directory via the flag `-D__DATA_DIR`.
  If you want to compile in a different path, set the variable `DATA_DIR` in your arch-file.
- `-D__NO_STATM_ACCESS` - Do not try to read from /proc/self/statm to get memory usage information.
  This is otherwise attempted on several. Linux-based architectures or using with the NAG, gfortran,
  compilers.
- `-D__CHECK_DIAG` Debug option which activates an orthonormality check of the eigenvectors
  calculated by the selected eigensolver

### 3c. Building CP2K as a library

You can build CP2K for use as a library by adding `libcp2k` as an option to your `make` command,
e.g.

```shell
make -j N ARCH=Linux-gnu-x86_64 VERSION=psmp libcp2k
```

This will create `libcp2k.a` in the relevant subdirectory of `./lib/`. You will need to add this
subdirectory to the library search path of your compiler (typically via the `LD_LIBRARY_PATH`
environment variable or the `-L` option to your compiler) and link to the library itself with
`-lcp2k`.

In order to use the functions in the library you will also require the `libcp2k.h` header file. This
can be found in `./src/start/` directory. You should add this directory to the header search path of
your compiler (typically via the `CPATH` environment variable or the `-I` option to your compiler).

For Fortran users, you will require the module interface file (`.mod` file) for every MODULE
encountered in the source. These are compiler specific and are to be found in the subdirectory of
`./obj/` that corresponds to your build, e.g.,

```shell
./obj/Linux-gnu-x86_64/psmp/
```

In order for your compiler to find these, you will need to indicate their location to the compiler
as is done for header files (typically via the `CPATH` environment variable or the `-I` option to
your compiler).

## 4. If it doesn't work

If things fail, take a break... go back to 2a (or skip to step 6).

## 5. Regtesting

If compilation works fine, it is recommended to test the generated binary, to exclude errors in
libraries, or miscompilations, etc.

```shell
make -j ARCH=... VERSION=... test
```

should work if you can locally execute CP2K without the need for, e.g., batch submission.

In the other case, you might need to configure the underlying testing script as described more
systematically at <https://www.cp2k.org/dev:regtesting>

## 6. Talk to us

In any case please tell us your comments, praise, criticism, thanks, etc. see
<https://www.cp2k.org>.

## 7. Manual

A reference manual of CP2K can be found on the web: <https://manual.cp2k.org> or can be generated
using the cp2k binary, see <https://manual.cp2k.org/trunk/generate_manual_howto.html>

## 8. Happy computing

The CP2K team.
