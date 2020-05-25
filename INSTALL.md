# How to compile the CP2K code

##  1. Acquire the code


For users, the preferred method is to [download a release](https://github.com/cp2k/cp2k/releases/).
For developers, the preferred method is to [download from Git](./README.md#downloading-cp2k-source-code).

For more details on downloading CP2K, see  https://www.cp2k.org/download.

## 2. Install prerequisites

The most convenient way to install pre-requisites is by using the [toolchain script](./tools/toolchain/install_cp2k_toolchain.sh).

For a complete introduction to the toolchain script, see the [README for users](./tools/toolchain/README_FOR_USERS.md) or the [README for developers](./tools/toolchain/README_FOR_DEVELOPERS.md).

The basic steps are:

- Read toolchain installation options:

```
> cd tools/toolchain/
> ./install_cp2k_toolchain.sh --help
```

- Launch toolchain script (example option choice)

```
> ./install_cp2k_toolchain.sh --with-libxsmm=install --with-openblas=system \
     --with-fftw=system --with-reflapack=no  --enable-cuda --enable-omp
```

- Once the script has completed successfully, follow the instructions given at the end of its output.
Note that the pre-built arch files provided by the toolchain are for the GNU compiler, users have to adapt them for other compilers. It is possible to use the provided [arch files](./arch) as guidance.

Sub-points here discuss prerequisites needed to build CP2K. Copies of the recommended versions of 3rd party software can be downloaded from https://www.cp2k.org/static/downloads/.

### 2a. GNU make (required, build system)

GNU make should be on your system (gmake or make on linux) and used for the build, go to https://www.gnu.org/software/make/make.html download from https://ftp.gnu.org/pub/gnu/make/

### 2b. Python (required, build system)
Python 3.5+ is needed to run the dependency generator. On most system Python is already installed. For more information visit: https://www.python.org/

### 2c. Fortran and C Compiler (required, build system)
A Fortran 2008 compiler and matching C99 compiler should be installed on your system. We have good experience with gcc/gfortran (gcc >=4.6 works, later version recommended). Be aware that some compilers have bugs that might cause them to fail (internal compiler errors, segfaults) or, worse, yield a mis-compiled CP2K. Report bugs to compiler vendors; they (and we) have an interest in fixing them. A list of tested compiler can be found [here](https://www.cp2k.org/dev:compiler_support). Always run a `make -j test` (See point 5.) after compilation to identify these problems.

### 2d. BLAS and LAPACK (required, base functionality)
BLAS and LAPACK should be installed.  Using vendor-provided libraries can make a very significant difference (up to 100%, e.g., ACML, MKL, ESSL), not all optimized libraries are bug free. Use the latest versions available, use the interfaces matching your compiler, and download all patches!

  * The canonical BLAS and LAPACK can be obtained from the Netlib repository:
    * http://www.netlib.org/blas/
    * http://www.netlib.org/lapack/ and see also
    * http://www.netlib.org/lapack-dev/
  * Open fast alternatives, include:
    * http://www.openblas.net/
    * http://math-atlas.sourceforge.net/
    * https://www.tacc.utexas.edu/research-development/tacc-software/gotoblas2

If compiling with OpenMP support then it is recommended to use a non-threaded version of BLAS. In particular if compiling with MKL and using OpenMP you must define `-D__MKL` to ensure the code is thread-safe. MKL with multiple OpenMP threads in CP2K requires that CP2K was compiled with the Intel compiler. If the `cpp` precompiler is used in a separate precompilation step in combination with the Intel Fortran compiler, `-D__INTEL_COMPILER` must be added explicitly (the Intel compiler sets `__INTEL_COMPILER` otherwise automatically).

On the Mac, BLAS and LAPACK may be provided by Apple's Accelerate framework. If using this framework, `-D__ACCELERATE` must be defined to account for some interface incompatibilities between Accelerate and reference BLAS/LAPACK.

When building on/for Windows using the Minimalist GNU for Windows (MinGW) environment, you must set `-D__MINGW`,  `-D__NO_STATM_ACCESS` and `-D__NO_IPI_DRIVER` to avoid undefined references during linking, respectively errors while printing the statistics.

### 2e. MPI and SCALAPACK (optional, required for MPI parallel builds)
MPI (version 2) and SCALAPACK are needed for parallel code. (Use the latest versions available and download all patches!).

:warning: Note that your MPI installation must match the used Fortran compiler. If your computing platform does not provide MPI, there are several freely available alternatives:

  * MPICH2 MPI:  http://www-unix.mcs.anl.gov/mpi/mpich/
  * OpenMPI MPI: http://www.open-mpi.org/
  * ScaLAPACK:
    * http://www.netlib.org/scalapack/
    * http://www.netlib.org/lapack-dev/
    * ScaLAPACK can be part of ACML or cluster MKL. These libraries are recommended if available.
    * Recently a [ScaLAPACK installer](http://www.netlib.org/scalapack/scalapack_installer.tgz) has been added that simplifies the installation.

CP2K assumes that the MPI library implements MPI version 3. If you have an older version of MPI (e.g. MPI 2.0) available you must define `-D__MPI_VERSION=2` in the arch file.

### 2f. FFTW (optional, improved performance of FFTs)
FFTW can be used to improve FFT speed on a wide range of architectures. It is strongly recommended to install and use FFTW3. The current version of CP2K works with FFTW 3.X (use `-D__FFTW3`). It can be downloaded from http://www.fftw.org/

:warning: Note that FFTW must know the Fortran compiler you will use in order to install properly (e.g., `export F77=gfortran` before configure if you intend to use gfortran).

:warning: Note that on machines and compilers which support SSE you can configure FFTW3 with `--enable-sse2`. Compilers/systems that do not align memory (NAG f95, Intel IA32/gfortran) should either not use `--enable-sse2` or otherwise set the define `-D__FFTW3_UNALIGNED` in the arch file. When building an OpenMP parallel version of CP2K (ssmp or psmp), the FFTW3 threading library libfftw3_threads (or libfftw3_omp) is required.

### 2g. LIBINT (optional, enables methods including HF exchange)
  * Hartree-Fock exchange (optional, use `-D__LIBINT`) requires the libint package to be installed.
  * Recommended way to build libint: Download a CP2K-configured libint library from [libint-cp2k](https://github.com/cp2k/libint-cp2k). Build and install libint by following the instructions provided there. Note that using a library configured for higher maximum angular momentum will increase build time and binary size of CP2K executable (assuming static linking).
  * CP2K is not hardwired to these provided libraries and any other libint library (version >= 2.5.0) should be compatible as long as it was compiled with `--enable-eri=1` and default ordering.
  * Avoid debugging information (`-g` flag) for compiling libint since this will increase library size by a large factor.
  * In the arch file of CP2K: add `-D__LIBINT` to the `DFLAGS`. Add `-L$(LIBINT_DIR)/lib -lint2 -lstdc++` to `LIBS` and `-I$(LIBINT_DIR)/include` to `FCFLAGS`. `lstdc++` is needed if you use the GNU C++ compiler.
  * Libint 1 is no longer supported and the previously needed flags `-D__LIBINT_MAX_AM` and `-D__LIBDERIV_MAX_AM1` are ignored.
  * `-D__MAX_CONTR=4` (default=2) can be used to compile efficient contraction kernels up to l=4, but the build time will increase accordingly.

### 2h. libsmm (optional, improved performance for matrix multiplication)
  * A library for small matrix multiplies can be built from the included source (see exts/dbcsr/tools/build_libsmm/README).  Usually only the double precision real and perhaps complex is needed.  Link to the generated libraries. For a couple of architectures prebuilt libsmm are available at https://www.cp2k.org/static/downloads/libsmm/.
  * Add `-D__HAS_smm_dnn` to the defines to make the code use the double precision real library.  Similarly use `-D__HAS_smm_snn` for single precision real and `-D__HAS_smm_znn` / `-D__HAS_smm_cnn` for double / single precision complex.
  * Add `-D__HAS_smm_vec` to enable the new vectorized interfaces of libsmm.

### 2i. libxsmm (optional, improved performance for matrix multiplication)
  * A library for matrix operations and deep learning primitives: https://github.com/hfp/libxsmm/
  * Add `-D__LIBXSMM` to enable it, with suitable include and library paths, e.g. `FCFLAGS += -I${LIBXSMM_DIR}/include -D__LIBXSMM` and `LIBS += -L${LIBXSMM_DIR}/lib -lxsmmf -lxsmm -ldl`

### 2j. CUDA (optional, improved performance on GPU systems)
  * Specify NVCC (e.g. `NVCC = nvcc`) and NVFLAGS (e.g. `NVFLAGS = -O3 -g -w --std=c++11`) variables.
  * `-D__ACC` needed to enable accelerator support.
  * Use the `-D__DBCSR_ACC` to enable accelerator support for matrix multiplications.
  * Add `-lstdc++ -lcudart -lnvrtc -lcuda -lcublas` to LIBS.
  * Specify the GPU type (e.g. `GPUVER   = P100`), possible values are K20X, K40, K80, P100, V100.
  * Specify the C++ compiler (e.g. `CXX = g++`). Remember to set the flags to support C++11 standard.
  * Use `-D__PW_CUDA` for CUDA support for PW (gather/scatter/fft) calculations.
  * CUFFT 7.0 has a known bug and is therefore disabled by default. NVIDIA's webpage list a patch (an upgraded version cufft i.e. >= 7.0.35) - use this together with `-D__HAS_PATCHED_CUFFT_70`.
  * Use `-D__CUDA_PROFILING` to turn on Nvidia Tools Extensions. It requires to link `-lnvToolsExt`.
  * Link to a blas/scalapack library that accelerates large DGEMMs (e.g. libsci_acc)

### 2k. libxc (optional, wider choice of xc functionals)
  * The version 4.0.3 (or later) of libxc can be downloaded from http://www.tddft.org/programs/octopus/wiki/index.php/Libxc.
  * During the installation, the directories `$(LIBXC_DIR)/lib` and `$(LIBXC_DIR)/include` are created.
  * Add `-D__LIBXC` to DFLAGS, `-I$(LIBXC_DIR)/include` to FCFLAGS and `-L$(LIBXC_DIR)/lib -lxcf03 -lxc` to LIBS.
  * :warning: Note that the deprecated flags `-D__LIBXC2` and `-D__LIBXC3` are ignored.

### 2l. ELPA (optional, improved performance for diagonalization)
Library ELPA for the solution of the eigenvalue problem
  * ELPA replaces the ScaLapack `SYEVD` to improve the performance of the diagonalization
  * A version of ELPA can to be downloaded from http://elpa.rzg.mpg.de/software.
  * During the installation the `libelpa.a` (or `libelpa_openmp.a` if OpenMP is enabled) is created.
  * Minimal supported version of ELPA is 2018.05.001.
  * Add `-D__ELPA` to `DFLAGS`
  * Add `-I$(ELPA_INCLUDE_DIR)/modules` to `FCFLAGS`
  * Add `-I$(ELPA_INCLUDE_DIR)/elpa` to `FCFLAGS`
  * Add `-L$(ELPA_DIR)` to `LDFLAGS`
  * Add `-lelpa` to `LIBS`
  * For specific architectures it can be better to install specifically
    optimized kernels (see BG) and/or employ a higher optimization level to compile it.

### 2m. PEXSI (optional, low scaling SCF method)
The Pole EXpansion and Selected Inversion (PEXSI) method requires the PEXSI library and two dependencies (ParMETIS or PT-Scotch and SuperLU_DIST).
  * Download PEXSI (www.pexsi.org) and install it and its dependencies by following its README.md.
  * PEXSI versions 0.10.x have been tested with CP2K. Older versions are not supported.
  * PEXSI needs to be built with `make finstall`.

In the arch file of CP2K:
  * Add `-lpexsi_${SUFFIX} -llapack -lblas -lsuperlu_dist_3.3 -lparmetis -lmetis`, and their paths (with `-L$(LIB_DIR)`) to LIBS.
  * It is important that a copy of LAPACK and BLAS is placed before and after these libraries  (replace `-llapack` and `-lblas` with the optimized versions as needed).
  * In order to link in PT-Scotch instead of ParMETIS replace `-lparmetis -lmetis` with: `-lptscotchparmetis -lptscotch -lptscotcherr -lscotchmetis -lscotch -lscotcherr`
  * Add `-I$(PEXSI_DIR)/fortran/` to FCFLAGS.
  * Add `-D__LIBPEXSI` to DFLAGS.

Below are some additional hints that may help in the compilation process:
  * For building PT-Scotch, the flag `-DSCOTCH_METIS_PREFIX` in `Makefile.inc` must not be set and the flag `-DSCOTCH_PTHREAD` must be removed.
  * For building SuperLU_DIST with PT-Scotch, you must set the following in `make.inc`:

```
METISLIB = -lscotchmetis -lscotch -lscotcherr
PARMETISLIB = -lptscotchparmetis -lptscotch -lptscotcherr
```

### 2n. QUIP (optional, wider range of interaction potentials)
QUIP - QUantum mechanics and Interatomic Potentials Support for QUIP can be enabled via the flag `-D__QUIP`.

For more information see http://www.libatoms.org/ .

### 2o. PLUMED (optional, enables various enhanced sampling methods)
CP2K can be compiled with PLUMED 2.x (`-D__PLUMED2`).

See https://cp2k.org/howto:install_with_plumed for full instructions.

### 2p. spglib (optional, crystal symmetries tools)
A library for finding and handling crystal symmetries
  * The spglib can be downloaded from https://github.com/atztogo/spglib
  * For building CP2K with the spglib add `-D__SPGLIB` to DFLAGS

### 2q. SIRIUS (optional, plane wave calculations)
SIRIUS is a domain specific library for electronic structure calculations.
  * The code is available at https://github.com/electronic-structure/SIRIUS
  * For building CP2K with SIRIUS add `-D__SIRIUS` to DFLAGS.
  * See https://electronic-structure.github.io/SIRIUS/ for more information.

### 2r. FPGA (optional, plane wave FFT calculations)
  * Use `-D__PW_FPGA` to enable FPGA support for PW (fft) calculations. Currently tested only for Intel Stratix 10 and Arria 10 GX1150 FPGAs.
  * Supports single precision and double precision fft calculations with the use of dedicated APIs.
  * Double precision is the default API chosen when set using the `-D__PW_FPGA` flag.
  * Single precision can be set using an additional `-D__PW_FPGA_SP` flag along with the `-D__PW_FPGA` flag.
  * Kernel code has to be synthesized separately and copied to a specific location.
  * See https://github.com/pc2/fft3d-fpga for the kernel code and instructions for synthesis.
  * Read `src/pw/fpga/README.md` for information on the specific location to copy the binaries to.
  * Currently supported FFT3d sizes - 16^3, 32^3, 64^3.
  * Include aocl compile flags and `-D__PW_FPGA -D__PW_FPGA_SP` to `CFLAGS`, aocl linker flags to `LDFLAGS` and aocl libs to `LIBS`.
  * CUDA and FPGA are mutually exclusive. Building with both `__PW_CUDA` and  `__PW_FPGA` will throw a compilation error.

### 2s. COSMA (Distributed Communication-Optimal Matrix-Matrix Multiplication Algorithm)
  * COSMA is a replacement of the pdgemm routine included in scalapack. The
    library supports both CPU and GPUs. No specific flag during compilation is
    needed to use the library in cp2k, excepted during linking time where the
    library should be placed in front of the scalapack library.
  * see https://github.com/eth-cscs/COSMA for more information.

## 3. Compile

### 3a. ARCH files
The location of compiler and libraries needs to be specified. Examples for a number of common architectures examples can be found in [arch folder](./arch/). The names of these files match `architecture.version` e.g., [Linux-x86-64-gfortran.sopt](./arch/Linux-x86-64-gfortran.sopt). Alternatively https://dashboard.cp2k.org/ provides sample arch files as part of the testing reports (click on the status field, search for 'ARCH-file').
  * With -DNDEBUG assertions may be stripped ("compiled out").
  * NDEBUG is the ANSI-conforming symbol name (not __NDEBUG).
  * Regular release builds may carry assertions for safety.

Conventionally, there are six versions:

| Acronym |        Meaning          |         Recommended for            |
|---------|-------------------------|------------------------------------|
| sdbg    | serial                  | single core testing and debugging  |
| sopt    | serial                  | general single core usage          |
| ssmp    | parallel (only OpenMP)  | optimized, single node, multi core |
| pdbg    | parallel (only MPI)     | multi-node testing and debugging   |
| popt    | parallel (only MPI)     | general usage, no threads          |
| psmp    | parallel (MPI + OpenMP) | general usage, threading might improve scalability and memory usage |

You'll need to modify one of these files to match your system's settings.

You can now build CP2K using these settings (where -j N allows for a parallel build using N processes):
```
> make -j N ARCH=architecture VERSION=version
```
e.g.
```
> make -j N ARCH=Linux-x86-64-gfortran VERSION=sopt
```
as a short-cut, you can build several version of the code at once
```
> make -j N ARCH=Linux-x86-64-gfortran VERSION="sopt popt ssmp psmp"
```
An executable should appear in the `./exe/` folder.

All compiled files, libraries, executables, .. of all architectures and versions can be removed with
```
> make distclean
```
To remove only objects and mod files (i.e., keep exe) for a given ARCH/VERSION use, e.g.,
```
> make ARCH=Linux-x86-64-gfortran VERSION=sopt clean
```
to remove everything for a given ARCH/VERSION use, e.g.,
```
> make ARCH=Linux-x86-64-gfortran VERSION=sopt realclean
```

### 3b. Compilation Flags

The following flags should be present (or not) in the arch file, partially depending on installed libraries (see 2.)
  * `-D__parallel -D__SCALAPACK` parallel runs
  * `-D__LIBINT` use libint (needed for HF exchange)
  * `-D__LIBXC` use libxc
  * `-D__ELPA` use ELPA in place of SYEVD  to solve the eigenvalue problem
  * `-D__FFTW3` FFTW version 3 is recommended
  * `-D__PW_CUDA` CUDA FFT and associated gather/scatter on the GPU
  * `-D__MKL` link the MKL library for linear algebra and/or FFT

  * with `-D__GRID_CORE=X` (with X=1..6) specific optimized core routines can be selected.  Reasonable defaults are [provided](./src/grid/collocate_fast.f90) but trial-and-error might yield (a small ~10%) speedup.
  * with `-D__HAS_LIBGRID` (and `-L/path/to/libgrid.a` in LIBS) tuned versions of integrate and collocate routines can be [generated](./tools/autotune_grid/README).
  * `-D__PILAENV_BLOCKSIZE`: can be used to specify the blocksize (e.g. `-D__PILAENV_BLOCKSIZE=1024`), which is a hack to overwrite (if the linker allows this) the PILAENV function provided by Scalapack. This can lead to much improved PDGEMM performance. The optimal value depends on hardware (GPU?) and precise problem. Alternatively, Cray provides an environment variable to this effect (e.g. `export LIBSCI_ACC_PILAENV=4000`)
  * `-D__STATM_RESIDENT` or `-D__STATM_TOTAL` toggles memory usage reporting between resident memory and total memory
  * `-D__CRAY_PM_ACCEL_ENERGY` or `-D__CRAY_PM_ENERGY` switch on energy profiling on Cray systems
  * `-D__NO_ABORT` to avoid calling abort, but STOP instead (useful for coverage testing, and to avoid core dumps on some systems)

Features useful to deal with legacy systems
  * `-D__NO_MPI_THREAD_SUPPORT_CHECK`  - Workaround for MPI libraries that do not declare they are thread safe (funneled) but you want to use them with OpenMP code anyways.
  * `-D__NO_IPI_DRIVER` disables the socket interface in case of troubles compiling on systems that do not support POSIX sockets.
  * `-D__HAS_IEEE_EXCEPTIONS` disables trapping temporarily for libraries like scalapack.
  * The Makefile automatically compiles in the path to the data directory via the flag `-D__DATA_DIR`. If you want to compile in a different path, set the variable `DATA_DIR` in your arch-file.
  * `-D__NO_STATM_ACCESS` - Do not try to read from /proc/self/statm to get memory usage information. This is otherwise attempted on several. Linux-based architectures or using with the NAG, gfortran, compilers.
  * `-D__CHECK_DIAG` Debug option which activates an orthonormality check of the eigenvectors calculated by the selected eigensolver

### 3c. Building CP2K as a library

You can build CP2K for use as a library by adding `libcp2k` as an option to your `make` command, e.g.
```
> make -j N ARCH=Linux-x86-64-gfortran VERSION=sopt libcp2k
```
This will create `libcp2k.a` in the relevant subdirectory of `./lib/`. You will need to add this subdirectory to the library search path of your compiler (typically via the `LD_LIBRARY_PATH` environment variable or the `-L` option to your compiler) and link to the library itself with `-lcp2k`.

In order to use the functions in the library you will also require the `libcp2k.h` header file. This can be found in `./src/start/` directory. You should add this directory to the header search path of your compiler (typically via the `CPATH` environment variable or the `-I` option to your compiler).

For Fortran users, you will require the module interface file (`.mod` file) for every MODULE encountered in the source. These are compiler specific and are to be found in the subdirectory of `./obj/` that corresponds to your build, e.g.,
```
./obj/Linux-x86-64-gfortran/sopt/
```
In order for your compiler to find these, you will need to indicate their location to the compiler as is done for header files  (typically via the `CPATH` environment variable or the `-I` option to your compiler).

## 4. If it doesn't work?
If things fail, take a break... go back to 2a (or skip to step 6).

## 5. Regtesting

If compilation works fine, it is recommended to test the generated binary, to exclude errors in libraries, or miscompilations, etc.
```
make -j ARCH=... VERSION=... test
```

should work if you can locally execute CP2K without the need for e.g. batch submission.

In the other case, you might need to configure the underlying testing script as described more systematically at https://www.cp2k.org/dev:regtesting

## 6. Talk to us
In any case please tell us your comments, praise, criticism, thanks,... see https://www.cp2k.org/

## 7. Manual
A reference manual of CP2K can be found on the web: https://manual.cp2k.org/ or can be generated using the cp2k executable, see https://manual.cp2k.org/trunk/generate_manual_howto.html

## 8. Happy computing!

 The CP2K team.
