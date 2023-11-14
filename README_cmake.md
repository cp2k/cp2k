# Build CP2K with CMake

This document regroups information about the CP2K CMake system. CMake is used to detect CP2K
dependencies and configure the compilation process. Dependencies should be installed independently
either with a distribution package manager, easybuild, or spack to name a few.

It is easier to build and install all manually built dependencies in a single directory ideally
where CP2K will also be installed. CMake will have less difficulties to find the `FindPACKAGE.cmake`
files and dependent libraries. CMake will also use environment variables such as `ORNL_FFTW3_ROOT`,
etc. Usually a standard prefix is used in HPC environments. If known just add it in the
`cmake/cp2k_utils.cmake` file.

The CMake build system requires a minimum set of dependencies:

- a c, C++, and fortran compiler (gcc, intel oneapi, AMD or nvidia SDK, xlf, etc...)
- an MPI implementation
- [DBCSR](https://cp2k.github.io/dbcsr/develop/).
- openmp
- any flavor of BLAS, LAPACK, SCALAPACK.
- fftw3
- CMake

Major vendors implementations of BLAS, LAPACK, and scalapack are supported. The build system was
tested with MKL, cray libsci, OpenBLAS, flexiblas, blis or ATLAS.

All optional dependencies are turned off by default. DBCSR is the only dependency that can be built
at the same time than CP2K. It is OFF by default. We assume that mpi, blas, lapack, and scalapack
are installed. To compile `CP2K and `DBCSR\` open a terminal and write

```shell
cmake -DCMAKE_INSTALL_PREFIX=/pyprefix -DCP2K_BUILD_DBCSR=ON ..
```

This command will build cp2k and dbcsr dependencies in cpu only mode. Openmp is turned on by default
for dbcsr (see dedicated session in this README).

If MKL is present on your system then add this

```shell
   cmake  -DCMAKE_INSTALL_PREFIX=/pyprefix \
      -DCP2K_BUILD_DBCSR=ON \    
      -DCP2K_BLAS_VENDOR=MKL \
      -DCP2K_SCALAPACK_VENDOR=MKL ..
make -j
```

CP2K can be compiled on CRAY systems with the following command line

```shell
   cmake  -DCMAKE_INSTALL_PREFIX=/pyprefix \
 -DCP2K_BUILD_DBCSR=ON \      
-DCP2K_BLAS_VENDOR=SCI \
      -DCP2K_SCALAPACK_VENDOR=SCI ..
make -j
```

MKL (or openblas, etc) can be combined with CUDA or HIP after adding the additional flags enabling
GPU support to the one of the previous cmake commands.

To get CUDA support (V100 gpus) write this command

```shell
cmake -DCMAKE_INSTALL_PREFIX=/pyprefix \
 -DCP2K_BUILD_DBCSR=ON \      
-DCP2K_USE_ACCEL=CUDA \
      -DCP2K_WUTH_GPU=V100 ..
make -j
```

HIP is also supported. To enable it use this command

```shell
cmake -DCMAKE_INSTALL_PREFIX=/pyprefix \
 -DCP2K_BUILD_DBCSR=ON \     
 -DCP2K_USE_ACCEL=HIP \
      -DCP2K_WUTH_GPU=Mi250 ..
make -j
```

The option `CP2K_ENABLE_REGTESTS=ON` will configure the build system such that the regtests are ran
as usuall. `CMake` will create the binaries in the `exe/cmake-build-{cpu,cuda,hip}` directory
located in the cp2k root directory. Users will have to set the environment variabble `CP2K_DATA_DIR`
accordingly.

For additional options please read the section below.

## Use of spack to build cp2k

The \[spack\]{https://github.com/spack/spack.git} tool can build CP2K and all its dependencies with
relative ease. We assume that you have spack available on your system. Enter the command

```shell
spack install cp2k@master build_system=cmake ^openblas
```

To build cp2k master with MPI, openblas, scalapack, fftw and openmp. CP2K with CUDA support is
compiled with

```shell
spack install cp2k@master build_system=cmake +cuda cuda_arch=70 ^openblas ^dbcsr+cuda cuda_arch=70
```

For hip support run this command line

```shell
spack install cp2k@master build_system=cmake +rocm amdgpu_target=gfx1030 ^openblas ^dbcsr+rocm amdgpu_target=gfx1030
```

To build the other depednencies add them to the command line. A full cp2k build can be obtained with

```shell
spack install \
cp2k@master \
build_system=cmake \
+cuda cuda_arch=70 \
+libint \
+libxc \
+spglib \
+sirius \
+elpa \
+pexsi \
+plumed \
+cosma \
^openblas \
^dbcsr+cuda cuda_arch=70 \
^sirius+scalapack \
^cosma+scalapack
```

## Additional dependencies

The CP2K cmake build system also supports the following dependencies :

- `CP2K_USE_SIRIUS = OFF`: add [SIRIUS](https://github.com/electronic-structure/SIRIUS) support to
  CP2K

- `CP2K_USE_FFTW3 = ON`: add support of [fftw3](https://www.fftw.org)

- `CP2K_USE_ELPA = OFF`: add [elpa](https://elpa.mpcdf.mpg.de) support. WARNING: Expect the
  detection to fail at that stage

- `CP2K_USE_PEXSI = OFF`

- `CP2K_USE_SUPERLU = OFF`: detection should work but needs improvement

- `CP2K_USE_COSMA = OFF` : add [cosma](https://github.com/eth-cscs/COSMA) drop-in replacement for
  sclapack pdgemnm

- `CP2K_USE_LIBINT2 = OFF`: add [libint2](https://github.com/evaleev/libint) support (detection
  works ok, module files may not be found at compilation time though)

- `CP2K_USE_VORI = OFF`: detection is fine compilation might fail at linking time (investigating
  why)

- `CP2K_USE_QUIP = OFF` :

- `CP2K_USE_SPGLIB = OFF`: everything alright

- `CP2K_USE_LIBXC = OFF`: Use `pkg-config` by default (ideally the library should be built with
  CMake, if so we can get rid of the `FindLibXC.cmake`). If you installed LIBXC in a non-standard
  location, modify the `PKG_CONFIG_PATH` variable accordingly.

- `CP2K_USE_SPLA = OFF`: enable spla off-loading capabilities (use CMake modules to detect it)

- `CP2K_USE_METIS = OFF`:

- `CP2K_USE_LIBXSMM = OFF`: use [libxsmm](https://libxsmm.readthedocs.io/en/latest/) library for
  small matrices operations. Detection based on `pkg-config`. If you installed libxsmm in a
  non-standard location, modify the `PKG_CONFIG_PATH` variable accordingly.

- `CP2K_USE_ACCEL = NONE, CUDA, HIP`: enable GPU support

- `CP2K_BLAS_VENDOR = auto`: CMake will search for the most common blas / lapack implementations. If
  possible indicate which implementation you are using. Supported values are: `auto` (default),
  `MKL`, `SCI`, `OpenBLAS`, `FlexiBLAS`, `Armpl`.

- `CP2K_SCALAPACK_VENDOR = MKL, SCI, GENERIC`: similar to the previous option but for scalapack

- `CP2K_BLAS_THREADING = sequential, openmp, etc...`: leave the default value (or use it at your own
  peril)

- `CP2K_BLAS_INTERFACE = 32 bits, 64 bits`: size of the integers for the matrices and vectors sizes.
  Default: 32 bits

- `CP2K_DEV_OPTIONS = OFF`: enable developer options. The main purpose is for debugging

  - `CP2K_USE_GRID_GPU = ON`: turn on of gpu support for collocate integrate
  - `CP2K_USE_PW_GPU = ON`: turn on or off gpu fft support
  - `CP2K_USE_DBM_GPU = ON`: turn on or off dbm gpu support
  - `CP2K_CHECK_CONSISTENCY = OFF` : Compare the list of compiled files to the files found in the
    tree. This function is only relevant to developers and should be left to its default value. It
    has no incidence on the binaries generated by the build.

ROCM 5.0.x is known to have a bug in the CMake configuration files. It is possible to go around this
but at the expense of time. The build system was not tested with ROCM 5.1.x but this version shows
performance regression and should be avoided. The Jiting capabilities of ROCM 5.2.x do not work
properly which affects DBCSR. It is highly recommended to update ROCM to the latest version to avoid
all these issues. CP2K can be built with ROCM 5.2.x but GPU support in DBCSR should be turned off
otherwise a crash is expected. ROCM 5.3.x and later seems to work fine.

## Threading with blas and lapack

CP2K expect by default a single threaded version of blas and lapack. The option
`-DCP2K_BLAS_THREADING` can change this behavior. Be careful when tweaking this specific option as
many implementations of blas / lapack are either threaded or (exclusive) sequential but not both. I
think the only exception to this is MKL. Also note that CP2K dependencies will most likely have the
same issue (COSMA with cray-libsci for instance).

## Typical examples of CMake use

The following list gives several examples of CMake command lines. Just add `-DCP2K_USE_SIRIUS=ON` to
add support of SIRIUS in CP2K

```shell
cmake -DCP2K_INSTALL_PREFIX=/myprefix -DCP2K_USE_SIRIUS=ON ..
```

then

```shell
make
```

- MKL

```shell
cmake -DCP2K_INSTALL_PREFIX=/myprefix -DCP2K_BLAS_VENDOR=MKL
-DCP2K_SCALAPACK_VENDOR=MKL ..
```

- Cray environments (with cray-libsci)

```shell
MPICC=cc MPICXX=CC cmake -DCP2K_INSTALL_PREFIX=/myprefix
-DCP2K_BLAS_VENDOR=SCI -DCP2K_SCALAPACK_VENDOR=SCI ..
```

## CUDA / HIP

Let us consider the case where OpenBLAS and netlib scalapack are installed (openmpi or mpich)

```shell
cmake -DCP2K_INSTALL_PREFIX=/myprefix -DCP2K_BLAS_VENDOR=OpenBLAS
-DCP2K_SCALAPACK_VENDOR=GENERIC -DCP2K_USE_ACCEL=CUDA -DCP2K_WITH_GPU=A100 ..
```

If HIP is needed then

```shell
cmake -DCP2K_INSTALL_PREFIX=/myprefix -DCP2K_BLAS_VENDOR=OpenBLAS
-DCP2K_SCALAPACK_VENDOR=GENERIC -DCP2K_USE_ACCEL=HIP -DCP2K_WITH_GPU=Mi250 ..
```

## Troubleshooting

This build system is relatevily stable and was tested on Cray, IBM, and redhat like distributions.
However it is not perfect and problems will show up, that's why the two build systems will be
available. We encourage the user to test the build system just reporting the output of `cmake ..` is
already beneficial.

The best way to report these problems is to open an issue including the CMake command line, error
message, and operating systems.

What is known to fail sometimes

- Nvidia HPC SDK: The location of the cuda maths libraries has changed recently. While CUDA support
  will be detected, the CUDA maths libraries may not.

- HIP : CMake support of ROCM is still under development and is known to fail from time to time.
  Update to ROCM 5.3.x or above to solve the issue.

- BLAS / LAPACK / SCALAPACK : use the options `CP2K_BLAS_VENDOR` and `CP2K_SCALPACK_VENDOR` if you
  know that `MKL` or `SCI` (cray libsci) are present. `-DCP2k_BLAS_VENDOR=OpenBLAS` will also help
  CMake to find OpenBLAS if it is used. Detecting the scalapack library might also fail if the user
  environment is not properly set up.

- BLAS / LAPACK / SCALAPACK: It is possible to set up BLAS / LAPACK / SCALAPACK libraries manually
  with the command

```shell
cmake -DCP2K_BLAS_LINK_LIBRARIES=libmyblas.so -DCP2K_BLAS_VENDOR=CUSTOM 
-DCP2K_LAPACK_LINK_LIBRARIES=libmylapack.so -DCP2K_SCALAPACK_VENDOR=GENERIC 
-DCP2K_SCALAPACK_LIBRARIES=libscalapack.so
```
