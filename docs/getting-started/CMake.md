# CMake

CP2K can be compiled with the [CMake] build system. [CMake] is used to detect CP2K dependencies and
configure the compilation process.

CP2K dependencies should be installed independently, either manually or using a package manager. See
\[CP2K's Spack Documentation\] for documentation on how to use the \[Spack\] package manager to
build CP2K's dependencies.

## Using CMake

The most common [CMake] workflow is to create a `build/` directory in the root directory of the
source tree, and run [CMake] as follows:

```bash
cd <CP2K_REPOSITORY>
mkdir build/
cmake -S . -B build             # -S <SOURCE_DIR> -B <BUILD_DIR>
cmake --build build -- -j 32    # -build <BUILD_DIR> -- <BUILD_TOOL_OPTIONS>
```

`<CP2K_REPOSITORY>` is a placeholder for the path to the CP2K repository, containing the source code
to be compiled.

### Generators

By default, [CMake] generates GNU Makefiles (on Linux). With [CMake], it is possible to generate
files for other build systems such as \[Ninja\]:

```bash
cmake -S . -B build -GNinja 
```

### Release and Debug Builds

[CMake] allows to specify a build type. `-DCMAKE_BUILD_TYPE=Release` turns on optimizations, while
`-DCMAKE_BUILD_TYPE=Debug` turns on debug options.

## Requirements

The minimum requirements to use the [CMake] build system are the following:

- [CMake]
- C, C++, and Fortran compilers
- \[DBCSR\]
- OpenMP\]
- BLAS and LAPACK

For an MPI build, the following dependencies are also required:

- MPI
- ScaLAPACK

### DBCSR

\[DBCSR\] is a required dependency of CP2K, and [CMake] expects \[DBCSR\] as an external dependency.

### BLAS, LAPACK, and ScaLAPACK

The major vendors' implementations of BLAS and LAPACK are supported. If CP2K is compiled with MPI
support (`-DCP2K_USE_MPI=ON`), ScaLAPACK is also required.

`-DCP2K_BLAS_VENDOR` and `-DCP2K_SCALAPACK_VENDOR` can be used to control which vendor library to
use.

To build with Intel OneAPI MKL:

```bash
cmake -S . -B build -DCP2K_BLAS_VENDOR=MKL -DCP2K_SCALAPACK_VENDOR=MKL
cmake --build build
```

To build with Cray LibSci:

```bash
cmake -S . -B build -DCP2K_BLAS_VENDOR=SCI -DCP2K_SCALAPACK_VENDOR=SCI
cmake --build build
```

## Build Customization

### Default Dependencies

All optional dependencies are turned off by default, with the exception of MPI.

For simplicity, `-DCP2K_BUILD_OPTIONS` is provided to turn on some of the optional dependencies.
`-DCP2K_BUILD_OPTIONS` can take the following values:

- `CUSTOM` (default)
- `DEFAULT`
- `MINIMAL`
- `FULL`
- `SERIAL`

By default, `-DCP2K_BUILD_OPTIONS=CUSTOM`, meaning that each dependency must be turned on explicitly
with `-DCP2K_USE_<LIBRARY>=ON`. `<LIBRARY>` is the name of the optional dependency (for example:
`-DCP2K_USE_COSMA=ON`, `-DCP2K_USE_LIBXC=ON`, ...).

Please refer to the `CMakeLists.txt` file for an up-to-date list of the dependencies enabled by each
option.

### GPUs

CP2K is GPU-accelerated. In order to enable GPU acceleration with \[CUDA\] or \[HIP\],
`-DCP2K_USE_ACCEL` can be used:

- `-DCP2K_USE_ACCEL=CUDA`: enables \[CUDA\]
- `-DCP2K_USE_ACCEL=HIP`: enables \[HIP\]

The target architecture can be selected with `-DCP2K_WITH_GPU`.

## Example

Build CP2K with CUDA acceleration for Nvidia A100 GPUs, with multiple optional dependencies:

```bash
cd <CP2K_REPOSITORY> && make build/
cmake -S . -B build \
    -GNinja \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_SPGLIB=ON \
    -DCP2K_USE_ELPA=ON \
    -DCP2K_USE_SPLA=ON \
    -DCP2K_USE_SIRIUS=ON \
    -DCP2K_USE_COSMA=ON \
    -DCP2K_USE_ACCEL=CUDA -DCP2K_WITH_GPU=A100

cmake --build build -j 32
```

\[CP2K's Spack Documentation\]: \[Spack\]: https://spack.readthedocs.io/en/latest/ \[DBCSR\]:
https://cp2k.github.io/dbcsr/develop/ \[CUDA\]: https://developer.nvidia.com/cuda-toolkit \[HIP\]:
https://rocm.docs.amd.com/projects/HIP/en/latest/

[cmake]: https://cmake.org
