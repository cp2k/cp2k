# Build from Source

CP2K uses the [CMake](https://cmake.org) build system, which detects dependencies and controls the
compilation process. The dependencies have to be installed in advance, either manually or through a
package manager like [Spack](https://spack.readthedocs.io).

## Dependencies

At a minimum CP2K requirements a modern C and Fortran compiler,
[DBCSR](https://github.com/cp2k/dbcsr/), BLAS, and LAPACK. For parallel builds it also needs at
least MPI and ScaLAPACK. Descriptions of all available build options can be found in the sections on
[](../technologies/eigensolvers/index), [](../technologies/accelerators/index), and
[](../technologies/libraries).

## Example

CMake is typically run *out-of-tree* in a seperate `build/` directory. The following example builds
CP2K with CUDA acceleration for Nvidia A100 GPUs and a few optional dependencies:

```bash
cd <CP2K_REPOSITORY>
mkdir build/
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

## Others Build Options

- `-GNinja` Generates [Ninja](https://ninja-build.org/) build files instead of GNU Makefiles.
- `-DCMAKE_BUILD_TYPE=Debug` Enables debug settings, recommended for development.
- `-DBUILD_SHARED_LIBS=OFF` Disables shared libraries.
- `-DCMAKE_POSITION_INDEPENDENT_CODE=OFF` Disables position-independent code.
- `-DCP2K_USE_EVERYTHING=ON` Enables all dependencies.
- `-DCP2K_ENABLE_CONSISTENCY_CHECKS=ON` Only used for
  [testing](https://dashboard.cp2k.org/archive/misc/index.html).
- `-DCP2K_USE_GRPP=ON` No used anymore.
