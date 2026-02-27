# Build from Source

CP2K uses the [CMake](https://cmake.org) build system, which detects dependencies and controls the
compilation process. The dependencies have to be installed in advance.

Currently, CP2K offers two convenient methods for installing CP2K and its dependencies. One is to
use the classic toolchain scripts to customize and install the required dependencies in a single
step, thereby preparing the environment for compiling CP2K. A more modern and automated approach is
via `make_cp2k.sh`, which leverages [Spack](https://spack.readthedocs.io) to install the
dependencies and subsequently build CP2K.

## Dependencies and build options

At a minimum CP2K requirements a modern C and Fortran compiler,
[DBCSR](https://github.com/cp2k/dbcsr/), BLAS, and LAPACK. For parallel builds it also needs at
least MPI and ScaLAPACK.

GCC is the most tested compiler from our side, while most of other compilers should work. It is not
recommended to use Intel oneAPI for now; although it could work (up to version 2025.2.1), there are
still some problems remaining unresolved.

Detailed descriptions of most build options can be found in the technologies section:

- **[](../technologies/eigensolvers/index)**
- **[](../technologies/accelerators/index)**
- **[](../technologies/libraries)**

There are some other important general options you may want to know:

- `-G <Generator>` Specifies which type of build files would be generated. Default is
  `Unix Makefiles`, which generates a GNU Makefile and allow you to build with running `make` in the
  build directory. For GPU-accelerated builds, it is strongly advised to use `Ninja` as generator,
  which is also used by `make_cp2k.sh`; in this case, please ensure that Ninja is installed on your
  host system.
- `-DCMAKE_INSTALL_PREFIX` Specifies the installation path of CP2K. Assuming it is set to
  `/path/to/installation`, there will be several subdirectories: `bin` for binaries like
  `cp2k.psmp`, `include` for module files and headers, `lib/lib64` for libraries, and `share` for
  some other files such as basis data. Default is `/usr/local`.
- `-DCP2K_DATA_DIR` Specifies the location of the data of basis and potentials. Default is
  `/path/to/installation/share/cp2k/data`.
- `-DBUILD_SHARED_LIBS` Specifies if shared libraries are built. Default is `ON`; if set `OFF`, a
  static library will be built instead.

Note that CMake is typically run *out-of-tree* in a seperate `build/` directory. We don't allow
in-source build; if you run CMake with `-B .` (default) in the root directory with `CMakeList.txt`,
it will give error.

### Example

The following example builds CP2K with CUDA acceleration for Nvidia A100 GPUs and a few optional
dependencies :

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
cmake --install build
```

- `-GNinja` Generates [Ninja](https://ninja-build.org/) build files instead of GNU Makefiles.
- `-DCMAKE_BUILD_TYPE=Debug` Enables debug settings, recommended for development.
- `-DBUILD_SHARED_LIBS=OFF` Disables shared libraries.
- `-DCMAKE_POSITION_INDEPENDENT_CODE=OFF` Disables position-independent code.
- `-DCP2K_USE_EVERYTHING=ON` Enables all dependencies.
- `-DCP2K_ENABLE_CONSISTENCY_CHECKS=ON` Only used for
  [testing](https://dashboard.cp2k.org/archive/misc/index.html).
- `-DCP2K_USE_CRAY_PM_ENERGY` Enables power monitoring on Cray systems.
- `-DCP2K_USE_CRAY_PM_ACCEL_ENERGY` Enables power monitoring of accelerators on Cray systems.
- `-DCP2K_USE_DBCSR_CONFIG` Make dbcsr cmake options (`DBCSR_USE_BLA`) available.
- The commands `cmake --build build -j 32` and `cmake --install build` can be replaced by a single
  command `cmake --build build --target install -j 32`

## Cleaning build cache

If you want to clean your build cache after installing in order to save space, simply run:

```
cmake --build build --target clean
```
