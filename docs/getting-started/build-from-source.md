# Build from Source

CP2K uses the [CMake](https://cmake.org) build system, which detects dependencies and controls the
compilation process. The dependencies have to be installed in advance.

Currently, CP2K offers three convenient methods for building CP2K from source with its dependencies.
One is to use the classic toolchain scripts to customize and install the required dependencies in a
single step, thereby preparing the environment for compiling CP2K. A more modern and automated
approach is via `make_cp2k.sh`, which leverages [Spack](https://spack.readthedocs.io) to install the
dependencies and subsequently build CP2K. Alternatively, one can use Spack directly to install all
the dependencies and setup the build environment, as described in
[](./build-with-spack.md#developer-workflow).

## Dependencies and build options

At a minimum CP2K requirements a modern C and Fortran compiler,
[DBCSR](https://github.com/cp2k/dbcsr/), BLAS, and LAPACK. For parallel builds it also needs at
least MPI and ScaLAPACK.

For currently supported compilers, see [here](https://www.cp2k.org/dev:compiler_support).

Detailed descriptions of most build options can be found in the technologies section:

- **[](../technologies/eigensolvers/index)**
- **[](../technologies/accelerators/index)**
- **[](../technologies/libraries)**

There are some other important general options you may want to know:

- `-G <Generator>` Specifies which type of build files would be generated. Default is
  `Unix Makefiles`, which generates a GNU Makefile and allows you to build with running `make` in
  the build directory. For GPU-accelerated builds, it is strongly advised to use `Ninja` as
  generator, which is also used by `make_cp2k.sh`; in this case, please ensure that Ninja is
  installed on your host system.
- `-DCMAKE_BUILD_TYPE` Valid vaules are `Release` (default) and `Debug` (enables debug settings and
  generates `pdbg` or `sdbg` instead of `psmp` or `ssmp`; recommended for development).
- `-DCMAKE_INSTALL_PREFIX` Specifies the installation path of CP2K. Assuming it is set to
  `/path/to/installation`, there will be several subdirectories: `bin` for binaries like
  `cp2k.psmp`, `include` for module files and headers, `lib/lib64` for libraries, and `share` for
  some other files such as basis data. Default is `/usr/local`.
- `-DBUILD_SHARED_LIBS` Specifies if shared libraries are built. Default is `ON`; if set `OFF`, a
  static library will be built instead.
- `-DCMAKE_POSITION_INDEPENDENT_CODE` Specifies if position-independent code is enabled.

Along with some options with CP2K:

- `-DCP2K_USE_EVERYTHING` Enables all dependencies or not.
- `-DCP2K_DATA_DIR` Specifies the location of the data of basis and potentials. Default is
  `/path/to/installation/share/cp2k/data`.
- `-DCP2K_ENABLE_CONSISTENCY_CHECKS` Only used for
  [testing](https://dashboard.cp2k.org/archive/misc/index.html).
- `-DCP2K_USE_CRAY_PM_ENERGY` Enables power monitoring on Cray systems.
- `-DCP2K_USE_CRAY_PM_ACCEL_ENERGY` Enables power monitoring of accelerators on Cray systems.
- `-DCP2K_USE_DBCSR_CONFIG` Make dbcsr cmake options (`DBCSR_USE_BLA`) available.

Note that CMake is typically run *out-of-tree* in a seperate `build/` directory. We don't allow
in-source builds; if you run CMake in the root directory, it will give error.

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

- The commands `cmake --build build -j 32` and `cmake --install build` can be replaced by a single
  command `cmake --build build --target install -j 32`

## Cleaning build cache

If you want to clean your build cache after installing in order to save space, simply run:

```
cmake --build build --target clean
```
