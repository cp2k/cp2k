# Build from Source

CP2K uses the [CMake](https://cmake.org) build system, which detects dependencies and controls the
compilation process. This page describes how to obtain a complete CP2K source tree and build it with
its dependencies.

## Obtaining source code

CP2K is available as a versioned release tarball or through the Git repository. The directory
produced by either method is the root of the CP2K source tree, referred to below as `CP2K_ROOT`.

### Release tarballs

For a stable released version, download the versioned `cp2k-<version>.tar.bz2` asset from the
[CP2K releases](https://github.com/cp2k/cp2k/releases) page and unpack it:

```shell
tar -xjf cp2k-<version>.tar.bz2
cd cp2k-<version>
```

```{tip}
It's strongly recommended to use the versioned release tarball rather than GitHub's automatically
generated `Source code` archive, especially for version <=2025.2. The versioned tarball is the
release artifact and contains any source components bundled for that release.
```

### Git checkout

A Git checkout is appropriate for development builds or when a particular branch is required:

```shell
git clone --recursive https://github.com/cp2k/cp2k.git cp2k
cd cp2k
```

To check out a supported release branch, replace the branch name as appropriate:

```shell
git clone --recursive -b support/v<version> https://github.com/cp2k/cp2k.git cp2k
cd cp2k
```

The `--recursive` is important for versions \<=2025.2, since it includes the DBCSR submodule. If the
repository was cloned without `--recursive`, initialize the required submodule before building:

```shell
git submodule update --init --recursive
```

```{warning}
When configuring CP2K with CMake after building the toolchain inside a CP2K Git checkout with tags,
DBCSR revision detection can use the enclosing CP2K repository and generate incorrect DBCSR version
metadata, causing CMake report an error regarding the DBCSR package discovery; see
[issue #5184](https://github.com/cp2k/cp2k/issues/5184). Therefore, use the versioned release
artifact rather than a Git tag for released CP2K versions.
```

## Setting up dependencies and building CP2K

At a minimum, CP2K requires a modern C and Fortran compiler,
[DBCSR](https://github.com/cp2k/dbcsr/), BLAS, and LAPACK. MPI builds additionally require MPI and
ScaLAPACK. For currently supported compilers, see the
[GitHub Wiki page](https://github.com/cp2k/cp2k/wiki/Compiler-Support).

Detailed descriptions of available dependencies can be found in the technologies section:

- **[](../technologies/eigensolvers/index)**
- **[](../technologies/accelerators/index)**
- **[](../technologies/libraries)**

The following two methods provide a CP2K-managed dependency stack. For a manually managed
environment, use the CMake configuration described below.

### Toolchain-based build

The toolchain scripts under `tools/toolchain` build a CP2K-compatible dependency stack and prepare
the environment for a subsequent CP2K build. Use `./install_cp2k_toolchain.sh --help` to display the
help message and the complete list of options.

To build CP2K with the toolchain, run `install_cp2k_toolchain.sh` from `tools/toolchain` with the
desired toolchain options to configure and install the requested dependencies, then use
`build_cp2k.sh` to build and install CP2K accordingly.

```{note}
The toolchain does not cover every optional dependency or feature combination, such as DLA-Future,
PEXSI, and optional SIRIUS features including NLCG. If these features are needed, it's recommended
to choose the Spack-based method.
```

### Spack-based build via `make_cp2k.sh`

`make_cp2k.sh` installs the selected dependency stack with [Spack](https://spack.readthedocs.io) and
then configures, builds, and installs CP2K with
[CMake](https://cmake.org/cmake/help/latest/index.html). It operates in `CP2K_ROOT`, the root of the
CP2K source tree.

```{note}
This build path is only available on the current CP2K master branch and CP2K release versions
`2026.2` and newer.
```

Run the script with its default options:

```console
./make_cp2k.sh
```

Use `./make_cp2k.sh --help` to display the complete list of options:

<details>

<summary>Click to see all options (version 1.9)</summary>

```
Usage: make_cp2k.sh [-bd | --build_deps]
                    [-bd_only | --build_deps_only]
                    [-bp | --build_path PATH]
                    [-bsl | --build_static_libcp2k]
                    [-bt | --build_type (Debug | Release | RelWithDebInfo)]
                    [-cray]
                    [-cv | --cp2k_version (pdbg | psmp | sdbg | ssmp | ssmp-static)]
                    [-df | --disable | --disable_feature (all | FEATURE | PACKAGE | none)
                    [-ef | --enable | --enable_feature (all | FEATURE | PACKAGE | none)
                    [-gm | -gpu  | --gpu_model (<CUDA SM code> | P100 | V100 | T400 | A100 | H100 | H200 | GH200 | none)]
                    [-gv | --gcc_version (10 | 11 | 12 | 13 | 14 | 15 | 16)]
                    [-h | --help]
                    [-ip | --install_path PATH]
                    [-j #PROCESSES]
                    [-mpi | --mpi_mode (mpich | no | openmpi)]
                    [-np | --num_packages #PACKAGES]
                    [-rc | --rebuild_cp2k]
                    [-t | --test "TESTOPTS"]
                    [-uc | --use_cache (folder | minio | no | none)]
                    [-ue | --use_externals]
                    [-v | --verbose]

Flags:
 --build_deps          : Force a rebuild of all CP2K dependencies from scratch (removes the spack folder)
 --build_deps_only     : Rebuild ONLY the CP2K dependencies from scratch (removes the spack folder)
 --build_path          : Define the CP2K build path (default: ${CP2K_ROOT})
 --build_static_libcp2k: Build a static CP2K library libcp2k.a instead of the default shared one libcp2k.so
 --build_type          : Set preferred CMake build type for CP2K (default: "Release")
 --cp2k_version        : CP2K version to be built (default: "psmp")
 -cray                 : Use Cray specific spack configuration
 --enable_feature      : Enable feature or package (default: all)
 --disable_feature     : Disable feature or package
 --help                : Print this help information
 --gcc_version         : Use the specified GCC version (default: automatically decided by spack)
 --gpu_model           : Select GPU model (default: none)
 --install_path        : Define the CP2K installation path (default: ./install)
 -j                    : Maximum number of processes used in parallel
 --mpi_mode            : Set preferred MPI mode (default: "mpich")
 --num_packages        : Maximum number of packages built by spack in parallel (default: 4)
 --rebuild_cp2k        : Rebuild CP2K: removes the build folder (default: no)
 --test                : Perform a regression test run after a successful build
 --use_cache           : Use a "folder", a "MinIO" object storage container (requires podman) or "no" cache
                         Set the environment variable SPACK_CACHE to specify the folder name, e.g.
                         SPACK_CACHE="file://${CP2K_ROOT}/spack_cache" (default)
 --use_externals       : Use external packages installed on the host system. This results in much
                         faster build times, but it can also cause conflicts with outdated packages
                         pulled in from the host system, e.g. old python or gcc versions
 --verbose             : Write verbose output

Hints:
 - Remove the folder ${CP2K_ROOT}/build to (re)build CP2K from scratch
   (see also --rebuild_cp2k flag)
 - Remove the folder ${CP2K_ROOT}/spack to (re)build CP2K and all its dependencies from scratch
   (see also --build_deps flag)
 - The folder ${CP2K_ROOT}/install is updated after each successful run

Packages: all | ace | cosma | deepmd | dftd4 | dlaf | elpa | fftw3 | gauxc | greenx | hdf5 | libfci |
          libint | libsmeagol | libtorch | libvdwxc | libxs | mimic | openpmd | pexsi | plumed |
          sirius | spfft | spglib | spla | tblite | trexio | vori

Features: cray_pm_accel_energy | cusolver_mp | dbm_gpu | elpa_gpu | grid_gpu | pw_gpu |
          spla_gemm_offloading | unified_memory
```

</details><br>

`make_cp2k.sh` creates and reuses the following directories below `CP2K_ROOT`:

- `spack/` contains the local Spack installation and dependency stack. Remove or rename it to
  rebuild all dependencies from scratch; `--build_deps` and `--build_deps_only` provide the same
  behavior from the script.
- `build/` contains the CMake build tree. Remove or rename it, or use `--rebuild_cp2k`, to
  reconfigure and rebuild CP2K from scratch.
- `install/` contains the installed CP2K files and is updated after each successful build.

By default, compiled packages are also stored in a local cache. This significantly accelerates later
dependency builds; see `--use_cache` for the available cache backends.

#### Testing

Add `-t` or `--test` followed by `TESTOPTS` to run a regression test after a successful build:

```console
./make_cp2k.sh --test "--maxtasks 16 --flagslow"
```

Alternatively, run `install/bin/run_tests` after a successful build. The script prints usage
examples at the end of a successful run.

## CMake configuration options

Detailed descriptions of most build options can be found in the technologies section together with
description of available dependencies.

Here are some other important general options you may want to know:

- `-G <Generator>` Specifies which type of build files would be generated. Default is
  `Unix Makefiles`, which generates a GNU Makefile and allows you to build with running `make` in
  the build directory. For GPU-accelerated builds, it is strongly advised to use `Ninja` as
  generator, which is also used by `make_cp2k.sh`; in this case, please ensure that Ninja is
  installed on your host system.
- `-DCMAKE_BUILD_TYPE` Valid vaules are `Release` (default) and `Debug` (enables debug settings and
  generates `pdbg` or `sdbg` instead of `psmp` or `ssmp`; recommended for development).
- `-DCMAKE_INSTALL_PREFIX` Specifies the installation path of CP2K. Assuming it is set to
  `/path/to/installation`, there will be several subdirectories: `bin` for binaries like
  `cp2k.psmp`, `include` for module files and headers, `lib` or `lib64` for libraries, and `share`
  for some other files such as basis data. Default is `/usr/local`.
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
    -DCP2K_USE_MPI=ON \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_SPGLIB=ON \
    -DCP2K_USE_ELPA=ON \
    -DCP2K_USE_SPLA=ON \
    -DCP2K_USE_SIRIUS=ON \
    -DCP2K_USE_COSMA=ON \
    -DCP2K_USE_ACCEL=CUDA \
    -DCP2K_WITH_GPU=A100

cmake --build build -j 32
cmake --install build
```

- The commands `cmake --build build -j 32` and `cmake --install build` can be replaced by a single
  command `cmake --build build --target install -j 32`.
- If you want to clean your build cache after installing in order to save space, simply run
  `cmake --build build --target clean`.
