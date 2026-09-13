# Build from Source

CP2K uses the [CMake](https://cmake.org) build system, which detects dependencies and controls the
compilation process. This page describes how to obtain a complete CP2K source tree, prepare the
dependencies and build and install CP2K.

```{note}
Historical releases from 2025 or prior used a custom system based on GNU Make and ARCH files, which
has been deprecated due to portability and maintenance difficulties. There is no longer interest in
providing tech support and documentation for it outside of legacy materials.

The 2026.1 release was the first to exclusively rely on the CMake build system, but the appropriate
[](./build-from-source.md#cmake-configuration-options) had to be assembled manually due to a lack
of automation. The inconvenience has been addressed since the subsequent 2026.2 release, which
introduced the scripts `build_cp2k.sh` (for a [](./build-from-source.md#toolchain-based-build)) and
`make_cp2k.sh` (for a [](./build-from-source.md#spack-based-build-via-make_cp2ksh)).

As such, the following documentation applies to the current CP2K master branch and release versions
`2026.2` and newer.
```

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
It is strongly recommended to use the versioned release tarball rather than GitHub's automatically
generated `Source code` archive, especially for version <=2025.2. The versioned tarball is the
release artifact and contains any source components bundled for that release.
```

### Git checkout

A Git checkout is appropriate for development builds or when a particular branch is required:

```shell
git clone https://github.com/cp2k/cp2k.git cp2k
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

At a minimum, CP2K requires a modern suite of C and Fortran compiler compliant with the C99 and the
F2008 standard respectively. For currently supported compiler versions, see the GitHub Wiki page on
[Compiler Support](https://github.com/cp2k/cp2k/wiki/Compiler-Support).

In addition, CP2K requires [DBCSR](https://github.com/cp2k/dbcsr/), BLAS, and LAPACK; on top of
these, MPI builds require MPI and ScaLAPACK.

Detailed descriptions of available dependencies can be found in the technologies section:

- **[](../technologies/eigensolvers/index)**
- **[](../technologies/accelerators/index)**
- **[](../technologies/libraries)**

Utilities available from package managers `apt-get`, `dnf` and the like, such as `python3`, `pip`,
`git`, `less`, `make`, `sed`, `tar`, `wget`, `zlib`, `unzip`, `bzip2`, and `xz` (this is not an
exhaustive list) are assumed to be readily available infrastructure.

The following two methods provide a CP2K-managed dependency stack. For a manually managed
environment, use the [CMake configuration](#cmake-configuration-options) described at the end.

### Toolchain-based build

The toolchain scripts under the `tools/toolchain` directory build a CP2K-compatible, self-contained
dependency stack and prepare the environment for a subsequent CP2K build. It is operated from the
`install_cp2k_toolchain.sh` script and accompanied by `build_cp2k.sh`. To enter the directory and
read help messages and the complete list of options, run:

```shell
cd ./tools/toolchain/
./install_cp2k_toolchain.sh --help
./build_cp2k.sh --help
```

To configure and install the dependencies requested by default options, run:

```shell
./install_cp2k_toolchain.sh
```

After that, to build and install CP2K linked against them, run:

```shell
./build_cp2k.sh
```

Everything is under the `CP2K_ROOT` directory mentioned above by default: the binaries and libraries
of the dependencies are in `tools/toolchain/install/`, the CP2K build tree in `build/`, and the
headers, modules, binary executables and a `cp2k_env` file in `install/`. Options are also available
to make installed dependencies and program outside of the source tree, and to specify a CMake preset
based on the target architecture.

```shell
./install_cp2k_toolchain.sh --install-dir=/opt/cp2k/toolchain
./build_cp2k.sh --prefix /opt/cp2k --preset native-gnu-x86_64
```

Remember to always source the `cp2k_env` file before starting the program in the current shell and
session. (This is easy to miss on a HPC server where one submits jobs from a login node to another
computing node with job scripts!) For instance, if `CP2K_ROOT` is `/opt/cp2k`, the `--prefix` option
of `build_cp2k.sh` is the default (i.e. `install`), then a quick look at the version information
would use the commands below (as input to the command line prompt or as part of the job script).

```shell
source /opt/cp2k/install/cp2k_env
cp2k.psmp -v
```

The intended output is something like this, with the block of compiler options omitted for brevity.

```text
 CP2K version 2026.2
 Source code revision git:c92cc08
 cp2kflags: omp libint fftw3 libxc elpa parallel scalapack mpi_f08 cosma libxs spglib openblas libdftd4 s_dftd3 mctc-lib tblite libvori libbqb
 compiler: GCC version 14.3.1 20251022 (Red Hat 14.3.1-4)
 compiler target: cpuid   1002 (x86_avx2)
 compiler options:
   [...]
```

#### Pros and Cons

The toolchain-based build is suitable for cases where only the (near-)minimal essential dependencies
are desired. Its reliance and interference with external package managers and internet connection is
also minimal, and as such it comes in more handy in an offline scenario.

The toolchain does not cover every optional dependency or feature combination, such as DLA-Future,
PEXSI, and optional SIRIUS features including NLCG. Its support for GPU-accelerated builds is also
limited. If these features are needed, the Spack-based build as detailed in the next section is the
more recommended method. In fact, due to difficulties in extending the toolchain for increasingly
sophisticated configurations involving nested dependencies, it _may_ be slimmed down or retired in
favor of the modern, well-maintained Spack workflow in the future.

### Spack-based build via `make_cp2k.sh`

`make_cp2k.sh` installs the selected dependency stack with [Spack](https://spack.readthedocs.io) and
then configures, builds, and installs CP2K with CMake. It operates in `CP2K_ROOT`, the root of the
CP2K source tree.

Run the script with its default options to perform a full build including (almost) all features:

```shell
./make_cp2k.sh
```

For a minimal build from scratch, run:

```shell
./make_cp2k.sh -bd -df all
```

Based on which desired features can be added explicitly for building a tailored CP2K binary.

```shell
./make_cp2k.sh -bd -df all -ef libint -ef libxc -ef spglib -ef tblite
```

Use `./make_cp2k.sh --help` to display the complete list of options:

<details>

<summary>Click to see all options (version 2.4)</summary>

```
Usage: make_cp2k.sh [-bd | --build_deps]
                    [-ase ASE_VERSION]
                    [-bd_only | --build_deps_only]
                    [-bp | --build_path PATH]
                    [-bsl | --build_static_libcp2k]
                    [-bt | --build_type (Debug | Release | RelWithDebInfo)]
                    [-cc | --check_conventions]
                    [-cray]
                    [-cv | --cp2k_version (pdbg | psmp | sdbg | ssmp | ssmp-static)]
                    [-df | --disable | --disable_feature (all | FEATURE | PACKAGE | none)]
                    [-ef | --enable | --enable_feature (all | FEATURE | PACKAGE | none)]
                    [-gm | -gpu  | --gpu_model (<CUDA SM code> | P100 | V100 | T400 | A100 | H100 | H200 | GH200 | B200 | none)]
                    [-gromacs GROMACS_VERSION]
                    [-gv | --gcc_version (10 | 11 | 12 | 13 | 14 | 15 | 16)]
                    [-h | --help]
                    [-ip | --install_path PATH]
                    [-j #PROCESSES]
                    [-mpi | --mpi_mode (mpich | no | openmpi)]
                    [-np | --num_packages #PACKAGES]
                    [-opencl]
                    [-preset (native-gnu-x86_64 | native-gnu-arm64 | native-intel | none)]
                    [-rc | --rebuild_cp2k]
                    [-t | --test "TESTOPTS"]
                    [-ta | --test_ase]
                    [-tc | --test_coverage]
                    [-tg | --test_gromacs]
                    [-tp | --test_performance "BENCHMARK_PROFILE"]
                    [-uc | --use_cache (folder | minio | no | none)]
                    [-ue | --use_externals]
                    [-v | --verbose]

Flags:
 -ase                  : Build CP2K with ASE support
 --build_deps          : Force a rebuild of all CP2K dependencies from scratch (removes the spack folder)
 --build_deps_only     : Rebuild ONLY the CP2K dependencies from scratch (removes the spack folder)
 --build_path          : Define the CP2K build path (default: ${CP2K_ROOT})
 --build_static_libcp2k: Build a static CP2K library libcp2k.a instead of the default shared one libcp2k.so
 --build_type          : Set preferred CMake build type for CP2K (default: "Release")
 --check_conventions   : Check compliance with CP2K's coding conventions
 --cp2k_version        : CP2K version to be built (default: "psmp")
 -cray                 : Use Cray specific spack configuration
 --enable_feature      : Enable feature or package (default: all)
 --disable_feature     : Disable feature or package
 -gromacs              : Build GROMACS with CP2K support
 --help                : Print this help information
 --gcc_version         : Use the specified GCC version (default: automatically decided by spack)
 --gpu_model           : Select GPU model (default: none)
 --install_path        : Define the CP2K installation path (default: ./install)
 -j                    : Maximum number of processes used in parallel
 --mpi_mode            : Set preferred MPI mode (default: "mpich")
 --num_packages        : Maximum number of packages built by spack in parallel (default: 4)
 -opencl               : Enable the use of the Open Computing Language (OpenCL)
 -preset               : Use a CMake configure preset, see "cmake --list-presets" (default: native-gnu-x86_64)
 --rebuild_cp2k        : Rebuild CP2K: removes the build folder (default: no)
 --test                : Perform a regression test run after a successful build
 --test_ase            : Build and test CP2K with ASE support
 --test_coverage       : Analyse the code coverage and generate a coverage report
 --test_gromacs        : Build and test GROMACS with CP2K support
 --test_performance    : Perform a benchmark run after a successful build
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
          libgint | libint | libsmeagol | libtorch | libvdwxc | libxs | mimic | openpmd | pexsi | plumed |
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

```{note}
The way Spack resolves a dependency stack is very different from that in a toolchain-based build;
even with the same set of libraries intended to be freshly installed and linked to CP2K, the package
download and disk usage can still make an overall difference.
```

The CP2K built with Spack can be started with the launcher script `install/bin/launch`. Suppose
`CP2K_ROOT` is `/opt/cp2k`, then a version check goes as follows.

```shell
/opt/cp2k/install/bin/launch cp2k.psmp -v
```

#### Testing

To run a regression test immediately after a successful build, add `-t` or `--test` followed by test
options surrounded by double quotes (`"TESTOPTS"`). The test options will be passed to the script
`tests/do_regtest.py` as arguments.

```shell
./make_cp2k.sh --test "--maxtasks 16 --flagslow"
# Alternatively: in case no options are needed, use an empty quote string
./make_cp2k.sh --test ""
```

Alternatively, the script `install/bin/run_tests` produced after a successful build can be used to
start a regression test later. The script prints usage examples at the end of a successful run.

(build-gromacs-cp2k)=

#### GROMACS/CP2K QM/MM

The latest supported GROMACS release (currently v2026.3) for GROMACS/CP2K QM/MM simulations can be
built and tested with

```shell
./make_cp2k.sh -bd --test_gromacs
```

for CP2K versions newer than v2026.2. Other (older) GROMACS versions can be built and tested with

```shell
./make_cp2k.sh -bd -gromacs v2025.2 --test_gromacs
```

A Dockerfile for building GROMACS/CP2K within a container with `podman` is also available. A usage
example is given in the header of that
[Dockerfile](https://raw.githubusercontent.com/cp2k/cp2k/refs/heads/master/tools/docker/Dockerfile.test_spack_gromacs).

## CMake configuration options

Both toolchain and Spack utilize CMake configurations automatically for convenience after preparing
the dependencies. Many of these options allow for CP2K to be built with support of a linked library;
refer to the technologies section for details together with description of available dependencies.

Here are some other important general options you may want to know:

- `-S <Source>` Specifies the path to source tree that contains the `src` directory and the
  `CMakeLists.txt` file. With `CP2K_ROOT` as the working directory, simply use `-S .` for it.
- `-B <Build>` Specifies the path to build. CMake is typically run *out-of-tree* in a separate
  `build/` directory under the `CP2K_ROOT`, corresponding to the `-B build` usage. Note that any
  in-source build or build in any directory with a `CMakeLists.txt` file is strictly *forbidden*.
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
- `-DCP2K_DATA_DIR` Specifies the location of the data for basis sets, pseudopotentials, etc. Any
  data filename specified in the CP2K input without absolute path will be interpreted as under this
  directory, which CP2K prioritizes when attempting to retrieve the data. It is created by copying
  over the existing `data` directory to the location during installation. Default is
  `/path/to/installation/share/cp2k/data`.
- `-DCP2K_ENABLE_CONSISTENCY_CHECKS` Only used for
  [testing](https://dashboard.cp2k.org/archive/misc/index.html).
- `-DCP2K_USE_CRAY_PM_ENERGY` Enables power monitoring on Cray systems.
- `-DCP2K_USE_CRAY_PM_ACCEL_ENERGY` Enables power monitoring of accelerators on Cray systems.
- `-DCP2K_USE_DBCSR_CONFIG` Make dbcsr cmake options (`DBCSR_USE_BLA`) available.

The options above do not cover architecture- and compiler-specific optimization profiles. Instead,
they are handled by [CMake presets](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html)
as defined in `CMakePresets.json`. For instance, the preset `native-gnu-x86_64` applies the flag
`-march=native` to the GNU compilers (namely `gfortran`, `gcc`, `g++`; use `--help=target` on any of
these for more information on the flag) via `CMAKE_<LANG>_FLAGS` variables. It will produce binary
executables that are optimized against the host machine by detecting and utilizing the appropriate
native CPU instruction sets on the target x86_64 architecture. CMake presets are listed with the
command `cmake --list-presets`, where the `native-*` ones are more relevant to everyday use.

### Example

The following example builds CP2K on a x86_64 machine with native optimization in GNU compilers,
also enabling CUDA acceleration for Nvidia A100 GPUs and a few optional dependencies:

```bash
cd <CP2K_REPOSITORY>
mkdir build/
cmake -S . -B build --preset native-gnu-x86_64 \
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
