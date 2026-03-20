# Build from Source

CP2K uses the [CMake](https://cmake.org) build system, which detects dependencies and controls the
compilation process. The dependencies have to be installed in advance.

Currently, CP2K offers two convenient ways for building CP2K from source with its dependencies.
One is to use the classic toolchain scripts to customize and install the required dependencies in a
single step, thereby preparing the environment for compiling CP2K.

A more modern and automated approach is via the `make_cp2k.sh` script, which
leverages [Spack](https://spack.readthedocs.io) to install the dependencies and
subsequently build CP2K. The latter approach is only available for the current
[CP2K master branch](https://github.com/cp2k/cp2k) and CP2K release versions `2026.02`
and newer.

## make_cp2k.sh

The `make_cp2k.sh` script supports and facilitates the build of CP2K and its dependencies from scratch
using [Spack](https://spack.readthedocs.io) and [CMake](https://cmake.org/cmake/help/latest/index.html).
The build is performed locally within the folder `CP2K_ROOT` which defaults to the current working directory.
This should be a `cp2k` folder containing the CP2K source tree which can be download with

```
> git clone https://github.com/cp2k/cp2k.git cp2k
```

and after

```  
> cd cp2k
```

this script can be run using the default options with

```
> ./make_cp2k.sh
```

The flags `-h` or `--help` print the available options.
<details>

  <summary>Click to see all options (of version 1.6)</summary>

```
Usage: make_cp2k.sh [-bd | --build_deps]
                    [-bd_only | --build_deps_only]
                    [-bt | --build_type (Debug | Release | RelWithDebInfo)]
                    [-cray]
                    [-cv | --cp2k_version (pdbg | psmp | sdbg | ssmp | ssmp-static)]
                    [-df | --disable | --disable_feature (all | FEATURE | PACKAGE | none)
                    [-dlc | --disable_local_cache]
                    [-ef | --enable | --enable_feature (all | FEATURE | PACKAGE | none)
                    [-gm | -gpu  | --gpu_model (<CUDA SM code> | P100 | V100 | T400 | A100 | H100 | H200 | GH200 | none)]
                    [-gv | --gcc_version (10 | 11 | 12 | 13 | 14 | 15 | 16)]
                    [-h | --help]
                    [-ip | --install_path PATH]
                    [-j #PROCESSES]
                    [-mpi | --mpi_mode (mpich | no | openmpi)]
                    [-np | --num_packages #PACKAGES]
                    [-rc | --rebuild_cp2k]
                    [-t | -test "TESTOPTS"]
                    [-ue | --use_externals]
                    [-v | --verbose]

Flags:
 --build_deps         : Force a rebuild of all CP2K dependencies from scratch (removes the spack folder)
 --build_deps_only    : Rebuild ONLY the CP2K dependencies from scratch (removes the spack folder)
 --build_type         : Set preferred CMake build type for CP2K (default: Release)
 --cp2k_version       : CP2K version to be built (default: psmp)
 -cray                : Use Cray specific spack configuration
 --disable_local_cache: Don't add local spack cache
 --enable_feature     : Enable feature or package (default: all)
 --disable_feature    : Disable feature or package
 --help               : Print this help information
 --gcc_version        : Use the specified GCC version (default: automatically decided by spack)
 --gpu_model          : Select GPU model (default: none)
 --install_path       : Define the CP2K installation path (default: ./install)
 -j                   : Maximum number of processes used in parallel
 --mpi_mode           : Set preferred MPI mode (default: mpich)
 --num_packages       : Maximum number of packages built by spack in parallel (default: 4)
 --rebuild_cp2k       : Rebuild CP2K: removes the build folder (default: no)
 --test               : Perform a regression test run after a successful build
 --use_externals      : Use external packages installed on the host system. This results in much
                        faster build times, but it can also cause conflicts with outdated packages
                        pulled in from the host system, e.g. old python or gcc versions
 --verbose            : Write verbose output

Hints:
 - Remove the folder ${CP2K_ROOT}/build to (re)build CP2K from scratch
   (see also --rebuild_cp2k flag)
 - Remove the folder ${CP2K_ROOT}/spack to (re)build CP2K and all its dependencies from scratch
   (see also --build_deps flag)
 - The folder ${CP2K_ROOT}/install is updated after each successful run

Packages: all | ace | cosma | deepmd | dftd4 | dlaf | elpa | fftw3 | greenx | hdf5 | libint2 |
          libsmeagol | libtorch | libvdwxc | libxsmm | mimic | openpmd | pexsi | plumed | sirius |
          spfft | spglib | spla | tblite | trexio | vori

Features: cray_pm_accel_energy | cusolver_mp | dbm_gpu | elpa_gpu | grid_gpu | pw_gpu |
          spla_gemm_offloading | unified_memory
```

</details>

The first run will take longer as it will build all CP2K dependencies with Spack. The Spack installation
is kept fully local in the subfolder `cp2k/spack` which corresponds to the `tools/toolchain/install` folder
created by the CP2K toolchain.

Subsequent runs of the `make_cp2k.sh` script will use the software stack from that `cp2k/spack` folder.
A rebuild of all CP2K dependencies can be enforced simply by removing or renaming the folder `cp2k/spack`.
The latter allows for keeping different software stacks (see also `-bd` and `-bd_only` flags).

It is recommended to install podman to take advantage of a local cache. This will accelerate the
(re)build of the CP2K dependencies with Spack significantly.

After the CP2K dependencies are built with Spack, CP2K itself is built and installed using `CMake`
in the subfolders `cp2k/build` and `cp2k/install`, respectively.

Subsequent runs of the script will use the CMake configuration in the subfolder `cp2k/build`.
A rebuild of CP2K from scratch can be enforced by removing or renaming that subfolder.
  
A CP2K regression run can be launched automatically by adding the flag `-t ""` (or `--test ""`).
This flag expects a string with the `TESTOPTS`, e.g.
```
> ./make_cp2k.sh -t "--maxtasks 8 --restrictdir QS/regtest-gpw-1"
```
Alternatively, the script `cp2k/install/bin/run_tests` can be launched after a successful CP2K build.
Usage examples are printed at the end of a (successful) `make_cp2k.sh` run.

### Docker containers

Check the folder [`cp2k/tools/docker`](https://github.com/cp2k/cp2k/tree/master/tools/docker) for Docker files
with the name pattern `Dockerfile.test_spack_*` like `Dockerfile.test_spack_psmp` to build Docker containers
using `make_cp2k.sh`. Each Docker file provides a usage example in its header
which employs [`podman`](https://docs.podman.io/en/latest/) for building the the CP2K docker containers.

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
