# The CP2K Toolchain

## Options

To use the CP2K toolchain installer, you may want to first follow
the instructions given in installer help message:

```shell
> ./install_cp2k_toolchain.sh --help
```

## Basic usage

If you are new to CP2K, and want a basic CP2K binary, then just calling

```shell
> ./install_cp2k_toolchain.sh
```

may be enough. This will use your system gcc, and mpi library (if
existing) and build libint, libxc, fftw and openblas (MKL will be used
instead if MKLROOT env variable is found) from scratch, and give you
a set of arch files that allow you to compile CP2K.

## Complete toolchain build

For a complete toolchain build, with everything installed from
scratch, use:

```shell
> ./install_cp2k_toolchain.sh --install-all
```

### Package settings

One can then change settings for some packages, by setting
`--with-PKG` options after the `--install-all` option. e.g.:

```shell
> ./install_CP2K_toolchain.sh --install-all --with-mkl=system
```

will set the script to look for a system MKL library to link, while
compile other packages from scratch.

### MPI implementation choice

If you do not have a MPI installation, by default the `--install-all`
option will install MPICH for you.  You can change this default
behavior by setting `--mpi-mode` after the `--install-all` option.

## Trouble Shooting

Below are solutions to some of the common problems you may encounter when
running this script.

### The script terminated with an error message

Look at the error message. If it does not indicate the reason for
failure then it is likely that some error occurred during
compilation of the package.  You can look at the compiler log in
the file make.log in the source directory of the package in
`./build`.

One of the causes on some systems may be the fact that too many
parallel make processes were initiated.  By default the script
tries to use all of the processors on you node. You can override
this behavior using `-j` option.

### The script failed at a tarball downloading stage

Try run again with `--no-check-certificate` option. See the help
section for this option for details.

### I've used --with-XYZ=system cannot find the XYZ library

The installation script in "system" mode will try to find a library
in the following system PATHS: `LD_LIBRARY_PATH`, `LD_RUN_PATH`,
`LIBRARY_PATH`, `/usr/local/lib64`, `/usr/local/lib`, `/usr/lib64`,
`/usr/lib`.

For MKL libraries, the installation script will try to look for
MKLROOT environment variable.

You can use:

```shell
> module show XYZ
```

to see exactly what happens when the module XYZ is loaded into your
system. Sometimes a module will define its own PATHS and
environment variables that is not in the default installation
script search path. And as a result the given library will likely
not be found.

The simplest solution is perhaps to find where the root
installation directory of the library or package is, and then use
`--with-XYZ=/some/location/to/XYZ` to tell the script exactly where
to look for the library.

## Licenses

The toolchain only downloads and installs packages that are
[compatible with the GPL](https://www.gnu.org/licenses/gpl-faq.html#WhatDoesCompatMean).
The following table list the licenses of all those packages. While the toolchain
does support linking proprietary software packages, like e.g. MKL, these have to
be installed separately by the user.

<!-- markdownlint-disable MD013 -->
| Package   | License                                                                                 | GPL Compatible |
| --------- | --------------------------------------------------------------------------------------- | -------------- |
| cmake     | [BSD 3-Clause](https://gitlab.kitware.com/cmake/cmake/raw/master/Copyright.txt)         | Yes            |
| cosma     | [BSD 3-Clause](https://github.com/eth-cscs/COSMA/blob/master/LICENCE)                   | Yes            |
| elpa      | [LGPL](https://gitlab.mpcdf.mpg.de/elpa/elpa/blob/master/LICENSE)                       | Yes            |
| fftw      | [GPL](http://www.fftw.org/doc/License-and-Copyright.html)                               | Yes            |
| gcc       | [GPL](https://gcc.gnu.org/git/?p=gcc.git;a=blob_plain;f=COPYING;hb=HEAD)                | Yes            |
| gsl       | [GPL](https://www.gnu.org/software/gsl/doc/html/gpl.html)                               | Yes            |
| hdf5      | [BSD 3-Clause](https://support.hdfgroup.org/ftp/HDF5/releases/COPYING)                  | Yes            |
| libint    | [GPL](https://github.com/evaleev/libint/blob/master/LICENSE)                            | Yes            |
| libsmm    | [GPL](https://github.com/cp2k/cp2k/blob/master/LICENSE)                                 | Yes            |
| libvdwxc  | [GPL](https://gitlab.com/libvdwxc/libvdwxc/blob/master/LICENSE)                         | Yes            |
| libxc     | [MPL](https://gitlab.com/libxc/libxc/blob/master/COPYING)                               | Yes            |
| libxsmm   | [BSD 3-Clause](https://github.com/hfp/libxsmm/blob/master/LICENSE.md)                   | Yes            |
| mpich     | [MPICH](https://github.com/pmodels/mpich/blob/master/COPYRIGHT)                         | [Yes](https://enterprise.dejacode.com/licenses/public/mpich/#license-conditions) |
| openblas  | [BSD 3-Clause](https://github.com/xianyi/OpenBLAS/blob/develop/LICENSE)                 | Yes            |
| openmpi   | [BSD 3-Clause](https://github.com/open-mpi/ompi/blob/master/LICENSE)                    | Yes            |
| pexsi     | [BSD 3-Clause](https://bitbucket.org/berkeleylab/pexsi/src/master/LICENSE)              | Yes            |
| plumed    | [LGPL](https://github.com/plumed/plumed2/blob/master/COPYING.LESSER)                    | Yes            |
| ptscotch  | [CeCILL-C](https://www.labri.fr/perso/pelegrin/scotch/)                                 | [Yes](https://cecill.info/faq.en.html#gpl) |
| quip      | [GPL](https://github.com/libAtoms/QUIP/blob/public/src/libAtoms/COPYRIGHT)              | Yes            |
| fox       | [BSD 3-Clause](https://github.com/andreww/fox/blob/master/LICENSE)                      | Yes            |
| reflapack | [BSD 3-Clause](http://www.netlib.org/lapack/LICENSE.txt)                                | Yes            |
| scalapack | [BSD 3-Clause](http://www.netlib.org/scalapack/LICENSE)                                 | Yes            |
| sirius    | [BSD 2-Clause](https://github.com/electronic-structure/SIRIUS/blob/master/LICENSE)      | Yes            |
| spfft     | [BSD 3-Clause](https://github.com/eth-cscs/SpFFT/blob/master/LICENSE)                   | Yes            |
| spglib    | [BSD 3-Clause](https://github.com/atztogo/spglib/blob/master/COPYING)                   | Yes            |
| superlu   | [BSD 3-Clause](https://github.com/xiaoyeli/superlu/blob/master/License.txt)             | Yes            |
| valgrind  | [GPL](https://sourceware.org/git/?p=valgrind.git;a=blob_plain;f=COPYING;hb=HEAD)        | Yes            |
<!-- markdownlint-enable MD013 -->

## For Developers

### Structure of the toolchain scripts

- `install_cp2k_toolchain.sh` is the main script that will call all
  other scripts.  It contains default flag settings, user input
  parser, calls to each package installation scripts and the
  generator of the CP2K arch files.

- `script/install_*.sh` are the installation scripts for individual
  packages. They are relatively independent, in the sense that by
  running `script/install_PKG.sh` it should install the package on its
  own. However, in practice due to dependencies to other libraries,
  sometimes for a package to be installed this way, it will depend
  on other libraries being already installed and the correct
  environment variables set. At the end of each script, it should
  write to __two__ files: `build/setup_PKG` and `install/setup`.
  - The `build/setup_PKG` file contains all the instructions to set
    the variables used by the `install_cp2k_toolchain.sh` and other
    `script/install_PKG.sh` scripts in order for them to correctly
    compile the toolchain and set the correct library flags for the
    arch files.
  - The `install/setup` file contains all the instructions for setting
    up the correct environment before the user can compile and/or
    run CP2K.

- `script/toolkit.sh` contains all the macros that may be used by all
  of the scripts, and provides functionalities such as prepending a
  path, checking if a library exists etc.

- `script/common_var.sh` contains all of the common variables used by
  each installation scripts. All of the variables in the file should
  have a default value, but allow the environment to set the values,
  using: `VAR=${VAR:-default_value}`.

- `script/parse_if.py` is a python code for parsing the `IF_XYZ(A|B)`
  constructs in the script. Nested structures will be parsed
  correctly. See
  [`IF_XYZ` constructs](./README_FOR_DEVELOPERS.md#the-if_xyz-constructs) below.

- `checksums.sha256` contains the pre-calculated SHA256 checksums for
  the tar balls of all of the packages. This is used by the
  `download_pkg` macro in `script/toolkit.sh`.

- `arch_base.tmpl` contains the template skeleton structure for the
  arch files. The `install_cp2k_toolchain` script will set all the
  variables used in the template file, and then do an eval to expand
  all of `${VARIABLE}` items in `arch_base.tmpl` to give the cp2k arch
  files.

### `enable-FEATURE` options

The `enable-FEATURE` options control whether a FEATURE is enabled or disabled.
Possible values are:

- `yes` (equivalent to using the option-keyword alone)
- `no`

### `with_PKG` and `PKG_MODE` variables

The `with_PKG` options controls how a package is going to be installed:

- either compiled and installed from source downloaded
  (`install`, or the option-keyword alone),
- or linked to locations provided by system search paths (`system`),
- or linked to locations provided by the user (`<path>`, path to some directory),
- or that the installer won't be used (`no`).

For most packages the `with_pkg` variables will act like a switch for
turning on or off the support for this package. However, for
packages serving the same purpose, with the installer needing only
one, an extra variable `PKG_MODE` (e.g. `MPI_MODE`) are used as a
selector.  In this case, while `with_PKG` controls the installation
method, the `PKG_MODE` variable picks which package to actually use.
This provides more flexibility.

### The IF_XYZ constructs

Due to the fact that `install_cp2k_toolchain.sh` needs to produce
several different versions of the arch files: `psmp`, `pdbg`,
`ssmp`, `sdbg`, etc, it will have to resolve different flags for
different arch file versions.

The solution used by this script is to use a syntax construct:

```shell
IF_XYZ(A | B)
```

A parser will then parse this expression to *A* if *XYZ* is passed
to the parser (python `parse_if.py` filename XYZ); and to *B* if *XYZ*
is not passed as command line option (python `parse_if.py` filename).

The `IF_XYZ(A|B)` construct can be nested, so things like:

```shell
IF_XYZ(IF_ABC(flag1|flag2) | flag3)
```

will parse to *flag1* if both *XYZ* and *ABC* are present in the command
line arguments of `parser_if.py`, to *flag2* if only *XYZ* is present,
and *flag3* if nothing is present.

### To ensure portability

- one should always pass compiler flags through the
  `allowed_gcc_flags` and `allowed_gfortran_flags` filters in
  `scripts/toolkit.sh` to omit any flags that are not supported by
  the gcc version used (or installed by this script).

- note that `allowed_gcc_flags` and `allowed_gfortran_flags` do not work
  with `IF_XYZ` constructs. So if you have something like:

```shell
FCFLAGS="IF_XYZ(flag1 flag2 | flag3 flag4)"
```

Then you should break this into:

```shell
XYZ_TRUE_FLAGS="flags1 flags2"
XYZ_FALSE_FLAGS="flags3 flags4"
# do filtering
XYZ_TRUE_FLAGS="$(allowed_gcc_flags $XYZ_TRUE_FLAGS)"
XYZ_FALSE_FLAGS="$(allowed_gcc_flags $XYZ_FALSE_FLAGS)"
```

So that:

```shell
FCFLAGS="IF_XYZ($XYZ_TRUE_FLAGS | $XYZ_FALSE_FLAGS)"
```

- For any intrinsic fortran modules that may be used, it is best to
  check with `check_gfortran_module` macro defined in
  `script/tool_kit.sh`. Depending on the gcc version, some intrinsic
  modules may not exist.

- Try to avoid as much hard coding as possible:
  e.g. instead of setting:

```shell
./configure --prefix=some_dir CC=mpicc FC=mpif90
```

use the common variables:

```shell
./configure --prefix=some_dir CC=${MPICC} FC=${MPIFC}
```

## To keep maintainability it is recommended that we follow these practices

- Reuse as much functionality from the macros defined in the
  `script/toolkit.sh` as possible

- When the existing macros in `script/toolkit.sh` do not provide the
  functionalities you want, it is better to write the new
  functionality as a macro in `script/toolkit.sh`, and then use the
  macro (repeatedly if required) in the actual installation
  script. This keeps the installation scripts uncluttered and more
  readable.

- All packages should install into their own directories, and with a
  lock file created in their respective directory to indicate
  installation has been successful. This allows the script to skip
  over the compilation stages of already installed packages if the
  user terminated the toolchain script at the middle of a run and
  then restarted the script.
