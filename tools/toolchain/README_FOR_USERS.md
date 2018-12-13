## Toolchain script usage

#### Options
To use the CP2K toolchain installer, you may want to first follow
the instructions given in installer help message:

```console
> ./install_cp2k_toolchain.sh --help
```

#### Basic usage

If you are new to CP2K, and want a basic CP2K binary, then just calling

```console
> ./install_cp2k_toolchain.sh
```

may be enough. This will use your system gcc, and mpi library (if
existing) and build libint, libxc, fftw and openblas (MKL will be used
instead if MKLROOT env variable is found) from scratch, and give you
a set of arch files that allow you to compile CP2K.

#### Complete toolchain build
For a complete toolchain build, with everything installed from
scratch, use:

```console
> ./install_cp2k_toolchain.sh --install-all
```
##### Package settings

One can then change settings for some packages, by setting
`--with-PKG` options after the `--install-all` option. e.g.:

```console
> ./install_CP2K_toolchain.sh --install-all --with-mkl=system
```

will set the script to look for a system MKL library to link, while
compile other packages from scratch.

##### MPI implementation choice

If you do not have a MPI installation, by default the `--install-all`
option will install MPICH for you.  You can change this default
behavior by setting `--mpi-mode` after the `--install-all` option.

## Trouble Shooting

Below are solutions to some of the common problems you may encounter when running
this script.

##### The script terminated with an error message while installing a particular package.

Look at the error message. If it does not indicate the reason for
failure then it is likely that some error occurred during
compilation of the package.  You can look at the compiler log in
the file make.log in the source directory of the package in
`./build`.

One of the causes on some systems may be the fact that too many
parallel make processes were initiated.  By default the script
tries to use all of the processors on you node. You can override
this behavior using `-j` option.

##### The script failed at a tarball downloading stage

Try run again with `--no-check-certificate` option. See the help
section for this option for details.

##### I have loaded libraries using module load XYZ, but --with-XYZ=system cannot find the XYZ library

The installation script in "system" mode will try to find a library
in the following system PATHS: `LD_LIBRARY_PATH`, `LD_RUN_PATH`,
`LIBRARY_PATH`, `/usr/local/lib64`, `/usr/local/lib`, `/usr/lib64`,
`/usr/lib`.

For MKL libraries, the installation script will try to look for
MKLROOT environment variable.

You can use:

```console
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
