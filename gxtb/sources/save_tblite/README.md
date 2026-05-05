# Light-weight tight-binding framework

[![License](https://img.shields.io/github/license/tblite/tblite)](https://github.com/tblite/tblite/blob/HEAD/COPYING.LESSER)
[![Build Status](https://github.com/tblite/tblite/workflows/CI/badge.svg)](https://github.com/tblite/tblite/actions)
[![Doxygen Status](https://github.com/tblite/tblite/workflows/docs/badge.svg)](https://tblite.github.io/tblite)
[![Documentation Status](https://readthedocs.org/projects/tblite/badge/?version=latest)](https://tblite.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/tblite/tblite/branch/main/graph/badge.svg?token=JXIE6myqNH)](https://codecov.io/gh/tblite/tblite)

This project is an effort to create a library implementation of the
extended tight binding (xTB) Hamiltonian which can be shared between
[``xtb``](https://github.com/grimme-lab/xtb) and
[``dftb+``](https://github.com/dftbplus/dftbplus).

Goals of this project are

- create a high-level interface to the extended tight binding methods
- allow low-level access to the components forming the actual energy expression
- provide a framework to handle and manipulate parametrization data

Explicit non-goals are

- provide functionality beyond singlepoint calculations in this library
  (like geometry optimization or molecular dynamics)


## Installation

### Conda package

[![Conda Version](https://img.shields.io/conda/vn/conda-forge/tblite.svg?label=tblite)](https://anaconda.org/conda-forge/tblite)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/tblite-python.svg?label=tblite-python)](https://anaconda.org/conda-forge/tblite-python)

This project is packaged for the *conda* package manager and available on the *conda-forge* channel.
To install the *mamba* package manager we recommend the [mambaforge](https://github.com/conda-forge/miniforge/releases) installer.
If the *conda-forge* channel is not yet enabled, add it to your channels with

```
mamba config --add channels conda-forge
mamba config --set channel_priority strict
```

Once the *conda-forge* channel has been enabled, this project can be installed with:

```
mamba install tblite
```

If you want to enable the Python API as well install

```
mamba install tblite-python
```

It is possible to list all of the versions available on your platform with:

```
mamba repoquery search tblite --channel conda-forge
```

Now you are ready to use ``tblite``.


### FreeBSD Port

[![FreeBSD port](https://repology.org/badge/version-for-repo/freebsd/tblite.svg)](https://www.freshports.org/science/tblite/)

A port for FreeBSD is available

```
pkg install science/tblite
```

In case no package is available build the port using

```
cd /usr/ports/science/tblite
make install clean
```

For more information see the [tblite port details](https://www.freshports.org/science/tblite/).


### Building from source

To build *tblite* from the source code in this repository you need to have
a Fortran compiler supporting Fortran 2008 and one of the supported build systems:

- [meson](https://mesonbuild.com) version 0.57.2 or newer, with
  a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.10 or newer
- [cmake](https://cmake.org) version 3.14 or newer, with
  a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.10 or newer
- [fpm](https://github.com/fortran-lang/fpm) version 0.3.0 or newer

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2008
- [meson](https://mesonbuild.com) version 0.57.2 or newer
- a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.10 or newer
- a LAPACK / BLAS provider, like MKL or OpenBLAS

Meson is the primary build system and provides feature-complete functionality of this project.
CMake and fpm support are available but the functionality of the project is limited.
Currently, *tblite* support GCC 8 and newer or Intel 18 and newer.

Detailed installation instruction are available in the project documentation under the [installation category](https://tblite.readthedocs.io/en/latest/installation.html).


#### Building with meson

Optional dependencies are
- asciidoctor to build the manual page
- C compiler to test the C-API and compile the Python extension module
- Python 3.6 or newer with the CFFI package installed to build the Python API

Setup a build with

```sh
meson setup _build
```

You can select the Fortran compiler by the `FC` environment variable.
To compile and run the projects testsuite use

```sh
meson test -C _build --print-errorlogs
```

To run the more extensive testing for the available parametrizations use

```sh
meson test -C _build --print-errorlogs --benchmark
```

If the testsuites pass you can install with

```sh
meson configure _build --prefix=/path/to/install
meson install -C _build
```

This might require administrator access depending on the chosen install prefix.
For more details see the [meson installation instructions](https://tblite.readthedocs.io/en/latest/installation.html#meson-based-build).


## Usage

This project provides multiple entry points for different usage scenarios.
The simplest way to check out this project is by using the command line driver.


### Command line interface

The ``tblite`` runner executable provides full access to the implemented Hamiltonians, with the [``tblite-run``](man/tblite-run.1.adoc) subcommand.
You can run a single point calculation by providing a geometry input with

```
tblite run --method gfn2 coord
```

To export a parametrization use the [``tblite-param``](man/tblite-param.1.adoc) subcommand

```
tblite param --method gfn2 --output gfn2-xtb.toml
```

The parameter file can be inspected or modified and than used to perform single point calculations with

```
tblite run --param gfn2-xtb.toml coord
```

A preliminary interfaces for the parameter optimization is provided by the [``tblite-fit``](man/tblite-fit.1.adoc) subcommand.
By providing a external command to evaluate the data set in the input file and setting the parameters to relax the fit can be started with

```
tblite fit gfn2-xtb.toml input.toml
```

The provided external program can callback to the main program to evaluate single points or create differences between data outputs using the [``tblite-tagdiff``](man/tblite-tagdiff.1.adoc) subcommand.
By adding the ``--dry-run`` option the setup of the parameter optimization can be inspected and with ``--copy copy.toml`` the input settings can be dumped for user inspection and tweaking.

For more details on all available subcommands checkout the [``tblite(1)``](man/tblite.1.adoc) man page and the respective subcommand man pages.


## Documentation

The user documentation is available at [readthedocs](https://tblite.readthedocs.io).
Additionally, the [doxygen](https://doxygen.nl) generated API documentation is available [here](https://tblite.github.io/tblite).

To build the user documentation locally we use sphinx, install the dependencies you can use the *mamba* package manager

```
mamba create -n sphinx --file doc/requirements.txt
mamba activate sphinx
```

The documentation is build with

```
sphinx-build doc _doc
```

You can inspect the generated documentation by starting a webserver

```
python3 -m http.server -d _doc
```

And open the down URL in a browser.


## Contributing

This is a volunteer open source projects and contributions are always welcome.
Please, take a moment to read the [contributing guidelines](CONTRIBUTING.md) on how to get involved in tblite.


## License

This project is free software: you can redistribute it and/or modify it under
the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This project is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
Lesser GNU General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this project by you, as defined in the
Lesser GNU General Public license, shall be licensed as above, without any
additional terms or conditions.
