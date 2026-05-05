# The D3 dispersion model

[![Latest Version](https://img.shields.io/github/v/release/dftd3/simple-dftd3)](https://github.com/dftd3/simple-dftd3/releases/latest)
[![LGPL-3.0-or-later](https://img.shields.io/github/license/dftd3/simple-dftd3)](COPYING)
[![JOSS](https://joss.theoj.org/papers/1a0f4b4571b8a362d596bd5759572d7f/status.svg)](https://joss.theoj.org/papers/1a0f4b4571b8a362d596bd5759572d7f)
[![CI](https://github.com/dftd3/simple-dftd3/workflows/CI/badge.svg)](https://github.com/dftd3/simple-dftd3/actions)
[![Documentation](https://readthedocs.org/projects/dftd3/badge/?version=latest)](https://dftd3.readthedocs.io/en/latest/)
[![docs](https://github.com/dftd3/simple-dftd3/actions/workflows/docs.yml/badge.svg)](https://dftd3.github.io/simple-dftd3/)
[![codecov](https://codecov.io/gh/dftd3/simple-dftd3/branch/main/graph/badge.svg)](https://codecov.io/gh/dftd3/simple-dftd3)

This package provides a library first implementation of the DFT-D3 dispersion correction
(see [*JCP* **132**, 154104 (2010)](https://dx.doi.org/10.1063/1.3382344)
and [*JCC* **32**, 1456 (2011)](https://dx.doi.org/10.1002/jcc.21759) for details).
Usable via the command line interface of the [`s-dftd3` executable](man/s-dftd3.1.adoc),
in Fortran via the [`dftd3` module](https://dftd3.readthedocs.io/en/latest/api/fortran.html),
with the C interface by including the [`dftd3.h` header](https://dftd3.readthedocs.io/en/latest/api/c.html),
or in Python by using the [`dftd3` module](https://dftd3.readthedocs.io/en/latest/api/python.html).
Additionally, the geometric counter-poise correction (see [*JCP* **136**, 154101 (2012)](https://dx.doi.org/10.1063/1.3700154))
is available to correct for basis set superposition errors and construct composite electronic structure methods of the 3c family.


## Installation

A full guide for installing this package or building from source can be found in the
[installation guide](https://dftd3.readthedocs.io/en/latest/installation.html).

### Conda package

[![Conda Version](https://img.shields.io/conda/vn/conda-forge/simple-dftd3.svg?label=simple-dftd3)](https://anaconda.org/conda-forge/simple-dftd3)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/dftd3-python.svg?label=dftd3-python)](https://anaconda.org/conda-forge/dftd3-python)

The preferred way of installing this package is from the *conda-forge* distribution via the *mamba* package manager.
To install the *mamba* package manager the [miniforge](https://github.com/conda-forge/miniforge/releases) installer can be used.
If you already have the *mamba* package manager installed add the *conda-forge* channel if it is not yet enabled:

```
mamba config --add channels conda-forge
```

Once the *conda-forge* channel has been enabled, this project can be installed with:

```
mamba install simple-dftd3
```

If you want to enable the Python API as well install

```
mamba install dftd3-python
```

Now you are ready to use ``s-dftd3``.


### Building from Source

To build this project from the source code in this repository you need a Fortran compiler supporting
at least the Fortran 2008 standard.
If a C compiler is available the C interface tests will be build with this compiler.
For building the Python bindings a C compiler is required to build the Python extension module.

For building from source the

- [meson](https://mesonbuild.com) version 0.58 or newer, with
  a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer
- [asciidoctor](https://asciidoctor.org/) to build the manual page (optional)
- [FORD](https://forddocs.readthedocs.io/en/stable/) to build the developer documentation (optional)
- Python 3.7 or newer with the [CFFI](https://cffi.readthedocs.io/en/latest/), [NumPy](https://numpy.org/) and [setuptools](https://setuptools.pypa.io/en/latest/index.html) package installed to build the Python API

This project builds on many existing Fortran packages, the following ones are used as dependencies

- [mctc-lib](https://github.com/grimme-lab/mctc-lib):
  For reading geometry files in a wide range of formats
  (like xyz, PBD, mol, SDF, Turbomole, DFTB+, Vasp, FHI-aims, Gaussian, QChem, QCSchema JSON, Chemical JSON)
- [toml-f](https://toml-f.readthedocs.org)
  For reading the parameter file with all damping parameters
  (see [`assets/`](assets/parameters.toml))
- [mstore](https://github.com/grimme-lab/mstore)
  For molecule fixtures to use when defining unit tests
  (testing only)
- [test-drive](https://github.com/fortran-lang/test-drive)
  Unit testing framework
  (testing only)

If these Fortran packages are not available, the build system will automatically fetch them and build
them as part of this project.

To setup a build use

```
meson setup _build --prefix=/path/to/install
```

You can select the Fortran compiler by the `FC` environment variable, similarly the C compiler is selected from the `CC` environment variable.
To compile and run the projects testsuite use

```
meson test -C _build --print-errorlogs
```

If the testsuite passes you can install with

```
meson install -C _build
```

This might require administrator access depending on the chosen install prefix.

Meson is the primary build system and provides feature-complete functionality of this project.
If you for any reason cannot use meson, this project also supports CMake and fpm as build systems.
For more details checkout the [installation guide](https://dftd3.readthedocs.io/en/latest/installation.html#building-from-source).


## Generating docs

The documentation is generated using [sphinx](https://www.sphinx-doc.org) for the general documentation and the Python API,
[FORD](https://forddocs.readthedocs.io) for the Fortran API documentation,
and [asciidoctor](https://asciidoctor.org/) for the command line interface manpage.


### Documentation pages (Sphinx)

[![Documentation](https://readthedocs.org/projects/dftd3/badge/?version=latest)](https://dftd3.readthedocs.io/en/latest/)

For generating the main documentation pages, install the documentation dependencies with

```
pip install -r doc/requirements.txt
```

The pages can be built with

```
sphinx-build doc _docs -b html
```

To view the final pages you can start a HTTP server via

```
python -m http.server -d _docs
```

And open the shown URL in a browser (usually this is https://localhost:8000).
The documentation is automatically deployed from the main branch and can be viewed on [readthedocs](https://dftd3.readthedocs.io/).


### API documentation (Ford)

[![docs](https://github.com/dftd3/simple-dftd3/actions/workflows/docs.yml/badge.svg)](https://dftd3.github.io/simple-dftd3/)

To generate the API documentation of the Fortran library install ford via

```
mamba install ford
```

The API documentation can be generated with

```
ford docs.md -o _api
```

To view the final pages you can start a HTTP server via

```
python -m http.server -d _docs
```

And open the shown URL in a browser.
The API documentation is automatically deployed from the main branch and can be viewed on [GitHub pages](https://dftd3.github.io/simple-dftd3).


## Usage

DFT-D3 calculations can be performed with the ``s-dftd3`` executable.
To calculate the dispersion correction for PBE0-D3(BJ)-ATM run:

```
s-dftd3 --bj pbe0 --atm coord
```

In case you want to access the DFT-D3 results from other programs, dump the results to JSON with
(the ``--noedisp`` flag prevents the ``.EDISP`` file generation):

```
s-dftd3 --bj pbe0 --atm --json --grad --noedisp struct.xyz
```

Dispersion related properties can be calculated as well:

```
s-dftd3 --property geo.gen
```

For an overview over all command line arguments use the ``--help`` argument or checkout the [``s-dftd3(1)``](man/s-dftd3.1.adoc) manpage.


## API access

This DFT-D3 implementation provides first class API support Fortran, C and Python.
Other programming languages should try to interface via one of those three APIs.
To provide first class API support for a new language the interface specification should be available as meson build files.


### Fortran API

The recommended way to access the Fortran module API is by using ``dftd3`` as a meson subproject.
Alternatively, the project is accessible by the Fortran package manager ([fpm](https://github.com/fortran-lang/fpm)) or as CMake subproject as explained above.

The complete API is available from ``dftd3`` module, the individual modules are available to the user as well but are not part of the public API and therefore not guaranteed to remain stable.
API compatibility is only guaranteed for the same minor version, while ABI compatibility cannot be guaranteed in a pre 1.0 stage.

The communication with the Fortran API uses the ``error_type`` and ``structure_type`` of the modular computation tool chain library (mctc-lib) to handle errors and represent geometries, respectively.


### C API

The C API provides access to the basic Fortran objects and their most important methods to interact with them.
All Fortran objects are available as opaque ``void*`` in C and can only be manipulated with the correct API calls.
To evaluate a dispersion correction in C four objects are available:

1. the error handler:

   Simple error handler to carry runtime exceptions created by the library.
   Exceptions can be handled and/or transfered to the downstream error handling system by this means.

2. the molecular structure data:

   Provides a representation of the molecular structure with immutable number of atoms, atomic species, total charge and boundary conditions.
   The object provides a way to update coordinates and lattice parameters, to update immutable quantities the object has to be recreated.

3. the dispersion model:

   Instantiated for a given molecular structure type, it carries no information on the geometry but relies on the atomic species of the structure object.
   Recreating a structure object requires to recreate the dispersion model as well.

4. the damping parameters:

   Damping parameter object determining the short-range behaviour of the dispersion correction.
   Standard damping parameters like the rational damping are independent of the molecular structure and can easily be reused for several structures or easily exchanged.

5. the counter-poise parameters:

   Counter-poise parameter object determining the basis set specific correction for basis set superposition error.
   Recreating a structure object requires to recreate the counter-poise parameters as well as they are dependent on the basis definition for each element type.

The user is responsible for creating and deleting the objects to avoid memory leaks.


### Python API

The Python API is disabled by default and can be built in-tree or out-of-tree.
The in-tree build is mainly meant for end users and packages.
To build the Python API with the normal project set the ``python`` option in the configuration step with

```sh
meson setup _build -Dpython=true -Dpython_version=$(which python3)
```

The Python version can be used to select a different Python version, it defaults to `'python3'`.
Python 2 is not supported with this project, the Python version key is meant to select between several local Python 3 versions.

Proceed with the build as described before and install the projects to make the Python API available in the selected prefix.

For the out-of-tree build see the instructions in the [``python``](./python) directory.


## Contributing

This is a volunteer open source projects and contributions are always welcome. 
Please, take a moment to read the [contributing guidelines](CONTRIBUTING.md).


## License

This project is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This project is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
GNU Lesser General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this project by you, as defined in the
GNU Lesser General Public license, shall be licensed as above, without any
additional terms or conditions.
