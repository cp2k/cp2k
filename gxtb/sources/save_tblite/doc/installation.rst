.. _install:

Installation
============

This project is currently in a highly experimental stage, however we can already offer binary distributions via some channels.


Installing from conda
---------------------

This project is packaged for the *conda* package manager and available on the *conda-forge* channel.
To install the *conda* package manager we recommend the `miniforge <https://github.com/conda-forge/miniforge/releases>`_ installer.
If the *conda-forge* channel is not yet enabled, add it to your channels with

.. code-block:: bash

    conda config --add channels conda-forge

Once the *conda-forge* channel has been enabled, this project can be installed with:

.. code-block:: bash

   conda install tblite

It is possible to list all of the versions available on your platform with:

.. code-block:: bash

   conda search tblite --channel conda-forge

Now you are ready to use the ``tblite`` executable, find the ``tblite.h`` header and link against the ``tblite`` library.


FreeBSD Port
------------

A port for FreeBSD is available

.. code-block:: bash

   pkg install science/tblite

In case no package is available build the port using

.. code-block:: bash

   cd /usr/ports/science/tblite
   make install clean

For more information see the `tblite port details <https://www.freshports.org/science/tblite/>`_.


Building from source
--------------------

This library depends on several Fortran modules to provide the desired functionality

- `mctc-lib`_: Modular computation tool chain library
- `dftd4`_: Reference implementation of the generally applicable charge-dependent London-dispersion correction, DFT-D4
- `s-dftd3`_: Reimplementation of the DFT-D3 dispersion correction
- `mstore`_: Molecular structure store (testing only)
- `toml-f`_: Library for processing and emitting TOML data

.. _dftd4: https://github.com/dftd4/dftd4
.. _s-dftd3: https://github.com/awvwgk/simple-dftd3
.. _multicharge: https://github.com/grimme-lab/multicharge
.. _mctc-lib: https://github.com/grimme-lab/mctc-lib
.. _mstore: https://github.com/grimme-lab/mstore
.. _toml-f: https://github.com/toml-f/toml-f

.. _meson: https://mesonbuild.com
.. _ninja: https://ninja-build.org
.. _asciidoctor: https://asciidoctor.org
.. _cmake: https://cmake.org
.. _fpm: https://fpm.fortran-lang.org
.. _cffi: https://cffi.readthedocs.io/
.. _numpy: https://numpy.org/
.. _pkgconfig: https://pypi.org/project/pkgconfig/


.. _meson-build:

Meson based build
~~~~~~~~~~~~~~~~~

The primary build system of this project is `meson`_.
For the full feature-complete build it is highly recommended to perform the build and development with meson.
To setup a build the following software is required

- A Fortran 2008 compliant compiler, like GCC Fortran and Intel Fortran classic
- `meson`_, version 0.57.2 or newer
- `ninja`_, version 1.10 or newer
- a linear algebra backend, like MKL or OpenBLAS

Optional dependencies are

- `asciidoctor`_, to build the manual pages
- A C compiler to test the C API and compile the Python extension module
- Python 3.6 or newer with the `CFFI`_ package installed to build the Python API

To setup a new build run

.. code:: text

   meson setup _build --prefix=$HOME/.local

The Fortran and C compiler can be selected with the ``FC`` and ``CC`` environment variable, respectively.
The installation location is selected using the ``--prefix`` option.
The required Fortran modules will be fetched automatically from the upstream repositories and checked out in the *subprojects* directory.

.. note::

   For Intel Fortran oneAPI (2021 or newer) builds with MKL backend the ``-Dfortran_link_args=-qopenmp`` option has to be added.

.. tip::

   To produce statically linked binaries set ``--default-library=static`` and add ``-Dfortran_link_args=-static`` as well.

To compile the project run

.. code:: text

   meson compile -C _build

Verify that the resulting projects is working correctly by running the testsuite with

.. code:: text

   meson test -C _build --print-errorlogs

In case meson is spawning too many concurrent test jobs limit the number of processes with the ``--num-processes`` option when running the test command.
By default the whole library and its subprojects are tested, to limit the testing to the project itself add ``--suite tblite`` as option.

To verify the included parametrizations are working correctly run the extra testsuite by passing the ``--benchmark`` argument

.. code:: text

   meson test -C _build --print-errorlogs --benchmark

Finally, you can make the project available by installation with

.. code:: text

   meson install -C _build


CMake based build
~~~~~~~~~~~~~~~~~

This project also provides support for `CMake`_ to give projects using it as build system an easier way to interface.
The CMake build files usually do not provide a feature-complete build, but contributions are more than welcome.
To setup a build the following software is required

- A Fortran 2008 compliant compiler, like GCC Fortran and Intel Fortran classic
- `cmake`_, version 3.14 or newer
- `ninja`_, version 1.10 or newer
- a linear algebra backend, like MKL or OpenBLAS

Configure a new build with

.. code:: text

   cmake -B _build -G Ninja -DCMAKE_INSTALL_PREFIX=$HOME/.local

You can set the Fortran compiler in the ``FC`` environment variable.
The installation location can be selected with the ``CMAKE_INSTALL_PREFIX``, GNU install directories are supported by default.
CMake will automatically fetch the required Fortran modules, you can provide specific version in the *subprojects* directory which will be used instead.

To run a build use

.. code:: text

   cmake --build _build

After a successful build make sure the testsuite passes

.. code:: text

   pushd _build && ctest --output-on-failure && popd

To make the project available install it with

.. code:: text

   cmake --install _build


Fpm based build
~~~~~~~~~~~~~~~

This projects supports building with the Fortran package manager (`fpm`_).
Create a new build by running

.. code:: text

   fpm build

You can adjust the Fortran compiler with the ``--compiler`` option and select the compilation profile with ``--profile release``.
To test the resulting build run the testsuite with

.. code:: text

   fpm test

The command line driver can be directly used from fpm wih

.. code:: text

   fpm run --profile release -- --help

To make the installation accessible install the project with

.. code:: text

   fpm install --profile release --prefix $HOME/.local


.. _python-build:

Python extension module
-----------------------

The Python API is available as Python extension module.
The easiest way to setup is to add ``-Dpython=true`` to a meson tree build and follow the :ref:`meson installation instructions <meson-build>`.
The extension module will become available once the project is installed.

.. important::

   When building with Intel compilers make sure to use the real-time version of the MKL.
   Add ``-Dlapack=mkl-rt`` when configuring the build.
   Otherwise, when using the normal MKL libraries dynamically loading the *tblite* library from Python will fail.

This section describes alternative ways to build the Python API


Using pip
~~~~~~~~~

This project support installation with pip as an easy way to build the Python API.

- C compiler to build the C-API and compile the extension module (the compiler name should be exported in the ``CC`` environment variable)
- Python 3.6 or newer
- The following Python packages are required additionally

  - `cffi`_
  - `numpy`_
  - `pkgconfig`_ (setup only)

Make sure to have your C compiler set to the ``CC`` environment variable

.. code:: sh

   export CC=gcc

Install the project with pip

.. code:: sh

   pip install .

To install extra dependencies as well use

.. code:: sh

   pip install '.[ase]'


Using meson
~~~~~~~~~~~

The Python extension module can be built on-top of an existing installation, either provided by meson or CMake.

Building requires against an existing *tblite* installation requires

- C compiler to build the C-API and compile the extension module
- `meson`_ version 0.55 or newer
- a build-system backend, *i.e.* `ninja`_ version 1.7 or newer
- Python 3.6 or newer with the `CFFI`_ package installed

Setup a build with

.. code:: sh

   meson setup _build_python python -Dpython_version=$(which python3)

The Python version can be used to select a different Python version, it defaults to ``'python3'``.
Python 2 is not supported with this project, the Python version key is meant to select between several local Python 3 versions.

Compile the project with

.. code:: sh

   meson compile -C _build

The extension module is now available in ``_build_python/tblite/_libtblite.*.so``.
You can install as usual with

.. code:: sh

   meson configure _build --prefix=/path/to/install
   meson install -C _build


Supported Compilers
-------------------

This is a non-comprehensive list of tested compilers for *tblite*.
Compilers with the label *latest* are tested with continuous integration for each commit.

=========== ================= ================ ================== =============================
 Compiler    Version           Platform         Architecture       *tblite*
=========== ================= ================ ================== =============================
 GCC         11.1, 10.3        Ubuntu 20.04     x86_64             0.2.0, 0.2.1, 0.3.0, latest
 GCC               10.3, 9.4   MacOS 11.6.5     x86_64             0.2.0, 0.2.1, 0.3.0, latest
 GCC                     9.4   MacOS 10.15.7    x86_64             0.2.0, 0.2.1
 GCC         11.0              MacOS 11.0       arm64              0.2.0, 0.2.1
 GCC               10.3        CentOS 7         aarch64, ppc64le   0.2.0, 0.2.1
 GCC/MinGW   11.2              Windows 2022     x86_64                    0.2.1  0.3.0, latest
 Intel       2021.2            Ubuntu 20.04     x86_64             0.2.0, 0.2.1, 0.3.0, latest
 NAG         7.1               AlmaLinux 8.5    x86_64             0.2.0, 0.2.1
=========== ================= ================ ================== =============================

Compiler known to fail are documented here, together with the last commit where this behaviour was encountered.
If available an issue in on the projects issue tracker or the issue tracker of the dependencies is linked.
Usually, it safe to assume that older versions of the same compiler will fail to compile as well and this failure is consistent over platforms and/or architectures.

========== ============= =============== ============== ==========================
 Compiler   Version       Platform        Architecture   Reference
========== ============= =============== ============== ==========================
 GCC        6.4.0         MacOS 10.15.7   x86_64         `abb17c3`_
 Intel      19.0.5        AlmaLinux 8.5   x86_64         `0542ce7`_, `tblite#45`_
 Intel      17.0.1        OpenSuse 42.1   x86_64         `abb17c3`_, `tblite#2`_
 Intel      16.0.3        CentOS 7.3      x86_64         `abb17c3`_, `dftd4#112`_
 Flang      20190329      Ubuntu 20.04    x86_64         `abb17c3`_, `toml-f#28`_
 NVHPC      20.9          Manjaro Linux   x86_64         `abb17c3`_, `toml-f#27`_
========== ============= =============== ============== ==========================

.. _0542ce7: https://github.com/tblite/tblite/tree/0542ce7ae0e323941156949a0620ca260bc0ce7f
.. _abb17c3: https://github.com/tblite/tblite/tree/abb17c3a8ea8e0336dde84ed78bdab8033144a0a
.. _tblite#2: https://github.com/tblite/tblite/issues/2
.. _dftd4#112: https://github.com/dftd4/dftd4/issues/112
.. _toml-f#28: https://github.com/toml-f/toml-f/issues/28
.. _toml-f#27: https://github.com/toml-f/toml-f/issues/27
.. _tblite#45: https://github.com/tblite/tblite/issues/45
