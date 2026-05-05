.. _install:

Installation
============

The preferred method for installing this package is to use a binary release distributed with a package manager.
If your preferred package manager is not available, you can also download the source code and build it yourself.


Installing from conda-forge
---------------------------

This project is packaged for the *mamba* package manager and available on the *conda-forge* channel.
To install the *mamba* package manager we recommend the `miniforge <https://github.com/conda-forge/miniforge/releases>`_ installer.
If the *conda-forge* channel is not yet enabled, add it to your channels with

.. code-block:: bash

    mamba config --add channels conda-forge

Once the *conda-forge* channel has been enabled, this project can be installed with:

.. code-block:: bash

   mamba install simple-dftd3

If you want to enable the Python API as well install

.. code-block:: bash

   mamba install dftd3-python

It is possible to list all of the versions available on your platform with:

.. code-block:: bash

   mamba search simple-dftd3 --channel conda-forge

Now you are ready to use the ``s-dftd3`` executable, find the ``dftd3.h`` header or import ``dftd3`` in your Python projects.


Building from source
--------------------

This library depends on few Fortran modules to provide the desired functionality

- `mctc-lib`_: Modular computation tool chain library
- `mstore`_: Molecular structure store (testing only)
- `toml-f`_: TOML parser for Fortran

When building from source the dependencies will be downloaded automatically if not available on the system.
To build offline install the dependencies before building the library or provide the source of the dependencies in the *subprojects* directory.

.. note::

   Using the released archive is recommended because it already bundles the source code of all dependencies.
   The complete source archive in the release page is named ``s-dftd3-<version>.tar.xz``.

.. _mctc-lib: https://github.com/grimme-lab/mctc-lib
.. _mstore: https://github.com/grimme-lab/mstore
.. _toml-f: https://toml-f.readthedocs.io

.. _meson: https://mesonbuild.com
.. _ninja: https://ninja-build.org
.. _asciidoctor: https://asciidoctor.org
.. _cmake: https://cmake.org
.. _fpm: https://github.com/fortran-lang/fpm
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
- `meson`_, version 0.55 or newer
- `ninja`_, version 1.7 or newer

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


.. tip::

   To produce statically linked binaries set ``--default-library=static`` and add ``-Dfortran_link_args=-static`` as well.

To compile the project run

.. code:: text

   meson compile -C _build

Verify that the resulting projects is working correctly by running the testsuite with

.. code:: text

   meson test -C _build --print-errorlogs

In case meson is spawning too many concurrent test jobs limit the number of processes with the ``--num-processes`` option when running the test command.
By default the whole library and its subprojects are tested, to limit the testing to the project itself add ``--suite s-dftd3`` as option.

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

   ctest --output-on-failure --test-dir _build

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

Building requires against an existing *s-dftd3* installation requires

- C compiler to build the C-API and compile the extension module
- `meson`_ version 0.55 or newer
- a build-system backend, *i.e.* `ninja`_ version 1.7 or newer
- Python 3.6 or newer with the `CFFI`_ package installed

Setup a build with

.. code:: sh

   meson setup _build_python python -Dpython_version=3

The Python version can be used to select a different Python version, it defaults to ``'3'``.
Python 2 is not supported with this project, the Python version key is meant to select between several local Python 3 versions.

Compile the project with

.. code:: sh

   meson compile -C _build

The extension module is now available in ``_build_python/dftd3/_libdftd3.*.so``.
You can install as usual with

.. code:: sh

   meson configure _build --prefix=/path/to/install
   meson install -C _build
