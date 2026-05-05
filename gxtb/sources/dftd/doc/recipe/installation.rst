.. _installation:

Installing DFT-D4
=================

This guide will walk you through installing the latest version of DFT-D4.


:fab:`apple` :fab:`linux` :fab:`windows` Installing from conda-forge
--------------------------------------------------------------------

.. image:: https://img.shields.io/conda/vn/conda-forge/dftd4
   :alt: Conda
   :target: https://github.com/conda-forge/dftd4-feedstock

.. image:: https://img.shields.io/conda/pn/conda-forge/dftd4
   :alt: Conda
   :target: https://github.com/conda-forge/dftd4-feedstock


This project is packaged for the *conda* package manager and available on the *conda-forge* channel.
To install the *conda* package manager we recommend the `miniforge <https://github.com/conda-forge/miniforge/releases>`_ installer.
If the *conda-forge* channel is not yet enabled, add it to your channels with

.. code-block:: bash

   conda config --add channels conda-forge
   conda config --set channel_priority strict

Once the *conda-forge* channel has been enabled, DFT-D4 can be installed with *conda*:

.. code-block:: shell

   conda install dftd4

or with *mamba*:

.. code-block:: shell

   mamba install dftd4

It is possible to list all of the versions of DFT-D4 available on your platform with *conda*:

.. code-block:: shell

   conda search dftd4 --channel conda-forge

or with *mamba*:

.. code-block:: shell

   mamba search dftd4 --channel conda-forge

Alternatively, *mamba repoquery* may provide more information:

.. code-block:: shell

   # Search all versions available on your platform:
   mamba repoquery search dftd4 --channel conda-forge

   # List packages depending on `dftd4`:
   mamba repoquery whoneeds dftd4 --channel conda-forge

   # List dependencies of `dftd4`:
   mamba repoquery depends dftd4 --channel conda-forge


:fab:`freebsd` FreeBSD ports
----------------------------

.. image:: https://repology.org/badge/version-for-repo/freebsd/dftd4.svg
   :alt: FreeBSD
   :target: https://www.freshports.org/science/dftd4/

A port for FreeBSD is available

.. code-block:: bash

   pkg install science/dftd4

In case no package is available build the port using

.. code-block:: bash

   cd /usr/ports/science/dftd4
   make install clean

For more information see the `dftd4 port details <https://www.freshports.org/science/dftd4/>`_.


Building from source
--------------------

To build this project from the source code in this repository you need to have

- a Fortran compiler supporting Fortran 2008
- One of the supported build systems

  - `meson <https://mesonbuild.com>`_ version 0.55 or newer
  - `CMake <https://cmake.org/>`_ version 3.9 or newer

.. note::
  
   GCC versions 15.0.x–15.1.x contain a `compiler bug <https://gcc.gnu.org/bugzilla/show_bug.cgi?id=119928>`__ that can lead to "Interface mismatch" errors during compilation.


First, get the source by cloning the repository

.. code-block:: bash

   git clone https://github.com/dftd4/dftd4
   cd dftd4


Using Meson
^^^^^^^^^^^

To build this project with meson a build-system backend is required, *i.e.* `ninja <https://ninja-build.org>`_ version 1.7 or newer.
Setup a build with

.. code-block:: bash

   meson setup _build --prefix=/path/to/installation

You can select the Fortran compiler by the ``FC`` environment variable.
To compile the project run

.. code-block:: bash

   meson compile -C _build

DFT-D4 comes with a comprehensive test suite.
Run the tests with

.. code-block:: bash

   meson test -C _build --print-errorlogs

Finally, you can install DFT-D4 with

.. code-block:: bash

   meson install -C _build


Using CMake
^^^^^^^^^^^

While meson is the preferred way to build this project it also offers CMake support.
Configure the CMake build with

.. code-block:: bash

   cmake -B_build -GNinja -DCMAKE_INSTALL_PREFIX=/path/to/installation

Similar to meson the compiler can be selected with the ``FC`` environment variable.
You can build the project using

.. code-block:: bash

   cmake --build _build

DFT-D4 comes with a comprehensive test suite.
Run the tests with

.. code-block:: bash

   ctest --test-dir _build --output-on-failure

Finally, you can install DFT-D4 with

.. code-block:: bash

   cmake --install _build
