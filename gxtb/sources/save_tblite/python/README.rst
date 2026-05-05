Python binding for tblite
=========================

This directory contains the Python API of the *tblite* project.

This interface provides access to the C-API of ``tblite`` via the CFFI module.
The low-level CFFI interface is available in the ``tblite.library`` module and only required for implementing other interfaces.
A more pythonic interface is provided in the ``tblite.interface`` module which can be used to build more specific interfaces.

.. code:: python

   from tblite.interface import Calculator
   import numpy as np
   numbers = np.array([1, 1, 6, 5, 1, 15, 8, 17, 13, 15, 5, 1, 9, 15, 1, 15])
   positions = np.array([  # Coordinates in Bohr
       [+2.79274810283778, +3.82998228828316, -2.79287054959216],
       [-1.43447454186833, +0.43418729987882, +5.53854345129809],
       [-3.26268343665218, -2.50644032426151, -1.56631149351046],
       [+2.14548759959147, -0.88798018953965, -2.24592534506187],
       [-4.30233097423181, -3.93631518670031, -0.48930754109119],
       [+0.06107643564880, -3.82467931731366, -2.22333344469482],
       [+0.41168550401858, +0.58105573172764, +5.56854609916143],
       [+4.41363836635653, +3.92515871809283, +2.57961724984000],
       [+1.33707758998700, +1.40194471661647, +1.97530004949523],
       [+3.08342709834868, +1.72520024666801, -4.42666116106828],
       [-3.02346932078505, +0.04438199934191, -0.27636197425010],
       [+1.11508390868455, -0.97617412809198, +6.25462847718180],
       [+0.61938955433011, +2.17903547389232, -6.21279842416963],
       [-2.67491681346835, +3.00175899761859, +1.05038813614845],
       [-4.13181080289514, -2.34226739863660, -3.44356159392859],
       [+2.85007173009739, -2.64884892757600, +0.71010806424206],
   ])
   calc = Calculator("GFN2-xTB", numbers, positions)
   res = calc.singlepoint()
   print(res.get("energy"))  # Results in atomic units
   # => -31.716159156026254


Building the extension module
-----------------------------

The Python bindings can be built against an existing installation of ``tblite`` or free-standing.
The free-standing implementation will select a matching version of the shared library, when building against an existing ``tblite`` library the API version of the two parts must match.


Setuptools build
~~~~~~~~~~~~~~~~

This project support installation with pip as an easy way to build the Python API.

- C compiler to build the C-API and compile the extension module (the compiler name should be exported in the ``CC`` environment variable)
- Python 3.6 or newer
- The following Python packages are required additionally

  - `cffi <https://cffi.readthedocs.io/>`_
  - `numpy <https://numpy.org/>`_
  - `pkgconfig <https://pypi.org/project/pkgconfig/>`_ (setup only)

Ensure that you can find ``tblite`` via

.. code:: sh

   pkg-config --modversion tblite

Adjust the ``PKG_CONFIG_PATH`` environment variable to include the correct directories to find the installation if necessary.
Alternatively, you can set the ``TBLITE_PREFIX`` environment variable to point to the installation of the library.

Make sure to have your C compiler set to the ``CC`` environment variable

.. code:: sh

   export CC=gcc

Install the project with pip

.. code:: sh

   pip install .


Using meson
~~~~~~~~~~~

This directory contains a separate meson build file to allow the out-of-tree build of the CFFI extension module.
The out-of-tree build requires

- C compiler to build the C-API and compile the extension module
- `meson <https://mesonbuild.com>`_ version 0.57.2 or newer
- a build-system backend, *i.e.* `ninja <https://ninja-build.org>`_ version 1.7 or newer
- Python 3.6 or newer with the `CFFI <https://cffi.readthedocs.io/>`_ package installed

To make a free-standing build you can provide the main repository as subproject to the Python bindings *without* having to build the shared library first.
This can be done for example by symlinking the main repository to the subprojects directory.

.. code:: sh

   mkdir subprojects
   ln -s $(realpath ..) subprojects/tblite

Note that this step is not needed if you built against an existing ``tblite`` installation.

Setup a build with

.. code:: sh

   meson setup _build -Dpython_version=$(which python3) --prefix=/path/to/install

The Python version can be used to select a different Python version, it defaults to ``'python3'``.
Python 2 is not supported with this project, the Python version key is meant to select between several local Python 3 versions.

Compile the project with

.. code:: sh

   meson compile -C _build

The extension module is now available in ``_build/tblite/_libtblite.*.so``.
You can install as usual with

.. code:: sh

   meson install -C _build
