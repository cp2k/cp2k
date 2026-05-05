DFT-D3 Python API
=================

Python interface for the D3 dispersion model.
This Python project is targeted at developers who want to interface their project via Python with ``s-dftd3``.

This interface provides access to the C-API of ``s-dftd3`` via the CFFI module.
The low-level CFFI interface is available in the ``dftd3.libdftd3`` module and only required for implementing other interfaces.
A more pythonic interface is provided in the ``dftd3.interface`` module which can be used to build more specific interfaces.

.. code:: python

   from dftd3.interface import RationalDampingParam, DispersionModel
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
   model = DispersionModel(numbers, positions)
   res = model.get_dispersion(RationalDampingParam(method="pbe0"), grad=False)
   print(res.get("energy"))  # Results in atomic units
   # => -0.029489232932494884


QCSchema Integration
--------------------

This Python API natively understands QCSchema and the `QCArchive infrastructure <http://docs.qcarchive.molssi.org>`_.
If the QCElemental package is installed the ``dftd3.qcschema`` module becomes importable and provides the ``run_qcschema`` function.

.. code:: python

   from dftd3.qcschema import run_qcschema
   import qcelemental as qcel
   atomic_input = qcel.models.AtomicInput(
       molecule = qcel.models.Molecule(
           symbols = ["O", "H", "H"],
           geometry = [
               0.00000000000000,  0.00000000000000, -0.73578586109551,
               1.44183152868459,  0.00000000000000,  0.36789293054775,
              -1.44183152868459,  0.00000000000000,  0.36789293054775
           ],
       ),
       driver = "energy",
       model = {
           "method": "tpss",
       },
       keywords = {
           "level_hint": "d3bj",
       },
   )

   atomic_result = run_qcschema(atomic_input)
   print(atomic_result.return_result)
   # => -0.0004204244108151285


ASE Integration
---------------

To integrate with `ASE <https://wiki.fysik.dtu.dk/ase/>`_ this interface implements an ASE Calculator.
The ``DFTD3`` calculator becomes importable if an ASE installation is available.

.. code:: python

   >>> from ase.build import molecule
   >>> from dftd3.ase import DFTD3
   >>> atoms = molecule('H2O')
   >>> atoms.calc = DFTD3(method="TPSS", damping="d3bj")
   >>> atoms.get_potential_energy()
   -0.0114416338147162
   >>> atoms.calc.set(method="PBE")
   {'method': 'PBE'}
   >>> atoms.get_potential_energy()
   -0.009781913226281063
   >>> atoms.get_forces()
   array([[-0.00000000e+00 -0.00000000e+00  9.56568982e-05]
          [-0.00000000e+00 -4.06046858e-05 -4.78284491e-05]
          [-0.00000000e+00  4.06046858e-05 -4.78284491e-05]])

To use the ``DFTD3`` calculator as dispersion correction the calculator can be combined using the `SumCalculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/mixing.html>`_ from the ``ase.calculators.mixing`` module.

.. code:: python

   >>> from ase.build import molecule
   >>> from ase.calculators.mixing import SumCalculator
   >>> from ase.calculators.nwchem import NWChem
   >>> from dftd3.ase import DFTD3
   >>> atoms = molecule('H2O')
   >>> atoms.calc = SumCalculator([DFTD3(method="PBE", damping="d3bj"), NWChem(xc="PBE")])

For convenience ``DFTD3`` allows to combine itself with another calculator by using the ``add_calculator`` method which returns a SumCalculator:

.. code:: python

   >>> from ase.build import molecule
   >>> from ase.calculators.emt import EMT
   >>> from dftd3.ase import DFTD3
   >>> atoms = molecule("C60")
   >>> atoms.calc = DFTD3(method="pbe", damping="d3bj").add_calculator(EMT())
   >>> atoms.get_potential_energy()
   7.513593999944228
   >>> [calc.get_potential_energy() for calc in atoms.calc.calcs]
   [-4.850025823367818, 12.363619823312046]

The individual contributions are available by iterating over the list of calculators in ``calc.calcs``.
Note that ``DFTD3`` will always place itself as first calculator in the list.


PySCF support
-------------

Integration with `PySCF <https://pyscf.org>`_ is possible by using the ``dftd3.pyscf`` module.
The module provides a ``DFTD3Dispersion`` class to construct a PySCF compatible calculator for evaluating the dispersion energy and gradients.

.. code:: python

   >>> from pyscf import gto
   >>> import dftd3.pyscf as disp
   >>> mol = gto.M(
   ...     atom="""
   ...          C   -0.189833176  -0.645396435   0.069807761
   ...          C    1.121636324  -0.354065576   0.439096514
   ...          C    1.486520953   0.962572632   0.712107225
   ...          C    0.549329390   1.989209324   0.617868956
   ...          C   -0.757627135   1.681862630   0.246856908
   ...          C   -1.138190460   0.370551816  -0.028582325
   ...          Br  -2.038462778   3.070459841   0.115165429
   ...          H    1.852935245  -1.146434699   0.514119204
   ...          H    0.825048723   3.012176989   0.829385472
   ...          H    2.502259769   1.196433556   1.000317333
   ...          H   -2.157140187   0.151608161  -0.313181471
   ...          H   -0.480820487  -1.664983631  -0.142918416
   ...          S   -4.157443472   5.729584377  -0.878761129
   ...          H   -4.823791426   4.796089466  -1.563433338
   ...          C   -2.828338520   5.970593053  -2.091189515
   ...          H   -2.167577293   6.722356639  -1.668621815
   ...          H   -2.264954814   5.054835899  -2.240198499
   ...          H   -3.218524904   6.337447714  -3.035087058
   ...          """
   ... )
   >>> d3 = disp.DFTD3Dispersion(mol, xc="PW6B95", version="d3bj")
   >>> d3.kernel()[0]
   array(-0.01009386)
   >>> d3.version = "d3zero"  # Change to zero damping
   >>> d3.kernel()[0]
   array(-0.00574098)
   >>> d3.atm = True  # Activate three-body dispersion
   >>> d3.kernel()[0]
   array(-0.00574289)

To make use of the dispersion correction together with other calculators, the ``energy`` method allows to apply a dispersion correction to an existing calculator.

.. code:: python

   >>> from pyscf import gto, scf
   >>> import dftd3.pyscf as disp
   >>> mol = gto.M(
   ...     atom="""
   ...          O  -1.65542061  -0.12330038   0.00000000
   ...          O   1.24621244   0.10268870   0.00000000
   ...          H  -0.70409026   0.03193167   0.00000000
   ...          H  -2.03867273   0.75372294   0.00000000
   ...          H   1.57598558  -0.38252146  -0.75856129
   ...          H   1.57598558  -0.38252146   0.75856129
   ...          """
   ... )
   >>> grad = disp.energy(scf.RHF(mol)).run().nuc_grad_method()
   converged SCF energy = -149.947191000075
   >>> g = grad.kernel()
   --------------- DFTD3 gradients ---------------
            x                y                z
   0 O     0.0171886976     0.0506606246     0.0000000000
   1 O     0.0383596853    -0.0459057549     0.0000000000
   2 H    -0.0313133974    -0.0125865676    -0.0000000000
   3 H     0.0066705789    -0.0380501872     0.0000000000
   4 H    -0.0154527822     0.0229409425     0.0215141991
   5 H    -0.0154527822     0.0229409425    -0.0215141991
   ----------------------------------------------


Installing
----------

.. image:: https://img.shields.io/conda/vn/conda-forge/dftd3-python.svg
   :alt: Conda Version
   :target: https://anaconda.org/conda-forge/dftd3-python

This project is packaged for the *conda* package manager and available on the *conda-forge* channel.
To install the *conda* package manager we recommend the `miniforge <https://github.com/conda-forge/miniforge/releases>`_ installer.
If the *conda-forge* channel is not yet enabled, add it to your channels with

.. code:: sh

   conda config --add channels conda-forge

Once the *conda-forge* channel has been enabled, this project can be installed with:

.. code:: sh

   conda install dftd3-python

Now you are ready to use ``dftd3``, check if you can import it with

.. code:: python

   >>> import dftd3
   >>> from dftd3.library import get_api_version
   >>> get_api_version()
   '1.3.2'


Building the extension module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To perform an out-of-tree build some version of ``s-dftd3`` has to be available on your system and preferably findable by ``pkg-config``.
Try to find a ``s-dftd3`` installation you build against first with

.. code:: sh

   pkg-config --modversion s-dftd3

Adjust the ``PKG_CONFIG_PATH`` environment variable to include the correct directories to find the installation if necessary.


Using pip
^^^^^^^^^

.. image:: https://img.shields.io/pypi/v/dftd3
   :target: https://pypi.org/project/dftd3/
   :alt: PyPI

This project support installation with pip as an easy way to build the Python API.
Precompiled Python wheels for Linux are available on `pypi <https://pypi.org/project/dftd3/>`_ and can be installed with

.. code:: sh

   pip install dftd3

Other platforms need to build from source, the following dependencies are required to do so

- C compiler to build the C-API and compile the extension module (the compiler name should be exported in the ``CC`` environment variable)
- Python 3.6 or newer
- The following Python packages are required additionally

  - `cffi <https://cffi.readthedocs.io/>`_
  - `numpy <https://numpy.org/>`_
  - `pkgconfig <https://pypi.org/project/pkgconfig/>`_ (setup only)

Make sure to have your C compiler set to the ``CC`` environment variable

.. code:: sh

   export CC=gcc

Install the project with pip

.. code:: sh

   pip install .

If you already have a ``s-dftd3`` installation, *e.g.* from conda-forge, you can build the Python extension module directly without cloning this repository

.. code:: sh

   pip install "https://github.com/dftd3/simple-dftd3/archive/refs/heads/main.zip#egg=dftd3-python&subdirectory=python"



Using meson
^^^^^^^^^^^

This directory contains a separate meson build file to allow the out-of-tree build of the CFFI extension module.
The out-of-tree build requires

- C compiler to build the C-API and compile the extension module
- `meson <https://mesonbuild.com>`_ version 0.53 or newer
- a build-system backend, *i.e.* `ninja <https://ninja-build.org>`_ version 1.7 or newer
- Python 3.6 or newer with the `CFFI <https://cffi.readthedocs.io/>`_ package installed

Setup a build with

.. code:: sh

   meson setup _build -Dpython_version=$(which python3)

The Python version can be used to select a different Python version, it defaults to ``'python3'``.
Python 2 is not supported with this project, the Python version key is meant to select between several local Python 3 versions.

Compile the project with

.. code:: sh

   meson compile -C _build

The extension module is now available in ``_build/dftd3/_libdftd3.*.so``.
You can install as usual with

.. code:: sh

   meson configure _build --prefix=/path/to/install
   meson install -C _build
