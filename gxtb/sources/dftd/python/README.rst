DFT-D4 Python API
-----------------

Python interface for the generally applicable atomic-charge dependent London dispersion correction, DFT-D4.
This Python project is targeted at developers who want to interface their project via Python with ``dftd4``.

This interface provides access to the C-API of ``dftd4`` via the CFFI module.
The low-level CFFI interface is available in the ``dftd4.library`` module and only required for implementing other interfaces.
A more pythonic interface is provided in the ``dftd4.interface`` module which can be used to build more specific interfaces.

.. code:: python

   >>> from dftd4.interface import DampingFunction, DampingParam, DispersionModel
   >>> import numpy as np
   >>> numbers = np.array([1, 1, 6, 5, 1, 15, 8, 17, 13, 15, 5, 1, 9, 15, 1, 15])
   >>> positions = np.array([  # Coordinates in Bohr
   ...     [+2.79274810283778, +3.82998228828316, -2.79287054959216],
   ...     [-1.43447454186833, +0.43418729987882, +5.53854345129809],
   ...     [-3.26268343665218, -2.50644032426151, -1.56631149351046],
   ...     [+2.14548759959147, -0.88798018953965, -2.24592534506187],
   ...     [-4.30233097423181, -3.93631518670031, -0.48930754109119],
   ...     [+0.06107643564880, -3.82467931731366, -2.22333344469482],
   ...     [+0.41168550401858, +0.58105573172764, +5.56854609916143],
   ...     [+4.41363836635653, +3.92515871809283, +2.57961724984000],
   ...     [+1.33707758998700, +1.40194471661647, +1.97530004949523],
   ...     [+3.08342709834868, +1.72520024666801, -4.42666116106828],
   ...     [-3.02346932078505, +0.04438199934191, -0.27636197425010],
   ...     [+1.11508390868455, -0.97617412809198, +6.25462847718180],
   ...     [+0.61938955433011, +2.17903547389232, -6.21279842416963],
   ...     [-2.67491681346835, +3.00175899761859, +1.05038813614845],
   ...     [-4.13181080289514, -2.34226739863660, -3.44356159392859],
   ...     [+2.85007173009739, -2.64884892757600, +0.71010806424206],
   ... ])
   >>> model = DispersionModel(numbers, positions)
   >>> damp = DampingFunction(model="d4")
   >>> param = DampingParam(method="scan", model="d4")
   >>> res = model.get_dispersion(damp, param, grad=False)
   >>> res.get("energy")  # Results in atomic units
   -0.005328888532435093
   >>> res.update(**model.get_properties())  # also allows access to properties
   >>> res.get("c6 coefficients")[0, 0]
   1.5976689760849156
   >>> res.get("polarizabilities")
   array([ 1.97521745,  1.48512704,  7.33564674, 10.28920458,  1.99973802,
          22.85298573,  6.65877552, 15.39410319, 22.73119177, 22.86303028,
          14.56038118,  1.4815783 ,  3.91266859, 25.8236368 ,  1.93444627,
          23.02494331])


Additional features
~~~~~~~~~~~~~~~~~~~

The ``dftd4.parameters`` module becomes available if a TOML parser is available, either `tomlkit <https://github.com/sdispater/tomlkit>`_ or `toml <https://github.com/uiri/toml>`_ can be used here.
The returned dict can be used to supply parameters to the constructor of the ``DampingParam`` object, only the ``s6``, ``s8``, ``s9``, ``a1``, ``a2`` and ``alp`` entries will be used, the remaining entries are meta data describing the damping parameters.

.. code-block:: python

   >>> from dftd4.parameters import get_damping_param
   >>> get_damping_param("b97m")
   {'s6': 1.0, 's9': 1.0, 'alp': 16.0, 's8': 0.6633, 'a1': 0.4288, 'a2': 3.9935}
   >>> get_damping_param("r2scan", keep_meta=True)
   {'s6': 1.0, 's9': 1.0, 'alp': 16.0, 'damping': 'bj', 'mbd': 'approx-atm', 's8': 0.6018749, 'a1': 0.51559235, 'a2': 5.77342911, 'doi': '10.1063/5.0041008'}


QCSchema Integration
~~~~~~~~~~~~~~~~~~~~

This Python API natively understands QCSchema and the `QCArchive infrastructure <http://docs.qcarchive.molssi.org>`_.
If the QCElemental package is installed the ``dftd4.qcschema`` module becomes importable and provides the ``run_qcschema`` function.

.. code:: python

   >>> from dftd4.qcschema import run_qcschema
   >>> import qcelemental as qcel
   >>> atomic_input = qcel.models.AtomicInput(
   ...     molecule = qcel.models.Molecule(
   ...         symbols = ["O", "H", "H"],
   ...         geometry = [
   ...             0.00000000000000,  0.00000000000000, -0.73578586109551,
   ...             1.44183152868459,  0.00000000000000,  0.36789293054775,
   ...            -1.44183152868459,  0.00000000000000,  0.36789293054775
   ...         ],
   ...     ),
   ...     driver = "energy",
   ...     model = {
   ...         "method": "TPSS-D4",
   ...     },
   ...     keywords = {},
   ... )
   ...
   >>> atomic_result = run_qcschema(atomic_input)
   >>> atomic_result.return_result
   -0.0002667885779142513


ASE Integration
~~~~~~~~~~~~~~~

To integrate with `ASE <https://wiki.fysik.dtu.dk/ase/>`_ this interface implements an ASE Calculator.
The ``DFTD4`` calculator becomes importable if an ASE installation is available.

.. code:: python

   >>> from ase.build import molecule
   >>> from dftd4.ase import DFTD4
   >>> atoms = molecule('H2O')
   >>> atoms.calc = DFTD4(method="TPSS")
   >>> atoms.get_potential_energy()
   -0.007310393443152083
   >>> atoms.calc.set(method="PBE")
   {'method': 'PBE'}
   >>> atoms.get_potential_energy()
   -0.005358475432239303
   >>> atoms.get_forces()
   array([[-0.        , -0.        ,  0.00296845],
          [-0.        ,  0.00119152, -0.00148423],
          [-0.        , -0.00119152, -0.00148423]])

To use the ``DFTD4`` calculator as dispersion correction the calculator can be combined using the `SumCalculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/mixing.html>`_ from the ``ase.calculators.mixing`` module.

.. code:: python

   >>> from ase.build import molecule
   >>> from ase.calculators.mixing import SumCalculator
   >>> from ase.calculators.nwchem import NWChem
   >>> from dftd4.ase import DFTD4
   >>> atoms = molecule('H2O')
   >>> atoms.calc = SumCalculator([DFTD4(method="PBE"), NWChem(xc="PBE")])

For convenience ``DFTD4`` allows to combine itself with another calculator by using the ``add_calculator`` method which returns a SumCalculator:

.. code:: python

   >>> from ase.build import molecule
   >>> from ase.calculators.emt import EMT
   >>> from dftd4.ase import DFTD4
   >>> atoms = molecule("C60")
   >>> atoms.calc = DFTD4(method="pbe").add_calculator(EMT())
   >>> atoms.get_potential_energy()
   6.348142387048062
   >>> [calc.get_potential_energy() for calc in atoms.calc.calcs]
   [-6.015477436263984, 12.363619823312046]

The individual contributions are available by iterating over the list of calculators in ``calc.calcs``.
Note that ``DFTD4`` will always place itself as first calculator in the list.


PySCF support
~~~~~~~~~~~~~

Integration with `PySCF <https://pyscf.org>`_ is possible by using the ``dftd4.pyscf`` module.
The module provides a ``DFTD4Dispersion`` class to construct a PySCF compatible calculator for evaluating the dispersion energy and gradients.

.. code:: python

   >>> from pyscf import gto
   >>> import dftd4.pyscf as disp
   >>> mol = gto.M(
   ...     atom="""
   ...          C   -0.755422531  -0.796459123  -1.023590391
   ...          C    0.634274834  -0.880017014  -1.075233285
   ...          C    1.406955202   0.199695367  -0.653144334
   ...          C    0.798863737   1.361204515  -0.180597909
   ...          C   -0.593166787   1.434312023  -0.133597923
   ...          C   -1.376239198   0.359205222  -0.553258516
   ...          I   -1.514344238   3.173268101   0.573601106
   ...          H    1.110906949  -1.778801728  -1.440619836
   ...          H    1.399172302   2.197767355   0.147412751
   ...          H    2.486417780   0.142466525  -0.689380574
   ...          H   -2.454252250   0.422581120  -0.512807958
   ...          H   -1.362353593  -1.630564523  -1.348743149
   ...          S   -3.112683203   6.289227834   1.226984439
   ...          H   -4.328789697   5.797771251   0.973373089
   ...          C   -2.689135032   6.703163830  -0.489062886
   ...          H   -1.684433029   7.115457372  -0.460265708
   ...          H   -2.683867206   5.816530502  -1.115183775
   ...          H   -3.365330613   7.451201412  -0.890098894
   ...          """
   ... )
   >>> d4 = disp.DFTD4Dispersion(mol, xc="r2SCAN")
   >>> d4.kernel()[0]
   array(-0.0050011)

To make use of the dispersion correction together with other calculators, the ``energy`` method allows to apply a dispersion correction to an existing calculator.

.. code:: python

   >>> from pyscf import gto, scf
   >>> import dftd4.pyscf as disp
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
   >>> mf = disp.energy(scf.RHF(mol)).run()
   converged SCF energy = -149.939098424774
   >>> grad = mf.nuc_grad_method().kernel()
   --------------- DFTD4 gradients ---------------
            x                y                z
   0 O     0.0172438133     0.0508406920     0.0000000000
   1 O     0.0380018285    -0.0460223790    -0.0000000000
   2 H    -0.0305058266    -0.0126478132    -0.0000000000
   3 H     0.0069233858    -0.0382898692    -0.0000000000
   4 H    -0.0158316004     0.0230596847     0.0218908543
   5 H    -0.0158316004     0.0230596847    -0.0218908543
   ----------------------------------------------


Installing
~~~~~~~~~~

.. image:: https://img.shields.io/conda/vn/conda-forge/dftd4-python.svg
   :alt: Conda Version
   :target: https://anaconda.org/conda-forge/dftd4-python

This project is packaged for the *conda* package manager and available on the *conda-forge* channel.
To install the *conda* package manager we recommend the `miniforge <https://github.com/conda-forge/miniforge/releases>`_ installer.
If the *conda-forge* channel is not yet enabled, add it to your channels with

.. code:: sh

   conda config --add channels conda-forge

Once the *conda-forge* channel has been enabled, this project can be installed with:

.. code:: sh

   conda install dftd4-python

Now you are ready to use ``dftd4``, check if you can import it with

.. code:: python

   >>> import dftd4
   >>> from dftd4.libdftd4 import get_api_version
   >>> get_api_version()
   '4.0.2'


Building the extension module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To perform an out-of-tree build some version of ``dftd4`` has to be available on your system and preferably findable by ``pkg-config``.
Try to find a ``dftd4`` installation you build against first with

.. code:: sh

   pkg-config --modversion dftd4

Adjust the ``PKG_CONFIG_PATH`` environment variable to include the correct directories to find the installation if necessary.


Using pip
^^^^^^^^^

This project support installation with pip as an easy way to build the Python API.

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

To install extra dependencies as well use

.. code:: sh

   pip install '.[parameters,ase,qcschema]'

If you already have a ``dftd4`` installation, *e.g.* from conda-forge, you can build the Python extension module directly without cloning this repository

.. code:: sh

   pip install "https://github.com/dftd4/dftd4/archive/refs/heads/main#egg=dftd4-python&subdirectory=python"



Using meson
^^^^^^^^^^^

This directory contains a separate meson build file to allow the out-of-tree build of the CFFI extension module.
The out-of-tree build requires

- C compiler to build the C-API and compile the extension module
- `meson <https://mesonbuild.com>`_ version 0.55 or newer
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

The extension module is now available in ``_build/dftd4/_libdftd4.*.so``.
You can install as usual with

.. code:: sh

   meson configure _build --prefix=/path/to/install
   meson install -C _build
