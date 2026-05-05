Atomistic features based on Mulliken partitioning
=================================================

Based on the electronic structure defined by the implemented xTB Hamiltonians in ``tblite``, 
we can compute atomistic properties that can be used for example in ML models or in general as chemical descriptors.

Currently, we feature the so-called ``xtbml`` properties, which can be organized into 4 categories:
  * :doc:`geometry-based properties <geometry>`
  * :doc:`density-based properties <density>`
  * :doc:`energy-based properties <energy>`
  * :doc:`orbital energy-based properties <orbital-energy>`

The details about the implemented properties can be found in the according sections.
We introduce the properties, with the equations defining them, and what key is used to store them in the output dictionary.

To capture non-local effects, we use a convolution kernel based on the coordination number of the surrounding atoms.
The convolution kernel is defined as a function of the distance between the atoms, and a scaling parameter ``a``.

.. math:: 

  f_{\text{log}}(R_A,R_B) = \frac{1}{1+\exp\left(\frac{-16a\cdot4(R_{A,\text{cov}}+R_{B,\text{cov}})}{3R_{AB}}-1\right)}

We refer to the properties computed using this convolution as "extended" properties. 
Such extended versions of properties are available for the geometry-, density- and orbital-energy-based properties.
The dictionary keys for those properties start with ``ext_``.

By default :math:`a = 1.0`, but it can be adjusted by the user, through a ``toml`` input file.
If multiple ``a`` values are provided, the value of a is added to the dictionary key, e.g. ``ext_CN_A_1.2``.

.. important::

    The properties are computed for each atom in the order of the input geometry.
    The properties are stored in the output file in the same order, which means that the properties are not permutation invariant as returned by the toml file.
    This can be achieved by sorting the properties according to a specific key.

How to compute properties
-------------------------

1. Using the CLI argument ``--post-processing``

The ``xtbml`` properties can be computed using CLI in the ``run`` mode.:

.. code-block:: bash

   tblite run --post-processing xtbml coord.xyz

This will return a ``toml`` file with the computed properties named ``post_processing.toml``.
In TOML files, the data is stored in key value pairs, the keys are defined in the detailed sections of the properties.
The values are returned as floating point arrays, since the properties are computed per atom.

Additionally, the properties can be computed using the ``xtbml_xyz`` post-processing string.
This means that the dipole moments are returned as three entries in the result dictionary, suffixed by ``_x``, ``_y``, and ``_z``.
The quadrupole moments are returned as six entries, with the entries suffixed by ``_xx``, ``_xy``, ``_xz``, ``_yy``, ``_yz``, and ``_zz``.
Of course, this breaks the rotation invariance of the properties, but can be useful for some applications.

2. Using a ``toml`` input file
   
To select individual properties to be computed, a ``toml`` file can be used as input.
Currently, the properties can be selected per category only.
The properties are selected by setting the corresponding key to ``true`` in the input file.

.. code-block:: toml

  [post-processing.xtbml]
  geometry = false
  density = true
  orbital = true
  energy = false
  convolution = true
  a = [1.0, 1.2]
  tensorial-output = true

``a`` is the scaling parameter for the convolution kernel, and ``tensorial-output`` is a flag to return the multipole moments as tensors instead of norms (see above).

Similar to the ``xtbml`` and ``xtbml_xyz`` post-processing strings, we can use the name of the ``toml`` file as an argument for the CLI:

.. code-block:: bash

   tblite run --post-processing xtbml.toml coord.xyz

The properties are written to the ``post_processing.toml`` file again.


3. Using the Python API

The properties can be computed using the Python API.
Similar to the CLI, the properties can be computed using the ``xtbml`` or ``xtbml_xyz`` post-processing strings or by providing a ``toml`` file as input.
Those strings are handed to the calculators ``add()`` method.
An example is shown below:

.. code-block:: python

  from tblite.interface import Calculator
  import numpy as np

  calc = Calculator(
      method="GFN2-xTB",
      numbers=np.array([1, 1]),
      positions=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
  )
  calc.add('xtbml')
  #calc.add('xtbml_xyz') with tensorial output
  #calc.add('xtbml.toml') use toml as input
  res = calc.singlepoint()

  dict_xtbml = res.get('post-processing-dict')
  
  
.. toctree:: 
    :maxdepth: 2
    
    geometry
    density
    energy
    orbital-energy