How to use this library?
========================

This section contains a few self-contained examples on how to use D3.


Compute energy with rational damping
------------------------------------

This example shows how to compute the dispersion energy with the rational damping function.

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. literalinclude:: minimal-example/energy.f90
         :language: fortran
         :caption: energy.f90

   .. tab-item:: C
      :sync: c

      .. literalinclude:: minimal-example/energy.c
         :language: c
         :caption: energy.c

   .. tab-item:: Python
      :sync: python

      .. literalinclude:: minimal-example/energy.py
         :language: python
         :caption: energy.py

To test this example you can install the dependencies with

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: text

         mamba create d3 simple-dftd3 fortran-compiler pkg-config
         mamba activate d3

   .. tab-item:: C
      :sync: c

      .. code-block:: text

         mamba create d3 simple-dftd3 c-compiler pkg-config
         mamba activate d3

   .. tab-item:: Python
      :sync: python

      .. code-block:: text

         mamba create d3 dftd3-python
         mamba activate d3

You can run the example code with

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: shell

         ❯ $FC energy.f90 $(pkg-config s-dftd3 mctc-lib --cflags --libs) && ./a.out
         Dispersion energy for PBE0-D3(BJ) is -0.0009218696 Hartree

   .. tab-item:: C
      :sync: c

      .. code-block:: shell

         ❯ $CC energy.c $(pkg-config s-dftd3 mctc-lib --cflags --libs) && ./a.out
         Dispersion energy for PBE0-D3(BJ) is -0.0009218696 Hartree

   .. tab-item:: Python
      :sync: python

      .. code-block:: shell

         ❯ python energy.py
         Dispersion energy for PBE0-D3(BJ) is -0.0009219059 Hartree

