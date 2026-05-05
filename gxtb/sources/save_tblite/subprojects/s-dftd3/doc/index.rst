Reimplementation of D3 dispersion correction
============================================

This pages describe the library first implementation of the DFT-D3 dispersion correction, in the `s-dftd3`_ software package.\ :footcite:`ehlert2024`.
This project aims to provide a user-friendly and uniform interface to the D3 dispersion model and for the calculation of DFT-D3 dispersion corrections.
Additionally, it provides the geometric counter-poise correction to create composite electronic structure methods of the 3c family.

Supported features include:

- Rational damping function, D3(BJ).\ :footcite:`grimme2011`
- Zero damping function D3(0).\ :footcite:`grimme2010`
- Modified rational and zero damping functions, D3M(BJ) and D3M(0).\ :footcite:`smith2016`
- Optimized power damping function, D3(op).\ :footcite:`witte2017`
- CSO (C6-scaled only) damping function, D3(CSO).\ :footcite:`schroeder2015`
- Axilrod-Teller-Muto three-center contribution.
- Pairwise analysis of dispersion contributions.
- Extensive parameter support for (almost) all published D3 parameters.
- Geometric counter-poise correction and short-range bond correction.\ :footcite:`kruse2012`
- Readily available in Fortran (:ref:`dftd3 module <fortran-api>`),
  C (:ref:`dftd3.h header <c-api>`),
  Python (:ref:`dftd3 package <python>`),
  and via command line (`s-dftd3 executable <https://github.com/dftd3/simple-dftd3/blob/main/man/s-dftd3.1.adoc>`__)


.. footbibliography::


.. _s-dftd3: https://github.com/dftd3/simple-dftd3

.. toctree::
   :maxdepth: 2

   Installation <installation>
   Tutorial <tutorial/index>
   How-to <guide/index>
   Comparison <comparison>
   API <api/index>
