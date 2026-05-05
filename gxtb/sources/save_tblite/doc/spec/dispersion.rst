.. _dispersion:

London dispersion interactions
==============================

.. contents::


D3 dispersion model
-------------------

The D3 dispersion model is provided by the `s-dftd3 <https://github.com/awvwgk/simple-dftd3>`_ library.
The DFT-D3 dispersion correction\ :footcite:`grimme2010` is used together with the rational damping function.\ :footcite:`grimme2011`


D4 dispersion model
-------------------

The D4 dispersion model is provided by the `dftd4 <https://github.com/dftd4/dftd4>`_ project and extended in this project as a self-consistent dispersion model.\ :footcite:`caldeweyher2019` The pairwise smoothed D4S interpolation between reference systems is also available.\ :footcite:`tkachenko2024`

The self-consistent D4 dispersion model can be classified as infinite order charge-dependent contribution in the density expansion, due to the non-linear dependency on the atomic partial charges in the charge scaling function.


Literature
----------

.. footbibliography::
