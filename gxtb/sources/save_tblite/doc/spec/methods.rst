.. _methods:

Built-in methods
================

The library supports GFN1-xTB\ :footcite:`grimme2017`, GFN2-xTB\ :footcite:`bannwarth2019`, and IPEA1-xTB\ :footcite:`asgeirsson2017` as built-in parametrizations.
The parameters can be loaded in the standalone using the ``--method`` keyword and can be loaded directly from the API without requiring to create a parametrization first.

=============== ===========
 method          keyword
=============== ===========
 GFN1-xTB        *gfn1*
 GFN2-xTB        *gfn2*
 IPEA1-xTB       *ipea1*
=============== ===========

The built-in parameters can be exported using the `tblite-param`_ command.

.. code-block:: none

   tblite param --method gfn2 --output gfn2-xtb.toml

For more information on the format of the parameter file see the :ref:`parameter specification <parameter>`.

.. _tblite-param: https://github.com/awvwgk/tblite/blob/main/man/tblite-param.1.adoc
