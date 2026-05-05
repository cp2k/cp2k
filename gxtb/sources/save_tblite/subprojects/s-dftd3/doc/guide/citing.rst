Suggesting references for D3 dispersion corrections
---------------------------------------------------

The D3 dispersion correction was originally introduced as D3(0)\ :footcite:`grimme2010` and later extended to D3(BJ).\ :footcite:`grimme2011`
Additional damping functions and parametrizations for new functionals have been introduced in the literature over time.
Providing proper attribution to the method, parameters and software can be challenging to keep up with.
The *s-dftd3* library keeps track of the references for each method and parameter set, and provides an option to suggest the used references.

For example when running a computation with SCAN-D3(BJ), we can use the following command to get the suggested references:

.. code-block:: text

   s-dftd3 structure.xyz --bj scan --citation

In the output we will see

.. code-block:: text

   -----------------------------------
    s i m p l e   D F T - D 3  v1.2.1
   -----------------------------------

   Rational (Becke-Johnson) damping: scan-D3(BJ)
   ---------------------
     s6         1.0000
     s8         0.0000
     s9         0.0000
     a1         0.5380
     a2         5.4200
    alp        14.0000
   --------------------

   Dispersion energy:      -1.0554287327920E-02 Eh

   [Info] Writing Dispersion energy to '.EDISP'
   [Info] Citation information written to 'dftd3.bib'

Inspecting this file shows three suggested references, here first the citation for the library itself used to perform the calculation, second the original publication of the D3(BJ) method and finally the publication introducing parameters for the SCAN functional.

.. code-block:: bib
   :caption: dftd3.bib

   @article{10.21105/joss.07169,
     title = {{Simple DFT-D3: Library first implementation of the D3 dispersion correction}},
     author = {Sebastian Ehlert},
     issue = {103},
     volume = {9},
     pages = {7169},
     doi = {10.21105/joss.07169},
     url = {https://doi.org/10.21105/joss.07169}
   }
   @article{10.1002/jcc.21759,
     title = {{Effect of the damping function in dispersion corrected density functional theory}},
     author = {Stefan Grimme
       and Stephan Ehrlich
       and Lars Goerigk},
     issue = {7},
     volume = {32},
     pages = {1456-1465},
     doi = {10.1002/jcc.21759},
     url = {https://doi.org/10.1002/jcc.21759}
   }
   @article{10.1103/physrevb.94.115144,
     title = {{Benchmark tests of a strongly constrained semilocal functional with a  long-range dispersion correction}},
     author = {J. G. Brandenburg
       and J. E. Bates
       and J. Sun
       and J. P. Perdew},
     volume = {94},
     pages = {115144},
     doi = {10.1103/physrevb.94.115144},
     url = {https://doi.org/10.1103/physrevb.94.115144}
   }

In the computational details section in a publication using SCAN-D3(BJ) we can now refer to these publications for example as

.. code-block:: tex

   The D3(BJ) dispersion correction\cite{10.1002/jcc.21759} for
   SCAN-D3(BJ)\cite{10.1002/jcc.21759} was calculated
   using the s-dftd3 library (version 1.2.1).\cite{10.21105/joss.07169}

.. important::

   The suggested references are based on the information available in the library.
   Always check the original publications for the most recent references and the correct citation.
