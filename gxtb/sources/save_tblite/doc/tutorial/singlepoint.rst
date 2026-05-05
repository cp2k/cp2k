Single point calculations
=========================

Evaluating the tight-binding methods is possible using the `tblite-run`_ command.

.. _tblite-run: https://github.com/awvwgk/tblite/blob/main/man/tblite-run.1.adoc
.. _tblite-param: https://github.com/awvwgk/tblite/blob/main/man/tblite-param.1.adoc
.. _tblite-tag: https://github.com/awvwgk/tblite/blob/main/man/tblite-tag.5.adoc
.. _dftb+: https://dftbplus.org


Calculating with the standalone
-------------------------------

The *tblite* standalone allows to perform single point calculations with the :ref:`built-in methods <methods>` and :ref:`parameter files <parameter>`.


Energy evaluation with different methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run a first calculation we only need a geometry input in one of the supported formats.
As a simple example we choose the caffeine molecule and save it as ``struc.xyz``:

.. code-block:: text
   :caption: struc.xyz

   24
   # caffeine
   C    1.07317    0.04885   -0.07573
   N    2.51365    0.01256   -0.07580
   C    3.35199    1.09592   -0.07533
   N    4.61898    0.73028   -0.07549
   C    4.57907   -0.63144   -0.07531
   C    3.30131   -1.10256   -0.07524
   C    2.98068   -2.48687   -0.07377
   O    1.82530   -2.90038   -0.07577
   N    4.11440   -3.30433   -0.06936
   C    5.45174   -2.85618   -0.07235
   O    6.38934   -3.65965   -0.07232
   N    5.66240   -1.47682   -0.07487
   C    7.00947   -0.93648   -0.07524
   C    3.92063   -4.74093   -0.06158
   H    0.73398    1.08786   -0.07503
   H    0.71239   -0.45698    0.82335
   H    0.71240   -0.45580   -0.97549
   H    2.99301    2.11762   -0.07478
   H    7.76531   -1.72634   -0.07591
   H    7.14864   -0.32182    0.81969
   H    7.14802   -0.32076   -0.96953
   H    2.86501   -5.02316   -0.05833
   H    4.40233   -5.15920    0.82837
   H    4.40017   -5.16929   -0.94780

To perform a single point calculation we use the `tblite-run`_ subcommand:

.. code-block:: text

   tblite run struc.xyz

The output is written to the standard output, no further files are produced at this point.
By default the GFN2-xTB method is used.

.. code-block:: text
   :caption: standard output for GFN2-xTB

    repulsion energy        4.9222938485600E-01 Eh
    dispersion energy       9.6582008183032E-05 Eh
    number of electrons     7.4000000000000E+01 e
    integral cutoff         1.7347787504999E+01 bohr

   ------------------------------------------------------------
     cycle        total energy    energy error   density error
   ------------------------------------------------------------
         1     -41.75569301832  -4.2248019E+01   1.9479952E-01
         2     -42.10609780052  -3.5040478E-01   7.5972187E-02
         3     -42.14069903418  -3.4601234E-02   4.6343369E-02
         4     -42.14507107704  -4.3720429E-03   1.2550655E-02
         5     -42.14724061029  -2.1695333E-03   6.1305110E-03
         6     -42.14770623967  -4.6562938E-04   2.4358040E-03
         7     -42.14738299681   3.2324286E-04   1.0726483E-03
         8     -42.14747017326  -8.7176447E-05   3.6117531E-04
         9     -42.14746157709   8.5961684E-06   1.6698164E-04
        10     -42.14746263709  -1.0599990E-06   5.4777834E-05
        11     -42.14746337487  -7.3778077E-07   2.5306995E-05
        12     -42.14746403618  -6.6130801E-07   8.2333411E-06
   ------------------------------------------------------------

    electronic energy      -4.2639790003039E+01 Eh
    total energy           -4.2147464036175E+01 Eh

   ------------------------------------------------------------------
         #    Occupation            Energy/Eh            Energy/eV
   ------------------------------------------------------------------
         1        2.0000           -0.7201238             -19.5956
       ...           ...                  ...                  ...
        29        2.0000           -0.4792385             -13.0407
        30        2.0000           -0.4764138             -12.9639
        31        2.0000           -0.4564304             -12.4201
        32        2.0000           -0.4379226             -11.9165
        33        2.0000           -0.4264752             -11.6050
        34        2.0000           -0.4193573             -11.4113
        35        2.0000           -0.4148063             -11.2875
        36        2.0000           -0.4035253             -10.9805
        37        2.0000           -0.3892737             -10.5927 (HOMO)
        38                         -0.2657013              -7.2301 (LUMO)
        39                         -0.2265846              -6.1657
        40                         -0.2098114              -5.7093
        41                         -0.1225595              -3.3350
        42                         -0.0559053              -1.5213
        43                         -0.0080016              -0.2177
        44                          0.0133961               0.3645
       ...                                ...                  ...
        66                          0.6592542              17.9392
   ------------------------------------------------------------------
                  HL-Gap            0.1235725 Eh            3.3626 eV
   ------------------------------------------------------------------

We can select a method by adding the ``--method`` argument to the command line and rerun

.. code-block:: text

   tblite run --method gfn1 struc.xyz

This will perform the calculation again using the GFN1-xTB method.

We can add a parametrized solvation model with the ``--alpb`` or ``--gbsa`` flag followed by the solvent name. Alternatively, the purely electrostatic generalized Born model (``--gb``/ ``-gbe``) or the polarizable continuum model (``--cpcm``) can be specified with the solvent name or the dielectric constant. The reference state can be corrected with the ``--solv-state`` flag (options: ``gsolv`` (default), ``bar1mol``, or ``reference``).

.. code-block:: text

   tblite run --method gfn1 --alpb water struc.xyz

Calculating derivatives
~~~~~~~~~~~~~~~~~~~~~~~

So far only energies are calculated, to enable the evaluation of molecular gradients we add the ``--grad`` keyword to the command line for our GFN1-xTB calculation

.. code-block:: text

   tblite run --grad --method gfn1 struc.xyz

.. note::

   The gradient option takes one optional argument, make sure to separate the geometry input file by the ``--`` separator in case you do not provide the output file for the derivatives:

   .. code-block:: text

      tblite run --grad -- struc.xyz

The output will show that the gradient output has been written to a file

.. code-block:: text

   [Info] Tight-binding results written to 'tblite.txt'

.. tip::

   Generally, *tblite* will never write to a file without explicitly stating it in the standard output.

We can inspect the result using a pager like *less* or any text editor of your choice, the output should look similar to this:

.. code-block:: text
   :caption: tblite.txt

   energy             :real:0:
    -4.4509702550094509E+01
   gradient           :real:2:3,24
     3.8751207233845646E-03  1.0589576505981196E-03  4.0101833008853537E-06
    -6.8197644662543388E-03 -2.5469121433516095E-02 -4.7470927704034663E-05
     1.0648728891662693E-02  8.8366335898543858E-03  4.7996917287539728E-05
     3.1902722527341862E-04  7.4242248622642722E-03 -2.8496918953763304E-05
    -2.6208973446986007E-02 -9.0131378921360604E-03 -3.8669871255270171E-05
    -1.7509318247989529E-03  3.6705299785918994E-03  3.4118097773531792E-05
     1.9484148290568090E-02  4.7415207908131207E-03  4.0471496678204742E-05
    -1.1722356880465069E-02  7.4828408528079754E-03 -1.1612153192208481E-04
    -1.2241607563977187E-03 -1.9132992220688570E-02  1.2461589545076379E-04
    -1.6766717498318165E-02  1.2103660828186829E-02 -3.4573645420344302E-05
     2.4844203069064843E-02 -1.9287321306216326E-02 -4.9197479201307158E-05
     1.1084471631753511E-02  1.7004159137354469E-02 -3.2214753172841759E-06
    -6.4592191454479423E-03 -4.2860776915391351E-03 -4.0858618678683595E-06
    -6.5899183170761598E-03  1.3221070732562403E-02 -9.0165446293088262E-05
    -1.3268732296168747E-03  2.1245633439957766E-03  4.0929253563542303E-06
     1.7227511511327858E-03  1.1717306153688867E-03  2.6844186750092148E-03
     1.7254183949399506E-03  1.1710614866087398E-03 -2.6879103256391306E-03
    -1.5499406234334675E-03  3.9313477434577580E-03  2.5224714308249690E-05
     6.3570530436816057E-03  2.2899202031480908E-03  3.4807215781027190E-06
     5.1739877369726952E-04 -5.0616912429584586E-04  1.9423398812888478E-03
     5.4061192657724807E-04 -4.8581923304342556E-04 -1.9455409758072379E-03
    -1.2675930676672822E-03 -7.6587507428546708E-03  8.4916914284767576E-05
     2.7792577845120209E-04 -1.8293010741325709E-04  1.4003154522332375E-03
     2.8959035627490622E-04 -2.0990206390925764E-04 -1.3505474151682963E-03
   virial             :real:2:3,3
     7.8713512660909479E-02 -3.0635463233997452E-02 -1.6134089053178151E-04
    -3.0635463233997452E-02  8.4807095398895124E-02 -1.1909491256157945E-04
    -1.6134089053176757E-04 -1.1909491256159344E-04  2.0323518845924864E-02

The ASCII output format is adopted from the `DFTB+`_ project and called tagged output (see `tblite-tag`_).
It has the advantage of being easily human-readable and a custom parser is straight-forward to write.

Instead of relying on a custom parser, a more standardized format like JSON might be preferable for automating calculations with *tblite*.
We can request to create to use JSON to store the results with the ``--json`` argument.

.. code-block:: text

   tblite run --grad --json --method gfn1 struc.xyz

In the output we will now find

.. code-block:: text

   [Info] JSON dump of results written to 'tblite.json'

.. note::

   The tagged output file will not be written by default when requesting a JSON output, to enable both provide the name of the output files in the command line

   .. code-block:: text

      tblite run --grad caffeine.txt --json caffeine.json struc.xyz

   The output will show that both output are now available

   .. code-block:: text

      [Info] Tight-binding results written to 'caffeine.txt'
      [Info] JSON dump of results written to 'caffeine.json'

We find that *tblite* also included meta data like the version number of the program in the JSON dump, which can be useful when automatically processing the data later.
