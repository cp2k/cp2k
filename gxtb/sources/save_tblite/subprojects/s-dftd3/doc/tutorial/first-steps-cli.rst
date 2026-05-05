Computing interaction energies
==============================

In this tutorial we want to compute the interaction energy of a non-covalently bound complex using DFT with D3 dispersion correction.
For this purpose we are using the forth system from the S66x8 data set.\ :footcite:`golokesh2022`

We start by providing the coordinates of the dimer system in xyz format.

.. code-block:: text
   :caption: dimer.xyz

   15
   water-peptide dimer, PBE0/def2-QZVP energy: -324.751193159385
   O    -3.2939688    0.4402024    0.1621802 
   H    -3.8134112    1.2387332    0.2637577 
   H    -2.3770466    0.7564365    0.1766203 
   C    -0.6611637   -1.4159110   -0.1449409 
   H    -0.0112009   -2.2770229   -0.2778563 
   H    -1.3421397   -1.3384389   -0.9888061 
   H    -1.2741806   -1.5547070    0.7420675 
   C     0.0935684   -0.1178981   -0.0123474 
   O    -0.4831471    0.9573968    0.1442414 
   N     1.4442015   -0.2154008   -0.0769653 
   H     1.8451531   -1.1259348   -0.2064804 
   C     2.3124436    0.9365697    0.0324778 
   H     1.6759495    1.8048701    0.1672624 
   H     2.9780331    0.8451145    0.8885706 
   H     2.9069093    1.0659902   -0.8697814 

Since D3 is a post SCF correction, we need to compute the SCF interaction energy with another program first.
For this system the PBE0/def2-QZVP interaction energy is computed with Orca 5.0.4 and already included in the comment line of the xyz file.

Our next step is to compute the D3 contribution to the PBE0-D3 energy, here we choose the rational damping version of D3\ :footcite:`grimme2011` denoted as Becke--Johnson (BJ) damping.
For this we call the ``s-dftd3`` command with the xyz file, here ``dimer.xyz`` and the option ``--bj`` with the functional name, here ``PBE0``.
The output should look like:

.. code-block:: text

    ❯ s-dftd3 dimer.xyz --bj PBE0
    -----------------------------------
    s i m p l e   D F T - D 3  v1.0.0
    -----------------------------------

    Rational (Becke-Johnson) damping: PBE0-D3(BJ)
    ---------------------
    s6         1.0000
    s8         1.2177
    s9         0.0000
    a1         0.4145
    a2         4.8593
    alp        14.0000
    --------------------

    Dispersion energy:      -8.2944752821052E-03 Eh

    [Info] Dispersion energy written to .EDISP

The output reports a few information for us, first it reports the program version we are using, second the parameters used in the dispersion calculation, we find the rational damping reported here as *PBE0-D3(BJ)*.
Finally the dispersion energy is writtem to the screen in the output but also stored in the ``.EDISP`` file for easily accessing it without having to search the output.

To compute the interaction energy for the system, we also need to compute the individual monomers.
Since we are computing an interaction energy, the monomers are in the same geometry as in the dimer, i.e. there is no geometry change.
For our system we copy the first three atoms and the last twelve atoms to separate files.

The water monomer is given by

.. code-block:: text
   :caption: monomer-a.xyz

   3
   water monomer, PBE0/def2-QZVP energy: -76.386381675761
   O    -3.2939688    0.4402024    0.1621802 
   H    -3.8134112    1.2387332    0.2637577 
   H    -2.3770466    0.7564365    0.1766203 

The peptide monomer is given by

.. code-block:: text
   :caption: monomer-b.xyz

   12
   peptide monomer, PBE0/def2-QZVP energy: -248.352741874853
   C    -0.6611637   -1.4159110   -0.1449409 
   H    -0.0112009   -2.2770229   -0.2778563 
   H    -1.3421397   -1.3384389   -0.9888061 
   H    -1.2741806   -1.5547070    0.7420675 
   C     0.0935684   -0.1178981   -0.0123474 
   O    -0.4831471    0.9573968    0.1442414 
   N     1.4442015   -0.2154008   -0.0769653 
   H     1.8451531   -1.1259348   -0.2064804 
   C     2.3124436    0.9365697    0.0324778 
   H     1.6759495    1.8048701    0.1672624 
   H     2.9780331    0.8451145    0.8885706 
   H     2.9069093    1.0659902   -0.8697814 

Again we report the PBE0/def2-QZVP total energies computed with Orca 5.0.4 in the comment lines of the xyz files.
For the monomers we now run again the ``s-dftd3`` program keeping the ``--bj`` option with the ``PBE0`` functional name, but changing the xyz file to ``monomer-a.xyz`` with the coordinates of the water monomer and afterwards running the command for the peptide monomer stored in ``monomer-b.xyz``.

.. code-block:: text

   ❯ s-dftd3 monomer-a.xyz --bj pbe0
   -----------------------------------
   s i m p l e   D F T - D 3  v1.0.0
   -----------------------------------

   Rational (Becke-Johnson) damping: pbe0-D3(BJ)
   ---------------------
   s6         1.0000
   s8         1.2177
   s9         0.0000
   a1         0.4145
   a2         4.8593
   alp        14.0000
   --------------------

   Dispersion energy:      -2.7688584084890E-04 Eh

   [Info] Dispersion energy written to .EDISP
   ❯ s-dftd3 monomer-b.xyz --bj pbe0
   -----------------------------------
   s i m p l e   D F T - D 3  v1.0.0
   -----------------------------------

   Rational (Becke-Johnson) damping: pbe0-D3(BJ)
   ---------------------
   s6         1.0000
   s8         1.2177
   s9         0.0000
   a1         0.4145
   a2         4.8593
   alp        14.0000
   --------------------

   Dispersion energy:      -6.5799988724592E-03 Eh

   [Info] Dispersion energy written to .EDISP

We can verify that the damping parameters in the output are identical with the ones before and we can record the total energies from the output or the ``.EDISP`` file.
Finally, we collect the calculation results in a table and compute the interaction energy:

============= ======================== ====================== ===============================
 system        E(PBE0/def2-QZVP) [Eh]   E(D3(BJ)) [Eh]         E(PBE0-D3(BJ)/def2-QZVP) [Eh]
============= ======================== ====================== ===============================
 dimer         -324.751193159385        -0.0082944752821052    -324.75948763466715
 monomer a      -76.386381675762        -0.0002768858408489     -76.38665856160284
 monomer b     -248.352741874853        -0.0065799988724592    -248.35932187372546
 interaction     -0.012069608770        -0.0014375905687970      -0.01350719933885
============= ======================== ====================== ===============================

The total interaction energy computed with PBE0-D3(BJ)/def2-QZVP is therefore -0.012507 Hartree, or -35.5 kJ/mol, from this interaction energy -3.8 kJ/mol or 10% are from the dispersion correction.

The ``s-dftd3`` binary provides access to other versions of D3 as well.
One contribution which is important to consider in dispersion is the non-additivity in the pairwise dispersion energy, in D3 this is possible by including three-body contributions from the Axilrod-Teller-Muto term using the ``--atm`` option.

We reevaluate our dispersion energy again using the D3(BJ)-ATM model this time:

.. code-block:: text

   ❯ s-dftd3 dimer.xyz --bj PBE0 --atm
   -----------------------------------
    s i m p l e   D F T - D 3  v1.0.0
   -----------------------------------

   Rational (Becke-Johnson) damping: PBE0-D3(BJ)-ATM
   ---------------------
     s6         1.0000
     s8         1.2177
     s9         1.0000
     a1         0.4145
     a2         4.8593
    alp        14.0000
   --------------------

   Dispersion energy:      -8.2920550838084E-03 Eh

   [Info] Dispersion energy written to .EDISP
   ❯ s-dftd3 monomer-a.xyz --bj PBE0 --atm
   -----------------------------------
    s i m p l e   D F T - D 3  v1.0.0
   -----------------------------------

   Rational (Becke-Johnson) damping: PBE0-D3(BJ)-ATM
   ---------------------
     s6         1.0000
     s8         1.2177
     s9         1.0000
     a1         0.4145
     a2         4.8593
    alp        14.0000
   --------------------

   Dispersion energy:      -2.7688559086715E-04 Eh

   [Info] Dispersion energy written to .EDISP
   ❯ s-dftd3 monomer-b.xyz --bj --atm
   PBE0 --atm
   -----------------------------------
    s i m p l e   D F T - D 3  v1.0.0
   -----------------------------------

   Rational (Becke-Johnson) damping: PBE0-D3(BJ)-ATM
   ---------------------
     s6         1.0000
     s8         1.2177
     s9         1.0000
     a1         0.4145
     a2         4.8593
    alp        14.0000
   --------------------

   Dispersion energy:      -6.5714439103566E-03 Eh

   [Info] Dispersion energy written to .EDISP

We see for the parameters in the output change, the *s9* value is now reported as one instead of zero as before.
Let's check the impact of the additional contribution for this system by collecting all calculations

============= ======================== ====================== ===================================
 system        E(PBE0/def2-QZVP) [Eh]   E(D3(BJ)-ATM) [Eh]     E(PBE0-D3(BJ)-ATM/def2-QZVP) [Eh]
============= ======================== ====================== ===================================
 dimer         -324.751193159385        -0.0082920550838084     -324.75948521446884
 monomer a      -76.386381675762        -0.0002768855908672      -76.38665856135286
 monomer b     -248.352741874853        -0.0065714439103566     -248.35931331876336
 interaction     -0.012069608770        -0.0014437255825847       -0.01351333435263
============= ======================== ====================== ===================================

In this case, the interaction energy is still -35.5 kJ/mol, as our system is rather small the impact of the non-additivity here is not large yet.
Since it is little effort to include this additional effect with D3 it is useful to include it by default.

In summary we learned in this session on how to *compute dispersion corrections for DFT* calculations with ``s-dftd3``, how to use the *rational damping scheme* for D3, and how to *include non-additive contributions* in the computed dispersion energies.

Literature
----------

.. footbibliography::
