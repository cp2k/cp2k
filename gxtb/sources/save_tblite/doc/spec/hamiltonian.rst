.. _hamiltonian:

Hamiltonian Specification
=========================

.. contents::


Basis set definition
--------------------

The Hamiltonian is defined in a minimal basis set using the linear combination of atomic orbitals (LCAO) approach.
Atomic orbitals are represented by contracted Gaussian basis functions of the STO-NG type basis sets.
This implementation support overlap integrals up to g-functions and one moment lower, f-functions, for derivatives.


Spherical harmonics ordering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Standard spherical harmonic sorting is used for the basis functions.

================== ===============================
 angular momentum   order
================== ===============================
 0 (s)              0
 1 (p)              -1, 0, 1
 2 (d)              -2, -1, 0, 1, 2
 3 (f)              -3, -2, -1, 0, 1, 2, 3
 4 (g)              -4, -3, -2, -1, 0, 1, 2, 3, 4
================== ===============================


Potential definition
--------------------

The potential shift vectors are defined with respect the partial charge and cumulative atomic moments.
This is the negative of the potential with respect to the populations.
The sign is applied when adding the potential shifts to the Hamiltonian.

.. note::

   This convention is in line with the sign convention for partial charges, but is unintuitive from a population / density perspective.
   The sign convention might change to the latter case as the latter convention is more in line with usual density functional theory.


Overlap integrals
-----------------

The overlap integral is defined as

.. math::

   S_{\mu\nu} = \langle \mu | \nu \rangle

Normalisation of the basis function in practise is accurate to 10:sup:`-10`, due to the precision of the tabulated STO-NG contraction coefficients.\ :footcite:`stewart1970`


Atomic partial charges
~~~~~~~~~~~~~~~~~~~~~~

Orbital populations are obtained by Mulliken population analysis

.. math::

   q_\nu = n^0_\nu - \sum_\kappa^{N_\text{ao}} P_{\nu\kappa} S_{\kappa\nu}

and summed up to respective shell-resolved or atomic partial charges.


Dipole moment integrals
-----------------------

The dipole moment integral is defined origin independent as

.. math::

   D_{\mu\nu} = \langle \mu | \vec r - \vec R_\nu | \nu \rangle

The dipole operator is always placed on the aufpunkt of the basis function in ket.


Cumulative dipole moments
~~~~~~~~~~~~~~~~~~~~~~~~~

The dipole moments are evaluated from the dipole moment integrals using Mulliken population analysis

.. math::

   \vec\mu_\nu = -\sum_\kappa^{N_\text{ao}} P_{\nu\kappa} D_{\kappa\nu}

and summed up to respective cumulative shell-resolved or atomic dipole moments.


Quadrupole moment integrals
---------------------------

The quadrupole moment integral is defined origin independent as

.. math::

   Q_{\mu\nu} = \langle \mu | (\vec r - \vec R_\nu) \times (\vec r - \vec R_\nu) | \nu \rangle

The quadrupole operator is always placed on the aufpunkt of the basis function in ket.
Only the lower triangle of the traceless quadrupole integrals are saved packed (ordering: xx, xy, yy, xz, yz, zz).


Cumulative quadrupole moments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The quadrupole moments are evaluated from the quadrupole moment integrals using Mulliken population analysis

.. math::

   \vec\theta_\nu = -\sum_\kappa^{N_\text{ao}} P_{\nu\kappa} Q_{\kappa\nu}

and summed up to respective cumulative shell-resolved or atomic quadrupole moments.


Literature
----------

.. footbibliography::
