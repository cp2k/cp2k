.. _coulomb:

Electrostatic interactions
==========================

.. contents::


Second order
------------


Isotropic electrostatics
~~~~~~~~~~~~~~~~~~~~~~~~

The isotropic electrostatic in a shell-resolved formulation is given by the parametrized Coulomb interaction between shellwise partial charges

.. math::

   E_\text{IES} =
   \frac12 \sum_{\text{A},\text{B}} \sum_{l,l'}^{s,p,d}
   q^{l}_\text{A} \gamma^{ll'}_\text{AB} q^{l'}_\text{B}

The interaction potential is parametrized by a Klopman–Ohno type potential in the xTB Hamiltonian or the γ-functional as used in the DFTB Hamiltonian.

Klopman–Ohno kernel
^^^^^^^^^^^^^^^^^^^

The interaction kernel for the Klopman–Ohno electrostatic is given by

.. math::

   \gamma^{ll'}_\text{AB} =
   \left(
   R_\text{AB}^g + f_\text{av}(\eta_A^l, \eta_B^{l'})^{-g}
   \right)^{-\frac1g}

where η:sub:`A/B` are the chemical hardness parameters of the respective shells and *g* is the exponent to manipulate the potential shape.


γ-functional kernel
^^^^^^^^^^^^^^^^^^^

The interaction kernel for the DFTB γ-functional is derived from the integral of two exponential densities

.. math::

   \begin{split}
   \gamma^{ll'}_\text{AB} =
   \frac1{R_\text{AB}}
   - \exp[-\tau_\text{A}R]
     \left(
     \frac{\tau_\text{B}^4\tau_\text{A}}{2(\tau_\text{A}^2-\tau_\text{B}^2)^2}
     - \frac{\tau_\text{B}^6\tau_\text{A} - 3\tau_\text{B}^4\tau_\text{A}^2}
       {(\tau_\text{A}^2-\tau_\text{B}^2)^3 R_\text{AB}}
     \right)
     \\
   - \exp[-\tau_\text{B}R]
     \left(
     \frac{\tau_\text{A}^4\tau_\text{B}}{2(\tau_\text{B}^2-\tau_\text{A}^2)^2}
     - \frac{\tau_\text{A}^6\tau_\text{B} - 3\tau_\text{A}^4\tau_\text{B}^2}
       {(\tau_\text{B}^2-\tau_\text{A}^2)^3 R_\text{AB}}
     \right)
   \end{split}

where τ:sub:`A/B` are scaled Hubbard parameters of the respective shells and *R* is the distance between the atomic sides.


Anisotropic electrostatics
~~~~~~~~~~~~~~~~~~~~~~~~~~

The anisotropic electrostatic in an atom-resolved formulation is given by the multipole interactions between the different moments:

.. math::

   E_\text{AES} =
   \sum_{\text{A},\text{B}} \sum_{k}^{x,y,z}
   q_\text{A} \gamma^{k}_\text{AB} \mu^{k}_\text{B}
   + \frac12 \sum_{\text{A},\text{B}} \sum_{k,k'}^{x,y,z}
   \mu^{k}_\text{A} \gamma^{kk'}_\text{AB} \mu^{k'}_\text{B}
   + \sum_{\text{A},\text{B}} \sum_{k,k'}^{x,y,z}
   q_\text{A} \gamma^{kk'}_\text{AB} \theta^{kk'}_\text{B}


Third order
-----------

The isotropic third-order contributions are included as the trace of the on-site shell-resolved Hubbard derivatives.

.. math::

   E_\text{IXC} =
   \frac13 \sum_\text{A} \sum_{l}
   \Gamma^l_\text{A} (q^l_\text{A})^3
