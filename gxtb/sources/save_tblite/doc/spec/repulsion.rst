.. _repulsion:

Repulsive contributions
=======================

.. contents::


Effective coulombic repulsion
-----------------------------

The xTB Hamiltonian is using an effective coulomb-like repulsion term

.. math::

   E_\text{rep} =
   \frac12
   \sum_\text{A}^{N_\text{at}}
   \sum_\text{B}^{N_\text{at}}
   \frac{Z_\text{AB}}{R_\text{AB}^{r_\text{AB}}}
   \exp[-\alpha_\text{AB}\cdot R_\text{AB}^{k_\text{AB}}]
   ,

where *R*\ :sub:`AB` is the interatomic distance and *Z*\ :sub:`AB`, *α*\ :sub:`AB`, *r*\ :sub:`AB` and *k*\ :sub:`AB` are pairwise parameters for each species pair.
The pairwise parameters are usually formed from element-specific parameters to limit the actual amount of free parameters.
The effective nuclear charges are formed by

.. math::

   Z_\text{AB} =
   Z_\text{A} \cdot Z_\text{B}
   .

Repulsion exponents are obtained as geometric mean from element specific exponents by

.. math::

   \alpha_\text{AB} =
   \sqrt{\alpha_\text{A}\alpha_\text{B}}
   .

Finally, the exponents, *r*\ :sub:`AB` and *k*\ :sub:`AB`, of the distance dependent contributions are usually set globally for all element-pairs.
