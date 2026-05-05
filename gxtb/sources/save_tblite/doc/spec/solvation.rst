.. _solvation:


Solvation
=========

Implicit solvation models are available for the calculation of the solvation free energies partitioned as 

.. math::
   \Delta G_{\text{solv}} = \Delta G_{\text{polar}} + \Delta G_{\text{npol}} + \Delta G_{\text{shift}}

including the polar contribution :math:`{\Delta G_{\text{polar}}}` based on electrostatics, the non-polar contribution :math:`{\Delta G_{\text{npol}}}` based on cavity formation and dispersion, and a constant shift :math:`{\Delta G_{\text{shift}}}` depending on the thermodynamic state of initial gas and final liquid solution.

ALPB solvation model
--------------------

The analytical linearized Poisson-Boltzmann (ALPB) model evaluates the polar contribution

.. math::
   \Delta G^{\text{ALPB}}_{\text{polar}} = 
   - \frac{1}{2} \left(\frac{1}{\epsilon_{\text{in}}} - \frac{1}{\epsilon_{\text{out}}}\right) 
   \frac{1}{1+\alpha\beta}
   \sum_{A,B} q_{A} q_{B} \left( \frac{1}{f(R_{AB, a_{A}, a_{B}})} + \frac{\alpha\beta}{\mathcal{A}_{\text{det}}} \right)

based on the ALPB constant :math:`{\alpha}` (set to 0.571214), the solute (:math:`{\epsilon_{\text{in}}=1}`) and solvent (:math:`{\epsilon_{\text{out}}}`) dielectric constants combined in :math:`{\beta=\frac{\epsilon_{\text{in}}}{\epsilon_{\text{out}}}}`, atomic partial charges :math:`{q_{A/B}}`, and the electrostatic size of the solute :math:`{\mathcal{A}_{\text{det}}}`. \ :footcite:`ehlert2021`
:math:`{f(R_{AB}, a_{A}, a_{B})}` is the interaction kernel with the Born radii :math:`{a_{A/B}}` and can take two forms, either 

.. math::
   f(R_{AB}, a_{A}, a_{B}) = \left( R_{AB}^2 + a_{A} a_{B} \exp\left[-\frac{R_{AB}^2}{4 a_{A} a_{B}} \right] \right)^{\frac{1}{2}}

proposed by Still (default in GBSA), or the more recent P16 kernel (default for ALPB): 

.. math::
   f(R_{AB}, a_{A}, a_{B}) = R_{AB} + \sqrt{a_{A} a_{B}} \left(1+\frac{1.028 R_{AB}}{16 \sqrt{a_{A} a_{B}}} \right)^{16}

For specific polar interactions, an atom-wise hydrogen bonding correction is introduced:

.. math::
   \Delta G^{\text{HB}}_{\text{polar}} = \sum \Delta G^{\text{HB}}_{\text{A}}

In addition to the polar contribution, the non-polar contribution is included with a cavity dispersion solvation term (CDS) based on the atomic surface tension :math:`\gamma_{A}` and the solvent-accessible surface area (SASA) :math:`\sigma_{A}`: 

.. math::
   \Delta G^{\text{CDS}}_{\text{npol}} = \sum_{A} \gamma_{A} \sigma_{A}

An additional empirical constant shift is applied to the solvation free energy.
A solution state correction can be activated but is not included by default.

GBSA solvation model
--------------------

The generalized Born solvation model (GBSA) is a simplified version of ALPB in the limit of an ideal conductor environment (:math:`{\epsilon_{\text{out}}}\rightarrow \infty` and :math:`{\beta\rightarrow 0}`).
As for ALPB, CDS and a constant shift shift are applied, while a solution state correction can be activated (only if the solvent is specified by name).

CPCM solvation model 
--------------------

The conductor-like polarizable continuum solvation model is implemented based on the domain-decomposition approach and is currently available only for the polar part :math:`{\Delta G_{\text{polar}}}`.

Solution state correction
-------------------------

For solvation free energies, the state of the inital gas and final liquid solution can be changed with a solution state correction.
By default no solution state correction is applied (gsolv, default), which is comparable with most other solvation models (SMD, COSMO-RS, ...).
For normal production runs, the option ``bar1mol`` should be used. For explicit comparisons with ``reference`` state corrected COSMO-RS, the ``reference`` option should be used (includes solvent-specific correction for infinite dilution).
Solution state correction is available for the ALPB and GBSA solvation models.

================== ====================================================================
 Name               Definition
================== ====================================================================
 gsolv (default)    1 L of ideal gas and 1 L of liquid solution 
 bar1mol            1 bar of ideal gas and 1 mol/L liquid solution 
 reference          1 bar of ideal gas and 1 mol/L liquid solution at infinite dilution
================== ====================================================================


Literature
----------

.. footbibliography::
