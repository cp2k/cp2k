# Constrained molecular dynamics

CP2K can constrain a collective variable (CV) during molecular dynamics with
[MOTION/CONSTRAINT/COLLECTIVE](#CP2K_INPUT.MOTION.CONSTRAINT.COLLECTIVE). The CV is defined in
[FORCE_EVAL/SUBSYS/COLVAR](#CP2K_INPUT.FORCE_EVAL.SUBSYS.COLVAR), and `COLVAR` in the constraint
section selects it by input order. Use
[INTERMOLECULAR](#CP2K_INPUT.MOTION.CONSTRAINT.COLLECTIVE.INTERMOLECULAR) for a CV whose atoms are
not all in the same molecular object.

A fixed-distance window between two atoms can be set up as follows:

```none
&FORCE_EVAL
  ...
  &SUBSYS
    ...
    &COLVAR
      &DISTANCE
        ATOMS 1 2
      &END DISTANCE
    &END COLVAR
  &END SUBSYS
&END FORCE_EVAL

&MOTION
  &CONSTRAINT
    &COLLECTIVE
      COLVAR 1
      INTERMOLECULAR TRUE
      TARGET [angstrom] 2.0
    &END COLLECTIVE
    &LAGRANGE_MULTIPLIERS ON
      FILENAME constraint_force
      COMMON_ITERATION_LEVELS 1
    &END LAGRANGE_MULTIPLIERS
  &END CONSTRAINT
  &MD
    ...
  &END MD
&END MOTION
```

## Lagrange-multiplier output

[LAGRANGE_MULTIPLIERS](#CP2K_INPUT.MOTION.CONSTRAINT.LAGRANGE_MULTIPLIERS) writes two records during
a velocity-Verlet step:

- `Shake Lagrangian Multipliers` contains the position-constraint multipliers. These are the values
  relevant to a configurational constraint force and constrained thermodynamic integration (TI).
- `Rattle Lagrangian Multipliers` contains the velocity-constraint multipliers. They enforce the
  time derivative of the constraint and are not a second configurational-force sample.

The values are raw CP2K internal quantities; specifying `TARGET` in another unit does not convert
the printed multipliers. A distance constraint therefore produces a multiplier in hartree/bohr,
while an angular constraint uses the corresponding internal angular unit. The file contains the
multipliers of all active constraints, first intramolecular constraints and then intermolecular
constraints. The ordering within each group is collective-variable, 3-by-3, and 4-by-6 constraints.

## Blue-moon ensemble correction

Constrained molecular dynamics generally produces biased statistical distributions. The blue-moon
ensemble average is used to correct these biased outputs and retrieve the statistical properties
corresponding to unconstrained molecular-dynamics conditions. The printed SHAKE multiplier is not,
in general, a complete blue-moon estimator. For one constrained reaction coordinate $\xi$, the
free-energy gradient then follows

$$
\frac{\mathrm d A}{\mathrm d\xi} 
= 
\frac{\left\langle Z^{-1/2}
\left(-\lambda + k_\mathrm{B} T G\right)
\right\rangle_\xi}
{\left\langle Z^{-1/2}\right\rangle_\xi}
$$

where $k_\mathrm{B}$ is the Boltzmann constant, $T$ is the temperature,
$\left\langle \cdots \right\rangle_\xi$ denotes the time average of reaction coordinate
$\xi(\mathbf r_1,\dots,\mathbf r_N)$ by MD simulation.

Here $Z$ defines the scalar mass metric

$$
Z = \sum_i \frac{1}{m_i}
\left|\nabla_i\xi\right|^2.
$$

With the CP2K convention that the constraint force is $-\lambda\nabla\xi$, the free-energy gradient
contains a $Z^{-1/2}$ reweighting and, for a general coordinate, an additional metric-derivative
term ($G$). With Cartesian gradient $\mathbf g=\nabla\xi$, Hessian $H=\nabla\nabla\xi$, and diagonal
Cartesian mass matrix $M$, this term is

$$
G = \frac{(M^{-1}\mathbf g)^T H (M^{-1}\mathbf g)}{Z^2}
  = \frac{(M^{-1}\mathbf g)\mathbin{\cdot}\nabla Z}{2Z^2}.
$$

CP2K currently writes $\lambda$ but does not evaluate or print the complete corrected blue-moon
estimator. The required metric terms therefore have to be evaluated during postprocessing for the
chosen reaction coordinate. The standalone
[blue-moon postprocessor](https://github.com/cp2k/cp2k/tree/master/tools/blue_moon) evaluates both
terms for **one fixed collective constraint**, including distances, angles, torsions, signed
point-plane distances, point-to-bond-center distances and their linear combinations. It requires the
matching XYZ trajectory, SHAKE output, actual atomic masses, fixed cell, temperature and exact CV
definition. Its README describes supported conventions, data alignment and limitations; it is not a
general multiple-constraint estimator. See
[Komeiji, Chem-Bio Informatics Journal 7, 12 (2007)](https://doi.org/10.1273/cbij.7.12) for the
general expression and explicit algorithms for two common coordinates.

For the distance between two atoms,

$$
\xi = |\mathbf r_i-\mathbf r_j|,
\qquad
Z = m_i^{-1}+m_j^{-1}.
$$

Here $Z$ is constant and the metric-derivative term is zero. The reweighting cancels, so the
free-energy gradient reduces to $-\langle\lambda\rangle_\xi$ with the sign convention above.

For the three-atom distance difference

$$
\xi = |\mathbf r_i-\mathbf r_j|-|\mathbf r_k-\mathbf r_j|,
$$

the metric is

$$
Z = m_i^{-1}+m_k^{-1}
    +2m_j^{-1}
    \left(1-\boldsymbol\rho_{ij}\mathbin{\cdot}\boldsymbol\rho_{kj}
    \right)
$$

and depends on the instantaneous angle. Here $r_{ij}=|\mathbf r_i-\mathbf r_j|$,
$r_{kj}=|\mathbf r_k-\mathbf r_j|$, $\boldsymbol\rho_{ij}=(\mathbf r_i-\mathbf r_j)/r_{ij}$ and
similarly for $\boldsymbol\rho_{kj}$. Applying the general Hessian expression above gives

$$
G = \frac{1-(\boldsymbol\rho_{ij}\mathbin{\cdot}\boldsymbol\rho_{kj})^2}
         {m_j^2 Z^2}
    \left(\frac{1}{r_{ij}}-\frac{1}{r_{kj}}\right).
$$

Thus $G$ vanishes for equal distances or collinear bonds, but not for a general configuration. Both
reweighting and the $k_\mathrm{B}TG$ term must otherwise be retained. In particular, the $G=0$
simplification stated in Eq. (35) of Komeiji does not follow from the general Eq. (8) for unequal,
noncollinear bonds; the postprocessor tests the general expression against independent analytic and
finite-difference evaluations.

Simplifications for a particular coordinate must not be assumed for an arbitrary
[COMBINE_COLVAR](#CP2K_INPUT.FORCE_EVAL.SUBSYS.COLVAR.COMBINE_COLVAR), coordination number, or
multiple simultaneous constraints. Additional constraints, including fixed atoms or rigid molecules,
change the accessible degrees of freedom and require an appropriate estimator.

## Fixed windows and moving constraints

For equilibrium constrained thermodynamic integration, run independently equilibrated trajectories
at a series of fixed `TARGET` values, calculate the corrected free-energy gradient in every window,
and integrate it over the reaction coordinate. Check the sampling length, correlation time, window
spacing, integration direction, and unit conversion.

[TARGET_GROWTH](#CP2K_INPUT.MOTION.CONSTRAINT.COLLECTIVE.TARGET_GROWTH) instead changes `TARGET`
linearly by `TARGET_GROWTH * TIMESTEP` at every MD step, optionally stopping at
[TARGET_LIMIT](#CP2K_INPUT.MOTION.CONSTRAINT.COLLECTIVE.TARGET_LIMIT). This is a moving-constraint
or slow-growth protocol. A simple example for the performance of CV as the distance between two
atoms:

```none
&FORCE_EVAL
  ...
  &SUBSYS
    ...
    &COLVAR
      &DISTANCE
        ATOMS 1 2
      &END DISTANCE
    &END COLVAR
  &END SUBSYS
&END FORCE_EVAL

&MOTION
  &CONSTRAINT
    &COLLECTIVE
      COLVAR 1
      INTERMOLECULAR TRUE
      TARGET [angstrom] 2.0
      TARGET_GROWTH [angstrom*fs^-1] 0.0008
      TARGET_LIMIT [angstrom] 3.0
    &END COLLECTIVE
    &LAGRANGE_MULTIPLIERS ON
      FILENAME constraint_force
      COMMON_ITERATION_LEVELS 1
    &END LAGRANGE_MULTIPLIERS
  &END CONSTRAINT
  &MD
    ...
  &END MD
&END MOTION
```

CP2K does not integrate the work or turn the resulting trajectory into an equilibrium free-energy
profile automatically. A finite pulling rate can cause lag, dissipation, and direction-dependent
hysteresis, so such a trajectory must be analysed with a method appropriate to the intended
nonequilibrium protocol. In the limit of infinitesimal change of $\xi$, the irreversible work ($W$)
describes the energy change between the initial and final state. The resulting work from slow growth
(typically the irreversible work) can be related to the free-energy change $(\Delta A)$ via
Jarzynski's identity

$$
\exp\left(-\frac{\Delta A}{k_\mathrm{B} T}\right) 
= 
\left\langle \exp\left(-\frac{W}{k_\mathrm{B} T}\right) 
\right\rangle
$$
