# Real-Time Bethe-Salpeter Propagation

Instead of solving the Casida equation in the linear response regime, an explicit time-integration
of the equation of motion of electrons can be carried out to determine the excitation frequencies.
In the real-time Bethe-Salpeter propagation (RTBSE) method, the equation of motion is the von
Neumann equation for the single particle density matrix $\hat{\rho}$ with an effective Hamiltonian
$\hat{H}$

$$ \frac{\mathrm{d} \hat{\rho}}{\mathrm{d} t} = -\mathrm{i} [\hat{H}(t), \hat{\rho}(t)] $$

The accuracy of such method is mainly determined by the choice of interaction model in the effective
Hamiltonian. Instead of using TDDFT functionals, the COHSEX approximation to the self-energy is
employed to calculate the time dependent behaviour of the density matrix
\[[Attaccalite2011](http://dx.doi.org/10.1103/PhysRevB.84.245110)\]. This requires a previous
determination of the screened Coulomb potential, done via the bandstructure
[GW](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW) calculation.

The equation of motion is solved in steps

$$ \hat{\rho} (t + \Delta t) = \mathrm{e} ^ {- i \hat{H} (t+\Delta t) \Delta t/2} \mathrm{e} ^ {-i \hat{H}(t) \Delta t/2}
\hat{\rho} (t) \mathrm{e} ^ {i \hat{H}(t) \Delta t/2} \mathrm{e} ^ {i \hat{H} (t + \Delta t) \Delta t/2}$$

which is called the _enforced time reversal
scheme_\[[Castro2004](https://doi.org/10.1063/1.1774980)\]. The effective Hamiltonian is given as

$$ \hat{H}(t) = \hat{h}^{G0W0} + \hat{U} (t) +
\hat{V}^{\mathrm{Hartree}} [\hat{\rho}(t)] - \hat{V}^{\mathrm{Hartree}} [\hat{\rho}_0] +
\hat{\Sigma}^{\mathrm{COHSEX}}[\hat{\rho}(t)] - \hat{\Sigma}^{\mathrm{COHSEX}}[\hat{\rho}_0]
$$

where $\hat{\rho}_0$ is the density matrix determined from the molecular orbitals used in
[GW](#CP2K_INPUT.FORCE_EVAL.PROPERTIES.BANDSTRUCTURE.GW) and $\hat{U}(t)$ is the external applied
field.

## Excitation scheme

Without the external field $\hat{U}(t)$, the density matrix only rotates in phase but does not
produce any measurable dynamics. The excitation of the dynamics can be done either by a real time
pulse (i.e. at each point, $\hat{U}(t)$ follows form due to some finite time dependent field
$\vec{E}(t)$) or by an infinitely sharp delta pulse, which we can understand as the limit of
$\vec{E}(t) \to I \vec{e} \delta(t)$, where $I$ is the delta pulse intensity and $\vec{e}$ its
direction.

## Observables

The dynamics can be traced through time with electric dipole moment associated with the density
matrix

$$ \mu_i(t) = \mathrm{Tr} (\hat{\rho}(t) (\hat{x}_i - x_{i,\mathrm{CC}})) \: ,
$$

where $x_{i,\mathrm{CC}}$ are the coordinates of center of molecular charge and $\hat{x}_i$ is the
position operator.

The electric polarizability (which is related to the photon absorption spectrum) is then determined
as

$$ \alpha_{ij} (\omega) = \frac{\mu_i(\omega)}{E_j(\omega)}
$$

where we Fourier transformed to the frequency domain. In order to stabilise the Fourier transform of
infinitely oscillating dipole moments, we introduce a damping factor
$\gamma$\[[Müller2020](https://doi.org/10.1002/jcc.26412)\]

$$ \mu_i(\omega) = \int _ 0 ^ T \mathrm{d}t \mathrm{e}^{-\gamma t} \mathrm{e} ^ {i \omega t} \mu_i(t) =
    \int _ 0 ^ T \mathrm{d}t \mathrm{e}^{i(\omega + i \gamma) t} \mu_i (t)
$$

One can easily verify that for real FT of the applied field, this leads Lorentzian peaks at the
frequencies of the oscillations of the moments present in the imaginary part of the corresponding
polarizability element.

## Running the Propagation

To run the RTBSE propagation, include the
[RTBSE](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.RTBSE) section in the input file.
[SECTION_PARAMETERS](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.RTBSE.SECTION_PARAMETERS) can
be used to choose TDDFT method as a debug method, but including the section alone leads to a RTBSE
run.

Furthermore, the [TIMESTEP](#CP2K_INPUT.MOTION.MD.TIMESTEP) and [STEPS](#CP2K_INPUT.MOTION.MD.STEPS)
influence the size of each timestep and the total time of propagation. From the properties of the
Fourier transform, one can determine that smaller [TIMESTEP](#CP2K_INPUT.MOTION.MD.TIMESTEP)
increases the maximum energy $\omega$ that is captured by the transform, while larger total
propagation time (influenced by [STEPS](#CP2K_INPUT.MOTION.MD.STEPS)) leads to a better energy
resolution (smaller $\Delta \omega$).

For gas phase/isotropic calculation of polarizability, one needs to run 3 calculations to determine
the trace of the polarizability tensor.

### ETRS Precision

The precision of the self-consistency in the ETRS loop is controlled by the
[EPS_ITER](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.EPS_ITER) keyword. Smaller threshold
(larger precision) lead to more stable propagation, but might require smaller timestep/more
self-consistent iterations.

[MAX_ITER](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.MAX_ITER) keyword is used to determine
the maximum number of self-consistent iterations for a single time step before the cycle is broken
and non-convergence is reported.

If the propagation is converging poorly (>50 ETRS iterations), smaller
[TIMESTEP](#CP2K_INPUT.MOTION.MD.TIMESTEP) may stabilize the propagation. A typical setup prints a
status after each successful ETRS iteration, similar to the following

```
 RTBSE| Simulation step         Convergence     Electron number  ETRS Iterations
 RTBSE|               0     0.55891101E-008     0.16000000E+002                5
 RTBSE| Simulation step         Convergence     Electron number  ETRS Iterations
 RTBSE|               1     0.31847656E-008     0.16000000E+002                5
 RTBSE| Simulation step         Convergence     Electron number  ETRS Iterations
 RTBSE|               2     0.38793291E-008     0.16000000E+002                5
```

### Exponential Method

The method used for the exponentiation of the Hamiltonian is set in the
[MAT_EXP](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.MAT_EXP) keyword, with the following
methods implemented for both TDDFT and RTBSE

- `BCH` - calculates the effect of matrix exponential by series of commutators using
  Baker-Campbell-Hausdorff expansion

and the following methods implemented only for RTBSE

- `EXACT` - Diagonalizes the instantaneous Hamiltonian to determine the exponential exactly

For inexact methods, a threshold for the cutoff of exponential series is provided by the
[EXP_ACCURACY](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.EXP_ACCURACY) keyword. For these,
the [MAX_ITER](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.MAX_ITER) keyword also sets the
maximum number of iterations before the program is stopped and non-convergence is reported.

### Excitation Method

The real time pulse can be specified in the [EFIELD](#CP2K_INPUT.FORCE_EVAL.DFT.EFIELD) section.

If delta pulse is required instead, use
[APPLY_DELTA_PULSE](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.APPLY_DELTA_PULSE) with
[DELTA_PULSE_DIRECTION](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.DELTA_PULSE_DIRECTION) used
for defining the $\vec{e}$ vector and
[DELTA_PULSE_SCALE](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.DELTA_PULSE_SCALE) setting the
$I$ scale of the delta pulse (in atomic units). Note that the definition of the vector is different
from the definition used in the TDDFT method.

The actual value of $I \vec{e}$ is printed out in atomic units, as well as the absolute value of the
maximum element difference between the density matrix before and after the application of the delta
pulse - so called metric difference after delta kick.

```
 RTBSE| Applying delta puls
 RTBSE| Delta pulse elements (a.u.) :   -0.1000E-003  -0.0000E+000  -0.0000E+000
 RTBSE| Metric difference after delta kick                       0.61399576E-004
```

If this metric difference is approaching 1.0, the ETRS cycle might have trouble converging - we
recommend reducing the
[DELTA_PULSE_SCALE](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.DELTA_PULSE_SCALE).

### Printing observables

The code is so far optimised for printing the polarizability elements, which are linked to the
absorption spectrum. The printing of all available properties is controlled in the
[PRINT](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT) section - the ones relevant for
RTBSE propagation are listed here

- [DENSITY_MATRIX](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.DENSITY_MATRIX) - Prints
  the elements of the density matrix in the MO basis into a file at every timestep
- [FIELD](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.FIELD) - Prints the elements of the
  electric field applied at every time step
- [MOMENTS](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.MOMENTS) - Prints the electric
  dipole moment elements reported at every time step
- [MOMENTS_FT](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.MOMENTS_FT) - Prints the
  Fourier transform of the dipole moment elements time series
- [POLARIZABILITY](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.POLARIZABILITY) - Prints
  an element of the Fourier transform of polarizability.
- [RESTART](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.RESTART) - Controls the name of
  the restart file

When [RESTART](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.RESTART),
[MOMENTS](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.MOMENTS) and
[FIELD](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.FIELD) are saved into files, one can
continue running the calculation in the same directory for longer time without rerunning the already
calculated time steps. Note that total length of the propagation time controls the energy/frequency
precision, while timestep size controls the energy/frequency range.

### Example Input

A typical input file which runs the RTBSE propagation will have the
[REAL_TIME_PROPAGATION](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION) section similar to this
one

```
&REAL_TIME_PROPAGATION
    &RTBSE ! Start the RTBSE method
    &END RTBSE
    EPS_ITER 1.0E-8 ! Check convergence
    MAT_EXP BCH
    EXP_ACCURACY 1.0E-14 ! Less than EPS_ITER
    INITIAL_WFN RT_RESTART
    APPLY_DELTA_PULSE
    DELTA_PULSE_DIRECTION 1 0 0
    DELTA_PULSE_SCALE 0.0001 ! Small
    &PRINT
        &MOMENTS
            FILENAME MOMENTS
        &END MOMENTS
        &MOMENTS_FT
            FILENAME MOMENTS-FT
            DAMPING 0.1 ! Exponential damping
            START_TIME ! Fourier transform offset
        &END MOMENTS_FT
        &FIELD
            FILENAME FIELD
        &END FIELD
        &POLARIZABILITY
            FILENAME POLARIZABILITY
            ELEMENT 1 1
        &END POLARIZABILITY
    &END PRINT
&END REAL_TIME_PROPAGATION
```

A complete example input file is available in the
[cp2k-examples](https://github.com/cp2k/cp2k-examples) repository.
