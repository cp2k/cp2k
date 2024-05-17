# X-Ray Absorption from δ-Kick

This tutorial shows how to run a Real-Time Time-Dependent DFT calculation using the so-called δ-kick
approach to compute absorption electronic spectra.

We warmly recommend having a look first at the Linear-Response approach, for instance at
[the X-Ray frequencies](./tddft). See also this article and the supplementary information: LINK

This tutorial is structured as follows:

- Brief theory overview
- CP2K input file overview

Then, the link between the CP2K output (time-dependent dipole moment) and the quantity targetted
(the absorption spectra in the frequency domain) is performed in this jupyter notebook:
[download this file](https://www.cp2k.org/_media/howto:delta_kick_analyse.zip). Of course, one can
use its own script and method. Yet, as these analyses present some subtilities, we propose our which
shows quantitative results with the Linear Response one.

## Theory Overview

### Response theory reminders

In the perturbative regime and at the linear order, the time-dependent response of an electronic
cloud can be decomposed in the Fourier space: the response at a given frequency is proportional to
the perturbation acting on the system at the same frequency. In particular, the induced dipole
moment (or equivalently the induced current) at a given frequency, $\mu^\omega$, is given by the
product of the polarizability matrix, $\alpha(\omega, \omega)$, with the perturbative electric field
at the same frequency, $F^\omega$, within the dipolar approximation.

$$
\mu^\omega = \alpha(\omega, \omega) \cdot F^\omega
$$

This quantity involved a dot product between the field vector and the polarizability matrix. For the
next paragraphs, we will drop the vectorial behavior and discuss it later.

In Linear-Response-TDDFT approaches, one computes the excited state promoted by an electric field
oscillating at a specific frequency $\omega$ and then computes the transition dipole moment
associated with such electronic transition. This transition dipole moment is then used to compute
the polarizability for resonant or non-resonant frequency.

In the Real-Time approach presented here, we will excite **all** the possible electronic transitions
at the beginning of the simulation and then propagate this excited electronic cloud to observe the
induced dipole moment. We will deduce the polarizability tensor for any frequency from this induced
dipole moment.

### The δ-kick technic

In the δ-kick technic, the electronic structure is first obtained at the ground state, for instance
using DFT, and then perturbed with an instantaneous electric field named a δ-kick:

$$
F(t) = F^0 \delta(t)
$$

This δ-kick perturbes the ground state wave-function at $t=0^-$ using either the dipole moment
operator (length gauge) or using the momentum one (velocity gauge). Then, this excited wave-function
is propagated in real time by numerically integrating the DFT equivalent of the Schrodinger
equation. 

The electric δ-kick can be written in the frequency domain:

$$
F^\omega = \frac{F^0}{2 \pi} \int_{-\infty}^{+ \infty} \delta(t) e^{i \omega t} dt = \frac{F^0}{2 \pi} e^{i \omega \times 0} = \frac{F^0}{2 \pi}
$$

The field amplitude in the Fourier space is $F^0 / 2 \pi$ for all frequencies: this instantaneous
perturbation does indeed contain all the frequencies. The time-dependent wave-function can be thus
described to be the one at the ground state plus all the possible excited states. During the time
propagation, the dipole moment of the molecule will be given by the superposition of all possible
oscillations related to all the excited states. This complex behavior in the time domain becomes
simple in the frequency domain, since we know the amplitude of the perturbation applied at each
frequency:

$$
\mu^\omega = \frac{1}{2 \pi} \alpha(\omega, \omega)  F^0
$$

Therefore, in order to get the polarizability $\alpha(\omega, \omega)$, one prepares the perturbed
state at $t=0^-$ for a given field amplitude, then propagates the wave-function in real time. During the propagation the dipole moment is computed at each time step. The dipole moment is finally Fourier transformed and the
polarizability is extracted using:

$$
\text{Re} \left[ \alpha(\omega, \omega) \right]  = 2 \pi \frac{\text{Re} \left[ \mu^\omega \right] }{F^0} \\
\text{Im} \left[ \alpha(\omega, \omega) \right] = 2 \pi \frac{\text{Im}  \left[ \mu^\omega \right] }{F^0}
$$

The amplitude of the dipole moment in the Fourier space may be very small if the frequencies are far
from resonances. But, in principle, one can extract the whole spectrum from one Real-Time
propagation. As a matter of fact, the whole spectrum from core to valence excitations is very broad and for numerical efficiency the propagation parameters are set to focus only on one specific part of the spectrum.

### Absorption spectrum

Using the time-dependent dipole moment, we got the frequency-dependent polarizability $\alpha$
assuming a perturbative and linear response regime. The real and imaginary part of the
polarizability defines the nature of the response. If the frequency is at an electronic resonance,
then the polarizability has an imaginary part. The polarizability is pure real off resonances.

From one Real-Time propagation, one can obtain 3 components of the polarizability tensor. For
instance, applying a δ-kick along the $x$ direction will provides the components $\alpha_{xx}$,
$\alpha_{yx}$ and $\alpha_{zx}$ by Fourier Transform the dipole moment along the $x$, $y$ and $z$
direction respectively. Therefore, to get the full polarizability tensor, one should run 3 Real-Time
propagations using a perturbation along 3 orthogonal directions, for instance, $x$, $y$, and $z$.

To compare with experiments, one often assumes that the system is averaged over all possible
orientations in space. Then, the absorption spectrum $I(\omega)$ is proportional to:

$$
I(\omega)  \propto \sum_{i=x,y,z} \text{Im} \left[ \alpha_{ii}(\omega, \omega) \right]
$$

## CP2K Input

In this tutorial, we will study the response of a carbon-monoxide in the gas phase upon a δ-kick
along its perpendicular direction. We will focus on the X-Ray range, but in principle, other ranges
of frequency can be sampled using this approach.

The input file RTP.inp presents an example to propagate the electronic density of an isolated carbon-monoxide molecule, upon 
the application of a predefined  δ-kick. Ste simulation is performedd at the DFT/PBEh level using a PCSEG-2 basis set. 

```none
&GLOBAL
  PROJECT RTP
  RUN_TYPE RT_PROPAGATION
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE NVE
    STEPS 50000
    TIMESTEP [fs] 0.00078
    TEMPERATURE [K] 0.0
  &END MD
&END MOTION

&FORCE_EVAL
  METHOD QS
  &DFT
    &REAL_TIME_PROPAGATION
      APPLY_DELTA_PULSE .TRUE.
      DELTA_PULSE_DIRECTION 1 0 0
      DELTA_PULSE_SCALE 0.001
      MAX_ITER 100
      MAT_EXP ARNOLDI
      EPS_ITER 1.0E-11
      INITIAL_WFN SCF_WFN
      PERIODIC .FALSE.
    &END REAL_TIME_PROPAGATION
    BASIS_SET_FILE_NAME BASIS_PCSEG2
    POTENTIAL_FILE_NAME POTENTIAL
    &MGRID
      CUTOFF 1000
      NGRIDS 5
      REL_CUTOFF 60
    &END MGRID
    &QS
      METHOD GAPW
      EPS_FIT 1.0E-6
    &END QS
    &SCF
      MAX_SCF 500
      SCF_GUESS RESTART
      EPS_SCF 1.0E-8
    &END SCF
    &POISSON
      POISSON_SOLVER WAVELET
      PERIODIC NONE
    &END POISSON
    &XC
      &XC_FUNCTIONAL PBE
        &PBE
          SCALE_X 0.55
        &END
      &END XC_FUNCTIONAL
      &HF
        FRACTION 0.45
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE TRUNCATED
          CUTOFF_RADIUS 7.0
        &END INTERACTION_POTENTIAL
      &END HF
    &END XC
    &PRINT
      &MULLIKEN OFF
      &END MULLIKEN
      &HIRSHFELD OFF
      &END HIRSHFELD
      &MOMENTS
       PERIODIC .FALSE.
         FILENAME =dipole
         COMMON_ITERATION_LEVELS 100000
         &EACH
            MD 1
         &END EACH
      &END MOMENTS
    &END PRINT
  &END DFT

  &SUBSYS
    &CELL
      ABC 10 10 10
      ALPHA_BETA_GAMMA 90 90 90
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME carbon-monoxide_opt.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND C
      BASIS_SET pcseg-2
      POTENTIAL ALL
    &END KIND
    &KIND O
      BASIS_SET pcseg-2
      POTENTIAL ALL
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
```

### Time step

The time step used to propagate the wave-function is adapted to optimally resolve the frequency of interest.
Typically, the time-step cjoice results from a trade-off between the resolution of the highest frequency of interest  and the computationl cost to sufficiently extend the sampling, such to achieve sufficient accuracy.
As a rule of thumb, the time step should be about one order of magnitude smaller than the period corresponding to the maximum frequency to be resolved. For core excitation from the Oxygen 1s, the K-edge is arounf 530 eV, corresponding to a period of 0.0078 fs.

According to Linear Response TDDFT calculations (using the XAS_TDP module, see
[this tutorial](./tddft)), the first transition occurs at 529 eV with an oscillator strength of
0.044 a.u. For the δ-kick approach, this means that the induced dipole moment should have an
oscillatory component around the frequency corresponding to 529 eV. Hence, the
[TIMESTEP](#CP2K_INPUT.MOTION.MD.TIMESTEP) is set ten times smaller than the maximum frequency, i.e., to `[fs] 0.00078`, which allawos the sampling of  the fasted oscillation by means of at 10 propagation steps.

The total time of the simulation is then given by the number of time steps, here 500000, and it affects the 
smoothness and sharpness of the resulting spectrum. Longer  simulation times give a better resolutin of the spectral peaks, i.e., a more accurate determination of the resonances.
There is no always valid rule to decide how long the propagation needs to be extended, because the fluctuations might be strongly system dependent. It is always a good idea to check whether the obtained spectrum has converged with respect to further extensiont of the sampling. 

### Field properties

The field perturbation (the δ-kick) is defined by its amplitude and polarization.

The amplitude depends on the system of interest, a typical value is $10^{-3}$. Also in this case a convergence check by running more than one propagation at different amplitudes is the best way to asses the parameter. On the other hand, this brute force approach might request too large computational respurces for large systems. 
Within the linear regime, the response of the system should double by dubling the field amplitude. When this is not the
case, the response is not linear and one should reduce the field amplitude. Note that applying
a too-low field may lead to numerical noise. For isolated carbon monoxide, we have checked that
$10^{-3}$ is within the linear regime.

**Please note that the actual perturbation applied in CP2K is not
[DELTA_PULSE_SCALE](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.DELTA_PULSE_SCALE) and depends
on the cell size. Look at the output file just before the Real Time propagation part to know the
exact perturbation amplitude used**

The polarization of the field determines which excited state is triggered: the electric field and
the transition dipole moment should be non-perpendicular. One way to know how to select the
polarization is to run a Linear-Response TDDFT calculation and to look at the transition dipole
moment decomposition along the laboratory axis (or the oscillator strength decomposition). For
instance, have a look at the first part of [this tutorial](./resonant_excitation) on Real Time TDDFT
with a time-dependent field for isolated carbon monoxide. Relying on this calculation, the
electronic transition at 529 eV is perpendicular to the CO bond. For the geometry we are using, it
means that the field polarization should be along x (or y) by setting
[DELTA_PULSE_DIRECTION](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.DELTA_PULSE_DIRECTION) to
`1 0 0` and [DELTA_PULSE_SCALE](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.DELTA_PULSE_SCALE)
to `0.001`.

### Note about the choice of the Gauge

For an isolated system, use the length gauge implementation because the perturbation is applied up to all
perturbative orders by setting [PERIODIC](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PERIODIC)
to `.FALSE.`.

For condensed phase systems, use the velocity gauge instead by setting
[PERIODIC](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PERIODIC) to `.TRUE.`. In this case, the
perturbation will be applied only within the first order. 

## Analyzing the results

To analyze the Real Time propagation simulation, we propose to go through a jupyter notebook:
[download this file](https://www.cp2k.org/_media/howto:delta_kick_analyse.zip). Note that this zip
file also contains the data needed to perform the analyses (the output of the CP2K run).
