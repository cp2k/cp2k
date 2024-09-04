# X-Ray Absorption from RTP and $\delta$-Kick perturbation

This tutorial shows how to run a Real-Time Time-Dependent DFT calculation using the so-called
$\delta$-kick approach to compute absorption electronic spectra.

For the corresponding Linear-Response approach, read the tutorial [the X-Ray frequencies](./tddft).

This tutorial is structured as follows:

- Brief theory overview
- CP2K input file overview

Along the RTP run, the time dependent dipole moment is sampled and it is then Fourier transformed to
produce the absorption spectrum in the frequency domain. An example script is provided
[here](https://github.com/cp2k/cp2k-examples/tree/master/x-ray/delta_kick_analyse).

## Theory Overview

### Response theory

In the perturbative regime and at the linear order, the time-dependent response of an electronic
cloud can be decomposed in the Fourier space: the response at a given frequency is proportional to
the perturbation acting on the system at the same frequency. In particular, the induced dipole
moment (or equivalently the induced current) at a given frequency, $\mu^\omega$, is given by the
product of the polarizability matrix, $\alpha(\omega, \omega)$, with the perturbative electric field
at the same frequency, $F^\omega$, within the dipolar approximation.

$$
\mu^\omega = \alpha(\omega, \omega) \cdot F^\omega
$$

This quantity involves a dot product between the field vector and the polarizability matrix. For the
next paragraphs, we will drop the vectorial behavior and discuss it later.

In Linear-Response-TDDFT approaches, one computes the excited state promoted by an electric field
oscillating at a specific frequency $\omega$ and then the transition dipole moment associated with
such electronic transition. This transition dipole moment is then used to derive the polarizability
for resonant or non-resonant frequency.

In the Real-Time approach presented here, we will excite **all** the possible electronic transitions
at the beginning of the simulation by providing an istantaneous pulse containing all the frequencies
The perturbed electronic wavefunction is then propagated and the fluctuations of the induced time
dependent dipole moment are recorded . The polarizability tensor at any frequency is finally derived
from the induced dipole moment.

### The $delta$-kick perturbation

To apply the $delta$-kick, i.e., an intense electric field in a very sort time, the electronic
structure is first obtained at the ground state, for instance using DFT, and then perturbed with an
instantaneous electric field:

$$
F(t) = F^0 \delta(t)
$$

The same result is obtained by applying a constant field with a very narrow Gaussian envelope. This
field perturbes the ground state wave-function at $t=0^-$. Then, this excited wave-function is
propagated in real time by numerically integrating the time dependent Schroedinger equation.

Theinstantaneaous field can be written in the frequency domain as

$$
F^\omega = \frac{F^0}{2 \pi} \int_{-\infty}^{+ \infty} \delta(t) e^{i \omega t} dt = \frac{F^0}{2 \pi} e^{i \omega \times 0} = \frac{F^0}{2 \pi}
$$

The field amplitude in the Fourier space is $F^0 / 2 \pi$ for all frequencies: this instantaneous
perturbation does indeed contain all the frequencies. The time-dependent wave-function can be
described to be the ground state one plus all possible excited states. The total dipole moment
results from the superposition of all possible oscillations related to all the excited states. This
complex behavior in the time domain becomes simple in the frequency domain, since we know the
amplitude of the perturbation applied at each frequency:

$$
\mu^\omega = \frac{1}{2 \pi} \alpha(\omega, \omega)  F^0
$$

Therefore, in order to get the polarizability $\alpha(\omega, \omega)$, one prepares the perturbed
state at $t=0^-$ for a given field amplitude of the field, then propagates the wave-function in real
time. The time dependent dipole moment is extracted along the propagation and it is finally Fourier
transformed to calculate the polarizability as

$$
\text{Re} \left[ \alpha(\omega, \omega) \right]  = 2 \pi \frac{\text{Re} \left[ \mu^\omega \right] }{F^0} \\
\text{Im} \left[ \alpha(\omega, \omega) \right] = 2 \pi \frac{\text{Im}  \left[ \mu^\omega \right] }{F^0}
$$

The amplitude of the dipole moment in the Fourier space may be very small, if the frequencies are
far from resonances. But, in principle, one can extract the whole spectrum from one RTP run. As a
matter of fact, the whole spectrum from core to valence excitations is very broad and for numerical
efficiency the propagation parameters are set to focus only on one specific part of the spectrum.

### Absorption spectrum

Under the assumption of linear response to the perturbation,\
the real and imaginary part of the
calculated frequency-dependent polarizability define the nature of the response. If the frequency is
at an electronic resonance, then the polarizability has an imaginary part. The polarizability is
pure real off resonances.

From one Real-Time propagation, one can obtain three components of the polarizability tensor. For
instance, applying a $\delta$-kick along $x$ provides the components $\alpha_{xx}$, $\alpha_{yx}$
and $\alpha_{zx}$, corresponding to the Fourier transform in $x$, $y$ and $z$, respectively.
Therefore, to get the full polarizability tensor, three RTP runs are needed, one per Cartesian axis.

To compare with experiments, one often assumes that the system is averaged over all possible
orientations in space. Then, the absorption spectrum $I(\omega)$ is proportional to:

$$
I(\omega)  \propto \sum_{i=x,y,z} \text{Im} \left[ \alpha_{ii}(\omega, \omega) \right]
$$

## CP2K Input

The following example is to simulate the response of carbon-monoxide in the gas phase when applying
the $\delta$-kick along the perpendicular to the CO bond. The analysis is done within the X-Ray
range, but in principle, any other range of frequency can be sampled using the same approach,
providing that fime step and length of the propagation are opportunely adjusted.

The ifollowing input file `RTP.inp` is for a simulation at the DFT/PBEh level of theory, it uses the
GAPW approah and the density is expandedn in the all-electron PCSEG-2 basis sets.

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

The time step used to propagate the wave-function is adapted to optimally resolve the frequency of
interest. Typically, the time-step choice results from a trade-off between the resolution of the
highest frequency of interest and the computational cost to sufficiently extend the sampling. As a
rule of thumb, the time step should be about one order of magnitude smaller than the period
corresponding to the maximum frequency to be resolved. For core excitation from the Oxygen 1s, the
K-edge is around 530 eV, corresponding to a period of 0.0078 fs.

According to Linear Response TDDFT calculations (using the XAS_TDP module, see
[this tutorial](./tddft)), the first transition occurs at 529 eV, with an oscillator strength of
0.044 a.u. For the $\delta$-kick approach, this means that the induced dipole moment should have an
oscillatory component around the frequency corresponding to 529 eV. Hence, the
[TIMESTEP](#CP2K_INPUT.MOTION.MD.TIMESTEP) is set ten times smaller than the maximum frequency,
i.e., to `[fs] 0.00078`, which allows the sampling of the fasted oscillation by means of at 10
propagation steps.

The total time of the simulation is then determined by the maximum number of time steps, here
500000\. The longer is the sampling, the smoother and sharper is the resulting spectrum. Longer
simulation times give a better resolution of the spectral peaks, i.e., a more accurate determination
of the resonances. There is no always valid rule to decide how long the propagation needs to be
extended, because the fluctuations might be strongly system dependent. It is always a good idea to
check whether the obtained spectrum has converged with respect to further extension of the sampling.

### Field properties

The field perturbation (the $\delta$-kick) is defined by its amplitude and polarization. The
amplitude depends on the system of interest, a typical value is $10^{-3}$. Also in this case a
convergence check by running more than one propagation at different amplitudes is the best way to
asses the parameter. Within the linear regime, the response of the system should double by twice the
field's amplitude. Note that when applying a too-low field, numerical noise might become dominant.
For isolated CO, $10^{-3}$ is a good value.

**Please note that the actual perturbation applied in CP2K is not
[DELTA_PULSE_SCALE](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.DELTA_PULSE_SCALE) and depends
on the cell size. The corresponding amplitude value is written in the output file, as part of the
RTP initialization**

The field's polarization determines what excited states are most probably triggered, according to
the selection rules. For the present example, we know that the electronic transition at 529 eV is
perpendicular to the CO bond. For the geometry we are using, it means that the field polarization
should be along $x$ (or $y$) , i.e.,
[DELTA_PULSE_DIRECTION](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.DELTA_PULSE_DIRECTION) to
`1 0 0` and [DELTA_PULSE_SCALE](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.DELTA_PULSE_SCALE)
to `0.001`.

### Note about the choice of the Gauge

For an isolated system, the length gauge form coupling coupling the position and the electric field
operators can be employed to all perturbative orders by setting
[PERIODIC](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PERIODIC) to `.FALSE.`.

For condensed phase systems, the velocity gauge form, implying a gauge transformation involving the
vector potential is instead required, and in CP2K it is activated by setting
[PERIODIC](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PERIODIC) to `.TRUE.`. In this case, the
perturbation will be applied only within the first order.
