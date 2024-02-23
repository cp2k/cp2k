# X-Ray Resonant Excitation

In this tutorial, we will present a simulation of resonant X-Ray excitation of an isolated carbon
monoxide in real-time using a time-dependent field. On this page, you will find an overview of the
method, some equations, and the CP2K input file. A longer version is available in the form of a
jupyter notebook file in [this zip file](https://www.cp2k.org/_media/howto:rtp_field_xas.zip) along
with the simulated data so that you do not need to run this calculation yourself to perform the
analysis. This kind of calculation is not easy to grasp: do not hesitate to have a first look before
diving into the equations and details! This tutorial is connected to this article REF where you can
find complementary information.

## Real-Time TDDFT for resonant excitation

Such real-time calculation aims to promote a small number of electrons from one specific state to
another using an electromagnetic field. The electronic dynamics are propagated as the field is
applied, leading to an electronic excitation increasing with time if the field frequency matches the
targetted electronic transition energy.

### Real-Time TDDFT reminders

The light is treated classically using either the length or velocity gauge to interact with the
electrons. The simulation starts usually from the electronic ground state, but one can also
construct another starting state in CP2K. The field is applied as the wave function is propagated in
real-time. This time-dependent field $\textbf{F}(t)$ has typically an envelope, $E_{\text{env}}$,
and is defined as:

$$
\textbf{F}(t) = \textbf{F} E_{\text{env}}(t) \cos \left( \omega_{0} t + \phi \right).
$$

Where $\textbf{F}$ is the field polarization, $\omega_{0}$ is the carrying frequency at which the
field oscillates, and $\phi$ its initial phase.

The DFT equivalent of the Schrodinger equation is solved numerically time step per time step:

$$
i \frac{\text{d}}{\text{d} t} |\psi_i(t)> = \hat{H}^{KS}(t) |\psi_i(t)>.
$$

Where $|\psi_i>$ is the $i^{\text{th}}$ time-dependent Molecular Orbital (MO) describing the
$i^{\text{th}}$ electrons and $\hat{H}^{KS}$ the time-dependent Kohn-Sham Hamiltonian. This
Hamiltonian is time-dependent because of the evolution of the electronic degrees of freedom and the
field evolution. In the case of Ehrenfest dynamics, the nuclei also move in real-time.

This equation is propagated for each $i$ electron with a time step $\Delta t$: the MOs evolved
collectively over time in a **continuous manner**. Especially, the energy evolution is continuous:
the electrons cannot absorb exactly one photon at a given frequency to instantaneously go from one
electronic state to another.\
If one applies a field resonant with a given electronic transition,
the MOs involved in this transition will gradually be transformed from their initial state to the
corresponding excited state. This will happen in real-time without defining any specific electronic
transition before starting the simulation but only the field parameters. Yet, it is better to
characterize first what electronic transition we want to address before running such Real-Time
Propagation (RTP) calculation. To do so, we will use the Linear Response TDDFT scheme.

### Connection between Real-Time and Linear Response

For X-Rays, we use the XAS_TDP scheme implemented in CP2K, see for instance [this tutorial](./tddft)
for more information. For isolated carbon monoxide and with our calculation parameters (PBEh
functional and PCSEG2 basis set), the first available excited state for the Oxygen 1s has a resonant
frequency of $\omega_{res} = \omega = 529$ eV with an associated oscillator strength of 0.044 a.u.
Note that this result is also confirmed using a real-time approach,
[the $\delta$-kick technic](./delta-kick).

Hence, applying an electric field with a carrying frequency of 529 eV should promote electrons from
the oxygen 1s toward this excited state, noted $|\omega>$. Within the RT-TDDFT approach, it means
that the Molecular Orbital corresponding to the Oxygen 1s, $|1s(t)>$, will be transformed with time
as:

$$
|1s(t)>  = \sqrt{\rho_{GS}(t)} |1s> + \sqrt{\rho_{Exc}(t)} |\omega>.
$$

$\rho_{GS}$ is the time-dependent ground state population and $\rho_{Exc}$ the excited state one.
Note that this equation leaves aside the time-dependent phase factor between these two states.

If the field is resonant with only one electronic transition and if we assume that the rest of the
electrons do not evolve we can use the generalized version of the Rabi oscillation formula to get a
prediction of the excited state population with time:

$$
\rho_{\text{Exc}}(t) = \sin  \left( \frac{ | \mathbf{\mu}^{res} \cdot \mathbf{F}| \times A(t)}{2} \right)^2; \ \ \ 
A(t) = \int_{- \infty}^t E_{env}(t') dt'
$$

Where $\mathbf{\mu}^{res}$ is the transition dipole moment between the ground and excited state,
$\mathbf{F}$ is the polarization of the field, and $A(t)$ is the integral of the field envelope. In
today's case where the amount of electron promoted is small, one can use this approximate formula:

$$
\rho_{\text{Exc}}(t) \approx    \tau \left( \frac{A(t)}{2} \right)^2.
$$

Where $\tau  = |\mathbf{\mu}^{res} \cdot \mathbf{F}|^2 $ can be seen as an instantaneous promoting
rate for the core electron upon the resonant field.

### Real-Time evolution expectations

In the perturbative regime, one can expect the above-mentioned excited population prediction to be
correct as the oscillator strength is provided by Linear Response. Thus, one can track in real-time
these ground and excited population using projection of the time-dependent molecular orbital:

$$
\rho_{GS}(t) = |<1s | 1s(t)>|^2 ; \ \ \ \rho_{Exc}(t) = |<\omega | 1s(t)>|^2
$$

Moreover, we should also have: $\rho_{Exc}(t) = 1 -  \rho_{GS}(t)$. These projections are written
for the Molecular Orbital originally corresponding to the Oxygen 1s, but we can also write them for
the other electrons. Yet, in the case of core X-Ray excitation, we do not expect the other electrons
to be affected within the perturbative regime: for any Molecular Orbital $i$ which is not the one
corresponding to the Oxygen 1s we should have:

$$
|\psi_i(t)> \approx e^{ -i E_i t} |\psi_i(t=0)>,
$$

with $E_i$ the energy associated to the \$i^\[\\text\{th}} molecular orbital. Let us define the
projection for all the molecular orbital either toward their ground state or the targetted excited
state:

$$
\rho^i_{GS}(t) = |<\psi_i(t=0) | \psi_i(t)>|^2 ; \ \ \ \rho^i_{Exc}(t) = |<\omega | \psi_i(t)>|^2
$$

For all MOs except the ones corresponding to the Oxygen 1s we should have $\rho^i_{GS}(t) = 1$ and
$\rho^i_{Exc}(t) \approx 0$. Therefore, here is the generalized formula for the ground state
population and the excited state one for the whole wave-function:

$$
\rho_{GS}(t) = \sum_{i=1}^{N_e} |<\psi_i(t=0) | \psi_i(t)>|^2 ; \ \ \  \rho_{Exc}(t) = \sum_{i=1}^{N_e} |<\omega| \psi_i(t)>|^2.
$$

We should have in the perturbative regime $\rho_{GS}(t)$ close to the number of electron $N_e$:
$\rho_{GS}(t) = N_e - \rho_{Exc}(t)$.

Finally, if we neglect the off-resonance interaction with the field, the evolution of the energy is:

$$
\Delta E(t) = E(t) - E_{GS} = \hbar \omega \rho_{Exc}(t)
$$

## CP2K input

RTP.inp is the input file used for such simulation:

```none
@set EXC_STATE_1     LR-xasat2_1s_singlet_idx1-1.wfn
@set EXC_STATE_2     LR-xasat2_1s_singlet_idx2-1.wfn

&GLOBAL
  PROJECT RTP
  RUN_TYPE RT_PROPAGATION
  PRINT_LEVEL MEDIUM
&END GLOBAL

&MOTION
  &MD
    ENSEMBLE NVE
    STEPS 5000
    TIMESTEP [fs] 0.00078
    TEMPERATURE [K] 0.0
  &END MD
&END MOTION

&FORCE_EVAL
  METHOD QS
  &DFT
    &REAL_TIME_PROPAGATION
      MAX_ITER 100
      MAT_EXP ARNOLDI
      EPS_ITER 1.0E-11
      INITIAL_WFN SCF_WFN
      &PRINT
        &FIELD
          FILENAME =applied_field
        &END FIELD
        &PROJECTION_MO
          REFERENCE_TYPE SCF
          REF_MO_FILE_NAME RTP-RESTART.wfn
          REF_MO_INDEX -1
          SUM_ON_ALL_REF .FALSE.
          TD_MO_INDEX -1
          SUM_ON_ALL_TD .FALSE.
          &PRINT
            &EACH
              MD 1
            &END EACH
          &END PRINT
        &END PROJECTION_MO
        &PROJECTION_MO
          REFERENCE_TYPE XAS_TDP
          REF_MO_FILE_NAME ${EXC_STATE_1}
          TD_MO_INDEX -1
          SUM_ON_ALL_TD .FALSE.
          &PRINT
            &EACH
              MD 1
            &END EACH
          &END PRINT
        &END PROJECTION_MO
        &PROJECTION_MO
          REFERENCE_TYPE XAS_TDP
          REF_MO_FILE_NAME ${EXC_STATE_2}
          TD_MO_INDEX -1
          SUM_ON_ALL_TD .FALSE.
          &PRINT
            &EACH
              MD 1
            &END EACH
          &END PRINT
        &END PROJECTION_MO
      &END PRINT
    &END REAL_TIME_PROPAGATION
    &EFIELD
      ENVELOP GAUSSIAN
      ! gaussian env in fs units
      &GAUSSIAN_ENV
        SIGMA [fs] 0.3073
        T0 [fs] 1.3190
      &END GAUSSIAN_ENV
      ! in W.cm-2
      INTENSITY 4.08E+13
      PHASE 0.0
      POLARISATION 1 0 0
      ! this is 529 eV:
      WAVELENGTH 2.34374655955
    &END EFIELD
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

We will not detail all the parameters and instead focus on the one related to the real-time
propagation, the field, and the time-dependent projection part.

### Real-Time Propagation parameters

In this simulation, we start from an SCF-optimized state which corresponds to the ground state by
setting [INITIAL_WFN](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.INITIAL_WFN) to `SCF_WFN`.

Before the Real-Time Propagation starts, a ground state calculation will thus be performed. We use a
Molecular Orbital-based description of the wave function for the propagation along with the Arnoldi
approach to compute the exponential by setting
[MAT_EXP](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.MAT_EXP) to `ARNOLDI`.

Note that for extensive systems, the density-based method can be used for linear scaling (see the
[DENSITY_PROPAGATION](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.DENSITY_PROPAGATION)
keyword). For each time step, the AERTS algorithm is used with a convergence threshold defined by
[EPS_ITER](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.EPS_ITER), which we set to `1.0E-11`.

Note that the smaller this threshold, the more iterations per time step will be needed to converge.
Hence, we have set the maximal iteration number
[MAX_ITER](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.MAX_ITER) at quite high value of 100.

The time step used is rather small since we have to describe a core-hole excitation process that
takes place with a typical frequency of 529 eV. Following the rule of thumb to set the time step to
about 10 times the field frequency, we set [TIMESTEP](#CP2K_INPUT.MOTION.MD.TIMESTEP) to
`[fs] 0.00078`

### Field parameters

The interaction between the field and the electron is done through the length gauge by default. For
periodic systems (or isolated ones), use the velocity gauge by setting the
[VELOCITY_GAUGE](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.VELOCITY_GAUGE) keyword to True.
The time-dependent electric field is defined for both gauges in the
[EFIELD](#CP2K_INPUT.FORCE_EVAL.DFT.EFIELD) section. It should be noted that one can define several
[EFIELD](#CP2K_INPUT.FORCE_EVAL.DFT.EFIELD) sections to apply several fields within the same
simulation.

The field is defined by its envelope, its intensity, its polarization along the laboratory x, y, and
z-axis, its wavelength, and the original phase. Several types of field envelopes can be used. Here
we use a Gaussian one with a width of $\sigma=0.3073$ fs and centered at $T0=1.3190$ fs. The
intensity used is 4.08E+13 W.cm$^{-2}$. Along with a carrying frequency of 529 eV (approximately
2.34374655955 nm), it should promote about $10^{-3}$ from the Oxygen 1s to the first available
excited state. See
[the jupyter notebook tutorial](https://www.cp2k.org/_media/howto:rtp_field_xas.zip) for more
information about how to choose these parameters depending on your system.

```none
&EFIELD
  ENVELOP GAUSSIAN
  &GAUSSIAN_ENV
    SIGMA [fs] 0.3073
    T0 [fs] 1.3190
  &END GAUSSIAN_ENV
  INTENSITY 4.08E+13
  PHASE 0.0
  POLARISATION 1 0 0
  WAVELENGTH 2.34374655955
&END EFIELD
```

The total number of [STEPS](#CP2K_INPUT.MOTION.MD.STEPS) is set to `5000` so that the field envelope
fits in the period of time simulated.

### Time-Dependent Projection parameters

Time-dependent projections can be defined in the
[REAL_TIME_PROPAGATION/PRINT](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT) section. One
can define several
[PROJECTION_MO](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.PROJECTION_MO) sections which
will be done consequitively and the results stored in different files.

First, let us have a look at the projection toward the ground state:

```none
&PROJECTION_MO
   REFERENCE_TYPE SCF
   REF_MO_FILE_NAME RTP-RESTART.wfn
   REF_MO_INDEX -1
   SUM_ON_ALL_REF .FALSE.
   TD_MO_INDEX -1
   SUM_ON_ALL_TD .FALSE.
   &PRINT
     &EACH
       MD 1
     &END EACH
   &END PRINT
&END PROJECTION_MO
```

All the time-dependent Molecular Orbitals are projected by setting
[TD_MO_INDEX](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.PROJECTION_MO.TD_MO_INDEX) to
`-1` and these projections are stored separately by setting
[SUM_ON_ALL_TD](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.PROJECTION_MO.SUM_ON_ALL_TD)
to `.FALSE.`.

This calculation is spin-independent so one does not have to define the spin of the MO to project.
The reference to projected to is loaded from the file `RTP-RESTART.wfn`, which is the ground state
obtained after the SCF cycles. Note that you can define a wave-function that is not the ground state
as long as the wave-function description (basis set, number of electrons...) is the same as the one
used for the real-time propagation. The
[REFERENCE_TYPE](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.PROJECTION_MO.REFERENCE_TYPE)
defaults to `SCF`.

All the molecular orbitals available in this reference wave-function will be used as a reference for
the projection by setting
[REF_MO_INDEX](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.PROJECTION_MO.REF_MO_INDEX) to
`-1`, and stored in separate files by setting
[SUM_ON_ALL_REF](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.PROJECTION_MO.SUM_ON_ALL_REF)
to `.FALSE.`. The projection will be computed at each time step.

There are $N_e/2 = 7$ molecular orbitals for carbon monoxide. There are thus 7 time-dependent MOs
(the $i$s) and 7 reference MO (the $j$s), so that there will be $7x7=49$ projection per time step:

$$
\rho_j^i(t) = |<\psi_j(t=0) | \psi_i(t)>|^2> = |\sum_{ab} c_j^a(t=0)^* c_i^b(t) S_{ab} |^2
$$

Where $c_j^a(t=0)$ is the $a^{\text{th}}$ atomic coefficient of the $j^{\text{th}}$ reference MO,
$c_i^b(t)$ the $b^{\text{th}}$ atomic coefficient of the $i^{\text{th}}$ time-dependent MO and
$S_{ab}$ the overlap matrix between atomic orbital $a$ and $b$. Using this definition, the
time-dependent ground state population is:

$$
\rho_{GS} = \sum_{i,j} \rho_j^i(t)
$$

Now, let us focus on the second and third projections which can be viewed as projections toward
excited-states:

```none
&PROJECTION_MO
   REFERENCE_TYPE XAS_TDP
   REF_MO_FILE_NAME ${EXC_STATE_1}
   TD_MO_INDEX -1
   SUM_ON_ALL_TD .FALSE.
   &PRINT
      &EACH
         MD 1
      &END EACH
   &END PRINT
&END PROJECTION_MO
```

In this case, all the time-dependent MO are involved in the projection and stored separately. This
time, the reference wave function is supposed to be from an XAS_TDP calculation, see
XAS_TDP%PRINT%RESTART_WFN section or the Linear Response input file in the
[tutorial archive](https://www.cp2k.org/_media/howto:rtp_field_xas.zip). Running with the proper
parameters, this XAS_TDP run saves one .wfn file per excited state: it is the $|\omega>$ presented
before. Using REFERENCE_TYPE XAS_TDP, the projection uses automatically the $|\omega>$ saved in the
.wfn file as the reference file:

$$
\rho_{\omega}^i(t) = |<\omega | \psi_i(t)>|^2 = |\sum_{ab} \left( c_\omega^a \right)^* c_i^b(t) S_{ab} |^2
$$

Where $c_\omega^a$ is the $a^{\text{th}}$ atomic coefficient of the excited state found in the
XAS_TDP module. There will be 7 projections per time step in this case: one for each time-dependent
MO. The excited state population associated with this excited state is:

$$
\rho_{\omega}(t) = \sum_i \rho_{\omega}^i(t)
$$

The carbon monoxide molecule has a rotational symmetry along its CO bond. If one notes this axis
$z$, then the $x$ and $y$-axis are equivalent by symmetry. It happens that the first available
excited state for the Oxygen 1s is degenerate: there are two available excited states orthogonal in
the $xy$-plane. That is why a third projection is requested which uses the other excited state
proposed by the XAD_TDP module.

Therefore, the excited state should be understood as the sum over the two equivalent excited states:

$$
\rho_{Exc}(t) = \rho_{\omega}(t) + \rho_{\omega'}(t)
$$

Where $\omega'$ stands for the other equivalent excited state, with $\omega=\omega'=529$ eV.

## Conclusion

In this tutorial, we have presented the ingredients needed to perform a Real-Time simulation with an
explicit field in order to trigger a specific electronic transition. Such approaches can be used for
other frequencies without much changes in the parameters and also extended to Ehrenfest dynamics.
Note that in this case, the projection toward the reference state is less trivial to define, see the
[PROPAGATE_REF](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.PRINT.PROJECTION_MO.PROPAGATE_REF)
keyword.

This has been done within the perturbative regime (the amount of excited electron is small compared
to 1) and [the jupyter notebook](https://www.cp2k.org/_media/howto:rtp_field_xas.zip) presents an
extensive connection to the Linear Response prediction and also discusses the obtained results.
Please note that applying a field with the intention to promote about one electron can bring more
challenges as well as applying a field with a large amplitude.
