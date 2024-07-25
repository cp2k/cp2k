# Real-Time Propagation and Ehrenfest MD

In real time time-dependent DFT, instead of solving the static Schroedinger equation (SE), the aim
is to solve the time dependent SE (TDSE)

$$
i \frac{\partial}{\partial t} \Psi({\bf r},t) = \hat{H}({\bf r},t) \Psi({\bf r},t),
$$

Contrary to the time independent case, there is no variational principle for the total energy, and
the total energy has to be replaced by the quantum mechanical action,

$$
A[\Psi] =  \int^{t_f}_{t_0} dt \langle \Psi(t) | i\frac{\partial}{\partial t} - \hat{H}(t)| \Psi(t)\rangle
$$

for which it is valid that the function $\Phi(t)$ that makes the action stationary will be the
solution. The resulting time-dependent Kohn-Sham (TD-KS) equations read as

$$
i \frac{\partial}{\partial t} \psi_i({\bf r},t) =  \left[  -\frac{\nabla^2}{2}+v_{\text{KS}}({\bf r},t) + \textbf{F}(t) \right] \psi_i({\bf r},t),
$$

with $\phi_i$ the KS orbitals and the time dependent density

$$
\rho({\bf r},t)  = \sum_i f_i |\psi_i({\bf r},t)|^2.
$$

$\textbf{F}(t)$ is a time dependent external field, if it has been specified. By applying an
external electric field, we intend to simulate the interaction between the electron and an external
radiation. The light is then treated classically using, either the length gauge for non periodic
systems, or the velocity gauge for periodic systems.

A rather general derivation for Ehrenfest (non-adiabatic quantum-classical molecular) dynamics can
be obtained starting from the action of a system. For a situation in which the the electrons are
treated quantum mechanically while the nuclei are treated classically the total action can be
written as the sum of the two environments

$$
A = A_c + A_q \qquad A_c = \int_{t_0}^{t_f} \left[ \sum_A \frac{M_A}{2}\dot{\bf R}_A -U({\bf R},t)\right]dt
$$

and the quantum mechanical action as above defined. The equations of motion are then derived by
making the action stationary

$$
\frac{\delta A}{\delta\langle \Psi({\bf r}, t)|} = 0 \qquad \frac{\delta A}{\delta\langle {\bf R}( t)|} = 0
$$

Evaluating these expressions in the framework of TDDFT, the equations of motion become

$$
M\ddot{\bf R}  = -\frac{\partial}{\partial {\bf R}} U({\bf R},t) - \sum_j\left\langle \Psi^j \left| \frac{\partial}{\partial{\bf R}} V_{\text{int}}({\bf r},{\bf R})\right| \Psi^j \right\rangle
$$

for the nuclear motion, wile for the electrons the time dependent SE as given above is used. These
equations are valid for the exact wavefunctions. In the Kohn-Sham approach, the wavefunctions are
replaced by a linear combination of basis functions. For plane waves the equations remain the same,
since they do not depend on the nuclear coordinates and therefore are independent of the change of
the ionic positions during MD. The case of a Gaussian basis set requires a more detailed analysis.
Using Gaussian basis sets for the expansion of the wavefunctions, an implicit dependence of the
wavefunctions on the nuclear positions is introduced. Since in an MD approach the nuclear positions
are function of time, the time derivative in the quantum mechanical action has to be replaced by the
total time derivative

$$
\frac{d}{dt} = \frac{\partial}{\partial t} + \sum_A \frac{\partial {\bf R}_A}{\partial t}\frac{\partial}{\partial {\bf R}_A}.
$$

Due to the introduction of the basis, the independent variables for making the action constant
become the expansion coefficients $C_{j\alpha}(t)$ of the molecular orbital $j$ in the contracted
basis functions $\phi_\alpha({\bf r},t)$

$$
\dot{C}_{j\alpha} = -\sum_{\alpha\beta}S^{-1}_{\beta\gamma}\left(  i H_{\beta\gamma} + B_{\beta\gamma} \right) C_{j\gamma}
$$

where

$$
S_{\alpha\beta} = \langle \phi_\alpha| \phi_\beta \rangle \qquad B_{\alpha\beta} = \langle \phi_\alpha | \frac{d}{dt} \phi_{\beta} \rangle
$$

Hence, in Ehrenfest dynamics an additional contribution due to the derivative of the basis functions
becomes part of the Hamiltonian used in the exponential operator. Instead of being purely imaginary,
the matrix in the exponential of the propagator becomes fully complex.

## Running real time time dependend DFT in CP2K

To run RTP or Ehrenfest MD with CP2K, the corresponding [RUN_TYPE](#CP2K_INPUT.GLOBAL.RUN_TYPE) in
the \[GLOBAL\] (#CP2K_INPUT.GLOBAL) section has to be selected, and the [MD](#CP2K_INPUT.MOTION.MD)
section has to be present to specify the time step length `TIMESTEP` and the desired number of steps
(`STEPS`). It is crucial to set an appropriate time step to propagate the electronic degrees of
freedom, i.e., in the order of atto-seconds. All other input parameters related to the RT-TDDFT run
are specified from the section
[REAL_TIME_PROPAGATION](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION). The simulation needs to
start from the electronic density at $t_0$ (`INITIAL_WFN`). This can be the ground state obtained by
means of an initial SCF optimization (`SCF_WFN`) or providing the restart file of a previously
computed state (`RESTART_WFN`). The propagation can also be restarted setting the initial
wavefunction as `RT_RESTART` and providing the correct restart file, which has a different format
with respect to the SCF restart file. The path to the restart file is in both cases provided by the
usual [WFN_RESTART_FILE_NAME](#CP2K_INPUT.FORCE_EVAL.DFT.WFN_RESTART_FILE_NAME) key.

Three different propagators are available in CP2K, the enforced time-reversible symmetry propagator
(ETRS), the exponential midpoint (EM) propagator, and the Crank-Nicholson propagator which can be
seen as a first order Pad'e approximation of the EM propagator. The ETRS approach starts with an
exponential approximation to the evolution operator $\hat{U}(t,0)=\exp{-it\hat{H}(0)}$, to then
compute the final time-reversible and unitary propagator self-consistently. In the real time
propagation scheme (fixed ionic positions) the self consistent solution only involves the
calculation of the new Kohn-Sham matrix for the propagated coefficients. For Ehrenfest MD, the
iterative procedure is embedded into the integrator of the nuclear EOM. Hence, each iteration step
for Ehrenfest dynamics involves a real time propagation step and a complete evaluation of the
forces. Due to this, more iterations are needed to reach self consistency for the propagator. The
convergence criterion implemented in CP2K is defined as

$$
||  \Delta C^T S \Delta C||_{\text{}max} < \epsilon
$$

with $\Delta C$ as the difference of coefficient matrices in two successive steps, and $\epsilon$
given by `EPS_ITER`.

In terms of computational cost, the most expensive part is the evaluation of the matrix exponential
in the propagator. Four different methods among those listed in have been implemented in CP2K, i.e.,
the Taylor expansion, he diagonalization, the Pad'e approximation, and the Arnoldi subspace
iteration, to be selected via the `MAT_ESP` input key. The Arnoldi method often provides a superior
performance. Comparing the theoretical scaling, the Arnoldi method is expected to be about 5 times
as fast as Pad'e or Taylor. However, the Pade approximation can be sometime the faster and more
stable choice than the Arnoldi method (e.g., large time step). Since propagation schemes require an
iterative procedure, it is convenient to apply an extrapolation scheme in order to speed up
convergence. For Ehrenfest dynamics, an extrapolation based on the product of the density and the
overlap matrix turns out to be the a good choice. Nevertheless, the quality of the different
extrapolation strongly depends on the system applied to and the settings of the simulation.
Therefore, which method to use and the extrapolation order needs to be tested on the system of
interest. For very large systems, the density-based method can be used to achieve linear scaling
performance (see the
[DENSITY_PROPAGATION](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.DENSITY_PROPAGATION)
keyword).

## The external electric field

One of the most relevant application domains for RT-TDDFT is the study of light?matter interactions,
e.g., in the field of spectroscopy, excited state dynamics and radiation damage. To mimic these
phenomena, at any time during the propagation, it is possible to apply a time dependent electric
field ${\bf E}(t)$. The applied field is in general modulated by an envelope function,
$E_{\text{env}}(t)$ and is defined as:

$$
\textbf{E}(t) = \textbf{P} E_{\text{env}}(t) \cos \left( \omega_{0} t + \phi \right).
$$

Where $\textbf{P}$ is the field polarization, $\omega_{0}$ is the carrying frequency at which the
field oscillates, and $\phi$ its initial phase. The characteristic of the applied field as well as
its time extension are provided through the section [EFIELD](#CP2K_INPUT.FORCE_EVAL.DFT.EFIELD). The
time dependent electric field defined within this section can only be used in combination with
RT-TDDFT. By default the coupling between the electric field and the electronic degrees of freedom
is described within the length gauge, by adding to the Hamiltonian the dipole coupling term
$ e \textbf{E}(t) \cdot {\bf r} $. This approach is only valid for isolated molecular systems, and
not when periodic boundary conditions are applied. The velocity-gauge form of the equations suitable
for periodic systems is obtained through a gauge transformation involving the vector potential

$$
{\bf A}(t) = c \int^t {\bf E}(t') dt' .
$$

In the time dependent KS equations, the vector potential ${\bf A}(t)$ appears in the kinetic energy
term and, in the case where non-local pseudopotentials are used, the gauge field also transforms the
electron?ion interaction. To use this representation the `VELOCITY_GAUGE` input keyword has to be
activated. The total energy varies with time as the applied external field interacts with the
system. To monitor the time evolution of the field as well as of the various terms contributing to
the total electronic energy, the corresponding print key can be activated from the
[REAL_TIME_PROPAGATION](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION) section.

with $\Delta C$ as the difference of coefficient matrices in two successive steps, and $\epsilon$
given by `EPS_ITER`.

## Resonant Excitation

The input file provided below starts a RT--TDDFT simulation for an isolated carbon monoxide molecule
perturbed by a resonant X-ray field. Such real-time calculation aims to promote electrons from some
core occupied state to empty. What states are going to be addressed depends on the frequency of the
field, which has to be resonant to some core state transition. The electronic dynamics are
propagated as the field is applied, leading to an electronic excitation increasing with time, if the
field frequency matches the aimed electronic transition energy.

This equation is propagated for each $i$ electron with a time step $\delta t$: the MOs evolved
collectively over time in a **continuous manner**. Especially, the energy evolution is continuous:
the electrons cannot absorb exactly one photon at a given frequency to instantaneously go from one
electronic state to another.\
If one applies a field resonant with a given electronic transition,
the MOs involved in this transition will gradually be transformed from their initial state to the
corresponding excited state. This will happen in real-time without defining any specific electronic
transition before starting the simulation. Some preliminary calculations to determine the spectral
features of the system under study are in general useful to better define the desired properties of
the perturbing field. For core state excitations these information can be obtained by XAS
simulations.

Once the RTP has started, the evolution of the electronic structure can be monitored by means of
several descriptors, like the time dependent dipole moment, current density, total electronic
density, as well as the projection of the time dependent orbitals onto some preselected reference
states. This latter is particularly useful to verify the variation in population of specific states,
while monitoring density differences is usually helpful to detect charge transfer processes.

By projecting the time-dependent molecular orbitals onto some reference states (e.g., the initial
MOs) means computing the overlap between the propagated orbital
$\psi_i({\bf r}, t) = \sum_\alpha C_{i\alpha}(t)\phi({\bf r})$ and any reference orbital
$\psi^{\text{ref}}_m$. For instance, considering as reference orbitals the static unoccupied ones,
the quantity

$$
N_{\text{exc}}(t)= \sum_m^{\text{unocc}}\sum_i^{\text{occ}} \left\|  \langle  \psi^{\text{ref}}_m |\psi_i(t)\rangle \right\|^2
$$

is an estimate of the number of excited electrons.

### CP2K input

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

In the following the chosen parameters for real-time propagation, electric field, and the
time-dependent projection are discussed.

### Real-Time Propagation parameters

In this simulation, we start from an SCF-optimized state which corresponds to the ground state by
setting [INITIAL_WFN](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.INITIAL_WFN) to `SCF_WFN`.
Before the Real-Time Propagation starts, a ground state calculation will thus be performed. We use a
Molecular Orbital-based description of the wave function for the propagation along with the Arnoldi
approach to compute the exponential by setting
[MAT_EXP](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.MAT_EXP) to `ARNOLDI`. For each time
step, the ERTS algorithm is used with a convergence threshold defined by
[EPS_ITER](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.EPS_ITER), which we set to `1.0E-11`.
Note that the smaller this threshold, the more iterations per time step will be needed to converge.
Hence, we have set the maximal iteration number
[MAX_ITER](#CP2K_INPUT.FORCE_EVAL.DFT.REAL_TIME_PROPAGATION.MAX_ITER) at quite high value of 100.
The time step used is rather small since we have to describe a core-hole excitation process that
takes place with a typical frequency of 529 eV. Following the rule of thumb to set the time step to
about 10 times the field frequency, we set [TIMESTEP](#CP2K_INPUT.MOTION.MD.TIMESTEP) to
`[fs] 0.00078`

### Field parameters

The field is defined by its envelope, its intensity, its polarization along the laboratory x, y, and
z-axis, its wavelength, and the original phase. Several types of field envelopes can be used. Here
we use a Gaussian one with a width of $\sigma=0.3073$ fs and centered at $T0=1.3190$ fs. The
intensity used is 4.08E+13 W.cm$^{-2}$. Along with a carrying frequency of 529 eV (approximately
2.34374655955 nm), it should promote about $10^{-3}$ from the Oxygen 1s to the first available
excited state.

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
(the $i$s) and 7 reference MO (the $m$s), so that there will be $7x7=49$ projection per time step:

$$
n_{mi} (t) = \left|\langle \psi_m(t=0) | \psi_i(t) \rangle\right|^2  = |\sum_{\alpha\beta } C_{m\alpha}(t=0)^* C_{i\beta}(t) S_{\alpha\beta} |^2
$$

where $S_{\alpha \beta}$ is the overlap matrix between basis set orbitals.

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
time, the reference wave function is supposed to be from an XAS_TDP calculation, see Linear Response
input file in the
[tutorial archive](https://github.com/cp2k/cp2k-examples/tree/master/rtp_field_xas). Running with
the proper parameters, this XAS_TDP run saves one .wfn file per excited state. Using
`REFERENCE_TYPE XAS_TDP`, the projection uses automatically the state saved in the produced wave
function file as the reference:

$$
n_{\omega}^i(t) = |<\omega | \psi_i(t)>|^2 = |\sum_{ab} \left( C_\omega^a \right)^* c_i^b(t) S_{ab} |^2
$$

where $C_{\omega\alpha}$ is the $a^{\text{th}}$ atomic coefficient of the excited state found in the
XAS_TDP module. There will be 7 projections per time step in this case: one for each time-dependent
MO. The excited state population associated with this excited state is:

$$
\rho_{\omega}(t) = \sum_i \rho_{\omega}^i(t)
$$

The carbon monoxide molecule has a rotational symmetry along its CO bond. If one notes this axis
$z$, then the $x$ and $y$-axis are equivalent by symmetry. It happens that the first available
excited state for the Oxygen 1s is degenerate: there are two available excited states orthogonal in
the $xy$-plane. That is why a third projection is requested which uses the other excited state
proposed by the XAD_TDP module. Therefore, the excited state should be understood as the sum over
the two equivalent excited states:

$$
n_{Exc}(t) = \rho_{\omega}(t) + \rho_{\omega'}(t)
$$

where $\omega'$ stands for the other equivalent excited state, with $\omega=\omega'=529$ eV.
