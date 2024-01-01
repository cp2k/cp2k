!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE pint_types

   USE cp_log_handling,                 ONLY: cp_logger_type
   USE gle_system_types,                ONLY: gle_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE parallel_rng_types,              ONLY: rng_record_length,&
                                              rng_stream_type
   USE replica_types,                   ONLY: replica_env_type
   USE simpar_types,                    ONLY: simpar_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Energy contributions - symbolic names for indexing energy arrays
   INTEGER, PARAMETER, PUBLIC :: e_conserved_id = 1, &
                                 e_potential_id = 2, &
                                 e_kin_thermo_id = 3, &
                                 e_kin_virial_id = 4

   ! Number of energy contributions for static array allocation
   INTEGER, PARAMETER, PUBLIC :: e_num_ids = 4

   INTEGER, PARAMETER, PUBLIC :: thermostat_none = 0, &
                                 thermostat_nose = 1, &
                                 thermostat_gle = 2, &
                                 thermostat_pile = 3, &
                                 thermostat_piglet = 4, &
                                 thermostat_qtb = 5

   PUBLIC :: pint_env_type
   PUBLIC :: normalmode_env_type
   PUBLIC :: staging_env_type
   PUBLIC :: pile_therm_type
   PUBLIC :: piglet_therm_type
   PUBLIC :: qtb_therm_type

   ! ***************************************************************************
   !> \brief environment for a path integral run
   !> \param p number of replicas/beads
   !> \param nnos nose hoover chain length
   !> \param nrespa number of respa steps
   !> \param nsteps - number of PIMD steps to be performed
   !> \param iter current iteration number
   !> \param ndim number of coordinates per replica/bead
   !> \param transform type of transform (normalmode or staging)
   !> \param t_tol temperature tolerance for rescaling
   !> \param v_tol velocity tolerance for rescaling
   !> \param kT boltzmann factor times temperature (simulation temperature
   !> \param    not necessarily the physical temperature)
   !> \param beta 1/kT (physical temperature)
   !> \param dt time step for dynamic
   !> \param e_pot_h potential energy in harmonic springs
   !> \param e_kin_beads (fictitious) kinetic energy of the beads
   !> \param e_pot_t potential energy of thermostats
   !> \param e_kin_t kinetic energy of thermostats
   !> \param energy - energy contributions updated every step REAL(e_num_ids)
   !> \param    e_kin_virial_id - virial estimator of the (real) kinetic energy
   !> \param t current simulation time
   !> \param replicas replica environment for force calculations
   !> \param input input data structure
   !> \param staging_env description for the staging transformation
   !> \param normalmode_env description for the normal mode transformation
   !> \param randomG random number stream descriptor
   !> \param mass real masses
   !> \param e_pot_bead array with last energies from QS per replica
   !> \param x array with real space coordinates (P, 3*N)
   !> \param v array with real space velocities
   !> \param f array with real space forces
   !> \param mass_beads masses of the beads for harmonic forces (harmonic mass)
   !> \param mass_fict fictitious mass of the beads for dynamics (kinetic mass)
   !> \param ux array with transformed space coordinates (P, 3*N)
   !> \param uv array with transformed velocities
   !> \param uv_t array with temporary transformed velocities
   !> \param uv_new array with new transformed velocities
   !> \param uf array with transformed accelerations (QS part)
   !> \param uf_h array with harmonic part transformed forces
   !> \param tx nose hoover chain positions (pint_env%nnos,pint_env%p,pint_env%ndim)
   !> \param tv nose hoover chain velocities
   !> \param tv_t nose hoover chain velocities (temporary)
   !> \param tv_old nose hoover chain velocities (older)
   !> \param tv_new nose hoover chain velocities (newer)
   !> \param tf nose hoover chain forces (?)
   !> \param Q nose hoover chain masses
   !> \param time_per_step - time per step in seconds (updated every step)
   !> \param pile_therm data used for the pile thermostat
   !> \param wsinex omega*sin(omega*deltat) for exact harminic integrator
   !> \param iwsinex 1/omega*sin(omega*deltat) for exact harminic integrator
   !> \param cosex cos(omega*deltat) for exact harminic integrator
   !> \param propagator contains propagator related constants
   !> \param harm_integrator selects between numeric and exact harmonic integrator scheme
   !> \param first_propagated_mode if 1 - propagate all normal modes,
   !>                              if 2 - keep centoid fixed
   !> \author fawzi
   !> \par History
   !>      Added some comments - hforbert
   !>      Added normal mode transformation - hforbert
   !>      2009-06-15 helium_solvent_type object is no longer a member of
   !>                   pint_env_type [lwalewski]
   !>      2014-10-23 added pile_therm [Felix Uhl]
   !>      2018-02-13 added qtb_therm [Fabien Brieuc]
   ! ***************************************************************************
   TYPE pint_env_type
      INTEGER :: p = 0, nnos = 0, nrespa = 0, iter = 0, ndim = 0, transform = 0
      INTEGER :: first_step = 0, last_step = 0, num_steps = 0, first_propagated_mode = 0
      INTEGER :: pimd_thermostat = 0, harm_integrator = 0, thermostat_rng_seed = 0
      REAL(KIND=dp) :: t_tol = 0.0_dp, v_tol = 0.0_dp, kT = 0.0_dp, beta = 0.0_dp, dt = 0.0_dp, &
                       e_gle = 0.0_dp, e_pile = 0.0_dp, e_piglet = 0.0_dp, e_qtb = 0.0_dp, e_pot_h = 0.0_dp, &
                       e_kin_beads = 0.0_dp, e_pot_t = 0.0_dp, e_kin_t = 0.0_dp, t = 0.0_dp, time_per_step = 0.0_dp
      REAL(KIND=dp) :: link_action = 0.0_dp, pot_action = 0.0_dp
      TYPE(cp_logger_type), POINTER :: logger => NULL()
      TYPE(replica_env_type), POINTER :: replicas => NULL()
      TYPE(section_vals_type), POINTER :: input => NULL()
      TYPE(staging_env_type), POINTER :: staging_env => NULL()
      TYPE(normalmode_env_type), POINTER :: normalmode_env => NULL()
      TYPE(rng_stream_type) :: randomG = rng_stream_type()
      TYPE(gle_type), POINTER        :: gle => NULL()
      REAL(KIND=dp), DIMENSION(e_num_ids) :: energy = 0.0_dp
      REAL(KIND=dp), DIMENSION(:), POINTER :: mass => NULL(), e_pot_bead => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: x => NULL(), v => NULL(), f => NULL(), mass_beads => NULL(), &
                                                 mass_fict => NULL(), ux => NULL(), ux_t => NULL(), uv => NULL(), uv_t => NULL(), &
                                                 uv_new => NULL(), uf => NULL(), uf_h => NULL(), external_f => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: centroid => NULL()
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: tx => NULL(), tv => NULL(), tv_t => NULL(), &
                                                    tv_old => NULL(), tv_new => NULL(), tf => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: Q => NULL() ! dim p, make it (p,ndim)?
      REAL(KIND=dp), DIMENSION(:), POINTER :: rtmp_ndim => NULL(), rtmp_natom => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: iwsinex => NULL(), wsinex => NULL(), cosex => NULL()
      TYPE(pile_therm_type), POINTER       :: pile_therm => NULL()
      TYPE(piglet_therm_type), POINTER     :: piglet_therm => NULL()
      TYPE(qtb_therm_type), POINTER        :: qtb_therm => NULL()
      TYPE(pint_propagator_type), POINTER  :: propagator => NULL()
      TYPE(simpar_type), POINTER           :: simpar => NULL()
      INTEGER                              :: n_atoms_constraints = 0
      INTEGER, DIMENSION(:), POINTER       :: atoms_constraints => NULL()
      LOGICAL                              :: beadwise_constraints = .FALSE.
      REAL(KIND=dp)                        :: kTcorr = 0.0_dp

   END TYPE pint_env_type

   ! ***************************************************************************
   !> \brief data to perform the normalmode transformation
   !> \note
   !>    p             - number of beads
   !>    Q_bead        - thermostat mass for a non-centroid bead
   !>    Q_centroid    - thermostat mass for a centroid degree of freedom
   !>    modefactor    - mass scale factor for non-centroid degrees of freedom
   !>    harm          - factor for harmonic potential ( w_p^2/modefactor )
   !>    x2u           - transformation matrix real coord to normal mode space
   !>    u2x           - transformation matrix normal mode coord to real space
   !>    lambda        - propagator frequencies of the ring polymer
   !>
   !>    This could be done via FFT calls as well, but for now...
   !> \author hforbert
   ! ***************************************************************************
   TYPE normalmode_env_type
      INTEGER :: p = 0
      REAL(KIND=dp) :: Q_bead = 0.0_dp, Q_centroid = 0.0_dp, modefactor = 0.0_dp, harm = 0.0_dp
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: x2u => NULL(), u2x => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: lambda => NULL()
   END TYPE normalmode_env_type

   ! ***************************************************************************
   !> \brief data to perform the staging transformation
   !> \note
   !>    nseg
   !>    j
   !>    p
   !>    w_p
   !>    w_j
   !>    Q_stage
   !>    Q_end
   !> \author fawzi
   ! ***************************************************************************
   TYPE staging_env_type
      INTEGER :: nseg = 0, j = 0, p = 0
      REAL(KIND=dp) :: w_p = 0.0_dp, w_j = 0.0_dp, Q_stage = 0.0_dp, Q_end = 0.0_dp
   END TYPE staging_env_type

   ! ***************************************************************************
   !> \brief data to use the pile thermostat
   !> \note
   !>    lamb          - coupling constant of pile to the normal modes
   !>    tau           - time constant for centroid mode
   !>    thermostat_energy        - energy difference for conxerved quantity
   !>    c1            - scaling of the old momenta
   !>    c2            - scaling of the friction term
   !>    g_fric        - mode specific friction
   !>    massfact      - Mass prefactor to get units right
   !>    gaussian_rng_stream - random number generator
   !> \author Felix Uhl
   ! ***************************************************************************
   TYPE pile_therm_type
      REAL(KIND=dp)                              :: lamb = 0.0_dp, tau = 0.0_dp, thermostat_energy = 0.0_dp
      REAL(KIND=dp), DIMENSION(:), POINTER       :: c1 => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER       :: c2 => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER       :: g_fric => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: massfact => NULL()
      TYPE(rng_stream_type)                      :: gaussian_rng_stream = rng_stream_type()
   END TYPE pile_therm_type

   ! ***************************************************************************
   !> \brief data to use the piglet thermostat
   !> \note
   !>    ndim          - number of degrees of freedom
   !>    p             - trotter number
   !>    nsp1          - number of additional degrees of freedom for Markovian
   !dynamics + 1
   !>    thermostat_energy        - energy difference for conxerved quantity
   !>    a_mat         - A matrices (9,9,P)
   !>    c_mat         - C matrices (9,9,P)
   !>    gle_t         - Deterministic part of propagator
   !>    gle_s         - Stochastic part of propagator
   !>    smalls        - Keeps a copy of momenta and additional degrees of
   !freedom
   !>                    to ensure Markovian dynamics
   !>    temp1         - Big storage array that is needed on the way
   !>    temp2         - vector to store the random numbers
   !>    sqrtmass      - contains the squareroot of the dynamical masses
   !>    gaussian_rng_stream - random number generator
   !> \author Felix Uhl
   ! ***************************************************************************
   TYPE piglet_therm_type
      INTEGER                              :: ndim = 0, p = 0, nsp1 = 0
      REAL(KIND=dp)                        :: thermostat_energy = 0.0_dp
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: a_mat => NULL(), c_mat => NULL()
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: gle_s => NULL(), gle_t => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: smalls => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: temp1 => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: temp2 => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: sqrtmass => NULL()
      TYPE(rng_stream_type)                :: gaussian_rng_stream = rng_stream_type()
   END TYPE piglet_therm_type

   ! ***************************************************************************
   !> \brief data to use the qtb thermostat
   !> \note
   !>    tau           - time constant (1/friction) for centroid mode
   !>    lamb          - scaling of time constants to the ring polymer NM freq.
   !>    taucut        - inverse of frequency cutoff for QTB forces
   !>    lambcut       - scaling of the cutoff angular freq. to the ring polymer
   !>    c1            - scaling of the old momenta
   !>    c2            - scaling of the friction term
   !>    g_fric        - mode specific friction
   !>    massfact      - Mass prefactor to get units right
   !>    rf            -  stores the QTB forces
   !>    h             - filter for computation of QTB forces
   !>    r             - store random numbers for computation of QTB forces
   !>                  - NM freq.
   !>    step          - update QTB forces every qtb_step
   !>    cpt           - to know when to draw new random forces (every qtb_step)
   !>    fp            - defines if we use f_P^(0) or f_P^(1)
   !>    nf            - nb of points used for the convolution product (memory)
   !>    gaussian_rng_stream - random number generator
   !>    rng_status          - keep track of rng status for restart purposes
   !>    thermostat_energy   - energy difference for conserved quantity
   !> \author Fabien Brieuc
   ! ***************************************************************************
   TYPE qtb_therm_type
      REAL(KIND=dp)                              :: tau = 0.0_dp, lamb = 0.0_dp
      REAL(KIND=dp)                              :: taucut = 0.0_dp, lambcut = 0.0_dp
      REAL(KIND=dp), DIMENSION(:), POINTER       :: c1 => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER       :: c2 => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER       :: g_fric => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: massfact => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: rf => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER    :: h => NULL()
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: r => NULL()
      INTEGER, DIMENSION(:), POINTER             :: step => NULL()
      INTEGER, DIMENSION(:), POINTER             :: cpt => NULL()
      INTEGER                                    :: fp = 0
      INTEGER                                    :: nf = 0
      REAL(KIND=dp)                              :: thermostat_energy = 0.0_dp
      TYPE(rng_stream_type)                      :: gaussian_rng_stream = rng_stream_type()
      CHARACTER(LEN=rng_record_length), DIMENSION(:), POINTER :: rng_status => NULL()
   END TYPE qtb_therm_type

   ! ***************************************************************************
   !> \brief data for the use of different Path Integral propagators
   !> \note
   !>    prop_kind     - selects a hamiltonian for the equations of motion
   !>    temp_sim2phys - conversion factor for simulation to physical temperature
   !>    temp_phys2sim - conversion factor for physical to simulation temperature
   !>    physpotscale  - factor to scale the physical interaction potential
   !> \author Felix Uhl
   ! ***************************************************************************
   TYPE pint_propagator_type
      INTEGER                              :: prop_kind = 0
      REAL(KIND=dp)                        :: temp_phys2sim = 0.0_dp
      REAL(KIND=dp)                        :: temp_sim2phys = 0.0_dp
      REAL(KIND=dp)                        :: physpotscale = 0.0_dp
   END TYPE pint_propagator_type

END MODULE pint_types
