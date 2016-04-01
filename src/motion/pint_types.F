!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE pint_types

   USE cp_log_handling,                 ONLY: cp_logger_type
   USE gle_system_types,                ONLY: gle_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE parallel_rng_types,              ONLY: rng_stream_type
   USE replica_types,                   ONLY: replica_env_type
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
                                 thermostat_piglet = 4

   PUBLIC :: pint_env_type
   PUBLIC :: normalmode_env_type
   PUBLIC :: staging_env_type
   PUBLIC :: pile_therm_type
   PUBLIC :: piglet_therm_type
   PUBLIC :: pint_propagator_type

   ! ***************************************************************************
   !> \brief environment for a path integral run
   !> \param ref_count reference count of this data structure
   !> \param id_nr identification number of this data structure
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
   !> \param e_kin_beads (fictious) kinetic energy of the beads
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
   !> \param mass_fict fictious mass of the beads for dynamics (kinetic mass)
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
   ! ***************************************************************************
   TYPE pint_env_type
      INTEGER :: ref_count, id_nr, p, nnos, nrespa, iter, ndim, transform
      INTEGER :: first_step, last_step, num_steps, first_propagated_mode
      INTEGER :: pimd_thermostat, harm_integrator, thermostat_rng_seed
      REAL(KIND=dp) :: t_tol, v_tol, kT, beta, dt, &
                       e_gle, e_pile, e_piglet, e_pot_h, e_kin_beads, e_pot_t, e_kin_t, t, time_per_step
      TYPE(cp_logger_type), POINTER :: logger
      TYPE(replica_env_type), POINTER :: replicas
      TYPE(section_vals_type), POINTER :: input
      TYPE(staging_env_type), POINTER :: staging_env
      TYPE(normalmode_env_type), POINTER :: normalmode_env
      TYPE(rng_stream_type), POINTER :: randomG
      TYPE(gle_type), POINTER        :: gle
      REAL(KIND=dp), DIMENSION(e_num_ids) :: energy
      REAL(KIND=dp), DIMENSION(:), POINTER :: mass, e_pot_bead
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: x, v, f, mass_beads, &
                                                 mass_fict, ux, ux_t, uv, uv_t, uv_new, uf, uf_h, external_f
      REAL(KIND=dp), DIMENSION(:), POINTER :: centroid
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: tx, tv, tv_t, tv_old, tv_new, tf
      REAL(KIND=dp), DIMENSION(:), POINTER :: Q ! dim p, make it (p,ndim)?
      REAL(KIND=dp), DIMENSION(:), POINTER :: rtmp_ndim, rtmp_natom
      REAL(KIND=dp), DIMENSION(:), POINTER :: iwsinex, wsinex, cosex
      TYPE(pile_therm_type), POINTER :: pile_therm
      TYPE(piglet_therm_type), POINTER :: piglet_therm
      TYPE(pint_propagator_type), POINTER :: propagator
   END TYPE pint_env_type

   ! ***************************************************************************
   !> \brief data to perform the normalmode transformation
   !> \note
   !>    ref_count     - reference count of this data structure
   !>    id_nr         - identification number of this data structure
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
      INTEGER :: id_nr, ref_count, p
      REAL(KIND=dp) :: Q_bead, Q_centroid, modefactor, harm
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: x2u, u2x
      REAL(KIND=dp), DIMENSION(:), POINTER :: lambda
   END TYPE normalmode_env_type

   ! ***************************************************************************
   !> \brief data to perform the staging transformation
   !> \note
   !>    ref_count     - reference count of this data structure
   !>    id_nr         - identification number of this data structure
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
      INTEGER :: id_nr, ref_count, nseg, j, p
      REAL(KIND=dp) w_p, w_j, Q_stage, Q_end
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
      INTEGER                              :: ref_count
      REAL(KIND=dp)                        :: lamb, tau, thermostat_energy
      REAL(KIND=dp), DIMENSION(:), POINTER :: c1
      REAL(KIND=dp), DIMENSION(:), POINTER :: c2
      REAL(KIND=dp), DIMENSION(:), POINTER :: g_fric
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: massfact
      TYPE(rng_stream_type), POINTER       :: gaussian_rng_stream
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
      INTEGER                              :: ref_count
      INTEGER                              :: ndim, p, nsp1
      REAL(KIND=dp)                        :: thermostat_energy
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: a_mat, c_mat
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: gle_s, gle_t
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: smalls
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: temp1
      REAL(KIND=dp), DIMENSION(:), POINTER   :: temp2
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: sqrtmass
      TYPE(rng_stream_type), POINTER       :: gaussian_rng_stream
   END TYPE piglet_therm_type

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
      INTEGER                              :: prop_kind
      REAL(KIND=dp)                        :: temp_phys2sim
      REAL(KIND=dp)                        :: temp_sim2phys
      REAL(KIND=dp)                        :: physpotscale
   END TYPE pint_propagator_type

END MODULE pint_types
