!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE pint_types

  USE gle_system_types,                ONLY: gle_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE parallel_rng_types,              ONLY: rng_stream_type
  USE replica_types,                   ONLY: replica_env_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! Energy contributions - symbolic names for indexing energy arrays
  INTEGER, PARAMETER, PUBLIC :: e_conserved_id  = 1, &
                                e_potential_id  = 2, &
                                e_kin_thermo_id = 3, &
                                e_kin_virial_id = 4

  ! Number of energy contributions for static array allocation
  INTEGER, PARAMETER, PUBLIC :: e_num_ids = 4

  INTEGER, PARAMETER, PUBLIC :: thermostat_none  = 0, &
                                thermostat_nose  = 1, &
                                thermostat_gle   = 2

  PUBLIC :: pint_env_type
  PUBLIC :: normalmode_env_type
  PUBLIC :: staging_env_type

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
  !> \param kT boltzmann factor times temperature
  !> \param beta 1/kT
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
  !> \param uf_h array with harmonic part transformed forces (sometimes with uf added)
  !> \param tx nose hoover chain positions (pint_env%nnos,pint_env%p,pint_env%ndim)
  !> \param tv nose hoover chain velocities
  !> \param tv_t nose hoover chain velocities (temporary)
  !> \param tv_old nose hoover chain velocities (older)
  !> \param tv_new nose hoover chain velocities (newer)
  !> \param tf nose hoover chain forces (?)
  !> \param Q nose hoover chain masses
  !> \param time_per_step - time per step in seconds (updated every step)
  !> \param first_propagated_mode if 1 - propagate all normal modes,
  !>                              if 2 - keep centoid fixed
  !> \author fawzi
  !> \par History
  !>      Added some comments - hforbert
  !>      Added normal mode transformation - hforbert
  !>      2009-06-15 helium_solvent_type object is no longer a member of
  !>                   pint_env_type [lwalewski]
  ! ***************************************************************************
  TYPE pint_env_type
     INTEGER :: ref_count, id_nr, p, nnos, nrespa, iter, ndim, transform
     INTEGER :: first_step, last_step, num_steps, first_propagated_mode
     INTEGER :: pimd_thermostat
     REAL(kind=dp) :: t_tol, v_tol, kT, beta, dt, &
          e_gle, e_pot_h, e_kin_beads, e_pot_t, e_kin_t, t, time_per_step
     TYPE(cp_logger_type), POINTER :: logger
     TYPE(replica_env_type), POINTER :: replicas
     TYPE(section_vals_type), POINTER :: input
     TYPE(staging_env_type), POINTER :: staging_env
     TYPE(normalmode_env_type), POINTER :: normalmode_env
     TYPE(rng_stream_type), POINTER :: randomG
     TYPE(gle_type), POINTER        :: gle
     REAL(kind=dp), DIMENSION(e_num_ids) :: energy
     REAL(kind=dp), DIMENSION(:), POINTER :: mass,e_pot_bead
     REAL(kind=dp), DIMENSION(:,:), POINTER :: x,v,f,mass_beads,&
          mass_fict,ux,uv,uv_t,uv_new,uf,uf_h, external_f
     REAL(kind=dp), DIMENSION(:,:,:), POINTER :: tx,tv,tv_t,tv_old,tv_new,tf
     REAL(kind=dp), DIMENSION(:), POINTER :: Q ! dim p, make it (p,ndim)?
     REAL(kind=dp), DIMENSION(:), POINTER :: rtmp_ndim, rtmp_natom
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
  !>
  !>    This could be done via FFT calls as well, but for now...
  !> \author hforbert
  ! ***************************************************************************
  TYPE normalmode_env_type
     INTEGER :: id_nr,ref_count,p
     REAL(kind=dp) :: Q_bead, Q_centroid, modefactor, harm
     REAL(kind=dp), DIMENSION(:,:), POINTER :: x2u, u2x
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
     INTEGER :: id_nr,ref_count,nseg,j,p
     REAL(kind=dp) w_p,w_j,Q_stage,Q_end
  END TYPE staging_env_type

END MODULE pint_types
