!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Debug routines for multipoles 
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
!> \date   05.2008
! *****************************************************************************
INTERFACE
  SUBROUTINE debug_ewald_multipoles(ewald_env, ewald_pw, nonbond_env, cell, &
       particle_set, local_particles, iw, debug_r_space, error )
    USE cell_types,                      ONLY: cell_type
    USE distribution_1d_types,           ONLY: distribution_1d_type
    USE ewald_environment_types,         ONLY: ewald_environment_type
    USE ewald_pw_types,                  ONLY: ewald_pw_type
    USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type
    USE kinds,                           ONLY: dp
    USE parallel_rng_types,              ONLY: UNIFORM,&
                                               create_rng_stream,&
                                               delete_rng_stream,&
                                               next_random_number,&
                                               rng_stream_type
    USE particle_types,                  ONLY: particle_type
#include "cp_common_uses.h"
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    INTEGER, INTENT(IN)                      :: iw
    LOGICAL, INTENT(IN)                      :: debug_r_space
    TYPE(cp_error_type), INTENT(inout)       :: error

  END SUBROUTINE debug_ewald_multipoles
END INTERFACE
