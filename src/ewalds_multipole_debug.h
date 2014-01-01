!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Interfaces for MULTIPOLES debug routines
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

INTERFACE
  SUBROUTINE debug_ewald_multipoles_fields(ewald_env, ewald_pw, nonbond_env, cell,&
       particle_set, local_particles, radii, charges, dipoles, quadrupoles, task, iw,&
       atomic_kind_set, force_env_section, error)
    USE atomic_kind_types,               ONLY: atomic_kind_type
    USE cell_types,                      ONLY: cell_type
    USE distribution_1d_types,           ONLY: distribution_1d_type
    USE ewald_environment_types,         ONLY: ewald_environment_type
    USE ewald_pw_types,                  ONLY: ewald_pw_type
    USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type
    USE input_section_types,             ONLY: section_vals_type
    USE kinds,                           ONLY: dp
    USE particle_types,                  ONLY: particle_type
#include "cp_common_uses.h"
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), POINTER             :: particle_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), DIMENSION(:), &
         POINTER, OPTIONAL                   :: radii, charges
    REAL(KIND=dp), DIMENSION(:, :), &
         POINTER, OPTIONAL                   :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
         POINTER, OPTIONAL                   :: quadrupoles
    LOGICAL, DIMENSION(3), INTENT(IN)        :: task
    INTEGER, INTENT(IN)                      :: iw
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

  END SUBROUTINE debug_ewald_multipoles_fields
END INTERFACE

INTERFACE
  SUBROUTINE debug_ewald_multipoles_fields2(ewald_env, ewald_pw, nonbond_env, cell,&
       particle_set, local_particles, radii, charges, dipoles, quadrupoles, task, iw,&
       error)
    USE cell_types,                      ONLY: cell_type
    USE distribution_1d_types,           ONLY: distribution_1d_type
    USE ewald_environment_types,         ONLY: ewald_environment_type
    USE ewald_pw_types,                  ONLY: ewald_pw_type
    USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type
    USE kinds,                           ONLY: dp
    USE particle_types,                  ONLY: particle_type
#include "cp_common_uses.h"
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), POINTER             :: particle_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), DIMENSION(:), &
         POINTER, OPTIONAL                   :: radii, charges
    REAL(KIND=dp), DIMENSION(:, :), &
         POINTER, OPTIONAL                   :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
         POINTER, OPTIONAL                   :: quadrupoles
    LOGICAL, DIMENSION(3), INTENT(IN)        :: task
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

  END SUBROUTINE debug_ewald_multipoles_fields2
END INTERFACE
