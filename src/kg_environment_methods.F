!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004 CP2K developers group                           !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/kg_environment_methods [1.0] *
!!
!!   NAME
!!     kg_environment_methods
!!
!!   FUNCTION
!!    to complete the setup of the kg_env      
!!   AUTHOR
!!     gt 15-10-2004    
!!   HISTORY 
!!*****
!******************************************************************************

MODULE kg_environment_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_type
  USE distribution_methods,            ONLY: distribute_molecules_2d
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type,&
                                             set_kg_env
  USE kg_rho_methods,                  ONLY: kg_rho_create
  USE kg_rho_types,                    ONLY: kg_rho_release,&
                                             kg_rho_type
  USE kg_rspw_methods,                 ONLY: kg_rspw_create,&
                                             kg_rspw_rebuild
  USE kg_rspw_types,                   ONLY: kg_rspw_release,&
                                             kg_rspw_type 
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE orbital_pointers,                ONLY: init_orbital_pointers
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
  USE particle_types,                  ONLY: particle_type  
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: topology_control

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: kg_env_setup

  CHARACTER(LEN=*), PARAMETER :: moduleN = "kg_environment_methods"
  
!-----------------------------------------------------------------------------!

CONTAINS

!***************************************************************************

  SUBROUTINE kg_env_setup(kg_env,error)

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_env_setup', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(kg_rho_type), POINTER               :: rho
    TYPE(kg_rspw_type), POINTER              :: rspw
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!***************************************************************************
  failure=.FALSE.

  NULLIFY(rho,rspw)
  NULLIFY(atomic_kind_set, distribution_2d,&
          molecule_kind_set, molecule_set, particle_set)


  CALL get_kg_env(kg_env=kg_env,&
                  atomic_kind_set=atomic_kind_set,&
                  molecule_kind_set=molecule_kind_set,&
                  molecule_set=molecule_set,&
                  particle_set=particle_set,&
                  rho=rho, rspw=rspw)

  CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)   
  IF (.NOT. failure) THEN

    !set up pw_pools, rs_pools
    CALL get_kg_env(kg_env=kg_env,rspw=rspw,error=error)
    IF (.NOT.ASSOCIATED(rspw)) THEN
       CALL kg_rspw_create(rspw,kg_env=kg_env,error=error)
       CALL set_kg_env(kg_env,rspw=rspw)
       CALL kg_rspw_release(rspw,error=error)
    ELSE
       CALL kg_rspw_rebuild(rspw,kg_env=kg_env,error=error)
    END IF

   ! allocate rho structures
    CALL kg_rho_create(kg_rho=rho,kg_env=kg_env,error=error)
    CALL set_kg_env(kg_env,rho=rho)
    CALL kg_rho_release(kg_rho=rho,error=error)
     ! create 2d distribution

    CALL distribute_molecules_2d(particle_kind_set=atomic_kind_set,&
                                  particle_set=particle_set,&
                                  molecule_kind_set=molecule_kind_set,&
                                  molecule_set=molecule_set,&
                                  distribution_2d=distribution_2d,&
                                  blacs_env=kg_env%blacs_env,&
                                  error=error)

    CALL set_kg_env(kg_env,distribution_2d=distribution_2d,error=error)
    CALL distribution_2d_release(distribution_2d, error=error)
  END IF

  END SUBROUTINE kg_env_setup

!******************************************************************************
END MODULE kg_environment_methods
