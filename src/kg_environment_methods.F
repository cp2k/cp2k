!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004 CP2K developers group                           !
!-----------------------------------------------------------------------------!
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
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_type
  USE distribution_methods,            ONLY: distribute_molecules_2d
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type,&
                                             set_kg_env
  USE kg_gpw_pw_env_methods,           ONLY: update_rho_mol_blocks
  USE kg_gpw_pw_env_types,             ONLY: kg_sub_pw_env_type
  USE kg_rho_methods,                  ONLY: kg_rho_create
  USE kg_rho_types,                    ONLY: kg_rho_release,&
                                             kg_rho_type
  USE kg_rspw_methods,                 ONLY: kg_rspw_create,&
                                             kg_rspw_rebuild
  USE kg_rspw_types,                   ONLY: kg_rspw_release,&
                                             kg_rspw_type
  USE kg_scf,                          ONLY: kg_scf_env_did_change
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE qs_collocate_density,            ONLY: calculate_rho_core
  USE qs_environment_methods,          ONLY: qs_env_rebuild_rho
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_did_change
  USE qs_matrix_pools,                 ONLY: mpools_rebuild_s_sm_pools
  USE sparse_matrix_types,             ONLY: deallocate_matrix_set,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: kg_env_setup, kg_qs_env_update

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_environment_methods'
  
!-----------------------------------------------------------------------------!

CONTAINS

!***************************************************************************

  SUBROUTINE kg_env_setup(kg_env,kg_gpw,error)

    TYPE(kg_environment_type), POINTER       :: kg_env
    LOGICAL, INTENT(IN), OPTIONAL            :: kg_gpw
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_env_setup', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_kg_gpw
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

  failure=.FALSE.
  my_kg_gpw = .FALSE.
  IF(PRESENT(kg_gpw)) my_kg_gpw=kg_gpw

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

    IF(.NOT. my_kg_gpw) THEN
      !set up pw_pools, rs_pools
      CALL get_kg_env(kg_env=kg_env,rspw=rspw,error=error)
      IF (.NOT.ASSOCIATED(rspw)) THEN
         CALL kg_rspw_create(rspw,kg_env=kg_env,error=error)
         CALL set_kg_env(kg_env,rspw=rspw)
         CALL kg_rspw_release(rspw,force_env_section=kg_env%input,error=error)
      ELSE
         CALL kg_rspw_rebuild(rspw,kg_env=kg_env,error=error)
      END IF
  
     ! allocate rho structures
      CALL kg_rho_create(kg_rho=rho,kg_env=kg_env,error=error)
      CALL set_kg_env(kg_env,rho=rho)
      CALL kg_rho_release(kg_rho=rho,error=error)
    END IF

    ! create 2d distribution
    CALL distribute_molecules_2d(particle_kind_set=atomic_kind_set,&
                                  particle_set=particle_set,&
                                  molecule_kind_set=molecule_kind_set,&
                                  molecule_set=molecule_set,&
                                  distribution_2d=distribution_2d,&
                                  blacs_env=kg_env%blacs_env,&
                                  force_env_section=kg_env%input,&
                                  error=error)

    CALL set_kg_env(kg_env,distribution_2d=distribution_2d,error=error)
    CALL distribution_2d_release(distribution_2d, error=error)
  END IF

  END SUBROUTINE kg_env_setup

!******************************************************************************

  SUBROUTINE kg_qs_env_update(kg_env, qs_env, error)

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_qs_env_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nspins
    LOGICAL                                  :: failure, use_tau
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(kg_sub_pw_env_type), POINTER        :: kg_sub_pw_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks

!------------------------------------------------------------------------------  

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
 
    NULLIFY(para_env, cell, dft_control,particle_set)
  
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
  
    CALL get_qs_env(qs_env, cell=cell, para_env=para_env, dft_control=dft_control,&
         particle_set=particle_set,matrix_ks=matrix_ks, error=error)
  
    nspins = dft_control%nspins
    use_tau = dft_control%use_kinetic_energy_density
    ! *** rebuilds the S_sm_pools ***
    CALL mpools_rebuild_s_sm_pools(qs_env%mpools,s=qs_env%matrix_s,&
         nspins=nspins,particle_set=particle_set,&
         para_env=para_env, error=error)
  
    CALL get_qs_env(qs_env, matrix_ks=matrix_ks, error=error)
  
    ! *** deallocate matrices that will have the wrong structure ***
    CALL deallocate_matrix_set(matrix_ks)
    CALL set_qs_env(qs_env,matrix_ks=matrix_ks)
  
    ! *** updates rho core ***
    CPPrecondition(ASSOCIATED(qs_env%rho_core),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env%pw_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL calculate_rho_core(qs_env%rho_core, &
            qs_env%qs_charges%total_rho_core_rspace,qs_env,error=error)
    END IF
   
    ! *** tell ks_env ***
    IF (ASSOCIATED(qs_env%ks_env)) THEN
       CALL qs_ks_did_change(qs_env%ks_env,s_mstruct_changed=.TRUE., error=error)
    END IF

    !   *** Updates rho structure ***
    CALL qs_env_rebuild_rho(qs_env=qs_env)

    !   *** Updates the molecular rho structures ***
    CALL get_qs_env(qs_env=qs_env, kg_sub_pw_env=kg_sub_pw_env)
    CALL update_rho_mol_blocks(kg_sub_pw_env, particle_set, cell, nspins, use_tau, error=error)
 
    ! *** tell scf_env ***
    IF (ASSOCIATED(kg_env%scf_env)) THEN
       CALL kg_scf_env_did_change(kg_env%scf_env,&
            error=error)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE kg_qs_env_update

!******************************************************************************
END MODULE kg_environment_methods
