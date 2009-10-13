!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Set of routines to apply restraints to the KS hamiltonian
! *****************************************************************************
MODULE qs_ks_apply_restraints

  USE cp_control_types,                ONLY: becke_restraint_type,&
                                             dft_control_type
  USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                             dbcsr_deallocate_matrix_set,&
                                             dbcsr_from_sm,&
                                             dbcsr_obj_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_methods,                   ONLY: dbcsr_init
  USE dbcsr_operations,                ONLY: dbcsr_copy,&
                                             dbcsr_set
  USE et_coupling,                     ONLY: becke_restraint
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mulliken,                        ONLY: mulliken_restraint
  USE pw_methods,                      ONLY: pw_scale
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_mo_types,                     ONLY: mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE s_square_methods,                ONLY: s2_restraint
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_apply_restraints'
  INTEGER, SAVE, PRIVATE :: last_ks_id_nr=0

  PUBLIC :: qs_ks_becke_restraint, qs_ks_mulliken_restraint, qs_ks_s2_restraint

CONTAINS

  SUBROUTINE qs_ks_becke_restraint(qs_env,auxbas_pw_pool,calculate_forces,matrix_s,becke,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(becke_restraint_type), POINTER      :: becke
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_becke_restraint', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: inv_vol

    IF(qs_env%dft_control%qs_control%becke_restraint)THEN
    !***** Check if becke potential is needed to constrain charges *****
       becke => qs_env%dft_control%qs_control%becke_control
       IF(becke%need_pot.OR.calculate_forces)THEN
          CALL pw_pool_create_pw(auxbas_pw_pool,becke%becke_pot%pw,use_data=REALDATA3D,&
               in_space=REALSPACE,error=error) 
          CALL becke_restraint(qs_env,becke_const=becke%becke_pot,calc_pot=.TRUE.,&
               calculate_forces=calculate_forces,error=error)
          CALL pw_scale(becke%becke_pot%pw,becke%becke_pot%pw%pw_grid%dvol,error=error)
          becke%need_pot=.FALSE.
       ELSE
          inv_vol=1.0_dp/becke%becke_pot%pw%pw_grid%dvol
          CALL pw_scale(becke%becke_pot%pw,inv_vol,error=error)
          CALL becke_restraint(qs_env,becke%becke_pot,calc_pot=.FALSE.,&
               calculate_forces=calculate_forces,error=error)
          CALL pw_scale(becke%becke_pot%pw,becke%becke_pot%pw%pw_grid%dvol,error=error)
       ENDIF

       IF(qs_env%dft_control%qs_control%et_coupling_calc)THEN
          IF(qs_env%et_coupling%keep_matrix)THEN
             IF(qs_env%et_coupling%first_run)THEN
                NULLIFY(qs_env%et_coupling%rest_mat(1)%matrix)                  
                ALLOCATE(qs_env%et_coupling%rest_mat(1)%matrix)
                CALL dbcsr_init(qs_env%et_coupling%rest_mat(1)%matrix)
                CALL dbcsr_copy(qs_env%et_coupling%rest_mat(1)%matrix, matrix_s(1)%matrix, &
                     name="ET_RESTRAINT_MATRIX_A", error=error)
                CALL dbcsr_set(qs_env%et_coupling%rest_mat(1)%matrix,0.0_dp,error=error)
                !CALL dbcsr_from_sm(qs_env%et_coupling%rest_mat(1)%matrix, matrix_s(1)%matrix, error)
                !CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                !     TARGET=qs_env%et_coupling%rest_mat(1)%matrix,&
                !     target_name="ET_RESTRAINT_MATRIX_A",error=error)
                CALL integrate_v_rspace(becke%becke_pot,&
                     h=qs_env%et_coupling%rest_mat(1),&
                     qs_env=qs_env,calculate_forces=.FALSE.,error=error)
                qs_env%et_coupling%order_p=qs_env%dft_control%qs_control%becke_control%becke_order_p
                qs_env%et_coupling%e1=qs_env%dft_control%qs_control%becke_control%strength
                qs_env%et_coupling%keep_matrix=.FALSE.
             ELSE
                NULLIFY(qs_env%et_coupling%rest_mat(2)%matrix)                  
                ALLOCATE(qs_env%et_coupling%rest_mat(2)%matrix)
                CALL dbcsr_init(qs_env%et_coupling%rest_mat(2)%matrix)
                CALL dbcsr_copy(qs_env%et_coupling%rest_mat(2)%matrix, matrix_s(1)%matrix, &
                     name="ET_RESTRAINT_MATRIX_B", error=error)
                CALL dbcsr_set(qs_env%et_coupling%rest_mat(2)%matrix,0.0_dp,error=error)
                !CALL dbcsr_from_sm(qs_env%et_coupling%rest_mat(2)%matrix, matrix_s(1)%matrix, error)
                !CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                !     TARGET=qs_env%et_coupling%rest_mat(2)%matrix,&
                !     target_name="ET_RESTRAINT_MATRIX_B",error=error)
                CALL integrate_v_rspace(becke%becke_pot,&
                     h=qs_env%et_coupling%rest_mat(2),&
                     qs_env=qs_env,calculate_forces=.FALSE.,error=error)

             END IF
          END IF
       END IF
    END IF
  END SUBROUTINE qs_ks_becke_restraint


  SUBROUTINE qs_ks_mulliken_restraint(energy,dft_control,just_energy,para_env,&
       ks_matrix, matrix_s,rho,mulliken_order_p,error)

    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(dft_control_type), POINTER          :: dft_control
    LOGICAL, INTENT(in)                      :: just_energy
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_s
    TYPE(qs_rho_type), POINTER               :: rho
    REAL(KIND=dp)                            :: mulliken_order_p
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_mulliken_restraint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao_sm

    energy%mulliken = 0.0_dp
    IF (dft_control%qs_control%mulliken_restraint) THEN

    rho_ao_sm=>rho%rho_ao
    NULLIFY(rho_ao)
    CALL dbcsr_allocate_matrix_set( rho_ao, SIZE(rho_ao_sm), error )!sm->dbcsr
    DO ispin=1,SIZE(rho_ao)!sm->dbcsr
       ALLOCATE(rho_ao(ispin)%matrix)!sm->dbcsr
       CALL dbcsr_init(rho_ao(ispin)%matrix)
       CALL dbcsr_from_sm(rho_ao(ispin)%matrix, rho_ao_sm(ispin)%matrix, error)!sm->dbcsr
    ENDDO!sm->dbcsr

       IF (just_energy) THEN
          CALL mulliken_restraint(dft_control%qs_control%mulliken_restraint_control, &
               para_env,matrix_s(1)%matrix, rho_ao,energy=energy%mulliken, &
               order_p=mulliken_order_p,error=error)
       ELSE
          CALL mulliken_restraint(dft_control%qs_control%mulliken_restraint_control,&
               para_env,matrix_s(1)%matrix,rho_ao,energy=energy%mulliken, &
               ks_matrix=ks_matrix, order_p=mulliken_order_p,error=error)
       ENDIF

   CALL dbcsr_deallocate_matrix_set( rho_ao, error=error )!sm->dbcsr

    ENDIF

  END SUBROUTINE qs_ks_mulliken_restraint


  SUBROUTINE qs_ks_s2_restraint(dft_control,qs_env,matrix_s,&
       mo_array,energy,calculate_forces, just_energy,error)

    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_s2_restraint', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs

    IF (dft_control%qs_control%s2_restraint) THEN
    ! adds s2_restraint energy and orbital derivatives
       CPPrecondition(dft_control%nspins == 2,cp_failure_level,routineP,error,failure)
       CPPrecondition(qs_env%requires_mo_derivs,cp_failure_level,routineP,error,failure)
       ! forces are not implemented (not difficult, but ... )
       CPPrecondition(.NOT. calculate_forces, cp_failure_level,routineP,error,failure)
       CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array,error=error)
       CALL s2_restraint(mo_array, matrix_s, mo_derivs, energy%s2_restraint, &
            dft_control%qs_control%s2_restraint_control, just_energy, error)
    ELSE
       energy%s2_restraint = 0.0_dp
    ENDIF
  END SUBROUTINE qs_ks_s2_restraint

END MODULE qs_ks_apply_restraints
