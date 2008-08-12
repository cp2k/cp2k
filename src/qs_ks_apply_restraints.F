!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that build the Kohn-Sham matrix (i.e calculate the coulomb
!>      and xc parts
!> \par History
!>      05.2002 moved from qs_scf (see there the history) [fawzi]
!>      JGH [30.08.02] multi-grid arrays independent from density and potential
!>      10.2002 introduced pools, uses updated rho as input,
!>              removed most temporary variables, renamed may vars,
!>              began conversion to LSD [fawzi]
!>      10.2004 moved calculate_w_matrix here [Joost VandeVondele]
!>              introduced energy derivative wrt MOs [Joost VandeVondele]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_ks_apply_restraints

  USE cp_control_types,                ONLY: becke_constraint_type,&
                                             ddapc_restraint_type,&
                                             dft_control_type
  USE cp_ddapc,                        ONLY: cp_ddapc_apply_CD,&
                                             cp_ddapc_apply_RF,&
                                             cp_ddapc_apply_RS
  USE cp_ddapc_forces,                 ONLY: reset_ch_pulay
  USE cp_fm_types,                     ONLY: cp_fm_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE et_coupling,                     ONLY: becke_constraint
  USE f77_blas
  USE input_constants,                 ONLY: do_spin_density
  USE kinds,                           ONLY: dp
  USE mulliken,                        ONLY: mulliken_restraint
  USE pw_methods,                      ONLY: pw_scale,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE s_square_methods,                ONLY: s2_restraint
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             replicate_matrix_structure
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_apply_restraints'
  INTEGER, SAVE, PRIVATE :: last_ks_id_nr=0

  PUBLIC :: qs_ks_ddapc_restraint, qs_ks_becke_restraint, qs_ks_mulliken_restraint,&
            qs_ks_s2_restraint
  !***

CONTAINS

  ! *****************************************************************************

  ! *****************************************************************************
  SUBROUTINE qs_ks_ddapc_restraint(qs_env,auxbas_pw_pool,rho_tot_gspace,v_hartree_gspace,&
       v_spin_ddapc_rest_r,energy,calculate_forces, ks_env, ks_matrix,just_energy,&
       ddapc_restraint_is_spin,explicit_potential,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_p_type)                          :: rho_tot_gspace, &
                                                v_hartree_gspace, &
                                                v_spin_ddapc_rest_r
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    LOGICAL, INTENT(in)                      :: just_energy
    LOGICAL                                  :: ddapc_restraint_is_spin, &
                                                explicit_potential
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_ddapc_restraint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ddapc_size, handle, i, my_id
    LOGICAL                                  :: et_coupling_calc, failure
    REAL(KIND=dp)                            :: ddapc_order_p
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(ddapc_restraint_type), POINTER      :: ddapc_restraint_control
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_p_type)                          :: v_spin_ddapc_rest_g
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: rest_mat

!   -----------------------------------------------------------------------
! In case decouple periodic images and/or apply restraints to charges

    CALL timeset(routineN,"I","",handle)
    failure        = .FALSE.
    logger         => cp_error_get_logger(error)
    ddapc_restraint_is_spin=.FALSE.
    et_coupling_calc=.FALSE.
    IF(qs_env%dft_control%qs_control%ddapc_restraint)THEN
       ddapc_size=SIZE( qs_env%dft_control%qs_control%ddapc_restraint_control )
       IF(SIZE(energy%ddapc_restraint).NE.ddapc_size)THEN
          DEALLOCATE(energy%ddapc_restraint)
          ALLOCATE(energy%ddapc_restraint(ddapc_size))
       END IF

       DO i=1,SIZE( qs_env%dft_control%qs_control%ddapc_restraint_control )
          my_id= qs_env%dft_control%qs_control%ddapc_restraint_control(i)%ddapc_restraint_control%density_type
          IF (my_id==do_spin_density.OR.ddapc_restraint_is_spin) &
                                      ddapc_restraint_is_spin=.TRUE.
       END DO
       et_coupling_calc= qs_env%dft_control%qs_control%et_coupling_calc
    ENDIF

    explicit_potential=ddapc_restraint_is_spin.OR.et_coupling_calc
    IF (explicit_potential) THEN
       CALL pw_pool_create_pw(auxbas_pw_pool,&
            v_spin_ddapc_rest_g%pw, &
            use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE,error=error)
       CALL pw_zero(v_spin_ddapc_rest_g%pw)
       CALL pw_pool_create_pw(auxbas_pw_pool,&
            v_spin_ddapc_rest_r%pw, &
            use_data=REALDATA3D,in_space=REALSPACE,error=error)
    ENDIF

    IF (.NOT.failure) THEN
       IF (calculate_forces) CALL reset_ch_pulay(qs_env, error=error)
       ! Decoupling/Recoupling
       CALL cp_ddapc_apply_CD(qs_env, rho_tot_gspace, energy%hartree, v_hartree_gspace,&
            calculate_forces, Itype_of_density="FULL DENSITY", error=error)
       IF(qs_env%dft_control%qs_control%ddapc_restraint)THEN
          ! Restraints/Constraints
          DO i=1,ddapc_size
             
             NULLIFY(ddapc_restraint_control)
             ddapc_restraint_control=> &
                  qs_env%dft_control%qs_control%ddapc_restraint_control(i)%ddapc_restraint_control

             CALL cp_ddapc_apply_RS(qs_env, energy%hartree, energy%ddapc_restraint(i),&
                                    v_hartree_gspace, v_spin_ddapc_rest_g, &
                                    i,ddapc_restraint_control,calculate_forces, error)
          END DO
       END IF


       CALL cp_ddapc_apply_RF(qs_env, rho_tot_gspace, energy%hartree, v_hartree_gspace,&
            calculate_forces, Itype_of_density="FULL DENSITY", error=error)
    END IF



    ! CJM Copying the real-space Hartree potential to KS_ENV
    IF ((.NOT. just_energy).OR.et_coupling_calc) THEN
       CALL pw_transfer(v_hartree_gspace%pw, ks_env%v_hartree_rspace%pw)
       CALL pw_scale(ks_env%v_hartree_rspace%pw, ks_env%v_hartree_rspace%pw%pw_grid%dvol)
       IF (explicit_potential) THEN
          CALL pw_transfer(v_spin_ddapc_rest_g%pw,v_spin_ddapc_rest_r%pw)
          CALL pw_scale(v_spin_ddapc_rest_r%pw,v_spin_ddapc_rest_r%pw%pw_grid%dvol)
          IF(et_coupling_calc)THEN
             IF(qs_env%et_coupling%keep_matrix)THEN
                IF(qs_env%et_coupling%first_run)THEN
                   NULLIFY(qs_env%et_coupling%rest_mat(1)%matrix)                  
                   CALL replicate_matrix_structure(source=ks_matrix(1)%matrix,&
                        TARGET=qs_env%et_coupling%rest_mat(1)%matrix,&
                        target_name="ET_RESTRAINT_MATRIX_B",error=error)
                   CALL integrate_v_rspace(v_spin_ddapc_rest_r,&
                        h=qs_env%et_coupling%rest_mat(1),&
                        qs_env=qs_env,calculate_forces=.FALSE.,error=error)
                   qs_env%et_coupling%order_p=&
                        qs_env%dft_control%qs_control%ddapc_restraint_control(1)%ddapc_restraint_control%ddapc_order_p
                   qs_env%et_coupling%e1=qs_env%dft_control%qs_control%ddapc_restraint_control(1)%ddapc_restraint_control%strength
                   qs_env%et_coupling%keep_matrix=.FALSE.
                ELSE
                   NULLIFY(qs_env%et_coupling%rest_mat(2)%matrix)                  
                   CALL replicate_matrix_structure(source=ks_matrix(1)%matrix,&
                        TARGET=qs_env%et_coupling%rest_mat(2)%matrix,&
                        target_name="ET_RESTRAINT_MATRIX_B",error=error)
                   CALL integrate_v_rspace(v_spin_ddapc_rest_r,&
                        h=qs_env%et_coupling%rest_mat(2),&
                        qs_env=qs_env,calculate_forces=.FALSE.,error=error)
                END IF
             END IF
          END IF
       ENDIF
    ENDIF

    IF (explicit_potential) THEN
       CALL pw_pool_give_back_pw(auxbas_pw_pool,v_spin_ddapc_rest_g%pw,&
            error=error)
    ENDIF

    CALL pw_pool_give_back_pw(auxbas_pw_pool,v_hartree_gspace%pw,&
         error=error)
    CALL timestop(0.0_dp, handle)

  END SUBROUTINE qs_ks_ddapc_restraint


  !   -----------------------------------------------------------------------
  SUBROUTINE qs_ks_becke_restraint(qs_env,auxbas_pw_pool,calculate_forces,matrix_s,becke,error)


    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(becke_constraint_type), POINTER     :: becke
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_becke_restraint', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: et_coupling_calc
    REAL(KIND=dp)                            :: inv_vol
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: rest_mat

!   -----------------------------------------------------------------------
!***** Check if becke potential is needed to constrain charges *****

    IF(qs_env%dft_control%qs_control%becke_constraint)THEN
       becke => qs_env%dft_control%qs_control%becke_control
       IF(becke%need_pot.OR.calculate_forces)THEN
          CALL pw_pool_create_pw(auxbas_pw_pool,becke%becke_pot%pw,use_data=REALDATA3D,&
               in_space=REALSPACE,error=error) 
          CALL becke_constraint(qs_env,becke_const=becke%becke_pot,calc_pot=.TRUE.,&
               calculate_forces=calculate_forces,error=error)
          CALL pw_scale(becke%becke_pot%pw,becke%becke_pot%pw%pw_grid%dvol)
          becke%need_pot=.FALSE.
       ELSE
          inv_vol=1.0_dp/becke%becke_pot%pw%pw_grid%dvol
          CALL pw_scale(becke%becke_pot%pw,inv_vol)
          CALL becke_constraint(qs_env,becke%becke_pot,calc_pot=.FALSE.,&
               calculate_forces=calculate_forces,error=error)
          CALL pw_scale(becke%becke_pot%pw,becke%becke_pot%pw%pw_grid%dvol)
       ENDIF

       IF(qs_env%dft_control%qs_control%et_coupling_calc)THEN
          IF(qs_env%et_coupling%keep_matrix)THEN
             IF(qs_env%et_coupling%first_run)THEN
                NULLIFY(qs_env%et_coupling%rest_mat(1)%matrix)                  
                CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                     TARGET=qs_env%et_coupling%rest_mat(1)%matrix,&
                     target_name="ET_RESTRAINT_MATRIX_A",error=error)
                CALL integrate_v_rspace(becke%becke_pot,&
                     h=qs_env%et_coupling%rest_mat(1),&
                     qs_env=qs_env,calculate_forces=.FALSE.,error=error)
                qs_env%et_coupling%order_p=qs_env%dft_control%qs_control%becke_control%becke_order_p
                qs_env%et_coupling%e1=qs_env%dft_control%qs_control%becke_control%strength
                qs_env%et_coupling%keep_matrix=.FALSE.
             ELSE
                NULLIFY(qs_env%et_coupling%rest_mat(2)%matrix)                  
                CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                     TARGET=qs_env%et_coupling%rest_mat(2)%matrix,&
                     target_name="ET_RESTRAINT_MATRIX_B",error=error)
                CALL integrate_v_rspace(becke%becke_pot,&
                     h=qs_env%et_coupling%rest_mat(2),&
                     qs_env=qs_env,calculate_forces=.FALSE.,error=error)

             END IF
          END IF
       END IF
    END IF
  END SUBROUTINE qs_ks_becke_restraint

  !   -----------------------------------------------------------------------
  SUBROUTINE qs_ks_mulliken_restraint(energy,dft_control,just_energy,para_env,&
       ks_matrix, matrix_s,rho,mulliken_order_p,error)

    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(dft_control_type), POINTER          :: dft_control
    LOGICAL, INTENT(in)                      :: just_energy
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_s
    TYPE(qs_rho_type), POINTER               :: rho
    REAL(KIND=dp)                            :: mulliken_order_p
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_mulliken_restraint', &
      routineP = moduleN//':'//routineN

!   -----------------------------------------------------------------------

    energy%mulliken = 0.0_dp
    IF (dft_control%qs_control%mulliken_restraint) THEN
       IF (just_energy) THEN
          CALL mulliken_restraint(dft_control%qs_control%mulliken_restraint_control, &
               para_env,matrix_s(1)%matrix, rho%rho_ao,energy=energy%mulliken, &
               order_p=mulliken_order_p)
       ELSE
          CALL mulliken_restraint(dft_control%qs_control%mulliken_restraint_control,&
               para_env,matrix_s(1)%matrix,rho%rho_ao,energy=energy%mulliken, &
               ks_matrix=ks_matrix, order_p=mulliken_order_p)
       ENDIF
    ENDIF

  END SUBROUTINE qs_ks_mulliken_restraint

  !   -----------------------------------------------------------------------
  SUBROUTINE qs_ks_s2_restraint(dft_control,qs_env,matrix_s,&
       mo_array,energy,calculate_forces, just_energy,error)

    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_s2_restraint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nspins
    LOGICAL                                  :: failure
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

!   -----------------------------------------------------------------------
! adds s2_restraint energy and orbital derivatives

    IF (dft_control%qs_control%s2_restraint) THEN
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
