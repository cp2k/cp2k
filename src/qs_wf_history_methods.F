!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_wf_history_methods [1.0] *
!!
!!   NAME
!!     qs_wf_history_methods
!!
!!   FUNCTION
!!     Storage of past states of the qs_env.
!!     Methods to interpolate (or actually normally extrapolate) the
!!     new guess for density and wavefunctions.
!!
!!   NOTES
!!     Most of the last snapshot should actually be in qs_env, but taking
!!     advantage of it would make the programming much convoluted
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!     11.2003 Joost VandeVondele : Implemented Nth order PS extrapolation 
!!     02.2005 modified for KG_GPW [MI]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_wf_history_methods
  USE coefficient_types,               ONLY: coeff_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_check,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale,&
                                             cp_fm_scale_and_add
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type, &
                                             cp_fm_to_fm
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type,&
                                             sm_pool_create_matrix,&
                                             sm_pools_create_matrix_vect
  USE global_types,                    ONLY: global_environment_type
  USE kg_gpw_fm_mol_types,             ONLY: kg_fm_mol_set_create,&
                                             kg_fm_mol_set_release,&
                                             kg_fm_mol_set_type,&
                                             kg_gpw_fm_mol_to_fm_mol, &
                                             duplicate_kg_fm_p_type
  USE kg_gpw_wf_history,               ONLY: kg_gpw_fm_mol_linear,&
                                             kg_gpw_fm_mol_ps,&
                                             kg_gpw_prev_wf
  USE kg_gpw_pw_env_methods,           ONLY: kg_rho_update_rho_mol
  USE kg_gpw_pw_env_types,             ONLY: kg_sub_pw_env_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: binomial
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_copy
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_reorthogonalize_vectors,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho, &
                                             duplicate_rho_type, &
                                             copy_rho_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type, &
                                             qs_rho_release
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             qs_wf_snapshot_type,&
                                             wfi_get_snapshot,&
                                             wfi_linear_p_method_nr,&
                                             wfi_linear_ps_method_nr,&
                                             wfi_linear_wf_method_nr,&
                                             wfi_ps_method_nr,&
                                             wfi_use_guess_method_nr,&
                                             wfi_use_prev_p_method_nr,&
                                             wfi_use_prev_rho_r_method_nr,&
                                             wfi_use_prev_wf_method_nr,&
                                             wfi_frozen_method_nr,&
                                             wfi_aspc, &
                                             wfi_clear_history, &
                                             wfs_release
  USE sparse_matrix_types,             ONLY: copy_matrix,&
                                             cp_sm_scale_and_add,&
                                             deallocate_matrix,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type,&
                                             transfer_matrix, &
                                             replicate_matrix, &
                                             allocate_matrix_set
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_wf_history_methods'
  INTEGER, SAVE, PRIVATE :: last_wfs_id=0, last_wfi_id=0

  PUBLIC :: wfs_create, wfi_create, wfi_update, &
            wfi_extrapolate, wfi_get_method_label, &
            wfi_duplicate_history, wfi_copy_history, &
            wfs_duplicate_snapshot, wfs_copy_snapshot

!!***
!****************************************************************************
CONTAINS

!!****f* qs_wf_history_methods/wfs_create [1.0] *
!!
!!   NAME
!!     wfs_create
!!
!!   SYNOPSIS
!!     Subroutine wfs_create(snapshot, error)
!!       Type(qs_wf_snapshot_type), Pointer:: snapshot
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine wfs_create
!!
!!   FUNCTION
!!     allocates and initialize a wavefunction snapshot
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - snapshot: the snapshot to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!     02.2005 added wf_mol [MI]
!!
!!*** **********************************************************************
SUBROUTINE wfs_create(snapshot, error)
    TYPE(qs_wf_snapshot_type), POINTER       :: snapshot
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfs_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  ALLOCATE(snapshot, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     last_wfs_id=last_wfs_id+1
     snapshot%id_nr=last_wfs_id
     NULLIFY(snapshot%wf, snapshot%wf_mol, snapshot%rho_r, &
             snapshot%rho_g, snapshot%rho_ao,&
             snapshot%overlap, snapshot%rho_frozen)
     snapshot%dt=1.0_dp
     snapshot%ref_count=1
  END IF
END SUBROUTINE wfs_create
!***************************************************************************

!!****f* qs_wf_history_methods/wfs_duplicate_snapshot [1.0] *
!!
!!   NAME
!!     wfs_duplicate_snapshot
!!
!!
!!   FUNCTION
!!     duplicates a single qs_wf_snapshot_type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - input_snapshot: the snapshot to be duplicated
!!     - output_snapshot: the out coming duplicate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [TdK]
!!
!!*** **********************************************************************
  SUBROUTINE wfs_duplicate_snapshot(input_snapshot, output_snapshot, qs_env, error)
    TYPE(qs_wf_snapshot_type), POINTER       :: input_snapshot, &
                                                output_snapshot
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfs_duplicate_snapshot', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: handle, &
                                                i, &
                                                stat, &
                                                nspins
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(dft_control_type), POINTER          :: dft_control

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)
    failure = .FALSE.
    
    NULLIFY(auxbas_pw_pool, pw_env, dft_control)

    CPPrecondition(ASSOCIATED(input_snapshot), cp_failure_level, routineP, error, failure)
    CPPrecondition(input_snapshot%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, dft_control=dft_control, error=error)
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)

      nspins = dft_control%nspins

      IF (ASSOCIATED(output_snapshot)) THEN
        CALL wfs_release(input_snapshot, error=error)
      ELSE
        CALL wfs_create(output_snapshot, error=error)
      END IF

      ! wf
      IF (ASSOCIATED(input_snapshot%wf)) THEN
        DO i = 1,SIZE(input_snapshot%wf)
          CALL cp_fm_create(matrix=output_snapshot%wf(i)%matrix, &
                            matrix_struct=input_snapshot%wf(i)%matrix%matrix_struct, &
                            error=error)
          CALL cp_fm_to_fm(source=input_snapshot%wf(i)%matrix, &
                           destination=output_snapshot%wf(i)%matrix, &
                           error=error)
        END DO
      END IF

      ! wf_mol
      IF (ASSOCIATED(input_snapshot%wf_mol)) THEN
        CALL duplicate_kg_fm_p_type(input_kg_fm_p_type=input_snapshot%wf_mol, &
                                    output_kg_fm_p_type=output_snapshot%wf_mol, &
                                    error=error)
      END IF

      ! rho_r
      IF (ASSOCIATED(input_snapshot%rho_r)) THEN
        DO i = 1,SIZE(input_snapshot%rho_r)
          ALLOCATE(output_snapshot%rho_r(nspins), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          IF (.NOT. failure) THEN
            CALL pw_pool_create_pw(auxbas_pw_pool, output_snapshot%rho_r(i)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE, &
                                   error=error)
            output_snapshot%rho_r(i)%pw%cr3d(:,:,:) = input_snapshot%rho_r(i)%pw%cr3d(:,:,:)
            IF (ASSOCIATED(input_snapshot%rho_r(i)%pw%mass_cr)) THEN
              output_snapshot%rho_r(i)%pw%mass_cr(:) = input_snapshot%rho_r(i)%pw%mass_cr(:)
            END IF
          END IF
        END DO
      END IF

      ! rho_g
      IF (ASSOCIATED(input_snapshot%rho_g)) THEN
        DO i = 1,SIZE(input_snapshot%rho_g)
          ALLOCATE(output_snapshot%rho_g(nspins), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          IF (.NOT. failure) THEN
            CALL pw_pool_create_pw(auxbas_pw_pool, output_snapshot%rho_g(i)%pw, &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                                   error=error)
            output_snapshot%rho_g(i)%pw%cc(:) = input_snapshot%rho_g(i)%pw%cc(:)
            IF (ASSOCIATED(input_snapshot%rho_g(i)%pw%mass_cc)) THEN
              output_snapshot%rho_g(i)%pw%mass_cc(:) = input_snapshot%rho_g(i)%pw%mass_cc(:)
            END IF
          END IF
        END DO
      END IF

      ! rho_ao
      IF (ASSOCIATED(input_snapshot%rho_ao)) THEN
        CALL allocate_matrix_set(output_snapshot%rho_ao, SIZE(input_snapshot%rho_ao))
        DO i = 1,SIZE(input_snapshot%rho_ao)
          CALL replicate_matrix(source=input_snapshot%rho_ao(i)%matrix, &
                                TARGET=output_snapshot%rho_ao(i)%matrix, &
                                target_name="myDensityMatrix_for_Spin_"//&
                                TRIM(ADJUSTL(cp_to_string(i))), &
                                allocate_blocks=.TRUE.)
        END DO
      END IF

      ! overlap
      IF (ASSOCIATED(input_snapshot%overlap)) THEN
        CALL replicate_matrix(source=input_snapshot%overlap, &
                              TARGET=output_snapshot%overlap, &
                              target_name="myOverlapMatrix", &
                              allocate_blocks=.TRUE.)
      END IF

      ! rho_frozen
      IF (ASSOCIATED(input_snapshot%rho_frozen)) THEN
        CALL duplicate_rho_type(rho_input=input_snapshot%rho_frozen, &
                                rho_output=output_snapshot%rho_frozen, &
                                qs_env=qs_env, error=error)
      END IF

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE wfs_duplicate_snapshot
!***************************************************************************

!!****f* qs_wf_history_methods/wfs_copy_snapshot [1.0] *
!!
!!   NAME
!!     wfs_copy_snapshot
!!
!!
!!   FUNCTION
!!     copies a single qs_wf_snapshot_type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - input_snapshot: the snapshot to be copied
!!     - output_snapshot: the out coming copy
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [TdK]
!!
!!*** **********************************************************************
  SUBROUTINE wfs_copy_snapshot(input_snapshot, output_snapshot, qs_env, error)
    TYPE(qs_wf_snapshot_type), POINTER       :: input_snapshot, &
                                                output_snapshot
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfs_copy_snapshot', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: handle, &
                                                i, &
                                                stat, &
                                                nspins, &
                                                ispin
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(dft_control_type), POINTER          :: dft_control

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)
    failure = .FALSE.
    
    NULLIFY(auxbas_pw_pool, pw_env, dft_control)

    CPPrecondition(ASSOCIATED(input_snapshot), cp_failure_level, routineP, error, failure)
    CPPrecondition(input_snapshot%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, dft_control=dft_control, error=error)
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)

      nspins = dft_control%nspins

      ! wf
      IF (ASSOCIATED(input_snapshot%wf)) THEN
        DO i = 1,SIZE(input_snapshot%wf)
          IF (ASSOCIATED(input_snapshot%wf(i)%matrix)) THEN
            CALL cp_fm_to_fm(source=input_snapshot%wf(i)%matrix, &
                             destination=output_snapshot%wf(i)%matrix, &
                             error=error)
          END IF
        END DO
      END IF

      ! wf_mol
      IF (ASSOCIATED(input_snapshot%wf_mol)) THEN
        CALL duplicate_kg_fm_p_type(input_kg_fm_p_type=input_snapshot%wf_mol, &
                                    output_kg_fm_p_type=output_snapshot%wf_mol, &
                                    error=error)
      END IF

      ! rho_r
      IF (ASSOCIATED(input_snapshot%rho_r)) THEN
        DO i = 1,SIZE(input_snapshot%rho_r)
          IF (ASSOCIATED(input_snapshot%rho_r(i)%pw)) THEN
            !CALL pw_pool_create_pw(auxbas_pw_pool, output_snapshot%rho_r(i)%pw, &
            !                       use_data=REALDATA3D, in_space=REALSPACE, &
            !                       error=error)
            output_snapshot%rho_r(i)%pw%cr3d(:,:,:) = input_snapshot%rho_r(i)%pw%cr3d(:,:,:)
            IF (ASSOCIATED(input_snapshot%rho_r(i)%pw%mass_cr)) THEN
              output_snapshot%rho_r(i)%pw%mass_cr(:) = input_snapshot%rho_r(i)%pw%mass_cr(:)
            END IF
          END IF
        END DO
      END IF

      ! rho_g
      IF (ASSOCIATED(input_snapshot%rho_g)) THEN
        DO i = 1,SIZE(input_snapshot%rho_g)
          IF (ASSOCIATED(input_snapshot%rho_g(i)%pw)) THEN
            !CALL pw_pool_create_pw(auxbas_pw_pool, output_snapshot%rho_g(i)%pw, &
            !                       use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
            !                       error=error)
            output_snapshot%rho_g(i)%pw%cc(:) = input_snapshot%rho_g(i)%pw%cc(:)
            IF (ASSOCIATED(input_snapshot%rho_g(i)%pw%mass_cc)) THEN
              output_snapshot%rho_g(i)%pw%mass_cc(:) = input_snapshot%rho_g(i)%pw%mass_cc(:)
            END IF
          END IF
        END DO
      END IF

      ! rho_ao
      IF (ASSOCIATED(input_snapshot%rho_ao)) THEN
        !CALL allocate_matrix_set(output_snapshot%rho_ao, SIZE(input_snapshot%rho_ao))
        DO i = 1,SIZE(input_snapshot%rho_ao)
          IF (ASSOCIATED(input_snapshot%rho_ao(i)%matrix)) THEN
            CALL replicate_matrix(source=input_snapshot%rho_ao(i)%matrix, &
                                  TARGET=output_snapshot%rho_ao(i)%matrix, &
                                  target_name="myDensityMatrix_for_Spin_"//&
                                  TRIM(ADJUSTL(cp_to_string(i))), &
                                  allocate_blocks=.FALSE.)
          END IF
        END DO
      END IF

      ! overlap
      IF (ASSOCIATED(input_snapshot%overlap)) THEN
        CALL replicate_matrix(source=input_snapshot%overlap, &
                              TARGET=output_snapshot%overlap, &
                              target_name="myOverlapMatrix", &
                              allocate_blocks=.FALSE.)
      END IF

      ! rho_frozen
      IF (ASSOCIATED(input_snapshot%rho_frozen)) THEN
        CALL copy_rho_type(rho_input=input_snapshot%rho_frozen, &
                           rho_output=output_snapshot%rho_frozen, &
                           qs_env=qs_env, error=error)
      END IF

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE wfs_copy_snapshot
!***************************************************************************

!!****f* qs_wf_history_methods/wfs_update [1.0] *
!!
!!   NAME
!!     wfs_update
!!
!!   FUNCTION
!!     updates the given snapshot
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - snapshot: the snapshot to be updated
!!     - wf_history: the history
!!     - qs_env: the qs_env that should be snapshotted
!!     - dt:the time of the snapshot (wrt. to the previous snapshot)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!     02.2005 added kg_fm_mol_set for KG_GPW [MI]
!!
!!*** **********************************************************************
SUBROUTINE wfs_update(snapshot,wf_history,qs_env,dt,kg_fm_mol_set,error)
    TYPE(qs_wf_snapshot_type), POINTER       :: snapshot
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: dt
    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER, OPTIONAL                      :: kg_fm_mol_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfs_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icount, ispin, &
                                                nmolecule_kind,&
                                                nspins, stat
    LOGICAL                                  :: failure, kg_gpw
    TYPE(coeff_type), DIMENSION(:), POINTER  :: rho_g, rho_r
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_pools
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao, s_array

  CALL timeset(routineN,"I","",handle)

  failure=.FALSE.
  NULLIFY(pw_env, auxbas_pw_pool, ao_mo_pools, dft_control, mos, mo_coeff,&
       rho, rho_r,rho_g,rho_ao, S_sm_pools,s_array)
  CALL get_qs_env(qs_env, pw_env=pw_env,&
       dft_control=dft_control,rho=rho,matrix_s=s_array,&
       error=error)
  CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_pools, &
       S_sm_pools=S_sm_pools, error=error)
  CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)

  CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
  IF (.not.ASSOCIATED(snapshot)) THEN
     CALL wfs_create(snapshot,error=error)
     CALL cp_error_check(error,failure)
  END IF
  IF (.not.failure) THEN
     CPPrecondition(wf_history%ref_count>0,cp_failure_level,routineP,error,failure)
     CPPrecondition(snapshot%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     nspins=dft_control%nspins
     snapshot%dt=1.0_dp
     IF (PRESENT(dt)) snapshot%dt=dt
     kg_gpw = .FALSE.
     IF(PRESENT(kg_fm_mol_set)) THEN
       IF(ASSOCIATED(kg_fm_mol_set)) kg_gpw = .TRUE.
     END IF
     IF (wf_history%store_wf) THEN
        IF(kg_gpw) THEN
          nmolecule_kind = SIZE(kg_fm_mol_set,1)
          IF (.NOT.ASSOCIATED(snapshot%wf_mol)) THEN
            CALL kg_fm_mol_set_create(snapshot%wf_mol,nmolecule_kind,error)
          END IF

          CALL kg_gpw_fm_mol_to_fm_mol(kg_fm_mol_set,snapshot%wf_mol%kg_fm_mol_set,error)
        ELSE
          CALL get_qs_env(qs_env,mos=mos,error=error)
          IF (.NOT.ASSOCIATED(snapshot%wf)) THEN
             CALL fm_pools_create_fm_vect(ao_mo_pools,snapshot%wf,&
                  name="ws_snap"//TRIM(ADJUSTL(cp_to_string(snapshot%id_nr)))//&
                  "ws",error=error)
             CPPostcondition(nspins==SIZE(snapshot%wf),cp_failure_level,routineP,error,failure)
          END IF
          DO ispin=1,nspins
             CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
             CALL cp_fm_to_fm(mo_coeff,snapshot%wf(ispin)%matrix)
          END DO
        END IF
     ELSE IF (ASSOCIATED(snapshot%wf)) THEN
        CALL fm_pools_give_back_fm_vect(ao_mo_pools,snapshot%wf,&
             error=error)
     ELSE IF (ASSOCIATED(snapshot%wf_mol)) THEN
        CALL kg_fm_mol_set_release(snapshot%wf_mol,error)
     END IF

     IF (wf_history%store_rho_r) THEN
        CALL qs_rho_get(rho, rho_r=rho_r)
        CPPrecondition(ASSOCIATED(rho_r),cp_failure_level,routineP,error,failure)
        IF (.NOT.ASSOCIATED(snapshot%rho_r)) THEN
           ALLOCATE(snapshot%rho_r(nspins),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO ispin=1,nspins
              NULLIFY(snapshot%rho_r(ispin)%pw)
              CALL pw_pool_create_pw(auxbas_pw_pool,snapshot%rho_r(ispin)%pw,&
                   in_space=REALSPACE, use_data=REALDATA3D,error=error)
           END DO
        END IF
        DO ispin=1,nspins
           CALL pw_copy(rho_r(ispin)%pw,snapshot%rho_r(ispin)%pw)
        END DO
     ELSE IF (ASSOCIATED(snapshot%rho_r)) THEN
        DO ispin=1,SIZE(snapshot%rho_r)
           CALL pw_pool_give_back_pw(auxbas_pw_pool,snapshot%rho_r(ispin)%pw,&
                error=error)
        END DO
        DEALLOCATE(snapshot%rho_r,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF

     IF (wf_history%store_rho_g) THEN
        CALL qs_rho_get(rho, rho_g=rho_g)
        CPPrecondition(ASSOCIATED(rho_g),cp_failure_level,routineP,error,failure)
        IF (.NOT.ASSOCIATED(snapshot%rho_g)) THEN
           ALLOCATE(snapshot%rho_g(nspins),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO ispin=1,nspins
              NULLIFY(snapshot%rho_g(ispin)%pw)
              CALL pw_pool_create_pw(auxbas_pw_pool,snapshot%rho_g(ispin)%pw,&
                   in_space=RECIPROCALSPACE,use_data=COMPLEXDATA1D,&
                   error=error)
           END DO
        END IF
        DO ispin=1,nspins
           CALL pw_copy(rho_g(ispin)%pw,snapshot%rho_g(ispin)%pw)
        END DO
     ELSE IF (ASSOCIATED(snapshot%rho_g)) THEN
        DO ispin=1,SIZE(snapshot%rho_g)
           CALL pw_pool_give_back_pw(auxbas_pw_pool,snapshot%rho_g(ispin)%pw,&
                error=error)
        END DO
        DEALLOCATE(snapshot%rho_g,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF

     IF (ASSOCIATED(snapshot%rho_ao)) THEN ! the sparsity might be different
        ! (future struct:check)
        CALL deallocate_matrix_set(snapshot%rho_ao)
     END IF
     IF (wf_history%store_rho_ao) THEN
        CALL qs_rho_get(rho, rho_ao=rho_ao)
        CPPrecondition(ASSOCIATED(rho_ao),cp_failure_level,routineP,error,failure)
        CALL sm_pools_create_matrix_vect(S_sm_pools,&
             snapshot%rho_ao,error=error)
        DO ispin=1,nspins
           CALL copy_matrix(rho_ao(ispin)%matrix,&
                snapshot%rho_ao(ispin)%matrix)
        END DO
     END IF

     IF (ASSOCIATED(snapshot%overlap)) THEN ! the sparsity might be different
        ! (future struct:check)
        CALL deallocate_matrix(snapshot%overlap)
     END IF
     IF (wf_history%store_overlap) THEN
        CPPrecondition(ASSOCIATED(s_array),cp_failure_level,routineP,error,failure)
        CPPrecondition(ASSOCIATED(s_array(1)%matrix),cp_failure_level,routineP,error,failure)
        CALL sm_pool_create_matrix(S_sm_pools(1)%pool,&
             snapshot%overlap,"snapshot%overlap",error=error)
        CALL copy_matrix(s_array(1)%matrix,snapshot%overlap)
     END IF

     IF (wf_history%store_frozen_density) THEN
       ! do nothing
       ! CALL deallocate_matrix_set(snapshot%rho_frozen%rho_ao)
     END IF

  END IF
  CALL timestop(0.0_dp,handle)

END SUBROUTINE wfs_update
!***************************************************************************

!!****f* qs_wf_history_methods/wfi_create [1.0] *
!!
!!   NAME
!!     wfi_create
!!
!!   SYNOPSIS
!!     Subroutine wfi_create(wf_history, interpolation_method_nr, extrapolation_order, error)
!!       Type(qs_wf_history_type), Pointer:: wf_history
!!       Integer, Intent (IN):: interpolation_method_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine wfi_create
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - interpolation_method_nr: the tag of the method used for 
!!       the extrapolation of the intial density for the next md step
!!       (see qs_wf_history_types:wfi_*_method_nr)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE wfi_create(wf_history, interpolation_method_nr, extrapolation_order, error)
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    INTEGER, INTENT(in)                      :: interpolation_method_nr
    INTEGER, INTENT(in), OPTIONAL            :: extrapolation_order
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, my_extrapolation_order, &
                                                stat
    LOGICAL                                  :: failure

  my_extrapolation_order = 1
  IF (PRESENT(extrapolation_order)) my_extrapolation_order=extrapolation_order

  failure=.FALSE.
  
  ALLOCATE(wf_history, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     last_wfi_id=last_wfi_id+1
     wf_history%id_nr=last_wfi_id
     wf_history%ref_count=1
     wf_history%memory_depth=0
     wf_history%snapshot_count=0
     wf_history%last_state_index=1
     wf_history%store_wf=.FALSE.
     wf_history%store_rho_r=.FALSE.
     wf_history%store_rho_g=.FALSE.
     wf_history%store_rho_ao=.FALSE.
     wf_history%store_overlap=.FALSE.
     wf_history%store_frozen_density=.FALSE.
     NULLIFY(wf_history%past_states)

     wf_history%interpolation_method_nr=interpolation_method_nr

     SELECT CASE(wf_history%interpolation_method_nr)
     CASE(wfi_use_guess_method_nr)
        wf_history%memory_depth=0
     CASE(wfi_use_prev_wf_method_nr)
        wf_history%memory_depth=0
     CASE(wfi_use_prev_p_method_nr)
        wf_history%memory_depth=1
        wf_history%store_rho_ao=.TRUE.
     CASE(wfi_use_prev_rho_r_method_nr)
        wf_history%memory_depth=1
        wf_history%store_rho_ao=.TRUE.
     CASE(wfi_linear_wf_method_nr)
        wf_history%memory_depth=2
        wf_history%store_wf=.TRUE.
     CASE(wfi_linear_p_method_nr)
        wf_history%memory_depth=2
        wf_history%store_rho_ao=.TRUE.
     CASE(wfi_linear_ps_method_nr)
        wf_history%memory_depth=2
        wf_history%store_wf=.TRUE.
        wf_history%store_overlap=.TRUE.
     CASE(wfi_ps_method_nr)
        wf_history%memory_depth=my_extrapolation_order+1
        wf_history%store_wf=.TRUE.
        wf_history%store_overlap=.TRUE.
     CASE(wfi_frozen_method_nr)
        wf_history%memory_depth=1
        wf_history%store_frozen_density=.TRUE.
     CASE (wfi_aspc)
       wf_history%memory_depth = my_extrapolation_order + 2
       wf_history%store_wf = .TRUE.
       wf_history%store_overlap = .TRUE.
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"Unknown interpolation method: "//&
             TRIM(ADJUSTL(cp_to_string(interpolation_method_nr)))//" in "//&
             CPSourceFileRef,&
             error, failure)
        wf_history%interpolation_method_nr=wfi_use_prev_rho_r_method_nr
     END SELECT
     ALLOCATE(wf_history%past_states(wf_history%memory_depth),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF (.NOT. failure) THEN
        DO i=1,SIZE(wf_history%past_states)
           NULLIFY(wf_history%past_states(i)%snapshot)
        END DO
     END IF
  END IF
END SUBROUTINE wfi_create
!***************************************************************************

!!****f* qs_wf_history_methods/wfi_get_method_label [1.0] *
!!
!!   NAME
!!     wfi_get_method_label
!!
!!   SYNOPSIS
!!     Function wfi_get_method_label(method_nr, error) Result(res)
!!       Integer, Intent (IN):: method_nr
!!       Character(Len=30):: res
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function wfi_get_method_label
!!
!!   FUNCTION
!!     returns a string describing the interpolation method
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - wf_history: the wf history object you want information about
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
FUNCTION wfi_get_method_label(method_nr,error) RESULT(res)
    INTEGER, INTENT(in)                      :: method_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    CHARACTER(len=30)                        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_get_method_label', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  res="unknown"
  IF (.NOT. failure) THEN
     SELECT CASE(method_nr)
     CASE(wfi_use_prev_p_method_nr)
        res="previous_p"
     CASE(wfi_use_prev_wf_method_nr)
        res="previous_wf"
     CASE(wfi_use_prev_rho_r_method_nr)
        res="previous_rho_r"
     CASE(wfi_use_guess_method_nr)
        res="initial_guess"
     CASE(wfi_linear_wf_method_nr)
        res="mo linear"
     CASE(wfi_linear_p_method_nr)
        res="P linear"
     CASE(wfi_linear_ps_method_nr)
        res="PS linear"
     CASE(wfi_ps_method_nr)
        res="PS Nth order"
     CASE(wfi_frozen_method_nr)
        res="frozen density approximation"
     CASE(wfi_aspc)
       res = "ASPC"
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"Unknown interpolation method: "//&
             TRIM(ADJUSTL(cp_to_string(method_nr)))//&
             " in "//&
             CPSourceFileRef,&
             error, failure)
     END SELECT
  END IF
END FUNCTION wfi_get_method_label
!***************************************************************************

!!****f* qs_wf_history_methods/wfi_extrapolate [1.0] *
!!
!!   NAME
!!     wfi_extrapolate
!!
!!   SYNOPSIS
!!     Subroutine wfi_extrapolate(wf_history, qs_env, dt,&
!!         extrapolation_method_nr, orthogonal_wf, error)
!!       Type(qs_wf_history_type), Pointer:: wf_history
!!       Type(qs_environment_type), Pointer:: qs_env
!!       Real(KIND = dp), Intent (IN):: dt
!!       Integer, Intent (OUT), Optional:: extrapolation_method_nr
!!       Logical, Intent (OUT), Optional:: orthogonal_wf
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine wfi_extrapolate
!!
!!   FUNCTION
!!     calculates the new starting state for the scf for the next
!!     wf optimization
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - wf_history: the previous history needed to extrapolate
!!     - qs_env: the qs env with the latest result, and that will contain
!!       the new starting state
!!     - dt: the time at which to extrapolate (wrt. to the last snapshot)
!!     - extrapolation_method_nr: returns the extrapolation method used
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!     11.2003 Joost VandeVondele : Implemented Nth order PS extrapolation 
!!
!!*** **********************************************************************
SUBROUTINE wfi_extrapolate(globenv,wf_history, qs_env, dt, extrapolation_method_nr, &
     orthogonal_wf, kg_gpw, gapw, gapw_xc, kg_fm_mol_set, error)
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: globenv
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(in)                :: dt
    INTEGER, INTENT(out), OPTIONAL           :: extrapolation_method_nr
    LOGICAL, INTENT(out), OPTIONAL           :: orthogonal_wf
    LOGICAL, INTENT(IN), OPTIONAL            :: gapw
    LOGICAL, INTENT(IN), OPTIONAL            :: gapw_xc
    LOGICAL, INTENT(IN), OPTIONAL            :: kg_gpw
    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER, OPTIONAL                      :: kg_fm_mol_set

    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_extrapolate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: actual_extrapolation_method_nr&
                                                , handle, i, ispin, k, n, &
                                                nmo, nvec
    LOGICAL                                  :: failure, my_gapw, my_gapw_xc, my_kg_gpw,&
                                                my_orthogonal_wf, ionode
    REAL(KIND=dp)                            :: alpha, beta, t0, t1, t2
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct, &
                                                matrix_struct_new
    TYPE(cp_fm_type), POINTER                :: csc, mo_coeff, fm_tmp
    TYPE(kg_sub_pw_env_type), POINTER        :: kg_sub_pw_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho, rho_xc
    TYPE(qs_wf_snapshot_type), POINTER       :: t0_state, t1_state
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s

  NULLIFY(mos, ao_mo_fm_pools, t0_state, t1_state, matrix_s, mo_coeff, rho, rho_xc)
  failure=.FALSE.
  my_gapw = .FALSE. 
  my_gapw_xc = .FALSE. 
  IF(PRESENT(gapw)) my_gapw = gapw
  IF(PRESENT(gapw_xc)) my_gapw_xc = gapw_xc
  my_kg_gpw = .FALSE. 
  IF(PRESENT(kg_gpw)) my_kg_gpw = kg_gpw

  CALL timeset(routineN,"I","",handle)
  
  CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,error,failure)
  CPPrecondition(wf_history%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)  
  IF (.NOT. failure) THEN
     CALL get_qs_env(qs_env, mos=mos,&
          rho=rho,rho_xc=rho_xc,matrix_s=matrix_s,error=error)
     CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools,&
          error=error)

     ! chooses the method for this extrapolation
     IF (wf_history%snapshot_count<1) THEN
        actual_extrapolation_method_nr=wfi_use_guess_method_nr
     ELSE
        actual_extrapolation_method_nr=wf_history%interpolation_method_nr
     END IF
     SELECT CASE(actual_extrapolation_method_nr)
      CASE(wfi_linear_wf_method_nr)
         IF (wf_history%snapshot_count<2) THEN
           actual_extrapolation_method_nr=wfi_use_prev_wf_method_nr
        END IF
     CASE(wfi_linear_p_method_nr)
        IF (wf_history%snapshot_count<2) THEN
           actual_extrapolation_method_nr=wfi_use_prev_wf_method_nr
        END IF
     CASE(wfi_linear_ps_method_nr)
        IF (wf_history%snapshot_count<2) THEN
           actual_extrapolation_method_nr=wfi_use_prev_wf_method_nr
        END IF
     END SELECT

     IF (PRESENT(extrapolation_method_nr)) &
          extrapolation_method_nr=actual_extrapolation_method_nr
     my_orthogonal_wf=.FALSE.
     ionode = (qs_env%para_env%mepos == qs_env%para_env%source)

     SELECT CASE (actual_extrapolation_method_nr)
     CASE(wfi_frozen_method_nr)
        t0_state => wfi_get_snapshot(wf_history, index=1, error=error)
        CPPrecondition(ASSOCIATED(t0_state%rho_frozen),cp_failure_level,routineP,error,failure)
        !DO ispin=1,SIZE(t0_state%rho_ao)
        !   CALL transfer_matrix(source_matrix=t0_state%rho_ao(ispin)%matrix,&
        !        target_matrix=rho%rho_ao(ispin)%matrix, error=error)
        !END DO

        DO ispin=1,SIZE(t0_state%rho_frozen%rho_ao)
           CALL transfer_matrix(source_matrix=t0_state%rho_frozen%rho_ao(ispin)%matrix,&
                target_matrix=rho%rho_ao(ispin)%matrix, error=error)
        END DO
        !FM updating rho_ao directly with t0_state%rho_ao would have the
        !FM wrong matrix structure
        CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
        CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE., error=error)

        my_orthogonal_wf=.TRUE.
     CASE(wfi_use_prev_rho_r_method_nr)
        t0_state => wfi_get_snapshot(wf_history, index=1, error=error)
        CPPrecondition(ASSOCIATED(t0_state%rho_ao),cp_failure_level,routineP,error,failure)
        DO ispin=1,SIZE(t0_state%rho_ao)
           CALL transfer_matrix(source_matrix=t0_state%rho_ao(ispin)%matrix,&
                target_matrix=rho%rho_ao(ispin)%matrix, error=error)
        END DO
        rho%rho_g_valid=.TRUE.
        IF(my_kg_gpw) THEN
           ! ** Updates the molecular rho for each molecule independently ***
           CALL get_qs_env(qs_env=qs_env, kg_sub_pw_env=kg_sub_pw_env)
           CALL kg_rho_update_rho_mol(kg_sub_pw_env, qs_env, kg_fm_mol_set, error=error)
        END IF

        ! does nothing
     CASE(wfi_use_prev_wf_method_nr)
        my_orthogonal_wf=.TRUE.

        IF(my_kg_gpw) THEN
           CPPrecondition(PRESENT(kg_fm_mol_set),cp_failure_level,routineP,error,failure)
           CALL kg_gpw_prev_wf(kg_fm_mol_set,rho%rho_ao,matrix_s,error=error)
        ELSE
          DO ispin=1,SIZE(mos)
             CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,&
                  nmo=nmo)
             CALL qs_env_reorthogonalize_vectors(qs_env,&
                  v_matrix=mo_coeff,ispin=ispin,&
                  n_col=nmo, error=error)
             CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
                  density_matrix=rho%rho_ao(ispin)%matrix,error=error)
          END DO
        END IF
        CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw = my_gapw, error=error)
        IF(my_gapw) THEN
           CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
        ENDIF
        IF(my_gapw_xc) THEN
            CALL qs_rho_update_rho(rho_xc,qs_env=qs_env, gapw=my_gapw_xc, error=error)
            CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
        END IF

        CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)
        IF(my_kg_gpw) THEN
           ! ** Updates the molecular rho for each molecule independently ***
           CALL get_qs_env(qs_env=qs_env, kg_sub_pw_env=kg_sub_pw_env)
           CALL kg_rho_update_rho_mol( kg_sub_pw_env, qs_env, kg_fm_mol_set, error=error)
        END IF

     CASE(wfi_use_prev_p_method_nr)
        t0_state => wfi_get_snapshot(wf_history, index=1, error=error)
        CPPrecondition(ASSOCIATED(t0_state%rho_ao),cp_failure_level,routineP,error,failure)
        DO ispin=1,SIZE(t0_state%rho_ao)
           CALL transfer_matrix(source_matrix=t0_state%rho_ao(ispin)%matrix,&
                target_matrix=rho%rho_ao(ispin)%matrix, error=error)
        END DO
        !FM updating rho_ao directly with t0_state%rho_ao would have the
        !FM wrong matrix structure
        CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw = my_gapw, error=error)
        IF(my_gapw) THEN
           CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
        ENDIF
        IF(my_gapw_xc) THEN
            CALL qs_rho_update_rho(rho_xc,qs_env=qs_env, gapw=my_gapw_xc, error=error)
            CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
        END IF

        CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE., error=error)
        IF(my_kg_gpw) THEN
           ! ** Updates the molecular rho for each molecule independently ***
           CALL get_qs_env(qs_env=qs_env, kg_sub_pw_env=kg_sub_pw_env)
           CALL kg_rho_update_rho_mol( kg_sub_pw_env, qs_env, kg_fm_mol_set, error=error)
        END IF
     CASE(wfi_use_guess_method_nr)
        !FM more clean to do it here, but it
        !FM might need to read a file (restart) and thus globenv
        !FM I do not want globenv here, thus done by the caller
        !FM (btw. it also needs the eigensolver, and unless you relocate it
        !FM gives circular dependencies)
     CASE(wfi_linear_wf_method_nr)
        t0_state => wfi_get_snapshot(wf_history, index=2, error=error)
        t1_state => wfi_get_snapshot(wf_history, index=1, error=error)
        CPPrecondition(ASSOCIATED(t0_state),cp_warning_level,routineP,error,failure)
        CPPrecondition(ASSOCIATED(t1_state),cp_warning_level,routineP,error,failure)
        IF(.not. failure .AND. my_kg_gpw) THEN
           CPPrecondition(ASSOCIATED(t0_state%wf_mol),cp_warning_level,routineP,error,failure)
           CPPrecondition(ASSOCIATED(t1_state%wf_mol),cp_warning_level,routineP,error,failure)           
        ELSE IF (.not.failure) THEN
           CPPrecondition(ASSOCIATED(t0_state%wf),cp_warning_level,routineP,error,failure)
           CPPrecondition(ASSOCIATED(t1_state%wf),cp_warning_level,routineP,error,failure)           
        END IF
        IF (.NOT.failure) THEN
           my_orthogonal_wf=.TRUE.
           t0=0.0_dp
           t1=t1_state%dt
           t2=t1+dt
           IF(my_kg_gpw) THEN
             CPPrecondition(PRESENT(kg_fm_mol_set),cp_failure_level,routineP,error,failure)
             CALL kg_gpw_fm_mol_linear(alpha=0.0_dp,set_a=kg_fm_mol_set,&
                  beta=(t2-t0)/(t1-t0),set_b=t1_state%wf_mol%kg_fm_mol_set,error=error)
             CALL kg_gpw_fm_mol_linear(alpha=1.0_dp,set_a=kg_fm_mol_set,&
                  beta=(t1-t2)/(t1-t0),set_b=t0_state%wf_mol%kg_fm_mol_set,&
                  matrix_p=rho%rho_ao,matrix_s=matrix_s,error=error)
           ELSE
             DO ispin=1,SIZE(mos)
                CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,&
                     nmo=nmo)
                CALL cp_fm_scale_and_add(alpha=0.0_dp,&
                     matrix_a=mo_coeff,&
                     matrix_b=t1_state%wf(ispin)%matrix,&
                     beta=(t2-t0)/(t1-t0),&
                     error=error) ! this copy should be unnecessary
                CALL cp_fm_scale_and_add(alpha=1.0_dp,&
                     matrix_a=mo_coeff,&
                     beta=(t1-t2)/(t1-t0), matrix_b=t0_state%wf(ispin)%matrix,&
                     error=error)
                CALL qs_env_reorthogonalize_vectors(qs_env,&
                     v_matrix=mo_coeff,ispin=ispin,&
                     n_col=nmo, error=error)
                CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
                     density_matrix=rho%rho_ao(ispin)%matrix,error=error)
             END DO
           END IF
           CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw=my_gapw, error=error)
           IF(my_gapw) THEN
             CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
           ENDIF
           IF(my_gapw_xc) THEN
              CALL qs_rho_update_rho(rho_xc,qs_env=qs_env, gapw=my_gapw_xc, error=error)
              CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
           END IF

           CALL qs_ks_did_change(qs_env%ks_env,&
                rho_changed=.TRUE., error=error)
           IF(my_kg_gpw) THEN
             ! ** Updates the molecular rho for each molecule independently ***
             CALL get_qs_env(qs_env=qs_env, kg_sub_pw_env=kg_sub_pw_env)
             CALL kg_rho_update_rho_mol( kg_sub_pw_env, qs_env, kg_fm_mol_set, error=error)
           END IF
        END IF 
     CASE(wfi_linear_p_method_nr)
        t0_state => wfi_get_snapshot(wf_history, index=2, error=error)
        t1_state => wfi_get_snapshot(wf_history, index=1, error=error)
        CPPrecondition(ASSOCIATED(t0_state),cp_warning_level,routineP,error,failure)
        CPPrecondition(ASSOCIATED(t1_state),cp_warning_level,routineP,error,failure)
        IF (.not.failure) THEN
           CPPrecondition(ASSOCIATED(t0_state%rho_ao),cp_warning_level,routineP,error,failure)
           CPPrecondition(ASSOCIATED(t1_state%rho_ao),cp_warning_level,routineP,error,failure)           
        END IF
        IF (.not.failure) THEN
          t0=0.0_dp
          t1=t1_state%dt
          t2=t1+dt
          DO ispin=1,SIZE(rho%rho_ao)
!             CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
             CALL cp_sm_scale_and_add(alpha=0.0_dp,&
                  matrix_a=rho%rho_ao(ispin)%matrix,&
                  matrix_b=t1_state%rho_ao(ispin)%matrix,&
                  beta=(t2-t0)/(t1-t0),&
                  error=error) ! this copy should be unnecessary
             CALL cp_sm_scale_and_add(alpha=1.0_dp,&
                  matrix_a=rho%rho_ao(ispin)%matrix,&
                  beta=(t1-t2)/(t1-t0), &
                  matrix_b=t0_state%rho_ao(ispin)%matrix,&
                  error=error)
                ! purify ?
          END DO
          CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw=my_gapw, error=error)
          IF(my_gapw) THEN
              CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
          ENDIF
          IF(my_gapw_xc) THEN
             CALL qs_rho_update_rho(rho_xc,qs_env=qs_env, gapw=my_gapw_xc, error=error)
             CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
          END IF

          CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE., error=error)
          IF(my_kg_gpw) THEN
             ! ** Updates the molecular rho for each molecule independently ***
             CALL get_qs_env(qs_env=qs_env, kg_sub_pw_env=kg_sub_pw_env)
             CALL kg_rho_update_rho_mol( kg_sub_pw_env, qs_env, kg_fm_mol_set, error=error)
          END IF
          ! wf not calculated, extract with PSC renormalized?
          ! use wf_linear?
        END IF
     CASE(wfi_linear_ps_method_nr)
        t0_state => wfi_get_snapshot(wf_history, index=2, error=error)
        t1_state => wfi_get_snapshot(wf_history, index=1, error=error)
        CPPrecondition(ASSOCIATED(t0_state),cp_warning_level,routineP,error,failure)
        CPPrecondition(ASSOCIATED(t1_state),cp_warning_level,routineP,error,failure)
        IF (.not.failure) THEN
           CPPrecondition(ASSOCIATED(t0_state%wf),cp_warning_level,routineP,error,failure)
           CPPrecondition(ASSOCIATED(t1_state%wf),cp_warning_level,routineP,error,failure)           
           CPPrecondition(ASSOCIATED(t0_state%overlap),cp_warning_level,routineP,error,failure)
           CPPrecondition(ASSOCIATED(t1_state%overlap),cp_warning_level,routineP,error,failure)           
        END IF
        IF (.not.failure) THEN
           my_orthogonal_wf=.TRUE.
           ! use PS_2=2 PS_1-PS_0
           ! C_2 comes from using PS_2 as a projector acting on C_1
           IF(my_kg_gpw) THEN
             stop 'kg_gpw wfi_linear_ps_method_nr and not implemented yet' 
             CPPrecondition(PRESENT(kg_fm_mol_set),cp_failure_level,routineP,error,failure)
           ELSE
             DO ispin=1,SIZE(mos)
                NULLIFY(mo_coeff,matrix_struct,matrix_struct_new,csc)
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
                CALL cp_fm_get_info(mo_coeff,nrow_global=n,ncol_global=k, &
                     matrix_struct=matrix_struct)
                CALL cp_fm_struct_create(matrix_struct_new,template_fmstruct=matrix_struct, &
                                         nrow_global=k,ncol_global=k)
                CALL cp_fm_create(csc,matrix_struct_new)
                CALL cp_fm_struct_release(matrix_struct_new)

                CALL cp_sm_fm_multiply(t0_state%overlap,t1_state%wf(ispin)%matrix,mo_coeff, k)
                CALL cp_fm_gemm('T','N',k,k,n,1.0_dp,t0_state%wf(ispin)%matrix,mo_coeff,0.0_dp,csc)
                CALL cp_fm_gemm('N','N',n,k,k,1.0_dp,t0_state%wf(ispin)%matrix,csc,0.0_dp,mo_coeff)
                CALL cp_fm_release(csc)
                CALL cp_fm_scale_and_add(-1.0_dp,mo_coeff,2.0_dp,t1_state%wf(ispin)%matrix)
                CALL qs_env_reorthogonalize_vectors(qs_env,&
                     v_matrix=mo_coeff,ispin=ispin,&
                     n_col=k, error=error)
                CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
                     density_matrix=rho%rho_ao(ispin)%matrix,error=error)
             END DO
           END IF
           CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw=my_gapw, error=error)
           IF(my_gapw) THEN
              CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
           ENDIF
           IF(my_gapw_xc) THEN
              CALL qs_rho_update_rho(rho_xc,qs_env=qs_env, gapw=my_gapw_xc, error=error)
              CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
           END IF
           CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)
           IF(my_kg_gpw) THEN
             ! ** Updates the molecular rho for each molecule independently ***
             CALL get_qs_env(qs_env=qs_env, kg_sub_pw_env=kg_sub_pw_env)
             CALL kg_rho_update_rho_mol( kg_sub_pw_env, qs_env, kg_fm_mol_set, error=error)
           END IF
        END IF
    CASE(wfi_ps_method_nr)
        ! figure out the actual number of vectors to use in the extrapolation:
        nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
        CPPrecondition(nvec .GT. 0,cp_failure_level,routineP,error,failure)
        my_orthogonal_wf=.TRUE.
        IF(my_kg_gpw) THEN
           CPPrecondition(PRESENT(kg_fm_mol_set),cp_failure_level,routineP,error,failure)
           ! first the most recent
           t1_state => wfi_get_snapshot(wf_history, index=1, error=error)
           beta = REAL(nvec,dp)
           IF(nvec==1) THEN
             CALL kg_gpw_fm_mol_linear(alpha=0.0_dp, set_a=kg_fm_mol_set,&
                  beta=beta, set_b=t1_state%wf_mol%kg_fm_mol_set,&
                  matrix_p=rho%rho_ao, matrix_s=matrix_s, error=error)
           ELSE
             CALL kg_gpw_fm_mol_linear(alpha=0.0_dp, set_a=kg_fm_mol_set,&
                  beta=beta, set_b=t1_state%wf_mol%kg_fm_mol_set, error=error)

             DO i=2,nvec-1
               t0_state => wfi_get_snapshot(wf_history, index=i, error=error)
               beta = -1.0_dp * beta * REAL(nvec - i + 1 , dp ) / REAL(i , dp )
               CALL kg_gpw_fm_mol_ps(set_a=kg_fm_mol_set,&
                    beta=beta,set_b=t0_state%wf_mol%kg_fm_mol_set,ov_b=t0_state%overlap,&
                    set_last=t1_state%wf_mol%kg_fm_mol_set,error=error)
             END DO
             ! finally calculate the density
             t0_state => wfi_get_snapshot(wf_history, index=nvec, error=error)
             beta = -1.0_dp * beta  / REAL(nvec , dp ) 
             CALL kg_gpw_fm_mol_ps(set_a=kg_fm_mol_set,&
                  beta=beta,set_b=t0_state%wf_mol%kg_fm_mol_set,ov_b=t0_state%overlap,&
                  set_last=t1_state%wf_mol%kg_fm_mol_set, matrix_p=rho%rho_ao,matrix_s=matrix_s,&
                  error=error)
           END IF
        ELSE
           DO ispin=1,SIZE(mos)
              NULLIFY(mo_coeff,matrix_struct,matrix_struct_new,csc,fm_tmp)
              CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
              CALL cp_fm_get_info(mo_coeff,nrow_global=n,ncol_global=k, &
                   matrix_struct=matrix_struct)
              CALL cp_fm_create(fm_tmp,matrix_struct)
              CALL cp_fm_struct_create(matrix_struct_new,template_fmstruct=matrix_struct, &
                                       nrow_global=k,ncol_global=k)
              CALL cp_fm_create(csc,matrix_struct_new)
              CALL cp_fm_struct_release(matrix_struct_new)
              ! first the most recent
              t1_state => wfi_get_snapshot(wf_history, index=1, error=error)
              CALL cp_fm_to_fm(t1_state%wf(ispin)%matrix,mo_coeff)
              alpha = nvec
              CALL cp_fm_scale(alpha,mo_coeff)
              DO i=2,nvec 
                 t0_state => wfi_get_snapshot(wf_history, index=i, error=error)
                 CALL cp_sm_fm_multiply(t0_state%overlap,t1_state%wf(ispin)%matrix,fm_tmp, k)
                 CALL cp_fm_gemm('T','N',k,k,n,1.0_dp,t0_state%wf(ispin)%matrix,fm_tmp,0.0_dp,csc)
                 CALL cp_fm_gemm('N','N',n,k,k,1.0_dp,t0_state%wf(ispin)%matrix,csc,0.0_dp,fm_tmp)
                 alpha = -1.0_dp * alpha * REAL(nvec - i + 1 , dp ) / REAL(i , dp ) 
                 CALL cp_fm_scale_and_add(1.0_dp,mo_coeff,alpha,fm_tmp)
              ENDDO

              CALL cp_fm_release(csc)
              CALL cp_fm_release(fm_tmp)
              CALL qs_env_reorthogonalize_vectors(qs_env,&
                   v_matrix=mo_coeff,ispin=ispin,&
                   n_col=k, error=error)
              CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
                   density_matrix=rho%rho_ao(ispin)%matrix,error=error)
           END DO
        END IF
        CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw=my_gapw, error=error)
        IF(my_gapw) THEN
            CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
        ENDIF
        IF(my_gapw_xc) THEN
            CALL qs_rho_update_rho(rho_xc,qs_env=qs_env, gapw=my_gapw_xc, error=error)
            CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
        END IF

        CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)
        IF(my_kg_gpw) THEN
           ! ** Updates the molecular rho for each molecule independently ***
           CALL get_qs_env(qs_env=qs_env, kg_sub_pw_env=kg_sub_pw_env)
           CALL kg_rho_update_rho_mol( kg_sub_pw_env, qs_env, kg_fm_mol_set, error=error)
        END IF
    CASE (wfi_aspc)
      ! figure out the actual number of vectors to use in the extrapolation:
      nvec = MIN(wf_history%memory_depth,wf_history%snapshot_count)
      CPPrecondition(nvec.GT.0,cp_failure_level,routineP,error,failure)
      ! Stable relaxation parameter of the corrector (i.e. Hartree damping)
      qs_env%scf_control%p_mix = REAL(nvec,KIND=dp)/REAL(2*nvec - 1,KIND=dp)
      IF (nvec >= wf_history%memory_depth) THEN
        qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_aspc
      END IF
      my_orthogonal_wf = .TRUE.
      CPPrecondition(.NOT.my_kg_gpw,cp_failure_level,routineP,error,failure)
      DO ispin=1,SIZE(mos)
        NULLIFY (mo_coeff,matrix_struct,matrix_struct_new,csc,fm_tmp)
        CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
        CALL cp_fm_get_info(mo_coeff,&
                            nrow_global=n,&
                            ncol_global=k,&
                            matrix_struct=matrix_struct)
        CALL cp_fm_create(fm_tmp,matrix_struct)
        CALL cp_fm_struct_create(matrix_struct_new,&
                                 template_fmstruct=matrix_struct,&
                                 nrow_global=k,&
                                 ncol_global=k)
        CALL cp_fm_create(csc,matrix_struct_new)
        CALL cp_fm_struct_release(matrix_struct_new)
        ! first the most recent
        t1_state => wfi_get_snapshot(wf_history,&
                                     index=1,&
                                     error=error)
        CALL cp_fm_to_fm(t1_state%wf(ispin)%matrix,mo_coeff)
        alpha = REAL(4*nvec - 2,KIND=dp)/REAL(nvec + 1,KIND=dp)
        IF (ionode) PRINT*," B( 1 ) =",alpha
        CALL cp_fm_scale(alpha,mo_coeff)

        DO i=2,nvec
          t0_state => wfi_get_snapshot(wf_history,index=i,error=error)
          CALL cp_sm_fm_multiply(t0_state%overlap,t1_state%wf(ispin)%matrix,fm_tmp,k)
          CALL cp_fm_gemm('T','N',k,k,n,1.0_dp,t0_state%wf(ispin)%matrix,fm_tmp,0.0_dp,csc)
          CALL cp_fm_gemm('N','N',n,k,k,1.0_dp,t0_state%wf(ispin)%matrix,csc,0.0_dp,fm_tmp)
          alpha = (-1.0_dp)**(i + 1)*REAL(i,KIND=dp)*&
                  binomial(2*nvec,nvec - i)/binomial(2*nvec - 2,nvec -1)
          IF (ionode) PRINT*," B(",i,") =",alpha
          CALL cp_fm_scale_and_add(1.0_dp,mo_coeff,alpha,fm_tmp)
        END DO
        CALL cp_fm_release(csc)
        CALL cp_fm_release(fm_tmp)
        CALL qs_env_reorthogonalize_vectors(qs_env,&
                                            v_matrix=mo_coeff,&
                                            ispin=ispin,&
                                            n_col=k,&
                                            error=error)
        CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
                                      density_matrix=rho%rho_ao(ispin)%matrix,&
                                      error=error)
      END DO
      CALL qs_rho_update_rho(rho,qs_env=qs_env,gapw=my_gapw,error=error)
      IF (my_gapw) CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
      IF (my_gapw_xc) THEN
        CALL qs_rho_update_rho(rho_xc,qs_env=qs_env, gapw=my_gapw_xc, error=error)
        CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
      END IF
      CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,error=error)
      !MK
      ! IF (my_kg_gpw) THEN
      !   ! Updates the molecular rho for each molecule independently
      !   CALL get_qs_env(qs_env=qs_env, kg_sub_pw_env=kg_sub_pw_env)
      !   CALL kg_rho_update_rho_mol( kg_sub_pw_env, qs_env, kg_fm_mol_set, error=error)
      ! END IF
      !MK
    CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"Unknown interpolation method: "//&
             TRIM(ADJUSTL(cp_to_string(wf_history%interpolation_method_nr)))//&
             " in "//&
             CPSourceFileRef,&
             error, failure)
     END SELECT
     IF (PRESENT(orthogonal_wf)) orthogonal_wf=my_orthogonal_wf
  END IF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE wfi_extrapolate
!***************************************************************************

!!****f* qs_wf_history_methods/wfi_update [1.0] *
!!
!!   NAME
!!     wfi_update
!!
!!   SYNOPSIS
!!     Subroutine wfi_update(wf_history, qs_env, dt, error)
!!       Type(qs_wf_history_type), Pointer:: wf_history
!!       Type(qs_environment_type), Pointer:: qs_env
!!       Real(KIND = dp), Intent (IN):: dt
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine wfi_update
!!
!!   FUNCTION
!!     updates the snapshot buffer, taking a new snapshot
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - wf_history: the history buffer to update
!!     - qs_env: the qs_env we get the info from
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE wfi_update(wf_history, qs_env, dt, kg_fm_mol_set, error)
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(in)                :: dt
    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER, OPTIONAL                      :: kg_fm_mol_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_update', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,error,failure)
  CPPrecondition(wf_history%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     wf_history%snapshot_count=wf_history%snapshot_count+1
     IF (wf_history%memory_depth>0) THEN
        wf_history%last_state_index=MODULO(wf_history%snapshot_count,&
             wf_history%memory_depth)+1
        CALL wfs_update(snapshot=wf_history%past_states &
             (wf_history%last_state_index)%snapshot,wf_history=wf_history,&
             qs_env=qs_env,dt=dt,kg_fm_mol_set=kg_fm_mol_set,error=error)
     END IF
  END IF
END SUBROUTINE wfi_update
!***************************************************************************

!!****f* qs_wf_history_methods/wfi_duplicate_history [1.0] *
!!
!!   NAME
!!     wfi_duplicate_history
!!
!!
!!   FUNCTION
!!     duplicates the whole qs_wf_history_type structure
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - input_wf_history: the history buffer to be duplicated
!!     - output_wf_history: the out coming duplicate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [TdK]
!!
!!*** **********************************************************************
  SUBROUTINE wfi_duplicate_history(input_wf_history, output_wf_history, qs_env, &
                                   error)
    TYPE(qs_wf_history_type), POINTER        :: input_wf_history, &
                                                output_wf_history    
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_duplicate_history', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: handle, &
                                                i, &
                                                stat

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)
    failure = .FALSE.

    CPPrecondition(ASSOCIATED(input_wf_history), cp_failure_level, routineP, error, failure)
    CPPrecondition(input_wf_history%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      IF (ASSOCIATED(output_wf_history)) THEN
        CALL wfi_clear_history(wf_history=input_wf_history, flush_buffer=.FALSE., &
                               error=error)
      ELSE
        ALLOCATE(output_wf_history, stat=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

        ALLOCATE(output_wf_history%past_states(input_wf_history%memory_depth), stat=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
        IF (.NOT. failure) THEN
          DO i = 1,SIZE(input_wf_history%past_states)
            NULLIFY(output_wf_history%past_states(i)%snapshot)
          END DO
        END IF
      END IF

      output_wf_history%memory_depth = input_wf_history%memory_depth
      output_wf_history%last_state_index = input_wf_history%last_state_index
      output_wf_history%interpolation_method_nr = input_wf_history%interpolation_method_nr
      output_wf_history%snapshot_count = input_wf_history%snapshot_count

      output_wf_history%store_wf = input_wf_history%store_wf
      output_wf_history%store_rho_r = input_wf_history%store_rho_r
      output_wf_history%store_rho_g = input_wf_history%store_rho_g
      output_wf_history%store_rho_ao = input_wf_history%store_rho_ao
      output_wf_history%store_overlap = input_wf_history%store_overlap
      output_wf_history%store_frozen_density = input_wf_history%store_frozen_density

      last_wfi_id=last_wfi_id+1
      output_wf_history%id_nr = last_wfi_id 
      output_wf_history%ref_count = 1

      DO i = 1,SIZE(input_wf_history%past_states)
        IF (ASSOCIATED(input_wf_history%past_states(i)%snapshot)) THEN
          CALL wfs_duplicate_snapshot(input_snapshot=input_wf_history%past_states(i)%snapshot, &
                                      output_snapshot=output_wf_history%past_states(i)%snapshot, &
                                      qs_env=qs_env, error=error)
        END IF
      END DO
    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE wfi_duplicate_history
!***************************************************************************

!!****f* qs_wf_history_methods/wfi_copy_history [1.0] *
!!
!!   NAME
!!     wfi_copy_history
!!
!!
!!   FUNCTION
!!     copies the whole qs_wf_history_type structure
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - input_wf_history: the history buffer to be copied
!!     - output_wf_history: the out coming copy
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     10.2005 created [TdK]
!!
!!*** **********************************************************************
  SUBROUTINE wfi_copy_history(input_wf_history, output_wf_history, qs_env, &
                              error)
    TYPE(qs_wf_history_type), POINTER        :: input_wf_history, &
                                                output_wf_history    
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfi_copy_history', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: handle, &
                                                i, &
                                                stat

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)
    failure = .FALSE.

    CPPrecondition(ASSOCIATED(input_wf_history), cp_failure_level, routineP, error, failure)
    CPPrecondition(input_wf_history%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      output_wf_history%memory_depth = input_wf_history%memory_depth
      output_wf_history%last_state_index = input_wf_history%last_state_index
      output_wf_history%interpolation_method_nr = input_wf_history%interpolation_method_nr
      output_wf_history%snapshot_count = input_wf_history%snapshot_count

      output_wf_history%store_wf = input_wf_history%store_wf
      output_wf_history%store_rho_r = input_wf_history%store_rho_r
      output_wf_history%store_rho_g = input_wf_history%store_rho_g
      output_wf_history%store_rho_ao = input_wf_history%store_rho_ao
      output_wf_history%store_overlap = input_wf_history%store_overlap
      output_wf_history%store_frozen_density = input_wf_history%store_frozen_density

      last_wfi_id=last_wfi_id+1
      output_wf_history%id_nr = last_wfi_id 
      output_wf_history%ref_count = 1

      DO i = 1,SIZE(input_wf_history%past_states)
        IF (ASSOCIATED(input_wf_history%past_states(i)%snapshot)) THEN
          CALL wfs_copy_snapshot(input_snapshot=input_wf_history%past_states(i)%snapshot, &
                                 output_snapshot=output_wf_history%past_states(i)%snapshot, &
                                 qs_env=qs_env, error=error)
        END IF
      END DO
    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE wfi_copy_history
!***************************************************************************

END MODULE qs_wf_history_methods
