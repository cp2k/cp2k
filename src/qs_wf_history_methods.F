!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_wf_history_methods [1.0] *
!!
!!   NAME
!!     qs_wf_history_methods
!!
!!   FUNCTION
!!     methods to interpolate (or actually noramlly extrapolate) the
!!     new guess for 
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_wf_history_methods
  USE coefficient_types,               ONLY: coeff_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_fm_basic_linalg,              ONLY: cp_fm_add
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pools_create_matrix_vect,&
                                             fm_pools_give_back_matrix_vect
  USE cp_fm_types,                     ONLY: cp_full_matrix_p_type,&
                                             cp_full_matrix_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dft_types,                       ONLY: dft_control_type
  USE kinds,                           ONLY: wp=>dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             NOSPACE,&
                                             REALDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_copy,&
                                             pw_p_type,&
                                             pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             qs_wf_snapshot_p_type,&
                                             qs_wf_snapshot_type,&
                                             wfi_get_snapshot
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_wf_history_methods'
  INTEGER, SAVE, PRIVATE :: last_wfs_id=0, last_wfi_id=0

!!***
!****************************************************************************
CONTAINS

!!****f* qs_wf_history_methods/wfs_create [1.0] *
!!
!!   NAME
!!     wfs_create
!!
!!   FUNCTION
!!     allocates and initialize a wavefunction snapshot
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
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
SUBROUTINE wfs_create(snapshot, error)
  TYPE(qs_wf_snapshot_type), POINTER :: snapshot
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='wfs_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: i, stat

  failure=.FALSE.
  
  ALLOCATE(snapshot, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     last_wfs_id=last_wfs_id+1
     snapshot%id_nr=last_wfs_id
     NULLIFY(snapshot%wf, snapshot%rho_r, snapshot%rho_g)
     snapshot%time=0.0_wp
     snapshot%ref_count=1
  END IF
END SUBROUTINE wfs_create
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
!!     - the time of the snapshot
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
SUBROUTINE wfs_update(snapshot,wf_history,qs_env,time,error)
  TYPE(qs_wf_snapshot_type), POINTER :: snapshot
  TYPE(qs_wf_history_type), POINTER :: wf_history
  TYPE(qs_environment_type), INTENT(inout),TARGET :: qs_env
  REAL(kind=wp), INTENT(in) :: time
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='wfs_update',&
        routineP=moduleN//':'//routineN
  TYPE(pw_env_type), POINTER :: pw_env
  TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
  TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER :: ao_mo_pools
  TYPE(dft_control_type), POINTER :: dft_control
  TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mos
  TYPE(cp_full_matrix_type), POINTER :: eigenvectors
  TYPE(qs_rho_type), POINTER :: rho
  TYPE(coeff_type), DIMENSION(:), POINTER :: rho_r, rho_g
  INTEGER :: ispin, nspins, stat

  failure=.FALSE.
  NULLIFY(pw_env, auxbas_pw_pool, ao_mo_pools, dft_control, mos, eigenvectors,&
       rho, rho_r,rho_g)
  CALL get_qs_env(qs_env, pw_env=pw_env, ao_mo_fm_pools=ao_mo_pools,&
       dft_control=dft_control,rho=rho,error=error)
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
     snapshot%time=time

     IF (wf_history%store_wf) THEN
        CALL get_qs_env(qs_env,c=mos,error=error)
        IF (.not.ASSOCIATED(snapshot%wf)) THEN
           CALL fm_pools_create_matrix_vect(ao_mo_pools,snapshot%wf,&
                name="ws_snap"//TRIM(ADJUSTL(cp_to_string(snapshot%id_nr)))//&
                "ws",error=error)
           CPPostcondition(nspins==SIZE(snapshot%wf),cp_failure_level,routineP,error,failure)
        END IF
        DO ispin=1,nspins
           CALL get_mo_set(mos(ispin)%mo_set,eigenvectors=eigenvectors)
           CALL cp_fm_add(alpha=0.0_wp,&
                matrix_a=snapshot%wf(ispin)%matrix,&
                matrix_b=eigenvectors,error=error)
        END DO
     ELSE IF (ASSOCIATED(snapshot%wf)) THEN
        CALL fm_pools_give_back_matrix_vect(ao_mo_pools,snapshot%wf,&
             error=error)
     END IF

     IF (wf_history%store_rho_r) THEN
        CALL qs_rho_get(rho, rho_r=rho_r)
        CPPrecondition(ASSOCIATED(rho_r),cp_failure_level,routineP,error,failure)
         IF (.not.ASSOCIATED(snapshot%rho_r)) THEN
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
         IF (.not.ASSOCIATED(snapshot%rho_g)) THEN
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
      
  END IF
END SUBROUTINE wfs_update
!***************************************************************************

!!****f* qs_wf_history_methods/wfi_create [1.0] *
!!
!!   NAME
!!     wfi_create
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - interpolation_method: the method used for the extrapolation of the
!!       intial density for the next md step
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
SUBROUTINE wfi_create(wf_history, interpolation_method, error)
  TYPE(qs_wf_history_type), POINTER :: wf_history
  CHARACTER(len=60), INTENT(in) :: interpolation_method
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='wfi_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat, i

  failure=.FALSE.
  
  ALLOCATE(wf_history, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     last_wfi_id=last_wfi_id+1
     wf_history%id_nr=last_wfi_id
     wf_history%ref_count=1
     wf_history%memory_depth=0
     wf_history%last_state_index=1
     wf_history%store_wf=.FALSE.
     wf_history%store_rho_r=.FALSE.
     wf_history%store_rho_g=.FALSE.
     NULLIFY(wf_history%past_states)

     wf_history%interpolation_method=interpolation_method
     CALL uppercase(wf_history%interpolation_method)
     SELECT CASE(wf_history%interpolation_method)
     CASE("NONE")
        wf_history%memory_depth=0
     CASE("LINEAR_WF")
        wf_history%memory_depth=1
        wf_history%store_wf=.TRUE.
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"Unknown interpolation method: "//&
             TRIM(interpolation_method)//" in "//&
             CPSourceFileRef,&
             error, failure)
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

!!****f* qs_wf_history_methods/wfi_extrapolate [1.0] *
!!
!!   NAME
!!     wfi_extrapolate
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
!!     - time: the time at which to extrapolate
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
SUBROUTINE wfi_extrapolate(wf_history, qs_env, time, error)
  TYPE(qs_wf_history_type), POINTER :: wf_history
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  REAL(kind=wp), INTENT(in) :: time
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='wfi_extrapolate',&
        routineP=moduleN//':'//routineN
  TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: tmp_matrix
  TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mos
  TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER :: ao_mo_fm_pools
  TYPE(qs_wf_snapshot_type), POINTER :: old_state
  TYPE(cp_full_matrix_type), POINTER :: eigenvectors
  REAL(kind=wp) :: t1,t2,t3
  INTEGER :: ispin

  NULLIFY(tmp_matrix, mos, ao_mo_fm_pools, old_state, eigenvectors)
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(wf_history),cp_failure_level,routineP,error,failure)
  CPPrecondition(wf_history%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL get_qs_env(qs_env,ao_mo_fm_pools=ao_mo_fm_pools, c=mos,&
          error=error)

     SELECT CASE (wf_history%interpolation_method)
     CASE("NONE")
        ! does nothing
     CASE("LINEAR_WF")
        ! with LSD one could spare one temporary (but ugly code)
        old_state => wfi_get_snapshot(wf_history, index=1, error=error)
        CPPrecondition(ASSOCIATED(old_state),cp_warning_level,routineP,error,failure)
        IF (.not.failure) THEN
           CPPrecondition(ASSOCIATED(old_state%wf),cp_warning_level,routineP,error,failure)
        END IF
        IF (.not.failure) THEN
           t1=0.0_wp
           t2=wf_history%actual_time
           t3=wf_history%actual_time+time
           CALL fm_pools_create_matrix_vect(ao_mo_fm_pools, &
                tmp_matrix, error=error)
           DO ispin=1,SIZE(tmp_matrix)
              CALL cp_fm_add(alpha=0.0_wp,matrix_a=tmp_matrix(ispin)%matrix,&
                   matrix_b=old_state%wf(ispin)%matrix,&
                   error=error)
           END DO
        END IF
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"Unknown interpolation method: "//&
             TRIM(wf_history%interpolation_method)//" in "//&
             CPSourceFileRef,&
             error, failure)
     END SELECT
     
     wf_history%last_state_index=MOD(wf_history%last_state_index+1,&
          wf_history%memory_depth)+1
     CALL wfs_update(snapshot=wf_history%past_states &
          (wf_history%last_state_index)%snapshot,wf_history=wf_history,&
          qs_env=qs_env,time=wf_history%actual_time,error=error)
        
     SELECT CASE (wf_history%interpolation_method)
     CASE("NONE")
        ! does nothing
     CASE("LINEAR_WF")
        IF (.not.failure) THEN
           DO ispin=1,SIZE(mos)
              CALL get_mo_set(mos(ispin)%mo_set,eigenvectors=eigenvectors)
              CALL cp_fm_add(alpha=0.0_wp,&
                   matrix_a=tmp_matrix(ispin)%matrix,&
                   matrix_b=eigenvectors,error=error)
              CALL cp_fm_add(alpha=(t3-t1)/(t2-t1),matrix_a=eigenvectors,&
                   beta=(t2-t3)/(t2-t1), matrix_b=tmp_matrix(ispin)%matrix,&
                   error=error)
           END DO
        END IF

        IF (ASSOCIATED(tmp_matrix)) THEN
           CALL fm_pools_give_back_matrix_vect(ao_mo_fm_pools,tmp_matrix,&
                error=error)
        END IF

     CASE default
        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
     END SELECT

     wf_history%actual_time=time
  END IF
END SUBROUTINE wfi_extrapolate
!***************************************************************************

END MODULE qs_wf_history_methods
