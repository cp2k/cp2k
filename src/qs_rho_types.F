!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief superstucture that hold various representations of the density and
!>      keeps track of which ones are valid
!> \par History
!>      08.2002 created [fawzi]
!>      08.2014 kpoints [JGH]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_rho_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_release
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho_types'

  PUBLIC :: qs_rho_p_type, qs_rho_type
  PUBLIC :: qs_rho_retain, qs_rho_release,&
            qs_rho_get

! *****************************************************************************
!> \brief keeps the density in various representations, keeping track of
!>      which ones are valid.
!> \param most attributes are array with either lda or lsd_alpha,lsd_beta.
!> \param rho_ao the filtered rho in the localized atom basis (to have rho(r)
!>        the filtred matrix is enough, but rho(r,r') is lost).
!> \param rho_ao_kp the filtered rho in the localized atom basis (to have rho(r)
!>        the filtered matrix is enough, but rho(r,r') is lost).
!>        for kpoints, in real space index form
!> \param rho_r grids with rho in the real space
!> \param tau_r grids with the kinetic energy density in real space
!> \param rho_g grids with rho in the g space
!> \param tau_g grids with the kinetic energy density in g space
!> \param rho_g_valid , rho_r_valid, tau_r_valid, tau_g_valid: if the
!>        corresponding component is valid
!> \param ref_count the reference count, when it becomes 0 the type
!>        is deallocated.
!> \param rebuild_each how often a rebuild should be done by default
!> \param tot_rho_r the total charge in r space (valid only if rho_r is)
!> \note
!>      If pw_p_type would implement retain/release it would be nice to
!>      store also the core charge and the qs_charges in this structure...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_rho_type
     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: rho_ao
     TYPE(cp_dbcsr_p_type), DIMENSION(:,:), POINTER :: rho_ao_kp
     TYPE(pw_p_type), DIMENSION(:), POINTER :: rho_g, rho_r
     TYPE(pw_p_type), DIMENSION(:), POINTER :: drho_g, drho_r
     TYPE(pw_p_type), DIMENSION(:), POINTER :: tau_g, tau_r
     LOGICAL ::  rho_g_valid, &
          rho_r_valid, &
          drho_r_valid, &
          drho_g_valid, &
          tau_r_valid, &
          tau_g_valid, &
          soft_valid
     INTEGER :: ref_count, &
          id_nr, &
          rebuild_each
     REAL(KIND = dp), DIMENSION(:), POINTER :: tot_rho_r, tot_rho_g
  END TYPE qs_rho_type

! *****************************************************************************
  TYPE qs_rho_p_type
    TYPE(qs_rho_type), POINTER         :: rho
  END TYPE qs_rho_p_type

CONTAINS

! *****************************************************************************
!> \brief retains a rho_struct by increasing the reference count by one
!>      (to be called when you want to keep a shared copy)
!> \param rho_struct the structure to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_rho_retain(rho_struct,error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
       rho_struct%ref_count=rho_struct%ref_count+1
    END IF
  END SUBROUTINE qs_rho_retain

! *****************************************************************************
!> \brief releases a rho_struct by decreasing the reference count by one
!>      and deallocating if it reaches 0 (to be called when you don't want
!>      anymore a shared copy)
!> \param rho_struct the structure to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_rho_release(rho_struct,error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(rho_struct)) THEN
       CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
       rho_struct%ref_count=rho_struct%ref_count-1
       IF (rho_struct%ref_count<1) THEN
          IF (ASSOCIATED(rho_struct%rho_ao))&
             CALL cp_dbcsr_deallocate_matrix_set(rho_struct%rho_ao,error=error)
          IF (ASSOCIATED(rho_struct%rho_ao_kp))&
             CALL cp_dbcsr_deallocate_matrix_set(rho_struct%rho_ao_kp,error=error)
          IF (ASSOCIATED(rho_struct%rho_r)) THEN
             DO i=1,SIZE(rho_struct%rho_r)
                CALL pw_release(rho_struct%rho_r(i)%pw,error=error)
             END DO
             DEALLOCATE(rho_struct%rho_r,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(rho_struct%drho_r)) THEN
             DO i=1,SIZE(rho_struct%drho_r)
                CALL pw_release(rho_struct%drho_r(i)%pw,error=error)
             END DO
             DEALLOCATE(rho_struct%drho_r,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(rho_struct%drho_g)) THEN
             DO i=1,SIZE(rho_struct%drho_g)
                CALL pw_release(rho_struct%drho_g(i)%pw,error=error)
             END DO
             DEALLOCATE(rho_struct%drho_g,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(rho_struct%tau_r)) THEN
             DO i=1,SIZE(rho_struct%tau_r)
                CALL pw_release(rho_struct%tau_r(i)%pw,error=error)
             END DO
             DEALLOCATE(rho_struct%tau_r,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(rho_struct%rho_g)) THEN
             DO i=1,SIZE(rho_struct%rho_g)
                CALL pw_release(rho_struct%rho_g(i)%pw,error=error)
             END DO
             DEALLOCATE(rho_struct%rho_g,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(rho_struct%tau_g)) THEN
             DO i=1,SIZE(rho_struct%tau_g)
                CALL pw_release(rho_struct%tau_g(i)%pw,error=error)
             END DO
             DEALLOCATE(rho_struct%tau_g,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(rho_struct%tot_rho_r)) THEN
             DEALLOCATE(rho_struct%tot_rho_r,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          ENDIF
          IF (ASSOCIATED(rho_struct%tot_rho_g)) THEN
             DEALLOCATE(rho_struct%tot_rho_g,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          ENDIF
          DEALLOCATE(rho_struct,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(rho_struct)
  END SUBROUTINE qs_rho_release

! *****************************************************************************
!> \brief this function should be called wehn something changes to update the
!>      internal status
!> \param rho_struct the density structure that changed
!> \param rho_r_valid , rho_g_valid, tau_r_valid, tau_g_valid: actual valid
!>        representations
!> \param drho_r_valid ...
!> \param tau_r_valid ...
!> \param rho_g_valid ...
!> \param drho_g_valid ...
!> \param tau_g_valid ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_rho_did_change(rho_struct, rho_r_valid, drho_r_valid, tau_r_valid,&
       rho_g_valid, drho_g_valid, tau_g_valid, error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    LOGICAL, INTENT(in), OPTIONAL            :: rho_r_valid, drho_r_valid, &
                                                tau_r_valid, rho_g_valid, &
                                                drho_g_valid, tau_g_valid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_did_change', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(rho_r_valid)) rho_struct%rho_r_valid=rho_r_valid
       IF (PRESENT(rho_g_valid)) rho_struct%rho_g_valid=rho_g_valid
       IF (PRESENT(drho_r_valid)) rho_struct%drho_r_valid=drho_r_valid
       IF (PRESENT(drho_g_valid)) rho_struct%drho_g_valid=drho_g_valid
       IF (PRESENT(tau_r_valid)) rho_struct%tau_r_valid=tau_r_valid
       IF (PRESENT(tau_g_valid)) rho_struct%tau_g_valid=tau_g_valid
    END IF
  END SUBROUTINE qs_rho_did_change

! *****************************************************************************
!> \brief returns info about the density described by this object.
!>      If some representation is not available an error is issued
!> \param rho_struct ...
!> \param rho_ao ...
!> \param rho_ao_kp ...
!> \param rho_r ...
!> \param drho_r ...
!> \param rho_g ...
!> \param drho_g ...
!> \param tau_r ...
!> \param tau_g ...
!> \param rho_r_valid ...
!> \param drho_r_valid ...
!> \param rho_g_valid ...
!> \param drho_g_valid ...
!> \param tau_r_valid ...
!> \param tau_g_valid ...
!> \param rebuild_each ...
!> \param tot_rho_r ...
!> \param soft_valid ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>      For the other inputs the qs_rho_type attributes
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_rho_get(rho_struct, rho_ao, rho_ao_kp, rho_r, drho_r, &
       rho_g, drho_g, tau_r, tau_g, rho_r_valid, drho_r_valid, rho_g_valid, &
       drho_g_valid, tau_r_valid, tau_g_valid, rebuild_each, tot_rho_r, soft_valid, &
       error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_ao
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: rho_ao_kp
    TYPE(pw_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_r, drho_r, rho_g, drho_g, &
                                                tau_r, tau_g
    LOGICAL, INTENT(out), OPTIONAL           :: rho_r_valid, drho_r_valid, &
                                                rho_g_valid, drho_g_valid, &
                                                tau_r_valid, tau_g_valid
    INTEGER, INTENT(out), OPTIONAL           :: rebuild_each
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: tot_rho_r
    LOGICAL, INTENT(out), OPTIONAL           :: soft_valid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(rho_ao)) THEN
             rho_ao => rho_struct%rho_ao
       END IF
       IF (PRESENT(rho_ao_kp)) THEN
             rho_ao_kp => rho_struct%rho_ao_kp
       END IF
       IF (PRESENT(rho_r)) THEN
          IF (PRESENT(rho_r_valid).or.rho_struct%rho_r_valid) THEN
             rho_r => rho_struct%rho_r
          ELSE
             CALL cp_unimplemented_error(routineP,"to do",error=error)
             NULLIFY(rho_r)
          END IF
       END IF
       IF (PRESENT(drho_r)) THEN
          IF (PRESENT(drho_r_valid).or.rho_struct%drho_r_valid) THEN
             drho_r => rho_struct%drho_r
          ELSE
             CALL cp_unimplemented_error(routineP,"to do",error=error)
             NULLIFY(drho_r)
          END IF
       END IF
       IF (PRESENT(rho_g)) THEN
          IF (PRESENT(rho_g_valid).or.rho_struct%rho_g_valid) THEN
             rho_g => rho_struct%rho_g
          ELSE
             CALL cp_unimplemented_error(routineP,"to do",error=error)
             NULLIFY(rho_g)
          END IF
       END IF
       IF (PRESENT(drho_g)) THEN
          IF (PRESENT(drho_g_valid).or.rho_struct%drho_g_valid) THEN
             drho_g => rho_struct%drho_g
          ELSE
             CALL cp_unimplemented_error(routineP,"to do",error=error)
             NULLIFY(drho_g)
          END IF
       END IF
       IF (PRESENT(tau_r)) THEN
          IF (PRESENT(tau_r_valid).or.rho_struct%tau_r_valid) THEN
             tau_r => rho_struct%tau_r
          ELSE
             CALL cp_unimplemented_error(routineP,"to do",error=error)
             NULLIFY(tau_r)
          END IF
       END IF
       IF (PRESENT(tau_g)) THEN
          IF (PRESENT(tau_g_valid).or.rho_struct%tau_g_valid) THEN
             tau_g => rho_struct%tau_g
          ELSE
             CALL cp_unimplemented_error(routineP,"to do",error=error)
             NULLIFY(tau_g)
          END IF
       END IF
       IF (PRESENT(rho_r_valid)) rho_r_valid=rho_struct%rho_r_valid
       IF (PRESENT(rho_g_valid)) rho_g_valid=rho_struct%rho_g_valid
       IF (PRESENT(drho_r_valid)) drho_r_valid=rho_struct%drho_r_valid
       IF (PRESENT(drho_g_valid)) drho_g_valid=rho_struct%drho_g_valid
       IF (PRESENT(tau_r_valid)) tau_r_valid=rho_struct%tau_r_valid
       IF (PRESENT(tau_g_valid)) tau_g_valid=rho_struct%tau_g_valid
       IF (PRESENT(soft_valid))  soft_valid=rho_struct%soft_valid
       IF (PRESENT(rebuild_each)) rebuild_each=rho_struct%rebuild_each
       ! give error if not valid?
       IF (PRESENT(tot_rho_r)) tot_rho_r => rho_struct%tot_rho_r
    END IF
  END SUBROUTINE qs_rho_get

END MODULE qs_rho_types

