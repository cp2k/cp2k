!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief buffer for the diis of the scf
!> \par History
!>      02.2003 rewamped [fawzi]
!> \author Matthias Krack
! *****************************************************************************
MODULE qs_diis_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_release
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_diis_types'

  PUBLIC :: qs_diis_buffer_type
  PUBLIC :: qs_diis_b_release
  PUBLIC :: qs_diis_buffer_type_sparse, &
            qs_diis_b_release_sparse

! *****************************************************************************
!> \brief keeps a buffer with the previous values of s,p,k
!> \par History
!>      02.2003 rewamped [fawzi]
!> \author Matthias Krack
! *****************************************************************************
  TYPE qs_diis_buffer_type
    INTEGER                                          :: nbuffer,ncall,&
         id_nr,ref_count
    TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER :: error,PARAMETER
    REAL(KIND = dp), DIMENSION(:,:), POINTER                :: b_matrix
  END TYPE qs_diis_buffer_type

! *****************************************************************************
!> \brief build array of pointers to diis buffers
!> \param diis_buffer the diis buffer pointer
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  TYPE qs_diis_buffer_p_type
     TYPE(qs_diis_buffer_type), POINTER :: diis_buffer
  END TYPE qs_diis_buffer_p_type

! *****************************************************************************
!> \brief build array of pointers to diis buffers for sparse matrix case
!> \param diis_buffer the diis buffer pointer
!> \par History
!>      10.2014 Modified from non-sparse case by Fredy W. Aquino
!> \author fwaq
! *****************************************************************************
  TYPE qs_diis_buffer_type_sparse
    INTEGER                                          :: nbuffer,ncall,&
                                                        id_nr,ref_count
    TYPE(cp_dbcsr_p_type), DIMENSION(:,:), POINTER   :: error,PARAMETER
    REAL(KIND = dp), DIMENSION(:,:), POINTER         :: b_matrix
  END TYPE qs_diis_buffer_type_sparse

  TYPE qs_diis_buffer_p_type_sparse
     TYPE(qs_diis_buffer_type_sparse), POINTER       :: diis_buffer
  END TYPE qs_diis_buffer_p_type_sparse

CONTAINS

! *****************************************************************************
!> \brief retains a diis buffer (see doc/ReferenceCounting.html)
!> \param diis_buffer the buffer to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE qs_diis_b_retain(diis_buffer,error)
    TYPE(qs_diis_buffer_type), POINTER       :: diis_buffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_diis_b_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(diis_buffer),cp_failure_level,routineP,error,failure)
  CPPreconditionNoFail(diis_buffer%ref_count>0,cp_failure_level,routineP,error)
  diis_buffer%ref_count=diis_buffer%ref_count+1
END SUBROUTINE qs_diis_b_retain

! *****************************************************************************
!> \brief releases the given diis buffer (see doc/ReferenceCounting.html)
!> \param diis_buffer the buffer to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE qs_diis_b_release(diis_buffer,error)
    TYPE(qs_diis_buffer_type), POINTER       :: diis_buffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_diis_b_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(diis_buffer)) THEN
     CPPreconditionNoFail(diis_buffer%ref_count>0,cp_failure_level,routineP,error)
     diis_buffer%ref_count=diis_buffer%ref_count-1
     IF (diis_buffer%ref_count<1) THEN
        IF (ASSOCIATED(diis_buffer%b_matrix)) THEN
           DEALLOCATE(diis_buffer%b_matrix,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(diis_buffer%error)) THEN
           DO j=1,SIZE(diis_buffer%error,2)
              DO i=1,SIZE(diis_buffer%error,1)
                 CALL cp_fm_release(diis_buffer%error(i,j)%matrix,error=error)
              END DO
           END DO
           DEALLOCATE(diis_buffer%error,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(diis_buffer%parameter)) THEN
           DO j=1,SIZE(diis_buffer%parameter,2)
              DO i=1,SIZE(diis_buffer%parameter,1)
                 CALL cp_fm_release(diis_buffer%parameter(i,j)%matrix,&
                      error=error)
              END DO
           END DO
           DEALLOCATE(diis_buffer%parameter,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        DEALLOCATE(diis_buffer,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE qs_diis_b_release

! *****************************************************************************
!> \brief releases the given diis buffer (see doc/ReferenceCounting.html)
!> \param diis_buffer the buffer to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10-11-14 created [FA] modified from qs_diis_b_release
!> \author Fredy W. Aquino
! *****************************************************************************
SUBROUTINE qs_diis_b_release_sparse(diis_buffer,error)

    TYPE(qs_diis_buffer_type_sparse), &
      POINTER                                :: diis_buffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_diis_b_release_sparse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(diis_buffer)) THEN
        IF (ASSOCIATED(diis_buffer%b_matrix)) THEN
           DEALLOCATE(diis_buffer%b_matrix,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(diis_buffer%error)) THEN
           DO j=1,SIZE(diis_buffer%error,2)
              DO i=1,SIZE(diis_buffer%error,1)
               CALL cp_dbcsr_release(diis_buffer%error(i,j)%matrix, &
                                     error=error)
               DEALLOCATE(diis_buffer%error(i,j)%matrix)
              END DO
           END DO
           DEALLOCATE(diis_buffer%error,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(diis_buffer%parameter)) THEN
           DO j=1,SIZE(diis_buffer%parameter,2)
              DO i=1,SIZE(diis_buffer%parameter,1)
               CALL cp_dbcsr_release(diis_buffer%parameter(i,j)%matrix, &
                                     error=error)
               DEALLOCATE(diis_buffer%parameter(i,j)%matrix)
              END DO
           END DO
           DEALLOCATE(diis_buffer%parameter,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        DEALLOCATE(diis_buffer,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
END SUBROUTINE qs_diis_b_release_sparse

END MODULE qs_diis_types
