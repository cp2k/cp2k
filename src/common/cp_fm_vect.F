!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routine to handle vectors of full matrixes
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_fm_vect
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_retain,&
                                             cp_fm_set_all,&
                                             cp_fm_write
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_vect'

  PUBLIC :: cp_fm_vect_dealloc, cp_fm_vect_copy,&
            cp_fm_vect_set_all, cp_fm_vect_write
!***
CONTAINS

! *****************************************************************************
!> \brief deallocate an array of pointers to blacs matrixes
!> \param matrixes the array of matrixes to deallocate
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_fm_vect_dealloc(matrixes, error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vect_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(matrixes)) THEN
       DO i=1,SIZE(matrixes)
          CALL cp_fm_release(matrixes(i)%matrix,error=error)
       END DO
       DEALLOCATE(matrixes,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
  END SUBROUTINE cp_fm_vect_dealloc

! *****************************************************************************
!> \brief Does a shallow copy of an array of full matrices (i.e. just retains
!>      the matrices)
!> \param matrixes the matrixes to copy
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_fm_vect_copy(matrixes, copy, error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(in)                             :: matrixes
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: copy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vect_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  ALLOCATE(copy(SIZE(matrixes)),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(matrixes)
        copy(i)%matrix => matrixes(i)%matrix
        CALL cp_fm_retain(matrixes(i)%matrix,error=error)
     END DO
  END IF
END SUBROUTINE cp_fm_vect_copy

! *****************************************************************************
!> \brief writes an array of full matrixes to the given output unit
!> \param matrixes the array of matrixes to output
!> \param unit_nr the unit where to output
!> \param long_description if a long description should be written
!>        (defaults to false)
!> \param local if the unit is a local unit or a global unit
!>        (defaults to false, i.e. global)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      routine to help debugging
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
   SUBROUTINE cp_fm_vect_write(matrixes, unit_nr, long_description,&
        local, error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    INTEGER, INTENT(in)                      :: unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: long_description, local
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vect_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, my_local, &
                                                my_long_description, should_w
    TYPE(cp_logger_type), POINTER            :: logger

   failure=.FALSE.; my_local=.FALSE.; my_long_description=.FALSE.
   logger => cp_error_get_logger(error)

   IF (PRESENT(local)) my_local=local
   IF (PRESENT(long_description)) my_long_description=long_description

   should_w=my_local .OR. logger%para_env%mepos==logger%para_env%source

   IF (ASSOCIATED(matrixes)) THEN
      IF (should_w) WRITE (unit=unit_nr,fmt="(a)") '('
      DO i=1,SIZE(matrixes)
         CALL cp_fm_write(matrix=matrixes(i)%matrix, unit_nr=unit_nr,&
              long_description=my_long_description,&
              local=my_local,error=error)
         IF (i/=SIZE(matrixes).and.should_w) WRITE(unit=unit_nr,fmt="(a)")','
      END DO
      IF (should_w) WRITE (unit=unit_nr,fmt="(a)") ')'
   ELSE
      IF (should_w) THEN
         WRITE (unit=unit_nr,fmt="(a)") ' <cp_fm_p array>:*null*'
      END IF
   END IF
 END SUBROUTINE cp_fm_vect_write

! *****************************************************************************
!> \brief sets all the elements of the matrixes to the given value
!> \param matrixes the matrixes to change
!> \param value the value that the element should receive (defaults to 0)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE cp_fm_vect_set_all(matrixes,value,error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: value
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_vect_set_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_value

  failure=.FALSE.
  my_value=0.0_dp

  IF (PRESENT(value)) my_value=value

  CPPrecondition(ASSOCIATED(matrixes),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(matrixes)
        CPPrecondition(ASSOCIATED(matrixes(i)%matrix),cp_failure_level,routineP,error,failure)
        CALL cp_fm_set_all(matrixes(i)%matrix,my_value,error=error)
     END DO
  END IF
END SUBROUTINE cp_fm_vect_set_all

END MODULE cp_fm_vect
