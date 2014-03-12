!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief a module to allow simple internal preprocessing in input files.
!> \par History
!>      - standalone proof-of-concept implementation (20.02.2008,AK)
!>      - integration into cp2k (22.02.2008,tlaino)
!>      - variables added (25.02.2008,AK)
!> \author Axel Kohlmeyer [AK] - CMM/UPenn Philadelphia
!> \date 25.02.2008
! *****************************************************************************
MODULE cp_parser_ilist_methods
  USE cp_parser_ilist_types,           ONLY: ilist_type
  USE f77_blas
#include "cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ilist_setup, ilist_update, ilist_reset
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_parser_ilist_methods'

CONTAINS

! ****************************************************************************
!> \brief setup the integer listing type
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date  08.2008
! *****************************************************************************
  SUBROUTINE ilist_setup(ilist, token, error)
    TYPE(ilist_type), POINTER                :: ilist
    CHARACTER(LEN=*)                         :: token
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ilist_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ind
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(ilist),cp_failure_level,routineP,error,failure)
    ind = INDEX(token,"..")
    READ (UNIT=token(:ind-1),FMT=*) ilist%istart
    READ (UNIT=token(ind+2:),FMT=*) ilist%iend
    CALL cp_assert((ilist%istart <= ilist%iend),cp_failure_level,cp_assertion_failed,routineP,&
                   "Invalid list range specified: "//&
                   TRIM(ADJUSTL(cp_to_string(ilist%istart)))//".."//&
                   TRIM(ADJUSTL(cp_to_string(ilist%iend))),&
                   only_ionode=.TRUE.)
    ilist%nel_list = ilist%iend - ilist%istart + 1
    ilist%ipresent = ilist%istart
    ilist%in_use   = .TRUE.

  END SUBROUTINE ilist_setup

! ****************************************************************************
!> \brief updates the integer listing type
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date  08.2008
! *****************************************************************************
  SUBROUTINE ilist_update(ilist, error)
    TYPE(ilist_type), POINTER                :: ilist
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ilist_update', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(ilist),cp_failure_level,routineP,error,failure)
    ilist%ipresent = ilist%ipresent + 1
    IF (ilist%ipresent>ilist%iend) THEN
       CALL ilist_reset(ilist, error)
    END IF
  END SUBROUTINE ilist_update

! ****************************************************************************
!> \brief updates the integer listing type
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date  08.2008
! *****************************************************************************
  SUBROUTINE ilist_reset(ilist, error)
    TYPE(ilist_type), POINTER                :: ilist
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ilist_reset', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(ilist),cp_failure_level,routineP,error,failure)
    IF (ilist%ipresent==ilist%iend) THEN
       ilist%istart   = HUGE(0)
       ilist%iend     = HUGE(0)
       ilist%nel_list = HUGE(0)
       ilist%ipresent = HUGE(0)
       ilist%in_use   = .FALSE.
    END IF
  END SUBROUTINE ilist_reset

END MODULE cp_parser_ilist_methods
