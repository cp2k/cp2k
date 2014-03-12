!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Scratch space for multiple threads writing to rs grids (see
!>        qs_collocate_density.F for an example
!> \par History
!>      IAB 26-Apr-2010 : initial version - moved out of qs_collocate_density.F
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2010 on behalf of the HECToR project
!> \author IAB
! *****************************************************************************

MODULE lgrid_types

  USE kinds,                           ONLY: dp
  USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                             rs_grid_max_ngpts
#include "../cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lgrid_type, lgrid_p_type, lgrid_release, lgrid_create, lgrid_retain, lgrid_allocate_grid

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lgrid_types'

  TYPE lgrid_type
     INTEGER :: ldim, ref_count
     REAL(dp), DIMENSION(:), POINTER :: r
  END TYPE lgrid_type

  TYPE lgrid_p_type
     TYPE(lgrid_type), POINTER :: l
  END TYPE lgrid_p_type

CONTAINS

! *****************************************************************************
!> \brief creates an lgrid, ldim set based on the rs_grid_descriptors.
!>        The grid is not allocated
!> \param lgrid the lgrid that gets created
!> \param rs_descs the rs grid descriptors used to set the lgrid size
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2011 created [IAB]
!> \author Iain Bethune
! *****************************************************************************
SUBROUTINE lgrid_create(lgrid,rs_descs,error)
    TYPE(lgrid_type), POINTER                :: lgrid
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lgrid_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ngpts, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(.NOT.ASSOCIATED(lgrid),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(lgrid, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     NULLIFY(lgrid%r)
     lgrid%ref_count=1
     ! Find the maximum number of grid points needed
     ngpts = 0
     DO i=1,SIZE(rs_descs)
       ngpts = MAX(ngpts, rs_grid_max_ngpts(rs_descs(i)%rs_desc))
     END DO
     lgrid%ldim = ngpts
  END IF
END SUBROUTINE lgrid_create

! *****************************************************************************
!> \brief retains the lgrid (see doc/ReferenceCounting.html)
!> \param lgrid the lgrid_type to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2011 created [IAB]
!> \author Iain Bethune
! *****************************************************************************
SUBROUTINE lgrid_retain(lgrid,error)
    TYPE(lgrid_type), POINTER                :: lgrid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lgrid_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(lgrid),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(lgrid%ref_count>0,cp_failure_level,routineP,error,failure)
     lgrid%ref_count=lgrid%ref_count+1
  END IF
END SUBROUTINE lgrid_retain

! *****************************************************************************
!> \brief releases the given lgrid (see doc/ReferenceCounting.html)
!> \param lgrid the lgrid_type to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2011 created [IAB]
!> \author Iain Bethune
! *****************************************************************************
SUBROUTINE lgrid_release(lgrid, error)
    TYPE(lgrid_type), POINTER                :: lgrid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lgrid_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(lgrid)) THEN
     CPPrecondition(lgrid%ref_count>0,cp_failure_level,routineP,error,failure)
     lgrid%ref_count=lgrid%ref_count-1
     IF (lgrid%ref_count<1) THEN
        IF (ASSOCIATED(lgrid%r)) THEN
           DEALLOCATE (lgrid%r,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        DEALLOCATE (lgrid,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        NULLIFY(lgrid)
     END IF
  END IF
END SUBROUTINE

! *****************************************************************************
!> \brief allocates the lgrid for a given number of threads
!> \param lgrid the lgrid_type for which the grid will be allocated
!> \param nthreads how many threads to allocate for
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2011 created [IAB]
!> \author Iain Bethune
! *****************************************************************************
SUBROUTINE lgrid_allocate_grid(lgrid, nthreads, error)
    TYPE(lgrid_type), POINTER                :: lgrid
    INTEGER, INTENT(in)                      :: nthreads
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lgrid_allocate_grid', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(lgrid),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(.NOT. ASSOCIATED(lgrid%r),cp_failure_level,routineP,error)
       ALLOCATE(lgrid%r(lgrid%ldim*nthreads))
    END IF
END SUBROUTINE

END MODULE

