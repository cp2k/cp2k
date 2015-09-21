!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lgrid_type, lgrid_release, lgrid_create, lgrid_allocate_grid

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
!> \par History
!>      10.2011 created [IAB]
!> \author Iain Bethune
! *****************************************************************************
SUBROUTINE lgrid_create(lgrid,rs_descs)
    TYPE(lgrid_type), POINTER                :: lgrid
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs

    CHARACTER(len=*), PARAMETER :: routineN = 'lgrid_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ngpts
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(.NOT.ASSOCIATED(lgrid))
  ALLOCATE(lgrid)
  NULLIFY(lgrid%r)
  lgrid%ref_count=1
  ! Find the maximum number of grid points needed
  ngpts = 0
  DO i=1,SIZE(rs_descs)
    ngpts = MAX(ngpts, rs_grid_max_ngpts(rs_descs(i)%rs_desc))
  END DO
  lgrid%ldim = ngpts
END SUBROUTINE lgrid_create

! *****************************************************************************
!> \brief retains the lgrid (see doc/ReferenceCounting.html)
!> \param lgrid the lgrid_type to retain
!> \par History
!>      10.2011 created [IAB]
!> \author Iain Bethune
! *****************************************************************************
SUBROUTINE lgrid_retain(lgrid)
    TYPE(lgrid_type), POINTER                :: lgrid

    CHARACTER(len=*), PARAMETER :: routineN = 'lgrid_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(lgrid))
  CPASSERT(lgrid%ref_count>0)
  lgrid%ref_count=lgrid%ref_count+1
END SUBROUTINE lgrid_retain

! *****************************************************************************
!> \brief releases the given lgrid (see doc/ReferenceCounting.html)
!> \param lgrid the lgrid_type to release
!> \par History
!>      10.2011 created [IAB]
!> \author Iain Bethune
! *****************************************************************************
SUBROUTINE lgrid_release(lgrid)
    TYPE(lgrid_type), POINTER                :: lgrid

    CHARACTER(len=*), PARAMETER :: routineN = 'lgrid_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(lgrid)) THEN
     CPASSERT(lgrid%ref_count>0)
     lgrid%ref_count=lgrid%ref_count-1
     IF (lgrid%ref_count<1) THEN
        IF (ASSOCIATED(lgrid%r)) THEN
           DEALLOCATE (lgrid%r)
        END IF
        DEALLOCATE (lgrid)
        NULLIFY(lgrid)
     END IF
  END IF
END SUBROUTINE

! *****************************************************************************
!> \brief allocates the lgrid for a given number of threads
!> \param lgrid the lgrid_type for which the grid will be allocated
!> \param nthreads how many threads to allocate for
!> \par History
!>      10.2011 created [IAB]
!> \author Iain Bethune
! *****************************************************************************
SUBROUTINE lgrid_allocate_grid(lgrid, nthreads)
    TYPE(lgrid_type), POINTER                :: lgrid
    INTEGER, INTENT(in)                      :: nthreads

    CHARACTER(len=*), PARAMETER :: routineN = 'lgrid_allocate_grid', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(ASSOCIATED(lgrid))
    CPASSERT(.NOT. ASSOCIATED(lgrid%r))
    ALLOCATE(lgrid%r(lgrid%ldim*nthreads))
END SUBROUTINE

END MODULE

