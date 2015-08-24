!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!>   \brief
!>     Routines to efficently collocate and integrate gaussians on a grid
!>     These use most of Joost's tricks and a couple more...
!>     result is *speed* and genericity
!>   \author Fawzi Mohamed, 2007
!>   \notes original available with BSD style license
! *****************************************************************************
MODULE gauss_colloc
  USE d3_poly
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE lgrid_types,                     ONLY: lgrid_type

  !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
#include "./common/cp_common_uses.f90"

IMPLICIT NONE
PRIVATE

  PUBLIC :: collocGauss,&
            integrateGaussFull

#ifdef FD_DEBUG
#define IF_CHECK(x,y) x
#else
#define IF_CHECK(x,y) y
#endif

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gauss_colloc'

    REAL(dp), PARAMETER :: small=TINY(1.0_dp)

  ! keep prettify happy (does not see the include)
  INTEGER(KIND=int_8), PARAMETER           :: unused_import_of_int_8=1


CONTAINS


! *****************************************************************************
!> \brief collocate a periodically repeated gaussian on a non orthormbic grid
!>
!> this routine has been tested and works well with cells with
!> det(h)/sqrt(tr(dot(h^T,h)))>0.2 (2 angles bigger than 24 deg or one angle
!> bigger than 11 deg).
!> Because of its numerics it might fail badly (infinity or NaN) with
!> with more deformed cells. Avoiding this would be bossible only using
!> IEEE numerics controls, which would also make everything slower and
!> less supported.
!> Still the actual numeric has been carefully tuned, and in normal cases
!> and most anormal it should work.
!> With det(h)/sqrt(tr(dot(h^T,h)))>0.2 I could not find any failure.
!>
!> \param h cell matrix
!> \param h_inv inverse of the cell matrix
!> \param grid the grid
!> \param poly polynomial (d3_poly format)
!> \param alphai exponential coeff
!> \param posi position of the gaussian
!> \param max_r2 maximum radius of collocation squared
!> \param periodic array of 0 or 1 that says which dimensions have pbc (1=pbc)
!> \param gdim dimension of the grid (grid might be a subset)
!> \param local_bounds local bounds of the grid piece that is kept locally
!>   (i.e. of grid) the global grid is assumed to atart at 0,0,0
!> \param local_shift start indexes of the local slice (i.e. of grid)
!> \param poly_shift position of posi in the polynomial reference system.
!>  Set it to posi to use the global reference system.
!> \param scale a global scale factor
!> \param lgrid ...
!> \param error type to control the error handling
! *****************************************************************************
SUBROUTINE collocGauss(h,h_inv,grid,poly,alphai,posi,max_r2,&
        periodic,gdim,local_bounds,local_shift,poly_shift,scale,lgrid,error)
    REAL(dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: h, h_inv
    REAL(dp), DIMENSION(0:, 0:, 0:), &
      INTENT(inout)                          :: grid
    REAL(dp), DIMENSION(:), INTENT(in)       :: poly
    REAL(dp), INTENT(in)                     :: alphai
    REAL(dp), DIMENSION(0:2), INTENT(in)     :: posi
    REAL(dp), INTENT(in)                     :: max_r2
    INTEGER, DIMENSION(0:2), INTENT(in)      :: periodic
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: gdim
    INTEGER, DIMENSION(2, 0:2), INTENT(in), &
      OPTIONAL                               :: local_bounds
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: local_shift
    REAL(dp), DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: poly_shift
    REAL(dp), INTENT(in), OPTIONAL           :: scale
    TYPE(lgrid_type), INTENT(inout), &
      OPTIONAL                               :: lgrid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'collocGauss', &
      routineP = moduleN//':'//routineN

#include "colloc_int_body.f90"

END SUBROUTINE


! *****************************************************************************
!> \brief integrates a gaussian times any polynomial up to a give order.
!>
!> Most things are the same as for collocGauss (see its comments).
!> Returns the integrals of all the monomials in d3 format into poly
!> \param h ...
!> \param h_inv ...
!> \param grid ...
!> \param poly ...
!> \param alphai ...
!> \param posi ...
!> \param max_r2 ...
!> \param periodic ...
!> \param gdim ...
!> \param local_bounds ...
!> \param local_shift ...
!> \param poly_shift ...
!> \param scale ...
!> \param error ...
! *****************************************************************************
SUBROUTINE integrateGaussFull(h,h_inv,grid,poly,alphai,posi,max_r2,&
        periodic,gdim,local_bounds,local_shift,poly_shift,scale,error)
    REAL(dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: h, h_inv
    REAL(dp), DIMENSION(0:, 0:, 0:), &
      INTENT(inout)                          :: grid
    REAL(dp), DIMENSION(:), INTENT(out)      :: poly
    REAL(dp), INTENT(in)                     :: alphai
    REAL(dp), DIMENSION(0:2), INTENT(in)     :: posi
    REAL(dp), INTENT(in)                     :: max_r2
    INTEGER, DIMENSION(0:2), INTENT(in)      :: periodic
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: gdim
    INTEGER, DIMENSION(2, 0:2), INTENT(in), &
      OPTIONAL                               :: local_bounds
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: local_shift
    REAL(dp), DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: poly_shift
    REAL(dp), INTENT(in), OPTIONAL           :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrateGaussFull', &
        routineP=moduleN//':'//routineN

#define FMG_INTEGRATE_FULL
#include "colloc_int_body.f90"
#undef FMG_INTEGRATE_FULL
END SUBROUTINE

END MODULE gauss_colloc
