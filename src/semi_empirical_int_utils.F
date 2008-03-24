!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for Integrals for semi-empiric methods
!> \author Teodoro Laino (03.2008) [tlaino] - University of Zurich 
! *****************************************************************************
MODULE semi_empirical_int_utils

  USE f77_blas
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_int_utils'

  PUBLIC :: set_taper_fn, taper, dtaper_ana

  REAL(KIND=dp), SAVE, PUBLIC          :: r0, al
  LOGICAL, SAVE, PUBLIC                :: taper_fn_init = .FALSE.
  REAL(KINd=dp), PARAMETER, PUBLIC     :: rij_threshold = 0.00002_dp
  
  ! Parameters for multipole sums
  REAL(KIND=dp), DIMENSION(72), PARAMETER, PUBLIC  :: pxy = (/&
  !    1          2          3          4         5         6           7          8          9         10
     1.0_dp,    0.5_dp,   -0.5_dp,   0.25_dp,  0.25_dp,   -0.5_dp,    0.5_dp,    0.5_dp,   -0.5_dp,   0.25_dp,& ! 0
    0.25_dp,   -0.5_dp,    0.5_dp,    0.5_dp,  -0.5_dp,   0.25_dp,   0.25_dp,  -0.25_dp,  -0.25_dp,   0.25_dp,& ! 1
   -0.25_dp,  -0.25_dp,   0.25_dp,   0.25_dp, -0.25_dp,  -0.25_dp,   0.25_dp, -0.125_dp,  0.125_dp, -0.125_dp,& ! 2
   0.125_dp, -0.125_dp,  0.125_dp, -0.125_dp, 0.125_dp,   0.25_dp,  0.125_dp,  0.125_dp,  -0.25_dp,  -0.25_dp,& ! 3
    0.25_dp, -0.125_dp,  0.125_dp, -0.125_dp, 0.125_dp, -0.125_dp,  0.125_dp, -0.125_dp,  0.125_dp, 0.0625_dp,& ! 4
  0.0625_dp, 0.0625_dp, 0.0625_dp,  -0.25_dp,  0.25_dp,   0.25_dp,  -0.25_dp,  -0.25_dp,   0.25_dp,   0.25_dp,& ! 5
   -0.25_dp,    0.0_dp,    0.0_dp,    0.0_dp, 0.125_dp, -0.125_dp, -0.125_dp,  0.125_dp, -0.125_dp,  0.125_dp,& ! 6
   0.125_dp, -0.125_dp/)

CONTAINS

! *****************************************************************************
!> \brief Initialization of taper functions
! *****************************************************************************
  SUBROUTINE set_taper_fn ( rc, range, error)
    REAL(KIND=dp), INTENT(IN)                :: rc, range
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_taper_fn', &
      routineP = moduleN//':'//routineN

    taper_fn_init = .TRUE.
    r0 = rc - 0.5_dp * range
    al = 24._dp/range

  END SUBROUTINE set_taper_fn

! *****************************************************************************
!> \brief Taper functions
! *****************************************************************************
  FUNCTION taper (rij, error) RESULT(ft)
    REAL(KIND=dp), INTENT(IN)                :: rij
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: ft

    CHARACTER(len=*), PARAMETER :: routineN = 'taper', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: dr, ep

    ft = 1._dp
    IF ( taper_fn_init ) THEN
       dr = -al*(rij-r0)
       IF ( dr < 20._dp ) THEN
          ep = EXP(dr)
          ft = ep/(1._dp+ep)
       END IF
    END IF
  END FUNCTION taper

! *****************************************************************************
!> \brief Analytical derivatives for taper function
! *****************************************************************************
  FUNCTION dtaper_ana ( rij, error ) RESULT(dft)
    REAL(KIND=dp), INTENT(IN)                :: rij
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: dft

    CHARACTER(len=*), PARAMETER :: routineN = 'dtaper_ana', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: dr, ep

    dft = 0.0_dp
    IF ( taper_fn_init ) THEN
       dr = -al*(rij-r0)
       IF ( dr < 20._dp ) THEN
          ep = EXP(dr)
          dft = 1._dp/(1._dp+ep)**2 * ep * (-al)
       END IF
    END IF

  END FUNCTION dtaper_ana

END MODULE semi_empirical_int_utils
