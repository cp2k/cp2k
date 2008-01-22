!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/dimer_utils [1.0] *
!!
!!   NAME
!!     dimer_utils
!!
!!   FUNCTION
!!     Contains utilities for a Dimer Method calculations
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
MODULE dimer_utils
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dimer_utils'

  PUBLIC :: rotate_dimer, dimer_thrs
  REAL(KIND=dp), PARAMETER                 :: dimer_thrs=EPSILON(0.0_dp)*1.0E4_dp

CONTAINS

!!****f* dimer_utils/rotate_dimer [1.0] *
!!
!!   NAME
!!     rotate_dimer
!!
!!   FUNCTION
!!     Performs a rotation of the unit dimer vector
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     -
!!
!!   AUTHOR
!!     Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE rotate_dimer(nvec,theta,dt,error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: nvec, theta
    REAL(KIND=dp)                            :: dt
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotate_dimer', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check, failure

    failure = .FALSE.
    ! Orthogonality check for the rotation..
    check = ABS(DOT_PRODUCT(nvec,theta))<dimer_thrs
    CPPrecondition(check,cp_failure_level,routineP,error,failure)
    nvec  = nvec * COS(dt) + theta * SIN(dt)

  END SUBROUTINE rotate_dimer

END MODULE dimer_utils
