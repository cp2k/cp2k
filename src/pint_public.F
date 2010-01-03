!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Public path integral routines that can be called from other modules
!> \author Lukasz Walewski
!> \date   2009-07-24
!> \note   Avoiding circular dependencies: please design new members of this
!>         module in such a way that they use pint_types module only.
! *****************************************************************************
MODULE pint_public

  USE kinds,                           ONLY: dp
  USE pint_types,                      ONLY: pint_env_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_public'

  PUBLIC :: pint_com_pos

  CONTAINS

  ! ***************************************************************************
  !> \brief  Return the center of mass of the PI system
  !> \author Lukasz Walewski
  !> \date   2009-07-24
  !> \par    History
  !>           2009-11-30 fixed serious bug in pint_env%x indexing [lwalewski]
  ! ***************************************************************************
  FUNCTION pint_com_pos(pint_env,error) RESULT(com_r)

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    REAL(kind=dp), DIMENSION(3)              :: com_r

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_com_pos', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib, ic
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: tmass

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)

    tmass = 0.0_dp
    com_r(:) = 0.0_dp
    DO ia = 1, pint_env%ndim/3
      DO ib = 1, pint_env%p
        DO ic = 1, 3
          com_r(ic) = com_r(ic) + &
            pint_env%x(ib,(ia-1)*3+ic) * pint_env%mass((ia-1)*3+ic)
          tmass = tmass + pint_env%mass((ia-1)*3+ic)
        END DO
      END DO
    END DO
    ! pint_env%mass is REAL, DIMENSION(NDIM) which means that each atom
    ! has it's mass defined three times - here we hope that all three
    ! values are equal
    tmass = tmass / 3.0_dp
    com_r(:) = com_r(:) / tmass
    RETURN
  END FUNCTION

  ! ***************************************************************************
  !> \brief  Return the center of geometry of the PI system
  !> \author Lukasz Walewski
  !> \date   2009-11-30
  ! ***************************************************************************
  FUNCTION pint_cog_pos(pint_env,error) RESULT(cntrd_r)

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    REAL(kind=dp), DIMENSION(3)              :: cntrd_r

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_cog_pos', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib, ic, natoms
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)

    cntrd_r(:) = 0.0_dp
    natoms = pint_env%ndim/3
    DO ia = 1, natoms
      DO ib = 1, pint_env%p
        DO ic = 1, 3
          cntrd_r(ic) = cntrd_r(ic) + pint_env%x(ib,(ia-1)*3+ic)
        END DO
      END DO
    END DO
    cntrd_r(:) = cntrd_r(:) / REAL(pint_env%p,dp) / REAL(natoms,dp)
    RETURN
  END FUNCTION pint_cog_pos

END MODULE pint_public
