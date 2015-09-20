!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE pint_transformations
  USE input_constants,                 ONLY: transformation_stage
  USE kinds,                           ONLY: dp
  USE pint_normalmode,                 ONLY: normalmode_f2uf,&
                                             normalmode_u2x,&
                                             normalmode_x2u
  USE pint_staging,                    ONLY: staging_f2uf,&
                                             staging_u2x,&
                                             staging_x2u
  USE pint_types,                      ONLY: pint_env_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_transformations'

  PUBLIC :: pint_x2u, &
            pint_u2x, &
            pint_f2uf

CONTAINS

! ***************************************************************************
!> \brief Transforms from the x into the u variables
!>      (at the moment a staging transformation for the positions)
!> \param pint_env the path integral environment
!> \param ux will contain the u variable (defaults to pint_env%ux)
!> \param x the positions to transform (defaults to pint_env%x)
!> \par History
!>      Added normal mode transformation [hforbert]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_x2u(pint_env,ux,x)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL, TARGET          :: ux
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: x

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_x2u', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_ux, my_x

    failure=.FALSE.

    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
    my_x => pint_env%x
    my_ux => pint_env%ux
    IF (PRESENT(x)) my_x => x
    IF (PRESENT(ux)) my_ux => ux
    CPASSERT(ASSOCIATED(my_ux))
    CPASSERT(ASSOCIATED(my_x))

    IF (pint_env%transform == transformation_stage) THEN
      CALL staging_x2u(pint_env%staging_env,ux=my_ux,x=my_x)
    ELSE
      CALL normalmode_x2u(pint_env%normalmode_env,ux=my_ux,x=my_x)
    END IF
    RETURN
  END SUBROUTINE pint_x2u

! ***************************************************************************
!> \brief transform from the u variable to the x (inverse of x2u)
!> \param pint_env path integral environment
!> \param ux the u variable (positions to be backtransformed)
!> \param x will contain the positions
!> \par History
!>      Added normal mode transformation by hforbert
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_u2x(pint_env,ux,x)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: ux
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL, TARGET          :: x

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_u2x', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_ux, my_x

    failure=.FALSE.

    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
    my_x => pint_env%x
    my_ux => pint_env%ux
    IF (PRESENT(x)) my_x => x
    IF (PRESENT(ux)) my_ux => ux
    CPASSERT(ASSOCIATED(my_ux))
    CPASSERT(ASSOCIATED(my_x))

    IF (pint_env%transform == transformation_stage) THEN
      CALL staging_u2x(pint_env%staging_env,ux=my_ux,x=my_x)
    ELSE
      CALL normalmode_u2x(pint_env%normalmode_env,ux=my_ux,x=my_x)
    END IF
    RETURN
  END SUBROUTINE pint_u2x

! ***************************************************************************
!> \brief transformation x to u for the forces
!> \param pint_env the path integral environment
!> \param uf will contain the accelerations for the transformed variables
!>        afterwards
!> \param f the forces to transform
!> \par History
!>      Added normal mode transformation [hforbert]
!>      Divide forces by the number of beads, since the replication
!>        environment (should) give raw forces [hforbert]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_f2uf(pint_env,uf,f)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL, TARGET          :: uf
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: f

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_f2uf', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_f, my_uf

    failure=.FALSE.

    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
    my_f => pint_env%f
    my_uf => pint_env%uf
    IF (PRESENT(f)) my_f => f
    IF (PRESENT(uf)) my_uf => uf
    CPASSERT(ASSOCIATED(my_uf))
    CPASSERT(ASSOCIATED(my_f))

    IF (pint_env%transform == transformation_stage) THEN
      CALL staging_f2uf(pint_env%staging_env,uf=my_uf,f=my_f)
    ELSE
      CALL normalmode_f2uf(pint_env%normalmode_env,uf=my_uf,f=my_f)
    END IF
    my_uf=my_uf/pint_env%mass_fict/REAL(pint_env%p,dp)
    RETURN
  END SUBROUTINE pint_f2uf

END MODULE pint_transformations
