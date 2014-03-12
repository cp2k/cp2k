!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition of the type to handle the 1/R^3 residual integral part
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
MODULE semi_empirical_expns3_types
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

! *****************************************************************************
!> \brief 1/R^3 expansion type
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
  TYPE semi_empirical_expns3_type
     REAL(KIND=dp)                                  :: core_core
     REAL(KIND=dp), DIMENSION(9)                    :: e1b, e2a
     REAL(KIND=dp), DIMENSION(81)                   :: w
  END TYPE semi_empirical_expns3_type

! *****************************************************************************
!> \brief 1/R^3 expansion type: array of pointers
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
  TYPE semi_empirical_expns3_p_type
     TYPE(semi_empirical_expns3_type), POINTER      :: expns3
  END TYPE semi_empirical_expns3_p_type

  ! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_expns3_types'

  PUBLIC :: semi_empirical_expns3_type,&
            semi_empirical_expns3_p_type,&
            semi_empirical_expns3_create,&
            semi_empirical_expns3_release

CONTAINS

! *****************************************************************************
!> \brief Allocate semi-empirical 1/R^3 expansion type
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
  SUBROUTINE semi_empirical_expns3_create(expns3, error)
    TYPE(semi_empirical_expns3_type), &
      POINTER                                :: expns3
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_expns3_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(expns3),cp_failure_level,routineP,error,failure)
    ALLOCATE (expns3,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    expns3%core_core = 0.0_dp
    expns3%e1b       = 0.0_dp
    expns3%e2a       = 0.0_dp
    expns3%w         = 0.0_dp
  END SUBROUTINE semi_empirical_expns3_create

! *****************************************************************************
!> \brief Deallocate the semi-empirical type
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
  SUBROUTINE semi_empirical_expns3_release(expns3, error)
    TYPE(semi_empirical_expns3_type), &
      POINTER                                :: expns3
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'semi_empirical_expns3_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(expns3)) THEN
       DEALLOCATE (expns3,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE semi_empirical_expns3_release

END MODULE semi_empirical_expns3_types
