!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition of the semi empirical parameter types.
!> \author Teodoro Laino [tlaino] - 10.2008 University of Zurich
! *****************************************************************************
MODULE taper_types
  
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'taper_types'

! *****************************************************************************
!> \brief Taper type
! *****************************************************************************
  TYPE taper_type
     LOGICAL                               :: apply_taper
     REAL(KIND=dp)                         :: r0, rscale
  END TYPE taper_type

  PUBLIC :: taper_type, taper_create, taper_release, taper_eval, dtaper_eval

CONTAINS

! *****************************************************************************
!> \brief Creates taper type
!> \param taper ...
!> \param rc ...
!> \param range ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE taper_create(taper, rc, range, error)
    TYPE(taper_type), POINTER                :: taper
    REAL(KIND=dp), INTENT(IN)                :: rc, range
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'taper_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(taper),cp_failure_level,routineP,error,failure)
    ALLOCATE (taper,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (range > EPSILON(0.0_dp)) THEN
       taper%apply_taper = .TRUE.
       CPPrecondition(range>0.0_dp,cp_failure_level,routineP,error,failure)
       taper%r0     = 2.0_dp*rc - 20.0_dp * range
       taper%rscale = 1.0_dp/range
    ELSE
       taper%apply_taper = .FALSE.
    END IF

  END SUBROUTINE taper_create

! *****************************************************************************
!> \brief Releases taper type
!> \param taper ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE taper_release(taper, error)
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'taper_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(taper)) THEN
       DEALLOCATE (taper,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE taper_release

! *****************************************************************************
!> \brief Taper functions
!> \param taper ...
!> \param rij ...
!> \param error ...
!> \retval ft ...
! *****************************************************************************
  FUNCTION taper_eval (taper, rij, error) RESULT(ft)
    TYPE(taper_type), POINTER                :: taper
    REAL(KIND=dp), INTENT(IN)                :: rij
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: ft

    CHARACTER(len=*), PARAMETER :: routineN = 'taper_eval', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: dr

    ft = 1._dp
    IF (taper%apply_taper) THEN
       dr = taper%rscale*(rij-taper%r0)
       ft = 0.5_dp*(1.0_dp-TANH(dr))
    END IF
  END FUNCTION taper_eval

! *****************************************************************************
!> \brief Analytical derivatives for taper function
!> \param taper ...
!> \param rij ...
!> \param error ...
!> \retval dft ...
! *****************************************************************************
  FUNCTION dtaper_eval (taper, rij, error ) RESULT(dft)
    TYPE(taper_type), POINTER                :: taper
    REAL(KIND=dp), INTENT(IN)                :: rij
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: dft

    CHARACTER(len=*), PARAMETER :: routineN = 'dtaper_eval', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: dr

    dft = 0.0_dp
    IF (taper%apply_taper) THEN
       dr  = taper%rscale*(rij-taper%r0)
       dft = -0.5_dp*(1.0_dp-TANH(dr)**2)*taper%rscale
    END IF
  END FUNCTION dtaper_eval

END MODULE taper_types
