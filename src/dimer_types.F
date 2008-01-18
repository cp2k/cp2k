!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/dimer_types [1.0] *
!!
!!   NAME
!!     dimer_types
!!
!!   FUNCTION
!!     Contains types used for a Dimer Method calculations
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
MODULE dimer_types
  USE f77_blas
  USE global_types,                    ONLY: global_environment_type
  USE gopt_param_types,                ONLY: gopt_param_type
  USE kinds,                           ONLY: dp
  USE parallel_rng_types,              ONLY: random_numbers
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dimer_types'
  INTEGER, PRIVATE, SAVE :: last_dimer_id=0

  PUBLIC :: dimer_env_type
  PUBLIC :: dimer_env_create, dimer_env_retain,&
            dimer_env_release

!!****s* dimer_types/dimer_env_type [1.0] *
!!
!!   NAME
!!     dimer_env_type
!!
!!   FUNCTION
!!     Defines the environment for a Dimer Method calculation
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
  !***************************************************************************
  TYPE dimer_env_type
     INTEGER :: ref_count, id_nr
     LOGICAL                                    :: rotate_dimer
     LOGICAL                                    :: interpolate_gradient
     REAL(KIND=dp)                              :: dr
     REAL(KIND=dp), POINTER, DIMENSION(:)       :: nvec, g0, g1, g1p
     REAL(KIND=dp)                              :: angle, dCdp, curvature
  END TYPE dimer_env_type

CONTAINS
  
!!****f* dimer_types/dimer_env_create [1.0] *
!!
!!   NAME
!!     dimer_env_create
!!
!!   FUNCTION
!!     -
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
  SUBROUTINE dimer_env_create(dimer_env, gopt_param, natom, globenv, error)
    TYPE(dimer_env_type), POINTER            :: dimer_env
    TYPE(gopt_param_type), POINTER           :: gopt_param
    INTEGER, INTENT(IN)                      :: natom
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dimer_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: norm

    failure=.FALSE.
    CPPostcondition(.NOT.ASSOCIATED(dimer_env),cp_failure_level,routineP,error,failure)
    ALLOCATE(dimer_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       dimer_env%ref_count=1
       last_dimer_id=last_dimer_id+1
       dimer_env%id_nr=last_dimer_id
       ! Setup NVEC
       NULLIFY(dimer_env%nvec)
       ! Allocate the working arrays
       ALLOCATE(dimer_env%nvec(natom*3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(dimer_env%g0(natom*3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(dimer_env%g1(natom*3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(dimer_env%g1p(natom*3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL random_numbers(dimer_env%nvec, globenv%gaussian_rng_stream, error)
       norm = SQRT(SUM(dimer_env%nvec**2))
       CPPostcondition(norm>0.0_dp,cp_failure_level,routineP,error,failure)
       dimer_env%nvec         = dimer_env%nvec / norm
       dimer_env%rotate_dimer = .TRUE.
    END IF
  END SUBROUTINE dimer_env_create
!***************************************************************************

!!****f* dimer_types/dimer_env_retain [1.0] *
!!
!!   NAME
!!     dimer_env_retain
!!
!!   FUNCTION
!!     -
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
  SUBROUTINE dimer_env_retain(dimer_env, error)
    TYPE(dimer_env_type), POINTER            :: dimer_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dimer_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(dimer_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPreconditionNoFail(dimer_env%ref_count>0,cp_failure_level,routineP,error)
       dimer_env%ref_count=dimer_env%ref_count+1
    END IF
  END SUBROUTINE dimer_env_retain
!***************************************************************************


!!****f* dimer_types/dimer_env_release [1.0] *
!!
!!   NAME
!!     dimer_env_release
!!
!!   FUNCTION
!!     -
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
  SUBROUTINE dimer_env_release(dimer_env, error)
    TYPE(dimer_env_type), POINTER            :: dimer_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dimer_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(dimer_env)) THEN
       CPPreconditionNoFail(dimer_env%ref_count>0,cp_failure_level,routineP,error)
       dimer_env%ref_count=dimer_env%ref_count-1
       IF (dimer_env%ref_count==0) THEN
          IF (ASSOCIATED(dimer_env%nvec)) THEN
             DEALLOCATE(dimer_env%nvec, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(dimer_env%g0)) THEN
             DEALLOCATE(dimer_env%g0, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(dimer_env%g1)) THEN
             DEALLOCATE(dimer_env%g1, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(dimer_env%g1p)) THEN
             DEALLOCATE(dimer_env%g1p, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          DEALLOCATE(dimer_env, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE dimer_env_release
!***************************************************************************

END MODULE dimer_types
