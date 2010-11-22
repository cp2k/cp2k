!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Splining replica for band calculations
!> \author Teodoro Laino 10.2006
! *****************************************************************************
MODULE spline_env_replica_types
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'spline_env_replica_types'

  PUBLIC :: spline_env_replica_type,&
            create_spline_env_replica,&
            release_spline_env_replica


! *****************************************************************************
  TYPE spline_env_replica_type
     LOGICAL                                   :: periodic
     INTEGER                                   :: order
     INTEGER                                   :: nt! This is the dimension of the
                                                    ! knot-vector 
     INTEGER                                   :: np! This is the number data points
     INTEGER                                   :: ns! Dimension of the Spline family or
                                                    ! number of control points
     REAL(KIND=dp), DIMENSION(:, :),  POINTER  :: C ! Control points
     REAL(KIND=dp), DIMENSION(:),  POINTER     :: t ! knot-vector
  END TYPE spline_env_replica_type

CONTAINS

! *****************************************************************************
!> \param  Creates a spline_env environment for splining replicas
!> \author Teodoro Laino [tlaino] 
!> \date   06.2009
! *****************************************************************************
  SUBROUTINE create_spline_env_replica(spline_env_replica, spline_order, nreplica,&
       nxyz, error)
    TYPE(spline_env_replica_type), POINTER   :: spline_env_replica
    INTEGER, INTENT(IN)                      :: spline_order, nreplica, nxyz
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_spline_env_replica', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: d1, d2, stat
    LOGICAL                                  :: check, failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       check = .NOT.ASSOCIATED(spline_env_replica)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       ALLOCATE(spline_env_replica, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       spline_env_replica%periodic = .FALSE. ! For the time being hardcoded to non-periodic
                                             ! The periodic case will be easy to handle
       ! Spline Order > 0
       CALL cp_assert(spline_order>0,&
            cp_fatal_level,cp_assertion_failed,routineP,&
            "Spline order for band methods must be larger than 0!"//&
 CPSourceFileRef,&
            error=error,failure=failure)
       spline_env_replica%order    = spline_order
       IF (spline_env_replica%periodic) THEN
          d1 = nreplica
          d2 = 0
       ELSE
          d1 = nreplica + 2*spline_order - 1
          d2 = spline_order - 1
       END IF
       spline_env_replica%nt       = d1
       spline_env_replica%np       = nreplica
       spline_env_replica%ns       = (nreplica-1)+d2
       ALLOCATE(spline_env_replica%C(nxyz, spline_env_replica%ns),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(spline_env_replica%t(spline_env_replica%nt),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF    
  END SUBROUTINE create_spline_env_replica

! *****************************************************************************
!> \param  Releases a spline_env environment for splining replicas
!> \author Teodoro Laino [tlaino] 
!> \date   06.2009
! *****************************************************************************
  SUBROUTINE release_spline_env_replica(spline_env_replica, error)
    TYPE(spline_env_replica_type), POINTER   :: spline_env_replica
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_spline_env_replica', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: check, failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       check = ASSOCIATED(spline_env_replica)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       DEALLOCATE(spline_env_replica%C,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(spline_env_replica%t,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(spline_env_replica, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF    
  END SUBROUTINE release_spline_env_replica

END MODULE spline_env_replica_types
