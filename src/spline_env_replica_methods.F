!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Methods working on spline replica for bands calculations
!> \author Teodoro Laino 10.2006
! *****************************************************************************
MODULE spline_env_replica_methods
  USE band_replica_utils,              ONLY: neb_replica_distance
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: invert_matrix
  USE neb_types,                       ONLY: neb_var_type
  USE spline_env_replica_types,        ONLY: spline_env_replica_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'spline_env_replica_methods'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

  PUBLIC :: eval_spline_replica, reparametrize_spline_replica

CONTAINS

! *****************************************************************************
!> \brief  Evaluates the spline env for the replica when performing the string
!>         method
!> \author Teodoro Laino [tlaino] - 06.2009
! *****************************************************************************
  SUBROUTINE eval_spline_replica(spline_env_replica, coords, nreplica, nxyz, error)
    TYPE(spline_env_replica_type), POINTER   :: spline_env_replica
    TYPE(neb_var_type), POINTER              :: coords
    INTEGER, INTENT(IN)                      :: nreplica, nxyz
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eval_spline_replica', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dim, i, j, k, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: inv_error, l_tot, xmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: distances
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: fmp, imp, qspl, tmp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: t

    failure = .FALSE.
    ! Allocate working arrays
    ALLOCATE(distances(nreplica-1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)    
    t => spline_env_replica%t
    ! First evaluate chord length along the B-spline
    DO i = 2, nreplica
       CALL neb_replica_distance(coords=coords, i0=i-1, i=i, distance=distances(i-1),&
            iw=-1, rotate=.FALSE., error=error)
    END DO
    l_tot = SUM(distances)

    dim   = 0
    IF (.NOT.spline_env_replica%periodic) dim = spline_env_replica%order-1 

    ! The parameter value t goes from 0 to 1.. (keeping into account boundary 
    ! conditions for non-periodic m-splines)
    t(1:dim+1) = 0.0_dp
    DO j = 2, nreplica
       i = j + dim
       t(i) = t(i-1) + distances(j-1)*(1.0_dp/l_tot)
    END DO
    ! Fix the last value to be 1.0_dp to avoid possible numerical errors..
    t(i) = 1.0_dp
    IF (.NOT.spline_env_replica%periodic) t(i:i+spline_env_replica%order) = t(i)
    
    ! Evaluate B-Spline Matrix for computing control points
    ALLOCATE(qspl(spline_env_replica%np, spline_env_replica%ns),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i =1, spline_env_replica%np
       qspl(i,:) = spline_low(spline_env_replica%periodic,spline_env_replica%nt,&
                              spline_env_replica%ns,t(dim+i),t,spline_env_replica%order,&
                              error)
    END DO

    ! Debug statements
    IF (debug_this_module) THEN
       WRITE(*,*)spline_env_replica%nt,t
       l_tot = 0.0_dp
       DO i =1, 101
          qspl(1,:) = spline_low(spline_env_replica%periodic,spline_env_replica%nt,&
                                 spline_env_replica%ns,l_tot,t,spline_env_replica%order,&
                                 error)
          WRITE(10,'(20F9.5)')l_tot,qspl(1,1:spline_env_replica%ns)
          l_tot=MIN(1.0_dp,l_tot+0.01_dp)
       END DO
       STOP
    END IF
    
    ! Find control points
    ALLOCATE(tmp(spline_env_replica%ns,spline_env_replica%ns),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(imp(spline_env_replica%ns,spline_env_replica%ns),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(fmp(spline_env_replica%ns,spline_env_replica%np),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    WRITE(*,'(7f12.5)')((qspl(i,j),j=1,spline_env_replica%ns),i=1,spline_env_replica%np)
    WRITE(*,*)
    tmp = MATMUL(TRANSPOSE(qspl),qspl)
    WRITE(*,*)spline_env_replica%ns
    WRITE(*,'(7f12.5)')((tmp(i,j),j=1,spline_env_replica%ns),i=1,spline_env_replica%ns)
    CALL invert_matrix(tmp, imp, inv_error, error=error)
    WRITE(*,*)inv_error
    fmp=MATMUL(tmp,TRANSPOSE(qspl))

    WRITE(*,*)
    WRITE(*,'(7f12.5)')MATMUL(imp,tmp)
    WRITE(*,*)
    DO i = 1, spline_env_replica%ns
       DO j = 1, nxyz
          xmp = 0.0_dp
          DO k = 1, spline_env_replica%np
             xmp = xmp + fmp(i,k)*coords%wrk(j,k)
          END DO
          spline_env_replica%C(j, i) = xmp
       END DO
    END DO
    STOP
    !TL stop here..

    DEALLOCATE(imp,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tmp,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Release working arrays
    DEALLOCATE(qspl,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(distances,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE eval_spline_replica

! *****************************************************************************
!> \brief  Evaluate spline (bspline or mspline) of order k along at point t 
!>         along all knots tarray
!> \author Teodoro Laino [tlaino] - 06.2009
! *****************************************************************************
  FUNCTION spline_low(periodic, nt, ns, t, tarray, k, error) RESULT(bspline)
    LOGICAL, INTENT(IN)                      :: periodic
    INTEGER, INTENT(IN)                      :: nt, ns
    REAL(KIND=dp), INTENT(IN)                :: t
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: tarray
    INTEGER, INTENT(IN)                      :: k
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp), DIMENSION(ns)             :: bspline

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure      = .FALSE.
    ! Control that t is contained in the tarray
    CPPostcondition(t>=tarray( 1),cp_failure_level,routineP,error,failure)
    CPPostcondition(t<=tarray(nt),cp_failure_level,routineP,error,failure)

    IF (periodic) THEN
       ! Not yet implemented
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ELSE
       DO i = 1, ns
          bspline(i) = mspline(t,k,i,tarray,error)
       END DO
    END IF

  END FUNCTION spline_low

! *****************************************************************************
!> \brief  Evaluate the kth order m-spline at t
!> \Notes  Arguments:
!>         t : the point at which the spline is evaluated
!>         k : the order of the spline (k>=1)
!>         i : the interval on which the spline is defined (i<LEN(tarray)-k)
!>         tarray: the knots
!>
!> \author Teodoro Laino [tlaino] - 07.2009
! *****************************************************************************
  RECURSIVE FUNCTION mspline(t,k,i,tarray,error) RESULT(res)
    REAL(KIND=dp), INTENT(IN)                :: t
    INTEGER, INTENT(IN)                      :: k, i
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: tarray
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'mspline', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: d1, d2, d3, spl

    failure = .FALSE.

    IF (t< tarray(i)  ) THEN 
       res = 0.0_dp
       RETURN
    END IF
    IF (t>tarray(i+k)) THEN
       res = 0.0_dp
       RETURN
    END IF

    IF (k==1) THEN
       d1 = tarray(i+1)-tarray(i)
       res = 0.0_dp
       IF (d1/=0.0_dp) res = 1.0_dp/(tarray(i+1)-tarray(i))
       RETURN
    END IF
    
    d3 = tarray(i+k)-tarray(i)
    IF (d3 /= 0.0_dp) THEN
       d1 = t-tarray(i)
       d2 = tarray(i+k)-t
       spl= d1*mspline(t,k-1,i,tarray,error)+d2*mspline(t,k-1,i+1,tarray,error)
       res= spl*REAL(k,KIND=dp)/(REAL(k-1,KIND=dp)*d3)
    ELSE
       res= 0.0_dp
    END IF
  END FUNCTION mspline

! *****************************************************************************
!> \brief  Reparametrize replica (for string method) using splines order N
!> \author Teodoro Laino [tlaino] - 06.2009
! *****************************************************************************
  SUBROUTINE reparametrize_spline_replica(spline_env_replica, coords, distances,&
       nreplica, nxyz, error)
    TYPE(spline_env_replica_type), POINTER   :: spline_env_replica
    TYPE(neb_var_type), POINTER              :: coords
    REAL(KIND=dp), DIMENSION(:)              :: distances
    INTEGER, INTENT(IN)                      :: nreplica, nxyz
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reparametrize_spline_replica', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: avg_distance, xtmp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: tmp_coords

    failure = .FALSE.

    ALLOCATE(tmp_coords(nxyz,nreplica),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    tmp_coords = coords%wrk

    ! Reparametrization
    SELECT CASE(spline_env_replica%order)
    CASE(1)
       ! Compute distances
       DO i = 2, nreplica
          CALL neb_replica_distance(coords=coords, i0=i-1, i=i, distance=distances(i-1),&
               iw=-1, rotate=.FALSE., error=error)
       END DO
       avg_distance = SUM(distances)/REAL(nreplica-1,KIND=dp)
       ! Redistribute frames
       DO i = 2, nreplica-1
          xtmp  = 0.0_dp
          DO j = 1, nreplica-1
             xtmp = xtmp + distances(j)
             IF (xtmp>avg_distance*REAL(i-1,KIND=dp)) THEN
                xtmp = (xtmp - avg_distance*REAL(i-1,KIND=dp))/distances(j)
                coords%wrk(:,i) = (1.0_dp-xtmp)*tmp_coords(:,j+1) + xtmp*tmp_coords(:,j)
                EXIT
             END IF
          END DO
       END DO
       ! Re-compute distances
       DO i = 2, nreplica
          CALL neb_replica_distance(coords=coords, i0=i-1, i=i, distance=distances(i-1),&
               iw=-1, rotate=.FALSE., error=error)
       END DO
    CASE DEFAULT
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="String Method: Spline order greater than 1 not implemented.", &
            error=error, error_level=cp_warning_level)
    END SELECT

    DEALLOCATE(tmp_coords,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE reparametrize_spline_replica

END MODULE spline_env_replica_methods
