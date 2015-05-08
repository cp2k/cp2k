!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief helper routines for cp_dbcsr
!> \par History
!>       2012.07 created [Florian Thoele]
!> \author Florian Thoele
! *****************************************************************************

MODULE cp_dbcsr_util
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_distribution,&
                                             cp_dbcsr_get_data_type,&
                                             cp_dbcsr_get_info,&
                                             cp_dbcsr_type,&
                                             dbcsr_distribution_mp,&
                                             dbcsr_mp_group,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_multiply_local
  USE kinds,                           ONLY: dp,&
                                             sp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_environ,&
                                             mp_max,&
                                             mp_min,&
                                             mp_sum
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_util'

  PUBLIC :: lanczos_alg_serial

CONTAINS

! *****************************************************************************
!> \brief Computes the highest and lowest eigenvalue of matrix_a (assumed to be symmetric)
!>        with the Lanczos algorithm. The 3 Lanczos vectors are replicated.
!>        speedups could be gained from an improved dbcsr x vector multiplication (distributed vector).
!>        if svd is true the square of the matrix is employed.
!> \param matrix_a ...
!> \param max_ev ...
!> \param min_ev ...
!> \param converged ...
!> \param threshold ...
!> \param max_iter ...
!> \param svd ...
!> \param error ...
!> \par History
!>       2012.07 created [Florian Thoele]
!> \author Florian Thoele 
! *****************************************************************************
  SUBROUTINE lanczos_alg_serial(matrix_a, max_ev, min_ev, converged, threshold, max_iter, svd, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(KIND=dp), INTENT(OUT)               :: max_ev, min_ev
    LOGICAL, INTENT(OUT)                     :: converged
    REAL(KIND=dp), INTENT(IN)                :: threshold
    INTEGER, INTENT(IN)                      :: max_iter
    LOGICAL, INTENT(IN), OPTIONAL            :: svd
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_alg_serial', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, mp_group, unit_nr
    LOGICAL                                  :: do_svd
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    mp_group = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_a)))

    do_svd = .FALSE.
    IF(PRESENT(svd)) do_svd = svd

    SELECT CASE(cp_dbcsr_get_data_type(matrix_a))
      CASE (dbcsr_type_real_8)
        CALL lanczos_alg_serial_d(matrix_a, max_ev, min_ev, converged, threshold, max_iter, do_svd, mp_group, error)
      CASE (dbcsr_type_real_4)
        CALL lanczos_alg_serial_s(matrix_a, max_ev, min_ev, converged, threshold, max_iter, do_svd, mp_group, error)
    !  CASE DEFAULT
    !    CALL stop_program(routineN,moduleN,__LINE__,"NYI")
    END SELECT

    CALL timestop(handle)
  END SUBROUTINE lanczos_alg_serial

! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param max_ev ...
!> \param min_ev ...
!> \param converged ...
!> \param threshold ...
!> \param max_iter ...
!> \param svd ...
!> \param mp_group ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE lanczos_alg_serial_d(matrix_a, max_ev, min_ev, converged, threshold, max_iter, svd, mp_group, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(KIND=dp), INTENT(OUT)               :: max_ev, min_ev
    LOGICAL, INTENT(OUT)                     :: converged
    REAL(KIND=dp), INTENT(IN)                :: threshold
    INTEGER, INTENT(IN)                      :: max_iter
    LOGICAL, INTENT(IN)                      :: svd
    INTEGER                                  :: mp_group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_alg_serial_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, info, lwork, n, ntask, rank
    REAL(KIND=dp)                            :: alpha, beta, DDOT, DNRM2, &
                                                nrm_v
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eval, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: evs, evs_tmp, u, u_tmp, v, vn
    TYPE(rng_stream_type), POINTER           :: rng_stream

    converged = .FALSE.

    CALL mp_environ(ntask,rank,mp_group)
    CALL cp_dbcsr_get_info(matrix_a, nfullcols_total=n)
    ALLOCATE(v(n,1), u(n,1), vn(n,1))
    IF(svd) ALLOCATE(u_tmp(n,1))

    ALLOCATE(evs(max_iter+1, max_iter+1), evs_tmp(max_iter+1, max_iter+1))
    lwork = max_iter*3+100
    ALLOCATE(eval(max_iter), work(lwork))

    NULLIFY (rng_stream)
    CALL create_rng_stream(rng_stream,name="lanczos_alg_serial",distribution_type=UNIFORM,&
                           extended_precision=.TRUE.,error=error)
    DO i=1,n
      v(i,1) = next_random_number(rng_stream,error=error) - 0.5_dp
    END DO
    CALL delete_rng_stream(rng_stream,error=error)
    CALL mp_sum(v, mp_group)
    nrm_v = DDOT(n,v, 1, v, 1)
    v(:,:) = v(:,:)/SQRT(nrm_v)

    u(:,:) = 0.0_dp
    IF(svd) THEN
      u_tmp(:,:) = 0.0_dp
      CALL cp_dbcsr_multiply_local(matrix_a, v, u_tmp, error=error)
      CALL mp_sum(u_tmp, mp_group)
      CALL cp_dbcsr_multiply_local(matrix_a, u_tmp, u, error=error)
      CALL mp_sum(u, mp_group)
    ELSE
      CALL cp_dbcsr_multiply_local(matrix_a, v, u, error=error)
      CALL mp_sum(u, mp_group)
    ENDIF

    evs(:,:) = 0.0_dp
    max_ev = 0.0_dp
    min_ev = 0.0_dp

    DO i=1,max_iter
      !Perform computation in serial to avoid numerical divergence
      !of ill-conditioned operations and then broadcast results
      IF (rank .EQ. 0) THEN
        alpha = DDOT(n, u, 1, v, 1)        ! alpha_i = ui^T*vi
        CALL DAXPY(n, -alpha, v, 1, u, 1)  ! ui = ui-alpha*vi   
        beta = DNRM2(n, u, 1)              ! beta_i = sqrt(ri^T*ri)

        evs(i,i) = alpha
        evs(i,i+1) = beta
        evs(i+1,i) = beta
        evs_tmp(:,:) = evs(:,:)
        CALL DSYEV('N', 'U', i, evs_tmp, max_iter+1, eval, work, lwork, info)
      END IF

      CALL mp_bcast(alpha,0,mp_group)
      CALL mp_bcast(beta,0,mp_group)
      CALL mp_bcast(eval,0,mp_group)
      
      IF (ABS(MAXVAL(eval(1:i))-max_ev) < threshold .AND. ABS(MINVAL(eval(1:i))-min_ev) < threshold .AND. i>=MIN(3,n)) THEN
        max_ev = MAXVAL(eval(1:i))
        min_ev = MINVAL(eval(1:i))
        converged = .TRUE.
        EXIT
      ENDIF

      ! update the current best value of max/min ev
      max_ev = MAXVAL(eval(1:i))
      min_ev = MINVAL(eval(1:i))

      ! if we have lost orthonormality, or the initial guess is an eigenstate (for example A is unit matrix), quit
      IF (ABS(beta)<EPSILON(beta)*5*SQRT(REAL(n,KIND=dp))) EXIT

      IF (rank .EQ. 0) THEN
        CALL DSCAL(n, 1/beta, u, 1)        ! 
        CALL DCOPY(n, u, 1, vn, 1)         ! vi+1 = ui/beta_i
      END IF
      CALL mp_bcast(vn, 0, mp_group)
  
      u(:,:) = 0.0_dp
      IF(svd) THEN
        u_tmp(:,:) = 0.0_dp
        CALL cp_dbcsr_multiply_local(matrix_a, vn, u_tmp, error=error)  !  
        CALL mp_sum(u_tmp, mp_group)
        CALL cp_dbcsr_multiply_local(matrix_a, u_tmp, u, error=error)  !  
        CALL mp_sum(u, mp_group)
      ELSE
        CALL cp_dbcsr_multiply_local(matrix_a, vn, u, error=error)  !  
        CALL mp_sum(u, mp_group)
      ENDIF

      IF (rank .EQ. 0) THEN
        CALL DAXPY(n, -beta, v, 1, u, 1)                            ! ui+1 = A*vi+1 - beta_i*vi
        CALL DCOPY(n, vn, 1, v, 1)
      END IF
      CALL mp_bcast(u, 0, mp_group)
      CALL mp_bcast(v, 0, mp_group)
    ENDDO

    IF(svd) THEN
      max_ev = SQRT(max_ev)
      min_ev = SQRT(min_ev)
    ENDIF

    ! make absolutely sure we have the same values on all processes
    CALL mp_min(min_ev,mp_group)
    CALL mp_max(max_ev,mp_group)

    DEALLOCATE(v, u, vn, evs, evs_tmp, eval, work)
    IF(svd) DEALLOCATE(u_tmp)
  END SUBROUTINE lanczos_alg_serial_d


! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param max_ev ...
!> \param min_ev ...
!> \param converged ...
!> \param threshold ...
!> \param max_iter ...
!> \param svd ...
!> \param mp_group ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE lanczos_alg_serial_s(matrix_a, max_ev, min_ev, converged, threshold, max_iter, svd, mp_group, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(KIND=dp), INTENT(OUT)               :: max_ev, min_ev
    LOGICAL, INTENT(OUT)                     :: converged
    REAL(KIND=dp), INTENT(IN)                :: threshold
    INTEGER, INTENT(IN)                      :: max_iter
    LOGICAL, INTENT(IN)                      :: svd
    INTEGER                                  :: mp_group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_alg_serial_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, info, lwork, n, ntask, rank
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eval, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: evs, evs_tmp
    REAL(kind=sp)                            :: alpha, beta, nrm_v, SDOT, &
                                                SNRM2
    REAL(KIND=sp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: u, u_tmp, v, vn
    TYPE(rng_stream_type), POINTER           :: rng_stream

    converged = .FALSE.

    CALL mp_environ(ntask,rank,mp_group)
    CALL cp_dbcsr_get_info(matrix_a, nfullcols_total=n)
    ALLOCATE(v(n,1), u(n,1), vn(n,1))
    IF(svd) ALLOCATE(u_tmp(n,1))

    ALLOCATE(evs(max_iter+1, max_iter+1), evs_tmp(max_iter+1, max_iter+1))
    lwork = max_iter*3+100
    ALLOCATE(eval(max_iter), work(lwork))

    NULLIFY (rng_stream)
    CALL create_rng_stream(rng_stream,name="lanczos_alg_serial",distribution_type=UNIFORM,&
                           extended_precision=.TRUE.,error=error)
    DO i=1,n
      v(i,1) = REAL(next_random_number(rng_stream,error=error) - 0.5_dp,KIND=sp)
    END DO
    CALL delete_rng_stream(rng_stream,error=error)
    CALL mp_sum(v, mp_group)
    nrm_v = SDOT(n,v, 1, v, 1)
    v(:,:) = v(:,:)/SQRT(nrm_v)

    u(:,:) = 0.0
    IF(svd) THEN
      u_tmp(:,:) = 0.0
      CALL cp_dbcsr_multiply_local(matrix_a, v, u_tmp, error=error)
      CALL mp_sum(u_tmp, mp_group)
      CALL cp_dbcsr_multiply_local(matrix_a, u_tmp, u, error=error)
      CALL mp_sum(u, mp_group)
    ELSE
      CALL cp_dbcsr_multiply_local(matrix_a, v, u, error=error)
      CALL mp_sum(u, mp_group)
    ENDIF

    evs(:,:) = 0.0_dp
    max_ev = 0.0_dp
    min_ev = 0.0_dp

    DO i=1,max_iter
      !Perform computation in serial to avoid numerical divergence
      !of ill-conditioned operations and then broadcast results
      IF (rank .EQ. 0) THEN
        alpha = SDOT(n, u, 1, v, 1)        ! alpha_i = ui^T*vi
        CALL SAXPY(n, -alpha, v, 1, u, 1)  ! ui = ui-alpha*vi   
        beta = SNRM2(n, u, 1)              ! beta_i = sqrt(ri^T*ri)

        evs(i,i) = REAL(alpha, dp)
        evs(i,i+1) = REAL(beta, dp)
        evs(i+1,i) = REAL(beta, dp)
        evs_tmp(:,:) = evs(:,:)
        CALL DSYEV('N', 'U', i, evs_tmp, max_iter+1, eval, work, lwork, info)
      END IF

      CALL mp_bcast(alpha,0,mp_group)
      CALL mp_bcast(beta,0,mp_group)
      CALL mp_bcast(eval,0,mp_group)

      IF (ABS(MAXVAL(eval(1:i))-max_ev) < threshold .AND. ABS(MINVAL(eval(1:i))-min_ev) < threshold .AND. i>=MIN(3,n)) THEN
        max_ev = MAXVAL(eval(1:i))
        min_ev = MINVAL(eval(1:i))
        converged = .TRUE.
        EXIT
      ENDIF

      ! update the current best value of max/min ev
      max_ev = MAXVAL(eval(1:i))
      min_ev = MINVAL(eval(1:i))

      ! if we have lost orthonormality, or the initial guess is an eigenstate (for example A is unit matrix), quit
      IF (ABS(beta)<EPSILON(beta)*5*SQRT(REAL(n,KIND=dp))) EXIT

      IF (rank .EQ. 0) THEN
        CALL SSCAL(n, 1/beta, u, 1)        ! 
        CALL SCOPY(n, u, 1, vn, 1)         ! vi+1 = ui/beta_i
      END IF
      CALL mp_bcast(vn, 0, mp_group)
  
      u(:,:) = 0.0
      IF(svd) THEN
        u_tmp(:,:) = 0.0
        CALL cp_dbcsr_multiply_local(matrix_a, vn, u_tmp, error=error)  !  
        CALL mp_sum(u_tmp, mp_group)
        CALL cp_dbcsr_multiply_local(matrix_a, u_tmp, u, error=error)  !  
        CALL mp_sum(u, mp_group)
      ELSE
        CALL cp_dbcsr_multiply_local(matrix_a, vn, u, error=error)  !  
        CALL mp_sum(u, mp_group)
      ENDIF

      IF (rank .EQ. 0) THEN
        CALL SAXPY(n, -beta, v, 1, u, 1)                            ! ui+1 = A*vi+1 - beta_i*vi
        CALL SCOPY(n, vn, 1, v, 1)
      END IF
      CALL mp_bcast(u, 0, mp_group)
      CALL mp_bcast(v, 0, mp_group)
    ENDDO

    IF(svd) THEN
      max_ev = SQRT(max_ev)
      min_ev = SQRT(min_ev)
    ENDIF

    ! make absolutely sure we have the same values on all processes
    CALL mp_min(min_ev,mp_group)
    CALL mp_max(max_ev,mp_group)

    DEALLOCATE(v, u, vn, evs, evs_tmp, eval, work)
    IF(svd) DEALLOCATE(u_tmp)
  END SUBROUTINE lanczos_alg_serial_s
END MODULE cp_dbcsr_util
