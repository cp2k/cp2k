!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief helper routines for cp_dbcsr
!> \par History
!>       2012.07 created [Florian Thoele]
!> \author Florian Thoele
! *****************************************************************************

MODULE cp_dbcsr_util
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_info
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_multiply_local
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE dbcsr_types,                     ONLY: dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_max,&
                                             mp_min,&
                                             mp_sum
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_util'

  PUBLIC lanczos_alg_serial

CONTAINS

! *****************************************************************************
!> \brief Computes the highest and lowest eigenvalue of matrix_a (assumed to be symmetric)
!>        with the Lanczos algorithm. The 3 Lanczos vectors are replicated.
!>        speedups could be gained from an improved dbcsr x vector multiplication (distributed vector).
!> \par History
!>       2012.07 created [Florian Thoele]
!> \author Florian Thoele 
! *****************************************************************************
  SUBROUTINE lanczos_alg_serial(matrix_a, max_ev, min_ev, threshold, max_iter, svd, converged, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(KIND=dp), INTENT(OUT)               :: max_ev, min_ev
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: threshold
    INTEGER, INTENT(IN), OPTIONAL            :: max_iter
    LOGICAL, INTENT(IN), OPTIONAL            :: svd
    LOGICAL, INTENT(OUT), OPTIONAL           :: converged
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_alg_serial', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: default_max_iter = 40
    REAL, PARAMETER                          :: default_threshold = 1.0E-6_dp

    INTEGER                                  :: handle, mp_group, &
                                                my_max_iter, unit_nr
    LOGICAL                                  :: do_svd
    REAL(KIND=dp)                            :: my_threshold
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    mp_group = matrix_a%matrix%m%dist%d%mp_env%mp%mp_group

    IF(PRESENT(converged)) converged = .FALSE.
    my_threshold = default_threshold
    IF(PRESENT(threshold)) my_threshold = threshold
    my_max_iter = default_max_iter
    IF(PRESENT(max_iter)) my_max_iter = max_iter
    do_svd = .FALSE.
    IF(PRESENT(svd)) do_svd = svd

    SELECT CASE(matrix_a%matrix%m%data_type)
      CASE (dbcsr_type_real_8)
        CALL lanczos_alg_serial_d(matrix_a, max_ev, min_ev, my_threshold, my_max_iter, do_svd, converged, mp_group, error)
      CASE (dbcsr_type_real_4)
        CALL lanczos_alg_serial_s(matrix_a, max_ev, min_ev, my_threshold, my_max_iter, do_svd, converged, mp_group, error)
    !  CASE DEFAULT
    !    CALL stop_program(routineN,moduleN,__LINE__,"NYI")
    END SELECT

    CALL timestop(handle)
  END SUBROUTINE lanczos_alg_serial

  SUBROUTINE lanczos_alg_serial_d(matrix_a, max_ev, min_ev, threshold, max_iter, svd, converged, mp_group, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(KIND=dp), INTENT(OUT)               :: max_ev, min_ev
    REAL(KIND=dp), INTENT(IN)                :: threshold
    INTEGER, INTENT(IN)                      :: max_iter
    LOGICAL, INTENT(IN)                      :: svd
    LOGICAL, INTENT(OUT), OPTIONAL           :: converged
    INTEGER                                  :: mp_group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_alg_serial_d', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: default_max_iter = 40
    REAL(KIND=dp), PARAMETER                 :: default_threshold = 1.0E-6_dp

    INTEGER                                  :: i, info, lwork, n
    REAL(KIND=dp)                            :: alpha, beta, DDOT, DNRM2, &
                                                nrm_v
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eval, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: evs, evs_tmp, u, u_tmp, v, vn
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CALL cp_dbcsr_get_info(matrix_a, nfullcols_total=n)
    ALLOCATE(v(n,1), u(n,1), vn(n,1))
    IF(svd) ALLOCATE(u_tmp(n,1))

    ALLOCATE(evs(max_iter+1, max_iter+1), evs_tmp(max_iter+1, max_iter+1))
    lwork = max_iter*3+100
    ALLOCATE(eval(max_iter), work(lwork))


!MK    CALL RANDOM_SEED()
!MK    DO i=1,n 
!MK      CALL RANDOM_NUMBER(v(i,1)); v(i,1 ) = v(i,1) - 0.5_dp
!MK    ENDDO
!MK Use the internal CP2K random number generator
    NULLIFY (rng_stream)
    CALL create_rng_stream(rng_stream,name="lanczos_alg_serial",distribution_type=UNIFORM,&
                           extended_precision=.TRUE.,error=error)
    DO i=1,n
      v(i,1) = next_random_number(rng_stream,error=error) - 0.5_dp
    END DO
    CALL delete_rng_stream(rng_stream,error=error)
    CALL mp_sum(v, mp_group)
    nrm_v = DDOT(n,v, 1, v, 1)
    v = v/SQRT(nrm_v)

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
      alpha = DDOT(n, u, 1, v, 1)        ! alpha_i = ui^T*vi
      CALL DAXPY(n, -alpha, v, 1, u, 1)  ! ui = ui-alpha*vi   
      beta = DNRM2(n, u, 1)              ! beta_i = sqrt(ri^T*ri)

      evs(i,i) = alpha
      evs(i,i+1) = beta
      evs(i+1,i) = beta
      evs_tmp(:,:) = evs(:,:)
      CALL DSYEV('N', 'U', i, evs_tmp, max_iter+1, eval, work, lwork, info)

      IF (ABS(MAXVAL(eval(1:i))-max_ev) < threshold .AND. ABS(MINVAL(eval(1:i))-min_ev) < threshold) THEN
        max_ev = MAXVAL(eval(1:i))
        min_ev = MINVAL(eval(1:i))
        IF(PRESENT(converged)) converged = .TRUE.
        EXIT
      ENDIF

      ! update the current best value of max/min ev
      max_ev = MAXVAL(eval(1:i))
      min_ev = MINVAL(eval(1:i))

      ! if we have lost orthonormality, or the initial guess is an eigenstate (for example A is unit matrix), quit
      IF (ABS(beta)<EPSILON(beta)*5*SQRT(REAL(n,KIND=dp))) EXIT

      CALL DSCAL(n, 1/beta, u, 1)        ! 
      CALL DCOPY(n, u, 1, vn, 1)         ! vi+1 = ui/beta_i
  
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

      CALL DAXPY(n, -beta, v, 1, u, 1)                            ! ui+1 = A*vi+1 - beta_i*vi
      CALL DCOPY(n, vn, 1, v, 1)
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


  SUBROUTINE lanczos_alg_serial_s(matrix_a, max_ev, min_ev, threshold, max_iter, svd, converged, mp_group, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(KIND=dp), INTENT(OUT)               :: max_ev, min_ev
    REAL(KIND=dp), INTENT(IN)                :: threshold
    INTEGER, INTENT(IN)                      :: max_iter
    LOGICAL, INTENT(IN)                      :: svd
    LOGICAL, INTENT(OUT), OPTIONAL           :: converged
    INTEGER                                  :: mp_group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_alg_serial_s', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: default_max_iter = 40
    REAL(KIND=dp), PARAMETER                 :: default_threshold = 1.0E-6_dp

    INTEGER                                  :: i, info, lwork, n
    REAL                                     :: alpha, beta, nrm_v, SDOT, &
                                                SNRM2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eval, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: evs, evs_tmp
    REAL, ALLOCATABLE, DIMENSION(:, :)       :: u, u_tmp, v, vn
    TYPE(rng_stream_type), POINTER           :: rng_stream

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
    nrm_v = SDOT(n,v, 1, v, 1)
    v = v/SQRT(nrm_v)

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
      alpha = SDOT(n, u, 1, v, 1)        ! alpha_i = ui^T*vi
      CALL SAXPY(n, -alpha, v, 1, u, 1)  ! ui = ui-alpha*vi   
      beta = SNRM2(n, u, 1)              ! beta_i = sqrt(ri^T*ri)

      evs(i,i) = REAL(alpha, dp)
      evs(i,i+1) = REAL(beta, dp)
      evs(i+1,i) = REAL(beta, dp)
      evs_tmp(:,:) = evs(:,:)
      CALL DSYEV('N', 'U', i, evs_tmp, max_iter+1, eval, work, lwork, info)

      IF (ABS(MAXVAL(eval(1:i))-max_ev) < threshold .AND. ABS(MINVAL(eval(1:i))-min_ev) < threshold) THEN
        max_ev = MAXVAL(eval(1:i))
        min_ev = MINVAL(eval(1:i))
        IF(PRESENT(converged)) converged = .TRUE.
        EXIT
      ENDIF

      ! update the current best value of max/min ev
      max_ev = MAXVAL(eval(1:i))
      min_ev = MINVAL(eval(1:i))

      ! if we have lost orthonormality, or the initial guess is an eigenstate (for example A is unit matrix), quit
      IF (ABS(beta)<EPSILON(beta)*5*SQRT(REAL(n,KIND=dp))) EXIT

      CALL SSCAL(n, 1/beta, u, 1)        ! 
      CALL SCOPY(n, u, 1, vn, 1)         ! vi+1 = ui/beta_i
  
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

      CALL SAXPY(n, -beta, v, 1, u, 1)                            ! ui+1 = A*vi+1 - beta_i*vi
      CALL SCOPY(n, vn, 1, v, 1)
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
