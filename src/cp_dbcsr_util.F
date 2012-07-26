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
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_util'

  PUBLIC lanczos_alg_serial

CONTAINS

! *****************************************************************************
!> \brief Computes the highest and lowest eigenvalue of matrix_a 
!>        with the Lanczos algorithm. The Lanczos vectors are replicated.
!> \par History
!>       2012.07 created [Florian Thoele]
!> \author Florian Thoele 
! *****************************************************************************
  SUBROUTINE lanczos_alg_serial(matrix_a, max_ev, min_ev, threshold, max_iter, converged, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(KIND=dp), INTENT(OUT)               :: max_ev, min_ev
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: threshold
    INTEGER, INTENT(IN), OPTIONAL            :: max_iter
    LOGICAL, INTENT(OUT), OPTIONAL           :: converged
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_alg_serial', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: default_max_iter = 40
    REAL, PARAMETER                          :: default_threshold = 1.0E-6_dp

    INTEGER                                  :: handle, i, info, lwork, &
                                                mp_group, my_max_iter, n, &
                                                unit_nr
    REAL(KIND=dp)                            :: alpha, beta, DDOT, DNRM2, &
                                                my_threshold, nrm_v
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eval, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: evs, evs_tmp, u, v, vn
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

    CALL cp_dbcsr_get_info(matrix_a, nfullcols_total=n)
    ALLOCATE(v(n,1), u(n,1), vn(n,1))
    ALLOCATE(evs(my_max_iter+1, my_max_iter+1), evs_tmp(my_max_iter+1, my_max_iter+1))
    lwork = my_max_iter*3+100
    ALLOCATE(eval(my_max_iter), work(lwork))

    CALL RANDOM_SEED()
    DO i=1,n 
      CALL RANDOM_NUMBER(v(i,1)); v(i,1 ) = v(i,1) - 0.5_dp
    ENDDO
    nrm_v = DDOT(n,v, 1, v, 1)
    v = v/SQRT(nrm_v)

    u(:,:) = 0.0_dp
    CALL cp_dbcsr_multiply_local(matrix_a, v, u, error=error)
    CALL mp_sum(u, mp_group)

    evs(:,:) = 0.0_dp
    max_ev = 0.0_dp
    min_ev = 0.0_dp

    DO i=1,my_max_iter
      alpha = DDOT(n, u, 1, v, 1)        ! alpha_i = ui^T*vi
      CALL DAXPY(n, -alpha, v, 1, u, 1)  ! ui = ui-alpha*vi   
      beta = DNRM2(n, u, 1)              ! beta_i = sqrt(ri^T*ri)
      CALL DSCAL(n, 1/beta, u, 1)        ! 
      CALL DCOPY(n, u, 1, vn, 1)         ! vi+1 = ui/beta_i
  
      u(:,:) = 0.0_dp
      CALL cp_dbcsr_multiply_local(matrix_a, vn, u, error=error)  !  
      CALL mp_sum(u, mp_group)
      CALL DAXPY(n, -beta, v, 1, u, 1)                            ! ui+1 = A*vi+1 - beta_i*vi
  
      evs(i,i) = alpha
      evs(i,i+1) = beta
      evs(i+1,i) = beta
      evs_tmp(:,:) = evs(:,:)
      CALL DSYEV('N', 'U', i, evs_tmp, my_max_iter+1, eval, work, lwork, info)

      IF (ABS(MAXVAL(eval(1:i))-max_ev) < my_threshold .AND. ABS(MINVAL(eval(1:i))-min_ev) < my_threshold) THEN
        max_ev = MAXVAL(eval(1:i))
        min_ev = MINVAL(eval(1:i))
        IF(PRESENT(converged)) converged = .TRUE.
        EXIT
      ENDIF

      max_ev = MAXVAL(eval(1:i))
      min_ev = MINVAL(eval(1:i))

      CALL DCOPY(n, vn, 1, v, 1)
    ENDDO

    DEALLOCATE(v, u, vn, evs, evs_tmp, eval, work)
    CALL timestop(handle)
  END SUBROUTINE lanczos_alg_serial
END MODULE cp_dbcsr_util
