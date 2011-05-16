!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Routines useful for iterative matrix calculations
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE iterate_matrix
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_filter, cp_dbcsr_frobenius_norm, cp_dbcsr_gershgorin_norm, &
       cp_dbcsr_get_occupation, cp_dbcsr_init, cp_dbcsr_multiply, &
       cp_dbcsr_release, cp_dbcsr_scale, cp_dbcsr_set
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'iterate_matrix'

  PUBLIC ::  invert_Hotelling, matrix_sign_Newton_Schulz, matrix_sqrt_Newton_Schulz

CONTAINS

! *****************************************************************************
!> \brief invert a symmetric positive definite matrix by Hotelling's method 
!>        explicit symmetrization makes this code not suitable for other matrix types
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE invert_Hotelling(matrix_inverse,matrix,threshold,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_inverse, matrix
    REAL(KIND=dp), INTENT(IN)                :: threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'invert_Hotelling', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, unit_nr
    INTEGER(KIND=int_8)                      :: flop1, flop2
    REAL(KIND=dp)                            :: frob_matrix, &
                                                frob_matrix_base, &
                                                gershgorin_norm, occ_matrix, &
                                                t1, t2
    TYPE(cp_dbcsr_type)                      :: tmp1, tmp2
    TYPE(cp_logger_type), POINTER            :: logger

! turn this off for the time being

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    gershgorin_norm=cp_dbcsr_gershgorin_norm(matrix)
    frob_matrix=cp_dbcsr_frobenius_norm(matrix)
    CALL cp_dbcsr_set(matrix_inverse,0.0_dp,error=error)
    CALL cp_dbcsr_add_on_diag(matrix_inverse,1/MIN(gershgorin_norm,frob_matrix),error=error)

    CALL cp_dbcsr_init(tmp1,error=error)
    CALL cp_dbcsr_create(tmp1,template=matrix_inverse,error=error)
    CALL cp_dbcsr_init(tmp2,error=error)
    CALL cp_dbcsr_create(tmp2,template=matrix_inverse,error=error)

    IF (unit_nr>0) WRITE(unit_nr,*)

    DO i=1,100

       t1 = m_walltime()
       ! tmp1 = S^-1 S
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_inverse, matrix,&
                              0.0_dp, tmp1,flop=flop1, error=error)

       ! for the convergence check
       frob_matrix_base=cp_dbcsr_frobenius_norm(tmp1)
       CALL cp_dbcsr_add_on_diag(tmp1,-1.0_dp,error=error)
       frob_matrix=cp_dbcsr_frobenius_norm(tmp1)
       CALL cp_dbcsr_add_on_diag(tmp1,+1.0_dp,error=error)
       occ_matrix=cp_dbcsr_get_occupation(matrix_inverse)

       ! tmp2 = S^-1 S S^-1
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, tmp1, matrix_inverse, 0.0_dp, tmp2,&
                              flop=flop2, error=error)
       ! S^-1_{n+1} = 2 S^-1 - S^-1 S S^-1
       CALL cp_dbcsr_add(matrix_inverse, tmp2, 2.0_dp, -1.0_dp, error=error)

       CALL cp_dbcsr_filter(matrix_inverse, threshold, error=error) 
       t2 = m_walltime()

       IF (unit_nr>0) THEN
          WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)') "Hotelling iter",i,occ_matrix, &
                                                             frob_matrix/frob_matrix_base,t2-t1,&
                                                                    (flop1+flop2)/(1.0E6_dp*(t2-t1))
          CALL m_flush(unit_nr)
       ENDIF
 
       ! convergence.... convergence is quadratic, so the current iteration made is below the threshold
       IF (frob_matrix/frob_matrix_base<SQRT(threshold)) EXIT

    ENDDO

    ! this check is not really needed
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_inverse, matrix, 0.0_dp, tmp1,error=error)
    frob_matrix_base=cp_dbcsr_frobenius_norm(tmp1)
    CALL cp_dbcsr_add_on_diag(tmp1,-1.0_dp,error=error)
    frob_matrix=cp_dbcsr_frobenius_norm(tmp1)
    occ_matrix=cp_dbcsr_get_occupation(matrix_inverse)
    IF (unit_nr>0) THEN
       WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.8,E12.3)') "Final Hotelling ",i,occ_matrix,frob_matrix/frob_matrix_base

       WRITE(unit_nr,'()') 
       CALL m_flush(unit_nr)
    ENDIF

    CALL cp_dbcsr_release(tmp1,error=error) 
    CALL cp_dbcsr_release(tmp2,error=error) 

    CALL timestop(handle)

  END SUBROUTINE invert_Hotelling

! *****************************************************************************
!> \brief compute the sign a matrix using Newton-Schulz iterations
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE matrix_sign_Newton_Schulz(matrix_sign,matrix,threshold,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_sign, matrix
    REAL(KIND=dp), INTENT(IN)                :: threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_sign_Newton_Schulz', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: count, handle, i, unit_nr
    INTEGER(KIND=int_8)                      :: flop1, flop2
    REAL(KIND=dp)                            :: frob_matrix, &
                                                frob_matrix_base, &
                                                gersh_matrix, occ_matrix, t1, &
                                                t2
    TYPE(cp_dbcsr_type)                      :: tmp1, tmp2
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    CALL cp_dbcsr_init(tmp1,error=error)
    CALL cp_dbcsr_create(tmp1,template=matrix_sign,error=error)

    CALL cp_dbcsr_init(tmp2,error=error)
    CALL cp_dbcsr_create(tmp2,template=matrix_sign,error=error)

    CALL cp_dbcsr_copy(matrix_sign,matrix,error=error)
    CALL cp_dbcsr_filter(matrix_sign,threshold,error=error)

    ! scale the matrix to get into the convergence range
    frob_matrix=cp_dbcsr_frobenius_norm(matrix_sign)
    gersh_matrix=cp_dbcsr_gershgorin_norm(matrix_sign)
    CALL cp_dbcsr_scale(matrix_sign,1/MIN(frob_matrix,gersh_matrix),error=error)

    IF (unit_nr>0) WRITE(unit_nr,*)

    count=0
    DO i=1,100

       t1 = m_walltime()
       ! tmp1 = X * X
       CALL cp_dbcsr_multiply("N", "N", -1.0_dp, matrix_sign, matrix_sign, 0.0_dp, tmp1,&
                              filter_eps=threshold, flop=flop1, error=error)

       ! check convergence (frob norm of what should be the identity matrix minus identity matrix)
       frob_matrix_base=cp_dbcsr_frobenius_norm(tmp1)
       CALL cp_dbcsr_add_on_diag(tmp1,+1.0_dp,error=error)
       frob_matrix=cp_dbcsr_frobenius_norm(tmp1)
        
       ! update the above to 3*I-X*X
       CALL cp_dbcsr_add_on_diag(tmp1,+2.0_dp,error=error)
       occ_matrix=cp_dbcsr_get_occupation(matrix_sign)

       ! tmp2 = 0.5 * X * (3*I-X*X)
       CALL cp_dbcsr_multiply("N", "N", 0.5_dp, matrix_sign, tmp1, 0.0_dp, tmp2, &
                              filter_eps=threshold, flop=flop2, error=error)

       ! done iterating
       ! CALL cp_dbcsr_filter(tmp2,threshold,error=error)
       CALL cp_dbcsr_copy(matrix_sign,tmp2,error=error)
       t2 = m_walltime()

       IF (unit_nr>0) THEN
          WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)') "NS sign iter ",i,occ_matrix,&
                                                                    frob_matrix/frob_matrix_base,t2-t1,&
                                                                    (flop1+flop2)/(1.0E6_dp*(t2-t1))
          CALL m_flush(unit_nr)
       ENDIF

       IF (frob_matrix/frob_matrix_base<SQRT(threshold)) EXIT

    ENDDO

    ! this check is not really needed
    CALL cp_dbcsr_multiply("N", "N", +1.0_dp, matrix_sign, matrix_sign, 0.0_dp, tmp1,&
                           filter_eps=threshold, error=error)
    frob_matrix_base=cp_dbcsr_frobenius_norm(tmp1)
    CALL cp_dbcsr_add_on_diag(tmp1,-1.0_dp,error=error)
    frob_matrix=cp_dbcsr_frobenius_norm(tmp1)
    occ_matrix=cp_dbcsr_get_occupation(matrix_sign)
    IF (unit_nr>0) THEN
       WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.8,E12.3)') "Final NS sign iter",i,occ_matrix,&
                                                    frob_matrix/frob_matrix_base
       WRITE(unit_nr,'()') 
       CALL m_flush(unit_nr)
    ENDIF

    CALL cp_dbcsr_release(tmp1,error=error) 
    CALL cp_dbcsr_release(tmp2,error=error) 

    CALL timestop(handle)

  END SUBROUTINE matrix_sign_Newton_Schulz

! *****************************************************************************
!> \brief compute the sqrt of a matrix via the sign function and the corresponding Newton-Schulz iterations
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE matrix_sqrt_Newton_Schulz(matrix_sqrt,matrix_sqrt_inv,matrix,threshold,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_sqrt, matrix_sqrt_inv, &
                                                matrix
    REAL(KIND=dp), INTENT(IN)                :: threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_sqrt_Newton_Schulz', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, unit_nr
    INTEGER(KIND=int_8)                      :: flop1, flop2, flop3
    REAL(KIND=dp)                            :: frob_matrix, &
                                                frob_matrix_base, &
                                                gersh_matrix, occ_matrix, &
                                                scaling, t1, t2
    TYPE(cp_dbcsr_type)                      :: tmp1, tmp2
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    CALL cp_dbcsr_init(tmp1,error=error)
    CALL cp_dbcsr_create(tmp1,template=matrix,error=error)
    CALL cp_dbcsr_init(tmp2,error=error)
    CALL cp_dbcsr_create(tmp2,template=matrix,error=error)

    CALL cp_dbcsr_set(matrix_sqrt_inv,0.0_dp,error=error)
    CALL cp_dbcsr_add_on_diag(matrix_sqrt_inv,1.0_dp,error=error)
    CALL cp_dbcsr_filter(matrix_sqrt_inv,threshold,error=error)

    ! scale the matrix to get into the convergence range
    CALL cp_dbcsr_copy(matrix_sqrt,matrix,error=error)
    frob_matrix=cp_dbcsr_frobenius_norm(matrix_sqrt)
    gersh_matrix=cp_dbcsr_gershgorin_norm(matrix_sqrt)
    scaling=1/MIN(frob_matrix,gersh_matrix)
    CALL cp_dbcsr_scale(matrix_sqrt,scaling,error=error)
    CALL cp_dbcsr_filter(matrix_sqrt,threshold,error=error)

    IF (unit_nr>0) WRITE(unit_nr,*)

    DO i=1,100

       t1 = m_walltime()

       ! tmp1 = I - Zk * Yk
       CALL cp_dbcsr_multiply("N", "N", -1.0_dp, matrix_sqrt_inv, matrix_sqrt, 0.0_dp, tmp1,&
                              filter_eps=threshold, flop=flop1, error=error)
       frob_matrix_base=cp_dbcsr_frobenius_norm(tmp1)
       CALL cp_dbcsr_add_on_diag(tmp1,+1.0_dp,error=error)
 
       ! check convergence (frob norm of what should be the identity matrix minus identity matrix)
       frob_matrix=cp_dbcsr_frobenius_norm(tmp1)

       ! update the above to 3*I-Zk*Yk
       CALL cp_dbcsr_add_on_diag(tmp1,+2.0_dp,error=error)

       ! tmp2 = 0.5 * Yk * tmp1 = Y(k+1)
       CALL cp_dbcsr_multiply("N", "N",  0.5_dp, matrix_sqrt, tmp1,  0.0_dp, tmp2,&
                              filter_eps=threshold, flop=flop2, error=error)
       ! CALL cp_dbcsr_filter(tmp2,threshold,error=error)
       CALL cp_dbcsr_copy(matrix_sqrt, tmp2, error=error)

       ! tmp2 = 0.5 * tmp1 * Zk = Z(k+1)
       CALL cp_dbcsr_multiply("N", "N",  0.5_dp, tmp1, matrix_sqrt_inv,  0.0_dp, tmp2,&
                              filter_eps=threshold,flop=flop3, error=error)
       ! CALL cp_dbcsr_filter(tmp2,threshold,error=error)
       CALL cp_dbcsr_copy(matrix_sqrt_inv, tmp2, error=error)

       occ_matrix=cp_dbcsr_get_occupation(matrix_sqrt_inv)

       ! done iterating
       t2 = m_walltime()

       IF (unit_nr>0) THEN
          WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)') "NS sqrt iter ",i,occ_matrix,&
                                                             frob_matrix/frob_matrix_base,t2-t1,&
                                                             (flop1+flop2+flop3)/(1.0E6_dp*(t2-t1))
          CALL m_flush(unit_nr)
       ENDIF

       IF (frob_matrix/frob_matrix_base<SQRT(threshold)) EXIT

    ENDDO

    ! this check is not really needed
    CALL cp_dbcsr_multiply("N", "N", +1.0_dp, matrix_sqrt_inv, matrix_sqrt, 0.0_dp, tmp1,&
                           filter_eps=threshold,error=error)
    frob_matrix_base=cp_dbcsr_frobenius_norm(tmp1)
    CALL cp_dbcsr_add_on_diag(tmp1,-1.0_dp,error=error)
    frob_matrix=cp_dbcsr_frobenius_norm(tmp1)
    occ_matrix=cp_dbcsr_get_occupation(matrix_sqrt_inv)
    IF (unit_nr>0) THEN
       WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.8,E12.3)') "Final NS sqrt iter ",i,occ_matrix,&
                                                    frob_matrix/frob_matrix_base
       WRITE(unit_nr,'()') 
       CALL m_flush(unit_nr)
    ENDIF

    ! scale to proper end results
    CALL cp_dbcsr_scale(matrix_sqrt,1/SQRT(scaling),error=error)
    CALL cp_dbcsr_scale(matrix_sqrt_inv,SQRT(scaling),error=error)

    CALL cp_dbcsr_release(tmp1,error=error) 
    CALL cp_dbcsr_release(tmp2,error=error) 

    CALL timestop(handle)

  END SUBROUTINE matrix_sqrt_Newton_Schulz

END MODULE iterate_matrix
