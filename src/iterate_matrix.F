!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Routines useful for iterative matrix calculations
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE iterate_matrix
!  USE cp_dbcsr_arnoldi,                ONLY: cp_dbcsr_eig_arnoldi
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_arnoldi_ev, &
       cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_filter, &
       cp_dbcsr_frobenius_norm, cp_dbcsr_gershgorin_norm, cp_dbcsr_get_info, &
       cp_dbcsr_get_occupation, cp_dbcsr_init, cp_dbcsr_multiply, &
       cp_dbcsr_norm, cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_scale, &
       cp_dbcsr_set, cp_dbcsr_setup_arnoldi_data, cp_dbcsr_transposed, &
       cp_dbcsr_type, dbcsr_arnoldi_data, dbcsr_norm_maxabsnorm, &
       dbcsr_type_no_symmetry, deallocate_arnoldi_data, get_selected_ritz_val
  USE cp_dbcsr_util,                   ONLY: lanczos_alg_serial
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE mathconstants,                   ONLY: ifac
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'iterate_matrix'

  INTERFACE purify_mcweeny
     MODULE PROCEDURE purify_mcweeny_orth, purify_mcweeny_nonorth
  END INTERFACE

  PUBLIC :: invert_Hotelling, matrix_sign_Newton_Schulz, matrix_sqrt_Newton_Schulz,&
            purify_mcweeny

CONTAINS

! *****************************************************************************
!> \brief invert a symmetric positive definite matrix by Hotelling's method
!>        explicit symmetrization makes this code not suitable for other matrix types
!> \param matrix_inverse ...
!> \param matrix ...
!> \param threshold ...
!> \param use_inv_as_guess ...
!> \param norm_convergence ...
!> \param error ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!>       2011.10 guess option added [Rustam Z Khaliullin]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE invert_Hotelling(matrix_inverse,matrix,threshold,use_inv_as_guess,&
                              norm_convergence,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT), &
      TARGET                                 :: matrix_inverse, matrix
    REAL(KIND=dp), INTENT(IN)                :: threshold
    LOGICAL, INTENT(IN), OPTIONAL            :: use_inv_as_guess
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: norm_convergence
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'invert_Hotelling', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nrow, unit_nr
    INTEGER(KIND=int_8)                      :: flop1, flop2
    LOGICAL                                  :: prepare_to_exit, use_inv_guess
    REAL(KIND=dp) :: convergence, frob_matrix, gershgorin_norm, max_abs_eig, &
      max_ev, maxnorm_matrix, min_ev, occ_matrix, t1, t2
    TYPE(cp_dbcsr_p_type), DIMENSION(2)      :: mymat
    TYPE(cp_dbcsr_type)                      :: tmp1, tmp2
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_arnoldi_data)                 :: my_arnoldi

! turn this off for the time being

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    convergence=threshold
    IF(PRESENT(norm_convergence))convergence=norm_convergence

    use_inv_guess=.FALSE.
    IF (PRESENT(use_inv_as_guess)) use_inv_guess=use_inv_as_guess

    IF (.NOT.use_inv_guess) THEN
       gershgorin_norm=cp_dbcsr_gershgorin_norm(matrix)
       frob_matrix=cp_dbcsr_frobenius_norm(matrix)
       CALL cp_dbcsr_set(matrix_inverse,0.0_dp,error=error)
       CALL cp_dbcsr_add_on_diag(matrix_inverse,1/MIN(gershgorin_norm,frob_matrix),error=error)
    ENDIF
    CALL cp_dbcsr_get_info(matrix, nfullrows_total=nrow )

    CALL cp_dbcsr_init(tmp1,error=error)
    CALL cp_dbcsr_create(tmp1,template=matrix_inverse,error=error)
    CALL cp_dbcsr_init(tmp2,error=error)
    CALL cp_dbcsr_create(tmp2,template=matrix_inverse,error=error)

    IF (unit_nr>0) WRITE(unit_nr,*)

    prepare_to_exit=.FALSE.

    ! scale the approximate inverse to be within the convergence radius
    IF(use_inv_guess)THEN
       mymat(1)%matrix=>matrix
       mymat(2)%matrix=>matrix_inverse
       CALL cp_dbcsr_setup_arnoldi_data(my_arnoldi,mymat,max_iter=30,threshold=1.0E-4_dp,&
                                        selection_crit=1,nval_request=2)
       CALL cp_dbcsr_arnoldi_ev(mymat,my_arnoldi,error)
       max_eV=REAL(get_selected_ritz_val(my_arnoldi,2),dp)
       min_eV=REAL(get_selected_ritz_val(my_arnoldi,1),dp)
       max_abs_eig=MAX(max_eV,ABS(min_eV))
       CALL deallocate_arnoldi_data(my_arnoldi)
       CALL cp_dbcsr_scale(matrix_inverse,1.0_dp/max_abs_eig,error=error)
    END IF

    DO i=1,100

       t1 = m_walltime()
       ! tmp1 = S^-1 S
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_inverse, matrix,&
                              0.0_dp, tmp1,flop=flop1, error=error)

       ! for the convergence check
       !frob_matrix_base=cp_dbcsr_frobenius_norm(tmp1)
       CALL cp_dbcsr_add_on_diag(tmp1,-1.0_dp,error=error)
       frob_matrix=cp_dbcsr_frobenius_norm(tmp1)

       CALL cp_dbcsr_norm(tmp1,&
               dbcsr_norm_maxabsnorm, norm_scalar=maxnorm_matrix, error=error)

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
          WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)')  "Hotelling iter",i,occ_matrix, &
                                                                     maxnorm_matrix,t2-t1,&
                                                                    (flop1+flop2)/(1.0E6_dp*MAX(0.001_dp,t2-t1))
          CALL m_flush(unit_nr)
       ENDIF

       !use the eigenvalues of the matrix to check conevrgence
       IF(PRESENT(norm_convergence))THEN
          mymat(1)%matrix=>matrix
          mymat(2)%matrix=>matrix_inverse
          CALL cp_dbcsr_setup_arnoldi_data(my_arnoldi,mymat,max_iter=30,threshold=1.0E-4_dp,&
                                           selection_crit=1,nval_request=2)
          CALL cp_dbcsr_arnoldi_ev(mymat,my_arnoldi,error)
          max_eV=REAL(get_selected_ritz_val(my_arnoldi,2),dp)
          min_eV=REAL(get_selected_ritz_val(my_arnoldi,1),dp)
          max_abs_eig=MAX(max_eV-1.0_dp,ABS(min_eV-1.0_dp))
          CALL deallocate_arnoldi_data(my_arnoldi)
          IF( ABS(max_abs_eig) .LT. norm_convergence)EXIT
       ELSE
!       ! use max_norm 
          IF (maxnorm_matrix<convergence) EXIT 
       END IF

    ENDDO
!    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_inverse, matrix, 0.0_dp, tmp1,error=error)
    !!frob_matrix_base=cp_dbcsr_frobenius_norm(tmp1)
    !CALL cp_dbcsr_add_on_diag(tmp1,-1.0_dp,error=error)
    !frob_matrix=cp_dbcsr_frobenius_norm(tmp1)
    !CALL cp_dbcsr_norm(tmp1,&
    !        dbcsr_norm_maxabsnorm, norm_scalar=maxnorm_matrix, error=error)
    !occ_matrix=cp_dbcsr_get_occupation(matrix_inverse)
    IF (unit_nr>0) THEN
    !   WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.8,E12.3)') "Final Hotelling ",i,occ_matrix,&
    !      !frob_matrix/frob_matrix_base
    !      maxnorm_matrix
       WRITE(unit_nr,'()')
       CALL m_flush(unit_nr)
    ENDIF

    CALL cp_dbcsr_release(tmp1,error=error)
    CALL cp_dbcsr_release(tmp2,error=error)

    CALL timestop(handle)

  END SUBROUTINE invert_Hotelling

! *****************************************************************************
!> \brief compute the sign a matrix using Newton-Schulz iterations
!> \param matrix_sign ...
!> \param matrix ...
!> \param threshold ...
!> \param error ...
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
                                                                    (flop1+flop2)/(1.0E6_dp*MAX(0.001_dp,t2-t1))
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
!>        the order of the algorithm should be 2..5, 3 or 5 is recommended
!> \param matrix_sqrt ...
!> \param matrix_sqrt_inv ...
!> \param matrix ...
!> \param threshold ...
!> \param order ...
!> \param eps_lanczos ...
!> \param max_iter_lanczos ...
!> \param error ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE matrix_sqrt_Newton_Schulz(matrix_sqrt,matrix_sqrt_inv,matrix,threshold, order, eps_lanczos, max_iter_lanczos, error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_sqrt, matrix_sqrt_inv, &
                                                matrix
    REAL(KIND=dp), INTENT(IN)                :: threshold
    INTEGER, INTENT(IN)                      :: order
    REAL(KIND=dp), INTENT(IN)                :: eps_lanczos
    INTEGER, INTENT(IN)                      :: max_iter_lanczos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_sqrt_Newton_Schulz', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, unit_nr
    INTEGER(KIND=int_8)                      :: flop1, flop2, flop3, flop4, &
                                                flop5
    LOGICAL                                  :: converged, failure
    REAL(KIND=dp)                            :: a, b, c, d, frob_matrix, &
                                                frob_matrix_base, max_ev, &
                                                min_ev, oa, ob, oc, &
                                                occ_matrix, od, scaling, t1, &
                                                t2
    TYPE(cp_dbcsr_type)                      :: tmp1, tmp2, tmp3
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure=.FALSE.

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! for stability symmetry can not be assumed
    CALL cp_dbcsr_init(tmp1,error=error)
    CALL cp_dbcsr_create(tmp1,template=matrix,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(tmp2,error=error)
    CALL cp_dbcsr_create(tmp2,template=matrix,matrix_type=dbcsr_type_no_symmetry,error=error)
    IF (order.GE.4) THEN
       CALL cp_dbcsr_init(tmp3,error=error)
       CALL cp_dbcsr_create(tmp3,template=matrix,matrix_type=dbcsr_type_no_symmetry,error=error)
    ENDIF

    CALL cp_dbcsr_set(matrix_sqrt_inv,0.0_dp,error=error)
    CALL cp_dbcsr_add_on_diag(matrix_sqrt_inv,1.0_dp,error=error)
    CALL cp_dbcsr_filter(matrix_sqrt_inv,threshold,error=error)
    CALL cp_dbcsr_copy(matrix_sqrt,matrix,error=error)

    ! scale the matrix to get into the convergence range
    CALL lanczos_alg_serial(matrix_sqrt, max_ev, min_ev, threshold=eps_lanczos, &
                            max_iter=max_iter_lanczos, converged=converged, error=error)
    IF (unit_nr>0) THEN
       WRITE(unit_nr,*)
       WRITE(unit_nr,'(T6,A,1X,L1,A,E12.3)') "Lanczos converged: ",converged," threshold:",eps_lanczos
       WRITE(unit_nr,'(T6,A,1X,E12.3,E12.3)') "Est. extremal eigenvalues:",max_ev,min_ev
       WRITE(unit_nr,'(T6,A,1X,E12.3)') "Est. condition number :", max_ev/MAX(min_ev,EPSILON(min_ev))
    ENDIF
    ! conservatively assume we get a relatively large error (100*threshold_lanczos) in the estimates
    ! and adjust the scaling to be on the safe side
    scaling=2/(max_ev+min_ev+100*eps_lanczos)

    CALL cp_dbcsr_scale(matrix_sqrt,scaling,error=error)
    CALL cp_dbcsr_filter(matrix_sqrt,threshold,error=error)

    DO i=1,100

       t1 = m_walltime()

       ! tmp1 = Zk * Yk - I
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_sqrt_inv, matrix_sqrt, 0.0_dp, tmp1,&
                              filter_eps=threshold, flop=flop1, error=error)
       frob_matrix_base=cp_dbcsr_frobenius_norm(tmp1)
       CALL cp_dbcsr_add_on_diag(tmp1,-1.0_dp,error=error)

       ! check convergence (frob norm of what should be the identity matrix minus identity matrix)
       frob_matrix=cp_dbcsr_frobenius_norm(tmp1)

       flop4=0 ; flop5=0
       SELECT CASE(order)
       CASE(2)
          ! update the above to 0.5*(3*I-Zk*Yk)
          CALL cp_dbcsr_add_on_diag(tmp1,-2.0_dp,error=error)
          CALL cp_dbcsr_scale(tmp1,-0.5_dp,error=error)
       CASE(3)
          ! tmp2 = tmp1 ** 2
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, tmp1, tmp1, 0.0_dp, tmp2,&
                                 filter_eps=threshold, flop=flop4, error=error)
          ! tmp1 = 1/16 * (16*I-8*tmp1+6*tmp1**2-5*tmp1**3)
          CALL cp_dbcsr_add(tmp1, tmp2, -4.0_dp, 3.0_dp, error=error)
          CALL cp_dbcsr_add_on_diag(tmp1,8.0_dp,error=error)
          CALL cp_dbcsr_scale(tmp1,0.125_dp,error=error)
       CASE(4) ! as expensive as case(5), so little need to use it
          ! tmp2 = tmp1 ** 2
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, tmp1, tmp1, 0.0_dp, tmp2,&
                                 filter_eps=threshold, flop=flop4, error=error)
          ! tmp3 = tmp2 * tmp1
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp1, 0.0_dp, tmp3,&
                                 filter_eps=threshold, flop=flop5, error=error)
          CALL cp_dbcsr_scale(tmp1,-8.0_dp,error=error)
          CALL cp_dbcsr_add_on_diag(tmp1,16.0_dp,error=error)
          CALL cp_dbcsr_add(tmp1, tmp2, 1.0_dp, 6.0_dp, error=error)
          CALL cp_dbcsr_add(tmp1, tmp3, 1.0_dp,-5.0_dp, error=error)
          CALL cp_dbcsr_scale(tmp1,1/16.0_dp,error=error)
       CASE(5)
          ! Knuth's reformulation to evaluate the polynomial of 4th degree in 2 multiplications
          ! p = y4+A*y3+B*y2+C*y+D
          ! z := y * (y+a); P := (z+y+b) * (z+c) + d.
          ! a=(A-1)/2 ; b=B*(a+1)-C-a*(a+1)*(a+1)
          ! c=B-b-a*(a+1)
          ! d=D-bc
          oa=-40.0_dp/35.0_dp 
          ob= 48.0_dp/35.0_dp
          oc=-64.0_dp/35.0_dp
          od=128.0_dp/35.0_dp
          a=(oa-1)/2
          b=ob*(a+1)-oc-a*(a+1)**2 
          c=ob-b-a*(a+1)
          d=od-b*c
          ! tmp2 = tmp1 ** 2 + a * tmp1
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, tmp1, tmp1, 0.0_dp, tmp2,&
                                 filter_eps=threshold, flop=flop4, error=error)
          CALL cp_dbcsr_add(tmp2, tmp1, 1.0_dp, a , error=error)
          ! tmp3 = tmp2 + tmp1 + b
          CALL cp_dbcsr_copy(tmp3,tmp2,error=error)
          CALL cp_dbcsr_add(tmp3, tmp1, 1.0_dp, 1.0_dp , error=error)
          CALL cp_dbcsr_add_on_diag(tmp3,b,error=error)
          ! tmp2 = tmp2 + c
          CALL cp_dbcsr_add_on_diag(tmp2,c,error=error)
          ! tmp1 = tmp2 * tmp3
          CALL cp_dbcsr_multiply("N", "N", 1.0_dp, tmp2, tmp3, 0.0_dp, tmp1,&
                                 filter_eps=threshold, flop=flop5, error=error)
          ! tmp1 = tmp1 + d
          CALL cp_dbcsr_add_on_diag(tmp1,d,error=error)
          ! final scale
          CALL cp_dbcsr_scale(tmp1,35.0_dp/128.0_dp,error=error)    
       CASE DEFAULT
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT

       ! tmp2 = Yk * tmp1 = Y(k+1)
       CALL cp_dbcsr_multiply("N", "N",  1.0_dp, matrix_sqrt, tmp1,  0.0_dp, tmp2,&
                              filter_eps=threshold, flop=flop2, error=error)
       ! CALL cp_dbcsr_filter(tmp2,threshold,error=error)
       CALL cp_dbcsr_copy(matrix_sqrt, tmp2, error=error)

       ! tmp2 = tmp1 * Zk = Z(k+1)
       CALL cp_dbcsr_multiply("N", "N",  1.0_dp, tmp1, matrix_sqrt_inv,  0.0_dp, tmp2,&
                              filter_eps=threshold,flop=flop3, error=error)
       ! CALL cp_dbcsr_filter(tmp2,threshold,error=error)
       CALL cp_dbcsr_copy(matrix_sqrt_inv, tmp2, error=error)

       occ_matrix=cp_dbcsr_get_occupation(matrix_sqrt_inv)

       ! done iterating
       t2 = m_walltime()

       IF (unit_nr>0) THEN
          WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.8,E12.3,F12.3,F13.3)') "NS sqrt iter ",i,occ_matrix,&
                                                             frob_matrix/frob_matrix_base,t2-t1,&
                                                 (flop1+flop2+flop3+flop4+flop5)/(1.0E6_dp*MAX(0.001_dp,t2-t1))
          CALL m_flush(unit_nr)
       ENDIF

       IF (frob_matrix/frob_matrix_base<SQRT(threshold)) EXIT

    ENDDO

    ! symmetrize the matrices as this is not guaranteed by the algorithm
    CALL cp_dbcsr_transposed(tmp1,matrix_sqrt_inv,error=error)
    CALL cp_dbcsr_add(matrix_sqrt_inv,tmp1,0.5_dp,0.5_dp,error=error)
    CALL cp_dbcsr_transposed(tmp1,matrix_sqrt,error=error)
    CALL cp_dbcsr_add(matrix_sqrt,tmp1,0.5_dp,0.5_dp,error=error)

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
    IF (order.GE.4) THEN
       CALL cp_dbcsr_release(tmp3,error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE matrix_sqrt_Newton_Schulz

! *****************************************************************************
!> \brief ...
!> \param matrix_exp ...
!> \param matrix ...
!> \param omega ...
!> \param alpha ...
!> \param threshold ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE matrix_exponential(matrix_exp,matrix,omega,alpha,threshold,error)
    ! compute matrix_exp=omega*exp(alpha*matrix)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_exp, matrix
    REAL(KIND=dp), INTENT(IN)                :: omega, alpha, threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_exponential', &
      routineP = moduleN//':'//routineN
    REAL(dp), PARAMETER                      :: one = 1.0_dp, &
                                                toll = 1.E-17_dp, &
                                                zero = 0.0_dp

    INTEGER                                  :: handle, i, k, unit_nr
    REAL(dp)                                 :: factorial, norm_C, norm_D, &
                                                norm_scalar
    TYPE(cp_dbcsr_type)                      :: B, B_square, C, D, D_product
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! Calculate the norm of the matrix alpha*matrix, and scale it until it is less than 1.0
    norm_scalar=ABS(alpha)*cp_dbcsr_frobenius_norm(matrix)

    ! k=scaling parameter
    k=1
    DO
      IF((norm_scalar/2.0_dp**k)<=one) EXIT
      k=k+1
    END DO

    ! copy and scale the input matrix in matrix C and in matrix D
    CALL cp_dbcsr_init(C,error=error)
    CALL cp_dbcsr_create(C,template=matrix,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_copy(C,matrix,error=error)
    CALL cp_dbcsr_scale(C, alpha_scalar=alpha/2.0_dp**k, error=error)

    CALL cp_dbcsr_init(D,error=error)
    CALL cp_dbcsr_create(D,template=matrix,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_copy(D,C,error=error)

    !   write(*,*)
    !   write(*,*)
    !   CALL cp_dbcsr_print(D, nodata=.FALSE., matlab_format=.TRUE., variable_name="D", unit_nr=6, error=error)

    ! set the B matrix as B=Identity+D
    CALL cp_dbcsr_init(B,error=error)
    CALL cp_dbcsr_create(B,template=matrix,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_copy(B,D,error=error)
    CALL cp_dbcsr_add_on_diag(B, alpha_scalar=one,  error=error)

    !   CALL cp_dbcsr_print(B, nodata=.FALSE., matlab_format=.TRUE., variable_name="B", unit_nr=6, error=error)

    ! Calculate the norm of C and moltiply by toll to be used as a threshold
    norm_C=toll*cp_dbcsr_frobenius_norm(matrix)

    ! iteration for the trucated taylor series expansion
    CALL cp_dbcsr_init(D_product,error=error)
    CALL cp_dbcsr_create(D_product,template=matrix,matrix_type=dbcsr_type_no_symmetry,error=error)
    i=1
    DO
      i=i+1
      ! compute D_product=D*C
      CALL cp_dbcsr_multiply("N", "N", one, D, C, &
                             zero, D_product, filter_eps=threshold,error=error)

      ! copy D_product in D
      CALL cp_dbcsr_copy(D,D_product,error=error)

      ! calculate B=B+D_product/fat(i)
      factorial=ifac(i)
      CALL cp_dbcsr_add(B, D_product, one, factorial, error=error)

      ! check for convergence using the norm of D (copy of the matrix D_product) and C
      norm_D=factorial*cp_dbcsr_frobenius_norm(D)
      IF(norm_D<norm_C) EXIT
    END DO

    ! start the k iteration for the squaring of the matrix
    CALL cp_dbcsr_init(B_square,error=error)
    CALL cp_dbcsr_create(B_square,template=matrix,matrix_type=dbcsr_type_no_symmetry,error=error)
    DO i=1, k
      !compute B_square=B*B
      CALL cp_dbcsr_multiply("N", "N", one, B, B, &
                             zero, B_square, filter_eps=threshold,error=error)
      ! copy Bsquare in B to iterate
      CALL cp_dbcsr_copy(B,B_square,error=error)
    END DO

    ! copy B_square in matrix_exp and
    CALL cp_dbcsr_copy(matrix_exp,B_square,error=error)

    ! scale matrix_exp by omega, matrix_exp=omega*B_square
    CALL cp_dbcsr_scale(matrix_exp, alpha_scalar=omega, error=error)
    ! write(*,*) alpha,omega

    CALL timestop(handle)

  END SUBROUTINE matrix_exponential

! *****************************************************************************
!> \brief McWeeny purification of a matrix in the orthonormal basis
!> \param matrix_p Matrix to purify (needs to be almost idempotent already)
!> \param threshold Threshold used as filter_eps and convergence criteria
!> \param max_steps Max number of iterations
!> \param error The error-handler
!> \par History
!>       2013.01 created [Florian Schiffmann]
!>       2014.07 slightly refactored [Ole Schuett]
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE purify_mcweeny_orth(matrix_p,threshold,max_steps,error)
    TYPE(cp_dbcsr_type), DIMENSION(:)        :: matrix_p
    REAL(KIND=dp)                            :: threshold
    INTEGER                                  :: max_steps
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mcweeny_orth', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, unit_nr
    REAL(KIND=dp)                            :: frob_norm
    TYPE(cp_dbcsr_type)                      :: matrix_pp, matrix_tmp

    CALL timeset(routineN,handle)
    unit_nr = cp_logger_get_default_unit_nr()

    CALL cp_dbcsr_init(matrix_pp,error=error)
    CALL cp_dbcsr_create(matrix_pp,template=matrix_p(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_tmp,error=error)
    CALL cp_dbcsr_create(matrix_tmp,template=matrix_p(1),matrix_type=dbcsr_type_no_symmetry,error=error)

    DO ispin=1,SIZE(matrix_p)
       DO i=1,max_steps
           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p(ispin), matrix_p(ispin),&
                0.0_dp, matrix_pp, filter_eps=threshold, error=error)

           ! test convergence
           CALL cp_dbcsr_copy(matrix_tmp, matrix_pp, error=error)
           CALL cp_dbcsr_add(matrix_tmp, matrix_p(ispin), 1.0_dp, -1.0_dp, error=error)
           frob_norm = cp_dbcsr_frobenius_norm(matrix_tmp) ! tmp = PP - P
           IF(unit_nr > 0) WRITE(unit_nr,'(t3,a,f16.8)')"McWeeny: Deviation of idempotency",frob_norm
           IF(frob_norm<1000_dp*threshold .AND. i>1) EXIT

           ! construct new P
           CALL cp_dbcsr_copy(matrix_tmp, matrix_pp, error=error)
           CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_pp, matrix_p(ispin),&
                3.0_dp, matrix_tmp, filter_eps=threshold, error=error)
           CALL cp_dbcsr_copy(matrix_p(ispin), matrix_tmp, error=error) ! tmp = 3PP - 2PPP
       END DO
    END DO

    CALL cp_dbcsr_release(matrix_pp,error)
    CALL cp_dbcsr_release(matrix_tmp,error)
    CALL timestop(handle)
  END SUBROUTINE purify_mcweeny_orth

! *****************************************************************************
!> \brief McWeeny purification of a matrix in the non-orthonormal basis
!> \param matrix_p Matrix to purify (needs to be almost idempotent already)
!> \param matrix_s Overlap-Matrix
!> \param threshold Threshold used as filter_eps and convergence criteria
!> \param max_steps Max number of iterations
!> \param error The error-handler
!> \par History
!>       2013.01 created [Florian Schiffmann]
!>       2014.07 slightly refactored [Ole Schuett]
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE purify_mcweeny_nonorth(matrix_p,matrix_s,threshold,max_steps,error)
    TYPE(cp_dbcsr_type), DIMENSION(:)        :: matrix_p
    TYPE(cp_dbcsr_type)                      :: matrix_s
    REAL(KIND=dp)                            :: threshold
    INTEGER                                  :: max_steps
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mcweeny_nonorth', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, unit_nr
    REAL(KIND=dp)                            :: frob_norm
    TYPE(cp_dbcsr_type)                      :: matrix_ps, matrix_psp, &
                                                matrix_test

    CALL timeset(routineN,handle)
    unit_nr = cp_logger_get_default_unit_nr()

    CALL cp_dbcsr_init(matrix_ps,error=error)
    CALL cp_dbcsr_create(matrix_ps,template=matrix_p(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_psp,error=error)
    CALL cp_dbcsr_create(matrix_psp,template=matrix_p(1),matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(matrix_test,error=error)
    CALL cp_dbcsr_create(matrix_test,template=matrix_p(1),matrix_type=dbcsr_type_no_symmetry,error=error)

    DO ispin=1,SIZE(matrix_p)
       DO i=1,max_steps
           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p(ispin), matrix_s,&
                0.0_dp, matrix_ps, filter_eps=threshold,error=error)
           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_ps, matrix_p(ispin),&
                0.0_dp, matrix_psp, filter_eps=threshold,error=error)

           ! test convergence
           CALL cp_dbcsr_copy(matrix_test, matrix_psp, error=error)
           CALL cp_dbcsr_add(matrix_test, matrix_p(ispin), 1.0_dp,-1.0_dp, error=error)
           frob_norm = cp_dbcsr_frobenius_norm(matrix_test) ! test = PSP - P
           IF(unit_nr>0) WRITE(unit_nr,'(t3,a,2f16.8)')"McWeeny: Deviation of idempotency",frob_norm
           IF(frob_norm<1000_dp*threshold .AND. i>1) EXIT

           ! construct new P
           CALL cp_dbcsr_copy(matrix_p(ispin), matrix_psp, error=error)
           CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_ps, matrix_psp,&
                3.0_dp, matrix_p(ispin), filter_eps=threshold, error=error)
       END DO
    END DO

    CALL cp_dbcsr_release(matrix_ps,error)
    CALL cp_dbcsr_release(matrix_psp,error)
    CALL cp_dbcsr_release(matrix_test,error)
    CALL timestop(handle)
  END SUBROUTINE purify_mcweeny_nonorth

END MODULE iterate_matrix
