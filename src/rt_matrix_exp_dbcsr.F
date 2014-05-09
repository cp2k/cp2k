!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for calculating a complex matrix exponential with dbcsr matrices.
!>        Based on the code in rt_matrix_exp.F from Florian Schiffmann
!> \author Samuel Andermatt (02.14)
! *****************************************************************************

MODULE rt_matrix_exp_dbcsr

  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_solve
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_filter, cp_dbcsr_frobenius_norm, cp_dbcsr_init, &
       cp_dbcsr_multiply, cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_transposed
  USE cp_dbcsr_operations,             ONLY: copy_cfm_to_dbcsr,&
                                             copy_dbcsr_to_cfm,&
                                             cp_dbcsr_deallocate_matrix
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_type
  USE dbcsr_api,                       ONLY: dbcsr_type_complex_8
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fac
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_matrix_exp_dbcsr'

  PUBLIC :: taylor_only_imaginary_dbcsr,&
            taylor_full_complex_dbcsr,&
            exp_pade_full_complex_dbcsr,&
            exp_pade_only_imaginary_dbcsr,&
            cp_complex_dbcsr_gemm,&
            cp_complex_dbcsr_gemm_3,&
            bch_expansion_imaginary_propagator,&
            bch_expansion_complex_propagator

CONTAINS

! *****************************************************************************
!> \brief Convenience function. Computes the matrix multiplications needed
!>        for the multiplication of complex sparse matrices.
!>        C = beta * C + alpha * ( A  ** transa ) * ( B ** transb )
!> \param A : m x k matrix ( ! for transa = 'N'), A_re is the real part, A_im the imaginary part
!> \param B : k x n matrix ( ! for transb = 'N')
!> \param C : m x n matrix
!> \param transa : 'N' -> normal   'T' -> transpose
!>      alpha,beta :: can be 0.0_dp and 1.0_dp
!> \param b_first_col : the k x n matrix starts at col b_first_col of matrix_b (avoid usage)
!> \note
!>      C should have no overlap with A, B
!>      This subroutine uses three real matrix multiplications instead of two complex
!>      This reduces the amount of flops and memory bandwidth by 25%, but for memory bandwidth
!>      true complex algebra is still superior (one third less bandwidth needed)
!>      limited cases matrix multiplications
!> \author Samuel Andermatt
! *****************************************************************************

SUBROUTINE cp_complex_dbcsr_gemm_3(transa, transb, alpha, A_re, A_im,&
           B_re, B_im, beta, C_re, C_im, filter_eps, error)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_dbcsr_type), INTENT(IN)          :: A_re, A_im, B_re, B_im
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: C_re, C_im
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: filter_eps
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_complex_dbcsr_gemm_3', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    CHARACTER                                :: matrix_type
    CHARACTER(LEN=1)                         :: transa2, transb2
    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha2, alpha3, alpha4
    TYPE(cp_dbcsr_type), POINTER             :: a_plus_b, ac, bd, c_plus_d

    CALL timeset(routineN,handle)
    !A complex matrix matrix multiplication can be done with only three multiplications
    !(a+ib)*(c+id)=ac-bd+i((a+b)*(c+d) - ac - bd)
    !A_re=a, A_im=b, B_re=c, B_im=d

    failure=.FALSE.

    alpha2=-alpha
    alpha3=alpha
    alpha4=alpha

    IF(transa=="C") THEN
       alpha2=-alpha2
       alpha3=-alpha3
       transa2="T"
    ELSE
       transa2=transa
    END IF
    IF(transb=="C") THEN
       alpha2=-alpha2
       alpha4=-alpha4
       transb2="T"
    ELSE
       transb2=transb
    END IF

    !create the work matrices
    NULLIFY(ac)
    ALLOCATE(ac,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(ac,error=error)
    CALL cp_dbcsr_create(ac,template=A_re,matrix_type="N",error=error)
    NULLIFY(bd)
    ALLOCATE(bd,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(bd,error=error)
    CALL cp_dbcsr_create(bd,template=A_re,matrix_type="N",error=error)
    NULLIFY(a_plus_b)
    ALLOCATE(a_plus_b,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(a_plus_b,error=error)
    CALL cp_dbcsr_create(a_plus_b,template=A_re,matrix_type="N",error=error)
    NULLIFY(c_plus_d)
    ALLOCATE(c_plus_d,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(c_plus_d,error=error)
    CALL cp_dbcsr_create(c_plus_d,template=A_re,matrix_type="N",error=error)
    
    !Do the neccesarry multiplications
    CALL cp_dbcsr_multiply(transa2,transb2,alpha,A_re,B_re,zero,ac,filter_eps=filter_eps,error=error)
    CALL cp_dbcsr_multiply(transa2,transb2,alpha2,A_im,B_im,zero,bd,filter_eps=filter_eps,error=error)

    CALL cp_dbcsr_add(a_plus_b,A_re,zero,alpha,error=error)
    CALL cp_dbcsr_add(a_plus_b,A_im,one,alpha3,error=error)
    CALL cp_dbcsr_add(c_plus_d,B_re,zero,alpha,error=error)
    CALL cp_dbcsr_add(c_plus_d,B_im,one,alpha4,error=error)

    !this can already be written into C_im
    !now both matrixes have been scaled which means we currently multiplied by alpha squared
    CALL cp_dbcsr_multiply(transa2,transb2,one/alpha,a_plus_b,c_plus_d,beta,C_im,filter_eps=filter_eps,error=error)

    !now add up all the terms into the result
    CALL cp_dbcsr_add(C_re,ac,beta,one,error=error)
    !the minus sign was already taken care of at the definition of alpha2
    CALL cp_dbcsr_add(C_re,bd,one,one,error=error)
    CALL cp_dbcsr_filter(C_re,filter_eps,error=error)

    CALL cp_dbcsr_add(C_im,ac,one,-one,error=error)
    !the minus sign was already taken care of at the definition of alpha2
    CALL cp_dbcsr_add(C_im,bd,one,one,error=error)
    CALL cp_dbcsr_filter(C_im,filter_eps,error=error)

    !Deallocate the work matrices
    CALL cp_dbcsr_deallocate_matrix(ac,error=error)
    CALL cp_dbcsr_deallocate_matrix(bd,error=error)
    CALL cp_dbcsr_deallocate_matrix(a_plus_b,error=error)
    CALL cp_dbcsr_deallocate_matrix(c_plus_d,error=error)

    CALL timestop(handle)

END SUBROUTINE

! *****************************************************************************
!> \brief Convenience function. Computes the matrix multiplications needed
!>        for the multiplication of complex sparse matrices.
!>        C = beta * C + alpha * ( A  ** transa ) * ( B ** transb )
!> \param A : m x k matrix ( ! for transa = 'N'), A_re is the real part, A_im the imaginary part
!> \param B : k x n matrix ( ! for transb = 'N')
!> \param C : m x n matrix
!> \param transa : 'N' -> normal   'T' -> transpose
!>      alpha,beta :: can be 0.0_dp and 1.0_dp
!> \param b_first_col : the k x n matrix starts at col b_first_col of matrix_b (avoid usage)
!> \note
!>      C should have no overlap with A, B
!> \author Samuel Andermatt
! *****************************************************************************

SUBROUTINE cp_complex_dbcsr_gemm(transa, transb, alpha, A_re, A_im,&
           B_re, B_im, beta, C_re, C_im,&
           first_row, last_row, first_column, last_column, first_k, last_k,&
           retain_sparsity, filter_eps, error)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_dbcsr_type), INTENT(IN)          :: A_re, A_im, B_re, B_im
    REAL(KIND=dp), INTENT(IN)                :: beta
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: C_re, C_im
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: filter_eps
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_complex_dbcsr_gemm', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    CHARACTER(LEN=1)                         :: transa2, transb2
    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha2, alpha3, alpha4

    CALL timeset(routineN,handle)
    !A complex matrix matrix multiplication can be done with only three multiplications
    !(a+ib)*(c+id)=ac-bd+i((a+b)*(c+d) - ac - bd)
    !A_re=a, A_im=b, B_re=c, B_im=d

    failure=.FALSE.

    alpha2=-alpha
    alpha3=alpha
    alpha4=alpha

    IF(transa=="C") THEN
       alpha2=-alpha2
       alpha3=-alpha3
       transa2="T"
    ELSE
       transa2=transa
    END IF
    IF(transb=="C") THEN
       alpha2=-alpha2
       alpha4=-alpha4
       transb2="T"
    ELSE
       transb2=transb
    END IF
    
    CALL cp_dbcsr_multiply(transa2,transb2,alpha,A_re,B_re,beta,C_re, filter_eps=filter_eps,error=error)
    CALL cp_dbcsr_multiply(transa2,transb2,alpha2,A_im,B_im,one,C_re, filter_eps=filter_eps,error=error)
    CALL cp_dbcsr_multiply(transa2,transb2,alpha3,A_im,B_re,beta,C_im,filter_eps=filter_eps,error=error)
    CALL cp_dbcsr_multiply(transa2,transb2,alpha4,A_re,B_im,one,C_im,filter_eps=filter_eps,error=error)


    CALL timestop(handle)

END SUBROUTINE

! *****************************************************************************
!> \brief specialized subroutine for purely imaginary matrix exponentials
!> \author Samuel Andermatt (02.2014)
! *****************************************************************************

  SUBROUTINE taylor_only_imaginary_dbcsr(exp_H,im_matrix,nsquare,ntaylor,filter_eps,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(2)      :: exp_H
    TYPE(cp_dbcsr_type), POINTER             :: im_matrix
    INTEGER                                  :: nsquare, ntaylor
    REAL(KIND=dp)                            :: filter_eps
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'taylor_only_imaginary_dbcsr', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, n, nloop, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: square_fac, Tfac, tmp
    TYPE(cp_dbcsr_type), POINTER             :: T1, T2, Tres_im, Tres_re

    CALL timeset(routineN,handle)

    failure=.FALSE.

    !The divider that comes from the scaling and squaring procedure
    square_fac=1.0_dp/(2.0_dp**REAL(nsquare,dp))

    !Allocate work matrices
    NULLIFY(T1)
    ALLOCATE(T1,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(T1,error=error)
    CALL cp_dbcsr_create(T1,template=im_matrix,matrix_type="N",error=error)
    NULLIFY(T2)
    ALLOCATE(T2,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(T2,error=error)
    CALL cp_dbcsr_create(T2,template=im_matrix,matrix_type="N",error=error)
    NULLIFY(Tres_re)
    ALLOCATE(Tres_re,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(Tres_re,error=error)
    CALL cp_dbcsr_create(Tres_re,template=im_matrix,matrix_type="N",error=error)
    NULLIFY(Tres_im)
    ALLOCATE(Tres_im,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(Tres_im,error=error)
    CALL cp_dbcsr_create(Tres_im,template=im_matrix,matrix_type="N",error=error)

    !Create the unit matrices
    CALL cp_dbcsr_set(T1,zero,error=error)
    CALL cp_dbcsr_add_on_diag(T1,one,error=error)
    CALL cp_dbcsr_set(Tres_re,zero,error=error)
    CALL cp_dbcsr_add_on_diag(Tres_re,one,error=error)
    CALL cp_dbcsr_set(Tres_im,zero,error=error)

    nloop=CEILING(REAL(ntaylor,dp)/2.0_dp)
    !the inverse of the prefactor in the taylor series
    tmp=1.0_dp
    DO i=1,nloop
       CALL cp_dbcsr_scale(T1,1/(REAL(i,dp)*2.0_dp-1.0_dp),error=error)
       CALL cp_dbcsr_filter(T1,filter_eps,error=error)
       CALL cp_dbcsr_multiply("N","N",square_fac,im_matrix,T1,zero,&
            T2,filter_eps=filter_eps,error=error)
       Tfac=one
       IF(MOD(i,2)==0) Tfac=-Tfac
       CALL cp_dbcsr_add(Tres_im,T2,one,Tfac,error=error)
       CALL cp_dbcsr_scale(T2,1/(REAL(i,dp)*2.0_dp),error=error)
       CALL cp_dbcsr_filter(T2,filter_eps,error=error)
       CALL cp_dbcsr_multiply("N","N",square_fac,im_matrix,T2,zero,&
            T1,filter_eps=filter_eps,error=error)
       Tfac=one
       IF(MOD(i,2)==1) Tfac=-Tfac
       CALL cp_dbcsr_add(Tres_re,T1,one,Tfac,error=error)
    END DO

    !Square the matrices, due to the scaling and squaring procedure
    IF(nsquare.GT.0)THEN
       DO i=1,nsquare
          CALL cp_complex_dbcsr_gemm("N","N",one,Tres_re,Tres_im,&
               Tres_re,Tres_im,zero,exp_H(1)%matrix,exp_H(2)%matrix,&
               filter_eps=filter_eps,error=error)
          CALL cp_dbcsr_copy(Tres_re,exp_H(1)%matrix,error=error)
          CALL cp_dbcsr_copy(Tres_im,exp_H(2)%matrix,error=error)
       END DO
    ELSE
       CALL cp_dbcsr_copy(exp_H(1)%matrix,Tres_re,error=error)
       CALL cp_dbcsr_copy(exp_H(2)%matrix,Tres_im,error=error)
    END IF
    CALL cp_dbcsr_deallocate_matrix(T1,error=error)
    CALL cp_dbcsr_deallocate_matrix(T2,error=error)
    CALL cp_dbcsr_deallocate_matrix(Tres_re,error=error)
    CALL cp_dbcsr_deallocate_matrix(Tres_im,error=error)

    CALL timestop(handle)


  END SUBROUTINE taylor_only_imaginary_dbcsr


! *****************************************************************************
!> \brief subroutine for general complex matrix exponentials
!>        on input a separate cp_fm_type for real and complex part
!>        on output a size 2 cp_fm_p_type, first element is the real part of
!>        the exponential second the imaginary
!> \author Samuel Andermatt (02.2014)
! *****************************************************************************
  SUBROUTINE taylor_full_complex_dbcsr(exp_H,re_part,im_part,nsquare,ntaylor,filter_eps,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(2)      :: exp_H
    TYPE(cp_dbcsr_type), POINTER             :: re_part, im_part
    INTEGER                                  :: nsquare, ntaylor
    REAL(KIND=dp)                            :: filter_eps
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'taylor_full_complex_dbcsr', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp)                         :: square_fac
    INTEGER                                  :: handle, i, stat
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type), POINTER             :: T1, T2, T3, Tres

    CALL timeset(routineN,handle)

    failure=.FALSE.

    !The divider that comes from the scaling and squaring procedure
    square_fac=CMPLX(1.0_dp/(2.0_dp**REAL(nsquare,dp)),0.0_dp,KIND=dp)

    !Allocate work matrices
    NULLIFY(T1)
    ALLOCATE(T1,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(T1,error=error)
    CALL cp_dbcsr_create(T1,template=re_part,matrix_type="N",&
         data_type=dbcsr_type_complex_8,error=error)
    NULLIFY(T2)
    ALLOCATE(T2,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(T2,error=error)
    CALL cp_dbcsr_create(T2,template=re_part,matrix_type="N",&
         data_type=dbcsr_type_complex_8,error=error)
    NULLIFY(T3)
    ALLOCATE(T3,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(T3,error=error)
    CALL cp_dbcsr_create(T3,template=re_part,matrix_type="N",&
         data_type=dbcsr_type_complex_8,error=error)
    NULLIFY(Tres)
    ALLOCATE(Tres,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(Tres,error=error)
    CALL cp_dbcsr_create(Tres,template=re_part,matrix_type="N",&
         data_type=dbcsr_type_complex_8,error=error)

    !Fuse the input matrices to a single complex matrix
    CALL cp_dbcsr_copy(T3,re_part,error=error)
    CALL cp_dbcsr_copy(Tres,im_part,error=error) !will later on be set back to zero
    CALL cp_dbcsr_scale(Tres,CMPLX(0.0_dp,1.0_dp,KIND=dp),error=error)
    CALL cp_dbcsr_add(T3,Tres,one,one,error=error)

    !Create the unit matrices
    CALL cp_dbcsr_set(T1,zero,error=error)
    CALL cp_dbcsr_add_on_diag(T1,1.0_dp,error=error)
    CALL cp_dbcsr_set(Tres,zero,error=error)
    CALL cp_dbcsr_add_on_diag(Tres,1.0_dp,error=error)

    DO i=1,ntaylor
       CALL cp_dbcsr_scale(T1,one/REAL(i,dp),error=error)
       CALL cp_dbcsr_filter(T1,filter_eps,error=error)
       CALL cp_dbcsr_multiply("N","N",square_fac,T1,T3,&
            zero,T2,filter_eps=filter_eps,error=error)
       CALL cp_dbcsr_add(Tres,T2,one,one,error=error)
       CALL cp_dbcsr_copy(T1,T2,error=error)
    END DO

    IF(nsquare.GT.0)THEN
       DO i=1,nsquare
          CALL cp_dbcsr_multiply("N","N",one,Tres,Tres,zero,&
               T2,filter_eps=filter_eps,error=error)
          CALL cp_dbcsr_copy(Tres,T2,error=error)
       END DO
    END IF

    CALL cp_dbcsr_copy(exp_H(1)%matrix,Tres,keep_imaginary=.FALSE.,error=error)
    CALL cp_dbcsr_scale(Tres,CMPLX(0.0_dp,-1.0_dp,KIND=dp),error=error)
    CALL cp_dbcsr_copy(exp_H(2)%matrix,Tres,keep_imaginary=.FALSE.,error=error)


    CALL cp_dbcsr_deallocate_matrix(T1,error=error)
    CALL cp_dbcsr_deallocate_matrix(T2,error=error)
    CALL cp_dbcsr_deallocate_matrix(T3,error=error)
    CALL cp_dbcsr_deallocate_matrix(Tres,error=error)

    CALL timestop(handle)

  END SUBROUTINE taylor_full_complex_dbcsr

! *****************************************************************************
!> \brief exponential of a complex matrix,
!>        calculated using pade approximation together with scaling and squaring
!> \author Samuel Andermatt (02.2014)
! *****************************************************************************

  SUBROUTINE exp_pade_full_complex_dbcsr(exp_H,re_part,im_part,nsquare,npade,filter_eps,fm_struct,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(2)      :: exp_H
    TYPE(cp_dbcsr_type), POINTER             :: re_part, im_part
    INTEGER                                  :: nsquare, npade
    REAL(KIND=dp)                            :: filter_eps
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'exp_pade_full_complex_dbcsr', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp)                         :: scaleD, scaleN
    INTEGER                                  :: handle, i, p, q, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: square_fac
    TYPE(cp_cfm_type), POINTER               :: Dpq_fm, Npq_fm
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: mult_p
    TYPE(cp_dbcsr_type), POINTER             :: Dpq, fin_p, Npq, T1, T2, tmp, &
                                                Tres

    CALL timeset(routineN,handle)

    failure = .FALSE.

    p=npade
    q=npade

    !Allocate work matrices
    NULLIFY(T1)
    ALLOCATE(T1,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(T1,error=error)
    CALL cp_dbcsr_create(T1,template=re_part,matrix_type="N",data_type=dbcsr_type_complex_8,error=error)
    NULLIFY(T2)
    ALLOCATE(T2,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(T2,error=error)
    CALL cp_dbcsr_create(T2,template=re_part,matrix_type="N",data_type=dbcsr_type_complex_8,error=error)
    NULLIFY(Tres)
    ALLOCATE(Tres,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(Tres,error=error)
    CALL cp_dbcsr_create(Tres,template=re_part,matrix_type="N",data_type=dbcsr_type_complex_8,error=error)
    NULLIFY(Dpq)
    ALLOCATE(Dpq,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(Dpq,error=error)
    CALL cp_dbcsr_create(Dpq,template=re_part,matrix_type="N",data_type=dbcsr_type_complex_8,error=error)
    NULLIFY(Npq)
    ALLOCATE(Npq,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(Npq,error=error)
    CALL cp_dbcsr_create(Npq,template=re_part,matrix_type="N",data_type=dbcsr_type_complex_8,error=error)
    NULLIFY(tmp)
    ALLOCATE(tmp,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(tmp,error=error)
    CALL cp_dbcsr_create(tmp,template=re_part,matrix_type="N",data_type=dbcsr_type_complex_8,error=error)


    NULLIFY(mult_p)
    ALLOCATE(mult_p(2),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL cp_dbcsr_set(T1,zero,error=error)
    CALL cp_dbcsr_set(T2,zero,error=error)

    square_fac=2.0_dp**REAL(nsquare,dp)

    CALL cp_dbcsr_set(tmp,zero,error=error)
    CALL cp_dbcsr_set(T2,zero,error=error)
    CALL cp_dbcsr_copy(T2,re_part,error=error)
    CALL cp_dbcsr_copy(tmp,im_part,error=error)
    CALL cp_dbcsr_scale(tmp,CMPLX(0.0_dp,1.0_dp,KIND=dp),error=error)
    CALL cp_dbcsr_add(T2,tmp,one,one,error=error)
    CALL cp_dbcsr_scale(T2,CMPLX(1/square_fac,0.0_dp,KIND=dp),error=error)

    CALL cp_dbcsr_copy(T1,T2,error=error)

    mult_p(1)%matrix=>T2
    mult_p(2)%matrix=>Tres

    !Create the unit matrices
    CALL cp_dbcsr_set(Dpq,zero,error=error)
    CALL cp_dbcsr_add_on_diag(Dpq,1.0_dp,error=error)
    CALL cp_dbcsr_set(Npq,zero,error=error)
    CALL cp_dbcsr_add_on_diag(Npq,1.0_dp,error=error)

    CALL cp_dbcsr_add(Npq,T2,one,one*0.5_dp,error=error)
    CALL cp_dbcsr_add(Dpq,T2,one,-one*0.5_dp,error=error)

    npade=MAX(p,q) !p will always be equal q here
    IF(npade.GT.2)THEN
       DO i=2,npade
          IF(i.LE.p)scaleN=CMPLX(fac(p +q -i )*fac(p )/(fac(p +q )*fac(i )*fac(p -i )),0.0_dp,kind=dp)
          scaleD=CMPLX((-1.0_dp)**i  * fac(p +q -i )*fac(q)/(fac(p +q )*fac(i )*fac(q -i)),0.0_dp,kind=dp)
          CALL cp_dbcsr_multiply("N","N",one,T1,mult_p(MOD(i,2)+1)%matrix,zero,&
               mult_p(MOD(i+1,2)+1)%matrix,filter_eps=filter_eps,error=error)
          IF(i.LE.p)CALL cp_dbcsr_add(Npq,mult_p(MOD(i+1,2)+1)%matrix,one,scaleN,error=error)
          IF(i.LE.q)CALL cp_dbcsr_add(Dpq,mult_p(MOD(i+1,2)+1)%matrix,one,scaleD,error=error)
       END DO
    END IF


!   Inversion of Dpq using full matrix algebra
    CALL cp_cfm_create(Dpq_fm, matrix_struct=fm_struct,&
         name="Dpq_fm",error=error)
    CALL cp_cfm_create(Npq_fm, matrix_struct=fm_struct,&
         name="Npq_fm",error=error)

    CALL copy_dbcsr_to_cfm(Dpq,Dpq_fm,error=error)
    CALL copy_dbcsr_to_cfm(Npq,Npq_fm,error=error)

    CALL cp_cfm_solve(Dpq_fm,Npq_fm,error=error)

    CALL copy_cfm_to_dbcsr(Dpq_fm,Dpq,error=error)
    CALL copy_cfm_to_dbcsr(Npq_fm,Npq,error=error)
    CALL cp_dbcsr_filter(Dpq,filter_eps,error=error)
    CALL cp_dbcsr_filter(Npq,filter_eps,error=error)

    CALL cp_cfm_release(Npq_fm,error=error)
    CALL cp_cfm_release(Dpq_fm,error=error)

    mult_p(2)%matrix=>Npq
    mult_p(1)%matrix=>Tres
    IF(nsquare.GT.0)THEN
       DO i=1,nsquare
          CALL cp_dbcsr_multiply("N","N",one,mult_p(MOD(i,2)+1)%matrix,mult_p(MOD(i,2)+1)%matrix,zero,&
               mult_p(MOD(i+1,2)+1)%matrix,filter_eps=filter_eps,error=error)
          fin_p=> mult_p(MOD(i+1,2)+1)%matrix
       END DO
    ELSE
       fin_p=>Npq
    END IF

    CALL cp_dbcsr_copy(exp_H(1)%matrix,fin_p,keep_imaginary=.FALSE.,error=error)
    CALL cp_dbcsr_scale(fin_p,CMPLX(0.0_dp,-1.0_dp,KIND=dp),error=error)
    CALL cp_dbcsr_copy(exp_H(2)%matrix,fin_p,keep_imaginary=.FALSE.,error=error)

    CALL cp_dbcsr_deallocate_matrix(tmp,error=error)
    CALL cp_dbcsr_deallocate_matrix(T1,error=error)
    CALL cp_dbcsr_deallocate_matrix(T2,error=error)
    CALL cp_dbcsr_deallocate_matrix(Tres,error=error)
    CALL cp_dbcsr_deallocate_matrix(Npq,error=error)
    CALL cp_dbcsr_deallocate_matrix(Dpq,error=error)
    DEALLOCATE(mult_p)
    CALL timestop(handle)

  END SUBROUTINE exp_pade_full_complex_dbcsr

! *****************************************************************************
!> \brief exponential of a complex matrix,
!>        calculated using pade approximation together with scaling and squaring
!> \author Samuel Andermatt (02.2014)
! *****************************************************************************

  SUBROUTINE exp_pade_only_imaginary_dbcsr(exp_H,im_part,nsquare,npade,filter_eps,fm_struct,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(2)      :: exp_H
    TYPE(cp_dbcsr_type), POINTER             :: im_part
    INTEGER                                  :: nsquare, npade
    REAL(KIND=dp)                            :: filter_eps
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'exp_pade_only_imaginary_dbcsr', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)
    REAL(KIND=dp), PARAMETER                 :: rone = (1.0_dp,0.0_dp) , &
                                                rzero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp)                         :: scaleD, scaleN
    INTEGER                                  :: handle, i, j, p, q, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_fac, square_fac
    TYPE(cp_cfm_type), POINTER               :: Dpq_fm, Npq_fm
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: cmult_p
    TYPE(cp_dbcsr_type), POINTER             :: Dpq, fin_p, Npq, T1, T2, tmp, &
                                                Tres

    CALL timeset(routineN,handle)
    p=npade 
    q=npade !p==q seems to be neccesary for the rest of the code

    ALLOCATE(cmult_p(2))

    !Allocate work matrices
    NULLIFY(T1)
    ALLOCATE(T1,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(T1,error=error)
    CALL cp_dbcsr_create(T1,template=im_part,matrix_type="N",data_type=dbcsr_type_complex_8,error=error)
    NULLIFY(T2)
    ALLOCATE(T2,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(T2,error=error)
    CALL cp_dbcsr_create(T2,template=im_part,matrix_type="N",error=error)
    NULLIFY(Tres)
    ALLOCATE(Tres,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(Tres,error=error)
    CALL cp_dbcsr_create(Tres,template=im_part,matrix_type="N",error=error)
    NULLIFY(Dpq)
    ALLOCATE(Dpq,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(Dpq,error=error)
    CALL cp_dbcsr_create(Dpq,template=im_part,matrix_type="N",data_type=dbcsr_type_complex_8,error=error)
    NULLIFY(Npq)
    ALLOCATE(Npq,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(Npq,error=error)
    CALL cp_dbcsr_create(Npq,template=im_part,matrix_type="N",data_type=dbcsr_type_complex_8,error=error)
    NULLIFY(tmp)
    ALLOCATE(tmp,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(tmp,error=error)
    CALL cp_dbcsr_create(tmp,template=im_part,matrix_type="N",data_type=dbcsr_type_complex_8,error=error)

    CALL cp_dbcsr_copy(T2,im_part,error=error)

    !Create the unit matrices
    CALL cp_dbcsr_set(T1,zero,error=error)
    CALL cp_dbcsr_set(Dpq,zero,error=error)
    CALL cp_dbcsr_add_on_diag(Dpq,1.0_dp,error=error)
    CALL cp_dbcsr_set(Npq,zero,error=error)
    CALL cp_dbcsr_add_on_diag(Npq,1.0_dp,error=error)

    square_fac=1/(2.0_dp**REAL(nsquare,dp))

    CALL cp_dbcsr_copy(tmp,im_part,error=error)
    CALL cp_dbcsr_scale(tmp,CMPLX(0.0_dp,square_fac,KIND=dp),error=error)

    CALL cp_dbcsr_add(Npq,tmp,one,one*0.5_dp,error=error)
    CALL cp_dbcsr_add(Dpq,tmp,one,-one*0.5_dp,error=error)

    IF(npade.GT.2)THEN
       DO j=1,FLOOR(npade/2.0_dp)
          i=2*j
          my_fac=(-1.0_dp)**j
          IF(i.LE.p)scaleN=my_fac*fac(p +q -i )*fac(p )/(fac(p +q )*fac(i )*fac(p -i ))
          scaleD=my_fac* fac(p +q -i )*fac(q)/(fac(p +q )*fac(i )*fac(q -i))
          CALL cp_dbcsr_multiply("N","N",square_fac,im_part,T2,rzero,Tres,filter_eps=filter_eps,error=error)
          CALL cp_dbcsr_copy(tmp,Tres,error=error) !This makes it a complex matrix so it can be added
          CALL cp_dbcsr_add(Npq,tmp,one,scaleN,error=error)
          CALL cp_dbcsr_add(Dpq,tmp,one,scaleD,error=error)

          IF(2*j+1.le.q)THEN
             i=2*j+1
             IF(i.LE.p)scaleN=my_fac*fac(p +q -i )*fac(p )/(fac(p +q )*fac(i )*fac(p -i ))
             scaleD=-my_fac*fac(p +q -i )*fac(q)/(fac(p +q )*fac(i )*fac(q -i))
             CALL cp_dbcsr_multiply("N","N",square_fac,im_part,Tres,rzero,T2,filter_eps=filter_eps,error=error)
             CALL cp_dbcsr_copy(tmp,T2,error=error) !This makes it a complex matrix so it can be added
             CALL cp_dbcsr_add(Npq,tmp,one,scaleN*CMPLX(0.0_dp,1.0_dp,KIND=dp),error=error)
             CALL cp_dbcsr_add(Dpq,tmp,one,scaleD*CMPLX(0.0_dp,1.0_dp,KIND=dp),error=error)
          ENDIF
       END DO
    END IF

!**********************************************************************************************************
    !Inversion currently fm

    CALL cp_cfm_create(Dpq_fm, matrix_struct=fm_struct,&
         name="Dpq_fm",error=error)
    CALL cp_cfm_create(Npq_fm, matrix_struct=fm_struct,&
         name="Npq_fm",error=error)

    CALL copy_dbcsr_to_cfm(Dpq,Dpq_fm,error=error)
    CALL copy_dbcsr_to_cfm(Npq,Npq_fm,error=error)

    CALL cp_cfm_solve(Dpq_fm,Npq_fm,error=error)

    CALL copy_cfm_to_dbcsr(Dpq_fm,Dpq,error=error)
    CALL copy_cfm_to_dbcsr(Npq_fm,Npq,error=error)
    CALL cp_dbcsr_filter(Dpq,filter_eps,error=error)
    CALL cp_dbcsr_filter(Npq,filter_eps,error=error)

    CALL cp_cfm_release(Npq_fm,error=error)
    CALL cp_cfm_release(Dpq_fm,error=error)
!***********************************************************************************************************************


    cmult_p(2)%matrix=>Npq
    cmult_p(1)%matrix=>T1
    IF(nsquare.GT.0)THEN
       DO i=1,nsquare
          CALL cp_dbcsr_multiply("N","N",one,cmult_p(MOD(i,2)+1)%matrix,cmult_p(MOD(i,2)+1)%matrix,zero,&
               cmult_p(MOD(i+1,2)+1)%matrix,filter_eps=filter_eps,error=error)
          fin_p=> cmult_p(MOD(i+1,2)+1)%matrix
       END DO
    ELSE
       fin_p=>Npq
    END IF

    CALL cp_dbcsr_copy(exp_H(1)%matrix,fin_p,keep_imaginary=.FALSE.,error=error)
    CALL cp_dbcsr_scale(fin_p,CMPLX(0.0_dp,-1.0_dp,KIND=dp),error=error)
    CALL cp_dbcsr_copy(exp_H(2)%matrix,fin_p,keep_imaginary=.FALSE.,error=error)

    CALL cp_dbcsr_deallocate_matrix(tmp,error=error)
    CALL cp_dbcsr_deallocate_matrix(T1,error=error)
    CALL cp_dbcsr_deallocate_matrix(T2,error=error)
    CALL cp_dbcsr_deallocate_matrix(Tres,error=error)
    CALL cp_dbcsr_deallocate_matrix(Npq,error=error)
    CALL cp_dbcsr_deallocate_matrix(Dpq,error=error)
    DEALLOCATE(cmult_p)
    CALL timestop(handle)

  END SUBROUTINE exp_pade_only_imaginary_dbcsr

! *****************************************************************************
!> \brief  The Baker–Campbell–Hausdorff expansion for a purely imaginary propagator (rtp)
!> \author Samuel Andermatt (02.2014)
! *****************************************************************************

   SUBROUTINE bch_expansion_imaginary_propagator(propagator,density_re,density_im,filter_eps,eps_exp,error)
    TYPE(cp_dbcsr_type), POINTER             :: propagator, density_re, &
                                                density_im
    REAL(KIND=dp)                            :: filter_eps, eps_exp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'bch_expansion_imaginary_propagator', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, stat, unit_nr
    LOGICAL                                  :: convergence, failure
    REAL(KIND=dp)                            :: alpha, max_alpha, prefac
    TYPE(cp_dbcsr_type), POINTER             :: comm, comm2, tmp, tmp2
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    unit_nr = cp_logger_get_default_io_unit(logger)

    failure=.FALSE.
    NULLIFY(tmp)
    ALLOCATE(tmp,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(tmp,error=error)
    CALL cp_dbcsr_create(tmp,template=propagator,error=error)
    NULLIFY(tmp2)
    ALLOCATE(tmp2,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(tmp2,error=error)
    CALL cp_dbcsr_create(tmp2,template=propagator,error=error)
    NULLIFY(comm)
    ALLOCATE(comm,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(comm,error=error)
    CALL cp_dbcsr_create(comm,template=propagator,error=error)
    NULLIFY(comm2)
    ALLOCATE(comm2,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(comm2,error=error)
    CALL cp_dbcsr_create(comm2,template=propagator,error=error)

    CALL cp_dbcsr_copy(tmp,density_re,error=error)
    CALL cp_dbcsr_copy(tmp2,density_im,error=error)

    convergence=.FALSE.
    DO i=1,20
       prefac=one/i
       CALL cp_dbcsr_multiply("N","N",-prefac,propagator,tmp2,zero,comm,&
            filter_eps=filter_eps,error=error)
       CALL cp_dbcsr_multiply("N","N",prefac,propagator,tmp,zero,comm2,&
            filter_eps=filter_eps,error=error)
       CALL cp_dbcsr_transposed(tmp,comm,error=error)
       CALL cp_dbcsr_transposed(tmp2,comm2,error=error)
       CALL cp_dbcsr_add(comm,tmp,one,one,error=error)
       CALL cp_dbcsr_add(comm2,tmp2,one,-one,error=error)
       CALL cp_dbcsr_add(density_re,comm,one,one,error=error)
       CALL cp_dbcsr_add(density_im,comm2,one,one,error=error)
       CALL cp_dbcsr_copy(tmp,comm,error=error)
       CALL cp_dbcsr_copy(tmp2,comm2,error=error)
       !check for convergence
       max_alpha=zero
       alpha = cp_dbcsr_frobenius_norm(comm)
       IF(alpha>max_alpha) max_alpha=alpha
       alpha = cp_dbcsr_frobenius_norm(comm2)
       IF(alpha>max_alpha) max_alpha=alpha
       IF(max_alpha<eps_exp) convergence = .TRUE.
       IF(convergence) THEN
          IF(unit_nr>0) WRITE (UNIT=unit_nr,FMT="((T3,A,I2,A))")&
                "BCH converged after ",i," steps"
          EXIT
       ENDIF
    ENDDO

    CALL cp_dbcsr_filter(density_re,filter_eps,error=error)
    CALL cp_dbcsr_filter(density_im,filter_eps,error=error)

    CALL cp_assert(convergence,cp_warning_level,cp_assertion_failed,routineP,&
         "BCH method did not converge"//&
         CPSourceFileRef,&
         only_ionode=.TRUE.)

    CALL cp_dbcsr_deallocate_matrix(tmp,error=error)
    CALL cp_dbcsr_deallocate_matrix(tmp2,error=error)
    CALL cp_dbcsr_deallocate_matrix(comm,error=error)
    CALL cp_dbcsr_deallocate_matrix(comm2,error=error)

    CALL timestop(handle)

   END SUBROUTINE

! *****************************************************************************
!> \brief  The Baker-Campbell-Hausdorff expansion for a purely imaginary propagator (rtp)
!>         Works for a non unitary propagator, because the density matrix is hermitian
!> \author Samuel Andermatt (02.2014)
! *****************************************************************************

   SUBROUTINE bch_expansion_complex_propagator(propagator_re,propagator_im,density_re,density_im,filter_eps,eps_exp,error)
    TYPE(cp_dbcsr_type), POINTER             :: propagator_re, propagator_im, &
                                                density_re, density_im
    REAL(KIND=dp)                            :: filter_eps, eps_exp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'bch_expansion_complex_propagator', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, stat, unit_nr
    LOGICAL                                  :: convergence, failure
    REAL(KIND=dp)                            :: alpha, max_alpha, prefac
    TYPE(cp_dbcsr_type), POINTER             :: comm, comm2, tmp, tmp2
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    unit_nr = cp_logger_get_default_io_unit(logger)

    failure=.FALSE.
    NULLIFY(tmp)
    ALLOCATE(tmp,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(tmp,error=error)
    CALL cp_dbcsr_create(tmp,template=propagator_re,error=error)
    NULLIFY(tmp2)
    ALLOCATE(tmp2,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(tmp2,error=error)
    CALL cp_dbcsr_create(tmp2,template=propagator_re,error=error)
    NULLIFY(comm)
    ALLOCATE(comm,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(comm,error=error)
    CALL cp_dbcsr_create(comm,template=propagator_re,error=error)
    NULLIFY(comm2)
    ALLOCATE(comm2,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(comm2,error=error)
    CALL cp_dbcsr_create(comm2,template=propagator_re,error=error)

    CALL cp_dbcsr_copy(tmp,density_re,error=error)
    CALL cp_dbcsr_copy(tmp2,density_im,error=error)

    convergence=.FALSE.

    DO i=1,20
       prefac=one/i
       CALL cp_complex_dbcsr_gemm_3("N", "N", prefac, propagator_re, propagator_im,&
            tmp, tmp2, zero, comm, comm2, filter_eps=filter_eps, error=error)
       CALL cp_dbcsr_transposed(tmp,comm,error=error)
       CALL cp_dbcsr_transposed(tmp2,comm2,error=error)
       CALL cp_dbcsr_add(comm,tmp,one,one,error=error)
       CALL cp_dbcsr_add(comm2,tmp2,one,-one,error=error)
       CALL cp_dbcsr_add(density_re,comm,one,one,error=error)
       CALL cp_dbcsr_add(density_im,comm2,one,one,error=error)
       CALL cp_dbcsr_copy(tmp,comm,error=error)
       CALL cp_dbcsr_copy(tmp2,comm2,error=error)
       !check for convergence
       max_alpha=zero
       alpha = cp_dbcsr_frobenius_norm(comm)
       IF(alpha>max_alpha) max_alpha=alpha
       alpha = cp_dbcsr_frobenius_norm(comm2)
       IF(alpha>max_alpha) max_alpha=alpha
       IF(max_alpha<eps_exp) convergence = .TRUE.
       IF(convergence) THEN
          IF(unit_nr>0) WRITE (UNIT=unit_nr,FMT="((T3,A,I2,A))")&
                "BCH converged after ",i," steps"
          EXIT
       ENDIF
    ENDDO

    CALL cp_dbcsr_filter(density_re,filter_eps,error=error)
    CALL cp_dbcsr_filter(density_im,filter_eps,error=error)

    CALL cp_assert(convergence,cp_warning_level,cp_assertion_failed,routineP,&
         "BCH method did not converge "//&
         CPSourceFileRef,&
         only_ionode=.TRUE.)

    CALL cp_dbcsr_deallocate_matrix(tmp,error=error)
    CALL cp_dbcsr_deallocate_matrix(tmp2,error=error)
    CALL cp_dbcsr_deallocate_matrix(comm,error=error)
    CALL cp_dbcsr_deallocate_matrix(comm2,error=error)

    CALL timestop(handle)

   END SUBROUTINE

END MODULE rt_matrix_exp_dbcsr
