!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for calculating a complex matrix exponential.
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE rt_matrix_exp

  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_add,&
                                             cp_cfm_gemm
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_set_all,&
                                             cp_cfm_to_cfm,&
                                             cp_cfm_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale,&
                                             cp_fm_scale_and_add
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_maxabsrownorm,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE kinds,                           ONLY: dp
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
#include "cp_common_uses.h"

  IMPLICIT NONE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_matrix_exp'

  PUBLIC :: exp_taylor,&
            exp_diag,&
            backtransform_matrix

  INTERFACE exp_taylor
    MODULE PROCEDURE taylor_only_imaginary,taylor_full_complex
  END INTERFACE


CONTAINS

! *****************************************************************************
!> \brief specialized subroutine for purely imaginary matrix exponentials
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE taylor_only_imaginary(exp_H,im_matrix,eps_exp,error)
    TYPE(cp_fm_p_type), DIMENSION(2)         :: exp_H
    TYPE(cp_fm_type), POINTER                :: im_matrix
    INTEGER                                  :: eps_exp
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: i, ndim, nloop, nsquare, &
                                                ntaylor
    REAL(KIND=dp)                            :: square_fac, Tfac, tmp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data_im
    TYPE(cp_fm_type), POINTER                :: T1, T2, T3, Tres_im, Tres_re

    CALL get_nsquare_ntaylor(im_part=im_matrix,nsquare=nsquare,ntaylor=ntaylor,&
                             eps_exp=eps_exp,error=error)

    CALL cp_fm_get_info(im_matrix,local_data=local_data_im,error=error)
    ndim=im_matrix%matrix_struct%nrow_global   

    square_fac=1.0_dp/(2.0_dp**REAL(nsquare,dp))
    CALL cp_fm_scale(square_fac,im_matrix,error)
    CALL cp_fm_create(T1,&
                       matrix_struct=im_matrix%matrix_struct,&
                       name="T1",&
                       error=error)    

    CALL cp_fm_create(T2,&
                      matrix_struct=T1%matrix_struct,&
                      name="T2",&
                      error=error)
    CALL cp_fm_create(Tres_im,&
                      matrix_struct=T1%matrix_struct,&
                      name="T3",&
                      error=error)
    CALL cp_fm_create(Tres_re,&
                      matrix_struct=T1%matrix_struct,&
                      name="Tres",&
                      error=error)
    tmp=1.0_dp

    CALL cp_fm_set_all(Tres_re,zero,one,error)
    CALL cp_fm_set_all(Tres_im,zero,zero,error)
    CALL cp_fm_set_all(T1,zero,one,error)

    Tfac=one
    nloop=CEILING(REAL(ntaylor,dp)/2.0_dp)

    DO i=1,nloop
       tmp=tmp*(REAL(i,dp)*2.0_dp-1.0_dp)
       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,im_matrix,T1,zero,&
            T2,error)
       Tfac=1._dp/tmp
       IF(MOD(i,2)==0) Tfac=-Tfac
       CALL cp_fm_scale_and_add(one,Tres_im,Tfac,T2,error)
       tmp=tmp*REAL(i,dp)*2.0_dp
       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,im_matrix,T2,zero,&
            T1,error)
       Tfac=1._dp/tmp
       IF(MOD(i,2)==1) Tfac=-Tfac
       CALL cp_fm_scale_and_add(one,Tres_re,Tfac,T1,error)

    END DO


    DO i=1,nsquare   
       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,Tres_re,Tres_re,zero,&
             exp_H(1)%matrix,error)
       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,-one,Tres_im,Tres_im,one,&
             exp_H(1)%matrix,error)

       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,Tres_re,Tres_im,zero,&
             exp_H(2)%matrix,error)
       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,Tres_im,Tres_re,one,&
             exp_H(2)%matrix,error)

       CALL cp_fm_to_fm(exp_H(1)%matrix,Tres_re,error)
       CALL cp_fm_to_fm(exp_H(2)%matrix,Tres_im,error)      
    END DO
    
    CALL cp_fm_release(T1,error)
    CALL cp_fm_release(T2,error)
    CALL cp_fm_release(Tres_re,error)
    CALL cp_fm_release(Tres_im,error)


 
  END SUBROUTINE taylor_only_imaginary


! *****************************************************************************
!> \brief subroutine for general complex matrix exponentials
!>        on input a separate cp_fm_type for real and complex part 
!>        on output a size 2 cp_fm_p_type, frst element is the real part of
!>        the exponential second the imaginary
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE taylor_full_complex(exp_H,re_part,im_part,eps_exp,error)
    TYPE(cp_fm_p_type), DIMENSION(2)         :: exp_H
    TYPE(cp_fm_type), POINTER                :: re_part, im_part
    INTEGER                                  :: eps_exp
    TYPE(cp_error_type), INTENT(inout)       :: error

    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp)                         :: Tfac
    INTEGER                                  :: i, ndim, nsquare, ntaylor
    REAL(KIND=dp)                            :: square_fac, tmp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data_im, local_data_re
    TYPE(cp_cfm_type), POINTER               :: T1, T2, T3, Tres

    CALL get_nsquare_ntaylor(re_part,im_part,nsquare,ntaylor,eps_exp,error)

    CALL cp_fm_get_info(re_part,local_data=local_data_re,error=error)
    CALL cp_fm_get_info(im_part,local_data=local_data_im,error=error)
    ndim=re_part%matrix_struct%nrow_global

    CALL cp_cfm_create(T1,&
                       matrix_struct=re_part%matrix_struct,&
                       name="T1",&
                       error=error)   
    
    square_fac=2.0_dp**REAL(nsquare,dp)
    
    T1%local_data=CMPLX(local_data_re/square_fac,local_data_im/square_fac,KIND=dp)

    CALL cp_cfm_create(T2,&
                      matrix_struct=T1%matrix_struct,&
                      name="T2",&
                      error=error)
    CALL cp_cfm_create(T3,&
                      matrix_struct=T1%matrix_struct,&
                      name="T3",&
                      error=error)
    CALL cp_cfm_create(Tres,&
                      matrix_struct=T1%matrix_struct,&
                      name="Tres",&
                      error=error)
    tmp=1.0_dp
    CALL cp_cfm_set_all(Tres,zero,one,error)
    CALL cp_cfm_set_all(T2,zero,one,error)
    Tfac=one
    DO i=1,ntaylor
       tmp=tmp*REAL(i,dp)
       CALL cp_cfm_gemm("N","N",ndim,ndim,ndim,one,T1,T2,zero,&
            T3,error)
       Tfac=CMPLX(1._dp/tmp,0.0_dp)
       CALL cp_cfm_add(one,Tres,Tfac,T3,error)
       CALL cp_cfm_to_cfm(T3,T2,error)
    END DO


    DO i=1,nsquare   
       CALL cp_cfm_gemm("N","N",ndim,ndim,ndim,one,Tres,Tres,zero,&
            T2,error)
       CALL cp_cfm_to_cfm(T2,Tres,error)
    END DO
    

    exp_H(1)%matrix%local_data=REAL(Tres%local_data,KIND=dp)
    exp_H(2)%matrix%local_data=AIMAG(Tres%local_data)

    CALL cp_cfm_release(T1,error)
    CALL cp_cfm_release(T2,error)
    CALL cp_cfm_release(T3,error)
    CALL cp_cfm_release(Tres,error)
 
  END SUBROUTINE taylor_full_complex


! *****************************************************************************
!> \brief exponential of a purely imaginary matrix,
!>        calculated using diagonalization
!> \author Florian Schiffmann (02.09)
! *****************************************************************************
  SUBROUTINE get_nsquare_ntaylor(re_part,im_part,nsquare,ntaylor,eps_exp,error)
    
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: re_part
    TYPE(cp_fm_type), POINTER                :: im_part
    INTEGER                                  :: nsquare, ntaylor, eps_exp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_nsquare_ntaylor', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol_local, ind_n, &
                                                irow_local, ncol_local, &
                                                nrow_local
    INTEGER, DIMENSION(3, 6)                 :: squares, taylors
    REAL(dp)                                 :: lognorm, norm
    TYPE(cp_fm_type), POINTER                :: tmp_fm

    squares(1,:)=(/1,2,3,5,9,13/)
    squares(2,:)=(/1,4,3,7,9,13/)
    squares(3,:)=(/1,4,5,7,10,14/)
    taylors(1,:)=(/3,4,6,8,7,6/)
    taylors(2,:)=(/4,4,8,7,9,8/)
    taylors(3,:)=(/5,5,7,9,10,8/)

    
    IF(PRESENT(re_part))THEN
       CALL cp_fm_create(tmp_fm,&
            matrix_struct=im_part%matrix_struct,&
            name="tmp_fm",&
            error=error)
       CALL cp_fm_get_info(tmp_fm, nrow_local=nrow_local, ncol_local=ncol_local,error=error)
       
       tmp_fm%local_data=0.0_dp
       
       DO icol_local=1,ncol_local
          DO irow_local=1,nrow_local
             tmp_fm%local_data(irow_local,icol_local) = &
                  im_part%local_data(irow_local,icol_local)*im_part%local_data(irow_local,icol_local)
          END DO
       END DO

       DO icol_local=1,ncol_local
          DO irow_local=1,nrow_local
             tmp_fm%local_data(irow_local,icol_local) = tmp_fm%local_data(irow_local,icol_local)+ &
                  re_part%local_data(irow_local,icol_local)*re_part%local_data(irow_local,icol_local)
          END DO
       END DO
       tmp_fm%local_data=SQRT(tmp_fm%local_data)
       CALL  cp_fm_maxabsrownorm(tmp_fm, norm,error)
       CALL cp_fm_release(tmp_fm,error)
    ELSE
       CALL  cp_fm_maxabsrownorm(im_part, norm,error)    
    END IF
    
    lognorm=LOG10(norm)
    
    ind_n=CEILING(lognorm)+3
    IF(ind_n.LT.1)ind_n=1
    IF(ind_n.GT.6)THEN
       ntaylor=12
       nsquare=16
       CALL cp_assert(.TRUE.,cp_warning_level,cp_assertion_failed,routineP,&
               "Norm of the matrix exceeds the parametrized range, no warranty that the exponential"//&
               "has the specified accuracy"//&
               CPSourceFileRef,&
               only_ionode=.TRUE.)
    ELSE
       ntaylor=taylors(eps_exp,ind_n)
       nsquare=squares(eps_exp,ind_n)
    END IF
    
  END SUBROUTINE get_nsquare_ntaylor


! *****************************************************************************
!> \brief exponential of a purely imaginary matrix,
!>        calculated using diagonalization
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE exp_diag(rtp,exp_H,H_fm,error)

    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(cp_fm_p_type), DIMENSION(2)         :: exp_H
    TYPE(cp_fm_type), POINTER                :: H_fm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'exp_diag', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: ndim, stat
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eigval_H
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: exp_eig_val
    TYPE(cp_fm_type), POINTER                :: eigvec_H, S_half, &
                                                S_minus_half, tmp_mat

    CALL get_rtp(rtp=rtp,S_half=S_half,S_minus_half=S_minus_half,error=error)

    CALL cp_fm_create(eigvec_H,&
         matrix_struct=exp_H(1)%matrix%matrix_struct,&
         name="tmp_EVEC",&
         error=error)
    CALL cp_fm_create(tmp_mat,&
         matrix_struct=exp_H(1)%matrix%matrix_struct,&
         name="tmp_mat",&
         error=error)

    !H_fm,eigvec_H,eigval_H is used as tmp, for diagonalizing S

    ndim=H_fm%matrix_struct%nrow_global
    ALLOCATE(eigval_H(ndim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(exp_eig_val(ndim,2),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,H_fm,S_minus_half,zero,&
         tmp_mat,error)
    CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,S_minus_half,tmp_mat,zero,&
         H_fm,error)
          
    CALL cp_fm_syevd(H_fm,eigvec_H,eigval_H,error)

    exp_eig_val(:,1)=COS(eigval_H(:))
    exp_eig_val(:,2)=SIN(eigval_H(:))
          
    CALL backtransform_matrix(exp_eig_val(:,1),eigvec_H,exp_H(1)%matrix,error)       
    
    CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,exp_H(1)%matrix,S_half,zero,&
         tmp_mat,error)
    
    CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,S_minus_half,tmp_mat,zero,&
         exp_H(1)%matrix,error)
    
    
    CALL backtransform_matrix(exp_eig_val(:,2),eigvec_H,exp_H(2)%matrix,error)
    
    CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,exp_H(2)%matrix,S_half,zero,&
         tmp_mat,error)
    
    CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,S_minus_half,tmp_mat,zero,&
         exp_H(2)%matrix,error)


  END SUBROUTINE exp_diag


  SUBROUTINE backtransform_matrix(Eval,eigenvec,matrix,error)

    REAL(dp), DIMENSION(:)                   :: Eval
    TYPE(cp_fm_type), POINTER                :: eigenvec, matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: i, j, l, m, ncol_local, ndim, &
                                                nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    TYPE(cp_fm_type), POINTER                :: tmp

    CALL cp_fm_create(tmp,&
         matrix_struct=matrix%matrix_struct,&
         name="TMP_BT",&
         error=error)
    CALL cp_fm_get_info(matrix, nrow_local=nrow_local, ncol_local=ncol_local,&
         row_indices=row_indices,col_indices=col_indices,error=error)

    ndim=matrix%matrix_struct%nrow_global

    CALL cp_fm_set_all(tmp,zero,zero,error)
    DO i=1,ncol_local
       l=col_indices(i)
       DO j=1,nrow_local
          m=row_indices(j)
          tmp%local_data(m,l)=eigenvec%local_data(m,l)*Eval(l)
       END DO
    END DO
    CALL cp_fm_gemm("N","T",ndim,ndim,ndim,one,tmp,eigenvec,zero,&
         matrix ,error)

    CALL cp_fm_release(tmp,error)

  END SUBROUTINE backtransform_matrix
    
END MODULE rt_matrix_exp
