!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for calculating a complex matrix exponential.
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE rt_matrix_exp

  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_add,&
                                             cp_cfm_gemm,&
                                             cp_cfm_solve
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_p_type,&
                                             cp_cfm_release,&
                                             cp_cfm_set_all,&
                                             cp_cfm_to_cfm,&
                                             cp_cfm_type
  USE cp_fm_basic_linalg,              ONLY: cp_complex_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_solve
  USE cp_fm_struct,                    ONLY: cp_fm_struct_double,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fac
  USE message_passing,                 ONLY: mp_comm_free,&
                                             mp_max,&
                                             mp_reordering,&
                                             mp_sum
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_matrix_exp'

  PUBLIC :: taylor_only_imaginary,&
            backtransform_matrix,&
            taylor_full_complex,&
            exp_pade_full_complex,&
            exp_pade_only_imaginary,&
            get_nsquare_norder,&
            arnoldi, exp_pade_real

CONTAINS

  ! *****************************************************************************
  !> \brief Is still left from diagonalization, should be removed later but is
  !>  still needed for the guess for the padé/Taylor method
  !> \author Florian Schiffmann (02.09)
  ! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param Eval ...
!> \param eigenvec ...
!> \param matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE backtransform_matrix(Eval,eigenvec,matrix,error)

    REAL(dp), DIMENSION(:)                   :: Eval
    TYPE(cp_fm_type), POINTER                :: eigenvec, matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'backtransform_matrix', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, j, l, ncol_local, &
                                                ndim, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    TYPE(cp_fm_type), POINTER                :: tmp

    CALL timeset(routineN,handle)
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
          tmp%local_data(j,i)=eigenvec%local_data(j,i)*Eval(l)
       END DO
    END DO
    CALL cp_gemm("N","T",ndim,ndim,ndim,one,tmp,eigenvec,zero,&
         matrix ,error)

    CALL cp_fm_release(tmp,error)
    CALL timestop(handle)

  END SUBROUTINE backtransform_matrix

! *****************************************************************************
!> \brief specialized subroutine for purely imaginary matrix exponentials
!> \param exp_H ...
!> \param im_matrix ...
!> \param nsquare ...
!> \param ntaylor ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE taylor_only_imaginary(exp_H,im_matrix,nsquare,ntaylor,error)
    TYPE(cp_fm_p_type), DIMENSION(2)         :: exp_H
    TYPE(cp_fm_type), POINTER                :: im_matrix
    INTEGER                                  :: nsquare, ntaylor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'taylor_only_imaginary', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, ndim, nloop
    REAL(KIND=dp)                            :: square_fac, Tfac, tmp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data_im
    TYPE(cp_fm_type), POINTER                :: T1, T2, Tres_im, Tres_re

    CALL timeset(routineN,handle)

    CALL cp_fm_get_info(im_matrix,local_data=local_data_im,error=error)
    ndim=im_matrix%matrix_struct%nrow_global

    square_fac=1.0_dp/(2.0_dp**REAL(nsquare,dp))
!    CALL cp_fm_scale(square_fac,im_matrix,error)
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
       CALL cp_gemm("N","N",ndim,ndim,ndim,square_fac,im_matrix,T1,zero,&
!       CALL cp_gemm("N","N",ndim,ndim,ndim,one,im_matrix,T1,zero,&
            T2,error)
       Tfac=1._dp/tmp
       IF(MOD(i,2)==0) Tfac=-Tfac
       CALL cp_fm_scale_and_add(one,Tres_im,Tfac,T2,error)
       tmp=tmp*REAL(i,dp)*2.0_dp
       CALL cp_gemm("N","N",ndim,ndim,ndim,square_fac,im_matrix,T2,zero,&
!       CALL cp_gemm("N","N",ndim,ndim,ndim,one,im_matrix,T2,zero,&
            T1,error)
       Tfac=1._dp/tmp
       IF(MOD(i,2)==1) Tfac=-Tfac
       CALL cp_fm_scale_and_add(one,Tres_re,Tfac,T1,error)

    END DO

    IF(nsquare.GT.0)THEN
       DO i=1,nsquare
          CALL cp_complex_fm_gemm("N","N",ndim,ndim,ndim,one,Tres_re,Tres_im,&
               Tres_re,Tres_im,zero,exp_H(1)%matrix,exp_H(2)%matrix,error)

          CALL cp_fm_to_fm(exp_H(1)%matrix,Tres_re,error)
          CALL cp_fm_to_fm(exp_H(2)%matrix,Tres_im,error)
       END DO
    ELSE
       CALL cp_fm_to_fm(Tres_re,exp_H(1)%matrix,error)
       CALL cp_fm_to_fm(Tres_im,exp_H(2)%matrix,error)
    END IF

    CALL cp_fm_release(T1,error)
    CALL cp_fm_release(T2,error)
    CALL cp_fm_release(Tres_re,error)
    CALL cp_fm_release(Tres_im,error)

    CALL timestop(handle)


  END SUBROUTINE taylor_only_imaginary


! *****************************************************************************
!> \brief subroutine for general complex matrix exponentials
!>        on input a separate cp_fm_type for real and complex part
!>        on output a size 2 cp_fm_p_type, first element is the real part of
!>        the exponential second the imaginary
!> \param exp_H ...
!> \param re_part ...
!> \param im_part ...
!> \param nsquare ...
!> \param ntaylor ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE taylor_full_complex(exp_H,re_part,im_part,nsquare,ntaylor,error)
    TYPE(cp_fm_p_type), DIMENSION(2)         :: exp_H
    TYPE(cp_fm_type), POINTER                :: re_part, im_part
    INTEGER                                  :: nsquare, ntaylor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'taylor_full_complex', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp)                         :: Tfac
    INTEGER                                  :: handle, i, ndim
    REAL(KIND=dp)                            :: square_fac, tmp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data_im, local_data_re
    TYPE(cp_cfm_type), POINTER               :: T1, T2, T3, Tres

    CALL timeset(routineN,handle)
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
       Tfac=CMPLX(1._dp/tmp,0.0_dp,kind=dp)
       CALL cp_cfm_add(one,Tres,Tfac,T3,error)
       CALL cp_cfm_to_cfm(T3,T2,error)
    END DO

    IF(nsquare.GT.0)THEN
       DO i=1,nsquare
          CALL cp_cfm_gemm("N","N",ndim,ndim,ndim,one,Tres,Tres,zero,&
               T2,error)
          CALL cp_cfm_to_cfm(T2,Tres,error)
       END DO
    END IF


    exp_H(1)%matrix%local_data=REAL(Tres%local_data,KIND=dp)
    exp_H(2)%matrix%local_data=AIMAG(Tres%local_data)

    CALL cp_cfm_release(T1,error)
    CALL cp_cfm_release(T2,error)
    CALL cp_cfm_release(T3,error)
    CALL cp_cfm_release(Tres,error)
    CALL timestop(handle)

  END SUBROUTINE taylor_full_complex


! *****************************************************************************
!> \brief optimization function for pade/taylor order and number of squaring steps
!> \param norm ...
!> \param nsquare ...
!> \param norder ...
!> \param eps_exp ...
!> \param method ...
!> \param do_emd ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************
  SUBROUTINE get_nsquare_norder(norm,nsquare,norder,eps_exp,method,do_emd,error)

    REAL(dp)                                 :: norm
    INTEGER                                  :: nsquare, norder
    REAL(dp)                                 :: eps_exp
    INTEGER                                  :: method
    LOGICAL                                  :: do_emd
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_nsquare_norder', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cost, i, iscale, orders(3), &
                                                p, prev_cost, q
    LOGICAL                                  :: new_scale
    REAL(dp)                                 :: D, eval, myval, N, scaleD, &
                                                scaleN

    orders(:)=(/12,12,12/)
    IF(method==2)THEN
       DO  iscale=0,12
          new_scale=.FALSE.
          eval=norm/(2.0_dp**REAL(iscale,dp))
          DO q=1,12
             DO p=MAX(1,q-1),q
                IF(p>q)EXIT
                D=1.0_dp
                N=1.0_dp
                DO i=1,q
                   IF(i.LE.p)scaleN=fac(p+q-i)*fac(p)/(fac(p+q)*fac(i)*fac(p-i))
                   scaleD=(-1.0)**i * fac(p+q-i)*fac(q)/(fac(p+q)*fac(i)*fac(q-i))
                   IF(i.LE.p)N=N+scaleN*eval**i
                   D=D+scaleD*eval**i
                END DO
                IF(ABS((EXP(norm)-(N/D)**(2.0_dp**iscale))/MAX(1.0_dp,EXP(norm))).LE.eps_exp)THEN
                   IF(do_emd)THEN
                      cost=iscale+q
                      prev_cost=orders(1)+orders(2)
                   ELSE
                      cost=iscale+CEILING(REAL(q,dp)/3.0_dp)
                      prev_cost=orders(1)+CEILING(REAL(orders(2),dp)/3.0_dp)
                   END IF
                   IF(cost.LT.prev_cost)THEN
                      orders(:)=(/iscale,q,p/)
                      myval=(N/D)**(2.0_dp**iscale)
                   END IF
                   new_scale=.TRUE.
                   EXIT
                END IF
             END DO
             IF(new_scale)EXIT
          END DO
          IF(iscale.GE.orders(1)+orders(2).AND.new_scale)EXIT
       END DO
    ELSE IF(method==1)THEN
       q=0
       eval=norm
       DO  iscale=0,6
          new_scale=.FALSE.
          IF(iscale.GE.1)eval=norm/(2.0_dp**REAL(iscale,dp))
          DO p=1,20
             D=1.0_dp
             N=1.0_dp
             DO i=1,p
                scaleN=1.0_dp/fac(i)
                N=N+scaleN*(eval**REAL(i,dp))
             END DO
             IF(ABS((EXP(norm)-N**(2.0_dp**REAL(iscale,dp)))/MAX(1.0_dp,EXP(norm))).LE.eps_exp)THEN
                IF(do_emd)THEN
                   cost=iscale+p
                   prev_cost=orders(1)+orders(2)
                ELSE
                   cost=iscale+CEILING(REAL(p,dp)/3.0_dp)
                   prev_cost=orders(1)+CEILING(REAL(orders(2),dp)/3.0_dp)
                END IF
                IF(cost.LT.prev_cost)THEN
                   orders(:)=(/iscale,p,0/)
                   myval=(N)**(2.0_dp**iscale)
                END IF
                new_scale=.TRUE.
                EXIT
             END IF
          END DO
          IF(iscale.GE.orders(1)+orders(2).AND.new_scale)EXIT
       END DO
    END IF

    nsquare=orders(1)
    norder=orders(2)

  END SUBROUTINE get_nsquare_norder

! *****************************************************************************
!> \brief exponential of a complex matrix,
!>        calculated using pade approximation together with scaling and squaring
!> \param exp_H ...
!> \param re_part ...
!> \param im_part ...
!> \param nsquare ...
!> \param npade ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE exp_pade_full_complex(exp_H,re_part,im_part,nsquare,npade,error)
    TYPE(cp_fm_p_type), DIMENSION(2)         :: exp_H
    TYPE(cp_fm_type), POINTER                :: re_part, im_part
    INTEGER                                  :: nsquare, npade
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'exp_pade_full_complex', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp)                         :: scaleD, scaleN
    INTEGER                                  :: handle, i, ldim, ndim, p, q
    REAL(KIND=dp)                            :: square_fac, tmp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data_im, local_data_re
    TYPE(cp_cfm_p_type), DIMENSION(:), &
      POINTER                                :: mult_p
    TYPE(cp_cfm_type), POINTER               :: Dpq, fin_p, Npq, T1, T2, Tres

    p=npade
    q=npade !p==q seems to be neccesary for the rest of the code

    CALL timeset(routineN,handle)
    CALL cp_fm_get_info(re_part,local_data=local_data_re,ncol_local=ldim,&
         nrow_global=ndim,error=error)
    CALL cp_fm_get_info(im_part,local_data=local_data_im,error=error)

    ALLOCATE(mult_p(2))

    CALL cp_cfm_create(Dpq,&
                       matrix_struct=re_part%matrix_struct,&
                       name="Dpq",&
                       error=error)

    square_fac=2.0_dp**REAL(nsquare,dp)

    CALL cp_cfm_create(T1,&
                      matrix_struct=Dpq%matrix_struct,&
                      name="T1",&
                      error=error)

    CALL cp_cfm_create(T2,&
                      matrix_struct=T1%matrix_struct,&
                      name="T2",&
                      error=error)
    CALL cp_cfm_create(Npq,&
                      matrix_struct=T1%matrix_struct,&
                      name="Npq",&
                      error=error)
    CALL cp_cfm_create(Tres,&
                      matrix_struct=T1%matrix_struct,&
                      name="Tres",&
                      error=error)

    DO i=1,ldim
       T2%local_data(:,i)=CMPLX(local_data_re(:,i)/square_fac,local_data_im(:,i)/square_fac,KIND=dp)
    END DO
    CALL cp_cfm_to_cfm(T2,T1,error)
    mult_p(1)%matrix=>T2
    mult_p(2)%matrix=>Tres
    tmp=1.0_dp
    CALL cp_cfm_set_all(Npq,zero,one,error)
    CALL cp_cfm_set_all(Dpq,zero,one,error)

    CALL cp_cfm_add(one,Npq,one*0.5_dp,T2,error)
    CALL cp_cfm_add(one,Dpq,-one*0.5_dp,T2,error)

    npade=MAX(p,q)
    IF(npade.GT.2)THEN
       DO i=2,npade
          IF(i.LE.p)scaleN=CMPLX(fac(p +q -i )*fac(p )/(fac(p +q )*fac(i )*fac(p -i )),0.0_dp,kind=dp)
          scaleD=CMPLX((-1.0_dp)**i  * fac(p +q -i )*fac(q)/(fac(p +q )*fac(i )*fac(q -i)),0.0_dp,kind=dp)
          CALL cp_cfm_gemm("N","N",ndim,ndim,ndim,one,T1,mult_p(MOD(i,2)+1)%matrix,zero,&
               mult_p(MOD(i+1,2)+1)%matrix,error)
          IF(i.LE.p)CALL cp_cfm_add(one,Npq,scaleN, mult_p(MOD(i+1,2)+1)%matrix,error)
          IF(i.LE.q)CALL cp_cfm_add(one,Dpq,scaleD, mult_p(MOD(i+1,2)+1)%matrix,error)
       END DO
    END IF

    CALL cp_cfm_solve(Dpq,Npq,error=error)

    mult_p(2)%matrix=>Npq
    mult_p(1)%matrix=>Tres
    IF(nsquare.GT.0)THEN
       DO i=1,nsquare
          CALL cp_cfm_gemm("N","N",ndim,ndim,ndim,one,mult_p(MOD(i,2)+1)%matrix,mult_p(MOD(i,2)+1)%matrix,zero,&
               mult_p(MOD(i+1,2)+1)%matrix,error)
          fin_p=> mult_p(MOD(i+1,2)+1)%matrix
       END DO
    ELSE
       fin_p=>Npq
    END IF
    DO i=1,ldim
       exp_H(1)%matrix%local_data(:,i)=REAL(fin_p%local_data(:,i),KIND=dp)
       exp_H(2)%matrix%local_data(:,i)=AIMAG(fin_p%local_data(:,i))
    END DO

    CALL cp_cfm_release(Npq,error)
    CALL cp_cfm_release(Dpq,error)
    CALL cp_cfm_release(T1,error)
    CALL cp_cfm_release(T2,error)
    CALL cp_cfm_release(Tres,error)
    DEALLOCATE(mult_p)
    CALL timestop(handle)

  END SUBROUTINE exp_pade_full_complex

! *****************************************************************************
!> \brief exponential of a complex matrix,
!>        calculated using pade approximation together with scaling and squaring
!> \param exp_H ...
!> \param im_part ...
!> \param nsquare ...
!> \param npade ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE exp_pade_only_imaginary(exp_H,im_part,nsquare,npade,error)
    TYPE(cp_fm_p_type), DIMENSION(2)         :: exp_H
    TYPE(cp_fm_type), POINTER                :: im_part
    INTEGER                                  :: nsquare, npade
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'exp_pade_only_imaginary', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)
    REAL(KIND=dp), PARAMETER                 :: rone = (1.0_dp,0.0_dp) , &
                                                rzero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp)                         :: scaleD, scaleN
    INTEGER                                  :: handle, i, j, k, ldim, ndim, &
                                                p, q
    REAL(KIND=dp)                            :: my_fac, square_fac
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data_im
    TYPE(cp_cfm_p_type), DIMENSION(:), &
      POINTER                                :: cmult_p
    TYPE(cp_cfm_type), POINTER               :: Dpq, fin_p, Npq, T1
    TYPE(cp_fm_type), POINTER                :: T2, Tres

    CALL timeset(routineN,handle)
    p=npade 
    q=npade !p==q seems to be neccesary for the rest of the code

    CALL cp_fm_get_info(im_part,local_data=local_data_im,ncol_local=ldim,nrow_global=ndim,error=error)
    square_fac=1/(2.0_dp**REAL(nsquare,dp))
    ALLOCATE(cmult_p(2))

    CALL cp_cfm_create(Dpq,&
                       matrix_struct=im_part%matrix_struct,&
                       name="Dpq",&
                       error=error)

    CALL cp_cfm_create(Npq,&
                      matrix_struct=Dpq%matrix_struct,&
                      name="Npq",&
                      error=error)

    CALL cp_cfm_create(T1,&
                      matrix_struct=Dpq%matrix_struct,&
                      name="T1",&
                      error=error)

    CALL cp_fm_create(T2,&
                      matrix_struct=T1%matrix_struct,&
                      name="T2",&
                      error=error)

    CALL cp_fm_create(Tres,&
                      matrix_struct=T1%matrix_struct,&
                      name="Tres",&
                      error=error)

!    DO i=1,ldim
!       local_data_im(:,i)=local_data_im(:,i)/square_fac
!    END DO

    CALL cp_fm_to_fm(im_part,T2,error)

    CALL cp_cfm_set_all(Npq,zero,one,error)
    CALL cp_cfm_set_all(Dpq,zero,one,error)

    DO i=1,ldim
       Npq%local_data(:,i)= Npq%local_data(:,i)+0.5_dp*square_fac*CMPLX(0.0_dp,local_data_im(:,i),dp)
       Dpq%local_data(:,i)= Dpq%local_data(:,i)-0.5_dp*square_fac*CMPLX(0.0_dp,local_data_im(:,i),dp)
    END DO




    IF(npade.GT.2)THEN
       DO j=1,FLOOR(npade/2.0_dp)
          i=2*j
          my_fac=(-1.0_dp)**j
          IF(i.LE.p)scaleN=my_fac*fac(p +q -i )*fac(p )/(fac(p +q )*fac(i )*fac(p -i ))
          scaleD=my_fac* fac(p +q -i )*fac(q)/(fac(p +q )*fac(i )*fac(q -i))
          CALL cp_gemm("N","N",ndim,ndim,ndim,square_fac,im_part,T2,rzero,Tres,error)

          DO k=1,ldim
             Npq%local_data(:,k)= Npq%local_data(:,k)+scaleN*CMPLX( Tres%local_data(:,k),0.0_dp,dp)
             Dpq%local_data(:,k)= Dpq%local_data(:,k)+scaleD*CMPLX( Tres%local_data(:,k),0.0_dp,dp)
          END DO

          IF(2*j+1.le.q)THEN
          i=2*j+1
          IF(i.LE.p)scaleN=my_fac*fac(p +q -i )*fac(p )/(fac(p +q )*fac(i )*fac(p -i ))
          scaleD=-my_fac*fac(p +q -i )*fac(q)/(fac(p +q )*fac(i )*fac(q -i))
          CALL cp_gemm("N","N",ndim,ndim,ndim,square_fac,im_part,Tres,rzero,T2,error)

          DO k=1,ldim
             Npq%local_data(:,k)= Npq%local_data(:,k)+scaleN*CMPLX(0.0_dp, T2%local_data(:,k),dp)
             Dpq%local_data(:,k)= Dpq%local_data(:,k)+scaleD*CMPLX(0.0_dp, T2%local_data(:,k),dp)
          END DO
          ENDIF
       END DO
    END IF

    CALL cp_cfm_solve(Dpq,Npq,error=error)

    cmult_p(2)%matrix=>Npq
    cmult_p(1)%matrix=>T1
    IF(nsquare.GT.0)THEN
       DO i=1,nsquare
          CALL cp_cfm_gemm("N","N",ndim,ndim,ndim,one,cmult_p(MOD(i,2)+1)%matrix,cmult_p(MOD(i,2)+1)%matrix,zero,&
               cmult_p(MOD(i+1,2)+1)%matrix,error)
          fin_p=> cmult_p(MOD(i+1,2)+1)%matrix
       END DO
    ELSE
       fin_p=>Npq
    END IF

    DO k=1,ldim
       exp_H(1)%matrix%local_data(:,k)=REAL(fin_p%local_data(:,k),KIND=dp)
       exp_H(2)%matrix%local_data(:,k)=AIMAG(fin_p%local_data(:,k))
    END DO

    CALL cp_cfm_release(Npq,error)
    CALL cp_cfm_release(Dpq,error)
    CALL cp_cfm_release(T1,error)
    CALL cp_fm_release(T2,error)
    CALL cp_fm_release(Tres,error)
    DEALLOCATE(cmult_p)
    CALL timestop(handle)
  END SUBROUTINE exp_pade_only_imaginary

 ! *****************************************************************************
!> \brief exponential of a real matrix,
!>        calculated using pade approximation together with scaling and squaring
!> \param exp_H ...
!> \param matrix ...
!> \param nsquare ...
!> \param npade ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE exp_pade_real(exp_H,matrix,nsquare,npade,error)
    TYPE(cp_fm_type)                         :: exp_H
    TYPE(cp_fm_type), POINTER                :: matrix
    INTEGER                                  :: nsquare, npade
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'exp_pade_real', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    COMPLEX(KIND=dp)                         :: scaleD, scaleN
    INTEGER                                  :: handle, i, j, k, ldim, ndim, &
                                                p, q
    REAL(KIND=dp)                            :: my_fac, square_fac
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mult_p
    TYPE(cp_fm_type), POINTER                :: Dpq, fin_p, Npq, T1, T2, Tres

    CALL timeset(routineN,handle)
    p=npade
    q=npade !p==q seems to be neccesary for the rest of the code

    CALL cp_fm_get_info(matrix,local_data=local_data,ncol_local=ldim,nrow_global=ndim,error=error)
    square_fac=2.0_dp**REAL(nsquare,dp)
    ALLOCATE(mult_p(2))

    CALL cp_fm_create(Dpq,&
                       matrix_struct=matrix%matrix_struct,&
                       name="Dpq",&
                       error=error)

    CALL cp_fm_create(Npq,&
                      matrix_struct=Dpq%matrix_struct,&
                      name="Npq",&
                      error=error)

    CALL cp_fm_create(T1,&
                      matrix_struct=Dpq%matrix_struct,&
                      name="T1",&
                      error=error)

    CALL cp_fm_create(T2,&
                      matrix_struct=T1%matrix_struct,&
                      name="T2",&
                      error=error)

    CALL cp_fm_create(Tres,&
                      matrix_struct=T1%matrix_struct,&
                      name="Tres",&
                      error=error)

    DO i=1,ldim
       T2%local_data(:,i)=local_data(:,i)/square_fac
    END DO

    CALL cp_fm_to_fm(T2,T1,error)
    CALL cp_fm_set_all(Npq,zero,one,error)
    CALL cp_fm_set_all(Dpq,zero,one,error)

    DO i=1,ldim
       Npq%local_data(:,i)= Npq%local_data(:,i)+0.5_dp*local_data(:,i)
       Dpq%local_data(:,i)= Dpq%local_data(:,i)-0.5_dp*local_data(:,i)
    END DO

    mult_p(1)%matrix=>T2
    mult_p(2)%matrix=>Tres

    IF(npade.GE.2)THEN
       DO j=2,npade
          my_fac=(-1.0_dp)**j
          scaleN=fac(p +q -j )*fac(p )/(fac(p +q )*fac(j )*fac(p -j ))
          scaleD=my_fac* fac(p +q -j )*fac(q)/(fac(p +q )*fac(j )*fac(q -j))
          CALL cp_gemm("N","N",ndim,ndim,ndim,one,mult_p(MOD(j,2)+1)%matrix,T1,&
          zero,mult_p(MOD(j+1,2)+1)%matrix,error)

          DO k=1,ldim
             Npq%local_data(:,k)= Npq%local_data(:,k)+scaleN* mult_p(MOD(j+1,2)+1)%matrix%local_data(:,k)
             Dpq%local_data(:,k)= Dpq%local_data(:,k)+scaleD* mult_p(MOD(j+1,2)+1)%matrix%local_data(:,k)
          END DO
       END DO
    END IF

    CALL cp_fm_solve(Dpq,Npq,error)

    mult_p(2)%matrix=>Npq
    mult_p(1)%matrix=>T1
    IF(nsquare.GT.0)THEN
       DO i=1,nsquare
          CALL cp_gemm("N","N",ndim,ndim,ndim,one,mult_p(MOD(i,2)+1)%matrix,mult_p(MOD(i,2)+1)%matrix,zero,&
               mult_p(MOD(i+1,2)+1)%matrix,error)
          fin_p=> mult_p(MOD(i+1,2)+1)%matrix
       END DO
    ELSE
       fin_p=>Npq
    END IF

    DO k=1,ldim
       exp_H%local_data(:,k)=fin_p%local_data(:,k)
    END DO

    CALL cp_fm_release(Npq,error)
    CALL cp_fm_release(Dpq,error)
    CALL cp_fm_release(T1,error)
    CALL cp_fm_release(T2,error)
    CALL cp_fm_release(Tres,error)
    DEALLOCATE(mult_p)
    CALL timestop(handle)

  END SUBROUTINE exp_pade_real

! ******************************************************************************
!> \brief exponential of a complex matrix,
!>        calculated using arnoldi subspace method (directly apllies to the MOs)
!> \param mos_old ...
!> \param mos_new ...
!> \param eps_exp ...
!> \param Hre ...
!> \param Him ...
!> \param mos_next ...
!> \param narn_old ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE arnoldi(mos_old,mos_new,eps_exp,Hre,Him,mos_next,narn_old,error)

    TYPE(cp_fm_p_type), DIMENSION(2)         :: mos_old, mos_new
    REAL(KIND=dp)                            :: eps_exp
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: Hre
    TYPE(cp_fm_type), POINTER                :: Him
    TYPE(cp_fm_p_type), DIMENSION(2), &
      OPTIONAL                               :: mos_next
    INTEGER                                  :: narn_old
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'arnoldi', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)

    INTEGER :: col_group, handle, i, icol_local, idim, info, j, l, mydim, &
      nao, narnoldi, ncol_local, newdim, nmo, npade, pade_step, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipivot
    INTEGER, DIMENSION(:), POINTER           :: col_indices, col_procs
    LOGICAL                                  :: convergence, double_col, &
                                                double_row, failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: last_norm, norm1, results
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: D, mat1, mat2, mat3, N
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: H_approx, H_approx_save
    REAL(KIND=dp)                            :: conv_norm, prefac, scaleD, &
                                                scaleN
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: V_mats
    TYPE(cp_fm_struct_type), POINTER         :: mo_struct, newstruct
    TYPE(cp_para_env_type), POINTER          :: para_env

    failure=.FALSE.
    CALL timeset(routineN,handle)
    para_env=>mos_new(1)%matrix%matrix_struct%para_env

    CALL cp_fm_get_info(mos_old(1)%matrix,ncol_local=ncol_local,col_indices=col_indices,&
                        nrow_global=nao,ncol_global=nmo,matrix_struct=mo_struct,error=error)
    narnoldi=MIN(18,nao)

    ALLOCATE(results(ncol_local),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(norm1(ncol_local),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(V_mats(narnoldi+1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(last_norm(ncol_local),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(H_approx(narnoldi,narnoldi,ncol_local),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(H_approx_save(narnoldi,narnoldi,ncol_local),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    col_procs=>mo_struct%context%blacs2mpi(:,mo_struct%context%mepos(2))
    CALL mp_reordering(para_env%group,col_group,col_procs)
    

    double_col=.TRUE.
    double_row=.FALSE.
    CALL cp_fm_struct_double(newstruct,mo_struct,mo_struct%context,double_col,double_row,error)
    H_approx_save=zero

    DO i=1,narnoldi+1
       CALL cp_fm_create(V_mats(i)%matrix, matrix_struct=newstruct,&
                          name="V_mat"//cp_to_string(i),error=error)
    END DO
    CALL cp_fm_get_info(V_mats(1)%matrix,ncol_global=newdim,error=error)

    norm1=0.0_dp
!$OMP PARALLEL DO PRIVATE(icol_local) DEFAULT(NONE) SHARED(V_mats,norm1,mos_old,ncol_local)
    DO icol_local=1,ncol_local
       V_mats(1)%matrix%local_data(:,icol_local)=mos_old(1)%matrix%local_data(:,icol_local)
       V_mats(1)%matrix%local_data(:,icol_local+ncol_local)=mos_old(2)%matrix%local_data(:,icol_local)
       norm1(icol_local)=SUM(V_mats(1)%matrix%local_data(:,icol_local)**2)+SUM(V_mats(1)%matrix%local_data(:,icol_local)**2)
    END DO

    CALL mp_sum(norm1,col_group)
    !!! normalize the mo vectors
    norm1(:)=SQRT(norm1(:))

!$OMP PARALLEL DO PRIVATE(icol_local) DEFAULT(NONE) SHARED(V_mats,norm1,ncol_local)
    DO icol_local=1,ncol_local
       V_mats(1)%matrix%local_data(:,icol_local)=V_mats(1)%matrix%local_data(:,icol_local)/norm1(icol_local)
       V_mats(1)%matrix%local_data(:,icol_local+ncol_local)=V_mats(1)%matrix%local_data(:,icol_local+ncol_local)/norm1(icol_local)
    END DO

    ! arnoldi subspace procedure to get H_approx
    DO i=2,narnoldi+1
       !Be careful, imaginary matrix multiplied with complex. Unfortunately requires a swap of arrays afterwards
       CALL cp_gemm("N","N",nao,newdim,nao,1.0_dp,Him,V_mats(i-1)%matrix,0.0_dp,V_mats(i)%matrix,error)

!$OMP PARALLEL DO PRIVATE(icol_local) DEFAULT(NONE) SHARED(mos_new,V_mats,ncol_local,i)
       DO icol_local=1,ncol_local
          mos_new(1)%matrix%local_data(:,icol_local)=-V_mats(i)%matrix%local_data(:,icol_local)
          V_mats(i)%matrix%local_data(:,icol_local)=V_mats(i)%matrix%local_data(:,icol_local+ncol_local)
          V_mats(i)%matrix%local_data(:,icol_local+ncol_local)=mos_new(1)%matrix%local_data(:,icol_local)
       END DO

       IF(PRESENT(Hre))THEN
           CALL cp_gemm("N","N",nao,newdim,nao,1.0_dp,Hre,V_mats(i-1)%matrix,1.0_dp,V_mats(i)%matrix,error)
       END IF

       DO l=1,i-1
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(results,V_mats,ncol_local,l,i)
          DO icol_local=1,ncol_local
             results(icol_local)=SUM(V_mats(l)%matrix%local_data(:,icol_local)*V_mats(i)%matrix%local_data(:,icol_local))+&
                SUM(V_mats(l)%matrix%local_data(:,icol_local+ncol_local)*V_mats(i)%matrix%local_data(:,icol_local+ncol_local))
          END DO

          CALL mp_sum(results,col_group)

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(H_approx_save,V_mats,ncol_local,l,i,results)
          DO icol_local=1,ncol_local
             H_approx_save(l,i-1,icol_local)=results(icol_local)
             V_mats(i)%matrix%local_data(:,icol_local)=V_mats(i)%matrix%local_data(:,icol_local)-&
                  results(icol_local) * V_mats(l)%matrix%local_data(:,icol_local)
             V_mats(i)%matrix%local_data(:,icol_local+ncol_local)=V_mats(i)%matrix%local_data(:,icol_local+ncol_local)-&
                  results(icol_local) * V_mats(l)%matrix%local_data(:,icol_local+ncol_local)
          END DO
       END DO

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(ncol_local,V_mats,results,i)
       DO icol_local=1,ncol_local
          results(icol_local)=SUM(V_mats(i)%matrix%local_data(:,icol_local)**2)+&
                               SUM(V_mats(i)%matrix%local_data(:,icol_local+ncol_local)**2)
       END DO

       CALL mp_sum(results,col_group)


       IF(i.LE.narnoldi)THEN

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(H_approx_save,last_norm,V_mats,ncol_local,i,results)
          DO icol_local=1,ncol_local
             H_approx_save(i,i-1,icol_local)=SQRT(results(icol_local))
             last_norm(icol_local)=SQRT(results(icol_local))
             V_mats(i)%matrix%local_data(:,icol_local)=V_mats(i)%matrix%local_data(:,icol_local)/SQRT(results(icol_local))
             V_mats(i)%matrix%local_data(:,icol_local+ncol_local)=&
                  V_mats(i)%matrix%local_data(:,icol_local+ncol_local)/SQRT(results(icol_local))
          END DO
       ELSE
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(ncol_local,last_norm,results)
          DO icol_local=1,ncol_local
             last_norm(icol_local)=SQRT(results(icol_local))
          END DO
       END IF

       H_approx(:,:,:)=H_approx_save

       ! PADE approximation for exp(H_approx), everything is done locally

       convergence=.FALSE.
       IF(i.ge.narn_old)THEN
          npade=9
          mydim=MIN(i,narnoldi)
          ALLOCATE(ipivot(mydim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(mat1(mydim,mydim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(mat2(mydim,mydim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(mat3(mydim,mydim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(N(mydim,mydim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(D(mydim,mydim),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO icol_local=1,ncol_local
             DO idim=1,mydim
                DO j=1,mydim
                   mat1(idim,j)=H_approx(idim,j,icol_local)/16.0_dp
                   mat3(idim,j)=mat1(idim,j)
                END DO
             END DO
             N=0.0_dp
             D=0.0_dp
             DO idim=1,mydim
                N(idim,idim)=1.0_dp
                D(idim,idim)=one
             END DO
             N(:,:)=N+0.5_dp*mat1
             D(:,:)=D-0.5_dp*mat1
             pade_step=1
             DO idim=1,4
                pade_step=pade_step+1
                CALL dgemm("N",'N',mydim,mydim,mydim,one,mat1(1,1),&
                     mydim,mat3(1,1),mydim,zero,mat2(1,1),mydim)
                scaleN=REAL(fac(2*npade -pade_step )*fac(npade )/&
                     (fac(2*npade )*fac(pade_step )*fac(npade -pade_step )),dp)
                scaleD=REAL((-1.0_dp)**pade_step  * fac(2*npade -pade_step )*fac(npade)/&
                     (fac(2*npade )*fac(pade_step )*fac(npade -pade_step)),dp)
                N(:,:)=N+scaleN*mat2
                D(:,:)=D+scaleD*mat2
                pade_step=pade_step+1
                CALL dgemm("N",'N',mydim,mydim,mydim,one,mat2(1,1),&
                     mydim,mat1(1,1),mydim,zero,mat3(1,1),mydim)
                scaleN=REAL(fac(2*npade -pade_step )*fac(npade )/&
                     (fac(2*npade )*fac(pade_step )*fac(npade -pade_step )),dp)
                scaleD=REAL((-1.0_dp)**pade_step  * fac(2*npade -pade_step )*fac(npade)/&
                     (fac(2*npade )*fac(pade_step )*fac(npade -pade_step)),dp)
                N(:,:)=N+scaleN*mat3
                D(:,:)=D+scaleD*mat3
             END DO
   
             CALL dgetrf(mydim,mydim,D(1,1),mydim,ipivot,info)
             CALL dgetrs("N",mydim,mydim,D(1,1),mydim,ipivot,N,mydim,info)
             CALL dgemm("N",'N',mydim,mydim,mydim,one,N(1,1),mydim,N(1,1),mydim,zero,mat1(1,1),mydim)
             CALL dgemm("N",'N',mydim,mydim,mydim,one,mat1(1,1),mydim,mat1(1,1),mydim,zero,N(1,1),mydim)
             CALL dgemm("N",'N',mydim,mydim,mydim,one,N(1,1),mydim,N(1,1),mydim,zero,mat1(1,1),mydim)
             CALL dgemm("N",'N',mydim,mydim,mydim,one,mat1(1,1),mydim,mat1(1,1),mydim,zero,N(1,1),mydim)
             DO idim=1,mydim
                DO j=1,mydim
                   H_approx(idim,j,icol_local)=N(idim,j)
                END DO
             END DO
          END DO
          ! H_approx is exp(H_approx) right now, calculate new MOs and check for convergence
          conv_norm=0.0_dp
          results=0.0_dp
          DO icol_local=1,ncol_local
             results(icol_local)= last_norm(icol_local)* H_approx(i-1,1,icol_local)
             conv_norm=MAX(conv_norm,ABS(results(icol_local)))
          END DO
   
          CALL mp_max(conv_norm,para_env%group)
   
          IF(conv_norm.LT.eps_exp.OR.i.GT.narnoldi)THEN
   
             mos_new(1)%matrix%local_data=zero
             mos_new(2)%matrix%local_data=zero
             DO icol_local=1,ncol_local
                DO idim=1,mydim
                   prefac=H_approx(idim,1,icol_local)*norm1(icol_local)
                   mos_new(1)%matrix%local_data(:,icol_local)=mos_new(1)%matrix%local_data(:,icol_local)+&
                        V_mats(idim)%matrix%local_data(:,icol_local)*prefac
                   mos_new(2)%matrix%local_data(:,icol_local)=mos_new(2)%matrix%local_data(:,icol_local)+&
                        V_mats(idim)%matrix%local_data(:,icol_local+ncol_local)*prefac
                END DO
             END DO

             IF(PRESENT(mos_next))THEN
               DO icol_local=1,ncol_local
                 DO idim=1,mydim
                    DO j=1,mydim
                       N(idim,j)=H_approx(idim,j,icol_local)
                    END DO
                 END DO
                 CALL dgemm("N",'N',mydim,mydim,mydim,one,N(1,1),mydim,N(1,1),mydim,zero,mat1(1,1),mydim)
                 DO idim=1,mydim
                   DO j=1,mydim
                     H_approx(idim,j,icol_local)=mat1(idim,j)
                   END DO
                 END DO
               END DO
               mos_next(1)%matrix%local_data=zero
               mos_next(2)%matrix%local_data=zero
               DO icol_local=1,ncol_local
                  DO idim=1,mydim
                     prefac=H_approx(idim,1,icol_local)*norm1(icol_local)
                     mos_next(1)%matrix%local_data(:,icol_local)=mos_next(1)%matrix%local_data(:,icol_local)+&
                          V_mats(idim)%matrix%local_data(:,icol_local)*prefac
                     mos_next(2)%matrix%local_data(:,icol_local)=mos_next(2)%matrix%local_data(:,icol_local)+&
                          V_mats(idim)%matrix%local_data(:,icol_local+ncol_local)*prefac
                  END DO
               END DO
             END IF
             IF(conv_norm.LT.eps_exp)THEN
                convergence=.TRUE.
                narn_old=i-1
             END IF
          END IF
   
          DEALLOCATE(ipivot,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(mat1,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(mat2,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(mat3,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(N,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(D,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(convergence)EXIT

    END DO
    CALL cp_assert(convergence,cp_warning_level,cp_assertion_failed,routineP,&
         "ARNOLDI method did not converge"//&
         CPSourceFileRef,&
         only_ionode=.TRUE.)
    !deallocate all work matrices

    DO i=1,SIZE(V_mats)
       CALL cp_fm_release(V_mats(i)%matrix,error)
    END DO
    CALL cp_fm_struct_release(newstruct,error)
    CALL mp_comm_free(col_group)

    DEALLOCATE(V_mats,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(H_approx,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(H_approx_save,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(results,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(norm1,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(last_norm,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
  END SUBROUTINE arnoldi


END MODULE rt_matrix_exp
