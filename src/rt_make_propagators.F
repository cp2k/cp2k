
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for calculating a complex matrix exponential.
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE rt_make_propagators

  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_add,&
                                             cp_cfm_solve
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_set_all,&
                                             cp_cfm_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale,&
                                             cp_fm_upper_to_full
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE input_constants,                 ONLY: do_diag,&
                                             do_em,&
                                             do_etrs,&
                                             do_pade,&
                                             do_taylor
  USE kinds,                           ONLY: dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE rt_matrix_exp,                   ONLY: arnoldi,&
                                             exp_diag,&
                                             exp_pade_full_complex,&
                                             exp_pade_only_imaginary,&
                                             taylor_full_complex,&
                                             taylor_only_imaginary
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_make_propagators'

  PUBLIC :: make_etrs_exp,&
            make_cn_propagator,&
            make_em_propagator,&
            propagate_arnoldi


CONTAINS
! *****************************************************************************
!> \brief calculates the exponentials for the ETRS propagator
!>        depending on input exp(i*H(t)*dt/2) or   exp(i*H(t+dt)*dt/2)
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE make_etrs_exp(exp_H,qs_env,error)


    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_etrs_exp', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, ispin, lb, ndim, ub
    LOGICAL                                  :: failure
    REAL(dp)                                 :: t
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: SinvH
    TYPE(cp_fm_type), POINTER                :: B_mat, H_fm, S_inv, SinvB, &
                                                tmp_mat_B, tmp_mat_H
    TYPE(rt_prop_type), POINTER              :: rtp

    failure=.FALSE.
    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,matrix_ks=matrix_ks,rtp=rtp,error=error)
    CALL get_rtp(rtp=rtp,S_inv=S_inv,dt=t,error=error)
    IF(.NOT.qs_env%dft_control%rtp_control%fixed_ions)CALL get_rtp(rtp=rtp,SinvH=SinvH,error=error)
    CALL cp_fm_create(H_fm,&
         matrix_struct=exp_H(1)%matrix%matrix_struct,&
         name="RTP_H_FM",&
         error=error)
    CALL cp_fm_create(tmp_mat_H,&
         matrix_struct=exp_H(1)%matrix%matrix_struct,&
         name="TMP_H",&
         error=error)

    ndim=H_fm%matrix_struct%nrow_global

    t=-t/2.0_dp
    DO ispin=1,SIZE(matrix_ks)

       lb=(ispin-1)*2+1
       ub=(ispin-1)*2+2

       CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,H_fm,error=error)

       CALL cp_fm_upper_to_full(H_fm,tmp_mat_H,error)


       SELECT CASE(qs_env%dft_control%rtp_control%mat_exp)

       CASE(do_taylor,do_pade)
          CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,rtp%S_inv,H_fm,zero,tmp_mat_H,error)
          IF(qs_env%dft_control%rtp_control%fixed_ions)THEN
             CALL cp_fm_scale(t,tmp_mat_H,error)
             IF(qs_env%dft_control%rtp_control%mat_exp==do_taylor)&
                  CALL taylor_only_imaginary(exp_H(lb:ub),tmp_mat_H,rtp%orders(1,ispin),rtp%orders(2,ispin),error)

             IF(qs_env%dft_control%rtp_control%mat_exp==do_pade)&
                  CALL exp_pade_only_imaginary(exp_H(lb:ub),tmp_mat_H,rtp%orders(1,ispin),rtp%orders(2,ispin),error)
          ELSE
             CALL cp_fm_to_fm(tmp_mat_H,SinvH(ispin)%matrix,error)
             CALL cp_fm_scale(t,tmp_mat_H,error)
             IF(ispin==1)THEN
                CALL get_rtp(rtp=rtp,B_mat=B_mat,SinvB=SinvB,error=error)

                CALL cp_fm_create(tmp_mat_B,&
                     matrix_struct=exp_H(1)%matrix%matrix_struct,&
                     name="TMP_B",&
                     error=error)

                CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,rtp%S_inv,B_mat,zero,tmp_mat_B,error)
                CALL cp_fm_to_fm(tmp_mat_B,SinvB,error)

                CALL cp_fm_scale(t,tmp_mat_B,error)
             END IF
             IF(qs_env%dft_control%rtp_control%mat_exp==do_taylor)&
                  CALL taylor_full_complex(exp_H(lb:ub),tmp_mat_B,tmp_mat_H,rtp%orders(1,ispin),rtp%orders(2,ispin),error)
             IF(qs_env%dft_control%rtp_control%mat_exp==do_pade)&
                  CALL exp_pade_full_complex(exp_H(lb:ub),tmp_mat_B,tmp_mat_H,rtp%orders(1,ispin),rtp%orders(2,ispin),&
                  qs_env%sim_step,error)
          END IF

       CASE(do_diag)
          CALL cp_assert(qs_env%dft_control%rtp_control%fixed_ions,&
               cp_failure_level,cp_assertion_failed,&
               routineP,"diagonalization does not yet work for the Ehrenfest scheme",&
               error,failure)
          CALL exp_diag(rtp,exp_H(lb:ub),H_fm,error)
       END SELECT

    END DO

    IF(.NOT.qs_env%dft_control%rtp_control%fixed_ions)CALL cp_fm_release(tmp_mat_B,error)
    CALL cp_fm_release(H_fm,error)
    CALL cp_fm_release(tmp_mat_H,error)

    CALL timestop(handle)

  END SUBROUTINE make_etrs_exp

! *****************************************************************************
!> \brief computes the Crank-Nicolson propagator
!>        [1-idt/2H(t+dt/2)-dt/2B(t+dt/2)]/ (1+idt/2H(t+dt/2)+dt/2B(t+dt/2))
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE make_cn_propagator(propagator,rtp,error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: propagator
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_cn_propagator', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: a
    INTEGER                                  :: i, im, ispin, ldim, ndim, re
    REAL(dp)                                 :: dt
    TYPE(cp_cfm_type), POINTER               :: M1, M2, M_inverse
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: new_Ham, old_Ham
    TYPE(cp_fm_type), POINTER                :: U_im, U_re

    CALL cp_cfm_create(M1,&
                      matrix_struct=propagator(1)%matrix%matrix_struct,&
                      name="M1",&
                      error=error)
    CALL cp_cfm_create(M2,&
                      matrix_struct=propagator(1)%matrix%matrix_struct,&
                      name="M2",&
                      error=error)
    CALL cp_cfm_create(M_inverse,&
                      matrix_struct=propagator(1)%matrix%matrix_struct,&
                      name="M_inverse",&
                      error=error)


    a=> M_inverse%local_data
    CALL get_rtp(rtp=rtp,exp_H_old=old_Ham,exp_H_new=new_Ham,dt=dt,error=error)

    CALL cp_fm_get_info(propagator(1)%matrix,ncol_local=ldim,nrow_global=ndim,error=error)
    DO ispin=1,SIZE(propagator)/2
       re=2*ispin-1
       im=2*ispin

       U_re=>propagator(re)%matrix
       U_im=>propagator(im)%matrix

       CALL cp_cfm_set_all(M1,zero,one,error)
       CALL cp_cfm_set_all(M2,zero,one,error)

       dt=dt*0.5_dp

       !M_inverse is used right here as temporary matrix
       DO i=1,ldim
          a(:,i)=&
               CMPLX(0.5_dp*dt*old_ham(re)%matrix%local_data(:,i)+0.5_dp*dt*new_ham(re)%matrix%local_data(:,i),&
               0.5_dp*dt*old_ham(im)%matrix%local_data(:,i)+0.5_dp*dt*new_ham(im)%matrix%local_data(:,i),KIND=dp)
       END DO

       CALL cp_cfm_add(one,M1,-one,M_inverse,error)
       CALL cp_cfm_add(one,M2,one,M_inverse,error)

       CALL cp_cfm_solve(M2,M1,error=error)

       !now M2 is temporary

       U_re%local_data=REAL(M1%local_data,KIND=dp)
       U_im%local_data=AIMAG(M1%local_data)


    END DO

    CALL cp_cfm_release(M1,error)
    CALL cp_cfm_release(M2,error)
    CALL cp_cfm_release(M_inverse,error)

  END SUBROUTINE make_cn_propagator

! *****************************************************************************
!> \brief computes exponential midpoint propagator
!>        exp[-idt/2H(t+dt/2)-tB(t+dt/2)]
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE make_em_propagator(propagator,rtp,qs_env,error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: propagator
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_em_propagator', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)

    INTEGER                                  :: handle, i, im, ispin, ldim, &
                                                ndim, re
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dt, eps_exp
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: new_Ham, old_Ham
    TYPE(cp_fm_type), POINTER                :: M1, M2

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CALL cp_fm_create(M1,&
                      matrix_struct=propagator(1)%matrix%matrix_struct,&
                      name="M1",&
                      error=error)
    CALL cp_fm_create(M2,&
                      matrix_struct=propagator(1)%matrix%matrix_struct,&
                      name="M2",&
                      error=error)


    CALL cp_fm_get_info(M1,ncol_local=ldim,nrow_global=ndim,error=error)
    CALL get_rtp(rtp=rtp,exp_H_old=old_Ham,exp_H_new=new_Ham,dt=dt,error=error)


    DO ispin=1,SIZE(propagator)/2
       re=2*ispin-1
       im=2*ispin

       dt=-dt

       DO i=1,ldim
          M2%local_data(:,i)=0.5_dp*dt*old_ham(im)%matrix%local_data(:,i)+0.5_dp*dt*new_ham(im)%matrix%local_data(:,i)
       END DO
       IF(.NOT.qs_env%dft_control%rtp_control%fixed_ions)THEN
          DO i=1,ldim
             M1%local_data(:,i)=0.5_dp*dt*old_ham(re)%matrix%local_data(:,i)+0.5_dp*dt*new_ham(re)%matrix%local_data(:,i)
          END DO
       END IF
       SELECT CASE(qs_env%dft_control%rtp_control%mat_exp)

       CASE(do_taylor)
          eps_exp=qs_env%dft_control%rtp_control%eps_exp
          IF(qs_env%dft_control%rtp_control%fixed_ions)THEN
             CALL taylor_only_imaginary(propagator(re:im),M2,rtp%orders(1,ispin),rtp%orders(2,ispin),error)
          ELSE
             CALL taylor_full_complex(propagator(re:im),M1,M2,rtp%orders(1,ispin),rtp%orders(2,ispin),error)
          END IF
       CASE(do_pade)
          eps_exp=qs_env%dft_control%rtp_control%eps_exp
          IF(qs_env%dft_control%rtp_control%fixed_ions)THEN
             CALL exp_pade_only_imaginary(propagator(re:im),M2,rtp%orders(1,ispin),rtp%orders(2,ispin),error)
          ELSE
             CALL exp_pade_full_complex(propagator(re:im),M1,M2,rtp%orders(1,ispin),rtp%orders(2,ispin),qs_env%sim_step,error)
          END IF

       CASE(do_diag)
          CALL cp_assert(.FALSE.,&
               cp_failure_level,cp_assertion_failed,&
               routineP,"diagonalization does not work for EM propagator scheme",&
               error,failure)
       END SELECT

    END DO
    CALL cp_fm_release(M1,error)
    CALL cp_fm_release(M2,error)
    CALL timestop(handle)

  END SUBROUTINE make_em_propagator

! *****************************************************************************
!> \brief computes U_prop*MOs using arnoldi subspace algorithm
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE propagate_arnoldi(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'propagate_arnoldi', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, ispin, lb, ldim, &
                                                ndim, nspin, ub
    REAL(dp)                                 :: dt, eps_arnoldi, t
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, mos_next, mos_old, &
                                                new_Ham, old_Ham
    TYPE(cp_fm_type), POINTER                :: tmp_mat_B, tmp_mat_H
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,rtp=rtp,error=error)
    CALL get_rtp(rtp=rtp,dt=t,mos_new=mos_new,mos_old=mos_old,&
         mos_next=mos_next,exp_H_old=old_Ham,exp_H_new=new_Ham,error=error)

    CALL cp_fm_create(tmp_mat_H,&
         matrix_struct=old_Ham(1)%matrix%matrix_struct,&
         name="TMP_H",&
         error=error)
    CALL cp_fm_create(tmp_mat_B,&
         matrix_struct=old_Ham(1)%matrix%matrix_struct,&
         name="TMP_B",&
         error=error)
    nspin=SIZE(mos_new)/2
    eps_arnoldi=qs_env%dft_control%rtp_control%eps_exp
    SELECT CASE(qs_env%dft_control%rtp_control%propagator)
    CASE(do_etrs)
       IF(rtp%iter==1.AND.rtp%istep==1)THEN
          dt=-t/2.0_dp
          DO ispin=1,nspin

             lb=ispin*2-1
             ub=ispin*2
             CALL cp_fm_to_fm(old_ham(ub)%matrix,tmp_mat_H,error)
             CALL cp_fm_scale(dt,tmp_mat_H,error)

             IF(qs_env%dft_control%rtp_control%fixed_ions)THEN
                CALL arnoldi(mos_old(lb:ub),mos_new(lb:ub),&
                     eps_arnoldi,Him=tmp_mat_H,mos_next=mos_next(lb:ub),error=error)
             ELSE
                CALL cp_fm_to_fm(old_ham(lb)%matrix,tmp_mat_B,error)
                CALL cp_fm_scale(dt,tmp_mat_B,error)

                CALL arnoldi(mos_old(lb:ub),mos_new(lb:ub),&
                     eps_arnoldi,Hre=tmp_mat_B,Him=tmp_mat_H,mos_next=mos_next(lb:ub),error=error)
             END IF
          END DO
          DO i=1,SIZE(mos_new)
             CALL cp_fm_to_fm(mos_new(i)%matrix,mos_old(i)%matrix,error)
             CALL cp_fm_to_fm(mos_next(i)%matrix,mos_new(i)%matrix,error)

          END DO
       ELSE
          IF(rtp%iter==1)THEN
             DO i=1,SIZE(mos_new)
                CALL cp_fm_to_fm(mos_next(i)%matrix,mos_old(i)%matrix,error)
             END DO
          END IF
       dt=-t/2.0_dp

       DO ispin=1,nspin

          lb=ispin*2-1
          ub=ispin*2
          CALL cp_fm_to_fm(new_ham(ub)%matrix,tmp_mat_H,error)
          CALL cp_fm_scale(dt,tmp_mat_H,error)
          IF(qs_env%dft_control%rtp_control%fixed_ions)THEN
             CALL arnoldi(mos_old(lb:ub),mos_new(lb:ub),&
                  eps_arnoldi,Him=tmp_mat_H,mos_next=mos_next(lb:ub),error=error)
          ELSE
             CALL cp_fm_to_fm(new_ham(lb)%matrix,tmp_mat_B,error)
             CALL cp_fm_scale(dt,tmp_mat_B,error)
             CALL arnoldi(mos_old(lb:ub),mos_new(lb:ub),&
                  eps_arnoldi,Hre=tmp_mat_B,Him=tmp_mat_H,&
                  mos_next=mos_next(lb:ub),error=error)
          END IF

       END DO
       END IF
       CASE(do_em)
          CALL cp_fm_get_info(tmp_mat_H,ncol_local=ldim,nrow_global=ndim,error=error)
          DO ispin=1,nspin
             lb=2*ispin-1
             ub=2*ispin
             dt=-t
             DO i=1,ldim
                tmp_mat_H%local_data(:,i)=0.5_dp*dt*old_ham(ub)%matrix%local_data(:,i)+&
                     0.5_dp*dt*new_ham(ub)%matrix%local_data(:,i)
             END DO
             IF(.NOT.qs_env%dft_control%rtp_control%fixed_ions)THEN
                DO i=1,ldim
                   tmp_mat_B%local_data(:,i)=0.5_dp*dt*old_ham(lb)%matrix%local_data(:,i)+&
                        0.5_dp*dt*new_ham(lb)%matrix%local_data(:,i)
                END DO
                CALL arnoldi(mos_old(lb:ub),mos_new(lb:ub),&
                     eps_arnoldi,Hre=tmp_mat_B,Him=tmp_mat_H,error=error)
             ELSE
                CALL arnoldi(mos_old(lb:ub),mos_new(lb:ub),&
                     eps_arnoldi,Him=tmp_mat_H,error=error)
             END IF


          END DO
       END SELECT
       CALL cp_fm_release(tmp_mat_H,error)
       CALL cp_fm_release(tmp_mat_B,error)
       CALL timestop(handle)

     END SUBROUTINE propagate_arnoldi

END MODULE rt_make_propagators
