!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for the prepares teh matrices for the propagators
!> \author Florian Schiffmann (02.09)
! *****************************************************************************
MODULE rt_propagator_init

  USE cp_control_types,                ONLY: dft_control_type,&
                                             rtp_control_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_upper_to_full
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE input_constants,                 ONLY: do_arnoldi,&
                                             do_cn,&
                                             do_em,&
                                             do_etrs,&
                                             do_pade,&
                                             do_taylor
  USE kinds,                           ONLY: dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE rt_make_propagators,             ONLY: compute_exponential,&
                                             propagate_arnoldi
  USE rt_matrix_exp,                   ONLY: backtransform_matrix,&
                                             get_nsquare_norder
  USE rt_propagation_methods,          ONLY: calc_SinvH,&
                                             calc_update_rho,&
                                             put_data_to_history,&
                                             s_matrices_create
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE rt_propagation_utils,            ONLY: calc_S_derivs
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagator_init'

  PUBLIC :: init_propagators,&
            init_emd_propagators


CONTAINS


! *****************************************************************************
!> \brief prepares the initial matrices for the propagators
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE init_propagators(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_propagators', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, imat
    REAL(KIND=dp)                            :: dt, prefac
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: s_mat
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, &
                                                mos_new, mos_next, mos_old, &
                                                propagator_matrix
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control

    CALL get_qs_env(qs_env,&
         rtp=rtp,&
         dft_control=dft_control,&
         matrix_s=s_mat,&
         error=error)

    rtp_control=>dft_control%rtp_control
    CALL get_rtp (rtp=rtp,exp_H_old=exp_H_old,exp_H_new=exp_H_new,mos_old=mos_old,&
                  propagator_matrix=propagator_matrix,dt=dt,error=error)


    CALL s_matrices_create (s_mat,rtp,error)
    CALL calc_SinvH(exp_H_new,qs_env,error)
    DO i=1,SIZE(exp_H_old)
       CALL cp_fm_to_fm(exp_H_new(i)%matrix,exp_H_old(i)%matrix,error)
    END DO

    ! use the fact that CN propagator is a first order pade approximation on the EM propagator
    IF(rtp_control%propagator==do_cn)THEN
       rtp%orders(1,:)=0;rtp%orders(2,:)=1; rtp_control%mat_exp=do_pade; rtp_control%propagator=do_em
    ELSE
       IF(rtp_control%mat_exp==do_pade.OR.rtp_control%mat_exp==do_taylor)&
          CALL get_maxabs_eigval(rtp,qs_env,error)
    END IF

    ! We have no clue yet about next H so we use initial H for t and t+dt
    ! Due to different nature of the propagator the prefactor has to be adopted

    SELECT CASE(rtp_control%propagator)
    CASE(do_etrs)
      prefac=-0.5_dp*dt
    CASE(do_em)
      prefac=-1.0_dp*dt
    END SELECT 
    DO imat=1,SIZE(exp_H_new)
        CALL cp_fm_scale_and_add(0.0_dp,propagator_matrix(imat)%matrix,prefac,exp_H_new(imat)%matrix,error)
    END DO

    ! For ETRS this bit could be avoided but it drastically simplifies the workflow afterwards.
    ! If we compute the half propagated mos/exponential already here, we ensure everything is done
    ! with the correct S matrix and all information as during RTP/EMD are computed.
    ! Therefore we might accept to compute an unnesscesary expoenential but understand the code afterwards
    IF(rtp_control%propagator==do_etrs)THEN
       IF(rtp_control%mat_exp==do_arnoldi)THEN
          rtp%iter=0
          CALL propagate_arnoldi(rtp,rtp_control,error)
          CALL get_rtp(rtp=rtp,mos_new=mos_new,mos_next=mos_next,error=error)
          DO imat=1,SIZE(mos_new)
             CALL cp_fm_to_fm(mos_new(imat)%matrix,mos_next(imat)%matrix,error)
          END DO
       ELSE
          CALL compute_exponential(exp_H_new,propagator_matrix,rtp_control,rtp,error)
          DO imat=1,SIZE(exp_H_new)
             CALL cp_fm_to_fm(exp_H_new(imat)%matrix,exp_H_old(imat)%matrix,error)
          END DO
       END IF
    END IF 

    CALL put_data_to_history(rtp,mos=mos_old,s_mat=s_mat,ihist=1,error=error)

  END SUBROUTINE init_propagators

! *****************************************************************************
!> \brief gets an estimate for the 2-norm of KS (diagnaliztion of KS) and
!>        calculates the order and number of squaring steps for Taylor or
!>        Pade matrix exponential
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE get_maxabs_eigval(rtp,qs_env,error)
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_maxabs_eigval', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: ispin, method, ndim, stat
    LOGICAL                                  :: emd, failure
    REAL(dp)                                 :: max_eval, min_eval, norm2, &
                                                scale, t
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eigval_H
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, s_mat
    TYPE(cp_fm_type), POINTER                :: eigvec_H, H_fm, S_half, &
                                                S_inv, S_minus_half, tmp, &
                                                tmp_mat_H

    failure=.FALSE.

    CALL get_qs_env(qs_env=qs_env,matrix_ks=matrix_ks,matrix_s=s_mat,error=error)
    CALL get_rtp(rtp=rtp,S_inv=S_inv,dt=t,error=error)

    CALL cp_fm_create(S_half,&
          matrix_struct=S_inv%matrix_struct,&
          name="S_half",&
          error=error)

    CALL cp_fm_create(S_minus_half,&
          matrix_struct=S_inv%matrix_struct,&
          name="S_minus_half",&
          error=error)

    CALL cp_fm_create(H_fm,&
         matrix_struct=S_inv%matrix_struct,&
         name="RTP_H_FM",&
         error=error)

    CALL cp_fm_create(tmp_mat_H,&
         matrix_struct=S_inv%matrix_struct,&
         name="TMP_H",&
         error=error)

    ndim=S_inv%matrix_struct%nrow_global
    scale=1.0_dp
    IF(qs_env%dft_control%rtp_control%propagator==do_etrs)scale=2.0_dp
    t=-t/scale

    ! Create the overlap matrices

    CALL cp_fm_create(tmp,&
         matrix_struct=S_inv%matrix_struct,&
         name="tmp_mat",&
         error=error)

    CALL cp_fm_create(eigvec_H,&
         matrix_struct=S_inv%matrix_struct,&
         name="tmp_EVEC",&
         error=error)

    ALLOCATE(eigval_H(ndim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL copy_dbcsr_to_fm(s_mat(1)%matrix,tmp,error=error)
    CALL cp_fm_upper_to_full(tmp,eigvec_H,error)

    CALL cp_fm_syevd(tmp,eigvec_H,eigval_H,error=error)

    eigval_H(:)=one/eigval_H(:)
    CALL backtransform_matrix(eigval_H,eigvec_H,S_inv,error)
    eigval_H(:)=SQRT(eigval_H(:))
    CALL backtransform_matrix(eigval_H,eigvec_H,S_minus_half,error)
    eigval_H(:)=one/eigval_H(:)
    CALL backtransform_matrix(eigval_H,eigvec_H,S_half,error)
    CALL cp_fm_release(eigvec_H,error)
    CALL cp_fm_release(tmp,error)

    IF(qs_env%dft_control%rtp_control%mat_exp==do_taylor)method=1
    IF(qs_env%dft_control%rtp_control%mat_exp==do_pade)method=2
    emd=(.NOT.qs_env%dft_control%rtp_control%fixed_ions)

    DO ispin=1,SIZE(matrix_ks)

       CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,H_fm,error=error)
       CALL cp_fm_upper_to_full(H_fm,tmp_mat_H,error)
       CALL cp_fm_scale(t,H_fm,error)

       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,H_fm,S_minus_half,zero,&
            tmp_mat_H,error)
       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,S_minus_half,tmp_mat_H,zero,&
            H_fm,error)

       CALL cp_fm_syevd(H_fm,tmp_mat_H,eigval_H,error=error)
       min_eval=MINVAL(eigval_H)
       max_eval=MAXVAL(eigval_H)
       norm2=2.0_dp*MAX(ABS(min_eval),ABS(max_eval))
       CALL get_nsquare_norder(norm2,rtp%orders(1,ispin),rtp%orders(2,ispin),&
            qs_env%dft_control%rtp_control%eps_exp,method ,emd,error)
    END DO

    DEALLOCATE(eigval_H,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL cp_fm_release(S_half,error)
    CALL cp_fm_release(S_minus_half,error)
    CALL cp_fm_release(H_fm,error)
    CALL cp_fm_release(tmp_mat_H,error)

  END SUBROUTINE get_maxabs_eigval


! *****************************************************************************
!> \brief computs all matrices needed to initialize emd (S_derivs, SinvH,
!>        SinvB and forces)
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE init_emd_propagators(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_emd_propagators', &
      routineP = moduleN//':'//routineN

    TYPE(qs_ks_env_type), POINTER            :: ks_env

!density matrix

       CALL calc_update_rho(qs_env,error)
       CALL get_qs_env(qs_env=qs_env,ks_env=ks_env,error=error)
       !derivatives of the potential
       CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,calculate_forces=.FALSE.,&
            error=error)
       !derivativs of the overlap needed for EMD
       CALL calc_S_derivs(qs_env,error)
       ! a bit hidden, but computs SinvH and SinvB (calc_SinvH for CN,EM and ARNOLDI)
       ! make_etrs_exp in case of ETRS in combination with TAYLOR and PADE
       CALL init_propagators(qs_env,error)

   END SUBROUTINE init_emd_propagators

END MODULE rt_propagator_init
