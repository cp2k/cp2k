!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for propagating the orbitals
!> \author Florian Schiffmann (02.09)
! *****************************************************************************
MODULE rt_propagation_methods
  USE bibliography,                    ONLY: Kolafa2004,&
                                             cite_reference
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_decompose,&
                                             cp_cfm_gemm,&
                                             cp_cfm_triangular_multiply
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_type
  USE cp_control_types,                ONLY: rtp_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_copy,&
                                             cp_dbcsr_init
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_double,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE input_constants,                 ONLY: do_arnoldi,&
                                             do_em,&
                                             do_etrs,&
                                             do_pade,&
                                             do_taylor
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: binomial
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: mo_set_p_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE rt_make_propagators,             ONLY: propagate_arnoldi,&
                                             propagate_exp
  USE rt_propagation_output,           ONLY: rt_convergence
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_methods'

  PUBLIC :: propagation_step,&
       s_matrices_create,&
       calc_update_rho,&
       aspc_extrapolate,&
       calc_sinvH,&
       put_data_to_history


CONTAINS

! *****************************************************************************
!> \brief performes a single propagation step a(t+Dt)=U(t+Dt,t)*a(0)
!>        and calculates the new exponential (at the moment for ETRS)
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE propagation_step(qs_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'propagation_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: aspc_order, handle, i
    LOGICAL                                  :: my_force
    REAL(dp)                                 :: delta_iter
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: delta_mos, exp_H_new, mos_new
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control

    CALL timeset(routineN,handle)
    NULLIFY(rho,ks_env,rtp)
    ! get everything needed and set some values
    CALL get_qs_env(qs_env,&
         rho=rho,&
         ks_env=ks_env,&
         rtp=rtp,&
         error=error)
    CALL get_rtp(rtp=rtp,exp_H_new=exp_H_new, mos_new=mos_new, error=error)
    rtp_control=>qs_env%dft_control%rtp_control

    my_force=.NOT.rtp_control%fixed_ions
    delta_iter=1.0_dp

    ALLOCATE(delta_mos(SIZE(mos_new)))
    DO i=1,SIZE(mos_new)
       CALL cp_fm_create(delta_mos(i)%matrix,&
            matrix_struct=mos_new(i)%matrix%matrix_struct,&
            name="delta_mos"//TRIM(ADJUSTL(cp_to_string(i))),&
            error=error)
       CALL cp_fm_to_fm(mos_new(i)%matrix,delta_mos(i)%matrix,error)
    END DO

    IF(rtp%iter==1)THEN
       aspc_order=rtp_control%aspc_order
       CALL aspc_extrapolate(qs_env,rtp,aspc_order,error)
       IF(rtp%history%mos_or_H==1)THEN
          CALL calc_update_rho(qs_env,error)
          CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,calculate_forces=.FALSE.,&
               error=error)
       END IF
    END IF

    CALL compute_propagator_matrix(qs_env,error)
 
    SELECT CASE(rtp_control%mat_exp)
    CASE( do_pade, do_taylor)
       CALL propagate_exp(rtp,rtp_control,error)
    CASE(do_arnoldi)
       CALL propagate_arnoldi(rtp,rtp_control,error)
    END SELECT
    CALL calc_update_rho(qs_env,error)

    IF(rtp_control%fixed_ions)my_force=.FALSE.
    CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,calculate_forces=my_force,&
         error=error)

    CALL step_finalize(qs_env,rtp_control,delta_mos,error)
    CALL cp_fm_vect_dealloc(delta_mos,error)

    CALL timestop(handle)

  END SUBROUTINE propagation_step

! *****************************************************************************
!> \brief Performes all the stuff to finish the step:
!>        convergence checks
!>        compying stuff into right place for the next step 
!>        updating the history for extrapolation 
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE step_finalize(qs_env,rtp_control,delta_mos,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(rtp_control_type), POINTER          :: rtp_control
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: delta_mos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'step_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ihist
    REAL(KIND=dp)                            :: delta_iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: s_mat
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, &
                                                mos_new, mos_old
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL get_qs_env(qs_env=qs_env,rtp=rtp,matrix_s=s_mat,error=error)
    CALL get_rtp(rtp=rtp,exp_H_old=exp_H_old,exp_H_new=exp_H_new,&
                 mos_old=mos_old,mos_new=mos_new,error=error)   
 
    IF(rtp_control%sc_check_start.LT.rtp%iter)&
        CALL rt_convergence(qs_env,delta_mos,delta_iter,error)
    rtp%converged=(delta_iter.LT.rtp_control%eps_ener)

    IF(rtp%converged)THEN
       DO i=1,SIZE(mos_new)
          CALL cp_fm_to_fm(mos_new(i)%matrix,mos_old(i)%matrix,error)
          IF(rtp_control%propagator==do_etrs)THEN
              CALL cp_fm_to_fm(exp_H_new(i)%matrix,exp_H_old(i)%matrix,error)
          ELSE
              CALL calc_SinvH(exp_H_old,qs_env,error)
          END IF
       END DO
       ihist=MOD(rtp%istep,rtp_control%aspc_order)+1
       IF(rtp_control%fixed_ions)THEN
          CALL  put_data_to_history(rtp,exp_H=exp_H_new,mos=mos_new,ihist=ihist,error=error)
       ELSE
          CALL  put_data_to_history(rtp,exp_H=exp_H_new,mos=mos_new,s_mat=s_mat,ihist=ihist,error=error)
       END IF
    END IF

  END SUBROUTINE step_finalize

! *****************************************************************************
!> \brief computes the propagator matrix for EM/ETRS, RTP/EMD
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE compute_propagator_matrix(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_propagator_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imat
    REAL(KIND=dp)                            :: dt, prefac
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, &
                                                propagator_matrix
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL get_qs_env(qs_env=qs_env,rtp=rtp,matrix_ks=matrix_ks,error=error)
    CALL get_rtp(rtp=rtp,exp_H_new=exp_H_new,exp_H_old=exp_H_old,&
                 propagator_matrix=propagator_matrix,dt=dt,error=error)

    prefac=-0.5_dp*dt
    SELECT CASE(qs_env%dft_control%rtp_control%propagator)
    ! HMMM, not really clean to put the matrix on the exp_H variable
    CASE(do_em)
       IF(rtp%iter.NE.1.OR.qs_env%dft_control%rtp_control%mat_exp==do_arnoldi) THEN
          CALL calc_SinvH(exp_H_new,qs_env,error)
          DO imat=1,SIZE(exp_H_new)
             CALL cp_fm_scale_and_add(0.0_dp,propagator_matrix(imat)%matrix,prefac,exp_H_old(imat)%matrix,error)
             CALL cp_fm_scale_and_add(1.0_dp,propagator_matrix(imat)%matrix,prefac,exp_H_new(imat)%matrix,error)
          END DO
       END IF
    CASE(do_etrs)
      IF(rtp%iter.NE.1.OR.qs_env%dft_control%rtp_control%mat_exp==do_arnoldi)THEN
          CALL calc_SinvH(exp_H_new,qs_env,error)
          DO imat=1,SIZE(exp_H_new)
             CALL cp_fm_scale_and_add(0.0_dp,propagator_matrix(imat)%matrix,prefac,exp_H_new(imat)%matrix,error)
          END DO
      END IF
       
    END SELECT

   END SUBROUTINE compute_propagator_matrix  

! *****************************************************************************
!> \brief computes t*S_inv*H, if needed t*Sinv*B
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calc_SinvH(exp_H,qs_env,error)
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_SinvH', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, im, ispin, ndim, re
    LOGICAL                                  :: failure
    REAL(dp)                                 :: t
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: SinvH
    TYPE(cp_fm_type), POINTER                :: B_mat, H_fm, S_inv, SinvB, &
                                                tmp_mat_H
    TYPE(rt_prop_type), POINTER              :: rtp

    failure=.FALSE.
    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,matrix_ks=matrix_ks,rtp=rtp,error=error)
    CALL get_rtp(rtp=rtp,S_inv=S_inv,dt=t,error=error)

    CALL cp_fm_create(H_fm,&
         matrix_struct=exp_H(1)%matrix%matrix_struct,&
         name="RTP_H_FM",&
         error=error)
    CALL cp_fm_create(tmp_mat_H,&
         matrix_struct=exp_H(1)%matrix%matrix_struct,&
         name="RTP_H_FM",&
         error=error)
    ndim=H_fm%matrix_struct%nrow_global

    DO ispin=1,SIZE(matrix_ks)

       re=ispin*2-1
       im=ispin*2

       CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,H_fm,error=error)
       CALL cp_fm_upper_to_full(H_fm,exp_H(im)%matrix,error)
       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,rtp%S_inv,H_fm,zero,exp_H(im)%matrix,error)

       IF(.NOT.qs_env%dft_control%rtp_control%fixed_ions)THEN
          CALL get_rtp(rtp=rtp,SinvH=SinvH,error=error)
          CALL cp_fm_to_fm(exp_H(im)%matrix,SinvH(ispin)%matrix,error)
       END IF
    END DO
    IF(.NOT.qs_env%dft_control%rtp_control%fixed_ions)THEN
       CALL get_rtp(rtp=rtp,B_mat=B_mat,SinvB=SinvB,error=error)
       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,rtp%S_inv,B_mat,zero,exp_H(1)%matrix,error)
       CALL cp_fm_to_fm(exp_H(1)%matrix,SinvB,error)
       IF(SIZE(matrix_ks)==2)CALL cp_fm_to_fm(exp_H(1)%matrix,exp_H(3)%matrix,error)
    END IF

    CALL cp_fm_release(H_fm,error)
    CALL cp_fm_release(tmp_mat_H,error)

    CALL timestop(handle)
  END SUBROUTINE calc_SinvH

! *****************************************************************************
!> \brief calculates the needed overlaplike matrices
!>        depending on the way the exponential is calculated, only S^-1 is needed
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE s_matrices_create (s_mat,rtp,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: s_mat
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 's_matrices_create', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, ndim
    LOGICAL                                  :: failure
    TYPE(cp_fm_type), POINTER                :: S_inv, tmp

    CALL timeset(routineN,handle)

    failure=.FALSE.

       CALL get_rtp(rtp=rtp,S_inv=S_inv,error=error)

       CALL cp_fm_create(tmp,&
            matrix_struct=S_inv%matrix_struct,&
            name="tmp_mat",&
            error=error)
       ndim=S_inv%matrix_struct%nrow_global


       CALL copy_dbcsr_to_fm(s_mat(1)%matrix,S_inv,error=error)
       CALL cp_fm_cholesky_decompose(S_inv,error=error)
       CALL cp_fm_cholesky_invert(S_inv,error=error)
       CALL cp_fm_upper_to_full(S_inv,tmp,error=error)
    CALL cp_fm_release(tmp,error)
    CALL timestop(handle)
  END SUBROUTINE s_matrices_create

! *****************************************************************************
!> \brief calculates the density from the complex MOs and passes the density to
!>        qs_env.
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calc_update_rho(qs_env,error)


    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_update_rho', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, im, re
    TYPE(cp_dbcsr_type), POINTER             :: tmp_rho_im, tmp_rho_re
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)

    NULLIFY (rho,ks_env,mos,rtp,tmp_rho_im, tmp_rho_re)
    CALL get_qs_env(qs_env=qs_env,ks_env=ks_env,rho=rho,mos=mos,rtp=rtp,error=error)
    CALL get_rtp(rtp=rtp,mos_new=mos_new,error=error)

    ALLOCATE(tmp_rho_im, tmp_rho_re)
    CALL cp_dbcsr_init(tmp_rho_re,error=error)
    CALL cp_dbcsr_init(tmp_rho_im,error=error)
    CALL cp_dbcsr_copy(tmp_rho_re,rho%rho_ao(1)%matrix,"tmp_tho_rts_re",error=error)
    CALL cp_dbcsr_copy(tmp_rho_im,rho%rho_ao(1)%matrix,"tmp_tho_rts_im",error=error)

    DO i=1,SIZE(mos)
       re=2*i-1
       im=2*i
       CALL cp_fm_to_fm(mos_new(im)%matrix,mos(i)%mo_set%mo_coeff,error)
       CALL calculate_density_matrix(mos(i)%mo_set,tmp_rho_im,error=error)
       CALL cp_fm_to_fm(mos_new(re)%matrix,mos(i)%mo_set%mo_coeff,error)
       CALL calculate_density_matrix(mos(i)%mo_set,tmp_rho_re,error=error)

       CALL cp_dbcsr_add(rho%rho_ao(i)%matrix,tmp_rho_re,&
            alpha_scalar=zero,beta_scalar=one,error=error)
       CALL cp_dbcsr_add(rho%rho_ao(i)%matrix,tmp_rho_im,&
            alpha_scalar=one,beta_scalar=one,error=error)
    END DO

    CALL qs_rho_update_rho(rho, qs_env, error=error)
    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)

    CALL cp_dbcsr_deallocate_matrix(tmp_rho_re,error)
    CALL cp_dbcsr_deallocate_matrix(tmp_rho_im,error)

    CALL timestop(handle)

  END SUBROUTINE calc_update_rho

  SUBROUTINE aspc_extrapolate(qs_env,rtp,aspc_order,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(rt_prop_type), POINTER              :: rtp
    INTEGER                                  :: aspc_order
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'aspc_extrapolate', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: cone = (1.0_dp,0.0_dp) , &
                                                czero = (0.0_dp,0.0_dp)
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, iaspc, icol_local, &
                                                ihist, imat, k, kdbl, n, &
                                                naspc, ncol_local, nmat
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_cfm_type), POINTER               :: cfm_tmp, cfm_tmp1, csc
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, s_hist
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, mos_new
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: exp_hist, mo_hist
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct, &
                                                matrix_struct_new
    TYPE(cp_fm_type), POINTER                :: fm_tmp, fm_tmp1, fm_tmp2

    CALL timeset(routineN,handle)
    NULLIFY(exp_hist)
    CALL get_rtp(rtp=rtp,exp_H_new=exp_H_new,mos_new=mos_new,error=error)
    nmat=SIZE(rtp%exp_H_new)
    CALL get_qs_env(qs_env,&
         matrix_s=matrix_s,&
         error=error)

    naspc=MIN(rtp%istep,aspc_order)
    CALL cite_reference(Kolafa2004)
    IF(rtp%history%mos_or_H==2)THEN
       exp_hist=>rtp%history%exp_history
       DO imat=1,nmat
          DO iaspc=1,naspc
             alpha=(-1.0_dp)**(iaspc + 1)*REAL(iaspc,KIND=dp)*&
                  binomial(2*naspc,naspc - iaspc)/binomial(2*naspc - 2,naspc -1)

             ihist=MOD(rtp%istep-iaspc,aspc_order)+1
             IF(iaspc==1)THEN
                CALL cp_fm_scale_and_add(zero,exp_H_new(imat)%matrix,alpha,exp_hist(imat,ihist)%matrix,error)
             ELSE
                CALL cp_fm_scale_and_add(one,exp_H_new(imat)%matrix,alpha,exp_hist(imat,ihist)%matrix,error)
             END IF
          END DO
       END DO
    END IF

    IF(rtp%history%mos_or_H==1)THEN !.AND.naspc.GT.1)THEN
       mo_hist=>rtp%history%mo_history
       DO imat=1,nmat
          DO iaspc=1,naspc
             alpha=(-1.0_dp)**(iaspc + 1)*REAL(iaspc,KIND=dp)*&
                  binomial(2*naspc,naspc - iaspc)/binomial(2*naspc - 2,naspc -1)
             ihist=MOD(rtp%istep-iaspc,aspc_order)+1
             IF(iaspc==1)THEN
                CALL cp_fm_scale_and_add(zero,mos_new(imat)%matrix,alpha,mo_hist(imat,ihist)%matrix,error)
             ELSE
                CALL cp_fm_scale_and_add(one,mos_new(imat)%matrix,alpha,mo_hist(imat,ihist)%matrix,error)
             END IF
          END DO
       END DO


       mo_hist=>rtp%history%mo_history
       s_hist=>rtp%history%s_history
       DO i=1,SIZE(mos_new)/2
          NULLIFY (matrix_struct,matrix_struct_new,csc,fm_tmp,fm_tmp1,fm_tmp2,cfm_tmp,cfm_tmp1)


          CALL cp_fm_struct_double(matrix_struct,&
                                mos_new(2*i)%matrix%matrix_struct,&
                                mos_new(2*i)%matrix%matrix_struct%context,&
                                .TRUE.,.FALSE.,error)

          CALL cp_fm_create(fm_tmp,matrix_struct,error=error)
          CALL cp_fm_create(fm_tmp1,matrix_struct,error=error)
          CALL cp_fm_create(fm_tmp2,mos_new(2*i)%matrix%matrix_struct,error=error)
          CALL cp_cfm_create(cfm_tmp,mos_new(2*i)%matrix%matrix_struct,error=error)
          CALL cp_cfm_create(cfm_tmp1,mos_new(2*i)%matrix%matrix_struct,error=error)

          CALL cp_fm_get_info(fm_tmp,&
                              ncol_global=kdbl,&
                              error=error)

          CALL cp_fm_get_info(mos_new(2*i)%matrix,&
                              nrow_global=n,&
                              ncol_global=k,&
                              ncol_local=ncol_local,&
                              error=error)

          CALL cp_fm_struct_create(matrix_struct_new,&
                                   template_fmstruct=matrix_struct,&
                                   nrow_global=k,&
                                   ncol_global=k,error=error)
          CALL cp_cfm_create(csc,matrix_struct_new,error=error)


          CALL cp_fm_struct_release(matrix_struct_new,error=error)
          CALL cp_fm_struct_release(matrix_struct,error=error)

          ! first the most recent


! reorthogonalize vectors

             DO icol_local=1,ncol_local
                fm_tmp%local_data(:,icol_local)=mos_new(2*i-1)%matrix%local_data(:,icol_local)
                fm_tmp%local_data(:,icol_local+ncol_local)=mos_new(2*i)%matrix%local_data(:,icol_local)
             END DO

             CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix,fm_tmp,fm_tmp1,kdbl,error=error)

             DO icol_local=1,ncol_local
                cfm_tmp%local_data(:,icol_local)=CMPLX(fm_tmp1%local_data(:,icol_local),&
                     fm_tmp1%local_data(:,icol_local+ncol_local),dp)
                cfm_tmp1%local_data(:,icol_local)=CMPLX(mos_new(2*i-1)%matrix%local_data(:,icol_local),&
                     mos_new(2*i)%matrix%local_data(:,icol_local),dp)
             END DO
             CALL cp_cfm_gemm('C','N',k,k,n,cone,cfm_tmp1,cfm_tmp,czero,csc,error=error)
             CALL cp_cfm_cholesky_decompose(csc,error=error)
             CALL cp_cfm_triangular_multiply(csc,cfm_tmp1,n_cols=k,side='R',invert_tr=.TRUE.,error=error)
             DO icol_local=1,ncol_local
                mos_new(2*i-1)%matrix%local_data(:,icol_local)=REAL(cfm_tmp1%local_data(:,icol_local),dp)
                mos_new(2*i)%matrix%local_data(:,icol_local)=AIMAG(cfm_tmp1%local_data(:,icol_local))
             END DO

! deallocate work matrices
             CALL cp_cfm_release(csc,error=error)
             CALL cp_fm_release(fm_tmp,error=error)
             CALL cp_fm_release(fm_tmp,error)
             CALL cp_fm_release(fm_tmp1,error)
             CALL cp_fm_release(fm_tmp2,error)
             CALL cp_cfm_release(cfm_tmp,error)
             CALL cp_cfm_release(cfm_tmp1,error)
          END DO

       END IF

    CALL timestop(handle)

  END SUBROUTINE aspc_extrapolate



  SUBROUTINE put_data_to_history(rtp,mos,exp_H,s_mat,ihist,error)
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos, exp_H
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: s_mat
    INTEGER                                  :: ihist
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'put_data_to_history', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF(rtp%history%mos_or_H==2)THEN
       DO i=1,SIZE(exp_H)
          CALL cp_fm_to_fm(exp_H(i)%matrix,rtp%history%exp_history(i,ihist)%matrix,error)
       END DO
    END IF

    IF(rtp%history%mos_or_H==1)THEN
       DO i=1,SIZE(mos)
          CALL cp_fm_to_fm(mos(i)%matrix,rtp%history%mo_history(i,ihist)%matrix,error)
       END DO

       IF(PRESENT(s_mat))THEN
          IF (ASSOCIATED(rtp%history%s_history(ihist)%matrix)) THEN ! the sparsity might be different
             ! (future struct:check)
             CALL cp_dbcsr_deallocate_matrix(rtp%history%s_history(ihist)%matrix,error=error)
          END IF
          ALLOCATE(rtp%history%s_history(ihist)%matrix)
          CALL cp_dbcsr_init(rtp%history%s_history(ihist)%matrix,error=error)
          CALL cp_dbcsr_copy(rtp%history%s_history(ihist)%matrix,s_mat(1)%matrix,error=error)
       END IF
    END IF

  END SUBROUTINE put_data_to_history
END MODULE rt_propagation_methods
