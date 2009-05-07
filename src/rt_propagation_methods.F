!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for propagating the orbitals
!> \author Florian Schiffmann (02.09)
! *****************************************************************************
MODULE rt_propagation_methods

  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type,&
                                             sm_pools_create_matrix_vect
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_cn,&
                                             do_em,&
                                             do_etrs,&
                                             do_diag, do_pade, do_taylor,do_arnoldi
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: binomial
  USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_environment_methods,          ONLY: qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_create,&
                                             qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: mo_set_p_type
  USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE rt_make_propagators,             ONLY: make_cn_propagator,&
                                             make_em_propagator,&
                                             make_etrs_exp,&
                                             propagate_arnoldi
  USE rt_matrix_exp,                   ONLY: backtransform_matrix
  USE rt_propagation_output,           ONLY: rt_convergence
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE sparse_matrix_types,             ONLY: add_matrices,&
                                             deallocate_matrix,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_methods'

  PUBLIC :: propagation_step,&
       s_matrices_create,&
       update_core_and_matrices,&
       calc_update_rho,&
       aspc_extrapolate,&
       calc_sinvH


CONTAINS

! *****************************************************************************
!> \brief performes a single propagation step a(t+Dt)=U(t+Dt,t)*a(0)
!>        and calculates the new exponential (at the moment for ETRS)
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE propagation_step(qs_env,delta_iter, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp)                                 :: delta_iter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'propagation_step', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: delta_mos, exp_H_new
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)
    NULLIFY(matrix_ks,rho,mos,ks_env,rtp)

    CALL get_qs_env(qs_env,&
         matrix_ks=matrix_ks,&
         rho=rho,&
         mos=mos,&
         matrix_s=matrix_s,&
         ks_env=ks_env,&
         rtp=rtp,&
         error=error)

    CALL get_rtp(rtp=rtp,exp_H_new=exp_H_new,error=error)

    ALLOCATE(delta_mos(2*SIZE(mos)))
    DO i=1,2*SIZE(mos)
       CALL cp_fm_create(delta_mos(i)%matrix,&
            matrix_struct=rtp%mos_new(i)%matrix%matrix_struct,&
            name="delta_mos"//TRIM(ADJUSTL(cp_to_string(i))),&
            error=error)
       CALL cp_fm_to_fm(rtp%mos_new(i)%matrix,delta_mos(i)%matrix,error)
    END DO
    
    SELECT CASE(qs_env%dft_control%rtp_control%mat_exp)
    CASE( do_diag, do_pade, do_taylor)
       
       SELECT CASE(qs_env%dft_control%rtp_control%propagator)
       CASE(do_cn,do_em)
          IF(rtp%iter.NE.1)&
               CALL calc_SinvH(exp_H_new,qs_env,error)
          CALL propagate_cn_or_em(qs_env,error)
       CASE(do_etrs)
          IF(rtp%iter.NE.1)&
               CALL make_etrs_exp(exp_H_new,qs_env,error)
          CALL propagate_etrs(qs_env,error)
       END SELECT
    CASE(do_arnoldi)
       IF(rtp%iter.NE.1)&
            CALL calc_SinvH(exp_H_new,qs_env,error)
       CALL propagate_arnoldi(qs_env,error)
       CALL calc_update_rho(qs_env,error)
    END SELECT

    DO i=1,2*SIZE(mos)
       CALL cp_fm_scale_and_add(-1.0_dp,delta_mos(i)%matrix,1.0_dp,rtp%mos_new(i)%matrix,error)
    END DO

    CALL rt_convergence(qs_env,delta_mos,delta_iter,error)
    CALL cp_fm_vect_dealloc(delta_mos,error)
    CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,calculate_forces=.FALSE.,&
         error=error)
    CALL timestop(handle)

  END SUBROUTINE propagation_step



! *****************************************************************************
!> \brief performes a crank nicholson or exponential midpoint propagation
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE propagate_cn_or_em(qs_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'propagate_cn_or_em', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: i, im_part, nao, nmo, &
                                                re_part, stat
    LOGICAL                                  :: failure
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, &
                                                mos_new, mos_old, propagator
    TYPE(cp_fm_type), POINTER                :: tmp_im, tmp_re
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rt_prop_type), POINTER              :: rtp

    failure=.FALSE.

    CALL get_qs_env(qs_env=qs_env,rho=rho,mos=mos,rtp=rtp,error=error)    
    CALL get_rtp(rtp=rtp,mos_old=mos_old,mos_new=mos_new,&
         exp_H_old=exp_H_old,exp_H_new=exp_H_new,error=error)

    ALLOCATE(propagator(2*SIZE(mos)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO i=1,2*SIZE(mos)
       CALL cp_fm_create(propagator(i)%matrix,&
            matrix_struct=exp_H_new(1)%matrix%matrix_struct,&
            error=error) 
    END DO

    SELECT CASE (qs_env%dft_control%rtp_control%propagator)
    CASE(do_cn)
       CALL make_cn_propagator(propagator,rtp,error)
    CASE(do_em)
       CALL make_em_propagator(propagator,rtp,qs_env,error)
    END SELECT

    DO i=1,SIZE(mos)
       re_part=2*i-1
       im_part=2*i
       
       tmp_re=>propagator(re_part)%matrix
       tmp_im=>propagator(im_part)%matrix
       
       CALL cp_fm_get_info(mos(i)%mo_set%mo_coeff,nrow_global=nao,ncol_global=nmo,error=error)
 

       CALL cp_fm_gemm("N","N",nao,nmo,nao,one,tmp_re,&
            mos_old(re_part)%matrix,zero,mos_new(re_part)%matrix,error)
       CALL cp_fm_gemm("N","N",nao,nmo,nao,one,tmp_im,&
            mos_old(re_part)%matrix,zero,mos_new(im_part)%matrix,error)

       CALL cp_fm_gemm("N","N",nao,nmo,nao,-one,tmp_im,&
            mos_old(im_part)%matrix,one,mos_new(re_part)%matrix,error)
       CALL cp_fm_gemm("N","N",nao,nmo,nao,one,tmp_re,&
            mos_old(im_part)%matrix,one,mos_new(im_part)%matrix,error)


    END DO

    CALL cp_fm_vect_dealloc(propagator,error)

    CALL calc_update_rho(qs_env,error)

  END SUBROUTINE propagate_cn_or_em



  
! *****************************************************************************
!> \brief applys the propagator on the MO coefficients
!>        at the moment only the ETRS, might have to be splitted when different
!>        propagators are implemented 
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE propagate_etrs(qs_env,error)


    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'propagate_etrs', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, im_part, nao, nmo, &
                                                re_part
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, &
                                                mos_new, mos_old
    TYPE(cp_fm_type), POINTER                :: tmp_im, tmp_re
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)
    NULLIFY (rho, ks_env, mos, rtp)
    CALL get_qs_env(qs_env=qs_env,ks_env=ks_env,rho=rho,mos=mos,rtp=rtp,error=error)    
    CALL get_rtp(rtp=rtp,mos_old=mos_old,mos_new=mos_new,&
         exp_H_old=exp_H_old,exp_H_new=exp_H_new,error=error)

    !case ETRS

    DO i=1,SIZE(mos)
       CALL cp_fm_get_info(mos(i)%mo_set%mo_coeff,nrow_global=nao,ncol_global=nmo,error=error)

       CALL cp_fm_create(tmp_re,&
            matrix_struct=exp_H_new(1)%matrix%matrix_struct,&
            name="tmp_re",&
            error=error)

       CALL cp_fm_create(tmp_im,&
            matrix_struct=exp_H_new(1)%matrix%matrix_struct,&
            name="tmp_im",&
            error=error)

       re_part=2*i-1
       im_part=2*i
       !whole bunch of gemms, since exp_H_old, exp_H_new and mo_coeff are all complex
       !case LSD needs 2 cycles
       CALL cp_fm_gemm("N","N",nao,nao,nao,one,exp_H_new(re_part)%matrix,&
            exp_H_old(re_part)%matrix,zero,tmp_re,error)
       CALL cp_fm_gemm("N","N",nao,nao,nao,one,exp_H_new(im_part)%matrix,&
            exp_H_old(re_part)%matrix,zero,tmp_im,error)

       CALL cp_fm_gemm("N","N",nao,nao,nao,-one,exp_H_new(im_part)%matrix,&
            exp_H_old(im_part)%matrix,one,tmp_re,error)
       CALL cp_fm_gemm("N","N",nao,nao,nao,one,exp_H_new(re_part)%matrix,&
            exp_H_old(im_part)%matrix,one,tmp_im,error)

       CALL cp_fm_gemm("N","N",nao,nmo,nao,one,tmp_re,&
            mos_old(re_part)%matrix,zero,mos_new(re_part)%matrix,error)
       CALL cp_fm_gemm("N","N",nao,nmo,nao,one,tmp_im,&
            mos_old(re_part)%matrix,zero,mos_new(im_part)%matrix,error)

       CALL cp_fm_gemm("N","N",nao,nmo,nao,-one,tmp_im,&
            mos_old(im_part)%matrix,one,mos_new(re_part)%matrix,error)
       CALL cp_fm_gemm("N","N",nao,nmo,nao,one,tmp_re,&
            mos_old(im_part)%matrix,one,mos_new(im_part)%matrix,error)

       CALL cp_fm_release(tmp_re,error)
       CALL cp_fm_release(tmp_im,error)

    END DO

    CALL calc_update_rho(qs_env,error)

    CALL timestop(handle)
  END SUBROUTINE propagate_etrs



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

    INTEGER                                  :: handle, ispin, lb, ndim, ub
    LOGICAL                                  :: failure
    REAL(dp)                                 :: t
    TYPE(cp_fm_type), POINTER                :: B_mat, H_fm, S_inv, tmp_mat_H
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
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
       
       lb=(ispin-1)*2+1
       ub=(ispin-1)*2+2

       CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,H_fm,error=error)
       CALL cp_fm_upper_to_full(H_fm,tmp_mat_H,error)
       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,rtp%S_inv,H_fm,zero,exp_H(ub)%matrix,error) 
       IF(.NOT.qs_env%dft_control%rtp_control%fixed_ions)THEN
          CALL get_rtp(rtp=rtp,B_mat=B_mat,error=error)
          CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,rtp%S_inv,B_mat,zero,exp_H(lb)%matrix,error)
       END IF
    END DO
            
    CALL cp_fm_release(H_fm,error)
    CALL cp_fm_release(tmp_mat_H,error)

    CALL timestop(handle)
  END SUBROUTINE calc_SinvH

       



! *****************************************************************************
!> \brief calculates the needed overlaplike matrices
!>        depending on the way the exponential is calculated, only S^-1
!>        or S^-0.5, S^0.5 and S^-1 are needed
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE s_matrices_create (s_mat,rtp,diagonalize,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: s_mat
    TYPE(rt_prop_type), POINTER              :: rtp
    LOGICAL                                  :: diagonalize
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 's_matrices_create', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, ndim, stat
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eigval_H
    TYPE(cp_fm_type), POINTER                :: eigvec_H, S_half, S_inv, &
                                                S_minus_half, tmp

    CALL timeset(routineN,handle)

    failure=.FALSE.

       CALL get_rtp(rtp=rtp,S_inv=S_inv,error=error)

       CALL cp_fm_create(tmp,&
            matrix_struct=S_inv%matrix_struct,&
            name="tmp_mat",&
            error=error)
       ndim=S_inv%matrix_struct%nrow_global

    IF(diagonalize)THEN
       CALL get_rtp(rtp=rtp,S_half=S_half,S_minus_half=S_minus_half,error=error)
       
       CALL cp_fm_create(eigvec_H,&
            matrix_struct=S_inv%matrix_struct,&
            name="tmp_EVEC",&
            error=error)
       
       !H_fm,eigvec_H,eigval_H is used as tmp, for diagonalizing S

       ALLOCATE(eigval_H(ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       CALL copy_sm_to_fm(s_mat(1)%matrix,tmp,error=error)
       CALL cp_fm_upper_to_full(tmp,eigvec_H,error)
       
       CALL cp_fm_syevd(tmp,eigvec_H,eigval_H,error)
       
       eigval_H(:)=one/eigval_H(:)
       CALL backtransform_matrix(eigval_H,eigvec_H,S_inv,error)       
       eigval_H(:)=SQRT(eigval_H(:))
       CALL backtransform_matrix(eigval_H,eigvec_H,S_minus_half,error)
       eigval_H(:)=one/eigval_H(:)
       CALL backtransform_matrix(eigval_H,eigvec_H,S_half,error)
       CALL cp_fm_release(eigvec_H,error)       

    ELSE
       CALL copy_sm_to_fm(s_mat(1)%matrix,S_inv,error=error)
       CALL cp_fm_cholesky_decompose(S_inv,error=error)
       CALL cp_fm_cholesky_invert(S_inv,error=error)
       CALL cp_fm_upper_to_full(S_inv,tmp,error=error)
    END IF
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

    INTEGER                                  :: handle, i, im_part, re_part
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_type), POINTER          :: tmp_rho_im, tmp_rho_re
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)

    NULLIFY (rho,ks_env,mos,tmp_rho_re,tmp_rho_im,rtp)
    CALL get_qs_env(qs_env=qs_env,ks_env=ks_env,rho=rho,mos=mos,rtp=rtp,error=error)    
    CALL get_rtp(rtp=rtp,mos_new=mos_new,error=error)

    CALL replicate_matrix_structure(rho%rho_ao(1)%matrix,&
         tmp_rho_re,&
         "tmp_tho_rts_re",&
         error=error)

    CALL replicate_matrix_structure(rho%rho_ao(1)%matrix,&
         tmp_rho_im,&
         "tmp_tho_rts_im",&
         error=error)

    DO i=1,SIZE(mos)
       re_part=2*(i-1)+1       
       im_part=2*(i-1)+2   
       CALL cp_fm_to_fm(mos_new(im_part)%matrix,mos(i)%mo_set%mo_coeff,error)
       CALL calculate_density_matrix(mos(i)%mo_set,tmp_rho_im,error)    
       CALL cp_fm_to_fm(mos_new(re_part)%matrix,mos(i)%mo_set%mo_coeff,error)
       CALL calculate_density_matrix(mos(i)%mo_set,tmp_rho_re,error)

       CALL add_matrices(rho%rho_ao(i)%matrix,one,tmp_rho_re,one,tmp_rho_im)
    END DO

    CALL qs_rho_update_rho(rho, qs_env, error)
    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)

    CALL deallocate_matrix(tmp_rho_re,error)
    CALL deallocate_matrix(tmp_rho_im,error)

    CALL timestop(handle)

  END SUBROUTINE calc_update_rho

! *****************************************************************************
!> \brief needed only for EMD. Recalculates neighborlists and integrals, when 
!>        ionic positions change
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE update_core_and_matrices(qs_env,globenv,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_core_and_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)

    para_env=>qs_env%para_env

    CALL get_qs_env(qs_env=qs_env,input=input, dft_control=dft_control,error=error)

    CALL build_qs_neighbor_lists(qs_env,para_env,force_env_section=input,error=error)
    CALL build_core_hamiltonian_matrix(qs_env=qs_env, globenv=globenv,&
         calculate_forces=.FALSE.,error=error)
    CALL qs_env_update_s_mstruct(qs_env,error=error)    
    
    CALL calculate_ecore_self(qs_env,error=error)
    CALL calculate_ecore_overlap(qs_env, para_env, &
         calculate_forces=.FALSE.,error=error)
    CALL get_qs_env(qs_env=qs_env,matrix_ks=matrix_ks,ks_env=ks_env,error=error)
    CALL mpools_get(qs_env%mpools,&
         S_sm_pools=S_sm_pools,&
         error=error)
    
    IF (.NOT.ASSOCIATED(matrix_ks)) THEN
       CALL sm_pools_create_matrix_vect(S_sm_pools,matrix_ks,&
            name="KOHN-SHAM_MATRIX",&
            error=error)
       CALL set_qs_env(qs_env=qs_env,&
            matrix_ks=matrix_ks,error=error)
    END IF

    IF (.NOT.ASSOCIATED(ks_env)) THEN
       CALL qs_ks_create(ks_env,qs_env=qs_env,error=error)
       CALL set_qs_env(qs_env, ks_env=ks_env,error=error)
       CALL qs_ks_release(ks_env,error=error)
    END IF
   
    
    CALL timestop(handle)
  
  END SUBROUTINE update_core_and_matrices

  SUBROUTINE aspc_extrapolate(rtp,aspc_order,error)
    TYPE(rt_prop_type), POINTER              :: rtp
    INTEGER                                  :: aspc_order
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'aspc_extrapolate', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: iaspc, ihist, imat, naspc, &
                                                nmat
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: exp_hist

    NULLIFY(exp_hist)
    CALL get_rtp(rtp=rtp,exp_H_new=exp_H_new,error=error)
    nmat=SIZE(rtp%exp_H_new)
    exp_hist=>rtp%history%exp_history

    naspc=MIN(rtp%istep,aspc_order)

    IF(naspc.GT.1)THEN
       DO imat=1,nmat
          DO iaspc=1,naspc
             IF(iaspc==1)alpha = REAL(4*naspc - 2,KIND=dp)/REAL(naspc + 1,KIND=dp)
             IF(iaspc.GT.1) alpha=(-1.0_dp)**(iaspc + 1)*REAL(iaspc,KIND=dp)*&
                  binomial(2*naspc,naspc - iaspc)/binomial(2*naspc - 2,naspc -1)
             
             ihist=MOD(rtp%istep-iaspc,aspc_order)+1
             IF(iaspc==1)THEN
                CALL cp_fm_scale_and_add(zero,exp_H_new(imat)%matrix,alpha,exp_hist(imat,ihist)%matrix,error)
             ELSE
                CALL cp_fm_scale_and_add(one,exp_H_new(imat)%matrix,alpha,exp_hist(imat,ihist)%matrix,error)
             END IF
          END DO
       END DO
    ELSE IF(naspc==1)THEN
       DO imat=1,nmat
          CALL cp_fm_scale_and_add(zero,exp_H_new(imat)%matrix,one,exp_hist(imat,1)%matrix,error)
       END DO
    END IF

  END SUBROUTINE aspc_extrapolate


END MODULE rt_propagation_methods
