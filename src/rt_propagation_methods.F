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
  USE cp_control_types,                ONLY: dft_control_type,&
                                             rtp_control_type
  USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                             cp_dbcsr_cholesky_invert
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, &
       cp_dbcsr_copy_into_existing, cp_dbcsr_create, &
       cp_dbcsr_deallocate_matrix, cp_dbcsr_deallocate_matrix_set, &
       cp_dbcsr_desymmetrize, cp_dbcsr_filter, cp_dbcsr_frobenius_norm, &
       cp_dbcsr_get_block_p, cp_dbcsr_init, cp_dbcsr_init_p, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_transposed, &
       cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
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
                                             do_bch,&
                                             do_em,&
                                             do_etrs,&
                                             do_pade,&
                                             do_taylor
  USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: binomial
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE rt_make_propagators,             ONLY: propagate_arnoldi,&
                                             propagate_bch,&
                                             propagate_exp,&
                                             propagate_exp_density
  USE rt_matrix_exp_dbcsr,             ONLY: cp_complex_dbcsr_gemm_3
  USE rt_propagation_output,           ONLY: report_density_occupation,&
                                             rt_convergence,&
                                             rt_convergence_density
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             get_rtp_hfx,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_methods'

  PUBLIC :: propagation_step,&
            s_matrices_create,&
            calc_update_rho,&
            calc_update_rho_sparse,&
            calc_sinvH,&
            put_data_to_history


CONTAINS

! *****************************************************************************
!> \brief performes a single propagation step a(t+Dt)=U(t+Dt,t)*a(0)
!>        and calculates the new exponential
!> \param qs_env ...
!> \param rtp ...
!> \param rtp_control ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE propagation_step(qs_env, rtp, rtp_control, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'propagation_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: aspc_order, handle, i
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: delta_P, rho_new
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: delta_mos, mos_new
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)
    NULLIFY(rho,ks_env,delta_P, rho_new,delta_mos, mos_new)
    ! get everything needed and set some values
    CALL get_qs_env(qs_env,&
         rho=rho,&
         ks_env=ks_env,&
         error=error)
    IF(rtp%iter==1)THEN
       aspc_order=rtp_control%aspc_order
       CALL aspc_extrapolate(qs_env,rtp,aspc_order,error=error)
       IF(rtp%linear_scaling) THEN
          CALL calc_update_rho_sparse(qs_env,error=error)
       ELSE
          CALL calc_update_rho(qs_env,error=error)
       ENDIF
       CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., error=error)
    END IF

    rtp%delta_iter=1.0_dp
    IF(rtp%linear_scaling) THEN
       ! keep temporary copy of the starting density matrix to check for convergence  
       CALL get_rtp(rtp=rtp,rho_new=rho_new,error=error)
       NULLIFY(delta_P)
       CALL cp_dbcsr_allocate_matrix_set(delta_P,SIZE(rho_new),error=error)
       DO i=1,SIZE(rho_new)
          CALL cp_dbcsr_init_p(delta_P(i)%matrix,error=error)
          CALL cp_dbcsr_create(delta_P(i)%matrix,template=rho_new(i)%matrix,error=error)
          CALL cp_dbcsr_copy(delta_P(i)%matrix,rho_new(i)%matrix,error=error)
       END DO
    ELSE
       ! keep temporary copy of the starting mos to check for convergence   
       CALL get_rtp(rtp=rtp, mos_new=mos_new, error=error)
       ALLOCATE(delta_mos(SIZE(mos_new)))
       DO i=1,SIZE(mos_new)
          CALL cp_fm_create(delta_mos(i)%matrix,&
               matrix_struct=mos_new(i)%matrix%matrix_struct,&
               name="delta_mos"//TRIM(ADJUSTL(cp_to_string(i))),&
               error=error)
          CALL cp_fm_to_fm(mos_new(i)%matrix,delta_mos(i)%matrix,error)
       END DO
    ENDIF

    CALL compute_propagator_matrix(qs_env,error)

    SELECT CASE(rtp_control%mat_exp)
    CASE( do_pade, do_taylor)
       IF(rtp%linear_scaling) THEN
         CALL propagate_exp_density(rtp,rtp_control,error)
         CALL calc_update_rho_sparse(qs_env,error)
       ELSE
          CALL propagate_exp(rtp,rtp_control,error)
          CALL calc_update_rho(qs_env,error)
       END IF
    CASE(do_arnoldi)
       CALL propagate_arnoldi(rtp,rtp_control,error)
       CALL calc_update_rho(qs_env,error)
    CASE(do_bch)
       CALL propagate_bch(rtp,rtp_control,error)
       CALL calc_update_rho_sparse(qs_env,error)
    END SELECT
    CALL step_finalize(qs_env,rtp_control,delta_mos,delta_P,error)
    IF(rtp%linear_scaling) THEN
       CALL cp_dbcsr_deallocate_matrix_set(delta_P,error)
    ELSE
       CALL cp_fm_vect_dealloc(delta_mos,error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE propagation_step

! *****************************************************************************
!> \brief Performes all the stuff to finish the step:
!>        convergence checks
!>        copying stuff into right place for the next step 
!>        updating the history for extrapolation 
!> \param qs_env ...
!> \param rtp_control ...
!> \param delta_mos ...
!> \param delta_P ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE step_finalize(qs_env,rtp_control,delta_mos,delta_P,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(rtp_control_type), POINTER          :: rtp_control
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: delta_mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: delta_P
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'step_finalize', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, ihist
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, &
                                                rho_new, rho_old, s_mat
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, mos_old
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,rtp=rtp,matrix_s=s_mat,error=error)
    CALL get_rtp(rtp=rtp,exp_H_old=exp_H_old,exp_H_new=exp_H_new,error=error)

    IF(rtp_control%sc_check_start.LT.rtp%iter) THEN
        IF(rtp%linear_scaling) THEN
           CALL rt_convergence_density(qs_env,delta_P,rtp%delta_iter,error)
        ELSE
           CALL rt_convergence(qs_env,delta_mos,rtp%delta_iter,error)
        END IF
        rtp%converged=(rtp%delta_iter.LT.rtp_control%eps_ener)
    END IF
    IF(rtp%converged)THEN
       IF(rtp%linear_scaling) THEN
          CALL get_rtp(rtp=rtp,rho_old=rho_old,rho_new=rho_new,error=error) 
          IF(rtp_control%mcweeny_max_iter>0.OR.rtp_control%calc_idempotency) THEN
             IF(SIZE(rho_new)==2) THEN
                CALL cp_dbcsr_scale(rho_new(1)%matrix,one/2,error=error)
                CALL cp_dbcsr_scale(rho_new(2)%matrix,one/2,error=error)
             ENDIF
             IF(rtp_control%orthonormal) THEN
                CALL purify_mcweeny_complex_orth(rho_new,s_mat,rtp%filter_eps,rtp%filter_eps_small,&
                     rtp_control%mcweeny_max_iter,rtp_control%mcweeny_eps,error)
             ELSE
                CALL purify_mcweeny_complex_nonorth(rho_new,s_mat,rtp%filter_eps,rtp%filter_eps_small,&
                     rtp_control%mcweeny_max_iter,rtp_control%mcweeny_eps,error)
             ENDIF
             IF(SIZE(rho_new)==2) THEN
                CALL cp_dbcsr_scale(rho_new(1)%matrix,one*2,error=error)
                CALL cp_dbcsr_scale(rho_new(2)%matrix,one*2,error=error)
             ENDIF
             CALL calc_update_rho_sparse(qs_env,error)
          ENDIF
          CALL report_density_occupation(rtp,rho_new,error=error)
          DO i=1,SIZE(rho_new)
             CALL cp_dbcsr_copy(rho_old(i)%matrix,rho_new(i)%matrix,error=error)
          END DO
       ELSE
          CALL get_rtp(rtp=rtp,mos_old=mos_old,mos_new=mos_new,error=error)
          DO i=1,SIZE(mos_new)
             CALL cp_fm_to_fm(mos_new(i)%matrix,mos_old(i)%matrix,error)
          END DO
       ENDIF
       IF(rtp_control%propagator==do_etrs)THEN
          DO i=1,SIZE(exp_H_new)
             CALL cp_dbcsr_copy(exp_H_old(i)%matrix,exp_H_new(i)%matrix,error=error)
          END DO
       ELSE
          CALL calc_SinvH(exp_H_old,qs_env,error)
       END IF
       ihist=MOD(rtp%istep,rtp_control%aspc_order)+1
       IF(rtp_control%fixed_ions)THEN
          CALL  put_data_to_history(rtp,rho=rho_new,mos=mos_new,ihist=ihist,error=error)
       ELSE
          CALL  put_data_to_history(rtp,rho=rho_new,mos=mos_new,s_mat=s_mat,ihist=ihist,error=error)
       END IF
    END IF
    CALL timestop(handle)

  END SUBROUTINE step_finalize

! *****************************************************************************
!> \brief computes the propagator matrix for EM/ETRS, RTP/EMD
!> \param qs_env ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE compute_propagator_matrix(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_propagator_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    REAL(Kind=dp)                            :: dt, prefac
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, &
                                                propagator_matrix
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control

    NULLIFY(dft_control)
    CALL timeset(routineN,handle)
    CALL get_qs_env(qs_env, rtp=rtp, dft_control=dft_control, error=error)
    CALL get_rtp(rtp=rtp,exp_H_new=exp_H_new,exp_H_old=exp_H_old,&
                 propagator_matrix=propagator_matrix,dt=dt,error=error)
    rtp_control=>dft_control%rtp_control

    prefac=-0.5_dp*dt

    CALL calc_SinvH(exp_H_new,qs_env,error)
    DO i=1,SIZE(exp_H_new)
       CALL cp_dbcsr_add(propagator_matrix(i)%matrix,exp_H_new(i)%matrix,0.0_dp,prefac,error)
       IF(rtp_control%propagator==do_em)&
          CALL cp_dbcsr_add(propagator_matrix(i)%matrix,exp_H_old(i)%matrix,1.0_dp,prefac,error)
    END DO

    CALL timestop(handle)

   END SUBROUTINE compute_propagator_matrix  

! *****************************************************************************
!> \brief computes t*S_inv*H, if needed t*Sinv*B
!> \param exp_H ...
!> \param qs_env ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calc_SinvH(exp_H,qs_env,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: exp_H
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_SinvH', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, im, ispin, re
    REAL(dp)                                 :: t
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, SinvB, SinvH
    TYPE(cp_dbcsr_type)                      :: matrix_ks_nosym, tmp
    TYPE(cp_dbcsr_type), POINTER             :: B_mat, S_inv, S_minus_half
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control

    NULLIFY(dft_control)
    CALL timeset(routineN,handle)
    CALL get_qs_env(qs_env,&
                    matrix_ks=matrix_ks,&
                    rtp=rtp,&
                    dft_control=dft_control,&
                    error=error)
    rtp_control=>dft_control%rtp_control
    CALL get_rtp(rtp=rtp,S_inv=S_inv,S_minus_half=S_minus_half,dt=t,error=error)
    CALL cp_dbcsr_init(matrix_ks_nosym,error=error)
    CALL cp_dbcsr_create(matrix_ks_nosym,template=matrix_ks(1)%matrix,matrix_type="N",error=error)
    CALL cp_dbcsr_init(tmp,error=error)
    CALL cp_dbcsr_create(tmp,template=matrix_ks(1)%matrix,matrix_type="N",error=error)
    DO ispin=1,SIZE(matrix_ks)
       re=ispin*2-1
       im=ispin*2
       CALL cp_dbcsr_desymmetrize(matrix_ks(ispin)%matrix,matrix_ks_nosym,error=error)
       IF(rtp_control%orthonormal) THEN
          CALL cp_dbcsr_multiply("N","N",one,S_minus_half,matrix_ks_nosym,zero,tmp,&
               filter_eps=rtp%filter_eps,error=error)
          CALL cp_dbcsr_multiply("N","N",one,tmp,S_minus_half,zero,exp_H(im)%matrix,&
               filter_eps=rtp%filter_eps,error=error)
       ELSE
          CALL cp_dbcsr_multiply("N","N",one,S_inv,matrix_ks_nosym,zero,exp_H(im)%matrix,&
               filter_eps=rtp%filter_eps,error=error)
       ENDIF
       IF(.NOT.dft_control%rtp_control%fixed_ions)THEN
          CALL get_rtp(rtp=rtp,SinvH=SinvH,error=error)
          CALL cp_dbcsr_copy(SinvH(ispin)%matrix,exp_H(im)%matrix,error=error)
       END IF
    END DO
    IF(.NOT.dft_control%rtp_control%fixed_ions.OR.rtp%do_hfx)THEN
       CALL get_rtp(rtp=rtp,B_mat=B_mat,SinvB=SinvB,error=error)
       IF(rtp%do_hfx)THEN
          DO ispin=1,SIZE(matrix_ks)
             re=ispin*2-1
             im=ispin*2           
             CALL cp_dbcsr_set(matrix_ks_nosym,0.0_dp,error)
             CALL cp_dbcsr_desymmetrize (rtp%rtp_hfx%matrix_ks_im(ispin)%matrix, matrix_ks_nosym,&
                     error=error)

             ! take care of the EMD case and add the velocity scaled S_derivative
             IF(.NOT.dft_control%rtp_control%fixed_ions)&
                CALL cp_dbcsr_add(matrix_ks_nosym,B_mat,1.0_dp,-1.0_dp,error=error)
            
             IF(rtp_control%orthonormal) THEN
                CALL cp_dbcsr_multiply("N","N",-one,S_minus_half,matrix_ks_nosym,zero,tmp,&
                     filter_eps=rtp%filter_eps,error=error)
                CALL cp_dbcsr_multiply("N","N",one,tmp,S_minus_half,zero,exp_H(re)%matrix,&
                     filter_eps=rtp%filter_eps,error=error)
             ELSE
                CALL cp_dbcsr_multiply("N","N",-one,S_inv,matrix_ks_nosym,zero,exp_H(re)%matrix,&
                     filter_eps=rtp%filter_eps,error=error)
             ENDIF

             IF(.NOT.dft_control%rtp_control%fixed_ions)&
                CALL cp_dbcsr_copy(SinvB(ispin)%matrix,exp_H(re)%matrix,error=error)
          END DO
       ELSE
          ! in case of pure EMD its only needed once as B is the same for both spins
          IF(rtp_control%orthonormal) THEN
             CALL cp_dbcsr_multiply("N","N",one,S_minus_half,B_mat,zero,tmp,&
                  filter_eps=rtp%filter_eps,error=error)
             CALL cp_dbcsr_multiply("N","N",one,tmp,S_minus_half,zero,exp_H(1)%matrix,&
                  filter_eps=rtp%filter_eps,error=error)
          ELSE
             CALL cp_dbcsr_multiply("N","N",one,S_inv,B_mat,zero,exp_H(1)%matrix,filter_eps=rtp%filter_eps,error=error)
          ENDIF

          CALL cp_dbcsr_copy(SinvB(1)%matrix,exp_H(1)%matrix,error=error)

          IF(SIZE(matrix_ks)==2)CALL cp_dbcsr_copy(exp_H(3)%matrix,exp_H(1)%matrix,error=error)
          IF(SIZE(matrix_ks)==2)CALL cp_dbcsr_copy(SinvB(2)%matrix,SinvB(1)%matrix,error=error)
       END IF
    ELSE
       !set real part to zero
       DO ispin=1,SIZE(exp_H)/2
          re=ispin*2-1
          im=ispin*2    
          CALL cp_dbcsr_set(exp_H(re)%matrix,zero,error=error)
       ENDDO
    END IF
    CALL cp_dbcsr_release(matrix_ks_nosym,error)
    CALL cp_dbcsr_release(tmp,error)
    CALL timestop(handle)
  END SUBROUTINE calc_SinvH

! *****************************************************************************
!> \brief calculates the needed overlaplike matrices
!>        depending on the way the exponential is calculated, only S^-1 is needed
!> \param s_mat ...
!> \param rtp ...
!> \param error ...
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

    INTEGER                                  :: handle
    TYPE(cp_dbcsr_type), POINTER             :: S_half, S_inv, S_minus_half

    CALL timeset(routineN,handle)

    CALL get_rtp(rtp=rtp,S_inv=S_inv,error=error)

    IF(rtp%linear_scaling) THEN
       CALL get_rtp(rtp=rtp,S_half=S_half,S_minus_half=S_minus_half,error=error)
       CALL matrix_sqrt_Newton_Schulz(S_half,S_minus_half,s_mat(1)%matrix,rtp%filter_eps,&
            rtp%newton_schulz_order,rtp%lanzcos_threshold,rtp%lanzcos_max_iter,error=error)
       CALL cp_dbcsr_multiply("N","N",one,S_minus_half,S_minus_half,zero,S_inv,&
            filter_eps=rtp%filter_eps,error=error)
    ELSE
       CALL cp_dbcsr_copy(S_inv,s_mat(1)%matrix,error=error)
       CALL cp_dbcsr_cholesky_decompose(S_inv,para_env=rtp%ao_ao_fmstruct%para_env,&
            blacs_env=rtp%ao_ao_fmstruct%context,error=error)
       CALL cp_dbcsr_cholesky_invert(S_inv,para_env=rtp%ao_ao_fmstruct%para_env,&
            blacs_env=rtp%ao_ao_fmstruct%context,upper_to_full=.TRUE.,error=error)
    ENDIF

    CALL timestop(handle)
  END SUBROUTINE s_matrices_create

! *****************************************************************************
!> \brief Calculates the frobenius norm of a omplex matrix represented by two real matrices
!> \param frob_norm ...
!> \param mat_re ...
!> \param mat_im ...
!> \param error ...
!> \author Samuel Andermatt (04.14)
! *****************************************************************************

  SUBROUTINE complex_frobenius_norm(frob_norm,mat_re,mat_im,error)

    REAL(KIND=dp), INTENT(out)               :: frob_norm
    TYPE(cp_dbcsr_type), POINTER             :: mat_re, mat_im
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'complex_frobenius_norm', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: col_atom, handle, row_atom
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:), POINTER          :: block_values, block_values2
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type), POINTER             :: tmp

    CALL timeset(routineN,handle)


    NULLIFY(tmp)
    ALLOCATE(tmp)
    CALL cp_dbcsr_init(tmp,error=error)
    CALL cp_dbcsr_create(tmp,template=mat_re,error=error)
    !make sure the tmp has the same sparsity pattern as the real and the complex part combined
    CALL cp_dbcsr_add(tmp,mat_re,zero,one,error=error)
    CALL cp_dbcsr_add(tmp,mat_im,zero,one,error=error)
    CALL cp_dbcsr_set(tmp,zero,error=error)
    !calculate the hadamard product
    CALL cp_dbcsr_iterator_start(iter, tmp)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
       CALL cp_dbcsr_get_block_p(mat_re, row_atom, col_atom, block_values2, found=found)
       IF(found) THEN
          block_values=block_values2*block_values2
       ENDIF
       CALL cp_dbcsr_get_block_p(mat_im, row_atom, col_atom, block_values2, found=found)
       IF(found) THEN
          block_values=block_values+block_values2*block_values2
       ENDIF
       block_values=SQRT(block_values)
    END DO
    CALL cp_dbcsr_iterator_stop (iter)
    frob_norm=cp_dbcsr_frobenius_norm(tmp)

    CALL cp_dbcsr_deallocate_matrix(tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE complex_frobenius_norm

! *****************************************************************************
!> \brief Does McWeeny for complex matrices in the non-orthogonal basis
!> \param P ...
!> \param s_mat ...
!> \param eps ...
!> \param eps_small ...
!> \param max_iter ...
!> \param threshold ...
!> \param error ...
!> \author Samuel Andermatt (04.14)
! *****************************************************************************

  SUBROUTINE purify_mcweeny_complex_nonorth(P,s_mat,eps,eps_small,max_iter,threshold,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: P, s_mat
    REAL(KIND=dp), INTENT(in)                :: eps, eps_small
    INTEGER, INTENT(in)                      :: max_iter
    REAL(KIND=dp), INTENT(in)                :: threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'purify_mcweeny_complex_nonorth', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, im, imax, ispin, &
                                                re, unit_nr
    REAL(KIND=dp)                            :: frob_norm
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: PS, PSP, tmp
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger   => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    NULLIFY(tmp,PS,PSP)
    CALL cp_dbcsr_allocate_matrix_set(tmp,SIZE(P),error=error)
    CALL cp_dbcsr_allocate_matrix_set(PSP,SIZE(P),error=error)
    CALL cp_dbcsr_allocate_matrix_set(PS,SIZE(P),error=error)
    DO i=1,SIZE(P)
       CALL cp_dbcsr_init_p(PS(i)%matrix,error=error)
       CALL cp_dbcsr_create(PS(i)%matrix,template=P(1)%matrix,error=error)
       CALL cp_dbcsr_init_p(PSP(i)%matrix,error=error)
       CALL cp_dbcsr_create(PSP(i)%matrix,template=P(1)%matrix,error=error)
       CALL cp_dbcsr_init_p(tmp(i)%matrix,error=error)
       CALL cp_dbcsr_create(tmp(i)%matrix,template=P(1)%matrix,error=error)
    ENDDO
    DO ispin=1,SIZE(P)/2
       re=2*ispin-1
       im=2*ispin
       imax=MAX(max_iter,1) !if max_iter is 0 then only the deviation from idempotency needs to be calculated
       DO i=1,imax
          CALL cp_dbcsr_multiply("N", "N", one, P(re)%matrix,s_mat(1)%matrix,&
               zero, PS(re)%matrix, filter_eps=eps_small,error=error)
          CALL cp_dbcsr_multiply("N", "N", one, P(im)%matrix,s_mat(1)%matrix,&
               zero, PS(im)%matrix, filter_eps=eps_small,error=error)
          CALL cp_complex_dbcsr_gemm_3("N","N",one,PS(re)%matrix,PS(im)%matrix,&
               P(re)%matrix,P(im)%matrix,zero,PSP(re)%matrix,PSP(im)%matrix,&
               filter_eps=eps_small,error=error)
          CALL cp_dbcsr_copy(tmp(re)%matrix,PSP(re)%matrix,error=error)
          CALL cp_dbcsr_copy(tmp(im)%matrix,PSP(im)%matrix,error=error)
          CALL cp_dbcsr_add(tmp(re)%matrix,P(re)%matrix,1.0_dp,-1.0_dp,error=error)
          CALL cp_dbcsr_add(tmp(im)%matrix,P(im)%matrix,1.0_dp,-1.0_dp,error=error)
          CALL complex_frobenius_norm(frob_norm,tmp(re)%matrix,tmp(im)%matrix,error=error)
          IF(unit_nr.gt.0) WRITE(unit_nr,'(t3,a,2f16.8)')"Deviation from idempotency: ",frob_norm
          IF(frob_norm.GT.threshold.AND.max_iter>0)THEN
             CALL cp_dbcsr_copy(P(re)%matrix,PSP(re)%matrix,error=error)
             CALL cp_dbcsr_copy(P(im)%matrix,PSP(im)%matrix,error=error)
             CALL cp_complex_dbcsr_gemm_3("N", "N", -2.0_dp, PS(re)%matrix,PS(im)%matrix,&
                  PSP(re)%matrix,PSP(im)%matrix,3.0_dp,P(re)%matrix,P(im)%matrix,&
                  filter_eps=eps_small,error=error)
             CALL cp_dbcsr_filter(P(re)%matrix,eps,error=error)
             CALL cp_dbcsr_filter(P(im)%matrix,eps,error=error)
             !make sure P is exactly hermitian
             CALL cp_dbcsr_transposed(tmp(re)%matrix,P(re)%matrix,error=error)
             CALL cp_dbcsr_add(P(re)%matrix,tmp(re)%matrix,one/2,one/2,error=error)
             CALL cp_dbcsr_transposed(tmp(im)%matrix,P(im)%matrix,error=error)
             CALL cp_dbcsr_add(P(im)%matrix,tmp(im)%matrix,one/2,-one/2,error=error)
          ELSE
              EXIT
          END IF
       END DO
       !make sure P is hermitian
       CALL cp_dbcsr_transposed(tmp(re)%matrix,P(re)%matrix,error=error)
       CALL cp_dbcsr_add(P(re)%matrix,tmp(re)%matrix,one/2,one/2,error=error)
       CALL cp_dbcsr_transposed(tmp(im)%matrix,P(im)%matrix,error=error)
       CALL cp_dbcsr_add(P(im)%matrix,tmp(im)%matrix,one/2,-one/2,error=error)
    END DO
    CALL cp_dbcsr_deallocate_matrix_set(tmp,error)
    CALL cp_dbcsr_deallocate_matrix_set(PS,error)
    CALL cp_dbcsr_deallocate_matrix_set(PSP,error)

    CALL timestop(handle)

  END SUBROUTINE purify_mcweeny_complex_nonorth

! *****************************************************************************
!> \brief Does McWeeny for complex matrices in the orthonormal basis
!> \param P ...
!> \param s_mat ...
!> \param eps ...
!> \param eps_small ...
!> \param max_iter ...
!> \param threshold ...
!> \param error ...
!> \author Samuel Andermatt (04.14)
! *****************************************************************************

  SUBROUTINE purify_mcweeny_complex_orth(P,s_mat,eps,eps_small,max_iter,threshold,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: P, s_mat
    REAL(KIND=dp), INTENT(in)                :: eps, eps_small
    INTEGER, INTENT(in)                      :: max_iter
    REAL(KIND=dp), INTENT(in)                :: threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'purify_mcweeny_complex_orth', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, im, imax, ispin, &
                                                re, unit_nr
    REAL(KIND=dp)                            :: frob_norm
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: PP, tmp
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger   => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    NULLIFY(tmp,PP)
    CALL cp_dbcsr_allocate_matrix_set(tmp,SIZE(P),error=error)
    CALL cp_dbcsr_allocate_matrix_set(PP,SIZE(P),error=error)
    DO i=1,SIZE(P)
       CALL cp_dbcsr_init_p(PP(i)%matrix,error=error)
       CALL cp_dbcsr_create(PP(i)%matrix,template=P(1)%matrix,error=error)
       CALL cp_dbcsr_init_p(tmp(i)%matrix,error=error)
       CALL cp_dbcsr_create(tmp(i)%matrix,template=P(1)%matrix,error=error)
    ENDDO
    DO ispin=1,SIZE(P)/2
       re=2*ispin-1
       im=2*ispin
       imax=MAX(max_iter,1) !if max_iter is 0 then only the deviation from idempotency needs to be calculated
       DO i=1,imax
          CALL cp_complex_dbcsr_gemm_3("N","N",one,P(re)%matrix,P(im)%matrix,&
               P(re)%matrix,P(im)%matrix,zero,PP(re)%matrix,PP(im)%matrix,&
               filter_eps=eps_small,error=error)
          CALL cp_dbcsr_copy(tmp(re)%matrix,PP(re)%matrix,error=error)
          CALL cp_dbcsr_copy(tmp(im)%matrix,PP(im)%matrix,error=error)
          CALL cp_dbcsr_add(tmp(re)%matrix,P(re)%matrix,1.0_dp,-1.0_dp,error=error)
          CALL cp_dbcsr_add(tmp(im)%matrix,P(im)%matrix,1.0_dp,-1.0_dp,error=error)
          CALL complex_frobenius_norm(frob_norm,tmp(re)%matrix,tmp(im)%matrix,error=error)
          IF(unit_nr.gt.0) WRITE(unit_nr,'(t3,a,2f16.8)')"Deviation from idempotency: ",frob_norm
          IF(frob_norm.GT.threshold.AND.max_iter>0)THEN
             CALL cp_dbcsr_copy(tmp(re)%matrix,P(re)%matrix,error=error)
             CALL cp_dbcsr_copy(tmp(im)%matrix,P(im)%matrix,error=error)
             CALL cp_dbcsr_copy(P(re)%matrix,PP(re)%matrix,error=error)
             CALL cp_dbcsr_copy(P(im)%matrix,PP(im)%matrix,error=error)            
             CALL cp_complex_dbcsr_gemm_3("N", "N", -2.0_dp, tmp(re)%matrix,tmp(im)%matrix,&
                  PP(re)%matrix,PP(im)%matrix,3.0_dp,P(re)%matrix,P(im)%matrix,&
                  filter_eps=eps_small,error=error)
             CALL cp_dbcsr_filter(P(re)%matrix,eps,error=error)
             CALL cp_dbcsr_filter(P(im)%matrix,eps,error=error)
             !make sure P is exactly hermitian
             CALL cp_dbcsr_transposed(tmp(re)%matrix,P(re)%matrix,error=error)
             CALL cp_dbcsr_add(P(re)%matrix,tmp(re)%matrix,one/2,one/2,error=error)
             CALL cp_dbcsr_transposed(tmp(im)%matrix,P(im)%matrix,error=error)
             CALL cp_dbcsr_add(P(im)%matrix,tmp(im)%matrix,one/2,-one/2,error=error)
          ELSE
              EXIT
          END IF
       END DO
    END DO
    CALL cp_dbcsr_deallocate_matrix_set(tmp,error)
    CALL cp_dbcsr_deallocate_matrix_set(PP,error)

    CALL timestop(handle)

  END SUBROUTINE purify_mcweeny_complex_orth

! *****************************************************************************
!> \brief calculates the density from the complex MOs and passes the density to
!>        qs_env.
!> \param qs_env ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calc_update_rho(qs_env,error)


    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_update_rho', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, im, ncol, re
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)

    NULLIFY (rho,ks_env,mos,rtp,dft_control)
    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    rho=rho,&
                    rtp=rtp,&
                    dft_control=dft_control,&
                    error=error)
    CALL get_rtp(rtp=rtp,mos_new=mos,error=error)
    DO i=1,SIZE(mos)/2
       re=2*i-1 ; im =2*i
       alpha=3.0_dp-REAL(SIZE(mos)/2,dp)
       CALL cp_dbcsr_set(rho%rho_ao(i)%matrix,0.0_dp,error=error)
       CALL cp_fm_get_info(mos(re)%matrix,ncol_global=ncol,error=error)
       CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho%rho_ao(i)%matrix,&
                               matrix_v=mos(re)%matrix,&
                               ncol=ncol,&
                               alpha=alpha,error=error)
       ! It is actually complex conjugate but i*i=-1 therfore it must be added
       CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho%rho_ao(i)%matrix,&
                        matrix_v=mos(im)%matrix,&
                        ncol=ncol,&
                        alpha=alpha,error=error)
    END DO

    CALL qs_rho_update_rho(rho, qs_env, error=error)

    IF(rtp%do_hfx)THEN
       CALL calculate_P_imaginary(rtp,dft_control%do_admm,error)
    END IF

    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)

    CALL timestop(handle)

  END SUBROUTINE calc_update_rho


! *****************************************************************************
!> \brief Copies the density matrix back into the qs_env%rho%rho_ao
!> \param qs_env ...
!> \param error ...
!> \author Samuel Andermatt (3.14)
! *****************************************************************************

  SUBROUTINE calc_update_rho_sparse(qs_env,error)


    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_update_rho_sparse', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, ispin, re
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_new
    TYPE(cp_dbcsr_type), POINTER             :: S_minus_half, tmp, tmp2
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control

    NULLIFY(rho,ks_env,rtp,tmp,tmp2,dft_control)
    CALL timeset(routineN,handle)
    CALL get_qs_env(qs_env,&
                   ks_env=ks_env,&
                   rho=rho,&
                   rtp=rtp,&
                   dft_control=dft_control,&
                   error=error)
    rtp_control=>dft_control%rtp_control
    CALL get_rtp(rtp=rtp,rho_new=rho_new,error=error)
    IF(rtp_control%orthonormal) THEN
       CALL get_rtp(rtp=rtp,S_minus_half=S_minus_half,error=error)
       ALLOCATE(tmp)
       CALL cp_dbcsr_init(tmp,error=error)
       CALL cp_dbcsr_create(tmp,template=rho_new(1)%matrix,error=error)
       ALLOCATE(tmp2)
       CALL cp_dbcsr_init(tmp2,error=error)
       CALL cp_dbcsr_create(tmp2,template=rho_new(1)%matrix,error=error)
       DO ispin=1,SIZE(rho%rho_ao)
          re = 2*ispin-1
          CALL cp_dbcsr_multiply("N","N",one,S_minus_half,rho_new(re)%matrix,zero,tmp,filter_eps=rtp%filter_eps,error=error)
          CALL cp_dbcsr_multiply("N","N",one,tmp,S_minus_half,zero,tmp2,filter_eps=rtp%filter_eps,error=error)
          CALL cp_dbcsr_set(rho%rho_ao(ispin)%matrix,zero,error=error)
          CALL cp_dbcsr_copy_into_existing(rho%rho_ao(ispin)%matrix,tmp2,error=error)
       END DO
       CALL cp_dbcsr_deallocate_matrix(tmp,error=error)
       CALL cp_dbcsr_deallocate_matrix(tmp2,error=error)
    ELSE
       DO ispin=1,SIZE(rho%rho_ao)
          CALL cp_dbcsr_set(rho%rho_ao(ispin)%matrix,zero,error=error)
          CALL cp_dbcsr_copy_into_existing(rho%rho_ao(ispin)%matrix,rho_new(ispin*2-1)%matrix,error=error)
       END DO
    ENDIF

    CALL qs_rho_update_rho(rho, qs_env, error=error)
    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)


    CALL timestop(handle)

  END SUBROUTINE calc_update_rho_sparse

! *****************************************************************************
!> \brief ...
!> \param rtp ...
!> \param do_admm ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_P_imaginary(rtp,do_admm,error)
    TYPE(rt_prop_type), POINTER              :: rtp
    LOGICAL, INTENT(in)                      :: do_admm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_P_imaginary', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, im, ncol, re
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p_im
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos

    CALL get_rtp(rtp=rtp,mos_new=mos,error=error)
    CALL get_rtp_hfx(rtp_hfx=rtp%rtp_hfx,matrix_p_im=matrix_p_im,error=error)

    DO i=1,SIZE(mos)/2
       re=2*i-1 ; im =2*i
       alpha=3.0_dp-REAL(SIZE(matrix_p_im),dp)
       CALL cp_dbcsr_set(matrix_p_im(i)%matrix,0.0_dp,error=error)
       CALL cp_fm_get_info(mos(re)%matrix,ncol_global=ncol,error=error)
       CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=matrix_p_im(i)%matrix,&
                               matrix_v=mos(im)%matrix,&
                               matrix_g=mos(re)%matrix,&
                               ncol=ncol,&
                               alpha=alpha,error=error)
       ! It is actually complex conjugate not only transposed
       alpha=-alpha
       CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=matrix_p_im(i)%matrix,&
                        matrix_v=mos(re)%matrix,&
                        matrix_g=mos(im)%matrix,&
                        ncol=ncol,&
                        alpha=alpha,error=error)
       
       
    END DO
  
  END SUBROUTINE calculate_P_imaginary

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rtp ...
!> \param aspc_order ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE aspc_extrapolate(qs_env,rtp,aspc_order,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(rt_prop_type), POINTER              :: rtp
    INTEGER, INTENT(in)                      :: aspc_order
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
      POINTER                                :: matrix_s, rho_new, s_hist
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: rho_hist
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: mo_hist
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct, &
                                                matrix_struct_new
    TYPE(cp_fm_type), POINTER                :: fm_tmp, fm_tmp1, fm_tmp2
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(rtp_control_type), POINTER          :: rtp_control

    NULLIFY(rho_hist, dft_control)
    CALL timeset(routineN,handle)
    CALL cite_reference(Kolafa2004)
    CALL get_qs_env(qs_env,&
                    matrix_s=matrix_s,&
                    dft_control=dft_control,&
                    error=error)

    rtp_control=>dft_control%rtp_control
    IF(rtp%linear_scaling) THEN
       CALL get_rtp(rtp=rtp,rho_new=rho_new,error=error)
    ELSE
       CALL get_rtp(rtp=rtp,mos_new=mos_new,error=error)
    ENDIF

    naspc=MIN(rtp%istep,aspc_order)
    IF(rtp%linear_scaling)THEN
       nmat=SIZE(rho_new)
       rho_hist=>rtp%history%rho_history
       DO imat=1,nmat
          DO iaspc=1,naspc
             alpha=(-1.0_dp)**(iaspc + 1)*REAL(iaspc,KIND=dp)*&
                  binomial(2*naspc,naspc - iaspc)/binomial(2*naspc - 2,naspc -1)
             ihist=MOD(rtp%istep-iaspc,aspc_order)+1
             IF(iaspc==1)THEN
                CALL cp_dbcsr_add(rho_new(imat)%matrix,rho_hist(imat,ihist)%matrix,zero,alpha,error)
             ELSE
                CALL cp_dbcsr_add(rho_new(imat)%matrix,rho_hist(imat,ihist)%matrix,one,alpha,error)
             END IF
          END DO
       END DO
    ELSE
       mo_hist=>rtp%history%mo_history
       nmat=SIZE(mos_new)
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

! *****************************************************************************
!> \brief ...
!> \param rtp ...
!> \param mos ...
!> \param rho ...
!> \param s_mat ...
!> \param ihist ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE put_data_to_history(rtp,mos,rho,s_mat,ihist,error)
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: s_mat
    INTEGER                                  :: ihist
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'put_data_to_history', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF(rtp%linear_scaling) THEN
       DO i=1,SIZE(rho)
          CALL cp_dbcsr_copy(rtp%history%rho_history(i,ihist)%matrix,rho(i)%matrix,error=error)
       END DO
    ELSE
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
