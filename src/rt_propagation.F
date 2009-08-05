!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for the real time propagation.
!> \author Florian Schiffmann (02.09)
! *****************************************************************************
MODULE rt_propagation

  USE cp_control_types,                ONLY: rtp_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale,&
                                             cp_fm_upper_to_full
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_rm_iter_level
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm
  USE efield_utils,                    ONLY: calculate_ecore_efield
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       do_arnoldi, do_cn, do_diag, do_em, do_etrs, do_pade, do_taylor, &
       extr_s_ks, real_time_propagation, use_restart_wfn, use_rt_restart, &
       use_scf_wfn
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type
  USE rt_make_propagators,             ONLY: make_etrs_exp
  USE rt_matrix_exp,                   ONLY: get_nsquare_norder
  USE rt_propagation_methods,          ONLY: calc_SinvH,&
                                             calc_update_rho,&
                                             propagation_step,&
                                             put_data_to_history,&
                                             s_matrices_create,&
                                             update_core_and_matrices
  USE rt_propagation_output,           ONLY: rt_prop_output
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_create,&
                                             rt_prop_type,&
                                             rtp_history_create
  USE rt_propagation_utils,            ONLY: calc_S_derivs,&
                                             get_restart_wfn
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE termination,                     ONLY: external_control
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation'

  PUBLIC :: rt_prop_setup,&
       init_ehrenfest_md,&
       run_propagation


CONTAINS

! *****************************************************************************
!> \brief creates rtp_type, gets the initial state, either by reading MO's
!>        from file or calling SCF run
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE rt_prop_setup(force_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_prop_setup', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: aspc_order, homo, i, ispin, &
                                                stat
    LOGICAL                                  :: extrapolate_H, failure
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, mos_old
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control
    TYPE(section_vals_type), POINTER         :: md_section, motion_section

    failure=.FALSE.
    NULLIFY(qs_env,rtp_control)
    
    CALL force_env_get(force_env=force_env,qs_env=qs_env,globenv=globenv,error=error)
    
    rtp_control=>qs_env%dft_control%rtp_control
     


    SELECT CASE(rtp_control%initial_wfn)
    CASE(use_scf_wfn)
       rtp_control%calc_restart=.FALSE.
       IF(rtp_control%fixed_ions)THEN
          qs_env%sim_time=0.0_dp
          qs_env%sim_step=0
          CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,&
               consistent_energies=.TRUE.,error=error)
       ELSE
          CALL force_env_calc_energy_force(force_env,calc_force=.TRUE.,&
               consistent_energies=.TRUE.,error=error) 
       END IF
       qs_env%run_rtp=.TRUE.
       ALLOCATE(qs_env%rtp, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
       CALL rt_prop_create(qs_env%rtp,qs_env%mos,qs_env%mpools,qs_env%dft_control,error)
       aspc_order=rtp_control%aspc_order
       extrapolate_H=(rtp_control%extrapolation==extr_s_ks)
       CALL rtp_history_create(qs_env%rtp,extrapolate_H,aspc_order,error)

    CASE(use_restart_wfn,use_rt_restart)
       CALL update_core_and_matrices(qs_env,globenv,error)       
       DO ispin=1,SIZE(qs_env%mos)
          CALL get_mo_set(qs_env%mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo)
          IF (.NOT.ASSOCIATED(mo_coeff)) THEN
             CALL init_mo_set(qs_env%mos(ispin)%mo_set,&
                  qs_env%mpools%ao_mo_fm_pools(ispin)%pool,&
                  name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                  "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                  error=error)
          END IF
       END DO
       ALLOCATE(qs_env%rtp, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       
       CALL rt_prop_create(qs_env%rtp,qs_env%mos,qs_env%mpools,qs_env%dft_control,error)
       aspc_order=rtp_control%aspc_order
       extrapolate_H=(rtp_control%extrapolation==extr_s_ks)
       CALL rtp_history_create(qs_env%rtp,extrapolate_H,aspc_order,error)
       
       CALL get_restart_wfn(qs_env,error)
       rtp_control%calc_restart=.TRUE.
       qs_env%run_rtp=.TRUE.
    END SELECT


    motion_section => section_vals_get_subs_vals(force_env%root_section,"MOTION",error=error)
    md_section     => section_vals_get_subs_vals(motion_section,"MD",error=error)
    
    CALL section_vals_val_get(md_section,"TIMESTEP",r_val=qs_env%rtp%dt,error=error)      
    CALL section_vals_val_get(md_section,"STEP_START_VAL",i_val=qs_env%rtp%i_start,error=error)

    IF(rtp_control%fixed_ions)THEN
       
       IF(rtp_control%initial_wfn==use_restart_wfn.OR.rtp_control%initial_wfn==use_rt_restart)THEN
          CALL calc_update_rho(qs_env,error)
          CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,calculate_forces=.FALSE.,&
               error=error)
       END IF
       CALL section_vals_val_get(md_section,"STEPS",i_val=qs_env%rtp%nsteps,error=error)
       CALL init_propagation_run(qs_env,error)
       CALL run_propagation(qs_env,globenv,error)
    ELSE
       IF(rtp_control%initial_wfn==use_restart_wfn.OR.rtp_control%initial_wfn==use_rt_restart)THEN
          CALL calc_update_rho(qs_env,error)
          CALL force_env_calc_energy_force(force_env,calc_force=.TRUE.,&
               error=error) 
          CALL get_qs_env(qs_env,&
               mos=mos,&
               rtp=rtp,&
               error=error)
          CALL get_rtp (rtp=rtp,mos_old=mos_old,mos_new=mos_new,error=error)
          DO i=1,SIZE(mos_new)
             CALL cp_fm_to_fm(mos_new(i)%matrix,mos_old(i)%matrix,error)
          END DO
       END IF
    END IF


  END SUBROUTINE rt_prop_setup

! *****************************************************************************
!> \brief calculates the matrices needed in the first step of RTP
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE init_propagation_run(qs_env,error)   
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: i
    REAL(kind=dp)                            :: t
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, mos_old
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: s_mat
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL get_qs_env(qs_env,&
         mos=mos,&
         rtp=rtp,&
         matrix_s=s_mat,&
         error=error)

    CALL get_rtp (rtp=rtp,mos_old=mos_old,mos_new=mos_new,dt=t,error=error)
    IF(qs_env%dft_control%rtp_control%initial_wfn.NE.use_rt_restart)THEN
       DO i=1,SIZE(qs_env%mos)
          CALL cp_fm_to_fm(mos(i)%mo_set%mo_coeff,mos_old(2*i-1)%matrix,error)
          CALL cp_fm_set_all(mos_old(2*i)%matrix,zero,zero,error)
       END DO
    ELSE
       DO i=1,SIZE(mos_old)
          CALL cp_fm_to_fm(mos_new(i)%matrix,mos_old(i)%matrix,error)
       END DO
    END IF
    CALL init_propagators(qs_env,error)   

  END SUBROUTINE init_propagation_run

! *****************************************************************************
!> \brief calculates the matrices needed in the first step of EMD
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE init_ehrenfest_md(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: i
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H, mos_old
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: s_mat
    TYPE(rt_prop_type), POINTER              :: rtp

     CALL get_qs_env(qs_env,&
         mos=mos,&
         rtp=rtp,&
         energy=energy,&
         matrix_s=s_mat,&
         error=error)
     qs_env%run_rtp=.TRUE.
     qs_env%dft_control%rtp_control%calc_restart=.FALSE.
     rtp%energy_old=energy%total
     
     CALL get_rtp (rtp=rtp,exp_H_old=exp_H,mos_old=mos_old,error=error)   
     IF(qs_env%dft_control%rtp_control%initial_wfn.NE.use_rt_restart)THEN
        DO i=1,SIZE(mos)
           CALL cp_fm_to_fm(mos(i)%mo_set%mo_coeff,mos_old(2*(i-1)+1)%matrix,error)
           CALL cp_fm_set_all(mos_old(2*(i-1)+2)%matrix,zero,zero,error)
        END DO
     END IF
     CALL calc_S_derivs(qs_env,error)
     CALL init_propagators(qs_env,error)

    
  END SUBROUTINE init_ehrenfest_md

! *****************************************************************************
!> \brief performes the real RTP run, gets information from MD section
!>        uses MD as iteration level 
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE run_propagation(qs_env,globenv,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'run_propagation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: aspc_order, handle, i, &
                                                i_iter, i_step, ihist, &
                                                max_iter, max_steps
    LOGICAL                                  :: failure, should_stop
    REAL(Kind=dp)                            :: delta_iter, eps_ener, &
                                                time_iter_start, &
                                                time_iter_stop, used_time
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, &
                                                mos_new, mos_old
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(rt_prop_type), POINTER              :: rtp

    failure=.FALSE.
    should_stop=.FALSE.
    CALL timeset(routineN,handle)
    logger   => cp_error_get_logger(error) 
    max_steps=qs_env%rtp%nsteps
    max_iter=qs_env%dft_control%rtp_control%max_iter
    eps_ener=qs_env%dft_control%rtp_control%eps_ener
   
 
    CALL get_qs_env(qs_env=qs_env,rtp=rtp,energy=energy,error=error)
    aspc_order=qs_env%dft_control%rtp_control%aspc_order

    CALL get_rtp(rtp=rtp,&
                 exp_H_old=exp_H_old,&
                 exp_H_new=exp_H_new,&
                 mos_old=mos_old,&
                 mos_new=mos_new,&
                 error=error)
                 
    rtp%energy_old=energy%total
    time_iter_start=m_walltime()
    CALL cp_add_iter_level(logger%iter_info,"MD",error=error)
    CALL cp_iterate(logger%iter_info,iter_nr=0,error=error)
    DO i_step=rtp%i_start+1,max_steps
       qs_env%energy%efield_core=0.0_dp

       qs_env%sim_time=REAL(i_step,dp)*rtp%dt
       qs_env%sim_step=i_step
       rtp%istep=i_step
       CALL calculate_ecore_efield(qs_env,.FALSE.,error)
       CALL cp_iterate(logger%iter_info,last=(i_step==max_steps),iter_nr=i_step,error=error)
       rtp%converged=.FALSE.
       DO i_iter=1,max_iter
          rtp%iter=i_iter
          CALL propagation_step(qs_env,delta_iter, error)
          rtp%energy_new=energy%total
          rtp%converged=(delta_iter.LT.eps_ener)

          IF(rtp%converged)EXIT
          CALL rt_prop_output(qs_env,real_time_propagation,delta_iter,error=error)
          
       END DO

       IF(rtp%converged)THEN
          ihist=MOD(rtp%istep,aspc_order)+1
          DO i=1,SIZE(mos_new)
             CALL cp_fm_to_fm(mos_new(i)%matrix,mos_old(i)%matrix,error)
             CALL cp_fm_to_fm(exp_H_new(i)%matrix,exp_H_old(i)%matrix,error)
          END DO
          IF(rtp%history%mos_or_H==2)&
             CALL  put_data_to_history(qs_env,exp_H=exp_H_new,ihist=ihist,error=error)
          IF(rtp%history%mos_or_H==1)&
             CALL put_data_to_history(qs_env,mos=mos_new,ihist=ihist,error=error)
          CALL external_control(should_stop,"MD",globenv,error=error)
          IF (should_stop)CALL cp_iterate(logger%iter_info,last=.TRUE.,iter_nr=i_step,error=error)
          time_iter_stop=m_walltime()
          used_time= time_iter_stop - time_iter_start
          time_iter_start=time_iter_stop
          CALL rt_prop_output(qs_env,real_time_propagation,delta_iter,used_time,error=error)
          IF (should_stop)EXIT
       ELSE
          EXIT
       END IF


    END DO
    CALL cp_rm_iter_level(logger%iter_info,"MD",error=error)
    IF(.NOT.rtp%converged)&
         CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
         routineP,"propagation did not converge, either increase MAX_ITER or use a smaller TIMESTEP",&
         error,failure)
    CALL timestop(handle)

  END SUBROUTINE run_propagation

! *****************************************************************************
!> \brief prepares the initial matrices for the propagators
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE init_propagators(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_propagators', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: diagonalize
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, mos_old
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: s_mat
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL get_qs_env(qs_env,&
         mos=mos,&
         rtp=rtp,&
         matrix_s=s_mat,&
         error=error)

    CALL get_rtp (rtp=rtp,exp_H_old=exp_H_old,exp_H_new=exp_H_new,mos_old=mos_old,error=error)


    
    SELECT CASE(qs_env%dft_control%rtp_control%mat_exp)
    CASE(do_taylor,do_pade,do_diag)
       CALL get_maxabs_eigval(rtp,qs_env,error)
       SELECT CASE(qs_env%dft_control%rtp_control%propagator)
       CASE (do_cn,do_em)
          diagonalize=(qs_env%dft_control%rtp_control%mat_exp==do_diag)
          IF(qs_env%dft_control%rtp_control%propagator==do_cn)diagonalize=.FALSE.
          CALL s_matrices_create (s_mat,rtp,diagonalize,error)
          CALL calc_SinvH(exp_H_old,qs_env,error)
          DO i=1,SIZE(exp_H_old)
             CALL cp_fm_to_fm(exp_H_old(i)%matrix,rtp%exp_H_new(i)%matrix,error)
          END DO
          
       CASE (do_etrs)
          
          diagonalize=(qs_env%dft_control%rtp_control%mat_exp==do_diag)
          CALL s_matrices_create (s_mat,rtp,diagonalize,error)
          
          CALL make_etrs_exp(exp_H_old,qs_env,error)
          DO i=1,SIZE(exp_H_old)
             CALL cp_fm_to_fm(exp_H_old(i)%matrix,rtp%exp_H_new(i)%matrix,error)
          END DO
       END SELECT
    CASE(do_arnoldi)
       diagonalize=.FALSE.
       CALL s_matrices_create (s_mat,rtp,diagonalize,error)
       CALL calc_SinvH(exp_H_old,qs_env,error)
       DO i=1,SIZE(exp_H_old)
          CALL cp_fm_to_fm(exp_H_old(i)%matrix,rtp%exp_H_new(i)%matrix,error)
       END DO
    END SELECT

    IF(rtp%history%mos_or_H==2)&
         CALL put_data_to_history(qs_env,exp_H=exp_H_old,ihist=1,error=error)
    IF(rtp%history%mos_or_H==1)&
         CALL put_data_to_history(qs_env,mos=mos_old,s_mat=s_mat,ihist=1,error=error)


  END SUBROUTINE init_propagators

  SUBROUTINE get_maxabs_eigval(rtp,qs_env,error)
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_maxabs_eigval', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: ispin, method, ndim, stat
    LOGICAL                                  :: emd, failure, is_present
    REAL(dp)                                 :: max_eval, min_eval, norm2, &
                                                scale, t
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eigval_H
    TYPE(cp_fm_type), POINTER                :: H_fm, tmp_mat_H
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, s_mat

    failure=.FALSE.

    is_present=(ASSOCIATED(rtp%S_half))
    IF(.NOT.is_present)THEN
       CALL cp_fm_create(rtp%S_half,&
            matrix_struct=rtp%S_inv%matrix_struct,&
            name="S_half",&
            error=error)
       CALL cp_fm_create(rtp%S_minus_half,&
            matrix_struct=rtp%S_inv%matrix_struct,&
            name="S_minus_half",&
            error=error)
    END IF
    
    t=rtp%dt
    CALL get_qs_env(qs_env=qs_env,matrix_ks=matrix_ks,matrix_s=s_mat,error=error) 
    CALL  s_matrices_create (s_mat,rtp,.TRUE.,error)

    CALL cp_fm_create(H_fm,&
         matrix_struct=rtp%S_inv%matrix_struct,&
         name="RTP_H_FM",&
         error=error)
    
    CALL cp_fm_create(tmp_mat_H,&
         matrix_struct=rtp%S_inv%matrix_struct,&
         name="TMP_H",&
         error=error)
    
    ndim=H_fm%matrix_struct%nrow_global   
    scale=1.0_dp
    IF(qs_env%dft_control%rtp_control%propagator==do_etrs)scale=2.0_dp
    t=-t/scale
    ALLOCATE(eigval_H(ndim),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF(qs_env%dft_control%rtp_control%mat_exp==do_taylor)method=1
    IF(qs_env%dft_control%rtp_control%mat_exp==do_pade)method=2
    emd=(.NOT.qs_env%dft_control%rtp_control%fixed_ions)


    DO ispin=1,SIZE(matrix_ks)
       
       CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,H_fm,error=error)
       CALL cp_fm_upper_to_full(H_fm,tmp_mat_H,error)
       CALL cp_fm_scale(t,H_fm,error)
       ndim=H_fm%matrix_struct%nrow_global

       
       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,H_fm,rtp%S_minus_half,zero,&
            tmp_mat_H,error)
       CALL cp_fm_gemm("N","N",ndim,ndim,ndim,one,rtp%S_minus_half,tmp_mat_H,zero,&
            H_fm,error)
       
       CALL cp_fm_syevd(H_fm,tmp_mat_H,eigval_H,error)
       min_eval=MINVAL(eigval_H)
       max_eval=MAXVAL(eigval_H)
       norm2=2.0_dp*MAX(ABS(min_eval),ABS(max_eval))
       CALL get_nsquare_norder(norm2,rtp%orders(1,ispin),rtp%orders(2,ispin),&
            qs_env%dft_control%rtp_control%eps_exp,method ,emd,error)
    END DO
    
    DEALLOCATE(eigval_H,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(.NOT.is_present)THEN
       CALL cp_fm_release(rtp%S_half,error)
       CALL cp_fm_release(rtp%S_minus_half,error)
    END IF

    CALL cp_fm_release(H_fm,error)
    CALL cp_fm_release(tmp_mat_H,error)

  END SUBROUTINE get_maxabs_eigval


END MODULE rt_propagation



      

 


