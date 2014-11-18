!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for the real time propagation.
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE rt_propagation
  USE cp_control_types,                ONLY: dft_control_type,&
                                             rtp_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_external_control,             ONLY: external_control
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_rm_iter_level
  USE efield_utils,                    ONLY: calculate_ecore_efield
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: real_time_propagation,&
                                             use_restart_wfn,&
                                             use_rt_restart,&
                                             use_scf_wfn
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE md_environment_types,            ONLY: md_environment_type
  USE qs_dftb_matrices,                ONLY: build_dftb_overlap
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_energy_utils,                 ONLY: qs_energies_init
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_external_potential,           ONLY: external_c_potential,&
                                             external_e_potential
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type
  USE rt_delta_pulse,                  ONLY: apply_delta_pulse,&
                                             apply_delta_pulse_periodic
  USE rt_hfx_utils,                    ONLY: rtp_hfx_rebuild
  USE rt_propagation_methods,          ONLY: calc_update_rho,&
                                             calc_update_rho_sparse,&
                                             propagation_step
  USE rt_propagation_output,           ONLY: rt_prop_output
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_create,&
                                             rt_prop_type,&
                                             rtp_history_create
  USE rt_propagation_utils,            ONLY: get_restart_wfn
  USE rt_propagator_init,              ONLY: init_propagators,&
                                             rt_initialize_rho_from_ks
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation'

  PUBLIC :: rt_prop_setup,&
            rt_write_input_restart


CONTAINS

! *****************************************************************************
!> \brief creates rtp_type, gets the initial state, either by reading MO's
!>        from file or calling SCF run
!> \param force_env ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE rt_prop_setup(force_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_prop_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: aspc_order
    INTEGER, DIMENSION(2)                    :: nelectron_spin
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control
    TYPE(section_vals_type), POINTER         :: hfx_sections, input, &
                                                ls_scf_section, md_section, &
                                                motion_section

    failure=.FALSE.
    NULLIFY(qs_env,rtp_control,dft_control)

    CALL force_env_get(force_env=force_env,qs_env=qs_env,globenv=globenv,error=error)
    CALL get_qs_env(qs_env, dft_control=dft_control, error=error)
    rtp_control=>dft_control%rtp_control


    ! Takes care that an initial wavefunction/density is available
    ! Can either be by performing an scf loop or reading a restart
    CALL rt_initial_guess(qs_env,force_env,rtp_control,error)

    ! Initializes the extrapolation
    CALL get_qs_env(qs_env=qs_env,rtp=rtp,input=input,error=error)
    aspc_order=rtp_control%aspc_order
    CALL rtp_history_create(rtp,aspc_order,error=error)


    ! Reads the simulation parameters from the input
    motion_section => section_vals_get_subs_vals(force_env%root_section,"MOTION",error=error)
    md_section     => section_vals_get_subs_vals(motion_section,"MD",error=error)
    hfx_sections => section_vals_get_subs_vals(force_env%root_section,"FORCE_EVAL%DFT%XC%HF",error=error)
    CALL section_vals_val_get(md_section,"TIMESTEP",r_val=qs_env%rtp%dt,error=error)
    CALL section_vals_val_get(md_section,"STEP_START_VAL",i_val=qs_env%rtp%i_start,error=error)
    CALL section_vals_val_get(md_section,"STEPS",i_val=rtp%nsteps,error=error)

    ls_scf_section => section_vals_get_subs_vals(input,"DFT%LS_SCF",error=error)
    CALL section_vals_val_get(ls_scf_section,"EPS_FILTER",r_val=rtp%filter_eps,error=error)
    IF(.NOT.qs_env%rtp%linear_scaling) rtp%filter_eps = 0.0_dp
    IF(rtp_control%acc_ref<1) rtp_control%acc_ref=1
    rtp%filter_eps_small=rtp%filter_eps/rtp_control%acc_ref
    CALL section_vals_val_get(ls_scf_section,"EPS_LANCZOS",r_val=rtp%lanzcos_threshold,error=error)
    CALL section_vals_val_get(ls_scf_section,"MAX_ITER_LANCZOS",i_val=rtp%lanzcos_max_iter,error=error)
    CALL section_vals_val_get(ls_scf_section,"SIGN_SQRT_ORDER",i_val=rtp%newton_schulz_order,error=error)
    CALL section_vals_get(hfx_sections,explicit=rtp%do_hfx,error=error)

    ! Hmm, not really like to initialize with the structure of S but I reckon it is
    ! done everywhere like this
    IF(rtp%do_hfx)&
       CALL rtp_hfx_rebuild(qs_env,error)

    IF(qs_env%rtp%linear_scaling.AND.rtp_control%initial_wfn==use_scf_wfn) THEN
       CALL get_qs_env(qs_env,&
                       matrix_ks=matrix_ks,&
                       matrix_s=matrix_s,&
                       nelectron_spin=nelectron_spin,&
                       error=error)
       CALL rt_initialize_rho_from_ks(rtp,matrix_ks,matrix_s,nelectron_spin,rtp_control%orthonormal,error=error)
       CALL calc_update_rho_sparse(qs_env,error)
    ENDIF

    IF(rtp_control%fixed_ions) THEN
       CALL init_propagation_run(qs_env,error)
       CALL run_propagation(qs_env,force_env,globenv,error)
    ELSE
       CALL init_ehrenfest_md(force_env,qs_env,error)
    END IF


  END SUBROUTINE rt_prop_setup

! *****************************************************************************
!> \brief calculates the matrices needed in the first step of RTP
!> \param qs_env ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE init_propagation_run(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: i
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, mos_old
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control

    NULLIFY(matrix_s,dft_control)

    CALL get_qs_env(qs_env,&
                    rtp=rtp,&
                    matrix_s=matrix_s,&
                    dft_control=dft_control,&
                    error=error)
    rtp_control=>dft_control%rtp_control

    IF(.NOT.rtp%linear_scaling) THEN
       CALL get_qs_env(qs_env,mos=mos,error=error)
       CALL get_rtp(rtp=rtp,mos_old=mos_old,mos_new=mos_new,error=error)
       IF(rtp_control%initial_wfn==use_scf_wfn)THEN
         IF (rtp_control%apply_delta_pulse) THEN
            IF(dft_control%qs_control%dftb)&
               CALL build_dftb_overlap(qs_env,1,matrix_s,error)
            IF (rtp_control%periodic) THEN
              CALL apply_delta_pulse_periodic(qs_env,mos_old,mos_new,error)
            ELSE
              CALL apply_delta_pulse(qs_env,mos_old,mos_new,error)
            ENDIF
         ELSE
            DO i=1,SIZE(mos)
                CALL cp_fm_to_fm(mos(i)%mo_set%mo_coeff,mos_old(2*i-1)%matrix,error)
                CALL cp_fm_set_all(mos_old(2*i)%matrix,zero,zero,error)
            END DO
         ENDIF
       END IF

       DO i=1,SIZE(mos_old)
          CALL cp_fm_to_fm(mos_old(i)%matrix,mos_new(i)%matrix,error)
       END DO
       CALL calc_update_rho(qs_env,error)
    END IF

    CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., error=error)

    CALL init_propagators(qs_env,error)

  END SUBROUTINE init_propagation_run

! *****************************************************************************
!> \brief calculates the matrices needed in the first step of EMD
!> \param force_env ...
!> \param qs_env ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE init_ehrenfest_md(force_env,qs_env,error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: i
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, mos_old
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control

    NULLIFY(dft_control)

     CALL get_qs_env(qs_env,&
                     mos=mos,&
                     rtp=rtp,&
                     energy=energy,&
                     dft_control=dft_control,&
                     error=error)
     rtp_control=>dft_control%rtp_control

     IF(.NOT.rtp%linear_scaling) THEN 
        CALL get_rtp(rtp=rtp,mos_old=mos_old,mos_new=mos_new,error=error)

        IF(rtp_control%initial_wfn==use_scf_wfn)THEN
           DO i=1,SIZE(mos)
               CALL cp_fm_to_fm(mos(i)%mo_set%mo_coeff,mos_old(2*i-1)%matrix,error)
               CALL cp_fm_set_all(mos_old(2*i)%matrix,zero,zero,error)
           END DO
        END IF

        DO i=1,SIZE(mos_old)
           CALL cp_fm_to_fm(mos_old(i)%matrix,mos_new(i)%matrix,error)
        END DO
     ENDIF

     rtp_control%initial_step=.TRUE.
     CALL force_env_calc_energy_force(force_env,calc_force=.TRUE.,&
          error=error)

     rtp_control%initial_step=.FALSE.
     rtp%energy_old=energy%total


  END SUBROUTINE init_ehrenfest_md

! *****************************************************************************
!> \brief performes the real RTP run, gets information from MD section
!>        uses MD as iteration level
!> \param qs_env ...
!> \param force_env ...
!> \param globenv ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE run_propagation(qs_env,force_env,globenv,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'run_propagation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: aspc_order, handle, i_iter, &
                                                i_step, max_iter, max_steps
    LOGICAL                                  :: failure, should_stop
    REAL(Kind=dp)                            :: eps_ener, time_iter_start, &
                                                time_iter_stop, used_time
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control

    failure=.FALSE.
    should_stop=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(logger,dft_control,energy,rtp,rtp_control)
    logger   => cp_error_get_logger(error)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,rtp=rtp,energy=energy,error=error)

    rtp_control=>dft_control%rtp_control
    max_steps=rtp%nsteps
    max_iter=rtp_control%max_iter
    eps_ener=rtp_control%eps_ener

    aspc_order=rtp_control%aspc_order

    rtp%energy_old=energy%total
    time_iter_start=m_walltime()
    CALL cp_add_iter_level(logger%iter_info,"MD",error=error)
    CALL cp_iterate(logger%iter_info,iter_nr=0,error=error)
    DO i_step=rtp%i_start+1,max_steps
       energy%efield_core=0.0_dp
       qs_env%sim_time=REAL(i_step,dp)*rtp%dt
       qs_env%sim_step=i_step
       rtp%istep=i_step-rtp%i_start
       CALL calculate_ecore_efield(qs_env,.FALSE.,error=error)
       CALL external_c_potential(qs_env,calculate_forces=.FALSE.,error=error)
       CALL external_e_potential(qs_env,error=error)
       CALL cp_iterate(logger%iter_info,last=(i_step==max_steps),iter_nr=i_step,error=error)
       rtp%converged=.FALSE.
       DO i_iter=1,max_iter
          IF(i_step==rtp%i_start+1.AND.i_iter==2.AND.rtp_control%hfx_redistribute)&
              CALL qs_ks_did_change(qs_env%ks_env,s_mstruct_changed=.TRUE., error=error)
          rtp%iter=i_iter
          CALL propagation_step(qs_env,rtp, rtp_control, error=error)
          CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., error=error)
          rtp%energy_new=energy%total
          IF(rtp%converged)EXIT
          CALL rt_prop_output(qs_env,real_time_propagation,rtp%delta_iter,error=error)
       END DO
       IF(rtp%converged)THEN
          CALL external_control(should_stop,"MD",globenv=globenv,error=error)
          IF (should_stop)CALL cp_iterate(logger%iter_info,last=.TRUE.,iter_nr=i_step,error=error)
          time_iter_stop=m_walltime()
          used_time= time_iter_stop - time_iter_start
          time_iter_start=time_iter_stop
          CALL rt_prop_output(qs_env,real_time_propagation,delta_iter=rtp%delta_iter,used_time=used_time,error=error)
          CALL rt_write_input_restart(force_env=force_env,error=error)
          IF (should_stop)  EXIT
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
!> \brief overwrites some values in the input file such that the .restart
!>        file will contain the appropriate information
!> \param md_env ...
!> \param force_env ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE rt_write_input_restart(md_env,force_env,error)
    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: md_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_write_input_restart', &
      routineP = moduleN//':'//routineN

    TYPE(section_vals_type), POINTER         :: motion_section, root_section, &
                                                rt_section

    root_section => force_env%root_section
    motion_section => section_vals_get_subs_vals(root_section,"MOTION",error=error)
    rt_section => section_vals_get_subs_vals(root_section,"FORCE_EVAL%DFT%REAL_TIME_PROPAGATION",error=error)
    CALL section_vals_val_set(rt_section,"INITIAL_WFN",i_val=use_rt_restart,error=error)
    ! coming from RTP
    IF (.NOT. PRESENT(md_env)) THEN
     CALL section_vals_val_set(motion_section,"MD%STEP_START_VAL",i_val=force_env%qs_env%sim_step,error=error)
    ENDIF

    CALL write_restart(md_env=md_env,root_section=root_section,error=error)

  END SUBROUTINE rt_write_input_restart

! *****************************************************************************
!> \brief Creates the initial electronic states and allocates the necessary 
!>        matrices 
!> \param qs_env ...
!> \param force_env ...
!> \param rtp_control ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE rt_initial_guess(qs_env,force_env,rtp_control,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(rtp_control_type), POINTER          :: rtp_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_initial_guess', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: homo, ispin
    LOGICAL                                  :: energy_consistency, failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools_aux_fit
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(dft_control_type), POINTER          :: dft_control

    NULLIFY(matrix_s,dft_control)
    CALL get_qs_env(qs_env, dft_control=dft_control, error=error)

    SELECT CASE(rtp_control%initial_wfn)
    CASE(use_scf_wfn)
       qs_env%sim_time=0.0_dp
       qs_env%sim_step=0
       energy_consistency=.TRUE.
       !in the linear scaling case we need a correct kohn-sham matrix, which we cannot get with consistent energies
       IF(rtp_control%linear_scaling) energy_consistency=.FALSE.
       CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,&
            consistent_energies=energy_consistency,error=error)
       qs_env%run_rtp=.TRUE.
       ALLOCATE(qs_env%rtp)
       CALL get_qs_env(qs_env, matrix_s=matrix_s, error=error)
       CALL rt_prop_create(qs_env%rtp,qs_env%mos,qs_env%mpools,dft_control,matrix_s(1)%matrix,&
                           rtp_control%linear_scaling,rtp_control%write_restart,qs_env%mos_aux_fit,error)

    CASE(use_restart_wfn,use_rt_restart)
       CALL qs_energies_init(qs_env, .FALSE. , error)
       IF(.NOT.rtp_control%linear_scaling.OR.rtp_control%initial_wfn==use_restart_wfn) THEN
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
          IF(dft_control%do_admm) THEN
             CALL mpools_get(qs_env%mpools_aux_fit, ao_mo_fm_pools=ao_mo_fm_pools_aux_fit,&
                  error=error)
             CPPrecondition(ASSOCIATED(qs_env%mos_aux_fit),cp_failure_level,routineP,error,failure)
             DO ispin=1,SIZE(qs_env%mos_aux_fit)
                CALL get_mo_set(qs_env%mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit,homo=homo)
                IF (.NOT.ASSOCIATED(mo_coeff_aux_fit)) THEN
                   CALL init_mo_set(qs_env%mos_aux_fit(ispin)%mo_set,&
                        ao_mo_fm_pools_aux_fit(ispin)%pool,&
                        name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                        "%mo_aux_fit"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                        error=error)
                END IF
             END DO
          END IF
       ENDIF
       ALLOCATE(qs_env%rtp)
       CALL get_qs_env(qs_env, matrix_s=matrix_s, error=error)
       CALL rt_prop_create(qs_env%rtp,qs_env%mos,qs_env%mpools,dft_control,matrix_s(1)%matrix,&
                           rtp_control%linear_scaling,rtp_control%write_restart,qs_env%mos_aux_fit,error)
       CALL get_restart_wfn(qs_env,error)

       qs_env%run_rtp=.TRUE.
    END SELECT
 
  END SUBROUTINE rt_initial_guess

END MODULE rt_propagation
