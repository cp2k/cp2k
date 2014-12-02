!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for the Quickstep SCF run.
!> \par History
!>      - Joost VandeVondele (02.2002)
!>           added code for: incremental (pab and gvg) update
!>                            initialisation (init_cube, l_info)
!>      - Joost VandeVondele (02.2002)
!>           called the poisson code of the classical part
!>           this takes into account the spherical cutoff and allows for
!>           isolated systems
!>      - Joost VandeVondele (02.2002)
!>           added multiple grid feature
!>           changed to spherical cutoff consistently (?)
!>           therefore removed the gradient correct functionals
!>      - updated with the new QS data structures (10.04.02,MK)
!>      - copy_matrix replaced by transfer_matrix (11.04.02,MK)
!>      - nrebuild_rho and nrebuild_gvg unified (12.04.02,MK)
!>      - set_mo_occupation for smearing of the MO occupation numbers
!>        (17.04.02,MK)
!>      - MO level shifting added (22.04.02,MK)
!>      - Usage of TYPE mo_set_p_type
!>      - Joost VandeVondele (05.2002)
!>            added cholesky based diagonalisation
!>      - 05.2002 added pao method [fawzi]
!>      - parallel FFT (JGH 22.05.2002)
!>      - 06.2002 moved KS matrix construction to qs_build_KS_matrix.F [fawzi]
!>      - started to include more LSD (01.2003,Joost VandeVondele)
!>      - 02.2003 scf_env [fawzi]
!>      - got rid of nrebuild (01.2004, Joost VandeVondele)
!>      - 10.2004 removed pao [fawzi]
!>      - 03.2006 large cleaning action [Joost VandeVondele]
!>      - High-spin ROKS added (05.04.06,MK)
!>      - Mandes (10.2013)
!>        intermediate energy communication with external communicator added
!>      - kpoints (08.2014, JGH)
!>      - unified k-point and gamma-point code (2014.11) [Ole Schuett]
!> \author Matthias Krack (30.04.2001)
! *****************************************************************************
MODULE qs_scf
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_name,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_fm_types,                     ONLY: cp_fm_release,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_p_file,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_methods,               ONLY: get_results,&
                                             test_for_result
  USE cp_result_types,                 ONLY: cp_result_type
  USE harris_env_types,                ONLY: harris_env_type
  USE harris_functional,               ONLY: harris_postprocessing
  USE input_constants,                 ONLY: do_method_mndo,&
                                             history_guess,&
                                             ot_precond_full_all,&
                                             ot_precond_full_single,&
                                             ot_precond_full_single_inverse,&
                                             ot_precond_none,&
                                             ot_precond_s_inverse
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE kpoint_types,                    ONLY: kpoint_type
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE message_passing,                 ONLY: mp_send
  USE particle_types,                  ONLY: particle_type
  USE preconditioner,                  ONLY: prepare_preconditioner,&
                                             restart_preconditioner
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_type
  USE qs_block_davidson_types,         ONLY: block_davidson_deallocate
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_density_mixing_types,         ONLY: gspace_mixing_nr
  USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                             qs_diis_b_create
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_io,                        ONLY: write_mo_set
  USE qs_mo_methods,                   ONLY: make_basis_simple,&
                                             make_basis_sm
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             set_mo_occupation
  USE qs_ot,                           ONLY: qs_ot_new_preconditioner
  USE qs_ot_scf,                       ONLY: ot_scf_init,&
                                             ot_scf_read_input
  USE qs_outer_scf,                    ONLY: outer_loop_optimize,&
                                             outer_loop_update_qs_env
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_scf_initialization,           ONLY: qs_scf_compute_properties,&
                                             qs_scf_env_initialize
  USE qs_scf_loop_utils,               ONLY: &
       qs_scf_check_inner_exit, qs_scf_check_outer_exit, &
       qs_scf_density_mixing, qs_scf_harris_e_correct, qs_scf_inner_finalize, &
       qs_scf_new_mos, qs_scf_new_mos_kp, qs_scf_rho_update, &
       qs_scf_set_loop_flags
  USE qs_scf_output,                   ONLY: qs_scf_loop_info,&
                                             qs_scf_loop_print,&
                                             qs_scf_outer_loop_info,&
                                             qs_scf_write_mos
  USE qs_scf_types,                    ONLY: block_davidson_diag_method_nr,&
                                             block_krylov_diag_method_nr,&
                                             general_diag_method_nr,&
                                             ot_diag_method_nr,&
                                             ot_method_nr,&
                                             qs_scf_env_type,&
                                             scf_env_release,&
                                             special_diag_method_nr
  USE qs_wf_history_methods,           ONLY: wfi_update
  USE scf_control_types,               ONLY: scf_control_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf'

  PUBLIC :: scf, scf_env_cleanup, scf_env_do_scf

CONTAINS

! *****************************************************************************
!> \brief perform an scf procedure in the given qs_env
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 introduced scf_env, moved real work to scf_env_do_scf [fawzi]
!> \author fawzi
!> \note
! *****************************************************************************
  SUBROUTINE scf(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: max_scf_tmp
    LOGICAL                                  :: converged, failure, &
                                                outer_scf_loop, should_stop
    LOGICAL, SAVE                            :: first_step_flag = .TRUE.
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                scf_section

    NULLIFY(scf_env)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
      CALL get_qs_env(qs_env,scf_env=scf_env,error=error,input=input, &
                       dft_control=dft_control,scf_control=scf_control)
      IF(scf_control%max_scf > 0) THEN

        dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
        scf_section => section_vals_get_subs_vals(dft_section,"SCF",error=error)

        IF(.NOT. ASSOCIATED(scf_env)) THEN
          CALL qs_scf_env_initialize(qs_env,scf_env,error=error)
!       Moved here from qs_scf_env_initialize to be able to have more scf_env
          CALL set_qs_env(qs_env,scf_env=scf_env,error=error)
          CALL scf_env_release(scf_env,error=error)
          CALL get_qs_env(qs_env=qs_env,scf_env=scf_env,error=error)
        ELSE
          CALL qs_scf_env_initialize(qs_env,scf_env,error=error)
        ENDIF

        IF ((scf_control%density_guess .EQ. history_guess) .AND. (first_step_flag)) THEN
          max_scf_tmp = scf_control%max_scf
          scf_control%max_scf = 1
          outer_scf_loop = scf_control%outer_scf%have_scf
          scf_control%outer_scf%have_scf = .FALSE.
        END IF

        CALL scf_env_do_scf(scf_env=scf_env, scf_control=scf_control, qs_env=qs_env,  &
                          converged=converged, should_stop=should_stop, error=error)

        !   *** add the converged wavefunction to the wavefunction history
        IF ((ASSOCIATED(qs_env%wf_history)) .AND. &
            ((scf_control%density_guess .NE. history_guess) .OR. &
             (.NOT. first_step_flag))) THEN
            CALL wfi_update(qs_env%wf_history,qs_env=qs_env,dt=1.0_dp, error=error)
        ELSE IF ((scf_control%density_guess .EQ. history_guess) .AND. &
                 (first_step_flag)) THEN
          scf_control%max_scf = max_scf_tmp
          scf_control%outer_scf%have_scf = outer_scf_loop
          first_step_flag = .FALSE.
        END IF

        ! *** compute properties that depend on the converged wavefunction
        IF (.NOT.(should_stop)) CALL qs_scf_compute_properties(qs_env,dft_section,error)

        ! *** cleanup
        CALL scf_env_cleanup(scf_env,error=error)

      END IF
     END IF

  END SUBROUTINE scf

! *****************************************************************************
!> \brief perform an scf loop
!> \param scf_env the scf_env where to perform the scf procedure
!> \param scf_control ...
!> \param qs_env the qs_env, the scf_env lives in
!> \param converged will be true / false if converged is reached
!> \param should_stop ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      long history, see cvs and qs_scf module history
!>      02.2003 introduced scf_env [fawzi]
!>      09.2005 Frozen density approximation [TdK]
!>      06.2007 Check for SCF iteration count early [jgh]
!> \author Matthias Krack
!> \note
! *****************************************************************************
  SUBROUTINE scf_env_do_scf(scf_env,scf_control,qs_env,converged,should_stop,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(OUT)                     :: converged, should_stop
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scf_env_do_scf', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: description
    INTEGER :: ext_master_id, external_comm, handle, handle2, i_tmp, ic, &
      ispin, iter_count, output_unit, scf_energy_message_tag, total_steps
    LOGICAL :: diis_step, do_kpoints, energy_only, exit_inner_loop, &
      exit_outer_loop, failure, gapw, gapw_xc, harris_flag, has_unit_metric, &
      inner_loop_converged, just_energy, outer_loop_converged, scp_nddo
    REAL(KIND=dp)                            :: t1, t2
    REAL(KIND=dp), DIMENSION(3)              :: res_val_3
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: rho_ao_kp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_result_type), POINTER            :: results
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho, rho_xc
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER         :: dft_section, harris_section, &
                                                input, scf_section

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(dft_control,rho,rho_xc,energy, &
            logger,qs_charges,ks_env,mos,atomic_kind_set,qs_kind_set, &
            particle_set,harris_env,dft_section,input,&
            scf_section, para_env, results, kpoints, pw_env, rho_ao_kp)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)

    logger => cp_error_get_logger(error)
    t1 = m_walltime()

    CALL get_qs_env(qs_env=qs_env,&
                    energy=energy,&
                    particle_set=particle_set,&
                    qs_charges=qs_charges,&
                    ks_env=ks_env, &
                    harris_env=harris_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    rho=rho,&
                    rho_xc=rho_xc,&
                    mos=mos,&
                    input=input,&
                    dft_control=dft_control,&
                    do_kpoints=do_kpoints,&
                    kpoints=kpoints,&
                    results=results,&
                    pw_env=pw_env,&
                    para_env=para_env,&
                    error=error)

    CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp, error=error)
    ! Defining SCP logicals
    scp_nddo = dft_control%qs_control%se_control%scp

    dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
    scf_section => section_vals_get_subs_vals(dft_section,"SCF",error=error)
    harris_section => section_vals_get_subs_vals(dft_section,"QS%HARRIS", error=error)

    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".scfLog",error=error)

    IF (output_unit>0) WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") &
            "SCF WAVEFUNCTION OPTIMIZATION"

    ! short cut flags setting the different methods
    gapw = dft_control%qs_control%gapw
    gapw_xc = dft_control%qs_control%gapw_xc
    ! **** SCP
    IF ( scp_nddo ) THEN
       CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error=error )
       IF (dft_control%qs_control%method_id.EQ.do_method_mndo) THEN
          ! this combination is buggy since revision 12141
          CALL stop_program(routineN,moduleN,__LINE__,"Not Yet Implemented (buggy since rev12141)")
       ENDIF
    END IF

    harris_flag = qs_env%use_harris
    CALL get_qs_env(qs_env,has_unit_metric=has_unit_metric,error=error)

    IF ((output_unit > 0).AND.(.NOT.scf_control%use_ot)) THEN
       WRITE (UNIT=output_unit,&
              FMT="(/,T3,A,T12,A,T31,A,T39,A,T59,A,T75,A,/,T3,A)")&
         "Step","Update method","Time","Convergence","Total energy","Change",&
         REPEAT("-",78)
    END IF
    CALL cp_add_iter_level(logger%iter_info,"QS_SCF",error=error)

    ! check for external communicator and if the inermediate energy should be sended
    res_val_3(:) = -1.0_dp
    description  = "[EXT_SCF_ENER_COMM]"
    IF(test_for_result(results,description=description, error=error)) THEN
      CALL get_results(results, description=description,&
                       values=res_val_3, n_entries=i_tmp, error=error)
      CPPostcondition(i_tmp.EQ.3,cp_failure_level,routineP,error,failure)
      CALL cp_assert(ANY(res_val_3(:).GT.0.0),cp_failure_level,&
             cp_assertion_failed,routineP,&
             " Trying to access result ("//TRIM(description)//&
             ") which is not correctly stored. ",&
             error, only_ionode=.TRUE.)
    END IF
    external_comm          = NINT(res_val_3(1))
    ext_master_id          = NINT(res_val_3(2))
    scf_energy_message_tag = NINT(res_val_3(3))

    ! *** outer loop of the scf, can treat other variables,
    ! *** such as lagrangian multipliers
    scf_env%outer_scf%iter_count=0
    iter_count = 0
    total_steps = 0
    energy%tot_old = 0.0_dp

    scf_outer_loop: DO

       CALL init_scf_loop(scf_env=scf_env, qs_env=qs_env, &
            scf_section=scf_section, error=error)

       CALL qs_scf_set_loop_flags(scp_nddo,scf_control,scf_env,diis_step,&
                      energy_only,just_energy,exit_inner_loop,error=error)

       scf_loop: DO

          CALL timeset(routineN//"_inner_loop",handle2)

          scf_env%iter_count = scf_env%iter_count + 1
          iter_count =  iter_count + 1
          CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=iter_count,error=error)

          IF (output_unit > 0) CALL m_flush(output_unit)

          total_steps = total_steps + 1
          just_energy = energy_only
 
          CALL qs_ks_update_qs_env(qs_env, just_energy=just_energy,&
                                   calculate_forces=.FALSE., error=error)

          ! print 'heavy weight' or relatively expensive quantities
          CALL qs_scf_loop_print(qs_env,scf_env,para_env,error=error)

          IF(do_kpoints) THEN
             ! kpoints
             CALL qs_scf_new_mos_kp(qs_env,scf_env,scf_control,scf_section,diis_step,energy_only,error=error)
          ELSE
             ! Gamma points only
             CALL qs_scf_new_mos(qs_env,scf_env,scf_control,scf_section,diis_step,energy_only,error=error)
          END IF

          ! another heavy weight print object, print controlled by dft_section
          CALL qs_scf_write_mos(mos,atomic_kind_set,qs_kind_set,particle_set,dft_section,.FALSE.,error=error)
          
          CALL qs_scf_harris_e_correct(qs_env,harris_env,para_env,scf_env,energy,harris_section,&
                                       harris_flag,error)

          CALL qs_scf_density_mixing(scf_env,rho,para_env,diis_step,do_kpoints,error=error)

          t2 = m_walltime()

          CALL qs_scf_loop_info(scf_env,output_unit,just_energy,t1,t2,energy)

          IF (.NOT.just_energy) energy%tot_old = energy%total

          IF (harris_flag) CALL harris_postprocessing(harris_env, error=error)

          ! check for external communicator and if the intermediate energy should be sent
          IF(scf_energy_message_tag.GT.0)THEN
            CALL mp_send(energy%total, ext_master_id, scf_energy_message_tag, external_comm)
          END IF

          CALL qs_scf_check_inner_exit(qs_env,scf_env,scf_control,should_stop,exit_inner_loop,&
                               inner_loop_converged,output_unit,error)

          ! In case we decide to exit we perform few more check to see if this one
          ! is really the last SCF step
          IF (exit_inner_loop) THEN

             CALL qs_scf_inner_finalize(scf_env,scf_control,qs_env,harris_flag,diis_step,output_unit,error)

             CALL qs_scf_check_outer_exit(qs_env,scf_env,scf_control,should_stop,&
                                         outer_loop_converged,exit_outer_loop,error)

             ! Let's tag the last SCF cycle so we can print informations only of the last step
             IF (exit_outer_loop) CALL cp_iterate(logger%iter_info,last=.TRUE.,iter_nr=iter_count,error=error)

          END IF

          ! Write WaveFunction restart file
          CALL write_mo_set(mos,particle_set,dft_section=dft_section,&
                            atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,&
                            kpoints=kpoints,error=error)

          ! Exit if we have finished with the SCF inner loop
          IF (exit_inner_loop) THEN
             CALL timestop(handle2)
             EXIT scf_loop
          END IF

          IF (.NOT.BTEST(cp_print_key_should_output(logger%iter_info,&
               scf_section,"PRINT%ITERATION_INFO/TIME_CUMUL",error=error),cp_p_file)) &
               t1 = m_walltime()

          ! mixing methods have the new density matrix in p_mix_new
          IF (scf_env%mixing_method > 0) THEN
             DO ic=1,SIZE(rho_ao_kp,2)
                DO ispin=1,dft_control%nspins
                   CALL cp_dbcsr_copy(rho_ao_kp(ispin,ic)%matrix,scf_env%p_mix_new(ispin,ic)%matrix,&
                                   name=TRIM(cp_dbcsr_name(rho_ao_kp(ispin,ic)%matrix)),error=error)
                END DO
             END DO
          END IF

          CALL qs_scf_rho_update(rho,qs_env,scf_env,ks_env,&
                    mix_rho=scf_env%mixing_method>=gspace_mixing_nr,error=error)

          CALL timestop(handle2)

       END DO scf_loop

       IF (.NOT. scf_control%outer_scf%have_scf) EXIT scf_outer_loop

       ! In case we use the OUTER SCF loop let's print some info..
       CALL qs_scf_outer_loop_info(output_unit,scf_control,scf_env,&
                                   energy,total_steps,should_stop,outer_loop_converged)

       IF (exit_outer_loop) EXIT scf_outer_loop

       CALL outer_loop_optimize(scf_env,scf_control,error)
       CALL outer_loop_update_qs_env(qs_env,scf_env,error)
       CALL qs_ks_did_change(ks_env,potential_changed=.TRUE.,error=error)

    END DO scf_outer_loop

    converged = inner_loop_converged .AND. outer_loop_converged

    ! if needed copy mo_coeff dbcsr->fm for later use in post_scf!fm->dbcsr
    DO ispin=1,SIZE(mos)!fm -> dbcsr
       IF(mos(ispin)%mo_set%use_mo_coeff_b) THEN!fm->dbcsr
          IF (.NOT.ASSOCIATED(mos(ispin)%mo_set%mo_coeff_b)) &!fm->dbcsr
             CALL stop_program(routineN,moduleN,__LINE__,"mo_coeff_b is not allocated")!fm->dbcsr
          CALL copy_dbcsr_to_fm(mos(ispin)%mo_set%mo_coeff_b, &!fm->dbcsr
               mos(ispin)%mo_set%mo_coeff,error=error)!fm -> dbcsr
       ENDIF!fm->dbcsr
    ENDDO!fm -> dbcsr


    IF(dft_control%qs_control%becke_restraint)THEN
       CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,&
            dft_control%qs_control%becke_control%becke_pot%pw,error=error)
       dft_control%qs_control%becke_control%need_pot=.TRUE.
    END IF
    CALL cp_rm_iter_level(logger%iter_info,level_name="QS_SCF",error=error)
    CALL timestop(handle)

  END SUBROUTINE scf_env_do_scf

! *****************************************************************************
!> \brief inits those objects needed if you want to restart the scf with, say
!>        only a new initial guess, or different density functional or ...
!>        this will happen just before the scf loop starts
!> \param scf_env ...
!> \param qs_env ...
!> \param scf_section ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE init_scf_loop(scf_env,qs_env,scf_section,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_scf_loop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nmo, &
                                                number_of_OT_envs, stat
    LOGICAL                                  :: do_rotation, failure, &
                                                has_unit_metric, is_full_all, &
                                                scp_nddo
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: ks_scp, orthogonality_metric, &
                                                pscp
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(scp_environment_type), POINTER      :: scp_env

    CALL timeset(routineN,handle)

    NULLIFY(scf_control,matrix_s,matrix_ks,dft_control,mos,mo_coeff)
! **SCP
    NULLIFY (scp_env, pscp, ks_scp)
! **SCP

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,&
         scf_control=scf_control,&
         dft_control=dft_control,&
         mos=mos,error=error)

    scp_nddo = dft_control%qs_control%se_control%scp.AND.&
              (dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb)

    ! if using mo_coeff_b then copy to fm
    DO ispin=1,SIZE(mos)!fm->dbcsr
       IF(mos(1)%mo_set%use_mo_coeff_b)THEN!fm->dbcsr
          CALL copy_dbcsr_to_fm(mos(ispin)%mo_set%mo_coeff_b,mos(ispin)%mo_set%mo_coeff,error=error)!fm->dbcsr
       ENDIF!fm->dbcsr
    ENDDO!fm->dbcsr

    ! this just guarantees that all mo_occupations match the eigenvalues, if smear
    DO ispin=1,dft_control%nspins
       CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,&
                              smear=scf_control%smear, error=error)
    ENDDO

    SELECT CASE (scf_env%method)
    CASE DEFAULT

      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method:"//&
               cp_to_string(scf_env%method),error,failure)

    CASE (general_diag_method_nr,special_diag_method_nr,block_krylov_diag_method_nr)
       IF(.NOT.scf_env%skip_diis) THEN
         IF (.NOT.ASSOCIATED(scf_env%scf_diis_buffer)) THEN
            CALL qs_diis_b_create(scf_env%scf_diis_buffer,nbuffer=scf_control%max_diis,error=error)
         END IF
         CALL qs_diis_b_clear(scf_env%scf_diis_buffer,error=error)
       END IF

    CASE (ot_diag_method_nr)
       CALL get_qs_env(qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s, error=error)

       IF(.NOT.scf_env%skip_diis) THEN
         IF (.NOT.ASSOCIATED(scf_env%scf_diis_buffer)) THEN
            CALL qs_diis_b_create(scf_env%scf_diis_buffer,nbuffer=scf_control%max_diis,error=error)
         END IF
         CALL qs_diis_b_clear(scf_env%scf_diis_buffer,error=error)
       END IF

       ! disable DFTB and SE for now
       IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%semi_empirical) THEN
         CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"DFTB and SE not available with OT/DIAG",error,failure)
       END IF

       ! if an old preconditioner is still around (i.e. outer SCF is active),
       ! remove it if this could be worthwhile
       CALL restart_preconditioner(qs_env,scf_env%ot_preconditioner,&
            scf_control%diagonalization%ot_settings%preconditioner_type,&
            dft_control%nspins,error=error)

       CALL prepare_preconditioner(qs_env,mos,matrix_ks,matrix_s,scf_env%ot_preconditioner,&
            scf_control%diagonalization%ot_settings%preconditioner_type,&
            scf_control%diagonalization%ot_settings%precond_solver_type,&
            scf_control%diagonalization%ot_settings%energy_gap,dft_control%nspins,error=error)

    CASE (block_davidson_diag_method_nr)
        ! Preconditioner initialized within the loop, when required
    CASE (ot_method_nr)
       CALL get_qs_env(qs_env,&
                       has_unit_metric=has_unit_metric,&
                       matrix_s=matrix_s,&
                       matrix_ks=matrix_ks,&
                       error=error)

       ! reortho the wavefunctions if we are having an outer scf and
       ! this is not the first iteration
       ! this is useful to avoid the build-up of numerical noise
       ! however, we can not play this trick if restricted (don't mix non-equivalent orbs)
       IF(scf_control%do_outer_scf_reortho) THEN
          IF (scf_control%outer_scf%have_scf .AND. .NOT. dft_control%restricted) THEN
             IF (scf_env%outer_scf%iter_count>0) THEN
                DO ispin=1,dft_control%nspins
                   CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo)
                   IF (has_unit_metric) THEN
                      CALL make_basis_simple(mo_coeff,nmo,error=error)
                   ELSE
                      CALL make_basis_sm(mo_coeff,nmo,matrix_s(1)%matrix,error=error)
                   ENDIF
                ENDDO
             ENDIF
          ENDIF
       ELSE
          ! dont need any dirty trick for the numerically stable irac algorithm.
       ENDIF

       IF (.NOT.ASSOCIATED(scf_env%qs_ot_env)) THEN

          ! restricted calculations require just one set of OT orbitals
          number_of_OT_envs=dft_control%nspins
          IF (dft_control%restricted) number_of_OT_envs=1

          ALLOCATE(scf_env%qs_ot_env(number_of_OT_envs),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          ! XXX Joost XXX should disentangle reading input from this part
          CALL ot_scf_read_input(scf_env%qs_ot_env,scf_section,error)

          ! keep a note that we are restricted
          IF (dft_control%restricted) THEN
              scf_env%qs_ot_env(1)%restricted=.TRUE.
              ! requires rotation
              CALL cp_assert(scf_env%qs_ot_env(1)%settings%do_rotation,cp_fatal_level,&
                   cp_assertion_failed,routineP,&
                   "Restricted calculation with OT requires orbital rotation. Please "//&
                   "activate the OT%ROTATION keyword! "//&
CPSourceFileRef)
          ELSE
              scf_env%qs_ot_env(:)%restricted=.FALSE.
          ENDIF

          ! might need the KS matrix to init properly
          CALL qs_ks_update_qs_env(qs_env, just_energy=.FALSE.,&
                                   calculate_forces=.FALSE., error=error)

          ! if an old preconditioner is still around (i.e. outer SCF is active),
          ! remove it if this could be worthwhile
          CALL restart_preconditioner(qs_env,scf_env%ot_preconditioner,&
               scf_env%qs_ot_env(1)%settings%preconditioner_type,&
               dft_control%nspins,error=error)

          !
          ! preconditioning still needs to be done correctly with has_unit_metric
          ! notice that a big part of the preconditioning (S^-1) is fine anyhow
          !
          IF (has_unit_metric) THEN
             NULLIFY(orthogonality_metric)
          ELSE
             orthogonality_metric=>matrix_s(1)%matrix
          ENDIF
          IF(dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb) THEN
             scf_env%qs_ot_env(1)%settings%mixed_precision=.FALSE.
          END IF

          CALL prepare_preconditioner(qs_env,mos,matrix_ks,matrix_s,scf_env%ot_preconditioner,&
               scf_env%qs_ot_env(1)%settings%preconditioner_type,&
               scf_env%qs_ot_env(1)%settings%precond_solver_type,&
               scf_env%qs_ot_env(1)%settings%energy_gap,dft_control%nspins,&
               has_unit_metric=has_unit_metric,&
               mixed_precision=scf_env%qs_ot_env(1)%settings%mixed_precision,error=error)

          ! we have to treat the semi-empirical methods separately
          IF ( scp_nddo ) THEN
             CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error = error )
             CALL get_scp_env ( scp_env = scp_env, pscp=pscp, ks_scp=ks_scp, error = error )
             CALL ot_scf_init(mos,orthogonality_metric, pscp=pscp, fscp=ks_scp, &
                  broyden_adaptive_sigma=qs_env%broyden_adaptive_sigma,&
                  qs_ot_env=scf_env%qs_ot_env, matrix_ks=matrix_ks(1)%matrix, &
                  error=error)
          ELSE
             CALL ot_scf_init(mo_array = mos, matrix_s = orthogonality_metric, &
                  broyden_adaptive_sigma=qs_env%broyden_adaptive_sigma,&
                  qs_ot_env = scf_env%qs_ot_env,matrix_ks=matrix_ks(1)%matrix,&
                  error=error)
          END IF

          SELECT CASE(scf_env%qs_ot_env(1)%settings%preconditioner_type)
             CASE(ot_precond_none)
             CASE(ot_precond_full_all,ot_precond_full_single_inverse)
                ! this will rotate the MOs to be eigen states, which is not compatible with rotation
                ! e.g. mo_derivs here do not yet include potentially different occupations numbers
                do_rotation=scf_env%qs_ot_env(1)%settings%do_rotation
                ! only full all needs rotation
                is_full_all=scf_env%qs_ot_env(1)%settings%preconditioner_type==ot_precond_full_all
                CPPrecondition(.NOT.(do_rotation.AND.is_full_all),cp_failure_level,routineP,error,failure)
                DO ispin=1,SIZE(scf_env%qs_ot_env)
                   CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin),&
                        scf_env%ot_preconditioner(ispin)%preconditioner,error=error)
                ENDDO
             CASE(ot_precond_s_inverse,ot_precond_full_single)
                DO ispin=1,SIZE(scf_env%qs_ot_env)
                   CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin),&
                        scf_env%ot_preconditioner(1)%preconditioner,error=error)
                ENDDO
             CASE DEFAULT
                DO ispin=1,SIZE(scf_env%qs_ot_env)
                   CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin),&
                        scf_env%ot_preconditioner(1)%preconditioner,error=error)
                ENDDO
          END SELECT
       ENDIF

       ! if we have non-uniform occupations we should be using rotation
       do_rotation=scf_env%qs_ot_env(1)%settings%do_rotation
       DO ispin=1,SIZE(mos)
          IF (.NOT. mos(ispin)%mo_set%uniform_occupation) THEN
             CPPrecondition(do_rotation,cp_failure_level,routineP,error,failure)
          ENDIF
       ENDDO
    END SELECT

    ! another safety check
    IF (dft_control%low_spin_roks) THEN
       CPPrecondition(scf_env%method==ot_method_nr,cp_failure_level,routineP,error,failure)
       do_rotation=scf_env%qs_ot_env(1)%settings%do_rotation
       CPPrecondition(do_rotation,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE init_scf_loop


! *****************************************************************************
!> \brief perform cleanup operations (like releasing temporary storage)
!>      at the end of the scf
!> \param scf_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf_env_cleanup(scf_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN

!   *** Release SCF work storage ***

       IF (ASSOCIATED(scf_env%scf_work1)) THEN
          DO ispin=1,SIZE(scf_env%scf_work1)
             CALL cp_fm_release(scf_env%scf_work1(ispin)%matrix,error=error)
          ENDDO
          DEALLOCATE(scf_env%scf_work1)
       ENDIF
       IF (ASSOCIATED(scf_env%scf_work2)) CALL cp_fm_release(scf_env%scf_work2,error)
       IF (ASSOCIATED(scf_env%ortho)) CALL cp_fm_release(scf_env%ortho,error=error)
       IF (ASSOCIATED(scf_env%ortho_m1)) CALL cp_fm_release(scf_env%ortho_m1,error=error)

       IF (ASSOCIATED(scf_env%ortho_dbcsr)) THEN
           CALL cp_dbcsr_deallocate_matrix(scf_env%ortho_dbcsr,error=error)
       END IF
       IF (ASSOCIATED(scf_env%buf1_dbcsr)) THEN
           CALL cp_dbcsr_deallocate_matrix(scf_env%buf1_dbcsr,error=error)
       END IF
       IF (ASSOCIATED(scf_env%buf2_dbcsr)) THEN
           CALL cp_dbcsr_deallocate_matrix(scf_env%buf2_dbcsr,error=error)
       END IF

       IF (ASSOCIATED(scf_env%p_mix_new)) THEN
          CALL cp_dbcsr_deallocate_matrix_set(scf_env%p_mix_new,error=error)
       END IF

       IF (ASSOCIATED(scf_env%p_delta)) THEN
          CALL cp_dbcsr_deallocate_matrix_set(scf_env%p_delta,error=error)
       END IF

! *** method dependent cleanup
       SELECT CASE(scf_env%method)
       CASE(ot_method_nr)
          !
       CASE(ot_diag_method_nr)
          !
       CASE(general_diag_method_nr)
          !
       CASE(special_diag_method_nr)
          !
       CASE(block_krylov_diag_method_nr)
       CASE(block_davidson_diag_method_nr)
           CALL block_davidson_deallocate(scf_env%block_davidson_env,error=error)
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method:"//&
               cp_to_string(scf_env%method),error,failure)
       END SELECT

       IF (ASSOCIATED(scf_env%outer_scf%variables)) THEN
          DEALLOCATE(scf_env%outer_scf%variables,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
       IF (ASSOCIATED(scf_env%outer_scf%count)) THEN
          DEALLOCATE(scf_env%outer_scf%count,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
       IF (ASSOCIATED(scf_env%outer_scf%gradient)) THEN
          DEALLOCATE(scf_env%outer_scf%gradient,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF
       IF (ASSOCIATED(scf_env%outer_scf%energy)) THEN
          DEALLOCATE(scf_env%outer_scf%energy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ENDIF

       IF (scf_env%do_diag_sub) THEN
          IF(ASSOCIATED(scf_env%subspace_env%p_matrix_store)) THEN
            CALL cp_dbcsr_deallocate_matrix_set(scf_env%subspace_env%p_matrix_store,error=error)
          END IF
          IF(ASSOCIATED(scf_env%subspace_env%p_matrix_mix)) THEN
            CALL cp_dbcsr_deallocate_matrix_set(scf_env%subspace_env%p_matrix_mix,error=error)
          END IF
       END IF

    END IF
    CALL timestop(handle)

  END SUBROUTINE scf_env_cleanup

END MODULE qs_scf
