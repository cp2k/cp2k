!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
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
!> \author Matthias Krack (30.04.2001)
! *****************************************************************************
MODULE qs_scf

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_symm,&
                                             cp_fm_triangular_invert,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_power,&
                                             cp_fm_syevd
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_get_el_struct
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_add_to_element,&
                                             cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_fm_to_sm,&
                                             copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type,&
                                             sm_pools_create_matrix_vect,&
                                             sm_pools_flush_cache,&
                                             sm_pools_give_back_matrix_vect
  USE f77_blas
  USE global_types,                    ONLY: global_environment_type
  USE harris_env_types,                ONLY: harris_env_type
  USE harris_functional,               ONLY: harris_eigenvalue_calculation,&
                                             harris_eigenvalue_trace_KS_Pmix,&
                                             harris_energy_correction,&
                                             harris_postprocessing
  USE input_constants,                 ONLY: &
       broy_mix, broy_mix_new, cholesky_inverse, cholesky_off, &
       cholesky_reduce, cholesky_restore, core_guess, debug_run, &
       densities_guess, diag_ot, diag_standard, direct_p_mix, general_roks, &
       high_spin_roks, history_guess, kerker_mix, no_mix, &
       ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
       ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, &
       ot_precond_sparse_diag, ot_precond_sparse_firstorder, &
       ot_precond_sparse_kinetic, outer_scf_none, outer_scf_scp, &
       plus_u_lowdin, pulay_mix, restart_guess, wfi_frozen_method_nr, &
       wfi_use_guess_method_nr
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE message_passing,                 ONLY: mp_min
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: evolt,&
                                             kcalmol
  USE preconditioner,                  ONLY: make_preconditioner
  USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                             init_preconditioner
  USE pw_env_types,                    ONLY: pw_env_get
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_type
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                             qs_diis_b_create,&
                                             qs_diis_b_step
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_reorthogonalize_vectors,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_gspace_mixing,                ONLY: broyden_mixing,&
                                             broyden_mixing_new,&
                                             gmix_potential_only,&
                                             pulay_mixing,&
                                             self_consistency_check
  USE qs_initial_guess,                ONLY: calculate_first_density_matrix
  USE qs_ks_methods,                   ONLY: qs_ks_create,&
                                             qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_scp_methods,               ONLY: qs_ks_scp_did_change,&
                                             qs_ks_scp_update
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             calculate_magnitude,&
                                             calculate_orthonormality,&
                                             calculate_subspace_eigenvalues,&
                                             make_basis_simple,&
                                             make_basis_sm
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type,&
                                             set_mo_occupation,&
                                             set_mo_set,&
                                             write_mo_set
  USE qs_ot,                           ONLY: qs_ot_new_preconditioner
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_ot_scf,                       ONLY: ot_scf_destroy,&
                                             ot_scf_init,&
                                             ot_scf_mini,&
                                             ot_scf_read_input
  USE qs_outer_scf,                    ONLY: outer_loop_extrapolate,&
                                             outer_loop_gradient,&
                                             outer_loop_optimize,&
                                             outer_loop_update_qs_env,&
                                             outer_loop_variables_count
  USE qs_rho_methods,                  ONLY: duplicate_rho_type,&
                                             qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_methods,                  ONLY: combine_ks_matrices,&
                                             eigensolver,&
                                             scf_env_density_mixing,&
                                             simple_eigensolver
  USE qs_scf_post_gpw,                 ONLY: scf_post_calculation_gpw
  USE qs_scf_post_se,                  ONLY: scf_post_calculation_se
  USE qs_scf_types,                    ONLY: &
       broyden_mixing_new_nr, broyden_mixing_nr, direct_mixing_nr, &
       general_diag_method_nr, gspace_mixing_nr, no_mixing_nr, &
       ot_diag_method_nr, ot_method_nr, pulay_mixing_nr, qs_scf_env_type, &
       scf_env_create, scf_env_release, special_diag_method_nr
  USE qs_wf_history_methods,           ONLY: wfi_extrapolate,&
                                             wfi_get_method_label,&
                                             wfi_update
  USE scf_control_types,               ONLY: scf_control_type,&
                                             smear_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type
  USE scp_density_methods,             ONLY: update_rhoscp
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type,&
                                             set_scp_env
  USE sparse_matrix_output,            ONLY: write_fm_with_basis_info,&
                                             write_sparse_matrix
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             real_matrix_type,&
                                             transfer_matrix
  USE termination,                     ONLY: external_control,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
  USE xas_restart,                     ONLY: xas_initialize_rho
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf'

  PUBLIC :: scf, scf_env_cleanup, scf_env_do_scf,&
            init_scf_run, init_scf_loop, &
            qs_scf_print_summary

CONTAINS

! *****************************************************************************
!> \brief perform an scf procedure in the given qs_env
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      globenv should be removed
!> \par History
!>      02.2003 introduced scf_env, moved real work to scf_env_do_scf [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf(qs_env,globenv,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: max_scf_tmp, output_unit, stat
    LOGICAL                                  :: converged, failure, ionode, &
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
    ionode = qs_env%para_env%ionode
    output_unit = -1
    IF(ionode) output_unit = cp_logger_get_default_unit_nr(logger)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL get_qs_env(qs_env,scf_env=scf_env,error=error,input=input, &
                       dft_control=dft_control, scf_control=scf_control)

       dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
       scf_section => section_vals_get_subs_vals(dft_section,"SCF",error=error)

       IF (.NOT.ASSOCIATED(scf_env)) THEN ! i.e. for MD this is associated on the second step (it so seems)
          CALL scf_env_create(scf_env, error=error)
          CALL set_qs_env(qs_env,scf_env=scf_env,error=error)
          CALL scf_env_release(scf_env,error=error)

          ! set some of the methods that might be used in this SCF.
          ! this might not yet be the ideal place to set this kind of stuff (who knows?)
          CALL get_qs_env(qs_env,scf_env=scf_env,scf_control=scf_control, &
                          dft_control=dft_control, error=error)

          IF (scf_control%use_diag) THEN
            IF ( scf_control%diagonalization%method == diag_standard ) THEN
              scf_env%method=general_diag_method_nr
              IF (dft_control%qs_control%semi_empirical) THEN
                IF (dft_control%qs_control%se_control%orthogonal_basis) &
                  scf_env%method=special_diag_method_nr
              END IF
              IF (dft_control%qs_control%dftb) THEN
                IF (dft_control%qs_control%dftb_control%orthogonal_basis) &
                  scf_env%method=special_diag_method_nr
              END IF
            ELSEIF ( scf_control%diagonalization%method == diag_ot ) THEN
              scf_env%method=ot_diag_method_nr
              IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%semi_empirical) THEN
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"DFTB and SE not possible with OT diagonalization",error,failure)
              END IF
            ELSE
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"Unknown diagonalization method",error,failure)
            END IF
          ELSEIF (scf_control%use_ot) THEN
            scf_env%method=ot_method_nr
          ELSE
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"OT or DIAGONALIZATION have to be set",error,failure)
          END IF

          SELECT CASE(scf_control%mixing_method)
          CASE(no_mix)
             scf_env%mixing_method=no_mixing_nr
             scf_env%p_mix_alpha = 1.0_dp
          CASE(direct_p_mix) 
             scf_env%mixing_method=direct_mixing_nr
             CALL section_vals_val_get(scf_section,"MIXING%ALPHA",r_val=scf_env%p_mix_alpha,error=error)
          CASE(kerker_mix)
             scf_env%mixing_method=gspace_mixing_nr
             CALL section_vals_val_get(scf_section,"MIXING%ALPHA",r_val=scf_env%p_mix_alpha,error=error)
             CALL section_vals_val_get(scf_section,"MIXING%BETA",r_val=scf_env%p_mix_beta,error=error)
          CASE(pulay_mix)
             scf_env%mixing_method=pulay_mixing_nr
             CALL section_vals_val_get(scf_section,"MIXING%ALPHA",r_val=scf_env%p_mix_alpha,error=error)
             CALL section_vals_val_get(scf_section,"MIXING%BETA",r_val=scf_env%p_mix_beta,error=error)
             CALL section_vals_val_get(scf_section,"MIXING%NPULAY",i_val=scf_env%npulay,error=error)
             CALL section_vals_val_get(scf_section,"MIXING%PULAY_ALPHA",r_val=scf_env%pulay_alpha,error=error)
          CASE(broy_mix)
             scf_env%mixing_method=broyden_mixing_nr
             CALL section_vals_val_get(scf_section,"MIXING%ALPHA",r_val=scf_env%p_mix_alpha,error=error)
             CALL section_vals_val_get(scf_section,"MIXING%BETA",r_val=scf_env%p_mix_beta,error=error)
             CALL section_vals_val_get(scf_section,"MIXING%NBROYDEN",i_val=scf_env%nbroy,error=error)
             CALL section_vals_val_get(scf_section,"MIXING%BROY_W0",r_val=scf_env%broy_w0,error=error)
          CASE(broy_mix_new)
             scf_env%mixing_method=broyden_mixing_new_nr
             CALL section_vals_val_get(scf_section,"MIXING%ALPHA",r_val=scf_env%p_mix_alpha,error=error)
             CALL section_vals_val_get(scf_section,"MIXING%BETA",r_val=scf_env%p_mix_beta,error=error)
             CALL section_vals_val_get(scf_section,"MIXING%NBROYDEN",i_val=scf_env%nbroy,error=error)
             ALLOCATE (scf_env%broy_store,STAT=stat)
             scf_env%broy_store%nbuffer = scf_env%nbroy
             CALL section_vals_val_get(scf_section,"MIXING%BROY_WREF",r_val=scf_env%broy_store%wc,error=error)
             CALL section_vals_val_get(scf_section,"MIXING%BROY_WMAX",r_val=scf_env%broy_store%wmax,error=error)
             scf_env%broy_store%bconst = 20.0_dp
             scf_env%broy_store%p_metric_method=1
          CASE DEFAULT        
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"Unknown mixing method",error,failure)
          END SELECT
          IF( scf_env%method==ot_method_nr ) THEN
            ! No mixing is used with OT
             scf_env%mixing_method=no_mixing_nr
             scf_env%p_mix_alpha = 1.0_dp
             scf_env%skip_diis = .TRUE.
          END IF


          IF(scf_env%mixing_method>direct_mixing_nr) THEN
            scf_env%skip_diis = .TRUE.
            IF(scf_env%p_mix_beta==0.0_dp) THEN
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                     routineP,"Mixing employing the Kerker damping factor: choose BETA /= 0.0",error,failure)
            END IF
          END IF

          IF(scf_env%mixing_method==direct_mixing_nr) THEN
            IF(qs_env%scf_control%eps_diis<qs_env%scf_control%eps_scf) THEN
             scf_env%skip_diis = .TRUE.
             IF(output_unit>0)  WRITE(output_unit,*) "WARNING the diis scheme is disabled "//&
                 "since eps_diis<eps_scf"
            END IF
            CALL section_vals_val_get(scf_section,"MIXING%NSKIP",i_val=scf_env%nskip_mixing,error=error)
          END IF

          CALL section_vals_val_get(scf_section,"CHOLESKY",i_val=scf_env%cholesky_method,error=error)
       END IF

       IF ( qs_env%scf_control%max_scf > 0 ) THEN

         CALL init_scf_run(scf_env=scf_env, qs_env=qs_env, scf_section=scf_section, error=error)

         IF ((qs_env%scf_control%density_guess .EQ. history_guess) .AND. (first_step_flag)) THEN
           max_scf_tmp = qs_env%scf_control%max_scf
           qs_env%scf_control%max_scf = 1
           outer_scf_loop = qs_env%scf_control%outer_scf%have_scf
           qs_env%scf_control%outer_scf%have_scf = .FALSE.
         END IF
         CALL scf_env_do_scf(scf_env=scf_env, qs_env=qs_env, globenv=globenv, &
                             converged=converged, should_stop=should_stop, error=error)

         !   *** add the converged wavefunction to the wavefunction history
         IF ((ASSOCIATED(qs_env%wf_history)) .AND. &
             ((qs_env%scf_control%density_guess .NE. history_guess) .OR. &
              (.NOT. first_step_flag))) THEN
             CALL wfi_update(qs_env%wf_history,qs_env=qs_env,dt=1.0_dp, error=error)
         ELSE IF ((qs_env%scf_control%density_guess .EQ. history_guess) .AND. &
                  (first_step_flag)) THEN
           qs_env%scf_control%max_scf = max_scf_tmp
           qs_env%scf_control%outer_scf%have_scf = outer_scf_loop
           first_step_flag = .FALSE.
         END IF

         ! *** compute properties that depend on the converged wavefunction
         IF (.NOT.(should_stop)) THEN
            IF     (dft_control%qs_control%semi_empirical) THEN
               CALL scf_post_calculation_se (dft_section, scf_env, qs_env, error)
            ELSEIF (dft_control%qs_control%dftb) THEN
               ! Do Nothing
            ELSE
               CALL scf_post_calculation_gpw(dft_section, scf_env, qs_env, error)
            END IF
         END IF

       END IF

       ! *** cleanup
       CALL scf_env_cleanup(scf_env,qs_env=qs_env,error=error)

    END IF

  END SUBROUTINE scf

! *****************************************************************************
!> \brief perform an scf loop
!> \param scf_env the scf_env where to perform the scf procedure
!> \param qs_env the qs_env, the scf_env lives in
!> \param converged will be true / false if converged is reached
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      globenv should be removed
!> \par History
!>      long history, see cvs and qs_scf module history
!>      02.2003 introduced scf_env [fawzi]
!>      09.2005 Frozen density approximation [TdK]
!>      06.2007 Check for SCF iteration count early [jgh]
!> \author Matthias Krack
! *****************************************************************************
  SUBROUTINE scf_env_do_scf(scf_env,qs_env,globenv,converged,should_stop,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(OUT)                     :: converged, should_stop
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scf_env_do_scf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, ispin, &
                                                iter_count, output_unit, &
                                                total_steps
    LOGICAL :: diis_step, energy_only, exit_inner_loop, failure, gapw, &
      gapw_xc, harris_energy, harris_flag, has_unit_metric, &
      outer_loop_converged, scp_dft, scp_nddo, use_jacobi, use_virial
    REAL(KIND=dp)                            :: delta_rhog, delta_rhog_p, &
                                                outer_loop_eps, t1, t2
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho, rho_xc
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(real_matrix_type), POINTER          :: ks_scp, pscp
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER         :: dft_section, harris_section, &
                                                input, scf_section
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    failure=.FALSE.
    converged=.TRUE.

    NULLIFY(dft_control,matrix_s,matrix_ks,rho,rho_xc,energy, &
            scf_control,logger,qs_charges,ks_env,mos,atomic_kind_set, &
            particle_set,harris_env,dft_section,input,&
            scf_section)
    NULLIFY ( aux_coeff_set, scp_env, pscp, ks_scp)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    para_env=>qs_env%para_env

    logger => cp_error_get_logger(error)
    t1 = m_walltime()

    CALL get_qs_env(qs_env=qs_env,&
                    matrix_s=matrix_s,energy=energy,&
                    particle_set=particle_set,&
                    qs_charges=qs_charges,&
                    ks_env=ks_env, &
                    harris_env=harris_env,&
                    atomic_kind_set=atomic_kind_set,&
                    matrix_ks=matrix_ks,rho=rho,rho_xc=rho_xc,mos=mos, &
                    input=input, dft_control=dft_control, scf_control=scf_control, &
                    virial=virial, error=error)
    ! Defining SCP logicals
    scp_dft = dft_control%scp
    scp_nddo = dft_control%qs_control%se_control%scp 

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

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
    IF ( ( scp_dft ) .OR. ( scp_nddo ) ) THEN
       CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error=error )
    END IF
! Assign SCP-DFT structures
    IF ( scp_dft ) THEN
      CALL get_scp_env ( scp_env = scp_env, aux_coeff_set = aux_coeff_set, error = error )
    ENDIF
! Assign SCP-NDDO structures
    IF ( scp_nddo ) THEN
      CALL get_scp_env ( scp_env = scp_env, pscp=pscp, ks_scp=ks_scp, error = error )
    ENDIF
    ! **** SCP
    harris_flag = qs_env%use_harris
    has_unit_metric=.FALSE.
    IF (dft_control%qs_control%semi_empirical) &
       has_unit_metric=dft_control%qs_control%se_control%orthogonal_basis
    IF (dft_control%qs_control%dftb) &
       has_unit_metric=dft_control%qs_control%dftb_control%orthogonal_basis

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,&
            FMT="(/,T3,A,T9,A,T36,A,T49,A,T68,A,/,T3,A)")&
            "Step","Update method","Time","Convergence","Total energy",&
            REPEAT("-",77)
    END IF
    CALL cp_add_iter_level(logger%iter_info,"QS_SCF",error=error)
    ! *** outer loop of the scf, can treat other variables, 
    ! *** such as lagrangian multipliers
    scf_env%outer_scf%iter_count=0
    iter_count = 0
    total_steps = 0
    scf_outer_loop: DO

       CALL init_scf_loop(scf_env=scf_env, qs_env=qs_env, &
            scf_section=scf_section, error=error)

!   ****** Setting Some SCP flags*******
! check if we are doing SCP:
!      **** qs_ot_env%settings%scp_dft,qs_ot_env%settings%ks ARE BOTH .TRUE. here *****
       IF ( scp_dft  ) THEN
! check if we are doing outer SCF optimization of SCP coeffs:
         IF ( scf_control%outer_scf%type==outer_scf_none ) THEN
             scf_env%qs_ot_env(1)%settings%ks = .TRUE. 
            !*** qs_ot_env%settings%scp_dft*** may be set .FALSE. on input
            !    and we will preserve this if there is no outer_scf
         ELSEIF ( scf_control%outer_scf%type==outer_scf_scp ) THEN
! IF outer scf count is even, then only optimize SCP, else optimize KS
           IF ( scf_env%outer_scf%iter_count == 0 ) THEN
             scf_env%qs_ot_env(1)%settings%ks = .TRUE.
             scf_env%qs_ot_env(1)%settings%scp_dft = .FALSE.
           ELSEIF ( MOD ( scf_env%outer_scf%iter_count, 2 ) /= 0 ) THEN
             scf_env%qs_ot_env(1)%settings%ks = .FALSE.
             scf_env%qs_ot_env(1)%settings%scp_dft = .TRUE.
           ELSE
             scf_env%qs_ot_env(1)%settings%ks = .TRUE.
             scf_env%qs_ot_env(1)%settings%scp_dft = .FALSE.
           END IF
         ENDIF
       END IF
       IF ( scp_nddo  ) THEN
! check if we are doing outer SCF optimization of SCP coeffs:
         IF ( scf_control%outer_scf%type==outer_scf_none ) THEN
             scf_env%qs_ot_env(1)%settings%ks = .TRUE. 
            !*** qs_ot_env%settings%scp_nddo*** may be set .FALSE. on input
            !    and we will preserve this if there is no outer_scf
         ELSEIF ( scf_control%outer_scf%type==outer_scf_scp ) THEN
! IF outer scf count is even, then only optimize SCP, else optimize KS
           IF ( scf_env%outer_scf%iter_count == 0 ) THEN
             scf_env%qs_ot_env(1)%settings%ks = .TRUE.
             scf_env%qs_ot_env(1)%settings%scp_nddo = .FALSE.
           ELSEIF ( MOD ( scf_env%outer_scf%iter_count, 2 ) /= 0 ) THEN
             scf_env%qs_ot_env(1)%settings%ks = .FALSE.
             scf_env%qs_ot_env(1)%settings%scp_nddo = .TRUE.
           ELSE
             scf_env%qs_ot_env(1)%settings%ks = .TRUE.
             scf_env%qs_ot_env(1)%settings%scp_nddo = .FALSE.
           END IF
         ENDIF
       END IF


       ! Some flags needed to be set at the beginning of the loop

       diis_step = .FALSE.
       use_jacobi = .FALSE.
       energy_only = .FALSE.


       ! SCF loop, optimisation of the wfn coefficients
       ! qs_env%rho%rho_r and qs_env%rho%rho_g should be up to date here

       scf_env%iter_count = 0
       exit_inner_loop    = .FALSE.
       scf_loop: DO

          CALL timeset(routineN//"_inner_loop",handle2)

          scf_env%iter_count = scf_env%iter_count + 1
          iter_count =  iter_count + 1
          CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=iter_count,error=error)

          IF (output_unit > 0) CALL m_flush(output_unit)

          total_steps = total_steps + 1
          ! ***SCP
          IF ( scp_dft ) CALL qs_ks_scp_update ( qs_env, just_energy=energy_only, error=error )
          ! ***SCP
          CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,calculate_forces=.FALSE.,&
               just_energy=energy_only,error=error)

          ! print 'heavy weight' or relatively expensive quantities
          CALL qs_scf_loop_print(qs_env,scf_env,para_env,error)

          scf_env%iter_param = 0.0_dp

          ! this takes known energy and derivatives and produces 
          ! new wfns and or density matrix
          SELECT CASE (scf_env%method)
          CASE DEFAULT
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"unknown scf method: "//&
                  cp_to_string(scf_env%method),error,failure)
          CASE(general_diag_method_nr)
             IF (dft_control%roks) THEN
                CALL qs_scf_loop_do_roks_diag(scf_env,mos,matrix_ks,matrix_s,&
                                              scf_control,scf_section,diis_step,&
                                              has_unit_metric,error)
             ELSE
                CALL qs_scf_loop_do_general_diag(scf_env,mos,matrix_ks,&
                                                 matrix_s,scf_control,scf_section, &
                                                 diis_step,use_jacobi,error)
             END IF
          CASE(special_diag_method_nr)
             IF (dft_control%roks) THEN
                CALL qs_scf_loop_do_roks_diag(scf_env,mos,matrix_ks,matrix_s,&
                                              scf_control,scf_section,diis_step,&
                                              has_unit_metric,error)
             ELSE
                CALL qs_scf_loop_do_special_diag(scf_env,mos,matrix_ks,&
                                                 scf_control,scf_section, &
                                                 diis_step,error)
             END IF
          CASE(ot_diag_method_nr)
             CALL cp_assert(.NOT.(dft_control%roks),cp_failure_level,cp_assertion_failed,&
                  routineP,"ROKS with OT diagonalization not possible",error,failure)
             CALL qs_scf_loop_do_ot_diag(scf_env,mos,matrix_ks,matrix_s,&
                                         scf_control,scf_section,diis_step,error)
          CASE(ot_method_nr)
             ! ***SCP
             IF ( scp_dft ) THEN
                CALL qs_scf_loop_do_ot(scf_env,scf_control%smear,mos,rho,qs_env%mo_derivs,energy%total, &
                     matrix_s, aux_coeff_set=aux_coeff_set, energy_only=energy_only, &
                     has_unit_metric=has_unit_metric,error=error)
             ELSEIF ( scp_nddo ) THEN
                CALL qs_scf_loop_do_ot(scf_env,scf_control%smear,mos,rho,qs_env%mo_derivs,energy%total, &
                     matrix_s, pscp=pscp, fscp=ks_scp, energy_only=energy_only, &
                     has_unit_metric=has_unit_metric,error=error)
                CALL set_scp_env ( scp_env = scp_env, pscp = pscp, error = error )
             ELSE
                CALL qs_scf_loop_do_ot(scf_env,scf_control%smear,mos,rho,qs_env%mo_derivs,energy%total, &
                     matrix_s, energy_only=energy_only,has_unit_metric=has_unit_metric,error=error)
             ENDIF
             ! ***SCP
          END SELECT

          ! another heavy weight print object, print controlled by dft_section
          IF (SIZE(mos) > 1) THEN
             CALL write_mo_set(mos(1)%mo_set,atomic_kind_set,particle_set,4,&
                               dft_section,spin="ALPHA",last=.FALSE.,error=error)
             CALL write_mo_set(mos(2)%mo_set,atomic_kind_set,particle_set,4,&
                               dft_section,spin="BETA",last=.FALSE.,error=error)
          ELSE
             CALL write_mo_set(mos(1)%mo_set,atomic_kind_set,particle_set,4,&
                               dft_section,last=.FALSE.,error=error)
          END IF

          ! ** calculation of the harris energy correction ***
          IF (harris_flag) THEN
             CALL harris_energy_correction(qs_env, harris_env, para_env=para_env, fast=.TRUE., error=error)
             IF (scf_env%method .NE. ot_method_nr) THEN
                CALL harris_eigenvalue_trace_KS_Pmix(scf_env, qs_env, harris_env,error=error)
             ELSE
                CALL harris_eigenvalue_calculation(qs_env=qs_env, harris_env=harris_env,error=error)
             END IF

             CALL section_vals_val_get(harris_section, "HARRIS_ENERGY",l_val=harris_energy,error=error)
             IF ((globenv%run_type_id == debug_run) .OR. (harris_energy)) THEN
                energy%total = harris_env%harris_energy%Eharris
             END IF
          END IF

          SELECT CASE(scf_env%mixing_method)
          CASE(direct_mixing_nr)
             CALL scf_env_density_mixing(scf_env%p_mix_new,&
                  scf_env, scf_env%iter_delta, qs_env=qs_env, &
                  diis=diis_step, error=error)
          CASE(gspace_mixing_nr,pulay_mixing_nr,broyden_mixing_nr,broyden_mixing_new_nr)
!             !matrix formulation of gspace mixing using the Kerker damping
             ! Compute the difference p_out-p_in
             CALL self_consistency_check(qs_env,scf_env,scf_env%p_mix_new,&
                  delta= scf_env%iter_delta, error=error)
          CASE(no_mixing_nr)
          CASE DEFAULT
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"unknown scf mixing method: "//&
                  cp_to_string(scf_env%mixing_method),error,failure)
          END SELECT

          t2 = m_walltime()
          IF ((output_unit>0).AND.scf_env%print_iter_line) THEN
             WRITE (UNIT=output_unit,&
                  FMT="(T2,I5,2X,A,T22,E10.2,T33,F10.3,T43,F17.10,F20.10)")&
                  scf_env%iter_count,TRIM(scf_env%iter_method),&
                  scf_env%iter_param,t2 - t1,scf_env%iter_delta,energy%total
          END IF

          IF (harris_flag) THEN
             CALL harris_postprocessing(harris_env, error=error)
          END IF

          ! ** can we exit the SCF loop  ?
          CALL external_control(should_stop,"SCF",globenv,error)
          IF (scf_env%iter_delta < scf_control%eps_scf) THEN
             IF (output_unit>0) THEN
                WRITE(UNIT=output_unit,FMT="(/,T3,A,I5,A/)")&
                     "*** SCF run converged in ",scf_env%iter_count," steps ***"
             END IF
             exit_inner_loop = .TRUE.
          ELSE IF (should_stop.OR.  scf_env%iter_count >= scf_control%max_scf) THEN
             IF (output_unit>0) THEN
                WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
                     "*** SCF run NOT converged ***"
             END IF
             converged=.FALSE.
             exit_inner_loop = .TRUE.
          END IF
          ! ** In case we decide to exit we perform few more check to see if this one
          ! ** is really the last SCF step
          IF (exit_inner_loop) THEN
! *********  Need to reset qs_ot_env%settings%scp_dft = .TRUE. to ensure *********
!            consistent deallocation in cleanup_scf_loop
             IF ( scp_dft  .AND. scf_control%outer_scf%type==outer_scf_scp ) &
             scf_env%qs_ot_env(1)%settings%scp_dft=.TRUE.
! *********  Need to reset qs_ot_env%settings%scp_dft = .TRUE. to ensure *********
!            consistent deallocation in cleanup_scf_loop
             IF ( scp_nddo  .AND. scf_control%outer_scf%type==outer_scf_scp ) &
             scf_env%qs_ot_env(1)%settings%scp_nddo=.TRUE.

             CALL cleanup_scf_loop(scf_env,error)

             ! now, print out energies and charges corresponding to the obtained wfn
             ! (this actually is not 100% consistent at this point)!
             IF ( scp_dft ) CALL qs_scf_scp_print_summary ( output_unit, scp_env, qs_env%qmmm, error )
             CALL qs_scf_print_summary(output_unit,rho,qs_charges,energy,scf_env, &
                  dft_control,qs_env%qmmm,gapw,gapw_xc)

             IF (harris_flag) THEN
                energy%total = harris_env%harris_energy%Eharris
             END IF

             ! *** mixing methods need to undo mixing of the density matrix 
             !     (restore original density) ***
             IF (scf_env%mixing_method>0) THEN
                SELECT CASE(scf_env%mixing_method)
                CASE(direct_mixing_nr)
                  CALL scf_env_density_mixing(scf_env%p_mix_new,&
                                            scf_env,scf_env%iter_delta,&
                                            qs_env=qs_env, diis=diis_step,&
                                            invert=.TRUE.,error=error)
                  DO ispin=1,dft_control%nspins
                     CALL transfer_matrix(scf_env%p_mix_new(ispin)%matrix,&
                                        rho%rho_ao(ispin)%matrix,error=error)
                  END DO
                CASE(gspace_mixing_nr)
                   !
                CASE(pulay_mixing_nr)
                   !
                CASE(broyden_mixing_nr)
                CASE(broyden_mixing_new_nr)
                END SELECT
             ENDIF

             !   *** update rspace rho since the mo changed
             !   *** this might not always be needed (i.e. no post calculation / no forces )
             !   *** but guarantees that rho and wfn are consistent at this point
             CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
             CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)
             ! ***SCP
             IF ( scp_dft ) THEN
                CALL update_rhoscp ( qs_env = qs_env, error = error )
                CALL qs_ks_scp_did_change(qs_env, rho_changed=.TRUE.,error=error)
             ENDIF
             ! ***SCP

             outer_loop_converged=.TRUE.
             IF (scf_control%outer_scf%have_scf) THEN
                ! We have an outer SCF loop...
                scf_env%outer_scf%iter_count=scf_env%outer_scf%iter_count+1
                outer_loop_converged=.FALSE.

                CALL outer_loop_gradient(qs_env,scf_env,error)
                outer_loop_eps=SQRT(SUM(scf_env%outer_scf%gradient(:,scf_env%outer_scf%iter_count)**2))/ &
                     SIZE(scf_env%outer_scf%gradient,1)
                IF (outer_loop_eps<scf_control%outer_scf%eps_scf) outer_loop_converged=.TRUE.
             END IF
             ! ** Let's tag the last SCF cycle so we can print informations only of the last step
             IF (outer_loop_converged.OR.&
             scf_env%outer_scf%iter_count>scf_control%outer_scf%max_scf .OR.&
             should_stop) THEN
                CALL cp_iterate(logger%iter_info,last=.TRUE.,iter_nr=iter_count,error=error)
             END IF
          END IF

          !   *** Write WaveFunction restart file ***
          CALL write_mo_set(mos,particle_set,dft_section=dft_section,scp=scp_dft,&
                            scp_env=scp_env,&
                            atomic_kind_set=atomic_kind_set,error=error)

          !   *** Exit if we have finished with the SCF inner loop ***
          IF (exit_inner_loop) THEN
             CALL timestop(handle2)
             EXIT scf_loop
          END IF

          IF (.NOT.BTEST(cp_print_key_should_output(logger%iter_info,&
               scf_section,"PRINT%ITERATION_INFO/TIME_CUMUL",error=error),cp_p_file)) &
               t1 = m_walltime()

          ! ** mixing methods have the new density matrix in p_mix_new **
          IF (scf_env%mixing_method>0) THEN
              DO ispin=1,dft_control%nspins
                 CALL transfer_matrix(scf_env%p_mix_new(ispin)%matrix,&
                                          rho%rho_ao(ispin)%matrix,error=error)
              END DO
          ENDIF

          ! ** update qs_env%rho
          CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
          ! ** Density mixing through density matrix or on the reciprocal space grid (exclusive)
          IF(scf_env%mixing_method==gspace_mixing_nr) THEN
            CALL gmix_potential_only(scf_env, rho, error)
          ELSEIF(scf_env%mixing_method==pulay_mixing_nr) THEN
            CALL pulay_mixing(scf_env, rho, qs_env%para_env, error)
          ELSEIF(scf_env%mixing_method==broyden_mixing_nr) THEN
            CALL broyden_mixing(scf_env, rho, qs_env%para_env, error)
          ELSEIF(scf_env%mixing_method==broyden_mixing_new_nr) THEN
            CALL broyden_mixing_new(scf_env, rho, qs_env%para_env, delta_rhog, delta_rhog_p, error)
          END IF
          CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)

          ! *** SCP
          IF ( scp_dft ) THEN
             CALL update_rhoscp ( qs_env = qs_env, error = error )
             CALL qs_ks_scp_did_change(qs_env, rho_changed=.TRUE.,error=error)
          ENDIF
          ! *** SCP
          CALL timestop(handle2)

       END DO scf_loop

       IF (.NOT. scf_control%outer_scf%have_scf) EXIT scf_outer_loop

       ! ** In case we use the OUTER SCF loop let's print some info..
       IF (output_unit>0) WRITE(output_unit,'(/,T3,A,I4,A,E10.2,A,F21.10)') &
                            "outer SCF iter = ",scf_env%outer_scf%iter_count, &
                            " RMS gradient = ",outer_loop_eps," energy =",energy%total

       IF (outer_loop_converged) THEN
            IF (output_unit>0) WRITE(output_unit,'(T3,A,I4,A,I4,A,/)') &
                  "outer SCF loop converged in",scf_env%outer_scf%iter_count,&
                  " iterations or ",total_steps," steps"
            EXIT scf_outer_loop
       ENDIF

       IF (scf_env%outer_scf%iter_count>scf_control%outer_scf%max_scf &
           .OR. should_stop ) THEN
         IF (output_unit>0) WRITE(output_unit,'(T3,A,I4,A,I4,A,/)') &
           "outer SCF loop FAILED to converge after ", &
           scf_env%outer_scf%iter_count," iterations or ",total_steps," steps"
         EXIT scf_outer_loop
       ENDIF

       CALL outer_loop_optimize(scf_env,scf_control,error)
       CALL outer_loop_update_qs_env(qs_env,scf_env,error)
       CALL qs_ks_did_change(ks_env,potential_changed=.TRUE.,error=error)
       ! *** SCP
       IF ( scp_dft ) CALL qs_ks_scp_did_change(qs_env, potential_changed=.TRUE.,error=error)
       ! *** SCP

    END DO scf_outer_loop

    IF(qs_env%dft_control%qs_control%becke_restraint)THEN
       CALL pw_env_get(qs_env%pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,&
            qs_env%dft_control%qs_control%becke_control%becke_pot%pw,error=error)
       qs_env%dft_control%qs_control%becke_control%need_pot=.TRUE.
    END IF
    CALL cp_rm_iter_level(logger%iter_info,level_name="QS_SCF",error=error)
    CALL timestop(handle)

  END SUBROUTINE scf_env_do_scf

! *****************************************************************************
!> \brief the inner loop of scf, specific to diagonalization with S matrix
!>       basically, in goes the ks matrix out goes a new p matrix
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE qs_scf_loop_do_general_diag(scf_env,mos,matrix_ks,&
                                         matrix_s,scf_control,scf_section,&
                                         diis_step,use_jacobi,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step, use_jacobi
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ispin
    LOGICAL                                  :: do_level_shift
    REAL(kind=dp)                            :: diis_error

    DO ispin=1,SIZE(matrix_ks)
      CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,&
                         scf_env%scf_work1(ispin)%matrix,&
                         error=error)
    END DO

    IF (scf_env%iter_count > 1 .AND. .NOT. scf_env%skip_diis) THEN
      CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_env%scf_work1,&
                          scf_env%scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                          scf_control%eps_diis,scf_control%nmixing,&
                          s_matrix=matrix_s,&
                          scf_section=scf_section,error=error)
    ELSE
      diis_step = .FALSE.
    END IF

    do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
                      ((scf_control%density_guess == core_guess).OR.&
                       (scf_env%iter_count > 1)))

    IF ((scf_env%iter_count > 1).AND.(scf_env%iter_delta < scf_control%diagonalization%eps_jacobi)) THEN
      use_jacobi = .TRUE.
    ELSE
      use_jacobi = .FALSE.
    END IF

    IF (diis_step) THEN
      scf_env%iter_param = diis_error
      IF (use_jacobi) THEN
        scf_env%iter_method = "DIIS/Jacobi"
      ELSE
        scf_env%iter_method = "DIIS/Diag."
      END IF
    ELSE
      SELECT CASE(scf_env%mixing_method)
      CASE(direct_mixing_nr)
        scf_env%iter_param = scf_env%p_mix_alpha
        IF (use_jacobi) THEN
          scf_env%iter_method = "P_Mix/Jacobi"
        ELSE
          scf_env%iter_method = "P_Mix/Diag."
        END IF
      CASE(gspace_mixing_nr)
        scf_env%iter_param = scf_env%p_mix_alpha
        IF (use_jacobi) THEN
          scf_env%iter_method = "Kerker/Jacobi"
        ELSE
          scf_env%iter_method = "Kerker/Diag."
        END IF
      CASE(pulay_mixing_nr)
        scf_env%iter_param = scf_env%pulay_alpha
        IF (use_jacobi) THEN
          scf_env%iter_method = "Pulay/Jacobi"
        ELSE
          scf_env%iter_method = "Pulay/Diag."
        END IF
      CASE(broyden_mixing_nr)
        scf_env%iter_param = scf_env%p_mix_alpha
        IF (use_jacobi) THEN
          scf_env%iter_method = "Broy/Jacobi"
        ELSE
          scf_env%iter_method = "Broy/Diag."
        END IF
      END SELECT
    END IF

    scf_env%iter_delta = 0.0_dp

    DO ispin=1,SIZE(matrix_ks)
      CALL eigensolver(matrix_ks=scf_env%scf_work1(ispin)%matrix,&
                       mo_set=mos(ispin)%mo_set,&
                       ortho=scf_env%ortho,&
                       work=scf_env%scf_work2,&
                       do_level_shift=do_level_shift,&
                       level_shift=scf_control%level_shift,&
                       cholesky_method=scf_env%cholesky_method,&
                       use_jacobi=use_jacobi,&
                       jacobi_threshold=scf_control%diagonalization%jacobi_threshold,&
                       error=error)
      CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,&
                             smear=scf_control%smear,&
                             error=error)
      CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                    scf_env%p_mix_new(ispin)%matrix,&
                                    error=error)
    END DO ! ispin

  END SUBROUTINE qs_scf_loop_do_general_diag

! *****************************************************************************
!> \brief the inner loop of scf, specific to diagonalization without S matrix
!>       basically, in goes the ks matrix out goes a new p matrix
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE qs_scf_loop_do_special_diag(scf_env,mos,matrix_ks,scf_control,&
                                         scf_section,diis_step,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ispin
    LOGICAL                                  :: do_level_shift
    REAL(kind=dp)                            :: diis_error

    DO ispin=1,SIZE(matrix_ks)
      CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,scf_env%scf_work1(ispin)%matrix,error=error)
    END DO
    IF (scf_env%iter_count > 1 .AND. .NOT. scf_env%skip_diis) THEN
      CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_env%scf_work1,&
                          scf_env%scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                          scf_control%eps_diis,scf_control%nmixing,&
                          scf_section=scf_section,&
                          error=error)
    ELSE
      diis_step = .FALSE.
    END IF

    do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
                      ((scf_control%density_guess == core_guess).OR.(scf_env%iter_count > 1)))
    IF (diis_step) THEN
      scf_env%iter_param = diis_error
      scf_env%iter_method = "DIIS/Diag."
    ELSE
      SELECT CASE(scf_env%mixing_method)
      CASE(direct_mixing_nr)
        scf_env%iter_param = scf_env%p_mix_alpha
          scf_env%iter_method = "P_Mix/Diag."
      CASE(gspace_mixing_nr)
        scf_env%iter_param = scf_env%p_mix_alpha
          scf_env%iter_method = "Kerker/Diag."
      CASE(pulay_mixing_nr)
        scf_env%iter_param = scf_env%pulay_alpha
          scf_env%iter_method = "Pulay/Diag."
      CASE(broyden_mixing_nr)
        scf_env%iter_param = scf_env%p_mix_alpha
        scf_env%iter_method = "Broy/Diag."
      END SELECT
    END IF
    scf_env%iter_delta = 0.0_dp

    DO ispin=1,SIZE(matrix_ks)
      CALL simple_eigensolver(matrix_ks=scf_env%scf_work1(ispin)%matrix,&
                              mo_set=mos(ispin)%mo_set,&
                              work=scf_env%scf_work2,&
                              do_level_shift=do_level_shift,&
                              level_shift=scf_control%level_shift,&
                              error=error)
      CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,&
                             smear=scf_control%smear,&
                             error=error)
      CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                    scf_env%p_mix_new(ispin)%matrix,&
                                    error=error)

    END DO ! ispin

  END SUBROUTINE qs_scf_loop_do_special_diag

! *****************************************************************************
!> \brief the inner loop of scf, specific to iterative diagonalization using OT
!>        with S matrix; basically, in goes the ks matrix out goes a new p matrix
!> \par History
!>      10.2008 created [JGH]
! *****************************************************************************
  SUBROUTINE qs_scf_loop_do_ot_diag(scf_env,mos,matrix_ks,matrix_s,&
               scf_control,scf_section,diis_step,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: homo, ispin, nmo
    LOGICAL                                  :: do_level_shift
    REAL(kind=dp)                            :: diis_error, eps_iter
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    DO ispin=1,SIZE(matrix_ks)
      CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,scf_env%scf_work1(ispin)%matrix,error=error)
    END DO

    IF (scf_env%iter_count > 1 .AND. .NOT. scf_env%skip_diis) THEN
      CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_env%scf_work1,&
                          scf_env%scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                          scf_control%eps_diis,scf_control%nmixing,&
                          s_matrix=matrix_s,&
                          scf_section=scf_section,error=error)
    ELSE
      diis_step = .FALSE.
    END IF

    do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
                      ((scf_control%density_guess == core_guess).OR.&
                       (scf_env%iter_count > 1)))

    eps_iter = scf_control%diagonalization%eps_iter
    IF (diis_step) THEN
      scf_env%iter_param = diis_error
      scf_env%iter_method = "DIIS/OTdiag"
      DO ispin=1,SIZE(matrix_ks)
        CALL copy_fm_to_sm(scf_env%scf_work1(ispin)%matrix,matrix_ks(ispin)%matrix)
      END DO
      eps_iter = MAX(eps_iter,scf_control%diagonalization%eps_adapt*diis_error)
    ELSE
      SELECT CASE(scf_env%mixing_method)
      CASE(direct_mixing_nr)
        scf_env%iter_param = scf_env%p_mix_alpha
          scf_env%iter_method = "P_Mix/OTdiag."
      CASE(gspace_mixing_nr)
        scf_env%iter_param = scf_env%p_mix_alpha
          scf_env%iter_method = "Kerker/OTdiag."
      CASE(pulay_mixing_nr)
        scf_env%iter_param = scf_env%pulay_alpha
          scf_env%iter_method = "Pulay/OTdiag."
      CASE(broyden_mixing_nr)
        scf_env%iter_param = scf_env%p_mix_alpha
        scf_env%iter_method = "Broy/Diag."
      END SELECT
    END IF

    scf_env%iter_delta = 0.0_dp

    DO ispin=1,SIZE(matrix_ks)
       CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,nmo=nmo,homo=homo)
       CALL ot_eigensolver(matrix_h=matrix_ks(ispin)%matrix,&
                           matrix_s=matrix_s(1)%matrix,&
                           matrix_c=mo_coeff,&
                           preconditioner=scf_env%ot_preconditioner(1)%preconditioner,&
                           eps_gradient=eps_iter,&
                           iter_max=scf_control%diagonalization%max_iter,&
                           silent=.TRUE.,&
                           ot_settings=scf_control%diagonalization%ot_settings,&
                           error=error)

       CALL calculate_subspace_eigenvalues(mo_coeff,matrix_ks(ispin)%matrix,&
              para_env=mo_coeff%matrix_struct%para_env, do_rotation = .TRUE., error=error)

       CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,smear=scf_control%smear,error=error)
       CALL calculate_density_matrix(mos(ispin)%mo_set,scf_env%p_mix_new(ispin)%matrix,error=error)
    END DO

  END SUBROUTINE qs_scf_loop_do_ot_diag

! *****************************************************************************
!> \brief the inner loop of scf, specific to using to the orbital transformation method
!>       basically, in goes the ks matrix out goes a new p matrix
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE qs_scf_loop_do_ot(scf_env,smear,mos,rho,mo_derivs,total_energy,&
                               matrix_s,aux_coeff_set,pscp,fscp,energy_only, &
                               has_unit_metric,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(smear_type)                         :: smear
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    REAL(KIND=dp), INTENT(IN)                :: total_energy
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(aux_coeff_set_type), OPTIONAL, &
      POINTER                                :: aux_coeff_set
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: pscp, fscp
    LOGICAL, INTENT(INOUT)                   :: energy_only
    LOGICAL, INTENT(IN)                      :: has_unit_metric
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_loop_do_ot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    TYPE(real_matrix_type), POINTER          :: orthogonality_metric

    CALL timeset(routineN,handle)


    IF (has_unit_metric) THEN
      NULLIFY(orthogonality_metric)
    ELSE
      orthogonality_metric=>matrix_s(1)%matrix
    END IF

    ! in case of LSD the first spin qs_ot_env will drive the minimization
    ! in the case of a restricted calculation, it will make sure the spin orbitals are equal

    CALL ot_scf_mini(mos,mo_derivs,smear,orthogonality_metric, &
                    total_energy, aux_coeff_set,    &
                    pscp,fscp,energy_only,scf_env%iter_delta, & 
                    scf_env%qs_ot_env,error=error)
    DO ispin=1,SIZE(mos)
       CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,&
                              smear=smear,&
                              error=error)
    ENDDO

    DO ispin=1,SIZE(mos)
      CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                    rho%rho_ao(ispin)%matrix,&
                                    error=error)
    END DO


    scf_env%iter_method=scf_env%qs_ot_env(1)%OT_METHOD_FULL
    scf_env%iter_param=scf_env%qs_ot_env(1)%ds_min

    CALL timestop(handle)

  END SUBROUTINE qs_scf_loop_do_ot

! *****************************************************************************
!> \brief writes rather detailed summary of SCP info
!>      after the SCF
!> \par History
!>      03.2008 created [CJM]
! *****************************************************************************
  SUBROUTINE qs_scf_scp_print_summary(output_unit,scp_env,qmmm,error)
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(scp_environment_type), POINTER      :: scp_env
    LOGICAL, INTENT(IN)                      :: qmmm
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(scp_energy_type), POINTER           :: energy

    CALL get_scp_env ( scp_env=scp_env,&
         energy=energy, error=error)
    IF (output_unit>0) THEN
       IF ( qmmm )  WRITE (UNIT=output_unit,FMT="((T3,A,T55,F25.14))")&
            "SCP Hartree (SCP density, QMMM potential)       ",energy % e_scp_qmmm
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "SCP Hartree (SCP,SCP)                          ",energy % e_scp_hartree,&
            "SCP Hartree Self (SCP,SCP)                     ",energy%e_scp_self,&
            "SCP Hartree (KS,SCP)                           ",energy % e_scp_ks,&
            "SCP Hartree Self (KS,SCP)                      ",energy%e_scp_ks_self,&
            "SCP Hartree Self Core (KS core,SCP)            ",energy%e_scp_core,&
            "SCP Polarization Kernel                        ",energy % e_scp_kernel, &
            "SCP TOTAL                                      ",energy % e_scp_total

       CALL m_flush(output_unit)
    END IF
  END SUBROUTINE qs_scf_scp_print_summary

! *****************************************************************************
!> \brief writes rather detailed summary of densities and energies
!>      after the SCF
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
SUBROUTINE qs_scf_print_summary(output_unit,rho,qs_charges,energy,&
                                scf_env,dft_control,qmmm,gapw,gapw_xc)

    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(dft_control_type), POINTER          :: dft_control
    LOGICAL, INTENT(IN)                      :: qmmm, gapw, gapw_xc

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_print_summary', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    REAL(kind=dp)                            :: tot1_h, tot1_s

    CALL timeset(routineN,handle)

    IF (output_unit>0) THEN
       IF(.NOT.(dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb)) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
            "Total electronic density (r-space): ",&
            accurate_sum(rho%tot_rho_r),&
            accurate_sum(rho%tot_rho_r)+ REAL(scf_env%nelectron,dp),&
            "Total core charge density (r-space):",&
            qs_charges%total_rho_core_rspace,&
            qs_charges%total_rho_core_rspace - REAL(scf_env%nelectron+dft_control%charge,dp)
         IF(gapw) THEN
            tot1_h =  qs_charges%total_rho1_hard(1)
            tot1_s =  qs_charges%total_rho1_soft(1)
            DO ispin=2,dft_control%nspins
              tot1_h = tot1_h + qs_charges%total_rho1_hard(ispin)
              tot1_s = tot1_s + qs_charges%total_rho1_soft(ispin)
            END DO
            WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
                  "Hard and soft densities (Lebedev):",&
                  tot1_h, tot1_s
            WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                 "Total Rho_soft + Rho1_hard - Rho1_soft (r-space): ",&
                 accurate_sum(rho%tot_rho_r)+ tot1_h - tot1_s ,&
                 "Total charge density (r-space):      ",&
                 accurate_sum(rho%tot_rho_r)+ tot1_h - tot1_s &
                 + qs_charges%total_rho_core_rspace,&
                 "Total Rho_soft + Rho0_soft (g-space):",&
                 qs_charges%total_rho_gspace
         ELSE
           WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
              "Total charge density (r-space):     ",&
              accurate_sum(rho%tot_rho_r)+&
              qs_charges%total_rho_core_rspace,&
              "Total charge density (g-space):     ",qs_charges%total_rho_gspace
         END IF
       END IF
       IF (dft_control%qs_control%semi_empirical) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Core-core repulsion energy [eV]:               ",energy%core_overlap*evolt,&
            "Core Hamiltonian energy [eV]:                  ",energy%core*evolt,&
            "Two-electron integral energy [eV]:             ",energy%hartree*evolt,&
            "Electronic energy [eV]:                        ",&
                                              (energy%core+0.5_dp*energy%hartree)*evolt
       ELSEIF (dft_control%qs_control%dftb) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Core Hamiltonian energy:                       ",energy%core,&
            "Repulsive potential energy:                    ",energy%repulsive,&
            "Electronic energy:                             ",energy%hartree,&
            "Dispersion energy:                             ",energy%dispersion
       ELSE
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Overlap energy of the core charge distribution:",energy%core_overlap,&
            "Self energy of the core charge distribution:   ",energy%core_self,&
            "Core Hamiltonian energy:                       ",energy%core,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc
         IF (energy%e_hartree /= 0.0_dp) &
              WRITE (UNIT=output_unit,FMT="(T3,A,/,T3,A,T55,F25.14)")&
              "Coulomb Electron-Electron Interaction Energy ",&
              "- Already included in the total Hartree term ",energy%e_hartree
         IF (energy%ex /= 0.0_dp)&
               WRITE (UNIT=output_unit,FMT="(T3,A,T55,F25.14)")&
               "Hartree-Fock Exchange energy:                  ",energy%ex
         IF (energy%dispersion /= 0.0_dp)&
               WRITE (UNIT=output_unit,FMT="(T3,A,T55,F25.14)")&
               "Dispersion energy:                             ",energy%dispersion
         IF(gapw) THEN
           WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
              "GAPW| Exc from hard and soft atomic rho1:      ",energy%exc1,&
              "GAPW| local Eh = 1 center integrals:           ",energy%hartree_1c,&
              "GAPW| local Eh = 2 center integrals:           ",energy%hartree_2c,&
              "GAPW| local Eh = 3 centers integrals:          ",energy%hartree_3c
         END IF
         IF(gapw_xc) THEN
           WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
              "GAPW_XC| Exc from hard and soft atomic rho1:      ",energy%exc1
         END IF
       END IF
       IF (dft_control%smear) THEN
         WRITE (UNIT=output_unit,FMT="((T3,A,T55,F25.14))")&
          "Electronic entropic energy:",energy%kTS
       END IF
       IF (dft_control%dft_plus_u) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
          "DFT+U energy:",energy%dft_plus_u
       END IF
       IF  (qmmm) THEN
          WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "QM/MM Electrostatic energy:                    ",energy%qmmm_el
       END IF
       IF (dft_control%qs_control%mulliken_restraint) THEN
          WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                "Mulliken restraint energy: ",energy%mulliken
       END IF
       IF (dft_control%qs_control%semi_empirical) THEN
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Total energy [eV]:                             ",energy%total*evolt
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Atomic reference energy [eV]:                  ",energy%core_self*evolt,&
            "Heat of formation [kcal/mol]:                  ",&
            (energy%total+energy%core_self)*kcalmol
       ELSE
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Total energy:                                  ",energy%total
       END IF
       CALL m_flush(output_unit)
    END IF

    CALL timestop(handle)

END SUBROUTINE qs_scf_print_summary

! *****************************************************************************
!> \brief collects the 'heavy duty' printing tasks out of the SCF loop
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE qs_scf_loop_print(qs_env,scf_env,para_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_loop_print', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, iw
    REAL(KIND=dp)                            :: mo_mag_max, mo_mag_min, &
                                                orthonormality
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_p, matrix_s
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                scf_section

    logger => cp_error_get_logger(error)
    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,input=input, dft_control=dft_control, &
                    error=error)

    dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
    scf_section => section_vals_get_subs_vals(dft_section,"SCF",error=error)

    CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, &
                    matrix_ks=matrix_ks,rho=rho,mos=mos, &
                    error=error)

    matrix_p => rho%rho_ao

    DO ispin=1,dft_control%nspins

      IF (BTEST(cp_print_key_should_output(logger%iter_info,&
           dft_section,"PRINT%AO_MATRICES/DENSITY",error=error),cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%AO_MATRICES/DENSITY",&
              extension=".Log",error=error)
         CALL write_sparse_matrix(matrix_p(ispin)%matrix,4,6,qs_env,para_env,&
              output_unit=iw,error=error)
         CALL cp_print_key_finished_output(iw,logger,dft_section,&
              "PRINT%AO_MATRICES/DENSITY", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info,&
           dft_section,"PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",error=error),cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",&
              extension=".Log",error=error)
         IF (dft_control%qs_control%semi_empirical) THEN
            CALL write_sparse_matrix(matrix_ks(ispin)%matrix,4,6,qs_env,para_env,&
                 scale=evolt,output_unit=iw,error=error)
         ELSE
            CALL write_sparse_matrix(matrix_ks(ispin)%matrix,4,6,qs_env,para_env,&
                 output_unit=iw,error=error)
         END IF
         CALL cp_print_key_finished_output(iw,logger,dft_section,&
              "PRINT%AO_MATRICES/KOHN_SHAM_MATRIX", error=error)
      END IF

    ENDDO

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         scf_section,"PRINT%MO_ORTHONORMALITY",error=error),cp_p_file)) THEN
       IF(scf_env%method==special_diag_method_nr) THEN
         CALL calculate_orthonormality(orthonormality,mos,error=error)
       ELSE
         CALL calculate_orthonormality(orthonormality,mos,matrix_s(1)%matrix,error=error)
       END IF
       iw=cp_print_key_unit_nr(logger,scf_section,"PRINT%MO_ORTHONORMALITY",&
            extension=".scfLog",error=error)
       IF (iw>0) THEN
          WRITE(iw,'(T8,A,T60,E20.4)')  &
                 " Maximum deviation from MO S-orthonormality",orthonormality
       ENDIF
       CALL cp_print_key_finished_output(iw,logger,scf_section,&
            "PRINT%MO_ORTHONORMALITY", error=error)
    ENDIF

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         scf_section,"PRINT%MO_MAGNITUDE",error=error),cp_p_file)) THEN
       CALL calculate_magnitude(mos,mo_mag_min,mo_mag_max,error=error)
       iw=cp_print_key_unit_nr(logger,scf_section,"PRINT%MO_MAGNITUDE",&
            extension=".scfLog",error=error)
       IF (iw>0) THEN
          WRITE(iw,'(T8,A,T40,2E20.4)')  &
                 " Minimum/Maximum MO magnitude ",mo_mag_min,mo_mag_max
       ENDIF
       CALL cp_print_key_finished_output(iw,logger,scf_section,&
            "PRINT%MO_MAGNITUDE", error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_loop_print

! *****************************************************************************
!> \brief inits those objects needed if you want to restart the scf with, say
!>        only a new initial guess, or different density functional or ...
!>        this will happen just before the scf loop starts
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

    CHARACTER(LEN=80)                        :: msg
    INTEGER                                  :: handle, i, ispin, nmo, &
                                                number_of_OT_envs, stat
    LOGICAL                                  :: dftb, do_rotation, failure, &
                                                has_unit_metric, scp_dft, &
                                                scp_nddo, semi
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(real_matrix_type), POINTER          :: ks_scp, orthogonality_metric, &
                                                pscp
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(scp_environment_type), POINTER      :: scp_env

    CALL timeset(routineN,handle)

    NULLIFY(scf_control,matrix_s,matrix_ks,dft_control,mos,mo_coeff)
! **SCP
    NULLIFY ( scp_env, aux_coeff_set, pscp, ks_scp )
! **SCP

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    para_env=>qs_env%para_env

    CALL get_qs_env(qs_env=qs_env,&
         scf_control=scf_control,&
         dft_control=dft_control,&
         mos=mos,matrix_ks=matrix_ks,&
         matrix_s=matrix_s, error=error)

    scp_dft = dft_control%scp
    scp_nddo = dft_control%qs_control%se_control%scp 

    ! this just guarantees that all mo_occupations match the eigenvalues, if smear
    ! TO BE FIXED, right now, can't be done with xas
    IF (.NOT. (dft_control%do_xas_calculation .AND. dft_control%xas_estate>0)) THEN
       DO ispin=1,dft_control%nspins
          CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,&
                                 smear=scf_control%smear,&
                                 error=error)
       ENDDO
    ENDIF

    SELECT CASE (scf_env%method)
    CASE DEFAULT

      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method:"//&
               cp_to_string(scf_env%method),error,failure)

    CASE (general_diag_method_nr,special_diag_method_nr)
       IF(.NOT. scf_env%skip_diis) THEN
         IF (.NOT.ASSOCIATED(scf_env%scf_diis_buffer)) THEN
            CALL qs_diis_b_create(scf_env%scf_diis_buffer,nbuffer=scf_control%max_diis,error=error)
         END IF
         CALL qs_diis_b_clear(scf_env%scf_diis_buffer,error=error)
       END IF

    CASE (ot_diag_method_nr)
  
       IF(.NOT. scf_env%skip_diis) THEN
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
       IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
         SELECT CASE(scf_control%diagonalization%ot_settings%preconditioner_type)
           ! these depend on the ks matrix
           CASE(ot_precond_full_all,ot_precond_full_single, ot_precond_full_single_inverse) 
             DO ispin=1,SIZE(scf_env%ot_preconditioner)
                CALL destroy_preconditioner(scf_env%ot_preconditioner(ispin)%preconditioner,error=error)
                DEALLOCATE(scf_env%ot_preconditioner(ispin)%preconditioner)
             ENDDO
             DEALLOCATE(scf_env%ot_preconditioner,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             CASE(ot_precond_none,ot_precond_full_kinetic,ot_precond_s_inverse, &
                  ot_precond_sparse_diag,ot_precond_sparse_kinetic,&
                  ot_precond_sparse_firstorder) ! these are 'independent'
           ! do nothing
           CASE DEFAULT
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
         END SELECT
       END IF

       ! add an OT preconditioner if none is present
       IF (.NOT.ASSOCIATED(scf_env%ot_preconditioner)) THEN
         SELECT CASE(scf_control%diagonalization%ot_settings%preconditioner_type)
           CASE(ot_precond_full_all,ot_precond_full_single_inverse)
             ALLOCATE(scf_env%ot_preconditioner(dft_control%nspins), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           CASE DEFAULT
             ALLOCATE(scf_env%ot_preconditioner(1), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         END SELECT
         DO i=1,SIZE(scf_env%ot_preconditioner)
           ALLOCATE(scf_env%ot_preconditioner(i)%preconditioner)
           CALL init_preconditioner(scf_env%ot_preconditioner(i)%preconditioner,&
                  para_env=qs_env%para_env,blacs_env=qs_env%blacs_env,error=error)
         ENDDO
       END IF

       SELECT CASE(scf_control%diagonalization%ot_settings%preconditioner_type)
         CASE(ot_precond_none)
           !
         CASE(ot_precond_full_all,ot_precond_full_single_inverse)
           DO ispin=1,dft_control%nspins
             CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
             CALL calculate_subspace_eigenvalues(mo_coeff,matrix_ks(ispin)%matrix,&
                        para_env=mo_coeff%matrix_struct%para_env, &
                        do_rotation = .TRUE., co_rotate=qs_env%mo_derivs(ispin)%matrix,error=error)
             CALL make_preconditioner(scf_env%ot_preconditioner(ispin)%preconditioner, &
                        scf_control%diagonalization%ot_settings%preconditioner_type, &
                        scf_control%diagonalization%ot_settings%precond_solver_type, &
                        matrix_ks(ispin)%matrix,matrix_s(1)%matrix,&
                        qs_env%kinetic(1)%matrix, &
                        mo_coeff,scf_control%diagonalization%ot_settings%energy_gap,error=error)
           ENDDO
         CASE DEFAULT
           ! constuct the preconditioner based on the first spin set.
           ! this might or might not be reasonable. Presumably it is.
           CALL get_mo_set(mo_set=mos(1)%mo_set, mo_coeff=mo_coeff)
           CALL make_preconditioner(scf_env%ot_preconditioner(1)%preconditioner, &
                     scf_control%diagonalization%ot_settings%preconditioner_type, &
                     scf_control%diagonalization%ot_settings%precond_solver_type, &
                     matrix_ks(1)%matrix,matrix_s(1)%matrix,&
                     qs_env%kinetic(1)%matrix, &
                     mo_coeff,scf_control%diagonalization%ot_settings%energy_gap,error=error)
       END SELECT

    CASE (ot_method_nr)


       dftb = .FALSE.
       semi = .FALSE.
       has_unit_metric=.FALSE.
       IF (dft_control%qs_control%semi_empirical) THEN
          IF (dft_control%qs_control%se_control%orthogonal_basis) THEN
                has_unit_metric=.TRUE.
          ENDIF
          semi = .TRUE.
       END IF
       IF (dft_control%qs_control%dftb) THEN
          IF (dft_control%qs_control%dftb_control%orthogonal_basis) THEN
                has_unit_metric=.TRUE.
          ENDIF
          dftb = .TRUE.
       END IF

       ! reortho the wavefunctions if we are having an outer scf and 
       ! this is not the first iteration
       ! this is useful to avoid the build-up of numerical noise
       ! however, we can not play this trick if restricted (don't mix non-equivalent orbs)
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


       IF (.NOT.ASSOCIATED(scf_env%qs_ot_env)) THEN

          ! restricted calculations require just one set of OT orbitals
          number_of_OT_envs=dft_control%nspins
          IF (dft_control%restricted) number_of_OT_envs=1

          ALLOCATE(scf_env%qs_ot_env(number_of_OT_envs),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          ! XXX Joost XXX should disentangle reading input from this part
          CALL ot_scf_read_input(scf_env%qs_ot_env,para_env,scf_section,error)

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
! **SCP
! Updates density and potential due to SCP
          IF ( scp_dft )  CALL qs_ks_scp_update ( qs_env, just_energy=.FALSE., error=error )
! **SCP
          CALL qs_ks_update_qs_env(qs_env%ks_env,&
                                   qs_env=qs_env,&
                                   calculate_forces=.FALSE.,&
                                   just_energy=.FALSE.,&
                                   error=error)

          ! if an old preconditioner is still around (i.e. outer SCF is active),
          ! remove it if this could be worthwhile
          IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
             SELECT CASE(scf_env%qs_ot_env(1)%settings%preconditioner_type)
             CASE(ot_precond_full_all,ot_precond_full_single, ot_precond_full_single_inverse) ! these depend on the ks matrix
               DO ispin=1,SIZE(scf_env%ot_preconditioner)
                  CALL destroy_preconditioner(scf_env%ot_preconditioner(ispin)%preconditioner,error=error)
                  DEALLOCATE(scf_env%ot_preconditioner(ispin)%preconditioner)
               ENDDO
               DEALLOCATE(scf_env%ot_preconditioner,stat=stat)
               CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             CASE(ot_precond_none,ot_precond_full_kinetic,ot_precond_s_inverse, &
                  ot_precond_sparse_diag,ot_precond_sparse_kinetic,&
                  ot_precond_sparse_firstorder) ! these are 'independent'
                ! do nothing
             CASE DEFAULT
              CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
          END IF

          ! add an OT preconditioner if none is present
          IF (.NOT.ASSOCIATED(scf_env%ot_preconditioner)) THEN
               SELECT CASE(scf_env%qs_ot_env(1)%settings%preconditioner_type)
               CASE(ot_precond_full_all,ot_precond_full_single_inverse)
                  ALLOCATE(scf_env%ot_preconditioner(dft_control%nspins), stat=stat)
                  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
               CASE DEFAULT
                  ALLOCATE(scf_env%ot_preconditioner(1), stat=stat)
                  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
               END SELECT
               DO i=1,SIZE(scf_env%ot_preconditioner)
                  ALLOCATE(scf_env%ot_preconditioner(i)%preconditioner)
                  CALL init_preconditioner(scf_env%ot_preconditioner(i)%preconditioner,&
                                           para_env=qs_env%para_env,&
                                           blacs_env=qs_env%blacs_env,error=error)
               ENDDO
          END IF

          !
          ! preconditioning still needs to be done correctly with has_unit_metric
          ! notice that a big part of the preconditioning (S^-1) is fine anyhow
          !
          IF (has_unit_metric) THEN
             NULLIFY(orthogonality_metric)
          ELSE
             orthogonality_metric=>matrix_s(1)%matrix
          ENDIF

          IF(semi .OR. dftb) THEN
             ! we have to treat the semi-empirical methods separately
             SELECT CASE(scf_env%qs_ot_env(1)%settings%preconditioner_type)
             CASE(ot_precond_none)
                CALL ot_scf_init(mo_array = mos, matrix_s = orthogonality_metric, & 
                     qs_ot_env = scf_env%qs_ot_env,error=error)
             CASE(ot_precond_full_all,ot_precond_full_single_inverse)
                ! this will rotate the MOs to be eigen states, which is not compatible with rotation
                ! e.g. mo_derivs here do not yet include potentially different occupations numbers
                do_rotation=scf_env%qs_ot_env(1)%settings%do_rotation
                CPPrecondition(.NOT.do_rotation,cp_failure_level,routineP,error,failure)
                DO ispin=1,dft_control%nspins
                   CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
                   CALL calculate_subspace_eigenvalues(mo_coeff,matrix_ks(ispin)%matrix,&
                        para_env=mo_coeff%matrix_struct%para_env, &
                        do_rotation = .TRUE., co_rotate=qs_env%mo_derivs(ispin)%matrix,error=error)
                   IF (has_unit_metric) THEN
                     CALL make_preconditioner(scf_env%ot_preconditioner(ispin)%preconditioner, &
                          scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                          scf_env%qs_ot_env(1)%settings%precond_solver_type, &
                          matrix_h=matrix_ks(ispin)%matrix, &
                          mo_coeff=mo_coeff, &
                          energy_gap=scf_env%qs_ot_env(1)%settings%energy_gap,error=error)
                   ELSE
                     CALL make_preconditioner(scf_env%ot_preconditioner(ispin)%preconditioner, &
                          scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                          scf_env%qs_ot_env(1)%settings%precond_solver_type, &
                          matrix_h=matrix_ks(ispin)%matrix, &
                          matrix_s=matrix_s(1)%matrix, &
                          matrix_t=matrix_s(1)%matrix, &
                          mo_coeff=mo_coeff,&
                          energy_gap=scf_env%qs_ot_env(1)%settings%energy_gap,error=error)
                   END IF
                ENDDO
                IF ( scp_nddo ) THEN
                   CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error = error )
                   CALL get_scp_env ( scp_env = scp_env, pscp=pscp, ks_scp=ks_scp, error = error )
                   CALL ot_scf_init(mos,orthogonality_metric, pscp=pscp, fscp=ks_scp, &
                        qs_ot_env=scf_env%qs_ot_env, error=error)
                ELSE
                   CALL ot_scf_init(mo_array = mos, matrix_s = orthogonality_metric, & 
                        qs_ot_env = scf_env%qs_ot_env,error=error)
                END IF
                DO ispin=1,SIZE(scf_env%qs_ot_env)
                   CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin),&
                        scf_env%ot_preconditioner(ispin)%preconditioner,error=error)
                ENDDO
             CASE(ot_precond_full_kinetic)
                msg="Full_kinetic not available for semi-empirical methods"
                CPErrorMessage(cp_failure_level,routineP,TRIM(msg),error)
                scf_env%qs_ot_env(1)%settings%preconditioner_type = ot_precond_none
                IF ( scp_nddo ) THEN
                   CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error = error )
                   CALL get_scp_env ( scp_env = scp_env, pscp=pscp, ks_scp=ks_scp, error = error )
                   CALL ot_scf_init(mos,orthogonality_metric, pscp=pscp, fscp=ks_scp, &
                        qs_ot_env=scf_env%qs_ot_env, error=error)
                ELSE
                   CALL ot_scf_init(mo_array = mos, matrix_s = orthogonality_metric, & 
                        qs_ot_env = scf_env%qs_ot_env,error=error)
                END IF
             CASE(ot_precond_s_inverse,ot_precond_sparse_diag,ot_precond_full_single)
                CALL get_mo_set(mo_set=mos(1)%mo_set, mo_coeff=mo_coeff)
                IF (has_unit_metric) THEN
                  CALL make_preconditioner(scf_env%ot_preconditioner(1)%preconditioner, &
                     scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                     scf_env%qs_ot_env(1)%settings%precond_solver_type, &
                     matrix_ks(1)%matrix,mo_coeff=mo_coeff,&
                     energy_gap=scf_env%qs_ot_env(1)%settings%energy_gap,error=error)
                ELSE
                  CALL make_preconditioner(scf_env%ot_preconditioner(1)%preconditioner, &
                     scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                     scf_env%qs_ot_env(1)%settings%precond_solver_type, &
                     matrix_ks(1)%matrix,matrix_s(1)%matrix,matrix_s(1)%matrix,&
                     mo_coeff,scf_env%qs_ot_env(1)%settings%energy_gap,error=error)
                END IF
                IF ( scp_nddo ) THEN
                   CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error = error )
                   CALL get_scp_env ( scp_env = scp_env, pscp=pscp, ks_scp=ks_scp, error = error )
                   CALL ot_scf_init(mos,orthogonality_metric, pscp=pscp, fscp=ks_scp, &
                        qs_ot_env=scf_env%qs_ot_env, error=error)
                ELSE
                   CALL ot_scf_init(mo_array = mos, matrix_s = orthogonality_metric, & 
                     qs_ot_env = scf_env%qs_ot_env,error=error)
                END IF
                DO ispin=1,SIZE(scf_env%qs_ot_env)
                   CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin),&
                        scf_env%ot_preconditioner(1)%preconditioner,error=error)
                ENDDO
             CASE DEFAULT
                IF ( scp_nddo ) THEN
                   CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error = error )
                   CALL get_scp_env ( scp_env = scp_env, pscp=pscp, ks_scp=ks_scp, error = error )
                   CALL ot_scf_init(mos,orthogonality_metric, pscp=pscp, fscp=ks_scp, &
                        qs_ot_env=scf_env%qs_ot_env, error=error)
                ELSE
                   CALL ot_scf_init(mo_array = mos, matrix_s = orthogonality_metric, & 
                        qs_ot_env = scf_env%qs_ot_env,error=error)
                END IF
             END SELECT
          ELSE
             SELECT CASE(scf_env%qs_ot_env(1)%settings%preconditioner_type)
             CASE(ot_precond_none)
                ! **** SCP
                IF ( scp_dft ) THEN
                   CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error = error )
                   CALL get_scp_env ( scp_env = scp_env, aux_coeff_set = aux_coeff_set, error = error )
                   CALL ot_scf_init(mos,orthogonality_metric, aux_coeff_set=aux_coeff_set, &
                        qs_ot_env=scf_env%qs_ot_env, error=error)
                ELSE
                   CALL ot_scf_init(mo_array = mos, matrix_s = orthogonality_metric, & 
                        qs_ot_env = scf_env%qs_ot_env,error=error)
                ENDIF
                ! **** SCP
             CASE(ot_precond_full_all,ot_precond_full_single_inverse)
                CPPrecondition(.NOT. has_unit_metric,cp_failure_level,routineP,error,failure)
                ! this will rotate the MOs to be eigen states, which is not compatible with rotation
                ! e.g. mo_derivs here do not yet include potentially different occupations numbers
                do_rotation=scf_env%qs_ot_env(1)%settings%do_rotation
                CPPrecondition(.NOT.do_rotation,cp_failure_level,routineP,error,failure)
                DO ispin=1,dft_control%nspins
                   CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
                   CALL calculate_subspace_eigenvalues(mo_coeff,matrix_ks(ispin)%matrix,&
                        para_env=mo_coeff%matrix_struct%para_env, &
                        do_rotation = .TRUE., co_rotate=qs_env%mo_derivs(ispin)%matrix,error=error)
                   CALL make_preconditioner(scf_env%ot_preconditioner(ispin)%preconditioner, &
                        scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                        scf_env%qs_ot_env(1)%settings%precond_solver_type, &
                        matrix_ks(ispin)%matrix,matrix_s(1)%matrix,&
                        qs_env%kinetic(1)%matrix, &
                        mo_coeff,scf_env%qs_ot_env(1)%settings%energy_gap,error=error)
                ENDDO
                ! **** SCP
                IF ( scp_dft ) THEN
                   CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error = error )
                   CALL get_scp_env ( scp_env = scp_env, aux_coeff_set = aux_coeff_set, error = error )
                   CALL ot_scf_init(mos,orthogonality_metric, aux_coeff_set=aux_coeff_set, &
                                    qs_ot_env=scf_env%qs_ot_env, error=error)
                ELSE
                   CALL ot_scf_init(mo_array = mos, matrix_s = orthogonality_metric, & 
                        qs_ot_env = scf_env%qs_ot_env,error=error)
                ENDIF
                ! **** SCP
                DO ispin=1,SIZE(scf_env%qs_ot_env)
                   CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin),&
                        scf_env%ot_preconditioner(ispin)%preconditioner,error=error)
                ENDDO
             CASE DEFAULT
                ! constuct the preconditioner based on the first spin set.
                ! this might or might not be reasonable. Presumably it is.
                CPPrecondition(.NOT. has_unit_metric,cp_failure_level,routineP,error,failure)
                CALL get_mo_set(mo_set=mos(1)%mo_set, mo_coeff=mo_coeff)
                CALL make_preconditioner(scf_env%ot_preconditioner(1)%preconditioner, &
                     scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                     scf_env%qs_ot_env(1)%settings%precond_solver_type, &
                     matrix_ks(1)%matrix,matrix_s(1)%matrix,&
                     qs_env%kinetic(1)%matrix, &
                     mo_coeff,scf_env%qs_ot_env(1)%settings%energy_gap,error=error)
                ! **** SCP
                IF ( scp_dft ) THEN
                   CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error = error )
                   CALL get_scp_env ( scp_env = scp_env, aux_coeff_set = aux_coeff_set, error = error )
                   CALL ot_scf_init(mos,orthogonality_metric, aux_coeff_set=aux_coeff_set, &
                                   qs_ot_env=scf_env%qs_ot_env, error=error)
                ELSE
                   CALL ot_scf_init(mo_array = mos, matrix_s = orthogonality_metric, & 
                        qs_ot_env = scf_env%qs_ot_env,error=error)
                ENDIF
                ! **** SCP
                DO ispin=1,SIZE(scf_env%qs_ot_env)
                   CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin),&
                        scf_env%ot_preconditioner(1)%preconditioner,error=error)
                ENDDO
             END SELECT
          END IF
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
!> \brief performs those initialisations that need to be done only once
!>       (e.g. that only depend on the atomic positions)
!>       this will be called in scf
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE init_scf_run(scf_env, qs_env, scf_section,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_scf_run', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, i, ispin, iw, &
                                                nao, ndep, nelectron_spin, &
                                                output_unit
    LOGICAL :: cho_inverse, failure, gth_potential_present, id_equal, &
      my_transition_potential, s_minus_half_available
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset(routineN,handle)

    NULLIFY(scf_control, atomic_kind_set, matrix_h, matrix_s, matrix_ks, &
         dft_control, mos, rho)
    failure=.FALSE.
    my_transition_potential = .FALSE.
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    para_env=>qs_env%para_env

    s_minus_half_available = .FALSE.

    CALL scf_env_check_i_alloc(scf_env=scf_env, qs_env=qs_env,&
         error=error)

    CALL get_qs_env(qs_env=qs_env,&
         scf_control=scf_control,&
         dft_control=dft_control,&
         atomic_kind_set=atomic_kind_set,&
         mos=mos,matrix_ks=matrix_ks, rho=rho, &
         matrix_h=matrix_h,matrix_s=matrix_s, error=error)

    ! some special checks to eliminate the first problems with restricted
    ! should move earlier I think
    IF (dft_control%restricted) THEN
       IF (scf_env%method .NE. ot_method_nr) THEN
          IF (output_unit>0) WRITE(output_unit,*) "OT only for restricted"
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
       CPPrecondition(dft_control%nspins.EQ.2,cp_failure_level,routineP,error,failure)
    ENDIF

    ! do some assertions here on these matrices having the same structure
    id_equal=ALL((/((matrix_s(i)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id),i=1,SIZE(matrix_s))/))
    CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    id_equal=ALL((/((rho%rho_ao(i)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id),i=1,SIZE(rho%rho_ao))/))
    CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    id_equal=ALL((/((matrix_h(i)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id),i=1,SIZE(matrix_h))/))
    CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    id_equal=ALL((/((matrix_ks(i)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id),i=1,SIZE(matrix_ks))/))
    CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         gth_potential_present=gth_potential_present)

! updates the total number of electrons
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         nelectron=scf_env%nelectron)
    scf_env%nelectron = scf_env%nelectron - dft_control%charge
! print occupation numbers
    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".scfLog",error=error)
    IF (output_unit>0) THEN
       DO ispin=1,dft_control%nspins
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,&
               homo=homo, &
               nelectron=nelectron_spin, &
               nao=nao)
          IF (dft_control%nspins > 1) THEN
            WRITE (UNIT=output_unit,FMT="(/,T3,A,I2)") "Spin",ispin
          END IF
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T31,I10))")&
            "Number of electrons:",nelectron_spin,&
            "Number of occupied orbitals:",homo
       ENDDO
       WRITE (UNIT=output_unit,FMT="((T3,A,T31,I10))")&
         "Number of orbital functions:",nao
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)

! calc ortho matrix
    ndep = 0
    IF (scf_env%method /= ot_method_nr .AND. &
        scf_env%method /= ot_diag_method_nr .AND. &
        scf_env%method /= special_diag_method_nr ) THEN
       CALL copy_sm_to_fm(matrix_s(1)%matrix,scf_env%ortho,error=error)
       cho_inverse = .FALSE.
       IF (scf_env%cholesky_method>cholesky_off) THEN
           CALL cp_fm_cholesky_decompose(scf_env%ortho,error=error)
           IF(scf_env%cholesky_method==cholesky_inverse) THEN
             CALL cp_fm_triangular_invert(scf_env%ortho,error=error)
           END IF
       ELSE
           CALL cp_fm_power(scf_env%ortho,scf_env%scf_work2,-0.5_dp,&
                           scf_control%eps_eigval,ndep,error=error)
           s_minus_half_available = .TRUE.
       END IF
       
       IF (BTEST(cp_print_key_should_output(logger%iter_info,&
            qs_env%input,"DFT%PRINT%AO_MATRICES/ORTHO",error=error),cp_p_file)) THEN
          iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/ORTHO",&
               extension=".Log",error=error)
          CALL write_fm_with_basis_info(scf_env%ortho,4,6,qs_env,para_env,output_unit=iw,error=error)
          CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
               "DFT%PRINT%AO_MATRICES/ORTHO", error=error)
       END IF
    END IF

    CALL get_mo_set(mo_set=mos(1)%mo_set,nao=nao)

    ! DFT+U methods based on Lowdin charges need S^(1/2)
    IF (dft_control%dft_plus_u) THEN
      IF (dft_control%plus_u_method_id == plus_u_lowdin) THEN
        IF (s_minus_half_available) THEN
          CALL cp_sm_fm_multiply(matrix_s(1)%matrix,scf_env%ortho,scf_env%s_half,&
                                 nao,error=error)
        ELSE
          CALL copy_sm_to_fm(matrix_s(1)%matrix,scf_env%s_half,error=error)
          CALL cp_fm_power(scf_env%s_half,scf_env%scf_work2,0.5_dp,&
                           scf_control%eps_eigval,ndep,error=error)
        END IF
      END IF
    END IF

    ! The P_mixing in gspace does not work "yet"
    ! The matrix formulation of G-space mixing needs S^(-1)
    IF (scf_env%mixing_method==gspace_mixing_nr .AND. .FALSE.) THEN
      IF (s_minus_half_available) THEN
          CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp, &
               scf_env%ortho,scf_env%ortho,0.0_dp,scf_env%s_minus_one,&
               error=error)
      ELSE
        CALL copy_sm_to_fm(matrix_s(1)%matrix,scf_env%s_minus_one,error=error)
        CALL cp_fm_power(scf_env%s_minus_one,scf_env%scf_work2,-1.0_dp,&
                         scf_control%eps_eigval,ndep,error=error)
      END IF
    END IF

    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
               extension=".scfLog",error=error)
    IF (output_unit > 0) THEN
       WRITE (UNIT=output_unit,FMT="(/,T3,A,I10)")&
         "Number of independent orbital functions:",nao - ndep
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
          "PRINT%PROGRAM_RUN_INFO", error=error)

    ! extrapolate outer loop variables
    IF (scf_control%outer_scf%have_scf) THEN
       CALL outer_loop_extrapolate(qs_env,error)
    ENDIF

    my_transition_potential = (dft_control%do_xas_calculation .AND. dft_control%xas_estate>0)
    ! initializes rho and the mos
    IF( my_transition_potential) THEN
      ! initialize the density with the localized mos
        CALL xas_initialize_rho(qs_env,error=error)
    ELSE
      CALL scf_env_initial_rho_setup(scf_env,qs_env=qs_env,&
           scf_section=scf_section, error=error)
    END IF

    ! Frozen density approximation
    IF (ASSOCIATED(qs_env%wf_history)) THEN
      IF (qs_env%wf_history%interpolation_method_nr==wfi_frozen_method_nr) THEN
        IF (.NOT. ASSOCIATED(qs_env%wf_history%past_states(1)%snapshot)) THEN
          CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp, &
                error=error)
          CALL duplicate_rho_type(rho_input=rho, &
                rho_output=qs_env%wf_history%past_states(1)%snapshot%rho_frozen, &
                qs_env=qs_env, error=error)
        END IF
      END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE init_scf_run

! *****************************************************************************
!> \brief checks the allocation status of the needed matrixes, and if necessary
!>      allocate them
!> \param scf_env the scf_env to be checked
!> \param qs_env the qs_env, the scf_env lives in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf_env_check_i_alloc(scf_env,qs_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scf_env_check_i_alloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, i, ispin, nao, &
                                                nbroy, nhistory, nrow_block, &
                                                nvariables, stat
    LOGICAL                                  :: failure, &
                                                gth_potential_present, &
                                                my_transition_potential, &
                                                uniform_occupation
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs, mo_derivs_aux_fit
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools, &
                                                ao_mo_fm_pools_aux_fit
    TYPE(cp_fm_struct_type), POINTER         :: ao_ao_fmstruct, ao_mo_fmstruct
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools, S_sm_pools_aux_fit
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset(routineN,handle)

    NULLIFY(matrix_ks, ao_mo_fm_pools, matrix_s, ao_mo_fmstruct, ao_ao_fmstruct,&
            dft_control, mos, occupation_numbers, S_sm_pools, ks_env,&
            ao_mo_fm_pools_aux_fit, matrix_ks_aux_fit, S_sm_pools_aux_fit)
    NULLIFY(atomic_kind_set, mo_derivs, mo_coeff, scf_control, mo_coeff_aux_fit,&
            mo_derivs_aux_fit)

    failure=.FALSE.

    my_transition_potential = .FALSE.
    uniform_occupation = .TRUE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         mos=mos,&
         mos_aux_fit=mos_aux_fit,&
         matrix_ks=matrix_ks,&
         ks_env=ks_env,&
         atomic_kind_set=atomic_kind_set,&
         matrix_s=matrix_s,error=error)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         gth_potential_present=gth_potential_present)
    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools,&
         S_sm_pools=S_sm_pools,&
         error=error)
    IF( dft_control%do_wfn_fitting) THEN
      CALL mpools_get(qs_env%mpools_aux_fit, ao_mo_fm_pools=ao_mo_fm_pools_aux_fit,&
           S_sm_pools=S_sm_pools_aux_fit,&
           error=error)
    END IF

    ! very first tests for xas

    my_transition_potential = (dft_control%do_xas_calculation .AND. dft_control%xas_estate>0)
!   *** finish initialization of the MOs ***
    CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       DO ispin=1,SIZE(mos)
          CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo)
          IF (.NOT.ASSOCIATED(mo_coeff)) THEN
             CALL init_mo_set(mos(ispin)%mo_set,&
                  ao_mo_fm_pools(ispin)%pool,&
                  name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                  "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                  error=error)
          END IF
          IF(my_transition_potential .AND. ispin==1) THEN
            CALL get_mo_set(mos(ispin)%mo_set,&
                 occupation_numbers=occupation_numbers)
            occupation_numbers(dft_control%xas_estate) = &
                 dft_control%xas_control%occ_estate
            occupation_numbers(homo) = dft_control%xas_control%occ_homo
            uniform_occupation = .FALSE.

            CALL set_mo_set(mos(ispin)%mo_set,&
                 uniform_occupation=uniform_occupation,error=error)
          END IF
       END DO
    END IF

!   *** finish initialization of the MOs for wfn fitting ***
    IF(dft_control%do_wfn_fitting) THEN
      CPPrecondition(ASSOCIATED(mos_aux_fit),cp_failure_level,routineP,error,failure)
      IF (.NOT.failure) THEN
         DO ispin=1,SIZE(mos_aux_fit)
            CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit,homo=homo)
            IF (.NOT.ASSOCIATED(mo_coeff_aux_fit)) THEN
               CALL init_mo_set(mos_aux_fit(ispin)%mo_set,&
                    ao_mo_fm_pools_aux_fit(ispin)%pool,&
                    name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                    "%mo_aux_fit"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                    error=error)
            END IF
            IF(my_transition_potential .AND. ispin==1) THEN
              CALL get_mo_set(mos_aux_fit(ispin)%mo_set,&
                   occupation_numbers=occupation_numbers)
              occupation_numbers(dft_control%xas_estate) = &
                   dft_control%xas_control%occ_estate
              occupation_numbers(homo) = dft_control%xas_control%occ_homo
              uniform_occupation = .FALSE.

              CALL set_mo_set(mos_aux_fit(ispin)%mo_set,&
                   uniform_occupation=uniform_occupation,error=error)
            END IF
         END DO
      END IF
    END IF

!   *** get the mo_derivs OK if needed ***
    IF (qs_env%requires_mo_derivs) THEN
     CALL get_qs_env(qs_env,mo_derivs=mo_derivs,error=error)
     IF (.NOT.ASSOCIATED(mo_derivs)) THEN
       IF (dft_control%restricted) THEN ! right now, there might be more mos than needed derivs
          ALLOCATE(mo_derivs(1))
          CALL get_mo_set(mos(1)%mo_set,mo_coeff=mo_coeff)
          CALL cp_fm_create(mo_derivs(1)%matrix,mo_coeff%matrix_struct,error=error)
       ELSE
          ALLOCATE(mo_derivs(dft_control%nspins))
          DO ispin=1,dft_control%nspins
             CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
             CALL cp_fm_create(mo_derivs(ispin)%matrix,mo_coeff%matrix_struct,error=error)
          ENDDO
       ENDIF
       CALL set_qs_env(qs_env,mo_derivs=mo_derivs,error=error)
     ENDIF
    ELSE
       ! nothing should be done
    ENDIF

!   *** get the mo_derivs OK if needed ***
    IF( dft_control%do_wfn_fitting) THEN
      IF (qs_env%requires_mo_derivs) THEN
       CALL get_qs_env(qs_env,mo_derivs_aux_fit=mo_derivs_aux_fit,error=error)
       IF (.NOT.ASSOCIATED(mo_derivs_aux_fit)) THEN
         IF (dft_control%restricted) THEN ! right now, there might be more mos than needed derivs
            ALLOCATE(mo_derivs_aux_fit(1))
            CALL get_mo_set(mos_aux_fit(1)%mo_set,mo_coeff=mo_coeff_aux_fit)
            CALL cp_fm_create(mo_derivs_aux_fit(1)%matrix,mo_coeff_aux_fit%matrix_struct,error=error)
         ELSE
            ALLOCATE(mo_derivs_aux_fit(dft_control%nspins))
            DO ispin=1,dft_control%nspins
               CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)
               CALL cp_fm_create(mo_derivs_aux_fit(ispin)%matrix,mo_coeff_aux_fit%matrix_struct,error=error)
            ENDDO
         ENDIF
         CALL set_qs_env(qs_env,mo_derivs_aux_fit=mo_derivs_aux_fit,error=error)
       ENDIF
      ELSE
         ! nothing should be done
      ENDIF
    END IF
!   *** Allocate the distributed SCF matrices ***

    IF ((.NOT.ASSOCIATED(scf_env%scf_work1)).OR.&
        (.NOT.ASSOCIATED(scf_env%scf_work2)).OR.&
        (.NOT.ASSOCIATED(scf_env%ortho)).OR.&
        (.NOT.ASSOCIATED(scf_env%s_half))) THEN

       ao_mo_fmstruct => fm_pool_get_el_struct(ao_mo_fm_pools(1)%pool,&
            error=error)
       CALL cp_fm_struct_get(ao_mo_fmstruct, nrow_block=nrow_block,&
            error=error)
       CALL get_mo_set(mos(1)%mo_set,nao=nao)
       CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct,&
            nrow_block=nrow_block,ncol_block=nrow_block,&
            nrow_global=nao, ncol_global=nao,&
            template_fmstruct=ao_mo_fmstruct, error=error)

       IF (.NOT.ASSOCIATED(scf_env%scf_work1)) THEN
          IF (scf_env%method /= ot_method_nr) THEN
             ALLOCATE(scf_env%scf_work1(dft_control%nspins), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(scf_env%scf_work1)) THEN
             DO ispin=1,SIZE(scf_env%scf_work1)
                NULLIFY(scf_env%scf_work1(ispin)%matrix)
                CALL cp_fm_create(scf_env%scf_work1(ispin)%matrix,&
                     matrix_struct=ao_ao_fmstruct,&
                     name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                     "WORK_MATRIX-1-"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                     error=error)
             ENDDO
          ENDIF
       END IF
       IF ((.NOT.ASSOCIATED(scf_env%ortho)).AND.&
           ((scf_env%method /= ot_method_nr).AND.&
            (scf_env%method /= special_diag_method_nr))) THEN
          CALL cp_fm_create(scf_env%ortho,&
               matrix_struct=ao_ao_fmstruct,&
               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
               "ORTHO_MATRIX",&
               error=error)
       END IF
       IF ((.NOT.ASSOCIATED(scf_env%scf_work2)).AND.&
           (scf_env%method /= ot_method_nr)) THEN
          CALL cp_fm_create(scf_env%scf_work2,&
               matrix_struct=ao_ao_fmstruct,&
               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
               "WORK_MATRIX-2",&
               error=error)
       END IF
       IF (dft_control%dft_plus_u) THEN
         IF (dft_control%plus_u_method_id == plus_u_lowdin) THEN
           IF (.NOT.ASSOCIATED(scf_env%scf_work2)) THEN
             CALL cp_fm_create(scf_env%scf_work2,&
                               matrix_struct=ao_ao_fmstruct,&
                               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                               "WORK_MATRIX-2",&
                               error=error)
           END IF
           IF (.NOT.ASSOCIATED(scf_env%s_half)) THEN
             CALL cp_fm_create(scf_env%s_half,&
                               matrix_struct=ao_ao_fmstruct,&
                               name="S**(1/2) MATRIX",&
                               error=error)
           END IF
         END IF
       END IF

      ! The P_mixing in gspace does not work "yet"
       IF (scf_env%mixing_method==gspace_mixing_nr .AND. .FALSE.) THEN
         IF (.NOT.ASSOCIATED(scf_env%s_minus_one)) THEN
            CALL cp_fm_create(scf_env%s_minus_one,&
                              matrix_struct=ao_ao_fmstruct,&
                              name="S**(-1) MATRIX",&
                              error=error)
         END IF
       END IF

       CALL cp_fm_struct_release(ao_ao_fmstruct,error=error)

    END IF

!   *** Allocate matrix_ks and put it in the QS environment ***

    IF (.not.ASSOCIATED(matrix_ks)) THEN
       CALL sm_pools_create_matrix_vect(S_sm_pools,matrix_ks,&
            name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
            "KOHN-SHAM_MATRIX",&
            error=error)
       CALL set_qs_env(qs_env=qs_env,&
            matrix_ks=matrix_ks,error=error)
    END IF

!   *** Allocate matrix_ks_aux_fit if requested and put it in the QS environment ***

    IF( dft_control%do_wfn_fitting) THEN
      IF (.not.ASSOCIATED(matrix_ks_aux_fit)) THEN
         CALL sm_pools_create_matrix_vect(S_sm_pools_aux_fit,matrix_ks_aux_fit,&
              name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
              "KOHN-SHAM_MATRIX for wfn fitting",&
              error=error)
         CALL set_qs_env(qs_env=qs_env,&
              matrix_ks_aux_fit=matrix_ks_aux_fit,error=error)
      END IF
    END IF


    IF (scf_env%mixing_method>0) THEN
!   *** allocate p_mix_new ***
       IF (.NOT.ASSOCIATED(scf_env%p_mix_new)) THEN
          CALL sm_pools_create_matrix_vect(S_sm_pools,scf_env%p_mix_new,&
               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
               "DENSITY",&
               error=error)
       END IF
!   *** allocate p_delta ***
       IF (scf_env%mixing_method>=gspace_mixing_nr) THEN
         IF(.NOT.ASSOCIATED(scf_env%p_delta)) THEN
            CALL sm_pools_create_matrix_vect(S_sm_pools,scf_env%p_delta,&
                 name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                 "DENSITY DIFFERENCE",error=error)
         END IF
       END IF
!   *** allocate kerker buffer ***
       IF (scf_env%mixing_method==gspace_mixing_nr) THEN
         IF(.NOT. ASSOCIATED(scf_env%rhoin_buffer)) THEN
           ALLOCATE(scf_env%rhoin_buffer(1,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO ispin = 1,dft_control%nspins
            NULLIFY(scf_env%rhoin_buffer(1,ispin)%cc)
           END DO 
         END IF
       END IF
!   *** allocate pulay buffer ***
       IF (scf_env%mixing_method==pulay_mixing_nr) THEN
         scf_env%mixing_ncall = 0
         IF(.NOT. ASSOCIATED(scf_env%pulay_matrix)) THEN
           ALLOCATE(scf_env%pulay_matrix(scf_env%npulay,scf_env%npulay),STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%rhoin_buffer(scf_env%npulay,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%res_buffer(scf_env%npulay,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO ispin = 1,dft_control%nspins
           DO i = 1,scf_env%npulay
              NULLIFY(scf_env%rhoin_buffer(i,ispin)%cc)
              NULLIFY(scf_env%res_buffer(i,ispin)%cc)
           END DO
           END DO
         END IF
       END IF
!   *** allocate broyden buffer ***
       IF (scf_env%mixing_method==broyden_mixing_nr) THEN
         scf_env%mixing_ncall = 0
         IF(.NOT. ASSOCIATED(scf_env%rhoin_buffer)) THEN
           ALLOCATE(scf_env%rhoin_buffer(scf_env%nbroy,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%res_buffer(scf_env%nbroy,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%drho_buffer(scf_env%nbroy,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%broy_last_res(dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO ispin = 1,dft_control%nspins
           DO i = 1,scf_env%nbroy
              NULLIFY(scf_env%rhoin_buffer(i,ispin)%cc)
              NULLIFY(scf_env%res_buffer(i,ispin)%cc)
              NULLIFY(scf_env%drho_buffer(i,ispin)%cc)
           END DO
              NULLIFY(scf_env%broy_last_res(ispin)%cc)
           END DO
         END IF
       END IF
!   *** allocate broyden buffer ***
       IF (scf_env%mixing_method==broyden_mixing_new_nr) THEN
         scf_env%mixing_ncall = 0
         nbroy = scf_env%nbroy 
         IF(.NOT. ASSOCIATED(scf_env%rhoin_buffer)) THEN
           ALLOCATE(scf_env%rhoin_buffer(nbroy,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%broy_store%last_res(dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%broy_store%delta_res(nbroy-1,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%broy_store%u_vec(nbroy-1,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%broy_store%z_vec(nbroy-1,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%broy_store%weight(nbroy,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%broy_store%fmat(nbroy-1,nbroy-1,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%broy_store%smat(nbroy-1,nbroy-1,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           ALLOCATE(scf_env%broy_store%gmat(nbroy-1,nbroy-1,dft_control%nspins), STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO ispin = 1,dft_control%nspins
           DO i = 1,nbroy
              NULLIFY(scf_env%rhoin_buffer(i,ispin)%cc)
           END DO
           DO i = 1,nbroy-1
              NULLIFY(scf_env%broy_store%delta_res(i,ispin)%cc)
              NULLIFY(scf_env%broy_store%u_vec(i,ispin)%cc)
              NULLIFY(scf_env%broy_store%z_vec(i,ispin)%cc)
           END DO
           NULLIFY(scf_env%broy_store%last_res(ispin)%cc)
           END DO
         END IF
       END IF
    END IF

!   *** allocate the ks env **
    IF (.NOT.ASSOCIATED(ks_env)) THEN
       CALL qs_ks_create(ks_env,qs_env=qs_env,error=error)
       CALL set_qs_env(qs_env, ks_env=ks_env,error=error)
       CALL qs_ks_release(ks_env,error=error)
    END IF

    ! If there is an outer scf loop allocate the space for the variables
    CALL get_qs_env(qs_env=qs_env,scf_control=scf_control, dft_control=dft_control,error=error)
    IF (scf_control%outer_scf%have_scf) THEN
       nhistory = scf_control%outer_scf%max_scf + 1
       nvariables = outer_loop_variables_count(scf_control,error)
       ALLOCATE(scf_env%outer_scf%variables(nvariables,nhistory),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(scf_env%outer_scf%count(nhistory),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       scf_env%outer_scf%count=0
       ALLOCATE(scf_env%outer_scf%gradient(nvariables,nhistory),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(scf_env%outer_scf%energy(nhistory),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE scf_env_check_i_alloc

! *****************************************************************************
!> \brief Initializes rho and the mos, so that an scf cycle can start
!> \param scf_env the scf env in which to do the scf
!> \param qs_env the qs env the scf_env lives in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf_env_initial_rho_setup(scf_env, qs_env, scf_section, error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_initial_rho_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: extrapolation_method_nr, &
                                                handle, ib, ig, ispin, istat, &
                                                ng, nmo, nspin, output_unit
    LOGICAL                                  :: failure, gapw, gapw_xc, &
                                                orthogonal_wf
    REAL(dp)                                 :: bconst, beta, fdamp, g2min, &
                                                kmin
    REAL(dp), DIMENSION(:), POINTER          :: g2
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(mo_coeff )
    gapw = qs_env%dft_control%qs_control%gapw
    gapw_xc = qs_env%dft_control%qs_control%gapw_xc
    logger => cp_error_get_logger(error)
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       extrapolation_method_nr=wfi_use_guess_method_nr
       IF (ASSOCIATED(qs_env%wf_history)) THEN 
          CALL wfi_extrapolate(qs_env%wf_history, &
               qs_env=qs_env, dt=1.0_dp, &
               extrapolation_method_nr=extrapolation_method_nr,&
               orthogonal_wf=orthogonal_wf, error=error)
          ! wfi_use_guess_method_nr the wavefunctions are not yet initialized
          IF ((.NOT.orthogonal_wf).AND.&
              (scf_env%method == ot_method_nr).AND.&
              (.NOT.(extrapolation_method_nr == wfi_use_guess_method_nr))) THEN 
             DO ispin=1,SIZE(qs_env%mos)
                CALL get_mo_set(qs_env%mos(ispin)%mo_set, &
                     mo_coeff=mo_coeff, nmo=nmo)
                CALL qs_env_reorthogonalize_vectors(qs_env, &
                     v_matrix=mo_coeff, n_col=nmo,&
                     error=error)
                CALL set_mo_occupation(mo_set=qs_env%mos(ispin)%mo_set, smear=qs_env%scf_control%smear, error=error)
             END DO
          END IF
       END IF
       output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
            extension=".scfLog",error=error)
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,T3,A)")&
               "Extrapolation method: "//&
               TRIM(wfi_get_method_label(extrapolation_method_nr,error=error))
       END IF

       CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)
       IF (extrapolation_method_nr==wfi_use_guess_method_nr) THEN
          CALL calculate_first_density_matrix(scf_env=scf_env,qs_env=qs_env,error=error)
          IF (.NOT.(qs_env%scf_control%density_guess==densities_guess)) THEN
            CALL qs_rho_update_rho(qs_env%rho,qs_env=qs_env, error=error)
            CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,&
                 error=error)
          END IF
       END IF
            IF(qs_env%scf_env%mixing_method>1) THEN
               nspin = SIZE(qs_env%rho%rho_g) 
               ng = SIZE(qs_env%rho%rho_g(1)%pw%pw_grid%gsq,1)
               IF(.NOT. ASSOCIATED(qs_env%scf_env%kerker_factor)) THEN
                 ALLOCATE(qs_env%scf_env%kerker_factor(ng),STAT=istat)
                 CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
               END IF
               g2 =>qs_env%rho%rho_g(1)%pw%pw_grid%gsq
               beta = qs_env%scf_env%p_mix_beta
               kmin = 0.1_dp
               DO ig =1,ng
                 qs_env%scf_env%kerker_factor(ig) = MAX(g2(ig)/(g2(ig)+beta*beta),kmin)
               END DO
               IF (qs_env%scf_env%mixing_method==gspace_mixing_nr) THEN
                  DO ispin = 1,nspin
                  IF(.NOT. ASSOCIATED(qs_env%scf_env%rhoin_buffer(1,1)%cc)) THEN
                     ALLOCATE(qs_env%scf_env%rhoin_buffer(1,ispin)%cc(ng),STAT=istat)
                     CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                  END IF
                  qs_env%scf_env%rhoin_buffer(1,ispin)%cc(1:ng) = qs_env%rho%rho_g(ispin)%pw%cc(1:ng)
                  END DO
               ELSEIF(qs_env%scf_env%mixing_method==pulay_mixing_nr) THEN
                   DO ispin = 1,nspin
                      IF(.NOT. ASSOCIATED(qs_env%scf_env%rhoin_buffer(1,1)%cc)) THEN
                        DO ib = 1,qs_env%scf_env%npulay
                           ALLOCATE(qs_env%scf_env%rhoin_buffer(ib,ispin)%cc(ng),STAT=istat)  
                           CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                           ALLOCATE(qs_env%scf_env%res_buffer(ib,ispin)%cc(ng),STAT=istat)  
                           CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                        END DO
                      END IF
                      qs_env%scf_env%rhoin_buffer(1,ispin)%cc(1:ng) = qs_env%rho%rho_g(ispin)%pw%cc(1:ng)
                   END DO
               ELSEIF(qs_env%scf_env%mixing_method==broyden_mixing_nr) THEN
                  DO ispin = 1,nspin
                    IF(.NOT. ASSOCIATED(qs_env%scf_env%rhoin_buffer(1,ispin)%cc)) THEN
                    DO ib = 1,qs_env%scf_env%nbroy
                      ALLOCATE(qs_env%scf_env%rhoin_buffer(ib,ispin)%cc(ng),STAT=istat)  
                      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      ALLOCATE(qs_env%scf_env%res_buffer(ib,ispin)%cc(ng),STAT=istat)  
                       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      ALLOCATE(qs_env%scf_env%drho_buffer(ib,ispin)%cc(ng),STAT=istat)  
                      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                    END DO
                    ALLOCATE(qs_env%scf_env%broy_last_res(ispin)%cc(ng),STAT=istat)  
                    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                    END IF
                    qs_env%scf_env%rhoin_buffer(1,ispin)%cc(1:ng) = qs_env%rho%rho_g(ispin)%pw%cc(1:ng)
                  END DO
               ELSEIF(qs_env%scf_env%mixing_method==broyden_mixing_new_nr) THEN
       CPPrecondition(ASSOCIATED(qs_env%scf_env%broy_store),cp_failure_level,routineP,error,failure)
                  DO ispin = 1,nspin
                    IF(.NOT. ASSOCIATED(qs_env%scf_env%rhoin_buffer(1,ispin)%cc)) THEN
                       DO ib = 1,qs_env%scf_env%nbroy-1
                         ALLOCATE(qs_env%scf_env%broy_store%delta_res(ib,ispin)%cc(ng),STAT=istat)  
                         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                         ALLOCATE(qs_env%scf_env%broy_store%u_vec(ib,ispin)%cc(ng),STAT=istat)  
                          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                         ALLOCATE(qs_env%scf_env%broy_store%z_vec(ib,ispin)%cc(ng),STAT=istat)  
                         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                         ALLOCATE(qs_env%scf_env%rhoin_buffer(ib,ispin)%cc(ng),STAT=istat)  
                         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                       END DO
                       ALLOCATE(qs_env%scf_env%rhoin_buffer(qs_env%scf_env%nbroy,ispin)%cc(ng),STAT=istat)  
                       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                       ALLOCATE(qs_env%scf_env%broy_store%last_res(ispin)%cc(ng),STAT=istat)  
                       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                       ALLOCATE(qs_env%scf_env%broy_store%p_metric(ng),STAT=istat)
                       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                       IF(qs_env%scf_env%broy_store%p_metric_method==1) THEN
                         bconst = qs_env%scf_env%broy_store%bconst
                         g2min = 1.E30_dp
                         DO ig = 1,ng
                           IF(g2(ig)>1.E-10_dp) g2min =  MIN(g2min,g2(ig))
                         END DO
                         CALL mp_min(g2min,qs_env%para_env%group)
                         fdamp = (bconst-1.0_dp)*g2min
                       ELSEIF(qs_env%scf_env%broy_store%p_metric_method==2) THEN
                         fdamp = beta*beta
                       END IF
                       DO ig = 1,ng
                           qs_env%scf_env%broy_store%p_metric(ig) = (g2(ig)+fdamp)/MAX(g2(ig),1.E-10_dp)
                       END DO
                       IF(qs_env%rho%rho_g(1)%pw%pw_grid%have_g0) qs_env%scf_env%broy_store%p_metric(1) = 0.0_dp
                    END IF
                    qs_env%scf_env%rhoin_buffer(1,ispin)%cc(1:ng) = qs_env%rho%rho_g(ispin)%pw%cc(1:ng)
                  END DO
               END IF
            END IF
! *** SCP
       IF ( qs_env % dft_control % scp ) THEN
         CALL update_rhoscp ( qs_env = qs_env, error = error )
         CALL qs_ks_scp_did_change ( qs_env=qs_env, rho_changed = .TRUE., error = error )
       END IF
! *** SCP
    END IF

    CALL timestop(handle)

  END SUBROUTINE scf_env_initial_rho_setup

! *****************************************************************************
!> \brief perform cleanup operations (like releasing temporary storage)
!>      at the end of the scf
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf_env_cleanup(scf_env,qs_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, stat
    LOGICAL                                  :: failure
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(S_sm_pools)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL mpools_get(qs_env%mpools,S_sm_pools=S_sm_pools,error=error)

!   *** Release SCF work storage ***

       IF (ASSOCIATED(scf_env%scf_work1)) THEN
          DO ispin=1,SIZE(scf_env%scf_work1)
             CALL cp_fm_release(scf_env%scf_work1(ispin)%matrix,error=error)
          ENDDO
          DEALLOCATE(scf_env%scf_work1)
       ENDIF
       IF (ASSOCIATED(scf_env%scf_work2)) CALL cp_fm_release(scf_env%scf_work2,error)
       IF (ASSOCIATED(scf_env%ortho)) CALL cp_fm_release(scf_env%ortho,error=error)

       IF (ASSOCIATED(scf_env%p_mix_new)) THEN
          CALL sm_pools_give_back_matrix_vect(S_sm_pools,scf_env%p_mix_new,&
               error=error)
          CALL sm_pools_flush_cache(S_sm_pools,error=error)
       END IF

       IF (ASSOCIATED(scf_env%p_delta)) THEN
          CALL sm_pools_give_back_matrix_vect(S_sm_pools,scf_env%p_delta,&
               error=error)
          CALL sm_pools_flush_cache(S_sm_pools,error=error)
       END IF

! *** method dependent cleanup
       SELECT CASE(scf_env%method)
       CASE(ot_method_nr)
          DO ispin=1,SIZE(scf_env%ot_preconditioner)
             CALL destroy_preconditioner(scf_env%ot_preconditioner(ispin)%preconditioner,error=error)
             DEALLOCATE(scf_env%ot_preconditioner(ispin)%preconditioner)
          ENDDO
          DEALLOCATE(scf_env%ot_preconditioner,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       CASE(ot_diag_method_nr)
          !
       CASE(general_diag_method_nr)
          !
       CASE(special_diag_method_nr)
          !
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

    END IF
    CALL timestop(handle)

  END SUBROUTINE scf_env_cleanup

! *****************************************************************************
!> \brief perform cleanup operations at the end of an scf loop
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE cleanup_scf_loop(scf_env,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cleanup_scf_loop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
! *** method dependent cleanup
       SELECT CASE(scf_env%method)
       CASE(ot_method_nr)
          DO ispin=1,SIZE(scf_env%qs_ot_env)
             CALL ot_scf_destroy(scf_env%qs_ot_env(ispin),error=error)
          ENDDO
          DEALLOCATE(scf_env%qs_ot_env,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       CASE(ot_diag_method_nr)
          !
       CASE(general_diag_method_nr)
          !
       CASE(special_diag_method_nr)
          !
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown scf method method:"//&
               cp_to_string(scf_env%method),error,failure)
       END SELECT

    END IF

    CALL timestop(handle)

  END SUBROUTINE cleanup_scf_loop

! *****************************************************************************
  SUBROUTINE qs_scf_loop_do_roks_diag(scf_env,mos,matrix_ks,matrix_s,&
                                      scf_control,scf_section,diis_step,&
                                      orthogonal_basis,error)

    ! Solve a set restricted open Kohn-Sham (ROKS) equations based on the
    ! alpha and beta Kohn-Sham matrices from unrestricted Kohn-Sham.
    ! Note, this is only a high-spin ROKS.

    ! Literature: - C. C. J. Roothaan, Rev. Mod. Phys. 32, 179 (1960)
    !             - M. F. Guest and V. R. Saunders, Mol. Phys. 28(3), 819 (1974)
    !             - M. Filatov and S. Shaik, Chem. Phys. Lett. 288, 689 (1998)

    ! - Creation (05.04.06,MK)
    ! - Revised (01.05.06,MK)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    LOGICAL, INTENT(IN)                      :: orthogonal_basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_loop_do_roks_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homoa, homob, imo, &
                                                nalpha, nao, nbeta, nmo
    REAL(KIND=dp)                            :: diis_error, level_shift_loc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eiga, eigb, occa, occb
    TYPE(cp_fm_type), POINTER                :: ksa, ksb, mo2ao, moa, mob, &
                                                ortho, work

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    ortho => scf_env%ortho
    work  => scf_env%scf_work2

    ksa => scf_env%scf_work1(1)%matrix
    ksb => scf_env%scf_work1(2)%matrix

    CALL copy_sm_to_fm(matrix_ks(1)%matrix,ksa,error=error)
    CALL copy_sm_to_fm(matrix_ks(2)%matrix,ksb,error=error)

    ! Get MO information

    CALL get_mo_set(mo_set=mos(1)%mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    nelectron=nalpha,&
                    homo=homoa,&
                    eigenvalues=eiga,&
                    occupation_numbers=occa,&
                    mo_coeff=moa)

    CALL get_mo_set(mo_set=mos(2)%mo_set,&
                    nelectron=nbeta,&
                    homo=homob,&
                    eigenvalues=eigb,&
                    occupation_numbers=occb,&
                    mo_coeff=mob)

    ! Define the amount of level-shifting

    IF ((scf_control%level_shift /= 0.0_dp).AND.&
        ((scf_control%density_guess == core_guess).OR.&
         (scf_control%density_guess == restart_guess).OR.&
         (scf_env%iter_count > 1))) THEN
      level_shift_loc = scf_control%level_shift
    ELSE
      level_shift_loc = 0.0_dp
    END IF

    IF ((scf_env%iter_count > 1).OR.&
        (scf_control%density_guess == core_guess).OR.&
        (scf_control%density_guess == restart_guess)) THEN

      ! Transform the spin unrestricted alpha and beta Kohn-Sham matrices
      ! from AO basis to MO basis: K(MO) = C(T)*K(AO)*C

      CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ksa,moa,0.0_dp,work,error=error)
      CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,moa,work,0.0_dp,ksa,error=error)

      CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ksb,moa,0.0_dp,work,error=error)
      CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,moa,work,0.0_dp,ksb,error=error)

      ! Combine the spin unrestricted alpha and beta Kohn-Sham matrices
      ! in the MO basis

      IF (scf_control%roks_scheme == general_roks) THEN
        CALL combine_ks_matrices(ksa,ksb,occa,occb,scf_control%roks_f,&
                                 nalpha,nbeta,error=error)
      ELSE IF (scf_control%roks_scheme == high_spin_roks) THEN
        CALL combine_ks_matrices(ksa,ksb,occa,occb,scf_control%roks_parameter,&
                                 error=error)
      ELSE
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Unknown ROKS scheme requested")
      END IF

      ! Back-transform the restricted open Kohn-Sham matrix from MO basis
      ! to AO basis

      IF (orthogonal_basis) THEN
        ! Q = C
        mo2ao => moa
      ELSE
        ! Q = S*C
        mo2ao => mob
!MK     CALL copy_sm_to_fm(matrix_s(1)%matrix,work)
!MK     CALL cp_fm_symm("L","U",nao,nao,1.0_dp,work,moa,0.0_dp,mo2ao)
        CALL cp_sm_fm_multiply(matrix_s(1)%matrix,moa,mo2ao,nao,error=error)
      END IF

      ! K(AO) = Q*K(MO)*Q(T)

      CALL cp_fm_gemm("N","T",nao,nao,nao,1.0_dp,ksa,mo2ao,0.0_dp,work,error=error)
      CALL cp_fm_gemm("N","N",nao,nao,nao,1.0_dp,mo2ao,work,0.0_dp,ksa,error=error)

    ELSE

      ! No transformation matrix available, yet. The closed shell part,
      ! i.e. the beta Kohn-Sham matrix in AO basis, is taken.
      ! There might be better choices, anyhow.

      CALL cp_fm_to_fm(ksb,ksa,error=error)

    END IF

    ! Update DIIS buffer and possibly perform DIIS extrapolation step

    IF (scf_env%iter_count > 1) THEN
      IF (orthogonal_basis) THEN
        CALL qs_diis_b_step(diis_buffer=scf_env%scf_diis_buffer,&
                            mo_array=mos,&
                            kc=scf_env%scf_work1,&
                            sc=work,&
                            delta=scf_env%iter_delta,&
                            error_max=diis_error,&
                            diis_step=diis_step,&
                            eps_diis=scf_control%eps_diis,&
                            scf_section=scf_section,&
                            roks=.TRUE.,&
                            error=error)
      ELSE
        CALL qs_diis_b_step(diis_buffer=scf_env%scf_diis_buffer,&
                            mo_array=mos,&
                            kc=scf_env%scf_work1,&
                            sc=work,&
                            delta=scf_env%iter_delta,&
                            error_max=diis_error,&
                            diis_step=diis_step,&
                            eps_diis=scf_control%eps_diis,&
                            scf_section=scf_section,&
                            s_matrix=matrix_s,&
                            roks=.TRUE.,&
                            error=error)
      END IF
    END IF

    IF (diis_step) THEN
      scf_env%iter_param = diis_error
      scf_env%iter_method = "DIIS/Diag."
    ELSE
      SELECT CASE(scf_env%mixing_method)
      CASE(direct_mixing_nr)
          scf_env%iter_method = "P_Mix/Diag."
          scf_env%iter_param = scf_env%p_mix_alpha
      CASE(gspace_mixing_nr)
          scf_env%iter_method = "Kerker/Diag."
          scf_env%iter_param = scf_env%p_mix_alpha
      CASE(pulay_mixing_nr)
          scf_env%iter_method = "Pulay/Diag."
          scf_env%iter_param = scf_env%pulay_alpha
      CASE(broyden_mixing_nr)
          scf_env%iter_method = "Broy/Diag."
          scf_env%iter_param = scf_env%p_mix_alpha
      END SELECT
    END IF

    scf_env%iter_delta = 0.0_dp

    IF (level_shift_loc /= 0.0_dp) THEN

      ! Transform the current Kohn-Sham matrix from AO to MO basis
      ! for level-shifting using the current MO set

      CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ksa,moa,0.0_dp,work,error=error)
      CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,moa,work,0.0_dp,ksa,error=error)

      ! Apply level-shifting using 50:50 split of the shift (could be relaxed)

      DO imo=homob+1,homoa
        CALL cp_fm_add_to_element(ksa,imo,imo,0.5_dp*level_shift_loc,error)
      END DO
      DO imo=homoa+1,nmo
        CALL cp_fm_add_to_element(ksa,imo,imo,level_shift_loc,error)
      END DO

    ELSE IF (.NOT.orthogonal_basis) THEN

      ! Transform the current Kohn-Sham matrix to an orthogonal basis
      SELECT CASE(scf_env%cholesky_method)
      CASE(cholesky_reduce)
        CALL cp_fm_cholesky_reduce(ksa,ortho,error=error)
      CASE(cholesky_restore)
         CALL cp_fm_upper_to_full(ksa,work,error=error)
         CALL cp_fm_cholesky_restore(ksa,nao,ortho,work,&
                "SOLVE",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,ksa,&
               "SOLVE",pos="LEFT",transa="T",error=error)
      CASE(cholesky_inverse)
         CALL cp_fm_upper_to_full(ksa,work,error=error)
         CALL cp_fm_cholesky_restore(ksa,nao,ortho,work,&
                "MULTIPLY",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,ksa,&
               "MULTIPLY",pos="LEFT",transa="T",error=error)
      CASE(cholesky_off)
        CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ksa,ortho,0.0_dp,work,error=error)
        CALL cp_fm_gemm("N","N",nao,nao,nao,1.0_dp,ortho,work,0.0_dp,ksa,error=error)
      END SELECT

    END IF

    ! Diagonalization of the ROKS operator matrix

    CALL cp_fm_syevd(ksa,work,eiga,error=error)

    ! Back-transformation of the orthonormal eigenvectors if needed

    IF (level_shift_loc /= 0.0_dp) THEN
      ! Use old MO set for back-transformation if level-shifting was applied
      CALL cp_fm_to_fm(moa,ortho,error=error)
      CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,ortho,work,0.0_dp,moa,error=error)
    ELSE
      IF (orthogonal_basis) THEN
        CALL cp_fm_to_fm(work,moa,error=error)
      ELSE
        SELECT CASE(scf_env%cholesky_method)
        CASE(cholesky_reduce,cholesky_restore)
          CALL cp_fm_cholesky_restore(work,nmo,ortho,moa,"SOLVE",error=error)
        CASE(cholesky_inverse)
          CALL cp_fm_cholesky_restore(work,nmo,ortho,moa,"MULTIPLY",error=error)
        CASE(cholesky_off) 
          CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,ortho,work,0.0_dp,moa,error=error)
        END SELECT
      END IF
    END IF

    ! Correct MO eigenvalues, if level-shifting was applied

    IF (level_shift_loc /= 0.0_dp) THEN
      DO imo=homob+1,homoa
        eiga(imo) = eiga(imo) - 0.5_dp*level_shift_loc
      END DO
      DO imo=homoa+1,nmo
        eiga(imo) = eiga(imo) - level_shift_loc
      END DO
    END IF

    ! Update also the beta MO set

    eigb(:) = eiga(:)
    CALL cp_fm_to_fm(moa,mob,error=error)

    ! Calculate the new alpha and beta density matrix

    CALL calculate_density_matrix(mos(1)%mo_set,scf_env%p_mix_new(1)%matrix,&
                                  error=error)
    CALL calculate_density_matrix(mos(2)%mo_set,scf_env%p_mix_new(2)%matrix,&
                                  error=error)

    CALL timestop(handle)

  END SUBROUTINE qs_scf_loop_do_roks_diag

END MODULE qs_scf
