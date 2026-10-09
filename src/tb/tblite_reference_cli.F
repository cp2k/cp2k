!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Optional native tblite command-line diagnostics.
!> \note This module has no dependency on the tblite library.
! **************************************************************************************************

MODULE tblite_reference_cli

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_control_types,                ONLY: dft_control_type,&
                                              xtb_reference_cli_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: debug_print_level,&
                                              high_print_level,&
                                              silent_print_level
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE input_constants,                 ONLY: &
        gfn1xtb, gfn2xtb, ipea1xtb, smear_energy_window, smear_fermi_dirac, smear_gaussian, &
        smear_list, smear_mp, smear_mv, tblite_cli_born_kernel_auto, tblite_cli_born_kernel_p16, &
        tblite_cli_born_kernel_still, tblite_cli_solution_state_bar1mol, &
        tblite_cli_solution_state_gsolv, tblite_cli_solution_state_reference, &
        tblite_cli_solvation_alpb, tblite_cli_solvation_cpcm, tblite_cli_solvation_gb, &
        tblite_cli_solvation_gbe, tblite_cli_solvation_gbsa, tblite_guess_ceh, tblite_guess_eeq, &
        tblite_guess_sad, tblite_mixer_damping_default, tblite_mixer_max_weight_default, &
        tblite_mixer_min_weight_default, tblite_mixer_omega0_default, &
        tblite_mixer_weight_factor_default, tblite_scc_mixer_cp2k, tblite_solver_gvd, &
        tblite_solver_gvr
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type,&
                                              total_qs_force
   USE scf_control_types,               ONLY: scf_control_type
   USE virial_types,                    ONLY: virial_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_reference_cli'

   PUBLIC :: tb_reference_cli_compare

CONTAINS

! **************************************************************************************************
!> \brief Run native tblite CLI and compare against CP2K/tblite.
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE tb_reference_cli_compare(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tb_reference_cli_compare'

      CHARACTER(LEN=16)                                  :: solvation_model_name
      CHARACTER(LEN=32) :: acc_str, charge_str, efield_x_str, efield_y_str, efield_z_str, &
         etemp_guess_val_str, etemp_str, iter_str, spin_str, spinpol_str
      CHARACTER(LEN=4*default_path_length+128) :: efield_str, etemp_guess_str, param_str, &
         post_processing_output_str, post_processing_str, restart_str, solvation_str, verbosity_str
      CHARACTER(LEN=64*default_path_length)              :: command
      CHARACTER(LEN=8)                                   :: guess, method, solver
      CHARACTER(LEN=default_path_length)                 :: file_base, gen_file, grad_file, &
                                                            json_file, log_file, &
                                                            post_processing_output_file
      INTEGER                                            :: handle, iounit, n_periodic, natom, nkp, &
                                                            reference_iterations, spin
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: do_kpoints, have_energy, have_gradient, &
                                                            have_virial, success, too_large, &
                                                            unsupported_kpoints
      REAL(KIND=dp)                                      :: cli_energy, cp_energy, ediff, etemp, &
                                                            etemp_guess, fmax, fsum, vmax, vsum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cli_gradient, cli_virial, cp_gradient
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_reference_cli_type)                       :: ref

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, cell, dft_control, energy, force, kpoints, logger, para_env, particle_set, &
               scf_control, virial, xkp)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                      dft_control=dft_control, energy=energy, force=force, &
                      para_env=para_env, particle_set=particle_set, scf_control=scf_control, &
                      virial=virial, do_kpoints=do_kpoints, kpoints=kpoints)

      ref = dft_control%qs_control%xtb_control%reference_cli
      IF (.NOT. ref%enabled) THEN
         CALL timestop(handle)
         RETURN
      END IF
      IF (.NOT. para_env%is_source()) THEN
         CALL timestop(handle)
         RETURN
      END IF

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      verbosity_str = ""
      IF (logger%iter_info%print_level == silent_print_level) THEN
         verbosity_str = " --silent"
      ELSE IF (logger%iter_info%print_level == high_print_level .OR. &
               logger%iter_info%print_level == debug_print_level) THEN
         verbosity_str = " --verbose"
      END IF
      IF (ref%solvation_active) THEN
         periodic = 0
         IF (ASSOCIATED(cell)) CALL get_cell(cell=cell, periodic=periodic)
         n_periodic = COUNT(periodic == 1)
         IF (n_periodic == 3) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") &
               "tblite reference CLI implicit solvation is not supported for PERIODIC XYZ."
            WRITE (UNIT=iounit, FMT="(T2,A)") &
               "Use PERIODIC NONE for molecular solvation diagnostics, or remove IMPLICIT_SOLVATION."
            CPABORT("REFERENCE_CLI implicit solvation is incompatible with PERIODIC XYZ")
         ELSE IF (n_periodic > 0) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A)") &
               "WARNING: tblite reference CLI implicit solvation with finite periodicity is diagnostic only."
            WRITE (UNIT=iounit, FMT="(T2,A,I0,A)") &
               "The generated native tblite reference geometry has ", n_periodic, &
               " periodic direction(s); continuum-solvation conventions are primarily molecular."
         END IF
      END IF
      unsupported_kpoints = .FALSE.
      IF (do_kpoints .AND. ASSOCIATED(kpoints)) THEN
         nkp = 0
         CALL get_kpoint_info(kpoint=kpoints, nkp=nkp, xkp=xkp)
         unsupported_kpoints = nkp > 1
         IF (nkp == 1 .AND. ASSOCIATED(xkp)) unsupported_kpoints = ANY(ABS(xkp(:, 1)) > 1.0E-12_dp)
      END IF
      IF (unsupported_kpoints) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "tblite reference CLI check skipped: CP2K KPOINTS are active."
         WRITE (UNIT=iounit, FMT="(T2,A)") &
            "The native tblite CLI reference path does not reproduce CP2K multi-k-point sampling."
         IF (ref%stop_on_error) THEN
            CPABORT("tblite reference CLI cannot check CP2K k-point calculations")
         END IF
         CALL timestop(handle)
         RETURN
      END IF
      IF (dft_control%qs_control%xtb_control%tblite_scc_mixer == tblite_scc_mixer_cp2k) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "WARNING: tblite reference CLI cannot reproduce XTB/SCC_MIXER CP2K."
         WRITE (UNIT=iounit, FMT="(T2,A)") &
            "The external native tblite run uses tblite's own SCC mixer; only the converged result is compared."
         IF (ref%stop_on_error) CPABORT("tblite reference CLI cannot reproduce SCC_MIXER CP2K")
      END IF
      natom = SIZE(particle_set)
      method = tb_reference_method_name(dft_control%qs_control%xtb_control%tblite_method)
      guess = tb_reference_guess_name(ref%guess)
      solver = tb_reference_solver_name(dft_control%qs_control%xtb_control%tblite_mixer_solver)
      file_base = tb_join_path(ref%work_directory, ref%prefix)
      gen_file = TRIM(file_base)//".gen"
      grad_file = TRIM(file_base)//".grad"
      json_file = TRIM(file_base)//".json"
      log_file = TRIM(file_base)//".log"
      post_processing_output_file = ""
      IF (LEN_TRIM(ref%grad_file) > 0) grad_file = ref%grad_file
      IF (LEN_TRIM(ref%json_file) > 0) json_file = ref%json_file
      IF (LEN_TRIM(ref%post_processing_output_file) > 0) THEN
         post_processing_output_file = ref%post_processing_output_file
      END IF

      WRITE (charge_str, "(I0)") dft_control%charge
      spin = MAX(0, dft_control%multiplicity - 1)
      WRITE (spin_str, "(I0)") spin
      WRITE (acc_str, "(ES16.8)") dft_control%qs_control%xtb_control%tblite_accuracy
      reference_iterations = dft_control%qs_control%xtb_control%tblite_mixer_iterations
      WRITE (iter_str, "(I0)") reference_iterations
      efield_str = ""
      IF (ref%efield_active) THEN
         WRITE (efield_x_str, "(ES16.8)") ref%efield(1)
         WRITE (efield_y_str, "(ES16.8)") ref%efield(2)
         WRITE (efield_z_str, "(ES16.8)") ref%efield(3)
         efield_str = " --efield "//TRIM(ADJUSTL(efield_x_str))//","// &
                      TRIM(ADJUSTL(efield_y_str))//","//TRIM(ADJUSTL(efield_z_str))
      END IF
      solvation_str = ""
      solvation_model_name = ""
      IF (ref%solvation_active) THEN
         SELECT CASE (ref%solvation_model)
         CASE (tblite_cli_solvation_alpb)
            solvation_model_name = "ALPB"
            solvation_str = " --alpb "
         CASE (tblite_cli_solvation_gbsa)
            solvation_model_name = "GBSA"
            solvation_str = " --gbsa "
         CASE (tblite_cli_solvation_gbe)
            solvation_model_name = "GBE"
            solvation_str = " --gbe "
         CASE (tblite_cli_solvation_gb)
            solvation_model_name = "GB"
            solvation_str = " --gb "
         CASE (tblite_cli_solvation_cpcm)
            solvation_model_name = "CPCM"
            solvation_str = " --cpcm "
         CASE DEFAULT
            CPABORT("Unknown tblite reference CLI implicit-solvation model")
         END SELECT
         solvation_str = TRIM(solvation_str)//" "//TRIM(tb_shell_quote(ref%solvation_solvent))
         SELECT CASE (ref%solvation_born_kernel)
         CASE (tblite_cli_born_kernel_auto)
         CASE (tblite_cli_born_kernel_p16)
            solvation_str = TRIM(solvation_str)//" --born-kernel p16"
         CASE (tblite_cli_born_kernel_still)
            solvation_str = TRIM(solvation_str)//" --born-kernel still"
         CASE DEFAULT
            CPABORT("Unknown tblite reference CLI Born kernel")
         END SELECT
         SELECT CASE (ref%solvation_state)
         CASE (tblite_cli_solution_state_gsolv)
         CASE (tblite_cli_solution_state_bar1mol)
            solvation_str = TRIM(solvation_str)//" --solv-state bar1mol"
         CASE (tblite_cli_solution_state_reference)
            solvation_str = TRIM(solvation_str)//" --solv-state reference"
         CASE DEFAULT
            CPABORT("Unknown tblite reference CLI solution state")
         END SELECT
      END IF
      IF (dft_control%qs_control%xtb_control%tblite_mixer_memory /= reference_iterations) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "WARNING: tblite reference CLI cannot reproduce XTB/TBLITE_MIXER/MEMORY."
         WRITE (UNIT=iounit, FMT="(T2,A,I0,A,I0,A)") &
            "The native reference run uses tblite's internal mixer memory tied to --iterations (", &
            reference_iterations, "), while CP2K uses MEMORY ", &
            dft_control%qs_control%xtb_control%tblite_mixer_memory, "."
         IF (ref%stop_on_error) THEN
            CPABORT("tblite reference CLI cannot reproduce TBLITE_MIXER/MEMORY")
         END IF
      END IF
      IF (ABS(dft_control%qs_control%xtb_control%tblite_mixer_damping - &
              tblite_mixer_damping_default) > 10.0_dp*EPSILON(1.0_dp)) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "WARNING: tblite reference CLI cannot reproduce XTB/TBLITE_MIXER/DAMPING."
         WRITE (UNIT=iounit, FMT="(T2,A,F8.4,A,F8.4,A)") &
            "The native reference run uses tblite's library default ", tblite_mixer_damping_default, &
            ", while CP2K uses DAMPING ", dft_control%qs_control%xtb_control%tblite_mixer_damping, "."
         IF (ref%stop_on_error) THEN
            CPABORT("tblite reference CLI cannot reproduce TBLITE_MIXER/DAMPING")
         END IF
      END IF
      IF (ABS(dft_control%qs_control%xtb_control%tblite_mixer_omega0 - &
              tblite_mixer_omega0_default) > 10.0_dp*EPSILON(1.0_dp)) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "WARNING: tblite reference CLI cannot reproduce XTB/TBLITE_MIXER/OMEGA0."
         WRITE (UNIT=iounit, FMT="(T2,A,ES12.4,A,ES12.4,A)") &
            "The native reference run uses tblite's library default ", tblite_mixer_omega0_default, &
            ", while CP2K uses OMEGA0 ", dft_control%qs_control%xtb_control%tblite_mixer_omega0, "."
         IF (ref%stop_on_error) THEN
            CPABORT("tblite reference CLI cannot reproduce TBLITE_MIXER/OMEGA0")
         END IF
      END IF
      IF (ABS(dft_control%qs_control%xtb_control%tblite_mixer_min_weight - &
              tblite_mixer_min_weight_default) > 10.0_dp*EPSILON(1.0_dp)) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "WARNING: tblite reference CLI cannot reproduce XTB/TBLITE_MIXER/MIN_WEIGHT."
         WRITE (UNIT=iounit, FMT="(T2,A,ES12.4,A,ES12.4,A)") &
            "The native reference run uses tblite's library default ", tblite_mixer_min_weight_default, &
            ", while CP2K uses MIN_WEIGHT ", dft_control%qs_control%xtb_control%tblite_mixer_min_weight, "."
         IF (ref%stop_on_error) THEN
            CPABORT("tblite reference CLI cannot reproduce TBLITE_MIXER/MIN_WEIGHT")
         END IF
      END IF
      IF (ABS(dft_control%qs_control%xtb_control%tblite_mixer_max_weight - &
              tblite_mixer_max_weight_default) > 10.0_dp*EPSILON(1.0_dp)) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "WARNING: tblite reference CLI cannot reproduce XTB/TBLITE_MIXER/MAX_WEIGHT."
         WRITE (UNIT=iounit, FMT="(T2,A,ES12.4,A,ES12.4,A)") &
            "The native reference run uses tblite's library default ", tblite_mixer_max_weight_default, &
            ", while CP2K uses MAX_WEIGHT ", dft_control%qs_control%xtb_control%tblite_mixer_max_weight, "."
         IF (ref%stop_on_error) THEN
            CPABORT("tblite reference CLI cannot reproduce TBLITE_MIXER/MAX_WEIGHT")
         END IF
      END IF
      IF (ABS(dft_control%qs_control%xtb_control%tblite_mixer_weight_factor - &
              tblite_mixer_weight_factor_default) > 10.0_dp*EPSILON(1.0_dp)) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "WARNING: tblite reference CLI cannot reproduce XTB/TBLITE_MIXER/WEIGHT_FACTOR."
         WRITE (UNIT=iounit, FMT="(T2,A,ES12.4,A,ES12.4,A)") &
            "The native reference run uses tblite's library default ", tblite_mixer_weight_factor_default, &
            ", while CP2K uses WEIGHT_FACTOR ", &
            dft_control%qs_control%xtb_control%tblite_mixer_weight_factor, "."
         IF (ref%stop_on_error) THEN
            CPABORT("tblite reference CLI cannot reproduce TBLITE_MIXER/WEIGHT_FACTOR")
         END IF
      END IF
      etemp = 300.0_dp
      IF (ASSOCIATED(scf_control)) THEN
         IF (ASSOCIATED(scf_control%smear)) THEN
            IF (scf_control%smear%do_smear) THEN
               etemp = cp_unit_from_cp2k(scf_control%smear%electronic_temperature, "K")
               IF (scf_control%smear%method /= smear_fermi_dirac) THEN
                  WRITE (UNIT=iounit, FMT="(/,T2,A,A,A)") &
                     "WARNING: tblite reference CLI cannot reproduce CP2K smearing method ", &
                     TRIM(tb_reference_smear_method_name(scf_control%smear%method)), "."
                  WRITE (UNIT=iounit, FMT="(T2,A,F12.3,A)") &
                     "The native reference run uses Fermi-Dirac electronic temperature ", etemp, " K instead."
               END IF
            END IF
         END IF
      END IF
      WRITE (etemp_str, "(ES16.8)") etemp
      etemp_guess = 0.0_dp
      etemp_guess_str = ""
      IF (ref%electronic_temperature_guess > 0.0_dp) THEN
         etemp_guess = cp_unit_from_cp2k(ref%electronic_temperature_guess, "K")
         WRITE (etemp_guess_val_str, "(ES16.8)") etemp_guess
         etemp_guess_str = " --etemp-guess "//TRIM(ADJUSTL(etemp_guess_val_str))
      END IF
      param_str = ""
      IF (LEN_TRIM(dft_control%qs_control%xtb_control%tblite_param_file) > 0) THEN
         param_str = " --param "//TRIM(tb_shell_quote(dft_control%qs_control%xtb_control%tblite_param_file))
      END IF
      spinpol_str = ""
      IF (dft_control%lsd) spinpol_str = " --spin-polarized"
      post_processing_str = ""
      IF (LEN_TRIM(ref%post_processing) > 0) THEN
         post_processing_str = " --post-processing "//TRIM(tb_shell_quote(ref%post_processing))
      END IF
      post_processing_output_str = ""
      IF (LEN_TRIM(post_processing_output_file) > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "WARNING: tblite reference CLI POST_PROCESSING_OUTPUT was requested explicitly."
         WRITE (UNIT=iounit, FMT="(T2,A)") &
            "Some tblite 0.5.0 command-line builds document --post-processing-output but do not parse it."
         post_processing_output_str = " --post-processing-output "// &
                                      TRIM(tb_shell_quote(post_processing_output_file))
      END IF
      restart_str = " --no-restart"
      IF (LEN_TRIM(ref%restart_file) > 0) THEN
         restart_str = " --restart "//TRIM(tb_shell_quote(ref%restart_file))
      END IF

      CALL tb_write_reference_gen(qs_env, TRIM(gen_file))

      command = TRIM(tb_shell_quote(ref%program_name))//" run --method "//TRIM(method)// &
                TRIM(param_str)// &
                TRIM(spinpol_str)// &
                " --charge "//TRIM(ADJUSTL(charge_str))// &
                " --spin "//TRIM(ADJUSTL(spin_str))// &
                " --acc "//TRIM(ADJUSTL(acc_str))// &
                " --guess "//TRIM(guess)// &
                " --solver "//TRIM(solver)// &
                " --iterations "//TRIM(ADJUSTL(iter_str))// &
                " --etemp "//TRIM(ADJUSTL(etemp_str))// &
                TRIM(etemp_guess_str)// &
                TRIM(efield_str)// &
                TRIM(solvation_str)// &
                TRIM(post_processing_str)// &
                TRIM(post_processing_output_str)// &
                TRIM(restart_str)// &
                TRIM(verbosity_str)//" --input "//TRIM(tb_shell_quote(ref%input_format))// &
                " --grad "//TRIM(tb_shell_quote(grad_file))// &
                " --json "//TRIM(tb_shell_quote(json_file))//" "//TRIM(tb_shell_quote(gen_file))// &
                " > "//TRIM(tb_shell_quote(log_file))//" 2>&1"

      CALL tb_reference_cli_execute(ref, "run", command, log_file, iounit, success)
      IF (.NOT. success) THEN
         ! Retain the geometry and diagnostic output on failure, even without KEEP_FILES.
         CALL timestop(handle)
         RETURN
      END IF

      ALLOCATE (cli_gradient(3, natom), cli_virial(3, 3))
      CALL tb_read_reference_grad(TRIM(grad_file), natom, cli_energy, cli_gradient, cli_virial, &
                                  have_energy, have_gradient, have_virial)

      WRITE (UNIT=iounit, FMT="(/,T2,A)") "tblite reference CLI check"
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Executable: ", TRIM(ref%program_name)
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Method:     ", TRIM(method)
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Guess:      ", TRIM(guess)
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Solver:     ", TRIM(solver)
      WRITE (UNIT=iounit, FMT="(T2,A,L1)") "Spin-pol.:  ", dft_control%lsd
      IF (ref%efield_active) THEN
         WRITE (UNIT=iounit, FMT="(T2,A,3ES16.8,A)") "Efield:     ", ref%efield, " V/Angstrom"
      END IF
      IF (ref%solvation_active) THEN
         WRITE (UNIT=iounit, FMT="(T2,A,A,1X,A)") "Solvation:  ", TRIM(solvation_model_name), &
            TRIM(ref%solvation_solvent)
      END IF
      IF (LEN_TRIM(ref%post_processing) > 0) THEN
         WRITE (UNIT=iounit, FMT="(T2,A,A)") "Post proc.: ", TRIM(ref%post_processing)
      END IF
      IF (LEN_TRIM(post_processing_output_file) > 0) THEN
         WRITE (UNIT=iounit, FMT="(T2,A,A)") "PP output:  ", TRIM(post_processing_output_file)
      END IF
      IF (ref%electronic_temperature_guess > 0.0_dp) THEN
         WRITE (UNIT=iounit, FMT="(T2,A,F12.3,A)") "Guess etemp:", etemp_guess, " K"
      END IF
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Grad file:  ", TRIM(grad_file)
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "JSON file:  ", TRIM(json_file)
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Log file:   ", TRIM(log_file)

      too_large = .FALSE.
      IF (ref%check_energy) THEN
         IF (have_energy) THEN
            cp_energy = energy%total
            ediff = ABS(cp_energy - cli_energy)
            WRITE (UNIT=iounit, FMT="(T2,A,3ES22.12)") &
               "Energy CP2K/CLI/absdiff:", cp_energy, cli_energy, ediff
            too_large = too_large .OR. ediff > ref%error_limit
         ELSE
            WRITE (UNIT=iounit, FMT="(T2,A)") "Energy check skipped: no CLI energy found."
         END IF
      END IF

      IF (ref%check_forces) THEN
         IF (have_gradient .AND. ASSOCIATED(force)) THEN
            ALLOCATE (cp_gradient(3, natom))
            CALL total_qs_force(cp_gradient, force, atomic_kind_set)
            fsum = SUM(ABS(cp_gradient - cli_gradient))
            fmax = MAXVAL(ABS(cp_gradient - cli_gradient))
            WRITE (UNIT=iounit, FMT="(T2,A,2ES22.12)") "Gradient diff sum/max:", fsum, fmax
            too_large = too_large .OR. fmax > ref%error_limit
            DEALLOCATE (cp_gradient)
         ELSE
            WRITE (UNIT=iounit, FMT="(T2,A)") "Gradient check skipped: no CLI gradient or CP2K force found."
         END IF
      END IF

      IF (ref%check_virial) THEN
         IF (have_virial .AND. ASSOCIATED(virial)) THEN
            ! Native tblite prints the positive cell derivative; CP2K stores the PV virial
            ! with the opposite sign.
            vsum = SUM(ABS(-virial%pv_virial - cli_virial))
            vmax = MAXVAL(ABS(-virial%pv_virial - cli_virial))
            WRITE (UNIT=iounit, FMT="(T2,A,2ES22.12)") "Virial diff sum/max:", vsum, vmax
            too_large = too_large .OR. vmax > ref%error_limit
         ELSE
            WRITE (UNIT=iounit, FMT="(T2,A)") "Virial check skipped: no CLI virial or CP2K virial found."
         END IF
      END IF

      IF (too_large) THEN
         WRITE (UNIT=iounit, FMT="(T2,A,ES12.4)") &
            "tblite reference CLI deviation exceeded ERROR_LIMIT = ", ref%error_limit
         IF (ref%stop_on_error) CPABORT("tblite reference CLI deviation exceeded ERROR_LIMIT")
      END IF

      CALL tb_reference_cli_aux_commands(ref, dft_control, gen_file, file_base, verbosity_str, iounit, success)

      IF (success .AND. .NOT. too_large) THEN
         CALL tb_reference_cleanup(ref, gen_file, grad_file, json_file, log_file, post_processing_output_file)
      END IF
      DEALLOCATE (cli_gradient, cli_virial)

      CALL timestop(handle)

   END SUBROUTINE tb_reference_cli_compare

! **************************************************************************************************
!> \brief Map CP2K tblite method id to native tblite CLI method name.
!> \param method_id ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_reference_method_name(method_id) RESULT(method)
      INTEGER, INTENT(IN)                                :: method_id
      CHARACTER(LEN=8)                                   :: method

      SELECT CASE (method_id)
      CASE (gfn1xtb)
         method = "gfn1"
      CASE (gfn2xtb)
         method = "gfn2"
      CASE (ipea1xtb)
         method = "ipea1"
      CASE DEFAULT
         CPABORT("Unknown tblite reference CLI method")
      END SELECT

   END FUNCTION tb_reference_method_name

! **************************************************************************************************
!> \brief Map CP2K tblite reference CLI guess id to native tblite CLI guess name.
!> \param guess_id ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_reference_guess_name(guess_id) RESULT(guess)
      INTEGER, INTENT(IN)                                :: guess_id
      CHARACTER(LEN=8)                                   :: guess

      SELECT CASE (guess_id)
      CASE (tblite_guess_sad)
         guess = "sad"
      CASE (tblite_guess_eeq)
         guess = "eeq"
      CASE (tblite_guess_ceh)
         guess = "ceh"
      CASE DEFAULT
         CPABORT("Unknown tblite reference CLI guess")
      END SELECT

   END FUNCTION tb_reference_guess_name

! **************************************************************************************************
!> \brief Map CP2K tblite solver id to native tblite CLI solver name.
!> \param solver_id ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_reference_solver_name(solver_id) RESULT(solver)
      INTEGER, INTENT(IN)                                :: solver_id
      CHARACTER(LEN=8)                                   :: solver

      SELECT CASE (solver_id)
      CASE (tblite_solver_gvd)
         solver = "gvd"
      CASE (tblite_solver_gvr)
         solver = "gvr"
      CASE DEFAULT
         CPABORT("Unknown tblite reference CLI solver")
      END SELECT

   END FUNCTION tb_reference_solver_name

! **************************************************************************************************
!> \brief Map CP2K smearing method id to a diagnostic label.
!> \param method_id ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_reference_smear_method_name(method_id) RESULT(method)
      INTEGER, INTENT(IN)                                :: method_id
      CHARACTER(LEN=24)                                  :: method

      SELECT CASE (method_id)
      CASE (smear_fermi_dirac)
         method = "FERMI_DIRAC"
      CASE (smear_energy_window)
         method = "ENERGY_WINDOW"
      CASE (smear_list)
         method = "LIST"
      CASE (smear_gaussian)
         method = "GAUSSIAN"
      CASE (smear_mp)
         method = "METHFESSEL_PAXTON"
      CASE (smear_mv)
         method = "MARZARI_VANDERBILT"
      CASE DEFAULT
         method = "UNKNOWN"
      END SELECT

   END FUNCTION tb_reference_smear_method_name

! **************************************************************************************************
!> \brief Run optional native tblite auxiliary subcommands.
!> \param ref ...
!> \param dft_control ...
!> \param gen_file ...
!> \param file_base ...
!> \param verbosity_str ...
!> \param iounit ...
!> \param success Whether all requested commands completed successfully.
! **************************************************************************************************
   SUBROUTINE tb_reference_cli_aux_commands(ref, dft_control, gen_file, file_base, verbosity_str, iounit, success)

      TYPE(xtb_reference_cli_type), INTENT(IN)           :: ref
      TYPE(dft_control_type), INTENT(IN)                 :: dft_control
      CHARACTER(LEN=*), INTENT(IN)                       :: gen_file, file_base, verbosity_str
      INTEGER, INTENT(IN)                                :: iounit
      LOGICAL, INTENT(OUT)                               :: success

      CHARACTER(LEN=32)                                  :: charge_str, efield_x_str, efield_y_str, &
                                                            efield_z_str, etemp_guess_val_str, &
                                                            spin_str
      CHARACTER(LEN=4*default_path_length+128)           :: copy_str, dry_run_str, efield_str, &
                                                            etemp_guess_str, grad_str, json_str, &
                                                            method_str, output_str
      CHARACTER(LEN=64*default_path_length)              :: command
      CHARACTER(LEN=default_path_length)                 :: guess_input, log_file
      INTEGER                                            :: spin
      LOGICAL                                            :: command_success
      REAL(KIND=dp)                                      :: etemp_guess

      success = .TRUE.
      WRITE (charge_str, "(I0)") dft_control%charge
      spin = MAX(0, dft_control%multiplicity - 1)
      WRITE (spin_str, "(I0)") spin

      IF (ref%guess_cli%enabled) THEN
         guess_input = ref%guess_cli%input_file
         IF (LEN_TRIM(guess_input) == 0) guess_input = gen_file
         etemp_guess_str = ""
         IF (ref%guess_cli%electronic_temperature_guess > 0.0_dp) THEN
            etemp_guess = cp_unit_from_cp2k(ref%guess_cli%electronic_temperature_guess, "K")
            WRITE (etemp_guess_val_str, "(ES16.8)") etemp_guess
            etemp_guess_str = " --etemp-guess "//TRIM(ADJUSTL(etemp_guess_val_str))
         END IF
         efield_str = ""
         IF (ref%guess_cli%efield_active) THEN
            WRITE (efield_x_str, "(ES16.8)") ref%guess_cli%efield(1)
            WRITE (efield_y_str, "(ES16.8)") ref%guess_cli%efield(2)
            WRITE (efield_z_str, "(ES16.8)") ref%guess_cli%efield(3)
            efield_str = " --efield "//TRIM(ADJUSTL(efield_x_str))//","// &
                         TRIM(ADJUSTL(efield_y_str))//","//TRIM(ADJUSTL(efield_z_str))
         END IF
         grad_str = ""
         IF (ref%guess_cli%grad) grad_str = " --grad"
         json_str = ""
         IF (LEN_TRIM(ref%guess_cli%json_file) > 0) THEN
            json_str = " --json "//TRIM(tb_shell_quote(ref%guess_cli%json_file))
         END IF
         log_file = TRIM(file_base)//".guess.log"
         command = TRIM(tb_shell_quote(ref%program_name))// &
                   " guess --charge "//TRIM(ADJUSTL(charge_str))// &
                   " --spin "//TRIM(ADJUSTL(spin_str))// &
                   " --method "//TRIM(tb_reference_guess_name(ref%guess_cli%method))// &
                   " --solver "//TRIM(tb_reference_solver_name(ref%guess_cli%solver))// &
                   TRIM(etemp_guess_str)// &
                   TRIM(efield_str)// &
                   TRIM(grad_str)// &
                   TRIM(json_str)// &
                   TRIM(verbosity_str)//" --input "//TRIM(tb_shell_quote(ref%guess_cli%input_format))// &
                   " "//TRIM(tb_shell_quote(guess_input))// &
                   " > "//TRIM(tb_shell_quote(log_file))//" 2>&1"
         CALL tb_reference_cli_execute(ref, "guess", command, log_file, iounit, command_success)
         success = success .AND. command_success
         IF (command_success .AND. .NOT. ref%keep_files) CALL tb_delete_file(log_file)
         IF (command_success .AND. .NOT. ref%keep_files .AND. LEN_TRIM(ref%guess_cli%json_file) > 0) THEN
            CALL tb_delete_file(ref%guess_cli%json_file)
         END IF
      END IF

      IF (ref%param_cli%enabled) THEN
         method_str = ""
         IF (ref%param_cli%method_explicit .OR. LEN_TRIM(ref%param_cli%input_file) == 0) THEN
            method_str = " --method "// &
                         TRIM(tb_reference_method_name(MERGE(ref%param_cli%method, &
                                                             dft_control%qs_control%xtb_control%tblite_method, &
                                                             ref%param_cli%method_explicit)))
         END IF
         output_str = ""
         IF (LEN_TRIM(ref%param_cli%output_file) > 0) THEN
            output_str = " --output "//TRIM(tb_shell_quote(ref%param_cli%output_file))
         END IF
         log_file = TRIM(file_base)//".param.log"
         command = TRIM(tb_shell_quote(ref%program_name))//" param"// &
                   TRIM(method_str)// &
                   TRIM(output_str)
         IF (LEN_TRIM(ref%param_cli%input_file) > 0) THEN
            command = TRIM(command)//" "//TRIM(tb_shell_quote(ref%param_cli%input_file))
         END IF
         command = TRIM(command)//" > "//TRIM(tb_shell_quote(log_file))//" 2>&1"
         CALL tb_reference_cli_execute(ref, "param", command, log_file, iounit, command_success)
         success = success .AND. command_success
         IF (command_success .AND. .NOT. ref%keep_files) CALL tb_delete_file(log_file)
         IF (command_success .AND. .NOT. ref%keep_files .AND. LEN_TRIM(ref%param_cli%output_file) > 0) THEN
            CALL tb_delete_file(ref%param_cli%output_file)
         END IF
      END IF

      IF (ref%fit_cli%enabled) THEN
         dry_run_str = ""
         IF (ref%fit_cli%dry_run) dry_run_str = " --dry-run"
         copy_str = ""
         IF (LEN_TRIM(ref%fit_cli%copy_file) > 0) THEN
            copy_str = " --copy "//TRIM(tb_shell_quote(ref%fit_cli%copy_file))
         END IF
         log_file = TRIM(file_base)//".fit.log"
         command = TRIM(tb_shell_quote(ref%program_name))//" fit"// &
                   TRIM(dry_run_str)// &
                   TRIM(copy_str)// &
                   TRIM(verbosity_str)//" "//TRIM(tb_shell_quote(ref%fit_cli%param_file))// &
                   " "//TRIM(tb_shell_quote(ref%fit_cli%input_file))// &
                   " > "//TRIM(tb_shell_quote(log_file))//" 2>&1"
         CALL tb_reference_cli_execute(ref, "fit", command, log_file, iounit, command_success)
         success = success .AND. command_success
         IF (command_success .AND. .NOT. ref%keep_files) CALL tb_delete_file(log_file)
         IF (command_success .AND. .NOT. ref%keep_files .AND. LEN_TRIM(ref%fit_cli%copy_file) > 0) THEN
            CALL tb_delete_file(ref%fit_cli%copy_file)
         END IF
      END IF

      IF (ref%tagdiff_cli%enabled) THEN
         method_str = ""
         IF (ref%tagdiff_cli%fit) method_str = " --fit"
         log_file = TRIM(file_base)//".tagdiff.log"
         command = TRIM(tb_shell_quote(ref%program_name))//" tagdiff"// &
                   TRIM(method_str)//" "//TRIM(tb_shell_quote(ref%tagdiff_cli%actual_file))// &
                   " "//TRIM(tb_shell_quote(ref%tagdiff_cli%reference_file))// &
                   " > "//TRIM(tb_shell_quote(log_file))//" 2>&1"
         CALL tb_reference_cli_execute(ref, "tagdiff", command, log_file, iounit, command_success)
         success = success .AND. command_success
         IF (command_success .AND. .NOT. ref%keep_files) CALL tb_delete_file(log_file)
      END IF

   END SUBROUTINE tb_reference_cli_aux_commands

! **************************************************************************************************
!> \brief Execute one native tblite REFERENCE_CLI command.
!> \param ref ...
!> \param label ...
!> \param command ...
!> \param log_file ...
!> \param iounit ...
!> \param success Whether all requested commands completed successfully.
!> \note Keep the external-process exception confined to this optional diagnostic path.
! **************************************************************************************************
   SUBROUTINE tb_reference_cli_execute(ref, label, command, log_file, iounit, success)

      TYPE(xtb_reference_cli_type), INTENT(IN)           :: ref
      CHARACTER(LEN=*), INTENT(IN)                       :: label, command, log_file
      INTEGER, INTENT(IN)                                :: iounit
      LOGICAL, INTENT(OUT)                               :: success

      CHARACTER(LEN=default_path_length)                 :: cmdmsg
      INTEGER                                            :: cmdstat, exitstat

      cmdmsg = ""
      cmdstat = 0
      exitstat = 0
      CALL EXECUTE_COMMAND_LINE(TRIM(command), EXITSTAT=exitstat, CMDSTAT=cmdstat, CMDMSG=cmdmsg)
      success = cmdstat == 0 .AND. exitstat == 0
      IF (.NOT. success) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") "tblite reference CLI check failed to run."
         WRITE (UNIT=iounit, FMT="(T2,A,A)") "Command: ", TRIM(command)
         WRITE (UNIT=iounit, FMT="(T2,A,I0,T32,A,I0)") "cmdstat:", cmdstat, "exitstat:", exitstat
         IF (cmdstat /= 0) WRITE (UNIT=iounit, FMT="(T2,A,A)") "Execution error: ", TRIM(cmdmsg)
         WRITE (UNIT=iounit, FMT="(T2,A,A)") "Diagnostic files retained; log: ", TRIM(log_file)
         IF (ref%stop_on_error) CPABORT("tblite reference CLI command failed")
      ELSE
         WRITE (UNIT=iounit, FMT="(/,T2,A,A)") "tblite reference CLI command completed: ", &
            TRIM(label)
         WRITE (UNIT=iounit, FMT="(T2,A,A)") "Log file: ", TRIM(log_file)
      END IF

   END SUBROUTINE tb_reference_cli_execute

! **************************************************************************************************
!> \brief Join directory and filename.
!> \param directory ...
!> \param filename ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_join_path(directory, filename) RESULT(path)
      CHARACTER(LEN=*), INTENT(IN)                       :: directory, filename
      CHARACTER(LEN=default_path_length)                 :: path

      IF (LEN_TRIM(directory) == 0 .OR. TRIM(directory) == ".") THEN
         path = TRIM(filename)
      ELSE IF (directory(LEN_TRIM(directory):LEN_TRIM(directory)) == "/") THEN
         path = TRIM(directory)//TRIM(filename)
      ELSE
         path = TRIM(directory)//"/"//TRIM(filename)
      END IF

   END FUNCTION tb_join_path

! **************************************************************************************************
!> \brief Shell-quote a filename or executable path.
!> \param text ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_shell_quote(text) RESULT(quoted)
      CHARACTER(LEN=*), INTENT(IN)                       :: text
      CHARACTER(LEN=4*LEN_TRIM(text)+2)                  :: quoted

      INTEGER                                            :: i, pos

      quoted = "'"
      pos = 2
      DO i = 1, LEN_TRIM(text)
         IF (text(i:i) == "'") THEN
            ! End the quoted part, escape the apostrophe, then resume quoting.
            ! ACHAR avoids dependence on the compiler's backslash-literal mode.
            quoted(pos:pos + 3) = "'"//ACHAR(92)//"''"
            pos = pos + 4
         ELSE
            quoted(pos:pos) = text(i:i)
            pos = pos + 1
         END IF
      END DO
      quoted(pos:pos) = "'"

   END FUNCTION tb_shell_quote

! **************************************************************************************************
!> \brief Write current CP2K geometry as DFTB+ gen format for native tblite.
!> \param qs_env ...
!> \param filename ...
! **************************************************************************************************
   SUBROUTINE tb_write_reference_gen(qs_env, filename)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(LEN=*), INTENT(IN)                       :: filename

      CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:)        :: symbols, unique_symbols
      INTEGER                                            :: iatom, ikind, ios, natom, nuniq, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: species
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: to_angstrom
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env=qs_env, cell=cell, particle_set=particle_set)

      natom = SIZE(particle_set)
      to_angstrom = cp_unit_from_cp2k(1.0_dp, "angstrom")
      ALLOCATE (symbols(natom), unique_symbols(natom), species(natom))
      nuniq = 0
      DO iatom = 1, natom
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, element_symbol=symbols(iatom))
         found = .FALSE.
         DO ikind = 1, nuniq
            IF (TRIM(unique_symbols(ikind)) == TRIM(symbols(iatom))) THEN
               found = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT. found) THEN
            nuniq = nuniq + 1
            unique_symbols(nuniq) = symbols(iatom)
            ikind = nuniq
         END IF
         species(iatom) = ikind
      END DO

      OPEN (NEWUNIT=unit_nr, FILE=TRIM(filename), STATUS="REPLACE", ACTION="WRITE", &
            FORM="FORMATTED", IOSTAT=ios)
      IF (ios /= 0) CPABORT("Could not open tblite reference CLI geometry file")

      CALL get_cell(cell=cell, periodic=periodic)
      IF (ANY(periodic == 1)) THEN
         WRITE (UNIT=unit_nr, FMT="(I0,1X,A)") natom, "S"
      ELSE
         WRITE (UNIT=unit_nr, FMT="(I0,1X,A)") natom, "C"
      END IF
      WRITE (UNIT=unit_nr, FMT="(*(A,1X))") (TRIM(unique_symbols(ikind)), ikind=1, nuniq)
      DO iatom = 1, natom
         WRITE (UNIT=unit_nr, FMT="(I0,1X,I0,3(1X,ES24.16))") &
            iatom, species(iatom), particle_set(iatom)%r(:)*to_angstrom
      END DO
      IF (ANY(periodic == 1)) THEN
         WRITE (UNIT=unit_nr, FMT="(3(1X,ES24.16))") 0.0_dp, 0.0_dp, 0.0_dp
         DO ikind = 1, 3
            WRITE (UNIT=unit_nr, FMT="(3(1X,ES24.16))") cell%hmat(:, ikind)*to_angstrom
         END DO
      END IF
      CLOSE (unit_nr)

      DEALLOCATE (symbols, unique_symbols, species)

   END SUBROUTINE tb_write_reference_gen

! **************************************************************************************************
!> \brief Read native tblite gradient file.
!> \param filename ...
!> \param natom ...
!> \param energy ...
!> \param gradient ...
!> \param virial ...
!> \param have_energy ...
!> \param have_gradient ...
!> \param have_virial ...
! **************************************************************************************************
   SUBROUTINE tb_read_reference_grad(filename, natom, energy, gradient, virial, &
                                     have_energy, have_gradient, have_virial)

      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      INTEGER, INTENT(IN)                                :: natom
      REAL(KIND=dp), INTENT(OUT)                         :: energy
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: gradient, virial
      LOGICAL, INTENT(OUT)                               :: have_energy, have_gradient, have_virial

      CHARACTER(LEN=1024)                                :: line
      INTEGER                                            :: ios, nread, unit_nr
      LOGICAL                                            :: exists
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: values

      have_energy = .FALSE.
      have_gradient = .FALSE.
      have_virial = .FALSE.
      energy = 0.0_dp
      gradient = 0.0_dp
      virial = 0.0_dp

      INQUIRE (FILE=TRIM(filename), EXIST=exists)
      IF (.NOT. exists) RETURN

      OPEN (NEWUNIT=unit_nr, FILE=TRIM(filename), STATUS="OLD", ACTION="READ", &
            FORM="FORMATTED", IOSTAT=ios)
      IF (ios /= 0) RETURN

      DO
         READ (UNIT=unit_nr, FMT="(A)", IOSTAT=ios) line
         IF (ios /= 0) EXIT
         IF (INDEX(line, "energy             :real:0:") > 0) THEN
            READ (UNIT=unit_nr, FMT="(A)", IOSTAT=ios) line
            IF (ios == 0) THEN
               READ (line, *, IOSTAT=ios) energy
               have_energy = ios == 0
            END IF
         ELSE IF (INDEX(line, "gradient           :real:2:3,") > 0) THEN
            ALLOCATE (values(3*natom))
            CALL tb_read_real_values(unit_nr, values, nread)
            IF (nread == 3*natom) THEN
               CALL tb_values_to_matrix(values, gradient)
               have_gradient = .TRUE.
            END IF
            DEALLOCATE (values)
         ELSE IF (INDEX(line, "virial             :real:2:3,3") > 0) THEN
            ALLOCATE (values(9))
            CALL tb_read_real_values(unit_nr, values, nread)
            IF (nread == 9) THEN
               CALL tb_values_to_matrix(values, virial)
               have_virial = .TRUE.
            END IF
            DEALLOCATE (values)
         END IF
      END DO
      CLOSE (unit_nr)

   END SUBROUTINE tb_read_reference_grad

! **************************************************************************************************
!> \brief Read a fixed number of real values from following lines.
!> \param unit_nr ...
!> \param values ...
!> \param nread ...
! **************************************************************************************************
   SUBROUTINE tb_read_real_values(unit_nr, values, nread)

      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: values
      INTEGER, INTENT(OUT)                               :: nread

      CHARACTER(LEN=1024)                                :: line
      INTEGER                                            :: ios

      nread = 0
      DO WHILE (nread < SIZE(values))
         READ (UNIT=unit_nr, FMT="(A)", IOSTAT=ios) line
         IF (ios /= 0) EXIT
         CALL tb_parse_real_line(line, values, nread)
      END DO

   END SUBROUTINE tb_read_real_values

! **************************************************************************************************
!> \brief Parse real values from one text line.
!> \param line ...
!> \param values ...
!> \param nread ...
! **************************************************************************************************
   SUBROUTINE tb_parse_real_line(line, values, nread)

      CHARACTER(LEN=*), INTENT(IN)                       :: line
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: values
      INTEGER, INTENT(INOUT)                             :: nread

      CHARACTER(LEN=128)                                 :: token
      INTEGER                                            :: first, ios, last, pos

      pos = 1
      DO WHILE (pos <= LEN_TRIM(line) .AND. nread < SIZE(values))
         DO WHILE (pos <= LEN_TRIM(line) .AND. INDEX(" ,[]", line(pos:pos)) > 0)
            pos = pos + 1
         END DO
         IF (pos > LEN_TRIM(line)) EXIT
         first = pos
         DO WHILE (pos <= LEN_TRIM(line) .AND. INDEX(" ,[]", line(pos:pos)) == 0)
            pos = pos + 1
         END DO
         last = pos - 1
         token = line(first:last)
         READ (token, *, IOSTAT=ios) values(nread + 1)
         IF (ios == 0) nread = nread + 1
      END DO

   END SUBROUTINE tb_parse_real_line

! **************************************************************************************************
!> \brief Convert flat native tblite values to CP2K atom-major matrix layout.
!> \param values ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE tb_values_to_matrix(values, matrix)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: values
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: matrix

      INTEGER                                            :: i, j, n

      n = 0
      DO j = 1, SIZE(matrix, 2)
         DO i = 1, SIZE(matrix, 1)
            n = n + 1
            matrix(i, j) = values(n)
         END DO
      END DO

   END SUBROUTINE tb_values_to_matrix

! **************************************************************************************************
!> \brief Delete temporary files unless requested otherwise.
!> \param ref ...
!> \param gen_file ...
!> \param grad_file ...
!> \param json_file ...
!> \param log_file ...
!> \param post_processing_output_file ...
! **************************************************************************************************
   SUBROUTINE tb_reference_cleanup(ref, gen_file, grad_file, json_file, log_file, post_processing_output_file)

      TYPE(xtb_reference_cli_type), INTENT(IN)           :: ref
      CHARACTER(LEN=*), INTENT(IN)                       :: gen_file, grad_file, json_file, &
                                                            log_file, post_processing_output_file

      IF (ref%keep_files) RETURN
      CALL tb_delete_file(gen_file)
      CALL tb_delete_file(grad_file)
      CALL tb_delete_file(json_file)
      CALL tb_delete_file(log_file)
      IF (LEN_TRIM(post_processing_output_file) > 0) THEN
         CALL tb_delete_file(post_processing_output_file)
      END IF

   END SUBROUTINE tb_reference_cleanup

! **************************************************************************************************
!> \brief Delete a file if it exists.
!> \param filename ...
! **************************************************************************************************
   SUBROUTINE tb_delete_file(filename)

      CHARACTER(LEN=*), INTENT(IN)                       :: filename

      INTEGER                                            :: ios, unit_nr
      LOGICAL                                            :: exists

      INQUIRE (FILE=TRIM(filename), EXIST=exists)
      IF (.NOT. exists) RETURN
      OPEN (NEWUNIT=unit_nr, FILE=TRIM(filename), STATUS="OLD", IOSTAT=ios)
      IF (ios == 0) CLOSE (unit_nr, STATUS="DELETE")

   END SUBROUTINE tb_delete_file

END MODULE tblite_reference_cli
