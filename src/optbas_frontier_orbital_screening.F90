!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Automatic screening of frontier-orbital basis optimizations
!> \par History
!>      28.08.2026 created [Jan Wilhelm]
!> \author Jan Wilhelm
! **************************************************************************************************
MODULE optbas_frontier_orbital_screening
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_files,                        ONLY: close_file,&
                                              discover_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: debug_print_level,&
                                              silent_print_level
   USE global_types,                    ONLY: global_environment_type
   USE input_constants,                 ONLY: atomic_guess,&
                                              do_qs
   USE input_section_types,             ONLY: section_type,&
                                              section_vals_duplicate,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_release,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE machine,                         ONLY: m_getcwd
   USE message_passing,                 ONLY: mp_para_env_type
   USE optbas_fenv_manipulation,        ONLY: update_basis_set
   USE optbas_frontier_orbitals_utils,  ONLY: frontier_orbitals_reference_init,&
                                              frontier_orbitals_reference_release,&
                                              frontier_orbitals_reference_type,&
                                              print_frontier_orbital_basis_sets
   USE optimize_basis_types,            ONLY: basis_optimization_type,&
                                              deallocate_basis_optimization_type,&
                                              method_mo_fit_occ_virtual
   USE optimize_basis_utils,            ONLY: optimize_basis_init_read_input,&
                                              write_basis
   USE physcon,                         ONLY: evolt
   USE qs_condnum,                      ONLY: overlap_condnum
   USE qs_energy,                       ONLY: qs_energies
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment,                  ONLY: qs_init
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_create,&
                                              qs_env_release,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE string_utilities,                ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: run_frontier_orbital_screening

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optbas_frontier_orbital_screening'

   TYPE frontier_orbital_screening_candidate_type
      REAL(KIND=dp) :: condition_number = 0.0_dp
      REAL(KIND=dp) :: gap = 0.0_dp
      REAL(KIND=dp) :: energy = 0.0_dp
      TYPE(basis_optimization_type) :: opt_bas
   END TYPE frontier_orbital_screening_candidate_type

   ABSTRACT INTERFACE
! **************************************************************************************************
!> \brief ...
!> \param opt_bas ...
!> \param input_declaration ...
!> \param para_env ...
!> \param globenv ...
!> \param training_input ...
!> \param reference ...
! **************************************************************************************************
      SUBROUTINE basis_optimization_driver(opt_bas, input_declaration, para_env, globenv, training_input, reference)
         IMPORT :: basis_optimization_type, frontier_orbitals_reference_type, global_environment_type, &
            mp_para_env_type, &
            section_type, section_vals_type
         TYPE(basis_optimization_type)                   :: opt_bas
         TYPE(section_type), POINTER                     :: input_declaration
         TYPE(mp_para_env_type), POINTER                 :: para_env
         TYPE(global_environment_type), POINTER          :: globenv
         TYPE(section_vals_type), OPTIONAL, POINTER      :: training_input
         TYPE(frontier_orbitals_reference_type), &
            DIMENSION(:), INTENT(IN), OPTIONAL, TARGET   :: reference
      END SUBROUTINE basis_optimization_driver
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Run serial frontier-orbital optimizations and evaluate their HOMO-LUMO gaps with new SCFs.
!> \param input_declaration ...
!> \param section ...
!> \param para_env ...
!> \param globenv ...
!> \param frontier_screening_section ...
!> \param optimize_basis_driver ...
! **************************************************************************************************
   SUBROUTINE run_frontier_orbital_screening(input_declaration, section, para_env, globenv, frontier_screening_section, &
                                             optimize_basis_driver)
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: frontier_screening_section

      PROCEDURE(basis_optimization_driver)                :: optimize_basis_driver

      INTEGER, PARAMETER :: number_parameter_sets = 5
      REAL(KIND=dp), DIMENSION(6, number_parameter_sets), PARAMETER :: &
         pre_tabulated_parameters = RESHAPE( &
         [10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 0.001_dp, 0.001_dp, &
          30.0_dp, 30.0_dp, 30.0_dp, 10.0_dp, 0.001_dp, 0.001_dp, &
          10.0_dp, 30.0_dp, 30.0_dp, 10.0_dp, 0.001_dp, 0.001_dp, &
          10.0_dp, 30.0_dp, 30.0_dp, 3.0_dp, 0.001_dp, 0.001_dp, &
          10.0_dp, 30.0_dp, 30.0_dp, 10.0_dp, 0.001_dp, 0.01_dp], &
         [6, number_parameter_sets])

      CHARACTER(LEN=default_path_length)                 :: final_basis_file, main_dir
      INTEGER                                            :: i, n_candidates, unit_nr
      REAL(KIND=dp) :: initial_condition_number, initial_energy, initial_gap, &
                       reference_condition_number, reference_energy, reference_gap
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(frontier_orbitals_reference_type), &
         ALLOCATABLE, DIMENSION(:), TARGET               :: reference
      TYPE(frontier_orbital_screening_candidate_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: candidates

! Rows are occupied, virtual, empty, gap, condition, and coefficient weights.

      CALL check_frontier_orbital_screening_input(section)
      CALL section_vals_val_get(frontier_screening_section, "NUMBER_OF_OPTIMIZATIONS", i_val=n_candidates)
      IF (n_candidates < 1 .OR. n_candidates > number_parameter_sets) THEN
         CALL cp_abort(__LOCATION__, &
                       "NUMBER_OF_OPTIMIZATIONS must be between 1 and "// &
                       TRIM(cp_to_string(number_parameter_sets)))
      END IF

      CALL m_getcwd(main_dir)
      ALLOCATE (candidates(n_candidates), reference(1))
      logger => cp_get_default_logger()
      unit_nr = -1
      IF (para_env%is_source()) unit_nr = cp_logger_get_default_unit_nr(logger)

      CALL run_frontier_orbital_reference_scf(section, para_env, globenv, reference(1), &
                                              reference_gap, reference_energy, &
                                              reference_condition_number)

      DO i = 1, n_candidates
         CALL optimize_basis_init_read_input(candidates(i)%opt_bas, section, para_env, &
                                             quiet_output=.TRUE., embedded_training=.TRUE.)
         IF (candidates(i)%opt_bas%method /= method_mo_fit_occ_virtual) THEN
            CPABORT("Internal error while initializing FRONTIER_ORBITAL_SCREENING")
         END IF
         IF (candidates(i)%opt_bas%ntraining_sets /= 1 .OR. candidates(i)%opt_bas%ncombinations /= 1) THEN
            CALL cp_abort(__LOCATION__, &
                          "FRONTIER_ORBITAL_SCREENING currently supports exactly one training system "// &
                          "and one candidate basis assignment")
         END IF

         IF (i == 1) THEN
            final_basis_file = work_path(main_dir, candidates(i)%opt_bas%output_basis_file)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(A)') ""
               WRITE (unit_nr, '(1X,A)') REPEAT("=", 79)
               WRITE (unit_nr, '(1X,A)') "Optimize basis set to improve frontier orbitals:"
               WRITE (unit_nr, '(1X,A)') "Using pre-tabulated parameters in the optimization"
               WRITE (unit_nr, '(1X,A)') REPEAT("=", 79)
               WRITE (unit_nr, '(A)') " "
            END IF
         END IF

         candidates(i)%opt_bas%method = method_mo_fit_occ_virtual
         candidates(i)%opt_bas%use_condition_number = .TRUE.
         candidates(i)%opt_bas%write_frequency = 5000
         candidates(i)%opt_bas%fval_weight = 1.0_dp
         candidates(i)%opt_bas%occupied_weight = pre_tabulated_parameters(1, i)
         candidates(i)%opt_bas%virtual_weight = pre_tabulated_parameters(2, i)
         candidates(i)%opt_bas%empty_overlap_weight = pre_tabulated_parameters(3, i)
         candidates(i)%opt_bas%gap_weight = pre_tabulated_parameters(4, i)
         candidates(i)%opt_bas%condition_weight = pre_tabulated_parameters(5, i)
         candidates(i)%opt_bas%coefficient_weight = pre_tabulated_parameters(6, i)
         candidates(i)%opt_bas%gap_energy_scale = 1.0_dp/evolt
         candidates(i)%opt_bas%virtual_energy_cutoff = 3.0_dp/evolt
         candidates(i)%opt_bas%virtual_energy_smoothing = 0.2_dp/evolt
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(1X,A,I0,A,I0,A)') "Candidate ", i, " of ", n_candidates, &
               ": optimizing basis set ..."
         END IF
         CALL optimize_basis_driver(candidates(i)%opt_bas, input_declaration, para_env, globenv, &
                                    training_input=section, reference=reference)
      END DO

      CALL run_frontier_orbital_scf(candidates(1)%opt_bas, section, para_env, globenv, &
                                    optimized_basis=.FALSE., &
                                    gap=initial_gap, energy=initial_energy, &
                                    condition_number=initial_condition_number)

      DO i = 1, n_candidates
         CALL run_frontier_orbital_scf(candidates(i)%opt_bas, section, para_env, globenv, &
                                       optimized_basis=.TRUE., &
                                       gap=candidates(i)%gap, energy=candidates(i)%energy, &
                                       condition_number=candidates(i)%condition_number)
      END DO

      CALL write_frontier_orbital_basis_file(final_basis_file, candidates, &
                                             reference_gap, reference_energy, reference_condition_number, para_env)
      CALL print_frontier_orbital_summary(unit_nr, final_basis_file, section, candidates, &
                                          reference_gap, reference_energy, reference_condition_number, &
                                          initial_gap, initial_energy, initial_condition_number)

      DO i = 1, n_candidates
         CALL deallocate_basis_optimization_type(candidates(i)%opt_bas)
      END DO
      CALL frontier_orbitals_reference_release(reference(1))
      DEALLOCATE (candidates, reference)

   END SUBROUTINE run_frontier_orbital_screening

! **************************************************************************************************
!> \brief Check the single-file input contract of frontier-orbital screening.
!> \param section ...
! **************************************************************************************************
   SUBROUTINE check_frontier_orbital_screening_input(section)
      TYPE(section_vals_type), POINTER                   :: section

      INTEGER                                            :: method, nforce_eval, ntraining
      LOGICAL                                            :: dft_explicit, keyword_explicit
      TYPE(section_vals_type), POINTER                   :: dft_section, force_eval_section, &
                                                            optbas_section, training_section

      optbas_section => section_vals_get_subs_vals(section, "OPTIMIZE_BASIS")

      CALL reject_controlled_keyword(optbas_section, "BASIS_WORK_FILE")
      CALL reject_controlled_keyword(optbas_section, "WRITE_FREQUENCY")
      CALL reject_controlled_keyword(optbas_section, "USE_CONDITION_NUMBER")
      CALL reject_controlled_keyword(optbas_section, "RESIDUUM_WEIGHT")
      CALL reject_controlled_keyword(optbas_section, "CONDITION_WEIGHT")
      CALL reject_controlled_keyword(optbas_section, "BASIS_COMBINATIONS")
      CALL reject_controlled_keyword(optbas_section, "GROUP_PARTITION")
      CALL section_vals_val_get(optbas_section, "BASIS_TEMPLATE_FILE", explicit=keyword_explicit)
      IF (.NOT. keyword_explicit) THEN
         CPABORT("FRONTIER_ORBITAL_SCREENING requires BASIS_TEMPLATE_FILE")
      END IF

      training_section => section_vals_get_subs_vals(optbas_section, "TRAINING_FILES")
      CALL section_vals_get(training_section, n_repetition=ntraining)
      IF (ntraining /= 0) CALL controlled_parameter_error("TRAINING_FILES")

      force_eval_section => section_vals_get_subs_vals(section, "FORCE_EVAL")
      CALL section_vals_get(force_eval_section, n_repetition=nforce_eval)
      IF (nforce_eval /= 1) THEN
         CALL cp_abort(__LOCATION__, &
                       "FRONTIER_ORBITAL_SCREENING requires exactly one FORCE_EVAL in the same input")
      END IF
      CALL section_vals_val_get(force_eval_section, "METHOD", i_val=method)
      IF (method /= do_qs) THEN
         CPABORT("FRONTIER_ORBITAL_SCREENING requires a Quickstep FORCE_EVAL")
      END IF
      dft_section => section_vals_get_subs_vals(force_eval_section, "DFT")
      CALL section_vals_get(dft_section, explicit=dft_explicit)
      IF (.NOT. dft_explicit) THEN
         CALL cp_abort(__LOCATION__, &
                       "FRONTIER_ORBITAL_SCREENING requires an explicit FORCE_EVAL%DFT section")
      END IF
   END SUBROUTINE check_frontier_orbital_screening_input

! **************************************************************************************************
!> \brief Reject an OPTIMIZE_BASIS keyword controlled by frontier-orbital screening.
!> \param section ...
!> \param keyword_name ...
! **************************************************************************************************
   SUBROUTINE reject_controlled_keyword(section, keyword_name)
      TYPE(section_vals_type), POINTER                   :: section
      CHARACTER(LEN=*), INTENT(IN)                       :: keyword_name

      LOGICAL                                            :: explicit

      CALL section_vals_val_get(section, keyword_name, explicit=explicit)
      IF (explicit) CALL controlled_parameter_error("OPTIMIZE_BASIS%"//TRIM(keyword_name))

   END SUBROUTINE reject_controlled_keyword

! **************************************************************************************************
!> \brief Abort for an input parameter replaced by frontier-orbital screening.
!> \param parameter_name ...
! **************************************************************************************************
   SUBROUTINE controlled_parameter_error(parameter_name)
      CHARACTER(LEN=*), INTENT(IN)                       :: parameter_name

      CALL cp_abort(__LOCATION__, &
                    TRIM(parameter_name)//" is not allowed together with FRONTIER_ORBITAL_SCREENING "// &
                    "because it is overwritten automatically")

   END SUBROUTINE controlled_parameter_error

! **************************************************************************************************
!> \brief Run and retain the large-basis reference from the embedded FORCE_EVAL.
!> \param section ...
!> \param para_env ...
!> \param globenv ...
!> \param reference ...
!> \param gap ...
!> \param energy ...
!> \param condition_number ...
! **************************************************************************************************
   SUBROUTINE run_frontier_orbital_reference_scf(section, para_env, globenv, reference, &
                                                 gap, energy, condition_number)
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(frontier_orbitals_reference_type), &
         INTENT(OUT)                                     :: reference
      REAL(KIND=dp), INTENT(OUT)                         :: gap, energy, condition_number

      INTEGER                                            :: outer_print_level
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input_file

      NULLIFY (input_file, logger, qs_env)
      CALL section_vals_duplicate(section, input_file)
      CALL section_vals_val_set(input_file, "GLOBAL%PRINT_LEVEL", i_val=silent_print_level)
      CALL section_vals_val_set(input_file, "FORCE_EVAL%DFT%SCF%SCF_GUESS", i_val=atomic_guess)
      CALL section_vals_val_set(input_file, "GLOBAL%PROGRAM_RUN_INFO%_SECTION_PARAMETERS_", &
                                i_val=debug_print_level)
      CALL section_vals_val_set(input_file, "GLOBAL%REFERENCES%_SECTION_PARAMETERS_", &
                                i_val=debug_print_level)
      CALL section_vals_val_set(input_file, "GLOBAL%TIMINGS%_SECTION_PARAMETERS_", &
                                i_val=debug_print_level)
      CALL section_vals_val_set(input_file, &
                                "FORCE_EVAL%DFT%SCF%PRINT%RESTART%_SECTION_PARAMETERS_", &
                                i_val=debug_print_level + 1)
      CALL section_vals_val_set(input_file, &
                                "FORCE_EVAL%DFT%SCF%PRINT%RESTART_HISTORY%_SECTION_PARAMETERS_", &
                                i_val=debug_print_level + 1)
      logger => cp_get_default_logger()
      outer_print_level = logger%iter_info%print_level
      logger%iter_info%print_level = silent_print_level
      CALL create_silent_qs_env(input_file, para_env, globenv, qs_env)
      CALL qs_energies(qs_env)
      CALL get_frontier_orbital_scf_result(qs_env, gap, energy, condition_number)
      CALL frontier_orbitals_reference_init(qs_env, reference)
      CALL qs_env_release(qs_env)
      DEALLOCATE (qs_env)
      logger%iter_info%print_level = outer_print_level
      CALL section_vals_release(input_file)

   END SUBROUTINE run_frontier_orbital_reference_scf

! **************************************************************************************************
!> \brief Run one silent SCF using a compact basis held in memory.
!> \param opt_bas ...
!> \param section ...
!> \param para_env ...
!> \param globenv ...
!> \param optimized_basis ...
!> \param gap ...
!> \param energy ...
!> \param condition_number ...
! **************************************************************************************************
   SUBROUTINE run_frontier_orbital_scf(opt_bas, section, para_env, globenv, &
                                       optimized_basis, gap, energy, condition_number)
      TYPE(basis_optimization_type), INTENT(IN)          :: opt_bas
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv
      LOGICAL, INTENT(IN)                                :: optimized_basis
      REAL(KIND=dp), INTENT(OUT)                         :: gap, energy, condition_number

      INTEGER                                            :: outer_print_level
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input_file

      NULLIFY (input_file, logger, qs_env)
      gap = 0.0_dp
      energy = 0.0_dp
      condition_number = 0.0_dp
      CALL section_vals_duplicate(section, input_file)
      CALL prepare_frontier_orbital_scf_input(input_file, opt_bas)
      logger => cp_get_default_logger()
      outer_print_level = logger%iter_info%print_level
      logger%iter_info%print_level = silent_print_level
      CALL create_silent_qs_env(input_file, para_env, globenv, qs_env)
      IF (optimized_basis) CALL update_basis_set(opt_bas, 1, "ORB", qs_env)
      CALL qs_energies(qs_env)
      CALL get_frontier_orbital_scf_result(qs_env, gap, energy, condition_number)
      CALL qs_env_release(qs_env)
      DEALLOCATE (qs_env)
      logger%iter_info%print_level = outer_print_level
      CALL section_vals_release(input_file)

   END SUBROUTINE run_frontier_orbital_scf

! **************************************************************************************************
!> \brief Initialize Quickstep directly under the active CP2K environment without nested startup.
!> \param input_file ...
!> \param para_env ...
!> \param globenv ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE create_silent_qs_env(input_file, para_env, globenv, qs_env)
      TYPE(section_vals_type), POINTER                   :: input_file
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(qs_environment_type), POINTER                 :: qs_env

      TYPE(section_vals_type), POINTER                   :: force_env_section, subsys_section

      force_env_section => section_vals_get_subs_vals(input_file, "FORCE_EVAL")
      subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")
      ALLOCATE (qs_env)
      CALL qs_env_create(qs_env, globenv)
      CALL qs_init(qs_env, para_env, input_file, globenv=globenv, &
                   force_env_section=force_env_section, subsys_section=subsys_section, &
                   use_motion_section=.TRUE., silent=.TRUE.)

   END SUBROUTINE create_silent_qs_env

! **************************************************************************************************
!> \brief Extract the gap, energy, and overlap condition number from a completed SCF.
!> \param qs_env ...
!> \param gap ...
!> \param energy ...
!> \param condition_number ...
! **************************************************************************************************
   SUBROUTINE get_frontier_orbital_scf_result(qs_env, gap, energy, condition_number)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(OUT)                         :: gap, energy, condition_number

      INTEGER                                            :: homo, nmo
      REAL(KIND=dp), DIMENSION(2)                        :: condnum
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: smat
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_energy_type), POINTER                      :: qs_energy

      NULLIFY (blacs_env, eigenvalues, matrix_s, mos, qs_energy, smat)
      gap = 0.0_dp
      CALL get_qs_env(qs_env, blacs_env=blacs_env, energy=qs_energy, matrix_s=matrix_s, mos=mos)
      energy = qs_energy%total
      ALLOCATE (smat(1, 1))
      smat(1, 1)%matrix => matrix_s(1)%matrix
      CALL overlap_condnum(smat, condnum, -1, .FALSE., .TRUE., .FALSE., blacs_env)
      condition_number = condnum(2)
      DEALLOCATE (smat)
      IF (SIZE(mos) /= 1) THEN
         CPABORT("FRONTIER_ORBITAL_SCREENING requires a restricted calculation")
      END IF
      CALL get_mo_set(mos(1), homo=homo, nmo=nmo, eigenvalues=eigenvalues)
      IF (homo < 1 .OR. homo >= nmo) THEN
         CPABORT("FRONTIER_ORBITAL_SCREENING could not determine a HOMO-LUMO gap")
      END IF
      gap = eigenvalues(homo + 1) - eigenvalues(homo)

   END SUBROUTINE get_frontier_orbital_scf_result

! **************************************************************************************************
!> \brief Make a parsed reference input suitable for a quiet, independent screening SCF.
!> \param input_file ...
!> \param opt_bas ...
! **************************************************************************************************
   SUBROUTINE prepare_frontier_orbital_scf_input(input_file, opt_bas)
      TYPE(section_vals_type), POINTER                   :: input_file
      TYPE(basis_optimization_type), INTENT(IN)          :: opt_bas

      CHARACTER(LEN=default_string_length)               :: element
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: basis_info
      INTEGER                                            :: ifit, ikind, nbasis_files, nkinds
      LOGICAL                                            :: element_explicit, found
      TYPE(section_vals_type), POINTER                   :: dft_section, kind_section, subsys_section

      CALL section_vals_val_set(input_file, "GLOBAL%PRINT_LEVEL", i_val=silent_print_level)
      CALL section_vals_val_set(input_file, "FORCE_EVAL%DFT%SCF%SCF_GUESS", i_val=atomic_guess)
      CALL section_vals_val_set(input_file, &
                                "FORCE_EVAL%DFT%SCF%PRINT%RESTART%_SECTION_PARAMETERS_", &
                                i_val=debug_print_level + 1)
      CALL section_vals_val_set(input_file, &
                                "FORCE_EVAL%DFT%SCF%PRINT%RESTART_HISTORY%_SECTION_PARAMETERS_", &
                                i_val=debug_print_level + 1)

      dft_section => section_vals_get_subs_vals(input_file, "FORCE_EVAL%DFT")
      CALL section_vals_val_get(dft_section, "BASIS_SET_FILE_NAME", n_rep_val=nbasis_files)
      IF (nbasis_files < 1) THEN
         CPABORT("The screening reference input has no basis-set file")
      END IF
      CALL section_vals_val_set( &
         dft_section, "BASIS_SET_FILE_NAME", i_rep_val=1, &
         c_val=discover_file(opt_bas%template_basis_file))

      subsys_section => section_vals_get_subs_vals(input_file, "FORCE_EVAL%SUBSYS")
      kind_section => section_vals_get_subs_vals(subsys_section, "KIND")
      CALL section_vals_get(kind_section, n_repetition=nkinds)
      DO ikind = 1, nkinds
         CALL section_vals_val_get(kind_section, "ELEMENT", explicit=element_explicit, &
                                   i_rep_section=ikind)
         IF (element_explicit) THEN
            CALL section_vals_val_get(kind_section, "ELEMENT", c_val=element, i_rep_section=ikind)
         ELSE
            CALL section_vals_val_get(kind_section, "_SECTION_PARAMETERS_", &
                                      c_val=element, i_rep_section=ikind)
         END IF
         CALL uppercase(element)
         found = .FALSE.
         DO ifit = 1, opt_bas%nkind
            IF (TRIM(element) == TRIM(opt_bas%kind_basis(ifit)%element)) THEN
               NULLIFY (basis_info)
               CALL section_vals_val_get(kind_section, "BASIS_SET", c_vals=basis_info, &
                                         i_rep_val=1, i_rep_section=ikind)
               basis_info(SIZE(basis_info)) = opt_bas%kind_basis(ifit)%flex_basis(0)%basis_name
               found = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT. found) THEN
            CALL cp_abort(__LOCATION__, &
                          "Every KIND in the screening reference input must have a corresponding "// &
                          "FIT_KIND: "//TRIM(element))
         END IF
      END DO

   END SUBROUTINE prepare_frontier_orbital_scf_input

! **************************************************************************************************
!> \brief Write all optimized candidates to one CP2K basis file.
!> \param filename ...
!> \param candidates ...
!> \param reference_gap ...
!> \param reference_energy ...
!> \param reference_condition_number ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE write_frontier_orbital_basis_file(filename, candidates, &
                                                reference_gap, reference_energy, reference_condition_number, para_env)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(frontier_orbital_screening_candidate_type), &
         DIMENSION(:)                                    :: candidates
      REAL(KIND=dp), INTENT(IN)                          :: reference_gap, reference_energy, &
                                                            reference_condition_number
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: candidate, ibasis, ikind, unit_nr

      IF (para_env%is_source()) THEN
         CALL open_file(file_name=filename, file_status="UNKNOWN", file_action="WRITE", unit_number=unit_nr)
      ELSE
         unit_nr = -999
      END IF
      DO candidate = 1, SIZE(candidates)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(A)') "#"
            WRITE (unit_nr, '(A,I0)') "# Optimized basis set ", candidate
            WRITE (unit_nr, '(A,F12.6,A)') "# Reference HOMO-LUMO gap: ", reference_gap*evolt, " eV"
            WRITE (unit_nr, '(A,F12.6,A)') "# Basis HOMO-LUMO gap:     ", candidates(candidate)%gap*evolt, " eV"
            WRITE (unit_nr, '(A,F12.6,A)') "# Absolute gap error:      ", &
               ABS(candidates(candidate)%gap - reference_gap)*evolt, " eV"
            WRITE (unit_nr, '(A,ES15.6)') "# Reference condition number: ", reference_condition_number
            WRITE (unit_nr, '(A,ES15.6)') &
               "# Basis condition number:     ", candidates(candidate)%condition_number
            WRITE (unit_nr, '(A,F22.12,A)') "# Reference total energy:  ", reference_energy, " hartree"
            WRITE (unit_nr, '(A,F22.12,A)') "# Basis total energy:      ", candidates(candidate)%energy, " hartree"
            WRITE (unit_nr, '(A)') "#"
         END IF
         DO ikind = 1, candidates(candidate)%opt_bas%nkind
            DO ibasis = 0, candidates(candidate)%opt_bas%kind_basis(ikind)%nbasis_deriv
               candidates(candidate)%opt_bas%kind_basis(ikind)%flex_basis(ibasis)%basis_name = &
                  TRIM(candidates(candidate)%opt_bas%kind_basis(ikind)%flex_basis(ibasis)%basis_name)// &
                  "_"//TRIM(ADJUSTL(cp_to_string(candidate)))
               CALL write_basis(candidates(candidate)%opt_bas%kind_basis(ikind)%flex_basis(ibasis), &
                                candidates(candidate)%opt_bas%kind_basis(ikind)%element, unit_nr)
            END DO
         END DO
      END DO
      IF (para_env%is_source()) CALL close_file(unit_number=unit_nr)

   END SUBROUTINE write_frontier_orbital_basis_file

! **************************************************************************************************
!> \brief Print the two-line column headings shared by screening result tables.
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_frontier_orbital_table_header(unit_nr)
      INTEGER, INTENT(IN)                                :: unit_nr

      WRITE (unit_nr, '(1X,A,T14,A,T24,A,T46,A,T59,A)') &
         "Basis set", "Gap / eV", "Abs. error to", "Condition", "Total energy / Ha"
      WRITE (unit_nr, '(T24,A,T46,A)') "ref. basis / meV", "number"

   END SUBROUTINE print_frontier_orbital_table_header

! **************************************************************************************************
!> \brief Collect reference and initial basis-set names for every fitted atom kind.
!> \param section ...
!> \param opt_bas ...
!> \param element_symbols ...
!> \param reference_basis_names ...
!> \param initial_basis_names ...
! **************************************************************************************************
   SUBROUTINE get_frontier_orbital_basis_sets(section, opt_bas, element_symbols, &
                                              reference_basis_names, initial_basis_names)
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(basis_optimization_type), INTENT(IN)          :: opt_bas
      CHARACTER(LEN=default_string_length), &
         ALLOCATABLE, DIMENSION(:), INTENT(OUT)          :: element_symbols, reference_basis_names, &
                                                            initial_basis_names

      CHARACTER(LEN=default_string_length)               :: element
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: basis_info
      INTEGER                                            :: ifit, ikind, nkinds
      LOGICAL                                            :: element_explicit, found
      TYPE(section_vals_type), POINTER                   :: kind_section, subsys_section

      subsys_section => section_vals_get_subs_vals(section, "FORCE_EVAL%SUBSYS")
      kind_section => section_vals_get_subs_vals(subsys_section, "KIND")
      CALL section_vals_get(kind_section, n_repetition=nkinds)

      ALLOCATE (element_symbols(opt_bas%nkind), reference_basis_names(opt_bas%nkind), &
                initial_basis_names(opt_bas%nkind))
      DO ifit = 1, opt_bas%nkind
         found = .FALSE.
         DO ikind = 1, nkinds
            CALL section_vals_val_get(kind_section, "ELEMENT", explicit=element_explicit, &
                                      i_rep_section=ikind)
            IF (element_explicit) THEN
               CALL section_vals_val_get(kind_section, "ELEMENT", c_val=element, i_rep_section=ikind)
            ELSE
               CALL section_vals_val_get(kind_section, "_SECTION_PARAMETERS_", &
                                         c_val=element, i_rep_section=ikind)
            END IF
            CALL uppercase(element)
            IF (TRIM(element) == TRIM(opt_bas%kind_basis(ifit)%element)) THEN
               NULLIFY (basis_info)
               CALL section_vals_val_get(kind_section, "BASIS_SET", c_vals=basis_info, &
                                         i_rep_val=1, i_rep_section=ikind)
               element_symbols(ifit) = opt_bas%kind_basis(ifit)%element
               reference_basis_names(ifit) = basis_info(SIZE(basis_info))
               initial_basis_names(ifit) = opt_bas%kind_basis(ifit)%basis_name
               found = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT. found) THEN
            CALL cp_abort(__LOCATION__, &
                          "Could not find the reference basis set for fitted atom kind "// &
                          TRIM(opt_bas%kind_basis(ifit)%element))
         END IF
      END DO

   END SUBROUTINE get_frontier_orbital_basis_sets

! **************************************************************************************************
!> \brief Print the compact, aligned human-readable screening summary.
!> \param unit_nr ...
!> \param filename ...
!> \param section ...
!> \param candidates ...
!> \param reference_gap ...
!> \param reference_energy ...
!> \param reference_condition_number ...
!> \param initial_gap ...
!> \param initial_energy ...
!> \param initial_condition_number ...
! **************************************************************************************************
   SUBROUTINE print_frontier_orbital_summary(unit_nr, filename, section, candidates, &
                                             reference_gap, reference_energy, reference_condition_number, &
                                             initial_gap, initial_energy, initial_condition_number)
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(frontier_orbital_screening_candidate_type), &
         DIMENSION(:), INTENT(IN)                        :: candidates
      REAL(KIND=dp), INTENT(IN) :: reference_gap, reference_energy, reference_condition_number, &
         initial_gap, initial_energy, initial_condition_number

      CHARACTER(LEN=default_string_length), &
         ALLOCATABLE, DIMENSION(:)                       :: element_symbols, initial_basis_names, &
                                                            reference_basis_names
      INTEGER                                            :: candidate, separator

      IF (unit_nr <= 0) RETURN
      separator = SCAN(TRIM(filename), "/", BACK=.TRUE.)
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(1X,A)') REPEAT("-", 79)
      WRITE (unit_nr, '(1X,A)') &
         "SCF calcs. with ref. basis set, initial small and optimized basis sets"
      WRITE (unit_nr, '(1X,A)') REPEAT("-", 79)
      WRITE (unit_nr, '(A)') ""
      CALL get_frontier_orbital_basis_sets(section, candidates(1)%opt_bas, element_symbols, &
                                           reference_basis_names, initial_basis_names)
      CALL print_frontier_orbital_basis_sets(unit_nr, element_symbols, reference_basis_names, &
                                             initial_basis_names)
      DEALLOCATE (element_symbols, reference_basis_names, initial_basis_names)
      WRITE (unit_nr, '(A)') ""
      CALL print_frontier_orbital_table_header(unit_nr)
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(1X,A,T14,F8.3,T24,F16.1,T46,ES11.3,T59,F17.12)') &
         "Reference", reference_gap*evolt, 0.0_dp, reference_condition_number, reference_energy
      WRITE (unit_nr, '(1X,A,T14,F8.3,T24,F16.1,T46,ES11.3,T59,F17.12)') &
         "Initial", initial_gap*evolt, ABS(initial_gap - reference_gap)*evolt*1000.0_dp, &
         initial_condition_number, initial_energy
      DO candidate = 1, SIZE(candidates)
         WRITE (unit_nr, '(1X,A,I0,T14,F8.3,T24,F16.1,T46,ES11.3,T59,F17.12)') &
            "Basis ", candidate, candidates(candidate)%gap*evolt, &
            ABS(candidates(candidate)%gap - reference_gap)*evolt*1000.0_dp, &
            candidates(candidate)%condition_number, candidates(candidate)%energy
      END DO
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(1X,A)') REPEAT("-", 79)
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(1X,A,A)') &
         "Basis set file containing all optimized basis sets: ", TRIM(filename(separator + 1:))
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(1X,A)') REPEAT("=", 79)
      WRITE (unit_nr, '(A)') ""

   END SUBROUTINE print_frontier_orbital_summary

! **************************************************************************************************
!> \brief Join a screening working directory and filename.
!> \param directory ...
!> \param filename ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION work_path(directory, filename) RESULT(path)
      CHARACTER(LEN=*), INTENT(IN)                       :: directory, filename
      CHARACTER(LEN=default_path_length)                 :: path

      INTEGER                                            :: n

      n = LEN_TRIM(directory)
      IF (filename(1:1) == "/") THEN
         path = filename
      ELSE IF (n == 0 .OR. directory(n:n) == "/") THEN
         path = TRIM(directory)//filename
      ELSE
         path = TRIM(directory)//"/"//filename
      END IF

   END FUNCTION work_path

END MODULE optbas_frontier_orbital_screening
