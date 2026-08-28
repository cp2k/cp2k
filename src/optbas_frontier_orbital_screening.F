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
                                              get_unit_number,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: debug_print_level,&
                                              silent_print_level
   USE f77_interface,                   ONLY: create_force_env,&
                                              destroy_force_env,&
                                              f_env_add_defaults,&
                                              f_env_rm_defaults,&
                                              f_env_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
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
   USE machine,                         ONLY: default_output_unit,&
                                              m_getcwd
   USE message_passing,                 ONLY: mp_para_env_type
   USE optbas_frontier_orbitals_utils,  ONLY: frontier_orbitals_reference_init,&
                                              frontier_orbitals_reference_release,&
                                              frontier_orbitals_reference_type
   USE optimize_basis_types,            ONLY: basis_optimization_type,&
                                              deallocate_basis_optimization_type,&
                                              method_mo_fit_occ_virtual
   USE optimize_basis_utils,            ONLY: optimize_basis_init_read_input,&
                                              write_basis
   USE physcon,                         ONLY: evolt
   USE qs_condnum,                      ONLY: overlap_condnum
   USE qs_energy,                       ONLY: qs_energies
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_part_release,&
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
      SUBROUTINE basis_optimization_driver(opt_bas, input_declaration, para_env, training_input, reference)
         IMPORT :: basis_optimization_type, frontier_orbitals_reference_type, mp_para_env_type, &
            section_type, section_vals_type
         TYPE(basis_optimization_type)                   :: opt_bas
         TYPE(section_type), POINTER                     :: input_declaration
         TYPE(mp_para_env_type), POINTER                 :: para_env
         TYPE(section_vals_type), OPTIONAL, POINTER      :: training_input
         TYPE(frontier_orbitals_reference_type), &
            DIMENSION(:), INTENT(IN), OPTIONAL, TARGET   :: reference
      END SUBROUTINE basis_optimization_driver
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Run serial frontier-orbital optimizations and rank them by HOMO-LUMO gaps from new SCFs.
! **************************************************************************************************
   SUBROUTINE run_frontier_orbital_screening(input_declaration, section, para_env, frontier_screening_section, &
                                             optimize_basis_driver)
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(mp_para_env_type), POINTER                    :: para_env
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

      CHARACTER(LEN=default_path_length)                 :: candidate_file, final_basis_file, &
                                                            main_dir, frontier_screening_work_file
      CHARACTER(LEN=default_path_length), ALLOCATABLE, &
         DIMENSION(:)                                    :: candidate_files
      INTEGER                                            :: i, n_candidates, n_ranked, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rank_order
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
      frontier_screening_work_file = work_path(main_dir, ".frontier_screening_work.basis")
      IF (para_env%is_source()) THEN
         CALL delete_work_file(frontier_screening_work_file)
      END IF
      CALL para_env%sync()

      ALLOCATE (candidates(n_candidates), candidate_files(n_candidates), rank_order(n_candidates), reference(1))
      logger => cp_get_default_logger()
      unit_nr = -1
      IF (para_env%is_source()) unit_nr = cp_logger_get_default_unit_nr(logger)

      CALL run_frontier_orbital_reference_scf(input_declaration, section, para_env, reference(1), &
                                              reference_gap, reference_energy, &
                                              reference_condition_number)

      DO i = 1, n_candidates
         CALL optimize_basis_init_read_input(candidates(i)%opt_bas, section, para_env, &
                                             quiet_output=.TRUE., embedded_training=.TRUE., &
                                             screening_work_file=frontier_screening_work_file)
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
         candidate_file = work_path(main_dir, ".frontier_screening_candidate_"// &
                                    TRIM(ADJUSTL(cp_to_string(i)))//".basis")
         candidate_files(i) = candidate_file
         candidates(i)%opt_bas%output_basis_file = candidate_file

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(1X,A,I0,A,I0,A)') "Candidate ", i, " of ", n_candidates, &
               ": optimizing basis set ..."
         END IF
         CALL optimize_basis_driver(candidates(i)%opt_bas, input_declaration, para_env, &
                                    training_input=section, reference=reference)
      END DO

      CALL run_frontier_orbital_scf(candidates(1)%opt_bas, input_declaration, section, para_env, &
                                    basis_file=work_path( &
                                    main_dir, discover_file(candidates(1)%opt_bas%template_basis_file)), &
                                    gap=initial_gap, energy=initial_energy, &
                                    condition_number=initial_condition_number)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(A)') " "
         WRITE (unit_nr, '(1X,A)') "SCF results for optimized basis sets"
         WRITE (unit_nr, '(A)') " "
         CALL print_frontier_orbital_table_header(unit_nr)
         WRITE (unit_nr, '(A)') " "
      END IF
      DO i = 1, n_candidates
         CALL run_frontier_orbital_scf(candidates(i)%opt_bas, input_declaration, section, para_env, &
                                       basis_file=candidate_files(i), &
                                       gap=candidates(i)%gap, energy=candidates(i)%energy, &
                                       condition_number=candidates(i)%condition_number)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(1X,A,I0,T14,F8.3,T24,F16.1,T46,ES11.3,T59,F17.12)') &
               "Candidate ", i, candidates(i)%gap*evolt, &
               ABS(candidates(i)%gap - reference_gap)*evolt*1000.0_dp, &
               candidates(i)%condition_number, candidates(i)%energy
         END IF
      END DO
      IF (unit_nr > 0) WRITE (unit_nr, '(A)') " "

      CALL rank_frontier_orbital_candidates(candidates, reference_gap, rank_order, n_ranked)
      CALL write_frontier_orbital_basis_file(final_basis_file, candidates, rank_order, n_ranked, &
                                             reference_gap, reference_energy, reference_condition_number, para_env)
      CALL print_frontier_orbital_summary(unit_nr, final_basis_file, candidates, rank_order, n_ranked, &
                                          reference_gap, reference_energy, reference_condition_number, &
                                          initial_gap, initial_energy, initial_condition_number)

      DO i = 1, n_candidates
         IF (para_env%is_source()) CALL delete_work_file(candidate_files(i))
         CALL deallocate_basis_optimization_type(candidates(i)%opt_bas)
      END DO
      IF (para_env%is_source()) CALL delete_work_file(frontier_screening_work_file)
      CALL para_env%sync()
      CALL frontier_orbitals_reference_release(reference(1))
      DEALLOCATE (candidates, candidate_files, rank_order, reference)

   END SUBROUTINE run_frontier_orbital_screening

! **************************************************************************************************
!> \brief Check the single-file input contract of frontier-orbital screening.
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
! **************************************************************************************************
   SUBROUTINE controlled_parameter_error(parameter_name)
      CHARACTER(LEN=*), INTENT(IN)                       :: parameter_name

      CALL cp_abort(__LOCATION__, &
                    TRIM(parameter_name)//" is not allowed together with FRONTIER_ORBITAL_SCREENING "// &
                    "because it is overwritten automatically")

   END SUBROUTINE controlled_parameter_error

! **************************************************************************************************
!> \brief Run and retain the large-basis reference from the embedded FORCE_EVAL.
! **************************************************************************************************
   SUBROUTINE run_frontier_orbital_reference_scf(input_declaration, section, para_env, reference, &
                                                 gap, energy, condition_number)
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(frontier_orbitals_reference_type), &
         INTENT(OUT)                                     :: reference
      REAL(KIND=dp), INTENT(OUT)                         :: gap, energy, condition_number

      INTEGER                                            :: env_id, stat
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input_file

      NULLIFY (f_env, force_env, input_file, qs_env)
      CALL section_vals_duplicate(section, input_file)
      CALL section_vals_val_set(input_file, "FORCE_EVAL%DFT%SCF%SCF_GUESS", i_val=atomic_guess)
      CALL section_vals_val_set(input_file, "GLOBAL%PROGRAM_RUN_INFO%_SECTION_PARAMETERS_", &
                                i_val=debug_print_level)
      CALL section_vals_val_set(input_file, "GLOBAL%REFERENCES%_SECTION_PARAMETERS_", &
                                i_val=debug_print_level)
      CALL section_vals_val_set(input_file, "GLOBAL%TIMINGS%_SECTION_PARAMETERS_", &
                                i_val=debug_print_level)
      CALL create_force_env(env_id, input_declaration=input_declaration, &
                            input_path="embedded frontier-orbital screening input", input=input_file, &
                            output_unit=default_output_unit, owns_out_unit=.FALSE., &
                            mpi_comm=para_env, ierr=stat)
      IF (stat /= 0) CPABORT("Could not create the embedded reference DFT environment")
      CALL f_env_add_defaults(env_id, f_env)
      force_env => f_env%force_env
      CALL force_env_get(force_env, qs_env=qs_env)
      CALL qs_energies(qs_env)
      CALL get_frontier_orbital_scf_result(qs_env, gap, energy, condition_number)
      CALL frontier_orbitals_reference_init(qs_env, reference)
      CALL f_env_rm_defaults(f_env, stat)
      CALL qs_env_part_release(qs_env)
      CALL destroy_force_env(env_id, stat)
      CALL section_vals_release(input_file)

   END SUBROUTINE run_frontier_orbital_reference_scf

! **************************************************************************************************
!> \brief Run one silent SCF using a compact basis file.
! **************************************************************************************************
   SUBROUTINE run_frontier_orbital_scf(opt_bas, input_declaration, section, para_env, &
                                       basis_file, gap, energy, condition_number)
      TYPE(basis_optimization_type), INTENT(IN)          :: opt_bas
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(mp_para_env_type), POINTER                    :: para_env
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_file
      REAL(KIND=dp), INTENT(OUT)                         :: gap, energy, condition_number

      INTEGER                                            :: env_id, scf_output_unit, stat
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input_file

      NULLIFY (f_env, force_env, input_file, qs_env)
      gap = 0.0_dp
      energy = 0.0_dp
      condition_number = 0.0_dp
      CALL section_vals_duplicate(section, input_file)
      CALL prepare_frontier_orbital_scf_input(input_file, opt_bas, basis_file)
      scf_output_unit = -1
      IF (para_env%is_source()) THEN
         scf_output_unit = get_unit_number()
         IF (scf_output_unit < 1) THEN
            CPABORT("No logical unit is available for disposable screening SCF output")
         END IF
         OPEN (UNIT=scf_output_unit, STATUS="SCRATCH", ACTION="READWRITE", IOSTAT=stat)
         IF (stat /= 0) THEN
            CPABORT("Could not open disposable frontier-orbital screening SCF output")
         END IF
      END IF
      CALL create_force_env(env_id, input_declaration=input_declaration, &
                            input_path="embedded frontier-orbital screening input", input=input_file, &
                            output_unit=scf_output_unit, owns_out_unit=.FALSE., &
                            mpi_comm=para_env, ierr=stat)
      IF (stat /= 0) THEN
         CPABORT("Could not create a frontier-orbital screening SCF environment")
      END IF
      CALL f_env_add_defaults(env_id, f_env)
      force_env => f_env%force_env
      CALL force_env_get(force_env, qs_env=qs_env)
      CALL qs_energies(qs_env)
      CALL get_frontier_orbital_scf_result(qs_env, gap, energy, condition_number)
      CALL f_env_rm_defaults(f_env, stat)
      CALL qs_env_part_release(qs_env)
      CALL destroy_force_env(env_id, stat)
      CALL section_vals_release(input_file)
      IF (para_env%is_source()) CLOSE (scf_output_unit)

   END SUBROUTINE run_frontier_orbital_scf

! **************************************************************************************************
!> \brief Extract the gap, energy, and overlap condition number from a completed SCF.
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
! **************************************************************************************************
   SUBROUTINE prepare_frontier_orbital_scf_input(input_file, opt_bas, basis_file)
      TYPE(section_vals_type), POINTER                   :: input_file
      TYPE(basis_optimization_type), INTENT(IN)          :: opt_bas
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_file

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
      CALL section_vals_val_set(dft_section, "BASIS_SET_FILE_NAME", i_rep_val=1, c_val=basis_file)

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
!> \brief Rank candidates by absolute self-consistent gap error.
! **************************************************************************************************
   SUBROUTINE rank_frontier_orbital_candidates(candidates, reference_gap, rank_order, n_ranked)
      TYPE(frontier_orbital_screening_candidate_type), &
         DIMENSION(:), INTENT(IN)                        :: candidates
      REAL(KIND=dp), INTENT(IN)                          :: reference_gap
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: rank_order
      INTEGER, INTENT(OUT)                               :: n_ranked

      INTEGER                                            :: i, j, saved

      rank_order = 0
      n_ranked = 0
      DO i = 1, SIZE(candidates)
         n_ranked = n_ranked + 1
         rank_order(n_ranked) = i
         j = n_ranked
         DO WHILE (j > 1)
            IF (ABS(candidates(rank_order(j))%gap - reference_gap) >= &
                ABS(candidates(rank_order(j - 1))%gap - reference_gap)) EXIT
            saved = rank_order(j - 1)
            rank_order(j - 1) = rank_order(j)
            rank_order(j) = saved
            j = j - 1
         END DO
      END DO

   END SUBROUTINE rank_frontier_orbital_candidates

! **************************************************************************************************
!> \brief Write all ranked candidates to one CP2K basis file.
! **************************************************************************************************
   SUBROUTINE write_frontier_orbital_basis_file(filename, candidates, rank_order, n_ranked, &
                                                reference_gap, reference_energy, reference_condition_number, para_env)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(frontier_orbital_screening_candidate_type), &
         DIMENSION(:)                                    :: candidates
      INTEGER, DIMENSION(:), INTENT(IN)                  :: rank_order
      INTEGER, INTENT(IN)                                :: n_ranked
      REAL(KIND=dp), INTENT(IN)                          :: reference_gap, reference_energy, &
                                                            reference_condition_number
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: candidate, ibasis, ikind, rank, unit_nr

      IF (para_env%is_source()) THEN
         CALL open_file(file_name=filename, file_status="UNKNOWN", file_action="WRITE", unit_number=unit_nr)
      ELSE
         unit_nr = -999
      END IF
      DO rank = 1, n_ranked
         candidate = rank_order(rank)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(A)') "#"
            WRITE (unit_nr, '(A,I0)') "# Optimized basis-set rank ", rank
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
                  "_"//TRIM(ADJUSTL(cp_to_string(rank)))
               CALL write_basis(candidates(candidate)%opt_bas%kind_basis(ikind)%flex_basis(ibasis), &
                                candidates(candidate)%opt_bas%kind_basis(ikind)%element, unit_nr)
            END DO
         END DO
      END DO
      IF (para_env%is_source()) CALL close_file(unit_number=unit_nr)

   END SUBROUTINE write_frontier_orbital_basis_file

! **************************************************************************************************
!> \brief Print the two-line column headings shared by screening result tables.
! **************************************************************************************************
   SUBROUTINE print_frontier_orbital_table_header(unit_nr)
      INTEGER, INTENT(IN)                                :: unit_nr

      WRITE (unit_nr, '(1X,A,T14,A,T24,A,T46,A,T59,A)') &
         "Basis set", "Gap / eV", "Abs. error to", "Condition", "Total energy / Ha"
      WRITE (unit_nr, '(T24,A,T46,A)') "ref. basis / meV", "number"

   END SUBROUTINE print_frontier_orbital_table_header

! **************************************************************************************************
!> \brief Print the compact, aligned human-readable screening summary.
! **************************************************************************************************
   SUBROUTINE print_frontier_orbital_summary(unit_nr, filename, candidates, rank_order, n_ranked, &
                                             reference_gap, reference_energy, reference_condition_number, &
                                             initial_gap, initial_energy, initial_condition_number)
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(frontier_orbital_screening_candidate_type), &
         DIMENSION(:), INTENT(IN)                        :: candidates
      INTEGER, DIMENSION(:), INTENT(IN)                  :: rank_order
      INTEGER, INTENT(IN)                                :: n_ranked
      REAL(KIND=dp), INTENT(IN) :: reference_gap, reference_energy, reference_condition_number, &
         initial_gap, initial_energy, initial_condition_number

      INTEGER                                            :: candidate, rank, separator

      IF (unit_nr <= 0) RETURN
      separator = SCAN(TRIM(filename), "/", BACK=.TRUE.)
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(1X,A)') REPEAT("-", 79)
      WRITE (unit_nr, '(1X,A)') &
         "SCF calculations with reference basis set and initial small basis set"
      WRITE (unit_nr, '(1X,A)') REPEAT("-", 79)
      WRITE (unit_nr, '(A)') ""
      CALL print_frontier_orbital_table_header(unit_nr)
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(1X,A,T14,F8.3,T24,F16.1,T46,ES11.3,T59,F17.12)') &
         "Reference", reference_gap*evolt, 0.0_dp, reference_condition_number, reference_energy
      WRITE (unit_nr, '(1X,A,T14,F8.3,T24,F16.1,T46,ES11.3,T59,F17.12)') &
         "Initial", initial_gap*evolt, ABS(initial_gap - reference_gap)*evolt*1000.0_dp, &
         initial_condition_number, initial_energy
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(1X,A)') REPEAT("-", 79)
      WRITE (unit_nr, '(1X,A)') "Ranking of optimized basis sets"
      WRITE (unit_nr, '(1X,A)') REPEAT("-", 79)
      WRITE (unit_nr, '(A)') ""
      CALL print_frontier_orbital_table_header(unit_nr)
      WRITE (unit_nr, '(A)') ""
      DO rank = 1, n_ranked
         candidate = rank_order(rank)
         WRITE (unit_nr, '(1X,A,I0,T14,F8.3,T24,F16.1,T46,ES11.3,T59,F17.12)') &
            "Basis ", rank, candidates(candidate)%gap*evolt, &
            ABS(candidates(candidate)%gap - reference_gap)*evolt*1000.0_dp, &
            candidates(candidate)%condition_number, candidates(candidate)%energy
      END DO
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(1X,A)') REPEAT("-", 79)
      IF (ABS(initial_gap - reference_gap) + 1.0E-6_dp/evolt < &
          ABS(candidates(rank_order(1))%gap - reference_gap)) THEN
         WRITE (unit_nr, '(1X,A)') &
            "Note: the initial basis matches the reference gap more closely than Basis 1."
      END IF
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(1X,A,A)') "Ranked basis set file: ", TRIM(filename(separator + 1:))
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(1X,A)') REPEAT("=", 79)
      WRITE (unit_nr, '(A)') ""

   END SUBROUTINE print_frontier_orbital_summary

! **************************************************************************************************
!> \brief Delete a private screening work file if it exists.
! **************************************************************************************************
   SUBROUTINE delete_work_file(filename)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename

      INTEGER                                            :: ios, unit_nr
      LOGICAL                                            :: exists

      INQUIRE (FILE=filename, EXIST=exists)
      IF (.NOT. exists) RETURN
      OPEN (NEWUNIT=unit_nr, FILE=filename, STATUS="OLD", ACTION="READWRITE", IOSTAT=ios)
      IF (ios == 0) CLOSE (unit_nr, STATUS="DELETE")

   END SUBROUTINE delete_work_file

! **************************************************************************************************
!> \brief Join a screening working directory and filename.
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
