!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE optimize_basis
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              discover_file,&
                                              get_unit_number,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: debug_print_level,&
                                              silent_print_level
   USE f77_interface,                   ONLY: create_force_env,&
                                              destroy_force_env,&
                                              f_env_add_defaults,&
                                              f_env_get_from_id,&
                                              f_env_rm_defaults,&
                                              f_env_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_constants,                 ONLY: atomic_guess
   USE input_cp2k_read,                 ONLY: empty_initial_variables,&
                                              read_input
   USE input_section_types,             ONLY: section_type,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_release,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE machine,                         ONLY: m_chdir,&
                                              m_getcwd,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type
   USE optbas_fenv_manipulation,        ONLY: allocate_mo_sets,&
                                              calculate_ks_matrix,&
                                              calculate_overlap_inverse,&
                                              modify_input_settings,&
                                              update_basis_set
   USE optbas_frontier_orbitals_utils,  ONLY: evaluate_frontier_orbitals_objective,&
                                              frontier_orbitals_reference_init,&
                                              frontier_orbitals_reference_release,&
                                              frontier_orbitals_reference_type,&
                                              frontier_orbitals_result_type,&
                                              get_ao_angular_momenta,&
                                              max_diagnostic_angular_momentum
   USE optbas_opt_utils,                ONLY: evaluate_optvals,&
                                              fit_mo_coeffs,&
                                              optbas_build_neighborlist
   USE optimize_basis_types,            ONLY: basis_optimization_type,&
                                              deallocate_basis_optimization_type,&
                                              method_mo_fit_occ,&
                                              method_mo_fit_occ_virtual,&
                                              subset_type
   USE optimize_basis_utils,            ONLY: get_set_and_basis_id,&
                                              optimize_basis_init_read_input,&
                                              update_derived_basis_sets,&
                                              write_basis
   USE physcon,                         ONLY: evolt
   USE powell,                          ONLY: powell_optimize
   USE qs_energy,                       ONLY: qs_energies
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_part_release,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE string_utilities,                ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: run_optimize_basis

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_basis'

   INTEGER, PARAMETER :: number_screening_presets = 5
   CHARACTER(LEN=16), DIMENSION(number_screening_presets), PARAMETER :: screening_preset_names = &
                                                   [CHARACTER(LEN=16) :: "requested", "all30", "virtual30", "lowgap", "regularized"]
   ! Rows are occupied, virtual, empty, gap, condition, and coefficient weights.
   REAL(KIND=dp), DIMENSION(6, number_screening_presets), PARAMETER :: screening_weights = RESHAPE([ &
                                                                           10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 0.001_dp, 0.001_dp, &
                                                                           30.0_dp, 30.0_dp, 30.0_dp, 10.0_dp, 0.001_dp, 0.001_dp, &
                                                                           10.0_dp, 30.0_dp, 30.0_dp, 10.0_dp, 0.001_dp, 0.001_dp, &
                                                                            10.0_dp, 30.0_dp, 30.0_dp, 3.0_dp, 0.001_dp, 0.001_dp, &
                                                                           10.0_dp, 30.0_dp, 30.0_dp, 10.0_dp, 0.001_dp, 0.01_dp], &
                                                                                                   [6, number_screening_presets])

   TYPE screening_candidate_type
      INTEGER :: preset = 0
      LOGICAL :: scf_converged = .FALSE.
      REAL(KIND=dp) :: gap = 0.0_dp
      REAL(KIND=dp) :: energy = 0.0_dp
      TYPE(basis_optimization_type) :: opt_bas
   END TYPE screening_candidate_type

   INTEGER, PARAMETER :: frontier_idx_loss_occupied = 1
   INTEGER, PARAMETER :: frontier_idx_loss_virtual = 2
   INTEGER, PARAMETER :: frontier_idx_loss_gap = 3
   INTEGER, PARAMETER :: frontier_idx_overlap_occupied = 4
   INTEGER, PARAMETER :: frontier_idx_gap_reference = 5
   INTEGER, PARAMETER :: frontier_idx_gap_candidate = 6
   INTEGER, PARAMETER :: frontier_idx_condition_candidate = 7
   INTEGER, PARAMETER :: frontier_idx_overlap_empty = 8
   INTEGER, PARAMETER :: frontier_idx_time = 9
   INTEGER, PARAMETER :: frontier_idx_loss_empty = 10
   INTEGER, PARAMETER :: frontier_idx_loss_coefficient = 11
   INTEGER, PARAMETER :: frontier_idx_similarity_virtual = 12
   INTEGER, PARAMETER :: frontier_idx_condition_reference = 13
   INTEGER, PARAMETER :: frontier_idx_nao_reference = 14
   INTEGER, PARAMETER :: frontier_idx_nao_candidate = 15
   INTEGER, PARAMETER :: frontier_idx_reference_count = 16
   INTEGER, PARAMETER :: frontier_idx_candidate_count = &
                         frontier_idx_reference_count + max_diagnostic_angular_momentum + 1
   INTEGER, PARAMETER :: number_frontier_diagnostics = &
                         frontier_idx_candidate_count + max_diagnostic_angular_momentum

CONTAINS

! **************************************************************************************************
!> \brief main entry point for methods aimed at optimizing basis sets
!> \param input_declaration ...
!> \param root_section ...
!> \param para_env ...
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE run_optimize_basis(input_declaration, root_section, para_env)
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(section_vals_type), POINTER                   :: root_section
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'run_optimize_basis'

      INTEGER                                            :: handle
      LOGICAL                                            :: screening_explicit
      TYPE(basis_optimization_type)                      :: opt_bas
      TYPE(section_vals_type), POINTER                   :: optbas_section, screening_section

      CALL timeset(routineN, handle)

      optbas_section => section_vals_get_subs_vals(root_section, "OPTIMIZE_BASIS")
      screening_section => section_vals_get_subs_vals(optbas_section, "FRONTIER_ORBITAL_SCREENING")
      CALL section_vals_get(screening_section, explicit=screening_explicit)

      IF (screening_explicit) THEN
         CALL run_frontier_orbital_screening(input_declaration, root_section, para_env, screening_section)
      ELSE
         CALL optimize_basis_init_read_input(opt_bas, root_section, para_env)
         CALL driver_para_opt_basis(opt_bas, input_declaration, para_env)
         CALL deallocate_basis_optimization_type(opt_bas)
      END IF

      CALL timestop(handle)

   END SUBROUTINE run_optimize_basis

! **************************************************************************************************
!> \brief Run serial frontier-orbital optimizations and rank them using fresh SCF gaps.
! **************************************************************************************************
   SUBROUTINE run_frontier_orbital_screening(input_declaration, root_section, para_env, screening_section)
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(section_vals_type), POINTER                   :: root_section
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: screening_section

      CHARACTER(LEN=default_path_length)                 :: candidate_file, final_basis_file, &
                                                            main_dir
      CHARACTER(LEN=default_path_length), ALLOCATABLE, &
         DIMENSION(:)                                    :: candidate_files
      INTEGER                                            :: i, n_candidates, n_ranked, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rank_order
      LOGICAL                                            :: initial_converged, reference_converged
      REAL(KIND=dp)                                      :: initial_energy, initial_gap, &
                                                            reference_energy, reference_gap
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(screening_candidate_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: candidates

      CALL section_vals_val_get(screening_section, "NUMBER_OF_OPTIMIZATIONS", i_val=n_candidates)
      IF (n_candidates < 1 .OR. n_candidates > number_screening_presets) THEN
         CALL cp_abort(__LOCATION__, &
                       "NUMBER_OF_OPTIMIZATIONS must be between 1 and "// &
                       TRIM(cp_to_string(number_screening_presets)))
      END IF

      CALL m_getcwd(main_dir)

      ALLOCATE (candidates(n_candidates), candidate_files(n_candidates), rank_order(n_candidates))
      logger => cp_get_default_logger()
      unit_nr = -1
      IF (para_env%is_source()) unit_nr = cp_logger_get_default_unit_nr(logger)

      DO i = 1, n_candidates
         candidates(i)%preset = i
         CALL optimize_basis_init_read_input(candidates(i)%opt_bas, root_section, para_env, quiet_output=.TRUE.)
         IF (candidates(i)%opt_bas%method /= method_mo_fit_occ_virtual) THEN
            CPABORT("FRONTIER_ORBITAL_SCREENING requires FRONTIER_ORBITALS")
         END IF
         IF (candidates(i)%opt_bas%ntraining_sets /= 1 .OR. candidates(i)%opt_bas%ncombinations /= 1) THEN
            CALL cp_abort(__LOCATION__, &
                          "FRONTIER_ORBITAL_SCREENING currently requires one training system "// &
                          "and one basis combination")
         END IF

         IF (i == 1) THEN
            final_basis_file = screening_path(main_dir, candidates(i)%opt_bas%output_basis_file)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(A)') ""
               WRITE (unit_nr, '(A)') REPEAT("=", 92)
               WRITE (unit_nr, '(A)') "FRONTIER-ORBITAL BASIS SCREENING"
               WRITE (unit_nr, '(A)') REPEAT("=", 92)
               WRITE (unit_nr, '(A,I0,A)') "Running ", n_candidates, &
                  " built-in weight presets serially; validation SCFs use an atomic guess."
               WRITE (unit_nr, '(A)') "Detailed optimizer evaluations and validation SCF output are suppressed."
               WRITE (unit_nr, '(A)') ""
               WRITE (unit_nr, '(A)') "Computing reference and initial-basis SCF baselines ..."
            END IF
            CALL run_screening_scf(candidates(i)%opt_bas, input_declaration, para_env, &
                                   use_reference_basis=.TRUE., basis_file="", &
                                   converged=reference_converged, gap=reference_gap, energy=reference_energy)
            CALL run_screening_scf(candidates(i)%opt_bas, input_declaration, para_env, &
                                   use_reference_basis=.FALSE., &
                                   basis_file=discover_file(candidates(i)%opt_bas%template_basis_file), &
                                   converged=initial_converged, gap=initial_gap, energy=initial_energy)
            IF (.NOT. reference_converged) CPABORT("The reference-basis screening SCF did not converge")
            IF (.NOT. initial_converged) CPABORT("The initial-basis screening SCF did not converge")
         END IF

         candidates(i)%opt_bas%occupied_weight = screening_weights(1, i)
         candidates(i)%opt_bas%virtual_weight = screening_weights(2, i)
         candidates(i)%opt_bas%empty_overlap_weight = screening_weights(3, i)
         candidates(i)%opt_bas%gap_weight = screening_weights(4, i)
         candidates(i)%opt_bas%condition_weight = screening_weights(5, i)
         candidates(i)%opt_bas%coefficient_weight = screening_weights(6, i)
         candidate_file = screening_path(main_dir, ".frontier_screening_candidate_"// &
                                         TRIM(ADJUSTL(cp_to_string(i)))//".basis")
         candidate_files(i) = candidate_file
         candidates(i)%opt_bas%output_basis_file = candidate_file

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(A,I0,A,I0,A,A,A)') "Candidate ", i, " of ", n_candidates, &
               ": optimizing preset '", TRIM(screening_preset_names(i)), "' ..."
         END IF
         CALL driver_para_opt_basis(candidates(i)%opt_bas, input_declaration, para_env)
         CALL run_screening_scf(candidates(i)%opt_bas, input_declaration, para_env, &
                                use_reference_basis=.FALSE., basis_file=candidate_file, &
                                converged=candidates(i)%scf_converged, &
                                gap=candidates(i)%gap, energy=candidates(i)%energy)
         IF (unit_nr > 0) THEN
            IF (candidates(i)%scf_converged) THEN
               WRITE (unit_nr, '(A,F10.6,A,F9.3,A)') "  SCF gap ", candidates(i)%gap*evolt, &
                  " eV; absolute reference error ", ABS(candidates(i)%gap - reference_gap)*evolt*1000.0_dp, " meV"
            ELSE
               WRITE (unit_nr, '(A)') "  SCF validation did not converge; candidate will not be ranked."
            END IF
         END IF
      END DO

      CALL rank_screening_candidates(candidates, reference_gap, rank_order, n_ranked)
      IF (n_ranked == 0) CPABORT("No frontier-orbital screening candidate produced a converged SCF")
      CALL write_screening_basis_file(final_basis_file, candidates, rank_order, n_ranked, &
                                      reference_gap, reference_energy, para_env)
      CALL write_screening_csv(TRIM(final_basis_file)//".screening.csv", candidates, rank_order, n_ranked, &
                               reference_gap, reference_energy, initial_gap, initial_energy, para_env)
      CALL print_screening_summary(unit_nr, final_basis_file, candidates, rank_order, n_ranked, &
                                   reference_gap, reference_energy, initial_gap, initial_energy)

      DO i = 1, n_candidates
         IF (para_env%is_source()) CALL delete_screening_file(candidate_files(i))
         CALL deallocate_basis_optimization_type(candidates(i)%opt_bas)
      END DO
      CALL para_env%sync()
      DEALLOCATE (candidates, candidate_files, rank_order)

   END SUBROUTINE run_frontier_orbital_screening

! **************************************************************************************************
!> \brief Run one silent SCF using either the reference or a compact basis file.
! **************************************************************************************************
   SUBROUTINE run_screening_scf(opt_bas, input_declaration, para_env, use_reference_basis, &
                                basis_file, converged, gap, energy)
      TYPE(basis_optimization_type), INTENT(IN)          :: opt_bas
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: use_reference_basis
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_file
      LOGICAL, INTENT(OUT)                               :: converged
      REAL(KIND=dp), INTENT(OUT)                         :: gap, energy

      CHARACTER(LEN=default_path_length)                 :: main_dir
      INTEGER                                            :: env_id, homo, ierr, nmo, &
                                                            scf_output_unit, stat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_energy_type), POINTER                      :: qs_energy
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input_file

      NULLIFY (eigenvalues, f_env, force_env, input_file, mos, qs_energy, qs_env)
      converged = .FALSE.
      gap = 0.0_dp
      energy = 0.0_dp
      CALL m_getcwd(main_dir)
      CALL m_chdir(TRIM(opt_bas%training_dir(1)), ierr)
      IF (ierr /= 0) CPABORT("Could not enter the frontier-orbital screening training directory")

      input_file => read_input(input_declaration, opt_bas%training_input(1), &
                               initial_variables=empty_initial_variables, para_env=para_env)
      CALL prepare_screening_scf_input(input_file, opt_bas, use_reference_basis, basis_file)
      scf_output_unit = -1
      IF (para_env%is_source()) THEN
         scf_output_unit = get_unit_number()
         IF (scf_output_unit < 1) THEN
            CPABORT("No logical unit is available for disposable screening SCF output")
         END IF
         OPEN (UNIT=scf_output_unit, STATUS="SCRATCH", ACTION="READWRITE", IOSTAT=stat)
         IF (stat /= 0) CPABORT("Could not open disposable frontier-orbital screening SCF output")
      END IF
      CALL create_force_env(env_id, input_declaration=input_declaration, &
                            input_path=opt_bas%training_input(1), input=input_file, &
                            output_unit=scf_output_unit, owns_out_unit=.FALSE., &
                            mpi_comm=para_env, ierr=stat)
      IF (stat /= 0) CPABORT("Could not create a frontier-orbital screening SCF environment")
      CALL f_env_add_defaults(env_id, f_env)
      force_env => f_env%force_env
      CALL force_env_get(force_env, qs_env=qs_env)
      CALL qs_energies(qs_env, scf_converged=converged)
      CALL get_qs_env(qs_env, energy=qs_energy, mos=mos)
      energy = qs_energy%total
      IF (converged .AND. SIZE(mos) == 1) THEN
         CALL get_mo_set(mos(1), homo=homo, nmo=nmo, eigenvalues=eigenvalues)
         IF (homo >= 1 .AND. homo < nmo) THEN
            gap = eigenvalues(homo + 1) - eigenvalues(homo)
         ELSE
            converged = .FALSE.
         END IF
      ELSE
         converged = .FALSE.
      END IF
      CALL f_env_rm_defaults(f_env, stat)
      CALL qs_env_part_release(qs_env)
      CALL destroy_force_env(env_id, stat)
      CALL section_vals_release(input_file)
      IF (para_env%is_source()) CLOSE (scf_output_unit)
      CALL m_chdir(TRIM(main_dir), ierr)
      IF (ierr /= 0) CPABORT("Could not return from the frontier-orbital screening training directory")

   END SUBROUTINE run_screening_scf

! **************************************************************************************************
!> \brief Make a parsed reference input suitable for a quiet, independent screening SCF.
! **************************************************************************************************
   SUBROUTINE prepare_screening_scf_input(input_file, opt_bas, use_reference_basis, basis_file)
      TYPE(section_vals_type), POINTER                   :: input_file
      TYPE(basis_optimization_type), INTENT(IN)          :: opt_bas
      LOGICAL, INTENT(IN)                                :: use_reference_basis
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
      IF (use_reference_basis) RETURN

      dft_section => section_vals_get_subs_vals(input_file, "FORCE_EVAL%DFT")
      CALL section_vals_val_get(dft_section, "BASIS_SET_FILE_NAME", n_rep_val=nbasis_files)
      IF (nbasis_files < 1) CPABORT("The screening reference input has no basis-set file")
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

   END SUBROUTINE prepare_screening_scf_input

! **************************************************************************************************
!> \brief Rank converged candidates by absolute self-consistent gap error.
! **************************************************************************************************
   SUBROUTINE rank_screening_candidates(candidates, reference_gap, rank_order, n_ranked)
      TYPE(screening_candidate_type), DIMENSION(:), &
         INTENT(IN)                                      :: candidates
      REAL(KIND=dp), INTENT(IN)                          :: reference_gap
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: rank_order
      INTEGER, INTENT(OUT)                               :: n_ranked

      INTEGER                                            :: i, j, saved

      rank_order = 0
      n_ranked = 0
      DO i = 1, SIZE(candidates)
         IF (.NOT. candidates(i)%scf_converged) CYCLE
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

   END SUBROUTINE rank_screening_candidates

! **************************************************************************************************
!> \brief Write all ranked candidates to one CP2K basis file.
! **************************************************************************************************
   SUBROUTINE write_screening_basis_file(filename, candidates, rank_order, n_ranked, &
                                         reference_gap, reference_energy, para_env)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(screening_candidate_type), DIMENSION(:)       :: candidates
      INTEGER, DIMENSION(:), INTENT(IN)                  :: rank_order
      INTEGER, INTENT(IN)                                :: n_ranked
      REAL(KIND=dp), INTENT(IN)                          :: reference_gap, reference_energy
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
            WRITE (unit_nr, '(A,I0,A,A)') "# Frontier-orbital screening rank ", rank, &
               "; preset: ", TRIM(screening_preset_names(candidates(candidate)%preset))
            WRITE (unit_nr, '(A,F12.6,A)') "# Reference HOMO-LUMO gap: ", reference_gap*evolt, " eV"
            WRITE (unit_nr, '(A,F12.6,A)') "# Basis HOMO-LUMO gap:     ", candidates(candidate)%gap*evolt, " eV"
            WRITE (unit_nr, '(A,F12.6,A)') "# Absolute gap error:      ", &
               ABS(candidates(candidate)%gap - reference_gap)*evolt, " eV"
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

   END SUBROUTINE write_screening_basis_file

! **************************************************************************************************
!> \brief Write machine-readable screening results next to the ranked basis file.
! **************************************************************************************************
   SUBROUTINE write_screening_csv(filename, candidates, rank_order, n_ranked, reference_gap, &
                                  reference_energy, initial_gap, initial_energy, para_env)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(screening_candidate_type), DIMENSION(:), &
         INTENT(IN)                                      :: candidates
      INTEGER, DIMENSION(:), INTENT(IN)                  :: rank_order
      INTEGER, INTENT(IN)                                :: n_ranked
      REAL(KIND=dp), INTENT(IN)                          :: reference_gap, reference_energy, &
                                                            initial_gap, initial_energy
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: candidate, rank, unit_nr

      IF (.NOT. para_env%is_source()) RETURN
      CALL open_file(file_name=filename, file_status="UNKNOWN", file_action="WRITE", unit_number=unit_nr)
      WRITE (unit_nr, '(A,A)') &
         "rank,label,preset,gap_ev,absolute_gap_error_ev,total_energy_hartree,reference_gap_ev,", &
         "reference_total_energy_hartree"
      WRITE (unit_nr, '(A,F0.10,A,F0.14,A,F0.10,A,F0.14)') &
         "0,reference,reference,", reference_gap*evolt, ",0.0,", reference_energy, ",", &
         reference_gap*evolt, ",", reference_energy
      WRITE (unit_nr, '(A,F0.10,A,F0.10,A,F0.14,A,F0.10,A,F0.14)') &
         "0,initial,initial,", initial_gap*evolt, ",", ABS(initial_gap - reference_gap)*evolt, ",", &
         initial_energy, ",", reference_gap*evolt, ",", reference_energy
      DO rank = 1, n_ranked
         candidate = rank_order(rank)
         WRITE (unit_nr, '(I0,A,I0,A,A,A,F0.10,A,F0.10,A,F0.14,A,F0.10,A,F0.14)') &
            rank, ",basis_", rank, ",", TRIM(screening_preset_names(candidates(candidate)%preset)), ",", &
            candidates(candidate)%gap*evolt, ",", ABS(candidates(candidate)%gap - reference_gap)*evolt, ",", &
            candidates(candidate)%energy, ",", reference_gap*evolt, ",", reference_energy
      END DO
      CALL close_file(unit_number=unit_nr)

   END SUBROUTINE write_screening_csv

! **************************************************************************************************
!> \brief Print the compact, aligned human-readable screening summary.
! **************************************************************************************************
   SUBROUTINE print_screening_summary(unit_nr, filename, candidates, rank_order, n_ranked, &
                                      reference_gap, reference_energy, initial_gap, initial_energy)
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(screening_candidate_type), DIMENSION(:), &
         INTENT(IN)                                      :: candidates
      INTEGER, DIMENSION(:), INTENT(IN)                  :: rank_order
      INTEGER, INTENT(IN)                                :: n_ranked
      REAL(KIND=dp), INTENT(IN)                          :: reference_gap, reference_energy, &
                                                            initial_gap, initial_energy

      INTEGER                                            :: candidate, i, rank

      IF (unit_nr <= 0) RETURN
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(A)') REPEAT("-", 92)
      WRITE (unit_nr, '(A)') "SELF-CONSISTENT BASELINES"
      WRITE (unit_nr, '(A)') REPEAT("-", 92)
      WRITE (unit_nr, '(A24,2X,A14,2X,A22)') "Basis", "Gap / eV", "Total energy / Ha"
      WRITE (unit_nr, '(A24,2X,F14.6,2X,F22.12)') "Reference", reference_gap*evolt, reference_energy
      WRITE (unit_nr, '(A24,2X,F14.6,2X,F22.12)') "Initial compact", initial_gap*evolt, initial_energy
      WRITE (unit_nr, '(A)') ""
      WRITE (unit_nr, '(A)') "RANKED OPTIMIZED BASES (ranking uses only absolute SCF gap error)"
      WRITE (unit_nr, '(A)') REPEAT("-", 92)
      WRITE (unit_nr, '(A4,2X,A16,2X,A10,2X,A14,2X,A22)') &
         "Rank", "Preset", "Gap / eV", "Error / meV", "Total energy / Ha"
      DO rank = 1, n_ranked
         candidate = rank_order(rank)
         WRITE (unit_nr, '(I4,2X,A16,2X,F10.6,2X,F14.3,2X,F22.12)') &
            rank, TRIM(screening_preset_names(candidates(candidate)%preset)), &
            candidates(candidate)%gap*evolt, ABS(candidates(candidate)%gap - reference_gap)*evolt*1000.0_dp, &
            candidates(candidate)%energy
      END DO
      WRITE (unit_nr, '(A)') REPEAT("-", 92)
      WRITE (unit_nr, '(A)') "Best optimized basis: basis_1 (all definition names end in _1)."
      IF (ABS(initial_gap - reference_gap) + 1.0E-6_dp/evolt < &
          ABS(candidates(rank_order(1))%gap - reference_gap)) THEN
         WRITE (unit_nr, '(A)') "Note: the initial compact basis matches the reference gap more closely than basis_1."
      END IF
      WRITE (unit_nr, '(A,A)') "Ranked basis file: ", TRIM(filename)
      WRITE (unit_nr, '(A,A)') "Machine-readable results: ", TRIM(filename)//".screening.csv"
      DO i = 1, SIZE(candidates)
         IF (.NOT. candidates(i)%scf_converged) THEN
            WRITE (unit_nr, '(A)') "One or more candidate SCFs failed and were excluded from the ranking."
            EXIT
         END IF
      END DO
      WRITE (unit_nr, '(A)') REPEAT("=", 92)
      WRITE (unit_nr, '(A)') ""

   END SUBROUTINE print_screening_summary

! **************************************************************************************************
!> \brief Delete a private screening work file if it exists.
! **************************************************************************************************
   SUBROUTINE delete_screening_file(filename)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename

      INTEGER                                            :: ios, unit_nr
      LOGICAL                                            :: exists

      INQUIRE (FILE=filename, EXIST=exists)
      IF (.NOT. exists) RETURN
      OPEN (NEWUNIT=unit_nr, FILE=filename, STATUS="OLD", ACTION="READWRITE", IOSTAT=ios)
      IF (ios == 0) CLOSE (unit_nr, STATUS="DELETE")

   END SUBROUTINE delete_screening_file

! **************************************************************************************************
!> \brief Join a screening working directory and filename.
! **************************************************************************************************
   PURE FUNCTION screening_path(directory, filename) RESULT(path)
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

   END FUNCTION screening_path

! **************************************************************************************************
!> \brief driver routine for the parallel part of the method
!> \param opt_bas ...
!> \param input_declaration ...
!> \param para_env ...
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE driver_para_opt_basis(opt_bas, input_declaration, para_env)
      TYPE(basis_optimization_type)                      :: opt_bas
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'driver_para_opt_basis'

      INTEGER                                            :: handle, n_groups_created
      TYPE(mp_comm_type)                                 :: opt_group
      INTEGER, DIMENSION(:), POINTER                     :: group_distribution_p
      INTEGER, DIMENSION(0:para_env%num_pe-1), TARGET    :: group_distribution

      CALL timeset(routineN, handle)
      group_distribution_p => group_distribution
      CALL opt_group%from_split(para_env, n_groups_created, group_distribution_p, &
                                n_subgroups=SIZE(opt_bas%group_partition), group_partition=opt_bas%group_partition)
      opt_bas%opt_id = group_distribution(para_env%mepos) + 1
      opt_bas%n_groups_created = n_groups_created
      ALLOCATE (opt_bas%sub_sources(0:para_env%num_pe - 1))

      CALL driver_optimization_para_low(opt_bas, input_declaration, para_env, opt_group)

      CALL opt_group%free()
      CALL timestop(handle)

   END SUBROUTINE driver_para_opt_basis

! **************************************************************************************************
!> \brief low level optimization routine includes initialization of the subsytems
!>        powell optimizer and deallocation of the various force envs
!> \param opt_bas ...
!> \param input_declaration ...
!> \param para_env_top ...
!> \param mpi_comm_opt ...
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE driver_optimization_para_low(opt_bas, input_declaration, para_env_top, mpi_comm_opt)
      TYPE(basis_optimization_type)                      :: opt_bas
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(mp_para_env_type), POINTER                    :: para_env_top
      TYPE(mp_comm_type), INTENT(IN)                     :: mpi_comm_opt

      CHARACTER(len=*), PARAMETER :: routineN = 'driver_optimization_para_low'

      INTEGER                                            :: handle, icalc, iopt, is, mp_id, stat
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: f_env_id, f_env_output_unit
      LOGICAL                                            :: write_basis
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tot_time
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: diagnostics_initial, &
                                                            diagnostics_optimized
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: matrix_S_inv
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(frontier_orbitals_reference_type), &
         ALLOCATABLE, DIMENSION(:)                       :: frontier_orbitals_reference
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (f_env)

      CALL timeset(routineN, handle)

      ! ======  initialize the f_env and precompute some matrices =====
      mp_id = opt_bas%opt_id
      NULLIFY (para_env, f_env)
      ALLOCATE (f_env_id(SIZE(opt_bas%comp_group(mp_id)%member_list)))
      ALLOCATE (f_env_output_unit(SIZE(opt_bas%comp_group(mp_id)%member_list)), source=-1)
      ALLOCATE (tot_time(opt_bas%ncombinations*opt_bas%ntraining_sets))
      ALLOCATE (diagnostics_initial(opt_bas%ncombinations*opt_bas%ntraining_sets, &
                                    number_frontier_diagnostics), source=0.0_dp)
      ALLOCATE (diagnostics_optimized(opt_bas%ncombinations*opt_bas%ntraining_sets, &
                                      number_frontier_diagnostics), source=0.0_dp)
      ALLOCATE (matrix_s_inv(SIZE(opt_bas%comp_group(mp_id)%member_list)))
      ALLOCATE (frontier_orbitals_reference(SIZE(opt_bas%comp_group(mp_id)%member_list)))

      ALLOCATE (para_env)
      para_env = mpi_comm_opt

      is = -1
      IF (para_env%is_source()) is = para_env_top%mepos
      CALL para_env_top%allgather(is, opt_bas%sub_sources)

      CALL init_training_force_envs(opt_bas, f_env_id, f_env_output_unit, input_declaration, matrix_s_inv, &
                                    frontier_orbitals_reference, para_env, mpi_comm_opt)

      CALL init_free_vars(opt_bas)
      IF (opt_bas%method == method_mo_fit_occ_virtual) THEN
         ! Store the normalized coefficient representation used by the basis code.
         CALL update_free_vars(opt_bas)
         CALL init_free_vars(opt_bas)
         ALLOCATE (opt_bas%x_initial, source=opt_bas%x_opt)
      END IF
      tot_time = 0.0_dp

      ! ======= The real optimization loop  =======
      DO iopt = 0, opt_bas%powell_param%maxfun
         SELECT CASE (opt_bas%method)
         CASE (method_mo_fit_occ)
            CALL compute_residuum_vectors(opt_bas, f_env_id, matrix_S_inv, tot_time, &
                                          para_env_top, para_env, iopt)
         CASE (method_mo_fit_occ_virtual)
            IF (iopt == 0 .AND. .NOT. opt_bas%quiet_output) THEN
               CALL compute_frontier_orbitals_objective( &
                  opt_bas, f_env_id, frontier_orbitals_reference, tot_time, &
                  para_env_top, para_env, iopt, diagnostics_initial)
            ELSE
               CALL compute_frontier_orbitals_objective( &
                  opt_bas, f_env_id, frontier_orbitals_reference, tot_time, &
                  para_env_top, para_env, iopt)
            END IF
         CASE DEFAULT
            CPABORT("Unknown basis optimization method")
         END SELECT
         IF (para_env_top%is_source()) THEN
            CALL powell_optimize(opt_bas%powell_param%nvar, opt_bas%x_opt, opt_bas%powell_param)
         END IF
         CALL para_env_top%bcast(opt_bas%powell_param%state)
         CALL para_env_top%bcast(opt_bas%x_opt)
         CALL update_free_vars(opt_bas)
         write_basis = MOD(iopt, opt_bas%write_frequency) == 0
         CALL update_derived_basis_sets(opt_bas, write_basis, opt_bas%output_basis_file, &
                                        para_env_top)
         IF (opt_bas%powell_param%state == -1) EXIT
      END DO

      ! ======= Update the basis set and print the final basis  =======
      IF (para_env_top%is_source()) THEN
         opt_bas%powell_param%state = 8
         CALL powell_optimize(opt_bas%powell_param%nvar, opt_bas%x_opt, opt_bas%powell_param)
      END IF

      CALL para_env_top%bcast(opt_bas%x_opt)
      CALL update_free_vars(opt_bas)
      IF (opt_bas%method == method_mo_fit_occ_virtual .AND. .NOT. opt_bas%quiet_output) THEN
         ! Report the actual best Powell vector, rather than the final trial point.
         CALL compute_frontier_orbitals_objective( &
            opt_bas, f_env_id, frontier_orbitals_reference, tot_time, &
            para_env_top, para_env, opt_bas%powell_param%maxfun + 1, diagnostics_optimized)
         CALL print_frontier_orbitals_diagnostics( &
            opt_bas, frontier_orbitals_reference, diagnostics_initial, diagnostics_optimized, &
            para_env_top, para_env)
      END IF
      CALL update_derived_basis_sets(opt_bas, .TRUE., opt_bas%output_basis_file, &
                                     para_env_top)

      ! ======  get rid of the f_env again =====

      DO icalc = SIZE(opt_bas%comp_group(mp_id)%member_list), 1, -1
         IF (opt_bas%method == method_mo_fit_occ_virtual) THEN
            CALL frontier_orbitals_reference_release(frontier_orbitals_reference(icalc))
         END IF
         CALL f_env_get_from_id(f_env_id(icalc), f_env)
         CALL destroy_force_env(f_env_id(icalc), stat)
         IF (opt_bas%quiet_output .AND. para_env%is_source()) CLOSE (f_env_output_unit(icalc))
      END DO
      DEALLOCATE (f_env_id, f_env_output_unit); DEALLOCATE (tot_time)
      DEALLOCATE (diagnostics_initial, diagnostics_optimized)
      DEALLOCATE (frontier_orbitals_reference)
      CALL cp_fm_release(matrix_s_inv)
      CALL mp_para_env_release(para_env)
      CALL timestop(handle)

   END SUBROUTINE driver_optimization_para_low

! **************************************************************************************************
!> \brief compute all ingredients for powell optimizer. Rho_diff,
!>        condition number, energy,... for all ttraining sets in
!>        the computational group
!> \param opt_bas ...
!> \param f_env_id ...
!> \param matrix_S_inv ...
!> \param tot_time ...
!> \param para_env_top ...
!> \param para_env ...
!> \param iopt ...
! **************************************************************************************************

   SUBROUTINE compute_residuum_vectors(opt_bas, f_env_id, matrix_S_inv, tot_time, &
                                       para_env_top, para_env, iopt)
      TYPE(basis_optimization_type)                      :: opt_bas
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: f_env_id
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: matrix_S_inv
      REAL(KIND=dp), DIMENSION(:)                        :: tot_time
      TYPE(mp_para_env_type), POINTER                    :: para_env_top, para_env
      INTEGER                                            :: iopt

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_residuum_vectors'

      CHARACTER(len=8)                                   :: basis_type
      INTEGER                                            :: bas_id, handle, icalc, icomb, ispin, &
                                                            mp_id, my_id, nao, ncalc, nelectron, &
                                                            nmo, nspins, set_id
      REAL(KIND=dp)                                      :: flexible_electron_count, maxocc, n_el_f
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cond_vec, energy, f_vec, my_time, &
                                                            start_time
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gdata
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s_aux, matrix_s_aux_orb
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(mo_set_type), ALLOCATABLE, DIMENSION(:)       :: mos_aux
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux, sab_aux_orb
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      basis_type = "AUX_OPT"
      !
      ncalc = opt_bas%ncombinations*opt_bas%ntraining_sets
      ALLOCATE (gdata(ncalc, 4))
      f_vec => gdata(:, 1)
      my_time => gdata(:, 2)
      cond_vec => gdata(:, 3)
      energy => gdata(:, 4)
      !
      f_vec = 0.0_dp; cond_vec = 0.0_dp; my_time = 0.0_dp; energy = 0.0_dp
      mp_id = opt_bas%opt_id
      ALLOCATE (start_time(SIZE(opt_bas%comp_group(mp_id)%member_list)))
      !
      DO icalc = 1, SIZE(opt_bas%comp_group(mp_id)%member_list)
         my_id = opt_bas%comp_group(mp_id)%member_list(icalc) + 1
         ! setup timings
         start_time(icalc) = m_walltime()

         NULLIFY (matrix_s_aux_orb, matrix_s_aux)
         CALL get_set_and_basis_id(opt_bas%comp_group(mp_id)%member_list(icalc), opt_bas, set_id, bas_id)
         CALL f_env_get_from_id(f_env_id(icalc), f_env)
         force_env => f_env%force_env
         CALL force_env_get(force_env, qs_env=qs_env)
         CALL get_qs_env(qs_env, ks_env=ks_env)
         CALL update_basis_set(opt_bas, bas_id, basis_type, qs_env)
         NULLIFY (sab_aux, sab_aux_orb)
         CALL optbas_build_neighborlist(qs_env, sab_aux, sab_aux_orb, basis_type)
         CALL build_overlap_matrix(ks_env, matrix_s=matrix_s_aux, &
                                   basis_type_a=basis_type, &
                                   basis_type_b=basis_type, &
                                   sab_nl=sab_aux)
         CALL build_overlap_matrix(ks_env, matrix_s=matrix_s_aux_orb, &
                                   basis_type_a=basis_type, &
                                   basis_type_b="ORB", &
                                   sab_nl=sab_aux_orb)
         CALL release_neighbor_list_sets(sab_aux)
         CALL release_neighbor_list_sets(sab_aux_orb)
         CALL get_qs_env(qs_env, mos=mos, matrix_ks=matrix_ks)

         nspins = SIZE(mos)
         ALLOCATE (mos_aux(nspins))
         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
         CALL get_qs_kind_set(qs_kind_set, nsgf=nao, basis_type=basis_type)
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, maxocc=maxocc, nelectron=nelectron, &
                            n_el_f=n_el_f, nmo=nmo, flexible_electron_count=flexible_electron_count)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo, &
                                     context=mo_coeff%matrix_struct%context, &
                                     para_env=mo_coeff%matrix_struct%para_env)
            CALL allocate_mo_set(mos_aux(ispin), nao, nmo, nelectron, &
                                 n_el_f, maxocc, flexible_electron_count)
            CALL init_mo_set(mo_set=mos_aux(ispin), fm_struct=fm_struct, name="MO_AUX")
            CALL cp_fm_struct_release(fm_struct)
         END DO

         CALL fit_mo_coeffs(matrix_s_aux, matrix_s_aux_orb, mos, mos_aux)
         CALL evaluate_optvals(mos, mos_aux, matrix_ks, matrix_s_aux_orb(1)%matrix, &
                               matrix_s_aux(1)%matrix, matrix_S_inv(icalc), &
                               f_vec(my_id), energy(my_id), cond_vec(my_id))

         DO ispin = 1, nspins
            CALL deallocate_mo_set(mos_aux(ispin))
         END DO
         DEALLOCATE (mos_aux)
         IF (ASSOCIATED(matrix_s_aux)) CALL dbcsr_deallocate_matrix_set(matrix_s_aux)
         IF (ASSOCIATED(matrix_s_aux_orb)) CALL dbcsr_deallocate_matrix_set(matrix_s_aux_orb)

         my_time(my_id) = m_walltime() - start_time(icalc)
      END DO

      DEALLOCATE (start_time)

      IF (.NOT. para_env%is_source()) THEN
         f_vec = 0.0_dp; cond_vec = 0.0_dp; my_time = 0.0_dp; energy = 0.0_dp
      END IF
      ! collect date from all subgroup ionodes on the main ionode
      CALL para_env_top%sum(gdata)

      opt_bas%powell_param%f = 0.0_dp
      IF (para_env_top%is_source()) THEN
         DO icalc = 1, SIZE(f_vec)
            icomb = MOD(icalc - 1, opt_bas%ncombinations)
            opt_bas%powell_param%f = opt_bas%powell_param%f + &
                                     (f_vec(icalc) + energy(icalc))*opt_bas%fval_weight(icomb)
            IF (opt_bas%use_condition_number) THEN
               opt_bas%powell_param%f = opt_bas%powell_param%f + &
                                        LOG(cond_vec(icalc))*opt_bas%condition_weight(icomb)
            END IF
         END DO
      ELSE
         f_vec = 0.0_dp; cond_vec = 0.0_dp; my_time = 0.0_dp; energy = 0.0_dp
      END IF
      CALL para_env_top%bcast(opt_bas%powell_param%f)

      ! output info if required
      CALL output_opt_info(f_vec, cond_vec, my_time, tot_time, opt_bas, iopt, para_env_top)
      DEALLOCATE (gdata)

      CALL para_env_top%sync()

      CALL timestop(handle)

   END SUBROUTINE compute_residuum_vectors

! **************************************************************************************************
!> \brief Evaluate the frozen-density frontier-orbital objective.
!> \param opt_bas basis optimization settings and current parameters
!> \param f_env_id force environments for locally assigned calculations
!> \param reference cached full reference eigensystems
!> \param tot_time accumulated timings for output
!> \param para_env_top top-level communicator
!> \param para_env calculation-group communicator
!> \param iopt Powell evaluation index
! **************************************************************************************************
   SUBROUTINE compute_frontier_orbitals_objective(opt_bas, f_env_id, reference, tot_time, &
                                                  para_env_top, para_env, iopt, diagnostics)
      TYPE(basis_optimization_type)                      :: opt_bas
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: f_env_id
      TYPE(frontier_orbitals_reference_type), &
         DIMENSION(:), INTENT(IN)                        :: reference
      REAL(KIND=dp), DIMENSION(:)                        :: tot_time
      TYPE(mp_para_env_type), POINTER                    :: para_env_top, para_env
      INTEGER                                            :: iopt
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: diagnostics

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_frontier_orbitals_objective'

      CHARACTER(len=8)                                   :: basis_type
      INTEGER                                            :: bas_id, handle, icalc, icomb, icont, &
                                                            ikind, ipgf, iset, ix, l, mp_id, &
                                                            my_id, ncalc, set_id
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: candidate_angular_momentum
      REAL(KIND=dp)                                      :: coefficient_loss, contribution
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coefficients_current
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gdata
      REAL(KIND=dp), DIMENSION(:), POINTER               :: start_time
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s_aux, matrix_s_aux_orb
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(frontier_orbitals_result_type)                :: result
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux, sab_aux_orb
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      basis_type = "AUX_OPT"
      ncalc = opt_bas%ncombinations*opt_bas%ntraining_sets
      ALLOCATE (gdata(ncalc, number_frontier_diagnostics), source=0.0_dp)
      ALLOCATE (coefficients_current, source=opt_bas%x_opt)
      DO ikind = 1, opt_bas%nkind
         DO iset = 1, opt_bas%kind_basis(ikind)%flex_basis(0)%nsets
            DO ipgf = 1, opt_bas%kind_basis(ikind)%flex_basis(0)%subset(iset)%nexp
               DO icont = 1, opt_bas%kind_basis(ikind)%flex_basis(0)%subset(iset)%ncon_tot
                  IF (opt_bas%kind_basis(ikind)%flex_basis(0)%subset(iset)%opt_coeff(ipgf, icont)) THEN
                     ix = opt_bas%kind_basis(ikind)%flex_basis(0)%subset(iset)%coeff_x_ind(ipgf, icont)
                     coefficients_current(ix) = &
                        opt_bas%kind_basis(ikind)%flex_basis(0)%subset(iset)%coeff(ipgf, icont)
                  END IF
               END DO
            END DO
         END DO
      END DO
      coefficient_loss = 0.0_dp
      IF (SUM(opt_bas%x_initial**2) > TINY(1.0_dp)) THEN
         coefficient_loss = SUM((coefficients_current - opt_bas%x_initial)**2)/ &
                            SUM(opt_bas%x_initial**2)
      END IF
      mp_id = opt_bas%opt_id
      ALLOCATE (start_time(SIZE(opt_bas%comp_group(mp_id)%member_list)))

      DO icalc = 1, SIZE(opt_bas%comp_group(mp_id)%member_list)
         my_id = opt_bas%comp_group(mp_id)%member_list(icalc) + 1
         start_time(icalc) = m_walltime()

         NULLIFY (ks_env, matrix_s_aux, matrix_s_aux_orb, sab_aux, sab_aux_orb)
         CALL get_set_and_basis_id(opt_bas%comp_group(mp_id)%member_list(icalc), opt_bas, &
                                   bas_id=bas_id, set_id=set_id)
         CALL f_env_get_from_id(f_env_id(icalc), f_env)
         force_env => f_env%force_env
         CALL force_env_get(force_env, qs_env=qs_env)
         CALL update_basis_set(opt_bas, bas_id, basis_type, qs_env)
         CALL optbas_build_neighborlist(qs_env, sab_aux, sab_aux_orb, basis_type)
         CALL get_qs_env(qs_env, ks_env=ks_env)
         CALL build_overlap_matrix(ks_env, matrix_s=matrix_s_aux, &
                                   basis_type_a=basis_type, basis_type_b=basis_type, sab_nl=sab_aux)
         CALL build_overlap_matrix(ks_env, matrix_s=matrix_s_aux_orb, &
                                   basis_type_a=basis_type, basis_type_b="ORB", sab_nl=sab_aux_orb)
         CALL release_neighbor_list_sets(sab_aux)
         CALL release_neighbor_list_sets(sab_aux_orb)

         IF (PRESENT(diagnostics)) THEN
            CALL get_ao_angular_momenta(qs_env, basis_type, candidate_angular_momentum)
            CALL evaluate_frontier_orbitals_objective( &
               reference(icalc), matrix_s_aux(1)%matrix, matrix_s_aux_orb(1)%matrix, &
               opt_bas%virtual_energy_cutoff, opt_bas%virtual_energy_smoothing, &
               opt_bas%gap_energy_scale, RESULT, candidate_angular_momentum)
            DEALLOCATE (candidate_angular_momentum)
         ELSE
            CALL evaluate_frontier_orbitals_objective( &
               reference(icalc), matrix_s_aux(1)%matrix, matrix_s_aux_orb(1)%matrix, &
               opt_bas%virtual_energy_cutoff, opt_bas%virtual_energy_smoothing, &
               opt_bas%gap_energy_scale, RESULT)
         END IF
         gdata(my_id, 1) = result%loss_occupied_subspace
         gdata(my_id, 2) = result%loss_virtual_subspace
         gdata(my_id, 3) = result%loss_gap
         gdata(my_id, 4) = result%occupied_subspace_overlap
         gdata(my_id, 5) = result%gap_reference
         gdata(my_id, 6) = result%gap_candidate
         gdata(my_id, 7) = result%condition_number
         gdata(my_id, 8) = result%empty_subspace_overlap
         gdata(my_id, 9) = m_walltime() - start_time(icalc)
         gdata(my_id, 10) = result%loss_empty_subspace
         gdata(my_id, frontier_idx_loss_coefficient) = coefficient_loss
         gdata(my_id, frontier_idx_similarity_virtual) = result%virtual_subspace_similarity
         gdata(my_id, frontier_idx_condition_reference) = reference(icalc)%condition_number
         gdata(my_id, frontier_idx_nao_reference) = REAL(reference(icalc)%number_reference_orbitals, dp)
         gdata(my_id, frontier_idx_nao_candidate) = REAL(result%number_candidate_orbitals, dp)
         IF (PRESENT(diagnostics)) THEN
            DO l = 0, max_diagnostic_angular_momentum
               gdata(my_id, frontier_idx_reference_count + l) = &
                  REAL(result%angular_momentum%reference_function_count(l), dp)
               gdata(my_id, frontier_idx_candidate_count + l) = &
                  REAL(result%angular_momentum%candidate_function_count(l), dp)
            END DO
         END IF

         IF (ASSOCIATED(matrix_s_aux)) CALL dbcsr_deallocate_matrix_set(matrix_s_aux)
         IF (ASSOCIATED(matrix_s_aux_orb)) CALL dbcsr_deallocate_matrix_set(matrix_s_aux_orb)
      END DO
      DEALLOCATE (start_time)

      IF (.NOT. para_env%is_source()) gdata = 0.0_dp
      CALL para_env_top%sum(gdata)
      IF (PRESENT(diagnostics)) THEN
         CPASSERT(SIZE(diagnostics, 1) == ncalc)
         CPASSERT(SIZE(diagnostics, 2) == number_frontier_diagnostics)
         diagnostics = gdata
      END IF

      opt_bas%powell_param%f = 0.0_dp
      IF (para_env_top%is_source()) THEN
         DO icalc = 1, ncalc
            icomb = MOD(icalc - 1, opt_bas%ncombinations)
            contribution = opt_bas%occupied_weight*gdata(icalc, 1) + &
                           opt_bas%virtual_weight*gdata(icalc, 2) + &
                           opt_bas%empty_overlap_weight*gdata(icalc, 10) + &
                           opt_bas%gap_weight*gdata(icalc, 3) + &
                           opt_bas%coefficient_weight*coefficient_loss
            IF (opt_bas%use_condition_number) THEN
               contribution = contribution + opt_bas%condition_weight(icomb)*LOG10(gdata(icalc, 7))
            END IF
            opt_bas%powell_param%f = opt_bas%powell_param%f + &
                                     opt_bas%fval_weight(icomb)*contribution
            tot_time(icalc) = tot_time(icalc) + gdata(icalc, 9)
         END DO

      END IF
      CALL para_env_top%bcast(opt_bas%powell_param%f)
      DEALLOCATE (coefficients_current, gdata)

      CALL para_env_top%sync()
      CALL timestop(handle)

   END SUBROUTINE compute_frontier_orbitals_objective

! **************************************************************************************************
!> \brief Print initial-to-optimized frontier-orbital diagnostics for every training calculation.
!> \param opt_bas basis optimization settings
!> \param reference locally held reference data and basis metadata
!> \param initial initial-small-basis diagnostics
!> \param optimized optimized-small-basis diagnostics
!> \param para_env_top top-level communicator
!> \param para_env calculation-group communicator
! **************************************************************************************************
   SUBROUTINE print_frontier_orbitals_diagnostics(opt_bas, reference, initial, optimized, &
                                                  para_env_top, para_env)
      TYPE(basis_optimization_type), INTENT(IN)          :: opt_bas
      TYPE(frontier_orbitals_reference_type), &
         DIMENSION(:), INTENT(IN)                        :: reference
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: initial, optimized
      TYPE(mp_para_env_type), POINTER                    :: para_env_top, para_env

      CHARACTER(LEN=17)                                  :: basis_label
      CHARACTER(LEN=default_string_length)               :: candidate_basis, element, reference_basis
      INTEGER                                            :: bas_id, icalc, ikind, iset, l, &
                                                            max_kinds, my_id, ncalc, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_count
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: metadata
      TYPE(cp_logger_type), POINTER                      :: logger

      ncalc = opt_bas%ncombinations*opt_bas%ntraining_sets
      CPASSERT(SIZE(initial, 1) == ncalc)
      CPASSERT(SIZE(optimized, 1) == ncalc)

      max_kinds = 0
      DO icalc = 1, SIZE(reference)
         max_kinds = MAX(max_kinds, reference(icalc)%number_kinds)
      END DO
      CALL para_env_top%max(max_kinds)
      ALLOCATE (kind_count(ncalc), source=0)
      ALLOCATE (metadata(ncalc, max_kinds, 3, default_string_length), source=0)

      IF (para_env%is_source()) THEN
         DO icalc = 1, SIZE(reference)
            my_id = opt_bas%comp_group(opt_bas%opt_id)%member_list(icalc) + 1
            kind_count(my_id) = reference(icalc)%number_kinds
            DO ikind = 1, reference(icalc)%number_kinds
               CALL encode_string(reference(icalc)%element_symbol(ikind), metadata(my_id, ikind, 1, :))
               CALL encode_string(reference(icalc)%reference_basis_name(ikind), metadata(my_id, ikind, 2, :))
               CALL encode_string(reference(icalc)%candidate_basis_name(ikind), metadata(my_id, ikind, 3, :))
            END DO
         END DO
      END IF
      CALL para_env_top%sum(kind_count)
      CALL para_env_top%sum(metadata)

      logger => cp_get_default_logger()
      unit_nr = -1
      IF (para_env_top%is_source()) unit_nr = cp_logger_get_default_unit_nr(logger)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(A)') ""
         WRITE (unit_nr, '(A)') REPEAT("=", 78)
         WRITE (unit_nr, '(A)') "BASIS OPTIMIZATION SUMMARY"
         WRITE (unit_nr, '(A)') REPEAT("=", 78)

         DO my_id = 1, ncalc
            CALL get_set_and_basis_id(my_id - 1, opt_bas, set_id=iset, bas_id=bas_id)
            WRITE (unit_nr, '(A)') ""
            WRITE (unit_nr, '(A,A,A,A)') "Reference calculation: ", &
               TRIM(opt_bas%training_dir(iset)), "/", TRIM(opt_bas%training_input(iset))
            IF (opt_bas%ncombinations > 1) THEN
               WRITE (unit_nr, '(A,I0)') "Basis combination: ", bas_id
            END IF
            WRITE (unit_nr, '(A)') ""
            WRITE (unit_nr, '(A)') "Atomic kind basis sets"
            WRITE (unit_nr, '(A)') REPEAT("-", 78)
            DO ikind = 1, kind_count(my_id)
               CALL decode_string(metadata(my_id, ikind, 1, :), element)
               CALL decode_string(metadata(my_id, ikind, 2, :), reference_basis)
               CALL decode_string(metadata(my_id, ikind, 3, :), candidate_basis)
               WRITE (unit_nr, '(A,4X,A,A)') TRIM(element), "Reference basis:       ", &
                  TRIM(reference_basis)
               WRITE (unit_nr, '(A,4X,A,A)') REPEAT(" ", LEN_TRIM(element)), &
                  "Initial small basis:   ", TRIM(candidate_basis)
               WRITE (unit_nr, '(A,4X,A,A)') REPEAT(" ", LEN_TRIM(element)), &
                  "Optimized small basis: ", TRIM(candidate_basis)
            END DO
            WRITE (unit_nr, '(A)') ""
            WRITE (unit_nr, '(A,A)') "Optimized basis file: ", TRIM(opt_bas%output_basis_file)
            WRITE (unit_nr, '(A)') ""

            WRITE (unit_nr, '(A)') "Overall subspace quality and conditioning"
            WRITE (unit_nr, '(A)') REPEAT("-", 78)
            basis_label = "Basis"
            WRITE (unit_nr, '(A17,1X,A4,2(1X,A11),1X,A12,1X,A9,1X,A6)') basis_label, "N_AO", &
               "Occ.overlap", "Vir.overlap", "Emp.coverage", "κ(S)", "Gap/eV"
            WRITE (unit_nr, '(A)') REPEAT("-", 78)
            basis_label = "Ref. basis"
            WRITE (unit_nr, '(A17,1X,I4,2(1X,F10.1,A),1X,F11.1,A,1X,ES9.2,1X,F6.3)') &
               basis_label, NINT(initial(my_id, frontier_idx_nao_reference)), &
               100.0_dp, "%", 100.0_dp, "%", 100.0_dp, "%", &
               initial(my_id, frontier_idx_condition_reference), &
               initial(my_id, frontier_idx_gap_reference)*evolt
            basis_label = "Init. small basis"
            WRITE (unit_nr, '(A17,1X,I4,2(1X,F10.1,A),1X,F11.1,A,1X,ES9.2,1X,F6.3)') &
               basis_label, NINT(initial(my_id, frontier_idx_nao_candidate)), &
               100.0_dp*initial(my_id, frontier_idx_overlap_occupied), "%", &
               100.0_dp*initial(my_id, frontier_idx_similarity_virtual), "%", &
               100.0_dp*initial(my_id, frontier_idx_overlap_empty), "%", &
               initial(my_id, frontier_idx_condition_candidate), &
               initial(my_id, frontier_idx_gap_candidate)*evolt
            basis_label = "Opt. small basis"
            WRITE (unit_nr, '(A17,1X,I4,2(1X,F10.1,A),1X,F11.1,A,1X,ES9.2,1X,F6.3)') &
               basis_label, NINT(optimized(my_id, frontier_idx_nao_candidate)), &
               100.0_dp*optimized(my_id, frontier_idx_overlap_occupied), "%", &
               100.0_dp*optimized(my_id, frontier_idx_similarity_virtual), "%", &
               100.0_dp*optimized(my_id, frontier_idx_overlap_empty), "%", &
               optimized(my_id, frontier_idx_condition_candidate), &
               optimized(my_id, frontier_idx_gap_candidate)*evolt
            WRITE (unit_nr, '(A)') REPEAT("-", 78)

            WRITE (unit_nr, '(A)') ""
            WRITE (unit_nr, '(A)') "Spherical basis-function content"
            WRITE (unit_nr, '(A)') REPEAT("-", 78)
            basis_label = "Basis"
            WRITE (unit_nr, '(A17)', ADVANCE='NO') basis_label
            DO l = 0, max_diagnostic_angular_momentum
               IF (NINT(initial(my_id, frontier_idx_reference_count + l)) > 0 .OR. &
                   NINT(initial(my_id, frontier_idx_candidate_count + l)) > 0) THEN
                  WRITE (unit_nr, '(1X,A4)', ADVANCE='NO') angular_momentum_label(l)
               END IF
            END DO
            WRITE (unit_nr, '(1X,A5)') "Total"
            basis_label = "Ref. basis"
            WRITE (unit_nr, '(A17)', ADVANCE='NO') basis_label
            DO l = 0, max_diagnostic_angular_momentum
               IF (NINT(initial(my_id, frontier_idx_reference_count + l)) > 0 .OR. &
                   NINT(initial(my_id, frontier_idx_candidate_count + l)) > 0) THEN
                  WRITE (unit_nr, '(1X,I4)', ADVANCE='NO') NINT(initial(my_id, frontier_idx_reference_count + l))
               END IF
            END DO
            WRITE (unit_nr, '(1X,I5)') NINT(initial(my_id, frontier_idx_nao_reference))
            basis_label = "Init. small basis"
            WRITE (unit_nr, '(A17)', ADVANCE='NO') basis_label
            DO l = 0, max_diagnostic_angular_momentum
               IF (NINT(initial(my_id, frontier_idx_reference_count + l)) > 0 .OR. &
                   NINT(initial(my_id, frontier_idx_candidate_count + l)) > 0) THEN
                  WRITE (unit_nr, '(1X,I4)', ADVANCE='NO') NINT(initial(my_id, frontier_idx_candidate_count + l))
               END IF
            END DO
            WRITE (unit_nr, '(1X,I5)') NINT(initial(my_id, frontier_idx_nao_candidate))
            basis_label = "Opt. small basis"
            WRITE (unit_nr, '(A17)', ADVANCE='NO') basis_label
            DO l = 0, max_diagnostic_angular_momentum
               IF (NINT(initial(my_id, frontier_idx_reference_count + l)) > 0 .OR. &
                   NINT(initial(my_id, frontier_idx_candidate_count + l)) > 0) THEN
                  WRITE (unit_nr, '(1X,I4)', ADVANCE='NO') NINT(optimized(my_id, frontier_idx_candidate_count + l))
               END IF
            END DO
            WRITE (unit_nr, '(1X,I5)') NINT(optimized(my_id, frontier_idx_nao_candidate))
            WRITE (unit_nr, '(A)') REPEAT("-", 78)
         END DO

         WRITE (unit_nr, '(A)') REPEAT("=", 78)
         WRITE (unit_nr, '(A)') ""
      END IF

      DEALLOCATE (kind_count, metadata)

   END SUBROUTINE print_frontier_orbitals_diagnostics

! **************************************************************************************************
!> \brief Encode a string as integers so basis metadata can be reduced across optimizer groups.
! **************************************************************************************************
   PURE SUBROUTINE encode_string(string, encoded)
      CHARACTER(LEN=*), INTENT(IN)                       :: string
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: encoded

      INTEGER                                            :: i

      encoded = 0
      DO i = 1, MIN(LEN_TRIM(string), SIZE(encoded))
         encoded(i) = IACHAR(string(i:i))
      END DO

   END SUBROUTINE encode_string

! **************************************************************************************************
!> \brief Decode integer character codes used to communicate basis metadata.
! **************************************************************************************************
   PURE SUBROUTINE decode_string(encoded, string)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: encoded
      CHARACTER(LEN=*), INTENT(OUT)                      :: string

      INTEGER                                            :: i

      string = ""
      DO i = 1, MIN(LEN(string), SIZE(encoded))
         IF (encoded(i) == 0) EXIT
         string(i:i) = ACHAR(encoded(i))
      END DO

   END SUBROUTINE decode_string

! **************************************************************************************************
!> \brief Return the spectroscopic label of an angular-momentum channel.
! **************************************************************************************************
   PURE FUNCTION angular_momentum_label(l) RESULT(label)
      INTEGER, INTENT(IN)                                :: l
      CHARACTER(LEN=7)                                   :: label

      CHARACTER(LEN=*), PARAMETER                        :: labels = "spdfghiklmn"

      label = "unknown"
      IF (l >= 0 .AND. l < LEN(labels)) label = labels(l + 1:l + 1)

   END FUNCTION angular_momentum_label

! **************************************************************************************************
!> \brief create the force_envs for every input in the computational group
!> \param opt_bas ...
!> \param f_env_id ...
!> \param input_declaration ...
!> \param matrix_s_inv ...
!> \param para_env ...
!> \param mpi_comm_opt ...
! **************************************************************************************************

   SUBROUTINE init_training_force_envs(opt_bas, f_env_id, f_env_output_unit, input_declaration, matrix_s_inv, &
                                       frontier_orbitals_reference, para_env, mpi_comm_opt)

      TYPE(basis_optimization_type)                      :: opt_bas
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: f_env_id
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: f_env_output_unit
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(cp_fm_type), DIMENSION(:), INTENT(OUT)        :: matrix_S_inv
      TYPE(frontier_orbitals_reference_type), &
         DIMENSION(:), INTENT(OUT)                       :: frontier_orbitals_reference
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(mp_comm_type)                                 :: mpi_comm_opt

      CHARACTER(len=*), PARAMETER :: routineN = 'init_training_force_envs'

      CHARACTER(len=default_path_length)                 :: main_dir
      INTEGER                                            :: bas_id, handle, icalc, ierr, mp_id, &
                                                            set_id, stat
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: input_file

      CALL timeset(routineN, handle)

      NULLIFY (matrix_s, blacs_env, ks_env)

      mp_id = opt_bas%opt_id
      CALL m_getcwd(main_dir)

      ! ======= Create f_env for all calculations in MPI group =======
      DO icalc = 1, SIZE(opt_bas%comp_group(mp_id)%member_list)
         NULLIFY (input_file)
         ! parse the input of the training sets
         CALL get_set_and_basis_id(opt_bas%comp_group(mp_id)%member_list(icalc), opt_bas, set_id, bas_id)
         CALL m_chdir(TRIM(opt_bas%training_dir(set_id)), ierr)
         IF (ierr /= 0) THEN
            CALL cp_abort(__LOCATION__, &
                          "Could not change to directory <"//TRIM(opt_bas%training_dir(set_id))//">")
         END IF
         input_file => read_input(input_declaration, &
                                  opt_bas%training_input(set_id), &
                                  initial_variables=empty_initial_variables, &
                                  para_env=para_env)

         CALL modify_input_settings(opt_bas, bas_id, input_file)
         IF (opt_bas%quiet_output) THEN
            IF (para_env%is_source()) THEN
               f_env_output_unit(icalc) = get_unit_number()
               IF (f_env_output_unit(icalc) < 1) THEN
                  CPABORT("No logical unit is available for disposable optimizer output")
               END IF
               OPEN (UNIT=f_env_output_unit(icalc), STATUS="SCRATCH", ACTION="READWRITE", IOSTAT=stat)
               IF (stat /= 0) CPABORT("Could not open disposable optimizer output")
            END IF
            CALL create_force_env(f_env_id(icalc), &
                                  input_declaration=input_declaration, &
                                  input_path=opt_bas%training_input(set_id), &
                                  input=input_file, &
                                  output_unit=f_env_output_unit(icalc), &
                                  owns_out_unit=.FALSE., &
                                  mpi_comm=mpi_comm_opt, &
                                  ierr=stat)
         ELSE
            CALL create_force_env(f_env_id(icalc), &
                                  input_declaration=input_declaration, &
                                  input_path=opt_bas%training_input(set_id), &
                                  input=input_file, &
                                  output_path="scrap_information", &
                                  mpi_comm=mpi_comm_opt, &
                                  ierr=stat)
         END IF

         ! some weirdness with the default stacks defaults have to be addded to get the
         ! correct default program name this causes trouble with the timer stack if kept
         CALL f_env_add_defaults(f_env_id(icalc), f_env)
         force_env => f_env%force_env
         CALL force_env_get(force_env, qs_env=qs_env)
         CALL allocate_mo_sets(qs_env)
         CALL f_env_rm_defaults(f_env, stat)
         CALL get_qs_env(qs_env, ks_env=ks_env)
         CALL build_qs_neighbor_lists(qs_env, para_env, molecular=.FALSE., &
                                      force_env_section=qs_env%input)
         CALL get_ks_env(ks_env, &
                         matrix_s=matrix_s, &
                         sab_orb=sab_orb)
         CALL build_overlap_matrix(ks_env, matrix_s=matrix_s, &
                                   matrix_name="OVERLAP", &
                                   basis_type_a="ORB", &
                                   basis_type_b="ORB", &
                                   sab_nl=sab_orb)
         CALL set_ks_env(ks_env, matrix_s=matrix_s)
         CALL get_qs_env(qs_env, matrix_s=matrix_s, blacs_env=blacs_env)
         CALL calculate_overlap_inverse(matrix_s(1)%matrix, matrix_s_inv(icalc), &
                                        para_env, blacs_env)
         CALL calculate_ks_matrix(qs_env)
         IF (opt_bas%method == method_mo_fit_occ_virtual) THEN
            CALL frontier_orbitals_reference_init(qs_env, frontier_orbitals_reference(icalc))
         END IF

         CALL section_vals_release(input_file)

         CALL qs_env_part_release(qs_env)

         CALL m_chdir(TRIM(ADJUSTL(main_dir)), ierr)
      END DO

      CALL timestop(handle)

   END SUBROUTINE init_training_force_envs

! **************************************************************************************************
!> \brief variable update from the powell vector for all sets
!> \param opt_bas ...
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE update_free_vars(opt_bas)
      TYPE(basis_optimization_type)                      :: opt_bas

      CHARACTER(len=*), PARAMETER                        :: routineN = 'update_free_vars'

      INTEGER                                            :: handle, ikind, iset, ix

      CALL timeset(routineN, handle)
      ix = 0
      DO ikind = 1, opt_bas%nkind
         DO iset = 1, opt_bas%kind_basis(ikind)%flex_basis(0)%nsets
            CALL update_subset_freevars(opt_bas%kind_basis(ikind)%flex_basis(0)%subset(iset), ix, opt_bas%x_opt)
         END DO
      END DO
      CALL timestop(handle)

   END SUBROUTINE update_free_vars

! **************************************************************************************************
!> \brief low level update for the basis sets. Exponents are transformed according to constraint
!> \param subset ...
!> \param ix ...
!> \param x ...
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE update_subset_freevars(subset, ix, x)
      TYPE(subset_type)                                  :: subset
      INTEGER                                            :: ix
      REAL(KIND=dp), DIMENSION(:)                        :: x

      CHARACTER(len=*), PARAMETER :: routineN = 'update_subset_freevars'

      INTEGER                                            :: handle, icon1, icon2, icont, iexp, il, &
                                                            istart
      REAL(KIND=dp)                                      :: fermi_f, gs_scale

      CALL timeset(routineN, handle)
      DO iexp = 1, subset%nexp
         IF (subset%opt_exps(iexp)) THEN
            ix = ix + 1
            subset%exps(iexp) = ABS(x(ix))
            IF (subset%exp_has_const(iexp)) THEN
               !use a fermi function to keep exponents in a given range around their initial value
               fermi_f = 1.0_dp/(EXP((x(ix) - 1.0_dp)/0.5_dp) + 1.0_dp)
               subset%exps(iexp) = (2.0_dp*fermi_f - 1.0_dp)*subset%exp_const(iexp)%var_fac*subset%exp_const(iexp)%init + &
                                   subset%exp_const(iexp)%init
            ELSE

            END IF
         END IF
         DO icont = 1, subset%ncon_tot
            IF (subset%opt_coeff(iexp, icont)) THEN
               ix = ix + 1
               subset%coeff(iexp, icont) = x(ix)
            END IF
         END DO
      END DO

      ! orthonormalize contraction coefficients using gram schmidt
      istart = 1
      DO il = 1, subset%nl
         DO icon1 = istart, istart + subset%l(il) - 2
            DO icon2 = icon1 + 1, istart + subset%l(il) - 1
               gs_scale = DOT_PRODUCT(subset%coeff(:, icon2), subset%coeff(:, icon1))/ &
                          DOT_PRODUCT(subset%coeff(:, icon1), subset%coeff(:, icon1))
               subset%coeff(:, icon2) = subset%coeff(:, icon2) - gs_scale*subset%coeff(:, icon1)
            END DO
         END DO
         istart = istart + subset%l(il)
      END DO

      DO icon1 = 1, subset%ncon_tot
         subset%coeff(:, icon1) = subset%coeff(:, icon1)/NORM2(subset%coeff(:, icon1))
      END DO
      CALL timestop(handle)

   END SUBROUTINE update_subset_freevars

! **************************************************************************************************
!> \brief variable initialization for the powell vector for all sets
!> \param opt_bas ...
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE init_free_vars(opt_bas)
      TYPE(basis_optimization_type)                      :: opt_bas

      CHARACTER(len=*), PARAMETER                        :: routineN = 'init_free_vars'

      INTEGER                                            :: handle, ikind, iset, ix

      CALL timeset(routineN, handle)
      ix = 0
      DO ikind = 1, opt_bas%nkind
         DO iset = 1, opt_bas%kind_basis(ikind)%flex_basis(0)%nsets
            CALL init_subset_freevars(opt_bas%kind_basis(ikind)%flex_basis(0)%subset(iset), ix, opt_bas%x_opt)
         END DO
      END DO
      CALL timestop(handle)

   END SUBROUTINE init_free_vars

! **************************************************************************************************
!> \brief variable initialization for the powell vector from low level informations
!>        constraint exponents will be mapped on a fermi function
!> \param subset ...
!> \param ix ...
!> \param x ...
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE init_subset_freevars(subset, ix, x)
      TYPE(subset_type)                                  :: subset
      INTEGER                                            :: ix
      REAL(KIND=dp), DIMENSION(:)                        :: x

      CHARACTER(len=*), PARAMETER :: routineN = 'init_subset_freevars'

      INTEGER                                            :: handle, icont, iexp
      REAL(KIND=dp)                                      :: fract

      CALL timeset(routineN, handle)

      DO iexp = 1, subset%nexp
         IF (subset%opt_exps(iexp)) THEN
            ix = ix + 1
            x(ix) = subset%exps(iexp)
            IF (subset%exp_has_const(iexp)) THEN
               IF (subset%exp_const(iexp)%const_type == 0) THEN
                  fract = 1.0_dp + (subset%exps(iexp) - subset%exp_const(iexp)%init)/ &
                          (subset%exp_const(iexp)%init*subset%exp_const(iexp)%var_fac)
                  x(ix) = 0.5_dp*LOG((2.0_dp/fract - 1.0_dp)) + 1.0_dp
               END IF
               IF (subset%exp_const(iexp)%const_type == 1) THEN
                  x(ix) = 1.0_dp
               END IF
            END IF
         END IF
         DO icont = 1, subset%ncon_tot
            IF (subset%opt_coeff(iexp, icont)) THEN
               ix = ix + 1
               x(ix) = subset%coeff(iexp, icont)
            END IF
         END DO
      END DO
      CALL timestop(handle)

   END SUBROUTINE init_subset_freevars

! **************************************************************************************************
!> \brief commuticates all info to the master and assembles the output
!> \param f_vec ...
!> \param cond_vec ...
!> \param my_time ...
!> \param tot_time ...
!> \param opt_bas ...
!> \param iopt ...
!> \param para_env_top ...
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE output_opt_info(f_vec, cond_vec, my_time, tot_time, opt_bas, iopt, para_env_top)
      REAL(KIND=dp), DIMENSION(:)                        :: f_vec, cond_vec, my_time, tot_time
      TYPE(basis_optimization_type)                      :: opt_bas
      INTEGER                                            :: iopt
      TYPE(mp_para_env_type), POINTER                    :: para_env_top

      CHARACTER(len=*), PARAMETER                        :: routineN = 'output_opt_info'

      INTEGER                                            :: handle, ibasis, icalc, iset, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      tot_time = tot_time + my_time

      unit_nr = -1
      IF (para_env_top%is_source() .AND. (MOD(iopt, opt_bas%write_frequency) == 0 .OR. iopt == opt_bas%powell_param%maxfun)) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger)
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(1X,A,I8)') "BASOPT| Information at iteration number:", iopt
         WRITE (unit_nr, '(1X,A)') "BASOPT| Training set | Combination | Rho difference | Condition num. | Time"
         WRITE (unit_nr, '(1X,A)') "BASOPT| -----------------------------------------------------------------------"
         icalc = 0
         DO iset = 1, opt_bas%ntraining_sets
            DO ibasis = 1, opt_bas%ncombinations
               icalc = icalc + 1
               WRITE (unit_nr, '(1X,A,2(5X,I3,5X,A),2(1X,E14.8,1X,A),1X,F8.1)') &
                  'BASOPT| ', iset, "|", ibasis, "|", f_vec(icalc), "|", cond_vec(icalc), "|", tot_time(icalc)
            END DO
         END DO
         WRITE (unit_nr, '(1X,A)') "BASOPT| -----------------------------------------------------------------------"
         WRITE (unit_nr, '(1X,A,E14.8)') "BASOPT| Total residuum value: ", opt_bas%powell_param%f
         WRITE (unit_nr, '(A)') ""
      END IF
      CALL timestop(handle)
   END SUBROUTINE output_opt_info

END MODULE optimize_basis
