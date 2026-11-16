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
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: silent_print_level
   USE f77_interface,                   ONLY: create_force_env,&
                                              destroy_force_env,&
                                              f_env_add_defaults,&
                                              f_env_get_from_id,&
                                              f_env_rm_defaults,&
                                              f_env_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE global_types,                    ONLY: global_environment_type
   USE input_constants,                 ONLY: restart_guess
   USE input_cp2k_check,                ONLY: check_cp2k_input
   USE input_cp2k_read,                 ONLY: empty_initial_variables,&
                                              read_input
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
   USE optbas_frontier_orbital_screening,ONLY: run_frontier_orbital_screening
   USE optbas_frontier_orbitals_utils,  ONLY: evaluate_frontier_orbitals_objective,&
                                              frontier_orbitals_reference_init,&
                                              frontier_orbitals_reference_release,&
                                              frontier_orbitals_reference_type,&
                                              frontier_orbitals_result_type,&
                                              print_frontier_orbital_basis_sets
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
   USE qs_environment,                  ONLY: qs_init
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_create,&
                                              qs_env_part_release,&
                                              qs_env_release,&
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
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: run_optimize_basis

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_basis'

CONTAINS

! **************************************************************************************************
!> \brief main entry point for methods aimed at optimizing basis sets
!> \param input_declaration ...
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \author Florian Schiffmann
!>         Jan Wilhelm, 08/2026: Add run_frontier_orbital_screening
! **************************************************************************************************
   SUBROUTINE run_optimize_basis(input_declaration, root_section, para_env, globenv)
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(section_vals_type), POINTER                   :: root_section
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv

      CHARACTER(len=*), PARAMETER :: routineN = 'run_optimize_basis'

      INTEGER                                            :: handle
      LOGICAL                                            :: frontier_orbitals_explicit, &
                                                            frontier_screening_explicit
      TYPE(basis_optimization_type)                      :: opt_bas
      TYPE(section_vals_type), POINTER                   :: frontier_orbitals_section, &
                                                            frontier_screening_section, &
                                                            optbas_section

      CALL timeset(routineN, handle)

      optbas_section => section_vals_get_subs_vals(root_section, "OPTIMIZE_BASIS")
      frontier_orbitals_section => section_vals_get_subs_vals(optbas_section, "FRONTIER_ORBITALS")
      frontier_screening_section => &
         section_vals_get_subs_vals(optbas_section, "FRONTIER_ORBITAL_SCREENING")
      CALL section_vals_get(frontier_orbitals_section, explicit=frontier_orbitals_explicit)
      CALL section_vals_get(frontier_screening_section, explicit=frontier_screening_explicit)

      IF (frontier_orbitals_explicit .AND. frontier_screening_explicit) THEN
         CPABORT("FRONTIER_ORBITALS and FRONTIER_ORBITAL_SCREENING cannot be activated together")
      END IF

      IF (frontier_screening_explicit) THEN
         CALL run_frontier_orbital_screening( &
            input_declaration, root_section, para_env, globenv, frontier_screening_section, driver_para_opt_basis)
      ELSE
         CALL optimize_basis_init_read_input(opt_bas, root_section, para_env)
         CALL driver_para_opt_basis(opt_bas, input_declaration, para_env, globenv)
         CALL deallocate_basis_optimization_type(opt_bas)
      END IF

      CALL timestop(handle)

   END SUBROUTINE run_optimize_basis

! **************************************************************************************************
!> \brief driver routine for the parallel part of the method
!> \param opt_bas ...
!> \param input_declaration ...
!> \param para_env ...
!> \param globenv ...
!> \param training_input ...
!> \param reference ...
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE driver_para_opt_basis(opt_bas, input_declaration, para_env, globenv, training_input, reference)
      TYPE(basis_optimization_type)                      :: opt_bas
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), OPTIONAL, POINTER         :: training_input
      TYPE(frontier_orbitals_reference_type), &
         DIMENSION(:), INTENT(IN), OPTIONAL, TARGET      :: reference

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

      IF (PRESENT(training_input)) THEN
         CPASSERT(PRESENT(reference))
         CALL driver_optimization_para_low(opt_bas, input_declaration, para_env, opt_group, &
                                           globenv, training_input, reference)
      ELSE
         CALL driver_optimization_para_low(opt_bas, input_declaration, para_env, opt_group, globenv)
      END IF

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
!> \param globenv ...
!> \param training_input ...
!> \param reference ...
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE driver_optimization_para_low(opt_bas, input_declaration, para_env_top, mpi_comm_opt, &
                                           globenv, training_input, reference)
      TYPE(basis_optimization_type)                      :: opt_bas
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(mp_para_env_type), POINTER                    :: para_env_top
      TYPE(mp_comm_type), INTENT(IN)                     :: mpi_comm_opt
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), OPTIONAL, POINTER         :: training_input
      TYPE(frontier_orbitals_reference_type), &
         DIMENSION(:), INTENT(IN), OPTIONAL, TARGET      :: reference

      CHARACTER(len=*), PARAMETER :: routineN = 'driver_optimization_para_low'

      INTEGER                                            :: handle, icalc, iopt, is, mp_id, &
                                                            outer_print_level, stat
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: f_env_id
      LOGICAL                                            :: use_direct_qs, write_basis
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tot_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: matrix_S_inv
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(frontier_orbitals_reference_type), &
         ALLOCATABLE, DIMENSION(:), TARGET               :: local_reference
      TYPE(frontier_orbitals_reference_type), &
         DIMENSION(:), POINTER                           :: active_reference
      TYPE(frontier_orbitals_result_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: initial_info, optimized_info
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_environment_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: training_qs_env

      NULLIFY (active_reference, f_env, logger)

      CALL timeset(routineN, handle)

      ! ======  initialize the f_env and precompute some matrices =====
      mp_id = opt_bas%opt_id
      use_direct_qs = opt_bas%method == method_mo_fit_occ_virtual
      NULLIFY (para_env, f_env)
      ALLOCATE (f_env_id(SIZE(opt_bas%comp_group(mp_id)%member_list)), source=-1)
      ALLOCATE (tot_time(opt_bas%ncombinations*opt_bas%ntraining_sets))
      IF (opt_bas%method == method_mo_fit_occ) THEN
         ALLOCATE (matrix_s_inv(SIZE(opt_bas%comp_group(mp_id)%member_list)))
      ELSE
         ALLOCATE (initial_info(opt_bas%ncombinations*opt_bas%ntraining_sets))
         ALLOCATE (optimized_info(opt_bas%ncombinations*opt_bas%ntraining_sets))
         IF (PRESENT(reference)) THEN
            active_reference => reference
         ELSE
            ALLOCATE (local_reference(SIZE(opt_bas%comp_group(mp_id)%member_list)))
            active_reference => local_reference
         END IF
      END IF

      ALLOCATE (para_env)
      para_env = mpi_comm_opt

      is = -1
      IF (para_env%is_source()) is = para_env_top%mepos
      CALL para_env_top%allgather(is, opt_bas%sub_sources)

      IF (use_direct_qs) THEN
         ALLOCATE (training_qs_env(SIZE(opt_bas%comp_group(mp_id)%member_list)))
         logger => cp_get_default_logger()
         outer_print_level = logger%iter_info%print_level
         logger%iter_info%print_level = silent_print_level
         IF (PRESENT(training_input)) THEN
            CALL init_training_qs_envs(opt_bas, training_qs_env, input_declaration, para_env, globenv, &
                                       training_input=training_input)
         ELSE
            CALL init_training_qs_envs(opt_bas, training_qs_env, input_declaration, para_env, globenv, &
                                       reference=local_reference)
         END IF
      ELSE
         CALL init_training_force_envs(opt_bas, f_env_id, input_declaration, para_env, mpi_comm_opt, &
                                       matrix_s_inv=matrix_s_inv)
      END IF

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
            IF (iopt == 0) THEN
               CALL compute_frontier_orbitals_objective( &
                  opt_bas, f_env_id, active_reference, tot_time, &
                  para_env_top, para_env, initial_info, training_qs_env)
            ELSE
               CALL compute_frontier_orbitals_objective( &
                  opt_bas, f_env_id, active_reference, tot_time, &
                  para_env_top, para_env, qs_envs=training_qs_env)
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
         write_basis = opt_bas%method == method_mo_fit_occ .AND. &
                       MOD(iopt, opt_bas%write_frequency) == 0
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
      IF (opt_bas%method == method_mo_fit_occ_virtual) THEN
         ! Report the actual best Powell vector, rather than the final trial point.
         CALL compute_frontier_orbitals_objective( &
            opt_bas, f_env_id, active_reference, tot_time, &
            para_env_top, para_env, optimized_info, training_qs_env)
         IF (.NOT. opt_bas%quiet_output) THEN
            CALL print_frontier_orbitals_info( &
               opt_bas, active_reference, initial_info, optimized_info, &
               para_env_top, para_env)
         END IF
      END IF
      CALL update_derived_basis_sets(opt_bas,.NOT. opt_bas%quiet_output, opt_bas%output_basis_file, &
                                     para_env_top)

      ! ======  get rid of the f_env again =====

      DO icalc = SIZE(opt_bas%comp_group(mp_id)%member_list), 1, -1
         IF (opt_bas%method == method_mo_fit_occ_virtual .AND. ALLOCATED(local_reference)) THEN
            CALL frontier_orbitals_reference_release(local_reference(icalc))
         END IF
         IF (use_direct_qs) THEN
            CALL qs_env_release(training_qs_env(icalc))
         ELSE
            CALL f_env_get_from_id(f_env_id(icalc), f_env)
            CALL destroy_force_env(f_env_id(icalc), stat)
         END IF
      END DO
      IF (ALLOCATED(training_qs_env)) THEN
         DEALLOCATE (training_qs_env)
         logger%iter_info%print_level = outer_print_level
      END IF
      DEALLOCATE (f_env_id); DEALLOCATE (tot_time)
      IF (ALLOCATED(initial_info)) DEALLOCATE (initial_info)
      IF (ALLOCATED(optimized_info)) DEALLOCATE (optimized_info)
      IF (ALLOCATED(local_reference)) DEALLOCATE (local_reference)
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
!> \param info ...
!> \param qs_envs ...
! **************************************************************************************************
   SUBROUTINE compute_frontier_orbitals_objective(opt_bas, f_env_id, reference, tot_time, &
                                                  para_env_top, para_env, info, qs_envs)
      TYPE(basis_optimization_type)                      :: opt_bas
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: f_env_id
      TYPE(frontier_orbitals_reference_type), &
         DIMENSION(:), INTENT(IN)                        :: reference
      REAL(KIND=dp), DIMENSION(:)                        :: tot_time
      TYPE(mp_para_env_type), POINTER                    :: para_env_top, para_env
      TYPE(frontier_orbitals_result_type), &
         DIMENSION(:), INTENT(OUT), OPTIONAL             :: info
      TYPE(qs_environment_type), DIMENSION(:), &
         OPTIONAL, TARGET                                :: qs_envs

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_frontier_orbitals_objective'

      CHARACTER(len=8)                                   :: basis_type
      INTEGER                                            :: bas_id, handle, icalc, icomb, icont, &
                                                            ikind, ipgf, iset, ix, mp_id, my_id, &
                                                            ncalc, set_id
      REAL(KIND=dp)                                      :: coefficient_loss, contribution
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coefficients_current
      REAL(KIND=dp), DIMENSION(:), POINTER               :: start_time
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s_aux, matrix_s_aux_orb
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(frontier_orbitals_result_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: results
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux, sab_aux_orb
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      basis_type = "AUX_OPT"
      ncalc = opt_bas%ncombinations*opt_bas%ntraining_sets
      ALLOCATE (results(ncalc))
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
         IF (PRESENT(qs_envs)) THEN
            qs_env => qs_envs(icalc)
         ELSE
            CALL f_env_get_from_id(f_env_id(icalc), f_env)
            force_env => f_env%force_env
            CALL force_env_get(force_env, qs_env=qs_env)
         END IF
         CALL update_basis_set(opt_bas, bas_id, basis_type, qs_env)
         CALL optbas_build_neighborlist(qs_env, sab_aux, sab_aux_orb, basis_type)
         CALL get_qs_env(qs_env, ks_env=ks_env)
         CALL build_overlap_matrix(ks_env, matrix_s=matrix_s_aux, &
                                   basis_type_a=basis_type, basis_type_b=basis_type, sab_nl=sab_aux)
         CALL build_overlap_matrix(ks_env, matrix_s=matrix_s_aux_orb, &
                                   basis_type_a=basis_type, basis_type_b="ORB", sab_nl=sab_aux_orb)
         CALL release_neighbor_list_sets(sab_aux)
         CALL release_neighbor_list_sets(sab_aux_orb)

         CALL evaluate_frontier_orbitals_objective( &
            reference(icalc), matrix_s_aux(1)%matrix, matrix_s_aux_orb(1)%matrix, &
            opt_bas%virtual_energy_cutoff, opt_bas%virtual_energy_smoothing, &
            opt_bas%gap_energy_scale, results(my_id))
         results(my_id)%reference_condition_number = reference(icalc)%condition_number
         results(my_id)%number_reference_orbitals = reference(icalc)%number_reference_orbitals
         results(my_id)%calculation_time = m_walltime() - start_time(icalc)

         IF (ASSOCIATED(matrix_s_aux)) CALL dbcsr_deallocate_matrix_set(matrix_s_aux)
         IF (ASSOCIATED(matrix_s_aux_orb)) CALL dbcsr_deallocate_matrix_set(matrix_s_aux_orb)
      END DO
      DEALLOCATE (start_time)

      CALL reduce_frontier_orbitals_results(results, para_env_top, para_env)
      IF (PRESENT(info)) THEN
         CPASSERT(SIZE(info) == ncalc)
         info = results
      END IF

      opt_bas%powell_param%f = 0.0_dp
      IF (para_env_top%is_source()) THEN
         DO icalc = 1, ncalc
            icomb = MOD(icalc - 1, opt_bas%ncombinations)
            contribution = opt_bas%occupied_weight*results(icalc)%loss_occupied_subspace + &
                           opt_bas%virtual_weight*results(icalc)%loss_virtual_subspace + &
                           opt_bas%empty_overlap_weight*results(icalc)%loss_empty_subspace + &
                           opt_bas%gap_weight*results(icalc)%loss_gap + &
                           opt_bas%coefficient_weight*coefficient_loss
            IF (opt_bas%use_condition_number) THEN
               contribution = contribution + opt_bas%condition_weight(icomb)* &
                              LOG10(results(icalc)%condition_number)
            END IF
            opt_bas%powell_param%f = opt_bas%powell_param%f + &
                                     opt_bas%fval_weight(icomb)*contribution
            tot_time(icalc) = tot_time(icalc) + results(icalc)%calculation_time
         END DO
      END IF
      CALL para_env_top%bcast(opt_bas%powell_param%f)
      DEALLOCATE (coefficients_current, results)

      CALL para_env_top%sync()
      CALL timestop(handle)

   END SUBROUTINE compute_frontier_orbitals_objective

! **************************************************************************************************
!> \brief Collect frontier-orbital results from all calculation groups.
!> \param results frontier-orbital results for every training calculation
!> \param para_env_top top-level communicator
!> \param para_env calculation-group communicator
! **************************************************************************************************
   SUBROUTINE reduce_frontier_orbitals_results(results, para_env_top, para_env)
      TYPE(frontier_orbitals_result_type), &
         DIMENSION(:), INTENT(INOUT)                     :: results
      TYPE(mp_para_env_type), POINTER                    :: para_env_top, para_env

      IF (.NOT. para_env%is_source()) results = frontier_orbitals_result_type()
      CALL para_env_top%sum(results%loss_occupied_subspace)
      CALL para_env_top%sum(results%loss_virtual_subspace)
      CALL para_env_top%sum(results%loss_empty_subspace)
      CALL para_env_top%sum(results%loss_gap)
      CALL para_env_top%sum(results%occupied_subspace_overlap)
      CALL para_env_top%sum(results%virtual_subspace_similarity)
      CALL para_env_top%sum(results%empty_subspace_overlap)
      CALL para_env_top%sum(results%gap_reference)
      CALL para_env_top%sum(results%gap_candidate)
      CALL para_env_top%sum(results%condition_number)
      CALL para_env_top%sum(results%reference_condition_number)
      CALL para_env_top%sum(results%calculation_time)
      CALL para_env_top%sum(results%number_candidate_orbitals)
      CALL para_env_top%sum(results%number_reference_orbitals)

   END SUBROUTINE reduce_frontier_orbitals_results

! **************************************************************************************************
!> \brief Print initial-to-optimized frontier-orbital information for every training calculation.
!> \param opt_bas basis optimization settings
!> \param reference locally held reference data and basis metadata
!> \param initial initial-small-basis results
!> \param optimized optimized-small-basis results
!> \param para_env_top top-level communicator
!> \param para_env calculation-group communicator
! **************************************************************************************************
   SUBROUTINE print_frontier_orbitals_info(opt_bas, reference, initial, optimized, &
                                           para_env_top, para_env)
      TYPE(basis_optimization_type), INTENT(IN)          :: opt_bas
      TYPE(frontier_orbitals_reference_type), &
         DIMENSION(:), INTENT(IN)                        :: reference
      TYPE(frontier_orbitals_result_type), &
         DIMENSION(:), INTENT(IN)                        :: initial, optimized
      TYPE(mp_para_env_type), POINTER                    :: para_env_top, para_env

      CHARACTER(LEN=17)                                  :: basis_label
      CHARACTER(LEN=default_string_length), &
         ALLOCATABLE, DIMENSION(:)                       :: candidate_basis, element, reference_basis
      INTEGER                                            :: bas_id, icalc, ikind, iset, max_kinds, &
                                                            my_id, ncalc, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_count
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: metadata
      TYPE(cp_logger_type), POINTER                      :: logger

      ncalc = opt_bas%ncombinations*opt_bas%ntraining_sets
      CPASSERT(SIZE(initial) == ncalc)
      CPASSERT(SIZE(optimized) == ncalc)

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
            ALLOCATE (candidate_basis(kind_count(my_id)), element(kind_count(my_id)), &
                      reference_basis(kind_count(my_id)))
            DO ikind = 1, kind_count(my_id)
               CALL decode_string(metadata(my_id, ikind, 1, :), element(ikind))
               CALL decode_string(metadata(my_id, ikind, 2, :), reference_basis(ikind))
               CALL decode_string(metadata(my_id, ikind, 3, :), candidate_basis(ikind))
            END DO
            CALL print_frontier_orbital_basis_sets(unit_nr, element, reference_basis, candidate_basis)
            DEALLOCATE (candidate_basis, element, reference_basis)
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
               basis_label, initial(my_id)%number_reference_orbitals, &
               100.0_dp, "%", 100.0_dp, "%", 100.0_dp, "%", &
               initial(my_id)%reference_condition_number, &
               initial(my_id)%gap_reference*evolt
            basis_label = "Init. small basis"
            WRITE (unit_nr, '(A17,1X,I4,2(1X,F10.1,A),1X,F11.1,A,1X,ES9.2,1X,F6.3)') &
               basis_label, initial(my_id)%number_candidate_orbitals, &
               100.0_dp*initial(my_id)%occupied_subspace_overlap, "%", &
               100.0_dp*initial(my_id)%virtual_subspace_similarity, "%", &
               100.0_dp*initial(my_id)%empty_subspace_overlap, "%", &
               initial(my_id)%condition_number, initial(my_id)%gap_candidate*evolt
            basis_label = "Opt. small basis"
            WRITE (unit_nr, '(A17,1X,I4,2(1X,F10.1,A),1X,F11.1,A,1X,ES9.2,1X,F6.3)') &
               basis_label, optimized(my_id)%number_candidate_orbitals, &
               100.0_dp*optimized(my_id)%occupied_subspace_overlap, "%", &
               100.0_dp*optimized(my_id)%virtual_subspace_similarity, "%", &
               100.0_dp*optimized(my_id)%empty_subspace_overlap, "%", &
               optimized(my_id)%condition_number, optimized(my_id)%gap_candidate*evolt
            WRITE (unit_nr, '(A)') REPEAT("-", 78)

         END DO

         WRITE (unit_nr, '(A)') REPEAT("=", 78)
         WRITE (unit_nr, '(A)') ""
      END IF

      DEALLOCATE (kind_count, metadata)

   END SUBROUTINE print_frontier_orbitals_info

! **************************************************************************************************
!> \brief Encode a string as integers so basis metadata can be reduced across optimizer groups.
!> \param string ...
!> \param encoded ...
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
!> \param encoded ...
!> \param string ...
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
!> \brief Initialize Quickstep environments without starting nested CP2K runs.
!> \param opt_bas ...
!> \param qs_envs ...
!> \param input_declaration ...
!> \param para_env ...
!> \param globenv ...
!> \param training_input ...
!> \param reference ...
! **************************************************************************************************
   SUBROUTINE init_training_qs_envs(opt_bas, qs_envs, input_declaration, para_env, globenv, &
                                    training_input, reference)
      TYPE(basis_optimization_type)                      :: opt_bas
      TYPE(qs_environment_type), DIMENSION(:), &
         INTENT(OUT), TARGET                             :: qs_envs
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), OPTIONAL, POINTER         :: training_input
      TYPE(frontier_orbitals_reference_type), &
         DIMENSION(:), INTENT(OUT), OPTIONAL             :: reference

      CHARACTER(len=*), PARAMETER :: routineN = 'init_training_qs_envs'

      CHARACTER(len=default_path_length)                 :: main_dir, restart_file
      CHARACTER(len=default_string_length)               :: project_name
      INTEGER                                            :: bas_id, handle, icalc, ierr, mp_id, &
                                                            n_rep_val, set_id
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_environment_type), POINTER                 :: qs_env, reference_qs_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: dft_section, force_env_section, &
                                                            input_file, subsys_section

      CALL timeset(routineN, handle)

      mp_id = opt_bas%opt_id
      CALL m_getcwd(main_dir)
      DO icalc = 1, SIZE(opt_bas%comp_group(mp_id)%member_list)
         NULLIFY (dft_section, force_env_section, input_file, ks_env, matrix_s, qs_env, &
                  reference_qs_env, sab_orb, subsys_section)
         CALL get_set_and_basis_id(opt_bas%comp_group(mp_id)%member_list(icalc), opt_bas, &
                                   set_id=set_id, bas_id=bas_id)
         IF (PRESENT(training_input)) THEN
            CALL section_vals_duplicate(training_input, input_file)
         ELSE
            CALL m_chdir(TRIM(opt_bas%training_dir(set_id)), ierr)
            IF (ierr /= 0) THEN
               CALL cp_abort(__LOCATION__, &
                             "Could not change to directory <"//TRIM(opt_bas%training_dir(set_id))//">")
            END IF
            input_file => read_input(input_declaration, &
                                     opt_bas%training_input(set_id), &
                                     initial_variables=empty_initial_variables, &
                                     para_env=para_env)
            CALL check_cp2k_input(input_declaration, input_file, para_env=para_env, &
                                  output_unit=-1)
         END IF

         CALL section_vals_val_set(input_file, "GLOBAL%PRINT_LEVEL", i_val=silent_print_level)
         IF (PRESENT(reference)) THEN
            force_env_section => section_vals_get_subs_vals(input_file, "FORCE_EVAL")
            dft_section => section_vals_get_subs_vals(force_env_section, "DFT")
            CALL section_vals_val_get(dft_section, "WFN_RESTART_FILE_NAME", n_rep_val=n_rep_val)
            IF (n_rep_val == 0) THEN
               CALL section_vals_val_get(input_file, "GLOBAL%PROJECT_NAME", c_val=project_name)
               restart_file = TRIM(project_name)//"-RESTART.wfn"
               CALL section_vals_val_set(dft_section, "WFN_RESTART_FILE_NAME", &
                                         c_val=TRIM(restart_file))
            END IF
            CALL section_vals_val_set(dft_section, "SCF%SCF_GUESS", i_val=restart_guess)
            subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")
            ALLOCATE (reference_qs_env)
            CALL qs_env_create(reference_qs_env, globenv)
            CALL qs_init(reference_qs_env, para_env, input_file, globenv=globenv, &
                         force_env_section=force_env_section, subsys_section=subsys_section, &
                         use_motion_section=.TRUE., silent=.TRUE.)
            CALL allocate_mo_sets(reference_qs_env)
            CALL get_qs_env(reference_qs_env, ks_env=ks_env)
            CALL build_qs_neighbor_lists(reference_qs_env, para_env, molecular=.FALSE., &
                                         force_env_section=reference_qs_env%input)
            CALL get_ks_env(ks_env, matrix_s=matrix_s, sab_orb=sab_orb)
            CALL build_overlap_matrix(ks_env, matrix_s=matrix_s, &
                                      matrix_name="OVERLAP", &
                                      basis_type_a="ORB", &
                                      basis_type_b="ORB", &
                                      sab_nl=sab_orb)
            CALL set_ks_env(ks_env, matrix_s=matrix_s)
            CALL calculate_ks_matrix(reference_qs_env)
            CALL frontier_orbitals_reference_init(reference_qs_env, reference(icalc))
            CALL qs_env_release(reference_qs_env)
            DEALLOCATE (reference_qs_env)
         END IF

         CALL modify_input_settings(opt_bas, bas_id, input_file)
         force_env_section => section_vals_get_subs_vals(input_file, "FORCE_EVAL")
         subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")
         qs_env => qs_envs(icalc)
         CALL qs_env_create(qs_env, globenv)
         CALL qs_init(qs_env, para_env, input_file, globenv=globenv, &
                      force_env_section=force_env_section, subsys_section=subsys_section, &
                      use_motion_section=.TRUE., silent=.TRUE.)

         CALL section_vals_release(input_file)
         IF (.NOT. PRESENT(training_input)) CALL m_chdir(TRIM(ADJUSTL(main_dir)), ierr)
      END DO

      CALL timestop(handle)

   END SUBROUTINE init_training_qs_envs

! **************************************************************************************************
!> \brief create the force_envs for every input in the computational group
!> \param opt_bas ...
!> \param f_env_id ...
!> \param input_declaration ...
!> \param para_env ...
!> \param mpi_comm_opt ...
!> \param matrix_s_inv ...
!> \param training_input ...
!> \param reference ...
! **************************************************************************************************

   SUBROUTINE init_training_force_envs(opt_bas, f_env_id, input_declaration, para_env, mpi_comm_opt, &
                                       matrix_s_inv, training_input, reference)

      TYPE(basis_optimization_type)                      :: opt_bas
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: f_env_id
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(mp_comm_type)                                 :: mpi_comm_opt
      TYPE(cp_fm_type), DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: matrix_S_inv
      TYPE(section_vals_type), OPTIONAL, POINTER         :: training_input
      TYPE(frontier_orbitals_reference_type), &
         DIMENSION(:), INTENT(OUT), OPTIONAL             :: reference

      CHARACTER(len=*), PARAMETER :: routineN = 'init_training_force_envs'

      CHARACTER(len=default_path_length)                 :: input_path, main_dir
      INTEGER                                            :: bas_id, f_env_output_unit, handle, &
                                                            icalc, ierr, mp_id, set_id, stat
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
         ! Parse an external training input or duplicate the embedded screening input.
         CALL get_set_and_basis_id(opt_bas%comp_group(mp_id)%member_list(icalc), opt_bas, set_id, bas_id)
         IF (PRESENT(training_input)) THEN
            CALL section_vals_duplicate(training_input, input_file)
            input_path = "embedded frontier-orbital screening input"
         ELSE
            CALL m_chdir(TRIM(opt_bas%training_dir(set_id)), ierr)
            IF (ierr /= 0) THEN
               CALL cp_abort(__LOCATION__, &
                             "Could not change to directory <"//TRIM(opt_bas%training_dir(set_id))//">")
            END IF
            input_file => read_input(input_declaration, &
                                     opt_bas%training_input(set_id), &
                                     initial_variables=empty_initial_variables, &
                                     para_env=para_env)
            input_path = opt_bas%training_input(set_id)
         END IF

         CALL modify_input_settings(opt_bas, bas_id, input_file)
         IF (opt_bas%quiet_output .OR. opt_bas%method == method_mo_fit_occ_virtual) THEN
            IF (opt_bas%quiet_output) THEN
               CALL section_vals_val_set(input_file, "GLOBAL%PRINT_LEVEL", i_val=silent_print_level)
            END IF
            f_env_output_unit = -1
            IF (para_env%is_source()) f_env_output_unit = cp_logger_get_default_unit_nr()
            CALL create_force_env(f_env_id(icalc), &
                                  input_declaration=input_declaration, &
                                  input_path=input_path, &
                                  input=input_file, &
                                  output_unit=f_env_output_unit, &
                                  owns_out_unit=.FALSE., &
                                  mpi_comm=mpi_comm_opt, &
                                  ierr=stat)
         ELSE
            CALL create_force_env(f_env_id(icalc), &
                                  input_declaration=input_declaration, &
                                  input_path=input_path, &
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
         IF (.NOT. PRESENT(training_input)) CALL allocate_mo_sets(qs_env)
         CALL f_env_rm_defaults(f_env, stat)
         IF (.NOT. PRESENT(training_input)) THEN
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
            IF (PRESENT(matrix_s_inv)) THEN
               CALL get_qs_env(qs_env, matrix_s=matrix_s, blacs_env=blacs_env)
               CALL calculate_overlap_inverse(matrix_s(1)%matrix, matrix_s_inv(icalc), &
                                              para_env, blacs_env)
            END IF
            CALL calculate_ks_matrix(qs_env)
         END IF
         IF (opt_bas%method == method_mo_fit_occ_virtual .AND. PRESENT(reference)) THEN
            CALL frontier_orbitals_reference_init(qs_env, reference(icalc))
         END IF

         CALL section_vals_release(input_file)

         CALL qs_env_part_release(qs_env)

         IF (.NOT. PRESENT(training_input)) CALL m_chdir(TRIM(ADJUSTL(main_dir)), ierr)
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
