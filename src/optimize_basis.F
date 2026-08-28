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
   USE f77_interface,                   ONLY: create_force_env,&
                                              destroy_force_env,&
                                              f_env_add_defaults,&
                                              f_env_get_from_id,&
                                              f_env_rm_defaults,&
                                              f_env_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_cp2k_read,                 ONLY: empty_initial_variables,&
                                              read_input
   USE input_section_types,             ONLY: section_type,&
                                              section_vals_release,&
                                              section_vals_type
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
                                              update_derived_basis_sets
   USE physcon,                         ONLY: evolt
   USE powell,                          ONLY: powell_optimize
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
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: run_optimize_basis

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_basis'

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
      TYPE(basis_optimization_type)                      :: opt_bas

      CALL timeset(routineN, handle)

      CALL optimize_basis_init_read_input(opt_bas, root_section, para_env)

      CALL driver_para_opt_basis(opt_bas, input_declaration, para_env)

      CALL deallocate_basis_optimization_type(opt_bas)
      CALL timestop(handle)

   END SUBROUTINE run_optimize_basis

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
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: f_env_id
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

      CALL init_training_force_envs(opt_bas, f_env_id, input_declaration, matrix_s_inv, &
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
            IF (iopt == 0) THEN
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
      IF (opt_bas%method == method_mo_fit_occ_virtual) THEN
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
      END DO
      DEALLOCATE (f_env_id); DEALLOCATE (tot_time)
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
                                                            my_id, ncalc, set_id, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: candidate_angular_momentum
      REAL(KIND=dp)                                      :: coefficient_loss, contribution
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coefficients_current
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gdata
      REAL(KIND=dp), DIMENSION(:), POINTER               :: start_time
      TYPE(cp_logger_type), POINTER                      :: logger
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

         logger => cp_get_default_logger()
         unit_nr = cp_logger_get_default_unit_nr(logger)
         IF (unit_nr > 0) THEN
            DO icalc = 1, ncalc
               WRITE (unit_nr, &
                      '(1X,A,I0,1X,A,F0.3,1X,A,F0.3,1X,A,F0.1,A,F0.1,A,F0.1,A,1X,A,I0,1X,A,I0)') &
                  "i=", iopt, &
                  "gap=", gdata(icalc, 6)*evolt, &
                  "ref=", gdata(icalc, 5)*evolt, &
                  "occ/vir/emp=", 100.0_dp*gdata(icalc, frontier_idx_overlap_occupied), "/", &
                  100.0_dp*gdata(icalc, frontier_idx_similarity_virtual), "/", &
                  100.0_dp*gdata(icalc, frontier_idx_overlap_empty), "%", &
                  "κ(S)=", NINT(gdata(icalc, frontier_idx_condition_candidate)), &
                  "L=", NINT(1.0E4_dp*opt_bas%powell_param%f)
            END DO
         END IF
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

   SUBROUTINE init_training_force_envs(opt_bas, f_env_id, input_declaration, matrix_s_inv, &
                                       frontier_orbitals_reference, para_env, mpi_comm_opt)

      TYPE(basis_optimization_type)                      :: opt_bas
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: f_env_id
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
         CALL create_force_env(f_env_id(icalc), &
                               input_declaration=input_declaration, &
                               input_path=opt_bas%training_input(set_id), &
                               input=input_file, &
                               output_path="scrap_information", &
                               mpi_comm=mpi_comm_opt, &
                               ierr=stat)

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
