!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for a Harris type energy correction on top of a
!>        Kohn-Sham calculation
!> \par History
!>       03.2014 created
!>       09.2019 Moved from KG to Kohn-Sham
!> \author JGH
! **************************************************************************************************
MODULE energy_corrections
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE core_ppl,                        ONLY: build_core_ppl
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_triangular_invert
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_init_random,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_triangular,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_result_methods,               ONLY: cp_results_erase,&
                                              put_results
   USE cp_result_types,                 ONLY: cp_result_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, dbcsr_distribution_type, &
        dbcsr_dot, dbcsr_filter, dbcsr_get_info, dbcsr_init_p, dbcsr_multiply, dbcsr_p_type, &
        dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE dm_ls_scf,                       ONLY: calculate_w_matrix_ls,&
                                              post_scf_sparsities
   USE dm_ls_scf_methods,               ONLY: apply_matrix_preconditioner,&
                                              density_matrix_sign,&
                                              density_matrix_tc2,&
                                              density_matrix_trs4,&
                                              ls_scf_init_matrix_s
   USE dm_ls_scf_qs,                    ONLY: matrix_ls_create,&
                                              matrix_ls_to_qs,&
                                              matrix_qs_to_ls
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE ec_efield_local,                 ONLY: ec_efield_integrals,&
                                              ec_efield_local_operator
   USE ec_env_types,                    ONLY: energy_correction_type
   USE ec_methods,                      ONLY: create_kernel,&
                                              ec_mos_init
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE input_constants,                 ONLY: &
        ec_diagonalization, ec_functional_harris, ec_matrix_sign, ec_matrix_tc2, ec_matrix_trs4, &
        ec_ot_atomic, ec_ot_diag, ec_ot_gs, ot_precond_full_single_inverse, &
        ot_precond_solver_default, vdw_pairpot_dftd3, vdw_pairpot_dftd3bj, xc_vdw_fun_pairpot
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE mao_basis,                       ONLY: mao_generate_basis
   USE mathlib,                         ONLY: det_3x3
   USE message_passing,                 ONLY: mp_sum
   USE molecule_types,                  ONLY: molecule_type
   USE moments_utils,                   ONLY: get_reference_point
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: ptable
   USE physcon,                         ONLY: bohr,&
                                              debye,&
                                              pascal
   USE preconditioner,                  ONLY: make_preconditioner
   USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                              init_preconditioner,&
                                              preconditioner_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integral_ab,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ptrace
   USE qs_density_matrices,             ONLY: calculate_density_matrix,&
                                              calculate_w_matrix
   USE qs_dispersion_pairpot,           ONLY: calculate_dispersion_pairpot
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type,&
                                              total_qs_force,&
                                              zero_qs_force
   USE qs_initial_guess,                ONLY: calculate_atomic_block_dm
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues,&
                                              make_basis_sm
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              get_mo_set,&
                                              mo_set_p_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE qs_ot_eigensolver,               ONLY: ot_eigensolver
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE response_solver,                 ONLY: ks_ref_potential,&
                                              response_calculation,&
                                              response_force
   USE string_utilities,                ONLY: uppercase
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE task_list_types,                 ONLY: allocate_task_list,&
                                              deallocate_task_list
   USE virial_methods,                  ONLY: one_third_sum_diag,&
                                              write_stress_tensor,&
                                              write_stress_tensor_components
   USE virial_types,                    ONLY: cp_virial,&
                                              symmetrize_virial,&
                                              virial_create,&
                                              virial_release,&
                                              virial_type
   USE voronoi_interface,               ONLY: entry_voronoi_or_bqb
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'energy_corrections'

   LOGICAL, PARAMETER                   :: debug_forces = .TRUE.
   LOGICAL, PARAMETER                   :: debug_stress = .TRUE.

   PUBLIC :: energy_correction

CONTAINS

! **************************************************************************************************
!> \brief Energy correction to a KG simulation
!>
!> \param qs_env ...
!> \param ec_init ...
!> \param calculate_forces ...
!> \par History
!>       03.2014 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE energy_correction(qs_env, ec_init, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: ec_init, calculate_forces

      CHARACTER(len=*), PARAMETER                        :: routineN = 'energy_correction'

      INTEGER                                            :: handle, unit_nr
      LOGICAL                                            :: my_calc_forces
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: ks_force

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      NULLIFY (ec_env)
      CALL get_qs_env(qs_env, ec_env=ec_env)

      ! Skip Harris functional calculation if ground-state is NOT converged
      IF (.NOT. ec_env%do_skip) THEN

         ec_env%should_update = .TRUE.
         IF (PRESENT(ec_init)) ec_env%should_update = ec_init

         my_calc_forces = .FALSE.
         IF (PRESENT(calculate_forces)) my_calc_forces = calculate_forces

         IF (ec_env%should_update) THEN
            ec_env%etotal = 0.0_dp
            ec_env%eband = 0.0_dp
            ec_env%ehartree = 0.0_dp
            ec_env%exc = 0.0_dp
            ec_env%vhxc = 0.0_dp
            ec_env%edispersion = 0.0_dp
         END IF

         IF (my_calc_forces) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,A,A,A,A)') "!", REPEAT("-", 25), &
                  " Energy Correction Forces ", REPEAT("-", 26), "!"
            END IF
            CALL get_qs_env(qs_env, force=ks_force)
            CALL zero_qs_force(ks_force)
         ELSE
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,A,A,A,A)') "!", REPEAT("-", 29), &
                  " Energy Correction ", REPEAT("-", 29), "!"
            END IF

         END IF

         CALL get_qs_env(qs_env, energy=energy)
         SELECT CASE (ec_env%energy_functional)
         CASE (ec_functional_harris)
            !
            CALL harris_energy(qs_env, ec_env, my_calc_forces, unit_nr)
            !
            IF (ec_env%should_update) THEN
               energy%nonscf_correction = ec_env%etotal - energy%total
               CALL evaluate_ec_core_matrix_traces(qs_env, ec_env)
               energy%total = ec_env%etotal
            END IF
         CASE DEFAULT
            CPABORT("unknown energy correction")
         END SELECT

         IF (.NOT. my_calc_forces .AND. unit_nr > 0) THEN
            WRITE (unit_nr, '(T3,A,T56,F25.15)') "Energy Correction ", energy%nonscf_correction
         END IF
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A,A,A)') "!", REPEAT("-", 77), "!"
         END IF

      ELSE

         ! Ground-state energy calculation did not converge,
         ! do not calculate energy correction
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A,A,A)') "!", REPEAT("-", 77), "!"
            WRITE (unit_nr, '(T2,A,A,A,A,A)') "!", REPEAT("-", 26), &
               " Skip Energy Correction ", REPEAT("-", 27), "!"
            WRITE (unit_nr, '(T2,A,A,A)') "!", REPEAT("-", 77), "!"
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE energy_correction

! **************************************************************************************************
!> \brief Calculates the traces of the core matrices and the density matrix.
!> \param qs_env ...
!> \param ec_env ...
!> \author Ole Schuett
!>         adapted for energy correction fbelle
! **************************************************************************************************
   SUBROUTINE evaluate_ec_core_matrix_traces(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_ec_core_matrix_traces'

      INTEGER                                            :: handle
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy

      CALL timeset(routineN, handle)
      NULLIFY (energy)

      CALL get_qs_env(qs_env, dft_control=dft_control, energy=energy)

      ! Core hamiltonian energy
      CALL calculate_ptrace(ec_env%matrix_h, ec_env%matrix_p, energy%core, dft_control%nspins)

      ! kinetic energy
      CALL calculate_ptrace(ec_env%matrix_t, ec_env%matrix_p, energy%kinetic, dft_control%nspins)

      CALL timestop(handle)

   END SUBROUTINE evaluate_ec_core_matrix_traces

! **************************************************************************************************
!> \brief Harris Energy Correction to a Kohn-Sham simulation
!>
!> \param qs_env ...
!> \param ec_env ...
!> \param calculate_forces ...
!> \param unit_nr ...
!> \par History
!>       03.2014 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE harris_energy(qs_env, ec_env, calculate_forces, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN)                                :: unit_nr

      INTEGER                                            :: ispin, nspins
      REAL(KIND=dp)                                      :: exc
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      IF (ec_env%should_update) THEN
         CALL ec_build_neighborlist(qs_env, ec_env)
         !
         CALL ec_build_core_hamiltonian(qs_env, ec_env)
         CALL ks_ref_potential(qs_env, ec_env%vh_rspace, ec_env%vxc_rspace, ec_env%vtau_rspace, ec_env%vadmm_rspace, &
                               ec_env%ehartree, exc)
         CALL ec_build_ks_matrix(qs_env, ec_env)
         !
         IF (ec_env%mao) THEN
            ! MAO basis
            IF (ASSOCIATED(ec_env%mao_coef)) CALL dbcsr_deallocate_matrix_set(ec_env%mao_coef)
            NULLIFY (ec_env%mao_coef)
            CALL mao_generate_basis(qs_env, ec_env%mao_coef, ref_basis_set="HARRIS", molecular=.TRUE., &
                                    max_iter=ec_env%mao_max_iter, eps_grad=ec_env%mao_eps_grad, &
                                    eps1_mao=ec_env%mao_eps1, iolevel=ec_env%mao_iolevel, unit_nr=unit_nr)
         END IF
         !
         CALL ec_ks_solver(qs_env, ec_env)
         !
         ! dispersion through pairpotentials
         CALL ec_disp(qs_env, ec_env, calculate_forces=.FALSE.)
         CALL ec_energy(ec_env, unit_nr)

      END IF
      IF (calculate_forces) THEN
         CALL ec_disp(qs_env, ec_env, calculate_forces=.TRUE.)
         CALL ec_build_core_hamiltonian_force(qs_env, ec_env)
         CALL ec_build_ks_matrix_force(qs_env, ec_env)

         CALL response_calculation(qs_env, ec_env)

         ! Allocate response density on real space grid for use in properties
         ! Calculated in response_force
         CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, pw_env=pw_env)
         nspins = dft_control%nspins
         CPASSERT(ASSOCIATED(pw_env))
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         ALLOCATE (ec_env%rhoz_r(nspins))
         DO ispin = 1, nspins
            NULLIFY (ec_env%rhoz_r(nspins)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, ec_env%rhoz_r(ispin)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
         END DO

         CALL response_force(qs_env, &
                             vh_rspace=ec_env%vh_rspace, &
                             vxc_rspace=ec_env%vxc_rspace, &
                             vtau_rspace=ec_env%vtau_rspace, &
                             vadmm_rspace=ec_env%vadmm_rspace, &
                             matrix_hz=ec_env%matrix_hz, &
                             matrix_pz=ec_env%matrix_z, &
                             matrix_pz_admm=ec_env%z_admm, &
                             matrix_wz=ec_env%matrix_wz, &
                             rhopz_r=ec_env%rhoz_r, &
                             zehartree=ec_env%ehartree, &
                             p_env=ec_env%p_env, &
                             zexc=ec_env%exc)

         CALL ec_properties(qs_env, ec_env)

         ! Deallocate Harris density and response density on grid
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, ec_env%rhoout_r(ispin)%pw)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, ec_env%rhoz_r(ispin)%pw)
         END DO
         DEALLOCATE (ec_env%rhoout_r, ec_env%rhoz_r)
         IF (ASSOCIATED(ec_env%tauout_r)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, ec_env%tauout_r(ispin)%pw)
         END DO
         DEALLOCATE (ec_env%tauout_r)
         END IF

         ! Deallocate matrices
         IF (ASSOCIATED(ec_env%matrix_ks)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_ks)
         IF (ASSOCIATED(ec_env%matrix_h)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_h)
         IF (ASSOCIATED(ec_env%matrix_s)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_s)
         IF (ASSOCIATED(ec_env%matrix_t)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_t)
         IF (ASSOCIATED(ec_env%matrix_p)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_p)
         IF (ASSOCIATED(ec_env%matrix_w)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_w)
         IF (ASSOCIATED(ec_env%matrix_hz)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_hz)
         IF (ASSOCIATED(ec_env%matrix_wz)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_wz)
         IF (ASSOCIATED(ec_env%matrix_z)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_z)

      END IF

   END SUBROUTINE harris_energy

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ec_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE ec_disp(qs_env, ec_env, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      REAL(KIND=dp)                                      :: edisp

      CALL calculate_dispersion_pairpot(qs_env, ec_env%dispersion_env, edisp, calculate_forces)
      IF (.NOT. calculate_forces) ec_env%edispersion = ec_env%edispersion + edisp

   END SUBROUTINE ec_disp

! **************************************************************************************************
!> \brief Construction of the Core Hamiltonian Matrix
!>        Short version of qs_core_hamiltonian
!> \param qs_env ...
!> \param ec_env ...
!> \author Creation (03.2014,JGH)
! **************************************************************************************************
   SUBROUTINE ec_build_core_hamiltonian(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_build_core_hamiltonian'

      INTEGER                                            :: handle, iounit, nder, nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: calculate_forces, use_virial
      REAL(KIND=dp)                                      :: eps_ppnl
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      END IF

      ! no k-points possible
      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      ks_env=ks_env, &
                      para_env=para_env)
      nimages = dft_control%nimages
      IF (nimages /= 1) THEN
         CPABORT("K-points for Harris functional not implemented")
      END IF

      ! check for GAPW/GAPW_XC
      IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) THEN
         CPABORT("Harris functional for GAPW not implemented")
      END IF

      ! get neighbor lists, we need the full sab_orb list from the ec_env
      NULLIFY (sab_orb, sac_ppl, sap_ppnl)
      sab_orb => ec_env%sab_orb
      sac_ppl => ec_env%sac_ppl
      sap_ppnl => ec_env%sap_ppnl

      nder = 0
      ! Overlap and kinetic energy matrices
      CALL build_overlap_matrix(ks_env, matrixkp_s=ec_env%matrix_s, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="HARRIS", &
                                basis_type_b="HARRIS", &
                                sab_nl=sab_orb)
      CALL build_kinetic_matrix(ks_env, matrixkp_t=ec_env%matrix_t, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="HARRIS", &
                                sab_nl=sab_orb)

      ! initialize H matrix
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_h, 1, 1)
      ALLOCATE (ec_env%matrix_h(1, 1)%matrix)
      CALL dbcsr_create(ec_env%matrix_h(1, 1)%matrix, template=ec_env%matrix_s(1, 1)%matrix)
      CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_h(1, 1)%matrix, sab_orb)

      ! add kinetic energy
      CALL dbcsr_copy(ec_env%matrix_h(1, 1)%matrix, ec_env%matrix_t(1, 1)%matrix, &
                      keep_sparsity=.TRUE., name="CORE HAMILTONIAN MATRIX")

      ! compute the ppl contribution to the core hamiltonian
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)
      NULLIFY (cell_to_index, virial)
      use_virial = .FALSE.
      calculate_forces = .FALSE.
      IF (ASSOCIATED(sac_ppl)) THEN
         CALL build_core_ppl(ec_env%matrix_h, ec_env%matrix_p, force, &
                             virial, calculate_forces, use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ppl, &
                             nimages, cell_to_index, "HARRIS")
      END IF

      ! compute the ppnl contribution to the core hamiltonian ***
      eps_ppnl = dft_control%qs_control%eps_ppnl
      IF (ASSOCIATED(sap_ppnl)) THEN
         CALL build_core_ppnl(ec_env%matrix_h, ec_env%matrix_p, force, &
                              virial, calculate_forces, use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, sab_orb, sap_ppnl, &
                              eps_ppnl, nimages, cell_to_index, "HARRIS")
      END IF

      ! External field (nonperiodic case)
      ec_env%efield_nuclear = 0.0_dp
      CALL ec_efield_local_operator(qs_env, ec_env, calculate_forces)

      CALL timestop(handle)

   END SUBROUTINE ec_build_core_hamiltonian

! **************************************************************************************************
!> \brief Solve KS equation for a given matrix
!>        calculate the complete KS matrix
!> \param qs_env ...
!> \param ec_env ...
!> \par History
!>      03.2014 adapted from qs_ks_build_kohn_sham_matrix [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ec_build_ks_matrix(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_build_ks_matrix'

      CHARACTER(LEN=default_string_length)               :: headline
      INTEGER                                            :: handle, iounit, ispin, nspins
      LOGICAL                                            :: calculate_forces, use_virial
      REAL(dp)                                           :: eexc, evhxc
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, tau_r, v_rspace, v_tau_rspace
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      END IF

      calculate_forces = .FALSE.

      ! get all information on the electronic density
      NULLIFY (rho, ks_env)
      CALL get_qs_env(qs_env=qs_env, rho=rho, virial=virial, dft_control=dft_control, &
                      para_env=para_env, blacs_env=blacs_env, ks_env=ks_env)

      nspins = dft_control%nspins
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      use_virial = .FALSE.

      ! Kohn-Sham matrix
      IF (ASSOCIATED(ec_env%matrix_ks)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_ks)
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_ks, nspins, 1)
      DO ispin = 1, nspins
         headline = "KOHN-SHAM MATRIX"
         ALLOCATE (ec_env%matrix_ks(ispin, 1)%matrix)
         CALL dbcsr_create(ec_env%matrix_ks(ispin, 1)%matrix, name=TRIM(headline), &
                           template=ec_env%matrix_s(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_ks(ispin, 1)%matrix, ec_env%sab_orb)
         CALL dbcsr_set(ec_env%matrix_ks(ispin, 1)%matrix, 0.0_dp)
      END DO

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CPASSERT(ASSOCIATED(pw_env))

      ! v_rspace and v_tau_rspace are generated from the auxbas pool
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      NULLIFY (v_rspace, v_tau_rspace)
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=ec_env%xc_section, &
                         vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=eexc, just_energy=.FALSE.)
      IF (.NOT. ASSOCIATED(v_rspace)) THEN
         ALLOCATE (v_rspace(nspins))
         DO ispin = 1, nspins
            NULLIFY (v_rspace(ispin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, v_rspace(ispin)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_zero(v_rspace(ispin)%pw)
         END DO
      END IF

      evhxc = 0.0_dp
      CALL qs_rho_get(rho, rho_r=rho_r)
      IF (ASSOCIATED(v_tau_rspace)) THEN
         CALL qs_rho_get(rho, tau_r=tau_r)
      END IF
      DO ispin = 1, nspins
         ! Add v_hartree + v_xc = v_rspace
         CALL pw_scale(v_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(ec_env%vh_rspace%pw, v_rspace(ispin)%pw)
         ! integrate over potential <a|V|b>
         CALL integrate_v_rspace(v_rspace=v_rspace(ispin), &
                                 hmat=ec_env%matrix_ks(ispin, 1), &
                                 qs_env=qs_env, &
                                 calculate_forces=.FALSE., &
                                 basis_type="HARRIS", &
                                 task_list_external=ec_env%task_list)

         IF (ASSOCIATED(v_tau_rspace)) THEN
            ! integrate over Tau-potential <nabla.a|V|nabla.b>
            CALL pw_scale(v_tau_rspace(ispin)%pw, v_tau_rspace(ispin)%pw%pw_grid%dvol)
            CALL integrate_v_rspace(v_rspace=v_tau_rspace(ispin), hmat=ec_env%matrix_ks(ispin, 1), &
                                    qs_env=qs_env, calculate_forces=.FALSE., compute_tau=.TRUE., &
                                    basis_type="HARRIS", &
                                    task_list_external=ec_env%task_list)
         END IF

         ! calclulate Int(vhxc*rho)dr and Int(vtau*tau)dr
         evhxc = evhxc + pw_integral_ab(rho_r(ispin)%pw, v_rspace(ispin)%pw)/ &
                 v_rspace(1)%pw%pw_grid%dvol
         IF (ASSOCIATED(v_tau_rspace)) THEN
            evhxc = evhxc + pw_integral_ab(tau_r(ispin)%pw, v_tau_rspace(ispin)%pw)/ &
                    v_tau_rspace(ispin)%pw%pw_grid%dvol
         END IF

      END DO

      ! return pw grids
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
         IF (ASSOCIATED(v_tau_rspace)) THEN
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
         END IF
      END DO
      DEALLOCATE (v_rspace)
      IF (ASSOCIATED(v_tau_rspace)) DEALLOCATE (v_tau_rspace)

      ! energies
      ec_env%exc = eexc
      ec_env%vhxc = evhxc

      ! add the core matrix
      DO ispin = 1, nspins
         CALL dbcsr_add(ec_env%matrix_ks(ispin, 1)%matrix, ec_env%matrix_h(1, 1)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         CALL dbcsr_filter(ec_env%matrix_ks(ispin, 1)%matrix, &
                           dft_control%qs_control%eps_filter_matrix)
      END DO

      CALL timestop(handle)

   END SUBROUTINE ec_build_ks_matrix

! **************************************************************************************************
!> \brief Construction of the Core Hamiltonian Matrix
!>        Short version of qs_core_hamiltonian
!> \param qs_env ...
!> \param ec_env ...
!> \author Creation (03.2014,JGH)
! **************************************************************************************************
   SUBROUTINE ec_build_core_hamiltonian_force(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_build_core_hamiltonian_force'

      INTEGER                                            :: handle, iounit, nder, nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: calculate_forces, use_virial
      REAL(KIND=dp)                                      :: eps_ppnl, fconv
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: stdeb, sttot
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: scrm
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      END IF

      IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "ec_build_core_hamiltonian_force - START"

      calculate_forces = .TRUE.

      ! no k-points possible
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, dft_control=dft_control)
      nimages = dft_control%nimages
      IF (nimages /= 1) THEN
         CPABORT("K-points for Harris functional not implemented")
      END IF

      ! check for virial (currently no stress tensor available)
      CALL get_qs_env(qs_env=qs_env, cell=cell, virial=virial)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      fconv = 1.0E-9_dp*pascal/cell%deth
      IF (debug_stress .AND. use_virial) THEN
         sttot = virial%pv_virial
      END IF

      ! check for GAPW/GAPW_XC
      IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) THEN
         CPABORT("Harris functional for GAPW not implemented")
      END IF

      ! get neighbor lists, we need the full sab_orb list from the ec_env
      NULLIFY (sab_orb, sac_ppl, sap_ppnl)
      sab_orb => ec_env%sab_orb
      sac_ppl => ec_env%sac_ppl
      sap_ppnl => ec_env%sap_ppnl

      ! initialize src matrix
      NULLIFY (scrm)
      CALL dbcsr_allocate_matrix_set(scrm, 1, 1)
      ALLOCATE (scrm(1, 1)%matrix)
      CALL dbcsr_create(scrm(1, 1)%matrix, template=ec_env%matrix_s(1, 1)%matrix)
      CALL cp_dbcsr_alloc_block_from_nbl(scrm(1, 1)%matrix, sab_orb)

      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, force=force)
      nder = 1
      IF (SIZE(ec_env%matrix_p, 1) == 2) THEN
         CALL dbcsr_add(ec_env%matrix_p(1, 1)%matrix, ec_env%matrix_p(2, 1)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         CALL dbcsr_add(ec_env%matrix_w(1, 1)%matrix, ec_env%matrix_w(2, 1)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF

      ! Overlap and kinetic energy matrices
      IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_overlap
      CALL build_overlap_matrix(ks_env, matrixkp_s=scrm, &
                                !CALL build_overlap_matrix(ks_env, matrixkp_s=ec_env%matrix_s, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="HARRIS", &
                                basis_type_b="HARRIS", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrixkp_p=ec_env%matrix_w)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%overlap(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Wout*dS    ", fodeb
         fodeb(1:3) = force(1)%kinetic(1:3, 1)
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_overlap - stdeb)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| Wout*dS', one_third_sum_diag(stdeb), det_3x3(stdeb)
         stdeb = virial%pv_ekinetic
      END IF
      CALL build_kinetic_matrix(ks_env, matrixkp_t=scrm, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="HARRIS", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrixkp_p=ec_env%matrix_p)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%kinetic(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pout*dT    ", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_ekinetic - stdeb)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| Pout*dT', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF
      IF (SIZE(ec_env%matrix_p, 1) == 2) THEN
         CALL dbcsr_add(ec_env%matrix_p(1, 1)%matrix, ec_env%matrix_p(2, 1)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      END IF

      ! compute the ppl contribution to the core hamiltonian
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)

      IF (ASSOCIATED(sac_ppl)) THEN
         IF (calculate_forces .AND. debug_forces) fodeb(1:3) = force(1)%gth_ppl(1:3, 1)
         IF (debug_stress .AND. use_virial) stdeb = virial%pv_ppl
         CALL build_core_ppl(scrm, ec_env%matrix_p, force, &
                             virial, calculate_forces, use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ppl, &
                             nimages, cell_to_index, "HARRIS")
         IF (calculate_forces .AND. debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppl(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pout*dH_PPL ", fodeb
         END IF
         IF (debug_stress .AND. use_virial) THEN
            stdeb = fconv*(virial%pv_ppl - stdeb)
            CALL mp_sum(stdeb, para_env%group)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| Pout*dH_PPL', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF
      END IF

      ! compute the ppnl contribution to the core hamiltonian ***
      eps_ppnl = dft_control%qs_control%eps_ppnl
      IF (ASSOCIATED(sap_ppnl)) THEN
         IF (calculate_forces .AND. debug_forces) fodeb(1:3) = force(1)%gth_ppnl(1:3, 1)
         IF (debug_stress .AND. use_virial) stdeb = virial%pv_ppnl
         CALL build_core_ppnl(scrm, ec_env%matrix_p, force, &
                              virial, calculate_forces, use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, &
                              sab_orb, sap_ppnl, eps_ppnl, &
                              nimages, cell_to_index, "HARRIS")
         IF (calculate_forces .AND. debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppnl(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pout*dH_PPNL", fodeb
         END IF
         IF (debug_stress .AND. use_virial) THEN
            stdeb = fconv*(virial%pv_ppnl - stdeb)
            CALL mp_sum(stdeb, para_env%group)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| Pout*dH_PPNL', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF
      END IF

      ! External field (nonperiodic case)
      ec_env%efield_nuclear = 0.0_dp
      IF (calculate_forces .AND. debug_forces) fodeb(1:3) = force(1)%efield(1:3, 1)
      CALL ec_efield_local_operator(qs_env, ec_env, calculate_forces)
      IF (calculate_forces .AND. debug_forces) THEN
         fodeb(1:3) = force(1)%efield(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pout*dEfield", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_virial - sttot)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| Stress Pout*dHcore   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") ' '
      END IF

      ! delete scr matrix
      CALL dbcsr_deallocate_matrix_set(scrm)

      CALL timestop(handle)

   END SUBROUTINE ec_build_core_hamiltonian_force

! **************************************************************************************************
!> \brief Solve KS equation for a given matrix
!> \brief calculate the complete KS matrix
!> \param qs_env ...
!> \param ec_env ...
!> \par History
!>      03.2014 adapted from qs_ks_build_kohn_sham_matrix [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ec_build_ks_matrix_force(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_build_ks_matrix_force'

      INTEGER                                            :: handle, i, iounit, ispin, natom, nspins
      LOGICAL                                            :: use_virial
      REAL(dp)                                           :: dehartree, eexc, ehartree, eovrl, exc, &
                                                            fconv
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: ftot
      REAL(dp), DIMENSION(3)                             :: fodeb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_stress, pv_loc, stdeb, sttot
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type)                                 :: scrm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: dv_hartree_rspace, rho_tot_gspace, &
                                                            rhodn_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace, vtot_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r, rhoout_g, rhoout_r, tau_r, &
                                                            tauout_g, tauout_r, v_rspace, &
                                                            v_tau_rspace, v_xc, v_xc_tau
      TYPE(pw_p_type), POINTER                           :: rho_core
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(virial_type), POINTER                         :: virdeb, virial

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      END IF

      ! get all information on the electronic density
      NULLIFY (atomic_kind_set, blacs_env, cell, dft_control, force, ks_env, &
               matrix_p, matrix_s, para_env, rho, sab_orb, virial)
      CALL get_qs_env(qs_env=qs_env, &
                      blacs_env=blacs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      ks_env=ks_env, &
                      para_env=para_env, &
                      rho=rho, &
                      sab_orb=sab_orb, &
                      virial=virial)

      nspins = dft_control%nspins
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      fconv = 1.0E-9_dp*pascal/cell%deth
      IF (debug_stress .AND. use_virial) THEN
         sttot = virial%pv_virial
      END IF

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CPASSERT(ASSOCIATED(pw_env))

      NULLIFY (auxbas_pw_pool, poisson_env)
      ! gets the tmp grids
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      ! Calculate the Hartree potential
      NULLIFY (v_hartree_gspace%pw, rhodn_tot_gspace%pw, v_hartree_rspace%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rhodn_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL pw_transfer(ec_env%vh_rspace%pw, v_hartree_rspace%pw)

      CALL get_qs_env(qs_env=qs_env, force=force)
      ! calculate output density on grid
      ! rho_in(R):   CALL qs_rho_get(rho, rho_r=rho_r)
      ! rho_in(G):   CALL qs_rho_get(rho, rho_g=rho_g)
      CALL qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g, tau_r=tau_r)
      NULLIFY (rhoout_r, rhoout_g)
      ALLOCATE (rhoout_r(nspins), rhoout_g(nspins))
      DO ispin = 1, nspins
         NULLIFY (rhoout_r(ispin)%pw, rhoout_g(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoout_r(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoout_g(ispin)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      END DO
      CALL pw_pool_create_pw(auxbas_pw_pool, dv_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, vtot_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL pw_zero(rhodn_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL calculate_rho_elec(ks_env=ks_env, matrix_p=ec_env%matrix_p(ispin, 1)%matrix, &
                                 rho=rhoout_r(ispin), &
                                 rho_gspace=rhoout_g(ispin), &
                                 basis_type="HARRIS", &
                                 task_list_external=ec_env%task_list)
      END DO

      ! Save Harris on real space grid for use in properties
      ALLOCATE (ec_env%rhoout_r(nspins))
      DO ispin = 1, nspins
         NULLIFY (ec_env%rhoout_r(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, ec_env%rhoout_r(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_copy(rhoout_r(ispin)%pw, ec_env%rhoout_r(ispin)%pw)
      END DO

      NULLIFY (tauout_r, tauout_g)
      IF (dft_control%use_kinetic_energy_density) THEN
         IF (use_virial) &
            CALL cp_abort(__LOCATION__, "Stress tensor of the Harris functional with "// &
                          "tau-dependent functionals not implemented!")
         ALLOCATE (tauout_r(nspins), tauout_g(nspins))
         DO ispin = 1, nspins
            NULLIFY (tauout_r(ispin)%pw, tauout_g(ispin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, tauout_r(ispin)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_pool_create_pw(auxbas_pw_pool, tauout_g(ispin)%pw, &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         END DO

         DO ispin = 1, nspins
            CALL calculate_rho_elec(ks_env=ks_env, matrix_p=ec_env%matrix_p(ispin, 1)%matrix, &
                                    rho=tauout_r(ispin), &
                                    rho_gspace=tauout_g(ispin), &
                                    total_rho=total_rhoout, &
                                    compute_tau=.TRUE., &
                                    basis_type="HARRIS", &
                                    task_list_external=ec_env%task_list)
         END DO

         ! Save Harris on real space grid for use in properties
         ALLOCATE (ec_env%tauout_r(nspins))
         DO ispin = 1, nspins
            NULLIFY (ec_env%tauout_r(ispin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, ec_env%tauout_r(ispin)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_copy(tauout_r(ispin)%pw, ec_env%tauout_r(ispin)%pw)
         END DO
      END IF

      IF (use_virial) THEN

         ! Calculate the Hartree potential
         NULLIFY (rho_tot_gspace%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

         ! Get the total input density in g-space [ions + electrons]
         CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

         ! make rho_tot_gspace with output density
         CALL get_qs_env(qs_env=qs_env, rho_core=rho_core)
         CALL pw_copy(rho_core%pw, rhodn_tot_gspace%pw)
         DO ispin = 1, dft_control%nspins
            CALL pw_axpy(rhoout_g(ispin)%pw, rhodn_tot_gspace%pw)
         END DO

         ! Volume and Green function terms
         h_stress(:, :) = 0.0_dp
         CALL pw_poisson_solve(poisson_env, &
                               density=rho_tot_gspace%pw, &  ! n_in
                               ehartree=ehartree, &
                               vhartree=v_hartree_gspace%pw, & ! v_H[n_in]
                               h_stress=h_stress, &
                               aux_density=rhodn_tot_gspace%pw) ! n_out

         virial%pv_ehartree = virial%pv_ehartree + h_stress/REAL(para_env%num_pe, dp)
         virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env%num_pe, dp)

         IF (debug_stress) THEN
            stdeb = fconv*(h_stress/REAL(para_env%num_pe, dp))
            CALL mp_sum(stdeb, para_env%group)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| GREEN 1st v_H[n_in]*n_out  ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF

         ! activate stress calculation
         virial%pv_calculate = .TRUE.

         NULLIFY (v_rspace, v_tau_rspace)
         CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=ec_env%xc_section, &
                            vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE.)

         ! Stress tensor XC-functional GGA contribution
         virial%pv_exc = virial%pv_exc - virial%pv_xc
         virial%pv_virial = virial%pv_virial - virial%pv_xc

         IF (debug_stress) THEN
            stdeb = -1.0_dp*fconv*virial%pv_xc
            CALL mp_sum(stdeb, para_env%group)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| GGA 1st E_xc[Pin]   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF

         IF (ASSOCIATED(v_rspace)) THEN
            DO ispin = 1, nspins
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
            END DO
            DEALLOCATE (v_rspace)
         END IF
         IF (ASSOCIATED(v_tau_rspace)) THEN
            DO ispin = 1, nspins
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
            END DO
            DEALLOCATE (v_tau_rspace)
         END IF
         CALL pw_zero(rhodn_tot_gspace%pw)

      END IF

      ! rho_out - rho_in
      DO ispin = 1, nspins
         CALL pw_axpy(rho_r(ispin)%pw, rhoout_r(ispin)%pw, -1.0_dp)
         CALL pw_axpy(rho_g(ispin)%pw, rhoout_g(ispin)%pw, -1.0_dp)
         CALL pw_axpy(rhoout_g(ispin)%pw, rhodn_tot_gspace%pw)
         IF (dft_control%use_kinetic_energy_density) CALL pw_axpy(tau_r(ispin)%pw, tauout_r(ispin)%pw, -1.0_dp)
      END DO
      ! calculate associated hartree potential
      IF (use_virial) THEN

         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
         NULLIFY (v_hartree_gspace%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

         ! Stress tensor - 2nd derivative Volume and Green function contribution
         h_stress(:, :) = 0.0_dp
         CALL pw_poisson_solve(poisson_env, &
                               density=rhodn_tot_gspace%pw, &  ! delta_n
                               ehartree=dehartree, &
                               vhartree=v_hartree_gspace%pw, & ! v_H[delta_n]
                               h_stress=h_stress, &
                               aux_density=rho_tot_gspace%pw)  ! n_in

         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)

         virial%pv_ehartree = virial%pv_ehartree + h_stress/REAL(para_env%num_pe, dp)
         virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env%num_pe, dp)

         IF (debug_stress) THEN
            stdeb = fconv*(h_stress/REAL(para_env%num_pe, dp))
            CALL mp_sum(stdeb, para_env%group)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| GREEN 2nd V_H[dP]*n_in  ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF

      ELSE
         ! v_H[dn]
         CALL pw_poisson_solve(poisson_env, rhodn_tot_gspace%pw, dehartree, &
                               v_hartree_gspace%pw)
      END IF

      CALL pw_transfer(v_hartree_gspace%pw, dv_hartree_rspace%pw)
      CALL pw_scale(dv_hartree_rspace%pw, dv_hartree_rspace%pw%pw_grid%dvol)
      ! Getting nuclear force contribution from the core charge density
      ! Vh(rho_in + rho_c) + Vh(rho_out - rho_in)
      CALL pw_transfer(v_hartree_rspace%pw, vtot_rspace%pw)
      CALL pw_axpy(dv_hartree_rspace%pw, vtot_rspace%pw)
      IF (debug_forces) fodeb(1:3) = force(1)%rho_core(1:3, 1)
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_ehartree
      CALL integrate_v_core_rspace(vtot_rspace, qs_env)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_core(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Vtot*dncore", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_ehartree - stdeb)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| Vtot*dncore', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF
      !
      ! Pulay force from Tr P_in (V_H(drho)+ Fxc(rho_in)*drho)
      ! RHS of CPKS equations: (V_H(drho)+ Fxc(rho_in)*drho)*C0
      ! Fxc*drho term
      NULLIFY (v_xc)
      xc_section => ec_env%xc_section

      IF (use_virial) virial%pv_xc = 0.0_dp
      NULLIFY (v_xc, v_xc_tau)
      CALL create_kernel(qs_env, &
                         vxc=v_xc, &
                         vxc_tau=v_xc_tau, &
                         rho=rho, &
                         rho1_r=rhoout_r, &
                         rho1_g=rhoout_g, &
                         tau1_r=tauout_r, &
                         xc_section=xc_section, &
                         compute_virial=use_virial, &
                         virial_xc=virial%pv_xc)

      IF (use_virial) THEN
         ! Stress-tensor XC-functional 2nd GGA terms
         virial%pv_exc = virial%pv_exc + virial%pv_xc
         virial%pv_virial = virial%pv_virial + virial%pv_xc
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = 1.0_dp*fconv*virial%pv_xc
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| GGA 2nd f_Hxc[dP]*Pin   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF
      !
      CALL get_qs_env(qs_env=qs_env, rho=rho, matrix_s_kp=matrix_s)
      NULLIFY (ec_env%matrix_hz)
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_hz, nspins)
      DO ispin = 1, nspins
         ALLOCATE (ec_env%matrix_hz(ispin)%matrix)
         CALL dbcsr_create(ec_env%matrix_hz(ispin)%matrix, template=matrix_s(1, 1)%matrix)
         CALL dbcsr_copy(ec_env%matrix_hz(ispin)%matrix, matrix_s(1, 1)%matrix)
         CALL dbcsr_set(ec_env%matrix_hz(ispin)%matrix, 0.0_dp)
      END DO
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ! vtot = v_xc(ispin) + dv_hartree
      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_virial

      ! Stress-tensor 2nd derivative integral contribution
      IF (use_virial) THEN
         pv_loc = virial%pv_virial
      END IF

      DO ispin = 1, nspins
         CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(dv_hartree_rspace%pw, v_xc(ispin)%pw)
         CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                 hmat=ec_env%matrix_hz(ispin), &
                                 pmat=matrix_p(ispin, 1), &
                                 calculate_forces=.TRUE.)
      END DO

      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dKdrho", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_virial - stdeb)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| INT 2nd f_Hxc[dP]*Pin    ', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF

      IF (ASSOCIATED(v_xc_tau)) THEN
         IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
         IF (debug_stress .AND. use_virial) stdeb = virial%pv_virial

         DO ispin = 1, nspins
            CALL pw_scale(v_xc_tau(ispin)%pw, v_xc_tau(ispin)%pw%pw_grid%dvol)
            CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc_tau(ispin), &
                                    hmat=ec_env%matrix_hz(ispin), &
                                    pmat=matrix_p(ispin, 1), &
                                    compute_tau=.TRUE., &
                                    calculate_forces=.TRUE.)
         END DO

         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dKtaudtau", fodeb
         END IF
         IF (debug_stress .AND. use_virial) THEN
            stdeb = fconv*(virial%pv_virial - stdeb)
            CALL mp_sum(stdeb, para_env%group)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| INT 2nd f_xctau[dP]*Pin    ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF
      END IF
      ! Stress-tensor 2nd derivative integral contribution
      IF (use_virial) THEN
         virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
      END IF

      ! initialize src matrix
      NULLIFY (scrm%matrix)
      ALLOCATE (scrm%matrix)
      CALL dbcsr_create(scrm%matrix, template=ec_env%matrix_s(1, 1)%matrix)
      CALL cp_dbcsr_alloc_block_from_nbl(scrm%matrix, ec_env%sab_orb)

      ! v_rspace and v_tau_rspace are generated from the auxbas pool
      NULLIFY (v_rspace, v_tau_rspace)

      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=ec_env%xc_section, &
                         vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=eexc, just_energy=.FALSE.)

      IF (use_virial) THEN
         eexc = 0.0_dp
         IF (ASSOCIATED(v_rspace)) THEN
            DO ispin = 1, nspins
               ! 2nd deriv xc-volume term
               eexc = eexc + pw_integral_ab(rhoout_r(ispin)%pw, v_rspace(ispin)%pw)
            END DO
         END IF
         IF (ASSOCIATED(v_tau_rspace)) THEN
            DO ispin = 1, nspins
               ! 2nd deriv xc-volume term
               eexc = eexc + pw_integral_ab(tauout_r(ispin)%pw, v_tau_rspace(ispin)%pw)
            END DO
         END IF
      END IF

      IF (.NOT. ASSOCIATED(v_rspace)) THEN
         ALLOCATE (v_rspace(nspins))
         DO ispin = 1, nspins
            NULLIFY (v_rspace(ispin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, v_rspace(ispin)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_zero(v_rspace(ispin)%pw)
         END DO
      END IF

      ! Stress-tensor contribution derivative of integrand
      ! int v_Hxc[n^în]*n^out
      IF (use_virial) THEN
         pv_loc = virial%pv_virial
      END IF

      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_virial
      DO ispin = 1, nspins
         ! Add v_hartree + v_xc = v_rspace
         CALL pw_scale(v_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(v_hartree_rspace%pw, v_rspace(ispin)%pw)
         ! integrate over potential <a|V|b>
         CALL integrate_v_rspace(v_rspace=v_rspace(ispin), &
                                 hmat=scrm, &
                                 pmat=ec_env%matrix_p(ispin, 1), &
                                 qs_env=qs_env, &
                                 calculate_forces=.TRUE., &
                                 basis_type="HARRIS", &
                                 task_list_external=ec_env%task_list)

         IF (ASSOCIATED(v_tau_rspace)) THEN
            ! integrate over Tau-potential <nabla.a|V|nabla.b>
            CALL pw_scale(v_tau_rspace(ispin)%pw, v_tau_rspace(ispin)%pw%pw_grid%dvol)
            CALL integrate_v_rspace(v_rspace=v_tau_rspace(ispin), hmat=scrm, &
                                    pmat=ec_env%matrix_p(ispin, 1), &
                                    qs_env=qs_env, calculate_forces=.TRUE., compute_tau=.TRUE., &
                                    basis_type="HARRIS", &
                                    task_list_external=ec_env%task_list)
         END IF

      END DO
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pout*dVhxc ", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_virial - stdeb)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| INT Pout*dVhxc   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF

      ! Stress-tensor
      IF (use_virial) THEN
         virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
      END IF

      ! delete scr matrix
      CALL dbcsr_release(scrm%matrix)
      DEALLOCATE (scrm%matrix)

      ! return pw grids
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
         IF (ASSOCIATED(v_tau_rspace)) THEN
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
         END IF
      END DO
      IF (ASSOCIATED(v_tau_rspace)) DEALLOCATE (v_tau_rspace)

      ! Core overlap
      IF (debug_forces) fodeb(1:3) = force(1)%core_overlap(1:3, 1)
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_ecore_overlap
      CALL calculate_ecore_overlap(qs_env, para_env, .TRUE., E_overlap_core=eovrl)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%core_overlap(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: CoreOverlap", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(stdeb - virial%pv_ecore_overlap)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| CoreOverlap   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF

      IF (debug_forces) THEN
         CALL get_qs_env(qs_env, natom=natom, atomic_kind_set=atomic_kind_set)
         ALLOCATE (ftot(3, natom))
         CALL total_qs_force(ftot, force, atomic_kind_set)
         fodeb(1:3) = ftot(1:3, 1)
         DEALLOCATE (ftot)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Force Explicit", fodeb
      END IF

      DEALLOCATE (v_rspace)
      !
      CALL pw_pool_give_back_pw(auxbas_pw_pool, dv_hartree_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, vtot_rspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoout_r(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoout_g(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      END DO
      DEALLOCATE (rhoout_r, rhoout_g, v_xc)
      IF (ASSOCIATED(tauout_r)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, tauout_r(ispin)%pw)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, tauout_g(ispin)%pw)
         END DO
         DEALLOCATE (tauout_r, tauout_g)
      END IF
      IF (ASSOCIATED(v_xc_tau)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc_tau(ispin)%pw)
         END DO
         DEALLOCATE (v_xc_tau)
      END IF
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhodn_tot_gspace%pw)

      ! Stress tensor - volume terms need to be stored,
      ! for a sign correction in QS at the end of qs_force
      IF (use_virial) THEN
         IF (qs_env%energy_correction) THEN
            ec_env%ehartree = ehartree + dehartree
            ec_env%exc = exc + eexc
         END IF
      END IF

      IF (debug_stress .AND. use_virial) THEN
         ! In total: -1.0*E_H
         stdeb = -1.0_dp*fconv*ehartree
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| VOL 1st v_H[n_in]*n_out', one_third_sum_diag(stdeb), det_3x3(stdeb)

         stdeb = -1.0_dp*fconv*exc
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| VOL 1st E_XC[n_in]', one_third_sum_diag(stdeb), det_3x3(stdeb)

         stdeb = -1.0_dp*fconv*dehartree
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| VOL 2nd v_H[dP]*n_in', one_third_sum_diag(stdeb), det_3x3(stdeb)

         stdeb = -1.0_dp*fconv*eexc
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| VOL 2nd v_XC[n_in]*dP', one_third_sum_diag(stdeb), det_3x3(stdeb)

         ! For debugging, create a second virial environment,
         ! apply volume terms immediately
         CALL virial_create(virdeb)
         CALL cp_virial(virial, virdeb)

         CALL mp_sum(virdeb%pv_overlap, para_env%group)
         CALL mp_sum(virdeb%pv_ekinetic, para_env%group)
         CALL mp_sum(virdeb%pv_ppl, para_env%group)
         CALL mp_sum(virdeb%pv_ppnl, para_env%group)
         CALL mp_sum(virdeb%pv_ecore_overlap, para_env%group)
         CALL mp_sum(virdeb%pv_ehartree, para_env%group)
         CALL mp_sum(virdeb%pv_exc, para_env%group)
         CALL mp_sum(virdeb%pv_exx, para_env%group)
         CALL mp_sum(virdeb%pv_vdw, para_env%group)
         CALL mp_sum(virdeb%pv_mp2, para_env%group)
         CALL mp_sum(virdeb%pv_nlcc, para_env%group)
         CALL mp_sum(virdeb%pv_gapw, para_env%group)
         CALL mp_sum(virdeb%pv_lrigpw, para_env%group)
         CALL mp_sum(virdeb%pv_virial, para_env%group)
         CALL symmetrize_virial(virdeb)

         ! apply stress-tensor 1st and 2nd volume terms
         DO i = 1, 3
            virdeb%pv_ehartree(i, i) = virdeb%pv_ehartree(i, i) - 2.0_dp*(ehartree + dehartree)
            virdeb%pv_virial(i, i) = virdeb%pv_virial(i, i) - exc - eexc &
                                     - 2.0_dp*(ehartree + dehartree)
            virdeb%pv_exc(i, i) = virdeb%pv_exc(i, i) - exc - eexc
            ! The factor 2 is a hack. It compensates the plus sign in h_stress/pw_poisson_solve.
            ! The sign in pw_poisson_solve is correct for FIST, but not for QS.
            ! There should be a more elegant solution to that ...
         END DO

         CALL mp_sum(sttot, para_env%group)
         stdeb = fconv*(virdeb%pv_virial - sttot)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| Explicit electronic stress   ', one_third_sum_diag(stdeb), det_3x3(stdeb)

         stdeb = fconv*(virdeb%pv_virial)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| Explicit total stress   ', one_third_sum_diag(stdeb), det_3x3(stdeb)

         IF (iounit > 0) WRITE (iounit, *) "Harris_ks_force - END"
         CALL write_stress_tensor_components(virdeb, iounit, cell)
         CALL write_stress_tensor(virdeb%pv_virial, iounit, cell, .FALSE.)

         CALL virial_release(virdeb)

      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_build_ks_matrix_force

! **************************************************************************************************
!> \brief Solve KS equation for a given matrix
!> \param qs_env ...
!> \param ec_env ...
!> \par History
!>      03.2014 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ec_ks_solver(qs_env, ec_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ec_ks_solver'

      CHARACTER(LEN=default_string_length)               :: headline
      INTEGER                                            :: handle, ispin, nspins
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ksmat, pmat, smat, wmat
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      nspins = dft_control%nspins

      ! create density matrix
      IF (.NOT. ASSOCIATED(ec_env%matrix_p)) THEN
         headline = "DENSITY MATRIX"
         CALL dbcsr_allocate_matrix_set(ec_env%matrix_p, nspins, 1)
         DO ispin = 1, nspins
            ALLOCATE (ec_env%matrix_p(ispin, 1)%matrix)
            CALL dbcsr_create(ec_env%matrix_p(ispin, 1)%matrix, name=TRIM(headline), &
                              template=ec_env%matrix_s(1, 1)%matrix)
            CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_p(ispin, 1)%matrix, ec_env%sab_orb)
         END DO
      END IF
      ! create energy weighted density matrix
      IF (.NOT. ASSOCIATED(ec_env%matrix_w)) THEN
         headline = "ENERGY WEIGHTED DENSITY MATRIX"
         CALL dbcsr_allocate_matrix_set(ec_env%matrix_w, nspins, 1)
         DO ispin = 1, nspins
            ALLOCATE (ec_env%matrix_w(ispin, 1)%matrix)
            CALL dbcsr_create(ec_env%matrix_w(ispin, 1)%matrix, name=TRIM(headline), &
                              template=ec_env%matrix_s(1, 1)%matrix)
            CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_w(ispin, 1)%matrix, ec_env%sab_orb)
         END DO
      END IF

      IF (ec_env%mao) THEN
         CALL mao_create_matrices(ec_env, ksmat, smat, pmat, wmat)
      ELSE
         ksmat => ec_env%matrix_ks
         smat => ec_env%matrix_s
         pmat => ec_env%matrix_p
         wmat => ec_env%matrix_w
      END IF

      SELECT CASE (ec_env%ks_solver)
      CASE (ec_diagonalization)
         CALL ec_diag_solver(qs_env, ksmat, smat, pmat, wmat)
      CASE (ec_ot_diag)
         CALL ec_ot_diag_solver(qs_env, ec_env, ksmat, smat, pmat, wmat)
      CASE (ec_matrix_sign, ec_matrix_trs4, ec_matrix_tc2)
         CALL ec_ls_init(qs_env, ksmat, smat)
         CALL ec_ls_solver(qs_env, pmat, wmat, ec_ls_method=ec_env%ks_solver)
      CASE DEFAULT
         CPASSERT(.FALSE.)
      END SELECT

      IF (ec_env%mao) THEN
         CALL mao_release_matrices(ec_env, ksmat, smat, pmat, wmat)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_ks_solver

! **************************************************************************************************
!> \brief Create matrices with MAO sizes
!> \param ec_env ...
!> \param ksmat ...
!> \param smat ...
!> \param pmat ...
!> \param wmat ...
!> \par History
!>      08.2016 created [JGH]
!> \author JGH
! **************************************************************************************************

   SUBROUTINE mao_create_matrices(ec_env, ksmat, smat, pmat, wmat)

      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ksmat, smat, pmat, wmat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mao_create_matrices'

      INTEGER                                            :: handle, ispin, nspins
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coef
      TYPE(dbcsr_type)                                   :: cgmat

      CALL timeset(routineN, handle)

      mao_coef => ec_env%mao_coef

      NULLIFY (ksmat, smat, pmat, wmat)
      nspins = SIZE(ec_env%matrix_ks, 1)
      CALL dbcsr_get_info(mao_coef(1)%matrix, col_blk_size=col_blk_sizes, distribution=dbcsr_dist)
      CALL dbcsr_allocate_matrix_set(ksmat, nspins, 1)
      CALL dbcsr_allocate_matrix_set(smat, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (ksmat(ispin, 1)%matrix)
         CALL dbcsr_create(ksmat(ispin, 1)%matrix, dist=dbcsr_dist, name="MAO KS mat", &
                           matrix_type=dbcsr_type_symmetric, row_blk_size=col_blk_sizes, &
                           col_blk_size=col_blk_sizes, nze=0)
         ALLOCATE (smat(ispin, 1)%matrix)
         CALL dbcsr_create(smat(ispin, 1)%matrix, dist=dbcsr_dist, name="MAO S mat", &
                           matrix_type=dbcsr_type_symmetric, row_blk_size=col_blk_sizes, &
                           col_blk_size=col_blk_sizes, nze=0)
      END DO
      !
      CALL dbcsr_create(cgmat, name="TEMP matrix", template=mao_coef(1)%matrix)
      DO ispin = 1, nspins
         CALL dbcsr_multiply("N", "N", 1.0_dp, ec_env%matrix_s(1, 1)%matrix, mao_coef(ispin)%matrix, &
                             0.0_dp, cgmat)
         CALL dbcsr_multiply("T", "N", 1.0_dp, mao_coef(ispin)%matrix, cgmat, 0.0_dp, smat(ispin, 1)%matrix)
         CALL dbcsr_multiply("N", "N", 1.0_dp, ec_env%matrix_ks(1, 1)%matrix, mao_coef(ispin)%matrix, &
                             0.0_dp, cgmat)
         CALL dbcsr_multiply("T", "N", 1.0_dp, mao_coef(ispin)%matrix, cgmat, 0.0_dp, ksmat(ispin, 1)%matrix)
      END DO
      CALL dbcsr_release(cgmat)

      CALL dbcsr_allocate_matrix_set(pmat, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (pmat(ispin, 1)%matrix)
         CALL dbcsr_create(pmat(ispin, 1)%matrix, template=smat(1, 1)%matrix, name="MAO P mat")
         CALL cp_dbcsr_alloc_block_from_nbl(pmat(ispin, 1)%matrix, ec_env%sab_orb)
      END DO

      CALL dbcsr_allocate_matrix_set(wmat, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (wmat(ispin, 1)%matrix)
         CALL dbcsr_create(wmat(ispin, 1)%matrix, template=smat(1, 1)%matrix, name="MAO W mat")
         CALL cp_dbcsr_alloc_block_from_nbl(wmat(ispin, 1)%matrix, ec_env%sab_orb)
      END DO

      CALL timestop(handle)

   END SUBROUTINE mao_create_matrices

! **************************************************************************************************
!> \brief Release matrices with MAO sizes
!> \param ec_env ...
!> \param ksmat ...
!> \param smat ...
!> \param pmat ...
!> \param wmat ...
!> \par History
!>      08.2016 created [JGH]
!> \author JGH
! **************************************************************************************************

   SUBROUTINE mao_release_matrices(ec_env, ksmat, smat, pmat, wmat)

      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ksmat, smat, pmat, wmat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mao_release_matrices'

      INTEGER                                            :: handle, ispin, nspins
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coef
      TYPE(dbcsr_type)                                   :: cgmat

      CALL timeset(routineN, handle)

      mao_coef => ec_env%mao_coef
      nspins = SIZE(mao_coef, 1)

      ! save pmat/wmat in full basis format
      CALL dbcsr_create(cgmat, name="TEMP matrix", template=mao_coef(1)%matrix)
      DO ispin = 1, nspins
         CALL dbcsr_multiply("N", "N", 1.0_dp, mao_coef(ispin)%matrix, pmat(ispin, 1)%matrix, 0.0_dp, cgmat)
         CALL dbcsr_multiply("N", "T", 1.0_dp, mao_coef(ispin)%matrix, cgmat, 0.0_dp, &
                             ec_env%matrix_p(ispin, 1)%matrix, retain_sparsity=.TRUE.)
         CALL dbcsr_multiply("N", "N", 1.0_dp, mao_coef(ispin)%matrix, wmat(ispin, 1)%matrix, 0.0_dp, cgmat)
         CALL dbcsr_multiply("N", "T", 1.0_dp, mao_coef(ispin)%matrix, cgmat, 0.0_dp, &
                             ec_env%matrix_w(ispin, 1)%matrix, retain_sparsity=.TRUE.)
      END DO
      CALL dbcsr_release(cgmat)

      CALL dbcsr_deallocate_matrix_set(ksmat)
      CALL dbcsr_deallocate_matrix_set(smat)
      CALL dbcsr_deallocate_matrix_set(pmat)
      CALL dbcsr_deallocate_matrix_set(wmat)

      CALL timestop(handle)

   END SUBROUTINE mao_release_matrices

! **************************************************************************************************
!> \brief Solve KS equation using diagonalization
!> \param qs_env ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param matrix_p ...
!> \param matrix_w ...
!> \par History
!>      03.2014 created [JGH]
!> \author JGH
! **************************************************************************************************

   SUBROUTINE ec_diag_solver(qs_env, matrix_ks, matrix_s, matrix_p, matrix_w)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s, matrix_p, matrix_w

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ec_diag_solver'

      INTEGER                                            :: handle, ispin, nmo(2), nsize, nspins
      REAL(KIND=dp)                                      :: eps_filter, focc(2)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_ks, fm_mo, fm_ortho
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type), POINTER                          :: buf1_dbcsr, buf2_dbcsr, buf3_dbcsr, &
                                                            ortho_dbcsr, ref_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, para_env)
      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env, para_env=para_env)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      eps_filter = dft_control%qs_control%eps_filter_matrix
      nspins = dft_control%nspins

      nmo = 0
      CALL get_qs_env(qs_env=qs_env, nelectron_spin=nmo)
      focc = 1._dp
      IF (nspins == 1) THEN
         focc = 2._dp
         nmo(1) = nmo(1)/2
      END IF

      CALL dbcsr_get_info(matrix_ks(1, 1)%matrix, nfullrows_total=nsize)
      ALLOCATE (eigenvalues(nsize))

      NULLIFY (fm_ortho, fm_ks, fm_mo, fm_struct, ref_matrix)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nsize, &
                               ncol_global=nsize, para_env=para_env)
      CALL cp_fm_create(fm_ortho, fm_struct)
      CALL cp_fm_create(fm_ks, fm_struct)
      CALL cp_fm_create(fm_mo, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ! factorization
      ref_matrix => matrix_s(1, 1)%matrix
      NULLIFY (ortho_dbcsr, buf1_dbcsr, buf2_dbcsr, buf3_dbcsr)
      CALL dbcsr_init_p(ortho_dbcsr)
      CALL dbcsr_create(ortho_dbcsr, template=ref_matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(buf1_dbcsr)
      CALL dbcsr_create(buf1_dbcsr, template=ref_matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(buf2_dbcsr)
      CALL dbcsr_create(buf2_dbcsr, template=ref_matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(buf3_dbcsr)
      CALL dbcsr_create(buf3_dbcsr, template=ref_matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      ref_matrix => matrix_s(1, 1)%matrix
      CALL copy_dbcsr_to_fm(ref_matrix, fm_ortho)
      CALL cp_fm_cholesky_decompose(fm_ortho)
      CALL cp_fm_triangular_invert(fm_ortho)
      CALL cp_fm_set_all(fm_ks, 0.0_dp)
      CALL cp_fm_to_fm_triangular(fm_ortho, fm_ks, "U")
      CALL copy_fm_to_dbcsr(fm_ks, ortho_dbcsr)
      DO ispin = 1, nspins
         ! calculate ZHZ(T)
         CALL dbcsr_desymmetrize(matrix_ks(ispin, 1)%matrix, buf1_dbcsr)
         CALL dbcsr_multiply("N", "N", 1.0_dp, buf1_dbcsr, ortho_dbcsr, &
                             0.0_dp, buf2_dbcsr, filter_eps=eps_filter)
         CALL dbcsr_multiply("T", "N", 1.0_dp, ortho_dbcsr, buf2_dbcsr, &
                             0.0_dp, buf1_dbcsr, filter_eps=eps_filter)
         ! copy to fm format
         CALL copy_dbcsr_to_fm(buf1_dbcsr, fm_ks)
         CALL choose_eigv_solver(fm_ks, fm_mo, eigenvalues)
         ! back transform of mos c = Z(T)*c
         CALL copy_fm_to_dbcsr(fm_mo, buf1_dbcsr)
         CALL dbcsr_multiply("N", "N", 1.0_dp, ortho_dbcsr, buf1_dbcsr, &
                             0.0_dp, buf2_dbcsr, filter_eps=eps_filter)
         ! density matrix
         CALL dbcsr_set(matrix_p(ispin, 1)%matrix, 0.0_dp)
         CALL dbcsr_multiply("N", "T", focc(ispin), buf2_dbcsr, buf2_dbcsr, &
                             1.0_dp, matrix_p(ispin, 1)%matrix, &
                             retain_sparsity=.TRUE., last_k=nmo(ispin))

         ! energy weighted density matrix
         CALL dbcsr_set(matrix_w(ispin, 1)%matrix, 0.0_dp)
         CALL cp_fm_column_scale(fm_mo, eigenvalues)
         CALL copy_fm_to_dbcsr(fm_mo, buf1_dbcsr)
         CALL dbcsr_multiply("N", "N", 1.0_dp, ortho_dbcsr, buf1_dbcsr, &
                             0.0_dp, buf3_dbcsr, filter_eps=eps_filter)
         CALL dbcsr_multiply("N", "T", focc(ispin), buf2_dbcsr, buf3_dbcsr, &
                             1.0_dp, matrix_w(ispin, 1)%matrix, &
                             retain_sparsity=.TRUE., last_k=nmo(ispin))
      END DO

      CALL cp_fm_release(fm_ks)
      CALL cp_fm_release(fm_mo)
      CALL cp_fm_release(fm_ortho)
      CALL dbcsr_release(ortho_dbcsr)
      CALL dbcsr_release(buf1_dbcsr)
      CALL dbcsr_release(buf2_dbcsr)
      CALL dbcsr_release(buf3_dbcsr)
      DEALLOCATE (ortho_dbcsr, buf1_dbcsr, buf2_dbcsr, buf3_dbcsr)
      DEALLOCATE (eigenvalues)

      CALL timestop(handle)

   END SUBROUTINE ec_diag_solver

! **************************************************************************************************
!> \brief Calculate the energy correction
!> \param ec_env ...
!> \param unit_nr ...
!> \author Creation (03.2014,JGH)
! **************************************************************************************************
   SUBROUTINE ec_energy(ec_env, unit_nr)
      TYPE(energy_correction_type)                       :: ec_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ec_energy'

      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: eband, trace

      CALL timeset(routineN, handle)

      SELECT CASE (ec_env%energy_functional)
      CASE (ec_functional_harris)
         nspins = SIZE(ec_env%matrix_ks, 1)
         eband = 0.0_dp
         DO ispin = 1, nspins

            CALL dbcsr_dot(ec_env%matrix_p(ispin, 1)%matrix, ec_env%matrix_s(1, 1)%matrix, trace)
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T65,F16.10)') 'Tr[PS] ', trace

            CALL dbcsr_dot(ec_env%matrix_ks(ispin, 1)%matrix, ec_env%matrix_p(ispin, 1)%matrix, trace)
            eband = eband + trace
         END DO
         ec_env%eband = eband + ec_env%efield_nuclear
         ec_env%etotal = ec_env%eband + ec_env%ehartree + ec_env%exc - ec_env%vhxc + ec_env%edispersion
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T3,A,T56,F25.15)') "Eband    ", ec_env%eband
            WRITE (unit_nr, '(T3,A,T56,F25.15)') "Ehartree ", ec_env%ehartree
            WRITE (unit_nr, '(T3,A,T56,F25.15)') "Exc      ", ec_env%exc
            WRITE (unit_nr, '(T3,A,T56,F25.15)') "Evhxc    ", ec_env%vhxc
            WRITE (unit_nr, '(T3,A,T56,F25.15)') "Edisp    ", ec_env%edispersion
            WRITE (unit_nr, '(T3,A,T56,F25.15)') "Etotal Harris Functional   ", ec_env%etotal
         END IF

      CASE DEFAULT

         CPASSERT(.FALSE.)

      END SELECT

      CALL timestop(handle)

   END SUBROUTINE ec_energy

! **************************************************************************************************
!> \brief builds either the full neighborlist or neighborlists of molecular
!> \brief subsets, depending on parameter values
!> \param qs_env ...
!> \param ec_env ...
!> \par History
!>       2012.07 created [Martin Haeufel]
!>       2016.07 Adapted for Harris functional {JGH]
!> \author Martin Haeufel
! **************************************************************************************************
   SUBROUTINE ec_build_neighborlist(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_build_neighborlist'

      INTEGER                                            :: handle, ikind, nkind, zat
      LOGICAL                                            :: gth_potential_present, &
                                                            sgp_potential_present, &
                                                            skip_load_balance_distributed
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: default_present, orb_present, &
                                                            ppl_present, ppnl_present
      REAL(dp)                                           :: subcells
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: c_radius, orb_radius, ppl_radius, &
                                                            ppnl_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_cn, sab_vdw
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind_set(qs_kind_set, gth_potential_present=gth_potential_present, &
                           sgp_potential_present=sgp_potential_present)
      nkind = SIZE(qs_kind_set)
      ALLOCATE (c_radius(nkind), default_present(nkind))
      ALLOCATE (orb_radius(nkind), ppl_radius(nkind), ppnl_radius(nkind))
      ALLOCATE (orb_present(nkind), ppl_present(nkind), ppnl_present(nkind))
      ALLOCATE (pair_radius(nkind, nkind))
      ALLOCATE (atom2d(nkind))

      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      cell=cell, &
                      distribution_2d=distribution_2d, &
                      local_particles=distribution_1d, &
                      particle_set=particle_set, &
                      molecule_set=molecule_set)

      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, .FALSE., particle_set)

      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom2d(ikind)%list)
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="HARRIS")
         IF (ASSOCIATED(basis_set)) THEN
            orb_present(ikind) = .TRUE.
            CALL get_gto_basis_set(gto_basis_set=basis_set, kind_radius=orb_radius(ikind))
         ELSE
            orb_present(ikind) = .FALSE.
            orb_radius(ikind) = 0.0_dp
         END IF
         CALL get_qs_kind(qs_kind, gth_potential=gth_potential, sgp_potential=sgp_potential)
         IF (gth_potential_present .OR. sgp_potential_present) THEN
            IF (ASSOCIATED(gth_potential)) THEN
               CALL get_potential(potential=gth_potential, &
                                  ppl_present=ppl_present(ikind), &
                                  ppl_radius=ppl_radius(ikind), &
                                  ppnl_present=ppnl_present(ikind), &
                                  ppnl_radius=ppnl_radius(ikind))
            ELSE IF (ASSOCIATED(sgp_potential)) THEN
               CALL get_potential(potential=sgp_potential, &
                                  ppl_present=ppl_present(ikind), &
                                  ppl_radius=ppl_radius(ikind), &
                                  ppnl_present=ppnl_present(ikind), &
                                  ppnl_radius=ppnl_radius(ikind))
            ELSE
               ppl_present(ikind) = .FALSE.
               ppl_radius(ikind) = 0.0_dp
               ppnl_present(ikind) = .FALSE.
               ppnl_radius(ikind) = 0.0_dp
            END IF
         END IF
      END DO

      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)

      ! overlap
      CALL pair_radius_setup(orb_present, orb_present, orb_radius, orb_radius, pair_radius)
      CALL build_neighbor_lists(ec_env%sab_orb, particle_set, atom2d, cell, pair_radius, &
                                subcells=subcells, nlname="sab_orb")
      ! pseudopotential
      IF (gth_potential_present .OR. sgp_potential_present) THEN
         IF (ANY(ppl_present)) THEN
            CALL pair_radius_setup(orb_present, ppl_present, orb_radius, ppl_radius, pair_radius)
            CALL build_neighbor_lists(ec_env%sac_ppl, particle_set, atom2d, cell, pair_radius, &
                                      subcells=subcells, operator_type="ABC", nlname="sac_ppl")
         END IF

         IF (ANY(ppnl_present)) THEN
            CALL pair_radius_setup(orb_present, ppnl_present, orb_radius, ppnl_radius, pair_radius)
            CALL build_neighbor_lists(ec_env%sap_ppnl, particle_set, atom2d, cell, pair_radius, &
                                      subcells=subcells, operator_type="ABBA", nlname="sap_ppnl")
         END IF
      END IF

      ! Build the neighbor lists for the vdW pair potential
      c_radius(:) = 0.0_dp
      dispersion_env => ec_env%dispersion_env
      sab_vdw => dispersion_env%sab_vdw
      sab_cn => dispersion_env%sab_cn
      IF (dispersion_env%type == xc_vdw_fun_pairpot) THEN
         c_radius(:) = dispersion_env%rc_disp
         default_present = .TRUE. !include all atoms in vdW (even without basis)
         CALL pair_radius_setup(default_present, default_present, c_radius, c_radius, pair_radius)
         CALL build_neighbor_lists(sab_vdw, particle_set, atom2d, cell, pair_radius, &
                                   subcells=subcells, operator_type="PP", nlname="sab_vdw")
         dispersion_env%sab_vdw => sab_vdw
         IF (dispersion_env%pp_type == vdw_pairpot_dftd3 .OR. &
             dispersion_env%pp_type == vdw_pairpot_dftd3bj) THEN
            ! Build the neighbor lists for coordination numbers as needed by the DFT-D3 method
            DO ikind = 1, nkind
               CALL get_atomic_kind(atomic_kind_set(ikind), z=zat)
               c_radius(ikind) = 4._dp*ptable(zat)%covalent_radius*bohr
            END DO
            CALL pair_radius_setup(default_present, default_present, c_radius, c_radius, pair_radius)
            CALL build_neighbor_lists(sab_cn, particle_set, atom2d, cell, pair_radius, &
                                      subcells=subcells, operator_type="PP", nlname="sab_cn")
            dispersion_env%sab_cn => sab_cn
         END IF
      END IF

      ! Release work storage
      CALL atom2d_cleanup(atom2d)
      DEALLOCATE (atom2d)
      DEALLOCATE (orb_present, default_present, ppl_present, ppnl_present)
      DEALLOCATE (orb_radius, ppl_radius, ppnl_radius, c_radius)
      DEALLOCATE (pair_radius)

      ! Task list
      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control)
      skip_load_balance_distributed = dft_control%qs_control%skip_load_balance_distributed
      IF (ASSOCIATED(ec_env%task_list)) CALL deallocate_task_list(ec_env%task_list)
      CALL allocate_task_list(ec_env%task_list)
      CALL generate_qs_task_list(ks_env, ec_env%task_list, &
                                 reorder_rs_grid_ranks=.FALSE., soft_valid=.FALSE., &
                                 skip_load_balance_distributed=skip_load_balance_distributed, &
                                 basis_type="HARRIS", sab_orb_external=ec_env%sab_orb)

      CALL timestop(handle)

   END SUBROUTINE ec_build_neighborlist

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ec_env ...
! **************************************************************************************************
   SUBROUTINE ec_properties(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ec_properties'

      CHARACTER(LEN=8), DIMENSION(3)                     :: rlab
      CHARACTER(LEN=default_path_length)                 :: filename, my_pos_voro
      CHARACTER(LEN=default_string_length)               :: description
      INTEGER :: akind, handle, i, ia, iatom, idir, ikind, iounit, ispin, maxmom, nspins, &
         reference, should_print_bqb, should_print_voro, unit_nr, unit_nr_voro
      LOGICAL                                            :: append_voro, magnetic, periodic, &
                                                            voro_print_txt
      REAL(KIND=dp)                                      :: charge, dd, focc, tmp
      REAL(KIND=dp), DIMENSION(3)                        :: cdip, pdip, rcc, rdip, ria, tdip
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, moments
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_elec_rspace
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type), POINTER                             :: mb_rho
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: ec_section, print_key, print_key_bqb, &
                                                            print_key_voro

      CALL timeset(routineN, handle)

      rlab(1) = "X"
      rlab(2) = "Y"
      rlab(3) = "Z"

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      END IF

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)

      ec_section => section_vals_get_subs_vals(qs_env%input, "DFT%ENERGY_CORRECTION")
      print_key => section_vals_get_subs_vals(section_vals=ec_section, &
                                              subsection_name="PRINT%MOMENTS")

      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN

         maxmom = section_get_ival(section_vals=ec_section, &
                                   keyword_name="PRINT%MOMENTS%MAX_MOMENT")
         periodic = section_get_lval(section_vals=ec_section, &
                                     keyword_name="PRINT%MOMENTS%PERIODIC")
         reference = section_get_ival(section_vals=ec_section, &
                                      keyword_name="PRINT%MOMENTS%REFERENCE")
         magnetic = section_get_lval(section_vals=ec_section, &
                                     keyword_name="PRINT%MOMENTS%MAGNETIC")
         NULLIFY (ref_point)
         CALL section_vals_val_get(ec_section, "PRINT%MOMENTS%REF_POINT", r_vals=ref_point)
         unit_nr = cp_print_key_unit_nr(logger=logger, basis_section=ec_section, &
                                        print_key_path="PRINT%MOMENTS", extension=".dat", &
                                        middle_name="moments", log_filename=.FALSE.)

         IF (iounit > 0) THEN
            IF (unit_nr /= iounit .AND. unit_nr > 0) THEN
               INQUIRE (UNIT=unit_nr, NAME=filename)
               WRITE (UNIT=iounit, FMT="(/,T2,A,2(/,T3,A),/)") &
                  "MOMENTS", "The electric/magnetic moments are written to file:", &
                  TRIM(filename)
            ELSE
               WRITE (UNIT=iounit, FMT="(/,T2,A)") "ELECTRIC/MAGNETIC MOMENTS"
            END IF
         END IF

         IF (periodic) THEN
            CPABORT("Periodic moments not implemented with EC")
         ELSE
            CPASSERT(maxmom < 2)
            CPASSERT(.NOT. magnetic)
            IF (maxmom == 1) THEN
               CALL get_qs_env(qs_env=qs_env, cell=cell, para_env=para_env)
               ! reference point
               CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)
               ! nuclear contribution
               cdip = 0.0_dp
               CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, &
                               qs_kind_set=qs_kind_set, local_particles=local_particles)
               DO ikind = 1, SIZE(local_particles%n_el)
                  DO ia = 1, local_particles%n_el(ikind)
                     iatom = local_particles%list(ikind)%array(ia)
                     ! fold atomic positions back into unit cell
                     ria = pbc(particle_set(iatom)%r - rcc, cell) + rcc
                     ria = ria - rcc
                     atomic_kind => particle_set(iatom)%atomic_kind
                     CALL get_atomic_kind(atomic_kind, kind_number=akind)
                     CALL get_qs_kind(qs_kind_set(akind), core_charge=charge)
                     cdip(1:3) = cdip(1:3) - charge*ria(1:3)
                  END DO
               END DO
               CALL mp_sum(cdip, para_env%group)
               !
               ! direct density contribution
               CALL ec_efield_integrals(qs_env, ec_env, rcc)
               !
               nspins = SIZE(ec_env%matrix_p, 1)
               pdip = 0.0_dp
               DO ispin = 1, nspins
                  DO idir = 1, 3
                     CALL dbcsr_dot(ec_env%matrix_p(ispin, 1)%matrix, &
                                    ec_env%efield%dipmat(idir)%matrix, tmp)
                     pdip(idir) = pdip(idir) + tmp
                  END DO
               END DO
               !
               ! response contribution
               CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
               NULLIFY (moments)
               CALL dbcsr_allocate_matrix_set(moments, 4)
               DO i = 1, 4
                  ALLOCATE (moments(i)%matrix)
                  CALL dbcsr_copy(moments(i)%matrix, matrix_s(1)%matrix, "Moments")
                  CALL dbcsr_set(moments(i)%matrix, 0.0_dp)
               END DO
               CALL build_local_moment_matrix(qs_env, moments, 1, ref_point=rcc)
               !
               focc = 2.0_dp
               IF (nspins == 2) focc = 1.0_dp
               rdip = 0.0_dp
               DO ispin = 1, nspins
                  DO idir = 1, 3
                     CALL dbcsr_dot(ec_env%matrix_z(ispin)%matrix, moments(idir)%matrix, tmp)
                     rdip(idir) = rdip(idir) + tmp
                  END DO
               END DO
               CALL dbcsr_deallocate_matrix_set(moments)
               !
               tdip = -(rdip + pdip + cdip)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, "(T3,A)") "Dipoles are based on the traditional operator."
                  dd = SQRT(SUM(tdip(1:3)**2))*debye
                  WRITE (unit_nr, "(T3,A)") "Dipole moment [Debye]"
                  WRITE (unit_nr, "(T5,3(A,A,F14.8,1X),T60,A,T67,F14.8)") &
                     (TRIM(rlab(i)), "=", tdip(i)*debye, i=1, 3), "Total=", dd
               END IF
            END IF
         END IF

         CALL cp_print_key_finished_output(unit_nr=unit_nr, logger=logger, &
                                           basis_section=ec_section, print_key_path="PRINT%MOMENTS")
         CALL get_qs_env(qs_env=qs_env, results=results)
         description = "[DIPOLE]"
         CALL cp_results_erase(results=results, description=description)
         CALL put_results(results=results, description=description, values=tdip(1:3))
      END IF

      ! Do a Voronoi Integration or write a compressed BQB File
      print_key_voro => section_vals_get_subs_vals(ec_section, "PRINT%VORONOI")
      print_key_bqb => section_vals_get_subs_vals(ec_section, "PRINT%E_DENSITY_BQB")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key_voro), cp_p_file)) THEN
         should_print_voro = 1
      ELSE
         should_print_voro = 0
      END IF
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key_bqb), cp_p_file)) THEN
         should_print_bqb = 1
      ELSE
         should_print_bqb = 0
      END IF
      IF ((should_print_voro /= 0) .OR. (should_print_bqb /= 0)) THEN

         CALL get_qs_env(qs_env=qs_env, &
                         pw_env=pw_env)
         CALL pw_env_get(pw_env=pw_env, &
                         auxbas_pw_pool=auxbas_pw_pool, &
                         pw_pools=pw_pools)
         CALL pw_pool_create_pw(pool=auxbas_pw_pool, &
                                pw=rho_elec_rspace%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)

         IF (dft_control%nspins > 1) THEN

            ! add Pout and Pz
            CALL pw_copy(ec_env%rhoout_r(1)%pw, rho_elec_rspace%pw)
            CALL pw_axpy(ec_env%rhoout_r(2)%pw, rho_elec_rspace%pw)

            CALL pw_axpy(ec_env%rhoz_r(1)%pw, rho_elec_rspace%pw)
            CALL pw_axpy(ec_env%rhoz_r(2)%pw, rho_elec_rspace%pw)

            mb_rho => rho_elec_rspace%pw
         ELSE

            ! add Pout and Pz
            CALL pw_copy(ec_env%rhoout_r(1)%pw, rho_elec_rspace%pw)
            CALL pw_axpy(ec_env%rhoz_r(1)%pw, rho_elec_rspace%pw)

            mb_rho => rho_elec_rspace%pw
         END IF ! nspins

         IF (should_print_voro /= 0) THEN
            CALL section_vals_val_get(print_key_voro, "OUTPUT_TEXT", l_val=voro_print_txt)
            IF (voro_print_txt) THEN
               append_voro = section_get_lval(ec_section, "PRINT%VORONOI%APPEND")
               my_pos_voro = "REWIND"
               IF (append_voro) THEN
                  my_pos_voro = "APPEND"
               END IF
               unit_nr_voro = cp_print_key_unit_nr(logger, ec_section, "PRINT%VORONOI", extension=".voronoi", &
                                                   file_position=my_pos_voro, log_filename=.FALSE.)
            ELSE
               unit_nr_voro = 0
            END IF
         ELSE
            unit_nr_voro = 0
         END IF

         CALL entry_voronoi_or_bqb(should_print_voro, should_print_bqb, print_key_voro, print_key_bqb, &
                                   unit_nr_voro, qs_env, mb_rho)

         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_elec_rspace%pw)

         IF (unit_nr_voro > 0) THEN
            CALL cp_print_key_finished_output(unit_nr_voro, logger, ec_section, "PRINT%VORONOI")
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_properties

! **************************************************************************************************
!> \brief Solve the Harris functional by linear scaling density purification scheme,
!>        instead of the diagonalization performed in ec_diag_solver
!>
!> \param qs_env ...
!> \param matrix_ks Harris Kohn-Sham matrix
!> \param matrix_s Overlap matrix in Harris functional basis
!> \par History
!>       09.2020 created
!> \author F.Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_ls_init(qs_env, matrix_ks, matrix_s)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_ls_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(ls_scf_env_type), POINTER                     :: ls_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      ec_env=ec_env)
      nspins = dft_control%nspins
      ls_env => ec_env%ls_env

      ! create the matrix template for use in the ls procedures
      CALL matrix_ls_create(matrix_ls=ls_env%matrix_s, matrix_qs=matrix_s(1, 1)%matrix, &
                            ls_mstruct=ls_env%ls_mstruct)

      IF (ALLOCATED(ls_env%matrix_p)) THEN
         DO ispin = 1, SIZE(ls_env%matrix_p)
            CALL dbcsr_release(ls_env%matrix_p(ispin))
         END DO
      ELSE
         ALLOCATE (ls_env%matrix_p(nspins))
      END IF

      DO ispin = 1, nspins
         CALL dbcsr_create(ls_env%matrix_p(ispin), template=ls_env%matrix_s, &
                           matrix_type=dbcsr_type_no_symmetry)
      END DO

      ALLOCATE (ls_env%matrix_ks(nspins))
      DO ispin = 1, nspins
         CALL dbcsr_create(ls_env%matrix_ks(ispin), template=ls_env%matrix_s, &
                           matrix_type=dbcsr_type_no_symmetry)
      END DO

      ! Set up S matrix and needed functions of S
      CALL ls_scf_init_matrix_s(matrix_s(1, 1)%matrix, ls_env)

      ! Bring KS matrix from QS to LS form
      ! EC KS-matrix already calculated
      DO ispin = 1, nspins
         CALL matrix_qs_to_ls(matrix_ls=ls_env%matrix_ks(ispin), &
                              matrix_qs=matrix_ks(ispin, 1)%matrix, &
                              ls_mstruct=ls_env%ls_mstruct, &
                              covariant=.TRUE.)
      END DO

      CALL timestop(handle)

   END SUBROUTINE ec_ls_init

! **************************************************************************************************
!> \brief Solve the Harris functional by linear scaling density purification scheme,
!>        instead of the diagonalization performed in ec_diag_solver
!>
!> \param qs_env ...
!> \param matrix_p Harris dentiy matrix, calculated here
!> \param matrix_w Harris energy weighted density matrix, calculated here
!> \param ec_ls_method which purification scheme should be used
!> \par History
!>      12.2019 created [JGH]
!>      08.2020 refactoring [fbelle]
!> \author Fabian Belleflamme
! **************************************************************************************************

   SUBROUTINE ec_ls_solver(qs_env, matrix_p, matrix_w, ec_ls_method)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_w
      INTEGER, INTENT(IN)                                :: ec_ls_method

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_ls_solver', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nelectron_spin_real, &
                                                            nspins
      INTEGER, DIMENSION(2)                              :: nmo
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: wmat
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: matrix_ks_deviation
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(ls_scf_env_type), POINTER                     :: ls_env

      CALL timeset(routineN, handle)

      NULLIFY (para_env)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env)
      nspins = dft_control%nspins
      ec_env => qs_env%ec_env
      ls_env => ec_env%ls_env

      nmo = 0
      CALL get_qs_env(qs_env=qs_env, nelectron_spin=nmo)
      IF (nspins == 1) nmo(1) = nmo(1)/2
      ls_env%homo_spin(:) = 0.0_dp
      ls_env%lumo_spin(:) = 0.0_dp

      ALLOCATE (matrix_ks_deviation(nspins))
      DO ispin = 1, nspins
         CALL dbcsr_create(matrix_ks_deviation(ispin), template=ls_env%matrix_ks(ispin))
         CALL dbcsr_set(matrix_ks_deviation(ispin), 0.0_dp)
      END DO

      ! F = S^-1/2 * F * S^-1/2
      IF (ls_env%has_s_preconditioner) THEN
         DO ispin = 1, nspins
            CALL apply_matrix_preconditioner(ls_env%matrix_ks(ispin), "forward", &
                                             ls_env%matrix_bs_sqrt, ls_env%matrix_bs_sqrt_inv)

            CALL dbcsr_filter(ls_env%matrix_ks(ispin), ls_env%eps_filter)
         END DO
      END IF

      DO ispin = 1, nspins
         nelectron_spin_real = ls_env%nelectron_spin(ispin)
         IF (ls_env%nspins == 1) nelectron_spin_real = nelectron_spin_real/2

         SELECT CASE (ec_ls_method)
         CASE (ec_matrix_sign)
            CALL density_matrix_sign(ls_env%matrix_p(ispin), &
                                     ls_env%mu_spin(ispin), &
                                     ls_env%fixed_mu, &
                                     ls_env%sign_method, &
                                     ls_env%sign_order, &
                                     ls_env%matrix_ks(ispin), &
                                     ls_env%matrix_s, &
                                     ls_env%matrix_s_inv, &
                                     nelectron_spin_real, &
                                     ec_env%eps_default)

         CASE (ec_matrix_trs4)
            CALL density_matrix_trs4( &
               ls_env%matrix_p(ispin), &
               ls_env%matrix_ks(ispin), &
               ls_env%matrix_s_sqrt_inv, &
               nelectron_spin_real, &
               ec_env%eps_default, &
               ls_env%homo_spin(ispin), &
               ls_env%lumo_spin(ispin), &
               ls_env%mu_spin(ispin), &
               matrix_ks_deviation=matrix_ks_deviation(ispin), &
               dynamic_threshold=ls_env%dynamic_threshold, &
               eps_lanczos=ls_env%eps_lanczos, &
               max_iter_lanczos=ls_env%max_iter_lanczos)

         CASE (ec_matrix_tc2)
            CALL density_matrix_tc2( &
               ls_env%matrix_p(ispin), &
               ls_env%matrix_ks(ispin), &
               ls_env%matrix_s_sqrt_inv, &
               nelectron_spin_real, &
               ec_env%eps_default, &
               ls_env%homo_spin(ispin), &
               ls_env%lumo_spin(ispin), &
               non_monotonic=ls_env%non_monotonic, &
               eps_lanczos=ls_env%eps_lanczos, &
               max_iter_lanczos=ls_env%max_iter_lanczos)

         END SELECT

      END DO

      ! de-orthonormalize
      IF (ls_env%has_s_preconditioner) THEN
         DO ispin = 1, nspins
            ! P = S^-1/2 * P_tilde * S^-1/2 (forward)
            CALL apply_matrix_preconditioner(ls_env%matrix_p(ispin), "forward", &
                                             ls_env%matrix_bs_sqrt, ls_env%matrix_bs_sqrt_inv)

            CALL dbcsr_filter(ls_env%matrix_p(ispin), ls_env%eps_filter)
         END DO
      END IF

      ! Closed-shell
      IF (nspins == 1) CALL dbcsr_scale(ls_env%matrix_p(1), 2.0_dp)

      IF (ls_env%report_all_sparsities) CALL post_scf_sparsities(ls_env)

      ! ls_scf_dm_to_ks
      ! Density matrix from LS to EC
      DO ispin = 1, nspins
         CALL matrix_ls_to_qs(matrix_qs=matrix_p(ispin, 1)%matrix, &
                              matrix_ls=ls_env%matrix_p(ispin), &
                              ls_mstruct=ls_env%ls_mstruct, &
                              covariant=.FALSE.)
      END DO

      wmat => matrix_w(:, 1)
      CALL calculate_w_matrix_ls(wmat, ec_env%ls_env)

      ! clean up
      CALL dbcsr_release(ls_env%matrix_s)
      IF (ls_env%has_s_preconditioner) THEN
         CALL dbcsr_release(ls_env%matrix_bs_sqrt)
         CALL dbcsr_release(ls_env%matrix_bs_sqrt_inv)
      END IF
      IF (ls_env%needs_s_inv) THEN
         CALL dbcsr_release(ls_env%matrix_s_inv)
      END IF
      IF (ls_env%use_s_sqrt) THEN
         CALL dbcsr_release(ls_env%matrix_s_sqrt)
         CALL dbcsr_release(ls_env%matrix_s_sqrt_inv)
      END IF

      DO ispin = 1, SIZE(ls_env%matrix_ks)
         CALL dbcsr_release(ls_env%matrix_ks(ispin))
      END DO
      DEALLOCATE (ls_env%matrix_ks)

      DO ispin = 1, nspins
         CALL dbcsr_release(matrix_ks_deviation(ispin))
      END DO
      DEALLOCATE (matrix_ks_deviation)

      CALL timestop(handle)

   END SUBROUTINE ec_ls_solver

! **************************************************************************************************
!> \brief Use OT-diagonalziation to obtain density matrix from Harris Kohn-Sham matrix
!>        Initial guess of density matrix is either the atomic block initial guess from SCF
!>        or the ground-state density matrix. The latter only works if the same basis is used
!>
!> \param qs_env ...
!> \param ec_env ...
!> \param matrix_ks Harris Kohn-Sham matrix
!> \param matrix_s Overlap matrix in Harris functional basis
!> \param matrix_p Harris dentiy matrix, calculated here
!> \param matrix_w Harris energy weighted density matrix, calculated here
!>
!> \par History
!>       09.2020 created
!> \author F.Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_ot_diag_solver(qs_env, ec_env, matrix_ks, matrix_s, matrix_p, matrix_w)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN), &
         POINTER                                         :: matrix_ks, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), &
         INTENT(INOUT), POINTER                          :: matrix_p, matrix_w

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_ot_diag_solver', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, homo, ikind, iounit, ispin, &
                                                            max_iter, nao, nkind, nmo, nspins
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, sv
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_rmpv
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set, harris_basis
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(preconditioner_type), POINTER                 :: local_preconditioner
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_unit_nr(logger)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(ec_env))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_s))
      CPASSERT(ASSOCIATED(matrix_p))
      CPASSERT(ASSOCIATED(matrix_w))

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      blacs_env=blacs_env, &
                      dft_control=dft_control, &
                      nelectron_spin=nelectron_spin, &
                      para_env=para_env, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)
      nspins = dft_control%nspins

      ! Maximum number of OT iterations for diagonalization
      max_iter = 200

      ! If linear scaling, need to allocate and init MO set
      ! set it to qs_env%mos
      IF (dft_control%qs_control%do_ls_scf) THEN
         CALL ec_mos_init(qs_env, matrix_s(1, 1)%matrix)
      END IF

      CALL get_qs_env(qs_env, mos=mos)

      ! Inital guess to use
      NULLIFY (p_rmpv)

      ! Using ether ground-state DM or ATOMIC GUESS requires
      ! Harris functional to use the same basis set
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, nkind=nkind)
      CALL uppercase(ec_env%basis)
      ! Harris basis may differ from ground-state basis only if explicitly added
      ! thus only two cases that need to be tested
      ! 1) explicit Harris basis present?
      IF (ec_env%basis == "HARRIS") THEN
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            ! Basis sets of ground-state
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="ORB")
            ! Basis sets of energy correction
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=harris_basis, basis_type="HARRIS")

            IF (basis_set%name .NE. harris_basis%name) THEN
               CPABORT("OT-Diag initial guess: Harris and ground state need to use the same basis")
            END IF
         END DO
      END IF
      ! 2) Harris uses MAOs
      IF (ec_env%mao) THEN
         CPABORT("OT-Diag initial guess: not implemented for MAOs")
      END IF

      ! Initital guess obtained for OT Diagonalization
      SELECT CASE (ec_env%ec_initial_guess)
      CASE (ec_ot_atomic)

         p_rmpv => matrix_p(:, 1)

         CALL calculate_atomic_block_dm(p_rmpv, matrix_s(1, 1)%matrix, &
                                        particle_set, atomic_kind_set, qs_kind_set, &
                                        nspins, nelectron_spin, iounit, para_env)

      CASE (ec_ot_gs)

         CALL get_qs_env(qs_env, rho=rho)
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao)
         p_rmpv => rho_ao(:, 1)

      CASE DEFAULT
         CPABORT("Unknown inital guess for OT-Diagonalization (Harris functional)")
      END SELECT

      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                         mo_coeff=mo_coeff, &
                         nmo=nmo, &
                         nao=nao, &
                         homo=homo)

         ! Calculate first MOs
         CALL cp_fm_set_all(mo_coeff, 0.0_dp)
         CALL cp_fm_init_random(mo_coeff, nmo)

         NULLIFY (sv)
         CALL cp_fm_create(sv, mo_coeff%matrix_struct, "SV")
         ! multiply times PS
         ! PS*C(:,1:nomo)+C(:,nomo+1:nmo) (nomo=NINT(nelectron/maxocc))
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1, 1)%matrix, mo_coeff, sv, nmo)
         CALL cp_dbcsr_sm_fm_multiply(p_rmpv(ispin)%matrix, sv, mo_coeff, homo)
         CALL cp_fm_release(sv)
         ! and ortho the result
         ! If DFBT or SE, then needs has_unit_metrix option
         CALL make_basis_sm(mo_coeff, nmo, matrix_s(1, 1)%matrix)
      END DO

      ! Preconditioner
      NULLIFY (local_preconditioner)
      ALLOCATE (local_preconditioner)
      CALL init_preconditioner(local_preconditioner, para_env=para_env, &
                               blacs_env=blacs_env)
      DO ispin = 1, nspins
         CALL make_preconditioner(local_preconditioner, &
                                  precon_type=ot_precond_full_single_inverse, &
                                  solver_type=ot_precond_solver_default, &
                                  matrix_h=matrix_ks(ispin, 1)%matrix, &
                                  matrix_s=matrix_s(ispin, 1)%matrix, &
                                  convert_precond_to_dbcsr=.TRUE., &
                                  mo_set=mos(ispin)%mo_set, energy_gap=0.2_dp)

         CALL get_mo_set(mos(ispin)%mo_set, &
                         mo_coeff=mo_coeff, &
                         eigenvalues=eigenvalues, &
                         nmo=nmo, &
                         homo=homo)
         CALL ot_eigensolver(matrix_h=matrix_ks(ispin, 1)%matrix, &
                             matrix_s=matrix_s(1, 1)%matrix, &
                             matrix_c_fm=mo_coeff, &
                             preconditioner=local_preconditioner, &
                             eps_gradient=ec_env%eps_default, &
                             iter_max=max_iter, &
                             silent=.FALSE.)
         CALL calculate_subspace_eigenvalues(mo_coeff, matrix_ks(ispin, 1)%matrix, &
                                             evals_arg=eigenvalues, do_rotation=.TRUE.)

         ! Deallocate preconditioner
         CALL destroy_preconditioner(local_preconditioner)
         DEALLOCATE (local_preconditioner)

         !fm->dbcsr
         CALL copy_fm_to_dbcsr(mos(ispin)%mo_set%mo_coeff, &
                               mos(ispin)%mo_set%mo_coeff_b)
      END DO

      ! Calculate density matrix from MOs
      DO ispin = 1, nspins
         CALL calculate_density_matrix(mos(ispin)%mo_set, matrix_p(ispin, 1)%matrix)

         CALL calculate_w_matrix(mos(ispin)%mo_set, matrix_w(ispin, 1)%matrix)
      END DO

      ! Get rid of MO environment again
      IF (dft_control%qs_control%do_ls_scf) THEN
         DO ispin = 1, nspins
            CALL deallocate_mo_set(mos(ispin)%mo_set)
         END DO
         IF (ASSOCIATED(qs_env%mos)) THEN
            DO ispin = 1, SIZE(qs_env%mos)
               CALL deallocate_mo_set(qs_env%mos(ispin)%mo_set)
            END DO
            DEALLOCATE (qs_env%mos)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_ot_diag_solver

END MODULE energy_corrections

