!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for total energy and forces of excited states
!> \par History
!>       01.2020 created
!> \author JGH
! **************************************************************************************************
MODULE excited_states
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
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
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
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
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, dbcsr_distribution_type, &
        dbcsr_dot, dbcsr_filter, dbcsr_get_info, dbcsr_init_p, dbcsr_multiply, dbcsr_p_type, &
        dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE dm_ls_scf_methods,               ONLY: density_matrix_sign,&
                                              density_matrix_tc2,&
                                              density_matrix_trs4
   USE exstates_types,                  ONLY: excited_energy_type
   USE ec_efield_local,                 ONLY: ec_efield_integrals,&
                                              ec_efield_local_operator
   USE ec_env_types,                    ONLY: energy_correction_type
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE input_constants,                 ONLY: ec_curvy_steps,&
                                              ec_diagonalization,&
                                              ec_functional_harris,&
                                              ec_matrix_sign,&
                                              ec_matrix_tc2,&
                                              ec_matrix_trs4,&
                                              ls_s_sqrt_ns,&
                                              ls_s_sqrt_proot
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_lval,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE iterate_matrix,                  ONLY: invert_Hotelling,&
                                              matrix_sqrt_Newton_Schulz,&
                                              matrix_sqrt_proot
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE mao_basis,                       ONLY: mao_generate_basis
   USE message_passing,                 ONLY: mp_sum
   USE molecule_types,                  ONLY: molecule_type
   USE moments_utils,                   ONLY: get_reference_point
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: debye
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integral_ab,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_core_energies,                ONLY: calculate_ecore_overlap
   USE qs_dispersion_pairpot,           ONLY: calculate_dispersion_pairpot
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: allocate_qs_force,&
                                              qs_force_type,&
                                              sum_qs_force,&
                                              total_qs_force,&
                                              zero_qs_force
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE response_solver,                 ONLY: ks_ref_potential,&
                                              response_equation,&
                                              response_force
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE task_list_types,                 ONLY: allocate_task_list,&
                                              deallocate_task_list
   USE virial_types,                    ONLY: virial_type
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'excited_states'

   LOGICAL, PARAMETER                   :: debug_forces = .TRUE.

   PUBLIC :: excited_state_energy

CONTAINS

! **************************************************************************************************
!> \brief Excited state energy and forces
!>
!> \param qs_env ...
!> \param calculate_forces ...
!> \par History
!>       03.2014 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE excited_state_energy(qs_env, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'excited_state_energy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(excited_energy_type), POINTER                 :: ex_env

      CALL timeset(routineN, handle)

      ! Check for energy correction
      IF (qs_env%excited_state) THEN
         logger => cp_get_default_logger()
         IF (logger%para_env%ionode) THEN
            unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
         ELSE
            unit_nr = -1
         ENDIF

         CALL get_qs_env(qs_env, exstate_env=ex_env, energy=energy)

         energy%total = energy%total + ex_env%evalue
         IF (calculate_forces) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,A,A,A,A)') "!", REPEAT("-", 27), &
                  " Excited State Forces ", REPEAT("-", 28), "!"
            END IF
         ELSE
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,A,A,A,A)') "!", REPEAT("-", 27), &
                  " Excited State Energy ", REPEAT("-", 28), "!"
               WRITE (unit_nr, '(T3,A,T65,F16.10)') "Excitation Energy [Hartree] ", ex_env%evalue
               WRITE (unit_nr, '(T3,A,T65,F16.10)') "Total Energy [Hartree]", energy%total
            END IF
         END IF

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A,A,A)') "!", REPEAT("-", 77), "!"
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE excited_state_energy

! **************************************************************************************************

END MODULE excited_states
