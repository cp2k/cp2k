!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'energy_corrections'

   LOGICAL, PARAMETER                   :: debug_forces = .TRUE.

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

      INTEGER                                            :: handle, nkind, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: natom_of_kind
      LOGICAL                                            :: my_calc_forces
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: ec_force, ks_force

      CALL timeset(routineN, handle)

      ! Check for energy correction
      IF (qs_env%energy_correction) THEN
         NULLIFY (ec_env)
         CALL get_qs_env(qs_env, ec_env=ec_env)

         ec_env%should_update = .TRUE.
         IF (PRESENT(ec_init)) ec_env%should_update = ec_init

         my_calc_forces = .FALSE.
         IF (PRESENT(calculate_forces)) my_calc_forces = calculate_forces

         logger => cp_get_default_logger()
         IF (logger%para_env%ionode) THEN
            unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
         ELSE
            unit_nr = -1
         ENDIF

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
            CALL get_qs_env(qs_env, force=ks_force, atomic_kind_set=atomic_kind_set)
            nkind = SIZE(atomic_kind_set)
            ALLOCATE (natom_of_kind(nkind))
            CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom_of_kind=natom_of_kind)
            NULLIFY (ec_force)
            CALL allocate_qs_force(ec_force, natom_of_kind)
            DEALLOCATE (natom_of_kind)
            CALL zero_qs_force(ec_force)
            CALL set_qs_env(qs_env, force=ec_force)
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
            END IF
            IF (my_calc_forces) THEN
               CALL get_qs_env(qs_env, force=ec_force)
               ec_env%force => ec_force
               CALL zero_qs_force(ks_force)
               CALL sum_qs_force(ks_force, ec_force)
               CALL set_qs_env(qs_env, force=ks_force)
            END IF
            energy%total = energy%total + energy%nonscf_correction
         CASE DEFAULT
            CPABORT("unknown energy correction")
         END SELECT

         IF (.NOT. my_calc_forces .AND. unit_nr > 0) THEN
            WRITE (unit_nr, '(T3,A,T65,F16.10)') "Energy Correction ", energy%nonscf_correction
         END IF
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A,A,A)') "!", REPEAT("-", 77), "!"
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE energy_correction

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

      REAL(KIND=dp)                                      :: exc

      IF (ec_env%should_update) THEN
         CALL ec_build_neighborlist(qs_env, ec_env)
         !
         CALL ec_build_core_hamiltonian(qs_env, ec_env)
         CALL ks_ref_potential(qs_env, ec_env%vh_rspace, ec_env%vxc_rspace, ec_env%vtau_rspace, &
                               ec_env%ehartree, exc)
         CALL ec_build_ks_matrix(qs_env, ec_env)
         !
         IF (ec_env%mao) THEN
            ! MAO basis
            IF (ASSOCIATED(ec_env%mao_coef)) CALL dbcsr_deallocate_matrix_set(ec_env%mao_coef)
            NULLIFY (ec_env%mao_coef)
            CALL mao_generate_basis(qs_env, ec_env%mao_coef, ref_basis_set="HARRIS", molecular=.TRUE., &
                                    max_iter=ec_env%mao_max_iter, eps_grad=ec_env%mao_eps_grad, unit_nr=unit_nr)
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
         !
         CALL response_equation(qs_env, ec_env%p_env, ec_env%cpmos)
         !
         CALL response_force(qs_env, ec_env%p_env, &
                             ec_env%vh_rspace, ec_env%vxc_rspace, ec_env%vtau_rspace, &
                             ec_env%matrix_hz)
         !
         CALL ec_properties(qs_env, ec_env)
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
      iounit = cp_logger_get_default_unit_nr(logger)

      ! no k-points possible
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, dft_control=dft_control)
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

      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env)
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
                              qs_kind_set, atomic_kind_set, particle_set, sab_orb, sap_ppnl, eps_ppnl, &
                              nimages, cell_to_index, "HARRIS")
      END IF

      ! External field (nonperiodic case)
      ec_env%efield_nuclear = 0.0_dp
      CALL ec_efield_local_operator(qs_env, ec_env, calculate_forces)

      CALL timestop(handle)

   END SUBROUTINE ec_build_core_hamiltonian

! **************************************************************************************************
!> \brief Solve KS equation for a given matrix
!> \brief calculate the complete KS matrix
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
      iounit = cp_logger_get_default_unit_nr(logger)

      calculate_forces = .FALSE.

      ! get all information on the electronic density
      NULLIFY (rho, ks_env)
      CALL get_qs_env(qs_env=qs_env, rho=rho, virial=virial, dft_control=dft_control, &
                      para_env=para_env, blacs_env=blacs_env, ks_env=ks_env)

      nspins = dft_control%nspins
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      CPASSERT(.NOT. use_virial)

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
      ENDDO

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
         ENDDO
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
         evhxc = evhxc + pw_integral_ab(rho_r(ispin)%pw, v_rspace(ispin)%pw)/v_rspace(1)%pw%pw_grid%dvol
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
      ENDDO
      DEALLOCATE (v_rspace)

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
      REAL(KIND=dp)                                      :: eps_ppnl
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
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
      iounit = cp_logger_get_default_unit_nr(logger)

      calculate_forces = .TRUE.

      ! no k-points possible
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, dft_control=dft_control)
      nimages = dft_control%nimages
      IF (nimages /= 1) THEN
         CPABORT("K-points for Harris functional not implemented")
      END IF

      ! check for virial (currently no stress tensor available)
      CALL get_qs_env(qs_env=qs_env, virial=virial)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      IF (use_virial) THEN
         CPABORT("Stress Tensor for Harris functional not implemented")
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
      CALL build_overlap_matrix(ks_env, matrixkp_s=scrm, &
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
      IF (SIZE(ec_env%matrix_p, 1) == 2) THEN
         CALL dbcsr_add(ec_env%matrix_p(1, 1)%matrix, ec_env%matrix_p(2, 1)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      END IF

      ! compute the ppl contribution to the core hamiltonian
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)
      NULLIFY (cell_to_index, virial)
      use_virial = .FALSE.
      IF (ASSOCIATED(sac_ppl)) THEN
         IF (calculate_forces .AND. debug_forces) fodeb(1:3) = force(1)%gth_ppl(1:3, 1)
         CALL build_core_ppl(scrm, ec_env%matrix_p, force, &
                             virial, calculate_forces, use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ppl, &
                             nimages, cell_to_index, "HARRIS")
         IF (calculate_forces .AND. debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppl(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pout*dH_PPL ", fodeb
         END IF
      END IF

      ! compute the ppnl contribution to the core hamiltonian ***
      eps_ppnl = dft_control%qs_control%eps_ppnl
      IF (ASSOCIATED(sap_ppnl)) THEN
         IF (calculate_forces .AND. debug_forces) fodeb(1:3) = force(1)%gth_ppnl(1:3, 1)
         CALL build_core_ppnl(scrm, ec_env%matrix_p, force, &
                              virial, calculate_forces, use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, sab_orb, sap_ppnl, eps_ppnl, &
                              nimages, cell_to_index, "HARRIS")
         IF (calculate_forces .AND. debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppnl(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pout*dH_PPNL", fodeb
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

      INTEGER                                            :: handle, i, iounit, ispin, nao, natom, &
                                                            norb, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: lsd, use_virial
      REAL(dp)                                           :: dehartree, eexc, eovrl, focc, &
                                                            total_rhoout
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: ftot
      REAL(dp), DIMENSION(3)                             :: fodeb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: vmos
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s, matrix_v, scrm
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: dv_hartree_rspace, rho_tot_gspace, &
                                                            v_hartree_gspace, v_hartree_rspace, &
                                                            vtot_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r, rhoout_g, rhoout_r, tau_r, &
                                                            v_rspace, v_tau_rspace, v_xc
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input, xc_fun_section, xc_section
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho1_set, rho_set

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_unit_nr(logger)

      ! get all information on the electronic density
      NULLIFY (rho, ks_env)
      CALL get_qs_env(qs_env=qs_env, rho=rho, virial=virial, dft_control=dft_control, &
                      para_env=para_env, blacs_env=blacs_env, ks_env=ks_env)

      nspins = dft_control%nspins
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      CPASSERT(.NOT. use_virial)

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CPASSERT(ASSOCIATED(pw_env))

      NULLIFY (auxbas_pw_pool, poisson_env)
      ! gets the tmp grids
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      ! Calculate the Hartree potential
      NULLIFY (v_hartree_gspace%pw, rho_tot_gspace%pw, v_hartree_rspace%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      ! Get the total density in g-space [ions + electrons]
      CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

      CALL pw_transfer(ec_env%vh_rspace%pw, v_hartree_rspace%pw)

      CALL get_qs_env(qs_env=qs_env, force=force)
      ! calculate output density on grid
      ! rho_in(R):   CALL qs_rho_get(rho, rho_r=rho_r)
      ! rho_in(G):   CALL qs_rho_get(rho, rho_g=rho_g)
      CALL qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g, tau_r=tau_r)
      ALLOCATE (rhoout_r(nspins), rhoout_g(nspins))
      DO ispin = 1, nspins
         NULLIFY (rhoout_r(ispin)%pw, rhoout_g(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoout_r(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoout_g(ispin)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      ENDDO
      CALL pw_pool_create_pw(auxbas_pw_pool, dv_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, vtot_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL pw_zero(rho_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL calculate_rho_elec(ks_env=ks_env, matrix_p=ec_env%matrix_p(ispin, 1)%matrix, &
                                 rho=rhoout_r(ispin), &
                                 rho_gspace=rhoout_g(ispin), &
                                 total_rho=total_rhoout, &
                                 basis_type="HARRIS", &
                                 task_list_external=ec_env%task_list)
         ! rho_out - rho_in
         CALL pw_axpy(rho_r(ispin)%pw, rhoout_r(ispin)%pw, -1.0_dp)
         CALL pw_axpy(rho_g(ispin)%pw, rhoout_g(ispin)%pw, -1.0_dp)
         CALL pw_axpy(rhoout_g(ispin)%pw, rho_tot_gspace%pw)
      END DO
      ! calculate associated hartree potential
      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, dehartree, &
                            v_hartree_gspace%pw)
      CALL pw_transfer(v_hartree_gspace%pw, dv_hartree_rspace%pw)
      CALL pw_scale(dv_hartree_rspace%pw, dv_hartree_rspace%pw%pw_grid%dvol)
      ! Getting nuclear force contribution from the core charge density
      ! Vh(rho_in + rho_c) + Vh(rho_out - rho_in)
      CALL pw_transfer(v_hartree_rspace%pw, vtot_rspace%pw)
      CALL pw_axpy(dv_hartree_rspace%pw, vtot_rspace%pw)
      IF (debug_forces) fodeb(1:3) = force(1)%rho_core(1:3, 1)
      CALL integrate_v_core_rspace(vtot_rspace, qs_env)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_core(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Vtot*dncore", fodeb
      END IF
      !
      ! Pulay force from Tr P_in (V_H(drho)+ Fxc(rho_in)*drho)
      ! RHS of CPKS equations: (V_H(drho)+ Fxc(rho_in)*drho)*C0
      ! Fxc*drho term
      ALLOCATE (v_xc(nspins))
      DO ispin = 1, nspins
         NULLIFY (v_xc(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(v_xc(ispin)%pw)
      END DO
      CALL get_qs_env(qs_env, input=input)
      xc_section => ec_env%xc_section
      lsd = (nspins == 2)
      NULLIFY (deriv_set, rho_set, rho1_set)
      CALL xc_prep_2nd_deriv(deriv_set, rho_set, rho_r, auxbas_pw_pool, xc_section=xc_section)
      bo = rhoout_r(1)%pw%pw_grid%bounds_local
      CALL xc_rho_set_create(rho1_set, bo, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

      ! calculate the arguments needed by the functionals
      CALL xc_rho_set_update(rho1_set, rhoout_r, rhoout_g, tau_r, needs, &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                             auxbas_pw_pool)
      CALL xc_calc_2nd_deriv(v_xc, deriv_set, rho_set, &
                             rho1_set, auxbas_pw_pool, xc_section=xc_section)
      CALL xc_dset_release(deriv_set)
      CALL xc_rho_set_release(rho_set)
      CALL xc_rho_set_release(rho1_set)
      !
      CALL get_qs_env(qs_env=qs_env, rho=rho, matrix_s_kp=matrix_s)
      NULLIFY (matrix_v)
      CALL dbcsr_allocate_matrix_set(matrix_v, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (matrix_v(ispin, 1)%matrix)
         CALL dbcsr_create(matrix_v(ispin, 1)%matrix, template=matrix_s(1, 1)%matrix)
         CALL dbcsr_copy(matrix_v(ispin, 1)%matrix, matrix_s(1, 1)%matrix)
         CALL dbcsr_set(matrix_v(ispin, 1)%matrix, 0.0_dp)
      END DO
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ! vtot = v_xc(ispin) + dv_hartree
      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      DO ispin = 1, nspins
         CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(dv_hartree_rspace%pw, v_xc(ispin)%pw)
         CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                 hmat=matrix_v(ispin, 1), &
                                 pmat=matrix_p(ispin, 1), &
                                 calculate_forces=.TRUE.)
      END DO
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dKdrho", fodeb
      END IF
      ! CPKS vector
      !
      ! sign is changed in linres_solver!
      ! Projector Q applied in linres_solver!
      focc = 2.0_dp
      IF (nspins == 1) focc = 4.0_dp
      CALL get_qs_env(qs_env=qs_env, mos=mos)
      ALLOCATE (vmos(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, matrix_struct=fm_struct, nrow_global=nao, ncol_global=norb)
         NULLIFY (vmos(ispin)%matrix)
         CALL cp_fm_create(vmos(ispin)%matrix, fm_struct)
         CALL cp_dbcsr_sm_fm_multiply(matrix_v(ispin, 1)%matrix, mo_coeff, vmos(ispin)%matrix, &
                                      norb, alpha=focc, beta=0.0_dp)
      END DO
      !MK The Intel compiler has trouble with that pointer assignment
      !MK ec_env%matrix_hz => matrix_v(:, 1)
      !MK Begin explicit copy
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_hz, nspins)
      DO ispin = 1, nspins
         ALLOCATE (ec_env%matrix_hz(ispin)%matrix)
         CALL dbcsr_create(ec_env%matrix_hz(ispin)%matrix, template=matrix_v(ispin, 1)%matrix)
         CALL dbcsr_copy(matrix_v(ispin, 1)%matrix, ec_env%matrix_hz(ispin)%matrix)
      END DO
      CALL dbcsr_deallocate_matrix_set(matrix_v)
      !MK End
      IF (ASSOCIATED(ec_env%cpmos)) THEN
         DO i = 1, SIZE(ec_env%cpmos)
            CALL cp_fm_release(ec_env%cpmos(i)%matrix)
         END DO
         DEALLOCATE (ec_env%cpmos)
         NULLIFY (ec_env%cpmos)
      END IF
      ec_env%cpmos => vmos
      !
      CALL pw_pool_give_back_pw(auxbas_pw_pool, dv_hartree_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, vtot_rspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoout_r(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoout_g(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      END DO
      DEALLOCATE (rhoout_r, rhoout_g, v_xc)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)

      ! initialize src matrix
      NULLIFY (scrm)
      CALL dbcsr_allocate_matrix_set(scrm, 1, 1)
      ALLOCATE (scrm(1, 1)%matrix)
      CALL dbcsr_create(scrm(1, 1)%matrix, template=ec_env%matrix_s(1, 1)%matrix)
      CALL cp_dbcsr_alloc_block_from_nbl(scrm(1, 1)%matrix, ec_env%sab_orb)

      ! v_rspace and v_tau_rspace are generated from the auxbas pool
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
         ENDDO
      END IF

      CALL qs_rho_get(rho, rho_r=rho_r)
      IF (ASSOCIATED(v_tau_rspace)) THEN
         CALL qs_rho_get(rho, tau_r=tau_r)
      END IF
      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      DO ispin = 1, nspins
         ! Add v_hartree + v_xc = v_rspace
         CALL pw_scale(v_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(v_hartree_rspace%pw, v_rspace(ispin)%pw)
         ! integrate over potential <a|V|b>
         CALL integrate_v_rspace(v_rspace=v_rspace(ispin), &
                                 hmat=scrm(1, 1), &
                                 pmat=ec_env%matrix_p(ispin, 1), &
                                 qs_env=qs_env, &
                                 calculate_forces=.TRUE., &
                                 basis_type="HARRIS", &
                                 task_list_external=ec_env%task_list)

         IF (ASSOCIATED(v_tau_rspace)) THEN
            ! integrate over Tau-potential <nabla.a|V|nabla.b>
            CALL pw_scale(v_tau_rspace(ispin)%pw, v_tau_rspace(ispin)%pw%pw_grid%dvol)
            CALL integrate_v_rspace(v_rspace=v_tau_rspace(ispin), hmat=scrm(1, 1), &
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

      ! delete scr matrix
      CALL dbcsr_deallocate_matrix_set(scrm)

      ! return pw grids
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
         IF (ASSOCIATED(v_tau_rspace)) THEN
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
         END IF
      ENDDO

      ! Core overlap
      IF (debug_forces) fodeb(1:3) = force(1)%core_overlap(1:3, 1)
      CALL calculate_ecore_overlap(qs_env, para_env, .TRUE., E_overlap_core=eovrl)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%core_overlap(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: CoreOverlap", fodeb
      END IF

      IF (debug_forces) THEN
         CALL get_qs_env(qs_env, natom=natom, atomic_kind_set=atomic_kind_set)
         ALLOCATE (ftot(3, natom))
         CALL total_qs_force(ftot, force, atomic_kind_set)
         fodeb(1:3) = ftot(1:3, 1)
         DEALLOCATE (ftot)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Force Explicit", fodeb
      END IF

      DEALLOCATE (v_rspace)

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
      ENDIF

      SELECT CASE (ec_env%ks_solver)
      CASE (ec_diagonalization)
         CALL ec_diag_solver(qs_env, ksmat, smat, pmat, wmat)
      CASE (ec_curvy_steps, ec_matrix_sign, ec_matrix_trs4, ec_matrix_tc2)
         CALL ec_ls_solver(qs_env, ec_env%ks_solver, ksmat, smat, pmat, wmat)
      CASE DEFAULT
         CPASSERT(.FALSE.)
      END SELECT

      IF (ec_env%mao) THEN
         CALL mao_release_matrices(ec_env, ksmat, smat, pmat, wmat)
      ENDIF

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
         CALL dbcsr_create(pmat(ispin, 1)%matrix, template=smat(1, 1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(pmat(ispin, 1)%matrix, ec_env%sab_orb)
      END DO

      CALL dbcsr_allocate_matrix_set(wmat, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (wmat(ispin, 1)%matrix)
         CALL dbcsr_create(wmat(ispin, 1)%matrix, template=smat(1, 1)%matrix)
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
                             1.0_dp, matrix_p(ispin, 1)%matrix, retain_sparsity=.TRUE., last_k=nmo(ispin))
         ! energy weighted density matrix
         CALL dbcsr_set(matrix_w(ispin, 1)%matrix, 0.0_dp)
         CALL cp_fm_column_scale(fm_mo, eigenvalues)
         CALL copy_fm_to_dbcsr(fm_mo, buf1_dbcsr)
         CALL dbcsr_multiply("N", "N", 1.0_dp, ortho_dbcsr, buf1_dbcsr, &
                             0.0_dp, buf3_dbcsr, filter_eps=eps_filter)
         CALL dbcsr_multiply("N", "T", focc(ispin), buf2_dbcsr, buf3_dbcsr, &
                             1.0_dp, matrix_w(ispin, 1)%matrix, retain_sparsity=.TRUE., last_k=nmo(ispin))
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
!> \brief Solve KS equation using linear scaling methods
!> \param qs_env ...
!> \param ls_method ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param matrix_p ...
!> \param matrix_w ...
!> \par History
!>      12.2019 created [JGH]
!> \author JGH
! **************************************************************************************************

   SUBROUTINE ec_ls_solver(qs_env, ls_method, matrix_ks, matrix_s, matrix_p, matrix_w)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ls_method
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s, matrix_p, matrix_w

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ec_ls_solver'

      INTEGER :: cluster_type, extrapolation_order, handle, ispin, max_iter_lanczos, nelectron, &
         nspins, s_inversion_type, s_preconditioner_type, s_sqrt_method, s_sqrt_order, &
         sign_method, sign_order
      INTEGER, DIMENSION(2)                              :: nmo
      LOGICAL                                            :: check_s_inv, converged, &
                                                            dynamic_threshold, fixed_mu, &
                                                            non_monotonic, report_all_sparsities, &
                                                            single_precision
      REAL(KIND=dp)                                      :: eps_filter, eps_lanczos, mu
      REAL(KIND=dp), DIMENSION(2)                        :: e_homo, e_lumo, mu_spin
      TYPE(dbcsr_type)                                   :: matrix_s_inv, matrix_s_sqrt, &
                                                            matrix_s_sqrt_inv, matrix_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input, solver_section

      CALL timeset(routineN, handle)

      CPASSERT(ls_method > 0)
      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_s))
      CPASSERT(ASSOCIATED(matrix_p))
      CPASSERT(ASSOCIATED(matrix_w))

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, input=input)
      nspins = dft_control%nspins
      !
      solver_section => section_vals_get_subs_vals(input, "DFT%ENERGY_CORRECTION%LS_SOLVER")
      CALL section_vals_val_get(solver_section, "EPS_FILTER", r_val=eps_filter)
      CALL section_vals_val_get(solver_section, "MU", r_val=mu)
      CALL section_vals_val_get(solver_section, "FIXED_MU", l_val=fixed_mu)
      mu_spin = mu
      CALL section_vals_val_get(solver_section, "S_PRECONDITIONER", i_val=s_preconditioner_type)
      CALL section_vals_val_get(solver_section, "MATRIX_CLUSTER_TYPE", i_val=cluster_type)
      CALL section_vals_val_get(solver_section, "SINGLE_PRECISION_MATRICES", l_val=single_precision)
      CALL section_vals_val_get(solver_section, "S_INVERSION", i_val=s_inversion_type)
      CALL section_vals_val_get(solver_section, "CHECK_S_INV", l_val=check_s_inv)
      CALL section_vals_val_get(solver_section, "REPORT_ALL_SPARSITIES", l_val=report_all_sparsities)
      CALL section_vals_val_get(solver_section, "SIGN_METHOD", i_val=sign_method)
      CALL section_vals_val_get(solver_section, "SIGN_ORDER", i_val=sign_order)
      CALL section_vals_val_get(solver_section, "DYNAMIC_THRESHOLD", l_val=dynamic_threshold)
      CALL section_vals_val_get(solver_section, "NON_MONOTONIC", l_val=non_monotonic)
      CALL section_vals_val_get(solver_section, "S_SQRT_METHOD", i_val=s_sqrt_method)
      CALL section_vals_val_get(solver_section, "S_SQRT_ORDER", i_val=s_sqrt_order)
      CALL section_vals_val_get(solver_section, "EXTRAPOLATION_ORDER", i_val=extrapolation_order)
      CALL section_vals_val_get(solver_section, "EPS_LANCZOS", r_val=eps_lanczos)
      CALL section_vals_val_get(solver_section, "MAX_ITER_LANCZOS", i_val=max_iter_lanczos)
      !
      nmo = 0
      CALL get_qs_env(qs_env=qs_env, nelectron_spin=nmo)
      IF (nspins == 1) nmo(1) = nmo(1)/2
      e_homo = -0.1_dp
      e_lumo = 0.1_dp

      SELECT CASE (ls_method)
      CASE (ec_curvy_steps)
         CPABORT("Curvy Step Method not available")
      CASE (ec_matrix_sign)
         CALL dbcsr_create(matrix_s_inv, template=matrix_s(1, 1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL invert_Hotelling(matrix_s_inv, matrix_s(1, 1)%matrix, eps_filter)
      CASE (ec_matrix_trs4, ec_matrix_tc2)
         CALL dbcsr_create(matrix_s_sqrt, template=matrix_s(1, 1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_s_sqrt_inv, template=matrix_s(1, 1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         SELECT CASE (s_sqrt_method)
         CASE (ls_s_sqrt_proot)
            CALL matrix_sqrt_proot(matrix_s_sqrt, matrix_s_sqrt_inv, &
                                   matrix_s(1, 1)%matrix, eps_filter, &
                                   s_sqrt_order, eps_lanczos, max_iter_lanczos, symmetrize=.TRUE.)
         CASE (ls_s_sqrt_ns)
            CALL matrix_sqrt_Newton_Schulz(matrix_s_sqrt, matrix_s_sqrt_inv, &
                                           matrix_s(1, 1)%matrix, eps_filter, &
                                           s_sqrt_order, eps_lanczos, max_iter_lanczos)
         CASE DEFAULT
            CPABORT("Unknown sqrt method.")
         END SELECT
         CALL dbcsr_release(matrix_s_sqrt)
      CASE DEFAULT
         CPABORT("Wrong ls_method")
      END SELECT

      DO ispin = 1, nspins
         nelectron = nmo(ispin)

         SELECT CASE (ls_method)
         CASE (ec_curvy_steps)
            CPABORT("Curvy Step Method not available")
         CASE (ec_matrix_sign)
            CALL density_matrix_sign(matrix_p(ispin, 1)%matrix, mu_spin(ispin), fixed_mu, &
                                     sign_method, sign_order, matrix_ks(ispin, 1)%matrix, &
                                     matrix_s(1, 1)%matrix, matrix_s_inv, nelectron, eps_filter)
         CASE (ec_matrix_trs4)
            CALL density_matrix_trs4(matrix_p(ispin, 1)%matrix, matrix_ks(ispin, 1)%matrix, &
                                     matrix_s_sqrt_inv, nelectron, eps_filter, &
                                     e_homo(ispin), e_lumo(ispin), mu_spin(ispin), &
                                     max_iter_lanczos=max_iter_lanczos, &
                                     eps_lanczos=eps_lanczos, converged=converged)
         CASE (ec_matrix_tc2)
            CPABORT("TC2 Method not available")
            CALL density_matrix_tc2(matrix_p(ispin, 1)%matrix, matrix_ks(ispin, 1)%matrix, &
                                    matrix_s_sqrt_inv, nelectron, eps_filter, &
                                    e_homo(ispin), e_lumo(ispin), &
                                    non_monotonic, eps_lanczos, max_iter_lanczos)
         CASE DEFAULT
            CPABORT("Wrong ls_method")
         END SELECT
      END DO

      SELECT CASE (ls_method)
      CASE (ec_curvy_steps)
         CPABORT("Curvy Step Method not available")
      CASE (ec_matrix_sign)
         CALL dbcsr_release(matrix_s_inv)
      CASE (ec_matrix_trs4, ec_matrix_tc2)
         CALL dbcsr_release(matrix_s_sqrt_inv)
      CASE DEFAULT
         CPABORT("Wrong ls_method")
      END SELECT

      ! W matrix
      CALL dbcsr_create(matrix_tmp, template=matrix_s(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      DO ispin = 1, nspins
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p(ispin, 1)%matrix, matrix_ks(ispin, 1)%matrix, &
                             0.0_dp, matrix_tmp, filter_eps=eps_filter)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp, matrix_p(ispin, 1)%matrix, &
                             0.0_dp, matrix_w(ispin, 1)%matrix, filter_eps=eps_filter)
      ENDDO
      CALL dbcsr_release(matrix_tmp)

      IF (nspins == 1) THEN
         CALL dbcsr_scale(matrix_p(1, 1)%matrix, 2.0_dp)
         CALL dbcsr_scale(matrix_w(1, 1)%matrix, 2.0_dp)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_ls_solver

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
            CALL dbcsr_dot(ec_env%matrix_ks(ispin, 1)%matrix, ec_env%matrix_p(ispin, 1)%matrix, trace)
            eband = eband + trace
         END DO
         ec_env%eband = eband + ec_env%efield_nuclear
         ec_env%etotal = ec_env%eband + ec_env%ehartree + ec_env%exc - ec_env%vhxc + ec_env%edispersion
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T3,A,T65,F16.10)') "Eband    ", ec_env%eband
            WRITE (unit_nr, '(T3,A,T65,F16.10)') "Ehartree ", ec_env%ehartree
            WRITE (unit_nr, '(T3,A,T65,F16.10)') "Exc      ", ec_env%exc
            WRITE (unit_nr, '(T3,A,T65,F16.10)') "Evhxc    ", ec_env%vhxc
            WRITE (unit_nr, '(T3,A,T65,F16.10)') "Edisp    ", ec_env%edispersion
            WRITE (unit_nr, '(T3,A,T65,F16.10)') "Etotal Harris Functional   ", ec_env%etotal
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

      INTEGER                                            :: handle, ikind, nkind
      LOGICAL                                            :: gth_potential_present, &
                                                            sgp_potential_present, &
                                                            skip_load_balance_distributed
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: orb_present, ppl_present, ppnl_present
      REAL(dp)                                           :: subcells
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: orb_radius, ppl_radius, ppnl_radius
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
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind_set(qs_kind_set, gth_potential_present=gth_potential_present, &
                           sgp_potential_present=sgp_potential_present)
      nkind = SIZE(qs_kind_set)
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

      ! Release work storage
      CALL atom2d_cleanup(atom2d)
      DEALLOCATE (atom2d)
      DEALLOCATE (orb_present, ppl_present, ppnl_present)
      DEALLOCATE (orb_radius, ppl_radius, ppnl_radius)
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
      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: akind, handle, i, ia, iatom, idir, &
                                                            ikind, iounit, ispin, maxmom, nspins, &
                                                            reference, unit_nr
      LOGICAL                                            :: magnetic, periodic
      REAL(KIND=dp)                                      :: charge, dd, focc, tmp
      REAL(KIND=dp), DIMENSION(3)                        :: cdip, pdip, rcc, rdip, ria, tdip
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, moments
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      rlab(1) = "X"
      rlab(2) = "Y"
      rlab(3) = "Z"

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_unit_nr(logger)

      print_key => section_vals_get_subs_vals(section_vals=qs_env%input, &
                                              subsection_name="DFT%PRINT%MOMENTS")

      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN

         maxmom = section_get_ival(section_vals=qs_env%input, &
                                   keyword_name="DFT%PRINT%MOMENTS%MAX_MOMENT")
         periodic = section_get_lval(section_vals=qs_env%input, &
                                     keyword_name="DFT%PRINT%MOMENTS%PERIODIC")
         reference = section_get_ival(section_vals=qs_env%input, &
                                      keyword_name="DFT%PRINT%MOMENTS%REFERENCE")
         magnetic = section_get_lval(section_vals=qs_env%input, &
                                     keyword_name="DFT%PRINT%MOMENTS%MAGNETIC")
         NULLIFY (ref_point)
         CALL section_vals_val_get(qs_env%input, "DFT%PRINT%MOMENTS%REF_POINT", r_vals=ref_point)
         unit_nr = cp_print_key_unit_nr(logger=logger, basis_section=qs_env%input, &
                                        print_key_path="DFT%PRINT%MOMENTS", extension=".dat", &
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
                     CALL dbcsr_dot(ec_env%p_env%p1(ispin)%matrix, moments(idir)%matrix, tmp)
                     rdip(idir) = rdip(idir) + tmp
                  END DO
               END DO
               CALL dbcsr_deallocate_matrix_set(moments)
               !
               IF (unit_nr > 0) THEN
                  tdip = -(rdip + pdip + cdip)
                  WRITE (unit_nr, "(T3,A)") "Dipoles are based on the traditional operator."
                  dd = SQRT(SUM(tdip(1:3)**2))*debye
                  WRITE (unit_nr, "(T3,A)") "Dipole moment [Debye]"
                  WRITE (unit_nr, "(T5,3(A,A,F14.8,1X),T60,A,T67,F14.8)") &
                     (TRIM(rlab(i)), "=", tdip(i)*debye, i=1, 3), "Total=", dd
               END IF
            ENDIF
         END IF

         CALL cp_print_key_finished_output(unit_nr=unit_nr, logger=logger, &
                                           basis_section=qs_env%input, print_key_path="DFT%PRINT%MOMENTS")
      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_properties

! **************************************************************************************************

END MODULE energy_corrections
