!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for a Harris type energy correction on top of a
!>        Kim-Gordon calculation
!> \par History
!>       03.2014 created
!> \author JGH
! *****************************************************************************
MODULE kg_energy_corrections
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE core_ppl,                        ONLY: build_core_ppl
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_add, cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, cp_dbcsr_create, &
        cp_dbcsr_deallocate_matrix_set, cp_dbcsr_desymmetrize, cp_dbcsr_filter, cp_dbcsr_get_info, &
        cp_dbcsr_init, cp_dbcsr_init_p, cp_dbcsr_multiply, cp_dbcsr_p_type, cp_dbcsr_release, &
        cp_dbcsr_set, cp_dbcsr_trace, cp_dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_triangular_invert
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_triangular,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE input_constants,                 ONLY: kg_ec_diagonalization,&
                                              kg_ec_functional_harris
   USE kg_environment_types,            ONLY: energy_correction_type,&
                                              kg_environment_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_transfer
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
                                              pw_p_type
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ecore_self
   USE qs_dispersion_pairpot,           ONLY: calculate_dispersion_pairpot
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_energy_corrections'

   PUBLIC :: kg_energy_correction

CONTAINS

! *****************************************************************************
!> \brief Energy correction to a KG simulation
!>
!> \param qs_env ...
!> \param calculate_forces ...
!> \par History
!>       03.2014 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE kg_energy_correction(qs_env, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'kg_energy_correction', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, unit_nr
      LOGICAL                                            :: my_calc_forces
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(kg_environment_type), POINTER                 :: kg_env

      CALL timeset(routineN, handle)

      my_calc_forces = .FALSE.
      IF (PRESENT(calculate_forces)) my_calc_forces = calculate_forces

      NULLIFY (ec_env, kg_env)
      CALL get_qs_env(qs_env=qs_env, kg_env=kg_env)

      ! Check for energy correction
      IF (kg_env%energy_correction) THEN

         ec_env => kg_env%ec_env

         ec_env%etotal = 0.0_dp
         ec_env%eband = 0.0_dp
         ec_env%ehartree = 0.0_dp
         ec_env%exc = 0.0_dp
         ec_env%vhxc = 0.0_dp
         ec_env%edispersion = 0.0_dp

         logger => cp_get_default_logger()
         IF (logger%para_env%mepos == logger%para_env%source) THEN
            unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
         ELSE
            unit_nr = -1
         ENDIF

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(/,T2,A)') '!-----------------------------------------------------------------------------!'
            WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 30), " KG energy correction ", REPEAT("-", 29)
         END IF
         CALL build_core_hamiltonian(qs_env, kg_env, my_calc_forces)
         CALL build_ks_matrix(qs_env, kg_env, my_calc_forces, unit_nr)
         CALL ks_solver(qs_env, kg_env, unit_nr)
         CALL ec_energy(qs_env, kg_env%ec_env, unit_nr)
!          CALL kg_response_solver(qs_env,kg_env,calculate_forces,error)

      END IF

      CALL timestop(handle)

   END SUBROUTINE kg_energy_correction

! *****************************************************************************
!> \brief Construction of the Core Hamiltonian Matrix
!>        Short version of qs_core_hamiltonian
!> \param qs_env ...
!> \param kg_env ...
!> \param calculate_forces ...
!> \author Creation (03.2014,JGH)
! **************************************************************************************************
   SUBROUTINE build_core_hamiltonian(qs_env, kg_env, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kg_environment_type), POINTER                 :: kg_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_core_hamiltonian', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nder, nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: eps_filter, eps_ppnl
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: matrix_p
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ae, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial

      IF (calculate_forces) THEN
         CALL timeset(routineN//"_forces", handle)
      ELSE
         CALL timeset(routineN, handle)
      ENDIF

      NULLIFY (ec_env)
      ec_env => kg_env%ec_env
      NULLIFY (atomic_kind_set, dft_control, ks_env, particle_set, qs_kind_set)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, dft_control=dft_control, &
                      ks_env=ks_env, particle_set=particle_set, qs_kind_set=qs_kind_set)

      nimages = dft_control%nimages
      CPASSERT(nimages == 1)

      ! check for virial
      CALL get_qs_env(qs_env=qs_env, virial=virial)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      CPASSERT(.NOT. use_virial)

      ! get neighbor lists, we need the full sab_orb list from the ec_env
      ! the other lists are non-molecular and OK from the qs_env
      NULLIFY (sab_orb, sac_ae, sac_ppl, sap_ppnl)
!    sab_orb => ec_env%sab_orb
      sab_orb => kg_env%sab_orb_full
      CALL get_qs_env(qs_env=qs_env, &
                      sac_ae=sac_ae, sac_ppl=sac_ppl, sap_ppnl=sap_ppnl)

      NULLIFY (matrix_p, force)
      nder = 0
      eps_filter = dft_control%qs_control%eps_filter_matrix

      CALL build_overlap_matrix(ks_env, nderivative=nder, matrixkp_s=ec_env%matrix_s, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", &
                                basis_type_b="ORB", &
                                sab_nl=sab_orb, qs_input=qs_env%input)

      CALL build_kinetic_matrix(ks_env, matrixkp_t=ec_env%matrix_t, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="ORB", &
                                sab_nl=sab_orb, &
                                eps_filter=eps_filter)

      ! initialize H matrix
      CALL cp_dbcsr_allocate_matrix_set(ec_env%matrix_h, 1, 1)
      ALLOCATE (ec_env%matrix_h(1, 1)%matrix)
      CALL cp_dbcsr_init(ec_env%matrix_h(1, 1)%matrix)
      CALL cp_dbcsr_create(ec_env%matrix_h(1, 1)%matrix, template=ec_env%matrix_s(1, 1)%matrix)
      CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_h(1, 1)%matrix, sab_orb)

      ! add kinetic energy
      CALL cp_dbcsr_copy(ec_env%matrix_h(1, 1)%matrix, ec_env%matrix_t(1, 1)%matrix, &
                         keep_sparsity=.TRUE., name="CORE HAMILTONIAN MATRIX")

      ! compute the ppl contribution to the core hamiltonian
      nimages = 1
      NULLIFY (cell_to_index)
      use_virial = .FALSE.
      IF (ASSOCIATED(sac_ppl)) THEN
         CALL build_core_ppl(ec_env%matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ppl, nimages, cell_to_index)
      END IF

      ! compute the ppnl contribution to the core hamiltonian ***
      eps_ppnl = dft_control%qs_control%eps_ppnl
      CALL build_core_ppnl(ec_env%matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder, &
                           qs_kind_set, atomic_kind_set, particle_set, sab_orb, sap_ppnl, eps_ppnl, nimages, cell_to_index)

      CALL timestop(handle)

   END SUBROUTINE build_core_hamiltonian

! *****************************************************************************
!> \brief calculate the complete KS matrix
!> \param qs_env ...
!> \param kg_env ...
!> \param calculate_forces ...
!> \param unit_nr ...
!> \par History
!>      03.2014 adapted from qs_ks_build_kohn_sham_matrix [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE build_ks_matrix(qs_env, kg_env, calculate_forces, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kg_environment_type), POINTER                 :: kg_env
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_ks_matrix', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: headline
      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: use_virial
      REAL(dp)                                           :: eexc, ehartree, eovrl, eself, evhxc, &
                                                            trace
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: rho_ao
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_rspace, v_tau_rspace
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      ! get all information on the electronic density
      NULLIFY (rho, ks_env)
      CALL get_qs_env(qs_env=qs_env, rho=rho, virial=virial, dft_control=dft_control, &
                      para_env=para_env, ks_env=ks_env)

      nspins = dft_control%nspins
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      ! Kohn-Sham matrix
      ec_env => kg_env%ec_env

      IF (ASSOCIATED(ec_env%matrix_ks)) CALL cp_dbcsr_deallocate_matrix_set(ec_env%matrix_ks)

      CALL cp_dbcsr_allocate_matrix_set(ec_env%matrix_ks, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (ec_env%matrix_ks(ispin, 1)%matrix)
         CALL cp_dbcsr_init(ec_env%matrix_ks(ispin, 1)%matrix)

         IF (nspins > 1) THEN
            IF (ispin == 1) THEN
               headline = "KOHN-SHAM MATRIX FOR ALPHA SPIN"
            ELSE
               headline = "KOHN-SHAM MATRIX FOR BETA SPIN"
            END IF
         ELSE
            headline = "KOHN-SHAM MATRIX"
         END IF
         CALL cp_dbcsr_create(ec_env%matrix_ks(ispin, 1)%matrix, name=TRIM(headline), &
                              template=ec_env%matrix_s(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_ks(ispin, 1)%matrix, kg_env%sab_orb_full)
         CALL cp_dbcsr_set(ec_env%matrix_ks(ispin, 1)%matrix, 0.0_dp)
      ENDDO

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CPASSERT(ASSOCIATED(pw_env))

      NULLIFY (auxbas_pw_pool, poisson_env, pw_pools)
      ! gets the tmp grids
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      pw_pools=pw_pools, poisson_env=poisson_env)

      ! Calculate the Hartree potential
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             v_hartree_rspace%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      ! Get the total density in g-space [ions + electrons]
      CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree, &
                            v_hartree_gspace%pw)
      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)

      ! v_rspace and v_tau_rspace are generated from the auxbas pool
      NULLIFY (v_rspace, v_tau_rspace)
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=ec_env%xc_section, &
                         vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=eexc, just_energy=.FALSE.)

      ! Add v_hartree + v_xc = v_rspace
      ! integrate over potentials
      DO ispin = 1, nspins
         CALL pw_scale(v_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(v_hartree_rspace%pw, v_rspace(ispin)%pw)

         CALL integrate_v_rspace(v_rspace=v_rspace(ispin), &
                                 hmat=ec_env%matrix_ks(ispin, 1), &
                                 qs_env=qs_env, &
                                 calculate_forces=.FALSE., &
                                 task_list_external=ec_env%task_list)

         IF (ASSOCIATED(v_tau_rspace)) THEN
            CALL pw_scale(v_tau_rspace(ispin)%pw, v_tau_rspace(ispin)%pw%pw_grid%dvol)

            CALL integrate_v_rspace(v_rspace=v_tau_rspace(ispin), hmat=ec_env%matrix_ks(ispin, 1), &
                                    qs_env=qs_env, calculate_forces=.FALSE., compute_tau=.TRUE., &
                                    task_list_external=ec_env%task_list)
         END IF

      END DO

      ! return pw grids
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
         IF (ASSOCIATED(v_tau_rspace)) THEN
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
         END IF
      ENDDO

      ! matrix_ks contains Vhxc at this point
      ! calculate Tr(Vhxc P0)
      evhxc = 0.0_dp
      CALL qs_rho_get(rho, rho_ao=rho_ao)

      DO ispin = 1, nspins
         CALL cp_dbcsr_trace(ec_env%matrix_ks(ispin, 1)%matrix, rho_ao(ispin)%matrix, trace)
         evhxc = evhxc+trace
      END DO

      ! energies
      CALL calculate_ecore_self(qs_env, E_self_core=eself)
      CALL calculate_ecore_overlap(qs_env, para_env, calculate_forces=.FALSE., E_overlap_core=eovrl)
      ec_env%exc = eexc
      ec_env%ehartree = ehartree+eovrl+eself
      ec_env%vhxc = evhxc
!dbg
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A,T16,F16.10)') 'Tr[KP0] ', trace
         WRITE (unit_nr, '(T2,A,T16,F16.10)') 'eexc ', eexc
         WRITE (unit_nr, '(T2,A,T16,F16.10)') 'evhxc ', evhxc
         WRITE (unit_nr, '(T2,A,T16,F16.10)') 'eovrl ', eovrl
         WRITE (unit_nr, '(T2,A,T16,F16.10)') 'eself ', eself
         WRITE (unit_nr, '(T2,A,T16,F16.10)') 'ehartree ', ehartree
      END IF
!dbg

      ! add the core matrix
      DO ispin = 1, nspins
         CALL cp_dbcsr_add(ec_env%matrix_ks(ispin, 1)%matrix, ec_env%matrix_h(1, 1)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END DO

      ! At this point the ks matrix is up to date, filter it if requested
      DO ispin = 1, nspins
         CALL cp_dbcsr_filter(ec_env%matrix_ks(ispin, 1)%matrix, &
                              dft_control%qs_control%eps_filter_matrix)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE build_ks_matrix

! *****************************************************************************
!> \brief  Solve KS equation for a given matrix
!> \param qs_env ...
!> \param kg_env ...
!> \param calculate_forces ...
!> \param unit_nr ...
!> \par History
!>      03.2014 adapted from qs_ks_build_kohn_sham_matrix [JGH]
!> \author JGH
! **************************************************************************************************

   SUBROUTINE ks_solver(qs_env, kg_env, unit_nr)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kg_environment_type), POINTER                 :: kg_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ks_solver', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: headline
      INTEGER                                            :: handle, info, ispin, nmo(2), nsize, &
                                                            nspins
      REAL(KIND=dp)                                      :: eps_filter, focc(2)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_type), POINTER                       :: buf1_dbcsr, buf2_dbcsr, ortho_dbcsr, &
                                                            ref_matrix
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_ks, fm_mo, fm_ortho
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(energy_correction_type), POINTER              :: ec_env

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, dft_control, ec_env, para_env)
      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env, dft_control=dft_control, &
                      para_env=para_env)

      nspins = dft_control%nspins
      ec_env => kg_env%ec_env

      nmo = 0
      CALL get_qs_env(qs_env=qs_env, nelectron_spin=nmo)
      focc = 1._dp
      IF (nspins == 1) THEN
         focc = 2._dp
         nmo(1) = nmo(1)/2
      END IF

      IF (.NOT. ASSOCIATED(ec_env%matrix_p)) THEN
         CALL cp_dbcsr_allocate_matrix_set(ec_env%matrix_p, nspins, 1)
         DO ispin = 1, nspins
            ALLOCATE (ec_env%matrix_p(ispin, 1)%matrix)
            IF (nspins > 1) THEN
               IF (ispin == 1) THEN
                  headline = "DENSITY MATRIX FOR ALPHA SPIN"
               ELSE
                  headline = "DENSITY MATRIX FOR BETA SPIN"
               END IF
            ELSE
               headline = "DENSITY MATRIX"
            END IF
            CALL cp_dbcsr_init(ec_env%matrix_p(ispin, 1)%matrix)
            CALL cp_dbcsr_create(ec_env%matrix_p(ispin, 1)%matrix, name=TRIM(headline), template=ec_env%matrix_s(1, 1)%matrix)
            CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_p(ispin, 1)%matrix, kg_env%sab_orb_full)
         END DO
      END IF

      SELECT CASE (ec_env%ks_solver)
      CASE (kg_ec_diagonalization)
         eps_filter = dft_control%qs_control%eps_filter_matrix

         CALL cp_dbcsr_get_info(ec_env%matrix_ks(1, 1)%matrix, nfullrows_total=nsize)
         ALLOCATE (eigenvalues(nsize))

         NULLIFY (fm_ortho, fm_ks, fm_mo, fm_struct, ref_matrix)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nsize, &
                                  ncol_global=nsize, para_env=para_env)
         CALL cp_fm_create(fm_ortho, fm_struct)
         CALL cp_fm_create(fm_ks, fm_struct)
         CALL cp_fm_create(fm_mo, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         ! factorization
         ref_matrix => ec_env%matrix_s(1, 1)%matrix
         CALL cp_dbcsr_init_p(ortho_dbcsr)
         CALL cp_dbcsr_create(ortho_dbcsr, template=ref_matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_init_p(buf1_dbcsr)
         CALL cp_dbcsr_create(buf1_dbcsr, template=ref_matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_init_p(buf2_dbcsr)
         CALL cp_dbcsr_create(buf2_dbcsr, template=ref_matrix, &
                              matrix_type=dbcsr_type_no_symmetry)

         CALL copy_dbcsr_to_fm(ref_matrix, fm_ortho)
         CALL cp_fm_cholesky_decompose(fm_ortho)
         CALL cp_fm_triangular_invert(fm_ortho)
         CALL cp_fm_set_all(fm_ks, 0.0_dp)
         CALL cp_fm_to_fm_triangular(fm_ortho, fm_ks, "U")
         CALL copy_fm_to_dbcsr(fm_ks, ortho_dbcsr)
         CALL cp_fm_set_all(fm_ks, 0.0_dp)

         DO ispin = 1, nspins
            ! calculate ZHZ(T)
            ! calculate Z(T)HZ
            CALL cp_dbcsr_desymmetrize(ec_env%matrix_ks(ispin, 1)%matrix, buf1_dbcsr)
            CALL cp_dbcsr_multiply("N", "N", 1.0_dp, buf1_dbcsr, ortho_dbcsr, &
                                   0.0_dp, buf2_dbcsr, filter_eps=eps_filter)
            CALL cp_dbcsr_multiply("T", "N", 1.0_dp, ortho_dbcsr, buf2_dbcsr, &
                                   0.0_dp, buf1_dbcsr, filter_eps=eps_filter)
            ! copy to fm format
            CALL copy_dbcsr_to_fm(buf1_dbcsr, fm_ks)
            CALL choose_eigv_solver(fm_ks, fm_mo, eigenvalues, info)
            CPASSERT(info == 0)
            ! back transform of mos c = Z(T)*c
            CALL copy_fm_to_dbcsr(fm_mo, buf1_dbcsr)
            CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ortho_dbcsr, buf1_dbcsr, &
                                   0.0_dp, buf2_dbcsr, filter_eps=eps_filter)
            ! density matrix
            CALL cp_dbcsr_set(ec_env%matrix_p(ispin, 1)%matrix, 0.0_dp)
            CALL cp_dbcsr_multiply("N", "T", focc(ispin), buf2_dbcsr, buf2_dbcsr, &
                                   1.0_dp, ec_env%matrix_p(ispin, 1)%matrix, retain_sparsity=.TRUE., last_k=nmo(ispin))

         END DO
         CALL cp_fm_release(fm_ks)
         CALL cp_fm_release(fm_mo)
         CALL cp_fm_release(fm_ortho)
         CALL cp_dbcsr_release(ortho_dbcsr)
         CALL cp_dbcsr_release(buf1_dbcsr)
         CALL cp_dbcsr_release(buf2_dbcsr)
         DEALLOCATE (eigenvalues)

!dbg
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A)') ' Corrected-KG density matrix done'
         END IF
!dbg

      CASE DEFAULT
         CPASSERT(.FALSE.)
      END SELECT
      CALL timestop(handle)

   END SUBROUTINE ks_solver

! *****************************************************************************
!> \brief Calculate the energy correction
!> \param qs_env ...
!> \param ec_env ...
!> \param unit_nr ...
!> \author Creation (03.2014,JGH)
! **************************************************************************************************
   SUBROUTINE ec_energy(qs_env, ec_env, unit_nr)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type)                       :: ec_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_energy', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: eband, energy, trace

      CALL timeset(routineN, handle)

      ! dispersion through pairpotentials
      CALL calculate_dispersion_pairpot(qs_env, ec_env%dispersion_env, energy, .FALSE.)
      ec_env%edispersion = ec_env%edispersion+energy

      SELECT CASE (ec_env%energy_functional)
      CASE (kg_ec_functional_harris)

         nspins = SIZE(ec_env%matrix_ks, 1)

         eband = 0.0_dp
         DO ispin = 1, nspins
!dbg
            CALL cp_dbcsr_trace(ec_env%matrix_p(ispin, 1)%matrix, ec_env%matrix_s(ispin, 1)%matrix, trace)
            IF (unit_nr > 0) WRITE (unit_nr, '(T2,A,T16,F16.10)') 'Tr[PS] ', trace
!dbg
            CALL cp_dbcsr_trace(ec_env%matrix_ks(ispin, 1)%matrix, ec_env%matrix_p(ispin, 1)%matrix, trace)
            eband = eband+trace
         END DO
         ec_env%eband = eband
         ec_env%etotal = ec_env%eband+ec_env%ehartree+ec_env%exc-ec_env%vhxc+ec_env%edispersion
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A,T16,F16.10)') "HF Etotal   ", ec_env%etotal
            WRITE (unit_nr, '(T2,A,T16,F16.10)') "Eband    ", ec_env%eband
            WRITE (unit_nr, '(T2,A,T16,F16.10)') "Ehartree ", ec_env%ehartree
            WRITE (unit_nr, '(T2,A,T16,F16.10)') "Exc      ", ec_env%exc
            WRITE (unit_nr, '(T2,A,T16,F16.10)') "Evhxc    ", ec_env%vhxc
            WRITE (unit_nr, '(T2,A,T16,F16.10)') "Edisp    ", ec_env%edispersion
         END IF

      CASE DEFAULT

         CPASSERT(.FALSE.)

      END SELECT

      CALL timestop(handle)

   END SUBROUTINE ec_energy

END MODULE kg_energy_corrections
