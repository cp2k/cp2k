!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate the CPKS equation and the resulting forces
!> \par History
!>       03.2014 created
!>       09.2019 Moved from KG to Kohn-Sham
!>       11.2019 Moved from energy_correction
!>       08.2020 AO linear response solver [fbelle]
!> \author JGH
! **************************************************************************************************
MODULE response_solver
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type
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
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_init_random,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_distribution_type, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE ec_env_types,                    ONLY: energy_correction_type
   USE ec_methods,                      ONLY: create_kernel,&
                                              ec_mos_init
   USE ec_orth_solver,                  ONLY: ec_response_ao
   USE exstates_types,                  ONLY: excited_energy_type
   USE hfx_derivatives,                 ONLY: derivatives_four_center
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: &
        do_admm_aux_exch_func_none, ec_ls_solver, ec_mo_solver, kg_tnadd_atomic, kg_tnadd_embed, &
        kg_tnadd_embed_ri, ls_s_sqrt_ns, ls_s_sqrt_proot, ot_precond_full_all, &
        ot_precond_full_kinetic, ot_precond_full_single, ot_precond_full_single_inverse, &
        ot_precond_none, ot_precond_s_inverse, precond_mlp
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kg_correction,                   ONLY: kg_ekin_subset
   USE kg_environment_types,            ONLY: kg_environment_type
   USE kg_tnadd_mat,                    ONLY: build_tnadd_mat
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE mathlib,                         ONLY: det_3x3
   USE message_passing,                 ONLY: mp_sum
   USE mulliken,                        ONLY: ao_charges
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: pascal
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
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
                                              pw_create,&
                                              pw_p_type,&
                                              pw_release
   USE qs_2nd_kernel_ao,                ONLY: build_dm_response
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ecore_self
   USE qs_density_matrices,             ONLY: calculate_whz_matrix,&
                                              calculate_wz_matrix
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: add_qs_force,&
                                              qs_force_type,&
                                              total_qs_force
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_methods,               ONLY: linres_solver
   USE qs_linres_types,                 ONLY: linres_control_create,&
                                              linres_control_release,&
                                              linres_control_type
   USE qs_matrix_pools,                 ONLY: mpools_rebuild_fm_pools
   USE qs_mo_methods,                   ONLY: make_basis_sm
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_force,&
                                              build_overlap_matrix
   USE qs_p_env_methods,                ONLY: p_env_create,&
                                              p_env_psi0_changed
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE task_list_types,                 ONLY: task_list_type
   USE virial_methods,                  ONLY: one_third_sum_diag
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
   USE xtb_ehess,                       ONLY: xtb_coulomb_hessian
   USE xtb_ehess_force,                 ONLY: calc_xtb_ehess_force
   USE xtb_matrices,                    ONLY: xtb_hab_force
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'response_solver'

   LOGICAL, PARAMETER                   :: debug_stress = .FALSE.
   LOGICAL, PARAMETER                   :: debug_forces = .FALSE.

   PUBLIC :: response_calculation, response_equation, response_force, response_force_xtb, &
             ks_ref_potential, admm_projection_derivative

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Initializes solver of linear response equation for Harris functional Lagrangian
!>
!> \param qs_env ...
!> \param ec_env Environemnt of Harris energy correction
!>
!> \date    01.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE response_calculation(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'response_calculation'

      INTEGER                                            :: handle, homo, ispin, nao, nmo, nspins, &
                                                            solver_method, unit_nr
      LOGICAL                                            :: should_stop
      REAL(KIND=dp)                                      :: focc
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: cpmos
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, sv
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_s_aux, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input, solver_section

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, logger, para_env, sab_orb, solver_section)

      ! Get useful output unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      input=input, &
                      matrix_s=matrix_s, &
                      para_env=para_env, &
                      sab_orb=sab_orb)
      nspins = dft_control%nspins

      ! initialize linres_control
      NULLIFY (linres_control)
      CALL linres_control_create(linres_control)
      linres_control%do_kernel = .TRUE.
      linres_control%lr_triplet = .FALSE.
      linres_control%converged = .FALSE.
      linres_control%energy_gap = 0.02_dp

      ! Read input
      solver_section => section_vals_get_subs_vals(input, "DFT%ENERGY_CORRECTION%RESPONSE_SOLVER")
      CALL section_vals_val_get(solver_section, "EPS", r_val=linres_control%eps)
      CALL section_vals_val_get(solver_section, "EPS_FILTER", r_val=linres_control%eps_filter)
      CALL section_vals_val_get(solver_section, "MAX_ITER", i_val=linres_control%max_iter)
      CALL section_vals_val_get(solver_section, "METHOD", i_val=solver_method)
      CALL section_vals_val_get(solver_section, "PRECONDITIONER", i_val=linres_control%preconditioner_type)
      CALL section_vals_val_get(solver_section, "RESTART", l_val=linres_control%linres_restart)
      CALL section_vals_val_get(solver_section, "RESTART_EVERY", i_val=linres_control%restart_every)
      CALL set_qs_env(qs_env, linres_control=linres_control)

      ! Write input section of response solver
      CALL response_solver_write_input(solver_section, linres_control, unit_nr)

      ! Allocate and initialize Z matrix, and energy weighted z-matrix
      ! Template is ground-state overlap matrix
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_wz, nspins)
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_z, nspins)
      DO ispin = 1, nspins
         ALLOCATE (ec_env%matrix_wz(ispin)%matrix)
         ALLOCATE (ec_env%matrix_z(ispin)%matrix)
         CALL dbcsr_create(ec_env%matrix_wz(ispin)%matrix, name="Wz MATRIX", &
                           template=matrix_s(1)%matrix)
         CALL dbcsr_create(ec_env%matrix_z(ispin)%matrix, name="Z MATRIX", &
                           template=matrix_s(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_wz(ispin)%matrix, sab_orb)
         CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_z(ispin)%matrix, sab_orb)
         CALL dbcsr_set(ec_env%matrix_wz(ispin)%matrix, 0.0_dp)
         CALL dbcsr_set(ec_env%matrix_z(ispin)%matrix, 0.0_dp)
      END DO

      ! Choose between MO-solver and AO-solver
      SELECT CASE (solver_method)
      CASE (ec_mo_solver)

         ! MO solver requires MO's calculated during ground-state,
         ! which are not available if linear scaling schemes are used.
         ! If ground-state calculation was performed with linear scaling methods,
         ! introduce MOs here
         IF (dft_control%qs_control%do_ls_scf) THEN

            ! Allocate and initialize MO environment
            CALL ec_mos_init(qs_env, matrix_s(1)%matrix)
            CALL get_qs_env(qs_env, mos=mos, rho=rho)

            ! Get ground-state density matrix
            CALL qs_rho_get(rho, rho_ao=rho_ao)

            DO ispin = 1, nspins
               CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                               mo_coeff=mo_coeff, &
                               nmo=nmo, nao=nao, homo=homo)

               CALL cp_fm_set_all(mo_coeff, 0.0_dp)
               CALL cp_fm_init_random(mo_coeff, nmo)

               NULLIFY (sv)
               CALL cp_fm_create(sv, mo_coeff%matrix_struct, "SV")
               ! multiply times PS
               ! PS*C(:,1:nomo)+C(:,nomo+1:nmo) (nomo=NINT(nelectron/maxocc))
               CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, mo_coeff, sv, nmo)
               CALL cp_dbcsr_sm_fm_multiply(rho_ao(ispin)%matrix, sv, mo_coeff, homo)
               CALL cp_fm_release(sv)
               ! and ortho the result
               CALL make_basis_sm(mo_coeff, nmo, matrix_s(1)%matrix)

               ! rebuilds fm_pools
               ! originally done in qs_env_setup, only when mos associated
               CALL get_qs_env(qs_env, blacs_env=blacs_env)
               CALL mpools_rebuild_fm_pools(qs_env%mpools, mos=mos, &
                                            blacs_env=blacs_env, para_env=para_env)
            END DO
         END IF

         ! initialized p_env
         CALL p_env_release(ec_env%p_env)
         CALL p_env_create(ec_env%p_env, qs_env, orthogonal_orbitals=.TRUE., &
                           linres_control=linres_control)
         CALL set_qs_env(qs_env, linres_control=linres_control)
         CALL p_env_psi0_changed(ec_env%p_env, qs_env)
         !
         p_env => ec_env%p_env
         !
         CALL get_qs_env(qs_env, matrix_s=matrix_s, sab_orb=sab_orb)
         CALL dbcsr_allocate_matrix_set(p_env%p1, nspins)
         CALL dbcsr_allocate_matrix_set(p_env%w1, nspins)
         DO ispin = 1, nspins
            ALLOCATE (p_env%p1(ispin)%matrix, p_env%w1(ispin)%matrix)
            CALL dbcsr_create(matrix=p_env%p1(ispin)%matrix, template=matrix_s(1)%matrix)
            CALL dbcsr_create(matrix=p_env%w1(ispin)%matrix, template=matrix_s(1)%matrix)
            CALL cp_dbcsr_alloc_block_from_nbl(p_env%p1(ispin)%matrix, sab_orb)
            CALL cp_dbcsr_alloc_block_from_nbl(p_env%w1(ispin)%matrix, sab_orb)
         END DO
         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, matrix_s_aux_fit=matrix_s_aux)
            CALL dbcsr_allocate_matrix_set(p_env%p1_admm, nspins)
            DO ispin = 1, nspins
               ALLOCATE (p_env%p1_admm(ispin)%matrix)
               CALL dbcsr_create(p_env%p1_admm(ispin)%matrix, &
                                 template=matrix_s_aux(1)%matrix)
               CALL dbcsr_copy(p_env%p1_admm(ispin)%matrix, matrix_s_aux(1)%matrix)
               CALL dbcsr_set(p_env%p1_admm(ispin)%matrix, 0.0_dp)
            END DO
         END IF

         ! CPKS vector cpmos - RHS of response equation as Ax + b = 0 (sign of b)
         ! Sign is changed in linres_solver!
         ! Projector Q applied in linres_solver!
         CALL get_qs_env(qs_env, mos=mos)
         ALLOCATE (cpmos(nspins))
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
            CALL cp_fm_get_info(mo_coeff, matrix_struct=fm_struct)
            NULLIFY (cpmos(ispin)%matrix)
            CALL cp_fm_create(cpmos(ispin)%matrix, fm_struct)
            CALL cp_fm_set_all(cpmos(ispin)%matrix, 0.0_dp)
         END DO

         focc = 2.0_dp
         IF (nspins == 1) focc = 4.0_dp
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo)
            CALL cp_fm_set_all(cpmos(ispin)%matrix, 0.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(ec_env%matrix_hz(ispin)%matrix, mo_coeff, &
                                         cpmos(ispin)%matrix, nmo, &
                                         alpha=focc, beta=0.0_dp)
         END DO

         CALL response_equation(qs_env, p_env, cpmos, unit_nr)

         ! Get the response density matrix,
         ! and energy-weighted response density matrix
         DO ispin = 1, nspins
            CALL dbcsr_copy(ec_env%matrix_z(ispin)%matrix, p_env%p1(ispin)%matrix)
            CALL dbcsr_copy(ec_env%matrix_wz(ispin)%matrix, p_env%w1(ispin)%matrix)
         END DO
         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CALL dbcsr_copy(ec_env%z_admm(ispin)%matrix, p_env%p1_admm(ispin)%matrix)
         END IF

         DO ispin = 1, nspins
            CALL cp_fm_release(cpmos(ispin)%matrix)
         END DO
         DEALLOCATE (cpmos)

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

      CASE (ec_ls_solver)
         ! AO ortho solver
         CALL ec_response_ao(qs_env, &
                             ec_env%matrix_hz, &
                             ec_env%matrix_z, &
                             ec_env%matrix_wz, &
                             unit_nr, &
                             should_stop)
      CASE DEFAULT
         CPABORT("Unknown solver for response equation requested")
      END SELECT

      CALL linres_control_release(linres_control)

      CALL timestop(handle)

   END SUBROUTINE response_calculation

! **************************************************************************************************
!> \brief Parse the input section of the response solver
!> \param input Input section which controls response solver parameters
!> \param linres_control Environment for general setting of linear response calculation
!> \param unit_nr ...
!> \par History
!>       2020.05 created [Fabian Belleflamme]
!> \author Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE response_solver_write_input(input, linres_control, unit_nr)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(linres_control_type), POINTER                 :: linres_control
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=*), PARAMETER :: routineN = 'response_solver_write_input'

      INTEGER                                            :: handle, max_iter_lanczos, s_sqrt_method, &
                                                            s_sqrt_order, solver_method
      REAL(KIND=dp)                                      :: eps_lanczos

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN

         ! linres_control
         WRITE (unit_nr, '(/,T2,A)') &
            REPEAT("-", 30)//" Linear Response Solver "//REPEAT("-", 25)

         ! Which type of solver is used
         CALL section_vals_val_get(input, "METHOD", i_val=solver_method)

         SELECT CASE (solver_method)
         CASE (ec_ls_solver)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Solver: ", "AO-based CG-solver"
         CASE (ec_mo_solver)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Solver: ", "MO-based CG-solver"
         END SELECT

         WRITE (unit_nr, '(T2,A,T61,E20.3)') "eps:", linres_control%eps
         WRITE (unit_nr, '(T2,A,T61,E20.3)') "eps_filter:", linres_control%eps_filter
         WRITE (unit_nr, '(T2,A,T61,I20)') "Max iter:", linres_control%max_iter

         SELECT CASE (linres_control%preconditioner_type)
         CASE (ot_precond_full_all)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Preconditioner: ", "FULL_ALL"
         CASE (ot_precond_full_single_inverse)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Preconditioner: ", "FULL_SINGLE_INVERSE"
         CASE (ot_precond_full_single)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Preconditioner: ", "FULL_SINGLE"
         CASE (ot_precond_full_kinetic)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Preconditioner: ", "FULL_KINETIC"
         CASE (ot_precond_s_inverse)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Preconditioner: ", "FULL_S_INVERSE"
         CASE (precond_mlp)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Preconditioner: ", "MULTI_LEVEL"
         CASE (ot_precond_none)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Preconditioner: ", "NONE"
         END SELECT

         SELECT CASE (solver_method)
         CASE (ec_ls_solver)

            CALL section_vals_val_get(input, "S_SQRT_METHOD", i_val=s_sqrt_method)
            CALL section_vals_val_get(input, "S_SQRT_ORDER", i_val=s_sqrt_order)
            CALL section_vals_val_get(input, "EPS_LANCZOS", r_val=eps_lanczos)
            CALL section_vals_val_get(input, "MAX_ITER_LANCZOS", i_val=max_iter_lanczos)

            ! Response solver transforms P and KS into orthonormal basis,
            ! reuires matrx S sqrt and its inverse
            SELECT CASE (s_sqrt_method)
            CASE (ls_s_sqrt_ns)
               WRITE (unit_nr, '(T2,A,T61,A20)') "S sqrt method:", "NEWTONSCHULZ"
            CASE (ls_s_sqrt_proot)
               WRITE (unit_nr, '(T2,A,T61,A20)') "S sqrt method:", "PROOT"
            CASE DEFAULT
               CPABORT("Unknown sqrt method.")
            END SELECT
            WRITE (unit_nr, '(T2,A,T61,I20)') "S sqrt order:", s_sqrt_order

         CASE (ec_mo_solver)
         END SELECT

         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)

         CALL m_flush(unit_nr)
      END IF

      CALL timestop(handle)

   END SUBROUTINE response_solver_write_input

! **************************************************************************************************
!> \brief Initializes vectors for MO-coefficient based linear response solver
!>        and calculates response density, and energy-weighted response density matrix
!>
!> \param qs_env ...
!> \param p_env ...
!> \param cpmos RHS of equation as Ax + b = 0 (sign of b)
!> \param iounit ...
!> \param p1_mat ...
!> \param lr_section ...
! **************************************************************************************************
   SUBROUTINE response_equation(qs_env, p_env, cpmos, iounit, p1_mat, lr_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: cpmos
      INTEGER, INTENT(IN)                                :: iounit
      TYPE(cp_fm_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: p1_mat
      TYPE(section_vals_type), OPTIONAL, POINTER         :: lr_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'response_equation', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nao, nao_aux, nspins
      LOGICAL                                            :: should_stop
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: psi0, psi1
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, matrix_ks, mo_coeff, mos, psi0, psi1)

      CALL get_qs_env(qs_env, dft_control=dft_control, matrix_ks=matrix_ks, &
                      matrix_s=matrix_s, mos=mos)
      nspins = dft_control%nspins

      ! Initialize vectors:
      ! psi0 : The ground-state MO-coefficients
      ! psi1 : The "perturbed" linear response orbitals
      ALLOCATE (psi0(nspins), psi1(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
         psi0(ispin)%matrix => mo_coeff
         CALL cp_fm_get_info(mo_coeff, matrix_struct=fm_struct)
         NULLIFY (psi1(ispin)%matrix)
         CALL cp_fm_create(psi1(ispin)%matrix, fm_struct)
         CALL cp_fm_set_all(psi1(ispin)%matrix, 0.0_dp)
      END DO

      should_stop = .FALSE.
      ! The response solver
      CALL linres_solver(p_env, qs_env, psi1, cpmos, psi0, iounit, should_stop)

      ! Building the response density matrix
      DO ispin = 1, nspins
         CALL dbcsr_copy(p_env%p1(ispin)%matrix, matrix_s(1)%matrix)
      END DO
      CALL build_dm_response(psi0, psi1, p_env%p1)
      DO ispin = 1, nspins
         CALL dbcsr_scale(p_env%p1(ispin)%matrix, 0.5_dp)
      END DO
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         CPASSERT(ASSOCIATED(admm_env%work_orb_orb))
         CPASSERT(ASSOCIATED(admm_env%work_aux_orb))
         CPASSERT(ASSOCIATED(admm_env%work_aux_aux))
         nao = admm_env%nao_orb
         nao_aux = admm_env%nao_aux_fit
         DO ispin = 1, nspins
            CALL copy_dbcsr_to_fm(p_env%p1(ispin)%matrix, admm_env%work_orb_orb)
            CALL cp_gemm('N', 'N', nao_aux, nao, nao, &
                         1.0_dp, admm_env%A, admm_env%work_orb_orb, 0.0_dp, &
                         admm_env%work_aux_orb)
            CALL cp_gemm('N', 'T', nao_aux, nao_aux, nao, &
                         1.0_dp, admm_env%work_aux_orb, admm_env%A, 0.0_dp, &
                         admm_env%work_aux_aux)
            CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, p_env%p1_admm(ispin)%matrix, &
                                  keep_sparsity=.TRUE.)
         END DO
      END IF

      ! Calculate Wz = 0.5*(psi1*eps*psi0^T + psi0*eps*psi1^T)
      DO ispin = 1, nspins
         CALL calculate_wz_matrix(mos(ispin)%mo_set, psi1(ispin)%matrix, matrix_ks(ispin)%matrix, &
                                  p_env%w1(ispin)%matrix)
      END DO
      IF (PRESENT(p1_mat)) THEN
         CPASSERT(.NOT. ASSOCIATED(p1_mat))
         p1_mat => psi1
      ELSE
         DO ispin = 1, nspins
            CALL cp_fm_release(psi1(ispin)%matrix)
         ENDDO
         DEALLOCATE (psi1)
      END IF
      DEALLOCATE (psi0)
      CALL linres_control_release(linres_control)

      CALL timestop(handle)

   END SUBROUTINE response_equation

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param vh_rspace ...
!> \param vxc_rspace ...
!> \param vtau_rspace ...
!> \param vadmm_rspace ...
!> \param matrix_hz Right-hand-side of linear response equation
!> \param matrix_pz Linear response density matrix
!> \param matrix_wz Energy-weighted linear response density
!> \param zehartree Hartree volume response contribution to stress tensor
!> \param zexc XC volume response contribution to stress tensor
!> \param p_env ...
!> \param ex_env ...
! **************************************************************************************************
   SUBROUTINE response_force(qs_env, vh_rspace, vxc_rspace, vtau_rspace, vadmm_rspace, &
                             matrix_hz, matrix_pz, matrix_wz, zehartree, zexc, p_env, ex_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), POINTER                           :: vh_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rspace, vtau_rspace, vadmm_rspace
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz, matrix_pz, matrix_wz
      REAL(KIND=dp), OPTIONAL                            :: zehartree, zexc
      TYPE(qs_p_env_type), OPTIONAL, POINTER             :: p_env
      TYPE(excited_energy_type), OPTIONAL, POINTER       :: ex_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'response_force'

      INTEGER                                            :: handle, iounit, ispin, mspin, n_rep_hf, &
                                                            nao, nao_aux, natom, nder, nimages, &
                                                            nspins
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: distribute_fock_matrix, do_ex, do_hfx, &
                                                            hfx_treat_lsd_in_core, resp_only, &
                                                            s_mstruct_changed, use_virial
      REAL(KIND=dp)                                      :: eh1, ehartree, ekin_mol, eps_filter, &
                                                            eps_ppnl, exc, fconv, focc, &
                                                            total_rhoz, total_rhoz_aux
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ftot1, ftot2, ftot3
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_stress, pv_loc, stdeb, sttot, sttot2
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pd, matrix_pz, matrix_s, &
                                                            matrix_wz, mpa, scrm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p, mhd, mhx, mhy, mhz, &
                                                            mpd, mpz
      TYPE(dbcsr_type), POINTER                          :: dbwork
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_tot_gspace, rhoz_tot_gspace, &
                                                            vhxc_rspace, zv_hartree_gspace, &
                                                            zv_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_g_aux, rho_r, rho_r_aux, &
                                                            rhoz_g, rhoz_g_aux, rhoz_r, &
                                                            rhoz_r_aux, tau_r, tau_r_aux, v_xc
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(section_vals_type), POINTER                   :: hfx_section, xc_section
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xc_rho_set_type), POINTER                     :: rho_set

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_unit_nr(logger)

      do_ex = .FALSE.
      IF (PRESENT(ex_env)) do_ex = .TRUE.

      NULLIFY (ks_env, sab_orb, sac_ppl, sap_ppnl, virial)
      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      force=force, &
                      ks_env=ks_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      sab_orb=sab_orb, &
                      sac_ppl=sac_ppl, &
                      sap_ppnl=sap_ppnl, &
                      virial=virial)
      nspins = dft_control%nspins

      IF (debug_forces) THEN
         CALL get_qs_env(qs_env, natom=natom, atomic_kind_set=atomic_kind_set)
         ALLOCATE (ftot1(3, natom))
         CALL total_qs_force(ftot1, force, atomic_kind_set)
      END IF

      ! check for virial
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      fconv = 1.0E-9_dp*pascal/cell%deth
      IF (debug_stress) THEN
         sttot = virial%pv_virial
      END IF

      !     *** If LSD, then combine alpha density and beta density to
      !     *** total density: alpha <- alpha + beta   and
      matrix_pz => p_env%p1
      NULLIFY (mpa)
      IF (do_ex) THEN
         CALL dbcsr_allocate_matrix_set(mpa, nspins)
         DO ispin = 1, nspins
            ALLOCATE (mpa(ispin)%matrix)
            CALL dbcsr_create(mpa(ispin)%matrix, template=matrix_pz(ispin)%matrix)
            CALL dbcsr_copy(mpa(ispin)%matrix, matrix_pz(ispin)%matrix)
            CALL dbcsr_add(mpa(ispin)%matrix, ex_env%matrix_pe(ispin)%matrix, 1.0_dp, 1.0_dp)
            CALL dbcsr_set(matrix_hz(ispin)%matrix, 0.0_dp)
         END DO
      ELSE
         mpa => p_env%p1
      END IF

      !
      ! START OF Tr(P+Z)Hcore
      !
      IF (nspins == 2) THEN
         CALL dbcsr_add(mpa(1)%matrix, mpa(2)%matrix, 1.0_dp, 1.0_dp)
      END IF
      nimages = 1

      ! Kinetic energy matrix
      NULLIFY (scrm)
      IF (debug_forces) fodeb(1:3) = force(1)%kinetic(1:3, 1)
      IF (debug_stress) stdeb = virial%pv_ekinetic
      CALL build_kinetic_matrix(ks_env, matrix_t=scrm, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=mpa(1)%matrix)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%kinetic(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dT      ", fodeb
      END IF
      IF (debug_stress) THEN
         stdeb = fconv*(virial%pv_ekinetic - stdeb)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| Kinetic energy', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF

      IF (nspins == 2) THEN
         CALL dbcsr_add(mpa(1)%matrix, mpa(2)%matrix, 1.0_dp, -1.0_dp)
      END IF
      CALL dbcsr_deallocate_matrix_set(scrm)

      ! Initialize a matrix scrm, later used for scratch purposes
      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
      CALL dbcsr_allocate_matrix_set(scrm, nspins)
      DO ispin = 1, nspins
         ALLOCATE (scrm(ispin)%matrix)
         CALL dbcsr_create(scrm(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(scrm(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(scrm(ispin)%matrix, 0.0_dp)
      END DO

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)

      NULLIFY (cell_to_index)
      ALLOCATE (matrix_p(nspins, 1), matrix_h(nspins, 1))
      DO ispin = 1, nspins
!<<<<<<<< MERGE
!        matrix_p(ispin, 1)%matrix => matrix_pz(ispin)%matrix
!======
         matrix_p(ispin, 1)%matrix => mpa(ispin)%matrix
         matrix_h(ispin, 1)%matrix => scrm(ispin)%matrix
!>>>>>> 59c27e49f... Further work on TDDFPT forces (ADMM part)
      END DO
      matrix_h(1, 1)%matrix => scrm(1)%matrix

      IF (ASSOCIATED(sac_ppl)) THEN
         nder = 1
         IF (debug_forces) fodeb(1:3) = force(1)%gth_ppl(1:3, 1)
         IF (debug_stress) stdeb = virial%pv_ppl
         CALL build_core_ppl(matrix_h, matrix_p, force, &
                             virial, .TRUE., use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, &
                             sab_orb, sac_ppl, nimages, cell_to_index, "ORB")
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppl(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dHppl   ", fodeb
         END IF
         IF (debug_stress) THEN
            stdeb = fconv*(virial%pv_ppl - stdeb)
            CALL mp_sum(stdeb, para_env%group)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| Pz*dHppl   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF
      END IF
      eps_ppnl = dft_control%qs_control%eps_ppnl
      IF (ASSOCIATED(sap_ppnl)) THEN
         nder = 1
         IF (debug_forces) fodeb(1:3) = force(1)%gth_ppnl(1:3, 1)
         IF (debug_stress) stdeb = virial%pv_ppnl
         CALL build_core_ppnl(matrix_h, matrix_p, force, &
                              virial, .TRUE., use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, &
                              sab_orb, sap_ppnl, eps_ppnl, nimages, cell_to_index, "ORB")
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppnl(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dHppnl  ", fodeb
         END IF
         IF (debug_stress) THEN
            stdeb = fconv*(virial%pv_ppnl - stdeb)
            CALL mp_sum(stdeb, para_env%group)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| Pz*dHppnl   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF

      END IF
      ! Kim-Gordon subsystem DFT
      ! Atomic potential for nonadditive kinetic energy contribution
      IF (dft_control%qs_control%do_kg) THEN
         IF (qs_env%kg_env%tnadd_method == kg_tnadd_atomic) THEN
            CALL get_qs_env(qs_env=qs_env, kg_env=kg_env, dbcsr_dist=dbcsr_dist)

            IF (use_virial) THEN
               pv_loc = virial%pv_virial
            END IF

            IF (debug_forces) fodeb(1:3) = force(1)%kinetic(1:3, 1)
            IF (debug_stress) stdeb = virial%pv_virial
            CALL build_tnadd_mat(kg_env=kg_env, matrix_p=matrix_p, force=force, virial=virial, &
                                 calculate_forces=.TRUE., use_virial=use_virial, &
                                 qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set, &
                                 particle_set=particle_set, sab_orb=sab_orb, dbcsr_dist=dbcsr_dist)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%kinetic(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dTnadd  ", fodeb
            END IF
            IF (debug_stress) THEN
               stdeb = fconv*(virial%pv_virial - stdeb)
               CALL mp_sum(stdeb, para_env%group)
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| Pz*dTnadd   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            END IF

            ! Stress-tensor update components
            IF (use_virial) THEN
               virial%pv_ekinetic = virial%pv_ekinetic + (virial%pv_virial - pv_loc)
            END IF

         END IF
      END IF

      DEALLOCATE (matrix_p, matrix_h)
      CALL dbcsr_deallocate_matrix_set(scrm)

      ! initialize src matrix
      ! Necessary as build_kinetic_matrix will only allocate scrm(1)
      ! and not scrm(2) in open-shell case
      NULLIFY (scrm)
      CALL dbcsr_allocate_matrix_set(scrm, nspins)
      DO ispin = 1, nspins
         ALLOCATE (scrm(ispin)%matrix)
         CALL dbcsr_create(scrm(ispin)%matrix, template=matrix_pz(1)%matrix)
         CALL dbcsr_copy(scrm(ispin)%matrix, matrix_pz(ispin)%matrix)
         CALL dbcsr_set(scrm(ispin)%matrix, 0.0_dp)
      END DO

      IF (debug_forces) THEN
         ALLOCATE (ftot2(3, natom))
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1) - ftot1(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Force Pz*dHcore", fodeb
      END IF
      IF (debug_stress) THEN
         stdeb = fconv*(virial%pv_virial - sttot)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| Stress Pz*dHcore   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         ! save current total viral, does not contain volume terms yet
         sttot2 = virial%pv_virial
      END IF
      !
      ! END OF Tr(P+Z)Hcore
      !
      !
      ! Vhxc (KS potentials calculated externally)
      !
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      CALL pw_pool_create_pw(auxbas_pw_pool, vhxc_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      ! Stress-tensor: integration contribution direct term
      ! int v_Hxc[n^in]*n^z
      IF (use_virial) THEN
         pv_loc = virial%pv_virial
      END IF

      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      IF (debug_stress) stdeb = virial%pv_virial
      DO ispin = 1, nspins
         CALL pw_transfer(vh_rspace%pw, vhxc_rspace%pw)
         CALL pw_axpy(vxc_rspace(ispin)%pw, vhxc_rspace%pw)
         CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                 hmat=scrm(ispin), pmat=mpa(ispin), &
                                 qs_env=qs_env, calculate_forces=.TRUE.)
         IF (ASSOCIATED(vtau_rspace)) THEN
            CALL integrate_v_rspace(v_rspace=vtau_rspace(ispin), &
                                    hmat=scrm(ispin), pmat=mpa(ispin), &
                                    qs_env=qs_env, calculate_forces=.TRUE., compute_tau=.TRUE.)
         END IF
      END DO
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dVhxc   ", fodeb
      END IF
      IF (debug_stress) THEN
         stdeb = fconv*(virial%pv_virial - pv_loc)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| INT Pz*dVhxc   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF
      CALL pw_pool_give_back_pw(auxbas_pw_pool, vhxc_rspace%pw)

      ! Stress-tensor Pz*v_Hxc[Pin]
      IF (use_virial) THEN
         virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
      END IF

      ! KG Embedding
      ! calculate kinetic energy potential and integrate with response density
      IF (dft_control%qs_control%do_kg) THEN
         IF (qs_env%kg_env%tnadd_method == kg_tnadd_embed .OR. &
             qs_env%kg_env%tnadd_method == kg_tnadd_embed_ri) THEN

            ekin_mol = 0.0_dp
            IF (use_virial) THEN
               pv_loc = virial%pv_virial
            END IF

            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            CALL kg_ekin_subset(qs_env=qs_env, &
                                ks_matrix=scrm, &
                                ekin_mol=ekin_mol, &
                                calc_force=.TRUE., &
                                do_kernel=.FALSE., &
                                pmat_ext=mpa)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dVkg   ", fodeb
            END IF
            IF (debug_stress) THEN
               !IF (iounit > 0) WRITE(iounit, *) &
               !   "response_force | VOL 1st KG - v_KG[n_in]*n_z: ", ekin_mol
               stdeb = 1.0_dp*fconv*ekin_mol
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| VOL KG Pz*dVKG ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            END IF
            IF (debug_stress) THEN
               stdeb = fconv*(virial%pv_virial - pv_loc)
               CALL mp_sum(stdeb, para_env%group)
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| INT KG Pz*dVKG  ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            END IF
            IF (debug_stress) THEN
               stdeb = fconv*virial%pv_xc
               CALL mp_sum(stdeb, para_env%group)
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| GGA KG Pz*dVKG  ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            END IF
            IF (use_virial) THEN
               ! Direct integral contribution
               virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
            END IF

         END IF ! tnadd_method
      END IF ! do_kg

      CALL dbcsr_deallocate_matrix_set(scrm)

      !
      ! Hartree potential of response density
      !
      ALLOCATE (rhoz_r(nspins), rhoz_g(nspins))
      DO ispin = 1, nspins
         NULLIFY (rhoz_r(ispin)%pw, rhoz_g(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_r(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_g(ispin)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      END DO
      CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, zv_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, zv_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      CALL pw_zero(rhoz_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL calculate_rho_elec(ks_env=ks_env, matrix_p=mpa(ispin)%matrix, &
                                 rho=rhoz_r(ispin), rho_gspace=rhoz_g(ispin), &
                                 total_rho=total_rhoz)
         CALL pw_axpy(rhoz_g(ispin)%pw, rhoz_tot_gspace%pw)
      END DO

      ! Stress-tensor contribution second derivative
      ! Volume : int v_H[n^z]*n_in
      ! Volume : int epsilon_xc*n_z
      IF (use_virial) THEN

         NULLIFY (rho_tot_gspace%pw)
         CALL get_qs_env(qs_env, rho=rho)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

         ! Get the total input density in g-space [ions + electrons]
         CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

         h_stress(:, :) = 0.0_dp
         ! calculate associated hartree potential
         ! This term appears twice in the derivation of the equations
         ! v_H{n_in]*n_z and v_H[n_z]*n_in
         ! due to symmetry we only need to call this routine once,
         ! and count the Volume and Green function contribution
         ! which is stored in h_stress twice
         CALL pw_poisson_solve(poisson_env, &
                               density=rhoz_tot_gspace%pw, &     ! n_z
                               ehartree=ehartree, &
                               vhartree=zv_hartree_gspace%pw, &  ! v_H[n_z]
                               h_stress=h_stress, &
                               aux_density=rho_tot_gspace%pw)  ! n_in

         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)

         ! Stress tensor Green function contribution
         virial%pv_ehartree = virial%pv_ehartree + 2.0_dp*h_stress/REAL(para_env%num_pe, dp)
         virial%pv_virial = virial%pv_virial + 2.0_dp*h_stress/REAL(para_env%num_pe, dp)

         IF (debug_stress) THEN
            stdeb = -1.0_dp*fconv*ehartree
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| VOL 1st v_H[n_z]*n_in  ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            stdeb = -1.0_dp*fconv*ehartree
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| VOL 2nd v_H[n_in]*n_z  ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            stdeb = fconv*(h_stress/REAL(para_env%num_pe, dp))
            CALL mp_sum(stdeb, para_env%group)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| GREEN 1st v_H[n_z]*n_in  ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            stdeb = fconv*(h_stress/REAL(para_env%num_pe, dp))
            CALL mp_sum(stdeb, para_env%group)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| GREEN 2nd v_H[n_in]*n_z   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF

         ! Stress tensor volume term: \int v_xc[n_in]*n_z
         ! vxc_rspace already scaled, we need to unscale it!
         exc = 0.0_dp
         DO ispin = 1, nspins
            exc = exc + pw_integral_ab(rhoz_r(ispin)%pw, vxc_rspace(ispin)%pw)/ &
                  vxc_rspace(ispin)%pw%pw_grid%dvol
         END DO

         ! Add KG embedding correction
         IF (dft_control%qs_control%do_kg) THEN
            IF (qs_env%kg_env%tnadd_method == kg_tnadd_embed .OR. &
                qs_env%kg_env%tnadd_method == kg_tnadd_embed_ri) THEN
               exc = exc - ekin_mol
            END IF
         END IF

         IF (debug_stress) THEN
            stdeb = -1.0_dp*fconv*exc
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| VOL 1st eps_XC[n_in]*n_z', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF

      ELSE
         ! calculate associated hartree potential
         CALL pw_poisson_solve(poisson_env, rhoz_tot_gspace%pw, ehartree, &
                               zv_hartree_gspace%pw)
      END IF

      IF (debug_forces) fodeb(1:3) = force(1)%rho_core(1:3, 1)
      IF (debug_stress) stdeb = virial%pv_ehartree
      CALL pw_transfer(zv_hartree_gspace%pw, zv_hartree_rspace%pw)
      CALL pw_scale(zv_hartree_rspace%pw, zv_hartree_rspace%pw%pw_grid%dvol)
      ! Getting nuclear force contribution from the core charge density
      CALL integrate_v_core_rspace(zv_hartree_rspace, qs_env)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_core(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Vh(rhoz)*dncore ", fodeb
      END IF
      IF (debug_stress) THEN
         stdeb = fconv*(virial%pv_ehartree - stdeb)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| INT Vh(rhoz)*dncore   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF
      !
      CALL get_qs_env(qs_env=qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g, tau_r=tau_r)
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      END IF

      IF (use_virial) THEN
         virial%pv_xc = 0.0_dp
      END IF

      NULLIFY (v_xc)
      CALL create_kernel(qs_env, &
                         vxc=v_xc, &
                         rho=rho, &
                         rho1_r=rhoz_r, &
                         rho1_g=rhoz_g, &
                         xc_section=xc_section, &
                         compute_virial=use_virial, &
                         virial_xc=virial%pv_xc)

      ! Stress-tensor XC-kernel GGA contribution
      IF (use_virial) THEN
         virial%pv_exc = virial%pv_exc + virial%pv_xc
         virial%pv_virial = virial%pv_virial + virial%pv_xc
      END IF

      IF (debug_stress) THEN
         stdeb = 1.0_dp*fconv*virial%pv_xc
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| GGA 2nd Pin*dK*rhoz', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF

      ! Stress-tensor integral contribution of 2nd derivative terms
      IF (use_virial) THEN
         pv_loc = virial%pv_virial
      END IF

      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      IF (debug_stress) stdeb = virial%pv_virial
      DO ispin = 1, nspins
         CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(zv_hartree_rspace%pw, v_xc(ispin)%pw)
         CALL integrate_v_rspace(qs_env=qs_env, &
                                 v_rspace=v_xc(ispin), &
                                 hmat=matrix_hz(ispin), &
                                 pmat=matrix_p(ispin, 1), &
                                 calculate_forces=.TRUE.)
      END DO

      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dK*rhoz ", fodeb
      END IF
      IF (debug_stress) THEN
         stdeb = fconv*(virial%pv_virial - stdeb)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| INT 2nd f_Hxc[Pz]*Pin', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF

      ! Stress-tensor integral contribution of 2nd derivative terms
      IF (use_virial) THEN
         virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
      END IF

      ! KG Embedding
      ! calculate kinetic energy kernel, folded with response density for partial integration
      IF (dft_control%qs_control%do_kg) THEN
         IF (qs_env%kg_env%tnadd_method == kg_tnadd_embed) THEN
            ekin_mol = 0.0_dp
            IF (use_virial) THEN
               pv_loc = virial%pv_virial
            END IF

            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            IF (use_virial) virial%pv_xc = 0.0_dp
            CALL kg_ekin_subset(qs_env=qs_env, &
                                ks_matrix=matrix_hz, &
                                ekin_mol=ekin_mol, &
                                calc_force=.TRUE., &
                                do_kernel=.TRUE., &
                                pmat_ext=matrix_pz)

            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*d(Kkg)*rhoz ", fodeb
            END IF
            IF (debug_stress) THEN
               stdeb = fconv*(virial%pv_virial - pv_loc)
               CALL mp_sum(stdeb, para_env%group)
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| INT KG Pin*d(KKG)*rhoz    ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            END IF
            IF (debug_stress) THEN
               stdeb = fconv*(virial%pv_xc)
               CALL mp_sum(stdeb, para_env%group)
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| GGA KG Pin*d(KKG)*rhoz    ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            END IF

            ! Stress tensor
            IF (use_virial) THEN
               ! XC-kernel Integral contribution
               virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)

               ! XC-kernel GGA contribution
               virial%pv_exc = virial%pv_exc - virial%pv_xc
               virial%pv_virial = virial%pv_virial - virial%pv_xc
               virial%pv_xc = 0.0_dp
            END IF
         END IF
      END IF
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_tot_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, zv_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, zv_hartree_rspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_r(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_g(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      END DO
      DEALLOCATE (rhoz_r, rhoz_g, v_xc)
      IF (debug_forces) THEN
         ALLOCATE (ftot3(3, natom))
         CALL total_qs_force(ftot3, force, atomic_kind_set)
         fodeb(1:3) = ftot3(1:3, 1) - ftot2(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Force Pin*V(rhoz)", fodeb
      END IF
      CALL dbcsr_deallocate_matrix_set(scrm)

      ! Stress-tensor volume contributions
      ! These need to be applied at the end of qs_force
      IF (use_virial) THEN
         ! Adding mixed Hartree energy twice, due to symmetry
         zehartree = zehartree + 2.0_dp*ehartree
         zexc = zexc + exc
      END IF

      IF (dft_control%do_admm) THEN
         IF (dft_control%admm_control%aux_exch_func == do_admm_aux_exch_func_none) THEN
            ! nothing to do
            NULLIFY (mpz, mhx, mhy)
         ELSE
            ! add ADMM xc_section_aux terms: Pz*Vxc + P0*K0[rhoz]
            CALL get_qs_env(qs_env, admm_env=admm_env, rho_aux_fit=rho_aux_fit, &
                            matrix_s_aux_fit=scrm, task_list_aux_fit=task_list_aux_fit)
            !
            NULLIFY (mpz, mhx, mhy)
            CALL dbcsr_allocate_matrix_set(mhx, nspins, 1)
            CALL dbcsr_allocate_matrix_set(mhy, nspins, 1)
            CALL dbcsr_allocate_matrix_set(mpz, nspins, 1)
            DO ispin = 1, nspins
               ALLOCATE (mhx(ispin, 1)%matrix)
               CALL dbcsr_create(mhx(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_copy(mhx(ispin, 1)%matrix, scrm(1)%matrix)
               CALL dbcsr_set(mhx(ispin, 1)%matrix, 0.0_dp)
               ALLOCATE (mhy(ispin, 1)%matrix)
               CALL dbcsr_create(mhy(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_copy(mhy(ispin, 1)%matrix, scrm(1)%matrix)
               CALL dbcsr_set(mhy(ispin, 1)%matrix, 0.0_dp)
               ALLOCATE (mpz(ispin, 1)%matrix)
               CALL dbcsr_create(mpz(ispin, 1)%matrix, template=p_env%p1_admm(ispin)%matrix)
               CALL dbcsr_copy(mpz(ispin, 1)%matrix, p_env%p1_admm(ispin)%matrix)
               IF (do_ex) THEN
                  CALL dbcsr_add(mpz(ispin, 1)%matrix, ex_env%matrix_pe_admm(ispin)%matrix, &
                                 1.0_dp, 1.0_dp)
               END IF
            END DO
            !
            xc_section => admm_env%xc_section_aux
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            DO ispin = 1, nspins
               CALL dbcsr_set(mhz(ispin, 1)%matrix, 0.0_dp)
               CALL integrate_v_rspace(v_rspace=vadmm_rspace(ispin), &
                                       hmat=mhz(ispin, 1), pmat=mpz(ispin, 1), &
!tddft         v_rspace(ispin)%pw%cr3d = v_rspace(ispin)%pw%pw_grid%dvol*v_rspace(ispin)%pw%cr3d
!              CALL integrate_v_rspace(v_rspace=v_rspace(ispin), &
!tddft                                 hmat=mhx(ispin, 1), pmat=mpz(ispin, 1), &
                                       qs_env=qs_env, calculate_forces=.TRUE., &
                                       basis_type="AUX_FIT", task_list_external=task_list_aux_fit)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, vadmm_rspace(ispin)%pw)
            END DO
            DEALLOCATE (vadmm_rspace)
!tddft      DO ispin = 1, nspins
!              CALL pw_release(v_rspace(ispin)%pw)
!           END DO
!           DEALLOCATE (v_rspace)
!tddft      NULLIFY (rho_g_aux, rho_r_aux, tau_r_aux, rhoz_g_aux, rhoz_r_aux)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*Vxc(rho_admm)", fodeb
            END IF
            !
            CALL qs_rho_get(rho_aux_fit, rho_r=rho_r_aux, rho_g=rho_g_aux, tau_r=tau_r_aux)
            ! rhoz_aux
            NULLIFY (rhoz_g_aux, rhoz_r_aux)
            ALLOCATE (rhoz_r_aux(nspins), rhoz_g_aux(nspins))
            DO ispin = 1, nspins
               NULLIFY (rhoz_r_aux(ispin)%pw, rhoz_g_aux(ispin)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_r_aux(ispin)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_g_aux(ispin)%pw, &
                                      use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            END DO
            DO ispin = 1, nspins
               CALL calculate_rho_elec(ks_env=ks_env, matrix_p=mpz(ispin, 1)%matrix, &
                                       rho=rhoz_r_aux(ispin), rho_gspace=rhoz_g_aux(ispin), &
                                       total_rho=total_rhoz_aux, basis_type="AUX_FIT", &
                                       task_list_external=task_list_aux_fit)
            END DO
            !
            NULLIFY (deriv_set, rho_set)
            CALL xc_prep_2nd_deriv(deriv_set, rho_set, rho_r_aux, auxbas_pw_pool, xc_section=xc_section)

            CALL xc_calc_2nd_deriv(v_xc, deriv_set, rho_set, &
                                   rhoz_r_aux, rhoz_g_aux, auxbas_pw_pool, xc_section=xc_section, gapw=.FALSE.)
            CALL xc_dset_release(deriv_set)
            CALL xc_rho_set_release(rho_set)
            !
            CALL qs_rho_get(rho_aux_fit, rho_ao_kp=matrix_p)
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            DO ispin = 1, nspins
               CALL dbcsr_set(mhy(ispin, 1)%matrix, 0.0_dp)
               CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
               CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                       hmat=mhy(ispin, 1), pmat=matrix_p(ispin, 1), &
                                       calculate_forces=.TRUE., basis_type="AUX_FIT", &
                                       task_list_external=task_list_aux_fit)
            END DO
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dK*rhoz_admm ", fodeb
            END IF
            DO ispin = 1, nspins
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_r_aux(ispin)%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_g_aux(ispin)%pw)
            END DO
            DEALLOCATE (v_xc, rhoz_r_aux, rhoz_g_aux)
            !
            nao = admm_env%nao_orb
            nao_aux = admm_env%nao_aux_fit
            ALLOCATE (dbwork)
            CALL dbcsr_create(dbwork, template=matrix_hz(1)%matrix)
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(mhy(ispin, 1)%matrix, admm_env%A, &
                                            admm_env%work_aux_orb, nao)
               CALL cp_gemm('T', 'N', nao, nao, nao_aux, &
                            1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                            admm_env%work_orb_orb)
               CALL dbcsr_copy(dbwork, matrix_hz(ispin)%matrix)
               CALL dbcsr_set(dbwork, 0.0_dp)
               CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbwork, keep_sparsity=.TRUE.)
               CALL dbcsr_add(matrix_hz(ispin)%matrix, dbwork, 1.0_dp, 1.0_dp)
            END DO
            CALL dbcsr_release(dbwork)
            DEALLOCATE (dbwork)
            CALL dbcsr_deallocate_matrix_set(mpz)
         END IF
      END IF

      ! HFX
      hfx_section => section_vals_get_subs_vals(xc_section, "HF")
      CALL section_vals_get(hfx_section, explicit=do_hfx)
      IF (do_hfx) THEN
         CALL section_vals_get(hfx_section, n_repetition=n_rep_hf)
         CPASSERT(n_rep_hf == 1)
         CALL section_vals_val_get(hfx_section, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core, &
                                   i_rep_section=1)
         mspin = 1
         IF (hfx_treat_lsd_in_core) mspin = nspins
         !
         CALL get_qs_env(qs_env=qs_env, rho=rho, x_data=x_data, &
                         s_mstruct_changed=s_mstruct_changed)
         distribute_fock_matrix = .TRUE.
         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CALL get_qs_env(qs_env=qs_env, matrix_s_aux_fit=scrm)
            CALL get_qs_env(qs_env=qs_env, rho_aux_fit=rho_aux_fit)
            CALL qs_rho_get(rho_aux_fit, rho_ao_kp=matrix_p)
            NULLIFY (mpz, mhz, mpd, mhd)
            ALLOCATE (mpd(nspins, 1))
            CALL dbcsr_allocate_matrix_set(mpz, nspins, 1)
            CALL dbcsr_allocate_matrix_set(mhz, nspins, 1)
            CALL dbcsr_allocate_matrix_set(mhd, nspins, 1)
            DO ispin = 1, nspins
               ALLOCATE (mhz(ispin, 1)%matrix, mhd(ispin, 1)%matrix)
               CALL dbcsr_create(mhz(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_create(mhd(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_copy(mhz(ispin, 1)%matrix, scrm(1)%matrix)
               CALL dbcsr_copy(mhd(ispin, 1)%matrix, scrm(1)%matrix)
               CALL dbcsr_set(mhz(ispin, 1)%matrix, 0.0_dp)
               CALL dbcsr_set(mhd(ispin, 1)%matrix, 0.0_dp)
               ALLOCATE (mpz(ispin, 1)%matrix)
               CALL dbcsr_create(mpz(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_copy(mpz(ispin, 1)%matrix, p_env%p1_admm(ispin)%matrix)
               IF (do_ex) THEN
                  CALL dbcsr_add(mpz(ispin, 1)%matrix, ex_env%matrix_pe_admm(ispin)%matrix, &
                                 1.0_dp, 1.0_dp)
               END IF
               mpd(ispin, 1)%matrix => matrix_p(ispin, 1)%matrix
            END DO
            !
            DO ispin = 1, mspin
               eh1 = 0.0
               CALL integrate_four_center(qs_env, x_data, mhz, eh1, mpz, hfx_section, &
                                          para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                          ispin=ispin)
            END DO
            DO ispin = 1, mspin
               eh1 = 0.0
               CALL integrate_four_center(qs_env, x_data, mhd, eh1, mpd, hfx_section, &
                                          para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                          ispin=ispin)
            END DO
            !
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CPASSERT(ASSOCIATED(admm_env%work_aux_orb))
            CPASSERT(ASSOCIATED(admm_env%work_orb_orb))
            nao = admm_env%nao_orb
            nao_aux = admm_env%nao_aux_fit
            ALLOCATE (dbwork)
            CALL dbcsr_create(dbwork, template=matrix_hz(1)%matrix)
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(mhz(ispin, 1)%matrix, admm_env%A, &
                                            admm_env%work_aux_orb, nao)
               CALL cp_gemm('T', 'N', nao, nao, nao_aux, &
                            1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                            admm_env%work_orb_orb)
               CALL dbcsr_copy(dbwork, matrix_hz(ispin)%matrix)
               CALL dbcsr_set(dbwork, 0.0_dp)
               CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbwork, keep_sparsity=.TRUE.)
               CALL dbcsr_add(matrix_hz(ispin)%matrix, dbwork, 1.0_dp, 1.0_dp)
            END DO
            CALL dbcsr_release(dbwork)
            DEALLOCATE (dbwork)
            ! derivatives Tr (Pz [A(T)H dA/dR])
            IF (debug_forces) fodeb(1:3) = force(1)%overlap_admm(1:3, 1)
            IF (ASSOCIATED(mhx) .AND. ASSOCIATED(mhy)) THEN
               DO ispin = 1, nspins
                  CALL dbcsr_add(mhd(ispin, 1)%matrix, mhx(ispin, 1)%matrix, 1.0_dp, 1.0_dp)
                  CALL dbcsr_add(mhz(ispin, 1)%matrix, mhy(ispin, 1)%matrix, 1.0_dp, 1.0_dp)
               END DO
            END IF
            CALL qs_rho_get(rho, rho_ao=matrix_pd)
            CALL admm_projection_derivative(qs_env, admm_env, mhd, mhz, matrix_pd, mpa)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%overlap_admm(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*hfx*S' ", fodeb
            END IF
            CALL dbcsr_deallocate_matrix_set(mpz)
            CALL dbcsr_deallocate_matrix_set(mhz)
            CALL dbcsr_deallocate_matrix_set(mhd)
            CALL dbcsr_deallocate_matrix_set(mhx)
            CALL dbcsr_deallocate_matrix_set(mhy)
            DEALLOCATE (mpd)
         ELSE
            ALLOCATE (mpz(nspins, 1), mhz(nspins, 1))
            DO ispin = 1, nspins
               mhz(ispin, 1)%matrix => matrix_hz(ispin)%matrix
               mpz(ispin, 1)%matrix => mpa(ispin)%matrix
            END DO
            DO ispin = 1, mspin
               eh1 = 0.0
               CALL integrate_four_center(qs_env, x_data, mhz, eh1, mpz, hfx_section, &
                                          para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                          ispin=ispin)
            END DO
            DEALLOCATE (mhz, mpz)
         END IF

         resp_only = .TRUE.
         IF (debug_forces) fodeb(1:3) = force(1)%fock_4c(1:3, 1)
         IF (dft_control%do_admm) THEN
            CALL qs_rho_get(rho_aux_fit, rho_ao_kp=matrix_p)
            NULLIFY (matrix_pz)
            CALL dbcsr_allocate_matrix_set(matrix_pz, nspins)
            DO ispin = 1, nspins
               ALLOCATE (matrix_pz(ispin)%matrix)
               CALL dbcsr_create(matrix_pz(ispin)%matrix, template=p_env%p1_admm(ispin)%matrix)
               CALL dbcsr_copy(matrix_pz(ispin)%matrix, p_env%p1_admm(ispin)%matrix)
               IF (do_ex) THEN
                  CALL dbcsr_add(matrix_pz(ispin)%matrix, ex_env%matrix_pe_admm(ispin)%matrix, &
                                 1.0_dp, 1.0_dp)
               END IF
            END DO
            CALL derivatives_four_center(qs_env, matrix_p, matrix_pz, hfx_section, para_env, &
                                         1, use_virial, resp_only=resp_only)
            CALL dbcsr_deallocate_matrix_set(matrix_pz)
         ELSE
            CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
            CALL derivatives_four_center(qs_env, matrix_p, mpa, hfx_section, para_env, &
                                         1, use_virial, resp_only=resp_only)
         END IF
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%fock_4c(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*hfx ", fodeb
         END IF
!<<<<<<< MERGE
!         IF (dft_control%do_admm) THEN
!            CALL dbcsr_deallocate_matrix_set(mhz)
!            DEALLOCATE (mpz)
!         ELSE
!            DEALLOCATE (mpz, mhz)
!         END IF
!         matrix_pz => p_env%p1
!=======
!>>>>>>> 59c27e49f... Further work on TDDFPT forces (ADMM part)
      END IF

      ! Overlap matrix
      ! H(drho+dz) + Wz
      ! If ground-state density matrix solved by diagonalization, then use this
      IF (.NOT. dft_control%qs_control%do_ls_scf .AND. nspins == 1) THEN
         matrix_wz => p_env%w1
         focc = 1.0_dp
         IF (nspins == 1) focc = 2.0_dp
         CALL get_qs_env(qs_env, mos=mos)
         DO ispin = 1, nspins
            CALL calculate_whz_matrix(mos(ispin)%mo_set%mo_coeff, matrix_hz(ispin)%matrix, &
                                      matrix_wz(ispin)%matrix, focc)
         END DO
      ELSE
         ! Ground-state density has been calculated by LS
         eps_filter = dft_control%qs_control%eps_filter_matrix
         CALL calculate_whz_ao_matrix(qs_env, matrix_hz, matrix_wz, eps_filter)
      END IF
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_wz(1)%matrix, matrix_wz(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF

      IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
      IF (debug_stress) stdeb = virial%pv_overlap
      NULLIFY (scrm)
      CALL build_overlap_matrix(ks_env, matrix_s=scrm, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", basis_type_b="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=matrix_wz(1)%matrix)

      IF (SIZE(matrix_wz, 1) == 2) THEN
         CALL dbcsr_add(matrix_wz(1)%matrix, matrix_wz(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      END IF

      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%overlap(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Wz*dS ", fodeb
      END IF
      IF (debug_stress) THEN
         stdeb = fconv*(virial%pv_overlap - stdeb)
         CALL mp_sum(stdeb, para_env%group)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| WHz   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF
      CALL dbcsr_deallocate_matrix_set(scrm)

      IF (debug_forces) THEN
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1) - ftot1(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Response Force", fodeb
         DEALLOCATE (ftot1, ftot2, ftot3)
      END IF

      IF (do_ex) THEN
         CALL dbcsr_deallocate_matrix_set(mpa)
      END IF

      CALL timestop(handle)

   END SUBROUTINE response_force

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param matrix_hz ...
!> \param ex_env ...
! **************************************************************************************************
   SUBROUTINE response_force_xtb(qs_env, p_env, matrix_hz, ex_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz
      TYPE(excited_energy_type), OPTIONAL, POINTER       :: ex_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'response_force_xtb'

      INTEGER                                            :: atom_a, handle, iatom, ikind, iounit, &
                                                            is, ispin, na, natom, natorb, nimages, &
                                                            nkind, ns, nsgf, nspins
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(5)                              :: occ
      LOGICAL                                            :: do_ex, use_virial
      REAL(KIND=dp)                                      :: focc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mcharge, mcharge1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aocg, aocg1, charges, charges1, ftot1, &
                                                            ftot2
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pz, matrix_wz, mpa, p_matrix, scrm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dbcsr_type), POINTER                          :: s_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_unit_nr(logger)

      do_ex = .FALSE.
      IF (PRESENT(ex_env)) do_ex = .TRUE.

      NULLIFY (ks_env, sab_orb)
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, dft_control=dft_control, &
                      sab_orb=sab_orb)
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, force=force)
      nspins = dft_control%nspins

      IF (debug_forces) THEN
         CALL get_qs_env(qs_env, natom=natom, atomic_kind_set=atomic_kind_set)
         ALLOCATE (ftot1(3, natom))
         ALLOCATE (ftot2(3, natom))
         CALL total_qs_force(ftot1, force, atomic_kind_set)
      END IF

      matrix_pz => p_env%p1
      NULLIFY (mpa)
      IF (do_ex) THEN
         CALL dbcsr_allocate_matrix_set(mpa, nspins)
         DO ispin = 1, nspins
            ALLOCATE (mpa(ispin)%matrix)
            CALL dbcsr_create(mpa(ispin)%matrix, template=matrix_pz(ispin)%matrix)
            CALL dbcsr_copy(mpa(ispin)%matrix, matrix_pz(ispin)%matrix)
            CALL dbcsr_add(mpa(ispin)%matrix, ex_env%matrix_pe(ispin)%matrix, 1.0_dp, 1.0_dp)
            CALL dbcsr_set(matrix_hz(ispin)%matrix, 0.0_dp)
         END DO
      ELSE
         mpa => p_env%p1
      END IF
      !
      ! START OF Tr(P+Z)Hcore
      !
      IF (nspins == 2) THEN
         CALL dbcsr_add(mpa(1)%matrix, mpa(2)%matrix, 1.0_dp, 1.0_dp)
      END IF
      ! Hcore  matrix
      IF (debug_forces) fodeb(1:3) = force(1)%all_potential(1:3, 1)
      CALL xtb_hab_force(qs_env, mpa(1)%matrix)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%all_potential(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dHcore  ", fodeb
      END IF
      IF (nspins == 2) THEN
         CALL dbcsr_add(mpa(1)%matrix, mpa(2)%matrix, 1.0_dp, -1.0_dp)
      END IF
      !
      ! END OF Tr(P+Z)Hcore
      !
      use_virial = .FALSE.
      nimages = 1
      !
      ! Hartree potential of response density
      !
      IF (dft_control%qs_control%xtb_control%coulomb_interaction) THEN
         ! Mulliken charges
         CALL get_qs_env(qs_env, rho=rho, particle_set=particle_set, matrix_s_kp=matrix_s)
         natom = SIZE(particle_set)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
         ALLOCATE (mcharge(natom), charges(natom, 5))
         ALLOCATE (mcharge1(natom), charges1(natom, 5))
         charges = 0.0_dp
         charges1 = 0.0_dp
         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
         nkind = SIZE(atomic_kind_set)
         CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)
         ALLOCATE (aocg(nsgf, natom))
         aocg = 0.0_dp
         ALLOCATE (aocg1(nsgf, natom))
         aocg1 = 0.0_dp
         p_matrix => matrix_p(:, 1)
         s_matrix => matrix_s(1, 1)%matrix
         CALL ao_charges(p_matrix, s_matrix, aocg, para_env)
         CALL ao_charges(mpa, s_matrix, aocg1, para_env)
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=na)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, occupation=occ)
            DO iatom = 1, na
               atom_a = atomic_kind_set(ikind)%atom_list(iatom)
               charges(atom_a, :) = REAL(occ(:), KIND=dp)
               DO is = 1, natorb
                  ns = lao(is) + 1
                  charges(atom_a, ns) = charges(atom_a, ns) - aocg(is, atom_a)
                  charges1(atom_a, ns) = charges1(atom_a, ns) - aocg1(is, atom_a)
               END DO
            END DO
         END DO
         DEALLOCATE (aocg, aocg1)
         DO iatom = 1, natom
            mcharge(iatom) = SUM(charges(iatom, :))
            mcharge1(iatom) = SUM(charges1(iatom, :))
         END DO
         ! Coulomb Kernel
         CALL xtb_coulomb_hessian(qs_env, matrix_hz, charges1, mcharge1, mcharge)
         CALL calc_xtb_ehess_force(qs_env, p_matrix, mpa, charges, mcharge, charges1, &
                                   mcharge1, debug_forces)
         !
         DEALLOCATE (charges, mcharge, charges1, mcharge1)
      END IF
      ! Overlap matrix
      ! H(drho+dz) + Wz
      matrix_wz => p_env%w1
      focc = 0.5_dp
      IF (nspins == 1) focc = 1.0_dp
      CALL get_qs_env(qs_env, mos=mos)
      DO ispin = 1, nspins
         CALL calculate_whz_matrix(mos(ispin)%mo_set%mo_coeff, matrix_hz(ispin)%matrix, &
                                   matrix_wz(ispin)%matrix, focc)
      END DO
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_wz(1)%matrix, matrix_wz(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF
      IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
      NULLIFY (scrm)
      CALL build_overlap_matrix(ks_env, matrix_s=scrm, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", basis_type_b="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=matrix_wz(1)%matrix)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%overlap(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Wz*dS ", fodeb
      END IF
      CALL dbcsr_deallocate_matrix_set(scrm)

      IF (debug_forces) THEN
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1) - ftot1(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Response Force", fodeb
         DEALLOCATE (ftot1, ftot2)
      END IF

      IF (do_ex) THEN
         CALL dbcsr_deallocate_matrix_set(mpa)
      END IF

      CALL timestop(handle)

   END SUBROUTINE response_force_xtb

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param admm_env ...
!> \param matrix_hd ...
!> \param matrix_hz ...
!> \param matrix_pd ...
!> \param matrix_pz ...
! **************************************************************************************************
   SUBROUTINE admm_projection_derivative(qs_env, admm_env, matrix_hd, matrix_hz, matrix_pd, matrix_pz)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_hd, matrix_hz
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pd, matrix_pz

      CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_projection_derivative', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nao, natom, naux, nspins
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: admm_force
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s_aux_fit, matrix_s_aux_fit_vs_orb
      TYPE(dbcsr_type), POINTER                          :: matrix_w_q, matrix_w_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux_fit_asymm, sab_aux_fit_vs_orb
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(admm_env))
      CPASSERT(ASSOCIATED(matrix_hd))
      CPASSERT(ASSOCIATED(matrix_pd))
      CPASSERT(ASSOCIATED(matrix_hz))
      CPASSERT(ASSOCIATED(matrix_pz))

      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      matrix_s_aux_fit=matrix_s_aux_fit, &
                      matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb, &
                      sab_aux_fit_vs_orb=sab_aux_fit_vs_orb, &
                      sab_aux_fit_asymm=sab_aux_fit_asymm)

      ALLOCATE (matrix_w_q)
      CALL dbcsr_copy(matrix_w_q, matrix_s_aux_fit_vs_orb(1)%matrix, &
                      "W MATRIX AUX Q")
      CALL cp_dbcsr_alloc_block_from_nbl(matrix_w_q, sab_aux_fit_vs_orb)
      ALLOCATE (matrix_w_s)
      CALL dbcsr_create(matrix_w_s, template=matrix_s_aux_fit(1)%matrix, &
                        name='W MATRIX AUX S', &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(matrix_w_s, sab_aux_fit_asymm)

      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                      natom=natom, force=force)
      ALLOCATE (admm_force(3, natom))

      nspins = SIZE(matrix_pz)
      nao = admm_env%nao_orb
      naux = admm_env%nao_aux_fit
      DO ispin = 1, nspins
         CALL copy_dbcsr_to_fm(matrix_hd(ispin, 1)%matrix, admm_env%work_aux_aux)
         CALL cp_gemm("N", "N", naux, naux, naux, 1.0_dp, admm_env%s_inv, &
                      admm_env%work_aux_aux, 0.0_dp, admm_env%work_aux_aux2)
         CALL cp_gemm("N", "N", naux, nao, naux, 1.0_dp, admm_env%work_aux_aux2, &
                      admm_env%A, 0.0_dp, admm_env%work_aux_orb)
         CALL copy_dbcsr_to_fm(matrix_pz(ispin)%matrix, admm_env%work_orb_orb)
         ! admm_env%work_aux_orb2 = S-1*H(D)*A*PZ
         CALL cp_gemm("N", "N", naux, nao, nao, 1.0_dp, admm_env%work_aux_orb, &
                      admm_env%work_orb_orb, 0.0_dp, admm_env%work_aux_orb2)

         CALL copy_dbcsr_to_fm(matrix_hz(ispin, 1)%matrix, admm_env%work_aux_aux)
         CALL cp_gemm("N", "N", naux, naux, naux, 1.0_dp, admm_env%s_inv, &
                      admm_env%work_aux_aux, 0.0_dp, admm_env%work_aux_aux2)
         CALL cp_gemm("N", "N", naux, nao, naux, 1.0_dp, admm_env%work_aux_aux2, &
                      admm_env%A, 0.0_dp, admm_env%work_aux_orb)
         CALL copy_dbcsr_to_fm(matrix_pd(ispin)%matrix, admm_env%work_orb_orb)
         ! admm_env%work_aux_orb2 = S-1*H(Z)*A*PD
         CALL cp_gemm("N", "N", naux, nao, nao, 1.0_dp, admm_env%work_aux_orb, &
                      admm_env%work_orb_orb, 1.0_dp, admm_env%work_aux_orb2)

         ! admm_env%work_aux_aux = [S-1*H(D)*A*PZ+S-1*H(Z)*A*PD]*A(T)
         CALL cp_gemm("N", "T", naux, naux, nao, 1.0_dp, admm_env%work_aux_orb2, &
                      admm_env%A, 0.0_dp, admm_env%work_aux_aux)
         !
         CALL copy_fm_to_dbcsr(admm_env%work_aux_orb2, matrix_w_q, keep_sparsity=.TRUE.)
         CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, matrix_w_s, keep_sparsity=.TRUE.)
! focc = 2?
         CALL dbcsr_scale(matrix_w_q, -2.0_dp)
         CALL dbcsr_scale(matrix_w_s, 2.0_dp)
         !
         admm_force = 0.0_dp
         CALL build_overlap_force(ks_env, admm_force, &
                                  basis_type_a="AUX_FIT", basis_type_b="AUX_FIT", &
                                  sab_nl=sab_aux_fit_asymm, matrix_p=matrix_w_s)
         CALL build_overlap_force(ks_env, admm_force, &
                                  basis_type_a="AUX_FIT", basis_type_b="ORB", &
                                  sab_nl=sab_aux_fit_vs_orb, matrix_p=matrix_w_q)
         ! add forces
         CALL add_qs_force(admm_force, force, "overlap_admm", atomic_kind_set)
      END DO

      DEALLOCATE (admm_force)
      CALL dbcsr_deallocate_matrix(matrix_w_s)
      CALL dbcsr_deallocate_matrix(matrix_w_q)

      CALL timestop(handle)

   END SUBROUTINE admm_projection_derivative

! **************************************************************************************************
!> \brief calculate the Kohn-Sham reference potential
!> \param qs_env ...
!> \param vh_rspace ...
!> \param vxc_rspace ...
!> \param vtau_rspace ...
!> \param vadmm_rspace ...
!> \param ehartree ...
!> \param exc ...
!> \param h_stress container for the stress tensor of the Hartree term
!> \par History
!>      10.2019 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ks_ref_potential(qs_env, vh_rspace, vxc_rspace, vtau_rspace, vadmm_rspace, ehartree, exc, h_stress)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), POINTER                           :: vh_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rspace, vtau_rspace, vadmm_rspace
      REAL(KIND=dp), INTENT(OUT)                         :: ehartree, exc
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: h_stress

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ks_ref_potential'

      INTEGER                                            :: handle, iab, ispin, nspins
      REAL(dp)                                           :: eadmm, eovrl, eself
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_xc
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_admm_rspace, v_admm_tau_rspace, &
                                                            v_rspace, v_tau_rspace
      TYPE(pw_p_type), POINTER                           :: rho_core
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      ! get all information on the electronic density
      NULLIFY (rho, ks_env)
      CALL get_qs_env(qs_env=qs_env, rho=rho, dft_control=dft_control, &
                      para_env=para_env, blacs_env=blacs_env, ks_env=ks_env, rho_core=rho_core)

      nspins = dft_control%nspins

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
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      ! Get the total density in g-space [ions + electrons]
      CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree, &
                            v_hartree_gspace%pw, h_stress=h_stress, rho_core=rho_core)
      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
      !
      CALL calculate_ecore_self(qs_env, E_self_core=eself)
      CALL calculate_ecore_overlap(qs_env, para_env, PRESENT(h_stress), E_overlap_core=eovrl)
      ehartree = ehartree + eovrl + eself

      ! v_rspace and v_tau_rspace are generated from the auxbas pool
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      END IF
      NULLIFY (v_rspace, v_tau_rspace)
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=xc_section, &
                         vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE.)

      NULLIFY (v_admm_rspace, v_admm_tau_rspace)
      IF (dft_control%do_admm) THEN
         IF (dft_control%admm_control%aux_exch_func /= do_admm_aux_exch_func_none) THEN
            ! For the virial, we have to save the pv_xc component because it will be reset in qs_vxc_create
            IF (PRESENT(h_stress)) THEN
               CALL get_qs_env(qs_env, virial=virial)
               virial_xc = virial%pv_xc
            END IF
            CALL get_qs_env(qs_env, rho_aux_fit=rho)
            xc_section => admm_env%xc_section_aux
            CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=xc_section, &
                               vxc_rho=v_admm_rspace, vxc_tau=v_admm_tau_rspace, exc=eadmm, just_energy=.FALSE.)
            IF (PRESENT(h_stress)) virial%pv_xc = virial%pv_xc + virial_xc
         END IF
      END IF

      ! allocate potentials
      IF (ASSOCIATED(vh_rspace)) THEN
         CALL pw_release(vh_rspace%pw)
      ELSE
         ALLOCATE (vh_rspace)
         NULLIFY (vh_rspace%pw)
      END IF
      IF (ASSOCIATED(vxc_rspace)) THEN
         DO iab = 1, SIZE(vxc_rspace)
            CALL pw_release(vxc_rspace(iab)%pw)
         END DO
      ELSE
         ALLOCATE (vxc_rspace(nspins))
         DO iab = 1, nspins
            NULLIFY (vxc_rspace(iab)%pw)
         END DO
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
         IF (ASSOCIATED(vtau_rspace)) THEN
            DO iab = 1, SIZE(vtau_rspace)
               CALL pw_release(vtau_rspace(iab)%pw)
            END DO
         ELSE
            ALLOCATE (vtau_rspace(nspins))
            DO iab = 1, nspins
               NULLIFY (vtau_rspace(iab)%pw)
            END DO
         END IF
      ELSE
         NULLIFY (vtau_rspace)
      END IF
      IF (ASSOCIATED(v_admm_rspace)) THEN
         IF (ASSOCIATED(vadmm_rspace)) THEN
            DO iab = 1, SIZE(vadmm_rspace)
               CALL pw_release(vadmm_rspace(iab)%pw)
            END DO
         ELSE
            ALLOCATE (vadmm_rspace(nspins))
            DO iab = 1, nspins
               NULLIFY (vadmm_rspace(iab)%pw)
            END DO
         END IF
      ELSE
         NULLIFY (vadmm_rspace)
      END IF

      pw_grid => v_hartree_rspace%pw%pw_grid
      CALL pw_create(vh_rspace%pw, pw_grid, use_data=REALDATA3D, in_space=REALSPACE)
      DO ispin = 1, nspins
         NULLIFY (vxc_rspace(ispin)%pw)
         CALL pw_create(vxc_rspace(ispin)%pw, pw_grid, &
                        use_data=REALDATA3D, in_space=REALSPACE)
         IF (ASSOCIATED(vtau_rspace)) THEN
            NULLIFY (vtau_rspace(ispin)%pw)
            CALL pw_create(vtau_rspace(ispin)%pw, pw_grid, &
                           use_data=REALDATA3D, in_space=REALSPACE)
         END IF
         IF (ASSOCIATED(vadmm_rspace)) THEN
            NULLIFY (vadmm_rspace(ispin)%pw)
            CALL pw_create(vadmm_rspace(ispin)%pw, pw_grid, &
                           use_data=REALDATA3D, in_space=REALSPACE)
         END IF
      END DO
      !
      CALL pw_transfer(v_hartree_rspace%pw, vh_rspace%pw)
      IF (ASSOCIATED(v_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_transfer(v_rspace(ispin)%pw, vxc_rspace(ispin)%pw)
            CALL pw_scale(vxc_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
            IF (ASSOCIATED(v_tau_rspace)) THEN
               CALL pw_transfer(v_tau_rspace(ispin)%pw, vtau_rspace(ispin)%pw)
               CALL pw_scale(vtau_rspace(ispin)%pw, v_tau_rspace(ispin)%pw%pw_grid%dvol)
            END IF
         END DO
      ELSE
         DO ispin = 1, nspins
            CALL pw_zero(vxc_rspace(ispin)%pw)
         END DO
      END IF
      IF (ASSOCIATED(v_admm_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_transfer(v_admm_rspace(ispin)%pw, vadmm_rspace(ispin)%pw)
            CALL pw_scale(vadmm_rspace(ispin)%pw, vadmm_rspace(ispin)%pw%pw_grid%dvol)
         END DO
      END IF

      ! return pw grids
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      IF (ASSOCIATED(v_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
            IF (ASSOCIATED(v_tau_rspace)) THEN
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
            END IF
         END DO
         DEALLOCATE (v_rspace)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) DEALLOCATE (v_tau_rspace)
      IF (ASSOCIATED(v_admm_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_admm_rspace(ispin)%pw)
         END DO
         DEALLOCATE (v_admm_rspace)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ks_ref_potential

! **************************************************************************************************
!> \brief Win = focc*(P*(H[P_out - P_in] + H[Z] )*P)
!>        Langrange multiplier matrix with response and perturbation (Harris) kernel matrices
!>
!> \param qs_env ...
!> \param matrix_hz ...
!> \param matrix_whz ...
!> \param eps_filter ...
!> \param
!> \par History
!>       2020.2 created [Fabian Belleflamme]
!> \author Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE calculate_whz_ao_matrix(qs_env, matrix_hz, matrix_whz, eps_filter)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: matrix_hz
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: matrix_whz
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_whz_ao_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: scaling
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_type)                                   :: matrix_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(matrix_hz))
      CPASSERT(ASSOCIATED(matrix_whz))

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      rho=rho, &
                      para_env=para_env)
      nspins = dft_control%nspins
      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ! init temp matrix
      CALL dbcsr_create(matrix_tmp, template=matrix_hz(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      !Spin factors simplify to
      scaling = 1.0_dp
      IF (nspins == 1) scaling = 0.5_dp

      ! Operation in MO-solver :
      ! Whz = focc*(CC^T*Hz*CC^T)
      ! focc = 2.0_dp Closed-shell
      ! focc = 1.0_dp Open-shell

      ! Operation in AO-solver :
      ! Whz = (scaling*P)*(focc*Hz)*(scaling*P)
      ! focc see above
      ! scaling = 0.5_dp Closed-shell (P = 2*CC^T), WHz = (0.5*P)*(2*Hz)*(0.5*P)
      ! scaling = 1.0_dp Open-shell, WHz = P*Hz*P

      ! Spin factors from Hz and P simplify to
      scaling = 1.0_dp
      IF (nspins == 1) scaling = 0.5_dp

      DO ispin = 1, nspins

         ! tmp = H*CC^T
         CALL dbcsr_multiply("N", "N", scaling, matrix_hz(ispin)%matrix, rho_ao(ispin)%matrix, &
                             0.0_dp, matrix_tmp, filter_eps=eps_filter)
         ! WHz = CC^T*tmp
         ! WHz = Wz + (scaling*P)*(focc*Hz)*(scaling*P)
         ! WHz = Wz + scaling*(P*Hz*P)
         CALL dbcsr_multiply("N", "N", 1.0_dp, rho_ao(ispin)%matrix, matrix_tmp, &
                             1.0_dp, matrix_whz(ispin)%matrix, filter_eps=eps_filter, &
                             retain_sparsity=.TRUE.)

      END DO

      CALL dbcsr_release(matrix_tmp)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************

END MODULE response_solver
