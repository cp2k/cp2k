!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
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
   USE admm_methods,                    ONLY: admm_projection_derivative
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type
   USE core_ae,                         ONLY: build_core_ae
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
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_init_random,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE ec_env_types,                    ONLY: energy_correction_type
   USE ec_methods,                      ONLY: create_kernel,&
                                              ec_mos_init
   USE ec_orth_solver,                  ONLY: ec_response_ao
   USE exstates_types,                  ONLY: excited_energy_type
   USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals,&
                                              init_coulomb_local
   USE hartree_local_types,             ONLY: hartree_local_create,&
                                              hartree_local_release,&
                                              hartree_local_type
   USE hfx_derivatives,                 ONLY: derivatives_four_center
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_ri,                          ONLY: hfx_ri_update_forces,&
                                              hfx_ri_update_ks
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: &
        do_admm_aux_exch_func_none, ec_ls_solver, ec_mo_solver, kg_tnadd_atomic, kg_tnadd_embed, &
        kg_tnadd_embed_ri, ls_s_sqrt_ns, ls_s_sqrt_proot, ot_precond_full_all, &
        ot_precond_full_kinetic, ot_precond_full_single, ot_precond_full_single_inverse, &
        ot_precond_none, ot_precond_s_inverse, precond_mlp, xc_none
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kg_correction,                   ONLY: kg_ekin_subset
   USE kg_environment_types,            ONLY: kg_environment_type
   USE kg_tnadd_mat,                    ONLY: build_tnadd_mat
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush
   USE mathlib,                         ONLY: det_3x3
   USE message_passing,                 ONLY: mp_para_env_type
   USE mulliken,                        ONLY: ao_charges
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: pascal
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
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE qs_2nd_kernel_ao,                ONLY: build_dm_response
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_density_matrices,             ONLY: calculate_whz_matrix,&
                                              calculate_wz_matrix
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: qs_force_type,&
                                              total_qs_force
   USE qs_gapw_densities,               ONLY: prepare_gapw_den
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_atom,                      ONLY: update_ks_atom
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_methods,               ONLY: linres_solver
   USE qs_linres_types,                 ONLY: linres_control_type
   USE qs_local_rho_types,              ONLY: local_rho_set_create,&
                                              local_rho_set_release,&
                                              local_rho_type
   USE qs_matrix_pools,                 ONLY: mpools_rebuild_fm_pools
   USE qs_mo_methods,                   ONLY: make_basis_sm
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_oce_types,                    ONLY: oce_matrix_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_p_env_methods,                ONLY: p_env_create,&
                                              p_env_psi0_changed
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace,&
                                              rho0_s_grid_create
   USE qs_rho0_methods,                 ONLY: init_rho0
   USE qs_rho0_types,                   ONLY: get_rho0_mpole
   USE qs_rho_atom_methods,             ONLY: allocate_rho_atom_internals,&
                                              calculate_rho_atom_coeff
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc_atom,                     ONLY: calculate_vxc_atom,&
                                              calculate_xc_2nd_deriv_atom
   USE task_list_types,                 ONLY: task_list_type
   USE virial_methods,                  ONLY: one_third_sum_diag
   USE virial_types,                    ONLY: virial_type
   USE xtb_ehess,                       ONLY: xtb_coulomb_hessian
   USE xtb_ehess_force,                 ONLY: calc_xtb_ehess_force
   USE xtb_matrices,                    ONLY: xtb_hab_force
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'response_solver'

   PUBLIC :: response_calculation, response_equation, response_force, response_force_xtb, &
             response_equation_new

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Initializes solver of linear response equation for energy correction
!> \brief Call AO or MO based linear response solver for energy correction
!>
!> \param qs_env The quickstep environment
!> \param ec_env The energy correction environment
!>
!> \date    01.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE response_calculation(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'response_calculation'

      INTEGER                                            :: handle, homo, ispin, nao, nao_aux, nmo, &
                                                            nocc, nspins, solver_method, unit_nr
      LOGICAL                                            :: should_stop
      REAL(KIND=dp)                                      :: focc
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: sv
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: cpmos, mo_occ
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_s_aux, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input, solver_section

      CALL timeset(routineN, handle)

      NULLIFY (admm_env, dft_control, energy, logger, matrix_s, matrix_s_aux, mo_coeff, mos, para_env, &
               rho_ao, sab_orb, solver_section)

      ! Get useful output unit
      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
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
      ALLOCATE (linres_control)
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

      ! Allocate and initialize response density matrix Z,
      ! and the energy weighted response density matrix
      ! Template is the ground-state overlap matrix
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

      ! MO solver requires MO's of the ground-state calculation,
      ! The MOs environment is not allocated if LS-DFT has been used.
      ! Introduce MOs here
      ! Remark: MOS environment also required for creation of p_env
      IF (dft_control%qs_control%do_ls_scf) THEN

         ! Allocate and initialize MO environment
         CALL ec_mos_init(qs_env, matrix_s(1)%matrix)
         CALL get_qs_env(qs_env, mos=mos, rho=rho)

         ! Get ground-state density matrix
         CALL qs_rho_get(rho, rho_ao=rho_ao)

         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin), &
                            mo_coeff=mo_coeff, &
                            nmo=nmo, nao=nao, homo=homo)

            CALL cp_fm_set_all(mo_coeff, 0.0_dp)
            CALL cp_fm_init_random(mo_coeff, nmo)

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
            NULLIFY (blacs_env)
            CALL get_qs_env(qs_env, blacs_env=blacs_env)
            CALL mpools_rebuild_fm_pools(qs_env%mpools, mos=mos, &
                                         blacs_env=blacs_env, para_env=para_env)
         END DO
      END IF

      ! initialize p_env
      ! Remark: mos environment is needed for this
      IF (ASSOCIATED(ec_env%p_env)) THEN
         CALL p_env_release(ec_env%p_env)
         DEALLOCATE (ec_env%p_env)
         NULLIFY (ec_env%p_env)
      END IF
      ALLOCATE (ec_env%p_env)
      CALL p_env_create(ec_env%p_env, qs_env, orthogonal_orbitals=.TRUE., &
                        linres_control=linres_control)
      CALL p_env_psi0_changed(ec_env%p_env, qs_env)
      ! Total energy overwritten, replace with Etot from energy correction
      CALL get_qs_env(qs_env, energy=energy)
      energy%total = ec_env%etotal
      !
      p_env => ec_env%p_env
      !
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
         CALL get_admm_env(qs_env%admm_env, matrix_s_aux_fit=matrix_s_aux)
         CALL dbcsr_allocate_matrix_set(p_env%p1_admm, nspins)
         DO ispin = 1, nspins
            ALLOCATE (p_env%p1_admm(ispin)%matrix)
            CALL dbcsr_create(p_env%p1_admm(ispin)%matrix, &
                              template=matrix_s_aux(1)%matrix)
            CALL dbcsr_copy(p_env%p1_admm(ispin)%matrix, matrix_s_aux(1)%matrix)
            CALL dbcsr_set(p_env%p1_admm(ispin)%matrix, 0.0_dp)
         END DO
      END IF

      ! Choose between MO-solver and AO-solver
      SELECT CASE (solver_method)
      CASE (ec_mo_solver)

         ! CPKS vector cpmos - RHS of response equation as Ax + b = 0 (sign of b)
         ! Sign is changed in linres_solver!
         ! Projector Q applied in linres_solver!
         CALL get_qs_env(qs_env, mos=mos)
         ALLOCATE (cpmos(nspins), mo_occ(nspins))
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, homo=nocc)
            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, ncol_global=nocc, &
                                     template_fmstruct=mo_coeff%matrix_struct)
            CALL cp_fm_create(cpmos(ispin), fm_struct)
            CALL cp_fm_set_all(cpmos(ispin), 0.0_dp)
            CALL cp_fm_create(mo_occ(ispin), fm_struct)
            CALL cp_fm_to_fm(mo_coeff, mo_occ(ispin), nocc)
            CALL cp_fm_struct_release(fm_struct)
         END DO

         focc = 2.0_dp
         IF (nspins == 1) focc = 4.0_dp
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, homo=nocc)
            CALL cp_dbcsr_sm_fm_multiply(ec_env%matrix_hz(ispin)%matrix, mo_occ(ispin), &
                                         cpmos(ispin), nocc, &
                                         alpha=focc, beta=0.0_dp)
         END DO
         CALL cp_fm_release(mo_occ)

         CALL response_equation_new(qs_env, p_env, cpmos, unit_nr)

         ! Get the response density matrix,
         ! and energy-weighted response density matrix
         DO ispin = 1, nspins
            CALL dbcsr_copy(ec_env%matrix_z(ispin)%matrix, p_env%p1(ispin)%matrix)
            CALL dbcsr_copy(ec_env%matrix_wz(ispin)%matrix, p_env%w1(ispin)%matrix)
         END DO
         CALL cp_fm_release(cpmos)

      CASE (ec_ls_solver)

         ! AO ortho solver
         CALL ec_response_ao(qs_env=qs_env, &
                             p_env=p_env, &
                             matrix_hz=ec_env%matrix_hz, &
                             matrix_pz=ec_env%matrix_z, &
                             matrix_wz=ec_env%matrix_wz, &
                             iounit=unit_nr, &
                             should_stop=should_stop)

         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CPASSERT(ASSOCIATED(admm_env%work_orb_orb))
            CPASSERT(ASSOCIATED(admm_env%work_aux_orb))
            CPASSERT(ASSOCIATED(admm_env%work_aux_aux))
            nao = admm_env%nao_orb
            nao_aux = admm_env%nao_aux_fit
            DO ispin = 1, nspins
               CALL copy_dbcsr_to_fm(ec_env%matrix_z(ispin)%matrix, admm_env%work_orb_orb)
               CALL parallel_gemm('N', 'N', nao_aux, nao, nao, &
                                  1.0_dp, admm_env%A, admm_env%work_orb_orb, 0.0_dp, &
                                  admm_env%work_aux_orb)
               CALL parallel_gemm('N', 'T', nao_aux, nao_aux, nao, &
                                  1.0_dp, admm_env%work_aux_orb, admm_env%A, 0.0_dp, &
                                  admm_env%work_aux_aux)
               CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, p_env%p1_admm(ispin)%matrix, &
                                     keep_sparsity=.TRUE.)
            END DO
         END IF

      CASE DEFAULT
         CPABORT("Unknown solver for response equation requested")
      END SELECT

      IF (dft_control%do_admm) THEN
         CALL dbcsr_allocate_matrix_set(ec_env%z_admm, nspins)
         DO ispin = 1, nspins
            ALLOCATE (ec_env%z_admm(ispin)%matrix)
            CALL dbcsr_create(matrix=ec_env%z_admm(ispin)%matrix, template=matrix_s_aux(1)%matrix)
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CALL dbcsr_copy(ec_env%z_admm(ispin)%matrix, p_env%p1_admm(ispin)%matrix)
         END DO
      END IF

      ! Get rid of MO environment again
      IF (dft_control%qs_control%do_ls_scf) THEN
         DO ispin = 1, nspins
            CALL deallocate_mo_set(mos(ispin))
         END DO
         IF (ASSOCIATED(qs_env%mos)) THEN
            DO ispin = 1, SIZE(qs_env%mos)
               CALL deallocate_mo_set(qs_env%mos(ispin))
            END DO
            DEALLOCATE (qs_env%mos)
         END IF
      END IF

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
!> \param cpmos ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE response_equation_new(qs_env, p_env, cpmos, iounit)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT)      :: cpmos
      INTEGER, INTENT(IN)                                :: iounit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'response_equation_new'

      INTEGER                                            :: handle, ispin, nao, nao_aux, nocc, nspins
      LOGICAL                                            :: should_stop
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: psi0, psi1
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, matrix_ks, mo_coeff, mos)

      CALL get_qs_env(qs_env, dft_control=dft_control, matrix_ks=matrix_ks, &
                      matrix_s=matrix_s, mos=mos)
      nspins = dft_control%nspins

      ! Initialize vectors:
      ! psi0 : The ground-state MO-coefficients
      ! psi1 : The "perturbed" linear response orbitals
      ALLOCATE (psi0(nspins), psi1(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, homo=nocc)
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, ncol_global=nocc, &
                                  template_fmstruct=mo_coeff%matrix_struct)
         CALL cp_fm_create(psi0(ispin), fm_struct)
         CALL cp_fm_to_fm(mo_coeff, psi0(ispin), nocc)
         CALL cp_fm_create(psi1(ispin), fm_struct)
         CALL cp_fm_set_all(psi1(ispin), 0.0_dp)
         CALL cp_fm_struct_release(fm_struct)
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
            CALL parallel_gemm('N', 'N', nao_aux, nao, nao, &
                               1.0_dp, admm_env%A, admm_env%work_orb_orb, 0.0_dp, &
                               admm_env%work_aux_orb)
            CALL parallel_gemm('N', 'T', nao_aux, nao_aux, nao, &
                               1.0_dp, admm_env%work_aux_orb, admm_env%A, 0.0_dp, &
                               admm_env%work_aux_aux)
            CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, p_env%p1_admm(ispin)%matrix, &
                                  keep_sparsity=.TRUE.)
         END DO
      END IF

      ! Calculate Wz = 0.5*(psi1*eps*psi0^T + psi0*eps*psi1^T)
      DO ispin = 1, nspins
         CALL calculate_wz_matrix(mos(ispin), psi1(ispin), matrix_ks(ispin)%matrix, &
                                  p_env%w1(ispin)%matrix)
      END DO
      DO ispin = 1, nspins
         CALL cp_fm_release(cpmos(ispin))
      END DO
      CALL cp_fm_release(psi1)
      CALL cp_fm_release(psi0)

      CALL timestop(handle)

   END SUBROUTINE response_equation_new

! **************************************************************************************************
!> \brief Initializes vectors for MO-coefficient based linear response solver
!>        and calculates response density, and energy-weighted response density matrix
!>
!> \param qs_env ...
!> \param p_env ...
!> \param cpmos RHS of equation as Ax + b = 0 (sign of b)
!> \param iounit ...
!> \param lr_section ...
! **************************************************************************************************
   SUBROUTINE response_equation(qs_env, p_env, cpmos, iounit, lr_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: cpmos
      INTEGER, INTENT(IN)                                :: iounit
      TYPE(section_vals_type), OPTIONAL, POINTER         :: lr_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'response_equation'

      INTEGER                                            :: handle, ispin, nao, nao_aux, nocc, nspins
      LOGICAL                                            :: should_stop
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: psi0, psi1
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s, matrix_s_aux
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      CALL timeset(routineN, handle)

      ! initialized linres_control
      NULLIFY (linres_control)
      ALLOCATE (linres_control)
      linres_control%do_kernel = .TRUE.
      linres_control%lr_triplet = .FALSE.
      IF (PRESENT(lr_section)) THEN
         CALL section_vals_val_get(lr_section, "RESTART", l_val=linres_control%linres_restart)
         CALL section_vals_val_get(lr_section, "MAX_ITER", i_val=linres_control%max_iter)
         CALL section_vals_val_get(lr_section, "EPS", r_val=linres_control%eps)
         CALL section_vals_val_get(lr_section, "EPS_FILTER", r_val=linres_control%eps_filter)
         CALL section_vals_val_get(lr_section, "RESTART_EVERY", i_val=linres_control%restart_every)
         CALL section_vals_val_get(lr_section, "PRECONDITIONER", i_val=linres_control%preconditioner_type)
         CALL section_vals_val_get(lr_section, "ENERGY_GAP", r_val=linres_control%energy_gap)
      ELSE
         linres_control%linres_restart = .FALSE.
         linres_control%max_iter = 100
         linres_control%eps = 1.0e-10_dp
         linres_control%eps_filter = 1.0e-15_dp
         linres_control%restart_every = 50
         linres_control%preconditioner_type = ot_precond_full_single_inverse
         linres_control%energy_gap = 0.02_dp
      END IF

      ! initialized p_env
      CALL p_env_create(p_env, qs_env, orthogonal_orbitals=.TRUE., &
                        linres_control=linres_control)
      CALL set_qs_env(qs_env, linres_control=linres_control)
      CALL p_env_psi0_changed(p_env, qs_env)
      p_env%new_preconditioner = .TRUE.

      CALL get_qs_env(qs_env, dft_control=dft_control, mos=mos)
      !
      nspins = dft_control%nspins

      ! Initialize vectors:
      ! psi0 : The ground-state MO-coefficients
      ! psi1 : The "perturbed" linear response orbitals
      ALLOCATE (psi0(nspins), psi1(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, homo=nocc)
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, ncol_global=nocc, &
                                  template_fmstruct=mo_coeff%matrix_struct)
         CALL cp_fm_create(psi0(ispin), fm_struct)
         CALL cp_fm_to_fm(mo_coeff, psi0(ispin), nocc)
         CALL cp_fm_create(psi1(ispin), fm_struct)
         CALL cp_fm_set_all(psi1(ispin), 0.0_dp)
         CALL cp_fm_struct_release(fm_struct)
      END DO

      should_stop = .FALSE.
      ! The response solver
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
         CALL get_admm_env(qs_env%admm_env, matrix_s_aux_fit=matrix_s_aux)
         CALL dbcsr_allocate_matrix_set(p_env%p1_admm, nspins)
         DO ispin = 1, nspins
            ALLOCATE (p_env%p1_admm(ispin)%matrix)
            CALL dbcsr_create(p_env%p1_admm(ispin)%matrix, &
                              template=matrix_s_aux(1)%matrix)
            CALL dbcsr_copy(p_env%p1_admm(ispin)%matrix, matrix_s_aux(1)%matrix)
            CALL dbcsr_set(p_env%p1_admm(ispin)%matrix, 0.0_dp)
         END DO
      END IF

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
            CALL parallel_gemm('N', 'N', nao_aux, nao, nao, &
                               1.0_dp, admm_env%A, admm_env%work_orb_orb, 0.0_dp, &
                               admm_env%work_aux_orb)
            CALL parallel_gemm('N', 'T', nao_aux, nao_aux, nao, &
                               1.0_dp, admm_env%work_aux_orb, admm_env%A, 0.0_dp, &
                               admm_env%work_aux_aux)
            CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, p_env%p1_admm(ispin)%matrix, &
                                  keep_sparsity=.TRUE.)
         END DO
      END IF

      ! Calculate Wz = 0.5*(psi1*eps*psi0^T + psi0*eps*psi1^T)
      CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
      DO ispin = 1, nspins
         CALL calculate_wz_matrix(mos(ispin), psi1(ispin), matrix_ks(ispin)%matrix, &
                                  p_env%w1(ispin)%matrix)
      END DO
      CALL cp_fm_release(psi0)
      CALL cp_fm_release(psi1)

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
!> \param matrix_pz_admm Linear response density matrix in ADMM basis
!> \param matrix_wz Energy-weighted linear response density
!> \param zehartree Hartree volume response contribution to stress tensor
!> \param zexc XC volume response contribution to stress tensor
!> \param zexc_aux_fit ADMM XC volume response contribution to stress tensor
!> \param rhopz_r Response density on real space grid
!> \param p_env ...
!> \param ex_env ...
!> \param debug ...
! **************************************************************************************************
   SUBROUTINE response_force(qs_env, vh_rspace, vxc_rspace, vtau_rspace, vadmm_rspace, &
                             matrix_hz, matrix_pz, matrix_pz_admm, matrix_wz, &
                             zehartree, zexc, zexc_aux_fit, rhopz_r, p_env, ex_env, debug)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_type), INTENT(IN)                          :: vh_rspace
      TYPE(pw_type), DIMENSION(:), POINTER               :: vxc_rspace, vtau_rspace, vadmm_rspace
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz, matrix_pz, matrix_pz_admm, &
                                                            matrix_wz
      REAL(KIND=dp), OPTIONAL                            :: zehartree, zexc, zexc_aux_fit
      TYPE(pw_type), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: rhopz_r
      TYPE(qs_p_env_type), OPTIONAL                      :: p_env
      TYPE(excited_energy_type), OPTIONAL, POINTER       :: ex_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: debug

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'response_force'

      CHARACTER(LEN=default_string_length)               :: basis_type
      INTEGER                                            :: handle, iounit, ispin, mspin, myfun, &
                                                            n_rep_hf, nao, nao_aux, natom, nder, &
                                                            nimages, nocc, nspins
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL :: debug_forces, debug_stress, distribute_fock_matrix, do_ex, do_hfx, gapw, gapw_xc, &
         hfx_treat_lsd_in_core, resp_only, s_mstruct_changed, use_virial
      REAL(KIND=dp)                                      :: eh1, ehartree, ekin_mol, eps_filter, &
                                                            eps_ppnl, exc, exc_aux_fit, fconv, &
                                                            focc, hartree_gs, hartree_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ftot1, ftot2, ftot3
      REAL(KIND=dp), DIMENSION(2)                        :: total_rho_gs, total_rho_t
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_stress, pv_loc, stdeb, sttot, sttot2
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ht, matrix_pd, matrix_pza, &
                                                            matrix_s, mpa, scrm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p, mhd, mhx, mhy, mhz, &
                                                            mpd, mpz
      TYPE(dbcsr_type), POINTER                          :: dbwork
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hartree_local_type), POINTER                  :: hartree_local_gs, hartree_local_t
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(local_rho_type), POINTER                      :: local_rho_set_f, local_rho_set_gs, &
                                                            local_rho_set_t, local_rho_set_vxc, &
                                                            local_rhoz_set_admm
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux_fit, sab_orb, sac_ae, sac_ppl, &
                                                            sap_ppnl
      TYPE(oce_matrix_type), POINTER                     :: oce
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type) :: rho_tot_gspace, rho_tot_gspace_gs, rho_tot_gspace_t, rhoz_tot_gspace, &
         v_hartree_gspace_gs, v_hartree_gspace_t, v_hartree_rspace_gs, v_hartree_rspace_t, &
         vhxc_rspace, zv_hartree_gspace, zv_hartree_rspace
      TYPE(pw_type), DIMENSION(:), POINTER :: rho_g_aux, rho_g_gs, rho_g_t, rho_r_aux, rho_r_gs, &
         rho_r_t, rhoz_g, rhoz_g_aux, rhoz_g_xc, rhoz_r, rhoz_r_aux, rhoz_r_xc, tau_r_aux, tauz_r, &
         tauz_r_xc, v_xc, v_xc_tau
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit, rho_xc
      TYPE(section_vals_type), POINTER                   :: hfx_section, xc_fun_section, xc_section
      TYPE(task_list_type), POINTER                      :: task_list, task_list_aux_fit
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      IF (PRESENT(debug)) THEN
         debug_forces = debug
         debug_stress = debug
      ELSE
         debug_forces = .FALSE.
         debug_stress = .FALSE.
      END IF

      logger => cp_get_default_logger()
      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      END IF

      do_ex = .FALSE.
      IF (PRESENT(ex_env)) do_ex = .TRUE.
      IF (do_ex) THEN
         CPASSERT(PRESENT(p_env))
      END IF

      NULLIFY (ks_env, sab_orb, sac_ae, sac_ppl, sap_ppnl, virial)
      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      force=force, &
                      ks_env=ks_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      sab_orb=sab_orb, &
                      sac_ae=sac_ae, &
                      sac_ppl=sac_ppl, &
                      sap_ppnl=sap_ppnl, &
                      virial=virial)
      nspins = dft_control%nspins
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc

      IF (debug_forces) THEN
         CALL get_qs_env(qs_env, natom=natom, atomic_kind_set=atomic_kind_set)
         ALLOCATE (ftot1(3, natom))
         CALL total_qs_force(ftot1, force, atomic_kind_set)
      END IF

      ! check for virial
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      fconv = 1.0E-9_dp*pascal/cell%deth
      IF (debug_stress .AND. use_virial) THEN
         sttot = virial%pv_virial
      END IF

      !     *** If LSD, then combine alpha density and beta density to
      !     *** total density: alpha <- alpha + beta   and
      NULLIFY (mpa)
      NULLIFY (matrix_ht)
      IF (do_ex) THEN
         CALL dbcsr_allocate_matrix_set(mpa, nspins)
         DO ispin = 1, nspins
            ALLOCATE (mpa(ispin)%matrix)
            CALL dbcsr_create(mpa(ispin)%matrix, template=p_env%p1(ispin)%matrix)
            CALL dbcsr_copy(mpa(ispin)%matrix, p_env%p1(ispin)%matrix)
            CALL dbcsr_add(mpa(ispin)%matrix, ex_env%matrix_pe(ispin)%matrix, 1.0_dp, 1.0_dp)
            CALL dbcsr_set(matrix_hz(ispin)%matrix, 0.0_dp)
         END DO

         CALL dbcsr_allocate_matrix_set(matrix_ht, nspins)
         DO ispin = 1, nspins
            ALLOCATE (matrix_ht(ispin)%matrix)
            CALL dbcsr_create(matrix_ht(ispin)%matrix, template=matrix_hz(ispin)%matrix)
            CALL dbcsr_copy(matrix_ht(ispin)%matrix, matrix_hz(ispin)%matrix)
            CALL dbcsr_set(matrix_ht(ispin)%matrix, 0.0_dp)
         END DO
      ELSE
         mpa => matrix_pz
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
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_ekinetic
      CALL build_kinetic_matrix(ks_env, matrix_t=scrm, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=mpa(1)%matrix)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%kinetic(1:3, 1) - fodeb(1:3)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dT      ", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_ekinetic - stdeb)
         CALL para_env%sum(stdeb)
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
         matrix_p(ispin, 1)%matrix => mpa(ispin)%matrix
         matrix_h(ispin, 1)%matrix => scrm(ispin)%matrix
      END DO
      matrix_h(1, 1)%matrix => scrm(1)%matrix

      IF (ASSOCIATED(sac_ae)) THEN
         nder = 1
         IF (debug_forces) fodeb(1:3) = force(1)%all_potential(1:3, 1)
         IF (debug_stress .AND. use_virial) stdeb = virial%pv_ppl
         CALL build_core_ae(matrix_h, matrix_p, force, virial, .TRUE., use_virial, nder, &
                            qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ae, &
                            nimages, cell_to_index)
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%all_potential(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dHae    ", fodeb
         END IF
         IF (debug_stress .AND. use_virial) THEN
            stdeb = fconv*(virial%pv_ppl - stdeb)
            CALL para_env%sum(stdeb)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| Pz*dHae    ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF
      END IF

      IF (ASSOCIATED(sac_ppl)) THEN
         nder = 1
         IF (debug_forces) fodeb(1:3) = force(1)%gth_ppl(1:3, 1)
         IF (debug_stress .AND. use_virial) stdeb = virial%pv_ppl
         CALL build_core_ppl(matrix_h, matrix_p, force, &
                             virial, .TRUE., use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, &
                             sab_orb, sac_ppl, nimages, cell_to_index, "ORB")
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppl(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dHppl   ", fodeb
         END IF
         IF (debug_stress .AND. use_virial) THEN
            stdeb = fconv*(virial%pv_ppl - stdeb)
            CALL para_env%sum(stdeb)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| Pz*dHppl   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF
      END IF
      eps_ppnl = dft_control%qs_control%eps_ppnl
      IF (ASSOCIATED(sap_ppnl)) THEN
         nder = 1
         IF (debug_forces) fodeb(1:3) = force(1)%gth_ppnl(1:3, 1)
         IF (debug_stress .AND. use_virial) stdeb = virial%pv_ppnl
         CALL build_core_ppnl(matrix_h, matrix_p, force, &
                              virial, .TRUE., use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, &
                              sab_orb, sap_ppnl, eps_ppnl, nimages, cell_to_index, "ORB")
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppnl(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dHppnl  ", fodeb
         END IF
         IF (debug_stress .AND. use_virial) THEN
            stdeb = fconv*(virial%pv_ppnl - stdeb)
            CALL para_env%sum(stdeb)
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
            IF (debug_stress .AND. use_virial) stdeb = virial%pv_virial
            CALL build_tnadd_mat(kg_env=kg_env, matrix_p=matrix_p, force=force, virial=virial, &
                                 calculate_forces=.TRUE., use_virial=use_virial, &
                                 qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set, &
                                 particle_set=particle_set, sab_orb=sab_orb, dbcsr_dist=dbcsr_dist)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%kinetic(1:3, 1) - fodeb(1:3)
               CALL para_env%sum(fodeb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dTnadd  ", fodeb
            END IF
            IF (debug_stress .AND. use_virial) THEN
               stdeb = fconv*(virial%pv_virial - stdeb)
               CALL para_env%sum(stdeb)
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| Pz*dTnadd   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            END IF

            ! Stress-tensor update components
            IF (use_virial) THEN
               virial%pv_ekinetic = virial%pv_ekinetic + (virial%pv_virial - pv_loc)
            END IF

         END IF
      END IF

      DEALLOCATE (matrix_h)
      DEALLOCATE (matrix_p)
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
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dHcore", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_virial - sttot)
         CALL para_env%sum(stdeb)
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
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)
      ! at the moment restricted to equal xc_sections for GS and ES
      xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      CALL section_vals_val_get(xc_fun_section, "_SECTION_PARAMETERS_", i_val=myfun)
      !
      IF (gapw .OR. gapw_xc) THEN
         NULLIFY (oce, sab_orb)
         CALL get_qs_env(qs_env=qs_env, oce=oce, sab_orb=sab_orb)
         ! set up local_rho_set for GS density
         NULLIFY (local_rho_set_gs)
         CALL get_qs_env(qs_env=qs_env, rho=rho)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
         CALL local_rho_set_create(local_rho_set_gs)
         CALL allocate_rho_atom_internals(local_rho_set_gs%rho_atom_set, atomic_kind_set, &
                                          qs_kind_set, dft_control, para_env)
         CALL init_rho0(local_rho_set_gs, qs_env, dft_control%qs_control%gapw_control)
         CALL rho0_s_grid_create(pw_env, local_rho_set_gs%rho0_mpole)
         CALL calculate_rho_atom_coeff(qs_env, matrix_p(:, 1), local_rho_set_gs%rho_atom_set, &
                                       qs_kind_set, oce, sab_orb, para_env)
         CALL prepare_gapw_den(qs_env, local_rho_set_gs, do_rho0=gapw)
         ! set up local_rho_set for response density
         NULLIFY (local_rho_set_t)
         CALL local_rho_set_create(local_rho_set_t)
         CALL allocate_rho_atom_internals(local_rho_set_t%rho_atom_set, atomic_kind_set, &
                                          qs_kind_set, dft_control, para_env)
         CALL init_rho0(local_rho_set_t, qs_env, dft_control%qs_control%gapw_control, &
                        zcore=0.0_dp)
         CALL rho0_s_grid_create(pw_env, local_rho_set_t%rho0_mpole)
         CALL calculate_rho_atom_coeff(qs_env, mpa(:), local_rho_set_t%rho_atom_set, &
                                       qs_kind_set, oce, sab_orb, para_env)
         CALL prepare_gapw_den(qs_env, local_rho_set_t, do_rho0=gapw)

         ! compute soft GS potential
         ALLOCATE (rho_r_gs(nspins), rho_g_gs(nspins))
         DO ispin = 1, nspins
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_r_gs(ispin), &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_g_gs(ispin), &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         END DO
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace_gs, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         ! compute soft GS density
         total_rho_gs = 0.0_dp
         CALL pw_zero(rho_tot_gspace_gs)
         DO ispin = 1, nspins
            CALL calculate_rho_elec(ks_env=ks_env, matrix_p=matrix_p(ispin, 1)%matrix, &
                                    rho=rho_r_gs(ispin), &
                                    rho_gspace=rho_g_gs(ispin), &
                                    soft_valid=(gapw .OR. gapw_xc), &
                                    total_rho=total_rho_gs(ispin))
            CALL pw_axpy(rho_g_gs(ispin), rho_tot_gspace_gs)
         END DO
         IF (gapw) THEN
            CALL get_qs_env(qs_env, natom=natom)
            ! add rho0 contributions to GS density (only for Coulomb) only for gapw
            CALL get_rho0_mpole(local_rho_set_gs%rho0_mpole)
            CALL pw_axpy(local_rho_set_gs%rho0_mpole%rho0_s_gs, rho_tot_gspace_gs)
            ! compute GS potential
            CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace_gs, &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace_gs, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            NULLIFY (hartree_local_gs)
            CALL hartree_local_create(hartree_local_gs)
            CALL init_coulomb_local(hartree_local_gs, natom)
            CALL pw_poisson_solve(poisson_env, rho_tot_gspace_gs, hartree_gs, v_hartree_gspace_gs)
            CALL pw_transfer(v_hartree_gspace_gs, v_hartree_rspace_gs)
            CALL pw_scale(v_hartree_rspace_gs, v_hartree_rspace_gs%pw_grid%dvol)
         END IF
      END IF

      IF (gapw) THEN
         !    ! Hartree grid PAW term
         CPASSERT(do_ex)
         CPASSERT(.NOT. use_virial)
         IF (debug_forces) fodeb(1:3) = force(1)%g0s_Vh_elec(1:3, 1)
         CALL Vh_1c_gg_integrals(qs_env, hartree_gs, hartree_local_gs%ecoul_1c, local_rho_set_t, para_env, tddft=.TRUE., &
                                 local_rho_set_2nd=local_rho_set_gs, core_2nd=.FALSE.) ! n^core for GS potential
         ! 1st to define integral space, 2nd for potential, integral contributions stored on local_rho_set_gs
         CALL integrate_vhg0_rspace(qs_env, v_hartree_rspace_gs, para_env, calculate_forces=.TRUE., &
                                    local_rho_set=local_rho_set_t, local_rho_set_2nd=local_rho_set_gs)
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%g0s_Vh_elec(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: (T+Dz)*dVh[D^GS]PAWg0", fodeb
         END IF
      END IF
      IF (gapw .OR. gapw_xc) THEN
         IF (myfun /= xc_none) THEN
            ! add 1c hard and soft XC contributions
            !xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
            NULLIFY (local_rho_set_vxc)
            CALL local_rho_set_create(local_rho_set_vxc)
            CALL allocate_rho_atom_internals(local_rho_set_vxc%rho_atom_set, atomic_kind_set, &
                                             qs_kind_set, dft_control, para_env)
            CALL calculate_rho_atom_coeff(qs_env, matrix_p(:, 1), local_rho_set_vxc%rho_atom_set, &
                                          qs_kind_set, oce, sab_orb, para_env)
            CALL prepare_gapw_den(qs_env, local_rho_set_vxc, do_rho0=.FALSE.)
            ! compute hard and soft atomic contributions
            CALL calculate_vxc_atom(qs_env, .FALSE., exc1=hartree_gs, xc_section_external=xc_section, &
                                    rho_atom_set_external=local_rho_set_vxc%rho_atom_set)
         END IF ! myfun
      END IF ! gapw

      CALL pw_pool_create_pw(auxbas_pw_pool, vhxc_rspace, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      !
      ! Stress-tensor: integration contribution direct term
      ! int v_Hxc[n^in]*n^z
      IF (use_virial) THEN
         pv_loc = virial%pv_virial
      END IF

      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_virial
      IF (gapw .OR. gapw_xc) THEN
         ! vtot = v_xc + v_hartree
         DO ispin = 1, nspins
            CALL pw_zero(vhxc_rspace)
            IF (gapw) THEN
               CALL pw_transfer(v_hartree_rspace_gs, vhxc_rspace)
            ELSEIF (gapw_xc) THEN
               CALL pw_transfer(vh_rspace, vhxc_rspace)
            END IF
            CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                    hmat=scrm(ispin), pmat=mpa(ispin), &
                                    qs_env=qs_env, gapw=gapw, &
                                    calculate_forces=.TRUE.)
         END DO
         IF (myfun /= xc_none) THEN
            DO ispin = 1, nspins
               CALL pw_zero(vhxc_rspace)
               CALL pw_axpy(vxc_rspace(ispin), vhxc_rspace)
               CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                       hmat=scrm(ispin), pmat=mpa(ispin), &
                                       qs_env=qs_env, gapw=(gapw .OR. gapw_xc), &
                                       calculate_forces=.TRUE.)
            END DO
         END IF
      ELSE ! original GPW with Standard Hartree as Potential
         DO ispin = 1, nspins
            CALL pw_transfer(vh_rspace, vhxc_rspace)
            CALL pw_axpy(vxc_rspace(ispin), vhxc_rspace)
            CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                    hmat=scrm(ispin), pmat=mpa(ispin), &
                                    qs_env=qs_env, gapw=gapw, calculate_forces=.TRUE.)
         END DO
      END IF

      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: (T+Dz)*dVhxc[D^GS]   ", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_virial - pv_loc)
         CALL para_env%sum(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| INT Pz*dVhxc   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF

      IF (gapw .OR. gapw_xc) THEN
         CPASSERT(do_ex)
         !      ! HXC term
         IF (debug_forces) fodeb(1:3) = force(1)%Vhxc_atom(1:3, 1)
         IF (gapw) CALL update_ks_atom(qs_env, scrm, mpa, forces=.TRUE., tddft=.FALSE., &
                                       rho_atom_external=local_rho_set_gs%rho_atom_set)
         IF (myfun /= xc_none) CALL update_ks_atom(qs_env, scrm, mpa, forces=.TRUE., tddft=.FALSE., &
                                                   rho_atom_external=local_rho_set_vxc%rho_atom_set)
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%Vhxc_atom(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: (T+Dz)*dVhxc[D^GS]PAW ", fodeb
         END IF
         !   ! release local environments for GAPW
         IF (myfun /= xc_none) THEN
            IF (ASSOCIATED(local_rho_set_vxc)) CALL local_rho_set_release(local_rho_set_vxc)
         END IF
         IF (ASSOCIATED(local_rho_set_gs)) CALL local_rho_set_release(local_rho_set_gs)
         IF (ASSOCIATED(local_rho_set_t)) CALL local_rho_set_release(local_rho_set_t)
         IF (gapw) THEN
            IF (ASSOCIATED(hartree_local_gs)) CALL hartree_local_release(hartree_local_gs)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace_gs)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace_gs)
         END IF
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace_gs)
         IF (ASSOCIATED(rho_r_gs)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r_gs(ispin))
         END DO
         DEALLOCATE (rho_r_gs)
         END IF
         IF (ASSOCIATED(rho_g_gs)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g_gs(ispin))
         END DO
         DEALLOCATE (rho_g_gs)
         END IF
      END IF !gapw

      IF (ASSOCIATED(vtau_rspace)) THEN
         CPASSERT(.NOT. (gapw .OR. gapw_xc))
         IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
         IF (debug_stress .AND. use_virial) stdeb = virial%pv_virial
         DO ispin = 1, nspins
            CALL integrate_v_rspace(v_rspace=vtau_rspace(ispin), &
                                    hmat=scrm(ispin), pmat=mpa(ispin), &
                                    qs_env=qs_env, gapw=(gapw .OR. gapw_xc), &
                                    calculate_forces=.TRUE., compute_tau=.TRUE.)
         END DO
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dVxc_tau   ", fodeb
         END IF
         IF (debug_stress .AND. use_virial) THEN
            stdeb = fconv*(virial%pv_virial - pv_loc)
            CALL para_env%sum(stdeb)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| INT Pz*dVxc_tau   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF
      END IF
      CALL pw_pool_give_back_pw(auxbas_pw_pool, vhxc_rspace)

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
               CALL para_env%sum(fodeb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dVkg   ", fodeb
            END IF
            IF (debug_stress .AND. use_virial) THEN
               !IF (iounit > 0) WRITE(iounit, *) &
               !   "response_force | VOL 1st KG - v_KG[n_in]*n_z: ", ekin_mol
               stdeb = 1.0_dp*fconv*ekin_mol
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| VOL KG Pz*dVKG ', one_third_sum_diag(stdeb), det_3x3(stdeb)

               stdeb = fconv*(virial%pv_virial - pv_loc)
               CALL para_env%sum(stdeb)
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| INT KG Pz*dVKG  ', one_third_sum_diag(stdeb), det_3x3(stdeb)

               stdeb = fconv*virial%pv_xc
               CALL para_env%sum(stdeb)
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
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_r(ispin), &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_g(ispin), &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      END DO
      CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_tot_gspace, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, zv_hartree_rspace, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, zv_hartree_gspace, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      CALL pw_zero(rhoz_tot_gspace)
      DO ispin = 1, nspins
         CALL calculate_rho_elec(ks_env=ks_env, matrix_p=mpa(ispin)%matrix, &
                                 rho=rhoz_r(ispin), rho_gspace=rhoz_g(ispin), &
                                 soft_valid=gapw)
         CALL pw_axpy(rhoz_g(ispin), rhoz_tot_gspace)
      END DO
      IF (gapw_xc) THEN
         NULLIFY (tauz_r_xc)
         ALLOCATE (rhoz_r_xc(nspins), rhoz_g_xc(nspins))
         DO ispin = 1, nspins
            CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_r_xc(ispin), &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_g_xc(ispin), &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         END DO
         DO ispin = 1, nspins
            CALL calculate_rho_elec(ks_env=ks_env, matrix_p=mpa(ispin)%matrix, &
                                    rho=rhoz_r_xc(ispin), rho_gspace=rhoz_g_xc(ispin), &
                                    soft_valid=gapw_xc)
         END DO
      END IF

      IF (ASSOCIATED(vtau_rspace)) THEN
         CPASSERT(.NOT. (gapw .OR. gapw_xc))
         BLOCK
            TYPE(pw_type) :: work_g
            ALLOCATE (tauz_r(nspins))
            CALL pw_pool_create_pw(auxbas_pw_pool, work_g, &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            DO ispin = 1, nspins
               CALL pw_pool_create_pw(auxbas_pw_pool, tauz_r(ispin), &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL calculate_rho_elec(ks_env=ks_env, matrix_p=mpa(ispin)%matrix, &
                                       rho=tauz_r(ispin), rho_gspace=work_g, &
                                       compute_tau=.TRUE.)
            END DO
            CALL pw_pool_give_back_pw(auxbas_pw_pool, work_g)
         END BLOCK
      END IF

      !
      IF (PRESENT(rhopz_r)) THEN
         DO ispin = 1, nspins
            CALL pw_copy(rhoz_r(ispin), rhopz_r(ispin))
         END DO
      END IF

      ! Stress-tensor contribution second derivative
      ! Volume : int v_H[n^z]*n_in
      ! Volume : int epsilon_xc*n_z
      IF (use_virial) THEN

         CALL get_qs_env(qs_env, rho=rho)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

         ! Get the total input density in g-space [ions + electrons]
         CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

         h_stress(:, :) = 0.0_dp
         ! calculate associated hartree potential
         ! This term appears twice in the derivation of the equations
         ! v_H[n_in]*n_z and v_H[n_z]*n_in
         ! due to symmetry we only need to call this routine once,
         ! and count the Volume and Green function contribution
         ! which is stored in h_stress twice
         CALL pw_poisson_solve(poisson_env, &
                               density=rhoz_tot_gspace, &     ! n_z
                               ehartree=ehartree, &
                               vhartree=zv_hartree_gspace, &  ! v_H[n_z]
                               h_stress=h_stress, &
                               aux_density=rho_tot_gspace)  ! n_in

         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace)

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
            CALL para_env%sum(stdeb)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| GREEN 1st v_H[n_z]*n_in  ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            stdeb = fconv*(h_stress/REAL(para_env%num_pe, dp))
            CALL para_env%sum(stdeb)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| GREEN 2nd v_H[n_in]*n_z   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF

         ! Stress tensor volume term: \int v_xc[n_in]*n_z
         ! vxc_rspace already scaled, we need to unscale it!
         exc = 0.0_dp
         DO ispin = 1, nspins
            exc = exc + pw_integral_ab(rhoz_r(ispin), vxc_rspace(ispin))/ &
                  vxc_rspace(ispin)%pw_grid%dvol
         END DO
         IF (ASSOCIATED(vtau_rspace)) THEN
            DO ispin = 1, nspins
               exc = exc + pw_integral_ab(tauz_r(ispin), vtau_rspace(ispin))/ &
                     vtau_rspace(ispin)%pw_grid%dvol
            END DO
         END IF

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

      ELSE ! use_virial

         ! calculate associated hartree potential
         ! contribution for both T and D^Z
         CALL pw_poisson_solve(poisson_env, rhoz_tot_gspace, ehartree, zv_hartree_gspace)

      END IF ! use virial

      IF (debug_forces) fodeb(1:3) = force(1)%rho_core(1:3, 1)
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_ehartree
      CALL pw_transfer(zv_hartree_gspace, zv_hartree_rspace)
      CALL pw_scale(zv_hartree_rspace, zv_hartree_rspace%pw_grid%dvol)
      ! Getting nuclear force contribution from the core charge density (not for GAPW)
      CALL integrate_v_core_rspace(zv_hartree_rspace, qs_env)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_core(1:3, 1) - fodeb(1:3)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Vh(rhoz)*dncore ", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_ehartree - stdeb)
         CALL para_env%sum(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| INT Vh(rhoz)*dncore   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF

      !
      IF (gapw_xc) THEN
         CALL get_qs_env(qs_env=qs_env, rho_xc=rho_xc)
      ELSE
         CALL get_qs_env(qs_env=qs_env, rho=rho)
      END IF
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      END IF

      IF (use_virial) THEN
         virial%pv_xc = 0.0_dp
      END IF
      !
      NULLIFY (v_xc, v_xc_tau)
      IF (gapw_xc) THEN
         CALL create_kernel(qs_env, vxc=v_xc, vxc_tau=v_xc_tau, &
                            rho=rho_xc, rho1_r=rhoz_r_xc, rho1_g=rhoz_g_xc, tau1_r=tauz_r_xc, &
                            xc_section=xc_section, compute_virial=use_virial, virial_xc=virial%pv_xc)
      ELSE
         CALL create_kernel(qs_env, vxc=v_xc, vxc_tau=v_xc_tau, &
                            rho=rho, rho1_r=rhoz_r, rho1_g=rhoz_g, tau1_r=tauz_r, &
                            xc_section=xc_section, compute_virial=use_virial, virial_xc=virial%pv_xc)
      END IF

      IF (gapw .OR. gapw_xc) THEN
         !get local_rho_set for GS density and response potential / density
         NULLIFY (local_rho_set_t)
         CALL local_rho_set_create(local_rho_set_t)
         CALL allocate_rho_atom_internals(local_rho_set_t%rho_atom_set, atomic_kind_set, &
                                          qs_kind_set, dft_control, para_env)
         CALL init_rho0(local_rho_set_t, qs_env, dft_control%qs_control%gapw_control, &
                        zcore=0.0_dp)
         CALL rho0_s_grid_create(pw_env, local_rho_set_t%rho0_mpole)
         CALL calculate_rho_atom_coeff(qs_env, mpa(:), local_rho_set_t%rho_atom_set, &
                                       qs_kind_set, oce, sab_orb, para_env)
         CALL prepare_gapw_den(qs_env, local_rho_set_t, do_rho0=gapw)
         NULLIFY (local_rho_set_gs)
         CALL local_rho_set_create(local_rho_set_gs)
         CALL allocate_rho_atom_internals(local_rho_set_gs%rho_atom_set, atomic_kind_set, &
                                          qs_kind_set, dft_control, para_env)
         CALL init_rho0(local_rho_set_gs, qs_env, dft_control%qs_control%gapw_control)
         CALL rho0_s_grid_create(pw_env, local_rho_set_gs%rho0_mpole)
         CALL calculate_rho_atom_coeff(qs_env, matrix_p(:, 1), local_rho_set_gs%rho_atom_set, &
                                       qs_kind_set, oce, sab_orb, para_env)
         CALL prepare_gapw_den(qs_env, local_rho_set_gs, do_rho0=gapw)
         ! compute response potential
         ALLOCATE (rho_r_t(nspins), rho_g_t(nspins))
         DO ispin = 1, nspins
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_r_t(ispin), &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_g_t(ispin), &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         END DO
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace_t, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         total_rho_t = 0.0_dp
         CALL pw_zero(rho_tot_gspace_t)
         DO ispin = 1, nspins
            CALL calculate_rho_elec(ks_env=ks_env, matrix_p=mpa(ispin)%matrix, &
                                    rho=rho_r_t(ispin), &
                                    rho_gspace=rho_g_t(ispin), &
                                    soft_valid=gapw, &
                                    total_rho=total_rho_t(ispin))
            CALL pw_axpy(rho_g_t(ispin), rho_tot_gspace_t)
         END DO
         ! add rho0 contributions to response density (only for Coulomb) only for gapw
         IF (gapw) THEN
            CALL get_rho0_mpole(local_rho_set_t%rho0_mpole)
            CALL pw_axpy(local_rho_set_t%rho0_mpole%rho0_s_gs, rho_tot_gspace_t)
            ! compute response Coulomb potential
            CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace_t, &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace_t, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            NULLIFY (hartree_local_t)
            CALL hartree_local_create(hartree_local_t)
            CALL init_coulomb_local(hartree_local_t, natom)
            CALL pw_poisson_solve(poisson_env, rho_tot_gspace_t, hartree_t, v_hartree_gspace_t)
            CALL pw_transfer(v_hartree_gspace_t, v_hartree_rspace_t)
            CALL pw_scale(v_hartree_rspace_t, v_hartree_rspace_t%pw_grid%dvol)
            !
            IF (debug_forces) fodeb(1:3) = force(1)%g0s_Vh_elec(1:3, 1)
            CALL Vh_1c_gg_integrals(qs_env, hartree_t, hartree_local_t%ecoul_1c, local_rho_set_gs, para_env, tddft=.FALSE., &
                                    local_rho_set_2nd=local_rho_set_t, core_2nd=.TRUE.) ! n^core for GS potential
            CALL integrate_vhg0_rspace(qs_env, v_hartree_rspace_t, para_env, calculate_forces=.TRUE., &
                                       local_rho_set=local_rho_set_gs, local_rho_set_2nd=local_rho_set_t)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%g0s_Vh_elec(1:3, 1) - fodeb(1:3)
               CALL para_env%sum(fodeb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Vh(T)*dncore PAWg0", fodeb
            END IF
         END IF !gapw
      END IF !gapw

      IF (gapw .OR. gapw_xc) THEN
         !GAPW compute atomic fxc contributions
         IF (myfun /= xc_none) THEN
            ! local_rho_set_f
            NULLIFY (local_rho_set_f)
            CALL local_rho_set_create(local_rho_set_f)
            CALL allocate_rho_atom_internals(local_rho_set_f%rho_atom_set, atomic_kind_set, &
                                             qs_kind_set, dft_control, para_env)
            CALL calculate_rho_atom_coeff(qs_env, mpa, local_rho_set_f%rho_atom_set, &
                                          qs_kind_set, oce, sab_orb, para_env)
            CALL prepare_gapw_den(qs_env, local_rho_set_f, do_rho0=.FALSE.)
            ! add hard and soft atomic contributions
            CALL calculate_xc_2nd_deriv_atom(local_rho_set_gs%rho_atom_set, &
                                             local_rho_set_f%rho_atom_set, &
                                             qs_env, xc_section, para_env, &
                                             do_tddft=.FALSE., do_triplet=.FALSE.)
         END IF ! myfun
      END IF

      ! Stress-tensor XC-kernel GGA contribution
      IF (use_virial) THEN
         virial%pv_exc = virial%pv_exc + virial%pv_xc
         virial%pv_virial = virial%pv_virial + virial%pv_xc
      END IF

      IF (debug_stress .AND. use_virial) THEN
         stdeb = 1.0_dp*fconv*virial%pv_xc
         CALL para_env%sum(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| GGA 2nd Pin*dK*rhoz', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF

      ! Stress-tensor integral contribution of 2nd derivative terms
      IF (use_virial) THEN
         pv_loc = virial%pv_virial
      END IF

      CALL get_qs_env(qs_env=qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_virial

      DO ispin = 1, nspins
         CALL pw_scale(v_xc(ispin), v_xc(ispin)%pw_grid%dvol)
      END DO
      IF ((.NOT. (gapw)) .AND. (.NOT. gapw_xc)) THEN
         IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
         DO ispin = 1, nspins
            CALL pw_axpy(zv_hartree_rspace, v_xc(ispin)) ! Hartree potential of response density
            CALL integrate_v_rspace(qs_env=qs_env, &
                                    v_rspace=v_xc(ispin), &
                                    hmat=matrix_hz(ispin), &
                                    pmat=matrix_p(ispin, 1), &
                                    gapw=.FALSE., &
                                    calculate_forces=.TRUE.)
         END DO
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dKh*rhoz ", fodeb
         END IF
      ELSE
         IF (myfun /= xc_none) THEN
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            DO ispin = 1, nspins
               CALL integrate_v_rspace(qs_env=qs_env, &
                                       v_rspace=v_xc(ispin), &
                                       hmat=matrix_hz(ispin), &
                                       pmat=matrix_p(ispin, 1), &
                                       gapw=.TRUE., &
                                       calculate_forces=.TRUE.)
            END DO
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL para_env%sum(fodeb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dKxc*rhoz ", fodeb
            END IF
         END IF ! my_fun
         ! Coulomb T+Dz
         IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
         DO ispin = 1, nspins
            CALL pw_zero(v_xc(ispin))
            IF (gapw) THEN ! Hartree potential of response density
               CALL pw_axpy(v_hartree_rspace_t, v_xc(ispin))
            ELSEIF (gapw_xc) THEN
               CALL pw_axpy(zv_hartree_rspace, v_xc(ispin))
            END IF
            CALL integrate_v_rspace(qs_env=qs_env, &
                                    v_rspace=v_xc(ispin), &
                                    hmat=matrix_ht(ispin), &
                                    pmat=matrix_p(ispin, 1), &
                                    gapw=gapw, &
                                    calculate_forces=.TRUE.)
         END DO
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dKh*rhoz ", fodeb
         END IF
      END IF

      IF (gapw .OR. gapw_xc) THEN
         ! compute hard and soft atomic contributions
         IF (myfun /= xc_none) THEN
            IF (debug_forces) fodeb(1:3) = force(1)%Vhxc_atom(1:3, 1)
            CALL update_ks_atom(qs_env, matrix_hz, matrix_p, forces=.TRUE., tddft=.FALSE., &
                                rho_atom_external=local_rho_set_f%rho_atom_set)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%Vhxc_atom(1:3, 1) - fodeb(1:3)
               CALL para_env%sum(fodeb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P^GS*dKxc*(Dz+T) PAW", fodeb
            END IF
         END IF !myfun
         ! Coulomb contributions
         IF (gapw) THEN
            IF (debug_forces) fodeb(1:3) = force(1)%Vhxc_atom(1:3, 1)
            CALL update_ks_atom(qs_env, matrix_ht, matrix_p, forces=.TRUE., tddft=.FALSE., &
                                rho_atom_external=local_rho_set_t%rho_atom_set)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%Vhxc_atom(1:3, 1) - fodeb(1:3)
               CALL para_env%sum(fodeb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P^GS*dKh*(Dz+T) PAW", fodeb
            END IF
         END IF
         ! add Coulomb and XC
         DO ispin = 1, nspins
            CALL dbcsr_add(matrix_hz(ispin)%matrix, matrix_ht(ispin)%matrix, 1.0_dp, 1.0_dp)
         END DO

         ! release
         IF (myfun /= xc_none) THEN
            IF (ASSOCIATED(local_rho_set_f)) CALL local_rho_set_release(local_rho_set_f)
         END IF
         IF (ASSOCIATED(local_rho_set_t)) CALL local_rho_set_release(local_rho_set_t)
         IF (ASSOCIATED(local_rho_set_gs)) CALL local_rho_set_release(local_rho_set_gs)
         IF (gapw) THEN
            IF (ASSOCIATED(hartree_local_t)) CALL hartree_local_release(hartree_local_t)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace_t)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace_t)
         END IF
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace_t)
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r_t(ispin))
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g_t(ispin))
         END DO
         DEALLOCATE (rho_r_t, rho_g_t)
      END IF ! gapw

      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_virial - stdeb)
         CALL para_env%sum(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| INT 2nd f_Hxc[Pz]*Pin', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF
      !
      IF (ASSOCIATED(v_xc_tau)) THEN
         CPASSERT(.NOT. (gapw .OR. gapw_xc))
         IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
         IF (debug_stress .AND. use_virial) stdeb = virial%pv_virial
         DO ispin = 1, nspins
            CALL pw_scale(v_xc_tau(ispin), v_xc_tau(ispin)%pw_grid%dvol)
            CALL integrate_v_rspace(qs_env=qs_env, &
                                    v_rspace=v_xc_tau(ispin), &
                                    hmat=matrix_hz(ispin), &
                                    pmat=matrix_p(ispin, 1), &
                                    compute_tau=.TRUE., &
                                    gapw=(gapw .OR. gapw_xc), &
                                    calculate_forces=.TRUE.)
         END DO
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dKtau*tauz ", fodeb
         END IF
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_virial - stdeb)
         CALL para_env%sum(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| INT 2nd f_xctau[Pz]*Pin', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF
      !! Stress-tensor integral contribution of 2nd derivative terms
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
               CALL para_env%sum(fodeb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*d(Kkg)*rhoz ", fodeb
            END IF
            IF (debug_stress .AND. use_virial) THEN
               stdeb = fconv*(virial%pv_virial - pv_loc)
               CALL para_env%sum(stdeb)
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| INT KG Pin*d(KKG)*rhoz    ', one_third_sum_diag(stdeb), det_3x3(stdeb)

               stdeb = fconv*(virial%pv_xc)
               CALL para_env%sum(stdeb)
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
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_tot_gspace)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, zv_hartree_gspace)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, zv_hartree_rspace)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_r(ispin))
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_g(ispin))
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin))
      END DO
      DEALLOCATE (rhoz_r, rhoz_g, v_xc)
      IF (gapw_xc) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_r_xc(ispin))
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_g_xc(ispin))
         END DO
         DEALLOCATE (rhoz_r_xc, rhoz_g_xc)
      END IF
      IF (ASSOCIATED(v_xc_tau)) THEN
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, tauz_r(ispin))
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc_tau(ispin))
      END DO
      DEALLOCATE (tauz_r, v_xc_tau)
      END IF
      IF (debug_forces) THEN
         ALLOCATE (ftot3(3, natom))
         CALL total_qs_force(ftot3, force, atomic_kind_set)
         fodeb(1:3) = ftot3(1:3, 1) - ftot2(1:3, 1)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*V(rhoz)", fodeb
      END IF
      CALL dbcsr_deallocate_matrix_set(scrm)
      CALL dbcsr_deallocate_matrix_set(matrix_ht)

      ! -----------------------------------------
      ! Apply ADMM exchange correction
      ! -----------------------------------------

      IF (dft_control%do_admm) THEN
         ! volume term
         exc_aux_fit = 0.0_dp

         IF (qs_env%admm_env%aux_exch_func == do_admm_aux_exch_func_none) THEN
            ! nothing to do
            NULLIFY (mpz, mhz, mhx, mhy)
         ELSE
            ! add ADMM xc_section_aux terms: Pz*Vxc + P0*K0[rhoz]
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CALL get_admm_env(admm_env, rho_aux_fit=rho_aux_fit, matrix_s_aux_fit=scrm, &
                              task_list_aux_fit=task_list_aux_fit)
            !
            NULLIFY (mpz, mhz, mhx, mhy)
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
               IF (do_ex) THEN
                  CALL dbcsr_create(mpz(ispin, 1)%matrix, template=p_env%p1_admm(ispin)%matrix)
                  CALL dbcsr_copy(mpz(ispin, 1)%matrix, p_env%p1_admm(ispin)%matrix)
                  CALL dbcsr_add(mpz(ispin, 1)%matrix, ex_env%matrix_pe_admm(ispin)%matrix, &
                                 1.0_dp, 1.0_dp)
               ELSE
                  CALL dbcsr_create(mpz(ispin, 1)%matrix, template=matrix_pz_admm(ispin)%matrix)
                  CALL dbcsr_copy(mpz(ispin, 1)%matrix, matrix_pz_admm(ispin)%matrix)
               END IF
            END DO
            !
            xc_section => admm_env%xc_section_aux
            ! Stress-tensor: integration contribution direct term
            ! int Pz*v_xc[rho_admm]
            IF (use_virial) THEN
               pv_loc = virial%pv_virial
            END IF

            basis_type = "AUX_FIT"
            task_list => task_list_aux_fit
            IF (admm_env%do_gapw) THEN
               basis_type = "AUX_FIT_SOFT"
               task_list => admm_env%admm_gapw_env%task_list
            END IF
            !
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            IF (debug_stress .AND. use_virial) stdeb = virial%pv_virial
            DO ispin = 1, nspins
               CALL integrate_v_rspace(v_rspace=vadmm_rspace(ispin), &
                                       hmat=mhx(ispin, 1), pmat=mpz(ispin, 1), &
                                       qs_env=qs_env, calculate_forces=.TRUE., &
                                       basis_type=basis_type, task_list_external=task_list)
            END DO
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL para_env%sum(fodeb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*Vxc(rho_admm)", fodeb
            END IF
            IF (debug_stress .AND. use_virial) THEN
               stdeb = fconv*(virial%pv_virial - pv_loc)
               CALL para_env%sum(stdeb)
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| INT 1st Pz*dVxc(rho_admm)   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            END IF
            ! Stress-tensor Pz_admm*v_xc[rho_admm]
            IF (use_virial) THEN
               virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
            END IF
            !
            IF (admm_env%do_gapw) THEN
               CALL get_admm_env(admm_env, sab_aux_fit=sab_aux_fit)
               IF (debug_forces) fodeb(1:3) = force(1)%Vhxc_atom(1:3, 1)
               CALL update_ks_atom(qs_env, mhx(:, 1), mpz(:, 1), forces=.TRUE., tddft=.FALSE., &
                                   rho_atom_external=ex_env%local_rho_set_admm%rho_atom_set, &
                                   kind_set_external=admm_env%admm_gapw_env%admm_kind_set, &
                                   oce_external=admm_env%admm_gapw_env%oce, &
                                   sab_external=sab_aux_fit)
               IF (debug_forces) THEN
                  fodeb(1:3) = force(1)%Vhxc_atom(1:3, 1) - fodeb(1:3)
                  CALL para_env%sum(fodeb)
                  IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*Vxc(rho_admm)PAW", fodeb
               END IF
            END IF
            !
            NULLIFY (rho_g_aux, rho_r_aux, tau_r_aux)
            CALL qs_rho_get(rho_aux_fit, rho_r=rho_r_aux, rho_g=rho_g_aux, tau_r=tau_r_aux)
            ! rhoz_aux
            NULLIFY (rhoz_g_aux, rhoz_r_aux)
            ALLOCATE (rhoz_r_aux(nspins), rhoz_g_aux(nspins))
            DO ispin = 1, nspins
               CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_r_aux(ispin), &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_g_aux(ispin), &
                                      use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            END DO
            DO ispin = 1, nspins
               CALL calculate_rho_elec(ks_env=ks_env, matrix_p=mpz(ispin, 1)%matrix, &
                                       rho=rhoz_r_aux(ispin), rho_gspace=rhoz_g_aux(ispin), &
                                       basis_type="AUX_FIT", &
                                       task_list_external=task_list_aux_fit)
            END DO
            !
            ! Add ADMM volume contribution to stress tensor
            IF (use_virial) THEN

               ! Stress tensor volume term: \int v_xc[n_in_admm]*n_z_admm
               ! vadmm_rspace already scaled, we need to unscale it!
               DO ispin = 1, nspins
                  exc_aux_fit = exc_aux_fit + pw_integral_ab(rhoz_r_aux(ispin), vadmm_rspace(ispin))/ &
                                vadmm_rspace(ispin)%pw_grid%dvol
               END DO

               IF (debug_stress) THEN
                  stdeb = -1.0_dp*fconv*exc_aux_fit
                  IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T43,2(1X,ES19.11))") &
                     'STRESS| VOL 1st eps_XC[n_in_admm]*n_z_admm', one_third_sum_diag(stdeb), det_3x3(stdeb)
               END IF

            END IF
            !
            NULLIFY (v_xc)

            IF (use_virial) virial%pv_xc = 0.0_dp

            CALL create_kernel(qs_env=qs_env, &
                               vxc=v_xc, &
                               vxc_tau=v_xc_tau, &
                               rho=rho_aux_fit, &
                               rho1_r=rhoz_r_aux, &
                               rho1_g=rhoz_g_aux, &
                               tau1_r=tau_r_aux, &
                               xc_section=xc_section, &
                               compute_virial=use_virial, &
                               virial_xc=virial%pv_xc)

            ! Stress-tensor ADMM-kernel GGA contribution
            IF (use_virial) THEN
               virial%pv_exc = virial%pv_exc + virial%pv_xc
               virial%pv_virial = virial%pv_virial + virial%pv_xc
            END IF

            IF (debug_stress .AND. use_virial) THEN
               stdeb = 1.0_dp*fconv*virial%pv_xc
               CALL para_env%sum(stdeb)
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| GGA 2nd Pin_admm*dK*rhoz_admm', one_third_sum_diag(stdeb), det_3x3(stdeb)
            END IF
            !
            CALL qs_rho_get(rho_aux_fit, rho_ao_kp=matrix_p)
            ! Stress-tensor Pin*dK*rhoz_admm
            IF (use_virial) THEN
               virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
            END IF
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            IF (debug_stress .AND. use_virial) stdeb = virial%pv_virial
            DO ispin = 1, nspins
               CALL dbcsr_set(mhy(ispin, 1)%matrix, 0.0_dp)
               CALL pw_scale(v_xc(ispin), v_xc(ispin)%pw_grid%dvol)
               CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                       hmat=mhy(ispin, 1), pmat=matrix_p(ispin, 1), &
                                       calculate_forces=.TRUE., basis_type="AUX_FIT", &
                                       task_list_external=task_list_aux_fit)
            END DO
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL para_env%sum(fodeb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dK*rhoz_admm ", fodeb
            END IF
            IF (debug_stress .AND. use_virial) THEN
               stdeb = fconv*(virial%pv_virial - pv_loc)
               CALL para_env%sum(stdeb)
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| INT 2nd Pin*dK*rhoz_admm   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            END IF
            ! Stress-tensor Pin*dK*rhoz_admm
            IF (use_virial) THEN
               virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
            END IF
            DO ispin = 1, nspins
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin))
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_r_aux(ispin))
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_g_aux(ispin))
            END DO
            DEALLOCATE (v_xc, rhoz_r_aux, rhoz_g_aux)
            !
            IF (admm_env%do_gapw) THEN
               CALL local_rho_set_create(local_rhoz_set_admm)
               CALL allocate_rho_atom_internals(local_rhoz_set_admm%rho_atom_set, atomic_kind_set, &
                                                admm_env%admm_gapw_env%admm_kind_set, dft_control, para_env)
               CALL calculate_rho_atom_coeff(qs_env, mpz(:, 1), local_rhoz_set_admm%rho_atom_set, &
                                             admm_env%admm_gapw_env%admm_kind_set, &
                                             admm_env%admm_gapw_env%oce, sab_aux_fit, para_env)
               CALL prepare_gapw_den(qs_env, local_rho_set=local_rhoz_set_admm, &
                                     do_rho0=.FALSE., kind_set_external=admm_env%admm_gapw_env%admm_kind_set)
               !compute the potential due to atomic densities
               CALL calculate_xc_2nd_deriv_atom(ex_env%local_rho_set_admm%rho_atom_set, local_rhoz_set_admm%rho_atom_set, &
                                                qs_env, xc_section, para_env, do_tddft=.FALSE., do_triplet=.FALSE., &
                                                kind_set_external=admm_env%admm_gapw_env%admm_kind_set)
               !
               IF (debug_forces) fodeb(1:3) = force(1)%Vhxc_atom(1:3, 1)
               CALL update_ks_atom(qs_env, mhy(:, 1), matrix_p(:, 1), forces=.TRUE., tddft=.FALSE., &
                                   rho_atom_external=local_rhoz_set_admm%rho_atom_set, &
                                   kind_set_external=admm_env%admm_gapw_env%admm_kind_set, &
                                   oce_external=admm_env%admm_gapw_env%oce, &
                                   sab_external=sab_aux_fit)
               IF (debug_forces) THEN
                  fodeb(1:3) = force(1)%Vhxc_atom(1:3, 1) - fodeb(1:3)
                  CALL para_env%sum(fodeb)
                  IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dK*rhoz_admm[PAW] ", fodeb
               END IF
               CALL local_rho_set_release(local_rhoz_set_admm)
            END IF
            !
            nao = admm_env%nao_orb
            nao_aux = admm_env%nao_aux_fit
            ALLOCATE (dbwork)
            CALL dbcsr_create(dbwork, template=matrix_hz(1)%matrix)
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(mhy(ispin, 1)%matrix, admm_env%A, &
                                            admm_env%work_aux_orb, nao)
               CALL parallel_gemm('T', 'N', nao, nao, nao_aux, &
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
         END IF ! qs_env%admm_env%aux_exch_func == do_admm_aux_exch_func_none
      END IF ! do_admm

      ! -----------------------------------------
      !  HFX
      ! -----------------------------------------

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
         IF (use_virial) virial%pv_fock_4c = 0.0_dp
         !
         CALL get_qs_env(qs_env=qs_env, rho=rho, x_data=x_data, &
                         s_mstruct_changed=s_mstruct_changed)
         distribute_fock_matrix = .TRUE.

         ! -----------------------------------------
         !  HFX-ADMM
         ! -----------------------------------------
         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env=qs_env, admm_env=admm_env)
            CALL get_admm_env(admm_env, matrix_s_aux_fit=scrm, rho_aux_fit=rho_aux_fit)
            CALL qs_rho_get(rho_aux_fit, rho_ao_kp=matrix_p)
            NULLIFY (mpz, mhz, mpd, mhd)
            CALL dbcsr_allocate_matrix_set(mpz, nspins, 1)
            CALL dbcsr_allocate_matrix_set(mhz, nspins, 1)
            CALL dbcsr_allocate_matrix_set(mpd, nspins, 1)
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
               IF (do_ex) THEN
                  CALL dbcsr_create(mpz(ispin, 1)%matrix, template=scrm(1)%matrix)
                  CALL dbcsr_copy(mpz(ispin, 1)%matrix, p_env%p1_admm(ispin)%matrix)
                  CALL dbcsr_add(mpz(ispin, 1)%matrix, ex_env%matrix_pe_admm(ispin)%matrix, &
                                 1.0_dp, 1.0_dp)
               ELSE
                  CALL dbcsr_create(mpz(ispin, 1)%matrix, template=scrm(1)%matrix)
                  CALL dbcsr_copy(mpz(ispin, 1)%matrix, matrix_pz_admm(ispin)%matrix)
               END IF
               mpd(ispin, 1)%matrix => matrix_p(ispin, 1)%matrix
            END DO
            !
            IF (x_data(1, 1)%do_hfx_ri) THEN

               eh1 = 0.0_dp
               CALL hfx_ri_update_ks(qs_env, x_data(1, 1)%ri_data, mhz, eh1, rho_ao=mpz, &
                                     geometry_did_change=s_mstruct_changed, nspins=nspins, &
                                     hf_fraction=x_data(1, 1)%general_parameter%fraction)

               eh1 = 0.0_dp
               CALL hfx_ri_update_ks(qs_env, x_data(1, 1)%ri_data, mhd, eh1, rho_ao=mpd, &
                                     geometry_did_change=s_mstruct_changed, nspins=nspins, &
                                     hf_fraction=x_data(1, 1)%general_parameter%fraction)

            ELSE
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
            END IF
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
               CALL parallel_gemm('T', 'N', nao, nao, nao_aux, &
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
            CALL admm_projection_derivative(qs_env, mhd(:, 1), mpa)
            CALL admm_projection_derivative(qs_env, mhz(:, 1), matrix_pd)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%overlap_admm(1:3, 1) - fodeb(1:3)
               CALL para_env%sum(fodeb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*hfx*S' ", fodeb
            END IF
            CALL dbcsr_deallocate_matrix_set(mpz)
            CALL dbcsr_deallocate_matrix_set(mhz)
            CALL dbcsr_deallocate_matrix_set(mhd)
            IF (ASSOCIATED(mhx) .AND. ASSOCIATED(mhy)) THEN
               CALL dbcsr_deallocate_matrix_set(mhx)
               CALL dbcsr_deallocate_matrix_set(mhy)
            END IF
            DEALLOCATE (mpd)
         ELSE
            ! -----------------------------------------
            !  conventional HFX
            ! -----------------------------------------
            ALLOCATE (mpz(nspins, 1), mhz(nspins, 1))
            DO ispin = 1, nspins
               mhz(ispin, 1)%matrix => matrix_hz(ispin)%matrix
               mpz(ispin, 1)%matrix => mpa(ispin)%matrix
            END DO

            IF (x_data(1, 1)%do_hfx_ri) THEN

               eh1 = 0.0_dp
               CALL hfx_ri_update_ks(qs_env, x_data(1, 1)%ri_data, mhz, eh1, rho_ao=mpz, &
                                     geometry_did_change=s_mstruct_changed, nspins=nspins, &
                                     hf_fraction=x_data(1, 1)%general_parameter%fraction)
            ELSE
               DO ispin = 1, mspin
                  eh1 = 0.0
                  CALL integrate_four_center(qs_env, x_data, mhz, eh1, mpz, hfx_section, &
                                             para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                             ispin=ispin)
               END DO
            END IF
            DEALLOCATE (mhz, mpz)
         END IF

         ! -----------------------------------------
         !  HFX FORCES
         ! -----------------------------------------

         resp_only = .TRUE.
         IF (debug_forces) fodeb(1:3) = force(1)%fock_4c(1:3, 1)
         IF (dft_control%do_admm) THEN
            ! -----------------------------------------
            !  HFX-ADMM FORCES
            ! -----------------------------------------
            CALL qs_rho_get(rho_aux_fit, rho_ao_kp=matrix_p)
            NULLIFY (matrix_pza)
            CALL dbcsr_allocate_matrix_set(matrix_pza, nspins)
            DO ispin = 1, nspins
               ALLOCATE (matrix_pza(ispin)%matrix)
               IF (do_ex) THEN
                  CALL dbcsr_create(matrix_pza(ispin)%matrix, template=p_env%p1_admm(ispin)%matrix)
                  CALL dbcsr_copy(matrix_pza(ispin)%matrix, p_env%p1_admm(ispin)%matrix)
                  CALL dbcsr_add(matrix_pza(ispin)%matrix, ex_env%matrix_pe_admm(ispin)%matrix, &
                                 1.0_dp, 1.0_dp)
               ELSE
                  CALL dbcsr_create(matrix_pza(ispin)%matrix, template=matrix_pz_admm(ispin)%matrix)
                  CALL dbcsr_copy(matrix_pza(ispin)%matrix, matrix_pz_admm(ispin)%matrix)
               END IF
            END DO
            IF (x_data(1, 1)%do_hfx_ri) THEN

               CALL hfx_ri_update_forces(qs_env, x_data(1, 1)%ri_data, nspins, &
                                         x_data(1, 1)%general_parameter%fraction, &
                                         rho_ao=matrix_p, rho_ao_resp=matrix_pza, &
                                         use_virial=use_virial, resp_only=resp_only)
            ELSE
               CALL derivatives_four_center(qs_env, matrix_p, matrix_pza, hfx_section, para_env, &
                                            1, use_virial, resp_only=resp_only)
            END IF
            CALL dbcsr_deallocate_matrix_set(matrix_pza)
         ELSE
            ! -----------------------------------------
            !  conventional HFX FORCES
            ! -----------------------------------------
            CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
            IF (x_data(1, 1)%do_hfx_ri) THEN

               CALL hfx_ri_update_forces(qs_env, x_data(1, 1)%ri_data, nspins, &
                                         x_data(1, 1)%general_parameter%fraction, &
                                         rho_ao=matrix_p, rho_ao_resp=mpa, &
                                         use_virial=use_virial, resp_only=resp_only)
            ELSE
               CALL derivatives_four_center(qs_env, matrix_p, mpa, hfx_section, para_env, &
                                            1, use_virial, resp_only=resp_only)
            END IF
         END IF ! do_admm

         IF (use_virial) THEN
            virial%pv_exx = virial%pv_exx - virial%pv_fock_4c
            virial%pv_virial = virial%pv_virial - virial%pv_fock_4c
            virial%pv_calculate = .FALSE.
         END IF

         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%fock_4c(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*hfx ", fodeb
         END IF
         IF (debug_stress .AND. use_virial) THEN
            stdeb = -1.0_dp*fconv*virial%pv_fock_4c
            CALL para_env%sum(stdeb)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| Pz*hfx  ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF
      END IF ! do_hfx

      ! Stress-tensor volume contributions
      ! These need to be applied at the end of qs_force
      IF (use_virial) THEN
         ! Adding mixed Hartree energy twice, due to symmetry
         zehartree = zehartree + 2.0_dp*ehartree
         zexc = zexc + exc
         ! ADMM contribution handled differently in qs_force
         IF (dft_control%do_admm) THEN
            zexc_aux_fit = zexc_aux_fit + exc_aux_fit
         END IF
      END IF

      ! Overlap matrix
      ! H(drho+dz) + Wz
      ! If ground-state density matrix solved by diagonalization, then use this
      IF (dft_control%qs_control%do_ls_scf) THEN
         ! Ground-state density has been calculated by LS
         eps_filter = dft_control%qs_control%eps_filter_matrix
         CALL calculate_whz_ao_matrix(qs_env, matrix_hz, matrix_wz, eps_filter)
      ELSE
         IF (do_ex) THEN
            matrix_wz => p_env%w1
         END IF
         focc = 1.0_dp
         IF (nspins == 1) focc = 2.0_dp
         CALL get_qs_env(qs_env, mos=mos)
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin), homo=nocc)
            CALL calculate_whz_matrix(mos(ispin)%mo_coeff, matrix_hz(ispin)%matrix, &
                                      matrix_wz(ispin)%matrix, focc, nocc)
         END DO
      END IF
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_wz(1)%matrix, matrix_wz(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF

      IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
      IF (debug_stress .AND. use_virial) stdeb = virial%pv_overlap
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
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Wz*dS ", fodeb
      END IF
      IF (debug_stress .AND. use_virial) THEN
         stdeb = fconv*(virial%pv_overlap - stdeb)
         CALL para_env%sum(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| WHz   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF
      CALL dbcsr_deallocate_matrix_set(scrm)

      IF (debug_forces) THEN
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1) - ftot1(1:3, 1)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Response Force", fodeb
         fodeb(1:3) = ftot2(1:3, 1)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Total Force ", fodeb
         DEALLOCATE (ftot1, ftot2, ftot3)
      END IF

      IF (do_ex) THEN
         CALL dbcsr_deallocate_matrix_set(mpa)
         CALL dbcsr_deallocate_matrix_set(matrix_hz)
      END IF

      CALL timestop(handle)

   END SUBROUTINE response_force

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param matrix_hz ...
!> \param ex_env ...
!> \param debug ...
! **************************************************************************************************
   SUBROUTINE response_force_xtb(qs_env, p_env, matrix_hz, ex_env, debug)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz
      TYPE(excited_energy_type), OPTIONAL, POINTER       :: ex_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: debug

      CHARACTER(LEN=*), PARAMETER :: routineN = 'response_force_xtb'

      INTEGER                                            :: atom_a, handle, iatom, ikind, iounit, &
                                                            is, ispin, na, natom, natorb, nimages, &
                                                            nkind, nocc, ns, nsgf, nspins
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(5)                              :: occ
      LOGICAL                                            :: debug_forces, do_ex, use_virial
      REAL(KIND=dp)                                      :: focc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mcharge, mcharge1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aocg, aocg1, charges, charges1, ftot1, &
                                                            ftot2
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pz, matrix_wz, mpa, p_matrix, scrm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dbcsr_type), POINTER                          :: s_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      IF (PRESENT(debug)) THEN
         debug_forces = debug
      ELSE
         debug_forces = .FALSE.
      END IF

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      END IF

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
         CALL para_env%sum(fodeb)
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
         CALL get_mo_set(mo_set=mos(ispin), homo=nocc)
         CALL calculate_whz_matrix(mos(ispin)%mo_coeff, matrix_hz(ispin)%matrix, &
                                   matrix_wz(ispin)%matrix, focc, nocc)
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
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Wz*dS ", fodeb
      END IF
      CALL dbcsr_deallocate_matrix_set(scrm)

      IF (debug_forces) THEN
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1) - ftot1(1:3, 1)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Response Force", fodeb
         DEALLOCATE (ftot1, ftot2)
      END IF

      IF (do_ex) THEN
         CALL dbcsr_deallocate_matrix_set(mpa)
      END IF

      CALL timestop(handle)

   END SUBROUTINE response_force_xtb

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

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_whz_ao_matrix'

      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: scaling
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_type)                                   :: matrix_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
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
