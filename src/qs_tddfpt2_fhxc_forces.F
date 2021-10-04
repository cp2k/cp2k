!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_fhxc_forces
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type,&
                                              stda_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_row_scale,&
                                              cp_fm_schur_product
   USE cp_fm_pool_types,                ONLY: fm_pool_create_fm,&
                                              fm_pool_give_back_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_vectorssum
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_filter, &
        dbcsr_get_block_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_transposed, dbcsr_type, dbcsr_type_no_symmetry
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE ewald_methods_tb,                ONLY: tb_ewald_overlap,&
                                              tb_spme_evaluate
   USE ewald_pw_types,                  ONLY: ewald_pw_type
   USE exstates_types,                  ONLY: excited_energy_type
   USE hfx_derivatives,                 ONLY: derivatives_four_center
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none,&
                                              tddfpt_kernel_full,&
                                              xc_kernel_method_analytic,&
                                              xc_kernel_method_best,&
                                              xc_kernel_method_numeric,&
                                              xc_none
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: oorootpi
   USE message_passing,                 ONLY: mp_sum
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
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
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: add_qs_force,&
                                              qs_force_type
   USE qs_fxc,                          ONLY: qs_fgxc_create,&
                                              qs_fgxc_release
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kernel_types,                 ONLY: full_kernel_env_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_force,&
                                              build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_tddfpt2_stda_types,           ONLY: stda_env_type
   USE qs_tddfpt2_stda_utils,           ONLY: get_lowdin_x,&
                                              setup_gamma
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos,&
                                              tddfpt_work_matrices
   USE task_list_types,                 ONLY: task_list_type
   USE util,                            ONLY: get_limit
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_fhxc_forces'

   PUBLIC :: fhxc_force, stda_force

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Calculate direct tddft forces
!> \param qs_env ...
!> \param ex_env ...
!> \param gs_mos ...
!> \param full_kernel ...
!> \param debug_forces ...
!> \par History
!>    * 01.2020 screated [JGH]
! **************************************************************************************************
   SUBROUTINE fhxc_force(qs_env, ex_env, gs_mos, full_kernel, debug_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(full_kernel_env_type), INTENT(IN)             :: full_kernel
      LOGICAL, INTENT(IN)                                :: debug_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fhxc_force', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iounit, ispin, mspin, myfun, &
                                                            n_rep_hf, nao, nao_aux, norb, nspins
      LOGICAL :: distribute_fock_matrix, do_admm, do_analytic, do_hfx, do_numeric, &
         hfx_treat_lsd_in_core, is_rks_triplets, s_mstruct_changed, use_virial
      REAL(KIND=dp)                                      :: eh1, focc, fval, total_rhox, &
                                                            total_rhox_aux, xehartree
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: cpmos, evect
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, fm_struct_mat
      TYPE(cp_fm_type), POINTER                          :: cvcmat, mos, vcvec
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: matrix_fx, matrix_gx, matrix_hfx, &
         matrix_hfx_admm, matrix_hx, matrix_p, matrix_p_admm, matrix_px1, matrix_px1_admm, &
         matrix_s, matrix_s_aux_fit, matrix_wx1, mdum, mfx, mgx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mhe, mpe
      TYPE(dbcsr_type), POINTER                          :: dbwork
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhox_tot_gspace, xv_hartree_gspace, &
                                                            xv_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: fxc_rho, fxc_tau, gxc_rho, gxc_tau, &
                                                            rho_g_aux, rho_r_aux, rhox_g, &
                                                            rhox_g_aux, rhox_r, rhox_r_aux, &
                                                            tau_r_aux
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit, rhox, rhox_aux
      TYPE(section_vals_type), POINTER                   :: hfx_section, xc_section
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      ENDIF

      CALL get_qs_env(qs_env, dft_control=dft_control)
      tddfpt_control => dft_control%tddfpt2_control
      nspins = dft_control%nspins
      is_rks_triplets = tddfpt_control%rks_triplets .AND. (nspins == 1)
      CPASSERT(tddfpt_control%kernel == tddfpt_kernel_full)
      do_hfx = tddfpt_control%do_hfx
      do_admm = tddfpt_control%do_admm

      evect => ex_env%evect
      matrix_px1 => ex_env%matrix_px1
      matrix_px1_admm => ex_env%matrix_px1_admm

      focc = 1.0_dp
      IF (nspins == 2) focc = 0.5_dp
      DO ispin = 1, nspins
         CALL dbcsr_set(matrix_px1(ispin)%matrix, 0.0_dp)
         CALL cp_fm_get_info(evect(ispin)%matrix, ncol_global=norb)
         CALL cp_dbcsr_plus_fm_fm_t(matrix_px1(ispin)%matrix, &
                                    matrix_v=gs_mos(ispin)%mos_occ, &
                                    matrix_g=evect(ispin)%matrix, &
                                    ncol=norb, alpha=focc)
         CALL cp_dbcsr_plus_fm_fm_t(matrix_px1(ispin)%matrix, &
                                    matrix_v=evect(ispin)%matrix, &
                                    matrix_g=gs_mos(ispin)%mos_occ, &
                                    ncol=norb, alpha=focc)
      ENDDO
      !
      IF (do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         nao_aux = admm_env%nao_aux_fit
         nao = admm_env%nao_orb
         !
         DO ispin = 1, nspins
            CALL copy_dbcsr_to_fm(matrix_px1(ispin)%matrix, admm_env%work_orb_orb)
            CALL cp_gemm('N', 'N', nao_aux, nao, nao, &
                         1.0_dp, admm_env%A, admm_env%work_orb_orb, 0.0_dp, &
                         admm_env%work_aux_orb)
            CALL cp_gemm('N', 'T', nao_aux, nao_aux, nao, &
                         1.0_dp, admm_env%work_aux_orb, admm_env%A, 0.0_dp, &
                         admm_env%work_aux_aux)
            CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, matrix_px1_admm(ispin)%matrix, &
                                  keep_sparsity=.TRUE.)
         END DO
      END IF
      !
      CALL get_qs_env(qs_env, ks_env=ks_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      ALLOCATE (rhox_r(nspins), rhox_g(nspins))
      DO ispin = 1, nspins
         NULLIFY (rhox_r(ispin)%pw, rhox_g(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhox_r(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhox_g(ispin)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      ENDDO
      CALL pw_pool_create_pw(auxbas_pw_pool, rhox_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      CALL pw_zero(rhox_tot_gspace%pw)
      DO ispin = 1, nspins
         IF (nspins == 2) CALL dbcsr_scale(matrix_px1(ispin)%matrix, 2.0_dp)
         CALL calculate_rho_elec(ks_env=ks_env, matrix_p=matrix_px1(ispin)%matrix, &
                                 rho=rhox_r(ispin), rho_gspace=rhox_g(ispin), &
                                 total_rho=total_rhox)
         CALL pw_axpy(rhox_g(ispin)%pw, rhox_tot_gspace%pw)
         IF (nspins == 2) CALL dbcsr_scale(matrix_px1(ispin)%matrix, 0.5_dp)
      END DO

      CALL get_qs_env(qs_env, matrix_s=matrix_s, force=force, para_env=para_env)

      IF (.NOT. is_rks_triplets) THEN
         CALL pw_pool_create_pw(auxbas_pw_pool, xv_hartree_rspace%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, xv_hartree_gspace%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         ! calculate associated hartree potential
         CALL pw_poisson_solve(poisson_env, rhox_tot_gspace%pw, xehartree, &
                               xv_hartree_gspace%pw)
         CALL pw_transfer(xv_hartree_gspace%pw, xv_hartree_rspace%pw)
         CALL pw_scale(xv_hartree_rspace%pw, xv_hartree_rspace%pw%pw_grid%dvol)
         !
         IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
         NULLIFY (matrix_hx)
         CALL dbcsr_allocate_matrix_set(matrix_hx, nspins)
         DO ispin = 1, nspins
            ALLOCATE (matrix_hx(ispin)%matrix)
            CALL dbcsr_create(matrix_hx(ispin)%matrix, template=matrix_s(1)%matrix)
            CALL dbcsr_copy(matrix_hx(ispin)%matrix, matrix_s(1)%matrix)
            CALL dbcsr_set(matrix_hx(ispin)%matrix, 0.0_dp)
            CALL integrate_v_rspace(qs_env=qs_env, v_rspace=xv_hartree_rspace, &
                                    pmat=matrix_px1(ispin), &
                                    hmat=matrix_hx(ispin), calculate_forces=.TRUE.)
         END DO
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Px*dKh[X]   ", fodeb
         END IF
      END IF

      ! XC
      NULLIFY (fxc_rho, fxc_tau, gxc_rho, gxc_tau)
      xc_section => full_kernel%xc_section
      CALL section_vals_val_get(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                i_val=myfun)
      IF (myfun /= xc_none) THEN
         SELECT CASE (ex_env%xc_kernel_method)
         CASE (xc_kernel_method_best)
            do_analytic = .TRUE.
            do_numeric = .TRUE.
         CASE (xc_kernel_method_analytic)
            do_analytic = .TRUE.
            do_numeric = .FALSE.
         CASE (xc_kernel_method_numeric)
            do_analytic = .FALSE.
            do_numeric = .TRUE.
         CASE DEFAULT
            CPABORT("invalid xc_kernel_method")
         END SELECT

         NULLIFY (rho, ks_env)
         CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, rho=rho)
         CALL qs_rho_get(rho, rho_ao=matrix_p)
         NULLIFY (rhox)
         CALL qs_rho_create(rhox)
         CALL qs_rho_set(rho_struct=rhox, rho_ao=matrix_px1, rho_r=rhox_r, rho_g=rhox_g, &
                         rho_r_valid=.TRUE., rho_g_valid=.TRUE.)
         IF (do_analytic .AND. .NOT. do_numeric) THEN
            CPABORT("Analytic 3rd EXC derivatives not available")
         ELSE
            CALL qs_fgxc_create(ks_env, rho, rhox, xc_section, 6, is_rks_triplets, &
                                fxc_rho, fxc_tau, gxc_rho, gxc_tau)
         END IF
         CALL qs_rho_set(rho_struct=rhox, rho_ao=NULL(), rho_r=NULL(), rho_g=NULL())
         CALL qs_rho_release(rhox)
         IF (nspins == 2) THEN
            DO ispin = 1, nspins
               CALL pw_scale(gxc_rho(ispin)%pw, 0.5_dp)
               IF (ASSOCIATED(gxc_tau)) CALL pw_scale(gxc_tau(ispin)%pw, 0.5_dp)
            END DO
         END IF

         IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
         NULLIFY (matrix_fx)
         CALL dbcsr_allocate_matrix_set(matrix_fx, nspins)
         DO ispin = 1, nspins
            ALLOCATE (matrix_fx(ispin)%matrix)
            CALL dbcsr_create(matrix_fx(ispin)%matrix, template=matrix_s(1)%matrix)
            CALL dbcsr_copy(matrix_fx(ispin)%matrix, matrix_s(1)%matrix)
            CALL dbcsr_set(matrix_fx(ispin)%matrix, 0.0_dp)
            CALL pw_scale(fxc_rho(ispin)%pw, fxc_rho(ispin)%pw%pw_grid%dvol)
            CALL integrate_v_rspace(qs_env=qs_env, v_rspace=fxc_rho(ispin), &
                                    pmat=matrix_px1(ispin), &
                                    hmat=matrix_fx(ispin), calculate_forces=.TRUE.)
         END DO
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Px*dKf[X]   ", fodeb
         END IF

         IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
         NULLIFY (matrix_gx)
         CALL dbcsr_allocate_matrix_set(matrix_gx, nspins)
         DO ispin = 1, nspins
            ALLOCATE (matrix_gx(ispin)%matrix)
            CALL dbcsr_create(matrix_gx(ispin)%matrix, template=matrix_s(1)%matrix)
            CALL dbcsr_copy(matrix_gx(ispin)%matrix, matrix_s(1)%matrix)
            CALL dbcsr_set(matrix_gx(ispin)%matrix, 0.0_dp)
            CALL pw_scale(gxc_rho(ispin)%pw, gxc_rho(ispin)%pw%pw_grid%dvol)
            CALL pw_scale(gxc_rho(ispin)%pw, 0.5_dp)
            CALL integrate_v_rspace(qs_env=qs_env, v_rspace=gxc_rho(ispin), &
                                    pmat=matrix_p(ispin), &
                                    hmat=matrix_gx(ispin), calculate_forces=.TRUE.)
            CALL dbcsr_scale(matrix_gx(ispin)%matrix, 2.0_dp)
         END DO
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dKg[X]   ", fodeb
         END IF
         CALL qs_fgxc_release(ks_env, fxc_rho, fxc_tau, gxc_rho, gxc_tau)

         IF (do_admm) THEN
            IF (dft_control%admm_control%aux_exch_func == do_admm_aux_exch_func_none) THEN
               ! nothing to do
            ELSE
               IF (.NOT. tddfpt_control%admm_symm) THEN
                  CALL cp_warn(__LOCATION__, "Forces need symmetric ADMM kernel corrections")
                  CPABORT("ADMM KERNEL CORRECTION")
               END IF
               xc_section => admm_env%xc_section_aux
               CALL get_qs_env(qs_env, admm_env=admm_env, rho_aux_fit=rho_aux_fit, &
                               matrix_s_aux_fit=matrix_s_aux_fit, &
                               task_list_aux_fit=task_list_aux_fit)
               !
               NULLIFY (mfx, mgx)
               CALL dbcsr_allocate_matrix_set(mfx, nspins)
               CALL dbcsr_allocate_matrix_set(mgx, nspins)
               DO ispin = 1, nspins
                  ALLOCATE (mfx(ispin)%matrix, mgx(ispin)%matrix)
                  CALL dbcsr_create(mfx(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix)
                  CALL dbcsr_copy(mfx(ispin)%matrix, matrix_s_aux_fit(1)%matrix)
                  CALL dbcsr_set(mfx(ispin)%matrix, 0.0_dp)
                  CALL dbcsr_create(mgx(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix)
                  CALL dbcsr_copy(mgx(ispin)%matrix, matrix_s_aux_fit(1)%matrix)
                  CALL dbcsr_set(mgx(ispin)%matrix, 0.0_dp)
               END DO

               NULLIFY (rho_g_aux, rho_r_aux, tau_r_aux, rhox_g_aux, rhox_r_aux)
               CALL qs_rho_get(rho_aux_fit, rho_r=rho_r_aux, rho_g=rho_g_aux, tau_r=tau_r_aux)
               CALL qs_rho_get(rho_aux_fit, rho_ao=matrix_p_admm)
               ! rhox_aux
               ALLOCATE (rhox_r_aux(nspins), rhox_g_aux(nspins))
               DO ispin = 1, nspins
                  NULLIFY (rhox_r_aux(ispin)%pw, rhox_g_aux(ispin)%pw)
                  CALL pw_pool_create_pw(auxbas_pw_pool, rhox_r_aux(ispin)%pw, &
                                         use_data=REALDATA3D, in_space=REALSPACE)
                  CALL pw_pool_create_pw(auxbas_pw_pool, rhox_g_aux(ispin)%pw, &
                                         use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
               ENDDO
               DO ispin = 1, nspins
                  CALL calculate_rho_elec(ks_env=ks_env, matrix_p=matrix_px1_admm(ispin)%matrix, &
                                          rho=rhox_r_aux(ispin), rho_gspace=rhox_g_aux(ispin), &
                                          total_rho=total_rhox_aux, basis_type="AUX_FIT", &
                                          task_list_external=task_list_aux_fit)
               END DO
               !
               NULLIFY (rhox_aux)
               CALL qs_rho_create(rhox_aux)
               CALL qs_rho_set(rho_struct=rhox_aux, rho_ao=matrix_px1_admm, &
                               rho_r=rhox_r_aux, rho_g=rhox_g_aux, &
                               rho_r_valid=.TRUE., rho_g_valid=.TRUE.)
               IF (do_analytic .AND. .NOT. do_numeric) THEN
                  CPABORT("Analytic 3rd derivatives of EXC not available")
               ELSE
                  CALL qs_fgxc_create(ks_env, rho_aux_fit, rhox_aux, xc_section, &
                                      6, is_rks_triplets, fxc_rho, fxc_tau, gxc_rho, gxc_tau)
               END IF
               CALL qs_rho_set(rho_struct=rhox_aux, rho_ao=NULL(), rho_r=NULL(), rho_g=NULL())
               CALL qs_rho_release(rhox_aux)
               DO ispin = 1, nspins
                  CALL pw_pool_give_back_pw(auxbas_pw_pool, rhox_r_aux(ispin)%pw)
                  CALL pw_pool_give_back_pw(auxbas_pw_pool, rhox_g_aux(ispin)%pw)
               END DO
               DEALLOCATE (rhox_r_aux, rhox_g_aux)
               IF (nspins == 2) THEN
                  DO ispin = 1, nspins
                     CALL pw_scale(gxc_rho(ispin)%pw, 0.5_dp)
                     IF (ASSOCIATED(gxc_tau)) CALL pw_scale(gxc_tau(ispin)%pw, 0.5_dp)
                  END DO
               END IF
               !
               IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
               DO ispin = 1, nspins
                  CALL pw_scale(fxc_rho(ispin)%pw, fxc_rho(ispin)%pw%pw_grid%dvol)
                  CALL integrate_v_rspace(qs_env=qs_env, v_rspace=fxc_rho(ispin), &
                                          hmat=mfx(ispin), &
                                          pmat=matrix_px1_admm(ispin), &
                                          basis_type="AUX_FIT", &
                                          calculate_forces=.TRUE., &
                                          task_list_external=task_list_aux_fit)
               END DO
               IF (debug_forces) THEN
                  fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
                  CALL mp_sum(fodeb, para_env%group)
                  IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Px*dKf[X]ADMM", fodeb
               END IF

               IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
               DO ispin = 1, nspins
                  CALL pw_scale(gxc_rho(ispin)%pw, gxc_rho(ispin)%pw%pw_grid%dvol)
                  CALL pw_scale(gxc_rho(ispin)%pw, 0.5_dp)
                  CALL integrate_v_rspace(qs_env=qs_env, v_rspace=gxc_rho(ispin), &
                                          hmat=mgx(ispin), &
                                          pmat=matrix_p_admm(ispin), &
                                          basis_type="AUX_FIT", &
                                          calculate_forces=.TRUE., &
                                          task_list_external=task_list_aux_fit)
                  CALL dbcsr_scale(mgx(ispin)%matrix, 2.0_dp)
               END DO
               IF (debug_forces) THEN
                  fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
                  CALL mp_sum(fodeb, para_env%group)
                  IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dKg[X]ADMM", fodeb
               END IF
               CALL qs_fgxc_release(ks_env, fxc_rho, fxc_tau, gxc_rho, gxc_tau)
               !
               ! A' fx A - Forces
               !
               IF (debug_forces) fodeb(1:3) = force(1)%overlap_admm(1:3, 1)
               fval = 2.0_dp*REAL(nspins, KIND=dp)
               CALL admm_projection_derivative(qs_env, admm_env, mfx, matrix_px1, fval)
               IF (debug_forces) THEN
                  fodeb(1:3) = force(1)%overlap_admm(1:3, 1) - fodeb(1:3)
                  CALL mp_sum(fodeb, para_env%group)
                  IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dfXC(P)*S' ", fodeb
               END IF
               IF (debug_forces) fodeb(1:3) = force(1)%overlap_admm(1:3, 1)
               fval = 2.0_dp*REAL(nspins, KIND=dp)
               CALL admm_projection_derivative(qs_env, admm_env, mgx, matrix_p, fval)
               IF (debug_forces) THEN
                  fodeb(1:3) = force(1)%overlap_admm(1:3, 1) - fodeb(1:3)
                  CALL mp_sum(fodeb, para_env%group)
                  IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dgXC(P)*S' ", fodeb
               END IF
               !
               ! Add ADMM fx/gx to the full basis fx/gx
               nao = admm_env%nao_orb
               nao_aux = admm_env%nao_aux_fit
               ALLOCATE (dbwork)
               CALL dbcsr_create(dbwork, template=matrix_fx(1)%matrix)
               DO ispin = 1, nspins
                  ! fx
                  CALL cp_dbcsr_sm_fm_multiply(mfx(ispin)%matrix, admm_env%A, &
                                               admm_env%work_aux_orb, nao)
                  CALL cp_gemm('T', 'N', nao, nao, nao_aux, &
                               1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                               admm_env%work_orb_orb)
                  CALL dbcsr_copy(dbwork, matrix_fx(1)%matrix)
                  CALL dbcsr_set(dbwork, 0.0_dp)
                  CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbwork, keep_sparsity=.TRUE.)
                  CALL dbcsr_add(matrix_fx(ispin)%matrix, dbwork, 1.0_dp, 1.0_dp)
                  ! gx
                  CALL cp_dbcsr_sm_fm_multiply(mgx(ispin)%matrix, admm_env%A, &
                                               admm_env%work_aux_orb, nao)
                  CALL cp_gemm('T', 'N', nao, nao, nao_aux, &
                               1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                               admm_env%work_orb_orb)
                  CALL dbcsr_set(dbwork, 0.0_dp)
                  CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbwork, keep_sparsity=.TRUE.)
                  CALL dbcsr_add(matrix_gx(ispin)%matrix, dbwork, 1.0_dp, 1.0_dp)
               END DO
               CALL dbcsr_release(dbwork)
               DEALLOCATE (dbwork)
               CALL dbcsr_deallocate_matrix_set(mfx)
               CALL dbcsr_deallocate_matrix_set(mgx)
            END IF
         END IF

      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhox_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhox_r(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhox_g(ispin)%pw)
      END DO
      DEALLOCATE (rhox_r, rhox_g)
      IF (.NOT. is_rks_triplets) THEN
         CALL pw_pool_give_back_pw(auxbas_pw_pool, xv_hartree_rspace%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, xv_hartree_gspace%pw)
      END IF

      ! HFX
      IF (do_hfx) THEN
         NULLIFY (matrix_hfx)
         CALL dbcsr_allocate_matrix_set(matrix_hfx, nspins)
         DO ispin = 1, nspins
            ALLOCATE (matrix_hfx(ispin)%matrix)
            CALL dbcsr_create(matrix_hfx(ispin)%matrix, template=matrix_s(1)%matrix)
            CALL dbcsr_copy(matrix_hfx(ispin)%matrix, matrix_s(1)%matrix)
            CALL dbcsr_set(matrix_hfx(ispin)%matrix, 0.0_dp)
         END DO
         !
         xc_section => full_kernel%xc_section
         hfx_section => section_vals_get_subs_vals(xc_section, "HF")
         CALL section_vals_get(hfx_section, n_repetition=n_rep_hf)
         CPASSERT(n_rep_hf == 1)
         CALL section_vals_val_get(hfx_section, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core, &
                                   i_rep_section=1)
         mspin = 1
         IF (hfx_treat_lsd_in_core) mspin = nspins
         !
         CALL get_qs_env(qs_env=qs_env, x_data=x_data, s_mstruct_changed=s_mstruct_changed)
         distribute_fock_matrix = .TRUE.
         !
         IF (do_admm) THEN
            CALL get_qs_env(qs_env=qs_env, matrix_s_aux_fit=matrix_s_aux_fit)
            NULLIFY (matrix_hfx_admm)
            CALL dbcsr_allocate_matrix_set(matrix_hfx_admm, nspins)
            DO ispin = 1, nspins
               ALLOCATE (matrix_hfx_admm(ispin)%matrix)
               CALL dbcsr_create(matrix_hfx_admm(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix)
               CALL dbcsr_copy(matrix_hfx_admm(ispin)%matrix, matrix_s_aux_fit(1)%matrix)
               CALL dbcsr_set(matrix_hfx_admm(ispin)%matrix, 0.0_dp)
            END DO
            !
            NULLIFY (mpe, mhe)
            ALLOCATE (mpe(nspins, 1), mhe(nspins, 1))
            DO ispin = 1, nspins
               mhe(ispin, 1)%matrix => matrix_hfx_admm(ispin)%matrix
               mpe(ispin, 1)%matrix => matrix_px1_admm(ispin)%matrix
            END DO
            DO ispin = 1, mspin
               eh1 = 0.0
               CALL integrate_four_center(qs_env, x_data, mhe, eh1, mpe, hfx_section, &
                                          para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                          ispin=ispin)
            END DO
            !
            nao = admm_env%nao_orb
            nao_aux = admm_env%nao_aux_fit
            ALLOCATE (dbwork)
            CALL dbcsr_create(dbwork, template=matrix_hfx(1)%matrix)
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(matrix_hfx_admm(ispin)%matrix, admm_env%A, &
                                            admm_env%work_aux_orb, nao)
               CALL cp_gemm('T', 'N', nao, nao, nao_aux, &
                            1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                            admm_env%work_orb_orb)
               CALL dbcsr_copy(dbwork, matrix_hfx(1)%matrix)
               CALL dbcsr_set(dbwork, 0.0_dp)
               CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbwork, keep_sparsity=.TRUE.)
               CALL dbcsr_add(matrix_hfx(ispin)%matrix, dbwork, 1.0_dp, 1.0_dp)
            END DO
            CALL dbcsr_release(dbwork)
            DEALLOCATE (dbwork)
            ! forces
            ! ADMM Projection force
            IF (debug_forces) fodeb(1:3) = force(1)%overlap_admm(1:3, 1)
            fval = 4.0_dp*REAL(nspins, KIND=dp)
            CALL admm_projection_derivative(qs_env, admm_env, matrix_hfx_admm, matrix_px1, fval)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%overlap_admm(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*Hx(P)*S' ", fodeb
            END IF
            !
            use_virial = .FALSE.
            NULLIFY (mdum)
            fval = 2.0_dp*REAL(nspins, KIND=dp)
            IF (debug_forces) fodeb(1:3) = force(1)%fock_4c(1:3, 1)
            CALL derivatives_four_center(qs_env, mpe, mdum, hfx_section, para_env, 1, use_virial, &
                                         adiabatic_rescale_factor=fval)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%fock_4c(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Px*hfx'*Px ", fodeb
            END IF
            !
            DEALLOCATE (mpe, mhe)
            !
            CALL dbcsr_deallocate_matrix_set(matrix_hfx_admm)
         ELSE
            NULLIFY (mpe, mhe)
            ALLOCATE (mpe(nspins, 1), mhe(nspins, 1))
            DO ispin = 1, nspins
               mhe(ispin, 1)%matrix => matrix_hfx(ispin)%matrix
               mpe(ispin, 1)%matrix => matrix_px1(ispin)%matrix
            END DO
            DO ispin = 1, mspin
               eh1 = 0.0
               CALL integrate_four_center(qs_env, x_data, mhe, eh1, mpe, hfx_section, &
                                          para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                          ispin=ispin)
            END DO
            ! forces
            use_virial = .FALSE.
            NULLIFY (mdum)
            fval = 2.0_dp*REAL(nspins, KIND=dp)
            IF (debug_forces) fodeb(1:3) = force(1)%fock_4c(1:3, 1)
            CALL derivatives_four_center(qs_env, mpe, mdum, hfx_section, para_env, 1, use_virial, &
                                         adiabatic_rescale_factor=fval)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%fock_4c(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Px*hfx'*Px ", fodeb
            END IF
            !
            DEALLOCATE (mpe, mhe)
         END IF
         fval = 2.0_dp*REAL(nspins, KIND=dp)
         DO ispin = 1, nspins
            CALL dbcsr_scale(matrix_hfx(ispin)%matrix, fval)
         END DO
      END IF

      CALL get_qs_env(qs_env, sab_orb=sab_orb)
      NULLIFY (matrix_wx1)
      CALL dbcsr_allocate_matrix_set(matrix_wx1, nspins)
      cpmos => ex_env%cpmos
      focc = 2.0_dp
      IF (nspins == 2) focc = 1.0_dp
      DO ispin = 1, nspins
         mos => gs_mos(ispin)%mos_occ
         CALL cp_fm_get_info(evect(ispin)%matrix, ncol_global=norb)
         CALL cp_fm_create(vcvec, mos%matrix_struct, "vcvec")
         CALL cp_fm_get_info(vcvec, matrix_struct=fm_struct, nrow_global=nao)
         CALL cp_fm_struct_create(fm_struct_mat, context=fm_struct%context, nrow_global=norb, &
                                  ncol_global=norb, para_env=fm_struct%para_env)
         CALL cp_fm_create(cvcmat, fm_struct_mat)
         CALL cp_fm_struct_release(fm_struct_mat)
         !
         ALLOCATE (matrix_wx1(ispin)%matrix)
         CALL dbcsr_create(matrix=matrix_wx1(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_wx1(ispin)%matrix, sab_orb)
         CALL dbcsr_set(matrix_wx1(ispin)%matrix, 0.0_dp)
         !
         IF (.NOT. is_rks_triplets) THEN
            CALL cp_dbcsr_sm_fm_multiply(matrix_hx(ispin)%matrix, evect(ispin)%matrix, &
                                         cpmos(ispin)%matrix, norb, alpha=focc, beta=1.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(matrix_hx(ispin)%matrix, mos, vcvec, norb, alpha=1.0_dp, beta=0.0_dp)
            CALL cp_gemm("T", "N", norb, norb, nao, 1.0_dp, mos, vcvec, 0.0_dp, cvcmat)
            CALL cp_gemm("N", "N", nao, norb, norb, 1.0_dp, evect(ispin)%matrix, cvcmat, 0.0_dp, vcvec)
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, vcvec, cpmos(ispin)%matrix, &
                                         norb, alpha=-focc, beta=1.0_dp)
            !
            CALL cp_dbcsr_plus_fm_fm_t(matrix_wx1(ispin)%matrix, matrix_v=mos, matrix_g=vcvec, &
                                       ncol=norb, alpha=1.0_dp)
            CALL cp_dbcsr_plus_fm_fm_t(matrix_wx1(ispin)%matrix, matrix_v=vcvec, matrix_g=mos, &
                                       ncol=norb, alpha=1.0_dp)
         END IF
         !
         IF (myfun /= xc_none) THEN
            CALL cp_dbcsr_sm_fm_multiply(matrix_fx(ispin)%matrix, evect(ispin)%matrix, &
                                         cpmos(ispin)%matrix, norb, alpha=focc, beta=1.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(matrix_fx(ispin)%matrix, mos, vcvec, norb, alpha=1.0_dp, beta=0.0_dp)
            CALL cp_gemm("T", "N", norb, norb, nao, 1.0_dp, mos, vcvec, 0.0_dp, cvcmat)
            CALL cp_gemm("N", "N", nao, norb, norb, 1.0_dp, evect(ispin)%matrix, cvcmat, 0.0_dp, vcvec)
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, vcvec, cpmos(ispin)%matrix, &
                                         norb, alpha=-focc, beta=1.0_dp)
            !
            CALL cp_dbcsr_plus_fm_fm_t(matrix_wx1(ispin)%matrix, matrix_v=mos, matrix_g=vcvec, &
                                       ncol=norb, alpha=1.0_dp)
            CALL cp_dbcsr_plus_fm_fm_t(matrix_wx1(ispin)%matrix, matrix_v=vcvec, matrix_g=mos, &
                                       ncol=norb, alpha=1.0_dp)
            !
            CALL cp_dbcsr_sm_fm_multiply(matrix_gx(ispin)%matrix, mos, &
                                         cpmos(ispin)%matrix, norb, alpha=focc, beta=1.0_dp)
            !
            CALL cp_dbcsr_sm_fm_multiply(matrix_gx(ispin)%matrix, mos, vcvec, norb, alpha=1.0_dp, beta=0.0_dp)
            CALL cp_gemm("T", "N", norb, norb, nao, 1.0_dp, mos, vcvec, 0.0_dp, cvcmat)
            CALL cp_gemm("N", "N", nao, norb, norb, 1.0_dp, mos, cvcmat, 0.0_dp, vcvec)
            CALL cp_dbcsr_plus_fm_fm_t(matrix_wx1(ispin)%matrix, matrix_v=mos, matrix_g=vcvec, &
                                       ncol=norb, alpha=0.5_dp)
            CALL cp_dbcsr_plus_fm_fm_t(matrix_wx1(ispin)%matrix, matrix_v=vcvec, matrix_g=mos, &
                                       ncol=norb, alpha=0.5_dp)
         END IF
         !
         IF (do_hfx) THEN
            CALL cp_dbcsr_sm_fm_multiply(matrix_hfx(ispin)%matrix, evect(ispin)%matrix, &
                                         cpmos(ispin)%matrix, norb, alpha=focc, beta=1.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(matrix_hfx(ispin)%matrix, mos, vcvec, norb, alpha=1.0_dp, beta=0.0_dp)
            CALL cp_gemm("T", "N", norb, norb, nao, 1.0_dp, mos, vcvec, 0.0_dp, cvcmat)
            CALL cp_gemm("N", "N", nao, norb, norb, 1.0_dp, evect(ispin)%matrix, cvcmat, 0.0_dp, vcvec)
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, vcvec, cpmos(ispin)%matrix, &
                                         norb, alpha=-focc, beta=1.0_dp)
            !
            CALL cp_dbcsr_plus_fm_fm_t(matrix_wx1(ispin)%matrix, matrix_v=mos, matrix_g=vcvec, &
                                       ncol=norb, alpha=1.0_dp)
            CALL cp_dbcsr_plus_fm_fm_t(matrix_wx1(ispin)%matrix, matrix_v=vcvec, matrix_g=mos, &
                                       ncol=norb, alpha=1.0_dp)
         END IF
         !
         CALL cp_fm_release(vcvec)
         CALL cp_fm_release(cvcmat)
      END DO

      IF (.NOT. is_rks_triplets) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_hx)
      END IF
      IF (ASSOCIATED(ex_env%matrix_wx1)) CALL dbcsr_deallocate_matrix_set(ex_env%matrix_wx1)
      ex_env%matrix_wx1 => matrix_wx1
      IF (myfun /= xc_none) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_fx)
         CALL dbcsr_deallocate_matrix_set(matrix_gx)
      END IF
      IF (do_hfx) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_hfx)
      END IF

      CALL timestop(handle)

   END SUBROUTINE fhxc_force

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param admm_env ...
!> \param matrix_hfx ...
!> \param matrix_pe ...
!> \param fval ...
! **************************************************************************************************
   SUBROUTINE admm_projection_derivative(qs_env, admm_env, matrix_hfx, matrix_pe, fval)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hfx, matrix_pe
      REAL(KIND=dp), INTENT(IN)                          :: fval

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
      CPASSERT(ASSOCIATED(matrix_hfx))
      CPASSERT(ASSOCIATED(matrix_pe))

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

      nspins = SIZE(matrix_pe)
      nao = admm_env%nao_orb
      naux = admm_env%nao_aux_fit
      DO ispin = 1, nspins
         CALL copy_dbcsr_to_fm(matrix_hfx(ispin)%matrix, admm_env%work_aux_aux)

         CALL cp_gemm("N", "N", naux, naux, naux, 1.0_dp, admm_env%s_inv, &
                      admm_env%work_aux_aux, 0.0_dp, admm_env%work_aux_aux2)
         CALL cp_gemm("N", "N", naux, nao, naux, 1.0_dp, admm_env%work_aux_aux2, &
                      admm_env%A, 0.0_dp, admm_env%work_aux_orb)
         CALL copy_dbcsr_to_fm(matrix_pe(ispin)%matrix, admm_env%work_orb_orb)
         ! admm_env%work_aux_orb2 = S-1*H(P)*A*P
         CALL cp_gemm("N", "N", naux, nao, nao, 1.0_dp, admm_env%work_aux_orb, &
                      admm_env%work_orb_orb, 0.0_dp, admm_env%work_aux_orb2)
         ! admm_env%work_aux_aux = [S-1*H(P)*A*P]*A(T)
         CALL cp_gemm("N", "T", naux, naux, nao, 1.0_dp, admm_env%work_aux_orb2, &
                      admm_env%A, 0.0_dp, admm_env%work_aux_aux)
         !
         CALL copy_fm_to_dbcsr(admm_env%work_aux_orb2, matrix_w_q, keep_sparsity=.TRUE.)
         CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, matrix_w_s, keep_sparsity=.TRUE.)
         !
         CALL dbcsr_scale(matrix_w_q, -fval)
         CALL dbcsr_scale(matrix_w_s, fval)
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
!> \brief Simplified Tamm Dancoff approach (sTDA). Kernel contribution to forces
!> \param qs_env ...
!> \param ex_env ...
!> \param gs_mos ...
!> \param stda_env ...
!> \param sub_env ...
!> \param work ...
!> \param debug_forces ...
! **************************************************************************************************
   SUBROUTINE stda_force(qs_env, ex_env, gs_mos, stda_env, sub_env, work, debug_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(stda_env_type), POINTER                       :: stda_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work
      LOGICAL, INTENT(IN)                                :: debug_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'stda_force', routineP = moduleN//':'//routineN

      INTEGER :: atom_i, atom_j, blk, ewald_type, handle, i, ia, iatom, idimk, ikind, iounit, is, &
         ispin, jatom, jkind, jspin, nao, natom, norb, nsgf, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, first_sgf, kind_of, &
                                                            last_sgf
      INTEGER, DIMENSION(2)                              :: nactive, nlim
      LOGICAL                                            :: calculate_forces, do_coulomb, do_ewald, &
                                                            found, is_rks_triplets, use_virial
      REAL(KIND=dp)                                      :: alpha, bp, dgabr, dr, eta, fdim, gabr, &
                                                            hfx, rbeta, spinfac, xfac
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tcharge, tv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gtcharge
      REAL(KIND=dp), DIMENSION(3)                        :: fij, fodeb, rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gab, pblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: cpmos, X, xtransformed
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct, fmstruct_mat, fmstructjspin
      TYPE(cp_fm_type), POINTER                          :: ct, ctjspin, cvcmat, cvec, cvecjspin, &
                                                            t0matrix, t1matrix, ucmatrix, &
                                                            uxmatrix, vcvec, xt, xvec
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gamma_matrix, matrix_plo, matrix_s, &
                                                            matrix_wx1, scrm
      TYPE(dbcsr_type)                                   :: pdens, ptrans
      TYPE(dbcsr_type), POINTER                          :: tempmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list, sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(stda_control_type), POINTER                   :: stda_control
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(ex_env))
      CPASSERT(ASSOCIATED(gs_mos))

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      ENDIF

      CALL get_qs_env(qs_env, dft_control=dft_control)
      tddfpt_control => dft_control%tddfpt2_control
      stda_control => tddfpt_control%stda_control
      nspins = dft_control%nspins
      is_rks_triplets = tddfpt_control%rks_triplets .AND. (nspins == 1)

      X => ex_env%evect

      nactive(:) = stda_env%nactive(:)
      xfac = 2.0_dp
      spinfac = 2.0_dp
      IF (nspins == 2) spinfac = 1.0_dp
      NULLIFY (matrix_wx1)
      CALL dbcsr_allocate_matrix_set(matrix_wx1, nspins)
      NULLIFY (matrix_plo)
      CALL dbcsr_allocate_matrix_set(matrix_plo, nspins)

      IF (nspins == 1 .AND. is_rks_triplets) THEN
         do_coulomb = .FALSE.
      ELSE
         do_coulomb = .TRUE.
      END IF
      do_ewald = stda_control%do_ewald

      CALL get_qs_env(qs_env, para_env=para_env, force=force)

      CALL get_qs_env(qs_env, natom=natom, cell=cell, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set)
      ALLOCATE (first_sgf(natom))
      ALLOCATE (last_sgf(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf)

      CALL get_qs_env(qs_env, ks_env=ks_env, matrix_s=matrix_s, sab_orb=sab_orb)

      ! calculate Loewdin transformed Davidson trial vector tilde(X)=S^1/2*X
      ! and tilde(tilde(X))=S^1/2_A*tilde(X)_A
      ALLOCATE (xtransformed(nspins))
      DO ispin = 1, nspins
         NULLIFY (fmstruct)
         ct => work%ctransformed(ispin)%matrix
         CALL cp_fm_get_info(ct, matrix_struct=fmstruct)
         CALL cp_fm_create(matrix=xtransformed(ispin)%matrix, matrix_struct=fmstruct, name="XTRANSFORMED")
      ENDDO
      CALL get_lowdin_x(work%shalf, X, xtransformed)

      ALLOCATE (tcharge(natom), gtcharge(natom, 4))

      cpmos => ex_env%cpmos

      DO ispin = 1, nspins
         ct => work%ctransformed(ispin)%matrix
         CALL cp_fm_get_info(ct, matrix_struct=fmstruct, nrow_global=nsgf)
         ALLOCATE (tv(nsgf))
         CALL cp_fm_create(cvec, fmstruct)
         CALL cp_fm_create(xvec, fmstruct)
         !
         ALLOCATE (matrix_wx1(ispin)%matrix)
         CALL dbcsr_create(matrix=matrix_wx1(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_wx1(ispin)%matrix, sab_orb)
         CALL dbcsr_set(matrix_wx1(ispin)%matrix, 0.0_dp)
         ALLOCATE (matrix_plo(ispin)%matrix)
         CALL dbcsr_create(matrix=matrix_plo(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_plo(ispin)%matrix, sab_orb)
         CALL dbcsr_set(matrix_plo(ispin)%matrix, 0.0_dp)
         !
         ! *** Coulomb contribution
         !
         IF (do_coulomb) THEN
            !
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            !
            ALLOCATE (kind_of(natom), atom_of_kind(natom))
            CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
            CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of, &
                                     atom_of_kind=atom_of_kind)
            !
            tcharge(:) = 0.0_dp
            DO jspin = 1, nspins
               ctjspin => work%ctransformed(jspin)%matrix
               xt => xtransformed(jspin)%matrix
               CALL cp_fm_get_info(ctjspin, matrix_struct=fmstructjspin)
               CALL cp_fm_get_info(ctjspin, matrix_struct=fmstructjspin, nrow_global=nsgf)
               CALL cp_fm_create(cvecjspin, fmstructjspin)
               ! CV(mu,j) = CT(mu,j)*XT(mu,j)
               CALL cp_fm_schur_product(ctjspin, xt, cvecjspin)
               ! TV(mu) = SUM_j CV(mu,j)
               CALL cp_fm_vectorssum(cvecjspin, tv, "R")
               ! contract charges
               ! TC(a) = SUM_(mu of a) TV(mu)
               DO ia = 1, natom
                  DO is = first_sgf(ia), last_sgf(ia)
                     tcharge(ia) = tcharge(ia) + tv(is)
                  END DO
               END DO
               CALL cp_fm_release(cvecjspin)
            ENDDO !jspin
            ! Apply tcharge*gab -> gtcharge
            ! gT(b) = SUM_a g(a,b)*TC(a)
            ! gab = work%gamma_exchange(1)%matrix
            gtcharge = 0.0_dp
            ! short range contribution
            NULLIFY (gamma_matrix)
            CALL setup_gamma(qs_env, stda_env, sub_env, gamma_matrix, ndim=4)
            tempmat => gamma_matrix(1)%matrix
            CALL dbcsr_iterator_start(iter, tempmat)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, iatom, jatom, gab, blk)
               gtcharge(iatom, 1) = gtcharge(iatom, 1) + gab(1, 1)*tcharge(jatom)
               IF (iatom /= jatom) THEN
                  gtcharge(jatom, 1) = gtcharge(jatom, 1) + gab(1, 1)*tcharge(iatom)
               END IF
               DO idimk = 2, 4
                  fdim = -1.0_dp
                  CALL dbcsr_get_block_p(matrix=gamma_matrix(idimk)%matrix, &
                                         row=iatom, col=jatom, block=gab, found=found)
                  IF (found) THEN
                     gtcharge(iatom, idimk) = gtcharge(iatom, idimk) + gab(1, 1)*tcharge(jatom)
                     IF (iatom /= jatom) THEN
                        gtcharge(jatom, idimk) = gtcharge(jatom, idimk) + fdim*gab(1, 1)*tcharge(iatom)
                     END IF
                  END IF
               END DO
            END DO
            CALL dbcsr_iterator_stop(iter)
            CALL dbcsr_deallocate_matrix_set(gamma_matrix)
            ! Ewald long range contribution
            IF (do_ewald) THEN
               ewald_env => work%ewald_env
               ewald_pw => work%ewald_pw
               CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
               CALL get_qs_env(qs_env=qs_env, sab_orb=n_list, virial=virial, atprop=atprop)
               use_virial = .FALSE.
               calculate_forces = .FALSE.
               CALL tb_ewald_overlap(gtcharge, tcharge, alpha, n_list, virial, use_virial, atprop)
               CALL tb_spme_evaluate(ewald_env, ewald_pw, particle_set, cell, &
                                     gtcharge, tcharge, calculate_forces, virial, use_virial, atprop)
               ! add self charge interaction contribution
               IF (para_env%ionode) THEN
                  gtcharge(:, 1) = gtcharge(:, 1) - 2._dp*alpha*oorootpi*tcharge(:)
               END IF
            ELSE
               nlim = get_limit(natom, para_env%num_pe, para_env%mepos)
               DO iatom = nlim(1), nlim(2)
                  DO jatom = 1, iatom - 1
                     rij = particle_set(iatom)%r - particle_set(jatom)%r
                     rij = pbc(rij, cell)
                     dr = SQRT(SUM(rij(:)**2))
                     IF (dr > 1.e-6_dp) THEN
                        gtcharge(iatom, 1) = gtcharge(iatom, 1) + tcharge(jatom)/dr
                        gtcharge(jatom, 1) = gtcharge(jatom, 1) + tcharge(iatom)/dr
                        DO idimk = 2, 4
                           gtcharge(iatom, idimk) = gtcharge(iatom, idimk) + rij(idimk - 1)*tcharge(jatom)/dr**3
                           gtcharge(jatom, idimk) = gtcharge(jatom, idimk) - rij(idimk - 1)*tcharge(iatom)/dr**3
                        END DO
                     END IF
                  END DO
               END DO
            END IF
            CALL mp_sum(gtcharge(:, 1), para_env%group)
            ! expand charges
            ! TV(mu) = TC(a of mu)
            tv(1:nsgf) = 0.0_dp
            DO ia = 1, natom
               DO is = first_sgf(ia), last_sgf(ia)
                  tv(is) = gtcharge(ia, 1)
               END DO
            END DO
            !
            DO iatom = 1, natom
               ikind = kind_of(iatom)
               atom_i = atom_of_kind(iatom)
               DO i = 1, 3
                  fij(i) = spinfac*spinfac*gtcharge(iatom, i + 1)*tcharge(iatom)
               END DO
               force(ikind)%rho_elec(1, atom_i) = force(ikind)%rho_elec(1, atom_i) - fij(1)
               force(ikind)%rho_elec(2, atom_i) = force(ikind)%rho_elec(2, atom_i) - fij(2)
               force(ikind)%rho_elec(3, atom_i) = force(ikind)%rho_elec(3, atom_i) - fij(3)
            END DO
            !
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Coul[X]   ", fodeb
            END IF
            norb = nactive(ispin)
            ! forces from Lowdin charge derivative
            CALL cp_fm_get_info(work%S_C0_C0T(ispin)%matrix, matrix_struct=fmstruct)
            CALL cp_fm_create(t0matrix, matrix_struct=fmstruct, name="T0 SCRATCH")
            CALL cp_fm_create(t1matrix, matrix_struct=fmstruct, name="T1 SCRATCH")
            CALL fm_pool_create_fm(work%fm_pool_ao_mo_occ(ispin)%pool, ucmatrix)
            CALL fm_pool_create_fm(work%fm_pool_ao_mo_occ(ispin)%pool, uxmatrix)
            ct => work%ctransformed(ispin)%matrix
            CALL cp_fm_to_fm(ct, cvec)
            CALL cp_fm_row_scale(cvec, tv)
            CALL cp_gemm('T', 'N', nsgf, norb, nsgf, 1.0_dp, work%S_eigenvectors, &
                         cvec, 0.0_dp, ucmatrix)
            CALL cp_gemm('T', 'N', nsgf, norb, nsgf, 1.0_dp, work%S_eigenvectors, &
                         X(ispin)%matrix, 0.0_dp, uxmatrix)
            CALL cp_gemm('N', 'T', nsgf, nsgf, norb, 1.0_dp, uxmatrix, ucmatrix, 0.0_dp, t0matrix)
            xt => xtransformed(ispin)%matrix
            CALL cp_fm_to_fm(xt, cvec)
            CALL cp_fm_row_scale(cvec, tv)
            CALL cp_gemm('T', 'N', nsgf, norb, nsgf, 1.0_dp, work%S_eigenvectors, &
                         cvec, 0.0_dp, uxmatrix)
            CALL cp_gemm('T', 'N', nsgf, norb, nsgf, 1.0_dp, work%S_eigenvectors, &
                         gs_mos(ispin)%mos_occ, 0.0_dp, ucmatrix)
            CALL cp_gemm('N', 'T', nsgf, nsgf, norb, 1.0_dp, ucmatrix, uxmatrix, 1.0_dp, t0matrix)
            CALL cp_fm_schur_product(work%slambda, t0matrix, t1matrix)
            !
            CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, spinfac, work%S_eigenvectors, t1matrix, &
                         0.0_dp, t0matrix)
            CALL cp_dbcsr_plus_fm_fm_t(matrix_plo(ispin)%matrix, matrix_v=t0matrix, &
                                       matrix_g=work%S_eigenvectors, ncol=nsgf, alpha=1.0_dp)
            CALL cp_dbcsr_plus_fm_fm_t(matrix_plo(ispin)%matrix, matrix_v=work%S_eigenvectors, &
                                       matrix_g=t0matrix, ncol=nsgf, alpha=1.0_dp)
            CALL fm_pool_give_back_fm(work%fm_pool_ao_mo_occ(ispin)%pool, ucmatrix)
            CALL fm_pool_give_back_fm(work%fm_pool_ao_mo_occ(ispin)%pool, uxmatrix)
            CALL cp_fm_release(t0matrix)
            CALL cp_fm_release(t1matrix)
            !
            ! CV(mu,i) = TV(mu)*XT(mu,i)
            xt => xtransformed(ispin)%matrix
            CALL cp_fm_to_fm(xt, cvec)
            CALL cp_fm_row_scale(cvec, tv)
            CALL cp_dbcsr_sm_fm_multiply(work%shalf, cvec, cpmos(ispin)%matrix, norb, 2.0_dp*spinfac, 1.0_dp)
            ! CV(mu,i) = TV(mu)*CT(mu,i)
            ct => work%ctransformed(ispin)%matrix
            CALL cp_fm_to_fm(ct, cvec)
            CALL cp_fm_row_scale(cvec, tv)
            ! Shalf(nu,mu)*CV(mu,i)
            CALL cp_fm_get_info(cvec, matrix_struct=fmstruct, nrow_global=nao)
            CALL cp_fm_create(vcvec, fmstruct)
            CALL cp_dbcsr_sm_fm_multiply(work%shalf, cvec, vcvec, norb, 1.0_dp, 0.0_dp)
            CALL cp_fm_struct_create(fmstruct_mat, context=fmstruct%context, nrow_global=norb, &
                                     ncol_global=norb, para_env=fmstruct%para_env)
            CALL cp_fm_create(cvcmat, fmstruct_mat)
            CALL cp_fm_struct_release(fmstruct_mat)
            CALL cp_gemm("T", "N", norb, norb, nao, 1.0_dp, gs_mos(ispin)%mos_occ, vcvec, 0.0_dp, cvcmat)
            CALL cp_gemm("N", "N", nao, norb, norb, 1.0_dp, X(ispin)%matrix, cvcmat, 0.0_dp, vcvec)
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, vcvec, cpmos(ispin)%matrix, &
                                         nactive(ispin), alpha=-2.0_dp*spinfac, beta=1.0_dp)
            ! wx1
            alpha = 2.0_dp
            CALL cp_dbcsr_plus_fm_fm_t(matrix_wx1(ispin)%matrix, matrix_v=gs_mos(ispin)%mos_occ, &
                                       matrix_g=vcvec, ncol=norb, alpha=alpha)
            CALL cp_dbcsr_plus_fm_fm_t(matrix_wx1(ispin)%matrix, matrix_v=vcvec, &
                                       matrix_g=gs_mos(ispin)%mos_occ, ncol=norb, alpha=alpha)
            CALL cp_fm_release(vcvec)
            CALL cp_fm_release(cvcmat)
            DEALLOCATE (atom_of_kind, kind_of)
         END IF
         !
         ! *** Exchange contribution
         !
         IF (stda_env%hfx_fraction > 0.0_dp) THEN
            !
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            !
            norb = nactive(ispin)
            !
            ALLOCATE (kind_of(natom), atom_of_kind(natom))
            CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
            CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of, &
                                     atom_of_kind=atom_of_kind)
            !
            tempmat => work%shalf
            CALL dbcsr_create(pdens, template=tempmat, matrix_type=dbcsr_type_no_symmetry)
            ! P(nu,mu) = SUM_j XT(nu,j)*CT(mu,j)
            ct => work%ctransformed(ispin)%matrix
            xt => xtransformed(ispin)%matrix
            CALL dbcsr_set(pdens, 0.0_dp)
            CALL cp_dbcsr_plus_fm_fm_t(pdens, xt, ct, nactive(ispin), &
                                       1.0_dp, keep_sparsity=.FALSE.)
            CALL dbcsr_filter(pdens, stda_env%eps_td_filter)
            ! Apply PP*gab -> PP; gab = gamma_coulomb
            ! P(nu,mu) = P(nu,mu)*g(a of nu,b of mu)
            bp = stda_env%beta_param
            hfx = stda_env%hfx_fraction
            CALL dbcsr_iterator_start(iter, pdens)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, iatom, jatom, pblock, blk)
               rij = particle_set(iatom)%r - particle_set(jatom)%r
               rij = pbc(rij, cell)
               dr = SQRT(SUM(rij(:)**2))
               ikind = kind_of(iatom)
               jkind = kind_of(jatom)
               eta = (stda_env%kind_param_set(ikind)%kind_param%hardness_param + &
                      stda_env%kind_param_set(jkind)%kind_param%hardness_param)/2.0_dp
               rbeta = dr**bp
               gabr = (1._dp/(rbeta + (hfx*eta)**(-bp)))**(1._dp/bp)
               ! forces
               IF (dr > 1.0e-6_dp) THEN
                  dgabr = -(1._dp/bp)*(1._dp/(rbeta + (hfx*eta)**(-bp)))**(1._dp/bp + 1._dp)
                  dgabr = bp*rbeta/dr**2*dgabr
                  dgabr = SUM(pblock**2)*dgabr
                  atom_i = atom_of_kind(iatom)
                  atom_j = atom_of_kind(jatom)
                  DO i = 1, 3
                     fij(i) = dgabr*rij(i)
                  END DO
                  force(ikind)%rho_elec(1, atom_i) = force(ikind)%rho_elec(1, atom_i) - fij(1)
                  force(ikind)%rho_elec(2, atom_i) = force(ikind)%rho_elec(2, atom_i) - fij(2)
                  force(ikind)%rho_elec(3, atom_i) = force(ikind)%rho_elec(3, atom_i) - fij(3)
                  force(jkind)%rho_elec(1, atom_j) = force(jkind)%rho_elec(1, atom_j) + fij(1)
                  force(jkind)%rho_elec(2, atom_j) = force(jkind)%rho_elec(2, atom_j) + fij(2)
                  force(jkind)%rho_elec(3, atom_j) = force(jkind)%rho_elec(3, atom_j) + fij(3)
               END IF
               !
               pblock = gabr*pblock
            END DO
            CALL dbcsr_iterator_stop(iter)
            !
            ! Transpose pdens matrix
            CALL dbcsr_create(ptrans, template=pdens)
            CALL dbcsr_transposed(ptrans, pdens)
            !
            ! forces from Lowdin charge derivative
            CALL cp_fm_get_info(work%S_C0_C0T(ispin)%matrix, matrix_struct=fmstruct)
            CALL cp_fm_create(t0matrix, matrix_struct=fmstruct, name="T0 SCRATCH")
            CALL cp_fm_create(t1matrix, matrix_struct=fmstruct, name="T1 SCRATCH")
            CALL fm_pool_create_fm(work%fm_pool_ao_mo_occ(ispin)%pool, ucmatrix)
            CALL fm_pool_create_fm(work%fm_pool_ao_mo_occ(ispin)%pool, uxmatrix)
            ct => work%ctransformed(ispin)%matrix
            CALL cp_dbcsr_sm_fm_multiply(pdens, ct, cvec, norb, 1.0_dp, 0.0_dp)
            CALL cp_gemm('T', 'N', nsgf, norb, nsgf, 1.0_dp, work%S_eigenvectors, &
                         cvec, 0.0_dp, ucmatrix)
            CALL cp_gemm('T', 'N', nsgf, norb, nsgf, 1.0_dp, work%S_eigenvectors, &
                         X(ispin)%matrix, 0.0_dp, uxmatrix)
            CALL cp_gemm('N', 'T', nsgf, nsgf, norb, 1.0_dp, uxmatrix, ucmatrix, 0.0_dp, t0matrix)
            xt => xtransformed(ispin)%matrix
            CALL cp_dbcsr_sm_fm_multiply(ptrans, xt, cvec, norb, 1.0_dp, 0.0_dp)
            CALL cp_gemm('T', 'N', nsgf, norb, nsgf, 1.0_dp, work%S_eigenvectors, &
                         cvec, 0.0_dp, uxmatrix)
            CALL cp_gemm('T', 'N', nsgf, norb, nsgf, 1.0_dp, work%S_eigenvectors, &
                         gs_mos(ispin)%mos_occ, 0.0_dp, ucmatrix)
            CALL cp_gemm('N', 'T', nsgf, nsgf, norb, 1.0_dp, ucmatrix, uxmatrix, 1.0_dp, t0matrix)
            CALL cp_fm_schur_product(work%slambda, t0matrix, t1matrix)
            CALL cp_gemm('N', 'N', nsgf, nsgf, nsgf, -1.0_dp, work%S_eigenvectors, t1matrix, &
                         0.0_dp, t0matrix)
            CALL cp_dbcsr_plus_fm_fm_t(matrix_plo(ispin)%matrix, matrix_v=t0matrix, &
                                       matrix_g=work%S_eigenvectors, ncol=nsgf, alpha=1.0_dp)
            CALL cp_dbcsr_plus_fm_fm_t(matrix_plo(ispin)%matrix, matrix_v=work%S_eigenvectors, &
                                       matrix_g=t0matrix, ncol=nsgf, alpha=1.0_dp)
            CALL fm_pool_give_back_fm(work%fm_pool_ao_mo_occ(ispin)%pool, ucmatrix)
            CALL fm_pool_give_back_fm(work%fm_pool_ao_mo_occ(ispin)%pool, uxmatrix)
            CALL cp_fm_release(t0matrix)
            CALL cp_fm_release(t1matrix)

            ! RHS contribution to response matrix
            ! CV(nu,i) = P(nu,mu)*XT(mu,i)
            CALL cp_dbcsr_sm_fm_multiply(ptrans, xt, cvec, norb, 1.0_dp, 0.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(work%shalf, cvec, cpmos(ispin)%matrix, norb, &
                                         alpha=-xfac, beta=1.0_dp)
            !
            CALL cp_fm_get_info(cvec, matrix_struct=fmstruct, nrow_global=nao)
            CALL cp_fm_create(vcvec, fmstruct)
            ! CV(nu,i) = P(nu,mu)*CT(mu,i)
            CALL cp_dbcsr_sm_fm_multiply(ptrans, ct, cvec, norb, 1.0_dp, 0.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(work%shalf, cvec, vcvec, norb, 1.0_dp, 0.0_dp)
            CALL cp_fm_struct_create(fmstruct_mat, context=fmstruct%context, nrow_global=norb, &
                                     ncol_global=norb, para_env=fmstruct%para_env)
            CALL cp_fm_create(cvcmat, fmstruct_mat)
            CALL cp_fm_struct_release(fmstruct_mat)
            CALL cp_gemm("T", "N", norb, norb, nao, 1.0_dp, gs_mos(ispin)%mos_occ, vcvec, 0.0_dp, cvcmat)
            CALL cp_gemm("N", "N", nao, norb, norb, 1.0_dp, X(ispin)%matrix, cvcmat, 0.0_dp, vcvec)
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, vcvec, cpmos(ispin)%matrix, &
                                         norb, alpha=xfac, beta=1.0_dp)
            ! wx1
            IF (nspins == 2) THEN
               alpha = -2.0_dp
            ELSE
               alpha = -1.0_dp
            END IF
            CALL cp_dbcsr_plus_fm_fm_t(matrix_wx1(ispin)%matrix, matrix_v=gs_mos(ispin)%mos_occ, &
                                       matrix_g=vcvec, &
                                       ncol=norb, alpha=alpha)
            CALL cp_dbcsr_plus_fm_fm_t(matrix_wx1(ispin)%matrix, matrix_v=vcvec, &
                                       matrix_g=gs_mos(ispin)%mos_occ, &
                                       ncol=norb, alpha=alpha)
            CALL cp_fm_release(vcvec)
            CALL cp_fm_release(cvcmat)
            !
            CALL dbcsr_release(pdens)
            CALL dbcsr_release(ptrans)
            DEALLOCATE (kind_of, atom_of_kind)
            !
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Exch[X]   ", fodeb
            END IF
         END IF
         !
         CALL cp_fm_release(cvec)
         CALL cp_fm_release(xvec)
         DEALLOCATE (tv)
      ENDDO

      DO ispin = 1, nspins
         CALL cp_fm_release(matrix=xtransformed(ispin)%matrix)
      END DO
      DEALLOCATE (xtransformed)
      DEALLOCATE (tcharge, gtcharge)
      DEALLOCATE (first_sgf, last_sgf)

      ! Lowdin forces
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_plo(1)%matrix, matrix_plo(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF
      CALL dbcsr_scale(matrix_plo(1)%matrix, -1.0_dp)
      NULLIFY (scrm)
      IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
      CALL build_overlap_matrix(ks_env, matrix_s=scrm, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", basis_type_b="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=matrix_plo(1)%matrix)
      CALL dbcsr_deallocate_matrix_set(scrm)
      CALL dbcsr_deallocate_matrix_set(matrix_plo)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%overlap(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Lowdin ", fodeb
      END IF

      IF (ASSOCIATED(ex_env%matrix_wx1)) CALL dbcsr_deallocate_matrix_set(ex_env%matrix_wx1)
      ex_env%matrix_wx1 => matrix_wx1

      CALL timestop(handle)

   END SUBROUTINE stda_force

! **************************************************************************************************

END MODULE qs_tddfpt2_fhxc_forces
