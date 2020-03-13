!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_forces
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE core_ppl,                        ONLY: build_core_ppl
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_type, dbcsr_dot, dbcsr_p_type, &
        dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type
   USE exstates_types,                  ONLY: excited_energy_type
   USE hfx_derivatives,                 ONLY: derivatives_four_center
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none,&
                                              kg_tnadd_atomic,&
                                              kg_tnadd_embed,&
                                              tddfpt_kernel_full,&
                                              tddfpt_kernel_none,&
                                              tddfpt_kernel_stda
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kg_correction,                   ONLY: kg_ekin_subset
   USE kg_environment_types,            ONLY: kg_environment_type
   USE kg_tnadd_mat,                    ONLY: build_tnadd_mat
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
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
                                              pw_p_type,&
                                              pw_release
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_density_matrices,             ONLY: calculate_wx_matrix,&
                                              calculate_xwx_matrix
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: allocate_qs_force,&
                                              deallocate_qs_force,&
                                              qs_force_type,&
                                              sum_qs_force,&
                                              total_qs_force,&
                                              zero_qs_force
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_tddfpt2_fhxc_forces,          ONLY: fhxc_force,&
                                              stda_force
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos,&
                                              tddfpt_kernel_env_type,&
                                              tddfpt_work_matrices
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE response_solver,                 ONLY: admm_projection_derivative,&
                                              ks_ref_potential
   USE task_list_types,                 ONLY: task_list_type
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_forces'

   PUBLIC :: tddfpt_forces, tddfpt_resvec1, tddfpt_resvec1_admm, tddfpt_resvec2, tddfpt_resvec3

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Calculate direct tddft forces
!> \param qs_env ...
!> \param ex_env ...
!> \param gs_mos ...
!> \param kernel_env ...
!> \param sub_env ...
!> \param work_matrices ...
!> \par History
!>    * 01.2020 screated [JGH]
! **************************************************************************************************
   SUBROUTINE tddfpt_forces(qs_env, ex_env, gs_mos, kernel_env, sub_env, work_matrices)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(tddfpt_kernel_env_type), INTENT(IN)           :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work_matrices

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_forces', routineP = moduleN//':'//routineN
      LOGICAL, PARAMETER                                 :: debug_forces = .TRUE.

      INTEGER                                            :: handle, iab, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: natom_of_kind
      REAL(KIND=dp)                                      :: ehartree, exc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vtau_rspace, vxc_rspace
      TYPE(pw_p_type), POINTER                           :: vh_rspace
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: ks_force, td_force

      CALL timeset(routineN, handle)

      ! prepare force array
      CALL get_qs_env(qs_env, force=ks_force, atomic_kind_set=atomic_kind_set)
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (natom_of_kind(nkind))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom_of_kind=natom_of_kind)
      NULLIFY (td_force)
      CALL allocate_qs_force(td_force, natom_of_kind)
      DEALLOCATE (natom_of_kind)
      CALL zero_qs_force(td_force)
      CALL set_qs_env(qs_env, force=td_force)
      !
      vh_rspace => ex_env%vh_rspace
      vxc_rspace => ex_env%vxc_rspace
      vtau_rspace => ex_env%vtau_rspace
      IF (ASSOCIATED(vh_rspace)) THEN
         CALL pw_release(vh_rspace%pw)
         DEALLOCATE (vh_rspace)
      END IF
      IF (ASSOCIATED(vxc_rspace)) THEN
         DO iab = 1, SIZE(vxc_rspace)
            CALL pw_release(vxc_rspace(iab)%pw)
         END DO
         DEALLOCATE (vxc_rspace)
      END IF
      IF (ASSOCIATED(vtau_rspace)) THEN
         DO iab = 1, SIZE(vtau_rspace)
            CALL pw_release(vtau_rspace(iab)%pw)
         END DO
         DEALLOCATE (vtau_rspace)
      END IF
      !
      NULLIFY (vh_rspace, vxc_rspace, vtau_rspace)
      CALL ks_ref_potential(qs_env, vh_rspace, vxc_rspace, vtau_rspace, ehartree, exc)
      ex_env%vh_rspace => vh_rspace
      ex_env%vxc_rspace => vxc_rspace
      ex_env%vtau_rspace => vtau_rspace
      !
      CALL tddfpt_force_direct(qs_env, ex_env, gs_mos, kernel_env, sub_env, work_matrices, debug_forces)
      !
      ! add TD and KS forces
      CALL get_qs_env(qs_env, force=td_force)
      CALL sum_qs_force(ks_force, td_force)
      CALL set_qs_env(qs_env, force=ks_force)
      CALL deallocate_qs_force(td_force)
      !
      CALL timestop(handle)

   END SUBROUTINE tddfpt_forces

! **************************************************************************************************
!> \brief Calculate direct tddft forces
!> \param qs_env ...
!> \param ex_env ...
!> \param gs_mos ...
!> \param kernel_env ...
!> \param sub_env ...
!> \param work_matrices ...
!> \param debug_forces ...
!> \par History
!>    * 01.2020 screated [JGH]
! **************************************************************************************************
   SUBROUTINE tddfpt_force_direct(qs_env, ex_env, gs_mos, kernel_env, sub_env, work_matrices, &
                                  debug_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(tddfpt_kernel_env_type), INTENT(IN)           :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work_matrices
      LOGICAL                                            :: debug_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_force_direct', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iounit, ispin, mspin, n_rep_hf, &
                                                            nao, nao_aux, natom, nder, nimages, &
                                                            norb, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: distribute_fock_matrix, do_hfx, &
                                                            hfx_treat_lsd_in_core, lsd, resp_only, &
                                                            s_mstruct_changed, use_virial
      REAL(KIND=dp)                                      :: eh1, ekin_mol, eps_ppnl, exc, &
                                                            total_rhoz, total_rhoz_aux, zehartree
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ftot1, ftot2, ftot3
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: evect
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_pe, matrix_pe_admm, &
                                                            matrix_wx1, matrix_wz, scrm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p, mhz, mpz
      TYPE(dbcsr_type), POINTER                          :: dbwork
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoz_tot_gspace, vhxc_rspace, &
                                                            zv_hartree_gspace, zv_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER :: rho_g, rho_g_aux, rho_r, rho_r_aux, rhoz_g, &
         rhoz_g_aux, rhoz_r, rhoz_r_aux, tau_r, tau_r_aux, v_rspace, v_tau_rspace, v_xc, &
         vtau_rspace, vxc_rspace
      TYPE(pw_p_type), POINTER                           :: vh_rspace
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(section_vals_type), POINTER                   :: hfx_section, xc_fun_section, xc_section
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho1_set, rho_set

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_unit_nr(logger)

      evect => ex_env%evect
      matrix_pe => ex_env%matrix_pe
      matrix_pe_admm => ex_env%matrix_pe_admm
      vh_rspace => ex_env%vh_rspace
      vxc_rspace => ex_env%vxc_rspace
      vtau_rspace => ex_env%vtau_rspace

      NULLIFY (ks_env, sab_orb, sac_ppl, sap_ppnl)
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, para_env=para_env, &
                      sab_orb=sab_orb, sac_ppl=sac_ppl, sap_ppnl=sap_ppnl, &
                      force=force)
      IF (debug_forces) THEN
         CALL get_qs_env(qs_env, natom=natom, atomic_kind_set=atomic_kind_set)
         ALLOCATE (ftot1(3, natom))
         CALL total_qs_force(ftot1, force, atomic_kind_set)
      END IF

      nspins = SIZE(matrix_pe, 1)
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_pe(1)%matrix, matrix_pe(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF
      ! Kinetic energy matrix
      NULLIFY (scrm)
      IF (debug_forces) fodeb(1:3) = force(1)%kinetic(1:3, 1)
      CALL build_kinetic_matrix(ks_env, matrix_t=scrm, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=matrix_pe(1)%matrix)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%kinetic(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pe*dT      ", fodeb
      END IF
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_pe(1)%matrix, matrix_pe(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      END IF

      ! pseudopotential
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, dft_control=dft_control)
      NULLIFY (cell_to_index, virial)
      use_virial = .FALSE.
      nimages = 1
      ALLOCATE (matrix_p(nspins, 1), matrix_h(nspins, 1))
      DO ispin = 1, nspins
         matrix_p(ispin, 1)%matrix => matrix_pe(ispin)%matrix
         matrix_h(ispin, 1)%matrix => scrm(1)%matrix
      END DO
      IF (ASSOCIATED(sac_ppl)) THEN
         nder = 1
         IF (debug_forces) fodeb(1:3) = force(1)%gth_ppl(1:3, 1)
         CALL build_core_ppl(matrix_h, matrix_p, force, &
                             virial, .TRUE., use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ppl, &
                             nimages, cell_to_index, "ORB")
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppl(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pe*dHppl   ", fodeb
         END IF
      END IF
      eps_ppnl = dft_control%qs_control%eps_ppnl
      IF (ASSOCIATED(sap_ppnl)) THEN
         nder = 1
         IF (debug_forces) fodeb(1:3) = force(1)%gth_ppnl(1:3, 1)
         CALL build_core_ppnl(matrix_h, matrix_p, force, &
                              virial, .TRUE., use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, sab_orb, sap_ppnl, eps_ppnl, &
                              nimages, cell_to_index, "ORB")
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppnl(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pe*dHppnl  ", fodeb
         END IF
      END IF
      ! Kim-Gordon subsystem DFT
      ! Atomic potential for nonadditive kinetic energy contribution
      IF (dft_control%qs_control%do_kg) THEN
         IF (qs_env%kg_env%tnadd_method == kg_tnadd_atomic) THEN
            CALL get_qs_env(qs_env=qs_env, kg_env=kg_env, dbcsr_dist=dbcsr_dist)
            IF (debug_forces) fodeb(1:3) = force(1)%kinetic(1:3, 1)
            CALL build_tnadd_mat(kg_env=kg_env, matrix_p=matrix_p, force=force, virial=virial, &
                                 calculate_forces=.TRUE., use_virial=use_virial, &
                                 qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set, &
                                 particle_set=particle_set, sab_orb=sab_orb, dbcsr_dist=dbcsr_dist)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%kinetic(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dTnadd  ", fodeb
            END IF
         END IF
      END IF
      DEALLOCATE (matrix_p, matrix_h)
      IF (debug_forces) THEN
         ALLOCATE (ftot2(3, natom))
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1) - ftot1(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Force Pe*dHcore", fodeb
      END IF

      ! Vhxc
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      CALL pw_pool_create_pw(auxbas_pw_pool, vhxc_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      DO ispin = 1, nspins
         CALL pw_transfer(vh_rspace%pw, vhxc_rspace%pw)
         CALL pw_axpy(vxc_rspace(ispin)%pw, vhxc_rspace%pw)
         CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                 hmat=scrm(1), pmat=matrix_pe(ispin), &
                                 qs_env=qs_env, calculate_forces=.TRUE.)
         IF (ASSOCIATED(vtau_rspace)) THEN
            CALL integrate_v_rspace(v_rspace=vtau_rspace(ispin), &
                                    hmat=scrm(1), pmat=matrix_pe(ispin), &
                                    qs_env=qs_env, calculate_forces=.TRUE., compute_tau=.TRUE.)
         END IF
      END DO
      CALL pw_pool_give_back_pw(auxbas_pw_pool, vhxc_rspace%pw)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pe*dVhxc   ", fodeb
      END IF
      ! KG Embedding
      ! calculate kinetic energy potential and integrate with transition density
      IF (dft_control%qs_control%do_kg) THEN
         IF (qs_env%kg_env%tnadd_method == kg_tnadd_embed) THEN
            CALL get_qs_env(qs_env, kg_env=kg_env)
            ekin_mol = 0.0_dp
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            CALL kg_ekin_subset(qs_env=qs_env, ks_matrix=scrm, ekin_mol=ekin_mol, &
                                calc_force=.TRUE., do_kernel=.FALSE., pmat_ext=matrix_pe)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pe*dVkg   ", fodeb
            END IF
         END IF
      END IF

      ALLOCATE (rhoz_r(nspins), rhoz_g(nspins))
      DO ispin = 1, nspins
         NULLIFY (rhoz_r(ispin)%pw, rhoz_g(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_r(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_g(ispin)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      ENDDO
      CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, zv_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, zv_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      CALL pw_zero(rhoz_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL calculate_rho_elec(ks_env=ks_env, matrix_p=matrix_pe(ispin)%matrix, &
                                 rho=rhoz_r(ispin), rho_gspace=rhoz_g(ispin), &
                                 total_rho=total_rhoz)
         CALL pw_axpy(rhoz_g(ispin)%pw, rhoz_tot_gspace%pw)
      END DO
      ! calculate associated hartree potential
      CALL pw_poisson_solve(poisson_env, rhoz_tot_gspace%pw, zehartree, &
                            zv_hartree_gspace%pw)
      CALL pw_transfer(zv_hartree_gspace%pw, zv_hartree_rspace%pw)
      CALL pw_scale(zv_hartree_rspace%pw, zv_hartree_rspace%pw%pw_grid%dvol)
      ! Getting nuclear force contribution from the core charge density
      IF (debug_forces) fodeb(1:3) = force(1)%rho_core(1:3, 1)
      CALL integrate_v_core_rspace(zv_hartree_rspace, qs_env)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_core(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Vh(rhoe)*dncore ", fodeb
      END IF
      !
      CALL get_qs_env(qs_env=qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g, tau_r=tau_r)
      ALLOCATE (v_xc(nspins))
      DO ispin = 1, nspins
         NULLIFY (v_xc(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(v_xc(ispin)%pw)
      END DO
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      END IF
      lsd = (nspins == 2)
      NULLIFY (deriv_set, rho_set, rho1_set)
      CALL xc_prep_2nd_deriv(deriv_set, rho_set, rho_r, auxbas_pw_pool, xc_section=xc_section)
      bo = rhoz_r(1)%pw%pw_grid%bounds_local
      CALL xc_rho_set_create(rho1_set, bo, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

      ! calculate the arguments needed by the functionals
      CALL xc_rho_set_update(rho1_set, rhoz_r, rhoz_g, tau_r, needs, &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                             auxbas_pw_pool)
      CALL xc_calc_2nd_deriv(v_xc, deriv_set, rho_set, &
                             rho1_set, auxbas_pw_pool, xc_section=xc_section)
      CALL xc_dset_release(deriv_set)
      CALL xc_rho_set_release(rho_set)
      CALL xc_rho_set_release(rho1_set)
      !
      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      DO ispin = 1, nspins
         CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(zv_hartree_rspace%pw, v_xc(ispin)%pw)
         CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                 hmat=scrm(1), &
                                 pmat=matrix_p(ispin, 1), &
                                 calculate_forces=.TRUE.)
      END DO
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dK*rhoe ", fodeb
      END IF
      ! KG Embedding
      ! calculate kinetic energy kernel, folded with response density for partial integration0
      IF (dft_control%qs_control%do_kg) THEN
         IF (qs_env%kg_env%tnadd_method == kg_tnadd_embed) THEN
            CALL get_qs_env(qs_env, kg_env=kg_env)
            ekin_mol = 0.0_dp
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            CALL kg_ekin_subset(qs_env=qs_env, ks_matrix=scrm, ekin_mol=ekin_mol, &
                                calc_force=.TRUE., do_kernel=.TRUE., pmat_ext=matrix_pe)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*d(Kkg)*rhoz ", fodeb
            END IF
         END IF
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_tot_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, zv_hartree_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, zv_hartree_gspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_r(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_g(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      END DO
      DEALLOCATE (rhoz_r, rhoz_g, v_xc)

      CALL dbcsr_deallocate_matrix_set(scrm)

      IF (debug_forces) THEN
         ALLOCATE (ftot3(3, natom))
         CALL total_qs_force(ftot3, force, atomic_kind_set)
         fodeb(1:3) = ftot3(1:3, 1) - ftot2(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Force Pin*V(rhoe)", fodeb
      END IF

      IF (dft_control%do_admm) THEN
         IF (dft_control%admm_control%aux_exch_func == do_admm_aux_exch_func_none) THEN
            ! nothing to do
         ELSE
            ! add ADMM xc_section_aux terms: Pz*Vxc + P0*K0[rhoz]
            CALL get_qs_env(qs_env, admm_env=admm_env, rho_aux_fit=rho_aux_fit, &
                            matrix_s_aux_fit=scrm, task_list_aux_fit=task_list_aux_fit)
            !
            NULLIFY (mpz, mhz)
            ALLOCATE (mpz(nspins, 1))
            CALL dbcsr_allocate_matrix_set(mhz, nspins, 1)
            DO ispin = 1, nspins
               ALLOCATE (mhz(ispin, 1)%matrix)
               CALL dbcsr_create(mhz(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_copy(mhz(ispin, 1)%matrix, scrm(1)%matrix)
               mpz(ispin, 1)%matrix => matrix_pe_admm(ispin)%matrix
            END DO
            !
            xc_section => admm_env%xc_section_aux
            NULLIFY (v_rspace, v_tau_rspace)
            CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_aux_fit, xc_section=xc_section, &
                               vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE.)
            CPASSERT(.NOT. ASSOCIATED(v_tau_rspace))
            DO ispin = 1, nspins
               CALL dbcsr_set(mhz(ispin, 1)%matrix, 0.0_dp)
               v_rspace(ispin)%pw%cr3d = v_rspace(ispin)%pw%pw_grid%dvol*v_rspace(ispin)%pw%cr3d
               CALL integrate_v_rspace(v_rspace=v_rspace(ispin), &
                                       hmat=mhz(ispin, 1), pmat=mpz(ispin, 1), &
                                       qs_env=qs_env, calculate_forces=.TRUE., &
                                       basis_type="AUX_FIT", task_list_external=task_list_aux_fit)
            END DO
            NULLIFY (rho_g_aux, rho_r_aux, tau_r_aux, rhoz_g_aux, rhoz_r_aux)
            !
            CALL qs_rho_get(rho_aux_fit, rho_r=rho_r_aux, rho_g=rho_g_aux, tau_r=tau_r_aux)
            ALLOCATE (v_xc(nspins))
            DO ispin = 1, nspins
               NULLIFY (v_xc(ispin)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_zero(v_xc(ispin)%pw)
            END DO
            lsd = (nspins == 2)
            ! rhoz_aux
            ALLOCATE (rhoz_r_aux(nspins), rhoz_g_aux(nspins))
            DO ispin = 1, nspins
               NULLIFY (rhoz_r_aux(ispin)%pw, rhoz_g_aux(ispin)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_r_aux(ispin)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_g_aux(ispin)%pw, &
                                      use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            ENDDO
            DO ispin = 1, nspins
               CALL calculate_rho_elec(ks_env=ks_env, matrix_p=mpz(ispin, 1)%matrix, &
                                       rho=rhoz_r_aux(ispin), rho_gspace=rhoz_g_aux(ispin), &
                                       total_rho=total_rhoz_aux, basis_type="AUX_FIT", &
                                       task_list_external=task_list_aux_fit)
            END DO
            !
            NULLIFY (deriv_set, rho_set, rho1_set)
            CALL xc_prep_2nd_deriv(deriv_set, rho_set, rho_r_aux, auxbas_pw_pool, xc_section=xc_section)
            bo = rhoz_r_aux(1)%pw%pw_grid%bounds_local
            CALL xc_rho_set_create(rho1_set, bo, &
                                   rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                                   drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                                   tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

            xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
            needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

            ! calculate the arguments needed by the functionals
            CALL xc_rho_set_update(rho1_set, rhoz_r_aux, rhoz_g_aux, tau_r_aux, needs, &
                                   section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                                   section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                                   auxbas_pw_pool)
            CALL xc_calc_2nd_deriv(v_xc, deriv_set, rho_set, &
                                   rho1_set, auxbas_pw_pool, xc_section=xc_section)
            CALL xc_dset_release(deriv_set)
            CALL xc_rho_set_release(rho_set)
            CALL xc_rho_set_release(rho1_set)
            !
            CALL qs_rho_get(rho_aux_fit, rho_ao_kp=matrix_p)
            DO ispin = 1, nspins
               CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
               CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                       hmat=mhz(ispin, 1), pmat=matrix_p(ispin, 1), &
                                       calculate_forces=.TRUE., basis_type="AUX_FIT", &
                                       task_list_external=task_list_aux_fit)
            END DO
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
            CALL dbcsr_create(dbwork, template=scrm(1)%matrix)
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(mhz(ispin, 1)%matrix, admm_env%A, &
                                            admm_env%work_aux_orb, nao)
               CALL cp_gemm('T', 'N', nao, nao, nao_aux, &
                            1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                            admm_env%work_orb_orb)
               CALL dbcsr_copy(dbwork, scrm(ispin)%matrix)
               CALL dbcsr_set(dbwork, 0.0_dp)
               CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbwork, keep_sparsity=.TRUE.)
               CALL dbcsr_add(scrm(ispin)%matrix, dbwork, 1.0_dp, 1.0_dp)
            END DO
            CALL dbcsr_release(dbwork)
            DEALLOCATE (dbwork)
            CALL dbcsr_deallocate_matrix_set(mhz)
            DEALLOCATE (mpz)
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
            NULLIFY (mpz, mhz)
            ALLOCATE (mpz(nspins, 1))
            CALL dbcsr_allocate_matrix_set(mhz, nspins, 1)
            DO ispin = 1, nspins
               ALLOCATE (mhz(ispin, 1)%matrix)
               CALL dbcsr_create(mhz(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_copy(mhz(ispin, 1)%matrix, scrm(1)%matrix)
               CALL dbcsr_set(mhz(ispin, 1)%matrix, 0.0_dp)
               mpz(ispin, 1)%matrix => matrix_pe_admm(ispin)%matrix
            END DO
         ELSE
            CALL get_qs_env(qs_env=qs_env, matrix_s=scrm)
            NULLIFY (mpz, mhz)
            ALLOCATE (mpz(nspins, 1))
            CALL dbcsr_allocate_matrix_set(mhz, nspins, 1)
            DO ispin = 1, nspins
               ALLOCATE (mhz(ispin, 1)%matrix)
               CALL dbcsr_create(mhz(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_copy(mhz(ispin, 1)%matrix, scrm(1)%matrix)
               CALL dbcsr_set(mhz(ispin, 1)%matrix, 0.0_dp)
               mpz(ispin, 1)%matrix => matrix_pe(ispin)%matrix
            END DO
         END IF
         DO ispin = 1, mspin
            eh1 = 0.0
            CALL integrate_four_center(qs_env, x_data, mhz, eh1, mpz, hfx_section, &
                                       para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                       ispin=ispin)
         END DO
         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CPASSERT(ASSOCIATED(admm_env%work_aux_orb))
            CPASSERT(ASSOCIATED(admm_env%work_orb_orb))
            nao = admm_env%nao_orb
            nao_aux = admm_env%nao_aux_fit
            ALLOCATE (dbwork)
            CALL dbcsr_create(dbwork, template=scrm(1)%matrix)
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(mhz(ispin, 1)%matrix, admm_env%A, &
                                            admm_env%work_aux_orb, nao)
               CALL cp_gemm('T', 'N', nao, nao, nao_aux, &
                            1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                            admm_env%work_orb_orb)
               CALL dbcsr_copy(dbwork, scrm(ispin)%matrix)
               CALL dbcsr_set(dbwork, 0.0_dp)
               CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbwork, keep_sparsity=.TRUE.)
               CALL dbcsr_add(scrm(ispin)%matrix, dbwork, 1.0_dp, 1.0_dp)
            END DO
            CALL dbcsr_release(dbwork)
            DEALLOCATE (dbwork)
            ! derivatives Tr (Pz [A(T)H dA/dR])
            CALL admm_projection_derivative(qs_env, admm_env, mhz, matrix_pe)
         END IF

         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env=qs_env, rho_aux_fit=rho_aux_fit)
            CALL qs_rho_get(rho_aux_fit, rho_ao_kp=matrix_p)
            matrix_pe => matrix_pe_admm
         ELSE
            CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
         END IF
         resp_only = .TRUE.
         CALL derivatives_four_center(qs_env, matrix_p, matrix_pe, hfx_section, para_env, &
                                      1, use_virial, resp_only=resp_only)
         CALL dbcsr_deallocate_matrix_set(mhz)
         DEALLOCATE (mpz)
      END IF

      CALL tddfpt_kernel_force(qs_env, ex_env, gs_mos, kernel_env, sub_env, work_matrices, debug_forces)

      ! Overlap matrix
      matrix_wx1 => ex_env%matrix_wx1
      CALL get_qs_env(qs_env=qs_env, matrix_s=scrm, matrix_ks=matrix_ks)
      NULLIFY (matrix_wz)
      CALL dbcsr_allocate_matrix_set(matrix_wz, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_wz(ispin)%matrix)
         CALL dbcsr_create(matrix=matrix_wz(ispin)%matrix, template=scrm(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_wz(ispin)%matrix, sab_orb)
         CALL dbcsr_set(matrix_wz(ispin)%matrix, 0.0_dp)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, ncol_global=norb)
         CALL cp_dbcsr_plus_fm_fm_t(matrix_wz(ispin)%matrix, matrix_v=evect(ispin)%matrix, ncol=norb)
         CALL dbcsr_scale(matrix_wz(ispin)%matrix, ex_env%evalue)
         CALL calculate_wx_matrix(gs_mos(ispin)%mos_occ, evect(ispin)%matrix, matrix_ks(ispin)%matrix, &
                                  matrix_wz(ispin)%matrix)
      END DO
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_wz(1)%matrix, matrix_wz(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF
      NULLIFY (scrm)
      IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
      CALL build_overlap_matrix(ks_env, matrix_s=scrm, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", basis_type_b="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=matrix_wz(1)%matrix)
      CALL dbcsr_deallocate_matrix_set(scrm)
      CALL dbcsr_deallocate_matrix_set(matrix_wz)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%overlap(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Wx*dS ", fodeb
      END IF

      ! Overlap matrix
      CALL get_qs_env(qs_env=qs_env, matrix_s=scrm, matrix_ks=matrix_ks)
      NULLIFY (matrix_wz)
      CALL dbcsr_allocate_matrix_set(matrix_wz, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_wz(ispin)%matrix)
         CALL dbcsr_create(matrix=matrix_wz(ispin)%matrix, template=scrm(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_wz(ispin)%matrix, sab_orb)
         CALL dbcsr_set(matrix_wz(ispin)%matrix, 0.0_dp)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, ncol_global=norb)
         CALL calculate_xwx_matrix(gs_mos(ispin)%mos_occ, evect(ispin)%matrix, scrm(1)%matrix, &
                                   matrix_ks(ispin)%matrix, matrix_wz(ispin)%matrix, ex_env%evalue)
      END DO
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_wz(1)%matrix, matrix_wz(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF
      NULLIFY (scrm)
      IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
      CALL build_overlap_matrix(ks_env, matrix_s=scrm, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", basis_type_b="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=matrix_wz(1)%matrix)
      CALL dbcsr_deallocate_matrix_set(scrm)
      CALL dbcsr_deallocate_matrix_set(matrix_wz)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%overlap(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: xWx*dS ", fodeb
      END IF

      ! Overlap matrix
      IF (ASSOCIATED(matrix_wx1)) THEN
         IF (nspins == 2) THEN
            CALL dbcsr_add(matrix_wx1(1)%matrix, matrix_wx1(2)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         END IF
         NULLIFY (scrm)
         IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
         CALL build_overlap_matrix(ks_env, matrix_s=scrm, &
                                   matrix_name="OVERLAP MATRIX", &
                                   basis_type_a="ORB", basis_type_b="ORB", &
                                   sab_nl=sab_orb, calculate_forces=.TRUE., &
                                   matrix_p=matrix_wx1(1)%matrix)
         CALL dbcsr_deallocate_matrix_set(scrm)
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%overlap(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: WK*dS ", fodeb
         END IF
      END IF

      IF (debug_forces) THEN
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1) - ftot1(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Excitation Force", fodeb
         DEALLOCATE (ftot1, ftot2, ftot3)
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_force_direct

! **************************************************************************************************
!> \brief ...
!> \param evect ...
!> \param mos_occ ...
!> \param matrix_s ...
!> \param matrix_pe ...
! **************************************************************************************************
   SUBROUTINE tddfpt_resvec1(evect, mos_occ, matrix_s, matrix_pe)

      TYPE(cp_fm_type), POINTER                          :: evect, mos_occ
      TYPE(dbcsr_type), POINTER                          :: matrix_s, matrix_pe

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_resvec1', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iounit, nao, norb
      REAL(KIND=dp)                                      :: tmp
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct, fmstruct2
      TYPE(cp_fm_type), POINTER                          :: cxmat, xxmat
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      ! X*X^T
      CALL cp_fm_get_info(mos_occ, nrow_global=nao, ncol_global=norb)
      CALL cp_dbcsr_plus_fm_fm_t(matrix_pe, matrix_v=evect, ncol=norb)
      ! X^T*S*X
      CALL cp_fm_get_info(evect, matrix_struct=fmstruct)
      NULLIFY (fmstruct2)
      CALL cp_fm_struct_create(fmstruct=fmstruct2, template_fmstruct=fmstruct, &
                               nrow_global=norb, ncol_global=norb)
      CALL cp_fm_create(xxmat, matrix_struct=fmstruct2)
      CALL cp_fm_struct_release(fmstruct2)
      CALL cp_fm_create(cxmat, matrix_struct=fmstruct)
      CALL cp_dbcsr_sm_fm_multiply(matrix_s, evect, cxmat, norb, alpha=1.0_dp, beta=0.0_dp)
      CALL cp_gemm('T', 'N', norb, norb, nao, 1.0_dp, cxmat, evect, 0.0_dp, xxmat)
      CALL cp_gemm('N', 'N', nao, norb, norb, 1.0_dp, mos_occ, xxmat, 0.0_dp, cxmat)
      CALL cp_fm_release(xxmat)
      ! C*C^T*XX
      CALL cp_dbcsr_plus_fm_fm_t(matrix_pe, matrix_v=mos_occ, matrix_g=cxmat, &
                                 ncol=norb, alpha=-1.0_dp)
      CALL cp_fm_release(cxmat)
      !
      ! Test for Tr(Pe*S)=0
      CALL dbcsr_dot(matrix_pe, matrix_s, tmp)
      IF (ABS(tmp) > 1.e-08_dp) THEN
         logger => cp_get_default_logger()
         iounit = cp_logger_get_default_unit_nr(logger)
         CPWARN("Electron count of excitation density matrix is non-zero.")
         IF (iounit > 0) THEN
            WRITE (iounit, "(T2,A,T61,G20.10)") "Measured electron count is ", tmp
            WRITE (iounit, "(T2,A,/)") REPEAT("*", 79)
         END IF
      END IF
      !
      CALL timestop(handle)

   END SUBROUTINE tddfpt_resvec1

! **************************************************************************************************
!> \brief ...
!> \param evect ...
!> \param mos_occ ...
!> \param admm_env ...
!> \param matrix_pe_admm ...
! **************************************************************************************************
   SUBROUTINE tddfpt_resvec1_admm(evect, mos_occ, admm_env, matrix_pe_admm)

      TYPE(cp_fm_type), POINTER                          :: evect, mos_occ
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_type), POINTER                          :: matrix_pe_admm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_resvec1_admm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nao, nao_aux, norb
      TYPE(cp_fm_type), POINTER                          :: work_aux_aux, work_aux_nmo, &
                                                            work_aux_nmo2, work_nmo_nmo1, &
                                                            work_orb_nmo

      CALL timeset(routineN, handle)
      !
      CALL cp_fm_get_info(mos_occ, nrow_global=nao, ncol_global=norb)
      nao_aux = admm_env%nao_aux_fit
      !
      work_orb_nmo => admm_env%work_orb_nmo(1)%matrix
      work_nmo_nmo1 => admm_env%work_nmo_nmo1(1)%matrix
      work_aux_nmo => admm_env%work_aux_nmo(1)%matrix
      work_aux_nmo2 => admm_env%work_aux_nmo2(1)%matrix
      work_aux_aux => admm_env%work_aux_aux
      ! (A*X)*(A*X)^T
      CALL cp_gemm('N', 'N', nao_aux, norb, nao, 1.0_dp, admm_env%A, evect, 0.0_dp, work_aux_nmo)
      CALL cp_gemm('N', 'T', nao_aux, nao_aux, norb, 1.0_dp, work_aux_nmo, work_aux_nmo, 0.0_dp, work_aux_aux)
      ! X^T*B*X
      CALL cp_gemm('N', 'N', nao, norb, nao, 1.0_dp, admm_env%B, evect, 0.0_dp, work_orb_nmo)
      CALL cp_gemm('T', 'N', norb, norb, nao, 1.0_dp, work_orb_nmo, evect, 0.0_dp, work_nmo_nmo1)
      ! A*C
      CALL cp_gemm('N', 'N', nao_aux, norb, nao, 1.0_dp, admm_env%A, mos_occ, 0.0_dp, work_aux_nmo)
      ! A*C*XBX
      CALL cp_gemm('N', 'N', nao_aux, norb, norb, 1.0_dp, work_aux_nmo, work_nmo_nmo1, 0.0_dp, work_aux_nmo2)
      !
      CALL cp_gemm('N', 'T', nao_aux, nao_aux, norb, -1.0_dp, work_aux_nmo2, work_aux_nmo, 1.0_dp, work_aux_aux)
      !
      CALL copy_fm_to_dbcsr(work_aux_aux, matrix_pe_admm, keep_sparsity=.TRUE.)
      !
      CALL timestop(handle)

   END SUBROUTINE tddfpt_resvec1_admm

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_pe ...
!> \param matrix_pe_admm ...
!> \param gs_mos ...
!> \param matrix_hz ...
!> \param cpmos ...
! **************************************************************************************************
   SUBROUTINE tddfpt_resvec2(qs_env, matrix_pe, matrix_pe_admm, gs_mos, matrix_hz, cpmos)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pe, matrix_pe_admm
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: cpmos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_resvec2', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iounit, ispin, mspin, n_rep_hf, &
                                                            nao, nao_aux, norb, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: distribute_fock_matrix, do_hfx, &
                                                            hfx_treat_lsd_in_core, lsd, &
                                                            s_mstruct_changed
      REAL(KIND=dp)                                      :: eh1, focc, thartree, total_rho
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), POINTER                          :: mos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: msaux
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mhz, mpe
      TYPE(dbcsr_type), POINTER                          :: dbwork
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r, tau_r, trho_g, trho_r, v_xc
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: hfx_section, input, xc_fun_section, &
                                                            xc_section
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho1_set, rho_set

      CALL timeset(routineN, handle)

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, ks_env=ks_env, dft_control=dft_control)
      CPASSERT(ASSOCIATED(pw_env))
      nspins = dft_control%nspins

      NULLIFY (auxbas_pw_pool, poisson_env)
      ! gets the tmp grids
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      NULLIFY (v_hartree_gspace%pw, rho_tot_gspace%pw, v_hartree_rspace%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      ALLOCATE (trho_r(nspins), trho_g(nspins))
      DO ispin = 1, nspins
         NULLIFY (trho_r(ispin)%pw, trho_g(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, trho_r(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, trho_g(ispin)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      ENDDO
      CALL pw_zero(rho_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL calculate_rho_elec(ks_env=ks_env, matrix_p=matrix_pe(ispin)%matrix, &
                                 rho=trho_r(ispin), &
                                 rho_gspace=trho_g(ispin), &
                                 total_rho=total_rho)
         CALL pw_axpy(trho_g(ispin)%pw, rho_tot_gspace%pw)
         IF (ABS(total_rho) > 1.e-08_dp) THEN
            logger => cp_get_default_logger()
            iounit = cp_logger_get_default_unit_nr(logger)
            CPWARN("Real space electron count of excitation density is non-zero.")
            IF (iounit > 0) THEN
               WRITE (iounit, "(T2,A,T61,G20.10)") "Measured electron count is ", total_rho
               WRITE (iounit, "(T2,A,/)") REPEAT("*", 79)
            END IF
         END IF
      END DO
      ! calculate associated hartree potential
      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, thartree, &
                            v_hartree_gspace%pw)
      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

      ! Fxc*drho term
      CALL get_qs_env(qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g, tau_r=tau_r)
      ALLOCATE (v_xc(nspins))
      DO ispin = 1, nspins
         NULLIFY (v_xc(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(v_xc(ispin)%pw)
      END DO
      CALL get_qs_env(qs_env, input=input)
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      lsd = (nspins == 2)
      NULLIFY (deriv_set, rho_set, rho1_set)
      CALL xc_prep_2nd_deriv(deriv_set, rho_set, rho_r, auxbas_pw_pool, xc_section=xc_section)
      bo = trho_r(1)%pw%pw_grid%bounds_local
      CALL xc_rho_set_create(rho1_set, bo, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

      ! calculate the arguments needed by the functionals
      CALL xc_rho_set_update(rho1_set, trho_r, trho_g, tau_r, needs, &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                             auxbas_pw_pool)
      CALL xc_calc_2nd_deriv(v_xc, deriv_set, rho_set, &
                             rho1_set, auxbas_pw_pool, xc_section=xc_section)
      CALL xc_dset_release(deriv_set)
      CALL xc_rho_set_release(rho_set)
      CALL xc_rho_set_release(rho1_set)

      ! vtot = v_xc(ispin) + v_hartree
      DO ispin = 1, nspins
         CALL dbcsr_set(matrix_hz(ispin)%matrix, 0.0_dp)
         CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(v_hartree_rspace%pw, v_xc(ispin)%pw)
         CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                 hmat=matrix_hz(ispin), &
                                 calculate_forces=.FALSE.)
      END DO

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, trho_r(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, trho_g(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      END DO
      DEALLOCATE (trho_r, trho_g, v_xc)

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
         CALL get_qs_env(qs_env=qs_env, rho=rho, x_data=x_data, para_env=para_env, &
                         s_mstruct_changed=s_mstruct_changed)
         distribute_fock_matrix = .TRUE.
         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env, matrix_s_aux_fit=msaux)
            NULLIFY (mpe, mhz)
            ALLOCATE (mpe(nspins, 1))
            CALL dbcsr_allocate_matrix_set(mhz, nspins, 1)
            DO ispin = 1, nspins
               ALLOCATE (mhz(ispin, 1)%matrix)
               CALL dbcsr_create(mhz(ispin, 1)%matrix, template=msaux(1)%matrix)
               CALL dbcsr_copy(mhz(ispin, 1)%matrix, msaux(1)%matrix)
               CALL dbcsr_set(mhz(ispin, 1)%matrix, 0.0_dp)
               mpe(ispin, 1)%matrix => matrix_pe_admm(ispin)%matrix
            END DO
         ELSE
            ALLOCATE (mpe(nspins, 1), mhz(nspins, 1))
            DO ispin = 1, nspins
               mhz(ispin, 1)%matrix => matrix_hz(ispin)%matrix
               mpe(ispin, 1)%matrix => matrix_pe(ispin)%matrix
            END DO
         END IF
         DO ispin = 1, mspin
            eh1 = 0.0
            CALL integrate_four_center(qs_env, x_data, mhz, eh1, mpe, hfx_section, &
                                       para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                       ispin=ispin)
         END DO
         IF (dft_control%do_admm) THEN
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
         END IF
         IF (dft_control%do_admm) THEN
            CALL dbcsr_deallocate_matrix_set(mhz)
            DEALLOCATE (mpe)
         ELSE
            DEALLOCATE (mpe, mhz)
         END IF
      END IF

      focc = 4.0_dp
      IF (nspins == 2) focc = 2.0_dp
      DO ispin = 1, nspins
         mos => gs_mos(ispin)%mos_occ
         CALL cp_fm_get_info(mos, ncol_global=norb)
         CALL cp_dbcsr_sm_fm_multiply(matrix_hz(ispin)%matrix, mos, cpmos(ispin)%matrix, &
                                      norb, alpha=focc, beta=0.0_dp)
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt_resvec2

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cpmos ...
!> \param work ...
! **************************************************************************************************
   SUBROUTINE tddfpt_resvec3(qs_env, cpmos, work)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: cpmos
      TYPE(tddfpt_work_matrices)                         :: work

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_resvec3', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nao, norb, nspins
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), POINTER                          :: cvec, omos, rvecs, scmat, umat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, mos=mos, dft_control=dft_control)
      nspins = dft_control%nspins

      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=omos)
         scmat => work%S_C0(ispin)%matrix
         rvecs => cpmos(ispin)%matrix
         CALL cp_fm_get_info(rvecs, nrow_global=nao, ncol_global=norb)
         CALL cp_fm_create(cvec, rvecs%matrix_struct, "cvec")
         CALL cp_fm_struct_create(fmstruct, context=rvecs%matrix_struct%context, nrow_global=norb, &
                                  ncol_global=norb, para_env=rvecs%matrix_struct%para_env)
         CALL cp_fm_create(umat, fmstruct, "umat")
         CALL cp_fm_struct_release(fmstruct)
         !
         CALL cp_gemm("T", "N", norb, norb, nao, 1.0_dp, omos, scmat, 0.0_dp, umat)
         CALL cp_fm_copy_general(rvecs, cvec, rvecs%matrix_struct%para_env)
         CALL cp_gemm("N", "T", nao, norb, norb, 1.0_dp, cvec, umat, 0.0_dp, rvecs)
         CALL cp_fm_release(cvec)
         CALL cp_fm_release(umat)
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt_resvec3

! **************************************************************************************************
!> \brief Calculate direct tddft forces
!> \param qs_env ...
!> \param ex_env ...
!> \param gs_mos ...
!> \param kernel_env ...
!> \param sub_env ...
!> \param work_matrices ...
!> \param debug_forces ...
!> \par History
!>    * 01.2020 screated [JGH]
! **************************************************************************************************
   SUBROUTINE tddfpt_kernel_force(qs_env, ex_env, gs_mos, kernel_env, sub_env, work_matrices, debug_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(tddfpt_kernel_env_type), INTENT(IN)           :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work_matrices
      LOGICAL, INTENT(IN)                                :: debug_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_kernel_force', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      MARK_USED(work_matrices)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      tddfpt_control => dft_control%tddfpt2_control

      IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
         ! full Kernel
         CALL fhxc_force(qs_env, ex_env, gs_mos, kernel_env%full_kernel, debug_forces)
      ELSE IF (tddfpt_control%kernel == tddfpt_kernel_stda) THEN
         ! sTDA Kernel
         CALL stda_force(qs_env, ex_env, gs_mos, kernel_env%stda_kernel, sub_env, work_matrices, debug_forces)
      ELSE IF (tddfpt_control%kernel == tddfpt_kernel_none) THEN
         ! nothing to be done here
         ex_env%matrix_wx1 => NULL()
      ELSE
         CPABORT('Unknown kernel type')
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_kernel_force

END MODULE qs_tddfpt2_forces
