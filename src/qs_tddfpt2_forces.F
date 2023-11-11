!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_forces
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
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
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_complete_redistribute, dbcsr_copy, dbcsr_create, dbcsr_dot, dbcsr_p_type, &
        dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric
   USE exstates_types,                  ONLY: excited_energy_type,&
                                              exstate_potential_release
   USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals,&
                                              init_coulomb_local
   USE hartree_local_types,             ONLY: hartree_local_create,&
                                              hartree_local_release,&
                                              hartree_local_type
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_ri,                          ONLY: hfx_ri_update_ks
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none,&
                                              oe_shift,&
                                              tddfpt_kernel_full,&
                                              tddfpt_kernel_none,&
                                              tddfpt_kernel_stda
   USE input_section_types,             ONLY: section_get_lval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE mulliken,                        ONLY: ao_charges
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
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
   USE qs_fxc,                          ONLY: qs_fxc_analytic,&
                                              qs_fxc_fdiff
   USE qs_gapw_densities,               ONLY: prepare_gapw_den
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kernel_types,                 ONLY: kernel_env_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_atom,                      ONLY: update_ks_atom
   USE qs_ks_reference,                 ONLY: ks_ref_potential,&
                                              ks_ref_potential_atom
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_local_rho_types,              ONLY: local_rho_set_create,&
                                              local_rho_set_release,&
                                              local_rho_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_oce_types,                    ONLY: oce_matrix_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace,&
                                              rho0_s_grid_create
   USE qs_rho0_methods,                 ONLY: init_rho0
   USE qs_rho0_types,                   ONLY: get_rho0_mpole
   USE qs_rho_atom_methods,             ONLY: allocate_rho_atom_internals,&
                                              calculate_rho_atom_coeff
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_tddfpt2_fhxc_forces,          ONLY: fhxc_force,&
                                              stda_force
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos,&
                                              tddfpt_work_matrices
   USE qs_vxc_atom,                     ONLY: calculate_xc_2nd_deriv_atom
   USE task_list_types,                 ONLY: task_list_type
   USE xtb_ehess,                       ONLY: xtb_coulomb_hessian
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_forces'

   PUBLIC :: tddfpt_forces_main

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Perform TDDFPT gradient calculation.
!> \param qs_env  Quickstep environment
!> \param gs_mos ...
!> \param ex_env ...
!> \param kernel_env ...
!> \param sub_env ...
!> \param work_matrices ...
!> \par History
!>    * 10.2022 created JHU
! **************************************************************************************************
   SUBROUTINE tddfpt_forces_main(qs_env, gs_mos, ex_env, kernel_env, sub_env, work_matrices)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(kernel_env_type)                              :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work_matrices

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_forces_main'

      INTEGER                                            :: handle, ispin, nspins
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pe_asymm, matrix_pe_symm, &
                                                            matrix_s, matrix_s_aux_fit
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins
      tddfpt_control => dft_control%tddfpt2_control
      ! rhs of linres equation
      IF (ASSOCIATED(ex_env%cpmos)) THEN
         DO ispin = 1, SIZE(ex_env%cpmos)
            CALL cp_fm_release(ex_env%cpmos(ispin))
         END DO
         DEALLOCATE (ex_env%cpmos)
      END IF
      ALLOCATE (ex_env%cpmos(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_get_info(matrix=ex_env%evect(ispin), matrix_struct=matrix_struct)
         CALL cp_fm_create(ex_env%cpmos(ispin), matrix_struct)
         CALL cp_fm_set_all(ex_env%cpmos(ispin), 0.0_dp)
      END DO
      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
      NULLIFY (matrix_pe_asymm, matrix_pe_symm)
      CALL dbcsr_allocate_matrix_set(ex_env%matrix_pe, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_pe_symm, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_pe_asymm, nspins)
      DO ispin = 1, nspins
         ALLOCATE (ex_env%matrix_pe(ispin)%matrix)
         CALL dbcsr_create(ex_env%matrix_pe(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(ex_env%matrix_pe(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(ex_env%matrix_pe(ispin)%matrix, 0.0_dp)

         ALLOCATE (matrix_pe_symm(ispin)%matrix)
         CALL dbcsr_create(matrix_pe_symm(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(matrix_pe_symm(ispin)%matrix, ex_env%matrix_pe(ispin)%matrix)

         ALLOCATE (matrix_pe_asymm(ispin)%matrix)
         CALL dbcsr_create(matrix_pe_asymm(ispin)%matrix, template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL dbcsr_complete_redistribute(ex_env%matrix_pe(ispin)%matrix, matrix_pe_asymm(ispin)%matrix)

         CALL tddfpt_resvec1(ex_env%evect(ispin), gs_mos(ispin)%mos_occ, &
                             matrix_s(1)%matrix, ex_env%matrix_pe(ispin)%matrix)
      END DO
      !
      ! ground state ADMM!
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         CALL get_admm_env(admm_env, matrix_s_aux_fit=matrix_s_aux_fit)
         CALL dbcsr_allocate_matrix_set(ex_env%matrix_pe_admm, nspins)
         DO ispin = 1, nspins
            ALLOCATE (ex_env%matrix_pe_admm(ispin)%matrix)
            CALL dbcsr_create(ex_env%matrix_pe_admm(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix)
            CALL dbcsr_copy(ex_env%matrix_pe_admm(ispin)%matrix, matrix_s_aux_fit(1)%matrix)
            CALL dbcsr_set(ex_env%matrix_pe_admm(ispin)%matrix, 0.0_dp)
            CALL tddfpt_resvec1_admm(ex_env%matrix_pe(ispin)%matrix, &
                                     admm_env, ex_env%matrix_pe_admm(ispin)%matrix)
         END DO
      END IF
      !
      CALL dbcsr_allocate_matrix_set(ex_env%matrix_hz, nspins)
      DO ispin = 1, nspins
         ALLOCATE (ex_env%matrix_hz(ispin)%matrix)
         CALL dbcsr_create(ex_env%matrix_hz(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(ex_env%matrix_hz(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(ex_env%matrix_hz(ispin)%matrix, 0.0_dp)
      END DO
      IF (dft_control%qs_control%xtb) THEN
         CALL tddfpt_resvec2_xtb(qs_env, ex_env%matrix_pe, gs_mos, ex_env%matrix_hz, ex_env%cpmos)
      ELSE
         CALL tddfpt_resvec2(qs_env, ex_env%matrix_pe, ex_env%matrix_pe_admm, &
                             gs_mos, ex_env%matrix_hz, ex_env%cpmos)
      END IF
      !
      CALL dbcsr_allocate_matrix_set(ex_env%matrix_px1, nspins)
      CALL dbcsr_allocate_matrix_set(ex_env%matrix_px1_asymm, nspins)
      DO ispin = 1, nspins
         ALLOCATE (ex_env%matrix_px1(ispin)%matrix)
         CALL dbcsr_create(ex_env%matrix_px1(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(ex_env%matrix_px1(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(ex_env%matrix_px1(ispin)%matrix, 0.0_dp)

         ALLOCATE (ex_env%matrix_px1_asymm(ispin)%matrix)
         CALL dbcsr_create(ex_env%matrix_px1_asymm(ispin)%matrix, template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL dbcsr_complete_redistribute(ex_env%matrix_px1(ispin)%matrix, ex_env%matrix_px1_asymm(ispin)%matrix)
      END DO
      ! Kernel ADMM
      IF (tddfpt_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         CALL get_admm_env(admm_env, matrix_s_aux_fit=matrix_s_aux_fit)
         CALL dbcsr_allocate_matrix_set(ex_env%matrix_px1_admm, nspins)
         CALL dbcsr_allocate_matrix_set(ex_env%matrix_px1_admm_asymm, nspins)
         DO ispin = 1, nspins
            ALLOCATE (ex_env%matrix_px1_admm(ispin)%matrix)
            CALL dbcsr_create(ex_env%matrix_px1_admm(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix)
            CALL dbcsr_copy(ex_env%matrix_px1_admm(ispin)%matrix, matrix_s_aux_fit(1)%matrix)
            CALL dbcsr_set(ex_env%matrix_px1_admm(ispin)%matrix, 0.0_dp)

            ALLOCATE (ex_env%matrix_px1_admm_asymm(ispin)%matrix)
            CALL dbcsr_create(ex_env%matrix_px1_admm_asymm(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL dbcsr_complete_redistribute(ex_env%matrix_px1_admm(ispin)%matrix, &
                                             ex_env%matrix_px1_admm_asymm(ispin)%matrix)
         END DO
      END IF
      ! TDA forces
      CALL tddfpt_forces(qs_env, ex_env, gs_mos, kernel_env, sub_env, work_matrices)
      ! Rotate res vector cpmos into original frame of occupied orbitals
      CALL tddfpt_resvec3(qs_env, ex_env%cpmos, work_matrices)

      CALL dbcsr_deallocate_matrix_set(matrix_pe_symm)
      CALL dbcsr_deallocate_matrix_set(matrix_pe_asymm)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_forces_main

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
      TYPE(kernel_env_type), INTENT(IN)                  :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work_matrices

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt_forces'

      INTEGER                                            :: handle
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: natom_of_kind
      LOGICAL                                            :: debug_forces
      REAL(KIND=dp)                                      :: ehartree, exc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: ks_force, td_force

      CALL timeset(routineN, handle)

      ! for extended debug output
      debug_forces = ex_env%debug_forces
      ! prepare force array
      CALL get_qs_env(qs_env, dft_control=dft_control, force=ks_force, &
                      atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom_of_kind=natom_of_kind)
      NULLIFY (td_force)
      CALL allocate_qs_force(td_force, natom_of_kind)
      DEALLOCATE (natom_of_kind)
      CALL zero_qs_force(td_force)
      CALL set_qs_env(qs_env, force=td_force)
      !
      IF (dft_control%qs_control%xtb) THEN
         CALL tddfpt_force_direct(qs_env, ex_env, gs_mos, kernel_env, sub_env, &
                                  work_matrices, debug_forces)
      ELSE
         !
         CALL exstate_potential_release(ex_env)
         CALL ks_ref_potential(qs_env, ex_env%vh_rspace, ex_env%vxc_rspace, &
                               ex_env%vtau_rspace, ex_env%vadmm_rspace, ehartree, exc)
         CALL ks_ref_potential_atom(qs_env, ex_env%local_rho_set, ex_env%local_rho_set_admm, &
                                    ex_env%vh_rspace)
         CALL tddfpt_force_direct(qs_env, ex_env, gs_mos, kernel_env, sub_env, &
                                  work_matrices, debug_forces)
      END IF
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
      TYPE(kernel_env_type), INTENT(IN)                  :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work_matrices
      LOGICAL                                            :: debug_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_force_direct'

      INTEGER                                            :: handle, iounit, ispin, natom, norb, &
                                                            nspins
      REAL(KIND=dp)                                      :: evalue
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ftot1, ftot2
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: evect
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s, matrix_wx1, &
                                                            matrix_wz, scrm
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      END IF

      evect => ex_env%evect

      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, para_env=para_env, &
                      sab_orb=sab_orb, dft_control=dft_control, force=force)
      NULLIFY (tddfpt_control)
      tddfpt_control => dft_control%tddfpt2_control
      nspins = dft_control%nspins

      IF (debug_forces) THEN
         CALL get_qs_env(qs_env, natom=natom, atomic_kind_set=atomic_kind_set)
         ALLOCATE (ftot1(3, natom))
         CALL total_qs_force(ftot1, force, atomic_kind_set)
      END IF

      CALL tddfpt_kernel_force(qs_env, ex_env, gs_mos, kernel_env, sub_env, work_matrices, debug_forces)

      ! Overlap matrix
      matrix_wx1 => ex_env%matrix_wx1
      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, matrix_ks=matrix_ks)
      NULLIFY (matrix_wz)
      CALL dbcsr_allocate_matrix_set(matrix_wz, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_wz(ispin)%matrix)
         CALL dbcsr_create(matrix=matrix_wz(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_wz(ispin)%matrix, sab_orb)
         CALL dbcsr_set(matrix_wz(ispin)%matrix, 0.0_dp)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, ncol_global=norb)
         CALL cp_dbcsr_plus_fm_fm_t(matrix_wz(ispin)%matrix, matrix_v=evect(ispin), ncol=norb)
         evalue = ex_env%evalue
         IF (tddfpt_control%oe_corr == oe_shift) THEN
            evalue = ex_env%evalue - tddfpt_control%ev_shift
         END IF
         CALL dbcsr_scale(matrix_wz(ispin)%matrix, evalue)
         CALL calculate_wx_matrix(gs_mos(ispin)%mos_occ, evect(ispin), matrix_ks(ispin)%matrix, &
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
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Wx*dS ", fodeb
      END IF

      ! Overlap matrix
      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, matrix_ks=matrix_ks)
      NULLIFY (matrix_wz)
      CALL dbcsr_allocate_matrix_set(matrix_wz, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_wz(ispin)%matrix)
         CALL dbcsr_create(matrix=matrix_wz(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_wz(ispin)%matrix, sab_orb)
         CALL dbcsr_set(matrix_wz(ispin)%matrix, 0.0_dp)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, ncol_global=norb)
         evalue = ex_env%evalue
         IF (tddfpt_control%oe_corr == oe_shift) THEN
            evalue = ex_env%evalue - tddfpt_control%ev_shift
         END IF
         CALL calculate_xwx_matrix(gs_mos(ispin)%mos_occ, evect(ispin), matrix_s(1)%matrix, &
                                   matrix_ks(ispin)%matrix, matrix_wz(ispin)%matrix, evalue)
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
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: xWx*dS ", fodeb
      END IF

      ! Overlap matrix
      IF (ASSOCIATED(matrix_wx1)) THEN
         IF (nspins == 2) THEN
            CALL dbcsr_add(matrix_wx1(1)%matrix, matrix_wx1(2)%matrix, &
                           alpha_scalar=0.5_dp, beta_scalar=0.5_dp)
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
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: WK*dS ", fodeb
         END IF
      END IF

      IF (debug_forces) THEN
         ALLOCATE (ftot2(3, natom))
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1) - ftot1(1:3, 1)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Excitation Force", fodeb
         DEALLOCATE (ftot1, ftot2)
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

      TYPE(cp_fm_type), INTENT(IN)                       :: evect, mos_occ
      TYPE(dbcsr_type), POINTER                          :: matrix_s, matrix_pe

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt_resvec1'

      INTEGER                                            :: handle, iounit, nao, norb
      REAL(KIND=dp)                                      :: tmp
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct, fmstruct2
      TYPE(cp_fm_type)                                   :: cxmat, xxmat
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
      CALL parallel_gemm('T', 'N', norb, norb, nao, 1.0_dp, cxmat, evect, 0.0_dp, xxmat)
      CALL parallel_gemm('N', 'N', nao, norb, norb, 1.0_dp, mos_occ, xxmat, 0.0_dp, cxmat)
      CALL cp_fm_release(xxmat)
      ! C*C^T*XX
      CALL cp_dbcsr_plus_fm_fm_t(matrix_pe, matrix_v=mos_occ, matrix_g=cxmat, &
                                 ncol=norb, alpha=-1.0_dp, symmetry_mode=1)
      CALL cp_fm_release(cxmat)
      !
      ! Test for Tr(Pe*S)=0
      CALL dbcsr_dot(matrix_pe, matrix_s, tmp)
      IF (ABS(tmp) > 1.e-08_dp) THEN
         logger => cp_get_default_logger()
         IF (logger%para_env%is_source()) THEN
            iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
         ELSE
            iounit = -1
         END IF
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
!> \brief PA = A * P * A(T)
!> \param matrix_pe ...
!> \param admm_env ...
!> \param matrix_pe_admm ...
! **************************************************************************************************
   SUBROUTINE tddfpt_resvec1_admm(matrix_pe, admm_env, matrix_pe_admm)

      TYPE(dbcsr_type), POINTER                          :: matrix_pe
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_type), POINTER                          :: matrix_pe_admm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_resvec1_admm'

      INTEGER                                            :: handle, nao, nao_aux

      CALL timeset(routineN, handle)
      !
      nao_aux = admm_env%nao_aux_fit
      nao = admm_env%nao_orb
      !
      CALL copy_dbcsr_to_fm(matrix_pe, admm_env%work_orb_orb)
      CALL parallel_gemm('N', 'N', nao_aux, nao, nao, &
                         1.0_dp, admm_env%A, admm_env%work_orb_orb, 0.0_dp, &
                         admm_env%work_aux_orb)
      CALL parallel_gemm('N', 'T', nao_aux, nao_aux, nao, &
                         1.0_dp, admm_env%work_aux_orb, admm_env%A, 0.0_dp, &
                         admm_env%work_aux_aux)
      CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, matrix_pe_admm, keep_sparsity=.TRUE.)
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
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: cpmos

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt_resvec2'

      CHARACTER(LEN=default_string_length)               :: basis_type
      INTEGER                                            :: handle, iounit, ispin, mspin, n_rep_hf, &
                                                            nao, nao_aux, natom, norb, nspins
      LOGICAL                                            :: deriv2_analytic, distribute_fock_matrix, &
                                                            do_hfx, gapw, gapw_xc, &
                                                            hfx_treat_lsd_in_core, &
                                                            s_mstruct_changed
      REAL(KIND=dp)                                      :: eh1, focc, rhotot, thartree
      REAL(KIND=dp), DIMENSION(2)                        :: total_rho
      REAL(KIND=dp), DIMENSION(:), POINTER               :: Qlm_tot
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_type), POINTER                          :: mos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: msaux
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mhz, mpe
      TYPE(dbcsr_type), POINTER                          :: dbwork
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hartree_local_type), POINTER                  :: hartree_local
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(local_rho_type), POINTER                      :: local_rho_set, local_rho_set_admm
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab, sab_aux_fit
      TYPE(oce_matrix_type), POINTER                     :: oce
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_type), DIMENSION(:), POINTER               :: rho_g, rho_g_aux, rho_r, rho_r_aux, &
                                                            rhoz_g_aux, rhoz_r_aux, tau_r, trho_g, &
                                                            trho_r, trho_xc_g, trho_xc_r, v_xc, &
                                                            v_xc_tau
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit, rho_xc, rhoz_aux, trho
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho1_atom_set, rho_atom_set
      TYPE(section_vals_type), POINTER                   :: hfx_section, input, xc_section
      TYPE(task_list_type), POINTER                      :: task_list

      CALL timeset(routineN, handle)

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, ks_env=ks_env, &
                      dft_control=dft_control, para_env=para_env)
      CPASSERT(ASSOCIATED(pw_env))
      nspins = dft_control%nspins
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc

      NULLIFY (auxbas_pw_pool, poisson_env)
      ! gets the tmp grids
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      CALL auxbas_pw_pool%create_pw(v_hartree_gspace, &
                                    use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL auxbas_pw_pool%create_pw(rho_tot_gspace, &
                                    use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL auxbas_pw_pool%create_pw(v_hartree_rspace, &
                                    use_data=REALDATA3D, in_space=REALSPACE)

      ALLOCATE (trho_r(nspins), trho_g(nspins))
      DO ispin = 1, nspins
         CALL auxbas_pw_pool%create_pw(trho_r(ispin), &
                                       use_data=REALDATA3D, in_space=REALSPACE)
         CALL auxbas_pw_pool%create_pw(trho_g(ispin), &
                                       use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      END DO
      IF (gapw_xc) THEN
         ALLOCATE (trho_xc_r(nspins), trho_xc_g(nspins))
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%create_pw(trho_xc_r(ispin), &
                                          use_data=REALDATA3D, in_space=REALSPACE)
            CALL auxbas_pw_pool%create_pw(trho_xc_g(ispin), &
                                          use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         END DO
      END IF

      ! GAPW/GAPW_XC initializations
      NULLIFY (hartree_local, local_rho_set)
      IF (gapw) THEN
         CALL get_qs_env(qs_env, &
                         atomic_kind_set=atomic_kind_set, &
                         natom=natom, &
                         qs_kind_set=qs_kind_set)
         CALL local_rho_set_create(local_rho_set)
         CALL allocate_rho_atom_internals(local_rho_set%rho_atom_set, atomic_kind_set, &
                                          qs_kind_set, dft_control, para_env)
         CALL init_rho0(local_rho_set, qs_env, dft_control%qs_control%gapw_control, &
                        zcore=0.0_dp)
         CALL rho0_s_grid_create(pw_env, local_rho_set%rho0_mpole)
         CALL hartree_local_create(hartree_local)
         CALL init_coulomb_local(hartree_local, natom)
      ELSEIF (gapw_xc) THEN
         CALL get_qs_env(qs_env, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set)
         CALL local_rho_set_create(local_rho_set)
         CALL allocate_rho_atom_internals(local_rho_set%rho_atom_set, atomic_kind_set, &
                                          qs_kind_set, dft_control, para_env)
      END IF

      total_rho = 0.0_dp
      CALL pw_zero(rho_tot_gspace)
      DO ispin = 1, nspins
         CALL calculate_rho_elec(ks_env=ks_env, matrix_p=matrix_pe(ispin)%matrix, &
                                 rho=trho_r(ispin), &
                                 rho_gspace=trho_g(ispin), &
                                 soft_valid=gapw, &
                                 total_rho=total_rho(ispin))
         CALL pw_axpy(trho_g(ispin), rho_tot_gspace)
         IF (gapw_xc) THEN
            CALL calculate_rho_elec(ks_env=ks_env, matrix_p=matrix_pe(ispin)%matrix, &
                                    rho=trho_xc_r(ispin), &
                                    rho_gspace=trho_xc_g(ispin), &
                                    soft_valid=gapw_xc, &
                                    total_rho=rhotot)
         END IF
      END DO

      ! GAPW o GAPW_XC require the calculation of hard and soft local densities
      IF (gapw .OR. gapw_xc) THEN
         CALL get_qs_env(qs_env=qs_env, oce=oce, sab_orb=sab)
         CALL calculate_rho_atom_coeff(qs_env, matrix_pe, local_rho_set%rho_atom_set, &
                                       qs_kind_set, oce, sab, para_env)
         CALL prepare_gapw_den(qs_env, local_rho_set, do_rho0=gapw)
      END IF
      rhotot = SUM(total_rho)
      IF (gapw) THEN
         CALL get_rho0_mpole(local_rho_set%rho0_mpole, Qlm_tot=Qlm_tot)
         rhotot = rhotot + local_rho_set%rho0_mpole%total_rho0_h
         CALL pw_axpy(local_rho_set%rho0_mpole%rho0_s_gs, rho_tot_gspace)
      END IF

      IF (ABS(rhotot) > 1.e-05_dp) THEN
         logger => cp_get_default_logger()
         IF (logger%para_env%is_source()) THEN
            iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
         ELSE
            iounit = -1
         END IF
         CPWARN("Real space electron count of excitation density is non-zero.")
         IF (iounit > 0) THEN
            WRITE (iounit, "(T2,A,T61,G20.10)") "Measured electron count is ", rhotot
            WRITE (iounit, "(T2,A,/)") REPEAT("*", 79)
         END IF
      END IF

      ! calculate associated hartree potential
      CALL pw_poisson_solve(poisson_env, rho_tot_gspace, thartree, &
                            v_hartree_gspace)
      CALL pw_transfer(v_hartree_gspace, v_hartree_rspace)
      CALL pw_scale(v_hartree_rspace, v_hartree_rspace%pw_grid%dvol)
      IF (gapw) THEN
         CALL Vh_1c_gg_integrals(qs_env, thartree, hartree_local%ecoul_1c, &
                                 local_rho_set, para_env, tddft=.TRUE.)
         CALL integrate_vhg0_rspace(qs_env, v_hartree_rspace, para_env, &
                                    calculate_forces=.FALSE., &
                                    local_rho_set=local_rho_set)
      END IF

      ! Fxc*drho term
      CALL get_qs_env(qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g)
      !
      CALL get_qs_env(qs_env, input=input)
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      END IF
      !
      deriv2_analytic = section_get_lval(xc_section, "2ND_DERIV_ANALYTICAL")
      IF (deriv2_analytic) THEN
         NULLIFY (v_xc, v_xc_tau, tau_r)
         IF (gapw_xc) THEN
            CALL get_qs_env(qs_env=qs_env, rho_xc=rho_xc)
            CALL qs_fxc_analytic(rho_xc, trho_xc_r, tau_r, xc_section, auxbas_pw_pool, .FALSE., v_xc, v_xc_tau)
         ELSE
            CALL qs_fxc_analytic(rho, trho_r, tau_r, xc_section, auxbas_pw_pool, .FALSE., v_xc, v_xc_tau)
         END IF
         IF (gapw .OR. gapw_xc) THEN
            CALL get_qs_env(qs_env, rho_atom_set=rho_atom_set)
            rho1_atom_set => local_rho_set%rho_atom_set
            CALL calculate_xc_2nd_deriv_atom(rho_atom_set, rho1_atom_set, qs_env, xc_section, para_env, &
                                             do_tddft=.TRUE., do_triplet=.FALSE.)
         END IF
      ELSE
         CPABORT("NYA 00006")
         NULLIFY (v_xc, trho)
         ALLOCATE (trho)
         CALL qs_rho_create(trho)
         CALL qs_rho_set(trho, rho_r=trho_r, rho_g=trho_g)
         CALL qs_fxc_fdiff(ks_env, rho, trho, xc_section, 6, .FALSE., v_xc, v_xc_tau)
         DEALLOCATE (trho)
      END IF

      DO ispin = 1, nspins
         CALL dbcsr_set(matrix_hz(ispin)%matrix, 0.0_dp)
         CALL pw_scale(v_xc(ispin), v_xc(ispin)%pw_grid%dvol)
      END DO
      IF (gapw_xc) THEN
         DO ispin = 1, nspins
            CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_hartree_rspace, &
                                    hmat=matrix_hz(ispin), &
                                    calculate_forces=.FALSE.)
            CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                    hmat=matrix_hz(ispin), &
                                    gapw=gapw_xc, calculate_forces=.FALSE.)
         END DO
      ELSE
         ! vtot = v_xc(ispin) + v_hartree
         DO ispin = 1, nspins
            CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                    hmat=matrix_hz(ispin), &
                                    gapw=gapw, calculate_forces=.FALSE.)
            CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_hartree_rspace, &
                                    hmat=matrix_hz(ispin), &
                                    gapw=gapw, calculate_forces=.FALSE.)
         END DO
      END IF
      IF (gapw .OR. gapw_xc) THEN
         mhz(1:nspins, 1:1) => matrix_hz(1:nspins)
         mpe(1:nspins, 1:1) => matrix_pe(1:nspins)
         CALL update_ks_atom(qs_env, mhz, mpe, forces=.FALSE., &
                             rho_atom_external=local_rho_set%rho_atom_set)
      END IF

      CALL auxbas_pw_pool%give_back_pw(v_hartree_gspace)
      CALL auxbas_pw_pool%give_back_pw(v_hartree_rspace)
      CALL auxbas_pw_pool%give_back_pw(rho_tot_gspace)
      DO ispin = 1, nspins
         CALL auxbas_pw_pool%give_back_pw(trho_r(ispin))
         CALL auxbas_pw_pool%give_back_pw(trho_g(ispin))
         CALL auxbas_pw_pool%give_back_pw(v_xc(ispin))
      END DO
      DEALLOCATE (trho_r, trho_g, v_xc)
      IF (gapw_xc) THEN
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(trho_xc_r(ispin))
            CALL auxbas_pw_pool%give_back_pw(trho_xc_g(ispin))
         END DO
         DEALLOCATE (trho_xc_r, trho_xc_g)
      END IF
      IF (ASSOCIATED(v_xc_tau)) THEN
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(v_xc_tau(ispin))
         END DO
         DEALLOCATE (v_xc_tau)
      END IF
      IF (dft_control%do_admm) THEN
         IF (qs_env%admm_env%aux_exch_func == do_admm_aux_exch_func_none) THEN
            ! nothing to do
         ELSE
            ! add ADMM xc_section_aux terms: f_x[rhoz_ADMM]
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CALL get_admm_env(admm_env, rho_aux_fit=rho_aux_fit, matrix_s_aux_fit=msaux, &
                              task_list_aux_fit=task_list)
            basis_type = "AUX_FIT"
            !
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
            !
            ! GAPW/GAPW_XC initializations
            NULLIFY (local_rho_set_admm)
            IF (admm_env%do_gapw) THEN
               basis_type = "AUX_FIT_SOFT"
               task_list => admm_env%admm_gapw_env%task_list
               CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
               CALL get_admm_env(admm_env, sab_aux_fit=sab_aux_fit)
               CALL local_rho_set_create(local_rho_set_admm)
               CALL allocate_rho_atom_internals(local_rho_set_admm%rho_atom_set, atomic_kind_set, &
                                                admm_env%admm_gapw_env%admm_kind_set, dft_control, para_env)
               CALL calculate_rho_atom_coeff(qs_env, matrix_pe_admm, &
                                             rho_atom_set=local_rho_set_admm%rho_atom_set, &
                                             qs_kind_set=admm_env%admm_gapw_env%admm_kind_set, &
                                             oce=admm_env%admm_gapw_env%oce, sab=sab_aux_fit, para_env=para_env)
               CALL prepare_gapw_den(qs_env, local_rho_set=local_rho_set_admm, &
                                     do_rho0=.FALSE., kind_set_external=admm_env%admm_gapw_env%admm_kind_set)
            END IF
            !
            xc_section => admm_env%xc_section_aux
            !
            NULLIFY (rho_g_aux, rho_r_aux, rhoz_g_aux, rhoz_r_aux)
            CALL qs_rho_get(rho_aux_fit, rho_r=rho_r_aux, rho_g=rho_g_aux)
            ! rhoz_aux
            ALLOCATE (rhoz_r_aux(nspins), rhoz_g_aux(nspins))
            DO ispin = 1, nspins
               CALL auxbas_pw_pool%create_pw(rhoz_r_aux(ispin), &
                                             use_data=REALDATA3D, in_space=REALSPACE)
               CALL auxbas_pw_pool%create_pw(rhoz_g_aux(ispin), &
                                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            END DO
            DO ispin = 1, nspins
               CALL calculate_rho_elec(ks_env=ks_env, matrix_p=mpe(ispin, 1)%matrix, &
                                       rho=rhoz_r_aux(ispin), rho_gspace=rhoz_g_aux(ispin), &
                                       basis_type=basis_type, &
                                       task_list_external=task_list)
            END DO
            !
            NULLIFY (v_xc)
            deriv2_analytic = section_get_lval(xc_section, "2ND_DERIV_ANALYTICAL")
            IF (deriv2_analytic) THEN
               NULLIFY (tau_r)
               CALL qs_fxc_analytic(rho_aux_fit, rhoz_r_aux, tau_r, xc_section, auxbas_pw_pool, .FALSE., v_xc, v_xc_tau)
            ELSE
               CPABORT("NYA 00007")
               NULLIFY (rhoz_aux)
               ALLOCATE (rhoz_aux)
               CALL qs_rho_create(rhoz_aux)
               CALL qs_rho_set(rhoz_aux, rho_r=rhoz_r_aux, rho_g=rhoz_g_aux)
               CALL qs_fxc_fdiff(ks_env, rho_aux_fit, rhoz_aux, xc_section, 6, .FALSE., v_xc, v_xc_tau)
               DEALLOCATE (rhoz_aux)
            END IF
            !
            DO ispin = 1, nspins
               CALL pw_scale(v_xc(ispin), v_xc(ispin)%pw_grid%dvol)
               CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                       hmat=mhz(ispin, 1), basis_type=basis_type, &
                                       calculate_forces=.FALSE., &
                                       task_list_external=task_list)
            END DO
            DO ispin = 1, nspins
               CALL auxbas_pw_pool%give_back_pw(v_xc(ispin))
               CALL auxbas_pw_pool%give_back_pw(rhoz_r_aux(ispin))
               CALL auxbas_pw_pool%give_back_pw(rhoz_g_aux(ispin))
            END DO
            DEALLOCATE (v_xc, rhoz_r_aux, rhoz_g_aux)
            !
            IF (admm_env%do_gapw) THEN
               rho_atom_set => admm_env%admm_gapw_env%local_rho_set%rho_atom_set
               rho1_atom_set => local_rho_set_admm%rho_atom_set
               CALL calculate_xc_2nd_deriv_atom(rho_atom_set, rho1_atom_set, qs_env, xc_section, &
                                                para_env, kind_set_external=admm_env%admm_gapw_env%admm_kind_set)
               CALL update_ks_atom(qs_env, mhz(:, 1), matrix_pe_admm, forces=.FALSE., tddft=.FALSE., &
                                   rho_atom_external=rho1_atom_set, &
                                   kind_set_external=admm_env%admm_gapw_env%admm_kind_set, &
                                   oce_external=admm_env%admm_gapw_env%oce, &
                                   sab_external=sab_aux_fit)
            END IF
            !
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
               CALL dbcsr_copy(dbwork, matrix_hz(1)%matrix)
               CALL dbcsr_set(dbwork, 0.0_dp)
               CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbwork, keep_sparsity=.TRUE.)
               CALL dbcsr_add(matrix_hz(ispin)%matrix, dbwork, 1.0_dp, 1.0_dp)
            END DO
            CALL dbcsr_release(dbwork)
            DEALLOCATE (dbwork)
            CALL dbcsr_deallocate_matrix_set(mhz)
            DEALLOCATE (mpe)
            IF (admm_env%do_gapw) THEN
               IF (ASSOCIATED(local_rho_set_admm)) CALL local_rho_set_release(local_rho_set_admm)
            END IF
         END IF
      END IF
      IF (gapw .OR. gapw_xc) THEN
         IF (ASSOCIATED(local_rho_set)) CALL local_rho_set_release(local_rho_set)
         IF (ASSOCIATED(hartree_local)) CALL hartree_local_release(hartree_local)
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
         CALL get_qs_env(qs_env=qs_env, rho=rho, x_data=x_data, para_env=para_env, &
                         s_mstruct_changed=s_mstruct_changed)
         distribute_fock_matrix = .TRUE.
         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CALL get_admm_env(admm_env, matrix_s_aux_fit=msaux)
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
            IF (x_data(1, 1)%do_hfx_ri) THEN
               eh1 = 0.0_dp
               CALL hfx_ri_update_ks(qs_env, x_data(1, 1)%ri_data, mhz, eh1, rho_ao=mpe, &
                                     geometry_did_change=s_mstruct_changed, nspins=nspins, &
                                     hf_fraction=x_data(1, 1)%general_parameter%fraction)
            ELSE
               DO ispin = 1, mspin
                  eh1 = 0.0
                  CALL integrate_four_center(qs_env, x_data, mhz, eh1, mpe, hfx_section, &
                                             para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                             ispin=ispin)
               END DO
            END IF
            !
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
            CALL dbcsr_deallocate_matrix_set(mhz)
            DEALLOCATE (mpe)
         ELSE
            NULLIFY (mpe, mhz)
            ALLOCATE (mpe(nspins, 1), mhz(nspins, 1))
            DO ispin = 1, nspins
               mhz(ispin, 1)%matrix => matrix_hz(ispin)%matrix
               mpe(ispin, 1)%matrix => matrix_pe(ispin)%matrix
            END DO
            IF (x_data(1, 1)%do_hfx_ri) THEN
               eh1 = 0.0_dp
               CALL hfx_ri_update_ks(qs_env, x_data(1, 1)%ri_data, mhz, eh1, rho_ao=mpe, &
                                     geometry_did_change=s_mstruct_changed, nspins=nspins, &
                                     hf_fraction=x_data(1, 1)%general_parameter%fraction)
            ELSE
               DO ispin = 1, mspin
                  eh1 = 0.0
                  CALL integrate_four_center(qs_env, x_data, mhz, eh1, mpe, hfx_section, &
                                             para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                             ispin=ispin)
               END DO
            END IF
            DEALLOCATE (mpe, mhz)
         END IF
      END IF

      focc = 4.0_dp
      IF (nspins == 2) focc = 2.0_dp
      DO ispin = 1, nspins
         mos => gs_mos(ispin)%mos_occ
         CALL cp_fm_get_info(mos, ncol_global=norb)
         CALL cp_dbcsr_sm_fm_multiply(matrix_hz(ispin)%matrix, mos, cpmos(ispin), &
                                      norb, alpha=focc, beta=0.0_dp)
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt_resvec2

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_pe ...
!> \param gs_mos ...
!> \param matrix_hz ...
!> \param cpmos ...
! **************************************************************************************************
   SUBROUTINE tddfpt_resvec2_xtb(qs_env, matrix_pe, gs_mos, matrix_hz, cpmos)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pe
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: cpmos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_resvec2_xtb'

      INTEGER                                            :: atom_a, handle, iatom, ikind, is, ispin, &
                                                            na, natom, natorb, nkind, norb, ns, &
                                                            nsgf, nspins
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(5)                              :: occ
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: mcharge, mcharge1
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: aocg, aocg1, charges, charges1
      REAL(KIND=dp)                                      :: focc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_type), POINTER                          :: mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_matrix
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dbcsr_type), POINTER                          :: s_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(matrix_pe))

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      nspins = dft_control%nspins

      DO ispin = 1, nspins
         CALL dbcsr_set(matrix_hz(ispin)%matrix, 0.0_dp)
      END DO

      IF (dft_control%qs_control%xtb_control%coulomb_interaction) THEN
         ! Mulliken charges
         CALL get_qs_env(qs_env, rho=rho, particle_set=particle_set, &
                         matrix_s_kp=matrix_s, para_env=para_env)
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
         CALL ao_charges(matrix_pe, s_matrix, aocg1, para_env)
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
         !
         DEALLOCATE (charges, mcharge, charges1, mcharge1)
      END IF

      focc = 2.0_dp
      IF (nspins == 2) focc = 1.0_dp
      DO ispin = 1, nspins
         mos => gs_mos(ispin)%mos_occ
         CALL cp_fm_get_info(mos, ncol_global=norb)
         CALL cp_dbcsr_sm_fm_multiply(matrix_hz(ispin)%matrix, mos, cpmos(ispin), &
                                      norb, alpha=focc, beta=0.0_dp)
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt_resvec2_xtb

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cpmos ...
!> \param work ...
! **************************************************************************************************
   SUBROUTINE tddfpt_resvec3(qs_env, cpmos, work)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: cpmos
      TYPE(tddfpt_work_matrices)                         :: work

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt_resvec3'

      INTEGER                                            :: handle, ispin, nao, norb, nspins
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type)                                   :: cvec, umat
      TYPE(cp_fm_type), POINTER                          :: omos
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, mos=mos, dft_control=dft_control)
      nspins = dft_control%nspins

      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin), mo_coeff=omos)
         ASSOCIATE (rvecs => cpmos(ispin))
            CALL cp_fm_get_info(rvecs, nrow_global=nao, ncol_global=norb)
            CALL cp_fm_create(cvec, rvecs%matrix_struct, "cvec")
            CALL cp_fm_struct_create(fmstruct, context=rvecs%matrix_struct%context, nrow_global=norb, &
                                     ncol_global=norb, para_env=rvecs%matrix_struct%para_env)
            CALL cp_fm_create(umat, fmstruct, "umat")
            CALL cp_fm_struct_release(fmstruct)
            !
            CALL parallel_gemm("T", "N", norb, norb, nao, 1.0_dp, omos, work%S_C0(ispin), 0.0_dp, umat)
            CALL cp_fm_copy_general(rvecs, cvec, rvecs%matrix_struct%para_env)
            CALL parallel_gemm("N", "T", nao, norb, norb, 1.0_dp, cvec, umat, 0.0_dp, rvecs)
         END ASSOCIATE
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
      TYPE(kernel_env_type), INTENT(IN)                  :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work_matrices
      LOGICAL, INTENT(IN)                                :: debug_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_kernel_force'

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
