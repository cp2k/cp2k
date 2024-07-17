!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contains ADMM methods which require molecular orbitals
!> \par History
!>      04.2008 created [Manuel Guidon]
!>      12.2019 Made GAPW compatible [A. Bussy]
!> \author Manuel Guidon
! **************************************************************************************************
MODULE admm_methods
   USE admm_types,                      ONLY: admm_gapw_r3d_rs_type,&
                                              admm_type,&
                                              get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Merlot2014,&
                                              cite_reference
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_decompose,&
                                              cp_cfm_cholesky_invert,&
                                              cp_cfm_scale,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_scale_and_add_fm,&
                                              cp_cfm_transpose
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, &
        dbcsr_dot, dbcsr_get_block_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric, &
        dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_schur_product,&
                                              cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert,&
                                              cp_fm_cholesky_reduce,&
                                              cp_fm_cholesky_restore
   USE cp_fm_diag,                      ONLY: cp_fm_syevd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        copy_info_type, cp_fm_cleanup_copy_general, cp_fm_create, cp_fm_finish_copy_general, &
        cp_fm_get_info, cp_fm_release, cp_fm_set_all, cp_fm_set_element, cp_fm_start_copy_general, &
        cp_fm_to_fm, cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE input_constants,                 ONLY: do_admm_purify_cauchy,&
                                              do_admm_purify_cauchy_subspace,&
                                              do_admm_purify_mo_diag,&
                                              do_admm_purify_mo_no_diag,&
                                              do_admm_purify_none
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: kpoint_density_matrices,&
                                              kpoint_density_transform,&
                                              rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_env,&
                                              get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: add_qs_force,&
                                              qs_force_type
   USE qs_gapw_densities,               ONLY: prepare_gapw_den
   USE qs_ks_atom,                      ONLY: update_ks_atom
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_local_rho_types,              ONLY: local_rho_set_create,&
                                              local_rho_set_release,&
                                              local_rho_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_force
   USE qs_rho_atom_methods,             ONLY: allocate_rho_atom_internals,&
                                              calculate_rho_atom_coeff
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE qs_vxc_atom,                     ONLY: calculate_vxc_atom
   USE task_list_types,                 ONLY: task_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: admm_mo_calc_rho_aux, &
             admm_mo_calc_rho_aux_kp, &
             admm_mo_merge_ks_matrix, &
             admm_mo_merge_derivs, &
             admm_aux_response_density, &
             calc_mixed_overlap_force, &
             scale_dm, &
             admm_fit_mo_coeffs, &
             admm_update_ks_atom, &
             calc_admm_mo_derivatives, &
             calc_admm_ovlp_forces, &
             calc_admm_ovlp_forces_kp, &
             admm_projection_derivative, &
             kpoint_calc_admm_matrices

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_methods'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE admm_mo_calc_rho_aux(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'admm_mo_calc_rho_aux'

      CHARACTER(LEN=default_string_length)               :: basis_type
      INTEGER                                            :: handle, ispin
      LOGICAL                                            :: gapw, s_mstruct_changed
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r_aux
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_s_aux_fit, &
                                                            matrix_s_aux_fit_vs_orb, rho_ao, &
                                                            rho_ao_aux
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos, mos_aux_fit
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g_aux
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r_aux
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(task_list_type), POINTER                      :: task_list

      CALL timeset(routineN, handle)

      NULLIFY (ks_env, admm_env, mos, mos_aux_fit, matrix_s_aux_fit, &
               matrix_s_aux_fit_vs_orb, matrix_s, rho, rho_aux_fit, para_env)
      NULLIFY (rho_g_aux, rho_r_aux, rho_ao, rho_ao_aux, tot_rho_r_aux, task_list)

      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      admm_env=admm_env, &
                      dft_control=dft_control, &
                      mos=mos, &
                      matrix_s=matrix_s, &
                      para_env=para_env, &
                      s_mstruct_changed=s_mstruct_changed, &
                      rho=rho)
      CALL get_admm_env(admm_env, mos_aux_fit=mos_aux_fit, matrix_s_aux_fit=matrix_s_aux_fit, &
                        matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb, rho_aux_fit=rho_aux_fit)

      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CALL qs_rho_get(rho_aux_fit, &
                      rho_ao=rho_ao_aux, &
                      rho_g=rho_g_aux, &
                      rho_r=rho_r_aux, &
                      tot_rho_r=tot_rho_r_aux)

      gapw = admm_env%do_gapw

      ! convert mos from full to dbcsr matrices
      DO ispin = 1, dft_control%nspins
         IF (mos(ispin)%use_mo_coeff_b) THEN
            CALL copy_dbcsr_to_fm(mos(ispin)%mo_coeff_b, mos(ispin)%mo_coeff)
         END IF
      END DO

      ! fit mo coeffcients
      CALL admm_fit_mo_coeffs(admm_env, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, &
                              mos, mos_aux_fit, s_mstruct_changed)

      DO ispin = 1, dft_control%nspins
         IF (admm_env%block_dm) THEN
            CALL blockify_density_matrix(admm_env, &
                                         density_matrix=rho_ao(ispin)%matrix, &
                                         density_matrix_aux=rho_ao_aux(ispin)%matrix, &
                                         ispin=ispin, &
                                         nspins=dft_control%nspins)

         ELSE

            ! Here, the auxiliary DM gets calculated and is written into rho_aux_fit%...
            CALL calculate_dm_mo_no_diag(admm_env, &
                                         mo_set=mos(ispin), &
                                         overlap_matrix=matrix_s_aux_fit(1)%matrix, &
                                         density_matrix=rho_ao_aux(ispin)%matrix, &
                                         overlap_matrix_large=matrix_s(1)%matrix, &
                                         density_matrix_large=rho_ao(ispin)%matrix, &
                                         ispin=ispin)

         END IF

         IF (admm_env%purification_method == do_admm_purify_cauchy) &
            CALL purify_dm_cauchy(admm_env, &
                                  mo_set=mos_aux_fit(ispin), &
                                  density_matrix=rho_ao_aux(ispin)%matrix, &
                                  ispin=ispin, &
                                  blocked=admm_env%block_dm)

         !GPW is the default, PW density is computed using the AUX_FIT basis and task_list
         !If GAPW, the we use the AUX_FIT_SOFT basis and task list
         basis_type = "AUX_FIT"
         task_list => admm_env%task_list_aux_fit
         IF (gapw) THEN
            basis_type = "AUX_FIT_SOFT"
            task_list => admm_env%admm_gapw_env%task_list
         END IF

         CALL calculate_rho_elec(ks_env=ks_env, &
                                 matrix_p=rho_ao_aux(ispin)%matrix, &
                                 rho=rho_r_aux(ispin), &
                                 rho_gspace=rho_g_aux(ispin), &
                                 total_rho=tot_rho_r_aux(ispin), &
                                 soft_valid=.FALSE., &
                                 basis_type=basis_type, &
                                 task_list_external=task_list)

      END DO

      !If GAPW, also need to prepare the atomic densities
      IF (gapw) THEN

         CALL calculate_rho_atom_coeff(qs_env, rho_ao_aux, &
                                       rho_atom_set=admm_env%admm_gapw_env%local_rho_set%rho_atom_set, &
                                       qs_kind_set=admm_env%admm_gapw_env%admm_kind_set, &
                                       oce=admm_env%admm_gapw_env%oce, sab=admm_env%sab_aux_fit, para_env=para_env)

         CALL prepare_gapw_den(qs_env, local_rho_set=admm_env%admm_gapw_env%local_rho_set, &
                               do_rho0=.FALSE., kind_set_external=admm_env%admm_gapw_env%admm_kind_set)
      END IF

      IF (dft_control%nspins == 1) THEN
         admm_env%gsi(3) = admm_env%gsi(1)
      ELSE
         admm_env%gsi(3) = (admm_env%gsi(1) + admm_env%gsi(2))/2.0_dp
      END IF

      CALL qs_rho_set(rho_aux_fit, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)

      CALL timestop(handle)

   END SUBROUTINE admm_mo_calc_rho_aux

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE admm_mo_calc_rho_aux_kp(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'admm_mo_calc_rho_aux_kp'

      CHARACTER(LEN=default_string_length)               :: basis_type
      INTEGER                                            :: handle, i, igroup, ik, ikp, img, indx, &
                                                            ispin, kplocal, nao_aux_fit, nao_orb, &
                                                            natom, nkp, nkp_groups, nmo, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: gapw, my_kpgrp, pmat_from_rs, &
                                                            use_real_wfn
      REAL(dp)                                           :: maxval_mos, nelec_aux(2), nelec_orb(2), &
                                                            tmp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: occ_num, occ_num_aux, tot_rho_r_aux
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(copy_info_type), ALLOCATABLE, DIMENSION(:, :) :: info
      TYPE(cp_cfm_type)                                  :: cA, cmo_coeff, cmo_coeff_aux_fit, &
                                                            cpmatrix, cwork_aux_aux, cwork_aux_orb
      TYPE(cp_fm_struct_type), POINTER                   :: mo_struct, mo_struct_aux_fit, &
                                                            struct_aux_aux, struct_aux_orb, &
                                                            struct_orb_orb
      TYPE(cp_fm_type)                                   :: fmdummy, work_aux_orb, work_orb_orb, &
                                                            work_orb_orb2
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_aux_fit
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_s_aux_fit, rho_ao_aux, &
                                                            rho_ao_orb
      TYPE(dbcsr_type)                                   :: pmatrix_tmp
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: pmatrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos, mos_aux_fit
      TYPE(mo_set_type), DIMENSION(:, :), POINTER        :: mos_aux_fit_kp, mos_kp
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux_fit, sab_kp
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g_aux
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r_aux
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_aux_fit, rho_orb
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(task_list_type), POINTER                      :: task_list

      CALL timeset(routineN, handle)

      NULLIFY (ks_env, admm_env, mos, mos_aux_fit, matrix_s, rho_orb, &
               matrix_s_aux_fit, rho_aux_fit, rho_ao_orb, &
               para_env, rho_g_aux, rho_r_aux, rho_ao_aux, tot_rho_r_aux, &
               kpoints, sab_aux_fit, sab_kp, kp, &
               struct_orb_orb, struct_aux_orb, struct_aux_aux, mo_struct, mo_struct_aux_fit)

      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      admm_env=admm_env, &
                      dft_control=dft_control, &
                      kpoints=kpoints, &
                      natom=natom, &
                      scf_env=scf_env, &
                      matrix_s_kp=matrix_s, &
                      rho=rho_orb)
      CALL get_admm_env(admm_env, &
                        rho_aux_fit=rho_aux_fit, &
                        matrix_s_aux_fit_kp=matrix_s_aux_fit, &
                        sab_aux_fit=sab_aux_fit)
      gapw = admm_env%do_gapw

      CALL qs_rho_get(rho_aux_fit, &
                      rho_ao_kp=rho_ao_aux, &
                      rho_g=rho_g_aux, &
                      rho_r=rho_r_aux, &
                      tot_rho_r=tot_rho_r_aux)

      CALL qs_rho_get(rho_orb, rho_ao_kp=rho_ao_orb)
      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, use_real_wfn=use_real_wfn, kp_range=kp_range, &
                           nkp_groups=nkp_groups, kp_dist=kp_dist, &
                           cell_to_index=cell_to_index, sab_nl=sab_kp)

      ! the temporary DBCSR matrices for the rskp_transform we have to manually allocate
      ! index 1 => real, index 2 => imaginary
      ALLOCATE (pmatrix(2))
      CALL dbcsr_create(pmatrix(1), template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(pmatrix(2), template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL dbcsr_create(pmatrix_tmp, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(pmatrix(1), sab_kp)
      CALL cp_dbcsr_alloc_block_from_nbl(pmatrix(2), sab_kp)

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb
      nspins = dft_control%nspins

      !Create fm and cfm work matrices, for each KP subgroup
      CALL cp_fm_struct_create(struct_orb_orb, context=kpoints%blacs_env, para_env=kpoints%para_env_kp, &
                               nrow_global=nao_orb, ncol_global=nao_orb)
      CALL cp_fm_create(work_orb_orb, struct_orb_orb)
      CALL cp_fm_create(work_orb_orb2, struct_orb_orb)

      CALL cp_fm_struct_create(struct_aux_aux, context=kpoints%blacs_env, para_env=kpoints%para_env_kp, &
                               nrow_global=nao_aux_fit, ncol_global=nao_aux_fit)

      CALL cp_fm_struct_create(struct_aux_orb, context=kpoints%blacs_env, para_env=kpoints%para_env_kp, &
                               nrow_global=nao_aux_fit, ncol_global=nao_orb)
      CALL cp_fm_create(work_aux_orb, struct_orb_orb)

      IF (.NOT. use_real_wfn) THEN
         CALL cp_cfm_create(cpmatrix, struct_orb_orb)

         CALL cp_cfm_create(cwork_aux_aux, struct_aux_aux)

         CALL cp_cfm_create(cA, struct_aux_orb)
         CALL cp_cfm_create(cwork_aux_orb, struct_aux_orb)

         CALL get_kpoint_env(kpoints%kp_env(1)%kpoint_env, mos=mos_kp)
         mos => mos_kp(1, :)
         CALL get_mo_set(mos(1), mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, matrix_struct=mo_struct)
         CALL cp_cfm_create(cmo_coeff, mo_struct)

         CALL get_kpoint_env(kpoints%kp_aux_env(1)%kpoint_env, mos=mos_aux_fit_kp)
         mos => mos_aux_fit_kp(1, :)
         CALL get_mo_set(mos(1), mo_coeff=mo_coeff_aux_fit)
         CALL cp_fm_get_info(mo_coeff_aux_fit, matrix_struct=mo_struct_aux_fit)
         CALL cp_cfm_create(cmo_coeff_aux_fit, mo_struct_aux_fit)
      END IF

      CALL cp_fm_struct_release(struct_orb_orb)
      CALL cp_fm_struct_release(struct_aux_aux)
      CALL cp_fm_struct_release(struct_aux_orb)

      para_env => kpoints%blacs_env_all%para_env
      kplocal = kp_range(2) - kp_range(1) + 1

      !We querry the maximum absolute value of the KP MOs to see if they are populated at all. If not, we
      !need to get the KP Pmat from the RS ones (happens at first SCF step, for example)
      maxval_mos = 0.0_dp
      indx = 0
      DO ikp = 1, kplocal
         DO ispin = 1, nspins
            DO igroup = 1, nkp_groups
               ! number of current kpoint
               ik = kp_dist(1, igroup) + ikp - 1
               my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
               indx = indx + 1

               CALL get_kpoint_env(kpoints%kp_env(ikp)%kpoint_env, mos=mos_kp)
               mos => mos_kp(1, :)
               CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
               maxval_mos = MAX(maxval_mos, MAXVAL(ABS(mo_coeff%local_data)))

               IF (.NOT. use_real_wfn) THEN
                  mos => mos_kp(2, :)
                  CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
                  maxval_mos = MAX(maxval_mos, MAXVAL(ABS(mo_coeff%local_data)))
               END IF
            END DO
         END DO
      END DO
      CALL para_env%sum(maxval_mos) !I think para_env is the global one

      pmat_from_rs = .FALSE.
      IF (maxval_mos < EPSILON(0.0_dp)) pmat_from_rs = .TRUE.

      !TODO: issue a warning when doing ADMM with ATOMIC guess. If small number of K-points => leads to bad things

      ALLOCATE (info(kplocal*nspins*nkp_groups, 2))
      !Start communication: only P matrix, and only if required
      indx = 0
      IF (pmat_from_rs) THEN
         DO ikp = 1, kplocal
            DO ispin = 1, nspins
               DO igroup = 1, nkp_groups
                  ! number of current kpoint
                  ik = kp_dist(1, igroup) + ikp - 1
                  my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
                  indx = indx + 1

                  ! FT of matrices P if required, then transfer to FM type
                  IF (use_real_wfn) THEN
                     CALL dbcsr_set(pmatrix(1), 0.0_dp)
                     CALL rskp_transform(rmatrix=pmatrix(1), rsmat=rho_ao_orb, ispin=ispin, &
                                         xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_kp)
                     CALL dbcsr_desymmetrize(pmatrix(1), pmatrix_tmp)
                     CALL copy_dbcsr_to_fm(pmatrix_tmp, admm_env%work_orb_orb)
                  ELSE
                     CALL dbcsr_set(pmatrix(1), 0.0_dp)
                     CALL dbcsr_set(pmatrix(2), 0.0_dp)
                     CALL rskp_transform(rmatrix=pmatrix(1), cmatrix=pmatrix(2), rsmat=rho_ao_orb, ispin=ispin, &
                                         xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_kp)
                     CALL dbcsr_desymmetrize(pmatrix(1), pmatrix_tmp)
                     CALL copy_dbcsr_to_fm(pmatrix_tmp, admm_env%work_orb_orb)
                     CALL dbcsr_desymmetrize(pmatrix(2), pmatrix_tmp)
                     CALL copy_dbcsr_to_fm(pmatrix_tmp, admm_env%work_orb_orb2)
                  END IF

                  IF (my_kpgrp) THEN
                     CALL cp_fm_start_copy_general(admm_env%work_orb_orb, work_orb_orb, para_env, info(indx, 1))
                     IF (.NOT. use_real_wfn) THEN
                        CALL cp_fm_start_copy_general(admm_env%work_orb_orb2, work_orb_orb2, para_env, info(indx, 2))
                     END IF
                  ELSE
                     CALL cp_fm_start_copy_general(admm_env%work_orb_orb, fmdummy, para_env, info(indx, 1))
                     IF (.NOT. use_real_wfn) THEN
                        CALL cp_fm_start_copy_general(admm_env%work_orb_orb2, fmdummy, para_env, info(indx, 2))
                     END IF
                  END IF !my_kpgrp
               END DO
            END DO
         END DO
      END IF !pmat_from_rs

      indx = 0
      DO ikp = 1, kplocal
         DO ispin = 1, nspins
            DO igroup = 1, nkp_groups
               ! number of current kpoint
               ik = kp_dist(1, igroup) + ikp - 1
               my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
               indx = indx + 1
               IF (my_kpgrp .AND. pmat_from_rs) THEN
                  CALL cp_fm_finish_copy_general(work_orb_orb, info(indx, 1))
                  IF (.NOT. use_real_wfn) THEN
                     CALL cp_fm_finish_copy_general(work_orb_orb2, info(indx, 2))
                     CALL cp_fm_to_cfm(work_orb_orb, work_orb_orb2, cpmatrix)
                  END IF
               END IF
            END DO

            IF (use_real_wfn) THEN

               nmo = admm_env%nmo(ispin)
               !! Each kpoint group has now information on a kpoint for which to calculate the MOS_aux
               CALL get_kpoint_env(kpoints%kp_env(ikp)%kpoint_env, mos=mos_kp)
               CALL get_kpoint_env(kpoints%kp_aux_env(ikp)%kpoint_env, mos=mos_aux_fit_kp)
               mos => mos_kp(1, :)
               mos_aux_fit => mos_aux_fit_kp(1, :)

               CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, occupation_numbers=occ_num)
               CALL get_mo_set(mos_aux_fit(ispin), mo_coeff=mo_coeff_aux_fit, &
                               occupation_numbers=occ_num_aux)

               kp => kpoints%kp_aux_env(ikp)%kpoint_env
               CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nao_orb, 1.0_dp, kp%amat(1, 1), &
                                  mo_coeff, 0.0_dp, mo_coeff_aux_fit)

               occ_num_aux(1:nmo) = occ_num(1:nmo)

               IF (pmat_from_rs) THEN
                  !We project on the AUX basis: P_aux = A * P *A^T
                  CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb, 1.0_dp, kp%amat(1, 1), &
                                     work_orb_orb, 0.0_dp, work_aux_orb)
                  CALL parallel_gemm('N', 'T', nao_aux_fit, nao_aux_fit, nao_orb, 1.0_dp, work_aux_orb, &
                                     kp%amat(1, 1), 0.0_dp, kpoints%kp_aux_env(ikp)%kpoint_env%pmat(1, ispin))
               END IF

            ELSE !complex wfn

               !construct the ORB MOs in complex format
               nmo = admm_env%nmo(ispin)
               CALL get_kpoint_env(kpoints%kp_env(ikp)%kpoint_env, mos=mos_kp)
               mos => mos_kp(1, :) !real
               CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
               CALL cp_cfm_scale_and_add_fm(z_zero, cmo_coeff, z_one, mo_coeff)
               mos => mos_kp(2, :) !complex
               CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
               CALL cp_cfm_scale_and_add_fm(z_one, cmo_coeff, gaussi, mo_coeff)

               !project
               kp => kpoints%kp_aux_env(ikp)%kpoint_env
               CALL cp_fm_to_cfm(kp%amat(1, 1), kp%amat(2, 1), cA)
               CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nao_orb, &
                                  z_one, cA, cmo_coeff, z_zero, cmo_coeff_aux_fit)

               !write result back to KP MOs
               CALL get_kpoint_env(kpoints%kp_aux_env(ikp)%kpoint_env, mos=mos_aux_fit_kp)
               mos_aux_fit => mos_aux_fit_kp(1, :)
               CALL get_mo_set(mos_aux_fit(ispin), mo_coeff=mo_coeff_aux_fit)
               CALL cp_cfm_to_fm(cmo_coeff_aux_fit, mtargetr=mo_coeff_aux_fit)
               mos_aux_fit => mos_aux_fit_kp(2, :)
               CALL get_mo_set(mos_aux_fit(ispin), mo_coeff=mo_coeff_aux_fit)
               CALL cp_cfm_to_fm(cmo_coeff_aux_fit, mtargeti=mo_coeff_aux_fit)

               DO i = 1, 2
                  mos => mos_kp(i, :)
                  CALL get_mo_set(mos(ispin), occupation_numbers=occ_num)
                  mos_aux_fit => mos_aux_fit_kp(i, :)
                  CALL get_mo_set(mos_aux_fit(ispin), occupation_numbers=occ_num_aux)
                  occ_num_aux(:) = occ_num(:)
               END DO

               IF (pmat_from_rs) THEN
                  CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb, z_one, cA, &
                                     cpmatrix, z_zero, cwork_aux_orb)
                  CALL parallel_gemm('N', 'C', nao_aux_fit, nao_aux_fit, nao_orb, z_one, cwork_aux_orb, &
                                     cA, z_zero, cwork_aux_aux)

                  CALL cp_cfm_to_fm(cwork_aux_aux, mtargetr=kpoints%kp_aux_env(ikp)%kpoint_env%pmat(1, ispin), &
                                    mtargeti=kpoints%kp_aux_env(ikp)%kpoint_env%pmat(2, ispin))
               END IF
            END IF

         END DO
      END DO

      !Clean-up communication
      IF (pmat_from_rs) THEN
         indx = 0
         DO ikp = 1, kplocal
            DO ispin = 1, nspins
               DO igroup = 1, nkp_groups
                  ! number of current kpoint
                  ik = kp_dist(1, igroup) + ikp - 1
                  my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
                  indx = indx + 1

                  CALL cp_fm_cleanup_copy_general(info(indx, 1))
                  IF (.NOT. use_real_wfn) CALL cp_fm_cleanup_copy_general(info(indx, 2))
               END DO
            END DO
         END DO
      END IF

      DEALLOCATE (info)
      CALL dbcsr_release(pmatrix(1))
      CALL dbcsr_release(pmatrix(2))
      CALL dbcsr_release(pmatrix_tmp)

      CALL cp_fm_release(work_orb_orb)
      CALL cp_fm_release(work_orb_orb2)
      CALL cp_fm_release(work_aux_orb)
      IF (.NOT. use_real_wfn) THEN
         CALL cp_cfm_release(cpmatrix)
         CALL cp_cfm_release(cwork_aux_aux)
         CALL cp_cfm_release(cwork_aux_orb)
         CALL cp_cfm_release(cA)
         CALL cp_cfm_release(cmo_coeff)
         CALL cp_cfm_release(cmo_coeff_aux_fit)
      END IF

      IF (.NOT. pmat_from_rs) CALL kpoint_density_matrices(kpoints, for_aux_fit=.TRUE.)
      CALL kpoint_density_transform(kpoints, rho_ao_aux, .FALSE., &
                                    matrix_s_aux_fit(1, 1)%matrix, sab_aux_fit, &
                                    admm_env%scf_work_aux_fit, for_aux_fit=.TRUE.)

      !ADMMQ, ADMMP, ADMMS
      IF (admm_env%do_admmq .OR. admm_env%do_admmp .OR. admm_env%do_admms) THEN

         CALL cite_reference(Merlot2014)

         nelec_orb = 0.0_dp
         nelec_aux = 0.0_dp
         admm_env%n_large_basis = 0.0_dp
         !Note: we can take the trace of the symmetric-typed matrices as P_mu^0,nu^b = P_nu^0,mu^-b
         !      and because of the sum over all images, all atomic blocks are accounted for
         DO img = 1, dft_control%nimages
            DO ispin = 1, dft_control%nspins
               CALL dbcsr_dot(rho_ao_orb(ispin, img)%matrix, matrix_s(1, img)%matrix, tmp)
               nelec_orb(ispin) = nelec_orb(ispin) + tmp
               CALL dbcsr_dot(rho_ao_aux(ispin, img)%matrix, matrix_s_aux_fit(1, img)%matrix, tmp)
               nelec_aux(ispin) = nelec_aux(ispin) + tmp
            END DO
         END DO

         DO ispin = 1, dft_control%nspins
            admm_env%n_large_basis(ispin) = nelec_orb(ispin)
            admm_env%gsi(ispin) = nelec_orb(ispin)/nelec_aux(ispin)
         END DO

         IF (admm_env%charge_constrain) THEN
            DO img = 1, dft_control%nimages
               DO ispin = 1, dft_control%nspins
                  CALL dbcsr_scale(rho_ao_aux(ispin, img)%matrix, admm_env%gsi(ispin))
               END DO
            END DO
         END IF

         IF (dft_control%nspins == 1) THEN
            admm_env%gsi(3) = admm_env%gsi(1)
         ELSE
            admm_env%gsi(3) = (admm_env%gsi(1) + admm_env%gsi(2))/2.0_dp
         END IF
      END IF

      basis_type = "AUX_FIT"
      task_list => admm_env%task_list_aux_fit
      IF (gapw) THEN
         basis_type = "AUX_FIT_SOFT"
         task_list => admm_env%admm_gapw_env%task_list
      END IF

      DO ispin = 1, nspins
         rho_ao => rho_ao_aux(ispin, :)
         CALL calculate_rho_elec(ks_env=ks_env, &
                                 matrix_p_kp=rho_ao, &
                                 rho=rho_r_aux(ispin), &
                                 rho_gspace=rho_g_aux(ispin), &
                                 total_rho=tot_rho_r_aux(ispin), &
                                 soft_valid=.FALSE., &
                                 basis_type=basis_type, &
                                 task_list_external=task_list)
      END DO

      IF (gapw) THEN
         CALL calculate_rho_atom_coeff(qs_env, rho_ao_aux, &
                                       rho_atom_set=admm_env%admm_gapw_env%local_rho_set%rho_atom_set, &
                                       qs_kind_set=admm_env%admm_gapw_env%admm_kind_set, &
                                       oce=admm_env%admm_gapw_env%oce, &
                                       sab=admm_env%sab_aux_fit, para_env=para_env)

         CALL prepare_gapw_den(qs_env, local_rho_set=admm_env%admm_gapw_env%local_rho_set, &
                               do_rho0=.FALSE., kind_set_external=admm_env%admm_gapw_env%admm_kind_set)
      END IF

      CALL qs_rho_set(rho_aux_fit, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)

      CALL timestop(handle)

   END SUBROUTINE admm_mo_calc_rho_aux_kp

! **************************************************************************************************
!> \brief Adds the GAPW exchange contribution to the aux_fit ks matrices
!> \param qs_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE admm_update_ks_atom(qs_env, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'admm_update_ks_atom'

      INTEGER                                            :: handle, img, ispin
      REAL(dp)                                           :: force_fac(2)
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_aux_fit, &
                                                            matrix_ks_aux_fit_dft, &
                                                            matrix_ks_aux_fit_hfx, rho_ao_aux
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_rho_type), POINTER                         :: rho_aux_fit

      NULLIFY (matrix_ks_aux_fit, matrix_ks_aux_fit_dft, matrix_ks_aux_fit_hfx, rho_ao_aux, rho_aux_fit)
      NULLIFY (admm_env, dft_control)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, admm_env=admm_env, dft_control=dft_control)
      CALL get_admm_env(admm_env, rho_aux_fit=rho_aux_fit, matrix_ks_aux_fit_kp=matrix_ks_aux_fit, &
                        matrix_ks_aux_fit_dft_kp=matrix_ks_aux_fit_dft, &
                        matrix_ks_aux_fit_hfx_kp=matrix_ks_aux_fit_hfx)
      CALL qs_rho_get(rho_aux_fit, rho_ao_kp=rho_ao_aux)

      !In case of ADMMS or ADMMP, need to scale the forces stemming from DFT exchagne correction
      force_fac = 1.0_dp
      IF (admm_env%do_admms) THEN
         DO ispin = 1, dft_control%nspins
            force_fac(ispin) = admm_env%gsi(ispin)**(2.0_dp/3.0_dp)
         END DO
      ELSE IF (admm_env%do_admmp) THEN
         DO ispin = 1, dft_control%nspins
            force_fac(ispin) = admm_env%gsi(ispin)**2
         END DO
      END IF

      CALL update_ks_atom(qs_env, matrix_ks_aux_fit, rho_ao_aux, calculate_forces, tddft=.FALSE., &
                          rho_atom_external=admm_env%admm_gapw_env%local_rho_set%rho_atom_set, &
                          kind_set_external=admm_env%admm_gapw_env%admm_kind_set, &
                          oce_external=admm_env%admm_gapw_env%oce, &
                          sab_external=admm_env%sab_aux_fit, fscale=force_fac)

      !Following the logic of sum_up_and_integrate to recover the pure DFT exchange contribution
      DO img = 1, dft_control%nimages
         DO ispin = 1, dft_control%nspins
            CALL dbcsr_add(matrix_ks_aux_fit_dft(ispin, img)%matrix, matrix_ks_aux_fit(ispin, img)%matrix, &
                           0.0_dp, -1.0_dp)
            CALL dbcsr_add(matrix_ks_aux_fit_dft(ispin, img)%matrix, matrix_ks_aux_fit_hfx(ispin, img)%matrix, &
                           1.0_dp, 1.0_dp)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE admm_update_ks_atom

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE admm_mo_merge_ks_matrix(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_mo_merge_ks_matrix'

      INTEGER                                            :: handle
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)
      NULLIFY (admm_env)

      CALL get_qs_env(qs_env, admm_env=admm_env, dft_control=dft_control)

      SELECT CASE (admm_env%purification_method)
      CASE (do_admm_purify_cauchy)
         CALL merge_ks_matrix_cauchy(qs_env)

      CASE (do_admm_purify_cauchy_subspace)
         CALL merge_ks_matrix_cauchy_subspace(qs_env)

      CASE (do_admm_purify_none)
         IF (dft_control%nimages > 1) THEN
            CALL merge_ks_matrix_none_kp(qs_env)
         ELSE
            CALL merge_ks_matrix_none(qs_env)
         END IF

      CASE (do_admm_purify_mo_diag, do_admm_purify_mo_no_diag)
         !do nothing
      CASE DEFAULT
         CPABORT("admm_mo_merge_ks_matrix: unknown purification method")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE admm_mo_merge_ks_matrix

! **************************************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param mo_set ...
!> \param mo_coeff ...
!> \param mo_coeff_aux_fit ...
!> \param mo_derivs ...
!> \param mo_derivs_aux_fit ...
!> \param matrix_ks_aux_fit ...
! **************************************************************************************************
   SUBROUTINE admm_mo_merge_derivs(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, mo_derivs, &
                                   mo_derivs_aux_fit, matrix_ks_aux_fit)
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(mo_set_type), INTENT(IN)                      :: mo_set
      TYPE(cp_fm_type), INTENT(IN)                       :: mo_coeff, mo_coeff_aux_fit
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_derivs, mo_derivs_aux_fit
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks_aux_fit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_mo_merge_derivs'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (admm_env%purification_method)
      CASE (do_admm_purify_mo_diag)
         CALL merge_mo_derivs_diag(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, &
                                   mo_derivs, mo_derivs_aux_fit, matrix_ks_aux_fit)

      CASE (do_admm_purify_mo_no_diag)
         CALL merge_mo_derivs_no_diag(ispin, admm_env, mo_set, mo_derivs, matrix_ks_aux_fit)

      CASE (do_admm_purify_none, do_admm_purify_cauchy, do_admm_purify_cauchy_subspace)
         !do nothing
      CASE DEFAULT
         CPABORT("admm_mo_merge_derivs: unknown purification method")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE admm_mo_merge_derivs

! **************************************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_mixed ...
!> \param mos ...
!> \param mos_aux_fit ...
!> \param geometry_did_change ...
! **************************************************************************************************
   SUBROUTINE admm_fit_mo_coeffs(admm_env, matrix_s_aux_fit, matrix_s_mixed, &
                                 mos, mos_aux_fit, geometry_did_change)

      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s_aux_fit, matrix_s_mixed
      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos, mos_aux_fit
      LOGICAL, INTENT(IN)                                :: geometry_did_change

      CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_fit_mo_coeffs'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (geometry_did_change) THEN
         CALL fit_mo_coeffs(admm_env, matrix_s_aux_fit, matrix_s_mixed)
      END IF

      SELECT CASE (admm_env%purification_method)
      CASE (do_admm_purify_mo_no_diag, do_admm_purify_cauchy_subspace)
         CALL purify_mo_cholesky(admm_env, mos, mos_aux_fit)

      CASE (do_admm_purify_mo_diag)
         CALL purify_mo_diag(admm_env, mos, mos_aux_fit)

      CASE DEFAULT
         CALL purify_mo_none(admm_env, mos, mos_aux_fit)
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE admm_fit_mo_coeffs

! **************************************************************************************************
!> \brief Calculate S^-1, Q, B full-matrices given sparse S_tilde and Q
!> \param admm_env ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_mixed ...
! **************************************************************************************************
   SUBROUTINE fit_mo_coeffs(admm_env, matrix_s_aux_fit, matrix_s_mixed)
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s_aux_fit, matrix_s_mixed

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fit_mo_coeffs'

      INTEGER                                            :: blk, handle, iatom, jatom, nao_aux_fit, &
                                                            nao_orb
      REAL(dp), DIMENSION(:, :), POINTER                 :: sparse_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type), POINTER                          :: matrix_s_tilde

      CALL timeset(routineN, handle)

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb

      ! *** This part only depends on overlap matrices ==> needs only to be calculated if the geometry changed

      IF (.NOT. admm_env%block_fit) THEN
         CALL copy_dbcsr_to_fm(matrix_s_aux_fit(1)%matrix, admm_env%S_inv)
      ELSE
         NULLIFY (matrix_s_tilde)
         ALLOCATE (matrix_s_tilde)
         CALL dbcsr_create(matrix_s_tilde, template=matrix_s_aux_fit(1)%matrix, &
                           name='MATRIX s_tilde', &
                           matrix_type=dbcsr_type_symmetric)

         CALL dbcsr_copy(matrix_s_tilde, matrix_s_aux_fit(1)%matrix)

         CALL dbcsr_iterator_start(iter, matrix_s_tilde)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
            IF (admm_env%block_map(iatom, jatom) == 0) THEN
               sparse_block = 0.0_dp
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
         CALL copy_dbcsr_to_fm(matrix_s_tilde, admm_env%S_inv)
         CALL dbcsr_deallocate_matrix(matrix_s_tilde)
      END IF

      CALL cp_fm_upper_to_full(admm_env%S_inv, admm_env%work_aux_aux)
      CALL cp_fm_to_fm(admm_env%S_inv, admm_env%S)

      CALL copy_dbcsr_to_fm(matrix_s_mixed(1)%matrix, admm_env%Q)

      !! Calculate S'_inverse
      CALL cp_fm_cholesky_decompose(admm_env%S_inv)
      CALL cp_fm_cholesky_invert(admm_env%S_inv)
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(admm_env%S_inv, admm_env%work_aux_aux)

      !! Calculate A=S'^(-1)*Q
      IF (admm_env%block_fit) THEN
         CALL cp_fm_set_all(admm_env%A, 0.0_dp, 1.0_dp)
      ELSE
         CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                            1.0_dp, admm_env%S_inv, admm_env%Q, 0.0_dp, &
                            admm_env%A)

         ! this multiplication is apparent not need for purify_none
         !! B=Q^(T)*A
         CALL parallel_gemm('T', 'N', nao_orb, nao_orb, nao_aux_fit, &
                            1.0_dp, admm_env%Q, admm_env%A, 0.0_dp, &
                            admm_env%B)
      END IF

      CALL timestop(handle)

   END SUBROUTINE fit_mo_coeffs

! **************************************************************************************************
!> \brief Calculates the MO coefficients for the auxiliary fitting basis set
!>        by minimizing int (psi_i - psi_aux_i)^2 using Lagrangian Multipliers
!>
!> \param admm_env The ADMM env
!> \param mos the MO's of the orbital basis set
!> \param mos_aux_fit the MO's of the auxiliary fitting basis set
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE purify_mo_cholesky(admm_env, mos, mos_aux_fit)

      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos, mos_aux_fit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mo_cholesky'

      INTEGER                                            :: handle, ispin, nao_aux_fit, nao_orb, &
                                                            nmo, nspins
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_aux_fit

      CALL timeset(routineN, handle)

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb
      nspins = SIZE(mos)

      ! *** Calculate the mo_coeffs for the fitting basis
      DO ispin = 1, nspins
         nmo = admm_env%nmo(ispin)
         IF (nmo == 0) CYCLE
         !! Lambda = C^(T)*B*C
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
         CALL get_mo_set(mos_aux_fit(ispin), mo_coeff=mo_coeff_aux_fit)
         CALL parallel_gemm('N', 'N', nao_orb, nmo, nao_orb, &
                            1.0_dp, admm_env%B, mo_coeff, 0.0_dp, &
                            admm_env%work_orb_nmo(ispin))
         CALL parallel_gemm('T', 'N', nmo, nmo, nao_orb, &
                            1.0_dp, mo_coeff, admm_env%work_orb_nmo(ispin), 0.0_dp, &
                            admm_env%lambda(ispin))
         CALL cp_fm_to_fm(admm_env%lambda(ispin), admm_env%work_nmo_nmo1(ispin))

         CALL cp_fm_cholesky_decompose(admm_env%work_nmo_nmo1(ispin))
         CALL cp_fm_cholesky_invert(admm_env%work_nmo_nmo1(ispin))
         !! Symmetrize the guy
         CALL cp_fm_upper_to_full(admm_env%work_nmo_nmo1(ispin), admm_env%lambda_inv(ispin))
         CALL cp_fm_to_fm(admm_env%work_nmo_nmo1(ispin), admm_env%lambda_inv(ispin))

         !! ** C_hat = AC
         CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nao_orb, &
                            1.0_dp, admm_env%A, mo_coeff, 0.0_dp, &
                            admm_env%C_hat(ispin))
         CALL cp_fm_to_fm(admm_env%C_hat(ispin), mo_coeff_aux_fit)

      END DO

      CALL timestop(handle)

   END SUBROUTINE purify_mo_cholesky

! **************************************************************************************************
!> \brief Calculates the MO coefficients for the auxiliary fitting basis set
!>        by minimizing int (psi_i - psi_aux_i)^2 using Lagrangian Multipliers
!>
!> \param admm_env The ADMM env
!> \param mos the MO's of the orbital basis set
!> \param mos_aux_fit the MO's of the auxiliary fitting basis set
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE purify_mo_diag(admm_env, mos, mos_aux_fit)

      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos, mos_aux_fit

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'purify_mo_diag'

      INTEGER                                            :: handle, i, ispin, nao_aux_fit, nao_orb, &
                                                            nmo, nspins
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: eig_work
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_aux_fit

      CALL timeset(routineN, handle)

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb
      nspins = SIZE(mos)

      ! *** Calculate the mo_coeffs for the fitting basis
      DO ispin = 1, nspins
         nmo = admm_env%nmo(ispin)
         IF (nmo == 0) CYCLE
         !! Lambda = C^(T)*B*C
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
         CALL get_mo_set(mos_aux_fit(ispin), mo_coeff=mo_coeff_aux_fit)
         CALL parallel_gemm('N', 'N', nao_orb, nmo, nao_orb, &
                            1.0_dp, admm_env%B, mo_coeff, 0.0_dp, &
                            admm_env%work_orb_nmo(ispin))
         CALL parallel_gemm('T', 'N', nmo, nmo, nao_orb, &
                            1.0_dp, mo_coeff, admm_env%work_orb_nmo(ispin), 0.0_dp, &
                            admm_env%lambda(ispin))
         CALL cp_fm_to_fm(admm_env%lambda(ispin), admm_env%work_nmo_nmo1(ispin))

         CALL cp_fm_syevd(admm_env%work_nmo_nmo1(ispin), admm_env%R(ispin), &
                          admm_env%eigvals_lambda(ispin)%eigvals%data)
         ALLOCATE (eig_work(nmo))
         DO i = 1, nmo
            eig_work(i) = 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(i))
         END DO
         CALL cp_fm_to_fm(admm_env%R(ispin), admm_env%work_nmo_nmo1(ispin))
         CALL cp_fm_column_scale(admm_env%work_nmo_nmo1(ispin), eig_work)
         CALL parallel_gemm('N', 'T', nmo, nmo, nmo, &
                            1.0_dp, admm_env%work_nmo_nmo1(ispin), admm_env%R(ispin), 0.0_dp, &
                            admm_env%lambda_inv_sqrt(ispin))
         CALL parallel_gemm('N', 'N', nao_orb, nmo, nmo, &
                            1.0_dp, mo_coeff, admm_env%lambda_inv_sqrt(ispin), 0.0_dp, &
                            admm_env%work_orb_nmo(ispin))
         CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nao_orb, &
                            1.0_dp, admm_env%A, admm_env%work_orb_nmo(ispin), 0.0_dp, &
                            mo_coeff_aux_fit)

         CALL cp_fm_to_fm(mo_coeff_aux_fit, admm_env%C_hat(ispin))
         CALL cp_fm_set_all(admm_env%lambda_inv(ispin), 0.0_dp, 1.0_dp)
         DEALLOCATE (eig_work)
      END DO

      CALL timestop(handle)

   END SUBROUTINE purify_mo_diag

! **************************************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param mos ...
!> \param mos_aux_fit ...
! **************************************************************************************************
   SUBROUTINE purify_mo_none(admm_env, mos, mos_aux_fit)
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos, mos_aux_fit

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'purify_mo_none'

      INTEGER                                            :: handle, ispin, nao_aux_fit, nao_orb, &
                                                            nmo, nmo_mos, nspins
      REAL(KIND=dp), DIMENSION(:), POINTER               :: occ_num, occ_num_aux
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_aux_fit

      CALL timeset(routineN, handle)

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb
      nspins = SIZE(mos)

      DO ispin = 1, nspins
         nmo = admm_env%nmo(ispin)
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, occupation_numbers=occ_num, nmo=nmo_mos)
         CALL get_mo_set(mos_aux_fit(ispin), mo_coeff=mo_coeff_aux_fit, &
                         occupation_numbers=occ_num_aux)

         CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nao_orb, &
                            1.0_dp, admm_env%A, mo_coeff, 0.0_dp, &
                            mo_coeff_aux_fit)
         CALL cp_fm_to_fm(mo_coeff_aux_fit, admm_env%C_hat(ispin))

         occ_num_aux(1:nmo) = occ_num(1:nmo)
         ! XXXX should only be done first time XXXX
         CALL cp_fm_set_all(admm_env%lambda(ispin), 0.0_dp, 1.0_dp)
         CALL cp_fm_set_all(admm_env%lambda_inv(ispin), 0.0_dp, 1.0_dp)
         CALL cp_fm_set_all(admm_env%lambda_inv_sqrt(ispin), 0.0_dp, 1.0_dp)
      END DO

      CALL timestop(handle)

   END SUBROUTINE purify_mo_none

! **************************************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param mo_set ...
!> \param density_matrix ...
!> \param ispin ...
!> \param blocked ...
! **************************************************************************************************
   SUBROUTINE purify_dm_cauchy(admm_env, mo_set, density_matrix, ispin, blocked)

      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(mo_set_type), INTENT(IN)                      :: mo_set
      TYPE(dbcsr_type), POINTER                          :: density_matrix
      INTEGER                                            :: ispin
      LOGICAL, INTENT(IN)                                :: blocked

      CHARACTER(len=*), PARAMETER                        :: routineN = 'purify_dm_cauchy'

      INTEGER                                            :: handle, i, nao_aux_fit, nao_orb, nmo, &
                                                            nspins
      REAL(KIND=dp)                                      :: pole
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_aux_fit

      CALL timeset(routineN, handle)

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb
      nmo = admm_env%nmo(ispin)

      nspins = SIZE(admm_env%P_to_be_purified)

      CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff_aux_fit)

      !! * For the time beeing, get the P to be purified from the mo_coeffs
      !! * This needs to be replaced with the a block modified P

      IF (.NOT. blocked) THEN
         CALL parallel_gemm('N', 'T', nao_aux_fit, nao_aux_fit, nmo, &
                            1.0_dp, mo_coeff_aux_fit, mo_coeff_aux_fit, 0.0_dp, &
                            admm_env%P_to_be_purified(ispin))
      END IF

      CALL cp_fm_to_fm(admm_env%S, admm_env%work_aux_aux)
      CALL cp_fm_to_fm(admm_env%P_to_be_purified(ispin), admm_env%work_aux_aux2)

      CALL cp_fm_cholesky_decompose(admm_env%work_aux_aux)

      CALL cp_fm_cholesky_reduce(admm_env%work_aux_aux2, admm_env%work_aux_aux, itype=3)

      CALL cp_fm_syevd(admm_env%work_aux_aux2, admm_env%R_purify(ispin), &
                       admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data)

      CALL cp_fm_cholesky_restore(admm_env%R_purify(ispin), nao_aux_fit, admm_env%work_aux_aux, &
                                  admm_env%work_aux_aux3, op="MULTIPLY", pos="LEFT", transa="T")

      CALL cp_fm_to_fm(admm_env%work_aux_aux3, admm_env%R_purify(ispin))

      ! *** Construct Matrix M for Hadamard Product
      CALL cp_fm_set_all(admm_env%M_purify(ispin), 0.0_dp)
      pole = 0.0_dp
      DO i = 1, nao_aux_fit
         pole = Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i) - 0.5_dp)
         CALL cp_fm_set_element(admm_env%M_purify(ispin), i, i, pole)
      END DO
      CALL cp_fm_upper_to_full(admm_env%M_purify(ispin), admm_env%work_aux_aux)

      CALL copy_dbcsr_to_fm(density_matrix, admm_env%work_aux_aux3)
      CALL cp_fm_upper_to_full(admm_env%work_aux_aux3, admm_env%work_aux_aux)

      ! ** S^(-1)*R
      CALL parallel_gemm('N', 'N', nao_aux_fit, nao_aux_fit, nao_aux_fit, &
                         1.0_dp, admm_env%S_inv, admm_env%R_purify(ispin), 0.0_dp, &
                         admm_env%work_aux_aux)
      ! ** S^(-1)*R*M
      CALL parallel_gemm('N', 'N', nao_aux_fit, nao_aux_fit, nao_aux_fit, &
                         1.0_dp, admm_env%work_aux_aux, admm_env%M_purify(ispin), 0.0_dp, &
                         admm_env%work_aux_aux2)
      ! ** S^(-1)*R*M*R^T*S^(-1)
      CALL parallel_gemm('N', 'T', nao_aux_fit, nao_aux_fit, nao_aux_fit, &
                         1.0_dp, admm_env%work_aux_aux2, admm_env%work_aux_aux, 0.0_dp, &
                         admm_env%work_aux_aux3)

      CALL copy_fm_to_dbcsr(admm_env%work_aux_aux3, density_matrix, keep_sparsity=.TRUE.)

      IF (nspins == 1) THEN
         CALL dbcsr_scale(density_matrix, 2.0_dp)
      END IF

      CALL timestop(handle)

   END SUBROUTINE purify_dm_cauchy

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE merge_ks_matrix_cauchy(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_ks_matrix_cauchy'

      INTEGER                                            :: blk, handle, i, iatom, ispin, j, jatom, &
                                                            nao_aux_fit, nao_orb, nmo
      REAL(dp)                                           :: eig_diff, pole, tmp
      REAL(dp), DIMENSION(:, :), POINTER                 :: sparse_block
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_aux_fit
      TYPE(dbcsr_type), POINTER                          :: matrix_k_tilde
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)
      NULLIFY (admm_env, dft_control, matrix_ks, matrix_ks_aux_fit, mos, mo_coeff)

      CALL get_qs_env(qs_env, &
                      admm_env=admm_env, &
                      dft_control=dft_control, &
                      matrix_ks=matrix_ks, &
                      mos=mos)
      CALL get_admm_env(admm_env, matrix_ks_aux_fit=matrix_ks_aux_fit)

      DO ispin = 1, dft_control%nspins
         nao_aux_fit = admm_env%nao_aux_fit
         nao_orb = admm_env%nao_orb
         nmo = admm_env%nmo(ispin)
         CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff)

         IF (.NOT. admm_env%block_dm) THEN
            !** Get P from mo_coeffs, otherwise we have troubles with occupation numbers ...
            CALL parallel_gemm('N', 'T', nao_orb, nao_orb, nmo, &
                               1.0_dp, mo_coeff, mo_coeff, 0.0_dp, &
                               admm_env%work_orb_orb)

            !! A*P
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb, &
                               1.0_dp, admm_env%A, admm_env%work_orb_orb, 0.0_dp, &
                               admm_env%work_aux_orb2)
            !! A*P*A^T
            CALL parallel_gemm('N', 'T', nao_aux_fit, nao_aux_fit, nao_orb, &
                               1.0_dp, admm_env%work_aux_orb2, admm_env%A, 0.0_dp, &
                               admm_env%P_to_be_purified(ispin))

         END IF

         CALL cp_fm_to_fm(admm_env%S, admm_env%work_aux_aux)
         CALL cp_fm_to_fm(admm_env%P_to_be_purified(ispin), admm_env%work_aux_aux2)

         CALL cp_fm_cholesky_decompose(admm_env%work_aux_aux)

         CALL cp_fm_cholesky_reduce(admm_env%work_aux_aux2, admm_env%work_aux_aux, itype=3)

         CALL cp_fm_syevd(admm_env%work_aux_aux2, admm_env%R_purify(ispin), &
                          admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data)

         CALL cp_fm_cholesky_restore(admm_env%R_purify(ispin), nao_aux_fit, admm_env%work_aux_aux, &
                                     admm_env%work_aux_aux3, op="MULTIPLY", pos="LEFT", transa="T")

         CALL cp_fm_to_fm(admm_env%work_aux_aux3, admm_env%R_purify(ispin))

         ! *** Construct Matrix M for Hadamard Product
         pole = 0.0_dp
         DO i = 1, nao_aux_fit
            DO j = i, nao_aux_fit
               eig_diff = (admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i) - &
                           admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j))
               ! *** two eigenvalues could be the degenerated. In that case use 2nd order formula for the poles
               IF (ABS(eig_diff) == 0.0_dp) THEN
                  pole = delta(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i) - 0.5_dp)
                  CALL cp_fm_set_element(admm_env%M_purify(ispin), i, j, pole)
               ELSE
                  pole = 1.0_dp/(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i) - &
                                 admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j))
                  tmp = Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i) - 0.5_dp)
                  tmp = tmp - Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j) - 0.5_dp)
                  pole = tmp*pole
                  CALL cp_fm_set_element(admm_env%M_purify(ispin), i, j, pole)
               END IF
            END DO
         END DO
         CALL cp_fm_upper_to_full(admm_env%M_purify(ispin), admm_env%work_aux_aux)

         CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix, admm_env%K(ispin))
         CALL cp_fm_upper_to_full(admm_env%K(ispin), admm_env%work_aux_aux)

         !! S^(-1)*R
         CALL parallel_gemm('N', 'N', nao_aux_fit, nao_aux_fit, nao_aux_fit, &
                            1.0_dp, admm_env%S_inv, admm_env%R_purify(ispin), 0.0_dp, &
                            admm_env%work_aux_aux)
         !! K*S^(-1)*R
         CALL parallel_gemm('N', 'N', nao_aux_fit, nao_aux_fit, nao_aux_fit, &
                            1.0_dp, admm_env%K(ispin), admm_env%work_aux_aux, 0.0_dp, &
                            admm_env%work_aux_aux2)
         !! R^T*S^(-1)*K*S^(-1)*R
         CALL parallel_gemm('T', 'N', nao_aux_fit, nao_aux_fit, nao_aux_fit, &
                            1.0_dp, admm_env%work_aux_aux, admm_env%work_aux_aux2, 0.0_dp, &
                            admm_env%work_aux_aux3)
         !! R^T*S^(-1)*K*S^(-1)*R x M
         CALL cp_fm_schur_product(admm_env%work_aux_aux3, admm_env%M_purify(ispin), &
                                  admm_env%work_aux_aux)

         !! R^T*A
         CALL parallel_gemm('T', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                            1.0_dp, admm_env%R_purify(ispin), admm_env%A, 0.0_dp, &
                            admm_env%work_aux_orb)

         !! (R^T*S^(-1)*K*S^(-1)*R x M) * R^T*A
         CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                            1.0_dp, admm_env%work_aux_aux, admm_env%work_aux_orb, 0.0_dp, &
                            admm_env%work_aux_orb2)
         !! A^T*R*(R^T*S^(-1)*K*S^(-1)*R x M) * R^T*A
         CALL parallel_gemm('T', 'N', nao_orb, nao_orb, nao_aux_fit, &
                            1.0_dp, admm_env%work_aux_orb, admm_env%work_aux_orb2, 0.0_dp, &
                            admm_env%work_orb_orb)

         NULLIFY (matrix_k_tilde)
         ALLOCATE (matrix_k_tilde)
         CALL dbcsr_create(matrix_k_tilde, template=matrix_ks(ispin)%matrix, &
                           name='MATRIX K_tilde', &
                           matrix_type=dbcsr_type_symmetric)

         CALL cp_fm_to_fm(admm_env%work_orb_orb, admm_env%ks_to_be_merged(ispin))

         CALL dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix)
         CALL dbcsr_set(matrix_k_tilde, 0.0_dp)
         CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE.)

         IF (admm_env%block_dm) THEN
            ! ** now loop through the list and nullify blocks
            CALL dbcsr_iterator_start(iter, matrix_k_tilde)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
               IF (admm_env%block_map(iatom, jatom) == 0) THEN
                  sparse_block = 0.0_dp
               END IF
            END DO
            CALL dbcsr_iterator_stop(iter)
         END IF

         CALL dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp)

         CALL dbcsr_deallocate_matrix(matrix_k_tilde)

      END DO !spin-loop

      CALL timestop(handle)

   END SUBROUTINE merge_ks_matrix_cauchy

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE merge_ks_matrix_cauchy_subspace(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_ks_matrix_cauchy_subspace'

      INTEGER                                            :: handle, ispin, nao_aux_fit, nao_orb, nmo
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_aux_fit
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_aux_fit
      TYPE(dbcsr_type), POINTER                          :: matrix_k_tilde
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos, mos_aux_fit

      CALL timeset(routineN, handle)
      NULLIFY (admm_env, dft_control, matrix_ks, matrix_ks_aux_fit, mos, mos_aux_fit, &
               mo_coeff, mo_coeff_aux_fit)

      CALL get_qs_env(qs_env, &
                      admm_env=admm_env, &
                      dft_control=dft_control, &
                      matrix_ks=matrix_ks, &
                      mos=mos)
      CALL get_admm_env(admm_env, matrix_ks_aux_fit=matrix_ks_aux_fit, mos_aux_fit=mos_aux_fit)

      DO ispin = 1, dft_control%nspins
         nao_aux_fit = admm_env%nao_aux_fit
         nao_orb = admm_env%nao_orb
         nmo = admm_env%nmo(ispin)
         CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff)
         CALL get_mo_set(mo_set=mos_aux_fit(ispin), mo_coeff=mo_coeff_aux_fit)

         !! Calculate Lambda^{-2}
         CALL cp_fm_to_fm(admm_env%lambda(ispin), admm_env%work_nmo_nmo1(ispin))
         CALL cp_fm_cholesky_decompose(admm_env%work_nmo_nmo1(ispin))
         CALL cp_fm_cholesky_invert(admm_env%work_nmo_nmo1(ispin))
         !! Symmetrize the guy
         CALL cp_fm_upper_to_full(admm_env%work_nmo_nmo1(ispin), admm_env%lambda_inv2(ispin))
         !! Take square
         CALL parallel_gemm('N', 'T', nmo, nmo, nmo, &
                            1.0_dp, admm_env%work_nmo_nmo1(ispin), admm_env%work_nmo_nmo1(ispin), 0.0_dp, &
                            admm_env%lambda_inv2(ispin))

         !! ** C_hat = AC
         CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nao_orb, &
                            1.0_dp, admm_env%A, mo_coeff, 0.0_dp, &
                            admm_env%C_hat(ispin))

         !! calc P_tilde from C_hat
         CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nmo, &
                            1.0_dp, admm_env%C_hat(ispin), admm_env%lambda_inv(ispin), 0.0_dp, &
                            admm_env%work_aux_nmo(ispin))

         CALL parallel_gemm('N', 'T', nao_aux_fit, nao_aux_fit, nmo, &
                            1.0_dp, admm_env%C_hat(ispin), admm_env%work_aux_nmo(ispin), 0.0_dp, &
                            admm_env%P_tilde(ispin))

         !! ** C_hat*Lambda^{-2}
         CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nmo, &
                            1.0_dp, admm_env%C_hat(ispin), admm_env%lambda_inv2(ispin), 0.0_dp, &
                            admm_env%work_aux_nmo(ispin))

         !! ** C_hat*Lambda^{-2}*C_hat^T
         CALL parallel_gemm('N', 'T', nao_aux_fit, nao_aux_fit, nmo, &
                            1.0_dp, admm_env%work_aux_nmo(ispin), admm_env%C_hat(ispin), 0.0_dp, &
                            admm_env%work_aux_aux)

         !! ** S*C_hat*Lambda^{-2}*C_hat^T
         CALL parallel_gemm('N', 'N', nao_aux_fit, nao_aux_fit, nao_aux_fit, &
                            1.0_dp, admm_env%S, admm_env%work_aux_aux, 0.0_dp, &
                            admm_env%work_aux_aux2)

         CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix, admm_env%K(ispin))
         CALL cp_fm_upper_to_full(admm_env%K(ispin), admm_env%work_aux_aux)

         !! ** S*C_hat*Lambda^{-2}*C_hat^T*H_tilde
         CALL parallel_gemm('N', 'N', nao_aux_fit, nao_aux_fit, nao_aux_fit, &
                            1.0_dp, admm_env%work_aux_aux2, admm_env%K(ispin), 0.0_dp, &
                            admm_env%work_aux_aux)

         !! ** P_tilde*S
         CALL parallel_gemm('N', 'N', nao_aux_fit, nao_aux_fit, nao_aux_fit, &
                            1.0_dp, admm_env%P_tilde(ispin), admm_env%S, 0.0_dp, &
                            admm_env%work_aux_aux2)

         !! ** -S*C_hat*Lambda^{-2}*C_hat^T*H_tilde*P_tilde*S
         CALL parallel_gemm('N', 'N', nao_aux_fit, nao_aux_fit, nao_aux_fit, &
                            -1.0_dp, admm_env%work_aux_aux, admm_env%work_aux_aux2, 0.0_dp, &
                            admm_env%work_aux_aux3)

         !! ** -S*C_hat*Lambda^{-2}*C_hat^T*H_tilde*P_tilde*S+S*C_hat*Lambda^{-2}*C_hat^T*H_tilde
         CALL cp_fm_scale_and_add(1.0_dp, admm_env%work_aux_aux3, 1.0_dp, admm_env%work_aux_aux)

         !! first_part*A
         CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                            1.0_dp, admm_env%work_aux_aux3, admm_env%A, 0.0_dp, &
                            admm_env%work_aux_orb)

         !! + first_part^T*A
         CALL parallel_gemm('T', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                            1.0_dp, admm_env%work_aux_aux3, admm_env%A, 1.0_dp, &
                            admm_env%work_aux_orb)

         !! A^T*(first+seccond)=H
         CALL parallel_gemm('T', 'N', nao_orb, nao_orb, nao_aux_fit, &
                            1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                            admm_env%work_orb_orb)

         NULLIFY (matrix_k_tilde)
         ALLOCATE (matrix_k_tilde)
         CALL dbcsr_create(matrix_k_tilde, template=matrix_ks(ispin)%matrix, &
                           name='MATRIX K_tilde', &
                           matrix_type=dbcsr_type_symmetric)

         CALL cp_fm_to_fm(admm_env%work_orb_orb, admm_env%ks_to_be_merged(ispin))

         CALL dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix)
         CALL dbcsr_set(matrix_k_tilde, 0.0_dp)
         CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE.)

         CALL parallel_gemm('N', 'N', nao_orb, nmo, nao_orb, &
                            1.0_dp, admm_env%work_orb_orb, mo_coeff, 0.0_dp, &
                            admm_env%mo_derivs_tmp(ispin))

         CALL dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp)

         CALL dbcsr_deallocate_matrix(matrix_k_tilde)

      END DO !spin loop
      CALL timestop(handle)

   END SUBROUTINE merge_ks_matrix_cauchy_subspace

! **************************************************************************************************
!> \brief Calculates the product Kohn-Sham-Matrix x mo_coeff for the auxiliary
!>        basis set and transforms it into the orbital basis. This is needed
!>        in order to use OT
!>
!> \param ispin which spin to transform
!> \param admm_env The ADMM env
!> \param mo_set ...
!> \param mo_coeff the MO coefficients from the orbital basis set
!> \param mo_coeff_aux_fit the MO coefficients from the auxiliary fitting basis set
!> \param mo_derivs KS x mo_coeff from the orbital basis set to which we add the
!>        auxiliary basis set part
!> \param mo_derivs_aux_fit ...
!> \param matrix_ks_aux_fit the Kohn-Sham matrix from the auxiliary fitting basis set
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE merge_mo_derivs_diag(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, mo_derivs, &
                                   mo_derivs_aux_fit, matrix_ks_aux_fit)
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(mo_set_type), INTENT(IN)                      :: mo_set
      TYPE(cp_fm_type), INTENT(IN)                       :: mo_coeff, mo_coeff_aux_fit
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_derivs, mo_derivs_aux_fit
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks_aux_fit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_mo_derivs_diag'

      INTEGER                                            :: handle, i, j, nao_aux_fit, nao_orb, nmo
      REAL(dp)                                           :: eig_diff, pole, tmp32, tmp52, tmp72, &
                                                            tmp92
      REAL(dp), DIMENSION(:), POINTER                    :: occupation_numbers, scaling_factor

      CALL timeset(routineN, handle)

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb
      nmo = admm_env%nmo(ispin)

      CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix, admm_env%K(ispin))
      CALL cp_fm_upper_to_full(admm_env%K(ispin), admm_env%work_aux_aux)

      CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nao_aux_fit, &
                         1.0_dp, admm_env%K(ispin), mo_coeff_aux_fit, 0.0_dp, &
                         admm_env%H(ispin))

      CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
      ALLOCATE (scaling_factor(SIZE(occupation_numbers)))
      scaling_factor = 2.0_dp*occupation_numbers

      CALL cp_fm_column_scale(admm_env%H(ispin), scaling_factor)

      CALL cp_fm_to_fm(admm_env%H(ispin), mo_derivs_aux_fit(ispin))

      ! *** Add first term
      CALL parallel_gemm('N', 'T', nao_aux_fit, nmo, nmo, &
                         1.0_dp, admm_env%H(ispin), admm_env%lambda_inv_sqrt(ispin), 0.0_dp, &
                         admm_env%work_aux_nmo(ispin))
      CALL parallel_gemm('T', 'N', nao_orb, nmo, nao_aux_fit, &
                         1.0_dp, admm_env%A, admm_env%work_aux_nmo(ispin), 0.0_dp, &
                         admm_env%mo_derivs_tmp(ispin))

      ! *** Construct Matrix M for Hadamard Product
      pole = 0.0_dp
      DO i = 1, nmo
         DO j = i, nmo
            eig_diff = (admm_env%eigvals_lambda(ispin)%eigvals%data(i) - &
                        admm_env%eigvals_lambda(ispin)%eigvals%data(j))
            ! *** two eigenvalues could be the degenerated. In that case use 2nd order formula for the poles
            IF (ABS(eig_diff) < 0.0001_dp) THEN
               tmp32 = 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(j))**3
               tmp52 = tmp32/admm_env%eigvals_lambda(ispin)%eigvals%data(j)*eig_diff
               tmp72 = tmp52/admm_env%eigvals_lambda(ispin)%eigvals%data(j)*eig_diff
               tmp92 = tmp72/admm_env%eigvals_lambda(ispin)%eigvals%data(j)*eig_diff

               pole = -0.5_dp*tmp32 + 3.0_dp/8.0_dp*tmp52 - 5.0_dp/16.0_dp*tmp72 + 35.0_dp/128.0_dp*tmp92
               CALL cp_fm_set_element(admm_env%M(ispin), i, j, pole)
            ELSE
               pole = 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(i))
               pole = pole - 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(j))
               pole = pole/(admm_env%eigvals_lambda(ispin)%eigvals%data(i) - &
                            admm_env%eigvals_lambda(ispin)%eigvals%data(j))
               CALL cp_fm_set_element(admm_env%M(ispin), i, j, pole)
            END IF
         END DO
      END DO
      CALL cp_fm_upper_to_full(admm_env%M(ispin), admm_env%work_nmo_nmo1(ispin))

      ! *** 2nd term to be added to fm_H

      !! Part 1: B^(T)*C* R*[R^(T)*c^(T)*A^(T)*H_aux_fit*R x M]*R^(T)
      !! Part 2: B*C*(R*[R^(T)*c^(T)*A^(T)*H_aux_fit*R x M]*R^(T))^(T)

      ! *** H'*R
      CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nmo, &
                         1.0_dp, admm_env%H(ispin), admm_env%R(ispin), 0.0_dp, &
                         admm_env%work_aux_nmo(ispin))
      ! *** A^(T)*H'*R
      CALL parallel_gemm('T', 'N', nao_orb, nmo, nao_aux_fit, &
                         1.0_dp, admm_env%A, admm_env%work_aux_nmo(ispin), 0.0_dp, &
                         admm_env%work_orb_nmo(ispin))
      ! *** c^(T)*A^(T)*H'*R
      CALL parallel_gemm('T', 'N', nmo, nmo, nao_orb, &
                         1.0_dp, mo_coeff, admm_env%work_orb_nmo(ispin), 0.0_dp, &
                         admm_env%work_nmo_nmo1(ispin))
      ! *** R^(T)*c^(T)*A^(T)*H'*R
      CALL parallel_gemm('T', 'N', nmo, nmo, nmo, &
                         1.0_dp, admm_env%R(ispin), admm_env%work_nmo_nmo1(ispin), 0.0_dp, &
                         admm_env%work_nmo_nmo2(ispin))
      ! *** R^(T)*c^(T)*A^(T)*H'*R x M
      CALL cp_fm_schur_product(admm_env%work_nmo_nmo2(ispin), &
                               admm_env%M(ispin), admm_env%work_nmo_nmo1(ispin))
      ! *** R* (R^(T)*c^(T)*A^(T)*H'*R x M)
      CALL parallel_gemm('N', 'N', nmo, nmo, nmo, &
                         1.0_dp, admm_env%R(ispin), admm_env%work_nmo_nmo1(ispin), 0.0_dp, &
                         admm_env%work_nmo_nmo2(ispin))

      ! *** R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)
      CALL parallel_gemm('N', 'T', nmo, nmo, nmo, &
                         1.0_dp, admm_env%work_nmo_nmo2(ispin), admm_env%R(ispin), 0.0_dp, &
                         admm_env%R_schur_R_t(ispin))

      ! *** B^(T)*c
      CALL parallel_gemm('T', 'N', nao_orb, nmo, nao_orb, &
                         1.0_dp, admm_env%B, mo_coeff, 0.0_dp, &
                         admm_env%work_orb_nmo(ispin))

      ! *** Add first term to fm_H
      ! *** B^(T)*c* R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)
      CALL parallel_gemm('N', 'N', nao_orb, nmo, nmo, &
                         1.0_dp, admm_env%work_orb_nmo(ispin), admm_env%R_schur_R_t(ispin), 1.0_dp, &
                         admm_env%mo_derivs_tmp(ispin))

      ! *** Add second term to fm_H
      ! *** B*C *[ R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)]^(T)
      CALL parallel_gemm('N', 'T', nao_orb, nmo, nmo, &
                         1.0_dp, admm_env%work_orb_nmo(ispin), admm_env%R_schur_R_t(ispin), 1.0_dp, &
                         admm_env%mo_derivs_tmp(ispin))

      DO i = 1, SIZE(scaling_factor)
         scaling_factor(i) = 1.0_dp/scaling_factor(i)
      END DO

      CALL cp_fm_column_scale(admm_env%mo_derivs_tmp(ispin), scaling_factor)

      CALL cp_fm_scale_and_add(1.0_dp, mo_derivs(ispin), 1.0_dp, admm_env%mo_derivs_tmp(ispin))

      DEALLOCATE (scaling_factor)

      CALL timestop(handle)

   END SUBROUTINE merge_mo_derivs_diag

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE merge_ks_matrix_none(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_ks_matrix_none'

      INTEGER                                            :: blk, handle, iatom, ispin, jatom, &
                                                            nao_aux_fit, nao_orb, nmo
      REAL(dp), DIMENSION(:, :), POINTER                 :: sparse_block
      REAL(KIND=dp)                                      :: ener_k(2), ener_x(2), ener_x1(2), &
                                                            gsi_square, trace_tmp, trace_tmp_two
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: matrix_ks, matrix_ks_aux_fit, &
         matrix_ks_aux_fit_dft, matrix_ks_aux_fit_hfx, matrix_s, matrix_s_aux_fit, rho_ao, &
         rho_ao_aux
      TYPE(dbcsr_type), POINTER                          :: matrix_k_tilde, &
                                                            matrix_ks_aux_fit_admms_tmp, &
                                                            matrix_TtsT
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit

      CALL timeset(routineN, handle)
      NULLIFY (admm_env, dft_control, matrix_ks, matrix_ks_aux_fit, matrix_ks_aux_fit_dft, &
               matrix_ks_aux_fit_hfx, matrix_s, matrix_s_aux_fit, rho_ao, rho_ao_aux, matrix_k_tilde, &
               matrix_TtsT, matrix_ks_aux_fit_admms_tmp, rho, rho_aux_fit, sparse_block, para_env, energy)

      CALL get_qs_env(qs_env, &
                      admm_env=admm_env, &
                      dft_control=dft_control, &
                      matrix_ks=matrix_ks, &
                      rho=rho, &
                      matrix_s=matrix_s, &
                      energy=energy, &
                      para_env=para_env)
      CALL get_admm_env(admm_env, matrix_ks_aux_fit=matrix_ks_aux_fit, matrix_ks_aux_fit_dft=matrix_ks_aux_fit_dft, &
                        matrix_ks_aux_fit_hfx=matrix_ks_aux_fit_hfx, rho_aux_fit=rho_aux_fit, &
                        matrix_s_aux_fit=matrix_s_aux_fit)

      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CALL qs_rho_get(rho_aux_fit, &
                      rho_ao=rho_ao_aux)

      DO ispin = 1, dft_control%nspins
         IF (admm_env%block_dm) THEN
            CALL dbcsr_iterator_start(iter, matrix_ks_aux_fit(ispin)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
               IF (admm_env%block_map(iatom, jatom) == 0) THEN
                  sparse_block = 0.0_dp
               END IF
            END DO
            CALL dbcsr_iterator_stop(iter)
            CALL dbcsr_add(matrix_ks(ispin)%matrix, matrix_ks_aux_fit(ispin)%matrix, 1.0_dp, 1.0_dp)

         ELSE

            nao_aux_fit = admm_env%nao_aux_fit
            nao_orb = admm_env%nao_orb
            nmo = admm_env%nmo(ispin)

            ! ADMMS: different matrix for calculating A^(T)*K*A, see Eq. (37) Merlot
            IF (admm_env%do_admms) THEN
               NULLIFY (matrix_ks_aux_fit_admms_tmp)
               ALLOCATE (matrix_ks_aux_fit_admms_tmp)
               CALL dbcsr_create(matrix_ks_aux_fit_admms_tmp, template=matrix_ks_aux_fit(ispin)%matrix, &
                                 name='matrix_ks_aux_fit_admms_tmp', matrix_type='s')
               ! matrix_ks_aux_fit_admms_tmp = k(d_Q)
               CALL dbcsr_copy(matrix_ks_aux_fit_admms_tmp, matrix_ks_aux_fit_hfx(ispin)%matrix)

               ! matrix_ks_aux_fit_admms_tmp = k(d_Q) - gsi^2/3 x(d_Q)
               CALL dbcsr_add(matrix_ks_aux_fit_admms_tmp, matrix_ks_aux_fit_dft(ispin)%matrix, &
                              1.0_dp, -(admm_env%gsi(ispin))**(2.0_dp/3.0_dp))
               CALL copy_dbcsr_to_fm(matrix_ks_aux_fit_admms_tmp, admm_env%K(ispin))
               CALL dbcsr_deallocate_matrix(matrix_ks_aux_fit_admms_tmp)
            ELSE
               CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix, admm_env%K(ispin))
            END IF

            CALL cp_fm_upper_to_full(admm_env%K(ispin), admm_env%work_aux_aux)

            !! K*A
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                               1.0_dp, admm_env%K(ispin), admm_env%A, 0.0_dp, &
                               admm_env%work_aux_orb)
            !! A^T*K*A
            CALL parallel_gemm('T', 'N', nao_orb, nao_orb, nao_aux_fit, &
                               1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                               admm_env%work_orb_orb)

            NULLIFY (matrix_k_tilde)
            ALLOCATE (matrix_k_tilde)
            CALL dbcsr_create(matrix_k_tilde, template=matrix_ks(ispin)%matrix, &
                              name='MATRIX K_tilde', matrix_type='S')
            CALL dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix)
            CALL dbcsr_set(matrix_k_tilde, 0.0_dp)
            CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE.)

            ! Scale matrix_K_tilde here. Then, the scaling has to be done for forces separately
            ! Scale matrix_K_tilde by gsi for ADMMQ and ADMMS (Eqs. (27), (37) in Merlot, 2014)
            IF (admm_env%do_admmq .OR. admm_env%do_admms) THEN
               CALL dbcsr_scale(matrix_k_tilde, admm_env%gsi(ispin))
            END IF

            ! Scale matrix_K_tilde by gsi^2 for ADMMP (Eq. (35) in Merlot, 2014)
            IF (admm_env%do_admmp) THEN
               gsi_square = (admm_env%gsi(ispin))*(admm_env%gsi(ispin))
               CALL dbcsr_scale(matrix_k_tilde, gsi_square)
            END IF

            admm_env%lambda_merlot(ispin) = 0

            ! Calculate LAMBDA according to Merlot, 1. IF: ADMMQ, 2. IF: ADMMP, 3. IF: ADMMS,
            IF (admm_env%do_admmq) THEN
               CALL dbcsr_dot(matrix_ks_aux_fit(ispin)%matrix, rho_ao_aux(ispin)%matrix, trace_tmp)

               ! Factor of 2 is missing compared to Eq. 28 in Merlot due to
               ! Tr(ds) = N in the code \neq 2N in Merlot
               admm_env%lambda_merlot(ispin) = trace_tmp/(admm_env%n_large_basis(ispin))

            ELSE IF (admm_env%do_admmp) THEN
               IF (dft_control%nspins == 2) THEN
                  CALL calc_spin_dep_aux_exch_ener(qs_env=qs_env, admm_env=admm_env, ener_k_ispin=ener_k(ispin), &
                                                   ener_x_ispin=ener_x(ispin), ener_x1_ispin=ener_x1(ispin), &
                                                   ispin=ispin)
                  admm_env%lambda_merlot(ispin) = 2.0_dp*(admm_env%gsi(ispin))**2* &
                                                  (ener_k(ispin) + ener_x(ispin) + ener_x1(ispin))/ &
                                                  (admm_env%n_large_basis(ispin))

               ELSE
                  admm_env%lambda_merlot(ispin) = 2.0_dp*(admm_env%gsi(ispin))**2* &
                                                  (energy%ex + energy%exc_aux_fit + energy%exc1_aux_fit) &
                                                  /(admm_env%n_large_basis(ispin))
               END IF

            ELSE IF (admm_env%do_admms) THEN
               CALL dbcsr_dot(matrix_ks_aux_fit_hfx(ispin)%matrix, rho_ao_aux(ispin)%matrix, trace_tmp)
               CALL dbcsr_dot(matrix_ks_aux_fit_dft(ispin)%matrix, rho_ao_aux(ispin)%matrix, trace_tmp_two)
               ! For ADMMS open-shell case we need k and x (Merlot) separately since gsi(a)\=gsi(b)
               IF (dft_control%nspins == 2) THEN
                  CALL calc_spin_dep_aux_exch_ener(qs_env=qs_env, admm_env=admm_env, ener_k_ispin=ener_k(ispin), &
                                                   ener_x_ispin=ener_x(ispin), ener_x1_ispin=ener_x1(ispin), &
                                                   ispin=ispin)
                  admm_env%lambda_merlot(ispin) = &
                     (trace_tmp + 2.0_dp/3.0_dp*((admm_env%gsi(ispin))**(2.0_dp/3.0_dp))* &
                      (ener_x(ispin) + ener_x1(ispin)) - ((admm_env%gsi(ispin))**(2.0_dp/3.0_dp))* &
                      trace_tmp_two)/(admm_env%n_large_basis(ispin))

               ELSE
                  admm_env%lambda_merlot(ispin) = (trace_tmp + (admm_env%gsi(ispin))**(2.0_dp/3.0_dp)* &
                                                   (2.0_dp/3.0_dp*(energy%exc_aux_fit + energy%exc1_aux_fit) - &
                                                    trace_tmp_two))/(admm_env%n_large_basis(ispin))
               END IF
            END IF

            ! Calculate variational distribution to KS matrix according
            ! to Eqs. (27), (35) and (37) in Merlot, 2014

            IF (admm_env%do_admmp .OR. admm_env%do_admmq .OR. admm_env%do_admms) THEN

               !! T^T*s_aux*T in (27) Merlot (T=A), as calculating A^T*K*A few lines above
               CALL copy_dbcsr_to_fm(matrix_s_aux_fit(1)%matrix, admm_env%work_aux_aux4)
               CALL cp_fm_upper_to_full(admm_env%work_aux_aux4, admm_env%work_aux_aux5)

               ! s_aux*T
               CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                                  1.0_dp, admm_env%work_aux_aux4, admm_env%A, 0.0_dp, &
                                  admm_env%work_aux_orb3)
               ! T^T*s_aux*T
               CALL parallel_gemm('T', 'N', nao_orb, nao_orb, nao_aux_fit, &
                                  1.0_dp, admm_env%A, admm_env%work_aux_orb3, 0.0_dp, &
                                  admm_env%work_orb_orb3)

               NULLIFY (matrix_TtsT)
               ALLOCATE (matrix_TtsT)
               CALL dbcsr_create(matrix_TtsT, template=matrix_ks(ispin)%matrix, &
                                 name='MATRIX TtsT', matrix_type='S')
               CALL dbcsr_copy(matrix_TtsT, matrix_ks(ispin)%matrix)
               CALL dbcsr_set(matrix_TtsT, 0.0_dp)
               CALL copy_fm_to_dbcsr(admm_env%work_orb_orb3, matrix_TtsT, keep_sparsity=.TRUE.)

               !Add -(gsi)*Lambda*TtsT and Lambda*S to the KS matrix according to Merlot2014

               CALL dbcsr_add(matrix_ks(ispin)%matrix, matrix_TtsT, 1.0_dp, &
                              (-admm_env%lambda_merlot(ispin))*admm_env%gsi(ispin))

               CALL dbcsr_add(matrix_ks(ispin)%matrix, matrix_s(1)%matrix, 1.0_dp, admm_env%lambda_merlot(ispin))

               CALL dbcsr_deallocate_matrix(matrix_TtsT)

            END IF

            CALL dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp)

            CALL dbcsr_deallocate_matrix(matrix_k_tilde)

         END IF
      END DO !spin loop

      ! Scale energy for ADMMP and ADMMS
      IF (admm_env%do_admmp) THEN
         !       ener_k = ener_k*(admm_env%gsi(1))*(admm_env%gsi(1))
         !       ener_x = ener_x*(admm_env%gsi(1))*(admm_env%gsi(1))
         !        PRINT *, 'energy%ex = ', energy%ex
         IF (dft_control%nspins == 2) THEN
            energy%exc_aux_fit = 0.0_dp
            energy%exc1_aux_fit = 0.0_dp
            energy%ex = 0.0_dp
            DO ispin = 1, dft_control%nspins
               energy%exc_aux_fit = energy%exc_aux_fit + (admm_env%gsi(ispin))**2.0_dp*ener_x(ispin)
               energy%exc1_aux_fit = energy%exc1_aux_fit + (admm_env%gsi(ispin))**2.0_dp*ener_x1(ispin)
               energy%ex = energy%ex + (admm_env%gsi(ispin))**2.0_dp*ener_k(ispin)
            END DO
         ELSE
            energy%exc_aux_fit = (admm_env%gsi(1))**2.0_dp*energy%exc_aux_fit
            energy%exc1_aux_fit = (admm_env%gsi(1))**2.0_dp*energy%exc1_aux_fit
            energy%ex = (admm_env%gsi(1))**2.0_dp*energy%ex
         END IF

      ELSE IF (admm_env%do_admms) THEN
         IF (dft_control%nspins == 2) THEN
            energy%exc_aux_fit = 0.0_dp
            energy%exc1_aux_fit = 0.0_dp
            DO ispin = 1, dft_control%nspins
               energy%exc_aux_fit = energy%exc_aux_fit + (admm_env%gsi(ispin))**(2.0_dp/3.0_dp)*ener_x(ispin)
               energy%exc1_aux_fit = energy%exc1_aux_fit + (admm_env%gsi(ispin))**(2.0_dp/3.0_dp)*ener_x1(ispin)
            END DO
         ELSE
            energy%exc_aux_fit = (admm_env%gsi(1))**(2.0_dp/3.0_dp)*energy%exc_aux_fit
            energy%exc1_aux_fit = (admm_env%gsi(1))**(2.0_dp/3.0_dp)*energy%exc1_aux_fit
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE merge_ks_matrix_none

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE merge_ks_matrix_none_kp(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_ks_matrix_none_kp'

      COMPLEX(dp)                                        :: fac, fac2
      INTEGER                                            :: handle, i, igroup, ik, ikp, img, indx, &
                                                            ispin, kplocal, nao_aux_fit, nao_orb, &
                                                            natom, nkp, nkp_groups, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: my_kpgrp, use_real_wfn
      REAL(dp)                                           :: ener_k(2), ener_x(2), ener_x1(2), tmp, &
                                                            trace_tmp, trace_tmp_two
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(copy_info_type), ALLOCATABLE, DIMENSION(:, :) :: info
      TYPE(cp_cfm_type)                                  :: cA, cK, cS, cwork_aux_aux, &
                                                            cwork_aux_orb, cwork_orb_orb
      TYPE(cp_fm_struct_type), POINTER                   :: struct_aux_aux, struct_aux_orb, &
                                                            struct_orb_orb
      TYPE(cp_fm_type)                                   :: fmdummy, work_aux_aux, work_aux_aux2, &
                                                            work_aux_orb
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fmwork
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_ks
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER :: matrix_k_tilde, matrix_ks_aux_fit, &
         matrix_ks_aux_fit_dft, matrix_ks_aux_fit_hfx, matrix_ks_kp, matrix_s, matrix_s_aux_fit, &
         rho_ao_aux
      TYPE(dbcsr_type)                                   :: tmpmatrix_ks
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: ksmatrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux_fit, sab_kp
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho_aux_fit
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CALL timeset(routineN, handle)
      NULLIFY (admm_env, rho_ao_aux, rho_aux_fit, &
               matrix_s_aux_fit, energy, &
               para_env, kpoints, sab_aux_fit, &
               matrix_k_tilde, matrix_ks_kp, matrix_ks_aux_fit, scf_env, &
               struct_orb_orb, struct_aux_orb, struct_aux_aux, kp, &
               matrix_ks_aux_fit_hfx, matrix_ks_aux_fit_dft)

      CALL get_qs_env(qs_env, &
                      admm_env=admm_env, &
                      dft_control=dft_control, &
                      matrix_ks_kp=matrix_ks_kp, &
                      matrix_s_kp=matrix_s, &
                      para_env=para_env, &
                      scf_env=scf_env, &
                      natom=natom, &
                      kpoints=kpoints, &
                      energy=energy)

      CALL get_admm_env(admm_env, &
                        matrix_ks_aux_fit_kp=matrix_ks_aux_fit, &
                        matrix_ks_aux_fit_hfx_kp=matrix_ks_aux_fit_hfx, &
                        matrix_ks_aux_fit_dft_kp=matrix_ks_aux_fit_dft, &
                        matrix_s_aux_fit_kp=matrix_s_aux_fit, &
                        sab_aux_fit=sab_aux_fit, &
                        rho_aux_fit=rho_aux_fit)
      CALL qs_rho_get(rho_aux_fit, rho_ao_kp=rho_ao_aux)

      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, use_real_wfn=use_real_wfn, kp_range=kp_range, &
                           nkp_groups=nkp_groups, kp_dist=kp_dist, sab_nl=sab_kp, &
                           cell_to_index=cell_to_index)

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb
      nspins = dft_control%nspins

      !Case study on ADMMQ, ADMMS and ADMMP

      !ADMMQ: calculate lamda as in Merlot eq (28)
      IF (admm_env%do_admmq) THEN
         admm_env%lambda_merlot = 0.0_dp
         DO img = 1, dft_control%nimages
            DO ispin = 1, nspins
               CALL dbcsr_dot(matrix_ks_aux_fit(ispin, img)%matrix, rho_ao_aux(ispin, img)%matrix, trace_tmp)
               admm_env%lambda_merlot(ispin) = admm_env%lambda_merlot(ispin) + trace_tmp/admm_env%n_large_basis(ispin)
            END DO
         END DO
      END IF

      !ADMMP: calculate lamda as in Merlot eq (34)
      IF (admm_env%do_admmp) THEN
         IF (nspins == 1) THEN
            admm_env%lambda_merlot(1) = 2.0_dp*(admm_env%gsi(1))**2* &
                                        (energy%ex + energy%exc_aux_fit + energy%exc1_aux_fit) &
                                        /(admm_env%n_large_basis(1))
         ELSE
            DO ispin = 1, nspins
               CALL calc_spin_dep_aux_exch_ener(qs_env=qs_env, admm_env=admm_env, &
                                                ener_k_ispin=ener_k(ispin), ener_x_ispin=ener_x(ispin), &
                                                ener_x1_ispin=ener_x1(ispin), ispin=ispin)
               admm_env%lambda_merlot(ispin) = 2.0_dp*(admm_env%gsi(ispin))**2* &
                                               (ener_k(ispin) + ener_x(ispin) + ener_x1(ispin))/ &
                                               (admm_env%n_large_basis(ispin))
            END DO
         END IF
      END IF

      !ADMMS: calculate lambda as in Merlot eq (36)
      IF (admm_env%do_admms) THEN
         IF (nspins == 1) THEN
            trace_tmp = 0.0_dp
            trace_tmp_two = 0.0_dp
            DO img = 1, dft_control%nimages
               CALL dbcsr_dot(matrix_ks_aux_fit_hfx(1, img)%matrix, rho_ao_aux(1, img)%matrix, tmp)
               trace_tmp = trace_tmp + tmp
               CALL dbcsr_dot(matrix_ks_aux_fit_dft(1, img)%matrix, rho_ao_aux(1, img)%matrix, tmp)
               trace_tmp_two = trace_tmp_two + tmp
            END DO
            admm_env%lambda_merlot(1) = (trace_tmp + (admm_env%gsi(1))**(2.0_dp/3.0_dp)* &
                                         (2.0_dp/3.0_dp*(energy%exc_aux_fit + energy%exc1_aux_fit) - &
                                          trace_tmp_two))/(admm_env%n_large_basis(1))
         ELSE

            DO ispin = 1, nspins
               trace_tmp = 0.0_dp
               trace_tmp_two = 0.0_dp
               DO img = 1, dft_control%nimages
                  CALL dbcsr_dot(matrix_ks_aux_fit_hfx(ispin, img)%matrix, rho_ao_aux(ispin, img)%matrix, tmp)
                  trace_tmp = trace_tmp + tmp
                  CALL dbcsr_dot(matrix_ks_aux_fit_dft(ispin, img)%matrix, rho_ao_aux(ispin, img)%matrix, tmp)
                  trace_tmp_two = trace_tmp_two + tmp
               END DO

               CALL calc_spin_dep_aux_exch_ener(qs_env=qs_env, admm_env=admm_env, &
                                                ener_k_ispin=ener_k(ispin), ener_x_ispin=ener_x(ispin), &
                                                ener_x1_ispin=ener_x1(ispin), ispin=ispin)

               admm_env%lambda_merlot(ispin) = &
                  (trace_tmp + 2.0_dp/3.0_dp*((admm_env%gsi(ispin))**(2.0_dp/3.0_dp))* &
                   (ener_x(ispin) + ener_x1(ispin)) - ((admm_env%gsi(ispin))**(2.0_dp/3.0_dp))* &
                   trace_tmp_two)/(admm_env%n_large_basis(ispin))
            END DO
         END IF

         !Here we buld the KS matrix: KS_hfx = gsi^2/3*KS_dft, the we then pass as the ususal KS_aux_fit
         NULLIFY (matrix_ks_aux_fit)
         ALLOCATE (matrix_ks_aux_fit(nspins, dft_control%nimages))
         DO img = 1, dft_control%nimages
            DO ispin = 1, nspins
               NULLIFY (matrix_ks_aux_fit(ispin, img)%matrix)
               ALLOCATE (matrix_ks_aux_fit(ispin, img)%matrix)
               CALL dbcsr_create(matrix_ks_aux_fit(ispin, img)%matrix, template=matrix_s_aux_fit(1, 1)%matrix)
               CALL dbcsr_copy(matrix_ks_aux_fit(ispin, img)%matrix, matrix_ks_aux_fit_hfx(ispin, img)%matrix)
               CALL dbcsr_add(matrix_ks_aux_fit(ispin, img)%matrix, matrix_ks_aux_fit_dft(ispin, img)%matrix, &
                              1.0_dp, -admm_env%gsi(ispin)**(2.0_dp/3.0_dp))
            END DO
         END DO
      END IF

      ! the temporary DBCSR matrices for the rskp_transform we have to manually allocate
      ALLOCATE (ksmatrix(2))
      CALL dbcsr_create(ksmatrix(1), template=matrix_ks_aux_fit(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(ksmatrix(2), template=matrix_ks_aux_fit(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL dbcsr_create(tmpmatrix_ks, template=matrix_ks_aux_fit(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(ksmatrix(1), sab_aux_fit)
      CALL cp_dbcsr_alloc_block_from_nbl(ksmatrix(2), sab_aux_fit)

      kplocal = kp_range(2) - kp_range(1) + 1
      para_env => kpoints%blacs_env_all%para_env

      CALL cp_fm_struct_create(struct_aux_aux, context=kpoints%blacs_env, para_env=kpoints%para_env_kp, &
                               nrow_global=nao_aux_fit, ncol_global=nao_aux_fit)
      CALL cp_fm_create(work_aux_aux, struct_aux_aux)
      CALL cp_fm_create(work_aux_aux2, struct_aux_aux)

      CALL cp_fm_struct_create(struct_aux_orb, context=kpoints%blacs_env, para_env=kpoints%para_env_kp, &
                               nrow_global=nao_aux_fit, ncol_global=nao_orb)
      CALL cp_fm_create(work_aux_orb, struct_aux_orb)

      CALL cp_fm_struct_create(struct_orb_orb, context=kpoints%blacs_env, para_env=kpoints%para_env_kp, &
                               nrow_global=nao_orb, ncol_global=nao_orb)

      !Create cfm work matrices
      IF (.NOT. use_real_wfn) THEN
         CALL cp_cfm_create(cS, struct_aux_aux)
         CALL cp_cfm_create(cK, struct_aux_aux)
         CALL cp_cfm_create(cwork_aux_aux, struct_aux_aux)

         CALL cp_cfm_create(cA, struct_aux_orb)
         CALL cp_cfm_create(cwork_aux_orb, struct_aux_orb)

         CALL cp_cfm_create(cwork_orb_orb, struct_orb_orb)
      END IF

      !We create the fms in which we store the KS ORB matrix at each kp
      ALLOCATE (fm_ks(kplocal, 2, nspins))
      DO ispin = 1, nspins
         DO i = 1, 2
            DO ikp = 1, kplocal
               CALL cp_fm_create(fm_ks(ikp, i, ispin), struct_orb_orb)
            END DO
         END DO
      END DO

      CALL cp_fm_struct_release(struct_aux_aux)
      CALL cp_fm_struct_release(struct_aux_orb)
      CALL cp_fm_struct_release(struct_orb_orb)

      ALLOCATE (info(kplocal*nspins*nkp_groups, 2))
      indx = 0
      DO ikp = 1, kplocal
         DO ispin = 1, nspins
            DO igroup = 1, nkp_groups
               ! number of current kpoint
               ik = kp_dist(1, igroup) + ikp - 1
               my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
               indx = indx + 1

               IF (use_real_wfn) THEN
                  CALL dbcsr_set(ksmatrix(1), 0.0_dp)
                  CALL rskp_transform(rmatrix=ksmatrix(1), rsmat=matrix_ks_aux_fit, ispin=ispin, &
                                      xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_aux_fit)
                  CALL dbcsr_desymmetrize(ksmatrix(1), tmpmatrix_ks)
                  CALL copy_dbcsr_to_fm(tmpmatrix_ks, admm_env%work_aux_aux)
               ELSE
                  CALL dbcsr_set(ksmatrix(1), 0.0_dp)
                  CALL dbcsr_set(ksmatrix(2), 0.0_dp)
                  CALL rskp_transform(rmatrix=ksmatrix(1), cmatrix=ksmatrix(2), rsmat=matrix_ks_aux_fit, ispin=ispin, &
                                      xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_aux_fit)
                  CALL dbcsr_desymmetrize(ksmatrix(1), tmpmatrix_ks)
                  CALL copy_dbcsr_to_fm(tmpmatrix_ks, admm_env%work_aux_aux)
                  CALL dbcsr_desymmetrize(ksmatrix(2), tmpmatrix_ks)
                  CALL copy_dbcsr_to_fm(tmpmatrix_ks, admm_env%work_aux_aux2)
               END IF

               IF (my_kpgrp) THEN
                  CALL cp_fm_start_copy_general(admm_env%work_aux_aux, work_aux_aux, para_env, info(indx, 1))
                  IF (.NOT. use_real_wfn) &
                     CALL cp_fm_start_copy_general(admm_env%work_aux_aux2, work_aux_aux2, &
                                                   para_env, info(indx, 2))
               ELSE
                  CALL cp_fm_start_copy_general(admm_env%work_aux_aux, fmdummy, para_env, info(indx, 1))
                  IF (.NOT. use_real_wfn) &
                     CALL cp_fm_start_copy_general(admm_env%work_aux_aux2, fmdummy, para_env, info(indx, 2))
               END IF
            END DO
         END DO
      END DO

      indx = 0
      DO ikp = 1, kplocal
         DO ispin = 1, nspins
            DO igroup = 1, nkp_groups
               ! number of current kpoint
               ik = kp_dist(1, igroup) + ikp - 1
               my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
               indx = indx + 1
               IF (my_kpgrp) THEN
                  CALL cp_fm_finish_copy_general(work_aux_aux, info(indx, 1))
                  IF (.NOT. use_real_wfn) THEN
                     CALL cp_fm_finish_copy_general(work_aux_aux2, info(indx, 2))
                     CALL cp_fm_to_cfm(work_aux_aux, work_aux_aux2, cK)
                  END IF
               END IF
            END DO

            kp => kpoints%kp_aux_env(ikp)%kpoint_env
            IF (use_real_wfn) THEN

               !! K*A
               CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                                  1.0_dp, work_aux_aux, kp%amat(1, 1), 0.0_dp, &
                                  work_aux_orb)
               !! A^T*K*A
               CALL parallel_gemm('T', 'N', nao_orb, nao_orb, nao_aux_fit, &
                                  1.0_dp, kp%amat(1, 1), work_aux_orb, 0.0_dp, &
                                  fm_ks(ikp, 1, ispin))
            ELSE

               IF (admm_env%do_admmq .OR. admm_env%do_admms) THEN
                  CALL cp_fm_to_cfm(kp%smat(1, 1), kp%smat(2, 1), cS)

                  !Need to subdtract lambda* S_aux to K_aux, and scale the whole thing by gsi
                  fac = CMPLX(-admm_env%lambda_merlot(ispin), 0.0_dp, dp)
                  CALL cp_cfm_scale_and_add(z_one, cK, fac, cS)
                  CALL cp_cfm_scale(admm_env%gsi(ispin), cK)
               END IF

               IF (admm_env%do_admmp) THEN
                  CALL cp_fm_to_cfm(kp%smat(1, 1), kp%smat(2, 1), cS)

                  !Need to substract labda*gsi*S_aux to gsi**2*K_aux
                  fac = CMPLX(-admm_env%gsi(ispin)*admm_env%lambda_merlot(ispin), 0.0_dp, dp)
                  fac2 = CMPLX(admm_env%gsi(ispin)**2, 0.0_dp, dp)
                  CALL cp_cfm_scale_and_add(fac2, cK, fac, cS)
               END IF

               CALL cp_fm_to_cfm(kp%amat(1, 1), kp%amat(2, 1), cA)
               CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, &
                                  z_one, cK, cA, z_zero, cwork_aux_orb)

               CALL parallel_gemm('C', 'N', nao_orb, nao_orb, nao_aux_fit, &
                                  z_one, cA, cwork_aux_orb, z_zero, cwork_orb_orb)

               CALL cp_cfm_to_fm(cwork_orb_orb, mtargetr=fm_ks(ikp, 1, ispin), mtargeti=fm_ks(ikp, 2, ispin))
            END IF
         END DO
      END DO

      indx = 0
      DO ikp = 1, kplocal
         DO ispin = 1, nspins
            DO igroup = 1, nkp_groups
               ! number of current kpoint
               ik = kp_dist(1, igroup) + ikp - 1
               my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
               indx = indx + 1
               CALL cp_fm_cleanup_copy_general(info(indx, 1))
               IF (.NOT. use_real_wfn) CALL cp_fm_cleanup_copy_general(info(indx, 2))
            END DO
         END DO
      END DO

      DEALLOCATE (info)
      CALL dbcsr_release(ksmatrix(1))
      CALL dbcsr_release(ksmatrix(2))
      CALL dbcsr_release(tmpmatrix_ks)

      CALL cp_fm_release(work_aux_aux)
      CALL cp_fm_release(work_aux_aux2)
      CALL cp_fm_release(work_aux_orb)
      IF (.NOT. use_real_wfn) THEN
         CALL cp_cfm_release(cS)
         CALL cp_cfm_release(cK)
         CALL cp_cfm_release(cwork_aux_aux)
         CALL cp_cfm_release(cA)
         CALL cp_cfm_release(cwork_aux_orb)
         CALL cp_cfm_release(cwork_orb_orb)
      END IF

      NULLIFY (matrix_k_tilde)

      CALL dbcsr_allocate_matrix_set(matrix_k_tilde, dft_control%nspins, dft_control%nimages)

      DO ispin = 1, nspins
         DO img = 1, dft_control%nimages
            ALLOCATE (matrix_k_tilde(ispin, img)%matrix)
            CALL dbcsr_create(matrix=matrix_k_tilde(ispin, img)%matrix, template=matrix_ks_kp(1, 1)%matrix, &
                              name='MATRIX K_tilde '//TRIM(ADJUSTL(cp_to_string(ispin)))//'_'//TRIM(ADJUSTL(cp_to_string(img))), &
                              matrix_type=dbcsr_type_symmetric, nze=0)
            CALL cp_dbcsr_alloc_block_from_nbl(matrix_k_tilde(ispin, img)%matrix, sab_kp)
            CALL dbcsr_set(matrix_k_tilde(ispin, img)%matrix, 0.0_dp)
         END DO
      END DO

      CALL cp_fm_get_info(admm_env%work_orb_orb, matrix_struct=struct_orb_orb)
      ALLOCATE (fmwork(2))
      CALL cp_fm_create(fmwork(1), struct_orb_orb)
      CALL cp_fm_create(fmwork(2), struct_orb_orb)

      ! reuse the density transform to FT the KS matrix
      CALL kpoint_density_transform(kpoints, matrix_k_tilde, .FALSE., &
                                    matrix_k_tilde(1, 1)%matrix, sab_kp, &
                                    fmwork, for_aux_fit=.FALSE., pmat_ext=fm_ks)
      CALL cp_fm_release(fmwork(1))
      CALL cp_fm_release(fmwork(2))

      DO ispin = 1, nspins
         DO i = 1, 2
            DO ikp = 1, kplocal
               CALL cp_fm_release(fm_ks(ikp, i, ispin))
            END DO
         END DO
      END DO

      DO ispin = 1, nspins
         DO img = 1, dft_control%nimages
            CALL dbcsr_add(matrix_ks_kp(ispin, img)%matrix, matrix_k_tilde(ispin, img)%matrix, 1.0_dp, 1.0_dp)
            IF (admm_env%do_admmq .OR. admm_env%do_admmp .OR. admm_env%do_admms) THEN
               !In ADMMQ and ADMMP, need to add lambda*S_orb (Merlot eq 27)
               CALL dbcsr_add(matrix_ks_kp(ispin, img)%matrix, matrix_s(1, img)%matrix, &
                              1.0_dp, admm_env%lambda_merlot(ispin))
            END IF
         END DO
      END DO

      !Scale the energies
      IF (admm_env%do_admmp) THEN
         IF (nspins == 1) THEN
            energy%exc_aux_fit = (admm_env%gsi(1))**2.0_dp*energy%exc_aux_fit
            energy%exc1_aux_fit = (admm_env%gsi(1))**2.0_dp*energy%exc1_aux_fit
            energy%ex = (admm_env%gsi(1))**2.0_dp*energy%ex
         ELSE
            energy%exc_aux_fit = 0.0_dp
            energy%exc1_aux_fit = 0.0_dp
            energy%ex = 0.0_dp
            DO ispin = 1, dft_control%nspins
               energy%exc_aux_fit = energy%exc_aux_fit + (admm_env%gsi(ispin))**2.0_dp*ener_x(ispin)
               energy%exc1_aux_fit = energy%exc1_aux_fit + (admm_env%gsi(ispin))**2.0_dp*ener_x1(ispin)
               energy%ex = energy%ex + (admm_env%gsi(ispin))**2.0_dp*ener_k(ispin)
            END DO
         END IF
      END IF

      !Scale the energies and clean-up
      IF (admm_env%do_admms) THEN
         IF (nspins == 1) THEN
            energy%exc_aux_fit = (admm_env%gsi(1))**(2.0_dp/3.0_dp)*energy%exc_aux_fit
            energy%exc1_aux_fit = (admm_env%gsi(1))**(2.0_dp/3.0_dp)*energy%exc1_aux_fit
         ELSE
            energy%exc_aux_fit = 0.0_dp
            energy%exc1_aux_fit = 0.0_dp
            DO ispin = 1, nspins
               energy%exc_aux_fit = energy%exc_aux_fit + (admm_env%gsi(ispin))**(2.0_dp/3.0_dp)*ener_x(ispin)
               energy%exc1_aux_fit = energy%exc1_aux_fit + (admm_env%gsi(ispin))**(2.0_dp/3.0_dp)*ener_x1(ispin)
            END DO
         END IF

         CALL dbcsr_deallocate_matrix_set(matrix_ks_aux_fit)
      END IF

      CALL dbcsr_deallocate_matrix_set(matrix_k_tilde)

      CALL timestop(handle)

   END SUBROUTINE merge_ks_matrix_none_kp

! **************************************************************************************************
!> \brief Calculate exchange correction energy (Merlot2014 Eqs. 32, 33) for every spin, for KP
!> \param qs_env ...
!> \param admm_env ...
!> \param ener_k_ispin exact ispin (Fock) exchange in auxiliary basis
!> \param ener_x_ispin ispin DFT exchange in auxiliary basis
!> \param ener_x1_ispin ispin DFT exchange in auxiliary basis, due to the GAPW atomic contributions
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE calc_spin_dep_aux_exch_ener(qs_env, admm_env, ener_k_ispin, ener_x_ispin, &
                                          ener_x1_ispin, ispin)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(admm_type), POINTER                           :: admm_env
      REAL(dp), INTENT(INOUT)                            :: ener_k_ispin, ener_x_ispin, ener_x1_ispin
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_spin_dep_aux_exch_ener'

      CHARACTER(LEN=default_string_length)               :: basis_type
      INTEGER                                            :: handle, img, myspin, nimg
      LOGICAL                                            :: gapw
      REAL(dp)                                           :: tmp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(admm_gapw_r3d_rs_type), POINTER               :: admm_gapw_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_aux_fit_hfx, rho_ao_aux, &
                                                            rho_ao_aux_buffer
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(local_rho_type)                               :: local_rho_buffer
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r, v_rspace_dummy, v_tau_rspace_dummy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_aux_fit, rho_aux_fit_buffer
      TYPE(section_vals_type), POINTER                   :: xc_section_aux
      TYPE(task_list_type), POINTER                      :: task_list

      CALL timeset(routineN, handle)

      NULLIFY (ks_env, rho_aux_fit, rho_aux_fit_buffer, rho_ao, &
               xc_section_aux, v_rspace_dummy, v_tau_rspace_dummy, &
               rho_ao_aux, rho_ao_aux_buffer, dft_control, &
               matrix_ks_aux_fit_hfx, task_list, admm_gapw_env)

      NULLIFY (rho_g, rho_r, tot_rho_r)

      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control)
      CALL get_admm_env(admm_env, rho_aux_fit=rho_aux_fit, rho_aux_fit_buffer=rho_aux_fit_buffer, &
                        matrix_ks_aux_fit_hfx_kp=matrix_ks_aux_fit_hfx)

      CALL qs_rho_get(rho_aux_fit, &
                      rho_ao_kp=rho_ao_aux)

      CALL qs_rho_get(rho_aux_fit_buffer, &
                      rho_ao_kp=rho_ao_aux_buffer, &
                      rho_g=rho_g, &
                      rho_r=rho_r, &
                      tot_rho_r=tot_rho_r)

      gapw = admm_env%do_gapw
      nimg = dft_control%nimages

!   Calculate rho_buffer = rho_aux(ispin) to get exchange of ispin electrons
      DO img = 1, nimg
         CALL dbcsr_set(rho_ao_aux_buffer(1, img)%matrix, 0.0_dp)
         CALL dbcsr_set(rho_ao_aux_buffer(2, img)%matrix, 0.0_dp)
         CALL dbcsr_add(rho_ao_aux_buffer(ispin, img)%matrix, &
                        rho_ao_aux(ispin, img)%matrix, 0.0_dp, 1.0_dp)
      END DO

      ! By default use standard AUX_FIT basis and task_list. IF GAPW use the soft ones
      basis_type = "AUX_FIT"
      task_list => admm_env%task_list_aux_fit
      IF (gapw) THEN
         basis_type = "AUX_FIT_SOFT"
         task_list => admm_env%admm_gapw_env%task_list
      END IF

      ! integration for getting the spin dependent density has to done for both spins!
      DO myspin = 1, dft_control%nspins

         rho_ao => rho_ao_aux_buffer(myspin, :)
         CALL calculate_rho_elec(ks_env=ks_env, &
                                 matrix_p_kp=rho_ao, &
                                 rho=rho_r(myspin), &
                                 rho_gspace=rho_g(myspin), &
                                 total_rho=tot_rho_r(myspin), &
                                 soft_valid=.FALSE., &
                                 basis_type="AUX_FIT", &
                                 task_list_external=task_list)

      END DO

      ! Write changes in buffer density matrix
      CALL qs_rho_set(rho_aux_fit_buffer, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)

      xc_section_aux => admm_env%xc_section_aux

      ener_x_ispin = 0.0_dp

      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_aux_fit_buffer, xc_section=xc_section_aux, &
                         vxc_rho=v_rspace_dummy, vxc_tau=v_tau_rspace_dummy, exc=ener_x_ispin, &
                         just_energy=.TRUE.)

      !atomic contributions: use the atomic density as stored in admm_env%gapw_env
      ener_x1_ispin = 0.0_dp
      IF (gapw) THEN

         admm_gapw_env => admm_env%admm_gapw_env
         CALL get_qs_env(qs_env, &
                         atomic_kind_set=atomic_kind_set, &
                         para_env=para_env)

         CALL local_rho_set_create(local_rho_buffer)
         CALL allocate_rho_atom_internals(local_rho_buffer%rho_atom_set, atomic_kind_set, &
                                          admm_gapw_env%admm_kind_set, dft_control, para_env)

         CALL calculate_rho_atom_coeff(qs_env, rho_ao_aux_buffer, &
                                       rho_atom_set=local_rho_buffer%rho_atom_set, &
                                       qs_kind_set=admm_gapw_env%admm_kind_set, &
                                       oce=admm_gapw_env%oce, sab=admm_env%sab_aux_fit, &
                                       para_env=para_env)

         CALL prepare_gapw_den(qs_env, local_rho_set=local_rho_buffer, do_rho0=.FALSE., &
                               kind_set_external=admm_gapw_env%admm_kind_set)

         CALL calculate_vxc_atom(qs_env, energy_only=.TRUE., exc1=ener_x1_ispin, &
                                 kind_set_external=admm_env%admm_gapw_env%admm_kind_set, &
                                 xc_section_external=xc_section_aux, &
                                 rho_atom_set_external=local_rho_buffer%rho_atom_set)

         CALL local_rho_set_release(local_rho_buffer)
      END IF

      ener_k_ispin = 0.0_dp

      !! ** Calculate the exchange energy
      DO img = 1, nimg
         CALL dbcsr_dot(matrix_ks_aux_fit_hfx(ispin, img)%matrix, rho_ao_aux_buffer(ispin, img)%matrix, tmp)
         ener_k_ispin = ener_k_ispin + tmp
      END DO

      ! Divide exchange for indivivual spin by two, since the ener_k_ispin originally is total
      ! exchange of alpha and beta
      ener_k_ispin = ener_k_ispin/2.0_dp

      CALL timestop(handle)

   END SUBROUTINE calc_spin_dep_aux_exch_ener

! **************************************************************************************************
!> \brief Scale density matrix by gsi(ispin), is needed for force scaling in ADMMP
!> \param qs_env ...
!> \param rho_ao_orb ...
!> \param scale_back ...
!> \author Jan Wilhelm, 12/2014
! **************************************************************************************************
   SUBROUTINE scale_dm(qs_env, rho_ao_orb, scale_back)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_orb
      LOGICAL, INTENT(IN)                                :: scale_back

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'scale_dm'

      INTEGER                                            :: handle, img, ispin
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      NULLIFY (admm_env, dft_control)

      CALL get_qs_env(qs_env, &
                      admm_env=admm_env, &
                      dft_control=dft_control)

      ! only for ADMMP
      IF (admm_env%do_admmp) THEN
         DO ispin = 1, dft_control%nspins
            DO img = 1, dft_control%nimages
               IF (scale_back) THEN
                  CALL dbcsr_scale(rho_ao_orb(ispin, img)%matrix, 1.0_dp/admm_env%gsi(ispin))
               ELSE
                  CALL dbcsr_scale(rho_ao_orb(ispin, img)%matrix, admm_env%gsi(ispin))
               END IF
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE scale_dm

! **************************************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param mo_set ...
!> \param mo_coeff_aux_fit ...
! **************************************************************************************************
   SUBROUTINE calc_aux_mo_derivs_none(ispin, admm_env, mo_set, mo_coeff_aux_fit)
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(mo_set_type), INTENT(IN)                      :: mo_set
      TYPE(cp_fm_type), INTENT(IN)                       :: mo_coeff_aux_fit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_aux_mo_derivs_none'

      INTEGER                                            :: handle, nao_aux_fit, nao_orb, nmo
      REAL(dp), DIMENSION(:), POINTER                    :: occupation_numbers, scaling_factor
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks_aux_fit, &
                                                            matrix_ks_aux_fit_dft, &
                                                            matrix_ks_aux_fit_hfx
      TYPE(dbcsr_type)                                   :: dbcsr_work

      NULLIFY (matrix_ks_aux_fit, matrix_ks_aux_fit_dft, matrix_ks_aux_fit_hfx)

      CALL timeset(routineN, handle)

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb
      nmo = admm_env%nmo(ispin)

      CALL get_admm_env(admm_env, matrix_ks_aux_fit=matrix_ks_aux_fit, &
                        matrix_ks_aux_fit_hfx=matrix_ks_aux_fit_hfx, &
                        matrix_ks_aux_fit_dft=matrix_ks_aux_fit_dft)

      ! just calculate the mo derivs in the aux basis
      ! only needs to be done on the converged ks matrix for the force calc
      ! Note with OT and purification NONE, the merging of the derivs
      ! happens implicitly because the KS matrices have been already been merged
      ! and adding them here would be double counting.

      IF (admm_env%do_admms) THEN
         !In ADMMS, we use the K matrix defined as K_hf - gsi^2/3*K_dft
         CALL dbcsr_create(dbcsr_work, template=matrix_ks_aux_fit(ispin)%matrix)
         CALL dbcsr_copy(dbcsr_work, matrix_ks_aux_fit_hfx(ispin)%matrix)
         CALL dbcsr_add(dbcsr_work, matrix_ks_aux_fit_dft(ispin)%matrix, 1.0_dp, -admm_env%gsi(ispin)**(2.0_dp/3.0_dp))
         CALL copy_dbcsr_to_fm(dbcsr_work, admm_env%K(ispin))
         CALL dbcsr_release(dbcsr_work)
      ELSE
         CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix, admm_env%K(ispin))
      END IF
      CALL cp_fm_upper_to_full(admm_env%K(ispin), admm_env%work_aux_aux)

      CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nao_aux_fit, &
                         1.0_dp, admm_env%K(ispin), mo_coeff_aux_fit, 0.0_dp, &
                         admm_env%H(ispin))

      CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
      ALLOCATE (scaling_factor(SIZE(occupation_numbers)))

      scaling_factor = 2.0_dp*occupation_numbers

      CALL cp_fm_column_scale(admm_env%H(ispin), scaling_factor)

      DEALLOCATE (scaling_factor)

      CALL timestop(handle)

   END SUBROUTINE calc_aux_mo_derivs_none

! **************************************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param mo_set ...
!> \param mo_derivs ...
!> \param matrix_ks_aux_fit ...
! **************************************************************************************************
   SUBROUTINE merge_mo_derivs_no_diag(ispin, admm_env, mo_set, mo_derivs, matrix_ks_aux_fit)
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(mo_set_type), INTENT(IN)                      :: mo_set
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_derivs
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks_aux_fit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_mo_derivs_no_diag'

      INTEGER                                            :: handle, nao_aux_fit, nao_orb, nmo
      REAL(dp), DIMENSION(:), POINTER                    :: occupation_numbers, scaling_factor

      CALL timeset(routineN, handle)

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb
      nmo = admm_env%nmo(ispin)

      CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix, admm_env%K(ispin))
      CALL cp_fm_upper_to_full(admm_env%K(ispin), admm_env%work_aux_aux)

      CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
      ALLOCATE (scaling_factor(SIZE(occupation_numbers)))
      scaling_factor = 0.5_dp

      !! ** calculate first part
      CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nmo, &
                         1.0_dp, admm_env%C_hat(ispin), admm_env%lambda_inv(ispin), 0.0_dp, &
                         admm_env%work_aux_nmo(ispin))
      CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nao_aux_fit, &
                         1.0_dp, admm_env%K(ispin), admm_env%work_aux_nmo(ispin), 0.0_dp, &
                         admm_env%work_aux_nmo2(ispin))
      CALL parallel_gemm('T', 'N', nao_orb, nmo, nao_aux_fit, &
                         2.0_dp, admm_env%A, admm_env%work_aux_nmo2(ispin), 0.0_dp, &
                         admm_env%mo_derivs_tmp(ispin))
      !! ** calculate second part
      CALL parallel_gemm('T', 'N', nmo, nmo, nao_aux_fit, &
                         1.0_dp, admm_env%work_aux_nmo(ispin), admm_env%work_aux_nmo2(ispin), 0.0_dp, &
                         admm_env%work_orb_orb)
      CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nmo, &
                         1.0_dp, admm_env%C_hat(ispin), admm_env%work_orb_orb, 0.0_dp, &
                         admm_env%work_aux_orb)
      CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nao_aux_fit, &
                         1.0_dp, admm_env%S, admm_env%work_aux_orb, 0.0_dp, &
                         admm_env%work_aux_nmo(ispin))
      CALL parallel_gemm('T', 'N', nao_orb, nmo, nao_aux_fit, &
                         -2.0_dp, admm_env%A, admm_env%work_aux_nmo(ispin), 1.0_dp, &
                         admm_env%mo_derivs_tmp(ispin))

      CALL cp_fm_column_scale(admm_env%mo_derivs_tmp(ispin), scaling_factor)

      CALL cp_fm_scale_and_add(1.0_dp, mo_derivs(ispin), 1.0_dp, admm_env%mo_derivs_tmp(ispin))

      DEALLOCATE (scaling_factor)

      CALL timestop(handle)

   END SUBROUTINE merge_mo_derivs_no_diag

! **************************************************************************************************
!> \brief Calculate the derivative of the AUX_FIT mo, based on the ORB mo_derivs
!> \param qs_env ...
!> \param mo_derivs the MO derivatives in the orbital basis
! **************************************************************************************************
   SUBROUTINE calc_admm_mo_derivatives(qs_env, mo_derivs)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mo_derivs

      INTEGER                                            :: ispin, nspins
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: mo_derivs_fm
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mo_derivs_aux_fit
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_aux_fit
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks_aux_fit
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mo_array, mos_aux_fit

      NULLIFY (mo_array, mos_aux_fit, matrix_ks_aux_fit, mo_coeff_aux_fit, &
               mo_derivs_aux_fit, mo_coeff)

      CALL get_qs_env(qs_env, admm_env=admm_env, mos=mo_array)
      CALL get_admm_env(admm_env, mos_aux_fit=mos_aux_fit, mo_derivs_aux_fit=mo_derivs_aux_fit, &
                        matrix_ks_aux_fit=matrix_ks_aux_fit)

      nspins = SIZE(mo_derivs)
      ALLOCATE (mo_derivs_fm(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mo_array(ispin), mo_coeff=mo_coeff)
         CALL cp_fm_create(mo_derivs_fm(ispin), mo_coeff%matrix_struct)
      END DO

      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mo_array(ispin), mo_coeff=mo_coeff)
         CALL get_mo_set(mo_set=mos_aux_fit(ispin), mo_coeff=mo_coeff_aux_fit)

         CALL copy_dbcsr_to_fm(mo_derivs(ispin)%matrix, mo_derivs_fm(ispin))
         CALL admm_mo_merge_derivs(ispin, admm_env, mo_array(ispin), mo_coeff, mo_coeff_aux_fit, &
                                   mo_derivs_fm, mo_derivs_aux_fit, matrix_ks_aux_fit)
         CALL copy_fm_to_dbcsr(mo_derivs_fm(ispin), mo_derivs(ispin)%matrix)
      END DO

      CALL cp_fm_release(mo_derivs_fm)

   END SUBROUTINE calc_admm_mo_derivatives

! **************************************************************************************************
!> \brief Calculate the forces due to the AUX/ORB basis overlap in ADMM
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calc_admm_ovlp_forces(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: ispin
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_aux_fit
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s_aux_fit, matrix_s_aux_fit_vs_orb
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos, mos_aux_fit
      TYPE(mo_set_type), POINTER                         :: mo_set

      CALL get_qs_env(qs_env, dft_control=dft_control)

      IF (dft_control%do_admm_dm) THEN
         CPABORT("Forces with ADMM DM methods not implemented")
      END IF
      IF (dft_control%do_admm_mo .AND. .NOT. qs_env%run_rtp) THEN
         NULLIFY (matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, mos_aux_fit, mos, admm_env)
         CALL get_qs_env(qs_env=qs_env, &
                         mos=mos, &
                         admm_env=admm_env)
         CALL get_admm_env(admm_env, matrix_s_aux_fit=matrix_s_aux_fit, mos_aux_fit=mos_aux_fit, &
                           matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb)
         DO ispin = 1, dft_control%nspins
            mo_set => mos(ispin)
            CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff)
            ! if no purification we need to calculate the H matrix for forces
            IF (admm_env%purification_method == do_admm_purify_none) THEN
               CALL get_mo_set(mo_set=mos_aux_fit(ispin), mo_coeff=mo_coeff_aux_fit)
               CALL calc_aux_mo_derivs_none(ispin, qs_env%admm_env, mo_set, mo_coeff_aux_fit)
            END IF
         END DO
         CALL calc_mixed_overlap_force(qs_env)
      END IF

   END SUBROUTINE calc_admm_ovlp_forces

! **************************************************************************************************
!> \brief Calculate the forces due to the AUX/ORB basis overlap in ADMM, in the KP case
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calc_admm_ovlp_forces_kp(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_admm_ovlp_forces_kp'

      COMPLEX(dp)                                        :: fac, fac2
      INTEGER                                            :: handle, i, igroup, ik, ikp, img, indx, &
                                                            ispin, kplocal, nao_aux_fit, nao_orb, &
                                                            natom, nimg, nkp, nkp_groups, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: gapw, my_kpgrp, use_real_wfn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: admm_force
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(copy_info_type), ALLOCATABLE, DIMENSION(:, :) :: info
      TYPE(cp_cfm_type)                                  :: cA, ckmatrix, cpmatrix, cQ, cS, cS_inv, &
                                                            cwork_aux_aux, cwork_aux_orb, &
                                                            cwork_aux_orb2
      TYPE(cp_fm_struct_type), POINTER                   :: struct_aux_aux, struct_aux_orb, &
                                                            struct_orb_orb
      TYPE(cp_fm_type)                                   :: fmdummy, S_inv, work_aux_aux, &
                                                            work_aux_aux2, work_aux_aux3, &
                                                            work_aux_orb
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_skap, fm_skapa
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: fmwork
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER :: matrix_ks_aux_fit, matrix_ks_aux_fit_dft, &
         matrix_ks_aux_fit_hfx, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, matrix_skap, &
         matrix_skapa, rho_ao_orb
      TYPE(dbcsr_type)                                   :: kmatrix_tmp
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: kmatrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux_fit, sab_aux_fit_asymm, &
                                                            sab_aux_fit_vs_orb, sab_kp
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      !Note: we only treat the case with purification none, there the overlap forces read as:
      !F = 2*Tr[P * A^T * K_aux * S^-1_aux * Q^(x)] - 2*Tr[A * P * A^T * K_aux * S^-1_aux *S_aux^(x)]
      !where P is the density matrix in the ORB basis. As a strategy, we FT all relevant matrices
      !from real space to KP, calculate the matrix products, back FT to real space, and calculate the
      !overlap forces

      NULLIFY (ks_env, admm_env, matrix_ks_aux_fit, &
               matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, rho, force, &
               para_env, atomic_kind_set, kpoints, sab_aux_fit, &
               sab_aux_fit_vs_orb, sab_aux_fit_asymm, struct_orb_orb, &
               struct_aux_orb, struct_aux_aux)

      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      admm_env=admm_env, &
                      dft_control=dft_control, &
                      kpoints=kpoints, &
                      natom=natom, &
                      atomic_kind_set=atomic_kind_set, &
                      force=force, &
                      rho=rho)
      nimg = dft_control%nimages
      CALL get_admm_env(admm_env, &
                        matrix_s_aux_fit_kp=matrix_s_aux_fit, &
                        matrix_s_aux_fit_vs_orb_kp=matrix_s_aux_fit_vs_orb, &
                        sab_aux_fit=sab_aux_fit, &
                        sab_aux_fit_vs_orb=sab_aux_fit_vs_orb, &
                        sab_aux_fit_asymm=sab_aux_fit_asymm, &
                        matrix_ks_aux_fit_kp=matrix_ks_aux_fit, &
                        matrix_ks_aux_fit_dft_kp=matrix_ks_aux_fit_dft, &
                        matrix_ks_aux_fit_hfx_kp=matrix_ks_aux_fit_hfx)

      gapw = admm_env%do_gapw
      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb
      nspins = dft_control%nspins

      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, use_real_wfn=use_real_wfn, kp_range=kp_range, &
                           nkp_groups=nkp_groups, kp_dist=kp_dist, &
                           cell_to_index=cell_to_index, sab_nl=sab_kp)

      !Case study on ADMMQ, ADMMS and ADMMP
      IF (admm_env%do_admms) THEN
         !Here we buld the KS matrix: KS_hfx = gsi^2/3*KS_dft, the we then pass as the ususal KS_aux_fit
         NULLIFY (matrix_ks_aux_fit)
         ALLOCATE (matrix_ks_aux_fit(nspins, dft_control%nimages))
         DO img = 1, dft_control%nimages
            DO ispin = 1, nspins
               NULLIFY (matrix_ks_aux_fit(ispin, img)%matrix)
               ALLOCATE (matrix_ks_aux_fit(ispin, img)%matrix)
               CALL dbcsr_create(matrix_ks_aux_fit(ispin, img)%matrix, template=matrix_s_aux_fit(1, 1)%matrix)
               CALL dbcsr_copy(matrix_ks_aux_fit(ispin, img)%matrix, matrix_ks_aux_fit_hfx(ispin, img)%matrix)
               CALL dbcsr_add(matrix_ks_aux_fit(ispin, img)%matrix, matrix_ks_aux_fit_dft(ispin, img)%matrix, &
                              1.0_dp, -admm_env%gsi(ispin)**(2.0_dp/3.0_dp))
            END DO
         END DO
      END IF

      ! the temporary DBCSR matrices for the rskp_transform we have to manually allocate
      ! index 1 => real, index 2 => imaginary
      ALLOCATE (kmatrix(2))
      CALL dbcsr_create(kmatrix(1), template=matrix_ks_aux_fit(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(kmatrix(2), template=matrix_ks_aux_fit(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL dbcsr_create(kmatrix_tmp, template=matrix_ks_aux_fit(1, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(kmatrix(1), sab_aux_fit)
      CALL cp_dbcsr_alloc_block_from_nbl(kmatrix(2), sab_aux_fit)

      kplocal = kp_range(2) - kp_range(1) + 1
      para_env => kpoints%blacs_env_all%para_env
      ALLOCATE (info(kplocal*nspins*nkp_groups, 2))

      CALL cp_fm_struct_create(struct_aux_aux, context=kpoints%blacs_env, para_env=kpoints%para_env_kp, &
                               nrow_global=nao_aux_fit, ncol_global=nao_aux_fit)
      CALL cp_fm_create(work_aux_aux, struct_aux_aux)
      CALL cp_fm_create(work_aux_aux2, struct_aux_aux)
      CALL cp_fm_create(work_aux_aux3, struct_aux_aux)
      CALL cp_fm_create(s_inv, struct_aux_aux)

      CALL cp_fm_struct_create(struct_aux_orb, context=kpoints%blacs_env, para_env=kpoints%para_env_kp, &
                               nrow_global=nao_aux_fit, ncol_global=nao_orb)
      CALL cp_fm_create(work_aux_orb, struct_aux_orb)

      CALL cp_fm_struct_create(struct_orb_orb, context=kpoints%blacs_env, para_env=kpoints%para_env_kp, &
                               nrow_global=nao_orb, ncol_global=nao_orb)

      !Create cfm work matrices
      IF (.NOT. use_real_wfn) THEN
         CALL cp_cfm_create(cpmatrix, struct_orb_orb)

         CALL cp_cfm_create(cS_inv, struct_aux_aux)
         CALL cp_cfm_create(cS, struct_aux_aux)
         CALL cp_cfm_create(cwork_aux_aux, struct_aux_aux)
         CALL cp_cfm_create(ckmatrix, struct_aux_aux)

         CALL cp_cfm_create(cA, struct_aux_orb)
         CALL cp_cfm_create(cQ, struct_aux_orb)
         CALL cp_cfm_create(cwork_aux_orb, struct_aux_orb)
         CALL cp_cfm_create(cwork_aux_orb2, struct_aux_orb)
      END IF

      !We create the fms in which we store the KP matrix products
      ALLOCATE (fm_skap(kplocal, 2, nspins), fm_skapa(kplocal, 2, nspins))
      DO ispin = 1, nspins
         DO i = 1, 2
            DO ikp = 1, kplocal
               CALL cp_fm_create(fm_skap(ikp, i, ispin), struct_aux_orb)
               CALL cp_fm_create(fm_skapa(ikp, i, ispin), struct_aux_aux)
            END DO
         END DO
      END DO

      CALL cp_fm_struct_release(struct_aux_aux)
      CALL cp_fm_struct_release(struct_aux_orb)
      CALL cp_fm_struct_release(struct_orb_orb)

      indx = 0
      DO ikp = 1, kplocal
         DO ispin = 1, nspins
            DO igroup = 1, nkp_groups
               ! number of current kpoint
               ik = kp_dist(1, igroup) + ikp - 1
               my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
               indx = indx + 1

               ! FT of matrices KS, then transfer to FM type
               IF (use_real_wfn) THEN
                  CALL dbcsr_set(kmatrix(1), 0.0_dp)
                  CALL rskp_transform(rmatrix=kmatrix(1), rsmat=matrix_ks_aux_fit, ispin=ispin, &
                                      xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_aux_fit)
                  CALL dbcsr_desymmetrize(kmatrix(1), kmatrix_tmp)
                  CALL copy_dbcsr_to_fm(kmatrix_tmp, admm_env%work_aux_aux)
               ELSE
                  CALL dbcsr_set(kmatrix(1), 0.0_dp)
                  CALL dbcsr_set(kmatrix(2), 0.0_dp)
                  CALL rskp_transform(rmatrix=kmatrix(1), cmatrix=kmatrix(2), rsmat=matrix_ks_aux_fit, ispin=ispin, &
                                      xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_aux_fit)
                  CALL dbcsr_desymmetrize(kmatrix(1), kmatrix_tmp)
                  CALL copy_dbcsr_to_fm(kmatrix_tmp, admm_env%work_aux_aux)
                  CALL dbcsr_desymmetrize(kmatrix(2), kmatrix_tmp)
                  CALL copy_dbcsr_to_fm(kmatrix_tmp, admm_env%work_aux_aux2)
               END IF

               IF (my_kpgrp) THEN
                  CALL cp_fm_start_copy_general(admm_env%work_aux_aux, work_aux_aux, para_env, info(indx, 1))
                  IF (.NOT. use_real_wfn) &
                     CALL cp_fm_start_copy_general(admm_env%work_aux_aux2, work_aux_aux2, para_env, info(indx, 2))
               ELSE
                  CALL cp_fm_start_copy_general(admm_env%work_aux_aux, fmdummy, para_env, info(indx, 1))
                  IF (.NOT. use_real_wfn) &
                     CALL cp_fm_start_copy_general(admm_env%work_aux_aux2, fmdummy, para_env, info(indx, 2))
               END IF
            END DO
         END DO
      END DO

      indx = 0
      DO ikp = 1, kplocal
         DO ispin = 1, nspins
            DO igroup = 1, nkp_groups
               ! number of current kpoint
               ik = kp_dist(1, igroup) + ikp - 1
               my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
               indx = indx + 1
               IF (my_kpgrp) THEN
                  CALL cp_fm_finish_copy_general(work_aux_aux, info(indx, 1))
                  IF (.NOT. use_real_wfn) THEN
                     CALL cp_fm_finish_copy_general(work_aux_aux2, info(indx, 2))
                     CALL cp_fm_to_cfm(work_aux_aux, work_aux_aux2, ckmatrix)
                  END IF
               END IF
            END DO
            kp => kpoints%kp_aux_env(ikp)%kpoint_env

            IF (use_real_wfn) THEN

               !! Calculate S'_inverse
               CALL cp_fm_to_fm(kp%smat(1, 1), S_inv)
               CALL cp_fm_cholesky_decompose(S_inv)
               CALL cp_fm_cholesky_invert(S_inv)
               !! Symmetrize the guy
               CALL cp_fm_upper_to_full(S_inv, work_aux_aux3)

               !We need to calculate S^-1*K*A*P and S^-1*K*A*P*A^T
               CALL parallel_gemm('N', 'N', nao_aux_fit, nao_aux_fit, nao_aux_fit, 1.0_dp, S_inv, &
                                  work_aux_aux, 0.0_dp, work_aux_aux3) ! S^-1 * K
               CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, 1.0_dp, work_aux_aux3, &
                                  kp%amat(1, 1), 0.0_dp, work_aux_orb) ! S^-1 * K * A
               CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb, 1.0_dp, work_aux_orb, &
                                  kpoints%kp_env(ikp)%kpoint_env%pmat(1, ispin), 0.0_dp, &
                                  fm_skap(ikp, 1, ispin)) ! S^-1 * K * A * P
               CALL parallel_gemm('N', 'T', nao_aux_fit, nao_aux_fit, nao_orb, 1.0_dp, fm_skap(ikp, 1, ispin), &
                                  kp%amat(1, 1), 0.0_dp, fm_skapa(ikp, 1, ispin))

            ELSE !complex wfn

               IF (admm_env%do_admmq .OR. admm_env%do_admms) THEN
                  CALL cp_fm_to_cfm(kp%smat(1, 1), kp%smat(2, 1), cS)

                  !Need to subdtract lambda* S_aux to K_aux, and scale the whole thing by gsi
                  fac = CMPLX(-admm_env%lambda_merlot(ispin), 0.0_dp, dp)
                  CALL cp_cfm_scale_and_add(z_one, ckmatrix, fac, cS)
                  CALL cp_cfm_scale(admm_env%gsi(ispin), ckmatrix)
               END IF

               IF (admm_env%do_admmp) THEN
                  CALL cp_fm_to_cfm(kp%smat(1, 1), kp%smat(2, 1), cS)

                  !Need to substract labda*gsi*S_aux to gsi**2*K_aux
                  fac = CMPLX(-admm_env%gsi(ispin)*admm_env%lambda_merlot(ispin), 0.0_dp, dp)
                  fac2 = CMPLX(admm_env%gsi(ispin)**2, 0.0_dp, dp)
                  CALL cp_cfm_scale_and_add(fac2, ckmatrix, fac, cS)
               END IF

               CALL cp_fm_to_cfm(kp%smat(1, 1), kp%smat(2, 1), cS_inv)
               CALL cp_cfm_cholesky_decompose(cS_inv)
               CALL cp_cfm_cholesky_invert(cS_inv)
               CALL own_cfm_upper_to_full(cS_inv, cwork_aux_aux)

               !Take the ORB density matrix from the kp_env
               CALL cp_fm_to_cfm(kpoints%kp_env(ikp)%kpoint_env%pmat(1, ispin), &
                                 kpoints%kp_env(ikp)%kpoint_env%pmat(2, ispin), &
                                 cpmatrix)

               !Do the same thing as in the real case
               !We need to calculate S^-1*K*A*P and S^-1*K*A*P*A^T
               CALL cp_fm_to_cfm(kp%amat(1, 1), kp%amat(2, 1), cA)
               CALL parallel_gemm('N', 'N', nao_aux_fit, nao_aux_fit, nao_aux_fit, z_one, cS_inv, &
                                  ckmatrix, z_zero, cwork_aux_aux) ! S^-1 * K
               CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, z_one, cwork_aux_aux, &
                                  cA, z_zero, cwork_aux_orb) ! S^-1 * K * A
               CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb, z_one, cwork_aux_orb, &
                                  cpmatrix, z_zero, cwork_aux_orb2) ! S^-1 * K * A * P
               CALL parallel_gemm('N', 'C', nao_aux_fit, nao_aux_fit, nao_orb, z_one, cwork_aux_orb2, &
                                  cA, z_zero, cwork_aux_aux)

               IF (admm_env%do_admmq .OR. admm_env%do_admmp .OR. admm_env%do_admms) THEN
                  !In ADMMQ, ADMMS, and ADMMP, there is an extra lambda*Tq *P* Tq^T matrix to contract with S_aux^(x)
                  !we calculate it and add it to fm_skapa (aka cwork_aux_aux)

                  !factor 0.5 because later multiplied by 2
                  fac = CMPLX(0.5_dp*admm_env%lambda_merlot(ispin)*admm_env%gsi(ispin), 0.0_dp, dp)
                  CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb, z_one, cA, cpmatrix, &
                                     z_zero, cwork_aux_orb)
                  CALL parallel_gemm('N', 'C', nao_aux_fit, nao_aux_fit, nao_orb, fac, cwork_aux_orb, &
                                     cA, z_one, cwork_aux_aux)
               END IF

               CALL cp_cfm_to_fm(cwork_aux_orb2, mtargetr=fm_skap(ikp, 1, ispin), mtargeti=fm_skap(ikp, 2, ispin))
               CALL cp_cfm_to_fm(cwork_aux_aux, mtargetr=fm_skapa(ikp, 1, ispin), mtargeti=fm_skapa(ikp, 2, ispin))

            END IF

         END DO
      END DO

      indx = 0
      DO ikp = 1, kplocal
         DO ispin = 1, nspins
            DO igroup = 1, nkp_groups
               ! number of current kpoint
               ik = kp_dist(1, igroup) + ikp - 1
               my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
               indx = indx + 1
               CALL cp_fm_cleanup_copy_general(info(indx, 1))
               IF (.NOT. use_real_wfn) CALL cp_fm_cleanup_copy_general(info(indx, 2))
            END DO
         END DO
      END DO

      DEALLOCATE (info)
      CALL dbcsr_release(kmatrix(1))
      CALL dbcsr_release(kmatrix(2))
      CALL dbcsr_release(kmatrix_tmp)

      CALL cp_fm_release(work_aux_aux)
      CALL cp_fm_release(work_aux_aux2)
      CALL cp_fm_release(work_aux_aux3)
      CALL cp_fm_release(S_inv)
      CALL cp_fm_release(work_aux_orb)
      IF (.NOT. use_real_wfn) THEN
         CALL cp_cfm_release(ckmatrix)
         CALL cp_cfm_release(cpmatrix)
         CALL cp_cfm_release(cS_inv)
         CALL cp_cfm_release(cS)
         CALL cp_cfm_release(cwork_aux_aux)
         CALL cp_cfm_release(cwork_aux_orb)
         CALL cp_cfm_release(cwork_aux_orb2)
         CALL cp_cfm_release(cA)
         CALL cp_cfm_release(cQ)
      END IF

      !Back FT to real space
      ALLOCATE (matrix_skap(nspins, nimg), matrix_skapa(nspins, nimg))
      DO img = 1, nimg
         DO ispin = 1, nspins
            ALLOCATE (matrix_skap(ispin, img)%matrix)
            CALL dbcsr_create(matrix_skap(ispin, img)%matrix, template=matrix_s_aux_fit_vs_orb(1, 1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_alloc_block_from_nbl(matrix_skap(ispin, img)%matrix, sab_aux_fit_vs_orb)

            ALLOCATE (matrix_skapa(ispin, img)%matrix)
            CALL dbcsr_create(matrix_skapa(ispin, img)%matrix, template=matrix_s_aux_fit(1, 1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_alloc_block_from_nbl(matrix_skapa(ispin, img)%matrix, sab_aux_fit_asymm)
         END DO
      END DO

      ALLOCATE (fmwork(2))
      CALL cp_fm_get_info(admm_env%work_aux_orb, matrix_struct=struct_aux_orb)
      CALL cp_fm_create(fmwork(1), struct_aux_orb)
      CALL cp_fm_create(fmwork(2), struct_aux_orb)
      CALL kpoint_density_transform(kpoints, matrix_skap, .FALSE., &
                                    matrix_s_aux_fit_vs_orb(1, 1)%matrix, sab_aux_fit_vs_orb, &
                                    fmwork, for_aux_fit=.TRUE., pmat_ext=fm_skap)
      CALL cp_fm_release(fmwork(1))
      CALL cp_fm_release(fmwork(2))

      CALL cp_fm_get_info(admm_env%work_aux_aux, matrix_struct=struct_aux_aux)
      CALL cp_fm_create(fmwork(1), struct_aux_aux)
      CALL cp_fm_create(fmwork(2), struct_aux_aux)
      CALL kpoint_density_transform(kpoints, matrix_skapa, .FALSE., &
                                    matrix_s_aux_fit(1, 1)%matrix, sab_aux_fit_asymm, &
                                    fmwork, for_aux_fit=.TRUE., pmat_ext=fm_skapa)
      CALL cp_fm_release(fmwork(1))
      CALL cp_fm_release(fmwork(2))
      DEALLOCATE (fmwork)

      DO img = 1, nimg
         DO ispin = 1, nspins
            CALL dbcsr_scale(matrix_skap(ispin, img)%matrix, -2.0_dp)
            CALL dbcsr_scale(matrix_skapa(ispin, img)%matrix, 2.0_dp)
         END DO
         IF (nspins == 2) THEN
            CALL dbcsr_add(matrix_skap(1, img)%matrix, matrix_skap(2, img)%matrix, 1.0_dp, 1.0_dp)
            CALL dbcsr_add(matrix_skapa(1, img)%matrix, matrix_skapa(2, img)%matrix, 1.0_dp, 1.0_dp)
         END IF
      END DO

      ALLOCATE (admm_force(3, natom))
      admm_force = 0.0_dp

      IF (admm_env%do_admmq .OR. admm_env%do_admmp .OR. admm_env%do_admms) THEN
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao_orb)
         DO img = 1, nimg
            DO ispin = 1, nspins
               CALL dbcsr_scale(rho_ao_orb(ispin, img)%matrix, -admm_env%lambda_merlot(ispin))
            END DO
            IF (nspins == 2) CALL dbcsr_add(rho_ao_orb(1, img)%matrix, rho_ao_orb(2, img)%matrix, 1.0_dp, 1.0_dp)
         END DO

         !In ADMMQ, ADMMS and ADMMP, there is an extra contribution from lambda*P_orb*S^(x)
         CALL build_overlap_force(qs_env%ks_env, admm_force, basis_type_a="ORB", basis_type_b="ORB", &
                                  sab_nl=sab_kp, matrixkp_p=rho_ao_orb(1, :))
         DO img = 1, nimg
            IF (nspins == 2) CALL dbcsr_add(rho_ao_orb(1, img)%matrix, rho_ao_orb(2, img)%matrix, 1.0_dp, -1.0_dp)
            DO ispin = 1, nspins
               CALL dbcsr_scale(rho_ao_orb(ispin, img)%matrix, -1.0_dp/admm_env%lambda_merlot(ispin))
            END DO
         END DO
      END IF

      CALL build_overlap_force(qs_env%ks_env, admm_force, basis_type_a="AUX_FIT", basis_type_b="ORB", &
                               sab_nl=sab_aux_fit_vs_orb, matrixkp_p=matrix_skap(1, :))
      CALL build_overlap_force(qs_env%ks_env, admm_force, basis_type_a="AUX_FIT", basis_type_b="AUX_FIT", &
                               sab_nl=sab_aux_fit_asymm, matrixkp_p=matrix_skapa(1, :))

      CALL add_qs_force(admm_force, force, "overlap_admm", atomic_kind_set)
      DEALLOCATE (admm_force)

      DO ispin = 1, nspins
         DO i = 1, 2
            DO ikp = 1, kplocal
               CALL cp_fm_release(fm_skap(ikp, i, ispin))
               CALL cp_fm_release(fm_skapa(ikp, i, ispin))
            END DO
         END DO
      END DO
      CALL dbcsr_deallocate_matrix_set(matrix_skap)
      CALL dbcsr_deallocate_matrix_set(matrix_skapa)

      IF (admm_env%do_admms) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_ks_aux_fit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calc_admm_ovlp_forces_kp

! **************************************************************************************************
!> \brief Calculate derivatives terms from overlap matrices
!> \param qs_env ...
!> \param matrix_hz Fock matrix part using the response density in admm basis
!> \param matrix_pz response density in orbital basis
!> \param fval ...
! **************************************************************************************************
   SUBROUTINE admm_projection_derivative(qs_env, matrix_hz, matrix_pz, fval)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: matrix_hz, matrix_pz
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: fval

      CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_projection_derivative'

      INTEGER                                            :: handle, ispin, nao, natom, naux, nspins
      REAL(KIND=dp)                                      :: my_fval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: admm_force
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s_aux_fit, matrix_s_aux_fit_vs_orb
      TYPE(dbcsr_type), POINTER                          :: matrix_w_q, matrix_w_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux_fit_asymm, sab_aux_fit_vs_orb
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm_env)
      CALL get_admm_env(admm_env, matrix_s_aux_fit=matrix_s_aux_fit, sab_aux_fit_asymm=sab_aux_fit_asymm, &
                        matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb, sab_aux_fit_vs_orb=sab_aux_fit_vs_orb)

      my_fval = 2.0_dp
      IF (PRESENT(fval)) my_fval = fval

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
      admm_force = 0.0_dp

      nspins = SIZE(matrix_pz)
      nao = admm_env%nao_orb
      naux = admm_env%nao_aux_fit

      CALL cp_fm_set_all(admm_env%work_aux_orb2, 0.0_dp)

      DO ispin = 1, nspins
         CALL copy_dbcsr_to_fm(matrix_hz(ispin)%matrix, admm_env%work_aux_aux)
         CALL parallel_gemm("N", "T", naux, naux, naux, 1.0_dp, admm_env%s_inv, &
                            admm_env%work_aux_aux, 0.0_dp, admm_env%work_aux_aux2)
         CALL parallel_gemm("N", "N", naux, nao, naux, 1.0_dp, admm_env%work_aux_aux2, &
                            admm_env%A, 0.0_dp, admm_env%work_aux_orb)
         CALL copy_dbcsr_to_fm(matrix_pz(ispin)%matrix, admm_env%work_orb_orb)
         ! admm_env%work_aux_orb2 = S-1*H*A*P
         CALL parallel_gemm("N", "N", naux, nao, nao, 1.0_dp, admm_env%work_aux_orb, &
                            admm_env%work_orb_orb, 1.0_dp, admm_env%work_aux_orb2)
      END DO

      CALL copy_fm_to_dbcsr(admm_env%work_aux_orb2, matrix_w_q, keep_sparsity=.TRUE.)

      ! admm_env%work_aux_aux = S-1*H*A*P*A(T)
      CALL parallel_gemm("N", "T", naux, naux, nao, 1.0_dp, admm_env%work_aux_orb2, &
                         admm_env%A, 0.0_dp, admm_env%work_aux_aux)
      CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, matrix_w_s, keep_sparsity=.TRUE.)

      CALL dbcsr_scale(matrix_w_q, -my_fval)
      CALL dbcsr_scale(matrix_w_s, my_fval)

      CALL build_overlap_force(ks_env, admm_force, &
                               basis_type_a="AUX_FIT", basis_type_b="AUX_FIT", &
                               sab_nl=sab_aux_fit_asymm, matrix_p=matrix_w_s)
      CALL build_overlap_force(ks_env, admm_force, &
                               basis_type_a="AUX_FIT", basis_type_b="ORB", &
                               sab_nl=sab_aux_fit_vs_orb, matrix_p=matrix_w_q)

      ! add forces
      CALL add_qs_force(admm_force, force, "overlap_admm", atomic_kind_set)

      DEALLOCATE (admm_force)
      CALL dbcsr_deallocate_matrix(matrix_w_s)
      CALL dbcsr_deallocate_matrix(matrix_w_q)

      CALL timestop(handle)

   END SUBROUTINE admm_projection_derivative

! **************************************************************************************************
!> \brief Calculates contribution of forces due to basis transformation
!>
!>        dE/dR = dE/dC'*dC'/dR
!>        dE/dC = Ks'*c'*occ = H'
!>
!>        dC'/dR = - tr(A*lambda^(-1/2)*H'^(T)*S^(-1) * dS'/dR)
!>                 - tr(A*C*Y^(T)*C^(T)*Q^(T)*A^(T) * dS'/dR)
!>                 + tr(C*lambda^(-1/2)*H'^(T)*S^(-1) * dQ/dR)
!>                 + tr(A*C*Y^(T)*c^(T) * dQ/dR)
!>                 + tr(C*Y^(T)*C^(T)*A^(T) * dQ/dR)
!>
!>        where
!>
!>        A = S'^(-1)*Q
!>        lambda = C^(T)*B*C
!>        B = Q^(T)*A
!>        Y = R*[ (R^(T)*C^(T)*A^(T)*H'*R) xx M ]*R^(T)
!>        lambda = R*D*R^(T)
!>        Mij = Poles-Matrix (see above)
!>        xx = schur product
!>
!> \param qs_env the QS environment
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE calc_mixed_overlap_force(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_mixed_overlap_force'

      INTEGER                                            :: handle, ispin, iw, nao_aux_fit, nao_orb, &
                                                            natom, neighbor_list_id, nmo
      LOGICAL                                            :: omit_headers
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: admm_force
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_s_aux_fit, &
                                                            matrix_s_aux_fit_vs_orb, rho_ao, &
                                                            rho_ao_aux
      TYPE(dbcsr_type), POINTER                          :: matrix_rho_aux_desymm_tmp, matrix_w_q, &
                                                            matrix_w_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit

      CALL timeset(routineN, handle)

      NULLIFY (admm_env, logger, dft_control, para_env, mos, mo_coeff, matrix_w_q, matrix_w_s, &
               rho, rho_aux_fit, energy, sab_orb, ks_env, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, matrix_s)

      CALL get_qs_env(qs_env, &
                      admm_env=admm_env, &
                      ks_env=ks_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      neighbor_list_id=neighbor_list_id, &
                      rho=rho, &
                      energy=energy, &
                      sab_orb=sab_orb, &
                      mos=mos, &
                      para_env=para_env)
      CALL get_admm_env(admm_env, matrix_s_aux_fit=matrix_s_aux_fit, rho_aux_fit=rho_aux_fit, &
                        matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb)

      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CALL qs_rho_get(rho_aux_fit, &
                      rho_ao=rho_ao_aux)

      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb

      logger => cp_get_default_logger()

      ! *** forces are only implemented for mo_diag or none and basis_projection ***
      IF (admm_env%block_dm) THEN
         CPABORT("ADMM Forces not implemented for blocked projection methods!")
      END IF

      IF (.NOT. (admm_env%purification_method == do_admm_purify_mo_diag .OR. &
                 admm_env%purification_method == do_admm_purify_none)) THEN
         CPABORT("ADMM Forces only implemented without purification or for MO_DIAG.")
      END IF

      ! *** Create sparse work matrices

      ALLOCATE (matrix_w_s)
      CALL dbcsr_create(matrix_w_s, template=matrix_s_aux_fit(1)%matrix, &
                        name='W MATRIX AUX S', &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(matrix_w_s, admm_env%sab_aux_fit_asymm)

      ALLOCATE (matrix_w_q)
      CALL dbcsr_copy(matrix_w_q, matrix_s_aux_fit_vs_orb(1)%matrix, &
                      "W MATRIX AUX Q")

      DO ispin = 1, dft_control%nspins
         nmo = admm_env%nmo(ispin)
         CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff)

         ! *** S'^(-T)*H'
         IF (.NOT. admm_env%purification_method == do_admm_purify_none) THEN
            CALL parallel_gemm('T', 'N', nao_aux_fit, nmo, nao_aux_fit, &
                               1.0_dp, admm_env%S_inv, admm_env%mo_derivs_aux_fit(ispin), 0.0_dp, &
                               admm_env%work_aux_nmo(ispin))
         ELSE

            CALL parallel_gemm('T', 'N', nao_aux_fit, nmo, nao_aux_fit, &
                               1.0_dp, admm_env%S_inv, admm_env%H(ispin), 0.0_dp, &
                               admm_env%work_aux_nmo(ispin))
         END IF

         ! *** S'^(-T)*H'*Lambda^(-T/2)
         CALL parallel_gemm('N', 'T', nao_aux_fit, nmo, nmo, &
                            1.0_dp, admm_env%work_aux_nmo(ispin), admm_env%lambda_inv_sqrt(ispin), 0.0_dp, &
                            admm_env%work_aux_nmo2(ispin))

         ! *** C*Lambda^(-1/2)*H'^(T)*S'^(-1) minus sign due to force = -dE/dR
         CALL parallel_gemm('N', 'T', nao_aux_fit, nao_orb, nmo, &
                            -1.0_dp, admm_env%work_aux_nmo2(ispin), mo_coeff, 0.0_dp, &
                            admm_env%work_aux_orb)

         ! *** A*C*Lambda^(-1/2)*H'^(T)*S'^(-1), minus sign to recover from above
         CALL parallel_gemm('N', 'T', nao_aux_fit, nao_aux_fit, nao_orb, &
                            -1.0_dp, admm_env%work_aux_orb, admm_env%A, 0.0_dp, &
                            admm_env%work_aux_aux)

         IF (.NOT. (admm_env%purification_method == do_admm_purify_none)) THEN
            ! *** C*Y
            CALL parallel_gemm('N', 'N', nao_orb, nmo, nmo, &
                               1.0_dp, mo_coeff, admm_env%R_schur_R_t(ispin), 0.0_dp, &
                               admm_env%work_orb_nmo(ispin))
            ! *** C*Y^(T)*C^(T)
            CALL parallel_gemm('N', 'T', nao_orb, nao_orb, nmo, &
                               1.0_dp, mo_coeff, admm_env%work_orb_nmo(ispin), 0.0_dp, &
                               admm_env%work_orb_orb)
            ! *** A*C*Y^(T)*C^(T) Add to work aux_orb, minus sign due to force = -dE/dR
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb, &
                               -1.0_dp, admm_env%A, admm_env%work_orb_orb, 1.0_dp, &
                               admm_env%work_aux_orb)

            ! *** C*Y^(T)
            CALL parallel_gemm('N', 'T', nao_orb, nmo, nmo, &
                               1.0_dp, mo_coeff, admm_env%R_schur_R_t(ispin), 0.0_dp, &
                               admm_env%work_orb_nmo(ispin))
            ! *** C*Y*C^(T)
            CALL parallel_gemm('N', 'T', nao_orb, nao_orb, nmo, &
                               1.0_dp, mo_coeff, admm_env%work_orb_nmo(ispin), 0.0_dp, &
                               admm_env%work_orb_orb)
            ! *** A*C*Y*C^(T) Add to work aux_orb, minus sign due to -dE/dR
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb, &
                               -1.0_dp, admm_env%A, admm_env%work_orb_orb, 1.0_dp, &
                               admm_env%work_aux_orb)
         END IF

         ! Add derivative contribution matrix*dQ/dR in additional last term in
         ! Eq. (26,32, 33) in Merlot2014 to the force
         ! ADMMS
         IF (admm_env%do_admms) THEN
            ! *** scale admm_env%work_aux_orb by gsi due to inner derivative
            CALL cp_fm_scale(admm_env%gsi(ispin), admm_env%work_aux_orb)
            CALL parallel_gemm('N', 'T', nao_orb, nao_orb, nmo, &
                               4.0_dp*(admm_env%gsi(ispin))*admm_env%lambda_merlot(ispin)/dft_control%nspins, &
                               mo_coeff, mo_coeff, 0.0_dp, admm_env%work_orb_orb2)

            ! *** prefactor*A*C*C^(T) Add to work aux_orb
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb, &
                               1.0_dp, admm_env%A, admm_env%work_orb_orb2, 1.0_dp, &
                               admm_env%work_aux_orb)

            ! ADMMP
         ELSE IF (admm_env%do_admmp) THEN
            CALL cp_fm_scale(admm_env%gsi(ispin)**2, admm_env%work_aux_orb)
            ! *** prefactor*C*C^(T), nspins since 2/n_spin*C*C^(T)=P
            CALL parallel_gemm('N', 'T', nao_orb, nao_orb, nmo, &
                               4.0_dp*(admm_env%gsi(ispin))*admm_env%lambda_merlot(ispin)/dft_control%nspins, &
                               mo_coeff, mo_coeff, 0.0_dp, admm_env%work_orb_orb2)

            ! *** prefactor*A*C*C^(T) Add to work aux_orb
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb, &
                               1.0_dp, admm_env%A, admm_env%work_orb_orb2, 1.0_dp, &
                               admm_env%work_aux_orb)

            ! ADMMQ
         ELSE IF (admm_env%do_admmq) THEN
            ! *** scale admm_env%work_aux_orb by gsi due to inner derivative
            CALL cp_fm_scale(admm_env%gsi(ispin), admm_env%work_aux_orb)
            CALL parallel_gemm('N', 'T', nao_orb, nao_orb, nmo, &
                               4.0_dp*(admm_env%gsi(ispin))*admm_env%lambda_merlot(ispin)/dft_control%nspins, &
                               mo_coeff, mo_coeff, 0.0_dp, admm_env%work_orb_orb2)

            ! *** prefactor*A*C*C^(T) Add to work aux_orb
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb, &
                               1.0_dp, admm_env%A, admm_env%work_orb_orb2, 1.0_dp, &
                               admm_env%work_aux_orb)
         END IF

         ! *** copy to sparse matrix
         CALL copy_fm_to_dbcsr(admm_env%work_aux_orb, matrix_w_q, keep_sparsity=.TRUE.)

         IF (.NOT. (admm_env%purification_method == do_admm_purify_none)) THEN
            ! *** A*C*Y^(T)*C^(T)
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb, &
                               1.0_dp, admm_env%A, admm_env%work_orb_orb, 0.0_dp, &
                               admm_env%work_aux_orb)
            ! *** A*C*Y^(T)*C^(T)*A^(T) add to aux_aux, minus sign cancels
            CALL parallel_gemm('N', 'T', nao_aux_fit, nao_aux_fit, nao_orb, &
                               1.0_dp, admm_env%work_aux_orb, admm_env%A, 1.0_dp, &
                               admm_env%work_aux_aux)
         END IF

         ! *** copy to sparse matrix
         CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, matrix_w_s, keep_sparsity=.TRUE.)

         ! Add derivative of Eq. (33) with respect to s_aux Merlot2014 to the force
         IF (admm_env%do_admmp .OR. admm_env%do_admmq .OR. admm_env%do_admms) THEN

            !Create desymmetrized auxiliary density matrix
            NULLIFY (matrix_rho_aux_desymm_tmp)
            ALLOCATE (matrix_rho_aux_desymm_tmp)
            CALL dbcsr_create(matrix_rho_aux_desymm_tmp, template=matrix_s_aux_fit(1)%matrix, &
                              name='Rho_aux non-symm', &
                              matrix_type=dbcsr_type_no_symmetry)

            CALL dbcsr_desymmetrize(rho_ao_aux(ispin)%matrix, matrix_rho_aux_desymm_tmp)

            ! ADMMS/Q 1. scale original matrix_w_s by gsi due to inner deriv.
            !       2. add derivative of variational term with resp. to s
            IF (admm_env%do_admms .OR. admm_env%do_admmq) THEN
               CALL dbcsr_scale(matrix_w_s, admm_env%gsi(ispin))
               CALL dbcsr_add(matrix_w_s, matrix_rho_aux_desymm_tmp, 1.0_dp, &
                              -admm_env%lambda_merlot(ispin))

               ! ADMMP scale by gsi^2 and add derivative of variational term with resp. to s
            ELSE IF (admm_env%do_admmp) THEN

               CALL dbcsr_scale(matrix_w_s, admm_env%gsi(ispin)**2)
               CALL dbcsr_add(matrix_w_s, matrix_rho_aux_desymm_tmp, 1.0_dp, &
                              (-admm_env%gsi(ispin))*admm_env%lambda_merlot(ispin))

            END IF

            CALL dbcsr_deallocate_matrix(matrix_rho_aux_desymm_tmp)

         END IF

         ! allocate force vector
         CALL get_qs_env(qs_env=qs_env, natom=natom)
         ALLOCATE (admm_force(3, natom))
         admm_force = 0.0_dp
         CALL build_overlap_force(ks_env, admm_force, &
                                  basis_type_a="AUX_FIT", basis_type_b="AUX_FIT", &
                                  sab_nl=admm_env%sab_aux_fit_asymm, matrix_p=matrix_w_s)
         CALL build_overlap_force(ks_env, admm_force, &
                                  basis_type_a="AUX_FIT", basis_type_b="ORB", &
                                  sab_nl=admm_env%sab_aux_fit_vs_orb, matrix_p=matrix_w_q)

         ! Add contribution of original basis set for ADMMQ, P, S
         IF (admm_env%do_admmq .OR. admm_env%do_admmp .OR. admm_env%do_admms) THEN
            CALL dbcsr_scale(rho_ao(ispin)%matrix, -admm_env%lambda_merlot(ispin))
            CALL build_overlap_force(ks_env, admm_force, &
                                     basis_type_a="ORB", basis_type_b="ORB", &
                                     sab_nl=sab_orb, matrix_p=rho_ao(ispin)%matrix)
            CALL dbcsr_scale(rho_ao(ispin)%matrix, -1.0_dp/admm_env%lambda_merlot(ispin))
         END IF

         ! add forces
         CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                         force=force)
         CALL add_qs_force(admm_force, force, "overlap_admm", atomic_kind_set)
         DEALLOCATE (admm_force)

         CALL section_vals_val_get(qs_env%input, "DFT%PRINT%AO_MATRICES%OMIT_HEADERS", l_val=omit_headers)
         IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                              qs_env%input, "DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT"), cp_p_file)) THEN
            iw = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT", &
                                      extension=".Log")
            CALL cp_dbcsr_write_sparse_matrix(matrix_w_s, 4, 6, qs_env, &
                                              para_env, output_unit=iw, omit_headers=omit_headers)
            CALL cp_print_key_finished_output(iw, logger, qs_env%input, &
                                              "DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT")
         END IF
         IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                              qs_env%input, "DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT"), cp_p_file)) THEN
            iw = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT", &
                                      extension=".Log")
            CALL cp_dbcsr_write_sparse_matrix(matrix_w_q, 4, 6, qs_env, &
                                              para_env, output_unit=iw, omit_headers=omit_headers)
            CALL cp_print_key_finished_output(iw, logger, qs_env%input, &
                                              "DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT")
         END IF

      END DO !spin loop

      ! *** Deallocated weighted density matrices
      CALL dbcsr_deallocate_matrix(matrix_w_s)
      CALL dbcsr_deallocate_matrix(matrix_w_q)

      CALL timestop(handle)

   END SUBROUTINE calc_mixed_overlap_force

! **************************************************************************************************
!> \brief ...
!> \param admm_env environment of auxiliary DM
!> \param mo_set ...
!> \param density_matrix auxiliary DM
!> \param overlap_matrix auxiliary OM
!> \param density_matrix_large DM of the original basis
!> \param overlap_matrix_large overlap matrix of original basis
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE calculate_dm_mo_no_diag(admm_env, mo_set, density_matrix, overlap_matrix, &
                                      density_matrix_large, overlap_matrix_large, ispin)
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(mo_set_type), INTENT(IN)                      :: mo_set
      TYPE(dbcsr_type), POINTER                          :: density_matrix, overlap_matrix, &
                                                            density_matrix_large, &
                                                            overlap_matrix_large
      INTEGER                                            :: ispin

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_dm_mo_no_diag'

      INTEGER                                            :: handle, nao_aux_fit, nmo
      REAL(KIND=dp)                                      :: alpha, nel_tmp_aux

! Number of electrons in the aux. DM

      CALL timeset(routineN, handle)

      CALL dbcsr_set(density_matrix, 0.0_dp)
      nao_aux_fit = admm_env%nao_aux_fit
      nmo = admm_env%nmo(ispin)
      CALL cp_fm_to_fm(admm_env%C_hat(ispin), admm_env%work_aux_nmo(ispin))
      CALL cp_fm_column_scale(admm_env%work_aux_nmo(ispin), mo_set%occupation_numbers(1:mo_set%homo))

      CALL parallel_gemm('N', 'N', nao_aux_fit, nmo, nmo, &
                         1.0_dp, admm_env%work_aux_nmo(ispin), admm_env%lambda_inv(ispin), 0.0_dp, &
                         admm_env%work_aux_nmo2(ispin))

      ! The following IF doesn't do anything unless !alpha=mo_set%maxocc is uncommented.
      IF (.NOT. mo_set%uniform_occupation) THEN ! not all orbitals 1..homo are equally occupied
         alpha = 1.0_dp
         CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix, &
                                    matrix_v=admm_env%C_hat(ispin), &
                                    matrix_g=admm_env%work_aux_nmo2(ispin), &
                                    ncol=mo_set%homo, &
                                    alpha=alpha)
      ELSE
         alpha = 1.0_dp
         !alpha=mo_set%maxocc
         CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix, &
                                    matrix_v=admm_env%C_hat(ispin), &
                                    matrix_g=admm_env%work_aux_nmo2(ispin), &
                                    ncol=mo_set%homo, &
                                    alpha=alpha)
      END IF

      !  The following IF checks whether gsi needs to be calculated. This is the case if
      !   the auxiliary density matrix gets scaled
      !   according to Eq. 22 (Merlot) or a scaling of exchange_correction is employed, Eq. 35 (Merlot).
      IF (admm_env%do_admmp .OR. admm_env%do_admmq .OR. admm_env%do_admms) THEN

         CALL cite_reference(Merlot2014)

         admm_env%n_large_basis(3) = 0.0_dp

         ! Calculate number of electrons in the original density matrix, transposing doesn't matter
         ! since both matrices are symmetric
         CALL dbcsr_dot(density_matrix_large, overlap_matrix_large, admm_env%n_large_basis(ispin))
         admm_env%n_large_basis(3) = admm_env%n_large_basis(3) + admm_env%n_large_basis(ispin)
         ! Calculate number of electrons in the auxiliary density matrix
         CALL dbcsr_dot(density_matrix, overlap_matrix, nel_tmp_aux)
         admm_env%gsi(ispin) = admm_env%n_large_basis(ispin)/nel_tmp_aux

         IF (admm_env%do_admmq .OR. admm_env%do_admms) THEN
            ! multiply aux. DM with gsi to get the scaled DM (Merlot, Eq. 21)
            CALL dbcsr_scale(density_matrix, admm_env%gsi(ispin))
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_dm_mo_no_diag

! **************************************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param density_matrix ...
!> \param density_matrix_aux ...
!> \param ispin ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE blockify_density_matrix(admm_env, density_matrix, density_matrix_aux, &
                                      ispin, nspins)
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_type), POINTER                          :: density_matrix, density_matrix_aux
      INTEGER                                            :: ispin, nspins

      CHARACTER(len=*), PARAMETER :: routineN = 'blockify_density_matrix'

      INTEGER                                            :: blk, handle, iatom, jatom
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: sparse_block, sparse_block_aux
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      ! ** set blocked density matrix to 0
      CALL dbcsr_set(density_matrix_aux, 0.0_dp)

      ! ** now loop through the list and copy corresponding blocks
      CALL dbcsr_iterator_start(iter, density_matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
         IF (admm_env%block_map(iatom, jatom) == 1) THEN
            CALL dbcsr_get_block_p(density_matrix_aux, &
                                   row=iatom, col=jatom, BLOCK=sparse_block_aux, found=found)
            IF (found) THEN
               sparse_block_aux = sparse_block
            END IF

         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL copy_dbcsr_to_fm(density_matrix_aux, admm_env%P_to_be_purified(ispin))
      CALL cp_fm_upper_to_full(admm_env%P_to_be_purified(ispin), admm_env%work_orb_orb2)

      IF (nspins == 1) THEN
         CALL cp_fm_scale(0.5_dp, admm_env%P_to_be_purified(ispin))
      END IF

      CALL timestop(handle)
   END SUBROUTINE blockify_density_matrix

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION delta(x)
      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: delta

      IF (x == 0.0_dp) THEN !TODO: exact comparison of reals?
         delta = 1.0_dp
      ELSE
         delta = 0.0_dp
      END IF

   END FUNCTION delta

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION Heaviside(x)
      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: Heaviside

      IF (x < 0.0_dp) THEN
         Heaviside = 0.0_dp
      ELSE
         Heaviside = 1.0_dp
      END IF
   END FUNCTION Heaviside

! **************************************************************************************************
!> \brief Calculate ADMM auxiliary response density
!> \param qs_env ...
!> \param dm ...
!> \param dm_admm ...
! **************************************************************************************************
   SUBROUTINE admm_aux_response_density(qs_env, dm, dm_admm)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: dm
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: dm_admm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_aux_response_density'

      INTEGER                                            :: handle, ispin, nao, nao_aux, ncol, nspins
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, admm_env=admm_env, dft_control=dft_control)

      nspins = dft_control%nspins

      CPASSERT(ASSOCIATED(admm_env%A))
      CPASSERT(ASSOCIATED(admm_env%work_orb_orb))
      CPASSERT(ASSOCIATED(admm_env%work_aux_orb))
      CPASSERT(ASSOCIATED(admm_env%work_aux_aux))
      CALL cp_fm_get_info(admm_env%A, nrow_global=nao_aux, ncol_global=nao)

      ! P1 -> AUX BASIS
      CALL cp_fm_get_info(admm_env%work_orb_orb, nrow_global=nao, ncol_global=ncol)
      DO ispin = 1, nspins
         CALL copy_dbcsr_to_fm(dm(ispin)%matrix, admm_env%work_orb_orb)
         CALL parallel_gemm('N', 'N', nao_aux, ncol, nao, 1.0_dp, admm_env%A, &
                            admm_env%work_orb_orb, 0.0_dp, admm_env%work_aux_orb)
         CALL parallel_gemm('N', 'T', nao_aux, nao_aux, nao, 1.0_dp, admm_env%A, &
                            admm_env%work_aux_orb, 0.0_dp, admm_env%work_aux_aux)
         CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, dm_admm(ispin)%matrix, keep_sparsity=.TRUE.)
      END DO

      CALL timestop(handle)

   END SUBROUTINE admm_aux_response_density

! **************************************************************************************************
!> \brief Fill the ADMM overlp and basis change  matrices in the KP env based on the real-space array
!> \param qs_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE kpoint_calc_admm_matrices(qs_env, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: calculate_forces

      INTEGER                                            :: ic, igroup, ik, ikp, indx, kplocal, &
                                                            nao_aux_fit, nao_orb, nc, nkp, &
                                                            nkp_groups
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: my_kpgrp, use_real_wfn
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(copy_info_type), ALLOCATABLE, DIMENSION(:, :) :: info
      TYPE(cp_cfm_type)                                  :: cmat_aux_fit, cmat_aux_fit_vs_orb, &
                                                            cwork_aux_fit, cwork_aux_fit_vs_orb
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct_aux_fit, &
                                                            matrix_struct_aux_fit_vs_orb
      TYPE(cp_fm_type)                                   :: fmdummy, imat_aux_fit, &
                                                            imat_aux_fit_vs_orb, rmat_aux_fit, &
                                                            rmat_aux_fit_vs_orb, work_aux_fit
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fmwork
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_aux_fit, matrix_s_aux_fit_vs_orb
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: dbcsr_aux_fit, dbcsr_aux_fit_vs_orb
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env_global, para_env_local
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux_fit, sab_aux_fit_vs_orb

      NULLIFY (xkp, kp_dist, para_env_local, cell_to_index, admm_env, kp, &
               kpoints, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, sab_aux_fit, sab_aux_fit_vs_orb, &
               para_env_global, matrix_struct_aux_fit, matrix_struct_aux_fit_vs_orb)

      CALL get_qs_env(qs_env, kpoints=kpoints, admm_env=admm_env)

      CALL get_admm_env(admm_env, matrix_s_aux_fit_kp=matrix_s_aux_fit, &
                        matrix_s_aux_fit_vs_orb_kp=matrix_s_aux_fit_vs_orb, &
                        sab_aux_fit=sab_aux_fit, &
                        sab_aux_fit_vs_orb=sab_aux_fit_vs_orb)

      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, use_real_wfn=use_real_wfn, kp_range=kp_range, &
                           nkp_groups=nkp_groups, kp_dist=kp_dist, cell_to_index=cell_to_index)
      kplocal = kp_range(2) - kp_range(1) + 1
      nc = 1
      IF (.NOT. use_real_wfn) nc = 2

      ALLOCATE (dbcsr_aux_fit(3))
      CALL dbcsr_create(dbcsr_aux_fit(1), template=matrix_s_aux_fit(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(dbcsr_aux_fit(2), template=matrix_s_aux_fit(1, 1)%matrix, matrix_type=dbcsr_type_antisymmetric)
      CALL dbcsr_create(dbcsr_aux_fit(3), template=matrix_s_aux_fit(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(dbcsr_aux_fit(1), sab_aux_fit)
      CALL cp_dbcsr_alloc_block_from_nbl(dbcsr_aux_fit(2), sab_aux_fit)

      ALLOCATE (dbcsr_aux_fit_vs_orb(2))
      CALL dbcsr_create(dbcsr_aux_fit_vs_orb(1), template=matrix_s_aux_fit_vs_orb(1, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(dbcsr_aux_fit_vs_orb(2), template=matrix_s_aux_fit_vs_orb(1, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(dbcsr_aux_fit_vs_orb(1), sab_aux_fit_vs_orb)
      CALL cp_dbcsr_alloc_block_from_nbl(dbcsr_aux_fit_vs_orb(2), sab_aux_fit_vs_orb)

      !Create global work fm
      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb
      para_env_global => kpoints%blacs_env_all%para_env

      ALLOCATE (fmwork(4))
      CALL cp_fm_struct_create(matrix_struct_aux_fit, context=kpoints%blacs_env_all, para_env=para_env_global, &
                               nrow_global=nao_aux_fit, ncol_global=nao_aux_fit)
      CALL cp_fm_create(fmwork(1), matrix_struct_aux_fit)
      CALL cp_fm_create(fmwork(2), matrix_struct_aux_fit)
      CALL cp_fm_struct_release(matrix_struct_aux_fit)

      CALL cp_fm_struct_create(matrix_struct_aux_fit_vs_orb, context=kpoints%blacs_env_all, para_env=para_env_global, &
                               nrow_global=nao_aux_fit, ncol_global=nao_orb)
      CALL cp_fm_create(fmwork(3), matrix_struct_aux_fit_vs_orb)
      CALL cp_fm_create(fmwork(4), matrix_struct_aux_fit_vs_orb)
      CALL cp_fm_struct_release(matrix_struct_aux_fit_vs_orb)

      !Create fm local to the KP groups
      nao_aux_fit = admm_env%nao_aux_fit
      nao_orb = admm_env%nao_orb
      para_env_local => kpoints%blacs_env%para_env

      CALL cp_fm_struct_create(matrix_struct_aux_fit, context=kpoints%blacs_env, para_env=para_env_local, &
                               nrow_global=nao_aux_fit, ncol_global=nao_aux_fit)
      CALL cp_fm_create(rmat_aux_fit, matrix_struct_aux_fit)
      CALL cp_fm_create(imat_aux_fit, matrix_struct_aux_fit)
      CALL cp_fm_create(work_aux_fit, matrix_struct_aux_fit)
      CALL cp_cfm_create(cwork_aux_fit, matrix_struct_aux_fit)
      CALL cp_cfm_create(cmat_aux_fit, matrix_struct_aux_fit)

      CALL cp_fm_struct_create(matrix_struct_aux_fit_vs_orb, context=kpoints%blacs_env, para_env=para_env_local, &
                               nrow_global=nao_aux_fit, ncol_global=nao_orb)
      CALL cp_fm_create(rmat_aux_fit_vs_orb, matrix_struct_aux_fit_vs_orb)
      CALL cp_fm_create(imat_aux_fit_vs_orb, matrix_struct_aux_fit_vs_orb)
      CALL cp_cfm_create(cwork_aux_fit_vs_orb, matrix_struct_aux_fit_vs_orb)
      CALL cp_cfm_create(cmat_aux_fit_vs_orb, matrix_struct_aux_fit_vs_orb)

      ALLOCATE (info(kplocal*nkp_groups, 4))

      ! Steup and start all the communication
      indx = 0
      DO ikp = 1, kplocal
         DO igroup = 1, nkp_groups
            ik = kp_dist(1, igroup) + ikp - 1
            my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
            indx = indx + 1

            IF (use_real_wfn) THEN
               !AUX-AUX overlap
               CALL dbcsr_set(dbcsr_aux_fit(1), 0.0_dp)
               CALL rskp_transform(rmatrix=dbcsr_aux_fit(1), rsmat=matrix_s_aux_fit, ispin=1, &
                                   xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_aux_fit)
               CALL dbcsr_desymmetrize(dbcsr_aux_fit(1), dbcsr_aux_fit(3))
               CALL copy_dbcsr_to_fm(dbcsr_aux_fit(3), fmwork(1))

               !AUX-ORB overlap
               CALL dbcsr_set(dbcsr_aux_fit_vs_orb(1), 0.0_dp)
               CALL rskp_transform(rmatrix=dbcsr_aux_fit_vs_orb(1), rsmat=matrix_s_aux_fit_vs_orb, ispin=1, &
                                   xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_aux_fit_vs_orb)
               CALL copy_dbcsr_to_fm(dbcsr_aux_fit_vs_orb(1), fmwork(3))
            ELSE
               !AUX-AUX overlap
               CALL dbcsr_set(dbcsr_aux_fit(1), 0.0_dp)
               CALL dbcsr_set(dbcsr_aux_fit(2), 0.0_dp)
               CALL rskp_transform(rmatrix=dbcsr_aux_fit(1), cmatrix=dbcsr_aux_fit(2), rsmat=matrix_s_aux_fit, &
                                   ispin=1, xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_aux_fit)
               CALL dbcsr_desymmetrize(dbcsr_aux_fit(1), dbcsr_aux_fit(3))
               CALL copy_dbcsr_to_fm(dbcsr_aux_fit(3), fmwork(1))
               CALL dbcsr_desymmetrize(dbcsr_aux_fit(2), dbcsr_aux_fit(3))
               CALL copy_dbcsr_to_fm(dbcsr_aux_fit(3), fmwork(2))

               !AUX-ORB overlap
               CALL dbcsr_set(dbcsr_aux_fit_vs_orb(1), 0.0_dp)
               CALL dbcsr_set(dbcsr_aux_fit_vs_orb(2), 0.0_dp)
               CALL rskp_transform(rmatrix=dbcsr_aux_fit_vs_orb(1), cmatrix=dbcsr_aux_fit_vs_orb(2), &
                                   rsmat=matrix_s_aux_fit_vs_orb, ispin=1, xkp=xkp(1:3, ik), &
                                   cell_to_index=cell_to_index, sab_nl=sab_aux_fit_vs_orb)
               CALL copy_dbcsr_to_fm(dbcsr_aux_fit_vs_orb(1), fmwork(3))
               CALL copy_dbcsr_to_fm(dbcsr_aux_fit_vs_orb(2), fmwork(4))
            END IF

            IF (my_kpgrp) THEN
               CALL cp_fm_start_copy_general(fmwork(1), rmat_aux_fit, para_env_global, info(indx, 1))
               CALL cp_fm_start_copy_general(fmwork(3), rmat_aux_fit_vs_orb, para_env_global, info(indx, 3))
               IF (.NOT. use_real_wfn) THEN
                  CALL cp_fm_start_copy_general(fmwork(2), imat_aux_fit, para_env_global, info(indx, 2))
                  CALL cp_fm_start_copy_general(fmwork(4), imat_aux_fit_vs_orb, para_env_global, info(indx, 4))
               END IF
            ELSE
               CALL cp_fm_start_copy_general(fmwork(1), fmdummy, para_env_global, info(indx, 1))
               CALL cp_fm_start_copy_general(fmwork(3), fmdummy, para_env_global, info(indx, 3))
               IF (.NOT. use_real_wfn) THEN
                  CALL cp_fm_start_copy_general(fmwork(2), fmdummy, para_env_global, info(indx, 2))
                  CALL cp_fm_start_copy_general(fmwork(4), fmdummy, para_env_global, info(indx, 4))
               END IF
            END IF

         END DO
      END DO

      ! Finish communication and store
      indx = 0
      DO ikp = 1, kplocal
         DO igroup = 1, nkp_groups
            ik = kp_dist(1, igroup) + ikp - 1
            my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
            indx = indx + 1

            IF (my_kpgrp) THEN
               CALL cp_fm_finish_copy_general(rmat_aux_fit, info(indx, 1))
               CALL cp_fm_finish_copy_general(rmat_aux_fit_vs_orb, info(indx, 3))
               IF (.NOT. use_real_wfn) THEN
                  CALL cp_fm_finish_copy_general(imat_aux_fit, info(indx, 2))
                  CALL cp_fm_finish_copy_general(imat_aux_fit_vs_orb, info(indx, 4))
               END IF
            END IF
         END DO

         kp => kpoints%kp_aux_env(ikp)%kpoint_env

         !Allocate local KP matrices
         CALL cp_fm_release(kp%amat)
         ALLOCATE (kp%amat(nc, 1))
         DO ic = 1, nc
            CALL cp_fm_create(kp%amat(ic, 1), matrix_struct_aux_fit_vs_orb)
         END DO

         !Only need the overlap in case of ADMMP, ADMMQ or ADMMS, or for forces
         IF (admm_env%do_admmp .OR. admm_env%do_admmq .OR. admm_env%do_admms .OR. calculate_forces) THEN
            CALL cp_fm_release(kp%smat)
            ALLOCATE (kp%smat(nc, 1))
            DO ic = 1, nc
               CALL cp_fm_create(kp%smat(ic, 1), matrix_struct_aux_fit)
            END DO
            CALL cp_fm_to_fm(rmat_aux_fit, kp%smat(1, 1))
            IF (.NOT. use_real_wfn) CALL cp_fm_to_fm(imat_aux_fit, kp%smat(2, 1))
         END IF

         IF (use_real_wfn) THEN
            !Invert S_aux
            CALL cp_fm_cholesky_decompose(rmat_aux_fit)
            CALL cp_fm_cholesky_invert(rmat_aux_fit)
            CALL cp_fm_upper_to_full(rmat_aux_fit, work_aux_fit)

            !A = S^-1 * Q
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, 1.0_dp, &
                               rmat_aux_fit, rmat_aux_fit_vs_orb, 0.0_dp, kp%amat(1, 1))
         ELSE

            !Invert S_aux
            CALL cp_fm_to_cfm(rmat_aux_fit, imat_aux_fit, cmat_aux_fit)
            CALL cp_cfm_cholesky_decompose(cmat_aux_fit)
            CALL cp_cfm_cholesky_invert(cmat_aux_fit)
            CALL own_cfm_upper_to_full(cmat_aux_fit, cwork_aux_fit)

            !A = S^-1 * Q
            CALL cp_fm_to_cfm(rmat_aux_fit_vs_orb, imat_aux_fit_vs_orb, cmat_aux_fit_vs_orb)
            CALL parallel_gemm('N', 'N', nao_aux_fit, nao_orb, nao_aux_fit, z_one, &
                               cmat_aux_fit, cmat_aux_fit_vs_orb, z_zero, cwork_aux_fit_vs_orb)
            CALL cp_cfm_to_fm(cwork_aux_fit_vs_orb, kp%amat(1, 1), kp%amat(2, 1))
         END IF
      END DO

      ! Clean up communication
      indx = 0
      DO ikp = 1, kplocal
         DO igroup = 1, nkp_groups
            ik = kp_dist(1, igroup) + ikp - 1
            my_kpgrp = (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2))
            indx = indx + 1

            IF (my_kpgrp) THEN
               CALL cp_fm_cleanup_copy_general(info(indx, 1))
               CALL cp_fm_cleanup_copy_general(info(indx, 3))
               IF (.NOT. use_real_wfn) THEN
                  CALL cp_fm_cleanup_copy_general(info(indx, 2))
                  CALL cp_fm_cleanup_copy_general(info(indx, 4))
               END IF
            END IF

         END DO
      END DO

      CALL cp_fm_release(rmat_aux_fit)
      CALL cp_fm_release(imat_aux_fit)
      CALL cp_fm_release(work_aux_fit)
      CALL cp_cfm_release(cwork_aux_fit)
      CALL cp_cfm_release(cmat_aux_fit)
      CALL cp_fm_release(rmat_aux_fit_vs_orb)
      CALL cp_fm_release(imat_aux_fit_vs_orb)
      CALL cp_cfm_release(cwork_aux_fit_vs_orb)
      CALL cp_cfm_release(cmat_aux_fit_vs_orb)
      CALL cp_fm_struct_release(matrix_struct_aux_fit)
      CALL cp_fm_struct_release(matrix_struct_aux_fit_vs_orb)

      CALL cp_fm_release(fmwork(1))
      CALL cp_fm_release(fmwork(2))
      CALL cp_fm_release(fmwork(3))
      CALL cp_fm_release(fmwork(4))

      CALL dbcsr_release(dbcsr_aux_fit(1))
      CALL dbcsr_release(dbcsr_aux_fit(2))
      CALL dbcsr_release(dbcsr_aux_fit(3))
      CALL dbcsr_release(dbcsr_aux_fit_vs_orb(1))
      CALL dbcsr_release(dbcsr_aux_fit_vs_orb(2))

   END SUBROUTINE kpoint_calc_admm_matrices

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_Q ...
!> \param cfm_mat_work ...
! **************************************************************************************************
   SUBROUTINE own_cfm_upper_to_full(cfm_mat_Q, cfm_mat_work)

      TYPE(cp_cfm_type), INTENT(IN)                      :: cfm_mat_Q, cfm_mat_work

      CHARACTER(LEN=*), PARAMETER :: routineN = 'own_cfm_upper_to_full'

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      !TODO: routine cpoied from rpa_gw_kpoints_util.F => put it somewhere centralized?

      ! get info of fm_mat_Q
      CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global < i_global) THEN
               cfm_mat_Q%local_data(iiB, jjB) = z_zero
            END IF
            IF (j_global == i_global) THEN
               cfm_mat_Q%local_data(iiB, jjB) = cfm_mat_Q%local_data(iiB, jjB)/(2.0_dp, 0.0_dp)
            END IF
         END DO
      END DO

      CALL cp_cfm_transpose(cfm_mat_Q, 'C', cfm_mat_work)

      CALL cp_cfm_scale_and_add(z_one, cfm_mat_Q, z_one, cfm_mat_work)

      CALL timestop(handle)

   END SUBROUTINE

END MODULE admm_methods
