!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate and distribute 2c- and 3c- integrals for RI
!> \par History
!>      06.2012 created [Mauro Del Ben]
!>      03.2019 separated from mp2_ri_gpw [Frederick Stein]
! **************************************************************************************************
MODULE mp2_integrals
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE bibliography,                    ONLY: DelBen2013,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_m_by_n_from_template
   USE cp_eri_mme_interface,            ONLY: cp_eri_mme_param,&
                                              cp_eri_mme_set_params
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_para_env,                     ONLY: cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_get_info, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_release_p, dbcsr_scalar, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, &
        dbcsr_type_real_8
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_clear, dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_create, dbcsr_t_destroy, &
        dbcsr_t_distribution_destroy, dbcsr_t_distribution_new, dbcsr_t_distribution_type, &
        dbcsr_t_filter, dbcsr_t_get_block, dbcsr_t_get_info, dbcsr_t_get_num_blocks, &
        dbcsr_t_get_stored_coordinates, dbcsr_t_mp_environ_pgrid, dbcsr_t_pgrid_create, &
        dbcsr_t_pgrid_destroy, dbcsr_t_pgrid_type, dbcsr_t_put_block, dbcsr_t_reserve_blocks, &
        dbcsr_t_reserved_block_indices, dbcsr_t_split_blocks, dbcsr_t_type
   USE group_dist_types,                ONLY: create_group_dist,&
                                              get_group_dist,&
                                              group_dist_d1_type
   USE hfx_types,                       ONLY: alloc_containers,&
                                              block_ind_type,&
                                              hfx_compression_type
   USE input_constants,                 ONLY: do_eri_gpw,&
                                              do_eri_mme,&
                                              do_eri_os,&
                                              do_potential_coulomb,&
                                              do_potential_id,&
                                              do_potential_long,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: kpoint_type
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE machine,                         ONLY: m_flush
   USE message_passing,                 ONLY: mp_cart_create,&
                                              mp_environ,&
                                              mp_max,&
                                              mp_sendrecv,&
                                              mp_sum,&
                                              mp_sync
   USE mp2_eri,                         ONLY: mp2_eri_3c_integrate
   USE mp2_eri_gpw,                     ONLY: cleanup_gpw,&
                                              mp2_eri_3c_integrate_gpw,&
                                              prepare_gpw
   USE mp2_ri_2c,                       ONLY: get_2c_integrals
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_tensors,                      ONLY: build_3c_integrals,&
                                              build_3c_neighbor_lists,&
                                              compress_tensor,&
                                              get_tensor_occupancy,&
                                              neighbor_list_3c_destroy
   USE qs_tensors_types,                ONLY: create_3c_tensor,&
                                              create_tensor_batches,&
                                              distribution_3d_create,&
                                              distribution_3d_type,&
                                              neighbor_list_3c_type,&
                                              pgf_block_sizes
   USE task_list_types,                 ONLY: task_list_type
   USE util,                            ONLY: get_limit

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_integrals'

   PUBLIC :: mp2_ri_gpw_compute_in

CONTAINS

! **************************************************************************************************
!> \brief with ri mp2 gpw
!> \param BIb_C ...
!> \param BIb_C_gw ...
!> \param BIb_C_bse_ij ...
!> \param BIb_C_bse_ab ...
!> \param gd_array ...
!> \param gd_B_virtual ...
!> \param dimen_RI ...
!> \param dimen_RI_red ...
!> \param qs_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param color_sub ...
!> \param cell ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param mo_coeff ...
!> \param fm_matrix_L_RI_metric ...
!> \param nmo ...
!> \param homo ...
!> \param mat_munu ...
!> \param sab_orb_sub ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param mo_coeff_all ...
!> \param mo_coeff_gw ...
!> \param eps_filter ...
!> \param unit_nr ...
!> \param mp2_memory ...
!> \param calc_PQ_cond_num ...
!> \param calc_forces ...
!> \param blacs_env_sub ...
!> \param my_do_gw ...
!> \param do_bse ...
!> \param gd_B_all ...
!> \param starts_array_mc ...
!> \param ends_array_mc ...
!> \param starts_array_mc_block ...
!> \param ends_array_mc_block ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param do_im_time ...
!> \param do_kpoints_cubic_RPA ...
!> \param kpoints ...
!> \param t_3c_M ...
!> \param t_3c_O ...
!> \param t_3c_O_compressed ...
!> \param t_3c_O_ind ...
!> \param ri_metric ...
!> \param gd_B_occ_bse ...
!> \param gd_B_virt_bse ...
!> \param BIb_C_beta ...
!> \param BIb_C_gw_beta ...
!> \param gd_B_virtual_beta ...
!> \param homo_beta ...
!> \param mo_coeff_o_beta ...
!> \param mo_coeff_v_beta ...
!> \param mo_coeff_all_beta ...
!> \param mo_coeff_gw_beta ...
!> \author Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE mp2_ri_gpw_compute_in(BIb_C, BIb_C_gw, BIb_C_bse_ij, BIb_C_bse_ab, gd_array, gd_B_virtual, &
                                    dimen_RI, dimen_RI_red, qs_env, para_env, para_env_sub, color_sub, &
                                    cell, particle_set, atomic_kind_set, qs_kind_set, mo_coeff, &
                                    fm_matrix_L_RI_metric, nmo, homo, &
                                    mat_munu, &
                                    sab_orb_sub, mo_coeff_o, mo_coeff_v, mo_coeff_all, &
                                    mo_coeff_gw, eps_filter, unit_nr, &
                                    mp2_memory, calc_PQ_cond_num, calc_forces, blacs_env_sub, my_do_gw, do_bse, &
                                    gd_B_all, starts_array_mc, ends_array_mc, &
                                    starts_array_mc_block, ends_array_mc_block, &
                                    gw_corr_lev_occ, gw_corr_lev_virt, &
                                    do_im_time, do_kpoints_cubic_RPA, kpoints, &
                                    t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                                    ri_metric, gd_B_occ_bse, gd_B_virt_bse, BIb_C_beta, BIb_C_gw_beta, &
                                    gd_B_virtual_beta, homo_beta, mo_coeff_o_beta, mo_coeff_v_beta, &
                                    mo_coeff_all_beta, mo_coeff_gw_beta)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: BIb_C, BIb_C_gw, BIb_C_bse_ij, &
                                                            BIb_C_bse_ab
      TYPE(group_dist_d1_type), INTENT(OUT)              :: gd_array, gd_B_virtual
      INTEGER, INTENT(OUT)                               :: dimen_RI, dimen_RI_red
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, INTENT(IN)                                :: color_sub
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_RI_metric
      INTEGER, INTENT(IN)                                :: nmo, homo
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_munu
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: sab_orb_sub
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_o, mo_coeff_v, mo_coeff_all, &
                                                            mo_coeff_gw
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), INTENT(IN)                          :: mp2_memory
      LOGICAL, INTENT(IN)                                :: calc_PQ_cond_num, calc_forces
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      LOGICAL, INTENT(IN)                                :: my_do_gw, do_bse
      TYPE(group_dist_d1_type), INTENT(OUT)              :: gd_B_all
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: starts_array_mc, ends_array_mc, &
                                                            starts_array_mc_block, &
                                                            ends_array_mc_block
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt
      LOGICAL, INTENT(IN)                                :: do_im_time, do_kpoints_cubic_RPA
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_t_type), INTENT(OUT)                    :: t_3c_M
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: t_3c_O
      TYPE(hfx_compression_type), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: t_3c_O_compressed
      TYPE(block_ind_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: t_3c_O_ind
      TYPE(libint_potential_type), INTENT(IN)            :: ri_metric
      TYPE(group_dist_d1_type), INTENT(OUT)              :: gd_B_occ_bse, gd_B_virt_bse
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT), OPTIONAL                           :: BIb_C_beta, BIb_C_gw_beta
      TYPE(group_dist_d1_type), INTENT(OUT), OPTIONAL    :: gd_B_virtual_beta
      INTEGER, INTENT(IN), OPTIONAL                      :: homo_beta
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: mo_coeff_o_beta, mo_coeff_v_beta, &
                                                            mo_coeff_all_beta, mo_coeff_gw_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_gpw_compute_in'

      INTEGER :: cm, cut_memory, cut_memory_int, eri_method, gw_corr_lev_total, handle, handle2, &
         handle4, i, i_counter, i_mem, itmp(2), j, jcell, kcell, LLL, max_row_col_local, &
         max_row_col_local_beta, max_row_col_local_gw, max_row_col_local_occ_bse, &
         max_row_col_local_virt_bse, min_bsize, mp_comm_t3c_2, my_B_all_end, my_B_all_size, &
         my_B_all_start, my_B_occ_bse_end, my_B_occ_bse_size, my_B_occ_bse_start, my_B_size, &
         my_B_size_beta, my_B_virt_bse_end, my_B_virt_bse_size, my_B_virt_bse_start, &
         my_B_virtual_end, my_B_virtual_end_beta, my_B_virtual_start, my_B_virtual_start_beta, &
         my_group_L_end
      INTEGER :: my_group_L_size, my_group_L_start, natom, ngroup, nimg, nkind, num_small_eigen, &
         potential_type, ri_metric_type, virtual, virtual_beta
      INTEGER(int_8)                                     :: nze
      INTEGER, ALLOCATABLE, DIMENSION(:) :: dist_AO_1, dist_AO_2, dist_RI, &
         ends_array_mc_block_int, ends_array_mc_int, sizes_AO, sizes_AO_split, sizes_RI, &
         sizes_RI_split, starts_array_mc_block_int, starts_array_mc_int, sub_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: local_col_row_info, local_col_row_info_beta, &
         local_col_row_info_gw, local_col_row_info_occ_bse, local_col_row_info_virt_bse
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(3)                              :: bounds_3c, pcoord, pdims, pdims_t3c, &
                                                            periodic
      LOGICAL                                            :: do_alpha_beta, do_gpw, &
                                                            do_kpoints_from_Gamma, do_svd, &
                                                            memory_info
      REAL(KIND=dp) :: compression_factor, cond_num, cutoff_old, eps_svd, mem_for_iaK, memory_3c, &
         occ, omega_pot, rc_ang, relative_cutoff_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: e_cutoff_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: my_Lrows
      TYPE(cp_eri_mme_param), POINTER                    :: eri_param
      TYPE(cp_fm_type), POINTER                          :: fm_BIb_bse_ab, fm_BIb_bse_ij, fm_BIb_gw, &
                                                            fm_BIb_gw_beta, fm_BIb_jb, &
                                                            fm_BIb_jb_beta, fm_matrix_L
      TYPE(cp_para_env_type), POINTER                    :: para_env_L
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_munu_local_L
      TYPE(dbcsr_t_pgrid_type)                           :: pgrid_t3c_M, pgrid_t3c_overl
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_template, t_3c_tmp
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:, :)   :: t_3c_overl_int
      TYPE(dbcsr_type) :: matrix_bse_ab, matrix_bse_anu, matrix_bse_ij, matrix_bse_inu, &
         matrix_ia_jb, matrix_ia_jb_beta, matrix_ia_jnu, matrix_ia_jnu_beta, matrix_in_jm, &
         matrix_in_jm_beta, matrix_in_jnu, matrix_in_jnu_beta
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_ao, basis_set_ri_aux
      TYPE(neighbor_list_3c_type)                        :: nl_3c
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_p_type)                                    :: pot_g, psi_L, rho_g, rho_r
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(task_list_type), POINTER                      :: task_list_sub

      CALL timeset(routineN, handle)

      CALL cite_reference(DelBen2013)

      do_alpha_beta = .FALSE.
      IF (PRESENT(BIb_C_beta) .AND. &
          PRESENT(gd_B_virtual_beta) .AND. &
          PRESENT(homo_beta) .AND. &
          PRESENT(mo_coeff_o_beta) .AND. &
          PRESENT(mo_coeff_v_beta)) do_alpha_beta = .TRUE.

      IF ((ri_metric%potential_type == do_potential_short .OR. ri_metric%potential_type == do_potential_truncated) &
          .AND. .NOT. do_im_time) THEN
         CPABORT("TRUNCATED and SHORTRANGE RI metrics not yet implemented")
      ENDIF

      virtual = nmo - homo
      gw_corr_lev_total = gw_corr_lev_virt + gw_corr_lev_occ

      eri_method = qs_env%mp2_env%eri_method
      eri_param => qs_env%mp2_env%eri_mme_param
      do_svd = qs_env%mp2_env%do_svd
      eps_svd = qs_env%mp2_env%eps_svd
      potential_type = qs_env%mp2_env%potential_parameter%potential_type
      ri_metric_type = ri_metric%potential_type
      omega_pot = qs_env%mp2_env%potential_parameter%omega

      ! whether we need gpw integrals (plus pw stuff)
      do_gpw = (eri_method == do_eri_gpw) .OR. &
               ((potential_type == do_potential_long .OR. ri_metric_type == do_potential_long) &
                .AND. qs_env%mp2_env%eri_method == do_eri_os) &
               .OR. (ri_metric_type == do_potential_id .AND. qs_env%mp2_env%eri_method == do_eri_mme)

      IF (do_svd .AND. calc_forces) THEN
         CPABORT("SVD not implemented for forces.!")
      END IF

      do_kpoints_from_Gamma = SUM(qs_env%mp2_env%ri_rpa_im_time%kp_grid) > 0
      IF (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         kpoints=kpoints)
      END IF
      IF (do_kpoints_from_Gamma) THEN
         CALL compute_kpoints(qs_env, kpoints, unit_nr)
      END IF

      ! in case of imaginary time, we do not need the intermediate matrices
      IF (.NOT. do_im_time) THEN

         CALL create_intermediate_matrices(matrix_ia_jnu, matrix_ia_jb, mo_coeff_o, virtual, homo, &
                                           fm_BIb_jb, "fm_BIb_jb", max_row_col_local, &
                                           blacs_env_sub, para_env_sub, local_col_row_info)

         CALL create_group_dist(gd_B_virtual, para_env_sub%num_pe, virtual)
         CALL get_group_dist(gd_B_virtual, para_env_sub%mepos, my_B_virtual_start, my_B_virtual_end, my_B_size)

         IF (do_alpha_beta) THEN

            virtual_beta = nmo - homo_beta

            CALL create_intermediate_matrices(matrix_ia_jnu_beta, matrix_ia_jb_beta, mo_coeff_o_beta, &
                                              virtual_beta, homo_beta, &
                                              fm_BIb_jb_beta, "fm_BIb_jb_beta", &
                                              max_row_col_local_beta, &
                                              blacs_env_sub, para_env_sub, local_col_row_info_beta)

            CALL create_group_dist(gd_B_virtual_beta, para_env_sub%num_pe, virtual_beta)
            CALL get_group_dist(gd_B_virtual_beta, para_env_sub%mepos, my_B_virtual_start_beta, my_B_virtual_end_beta, &
                                my_B_size_beta)

         END IF

         ! in the case of G0W0, we need (K|nm), n,m may be occ or virt (m restricted to corrected levels)
         IF (my_do_gw) THEN

            CALL create_intermediate_matrices(matrix_in_jnu, matrix_in_jm, mo_coeff_gw, &
                                              nmo, gw_corr_lev_total, &
                                              fm_BIb_gw, "fm_BIb_gw", &
                                              max_row_col_local_gw, &
                                              blacs_env_sub, para_env_sub, local_col_row_info_gw)

            CALL create_group_dist(gd_B_all, para_env_sub%num_pe, nmo)
            CALL get_group_dist(gd_B_all, para_env_sub%mepos, my_B_all_start, my_B_all_end, my_B_all_size)

            IF (do_alpha_beta) THEN
               ! deallocate local_col_row_info_gw, otherwise it gets twice allocated in create_intermediate_m
               DEALLOCATE (local_col_row_info_gw)
               CALL create_intermediate_matrices(matrix_in_jnu_beta, matrix_in_jm_beta, mo_coeff_gw_beta, &
                                                 nmo, gw_corr_lev_total, &
                                                 fm_BIb_gw_beta, "fm_BIb_gw_beta", &
                                                 max_row_col_local_gw, &
                                                 blacs_env_sub, para_env_sub, local_col_row_info_gw)

               ! we don"t need parallelization arrays for beta since the matrix sizes of B_nm^P is the same
               ! for the beta case and therefore the parallelization of beta is the same than for alpha

            END IF
         END IF
      END IF ! not for imaginary time

      IF (do_bse) THEN

         CPASSERT(my_do_gw)
         CPASSERT(.NOT. do_im_time)
         ! GPW integrals have to be implemented later
         CPASSERT(.NOT. (eri_method == do_eri_gpw))

         ! virt x virt matrices
         CALL create_intermediate_matrices(matrix_bse_anu, matrix_bse_ab, mo_coeff_v, virtual, virtual, &
                                           fm_BIb_bse_ab, "fm_BIb_bse_ab", max_row_col_local_virt_bse, &
                                           blacs_env_sub, para_env_sub, local_col_row_info_virt_bse)

         CALL create_group_dist(gd_B_virt_bse, para_env_sub%num_pe, virtual)
         CALL get_group_dist(gd_B_virt_bse, para_env_sub%mepos, my_B_virt_bse_start, my_B_virt_bse_end, my_B_virt_bse_size)

         ! occ x occ matrices
         CALL create_intermediate_matrices(matrix_bse_inu, matrix_bse_ij, mo_coeff_o, homo, homo, &
                                           fm_BIb_bse_ij, "fm_BIb_bse_ij", max_row_col_local_occ_bse, &
                                           blacs_env_sub, para_env_sub, local_col_row_info_occ_bse)

         CALL create_group_dist(gd_B_occ_bse, para_env_sub%num_pe, homo)
         CALL get_group_dist(gd_B_occ_bse, para_env_sub%mepos, my_B_occ_bse_start, my_B_occ_bse_end, my_B_occ_bse_size)

      END IF

      ngroup = para_env%num_pe/para_env_sub%num_pe

      ! Preparations for MME method to compute ERIs
      IF (qs_env%mp2_env%eri_method .EQ. do_eri_mme) THEN
         ! cell might have changed, so we need to reset parameters
         CALL cp_eri_mme_set_params(eri_param, cell, qs_kind_set, basis_type_1="ORB", basis_type_2="RI_AUX", para_env=para_env)
      ENDIF

      CALL get_cell(cell=cell, periodic=periodic)
      ! for minimax Ewald summation, full periodicity is required
      IF (eri_method == do_eri_mme) THEN
         CPASSERT(periodic(1) == 1 .AND. periodic(2) == 1 .AND. periodic(3) == 1)
      END IF

      IF (do_svd .AND. (do_kpoints_from_Gamma .OR. do_kpoints_cubic_RPA)) THEN
         CPABORT("SVD with kpoints not implemented yet!")
      END IF

      CALL get_2c_integrals(qs_env, eri_method, eri_param, para_env, para_env_sub, para_env_L, mp2_memory, &
                            fm_matrix_L, ngroup, color_sub, dimen_RI, dimen_RI_red, kpoints, mo_coeff, &
                            my_group_L_size, my_group_L_start, my_group_L_end, &
                            gd_array, calc_PQ_cond_num .AND. .NOT. do_svd, cond_num, do_svd, eps_svd, &
                            num_small_eigen, qs_env%mp2_env%potential_parameter, ri_metric, &
                            fm_matrix_L_RI_metric, &
                            do_im_time, do_kpoints_from_Gamma .OR. do_kpoints_cubic_RPA, qs_env%mp2_env%mp2_gpw%eps_pgf_orb_S, &
                            qs_kind_set, sab_orb_sub)

      IF (unit_nr > 0) THEN
         ASSOCIATE (ri_metric=>qs_env%mp2_env%ri_metric)
            SELECT CASE (ri_metric%potential_type)
            CASE (do_potential_coulomb)
               WRITE (unit_nr, FMT="(/T3,A,T74,A)") &
                  "RI_INFO| RI metric: ", "COULOMB"
            CASE (do_potential_short)
               WRITE (unit_nr, FMT="(T3,A,T71,A)") &
                  "RI_INFO| RI metric: ", "SHORTRANGE"
               WRITE (unit_nr, '(T3,A,T61,F20.10)') &
                  "RI_INFO| Omega:     ", ri_metric%omega
               rc_ang = cp_unit_from_cp2k(ri_metric%cutoff_radius, "angstrom")
               WRITE (unit_nr, '(T3,A,T61,F20.10)') &
                  "RI_INFO| Cutoff Radius [angstrom]:     ", rc_ang
            CASE (do_potential_long)
               WRITE (unit_nr, FMT="(T3,A,T72,A)") &
                  "RI_INFO| RI metric: ", "LONGRANGE"
               WRITE (unit_nr, '(T3,A,T61,F20.10)') &
                  "RI_INFO| Omega:     ", ri_metric%omega
            CASE (do_potential_id)
               WRITE (unit_nr, FMT="(T3,A,T73,A)") &
                  "RI_INFO| RI metric: ", "IDENTITY"
            CASE (do_potential_truncated)
               WRITE (unit_nr, FMT="(T3,A,T72,A)") &
                  "RI_INFO| RI metric: ", "TRUNCATED"
               rc_ang = cp_unit_from_cp2k(ri_metric%cutoff_radius, "angstrom")
               WRITE (unit_nr, '(T3,A,T61,F20.2)') &
                  "RI_INFO| Cutoff Radius [angstrom]:     ", rc_ang
            END SELECT
         END ASSOCIATE
      ENDIF

      IF (unit_nr > 0 .AND. .NOT. do_im_time) THEN
         WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") "RI_INFO| Cholesky decomposition group size:", para_env_L%num_pe
         WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") "RI_INFO| Number of groups for auxiliary basis functions", ngroup
      END IF
      IF (calc_PQ_cond_num .OR. do_svd) THEN
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T67,ES14.5)") &
            "RI_INFO| Condition number of the (P|Q):", cond_num
         IF (do_svd) THEN
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,ES8.1,A,T75,i6)") &
               "RI_INFO| Number of neglected Eigenvalues of (P|Q) smaller than ", eps_svd, ":", num_small_eigen
         ELSE
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,ES14.1,A,T75,i6)") &
               "RI_INFO| Number of Eigenvalue of (P|Q) smaller ", eps_svd, ":", num_small_eigen
         END IF
      END IF

      IF (.NOT. do_im_time) THEN
         ! replicate the necessary row of the L^{-1} matrix on each proc
         CALL grep_Lcols(para_env_L, dimen_RI, fm_matrix_L, &
                         my_group_L_start, my_group_L_end, my_group_L_size, my_Lrows)
      END IF
      ! clean the L^{-1} matrix
      CALL cp_fm_release(fm_matrix_L)

      ! in case of imag. time we need the para_env_L later
      IF (.NOT. do_im_time) THEN
         CALL cp_para_env_release(para_env_L)
      END IF

      IF (calc_forces) THEN
         ! we need (P|Q)^(-1/2) for future use, just save it
         ! in a fully (home made) distributed way
         itmp = get_limit(dimen_RI, para_env_sub%num_pe, para_env_sub%mepos)
         lll = itmp(2) - itmp(1) + 1
         ALLOCATE (qs_env%mp2_env%ri_grad%PQ_half(lll, my_group_L_size))
         qs_env%mp2_env%ri_grad%PQ_half(:, :) = my_Lrows(itmp(1):itmp(2), 1:my_group_L_size)
      END IF

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
            "RI_INFO| Number of auxiliary basis functions:", dimen_RI, &
            "GENERAL_INFO| Number of basis functions:", homo + virtual, &
            "GENERAL_INFO| Number of occupied orbitals:", homo, &
            "GENERAL_INFO| Number of virtual orbitals:", virtual
         IF (do_svd) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "RI_INFO| Reduced auxiliary basis set size:", dimen_RI_red
         END IF

         mem_for_iaK = dimen_RI*REAL(homo, KIND=dp)*virtual*8.0_dp/(1024_dp**2)
         IF (do_alpha_beta) mem_for_iaK = mem_for_iaK + &
                                          dimen_RI*REAL(homo_beta, KIND=dp)*(nmo - homo_beta)*8.0_dp/(1024_dp**2)

         IF (.NOT. do_im_time) THEN
            WRITE (unit_nr, '(T3,A,T66,F11.2,A4)') 'RI_INFO| Total memory for (ia|K) integrals:', &
               mem_for_iaK, ' MiB'
            IF (my_do_gw .AND. .NOT. do_im_time) THEN
               mem_for_iaK = dimen_RI*REAL(nmo, KIND=dp)*gw_corr_lev_total*8.0_dp/(1024_dp**2)

               WRITE (unit_nr, '(T3,A,T66,F11.2,A4)') 'RI_INFO| Total memory for G0W0-(nm|K) integrals:', &
                  mem_for_iaK, ' MiB'
            END IF
         END IF
         CALL m_flush(unit_nr)
      ENDIF

      CALL mp_sync(para_env%group) ! sync to see memory output

      ! in case we do imaginary time, we need the overlap tensor (alpha beta P) or trunc. Coulomb tensor
      IF (.NOT. do_im_time) THEN

         ALLOCATE (sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe - 1))
         DO i = 0, para_env_sub%num_pe - 1
            sub_proc_map(i) = i
            sub_proc_map(-i - 1) = para_env_sub%num_pe - i - 1
            sub_proc_map(para_env_sub%num_pe + i) = i
         END DO

         ! array that will store the (ia|K) integrals
         ALLOCATE (BIb_C(my_group_L_size, my_B_size, homo))
         BIb_C = 0.0_dp

         IF (do_alpha_beta) THEN
            ALLOCATE (BIb_C_beta(my_group_L_size, my_B_size_beta, homo_beta))
            BIb_C_beta = 0.0_dp
         END IF

         ! in the case of GW, we also need (nm|K)
         IF (my_do_gw) THEN

            ALLOCATE (BIb_C_gw(my_group_L_size, my_B_all_size, gw_corr_lev_total))
            BIb_C_gw = 0.0_dp
            IF (do_alpha_beta) THEN
               ALLOCATE (BIb_C_gw_beta(my_group_L_size, my_B_all_size, gw_corr_lev_total))
               BIb_C_gw_beta = 0.0_dp
            END IF

         END IF

         IF (do_bse) THEN

            ALLOCATE (BIb_C_bse_ij(my_group_L_size, my_B_occ_bse_size, homo))
            BIb_C_bse_ij = 0.0_dp

            ALLOCATE (BIb_C_bse_ab(my_group_L_size, my_B_virt_bse_size, virtual))
            BIb_C_bse_ab = 0.0_dp

         END IF

         CALL timeset(routineN//"_loop", handle2)

         IF (eri_method == do_eri_mme .AND. &
             (ri_metric%potential_type == do_potential_coulomb .OR. ri_metric%potential_type == do_potential_long) .OR. &
             eri_method == do_eri_os .AND. ri_metric%potential_type == do_potential_coulomb) THEN

            NULLIFY (mat_munu_local_L)
            ALLOCATE (mat_munu_local_L(my_group_L_size))
            DO LLL = 1, my_group_L_size
               NULLIFY (mat_munu_local_L(LLL)%matrix)
               ALLOCATE (mat_munu_local_L(LLL)%matrix)
               CALL dbcsr_copy(mat_munu_local_L(LLL)%matrix, mat_munu%matrix)
               CALL dbcsr_set(mat_munu_local_L(LLL)%matrix, 0.0_dp)
            ENDDO
            CALL mp2_eri_3c_integrate(eri_param, ri_metric%potential_type, ri_metric%omega, para_env_sub, qs_env, &
                                      first_c=my_group_L_start, last_c=my_group_L_end, &
                                      mat_ab=mat_munu_local_L, &
                                      basis_type_a="ORB", basis_type_b="ORB", &
                                      basis_type_c="RI_AUX", &
                                      sab_nl=sab_orb_sub, eri_method=eri_method)

            DO LLL = 1, my_group_L_size
               CALL ao_to_mo_and_store_B(para_env_sub, mat_munu_local_L(LLL), matrix_ia_jnu, matrix_ia_jb, &
                                         fm_BIb_jb, BIb_C(LLL, 1:my_B_size, 1:homo), &
                                         mo_coeff_o, mo_coeff_v, eps_filter, max_row_col_local, sub_proc_map, &
                                         local_col_row_info, my_B_virtual_end, my_B_virtual_start, "alpha")
            ENDDO
            CALL contract_B_L(BIb_C, my_Lrows, gd_B_virtual%sizes, gd_array%sizes, qs_env%mp2_env%eri_blksize, &
                              ngroup, color_sub, para_env%group, para_env_sub)

            IF (do_alpha_beta) THEN

               DO LLL = 1, my_group_L_size
                  CALL ao_to_mo_and_store_B(para_env_sub, mat_munu_local_L(LLL), matrix_ia_jnu_beta, matrix_ia_jb_beta, &
                                            fm_BIb_jb_beta, BIb_C_beta(LLL, 1:my_B_size_beta, 1:homo_beta), &
                                            mo_coeff_o_beta, mo_coeff_v_beta, eps_filter, max_row_col_local_beta, sub_proc_map, &
                                            local_col_row_info_beta, my_B_virtual_end_beta, my_B_virtual_start_beta, "beta")
               ENDDO
               CALL contract_B_L(BIb_C_beta, my_Lrows, gd_B_virtual_beta%sizes, gd_array%sizes, qs_env%mp2_env%eri_blksize, &
                                 ngroup, color_sub, para_env%group, para_env_sub)

            ENDIF

            IF (my_do_gw) THEN

               DO LLL = 1, my_group_L_size
                  CALL ao_to_mo_and_store_B(para_env_sub, mat_munu_local_L(LLL), matrix_in_jnu, matrix_in_jm, &
                                            fm_BIb_gw, BIb_C_gw(LLL, 1:my_B_all_size, 1:gw_corr_lev_total), &
                                            mo_coeff_gw, mo_coeff_all, eps_filter, max_row_col_local_gw, sub_proc_map, &
                                            local_col_row_info_gw, my_B_all_end, my_B_all_start, "gw_alpha")
               ENDDO
               CALL contract_B_L(BIb_C_gw, my_Lrows, gd_B_all%sizes, gd_array%sizes, qs_env%mp2_env%eri_blksize, &
                                 ngroup, color_sub, para_env%group, para_env_sub)

               IF (do_alpha_beta) THEN

                  DO LLL = 1, my_group_L_size
                     CALL ao_to_mo_and_store_B(para_env_sub, mat_munu_local_L(LLL), matrix_in_jnu_beta, matrix_in_jm_beta, &
                                               fm_BIb_gw_beta, BIb_C_gw_beta(LLL, 1:my_B_all_size, 1:gw_corr_lev_total), &
                                               mo_coeff_gw_beta, mo_coeff_all_beta, eps_filter, max_row_col_local_gw, &
                                               sub_proc_map, local_col_row_info_gw, my_B_all_end, my_B_all_start, "gw_beta")
                  ENDDO
                  CALL contract_B_L(BIb_C_gw_beta, my_Lrows, gd_B_all%sizes, gd_array%sizes, qs_env%mp2_env%eri_blksize, &
                                    ngroup, color_sub, para_env%group, para_env_sub)

               ENDIF
            ENDIF

            IF (do_bse) THEN

               ! B^ab_P matrix elements for BSE
               DO LLL = 1, my_group_L_size
                  CALL ao_to_mo_and_store_B(para_env_sub, mat_munu_local_L(LLL), matrix_bse_anu, matrix_bse_ab, &
                                            fm_BIb_bse_ab, BIb_C_bse_ab(LLL, 1:my_B_virt_bse_size, 1:virtual), &
                                            mo_coeff_v, mo_coeff_v, eps_filter, max_row_col_local_virt_bse, &
                                            sub_proc_map, local_col_row_info_virt_bse, my_B_all_end, my_B_all_start, "bse_ab")
               ENDDO
               CALL contract_B_L(BIb_C_bse_ab, my_Lrows, gd_B_virt_bse%sizes, gd_array%sizes, qs_env%mp2_env%eri_blksize, &
                                 ngroup, color_sub, para_env%group, para_env_sub)

               ! B^ij_P matrix elements for BSE
               DO LLL = 1, my_group_L_size
                  CALL ao_to_mo_and_store_B(para_env_sub, mat_munu_local_L(LLL), matrix_bse_inu, matrix_bse_ij, &
                                            fm_BIb_bse_ij, BIb_C_bse_ij(LLL, 1:my_B_occ_bse_size, 1:homo), &
                                            mo_coeff_o, mo_coeff_o, eps_filter, max_row_col_local_occ_bse, &
                                            sub_proc_map, local_col_row_info_occ_bse, &
                                            my_B_occ_bse_end, my_B_occ_bse_start, "bse_ij")
               ENDDO
               CALL contract_B_L(BIb_C_bse_ij, my_Lrows, gd_B_occ_bse%sizes, gd_array%sizes, qs_env%mp2_env%eri_blksize, &
                                 ngroup, color_sub, para_env%group, para_env_sub)

            END IF

            DO LLL = 1, my_group_L_size
               CALL dbcsr_release_p(mat_munu_local_L(LLL)%matrix)
            ENDDO
            DEALLOCATE (mat_munu_local_L)

         ELSE IF (do_gpw) THEN

            CALL prepare_gpw(qs_env, dft_control, e_cutoff_old, cutoff_old, relative_cutoff_old, para_env_sub, pw_env_sub, &
                             auxbas_pw_pool, poisson_env, task_list_sub, rho_r, rho_g, pot_g, psi_L, sab_orb_sub)

            DO i_counter = 1, my_group_L_size

               CALL mp2_eri_3c_integrate_gpw(mo_coeff, psi_L, rho_g, atomic_kind_set, qs_kind_set, cell, dft_control, &
                                             particle_set, pw_env_sub, my_Lrows(:, i_counter), poisson_env, rho_r, pot_g, &
                                             ri_metric%potential_type, ri_metric%omega, mat_munu, qs_env, task_list_sub)

               CALL ao_to_mo_and_store_B(para_env_sub, mat_munu, matrix_ia_jnu, matrix_ia_jb, &
                                         fm_BIb_jb, BIb_C(i_counter, 1:my_B_size, 1:homo), &
                                         mo_coeff_o, mo_coeff_v, eps_filter, max_row_col_local, sub_proc_map, &
                                         local_col_row_info, my_B_virtual_end, my_B_virtual_start, "alpha")

               IF (do_alpha_beta) THEN
                  CALL ao_to_mo_and_store_B(para_env_sub, mat_munu, matrix_ia_jnu_beta, matrix_ia_jb_beta, &
                                            fm_BIb_jb_beta, BIb_C_beta(i_counter, 1:my_B_size_beta, 1:homo_beta), &
                                            mo_coeff_o_beta, mo_coeff_v_beta, eps_filter, max_row_col_local_beta, sub_proc_map, &
                                            local_col_row_info_beta, my_B_virtual_end_beta, my_B_virtual_start_beta, "beta")

               ENDIF

               IF (my_do_gw) THEN
                  ! transform (K|mu nu) to (K|nm), n corresponds to corrected GW levels, m is in nmo
                  CALL ao_to_mo_and_store_B(para_env_sub, mat_munu, matrix_in_jnu, matrix_in_jm, &
                                            fm_BIb_gw, BIb_C_gw(i_counter, 1:my_B_all_size, 1:gw_corr_lev_total), &
                                            mo_coeff_gw, mo_coeff_all, eps_filter, max_row_col_local_gw, sub_proc_map, &
                                            local_col_row_info_gw, my_B_all_end, my_B_all_start, "gw_alpha")

                  ! the same for beta
                  IF (do_alpha_beta) THEN
                     CALL ao_to_mo_and_store_B(para_env_sub, mat_munu, matrix_in_jnu_beta, matrix_in_jm_beta, &
                                               fm_BIb_gw_beta, BIb_C_gw_beta(i_counter, 1:my_B_all_size, 1:gw_corr_lev_total), &
                                               mo_coeff_gw_beta, mo_coeff_all_beta, eps_filter, max_row_col_local_gw, &
                                               sub_proc_map, local_col_row_info_gw, my_B_all_end, my_B_all_start, "gw_beta")

                  END IF
               END IF

            END DO

            CALL cleanup_gpw(qs_env, e_cutoff_old, cutoff_old, relative_cutoff_old, pw_env_sub, &
                             task_list_sub, auxbas_pw_pool, rho_r, rho_g, pot_g, psi_L)
         ELSE
            CPABORT("Integration method not implemented!")
         END IF

         CALL timestop(handle2)

         DEALLOCATE (my_Lrows)

         CALL cp_fm_release(fm_BIb_jb)
         DEALLOCATE (local_col_row_info)

         CALL dbcsr_release(matrix_ia_jnu)
         CALL dbcsr_release(matrix_ia_jb)
         IF (do_alpha_beta) THEN
            CALL dbcsr_release(matrix_ia_jnu_beta)
            CALL dbcsr_release(matrix_ia_jb_beta)
            CALL cp_fm_release(fm_BIb_jb_beta)
            DEALLOCATE (local_col_row_info_beta)
         END IF

         IF (my_do_gw) THEN
            CALL dbcsr_release(matrix_in_jnu)
            CALL dbcsr_release(matrix_in_jm)
            CALL cp_fm_release(fm_BIb_gw)
            DEALLOCATE (local_col_row_info_gw)
            IF (do_alpha_beta) THEN
               CALL dbcsr_release(matrix_in_jnu_beta)
               CALL dbcsr_release(matrix_in_jm_beta)
               CALL cp_fm_release(fm_BIb_gw_beta)
            END IF
         END IF

         IF (do_bse) THEN
            CALL dbcsr_release(matrix_bse_anu)
            CALL dbcsr_release(matrix_bse_ab)
            CALL cp_fm_release(fm_BIb_bse_ab)
            CALL dbcsr_release(matrix_bse_inu)
            CALL dbcsr_release(matrix_bse_ij)
            CALL cp_fm_release(fm_BIb_bse_ij)
         END IF

         DEALLOCATE (sub_proc_map)

         ! imag. time = low-scaling SOS-MP2, RPA, GW
      ELSE

         memory_info = qs_env%mp2_env%ri_rpa_im_time%memory_info

         ! we need 3 tensors:
         ! 1) t_3c_overl_int: 3c overlap integrals, optimized for easy access to integral blocks
         !                   (atomic blocks)
         ! 2) t_3c_O: 3c overlap integrals, optimized for contraction (split blocks)
         ! 3) t_3c_M: tensor M, optimized for contraction

         CALL get_qs_env(qs_env, natom=natom, nkind=nkind, dft_control=dft_control)

         pdims_t3c = 0
         CALL dbcsr_t_pgrid_create(para_env%group, pdims_t3c, pgrid_t3c_overl)

         ! set up basis
         ALLOCATE (sizes_RI(natom), sizes_AO(natom))
         ALLOCATE (basis_set_ri_aux(nkind), basis_set_ao(nkind))
         CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
         CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=basis_set_ri_aux)
         CALL basis_set_list_setup(basis_set_ao, "ORB", qs_kind_set)
         CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_AO, basis=basis_set_ao)

         cut_memory_int = qs_env%mp2_env%ri_rpa_im_time%cut_memory
         CALL create_tensor_batches(sizes_AO, cut_memory_int, starts_array_mc_int, ends_array_mc_int, &
                                    starts_array_mc_block_int, ends_array_mc_block_int)

         DEALLOCATE (starts_array_mc_int, ends_array_mc_int)

         CALL create_3c_tensor(t_3c_overl_int_template, dist_RI, dist_AO_1, dist_AO_2, pgrid_t3c_overl, &
                               sizes_RI, sizes_AO, sizes_AO, map1=[1, 2], map2=[3], &
                               name="O (RI AO | AO)")

         CALL get_qs_env(qs_env, nkind=nkind, particle_set=particle_set)
         CALL dbcsr_t_mp_environ_pgrid(pgrid_t3c_overl, pdims, pcoord)
         CALL mp_cart_create(pgrid_t3c_overl%mp_comm_2d, 3, pdims, pcoord, mp_comm_t3c_2)
         CALL distribution_3d_create(dist_3d, dist_RI, dist_AO_1, dist_AO_2, &
                                     nkind, particle_set, mp_comm_t3c_2, own_comm=.TRUE.)
         DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)

         CALL build_3c_neighbor_lists(nl_3c, basis_set_ri_aux, basis_set_ao, basis_set_ao, &
                                      dist_3d, ri_metric, "RPA_3c_nl", qs_env, &
                                      sym_jk=.NOT. do_kpoints_cubic_RPA, own_dist=.TRUE.)

         ! init k points
         IF (do_kpoints_cubic_RPA) THEN
            ! set up new kpoint type with periodic images according to eps_grid from MP2 section
            ! instead of eps_pgf_orb from QS section
            CALL kpoint_init_cell_index(kpoints, nl_3c%jk_list, para_env, dft_control)
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "3C_OVERLAP_INTEGRALS_INFO| Number of periodic images considered:", dft_control%nimages

            nimg = dft_control%nimages
         ELSE
            nimg = 1
         ENDIF

         ALLOCATE (t_3c_overl_int(nimg, nimg))

         DO i = 1, SIZE(t_3c_overl_int, 1)
            DO j = 1, SIZE(t_3c_overl_int, 2)
               CALL dbcsr_t_create(t_3c_overl_int_template, t_3c_overl_int(i, j))
            ENDDO
         ENDDO

         CALL dbcsr_t_destroy(t_3c_overl_int_template)

         ! split blocks to improve load balancing for tensor contraction
         min_bsize = qs_env%mp2_env%ri_rpa_im_time%min_bsize

         CALL pgf_block_sizes(atomic_kind_set, basis_set_ao, min_bsize, sizes_AO_split)
         CALL pgf_block_sizes(atomic_kind_set, basis_set_ri_aux, min_bsize, sizes_RI_split)

         pdims_t3c = 0
         CALL dbcsr_t_pgrid_create(para_env%group, pdims_t3c, pgrid_t3c_M)

         ASSOCIATE (cut_memory=>qs_env%mp2_env%ri_rpa_im_time%cut_memory)
            CALL create_tensor_batches(sizes_AO_split, cut_memory, starts_array_mc, ends_array_mc, &
                                       starts_array_mc_block, ends_array_mc_block)
         END ASSOCIATE
         cut_memory = qs_env%mp2_env%ri_rpa_im_time%cut_memory

         CALL create_3c_tensor(t_3c_M, dist_RI, dist_AO_1, dist_AO_2, pgrid_t3c_M, &
                               sizes_RI_split, sizes_AO_split, sizes_AO_split, &
                               map1=[1], map2=[2, 3], &
                               name="M (RI | AO AO)")

         DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)

         CALL dbcsr_t_pgrid_destroy(pgrid_t3c_M)

         ALLOCATE (t_3c_O(SIZE(t_3c_overl_int, 1), SIZE(t_3c_overl_int, 2)))
         ALLOCATE (t_3c_O_compressed(SIZE(t_3c_overl_int, 1), SIZE(t_3c_overl_int, 2), cut_memory))
         ALLOCATE (t_3c_O_ind(SIZE(t_3c_overl_int, 1), SIZE(t_3c_overl_int, 2), cut_memory))
         CALL create_3c_tensor(t_3c_O(1, 1), dist_RI, dist_AO_1, dist_AO_2, pgrid_t3c_overl, &
                               sizes_RI_split, sizes_AO_split, sizes_AO_split, &
                               map1=[1, 2], map2=[3], &
                               name="O (RI AO | AO)")
         DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)

         CALL dbcsr_t_pgrid_destroy(pgrid_t3c_overl)
         DO i = 1, SIZE(t_3c_O, 1)
            DO j = 1, SIZE(t_3c_O, 2)
               IF (i > 1 .OR. j > 1) CALL dbcsr_t_create(t_3c_O(1, 1), t_3c_O(i, j))
            ENDDO
         ENDDO

         ! build integrals in batches and copy to optimized format
         ! note: integrals are stored in terms of atomic blocks. To avoid a memory bottleneck,
         ! integrals are calculated in batches and copied to optimized format with subatomic blocks

         DO cm = 1, cut_memory_int
            CALL build_3c_integrals(t_3c_overl_int, &
                                    qs_env%mp2_env%ri_rpa_im_time%eps_filter/2, &
                                    qs_env, &
                                    nl_3c, &
                                    basis_i=basis_set_ri_aux, &
                                    basis_j=basis_set_ao, basis_k=basis_set_ao, &
                                    potential_parameter=ri_metric, &
                                    do_kpoints=do_kpoints_cubic_RPA, &
                                    bounds_j=[starts_array_mc_block_int(cm), ends_array_mc_block_int(cm)], desymmetrize=.FALSE.)
            CALL timeset(routineN//"_copy_3c", handle4)
            ! copy integral tensor t_3c_overl_int to t_3c_O tensor optimized for contraction
            DO i = 1, SIZE(t_3c_overl_int, 1)
               DO j = 1, SIZE(t_3c_overl_int, 2)

                  CALL dbcsr_t_copy(t_3c_overl_int(i, j), t_3c_O(i, j), order=[1, 3, 2], &
                                    summation=.TRUE., move_data=.TRUE.)
                  CALL dbcsr_t_clear(t_3c_overl_int(i, j))
                  CALL dbcsr_t_filter(t_3c_O(i, j), qs_env%mp2_env%ri_rpa_im_time%eps_filter/2)
               ENDDO
            ENDDO

            CALL timestop(handle4)
         ENDDO

         DO i = 1, SIZE(t_3c_overl_int, 1)
            DO j = 1, SIZE(t_3c_overl_int, 2)
               CALL dbcsr_t_destroy(t_3c_overl_int(i, j))
            ENDDO
         ENDDO
         DEALLOCATE (t_3c_overl_int)

         CALL timeset(routineN//"_copy_3c", handle4)
         ! desymmetrize
         CALL dbcsr_t_create(t_3c_O(1, 1), t_3c_tmp)
         DO jcell = 1, nimg
            DO kcell = 1, jcell
               CALL dbcsr_t_copy(t_3c_O(jcell, kcell), t_3c_tmp)
               CALL dbcsr_t_copy(t_3c_tmp, t_3c_O(kcell, jcell), order=[1, 3, 2], summation=.TRUE., move_data=.TRUE.)
               CALL dbcsr_t_filter(t_3c_O(kcell, jcell), qs_env%mp2_env%ri_rpa_im_time%eps_filter)
            ENDDO
         ENDDO
         DO jcell = 1, nimg
            DO kcell = jcell + 1, nimg
               CALL dbcsr_t_copy(t_3c_O(jcell, kcell), t_3c_tmp)
               CALL dbcsr_t_copy(t_3c_tmp, t_3c_O(kcell, jcell), order=[1, 3, 2], summation=.FALSE., move_data=.TRUE.)
               CALL dbcsr_t_filter(t_3c_O(kcell, jcell), qs_env%mp2_env%ri_rpa_im_time%eps_filter)
            ENDDO
         ENDDO

         CALL dbcsr_t_get_info(t_3c_O(1, 1), nfull_total=bounds_3c)
         CALL get_tensor_occupancy(t_3c_O(1, 1), nze, occ)
         memory_3c = 0.0_dp

         bounds(:, 1) = [1, bounds_3c(1)]
         bounds(:, 3) = [1, bounds_3c(3)]
         DO i = 1, SIZE(t_3c_O, 1)
            DO j = 1, SIZE(t_3c_O, 2)
               DO i_mem = 1, cut_memory
                  bounds(:, 2) = [starts_array_mc(i_mem), ends_array_mc(i_mem)]
                  CALL dbcsr_t_copy(t_3c_O(i, j), t_3c_tmp, bounds=bounds)

                  ALLOCATE (t_3c_O_ind(i, j, i_mem)%ind(dbcsr_t_get_num_blocks(t_3c_tmp), 3))
                  CALL dbcsr_t_reserved_block_indices(t_3c_tmp, t_3c_O_ind(i, j, i_mem)%ind)
                  CALL alloc_containers(t_3c_O_compressed(i, j, i_mem), 1)
                  CALL compress_tensor(t_3c_tmp, t_3c_O_compressed(i, j, i_mem), &
                                       qs_env%mp2_env%ri_rpa_im_time%eps_compress, memory_3c)
               ENDDO
               CALL dbcsr_t_clear(t_3c_O(i, j))
            ENDDO
         ENDDO

         CALL mp_sum(memory_3c, para_env%group)

         compression_factor = REAL(nze, dp)*1.0E-06*8.0_dp/memory_3c

         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="((T3,A,T66,F11.2,A4))") &
               "MEMORY_INFO| Memory for 3-center integrals (compressed):", memory_3c, ' MiB'

            WRITE (UNIT=unit_nr, FMT="((T3,A,T60,F21.2))") &
               "MEMORY_INFO| Compression factor:                  ", compression_factor
         END IF

         CALL dbcsr_t_destroy(t_3c_tmp)

         CALL timestop(handle4)

         DEALLOCATE (basis_set_ri_aux, basis_set_ao)

         CALL neighbor_list_3c_destroy(nl_3c)

         CALL cp_para_env_release(para_env_L)

      END IF

      CALL timestop(handle)

   END SUBROUTINE mp2_ri_gpw_compute_in

! **************************************************************************************************
!> \brief Contract (P|ai) = (R|P) x (R|ai)
!> \param BIb_C (R|ai)
!> \param my_Lrows (R|P)
!> \param sizes_B number of a (virtual) indices per subgroup process
!> \param sizes_L number of P / R (auxiliary) indices per subgroup
!> \param blk_size ...
!> \param ngroup how many subgroups (NG)
!> \param igroup subgroup color
!> \param mp_comm communicator
!> \param para_env_sub ...
! **************************************************************************************************
   SUBROUTINE contract_B_L(BIb_C, my_Lrows, sizes_B, sizes_L, blk_size, ngroup, igroup, mp_comm, para_env_sub)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: BIb_C
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: my_Lrows
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_B, sizes_L
      INTEGER, DIMENSION(2), INTENT(IN)                  :: blk_size
      INTEGER, INTENT(IN)                                :: ngroup, igroup, mp_comm
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_B_L'
      LOGICAL, PARAMETER                                 :: debug = .FALSE.

      INTEGER                                            :: check_proc, handle, i, iend, ii, ioff, &
                                                            iproc, iproc_glob, istart, loc_a, &
                                                            loc_P, nproc, nproc_glob
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: block_ind_L_P, block_ind_L_R
      INTEGER, DIMENSION(1)                              :: dist_B_i, map_B_1, map_L_1, map_L_2, &
                                                            sizes_i
      INTEGER, DIMENSION(2)                              :: map_B_2, pdims_L
      INTEGER, DIMENSION(3)                              :: pdims_B
      LOGICAL                                            :: found
      INTEGER, DIMENSION(ngroup)                         :: dist_L_P, dist_L_R
      INTEGER, DIMENSION(para_env_sub%num_pe)            :: dist_B_a
      TYPE(dbcsr_t_distribution_type)                    :: dist_B, dist_L
      TYPE(dbcsr_t_pgrid_type)                           :: mp_comm_B, mp_comm_L
      TYPE(dbcsr_t_type)                                 :: tB_in, tB_in_split, tB_out, &
                                                            tB_out_split, tL, tL_split

      CALL timeset(routineN, handle)

      sizes_i(1) = SIZE(BIb_C, 3)

      nproc = para_env_sub%num_pe ! number of processes per subgroup (Nw)
      iproc = para_env_sub%mepos ! subgroup-local process ID

      ! Total number of processes and global process ID
      CALL mp_environ(nproc_glob, iproc_glob, mp_comm)

      ! local block index for R/P and a
      loc_P = igroup + 1; loc_a = iproc + 1

      CPASSERT(SIZE(sizes_L) .EQ. ngroup)
      CPASSERT(SIZE(sizes_B) .EQ. nproc)
      CPASSERT(sizes_L(loc_P) .EQ. SIZE(BIb_C, 1))
      CPASSERT(sizes_L(loc_P) .EQ. SIZE(my_Lrows, 2))
      CPASSERT(sizes_B(loc_a) .EQ. SIZE(BIb_C, 2))

      ! Tensor distributions as follows:
      ! Process grid NG x Nw
      ! Each process has coordinates (np, nw)
      ! tB_in: (R|ai): R distributed (np), a distributed (nw)
      ! tB_out: (P|ai): P distributed (np), a distributed (nw)
      ! tL: (R|P): R distributed (nw), P distributed (np)

      ! define mappings between tensor index and matrix index:
      ! (R|ai) and (P|ai):
      map_B_1 = [1] ! index 1 (R or P) maps to 1st matrix index (np distributed)
      map_B_2 = [2, 3] ! indices 2, 3 (a, i) map to 2nd matrix index (nw distributed)
      ! (R|P):
      map_L_1 = [2] ! index 2 (P) maps to 1st matrix index (np distributed)
      map_L_2 = [1] ! index 1 (R) maps to 2nd matrix index (nw distributed)

      ! derive nd process grid that is compatible with distributions and 2d process grid
      ! (R|ai) / (P|ai) on process grid NG x Nw x 1
      ! (R|P) on process grid NG x Nw
      pdims_B = [ngroup, nproc, 1]
      pdims_L = [nproc, ngroup]

      CALL dbcsr_t_pgrid_create(mp_comm, pdims_B, mp_comm_B)
      CALL dbcsr_t_pgrid_create(mp_comm, pdims_L, mp_comm_L)

      ! setup distribution vectors such that distribution matches parallel data layout of BIb_C and my_Lrows
      dist_B_i = [0]
      dist_B_a = (/(i, i=0, nproc - 1)/)
      dist_L_R = (/(MODULO(i, nproc), i=0, ngroup - 1)/) ! R index is replicated in my_Lrows, we impose a cyclic distribution
      dist_L_P = (/(i, i=0, ngroup - 1)/)

      ! create distributions and tensors
      CALL dbcsr_t_distribution_new(dist_B, mp_comm_B, dist_L_P, dist_B_a, dist_B_i)
      CALL dbcsr_t_distribution_new(dist_L, mp_comm_L, dist_L_R, dist_L_P)

      CALL dbcsr_t_create(tB_in, "(R|ai)", dist_B, map_B_1, map_B_2, dbcsr_type_real_8, sizes_L, sizes_B, sizes_i)
      CALL dbcsr_t_create(tB_out, "(P|ai)", dist_B, map_B_1, map_B_2, dbcsr_type_real_8, sizes_L, sizes_B, sizes_i)
      CALL dbcsr_t_create(tL, "(R|P)", dist_L, map_L_1, map_L_2, dbcsr_type_real_8, sizes_L, sizes_L)

      IF (debug) THEN
         ! check that tensor distribution is correct
         CALL dbcsr_t_get_stored_coordinates(tB_in, [loc_P, loc_a, 1], check_proc)
         CPASSERT(check_proc .EQ. iproc_glob)
      ENDIF

      ! reserve (R|ai) block
      CALL dbcsr_t_reserve_blocks(tB_in, [loc_P], [loc_a], [1])

      ! reserve (R|P) blocks
      ! in my_Lrows, R index is replicated. For (R|P), we distribute quadratic blocks cyclically over
      ! the processes in a subgroup.
      ! There are NG blocks, so each process holds at most NG/Nw+1 blocks.
      ALLOCATE (block_ind_L_R(ngroup/nproc + 1))
      ALLOCATE (block_ind_L_P(ngroup/nproc + 1))
      block_ind_L_R(:) = 0; block_ind_L_P(:) = 0
      ii = 0
      DO i = 1, ngroup
         CALL dbcsr_t_get_stored_coordinates(tL, [i, loc_P], check_proc)
         IF (check_proc == iproc_glob) THEN
            ii = ii + 1
            block_ind_L_R(ii) = i
            block_ind_L_P(ii) = loc_P
         ENDIF
      ENDDO
      CALL dbcsr_t_reserve_blocks(tL, block_ind_L_R(1:ii), block_ind_L_P(1:ii))

      ! insert (R|ai) block
      CALL dbcsr_t_put_block(tB_in, [loc_P, loc_a, 1], SHAPE(BIb_C), BIb_C)

      ! insert (R|P) blocks
      ioff = 0
      DO i = 1, ngroup
         istart = ioff + 1; iend = ioff + sizes_L(i)
         ioff = ioff + sizes_L(i)
         CALL dbcsr_t_get_stored_coordinates(tL, [i, loc_P], check_proc)
         IF (check_proc == iproc_glob) THEN
            CALL dbcsr_t_put_block(tL, [i, loc_P], [sizes_L(i), sizes_L(loc_P)], my_Lrows(istart:iend, :))
         ENDIF
      ENDDO

      CALL dbcsr_t_split_blocks(tB_in, tB_in_split, [blk_size(2), blk_size(1), blk_size(1)])
      CALL dbcsr_t_split_blocks(tL, tL_split, [blk_size(2), blk_size(2)])
      CALL dbcsr_t_split_blocks(tB_out, tB_out_split, [blk_size(2), blk_size(1), blk_size(1)])

      ! contract
      CALL dbcsr_t_contract(alpha=dbcsr_scalar(1.0_dp), tensor_1=tB_in_split, tensor_2=tL_split, &
                            beta=dbcsr_scalar(0.0_dp), tensor_3=tB_out_split, &
                            contract_1=[1], notcontract_1=[2, 3], &
                            contract_2=[1], notcontract_2=[2], &
                            map_1=[2, 3], map_2=[1], optimize_dist=.TRUE.)

      ! retrieve local block of contraction result (P|ai)
      CALL dbcsr_t_copy(tB_out_split, tB_out)

      CALL dbcsr_t_get_block(tB_out, [loc_P, loc_a, 1], SHAPE(BIb_C), BIb_C, found)
      CPASSERT(found)

      ! cleanup
      CALL dbcsr_t_destroy(tB_in)
      CALL dbcsr_t_destroy(tB_in_split)
      CALL dbcsr_t_destroy(tB_out)
      CALL dbcsr_t_destroy(tB_out_split)
      CALL dbcsr_t_destroy(tL)
      CALL dbcsr_t_destroy(tL_split)

      CALL dbcsr_t_distribution_destroy(dist_B)
      CALL dbcsr_t_distribution_destroy(dist_L)

      CALL dbcsr_t_pgrid_destroy(mp_comm_B)
      CALL dbcsr_t_pgrid_destroy(mp_comm_L)

      CALL timestop(handle)

   END SUBROUTINE contract_B_L

! **************************************************************************************************
!> \brief Encapsulate building of intermediate matrices matrix_ia_jnu(_beta
!>         matrix_ia_jb(_beta),fm_BIb_jb(_beta),matrix_in_jnu(for G0W0) and
!>         fm_BIb_all(for G0W0)
!> \param matrix_ia_jnu ...
!> \param matrix_ia_jb ...
!> \param mo_coeff_templ ...
!> \param size_1 ...
!> \param size_2 ...
!> \param fm_BIb_jb ...
!> \param matrix_name_2 ...
!> \param max_row_col_local ...
!> \param blacs_env_sub ...
!> \param para_env_sub ...
!> \param local_col_row_info ...
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE create_intermediate_matrices(matrix_ia_jnu, matrix_ia_jb, mo_coeff_templ, size_1, size_2, &
                                           fm_BIb_jb, matrix_name_2, max_row_col_local, &
                                           blacs_env_sub, para_env_sub, local_col_row_info)

      TYPE(dbcsr_type), INTENT(OUT)                      :: matrix_ia_jnu, matrix_ia_jb
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mo_coeff_templ
      INTEGER, INTENT(IN)                                :: size_1, size_2
      TYPE(cp_fm_type), POINTER                          :: fm_BIb_jb
      CHARACTER(LEN=*), INTENT(IN)                       :: matrix_name_2
      INTEGER, INTENT(OUT)                               :: max_row_col_local
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: local_col_row_info

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_intermediate_matrices'

      INTEGER                                            :: handle, ncol_local, nfullcols_total, &
                                                            nfullrows_total, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      ! initialize and create the matrix (K|jnu)
      CALL dbcsr_create(matrix_ia_jnu, template=mo_coeff_templ)

      ! Allocate Sparse matrices: (K|jb)
      CALL cp_dbcsr_m_by_n_from_template(matrix_ia_jb, template=mo_coeff_templ, m=size_2, n=size_1, &
                                         sym=dbcsr_type_no_symmetry)

      ! set all to zero in such a way that the memory is actually allocated
      CALL dbcsr_set(matrix_ia_jnu, 0.0_dp)
      CALL dbcsr_set(matrix_ia_jb, 0.0_dp)

      ! create the analogous of matrix_ia_jb in fm type
      NULLIFY (fm_BIb_jb)
      NULLIFY (fm_struct)
      CALL dbcsr_get_info(matrix_ia_jb, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env_sub, nrow_global=nfullrows_total, &
                               ncol_global=nfullcols_total, para_env=para_env_sub)
      CALL cp_fm_create(fm_BIb_jb, fm_struct, name=matrix_name_2)

      CALL copy_dbcsr_to_fm(matrix_ia_jb, fm_BIb_jb)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_get_info(matrix=fm_BIb_jb, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      max_row_col_local = MAX(nrow_local, ncol_local)
      CALL mp_max(max_row_col_local, para_env_sub%group)

      ALLOCATE (local_col_row_info(0:max_row_col_local, 2))
      local_col_row_info = 0
      ! 0,1 nrows
      local_col_row_info(0, 1) = nrow_local
      local_col_row_info(1:nrow_local, 1) = row_indices(1:nrow_local)
      ! 0,2 ncols
      local_col_row_info(0, 2) = ncol_local
      local_col_row_info(1:ncol_local, 2) = col_indices(1:ncol_local)

      CALL timestop(handle)

   END SUBROUTINE create_intermediate_matrices

! **************************************************************************************************
!> \brief Encapsulate ERI postprocessing: AO to MO transformation and store in B matrix.
!> \param para_env ...
!> \param mat_munu ...
!> \param matrix_ia_jnu ...
!> \param matrix_ia_jb ...
!> \param fm_BIb_jb ...
!> \param BIb_jb ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param eps_filter ...
!> \param max_row_col_local ...
!> \param sub_proc_map ...
!> \param local_col_row_info ...
!> \param my_B_end ...
!> \param my_B_start ...
!> \param descr ...
! **************************************************************************************************
   SUBROUTINE ao_to_mo_and_store_B(para_env, mat_munu, matrix_ia_jnu, matrix_ia_jb, fm_BIb_jb, BIb_jb, &
                                   mo_coeff_o, mo_coeff_v, eps_filter, &
                                   max_row_col_local, sub_proc_map, local_col_row_info, &
                                   my_B_end, my_B_start, descr)
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_ia_jnu, matrix_ia_jb
      TYPE(cp_fm_type), POINTER                          :: fm_BIb_jb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: BIb_jb
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_o, mo_coeff_v
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      INTEGER, INTENT(IN)                                :: max_row_col_local
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: sub_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(IN)  :: local_col_row_info
      INTEGER, INTENT(IN)                                :: my_B_end, my_B_start
      CHARACTER(LEN=*), INTENT(IN)                       :: descr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ao_to_mo_and_store_B'

      INTEGER                                            :: handle1, handle2

      CALL timeset(routineN//"_mult_"//descr, handle1)

      CALL dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_o, &
                          0.0_dp, matrix_ia_jnu, filter_eps=eps_filter)
      CALL dbcsr_multiply("T", "N", 1.0_dp, matrix_ia_jnu, mo_coeff_v, &
                          0.0_dp, matrix_ia_jb, filter_eps=eps_filter)
      CALL timestop(handle1)

      CALL timeset(routineN//"_E_Ex_"//descr, handle2)
      CALL copy_dbcsr_to_fm(matrix_ia_jb, fm_BIb_jb)

      IF (.NOT. (descr .EQ. "bse_ab")) THEN

         CALL grep_my_integrals(para_env, fm_BIb_jb, BIb_jb, max_row_col_local, &
                                sub_proc_map, local_col_row_info, &
                                my_B_end, my_B_start)

      END IF

      CALL timestop(handle2)
   END SUBROUTINE ao_to_mo_and_store_B

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kpoints ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE compute_kpoints(qs_env, kpoints, unit_nr)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_kpoints'

      INTEGER                                            :: handle, i, i_dim, ix, iy, iz, nkp, &
                                                            num_dim
      INTEGER, DIMENSION(3)                              :: nkp_grid, periodic
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      CALL timeset(routineN, handle)

      NULLIFY (cell, dft_control, para_env)
      CALL get_qs_env(qs_env=qs_env, cell=cell, para_env=para_env, dft_control=dft_control, sab_orb=sab_orb)
      CALL get_cell(cell=cell, periodic=periodic)

      ! general because we augment a Monkhorst-Pack mesh by additional points in the BZ
      kpoints%kp_scheme = "GENERAL"
      kpoints%symmetry = .FALSE.
      kpoints%verbose = .FALSE.
      kpoints%full_grid = .FALSE.
      kpoints%use_real_wfn = .FALSE.
      kpoints%eps_geo = 1.e-6_dp
      kpoints%full_grid = .TRUE.
      nkp_grid(1:3) = qs_env%mp2_env%ri_rpa_im_time%kp_grid(1:3)
      kpoints%nkp_grid(1:3) = nkp_grid(1:3)

      num_dim = periodic(1) + periodic(2) + periodic(3)

      DO i_dim = 1, 3
         IF (periodic(i_dim) == 1) THEN
            CPASSERT(MODULO(kpoints%nkp_grid(i_dim), 2) == 0)
         END IF
      END DO

      IF (num_dim == 3) THEN
         nkp = nkp_grid(1)*nkp_grid(2)*nkp_grid(3)/2
      ELSE IF (num_dim == 2) THEN
         nkp = MAX(nkp_grid(1)*nkp_grid(2)/2, nkp_grid(1)*nkp_grid(3)/2, nkp_grid(2)*nkp_grid(3)/2)
      ELSE IF (num_dim == 1) THEN
         nkp = MAX(nkp_grid(1)/2, nkp_grid(2)/2, nkp_grid(3)/2)
      END IF

      kpoints%nkp = nkp

      IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
         "KPOINT_INFO| Number of kpoints for V and W:", nkp

      ALLOCATE (kpoints%xkp(3, nkp), kpoints%wkp(nkp))
      kpoints%wkp(:) = 1.0_dp/REAL(nkp, KIND=dp)

      i = 0
      DO ix = 1, nkp_grid(1)
         DO iy = 1, nkp_grid(2)
            DO iz = 1, nkp_grid(3)

               i = i + 1
               IF (i > nkp) CYCLE

               kpoints%xkp(1, i) = REAL(2*ix - nkp_grid(1) - 1, KIND=dp)/(2._dp*REAL(nkp_grid(1), KIND=dp))
               kpoints%xkp(2, i) = REAL(2*iy - nkp_grid(2) - 1, KIND=dp)/(2._dp*REAL(nkp_grid(2), KIND=dp))
               kpoints%xkp(3, i) = REAL(2*iz - nkp_grid(3) - 1, KIND=dp)/(2._dp*REAL(nkp_grid(3), KIND=dp))

            END DO
         END DO
      END DO

      CALL kpoint_init_cell_index(kpoints, sab_orb, para_env, dft_control)

      CALL set_qs_env(qs_env, kpoints=kpoints)

      CALL timestop(handle)

   END SUBROUTINE compute_kpoints

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param dimen_RI ...
!> \param fm_matrix_L ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param my_group_L_size ...
!> \param my_Lrows ...
! **************************************************************************************************
   SUBROUTINE grep_Lcols(para_env, dimen_RI, fm_matrix_L, &
                         my_group_L_start, my_group_L_end, my_group_L_size, my_Lrows)
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: dimen_RI
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L
      INTEGER, INTENT(IN)                                :: my_group_L_start, my_group_L_end, &
                                                            my_group_L_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: my_Lrows

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'grep_Lcols'

      INTEGER :: handle, handle2, i_global, iiB, j_global, jjB, max_row_col_local, ncol_local, &
         ncol_rec, nrow_local, nrow_rec, proc_receive, proc_receive_static, proc_send, &
         proc_send_static, proc_shift
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_col_row_info, rec_col_row_info
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, col_indices_rec, &
                                                            row_indices, row_indices_rec
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: local_L, local_L_internal, rec_L

      CALL timeset(routineN, handle)

      ALLOCATE (my_Lrows(dimen_RI, my_group_L_size))
      my_Lrows = 0.0_dp

      ! proc_map, vector that replicates the processor numbers also
      ! for negative and positive number > num_pe
      ! needed to know which is the processor, to respect to another one,
      ! for a given shift
      ALLOCATE (proc_map(-para_env%num_pe:2*para_env%num_pe - 1))
      DO iiB = 0, para_env%num_pe - 1
         proc_map(iiB) = iiB
         proc_map(-iiB - 1) = para_env%num_pe - iiB - 1
         proc_map(para_env%num_pe + iiB) = iiB
      END DO

      CALL cp_fm_get_info(matrix=fm_matrix_L, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices, &
                          local_data=local_L_internal)

      ALLOCATE (local_L(nrow_local, ncol_local))
      local_L = local_L_internal(1:nrow_local, 1:ncol_local)

      max_row_col_local = MAX(nrow_local, ncol_local)
      CALL mp_max(max_row_col_local, para_env%group)

      ALLOCATE (local_col_row_info(0:max_row_col_local, 2))
      local_col_row_info = 0
      ! 0,1 nrows
      local_col_row_info(0, 1) = nrow_local
      local_col_row_info(1:nrow_local, 1) = row_indices(1:nrow_local)
      ! 0,2 ncols
      local_col_row_info(0, 2) = ncol_local
      local_col_row_info(1:ncol_local, 2) = col_indices(1:ncol_local)

      ALLOCATE (rec_col_row_info(0:max_row_col_local, 2))

      ! accumulate data on my_Lrows starting from myself
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         IF (j_global >= my_group_L_start .AND. j_global <= my_group_L_end) THEN
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)
               my_Lrows(i_global, j_global - my_group_L_start + 1) = local_L(iiB, jjB)
            END DO
         END IF
      END DO

      proc_send_static = proc_map(para_env%mepos + 1)
      proc_receive_static = proc_map(para_env%mepos - 1)

      CALL timeset(routineN//"_comm", handle2)

      DO proc_shift = 1, para_env%num_pe - 1
         proc_send = proc_map(para_env%mepos + proc_shift)
         proc_receive = proc_map(para_env%mepos - proc_shift)

         ! first exchange information on the local data
         rec_col_row_info = 0
         CALL mp_sendrecv(local_col_row_info, proc_send_static, rec_col_row_info, proc_receive_static, para_env%group)
         nrow_rec = rec_col_row_info(0, 1)
         ncol_rec = rec_col_row_info(0, 2)

         ALLOCATE (row_indices_rec(nrow_rec))
         row_indices_rec = rec_col_row_info(1:nrow_rec, 1)

         ALLOCATE (col_indices_rec(ncol_rec))
         col_indices_rec = rec_col_row_info(1:ncol_rec, 2)

         ALLOCATE (rec_L(nrow_rec, ncol_rec))
         rec_L = 0.0_dp

         ! then send and receive the real data
         CALL mp_sendrecv(local_L, proc_send_static, rec_L, proc_receive_static, para_env%group)

         ! accumulate the received data on my_Lrows
         DO jjB = 1, ncol_rec
            j_global = col_indices_rec(jjB)
            IF (j_global >= my_group_L_start .AND. j_global <= my_group_L_end) THEN
               DO iiB = 1, nrow_rec
                  i_global = row_indices_rec(iiB)
                  my_Lrows(i_global, j_global - my_group_L_start + 1) = rec_L(iiB, jjB)
               END DO
            END IF
         END DO

         local_col_row_info(:, :) = rec_col_row_info
         DEALLOCATE (local_L)
         ALLOCATE (local_L(nrow_rec, ncol_rec))
         local_L = rec_L

         DEALLOCATE (col_indices_rec)
         DEALLOCATE (row_indices_rec)
         DEALLOCATE (rec_L)
      END DO
      CALL timestop(handle2)

      DEALLOCATE (local_col_row_info)
      DEALLOCATE (rec_col_row_info)
      DEALLOCATE (proc_map)
      DEALLOCATE (local_L)

      CALL timestop(handle)

   END SUBROUTINE grep_Lcols

! **************************************************************************************************
!> \brief ...
!> \param para_env_sub ...
!> \param fm_BIb_jb ...
!> \param BIb_jb ...
!> \param max_row_col_local ...
!> \param proc_map ...
!> \param local_col_row_info ...
!> \param my_B_virtual_end ...
!> \param my_B_virtual_start ...
! **************************************************************************************************
   SUBROUTINE grep_my_integrals(para_env_sub, fm_BIb_jb, BIb_jb, max_row_col_local, &
                                proc_map, local_col_row_info, &
                                my_B_virtual_end, my_B_virtual_start)
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      TYPE(cp_fm_type), POINTER                          :: fm_BIb_jb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: BIb_jb
      INTEGER, INTENT(IN)                                :: max_row_col_local
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(IN)  :: local_col_row_info
      INTEGER, INTENT(IN)                                :: my_B_virtual_end, my_B_virtual_start

      INTEGER                                            :: i_global, iiB, j_global, jjB, ncol_rec, &
                                                            nrow_rec, proc_receive, proc_send, &
                                                            proc_shift
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: rec_col_row_info
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_rec, row_indices_rec
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: local_BI, rec_BI

      ALLOCATE (rec_col_row_info(0:max_row_col_local, 2))

      rec_col_row_info(:, :) = local_col_row_info

      nrow_rec = rec_col_row_info(0, 1)
      ncol_rec = rec_col_row_info(0, 2)

      ALLOCATE (row_indices_rec(nrow_rec))
      row_indices_rec = rec_col_row_info(1:nrow_rec, 1)

      ALLOCATE (col_indices_rec(ncol_rec))
      col_indices_rec = rec_col_row_info(1:ncol_rec, 2)

      ! accumulate data on BIb_jb buffer starting from myself
      DO jjB = 1, ncol_rec
         j_global = col_indices_rec(jjB)
         IF (j_global >= my_B_virtual_start .AND. j_global <= my_B_virtual_end) THEN
            DO iiB = 1, nrow_rec
               i_global = row_indices_rec(iiB)
               BIb_jb(j_global - my_B_virtual_start + 1, i_global) = fm_BIb_jb%local_data(iiB, jjB)
            END DO
         END IF
      END DO

      DEALLOCATE (row_indices_rec)
      DEALLOCATE (col_indices_rec)

      IF (para_env_sub%num_pe > 1) THEN
         ALLOCATE (local_BI(nrow_rec, ncol_rec))
         local_BI(1:nrow_rec, 1:ncol_rec) = fm_BIb_jb%local_data(1:nrow_rec, 1:ncol_rec)

         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_send = proc_map(para_env_sub%mepos + proc_shift)
            proc_receive = proc_map(para_env_sub%mepos - proc_shift)

            ! first exchange information on the local data
            rec_col_row_info = 0
            CALL mp_sendrecv(local_col_row_info, proc_send, rec_col_row_info, proc_receive, para_env_sub%group)
            nrow_rec = rec_col_row_info(0, 1)
            ncol_rec = rec_col_row_info(0, 2)

            ALLOCATE (row_indices_rec(nrow_rec))
            row_indices_rec = rec_col_row_info(1:nrow_rec, 1)

            ALLOCATE (col_indices_rec(ncol_rec))
            col_indices_rec = rec_col_row_info(1:ncol_rec, 2)

            ALLOCATE (rec_BI(nrow_rec, ncol_rec))
            rec_BI = 0.0_dp

            ! then send and receive the real data
            CALL mp_sendrecv(local_BI, proc_send, rec_BI, proc_receive, para_env_sub%group)

            ! accumulate the received data on BIb_jb buffer
            DO jjB = 1, ncol_rec
               j_global = col_indices_rec(jjB)
               IF (j_global >= my_B_virtual_start .AND. j_global <= my_B_virtual_end) THEN
                  DO iiB = 1, nrow_rec
                     i_global = row_indices_rec(iiB)
                     BIb_jb(j_global - my_B_virtual_start + 1, i_global) = rec_BI(iiB, jjB)
                  END DO
               END IF
            END DO

            DEALLOCATE (col_indices_rec)
            DEALLOCATE (row_indices_rec)
            DEALLOCATE (rec_BI)
         END DO

         DEALLOCATE (local_BI)
      END IF

      DEALLOCATE (rec_col_row_info)

   END SUBROUTINE grep_my_integrals

END MODULE mp2_integrals
