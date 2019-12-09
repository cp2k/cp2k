!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate and distribute 2c- and 3c- integrals for RI
!> \par History
!>      06.2012 created [Mauro Del Ben]
!>      03.2019 separated from mp2_ri_gpw [Frederick Stein]
! **************************************************************************************************
MODULE mp2_integrals
   USE ai_contraction_sphi,             ONLY: ab_contract,&
                                              abc_contract
   USE ai_overlap3,                     ONLY: overlap3
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: DelBen2013,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_m_by_n_from_template,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
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
   USE dbcsr_api,                       ONLY: &
        dbcsr_complete_redistribute, dbcsr_copy, dbcsr_create, dbcsr_distribution_type, &
        dbcsr_filter, dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_get_stored_coordinates, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_release_p, dbcsr_reserve_blocks, &
        dbcsr_scalar, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_real_8
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_create, dbcsr_t_destroy, &
        dbcsr_t_distribution_destroy, dbcsr_t_distribution_new, dbcsr_t_distribution_type, &
        dbcsr_t_filter, dbcsr_t_get_block, dbcsr_t_get_stored_coordinates, &
        dbcsr_t_mp_environ_pgrid, dbcsr_t_pgrid_create, dbcsr_t_pgrid_destroy, dbcsr_t_pgrid_type, &
        dbcsr_t_put_block, dbcsr_t_reserve_blocks, dbcsr_t_split_blocks, dbcsr_t_type
   USE group_dist_types,                ONLY: create_group_dist,&
                                              get_group_dist,&
                                              group_dist_d1_type,&
                                              release_group_dist
   USE input_constants,                 ONLY: do_eri_gpw,&
                                              do_eri_mme,&
                                              do_eri_os,&
                                              do_potential_coulomb,&
                                              do_potential_id,&
                                              do_potential_long,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE machine,                         ONLY: m_flush
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_cart_create,&
                                              mp_dims_create,&
                                              mp_environ,&
                                              mp_max,&
                                              mp_sendrecv,&
                                              mp_sync
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE mp2_eri,                         ONLY: mp2_eri_3c_integrate
   USE mp2_eri_gpw,                     ONLY: cleanup_gpw,&
                                              mp2_eri_3c_integrate_gpw,&
                                              prepare_gpw
   USE mp2_ri_2c,                       ONLY: get_2c_integrals
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              one_dim_int_array,&
                                              two_dim_real_array
   USE orbital_pointers,                ONLY: ncoset
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
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_tensors,                      ONLY: build_3c_integrals,&
                                              build_3c_neighbor_lists,&
                                              contiguous_tensor_dist,&
                                              neighbor_list_3c_destroy
   USE qs_tensors_types,                ONLY: cyclic_tensor_dist,&
                                              distribution_3d_create,&
                                              distribution_3d_type,&
                                              neighbor_list_3c_type,&
                                              split_block_sizes
   USE rpa_communication,               ONLY: communicate_buffer
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
!> \param mat_munu_mao_occ_virt ...
!> \param mat_munu_mao_virt_occ ...
!> \param sab_orb_sub ...
!> \param sab_orb_all ...
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
!> \param starts_array_mc_t ...
!> \param ends_array_mc_t ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param do_im_time ...
!> \param do_mao ...
!> \param do_kpoints_cubic_RPA ...
!> \param kpoints ...
!> \param mat_3c_overl_int ...
!> \param do_dbcsr_t ...
!> \param t_3c_overl_int ...
!> \param t_3c_M ...
!> \param t_3c_O ...
!> \param mat_3c_overl_int_mao_for_occ ...
!> \param mat_3c_overl_int_mao_for_virt ...
!> \param mao_coeff_occ ...
!> \param mao_coeff_virt ...
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
                                    mat_munu, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, &
                                    sab_orb_sub, sab_orb_all, mo_coeff_o, mo_coeff_v, mo_coeff_all, &
                                    mo_coeff_gw, eps_filter, unit_nr, &
                                    mp2_memory, calc_PQ_cond_num, calc_forces, blacs_env_sub, my_do_gw, do_bse, &
                                    gd_B_all, starts_array_mc_t, ends_array_mc_t, &
                                    gw_corr_lev_occ, gw_corr_lev_virt, &
                                    do_im_time, do_mao, do_kpoints_cubic_RPA, kpoints, &
                                    mat_3c_overl_int, do_dbcsr_t, t_3c_overl_int, t_3c_M, t_3c_O, mat_3c_overl_int_mao_for_occ, &
                                    mat_3c_overl_int_mao_for_virt, mao_coeff_occ, mao_coeff_virt, &
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
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_munu, mat_munu_mao_occ_virt, &
                                                            mat_munu_mao_virt_occ
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: sab_orb_sub, sab_orb_all
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_o, mo_coeff_v, mo_coeff_all, &
                                                            mo_coeff_gw
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), INTENT(IN)                          :: mp2_memory
      LOGICAL, INTENT(IN)                                :: calc_PQ_cond_num, calc_forces
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      LOGICAL, INTENT(IN)                                :: my_do_gw, do_bse
      TYPE(group_dist_d1_type), INTENT(OUT)              :: gd_B_all
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: starts_array_mc_t, ends_array_mc_t
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt
      LOGICAL, INTENT(IN)                                :: do_im_time, do_mao, do_kpoints_cubic_RPA
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int
      LOGICAL, INTENT(IN)                                :: do_dbcsr_t
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:, :)   :: t_3c_overl_int
      TYPE(dbcsr_t_type), INTENT(OUT)                    :: t_3c_M
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: t_3c_O
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int_mao_for_occ, &
                                                            mat_3c_overl_int_mao_for_virt
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coeff_occ, mao_coeff_virt
      TYPE(libint_potential_type), INTENT(IN)            :: ri_metric
      TYPE(group_dist_d1_type), INTENT(OUT)              :: gd_B_occ_bse, gd_B_virt_bse
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT), OPTIONAL                           :: BIb_C_beta, BIb_C_gw_beta
      TYPE(group_dist_d1_type), INTENT(OUT), OPTIONAL    :: gd_B_virtual_beta
      INTEGER, INTENT(IN), OPTIONAL                      :: homo_beta
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: mo_coeff_o_beta, mo_coeff_v_beta, &
                                                            mo_coeff_all_beta, mo_coeff_gw_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_gpw_compute_in', &
         routineP = moduleN//':'//routineN

      INTEGER :: color_sub_3c, cut_memory, cut_RI, eri_method, gw_corr_lev_total, handle, handle2, &
         handle4, i, i_counter, itmp(2), j, LLL, max_row_col_local, max_row_col_local_beta, &
         max_row_col_local_gw, max_row_col_local_occ_bse, max_row_col_local_virt_bse, &
         my_B_all_end, my_B_all_size, my_B_all_start, my_B_occ_bse_end, my_B_occ_bse_size, &
         my_B_occ_bse_start, my_B_size, my_B_size_beta, my_B_virt_bse_end, my_B_virt_bse_size, &
         my_B_virt_bse_start, my_B_virtual_end, my_B_virtual_end_beta, my_B_virtual_start, &
         my_B_virtual_start_beta, my_group_L_end, my_group_L_size, my_group_L_start, natom, ngroup
      INTEGER :: ngroup_im_time, ngroup_im_time_P, nimg, nkind, num_diff_blk, num_small_eigen, &
         potential_type, ri_metric_type, sqrt_max_bsize, virtual, virtual_beta
      INTEGER, ALLOCATABLE, DIMENSION(:) :: local_atoms_for_mao_basis, my_group_L_ends_im_time, &
         my_group_L_sizes_im_time, my_group_L_starts_im_time, sizes_AO, sizes_AO_split, sizes_RI, &
         sizes_RI_split, sub_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: local_col_row_info, local_col_row_info_beta, &
         local_col_row_info_gw, local_col_row_info_occ_bse, local_col_row_info_virt_bse
      INTEGER, DIMENSION(3)                              :: pdims_t3c, periodic
      LOGICAL                                            :: do_alpha_beta, do_gpw, &
                                                            do_kpoints_from_Gamma, do_svd, &
                                                            impose_split, memory_info
      REAL(KIND=dp)                                      :: cond_num, cutoff_old, eps_svd, &
                                                            mem_for_iaK, omega_pot, &
                                                            relative_cutoff_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: e_cutoff_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: my_Lrows
      TYPE(cp_eri_mme_param), POINTER                    :: eri_param
      TYPE(cp_fm_type), POINTER                          :: fm_BIb_bse_ab, fm_BIb_bse_ij, fm_BIb_gw, &
                                                            fm_BIb_gw_beta, fm_BIb_jb, &
                                                            fm_BIb_jb_beta, fm_matrix_L
      TYPE(cp_para_env_type), POINTER                    :: para_env_L
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_munu_local_L
      TYPE(dbcsr_t_pgrid_type)                           :: pgrid_t3c_M, pgrid_t3c_overl
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_template
      TYPE(dbcsr_type) :: matrix_bse_ab, matrix_bse_anu, matrix_bse_ij, matrix_bse_inu, &
         matrix_ia_jb, matrix_ia_jb_beta, matrix_ia_jnu, matrix_ia_jnu_beta, matrix_in_jm, &
         matrix_in_jm_beta, matrix_in_jnu, matrix_in_jnu_beta
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_ao, basis_set_ri_aux
      TYPE(neighbor_list_3c_type)                        :: nl_3c
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_p_type)                                    :: pot_g, psi_L, rho_g, rho_r
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(task_list_type), POINTER                      :: task_list_sub
      TYPE(two_dim_real_array), ALLOCATABLE, &
         DIMENSION(:)                                    :: mao_coeff_repl_occ, mao_coeff_repl_virt

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
                            atomic_kind_set, qs_kind_set, sab_orb_sub)

      IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
         "RI_INFO| Cholesky decomposition group size:", para_env_L%num_pe
      IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
         "RI_INFO| Number of groups for auxiliary basis functions", ngroup
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
            "RI_INFO| Occupied  basis set size:", homo, &
            "RI_INFO| Virtual   basis set size:", virtual, &
            "RI_INFO| Auxiliary basis set size:", dimen_RI
         IF (do_svd) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "RI_INFO| Reduced auxiliary basis set size:", dimen_RI_red
         END IF

         mem_for_iaK = dimen_RI*REAL(homo, KIND=dp)*virtual*8.0_dp/(1024_dp**2)
         IF (do_alpha_beta) mem_for_iaK = mem_for_iaK + &
                                          dimen_RI*REAL(homo_beta, KIND=dp)*(nmo - homo_beta)*8.0_dp/(1024_dp**2)

         WRITE (unit_nr, '(T3,A,T66,F11.2,A4)') 'RI_INFO| Total memory for (ia|K) integrals:', &
            mem_for_iaK, ' MiB'
         IF (my_do_gw .AND. .NOT. do_im_time) THEN
            mem_for_iaK = dimen_RI*REAL(nmo, KIND=dp)*gw_corr_lev_total*8.0_dp/(1024_dp**2)

            WRITE (unit_nr, '(T3,A,T66,F11.2,A4)') 'RI_INFO| Total memory for G0W0-(nm|K) integrals:', &
               mem_for_iaK, ' MiB'
         END IF
         CALL m_flush(unit_nr)
      ENDIF

      CALL mp_sync(para_env%group) ! sync to see memory output

      ! in case we do imaginary time, we need the overlap matrix (alpha beta P)
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

      ELSE

         cut_memory = qs_env%mp2_env%ri_rpa_im_time%cut_memory

         ngroup_im_time = para_env%num_pe/qs_env%mp2_env%ri_rpa_im_time%group_size_3c
         ngroup_im_time_P = para_env%num_pe/qs_env%mp2_env%ri_rpa_im_time%group_size_p
         impose_split = .NOT. qs_env%mp2_env%ri_rpa_im_time%group_size_internal

         color_sub_3c = para_env%mepos/qs_env%mp2_env%ri_rpa_im_time%group_size_3c

         CALL setup_group_L_im_time(my_group_L_starts_im_time, my_group_L_ends_im_time, my_group_L_sizes_im_time, &
                                    dimen_RI, ngroup_im_time, cut_memory, &
                                    cut_RI, &
                                    color_sub_3c, qs_env)

         memory_info = qs_env%mp2_env%ri_rpa_im_time%memory_info

         IF (do_dbcsr_t) THEN
            ! we need 3 tensors:
            ! 1) t_3c_overl_int: 3c overlap integrals, optimized for easy access to integral blocks
            !                   (atomic blocks)
            ! 2) t_3c_O: 3c overlap integrals, optimized for contraction (split blocks)
            ! 3) t_3c_M: tensor M, optimized for contraction

            CALL get_qs_env(qs_env, natom=natom, nkind=nkind, dft_control=dft_control)

            IF (.NOT. impose_split) THEN
               pdims_t3c = 0
               CALL dbcsr_t_pgrid_create(para_env%group, pdims_t3c, pgrid_t3c_overl, map1_2d=[1, 2], map2_2d=[3])
            ELSE
               pgrid_t3c_overl = get_pgrid_from_ngroup(para_env, ngroup_im_time, map1_2d=[1, 2], map2_2d=[3])
            ENDIF

            ! set up basis
            ALLOCATE (sizes_RI(natom), sizes_AO(natom))
            ALLOCATE (basis_set_ri_aux(nkind), basis_set_ao(nkind))
            CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
            CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=basis_set_ri_aux)
            CALL basis_set_list_setup(basis_set_ao, "ORB", qs_kind_set)
            CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_AO, basis=basis_set_ao)

            CALL create_tensor_O_3c(t_3c_overl_int_template, pgrid_t3c_overl, &
                                    sizes_AO, sizes_AO, sizes_RI, &
                                    ri_metric=ri_metric, nl_3c=nl_3c, qs_env=qs_env, &
                                    basis_ao_1=basis_set_ao, basis_ao_2=basis_set_ao, basis_RI=basis_set_ri_aux, &
                                    do_kpoints=do_kpoints_cubic_RPA)

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
            sqrt_max_bsize = qs_env%mp2_env%ri_rpa_im_time%max_bsize_sqrt

            CALL split_block_sizes(sizes_AO, sizes_AO_split, sqrt_max_bsize)
            CALL split_block_sizes(sizes_RI, sizes_RI_split, sqrt_max_bsize)

            IF (.NOT. impose_split) THEN
               pdims_t3c = 0
               CALL dbcsr_t_pgrid_create(para_env%group, pdims_t3c, pgrid_t3c_M, map1_2d=[1], map2_2d=[2, 3])
            ELSE
               pgrid_t3c_M = get_pgrid_from_ngroup(para_env, ngroup_im_time_P, map1_2d=[1], map2_2d=[2, 3])
            ENDIF

            CALL create_tensor_M_3c(t_3c_M, pgrid_t3c_M, &
                                    cut_memory, &
                                    sizes_AO_split, sizes_RI, starts_array_mc_t, ends_array_mc_t)

            CALL dbcsr_t_pgrid_destroy(pgrid_t3c_M)

            ALLOCATE (t_3c_O(SIZE(t_3c_overl_int, 1), SIZE(t_3c_overl_int, 2)))
            CALL create_tensor_O_3c(t_3c_O(1, 1), pgrid_t3c_overl, sizes_AO_split, &
                                    sizes_AO, sizes_RI_split)

            CALL dbcsr_t_pgrid_destroy(pgrid_t3c_overl)
            DO i = 1, SIZE(t_3c_O, 1)
               DO j = 1, SIZE(t_3c_O, 2)
                  IF (i > 1 .OR. j > 1) CALL dbcsr_t_create(t_3c_O(1, 1), t_3c_O(i, j))
               ENDDO
            ENDDO

            CALL build_3c_integrals(t_3c_overl_int, &
                                    qs_env%mp2_env%mp2_gpw%eps_filter, &
                                    qs_env, &
                                    nl_3c, &
                                    basis_i=basis_set_ri_aux, &
                                    basis_j=basis_set_ao, basis_k=basis_set_ao, &
                                    potential_parameter=ri_metric, &
                                    do_kpoints=do_kpoints_cubic_RPA)

            DEALLOCATE (basis_set_ri_aux, basis_set_ao)
            ! copy integral tensor t_3c_overl_int to t_3c_O tensor optimized for contraction
            CALL timeset(routineN//"_copy_3c", handle4)
            DO i = 1, SIZE(t_3c_overl_int, 1)
               DO j = 1, SIZE(t_3c_overl_int, 2)

                  CALL dbcsr_t_copy(t_3c_overl_int(i, j), t_3c_O(i, j), order=[1, 3, 2], move_data=.FALSE.)
                  CALL dbcsr_t_filter(t_3c_O(i, j), qs_env%mp2_env%mp2_gpw%eps_filter)
               ENDDO
            ENDDO
            CALL timestop(handle4)

            CALL neighbor_list_3c_destroy(nl_3c)
         ELSE
            CALL reserve_blocks_3c(mat_3c_overl_int, &
                                   mat_munu, qs_env, &
                                   dimen_RI, cut_RI, unit_nr, &
                                   my_group_L_starts_im_time, my_group_L_ends_im_time, &
                                   my_group_L_sizes_im_time, &
                                   sab_orb_sub, sab_orb_all, para_env, num_diff_blk)

            ! if we do modified atomic orbitals for the primary basis (one MAO basis for D^occ one for D^virt, then
            ! we also need different 3-center overlap matrix elements!
            IF (do_mao) THEN

               natom = SIZE(particle_set)
               ALLOCATE (local_atoms_for_mao_basis(natom))
               local_atoms_for_mao_basis = 0

               CALL reserve_blocks_3c(mat_3c_overl_int_mao_for_occ, mat_munu_mao_occ_virt, qs_env, &
                                      dimen_RI, cut_RI, unit_nr, &
                                      my_group_L_starts_im_time, my_group_L_ends_im_time, &
                                      my_group_L_sizes_im_time, sab_orb_sub, sab_orb_all, para_env, num_diff_blk, &
                                      local_atoms_for_mao_basis=local_atoms_for_mao_basis)

               CALL reserve_blocks_3c(mat_3c_overl_int_mao_for_virt, mat_munu_mao_virt_occ, qs_env, &
                                      dimen_RI, cut_RI, unit_nr, &
                                      my_group_L_starts_im_time, my_group_L_ends_im_time, &
                                      my_group_L_sizes_im_time, sab_orb_sub, sab_orb_all, para_env, num_diff_blk, &
                                      local_atoms_for_mao_basis=local_atoms_for_mao_basis)

               CALL replicate_mao_coeff(mao_coeff_repl_occ, mao_coeff_occ, local_atoms_for_mao_basis, natom, para_env)

               CALL replicate_mao_coeff(mao_coeff_repl_virt, mao_coeff_virt, local_atoms_for_mao_basis, natom, para_env)

            END IF

            CALL compute_3c_overlap_int(mat_3c_overl_int, &
                                        mat_3c_overl_int_mao_for_occ, &
                                        mat_3c_overl_int_mao_for_virt, &
                                        mat_munu, mat_munu_mao_occ_virt, &
                                        qs_env, my_group_L_starts_im_time, &
                                        my_group_L_ends_im_time, my_group_L_sizes_im_time, &
                                        dimen_RI, cut_RI, sab_orb_sub, sab_orb_all, do_mao, &
                                        mao_coeff_repl_occ, mao_coeff_repl_virt, &
                                        num_diff_blk)
         ENDIF

         CALL cp_para_env_release(para_env_L)

         IF (do_mao) THEN
            CALL clean_up(mao_coeff_repl_occ, mao_coeff_repl_virt, local_atoms_for_mao_basis, natom)
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE mp2_ri_gpw_compute_in

! **************************************************************************************************
!> \brief heuristic to get a global 3d process grid that is consistent with a number of process subgroups
!>        such that balanced (ideally square) 2d grids on subgroups are obtained
!> \param para_env ...
!> \param ngroup number of process groups
!> \param map1_2d mapping of 3d grid to 2 dimensions (see dbcsr_t_pgrid_create)
!> \param map2_2d mapping of 3d grid to 2 dimensions (see dbcsr_t_pgrid_create)
!> \return process grid object compatible with DBCSR tensors
! **************************************************************************************************
   FUNCTION get_pgrid_from_ngroup(para_env, ngroup, map1_2d, map2_2d) RESULT(pgrid)
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, INTENT(IN)                                :: ngroup
      INTEGER, DIMENSION(:), INTENT(IN)                  :: map1_2d, map2_2d
      TYPE(dbcsr_t_pgrid_type)                           :: pgrid

      INTEGER                                            :: dimsplit, nproc_sub
      INTEGER, DIMENSION(2)                              :: pdims_2_2d, pdims_sub_2d
      INTEGER, DIMENSION(3)                              :: pdims_t3c

      nproc_sub = para_env%num_pe/ngroup
      pdims_sub_2d = 0
      CALL mp_dims_create(nproc_sub, pdims_sub_2d)
      pdims_2_2d = 0
      CALL mp_dims_create(pdims_sub_2d(2)*ngroup, pdims_2_2d)
      IF (SIZE(map1_2d) == 1) THEN
         dimsplit = 2
         pdims_t3c(map1_2d(1)) = pdims_sub_2d(1)
         pdims_t3c(map2_2d) = pdims_2_2d
      ELSEIF (SIZE(map2_2d) == 1) THEN
         dimsplit = 1
         pdims_t3c(map2_2d(1)) = pdims_sub_2d(1)
         pdims_t3c(map1_2d) = pdims_2_2d
      ELSE
         CPABORT("map1_2d and map2_2d not compatible with a 3d grid")
      ENDIF
      CALL dbcsr_t_pgrid_create(para_env%group, pdims_t3c, pgrid, map1_2d=map1_2d, map2_2d=map2_2d, &
                                nsplit=ngroup, dimsplit=dimsplit)
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param t_3c_M ...
!> \param pgrid_t3c_M ...
!> \param mem_cut ...
!> \param sizes_AO ...
!> \param sizes_RI ...
!> \param starts_array_mc ...
!> \param ends_array_mc ...
! **************************************************************************************************
   SUBROUTINE create_tensor_M_3c(t_3c_M, pgrid_t3c_M, mem_cut, sizes_AO, sizes_RI, &
                                 starts_array_mc, ends_array_mc)
      TYPE(dbcsr_t_type), INTENT(OUT)                    :: t_3c_M
      TYPE(dbcsr_t_pgrid_type), INTENT(IN)               :: pgrid_t3c_M
      INTEGER, INTENT(IN)                                :: mem_cut
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_AO, sizes_RI
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: starts_array_mc, ends_array_mc

      INTEGER                                            :: bsum, imem, size_AO, size_AO_cut, size_RI
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_ao_1, dist_ao_2, dist_RI, &
                                                            ends_array_mc_block, &
                                                            starts_array_mc_block
      INTEGER, DIMENSION(3)                              :: pcoord, pdims
      TYPE(dbcsr_t_distribution_type)                    :: dist

      CALL dbcsr_t_mp_environ_pgrid(pgrid_t3c_M, pdims, pcoord)

      size_RI = SIZE(sizes_RI)
      size_AO = SIZE(sizes_AO)

      ALLOCATE (dist_RI(size_RI))
      CALL cyclic_tensor_dist(size_RI, pdims(1), sizes_RI, dist_RI)

      ALLOCATE (starts_array_mc(mem_cut))
      ALLOCATE (ends_array_mc(mem_cut))
      ALLOCATE (starts_array_mc_block(mem_cut))
      ALLOCATE (ends_array_mc_block(mem_cut))

      IF (size_AO < mem_cut) THEN
         CPABORT("Use smaller MEMORY_CUT")
      ENDIF

      CALL contiguous_tensor_dist(size_AO, mem_cut, sizes_AO, limits_start=starts_array_mc_block, limits_end=ends_array_mc_block)

      bsum = 0
      DO imem = 1, mem_cut
         starts_array_mc(imem) = bsum + 1
         bsum = bsum + SUM(sizes_AO(starts_array_mc_block(imem):ends_array_mc_block(imem)))
         ends_array_mc(imem) = bsum
      ENDDO

      ALLOCATE (dist_ao_1(size_AO))
      ALLOCATE (dist_ao_2(size_AO))
      DO imem = 1, mem_cut

         size_AO_cut = ends_array_mc_block(imem) - starts_array_mc_block(imem) + 1

         IF (size_AO_cut < MIN(pdims(2), pdims(3))) THEN
            CPABORT("use smaller MEMORY_CUT or use less MPI ranks")
         ENDIF
         CALL cyclic_tensor_dist(size_AO_cut, pdims(2), sizes_AO(starts_array_mc_block(imem):ends_array_mc_block(imem)), &
                                 dist=dist_ao_1(starts_array_mc_block(imem):ends_array_mc_block(imem)))
         CALL cyclic_tensor_dist(size_AO_cut, pdims(3), sizes_AO(starts_array_mc_block(imem):ends_array_mc_block(imem)), &
                                 dist=dist_ao_2(starts_array_mc_block(imem):ends_array_mc_block(imem)))
      ENDDO

      CALL dbcsr_t_distribution_new(dist, pgrid_t3c_M, [1], [2, 3], dist_RI, dist_ao_1, dist_ao_2)
      CALL dbcsr_t_create(t_3c_M, "M (RI | AO AO)", dist, [1], [2, 3], dbcsr_type_real_8, sizes_RI, &
                          sizes_AO, sizes_AO)
      CALL dbcsr_t_distribution_destroy(dist)

   END SUBROUTINE create_tensor_M_3c

! **************************************************************************************************
!> \brief ...
!> \param t_3c_O Create tensor with overlap integrals, optionally create 3c neighbor list matching tensor
!>        distribution
!> \param mp_comm_t3c ...
!> \param ao_sizes_1 block sizes of first AO index
!> \param ao_sizes_2 block sizes of second AO index
!> \param sizes_RI block sizes of RI index
!> \param ri_metric ...
!> \param nl_3c 3c-neighborlist - if present, the sizes arguments must refer to atomic blocks
!> \param qs_env needed for nl_3c
!> \param basis_ao_1 needed for nl_3c
!> \param basis_ao_2 needed for nl_3c
!> \param basis_RI needed for nl_3c
!> \param do_kpoints ...
! **************************************************************************************************
   SUBROUTINE create_tensor_O_3c(t_3c_O, mp_comm_t3c, ao_sizes_1, ao_sizes_2, sizes_RI, &
                                 ri_metric, nl_3c, qs_env, basis_ao_1, basis_ao_2, basis_RI, &
                                 do_kpoints)
      TYPE(dbcsr_t_type), INTENT(OUT)                    :: t_3c_O
      TYPE(dbcsr_t_pgrid_type), INTENT(IN)               :: mp_comm_t3c
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_sizes_1, ao_sizes_2, sizes_RI
      TYPE(libint_potential_type), INTENT(IN), OPTIONAL  :: ri_metric
      TYPE(neighbor_list_3c_type), INTENT(OUT), OPTIONAL :: nl_3c
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: basis_ao_1, basis_ao_2, basis_RI
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_kpoints

      INTEGER                                            :: mp_comm_t3c_2, nkind, size_AO_1, &
                                                            size_AO_2, size_RI
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_ao_1, dist_ao_2, dist_RI
      INTEGER, DIMENSION(3)                              :: pcoor, pcoord, pdims
      LOGICAL                                            :: kp
      TYPE(dbcsr_t_distribution_type)                    :: dist
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      IF (PRESENT(do_kpoints)) THEN
         kp = do_kpoints
      ELSE
         kp = .FALSE.
      ENDIF

      CALL dbcsr_t_mp_environ_pgrid(mp_comm_t3c, pdims, pcoord)

      size_AO_1 = SIZE(ao_sizes_1)
      size_AO_2 = SIZE(ao_sizes_2)
      size_RI = SIZE(sizes_RI)

      ALLOCATE (dist_RI(size_RI))
      ALLOCATE (dist_ao_1(size_AO_1))
      ALLOCATE (dist_ao_2(size_AO_2))

      CALL cyclic_tensor_dist(size_RI, pdims(1), sizes_RI, dist_RI)
      CALL cyclic_tensor_dist(size_AO_1, pdims(2), ao_sizes_1, dist_ao_1)
      CALL cyclic_tensor_dist(size_AO_2, pdims(3), ao_sizes_2, dist_ao_2)

      IF (PRESENT(nl_3c)) THEN
         CPASSERT(PRESENT(qs_env))
         CPASSERT(PRESENT(ri_metric))
         CALL get_qs_env(qs_env, nkind=nkind, particle_set=particle_set)
         CPASSERT(PRESENT(basis_RI))
         CPASSERT(PRESENT(basis_ao_1))
         CPASSERT(PRESENT(basis_ao_2))
         CALL dbcsr_t_mp_environ_pgrid(mp_comm_t3c, pdims, pcoor)
         CALL mp_cart_create(mp_comm_t3c%mp_comm_2d, 3, pdims, pcoord, mp_comm_t3c_2)
         CALL distribution_3d_create(dist_3d, dist_RI, dist_ao_1, dist_ao_2, &
                                     nkind, particle_set, mp_comm_t3c_2, own_comm=.TRUE.)

         CALL build_3c_neighbor_lists(nl_3c, basis_RI, basis_ao_1, basis_ao_2, &
                                      dist_3d, ri_metric, "RPA_3c_nl", qs_env, &
                                      sym_jk=.NOT. kp, own_dist=.TRUE.)
      ENDIF

      CALL dbcsr_t_distribution_new(dist, mp_comm_t3c, [1, 2], [3], dist_RI, dist_ao_1, dist_ao_2)
      CALL dbcsr_t_create(t_3c_O, "O (RI AO | AO)", dist, [1, 2], [3], dbcsr_type_real_8, sizes_RI, &
                          ao_sizes_1, ao_sizes_2)
      CALL dbcsr_t_distribution_destroy(dist)

   END SUBROUTINE create_tensor_O_3c

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_B_L', routineP = moduleN//':'//routineN
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
      CALL dbcsr_t_distribution_new(dist_B, mp_comm_B, map_B_1, map_B_2, dist_L_P, dist_B_a, dist_B_i)
      CALL dbcsr_t_distribution_new(dist_L, mp_comm_L, map_L_1, map_L_2, dist_L_R, dist_L_P)

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_intermediate_matrices', &
         routineP = moduleN//':'//routineN

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_kpoints', &
         routineP = moduleN//':'//routineN

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'grep_Lcols', routineP = moduleN//':'//routineN

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'grep_my_integrals', &
         routineP = moduleN//':'//routineN

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

! **************************************************************************************************
!> \brief compute three center overlap integrals and write them to mat_3_overl_int
!> \param mat_3c_overl_int ...
!> \param mat_3c_overl_int_mao_for_occ ...
!> \param mat_3c_overl_int_mao_for_virt ...
!> \param mat_munu ...
!> \param mat_munu_mao_occ_virt ...
!> \param qs_env ...
!> \param my_group_L_starts_im_time ...
!> \param my_group_L_ends_im_time ...
!> \param my_group_L_sizes_im_time ...
!> \param dimen_RI ...
!> \param cut_RI ...
!> \param sab_orb_sub ...
!> \param sab_orb_all ...
!> \param do_mao ...
!> \param mao_coeff_repl_occ ...
!> \param mao_coeff_repl_virt ...
!> \param num_diff_blk ...
! **************************************************************************************************
   SUBROUTINE compute_3c_overlap_int(mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, &
                                     mat_3c_overl_int_mao_for_virt, mat_munu, mat_munu_mao_occ_virt, &
                                     qs_env, my_group_L_starts_im_time, &
                                     my_group_L_ends_im_time, my_group_L_sizes_im_time, dimen_RI, &
                                     cut_RI, sab_orb_sub, sab_orb_all, do_mao, mao_coeff_repl_occ, &
                                     mao_coeff_repl_virt, num_diff_blk)

      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int, &
                                                            mat_3c_overl_int_mao_for_occ, &
                                                            mat_3c_overl_int_mao_for_virt
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu, mat_munu_mao_occ_virt
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: my_group_L_starts_im_time, &
                                                            my_group_L_ends_im_time, &
                                                            my_group_L_sizes_im_time
      INTEGER, INTENT(IN)                                :: dimen_RI, cut_RI
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub, sab_orb_all
      LOGICAL, INTENT(IN)                                :: do_mao
      TYPE(two_dim_real_array), ALLOCATABLE, &
         DIMENSION(:), INTENT(IN)                        :: mao_coeff_repl_occ, mao_coeff_repl_virt
      INTEGER, INTENT(IN)                                :: num_diff_blk

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_3c_overlap_int', &
         routineP = moduleN//':'//routineN

      INTEGER :: acell, atom_RI, atom_RI_end, atom_RI_start, bcell, block_end_col, block_end_row, &
         block_start_col, block_start_row, col, col_basis_size_mao_occ, col_basis_size_mao_virt, &
         end_col_from_LLL_mao_occ, end_col_from_LLL_mao_virt, end_row_from_LLL_mao_occ, &
         end_row_from_LLL_mao_virt, handle, handle2, i_cut_RI, i_img, iatom, &
         iatom_basis_size_mao_occ, iatom_basis_size_mao_virt, iatom_outer, iblk, iblk_send, ikind, &
         img_col, img_row, iset, isgf_RI, isgfa, j_img, jatom, jatom_basis_size_mao_occ, &
         jatom_basis_size_mao_virt, jatom_outer, jkind, jset, jsgfb, kind_RI, LLL, LLL_set_end
      INTEGER :: LLL_set_start, maxval_1, maxval_2, maxval_3, mepos, minval_1, minval_2, minval_3, &
         my_group_L_end, my_group_L_size, my_group_L_start, natom, nco_RI, ncoa, ncob, nimg, &
         nkind, nset_RI, nseta, nsetb, nthread, offset_col_from_LLL, offset_row_from_LLL, row, &
         row_basis_size_mao_occ, row_basis_size_mao_virt, set_RI, sgf_RI, sgfa, sgfb, &
         start_col_from_LLL_mao_occ, start_col_from_LLL_mao_virt, start_row_from_LLL_mao_occ, &
         start_row_from_LLL_mao_virt
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_RI_index, kind_of
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: indices_for_uncommon_blocks
      INTEGER, DIMENSION(3)                              :: bcellvec, cell_vec, outer_cell_vec, &
                                                            size_cell_to_index
      INTEGER, DIMENSION(:), POINTER :: blk_sizes_mao_occ, blk_sizes_mao_virt, col_blk_sizes, &
         l_max_RI, l_min_RI, la_max, la_min, lb_max, lb_min, npgf_RI, npgfa, npgfb, nsgf_RI, &
         nsgfa, nsgfb, row_blk_end, row_blk_sizes, row_blk_start
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_RI, first_sgfa, first_sgfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints_cubic_RPA, found
      REAL(KIND=dp)                                      :: dab, daRI, dbRI
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: block, block_mao_for_occ, &
         block_mao_for_occ_transposed, block_mao_for_virt, block_mao_for_virt_transposed, &
         block_transposed, temp_mat_mao_occ_virt, temp_mat_mao_virt_occ
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: s_abRI, s_abRI_contr, &
                                                            s_abRI_contr_transposed
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rab_outer, raRI, rbRI
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b, set_radius_RI
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dummy_block, rpgf_RI, rpgfa, rpgfb, &
                                                            sphi_a, sphi_b, sphi_RI, zet_RI, zeta, &
                                                            zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_RI_tmp
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b, basis_set_RI
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator, nl_iterator_outer
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(two_dim_real_array), ALLOCATABLE, &
         DIMENSION(:)                                    :: blocks_to_send_around

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, atomic_kind_set, cell, particle_set, molecule_kind_set, basis_set_RI, basis_set_a, basis_set_b, &
               nl_iterator)

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      cell=cell, particle_set=particle_set, molecule_kind_set=molecule_kind_set, &
                      nkind=nkind, natom=natom, kpoints=kpoints, dft_control=dft_control)

      do_kpoints_cubic_RPA = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints

      IF (do_kpoints_cubic_RPA) THEN
         ! No MAOs kpoints
         CPASSERT(.NOT. do_mao)
         nimg = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
         size_cell_to_index(1) = SIZE(cell_to_index, 1)
         size_cell_to_index(2) = SIZE(cell_to_index, 2)
         size_cell_to_index(3) = SIZE(cell_to_index, 3)
         ! 1: row, 2: col, 3: i_img, 4: j_img, 5: cut_RI, 6: proc_to_send
         ALLOCATE (indices_for_uncommon_blocks(num_diff_blk, 6))
         indices_for_uncommon_blocks(:, :) = 0
         ALLOCATE (blocks_to_send_around(num_diff_blk))
         iblk_send = 0
         minval_1 = MINVAL(kpoints%index_to_cell(1, :))
         maxval_1 = MAXVAL(kpoints%index_to_cell(1, :))
         minval_2 = MINVAL(kpoints%index_to_cell(2, :))
         maxval_2 = MAXVAL(kpoints%index_to_cell(2, :))
         minval_3 = MINVAL(kpoints%index_to_cell(3, :))
         maxval_3 = MAXVAL(kpoints%index_to_cell(3, :))
      ELSE
         nimg = 1
      END IF

      ALLOCATE (row_blk_start(natom))
      ALLOCATE (row_blk_end(natom))

      ALLOCATE (basis_set_RI_tmp(nkind))
      CALL basis_set_list_setup(basis_set_RI_tmp, "RI_AUX", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=row_blk_start, last_sgf=row_blk_end, &
                            basis=basis_set_RI_tmp)
      DEALLOCATE (basis_set_RI_tmp)

      ALLOCATE (atom_from_RI_index(dimen_RI))

      DO LLL = 1, dimen_RI

         DO iatom = 1, natom

            IF (LLL >= row_blk_start(iatom) .AND. LLL <= row_blk_end(iatom)) THEN

               atom_from_RI_index(LLL) = iatom

            END IF

         END DO

      END DO

      CALL dbcsr_get_info(mat_munu%matrix, &
                          row_blk_size=row_blk_sizes, &
                          col_blk_size=col_blk_sizes)

      IF (do_mao) THEN

         CALL dbcsr_get_info(mat_munu_mao_occ_virt%matrix, &
                             row_blk_size=blk_sizes_mao_occ, &
                             col_blk_size=blk_sizes_mao_virt)

      END IF

      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      DO i_cut_RI = 1, cut_RI

         my_group_L_start = my_group_L_starts_im_time(i_cut_RI)
         my_group_L_end = my_group_L_ends_im_time(i_cut_RI)
         my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

         atom_RI_start = atom_from_RI_index(my_group_L_start)
         atom_RI_end = atom_from_RI_index(my_group_L_end)

         DO atom_RI = atom_RI_start, atom_RI_end

            kind_RI = kind_of(atom_RI)
            CALL get_qs_kind(qs_kind=qs_kind_set(kind_RI), basis_set=basis_set_RI, basis_type="RI_AUX")

            first_sgf_RI => basis_set_RI%first_sgf
            l_max_RI => basis_set_RI%lmax
            l_min_RI => basis_set_RI%lmin
            npgf_RI => basis_set_RI%npgf
            nset_RI = basis_set_RI%nset
            nsgf_RI => basis_set_RI%nsgf_set
            rpgf_RI => basis_set_RI%pgf_radius
            set_radius_RI => basis_set_RI%set_radius
            sphi_RI => basis_set_RI%sphi
            zet_RI => basis_set_RI%zet

            CALL neighbor_list_iterator_create(nl_iterator_outer, sab_orb_all)
            DO WHILE (neighbor_list_iterate(nl_iterator_outer) == 0)

               CALL get_iterator_info(nl_iterator_outer, &
                                      iatom=iatom_outer, jatom=jatom_outer, r=rab_outer, &
                                      cell=outer_cell_vec)

               IF (iatom_outer .NE. atom_RI .AND. jatom_outer .NE. atom_RI) CYCLE

               nthread = 1
!$             nthread = omp_get_max_threads()

               CALL neighbor_list_iterator_create(nl_iterator, sab_orb_sub, nthread=nthread)

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED (nthread,atom_RI,first_sgf_RI,l_max_RI,l_min_RI,npgf_RI,nset_RI,nsgf_RI,rpgf_RI,set_radius_RI,sphi_RI,zet_RI,&
!$OMP           iatom_outer,jatom_outer,rab_outer,kind_of,my_group_L_start,row_blk_start,&
!$OMP           mat_3c_overl_int,nl_iterator,qs_kind_set,ncoset,my_group_L_size,row_blk_sizes,&
!$OMP           col_blk_sizes, i_cut_RI, do_mao, blk_sizes_mao_occ, blk_sizes_mao_virt, outer_cell_vec, &
!$OMP           mao_coeff_repl_occ, mao_coeff_repl_virt, do_kpoints_cubic_RPA, &
!$OMP           mat_3c_overl_int_mao_for_occ, mat_3c_overl_int_mao_for_virt, kpoints, cell_to_index, nimg, &
!$OMP           minval_1, minval_2, minval_3, maxval_1, maxval_2, maxval_3) &
!$OMP   PRIVATE (mepos,iatom,jatom,rab,set_RI,raRI,rbRI,dab,daRI,dbRI,ikind,first_sgfa,la_max,la_min,npgfa,nseta,nsgfa,rpgfa,&
!$OMP            set_radius_a,sphi_a,zeta,jkind,first_sgfb,lb_max,lb_min,npgfb,nsetb,nsgfb,rpgfb,set_radius_b,sphi_b,zetb,&
!$OMP            iset,jset,LLL_set_start,nco_RI,ncoa,ncob,sgf_RI,sgfa,sgfb,isgf_RI,isgfa,jsgfb,LLL,row,col,&
!$OMP            block,block_transposed,LLL_set_end, start_row_from_LLL_mao_occ, start_col_from_LLL_mao_occ, &
!$OMP            start_row_from_LLL_mao_virt, start_col_from_LLL_mao_virt, end_row_from_LLL_mao_occ, &
!$OMP            end_col_from_LLL_mao_occ, end_row_from_LLL_mao_virt, end_col_from_LLL_mao_virt, &
!$OMP            row_basis_size_mao_occ, col_basis_size_mao_occ, row_basis_size_mao_virt, col_basis_size_mao_virt, &
!$OMP            offset_row_from_LLL,offset_col_from_LLL,block_start_row,block_end_row,block_start_col,block_end_col,&
!$OMP            handle2,s_abRI,s_abRI_contr,s_abRI_contr_transposed,basis_set_a,basis_set_b, jatom_basis_size_mao_occ, &
!$OMP            iatom_basis_size_mao_virt, jatom_basis_size_mao_virt, temp_mat_mao_occ_virt, block_mao_for_occ, &
!$OMP            block_mao_for_occ_transposed, block_mao_for_virt, block_mao_for_virt_transposed, temp_mat_mao_virt_occ, &
!$OMP            iatom_basis_size_mao_occ, found, dummy_block, cell_vec, acell, bcell, bcellvec, &
!$OMP            img_row, img_col, i_img, j_img)

               mepos = 0
!$             mepos = omp_get_thread_num()

               DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)

                  CALL get_iterator_info(nl_iterator, mepos=mepos, &
                                         iatom=iatom, jatom=jatom, r=rab, &
                                         cell=cell_vec)

                  DO i_img = 1, nimg

                     DO j_img = 1, nimg

                        IF (atom_RI .EQ. iatom_outer .AND. &
                            iatom .NE. jatom_outer .AND. &
                            jatom .NE. jatom_outer) &
                           CYCLE

                        IF (atom_RI .EQ. jatom_outer .AND. &
                            iatom .NE. iatom_outer .AND. &
                            jatom .NE. iatom_outer) &
                           CYCLE

                        IF (iatom .EQ. iatom_outer .AND. &
                            atom_RI .EQ. jatom_outer) THEN

                           IF (do_kpoints_cubic_RPA) THEN

                              IF (outer_cell_vec(1) < minval_1) CYCLE
                              IF (outer_cell_vec(1) > maxval_1) CYCLE
                              IF (outer_cell_vec(2) < minval_2) CYCLE
                              IF (outer_cell_vec(2) > maxval_2) CYCLE
                              IF (outer_cell_vec(3) < minval_3) CYCLE
                              IF (outer_cell_vec(3) > maxval_3) CYCLE

                              acell = cell_to_index(-outer_cell_vec(1), -outer_cell_vec(2), -outer_cell_vec(3))

                              IF (.NOT. (acell == i_img)) CYCLE

                              bcellvec = -outer_cell_vec + cell_vec

                              IF (bcellvec(1) < minval_1) CYCLE
                              IF (bcellvec(1) > maxval_1) CYCLE
                              IF (bcellvec(2) < minval_2) CYCLE
                              IF (bcellvec(2) > maxval_2) CYCLE
                              IF (bcellvec(3) < minval_3) CYCLE
                              IF (bcellvec(3) > maxval_3) CYCLE

                              bcell = cell_to_index(bcellvec(1), bcellvec(2), bcellvec(3))

                              IF (.NOT. (bcell == j_img)) CYCLE

                           END IF

                           raRI = rab_outer
                           rbRI = raRI - rab

                        ELSE IF (atom_RI .EQ. iatom_outer .AND. &
                                 iatom .EQ. jatom_outer) THEN

                           IF (do_kpoints_cubic_RPA) THEN
                              ! for kpoints we have the non-symmetric neighbor list
                              CYCLE
                           END IF

                           raRI = -rab_outer
                           rbRI = raRI - rab

                        ELSE IF (jatom .EQ. iatom_outer .AND. &
                                 atom_RI .EQ. jatom_outer) THEN

                           IF (do_kpoints_cubic_RPA) THEN
                              ! for kpoints we have the non-symmetric neighbor list
                              CYCLE
                           END IF

                           rbRI = rab_outer
                           raRI = rbRI + rab

                        ELSE IF (atom_RI .EQ. iatom_outer .AND. &
                                 jatom .EQ. jatom_outer) THEN

                           IF (do_kpoints_cubic_RPA) THEN
                              ! for kpoints we have the non-symmetric neighbor list
                              CYCLE
                           END IF

                           rbRI = -rab_outer
                           raRI = rbRI + rab

                        END IF

                        DO set_RI = 1, nset_RI

                           dab = SQRT(SUM(rab**2))
                           daRI = SQRT(SUM(raRI**2))
                           dbRI = SQRT(SUM(rbRI**2))

                           ikind = kind_of(iatom)
                           CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a)
                           first_sgfa => basis_set_a%first_sgf
                           la_max => basis_set_a%lmax
                           la_min => basis_set_a%lmin
                           npgfa => basis_set_a%npgf
                           nseta = basis_set_a%nset
                           nsgfa => basis_set_a%nsgf_set
                           rpgfa => basis_set_a%pgf_radius
                           set_radius_a => basis_set_a%set_radius
                           sphi_a => basis_set_a%sphi
                           zeta => basis_set_a%zet

                           jkind = kind_of(jatom)
                           CALL get_qs_kind(qs_kind=qs_kind_set(jkind), basis_set=basis_set_b)
                           first_sgfb => basis_set_b%first_sgf
                           lb_max => basis_set_b%lmax
                           lb_min => basis_set_b%lmin
                           npgfb => basis_set_b%npgf
                           nsetb = basis_set_b%nset
                           nsgfb => basis_set_b%nsgf_set
                           rpgfb => basis_set_b%pgf_radius
                           set_radius_b => basis_set_b%set_radius
                           sphi_b => basis_set_b%sphi
                           zetb => basis_set_b%zet

                           DO iset = 1, nseta

                              IF (set_radius_a(iset) + set_radius_RI(set_RI) < daRI) CYCLE

                              DO jset = 1, nsetb

                                 IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
                                 IF (set_radius_b(jset) + set_radius_RI(set_RI) < dbRI) CYCLE

                                 nco_RI = npgf_RI(set_RI)*ncoset(l_max_RI(set_RI))
                                 ncoa = npgfa(iset)*ncoset(la_max(iset))
                                 ncob = npgfb(jset)*ncoset(lb_max(jset))

                                 sgf_RI = first_sgf_RI(1, set_RI)
                                 sgfa = first_sgfa(1, iset)
                                 sgfb = first_sgfb(1, jset)

                                 LLL_set_start = 1 + sgf_RI - 1 - my_group_L_start + 1 + row_blk_start(atom_RI) - 1

                                 IF (LLL_set_start > my_group_L_size) CYCLE

                                 LLL_set_end = nsgf_RI(set_RI) + sgf_RI - 1 - my_group_L_start + 1 + row_blk_start(atom_RI) - 1

                                 IF (LLL_set_end < 1) CYCLE

                                 IF (ncoa*ncob*nco_RI > 0) THEN

                                    ALLOCATE (s_abRI(ncoa, ncob, nco_RI))
                                    s_abRI(:, :, :) = 0.0_dp

                                    CALL overlap3(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                                  lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                                  l_max_RI(set_RI), npgf_RI(set_RI), zet_RI(:, set_RI), rpgf_RI(:, set_RI), &
                                                  l_min_RI(set_RI), rab, dab, raRI, daRI, rbRI, dbRI, s_abRI)

                                    ALLOCATE (s_abRI_contr(nsgfa(iset), nsgfb(jset), nsgf_RI(set_RI)))

                                    CALL abc_contract(s_abRI_contr, s_abRI, &
                                                      sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_RI(:, sgf_RI:), &
                                                      ncoa, ncob, nco_RI, nsgfa(iset), nsgfb(jset), nsgf_RI(set_RI))

                                    ALLOCATE (s_abRI_contr_transposed(nsgfb(jset), nsgfa(iset), nsgf_RI(set_RI)))

                                    DO isgf_RI = 1, nsgf_RI(set_RI)
                                       DO isgfa = 1, nsgfa(iset)
                                          DO jsgfb = 1, nsgfb(jset)
                                             s_abRI_contr_transposed(jsgfb, isgfa, isgf_RI) = &
                                                s_abRI_contr(isgfa, jsgfb, isgf_RI)
                                          END DO
                                       END DO
                                    END DO

                                    DO isgf_RI = 1, nsgf_RI(set_RI)

                                       LLL = isgf_RI + sgf_RI - 1 - my_group_L_start + 1 + row_blk_start(atom_RI) - 1

                                       IF (LLL < 1) CYCLE
                                       IF (LLL > my_group_L_size) CYCLE

                                       IF (iatom < jatom) THEN

                                          row = iatom
                                          col = jatom
                                          img_row = i_img
                                          img_col = j_img

                                       ELSE IF (iatom >= jatom) THEN

                                          row = jatom
                                          col = iatom
                                          img_row = j_img
                                          img_col = i_img

                                       END IF

                                       ALLOCATE (block(row_blk_sizes(row), col_blk_sizes(col)*my_group_L_size))
                                       block = 0.0_dp

                                       ALLOCATE (block_transposed(col_blk_sizes(col), row_blk_sizes(row)*my_group_L_size))
                                       block_transposed = 0.0_dp

                                       offset_row_from_LLL = (LLL - 1)*row_blk_sizes(row)
                                       offset_col_from_LLL = (LLL - 1)*col_blk_sizes(col)

                                       IF (iatom < jatom) THEN

                                          block_start_row = sgfa
                                          block_end_row = sgfa + nsgfa(iset) - 1
                                          block_start_col = sgfb
                                          block_end_col = sgfb + nsgfb(jset) - 1

                                          ! factor 0.5 is necessary due to double filling due to double iterate loop
                                          block(block_start_row:block_end_row, &
                                                block_start_col + offset_col_from_LLL:block_end_col + offset_col_from_LLL) = &
                                             0.5_dp*s_abRI_contr(:, :, isgf_RI)

                                          block_transposed(block_start_col:block_end_col, &
                                                           block_start_row + offset_row_from_LLL: &
                                                           block_end_row + offset_row_from_LLL) = &
                                             0.5_dp*s_abRI_contr_transposed(:, :, isgf_RI)

                                       ELSE IF (iatom > jatom) THEN

                                          block_start_row = sgfb
                                          block_end_row = sgfb + nsgfb(jset) - 1
                                          block_start_col = sgfa
                                          block_end_col = sgfa + nsgfa(iset) - 1

                                          ! factor 0.5 is necessary due to double filling due to double iterate loop
                                          block(block_start_row:block_end_row, &
                                                block_start_col + offset_col_from_LLL:block_end_col + offset_col_from_LLL) = &
                                             0.5_dp*s_abRI_contr_transposed(:, :, isgf_RI)

                                          block_transposed(block_start_col:block_end_col, &
                                                           block_start_row + offset_row_from_LLL: &
                                                           block_end_row + offset_row_from_LLL) = &
                                             0.5_dp*s_abRI_contr(:, :, isgf_RI)

                                       ELSE IF (iatom .EQ. jatom) THEN

                                          block_start_row = sgfa
                                          block_end_row = sgfa + nsgfa(iset) - 1
                                          block_start_col = sgfb + offset_col_from_LLL
                                          block_end_col = sgfb + nsgfb(jset) - 1 + offset_col_from_LLL

                                          block(block_start_row:block_end_row, &
                                                block_start_col:block_end_col) = s_abRI_contr(:, :, isgf_RI)

                                       END IF

                                       CALL timeset(routineN//"_put_dbcsr", handle2)

!$OMP CRITICAL(putblock)
                                       NULLIFY (dummy_block)
                                       CALL dbcsr_get_block_p(matrix=mat_3c_overl_int(i_cut_RI, img_row, img_col)%matrix, &
                                                              row=row, col=col, block=dummy_block, found=found)
                                       IF (found) THEN
                                          CALL dbcsr_put_block(mat_3c_overl_int(i_cut_RI, img_row, img_col)%matrix, &
                                                               row=row, col=col, block=block, summation=.TRUE.)
                                       END IF
!$OMP END CRITICAL(putblock)

                                       IF (row .NE. col) THEN
!$OMP CRITICAL(putblock)
                                          NULLIFY (dummy_block)
                                          CALL dbcsr_get_block_p(matrix=mat_3c_overl_int(i_cut_RI, img_col, img_row)%matrix, &
                                                                 row=col, col=row, block=dummy_block, found=found)
                                          IF (found) THEN
                                             CALL dbcsr_put_block(mat_3c_overl_int(i_cut_RI, img_col, img_row)%matrix, &
                                                                  row=col, col=row, block=block_transposed, summation=.TRUE.)
                                          END IF
!$OMP END CRITICAL(putblock)
                                       END IF

                                       DEALLOCATE (block, block_transposed)

                                       IF (do_mao) THEN

                                          iatom_basis_size_mao_occ = blk_sizes_mao_occ(iatom)
                                          jatom_basis_size_mao_occ = blk_sizes_mao_occ(jatom)
                                          iatom_basis_size_mao_virt = blk_sizes_mao_virt(iatom)
                                          jatom_basis_size_mao_virt = blk_sizes_mao_virt(jatom)

                                          row_basis_size_mao_occ = blk_sizes_mao_occ(row)
                                          col_basis_size_mao_occ = blk_sizes_mao_occ(col)
                                          row_basis_size_mao_virt = blk_sizes_mao_virt(row)
                                          col_basis_size_mao_virt = blk_sizes_mao_virt(col)

                                          ALLOCATE (temp_mat_mao_occ_virt(iatom_basis_size_mao_occ, jatom_basis_size_mao_virt))
                                          ALLOCATE (temp_mat_mao_virt_occ(iatom_basis_size_mao_virt, jatom_basis_size_mao_occ))

                                          ! transform from the primary Gaussian basis into the MAO basis
                                          CALL ab_contract(temp_mat_mao_occ_virt(:, :), s_abRI_contr(:, :, isgf_RI), &
                                                           mao_coeff_repl_occ(iatom)%array(sgfa:sgfa + nsgfa(iset) - 1, :), &
                                                           mao_coeff_repl_virt(jatom)%array(sgfb:sgfb + nsgfb(jset) - 1, :), &
                                                           nsgfa(iset), nsgfb(jset), &
                                                           iatom_basis_size_mao_occ, jatom_basis_size_mao_virt)

                                          CALL ab_contract(temp_mat_mao_virt_occ(:, :), s_abRI_contr(:, :, isgf_RI), &
                                                           mao_coeff_repl_virt(iatom)%array(sgfa:sgfa + nsgfa(iset) - 1, :), &
                                                           mao_coeff_repl_occ(jatom)%array(sgfb:sgfb + nsgfb(jset) - 1, :), &
                                                           nsgfa(iset), nsgfb(jset), &
                                                           iatom_basis_size_mao_virt, jatom_basis_size_mao_occ)

                                          start_row_from_LLL_mao_occ = (LLL - 1)*row_basis_size_mao_occ + 1
                                          start_col_from_LLL_mao_occ = (LLL - 1)*col_basis_size_mao_occ + 1

                                          start_row_from_LLL_mao_virt = (LLL - 1)*row_basis_size_mao_virt + 1
                                          start_col_from_LLL_mao_virt = (LLL - 1)*col_basis_size_mao_virt + 1

                                          end_row_from_LLL_mao_occ = LLL*row_basis_size_mao_occ
                                          end_col_from_LLL_mao_occ = LLL*col_basis_size_mao_occ

                                          end_row_from_LLL_mao_virt = LLL*row_basis_size_mao_virt
                                          end_col_from_LLL_mao_virt = LLL*col_basis_size_mao_virt

                                          ALLOCATE (block_mao_for_occ(row_basis_size_mao_occ, &
                                                                      col_basis_size_mao_virt*my_group_L_size))
                                          block_mao_for_occ = 0.0_dp

                                          ALLOCATE (block_mao_for_occ_transposed(col_basis_size_mao_occ, &
                                                                                 row_basis_size_mao_virt*my_group_L_size))
                                          block_mao_for_occ_transposed = 0.0_dp

                                          ALLOCATE (block_mao_for_virt(row_basis_size_mao_virt, &
                                                                       col_basis_size_mao_occ*my_group_L_size))
                                          block_mao_for_virt = 0.0_dp

                                          ALLOCATE (block_mao_for_virt_transposed(col_basis_size_mao_virt, &
                                                                                  row_basis_size_mao_occ*my_group_L_size))
                                          block_mao_for_virt_transposed = 0.0_dp

                                          IF (iatom < jatom) THEN

                                             block_mao_for_occ(:, start_col_from_LLL_mao_virt:end_col_from_LLL_mao_virt) = &
                                                0.5_dp*temp_mat_mao_occ_virt(:, :)

                                             block_mao_for_occ_transposed(:, start_row_from_LLL_mao_virt: &
                                                                          end_row_from_LLL_mao_virt) = &
                                                0.5_dp*TRANSPOSE(temp_mat_mao_virt_occ(:, :))

                                             block_mao_for_virt(:, start_col_from_LLL_mao_occ:end_col_from_LLL_mao_occ) = &
                                                0.5_dp*temp_mat_mao_virt_occ(:, :)

                                             block_mao_for_virt_transposed(:, start_row_from_LLL_mao_occ: &
                                                                           end_row_from_LLL_mao_occ) = &
                                                0.5_dp*TRANSPOSE(temp_mat_mao_occ_virt(:, :))

                                          ELSE IF (iatom > jatom) THEN

                                             block_mao_for_occ(:, start_col_from_LLL_mao_virt:end_col_from_LLL_mao_virt) = &
                                                0.5_dp*TRANSPOSE(temp_mat_mao_virt_occ(:, :))

                                             block_mao_for_occ_transposed(:, start_row_from_LLL_mao_virt: &
                                                                          end_row_from_LLL_mao_virt) = &
                                                0.5_dp*temp_mat_mao_occ_virt(:, :)

                                             block_mao_for_virt(:, start_col_from_LLL_mao_occ:end_col_from_LLL_mao_occ) = &
                                                0.5_dp*TRANSPOSE(temp_mat_mao_occ_virt(:, :))

                                             block_mao_for_virt_transposed(:, start_row_from_LLL_mao_occ: &
                                                                           end_row_from_LLL_mao_occ) = &
                                                0.5_dp*temp_mat_mao_virt_occ(:, :)

                                          ELSE IF (iatom .EQ. jatom) THEN

                                             block_mao_for_occ(:, start_col_from_LLL_mao_virt:end_col_from_LLL_mao_virt) = &
                                                temp_mat_mao_occ_virt(:, :)

                                             block_mao_for_virt(:, start_col_from_LLL_mao_occ:end_col_from_LLL_mao_occ) = &
                                                temp_mat_mao_virt_occ(:, :)

                                          END IF

!$OMP CRITICAL(putblock)
                                          CALL dbcsr_put_block(mat_3c_overl_int_mao_for_occ(i_cut_RI, i_img, j_img)%matrix, &
                                                               row=row, col=col, block=block_mao_for_occ, summation=.TRUE.)
!$OMP END CRITICAL(putblock)

!$OMP CRITICAL(putblock)
                                          CALL dbcsr_put_block(mat_3c_overl_int_mao_for_virt(i_cut_RI, i_img, j_img)%matrix, &
                                                               row=row, col=col, block=block_mao_for_virt, summation=.TRUE.)
!$OMP END CRITICAL(putblock)

                                          IF (row .NE. col) THEN
!$OMP CRITICAL(putblock)
                                             CALL dbcsr_put_block(mat_3c_overl_int_mao_for_occ(i_cut_RI, j_img, i_img)%matrix, &
                                                                  row=col, col=row, block=block_mao_for_occ_transposed, &
                                                                  summation=.TRUE.)
!$OMP END CRITICAL(putblock)

!$OMP CRITICAL(putblock)
                                             CALL dbcsr_put_block(mat_3c_overl_int_mao_for_virt(i_cut_RI, j_img, i_img)%matrix, &
                                                                  row=col, col=row, block=block_mao_for_virt_transposed, &
                                                                  summation=.TRUE.)
!$OMP END CRITICAL(putblock)

                                          END IF

                                          DEALLOCATE (temp_mat_mao_occ_virt, temp_mat_mao_virt_occ, block_mao_for_occ, &
                                                      block_mao_for_occ_transposed, block_mao_for_virt, &
                                                      block_mao_for_virt_transposed)

                                       END IF

                                       CALL timestop(handle2)

                                    END DO ! single RI basis functions

                                    DEALLOCATE (s_abRI, s_abRI_contr, s_abRI_contr_transposed)

                                 END IF ! number of triples > 0

                              END DO ! jset

                           END DO ! iset

                        END DO ! set RI

                     END DO ! j_img

                  END DO ! i_img

               END DO ! inner neighbor list

!$OMP   END PARALLEL

               CALL neighbor_list_iterator_release(nl_iterator)

            END DO ! outer neighbor list

            CALL neighbor_list_iterator_release(nl_iterator_outer)

         END DO ! atom_RI

      END DO ! RI index ranges

      ! when doing kpoints, there might be some zero blocks
      IF (do_kpoints_cubic_RPA) THEN

         NULLIFY (mat_3c_overl_int_tmp)
         CALL dbcsr_allocate_matrix_set(mat_3c_overl_int_tmp, cut_RI, nimg, nimg)

         DO i_img = 1, nimg
            DO j_img = 1, nimg
               DO i_cut_RI = 1, cut_RI

                  CALL dbcsr_filter(mat_3c_overl_int(i_cut_RI, i_img, j_img)%matrix, 1.0E-30_dp)

                  ALLOCATE (mat_3c_overl_int_tmp(i_cut_RI, i_img, j_img)%matrix)
                  CALL dbcsr_create(mat_3c_overl_int_tmp(i_cut_RI, i_img, j_img)%matrix, &
                                    template=mat_3c_overl_int(i_cut_RI, i_img, j_img)%matrix)
                  CALL dbcsr_copy(mat_3c_overl_int_tmp(i_cut_RI, i_img, j_img)%matrix, &
                                  mat_3c_overl_int(i_cut_RI, i_img, j_img)%matrix)
                  CALL dbcsr_set(mat_3c_overl_int(i_cut_RI, i_img, j_img)%matrix, 0.0_dp)
                  CALL dbcsr_filter(mat_3c_overl_int(i_cut_RI, i_img, j_img)%matrix, 1.0E-30_dp)
                  CALL dbcsr_complete_redistribute(mat_3c_overl_int_tmp(i_cut_RI, i_img, j_img)%matrix, &
                                                   mat_3c_overl_int(i_cut_RI, i_img, j_img)%matrix, summation=.TRUE.)

               END DO
            END DO
         END DO

         DO iblk = 1, num_diff_blk
            DEALLOCATE (blocks_to_send_around(iblk)%array)
         END DO

         DEALLOCATE (indices_for_uncommon_blocks, blocks_to_send_around)
         CALL dbcsr_deallocate_matrix_set(mat_3c_overl_int_tmp)

      END IF

      DEALLOCATE (atom_from_RI_index)

      DEALLOCATE (row_blk_start, row_blk_end)

      CALL timestop(handle)

   END SUBROUTINE compute_3c_overlap_int

! **************************************************************************************************
!> \brief compute three center overlap integrals and write them to mat_3_overl_int
!> \param mat_3c_overl_int ...
!> \param mat_munu ...
!> \param qs_env ...
!> \param dimen_RI ...
!> \param cut_RI ...
!> \param unit_nr ...
!> \param my_group_L_starts_im_time ...
!> \param my_group_L_ends_im_time ...
!> \param my_group_L_sizes_im_time ...
!> \param sab_orb_sub ...
!> \param sab_orb_all ...
!> \param para_env ...
!> \param num_diff_blk ...
!> \param local_atoms_for_mao_basis ...
! **************************************************************************************************
   SUBROUTINE reserve_blocks_3c(mat_3c_overl_int, mat_munu, qs_env, dimen_RI, &
                                cut_RI, unit_nr, my_group_L_starts_im_time, my_group_L_ends_im_time, &
                                my_group_L_sizes_im_time, sab_orb_sub, sab_orb_all, para_env, &
                                num_diff_blk, local_atoms_for_mao_basis)

      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: dimen_RI, cut_RI, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: my_group_L_starts_im_time, &
                                                            my_group_L_ends_im_time, &
                                                            my_group_L_sizes_im_time
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub, sab_orb_all
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(OUT)                               :: num_diff_blk
      INTEGER, ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: local_atoms_for_mao_basis

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reserve_blocks_3c', &
         routineP = moduleN//':'//routineN

      INTEGER :: acell, atom_RI, atom_RI_end, atom_RI_start, bcell, blk, blk_cnt, col, handle, &
         handle2, i_cut_RI, i_img, i_img_outer, i_loop, iatom, iatom_outer, ikind, iset, isgf_RI, &
         j_img, j_img_outer, jatom, jatom_outer, jkind, jset, kind_RI, LLL, LLL_set_end, &
         LLL_set_start, maxval_1, maxval_2, maxval_3, minval_1, minval_2, minval_3, &
         my_group_L_end, my_group_L_size, my_group_L_start, natom, nblkrows_total, nco_RI, ncoa, &
         ncob, nimg, nkind, nset_RI, nseta, nsetb, num_loops, row, set_RI, sgf_RI, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_RI_index, kind_of, tmp
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: old_col, old_row
      INTEGER, DIMENSION(3)                              :: cell_vec, outer_cell_vec
      INTEGER, DIMENSION(:), POINTER :: col_blk_sizes, col_blk_sizes_scaled, l_max_RI, l_min_RI, &
         la_max, la_min, lb_max, lb_min, npgf_RI, npgfa, npgfb, nsgf_RI, nsgfa, nsgfb, &
         row_blk_end, row_blk_offset, row_blk_sizes, row_blk_sizes_scaled, row_blk_start
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_RI, first_sgfa, first_sgfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints_cubic_RPA, new_block
      REAL(KIND=dp)                                      :: dab, daRI, dbRI
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rab_outer, raRI, rbRI
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b, set_radius_RI
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_RI, rpgfa, rpgfb, sphi_a, sphi_b, &
                                                            sphi_RI, zet_RI, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type)                      :: dist_sub
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_RI_tmp
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b, basis_set_RI
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator, nl_iterator_outer
      TYPE(one_dim_int_array), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: cols_to_alloc, rows_to_alloc
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, atomic_kind_set, cell, particle_set, molecule_kind_set, basis_set_RI, basis_set_a, basis_set_b, &
               nl_iterator)

      do_kpoints_cubic_RPA = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, dft_control=dft_control, &
                      cell=cell, particle_set=particle_set, molecule_kind_set=molecule_kind_set, &
                      nkind=nkind, natom=natom, kpoints=kpoints)

      IF (do_kpoints_cubic_RPA) THEN
         ! set up new kpoint type with periodic images according to eps_grid from MP2 section
         ! instead of eps_pgf_orb from QS section
         CALL kpoint_init_cell_index(kpoints, sab_orb_sub, para_env, dft_control)

         nimg = dft_control%nimages
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)

         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
            "3C_OVERLAP_INTEGRALS_INFO| Number of periodic images considered:", nimg

         num_diff_blk = 0
      ELSE
         nimg = 1
      END IF

      ! allocate matrix for storing 3c overlap integrals
      NULLIFY (mat_3c_overl_int)
      CALL dbcsr_allocate_matrix_set(mat_3c_overl_int, cut_RI, nimg, nimg)

      ALLOCATE (row_blk_start(natom))
      ALLOCATE (row_blk_end(natom))

      ALLOCATE (basis_set_RI_tmp(nkind))
      CALL basis_set_list_setup(basis_set_RI_tmp, "RI_AUX", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=row_blk_start, last_sgf=row_blk_end, &
                            basis=basis_set_RI_tmp)
      DEALLOCATE (basis_set_RI_tmp)

      ALLOCATE (atom_from_RI_index(dimen_RI))

      DO LLL = 1, dimen_RI

         DO iatom = 1, natom

            IF (LLL >= row_blk_start(iatom) .AND. LLL <= row_blk_end(iatom)) THEN

               atom_from_RI_index(LLL) = iatom

            END IF

         END DO

      END DO

      CALL dbcsr_get_info(mat_munu%matrix, &
                          row_blk_offset=row_blk_offset, &
                          row_blk_size=row_blk_sizes, &
                          nblkrows_total=nblkrows_total)

      IF (do_kpoints_cubic_RPA) THEN

         minval_1 = MINVAL(kpoints%index_to_cell(1, :))
         maxval_1 = MAXVAL(kpoints%index_to_cell(1, :))
         minval_2 = MINVAL(kpoints%index_to_cell(2, :))
         maxval_2 = MAXVAL(kpoints%index_to_cell(2, :))
         minval_3 = MINVAL(kpoints%index_to_cell(3, :))
         maxval_3 = MAXVAL(kpoints%index_to_cell(3, :))

      END IF

      ALLOCATE (kind_of(natom))
      ALLOCATE (rows_to_alloc(nimg, nimg))
      ALLOCATE (cols_to_alloc(nimg, nimg))

      ALLOCATE (old_row(nimg, nimg))
      old_row = -1
      ALLOCATE (old_col(nimg, nimg))
      old_col = -1

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      DO i_cut_RI = 1, cut_RI

         DO i_img = 1, nimg
         DO j_img = 1, nimg
            ALLOCATE (rows_to_alloc(i_img, j_img)%array(1))
            rows_to_alloc(i_img, j_img)%array(1) = 0
            ALLOCATE (cols_to_alloc(i_img, j_img)%array(1))
            cols_to_alloc(i_img, j_img)%array(1) = 0
         END DO
         END DO

         my_group_L_start = my_group_L_starts_im_time(i_cut_RI)
         my_group_L_end = my_group_L_ends_im_time(i_cut_RI)
         my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

         atom_RI_start = atom_from_RI_index(my_group_L_start)
         atom_RI_end = atom_from_RI_index(my_group_L_end)

         DO atom_RI = atom_RI_start, atom_RI_end

            kind_RI = kind_of(atom_RI)
            CALL get_qs_kind(qs_kind=qs_kind_set(kind_RI), basis_set=basis_set_RI, basis_type="RI_AUX")

            first_sgf_RI => basis_set_RI%first_sgf

            l_max_RI => basis_set_RI%lmax
            l_min_RI => basis_set_RI%lmin
            npgf_RI => basis_set_RI%npgf
            nset_RI = basis_set_RI%nset
            nsgf_RI => basis_set_RI%nsgf_set
            rpgf_RI => basis_set_RI%pgf_radius
            set_radius_RI => basis_set_RI%set_radius
            sphi_RI => basis_set_RI%sphi
            zet_RI => basis_set_RI%zet

            CALL neighbor_list_iterator_create(nl_iterator_outer, sab_orb_all)
            DO WHILE (neighbor_list_iterate(nl_iterator_outer) == 0)

               CALL get_iterator_info(nl_iterator_outer, &
                                      iatom=iatom_outer, jatom=jatom_outer, r=rab_outer, cell=outer_cell_vec)

               IF (iatom_outer .NE. atom_RI .AND. jatom_outer .NE. atom_RI) CYCLE

               CALL neighbor_list_iterator_create(nl_iterator, sab_orb_sub)
               DO WHILE (neighbor_list_iterate(nl_iterator) == 0)

                  CALL get_iterator_info(nl_iterator, &
                                         iatom=iatom, jatom=jatom, r=rab, cell=cell_vec)

                  DO i_img_outer = 1, nimg

                     DO j_img_outer = 1, nimg

                        IF (i_img_outer > j_img_outer) CYCLE
                        IF (i_img_outer < j_img_outer) THEN
                           num_loops = 2
                        ELSE IF (i_img_outer == j_img_outer) THEN
                           num_loops = 1
                        END IF

                        DO i_loop = 1, num_loops

                           IF (i_loop == 1) THEN
                              i_img = i_img_outer
                              j_img = j_img_outer
                           ELSE IF (i_loop == 2) THEN
                              i_img = j_img_outer
                              j_img = i_img_outer
                           END IF

                           IF (atom_RI .EQ. iatom_outer .AND. &
                               iatom .NE. jatom_outer .AND. &
                               jatom .NE. jatom_outer) &
                              CYCLE

                           IF (atom_RI .EQ. jatom_outer .AND. &
                               iatom .NE. iatom_outer .AND. &
                               jatom .NE. iatom_outer) &
                              CYCLE

                           IF (iatom .EQ. iatom_outer .AND. &
                               atom_RI .EQ. jatom_outer) THEN

                              IF (do_kpoints_cubic_RPA) THEN

                                 IF (outer_cell_vec(1) < minval_1) CYCLE
                                 IF (outer_cell_vec(1) > maxval_1) CYCLE
                                 IF (outer_cell_vec(2) < minval_2) CYCLE
                                 IF (outer_cell_vec(2) > maxval_2) CYCLE
                                 IF (outer_cell_vec(3) < minval_3) CYCLE
                                 IF (outer_cell_vec(3) > maxval_3) CYCLE

                                 IF (-outer_cell_vec(1) + cell_vec(1) < minval_1) CYCLE
                                 IF (-outer_cell_vec(1) + cell_vec(1) > maxval_1) CYCLE
                                 IF (-outer_cell_vec(2) + cell_vec(2) < minval_2) CYCLE
                                 IF (-outer_cell_vec(2) + cell_vec(2) > maxval_2) CYCLE
                                 IF (-outer_cell_vec(3) + cell_vec(3) < minval_3) CYCLE
                                 IF (-outer_cell_vec(3) + cell_vec(3) > maxval_3) CYCLE

                                 acell = cell_to_index(-outer_cell_vec(1), -outer_cell_vec(2), -outer_cell_vec(3))
                                 IF (.NOT. (acell == i_img)) CYCLE
                                 bcell = cell_to_index(-outer_cell_vec(1) + cell_vec(1), &
                                                       -outer_cell_vec(2) + cell_vec(2), &
                                                       -outer_cell_vec(3) + cell_vec(3))
                                 IF (.NOT. (bcell == j_img)) CYCLE

                              END IF

                              raRI = rab_outer
                              rbRI = raRI - rab

                           ELSE IF (atom_RI .EQ. iatom_outer .AND. &
                                    iatom .EQ. jatom_outer) THEN

                              IF (do_kpoints_cubic_RPA) THEN
                                 ! for kpoints we have the non-symmetric neighbor list
                                 CYCLE
                              END IF

                              raRI = -rab_outer
                              rbRI = raRI - rab

                           ELSE IF (jatom .EQ. iatom_outer .AND. &
                                    atom_RI .EQ. jatom_outer) THEN

                              IF (do_kpoints_cubic_RPA) THEN
                                 ! for kpoints we have the non-symmetric neighbor list
                                 CYCLE
                              END IF

                              rbRI = rab_outer
                              raRI = rbRI + rab

                           ELSE IF (atom_RI .EQ. iatom_outer .AND. &
                                    jatom .EQ. jatom_outer) THEN

                              IF (do_kpoints_cubic_RPA) THEN
                                 ! for kpoints we have the non-symmetric neighbor list
                                 CYCLE
                              END IF

                              rbRI = -rab_outer
                              raRI = rbRI + rab

                           END IF

                           IF (PRESENT(local_atoms_for_mao_basis)) THEN
                              local_atoms_for_mao_basis(iatom) = iatom
                              local_atoms_for_mao_basis(jatom) = jatom
                           END IF

                           dab = SQRT(SUM(rab**2))
                           daRI = SQRT(SUM(raRI**2))
                           dbRI = SQRT(SUM(rbRI**2))

                           ikind = kind_of(iatom)
                           CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a)
                           first_sgfa => basis_set_a%first_sgf
                           la_max => basis_set_a%lmax
                           la_min => basis_set_a%lmin
                           npgfa => basis_set_a%npgf
                           nseta = basis_set_a%nset
                           nsgfa => basis_set_a%nsgf_set
                           rpgfa => basis_set_a%pgf_radius
                           set_radius_a => basis_set_a%set_radius
                           sphi_a => basis_set_a%sphi
                           zeta => basis_set_a%zet

                           jkind = kind_of(jatom)
                           CALL get_qs_kind(qs_kind=qs_kind_set(jkind), basis_set=basis_set_b)
                           first_sgfb => basis_set_b%first_sgf
                           lb_max => basis_set_b%lmax
                           lb_min => basis_set_b%lmin
                           npgfb => basis_set_b%npgf
                           nsetb = basis_set_b%nset
                           nsgfb => basis_set_b%nsgf_set
                           rpgfb => basis_set_b%pgf_radius
                           set_radius_b => basis_set_b%set_radius
                           sphi_b => basis_set_b%sphi
                           zetb => basis_set_b%zet

                           DO set_RI = 1, nset_RI

                              DO iset = 1, nseta

                                 IF (set_radius_a(iset) + set_radius_RI(set_RI) < daRI) CYCLE

                                 DO jset = 1, nsetb

                                    IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
                                    IF (set_radius_b(jset) + set_radius_RI(set_RI) < dbRI) CYCLE

                                    nco_RI = npgf_RI(set_RI)*ncoset(l_max_RI(set_RI))
                                    ncoa = npgfa(iset)*ncoset(la_max(iset))
                                    ncob = npgfb(jset)*ncoset(lb_max(jset))

                                    sgf_RI = first_sgf_RI(1, set_RI)
                                    sgfa = first_sgfa(1, iset)
                                    sgfb = first_sgfb(1, jset)

                                    LLL_set_start = 1 + sgf_RI - 1 - my_group_L_start + 1 + row_blk_start(atom_RI) - 1

                                    IF (LLL_set_start > my_group_L_size) CYCLE

                                    LLL_set_end = nsgf_RI(set_RI) + sgf_RI - 1 - my_group_L_start + 1 + row_blk_start(atom_RI) - 1

                                    IF (LLL_set_end < 1) CYCLE

                                    IF (ncoa*ncob*nco_RI > 0) THEN

                                       DO isgf_RI = 1, nsgf_RI(set_RI)

                                          LLL = isgf_RI + sgf_RI - 1 - my_group_L_start + 1 + row_blk_start(atom_RI) - 1

                                          IF (LLL < 1) CYCLE
                                          IF (LLL > my_group_L_size) CYCLE

                                          IF (iatom < jatom) THEN

                                             row = iatom
                                             col = jatom

                                          ELSE IF (iatom >= jatom) THEN

                                             row = jatom
                                             col = iatom

                                          END IF

                                          IF (SIZE(rows_to_alloc(i_img_outer, j_img_outer)%array) == 1 .AND. &
                                              rows_to_alloc(i_img_outer, j_img_outer)%array(1) == 0) THEN

                                             IF (row == col) THEN

                                                rows_to_alloc(i_img_outer, j_img_outer)%array = row
                                                cols_to_alloc(i_img_outer, j_img_outer)%array = col

                                             ELSE

                                                DEALLOCATE (rows_to_alloc(i_img_outer, j_img_outer)%array)
                                                DEALLOCATE (cols_to_alloc(i_img_outer, j_img_outer)%array)

                                                ALLOCATE (rows_to_alloc(i_img_outer, j_img_outer)%array(2))
                                                ALLOCATE (cols_to_alloc(i_img_outer, j_img_outer)%array(2))

                                                rows_to_alloc(i_img_outer, j_img_outer)%array(1) = row
                                                rows_to_alloc(i_img_outer, j_img_outer)%array(2) = col

                                                cols_to_alloc(i_img_outer, j_img_outer)%array(1) = col
                                                cols_to_alloc(i_img_outer, j_img_outer)%array(2) = row

                                             END IF

                                             old_row(i_img_outer, j_img_outer) = row
                                             old_col(i_img_outer, j_img_outer) = col

                                          ELSE

                                             new_block = .FALSE.

                                             IF ((old_row(i_img_outer, j_img_outer) .NE. row) .OR. &
                                                 (old_col(i_img_outer, j_img_outer) .NE. col)) THEN

                                                blk_cnt = SIZE(rows_to_alloc(i_img_outer, j_img_outer)%array)

                                                new_block = .TRUE.

                                                DO blk = 1, blk_cnt

                                                   IF (rows_to_alloc(i_img_outer, j_img_outer)%array(blk) == row .AND. &
                                                       cols_to_alloc(i_img_outer, j_img_outer)%array(blk) == col) THEN

                                                      new_block = .FALSE.
                                                      EXIT

                                                   END IF

                                                END DO

                                                IF (new_block .AND. row == col) THEN

                                                   ALLOCATE (tmp(blk_cnt))
                                                   tmp(1:blk_cnt) = rows_to_alloc(i_img_outer, j_img_outer)%array(1:blk_cnt)
                                                   DEALLOCATE (rows_to_alloc(i_img_outer, j_img_outer)%array)
                                                   ALLOCATE (rows_to_alloc(i_img_outer, j_img_outer)%array(blk_cnt + 1))
                                                   rows_to_alloc(i_img_outer, j_img_outer)%array(1:blk_cnt) = tmp(1:blk_cnt)
                                                   rows_to_alloc(i_img_outer, j_img_outer)%array(blk_cnt + 1) = row
                                                   tmp(1:blk_cnt) = cols_to_alloc(i_img_outer, j_img_outer)%array(1:blk_cnt)
                                                   DEALLOCATE (cols_to_alloc(i_img_outer, j_img_outer)%array)
                                                   ALLOCATE (cols_to_alloc(i_img_outer, j_img_outer)%array(blk_cnt + 1))
                                                   cols_to_alloc(i_img_outer, j_img_outer)%array(1:blk_cnt) = tmp(1:blk_cnt)
                                                   cols_to_alloc(i_img_outer, j_img_outer)%array(blk_cnt + 1) = col
                                                   DEALLOCATE (tmp)

                                                   old_row(i_img_outer, j_img_outer) = row
                                                   old_col(i_img_outer, j_img_outer) = col

                                                ELSE IF (new_block .AND. row .NE. col) THEN

                                                   ALLOCATE (tmp(blk_cnt))
                                                   tmp(1:blk_cnt) = rows_to_alloc(i_img_outer, j_img_outer)%array(1:blk_cnt)
                                                   DEALLOCATE (rows_to_alloc(i_img_outer, j_img_outer)%array)
                                                   ALLOCATE (rows_to_alloc(i_img_outer, j_img_outer)%array(blk_cnt + 2))
                                                   rows_to_alloc(i_img_outer, j_img_outer)%array(1:blk_cnt) = tmp(1:blk_cnt)
                                                   rows_to_alloc(i_img_outer, j_img_outer)%array(blk_cnt + 1) = row
                                                   rows_to_alloc(i_img_outer, j_img_outer)%array(blk_cnt + 2) = col
                                                   tmp(1:blk_cnt) = cols_to_alloc(i_img_outer, j_img_outer)%array(1:blk_cnt)
                                                   DEALLOCATE (cols_to_alloc(i_img_outer, j_img_outer)%array)
                                                   ALLOCATE (cols_to_alloc(i_img_outer, j_img_outer)%array(blk_cnt + 2))
                                                   cols_to_alloc(i_img_outer, j_img_outer)%array(1:blk_cnt) = tmp(1:blk_cnt)
                                                   cols_to_alloc(i_img_outer, j_img_outer)%array(blk_cnt + 1) = col
                                                   cols_to_alloc(i_img_outer, j_img_outer)%array(blk_cnt + 2) = row
                                                   DEALLOCATE (tmp)

                                                   old_row(i_img_outer, j_img_outer) = row
                                                   old_col(i_img_outer, j_img_outer) = col

                                                END IF

                                             END IF

                                          END IF

                                       END DO ! single RI basis functions

                                    END IF ! number of triples > 0

                                 END DO ! jset

                              END DO ! iset

                           END DO ! i_loop for i_img and j_img for kpoints

                        END DO ! j_img

                     END DO ! i_img

                  END DO ! set RI

               END DO ! inner neighbor list

               CALL neighbor_list_iterator_release(nl_iterator)

            END DO ! outer neighbor list

            CALL neighbor_list_iterator_release(nl_iterator_outer)

         END DO ! atom_RI

         DO i_img = 1, nimg

            DO j_img = 1, nimg

               IF (i_img > j_img) CYCLE

               CALL timeset(routineN//"_reserve_dbcsr", handle2)

               blk_cnt = SIZE(rows_to_alloc(i_img, j_img)%array)

               CALL dbcsr_get_info(mat_munu%matrix, &
                                   row_blk_size=row_blk_sizes, &
                                   col_blk_size=col_blk_sizes, &
                                   distribution=dist_sub)

               ALLOCATE (col_blk_sizes_scaled(SIZE(col_blk_sizes)))
               col_blk_sizes_scaled(:) = col_blk_sizes(:)*my_group_L_size

               ALLOCATE (mat_3c_overl_int(i_cut_RI, i_img, j_img)%matrix)

               CALL dbcsr_create(matrix=mat_3c_overl_int(i_cut_RI, i_img, j_img)%matrix, &
                                 name="(munuP)", &
                                 dist=dist_sub, &
                                 matrix_type=dbcsr_type_no_symmetry, &
                                 row_blk_size=row_blk_sizes, &
                                 col_blk_size=col_blk_sizes_scaled, &
                                 nze=0)

               IF (rows_to_alloc(i_img, j_img)%array(1) .NE. 0) THEN

                  ! that should be okay for the kpoints sinceif block (row, col) is allocated, also (col, row)
                  ! afterwards, the matrix will be filtered anyway
                  CALL dbcsr_reserve_blocks(mat_3c_overl_int(i_cut_RI, i_img, j_img)%matrix, &
                                            rows=rows_to_alloc(i_img, j_img)%array(1:blk_cnt), &
                                            cols=cols_to_alloc(i_img, j_img)%array(1:blk_cnt))

               END IF

               CALL dbcsr_finalize(mat_3c_overl_int(i_cut_RI, i_img, j_img)%matrix)

               IF (i_img .NE. j_img) THEN

                  ALLOCATE (row_blk_sizes_scaled(SIZE(row_blk_sizes)))
                  row_blk_sizes_scaled(:) = row_blk_sizes(:)*my_group_L_size

                  ALLOCATE (mat_3c_overl_int(i_cut_RI, j_img, i_img)%matrix)

                  CALL dbcsr_create(matrix=mat_3c_overl_int(i_cut_RI, j_img, i_img)%matrix, &
                                    name="(munuP)", &
                                    dist=dist_sub, &
                                    matrix_type=dbcsr_type_no_symmetry, &
                                    row_blk_size=col_blk_sizes, &
                                    col_blk_size=row_blk_sizes_scaled, &
                                    nze=0)

                  DEALLOCATE (row_blk_sizes_scaled)

                  IF (rows_to_alloc(i_img, j_img)%array(1) .NE. 0) THEN
                     CALL dbcsr_reserve_blocks(mat_3c_overl_int(i_cut_RI, j_img, i_img)%matrix, &
                                               rows=cols_to_alloc(i_img, j_img)%array(1:blk_cnt), &
                                               cols=rows_to_alloc(i_img, j_img)%array(1:blk_cnt))
                  END IF
                  CALL dbcsr_finalize(mat_3c_overl_int(i_cut_RI, j_img, i_img)%matrix)
               END IF

               DEALLOCATE (col_blk_sizes_scaled)

               CALL timestop(handle2)

            END DO ! j_img

         END DO ! i_img

         DO i_img = 1, nimg
         DO j_img = 1, nimg
            DEALLOCATE (rows_to_alloc(i_img, j_img)%array)
            DEALLOCATE (cols_to_alloc(i_img, j_img)%array)
         END DO
         END DO

      END DO ! cut_RI

      DEALLOCATE (rows_to_alloc, cols_to_alloc)
      DEALLOCATE (atom_from_RI_index)
      DEALLOCATE (row_blk_start, row_blk_end)

      CALL timestop(handle)

   END SUBROUTINE reserve_blocks_3c

! **************************************************************************************************
!> \brief ...
!> \param mao_coeff_repl ...
!> \param mao_coeff ...
!> \param local_atoms_for_mao_basis ...
!> \param natom ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE replicate_mao_coeff(mao_coeff_repl, mao_coeff, local_atoms_for_mao_basis, natom, para_env)
      TYPE(two_dim_real_array), ALLOCATABLE, &
         DIMENSION(:)                                    :: mao_coeff_repl
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coeff
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: local_atoms_for_mao_basis
      INTEGER, INTENT(IN)                                :: natom
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_mao_coeff', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_rec, block_size, col, col_size, end_msg, handle, iatom, imepos, mao_size, &
         num_atom_from_imepos, prim_size, row, row_size, start_msg
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_counter, entry_counter, &
                                                            num_entries_atoms_rec, &
                                                            num_entries_atoms_send, sizes_rec, &
                                                            sizes_send
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes_mao, row_blk_sizes_prim
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mao_coeff(1)%matrix, &
                          row_blk_size=row_blk_sizes_prim, &
                          col_blk_size=col_blk_sizes_mao)

      ALLOCATE (mao_coeff_repl(1:natom))

      ALLOCATE (num_entries_atoms_send(0:2*para_env%num_pe - 1))
      num_entries_atoms_send(:) = 0

      ALLOCATE (num_entries_atoms_rec(0:2*para_env%num_pe - 1))
      num_entries_atoms_rec(:) = 0

      DO iatom = 1, natom

         ! check whether local atom is there
         IF (local_atoms_for_mao_basis(iatom) == iatom) THEN
            ! row index belongs to primary basis, col index belongs to mao basis
            prim_size = row_blk_sizes_prim(iatom)
            mao_size = col_blk_sizes_mao(iatom)
            ALLOCATE (mao_coeff_repl(iatom)%array(prim_size, mao_size))

            CALL dbcsr_get_stored_coordinates(mao_coeff(1)%matrix, iatom, iatom, imepos)

            num_entries_atoms_rec(2*imepos) = num_entries_atoms_rec(2*imepos) + prim_size*mao_size
            num_entries_atoms_rec(2*imepos + 1) = num_entries_atoms_rec(2*imepos + 1) + 1

         END IF

      END DO

      IF (para_env%num_pe > 1) THEN
         CALL mp_alltoall(num_entries_atoms_rec, num_entries_atoms_send, 2, para_env%group)
      ELSE
         num_entries_atoms_send(0:1) = num_entries_atoms_rec(0:1)
      END IF

      ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env%num_pe - 1))

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe - 1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_atoms_rec(2*imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_atoms_send(2*imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_entries_atoms_rec(2*imepos + 1), 6))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_entries_atoms_send(2*imepos + 1), 6))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = 0

      ALLOCATE (atom_counter(0:para_env%num_pe - 1))
      atom_counter = 0

      DO iatom = 1, natom

         ! check whether local atom is there
         IF (local_atoms_for_mao_basis(iatom) == iatom) THEN

            CALL dbcsr_get_stored_coordinates(mao_coeff(1)%matrix, iatom, iatom, imepos)

            atom_counter(imepos) = atom_counter(imepos) + 1

            buffer_rec(imepos)%indx(atom_counter(imepos), 1) = iatom

         END IF

      END DO

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      ALLOCATE (sizes_rec(0:para_env%num_pe - 1))
      ALLOCATE (sizes_send(0:para_env%num_pe - 1))

      DO imepos = 0, para_env%num_pe - 1

         sizes_send(imepos) = num_entries_atoms_send(2*imepos)
         sizes_rec(imepos) = num_entries_atoms_rec(2*imepos)

      END DO

      ! change rec and send!
      CALL communicate_buffer(para_env, sizes_send, sizes_rec, buffer_send, buffer_rec, req_array, do_msg=.FALSE.)

      atom_counter(:) = 0

      CALL dbcsr_iterator_start(iter, mao_coeff(1)%matrix)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         ! for the mao trafo matrix, only the diagonal blocks should be allocated
         CPASSERT(row == col)

         DO imepos = 0, para_env%num_pe - 1

            IF (ANY(buffer_send(imepos)%indx(:, 1) == row)) THEN

               block_size = row_size*col_size
               start_msg = entry_counter(imepos) + 1
               end_msg = entry_counter(imepos) + block_size

               buffer_send(imepos)%msg(start_msg:end_msg) = RESHAPE(data_block(1:row_size, 1:col_size), (/block_size/))

               entry_counter(imepos) = entry_counter(imepos) + block_size
               atom_counter(imepos) = atom_counter(imepos) + 1

               buffer_send(imepos)%indx(atom_counter(imepos), 2) = row
               buffer_send(imepos)%indx(atom_counter(imepos), 3) = row_size
               buffer_send(imepos)%indx(atom_counter(imepos), 4) = col_size
               buffer_send(imepos)%indx(atom_counter(imepos), 5) = start_msg
               buffer_send(imepos)%indx(atom_counter(imepos), 6) = end_msg

            END IF

         END DO

      END DO ! blocks

      CALL dbcsr_iterator_stop(iter)

      ! communicate the mao coefficients
      CALL communicate_buffer(para_env, sizes_rec, sizes_send, buffer_rec, buffer_send, req_array)

      DEALLOCATE (req_array, sizes_rec, sizes_send)

      ! fill mao_coeff_repl
      DO imepos = 0, para_env%num_pe - 1

         num_atom_from_imepos = num_entries_atoms_rec(2*imepos + 1)

         DO atom_rec = 1, num_atom_from_imepos

            iatom = buffer_rec(imepos)%indx(atom_rec, 2)
            row_size = buffer_rec(imepos)%indx(atom_rec, 3)
            col_size = buffer_rec(imepos)%indx(atom_rec, 4)
            start_msg = buffer_rec(imepos)%indx(atom_rec, 5)
            end_msg = buffer_rec(imepos)%indx(atom_rec, 6)

            mao_coeff_repl(iatom)%array(1:row_size, 1:col_size) = &
               RESHAPE(buffer_rec(imepos)%msg(start_msg:end_msg), (/row_size, col_size/))

         END DO

      END DO

      DEALLOCATE (num_entries_atoms_send, num_entries_atoms_rec, entry_counter, atom_counter)

      DO imepos = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)

      CALL timestop(handle)

   END SUBROUTINE replicate_mao_coeff

!***************************************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param mao_coeff_repl_occ ...
!> \param mao_coeff_repl_virt ...
!> \param local_atoms_for_mao_basis ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE clean_up(mao_coeff_repl_occ, mao_coeff_repl_virt, local_atoms_for_mao_basis, natom)
      TYPE(two_dim_real_array), ALLOCATABLE, &
         DIMENSION(:), INTENT(INOUT)                     :: mao_coeff_repl_occ, mao_coeff_repl_virt
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: local_atoms_for_mao_basis
      INTEGER, INTENT(IN)                                :: natom

      INTEGER                                            :: iatom

      DO iatom = 1, natom
         ! check whether local atom is there
         IF (local_atoms_for_mao_basis(iatom) == iatom) THEN
            DEALLOCATE (mao_coeff_repl_occ(iatom)%array)
         END IF
      END DO

      DO iatom = 1, natom
         ! check whether local atom is there
         IF (local_atoms_for_mao_basis(iatom) == iatom) THEN
            DEALLOCATE (mao_coeff_repl_virt(iatom)%array)
         END IF
      END DO

      DEALLOCATE (mao_coeff_repl_occ, mao_coeff_repl_virt)

      DEALLOCATE (local_atoms_for_mao_basis)

   END SUBROUTINE clean_up

! **************************************************************************************************
!> \brief ...
!> \param my_group_L_starts_im_time ...
!> \param my_group_L_ends_im_time ...
!> \param my_group_L_sizes_im_time ...
!> \param dimen_RI ...
!> \param ngroup ...
!> \param cut_memory ...
!> \param cut_RI ...
!> \param color_sub ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE setup_group_L_im_time(my_group_L_starts_im_time, my_group_L_ends_im_time, my_group_L_sizes_im_time, &
                                    dimen_RI, ngroup, cut_memory, cut_RI, color_sub, qs_env)

      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: my_group_L_starts_im_time, &
                                                            my_group_L_ends_im_time, &
                                                            my_group_L_sizes_im_time
      INTEGER, INTENT(IN)                                :: dimen_RI, ngroup, cut_memory
      INTEGER, INTENT(OUT)                               :: cut_RI
      INTEGER, INTENT(IN)                                :: color_sub
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_group_L_im_time', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: cut_RI_old, handle, iblk_RI, icut, &
                                                            itmp(2), natom, nblks_RI, nkind, &
                                                            size_RI, start_RI
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_group_L_ends_im_time_tmp, &
                                                            my_group_L_starts_im_time_tmp, &
                                                            row_blk_end_RI, row_blk_offset_RI
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes_RI
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(group_dist_d1_type)                           :: gd_array_mem_cut
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_ri_aux
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL create_group_dist(gd_array_mem_cut, cut_memory, dimen_RI)

      ALLOCATE (my_group_L_starts_im_time(1:cut_memory))
      my_group_L_starts_im_time = 0
      ALLOCATE (my_group_L_ends_im_time(1:cut_memory))
      my_group_L_ends_im_time = 0

      DO icut = 0, cut_memory - 1

         CALL get_group_dist(gd_array_mem_cut, icut, sizes=size_RI, starts=start_RI)

         itmp = get_limit(size_RI, ngroup, color_sub)

         my_group_L_starts_im_time(icut + 1) = itmp(1) + start_RI - 1
         my_group_L_ends_im_time(icut + 1) = itmp(2) + start_RI - 1

      END DO

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set)

      natom = SIZE(particle_set)
      nkind = SIZE(atomic_kind_set)

      ALLOCATE (row_blk_sizes_RI(natom))

      ALLOCATE (basis_set_ri_aux(nkind))
      CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes_RI, basis=basis_set_ri_aux)
      DEALLOCATE (basis_set_ri_aux)

      ! check whether my_group_L_im_time covers different atoms. In this case, cut my_group_L_im_time
      ! especially improves memory a bit
      nblks_RI = SIZE(row_blk_sizes_RI)

      ALLOCATE (row_blk_offset_RI(nblks_RI))
      ALLOCATE (row_blk_end_RI(nblks_RI))

      row_blk_offset_RI(1) = 1
      DO iblk_RI = 2, nblks_RI
         row_blk_offset_RI(iblk_RI) = row_blk_offset_RI(iblk_RI - 1) + row_blk_sizes_RI(iblk_RI - 1)
      END DO

      row_blk_end_RI(nblks_RI) = dimen_RI
      DO iblk_RI = 1, nblks_RI - 1
         row_blk_end_RI(iblk_RI) = row_blk_offset_RI(iblk_RI + 1) - 1
      END DO

      cut_RI = cut_memory
      cut_RI_old = 0

      DO WHILE (cut_RI .NE. cut_RI_old)

         cut_RI_old = cut_RI

         DO iblk_RI = 1, nblks_RI

            DO icut = 1, cut_RI_old

               IF (row_blk_offset_RI(iblk_RI) >= my_group_L_starts_im_time(icut) .AND. &
                   row_blk_offset_RI(iblk_RI) <= my_group_L_ends_im_time(icut) .AND. &
                   row_blk_end_RI(iblk_RI) < my_group_L_ends_im_time(icut)) THEN

                  cut_RI = cut_RI + 1

                  ALLOCATE (my_group_L_starts_im_time_tmp(cut_RI - 1))
                  ALLOCATE (my_group_L_ends_im_time_tmp(cut_RI - 1))

                  my_group_L_starts_im_time_tmp(:) = my_group_L_starts_im_time(:)
                  my_group_L_ends_im_time_tmp(:) = my_group_L_ends_im_time(:)

                  DEALLOCATE (my_group_L_starts_im_time)
                  DEALLOCATE (my_group_L_ends_im_time)

                  ALLOCATE (my_group_L_starts_im_time(cut_RI))
                  ALLOCATE (my_group_L_ends_im_time(cut_RI))

                  my_group_L_starts_im_time(1:icut) = my_group_L_starts_im_time_tmp(1:icut)
                  my_group_L_starts_im_time(icut + 1) = row_blk_offset_RI(iblk_RI + 1)
                  IF (cut_RI >= icut + 2) THEN
                     my_group_L_starts_im_time(icut + 2:cut_RI) = my_group_L_starts_im_time_tmp(icut + 1:cut_RI - 1)
                  END IF

                  IF (icut - 1 >= 1) THEN
                     my_group_L_ends_im_time(1:icut - 1) = my_group_L_ends_im_time_tmp(1:icut - 1)
                  END IF
                  my_group_L_ends_im_time(icut) = row_blk_offset_RI(iblk_RI + 1) - 1
                  my_group_L_ends_im_time(icut + 1:cut_RI) = my_group_L_ends_im_time_tmp(icut:cut_RI - 1)

                  DEALLOCATE (my_group_L_starts_im_time_tmp)
                  DEALLOCATE (my_group_L_ends_im_time_tmp)

               END IF

            END DO

         END DO

      END DO

      ALLOCATE (my_group_L_sizes_im_time(cut_RI))
      my_group_L_sizes_im_time(:) = my_group_L_ends_im_time(:) - my_group_L_starts_im_time(:) + 1

      ALLOCATE (qs_env%mp2_env%ri_rpa_im_time_util(cut_memory))

      qs_env%mp2_env%ri_rpa_im_time_util(1)%cut_RI = cut_RI

      ALLOCATE (qs_env%mp2_env%ri_rpa_im_time_util(1)%my_group_L_starts_im_time(1:cut_RI))
      qs_env%mp2_env%ri_rpa_im_time_util(1)%my_group_L_starts_im_time(:) = my_group_L_starts_im_time(:)
      ALLOCATE (qs_env%mp2_env%ri_rpa_im_time_util(1)%my_group_L_ends_im_time(1:cut_RI))
      qs_env%mp2_env%ri_rpa_im_time_util(1)%my_group_L_ends_im_time(:) = my_group_L_ends_im_time(:)
      ALLOCATE (qs_env%mp2_env%ri_rpa_im_time_util(1)%my_group_L_sizes_im_time(1:cut_RI))
      qs_env%mp2_env%ri_rpa_im_time_util(1)%my_group_L_sizes_im_time(:) = my_group_L_sizes_im_time(:)

      CALL release_group_dist(gd_array_mem_cut)
      DEALLOCATE (row_blk_sizes_RI)

      CALL timestop(handle)

   END SUBROUTINE setup_group_L_im_time

END MODULE mp2_integrals
