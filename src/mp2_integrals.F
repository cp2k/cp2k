!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
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
   USE OMP_LIB,                         ONLY: omp_get_num_threads,&
                                              omp_get_thread_num
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
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
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_get_info, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_release_p, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE dbt_api,                         ONLY: &
        dbt_clear, dbt_contract, dbt_copy, dbt_create, dbt_destroy, dbt_distribution_destroy, &
        dbt_distribution_new, dbt_distribution_type, dbt_filter, dbt_get_block, dbt_get_info, &
        dbt_get_stored_coordinates, dbt_mp_environ_pgrid, dbt_pgrid_create, dbt_pgrid_destroy, &
        dbt_pgrid_type, dbt_put_block, dbt_reserve_blocks, dbt_split_blocks, dbt_type
   USE group_dist_types,                ONLY: create_group_dist,&
                                              get_group_dist,&
                                              group_dist_d1_type
   USE hfx_types,                       ONLY: alloc_containers,&
                                              block_ind_type,&
                                              hfx_compression_type
   USE input_constants,                 ONLY: &
        do_eri_gpw, do_eri_mme, do_eri_os, do_potential_coulomb, do_potential_id, &
        do_potential_long, do_potential_short, do_potential_truncated, kp_weights_W_auto, &
        kp_weights_W_tailored, kp_weights_W_uniform
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: kpoint_type
   USE libint_2c_3c,                    ONLY: compare_potential_types,&
                                              libint_potential_type
   USE machine,                         ONLY: m_flush
   USE message_passing,                 ONLY: mp_cart_type,&
                                              mp_comm_type,&
                                              mp_para_env_type
   USE mp2_eri,                         ONLY: mp2_eri_3c_integrate
   USE mp2_eri_gpw,                     ONLY: cleanup_gpw,&
                                              mp2_eri_3c_integrate_gpw,&
                                              prepare_gpw
   USE mp2_ri_2c,                       ONLY: get_2c_integrals
   USE mp2_types,                       ONLY: three_dim_real_array
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
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
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_integrals'

   PUBLIC :: mp2_ri_gpw_compute_in, compute_kpoints

   TYPE intermediate_matrix_type
      TYPE(dbcsr_type) :: matrix_ia_jnu, matrix_ia_jb
      INTEGER :: max_row_col_local = 0
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: local_col_row_info
      TYPE(cp_fm_type) :: fm_BIb_jb = cp_fm_type()
      CHARACTER(LEN=default_string_length) :: descr = ""
   END TYPE intermediate_matrix_type

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
!> \param fm_matrix_PQ ...
!> \param fm_matrix_L_kpoints ...
!> \param fm_matrix_Minv_L_kpoints ...
!> \param fm_matrix_Minv ...
!> \param fm_matrix_Minv_Vtrunc_Minv ...
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
!> \author Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE mp2_ri_gpw_compute_in(BIb_C, BIb_C_gw, BIb_C_bse_ij, BIb_C_bse_ab, gd_array, gd_B_virtual, &
                                    dimen_RI, dimen_RI_red, qs_env, para_env, para_env_sub, color_sub, &
                                    cell, particle_set, atomic_kind_set, qs_kind_set, mo_coeff, &
                                    fm_matrix_PQ, fm_matrix_L_kpoints, fm_matrix_Minv_L_kpoints, &
                                    fm_matrix_Minv, fm_matrix_Minv_Vtrunc_Minv, &
                                    nmo, homo, mat_munu, &
                                    sab_orb_sub, mo_coeff_o, mo_coeff_v, mo_coeff_all, &
                                    mo_coeff_gw, eps_filter, unit_nr, &
                                    mp2_memory, calc_PQ_cond_num, calc_forces, blacs_env_sub, my_do_gw, do_bse, &
                                    gd_B_all, starts_array_mc, ends_array_mc, &
                                    starts_array_mc_block, ends_array_mc_block, &
                                    gw_corr_lev_occ, gw_corr_lev_virt, &
                                    do_im_time, do_kpoints_cubic_RPA, kpoints, &
                                    t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                                    ri_metric, gd_B_occ_bse, gd_B_virt_bse)

      TYPE(three_dim_real_array), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: BIb_C, BIb_C_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: BIb_C_bse_ij, BIb_C_bse_ab
      TYPE(group_dist_d1_type), INTENT(OUT)              :: gd_array
      TYPE(group_dist_d1_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: gd_B_virtual
      INTEGER, INTENT(OUT)                               :: dimen_RI, dimen_RI_red
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, INTENT(IN)                                :: color_sub
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_matrix_PQ
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_matrix_L_kpoints, &
                                                            fm_matrix_Minv_L_kpoints, &
                                                            fm_matrix_Minv, &
                                                            fm_matrix_Minv_Vtrunc_Minv
      INTEGER, INTENT(IN)                                :: nmo
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_munu
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: sab_orb_sub
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: mo_coeff_o, mo_coeff_v, mo_coeff_all, &
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
      TYPE(dbt_type), INTENT(OUT)                        :: t_3c_M
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: t_3c_O
      TYPE(hfx_compression_type), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: t_3c_O_compressed
      TYPE(block_ind_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: t_3c_O_ind
      TYPE(libint_potential_type), INTENT(IN)            :: ri_metric
      TYPE(group_dist_d1_type), INTENT(OUT)              :: gd_B_occ_bse, gd_B_virt_bse

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_gpw_compute_in'

      INTEGER :: cm, cut_memory, cut_memory_int, eri_method, gw_corr_lev_total, handle, handle2, &
         handle4, i, i_counter, i_mem, ibasis, ispin, itmp(2), j, jcell, kcell, LLL, min_bsize, &
         my_B_all_end, my_B_all_size, my_B_all_start, my_B_occ_bse_end, my_B_occ_bse_size, &
         my_B_occ_bse_start, my_B_virt_bse_end, my_B_virt_bse_size, my_B_virt_bse_start, &
         my_group_L_end, my_group_L_size, my_group_L_start, n_rep, natom, ngroup, nimg, nkind, &
         nspins, potential_type, ri_metric_type
      INTEGER(int_8)                                     :: nze
      INTEGER, ALLOCATABLE, DIMENSION(:) :: dist_AO_1, dist_AO_2, dist_RI, &
         ends_array_mc_block_int, ends_array_mc_int, my_B_size, my_B_virtual_end, &
         my_B_virtual_start, sizes_AO, sizes_AO_split, sizes_RI, sizes_RI_split, &
         starts_array_mc_block_int, starts_array_mc_int, virtual
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(3)                              :: bounds_3c, pcoord, pdims, pdims_t3c, &
                                                            periodic
      LOGICAL                                            :: do_gpw, do_kpoints_from_Gamma, do_svd, &
                                                            memory_info
      REAL(KIND=dp) :: compression_factor, cutoff_old, eps_pgf_orb, eps_pgf_orb_old, mem_for_iaK, &
         memory_3c, occ, omega_pot, rc_ang, relative_cutoff_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: e_cutoff_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: my_Lrows, my_Vrows
      TYPE(cp_eri_mme_param), POINTER                    :: eri_param
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_munu_local_L
      TYPE(dbt_pgrid_type)                               :: pgrid_t3c_M, pgrid_t3c_overl
      TYPE(dbt_type)                                     :: t_3c_overl_int_template, t_3c_tmp
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_3c_overl_int
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_ao, basis_set_ri_aux
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis, ri_basis
      TYPE(intermediate_matrix_type)                     :: intermed_mat_bse_ab, intermed_mat_bse_ij
      TYPE(intermediate_matrix_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: intermed_mat, intermed_mat_gw
      TYPE(mp_cart_type)                                 :: mp_comm_t3c_2
      TYPE(neighbor_list_3c_type)                        :: nl_3c
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: pot_g, psi_L, rho_g, rho_r
      TYPE(section_vals_type), POINTER                   :: qs_section
      TYPE(task_list_type), POINTER                      :: task_list_sub

      CALL timeset(routineN, handle)

      CALL cite_reference(DelBen2013)

      nspins = SIZE(homo)

      ALLOCATE (virtual(nspins))
      virtual(:) = nmo - homo(:)
      gw_corr_lev_total = gw_corr_lev_virt + gw_corr_lev_occ

      eri_method = qs_env%mp2_env%eri_method
      eri_param => qs_env%mp2_env%eri_mme_param
      do_svd = qs_env%mp2_env%do_svd
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

      do_kpoints_from_Gamma = qs_env%mp2_env%ri_rpa_im_time%do_kpoints_from_Gamma
      IF (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         kpoints=kpoints)
      END IF
      IF (do_kpoints_from_Gamma) THEN
         CALL compute_kpoints(qs_env, kpoints, unit_nr)
      END IF

      IF (do_bse) THEN
         CPASSERT(my_do_gw)
         CPASSERT(.NOT. do_im_time)
         ! GPW integrals have to be implemented later
         CPASSERT(.NOT. (eri_method == do_eri_gpw))
      END IF

      ngroup = para_env%num_pe/para_env_sub%num_pe

      ! Preparations for MME method to compute ERIs
      IF (qs_env%mp2_env%eri_method .EQ. do_eri_mme) THEN
         ! cell might have changed, so we need to reset parameters
         CALL cp_eri_mme_set_params(eri_param, cell, qs_kind_set, basis_type_1="ORB", basis_type_2="RI_AUX", para_env=para_env)
      END IF

      CALL get_cell(cell=cell, periodic=periodic)
      ! for minimax Ewald summation, full periodicity is required
      IF (eri_method == do_eri_mme) THEN
         CPASSERT(periodic(1) == 1 .AND. periodic(2) == 1 .AND. periodic(3) == 1)
      END IF

      IF (do_svd .AND. (do_kpoints_from_Gamma .OR. do_kpoints_cubic_RPA)) THEN
         CPABORT("SVD with kpoints not implemented yet!")
      END IF

      CALL get_2c_integrals(qs_env, eri_method, eri_param, para_env, para_env_sub, mp2_memory, &
                            my_Lrows, my_Vrows, fm_matrix_PQ, ngroup, color_sub, dimen_RI, dimen_RI_red, &
                            kpoints, my_group_L_size, my_group_L_start, my_group_L_end, &
                            gd_array, calc_PQ_cond_num .AND. .NOT. do_svd, do_svd, &
                            qs_env%mp2_env%potential_parameter, ri_metric, &
                            fm_matrix_L_kpoints, fm_matrix_Minv_L_kpoints, fm_matrix_Minv, fm_matrix_Minv_Vtrunc_Minv, &
                            do_im_time, do_kpoints_from_Gamma .OR. do_kpoints_cubic_RPA, qs_env%mp2_env%mp2_gpw%eps_pgf_orb_S, &
                            qs_kind_set, sab_orb_sub, calc_forces, unit_nr)

      IF (unit_nr > 0) THEN
         ASSOCIATE (ri_metric => qs_env%mp2_env%ri_metric)
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
               WRITE (unit_nr, FMT="(T3,A,T74,A)") &
                  "RI_INFO| RI metric: ", "OVERLAP"
            CASE (do_potential_truncated)
               WRITE (unit_nr, FMT="(T3,A,T64,A)") &
                  "RI_INFO| RI metric: ", "TRUNCATED COULOMB"
               rc_ang = cp_unit_from_cp2k(ri_metric%cutoff_radius, "angstrom")
               WRITE (unit_nr, '(T3,A,T61,F20.2)') &
                  "RI_INFO| Cutoff Radius [angstrom]:     ", rc_ang
            END SELECT
         END ASSOCIATE
      END IF

      IF (calc_forces .AND. .NOT. do_im_time) THEN
         ! we need (P|Q)^(-1/2) for future use, just save it
         ! in a fully (home made) distributed way
         itmp = get_limit(dimen_RI, para_env_sub%num_pe, para_env_sub%mepos)
         lll = itmp(2) - itmp(1) + 1
         ALLOCATE (qs_env%mp2_env%ri_grad%PQ_half(lll, my_group_L_size))
         qs_env%mp2_env%ri_grad%PQ_half(:, :) = my_Lrows(itmp(1):itmp(2), 1:my_group_L_size)
         IF (.NOT. compare_potential_types(qs_env%mp2_env%ri_metric, qs_env%mp2_env%potential_parameter)) THEN
            ALLOCATE (qs_env%mp2_env%ri_grad%operator_half(lll, my_group_L_size))
            qs_env%mp2_env%ri_grad%operator_half(:, :) = my_Vrows(itmp(1):itmp(2), 1:my_group_L_size)
            DEALLOCATE (my_Vrows)
         END IF
      END IF

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
            "RI_INFO| Number of auxiliary basis functions:", dimen_RI, &
            "GENERAL_INFO| Number of basis functions:", nmo, &
            "GENERAL_INFO| Number of occupied orbitals:", homo(1), &
            "GENERAL_INFO| Number of virtual orbitals:", virtual(1)
         IF (do_svd) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "RI_INFO| Reduced auxiliary basis set size:", dimen_RI_red
         END IF

         mem_for_iaK = dimen_RI*REAL(SUM(homo*virtual), KIND=dp)*8.0_dp/(1024_dp**2)

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
      END IF

      CALL para_env%sync() ! sync to see memory output

      ! in case we do imaginary time, we need the overlap tensor (alpha beta P) or trunc. Coulomb tensor
      IF (.NOT. do_im_time) THEN

         ALLOCATE (gd_B_virtual(nspins), intermed_mat(nspins))
         ALLOCATE (my_B_virtual_start(nspins), my_B_virtual_end(nspins), my_B_size(nspins))
         DO ispin = 1, nspins

            CALL create_intermediate_matrices(intermed_mat(ispin), mo_coeff_o(ispin)%matrix, virtual(ispin), homo(ispin), &
                                              TRIM(ADJUSTL(cp_to_string(ispin))), blacs_env_sub, para_env_sub)

            CALL create_group_dist(gd_B_virtual(ispin), para_env_sub%num_pe, virtual(ispin))
            CALL get_group_dist(gd_B_virtual(ispin), para_env_sub%mepos, my_B_virtual_start(ispin), my_B_virtual_end(ispin), &
                                my_B_size(ispin))

         END DO

         ! in the case of G0W0, we need (K|nm), n,m may be occ or virt (m restricted to corrected levels)
         IF (my_do_gw) THEN

            ALLOCATE (intermed_mat_gw(nspins))
            DO ispin = 1, nspins
               CALL create_intermediate_matrices(intermed_mat_gw(ispin), mo_coeff_gw(ispin)%matrix, &
                                                 nmo, gw_corr_lev_total, &
                                                 "gw_"//TRIM(ADJUSTL(cp_to_string(ispin))), &
                                                 blacs_env_sub, para_env_sub)

            END DO

            CALL create_group_dist(gd_B_all, para_env_sub%num_pe, nmo)
            CALL get_group_dist(gd_B_all, para_env_sub%mepos, my_B_all_start, my_B_all_end, my_B_all_size)

            IF (do_bse) THEN

               ! virt x virt matrices
               CALL create_intermediate_matrices(intermed_mat_bse_ab, mo_coeff_v(1)%matrix, virtual(1), virtual(1), &
                                                 "bse_ab", blacs_env_sub, para_env_sub)

               CALL create_group_dist(gd_B_virt_bse, para_env_sub%num_pe, virtual(1))
               CALL get_group_dist(gd_B_virt_bse, para_env_sub%mepos, my_B_virt_bse_start, my_B_virt_bse_end, my_B_virt_bse_size)

               ! occ x occ matrices
               CALL create_intermediate_matrices(intermed_mat_bse_ij, mo_coeff_o(1)%matrix, homo(1), homo(1), &
                                                 "bse_ij", blacs_env_sub, para_env_sub)

               CALL create_group_dist(gd_B_occ_bse, para_env_sub%num_pe, homo(1))
               CALL get_group_dist(gd_B_occ_bse, para_env_sub%mepos, my_B_occ_bse_start, my_B_occ_bse_end, my_B_occ_bse_size)

            END IF
         END IF

         ! array that will store the (ia|K) integrals
         ALLOCATE (BIb_C(nspins))
         DO ispin = 1, nspins
            ALLOCATE (BIb_C(ispin)%array(my_group_L_size, my_B_size(ispin), homo(ispin)))
            BIb_C(ispin)%array = 0.0_dp
         END DO

         ! in the case of GW, we also need (nm|K)
         IF (my_do_gw) THEN

            ALLOCATE (BIb_C_gw(nspins))
            DO ispin = 1, nspins
               ALLOCATE (BIb_C_gw(ispin)%array(my_group_L_size, my_B_all_size, gw_corr_lev_total))
               BIb_C_gw(ispin)%array = 0.0_dp
            END DO

         END IF

         IF (do_bse) THEN

            ALLOCATE (BIb_C_bse_ij(my_group_L_size, my_B_occ_bse_size, homo(1)))
            BIb_C_bse_ij = 0.0_dp

            ALLOCATE (BIb_C_bse_ab(my_group_L_size, my_B_virt_bse_size, virtual(1)))
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
            END DO
            CALL mp2_eri_3c_integrate(eri_param, ri_metric, para_env_sub, qs_env, &
                                      first_c=my_group_L_start, last_c=my_group_L_end, &
                                      mat_ab=mat_munu_local_L, &
                                      basis_type_a="ORB", basis_type_b="ORB", &
                                      basis_type_c="RI_AUX", &
                                      sab_nl=sab_orb_sub, eri_method=eri_method)

            DO ispin = 1, nspins
               DO LLL = 1, my_group_L_size
                  CALL ao_to_mo_and_store_B(para_env_sub, mat_munu_local_L(LLL), intermed_mat(ispin), &
                                            BIb_C(ispin)%array(LLL, :, :), &
                                            mo_coeff_o(ispin)%matrix, mo_coeff_v(ispin)%matrix, &
                                            eps_filter, &
                                            my_B_virtual_end(ispin), my_B_virtual_start(ispin))
               END DO
               CALL contract_B_L(BIb_C(ispin)%array, my_Lrows, gd_B_virtual(ispin)%sizes, &
                                 gd_array%sizes, qs_env%mp2_env%eri_blksize, &
                                 ngroup, color_sub, para_env, para_env_sub)
            END DO

            IF (my_do_gw) THEN

               DO ispin = 1, nspins
                  DO LLL = 1, my_group_L_size
                     CALL ao_to_mo_and_store_B(para_env_sub, mat_munu_local_L(LLL), intermed_mat_gw(ispin), &
                                               BIb_C_gw(ispin)%array(LLL, :, :), &
                                               mo_coeff_gw(ispin)%matrix, mo_coeff_all(ispin)%matrix, eps_filter, &
                                               my_B_all_end, my_B_all_start)
                  END DO
                  CALL contract_B_L(BIb_C_gw(ispin)%array, my_Lrows, gd_B_all%sizes, gd_array%sizes, qs_env%mp2_env%eri_blksize, &
                                    ngroup, color_sub, para_env, para_env_sub)
               END DO
            END IF

            IF (do_bse) THEN

               ! B^ab_P matrix elements for BSE
               DO LLL = 1, my_group_L_size
                  CALL ao_to_mo_and_store_B(para_env_sub, mat_munu_local_L(LLL), intermed_mat_bse_ab, &
                                            BIb_C_bse_ab(LLL, :, :), &
                                            mo_coeff_v(1)%matrix, mo_coeff_v(1)%matrix, eps_filter, &
                                            my_B_all_end, my_B_all_start)
               END DO
               CALL contract_B_L(BIb_C_bse_ab, my_Lrows, gd_B_virt_bse%sizes, gd_array%sizes, qs_env%mp2_env%eri_blksize, &
                                 ngroup, color_sub, para_env, para_env_sub)

               ! B^ij_P matrix elements for BSE
               DO LLL = 1, my_group_L_size
                  CALL ao_to_mo_and_store_B(para_env_sub, mat_munu_local_L(LLL), intermed_mat_bse_ij, &
                                            BIb_C_bse_ij(LLL, :, :), &
                                            mo_coeff_o(1)%matrix, mo_coeff_o(1)%matrix, eps_filter, &
                                            my_B_occ_bse_end, my_B_occ_bse_start)
               END DO
               CALL contract_B_L(BIb_C_bse_ij, my_Lrows, gd_B_occ_bse%sizes, gd_array%sizes, qs_env%mp2_env%eri_blksize, &
                                 ngroup, color_sub, para_env, para_env_sub)

            END IF

            DO LLL = 1, my_group_L_size
               CALL dbcsr_release_p(mat_munu_local_L(LLL)%matrix)
            END DO
            DEALLOCATE (mat_munu_local_L)

         ELSE IF (do_gpw) THEN

            CALL prepare_gpw(qs_env, dft_control, e_cutoff_old, cutoff_old, relative_cutoff_old, para_env_sub, pw_env_sub, &
                             auxbas_pw_pool, poisson_env, task_list_sub, rho_r, rho_g, pot_g, psi_L, sab_orb_sub)

            DO i_counter = 1, my_group_L_size

               CALL mp2_eri_3c_integrate_gpw(mo_coeff(1), psi_L, rho_g, atomic_kind_set, qs_kind_set, cell, dft_control, &
                                             particle_set, pw_env_sub, my_Lrows(:, i_counter), poisson_env, rho_r, pot_g, &
                                             ri_metric, mat_munu, qs_env, task_list_sub)

               DO ispin = 1, nspins
                  CALL ao_to_mo_and_store_B(para_env_sub, mat_munu, intermed_mat(ispin), &
                                            BIb_C(ispin)%array(i_counter, :, :), &
                                            mo_coeff_o(ispin)%matrix, mo_coeff_v(ispin)%matrix, eps_filter, &
                                            my_B_virtual_end(ispin), my_B_virtual_start(ispin))

               END DO

               IF (my_do_gw) THEN
                  ! transform (K|mu nu) to (K|nm), n corresponds to corrected GW levels, m is in nmo
                  DO ispin = 1, nspins
                     CALL ao_to_mo_and_store_B(para_env_sub, mat_munu, intermed_mat_gw(ispin), &
                                               BIb_C_gw(ispin)%array(i_counter, :, :), &
                                               mo_coeff_gw(ispin)%matrix, mo_coeff_all(ispin)%matrix, eps_filter, &
                                               my_B_all_end, my_B_all_start)

                  END DO
               END IF

            END DO

            CALL cleanup_gpw(qs_env, e_cutoff_old, cutoff_old, relative_cutoff_old, para_env_sub, pw_env_sub, &
                             task_list_sub, auxbas_pw_pool, rho_r, rho_g, pot_g, psi_L)
         ELSE
            CPABORT("Integration method not implemented!")
         END IF

         CALL timestop(handle2)

         DEALLOCATE (my_Lrows)

         DO ispin = 1, nspins
            CALL release_intermediate_matrices(intermed_mat(ispin))
         END DO
         DEALLOCATE (intermed_mat)

         IF (my_do_gw) THEN
            DO ispin = 1, nspins
               CALL release_intermediate_matrices(intermed_mat_gw(ispin))
            END DO
            DEALLOCATE (intermed_mat_gw)
         END IF

         IF (do_bse) THEN
            CALL release_intermediate_matrices(intermed_mat_bse_ab)
            CALL release_intermediate_matrices(intermed_mat_bse_ij)
         END IF

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
         CALL dbt_pgrid_create(para_env, pdims_t3c, pgrid_t3c_overl)

         ! set up basis
         ALLOCATE (sizes_RI(natom), sizes_AO(natom))
         ALLOCATE (basis_set_ri_aux(nkind), basis_set_ao(nkind))
         CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
         CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=basis_set_ri_aux)
         CALL basis_set_list_setup(basis_set_ao, "ORB", qs_kind_set)
         CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_AO, basis=basis_set_ao)

         ! make sure we use the QS%EPS_PGF_ORB
         qs_section => section_vals_get_subs_vals(qs_env%input, "DFT%QS")
         CALL section_vals_val_get(qs_section, "EPS_PGF_ORB", n_rep_val=n_rep)
         IF (n_rep /= 0) THEN
            CALL section_vals_val_get(qs_section, "EPS_PGF_ORB", r_val=eps_pgf_orb)
         ELSE
            CALL section_vals_val_get(qs_section, "EPS_DEFAULT", r_val=eps_pgf_orb)
            eps_pgf_orb = SQRT(eps_pgf_orb)
         END IF
         eps_pgf_orb_old = dft_control%qs_control%eps_pgf_orb

         DO ibasis = 1, SIZE(basis_set_ao)
            orb_basis => basis_set_ao(ibasis)%gto_basis_set
            CALL init_interaction_radii_orb_basis(orb_basis, eps_pgf_orb)
            ri_basis => basis_set_ri_aux(ibasis)%gto_basis_set
            CALL init_interaction_radii_orb_basis(ri_basis, eps_pgf_orb)
         END DO

         cut_memory_int = qs_env%mp2_env%ri_rpa_im_time%cut_memory
         CALL create_tensor_batches(sizes_RI, cut_memory_int, starts_array_mc_int, ends_array_mc_int, &
                                    starts_array_mc_block_int, ends_array_mc_block_int)

         DEALLOCATE (starts_array_mc_int, ends_array_mc_int)

         CALL create_3c_tensor(t_3c_overl_int_template, dist_RI, dist_AO_1, dist_AO_2, pgrid_t3c_overl, &
                               sizes_RI, sizes_AO, sizes_AO, map1=[1, 2], map2=[3], &
                               name="O (RI AO | AO)")

         CALL get_qs_env(qs_env, nkind=nkind, particle_set=particle_set)
         CALL dbt_mp_environ_pgrid(pgrid_t3c_overl, pdims, pcoord)
         CALL mp_comm_t3c_2%create(pgrid_t3c_overl%mp_comm_2d, 3, pdims)
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
         END IF

         ALLOCATE (t_3c_overl_int(nimg, nimg))

         DO i = 1, SIZE(t_3c_overl_int, 1)
            DO j = 1, SIZE(t_3c_overl_int, 2)
               CALL dbt_create(t_3c_overl_int_template, t_3c_overl_int(i, j))
            END DO
         END DO

         CALL dbt_destroy(t_3c_overl_int_template)

         ! split blocks to improve load balancing for tensor contraction
         min_bsize = qs_env%mp2_env%ri_rpa_im_time%min_bsize

         CALL pgf_block_sizes(atomic_kind_set, basis_set_ao, min_bsize, sizes_AO_split)
         CALL pgf_block_sizes(atomic_kind_set, basis_set_ri_aux, min_bsize, sizes_RI_split)

         pdims_t3c = 0
         CALL dbt_pgrid_create(para_env, pdims_t3c, pgrid_t3c_M)

         ASSOCIATE (cut_memory => qs_env%mp2_env%ri_rpa_im_time%cut_memory)
            CALL create_tensor_batches(sizes_AO_split, cut_memory, starts_array_mc, ends_array_mc, &
                                       starts_array_mc_block, ends_array_mc_block)
            CALL create_tensor_batches(sizes_RI_split, cut_memory, &
                                       qs_env%mp2_env%ri_rpa_im_time%starts_array_mc_RI, &
                                       qs_env%mp2_env%ri_rpa_im_time%ends_array_mc_RI, &
                                       qs_env%mp2_env%ri_rpa_im_time%starts_array_mc_block_RI, &
                                       qs_env%mp2_env%ri_rpa_im_time%ends_array_mc_block_RI)

         END ASSOCIATE
         cut_memory = qs_env%mp2_env%ri_rpa_im_time%cut_memory

         CALL create_3c_tensor(t_3c_M, dist_RI, dist_AO_1, dist_AO_2, pgrid_t3c_M, &
                               sizes_RI_split, sizes_AO_split, sizes_AO_split, &
                               map1=[1], map2=[2, 3], &
                               name="M (RI | AO AO)")
         DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)
         CALL dbt_pgrid_destroy(pgrid_t3c_M)

         ALLOCATE (t_3c_O(SIZE(t_3c_overl_int, 1), SIZE(t_3c_overl_int, 2)))
         ALLOCATE (t_3c_O_compressed(SIZE(t_3c_overl_int, 1), SIZE(t_3c_overl_int, 2), cut_memory))
         ALLOCATE (t_3c_O_ind(SIZE(t_3c_overl_int, 1), SIZE(t_3c_overl_int, 2), cut_memory))
         CALL create_3c_tensor(t_3c_O(1, 1), dist_RI, dist_AO_1, dist_AO_2, pgrid_t3c_overl, &
                               sizes_RI_split, sizes_AO_split, sizes_AO_split, &
                               map1=[1, 2], map2=[3], &
                               name="O (RI AO | AO)")
         DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)
         CALL dbt_pgrid_destroy(pgrid_t3c_overl)

         DO i = 1, SIZE(t_3c_O, 1)
            DO j = 1, SIZE(t_3c_O, 2)
               IF (i > 1 .OR. j > 1) CALL dbt_create(t_3c_O(1, 1), t_3c_O(i, j))
            END DO
         END DO

         ! build integrals in batches and copy to optimized format
         ! note: integrals are stored in terms of atomic blocks. To avoid a memory bottleneck,
         ! integrals are calculated in batches and copied to optimized format with subatomic blocks

         DO cm = 1, cut_memory_int
            CALL build_3c_integrals(t_3c_overl_int, &
                                    qs_env%mp2_env%ri_rpa_im_time%eps_filter/2, &
                                    qs_env, &
                                    nl_3c, &
                                    int_eps=qs_env%mp2_env%ri_rpa_im_time%eps_filter/2, &
                                    basis_i=basis_set_ri_aux, &
                                    basis_j=basis_set_ao, basis_k=basis_set_ao, &
                                    potential_parameter=ri_metric, &
                                    do_kpoints=do_kpoints_cubic_RPA, &
                                    bounds_i=[starts_array_mc_block_int(cm), ends_array_mc_block_int(cm)], desymmetrize=.FALSE.)
            CALL timeset(routineN//"_copy_3c", handle4)
            ! copy integral tensor t_3c_overl_int to t_3c_O tensor optimized for contraction
            DO i = 1, SIZE(t_3c_overl_int, 1)
               DO j = 1, SIZE(t_3c_overl_int, 2)

                  CALL dbt_copy(t_3c_overl_int(i, j), t_3c_O(i, j), order=[1, 3, 2], &
                                summation=.TRUE., move_data=.TRUE.)
                  CALL dbt_clear(t_3c_overl_int(i, j))
                  CALL dbt_filter(t_3c_O(i, j), qs_env%mp2_env%ri_rpa_im_time%eps_filter/2)
               END DO
            END DO
            CALL timestop(handle4)
         END DO

         DO i = 1, SIZE(t_3c_overl_int, 1)
            DO j = 1, SIZE(t_3c_overl_int, 2)
               CALL dbt_destroy(t_3c_overl_int(i, j))
            END DO
         END DO
         DEALLOCATE (t_3c_overl_int)

         CALL timeset(routineN//"_copy_3c", handle4)
         ! desymmetrize
         CALL dbt_create(t_3c_O(1, 1), t_3c_tmp)
         DO jcell = 1, nimg
            DO kcell = 1, jcell
               CALL dbt_copy(t_3c_O(jcell, kcell), t_3c_tmp)
               CALL dbt_copy(t_3c_tmp, t_3c_O(kcell, jcell), order=[1, 3, 2], summation=.TRUE., move_data=.TRUE.)
               CALL dbt_filter(t_3c_O(kcell, jcell), qs_env%mp2_env%ri_rpa_im_time%eps_filter)
            END DO
         END DO
         DO jcell = 1, nimg
            DO kcell = jcell + 1, nimg
               CALL dbt_copy(t_3c_O(jcell, kcell), t_3c_tmp)
               CALL dbt_copy(t_3c_tmp, t_3c_O(kcell, jcell), order=[1, 3, 2], summation=.FALSE., move_data=.TRUE.)
               CALL dbt_filter(t_3c_O(kcell, jcell), qs_env%mp2_env%ri_rpa_im_time%eps_filter)
            END DO
         END DO

         CALL dbt_get_info(t_3c_O(1, 1), nfull_total=bounds_3c)
         CALL get_tensor_occupancy(t_3c_O(1, 1), nze, occ)
         memory_3c = 0.0_dp

         bounds(:, 1) = [1, bounds_3c(1)]
         bounds(:, 3) = [1, bounds_3c(3)]
         DO i = 1, SIZE(t_3c_O, 1)
            DO j = 1, SIZE(t_3c_O, 2)
               DO i_mem = 1, cut_memory
                  bounds(:, 2) = [starts_array_mc(i_mem), ends_array_mc(i_mem)]
                  CALL dbt_copy(t_3c_O(i, j), t_3c_tmp, bounds=bounds)

                  CALL alloc_containers(t_3c_O_compressed(i, j, i_mem), 1)
                  CALL compress_tensor(t_3c_tmp, t_3c_O_ind(i, j, i_mem)%ind, &
                                       t_3c_O_compressed(i, j, i_mem), &
                                       qs_env%mp2_env%ri_rpa_im_time%eps_compress, memory_3c)
               END DO
               CALL dbt_clear(t_3c_O(i, j))
            END DO
         END DO

         CALL para_env%sum(memory_3c)

         compression_factor = REAL(nze, dp)*1.0E-06*8.0_dp/memory_3c

         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="((T3,A,T66,F11.2,A4))") &
               "MEMORY_INFO| Memory for 3-center integrals (compressed):", memory_3c, ' MiB'

            WRITE (UNIT=unit_nr, FMT="((T3,A,T60,F21.2))") &
               "MEMORY_INFO| Compression factor:                  ", compression_factor
         END IF

         CALL dbt_destroy(t_3c_tmp)

         CALL timestop(handle4)

         DO ibasis = 1, SIZE(basis_set_ao)
            orb_basis => basis_set_ao(ibasis)%gto_basis_set
            CALL init_interaction_radii_orb_basis(orb_basis, eps_pgf_orb_old)
            ri_basis => basis_set_ri_aux(ibasis)%gto_basis_set
            CALL init_interaction_radii_orb_basis(ri_basis, eps_pgf_orb_old)
         END DO

         DEALLOCATE (basis_set_ri_aux, basis_set_ao)

         CALL neighbor_list_3c_destroy(nl_3c)

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
      INTEGER, INTENT(IN)                                :: ngroup, igroup

      CLASS(mp_comm_type), INTENT(IN)                    :: mp_comm
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env_sub

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_B_L'
      LOGICAL, PARAMETER                                 :: debug = .FALSE.

      INTEGER                                            :: check_proc, handle, i, iend, ii, ioff, &
                                                            istart, loc_a, loc_P, nblk_per_thread
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: block_ind_L_P, block_ind_L_R
      INTEGER, DIMENSION(1)                              :: dist_B_i, map_B_1, map_L_1, map_L_2, &
                                                            sizes_i
      INTEGER, DIMENSION(2)                              :: map_B_2, pdims_L
      INTEGER, DIMENSION(3)                              :: pdims_B
      LOGICAL                                            :: found
      INTEGER, DIMENSION(ngroup)                         :: dist_L_P, dist_L_R
      INTEGER, DIMENSION(para_env_sub%num_pe)            :: dist_B_a
      TYPE(dbt_distribution_type)                        :: dist_B, dist_L
      TYPE(dbt_pgrid_type)                               :: mp_comm_B, mp_comm_L
      TYPE(dbt_type)                                     :: tB_in, tB_in_split, tB_out, &
                                                            tB_out_split, tL, tL_split

      CALL timeset(routineN, handle)

      sizes_i(1) = SIZE(BIb_C, 3)

      ASSOCIATE (nproc => para_env_sub%num_pe, iproc => para_env_sub%mepos, iproc_glob => mp_comm%mepos)

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

         CALL dbt_pgrid_create(mp_comm, pdims_B, mp_comm_B)
         CALL dbt_pgrid_create(mp_comm, pdims_L, mp_comm_L)

         ! setup distribution vectors such that distribution matches parallel data layout of BIb_C and my_Lrows
         dist_B_i = [0]
         dist_B_a = (/(i, i=0, nproc - 1)/)
         dist_L_R = (/(MODULO(i, nproc), i=0, ngroup - 1)/) ! R index is replicated in my_Lrows, we impose a cyclic distribution
         dist_L_P = (/(i, i=0, ngroup - 1)/)

         ! create distributions and tensors
         CALL dbt_distribution_new(dist_B, mp_comm_B, dist_L_P, dist_B_a, dist_B_i)
         CALL dbt_distribution_new(dist_L, mp_comm_L, dist_L_R, dist_L_P)

         CALL dbt_create(tB_in, "(R|ai)", dist_B, map_B_1, map_B_2, sizes_L, sizes_B, sizes_i)
         CALL dbt_create(tB_out, "(P|ai)", dist_B, map_B_1, map_B_2, sizes_L, sizes_B, sizes_i)
         CALL dbt_create(tL, "(R|P)", dist_L, map_L_1, map_L_2, sizes_L, sizes_L)

         IF (debug) THEN
            ! check that tensor distribution is correct
            CALL dbt_get_stored_coordinates(tB_in, [loc_P, loc_a, 1], check_proc)
            CPASSERT(check_proc .EQ. iproc_glob)
         END IF

         ! reserve (R|ai) block
!$OMP PARALLEL DEFAULT(NONE) SHARED(tB_in,loc_P,loc_a)
         CALL dbt_reserve_blocks(tB_in, [loc_P], [loc_a], [1])
!$OMP END PARALLEL

         ! reserve (R|P) blocks
         ! in my_Lrows, R index is replicated. For (R|P), we distribute quadratic blocks cyclically over
         ! the processes in a subgroup.
         ! There are NG blocks, so each process holds at most NG/Nw+1 blocks.
         ALLOCATE (block_ind_L_R(ngroup/nproc + 1))
         ALLOCATE (block_ind_L_P(ngroup/nproc + 1))
         block_ind_L_R(:) = 0; block_ind_L_P(:) = 0
         ii = 0
         DO i = 1, ngroup
            CALL dbt_get_stored_coordinates(tL, [i, loc_P], check_proc)
            IF (check_proc == iproc_glob) THEN
               ii = ii + 1
               block_ind_L_R(ii) = i
               block_ind_L_P(ii) = loc_P
            END IF
         END DO

!TODO: Parallelize creation of block list.
!$OMP PARALLEL DEFAULT(NONE) SHARED(tL,block_ind_L_R,block_ind_L_P,ii) &
!$OMP PRIVATE(nblk_per_thread,istart,iend)
         nblk_per_thread = ii/omp_get_num_threads() + 1
         istart = omp_get_thread_num()*nblk_per_thread + 1
         iend = MIN(istart + nblk_per_thread, ii)
         CALL dbt_reserve_blocks(tL, block_ind_L_R(istart:iend), block_ind_L_P(istart:iend))
!$OMP END PARALLEL

         ! insert (R|ai) block
         CALL dbt_put_block(tB_in, [loc_P, loc_a, 1], SHAPE(BIb_C), BIb_C)

         ! insert (R|P) blocks
         ioff = 0
         DO i = 1, ngroup
            istart = ioff + 1; iend = ioff + sizes_L(i)
            ioff = ioff + sizes_L(i)
            CALL dbt_get_stored_coordinates(tL, [i, loc_P], check_proc)
            IF (check_proc == iproc_glob) THEN
               CALL dbt_put_block(tL, [i, loc_P], [sizes_L(i), sizes_L(loc_P)], my_Lrows(istart:iend, :))
            END IF
         END DO
      END ASSOCIATE

      CALL dbt_split_blocks(tB_in, tB_in_split, [blk_size(2), blk_size(1), blk_size(1)])
      CALL dbt_split_blocks(tL, tL_split, [blk_size(2), blk_size(2)])
      CALL dbt_split_blocks(tB_out, tB_out_split, [blk_size(2), blk_size(1), blk_size(1)])

      ! contract
      CALL dbt_contract(alpha=1.0_dp, tensor_1=tB_in_split, tensor_2=tL_split, &
                        beta=0.0_dp, tensor_3=tB_out_split, &
                        contract_1=[1], notcontract_1=[2, 3], &
                        contract_2=[1], notcontract_2=[2], &
                        map_1=[2, 3], map_2=[1], optimize_dist=.TRUE.)

      ! retrieve local block of contraction result (P|ai)
      CALL dbt_copy(tB_out_split, tB_out)

      CALL dbt_get_block(tB_out, [loc_P, loc_a, 1], SHAPE(BIb_C), BIb_C, found)
      CPASSERT(found)

      ! cleanup
      CALL dbt_destroy(tB_in)
      CALL dbt_destroy(tB_in_split)
      CALL dbt_destroy(tB_out)
      CALL dbt_destroy(tB_out_split)
      CALL dbt_destroy(tL)
      CALL dbt_destroy(tL_split)

      CALL dbt_distribution_destroy(dist_B)
      CALL dbt_distribution_destroy(dist_L)

      CALL dbt_pgrid_destroy(mp_comm_B)
      CALL dbt_pgrid_destroy(mp_comm_L)

      CALL timestop(handle)

   END SUBROUTINE contract_B_L

! **************************************************************************************************
!> \brief Encapsulate building of intermediate matrices matrix_ia_jnu(_beta
!>         matrix_ia_jb(_beta),fm_BIb_jb(_beta),matrix_in_jnu(for G0W0) and
!>         fm_BIb_all(for G0W0)
!> \param intermed_mat ...
!> \param mo_coeff_templ ...
!> \param size_1 ...
!> \param size_2 ...
!> \param matrix_name_2 ...
!> \param blacs_env_sub ...
!> \param para_env_sub ...
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE create_intermediate_matrices(intermed_mat, mo_coeff_templ, size_1, size_2, &
                                           matrix_name_2, blacs_env_sub, para_env_sub)

      TYPE(intermediate_matrix_type), INTENT(OUT)        :: intermed_mat
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mo_coeff_templ
      INTEGER, INTENT(IN)                                :: size_1, size_2
      CHARACTER(LEN=*), INTENT(IN)                       :: matrix_name_2
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(mp_para_env_type), POINTER                    :: para_env_sub

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_intermediate_matrices'

      INTEGER                                            :: handle, ncol_local, nfullcols_total, &
                                                            nfullrows_total, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      ! initialize and create the matrix (K|jnu)
      CALL dbcsr_create(intermed_mat%matrix_ia_jnu, template=mo_coeff_templ)

      ! Allocate Sparse matrices: (K|jb)
      CALL cp_dbcsr_m_by_n_from_template(intermed_mat%matrix_ia_jb, template=mo_coeff_templ, m=size_2, n=size_1, &
                                         sym=dbcsr_type_no_symmetry)

      ! set all to zero in such a way that the memory is actually allocated
      CALL dbcsr_set(intermed_mat%matrix_ia_jnu, 0.0_dp)
      CALL dbcsr_set(intermed_mat%matrix_ia_jb, 0.0_dp)

      ! create the analogous of matrix_ia_jb in fm type
      NULLIFY (fm_struct)
      CALL dbcsr_get_info(intermed_mat%matrix_ia_jb, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env_sub, nrow_global=nfullrows_total, &
                               ncol_global=nfullcols_total, para_env=para_env_sub)
      CALL cp_fm_create(intermed_mat%fm_BIb_jb, fm_struct, name="fm_BIb_jb_"//matrix_name_2)

      CALL copy_dbcsr_to_fm(intermed_mat%matrix_ia_jb, intermed_mat%fm_BIb_jb)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_get_info(matrix=intermed_mat%fm_BIb_jb, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      intermed_mat%max_row_col_local = MAX(nrow_local, ncol_local)
      CALL para_env_sub%max(intermed_mat%max_row_col_local)

      ALLOCATE (intermed_mat%local_col_row_info(0:intermed_mat%max_row_col_local, 2))
      intermed_mat%local_col_row_info = 0
      ! 0,1 nrows
      intermed_mat%local_col_row_info(0, 1) = nrow_local
      intermed_mat%local_col_row_info(1:nrow_local, 1) = row_indices(1:nrow_local)
      ! 0,2 ncols
      intermed_mat%local_col_row_info(0, 2) = ncol_local
      intermed_mat%local_col_row_info(1:ncol_local, 2) = col_indices(1:ncol_local)

      intermed_mat%descr = matrix_name_2

      CALL timestop(handle)

   END SUBROUTINE create_intermediate_matrices

! **************************************************************************************************
!> \brief Encapsulate ERI postprocessing: AO to MO transformation and store in B matrix.
!> \param para_env ...
!> \param mat_munu ...
!> \param intermed_mat ...
!> \param BIb_jb ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param eps_filter ...
!> \param my_B_end ...
!> \param my_B_start ...
! **************************************************************************************************
   SUBROUTINE ao_to_mo_and_store_B(para_env, mat_munu, intermed_mat, BIb_jb, &
                                   mo_coeff_o, mo_coeff_v, eps_filter, &
                                   my_B_end, my_B_start)
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu
      TYPE(intermediate_matrix_type), INTENT(INOUT)      :: intermed_mat
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: BIb_jb
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_o, mo_coeff_v
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      INTEGER, INTENT(IN)                                :: my_B_end, my_B_start

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ao_to_mo_and_store_B'

      INTEGER                                            :: handle

      CALL timeset(routineN//"_mult_"//TRIM(intermed_mat%descr), handle)

      CALL dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_o, &
                          0.0_dp, intermed_mat%matrix_ia_jnu, filter_eps=eps_filter)
      CALL dbcsr_multiply("T", "N", 1.0_dp, intermed_mat%matrix_ia_jnu, mo_coeff_v, &
                          0.0_dp, intermed_mat%matrix_ia_jb, filter_eps=eps_filter)
      CALL timestop(handle)

      CALL timeset(routineN//"_E_Ex_"//TRIM(intermed_mat%descr), handle)
      CALL copy_dbcsr_to_fm(intermed_mat%matrix_ia_jb, intermed_mat%fm_BIb_jb)

      IF (.NOT. (TRIM(intermed_mat%descr) .EQ. "bse_ab")) THEN

         CALL grep_my_integrals(para_env, intermed_mat%fm_BIb_jb, BIb_jb, intermed_mat%max_row_col_local, &
                                intermed_mat%local_col_row_info, &
                                my_B_end, my_B_start)

      END IF

      CALL timestop(handle)
   END SUBROUTINE ao_to_mo_and_store_B

! **************************************************************************************************
!> \brief ...
!> \param intermed_mat ...
! **************************************************************************************************
   SUBROUTINE release_intermediate_matrices(intermed_mat)
      TYPE(intermediate_matrix_type), INTENT(INOUT)      :: intermed_mat

      CALL dbcsr_release(intermed_mat%matrix_ia_jnu)
      CALL dbcsr_release(intermed_mat%matrix_ia_jb)
      CALL cp_fm_release(intermed_mat%fm_BIb_jb)
      DEALLOCATE (intermed_mat%local_col_row_info)

   END SUBROUTINE

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
                                                            nkp_extra, nkp_orig
      INTEGER, DIMENSION(3)                              :: nkp_grid, nkp_grid_extra, periodic
      LOGICAL                                            :: do_extrapolate_kpoints
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
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
      kpoints%full_grid = .TRUE.
      kpoints%use_real_wfn = .FALSE.
      kpoints%eps_geo = 1.e-6_dp
      nkp_grid(1:3) = qs_env%mp2_env%ri_rpa_im_time%kp_grid(1:3)
      do_extrapolate_kpoints = qs_env%mp2_env%ri_rpa_im_time%do_extrapolate_kpoints

      DO i_dim = 1, 3
         IF (periodic(i_dim) == 1) THEN
            CPASSERT(MODULO(nkp_grid(i_dim), 2) == 0)
         END IF
         IF (periodic(i_dim) == 0) THEN
            CPASSERT(nkp_grid(i_dim) == 1)
         END IF
      END DO

      nkp_orig = nkp_grid(1)*nkp_grid(2)*nkp_grid(3)/2

      IF (do_extrapolate_kpoints) THEN

         CPASSERT(qs_env%mp2_env%ri_rpa_im_time%kpoint_weights_W_method == kp_weights_W_uniform)

         DO i_dim = 1, 3
            IF (periodic(i_dim) == 1) nkp_grid_extra(i_dim) = nkp_grid(i_dim) + 2
            IF (periodic(i_dim) == 0) nkp_grid_extra(i_dim) = 1
         END DO

         qs_env%mp2_env%ri_rpa_im_time%kp_grid_extra(1:3) = nkp_grid_extra(1:3)

         nkp_extra = nkp_grid_extra(1)*nkp_grid_extra(2)*nkp_grid_extra(3)/2

      ELSE

         nkp_grid_extra(1:3) = 0
         nkp_extra = 0

      END IF

      nkp = nkp_orig + nkp_extra

      qs_env%mp2_env%ri_rpa_im_time%nkp_orig = nkp_orig
      qs_env%mp2_env%ri_rpa_im_time%nkp_extra = nkp_extra

      ALLOCATE (kpoints%xkp(3, nkp), kpoints%wkp(nkp))

      kpoints%nkp_grid(1:3) = nkp_grid(1:3)
      kpoints%nkp = nkp

      ALLOCATE (qs_env%mp2_env%ri_rpa_im_time%wkp_V(nkp))
      IF (do_extrapolate_kpoints) THEN
         kpoints%wkp(1:nkp_orig) = 1.0_dp/REAL(nkp_orig, KIND=dp) &
                                   /(1.0_dp - SQRT(REAL(nkp_extra, KIND=dp)/REAL(nkp_orig, KIND=dp)))
         kpoints%wkp(nkp_orig + 1:nkp) = 1.0_dp/REAL(nkp_extra, KIND=dp) &
                                         /(1.0_dp - SQRT(REAL(nkp_orig, KIND=dp)/REAL(nkp_extra, KIND=dp)))
         qs_env%mp2_env%ri_rpa_im_time%wkp_V(1:nkp_orig) = 0.0_dp
         qs_env%mp2_env%ri_rpa_im_time%wkp_V(nkp_orig + 1:nkp) = 1.0_dp/REAL(nkp_extra, KIND=dp)
      ELSE
         kpoints%wkp(:) = 1.0_dp/REAL(nkp, KIND=dp)
         qs_env%mp2_env%ri_rpa_im_time%wkp_V(:) = kpoints%wkp(:)
      END IF

      i = 0
      DO ix = 1, nkp_grid(1)
         DO iy = 1, nkp_grid(2)
            DO iz = 1, nkp_grid(3)

               IF (i == nkp_orig) CYCLE
               i = i + 1

               kpoints%xkp(1, i) = REAL(2*ix - nkp_grid(1) - 1, KIND=dp)/(2._dp*REAL(nkp_grid(1), KIND=dp))
               kpoints%xkp(2, i) = REAL(2*iy - nkp_grid(2) - 1, KIND=dp)/(2._dp*REAL(nkp_grid(2), KIND=dp))
               kpoints%xkp(3, i) = REAL(2*iz - nkp_grid(3) - 1, KIND=dp)/(2._dp*REAL(nkp_grid(3), KIND=dp))

            END DO
         END DO
      END DO

      DO ix = 1, nkp_grid_extra(1)
         DO iy = 1, nkp_grid_extra(2)
            DO iz = 1, nkp_grid_extra(3)

               i = i + 1
               IF (i > nkp) CYCLE

               kpoints%xkp(1, i) = REAL(2*ix - nkp_grid_extra(1) - 1, KIND=dp)/(2._dp*REAL(nkp_grid_extra(1), KIND=dp))
               kpoints%xkp(2, i) = REAL(2*iy - nkp_grid_extra(2) - 1, KIND=dp)/(2._dp*REAL(nkp_grid_extra(2), KIND=dp))
               kpoints%xkp(3, i) = REAL(2*iz - nkp_grid_extra(3) - 1, KIND=dp)/(2._dp*REAL(nkp_grid_extra(3), KIND=dp))

            END DO
         END DO
      END DO

      CALL kpoint_init_cell_index(kpoints, sab_orb, para_env, dft_control)

      CALL set_qs_env(qs_env, kpoints=kpoints)

      IF (unit_nr > 0) THEN

         IF (do_extrapolate_kpoints) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T69,3I4)") "KPOINT_INFO| K-point mesh for V (leading to Sigma^x):", nkp_grid(1:3)
            WRITE (UNIT=unit_nr, FMT="(T3,A,T69)") "KPOINT_INFO| K-point extrapolation for W^c is used (W^c leads to Sigma^c):"
            WRITE (UNIT=unit_nr, FMT="(T3,A,T69,3I4)") "KPOINT_INFO| K-point mesh 1 for W^c:", nkp_grid(1:3)
            WRITE (UNIT=unit_nr, FMT="(T3,A,T69,3I4)") "KPOINT_INFO| K-point mesh 2 for W^c:", nkp_grid_extra(1:3)
         ELSE
            WRITE (UNIT=unit_nr, FMT="(T3,A,T69,3I4)") "KPOINT_INFO| K-point mesh for V and W:", nkp_grid(1:3)
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,I6)") "KPOINT_INFO| Number of kpoints for V and W:", nkp
         END IF

         SELECT CASE (qs_env%mp2_env%ri_rpa_im_time%kpoint_weights_W_method)
         CASE (kp_weights_W_tailored)
            WRITE (UNIT=unit_nr, FMT="(T3,A,T81)") &
               "KPOINT_INFO| K-point weights for W:                                   TAILORED"
         CASE (kp_weights_W_auto)
            WRITE (UNIT=unit_nr, FMT="(T3,A,T81)") &
               "KPOINT_INFO| K-point weights for W:                                       AUTO"
         CASE (kp_weights_W_uniform)
            WRITE (UNIT=unit_nr, FMT="(T3,A,T81)") &
               "KPOINT_INFO| K-point weights for W:                                    UNIFORM"
         END SELECT

      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_kpoints

! **************************************************************************************************
!> \brief ...
!> \param para_env_sub ...
!> \param fm_BIb_jb ...
!> \param BIb_jb ...
!> \param max_row_col_local ...
!> \param local_col_row_info ...
!> \param my_B_virtual_end ...
!> \param my_B_virtual_start ...
! **************************************************************************************************
   SUBROUTINE grep_my_integrals(para_env_sub, fm_BIb_jb, BIb_jb, max_row_col_local, &
                                local_col_row_info, &
                                my_B_virtual_end, my_B_virtual_start)
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env_sub
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_BIb_jb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: BIb_jb
      INTEGER, INTENT(IN)                                :: max_row_col_local
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
            proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
            proc_receive = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)

            ! first exchange information on the local data
            rec_col_row_info = 0
            CALL para_env_sub%sendrecv(local_col_row_info, proc_send, rec_col_row_info, proc_receive)
            nrow_rec = rec_col_row_info(0, 1)
            ncol_rec = rec_col_row_info(0, 2)

            ALLOCATE (row_indices_rec(nrow_rec))
            row_indices_rec = rec_col_row_info(1:nrow_rec, 1)

            ALLOCATE (col_indices_rec(ncol_rec))
            col_indices_rec = rec_col_row_info(1:ncol_rec, 2)

            ALLOCATE (rec_BI(nrow_rec, ncol_rec))
            rec_BI = 0.0_dp

            ! then send and receive the real data
            CALL para_env_sub%sendrecv(local_BI, proc_send, rec_BI, proc_receive)

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
