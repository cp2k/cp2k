!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Rountines to calculate RI-GPW-MP2 energy using pw
!> \par History
!>      06.2012 created [Mauro Del Ben]
! **************************************************************************************************
MODULE mp2_ri_gpw
   USE ai_contraction_sphi,             ONLY: abc_contract
   USE ai_coulomb,                      ONLY: coulomb2_new
   USE ai_overlap3,                     ONLY: overlap3
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: DelBen2013,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_allocate_matrix_set, cp_dbcsr_col_block_sizes, cp_dbcsr_create, &
        cp_dbcsr_deallocate_matrix, cp_dbcsr_desymmetrize, cp_dbcsr_distribution, &
        cp_dbcsr_finalize, cp_dbcsr_get_data_size, cp_dbcsr_get_data_type, cp_dbcsr_get_info, &
        cp_dbcsr_init, cp_dbcsr_multiply, cp_dbcsr_p_type, cp_dbcsr_put_block, cp_dbcsr_release, &
        cp_dbcsr_reserve_blocks, cp_dbcsr_row_block_sizes, cp_dbcsr_set, cp_dbcsr_type, &
        dbcsr_distribution_obj, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_m_by_n_from_template
   USE cp_eri_mme_interface,            ONLY: eri_mme_2c_integrate,&
                                              eri_mme_set_params
   USE cp_fm_basic_linalg,              ONLY: cp_fm_triangular_invert
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: cp_fm_syevx
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE eri_mme_types,                   ONLY: eri_mme_param
   USE gaussian_gridlevels,             ONLY: gaussian_gridlevel
   USE input_constants,                 ONLY: do_eri_gpw,&
                                              do_eri_mme,&
                                              ri_coulomb,&
                                              ri_overlap
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE machine,                         ONLY: m_flush,&
                                              m_memory,&
                                              m_walltime
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_allgather,&
                                              mp_alltoall,&
                                              mp_comm_split_direct,&
                                              mp_max,&
                                              mp_min,&
                                              mp_sendrecv,&
                                              mp_sum,&
                                              mp_sync
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types_new,              ONLY: molecule_type
   USE mp2_ri_grad_util,                ONLY: complete_gamma
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type,&
                                              one_dim_int_array
   USE orbital_pointers,                ONLY: ncoset
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_copy,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_integrate_potential,          ONLY: integrate_pgf_product_rspace,&
                                              integrate_v_rspace
   USE qs_interactions,                 ONLY: init_interaction_radii
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: set_ks_env
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                              realspace_grid_p_type,&
                                              rs_grid_release,&
                                              rs_grid_retain
   USE rpa_im_time,                     ONLY: communicate_buffer
   USE rs_pw_interface,                 ONLY: potential_pw2rs
   USE task_list_types,                 ONLY: task_list_type
   USE util,                            ONLY: get_limit

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_ri_gpw'

   PUBLIC :: mp2_ri_gpw_compute_in, &
             mp2_ri_gpw_compute_en

CONTAINS

! **************************************************************************************************
!> \brief with ri mp2 gpw
!> \param BIb_C ...
!> \param BIb_C_gw ...
!> \param ends_array ...
!> \param ends_B_virtual ...
!> \param sizes_array ...
!> \param sizes_B_virtual ...
!> \param starts_array ...
!> \param starts_B_virtual ...
!> \param dimen_RI ...
!> \param qs_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param color_sub ...
!> \param dft_control ...
!> \param cell ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param mo_coeff ...
!> \param nmo ...
!> \param homo ...
!> \param rho_r ...
!> \param rho_g ...
!> \param pot_g ...
!> \param mat_munu ...
!> \param sab_orb_sub ...
!> \param sab_orb_all ...
!> \param pw_env_sub ...
!> \param poisson_env ...
!> \param auxbas_pw_pool ...
!> \param task_list_sub ...
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
!> \param starts_B_all ...
!> \param sizes_B_all ...
!> \param ends_B_all ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param do_im_time ...
!> \param mat_3c_overl_int ...
!> \param ri_metric ...
!> \param ri_metric_gw ...
!> \param my_Lrows ...
!> \param BIb_C_beta ...
!> \param BIb_C_gw_beta ...
!> \param ends_B_virtual_beta ...
!> \param sizes_B_virtual_beta ...
!> \param starts_B_virtual_beta ...
!> \param homo_beta ...
!> \param mo_coeff_o_beta ...
!> \param mo_coeff_v_beta ...
!> \param mo_coeff_all_beta ...
!> \param mo_coeff_gw_beta ...
!> \author Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE mp2_ri_gpw_compute_in(BIb_C, BIb_C_gw, ends_array, ends_B_virtual, &
                                    sizes_array, sizes_B_virtual, starts_array, starts_B_virtual, &
                                    dimen_RI, qs_env, para_env, para_env_sub, color_sub, dft_control, &
                                    cell, particle_set, atomic_kind_set, qs_kind_set, mo_coeff, nmo, homo, &
                                    rho_r, rho_g, pot_g, mat_munu, sab_orb_sub, sab_orb_all, pw_env_sub, &
                                    poisson_env, auxbas_pw_pool, task_list_sub, mo_coeff_o, mo_coeff_v, mo_coeff_all, &
                                    mo_coeff_gw, eps_filter, unit_nr, &
                                    mp2_memory, calc_PQ_cond_num, calc_forces, blacs_env_sub, my_do_gw, &
                                    starts_B_all, sizes_B_all, ends_B_all, gw_corr_lev_occ, gw_corr_lev_virt, &
                                    do_im_time, mat_3c_overl_int, ri_metric, ri_metric_gw, my_Lrows, &
                                    BIb_C_beta, BIb_C_gw_beta, ends_B_virtual_beta, sizes_B_virtual_beta, &
                                    starts_B_virtual_beta, homo_beta, mo_coeff_o_beta, mo_coeff_v_beta, &
                                    mo_coeff_all_beta, mo_coeff_gw_beta)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C, BIb_C_gw
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ends_array, ends_B_virtual, sizes_array, &
                                                            sizes_B_virtual, starts_array, &
                                                            starts_B_virtual
      INTEGER                                            :: dimen_RI
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER                                            :: color_sub
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      INTEGER                                            :: nmo, homo
      TYPE(pw_p_type)                                    :: rho_r, rho_g, pot_g
      TYPE(cp_dbcsr_p_type)                              :: mat_munu
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub, sab_orb_all
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(task_list_type), POINTER                      :: task_list_sub
      TYPE(cp_dbcsr_type), POINTER                       :: mo_coeff_o, mo_coeff_v, mo_coeff_all, &
                                                            mo_coeff_gw
      REAL(KIND=dp)                                      :: eps_filter
      INTEGER                                            :: unit_nr
      REAL(KIND=dp)                                      :: mp2_memory
      LOGICAL                                            :: calc_PQ_cond_num, calc_forces
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      LOGICAL                                            :: my_do_gw
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: starts_B_all, sizes_B_all, ends_B_all
      INTEGER                                            :: gw_corr_lev_occ, gw_corr_lev_virt
      LOGICAL                                            :: do_im_time
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mat_3c_overl_int
      INTEGER                                            :: ri_metric, ri_metric_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: my_Lrows
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         OPTIONAL                                        :: BIb_C_beta, BIb_C_gw_beta
      INTEGER, ALLOCATABLE, DIMENSION(:), OPTIONAL       :: ends_B_virtual_beta, &
                                                            sizes_B_virtual_beta, &
                                                            starts_B_virtual_beta
      INTEGER, OPTIONAL                                  :: homo_beta
      TYPE(cp_dbcsr_type), OPTIONAL, POINTER             :: mo_coeff_o_beta, mo_coeff_v_beta, &
                                                            mo_coeff_all_beta, mo_coeff_gw_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_gpw_compute_in', &
         routineP = moduleN//':'//routineN

      INTEGER :: add_cut_RI, color_sub_3c, cut_memory, cut_RI, eri_method, gw_corr_lev_total, &
         handle, handle2, handle3, i, i_counter, itmp(2), LLL, max_row_col_local, &
         max_row_col_local_beta, max_row_col_local_gw, my_B_all_end, my_B_all_size, &
         my_B_all_start, my_B_size, my_B_size_beta, my_B_virtual_end, my_B_virtual_end_beta, &
         my_B_virtual_start, my_B_virtual_start_beta, my_group_L_end, my_group_L_size, &
         my_group_L_start, ncol_local, nfullcols_total, nfullrows_total, ngroup, ngroup_im_time, &
         nrow_local, num_small_eigen, virtual, virtual_beta
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_group_L_ends_im_time, &
                                                            my_group_L_sizes_im_time, &
                                                            my_group_L_starts_im_time, sub_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_col_row_info, &
                                                            local_col_row_info_beta, &
                                                            local_col_row_info_gw
      LOGICAL                                            :: do_alpha_beta, do_mult_2c, do_os, &
                                                            memory_info
      REAL(KIND=dp)                                      :: cond_num, mem_for_iaK, pair_energy, &
                                                            wfn_size
      TYPE(cp_dbcsr_type) :: matrix_ia_jb, matrix_ia_jb_beta, matrix_ia_jnu, matrix_ia_jnu_beta, &
         matrix_in_jm, matrix_in_jm_beta, matrix_in_jnu, matrix_in_jnu_beta
      TYPE(cp_fm_type), POINTER                          :: fm_BIb_gw, fm_BIb_gw_beta, fm_BIb_jb, &
                                                            fm_BIb_jb_beta, fm_matrix_L, &
                                                            fm_matrix_L_RI_metric
      TYPE(cp_para_env_type), POINTER                    :: para_env_L
      TYPE(eri_mme_param), POINTER                       :: eri_param
      TYPE(pw_p_type)                                    :: psi_L

      CALL timeset(routineN, handle)

      CALL cite_reference(DelBen2013)

      do_alpha_beta = .FALSE.
      IF (PRESENT(BIb_C_beta) .AND. &
          PRESENT(ends_B_virtual_beta) .AND. &
          PRESENT(sizes_B_virtual_beta) .AND. &
          PRESENT(starts_B_virtual_beta) .AND. &
          PRESENT(homo_beta) .AND. &
          PRESENT(mo_coeff_o_beta) .AND. &
          PRESENT(mo_coeff_v_beta)) do_alpha_beta = .TRUE.

      do_os = qs_env%mp2_env%mp2_gpw%do_os

      ! imag. time only without Coulomb metric implemented
      CPASSERT(.NOT. (ri_metric == ri_coulomb .AND. do_im_time))

      virtual = nmo-homo

      ! in case of imaginary time, we do not need the intermediate matrices
      IF (.NOT. do_im_time) THEN

         CALL create_intermediate_matrices(matrix_ia_jnu, matrix_ia_jb, mo_coeff_o, virtual, homo, &
                                           fm_BIb_jb, "fm_BIb_jb", max_row_col_local, &
                                           nfullrows_total, nfullcols_total, nrow_local, &
                                           ncol_local, blacs_env_sub, &
                                           para_env_sub, local_col_row_info)

         CALL create_parallelization_arrays(para_env_sub, ends_B_virtual, sizes_B_virtual, &
                                            starts_B_virtual, my_B_virtual_start, &
                                            my_B_virtual_end, my_B_size, &
                                            virtual)

         IF (do_alpha_beta) THEN

            virtual_beta = nmo-homo_beta

            CALL create_intermediate_matrices(matrix_ia_jnu_beta, matrix_ia_jb_beta, mo_coeff_o_beta, &
                                              virtual_beta, homo_beta, &
                                              fm_BIb_jb_beta, "fm_BIb_jb_beta", &
                                              max_row_col_local_beta, &
                                              nfullrows_total, nfullcols_total, nrow_local, &
                                              ncol_local, blacs_env_sub, &
                                              para_env_sub, local_col_row_info_beta)

            CALL create_parallelization_arrays(para_env_sub, ends_B_virtual_beta, sizes_B_virtual_beta, &
                                               starts_B_virtual_beta, my_B_virtual_start_beta, &
                                               my_B_virtual_end_beta, my_B_size_beta, &
                                               virtual_beta)

         END IF

      END IF ! not for imaginary time

      ! in the case of G0W0, we need (K|nm), n,m may be occ or virt (m restricted to corrected levels)
      IF (my_do_gw) THEN

         gw_corr_lev_total = gw_corr_lev_virt+gw_corr_lev_occ

         CALL create_intermediate_matrices(matrix_in_jnu, matrix_in_jm, mo_coeff_gw, &
                                           nmo, gw_corr_lev_total, &
                                           fm_BIb_gw, "fm_BIb_gw", &
                                           max_row_col_local_gw, &
                                           nfullrows_total, nfullcols_total, nrow_local, &
                                           ncol_local, blacs_env_sub, &
                                           para_env_sub, local_col_row_info_gw)

         CALL create_parallelization_arrays(para_env_sub, ends_B_all, sizes_B_all, &
                                            starts_B_all, my_B_all_start, &
                                            my_B_all_end, my_B_all_size, &
                                            nmo)

         IF (do_alpha_beta) THEN
            ! deallocate local_col_row_info_gw, otherwise it gets twice allocated in create_intermediate_m
            DEALLOCATE (local_col_row_info_gw)
            CALL create_intermediate_matrices(matrix_in_jnu_beta, matrix_in_jm_beta, mo_coeff_gw_beta, &
                                              nmo, gw_corr_lev_total, &
                                              fm_BIb_gw_beta, "fm_BIb_gw_beta", &
                                              max_row_col_local_gw, &
                                              nfullrows_total, nfullcols_total, nrow_local, &
                                              ncol_local, blacs_env_sub, &
                                              para_env_sub, local_col_row_info_gw)

            ! we don"t need parallelization arrays for beta since the matrix sizes of B_nm^P is the same
            ! for the beta case and therefore the parallelization of beta is the same than for alpha

         END IF
      END IF

      ngroup = para_env%num_pe/para_env_sub%num_pe

      ALLOCATE (sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe-1))
      DO i = 0, para_env_sub%num_pe-1
         sub_proc_map(i) = i
         sub_proc_map(-i-1) = para_env_sub%num_pe-i-1
         sub_proc_map(para_env_sub%num_pe+i) = i
      END DO

      IF (.NOT. do_im_time) THEN

         ! maybe remove wfn_size??
         wfn_size = REAL(SIZE(rho_r%pw%cr3d), KIND=dp)
         CALL mp_max(wfn_size, para_env%group)

         ! start real calculation
         NULLIFY (psi_L%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, psi_L%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)

      END IF

      do_mult_2c = qs_env%mp2_env%mp2_gpw%do_mult_2c
      eri_method = qs_env%mp2_env%eri_method
      eri_param => qs_env%mp2_env%eri_mme_param

      ! Preparations for ERI MME method to compute 2c ERIs
      IF (qs_env%mp2_env%eri_method .EQ. do_eri_mme) THEN
         CALL eri_mme_set_params(eri_param, cell, qs_kind_set, basis_type="RI_AUX", para_env=para_env)
      ENDIF

      CALL calculate_Lmin1(qs_env, eri_method, eri_param, para_env, para_env_sub, para_env_L, mp2_memory, &
                           fm_matrix_L, ngroup, color_sub, dimen_RI, &
                           mo_coeff, dft_control, psi_L, rho_r, rho_g, pot_g, pw_env_sub, poisson_env, &
                           my_group_L_size, my_group_L_start, my_group_L_end, sab_orb_sub, &
                           sizes_array, starts_array, ends_array, calc_PQ_cond_num, cond_num, &
                           num_small_eigen, auxbas_pw_pool, do_mult_2c, ri_metric, &
                           fm_matrix_L_RI_metric, do_im_time, do_os, qs_env%mp2_env%mp2_gpw%eps_pgf_orb_S, &
                           atomic_kind_set, qs_kind_set)

      IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
         "RI_INFO| Cholesky decomposition group size:", para_env_L%num_pe
      IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
         "RI_INFO| Number of groups for auxiliary basis functions", ngroup
      IF (calc_PQ_cond_num) THEN
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T67,ES14.5)") &
            "RI_INFO| Condition number of the (P|Q):", cond_num
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
            "RI_INFO| Number of Eigenvalue of (P|Q) smaller than 10^(-3):", num_small_eigen
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
         lll = itmp(2)-itmp(1)+1
         ALLOCATE (qs_env%mp2_env%ri_grad%PQ_half(lll, my_group_L_size))
         qs_env%mp2_env%ri_grad%PQ_half(:, :) = my_Lrows(itmp(1):itmp(2), 1:my_group_L_size)
      END IF

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
            "RI_INFO| Occupied  basis set size:", homo, &
            "RI_INFO| Virtual   basis set size:", virtual, &
            "RI_INFO| Auxiliary basis set size:", dimen_RI

         mem_for_iaK = dimen_RI*REAL(homo, KIND=dp)*virtual*8.0_dp/(1024_dp**2)
         IF (do_alpha_beta) mem_for_iaK = mem_for_iaK+ &
                                          dimen_RI*REAL(homo_beta, KIND=dp)*(nmo-homo_beta)*8.0_dp/(1024_dp**2)

         WRITE (unit_nr, '(T3,A,T66,F11.2,A4)') 'RI_INFO| Total memory for (ia|K) integrals:', &
            mem_for_iaK, ' MiB'
         IF (my_do_gw) THEN
            mem_for_iaK = dimen_RI*REAL(nmo, KIND=dp)*gw_corr_lev_total*8.0_dp/(1024_dp**2)

            WRITE (unit_nr, '(T3,A,T66,F11.2,A4)') 'RI_INFO| Total memory for G0W0-(nm|K) integrals:', &
               mem_for_iaK, ' MiB'
         END IF
         CALL m_flush(unit_nr)
      ENDIF

      CALL timeset(routineN//"_alloc", handle2)
      CALL mp_sync(para_env%group) ! sync to see memory output

      ! in case we do imaginary time, we need the overlap matrix (alpha beta P)
      IF (.NOT. do_im_time) THEN

         ! array that will store the (ia|K) integrals
         ALLOCATE (BIb_C(my_group_L_size, my_B_size, homo))
         BIb_C = 0.0_dp

         IF (do_alpha_beta) THEN
            ALLOCATE (BIb_C_beta(my_group_L_size, my_B_size_beta, homo_beta))
            BIb_C_beta = 0.0_dp
         END IF

      END IF ! imag. time

      ! in the case of GW, we also need (nm|K)
      IF (my_do_gw) THEN

         ALLOCATE (BIb_C_gw(my_group_L_size, my_B_all_size, gw_corr_lev_total))
         BIb_C_gw = 0.0_dp
         IF (do_alpha_beta) THEN
            ALLOCATE (BIb_C_gw_beta(my_group_L_size, my_B_all_size, gw_corr_lev_total))
            BIb_C_gw_beta = 0.0_dp
         END IF

      END IF

      CALL timestop(handle2)

      CALL timeset(routineN//"_loop", handle2)
      i_counter = 0
      DO LLL = my_group_L_start, my_group_L_end
         i_counter = i_counter+1

         ! if we do Obara Saika for the three-center overlap integrals, exit
         IF (do_im_time) EXIT

         ! pseudo psi_L
         CALL calculate_wavefunction(mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                     qs_kind_set, cell, dft_control, particle_set, pw_env_sub, &
                                     basis_type="RI_AUX", &
                                     external_vector=my_Lrows(:, LLL-my_group_L_start+1))

         CALL timeset(routineN//"_pot", handle3)
         rho_r%pw%cr3d = psi_L%pw%cr3d

         ! in case we do Coulomb metric for RI, we need the Coulomb operator, but for RI with the
         ! overlap metric, we do not need the Coulomb operator
         IF (ri_metric == ri_coulomb) THEN
            CALL pw_transfer(rho_r%pw, rho_g%pw)
            CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw)
            CALL pw_transfer(pot_g%pw, rho_r%pw)
         END IF

         CALL pw_scale(rho_r%pw, rho_r%pw%pw_grid%dvol)
         CALL timestop(handle3)

         ! and finally (K|mu nu)
         CALL timeset(routineN//"_int", handle3)
         CALL cp_dbcsr_set(mat_munu%matrix, 0.0_dp)
         CALL integrate_v_rspace(rho_r, hmat=mat_munu, qs_env=qs_env, calculate_forces=.FALSE., compute_tau=.FALSE., gapw=.FALSE., &
                                 pw_env_external=pw_env_sub, task_list_external=task_list_sub)
         CALL timestop(handle3)

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

         IF (my_do_gw .AND. (ri_metric_gw == ri_metric)) THEN
            CALL ao_to_mo_and_store_B(para_env_sub, mat_munu, matrix_in_jnu, matrix_in_jm, &
                                      fm_BIb_gw, BIb_C_gw(i_counter, 1:my_B_all_size, 1:gw_corr_lev_total), &
                                      mo_coeff_gw, mo_coeff_all, eps_filter, max_row_col_local_gw, sub_proc_map, &
                                      local_col_row_info_gw, my_B_all_end, my_B_all_start, "gw_alpha")

            IF (do_alpha_beta) THEN
               CALL ao_to_mo_and_store_B(para_env_sub, mat_munu, matrix_in_jnu_beta, matrix_in_jm_beta, &
                                         fm_BIb_gw_beta, BIb_c_gw_beta(i_counter, 1:my_B_all_size, 1:gw_corr_lev_total), &
                                         mo_coeff_gw_beta, mo_coeff_all_beta, eps_filter, max_row_col_local_gw, sub_proc_map, &
                                         local_col_row_info_gw, my_B_all_end, my_B_all_start, "gw_beta")
            ENDIF
         ENDIF

         IF (my_do_gw .AND. (ri_metric_gw == ri_metric)) THEN
            ! transform (K|mu nu) to (K|nm), n corresponds to corrected GW levels, m is in nmo
            CALL ao_to_mo_and_store_B(para_env_sub, mat_munu, matrix_in_jnu, matrix_in_jm, &
                                      fm_BIb_gw, BIb_C_gw(i_counter, 1:my_B_all_size, 1:gw_corr_lev_total), &
                                      mo_coeff_gw, mo_coeff_all, eps_filter, max_row_col_local_gw, sub_proc_map, &
                                      local_col_row_info_gw, my_B_all_end, my_B_all_start, "gw_alpha")

            ! the same for beta
            IF (do_alpha_beta) THEN
               CALL ao_to_mo_and_store_B(para_env_sub, mat_munu, matrix_in_jnu_beta, matrix_in_jm_beta, &
                                         fm_BIb_gw_beta, BIb_C_gw_beta(i_counter, 1:my_B_all_size, 1:gw_corr_lev_total), &
                                         mo_coeff_gw_beta, mo_coeff_all_beta, eps_filter, max_row_col_local_gw, sub_proc_map, &
                                         local_col_row_info_gw, my_B_all_end, my_B_all_start, "gw_beta")

            END IF
         END IF

      END DO

      ! use possibly different RI metric for GW matrix elements B^nn_P
      IF (my_do_gw .AND. ri_metric_gw /= ri_metric) THEN

         DEALLOCATE (my_Lrows, sizes_array, starts_array, ends_array)

         ! difference to call of calculate_Lmin1 above: other RI metric for GW
         CALL calculate_Lmin1(qs_env, eri_method, eri_param, para_env, para_env_sub, para_env_L, mp2_memory, &
                              fm_matrix_L, ngroup, color_sub, dimen_RI, &
                              mo_coeff, dft_control, psi_L, rho_r, rho_g, pot_g, pw_env_sub, poisson_env, &
                              my_group_L_size, my_group_L_start, my_group_L_end, sab_orb_sub, &
                              sizes_array, starts_array, ends_array, calc_PQ_cond_num, cond_num, &
                              num_small_eigen, auxbas_pw_pool, do_mult_2c, ri_metric_gw, &
                              fm_matrix_L_RI_metric, do_im_time, do_os, qs_env%mp2_env%mp2_gpw%eps_pgf_orb_S, &
                              atomic_kind_set, qs_kind_set)

         ! replicate the necessary row of the L^{-1} matrix on each proc
         CALL grep_Lcols(para_env_L, dimen_RI, fm_matrix_L, &
                         my_group_L_start, my_group_L_end, my_group_L_size, my_Lrows)

         ! clean the L^{-1} matrix
         CALL cp_fm_release(fm_matrix_L)
         CALL cp_para_env_release(para_env_L)

         i_counter = 0
         DO LLL = my_group_L_start, my_group_L_end
            i_counter = i_counter+1

            ! pseudo psi_L
            CALL calculate_wavefunction(mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                        qs_kind_set, cell, dft_control, particle_set, pw_env_sub, &
                                        basis_type="RI_AUX", &
                                        external_vector=my_Lrows(:, LLL-my_group_L_start+1))

            CALL timeset(routineN//"_pot_trunc", handle3)
            rho_r%pw%cr3d = psi_L%pw%cr3d

            ! in case we do Coulomb metric for RI, we need the Coulomb operator, but for RI with the
            ! overlap metric, we do not need the Coulomb operator
            IF (ri_metric_gw == ri_coulomb) THEN
               CALL pw_transfer(rho_r%pw, rho_g%pw)
               CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw)
               CALL pw_transfer(pot_g%pw, rho_r%pw)
            END IF

            CALL pw_scale(rho_r%pw, rho_r%pw%pw_grid%dvol)

            CALL timestop(handle3)

            ! and finally (K|mu nu)
            CALL timeset(routineN//"_int_trunc", handle3)
            CALL cp_dbcsr_set(mat_munu%matrix, 0.0_dp)
            CALL integrate_v_rspace(rho_r, hmat=mat_munu, qs_env=qs_env, &
                                    calculate_forces=.FALSE., compute_tau=.FALSE., gapw=.FALSE., &
                                    pw_env_external=pw_env_sub, task_list_external=task_list_sub)
            CALL timestop(handle3)

            CALL ao_to_mo_and_store_B(para_env_sub, mat_munu, matrix_in_jnu, matrix_in_jm, &
                                      fm_BIb_gw, BIb_C_gw(i_counter, 1:my_B_all_size, 1:gw_corr_lev_total), &
                                      mo_coeff_gw, mo_coeff_all, eps_filter, max_row_col_local_gw, sub_proc_map, &
                                      local_col_row_info_gw, my_B_all_end, my_B_all_start, "gw_trunc")

         END DO

         ! set this flag so that the truncation of above will be overwritten the next time the poisson_env is used
         IF (my_group_L_start <= my_group_L_end) &
            poisson_env%rebuild = .TRUE.

      END IF ! use truncated operator for GW
      CALL timestop(handle2)

      IF (do_im_time) THEN

         cut_memory = qs_env%mp2_env%ri_rpa_im_time%cut_memory

         add_cut_RI = qs_env%mp2_env%ri_rpa_im_time%add_cut_RI

         ngroup_im_time = para_env%num_pe/qs_env%mp2_env%ri_rpa_im_time%group_size_3c

         color_sub_3c = para_env%mepos/qs_env%mp2_env%ri_rpa_im_time%group_size_3c

         CALL setup_group_L_im_time(my_group_L_starts_im_time, my_group_L_ends_im_time, my_group_L_sizes_im_time, &
                                    dimen_RI, ngroup_im_time, cut_memory*add_cut_RI, cut_RI, color_sub_3c, qs_env)

         qs_env%mp2_env%ri_rpa_im_time_util(1)%n_group_RI_orig = ngroup

         memory_info = qs_env%mp2_env%ri_rpa_im_time%memory_info

         CALL reserve_blocks_3c(mat_3c_overl_int, mat_munu, qs_env, nmo, dimen_RI, cut_RI, &
                                my_group_L_starts_im_time, my_group_L_ends_im_time, &
                                my_group_L_sizes_im_time, sab_orb_sub, sab_orb_all)

         CALL compute_3c_overlap_int_new(mat_3c_overl_int, mat_munu, qs_env, my_group_L_starts_im_time, &
                                         my_group_L_ends_im_time, my_group_L_sizes_im_time, &
                                         nmo, dimen_RI, cut_RI, sab_orb_sub, sab_orb_all)

         ! transfer entries of fm_matrix_L_RI_metric to RPA
         CALL grep_Lcols_nonblocking(para_env_L, dimen_RI, fm_matrix_L_RI_metric, &
                                     my_group_L_start, my_group_L_end, my_group_L_size, my_Lrows)
         CALL cp_fm_release(fm_matrix_L_RI_metric)
         CALL cp_para_env_release(para_env_L)

      ELSE

         DEALLOCATE (my_Lrows)

         CALL cp_fm_release(fm_BIb_jb)
         DEALLOCATE (local_col_row_info)

         CALL cp_dbcsr_release(matrix_ia_jnu)
         CALL cp_dbcsr_release(matrix_ia_jb)
         IF (do_alpha_beta) THEN
            CALL cp_dbcsr_release(matrix_ia_jnu_beta)
            CALL cp_dbcsr_release(matrix_ia_jb_beta)
            CALL cp_fm_release(fm_BIb_jb_beta)
            DEALLOCATE (local_col_row_info_beta)
         END IF

         CALL pw_pool_give_back_pw(auxbas_pw_pool, psi_L%pw)

      END IF

      IF (my_do_gw) THEN
         CALL cp_dbcsr_release(matrix_in_jnu)
         CALL cp_dbcsr_release(matrix_in_jm)
         CALL cp_fm_release(fm_BIb_gw)
         DEALLOCATE (local_col_row_info_gw)
         IF (do_alpha_beta) THEN
            CALL cp_dbcsr_release(matrix_in_jnu_beta)
            CALL cp_dbcsr_release(matrix_in_jm_beta)
            CALL cp_fm_release(fm_BIb_gw_beta)
         END IF
      END IF

      DEALLOCATE (sub_proc_map)

      CALL timestop(handle)

   END SUBROUTINE mp2_ri_gpw_compute_in

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
!> \param nfullrows_total ...
!> \param nfullcols_total ...
!> \param nrow_local ...
!> \param ncol_local ...
!> \param blacs_env_sub ...
!> \param para_env_sub ...
!> \param local_col_row_info ...
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE create_intermediate_matrices(matrix_ia_jnu, matrix_ia_jb, mo_coeff_templ, size_1, size_2, &
                                           fm_BIb_jb, matrix_name_2, max_row_col_local, &
                                           nfullrows_total, nfullcols_total, nrow_local, &
                                           ncol_local, blacs_env_sub, &
                                           para_env_sub, local_col_row_info)

      TYPE(cp_dbcsr_type)                                :: matrix_ia_jnu, matrix_ia_jb, &
                                                            mo_coeff_templ
      INTEGER                                            :: size_1, size_2
      TYPE(cp_fm_type), POINTER                          :: fm_BIb_jb
      CHARACTER(LEN=*)                                   :: matrix_name_2
      INTEGER                                            :: max_row_col_local, nfullrows_total, &
                                                            nfullcols_total, nrow_local, ncol_local
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_col_row_info

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_intermediate_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      ! initialize and create the matrix (K|jnu)
      CALL cp_dbcsr_init(matrix_ia_jnu)
      CALL cp_dbcsr_create(matrix_ia_jnu, template=mo_coeff_templ)

      ! Allocate Sparse matrices: (K|jb)
      CALL cp_dbcsr_init(matrix_ia_jb)
      CALL cp_dbcsr_m_by_n_from_template(matrix_ia_jb, template=mo_coeff_templ, m=size_2, n=size_1, &
                                         sym=dbcsr_type_no_symmetry)

      ! set all to zero in such a way that the memory is actually allocated
      CALL cp_dbcsr_set(matrix_ia_jnu, 0.0_dp)
      CALL cp_dbcsr_set(matrix_ia_jb, 0.0_dp)

      ! create the analogous of matrix_ia_jb in fm type
      NULLIFY (fm_BIb_jb)
      NULLIFY (fm_struct)
      CALL cp_dbcsr_get_info(matrix_ia_jb, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
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
!> \brief Encapsulate building of parallelization arrays
!> \param para_env_sub ...
!> \param ends_B ...
!> \param sizes_B ...
!> \param starts_B ...
!> \param my_B_start ...
!> \param my_B_end ...
!> \param my_B_size ...
!> \param size_1 ...
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE create_parallelization_arrays(para_env_sub, ends_B, sizes_B, starts_B, &
                                            my_B_start, my_B_end, my_B_size, size_1)

      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ends_B, sizes_B, starts_B
      INTEGER                                            :: my_B_start, my_B_end, my_B_size, size_1

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_parallelization_arrays', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iproc, itmp(2)

      CALL timeset(routineN, handle)

      ! divide the b states in the sub_group in such a way to create
      ! b_start and b_end for each proc inside the sub_group
      ALLOCATE (sizes_B(0:para_env_sub%num_pe-1))
      sizes_B = 0
      ALLOCATE (starts_B(0:para_env_sub%num_pe-1))
      starts_B = 0
      ALLOCATE (ends_B(0:para_env_sub%num_pe-1))
      ends_B = 0

      DO iproc = 0, para_env_sub%num_pe-1
         itmp = get_limit(size_1, para_env_sub%num_pe, iproc)
         !IF (itmp(2) < itmp(1)) THEN
         !   itmp(2) = itmp(1)
         !ENDIF
         starts_B(iproc) = itmp(1)
         ends_B(iproc) = itmp(2)
         sizes_B(iproc) = itmp(2)-itmp(1)+1
      END DO

      my_B_start = starts_B(para_env_sub%mepos)
      my_B_end = ends_B(para_env_sub%mepos)
      my_B_size = sizes_B(para_env_sub%mepos)

      CALL timestop(handle)

   END SUBROUTINE create_parallelization_arrays

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
      TYPE(cp_dbcsr_p_type)                              :: mat_munu
      TYPE(cp_dbcsr_type)                                :: matrix_ia_jnu, matrix_ia_jb
      TYPE(cp_fm_type), POINTER                          :: fm_BIb_jb
      REAL(KIND=dp), DIMENSION(:, :)                     :: BIb_jb
      TYPE(cp_dbcsr_type), POINTER                       :: mo_coeff_o, mo_coeff_v
      REAL(KIND=dp)                                      :: eps_filter
      INTEGER                                            :: max_row_col_local
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sub_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_col_row_info
      INTEGER                                            :: my_B_end, my_B_start
      CHARACTER(LEN=*)                                   :: descr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ao_to_mo_and_store_B'

      INTEGER                                            :: handle1, handle2

      CALL timeset(routineN//"_mult_"//descr, handle1)

      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_o, &
                             0.0_dp, matrix_ia_jnu, filter_eps=eps_filter)
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, matrix_ia_jnu, mo_coeff_v, &
                             0.0_dp, matrix_ia_jb, filter_eps=eps_filter)
      CALL timestop(handle1)

      CALL timeset(routineN//"_E_Ex_"//descr, handle2)
      CALL copy_dbcsr_to_fm(matrix_ia_jb, fm_BIb_jb)

      CALL grep_my_integrals(para_env, fm_BIb_jb, BIb_jb, max_row_col_local, &
                             sub_proc_map, local_col_row_info, &
                             my_B_end, my_B_start)

      CALL timestop(handle2)
   END SUBROUTINE ao_to_mo_and_store_B

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param eri_method ...
!> \param eri_param ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param para_env_L ...
!> \param mp2_memory ...
!> \param fm_matrix_L ...
!> \param ngroup ...
!> \param color_sub ...
!> \param dimen_RI ...
!> \param mo_coeff ...
!> \param dft_control ...
!> \param psi_L ...
!> \param rho_r ...
!> \param rho_g ...
!> \param pot_g ...
!> \param pw_env_sub ...
!> \param poisson_env ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param sab_orb_sub ...
!> \param sizes_array ...
!> \param starts_array ...
!> \param ends_array ...
!> \param calc_PQ_cond_num ...
!> \param cond_num ...
!> \param num_small_eigen ...
!> \param auxbas_pw_pool ...
!> \param do_mult_2c ...
!> \param ri_metric ...
!> \param fm_matrix_L_RI_metric ...
!> \param do_im_time ...
!> \param do_os ...
!> \param mp2_eps_pgf_orb_S ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
! **************************************************************************************************
   SUBROUTINE calculate_Lmin1(qs_env, eri_method, eri_param, para_env, para_env_sub, para_env_L, mp2_memory, &
                              fm_matrix_L, ngroup, color_sub, dimen_RI, &
                              mo_coeff, dft_control, psi_L, rho_r, rho_g, pot_g, pw_env_sub, poisson_env, &
                              my_group_L_size, my_group_L_start, my_group_L_end, sab_orb_sub, &
                              sizes_array, starts_array, ends_array, calc_PQ_cond_num, cond_num, &
                              num_small_eigen, auxbas_pw_pool, do_mult_2c, &
                              ri_metric, fm_matrix_L_RI_metric, do_im_time, do_os, mp2_eps_pgf_orb_S, &
                              atomic_kind_set, qs_kind_set)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: eri_method
      TYPE(eri_mme_param), POINTER                       :: eri_param
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub, para_env_L
      REAL(KIND=dp)                                      :: mp2_memory
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L
      INTEGER                                            :: ngroup, color_sub, dimen_RI
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_p_type)                                    :: psi_L, rho_r, rho_g, pot_g
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      INTEGER                                            :: my_group_L_size, my_group_L_start, &
                                                            my_group_L_end
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_array, starts_array, ends_array
      LOGICAL                                            :: calc_PQ_cond_num
      REAL(KIND=dp)                                      :: cond_num
      INTEGER                                            :: num_small_eigen
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      LOGICAL                                            :: do_mult_2c
      INTEGER                                            :: ri_metric
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L_RI_metric
      LOGICAL                                            :: do_im_time, do_os
      REAL(KIND=dp)                                      :: mp2_eps_pgf_orb_S
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_Lmin1', &
         routineP = moduleN//':'//routineN

      LOGICAL                                            :: do_inversion
      REAL(KIND=dp)                                      :: eps_pgf_orb_old
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_S_inv_work

      ! we have to invert the Coulomb matrix V_PQ only when doing RI with the Coulomb metric

      IF (ri_metric == ri_coulomb) THEN
         do_inversion = .TRUE.
      ELSE
         do_inversion = .FALSE.
      END IF

      ! calculate L^{-1} or L (for overlap and trunc. Coulomb metric)
      CALL compute_2c_in_and_cholesky(qs_env, eri_method, eri_param, para_env, para_env_sub, para_env_L, mp2_memory, &
                                      fm_matrix_L, ngroup, color_sub, dimen_RI, &
                                      mo_coeff, dft_control, psi_L, rho_r, rho_g, pot_g, pw_env_sub, poisson_env, &
                                      my_group_L_size, my_group_L_start, my_group_L_end, sab_orb_sub, &
                                      sizes_array, starts_array, ends_array, calc_PQ_cond_num, cond_num, &
                                      num_small_eigen, auxbas_pw_pool, do_mult_2c, do_inversion)

      IF (ri_metric == ri_overlap) THEN

         ! Obara-Saika overlap matrix
         IF (do_os) THEN

            ! re-init the radii to be able to generate pair lists with appropriate screening for overlap matrix
            eps_pgf_orb_old = dft_control%qs_control%eps_pgf_orb
            dft_control%qs_control%eps_pgf_orb = mp2_eps_pgf_orb_S
            CALL init_interaction_radii(dft_control%qs_control, atomic_kind_set, qs_kind_set)

            CALL OS_overlap_mat_and_cholesky(qs_env, fm_matrix_L_RI_metric, fm_matrix_L, para_env_L, dimen_RI)

            ! re-init the radii to previous values
            dft_control%qs_control%eps_pgf_orb = eps_pgf_orb_old
            CALL init_interaction_radii(dft_control%qs_control, atomic_kind_set, qs_kind_set)

            ! GPW overlap matrix
         ELSE

            CALL cp_para_env_release(para_env_L)
            DEALLOCATE (sizes_array, starts_array, ends_array)

            ! Calculate Chol. decomp. of matrix of RI operator, store it in fm_matrix_L_RI_metric
            CALL compute_2c_in_and_cholesky(qs_env, eri_method, eri_param, para_env, para_env_sub, para_env_L, mp2_memory, &
                                            fm_matrix_L_RI_metric, ngroup, color_sub, dimen_RI, &
                                            mo_coeff, dft_control, psi_L, rho_r, rho_g, pot_g, pw_env_sub, poisson_env, &
                                            my_group_L_size, my_group_L_start, my_group_L_end, sab_orb_sub, &
                                            sizes_array, starts_array, ends_array, calc_PQ_cond_num, cond_num, &
                                            num_small_eigen, auxbas_pw_pool, do_mult_2c=.FALSE., do_inversion=.TRUE., &
                                            external_matrix_struc=.TRUE., fm_matrix_L_extern=fm_matrix_L, &
                                            ri_metric=ri_metric)

         END IF

         NULLIFY (fm_matrix_S_inv_work)
         CALL cp_fm_create(fm_matrix_S_inv_work, fm_matrix_L_RI_metric%matrix_struct)
         CALL cp_fm_set_all(fm_matrix_S_inv_work, 0.0_dp)

         CALL cp_gemm('N', 'T', dimen_RI, dimen_RI, dimen_RI, 1.0_dp, fm_matrix_L_RI_metric, fm_matrix_L_RI_metric, &
                      0.0_dp, fm_matrix_S_inv_work)

         ! L (=fm_matrix_L_RI_metric) = S^(-1)*chol_dec(V)
         CALL cp_gemm('T', 'T', dimen_RI, dimen_RI, dimen_RI, 1.0_dp, fm_matrix_S_inv_work, fm_matrix_L, &
                      0.0_dp, fm_matrix_L_RI_metric)

         IF (do_im_time) THEN

            ! set L to identity matrix to compute the 3-center matrix elements (alpha beta | aux) with
            ! a single auxiliary basis function aux
            CALL set_identity(fm_matrix_L)

         ELSE

            CALL cp_fm_to_fm(fm_matrix_L_RI_metric, fm_matrix_L)
            CALL cp_fm_release(fm_matrix_L_RI_metric)

         END IF

         ! clean the S_inv matrix
         CALL cp_fm_release(fm_matrix_S_inv_work)

      END IF

   END SUBROUTINE calculate_Lmin1

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param fm_matrix_L_RI_metric ...
!> \param fm_matrix_L ...
!> \param para_env_L ...
!> \param dimen_RI ...
! **************************************************************************************************
   SUBROUTINE OS_overlap_mat_and_cholesky(qs_env, fm_matrix_L_RI_metric, fm_matrix_L, para_env_L, dimen_RI)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L_RI_metric, fm_matrix_L
      TYPE(cp_para_env_type), POINTER                    :: para_env_L
      INTEGER                                            :: dimen_RI

      CHARACTER(LEN=*), PARAMETER :: routineN = 'OS_overlap_mat_and_cholesky', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: comm_exchange, handle, i_global, iiB, &
                                                            info_chol, j_global, jjB, ncol_local, &
                                                            nrow_local, sub_sub_color
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s_RI_aux
      TYPE(cp_dbcsr_type), POINTER                       :: matrix_s_RI_aux_desymm
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_dummy, fm_matrix_S_global
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_exchange
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      CALL timeset(routineN, handle)

      NULLIFY (sab_orb, matrix_s_RI_aux)

      CALL get_qs_env(qs_env=qs_env, &
                      sab_orb=sab_orb, &
                      para_env=para_env, &
                      blacs_env=blacs_env, &
                      matrix_s_RI_aux=matrix_s_RI_aux)

      CALL build_overlap_matrix(qs_env%ks_env, &
                                matrix_s=matrix_s_RI_aux, &
                                basis_type_a='RI_AUX', &
                                basis_type_b='RI_AUX', &
                                sab_nl=sab_orb, &
                                qs_input=qs_env%input)

      CALL set_ks_env(qs_env%ks_env, matrix_s_RI_aux=matrix_s_RI_aux)

      NULLIFY (matrix_s_RI_aux_desymm)
      ALLOCATE (matrix_s_RI_aux_desymm)
      CALL cp_dbcsr_init(matrix_s_RI_aux_desymm)
      CALL cp_dbcsr_create(matrix_s_RI_aux_desymm, 'S_RI non_ysmm', &
                           cp_dbcsr_distribution(matrix_s_RI_aux(1)%matrix), dbcsr_type_no_symmetry, &
                           cp_dbcsr_row_block_sizes(matrix_s_RI_aux(1)%matrix), &
                           cp_dbcsr_col_block_sizes(matrix_s_RI_aux(1)%matrix), &
                           cp_dbcsr_get_data_size(matrix_s_RI_aux(1)%matrix), &
                           cp_dbcsr_get_data_type(matrix_s_RI_aux(1)%matrix))

      CALL cp_dbcsr_desymmetrize(matrix_s_RI_aux(1)%matrix, matrix_s_RI_aux_desymm)

      CALL cp_fm_create(fm_matrix_L_RI_metric, fm_matrix_L%matrix_struct)
      CALL cp_fm_set_all(fm_matrix_L_RI_metric, 0.0_dp)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=dimen_RI, &
                               ncol_global=dimen_RI, para_env=para_env)

      NULLIFY (fm_matrix_S_global)

      CALL cp_fm_create(fm_matrix_S_global, fm_struct)
      CALL cp_fm_set_all(fm_matrix_S_global, 0.0_dp)

      CALL copy_dbcsr_to_fm(matrix_s_RI_aux_desymm, fm_matrix_S_global)

      CALL cp_dbcsr_deallocate_matrix(matrix_s_RI_aux_desymm)

      CALL cp_fm_struct_release(fm_struct)

      ! create the new group for the mp_sum of the local data
      sub_sub_color = para_env_L%mepos
      CALL mp_comm_split_direct(para_env%group, comm_exchange, sub_sub_color)

      NULLIFY (para_env_exchange)
      CALL cp_para_env_create(para_env_exchange, comm_exchange)

      NULLIFY (fm_dummy)

      IF (para_env_exchange%mepos == 0) THEN
         CALL cp_fm_copy_general(fm_matrix_S_global, fm_matrix_L_RI_metric, blacs_env%group)
      ELSE
         CALL cp_fm_copy_general(fm_matrix_S_global, fm_dummy, blacs_env%group)
      END IF

      CALL cp_fm_release(fm_matrix_S_global)

      CALL mp_sum(fm_matrix_L_RI_metric%local_data, para_env_exchange%group)

      CALL cp_para_env_release(para_env_exchange)

      ! do cholesky decomposition
      CALL cp_fm_cholesky_decompose(matrix=fm_matrix_L_RI_metric, n=dimen_RI, info_out=info_chol)
      CPASSERT(info_chol == 0)

      CALL cp_fm_triangular_invert(matrix_a=fm_matrix_L_RI_metric, uplo_tr='U')

      ! clean the lower part of the L^{-1} matrix (just to not have surprises afterwards)
      CALL cp_fm_get_info(matrix=fm_matrix_L_RI_metric, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)
      DO iiB = 1, nrow_local
         i_global = row_indices(iiB)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            IF (j_global < i_global) fm_matrix_L_RI_metric%local_data(iiB, jjB) = 0.0_dp
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param eri_method ...
!> \param eri_param ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param para_env_L ...
!> \param mp2_memory ...
!> \param fm_matrix_L ...
!> \param ngroup ...
!> \param color_sub ...
!> \param dimen_RI ...
!> \param mo_coeff ...
!> \param dft_control ...
!> \param psi_L ...
!> \param rho_r ...
!> \param rho_g ...
!> \param pot_g ...
!> \param pw_env_sub ...
!> \param poisson_env ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param sab_orb_sub ...
!> \param sizes_array ...
!> \param starts_array ...
!> \param ends_array ...
!> \param calc_PQ_cond_num ...
!> \param cond_num ...
!> \param num_small_eigen ...
!> \param auxbas_pw_pool ...
!> \param do_mult_2c ...
!> \param do_inversion ...
!> \param external_matrix_struc ...
!> \param fm_matrix_L_extern ...
!> \param ri_metric ...
! **************************************************************************************************
   SUBROUTINE compute_2c_in_and_cholesky(qs_env, eri_method, eri_param, para_env, para_env_sub, para_env_L, mp2_memory, &
                                         fm_matrix_L, ngroup, color_sub, dimen_RI, &
                                         mo_coeff, dft_control, psi_L, rho_r, rho_g, pot_g, pw_env_sub, poisson_env, &
                                         my_group_L_size, my_group_L_start, my_group_L_end, sab_orb_sub, &
                                         sizes_array, starts_array, ends_array, calc_PQ_cond_num, cond_num, &
                                         num_small_eigen, auxbas_pw_pool, do_mult_2c, do_inversion, &
                                         external_matrix_struc, fm_matrix_L_extern, ri_metric)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: eri_method
      TYPE(eri_mme_param), POINTER                       :: eri_param
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub, para_env_L
      REAL(KIND=dp)                                      :: mp2_memory
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L
      INTEGER                                            :: ngroup, color_sub, dimen_RI
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_p_type)                                    :: psi_L, rho_r, rho_g, pot_g
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      INTEGER                                            :: my_group_L_size, my_group_L_start, &
                                                            my_group_L_end
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_array, starts_array, ends_array
      LOGICAL                                            :: calc_PQ_cond_num
      REAL(KIND=dp)                                      :: cond_num
      INTEGER                                            :: num_small_eigen
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      LOGICAL                                            :: do_mult_2c, do_inversion
      LOGICAL, OPTIONAL                                  :: external_matrix_struc
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_matrix_L_extern
      INTEGER, OPTIONAL                                  :: ri_metric

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_2c_in_and_cholesky', &
         routineP = moduleN//':'//routineN

      INTEGER :: best_group_size, color_exc, color_L, comm_exc, comm_exchange, comm_L, dir, &
         group_size, handle, handle2, handle3, i, i_bf, i_counter, i_counter_diff_bf, &
         i_counter_diff_bf_at_any_atom, i_counter_diff_sets, i_global, iatom, igrid_level, igroup, &
         iiB, ikind, info_chol, inode, ipgf, iproc, iset, itmp(2), j, j_global, j_ref_atom, &
         j_ref_bf, j_ref_kind, j_ref_set, jatom, jjB, jkind, jset, lb(3), LLL, location(3), m, &
         my_group_L_q_end, my_group_L_q_size, my_group_L_q_start, my_group_pot_end, &
         my_group_pot_size, my_group_pot_start, my_ri_metric, na1, na2, natom, ncoa, ncoa_a, ncoa_b
      INTEGER :: ncol_local, nkind, nrow_local, nseta, nsetb, nsgf, offset, offsetb, proc_receive, &
         proc_receive_static, proc_send, proc_send_static, proc_shift, rec_L_end, rec_L_size, &
         rec_L_start, sgfa, sgfb, strat_group_size, sub_sub_color, tp(3), ub(3)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: diffbf_atom_arbitr_atom, diffbf_noinbas, &
         diffbfatanyatom_noinbas, ends_array_pot, kind_of, proc_map, sizes_array_pot, &
         starts_array_pot, sub_ends_array, sub_sizes_array, sub_starts_array
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: atom_set_afo_diffbf, &
                                                            atom_set_afo_diffbf_atanyatom, &
                                                            diffsets_atom_set_bf, ref_matrix
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: naolimits, ref_at_se_sgf_ofs
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, la_max, la_min, lb_max, &
                                                            lb_min, npgfa, npgfb, nsgfa, nsgfb, &
                                                            row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: check_i, check_j, map_it_here, &
                                                            my_external_matrix_struc
      REAL(KIND=dp)                                      :: min_mem_for_QK, pair_energy, rab2, rac2, &
                                                            scaling, zet_tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: egen_L, f, wf_vector
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: I_ab, I_tmp, L_external_col, &
                                                            L_local_col, L_local_col_q, vac
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: v
      REAL(KIND=dp), DIMENSION(1)                        :: exp_q
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rac
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: I_tmp2, rpgfa, rpgfb, sphi_a, sphi_b, &
                                                            zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_L
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L_diag
      TYPE(cp_para_env_type), POINTER                    :: para_env_exc, para_env_exchange
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_p_type), ALLOCATABLE, DIMENSION(:)         :: rho_r_array
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v

      CALL timeset(routineN, handle)

      ! default: Coulomb metric
      my_ri_metric = ri_coulomb
      IF (PRESENT(ri_metric)) my_ri_metric = ri_metric

      my_external_matrix_struc = .FALSE.
      IF (PRESENT(external_matrix_struc)) THEN
         my_external_matrix_struc = external_matrix_struc
      END IF

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      cell=cell, molecule_set=molecule_set, particle_set=particle_set, &
                      molecule_kind_set=molecule_kind_set)

      ! blacs_env => qs_env%blacs_env

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type="RI_AUX")
         CPASSERT(ASSOCIATED(basis_set_a))
      END DO

      dimen_RI = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), nsgf=nsgf, basis_type="RI_AUX")
         dimen_RI = dimen_RI+nsgf
      END DO

      ! calculate which rows of L^{-1} to have
      ALLOCATE (sizes_array(0:ngroup-1))
      sizes_array = 0
      ALLOCATE (starts_array(0:ngroup-1))
      starts_array = 0
      ALLOCATE (ends_array(0:ngroup-1))
      ends_array = 0

      DO igroup = 0, ngroup-1
         itmp = get_limit(dimen_RI, ngroup, igroup)
         starts_array(igroup) = itmp(1)
         ends_array(igroup) = itmp(2)
         sizes_array(igroup) = itmp(2)-itmp(1)+1
      ENDDO

      my_group_L_size = sizes_array(color_sub)
      my_group_L_start = starts_array(color_sub)
      my_group_L_end = ends_array(color_sub)

      CALL timeset(routineN//"_loop_lm", handle2)

      IF (do_mult_2c) THEN

         ! write arrays with number(=label) of every basis function, its quantum number (l,m)
         ! and the reference to the first basis function with (l,m)
         ALLOCATE (ref_matrix(dimen_RI, 10))

         ! naolimits(iatom,iset,1/2)=begin/end reference entries in L_local_q
         ALLOCATE (naolimits(natom, 100, 2))

         ! array with reference (atom,set) for arbitrary (atom,set);
         ALLOCATE (ref_at_se_sgf_ofs(natom, 100, 4))

         i_counter = 0
         i_counter_diff_bf = 0
         i_counter_diff_sets = 0
         i_counter_diff_bf_at_any_atom = 0

         ! 1) build up reference matrices to get the reference q_alpha for every basis function alpha

         offset = 0

         DO iatom = 1, natom
            ikind = kind_of(iatom)
            CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type="RI_AUX")

            la_max => basis_set_a%lmax
            la_min => basis_set_a%lmin
            npgfa => basis_set_a%npgf
            nseta = basis_set_a%nset
            nsgfa => basis_set_a%nsgf_set
            zeta => basis_set_a%zet
            first_sgfa => basis_set_a%first_sgf

            DO iset = 1, nseta

               DO ipgf = 1, npgfa(iset)

                  !pseudo-m quantum number runs from 1 to 2l+1
                  DO m = 1, nsgfa(iset)
                     i_counter = i_counter+1
                     !first entry: l; only in case of la_max(iset)==la_min(iset) the shortcut will work, see IF below
                     ref_matrix(i_counter, 1) = la_max(iset)
                     !second entry: pseudo-m quantum number runs from 1 to 2l+1
                     ref_matrix(i_counter, 2) = m
                     ref_matrix(i_counter, 3) = iatom
                     ref_matrix(i_counter, 4) = iset
                     ref_matrix(i_counter, 5) = ipgf
                     ref_matrix(i_counter, 8) = first_sgfa(1, iset)

                     !Check, whether (l,m) combination already exists
                     IF (i_counter >= 2) THEN
                        DO i_bf = 1, i_counter-1
                           ! 6th row of ref_matrix contains number of first basis function with identical (l,m)
                           ! for set with la_min=/la_max at the same atom
                           IF (ref_matrix(i_bf, 1) == la_max(iset) .AND. ref_matrix(i_bf, 2) == m .AND. &
                               la_max(iset) == la_min(iset) .AND. ref_matrix(i_bf, 3) == iatom) THEN
                              ref_matrix(i_counter, 6) = i_bf

                              ! we are done with checking the basis functions in DO when a (l,m) match was found.
                              EXIT

                           END IF
                           ! 0 if no basis function has this (l,m) combination
                           ref_matrix(i_counter, 6) = 0

                        END DO !loop already existing basis functions

                        DO i_bf = 1, i_counter-1
                           ! 6th row of ref_matrix contains number of first basis function with identical (l,m)
                           ! for set with la_min=/la_max AT ANY ATOM (!!DIFFERENCE!!)
                           IF (ref_matrix(i_bf, 1) == la_max(iset) .AND. ref_matrix(i_bf, 2) == m .AND. &
                               la_max(iset) == la_min(iset)) THEN
                              ref_matrix(i_counter, 9) = i_bf

                              ! we are done with checking the basis functions in DO when a (l,m) match was found.
                              EXIT

                           END IF
                           ! 0 if no basis function has this (l,m) combination
                           ref_matrix(i_counter, 9) = 0

                        END DO !loop already existing basis functions

                     ELSE
                        ref_matrix(1, 6) = 0
                        ref_matrix(1, 9) = 0

                     END IF
                     ! for new (l,m) combination increase i_counter_diff_bf by 1 and set reference to itself
                     IF (ref_matrix(i_counter, 6) == 0) THEN
                        i_counter_diff_bf = i_counter_diff_bf+1
                        ref_matrix(i_counter, 7) = i_counter_diff_bf
                        ! for already existing (l,m) combination, 7th entry refers to number/label of first (l,m) combination
                     ELSE
                        ref_matrix(i_counter, 7) = &
                           ref_matrix(ref_matrix(i_counter, 6), 7)
                     END IF

                     ! for new (l,m) at any atom combination increase i_counter_diff_bf_at_any_atom by 1
                     !  and set reference to itself
                     IF (ref_matrix(i_counter, 9) == 0) THEN
                        i_counter_diff_bf_at_any_atom = i_counter_diff_bf_at_any_atom+1
                        ref_matrix(i_counter, 10) = i_counter_diff_bf_at_any_atom
                        ! for already existing (l,m) combination, 7th entry refers to number/label of first (l,m) combination
                     ELSE
                        ref_matrix(i_counter, 10) = &
                           ref_matrix(ref_matrix(i_counter, 9), 10)
                     END IF

                  END DO !m

               END DO !prim. gf.

               IF (ref_matrix(i_counter, 6) == 0) THEN
                  i_counter_diff_sets = i_counter_diff_sets+1
                  ref_at_se_sgf_ofs(iatom, iset, 1) = iatom
                  ref_at_se_sgf_ofs(iatom, iset, 2) = iset
                  ref_at_se_sgf_ofs(iatom, iset, 3) = first_sgfa(1, iset)

               ELSE
                  ref_at_se_sgf_ofs(iatom, iset, 1) = ref_matrix(ref_matrix(i_counter, 6), 3)
                  ref_at_se_sgf_ofs(iatom, iset, 2) = ref_matrix(ref_matrix(i_counter, 6), 4)
                  ref_at_se_sgf_ofs(iatom, iset, 3) = ref_matrix(ref_matrix(i_counter, 6), 8)
               END IF

               naolimits(iatom, iset, 1) = ref_matrix(i_counter, 7)-nsgfa(iset)+1
               naolimits(iatom, iset, 2) = ref_matrix(i_counter, 7)

               ref_at_se_sgf_ofs(iatom, iset, 4) = offset
               offset = offset+nsgfa(iset)

            END DO !sets

         END DO !atoms

         ! array with number of q_a basis function at an atom
         ALLOCATE (diffbf_noinbas(i_counter_diff_bf))

         ! array with number of a basis function with first unique (l,m) combination
         ALLOCATE (diffbfatanyatom_noinbas(i_counter_diff_bf_at_any_atom))

         ! array with number of q_a (at arbitrary atom) with q_a at specific atom
         ALLOCATE (diffbf_atom_arbitr_atom(i_counter_diff_bf))

         ! array with atom and set number for every q_a bf
         ALLOCATE (atom_set_afo_diffbf(i_counter_diff_bf, 2))

         ! array with atom and set number for every q_a bf
         ALLOCATE (atom_set_afo_diffbf_atanyatom(i_counter_diff_bf, 2))

         j = 1
         DO i = 1, dimen_RI
            IF (ref_matrix(i, 9) == 0) THEN
               diffbfatanyatom_noinbas(j) = i
               atom_set_afo_diffbf_atanyatom(j, 1) = ref_matrix(i, 3)
               atom_set_afo_diffbf_atanyatom(j, 2) = ref_matrix(i, 4)
               j = j+1
            END IF
         END DO

         j = 1
         DO i = 1, dimen_RI
            IF (ref_matrix(i, 6) == 0) THEN
               diffbf_noinbas(j) = i
               atom_set_afo_diffbf(j, 1) = ref_matrix(i, 3)
               atom_set_afo_diffbf(j, 2) = ref_matrix(i, 4)
               diffbf_atom_arbitr_atom(j) = ref_matrix(i, 10)
               j = j+1
            END IF
         END DO

         ! array with number of first set with quantum number l
         ALLOCATE (diffsets_atom_set_bf(i_counter_diff_sets, 3))

         ! first set with new quantum number l is on atom 1, set 1
         diffsets_atom_set_bf(1, 1) = 1
         diffsets_atom_set_bf(1, 2) = 1
         diffsets_atom_set_bf(1, 3) = 1

         j = 2
         DO i = 2, dimen_RI
            IF (ref_matrix(i, 6) == 0 .AND. (diffsets_atom_set_bf(j-1, 1) /= ref_matrix(i, 3) .OR. &
                                             diffsets_atom_set_bf(j-1, 2) /= ref_matrix(i, 4))) THEN
               diffsets_atom_set_bf(j, 1) = ref_matrix(i, 3)
               diffsets_atom_set_bf(j, 2) = ref_matrix(i, 4)
               diffsets_atom_set_bf(j, 3) = ref_matrix(i, 7)
               j = j+1
            END IF
         END DO

         ! For parallelization: Check, which elements of L_local_col_q
         ! are needed for setting up L_local_q

         my_group_L_q_start = naolimits(ref_matrix(my_group_L_start, 3), 1, 1)
         IF (ref_matrix(my_group_L_end, 3) /= natom) THEN
            my_group_L_q_end = naolimits(ref_matrix(my_group_L_end, 3)+1, 1, 1)-1
         ELSE
            nkind = kind_of(natom)
            CALL get_qs_kind(qs_kind=qs_kind_set(nkind), basis_set=basis_set_a, basis_type="RI_AUX")
            nseta = basis_set_a%nset
            my_group_L_q_end = naolimits(ref_matrix(my_group_L_end, 3), nseta, 2)
         END IF
         my_group_L_q_size = my_group_L_q_end-my_group_L_q_start+1

         ! calculate which potentials to calculate
         ALLOCATE (sizes_array_pot(0:ngroup-1))
         sizes_array_pot = 0
         ALLOCATE (starts_array_pot(0:ngroup-1))
         starts_array_pot = 0
         ALLOCATE (ends_array_pot(0:ngroup-1))
         ends_array_pot = 0

         DO igroup = 0, ngroup-1
            itmp = get_limit(i_counter_diff_bf_at_any_atom, ngroup, igroup)
            starts_array_pot(igroup) = itmp(1)
            ends_array_pot(igroup) = itmp(2)
            sizes_array_pot(igroup) = itmp(2)-itmp(1)+1
         ENDDO

         my_group_pot_size = sizes_array_pot(color_sub)
         my_group_pot_start = starts_array_pot(color_sub)
         my_group_pot_end = ends_array_pot(color_sub)

         DEALLOCATE (ref_matrix)

         ! 2.0) calculate potential (q_alpha|_p for first (l,m) combination at an atom
         i_counter = 0
         exp_q = 0.5_dp

         ALLOCATE (rho_r_array(i_counter_diff_bf_at_any_atom))

         ALLOCATE (wf_vector(dimen_RI))

         DO j = my_group_pot_start, my_group_pot_end
            wf_vector = 0.0_dp
            ! only take one basis function unique l,m-combination at any atom
            wf_vector(diffbfatanyatom_noinbas(j)) = 1.0_dp

            ! set exponent of aux basis set temporarily to exp_q since reference wavefunction
            ! q_a/b should have exponent exp_q. Then integration (q_a|q_b) is numerically
            ! feasible in reciprocal space
            jatom = atom_set_afo_diffbf_atanyatom(j, 1)
            jset = atom_set_afo_diffbf_atanyatom(j, 2)
            jkind = kind_of(jatom)
            CALL get_qs_kind(qs_kind=qs_kind_set(jkind), basis_set=basis_set_b, basis_type="RI_AUX")
            zetb => basis_set_b%zet
            lb_max => basis_set_b%lmax

            zet_tmp = zetb(1, jset)
            zetb(1, jset) = exp_q(1)

            ! pseudo psi_L
            CALL calculate_wavefunction(mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                        qs_kind_set, cell, dft_control, particle_set, pw_env_sub, &
                                        basis_type="RI_AUX", &
                                        external_vector=wf_vector)

            ! set back the exponent
            zetb(1, jset) = zet_tmp

            rho_r%pw%cr3d = psi_L%pw%cr3d
            CALL pw_transfer(rho_r%pw, rho_g%pw)
            CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw)
            CALL pw_transfer(pot_g%pw, rho_r%pw)
            CALL pw_scale(rho_r%pw, rho_r%pw%pw_grid%dvol)

            ! Writing array of potential copied from xc/xc_rho_set_types.F

            NULLIFY (rho_r_array(j)%pw)

            rho_r_array(j)%pw => rho_r%pw

            CALL pw_pool_create_pw(auxbas_pw_pool, rho_r_array(j)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_copy(rho_r%pw, rho_r_array(j)%pw)

         END DO

         ! Fill pw_grid and zero-components cr3d, aso. into rho_r_array for parallization
         ! If no potential was calculated because of my_group_pot_start>my_group_pot_end,
         ! calculate potential of my_group_pot_start-component to get a pw_grid
         IF (my_group_pot_start > my_group_pot_end) THEN
            wf_vector = 0.0_dp
            ! pseudo psi_L
            CALL calculate_wavefunction(mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                        qs_kind_set, cell, dft_control, particle_set, pw_env_sub, &
                                        basis_type="RI_AUX", &
                                        external_vector=wf_vector)

            rho_r%pw%cr3d = psi_L%pw%cr3d
            CALL pw_transfer(rho_r%pw, rho_g%pw)
            CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw)
            CALL pw_transfer(pot_g%pw, rho_r%pw)
            CALL pw_scale(rho_r%pw, rho_r%pw%pw_grid%dvol)

         END IF

         DEALLOCATE (wf_vector)

         ! Fill pw_grid and 0.0_dp potential in every component of rho_r_array
         DO j = 1, i_counter_diff_bf_at_any_atom
            IF (j < my_group_pot_start .OR. j > my_group_pot_end) THEN

               NULLIFY (rho_r_array(j)%pw)
               rho_r%pw%cr3d = 0.0_dp
               rho_r_array(j)%pw => rho_r%pw

               CALL pw_pool_create_pw(auxbas_pw_pool, rho_r_array(j)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_copy(rho_r%pw, rho_r_array(j)%pw)

            END IF
         END DO

         color_exc = para_env_sub%mepos
         CALL mp_comm_split_direct(para_env%group, comm_exc, color_exc)
         NULLIFY (para_env_exc)
         CALL cp_para_env_create(para_env_exc, comm_exc)

         ! collect all calculated potentials
         DO j = 1, i_counter_diff_bf_at_any_atom
            CALL mp_sum(rho_r_array(j)%pw%cr3d, para_env_exc%group)
         END DO

         CALL cp_para_env_release(para_env_exc)

         ALLOCATE (L_local_col_q(i_counter_diff_bf, my_group_L_q_size))

         L_local_col_q = 0.0_dp

         i_counter = 0

         ! 2.1) calculate (q_alpha|q_beta)_p for every (l_alpha,m_alpha), (l_beta,m_beta)
         !      at every atom. Use the potential (q_alpha|_p calculated in 2.0)

         DO j = my_group_L_q_start, my_group_L_q_end
            i_counter = i_counter+1

            jatom = atom_set_afo_diffbf(j, 1)
            jset = atom_set_afo_diffbf(j, 2)

            j_ref_bf = diffbf_atom_arbitr_atom(j)
            j_ref_atom = atom_set_afo_diffbf_atanyatom(j_ref_bf, 1)
            j_ref_set = atom_set_afo_diffbf_atanyatom(j_ref_bf, 2)

            j_ref_kind = kind_of(j_ref_atom)
            CALL get_qs_kind(qs_kind=qs_kind_set(j_ref_kind), basis_set=basis_set_b, basis_type="RI_AUX")
            zetb => basis_set_b%zet
            lb_max => basis_set_b%lmax

            NULLIFY (rs_v)
            NULLIFY (rs_descs)
            CALL pw_env_get(pw_env_sub, rs_descs=rs_descs, rs_grids=rs_v)
            DO i = 1, SIZE(rs_v)
               ! allocation
               CALL rs_grid_retain(rs_v(i)%rs_grid)
            END DO
            CALL potential_pw2rs(rs_v, rho_r_array(j_ref_bf), pw_env_sub)

            ! integrate the little bastards
            offset = 0
            ! only take one single set once for every (l,m) combination
            DO i = 1, i_counter_diff_sets
               iatom = diffsets_atom_set_bf(i, 1)
               iset = diffsets_atom_set_bf(i, 2)

               ikind = kind_of(iatom)
               CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type="RI_AUX")

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

               ra(:) = pbc(particle_set(iatom)%r, cell)- &
                       pbc(particle_set(jatom)%r, cell)+ &
                       pbc(particle_set(j_ref_atom)%r, cell)
               rab = 0.0_dp
               rab2 = 0.0_dp

               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)

               ALLOCATE (I_tmp2(ncoa, 1))
               I_tmp2 = 0.0_dp
               ALLOCATE (I_ab(nsgfa(iset), 1))
               I_ab = 0.0_dp

               igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info, MINVAL(exp_q(:)))
               map_it_here = .FALSE.
               IF (.NOT. ALL(rs_v(igrid_level)%rs_grid%desc%perd == 1)) THEN
                  DO dir = 1, 3
                     ! bounds of local grid (i.e. removing the 'wings'), if periodic
                     tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir, :), ra)*rs_v(igrid_level)%rs_grid%desc%npts(dir))
                     tp(dir) = MODULO(tp(dir), rs_v(igrid_level)%rs_grid%desc%npts(dir))
                     IF (rs_v(igrid_level)%rs_grid%desc%perd(dir) .NE. 1) THEN
                        lb(dir) = rs_v(igrid_level)%rs_grid%lb_local(dir)+rs_v(igrid_level)%rs_grid%desc%border
                        ub(dir) = rs_v(igrid_level)%rs_grid%ub_local(dir)-rs_v(igrid_level)%rs_grid%desc%border
                     ELSE
                        lb(dir) = rs_v(igrid_level)%rs_grid%lb_local(dir)
                        ub(dir) = rs_v(igrid_level)%rs_grid%ub_local(dir)
                     ENDIF
                     ! distributed grid, only map if it is local to the grid
                     location(dir) = tp(dir)+rs_v(igrid_level)%rs_grid%desc%lb(dir)
                  ENDDO
                  IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                      lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                      lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                     map_it_here = .TRUE.
                  ENDIF
               ELSE
                  ! not distributed, just a round-robin distribution over the full set of CPUs
                  IF (MODULO(offset, para_env_sub%num_pe) == para_env_sub%mepos) map_it_here = .TRUE.
               ENDIF

               offset = offset+nsgfa(iset)

               IF (map_it_here) THEN
                  DO ipgf = 1, npgfa(iset)
                     na1 = (ipgf-1)*ncoset(la_max(iset))+1
                     na2 = ipgf*ncoset(la_max(iset))
                     igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info, exp_q(1))

                     CALL integrate_pgf_product_rspace(la_max=la_max(iset), zeta=exp_q(1)/2.0_dp, la_min=la_min(iset), &
                                                       lb_max=0, zetb=exp_q(1)/2.0_dp, lb_min=0, &
                                                       ra=ra, rab=rab, rab2=rab2, &
                                                       rsgrid=rs_v(igrid_level)%rs_grid, &
                                                       cell=cell, &
                                                       cube_info=pw_env_sub%cube_info(igrid_level), &
                                                       hab=I_tmp2, &
                                                       o1=na1-1, &
                                                       o2=0, &
                                                       map_consistent=.TRUE., &
                                                       eps_gvg_rspace=dft_control%qs_control%eps_gvg_rspace, &
                                                       calculate_forces=.FALSE.)
                  END DO

                  CALL dgemm("T", "N", nsgfa(iset), 1, ncoa, &
                             1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             I_tmp2(1, 1), SIZE(I_tmp2, 1), &
                             1.0_dp, I_ab(1, 1), SIZE(I_ab, 1))

                  ! integration result has to be scaled since from basis set differing exponent exp_q
                  ! causes different normalization
                  L_local_col_q(offset-nsgfa(iset)+1:offset, i_counter) = &
                     I_ab(1:nsgfa(iset), 1)* &
                     (exp_q(1)/zeta(1, iset))**(la_max(iset)*0.5_dp+3.0_dp/4.0_dp)* &
                     (exp_q(1)/zetb(1, j_ref_set))**(lb_max(j_ref_set)*0.5_dp+3.0_dp/4.0_dp)
               END IF

               DEALLOCATE (I_tmp2)
               DEALLOCATE (I_ab)

            END DO ! q_beta)-set loop

            ! deallocation
            DO i = 1, SIZE(rs_v)
               CALL rs_grid_release(rs_v(i)%rs_grid)
            END DO

         END DO !(q_alpha| loop

         ! post cleanup
         DO j = 1, i_counter_diff_bf_at_any_atom
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r_array(j)%pw, &
                                      accept_non_compatible=.TRUE.)
         END DO

         DEALLOCATE (rho_r_array)

         DEALLOCATE (diffbf_noinbas)
         DEALLOCATE (diffbfatanyatom_noinbas)

         !  3) Distribute (q_alpha|q_beta)_p to L_local_col

         ALLOCATE (L_local_col(dimen_RI, my_group_L_size))
         L_local_col = 0.0_dp

         DO iatom = 1, natom
            ikind = kind_of(iatom)
            CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type="RI_AUX")
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

            DO iset = 1, nseta

               offsetb = 0

               DO jatom = 1, natom
                  jkind = kind_of(jatom)
                  CALL get_qs_kind(qs_kind=qs_kind_set(jkind), basis_set=basis_set_b, basis_type="RI_AUX")
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

                  DO jset = 1, nsetb

                     ! for parallelization: check, if jset is in paralleliz. range
                     IF (ref_at_se_sgf_ofs(jatom, jset, 4)+nsgfb(jset) >= my_group_L_start .AND. &
                         ref_at_se_sgf_ofs(jatom, jset, 4)+1 <= my_group_L_end) THEN

                        ! add contribution from (q_a|q_b)_p
                        scaling = 1.0_dp/(zeta(1, iset))**(la_max(iset)*0.5_dp+3.0_dp/4.0_dp)/ &
                                  (zetb(1, jset))**(lb_max(jset)*0.5_dp+3.0_dp/4.0_dp)* &
                                  exp_q(1)**(la_max(iset)*0.5_dp+3.0_dp/4.0_dp)* &
                                  exp_q(1)**(lb_max(jset)*0.5_dp+3.0_dp/4.0_dp)

                        DO i = 1, nsgfa(iset)
                           DO j = 1, nsgfb(jset)
                              ! for parallelization: check, if sgf_j is in paralleliz. range
                              IF (ref_at_se_sgf_ofs(jatom, jset, 4)+j >= my_group_L_start .AND. &
                                  ref_at_se_sgf_ofs(jatom, jset, 4)+j <= my_group_L_end) THEN

                                 ! Scale (q_a|q_b) by c_a*c_b, only implemented for uncontracted Gaussians and lmin=lmax
                                 L_local_col(ref_at_se_sgf_ofs(iatom, iset, 4)+i, &
                                             ref_at_se_sgf_ofs(jatom, jset, 4)+j-my_group_L_start+1) = &
                                    L_local_col_q(naolimits(iatom, iset, 1)+i-1, &
                                                  naolimits(jatom, jset, 1)+j-1-my_group_L_q_start+1)* &
                                    scaling

                              END IF
                           END DO
                        END DO

                        ! for s-functions alpha and beta, we need additional correction V_a(G=0)-V_q(G=0))*beta(G=0)
                        ! and V_b(G=0)-V_q(G=0))*alpha(G=0), only add once for every subgroup (for mepos==0)
                        IF (la_max(iset) == 0 .AND. lb_max(jset) == 0 .AND. para_env_sub%mepos == 0) THEN
                           L_local_col(ref_at_se_sgf_ofs(iatom, iset, 4)+1, &
                                       ref_at_se_sgf_ofs(jatom, jset, 4)-my_group_L_start+2) = &
                              L_local_col(ref_at_se_sgf_ofs(iatom, iset, 4)+1, &
                                          ref_at_se_sgf_ofs(jatom, jset, 4)-my_group_L_start+2)- &
                              Pi/(rho_r%pw%pw_grid%vol)*(4.0_dp*Pi**2.0_dp/(zeta(1, iset)*zetb(1, jset)))**0.75_dp* &
                              (2.0_dp/exp_q(1)-1.0_dp/zeta(1, iset)-1.0_dp/zetb(1, jset))
                        END IF

                     END IF

                  END DO !jset
               END DO !jatom
            END DO !iset
         END DO !iatom

         DEALLOCATE (naolimits)
         DEALLOCATE (L_local_col_q)

         ! End of calculation (q_alpha|q_beta) for every (l_alpha,m_alpha), (l_beta,m_beta) combination

         ! 4) Calculate (a|b)_s - (q_alpha|q_beta)_s and add them to (q_a|q_b)_p

         ! Loop over neighbor list
         CALL neighbor_list_iterator_create(nl_iterator, sab_orb_sub)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                   iatom=iatom, jatom=jatom, r=rac)

            ikind = kind_of(iatom)
            CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type="RI_AUX")
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

            DO iset = 1, nseta

               jkind = kind_of(jatom)
               CALL get_qs_kind(qs_kind=qs_kind_set(jkind), basis_set=basis_set_b, basis_type="RI_AUX")
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

               DO jset = 1, nsetb

                  ! for parallelization: check, if (jatom,jset) is in paralleliz. range
                  check_j = ref_at_se_sgf_ofs(jatom, jset, 4)+nsgfb(jset) >= my_group_L_start .AND. &
                            ref_at_se_sgf_ofs(jatom, jset, 4)+1 <= my_group_L_end

                  check_i = ref_at_se_sgf_ofs(iatom, iset, 4)+nsgfa(iset) >= my_group_L_start .AND. &
                            ref_at_se_sgf_ofs(iatom, iset, 4)+1 <= my_group_L_end .AND. &
                            iatom /= jatom

                  IF (check_j .OR. check_i) THEN

                     ncoa_a = npgfa(iset)*ncoset(la_max(iset))
                     ncoa_b = npgfb(jset)*ncoset(lb_max(jset))
                     sgfa = first_sgfa(1, iset)
                     sgfb = first_sgfb(1, jset)

                     ALLOCATE (I_tmp(ncoa_a, nsgfb(jset)))
                     ALLOCATE (I_tmp2(nsgfa(iset), nsgfb(jset)))

                     ALLOCATE (f(0:la_max(iset)+lb_max(jset)+2), &
                               v(ncoa_a, ncoa_b, la_max(iset)+lb_max(jset)+1), &
                               vac(ncoa_a, ncoa_b))
                     f = 0.0_dp
                     v = 0.0_dp
                     vac = 0.0_dp
                     I_tmp = 0.0_dp
                     I_tmp2 = 0.0_dp

                     rac2 = rac(1)*rac(1)+rac(2)*rac(2)+rac(3)*rac(3)

                     ! calculate spherical matrix elements (a|b)_s
                     CALL coulomb2_new(la_max=la_max(iset), npgfa=npgfa(iset), zeta=zeta(:, iset), &
                                       la_min=la_min(iset), &
                                       lc_max=lb_max(jset), npgfc=npgfb(jset), zetc=zetb(:, jset), &
                                       lc_min=lb_min(jset), rac=rac, &
                                       rac2=rac2, vac=vac, v=v, f=f)

                     ! trafo to spherical and add (a|b)_s to L_local_col
                     I_tmp(1:ncoa_a, 1:nsgfb(jset)) = MATMUL(vac(1:ncoa_a, 1:ncoa_b), sphi_b(1:ncoa_b, sgfb:(sgfb+nsgfb(jset)-1)))
                     I_tmp2(1:nsgfa(iset), 1:nsgfb(jset)) = MATMUL(TRANSPOSE(sphi_a(1:ncoa_a, sgfa:(sgfa+nsgfa(iset)-1))), &
                                                                   I_tmp(1:ncoa_a, 1:nsgfb(jset)))

                     f = 0.0_dp
                     v = 0.0_dp
                     vac = 0.0_dp

                     ! calculate -(q_a|q_b)_s
                     CALL coulomb2_new(la_max=la_max(iset), npgfa=npgfa(iset), zeta=exp_q(:), &
                                       la_min=la_min(iset), &
                                       lc_max=lb_max(jset), npgfc=npgfb(jset), zetc=exp_q(:), &
                                       lc_min=lb_min(jset), rac=rac, &
                                       rac2=rac2, vac=vac, v=v, f=f)

                     ! trafo to spherical and scale by c_a*c_b, only implemented for uncontracted Gaussians and lmin=lmax
                     I_tmp(1:ncoa_a, 1:nsgfb(jset)) = MATMUL(vac(1:ncoa_a, 1:ncoa_b), &
                                                             sphi_b(1:ncoa_b, ref_at_se_sgf_ofs(jatom, jset, 3): &
                                                                    (ref_at_se_sgf_ofs(jatom, jset, 3)+nsgfb(jset)-1)))
                     I_tmp2(1:nsgfa(iset), 1:nsgfb(jset)) = &
                        I_tmp2(1:nsgfa(iset), 1:nsgfb(jset))- &
                        MATMUL(TRANSPOSE(sphi_a(1:ncoa_a, ref_at_se_sgf_ofs(iatom, iset, 3): &
                                                (ref_at_se_sgf_ofs(iatom, iset, 3)+nsgfa(iset)-1))), &
                               I_tmp(1:ncoa_a, 1:nsgfb(jset)))/ &
                        (zeta(1, iset))**(la_max(iset)*0.5_dp+3.0_dp/4.0_dp)/ &
                        (zetb(1, jset))**(lb_max(jset)*0.5_dp+3.0_dp/4.0_dp)* &
                        exp_q(1)**(la_max(iset)*1.0_dp+3.0_dp/2.0_dp)* &
                        exp_q(1)**(lb_max(jset)*1.0_dp+3.0_dp/2.0_dp)/ &
                        (zeta(1, ref_at_se_sgf_ofs(iatom, iset, 2)))**(la_max(iset)*0.5_dp+3.0_dp/4.0_dp)/ &
                        (zetb(1, ref_at_se_sgf_ofs(jatom, jset, 2)))**(lb_max(jset)*0.5_dp+3.0_dp/4.0_dp)

                     ! add contribution from (a|b)_s-(q_a|q_b)_s
                     DO i = 1, nsgfa(iset)
                        DO j = 1, nsgfb(jset)
                           ! for parallelization: check, if sgf_j is in paralleliz. range
                           IF (ref_at_se_sgf_ofs(jatom, jset, 4)+j >= my_group_L_start .AND. &
                               ref_at_se_sgf_ofs(jatom, jset, 4)+j <= my_group_L_end .AND. &
                               check_j) THEN

                              L_local_col(ref_at_se_sgf_ofs(iatom, iset, 4)+i, &
                                          ref_at_se_sgf_ofs(jatom, jset, 4)+j-my_group_L_start+1) = &
                                 L_local_col(ref_at_se_sgf_ofs(iatom, iset, 4)+i, &
                                             ref_at_se_sgf_ofs(jatom, jset, 4)+j-my_group_L_start+1)+ &
                                 I_tmp2(i, j)
                           END IF

                           IF (ref_at_se_sgf_ofs(iatom, iset, 4)+i >= my_group_L_start .AND. &
                               ref_at_se_sgf_ofs(iatom, iset, 4)+i <= my_group_L_end .AND. &
                               check_i) THEN
                              L_local_col(ref_at_se_sgf_ofs(jatom, jset, 4)+j, &
                                          ref_at_se_sgf_ofs(iatom, iset, 4)+i-my_group_L_start+1) = &
                                 L_local_col(ref_at_se_sgf_ofs(jatom, jset, 4)+j, &
                                             ref_at_se_sgf_ofs(iatom, iset, 4)+i-my_group_L_start+1)+ &
                                 I_tmp2(i, j)
                           END IF
                        END DO
                     END DO

                     DEALLOCATE (f, v, vac, I_tmp, I_tmp2)

                  END IF !parallelization

               END DO !jset
            END DO !iset

         END DO

         CALL neighbor_list_iterator_release(nl_iterator)

         DEALLOCATE (diffsets_atom_set_bf, ref_at_se_sgf_ofs)
         DEALLOCATE (atom_set_afo_diffbf)
         DEALLOCATE (sizes_array_pot, starts_array_pot, ends_array_pot)

         CALL mp_sum(L_local_col, para_env_sub%group)

         ! end if of new algorithm
      ELSE
         ! old algorithm

         ALLOCATE (L_local_col(dimen_RI, my_group_L_size))
         L_local_col = 0.0_dp

         IF (eri_method .EQ. do_eri_mme .AND. my_ri_metric .EQ. ri_coulomb) THEN
            CALL eri_mme_2c_integrate(eri_param, para_env_sub, qs_env, &
                                      basis_type_a="RI_AUX", basis_type_b="RI_AUX", &
                                      hab=L_local_col, first_b=my_group_L_start, last_b=my_group_L_end)

         ELSEIF (eri_method .EQ. do_eri_gpw .OR. my_ri_metric .NE. ri_coulomb) THEN

            ALLOCATE (wf_vector(dimen_RI))

            i_counter = 0
            DO LLL = my_group_L_start, my_group_L_end
               i_counter = i_counter+1

               wf_vector = 0.0_dp
               wf_vector(LLL) = 1.0_dp

               ! pseudo psi_L
               CALL calculate_wavefunction(mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                           qs_kind_set, cell, dft_control, particle_set, pw_env_sub, &
                                           basis_type="RI_AUX", &
                                           external_vector=wf_vector)

               CALL timeset(routineN//"_pot_lm", handle3)
               rho_r%pw%cr3d = psi_L%pw%cr3d

               ! in case we do the overlap metric, we have no potential and hence no FT has to be carried out
               IF (my_ri_metric == ri_coulomb) THEN
                  CALL pw_transfer(rho_r%pw, rho_g%pw)
                  CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw)
                  CALL pw_transfer(pot_g%pw, rho_r%pw)
               END IF

               CALL pw_scale(rho_r%pw, rho_r%pw%pw_grid%dvol)

               NULLIFY (rs_v)
               NULLIFY (rs_descs)
               CALL pw_env_get(pw_env_sub, rs_descs=rs_descs, rs_grids=rs_v)
               DO i = 1, SIZE(rs_v)
                  CALL rs_grid_retain(rs_v(i)%rs_grid)
               END DO
               CALL potential_pw2rs(rs_v, rho_r, pw_env_sub)

               CALL timestop(handle3)

               ! integrate the little bastards
               offset = 0
               DO iatom = 1, natom
                  ikind = kind_of(iatom)
                  CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, basis_type="RI_AUX")

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

                  ra(:) = pbc(particle_set(iatom)%r, cell)
                  rab = 0.0_dp
                  rab2 = 0.0_dp

                  DO iset = 1, nseta
                     ncoa = npgfa(iset)*ncoset(la_max(iset))
                     sgfa = first_sgfa(1, iset)

                     ALLOCATE (I_tmp2(ncoa, 1))
                     I_tmp2 = 0.0_dp
                     ALLOCATE (I_ab(nsgfa(iset), 1))
                     I_ab = 0.0_dp

                     igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info, MINVAL(zeta(:, iset)))
                     map_it_here = .FALSE.
                     IF (.NOT. ALL(rs_v(igrid_level)%rs_grid%desc%perd == 1)) THEN
                        DO dir = 1, 3
                           ! bounds of local grid (i.e. removing the 'wings'), if periodic
                           tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir, :), ra)*rs_v(igrid_level)%rs_grid%desc%npts(dir))
                           tp(dir) = MODULO(tp(dir), rs_v(igrid_level)%rs_grid%desc%npts(dir))
                           IF (rs_v(igrid_level)%rs_grid%desc%perd(dir) .NE. 1) THEN
                              lb(dir) = rs_v(igrid_level)%rs_grid%lb_local(dir)+rs_v(igrid_level)%rs_grid%desc%border
                              ub(dir) = rs_v(igrid_level)%rs_grid%ub_local(dir)-rs_v(igrid_level)%rs_grid%desc%border
                           ELSE
                              lb(dir) = rs_v(igrid_level)%rs_grid%lb_local(dir)
                              ub(dir) = rs_v(igrid_level)%rs_grid%ub_local(dir)
                           ENDIF
                           ! distributed grid, only map if it is local to the grid
                           location(dir) = tp(dir)+rs_v(igrid_level)%rs_grid%desc%lb(dir)
                        ENDDO
                        IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                            lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                            lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                           map_it_here = .TRUE.
                        ENDIF
                     ELSE
                        ! not distributed, just a round-robin distribution over the full set of CPUs
                        IF (MODULO(offset, para_env_sub%num_pe) == para_env_sub%mepos) map_it_here = .TRUE.
                     ENDIF

                     offset = offset+nsgfa(iset)

                     IF (map_it_here) THEN
                        DO ipgf = 1, npgfa(iset)
                           sgfa = first_sgfa(1, iset)

                           na1 = (ipgf-1)*ncoset(la_max(iset))+1
                           na2 = ipgf*ncoset(la_max(iset))
                           igrid_level = gaussian_gridlevel(pw_env_sub%gridlevel_info, zeta(ipgf, iset))

                           CALL integrate_pgf_product_rspace( &
                              la_max=la_max(iset), zeta=zeta(ipgf, iset)/2.0_dp, la_min=la_min(iset), &
                              lb_max=0, zetb=zeta(ipgf, iset)/2.0_dp, lb_min=0, &
                              ra=ra, rab=rab, rab2=rab2, &
                              rsgrid=rs_v(igrid_level)%rs_grid, &
                              cell=cell, &
                              cube_info=pw_env_sub%cube_info(igrid_level), &
                              hab=I_tmp2, &
                              o1=na1-1, &
                              o2=0, &
                              map_consistent=.TRUE., &
                              eps_gvg_rspace=dft_control%qs_control%eps_gvg_rspace, &
                              calculate_forces=.FALSE.)
                        END DO

                        CALL dgemm("T", "N", nsgfa(iset), 1, ncoa, &
                                   1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                   I_tmp2(1, 1), SIZE(I_tmp2, 1), &
                                   1.0_dp, I_ab(1, 1), SIZE(I_ab, 1))

                        L_local_col(offset-nsgfa(iset)+1:offset, i_counter) = I_ab(1:nsgfa(iset), 1)
                     END IF

                     DEALLOCATE (I_tmp2)
                     DEALLOCATE (I_ab)

                  END DO
               END DO

               DO i = 1, SIZE(rs_v)
                  CALL rs_grid_release(rs_v(i)%rs_grid)
               END DO

            END DO

            DEALLOCATE (wf_vector)

            CALL mp_sum(L_local_col, para_env_sub%group)

         ELSE
            CPABORT("unknown ERI method")
         ENDIF
      END IF !end old code

      CALL timestop(handle2)

      ! split the total number of proc in a subgroup of the size of ~1/10 of the
      ! total num of proc
      best_group_size = para_env%num_pe

      strat_group_size = MAX(1, para_env%num_pe/10)

      min_mem_for_QK = REAL(dimen_RI, KIND=dp)*dimen_RI*3.0_dp*8.0_dp/1024_dp/1024_dp

      group_size = strat_group_size-1
      DO iproc = strat_group_size, para_env%num_pe
         group_size = group_size+1
         ! check that group_size is a multiple of sub_group_size and a divisor of
         ! the total num of proc
         IF (MOD(para_env%num_pe, group_size) /= 0 .OR. MOD(group_size, para_env_sub%num_pe) /= 0) CYCLE

         ! check for memory
         IF (REAL(group_size, KIND=dp)*mp2_memory < min_mem_for_QK) CYCLE

         best_group_size = group_size
         EXIT
      END DO

      ! create the L group
      color_L = para_env%mepos/best_group_size
      CALL mp_comm_split_direct(para_env%group, comm_L, color_L)
      NULLIFY (para_env_L)
      CALL cp_para_env_create(para_env_L, comm_L)

      ! create the blacs_L
      NULLIFY (blacs_env_L)
      CALL cp_blacs_env_create(blacs_env=blacs_env_L, para_env=para_env_L)

      ! now create the exchange group (for communication only between members not belonging to the
      ! same group
      sub_sub_color = para_env_sub%mepos
      CALL mp_comm_split_direct(para_env_L%group, comm_exchange, sub_sub_color)
      NULLIFY (para_env_exchange)
      CALL cp_para_env_create(para_env_exchange, comm_exchange)

      ! crate the proc maps
      ALLOCATE (proc_map(-para_env_exchange%num_pe:2*para_env_exchange%num_pe-1))
      DO i = 0, para_env_exchange%num_pe-1
         proc_map(i) = i
         proc_map(-i-1) = para_env_exchange%num_pe-i-1
         proc_map(para_env_exchange%num_pe+i) = i
      END DO

      ! create the information array
      ALLOCATE (sub_sizes_array(0:para_env_exchange%num_pe-1))
      sub_sizes_array = 0
      ALLOCATE (sub_starts_array(0:para_env_exchange%num_pe-1))
      sub_starts_array = 0
      ALLOCATE (sub_ends_array(0:para_env_exchange%num_pe-1))
      sub_ends_array = 0

      sub_sizes_array(para_env_exchange%mepos) = my_group_L_size
      sub_starts_array(para_env_exchange%mepos) = my_group_L_start
      sub_ends_array(para_env_exchange%mepos) = my_group_L_end

      CALL mp_sum(sub_sizes_array, para_env_exchange%group)
      CALL mp_sum(sub_starts_array, para_env_exchange%group)
      CALL mp_sum(sub_ends_array, para_env_exchange%group)

      ! create the full matrix L defined in the L group
      NULLIFY (fm_matrix_L)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env_L, nrow_global=dimen_RI, &
                               ncol_global=dimen_RI, para_env=para_env_L)
      IF (my_external_matrix_struc) THEN
         CALL cp_fm_create(fm_matrix_L, fm_matrix_L_extern%matrix_struct, name="fm_matrix_L")
      ELSE
         CALL cp_fm_create(fm_matrix_L, fm_struct, name="fm_matrix_L")
      END IF

      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_set_all(matrix=fm_matrix_L, alpha=0.0_dp)

      CALL cp_fm_get_info(matrix=fm_matrix_L, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         IF (j_global >= my_group_L_start .AND. j_global <= my_group_L_end) THEN
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)
               fm_matrix_L%local_data(iiB, jjB) = L_local_col(i_global, j_global-my_group_L_start+1)
            END DO
         END IF
      END DO

      proc_send_static = proc_map(para_env_exchange%mepos+1)
      proc_receive_static = proc_map(para_env_exchange%mepos-1)

      DO proc_shift = 1, para_env_exchange%num_pe-1
         proc_send = proc_map(para_env_exchange%mepos+proc_shift)
         proc_receive = proc_map(para_env_exchange%mepos-proc_shift)

         rec_L_size = sub_sizes_array(proc_receive)
         rec_L_start = sub_starts_array(proc_receive)
         rec_L_end = sub_ends_array(proc_receive)

         ALLOCATE (L_external_col(dimen_RI, rec_L_size))
         L_external_col = 0.0_dp

         CALL mp_sendrecv(L_local_col, proc_send_static, L_external_col, proc_receive_static, para_env_exchange%group)

         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            IF (j_global >= rec_L_start .AND. j_global <= rec_L_end) THEN
               DO iiB = 1, nrow_local
                  i_global = row_indices(iiB)
                  fm_matrix_L%local_data(iiB, jjB) = L_external_col(i_global, j_global-rec_L_start+1)
               END DO
            END IF
         END DO

         DEALLOCATE (L_local_col)
         ALLOCATE (L_local_col(dimen_RI, rec_L_size))
         L_local_col(:, :) = L_external_col

         DEALLOCATE (L_external_col)

      END DO

      DEALLOCATE (L_local_col)

      ! free the old exchange group stuff
      DEALLOCATE (proc_map)
      CALL cp_para_env_release(para_env_exchange)

      DEALLOCATE (sub_sizes_array)
      DEALLOCATE (sub_starts_array)
      DEALLOCATE (sub_ends_array)

      ! create the new group for the mp_sum of the local data
      sub_sub_color = para_env_L%mepos
      CALL mp_comm_split_direct(para_env%group, comm_exchange, sub_sub_color)
      NULLIFY (para_env_exchange)
      CALL cp_para_env_create(para_env_exchange, comm_exchange)

      CALL mp_sum(fm_matrix_L%local_data, para_env_exchange%group)

      CALL cp_para_env_release(para_env_exchange)

      cond_num = 1.0_dp
      num_small_eigen = 0
      IF (calc_PQ_cond_num) THEN
         ! calculate the condition number of the (P|Q) matrix
         ! create a copy of the matrix
         NULLIFY (fm_matrix_L_diag)
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env_L, nrow_global=dimen_RI, &
                                  ncol_global=dimen_RI, para_env=para_env_L)
         IF (my_external_matrix_struc) THEN
            CALL cp_fm_create(fm_matrix_L_diag, fm_matrix_L_extern%matrix_struct, name="fm_matrix_L_diag")
         ELSE
            CALL cp_fm_create(fm_matrix_L_diag, fm_struct, name="fm_matrix_L_diag")
         END IF

         CALL cp_fm_struct_release(fm_struct)

         CALL cp_fm_set_all(matrix=fm_matrix_L_diag, alpha=0.0_dp)

         CALL cp_fm_to_fm(source=fm_matrix_L, destination=fm_matrix_L_diag)

         ALLOCATE (egen_L(dimen_RI))

         egen_L = 0.0_dp
         CALL cp_fm_syevx(matrix=fm_matrix_L_diag, eigenvalues=egen_L)

         num_small_eigen = 0
         DO iiB = 1, dimen_RI
            IF (ABS(egen_L(iiB)) < 0.001_dp) num_small_eigen = num_small_eigen+1
         END DO

         cond_num = MAXVAL(ABS(egen_L))/MINVAL(ABS(egen_L))

         CALL cp_fm_release(fm_matrix_L_diag)

         DEALLOCATE (egen_L)
      END IF

      ! do cholesky decomposition
      CALL cp_fm_cholesky_decompose(matrix=fm_matrix_L, n=dimen_RI, info_out=info_chol)
      CPASSERT(info_chol == 0)

      IF (do_inversion) THEN
         CALL cp_fm_triangular_invert(matrix_a=fm_matrix_L, uplo_tr='U')
      END IF

      ! clean the lower part of the L^{-1} matrix (just to not have surprises afterwards)
      CALL cp_fm_get_info(matrix=fm_matrix_L, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)
      DO iiB = 1, nrow_local
         i_global = row_indices(iiB)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            IF (j_global < i_global) fm_matrix_L%local_data(iiB, jjB) = 0.0_dp
         END DO
      END DO

      ! release blacs_env
      CALL cp_blacs_env_release(blacs_env_L)

      CALL timestop(handle)

   END SUBROUTINE compute_2c_in_and_cholesky

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
      INTEGER                                            :: dimen_RI
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L
      INTEGER                                            :: my_group_L_start, my_group_L_end, &
                                                            my_group_L_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: my_Lrows

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

      ! proc_map, vector that replicate the processor numbers also
      ! for negative and positive number > num_pe
      ! needed to know which is the processor, to respect to another one,
      ! for a given shift
      ALLOCATE (proc_map(-para_env%num_pe:2*para_env%num_pe-1))
      DO iiB = 0, para_env%num_pe-1
         proc_map(iiB) = iiB
         proc_map(-iiB-1) = para_env%num_pe-iiB-1
         proc_map(para_env%num_pe+iiB) = iiB
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
               my_Lrows(i_global, j_global-my_group_L_start+1) = local_L(iiB, jjB)
            END DO
         END IF
      END DO

      proc_send_static = proc_map(para_env%mepos+1)
      proc_receive_static = proc_map(para_env%mepos-1)

      CALL timeset(routineN//"_comm", handle2)

      DO proc_shift = 1, para_env%num_pe-1
         proc_send = proc_map(para_env%mepos+proc_shift)
         proc_receive = proc_map(para_env%mepos-proc_shift)

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
                  my_Lrows(i_global, j_global-my_group_L_start+1) = rec_L(iiB, jjB)
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

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param para_env_L ...
!> \param dimen_RI ...
!> \param fm_matrix_L ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param my_group_L_size ...
!> \param my_Lrows ...
!> \author Jan Wilhelm, replace old grep_Lcols by a new one with nonblocking MPI
! **************************************************************************************************
   SUBROUTINE grep_Lcols_nonblocking(para_env_L, dimen_RI, fm_matrix_L, &
                                     my_group_L_start, my_group_L_end, my_group_L_size, my_Lrows)
      TYPE(cp_para_env_type), POINTER                    :: para_env_L
      INTEGER                                            :: dimen_RI
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L
      INTEGER                                            :: my_group_L_start, my_group_L_end, &
                                                            my_group_L_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: my_Lrows

      CHARACTER(LEN=*), PARAMETER :: routineN = 'grep_Lcols_nonblocking', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: end_idx, handle, i_element, i_global, &
                                                            iiB, imepos_L, j_global, jjB, L_end, &
                                                            L_start, ncol_local, nrow_local, &
                                                            start_idx
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: offset_buffer_send, sizes_rec, sizes_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: L_start_end
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: local_L, local_L_internal
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      ALLOCATE (my_Lrows(dimen_RI, my_group_L_size))
      my_Lrows = 0.0_dp

      ALLOCATE (L_start_end(0:para_env_L%num_pe-1, 2))
      L_start_end = 0

      L_start_end(para_env_L%mepos, 1) = my_group_L_start
      L_start_end(para_env_L%mepos, 2) = my_group_L_end

      ! get the L ranges from the other members of the para_env_L
      CALL mp_sum(L_start_end, para_env_L%group)

      CALL cp_fm_get_info(matrix=fm_matrix_L, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices, &
                          local_data=local_L_internal)

      ALLOCATE (local_L(nrow_local, ncol_local))
      local_L = local_L_internal(1:nrow_local, 1:ncol_local)

      ! since fm_matrix_L is fully replicated in every para_env_L, we only communicate inside para_env_L
      ALLOCATE (sizes_rec(0:para_env_L%num_pe-1))
      ALLOCATE (sizes_send(0:para_env_L%num_pe-1))
      sizes_rec = 0
      sizes_send = 0

      ! write info how many elements a process receives
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO imepos_L = 0, para_env_L%num_pe-1
            L_start = L_start_end(imepos_L, 1)
            L_end = L_start_end(imepos_L, 2)
            IF (j_global >= L_start .AND. j_global <= L_end) THEN
               sizes_send(imepos_L) = sizes_send(imepos_L)+nrow_local
            END IF
         END DO
      END DO

      IF (para_env_L%num_pe > 1) THEN
         CALL mp_alltoall(sizes_send, sizes_rec, 1, para_env_L%group)
      ELSE
         sizes_rec(0) = sizes_send(0)
      END IF

      ALLOCATE (buffer_rec(0:para_env_L%num_pe-1))
      ALLOCATE (buffer_send(0:para_env_L%num_pe-1))

      ! allocate data message and corresponding indices
      DO imepos_L = 0, para_env_L%num_pe-1

         ALLOCATE (buffer_rec(imepos_L)%msg(sizes_rec(imepos_L)))
         buffer_rec(imepos_L)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos_L)%msg(sizes_send(imepos_L)))
         buffer_send(imepos_L)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos_L)%indx(sizes_rec(imepos_L), 2))
         buffer_rec(imepos_L)%indx = 0

         ALLOCATE (buffer_send(imepos_L)%indx(sizes_send(imepos_L), 2))
         buffer_send(imepos_L)%indx = 0

      END DO

      ALLOCATE (offset_buffer_send(0:para_env_L%num_pe-1))
      offset_buffer_send = 0

      ! write matrix elements and info
      DO jjB = 1, ncol_local

         j_global = col_indices(jjB)

         DO imepos_L = 0, para_env_L%num_pe-1

            L_start = L_start_end(imepos_L, 1)
            L_end = L_start_end(imepos_L, 2)

            IF (j_global >= L_start .AND. j_global <= L_end) THEN

               start_idx = offset_buffer_send(imepos_L)+1
               end_idx = offset_buffer_send(imepos_L)+nrow_local

               DO iiB = 1, nrow_local

                  i_global = row_indices(iiB)

                  buffer_send(imepos_L)%indx(start_idx+iiB-1, 1) = i_global
                  buffer_send(imepos_L)%indx(start_idx+iiB-1, 2) = j_global

               END DO

               buffer_send(imepos_L)%msg(start_idx:end_idx) = local_L(1:nrow_local, jjB)

               offset_buffer_send(imepos_L) = offset_buffer_send(imepos_L)+nrow_local

            END IF

         END DO

      END DO

      ALLOCATE (req_array(1:para_env_L%num_pe, 4))

      CALL communicate_buffer(para_env_L, sizes_rec, sizes_send, buffer_rec, buffer_send, req_array)

      ! accumulate data on L
      DO imepos_L = 0, para_env_L%num_pe-1
         DO i_element = 1, sizes_rec(imepos_L)
            i_global = buffer_rec(imepos_L)%indx(i_element, 1)
            j_global = buffer_rec(imepos_L)%indx(i_element, 2)

            my_Lrows(i_global, j_global-my_group_L_start+1) = buffer_rec(imepos_L)%msg(i_element)
         END DO
      END DO

      DEALLOCATE (offset_buffer_send, req_array)

      DO imepos_L = 0, para_env_L%num_pe-1
         DEALLOCATE (buffer_send(imepos_L)%msg)
         DEALLOCATE (buffer_send(imepos_L)%indx)
      END DO

      DO imepos_L = 0, para_env_L%num_pe-1
         DEALLOCATE (buffer_rec(imepos_L)%msg)
         DEALLOCATE (buffer_rec(imepos_L)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)

      DEALLOCATE (sizes_rec, sizes_send)
      DEALLOCATE (local_L)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param fm_matrix_L ...
! **************************************************************************************************
   SUBROUTINE set_identity(fm_matrix_L)
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L

      CHARACTER(LEN=*), PARAMETER :: routineN = 'set_identity', routineP = moduleN//':'//routineN

      INTEGER                                            :: i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL cp_fm_get_info(matrix=fm_matrix_L, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)
      DO iiB = 1, nrow_local
         i_global = row_indices(iiB)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            IF (j_global == i_global) THEN
               fm_matrix_L%local_data(iiB, jjB) = 1.0_dp
            ELSE
               fm_matrix_L%local_data(iiB, jjB) = 0.0_dp
            END IF
         END DO
      END DO

   END SUBROUTINE set_identity

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
      REAL(KIND=dp), DIMENSION(:, :)                     :: BIb_jb
      INTEGER                                            :: max_row_col_local
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_col_row_info
      INTEGER                                            :: my_B_virtual_end, my_B_virtual_start

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
               BIb_jb(j_global-my_B_virtual_start+1, i_global) = fm_BIb_jb%local_data(iiB, jjB)
            END DO
         END IF
      END DO

      DEALLOCATE (row_indices_rec)
      DEALLOCATE (col_indices_rec)

      IF (para_env_sub%num_pe > 1) THEN
         ALLOCATE (local_BI(nrow_rec, ncol_rec))
         local_BI(1:nrow_rec, 1:ncol_rec) = fm_BIb_jb%local_data(1:nrow_rec, 1:ncol_rec)

         DO proc_shift = 1, para_env_sub%num_pe-1
            proc_send = proc_map(para_env_sub%mepos+proc_shift)
            proc_receive = proc_map(para_env_sub%mepos-proc_shift)

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
                     BIb_jb(j_global-my_B_virtual_start+1, i_global) = rec_BI(iiB, jjB)
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
!> \brief ...
!> \param Emp2 ...
!> \param Emp2_Cou ...
!> \param Emp2_EX ...
!> \param BIb_C ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param color_sub ...
!> \param ends_array ...
!> \param ends_B_virtual ...
!> \param sizes_array ...
!> \param sizes_B_virtual ...
!> \param starts_array ...
!> \param starts_B_virtual ...
!> \param Eigenval ...
!> \param nmo ...
!> \param homo ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param calc_forces ...
!> \param calc_ex ...
!> \param open_shell_SS ...
!> \param BIb_C_beta ...
!> \param homo_beta ...
!> \param Eigenval_beta ...
!> \param ends_B_virtual_beta ...
!> \param sizes_B_virtual_beta ...
!> \param starts_B_virtual_beta ...
! **************************************************************************************************
   SUBROUTINE mp2_ri_gpw_compute_en(Emp2, Emp2_Cou, Emp2_EX, BIb_C, mp2_env, para_env, para_env_sub, color_sub, &
                                    ends_array, ends_B_virtual, sizes_array, sizes_B_virtual, starts_array, starts_B_virtual, &
                                    Eigenval, nmo, homo, dimen_RI, unit_nr, calc_forces, calc_ex, &
                                    open_shell_SS, BIb_C_beta, homo_beta, Eigenval_beta, &
                                    ends_B_virtual_beta, sizes_B_virtual_beta, starts_B_virtual_beta)
      REAL(KIND=dp)                                      :: Emp2, Emp2_Cou, Emp2_EX
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER                                            :: color_sub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ends_array, ends_B_virtual, sizes_array, &
                                                            sizes_B_virtual, starts_array, &
                                                            starts_B_virtual
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: nmo, homo, dimen_RI, unit_nr
      LOGICAL, INTENT(IN)                                :: calc_forces
      LOGICAL                                            :: calc_ex
      LOGICAL, OPTIONAL                                  :: open_shell_SS
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         OPTIONAL                                        :: BIb_C_beta
      INTEGER, OPTIONAL                                  :: homo_beta
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: Eigenval_beta
      INTEGER, ALLOCATABLE, DIMENSION(:), OPTIONAL       :: ends_B_virtual_beta, &
                                                            sizes_B_virtual_beta, &
                                                            starts_B_virtual_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_gpw_compute_en', &
         routineP = moduleN//':'//routineN

      INTEGER :: a, a_global, b, b_global, best_block_size, best_integ_group_size, block_size, &
         comm_exchange, comm_P, end_point, handle, handle2, handle3, iiB, ij_counter, &
         ij_counter_send, ij_index, integ_group_size, irep, jjB, Lend_pos, Lstart_pos, &
         max_ij_pairs, min_integ_group_size, my_B_size, my_B_size_beta, my_B_virtual_end, &
         my_B_virtual_end_beta, my_B_virtual_start, my_B_virtual_start_beta, my_block_size, &
         my_group_L_end, my_group_L_size, my_group_L_size_orig, my_group_L_start, my_homo_beta, &
         my_i, my_ij_pairs, my_j, my_new_group_L_size, my_num_dgemm_call, ngroup, num_IJ_blocks
      INTEGER :: num_integ_group, pos_integ_group, proc_receive, proc_send, proc_shift, &
         rec_B_size, rec_B_virtual_end, rec_B_virtual_start, rec_L_size, send_B_size, &
         send_B_virtual_end, send_B_virtual_start, send_block_size, send_i, send_ij_index, send_j, &
         start_point, sub_P_color, sub_sub_color, total_ij_pairs, virtual, virtual_beta
      INTEGER, ALLOCATABLE, DIMENSION(:) :: integ_group_pos2color_sub, new_sizes_array, &
         num_ij_pairs, proc_map, proc_map_rep, sizes_array_orig, sub_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ij_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: ranges_info_array
      LOGICAL                                            :: my_alpha_alpha_case, my_alpha_beta_case, &
                                                            my_beta_beta_case, my_open_shell_SS
      REAL(KIND=dp) :: actual_flop_rate, amp_fac, mem_for_aK, mem_for_comm, mem_for_iaK, &
         mem_for_rep, mem_min, mem_per_group, mem_real, my_flop_rate, null_mat_rec(2, 2, 2), &
         null_mat_send(2, 2, 2), sym_fac, t_end, t_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: external_ab, external_i_aL, local_ab, &
                                                            local_ba, t_ab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_ia_Q, B_ia_Q_beta, BI_C_rec, &
                                                            local_i_aL, local_j_aL, Y_i_aP, &
                                                            Y_j_aP, Y_j_aP_beta
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange, para_env_P, &
                                                            para_env_rep

      CALL timeset(routineN, handle)

      my_open_shell_SS = .FALSE.
      IF (PRESENT(open_shell_SS)) my_open_shell_SS = open_shell_SS

      ! t_ab = amp_fac*(:,a|:,b)-(:,b|:,a)
      IF (calc_forces) amp_fac = 2.0_dp
      ! If we calculate the gradient we need to distinguish
      ! between alpha-alpha and beta-beta cases for UMP2

      my_alpha_alpha_case = .FALSE.
      my_beta_beta_case = .FALSE.
      my_alpha_beta_case = .FALSE.
      IF (calc_forces) THEN
         IF (my_open_shell_SS) THEN
            amp_fac = 1.0_dp
            IF ((.NOT. ALLOCATED(mp2_env%ri_grad%P_ij)) &
                .AND. (.NOT. ALLOCATED(mp2_env%ri_grad%P_ab))) THEN
               my_alpha_alpha_case = .TRUE.
               amp_fac = 1.0_dp
            ELSE
               IF ((.NOT. ALLOCATED(mp2_env%ri_grad%P_ij_beta)) &
                   .AND. (.NOT. ALLOCATED(mp2_env%ri_grad%P_ab_beta))) THEN
                  my_beta_beta_case = .TRUE.
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      my_alpha_beta_case = .FALSE.
      IF (PRESENT(BIb_C_beta) .AND. &
          PRESENT(ends_B_virtual_beta) .AND. &
          PRESENT(sizes_B_virtual_beta) .AND. &
          PRESENT(starts_B_virtual_beta) .AND. &
          PRESENT(homo_beta) .AND. &
          PRESENT(Eigenval_beta)) THEN
         my_alpha_beta_case = .TRUE.
         my_alpha_alpha_case = .FALSE.
      ENDIF

      IF (my_alpha_beta_case) amp_fac = 1.0_dp

      virtual = nmo-homo
      IF (my_alpha_beta_case) virtual_beta = nmo-homo_beta

      CALL mp2_ri_get_sizes( &
         mp2_env, para_env, para_env_sub, ends_array, ends_B_virtual, sizes_array, sizes_B_virtual, &
         starts_array, starts_B_virtual, homo, dimen_RI, unit_nr, color_sub, best_block_size, best_integ_group_size, block_size, &
         integ_group_size, min_integ_group_size, my_B_size, my_B_virtual_end, my_B_virtual_start, my_group_L_size, &
         my_group_L_start, my_group_L_end, ngroup, num_IJ_blocks, num_integ_group, pos_integ_group, virtual, my_alpha_beta_case, &
         my_open_shell_SS, mem_for_aK, mem_for_comm, mem_for_iaK, mem_for_rep, mem_min, mem_per_group, mem_real)

      IF (my_alpha_beta_case) THEN
         my_B_virtual_start_beta = starts_B_virtual_beta(para_env_sub%mepos)
         my_B_virtual_end_beta = ends_B_virtual_beta(para_env_sub%mepos)
         my_B_size_beta = sizes_B_virtual_beta(para_env_sub%mepos)
         my_homo_beta = homo_beta
      ELSE
         my_B_virtual_start_beta = my_B_virtual_start
         my_B_virtual_end_beta = my_B_virtual_end
         my_B_size_beta = my_B_size
         my_homo_beta = homo
      END IF

      ! now create a group that contains all the proc that have the same virtual starting point
      ! in the integ group
      ! sub_sub_color=para_env_sub%mepos
      CALL mp2_ri_create_group( &
         BIb_C, para_env, para_env_sub, homo, color_sub, &
         sizes_array, calc_forces, &
         comm_exchange, integ_group_size, my_B_size, iiB, my_group_L_end, &
         my_group_L_size, my_group_L_size_orig, my_group_L_start, my_new_group_L_size, &
         sub_sub_color, integ_group_pos2color_sub, new_sizes_array, proc_map, proc_map_rep, sizes_array_orig, &
         sub_proc_map, ranges_info_array, para_env_exchange, para_env_rep, num_integ_group)

      ! *****************************************************************
      ! **********  REPLICATION-BLOCKED COMMUNICATION SCHEME  ***********
      ! *****************************************************************
      ! introduce block size, the number of occupied orbitals has to be a
      ! multiple of the block size

      ! Calculate the maximum number of ij pairs that have to be computed
      ! among groups
      CALL mp2_ri_communication(my_alpha_beta_case, total_ij_pairs, homo, homo_beta, num_IJ_blocks, &
                                block_size, ngroup, ij_map, color_sub, my_ij_pairs, my_open_shell_SS, unit_nr)

      ALLOCATE (num_ij_pairs(0:para_env_exchange%num_pe-1))
      num_ij_pairs = 0
      num_ij_pairs(para_env_exchange%mepos) = my_ij_pairs
      CALL mp_sum(num_ij_pairs, para_env_exchange%group)

      max_ij_pairs = MAXVAL(num_ij_pairs)

      ! start real stuff
      IF (.NOT. my_alpha_beta_case) THEN
         CALL mp2_ri_allocate(local_ab, t_ab, mp2_env, homo, virtual, dimen_RI, my_B_size, &
                              block_size, my_B_size_beta, my_group_L_size, local_i_aL, &
                              local_j_aL, calc_forces, Y_i_aP, Y_j_aP, &
                              my_alpha_beta_case, &
                              my_beta_beta_case)
      ELSE
         CALL mp2_ri_allocate(local_ab, t_ab, mp2_env, homo, virtual, dimen_RI, my_B_size, &
                              block_size, my_B_size_beta, my_group_L_size, local_i_aL, &
                              local_j_aL, calc_forces, Y_i_aP, Y_j_aP_beta, &
                              my_alpha_beta_case, &
                              my_beta_beta_case, local_ba, virtual_beta)
      ENDIF

      CALL timeset(routineN//"_RI_loop", handle2)
      null_mat_rec = 0.0_dp
      null_mat_send = 0.0_dp
      Emp2 = 0.0_dp
      Emp2_Cou = 0.0_dp
      Emp2_EX = 0.0_dp
      my_num_dgemm_call = 0
      my_flop_rate = 0.0_dp
      DO ij_index = 1, max_ij_pairs

         IF (ij_index <= my_ij_pairs) THEN
            ! We have work to do
            ij_counter = (ij_index-MIN(1, color_sub))*ngroup+color_sub
            my_i = ij_map(ij_counter, 1)
            my_j = ij_map(ij_counter, 2)
            my_block_size = ij_map(ij_counter, 3)

            local_i_aL = 0.0_dp
            ! local_i_aL(my_group_L_start:my_group_L_end,1:my_B_size)=BIb_C(1:my_group_L_size,1:my_B_size,my_i)
            DO irep = 0, num_integ_group-1
               Lstart_pos = ranges_info_array(1, irep, para_env_exchange%mepos)
               Lend_pos = ranges_info_array(2, irep, para_env_exchange%mepos)
               start_point = ranges_info_array(3, irep, para_env_exchange%mepos)
               end_point = ranges_info_array(4, irep, para_env_exchange%mepos)

               local_i_aL(Lstart_pos:Lend_pos, 1:my_B_size, 1:my_block_size) = &
                  BIb_C(start_point:end_point, 1:my_B_size, my_i:my_i+my_block_size-1)
            END DO

            local_j_aL = 0.0_dp
            ! local_j_aL(my_group_L_start:my_group_L_end,1:my_B_size)=BIb_C(1:my_group_L_size,1:my_B_size,my_j)
            DO irep = 0, num_integ_group-1
               Lstart_pos = ranges_info_array(1, irep, para_env_exchange%mepos)
               Lend_pos = ranges_info_array(2, irep, para_env_exchange%mepos)
               start_point = ranges_info_array(3, irep, para_env_exchange%mepos)
               end_point = ranges_info_array(4, irep, para_env_exchange%mepos)

               IF (.NOT. my_alpha_beta_case) THEN
                  local_j_aL(Lstart_pos:Lend_pos, 1:my_B_size, 1:my_block_size) = &
                     BIb_C(start_point:end_point, 1:my_B_size, my_j:my_j+my_block_size-1)
               ELSE
                  local_j_aL(Lstart_pos:Lend_pos, 1:my_B_size_beta, 1:my_block_size) = &
                     BIb_C_beta(start_point:end_point, 1:my_B_size_beta, my_j:my_j+my_block_size-1)
               END IF
            END DO

            ! collect data from other proc
            CALL timeset(routineN//"_comm", handle3)
            DO proc_shift = 1, para_env_exchange%num_pe-1
               proc_send = proc_map(para_env_exchange%mepos+proc_shift)
               proc_receive = proc_map(para_env_exchange%mepos-proc_shift)

               send_ij_index = num_ij_pairs(proc_send)

               rec_L_size = sizes_array(proc_receive)
               ALLOCATE (BI_C_rec(rec_L_size, MAX(my_B_size, my_B_size_beta), my_block_size))

               IF (ij_index <= send_ij_index) THEN
                  ! ij_counter_send=(ij_index-MIN(1,proc_send))*ngroup+proc_send
                  ij_counter_send = (ij_index-MIN(1, integ_group_pos2color_sub(proc_send)))*ngroup+ &
                                    integ_group_pos2color_sub(proc_send)
                  send_i = ij_map(ij_counter_send, 1)
                  send_j = ij_map(ij_counter_send, 2)
                  send_block_size = ij_map(ij_counter_send, 3)

                  ! occupied i
                  BI_C_rec = 0.0_dp
                  CALL mp_sendrecv(BIb_C(1:my_group_L_size, 1:my_B_size, send_i:send_i+send_block_size-1), proc_send, &
                                   BI_C_rec(1:rec_L_size, 1:my_B_size, 1:my_block_size), proc_receive, &
                                   para_env_exchange%group)
                  DO irep = 0, num_integ_group-1
                     Lstart_pos = ranges_info_array(1, irep, proc_receive)
                     Lend_pos = ranges_info_array(2, irep, proc_receive)
                     start_point = ranges_info_array(3, irep, proc_receive)
                     end_point = ranges_info_array(4, irep, proc_receive)

                     local_i_aL(Lstart_pos:Lend_pos, 1:my_B_size, 1:my_block_size) = &
                        BI_C_rec(start_point:end_point, 1:my_B_size, 1:my_block_size)

                  END DO

                  ! occupied j
                  BI_C_rec = 0.0_dp
                  IF (.NOT. my_alpha_beta_case) THEN
                     CALL mp_sendrecv(BIb_C(1:my_group_L_size, 1:my_B_size, send_j:send_j+send_block_size-1), proc_send, &
                                      BI_C_rec(1:rec_L_size, 1:my_B_size, 1:my_block_size), proc_receive, &
                                      para_env_exchange%group)
                  ELSE
                     CALL mp_sendrecv(BIb_C_beta(1:my_group_L_size, 1:my_B_size_beta, send_j:send_j+send_block_size-1), proc_send, &
                                      BI_C_rec(1:rec_L_size, 1:my_B_size_beta, 1:my_block_size), proc_receive, &
                                      para_env_exchange%group)
                  END IF

                  DO irep = 0, num_integ_group-1
                     Lstart_pos = ranges_info_array(1, irep, proc_receive)
                     Lend_pos = ranges_info_array(2, irep, proc_receive)
                     start_point = ranges_info_array(3, irep, proc_receive)
                     end_point = ranges_info_array(4, irep, proc_receive)

                     IF (.NOT. my_alpha_beta_case) THEN
                        local_j_aL(Lstart_pos:Lend_pos, 1:my_B_size, 1:my_block_size) = &
                           BI_C_rec(start_point:end_point, 1:my_B_size, 1:my_block_size)
                     ELSE
                        local_j_aL(Lstart_pos:Lend_pos, 1:my_B_size_beta, 1:my_block_size) = &
                           BI_C_rec(start_point:end_point, 1:my_B_size_beta, 1:my_block_size)
                     END IF

                  END DO

               ELSE
                  ! we send the null matrix while we know that we have to receive something

                  ! occupied i
                  BI_C_rec = 0.0_dp
                  CALL mp_sendrecv(null_mat_send, proc_send, &
                                   BI_C_rec(1:rec_L_size, 1:my_B_size, 1:my_block_size), proc_receive, &
                                   para_env_exchange%group)

                  DO irep = 0, num_integ_group-1
                     Lstart_pos = ranges_info_array(1, irep, proc_receive)
                     Lend_pos = ranges_info_array(2, irep, proc_receive)
                     start_point = ranges_info_array(3, irep, proc_receive)
                     end_point = ranges_info_array(4, irep, proc_receive)

                     local_i_aL(Lstart_pos:Lend_pos, 1:my_B_size, 1:my_block_size) = &
                        BI_C_rec(start_point:end_point, 1:my_B_size, 1:my_block_size)

                  END DO

                  ! occupied j
                  BI_C_rec = 0.0_dp
                  IF (.NOT. my_alpha_beta_case) THEN
                     CALL mp_sendrecv(null_mat_send, proc_send, &
                                      BI_C_rec(1:rec_L_size, 1:my_B_size, 1:my_block_size), proc_receive, &
                                      para_env_exchange%group)
                  ELSE
                     CALL mp_sendrecv(null_mat_send, proc_send, &
                                      BI_C_rec(1:rec_L_size, 1:my_B_size_beta, 1:my_block_size), proc_receive, &
                                      para_env_exchange%group)
                  END IF
                  DO irep = 0, num_integ_group-1
                     Lstart_pos = ranges_info_array(1, irep, proc_receive)
                     Lend_pos = ranges_info_array(2, irep, proc_receive)
                     start_point = ranges_info_array(3, irep, proc_receive)
                     end_point = ranges_info_array(4, irep, proc_receive)

                     IF (.NOT. my_alpha_beta_case) THEN
                        local_j_aL(Lstart_pos:Lend_pos, 1:my_B_size, 1:my_block_size) = &
                           BI_C_rec(start_point:end_point, 1:my_B_size, 1:my_block_size)
                     ELSE
                        local_j_aL(Lstart_pos:Lend_pos, 1:my_B_size_beta, 1:my_block_size) = &
                           BI_C_rec(start_point:end_point, 1:my_B_size_beta, 1:my_block_size)
                     END IF

                  END DO

               END IF

               DEALLOCATE (BI_C_rec)

            END DO
            CALL timestop(handle3)

            ! loop over the block elements
            DO iiB = 1, my_block_size
               DO jjB = 1, my_block_size
                  CALL timeset(routineN//"_expansion", handle3)
                  ! calculate the integrals (ia|jb) strating from my local data ...
                  local_ab = 0.0_dp
                  IF ((my_alpha_beta_case) .AND. (calc_forces)) THEN
                     local_ba = 0.0_dp
                  ENDIF
                  t_start = m_walltime()
                  CALL dgemm('T', 'N', my_B_size, my_B_size_beta, dimen_RI, 1.0_dp, &
                             local_i_aL(:, :, iiB), dimen_RI, local_j_aL(:, :, jjB), dimen_RI, &
                             0.0_dp, local_ab(my_B_virtual_start:my_B_virtual_end, 1:my_B_size_beta), my_B_size)
                  t_end = m_walltime()
                  actual_flop_rate = 2.0_dp*my_B_size*my_B_size_beta*REAL(dimen_RI, KIND=dp)/(MAX(0.01_dp, t_end-t_start))
                  my_flop_rate = my_flop_rate+actual_flop_rate
                  my_num_dgemm_call = my_num_dgemm_call+1
                  ! Additional integrals only for alpha_beta case and forces
                  IF ((my_alpha_beta_case) .AND. (calc_forces)) THEN
                     t_start = m_walltime()
                     CALL dgemm('T', 'N', my_B_size_beta, my_B_size, dimen_RI, 1.0_dp, &
                                local_j_aL(:, :, iiB), dimen_RI, local_i_aL(:, :, jjB), dimen_RI, &
                                0.0_dp, local_ba(my_B_virtual_start_beta:my_B_virtual_end_beta, 1:my_B_size), my_B_size_beta)
                     t_end = m_walltime()
                     actual_flop_rate = 2.0_dp*my_B_size*my_B_size_beta*REAL(dimen_RI, KIND=dp)/(MAX(0.01_dp, t_end-t_start))
                     my_flop_rate = my_flop_rate+actual_flop_rate
                     my_num_dgemm_call = my_num_dgemm_call+1
                  ENDIF
                  ! ... and from the other of my subgroup
                  DO proc_shift = 1, para_env_sub%num_pe-1
                     proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
                     proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)

                     rec_B_size = sizes_B_virtual(proc_receive)
                     rec_B_virtual_end = ends_B_virtual(proc_receive)
                     rec_B_virtual_start = starts_B_virtual(proc_receive)

                     ALLOCATE (external_i_aL(dimen_RI, rec_B_size))
                     external_i_aL = 0.0_dp

                     CALL mp_sendrecv(local_i_aL(:, :, iiB), proc_send, &
                                      external_i_aL, proc_receive, &
                                      para_env_sub%group)

                     ! local_ab(rec_B_virtual_start:rec_B_virtual_end,1:my_B_size)=MATMUL(TRANSPOSE(external_i_aL),local_j_aL)
                     t_start = m_walltime()
                     CALL dgemm('T', 'N', rec_B_size, my_B_size_beta, dimen_RI, 1.0_dp, &
                                external_i_aL, dimen_RI, local_j_aL(:, :, jjB), dimen_RI, &
                                0.0_dp, local_ab(rec_B_virtual_start:rec_B_virtual_end, 1:my_B_size_beta), rec_B_size)

                     t_end = m_walltime()
                     actual_flop_rate = 2.0_dp*rec_B_size*my_B_size_beta*REAL(dimen_RI, KIND=dp)/(MAX(0.01_dp, t_end-t_start))
                     my_flop_rate = my_flop_rate+actual_flop_rate
                     my_num_dgemm_call = my_num_dgemm_call+1

                     DEALLOCATE (external_i_aL)
                     ! Additional integrals only for alpha_beta case and forces
                     IF ((my_alpha_beta_case) .AND. (calc_forces)) THEN

                        rec_B_size = sizes_B_virtual_beta(proc_receive)
                        rec_B_virtual_end = ends_B_virtual_beta(proc_receive)
                        rec_B_virtual_start = starts_B_virtual_beta(proc_receive)

                        ALLOCATE (external_i_aL(dimen_RI, rec_B_size))
                        external_i_aL = 0.0_dp

                        CALL mp_sendrecv(local_j_aL(:, :, jjB), proc_send, &
                                         external_i_aL, proc_receive, &
                                         para_env_sub%group)

                        t_start = m_walltime()
                        CALL dgemm('T', 'N', rec_B_size, my_B_size, dimen_RI, 1.0_dp, &
                                   external_i_aL, dimen_RI, local_i_aL(:, :, iiB), dimen_RI, &
                                   0.0_dp, local_ba(rec_B_virtual_start:rec_B_virtual_end, 1:my_B_size), rec_B_size)
                        t_end = m_walltime()
                        actual_flop_rate = 2.0_dp*rec_B_size*my_B_size*REAL(dimen_RI, KIND=dp)/(MAX(0.01_dp, t_end-t_start))
                        my_flop_rate = my_flop_rate+actual_flop_rate
                        my_num_dgemm_call = my_num_dgemm_call+1

                        DEALLOCATE (external_i_aL)
                     ENDIF

                  END DO
                  CALL timestop(handle3)

                  !sample peak memory
                  CALL m_memory()

                  CALL timeset(routineN//"_ener", handle3)
                  ! calculate coulomb only MP2
                  sym_fac = 2.0_dp
                  IF (my_i == my_j) sym_fac = 1.0_dp
                  IF (.NOT. my_alpha_beta_case) THEN
                     ! IF(my_open_shell_SS) sym_fac=sym_fac/2.0_dp
                     DO b = 1, my_B_size
                        b_global = b+my_B_virtual_start-1
                        DO a = 1, virtual
                           Emp2_Cou = Emp2_Cou-sym_fac*2.0_dp*local_ab(a, b)**2/ &
                                      (Eigenval(homo+a)+Eigenval(homo+b_global)-Eigenval(my_i+iiB-1)-Eigenval(my_j+jjB-1))
                        END DO
                     END DO
                  ELSE
                     DO b = 1, my_B_size_beta
                        b_global = b+my_B_virtual_start_beta-1
                        DO a = 1, virtual
                           Emp2_Cou = Emp2_Cou-local_ab(a, b)**2/ &
                                      (Eigenval(homo+a)+Eigenval_beta(homo_beta+b_global)- &
                                       Eigenval(my_i+iiB-1)-Eigenval_beta(my_j+jjB-1))
                        END DO
                     END DO
                  END IF

                  IF (calc_ex) THEN
                     ! contract integrals with orbital energies for exchange MP2 energy
                     ! starting with local ...
                     ! IF(my_open_shell_SS) sym_fac=sym_fac*2.0_dp
                     IF (calc_forces .AND. (.NOT. my_alpha_beta_case)) t_ab = 0.0_dp
                     DO b = 1, my_B_size
                        b_global = b+my_B_virtual_start-1
                        DO a = 1, my_B_size
                           a_global = a+my_B_virtual_start-1
                           Emp2_Ex = Emp2_Ex+sym_fac*local_ab(a_global, b)*local_ab(b_global, a)/ &
                                     (Eigenval(homo+a_global)+Eigenval(homo+b_global)-Eigenval(my_i+iiB-1)-Eigenval(my_j+jjB-1))
                           IF (calc_forces .AND. (.NOT. my_alpha_beta_case)) &
                              t_ab(a_global, b) = -(amp_fac*local_ab(a_global, b)-local_ab(b_global, a))/ &
                                                  (Eigenval(homo+a_global)+Eigenval(homo+b_global)- &
                                                   Eigenval(my_i+iiB-1)-Eigenval(my_j+jjB-1))
                        END DO
                     END DO
                     ! ... and then with external data
                     DO proc_shift = 1, para_env_sub%num_pe-1
                        proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
                        proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)

                        rec_B_size = sizes_B_virtual(proc_receive)
                        rec_B_virtual_end = ends_B_virtual(proc_receive)
                        rec_B_virtual_start = starts_B_virtual(proc_receive)

                        send_B_size = sizes_B_virtual(proc_send)
                        send_B_virtual_end = ends_B_virtual(proc_send)
                        send_B_virtual_start = starts_B_virtual(proc_send)

                        ALLOCATE (external_ab(my_B_size, rec_B_size))
                        external_ab = 0.0_dp

                        CALL mp_sendrecv(local_ab(send_B_virtual_start:send_B_virtual_end, 1:my_B_size), proc_send, &
                                         external_ab(1:my_B_size, 1:rec_B_size), proc_receive, &
                                         para_env_sub%group)

                        DO b = 1, my_B_size
                           b_global = b+my_B_virtual_start-1
                           DO a = 1, rec_B_size
                              a_global = a+rec_B_virtual_start-1
                              Emp2_Ex = Emp2_Ex+sym_fac*local_ab(a_global, b)*external_ab(b, a)/ &
                                        (Eigenval(homo+a_global)+Eigenval(homo+b_global)-Eigenval(my_i+iiB-1)-Eigenval(my_j+jjB-1))
                              IF (calc_forces .AND. (.NOT. my_alpha_beta_case)) &
                                 t_ab(a_global, b) = -(amp_fac*local_ab(a_global, b)-external_ab(b, a))/ &
                                                     (Eigenval(homo+a_global)+Eigenval(homo+b_global)- &
                                                      Eigenval(my_i+iiB-1)-Eigenval(my_j+jjB-1))
                           END DO
                        END DO

                        DEALLOCATE (external_ab)
                     END DO
                  END IF
                  CALL timestop(handle3)

                  IF (calc_forces) THEN
                     ! update P_ab, Gamma_P_ia
                     IF (.NOT. my_alpha_beta_case) THEN
                        Y_i_aP = 0.0_dp
                        Y_j_aP = 0.0_dp
                        CALL mp2_update_P_gamma(mp2_env, para_env_sub, ends_B_virtual, sizes_B_virtual, &
                                                starts_B_virtual, Eigenval, homo, dimen_RI, iiB, jjB, my_B_size, &
                                                my_B_virtual_end, my_B_virtual_start, my_i, my_j, virtual, &
                                                sub_proc_map, local_ab, t_ab, local_i_aL, local_j_aL, &
                                                my_open_shell_ss, my_alpha_alpha_case, my_beta_beta_case, Y_i_aP, Y_j_aP)
                     ELSE
                        Y_i_aP = 0.0_dp
                        Y_j_aP_beta = 0.0_dp
                        CALL mp2_update_P_gamma(mp2_env, para_env_sub, ends_B_virtual, sizes_B_virtual, &
                                                starts_B_virtual, Eigenval, homo, dimen_RI, iiB, jjB, my_B_size, &
                                                my_B_virtual_end, my_B_virtual_start, my_i, my_j, virtual, sub_proc_map, &
                                                local_ab, t_ab, local_i_aL, local_j_aL, my_open_shell_ss, my_alpha_alpha_case, &
                                                my_beta_beta_case, Y_i_aP, Y_j_aP_beta, Eigenval_beta, homo_beta, my_B_size_beta, &
                                                ends_B_virtual_beta, sizes_B_virtual_beta, starts_B_virtual_beta, &
                                                my_B_virtual_start_beta, my_B_virtual_end_beta, virtual_beta, local_ba)
                     ENDIF

                  END IF

               END DO ! jjB
            END DO ! iiB

         ELSE
            ! No work to do and we know that we have to receive nothing, but send something
            ! send data to other proc
            DO proc_shift = 1, para_env_exchange%num_pe-1
               proc_send = proc_map(para_env_exchange%mepos+proc_shift)
               proc_receive = proc_map(para_env_exchange%mepos-proc_shift)

               send_ij_index = num_ij_pairs(proc_send)

               IF (ij_index <= send_ij_index) THEN
                  ! something to send
                  ! ij_counter_send=(ij_index-MIN(1,proc_send))*ngroup+proc_send
                  ij_counter_send = (ij_index-MIN(1, integ_group_pos2color_sub(proc_send)))*ngroup+ &
                                    integ_group_pos2color_sub(proc_send)
                  send_i = ij_map(ij_counter_send, 1)
                  send_j = ij_map(ij_counter_send, 2)
                  send_block_size = ij_map(ij_counter_send, 3)

                  ! occupied i
                  CALL mp_sendrecv(BIb_C(1:my_group_L_size, 1:my_B_size, send_i:send_i+send_block_size-1), proc_send, &
                                   null_mat_rec, proc_receive, &
                                   para_env_exchange%group)
                  ! occupied j
                  IF (.NOT. my_alpha_beta_case) THEN
                     CALL mp_sendrecv(BIb_C(1:my_group_L_size, 1:my_B_size, send_j:send_j+send_block_size-1), proc_send, &
                                      null_mat_rec, proc_receive, &
                                      para_env_exchange%group)
                  ELSE
                     CALL mp_sendrecv(BIb_C_beta(1:my_group_L_size, 1:my_B_size_beta, send_j:send_j+send_block_size-1), proc_send, &
                                      null_mat_rec, proc_receive, &
                                      para_env_exchange%group)
                  END IF

               ELSE
                  ! nothing to send
                  ! occupied i
                  CALL mp_sendrecv(null_mat_send, proc_send, &
                                   null_mat_rec, proc_receive, &
                                   para_env_exchange%group)
                  ! occupied j
                  CALL mp_sendrecv(null_mat_send, proc_send, &
                                   null_mat_rec, proc_receive, &
                                   para_env_exchange%group)

               END IF
            END DO
         END IF

         ! redistribute gamma
         IF (calc_forces) THEN
            ! Closed shell, alpha-alpha or beta-beta case
            IF ((.NOT. my_alpha_beta_case)) THEN
               CALL mp2_redistribute_gamma(mp2_env, ij_index, my_B_size, &
                                           my_block_size, my_group_L_size, my_i, my_ij_pairs, my_j, ngroup, &
                                           num_integ_group, integ_group_pos2color_sub, num_ij_pairs, proc_map, &
                                           ij_map, ranges_info_array, Y_i_aP, Y_j_aP, para_env_exchange, &
                                           null_mat_rec, null_mat_send, sizes_array, my_alpha_alpha_case, &
                                           my_beta_beta_case, my_alpha_beta_case, my_open_shell_ss)
            ELSE
               ! Alpha-beta case
               CALL mp2_redistribute_gamma(mp2_env, ij_index, my_B_size, &
                                           my_block_size, my_group_L_size, my_i, my_ij_pairs, my_j, ngroup, &
                                           num_integ_group, integ_group_pos2color_sub, num_ij_pairs, proc_map, &
                                           ij_map, ranges_info_array, Y_i_aP, Y_j_aP_beta, para_env_exchange, &
                                           null_mat_rec, null_mat_send, sizes_array, my_alpha_alpha_case, my_beta_beta_case, &
                                           my_alpha_beta_case, my_open_shell_ss, my_B_size_beta)
            ENDIF
         END IF

      END DO
      CALL timestop(handle2)

      DEALLOCATE (local_i_aL)
      DEALLOCATE (local_j_aL)
      DEALLOCATE (ij_map)
      DEALLOCATE (num_ij_pairs)

      IF (calc_forces) THEN
         DEALLOCATE (Y_i_aP)
         IF (.NOT. my_alpha_beta_case) THEN
            DEALLOCATE (Y_j_aP)
         ELSE
            DEALLOCATE (Y_j_aP_beta)
         ENDIF
         IF (ALLOCATED(t_ab)) THEN
            DEALLOCATE (t_ab)
         ENDIF
         ! Deallocate additional integrals: alpha_beta case with forces
         IF (ALLOCATED(local_ba)) THEN
            DEALLOCATE (local_ba)
         ENDIF

         ! here we check if there are almost degenerate ij
         ! pairs and we update P_ij with these contribution.
         ! If all pairs are degenerate with each other this step will scale O(N^6),
         ! if the number of degenerate pairs scales linearly with the system size
         ! this step will scale O(N^5).
         ! Start counting the number of almost degenerate ij pairs according
         ! to eps_canonical
         IF (.NOT. my_alpha_beta_case) THEN
            CALL quasi_degenerate_P_ij( &
               mp2_env, Eigenval, homo, virtual, my_open_shell_ss, &
               my_beta_beta_case, my_alpha_beta_case, Bib_C, unit_nr, dimen_RI, &
               my_B_size, ngroup, num_integ_group, my_group_L_size, &
               color_sub, ranges_info_array, para_env_exchange, para_env_sub, proc_map, &
               my_B_virtual_start, my_B_virtual_end, sizes_array, ends_B_virtual, sizes_B_virtual, &
               starts_B_virtual, sub_proc_map, integ_group_pos2color_sub, local_ab)
         ELSE
            CALL quasi_degenerate_P_ij( &
               mp2_env, Eigenval, homo, virtual, my_open_shell_ss, &
               my_beta_beta_case, my_alpha_beta_case, Bib_C, unit_nr, dimen_RI, &
               my_B_size, ngroup, num_integ_group, my_group_L_size, &
               color_sub, ranges_info_array, para_env_exchange, para_env_sub, proc_map, &
               my_B_virtual_start, my_B_virtual_end, sizes_array, ends_B_virtual, sizes_B_virtual, &
               starts_B_virtual, sub_proc_map, integ_group_pos2color_sub, local_ab, BIb_C_beta, my_B_size_beta, &
               ends_B_virtual_beta, sizes_B_virtual_beta, starts_B_virtual_beta, my_B_virtual_start_beta, &
               virtual_beta, homo_beta, Eigenval_beta, my_B_virtual_end_beta)
         ENDIF

      END IF

      DEALLOCATE (integ_group_pos2color_sub)
      DEALLOCATE (local_ab)

      CALL mp_sum(Emp2_Cou, para_env%group)
      CALL mp_sum(Emp2_Ex, para_env%group)

      IF (calc_forces) THEN
         ! sum P_ab
         IF (.NOT. my_open_shell_ss) THEN
            mp2_env%ri_grad%P_ab(:, :) = mp2_env%ri_grad%P_ab(:, :)*amp_fac
            IF (my_alpha_beta_case) mp2_env%ri_grad%P_ab_beta(:, :) = &
               mp2_env%ri_grad%P_ab_beta(:, :)*amp_fac
            sub_P_color = para_env_sub%mepos
            CALL mp_comm_split_direct(para_env%group, comm_P, sub_P_color)
            NULLIFY (para_env_P)
            CALL cp_para_env_create(para_env_P, comm_P)
            CALL mp_sum(mp2_env%ri_grad%P_ab, para_env_P%group)
            IF (my_alpha_beta_case) CALL mp_sum(mp2_env%ri_grad%P_ab_beta, para_env_P%group)
            ! release para_env_P
            CALL cp_para_env_release(para_env_P)
         ENDIF

         ! sum P_ij (later)
         ! mp2_env%ri_grad%P_ij=mp2_env%ri_grad%P_ij*2.0_dp
         ! CALL mp_sum(mp2_env%ri_grad%P_ij,para_env%group)
         ! WRITE(*,*) mp2_env%ri_mp2%eps_canonical

         ! XXXXXXXXXXXXXX
         ! Write P_ab
         ! ALLOCATE(external_ab(virtual,virtual))
         ! external_ab=0.0_dp
         ! external_ab(my_B_virtual_start:my_B_virtual_end,1:virtual)=mp2_env%ri_grad%P_ab
         ! CALL mp_sum(external_ab,para_env_sub%group)
         ! IF(para_env%mepos==0) CALL write_array(external_ab)
         ! DEALLOCATE(external_ab)
         ! ! Write P_ij
         ! IF(para_env%mepos==0) CALL write_array(mp2_env%ri_grad%P_ij)
         ! XXXXXXXXXXXXXX

         ! recover original information (before replication)
         DEALLOCATE (sizes_array)
         iiB = SIZE(sizes_array_orig)
         ALLOCATE (sizes_array(0:iiB-1))
         sizes_array(:) = sizes_array_orig
         DEALLOCATE (sizes_array_orig)

         ! make a copy of the original integrals (ia|Q)
         my_group_L_size = my_group_L_size_orig
         ALLOCATE (B_ia_Q(homo, my_B_size, my_group_L_size))
         B_ia_Q = 0.0_dp
         DO jjB = 1, homo
            DO iiB = 1, my_B_size
               B_ia_Q(jjB, iiB, 1:my_group_L_size) = BIb_C(1:my_group_L_size, iiB, jjB)
            END DO
         END DO
         DEALLOCATE (BIb_C)
         IF (my_alpha_beta_case) THEN
            ALLOCATE (B_ia_Q_beta(homo_beta, my_B_size_beta, my_group_L_size))
            B_ia_Q_beta = 0.0_dp
            DO jjB = 1, homo_beta
               DO iiB = 1, my_B_size_beta
                  B_ia_Q_beta(jjB, iiB, 1:my_group_L_size) = &
                     BIb_C_beta(1:my_group_L_size, iiB, jjB)
               END DO
            END DO
            DEALLOCATE (BIb_C_beta)
         ENDIF

         ! sum Gamma and dereplicate
         ALLOCATE (BIb_C(homo, my_B_size, my_group_L_size))
         IF (my_alpha_beta_case) ALLOCATE (BIb_C_beta(homo_beta, my_B_size_beta, my_group_L_size))
         DO proc_shift = 1, para_env_rep%num_pe-1
            ! invert order
            proc_send = proc_map_rep(para_env_rep%mepos-proc_shift)
            proc_receive = proc_map_rep(para_env_rep%mepos+proc_shift)

            start_point = ranges_info_array(3, proc_shift, para_env_exchange%mepos)
            end_point = ranges_info_array(4, proc_shift, para_env_exchange%mepos)

            BIb_C = 0.0_dp
            ! Closed shell, alpha-alpha, and alpha-alpha part in alpha-beta case
            IF (my_alpha_alpha_case .OR. (.NOT. my_open_shell_ss)) THEN
               CALL mp_sendrecv(mp2_env%ri_grad%Gamma_P_ia(1:homo, 1:my_B_size, start_point:end_point), &
                                proc_send, BIb_C, proc_receive, para_env_rep%group)
               mp2_env%ri_grad%Gamma_P_ia(1:homo, 1:my_B_size, 1:my_group_L_size) = &
                  mp2_env%ri_grad%Gamma_P_ia(1:homo, 1:my_B_size, 1:my_group_L_size)+BIb_C
            ENDIF
            ! Beta-beta
            IF (my_beta_beta_case) THEN
               CALL mp_sendrecv(mp2_env%ri_grad%Gamma_P_ia_beta(1:homo, 1:my_B_size, start_point:end_point), &
                                proc_send, BIb_C, proc_receive, para_env_rep%group)
               mp2_env%ri_grad%Gamma_P_ia_beta(1:homo, 1:my_B_size, 1:my_group_L_size) = &
                  mp2_env%ri_grad%Gamma_P_ia_beta(1:homo, 1:my_B_size, 1:my_group_L_size)+BIb_C
            ENDIF
            IF (my_alpha_beta_case) THEN ! Beta-beta part of alpha-beta case
               CALL mp_sendrecv(mp2_env%ri_grad%Gamma_P_ia_beta(1:homo_beta, 1:my_B_size_beta, start_point:end_point), &
                                proc_send, BIb_C_beta, proc_receive, para_env_rep%group)
               mp2_env%ri_grad%Gamma_P_ia_beta(1:homo_beta, 1:my_B_size_beta, 1:my_group_L_size) = &
                  mp2_env%ri_grad%Gamma_P_ia_beta(1:homo_beta, 1:my_B_size_beta, 1:my_group_L_size)+BIb_C_beta
            ENDIF

         END DO
         IF (.NOT. my_open_shell_ss) THEN
            BIb_C(:, :, :) = mp2_env%ri_grad%Gamma_P_ia(1:homo, 1:my_B_size, 1:my_group_L_size)
            DEALLOCATE (mp2_env%ri_grad%Gamma_P_ia)
            ALLOCATE (mp2_env%ri_grad%Gamma_P_ia(homo, my_B_size, my_group_L_size))
            mp2_env%ri_grad%Gamma_P_ia(:, :, :) = BIb_C
            DEALLOCATE (BIb_C)
            IF (my_alpha_beta_case) THEN
               BIb_C_beta(:, :, :) = mp2_env%ri_grad%Gamma_P_ia_beta(1:homo_beta, 1:my_B_size_beta, 1:my_group_L_size)
               DEALLOCATE (mp2_env%ri_grad%Gamma_P_ia_beta)
               ALLOCATE (mp2_env%ri_grad%Gamma_P_ia_beta(homo_beta, my_B_size_beta, my_group_L_size))
               mp2_env%ri_grad%Gamma_P_ia_beta(:, :, :) = BIb_C_beta
               DEALLOCATE (BIb_C_beta)
            ENDIF
         ENDIF
         ! For open shell systems, we need to pass Bib_C through the subroutine in alpha-alpha and beta-beta case.
         ! Only for forces, as IF suggests! Here we deallocate it, but restore after complete_gamma
         IF (my_open_shell_ss) THEN
            DEALLOCATE (BIb_C)
         ENDIF

         IF (.NOT. my_open_shell_ss) THEN
            CALL complete_gamma(mp2_env, B_ia_Q, dimen_RI, homo, virtual, para_env, para_env_sub, ngroup, &
                                my_group_L_size, my_group_L_start, my_group_L_end, &
                                my_B_size, my_B_virtual_start, &
                                ends_array, ends_B_virtual, sizes_array, sizes_B_virtual, starts_array, starts_B_virtual, &
                                sub_proc_map, .TRUE.)
            IF (my_alpha_beta_case) THEN
               CALL complete_gamma(mp2_env, B_ia_Q_beta, dimen_RI, homo_beta, virtual_beta, para_env, para_env_sub, &
                                   ngroup, my_group_L_size, my_group_L_start, my_group_L_end, &
                                   my_B_size_beta, my_B_virtual_start_beta, &
                                   ends_array, ends_B_virtual_beta, sizes_array, sizes_B_virtual_beta, starts_array, &
                                   starts_B_virtual_beta, sub_proc_map, .FALSE.)
            ENDIF
         ENDIF
         ! Here we restore BIb_C
         IF (my_open_shell_ss) THEN
            ALLOCATE (BIb_C(my_group_L_size, my_B_size, homo))
            BIb_C = 0.0_dp
            ! copy the integrals (ia|Q) back
            DO jjB = 1, homo
               DO iiB = 1, my_B_size
                  BIb_C(1:my_group_L_size, iiB, jjB) = &
                     B_ia_Q(jjB, iiB, 1:my_group_L_size)
               END DO
            END DO
         ENDIF

         ! ALLOCATE(Y_i_aP(homo,virtual,dimen_RI))
         ! Y_i_aP=0.0_dp
         ! Y_i_aP(1:homo,my_B_virtual_start:my_B_virtual_end,my_group_L_start:my_group_L_end)=mp2_env%ri_grad%Gamma_P_ia
         ! CALL mp_sum(Y_i_aP,para_env%group)
         ! IF(para_env%mepos==0) THEN
         !   DO iiB=1, dimen_RI
         !     WRITE(*,*) iiB
         !     CALL write_array(Y_i_aP(1:homo,1:virtual,iiB))
         !     WRITE(*,*)
         !   END DO
         ! END IF
         ! DEALLOCATE(Y_i_aP)

      END IF

      Emp2 = Emp2_Cou+Emp2_EX

      DEALLOCATE (proc_map)
      DEALLOCATE (sub_proc_map)
      DEALLOCATE (proc_map_rep)
      DEALLOCATE (ranges_info_array)

      IF (.NOT. my_open_shell_SS) THEN
         ! keep the array for the next calculations
         IF (ALLOCATED(BIb_C)) DEALLOCATE (BIb_C)
         DEALLOCATE (sizes_array)
         DEALLOCATE (starts_array)
         DEALLOCATE (ends_array)
         DEALLOCATE (starts_B_virtual)
         DEALLOCATE (ends_B_virtual)
         DEALLOCATE (sizes_B_virtual)
         IF (my_alpha_beta_case) THEN
            DEALLOCATE (starts_B_virtual_beta)
            DEALLOCATE (ends_B_virtual_beta)
            DEALLOCATE (sizes_B_virtual_beta)
         END IF
      END IF

      CALL cp_para_env_release(para_env_exchange)
      CALL cp_para_env_release(para_env_rep)

      my_flop_rate = my_flop_rate/REAL(MAX(my_num_dgemm_call, 1), KIND=dp)/1.0E9_dp
      CALL mp_sum(my_flop_rate, para_env%group)
      my_flop_rate = my_flop_rate/para_env%num_pe
      IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.2)") &
         "PERFORMANCE| DGEMM flop rate (Gflops / MPI rank):", my_flop_rate

      CALL timestop(handle)

   END SUBROUTINE mp2_ri_gpw_compute_en

! **************************************************************************************************
!> \brief ...
!> \param BIb_C ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param para_env_exchange ...
!> \param para_env_rep ...
!> \param homo ...
!> \param proc_map_rep ...
!> \param sizes_array ...
!> \param my_B_size ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param my_new_group_L_size ...
!> \param new_sizes_array ...
!> \param ranges_info_array ...
! **************************************************************************************************
   SUBROUTINE replicate_iaK_2intgroup(BIb_C, para_env, para_env_sub, para_env_exchange, para_env_rep, &
                                      homo, proc_map_rep, &
                                      sizes_array, &
                                      my_B_size, &
                                      my_group_L_size, my_group_L_start, my_group_L_end, &
                                      my_new_group_L_size, new_sizes_array, ranges_info_array)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub, &
                                                            para_env_exchange, para_env_rep
      INTEGER                                            :: homo
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: proc_map_rep, sizes_array
      INTEGER                                            :: my_B_size, my_group_L_size, &
                                                            my_group_L_start, my_group_L_end, &
                                                            my_new_group_L_size
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: new_sizes_array
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: ranges_info_array

      CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_iaK_2intgroup', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: comm_rep, end_point, handle, i, &
                                                            max_L_size, proc_receive, proc_send, &
                                                            proc_shift, start_point, sub_sub_color
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rep_ends_array, rep_sizes_array, &
                                                            rep_starts_array
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C_copy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: BIb_C_gather

      CALL timeset(routineN, handle)

      ! create the replication group
      sub_sub_color = para_env_sub%mepos*para_env_exchange%num_pe+para_env_exchange%mepos
      CALL mp_comm_split_direct(para_env%group, comm_rep, sub_sub_color)
      NULLIFY (para_env_rep)
      CALL cp_para_env_create(para_env_rep, comm_rep)

      ! crate the proc maps
      ALLOCATE (proc_map_rep(-para_env_rep%num_pe:2*para_env_rep%num_pe-1))
      DO i = 0, para_env_rep%num_pe-1
         proc_map_rep(i) = i
         proc_map_rep(-i-1) = para_env_rep%num_pe-i-1
         proc_map_rep(para_env_rep%num_pe+i) = i
      END DO

      ! create the new limits for K according to the size
      ! of the integral group
      ALLOCATE (new_sizes_array(0:para_env_exchange%num_pe-1))
      new_sizes_array = 0
      ALLOCATE (ranges_info_array(4, 0:para_env_rep%num_pe-1, 0:para_env_exchange%num_pe-1))
      ranges_info_array = 0

      ! info array for replication
      ALLOCATE (rep_ends_array(0:para_env_rep%num_pe-1))
      rep_ends_array = 0
      ALLOCATE (rep_starts_array(0:para_env_rep%num_pe-1))
      rep_starts_array = 0
      ALLOCATE (rep_sizes_array(0:para_env_rep%num_pe-1))
      rep_sizes_array = 0

      rep_sizes_array(para_env_rep%mepos) = my_group_L_size
      rep_starts_array(para_env_rep%mepos) = my_group_L_start
      rep_ends_array(para_env_rep%mepos) = my_group_L_end

      CALL mp_sum(rep_sizes_array, para_env_rep%group)
      CALL mp_sum(rep_starts_array, para_env_rep%group)
      CALL mp_sum(rep_ends_array, para_env_rep%group)

      ! calculate my_new_group_L_size according to sizes_array
      my_new_group_L_size = my_group_L_size
      ranges_info_array(1, 0, para_env_exchange%mepos) = my_group_L_start
      ranges_info_array(2, 0, para_env_exchange%mepos) = my_group_L_end
      ranges_info_array(3, 0, para_env_exchange%mepos) = 1
      ranges_info_array(4, 0, para_env_exchange%mepos) = my_group_L_size

      DO proc_shift = 1, para_env_rep%num_pe-1
         proc_send = proc_map_rep(para_env_rep%mepos+proc_shift)
         proc_receive = proc_map_rep(para_env_rep%mepos-proc_shift)

         my_new_group_L_size = my_new_group_L_size+rep_sizes_array(proc_receive)

         ranges_info_array(1, proc_shift, para_env_exchange%mepos) = rep_starts_array(proc_receive)
         ranges_info_array(2, proc_shift, para_env_exchange%mepos) = rep_ends_array(proc_receive)
         ranges_info_array(3, proc_shift, para_env_exchange%mepos) = ranges_info_array(4, proc_shift-1, para_env_exchange%mepos)+1
         ranges_info_array(4, proc_shift, para_env_exchange%mepos) = my_new_group_L_size

      END DO
      new_sizes_array(para_env_exchange%mepos) = my_new_group_L_size

      CALL mp_sum(new_sizes_array, para_env_exchange%group)
      CALL mp_sum(ranges_info_array, para_env_exchange%group)

      IF (.FALSE.) THEN
         ! replication scheme using mp_sendrecv
         ALLOCATE (BIb_C_copy(my_group_L_size, my_B_size, homo))
         BIb_C_copy(:, :, :) = BIb_C

         DEALLOCATE (BIb_C)

         ALLOCATE (BIb_C(my_new_group_L_size, my_B_size, homo))
         BIb_C = 0.0_dp

         start_point = ranges_info_array(3, 0, para_env_exchange%mepos)
         end_point = ranges_info_array(4, 0, para_env_exchange%mepos)

         BIb_C(start_point:end_point, 1:my_B_size, 1:homo) = BIb_C_copy

         DO proc_shift = 1, para_env_rep%num_pe-1
            proc_send = proc_map_rep(para_env_rep%mepos+proc_shift)
            proc_receive = proc_map_rep(para_env_rep%mepos-proc_shift)

            start_point = ranges_info_array(3, proc_shift, para_env_exchange%mepos)
            end_point = ranges_info_array(4, proc_shift, para_env_exchange%mepos)

            CALL mp_sendrecv(BIb_C_copy(1:my_group_L_size, 1:my_B_size, 1:homo), proc_send, &
                             BIb_C(start_point:end_point, 1:my_B_size, 1:homo), proc_receive, &
                             para_env_rep%group)

         END DO

         DEALLOCATE (BIb_C_copy)

      ELSE
         ! replication scheme using mp_allgather
         ! get the max L size of the
         max_L_size = MAXVAL(sizes_array)

         ALLOCATE (BIb_C_copy(max_L_size, my_B_size, homo))
         BIb_C_copy = 0.0_dp
         BIb_C_copy(1:my_group_L_size, 1:my_B_size, 1:homo) = BIb_C

         DEALLOCATE (BIb_C)

         ALLOCATE (BIb_C_gather(max_L_size, my_B_size, homo, 0:para_env_rep%num_pe-1))
         BIb_C_gather = 0.0_dp

         CALL mp_allgather(BIb_C_copy, BIb_C_gather, para_env_rep%group)

         DEALLOCATE (BIb_C_copy)

         ALLOCATE (BIb_C(my_new_group_L_size, my_B_size, homo))
         BIb_C = 0.0_dp

         ! reorder data
         DO proc_shift = 0, para_env_rep%num_pe-1
            proc_send = proc_map_rep(para_env_rep%mepos+proc_shift)
            proc_receive = proc_map_rep(para_env_rep%mepos-proc_shift)

            start_point = ranges_info_array(3, proc_shift, para_env_exchange%mepos)
            end_point = ranges_info_array(4, proc_shift, para_env_exchange%mepos)

            BIb_C(start_point:end_point, 1:my_B_size, 1:homo) = &
               BIb_C_gather(1:end_point-start_point+1, 1:my_B_size, 1:homo, proc_receive)

         END DO

         DEALLOCATE (BIb_C_gather)

      END IF

      ! DEALLOCATE(proc_map_rep)
      DEALLOCATE (rep_sizes_array)
      DEALLOCATE (rep_starts_array)
      DEALLOCATE (rep_ends_array)

      ! CALL  cp_para_env_release(para_env_rep)

      CALL timestop(handle)

   END SUBROUTINE replicate_iaK_2intgroup

! **************************************************************************************************
!> \brief ...
!> \param mat ...
! **************************************************************************************************
   SUBROUTINE write_array(mat)
      REAL(KIND=dp), DIMENSION(:, :)                     :: mat

      INTEGER                                            :: iii, jjj

      WRITE (*, *)
      DO iii = 1, SIZE(mat, 1)
         DO jjj = 1, SIZE(mat, 2), 20
            WRITE (*, '(1000F10.5)') mat(iii, jjj:MIN(SIZE(mat, 2), jjj+19))
         END DO
         WRITE (*, *)
      END DO
      WRITE (*, *)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param local_ab ...
!> \param t_ab ...
!> \param mp2_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param my_B_size ...
!> \param block_size ...
!> \param my_B_size_beta ...
!> \param my_group_L_size ...
!> \param local_i_aL ...
!> \param local_j_aL ...
!> \param calc_forces ...
!> \param Y_i_aP ...
!> \param Y_j_aP ...
!> \param alpha_beta ...
!> \param beta_beta ...
!> \param local_ba ...
!> \param virtual_beta ...
! **************************************************************************************************
   SUBROUTINE mp2_ri_allocate(local_ab, t_ab, mp2_env, homo, virtual, dimen_RI, my_B_size, &
                              block_size, my_B_size_beta, my_group_L_size, &
                              local_i_aL, local_j_aL, calc_forces, &
                              Y_i_aP, Y_j_aP, alpha_beta, &
                              beta_beta, local_ba, virtual_beta)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: local_ab, t_ab
      TYPE(mp2_type), POINTER                            :: mp2_env
      INTEGER                                            :: homo, virtual, dimen_RI, my_B_size, &
                                                            block_size, my_B_size_beta, &
                                                            my_group_L_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: local_i_aL, local_j_aL
      LOGICAL                                            :: calc_forces
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Y_i_aP, Y_j_aP
      LOGICAL                                            :: alpha_beta, beta_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         OPTIONAL                                        :: local_ba
      INTEGER, OPTIONAL                                  :: virtual_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_allocate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ALLOCATE (local_i_aL(dimen_RI, my_B_size, block_size))
      ALLOCATE (local_j_aL(dimen_RI, my_B_size_beta, block_size))
      ALLOCATE (local_ab(virtual, my_B_size_beta))

      IF (calc_forces) THEN
         ALLOCATE (Y_i_aP(my_B_size, dimen_RI, block_size))
         Y_i_aP = 0.0_dp
         ! For  closed-shell, alpha-alpha and beta-beta my_B_size_beta=my_b_size
         ! Not for alpha-beta case: Y_j_aP_beta is sent and received as Y_j_aP
         ALLOCATE (Y_j_aP(my_B_size_beta, dimen_RI, block_size))
         Y_j_aP = 0.0_dp
         ! Closed shell or alpha-alpha case
         IF (.NOT. (beta_beta .OR. alpha_beta)) THEN
            ALLOCATE (mp2_env%ri_grad%P_ij(homo, homo))
            ALLOCATE (mp2_env%ri_grad%P_ab(my_B_size, virtual))
            mp2_env%ri_grad%P_ij = 0.0_dp
            mp2_env%ri_grad%P_ab = 0.0_dp
            ALLOCATE (mp2_env%ri_grad%Gamma_P_ia(homo, my_B_size, my_group_L_size))
            mp2_env%ri_grad%Gamma_P_ia = 0.0_dp
         ELSE
            IF (beta_beta) THEN
               ALLOCATE (mp2_env%ri_grad%P_ij_beta(homo, homo))
               ALLOCATE (mp2_env%ri_grad%P_ab_beta(my_B_size, virtual))
               mp2_env%ri_grad%P_ij_beta = 0.0_dp
               mp2_env%ri_grad%P_ab_beta = 0.0_dp
               ALLOCATE (mp2_env%ri_grad%Gamma_P_ia_beta(homo, my_B_size_beta, my_group_L_size))
               mp2_env%ri_grad%Gamma_P_ia_beta = 0.0_dp
            ENDIF
         ENDIF
         IF (.NOT. alpha_beta) THEN
            ! For non-alpha-beta case we need amplitudes
            ALLOCATE (t_ab(virtual, my_B_size_beta))
         ELSE
            ! We need more integrals
            ALLOCATE (local_ba(virtual_beta, my_B_size))
         ENDIF
      END IF
      !

      CALL timestop(handle)

   END SUBROUTINE mp2_ri_allocate

! **************************************************************************************************
!> \brief ...
!> \param my_alpha_beta_case ...
!> \param total_ij_pairs ...
!> \param homo ...
!> \param homo_beta ...
!> \param num_IJ_blocks ...
!> \param block_size ...
!> \param ngroup ...
!> \param ij_map ...
!> \param color_sub ...
!> \param my_ij_pairs ...
!> \param my_open_shell_SS ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE mp2_ri_communication(my_alpha_beta_case, total_ij_pairs, homo, homo_beta, num_IJ_blocks, &
                                   block_size, ngroup, ij_map, color_sub, my_ij_pairs, my_open_shell_SS, unit_nr)
      LOGICAL                                            :: my_alpha_beta_case
      INTEGER                                            :: total_ij_pairs, homo, homo_beta, &
                                                            num_IJ_blocks, block_size, ngroup
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ij_map
      INTEGER                                            :: color_sub, my_ij_pairs
      LOGICAL                                            :: my_open_shell_SS
      INTEGER                                            :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_communication', &
         routineP = moduleN//':'//routineN

      INTEGER :: assigned_blocks, first_I_block, first_J_block, handle, iiB, ij_block_counter, &
         ij_counter, jjB, last_i_block, last_J_block, num_block_per_group, total_ij_block, &
         total_ij_pairs_blocks
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ij_marker

! Calculate the maximum number of ij pairs that have to be computed
! among groups

      CALL timeset(routineN, handle)

      IF (.NOT. my_alpha_beta_case) THEN
         total_ij_pairs = homo*(1+homo)/2
         num_IJ_blocks = homo/block_size-1

         first_I_block = 1
         last_i_block = block_size*(num_IJ_blocks-1)

         first_J_block = block_size+1
         last_J_block = block_size*(num_IJ_blocks+1)

         ij_block_counter = 0
         DO iiB = first_I_block, last_i_block, block_size
            DO jjB = iiB+block_size, last_J_block, block_size
               ij_block_counter = ij_block_counter+1
            END DO
         END DO

         total_ij_block = ij_block_counter
         num_block_per_group = total_ij_block/ngroup
         assigned_blocks = num_block_per_group*ngroup

         total_ij_pairs_blocks = assigned_blocks+(total_ij_pairs-assigned_blocks*(block_size**2))

         ALLOCATE (ij_marker(homo, homo))
         ij_marker = 0
         ALLOCATE (ij_map(total_ij_pairs_blocks, 3))
         ij_map = 0
         ij_counter = 0
         my_ij_pairs = 0
         DO iiB = first_I_block, last_i_block, block_size
            DO jjB = iiB+block_size, last_J_block, block_size
               IF (ij_counter+1 > assigned_blocks) EXIT
               ij_counter = ij_counter+1
               ij_marker(iiB:iiB+block_size-1, jjB:jjB+block_size-1) = 1
               ij_map(ij_counter, 1) = iiB
               ij_map(ij_counter, 2) = jjB
               ij_map(ij_counter, 3) = block_size
               IF (MOD(ij_counter, ngroup) == color_sub) my_ij_pairs = my_ij_pairs+1
            END DO
         END DO
         DO iiB = 1, homo
            DO jjB = iiB, homo
               IF (ij_marker(iiB, jjB) == 0) THEN
                  ij_counter = ij_counter+1
                  ij_map(ij_counter, 1) = iiB
                  ij_map(ij_counter, 2) = jjB
                  ij_map(ij_counter, 3) = 1
                  IF (MOD(ij_counter, ngroup) == color_sub) my_ij_pairs = my_ij_pairs+1
               END IF
            END DO
         END DO
         DEALLOCATE (ij_marker)

         IF ((.NOT. my_open_shell_SS)) THEN
            IF (unit_nr > 0) THEN
               IF (block_size == 1) THEN
                  WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.1)") &
                     "RI_INFO| Percentage of ij pairs communicated with block size 1:", 100.0_dp
               ELSE
                  WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.1)") &
                     "RI_INFO| Percentage of ij pairs communicated with block size 1:", &
                     100.0_dp*REAL((total_ij_pairs-assigned_blocks*(block_size**2)), KIND=dp)/REAL(total_ij_pairs, KIND=dp)
               END IF
               CALL m_flush(unit_nr)
            END IF
         END IF

      ELSE
         ! alpha-beta case no index symmetry
         total_ij_pairs = homo*homo_beta
         ALLOCATE (ij_map(total_ij_pairs, 3))
         ij_map = 0
         ij_counter = 0
         my_ij_pairs = 0
         DO iiB = 1, homo
            DO jjB = 1, homo_beta
               ij_counter = ij_counter+1
               ij_map(ij_counter, 1) = iiB
               ij_map(ij_counter, 2) = jjB
               ij_map(ij_counter, 3) = 1
               IF (MOD(ij_counter, ngroup) == color_sub) my_ij_pairs = my_ij_pairs+1
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE mp2_ri_communication

! **************************************************************************************************
!> \brief ...
!> \param BIb_C ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param homo ...
!> \param color_sub ...
!> \param sizes_array ...
!> \param calc_forces ...
!> \param comm_exchange ...
!> \param integ_group_size ...
!> \param my_B_size ...
!> \param iiB ...
!> \param my_group_L_end ...
!> \param my_group_L_size ...
!> \param my_group_L_size_orig ...
!> \param my_group_L_start ...
!> \param my_new_group_L_size ...
!> \param sub_sub_color ...
!> \param integ_group_pos2color_sub ...
!> \param new_sizes_array ...
!> \param proc_map ...
!> \param proc_map_rep ...
!> \param sizes_array_orig ...
!> \param sub_proc_map ...
!> \param ranges_info_array ...
!> \param para_env_exchange ...
!> \param para_env_rep ...
!> \param num_integ_group ...
! **************************************************************************************************
   SUBROUTINE mp2_ri_create_group(BIb_C, para_env, para_env_sub, homo, color_sub, &
                                  sizes_array, calc_forces, &
                                  comm_exchange, integ_group_size, my_B_size, iiB, my_group_L_end, &
                                  my_group_L_size, my_group_L_size_orig, my_group_L_start, my_new_group_L_size, &
                                  sub_sub_color, integ_group_pos2color_sub, new_sizes_array, &
                                  proc_map, proc_map_rep, sizes_array_orig, &
                                  sub_proc_map, ranges_info_array, para_env_exchange, para_env_rep, num_integ_group)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER                                            :: homo, color_sub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_array
      LOGICAL, INTENT(IN)                                :: calc_forces
      INTEGER :: comm_exchange, integ_group_size, my_B_size, iiB, my_group_L_end, my_group_L_size, &
         my_group_L_size_orig, my_group_L_start, my_new_group_L_size, sub_sub_color
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: integ_group_pos2color_sub, &
                                                            new_sizes_array, proc_map, &
                                                            proc_map_rep, sizes_array_orig, &
                                                            sub_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: ranges_info_array
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange, para_env_rep
      INTEGER                                            :: num_integ_group

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_create_group', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i

!

      CALL timeset(routineN, handle)
      !
      sub_sub_color = para_env_sub%mepos*num_integ_group+color_sub/integ_group_size
      CALL mp_comm_split_direct(para_env%group, comm_exchange, sub_sub_color)
      NULLIFY (para_env_exchange)
      CALL cp_para_env_create(para_env_exchange, comm_exchange)

      ! create the proc maps
      ALLOCATE (proc_map(-para_env_exchange%num_pe:2*para_env_exchange%num_pe-1))
      DO i = 0, para_env_exchange%num_pe-1
         proc_map(i) = i
         proc_map(-i-1) = para_env_exchange%num_pe-i-1
         proc_map(para_env_exchange%num_pe+i) = i
      END DO

      ALLOCATE (sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe-1))
      DO i = 0, para_env_sub%num_pe-1
         sub_proc_map(i) = i
         sub_proc_map(-i-1) = para_env_sub%num_pe-i-1
         sub_proc_map(para_env_sub%num_pe+i) = i
      END DO

      CALL replicate_iaK_2intgroup(BIb_C, para_env, para_env_sub, para_env_exchange, para_env_rep, &
                                   homo, proc_map_rep, &
                                   sizes_array, &
                                   my_B_size, &
                                   my_group_L_size, my_group_L_start, my_group_L_end, &
                                   my_new_group_L_size, new_sizes_array, ranges_info_array)

      ALLOCATE (integ_group_pos2color_sub(0:para_env_exchange%num_pe-1))
      integ_group_pos2color_sub = 0
      integ_group_pos2color_sub(para_env_exchange%mepos) = color_sub
      CALL mp_sum(integ_group_pos2color_sub, para_env_exchange%group)

      IF (calc_forces) THEN
         iiB = SIZE(sizes_array)
         ALLOCATE (sizes_array_orig(0:iiB-1))
         sizes_array_orig(:) = sizes_array
      END IF

      my_group_L_size_orig = my_group_L_size
      my_group_L_size = my_new_group_L_size
      DEALLOCATE (sizes_array)

      ALLOCATE (sizes_array(0:integ_group_size-1))
      sizes_array(:) = new_sizes_array

      DEALLOCATE (new_sizes_array)
      !
      CALL timestop(handle)

   END SUBROUTINE mp2_ri_create_group

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param ends_array ...
!> \param ends_B_virtual ...
!> \param sizes_array ...
!> \param sizes_B_virtual ...
!> \param starts_array ...
!> \param starts_B_virtual ...
!> \param homo ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param color_sub ...
!> \param best_block_size ...
!> \param best_integ_group_size ...
!> \param block_size ...
!> \param integ_group_size ...
!> \param min_integ_group_size ...
!> \param my_B_size ...
!> \param my_B_virtual_end ...
!> \param my_B_virtual_start ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param ngroup ...
!> \param num_IJ_blocks ...
!> \param num_integ_group ...
!> \param pos_integ_group ...
!> \param virtual ...
!> \param my_alpha_beta_case ...
!> \param my_open_shell_SS ...
!> \param mem_for_aK ...
!> \param mem_for_comm ...
!> \param mem_for_iaK ...
!> \param mem_for_rep ...
!> \param mem_min ...
!> \param mem_per_group ...
!> \param mem_real ...
! **************************************************************************************************
   SUBROUTINE mp2_ri_get_sizes(mp2_env, para_env, para_env_sub, ends_array, ends_B_virtual, sizes_array, sizes_B_virtual, &
                               starts_array, starts_B_virtual, homo, dimen_RI, unit_nr, color_sub, &
                               best_block_size, best_integ_group_size, block_size, &
                               integ_group_size, min_integ_group_size, my_B_size, &
                               my_B_virtual_end, my_B_virtual_start, my_group_L_size, &
                               my_group_L_start, my_group_L_end, ngroup, num_IJ_blocks, num_integ_group, &
                               pos_integ_group, virtual, my_alpha_beta_case, &
                               my_open_shell_SS, mem_for_aK, mem_for_comm, &
                               mem_for_iaK, mem_for_rep, mem_min, mem_per_group, mem_real)
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ends_array, ends_B_virtual, sizes_array, &
                                                            sizes_B_virtual, starts_array, &
                                                            starts_B_virtual
      INTEGER :: homo, dimen_RI, unit_nr, color_sub, best_block_size, best_integ_group_size, &
         block_size, integ_group_size, min_integ_group_size, my_B_size, my_B_virtual_end, &
         my_B_virtual_start, my_group_L_size, my_group_L_start, my_group_L_end, ngroup, &
         num_IJ_blocks, num_integ_group, pos_integ_group, virtual
      LOGICAL                                            :: my_alpha_beta_case, my_open_shell_SS
      REAL(KIND=dp)                                      :: mem_for_aK, mem_for_comm, mem_for_iaK, &
                                                            mem_for_rep, mem_min, mem_per_group, &
                                                            mem_real

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_get_sizes', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iiB
      INTEGER(KIND=int_8)                                :: mem

      CALL timeset(routineN, handle)

      ngroup = para_env%num_pe/para_env_sub%num_pe

      ! Calculate available memory and create integral group according to that
      ! mem_for_iaK is the memory needed for storing the 3 centre integrals
      mem_for_iaK = REAL(homo, KIND=dp)*virtual*dimen_RI*8.0_dp/(1024_dp**2)
      mem_for_aK = REAL(virtual, KIND=dp)*dimen_RI*8.0_dp/(1024_dp**2)

      CALL m_memory(mem)
      mem_real = (mem+1024*1024-1)/(1024*1024)
      ! mp_min .... a hack.. it should be mp_max, but as it turns out, on some processes the previously freed memory (hfx)
      ! has not been given back to the OS yet.
      CALL mp_min(mem_real, para_env%group)

      mem_min = 2.0_dp*REAL(homo, KIND=dp)*MAXVAL(sizes_B_virtual)*MAXVAL(sizes_array)*8.0_dp/(1024**2)
      mem_min = mem_min+3.0_dp*MAXVAL(sizes_B_virtual)*REAL(dimen_RI, KIND=dp)*8.0_dp/(1024**2)

      IF ((.NOT. my_open_shell_SS) .AND. (.NOT. my_alpha_beta_case)) THEN
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T68,F9.2,A4)') 'RI_INFO| Minimum required memory per MPI process:', &
            mem_min, ' MiB'
      END IF

      mem_real = mp2_env%mp2_memory

      mem_per_group = mem_real*para_env_sub%num_pe

      ! here we try to find the best block_size and integ_group_size
      best_integ_group_size = ngroup
      best_block_size = 1

      ! in the open shell case no replication and no block communication is done
      IF ((.NOT. my_open_shell_SS) .AND. (.NOT. my_alpha_beta_case)) THEN
         ! Here we split the memory half for the communication, half for replication
         IF (mp2_env%ri_mp2%block_size > 0) THEN
            best_block_size = mp2_env%ri_mp2%block_size
            mem_for_rep = MAX(mem_min, mem_per_group-2.0_dp*mem_for_aK*best_block_size)
         ELSE
            mem_for_rep = mem_per_group/2.0_dp
         END IF
         ! calculate the minimum replication group size according to the available memory
         min_integ_group_size = CEILING(2.0_dp*mem_for_iaK/mem_for_rep)

         integ_group_size = MIN(min_integ_group_size, ngroup)-1
         DO iiB = min_integ_group_size+1, ngroup
            integ_group_size = integ_group_size+1
            ! check that the ngroup is a multiple of  integ_group_size
            IF (MOD(ngroup, integ_group_size) /= 0) CYCLE
            ! check that the integ group size is not too small (10% is empirical for now)
            IF (REAL(integ_group_size, KIND=dp)/REAL(ngroup, KIND=dp) < 0.1_dp) CYCLE

            best_integ_group_size = integ_group_size
            EXIT
         END DO

         IF (.NOT. (mp2_env%ri_mp2%block_size > 0)) THEN
            mem_for_comm = mem_per_group-2.0_dp*mem_for_iaK/best_integ_group_size
            !k_size=MAX(INT(mem_for_comm/(2.0_dp*mem_for_aK)),1)
            DO
               num_IJ_blocks = (homo/best_block_size)
               num_IJ_blocks = (num_IJ_blocks*num_IJ_blocks-num_IJ_blocks)/2
               IF (num_IJ_blocks > ngroup .OR. best_block_size == 1) THEN
                  EXIT
               ELSE
                  best_block_size = best_block_size-1
               END IF
            END DO
         END IF

         ! check that best_block_size is not bigger than homo/2-1
         best_block_size = MIN(MAX(homo/2-1+MOD(homo, 2), 1), best_block_size)
      END IF

      integ_group_size = best_integ_group_size
      block_size = best_block_size

      IF ((.NOT. my_open_shell_SS) .AND. (.NOT. my_alpha_beta_case)) THEN
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "RI_INFO| Group size for integral replication:", integ_group_size*para_env_sub%num_pe
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "RI_INFO| Block size:", block_size
            CALL m_flush(unit_nr)
         END IF
      END IF

      num_integ_group = ngroup/integ_group_size

      pos_integ_group = MOD(color_sub, integ_group_size)

      my_group_L_size = sizes_array(color_sub)
      my_group_L_start = starts_array(color_sub)
      my_group_L_end = ends_array(color_sub)

      my_B_virtual_start = starts_B_virtual(para_env_sub%mepos)
      my_B_virtual_end = ends_B_virtual(para_env_sub%mepos)
      my_B_size = sizes_B_virtual(para_env_sub%mepos)

      CALL timestop(handle)

   END SUBROUTINE mp2_ri_get_sizes

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param para_env_sub ...
!> \param ends_B_virtual ...
!> \param sizes_B_virtual ...
!> \param starts_B_virtual ...
!> \param Eigenval ...
!> \param homo ...
!> \param dimen_RI ...
!> \param iiB ...
!> \param jjB ...
!> \param my_B_size ...
!> \param my_B_virtual_end ...
!> \param my_B_virtual_start ...
!> \param my_i ...
!> \param my_j ...
!> \param virtual ...
!> \param sub_proc_map ...
!> \param local_ab ...
!> \param t_ab ...
!> \param local_i_aL ...
!> \param local_j_aL ...
!> \param open_ss ...
!> \param alpha_alpha ...
!> \param beta_beta ...
!> \param Y_i_aP ...
!> \param Y_j_aP ...
!> \param eigenval_beta ...
!> \param homo_beta ...
!> \param my_B_size_beta ...
!> \param ends_B_virtual_beta ...
!> \param sizes_B_virtual_beta ...
!> \param starts_B_virtual_beta ...
!> \param my_B_virtual_start_beta ...
!> \param my_B_virtual_end_beta ...
!> \param virtual_beta ...
!> \param local_ba ...
! **************************************************************************************************
   SUBROUTINE mp2_update_P_gamma(mp2_env, para_env_sub, ends_B_virtual, sizes_B_virtual, starts_B_virtual, &
                                 Eigenval, homo, dimen_RI, iiB, jjB, my_B_size, &
                                 my_B_virtual_end, my_B_virtual_start, my_i, my_j, virtual, sub_proc_map, local_ab, &
                                 t_ab, local_i_aL, local_j_aL, open_ss, alpha_alpha, beta_beta, Y_i_aP, Y_j_aP, &
                                 eigenval_beta, homo_beta, my_B_size_beta, ends_B_virtual_beta, sizes_B_virtual_beta, &
                                 starts_B_virtual_beta, my_B_virtual_start_beta, my_B_virtual_end_beta, virtual_beta, local_ba)
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ends_B_virtual, sizes_B_virtual, &
                                                            starts_B_virtual
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: homo, dimen_RI, iiB, jjB, my_B_size, &
                                                            my_B_virtual_end, my_B_virtual_start, &
                                                            my_i, my_j, virtual
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sub_proc_map
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: local_ab, t_ab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: local_i_aL, local_j_aL
      LOGICAL                                            :: open_ss, alpha_alpha, beta_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Y_i_aP, Y_j_aP
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: Eigenval_beta
      INTEGER, OPTIONAL                                  :: homo_beta, my_B_size_beta
      INTEGER, ALLOCATABLE, DIMENSION(:), OPTIONAL       :: ends_B_virtual_beta, &
                                                            sizes_B_virtual_beta, &
                                                            starts_B_virtual_beta
      INTEGER, OPTIONAL                                  :: my_B_virtual_start_beta, &
                                                            my_B_virtual_end_beta, virtual_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         OPTIONAL                                        :: local_ba

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_update_P_gamma', &
         routineP = moduleN//':'//routineN

      INTEGER :: a, b, b_global, handle, proc_receive, proc_send, proc_shift, rec_B_size, &
         rec_B_virtual_end, rec_B_virtual_start, send_B_size, send_B_virtual_end, &
         send_B_virtual_start
      LOGICAL                                            :: alpha_beta
      REAL(KIND=dp)                                      :: factor, P_ij_diag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: external_ab, send_ab

!
!  In alpha-beta case Y_j_aP_beta is sent and is received as Y_j_aP
!

      CALL timeset(routineN//"_Pia", handle)
! Find out whether we have an alpha-beta case
      alpha_beta = .FALSE.
      IF (PRESENT(Eigenval_beta) .AND. PRESENT(homo_beta) .AND. PRESENT(my_B_size_beta)) &
         alpha_beta = .TRUE.
! update P_ab, Gamma_P_ia
! First, P_ab
      IF (open_ss) THEN
         factor = 1.0_dp
      ELSE
         factor = 2.0_dp
      ENDIF
      ! divide the (ia|jb) integrals by Delta_ij^ab
      IF (.NOT. alpha_beta) THEN
         DO b = 1, my_B_size
            b_global = b+my_B_virtual_start-1
            DO a = 1, virtual
               local_ab(a, b) = -local_ab(a, b)/ &
                                (Eigenval(homo+a)+Eigenval(homo+b_global)-Eigenval(my_i+iiB-1)-Eigenval(my_j+jjB-1))
            END DO
         END DO
         ! update diagonal part of P_ij
         P_ij_diag = -SUM(local_ab*t_ab)*factor
      ELSE
         DO b = 1, my_B_size_beta
            b_global = b+my_B_virtual_start_beta-1
            DO a = 1, virtual
               local_ab(a, b) = -local_ab(a, b)/ &
                                (Eigenval(homo+a)+Eigenval_beta(homo_beta+b_global)-Eigenval(my_i+iiB-1)-Eigenval_beta(my_j+jjB-1))
            END DO
         END DO
         ! update diagonal part of P_ij
         P_ij_diag = -SUM(local_ab*local_ab)
         ! More integrals needed only for alpha-beta case: local_ba
         DO b = 1, my_B_size
            b_global = b+my_B_virtual_start-1
            DO a = 1, virtual_beta
               local_ba(a, b) = -local_ba(a, b)/ &
                                (Eigenval_beta(homo_beta+a)+Eigenval(homo+b_global)-Eigenval(my_i+iiB-1)-Eigenval_beta(my_j+jjB-1))
            END DO
         END DO
      ENDIF

      ! P_ab and add diagonal part of P_ij

      ! Alpha_alpha or closed-shell case
      IF (((.NOT. open_ss) .AND. (.NOT. alpha_beta)) .OR. alpha_alpha) THEN
         CALL dgemm('T', 'N', my_B_size, my_B_size, virtual, 1.0_dp, &
                    t_ab(:, :), virtual, local_ab(:, :), virtual, &
                    1.0_dp, mp2_env%ri_grad%P_ab(1:my_B_size, my_B_virtual_start:my_B_virtual_end), my_B_size)
         mp2_env%ri_grad%P_ij(my_i+iiB-1, my_i+iiB-1) = &
            mp2_env%ri_grad%P_ij(my_i+iiB-1, my_i+iiB-1)+P_ij_diag
      ENDIF
      ! Beta_beta case
      IF (beta_beta) THEN
         CALL dgemm('T', 'N', my_B_size, my_B_size, virtual, 1.0_dp, &
                    t_ab(:, :), virtual, local_ab(:, :), virtual, &
                    1.0_dp, mp2_env%ri_grad%P_ab_beta(1:my_B_size, my_B_virtual_start:my_B_virtual_end), my_B_size)
         mp2_env%ri_grad%P_ij_beta(my_i+iiB-1, my_i+iiB-1) = &
            mp2_env%ri_grad%P_ij_beta(my_i+iiB-1, my_i+iiB-1)+P_ij_diag
      ENDIF
      ! Alpha_beta case
      IF (alpha_beta) THEN
         CALL dgemm('T', 'N', my_B_size, my_B_size, virtual_beta, 1.0_dp, &
                    local_ba(:, :), virtual_beta, local_ba(:, :), virtual_beta, 1.0_dp, &
                    mp2_env%ri_grad%P_ab(1:my_B_size, my_B_virtual_start: &
                                         my_B_virtual_end), my_B_size)
         mp2_env%ri_grad%P_ij(my_i+iiB-1, my_i+iiB-1) = &
            mp2_env%ri_grad%P_ij(my_i+iiB-1, my_i+iiB-1)+P_ij_diag
         CALL dgemm('T', 'N', my_B_size_beta, my_B_size_beta, virtual, 1.0_dp, &
                    local_ab(:, :), virtual, local_ab(:, :), virtual, 1.0_dp, &
                    mp2_env%ri_grad%P_ab_beta(1:my_B_size_beta, my_B_virtual_start_beta: &
                                              my_B_virtual_end_beta), my_B_size_beta)
         mp2_env%ri_grad%P_ij_beta(my_j+jjB-1, my_j+jjB-1) = &
            mp2_env%ri_grad%P_ij_beta(my_j+jjB-1, my_j+jjB-1)+P_ij_diag
      ENDIF
      ! The summation is over unique pairs. In alpha-beta case, all pairs are unique: subroutine is called for
      ! both i^alpha,j^beta and i^beta,j^alpha. Formally, my_i can be equal to my_j, but they are different
      ! due to spin in alpha-beta case.
      IF ((my_i /= my_j) .AND. (.NOT. alpha_beta)) THEN
         ! CALL dgemm('N','T',virtual,virtual,my_B_size,1.0_dp,&
         !                t_ab(:,:),virtual,local_ab(:,:),virtual,&
         !                1.0_dp,mp2_env%ri_grad%P_ab(:,:),virtual)

         ! Alpha_alpha or closed-shell case
         IF (((.NOT. open_ss) .AND. (.NOT. alpha_beta)) .OR. alpha_alpha) THEN
            CALL dgemm('N', 'T', my_B_size, virtual, my_B_size, 1.0_dp, &
                       t_ab(my_B_virtual_start:my_B_virtual_end, 1:my_B_size), my_B_size, &
                       local_ab(:, :), virtual, &
                       1.0_dp, mp2_env%ri_grad%P_ab(1:my_B_size, 1:virtual), my_B_size)
            ! CHECK CORRECTNESS IF BLOCK SIZE greater than 1
            mp2_env%ri_grad%P_ij(my_j+jjB-1, my_j+jjB-1) = &
               mp2_env%ri_grad%P_ij(my_j+jjB-1, my_j+jjB-1)+P_ij_diag
         ENDIF
         ! Beta_beta_case
         IF (beta_beta) THEN
            CALL dgemm('N', 'T', my_B_size, virtual, my_B_size, 1.0_dp, &
                       t_ab(my_B_virtual_start:my_B_virtual_end, 1:my_B_size), my_B_size, &
                       local_ab(:, :), virtual, &
                       1.0_dp, mp2_env%ri_grad%P_ab_beta(1:my_B_size, 1:virtual), my_B_size)
            ! CHECK CORRECTNESS IF BLOCK SIZE greater than 1
            mp2_env%ri_grad%P_ij_beta(my_j+jjB-1, my_j+jjB-1) = &
               mp2_env%ri_grad%P_ij_beta(my_j+jjB-1, my_j+jjB-1)+P_ij_diag
         ENDIF
      END IF
      DO proc_shift = 1, para_env_sub%num_pe-1
         proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
         proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)
         ! Alpha-alpha, beta-beta, closed shell
         IF (.NOT. alpha_beta) THEN
            rec_B_size = sizes_B_virtual(proc_receive)
            rec_B_virtual_end = ends_B_virtual(proc_receive)
            rec_B_virtual_start = starts_B_virtual(proc_receive)
            send_B_size = sizes_B_virtual(proc_send)
            send_B_virtual_end = ends_B_virtual(proc_send)
            send_B_virtual_start = starts_B_virtual(proc_send)
         ELSE ! Alpha-beta case
            rec_B_size = sizes_B_virtual_beta(proc_receive)
            rec_B_virtual_end = ends_B_virtual_beta(proc_receive)
            rec_B_virtual_start = starts_B_virtual_beta(proc_receive)
            send_B_size = sizes_B_virtual_beta(proc_send)
            send_B_virtual_end = ends_B_virtual_beta(proc_send)
            send_B_virtual_start = starts_B_virtual_beta(proc_send)
         ENDIF

         ALLOCATE (external_ab(virtual, rec_B_size))
         external_ab = 0.0_dp

         IF (.NOT. alpha_beta) THEN
            CALL mp_sendrecv(local_ab(1:virtual, 1:my_B_size), proc_send, &
                             external_ab(1:virtual, 1:rec_B_size), proc_receive, &
                             para_env_sub%group)
         ELSE
            CALL mp_sendrecv(local_ab(1:virtual, 1:my_B_size_beta), proc_send, &
                             external_ab(1:virtual, 1:rec_B_size), proc_receive, &
                             para_env_sub%group)
         ENDIF

         ! Alpha-alpha or closed-shell case
         IF (((.NOT. open_ss) .AND. (.NOT. alpha_beta)) .OR. alpha_alpha) &
            CALL dgemm('T', 'N', my_B_size, rec_B_size, virtual, 1.0_dp, &
                       t_ab(:, :), virtual, external_ab(:, :), virtual, &
                       1.0_dp, mp2_env%ri_grad%P_ab(1:my_B_size, rec_B_virtual_start:rec_B_virtual_end), my_B_size)
         ! Beta-beta case
         IF (beta_beta) &
            CALL dgemm('T', 'N', my_B_size, rec_B_size, virtual, 1.0_dp, &
                       t_ab(:, :), virtual, external_ab(:, :), virtual, &
                       1.0_dp, mp2_env%ri_grad%P_ab_beta(1:my_B_size, rec_B_virtual_start:rec_B_virtual_end), my_B_size)
         ! Alpha-beta case
         IF (alpha_beta) THEN
            !   CALL dgemm('N','T',my_B_size,virtual,rec_B_size,1.0_dp,&
            !             local_ab(1:my_B_size,rec_B_virtual_start:rec_B_virtual_end),my_B_size,external_ab(:,:),rec_B_size,&
            !             1.0_dp,mp2_env%ri_grad%P_ab(1:my_B_size,1:virtual),my_B_size)

            ! Alpha-beta part of beta-beta density
            CALL dgemm('T', 'N', my_B_size_beta, rec_B_size, virtual, 1.0_dp, &
                       local_ab(:, :), virtual, external_ab(:, :), virtual, &
                       1.0_dp, mp2_env%ri_grad%P_ab_beta(1:my_B_size_beta, rec_B_virtual_start:rec_B_virtual_end), &
                       my_B_size_beta)

            ! For alpha-beta part of alpha-density we need a new parallel code
            DEALLOCATE (external_ab)
            ! And new external_ab (of a different size)
            rec_B_size = sizes_B_virtual(proc_receive)
            rec_B_virtual_end = ends_B_virtual(proc_receive)
            rec_B_virtual_start = starts_B_virtual(proc_receive)
            send_B_size = sizes_B_virtual(proc_send)
            send_B_virtual_end = ends_B_virtual(proc_send)
            send_B_virtual_start = starts_B_virtual(proc_send)
            ALLOCATE (external_ab(virtual_beta, rec_B_size))
            external_ab = 0.0_dp
            CALL mp_sendrecv(local_ba(1:virtual_beta, 1:my_B_size), proc_send, &
                             external_ab(1:virtual_beta, 1:rec_B_size), proc_receive, &
                             para_env_sub%group)
            CALL dgemm('T', 'N', my_B_size, rec_B_size, virtual_beta, 1.0_dp, &
                       local_ba(:, :), virtual_beta, external_ab(:, :), virtual_beta, &
                       1.0_dp, mp2_env%ri_grad%P_ab(1:my_B_size, rec_B_virtual_start:rec_B_virtual_end), my_B_size)
         ENDIF

         DEALLOCATE (external_ab)

         IF ((my_i /= my_j) .AND. (.NOT. alpha_beta)) THEN
            ALLOCATE (external_ab(my_B_size, virtual))
            external_ab = 0.0_dp

            ALLOCATE (send_ab(send_B_size, virtual))
            send_ab = 0.0_dp

            CALL dgemm('N', 'T', send_B_size, virtual, my_B_size, 1.0_dp, &
                       t_ab(send_B_virtual_start:send_B_virtual_end, 1:my_B_size), send_B_size, &
                       local_ab(:, :), virtual, &
                       0.0_dp, send_ab(1:send_B_size, 1:virtual), send_B_size)

            CALL mp_sendrecv(send_ab, proc_send, &
                             external_ab, proc_receive, &
                             para_env_sub%group)

            ! Alpha_alpha or closed-shell case
            IF (((.NOT. open_ss) .AND. (.NOT. alpha_beta)) .OR. alpha_alpha) &
               mp2_env%ri_grad%P_ab(:, :) = mp2_env%ri_grad%P_ab+external_ab
            ! Beta_beta case
            IF (beta_beta) &
               mp2_env%ri_grad%P_ab_beta(:, :) = mp2_env%ri_grad%P_ab_beta+external_ab

            DEALLOCATE (external_ab)
            DEALLOCATE (send_ab)
         END IF

      END DO
      CALL timestop(handle)

      ! Now, Gamma_P_ia (made of Y_ia_P)

      CALL timeset(routineN//"_Gamma", handle)
      ! CALL dgemm('T','N',my_B_size,my_B_size_beta,dimen_RI,1.0_dp,&
      !                local_i_aL(:,:,iiB),dimen_RI,local_j_aL(:,:,jjB),dimen_RI,&
      !                0.0_dp,local_ab(my_B_virtual_start:my_B_virtual_end,1:my_B_size_beta),my_B_size)
      IF (.NOT. alpha_beta) THEN
         ! Alpha-alpha, beta-beta and closed shell
         CALL dgemm('N', 'T', my_B_size, dimen_RI, my_B_size, 1.0_dp, &
                    t_ab(my_B_virtual_start:my_B_virtual_end, 1:my_B_size), my_B_size, &
                    local_j_aL(1:dimen_RI, 1:my_B_size, jjB), dimen_RI, &
                    1.0_dp, Y_i_aP(1:my_B_size, 1:dimen_RI, iiB), my_B_size)
      ELSE ! Alpha-beta
         CALL dgemm('N', 'T', my_B_size, dimen_RI, my_B_size_beta, 1.0_dp, &
                    local_ab(my_B_virtual_start:my_B_virtual_end, 1:my_B_size_beta), my_B_size, &
                    local_j_aL(1:dimen_RI, 1:my_B_size_beta, jjB), dimen_RI, &
                    1.0_dp, Y_i_aP(1:my_B_size, 1:dimen_RI, iiB), my_B_size)
         CALL dgemm('T', 'T', my_B_size_beta, dimen_RI, my_B_size, 1.0_dp, &
                    local_ab(my_B_virtual_start:my_B_virtual_end, 1:my_B_size_beta), my_B_size, &
                    local_i_aL(1:dimen_RI, 1:my_B_size, jjB), dimen_RI, &
                    1.0_dp, Y_j_aP(1:my_B_size_beta, 1:dimen_RI, jjB), my_B_size_beta)

      ENDIF

      ALLOCATE (external_ab(my_B_size, dimen_RI))
      external_ab = 0.0_dp
      !
      DO proc_shift = 1, para_env_sub%num_pe-1
         proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
         proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)

         rec_B_size = sizes_B_virtual(proc_receive)
         rec_B_virtual_end = ends_B_virtual(proc_receive)
         rec_B_virtual_start = starts_B_virtual(proc_receive)

         send_B_size = sizes_B_virtual(proc_send)
         send_B_virtual_end = ends_B_virtual(proc_send)
         send_B_virtual_start = starts_B_virtual(proc_send)

         IF (.NOT. alpha_beta) THEN
            ALLOCATE (send_ab(send_B_size, dimen_RI))
            send_ab = 0.0_dp
            CALL dgemm('N', 'T', send_B_size, dimen_RI, my_B_size, 1.0_dp, &
                       t_ab(send_B_virtual_start:send_B_virtual_end, 1:my_B_size), send_B_size, &
                       local_j_aL(1:dimen_RI, 1:my_B_size, jjB), dimen_RI, &
                       0.0_dp, send_ab(1:send_B_size, 1:dimen_RI), send_B_size)
            CALL mp_sendrecv(send_ab, proc_send, external_ab, proc_receive, para_env_sub%group)
            ! Alpha-alpha, beta-beta and closed shell

            Y_i_aP(1:my_B_size, 1:dimen_RI, iiB) = &
               Y_i_aP(1:my_B_size, 1:dimen_RI, iiB)+external_ab

            DEALLOCATE (send_ab)
         ELSE ! Alpha-beta case
            ! Alpha-alpha part
            ALLOCATE (send_ab(send_B_size, dimen_RI))
            send_ab = 0.0_dp
            CALL dgemm('N', 'T', send_B_size, dimen_RI, my_B_size_beta, 1.0_dp, &
                       local_ab(send_B_virtual_start:send_B_virtual_end, 1:my_B_size_beta), send_B_size, &
                       local_j_aL(1:dimen_RI, 1:my_B_size_beta, jjB), dimen_RI, &
                       0.0_dp, send_ab(1:send_B_size, 1:dimen_RI), send_B_size)
            CALL mp_sendrecv(send_ab, proc_send, external_ab, proc_receive, para_env_sub%group)
            Y_i_aP(1:my_B_size, 1:dimen_RI, iiB) = &
               Y_i_aP(1:my_B_size, 1:dimen_RI, iiB)+external_ab
            DEALLOCATE (send_ab)
         ENDIF
      END DO
      DEALLOCATE (external_ab)

      IF (alpha_beta) THEN
         ! For beta-beta part (in alpha-beta case) we need a new parallel code
         ALLOCATE (external_ab(my_B_size_beta, dimen_RI))
         external_ab = 0.0_dp
         DO proc_shift = 1, para_env_sub%num_pe-1
            proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
            proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)

            send_B_size = sizes_B_virtual_beta(proc_send)
            send_B_virtual_end = ends_B_virtual_beta(proc_send)
            send_B_virtual_start = starts_B_virtual_beta(proc_send)
            ALLOCATE (send_ab(send_B_size, dimen_RI))
            send_ab = 0.0_dp

            CALL dgemm('N', 'T', send_B_size, dimen_RI, my_B_size, 1.0_dp, &
                       local_ba(send_B_virtual_start:send_B_virtual_end, 1:my_B_size), send_B_size, &
                       local_i_aL(1:dimen_RI, 1:my_B_size, jjB), dimen_RI, &
                       0.0_dp, send_ab(1:send_B_size, 1:dimen_RI), send_B_size)
            CALL mp_sendrecv(send_ab, proc_send, external_ab, proc_receive, para_env_sub%group)
            Y_j_aP(1:my_B_size_beta, 1:dimen_RI, jjB) = &
               Y_j_aP(1:my_B_size_beta, 1:dimen_RI, jjB)+external_ab
            DEALLOCATE (send_ab)

         END DO
         DEALLOCATE (external_ab)
      ENDIF

      IF ((my_i /= my_j) .AND. (.NOT. alpha_beta)) THEN
         ! Alpha-alpha, beta-beta and closed shell
         CALL dgemm('T', 'T', my_B_size, dimen_RI, my_B_size, 1.0_dp, &
                    t_ab(my_B_virtual_start:my_B_virtual_end, 1:my_B_size), my_B_size, &
                    local_i_aL(1:dimen_RI, 1:my_B_size, iiB), dimen_RI, &
                    1.0_dp, Y_j_aP(1:my_B_size, 1:dimen_RI, jjB), my_B_size)

         DO proc_shift = 1, para_env_sub%num_pe-1
            proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
            proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)

            rec_B_size = sizes_B_virtual(proc_receive)
            rec_B_virtual_end = ends_B_virtual(proc_receive)
            rec_B_virtual_start = starts_B_virtual(proc_receive)

            ALLOCATE (external_ab(dimen_RI, rec_B_size))
            external_ab = 0.0_dp

            CALL mp_sendrecv(local_i_aL(1:dimen_RI, 1:my_B_size, iiB), proc_send, &
                             external_ab, proc_receive, para_env_sub%group)

            ! Alpha-alpha, beta-beta and closed shell
            CALL dgemm('T', 'T', my_B_size, dimen_RI, rec_B_size, 1.0_dp, &
                       t_ab(rec_B_virtual_start:rec_B_virtual_end, 1:my_B_size), rec_B_size, &
                       external_ab(1:dimen_RI, 1:rec_B_size), dimen_RI, &
                       1.0_dp, Y_j_aP(1:my_B_size, 1:dimen_RI, jjB), my_B_size)

            DEALLOCATE (external_ab)
         END DO
      END IF

      ! Y_i_aP=0.0_dp
      ! Y_j_aP=0.0_dp

      CALL timestop(handle)

   END SUBROUTINE mp2_update_P_gamma

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param ij_index ...
!> \param my_B_size ...
!> \param my_block_size ...
!> \param my_group_L_size ...
!> \param my_i ...
!> \param my_ij_pairs ...
!> \param my_j ...
!> \param ngroup ...
!> \param num_integ_group ...
!> \param integ_group_pos2color_sub ...
!> \param num_ij_pairs ...
!> \param proc_map ...
!> \param ij_map ...
!> \param ranges_info_array ...
!> \param Y_i_aP ...
!> \param Y_j_aP ...
!> \param para_env_exchange ...
!> \param null_mat_rec ...
!> \param null_mat_send ...
!> \param sizes_array ...
!> \param alpha_alpha ...
!> \param beta_beta ...
!> \param alpha_beta ...
!> \param open_shell ...
!> \param my_b_size_beta ...
! **************************************************************************************************
   SUBROUTINE mp2_redistribute_gamma(mp2_env, ij_index, my_B_size, &
                                     my_block_size, my_group_L_size, my_i, my_ij_pairs, my_j, ngroup, &
                                     num_integ_group, integ_group_pos2color_sub, num_ij_pairs, proc_map, &
                                     ij_map, ranges_info_array, Y_i_aP, Y_j_aP, para_env_exchange, &
                                     null_mat_rec, null_mat_send, sizes_array, alpha_alpha, beta_beta, &
                                     alpha_beta, open_shell, my_b_size_beta)

      TYPE(mp2_type), POINTER                            :: mp2_env
      INTEGER                                            :: ij_index, my_B_size, my_block_size, &
                                                            my_group_L_size, my_i, my_ij_pairs, &
                                                            my_j, ngroup, num_integ_group
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: integ_group_pos2color_sub, num_ij_pairs, &
                                                            proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ij_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: ranges_info_array
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Y_i_aP, Y_j_aP
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange
      REAL(KIND=dp)                                      :: null_mat_rec(:, :, :), &
                                                            null_mat_send(:, :, :)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_array
      LOGICAL                                            :: alpha_alpha, beta_beta, alpha_beta, &
                                                            open_shell
      INTEGER, OPTIONAL                                  :: my_B_size_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_redistribute_gamma', &
         routineP = moduleN//':'//routineN

      INTEGER :: end_point, handle, handle2, iiB, ij_counter_rec, irep, jjb, kkk, Lend_pos, lll, &
         Lstart_pos, proc_receive, proc_send, proc_shift, rec_block_size, rec_i, rec_ij_index, &
         rec_j, send_L_size, start_point
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BI_C_rec, BI_C_rec_beta, BI_C_send, &
                                                            BI_C_send_beta

! In alpha-beta case Y_i_aP_beta is sent as Y_j_aP

      CALL timeset(routineN//"_comm2", handle)
      IF (ij_index <= my_ij_pairs) THEN
         ! somethig to send
         ! start with myself
         CALL timeset(routineN//"_comm2_w", handle2)
         DO irep = 0, num_integ_group-1
            Lstart_pos = ranges_info_array(1, irep, para_env_exchange%mepos)
            Lend_pos = ranges_info_array(2, irep, para_env_exchange%mepos)
            start_point = ranges_info_array(3, irep, para_env_exchange%mepos)
            end_point = ranges_info_array(4, irep, para_env_exchange%mepos)
            DO iiB = 1, my_block_size
!$OMP           PARALLEL DO DEFAULT(NONE) PRIVATE(kkk,lll) &
!$OMP                                     SHARED(start_point,end_point,Lstart_pos,Lend_pos,&
!$OMP                                            mp2_env,my_i,iiB,my_B_size,Y_i_aP,&
!$OMP                                     alpha_alpha,beta_beta,open_shell)
               DO kkk = start_point, end_point
                  lll = kkk-start_point+Lstart_pos
                  IF (alpha_alpha .OR. (.NOT. open_shell)) THEN
                     mp2_env%ri_grad%Gamma_P_ia(my_i+iiB-1, 1:my_B_size, kkk) = &
                        mp2_env%ri_grad%Gamma_P_ia(my_i+iiB-1, 1:my_B_size, kkk)+ &
                        Y_i_aP(1:my_B_size, lll, iiB)
                  ENDIF
                  IF (beta_beta) THEN
                     mp2_env%ri_grad%Gamma_P_ia_beta(my_i+iiB-1, 1:my_B_size, kkk) = &
                        mp2_env%ri_grad%Gamma_P_ia_beta(my_i+iiB-1, 1:my_B_size, kkk)+ &
                        Y_i_aP(1:my_B_size, lll, iiB)
                  ENDIF
               END DO
!$OMP           END PARALLEL DO
            END DO
            DO jjB = 1, my_block_size
               IF (.NOT. alpha_beta) THEN
!$OMP              PARALLEL DO DEFAULT(NONE) PRIVATE(kkk,lll) &
!$OMP                                        SHARED(start_point,end_point,Lstart_pos,Lend_pos,&
!$OMP                                               mp2_env,my_j,jjB,my_B_size,Y_j_aP,&
!$OMP                                        alpha_alpha,beta_beta,open_shell)
                  DO kkk = start_point, end_point
                     lll = kkk-start_point+Lstart_pos
                     IF (alpha_alpha .OR. (.NOT. open_shell)) THEN
                        mp2_env%ri_grad%Gamma_P_ia(my_j+jjB-1, 1:my_B_size, kkk) = &
                           mp2_env%ri_grad%Gamma_P_ia(my_j+jjB-1, 1:my_B_size, kkk)+ &
                           Y_j_aP(1:my_B_size, lll, jjB)
                     ENDIF
                     IF (beta_beta) THEN
                        mp2_env%ri_grad%Gamma_P_ia_beta(my_j+jjB-1, 1:my_B_size, kkk) = &
                           mp2_env%ri_grad%Gamma_P_ia_beta(my_j+jjB-1, 1:my_B_size, kkk)+ &
                           Y_j_aP(1:my_B_size, lll, jjB)
                     ENDIF
                  END DO
!$OMP              END PARALLEL DO
               ELSE
!$OMP              PARALLEL DO DEFAULT(NONE) PRIVATE(kkk,lll) &
!$OMP                                        SHARED(start_point,end_point,Lstart_pos,Lend_pos,&
!$OMP                                               mp2_env,my_j,jjB,my_B_size_beta,Y_j_aP)
                  DO kkk = start_point, end_point
                     lll = kkk-start_point+Lstart_pos
                     mp2_env%ri_grad%Gamma_P_ia_beta(my_j+jjB-1, 1:my_B_size_beta, kkk) = &
                        mp2_env%ri_grad%Gamma_P_ia_beta(my_j+jjB-1, 1:my_B_size_beta, kkk)+ &
                        Y_j_aP(1:my_B_size_beta, lll, jjB)
                  ENDDO
!$OMP              END PARALLEL DO
               ENDIF
            END DO
            ! DO iiB=1, my_block_size
            !   mp2_env%ri_grad%Gamma_P_ia(my_i+iiB-1,1:my_B_size,start_point:end_point)=&
            !   mp2_env%ri_grad%Gamma_P_ia(my_i+iiB-1,1:my_B_size,start_point:end_point)+&
            !   Y_i_aP(1:my_B_size,Lstart_pos:Lend_pos,iiB)
            ! END DO
            ! DO jjB=1, my_block_size
            !   mp2_env%ri_grad%Gamma_P_ia(my_j+jjB-1,1:my_B_size,start_point:end_point)=&
            !   mp2_env%ri_grad%Gamma_P_ia(my_j+jjB-1,1:my_B_size,start_point:end_point)+&
            !   Y_j_aP(1:my_B_size,Lstart_pos:Lend_pos,jjB)
            ! END DO
         END DO
         CALL timestop(handle2)

         ! Y_i_aP(my_B_size,dimen_RI,block_size)

         DO proc_shift = 1, para_env_exchange%num_pe-1
            proc_send = proc_map(para_env_exchange%mepos+proc_shift)
            proc_receive = proc_map(para_env_exchange%mepos-proc_shift)

            send_L_size = sizes_array(proc_send)
            IF (.NOT. alpha_beta) THEN
               ALLOCATE (BI_C_send(2*my_block_size, my_B_size, send_L_size))
            ELSE
               ALLOCATE (BI_C_send(my_block_size, my_B_size, send_L_size))
               ALLOCATE (BI_C_send_beta(my_block_size, my_B_size_beta, send_L_size))
            ENDIF
            CALL timeset(routineN//"_comm2_w", handle2)
            BI_C_send = 0.0_dp
            IF (alpha_beta) BI_C_send_beta = 0.0_dp
            DO irep = 0, num_integ_group-1
               Lstart_pos = ranges_info_array(1, irep, proc_send)
               Lend_pos = ranges_info_array(2, irep, proc_send)
               start_point = ranges_info_array(3, irep, proc_send)
               end_point = ranges_info_array(4, irep, proc_send)
               DO iiB = 1, my_block_size
!$OMP             PARALLEL DO DEFAULT(NONE) PRIVATE(kkk,lll) &
!$OMP                                       SHARED(start_point,end_point,Lstart_pos,Lend_pos,&
!$OMP                                              BI_C_send,iiB,my_B_size,Y_i_aP)
                  DO kkk = start_point, end_point
                     lll = kkk-start_point+Lstart_pos
                     BI_C_send(iiB, 1:my_B_size, kkk) = Y_i_aP(1:my_B_size, lll, iiB)
                  END DO
!$OMP             END PARALLEL DO
               END DO
               DO jjB = 1, my_block_size
                  IF (.NOT. alpha_beta) THEN
!$OMP                PARALLEL DO DEFAULT(NONE) PRIVATE(kkk,lll) &
!$OMP                                          SHARED(start_point,end_point,Lstart_pos,Lend_pos,my_block_size,&
!$OMP                                                 BI_C_send,jjB,my_B_size,Y_j_aP)
                     DO kkk = start_point, end_point
                        lll = kkk-start_point+Lstart_pos
                        BI_C_send(jjB+my_block_size, 1:my_B_size, kkk) = Y_j_aP(1:my_B_size, lll, jjB)
                     END DO
!$OMP                END PARALLEL DO
                  ELSE
!$OMP                PARALLEL DO DEFAULT(NONE) PRIVATE(kkk,lll) &
!$OMP                                          SHARED(start_point,end_point,Lstart_pos,Lend_pos,&
!$OMP                                                 BI_C_send_beta,jjB,my_B_size_beta,Y_j_aP)
                     DO kkk = start_point, end_point
                        lll = kkk-start_point+Lstart_pos
                        BI_C_send_beta(jjB, 1:my_B_size_beta, kkk) = Y_j_aP(1:my_B_size_beta, lll, jjB)
                     END DO
!$OMP                END PARALLEL DO
                  ENDIF
               END DO
               ! DO iiB=1, my_block_size
               !   BI_C_send(iiB,1:my_B_size,start_point:end_point)=&
               !      Y_i_aP(1:my_B_size,Lstart_pos:Lend_pos,iiB)
               ! END DO
               ! DO jjB=1, my_block_size
               !   BI_C_send(jjB+my_block_size,1:my_B_size,start_point:end_point)=&
               !      Y_j_aP(1:my_B_size,Lstart_pos:Lend_pos,jjB)
               ! END DO
            END DO
            CALL timestop(handle2)

            rec_ij_index = num_ij_pairs(proc_receive)

            IF (ij_index <= rec_ij_index) THEN
               ! we know that proc_receive has something to send for us, let's see what
               ij_counter_rec = &
                  (ij_index-MIN(1, integ_group_pos2color_sub(proc_receive)))*ngroup+integ_group_pos2color_sub(proc_receive)

               rec_i = ij_map(ij_counter_rec, 1)
               rec_j = ij_map(ij_counter_rec, 2)
               rec_block_size = ij_map(ij_counter_rec, 3)

               IF (.NOT. alpha_beta) THEN
                  ALLOCATE (BI_C_rec(2*rec_block_size, my_B_size, my_group_L_size))
               ELSE
                  ALLOCATE (BI_C_rec(rec_block_size, my_B_size, my_group_L_size))
                  ALLOCATE (BI_C_rec_beta(rec_block_size, my_B_size_beta, my_group_L_size))
               ENDIF

               BI_C_rec = 0.0_dp
               IF (alpha_beta) BI_C_rec_beta = 0.0_dp

               CALL mp_sendrecv(BI_C_send, proc_send, &
                                BI_C_rec, proc_receive, &
                                para_env_exchange%group)
               IF (alpha_beta) THEN
                  CALL mp_sendrecv(BI_C_send_beta, proc_send, &
                                   BI_C_rec_beta, proc_receive, &
                                   para_env_exchange%group)
               ENDIF

               CALL timeset(routineN//"_comm2_w", handle2)
               DO irep = 0, num_integ_group-1
                  Lstart_pos = ranges_info_array(1, irep, para_env_exchange%mepos)
                  Lend_pos = ranges_info_array(2, irep, para_env_exchange%mepos)
                  start_point = ranges_info_array(3, irep, para_env_exchange%mepos)
                  end_point = ranges_info_array(4, irep, para_env_exchange%mepos)
                  DO iiB = 1, rec_block_size
!$OMP               PARALLEL DO DEFAULT(NONE) PRIVATE(kkk) &
!$OMP                                         SHARED(start_point,end_point,&
!$OMP                                                mp2_env,rec_i,iiB,my_B_size,BI_C_rec,&
!$OMP                                         alpha_alpha,beta_beta,open_shell)
                     DO kkk = start_point, end_point
                        IF (alpha_alpha .OR. (.NOT. open_shell)) THEN
                           mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1, 1:my_B_size, kkk) = &
                              mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1, 1:my_B_size, kkk)+ &
                              BI_C_rec(iiB, 1:my_B_size, kkk)
                        ENDIF
                        IF (beta_beta) THEN
                           mp2_env%ri_grad%Gamma_P_ia_beta(rec_i+iiB-1, 1:my_B_size, kkk) = &
                              mp2_env%ri_grad%Gamma_P_ia_beta(rec_i+iiB-1, 1:my_B_size, kkk)+ &
                              BI_C_rec(iiB, 1:my_B_size, kkk)
                        ENDIF
                     END DO
!$OMP               END PARALLEL DO
                  END DO
                  DO jjB = 1, rec_block_size
                     IF (.NOT. alpha_beta) THEN
!$OMP                   PARALLEL DO DEFAULT(NONE) PRIVATE(kkk) &
!$OMP                                             SHARED(start_point,end_point,rec_block_size,&
!$OMP                                                    mp2_env,rec_j,jjB,my_B_size,BI_C_rec,&
!$OMP                                             alpha_alpha,beta_beta,open_shell)
                        DO kkk = start_point, end_point
                           IF (alpha_alpha .OR. (.NOT. open_shell)) THEN
                              mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1, 1:my_B_size, kkk) = &
                                 mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1, 1:my_B_size, kkk)+ &
                                 BI_C_rec(jjB+rec_block_size, 1:my_B_size, kkk)
                           ENDIF
                           IF (beta_beta) THEN
                              mp2_env%ri_grad%Gamma_P_ia_beta(rec_j+jjB-1, 1:my_B_size, kkk) = &
                                 mp2_env%ri_grad%Gamma_P_ia_beta(rec_j+jjB-1, 1:my_B_size, kkk)+ &
                                 BI_C_rec(jjB+rec_block_size, 1:my_B_size, kkk)
                           ENDIF
                        END DO
!$OMP                   END PARALLEL DO
                     ELSE
!$OMP                   PARALLEL DO DEFAULT(NONE) PRIVATE(kkk) &
!$OMP                                             SHARED(start_point,end_point,&
!$OMP                                           mp2_env,rec_j,jjB,my_B_size_beta,BI_C_rec_beta)
                        DO kkk = start_point, end_point
                           mp2_env%ri_grad%Gamma_P_ia_beta(rec_j+jjB-1, 1:my_B_size_beta, kkk) = &
                              mp2_env%ri_grad%Gamma_P_ia_beta(rec_j+jjB-1, 1:my_B_size_beta, kkk)+ &
                              BI_C_rec_beta(jjB, 1:my_B_size_beta, kkk)
                        END DO
!$OMP                   END PARALLEL DO
                     ENDIF
                  END DO
                  ! DO iiB=1, rec_block_size
                  !   mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1,1:my_B_size,start_point:end_point)=&
                  !   mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1,1:my_B_size,start_point:end_point)+&
                  !   BI_C_rec(iiB,1:my_B_size,start_point:end_point)
                  ! END DO
                  ! DO jjB=1, rec_block_size
                  !   mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1,1:my_B_size,start_point:end_point)=&
                  !   mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1,1:my_B_size,start_point:end_point)+&
                  !   BI_C_rec(jjB+rec_block_size,1:my_B_size,start_point:end_point)
                  ! END DO
               END DO
               CALL timestop(handle2)

               DEALLOCATE (BI_C_rec)
               IF (alpha_beta) DEALLOCATE (BI_C_rec_beta)

            ELSE
               ! we have something to send but nothing to receive

               CALL mp_sendrecv(BI_C_send, proc_send, &
                                null_mat_rec, proc_receive, &
                                para_env_exchange%group)
               IF (alpha_beta) THEN
                  CALL mp_sendrecv(BI_C_send_beta, proc_send, &
                                   null_mat_rec, proc_receive, &
                                   para_env_exchange%group)
               ENDIF

            END IF

            DEALLOCATE (BI_C_send)
            IF (alpha_beta) DEALLOCATE (BI_C_send_beta)
         END DO

      ELSE
         ! noting to send check if we have to receive
         DO proc_shift = 1, para_env_exchange%num_pe-1
            proc_send = proc_map(para_env_exchange%mepos+proc_shift)
            proc_receive = proc_map(para_env_exchange%mepos-proc_shift)
            rec_ij_index = num_ij_pairs(proc_receive)

            IF (ij_index <= rec_ij_index) THEN
               ! we know that proc_receive has something to send for us, let's see what
               ij_counter_rec = &
                  (ij_index-MIN(1, integ_group_pos2color_sub(proc_receive)))*ngroup+integ_group_pos2color_sub(proc_receive)

               rec_i = ij_map(ij_counter_rec, 1)
               rec_j = ij_map(ij_counter_rec, 2)
               rec_block_size = ij_map(ij_counter_rec, 3)

               IF (.NOT. alpha_beta) THEN
                  ALLOCATE (BI_C_rec(2*rec_block_size, my_B_size, my_group_L_size))
               ELSE
                  ALLOCATE (BI_C_rec(rec_block_size, my_B_size, my_group_L_size))
                  ALLOCATE (BI_C_rec_beta(rec_block_size, my_B_size_beta, my_group_L_size))
               ENDIF

               BI_C_rec = 0.0_dp
               IF (alpha_beta) BI_C_rec_beta = 0.0_dp

               CALL mp_sendrecv(null_mat_send, proc_send, &
                                BI_C_rec, proc_receive, &
                                para_env_exchange%group)
               IF (alpha_beta) THEN
                  CALL mp_sendrecv(null_mat_send, proc_send, &
                                   BI_C_rec_beta, proc_receive, &
                                   para_env_exchange%group)
               ENDIF

               CALL timeset(routineN//"_comm2_w", handle2)
               DO irep = 0, num_integ_group-1
                  Lstart_pos = ranges_info_array(1, irep, para_env_exchange%mepos)
                  Lend_pos = ranges_info_array(2, irep, para_env_exchange%mepos)
                  start_point = ranges_info_array(3, irep, para_env_exchange%mepos)
                  end_point = ranges_info_array(4, irep, para_env_exchange%mepos)
                  DO iiB = 1, rec_block_size
!$OMP               PARALLEL DO DEFAULT(NONE) PRIVATE(kkk) &
!$OMP                                         SHARED(start_point,end_point,&
!$OMP                                                mp2_env,rec_i,iiB,my_B_size,BI_C_rec,&
!$OMP                                         alpha_alpha,beta_beta,open_shell)
                     DO kkk = start_point, end_point
                        IF (alpha_alpha .OR. (.NOT. open_shell)) THEN
                           mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1, 1:my_B_size, kkk) = &
                              mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1, 1:my_B_size, kkk)+ &
                              BI_C_rec(iiB, 1:my_B_size, kkk)
                        ENDIF
                        IF (beta_beta) THEN
                           mp2_env%ri_grad%Gamma_P_ia_beta(rec_i+iiB-1, 1:my_B_size, kkk) = &
                              mp2_env%ri_grad%Gamma_P_ia_beta(rec_i+iiB-1, 1:my_B_size, kkk)+ &
                              BI_C_rec(iiB, 1:my_B_size, kkk)
                        ENDIF
                     END DO
!$OMP               END PARALLEL DO
                  END DO
                  DO jjB = 1, rec_block_size
                     IF (.NOT. alpha_beta) THEN
!$OMP                  PARALLEL DO DEFAULT(NONE) PRIVATE(kkk) &
!$OMP                                         SHARED(start_point,end_point,rec_block_size,&
!$OMP                                                mp2_env,rec_j,jjB,my_B_size,BI_C_rec,&
!$OMP                                            alpha_alpha,beta_beta,open_shell)
                        DO kkk = start_point, end_point
                           IF (alpha_alpha .OR. (.NOT. open_shell)) THEN
                              mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1, 1:my_B_size, kkk) = &
                                 mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1, 1:my_B_size, kkk)+ &
                                 BI_C_rec(jjB+rec_block_size, 1:my_B_size, kkk)
                           ENDIF
                           IF (beta_beta) THEN
                              mp2_env%ri_grad%Gamma_P_ia_beta(rec_j+jjB-1, 1:my_B_size, kkk) = &
                                 mp2_env%ri_grad%Gamma_P_ia_beta(rec_j+jjB-1, 1:my_B_size, kkk)+ &
                                 BI_C_rec(jjB+rec_block_size, 1:my_B_size, kkk)
                           ENDIF
                        END DO
!$OMP                  END PARALLEL DO
                     ELSE
!$OMP                  PARALLEL DO DEFAULT(NONE) PRIVATE(kkk) &
!$OMP                                         SHARED(start_point,end_point,&
!$OMP                                       mp2_env,rec_j,jjB,my_B_size_beta,BI_C_rec_beta)
                        DO kkk = start_point, end_point
                           mp2_env%ri_grad%Gamma_P_ia_beta(rec_j+jjB-1, 1:my_B_size_beta, kkk) = &
                              mp2_env%ri_grad%Gamma_P_ia_beta(rec_j+jjB-1, 1:my_B_size_beta, kkk)+ &
                              BI_C_rec_beta(jjB, 1:my_B_size_beta, kkk)
                        END DO
!$OMP                  END PARALLEL DO
                     ENDIF
                  END DO
                  ! DO iiB=1, rec_block_size
                  !   mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1,1:my_B_size,start_point:end_point)=&
                  !   mp2_env%ri_grad%Gamma_P_ia(rec_i+iiB-1,1:my_B_size,start_point:end_point)+&
                  !   BI_C_rec(iiB,1:my_B_size,start_point:end_point)
                  ! END DO
                  ! DO jjB=1, rec_block_size
                  !   mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1,1:my_B_size,start_point:end_point)=&
                  !   mp2_env%ri_grad%Gamma_P_ia(rec_j+jjB-1,1:my_B_size,start_point:end_point)+&
                  !   BI_C_rec(jjB+rec_block_size,1:my_B_size,start_point:end_point)
                  ! END DO
               END DO
               CALL timestop(handle2)

               DEALLOCATE (BI_C_rec)
               IF (alpha_beta) DEALLOCATE (BI_C_rec_beta)

            ELSE
               ! nothing to send nothing to receive
               CALL mp_sendrecv(null_mat_send, proc_send, &
                                null_mat_rec, proc_receive, &
                                para_env_exchange%group)
               IF (alpha_beta) THEN
                  CALL mp_sendrecv(null_mat_send, proc_send, &
                                   null_mat_rec, proc_receive, &
                                   para_env_exchange%group)
               ENDIF

            END IF
         END DO

      END IF
      CALL timestop(handle)

   END SUBROUTINE mp2_redistribute_gamma

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param Eigenval ...
!> \param homo ...
!> \param virtual ...
!> \param open_shell ...
!> \param beta_beta ...
!> \param alpha_beta ...
!> \param Bib_C ...
!> \param unit_nr ...
!> \param dimen_RI ...
!> \param my_B_size ...
!> \param ngroup ...
!> \param num_integ_group ...
!> \param my_group_L_size ...
!> \param color_sub ...
!> \param ranges_info_array ...
!> \param para_env_exchange ...
!> \param para_env_sub ...
!> \param proc_map ...
!> \param my_B_virtual_start ...
!> \param my_B_virtual_end ...
!> \param sizes_array ...
!> \param ends_B_virtual ...
!> \param sizes_B_virtual ...
!> \param starts_B_virtual ...
!> \param sub_proc_map ...
!> \param integ_group_pos2color_sub ...
!> \param local_ab ...
!> \param BIb_C_beta ...
!> \param my_B_size_beta ...
!> \param ends_B_virtual_beta ...
!> \param sizes_B_virtual_beta ...
!> \param starts_B_virtual_beta ...
!> \param my_B_virtual_start_beta ...
!> \param virtual_beta ...
!> \param homo_beta ...
!> \param Eigenval_beta ...
!> \param my_B_virtual_end_beta ...
! **************************************************************************************************
   SUBROUTINE quasi_degenerate_P_ij(mp2_env, Eigenval, homo, virtual, open_shell, &
                                    beta_beta, alpha_beta, Bib_C, unit_nr, dimen_RI, &
                                    my_B_size, ngroup, num_integ_group, my_group_L_size, &
                                    color_sub, ranges_info_array, para_env_exchange, para_env_sub, proc_map, &
                                    my_B_virtual_start, my_B_virtual_end, sizes_array, ends_B_virtual, sizes_B_virtual, &
                                    starts_B_virtual, sub_proc_map, integ_group_pos2color_sub, &
                                    local_ab, BIb_C_beta, my_B_size_beta, &
                                    ends_B_virtual_beta, sizes_B_virtual_beta, starts_B_virtual_beta, my_B_virtual_start_beta, &
                                    virtual_beta, homo_beta, Eigenval_beta, my_B_virtual_end_beta)
      TYPE(mp2_type), POINTER                            :: mp2_env
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: homo, virtual
      LOGICAL                                            :: open_shell, beta_beta, alpha_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C
      INTEGER                                            :: unit_nr, dimen_RI, my_B_size, ngroup, &
                                                            num_integ_group, my_group_L_size, &
                                                            color_sub
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: ranges_info_array
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange, para_env_sub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: proc_map
      INTEGER                                            :: my_B_virtual_start, my_B_virtual_end
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_array, ends_B_virtual, &
                                                            sizes_B_virtual, starts_B_virtual, &
                                                            sub_proc_map, integ_group_pos2color_sub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: local_ab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         OPTIONAL                                        :: BIb_C_beta
      INTEGER, OPTIONAL                                  :: my_B_size_beta
      INTEGER, ALLOCATABLE, DIMENSION(:), OPTIONAL       :: ends_B_virtual_beta, &
                                                            sizes_B_virtual_beta, &
                                                            starts_B_virtual_beta
      INTEGER, OPTIONAL                                  :: my_B_virtual_start_beta, virtual_beta, &
                                                            homo_beta
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: Eigenval_beta
      INTEGER, OPTIONAL                                  :: my_B_virtual_end_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'quasi_degenerate_P_ij', &
         routineP = moduleN//':'//routineN

      INTEGER :: a, a_global, b, b_global, end_point, handle, handle2, ijk_counter, &
         ijk_counter_send, ijk_index, iloops, irep, Lend_pos, Lstart_pos, max_ijk, max_ijk_beta, &
         max_ijk_loop, my_i, my_ijk, my_ijk_beta, my_j, my_k, my_virtual, nloops, proc_receive, &
         proc_send, proc_shift, rec_B_size, rec_B_virtual_end, rec_B_virtual_start, rec_L_size, &
         send_B_size, send_B_virtual_end, send_B_virtual_start, send_i, send_ijk_index, send_j, &
         send_k, size_B_i, size_B_j, size_B_k, start_point
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_ijk, num_ijk_beta
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ijk_map, ijk_map_beta
      REAL(KIND=dp)                                      :: amp_fac, null_mat_rec(2, 2, 2), &
                                                            null_mat_send(2, 2, 2), P_ij_elem
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: external_ab, external_i_aL, t_ab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BI_C_rec, local_i_aL, local_j_aL, &
                                                            local_k_aL

!

      CALL timeset(routineN//"_ij_sing", handle)
! Define the number of loops over orbital triplets

      nloops = 1
      IF (alpha_beta) nloops = 2

      ! Find the number of quasi-degenerate orbitals and orbital triplets

      IF (.NOT. alpha_beta) THEN
         CALL Find_quasi_degenerate_ij(my_ijk, homo, Eigenval, mp2_env, ijk_map, unit_nr, ngroup, &
                                       beta_beta, alpha_beta, para_env_exchange, num_ijk, max_ijk, color_sub)
      ELSE
         CALL Find_quasi_degenerate_ij(my_ijk, homo, Eigenval, mp2_env, ijk_map, unit_nr, ngroup, &
                                       beta_beta, alpha_beta, para_env_exchange, num_ijk, max_ijk, color_sub, &
                                       Eigenval_beta, homo_beta, ijk_map_beta, num_ijk_beta, max_ijk_beta, my_ijk_beta)
      ENDIF

      ! Set amplitude factor
      amp_fac = 2.0_dp
      IF (open_shell) amp_fac = 1.0_dp

      ! Loop(s) over orbital triplets
      DO iloops = 1, nloops
         IF (iloops .EQ. 1) THEN
            size_B_i = my_B_size
            size_B_j = my_B_size
            max_ijk_loop = max_ijk
            my_virtual = virtual
            IF (alpha_beta) THEN
               size_B_k = my_B_size_beta
            ELSE
               size_B_k = my_B_size
            ENDIF
         ELSE
            size_B_i = my_B_size_beta
            size_B_j = my_B_size_beta
            size_B_k = my_B_size
            my_virtual = virtual_beta
            max_ijk_loop = max_ijk_beta
         ENDIF

         ALLOCATE (local_i_aL(dimen_RI, size_B_i, 1))
         ALLOCATE (local_j_aL(dimen_RI, size_B_j, 1))
         ALLOCATE (local_k_aL(dimen_RI, size_B_k, 1))
         ALLOCATE (t_ab(my_virtual, size_B_k))

         DO ijk_index = 1, max_ijk_loop
            IF (iloops .EQ. 2) my_ijk = my_ijk_beta
            IF (ijk_index <= my_ijk) THEN
               ! work to be done
               ijk_counter = (ijk_index-MIN(1, color_sub))*ngroup+color_sub
               IF (iloops .EQ. 1) THEN
                  my_i = ijk_map(ijk_counter, 1)
                  my_j = ijk_map(ijk_counter, 2)
                  my_k = ijk_map(ijk_counter, 3)
               ELSE
                  my_i = ijk_map_beta(ijk_counter, 1)
                  my_j = ijk_map_beta(ijk_counter, 2)
                  my_k = ijk_map_beta(ijk_counter, 3)
               ENDIF

               local_i_aL = 0.0_dp
               local_j_al = 0.0_dp
               local_k_al = 0.0_dp
               DO irep = 0, num_integ_group-1
                  Lstart_pos = ranges_info_array(1, irep, para_env_exchange%mepos)
                  Lend_pos = ranges_info_array(2, irep, para_env_exchange%mepos)
                  start_point = ranges_info_array(3, irep, para_env_exchange%mepos)
                  end_point = ranges_info_array(4, irep, para_env_exchange%mepos)

                  IF (.NOT. alpha_beta) THEN
                     local_i_aL(Lstart_pos:Lend_pos, 1:size_B_i, 1) = BIb_C(start_point:end_point, 1:size_B_i, my_i)
                     local_j_aL(Lstart_pos:Lend_pos, 1:size_B_j, 1) = BIb_C(start_point:end_point, 1:size_B_j, my_j)
                     local_k_aL(Lstart_pos:Lend_pos, 1:size_B_k, 1) = BIb_C(start_point:end_point, 1:size_B_k, my_k)
                  ELSE
                     IF (iloops .EQ. 1) THEN ! For alpha-alpha density
                        local_i_aL(Lstart_pos:Lend_pos, 1:size_B_i, 1) = BIb_C(start_point:end_point, 1:size_B_i, my_i)
                        local_j_aL(Lstart_pos:Lend_pos, 1:size_B_j, 1) = BIb_C(start_point:end_point, 1:size_B_j, my_j)
                        local_k_aL(Lstart_pos:Lend_pos, 1:size_B_k, 1) = BIb_C_beta(start_point:end_point, 1:size_B_k, my_k)
                     ELSE ! For beta-beta density
                        local_i_aL(Lstart_pos:Lend_pos, 1:size_B_i, 1) = BIb_C_beta(start_point:end_point, 1:size_B_i, my_i)
                        local_j_aL(Lstart_pos:Lend_pos, 1:size_B_j, 1) = BIb_C_beta(start_point:end_point, 1:size_B_j, my_j)
                        local_k_aL(Lstart_pos:Lend_pos, 1:size_B_k, 1) = BIb_C(start_point:end_point, 1:size_B_k, my_k)
                     ENDIF
                  ENDIF
               END DO

               DO proc_shift = 1, para_env_exchange%num_pe-1
                  proc_send = proc_map(para_env_exchange%mepos+proc_shift)
                  proc_receive = proc_map(para_env_exchange%mepos-proc_shift)

                  send_ijk_index = num_ijk(proc_send)
                  IF (iloops .EQ. 2) send_ijk_index = num_ijk_beta(proc_send)

                  rec_L_size = sizes_array(proc_receive)
                  ALLOCATE (BI_C_rec(rec_L_size, size_B_i, 1))

                  IF (ijk_index <= send_ijk_index) THEN
                     ! something to send
                     ijk_counter_send = (ijk_index-MIN(1, integ_group_pos2color_sub(proc_send)))* &
                                        ngroup+integ_group_pos2color_sub(proc_send)
                     IF (iloops .EQ. 1) THEN
                        send_i = ijk_map(ijk_counter_send, 1)
                        send_j = ijk_map(ijk_counter_send, 2)
                        send_k = ijk_map(ijk_counter_send, 3)
                     ELSE
                        send_i = ijk_map_beta(ijk_counter_send, 1)
                        send_j = ijk_map_beta(ijk_counter_send, 2)
                        send_k = ijk_map_beta(ijk_counter_send, 3)
                     ENDIF
                  END IF

                  ! occupied i
                  BI_C_rec = 0.0_dp
                  IF (ijk_index <= send_ijk_index) THEN
                     IF (iloops .EQ. 1) THEN
                        CALL mp_sendrecv(BIb_C(1:my_group_L_size, 1:size_B_i, send_i), proc_send, &
                                         BI_C_rec(1:rec_L_size, 1:size_B_i, 1), proc_receive, &
                                         para_env_exchange%group)
                     ELSE
                        CALL mp_sendrecv(BIb_C_beta(1:my_group_L_size, 1:size_B_i, send_i), proc_send, &
                                         BI_C_rec(1:rec_L_size, 1:size_B_i, 1), proc_receive, &
                                         para_env_exchange%group)
                     ENDIF
                  ELSE
                     ! nothing to send
                     CALL mp_sendrecv(null_mat_send, proc_send, &
                                      BI_C_rec(1:rec_L_size, 1:size_B_i, 1:1), proc_receive, &
                                      para_env_exchange%group)
                  END IF
                  DO irep = 0, num_integ_group-1
                     Lstart_pos = ranges_info_array(1, irep, proc_receive)
                     Lend_pos = ranges_info_array(2, irep, proc_receive)
                     start_point = ranges_info_array(3, irep, proc_receive)
                     end_point = ranges_info_array(4, irep, proc_receive)

                     local_i_aL(Lstart_pos:Lend_pos, 1:size_B_i, 1) = BI_C_rec(start_point:end_point, 1:size_B_i, 1)
                  END DO

                  ! occupied j
                  BI_C_rec = 0.0_dp
                  IF (ijk_index <= send_ijk_index) THEN
                     IF (iloops .EQ. 1) THEN
                        CALL mp_sendrecv(BIb_C(1:my_group_L_size, 1:size_B_j, send_j), proc_send, &
                                         BI_C_rec(1:rec_L_size, 1:size_B_j, 1), proc_receive, &
                                         para_env_exchange%group)
                     ELSE ! For beta_beta density, the size is different now
                        !DEALLOCATE(BI_C_rec)
                        !ALLOCATE(BI_C_rec(rec_L_size,size_B_j,1))
                        !BI_C_rec=0.0_dp
                        CALL mp_sendrecv(BIb_C_beta(1:my_group_L_size, 1:size_B_j, send_j), proc_send, &
                                         BI_C_rec(1:rec_L_size, 1:size_B_j, 1), proc_receive, &
                                         para_env_exchange%group)
                     ENDIF
                  ELSE
                     ! nothing to send
                     CALL mp_sendrecv(null_mat_send, proc_send, &
                                      BI_C_rec(1:rec_L_size, 1:size_B_j, 1:1), proc_receive, &
                                      para_env_exchange%group)
                  END IF
                  DO irep = 0, num_integ_group-1
                     Lstart_pos = ranges_info_array(1, irep, proc_receive)
                     Lend_pos = ranges_info_array(2, irep, proc_receive)
                     start_point = ranges_info_array(3, irep, proc_receive)
                     end_point = ranges_info_array(4, irep, proc_receive)

                     local_j_aL(Lstart_pos:Lend_pos, 1:size_B_j, 1) = BI_C_rec(start_point:end_point, 1:size_B_j, 1)
                  END DO

                  ! occupied k
                  BI_C_rec = 0.0_dp
                  DEALLOCATE (BI_C_rec)
                  ALLOCATE (BI_C_rec(rec_L_size, size_B_k, 1))
                  BI_C_rec = 0.0_dp
                  IF (ijk_index <= send_ijk_index) THEN
                     IF (iloops .EQ. 1) THEN
                        IF (.NOT. alpha_beta) THEN
                           CALL mp_sendrecv(BIb_C(1:my_group_L_size, 1:size_B_k, send_k), proc_send, &
                                            BI_C_rec(1:rec_L_size, 1:size_B_k, 1), proc_receive, &
                                            para_env_exchange%group)
                        ELSE
                           CALL mp_sendrecv(BIb_C_beta(1:my_group_L_size, 1:size_B_k, send_k), proc_send, &
                                            BI_C_rec(1:rec_L_size, 1:size_B_k, 1), proc_receive, &
                                            para_env_exchange%group)
                        ENDIF
                     ELSE ! For beta_beta density, the size is different now
                        BI_C_rec = 0.0_dp
                        CALL mp_sendrecv(BIb_C(1:my_group_L_size, 1:size_B_k, send_k), proc_send, &
                                         BI_C_rec(1:rec_L_size, 1:size_B_k, 1), proc_receive, &
                                         para_env_exchange%group)
                     ENDIF
                  ELSE
                     ! nothing to send
                     BI_C_rec = 0.0_dp
                     CALL mp_sendrecv(null_mat_send, proc_send, &
                                      BI_C_rec(1:rec_L_size, 1:size_B_k, 1:1), proc_receive, &
                                      para_env_exchange%group)
                  END IF
                  DO irep = 0, num_integ_group-1
                     Lstart_pos = ranges_info_array(1, irep, proc_receive)
                     Lend_pos = ranges_info_array(2, irep, proc_receive)
                     start_point = ranges_info_array(3, irep, proc_receive)
                     end_point = ranges_info_array(4, irep, proc_receive)

                     local_k_aL(Lstart_pos:Lend_pos, 1:size_B_k, 1) = BI_C_rec(start_point:end_point, 1:size_B_k, 1)
                  END DO

                  DEALLOCATE (BI_C_rec)
               END DO

               ! expand integrals
               CALL timeset(routineN//"_exp_ik", handle2)
               local_ab = 0.0_dp
               IF (iloops .EQ. 2) THEN ! For alpha-beta case for beta-beta density the dimensions are different
                  DEALLOCATE (local_ab)
                  ALLOCATE (local_ab(virtual_beta, size_B_k))
                  local_ab = 0.0_dp
                  CALL dgemm('T', 'N', size_B_i, size_B_k, dimen_RI, 1.0_dp, &
                             local_i_aL(:, :, 1), dimen_RI, local_k_aL(:, :, 1), dimen_RI, &
                             0.0_dp, local_ab(my_B_virtual_start_beta:my_B_virtual_end_beta, 1:size_B_k), size_B_i)
               ELSE
                  CALL dgemm('T', 'N', size_B_i, size_B_k, dimen_RI, 1.0_dp, &
                             local_i_aL(:, :, 1), dimen_RI, local_k_aL(:, :, 1), dimen_RI, &
                             0.0_dp, local_ab(my_B_virtual_start:my_B_virtual_end, 1:size_B_k), size_B_i)
               ENDIF
               DO proc_shift = 1, para_env_sub%num_pe-1
                  proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
                  proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)

                  IF (iloops .EQ. 1) THEN
                     rec_B_size = sizes_B_virtual(proc_receive)
                     rec_B_virtual_end = ends_B_virtual(proc_receive)
                     rec_B_virtual_start = starts_B_virtual(proc_receive)
                  ELSE
                     rec_B_size = sizes_B_virtual_beta(proc_receive)
                     rec_B_virtual_end = ends_B_virtual_beta(proc_receive)
                     rec_B_virtual_start = starts_B_virtual_beta(proc_receive)
                  ENDIF

                  ALLOCATE (external_i_aL(dimen_RI, rec_B_size))
                  external_i_aL = 0.0_dp

                  CALL mp_sendrecv(local_i_aL(:, :, 1), proc_send, &
                                   external_i_aL, proc_receive, &
                                   para_env_sub%group)

                  CALL dgemm('T', 'N', rec_B_size, size_B_k, dimen_RI, 1.0_dp, &
                             external_i_aL, dimen_RI, local_k_aL(:, :, 1), dimen_RI, &
                             0.0_dp, local_ab(rec_B_virtual_start:rec_B_virtual_end, 1:size_B_k), rec_B_size)

                  DEALLOCATE (external_i_aL)
               END DO
               CALL timestop(handle2)

               ! Amplitudes
               CALL timeset(routineN//"_tab", handle2)
               t_ab = 0.0_dp
               ! Alpha-alpha, beta-beta and closed shell
               IF (.NOT. alpha_beta) THEN
                  DO b = 1, size_B_k
                     b_global = b+my_B_virtual_start-1
                     DO a = 1, my_B_size
                        a_global = a+my_B_virtual_start-1
                        t_ab(a_global, b) = (amp_fac*local_ab(a_global, b)-local_ab(b_global, a))/ &
                                            (Eigenval(my_i)+Eigenval(my_k)-Eigenval(homo+a_global)-Eigenval(homo+b_global))
                     END DO
                  END DO
               ELSE
                  IF (iloops .EQ. 1) THEN ! Alpha-beta for alpha-alpha density
                     DO b = 1, size_B_k
                        b_global = b+my_B_virtual_start_beta-1
                        DO a = 1, my_B_size
                           a_global = a+my_B_virtual_start-1
                           t_ab(a_global, b) = local_ab(a_global, b)/ &
                                               (Eigenval(my_i)+Eigenval_beta(my_k)-Eigenval(homo+a_global)- &
                                                Eigenval_beta(homo_beta+b_global))
                        END DO
                     END DO
                  ELSE ! Alpha-beta for beta-beta density
                     DO b = 1, size_B_k
                        b_global = b+my_B_virtual_start-1
                        DO a = 1, my_B_size_beta
                           a_global = a+my_B_virtual_start_beta-1
                           t_ab(a_global, b) = local_ab(a_global, b)/ &
                                               (Eigenval_beta(my_i)+Eigenval(my_k)-Eigenval_beta(homo_beta+a_global)- &
                                                Eigenval(homo+b_global))
                        END DO
                     END DO
                  ENDIF
               ENDIF

               IF (.NOT. alpha_beta) THEN
                  DO proc_shift = 1, para_env_sub%num_pe-1
                     proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
                     proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)
                     rec_B_size = sizes_B_virtual(proc_receive)
                     rec_B_virtual_end = ends_B_virtual(proc_receive)
                     rec_B_virtual_start = starts_B_virtual(proc_receive)
                     send_B_size = sizes_B_virtual(proc_send)
                     send_B_virtual_end = ends_B_virtual(proc_send)
                     send_B_virtual_start = starts_B_virtual(proc_send)

                     ALLOCATE (external_ab(size_B_i, rec_B_size))
                     external_ab = 0.0_dp
                     CALL mp_sendrecv(local_ab(send_B_virtual_start:send_B_virtual_end, 1:size_B_k), proc_send, &
                                      external_ab(1:size_B_i, 1:rec_B_size), proc_receive, para_env_sub%group)

                     DO b = 1, my_B_size
                        b_global = b+my_B_virtual_start-1
                        DO a = 1, rec_B_size
                           a_global = a+rec_B_virtual_start-1
                           t_ab(a_global, b) = (amp_fac*local_ab(a_global, b)-external_ab(b, a))/ &
                                               (Eigenval(my_i)+Eigenval(my_k)-Eigenval(homo+a_global)-Eigenval(homo+b_global))
                        END DO
                     END DO

                     DEALLOCATE (external_ab)
                  END DO
               ENDIF
               CALL timestop(handle2)

               ! Expand the second set of integrals
               CALL timeset(routineN//"_exp_jk", handle2)
               local_ab = 0.0_dp

               IF (iloops .EQ. 2) THEN ! In alpha-beta case for beta-beta density the dimensions are different
                  CALL dgemm('T', 'N', size_B_j, size_B_k, dimen_RI, 1.0_dp, &
                             local_j_aL(:, :, 1), dimen_RI, local_k_aL(:, :, 1), dimen_RI, &
                             0.0_dp, local_ab(my_B_virtual_start_beta:my_B_virtual_end_beta, 1:size_B_k), size_B_j)
               ELSE
                  CALL dgemm('T', 'N', size_B_j, size_B_k, dimen_RI, 1.0_dp, &
                             local_j_aL(:, :, 1), dimen_RI, local_k_aL(:, :, 1), dimen_RI, &
                             0.0_dp, local_ab(my_B_virtual_start:my_B_virtual_end, 1:size_B_k), size_B_j)
               ENDIF

               DO proc_shift = 1, para_env_sub%num_pe-1
                  proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
                  proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)

                  IF (iloops .EQ. 1) THEN
                     rec_B_size = sizes_B_virtual(proc_receive)
                     rec_B_virtual_end = ends_B_virtual(proc_receive)
                     rec_B_virtual_start = starts_B_virtual(proc_receive)
                  ELSE
                     rec_B_size = sizes_B_virtual_beta(proc_receive)
                     rec_B_virtual_end = ends_B_virtual_beta(proc_receive)
                     rec_B_virtual_start = starts_B_virtual_beta(proc_receive)
                  ENDIF

                  ALLOCATE (external_i_aL(dimen_RI, rec_B_size))
                  external_i_aL = 0.0_dp

                  CALL mp_sendrecv(local_j_aL(:, :, 1), proc_send, &
                                   external_i_aL, proc_receive, &
                                   para_env_sub%group)

                  CALL dgemm('T', 'N', rec_B_size, size_B_k, dimen_RI, 1.0_dp, &
                             external_i_aL, dimen_RI, local_k_aL(:, :, 1), dimen_RI, &
                             0.0_dp, local_ab(rec_B_virtual_start:rec_B_virtual_end, 1:size_B_k), rec_B_size)

                  DEALLOCATE (external_i_aL)
               END DO
               CALL timestop(handle2)

               CALL timeset(routineN//"_Pij", handle2)
               ! Alpha-alpha, beta-beta and closed shell
               IF (.NOT. alpha_beta) THEN
                  DO b = 1, size_B_k
                     b_global = b+my_B_virtual_start-1
                     DO a = 1, my_B_size
                        a_global = a+my_B_virtual_start-1
                        local_ab(a_global, b) = local_ab(a_global, b)/ &
                                                (Eigenval(my_j)+Eigenval(my_k)-Eigenval(homo+a_global)-Eigenval(homo+b_global))
                     END DO
                  END DO
               ELSE
                  IF (iloops .EQ. 1) THEN ! Alpha-beta for alpha-alpha density
                     DO b = 1, size_B_k
                        b_global = b+my_B_virtual_start_beta-1
                        DO a = 1, my_B_size
                           a_global = a+my_B_virtual_start-1
                           local_ab(a_global, b) = local_ab(a_global, b)/ &
                                                   (Eigenval(my_j)+Eigenval_beta(my_k)-Eigenval(homo+a_global)- &
                                                    Eigenval_beta(homo_beta+b_global))
                        END DO
                     END DO
                  ELSE ! Alpha-beta for beta-beta density
                     DO b = 1, size_B_k
                        b_global = b+my_B_virtual_start-1
                        DO a = 1, my_B_size_beta
                           a_global = a+my_B_virtual_start_beta-1
                           local_ab(a_global, b) = local_ab(a_global, b)/ &
                                                   (Eigenval_beta(my_j)+Eigenval(my_k)-Eigenval_beta(homo_beta+a_global)- &
                                                    Eigenval(homo+b_global))
                        END DO
                     END DO
                  ENDIF
               ENDIF
               !
               P_ij_elem = SUM(local_ab*t_ab)
               IF ((.NOT. open_shell) .AND. (.NOT. alpha_beta)) THEN
                  P_ij_elem = P_ij_elem*2.0_dp
               ENDIF
               IF ((beta_beta) .OR. (iloops .EQ. 2)) THEN
                  mp2_env%ri_grad%P_ij_beta(my_i, my_j) = mp2_env%ri_grad%P_ij_beta(my_i, my_j)-P_ij_elem
                  mp2_env%ri_grad%P_ij_beta(my_j, my_i) = mp2_env%ri_grad%P_ij_beta(my_j, my_i)-P_ij_elem
               ELSE
                  mp2_env%ri_grad%P_ij(my_i, my_j) = mp2_env%ri_grad%P_ij(my_i, my_j)-P_ij_elem
                  mp2_env%ri_grad%P_ij(my_j, my_i) = mp2_env%ri_grad%P_ij(my_j, my_i)-P_ij_elem
               ENDIF
               CALL timestop(handle2)
            ELSE
               ! no work to be done, possible messeges to be exchanged
               DO proc_shift = 1, para_env_exchange%num_pe-1
                  proc_send = proc_map(para_env_exchange%mepos+proc_shift)
                  proc_receive = proc_map(para_env_exchange%mepos-proc_shift)

                  send_ijk_index = num_ijk(proc_send)
                  IF (iloops .EQ. 2) send_ijk_index = num_ijk_beta(proc_send)

                  IF (ijk_index <= send_ijk_index) THEN
                     ! somethig to send
                     ijk_counter_send = (ijk_index-MIN(1, integ_group_pos2color_sub(proc_send)))*ngroup+ &
                                        integ_group_pos2color_sub(proc_send)
                     IF (iloops .EQ. 1) THEN
                        send_i = ijk_map(ijk_counter_send, 1)
                        send_j = ijk_map(ijk_counter_send, 2)
                        send_k = ijk_map(ijk_counter_send, 3)
                     ELSE
                        send_i = ijk_map_beta(ijk_counter_send, 1)
                        send_j = ijk_map_beta(ijk_counter_send, 2)
                        send_k = ijk_map_beta(ijk_counter_send, 3)
                     ENDIF
                     ! occupied i
                     IF (iloops .EQ. 1) THEN
                        CALL mp_sendrecv(BIb_C(1:my_group_L_size, 1:size_B_i, send_i:send_i), proc_send, &
                                         null_mat_rec, proc_receive, para_env_exchange%group)
                     ELSE
                        CALL mp_sendrecv(BIb_C_beta(1:my_group_L_size, 1:size_B_i, send_i:send_i), proc_send, &
                                         null_mat_rec, proc_receive, para_env_exchange%group)
                     ENDIF
                     ! occupied j
                     IF (iloops .EQ. 1) THEN
                        CALL mp_sendrecv(BIb_C(1:my_group_L_size, 1:size_B_j, send_j:send_j), proc_send, &
                                         null_mat_rec, proc_receive, para_env_exchange%group)
                     ELSE ! For beta_beta density, the size is different now
                        CALL mp_sendrecv(BIb_C_beta(1:my_group_L_size, 1:size_B_j, send_j:send_j), proc_send, &
                                         null_mat_rec, proc_receive, para_env_exchange%group)
                     ENDIF
                     ! occupied k
                     IF (iloops .EQ. 1) THEN
                        IF (.NOT. alpha_beta) THEN
                           CALL mp_sendrecv(BIb_C(1:my_group_L_size, 1:size_B_k, send_k:send_k), proc_send, &
                                            null_mat_rec, proc_receive, para_env_exchange%group)
                        ELSE
                           CALL mp_sendrecv(BIb_C_beta(1:my_group_L_size, 1:size_B_k, send_k:send_k), proc_send, &
                                            null_mat_rec, proc_receive, para_env_exchange%group)
                        ENDIF
                     ELSE ! For beta_beta density, the size is different now
                        CALL mp_sendrecv(BIb_C(1:my_group_L_size, 1:size_B_k, send_k:send_k), proc_send, &
                                         null_mat_rec, proc_receive, para_env_exchange%group)
                     ENDIF

                  ELSE
                     ! nothing to send
                     ! occupied i
                     CALL mp_sendrecv(null_mat_send, proc_send, &
                                      null_mat_rec, proc_receive, &
                                      para_env_exchange%group)
                     ! occupied j
                     CALL mp_sendrecv(null_mat_send, proc_send, &
                                      null_mat_rec, proc_receive, &
                                      para_env_exchange%group)
                     ! occupied k
                     CALL mp_sendrecv(null_mat_send, proc_send, &
                                      null_mat_rec, proc_receive, &
                                      para_env_exchange%group)
                  END IF

               END DO ! proc loop
            END IF
         END DO ! ijk_index loop
         DEALLOCATE (local_i_aL)
         DEALLOCATE (local_j_aL)
         DEALLOCATE (local_k_aL)
         DEALLOCATE (t_ab)
      ENDDO ! over number of loops (iloop)
      !
      DEALLOCATE (ijk_map)
      DEALLOCATE (num_ijk)
      IF (alpha_beta) THEN
         DEALLOCATE (ijk_map_beta, num_ijk_beta)
      ENDIF
      CALL timestop(handle)

   END SUBROUTINE Quasi_degenerate_P_ij

! **************************************************************************************************
!> \brief ...
!> \param my_ijk ...
!> \param homo ...
!> \param Eigenval ...
!> \param mp2_env ...
!> \param ijk_map ...
!> \param unit_nr ...
!> \param ngroup ...
!> \param beta_beta ...
!> \param alpha_beta ...
!> \param para_env_exchange ...
!> \param num_ijk ...
!> \param max_ijk ...
!> \param color_sub ...
!> \param Eigenval_beta ...
!> \param homo_beta ...
!> \param ijk_map_beta ...
!> \param num_ijk_beta ...
!> \param max_ijk_beta ...
!> \param my_ijk_beta ...
! **************************************************************************************************
   SUBROUTINE Find_quasi_degenerate_ij(my_ijk, homo, Eigenval, mp2_env, ijk_map, unit_nr, ngroup, &
                                       beta_beta, alpha_beta, para_env_exchange, num_ijk, max_ijk, color_sub, Eigenval_beta, &
                                       homo_beta, ijk_map_beta, num_ijk_beta, max_ijk_beta, my_ijk_beta)

      INTEGER                                            :: my_ijk, homo
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      TYPE(mp2_type), POINTER                            :: mp2_env
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ijk_map
      INTEGER                                            :: unit_nr, ngroup
      LOGICAL                                            :: beta_beta, alpha_beta
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_ijk
      INTEGER                                            :: max_ijk, color_sub
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: Eigenval_beta
      INTEGER, OPTIONAL                                  :: homo_beta
      INTEGER, ALLOCATABLE, DIMENSION(:, :), OPTIONAL    :: ijk_map_beta
      INTEGER, ALLOCATABLE, DIMENSION(:), OPTIONAL       :: num_ijk_beta
      INTEGER, OPTIONAL                                  :: max_ijk_beta, my_ijk_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'Find_quasi_degenerate_ij', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iib, ijk_counter, jjb, kkb, my_homo, &
                                                            num_sing_ij, total_ijk

      IF (alpha_beta) THEN
         my_homo = homo_beta
      ELSE
         my_homo = homo
      ENDIF

      ! General case
      num_sing_ij = 0
      DO iiB = 1, homo
         ! diagonal elements already updated
         DO jjB = iiB+1, homo
            IF (ABS(Eigenval(jjB)-Eigenval(iiB)) < mp2_env%ri_mp2%eps_canonical) &
               num_sing_ij = num_sing_ij+1
         END DO
      END DO
      IF (.NOT. beta_beta) THEN
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MO_INFO| Number of ij pairs below EPS_CANONICAL:", num_sing_ij
         END IF
      ELSE
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MO_INFO| Number of ij pairs (spin beta) below EPS_CANONICAL:", num_sing_ij
         END IF
      ENDIF
      total_ijk = my_homo*num_sing_ij
      ALLOCATE (ijk_map(total_ijk, 3))
      ijk_map = 0

      my_ijk = 0
      ijk_counter = 0
      DO iiB = 1, homo
         ! diagonal elements already updated
         DO jjB = iiB+1, homo
            IF (ABS(Eigenval(jjB)-Eigenval(iiB)) >= mp2_env%ri_mp2%eps_canonical) CYCLE
            DO kkB = 1, my_homo
               ijk_counter = ijk_counter+1
               ijk_map(ijk_counter, 1) = iiB
               ijk_map(ijk_counter, 2) = jjB
               ijk_map(ijk_counter, 3) = kkB
               IF (MOD(ijk_counter, ngroup) == color_sub) my_ijk = my_ijk+1
            END DO
         END DO
      END DO

      ALLOCATE (num_ijk(0:para_env_exchange%num_pe-1))
      num_ijk = 0
      num_ijk(para_env_exchange%mepos) = my_ijk
      CALL mp_sum(num_ijk, para_env_exchange%group)
      max_ijk = MAXVAL(num_ijk)

      ! Alpha-beta case: we need a second map
      IF (alpha_beta) THEN
         num_sing_ij = 0
         DO iiB = 1, homo_beta
            ! diagonal elements already updated
            DO jjB = iiB+1, homo_beta
               IF (ABS(Eigenval_beta(jjB)-Eigenval_beta(iiB)) < mp2_env%ri_mp2%eps_canonical) &
                  num_sing_ij = num_sing_ij+1
            END DO
         END DO
         ! total number of elemets that have to be computed
         total_ijk = homo*num_sing_ij
         ALLOCATE (ijk_map_beta(total_ijk, 3))
         ijk_map_beta = 0
         my_ijk_beta = 0
         ijk_counter = 0
         DO iiB = 1, homo_beta
            ! diagonal elements already updated
            DO jjB = iiB+1, homo_beta
               IF (ABS(Eigenval_beta(jjB)-Eigenval_beta(iiB)) >= mp2_env%ri_mp2%eps_canonical) CYCLE
               DO kkB = 1, homo
                  ijk_counter = ijk_counter+1
                  ijk_map_beta(ijk_counter, 1) = iiB
                  ijk_map_beta(ijk_counter, 2) = jjB
                  ijk_map_beta(ijk_counter, 3) = kkB
                  IF (MOD(ijk_counter, ngroup) == color_sub) my_ijk_beta = my_ijk_beta+1
               END DO
            END DO
         END DO
         ALLOCATE (num_ijk_beta(0:para_env_exchange%num_pe-1))
         num_ijk_beta = 0
         num_ijk_beta(para_env_exchange%mepos) = my_ijk_beta
         CALL mp_sum(num_ijk_beta, para_env_exchange%group)
         max_ijk_beta = MAXVAL(num_ijk_beta)
      ENDIF

   END SUBROUTINE Find_quasi_degenerate_ij

! **************************************************************************************************
!> \brief compute three center overlap integrals and write them to mat_B_munu
!> \param mat_B_munu ...
!> \param qs_env ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param nmo ...
!> \param dimen_RI ...
!> \param sab_orb_sub ...
!> \param sab_orb_all ...
!> \param reserve_blocks ...
!> \param compute_3c ...
! **************************************************************************************************
   SUBROUTINE compute_3c_overlap_int(mat_B_munu, qs_env, my_group_L_size, my_group_L_start, &
                                     my_group_L_end, nmo, dimen_RI, sab_orb_sub, sab_orb_all, &
                                     reserve_blocks, compute_3c)

      TYPE(cp_dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: mat_B_munu
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: my_group_L_size, my_group_L_start, &
                                                            my_group_L_end, nmo, dimen_RI
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub, sab_orb_all
      LOGICAL                                            :: reserve_blocks, compute_3c

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_3c_overlap_int', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_RI, atom_RI_end, atom_RI_start, blk, blk_cnt, block_end_col, block_end_row, &
         block_start_col, block_start_row, col, handle, handle2, i_bf, iatom, iatom_outer, ikind, &
         iset, isgf_RI, isgfa, jatom, jatom_outer, jkind, jset, jsgfb, kind_RI, LLL, LLL_set_end, &
         LLL_set_start, natom, nblkrows_total, nco_RI, ncoa, ncob, nkind, nset_RI, nseta, nsetb, &
         old_row, row, set_RI, sgf_RI, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_RI_index, kind_of, &
                                                            offset_in_block, old_col_of_LLL, &
                                                            old_row_of_LLL, &
                                                            row_col_from_prim_index, tmp
      INTEGER, DIMENSION(:), POINTER :: col_blk_sizes, l_max_RI, l_min_RI, la_max, la_min, lb_max, &
         lb_min, npgf_RI, npgfa, npgfb, nsgf_RI, nsgfa, nsgfb, row_blk_end, row_blk_offset, &
         row_blk_sizes, row_blk_start
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_RI, first_sgfa, first_sgfb
      LOGICAL                                            :: new_block
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: first_cycle
      REAL(KIND=dp)                                      :: dab, daRI, dbRI
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: block
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: s_abRI, s_abRI_contr, &
                                                            s_abRI_contr_transposed
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rab_outer, raRI, rbRI
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b, set_radius_RI
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_RI, rpgfa, rpgfb, sphi_a, sphi_b, &
                                                            sphi_RI, zet_RI, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_RI_tmp
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b, basis_set_RI
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator, nl_iterator_outer
      TYPE(one_dim_int_array), ALLOCATABLE, DIMENSION(:) :: cols_to_alloc, rows_to_alloc
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      ! either reserve blocks or compute the three-center integrals
      CPASSERT(reserve_blocks .NEQV. compute_3c)

      NULLIFY (qs_kind_set, atomic_kind_set, cell, particle_set, molecule_kind_set, basis_set_RI, basis_set_a, basis_set_b, &
               nl_iterator)

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      cell=cell, particle_set=particle_set, molecule_kind_set=molecule_kind_set, &
                      nkind=nkind, natom=natom)

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

      ! maybe remove row_col_from_prim_index
      ALLOCATE (row_col_from_prim_index(nmo))
      ALLOCATE (offset_in_block(nmo))

      CALL cp_dbcsr_get_info(mat_B_munu(1)%matrix, &
                             row_blk_offset=row_blk_offset, &
                             row_blk_size=row_blk_sizes, &
                             col_blk_size=col_blk_sizes, &
                             nblkrows_total=nblkrows_total)

      old_row = 0

      DO i_bf = 1, nmo

         DO row = 1, nblkrows_total

            IF (i_bf >= row_blk_offset(row) .AND. i_bf < row_blk_offset(row)+row_blk_sizes(row)) THEN

               row_col_from_prim_index(i_bf) = row
               IF (old_row == row) THEN
                  offset_in_block(i_bf) = offset_in_block(i_bf)+1
               ELSE
                  offset_in_block(i_bf) = 1
               END IF

               old_row = row

            END IF

         END DO

      END DO

      IF (reserve_blocks) THEN

         ALLOCATE (rows_to_alloc(my_group_L_size))
         ALLOCATE (cols_to_alloc(my_group_L_size))

         ALLOCATE (first_cycle(my_group_L_size))
         first_cycle = .TRUE.

         ALLOCATE (old_row_of_LLL(my_group_L_size))
         ALLOCATE (old_col_of_LLL(my_group_L_size))

         DO LLL = 1, my_group_L_size

            ALLOCATE (rows_to_alloc(LLL)%array(1))
            rows_to_alloc(LLL)%array(1) = 0
            ALLOCATE (cols_to_alloc(LLL)%array(1))
            cols_to_alloc(LLL)%array(1) = 0

         END DO

      END IF

      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

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
                                   iatom=iatom_outer, jatom=jatom_outer, r=rab_outer)

            IF (iatom_outer .NE. atom_RI .AND. jatom_outer .NE. atom_RI) CYCLE

            CALL neighbor_list_iterator_create(nl_iterator, sab_orb_sub)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)

               CALL get_iterator_info(nl_iterator, &
                                      iatom=iatom, jatom=jatom, r=rab)

               DO set_RI = 1, nset_RI

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

                     raRI = rab_outer
                     rbRI = raRI-rab

                  ELSE IF (atom_RI .EQ. iatom_outer .AND. &
                           iatom .EQ. jatom_outer) THEN

                     raRI = -rab_outer
                     rbRI = raRI-rab

                  ELSE IF (jatom .EQ. iatom_outer .AND. &
                           atom_RI .EQ. jatom_outer) THEN

                     rbRI = rab_outer
                     raRI = rbRI+rab

                  ELSE IF (atom_RI .EQ. iatom_outer .AND. &
                           jatom .EQ. jatom_outer) THEN

                     rbRI = -rab_outer
                     raRI = rbRI+rab

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

                  DO iset = 1, nseta

                     IF (set_radius_a(iset)+set_radius_RI(set_RI) < daRI) CYCLE

                     DO jset = 1, nsetb

                        IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE
                        IF (set_radius_b(jset)+set_radius_RI(set_RI) < dbRI) CYCLE

                        nco_RI = npgf_RI(set_RI)*ncoset(l_max_RI(set_RI))
                        ncoa = npgfa(iset)*ncoset(la_max(iset))
                        ncob = npgfb(jset)*ncoset(lb_max(jset))

                        sgf_RI = first_sgf_RI(1, set_RI)
                        sgfa = first_sgfa(1, iset)
                        sgfb = first_sgfb(1, jset)

                        LLL_set_start = 1+sgf_RI-1-my_group_L_start+1+row_blk_start(atom_RI)-1

                        IF (LLL_set_start > my_group_L_size) CYCLE

                        LLL_set_end = nsgf_RI(set_RI)+sgf_RI-1-my_group_L_start+1+row_blk_start(atom_RI)-1

                        IF (LLL_set_end < 1) CYCLE

                        IF (ncoa*ncob*nco_RI > 0) THEN

                           IF (compute_3c) THEN

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

                                       s_abRI_contr_transposed(jsgfb, isgfa, isgf_RI) = s_abRI_contr(isgfa, jsgfb, isgf_RI)

                                    END DO

                                 END DO

                              END DO

                           END IF

                           DO isgf_RI = 1, nsgf_RI(set_RI)

                              LLL = isgf_RI+sgf_RI-1-my_group_L_start+1+row_blk_start(atom_RI)-1

                              IF (LLL < 1) CYCLE
                              IF (LLL > my_group_L_size) CYCLE

                              IF (iatom < jatom) THEN

                                 row = iatom
                                 col = jatom

                              ELSE IF (iatom >= jatom) THEN

                                 row = jatom
                                 col = iatom

                              END IF

                              IF (reserve_blocks) THEN

                                 IF (first_cycle(LLL)) THEN

                                    rows_to_alloc(LLL)%array(1) = row
                                    cols_to_alloc(LLL)%array(1) = col

                                    first_cycle(LLL) = .FALSE.

                                    old_row_of_LLL(LLL) = row
                                    old_col_of_LLL(LLL) = col

                                 ELSE

                                    IF ((old_row_of_LLL(LLL) .NE. row) .OR. (old_col_of_LLL(LLL) .NE. col)) THEN

                                       blk_cnt = SIZE(rows_to_alloc(LLL)%array)

                                       new_block = .TRUE.

                                       DO blk = 1, blk_cnt

                                          IF (rows_to_alloc(LLL)%array(blk) == row .AND. &
                                              cols_to_alloc(LLL)%array(blk) == col) THEN

                                             new_block = .FALSE.
                                             EXIT

                                          END IF

                                       END DO

                                       IF (new_block) THEN

                                          ALLOCATE (tmp(blk_cnt+1))
                                          tmp(1:blk_cnt) = rows_to_alloc(LLL)%array(1:blk_cnt)
                                          DEALLOCATE (rows_to_alloc(LLL)%array)
                                          ALLOCATE (rows_to_alloc(LLL)%array(blk_cnt+1))
                                          rows_to_alloc(LLL)%array(1:blk_cnt) = tmp(1:blk_cnt)
                                          rows_to_alloc(LLL)%array(blk_cnt+1) = row
                                          tmp(1:blk_cnt) = cols_to_alloc(LLL)%array(1:blk_cnt)
                                          DEALLOCATE (cols_to_alloc(LLL)%array)
                                          ALLOCATE (cols_to_alloc(LLL)%array(blk_cnt+1))
                                          cols_to_alloc(LLL)%array(1:blk_cnt) = tmp(1:blk_cnt)
                                          cols_to_alloc(LLL)%array(blk_cnt+1) = col
                                          DEALLOCATE (tmp)

                                          old_row_of_LLL(LLL) = row
                                          old_col_of_LLL(LLL) = col

                                       END IF

                                    END IF

                                 END IF

                              END IF

                              IF (compute_3c) THEN

                                 ALLOCATE (block(row_blk_sizes(row), col_blk_sizes(col)))
                                 block = 0.0_dp

                                 IF (iatom < jatom) THEN

                                    block_start_row = sgfa
                                    block_end_row = sgfa+nsgfa(iset)-1
                                    block_start_col = sgfb
                                    block_end_col = sgfb+nsgfb(jset)-1

                                    ! factor 0.5 is necessary due to double filling due to double iterate loop
                                    block(block_start_row:block_end_row, &
                                          block_start_col:block_end_col) = 0.5_dp*s_abRI_contr(:, :, isgf_RI)

                                 ELSE IF (iatom > jatom) THEN

                                    block_start_row = sgfb
                                    block_end_row = sgfb+nsgfb(jset)-1
                                    block_start_col = sgfa
                                    block_end_col = sgfa+nsgfa(iset)-1

                                    ! factor 0.5 is necessary due to double filling due to double iterate loop
                                    block(block_start_row:block_end_row, &
                                          block_start_col:block_end_col) = 0.5_dp*s_abRI_contr_transposed(:, :, isgf_RI)

                                 ELSE IF (iatom .EQ. jatom) THEN

                                    block_start_row = sgfa
                                    block_end_row = sgfa+nsgfa(iset)-1
                                    block_start_col = sgfb
                                    block_end_col = sgfb+nsgfb(jset)-1

                                    block(block_start_row:block_end_row, &
                                          block_start_col:block_end_col) = s_abRI_contr(:, :, isgf_RI)

                                 END IF

                                 CALL timeset(routineN//"_put_dbcsr", handle2)

                                 CALL cp_dbcsr_put_block(mat_B_munu(LLL)%matrix, row=row, col=col, block=block, summation=.TRUE.)
                                 DEALLOCATE (block)

                                 CALL timestop(handle2)

                              END IF

                           END DO ! single RI basis functions

                           IF (compute_3c) THEN

                              DEALLOCATE (s_abRI, s_abRI_contr, s_abRI_contr_transposed)

                           END IF

                        END IF ! number of triples > 0

                     END DO ! jset

                  END DO ! iset

               END DO ! set RI

            END DO ! inner neighbor list

            CALL neighbor_list_iterator_release(nl_iterator)

         END DO ! outer neighbor list

         CALL neighbor_list_iterator_release(nl_iterator_outer)

      END DO ! atom_RI

      IF (reserve_blocks) THEN

         CALL timeset(routineN//"_reserve_dbcsr", handle2)

         DO LLL = 1, my_group_L_size

            blk_cnt = SIZE(rows_to_alloc(LLL)%array)

            CALL cp_dbcsr_reserve_blocks(mat_B_munu(LLL)%matrix, &
                                         rows=rows_to_alloc(LLL)%array(1:blk_cnt), &
                                         cols=cols_to_alloc(LLL)%array(1:blk_cnt))

            CALL cp_dbcsr_finalize(mat_B_munu(LLL)%matrix)

            DEALLOCATE (rows_to_alloc(LLL)%array)
            DEALLOCATE (cols_to_alloc(LLL)%array)

         END DO

         CALL timestop(handle2)

      END IF

      DEALLOCATE (atom_from_RI_index, row_col_from_prim_index)

      DEALLOCATE (row_blk_start, row_blk_end)

      IF (reserve_blocks) THEN

         DEALLOCATE (rows_to_alloc, cols_to_alloc, old_row_of_LLL, old_col_of_LLL)

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief compute three center overlap integrals and write them to mat_3_overl_int
!> \param mat_3c_overl_int ...
!> \param mat_munu ...
!> \param qs_env ...
!> \param my_group_L_starts_im_time ...
!> \param my_group_L_ends_im_time ...
!> \param my_group_L_sizes_im_time ...
!> \param nmo ...
!> \param dimen_RI ...
!> \param cut_RI ...
!> \param sab_orb_sub ...
!> \param sab_orb_all ...
! **************************************************************************************************
   SUBROUTINE compute_3c_overlap_int_new(mat_3c_overl_int, mat_munu, qs_env, my_group_L_starts_im_time, &
                                         my_group_L_ends_im_time, my_group_L_sizes_im_time, nmo, dimen_RI, &
                                         cut_RI, sab_orb_sub, sab_orb_all)

      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mat_3c_overl_int
      TYPE(cp_dbcsr_p_type)                              :: mat_munu
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_group_L_starts_im_time, &
                                                            my_group_L_ends_im_time, &
                                                            my_group_L_sizes_im_time
      INTEGER                                            :: nmo, dimen_RI, cut_RI
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub, sab_orb_all

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_3c_overlap_int_new', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_RI, atom_RI_end, atom_RI_start, block_end_col, block_end_row, &
         block_start_col, block_start_row, col, handle, handle2, i_bf, i_cut_RI, iatom, &
         iatom_outer, ikind, iset, isgf_RI, isgfa, jatom, jatom_outer, jkind, jset, jsgfb, &
         kind_RI, LLL, LLL_set_end, LLL_set_start, mepos, my_group_L_end, my_group_L_size, &
         my_group_L_start, natom, nblkrows_total, nco_RI, ncoa, ncob, nkind, nset_RI, nseta, &
         nsetb, nthread, offset_col_due_to_LLL, offset_row_due_to_LLL, old_row, row, set_RI, &
         sgf_RI, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_RI_index, kind_of, &
                                                            offset_in_block, &
                                                            row_col_from_prim_index
      INTEGER, DIMENSION(:), POINTER :: col_blk_sizes, l_max_RI, l_min_RI, la_max, la_min, lb_max, &
         lb_min, npgf_RI, npgfa, npgfb, nsgf_RI, nsgfa, nsgfb, row_blk_end, row_blk_offset, &
         row_blk_sizes, row_blk_start
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_RI, first_sgfa, first_sgfb
      REAL(KIND=dp)                                      :: dab, daRI, dbRI
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: block, block_transposed
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: s_abRI, s_abRI_contr, &
                                                            s_abRI_contr_transposed
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rab_outer, raRI, rbRI
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b, set_radius_RI
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_RI, rpgfa, rpgfb, sphi_a, sphi_b, &
                                                            sphi_RI, zet_RI, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_RI_tmp
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b, basis_set_RI
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator, nl_iterator_outer
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, atomic_kind_set, cell, particle_set, molecule_kind_set, basis_set_RI, basis_set_a, basis_set_b, &
               nl_iterator)

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      cell=cell, particle_set=particle_set, molecule_kind_set=molecule_kind_set, &
                      nkind=nkind, natom=natom)

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

      ! maybe remove row_col_from_prim_index
      ALLOCATE (row_col_from_prim_index(nmo))
      ALLOCATE (offset_in_block(nmo))

      CALL cp_dbcsr_get_info(mat_munu%matrix, &
                             row_blk_offset=row_blk_offset, &
                             row_blk_size=row_blk_sizes, &
                             col_blk_size=col_blk_sizes, &
                             nblkrows_total=nblkrows_total)

      old_row = 0

      DO i_bf = 1, nmo

         DO row = 1, nblkrows_total

            IF (i_bf >= row_blk_offset(row) .AND. i_bf < row_blk_offset(row)+row_blk_sizes(row)) THEN

               row_col_from_prim_index(i_bf) = row
               IF (old_row == row) THEN
                  offset_in_block(i_bf) = offset_in_block(i_bf)+1
               ELSE
                  offset_in_block(i_bf) = 1
               END IF

               old_row = row

            END IF

         END DO

      END DO

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
                                      iatom=iatom_outer, jatom=jatom_outer, r=rab_outer)

               IF (iatom_outer .NE. atom_RI .AND. jatom_outer .NE. atom_RI) CYCLE

               nthread = 1
!$             nthread = omp_get_max_threads()

               CALL neighbor_list_iterator_create(nl_iterator, sab_orb_sub)

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED (nthread,atom_RI,first_sgf_RI,l_max_RI,l_min_RI,npgf_RI,nset_RI,nsgf_RI,rpgf_RI,set_radius_RI,sphi_RI,zet_RI,&
!$OMP           iatom_outer,jatom_outer,rab_outer,kind_of,my_group_L_start,row_blk_start,LLL_set_end,&
!$OMP           mat_3c_overl_int,nl_iterator,qs_kind_set,ncoset,my_group_L_size,row_blk_sizes,&
!$OMP           col_blk_sizes, sab_orb_sub,i_cut_RI) &
!$OMP   PRIVATE (mepos,iatom,jatom,rab,set_RI,raRI,rbRI,dab,daRI,dbRI,ikind,first_sgfa,la_max,la_min,npgfa,nseta,nsgfa,rpgfa,&
!$OMP            set_radius_a,sphi_a,zeta,jkind,first_sgfb,lb_max,lb_min,npgfb,nsetb,nsgfb,rpgfb,set_radius_b,sphi_b,zetb,&
!$OMP            iset,jset,LLL_set_start,nco_RI,ncoa,ncob,sgf_RI,sgfa,sgfb,isgf_RI,isgfa,jsgfb,LLL,row,col,old_row,&
!$OMP            block,block_transposed,&
!$OMP            offset_row_due_to_LLL,offset_col_due_to_LLL,block_start_row,block_end_row,block_start_col,block_end_col,&
!$OMP            handle2,s_abRI,s_abRI_contr,s_abRI_contr_transposed,basis_set_a,basis_set_b)

               mepos = 0
!$             mepos = omp_get_thread_num()

               DO WHILE (neighbor_list_iterate(nl_iterator) == 0)

                  CALL get_iterator_info(nl_iterator, &
                                         iatom=iatom, jatom=jatom, r=rab)

                  DO set_RI = 1, nset_RI

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

                        raRI = rab_outer
                        rbRI = raRI-rab

                     ELSE IF (atom_RI .EQ. iatom_outer .AND. &
                              iatom .EQ. jatom_outer) THEN

                        raRI = -rab_outer
                        rbRI = raRI-rab

                     ELSE IF (jatom .EQ. iatom_outer .AND. &
                              atom_RI .EQ. jatom_outer) THEN

                        rbRI = rab_outer
                        raRI = rbRI+rab

                     ELSE IF (atom_RI .EQ. iatom_outer .AND. &
                              jatom .EQ. jatom_outer) THEN

                        rbRI = -rab_outer
                        raRI = rbRI+rab

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

                     DO iset = 1, nseta

                        IF (set_radius_a(iset)+set_radius_RI(set_RI) < daRI) CYCLE

                        DO jset = 1, nsetb

                           IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE
                           IF (set_radius_b(jset)+set_radius_RI(set_RI) < dbRI) CYCLE

                           nco_RI = npgf_RI(set_RI)*ncoset(l_max_RI(set_RI))
                           ncoa = npgfa(iset)*ncoset(la_max(iset))
                           ncob = npgfb(jset)*ncoset(lb_max(jset))

                           sgf_RI = first_sgf_RI(1, set_RI)
                           sgfa = first_sgfa(1, iset)
                           sgfb = first_sgfb(1, jset)

                           LLL_set_start = 1+sgf_RI-1-my_group_L_start+1+row_blk_start(atom_RI)-1

                           IF (LLL_set_start > my_group_L_size) CYCLE

                           LLL_set_end = nsgf_RI(set_RI)+sgf_RI-1-my_group_L_start+1+row_blk_start(atom_RI)-1

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
                                       s_abRI_contr_transposed(jsgfb, isgfa, isgf_RI) = s_abRI_contr(isgfa, jsgfb, isgf_RI)
                                    END DO
                                 END DO
                              END DO

                              DO isgf_RI = 1, nsgf_RI(set_RI)

                                 LLL = isgf_RI+sgf_RI-1-my_group_L_start+1+row_blk_start(atom_RI)-1

                                 IF (LLL < 1) CYCLE
                                 IF (LLL > my_group_L_size) CYCLE

                                 IF (iatom < jatom) THEN

                                    row = iatom
                                    col = jatom

                                 ELSE IF (iatom >= jatom) THEN

                                    row = jatom
                                    col = iatom

                                 END IF

                                 ALLOCATE (block(row_blk_sizes(row), col_blk_sizes(col)*my_group_L_size))
                                 block = 0.0_dp

                                 ALLOCATE (block_transposed(col_blk_sizes(col), row_blk_sizes(row)*my_group_L_size))
                                 block_transposed = 0.0_dp

                                 offset_row_due_to_LLL = (LLL-1)*row_blk_sizes(row)
                                 offset_col_due_to_LLL = (LLL-1)*col_blk_sizes(col)

                                 IF (iatom < jatom) THEN

                                    block_start_row = sgfa
                                    block_end_row = sgfa+nsgfa(iset)-1
                                    block_start_col = sgfb
                                    block_end_col = sgfb+nsgfb(jset)-1

                                    ! factor 0.5 is necessary due to double filling due to double iterate loop
                                    block(block_start_row:block_end_row, &
                                          block_start_col+offset_col_due_to_LLL:block_end_col+offset_col_due_to_LLL) = &
                                       0.5_dp*s_abRI_contr(:, :, isgf_RI)

                                    block_transposed(block_start_col:block_end_col, &
                                                     block_start_row+offset_row_due_to_LLL:block_end_row+offset_row_due_to_LLL) = &
                                       0.5_dp*s_abRI_contr_transposed(:, :, isgf_RI)

                                 ELSE IF (iatom > jatom) THEN

                                    block_start_row = sgfb
                                    block_end_row = sgfb+nsgfb(jset)-1
                                    block_start_col = sgfa
                                    block_end_col = sgfa+nsgfa(iset)-1

                                    ! factor 0.5 is necessary due to double filling due to double iterate loop
                                    block(block_start_row:block_end_row, &
                                          block_start_col+offset_col_due_to_LLL:block_end_col+offset_col_due_to_LLL) = &
                                       0.5_dp*s_abRI_contr_transposed(:, :, isgf_RI)

                                    block_transposed(block_start_col:block_end_col, &
                                                     block_start_row+offset_row_due_to_LLL:block_end_row+offset_row_due_to_LLL) = &
                                       0.5_dp*s_abRI_contr(:, :, isgf_RI)

                                 ELSE IF (iatom .EQ. jatom) THEN

                                    block_start_row = sgfa
                                    block_end_row = sgfa+nsgfa(iset)-1
                                    block_start_col = sgfb+offset_col_due_to_LLL
                                    block_end_col = sgfb+nsgfb(jset)-1+offset_col_due_to_LLL

                                    block(block_start_row:block_end_row, &
                                          block_start_col:block_end_col) = s_abRI_contr(:, :, isgf_RI)

                                 END IF

!                         CALL cp_dbcsr_get_info(mat_3c_overl_int(i_cut_RI)%matrix, &
!                                                local_rows=local_rows,&
!                                                local_cols=local_cols)
!
!                         is_local_block = .FALSE.
!
!                         num_local_blocks = SIZE(local_rows)
!
!                         DO iblk_local = 1, num_local_blocks
!
!                           IF(local_rows(iblk_local)==row .AND. local_cols(iblk_local)==col) THEN
!
!                             is_local_block = .TRUE.
!                             EXIT
!
!                           END IF
!
!                         END DO

                                 CALL timeset(routineN//"_put_dbcsr", handle2)

                                 CALL cp_dbcsr_put_block(mat_3c_overl_int(i_cut_RI)%matrix, &
                                                         row=row, col=col, block=block, summation=.TRUE.)

                                 IF (row .NE. col) THEN

                                    CALL cp_dbcsr_put_block(mat_3c_overl_int(i_cut_RI)%matrix, &
                                                            row=col, col=row, block=block_transposed, summation=.TRUE.)

                                 END IF

                                 DEALLOCATE (block, block_transposed)

                                 CALL timestop(handle2)

                              END DO ! single RI basis functions

                              DEALLOCATE (s_abRI, s_abRI_contr, s_abRI_contr_transposed)

                           END IF ! number of triples > 0

                        END DO ! jset

                     END DO ! iset

                  END DO ! set RI

               END DO ! inner neighbor list

!$OMP   END PARALLEL

               CALL neighbor_list_iterator_release(nl_iterator)

            END DO ! outer neighbor list

            CALL neighbor_list_iterator_release(nl_iterator_outer)

         END DO ! atom_RI

      END DO ! RI index ranges

      DEALLOCATE (atom_from_RI_index, row_col_from_prim_index)

      DEALLOCATE (row_blk_start, row_blk_end)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief compute three center overlap integrals and write them to mat_3_overl_int
!> \param mat_3c_overl_int ...
!> \param mat_munu ...
!> \param qs_env ...
!> \param nmo ...
!> \param dimen_RI ...
!> \param cut_RI ...
!> \param my_group_L_starts_im_time ...
!> \param my_group_L_ends_im_time ...
!> \param my_group_L_sizes_im_time ...
!> \param sab_orb_sub ...
!> \param sab_orb_all ...
! **************************************************************************************************
   SUBROUTINE reserve_blocks_3c(mat_3c_overl_int, mat_munu, qs_env, nmo, dimen_RI, &
                                cut_RI, my_group_L_starts_im_time, my_group_L_ends_im_time, &
                                my_group_L_sizes_im_time, sab_orb_sub, sab_orb_all)

      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mat_3c_overl_int
      TYPE(cp_dbcsr_p_type)                              :: mat_munu
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: nmo, dimen_RI, cut_RI
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_group_L_starts_im_time, &
                                                            my_group_L_ends_im_time, &
                                                            my_group_L_sizes_im_time
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub, sab_orb_all

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reserve_blocks_3c', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_RI, atom_RI_end, atom_RI_start, blk, blk_cnt, col, handle, handle2, i_bf, &
         i_cut_RI, iatom, iatom_outer, ikind, iset, isgf_RI, jatom, jatom_outer, jkind, jset, &
         kind_RI, LLL, LLL_set_end, LLL_set_start, my_group_L_end, my_group_L_size, &
         my_group_L_start, natom, nblkrows_total, nco_RI, ncoa, ncob, nkind, nset_RI, nseta, &
         nsetb, old_col, old_row, row, set_RI, sgf_RI, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_RI_index, cols_to_alloc, &
                                                            kind_of, offset_in_block, &
                                                            row_col_from_prim_index, &
                                                            rows_to_alloc, tmp
      INTEGER, DIMENSION(:), POINTER :: col_blk_sizes, col_blk_sizes_scaled, l_max_RI, l_min_RI, &
         la_max, la_min, lb_max, lb_min, npgf_RI, npgfa, npgfb, nsgf_RI, nsgfa, nsgfb, &
         row_blk_end, row_blk_offset, row_blk_sizes, row_blk_start
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_RI, first_sgfa, first_sgfb
      LOGICAL                                            :: first_cycle, new_block
      REAL(KIND=dp)                                      :: dab, daRI, dbRI
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rab_outer, raRI, rbRI
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b, set_radius_RI
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_RI, rpgfa, rpgfb, sphi_a, sphi_b, &
                                                            sphi_RI, zet_RI, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_obj)                       :: dist_sub
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_RI_tmp
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b, basis_set_RI
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator, nl_iterator_outer
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, atomic_kind_set, cell, particle_set, molecule_kind_set, basis_set_RI, basis_set_a, basis_set_b, &
               nl_iterator)

      ! allocate matrix for storing 3c overlap integrals
      NULLIFY (mat_3c_overl_int)
      CALL cp_dbcsr_allocate_matrix_set(mat_3c_overl_int, cut_RI)

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      cell=cell, particle_set=particle_set, molecule_kind_set=molecule_kind_set, &
                      nkind=nkind, natom=natom)

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

      ALLOCATE (row_col_from_prim_index(nmo))
      ALLOCATE (offset_in_block(nmo))

      CALL cp_dbcsr_get_info(mat_munu%matrix, &
                             row_blk_offset=row_blk_offset, &
                             row_blk_size=row_blk_sizes, &
                             col_blk_size=col_blk_sizes, &
                             nblkrows_total=nblkrows_total)

      old_row = 0

      DO i_bf = 1, nmo

         DO row = 1, nblkrows_total

            IF (i_bf >= row_blk_offset(row) .AND. i_bf < row_blk_offset(row)+row_blk_sizes(row)) THEN

               row_col_from_prim_index(i_bf) = row
               IF (old_row == row) THEN
                  offset_in_block(i_bf) = offset_in_block(i_bf)+1
               ELSE
                  offset_in_block(i_bf) = 1
               END IF

               old_row = row

            END IF

         END DO

      END DO

      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      DO i_cut_RI = 1, cut_RI

         my_group_L_start = my_group_L_starts_im_time(i_cut_RI)
         my_group_L_end = my_group_L_ends_im_time(i_cut_RI)
         my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

         ALLOCATE (rows_to_alloc(1))
         rows_to_alloc(1) = 0
         ALLOCATE (cols_to_alloc(1))
         cols_to_alloc(1) = 0

         first_cycle = .TRUE.

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
                                      iatom=iatom_outer, jatom=jatom_outer, r=rab_outer)

               IF (iatom_outer .NE. atom_RI .AND. jatom_outer .NE. atom_RI) CYCLE

               CALL neighbor_list_iterator_create(nl_iterator, sab_orb_sub)
               DO WHILE (neighbor_list_iterate(nl_iterator) == 0)

                  CALL get_iterator_info(nl_iterator, &
                                         iatom=iatom, jatom=jatom, r=rab)

                  DO set_RI = 1, nset_RI

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

                        raRI = rab_outer
                        rbRI = raRI-rab

                     ELSE IF (atom_RI .EQ. iatom_outer .AND. &
                              iatom .EQ. jatom_outer) THEN

                        raRI = -rab_outer
                        rbRI = raRI-rab

                     ELSE IF (jatom .EQ. iatom_outer .AND. &
                              atom_RI .EQ. jatom_outer) THEN

                        rbRI = rab_outer
                        raRI = rbRI+rab

                     ELSE IF (atom_RI .EQ. iatom_outer .AND. &
                              jatom .EQ. jatom_outer) THEN

                        rbRI = -rab_outer
                        raRI = rbRI+rab

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

                     DO iset = 1, nseta

                        IF (set_radius_a(iset)+set_radius_RI(set_RI) < daRI) CYCLE

                        DO jset = 1, nsetb

                           IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE
                           IF (set_radius_b(jset)+set_radius_RI(set_RI) < dbRI) CYCLE

                           nco_RI = npgf_RI(set_RI)*ncoset(l_max_RI(set_RI))
                           ncoa = npgfa(iset)*ncoset(la_max(iset))
                           ncob = npgfb(jset)*ncoset(lb_max(jset))

                           sgf_RI = first_sgf_RI(1, set_RI)
                           sgfa = first_sgfa(1, iset)
                           sgfb = first_sgfb(1, jset)

                           LLL_set_start = 1+sgf_RI-1-my_group_L_start+1+row_blk_start(atom_RI)-1

                           IF (LLL_set_start > my_group_L_size) CYCLE

                           LLL_set_end = nsgf_RI(set_RI)+sgf_RI-1-my_group_L_start+1+row_blk_start(atom_RI)-1

                           IF (LLL_set_end < 1) CYCLE

                           IF (ncoa*ncob*nco_RI > 0) THEN

                              DO isgf_RI = 1, nsgf_RI(set_RI)

                                 LLL = isgf_RI+sgf_RI-1-my_group_L_start+1+row_blk_start(atom_RI)-1

                                 IF (LLL < 1) CYCLE
                                 IF (LLL > my_group_L_size) CYCLE

                                 IF (iatom < jatom) THEN

                                    row = iatom
                                    col = jatom

                                 ELSE IF (iatom >= jatom) THEN

                                    row = jatom
                                    col = iatom

                                 END IF

                                 IF (first_cycle) THEN

                                    IF (row == col) THEN

                                       rows_to_alloc = row
                                       cols_to_alloc = col

                                    ELSE

                                       DEALLOCATE (rows_to_alloc)
                                       DEALLOCATE (cols_to_alloc)

                                       ALLOCATE (rows_to_alloc(2))
                                       ALLOCATE (cols_to_alloc(2))

                                       rows_to_alloc(1) = row
                                       rows_to_alloc(2) = col

                                       cols_to_alloc(1) = col
                                       cols_to_alloc(2) = row

                                    END IF

                                    old_row = row
                                    old_col = col

                                    first_cycle = .FALSE.

                                 ELSE

                                    new_block = .FALSE.

                                    IF ((old_row .NE. row) .OR. (old_col .NE. col)) THEN

                                       blk_cnt = SIZE(rows_to_alloc)

                                       new_block = .TRUE.

                                       DO blk = 1, blk_cnt

                                          IF (rows_to_alloc(blk) == row .AND. &
                                              cols_to_alloc(blk) == col) THEN

                                             new_block = .FALSE.
                                             EXIT

                                          END IF

                                       END DO

                                       IF (new_block .AND. row == col) THEN

                                          ALLOCATE (tmp(blk_cnt))
                                          tmp(1:blk_cnt) = rows_to_alloc(1:blk_cnt)
                                          DEALLOCATE (rows_to_alloc)
                                          ALLOCATE (rows_to_alloc(blk_cnt+1))
                                          rows_to_alloc(1:blk_cnt) = tmp(1:blk_cnt)
                                          rows_to_alloc(blk_cnt+1) = row
                                          tmp(1:blk_cnt) = cols_to_alloc(1:blk_cnt)
                                          DEALLOCATE (cols_to_alloc)
                                          ALLOCATE (cols_to_alloc(blk_cnt+1))
                                          cols_to_alloc(1:blk_cnt) = tmp(1:blk_cnt)
                                          cols_to_alloc(blk_cnt+1) = col
                                          DEALLOCATE (tmp)

                                          old_row = row
                                          old_col = col

                                       ELSE IF (new_block .AND. row .NE. col) THEN

                                          ALLOCATE (tmp(blk_cnt))
                                          tmp(1:blk_cnt) = rows_to_alloc(1:blk_cnt)
                                          DEALLOCATE (rows_to_alloc)
                                          ALLOCATE (rows_to_alloc(blk_cnt+2))
                                          rows_to_alloc(1:blk_cnt) = tmp(1:blk_cnt)
                                          rows_to_alloc(blk_cnt+1) = row
                                          rows_to_alloc(blk_cnt+2) = col
                                          tmp(1:blk_cnt) = cols_to_alloc(1:blk_cnt)
                                          DEALLOCATE (cols_to_alloc)
                                          ALLOCATE (cols_to_alloc(blk_cnt+2))
                                          cols_to_alloc(1:blk_cnt) = tmp(1:blk_cnt)
                                          cols_to_alloc(blk_cnt+1) = col
                                          cols_to_alloc(blk_cnt+2) = row
                                          DEALLOCATE (tmp)

                                          old_row = row
                                          old_col = col

                                       END IF

                                    END IF

                                 END IF

                              END DO ! single RI basis functions

                           END IF ! number of triples > 0

                        END DO ! jset

                     END DO ! iset

                  END DO ! set RI

               END DO ! inner neighbor list

               CALL neighbor_list_iterator_release(nl_iterator)

            END DO ! outer neighbor list

            CALL neighbor_list_iterator_release(nl_iterator_outer)

         END DO ! atom_RI

         CALL timeset(routineN//"_reserve_dbcsr", handle2)

         blk_cnt = SIZE(rows_to_alloc)

         CALL cp_dbcsr_get_info(mat_munu%matrix, &
                                row_blk_size=row_blk_sizes, &
                                col_blk_size=col_blk_sizes, &
                                distribution=dist_sub)

         ALLOCATE (col_blk_sizes_scaled(SIZE(col_blk_sizes)))
         col_blk_sizes_scaled(:) = col_blk_sizes(:)*my_group_L_size

         ALLOCATE (mat_3c_overl_int(i_cut_RI)%matrix)
         CALL cp_dbcsr_init(mat_3c_overl_int(i_cut_RI)%matrix)

         CALL cp_dbcsr_create(matrix=mat_3c_overl_int(i_cut_RI)%matrix, &
                              name="(munuP)", &
                              dist=dist_sub, &
                              matrix_type=dbcsr_type_no_symmetry, &
                              row_blk_size=row_blk_sizes, &
                              col_blk_size=col_blk_sizes_scaled, &
                              nze=0)

         IF (rows_to_alloc(1) .NE. 0) THEN

            CALL cp_dbcsr_reserve_blocks(mat_3c_overl_int(i_cut_RI)%matrix, &
                                         rows=rows_to_alloc(1:blk_cnt), &
                                         cols=cols_to_alloc(1:blk_cnt))

         END IF

         CALL cp_dbcsr_finalize(mat_3c_overl_int(i_cut_RI)%matrix)

         DEALLOCATE (col_blk_sizes_scaled)
         DEALLOCATE (rows_to_alloc)
         DEALLOCATE (cols_to_alloc)

         CALL timestop(handle2)

      END DO

      DEALLOCATE (atom_from_RI_index, row_col_from_prim_index)

      DEALLOCATE (row_blk_start, row_blk_end)

      CALL timestop(handle)

   END SUBROUTINE

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

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_group_L_starts_im_time, &
                                                            my_group_L_ends_im_time, &
                                                            my_group_L_sizes_im_time
      INTEGER                                            :: dimen_RI, ngroup, cut_memory, cut_RI, &
                                                            color_sub
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_group_L_im_time', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: cut_RI_old, handle, iblk_RI, icut, &
                                                            itmp(2), natom, nblks_RI, nkind, &
                                                            size_RI, start_RI
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ends_array_mem_cut, my_group_L_ends_im_time_tmp, &
         my_group_L_starts_im_time_tmp, row_blk_end_RI, row_blk_offset_RI, sizes_array_mem_cut, &
         starts_array_mem_cut
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes_RI
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_ri_aux
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      ALLOCATE (sizes_array_mem_cut(0:cut_memory-1))
      sizes_array_mem_cut = 0
      ALLOCATE (starts_array_mem_cut(0:cut_memory-1))
      starts_array_mem_cut = 0
      ALLOCATE (ends_array_mem_cut(0:cut_memory-1))
      ends_array_mem_cut = 0

      DO icut = 0, cut_memory-1

         itmp = get_limit(dimen_RI, cut_memory, icut)
         starts_array_mem_cut(icut) = itmp(1)
         ends_array_mem_cut(icut) = itmp(2)
         sizes_array_mem_cut(icut) = itmp(2)-itmp(1)+1

      END DO

      ALLOCATE (my_group_L_starts_im_time(1:cut_memory))
      my_group_L_starts_im_time = 0
      ALLOCATE (my_group_L_ends_im_time(1:cut_memory))
      my_group_L_ends_im_time = 0

      DO icut = 0, cut_memory-1

         size_RI = sizes_array_mem_cut(icut)
         start_RI = starts_array_mem_cut(icut)

         itmp = get_limit(size_RI, ngroup, color_sub)

         my_group_L_starts_im_time(icut+1) = itmp(1)+start_RI-1
         my_group_L_ends_im_time(icut+1) = itmp(2)+start_RI-1

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
         row_blk_offset_RI(iblk_RI) = row_blk_offset_RI(iblk_RI-1)+row_blk_sizes_RI(iblk_RI-1)
      END DO

      row_blk_end_RI(nblks_RI) = dimen_RI
      DO iblk_RI = 1, nblks_RI-1
         row_blk_end_RI(iblk_RI) = row_blk_offset_RI(iblk_RI+1)-1
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

                  cut_RI = cut_RI+1

                  ALLOCATE (my_group_L_starts_im_time_tmp(cut_RI-1))
                  ALLOCATE (my_group_L_ends_im_time_tmp(cut_RI-1))

                  my_group_L_starts_im_time_tmp(:) = my_group_L_starts_im_time(:)
                  my_group_L_ends_im_time_tmp(:) = my_group_L_ends_im_time(:)

                  DEALLOCATE (my_group_L_starts_im_time)
                  DEALLOCATE (my_group_L_ends_im_time)

                  ALLOCATE (my_group_L_starts_im_time(cut_RI))
                  ALLOCATE (my_group_L_ends_im_time(cut_RI))

                  my_group_L_starts_im_time(1:icut) = my_group_L_starts_im_time_tmp(1:icut)
                  my_group_L_starts_im_time(icut+1) = row_blk_offset_RI(iblk_RI+1)
                  IF (cut_RI >= icut+2) THEN
                     my_group_L_starts_im_time(icut+2:cut_RI) = my_group_L_starts_im_time_tmp(icut+1:cut_RI-1)
                  END IF

                  IF (icut-1 >= 1) THEN
                     my_group_L_ends_im_time(1:icut-1) = my_group_L_ends_im_time_tmp(1:icut-1)
                  END IF
                  my_group_L_ends_im_time(icut) = row_blk_offset_RI(iblk_RI+1)-1
                  my_group_L_ends_im_time(icut+1:cut_RI) = my_group_L_ends_im_time_tmp(icut:cut_RI-1)

                  DEALLOCATE (my_group_L_starts_im_time_tmp)
                  DEALLOCATE (my_group_L_ends_im_time_tmp)

               END IF

            END DO

         END DO

      END DO

      ALLOCATE (my_group_L_sizes_im_time(cut_RI))
      my_group_L_sizes_im_time(:) = my_group_L_ends_im_time(:)-my_group_L_starts_im_time(:)+1

      ALLOCATE (qs_env%mp2_env%ri_rpa_im_time_util(cut_memory))

      qs_env%mp2_env%ri_rpa_im_time_util(1)%cut_RI = cut_RI

      ALLOCATE (qs_env%mp2_env%ri_rpa_im_time_util(1)%my_group_L_starts_im_time(1:cut_RI))
      qs_env%mp2_env%ri_rpa_im_time_util(1)%my_group_L_starts_im_time(:) = my_group_L_starts_im_time(:)
      ALLOCATE (qs_env%mp2_env%ri_rpa_im_time_util(1)%my_group_L_ends_im_time(1:cut_RI))
      qs_env%mp2_env%ri_rpa_im_time_util(1)%my_group_L_ends_im_time(:) = my_group_L_ends_im_time(:)
      ALLOCATE (qs_env%mp2_env%ri_rpa_im_time_util(1)%my_group_L_sizes_im_time(1:cut_RI))
      qs_env%mp2_env%ri_rpa_im_time_util(1)%my_group_L_sizes_im_time(:) = my_group_L_sizes_im_time(:)

      DEALLOCATE (sizes_array_mem_cut, starts_array_mem_cut, ends_array_mem_cut)
      DEALLOCATE (row_blk_sizes_RI)

      CALL timestop(handle)

   END SUBROUTINE

END MODULE mp2_ri_gpw
