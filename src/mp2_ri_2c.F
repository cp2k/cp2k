!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Framework for 2c-integrals for RI
!> \par History
!>      06.2012 created [Mauro Del Ben]
!>      03.2019 separated from mp2_ri_gpw [Frederick Stein]
! **************************************************************************************************
MODULE mp2_ri_2c
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_decompose,&
                                              cp_cfm_gemm,&
                                              cp_cfm_scale_and_add_fm,&
                                              cp_cfm_triangular_invert
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_dist2d_to_dist,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_eri_mme_interface,            ONLY: cp_eri_mme_param
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_triangular_invert
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                              cp_fm_syevx
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_copy_general, cp_fm_create, cp_fm_get_info, cp_fm_indxg2p, cp_fm_p_type, &
        cp_fm_release, cp_fm_set_all, cp_fm_to_fm, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_p_type, dbcsr_reserve_all_blocks, dbcsr_set, dbcsr_type, &
        dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE group_dist_types,                ONLY: create_group_dist,&
                                              get_group_dist,&
                                              group_dist_d1_type,&
                                              release_group_dist
   USE input_constants,                 ONLY: do_eri_gpw,&
                                              do_eri_mme,&
                                              do_eri_os,&
                                              do_potential_coulomb,&
                                              do_potential_id,&
                                              do_potential_long
   USE kinds,                           ONLY: dp
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE kpoint_methods,                  ONLY: rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_bcast,&
                                              mp_comm_split_direct,&
                                              mp_sendrecv,&
                                              mp_sum
   USE mp2_eri,                         ONLY: mp2_eri_2c_integrate
   USE mp2_eri_gpw,                     ONLY: mp2_eri_2c_integrate_gpw
   USE mp2_types,                       ONLY: integ_mat_buffer_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: set_ks_env
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_tensors,                      ONLY: build_2c_integrals,&
                                              build_2c_neighbor_lists
   USE rpa_communication,               ONLY: communicate_buffer

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_ri_2c'

   PUBLIC :: get_2c_integrals

CONTAINS

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
!> \param dimen_RI_red reduced RI dimension,  needed if we perform SVD instead of Cholesky
!> \param kpoints ...
!> \param mo_coeff ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param gd_array ...
!> \param calc_PQ_cond_num ...
!> \param cond_num ...
!> \param do_svd ...
!> \param eps_svd ...
!> \param num_small_eigen ...
!> \param potential ...
!> \param ri_metric ...
!> \param fm_matrix_L_RI_metric ...
!> \param do_im_time ...
!> \param do_kpoints ...
!> \param mp2_eps_pgf_orb_S ...
!> \param qs_kind_set ...
!> \param sab_orb_sub ...
! **************************************************************************************************
   SUBROUTINE get_2c_integrals(qs_env, eri_method, eri_param, para_env, para_env_sub, para_env_L, mp2_memory, &
                               fm_matrix_L, ngroup, color_sub, dimen_RI, dimen_RI_red, kpoints, mo_coeff, &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               gd_array, calc_PQ_cond_num, cond_num, do_svd, eps_svd, &
                               num_small_eigen, potential, ri_metric, &
                               fm_matrix_L_RI_metric, &
                               do_im_time, do_kpoints, mp2_eps_pgf_orb_S, qs_kind_set, sab_orb_sub)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: eri_method
      TYPE(cp_eri_mme_param), POINTER                    :: eri_param
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub, para_env_L
      REAL(KIND=dp), INTENT(IN)                          :: mp2_memory
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L
      INTEGER, INTENT(IN)                                :: ngroup, color_sub
      INTEGER, INTENT(OUT)                               :: dimen_RI, dimen_RI_red
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      INTEGER, INTENT(OUT)                               :: my_group_L_size, my_group_L_start, &
                                                            my_group_L_end
      TYPE(group_dist_d1_type), INTENT(OUT)              :: gd_array
      LOGICAL, INTENT(IN)                                :: calc_PQ_cond_num
      REAL(KIND=dp), INTENT(OUT)                         :: cond_num
      LOGICAL, INTENT(IN)                                :: do_svd
      REAL(KIND=dp), INTENT(IN)                          :: eps_svd
      INTEGER, INTENT(OUT)                               :: num_small_eigen
      TYPE(libint_potential_type)                        :: potential, ri_metric
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_RI_metric
      LOGICAL, INTENT(IN)                                :: do_im_time, do_kpoints
      REAL(KIND=dp), INTENT(IN)                          :: mp2_eps_pgf_orb_S
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_2c_integrals'

      INTEGER                                            :: handle, i_real_imag, i_size, j_size
      REAL(KIND=dp)                                      :: eps_pgf_orb_old
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_kpoints
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_S_inv_work
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      ! calculate V and store it in fm_matrix_L
      CALL compute_2c_integrals(qs_env, eri_method, eri_param, para_env, para_env_sub, para_env_L, mp2_memory, &
                                fm_matrix_L, ngroup, color_sub, dimen_RI, mo_coeff, &
                                my_group_L_size, my_group_L_start, my_group_L_end, &
                                gd_array, calc_PQ_cond_num, cond_num, &
                                num_small_eigen, potential, sab_orb_sub, do_im_time=do_im_time)

      dimen_RI_red = dimen_RI

      IF (ri_metric%potential_type == potential%potential_type .AND. .NOT. do_im_time) THEN
         CALL decomp_mat_L(fm_matrix_L, do_svd, eps_svd, num_small_eigen, cond_num, .TRUE., gd_array, ngroup, &
                           dimen_RI, dimen_RI_red, para_env)
      ELSE

         ! Obara-Saika overlap matrix
         IF (do_im_time) THEN
            CALL get_qs_env(qs_env, dft_control=dft_control)

            ! re-init the radii to be able to generate pair lists with appropriate screening for overlap matrix
            eps_pgf_orb_old = dft_control%qs_control%eps_pgf_orb
            dft_control%qs_control%eps_pgf_orb = mp2_eps_pgf_orb_S
            CALL init_interaction_radii(dft_control%qs_control, qs_kind_set)

            CALL RI_2c_integral_mat(qs_env, fm_matrix_L_RI_metric, fm_matrix_L, dimen_RI, ri_metric, &
                                    do_kpoints, kpoints)

            ! re-init the radii to previous values
            dft_control%qs_control%eps_pgf_orb = eps_pgf_orb_old
            CALL init_interaction_radii(dft_control%qs_control, qs_kind_set)

            ! GPW overlap matrix
         ELSE

            CALL cp_para_env_release(para_env_L)
            CALL release_group_dist(gd_array)

            ALLOCATE (fm_matrix_L_RI_metric(1, 1))

            ! Calculate matrix of RI operator (for overlap metric: S), store it in fm_matrix_L_RI_metric
            CALL compute_2c_integrals(qs_env, eri_method, eri_param, para_env, para_env_sub, para_env_L, mp2_memory, &
                                      fm_matrix_L_RI_metric(1, 1)%matrix, ngroup, color_sub, dimen_RI, mo_coeff, &
                                      my_group_L_size, my_group_L_start, my_group_L_end, &
                                      gd_array, calc_PQ_cond_num, cond_num, &
                                      num_small_eigen, ri_metric, sab_orb_sub, &
                                      fm_matrix_L_extern=fm_matrix_L)

         END IF

         IF (do_kpoints) THEN

            CALL compute_V_by_lattice_sum(qs_env, fm_matrix_L_kpoints, fm_matrix_L_RI_metric, kpoints)

            CALL inversion_of_S_and_mult_with_chol_dec_of_V(fm_matrix_L_RI_metric, fm_matrix_L_kpoints, &
                                                            dimen_RI, kpoints)

         ELSE

            CALL decomp_mat_L(fm_matrix_L, do_svd, eps_svd, num_small_eigen, cond_num, .FALSE., gd_array, ngroup, &
                              dimen_RI, dimen_RI_red, para_env)

            CALL decomp_mat_L(fm_matrix_L_RI_metric(1, 1)%matrix, .FALSE., 0.0_dp, num_small_eigen, cond_num, .TRUE., &
                              gd_array, ngroup, dimen_RI, dimen_RI_red, para_env)

            NULLIFY (fm_matrix_S_inv_work)
            CALL cp_fm_create(fm_matrix_S_inv_work, fm_matrix_L_RI_metric(1, 1)%matrix%matrix_struct)
            CALL cp_fm_set_all(fm_matrix_S_inv_work, 0.0_dp)

            CALL cp_gemm('N', 'T', dimen_RI, dimen_RI, dimen_RI, 1.0_dp, fm_matrix_L_RI_metric(1, 1)%matrix, &
                         fm_matrix_L_RI_metric(1, 1)%matrix, 0.0_dp, fm_matrix_S_inv_work)

            IF (do_svd) THEN
               ! We have to reset the size of fm_matrix_L_RI_metric
               CALL reset_size_matrix(fm_matrix_L_RI_metric(1, 1)%matrix, dimen_RI_red, fm_matrix_L%matrix_struct)

               ! L (=fm_matrix_L_RI_metric) = S^(-1)*chol_dec(V)
               CALL cp_gemm('T', 'N', dimen_RI, dimen_RI_red, dimen_RI, 1.0_dp, fm_matrix_S_inv_work, fm_matrix_L, &
                            0.0_dp, fm_matrix_L_RI_metric(1, 1)%matrix)
            ELSE
               ! L (=fm_matrix_L_RI_metric) = S^(-1)*chol_dec(V)
               CALL cp_gemm('T', 'T', dimen_RI, dimen_RI, dimen_RI, 1.0_dp, fm_matrix_S_inv_work, fm_matrix_L, &
                            0.0_dp, fm_matrix_L_RI_metric(1, 1)%matrix)
            END IF

            ! clean the S_inv matrix
            CALL cp_fm_release(fm_matrix_S_inv_work)
         END IF

         IF (.NOT. do_im_time) THEN

            CALL cp_fm_to_fm(fm_matrix_L_RI_metric(1, 1)%matrix, fm_matrix_L)
            DO j_size = 1, SIZE(fm_matrix_L_RI_metric, 1)
            DO i_size = 1, SIZE(fm_matrix_L_RI_metric, 2)
               CALL cp_fm_release(fm_matrix_L_RI_metric(i_size, j_size)%matrix)
            END DO
            END DO
            DEALLOCATE (fm_matrix_L_RI_metric)

         END IF

      END IF

      IF (do_kpoints) THEN

         DO i_size = 1, SIZE(fm_matrix_L_kpoints, 1)
         DO i_real_imag = 1, SIZE(fm_matrix_L_kpoints, 2)
            CALL cp_fm_release(fm_matrix_L_kpoints(i_size, i_real_imag)%matrix)
         END DO
         END DO

         DEALLOCATE (fm_matrix_L_kpoints)

      END IF

      ! If the group distribution changed because of an SVD, we get the new values here
      CALL get_group_dist(gd_array, color_sub, my_group_L_start, my_group_L_end, my_group_L_size)

      CALL timestop(handle)

   END SUBROUTINE get_2c_integrals

! **************************************************************************************************
!> \brief ...
!> \param fm_matrix_L ...
!> \param do_svd ...
!> \param eps_svd ...
!> \param num_small_eigen ...
!> \param cond_num ...
!> \param do_inversion ...
!> \param gd_array ...
!> \param ngroup ...
!> \param dimen_RI ...
!> \param dimen_RI_red ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE decomp_mat_L(fm_matrix_L, do_svd, eps_svd, num_small_eigen, cond_num, do_inversion, gd_array, ngroup, &
                           dimen_RI, dimen_RI_red, para_env)

      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L
      LOGICAL, INTENT(IN)                                :: do_svd
      REAL(KIND=dp), INTENT(IN)                          :: eps_svd
      INTEGER, INTENT(INOUT)                             :: num_small_eigen
      REAL(KIND=dp), INTENT(INOUT)                       :: cond_num
      LOGICAL, INTENT(IN)                                :: do_inversion
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: gd_array
      INTEGER, INTENT(IN)                                :: ngroup, dimen_RI
      INTEGER, INTENT(INOUT)                             :: dimen_RI_red
      TYPE(cp_para_env_type), POINTER                    :: para_env

      IF (do_svd) THEN
         CALL matrix_root_with_svd(fm_matrix_L, eps_svd, num_small_eigen, cond_num, do_inversion, para_env)

         dimen_RI_red = dimen_RI - num_small_eigen

         ! We changed the size of fm_matrix_L in matrix_root_with_svd.
         ! So, we have to get new group sizes
         CALL release_group_dist(gd_array)

         CALL create_group_dist(gd_array, ngroup, dimen_RI_red)
      ELSE

         ! calculate Cholesky decomposition L (V=LL^T)
         CALL cholesky_decomp(fm_matrix_L, dimen_RI, do_inversion=do_inversion)

         IF (do_inversion) CALL invert_L(fm_matrix_L)
      END IF

   END SUBROUTINE decomp_mat_L

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param fm_matrix_L_kpoints ...
!> \param fm_matrix_L_RI_metric ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE compute_V_by_lattice_sum(qs_env, fm_matrix_L_kpoints, fm_matrix_L_RI_metric, kpoints)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_kpoints, &
                                                            fm_matrix_L_RI_metric
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_V_by_lattice_sum'

      INTEGER                                            :: handle, i_dim, i_real_imag, ikp, nkp
      INTEGER, DIMENSION(3)                              :: periodic
      REAL(dp)                                           :: cutoff
      REAL(KIND=dp), DIMENSION(3)                        :: abc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_RI_aux_transl, matrix_v_RI_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (matrix_s_RI_aux_transl, particle_set, cell, qs_kind_set)

      CALL get_qs_env(qs_env=qs_env, &
                      matrix_s_RI_aux_kp=matrix_s_RI_aux_transl, &
                      particle_set=particle_set, &
                      cell=cell, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set)

      ! check that we have a 2n x 2m x 2k mesh to guarantee good convergence
      CALL get_cell(cell=cell, periodic=periodic)
      DO i_dim = 1, 3
         IF (periodic(i_dim) == 1) THEN
            CPASSERT(MODULO(kpoints%nkp_grid(i_dim), 2) == 0)
         END IF
      END DO

      nkp = kpoints%nkp

      NULLIFY (fm_matrix_L_kpoints)
      ALLOCATE (fm_matrix_L_kpoints(nkp, 2))
      DO ikp = 1, nkp
         DO i_real_imag = 1, 2
            CALL cp_fm_create(fm_matrix_L_kpoints(ikp, i_real_imag)%matrix, &
                              fm_matrix_L_RI_metric(1, i_real_imag)%matrix%matrix_struct)
            CALL cp_fm_set_all(fm_matrix_L_kpoints(ikp, i_real_imag)%matrix, 0.0_dp)
         END DO
      END DO

      NULLIFY (matrix_v_RI_kp)
      CALL dbcsr_allocate_matrix_set(matrix_v_RI_kp, nkp, 2)

      DO ikp = 1, nkp

         ALLOCATE (matrix_v_RI_kp(ikp, 1)%matrix)
         CALL dbcsr_create(matrix_v_RI_kp(ikp, 1)%matrix, template=matrix_s_RI_aux_transl(1, 1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_reserve_all_blocks(matrix_v_RI_kp(ikp, 1)%matrix)
         CALL dbcsr_set(matrix_v_RI_kp(ikp, 1)%matrix, 0.0_dp)

         ALLOCATE (matrix_v_RI_kp(ikp, 2)%matrix)
         CALL dbcsr_create(matrix_v_RI_kp(ikp, 2)%matrix, template=matrix_s_RI_aux_transl(1, 1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_reserve_all_blocks(matrix_v_RI_kp(ikp, 2)%matrix)
         CALL dbcsr_set(matrix_v_RI_kp(ikp, 2)%matrix, 0.0_dp)

      END DO

      CALL get_cell(cell=cell, abc=abc, periodic=periodic)
      IF (periodic(1) == 0) THEN
         abc(1) = 1.0E10_dp
      END IF
      IF (periodic(2) == 0) THEN
         abc(2) = 1.0E10_dp
      END IF
      IF (periodic(3) == 0) THEN
         abc(3) = 1.0E10_dp
      END IF
      cutoff = 1.0_dp/(qs_env%mp2_env%ri_rpa_im_time%cutoff*MIN(abc(1), abc(2), abc(3)))

      CALL build_2c_coulomb_matrix_kp(matrix_v_RI_kp, kpoints, basis_type="RI_AUX", &
                                      cell=cell, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                                      atomic_kind_set=atomic_kind_set, &
                                      size_lattice_sum=qs_env%mp2_env%mp2_gpw%size_lattice_sum, &
                                      operator_type=operator_coulomb, cutoff=cutoff)

      DO ikp = 1, nkp

         CALL copy_dbcsr_to_fm(matrix_v_RI_kp(ikp, 1)%matrix, fm_matrix_L_kpoints(ikp, 1)%matrix)
         CALL copy_dbcsr_to_fm(matrix_v_RI_kp(ikp, 2)%matrix, fm_matrix_L_kpoints(ikp, 2)%matrix)

      END DO

      CALL dbcsr_deallocate_matrix_set(matrix_v_RI_kp)

      CALL timestop(handle)

   END SUBROUTINE compute_V_by_lattice_sum

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param fm_matrix_L_RI_metric ...
!> \param fm_matrix_L ...
!> \param dimen_RI ...
!> \param ri_metric ...
!> \param do_kpoints ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE RI_2c_integral_mat(qs_env, fm_matrix_L_RI_metric, fm_matrix_L, dimen_RI, &
                                 ri_metric, do_kpoints, kpoints)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_RI_metric
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L
      INTEGER, INTENT(IN)                                :: dimen_RI
      TYPE(libint_potential_type)                        :: ri_metric
      LOGICAL, INTENT(IN)                                :: do_kpoints
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'RI_2c_integral_mat'

      INTEGER                                            :: handle, i_real_imag, i_size, ikp, img, &
                                                            n_real_imag, natom, nimg, nkind, nkp
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_RI
      INTEGER, DIMENSION(:), POINTER                     :: col_bsize, row_bsize
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_S_global
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_RI_aux_transl
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: mat_2c
      TYPE(dbcsr_type), POINTER                          :: cmatrix, matrix_s_RI_aux_desymm, &
                                                            rmatrix, tmpmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_RI
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_RI
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (sab_RI, matrix_s_RI_aux_transl)

      CALL get_qs_env(qs_env=qs_env, &
                      para_env=para_env, &
                      blacs_env=blacs_env, &
                      nkind=nkind, &
                      distribution_2d=dist_2d, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      dft_control=dft_control, &
                      natom=natom)

      ALLOCATE (sizes_RI(natom))
      ALLOCATE (basis_set_RI(nkind))

      CALL basis_set_list_setup(basis_set_RI, 'RI_AUX', qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=basis_set_RI)

      CALL build_2c_neighbor_lists(sab_RI, basis_set_RI, basis_set_RI, ri_metric, "RPA_2c_nl_RI", qs_env, &
                                   sym_ij=.TRUE., dist_2d=dist_2d)

      CALL cp_dbcsr_dist2d_to_dist(dist_2d, dbcsr_dist)
      ALLOCATE (row_bsize(SIZE(sizes_RI)))
      ALLOCATE (col_bsize(SIZE(sizes_RI)))
      row_bsize(:) = sizes_RI
      col_bsize(:) = sizes_RI

      IF (do_kpoints) THEN
         nimg = dft_control%nimages
      ELSE
         nimg = 1
      ENDIF

      ALLOCATE (mat_2c(nimg))
      CALL dbcsr_create(mat_2c(1), "(RI|RI)", dbcsr_dist, dbcsr_type_symmetric, &
                        row_bsize, col_bsize, reuse_arrays=.TRUE.)

      DO img = 2, nimg
         CALL dbcsr_create(mat_2c(img), template=mat_2c(1))
      ENDDO

      CALL build_2c_integrals(mat_2c, 0.0_dp, qs_env, sab_RI, basis_set_RI, basis_set_RI, &
                              ri_metric, do_kpoints=do_kpoints)

      CALL dbcsr_distribution_release(dbcsr_dist)
      DEALLOCATE (basis_set_RI)

      CALL dbcsr_allocate_matrix_set(matrix_s_RI_aux_transl, 1, nimg)
      DO img = 1, nimg
         ALLOCATE (matrix_s_RI_aux_transl(1, img)%matrix)
         matrix_s_RI_aux_transl(1, img)%matrix = mat_2c(img)
      ENDDO

      CALL set_ks_env(qs_env%ks_env, matrix_s_RI_aux_kp=matrix_s_RI_aux_transl)

      IF (do_kpoints) THEN
         CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, &
                              cell_to_index=cell_to_index)
         n_real_imag = 2
      ELSE
         nkp = 1
         n_real_imag = 1
      END IF

      ALLOCATE (fm_matrix_L_RI_metric(nkp, n_real_imag))
      DO i_size = 1, nkp
         DO i_real_imag = 1, n_real_imag
            CALL cp_fm_create(fm_matrix_L_RI_metric(i_size, i_real_imag)%matrix, fm_matrix_L%matrix_struct)
            CALL cp_fm_set_all(fm_matrix_L_RI_metric(i_size, i_real_imag)%matrix, 0.0_dp)
         END DO
      END DO

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=dimen_RI, &
                               ncol_global=dimen_RI, para_env=para_env)

      NULLIFY (fm_matrix_S_global)

      CALL cp_fm_create(fm_matrix_S_global, fm_struct)
      CALL cp_fm_set_all(fm_matrix_S_global, 0.0_dp)

      IF (do_kpoints) THEN

         ALLOCATE (rmatrix, cmatrix, tmpmat)
         CALL dbcsr_create(rmatrix, template=matrix_s_RI_aux_transl(1, 1)%matrix, &
                           matrix_type=dbcsr_type_symmetric)
         CALL dbcsr_create(cmatrix, template=matrix_s_RI_aux_transl(1, 1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL dbcsr_create(tmpmat, template=matrix_s_RI_aux_transl(1, 1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_RI)
         CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_RI)

         DO ikp = 1, nkp

            ! s matrix is not spin dependent, double the work
            CALL dbcsr_set(rmatrix, 0.0_dp)
            CALL dbcsr_set(cmatrix, 0.0_dp)

            CALL rskp_transform(rmatrix=rmatrix, cmatrix=cmatrix, rsmat=matrix_s_RI_aux_transl, ispin=1, &
                                xkp=xkp(1:3, ikp), cell_to_index=cell_to_index, sab_nl=sab_RI)

            CALL dbcsr_set(tmpmat, 0.0_dp)
            CALL dbcsr_desymmetrize(rmatrix, tmpmat)

            CALL cp_fm_set_all(fm_matrix_S_global, 0.0_dp)
            CALL copy_dbcsr_to_fm(tmpmat, fm_matrix_S_global)
            CALL cp_fm_copy_general(fm_matrix_S_global, fm_matrix_L_RI_metric(ikp, 1)%matrix, para_env)

            CALL dbcsr_set(tmpmat, 0.0_dp)
            CALL dbcsr_desymmetrize(cmatrix, tmpmat)

            CALL cp_fm_set_all(fm_matrix_S_global, 0.0_dp)
            CALL copy_dbcsr_to_fm(tmpmat, fm_matrix_S_global)
            CALL cp_fm_copy_general(fm_matrix_S_global, fm_matrix_L_RI_metric(ikp, 2)%matrix, para_env)

         END DO

         CALL dbcsr_deallocate_matrix(rmatrix)
         CALL dbcsr_deallocate_matrix(cmatrix)
         CALL dbcsr_deallocate_matrix(tmpmat)

      ELSE

         NULLIFY (matrix_s_RI_aux_desymm)
         ALLOCATE (matrix_s_RI_aux_desymm)
         CALL dbcsr_create(matrix_s_RI_aux_desymm, template=matrix_s_RI_aux_transl(1, 1)%matrix, &
                           name='S_RI non_symm', matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_desymmetrize(matrix_s_RI_aux_transl(1, 1)%matrix, matrix_s_RI_aux_desymm)

         CALL copy_dbcsr_to_fm(matrix_s_RI_aux_desymm, fm_matrix_S_global)

         CALL cp_fm_copy_general(fm_matrix_S_global, fm_matrix_L_RI_metric(1, 1)%matrix, para_env)

         CALL dbcsr_deallocate_matrix(matrix_s_RI_aux_desymm)

      END IF

      CALL release_neighbor_list_sets(sab_RI)

      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_release(fm_matrix_S_global)

      CALL timestop(handle)

   END SUBROUTINE RI_2c_integral_mat

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
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param gd_array ...
!> \param calc_PQ_cond_num ...
!> \param cond_num ...
!> \param num_small_eigen ...
!> \param potential ...
!> \param sab_orb_sub ...
!> \param do_im_time ...
!> \param fm_matrix_L_extern ...
! **************************************************************************************************
   SUBROUTINE compute_2c_integrals(qs_env, eri_method, eri_param, para_env, para_env_sub, para_env_L, mp2_memory, &
                                   fm_matrix_L, ngroup, color_sub, dimen_RI, mo_coeff, &
                                   my_group_L_size, my_group_L_start, my_group_L_end, &
                                   gd_array, calc_PQ_cond_num, cond_num, num_small_eigen, potential, &
                                   sab_orb_sub, do_im_time, fm_matrix_L_extern)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: eri_method
      TYPE(cp_eri_mme_param), POINTER                    :: eri_param
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub, para_env_L
      REAL(KIND=dp), INTENT(IN)                          :: mp2_memory
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L
      INTEGER, INTENT(IN)                                :: ngroup, color_sub
      INTEGER, INTENT(OUT)                               :: dimen_RI
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      INTEGER, INTENT(OUT)                               :: my_group_L_size, my_group_L_start, &
                                                            my_group_L_end
      TYPE(group_dist_d1_type), INTENT(OUT)              :: gd_array
      LOGICAL, INTENT(IN)                                :: calc_PQ_cond_num
      REAL(KIND=dp), INTENT(OUT)                         :: cond_num
      INTEGER, INTENT(OUT)                               :: num_small_eigen
      TYPE(libint_potential_type), INTENT(IN)            :: potential
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_im_time
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_matrix_L_extern

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_2c_integrals'

      INTEGER :: best_group_size, color_L, comm_exchange, comm_L, group_size, handle, handle2, i, &
         i_global, iatom, iiB, ikind, iproc, j_global, jjB, natom, ncol_local, nkind, nrow_local, &
         nsgf, potential_type, proc_receive, proc_receive_static, proc_send, proc_send_static, &
         proc_shift, rec_L_end, rec_L_size, rec_L_start, strat_group_size, sub_sub_color
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of, proc_map
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_do_im_time
      REAL(KIND=dp)                                      :: min_mem_for_QK
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: egen_L
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: L_external_col, L_local_col
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_L
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L_diag
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange
      TYPE(group_dist_d1_type)                           :: gd_sub_array
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      my_do_im_time = .FALSE.
      IF (PRESENT(do_im_time)) THEN
         my_do_im_time = do_im_time
      END IF

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      particle_set=particle_set)

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
         dimen_RI = dimen_RI + nsgf
      END DO

      ! check that very small systems are not running on too many processes
      CPASSERT(dimen_RI >= para_env%num_pe)

      CALL create_group_dist(gd_array, ngroup, dimen_RI)

      CALL get_group_dist(gd_array, color_sub, my_group_L_start, my_group_L_end, my_group_L_size)

      CALL timeset(routineN//"_loop_lm", handle2)

      ALLOCATE (L_local_col(dimen_RI, my_group_L_size))
      L_local_col = 0.0_dp

      potential_type = potential%potential_type

      IF ((eri_method .EQ. do_eri_mme .OR. eri_method .EQ. do_eri_os) &
          .AND. potential_type .EQ. do_potential_coulomb .OR. &
          (eri_method .EQ. do_eri_mme .AND. potential_type .EQ. do_potential_long)) THEN

         CALL mp2_eri_2c_integrate(eri_param, potential_type, potential%omega, para_env_sub, qs_env, &
                                   basis_type_a="RI_AUX", basis_type_b="RI_AUX", &
                                   hab=L_local_col, first_b=my_group_L_start, last_b=my_group_L_end, &
                                   eri_method=eri_method)

      ELSEIF (eri_method .EQ. do_eri_gpw .OR. &
              (potential_type == do_potential_long .AND. qs_env%mp2_env%eri_method == do_eri_os) &
              .OR. (potential_type == do_potential_id .AND. qs_env%mp2_env%eri_method == do_eri_mme)) THEN

         CALL mp2_eri_2c_integrate_gpw(qs_env, para_env_sub, dimen_RI, mo_coeff, my_group_L_start, my_group_L_end, &
                                       potential_type, &
                                       natom, potential%omega, sab_orb_sub, L_local_col, kind_of)

      ELSE
         CPABORT("unknown ERI method")
      ENDIF

      CALL timestop(handle2)

      ! split the total number of proc in a subgroup of the size of ~1/10 of the
      ! total num of proc
      best_group_size = para_env%num_pe

      strat_group_size = MAX(1, para_env%num_pe/10)

      min_mem_for_QK = REAL(dimen_RI, KIND=dp)*dimen_RI*3.0_dp*8.0_dp/1024_dp/1024_dp

      group_size = strat_group_size - 1
      DO iproc = strat_group_size, para_env%num_pe
         group_size = group_size + 1
         ! check that group_size is a multiple of sub_group_size and a divisor of
         ! the total num of proc
         IF (MOD(para_env%num_pe, group_size) /= 0 .OR. MOD(group_size, para_env_sub%num_pe) /= 0) CYCLE

         ! check for memory
         IF (REAL(group_size, KIND=dp)*mp2_memory < min_mem_for_QK) CYCLE

         best_group_size = group_size
         EXIT
      END DO

      IF (my_do_im_time) THEN
         ! para_env_L is the para_env for im. time to avoid bug
         best_group_size = para_env%num_pe
      END IF

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
      ALLOCATE (proc_map(-para_env_exchange%num_pe:2*para_env_exchange%num_pe - 1))
      DO i = 0, para_env_exchange%num_pe - 1
         proc_map(i) = i
         proc_map(-i - 1) = para_env_exchange%num_pe - i - 1
         proc_map(para_env_exchange%num_pe + i) = i
      END DO

      CALL create_group_dist(gd_sub_array, my_group_L_start, &
                             my_group_L_end, my_group_L_size, para_env_exchange)

      ! create the full matrix L defined in the L group

      CALL create_matrix_L(fm_matrix_L, blacs_env_L, dimen_RI, para_env_L, "fm_matrix_L", fm_matrix_L_extern)

      IF (my_do_im_time .AND. para_env%num_pe > 1) THEN

         CALL fill_fm_L_from_L_loc_non_blocking(fm_matrix_L, L_local_col, para_env, &
                                                my_group_L_start, my_group_L_end, &
                                                dimen_RI)

      ELSE

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
                  fm_matrix_L%local_data(iiB, jjB) = L_local_col(i_global, j_global - my_group_L_start + 1)
               END DO
            END IF
         END DO

         proc_send_static = proc_map(para_env_exchange%mepos + 1)
         proc_receive_static = proc_map(para_env_exchange%mepos - 1)

         DO proc_shift = 1, para_env_exchange%num_pe - 1
            proc_send = proc_map(para_env_exchange%mepos + proc_shift)
            proc_receive = proc_map(para_env_exchange%mepos - proc_shift)

            CALL get_group_dist(gd_sub_array, proc_receive, rec_L_start, rec_L_end, rec_L_size)

            ALLOCATE (L_external_col(dimen_RI, rec_L_size))
            L_external_col = 0.0_dp

            CALL mp_sendrecv(L_local_col, proc_send_static, L_external_col, proc_receive_static, para_env_exchange%group)

            DO jjB = 1, ncol_local
               j_global = col_indices(jjB)
               IF (j_global >= rec_L_start .AND. j_global <= rec_L_end) THEN
                  DO iiB = 1, nrow_local
                     i_global = row_indices(iiB)
                     fm_matrix_L%local_data(iiB, jjB) = L_external_col(i_global, j_global - rec_L_start + 1)
                  END DO
               END IF
            END DO

            DEALLOCATE (L_local_col)
            ALLOCATE (L_local_col(dimen_RI, rec_L_size))
            L_local_col(:, :) = L_external_col

            DEALLOCATE (L_external_col)

         END DO

      END IF

      DEALLOCATE (L_local_col)

      ! free the old exchange group stuff
      DEALLOCATE (proc_map)
      CALL cp_para_env_release(para_env_exchange)

      CALL release_group_dist(gd_sub_array)

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

         CALL create_matrix_L(fm_matrix_L_diag, blacs_env_L, dimen_RI, para_env_L, "fm_matrix_L_diag", fm_matrix_L_extern)

         CALL cp_fm_to_fm(source=fm_matrix_L, destination=fm_matrix_L_diag)

         ALLOCATE (egen_L(dimen_RI))

         egen_L = 0.0_dp
         CALL cp_fm_syevx(matrix=fm_matrix_L_diag, eigenvalues=egen_L)

         num_small_eigen = 0
         DO iiB = 1, dimen_RI
            IF (ABS(egen_L(iiB)) < 0.001_dp) num_small_eigen = num_small_eigen + 1
         END DO

         cond_num = MAXVAL(ABS(egen_L))/MINVAL(ABS(egen_L))

         CALL cp_fm_release(fm_matrix_L_diag)

         DEALLOCATE (egen_L)
      END IF

      ! release blacs_env
      CALL cp_blacs_env_release(blacs_env_L)

      CALL timestop(handle)

   END SUBROUTINE compute_2c_integrals

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eps_svd ...
!> \param num_small_evals ...
!> \param cond_num ...
!> \param do_inversion ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE matrix_root_with_svd(matrix, eps_svd, num_small_evals, cond_num, do_inversion, para_env)
      TYPE(cp_fm_type), POINTER                          :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: eps_svd
      INTEGER, INTENT(OUT)                               :: num_small_evals
      REAL(KIND=dp), INTENT(OUT)                         :: cond_num
      LOGICAL, INTENT(IN)                                :: do_inversion
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_root_with_svd'

      INTEGER                                            :: comm_exchange, group_size_L, handle, ii, &
                                                            needed_evals, nrow, pos_max
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_eval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: evals
      TYPE(cp_fm_type), POINTER                          :: evecs
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange

      CALL timeset(routineN, handle)

      CPASSERT(eps_svd >= 0)

      ! Create arrays carrying eigenvalues and eigenvectors later
      CALL cp_fm_get_info(matrix=matrix, nrow_global=nrow)

      ALLOCATE (evals(nrow))
      evals = 0

      NULLIFY (evecs)
      CALL cp_fm_create(evecs, matrix%matrix_struct)

      ! Perform the EVD (=SVD of a positive semidefinite matrix)
      CALL choose_eigv_solver(matrix, evecs, evals)

      ! Determine the number of neglectable eigenvalues assuming that the eigenvalues are in ascending order
      num_small_evals = 0
      DO ii = 1, nrow
         IF (evals(ii) > eps_svd) THEN
            num_small_evals = ii - 1
            EXIT
         END IF
      END DO
      needed_evals = nrow - num_small_evals

      ! Get the condition number w.r.t. considered eigenvalues
      cond_num = evals(nrow)/evals(num_small_evals + 1)

      ! Determine the eigenvalues of the request matrix root or its inverse
      evals(1:num_small_evals) = 0.0_dp
      IF (do_inversion) THEN
         evals(num_small_evals + 1:nrow) = 1.0_dp/SQRT(evals(num_small_evals + 1:nrow))
      ELSE
         evals(num_small_evals + 1:nrow) = SQRT(evals(num_small_evals + 1:nrow))
      END IF

      CALL cp_fm_column_scale(evecs, evals)

      ! As it turns out, the results in the subgroups differ.
      ! Thus, we have to equilibrate the results.
      ! Step 1: Get a communicator connecting processes with same id within subgroups
      ! use that group_size_L is divisible by the total number of processes (see above)
      group_size_L = para_env%num_pe/matrix%matrix_struct%para_env%num_pe
      CALL mp_comm_split_direct(para_env%group, comm_exchange, matrix%matrix_struct%para_env%mepos)
      NULLIFY (para_env_exchange)
      CALL cp_para_env_create(para_env_exchange, comm_exchange)

      ALLOCATE (num_eval(0:group_size_L - 1))
      num_eval = 0
      num_eval(para_env_exchange%mepos) = num_small_evals
      CALL mp_sum(num_eval, para_env_exchange%group)

      num_small_evals = MINVAL(num_eval)

      IF (num_small_evals /= MAXVAL(num_eval)) THEN
         ! Step 2: Get position of maximum value
         DO ii = 0, group_size_L - 1
            IF (num_eval(ii) == num_small_evals) THEN
               pos_max = ii
               EXIT
            END IF
         END DO
         num_small_evals = num_eval(pos_max)
         needed_evals = nrow - num_small_evals

         ! Step 3: Broadcast your local data to all other processes
         CALL mp_bcast(evecs%local_data, pos_max, para_env_exchange%group)
         CALL mp_bcast(cond_num, pos_max, para_env_exchange%group)
      END IF

      DEALLOCATE (num_eval)

      CALL cp_para_env_release(para_env_exchange)

      CALL reset_size_matrix(matrix, needed_evals, matrix%matrix_struct)

      ! Copy the needed eigenvectors
      CALL cp_fm_to_fm(evecs, matrix, needed_evals, num_small_evals + 1)

      CALL cp_fm_release(evecs)

      CALL timestop(handle)

   END SUBROUTINE matrix_root_with_svd

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param new_size ...
!> \param fm_struct_template ...
! **************************************************************************************************
   SUBROUTINE reset_size_matrix(matrix, new_size, fm_struct_template)
      TYPE(cp_fm_type), POINTER                          :: matrix
      INTEGER, INTENT(IN)                                :: new_size
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_template

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'reset_size_matrix'

      INTEGER                                            :: handle
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      ! Choose the block sizes as large as in the template for the later steps
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, ncol_global=new_size, template_fmstruct=fm_struct_template, force_block=.TRUE.)

      CALL cp_fm_release(matrix)

      NULLIFY (matrix)
      CALL cp_fm_create(matrix, fm_struct)
      CALL cp_fm_set_all(matrix, 0.0_dp)

      CALL cp_fm_struct_release(fm_struct)

      CALL timestop(handle)

   END SUBROUTINE reset_size_matrix

! **************************************************************************************************
!> \brief ...
!> \param fm_matrix_L ...
!> \param L_local_col ...
!> \param para_env ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param dimen_RI ...
! **************************************************************************************************
   SUBROUTINE fill_fm_L_from_L_loc_non_blocking(fm_matrix_L, L_local_col, para_env, my_group_L_start, &
                                                my_group_L_end, dimen_RI)
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: L_local_col
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: my_group_L_start, my_group_L_end, &
                                                            dimen_RI

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_fm_L_from_L_loc_non_blocking'

      INTEGER :: dummy_proc, handle, handle2, i_entry_rec, i_row, i_row_global, iproc, j_col, &
         j_col_global, LLL, MMM, ncol_block, ncol_local, npcol, nprow, nrow_block, nrow_local, &
         proc_send, send_pcol, send_prow
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: entry_counter, num_entries_rec, &
                                                            num_entries_send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      CALL timeset(routineN//"_1", handle2)

      ! get info for the dest
      CALL cp_fm_get_info(matrix=fm_matrix_L, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices, &
                          nrow_block=nrow_block, &
                          ncol_block=ncol_block)

      nprow = fm_matrix_L%matrix_struct%context%num_pe(1)
      npcol = fm_matrix_L%matrix_struct%context%num_pe(2)

      ALLOCATE (num_entries_rec(0:para_env%num_pe - 1))
      ALLOCATE (num_entries_send(0:para_env%num_pe - 1))

      num_entries_rec(:) = 0
      num_entries_send(:) = 0

      dummy_proc = 0

      ! get the process, where the elements from L_local_col have to be sent
      DO LLL = 1, dimen_RI

         send_prow = cp_fm_indxg2p(LLL, nrow_block, dummy_proc, &
                                   fm_matrix_L%matrix_struct%first_p_pos(1), nprow)

         DO MMM = my_group_L_start, my_group_L_end

            send_pcol = cp_fm_indxg2p(MMM, ncol_block, dummy_proc, &
                                      fm_matrix_L%matrix_struct%first_p_pos(2), npcol)

            proc_send = fm_matrix_L%matrix_struct%context%blacs2mpi(send_prow, send_pcol)

            num_entries_send(proc_send) = num_entries_send(proc_send) + 1

         END DO

      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_2", handle2)

      CALL mp_alltoall(num_entries_send, num_entries_rec, 1, para_env%group)

      CALL timestop(handle2)

      CALL timeset(routineN//"_3", handle2)

      ! allocate buffers to send the entries and the information of the entries
      ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env%num_pe - 1))

      ! allocate data message and corresponding indices
      DO iproc = 0, para_env%num_pe - 1

         ALLOCATE (buffer_rec(iproc)%msg(num_entries_rec(iproc)))
         buffer_rec(iproc)%msg = 0.0_dp

      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_4", handle2)

      DO iproc = 0, para_env%num_pe - 1

         ALLOCATE (buffer_send(iproc)%msg(num_entries_send(iproc)))
         buffer_send(iproc)%msg = 0.0_dp

      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_5", handle2)

      DO iproc = 0, para_env%num_pe - 1

         ALLOCATE (buffer_rec(iproc)%indx(num_entries_rec(iproc), 2))
         buffer_rec(iproc)%indx = 0

      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_6", handle2)

      DO iproc = 0, para_env%num_pe - 1

         ALLOCATE (buffer_send(iproc)%indx(num_entries_send(iproc), 2))
         buffer_send(iproc)%indx = 0

      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_7", handle2)

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = 0

      ! get the process, where the elements from L_local_col have to be sent and
      ! write the data and the info to buffer_send
      DO LLL = 1, dimen_RI

         send_prow = cp_fm_indxg2p(LLL, nrow_block, dummy_proc, &
                                   fm_matrix_L%matrix_struct%first_p_pos(1), nprow)

         DO MMM = my_group_L_start, my_group_L_end

            send_pcol = cp_fm_indxg2p(MMM, ncol_block, dummy_proc, &
                                      fm_matrix_L%matrix_struct%first_p_pos(2), npcol)

            proc_send = fm_matrix_L%matrix_struct%context%blacs2mpi(send_prow, send_pcol)

            entry_counter(proc_send) = entry_counter(proc_send) + 1

            buffer_send(proc_send)%msg(entry_counter(proc_send)) = &
               L_local_col(LLL, MMM - my_group_L_start + 1)

            buffer_send(proc_send)%indx(entry_counter(proc_send), 1) = LLL
            buffer_send(proc_send)%indx(entry_counter(proc_send), 2) = MMM

         END DO

      END DO

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      CALL timestop(handle2)

      CALL timeset(routineN//"_8", handle2)

      ! communicate the buffer
      CALL communicate_buffer(para_env, num_entries_rec, num_entries_send, buffer_rec, &
                              buffer_send, req_array)

      fm_matrix_L%local_data = 0.0_dp

      CALL timestop(handle2)

      CALL timeset(routineN//"_9", handle2)

      ! fill fm_matrix_L with the entries from buffer_rec
      DO iproc = 0, para_env%num_pe - 1

         DO i_entry_rec = 1, num_entries_rec(iproc)

            DO i_row = 1, nrow_local

               i_row_global = row_indices(i_row)

               DO j_col = 1, ncol_local

                  j_col_global = col_indices(j_col)

                  IF (i_row_global == buffer_rec(iproc)%indx(i_entry_rec, 1) .AND. &
                      j_col_global == buffer_rec(iproc)%indx(i_entry_rec, 2)) THEN

                     fm_matrix_L%local_data(i_row, j_col) = buffer_rec(iproc)%msg(i_entry_rec)

                  END IF

               END DO

            END DO

         END DO

      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_10", handle2)

      DO iproc = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_rec(iproc)%msg)
         DEALLOCATE (buffer_rec(iproc)%indx)
         DEALLOCATE (buffer_send(iproc)%msg)
         DEALLOCATE (buffer_send(iproc)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)
      DEALLOCATE (req_array)
      DEALLOCATE (entry_counter)
      DEALLOCATE (num_entries_rec, num_entries_send)

      CALL timestop(handle2)

      CALL timestop(handle)

   END SUBROUTINE fill_fm_L_from_L_loc_non_blocking

! **************************************************************************************************
!> \brief ...
!> \param fm_matrix_L ...
!> \param dimen_RI ...
!> \param do_inversion ...
! **************************************************************************************************
   SUBROUTINE cholesky_decomp(fm_matrix_L, dimen_RI, do_inversion)

      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L
      INTEGER, INTENT(IN)                                :: dimen_RI
      LOGICAL, INTENT(IN)                                :: do_inversion

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cholesky_decomp'

      INTEGER                                            :: handle, i_global, iiB, info_chol, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      ! do cholesky decomposition
      CALL cp_fm_cholesky_decompose(matrix=fm_matrix_L, n=dimen_RI, info_out=info_chol)
      CPASSERT(info_chol == 0)

      IF (.NOT. do_inversion) THEN
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

      END IF

      CALL timestop(handle)

   END SUBROUTINE cholesky_decomp

! **************************************************************************************************
!> \brief ...
!> \param fm_matrix_L_RI_metric ...
!> \param fm_matrix_L_kpoints ...
!> \param dimen_RI ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE inversion_of_S_and_mult_with_chol_dec_of_V(fm_matrix_L_RI_metric, fm_matrix_L_kpoints, &
                                                         dimen_RI, kpoints)

      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_RI_metric, &
                                                            fm_matrix_L_kpoints
      INTEGER, INTENT(IN)                                :: dimen_RI
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'inversion_of_S_and_mult_with_chol_dec_of_V'
      COMPLEX(KIND=dp), PARAMETER :: cone = CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
         czero = CMPLX(0.0_dp, 0.0_dp, KIND=dp), ione = CMPLX(0.0_dp, 1.0_dp, KIND=dp)

      INTEGER                                            :: handle, i_global, iiB, ikp, info_chol, &
                                                            j_global, jjB, ncol_local, nkp, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_cfm_type), POINTER                         :: cfm_matrix_K_tmp, cfm_matrix_L_tmp, &
                                                            cfm_matrix_S_inv_tmp, cfm_matrix_S_tmp
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct

      CALL timeset(routineN, handle)

      NULLIFY (cfm_matrix_S_tmp, cfm_matrix_S_inv_tmp, cfm_matrix_L_tmp, cfm_matrix_K_tmp)

      CALL cp_fm_get_info(fm_matrix_L_RI_metric(1, 1)%matrix, matrix_struct=matrix_struct)

      CALL cp_cfm_create(cfm_matrix_S_tmp, matrix_struct)
      CALL cp_cfm_create(cfm_matrix_S_inv_tmp, matrix_struct)
      CALL cp_cfm_create(cfm_matrix_L_tmp, matrix_struct)
      CALL cp_cfm_create(cfm_matrix_K_tmp, matrix_struct)

      CALL get_kpoint_info(kpoints, nkp=nkp)

      DO ikp = 1, nkp

         CALL cp_cfm_scale_and_add_fm(czero, cfm_matrix_S_tmp, cone, fm_matrix_L_RI_metric(ikp, 1)%matrix)
         CALL cp_cfm_scale_and_add_fm(cone, cfm_matrix_S_tmp, ione, fm_matrix_L_RI_metric(ikp, 2)%matrix)

         CALL cp_cfm_scale_and_add_fm(czero, cfm_matrix_L_tmp, cone, fm_matrix_L_kpoints(ikp, 1)%matrix)
         CALL cp_cfm_scale_and_add_fm(cone, cfm_matrix_L_tmp, ione, fm_matrix_L_kpoints(ikp, 2)%matrix)

         ! Cholesky decomposition S=U(T)U
         CALL cp_cfm_cholesky_decompose(cfm_matrix_S_tmp, info_out=info_chol)
         CPASSERT(info_chol == 0)

         ! Invert to get U^(-1)
         CALL cp_cfm_triangular_invert(cfm_matrix_S_tmp, uplo='U')

         ! clean the lower part of the U^{-1} matrix (just to not have surprises afterwards)
         CALL cp_cfm_get_info(matrix=cfm_matrix_S_tmp, &
                              nrow_local=nrow_local, &
                              ncol_local=ncol_local, &
                              row_indices=row_indices, &
                              col_indices=col_indices)

         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            DO jjB = 1, ncol_local
               j_global = col_indices(jjB)
               IF (j_global < i_global) cfm_matrix_S_tmp%local_data(iiB, jjB) = czero
            END DO
         END DO

         ! get S^(-1) = U^(-1)U^(-H)
         CALL cp_cfm_gemm("N", "C", dimen_RI, dimen_RI, dimen_RI, cone, cfm_matrix_S_tmp, cfm_matrix_S_tmp, &
                          czero, cfm_matrix_S_inv_tmp)

         ! Cholesky decomposition V=L(T)L
         CALL cp_cfm_cholesky_decompose(cfm_matrix_L_tmp, info_out=info_chol)
         CPASSERT(info_chol == 0)

         ! clean the lower part of the L matrix (just to not have surprises afterwards)
         CALL cp_cfm_get_info(matrix=cfm_matrix_L_tmp, &
                              nrow_local=nrow_local, &
                              ncol_local=ncol_local, &
                              row_indices=row_indices, &
                              col_indices=col_indices)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            DO jjB = 1, ncol_local
               j_global = col_indices(jjB)
               IF (j_global < i_global) cfm_matrix_L_tmp%local_data(iiB, jjB) = czero
            END DO
         END DO

         ! get K = S^(-1)*L, where L is the Cholesky decomposition of V = L^T*L
         CALL cp_cfm_gemm("N", "C", dimen_RI, dimen_RI, dimen_RI, cone, cfm_matrix_S_inv_tmp, cfm_matrix_L_tmp, &
                          czero, cfm_matrix_K_tmp)

         ! move
         CALL cp_cfm_to_fm(cfm_matrix_K_tmp, fm_matrix_L_RI_metric(ikp, 1)%matrix, fm_matrix_L_RI_metric(ikp, 2)%matrix)

      END DO

      CALL cp_cfm_release(cfm_matrix_S_tmp)
      CALL cp_cfm_release(cfm_matrix_S_inv_tmp)
      CALL cp_cfm_release(cfm_matrix_L_tmp)
      CALL cp_cfm_release(cfm_matrix_K_tmp)

      CALL timestop(handle)

   END SUBROUTINE inversion_of_S_and_mult_with_chol_dec_of_V

! **************************************************************************************************
!> \brief ...
!> \param fm_matrix_L ...
! **************************************************************************************************
   SUBROUTINE invert_L(fm_matrix_L)

      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'invert_L'

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_fm_triangular_invert(matrix_a=fm_matrix_L, uplo_tr='U')

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

      CALL timestop(handle)

   END SUBROUTINE invert_L

! **************************************************************************************************
!> \brief ...
!> \param fm_matrix_L ...
!> \param blacs_env_L ...
!> \param dimen_RI ...
!> \param para_env_L ...
!> \param name ...
!> \param fm_matrix_L_extern ...
! **************************************************************************************************
   SUBROUTINE create_matrix_L(fm_matrix_L, blacs_env_L, dimen_RI, para_env_L, name, fm_matrix_L_extern)
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_L
      INTEGER, INTENT(IN)                                :: dimen_RI
      TYPE(cp_para_env_type), POINTER                    :: para_env_L
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_matrix_L_extern

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_matrix_L'

      INTEGER                                            :: handle
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      ! create the full matrix L defined in the L group
      NULLIFY (fm_matrix_L)
      IF (PRESENT(fm_matrix_L_extern)) THEN
         CALL cp_fm_create(fm_matrix_L, fm_matrix_L_extern%matrix_struct, name=name)
      ELSE
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env_L, nrow_global=dimen_RI, &
                                  ncol_global=dimen_RI, para_env=para_env_L)

         CALL cp_fm_create(fm_matrix_L, fm_struct, name=name)

         CALL cp_fm_struct_release(fm_struct)
      END IF

      CALL cp_fm_set_all(matrix=fm_matrix_L, alpha=0.0_dp)

      CALL timestop(handle)

   END SUBROUTINE create_matrix_L
END MODULE mp2_ri_2c
