!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for GW
!> \par History
!>      03.2019 created [Frederick Stein]
! **************************************************************************************************
MODULE rpa_gw
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_cfm_types,                    ONLY: cp_cfm_p_type,&
                                              cp_cfm_release
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_diag,                      ONLY: cp_fm_syevd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add_on_diag, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, dbcsr_filter, &
        dbcsr_get_info, dbcsr_init_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_release, dbcsr_release_p, dbcsr_scalar, dbcsr_scale, dbcsr_set, &
        dbcsr_type, dbcsr_type_no_symmetry
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_batched_contract_finalize, dbcsr_t_batched_contract_init, dbcsr_t_clear, &
        dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_copy_matrix_to_tensor, dbcsr_t_create, &
        dbcsr_t_destroy, dbcsr_t_get_block, dbcsr_t_get_info, dbcsr_t_iterator_blocks_left, &
        dbcsr_t_iterator_next_block, dbcsr_t_iterator_start, dbcsr_t_iterator_stop, &
        dbcsr_t_iterator_type, dbcsr_t_nblks_total, dbcsr_t_pgrid_create, dbcsr_t_pgrid_destroy, &
        dbcsr_t_pgrid_type, dbcsr_t_type
   USE hfx_types,                       ONLY: block_ind_type,&
                                              dealloc_containers,&
                                              hfx_compression_type
   USE input_constants,                 ONLY: gw_pade_approx,&
                                              gw_two_pole_model,&
                                              ri_rpa_g0w0_crossing_bisection,&
                                              ri_rpa_g0w0_crossing_newton,&
                                              ri_rpa_g0w0_crossing_z_shot
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_release,&
                                              kpoint_sym_create,&
                                              kpoint_type
   USE mathconstants,                   ONLY: fourpi,&
                                              gaussi,&
                                              pi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_sum,&
                                              mp_sync
   USE mp2_types,                       ONLY: mp2_type,&
                                              one_dim_real_array,&
                                              two_dim_int_array
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE qs_band_structure,               ONLY: calculate_kp_orbitals
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_release,&
                                              qs_environment_type
   USE qs_gamma2kp,                     ONLY: create_kp_from_gamma
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_moments,                      ONLY: build_berry_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: setup_neighbor_list
   USE qs_overlap,                      ONLY: build_overlap_matrix_simple
   USE qs_tensors,                      ONLY: decompress_tensor
   USE qs_tensors_types,                ONLY: create_2c_tensor
   USE rpa_gw_ic,                       ONLY: apply_ic_corr
   USE rpa_gw_im_time_util,             ONLY: get_tensor_3c_overl_int_gw
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw'

   PUBLIC :: allocate_matrices_gw_im_time, allocate_matrices_gw, GW_matrix_operations, compute_QP_energies, &
             deallocate_matrices_gw_im_time, deallocate_matrices_gw

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param nmo ...
!> \param num_integ_points ...
!> \param unit_nr ...
!> \param RI_blk_sizes ...
!> \param do_ic_model ...
!> \param para_env ...
!> \param fm_mat_W ...
!> \param fm_mat_Q ...
!> \param mo_coeff ...
!> \param t_3c_overl_int_ao_mo ...
!> \param t_3c_O_mo_compressed ...
!> \param t_3c_O_mo_ind ...
!> \param t_3c_overl_int_gw_RI ...
!> \param t_3c_overl_int_gw_AO ...
!> \param starts_array_mc ...
!> \param ends_array_mc ...
!> \param t_3c_overl_nnP_ic ...
!> \param t_3c_overl_nnP_ic_reflected ...
!> \param matrix_s ...
!> \param mat_W ...
!> \param t_3c_overl_int ...
!> \param t_3c_O_compressed ...
!> \param t_3c_O_ind ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE allocate_matrices_gw_im_time(gw_corr_lev_occ, gw_corr_lev_virt, homo, nmo, &
                                           num_integ_points, unit_nr, &
                                           RI_blk_sizes, do_ic_model, &
                                           para_env, fm_mat_W, fm_mat_Q, &
                                           mo_coeff, &
                                           t_3c_overl_int_ao_mo, t_3c_O_mo_compressed, t_3c_O_mo_ind, &
                                           t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                           starts_array_mc, ends_array_mc, &
                                           t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected, &
                                           matrix_s, mat_W, t_3c_overl_int, &
                                           t_3c_O_compressed, t_3c_O_ind, &
                                           qs_env)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: gw_corr_lev_occ, gw_corr_lev_virt, homo
      INTEGER, INTENT(IN)                                :: nmo, num_integ_points, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: RI_blk_sizes
      LOGICAL, INTENT(IN)                                :: do_ic_model
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: fm_mat_W
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: mo_coeff
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_ao_mo
      TYPE(hfx_compression_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: t_3c_O_mo_compressed
      TYPE(two_dim_int_array), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: t_3c_O_mo_ind
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: t_3c_overl_int_gw_RI, &
                                                            t_3c_overl_int_gw_AO
      INTEGER, DIMENSION(:), INTENT(IN)                  :: starts_array_mc, ends_array_mc
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: t_3c_overl_nnP_ic, &
                                                            t_3c_overl_nnP_ic_reflected
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: mat_W
      TYPE(dbcsr_t_type), DIMENSION(:, :)                :: t_3c_overl_int
      TYPE(hfx_compression_type), DIMENSION(:, :, :)     :: t_3c_O_compressed
      TYPE(block_ind_type), DIMENSION(:, :, :)           :: t_3c_O_ind
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_matrices_gw_im_time'

      INTEGER                                            :: handle, jquad, nspins
      LOGICAL                                            :: my_open_shell
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_ao_mo_beta

      CALL timeset(routineN, handle)

      nspins = SIZE(homo)
      my_open_shell = (nspins == 2)

      ALLOCATE (t_3c_O_mo_ind(nspins), t_3c_overl_int_gw_AO(nspins), t_3c_overl_int_gw_RI(nspins), &
                t_3c_overl_nnP_ic(nspins), t_3c_overl_nnP_ic_reflected(nspins), t_3c_O_mo_compressed(nspins))
      CALL get_tensor_3c_overl_int_gw(t_3c_overl_int, &
                                      t_3c_O_compressed, t_3c_O_ind, &
                                      t_3c_overl_int_ao_mo, t_3c_O_mo_compressed(1), t_3c_O_mo_ind(1)%array, &
                                      t_3c_overl_int_gw_RI(1), t_3c_overl_int_gw_AO(1), &
                                      starts_array_mc, ends_array_mc, &
                                      mo_coeff(1)%matrix, matrix_s, &
                                      gw_corr_lev_occ(1), gw_corr_lev_virt(1), homo(1), nmo, &
                                      para_env, &
                                      do_ic_model, &
                                      t_3c_overl_nnP_ic(1), t_3c_overl_nnP_ic_reflected(1), &
                                      qs_env, unit_nr)

      IF (my_open_shell) THEN

         CALL get_tensor_3c_overl_int_gw(t_3c_overl_int, &
                                         t_3c_O_compressed, t_3c_O_ind, &
                                         t_3c_overl_int_ao_mo_beta, t_3c_O_mo_compressed(2), t_3c_O_mo_ind(2)%array, &
                                         t_3c_overl_int_gw_RI(2), t_3c_overl_int_gw_AO(2), &
                                         starts_array_mc, ends_array_mc, &
                                         mo_coeff(2)%matrix, matrix_s, &
                                         gw_corr_lev_occ(2), gw_corr_lev_virt(2), homo(2), nmo, &
                                         para_env, &
                                         do_ic_model, &
                                         t_3c_overl_nnP_ic(2), t_3c_overl_nnP_ic_reflected(2), &
                                         qs_env, unit_nr)

         CALL dbcsr_t_destroy(t_3c_overl_int_ao_mo_beta)

      END IF

      ALLOCATE (fm_mat_W(num_integ_points))

      DO jquad = 1, num_integ_points

         NULLIFY (fm_mat_W(jquad)%matrix)
         CALL cp_fm_create(fm_mat_W(jquad)%matrix, fm_mat_Q%matrix_struct)
         CALL cp_fm_to_fm(fm_mat_Q, fm_mat_W(jquad)%matrix)
         CALL cp_fm_set_all(fm_mat_W(jquad)%matrix, 0.0_dp)

      END DO

      NULLIFY (mat_W)
      CALL dbcsr_init_p(mat_W)
      CALL dbcsr_create(matrix=mat_W, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=RI_blk_sizes, &
                        col_blk_size=RI_blk_sizes)

      CALL timestop(handle)

   END SUBROUTINE allocate_matrices_gw_im_time

! **************************************************************************************************
!> \brief ...
!> \param vec_Sigma_c_gw ...
!> \param color_rpa_group ...
!> \param dimen_nm_gw ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param nmo ...
!> \param num_integ_group ...
!> \param num_integ_points ...
!> \param unit_nr ...
!> \param gw_corr_lev_tot ...
!> \param num_fit_points ...
!> \param omega_max_fit ...
!> \param do_minimax_quad ...
!> \param do_periodic ...
!> \param do_ri_Sigma_x ...
!> \param my_do_gw ...
!> \param first_cycle_periodic_correction ...
!> \param a_scaling ...
!> \param Eigenval ...
!> \param tj ...
!> \param vec_omega_fit_gw ...
!> \param vec_Sigma_x_gw ...
!> \param delta_corr ...
!> \param Eigenval_last ...
!> \param Eigenval_scf ...
!> \param vec_W_gw ...
!> \param fm_mat_S_gw ...
!> \param fm_mat_S_gw_work ...
!> \param para_env ...
!> \param mp2_env ...
!> \param kpoints ...
!> \param nkp ...
!> \param nkp_self_energy ...
!> \param do_kpoints_cubic_RPA ...
!> \param do_kpoints_from_Gamma ...
! **************************************************************************************************
   SUBROUTINE allocate_matrices_gw(vec_Sigma_c_gw, color_rpa_group, dimen_nm_gw, &
                                   gw_corr_lev_occ, gw_corr_lev_virt, homo, &
                                   nmo, num_integ_group, num_integ_points, unit_nr, &
                                   gw_corr_lev_tot, num_fit_points, omega_max_fit, &
                                   do_minimax_quad, do_periodic, do_ri_Sigma_x, my_do_gw, &
                                   first_cycle_periodic_correction, &
                                   a_scaling, Eigenval, tj, vec_omega_fit_gw, vec_Sigma_x_gw, &
                                   delta_corr, Eigenval_last, Eigenval_scf, vec_W_gw, &
                                   fm_mat_S_gw, fm_mat_S_gw_work, &
                                   para_env, mp2_env, kpoints, nkp, nkp_self_energy, &
                                   do_kpoints_cubic_RPA, do_kpoints_from_Gamma)

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), INTENT(OUT)              :: vec_Sigma_c_gw
      INTEGER, INTENT(IN)                                :: color_rpa_group, dimen_nm_gw
      INTEGER, DIMENSION(:), INTENT(IN)                  :: gw_corr_lev_occ, gw_corr_lev_virt, homo
      INTEGER, INTENT(IN)                                :: nmo, num_integ_group, num_integ_points, &
                                                            unit_nr
      INTEGER, INTENT(INOUT)                             :: gw_corr_lev_tot, num_fit_points
      REAL(KIND=dp)                                      :: omega_max_fit
      LOGICAL, INTENT(IN)                                :: do_minimax_quad, do_periodic, &
                                                            do_ri_Sigma_x, my_do_gw
      LOGICAL, INTENT(OUT) :: first_cycle_periodic_correction
      REAL(KIND=dp), INTENT(IN)                          :: a_scaling
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: vec_omega_fit_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: vec_Sigma_x_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: delta_corr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: Eigenval_last, Eigenval_scf, vec_W_gw
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_S_gw
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: fm_mat_S_gw_work
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(OUT)                               :: nkp, nkp_self_energy
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA, &
                                                            do_kpoints_from_Gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_matrices_gw'

      INTEGER                                            :: handle, iquad, ispin, jquad, nspins
      LOGICAL                                            :: my_open_shell
      REAL(KIND=dp)                                      :: omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vec_omega_gw

      CALL timeset(routineN, handle)

      nspins = SIZE(Eigenval, 2)
      my_open_shell = (nspins == 2)

      gw_corr_lev_tot = gw_corr_lev_occ(1) + gw_corr_lev_virt(1)

      ! fill the omega_frequency vector
      ALLOCATE (vec_omega_gw(num_integ_points))
      vec_omega_gw = 0.0_dp

      DO jquad = 1, num_integ_points
         IF (do_minimax_quad) THEN
            omega = tj(jquad)
         ELSE
            omega = a_scaling/TAN(tj(jquad))
         END IF
         vec_omega_gw(jquad) = omega
      END DO

      ! determine number of fit points in the interval [0,w_max] for virt, or [-w_max,0] for occ
      num_fit_points = 0

      DO jquad = 1, num_integ_points
         IF (vec_omega_gw(jquad) < omega_max_fit) THEN
            num_fit_points = num_fit_points + 1
         END IF
      END DO

      IF (mp2_env%ri_g0w0%analytic_continuation == gw_pade_approx) THEN
         IF (mp2_env%ri_g0w0%nparam_pade > num_fit_points) THEN
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A)") &
               "Pade approximation: more parameters than data points. Reset # of parameters."
            mp2_env%ri_g0w0%nparam_pade = num_fit_points
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T74,I7)") &
               "Number of pade parameters:", mp2_env%ri_g0w0%nparam_pade
         ENDIF
      ENDIF

      ! create new arrays containing omega values at which we calculate vec_Sigma_c_gw
      ALLOCATE (vec_omega_fit_gw(num_fit_points))

      ! fill the omega vector with frequencies, where we calculate the self-energy
      iquad = 0
      DO jquad = 1, num_integ_points
         IF (vec_omega_gw(jquad) < omega_max_fit) THEN
            iquad = iquad + 1
            vec_omega_fit_gw(iquad) = vec_omega_gw(jquad)
         END IF
      END DO

      DEALLOCATE (vec_omega_gw)

      IF (do_kpoints_cubic_RPA) THEN
         CALL get_kpoint_info(kpoints, nkp=nkp)
         IF (mp2_env%ri_g0w0%do_gamma_only_sigma) THEN
            nkp_self_energy = 1
         ELSE
            nkp_self_energy = nkp
         END IF
      ELSE IF (do_kpoints_from_Gamma) THEN
         CALL get_kpoint_info(kpoints, nkp=nkp)
         nkp_self_energy = 1
      ELSE
         nkp = 1
         nkp_self_energy = 1
      END IF
      ALLOCATE (vec_Sigma_c_gw(gw_corr_lev_tot, num_fit_points, nkp_self_energy, nspins))
      vec_Sigma_c_gw = z_zero

      ALLOCATE (Eigenval_scf(nmo, nspins))
      Eigenval_scf(:, :) = Eigenval(:, :)

      ALLOCATE (Eigenval_last(nmo, nspins))
      Eigenval_last(:, :) = Eigenval(:, :)

      IF (do_periodic) THEN

         ALLOCATE (delta_corr(1 + homo(1) - gw_corr_lev_occ(1):homo(1) + gw_corr_lev_virt(1)))
         delta_corr(:) = 0.0_dp

         first_cycle_periodic_correction = .TRUE.

      END IF

      IF (do_ri_Sigma_x) THEN
         ALLOCATE (vec_Sigma_x_gw(nmo, nkp_self_energy, nspins))
         vec_Sigma_x_gw = 0.0_dp
      END IF

      IF (my_do_gw) THEN

         ! minimax grids not implemented for O(N^4) GW
         CPASSERT(.NOT. do_minimax_quad)

         ! create temporary matrix to store B*([1+Q(iw')]^-1-1), has the same size as B
         ALLOCATE (fm_mat_S_gw_work(nspins))
         DO ispin = 1, nspins
            NULLIFY (fm_mat_S_gw_work(ispin)%matrix)
            CALL cp_fm_create(fm_mat_S_gw_work(ispin)%matrix, fm_mat_S_gw(ispin)%matrix%matrix_struct)
            CALL cp_fm_set_all(matrix=fm_mat_S_gw_work(ispin)%matrix, alpha=0.0_dp)
         END DO

         ALLOCATE (vec_W_gw(dimen_nm_gw, nspins))
         vec_W_gw = 0.0_dp

         ! in case we do RI for Sigma_x, we calculate Sigma_x right here
         IF (do_ri_Sigma_x) THEN

            CALL get_vec_sigma_x(vec_Sigma_x_gw(:, :, 1), nmo, fm_mat_S_gw(1)%matrix, para_env, num_integ_group, color_rpa_group, &
                                 homo(1), gw_corr_lev_occ(1), mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1, 1))

            IF (my_open_shell) THEN
               CALL get_vec_sigma_x(vec_Sigma_x_gw(:, :, 2), nmo, fm_mat_S_gw(2)%matrix, para_env, num_integ_group, &
                                    color_rpa_group, homo(2), gw_corr_lev_occ(2), &
                                    mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2, 1))
            END IF

         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE allocate_matrices_gw

! **************************************************************************************************
!> \brief ...
!> \param vec_Sigma_x_gw ...
!> \param nmo ...
!> \param fm_mat_S_gw ...
!> \param para_env ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param homo ...
!> \param gw_corr_lev_occ ...
!> \param vec_Sigma_x_minus_vxc_gw11 ...
! **************************************************************************************************
   SUBROUTINE get_vec_sigma_x(vec_Sigma_x_gw, nmo, fm_mat_S_gw, para_env, num_integ_group, color_rpa_group, homo, &
                              gw_corr_lev_occ, vec_Sigma_x_minus_vxc_gw11)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: vec_Sigma_x_gw
      INTEGER, INTENT(IN)                                :: nmo
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S_gw
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: num_integ_group, color_rpa_group, homo, &
                                                            gw_corr_lev_occ
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: vec_Sigma_x_minus_vxc_gw11

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_vec_sigma_x'

      INTEGER                                            :: handle, iiB, jjB, m_global, n_global, &
                                                            ncol_local, nm_global, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: row_indices

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm_mat_S_gw, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices)

      CALL mp_sync(para_env%group)

      ! loop over (nm) index
      DO iiB = 1, nrow_local

         ! this is needed for correct values within parallelization
         IF (MODULO(1, num_integ_group) /= color_rpa_group) CYCLE

         nm_global = row_indices(iiB)

         ! transform the index nm to n and m, formulae copied from Mauro's code
         n_global = MAX(1, nm_global - 1)/nmo + 1
         m_global = nm_global - (n_global - 1)*nmo
         n_global = n_global + homo - gw_corr_lev_occ

         IF (m_global <= homo) THEN

            ! loop over auxiliary basis functions
            DO jjB = 1, ncol_local

               ! Sigma_x_n = -sum_m^occ sum_P (B_(nm)^P)^2
               vec_Sigma_x_gw(n_global, 1) = &
                  vec_Sigma_x_gw(n_global, 1) - &
                  fm_mat_S_gw%local_data(iiB, jjB)**2

            END DO

         END IF

      END DO

      CALL mp_sync(para_env%group)

      CALL mp_sum(vec_Sigma_x_gw, para_env%group)

      vec_Sigma_x_minus_vxc_gw11(:) = &
         vec_Sigma_x_minus_vxc_gw11(:) + &
         vec_Sigma_x_gw(:, 1)

      CALL timestop(handle)

   END SUBROUTINE get_vec_sigma_x

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S_gw_work ...
!> \param vec_W_gw ...
!> \param vec_Sigma_c_gw ...
!> \param vec_omega_fit_gw ...
!> \param vec_Sigma_x_minus_vxc_gw ...
!> \param Eigenval_last ...
!> \param Eigenval_scf ...
!> \param do_periodic ...
!> \param matrix_berry_re_mo_mo ...
!> \param matrix_berry_im_mo_mo ...
!> \param kpoints ...
!> \param do_ri_Sigma_x ...
!> \param vec_Sigma_x_gw ...
!> \param my_do_gw ...
! **************************************************************************************************
   SUBROUTINE deallocate_matrices_gw(fm_mat_S_gw_work, vec_W_gw, vec_Sigma_c_gw, vec_omega_fit_gw, &
                                     vec_Sigma_x_minus_vxc_gw, Eigenval_last, &
                                     Eigenval_scf, do_periodic, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, kpoints, &
                                     do_ri_Sigma_x, vec_Sigma_x_gw, my_do_gw)

      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: fm_mat_S_gw_work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: vec_W_gw
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), INTENT(INOUT)            :: vec_Sigma_c_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: vec_omega_fit_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: vec_Sigma_x_minus_vxc_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: Eigenval_last, Eigenval_scf
      LOGICAL, INTENT(IN)                                :: do_periodic
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      TYPE(kpoint_type), POINTER                         :: kpoints
      LOGICAL, INTENT(IN)                                :: do_ri_Sigma_x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: vec_Sigma_x_gw
      LOGICAL, INTENT(IN)                                :: my_do_gw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_matrices_gw'

      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: my_open_shell

      CALL timeset(routineN, handle)

      nspins = SIZE(Eigenval_last, 2)
      my_open_shell = (nspins == 2)

      IF (my_do_gw) THEN
         DO ispin = 1, nspins
            CALL cp_fm_release(fm_mat_S_gw_work(ispin)%matrix)
         END DO
         DEALLOCATE (fm_mat_S_gw_work)
         DEALLOCATE (vec_Sigma_x_minus_vxc_gw)
         DEALLOCATE (vec_W_gw)
      END IF

      DEALLOCATE (vec_Sigma_c_gw)
      DEALLOCATE (vec_omega_fit_gw)
      DEALLOCATE (Eigenval_last)
      DEALLOCATE (Eigenval_scf)

      IF (do_periodic) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_berry_re_mo_mo)
         CALL dbcsr_deallocate_matrix_set(matrix_berry_im_mo_mo)
         CALL kpoint_release(kpoints)
      END IF
      IF (do_ri_Sigma_x) THEN
         DEALLOCATE (vec_Sigma_x_gw)
      END IF

      CALL timestop(handle)

   END SUBROUTINE deallocate_matrices_gw

! **************************************************************************************************
!> \brief ...
!> \param weights_cos_tf_w_to_t ...
!> \param weights_sin_tf_t_to_w ...
!> \param do_ic_model ...
!> \param do_kpoints_cubic_RPA ...
!> \param fm_mat_W ...
!> \param t_3c_overl_int_ao_mo ...
!> \param t_3c_O_mo_compressed ...
!> \param t_3c_O_mo_ind ...
!> \param t_3c_overl_int_gw_RI ...
!> \param t_3c_overl_int_gw_AO ...
!> \param t_3c_overl_nnP_ic ...
!> \param t_3c_overl_nnP_ic_reflected ...
!> \param mat_W ...
!> \param ikp_local ...
!> \param cfm_mat_W_kp_tau ...
! **************************************************************************************************
   SUBROUTINE deallocate_matrices_gw_im_time(weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, do_ic_model, do_kpoints_cubic_RPA, &
                                             fm_mat_W, &
                                             t_3c_overl_int_ao_mo, t_3c_O_mo_compressed, t_3c_O_mo_ind, &
                                             t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                             t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected, mat_W, &
                                             ikp_local, cfm_mat_W_kp_tau)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w
      LOGICAL, INTENT(IN)                                :: do_ic_model, do_kpoints_cubic_RPA
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: fm_mat_W
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t_3c_overl_int_ao_mo
      TYPE(hfx_compression_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: t_3c_O_mo_compressed
      TYPE(two_dim_int_array), ALLOCATABLE, DIMENSION(:) :: t_3c_O_mo_ind
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: t_3c_overl_int_gw_RI, &
                                                            t_3c_overl_int_gw_AO, &
                                                            t_3c_overl_nnP_ic, &
                                                            t_3c_overl_nnP_ic_reflected
      TYPE(dbcsr_type), POINTER                          :: mat_W
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: ikp_local
      TYPE(cp_cfm_p_type), ALLOCATABLE, &
         DIMENSION(:, :), INTENT(INOUT)                  :: cfm_mat_W_kp_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_matrices_gw_im_time'

      INTEGER                                            :: handle, ikp, ispin, jquad, nspins, unused
      LOGICAL                                            :: my_open_shell

      CALL timeset(routineN, handle)

      nspins = SIZE(t_3c_overl_int_gw_RI)
      my_open_shell = (nspins == 2)

      IF (ALLOCATED(weights_cos_tf_w_to_t)) DEALLOCATE (weights_cos_tf_w_to_t)
      IF (ALLOCATED(weights_sin_tf_t_to_w)) DEALLOCATE (weights_sin_tf_t_to_w)

      IF (.NOT. do_kpoints_cubic_RPA) THEN

         DO jquad = 1, SIZE(fm_mat_W, 1)
            CALL cp_fm_release(fm_mat_W(jquad)%matrix)
         END DO

         DEALLOCATE (fm_mat_W)

         CALL dbcsr_release_P(mat_W)

      ELSE
         DO jquad = 1, SIZE(cfm_mat_W_kp_tau, 2)
            DO ikp = 1, SIZE(cfm_mat_W_kp_tau, 1)
               IF (.NOT. (ANY(ikp_local(:) == ikp))) CYCLE
               CALL cp_cfm_release(cfm_mat_W_kp_tau(ikp, jquad)%matrix)
            END DO
         END DO
         DEALLOCATE (cfm_mat_W_kp_tau)

      END IF

      DO ispin = 1, nspins
         CALL dbcsr_t_destroy(t_3c_overl_int_gw_RI(ispin))
         CALL dbcsr_t_destroy(t_3c_overl_int_gw_AO(ispin))
      END DO
      DEALLOCATE (t_3c_overl_int_gw_AO, t_3c_overl_int_gw_RI)
      IF (do_ic_model) THEN
         DO ispin = 1, nspins
            CALL dbcsr_t_destroy(t_3c_overl_nnP_ic(ispin))
            CALL dbcsr_t_destroy(t_3c_overl_nnP_ic_reflected(ispin))
         END DO
         DEALLOCATE (t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected)
      ENDIF

      DO ispin = 1, nspins
         DEALLOCATE (t_3c_O_mo_ind(ispin)%array)
         CALL dealloc_containers(t_3c_O_mo_compressed(ispin), unused)
      END DO
      DEALLOCATE (t_3c_O_mo_ind, t_3c_O_mo_compressed)

      CALL dbcsr_t_destroy(t_3c_overl_int_ao_mo)

      CALL timestop(handle)

   END SUBROUTINE deallocate_matrices_gw_im_time

! **************************************************************************************************
!> \brief ...
!> \param vec_Sigma_c_gw ...
!> \param dimen_nm_gw ...
!> \param dimen_RI ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param jquad ...
!> \param nmo ...
!> \param num_fit_points ...
!> \param num_integ_points ...
!> \param do_bse ...
!> \param do_im_time ...
!> \param do_periodic ...
!> \param first_cycle_periodic_correction ...
!> \param fermi_level_offset ...
!> \param fermi_level_offset_input ...
!> \param omega ...
!> \param Eigenval ...
!> \param delta_corr ...
!> \param tau_tj ...
!> \param tj ...
!> \param vec_omega_fit_gw ...
!> \param vec_W_gw ...
!> \param wj ...
!> \param weights_cos_tf_w_to_t ...
!> \param fm_mat_W ...
!> \param fm_mat_L ...
!> \param fm_mat_Q ...
!> \param fm_mat_Q_static_bse ...
!> \param fm_mat_R_gw ...
!> \param fm_mat_S_gw ...
!> \param fm_mat_S_gw_work ...
!> \param fm_mat_work ...
!> \param mo_coeff ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param matrix_berry_im_mo_mo ...
!> \param matrix_berry_re_mo_mo ...
!> \param kpoints ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param do_kpoints_cubic_RPA ...
!> \param do_kpoints_from_Gamma ...
! **************************************************************************************************
   SUBROUTINE GW_matrix_operations(vec_Sigma_c_gw, dimen_nm_gw, dimen_RI, gw_corr_lev_occ, &
                                   gw_corr_lev_virt, homo, jquad, nmo, num_fit_points, num_integ_points, &
                                   do_bse, do_im_time, do_periodic, first_cycle_periodic_correction, &
                                   fermi_level_offset, fermi_level_offset_input, &
                                   omega, Eigenval, delta_corr, tau_tj, tj, vec_omega_fit_gw, &
                                   vec_W_gw, wj, weights_cos_tf_w_to_t, fm_mat_W, fm_mat_L, &
                                   fm_mat_Q, fm_mat_Q_static_bse, fm_mat_R_gw, fm_mat_S_gw, &
                                   fm_mat_S_gw_work, fm_mat_work, mo_coeff, para_env, &
                                   para_env_RPA, matrix_berry_im_mo_mo, matrix_berry_re_mo_mo, &
                                   kpoints, qs_env, mp2_env, do_kpoints_cubic_RPA, do_kpoints_from_Gamma)

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), INTENT(INOUT)            :: vec_Sigma_c_gw
      INTEGER, INTENT(IN)                                :: dimen_nm_gw, dimen_RI
      INTEGER, DIMENSION(:), INTENT(IN)                  :: gw_corr_lev_occ, gw_corr_lev_virt, homo
      INTEGER, INTENT(IN)                                :: jquad, nmo, num_fit_points, &
                                                            num_integ_points
      LOGICAL, INTENT(IN)                                :: do_bse, do_im_time, do_periodic
      LOGICAL, INTENT(INOUT) :: first_cycle_periodic_correction
      REAL(KIND=dp), INTENT(INOUT)                       :: fermi_level_offset
      REAL(KIND=dp), INTENT(IN)                          :: fermi_level_offset_input
      REAL(KIND=dp), INTENT(INOUT)                       :: omega
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: delta_corr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tj, vec_omega_fit_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: vec_W_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: wj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: weights_cos_tf_w_to_t
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_W
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_mat_L
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q, fm_mat_Q_static_bse, &
                                                            fm_mat_R_gw
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_S_gw, fm_mat_S_gw_work
      TYPE(cp_fm_type), POINTER                          :: fm_mat_work, mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_im_mo_mo, &
                                                            matrix_berry_re_mo_mo
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA, &
                                                            do_kpoints_from_Gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'GW_matrix_operations'

      INTEGER                                            :: handle, i_global, iiB, iquad, ispin, &
                                                            j_global, jjB, ncol_local, nrow_local, &
                                                            nspins
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_open_shell
      REAL(KIND=dp)                                      :: tau, weight

      CALL timeset(routineN, handle)

      nspins = SIZE(fm_mat_S_gw)
      my_open_shell = (nspins == 2)

      ! Fermi level offset should have a maximum such that the Fermi level of occupied orbitals
      ! is always closer to occupied orbitals than to virtual orbitals and vice versa
      ! that means, the Fermi level offset is at most as big as half the bandgap
      fermi_level_offset = fermi_level_offset_input
      DO ispin = 1, nspins
         fermi_level_offset = MIN(fermi_level_offset, (Eigenval(homo(ispin) + 1, ispin) - Eigenval(homo(ispin), ispin))*0.5_dp)
      END DO

      CALL cp_fm_get_info(matrix=fm_mat_Q, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      IF (.NOT. do_im_time) THEN
         ! calculate [1+Q(iw')]^-1
         CALL cp_fm_cholesky_invert(fm_mat_Q)
         ! symmetrize the result, fm_mat_R_gw is only temporary work matrix
         CALL cp_fm_upper_to_full(fm_mat_Q, fm_mat_R_gw)

         IF (do_bse .AND. jquad == 1) THEN
            CALL cp_fm_to_fm(fm_mat_Q, fm_mat_Q_static_bse)
         END IF

         ! periodic correction for GW
         IF (do_periodic) THEN
            CALL calc_periodic_correction(delta_corr, qs_env, para_env, para_env_RPA, &
                                          mp2_env%ri_g0w0%kp_grid, homo(1), nmo, gw_corr_lev_occ(1), &
                                          gw_corr_lev_virt(1), omega, mo_coeff, Eigenval(:, 1), &
                                          matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                          first_cycle_periodic_correction, kpoints, &
                                          mp2_env%ri_g0w0%do_mo_coeff_gamma, &
                                          mp2_env%ri_g0w0%num_kp_grids, mp2_env%ri_g0w0%eps_kpoint, &
                                          mp2_env%ri_g0w0%do_extra_kpoints, &
                                          mp2_env%ri_g0w0%do_aux_bas_gw, mp2_env%ri_g0w0%frac_aux_mos)
         END IF

         CALL mp_sync(para_env%group)

         ! subtract 1 from the diagonal to get rid of exchange self-energy
!$OMP           PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
!$OMP                       SHARED(ncol_local,nrow_local,col_indices,row_indices,fm_mat_Q,dimen_RI)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)
               IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
                  fm_mat_Q%local_data(iiB, jjB) = fm_mat_Q%local_data(iiB, jjB) - 1.0_dp
               END IF
            END DO
         END DO

         CALL mp_sync(para_env%group)

         DO ispin = 1, nspins
            CALL compute_self_energy_gw(vec_Sigma_c_gw(:, :, :, ispin), dimen_nm_gw, dimen_RI, &
                                        gw_corr_lev_occ(ispin), homo(ispin), jquad, nmo, &
                                        num_fit_points, do_periodic, fermi_level_offset, omega, Eigenval(:, ispin), delta_corr, &
                                        vec_omega_fit_gw, vec_W_gw(:, ispin), wj, fm_mat_Q, &
                                        fm_mat_S_gw(ispin)%matrix, fm_mat_S_gw_work(ispin)%matrix)
         END DO

      END IF ! GW

      ! cubic scaling GW calculation for molecules
      IF (do_im_time .AND. .NOT. (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma)) THEN

         ! calculate [1+Q(iw')]^-1
         CALL cp_fm_cholesky_invert(fm_mat_Q)

         ! symmetrize the result
         CALL cp_fm_upper_to_full(fm_mat_Q, fm_mat_work)

         ! subtract 1 from the diagonal to get rid of exchange self-energy
!$OMP           PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
!$OMP                       SHARED(ncol_local,nrow_local,col_indices,row_indices,fm_mat_Q,dimen_RI)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)
               IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
                  fm_mat_Q%local_data(iiB, jjB) = fm_mat_Q%local_data(iiB, jjB) - 1.0_dp
               END IF
            END DO
         END DO

         ! multiply with L from the left and the right to get the screened Coulomb interaction
         CALL cp_gemm('T', 'N', dimen_RI, dimen_RI, dimen_RI, 1.0_dp, fm_mat_L(1, 1)%matrix, fm_mat_Q, &
                      0.0_dp, fm_mat_work)

         CALL cp_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, 1.0_dp, fm_mat_work, fm_mat_L(1, 1)%matrix, &
                      0.0_dp, fm_mat_Q)

         ! Fourier transform from w to t
         DO iquad = 1, num_integ_points

            omega = tj(jquad)
            tau = tau_tj(iquad)
            weight = weights_cos_tf_w_to_t(iquad, jquad)*COS(tau*omega)

            IF (jquad == 1) THEN

               CALL cp_fm_set_all(matrix=fm_mat_W(iquad)%matrix, alpha=0.0_dp)

            END IF

            CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_mat_W(iquad)%matrix, beta=weight, matrix_b=fm_mat_Q)

         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE GW_matrix_operations

! **************************************************************************************************
!> \brief ...
!> \param vec_Sigma_c_gw ...
!> \param dimen_nm_gw ...
!> \param dimen_RI ...
!> \param gw_corr_lev_occ ...
!> \param homo ...
!> \param jquad ...
!> \param nmo ...
!> \param num_fit_points ...
!> \param do_periodic ...
!> \param fermi_level_offset ...
!> \param omega ...
!> \param Eigenval ...
!> \param delta_corr ...
!> \param vec_omega_fit_gw ...
!> \param vec_W_gw ...
!> \param wj ...
!> \param fm_mat_Q ...
!> \param fm_mat_S_gw ...
!> \param fm_mat_S_gw_work ...
! **************************************************************************************************
   SUBROUTINE compute_self_energy_gw(vec_Sigma_c_gw, dimen_nm_gw, dimen_RI, gw_corr_lev_occ, homo, jquad, nmo, num_fit_points, &
                                     do_periodic, fermi_level_offset, omega, Eigenval, delta_corr, vec_omega_fit_gw, vec_W_gw, &
                                     wj, fm_mat_Q, fm_mat_S_gw, fm_mat_S_gw_work)

      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: vec_Sigma_c_gw
      INTEGER, INTENT(IN)                                :: dimen_nm_gw, dimen_RI, gw_corr_lev_occ, &
                                                            homo, jquad, nmo, num_fit_points
      LOGICAL, INTENT(IN)                                :: do_periodic
      REAL(KIND=dp), INTENT(IN)                          :: fermi_level_offset
      REAL(KIND=dp), INTENT(INOUT)                       :: omega
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: delta_corr, vec_omega_fit_gw
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: vec_W_gw
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wj
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q, fm_mat_S_gw, fm_mat_S_gw_work

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_self_energy_gw'

      INTEGER                                            :: handle, iiB, iquad, jjB, m_global, &
                                                            n_global, ncol_local, nm_global, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: delta_corr_nn, e_fermi, omega_i, &
                                                            sign_occ_virt

      CALL timeset(routineN, handle)

      ! S_work_(nm)Q = B_(nm)P * ([1+Q]^-1-1)_PQ
      CALL cp_gemm(transa="N", transb="N", m=dimen_nm_gw, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                   matrix_a=fm_mat_S_gw, matrix_b=fm_mat_Q, beta=0.0_dp, &
                   matrix_c=fm_mat_S_gw_work)

      CALL cp_fm_get_info(matrix=fm_mat_S_gw, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! vector W_(nm) = S_work_(nm)Q * [B_(nm)Q]^T

      vec_W_gw = 0.0_dp

      DO iiB = 1, nrow_local
         nm_global = row_indices(iiB)
         DO jjB = 1, ncol_local
            vec_W_gw(nm_global) = vec_W_gw(nm_global) + &
                                  fm_mat_S_gw_work%local_data(iiB, jjB)*fm_mat_S_gw%local_data(iiB, jjB)
         END DO

         ! transform the index nm of vec_W_gw back to n and m, formulae copied from Mauro's code
         n_global = MAX(1, nm_global - 1)/nmo + 1
         m_global = nm_global - (n_global - 1)*nmo
         n_global = n_global + homo - gw_corr_lev_occ

         ! compute self-energy for imaginary frequencies
         DO iquad = 1, num_fit_points

            ! for occ orbitals, we compute the self-energy for negative frequencies
            IF (n_global <= homo) THEN
               sign_occ_virt = -1.0_dp
            ELSE
               sign_occ_virt = 1.0_dp
            END IF

            omega_i = vec_omega_fit_gw(iquad)*sign_occ_virt

            ! set the Fermi energy for occ orbitals slightly above the HOMO and
            ! for virt orbitals slightly below the LUMO
            IF (n_global <= homo) THEN
               e_fermi = Eigenval(homo) + fermi_level_offset
            ELSE
               e_fermi = Eigenval(homo + 1) - fermi_level_offset
            END IF

            ! add here the periodic correction
            IF (do_periodic .AND. col_indices(1) == 1 .AND. n_global == m_global) THEN
               delta_corr_nn = delta_corr(n_global)
            ELSE
               delta_corr_nn = 0.0_dp
            END IF

            ! update the self-energy (use that vec_W_gw(iw) is symmetric), divide the integration
            ! weight by 2, because the integration is from -infty to +infty and not just 0 to +infty
            ! as for RPA, also we need for virtual orbitals a complex conjugate
            vec_Sigma_c_gw(n_global - homo + gw_corr_lev_occ, iquad, 1) = &
               vec_Sigma_c_gw(n_global - homo + gw_corr_lev_occ, iquad, 1) - &
               0.5_dp/pi*wj(jquad)/2.0_dp*(vec_W_gw(nm_global) + delta_corr_nn)* &
               (1.0_dp/(gaussi*(omega + omega_i) + e_fermi - Eigenval(m_global)) + &
                1.0_dp/(gaussi*(-omega + omega_i) + e_fermi - Eigenval(m_global)))
         END DO

      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_self_energy_gw

! **************************************************************************************************
!> \brief ...
!> \param vec_Sigma_c_gw ...
!> \param count_ev_sc_GW ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_tot ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param nmo ...
!> \param num_fit_points ...
!> \param num_integ_points ...
!> \param unit_nr ...
!> \param do_apply_ic_corr_to_gw ...
!> \param do_im_time ...
!> \param do_periodic ...
!> \param do_ri_Sigma_x ...
!> \param first_cycle_periodic_correction ...
!> \param e_fermi ...
!> \param eps_filter ...
!> \param fermi_level_offset ...
!> \param delta_corr ...
!> \param Eigenval ...
!> \param Eigenval_last ...
!> \param Eigenval_scf ...
!> \param iter_sc_GW0 ...
!> \param exit_ev_gw ...
!> \param tau_tj ...
!> \param tj ...
!> \param vec_omega_fit_gw ...
!> \param vec_Sigma_x_gw ...
!> \param ic_corr_list ...
!> \param weights_cos_tf_t_to_w ...
!> \param weights_sin_tf_t_to_w ...
!> \param fm_mo_coeff_occ_scaled ...
!> \param fm_mo_coeff_virt_scaled ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param mo_coeff ...
!> \param fm_mat_W ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param mat_dm ...
!> \param mat_SinvVSinv ...
!> \param t_3c_overl_int_ao_mo ...
!> \param t_3c_O_mo_compressed ...
!> \param t_3c_O_mo_ind ...
!> \param t_3c_overl_int_gw_RI ...
!> \param t_3c_overl_int_gw_AO ...
!> \param matrix_berry_im_mo_mo ...
!> \param matrix_berry_re_mo_mo ...
!> \param mat_W ...
!> \param matrix_s ...
!> \param kpoints ...
!> \param mp2_env ...
!> \param qs_env ...
!> \param nkp_self_energy ...
!> \param do_kpoints_cubic_RPA ...
!> \param Eigenval_kp ...
!> \param Eigenval_scf_kp ...
! **************************************************************************************************
   SUBROUTINE compute_QP_energies(vec_Sigma_c_gw, count_ev_sc_GW, gw_corr_lev_occ, &
                                  gw_corr_lev_tot, gw_corr_lev_virt, homo, &
                                  nmo, num_fit_points, num_integ_points, &
                                  unit_nr, do_apply_ic_corr_to_gw, do_im_time, &
                                  do_periodic, do_ri_Sigma_x, &
                                  first_cycle_periodic_correction, e_fermi, eps_filter, &
                                  fermi_level_offset, delta_corr, Eigenval, &
                                  Eigenval_last, Eigenval_scf, iter_sc_GW0, exit_ev_gw, tau_tj, tj, &
                                  vec_omega_fit_gw, vec_Sigma_x_gw, ic_corr_list, &
                                  weights_cos_tf_t_to_w, weights_sin_tf_t_to_w, &
                                  fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, fm_mo_coeff_occ, &
                                  fm_mo_coeff_virt, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, &
                                  mo_coeff, fm_mat_W, para_env, para_env_RPA, mat_dm, mat_SinvVSinv, &
                                  t_3c_overl_int_ao_mo, t_3c_O_mo_compressed, t_3c_O_mo_ind, &
                                  t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, matrix_berry_im_mo_mo, &
                                  matrix_berry_re_mo_mo, mat_W, matrix_s, &
                                  kpoints, mp2_env, qs_env, &
                                  nkp_self_energy, do_kpoints_cubic_RPA, Eigenval_kp, Eigenval_scf_kp)

      COMPLEX(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(OUT)                                     :: vec_Sigma_c_gw
      INTEGER, INTENT(IN)                                :: count_ev_sc_GW
      INTEGER, DIMENSION(:), INTENT(IN)                  :: gw_corr_lev_occ
      INTEGER, INTENT(IN)                                :: gw_corr_lev_tot
      INTEGER, DIMENSION(:), INTENT(IN)                  :: gw_corr_lev_virt, homo
      INTEGER, INTENT(IN)                                :: nmo, num_fit_points, num_integ_points, &
                                                            unit_nr
      LOGICAL, INTENT(IN)                                :: do_apply_ic_corr_to_gw, do_im_time, &
                                                            do_periodic, do_ri_Sigma_x
      LOGICAL, INTENT(INOUT) :: first_cycle_periodic_correction
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: e_fermi
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter, fermi_level_offset
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: delta_corr
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: Eigenval_last, Eigenval_scf
      INTEGER, INTENT(IN)                                :: iter_sc_GW0
      LOGICAL, INTENT(INOUT)                             :: exit_ev_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: tau_tj, tj, vec_omega_fit_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: vec_Sigma_x_gw
      TYPE(one_dim_real_array), DIMENSION(2), INTENT(IN) :: ic_corr_list
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: weights_cos_tf_t_to_w, &
                                                            weights_sin_tf_t_to_w
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff_occ_scaled, &
                                                            fm_mo_coeff_virt_scaled
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mo_coeff_occ, fm_mo_coeff_virt
      TYPE(cp_fm_type), POINTER                          :: fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau, mo_coeff
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_W
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_dm, mat_SinvVSinv
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_ao_mo
      TYPE(hfx_compression_type), DIMENSION(:)           :: t_3c_O_mo_compressed
      TYPE(two_dim_int_array), DIMENSION(:)              :: t_3c_O_mo_ind
      TYPE(dbcsr_t_type), DIMENSION(:)                   :: t_3c_overl_int_gw_RI, &
                                                            t_3c_overl_int_gw_AO
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_im_mo_mo, &
                                                            matrix_berry_re_mo_mo
      TYPE(dbcsr_type), POINTER                          :: mat_W
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: nkp_self_energy
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: Eigenval_kp, Eigenval_scf_kp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_QP_energies'

      INTEGER :: count_ev_sc_GW_print, count_sc_GW0, count_sc_GW0_print, crossing_search, handle, &
         ikp, ispin, n_level_gw, nspins, num_points_corr, num_poles
      LOGICAL                                            :: my_open_shell
      REAL(KIND=dp)                                      :: stop_crit
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: m_value, vec_gw_energ, z_value

      CALL timeset(routineN, handle)

      nspins = SIZE(homo)
      my_open_shell = (nspins == 2)

      DO count_sc_GW0 = 1, iter_sc_GW0

         ! postprocessing for cubic scaling GW calculation
         IF (do_im_time .AND. .NOT. do_kpoints_cubic_RPA) THEN
            num_points_corr = mp2_env%ri_g0w0%num_omega_points

            CALL compute_self_energy_cubic_gw(num_integ_points, nmo, tau_tj, tj, &
                                              matrix_s, fm_mo_coeff_occ(1)%matrix, &
                                              fm_mo_coeff_virt(1)%matrix, fm_mo_coeff_occ_scaled, &
                                              fm_mo_coeff_virt_scaled, fm_scaled_dm_occ_tau, &
                                              fm_scaled_dm_virt_tau, Eigenval(:, 1), eps_filter, &
                                              e_fermi(1), fm_mat_W, &
                                              gw_corr_lev_tot, gw_corr_lev_occ(1), gw_corr_lev_virt(1), homo(1), &
                                              count_ev_sc_GW, count_sc_GW0, &
                                              t_3c_overl_int_ao_mo, t_3c_O_mo_compressed(1), t_3c_O_mo_ind(1)%array, &
                                              t_3c_overl_int_gw_RI(1), t_3c_overl_int_gw_AO(1), &
                                              mat_W, mat_SinvVSinv, mat_dm, &
                                              weights_cos_tf_t_to_w, weights_sin_tf_t_to_w, vec_Sigma_c_gw(:, :, :, 1), &
                                              do_periodic, num_points_corr, delta_corr, qs_env, para_env, para_env_RPA, &
                                              mp2_env, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                              first_cycle_periodic_correction, kpoints, num_fit_points, mo_coeff, &
                                              do_ri_Sigma_x, vec_Sigma_x_gw(:, :, 1), unit_nr, 1)

            IF (my_open_shell) THEN

               CALL compute_self_energy_cubic_gw(num_integ_points, nmo, tau_tj, tj, &
                                                 matrix_s, fm_mo_coeff_occ(2)%matrix, &
                                                 fm_mo_coeff_virt(2)%matrix, fm_mo_coeff_occ_scaled, &
                                                 fm_mo_coeff_virt_scaled, fm_scaled_dm_occ_tau, &
                                                 fm_scaled_dm_virt_tau, Eigenval(:, 2), eps_filter, &
                                                 e_fermi(2), fm_mat_W, &
                                                 gw_corr_lev_tot, gw_corr_lev_occ(2), gw_corr_lev_virt(2), homo(2), &
                                                 count_ev_sc_GW, count_sc_GW0, &
                                                 t_3c_overl_int_ao_mo, t_3c_O_mo_compressed(2), t_3c_O_mo_ind(2)%array, &
                                                 t_3c_overl_int_gw_RI(2), t_3c_overl_int_gw_AO(2), &
                                                 mat_W, mat_SinvVSinv, mat_dm, &
                                                 weights_cos_tf_t_to_w, weights_sin_tf_t_to_w, vec_Sigma_c_gw(:, :, :, 2), &
                                                 do_periodic, num_points_corr, delta_corr, qs_env, para_env, para_env_RPA, &
                                                 mp2_env, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                                 first_cycle_periodic_correction, kpoints, num_fit_points, mo_coeff, &
                                                 do_ri_Sigma_x, vec_Sigma_x_gw(:, :, 2), unit_nr, 2)

            END IF

         END IF

         IF (do_periodic .AND. mp2_env%ri_g0w0%do_average_deg_levels) THEN

            DO ispin = 1, nspins
               CALL average_degenerate_levels(vec_Sigma_c_gw(:, :, :, ispin), &
                                              Eigenval(1 + homo(ispin) - gw_corr_lev_occ(ispin): &
                                                       homo(ispin) + gw_corr_lev_virt(ispin), ispin), &
                                              mp2_env%ri_g0w0%eps_eigenval)
            END DO
         END IF

         IF (.NOT. do_im_time) THEN
            CALL mp_sum(vec_Sigma_c_gw, para_env%group)
         END IF

         CALL mp_sync(para_env%group)

         stop_crit = 1.0e-7
         num_poles = mp2_env%ri_g0w0%num_poles
         crossing_search = mp2_env%ri_g0w0%crossing_search

         ! arrays storing the correlation self-energy, stat. error and z-shot value
         ALLOCATE (vec_gw_energ(gw_corr_lev_tot, nspins))
         vec_gw_energ = 0.0_dp
         ALLOCATE (z_value(gw_corr_lev_tot, nspins))
         z_value = 0.0_dp
         ALLOCATE (m_value(gw_corr_lev_tot, nspins))
         m_value = 0.0_dp

         ! for the normal code for molecules or Gamma only: nkp = 1
         DO ikp = 1, nkp_self_energy

            IF (do_kpoints_cubic_RPA) THEN

               vec_gw_energ = 0.0_dp
               z_value = 0.0_dp
               m_value = 0.0_dp

               CALL get_eigenval_for_conti(Eigenval(:, 1), Eigenval_scf(:, 1), Eigenval_kp, Eigenval_scf_kp, kpoints, &
                                           ikp, my_open_shell)
            END IF

            ! fit the self-energy on imaginary frequency axis and evaluate the fit on the MO energy of the SCF
            DO n_level_gw = 1, gw_corr_lev_tot
               ! processes perform different fits
               IF (MODULO(n_level_gw, para_env%num_pe) /= para_env%mepos) CYCLE

               SELECT CASE (mp2_env%ri_g0w0%analytic_continuation)
               CASE (gw_two_pole_model)
                  CALL fit_and_continuation_2pole(vec_gw_energ(:, 1), vec_omega_fit_gw, &
                                                  z_value(:, 1), m_value(:, 1), vec_Sigma_c_gw(:, :, ikp, 1), &
                                                  mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1, ikp), &
                                                  Eigenval(:, 1), Eigenval_scf(:, 1), n_level_gw, gw_corr_lev_occ(1), num_poles, &
                                                  num_fit_points, crossing_search, homo(1), stop_crit, &
                                                  fermi_level_offset, do_im_time)

               CASE (gw_pade_approx)
                  CALL continuation_pade(vec_gw_energ(:, 1), vec_omega_fit_gw, &
                                         z_value(:, 1), m_value(:, 1), vec_Sigma_c_gw(:, :, ikp, 1), &
                                         mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1, ikp), &
                                         Eigenval(:, 1), Eigenval_scf(:, 1), n_level_gw, &
                                         gw_corr_lev_occ(1), mp2_env%ri_g0w0%nparam_pade, &
                                         num_fit_points, crossing_search, homo(1), fermi_level_offset, &
                                         do_im_time, mp2_env%ri_g0w0%print_self_energy, count_ev_sc_GW)
               CASE DEFAULT
                  CPABORT("Only two-model and Pade approximation are implemented.")
               END SELECT

               IF (my_open_shell) THEN
                  SELECT CASE (mp2_env%ri_g0w0%analytic_continuation)
                  CASE (gw_two_pole_model)
                     CALL fit_and_continuation_2pole( &
                        vec_gw_energ(:, 2), vec_omega_fit_gw, &
                        z_value(:, 2), m_value(:, 2), vec_Sigma_c_gw(:, :, ikp, 2), &
                        mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2, ikp), &
                        Eigenval(:, 2), Eigenval_scf(:, 2), n_level_gw, &
                        gw_corr_lev_occ(2), num_poles, &
                        num_fit_points, crossing_search, homo(2), stop_crit, &
                        fermi_level_offset, do_im_time)
                  CASE (gw_pade_approx)
                     CALL continuation_pade(vec_gw_energ(:, 2), vec_omega_fit_gw, &
                                            z_value(:, 2), m_value(:, 2), vec_Sigma_c_gw(:, :, ikp, 2), &
                                            mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2, ikp), &
                                            Eigenval(:, 2), Eigenval_scf(:, 2), n_level_gw, &
                                            gw_corr_lev_occ(2), mp2_env%ri_g0w0%nparam_pade, &
                                            num_fit_points, crossing_search, homo(2), &
                                            fermi_level_offset, do_im_time, &
                                            mp2_env%ri_g0w0%print_self_energy, count_ev_sc_GW)
                  CASE DEFAULT
                     CPABORT("Only two-model and Pade approximation are implemented.")
                  END SELECT

               END IF

            END DO ! n_level_gw

            CALL mp_sum(vec_gw_energ, para_env%group)
            CALL mp_sum(z_value, para_env%group)
            CALL mp_sum(m_value, para_env%group)

            IF (do_im_time .OR. mp2_env%ri_g0w0%iter_sc_GW0 == 1) THEN
               count_ev_sc_GW_print = count_ev_sc_GW
               count_sc_GW0_print = count_sc_GW0
            ELSE
               count_ev_sc_GW_print = count_sc_GW0
               count_sc_GW0_print = count_ev_sc_GW
            END IF

            ! print the quasiparticle energies and update Eigenval in case you do eigenvalue self-consistent GW
            IF (my_open_shell) THEN

               CALL print_and_update_for_ev_sc( &
                  vec_gw_energ(:, 1), &
                  z_value(:, 1), m_value(:, 1), mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1, ikp), Eigenval(:, 1), &
                  Eigenval_last(:, 1), Eigenval_scf(:, 1), gw_corr_lev_occ(1), gw_corr_lev_tot, &
                  crossing_search, homo(1), unit_nr, count_ev_sc_GW_print, count_sc_GW0_print, &
                  ikp, nkp_self_energy, kpoints, 1)

               CALL print_and_update_for_ev_sc( &
                  vec_gw_energ(:, 2), &
                  z_value(:, 2), m_value(:, 2), mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2, ikp), Eigenval(:, 2), &
                  Eigenval_last(:, 2), Eigenval_scf(:, 2), gw_corr_lev_occ(2), gw_corr_lev_tot, &
                  crossing_search, homo(2), unit_nr, count_ev_sc_GW_print, count_sc_GW0_print, &
                  ikp, nkp_self_energy, kpoints, 2)

               IF (do_apply_ic_corr_to_gw .AND. count_ev_sc_GW == 1) THEN

                  CALL apply_ic_corr(Eigenval(:, 1), Eigenval_scf(:, 1), ic_corr_list(1)%array, &
                                     gw_corr_lev_occ(1), gw_corr_lev_virt(1), gw_corr_lev_tot, &
                                     homo(1), nmo, unit_nr, do_alpha=.TRUE.)

                  CALL apply_ic_corr(Eigenval(:, 2), Eigenval_scf(:, 2), ic_corr_list(2)%array, &
                                     gw_corr_lev_occ(2), gw_corr_lev_virt(2), gw_corr_lev_tot, &
                                     homo(2), nmo, unit_nr, do_beta=.TRUE.)

               END IF

            ELSE

               CALL print_and_update_for_ev_sc( &
                  vec_gw_energ(:, 1), &
                  z_value(:, 1), m_value(:, 1), mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1, ikp), Eigenval(:, 1), &
                  Eigenval_last(:, 1), Eigenval_scf(:, 1), gw_corr_lev_occ(1), gw_corr_lev_tot, &
                  crossing_search, homo(1), unit_nr, count_ev_sc_GW_print, count_sc_GW0_print, &
                  ikp, nkp_self_energy, kpoints, 0)

               IF (do_apply_ic_corr_to_gw .AND. count_ev_sc_GW == 1) THEN

                  CALL apply_ic_corr(Eigenval(:, 1), Eigenval_scf(:, 1), ic_corr_list(1)%array, &
                                     gw_corr_lev_occ(1), gw_corr_lev_virt(1), gw_corr_lev_tot, &
                                     homo(1), nmo, unit_nr)

               END IF

            END IF

         END DO ! ikp

         DEALLOCATE (z_value)
         DEALLOCATE (m_value)
         DEALLOCATE (vec_gw_energ)

         exit_ev_gw = .FALSE.

         ! if HOMO-LUMO gap differs by less than mp2_env%ri_g0w0%eps_sc_iter, exit ev sc GW loop
         IF (ABS(Eigenval(homo(1), 1) - Eigenval_last(homo(1), 1) - Eigenval(homo(1) + 1, 1) + Eigenval_last(homo(1) + 1, 1)) &
             < mp2_env%ri_g0w0%eps_iter) THEN
            IF (count_sc_GW0 == 1) exit_ev_gw = .TRUE.
            EXIT
         END IF

         DO ispin = 1, nspins
            CALL shift_unshifted_levels(Eigenval(:, ispin), Eigenval_last(:, ispin), gw_corr_lev_occ(ispin), &
                                        gw_corr_lev_virt(ispin), homo(ispin), nmo)
         END DO

         ! in case of N^4 scaling GW, the scGW0 cycle is the eigenvalue sc cycle
         IF (.NOT. do_im_time) EXIT

      END DO ! scGW0

      CALL timestop(handle)

   END SUBROUTINE compute_QP_energies

! **************************************************************************************************
!> \brief ...
!> \param delta_corr ...
!> \param qs_env ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param kp_grid ...
!> \param homo ...
!> \param nmo ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param omega ...
!> \param fm_mo_coeff ...
!> \param Eigenval ...
!> \param matrix_berry_re_mo_mo ...
!> \param matrix_berry_im_mo_mo ...
!> \param first_cycle_periodic_correction ...
!> \param kpoints ...
!> \param do_mo_coeff_Gamma_only ...
!> \param num_kp_grids ...
!> \param eps_kpoint ...
!> \param do_extra_kpoints ...
!> \param do_aux_bas ...
!> \param frac_aux_mos ...
! **************************************************************************************************
   SUBROUTINE calc_periodic_correction(delta_corr, qs_env, para_env, para_env_RPA, kp_grid, homo, nmo, &
                                       gw_corr_lev_occ, gw_corr_lev_virt, omega, fm_mo_coeff, Eigenval, &
                                       matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                       first_cycle_periodic_correction, kpoints, do_mo_coeff_Gamma_only, &
                                       num_kp_grids, eps_kpoint, do_extra_kpoints, do_aux_bas, frac_aux_mos)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: delta_corr
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      INTEGER, DIMENSION(:), POINTER                     :: kp_grid
      INTEGER, INTENT(IN)                                :: homo, nmo, gw_corr_lev_occ, &
                                                            gw_corr_lev_virt
      REAL(KIND=dp), INTENT(IN)                          :: omega
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      LOGICAL, INTENT(INOUT) :: first_cycle_periodic_correction
      TYPE(kpoint_type), POINTER                         :: kpoints
      LOGICAL, INTENT(IN)                                :: do_mo_coeff_Gamma_only
      INTEGER, INTENT(IN)                                :: num_kp_grids
      REAL(KIND=dp), INTENT(IN)                          :: eps_kpoint
      LOGICAL, INTENT(IN)                                :: do_extra_kpoints, do_aux_bas
      REAL(KIND=dp), INTENT(IN)                          :: frac_aux_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_periodic_correction'

      INTEGER                                            :: handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eps_head, eps_inv_head
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv

      CALL timeset(routineN, handle)

      IF (first_cycle_periodic_correction) THEN

         CALL get_kpoints(qs_env, kpoints, kp_grid, num_kp_grids, para_env, h_inv, nmo, do_mo_coeff_Gamma_only, &
                          do_extra_kpoints)

         CALL get_berry_phase(qs_env, kpoints, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, fm_mo_coeff, &
                              para_env, do_mo_coeff_Gamma_only, homo, nmo, gw_corr_lev_virt, eps_kpoint, do_aux_bas, &
                              frac_aux_mos)

      END IF

      CALL compute_eps_head_Berry(eps_head, kpoints, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, para_env_RPA, &
                                  qs_env, homo, Eigenval, omega)

      CALL compute_eps_inv_head(eps_inv_head, eps_head, kpoints)

      CALL kpoint_sum_for_eps_inv_head_Berry(delta_corr, eps_inv_head, kpoints, qs_env, &
                                             matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                             homo, gw_corr_lev_occ, gw_corr_lev_virt, para_env_RPA, &
                                             do_extra_kpoints)

      DEALLOCATE (eps_head, eps_inv_head)

      first_cycle_periodic_correction = .FALSE.

      CALL timestop(handle)

   END SUBROUTINE calc_periodic_correction

! **************************************************************************************************
!> \brief ...
!> \param eps_head ...
!> \param kpoints ...
!> \param matrix_berry_re_mo_mo ...
!> \param matrix_berry_im_mo_mo ...
!> \param para_env_RPA ...
!> \param qs_env ...
!> \param homo ...
!> \param Eigenval ...
!> \param omega ...
! **************************************************************************************************
   SUBROUTINE compute_eps_head_Berry(eps_head, kpoints, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, para_env_RPA, &
                                     qs_env, homo, Eigenval, omega)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: eps_head
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_RPA
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: homo
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: omega

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_eps_head_Berry'

      INTEGER :: col, col_end_in_block, col_offset, col_size, handle, i_col, i_row, ikp, nkp, row, &
         row_offset, row_size, row_start_in_block
      REAL(KIND=dp)                                      :: abs_k_square, cell_volume, &
                                                            correct_kpoint(3), cos_square, &
                                                            eigen_diff, relative_kpoint(3), &
                                                            sin_square
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: P_head
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, cell=cell)
      CALL get_cell(cell=cell, deth=cell_volume)

      NULLIFY (data_block)

      nkp = kpoints%nkp

      ALLOCATE (P_head(nkp))
      P_head(:) = 0.0_dp

      ALLOCATE (eps_head(nkp))
      eps_head(:) = 0.0_dp

      DO ikp = 1, nkp

         relative_kpoint(1:3) = MATMUL(cell%hmat, kpoints%xkp(1:3, ikp))

         correct_kpoint(1:3) = twopi*kpoints%xkp(1:3, ikp)

         abs_k_square = (correct_kpoint(1))**2 + (correct_kpoint(2))**2 + (correct_kpoint(3))**2

         ! real part of the Berry phase
         CALL dbcsr_iterator_start(iter, matrix_berry_re_mo_mo(ikp)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            IF (row_offset + row_size <= homo .OR. col_offset > homo) CYCLE

            IF (row_offset <= homo) THEN
               row_start_in_block = homo - row_offset + 2
            ELSE
               row_start_in_block = 1
            END IF

            IF (col_offset + col_size - 1 > homo) THEN
               col_end_in_block = homo - col_offset + 1
            ELSE
               col_end_in_block = col_size
            END IF

            DO i_row = row_start_in_block, row_size

               DO i_col = 1, col_end_in_block

                  eigen_diff = Eigenval(i_col + col_offset - 1) - Eigenval(i_row + row_offset - 1)

                  cos_square = (data_block(i_row, i_col))**2

                  P_head(ikp) = P_head(ikp) + 2.0_dp*eigen_diff/(omega**2 + eigen_diff**2)*cos_square/abs_k_square

               END DO

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter)

         ! imaginary part of the Berry phase
         CALL dbcsr_iterator_start(iter, matrix_berry_im_mo_mo(ikp)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            IF (row_offset + row_size <= homo .OR. col_offset > homo) CYCLE

            IF (row_offset <= homo) THEN
               row_start_in_block = homo - row_offset + 2
            ELSE
               row_start_in_block = 1
            END IF

            IF (col_offset + col_size - 1 > homo) THEN
               col_end_in_block = homo - col_offset + 1
            ELSE
               col_end_in_block = col_size
            END IF

            DO i_row = row_start_in_block, row_size

               DO i_col = 1, col_end_in_block

                  eigen_diff = Eigenval(i_col + col_offset - 1) - Eigenval(i_row + row_offset - 1)

                  sin_square = (data_block(i_row, i_col))**2

                  P_head(ikp) = P_head(ikp) + 2.0_dp*eigen_diff/(omega**2 + eigen_diff**2)*sin_square/abs_k_square

               END DO

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter)

      END DO

      CALL mp_sum(P_head, para_env_RPA%group)

      ! normalize eps_head
      ! 2.0_dp due to closed shell
      eps_head(:) = 1.0_dp - 2.0_dp*P_head(:)/cell_volume*fourpi

      DEALLOCATE (P_head)

      CALL timestop(handle)

   END SUBROUTINE compute_eps_head_Berry

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kpoints ...
!> \param matrix_berry_re_mo_mo ...
!> \param matrix_berry_im_mo_mo ...
!> \param fm_mo_coeff ...
!> \param para_env ...
!> \param do_mo_coeff_Gamma_only ...
!> \param homo ...
!> \param nmo ...
!> \param gw_corr_lev_virt ...
!> \param eps_kpoint ...
!> \param do_aux_bas ...
!> \param frac_aux_mos ...
! **************************************************************************************************
   SUBROUTINE get_berry_phase(qs_env, kpoints, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, fm_mo_coeff, para_env, &
                              do_mo_coeff_Gamma_only, homo, nmo, gw_corr_lev_virt, eps_kpoint, do_aux_bas, &
                              frac_aux_mos)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: do_mo_coeff_Gamma_only
      INTEGER, INTENT(IN)                                :: homo, nmo, gw_corr_lev_virt
      REAL(KIND=dp), INTENT(IN)                          :: eps_kpoint
      LOGICAL, INTENT(IN)                                :: do_aux_bas
      REAL(KIND=dp), INTENT(IN)                          :: frac_aux_mos

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_berry_phase'

      INTEGER                                            :: col_index, handle, i_col_local, ikind, &
                                                            ikp, nao_aux, ncol_local, nkind, nkp, &
                                                            nmo_for_aux_bas
      INTEGER, DIMENSION(:), POINTER                     :: col_indices
      REAL(dp)                                           :: abs_kpoint, correct_kpoint(3), &
                                                            scale_kpoint
      REAL(KIND=dp), DIMENSION(:), POINTER               :: evals_P, evals_P_sqrt_inv
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_aux_aux
      TYPE(cp_fm_type), POINTER :: fm_mat_eigv_P, fm_mat_P, fm_mat_P_sqrt_inv, &
         fm_mat_s_aux_aux_inv, fm_mat_scaled_eigv_P, fm_mat_work_aux_aux
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_s_aux_aux, &
                                                            matrix_s_aux_orb
      TYPE(dbcsr_type), POINTER :: cosmat, cosmat_desymm, mat_mo_coeff_aux, mat_mo_coeff_aux_2, &
         mat_mo_coeff_Gamma_all, mat_mo_coeff_Gamma_occ_and_GW, mat_mo_coeff_im, mat_mo_coeff_re, &
         mat_work_aux_orb, mat_work_aux_orb_2, matrix_P, matrix_P_sqrt, matrix_P_sqrt_inv, &
         matrix_s_inv_aux_aux, sinmat, sinmat_desymm, tmp
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: gw_aux_basis_set_list, orb_basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_gw_aux
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sab_orb_mic, sgwgw_list, &
                                                            sgworb_list
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      nkp = kpoints%nkp

      NULLIFY (matrix_berry_re_mo_mo, matrix_s, cell, matrix_berry_im_mo_mo, sinmat, cosmat, tmp, &
               cosmat_desymm, sinmat_desymm, qs_kind_set, orb_basis_set_list, sab_orb_mic)

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      matrix_s=matrix_s, &
                      qs_kind_set=qs_kind_set, &
                      nkind=nkind, &
                      ks_env=ks_env, &
                      sab_orb=sab_orb)

      ALLOCATE (orb_basis_set_list(nkind))
      CALL basis_set_list_setup(orb_basis_set_list, "ORB", qs_kind_set)

      CALL setup_neighbor_list(sab_orb_mic, orb_basis_set_list, qs_env=qs_env, mic=.FALSE.)

      ! create dbcsr matrix of mo_coeff for multiplcation
      NULLIFY (mat_mo_coeff_re)
      CALL dbcsr_init_p(mat_mo_coeff_re)
      CALL dbcsr_create(matrix=mat_mo_coeff_re, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      NULLIFY (mat_mo_coeff_im)
      CALL dbcsr_init_p(mat_mo_coeff_im)
      CALL dbcsr_create(matrix=mat_mo_coeff_im, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      NULLIFY (mat_mo_coeff_Gamma_all)
      CALL dbcsr_init_p(mat_mo_coeff_Gamma_all)
      CALL dbcsr_create(matrix=mat_mo_coeff_Gamma_all, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL copy_fm_to_dbcsr(fm_mo_coeff, mat_mo_coeff_Gamma_all, keep_sparsity=.FALSE.)

      NULLIFY (mat_mo_coeff_Gamma_occ_and_GW)
      CALL dbcsr_init_p(mat_mo_coeff_Gamma_occ_and_GW)
      CALL dbcsr_create(matrix=mat_mo_coeff_Gamma_occ_and_GW, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL copy_fm_to_dbcsr(fm_mo_coeff, mat_mo_coeff_Gamma_occ_and_GW, keep_sparsity=.FALSE.)

      IF (.NOT. do_aux_bas) THEN

         ! allocate intermediate matrices
         CALL dbcsr_init_p(cosmat)
         CALL dbcsr_init_p(sinmat)
         CALL dbcsr_init_p(tmp)
         CALL dbcsr_init_p(cosmat_desymm)
         CALL dbcsr_init_p(sinmat_desymm)
         CALL dbcsr_create(matrix=cosmat, template=matrix_s(1)%matrix)
         CALL dbcsr_create(matrix=sinmat, template=matrix_s(1)%matrix)
         CALL dbcsr_create(matrix=tmp, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix=cosmat_desymm, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix=sinmat_desymm, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(cosmat, matrix_s(1)%matrix)
         CALL dbcsr_copy(sinmat, matrix_s(1)%matrix)
         CALL dbcsr_set(cosmat, 0.0_dp)
         CALL dbcsr_set(sinmat, 0.0_dp)

         CALL dbcsr_allocate_matrix_set(matrix_berry_re_mo_mo, nkp)
         CALL dbcsr_allocate_matrix_set(matrix_berry_im_mo_mo, nkp)

      ELSE

         NULLIFY (gw_aux_basis_set_list)
         ALLOCATE (gw_aux_basis_set_list(nkind))

         DO ikind = 1, nkind

            NULLIFY (gw_aux_basis_set_list(ikind)%gto_basis_set)

            NULLIFY (basis_set_gw_aux)

            qs_kind => qs_kind_set(ikind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_gw_aux, basis_type="AUX_GW")
            CPASSERT(ASSOCIATED(basis_set_gw_aux))

            basis_set_gw_aux%kind_radius = orb_basis_set_list(ikind)%gto_basis_set%kind_radius

            gw_aux_basis_set_list(ikind)%gto_basis_set => basis_set_gw_aux

         END DO

         ! neighbor lists
         NULLIFY (sgwgw_list, sgworb_list)
         CALL setup_neighbor_list(sgwgw_list, gw_aux_basis_set_list, qs_env=qs_env)
         CALL setup_neighbor_list(sgworb_list, gw_aux_basis_set_list, orb_basis_set_list, qs_env=qs_env)

         NULLIFY (matrix_s_aux_aux, matrix_s_aux_orb)

         ! build overlap matrix in gw aux basis and the mixed gw aux basis-orb basis
         CALL build_overlap_matrix_simple(ks_env, matrix_s_aux_aux, &
                                          gw_aux_basis_set_list, gw_aux_basis_set_list, sgwgw_list)

         CALL build_overlap_matrix_simple(ks_env, matrix_s_aux_orb, &
                                          gw_aux_basis_set_list, orb_basis_set_list, sgworb_list)

         CALL dbcsr_get_info(matrix_s_aux_aux(1)%matrix, nfullrows_total=nao_aux)

         nmo_for_aux_bas = FLOOR(frac_aux_mos*REAL(nao_aux, KIND=dp))

         CALL cp_fm_struct_create(fm_struct_aux_aux, &
                                  context=fm_mo_coeff%matrix_struct%context, &
                                  nrow_global=nao_aux, &
                                  ncol_global=nao_aux, &
                                  para_env=para_env)

         NULLIFY (mat_work_aux_orb)
         CALL dbcsr_init_p(mat_work_aux_orb)
         CALL dbcsr_create(matrix=mat_work_aux_orb, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (mat_work_aux_orb_2)
         CALL dbcsr_init_p(mat_work_aux_orb_2)
         CALL dbcsr_create(matrix=mat_work_aux_orb_2, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (mat_mo_coeff_aux)
         CALL dbcsr_init_p(mat_mo_coeff_aux)
         CALL dbcsr_create(matrix=mat_mo_coeff_aux, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (mat_mo_coeff_aux_2)
         CALL dbcsr_init_p(mat_mo_coeff_aux_2)
         CALL dbcsr_create(matrix=mat_mo_coeff_aux_2, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (matrix_s_inv_aux_aux)
         CALL dbcsr_init_p(matrix_s_inv_aux_aux)
         CALL dbcsr_create(matrix=matrix_s_inv_aux_aux, &
                           template=matrix_s_aux_aux(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (matrix_P)
         CALL dbcsr_init_p(matrix_P)
         CALL dbcsr_create(matrix=matrix_P, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (matrix_P_sqrt)
         CALL dbcsr_init_p(matrix_P_sqrt)
         CALL dbcsr_create(matrix=matrix_P_sqrt, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (matrix_P_sqrt_inv)
         CALL dbcsr_init_p(matrix_P_sqrt_inv)
         CALL dbcsr_create(matrix=matrix_P_sqrt_inv, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (fm_mat_s_aux_aux_inv)
         CALL cp_fm_create(fm_mat_s_aux_aux_inv, fm_struct_aux_aux, name="inverse overlap mat")

         NULLIFY (fm_mat_work_aux_aux)
         CALL cp_fm_create(fm_mat_work_aux_aux, fm_struct_aux_aux, name="work mat")

         NULLIFY (fm_mat_P)
         CALL cp_fm_create(fm_mat_P, fm_mo_coeff%matrix_struct)

         NULLIFY (fm_mat_eigv_P)
         CALL cp_fm_create(fm_mat_eigv_P, fm_mo_coeff%matrix_struct)

         NULLIFY (fm_mat_scaled_eigv_P)
         CALL cp_fm_create(fm_mat_scaled_eigv_P, fm_mo_coeff%matrix_struct)

         NULLIFY (fm_mat_P_sqrt_inv)
         CALL cp_fm_create(fm_mat_P_sqrt_inv, fm_mo_coeff%matrix_struct)

         NULLIFY (evals_P)
         ALLOCATE (evals_P(nmo))

         NULLIFY (evals_P_sqrt_inv)
         ALLOCATE (evals_P_sqrt_inv(nmo))

         CALL copy_dbcsr_to_fm(matrix_s_aux_aux(1)%matrix, fm_mat_s_aux_aux_inv)
         ! Calculate S_inverse
         CALL cp_fm_cholesky_decompose(fm_mat_s_aux_aux_inv)
         CALL cp_fm_cholesky_invert(fm_mat_s_aux_aux_inv)
         ! Symmetrize the guy
         CALL cp_fm_upper_to_full(fm_mat_s_aux_aux_inv, fm_mat_work_aux_aux)

         CALL copy_fm_to_dbcsr(fm_mat_s_aux_aux_inv, matrix_s_inv_aux_aux, keep_sparsity=.FALSE.)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s_inv_aux_aux, matrix_s_aux_orb(1)%matrix, 0.0_dp, mat_work_aux_orb, &
                             filter_eps=1.0E-15_dp)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, mat_work_aux_orb, mat_mo_coeff_Gamma_all, 0.0_dp, mat_mo_coeff_aux_2, &
                             last_column=nmo_for_aux_bas, filter_eps=1.0E-15_dp)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s_aux_aux(1)%matrix, mat_mo_coeff_aux_2, 0.0_dp, mat_work_aux_orb, &
                             filter_eps=1.0E-15_dp)

         CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_mo_coeff_aux_2, mat_work_aux_orb, 0.0_dp, matrix_P, &
                             filter_eps=1.0E-15_dp)

         CALL copy_dbcsr_to_fm(matrix_P, fm_mat_P)

         CALL cp_fm_syevd(fm_mat_P, fm_mat_eigv_P, evals_P)

         ! only invert the eigenvalues which correspond to the MOs used in the aux. basis
         evals_P_sqrt_inv(1:nmo - nmo_for_aux_bas) = 0.0_dp
         evals_P_sqrt_inv(nmo - nmo_for_aux_bas + 1:nmo) = 1.0_dp/SQRT(evals_P(nmo - nmo_for_aux_bas + 1:nmo))

         CALL cp_fm_to_fm(fm_mat_eigv_P, fm_mat_scaled_eigv_P)

         CALL cp_fm_get_info(matrix=fm_mat_scaled_eigv_P, &
                             ncol_local=ncol_local, &
                             col_indices=col_indices)

         CALL mp_sync(para_env%group)

         ! multiply eigenvectors with inverse sqrt of eigenvalues
         DO i_col_local = 1, ncol_local

            col_index = col_indices(i_col_local)

            fm_mat_scaled_eigv_P%local_data(:, i_col_local) = &
               fm_mat_scaled_eigv_P%local_data(:, i_col_local)*evals_P_sqrt_inv(col_index)

         END DO

         CALL mp_sync(para_env%group)

         CALL cp_gemm(transa="N", transb="T", m=nmo, n=nmo, k=nmo, alpha=1.0_dp, &
                      matrix_a=fm_mat_eigv_P, matrix_b=fm_mat_scaled_eigv_P, beta=0.0_dp, &
                      matrix_c=fm_mat_P_sqrt_inv)

         CALL copy_fm_to_dbcsr(fm_mat_P_sqrt_inv, matrix_P_sqrt_inv, keep_sparsity=.FALSE.)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, mat_mo_coeff_aux_2, matrix_P_sqrt_inv, 0.0_dp, mat_mo_coeff_aux, &
                             filter_eps=1.0E-15_dp)

         ! allocate intermediate matrices
         CALL dbcsr_init_p(cosmat)
         CALL dbcsr_init_p(sinmat)
         CALL dbcsr_init_p(tmp)
         CALL dbcsr_init_p(cosmat_desymm)
         CALL dbcsr_init_p(sinmat_desymm)
         CALL dbcsr_create(matrix=cosmat, template=matrix_s_aux_aux(1)%matrix)
         CALL dbcsr_create(matrix=sinmat, template=matrix_s_aux_aux(1)%matrix)
         CALL dbcsr_create(matrix=tmp, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix=cosmat_desymm, &
                           template=matrix_s_aux_aux(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix=sinmat_desymm, &
                           template=matrix_s_aux_aux(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(cosmat, matrix_s_aux_aux(1)%matrix)
         CALL dbcsr_copy(sinmat, matrix_s_aux_aux(1)%matrix)
         CALL dbcsr_set(cosmat, 0.0_dp)
         CALL dbcsr_set(sinmat, 0.0_dp)

         CALL dbcsr_allocate_matrix_set(matrix_berry_re_mo_mo, nkp)
         CALL dbcsr_allocate_matrix_set(matrix_berry_im_mo_mo, nkp)

         ! allocate the new MO coefficients in the aux basis
         CALL dbcsr_release_p(mat_mo_coeff_Gamma_all)
         CALL dbcsr_release_p(mat_mo_coeff_Gamma_occ_and_GW)

         NULLIFY (mat_mo_coeff_Gamma_all)
         CALL dbcsr_init_p(mat_mo_coeff_Gamma_all)
         CALL dbcsr_create(matrix=mat_mo_coeff_Gamma_all, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_copy(mat_mo_coeff_Gamma_all, mat_mo_coeff_aux)

         NULLIFY (mat_mo_coeff_Gamma_occ_and_GW)
         CALL dbcsr_init_p(mat_mo_coeff_Gamma_occ_and_GW)
         CALL dbcsr_create(matrix=mat_mo_coeff_Gamma_occ_and_GW, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_copy(mat_mo_coeff_Gamma_occ_and_GW, mat_mo_coeff_aux)

         DEALLOCATE (evals_P, evals_P_sqrt_inv)

      END IF

      CALL remove_unnecessary_blocks(mat_mo_coeff_Gamma_occ_and_GW, homo, gw_corr_lev_virt)

      DO ikp = 1, nkp

         ALLOCATE (matrix_berry_re_mo_mo(ikp)%matrix)
         CALL dbcsr_init_p(matrix_berry_re_mo_mo(ikp)%matrix)
         CALL dbcsr_create(matrix_berry_re_mo_mo(ikp)%matrix, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_desymmetrize(matrix_s(1)%matrix, matrix_berry_re_mo_mo(ikp)%matrix)
         CALL dbcsr_set(matrix_berry_re_mo_mo(ikp)%matrix, 0.0_dp)

         ALLOCATE (matrix_berry_im_mo_mo(ikp)%matrix)
         CALL dbcsr_init_p(matrix_berry_im_mo_mo(ikp)%matrix)
         CALL dbcsr_create(matrix_berry_im_mo_mo(ikp)%matrix, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_desymmetrize(matrix_s(1)%matrix, matrix_berry_im_mo_mo(ikp)%matrix)
         CALL dbcsr_set(matrix_berry_im_mo_mo(ikp)%matrix, 0.0_dp)

         correct_kpoint(1:3) = -twopi*kpoints%xkp(1:3, ikp)

         abs_kpoint = SQRT(correct_kpoint(1)**2 + correct_kpoint(2)**2 + correct_kpoint(3)**2)

         IF (abs_kpoint < eps_kpoint) THEN

            scale_kpoint = eps_kpoint/abs_kpoint
            correct_kpoint(:) = correct_kpoint(:)*scale_kpoint

         END IF

         ! get the Berry phase
         IF (do_aux_bas) THEN
            CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, correct_kpoint, sab_orb_external=sab_orb_mic, &
                                           basis_type="AUX_GW")
         ELSE
            CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, correct_kpoint, sab_orb_external=sab_orb_mic, &
                                           basis_type="ORB")
         END IF

         IF (do_mo_coeff_Gamma_only) THEN

            CALL dbcsr_desymmetrize(cosmat, cosmat_desymm)

            CALL dbcsr_multiply('N', 'N', 1.0_dp, cosmat_desymm, mat_mo_coeff_Gamma_occ_and_GW, 0.0_dp, tmp, &
                                filter_eps=1.0E-15_dp)

            CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_mo_coeff_Gamma_all, tmp, 0.0_dp, &
                                matrix_berry_re_mo_mo(ikp)%matrix, filter_eps=1.0E-15_dp)

            CALL dbcsr_desymmetrize(sinmat, sinmat_desymm)

            CALL dbcsr_multiply('N', 'N', 1.0_dp, sinmat_desymm, mat_mo_coeff_Gamma_occ_and_GW, 0.0_dp, tmp, &
                                filter_eps=1.0E-15_dp)

            CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_mo_coeff_Gamma_all, tmp, 0.0_dp, &
                                matrix_berry_im_mo_mo(ikp)%matrix, filter_eps=1.0E-15_dp)

         ELSE

            ! get mo coeff at the ikp
            CALL copy_fm_to_dbcsr(kpoints%kp_env(ikp)%kpoint_env%mos(1, 1)%mo_set%mo_coeff, &
                                  mat_mo_coeff_re, keep_sparsity=.FALSE.)

            CALL copy_fm_to_dbcsr(kpoints%kp_env(ikp)%kpoint_env%mos(2, 1)%mo_set%mo_coeff, &
                                  mat_mo_coeff_im, keep_sparsity=.FALSE.)

            CALL dbcsr_desymmetrize(cosmat, cosmat_desymm)

            CALL dbcsr_desymmetrize(sinmat, sinmat_desymm)

            ! I.
            CALL dbcsr_multiply('N', 'N', 1.0_dp, cosmat_desymm, mat_mo_coeff_re, 0.0_dp, tmp)

            ! I.1
            CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_mo_coeff_Gamma_all, tmp, 0.0_dp, &
                                matrix_berry_re_mo_mo(ikp)%matrix)

            ! II.
            CALL dbcsr_multiply('N', 'N', 1.0_dp, sinmat_desymm, mat_mo_coeff_re, 0.0_dp, tmp)

            ! II.5
            CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_mo_coeff_Gamma_all, tmp, 0.0_dp, &
                                matrix_berry_im_mo_mo(ikp)%matrix)

            ! III.
            CALL dbcsr_multiply('N', 'N', 1.0_dp, cosmat_desymm, mat_mo_coeff_im, 0.0_dp, tmp)

            ! III.7
            CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_mo_coeff_Gamma_all, tmp, 1.0_dp, &
                                matrix_berry_im_mo_mo(ikp)%matrix)

            ! IV.
            CALL dbcsr_multiply('N', 'N', 1.0_dp, sinmat_desymm, mat_mo_coeff_im, 0.0_dp, tmp)

            ! IV.3
            CALL dbcsr_multiply('T', 'N', -1.0_dp, mat_mo_coeff_Gamma_all, tmp, 1.0_dp, &
                                matrix_berry_re_mo_mo(ikp)%matrix)

         END IF

         IF (abs_kpoint < eps_kpoint) THEN

            CALL dbcsr_scale(matrix_berry_im_mo_mo(ikp)%matrix, 1.0_dp/scale_kpoint)
            CALL dbcsr_set(matrix_berry_re_mo_mo(ikp)%matrix, 0.0_dp)
            CALL dbcsr_add_on_diag(matrix_berry_re_mo_mo(ikp)%matrix, 1.0_dp)

         END IF

      END DO

      CALL dbcsr_release_p(cosmat)
      CALL dbcsr_release_p(sinmat)
      CALL dbcsr_release_p(mat_mo_coeff_re)
      CALL dbcsr_release_p(mat_mo_coeff_im)
      CALL dbcsr_release_p(mat_mo_coeff_Gamma_all)
      CALL dbcsr_release_p(mat_mo_coeff_Gamma_occ_and_GW)
      CALL dbcsr_release_p(tmp)
      CALL dbcsr_release_p(cosmat_desymm)
      CALL dbcsr_release_p(sinmat_desymm)
      DEALLOCATE (orb_basis_set_list)

      CALL release_neighbor_list_sets(sab_orb_mic)

      IF (do_aux_bas) THEN

         DEALLOCATE (gw_aux_basis_set_list)
         CALL dbcsr_deallocate_matrix_set(matrix_s_aux_aux)
         CALL dbcsr_deallocate_matrix_set(matrix_s_aux_orb)
         CALL dbcsr_release_p(mat_work_aux_orb)
         CALL dbcsr_release_p(mat_work_aux_orb_2)
         CALL dbcsr_release_p(mat_mo_coeff_aux)
         CALL dbcsr_release_p(mat_mo_coeff_aux_2)
         CALL dbcsr_release_p(matrix_s_inv_aux_aux)
         CALL dbcsr_release_p(matrix_P)
         CALL dbcsr_release_p(matrix_P_sqrt)
         CALL dbcsr_release_p(matrix_P_sqrt_inv)

         CALL cp_fm_struct_release(fm_struct_aux_aux)

         CALL cp_fm_release(fm_mat_s_aux_aux_inv)
         CALL cp_fm_release(fm_mat_work_aux_aux)
         CALL cp_fm_release(fm_mat_P)
         CALL cp_fm_release(fm_mat_eigv_P)
         CALL cp_fm_release(fm_mat_scaled_eigv_P)
         CALL cp_fm_release(fm_mat_P_sqrt_inv)

         ! Deallocate the neighbor list structure
         CALL release_neighbor_list_sets(sgwgw_list)
         CALL release_neighbor_list_sets(sgworb_list)

      END IF

      CALL timestop(handle)

   END SUBROUTINE get_berry_phase

! **************************************************************************************************
!> \brief ...
!> \param mat_mo_coeff_Gamma_occ_and_GW ...
!> \param homo ...
!> \param gw_corr_lev_virt ...
! **************************************************************************************************
   SUBROUTINE remove_unnecessary_blocks(mat_mo_coeff_Gamma_occ_and_GW, homo, gw_corr_lev_virt)

      TYPE(dbcsr_type), POINTER                          :: mat_mo_coeff_Gamma_occ_and_GW
      INTEGER, INTENT(IN)                                :: homo, gw_corr_lev_virt

      INTEGER                                            :: col, col_offset, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_iterator_start(iter, mat_mo_coeff_Gamma_occ_and_GW)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        col_offset=col_offset)

         IF (col_offset > homo + gw_corr_lev_virt) THEN

            data_block = 0.0_dp

         END IF

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_filter(mat_mo_coeff_Gamma_occ_and_GW, 1.0E-15_dp)

   END SUBROUTINE remove_unnecessary_blocks

! **************************************************************************************************
!> \brief ...
!> \param delta_corr ...
!> \param eps_inv_head ...
!> \param kpoints ...
!> \param qs_env ...
!> \param matrix_berry_re_mo_mo ...
!> \param matrix_berry_im_mo_mo ...
!> \param homo ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param para_env_RPA ...
!> \param do_extra_kpoints ...
! **************************************************************************************************
   SUBROUTINE kpoint_sum_for_eps_inv_head_Berry(delta_corr, eps_inv_head, kpoints, qs_env, matrix_berry_re_mo_mo, &
                                                matrix_berry_im_mo_mo, homo, gw_corr_lev_occ, gw_corr_lev_virt, &
                                                para_env_RPA, do_extra_kpoints)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: delta_corr
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eps_inv_head
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      INTEGER, INTENT(IN)                                :: homo, gw_corr_lev_occ, gw_corr_lev_virt
      TYPE(cp_para_env_type), INTENT(IN), OPTIONAL       :: para_env_RPA
      LOGICAL, INTENT(IN)                                :: do_extra_kpoints

      INTEGER                                            :: col, col_offset, col_size, i_col, i_row, &
                                                            ikp, m_level, n_level_gw, nkp, row, &
                                                            row_offset, row_size
      REAL(KIND=dp)                                      :: abs_k_square, cell_volume, &
                                                            check_int_one_over_ksq, contribution, &
                                                            weight
      REAL(KIND=dp), DIMENSION(3)                        :: correct_kpoint
      REAL(KIND=dp), DIMENSION(:), POINTER               :: delta_corr_extra
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_iterator_type)                          :: iter, iter_new

      CALL get_qs_env(qs_env=qs_env, cell=cell)

      CALL get_cell(cell=cell, deth=cell_volume)

      nkp = kpoints%nkp

      delta_corr = 0.0_dp

      IF (do_extra_kpoints) THEN
         NULLIFY (delta_corr_extra)
         ALLOCATE (delta_corr_extra(1 + homo - gw_corr_lev_occ:homo + gw_corr_lev_virt))
         delta_corr_extra = 0.0_dp
      END IF

      check_int_one_over_ksq = 0.0_dp

      DO ikp = 1, nkp

         weight = kpoints%wkp(ikp)

         correct_kpoint(1:3) = twopi*kpoints%xkp(1:3, ikp)

         abs_k_square = (correct_kpoint(1))**2 + (correct_kpoint(2))**2 + (correct_kpoint(3))**2

         ! cos part of the Berry phase
         CALL dbcsr_iterator_start(iter, matrix_berry_re_mo_mo(ikp)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            DO i_col = 1, col_size

               DO n_level_gw = 1 + homo - gw_corr_lev_occ, homo + gw_corr_lev_virt

                  IF (n_level_gw == i_col + col_offset - 1) THEN

                     DO i_row = 1, row_size

                        contribution = weight*(eps_inv_head(ikp) - 1.0_dp)/abs_k_square*(data_block(i_row, i_col))**2

                        m_level = i_row + row_offset - 1

                        ! we only compute the correction for n=m
                        IF (m_level .NE. n_level_gw) CYCLE

                        IF (.NOT. do_extra_kpoints) THEN

                           delta_corr(n_level_gw) = delta_corr(n_level_gw) + contribution

                        ELSE

                           IF (ikp <= nkp*8/9) THEN

                              delta_corr(n_level_gw) = delta_corr(n_level_gw) + contribution

                           ELSE

                              delta_corr_extra(n_level_gw) = delta_corr_extra(n_level_gw) + contribution

                           END IF

                        END IF

                     END DO

                  END IF

               END DO

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter)

         ! the same for the im. part of the Berry phase
         CALL dbcsr_iterator_start(iter_new, matrix_berry_im_mo_mo(ikp)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter_new))

            CALL dbcsr_iterator_next_block(iter_new, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            DO i_col = 1, col_size

               DO n_level_gw = 1 + homo - gw_corr_lev_occ, homo + gw_corr_lev_virt

                  IF (n_level_gw == i_col + col_offset - 1) THEN

                     DO i_row = 1, row_size

                        m_level = i_row + row_offset - 1

                        contribution = weight*(eps_inv_head(ikp) - 1.0_dp)/abs_k_square*(data_block(i_row, i_col))**2

                        ! we only compute the correction for n=m
                        IF (m_level .NE. n_level_gw) CYCLE

                        IF (.NOT. do_extra_kpoints) THEN

                           delta_corr(n_level_gw) = delta_corr(n_level_gw) + contribution

                        ELSE

                           IF (ikp <= nkp*8/9) THEN

                              delta_corr(n_level_gw) = delta_corr(n_level_gw) + contribution

                           ELSE

                              delta_corr_extra(n_level_gw) = delta_corr_extra(n_level_gw) + contribution

                           END IF

                        END IF

                     END DO

                  END IF

               END DO

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter_new)

         check_int_one_over_ksq = check_int_one_over_ksq + weight/abs_k_square

      END DO

      ! normalize by the cell volume
      delta_corr = delta_corr/cell_volume*fourpi

      check_int_one_over_ksq = check_int_one_over_ksq/cell_volume

      CALL mp_sum(delta_corr, para_env_RPA%group)

      IF (do_extra_kpoints) THEN

         delta_corr_extra = delta_corr_extra/cell_volume*fourpi

         CALL mp_sum(delta_corr_extra, para_env_RPA%group)

         delta_corr(:) = delta_corr(:) + (delta_corr(:) - delta_corr_extra(:))

         DEALLOCATE (delta_corr_extra)

      END IF

   END SUBROUTINE kpoint_sum_for_eps_inv_head_Berry

! **************************************************************************************************
!> \brief ...
!> \param eps_inv_head ...
!> \param eps_head ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE compute_eps_inv_head(eps_inv_head, eps_head, kpoints)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: eps_inv_head
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eps_head
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_eps_inv_head'

      INTEGER                                            :: handle, ikp, nkp

      CALL timeset(routineN, handle)

      nkp = kpoints%nkp

      ALLOCATE (eps_inv_head(nkp))

      DO ikp = 1, nkp

         eps_inv_head(ikp) = 1.0_dp/eps_head(ikp)

      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_eps_inv_head

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kpoints ...
!> \param kp_grid ...
!> \param num_kp_grids ...
!> \param para_env ...
!> \param h_inv ...
!> \param nmo ...
!> \param do_mo_coeff_Gamma_only ...
!> \param do_extra_kpoints ...
! **************************************************************************************************
   SUBROUTINE get_kpoints(qs_env, kpoints, kp_grid, num_kp_grids, para_env, h_inv, nmo, &
                          do_mo_coeff_Gamma_only, do_extra_kpoints)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, DIMENSION(:), POINTER                     :: kp_grid
      INTEGER, INTENT(IN)                                :: num_kp_grids
      TYPE(cp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: h_inv
      INTEGER, INTENT(IN)                                :: nmo
      LOGICAL, INTENT(IN)                                :: do_mo_coeff_Gamma_only, do_extra_kpoints

      INTEGER                                            :: end_kp, i, i_grid_level, ix, iy, iz, &
                                                            nkp_inner_grid, nkp_outer_grid, &
                                                            npoints, start_kp
      INTEGER, DIMENSION(3)                              :: outer_kp_grid
      REAL(KIND=dp)                                      :: kpoint_weight_left, single_weight
      REAL(KIND=dp), DIMENSION(3)                        :: kpt_latt, reducing_factor
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_environment_type), POINTER                 :: qs_env_kp_Gamma_only

      NULLIFY (kpoints, cell, particle_set, qs_env_kp_Gamma_only)

      ! check whether kp_grid includes the Gamma point. If so, abort.
      CPASSERT(MOD(kp_grid(1)*kp_grid(2)*kp_grid(3), 2) == 0)
      IF (do_extra_kpoints) THEN
         CPASSERT(do_mo_coeff_Gamma_only)
      END IF

      IF (do_mo_coeff_Gamma_only) THEN

         outer_kp_grid(1) = kp_grid(1) - 1
         outer_kp_grid(2) = kp_grid(2) - 1
         outer_kp_grid(3) = kp_grid(3) - 1

         CALL get_qs_env(qs_env=qs_env, cell=cell, particle_set=particle_set)

         CALL get_cell(cell, h_inv=h_inv)

         CALL kpoint_create(kpoints)

         kpoints%kp_scheme = "GENERAL"
         kpoints%symmetry = .FALSE.
         kpoints%verbose = .FALSE.
         kpoints%full_grid = .FALSE.
         kpoints%use_real_wfn = .FALSE.
         kpoints%eps_geo = 1.e-6_dp
         npoints = kp_grid(1)*kp_grid(2)*kp_grid(3)/2 + &
                   (num_kp_grids - 1)*((outer_kp_grid(1) + 1)/2*outer_kp_grid(2)*outer_kp_grid(3) - 1)

         IF (do_extra_kpoints) THEN

            CPASSERT(num_kp_grids == 1)
            CPASSERT(MOD(kp_grid(1), 4) == 0)
            CPASSERT(MOD(kp_grid(2), 4) == 0)
            CPASSERT(MOD(kp_grid(3), 4) == 0)

         END IF

         IF (do_extra_kpoints) THEN

            npoints = kp_grid(1)*kp_grid(2)*kp_grid(3)/2 + kp_grid(1)*kp_grid(2)*kp_grid(3)/2/8

         END IF

         kpoints%full_grid = .TRUE.
         kpoints%nkp = npoints
         ALLOCATE (kpoints%xkp(3, npoints), kpoints%wkp(npoints))
         kpoints%xkp = 0.0_dp
         kpoints%wkp = 0.0_dp

         nkp_outer_grid = outer_kp_grid(1)*outer_kp_grid(2)*outer_kp_grid(3)
         nkp_inner_grid = kp_grid(1)*kp_grid(2)*kp_grid(3)

         i = 0
         reducing_factor(:) = 1.0_dp
         kpoint_weight_left = 1.0_dp

         ! the outer grids
         DO i_grid_level = 1, num_kp_grids - 1

            single_weight = kpoint_weight_left/REAL(nkp_outer_grid, KIND=dp)

            start_kp = i + 1

            DO ix = 1, outer_kp_grid(1)
               DO iy = 1, outer_kp_grid(2)
                  DO iz = 1, outer_kp_grid(3)

                     ! exclude Gamma
                     IF (2*ix - outer_kp_grid(1) - 1 == 0 .AND. 2*iy - outer_kp_grid(2) - 1 == 0 .AND. &
                         2*iz - outer_kp_grid(3) - 1 == 0) CYCLE

                     ! use time reversal symmetry k<->-k
                     IF (2*ix - outer_kp_grid(1) - 1 < 0) CYCLE

                     i = i + 1
                     kpt_latt(1) = REAL(2*ix - outer_kp_grid(1) - 1, KIND=dp)/(2._dp*REAL(outer_kp_grid(1), KIND=dp)) &
                                   *reducing_factor(1)
                     kpt_latt(2) = REAL(2*iy - outer_kp_grid(2) - 1, KIND=dp)/(2._dp*REAL(outer_kp_grid(2), KIND=dp)) &
                                   *reducing_factor(2)
                     kpt_latt(3) = REAL(2*iz - outer_kp_grid(3) - 1, KIND=dp)/(2._dp*REAL(outer_kp_grid(3), KIND=dp)) &
                                   *reducing_factor(3)
                     kpoints%xkp(1:3, i) = MATMUL(TRANSPOSE(h_inv), kpt_latt(:))

                     IF (2*ix - outer_kp_grid(1) - 1 == 0) THEN
                        kpoints%wkp(i) = single_weight
                     ELSE
                        kpoints%wkp(i) = 2._dp*single_weight
                     END IF

                  END DO
               END DO
            END DO

            end_kp = i

            kpoint_weight_left = kpoint_weight_left - SUM(kpoints%wkp(start_kp:end_kp))

            reducing_factor(1) = reducing_factor(1)/REAL(outer_kp_grid(1), KIND=dp)
            reducing_factor(2) = reducing_factor(2)/REAL(outer_kp_grid(2), KIND=dp)
            reducing_factor(3) = reducing_factor(3)/REAL(outer_kp_grid(3), KIND=dp)

         END DO

         single_weight = kpoint_weight_left/REAL(nkp_inner_grid, KIND=dp)

         ! the inner grid
         DO ix = 1, kp_grid(1)
            DO iy = 1, kp_grid(2)
               DO iz = 1, kp_grid(3)

                  ! use time reversal symmetry k<->-k
                  IF (2*ix - kp_grid(1) - 1 < 0) CYCLE

                  i = i + 1
                  kpt_latt(1) = REAL(2*ix - kp_grid(1) - 1, KIND=dp)/(2._dp*REAL(kp_grid(1), KIND=dp))*reducing_factor(1)
                  kpt_latt(2) = REAL(2*iy - kp_grid(2) - 1, KIND=dp)/(2._dp*REAL(kp_grid(2), KIND=dp))*reducing_factor(2)
                  kpt_latt(3) = REAL(2*iz - kp_grid(3) - 1, KIND=dp)/(2._dp*REAL(kp_grid(3), KIND=dp))*reducing_factor(3)

                  kpoints%xkp(1:3, i) = MATMUL(TRANSPOSE(h_inv), kpt_latt(:))

                  kpoints%wkp(i) = 2._dp*single_weight

               END DO
            END DO
         END DO

         IF (do_extra_kpoints) THEN

            single_weight = kpoint_weight_left/REAL(kp_grid(1)*kp_grid(2)*kp_grid(3)/8, KIND=dp)

            DO ix = 1, kp_grid(1)/2
               DO iy = 1, kp_grid(2)/2
                  DO iz = 1, kp_grid(3)/2

                     ! use time reversal symmetry k<->-k
                     IF (2*ix - kp_grid(1)/2 - 1 < 0) CYCLE

                     i = i + 1
                     kpt_latt(1) = REAL(2*ix - kp_grid(1)/2 - 1, KIND=dp)/(REAL(kp_grid(1), KIND=dp))
                     kpt_latt(2) = REAL(2*iy - kp_grid(2)/2 - 1, KIND=dp)/(REAL(kp_grid(2), KIND=dp))
                     kpt_latt(3) = REAL(2*iz - kp_grid(3)/2 - 1, KIND=dp)/(REAL(kp_grid(3), KIND=dp))

                     kpoints%xkp(1:3, i) = MATMUL(TRANSPOSE(h_inv), kpt_latt(:))

                     kpoints%wkp(i) = 2._dp*single_weight

                  END DO
               END DO
            END DO

         END IF

         ! default: no symmetry settings
         ALLOCATE (kpoints%kp_sym(kpoints%nkp))
         DO i = 1, kpoints%nkp
            NULLIFY (kpoints%kp_sym(i)%kpoint_sym)
            CALL kpoint_sym_create(kpoints%kp_sym(i)%kpoint_sym)
         END DO

      ELSE

         CALL create_kp_from_gamma(qs_env, qs_env_kp_Gamma_only)

         CALL get_qs_env(qs_env=qs_env, cell=cell, particle_set=particle_set)

         CALL calculate_kp_orbitals(qs_env_kp_Gamma_only, kpoints, "MONKHORST-PACK", nadd=nmo, mp_grid=kp_grid(1:3), &
                                    group_size_ext=para_env%num_pe)

         CALL qs_env_release(qs_env_kp_Gamma_only)

      END IF

   END SUBROUTINE get_kpoints

! **************************************************************************************************
!> \brief ...
!> \param vec_Sigma_c_gw ...
!> \param Eigenval_DFT ...
!> \param eps_eigenval ...
! **************************************************************************************************
   PURE SUBROUTINE average_degenerate_levels(vec_Sigma_c_gw, Eigenval_DFT, eps_eigenval)
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: vec_Sigma_c_gw
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval_DFT
      REAL(KIND=dp), INTENT(IN)                          :: eps_eigenval

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: avg_self_energy
      INTEGER :: degeneracy, first_degenerate_level, i_deg_level, i_level_gw, j_deg_level, jquad, &
         num_deg_levels, num_integ_points, num_levels_gw
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: list_degenerate_levels

      num_levels_gw = SIZE(vec_Sigma_c_gw, 1)

      ALLOCATE (list_degenerate_levels(num_levels_gw))
      list_degenerate_levels = 1

      num_integ_points = SIZE(vec_Sigma_c_gw, 2)

      ALLOCATE (avg_self_energy(num_integ_points))

      DO i_level_gw = 2, num_levels_gw

         IF (ABS(Eigenval_DFT(i_level_gw) - Eigenval_DFT(i_level_gw - 1)) < eps_eigenval) THEN

            list_degenerate_levels(i_level_gw) = list_degenerate_levels(i_level_gw - 1)

         ELSE

            list_degenerate_levels(i_level_gw) = list_degenerate_levels(i_level_gw - 1) + 1

         END IF

      END DO

      num_deg_levels = list_degenerate_levels(num_levels_gw)

      DO i_deg_level = 1, num_deg_levels

         degeneracy = 0

         DO i_level_gw = 1, num_levels_gw

            IF (degeneracy == 0 .AND. i_deg_level == list_degenerate_levels(i_level_gw)) THEN

               first_degenerate_level = i_level_gw

            END IF

            IF (i_deg_level == list_degenerate_levels(i_level_gw)) THEN

               degeneracy = degeneracy + 1

            END IF

         END DO

         DO jquad = 1, num_integ_points

            avg_self_energy(jquad) = SUM(vec_Sigma_c_gw(first_degenerate_level:first_degenerate_level + degeneracy - 1, jquad, 1)) &
                                     /REAL(degeneracy, KIND=dp)

         END DO

         DO j_deg_level = 0, degeneracy - 1

            vec_Sigma_c_gw(first_degenerate_level + j_deg_level, :, 1) = avg_self_energy(:)

         END DO

      END DO

   END SUBROUTINE average_degenerate_levels

! **************************************************************************************************
!> \brief ...
!> \param Eigenval ...
!> \param Eigenval_scf ...
!> \param Eigenval_kp ...
!> \param Eigenval_scf_kp ...
!> \param kpoints ...
!> \param ikp ...
!> \param my_open_shell ...
! **************************************************************************************************
   SUBROUTINE get_eigenval_for_conti(Eigenval, Eigenval_scf, Eigenval_kp, Eigenval_scf_kp, kpoints, ikp, my_open_shell)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval, Eigenval_scf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: Eigenval_kp, Eigenval_scf_kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: ikp
      LOGICAL, INTENT(IN)                                :: my_open_shell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_eigenval_for_conti'

      INTEGER                                            :: handle, ispin, jkp, nkp, nmo, nspin
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      TYPE(mo_set_type), POINTER                         :: mo_set

      CALL timeset(routineN, handle)

      CALL get_kpoint_info(kpoints, nkp=nkp)

      nmo = SIZE(Eigenval)

      IF (my_open_shell) THEN
         nspin = 2
      ELSE
         nspin = 1
      END IF

      IF (ikp == 1) THEN
         ALLOCATE (Eigenval_kp(SIZE(Eigenval), nkp))
         ALLOCATE (Eigenval_scf_kp(SIZE(Eigenval), nkp))

         DO jkp = 1, nkp

            DO ispin = 1, nspin

               mo_set => kpoints%kp_env(jkp)%kpoint_env%mos(1, ispin)%mo_set

               CALL get_mo_set(mo_set=mo_set, eigenvalues=mo_eigenvalues)

               Eigenval_kp(1:nmo, jkp) = mo_eigenvalues(1:nmo)
               Eigenval_scf_kp(1:nmo, jkp) = mo_eigenvalues(1:nmo)

            END DO

         END DO

      END IF

      Eigenval(1:nmo) = Eigenval_kp(1:nmo, ikp)
      Eigenval_scf(1:nmo) = Eigenval_scf_kp(1:nmo, ikp)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param vec_gw_energ ...
!> \param vec_omega_fit_gw ...
!> \param z_value ...
!> \param m_value ...
!> \param vec_Sigma_c_gw ...
!> \param vec_Sigma_x_minus_vxc_gw ...
!> \param Eigenval ...
!> \param Eigenval_scf ...
!> \param n_level_gw ...
!> \param gw_corr_lev_occ ...
!> \param num_poles ...
!> \param num_fit_points ...
!> \param crossing_search ...
!> \param homo ...
!> \param stop_crit ...
!> \param fermi_level_offset ...
!> \param do_gw_im_time ...
! **************************************************************************************************
   SUBROUTINE fit_and_continuation_2pole(vec_gw_energ, vec_omega_fit_gw, &
                                         z_value, m_value, vec_Sigma_c_gw, vec_Sigma_x_minus_vxc_gw, &
                                         Eigenval, Eigenval_scf, n_level_gw, gw_corr_lev_occ, num_poles, &
                                         num_fit_points, crossing_search, homo, stop_crit, &
                                         fermi_level_offset, do_gw_im_time)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: vec_gw_energ, vec_omega_fit_gw, z_value, &
                                                            m_value
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: vec_Sigma_c_gw
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vec_Sigma_x_minus_vxc_gw, Eigenval, &
                                                            Eigenval_scf
      INTEGER, INTENT(IN)                                :: n_level_gw, gw_corr_lev_occ, num_poles, &
                                                            num_fit_points, crossing_search, homo
      REAL(KIND=dp), INTENT(IN)                          :: stop_crit, fermi_level_offset
      LOGICAL, INTENT(IN)                                :: do_gw_im_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fit_and_continuation_2pole'

      COMPLEX(KIND=dp)                                   :: func_val, rho1
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: dLambda, dLambda_2, Lambda, &
                                                            Lambda_without_offset, vec_b_gw, &
                                                            vec_b_gw_copy
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: mat_A_gw, mat_B_gw
      INTEGER                                            :: handle4, ierr, iii, iiter, info, &
                                                            integ_range, jjj, jquad, kkk, &
                                                            max_iter_fit, n_level_gw_ref, num_var, &
                                                            xpos
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipiv
      LOGICAL                                            :: could_exit
      REAL(KIND=dp) :: chi2, chi2_old, delta, deriv_val_real, e_fermi, gw_energ, Ldown, &
         level_energ_GW, Lup, range_step, ScalParam, sign_occ_virt, stat_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Lambda_Im, Lambda_Re, stat_errors, &
                                                            vec_N_gw, vec_omega_fit_gw_sign
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_N_gw

      max_iter_fit = 10000

      num_var = 2*num_poles + 1
      ALLOCATE (Lambda(num_var))
      Lambda = z_zero
      ALLOCATE (Lambda_without_offset(num_var))
      Lambda_without_offset = z_zero
      ALLOCATE (Lambda_Re(num_var))
      Lambda_Re = 0.0_dp
      ALLOCATE (Lambda_Im(num_var))
      Lambda_Im = 0.0_dp

      ALLOCATE (vec_omega_fit_gw_sign(num_fit_points))

      IF (n_level_gw <= gw_corr_lev_occ) THEN
         sign_occ_virt = -1.0_dp
      ELSE
         sign_occ_virt = 1.0_dp
      END IF

      n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ

      DO jquad = 1, num_fit_points
         vec_omega_fit_gw_sign(jquad) = ABS(vec_omega_fit_gw(jquad))*sign_occ_virt
      END DO

      ! initial guess
      range_step = (vec_omega_fit_gw_sign(num_fit_points) - vec_omega_fit_gw_sign(1))/(num_poles - 1)
      DO iii = 1, num_poles
         Lambda_Im(2*iii + 1) = vec_omega_fit_gw_sign(1) + (iii - 1)*range_step
      END DO
      range_step = (vec_omega_fit_gw_sign(num_fit_points) - vec_omega_fit_gw_sign(1))/num_poles
      DO iii = 1, num_poles
         Lambda_Re(2*iii + 1) = ABS(vec_omega_fit_gw_sign(1) + (iii - 0.5_dp)*range_step)
      END DO

      DO iii = 1, num_var
         Lambda(iii) = Lambda_Re(iii) + gaussi*Lambda_Im(iii)
      END DO

      CALL calc_chi2(chi2_old, Lambda, vec_Sigma_c_gw, vec_omega_fit_gw_sign, num_poles, &
                     num_fit_points, n_level_gw)

      ALLOCATE (mat_A_gw(num_poles + 1, num_poles + 1))
      ALLOCATE (vec_b_gw(num_poles + 1))
      ALLOCATE (ipiv(num_poles + 1))
      mat_A_gw = z_zero
      vec_b_gw = 0.0_dp

      mat_A_gw(1:num_poles + 1, 1) = z_one
      integ_range = num_fit_points/num_poles
      DO kkk = 1, num_poles + 1
         xpos = (kkk - 1)*integ_range + 1
         xpos = MIN(xpos, num_fit_points)
         ! calculate coefficient at this point
         DO iii = 1, num_poles
            jjj = iii*2
            func_val = z_one/(gaussi*vec_omega_fit_gw_sign(xpos) - &
                              CMPLX(Lambda_Re(jjj + 1), Lambda_Im(jjj + 1), KIND=dp))
            mat_A_gw(kkk, iii + 1) = func_val
         END DO
         vec_b_gw(kkk) = vec_Sigma_c_gw(n_level_gw, xpos)
      END DO

      ! Solve system of linear equations
      CALL ZGETRF(num_poles + 1, num_poles + 1, mat_A_gw, num_poles + 1, ipiv, info)

      CALL ZGETRS('N', num_poles + 1, 1, mat_A_gw, num_poles + 1, ipiv, vec_b_gw, num_poles + 1, info)

      Lambda_Re(1) = REAL(vec_b_gw(1))
      Lambda_Im(1) = AIMAG(vec_b_gw(1))
      DO iii = 1, num_poles
         jjj = iii*2
         Lambda_Re(jjj) = REAL(vec_b_gw(iii + 1))
         Lambda_Im(jjj) = AIMAG(vec_b_gw(iii + 1))
      END DO

      DEALLOCATE (mat_A_gw)
      DEALLOCATE (vec_b_gw)
      DEALLOCATE (ipiv)

      ALLOCATE (mat_A_gw(num_var*2, num_var*2))
      ALLOCATE (mat_B_gw(num_fit_points, num_var*2))
      ALLOCATE (dLambda(num_fit_points))
      ALLOCATE (dLambda_2(num_fit_points))
      ALLOCATE (vec_b_gw(num_var*2))
      ALLOCATE (vec_b_gw_copy(num_var*2))
      ALLOCATE (ipiv(num_var*2))

      ScalParam = 0.01_dp
      Ldown = 1.5_dp
      Lup = 10.0_dp
      could_exit = .FALSE.

      ! iteration loop for fitting
      DO iiter = 1, max_iter_fit

         CALL timeset(routineN//"_fit_loop_1", handle4)

         ! calc delta lambda
         DO iii = 1, num_var
            Lambda(iii) = Lambda_Re(iii) + gaussi*Lambda_Im(iii)
         END DO
         dLambda = z_zero

         DO kkk = 1, num_fit_points
            func_val = Lambda(1)
            DO iii = 1, num_poles
               jjj = iii*2
               func_val = func_val + Lambda(jjj)/(vec_omega_fit_gw_sign(kkk)*gaussi - Lambda(jjj + 1))
            END DO
            dLambda(kkk) = vec_Sigma_c_gw(n_level_gw, kkk) - func_val
         END DO
         rho1 = SUM(dLambda*dLambda)

         ! fill matrix
         mat_B_gw = z_zero
         DO iii = 1, num_fit_points
            mat_B_gw(iii, 1) = 1.0_dp
            mat_B_gw(iii, num_var + 1) = gaussi
         END DO
         DO iii = 1, num_poles
            jjj = iii*2
            DO kkk = 1, num_fit_points
               mat_B_gw(kkk, jjj) = 1.0_dp/(gaussi*vec_omega_fit_gw_sign(kkk) - Lambda(jjj + 1))
               mat_B_gw(kkk, jjj + num_var) = gaussi/(gaussi*vec_omega_fit_gw_sign(kkk) - Lambda(jjj + 1))
               mat_B_gw(kkk, jjj + 1) = Lambda(jjj)/(gaussi*vec_omega_fit_gw_sign(kkk) - Lambda(jjj + 1))**2
               mat_B_gw(kkk, jjj + 1 + num_var) = (-Lambda_Im(jjj) + gaussi*Lambda_Re(jjj))/ &
                                                  (gaussi*vec_omega_fit_gw_sign(kkk) - Lambda(jjj + 1))**2
            END DO
         END DO

         CALL timestop(handle4)

         CALL timeset(routineN//"_fit_matmul_1", handle4)

         CALL zgemm('C', 'N', num_var*2, num_var*2, num_fit_points, z_one, mat_B_gw, num_fit_points, mat_B_gw, num_fit_points, &
                    z_zero, mat_A_gw, num_var*2)
         CALL timestop(handle4)

         CALL timeset(routineN//"_fit_zgemv_1", handle4)
         CALL zgemv('C', num_fit_points, num_var*2, z_one, mat_B_gw, num_fit_points, dLambda, 1, &
                    z_zero, vec_b_gw, 1)

         CALL timestop(handle4)

         ! scale diagonal elements of a_mat
         DO iii = 1, num_var*2
            mat_A_gw(iii, iii) = mat_A_gw(iii, iii) + ScalParam*mat_A_gw(iii, iii)
         END DO

         ! solve linear system
         ierr = 0
         ipiv = 0

         CALL timeset(routineN//"_fit_lin_eq_2", handle4)

         CALL ZGETRF(2*num_var, 2*num_var, mat_A_gw, 2*num_var, ipiv, info)

         CALL ZGETRS('N', 2*num_var, 1, mat_A_gw, 2*num_var, ipiv, vec_b_gw, 2*num_var, info)

         CALL timestop(handle4)

         DO iii = 1, num_var
            Lambda(iii) = Lambda_Re(iii) + gaussi*Lambda_Im(iii) + vec_b_gw(iii) + vec_b_gw(iii + num_var)
         END DO

         ! calculate chi2
         CALL calc_chi2(chi2, Lambda, vec_Sigma_c_gw, vec_omega_fit_gw_sign, num_poles, &
                        num_fit_points, n_level_gw)

         ! if the fit is already super accurate, exit. otherwise maybe issues when dividing by 0
         IF (chi2 < 1.0E-30_dp) EXIT

         IF (chi2 < chi2_old) THEN
            ScalParam = MAX(ScalParam/Ldown, 1E-12_dp)
            DO iii = 1, num_var
               Lambda_Re(iii) = Lambda_Re(iii) + REAL(vec_b_gw(iii) + vec_b_gw(iii + num_var))
               Lambda_Im(iii) = Lambda_Im(iii) + AIMAG(vec_b_gw(iii) + vec_b_gw(iii + num_var))
            END DO
            IF (chi2_old/chi2 - 1.0_dp < stop_crit) could_exit = .TRUE.
            chi2_old = chi2
         ELSE
            ScalParam = ScalParam*Lup
         END IF
         IF (ScalParam > 100.0_dp .AND. could_exit) EXIT

         IF (ScalParam > 1E+10_dp) ScalParam = 1E-4_dp

      END DO

      IF (.NOT. do_gw_im_time) THEN

         ! change a_0 [Lambda(1)], so that Sigma(i0) = Fit(i0)
         ! do not do this for imaginary time since we do not have many fit points and the fit should be perfect
         func_val = Lambda(1)
         DO iii = 1, num_poles
            jjj = iii*2
            ! calculate value of the fit function
            func_val = func_val + Lambda(jjj)/(-Lambda(jjj + 1))
         END DO

         Lambda_Re(1) = Lambda_Re(1) - REAL(func_val) + REAL(vec_Sigma_c_gw(n_level_gw, num_fit_points))
         Lambda_Im(1) = Lambda_Im(1) - AIMAG(func_val) + AIMAG(vec_Sigma_c_gw(n_level_gw, num_fit_points))

      END IF

      Lambda_without_offset(:) = Lambda(:)

      DO iii = 1, num_var
         Lambda(iii) = CMPLX(Lambda_Re(iii), Lambda_Im(iii), KIND=dp)
      END DO

      IF (do_gw_im_time) THEN
         ! for cubic-scaling GW, we have one Green's function for occ and virt states with the Fermi level
         ! in the middle of homo and lumo
         e_fermi = 0.5_dp*(Eigenval(homo) + Eigenval(homo + 1))
      ELSE
         ! in case of O(N^4) GW, we have the Fermi level differently for occ and virt states, see
         ! Fig. 1 in JCTC 12, 3623-3635 (2016)
         IF (n_level_gw <= gw_corr_lev_occ) THEN
            e_fermi = Eigenval(homo) + fermi_level_offset
         ELSE
            e_fermi = Eigenval(homo + 1) - fermi_level_offset
         END IF
      END IF

      ! either Z-shot or Newton/bisection crossing search for evaluating Sigma_c
      IF (crossing_search == ri_rpa_g0w0_crossing_z_shot .OR. &
          crossing_search == ri_rpa_g0w0_crossing_newton) THEN

         ! calculate Sigma_c_fit(e_n) and Z
         func_val = Lambda(1)
         z_value(n_level_gw) = 1.0_dp
         DO iii = 1, num_poles
            jjj = iii*2
            z_value(n_level_gw) = z_value(n_level_gw) + REAL(Lambda(jjj)/ &
                                                             (Eigenval(n_level_gw_ref) - e_fermi - Lambda(jjj + 1))**2)
            func_val = func_val + Lambda(jjj)/(Eigenval(n_level_gw_ref) - e_fermi - Lambda(jjj + 1))
         END DO
         ! m is the slope of the correl self-energy
         m_value(n_level_gw) = 1.0_dp - z_value(n_level_gw)
         z_value(n_level_gw) = 1.0_dp/z_value(n_level_gw)
         gw_energ = REAL(func_val)
         vec_gw_energ(n_level_gw) = gw_energ

         ! in case one wants to do Newton-Raphson on top of the Z-shot
         IF (crossing_search == ri_rpa_g0w0_crossing_newton) THEN

            level_energ_GW = (Eigenval_scf(n_level_gw_ref) - &
                              m_value(n_level_gw)*Eigenval(n_level_gw_ref) + &
                              vec_gw_energ(n_level_gw) + &
                              vec_Sigma_x_minus_vxc_gw(n_level_gw_ref))* &
                             z_value(n_level_gw)

            ! Newton-Raphson iteration
            DO kkk = 1, 1000

               ! calculate the value of the fit function for level_energ_GW
               func_val = Lambda(1)
               z_value(n_level_gw) = 1.0_dp
               DO iii = 1, num_poles
                  jjj = iii*2
                  func_val = func_val + Lambda(jjj)/(level_energ_GW - e_fermi - Lambda(jjj + 1))
               END DO

               ! calculate the derivative of the fit function for level_energ_GW
               deriv_val_real = -1.0_dp
               DO iii = 1, num_poles
                  jjj = iii*2
                  deriv_val_real = deriv_val_real + REAL(Lambda(jjj))/((ABS(level_energ_GW - e_fermi - Lambda(jjj + 1)))**2) &
                                   - (REAL(Lambda(jjj))*(level_energ_GW - e_fermi) - REAL(Lambda(jjj)*CONJG(Lambda(jjj + 1))))* &
                                   2.0_dp*(level_energ_GW - e_fermi - REAL(Lambda(jjj + 1)))/ &
                                   ((ABS(level_energ_GW - e_fermi - Lambda(jjj + 1)))**2)

               END DO

              delta = (Eigenval_scf(n_level_gw_ref) + vec_Sigma_x_minus_vxc_gw(n_level_gw_ref) + REAL(func_val) - level_energ_GW)/ &
                       deriv_val_real

               level_energ_GW = level_energ_GW - delta

               IF (ABS(delta) < 1.0E-08) EXIT

            END DO

            ! update the GW-energy by Newton-Raphson and set the Z-value to 1

            vec_gw_energ(n_level_gw) = REAL(func_val)
            z_value(n_level_gw) = 1.0_dp
            m_value(n_level_gw) = 0.0_dp

         END IF ! Newton-Raphson on top of Z-shot

      ELSE
         CPABORT("Only NONE, ZSHOT and NEWTON implemented for 2-pole model")
      END IF ! decision crossing search none, Z-shot

      !   --------------------------------------------
      !  | calculate statistical error due to fitting |
      !   --------------------------------------------

      ! estimate the statistical error of the calculated Sigma_c(i*omega)
      ! by sqrt(chi2/n), where n is the number of fit points

      CALL calc_chi2(chi2, Lambda_without_offset, vec_Sigma_c_gw, vec_omega_fit_gw_sign, num_poles, &
                     num_fit_points, n_level_gw)

      ! Estimate the statistical error of every fit point
      stat_error = SQRT(chi2/num_fit_points)

      ! allocate N array containing the second derivatives of chi^2
      ALLOCATE (vec_N_gw(num_var*2))
      vec_N_gw = 0.0_dp

      ALLOCATE (mat_N_gw(num_var*2, num_var*2))
      mat_N_gw = 0.0_dp

      DO iii = 1, num_var*2
         CALL calc_mat_N(vec_N_gw(iii), Lambda_without_offset, vec_Sigma_c_gw, vec_omega_fit_gw_sign, &
                         iii, iii, num_poles, num_fit_points, n_level_gw, 0.001_dp)
      END DO

      DO iii = 1, num_var*2
         DO jjj = 1, num_var*2
            CALL calc_mat_N(mat_N_gw(iii, jjj), Lambda_without_offset, vec_Sigma_c_gw, vec_omega_fit_gw_sign, &
                            iii, jjj, num_poles, num_fit_points, n_level_gw, 0.001_dp)
         END DO
      END DO

      CALL DGETRF(2*num_var, 2*num_var, mat_N_gw, 2*num_var, ipiv, info)

      ! vec_b_gw is only working array
      CALL DGETRI(2*num_var, mat_N_gw, 2*num_var, ipiv, vec_b_gw, 2*num_var, info)

      ALLOCATE (stat_errors(2*num_var))
      stat_errors = 0.0_dp

      DO iii = 1, 2*num_var
         stat_errors(iii) = SQRT(ABS(mat_N_gw(iii, iii)))*stat_error
      END DO

      DEALLOCATE (mat_N_gw)
      DEALLOCATE (vec_N_gw)
      DEALLOCATE (mat_A_gw)
      DEALLOCATE (mat_B_gw)
      DEALLOCATE (stat_errors)
      DEALLOCATE (dLambda)
      DEALLOCATE (dLambda_2)
      DEALLOCATE (vec_b_gw)
      DEALLOCATE (vec_b_gw_copy)
      DEALLOCATE (ipiv)
      DEALLOCATE (vec_omega_fit_gw_sign)
      DEALLOCATE (Lambda)
      DEALLOCATE (Lambda_without_offset)
      DEALLOCATE (Lambda_Re)
      DEALLOCATE (Lambda_Im)

   END SUBROUTINE fit_and_continuation_2pole

! **************************************************************************************************
!> \brief perform analytic continuation with pade approximation
!> \param vec_gw_energ real Sigma_c
!> \param vec_omega_fit_gw frequency points for Sigma_c(iomega)
!> \param z_value 1/(1-dev)
!> \param m_value derivative of real Sigma_c
!> \param vec_Sigma_c_gw complex Sigma_c(iomega)
!> \param vec_Sigma_x_minus_vxc_gw ...
!> \param Eigenval quasiparticle energy during ev self-consistent GW
!> \param Eigenval_scf KS/HF eigenvalue
!> \param n_level_gw ...
!> \param gw_corr_lev_occ ...
!> \param nparam_pade number of pade parameters
!> \param num_fit_points number of fit points for Sigma_c(iomega)
!> \param crossing_search type ofr cross search to find quasiparticle energies
!> \param homo ...
!> \param fermi_level_offset ...
!> \param do_gw_im_time ...
!> \param print_self_energy ...
!> \param count_ev_sc_GW ...
! **************************************************************************************************
   SUBROUTINE continuation_pade(vec_gw_energ, vec_omega_fit_gw, &
                                z_value, m_value, vec_Sigma_c_gw, vec_Sigma_x_minus_vxc_gw, &
                                Eigenval, Eigenval_scf, n_level_gw, gw_corr_lev_occ, nparam_pade, &
                                num_fit_points, crossing_search, homo, &
                                fermi_level_offset, do_gw_im_time, print_self_energy, count_ev_sc_GW)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: vec_gw_energ
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vec_omega_fit_gw
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: z_value, m_value
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: vec_Sigma_c_gw
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vec_Sigma_x_minus_vxc_gw, Eigenval, &
                                                            Eigenval_scf
      INTEGER, INTENT(IN)                                :: n_level_gw, gw_corr_lev_occ, &
                                                            nparam_pade, num_fit_points, &
                                                            crossing_search, homo
      REAL(KIND=dp), INTENT(IN)                          :: fermi_level_offset
      LOGICAL, INTENT(IN)                                :: do_gw_im_time, print_self_energy
      INTEGER, INTENT(IN)                                :: count_ev_sc_GW

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'continuation_pade'

      CHARACTER(len=default_path_length)                 :: filename
      COMPLEX(KIND=dp)                                   :: sigma_c_pade, sigma_c_pade_im_freq
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: coeff_pade, omega_points_pade, &
                                                            Sigma_c_gw_reorder
      INTEGER                                            :: handle, i_omega, iunit, jquad, &
                                                            n_level_gw_ref, num_omega
      REAL(KIND=dp)                                      :: e_fermi, energy_val, level_energ_GW, &
                                                            omega, sign_occ_virt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vec_omega_fit_gw_sign, &
                                                            vec_omega_fit_gw_sign_reorder

      CALL timeset(routineN, handle)

      ALLOCATE (vec_omega_fit_gw_sign(num_fit_points))

      IF (n_level_gw <= gw_corr_lev_occ) THEN
         sign_occ_virt = -1.0_dp
      ELSE
         sign_occ_virt = 1.0_dp
      END IF

      DO jquad = 1, num_fit_points
         vec_omega_fit_gw_sign(jquad) = ABS(vec_omega_fit_gw(jquad))*sign_occ_virt
      END DO

      IF (do_gw_im_time) THEN
         ! for cubic-scaling GW, we have one Green's function for occ and virt states with the Fermi level
         ! in the middle of homo and lumo
         e_fermi = 0.5_dp*(Eigenval(homo) + Eigenval(homo + 1))
      ELSE
         ! in case of O(N^4) GW, we have the Fermi level differently for occ and virt states, see
         ! Fig. 1 in JCTC 12, 3623-3635 (2016)
         IF (n_level_gw <= gw_corr_lev_occ) THEN
            e_fermi = Eigenval(homo) + fermi_level_offset
         ELSE
            e_fermi = Eigenval(homo + 1) - fermi_level_offset
         END IF
      END IF

      n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ

      !*** reorder, such that omega=i*0 is first entry
      ALLOCATE (Sigma_c_gw_reorder(num_fit_points))
      ALLOCATE (vec_omega_fit_gw_sign_reorder(num_fit_points))
      ! for cubic scaling GW fit points are ordered differently than in N^4 GW
      IF (do_gw_im_time) THEN
         DO jquad = 1, num_fit_points
            Sigma_c_gw_reorder(jquad) = vec_Sigma_c_gw(n_level_gw, jquad)
            vec_omega_fit_gw_sign_reorder(jquad) = vec_omega_fit_gw_sign(jquad)
         ENDDO
      ELSE
         DO jquad = 1, num_fit_points
            Sigma_c_gw_reorder(jquad) = vec_Sigma_c_gw(n_level_gw, num_fit_points - jquad + 1)
            vec_omega_fit_gw_sign_reorder(jquad) = vec_omega_fit_gw_sign(num_fit_points - jquad + 1)
         ENDDO
      ENDIF

      !*** evaluate parameters for pade approximation
      ALLOCATE (coeff_pade(nparam_pade))
      ALLOCATE (omega_points_pade(nparam_pade))
      coeff_pade = 0.0_dp
      CALL get_pade_parameters(Sigma_c_gw_reorder, vec_omega_fit_gw_sign_reorder, &
                               num_fit_points, nparam_pade, omega_points_pade, coeff_pade)

      !*** calculate start_value for iterative cross-searching methods
      IF ((crossing_search == ri_rpa_g0w0_crossing_bisection) .OR. &
          (crossing_search == ri_rpa_g0w0_crossing_newton)) THEN
         energy_val = Eigenval(n_level_gw_ref) - e_fermi
         CALL evaluate_pade_function(energy_val, nparam_pade, omega_points_pade, &
                                     coeff_pade, sigma_c_pade)
         CALL get_z_and_m_value_pade(energy_val, nparam_pade, omega_points_pade, &
                                     coeff_pade, z_value(n_level_gw), m_value(n_level_gw))
         level_energ_GW = (Eigenval_scf(n_level_gw_ref) - &
                           m_value(n_level_gw)*Eigenval(n_level_gw_ref) + &
                           REAL(sigma_c_pade) + &
                           vec_Sigma_x_minus_vxc_gw(n_level_gw_ref))* &
                          z_value(n_level_gw)
      ENDIF

      !*** perform crossing search
      SELECT CASE (crossing_search)
      CASE (ri_rpa_g0w0_crossing_z_shot)
         energy_val = Eigenval(n_level_gw_ref) - e_fermi
         CALL evaluate_pade_function(energy_val, nparam_pade, omega_points_pade, &
                                     coeff_pade, sigma_c_pade)
         vec_gw_energ(n_level_gw) = REAL(sigma_c_pade)

         CALL get_z_and_m_value_pade(energy_val, nparam_pade, omega_points_pade, &
                                     coeff_pade, z_value(n_level_gw), m_value(n_level_gw))

      CASE (ri_rpa_g0w0_crossing_bisection)
         CALL get_sigma_c_bisection_pade(vec_gw_energ(n_level_gw), Eigenval_scf(n_level_gw_ref), &
                                         vec_Sigma_x_minus_vxc_gw(n_level_gw_ref), e_fermi, &
                                         nparam_pade, omega_points_pade, coeff_pade, &
                                         n_level_gw_ref, start_val=level_energ_GW)
         z_value(n_level_gw) = 1.0_dp
         m_value(n_level_gw) = 0.0_dp

      CASE (ri_rpa_g0w0_crossing_newton)
         CALL get_sigma_c_newton_pade(vec_gw_energ(n_level_gw), Eigenval_scf(n_level_gw_ref), &
                                      vec_Sigma_x_minus_vxc_gw(n_level_gw_ref), e_fermi, &
                                      nparam_pade, omega_points_pade, coeff_pade, &
                                      n_level_gw_ref, start_val=level_energ_GW)
         z_value(n_level_gw) = 1.0_dp
         m_value(n_level_gw) = 0.0_dp

      CASE DEFAULT
         CPABORT("Only Z_SHOT, NEWTON, and BISECTION crossing search implemented.")
      END SELECT

      IF (print_self_energy) THEN

         IF (count_ev_sc_GW == 1) THEN

            IF (n_level_gw_ref < 10) THEN
               WRITE (filename, "(A26,I1)") "G0W0_self_energy_level_000", n_level_gw_ref
            ELSE IF (n_level_gw_ref < 100) THEN
               WRITE (filename, "(A25,I2)") "G0W0_self_energy_level_00", n_level_gw_ref
            ELSE IF (n_level_gw_ref < 1000) THEN
               WRITE (filename, "(A24,I3)") "G0W0_self_energy_level_0", n_level_gw_ref
            ELSE
               WRITE (filename, "(A23,I4)") "G0W0_self_energy_level_", n_level_gw_ref
            END IF

         ELSE

            IF (n_level_gw_ref < 10) THEN
               WRITE (filename, "(A11,I1,A22,I1)") "evGW_cycle_", count_ev_sc_GW, &
                  "_self_energy_level_000", n_level_gw_ref
            ELSE IF (n_level_gw_ref < 100) THEN
               WRITE (filename, "(A11,I1,A21,I2)") "evGW_cycle_", count_ev_sc_GW, &
                  "_self_energy_level_00", n_level_gw_ref
            ELSE IF (n_level_gw_ref < 1000) THEN
               WRITE (filename, "(A11,I1,A20,I3)") "evGW_cycle_", count_ev_sc_GW, &
                  "_self_energy_level_0", n_level_gw_ref
            ELSE
               WRITE (filename, "(A11,I1,A19,I4)") "evGW_cycle_", count_ev_sc_GW, &
                  "_self_energy_level_", n_level_gw_ref
            END IF

         END IF

         CALL open_file(TRIM(filename), unit_number=iunit, file_status="UNKNOWN", file_action="WRITE")

         num_omega = 10000

         WRITE (iunit, "(2A42)") " omega (eV)     Sigma(omega) (eV)  ", &
            "  omega - e_n^DFT - Sigma_n^x - v_n^xc (eV)"

         DO i_omega = 0, num_omega

            omega = -50.0_dp/evolt + REAL(i_omega, KIND=dp)/REAL(num_omega, KIND=dp)*100.0_dp/evolt

            CALL evaluate_pade_function(omega - e_fermi, nparam_pade, omega_points_pade, &
                                        coeff_pade, sigma_c_pade)

            WRITE (iunit, "(F12.2,2F17.5)") omega*evolt, REAL(sigma_c_pade)*evolt, &
               (omega - Eigenval_scf(n_level_gw_ref) - vec_Sigma_x_minus_vxc_gw(n_level_gw_ref))*evolt

         END DO

         WRITE (iunit, "(A51,A39)") " w (eV)  Re(Sigma(i*w)) (eV)   Im(Sigma(i*w)) (eV) ", &
            "  Re(Fit(i*w)) (eV)    Im(Fit(iw)) (eV)"

         DO jquad = 1, num_fit_points

            CALL evaluate_pade_function(vec_omega_fit_gw_sign_reorder(jquad), &
                                        nparam_pade, omega_points_pade, &
                                        coeff_pade, sigma_c_pade_im_freq, do_imag_freq=.TRUE.)

            WRITE (iunit, "(F12.2,4F17.5)") vec_omega_fit_gw_sign_reorder(jquad)*evolt, &
               REAL(Sigma_c_gw_reorder(jquad)*evolt), &
               AIMAG(Sigma_c_gw_reorder(jquad)*evolt), &
               REAL(sigma_c_pade_im_freq*evolt), &
               AIMAG(sigma_c_pade_im_freq*evolt)

         END DO

         CALL close_file(iunit)

      END IF

      DEALLOCATE (vec_omega_fit_gw_sign)
      DEALLOCATE (Sigma_c_gw_reorder)
      DEALLOCATE (vec_omega_fit_gw_sign_reorder)
      DEALLOCATE (coeff_pade, omega_points_pade)

      CALL timestop(handle)

   END SUBROUTINE continuation_pade

! **************************************************************************************************
!> \brief calculate pade parameter recursively as in  Eq. (A2) in J. Low Temp. Phys., Vol. 29,
!>          1977, pp. 179
!> \param y f(x), here: Sigma_c(iomega)
!> \param x the frequency points omega
!> \param num_fit_points ...
!> \param nparam number of pade parameters
!> \param xpoints set of points used in pade approximation, selection of x
!> \param coeff pade coefficients
! **************************************************************************************************
   PURE SUBROUTINE get_pade_parameters(y, x, num_fit_points, nparam, xpoints, coeff)

      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: y
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: x
      INTEGER, INTENT(IN)                                :: num_fit_points, nparam
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(INOUT)      :: xpoints, coeff

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_pade_parameters'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: ypoints
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: g_mat
      INTEGER                                            :: idat, iparam, nstep

      nstep = INT(num_fit_points/(nparam - 1))

      ALLOCATE (ypoints(nparam))
      !omega=i0 is in element x(1)
      idat = 1
      DO iparam = 1, nparam - 1
         xpoints(iparam) = gaussi*x(idat)
         ypoints(iparam) = y(idat)
         idat = idat + nstep
      ENDDO
      xpoints(nparam) = gaussi*x(num_fit_points)
      ypoints(nparam) = y(num_fit_points)

      !*** generate parameters recursively

      ALLOCATE (g_mat(nparam, nparam))
      g_mat(:, 1) = ypoints(:)
      DO iparam = 2, nparam
         DO idat = iparam, nparam
            g_mat(idat, iparam) = (g_mat(iparam - 1, iparam - 1) - g_mat(idat, iparam - 1))/ &
                                  ((xpoints(idat) - xpoints(iparam - 1))*g_mat(idat, iparam - 1))
         ENDDO
      ENDDO

      DO iparam = 1, nparam
         coeff(iparam) = g_mat(iparam, iparam)
      ENDDO

      DEALLOCATE (ypoints)
      DEALLOCATE (g_mat)

   END SUBROUTINE get_pade_parameters

! **************************************************************************************************
!> \brief evaluate pade function for a real value x_val
!> \param x_val real value
!> \param nparam number of pade parameters
!> \param xpoints selection of points of the original complex function, i.e. here of Sigma_c(iomega)
!> \param coeff pade coefficients
!> \param func_val function value
!> \param do_imag_freq ...
! **************************************************************************************************
   PURE SUBROUTINE evaluate_pade_function(x_val, nparam, xpoints, coeff, func_val, do_imag_freq)

      REAL(KIND=dp), INTENT(IN)                          :: x_val
      INTEGER, INTENT(IN)                                :: nparam
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: xpoints, coeff
      COMPLEX(KIND=dp), INTENT(OUT)                      :: func_val
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_imag_freq

      CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_pade_function'

      INTEGER                                            :: iparam
      LOGICAL                                            :: my_do_imag_freq

      my_do_imag_freq = .FALSE.
      IF (PRESENT(do_imag_freq)) my_do_imag_freq = do_imag_freq

      func_val = z_one
      DO iparam = nparam, 2, -1
         IF (my_do_imag_freq) THEN
            func_val = z_one + coeff(iparam)*(gaussi*x_val - xpoints(iparam - 1))/func_val
         ELSE
            func_val = z_one + coeff(iparam)*(x_val*z_one - xpoints(iparam - 1))/func_val
         END IF
      ENDDO

      func_val = coeff(1)/func_val

   END SUBROUTINE evaluate_pade_function

! **************************************************************************************************
!> \brief get the z-value and the m-value (derivative) of the pade function
!> \param x_val real value
!> \param nparam number of pade parameters
!> \param xpoints selection of points of the original complex function, i.e. here of Sigma_c(iomega)
!> \param coeff pade coefficients
!> \param z_value 1/(1-dev)
!> \param m_value derivative
! **************************************************************************************************
   PURE SUBROUTINE get_z_and_m_value_pade(x_val, nparam, xpoints, coeff, z_value, m_value)

      REAL(KIND=dp), INTENT(IN)                          :: x_val
      INTEGER, INTENT(IN)                                :: nparam
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: xpoints, coeff
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: z_value, m_value

      COMPLEX(KIND=dp)                                   :: denominator, dev_denominator, &
                                                            dev_numerator, dev_val, func_val, &
                                                            numerator
      INTEGER                                            :: iparam

      func_val = z_one
      dev_val = z_zero
      DO iparam = nparam, 2, -1
         numerator = coeff(iparam)*(x_val*z_one - xpoints(iparam - 1))
         dev_numerator = coeff(iparam)*z_one
         denominator = func_val
         dev_denominator = dev_val
         dev_val = dev_numerator/denominator - (numerator*dev_denominator)/(denominator**2)
         func_val = z_one + coeff(iparam)*(x_val*z_one - xpoints(iparam - 1))/func_val
      ENDDO

      dev_val = -1.0_dp*coeff(1)/(func_val**2)*dev_val
      func_val = coeff(1)/func_val

      IF (PRESENT(z_value)) THEN
         z_value = 1.0_dp - REAL(dev_val)
         z_value = 1.0_dp/z_value
      ENDIF
      IF (PRESENT(m_value)) m_value = REAL(dev_val)

   END SUBROUTINE get_z_and_m_value_pade

! **************************************************************************************************
!> \brief crossing search using the bisection method to find the quasiparticle energy
!> \param gw_energ real Sigma_c
!> \param Eigenval_scf Eigenvalue from the SCF
!> \param Sigma_x_minus_vxc_gw ...
!> \param e_fermi fermi level
!> \param nparam_pade number of pade parameters
!> \param omega_points_pade selection of frequency points of Sigma_c(iomega)
!> \param coeff_pade pade coefficients
!> \param n_level_gw_ref ...
!> \param start_val start value for the quasiparticle iteration
! **************************************************************************************************
   SUBROUTINE get_sigma_c_bisection_pade(gw_energ, Eigenval_scf, Sigma_x_minus_vxc_gw, e_fermi, &
                                         nparam_pade, omega_points_pade, coeff_pade, n_level_gw_ref, start_val)

      REAL(KIND=dp), INTENT(OUT)                         :: gw_energ
      REAL(KIND=dp), INTENT(IN)                          :: Eigenval_scf, Sigma_x_minus_vxc_gw, &
                                                            e_fermi
      INTEGER, INTENT(IN)                                :: nparam_pade
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: omega_points_pade, coeff_pade
      INTEGER, INTENT(IN)                                :: n_level_gw_ref
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: start_val

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_sigma_c_bisection_pade'

      CHARACTER(LEN=512)                                 :: error_msg
      CHARACTER(LEN=64)                                  :: n_level_gw_ref_char
      COMPLEX(KIND=dp)                                   :: sigma_c
      INTEGER                                            :: handle, icount
      REAL(KIND=dp)                                      :: delta, energy_val, my_start_val, &
                                                            qp_energy, qp_energy_old, threshold

      CALL timeset(routineN, handle)

      threshold = 1.0E-7_dp

      IF (PRESENT(start_val)) THEN
         my_start_val = start_val
      ELSE
         my_start_val = Eigenval_scf
      ENDIF

      qp_energy = my_start_val
      qp_energy_old = my_start_val
      delta = 1.0E-3_dp

      icount = 0
      DO WHILE (ABS(delta) > threshold)
         icount = icount + 1
         qp_energy = qp_energy_old + 0.5_dp*delta
         qp_energy_old = qp_energy
         energy_val = qp_energy - e_fermi
         CALL evaluate_pade_function(energy_val, nparam_pade, omega_points_pade, &
                                     coeff_pade, sigma_c)
         qp_energy = Eigenval_scf + REAL(sigma_c) + Sigma_x_minus_vxc_gw
         delta = qp_energy - qp_energy_old
         IF (icount > 500) THEN
            WRITE (n_level_gw_ref_char, '(I10)') n_level_gw_ref
            WRITE (error_msg, '(A,A,A)') " Self-consistent quasi-particle solution of "// &
               "MO ", TRIM(n_level_gw_ref_char), " has not been found."
            CPWARN(error_msg)
            EXIT
         ENDIF
      ENDDO

      gw_energ = REAL(sigma_c)

      CALL timestop(handle)

   END SUBROUTINE get_sigma_c_bisection_pade

! **************************************************************************************************
!> \brief crossing search using the Newton method to find the quasiparticle energy
!> \param gw_energ real Sigma_c
!> \param Eigenval_scf Eigenvalue from the SCF
!> \param Sigma_x_minus_vxc_gw ...
!> \param e_fermi fermi level
!> \param nparam_pade number of pade parameters
!> \param omega_points_pade selection of frequency points of Sigma_c(iomega)
!> \param coeff_pade pade coefficients
!> \param n_level_gw_ref ...
!> \param start_val start value for the quasiparticle iteration
! **************************************************************************************************
   SUBROUTINE get_sigma_c_newton_pade(gw_energ, Eigenval_scf, Sigma_x_minus_vxc_gw, e_fermi, &
                                      nparam_pade, omega_points_pade, coeff_pade, n_level_gw_ref, start_val)

      REAL(KIND=dp), INTENT(OUT)                         :: gw_energ
      REAL(KIND=dp), INTENT(IN)                          :: Eigenval_scf, Sigma_x_minus_vxc_gw, &
                                                            e_fermi
      INTEGER, INTENT(IN)                                :: nparam_pade
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: omega_points_pade, coeff_pade
      INTEGER, INTENT(IN)                                :: n_level_gw_ref
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: start_val

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_sigma_c_newton_pade'

      CHARACTER(LEN=512)                                 :: error_msg
      CHARACTER(LEN=64)                                  :: n_level_gw_ref_char
      COMPLEX(KIND=dp)                                   :: sigma_c
      INTEGER                                            :: handle, icount
      REAL(KIND=dp)                                      :: delta, energy_val, m_value, &
                                                            my_start_val, qp_energy, &
                                                            qp_energy_old, threshold

      CALL timeset(routineN, handle)

      threshold = 1.0E-7_dp

      IF (PRESENT(start_val)) THEN
         my_start_val = start_val
      ELSE
         my_start_val = Eigenval_scf
      ENDIF

      qp_energy = my_start_val
      qp_energy_old = my_start_val
      delta = 1.0E-3_dp

      icount = 0
      DO WHILE (ABS(delta) > threshold)
         icount = icount + 1
         energy_val = qp_energy - e_fermi
         CALL evaluate_pade_function(energy_val, nparam_pade, omega_points_pade, &
                                     coeff_pade, sigma_c)
         !get m_value --> derivative of function
         CALL get_z_and_m_value_pade(energy_val, nparam_pade, omega_points_pade, &
                                     coeff_pade, m_value=m_value)
         qp_energy_old = qp_energy
         qp_energy = qp_energy - (Eigenval_scf + Sigma_x_minus_vxc_gw + REAL(sigma_c) - qp_energy)/ &
                     (m_value - 1.0_dp)
         delta = qp_energy - qp_energy_old
         IF (icount > 500) THEN
            WRITE (n_level_gw_ref_char, '(I10)') n_level_gw_ref
            WRITE (error_msg, '(A,A,A)') " Self-consistent quasi-particle solution of "// &
               "MO ", TRIM(n_level_gw_ref_char), " has not been found."
            CPWARN(error_msg)
            EXIT
         ENDIF
      ENDDO

      gw_energ = REAL(sigma_c)

      CALL timestop(handle)

   END SUBROUTINE get_sigma_c_newton_pade

! **************************************************************************************************
!> \brief Prints the GW stuff to the output and optinally to an external file.
!>        Also updates the eigenvalues for eigenvalue-self-consistent GW
!> \param vec_gw_energ ...
!> \param z_value ...
!> \param m_value ...
!> \param vec_Sigma_x_minus_vxc_gw ...
!> \param Eigenval ...
!> \param Eigenval_last ...
!> \param Eigenval_scf ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_tot ...
!> \param crossing_search ...
!> \param homo ...
!> \param unit_nr ...
!> \param count_ev_sc_GW ...
!> \param count_sc_GW0 ...
!> \param ikp ...
!> \param nkp_self_energy ...
!> \param kpoints ...
!> \param ispin requested spin-state (1 for alpha, 2 for beta, else closed-shell)
! **************************************************************************************************
   SUBROUTINE print_and_update_for_ev_sc(vec_gw_energ, &
                                         z_value, m_value, vec_Sigma_x_minus_vxc_gw, Eigenval, &
                                         Eigenval_last, Eigenval_scf, gw_corr_lev_occ, gw_corr_lev_tot, &
                                         crossing_search, homo, unit_nr, count_ev_sc_GW, count_sc_GW0, &
                                         ikp, nkp_self_energy, kpoints, ispin)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vec_gw_energ, z_value, m_value
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: vec_Sigma_x_minus_vxc_gw, Eigenval, &
                                                            Eigenval_last, Eigenval_scf
      INTEGER, INTENT(IN) :: gw_corr_lev_occ, gw_corr_lev_tot, crossing_search, homo, unit_nr, &
         count_ev_sc_GW, count_sc_GW0, ikp, nkp_self_energy
      TYPE(kpoint_type), INTENT(IN), POINTER             :: kpoints
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_and_update_for_ev_sc'

      CHARACTER(4)                                       :: occ_virt
      INTEGER                                            :: handle, n_level_gw, n_level_gw_ref
      LOGICAL                                            :: do_alpha, do_beta, do_closed_shell, &
                                                            do_kpoints, is_energy_okay
      REAL(KIND=dp)                                      :: new_energy

      CALL timeset(routineN, handle)

      do_alpha = (ispin == 1)
      do_beta = (ispin == 2)
      do_closed_shell = .NOT. (do_alpha .OR. do_beta)
      do_kpoints = (nkp_self_energy > 1)

      Eigenval_last(:) = Eigenval(:)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, *) ' '

         IF (count_ev_sc_GW == 1 .AND. count_sc_GW0 == 1) THEN

            IF (do_closed_shell) THEN
               WRITE (unit_nr, *) ' '
               WRITE (unit_nr, '(T3,A)') '******************************************************************************'
               WRITE (unit_nr, '(T3,A)') '**                                                                          **'
               WRITE (unit_nr, '(T3,A)') '**                        GW QUASIPARTICLE ENERGIES                         **'
               WRITE (unit_nr, '(T3,A)') '**                                                                          **'
               WRITE (unit_nr, '(T3,A)') '******************************************************************************'

            ELSE IF (do_alpha) THEN
               WRITE (unit_nr, *) ' '
               WRITE (unit_nr, '(T3,A)') '---------------------------------------'
               WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies of alpha spins'
               WRITE (unit_nr, '(T3,A)') '----------------------------------------'
            ELSE IF (do_beta) THEN
               WRITE (unit_nr, *) ' '
               WRITE (unit_nr, '(T3,A)') '---------------------------------------'
               WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies of beta spins'
               WRITE (unit_nr, '(T3,A)') '---------------------------------------'
            END IF

            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A)') 'The GW quasiparticle energies are calculated according to: '
            IF (crossing_search == ri_rpa_g0w0_crossing_z_shot) THEN
               WRITE (unit_nr, '(T3,A)') 'E_GW = E_SCF + Z * ( Sigc(E_SCF) + Sigx - vxc )'
            ELSE
               WRITE (unit_nr, '(T3,A)') ' '
               WRITE (unit_nr, '(T3,A)') '                    E_GW = E_SCF + Sigc(E_GW) + Sigx - vxc '
               WRITE (unit_nr, '(T3,A)') ' '
               WRITE (unit_nr, '(T3,A)') 'Upper equation is solved self-consistently for E_GW, see Eq. (12) in J. Phys.'
               WRITE (unit_nr, '(T3,A)') 'Chem. Lett. 9, 306 (2018), doi: 10.1021/acs.jpclett.7b02740'
            END IF

            WRITE (unit_nr, *) ' '
            WRITE (unit_nr, *) ' '
            WRITE (unit_nr, '(T3,A)') '------------'
            WRITE (unit_nr, '(T3,A)') 'G0W0 results'
            WRITE (unit_nr, '(T3,A)') '------------'

         END IF

         IF (count_ev_sc_GW > 1) THEN
            WRITE (unit_nr, *) ' '
            WRITE (unit_nr, '(T3,A)') '---------------------------------------'
            WRITE (unit_nr, '(T3,A,I4)') 'Eigenvalue-selfconsistency cycle: ', count_ev_sc_GW
            WRITE (unit_nr, '(T3,A)') '---------------------------------------'
         END IF

         IF (count_sc_GW0 > 1) THEN
            WRITE (unit_nr, '(T3,A)') '----------------------------------'
            WRITE (unit_nr, '(T3,A,I4)') 'scGW0 selfconsistency cycle: ', count_sc_GW0
            WRITE (unit_nr, '(T3,A)') '----------------------------------'
         END IF

         IF (do_kpoints) THEN
            WRITE (unit_nr, *) ' '
            WRITE (unit_nr, '(T3,A7,I3,A3,I3,A8,3F7.3,A12,3F7.3)') 'Kpoint ', ikp, '  /', nkp_self_energy, &
               '   xkp =', kpoints%xkp(1, ikp), kpoints%xkp(2, ikp), kpoints%xkp(3, ikp), &
               '  and  xkp =', -kpoints%xkp(1, ikp), -kpoints%xkp(2, ikp), -kpoints%xkp(3, ikp)
            WRITE (unit_nr, '(T3,A72)') '(Relative Brillouin zone size: [-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5])'
         END IF

      END IF

      DO n_level_gw = 1, gw_corr_lev_tot

         n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ

         new_energy = (Eigenval_scf(n_level_gw_ref) - &
                       m_value(n_level_gw)*Eigenval(n_level_gw_ref) + &
                       vec_gw_energ(n_level_gw) + &
                       vec_Sigma_x_minus_vxc_gw(n_level_gw_ref))* &
                      z_value(n_level_gw)

         is_energy_okay = .TRUE.

         IF (n_level_gw_ref > homo .AND. new_energy < Eigenval(homo)) THEN
            is_energy_okay = .FALSE.
         END IF

         IF (is_energy_okay) THEN
            Eigenval(n_level_gw_ref) = new_energy
         END IF

      END DO

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A)') ' '
         IF (crossing_search == ri_rpa_g0w0_crossing_z_shot) THEN
            WRITE (unit_nr, '(T13,2A)') 'MO    E_SCF (eV)    Sigc (eV)   Sigx-vxc (eV)    Z         E_GW (eV)'
         ELSE
            WRITE (unit_nr, '(T3,2A)') 'Molecular orbital   E_SCF (eV)       Sigc (eV)   Sigx-vxc (eV)       E_GW (eV)'
         END IF
      END IF

      DO n_level_gw = 1, gw_corr_lev_tot
         n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
         IF (n_level_gw <= gw_corr_lev_occ) THEN
            occ_virt = 'occ'
         ELSE
            occ_virt = 'vir'
         END IF

         IF (unit_nr > 0) THEN
            IF (crossing_search == ri_rpa_g0w0_crossing_z_shot) THEN
               WRITE (unit_nr, '(T3,I4,3A,5F13.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ') ', &
                  Eigenval_last(n_level_gw_ref)*evolt, &
                  vec_gw_energ(n_level_gw)*evolt, &
                  vec_Sigma_x_minus_vxc_gw(n_level_gw_ref)*evolt, &
                  z_value(n_level_gw), &
                  Eigenval(n_level_gw_ref)*evolt
            ELSE
               WRITE (unit_nr, '(T3,I4,3A,4F16.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ')  ', &
                  Eigenval_last(n_level_gw_ref)*evolt, &
                  vec_gw_energ(n_level_gw)*evolt, &
                  vec_Sigma_x_minus_vxc_gw(n_level_gw_ref)*evolt, &
                  Eigenval(n_level_gw_ref)*evolt
            END IF
         END IF
      END DO

      IF (unit_nr > 0) THEN
         IF (do_closed_shell) THEN
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A,F57.2)') 'GW HOMO-LUMO gap (eV)', (Eigenval(homo + 1) - Eigenval(homo))*evolt
         ELSE IF (do_alpha) THEN
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A,F51.2)') 'Alpha GW HOMO-LUMO gap (eV)', (Eigenval(homo + 1) - Eigenval(homo))*evolt
         ELSE IF (do_beta) THEN
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A,F52.2)') 'Beta GW HOMO-LUMO gap (eV)', (Eigenval(homo + 1) - Eigenval(homo))*evolt
         END IF
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_and_update_for_ev_sc

! **************************************************************************************************
!> \brief ...
!> \param Eigenval ...
!> \param Eigenval_last ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param nmo ...
! **************************************************************************************************
   PURE SUBROUTINE shift_unshifted_levels(Eigenval, Eigenval_last, gw_corr_lev_occ, gw_corr_lev_virt, &
                                          homo, nmo)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval, Eigenval_last
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt, homo, &
                                                            nmo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'shift_unshifted_levels'

      INTEGER                                            :: n_level_gw, n_level_gw_ref
      REAL(KIND=dp)                                      :: eigen_diff

      ! for eigenvalue self-consistent GW, all eigenvalues have to be corrected
      ! 1) the occupied; check if there are occupied MOs not being corrected by GW
      IF (gw_corr_lev_occ < homo .AND. gw_corr_lev_occ > 0) THEN

         ! calculate average GW correction for occupied orbitals
         eigen_diff = 0.0_dp

         DO n_level_gw = 1, gw_corr_lev_occ
            n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
            eigen_diff = eigen_diff + Eigenval(n_level_gw_ref) - Eigenval_last(n_level_gw_ref)
         END DO
         eigen_diff = eigen_diff/gw_corr_lev_occ

         ! correct the eigenvalues of the occupied orbitals which have not been corrected by GW
         DO n_level_gw = 1, homo - gw_corr_lev_occ
            Eigenval(n_level_gw) = Eigenval(n_level_gw) + eigen_diff
         END DO

      END IF

      ! 2) the virtual: check if there are virtual orbitals not being corrected by GW
      IF (gw_corr_lev_virt < nmo - homo .AND. gw_corr_lev_virt > 0) THEN

         ! calculate average GW correction for virtual orbitals
         eigen_diff = 0.0_dp
         DO n_level_gw = 1, gw_corr_lev_virt
            n_level_gw_ref = n_level_gw + homo
            eigen_diff = eigen_diff + Eigenval(n_level_gw_ref) - Eigenval_last(n_level_gw_ref)
         END DO
         eigen_diff = eigen_diff/gw_corr_lev_virt

         ! correct the eigenvalues of the virtual orbitals which have not been corrected by GW
         DO n_level_gw = homo + gw_corr_lev_virt + 1, nmo
            Eigenval(n_level_gw) = Eigenval(n_level_gw) + eigen_diff
         END DO

      END IF

   END SUBROUTINE shift_unshifted_levels

! **************************************************************************************************
!> \brief Calculate the matrix mat_N_gw containing the second derivatives
!>        with respect to the fitting parameters. The second derivatives are
!>        calculated numerically by finite differences.
!> \param N_ij matrix element
!> \param Lambda fitting parameters
!> \param Sigma_c ...
!> \param vec_omega_fit_gw ...
!> \param i ...
!> \param j ...
!> \param num_poles ...
!> \param num_fit_points ...
!> \param n_level_gw ...
!> \param h  ...
! **************************************************************************************************
   SUBROUTINE calc_mat_N(N_ij, Lambda, Sigma_c, vec_omega_fit_gw, i, j, &
                         num_poles, num_fit_points, n_level_gw, h)
      REAL(KIND=dp), INTENT(OUT)                         :: N_ij
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: Lambda
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: Sigma_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: vec_omega_fit_gw
      INTEGER, INTENT(IN)                                :: i, j, num_poles, num_fit_points, &
                                                            n_level_gw
      REAL(KIND=dp), INTENT(IN)                          :: h

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calc_mat_N'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: Lambda_tmp
      INTEGER                                            :: handle, num_var
      REAL(KIND=dp)                                      :: chi2, chi2_sum

      CALL timeset(routineN, handle)

      num_var = 2*num_poles + 1
      ALLOCATE (Lambda_tmp(num_var))
      Lambda_tmp = z_zero
      chi2_sum = 0.0_dp

      !test
      Lambda_tmp(:) = Lambda(:)
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)

      ! Fitting parameters with offset h
      Lambda_tmp(:) = Lambda(:)
      IF (MODULO(i, 2) == 0) THEN
         Lambda_tmp(i/2) = Lambda_tmp(i/2) + h*z_one
      ELSE
         Lambda_tmp((i + 1)/2) = Lambda_tmp((i + 1)/2) + h*gaussi
      END IF
      IF (MODULO(j, 2) == 0) THEN
         Lambda_tmp(j/2) = Lambda_tmp(j/2) + h*z_one
      ELSE
         Lambda_tmp((j + 1)/2) = Lambda_tmp((j + 1)/2) + h*gaussi
      END IF
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)
      chi2_sum = chi2_sum + chi2

      IF (MODULO(i, 2) == 0) THEN
         Lambda_tmp(i/2) = Lambda_tmp(i/2) - 2.0_dp*h*z_one
      ELSE
         Lambda_tmp((i + 1)/2) = Lambda_tmp((i + 1)/2) - 2.0_dp*h*gaussi
      END IF
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)
      chi2_sum = chi2_sum - chi2

      IF (MODULO(j, 2) == 0) THEN
         Lambda_tmp(j/2) = Lambda_tmp(j/2) - 2.0_dp*h*z_one
      ELSE
         Lambda_tmp((j + 1)/2) = Lambda_tmp((j + 1)/2) - 2.0_dp*h*gaussi
      END IF
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)
      chi2_sum = chi2_sum + chi2

      IF (MODULO(i, 2) == 0) THEN
         Lambda_tmp(i/2) = Lambda_tmp(i/2) + 2.0_dp*h*z_one
      ELSE
         Lambda_tmp((i + 1)/2) = Lambda_tmp((i + 1)/2) + 2.0_dp*h*gaussi
      END IF
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)
      chi2_sum = chi2_sum - chi2

      ! Second derivative with symmetric difference quotient
      N_ij = 1.0_dp/2.0_dp*chi2_sum/(4.0_dp*h*h)

      DEALLOCATE (Lambda_tmp)

      CALL timestop(handle)

   END SUBROUTINE calc_mat_N

! **************************************************************************************************
!> \brief Calculate chi2
!> \param chi2 ...
!> \param Lambda fitting parameters
!> \param Sigma_c ...
!> \param vec_omega_fit_gw ...
!> \param num_poles ...
!> \param num_fit_points ...
!> \param n_level_gw ...
! **************************************************************************************************
   PURE SUBROUTINE calc_chi2(chi2, Lambda, Sigma_c, vec_omega_fit_gw, num_poles, &
                             num_fit_points, n_level_gw)
      REAL(KIND=dp), INTENT(OUT)                         :: chi2
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: Lambda
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: Sigma_c
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vec_omega_fit_gw
      INTEGER, INTENT(IN)                                :: num_poles, num_fit_points, n_level_gw

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calc_chi2'

      COMPLEX(KIND=dp)                                   :: func_val
      INTEGER                                            :: iii, jjj, kkk

      chi2 = 0.0_dp
      DO kkk = 1, num_fit_points
         func_val = Lambda(1)
         DO iii = 1, num_poles
            jjj = iii*2
            ! calculate value of the fit function
            func_val = func_val + Lambda(jjj)/(gaussi*vec_omega_fit_gw(kkk) - Lambda(jjj + 1))
         END DO
         chi2 = chi2 + (ABS(Sigma_c(n_level_gw, kkk) - func_val))**2
      END DO

   END SUBROUTINE calc_chi2

! **************************************************************************************************
!> \brief ...
!> \param num_integ_points ...
!> \param nmo ...
!> \param tau_tj ...
!> \param tj ...
!> \param matrix_s ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_mo_coeff_occ_scaled ...
!> \param fm_mo_coeff_virt_scaled ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param Eigenval ...
!> \param eps_filter ...
!> \param e_fermi ...
!> \param fm_mat_W ...
!> \param gw_corr_lev_tot ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param count_ev_sc_GW ...
!> \param count_sc_GW0 ...
!> \param t_3c_overl_int_ao_mo ...
!> \param t_3c_O_mo_compressed ...
!> \param t_3c_O_mo_ind ...
!> \param t_3c_overl_int_gw_RI ...
!> \param t_3c_overl_int_gw_AO ...
!> \param mat_W ...
!> \param mat_SinvVSinv ...
!> \param mat_dm ...
!> \param weights_cos_tf_t_to_w ...
!> \param weights_sin_tf_t_to_w ...
!> \param vec_Sigma_c_gw ...
!> \param do_periodic ...
!> \param num_points_corr ...
!> \param delta_corr ...
!> \param qs_env ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param mp2_env ...
!> \param matrix_berry_re_mo_mo ...
!> \param matrix_berry_im_mo_mo ...
!> \param first_cycle_periodic_correction ...
!> \param kpoints ...
!> \param num_fit_points ...
!> \param fm_mo_coeff ...
!> \param do_ri_Sigma_x ...
!> \param vec_Sigma_x_gw ...
!> \param unit_nr ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE compute_self_energy_cubic_gw(num_integ_points, nmo, tau_tj, tj, &
                                           matrix_s, fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
                                           fm_mo_coeff_virt_scaled, fm_scaled_dm_occ_tau, &
                                           fm_scaled_dm_virt_tau, Eigenval, eps_filter, &
                                           e_fermi, fm_mat_W, &
                                           gw_corr_lev_tot, gw_corr_lev_occ, gw_corr_lev_virt, homo, &
                                           count_ev_sc_GW, count_sc_GW0, &
                                           t_3c_overl_int_ao_mo, t_3c_O_mo_compressed, t_3c_O_mo_ind, &
                                           t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                           mat_W, mat_SinvVSinv, mat_dm, &
                                           weights_cos_tf_t_to_w, weights_sin_tf_t_to_w, vec_Sigma_c_gw, &
                                           do_periodic, num_points_corr, delta_corr, qs_env, para_env, para_env_RPA, &
                                           mp2_env, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                           first_cycle_periodic_correction, kpoints, num_fit_points, fm_mo_coeff, &
                                           do_ri_Sigma_x, vec_Sigma_x_gw, unit_nr, ispin)
      INTEGER, INTENT(IN)                                :: num_integ_points, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tj
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: matrix_s
      TYPE(cp_fm_type), POINTER :: fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
         fm_mo_coeff_virt_scaled, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      REAL(KIND=dp), INTENT(INOUT)                       :: e_fermi
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_W
      INTEGER, INTENT(IN)                                :: gw_corr_lev_tot, gw_corr_lev_occ, &
                                                            gw_corr_lev_virt, homo, &
                                                            count_ev_sc_GW, count_sc_GW0
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_ao_mo
      TYPE(hfx_compression_type)                         :: t_3c_O_mo_compressed
      INTEGER, DIMENSION(:, :)                           :: t_3c_O_mo_ind
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_gw_RI, &
                                                            t_3c_overl_int_gw_AO
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: mat_W
      TYPE(dbcsr_p_type)                                 :: mat_SinvVSinv, mat_dm
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: weights_cos_tf_t_to_w, &
                                                            weights_sin_tf_t_to_w
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: vec_Sigma_c_gw
      LOGICAL, INTENT(IN)                                :: do_periodic
      INTEGER, INTENT(IN)                                :: num_points_corr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: delta_corr
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      LOGICAL, INTENT(INOUT) :: first_cycle_periodic_correction
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: num_fit_points
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff
      LOGICAL, INTENT(IN)                                :: do_ri_Sigma_x
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: vec_Sigma_x_gw
      INTEGER, INTENT(IN)                                :: unit_nr, ispin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_self_energy_cubic_gw'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: delta_corr_omega
      INTEGER :: gw_lev_end, gw_lev_start, handle, handle3, i, iblk_mo, iquad, jquad, mo_end, &
         mo_start, n_level_gw, n_level_gw_ref, nblk_mo, unit_nr_prv
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: batch_range_mo, dist1, dist2, mo_bsizes, &
                                                            mo_offsets, sizes_AO, sizes_RI
      INTEGER, DIMENSION(2)                              :: mo_bounds, pdims_2d
      LOGICAL                                            :: memory_info
      REAL(KIND=dp)                                      :: ext_scaling, omega, omega_i, omega_sign, &
                                                            sign_occ_virt, t_i_Clenshaw, tau, &
                                                            weight_cos, weight_i, weight_sin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: vec_Sigma_c_gw_cos_omega, &
         vec_Sigma_c_gw_cos_tau, vec_Sigma_c_gw_neg_tau, vec_Sigma_c_gw_pos_tau, &
         vec_Sigma_c_gw_sin_omega, vec_Sigma_c_gw_sin_tau
      TYPE(dbcsr_t_pgrid_type)                           :: pgrid_2d
      TYPE(dbcsr_t_type)                                 :: t_3c_ctr_AO, t_3c_ctr_RI, t_AO_tmp, &
                                                            t_dm, t_greens_fct_occ, &
                                                            t_greens_fct_virt, t_RI_tmp, &
                                                            t_SinvVSinv, t_W
      TYPE(dbcsr_type), TARGET                           :: mat_greens_fct_occ, mat_greens_fct_virt

      CALL timeset(routineN, handle)

      CALL decompress_tensor(t_3c_overl_int_ao_mo, t_3c_O_mo_ind, t_3c_O_mo_compressed, &
                             mp2_env%ri_rpa_im_time%eps_compress)

      CALL dbcsr_t_copy(t_3c_overl_int_ao_mo, t_3c_overl_int_gw_RI)
      CALL dbcsr_t_copy(t_3c_overl_int_ao_mo, t_3c_overl_int_gw_AO, order=[2, 1, 3], move_data=.TRUE.)

      memory_info = mp2_env%ri_rpa_im_time%memory_info
      IF (memory_info) THEN
         unit_nr_prv = unit_nr
      ELSE
         unit_nr_prv = 0
      ENDIF

      mo_start = homo - gw_corr_lev_occ + 1
      mo_end = homo + gw_corr_lev_virt
      CPASSERT(mo_end - mo_start + 1 == gw_corr_lev_tot)

      vec_Sigma_c_gw = z_zero
      ALLOCATE (vec_Sigma_c_gw_pos_tau(gw_corr_lev_tot, num_integ_points))
      vec_Sigma_c_gw_pos_tau = 0.0_dp
      ALLOCATE (vec_Sigma_c_gw_neg_tau(gw_corr_lev_tot, num_integ_points))
      vec_Sigma_c_gw_neg_tau = 0.0_dp
      ALLOCATE (vec_Sigma_c_gw_cos_tau(gw_corr_lev_tot, num_integ_points))
      vec_Sigma_c_gw_cos_tau = 0.0_dp
      ALLOCATE (vec_Sigma_c_gw_sin_tau(gw_corr_lev_tot, num_integ_points))
      vec_Sigma_c_gw_sin_tau = 0.0_dp

      ALLOCATE (vec_Sigma_c_gw_cos_omega(gw_corr_lev_tot, num_integ_points))
      vec_Sigma_c_gw_cos_omega = 0.0_dp
      ALLOCATE (vec_Sigma_c_gw_sin_omega(gw_corr_lev_tot, num_integ_points))
      vec_Sigma_c_gw_sin_omega = 0.0_dp

      ALLOCATE (delta_corr_omega(1 + homo - gw_corr_lev_occ:homo + gw_corr_lev_virt, num_integ_points))
      delta_corr_omega(:, :) = z_zero

      CALL dbcsr_create(matrix=mat_greens_fct_occ, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_create(matrix=mat_greens_fct_virt, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      e_fermi = 0.5_dp*(Eigenval(homo) + Eigenval(homo + 1))

      nblk_mo = dbcsr_t_nblks_total(t_3c_overl_int_gw_AO, 3)
      ALLOCATE (mo_offsets(nblk_mo))
      ALLOCATE (mo_bsizes(nblk_mo))
      ALLOCATE (batch_range_mo(nblk_mo - 1))
      CALL dbcsr_t_get_info(t_3c_overl_int_gw_AO, blk_offset_3=mo_offsets, blk_size_3=mo_bsizes)

      pdims_2d = 0
      CALL dbcsr_t_pgrid_create(para_env%group, pdims_2d, pgrid_2d)
      ALLOCATE (sizes_RI(dbcsr_t_nblks_total(t_3c_overl_int_gw_RI, 1)))
      CALL dbcsr_t_get_info(t_3c_overl_int_gw_RI, blk_size_1=sizes_RI)

      CALL create_2c_tensor(t_W, dist1, dist2, pgrid_2d, sizes_RI, sizes_RI, name="(RI|RI)")

      DEALLOCATE (dist1, dist2)

      CALL dbcsr_t_create(mat_W, t_RI_tmp, name="(RI|RI)")

      CALL dbcsr_t_create(t_3c_overl_int_gw_RI, t_3c_ctr_RI)
      CALL dbcsr_t_create(t_3c_overl_int_gw_AO, t_3c_ctr_AO)

      ALLOCATE (sizes_AO(dbcsr_t_nblks_total(t_3c_overl_int_gw_AO, 1)))
      CALL dbcsr_t_get_info(t_3c_overl_int_gw_AO, blk_size_1=sizes_AO)
      CALL create_2c_tensor(t_greens_fct_occ, dist1, dist2, pgrid_2d, sizes_AO, sizes_AO, name="(AO|AO)")
      DEALLOCATE (dist1, dist2)
      CALL create_2c_tensor(t_greens_fct_virt, dist1, dist2, pgrid_2d, sizes_AO, sizes_AO, name="(AO|AO)")
      DEALLOCATE (dist1, dist2)

      DO jquad = 1, num_integ_points

         CALL compute_Greens_function_time(mat_greens_fct_occ, mat_greens_fct_virt, &
                                           fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                           fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                           fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, Eigenval, &
                                           nmo, eps_filter, e_fermi, tau_tj(jquad), para_env)

         CALL dbcsr_set(mat_W, 0.0_dp)
         CALL copy_fm_to_dbcsr(fm_mat_W(jquad)%matrix, mat_W, keep_sparsity=.FALSE.)

         IF (jquad == 1) CALL dbcsr_t_create(mat_greens_fct_occ, t_AO_tmp, name="(AO|AO)")

         CALL dbcsr_t_copy_matrix_to_tensor(mat_W, t_RI_tmp)
         CALL dbcsr_t_copy(t_RI_tmp, t_W)
         CALL dbcsr_t_copy_matrix_to_tensor(mat_greens_fct_occ, t_AO_tmp)
         CALL dbcsr_t_copy(t_AO_tmp, t_greens_fct_occ)
         CALL dbcsr_t_copy_matrix_to_tensor(mat_greens_fct_virt, t_AO_tmp)
         CALL dbcsr_t_copy(t_AO_tmp, t_greens_fct_virt)

         batch_range_mo(:) = [(i, i=2, nblk_mo)]
         CALL dbcsr_t_batched_contract_init(t_3c_overl_int_gw_AO, batch_range_3=batch_range_mo)
         CALL dbcsr_t_batched_contract_init(t_3c_overl_int_gw_RI, batch_range_3=batch_range_mo)
         CALL dbcsr_t_batched_contract_init(t_3c_ctr_AO, batch_range_3=batch_range_mo)
         CALL dbcsr_t_batched_contract_init(t_3c_ctr_RI, batch_range_3=batch_range_mo)
         CALL dbcsr_t_batched_contract_init(t_W)
         CALL dbcsr_t_batched_contract_init(t_greens_fct_occ)
         CALL dbcsr_t_batched_contract_init(t_greens_fct_virt)

         ! in iteration over MO blocks skip first and last block because they correspond to the MO s
         ! outside of the GW range of required MOs
         DO iblk_mo = 2, nblk_mo - 1
            mo_bounds = [mo_offsets(iblk_mo), mo_offsets(iblk_mo) + mo_bsizes(iblk_mo) - 1]
            CALL contract_cubic_gw(t_3c_overl_int_gw_AO, t_3c_overl_int_gw_RI, &
                                   t_greens_fct_occ, t_W, [1.0_dp, -1.0_dp], &
                                   mo_bounds, unit_nr_prv, &
                                   t_3c_ctr_RI, t_3c_ctr_AO, calculate_ctr_ri=.TRUE.)
            CALL trace_sigma_gw(t_3c_ctr_AO, t_3c_ctr_RI, vec_Sigma_c_gw_neg_tau(:, jquad), mo_start, mo_bounds, para_env)

            CALL contract_cubic_gw(t_3c_overl_int_gw_AO, t_3c_overl_int_gw_RI, &
                                   t_greens_fct_virt, t_W, [1.0_dp, 1.0_dp], &
                                   mo_bounds, unit_nr_prv, &
                                   t_3c_ctr_RI, t_3c_ctr_AO, calculate_ctr_ri=.FALSE.)

            CALL trace_sigma_gw(t_3c_ctr_AO, t_3c_ctr_RI, vec_Sigma_c_gw_pos_tau(:, jquad), mo_start, mo_bounds, para_env)
         ENDDO
         CALL dbcsr_t_batched_contract_finalize(t_3c_overl_int_gw_AO)
         CALL dbcsr_t_batched_contract_finalize(t_3c_overl_int_gw_RI)
         CALL dbcsr_t_batched_contract_finalize(t_3c_ctr_AO)
         CALL dbcsr_t_batched_contract_finalize(t_3c_ctr_RI)
         CALL dbcsr_t_batched_contract_finalize(t_W)
         CALL dbcsr_t_batched_contract_finalize(t_greens_fct_occ)
         CALL dbcsr_t_batched_contract_finalize(t_greens_fct_virt)

         CALL dbcsr_t_clear(t_3c_ctr_AO)
         CALL dbcsr_t_clear(t_3c_ctr_RI)

         vec_Sigma_c_gw_cos_tau(:, jquad) = 0.5_dp*(vec_Sigma_c_gw_pos_tau(:, jquad) + &
                                                    vec_Sigma_c_gw_neg_tau(:, jquad))

         vec_Sigma_c_gw_sin_tau(:, jquad) = 0.5_dp*(vec_Sigma_c_gw_pos_tau(:, jquad) - &
                                                    vec_Sigma_c_gw_neg_tau(:, jquad))

      END DO ! jquad (tau)
      CALL dbcsr_t_destroy(t_W)

      CALL dbcsr_t_destroy(t_greens_fct_occ)
      CALL dbcsr_t_destroy(t_greens_fct_virt)

      ! Fourier transform from time to frequency
      DO jquad = 1, num_fit_points

         DO iquad = 1, num_integ_points

            omega = tj(jquad)
            tau = tau_tj(iquad)
            weight_cos = weights_cos_tf_t_to_w(jquad, iquad)*COS(omega*tau)
            weight_sin = weights_sin_tf_t_to_w(jquad, iquad)*SIN(omega*tau)

            vec_Sigma_c_gw_cos_omega(:, jquad) = vec_Sigma_c_gw_cos_omega(:, jquad) + &
                                                 weight_cos*vec_Sigma_c_gw_cos_tau(:, iquad)

            vec_Sigma_c_gw_sin_omega(:, jquad) = vec_Sigma_c_gw_sin_omega(:, jquad) + &
                                                 weight_sin*vec_Sigma_c_gw_sin_tau(:, iquad)

         END DO

      END DO

      ! for occupied levels, we need the correlation self-energy for negative omega. Therefore, weight_sin
      ! should be computed with -omega, which results in an additional minus for vec_Sigma_c_gw_sin_omega:
      vec_Sigma_c_gw_sin_omega(1:gw_corr_lev_occ, :) = -vec_Sigma_c_gw_sin_omega(1:gw_corr_lev_occ, :)

      vec_Sigma_c_gw(:, 1:num_fit_points, 1) = vec_Sigma_c_gw_cos_omega(:, 1:num_fit_points) + &
                                               gaussi*vec_Sigma_c_gw_sin_omega(:, 1:num_fit_points)

      CALL dbcsr_release(mat_greens_fct_occ)
      CALL dbcsr_release(mat_greens_fct_virt)

      IF (do_ri_Sigma_x .AND. count_ev_sc_GW == 1 .AND. count_sc_GW0 == 1) THEN

         CALL timeset(routineN//"_RI_HFX_operation_1", handle3)

         ! get density matrix
         CALL cp_gemm(transa="N", transb="T", m=nmo, n=nmo, k=nmo, alpha=1.0_dp, &
                      matrix_a=fm_mo_coeff_occ, matrix_b=fm_mo_coeff_occ, beta=0.0_dp, &
                      matrix_c=fm_scaled_dm_occ_tau)

         CALL timestop(handle3)

         CALL timeset(routineN//"_RI_HFX_operation_2", handle3)

         CALL copy_fm_to_dbcsr(fm_scaled_dm_occ_tau, &
                               mat_dm%matrix, &
                               keep_sparsity=.FALSE.)

         CALL timestop(handle3)

         CALL create_2c_tensor(t_dm, dist1, dist2, pgrid_2d, sizes_AO, sizes_AO, name="(AO|AO)")
         DEALLOCATE (dist1, dist2)

         CALL dbcsr_t_copy_matrix_to_tensor(mat_dm%matrix, t_AO_tmp)
         CALL dbcsr_t_copy(t_AO_tmp, t_dm)

         CALL create_2c_tensor(t_SinvVSinv, dist1, dist2, pgrid_2d, sizes_RI, sizes_RI, name="(RI|RI)")
         DEALLOCATE (dist1, dist2)

         CALL dbcsr_t_copy_matrix_to_tensor(mat_SinvVSinv%matrix, t_RI_tmp)
         CALL dbcsr_t_copy(t_RI_tmp, t_SinvVSinv)

         CALL dbcsr_t_batched_contract_init(t_3c_overl_int_gw_AO, batch_range_3=batch_range_mo)
         CALL dbcsr_t_batched_contract_init(t_3c_overl_int_gw_RI, batch_range_3=batch_range_mo)
         CALL dbcsr_t_batched_contract_init(t_3c_ctr_RI, batch_range_3=batch_range_mo)
         CALL dbcsr_t_batched_contract_init(t_3c_ctr_AO, batch_range_3=batch_range_mo)
         CALL dbcsr_t_batched_contract_init(t_dm)
         CALL dbcsr_t_batched_contract_init(t_SinvVSinv)

         DO iblk_mo = 2, nblk_mo - 1
            mo_bounds = [mo_offsets(iblk_mo), mo_offsets(iblk_mo) + mo_bsizes(iblk_mo) - 1]

            CALL contract_cubic_gw(t_3c_overl_int_gw_AO, t_3c_overl_int_gw_RI, &
                                   t_dm, t_SinvVSinv, [1.0_dp, -1.0_dp], &
                                   mo_bounds, unit_nr_prv, &
                                   t_3c_ctr_RI, t_3c_ctr_AO, calculate_ctr_ri=.TRUE.)

            CALL trace_sigma_gw(t_3c_ctr_AO, t_3c_ctr_RI, vec_Sigma_x_gw(mo_start:mo_end, 1), mo_start, mo_bounds, para_env)
         ENDDO
         CALL dbcsr_t_batched_contract_finalize(t_3c_overl_int_gw_AO)
         CALL dbcsr_t_batched_contract_finalize(t_3c_overl_int_gw_RI)
         CALL dbcsr_t_batched_contract_finalize(t_dm)
         CALL dbcsr_t_batched_contract_finalize(t_SinvVSinv)
         CALL dbcsr_t_batched_contract_finalize(t_3c_ctr_RI)
         CALL dbcsr_t_batched_contract_finalize(t_3c_ctr_AO)

         CALL dbcsr_t_destroy(t_dm)
         CALL dbcsr_t_destroy(t_SinvVSinv)

         mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, ispin, 1) = &
            mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, ispin, 1) + &
            vec_Sigma_x_gw(:, 1)

      END IF

      CALL dbcsr_t_pgrid_destroy(pgrid_2d)

      CALL dbcsr_t_destroy(t_3c_ctr_RI)
      CALL dbcsr_t_destroy(t_3c_ctr_AO)
      CALL dbcsr_t_destroy(t_AO_tmp)
      CALL dbcsr_t_destroy(t_RI_tmp)

      ! compute and add the periodic correction
      IF (do_periodic) THEN

         ext_scaling = 0.2_dp

         ! loop over omega' (integration)
         DO iquad = 1, num_points_corr

            ! use the Clenshaw-grid
            t_i_Clenshaw = iquad*pi/(2.0_dp*num_points_corr)
            omega_i = ext_scaling/TAN(t_i_Clenshaw)

            IF (iquad < num_points_corr) THEN
               weight_i = ext_scaling*pi/(num_points_corr*SIN(t_i_Clenshaw)**2)
            ELSE
               weight_i = ext_scaling*pi/(2.0_dp*num_points_corr*SIN(t_i_Clenshaw)**2)
            END IF

            CALL calc_periodic_correction(delta_corr, qs_env, para_env, para_env_RPA, &
                                          mp2_env%ri_g0w0%kp_grid, homo, nmo, gw_corr_lev_occ, &
                                          gw_corr_lev_virt, omega_i, fm_mo_coeff, Eigenval, &
                                          matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                          first_cycle_periodic_correction, kpoints, &
                                          mp2_env%ri_g0w0%do_mo_coeff_gamma, &
                                          mp2_env%ri_g0w0%num_kp_grids, mp2_env%ri_g0w0%eps_kpoint, &
                                          mp2_env%ri_g0w0%do_extra_kpoints, &
                                          mp2_env%ri_g0w0%do_aux_bas_gw, mp2_env%ri_g0w0%frac_aux_mos)

            DO n_level_gw = 1, gw_corr_lev_tot

               n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ

               IF (n_level_gw <= gw_corr_lev_occ) THEN
                  sign_occ_virt = -1.0_dp
               ELSE
                  sign_occ_virt = 1.0_dp
               END IF

               DO jquad = 1, num_integ_points

                  omega_sign = tj(jquad)*sign_occ_virt

                  delta_corr_omega(n_level_gw_ref, jquad) = &
                     delta_corr_omega(n_level_gw_ref, jquad) - &
                     0.5_dp/pi*weight_i/2.0_dp*delta_corr(n_level_gw_ref)* &
                     (1.0_dp/(gaussi*(omega_i + omega_sign) + e_fermi - Eigenval(n_level_gw_ref)) + &
                      1.0_dp/(gaussi*(-omega_i + omega_sign) + e_fermi - Eigenval(n_level_gw_ref)))

               END DO

            END DO

         END DO

         gw_lev_start = 1 + homo - gw_corr_lev_occ
         gw_lev_end = homo + gw_corr_lev_virt

         ! add the periodic correction
         vec_Sigma_c_gw(1:gw_corr_lev_tot, :, 1) = vec_Sigma_c_gw(1:gw_corr_lev_tot, :, 1) + &
                                                   delta_corr_omega(gw_lev_start:gw_lev_end, 1:num_fit_points)

      END IF

      DEALLOCATE (vec_Sigma_c_gw_pos_tau)
      DEALLOCATE (vec_Sigma_c_gw_neg_tau)
      DEALLOCATE (vec_Sigma_c_gw_cos_tau)
      DEALLOCATE (vec_Sigma_c_gw_sin_tau)
      DEALLOCATE (vec_Sigma_c_gw_cos_omega)
      DEALLOCATE (vec_Sigma_c_gw_sin_omega)
      DEALLOCATE (delta_corr_omega)

      CALL timestop(handle)

   END SUBROUTINE compute_self_energy_cubic_gw

! **************************************************************************************************
!> \brief ...
!> \param t_3c_overl_int_gw_AO ...
!> \param t_3c_overl_int_gw_RI ...
!> \param t_AO ...
!> \param t_RI ...
!> \param prefac ...
!> \param mo_bounds ...
!> \param unit_nr ...
!> \param t_3c_ctr_RI ...
!> \param t_3c_ctr_AO ...
!> \param calculate_ctr_RI ...
! **************************************************************************************************
   SUBROUTINE contract_cubic_gw(t_3c_overl_int_gw_AO, t_3c_overl_int_gw_RI, &
                                t_AO, t_RI, prefac, &
                                mo_bounds, unit_nr, &
                                t_3c_ctr_RI, t_3c_ctr_AO, calculate_ctr_RI)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t_3c_overl_int_gw_AO, &
                                                            t_3c_overl_int_gw_RI, t_AO, t_RI
      REAL(dp), DIMENSION(2), INTENT(IN)                 :: prefac
      INTEGER, DIMENSION(2), INTENT(IN)                  :: mo_bounds
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t_3c_ctr_RI, t_3c_ctr_AO
      LOGICAL, INTENT(IN)                                :: calculate_ctr_RI

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_cubic_gw'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(2, 2)                           :: ctr_bounds_mo
      INTEGER, DIMENSION(3)                              :: bounds_3c

      CALL timeset(routineN, handle)

      IF (calculate_ctr_RI) THEN
         CALL dbcsr_t_get_info(t_3c_overl_int_gw_RI, nfull_total=bounds_3c)
         ctr_bounds_mo(:, 1) = [1, bounds_3c(2)]
         ctr_bounds_mo(:, 2) = mo_bounds

         CALL dbcsr_t_contract(dbcsr_scalar(prefac(1)), t_RI, t_3c_overl_int_gw_RI, dbcsr_scalar(0.0_dp), &
                               t_3c_ctr_RI, &
                               contract_1=[2], notcontract_1=[1], &
                               contract_2=[1], notcontract_2=[2, 3], &
                               map_1=[1], map_2=[2, 3], &
                               bounds_3=ctr_bounds_mo, &
                               unit_nr=unit_nr)

      ENDIF

      CALL dbcsr_t_get_info(t_3c_overl_int_gw_AO, nfull_total=bounds_3c)
      ctr_bounds_mo(:, 1) = [1, bounds_3c(2)]
      ctr_bounds_mo(:, 2) = mo_bounds

      CALL dbcsr_t_contract(dbcsr_scalar(prefac(2)), t_AO, t_3c_overl_int_gw_AO, dbcsr_scalar(0.0_dp), &
                            t_3c_ctr_AO, &
                            contract_1=[2], notcontract_1=[1], &
                            contract_2=[1], notcontract_2=[2, 3], &
                            map_1=[1], map_2=[2, 3], &
                            bounds_3=ctr_bounds_mo, &
                            unit_nr=unit_nr)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param t3c_1 ...
!> \param t3c_2 ...
!> \param vec_sigma ...
!> \param mo_offset ...
!> \param mo_bounds ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE trace_sigma_gw(t3c_1, t3c_2, vec_sigma, mo_offset, mo_bounds, para_env)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t3c_1, t3c_2
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: vec_Sigma
      INTEGER, INTENT(IN)                                :: mo_offset
      INTEGER, DIMENSION(2), INTENT(IN)                  :: mo_bounds
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'trace_sigma_gw'

      INTEGER                                            :: blk, handle, n, n_end, n_end_block, &
                                                            n_start, n_start_block
      INTEGER, DIMENSION(1)                              :: trace_shape
      INTEGER, DIMENSION(2)                              :: mo_bounds_off
      INTEGER, DIMENSION(3)                              :: boff, bsize, ind
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: block_1, block_2
      REAL(KIND=dp), &
         DIMENSION(mo_bounds(2)-mo_bounds(1)+1)          :: vec_Sigma_prv
      TYPE(dbcsr_t_iterator_type)                        :: iter
      TYPE(dbcsr_t_type)                                 :: t3c_1_redist

      CALL timeset(routineN, handle)

      CALL dbcsr_t_create(t3c_2, t3c_1_redist)
      CALL dbcsr_t_copy(t3c_1, t3c_1_redist, order=[2, 1, 3], move_data=.TRUE.)

      vec_Sigma_prv = 0.0_dp

      CALL dbcsr_t_iterator_start(iter, t3c_1_redist)
      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, ind, blk, blk_size=bsize, blk_offset=boff)
         CALL dbcsr_t_get_block(t3c_1_redist, ind, block_1, found)
         CPASSERT(found)
         CALL dbcsr_t_get_block(t3c_2, ind, block_2, found)
         IF (.NOT. found) CYCLE

         IF (boff(3) < mo_bounds(1)) THEN
            n_start_block = mo_bounds(1) - boff(3) + 1
            n_start = 1
         ELSE
            n_start_block = 1
            n_start = boff(3) - mo_bounds(1) + 1
         ENDIF

         IF (boff(3) + bsize(3) - 1 > mo_bounds(2)) THEN
            n_end_block = mo_bounds(2) - boff(3) + 1
            n_end = mo_bounds(2) - mo_bounds(1) + 1
         ELSE
            n_end_block = bsize(3)
            n_end = boff(3) + bsize(3) - mo_bounds(1)
         ENDIF

         trace_shape(1) = SIZE(block_1, 1)*SIZE(block_1, 2)
         vec_Sigma_prv(n_start:n_end) = &
            vec_Sigma_prv(n_start:n_end) + &
            (/(DOT_PRODUCT(RESHAPE(block_1(:, :, n), trace_shape), &
                           RESHAPE(block_2(:, :, n), trace_shape)), &
               n=n_start_block, n_end_block)/)
         DEALLOCATE (block_1, block_2)
      ENDDO
      CALL dbcsr_t_iterator_stop(iter)

      CALL dbcsr_t_destroy(t3c_1_redist)

      CALL mp_sum(vec_Sigma_prv, para_env%group)

      mo_bounds_off = mo_bounds - mo_offset + 1
      vec_Sigma(mo_bounds_off(1):mo_bounds_off(2)) = &
         vec_Sigma(mo_bounds_off(1):mo_bounds_off(2)) + vec_Sigma_prv

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_greens_fct_occ ...
!> \param mat_greens_fct_virt ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_mo_coeff_occ_scaled ...
!> \param fm_mo_coeff_virt_scaled ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param Eigenval ...
!> \param nmo ...
!> \param eps_filter ...
!> \param e_fermi ...
!> \param tau ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE compute_Greens_function_time(mat_greens_fct_occ, mat_greens_fct_virt, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                           fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                           fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, Eigenval, nmo, &
                                           eps_filter, e_fermi, tau, para_env)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_greens_fct_occ, mat_greens_fct_virt
      TYPE(cp_fm_type), POINTER :: fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
         fm_mo_coeff_virt_scaled, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: nmo
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter, e_fermi, tau
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Greens_function_time'

      INTEGER                                            :: handle, i_global, iiB, jjB, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: stabilize_exp

      CALL timeset(routineN, handle)

      CALL mp_sync(para_env%group)

      ! get info of fm_mo_coeff_occ
      CALL cp_fm_get_info(matrix=fm_mo_coeff_occ, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! Multiply the occupied and the virtual MO coefficients with the factor exp((-e_i-e_F)*tau/2).
      ! Then, we simply get the sum over all occ states and virt. states by a simple matrix-matrix
      ! multiplication.

      stabilize_exp = 70.0_dp

      ! first, the occ
      DO jjB = 1, nrow_local
         DO iiB = 1, ncol_local
            i_global = col_indices(iiB)

            IF (ABS(tau*0.5_dp*(Eigenval(i_global) - e_fermi)) < stabilize_exp) THEN
               fm_mo_coeff_occ_scaled%local_data(jjB, iiB) = &
                  fm_mo_coeff_occ%local_data(jjB, iiB)*EXP(tau*0.5_dp*(Eigenval(i_global) - e_fermi))
            ELSE
               fm_mo_coeff_occ_scaled%local_data(jjB, iiB) = 0.0_dp
            END IF

         END DO
      END DO

      ! the same for virt
      DO jjB = 1, nrow_local
         DO iiB = 1, ncol_local
            i_global = col_indices(iiB)

            IF (ABS(tau*0.5_dp*(Eigenval(i_global) - e_fermi)) < stabilize_exp) THEN
               fm_mo_coeff_virt_scaled%local_data(jjB, iiB) = &
                  fm_mo_coeff_virt%local_data(jjB, iiB)*EXP(-tau*0.5_dp*(Eigenval(i_global) - e_fermi))
            ELSE
               fm_mo_coeff_virt_scaled%local_data(jjB, iiB) = 0.0_dp
            END IF

         END DO
      END DO

      CALL mp_sync(para_env%group)

      CALL cp_gemm(transa="N", transb="T", m=nmo, n=nmo, k=nmo, alpha=1.0_dp, &
                   matrix_a=fm_mo_coeff_occ_scaled, matrix_b=fm_mo_coeff_occ_scaled, beta=0.0_dp, &
                   matrix_c=fm_scaled_dm_occ_tau)

      CALL cp_gemm(transa="N", transb="T", m=nmo, n=nmo, k=nmo, alpha=1.0_dp, &
                   matrix_a=fm_mo_coeff_virt_scaled, matrix_b=fm_mo_coeff_virt_scaled, beta=0.0_dp, &
                   matrix_c=fm_scaled_dm_virt_tau)

      CALL dbcsr_set(mat_greens_fct_occ, 0.0_dp)

      CALL copy_fm_to_dbcsr(fm_scaled_dm_occ_tau, &
                            mat_greens_fct_occ, &
                            keep_sparsity=.FALSE.)

      CALL dbcsr_filter(mat_greens_fct_occ, eps_filter)

      CALL dbcsr_set(mat_greens_fct_virt, 0.0_dp)

      CALL copy_fm_to_dbcsr(fm_scaled_dm_virt_tau, &
                            mat_greens_fct_virt, &
                            keep_sparsity=.FALSE.)

      CALL dbcsr_filter(mat_greens_fct_virt, eps_filter)

      CALL timestop(handle)

   END SUBROUTINE compute_Greens_function_time

END MODULE rpa_gw

