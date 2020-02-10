!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
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
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add_on_diag, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, dbcsr_filter, &
        dbcsr_get_info, dbcsr_init_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_release_p, dbcsr_scalar, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_copy_matrix_to_tensor, dbcsr_t_create, &
        dbcsr_t_destroy, dbcsr_t_get_block, dbcsr_t_iterator_blocks_left, &
        dbcsr_t_iterator_next_block, dbcsr_t_iterator_start, dbcsr_t_iterator_stop, &
        dbcsr_t_iterator_type, dbcsr_t_type
   USE input_constants,                 ONLY: gw_pade_approx,&
                                              gw_two_pole_model,&
                                              ri_rpa_g0w0_crossing_bisection,&
                                              ri_rpa_g0w0_crossing_newton,&
                                              ri_rpa_g0w0_crossing_none,&
                                              ri_rpa_g0w0_crossing_z_shot
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_release,&
                                              kpoint_sym_create,&
                                              kpoint_type
   USE mathconstants,                   ONLY: fourpi,&
                                              pi,&
                                              twopi
   USE message_passing,                 ONLY: mp_sum,&
                                              mp_sync
   USE mp2_types,                       ONLY: mp2_type
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
   USE rpa_gw_im_time_util,             ONLY: get_tensor_3c_overl_int_gw
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw'

   PUBLIC :: allocate_matrices_gw_im_time, allocate_matrices_gw, GW_matrix_operations, GW_postprocessing, &
             deallocate_matrices_gw_im_time, deallocate_matrices_gw

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param nmo ...
!> \param num_integ_points_freq ...
!> \param unit_nr ...
!> \param RI_blk_sizes ...
!> \param do_ic_model ...
!> \param para_env ...
!> \param fm_mat_W ...
!> \param fm_mat_Q ...
!> \param mo_coeff ...
!> \param t_3c_overl_int_gw_RI ...
!> \param t_3c_overl_int_gw_AO ...
!> \param t_3c_overl_nnP_ic ...
!> \param t_3c_overl_nnP_ic_reflected ...
!> \param matrix_s ...
!> \param mat_W ...
!> \param t_3c_overl_int ...
!> \param qs_env ...
!> \param gw_corr_lev_occ_beta ...
!> \param gw_corr_lev_virt_beta ...
!> \param homo_beta ...
!> \param mo_coeff_beta ...
!> \param t_3c_overl_int_gw_RI_beta ...
!> \param t_3c_overl_int_gw_AO_beta ...
!> \param t_3c_overl_nnP_ic_beta ...
!> \param t_3c_overl_nnP_ic_reflected_beta ...
! **************************************************************************************************
   SUBROUTINE allocate_matrices_gw_im_time(gw_corr_lev_occ, gw_corr_lev_virt, homo, nmo, &
                                           num_integ_points_freq, unit_nr, &
                                           RI_blk_sizes, do_ic_model, &
                                           para_env, fm_mat_W, fm_mat_Q, &
                                           mo_coeff, &
                                           t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                           t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected, &
                                           matrix_s, mat_W, t_3c_overl_int, &
                                           qs_env, &
                                           gw_corr_lev_occ_beta, gw_corr_lev_virt_beta, homo_beta, &
                                           mo_coeff_beta, &
                                           t_3c_overl_int_gw_RI_beta, t_3c_overl_int_gw_AO_beta, &
                                           t_3c_overl_nnP_ic_beta, t_3c_overl_nnP_ic_reflected_beta)

      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt, homo, &
                                                            nmo, num_integ_points_freq, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: RI_blk_sizes
      LOGICAL, INTENT(IN)                                :: do_ic_model
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fm_mat_W
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q, mo_coeff
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_gw_RI, &
                                                            t_3c_overl_int_gw_AO, &
                                                            t_3c_overl_nnP_ic, &
                                                            t_3c_overl_nnP_ic_reflected
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: mat_W
      TYPE(dbcsr_t_type), DIMENSION(:, :)                :: t_3c_overl_int
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN), OPTIONAL                      :: gw_corr_lev_occ_beta, &
                                                            gw_corr_lev_virt_beta, homo_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: mo_coeff_beta
      TYPE(dbcsr_t_type), OPTIONAL :: t_3c_overl_int_gw_RI_beta, t_3c_overl_int_gw_AO_beta, &
         t_3c_overl_nnP_ic_beta, t_3c_overl_nnP_ic_reflected_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_matrices_gw_im_time', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, jquad
      LOGICAL                                            :: my_open_shell

      CALL timeset(routineN, handle)

      my_open_shell = .FALSE.
      IF (PRESENT(gw_corr_lev_occ_beta) .AND. PRESENT(gw_corr_lev_virt_beta) .AND. PRESENT(homo_beta) .AND. &
          PRESENT(mo_coeff_beta) .AND. PRESENT(t_3c_overl_int_gw_RI_beta) .AND. PRESENT(t_3c_overl_int_gw_AO_beta) .AND. &
          PRESENT(t_3c_overl_nnP_ic_beta) &
          .AND. PRESENT(t_3c_overl_nnP_ic_reflected_beta)) THEN
         my_open_shell = .TRUE.

      END IF

      CALL get_tensor_3c_overl_int_gw(t_3c_overl_int, &
                                      t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                      mo_coeff, matrix_s, &
                                      gw_corr_lev_occ, gw_corr_lev_virt, homo, nmo, &
                                      para_env, &
                                      do_ic_model, &
                                      t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected, &
                                      qs_env, unit_nr)

      IF (my_open_shell) THEN

         CALL get_tensor_3c_overl_int_gw(t_3c_overl_int, &
                                         t_3c_overl_int_gw_RI_beta, t_3c_overl_int_gw_AO_beta, &
                                         mo_coeff_beta, matrix_s, &
                                         gw_corr_lev_occ_beta, gw_corr_lev_virt_beta, homo_beta, nmo, &
                                         para_env, &
                                         do_ic_model, &
                                         t_3c_overl_nnP_ic_beta, t_3c_overl_nnP_ic_reflected_beta, &
                                         qs_env, unit_nr)

      END IF

      NULLIFY (fm_mat_W)
      ALLOCATE (fm_mat_W(num_integ_points_freq))

      DO jquad = 1, num_integ_points_freq

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
!> \param num_integ_points_freq ...
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
!> \param vec_Sigma_c_gw_beta ...
!> \param gw_corr_lev_occ_beta ...
!> \param homo_beta ...
!> \param Eigenval_beta ...
!> \param vec_Sigma_x_gw_beta ...
!> \param Eigenval_last_beta ...
!> \param Eigenval_scf_beta ...
!> \param vec_W_gw_beta ...
!> \param fm_mat_S_gw_work_beta ...
!> \param fm_mat_S_gw_beta ...
! **************************************************************************************************
   SUBROUTINE allocate_matrices_gw(vec_Sigma_c_gw, color_rpa_group, dimen_nm_gw, &
                                   gw_corr_lev_occ, gw_corr_lev_virt, homo, &
                                   nmo, num_integ_group, num_integ_points_freq, unit_nr, &
                                   gw_corr_lev_tot, num_fit_points, omega_max_fit, &
                                   do_minimax_quad, do_periodic, do_ri_Sigma_x, my_do_gw, &
                                   first_cycle_periodic_correction, &
                                   a_scaling, Eigenval, tj, vec_omega_fit_gw, vec_Sigma_x_gw, &
                                   delta_corr, Eigenval_last, Eigenval_scf, vec_W_gw, &
                                   fm_mat_S_gw, fm_mat_S_gw_work, &
                                   para_env, mp2_env, kpoints, nkp, nkp_self_energy, &
                                   do_kpoints_cubic_RPA, do_kpoints_from_Gamma, &
                                   vec_Sigma_c_gw_beta, gw_corr_lev_occ_beta, homo_beta, Eigenval_beta, &
                                   vec_Sigma_x_gw_beta, Eigenval_last_beta, Eigenval_scf_beta, vec_W_gw_beta, &
                                   fm_mat_S_gw_work_beta, fm_mat_S_gw_beta)

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(OUT)                 :: vec_Sigma_c_gw
      INTEGER, INTENT(IN) :: color_rpa_group, dimen_nm_gw, gw_corr_lev_occ, gw_corr_lev_virt, &
         homo, nmo, num_integ_group, num_integ_points_freq, unit_nr
      INTEGER, INTENT(INOUT)                             :: gw_corr_lev_tot, num_fit_points
      REAL(KIND=dp)                                      :: omega_max_fit
      LOGICAL, INTENT(IN)                                :: do_minimax_quad, do_periodic, &
                                                            do_ri_Sigma_x, my_do_gw
      LOGICAL, INTENT(OUT) :: first_cycle_periodic_correction
      REAL(KIND=dp), INTENT(IN)                          :: a_scaling
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: vec_omega_fit_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: vec_Sigma_x_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: delta_corr, Eigenval_last, Eigenval_scf, &
                                                            vec_W_gw
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S_gw, fm_mat_S_gw_work
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(OUT)                               :: nkp, nkp_self_energy
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA, &
                                                            do_kpoints_from_Gamma
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(OUT), OPTIONAL       :: vec_Sigma_c_gw_beta
      INTEGER, INTENT(IN), OPTIONAL                      :: gw_corr_lev_occ_beta, homo_beta
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: Eigenval_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT), OPTIONAL                           :: vec_Sigma_x_gw_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: Eigenval_last_beta, Eigenval_scf_beta, &
                                                            vec_W_gw_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mat_S_gw_work_beta, fm_mat_S_gw_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_matrices_gw', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iquad, jquad, n_level_gw, &
                                                            n_level_gw_ref
      LOGICAL                                            :: my_open_shell
      REAL(KIND=dp)                                      :: omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vec_omega_gw

      CALL timeset(routineN, handle)

      my_open_shell = .FALSE.
      IF (PRESENT(vec_Sigma_c_gw_beta) .AND. PRESENT(gw_corr_lev_occ_beta) .AND. PRESENT(homo_beta) .AND. &
          PRESENT(Eigenval_beta) .AND. PRESENT(vec_Sigma_x_gw_beta) .AND. PRESENT(Eigenval_last_beta) .AND. &
          PRESENT(Eigenval_scf_beta) .AND. PRESENT(vec_W_gw_beta) .AND. PRESENT(fm_mat_S_gw_work_beta) .AND. &
          PRESENT(fm_mat_S_gw_beta)) THEN
         my_open_shell = .TRUE.
      END IF

      gw_corr_lev_tot = gw_corr_lev_occ + gw_corr_lev_virt

      ! fill the omega_frequency vector
      ALLOCATE (vec_omega_gw(num_integ_points_freq))
      vec_omega_gw = 0.0_dp

      DO jquad = 1, num_integ_points_freq
         IF (do_minimax_quad) THEN
            omega = tj(jquad)
         ELSE
            omega = a_scaling/TAN(tj(jquad))
         END IF
         vec_omega_gw(jquad) = omega
      END DO

      ! determine number of fit points in the interval [0,w_max] for virt, or [-w_max,0] for occ
      num_fit_points = 0

      DO jquad = 1, num_integ_points_freq
         IF (vec_omega_gw(jquad) < omega_max_fit) THEN
            num_fit_points = num_fit_points + 1
         END IF
      END DO

      IF (mp2_env%ri_g0w0%analytic_continuation == gw_pade_approx) THEN
         IF (mp2_env%ri_g0w0%nparam_pade > num_fit_points) THEN
            IF (unit_nr > 0) &
               CPWARN("Pade approximation: more parameters than data points. Reset # of parameters.")
            mp2_env%ri_g0w0%nparam_pade = num_fit_points
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T74,I7)") &
               "Number of pade parameters:", mp2_env%ri_g0w0%nparam_pade
         ENDIF
      ENDIF

      ! create new arrays containing omega values at which we calculate vec_Sigma_c_gw
      ALLOCATE (vec_omega_fit_gw(num_fit_points))

      ! fill the omega vector with frequencies, where we calculate the self-energy
      iquad = 0
      DO jquad = 1, num_integ_points_freq
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
      ALLOCATE (vec_Sigma_c_gw(gw_corr_lev_tot, num_fit_points, nkp_self_energy))
      vec_Sigma_c_gw = (0.0_dp, 0.0_dp)

      IF (my_open_shell) THEN
         ALLOCATE (vec_Sigma_c_gw_beta(gw_corr_lev_tot, num_fit_points, nkp_self_energy))
         vec_Sigma_c_gw_beta = (0.0_dp, 0.0_dp)
      END IF

      ALLOCATE (Eigenval_scf(nmo))
      Eigenval_scf(:) = Eigenval(:)

      ALLOCATE (Eigenval_last(nmo))
      Eigenval_last(:) = Eigenval(:)

      ! in the case of HF_diag approach of X. Blase (PRB 83, 115103 (2011), Sec. IV), we subtract the
      ! XC potential and add exact exchange
      IF (mp2_env%ri_g0w0%hf_like_ev_start) THEN
         DO n_level_gw = 1, gw_corr_lev_tot
            n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
            Eigenval(n_level_gw_ref) = Eigenval(n_level_gw_ref) + &
                                       mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(n_level_gw_ref, 1, 1)
         END DO
      END IF

      ! Eigenval for beta
      IF (my_open_shell) THEN
         ALLOCATE (Eigenval_scf_beta(nmo))
         Eigenval_scf_beta(:) = Eigenval_beta(:)

         ALLOCATE (Eigenval_last_beta(nmo))
         Eigenval_last_beta(:) = Eigenval_beta(:)

         ! in the case of HF_diag approach of X. Blase (PRB 83, 115103 (2011), Sec. IV), we subtract the
         ! XC potential and add exact exchange
         IF (mp2_env%ri_g0w0%hf_like_ev_start) THEN
            DO n_level_gw = 1, gw_corr_lev_tot
               n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
               Eigenval_beta(n_level_gw_ref) = Eigenval_beta(n_level_gw_ref) + &
                                               mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(n_level_gw_ref, 2, 1)
            END DO
         END IF
      END IF

      IF (do_periodic) THEN

         ALLOCATE (delta_corr(1 + homo - gw_corr_lev_occ:homo + gw_corr_lev_virt))
         delta_corr(:) = 0.0_dp

         first_cycle_periodic_correction = .TRUE.

      END IF

      IF (do_ri_Sigma_x) THEN
         ALLOCATE (vec_Sigma_x_gw(nmo, nkp_self_energy))
         vec_Sigma_x_gw = 0.0_dp

         IF (my_open_shell) THEN
            ALLOCATE (vec_Sigma_x_gw_beta(nmo, nkp_self_energy))
            vec_Sigma_x_gw_beta = 0.0_dp
         END IF
      END IF

      IF (my_do_gw) THEN

         ! minimax grids not implemented for O(N^4) GW
         CPASSERT(.NOT. do_minimax_quad)

         ! create temporary matrix to store B*([1+Q(iw')]^-1-1), has the same size as B
         NULLIFY (fm_mat_S_gw_work)
         CALL cp_fm_create(fm_mat_S_gw_work, fm_mat_S_gw%matrix_struct)
         CALL cp_fm_set_all(matrix=fm_mat_S_gw_work, alpha=0.0_dp)

         IF (my_open_shell) THEN
            NULLIFY (fm_mat_S_gw_work_beta)
            CALL cp_fm_create(fm_mat_S_gw_work_beta, fm_mat_S_gw%matrix_struct)
            CALL cp_fm_set_all(matrix=fm_mat_S_gw_work_beta, alpha=0.0_dp)
         END IF

         ALLOCATE (vec_W_gw(dimen_nm_gw))
         vec_W_gw = 0.0_dp

         IF (my_open_shell) THEN
            ALLOCATE (vec_W_gw_beta(dimen_nm_gw))
            vec_W_gw_beta = 0.0_dp
         END IF

         ! in case we do RI for Sigma_x, we calculate Sigma_x right here
         IF (do_ri_Sigma_x) THEN

            CALL get_vec_sigma_x(vec_Sigma_x_gw, nmo, fm_mat_S_gw, para_env, num_integ_group, color_rpa_group, &
                                 homo, gw_corr_lev_occ, mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1, 1))

            IF (my_open_shell) THEN
               CALL get_vec_sigma_x(vec_Sigma_x_gw_beta, nmo, fm_mat_S_gw_beta, para_env, num_integ_group, &
                                    color_rpa_group, homo_beta, gw_corr_lev_occ_beta, &
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_vec_sigma_x', &
         routineP = moduleN//':'//routineN

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
!> \param fm_mat_S_gw_work_beta ...
!> \param vec_W_gw_beta ...
!> \param vec_Sigma_c_gw_beta ...
!> \param Eigenval_last_beta ...
!> \param Eigenval_scf_beta ...
!> \param vec_Sigma_x_gw_beta ...
! **************************************************************************************************
   SUBROUTINE deallocate_matrices_gw(fm_mat_S_gw_work, vec_W_gw, vec_Sigma_c_gw, vec_omega_fit_gw, &
                                     vec_Sigma_x_minus_vxc_gw, Eigenval_last, &
                                     Eigenval_scf, do_periodic, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, kpoints, &
                                     do_ri_Sigma_x, vec_Sigma_x_gw, my_do_gw, &
                                     fm_mat_S_gw_work_beta, vec_W_gw_beta, &
                                     vec_Sigma_c_gw_beta, Eigenval_last_beta, Eigenval_scf_beta, vec_Sigma_x_gw_beta)

      TYPE(cp_fm_type), POINTER                          :: fm_mat_S_gw_work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: vec_W_gw
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: vec_Sigma_c_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: vec_omega_fit_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: vec_Sigma_x_minus_vxc_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: Eigenval_last, Eigenval_scf
      LOGICAL, INTENT(IN)                                :: do_periodic
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      TYPE(kpoint_type), POINTER                         :: kpoints
      LOGICAL, INTENT(IN)                                :: do_ri_Sigma_x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: vec_Sigma_x_gw
      LOGICAL, INTENT(IN)                                :: my_do_gw
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mat_S_gw_work_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: vec_W_gw_beta
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT), OPTIONAL     :: vec_Sigma_c_gw_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: Eigenval_last_beta, Eigenval_scf_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT), OPTIONAL                         :: vec_Sigma_x_gw_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_matrices_gw', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: my_open_shell

      CALL timeset(routineN, handle)

      my_open_shell = .FALSE.
      IF (PRESENT(fm_mat_S_gw_work_beta) .AND. PRESENT(vec_W_gw_beta) .AND. PRESENT(vec_Sigma_c_gw_beta) .AND. &
          PRESENT(Eigenval_last_beta) .AND. PRESENT(Eigenval_scf_beta) .AND. PRESENT(vec_Sigma_x_gw_beta)) THEN
         my_open_shell = .TRUE.
      END IF

      IF (my_do_gw) THEN
         CALL cp_fm_release(fm_mat_S_gw_work)
         DEALLOCATE (vec_Sigma_x_minus_vxc_gw)
         DEALLOCATE (vec_W_gw)
         IF (my_open_shell) THEN
            CALL cp_fm_release(fm_mat_S_gw_work_beta)
            DEALLOCATE (vec_W_gw_beta)
         END IF
      END IF

      DEALLOCATE (vec_Sigma_c_gw)
      DEALLOCATE (vec_omega_fit_gw)
      DEALLOCATE (Eigenval_last)
      DEALLOCATE (Eigenval_scf)
      IF (my_open_shell) THEN
         DEALLOCATE (vec_Sigma_c_gw_beta)
         DEALLOCATE (Eigenval_last_beta)
         DEALLOCATE (Eigenval_scf_beta)
      END IF

      IF (do_periodic) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_berry_re_mo_mo)
         CALL dbcsr_deallocate_matrix_set(matrix_berry_im_mo_mo)
         CALL kpoint_release(kpoints)
      END IF
      IF (do_ri_Sigma_x) THEN
         DEALLOCATE (vec_Sigma_x_gw)
         IF (my_open_shell) THEN
            DEALLOCATE (vec_Sigma_x_gw_beta)
         END IF
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
!> \param t_3c_overl_int_gw_RI ...
!> \param t_3c_overl_int_gw_AO ...
!> \param t_3c_overl_nnP_ic ...
!> \param t_3c_overl_nnP_ic_reflected ...
!> \param mat_W ...
!> \param ikp_local ...
!> \param cfm_mat_W_kp_tau ...
!> \param t_3c_overl_int_gw_RI_beta ...
!> \param t_3c_overl_int_gw_AO_beta ...
!> \param t_3c_overl_nnP_ic_beta ...
!> \param t_3c_overl_nnP_ic_reflected_beta ...
! **************************************************************************************************
   SUBROUTINE deallocate_matrices_gw_im_time(weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, do_ic_model, do_kpoints_cubic_RPA, &
                                             fm_mat_W, t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                             t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected, mat_W, &
                                             ikp_local, cfm_mat_W_kp_tau, &
                                             t_3c_overl_int_gw_RI_beta, t_3c_overl_int_gw_AO_beta, &
                                             t_3c_overl_nnP_ic_beta, t_3c_overl_nnP_ic_reflected_beta)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w
      LOGICAL, INTENT(IN)                                :: do_ic_model, do_kpoints_cubic_RPA
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fm_mat_W
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_gw_RI, &
                                                            t_3c_overl_int_gw_AO, &
                                                            t_3c_overl_nnP_ic, &
                                                            t_3c_overl_nnP_ic_reflected
      TYPE(dbcsr_type), POINTER                          :: mat_W
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: ikp_local
      TYPE(cp_cfm_p_type), DIMENSION(:, :), POINTER      :: cfm_mat_W_kp_tau
      TYPE(dbcsr_t_type), OPTIONAL :: t_3c_overl_int_gw_RI_beta, t_3c_overl_int_gw_AO_beta, &
         t_3c_overl_nnP_ic_beta, t_3c_overl_nnP_ic_reflected_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_matrices_gw_im_time', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikp, jquad
      LOGICAL                                            :: my_open_shell

      CALL timeset(routineN, handle)

      my_open_shell = .FALSE.
      IF (PRESENT(t_3c_overl_int_gw_RI_beta) .AND. PRESENT(t_3c_overl_int_gw_AO_beta) .AND. PRESENT(t_3c_overl_nnP_ic_beta) &
          .AND. PRESENT(t_3c_overl_nnP_ic_reflected_beta)) THEN
         my_open_shell = .TRUE.
      END IF

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

      CALL dbcsr_t_destroy(t_3c_overl_int_gw_RI)
      CALL dbcsr_t_destroy(t_3c_overl_int_gw_AO)
      IF (PRESENT(t_3c_overl_int_gw_RI_beta)) CALL dbcsr_t_destroy(t_3c_overl_int_gw_RI_beta)
      IF (PRESENT(t_3c_overl_int_gw_AO_beta)) CALL dbcsr_t_destroy(t_3c_overl_int_gw_AO_beta)
      IF (do_ic_model) THEN
         CALL dbcsr_t_destroy(t_3c_overl_nnP_ic)
         CALL dbcsr_t_destroy(t_3c_overl_nnP_ic_reflected)
         IF (PRESENT(t_3c_overl_nnP_ic_beta)) CALL dbcsr_t_destroy(t_3c_overl_nnP_ic_beta)
         IF (PRESENT(t_3c_overl_nnP_ic_reflected_beta)) CALL dbcsr_t_destroy(t_3c_overl_nnP_ic_reflected_beta)
      ENDIF

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
!> \param do_gw_im_time_mix_minimax_clenshaw ...
!> \param first_cycle_periodic_correction ...
!> \param fermi_level_offset ...
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
!> \param vec_Sigma_c_gw_beta ...
!> \param gw_corr_lev_occ_beta ...
!> \param homo_beta ...
!> \param Eigenval_beta ...
!> \param vec_W_gw_beta ...
!> \param fm_mat_S_gw_beta ...
!> \param fm_mat_S_gw_work_beta ...
! **************************************************************************************************
   SUBROUTINE GW_matrix_operations(vec_Sigma_c_gw, dimen_nm_gw, dimen_RI, gw_corr_lev_occ, &
                                   gw_corr_lev_virt, homo, jquad, nmo, num_fit_points, num_integ_points, &
                                   do_bse, do_im_time, do_periodic, do_gw_im_time_mix_minimax_clenshaw, &
                                   first_cycle_periodic_correction, fermi_level_offset, &
                                   omega, Eigenval, delta_corr, tau_tj, tj, vec_omega_fit_gw, &
                                   vec_W_gw, wj, weights_cos_tf_w_to_t, fm_mat_W, fm_mat_L, &
                                   fm_mat_Q, fm_mat_Q_static_bse, fm_mat_R_gw, fm_mat_S_gw, &
                                   fm_mat_S_gw_work, fm_mat_work, mo_coeff, para_env, &
                                   para_env_RPA, matrix_berry_im_mo_mo, matrix_berry_re_mo_mo, &
                                   kpoints, qs_env, mp2_env, do_kpoints_cubic_RPA, do_kpoints_from_Gamma, &
                                   vec_Sigma_c_gw_beta, gw_corr_lev_occ_beta, homo_beta, Eigenval_beta, &
                                   vec_W_gw_beta, fm_mat_S_gw_beta, fm_mat_S_gw_work_beta)

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: vec_Sigma_c_gw
      INTEGER, INTENT(IN)                                :: dimen_nm_gw, dimen_RI, gw_corr_lev_occ, &
                                                            gw_corr_lev_virt, homo, jquad, nmo, &
                                                            num_fit_points, num_integ_points
      LOGICAL, INTENT(IN) :: do_bse, do_im_time, do_periodic, do_gw_im_time_mix_minimax_clenshaw
      LOGICAL, INTENT(INOUT) :: first_cycle_periodic_correction
      REAL(KIND=dp), INTENT(IN)                          :: fermi_level_offset
      REAL(KIND=dp), INTENT(INOUT)                       :: omega
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: delta_corr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tj, vec_omega_fit_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: vec_W_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: wj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: weights_cos_tf_w_to_t
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fm_mat_W
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_mat_L
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q, fm_mat_Q_static_bse, &
                                                            fm_mat_R_gw, fm_mat_S_gw, &
                                                            fm_mat_S_gw_work, fm_mat_work, mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_im_mo_mo, &
                                                            matrix_berry_re_mo_mo
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA, &
                                                            do_kpoints_from_Gamma
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT), OPTIONAL     :: vec_Sigma_c_gw_beta
      INTEGER, INTENT(IN), OPTIONAL                      :: gw_corr_lev_occ_beta, homo_beta
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: Eigenval_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: vec_W_gw_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mat_S_gw_beta, fm_mat_S_gw_work_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'GW_matrix_operations', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, iquad, j_global, &
                                                            jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_open_shell
      REAL(KIND=dp)                                      :: tau, weight

      CALL timeset(routineN, handle)

      my_open_shell = .FALSE.
      IF (PRESENT(vec_Sigma_c_gw_beta) .AND. PRESENT(gw_corr_lev_occ_beta) .AND. PRESENT(homo_beta) &
          .AND. PRESENT(Eigenval_beta) .AND. PRESENT(vec_W_gw_beta) .AND. PRESENT(fm_mat_S_gw_beta) &
          .AND. PRESENT(fm_mat_S_gw_work_beta)) THEN
         my_open_shell = .TRUE.
      END IF

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
                                          mp2_env%ri_g0w0%kp_grid, homo, nmo, gw_corr_lev_occ, &
                                          gw_corr_lev_virt, omega, mo_coeff, Eigenval, &
                                          matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                          first_cycle_periodic_correction, kpoints, &
                                          mp2_env%ri_g0w0%do_mo_coeff_gamma, &
                                          mp2_env%ri_g0w0%num_kp_grids, mp2_env%ri_g0w0%eps_kpoint, &
                                          mp2_env%ri_g0w0%do_extra_kpoints, &
                                          mp2_env%ri_g0w0%do_aux_bas_gw, mp2_env%ri_g0w0%frac_aux_mos)
         END IF

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

         CALL calc_vec_W_gw(vec_Sigma_c_gw, dimen_nm_gw, dimen_RI, gw_corr_lev_occ, homo, jquad, nmo, num_fit_points, &
                            do_periodic, fermi_level_offset, omega, Eigenval, delta_corr, vec_omega_fit_gw, vec_W_gw, &
                            wj, fm_mat_Q, fm_mat_S_gw, fm_mat_S_gw_work)

         IF (my_open_shell) THEN
            CALL calc_vec_W_gw(vec_Sigma_c_gw_beta, dimen_nm_gw, dimen_RI, gw_corr_lev_occ_beta, homo_beta, jquad, nmo, &
                               num_fit_points, do_periodic, fermi_level_offset, omega, Eigenval_beta, delta_corr, &
                               vec_omega_fit_gw, vec_W_gw_beta, wj, fm_mat_Q, fm_mat_S_gw_beta, fm_mat_S_gw_work_beta)
         END IF

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

         IF (do_gw_im_time_mix_minimax_clenshaw) THEN

            CALL cp_fm_to_fm(fm_mat_Q, fm_mat_W(jquad)%matrix)

         ELSE

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
   SUBROUTINE calc_vec_W_gw(vec_Sigma_c_gw, dimen_nm_gw, dimen_RI, gw_corr_lev_occ, homo, jquad, nmo, num_fit_points, &
                            do_periodic, fermi_level_offset, omega, Eigenval, delta_corr, vec_omega_fit_gw, vec_W_gw, &
                            wj, fm_mat_Q, fm_mat_S_gw, fm_mat_S_gw_work)

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: vec_Sigma_c_gw
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_vec_W_gw', routineP = moduleN//':'//routineN
      COMPLEX(KIND=dp), PARAMETER                        :: im_unit = (0.0_dp, 1.0_dp)

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
               (1.0_dp/(im_unit*(omega + omega_i) + e_fermi - Eigenval(m_global)) + &
                1.0_dp/(im_unit*(-omega + omega_i) + e_fermi - Eigenval(m_global)))
         END DO

      END DO

      CALL timestop(handle)

   END SUBROUTINE calc_vec_W_gw

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
!> \param num_integ_points_freq ...
!> \param num_points_corr ...
!> \param unit_nr ...
!> \param do_apply_ic_corr_to_gw ...
!> \param do_im_time ...
!> \param do_periodic ...
!> \param do_gw_im_time_mix_minimax_clenshaw ...
!> \param do_ri_Sigma_x ...
!> \param first_cycle_periodic_correction ...
!> \param e_fermi ...
!> \param eps_filter ...
!> \param fermi_level_offset ...
!> \param stabilize_exp ...
!> \param delta_corr ...
!> \param Eigenval ...
!> \param Eigenval_last ...
!> \param Eigenval_scf ...
!> \param tau_tj ...
!> \param tj ...
!> \param wj ...
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
!> \param iter_ev_sc ...
!> \param vec_Sigma_c_gw_beta ...
!> \param gw_corr_lev_occ_beta ...
!> \param gw_corr_lev_virt_beta ...
!> \param homo_beta ...
!> \param e_fermi_beta ...
!> \param Eigenval_beta ...
!> \param Eigenval_last_beta ...
!> \param Eigenval_scf_beta ...
!> \param vec_Sigma_x_gw_beta ...
!> \param ic_corr_list_beta ...
!> \param fm_mo_coeff_occ_beta ...
!> \param fm_mo_coeff_virt_beta ...
!> \param t_3c_overl_int_gw_RI_beta ...
!> \param t_3c_overl_int_gw_AO_beta ...
! **************************************************************************************************
   SUBROUTINE GW_postprocessing(vec_Sigma_c_gw, count_ev_sc_GW, gw_corr_lev_occ, &
                                gw_corr_lev_tot, gw_corr_lev_virt, homo, &
                                nmo, num_fit_points, num_integ_points, num_integ_points_freq, &
                                num_points_corr, unit_nr, do_apply_ic_corr_to_gw, do_im_time, &
                                do_periodic, do_gw_im_time_mix_minimax_clenshaw, do_ri_Sigma_x, &
                                first_cycle_periodic_correction, e_fermi, eps_filter, &
                                fermi_level_offset, stabilize_exp, delta_corr, Eigenval, &
                                Eigenval_last, Eigenval_scf, tau_tj, tj, wj, &
                                vec_omega_fit_gw, vec_Sigma_x_gw, ic_corr_list, &
                                weights_cos_tf_t_to_w, weights_sin_tf_t_to_w, &
                                fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, fm_mo_coeff_occ, &
                                fm_mo_coeff_virt, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, &
                                mo_coeff, fm_mat_W, para_env, para_env_RPA, mat_dm, mat_SinvVSinv, &
                                t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, matrix_berry_im_mo_mo, &
                                matrix_berry_re_mo_mo, mat_W, matrix_s, &
                                kpoints, mp2_env, qs_env, &
                                nkp_self_energy, do_kpoints_cubic_RPA, Eigenval_kp, Eigenval_scf_kp, iter_ev_sc, &
                                vec_Sigma_c_gw_beta, gw_corr_lev_occ_beta, gw_corr_lev_virt_beta, homo_beta, &
                                e_fermi_beta, Eigenval_beta, Eigenval_last_beta, Eigenval_scf_beta, &
                                vec_Sigma_x_gw_beta, ic_corr_list_beta, fm_mo_coeff_occ_beta, fm_mo_coeff_virt_beta, &
                                t_3c_overl_int_gw_RI_beta, t_3c_overl_int_gw_AO_beta)

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: vec_Sigma_c_gw
      INTEGER, INTENT(IN) :: count_ev_sc_GW, gw_corr_lev_occ, gw_corr_lev_tot, gw_corr_lev_virt, &
         homo, nmo, num_fit_points, num_integ_points, num_integ_points_freq
      INTEGER                                            :: num_points_corr
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN) :: do_apply_ic_corr_to_gw, do_im_time, do_periodic, &
         do_gw_im_time_mix_minimax_clenshaw, do_ri_Sigma_x
      LOGICAL, INTENT(INOUT) :: first_cycle_periodic_correction
      REAL(KIND=dp), INTENT(IN)                          :: e_fermi, eps_filter, fermi_level_offset, &
                                                            stabilize_exp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: delta_corr
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: Eigenval_last, Eigenval_scf, tau_tj, tj, &
                                                            wj, vec_omega_fit_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: vec_Sigma_x_gw
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ic_corr_list
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: weights_cos_tf_t_to_w, &
                                                            weights_sin_tf_t_to_w
      TYPE(cp_fm_type), POINTER :: fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
         fm_mo_coeff_occ, fm_mo_coeff_virt, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, mo_coeff
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fm_mat_W
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_dm, mat_SinvVSinv
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_gw_RI, &
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
      INTEGER, INTENT(INOUT)                             :: iter_ev_sc
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT), OPTIONAL     :: vec_Sigma_c_gw_beta
      INTEGER, INTENT(IN), OPTIONAL                      :: gw_corr_lev_occ_beta, &
                                                            gw_corr_lev_virt_beta, homo_beta
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: e_fermi_beta
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: Eigenval_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: Eigenval_last_beta, Eigenval_scf_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT), OPTIONAL                         :: vec_Sigma_x_gw_beta
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: ic_corr_list_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mo_coeff_occ_beta, &
                                                            fm_mo_coeff_virt_beta
      TYPE(dbcsr_t_type), OPTIONAL                       :: t_3c_overl_int_gw_RI_beta, &
                                                            t_3c_overl_int_gw_AO_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'GW_postprocessing', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: crossing_search, handle, ikp, &
                                                            max_iter_fit, n_level_gw, num_poles
      LOGICAL                                            :: check_fit, my_open_shell, &
                                                            remove_neg_virt_energies
      REAL(KIND=dp)                                      :: stop_crit
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: m_value, m_value_beta, vec_gw_energ, &
         vec_gw_energ_beta, vec_gw_energ_error_fit, vec_gw_energ_error_fit_beta, z_value, &
         z_value_beta

      CALL timeset(routineN, handle)

      my_open_shell = .FALSE.
      IF (PRESENT(vec_Sigma_c_gw_beta) .AND. PRESENT(gw_corr_lev_occ_beta) .AND. PRESENT(gw_corr_lev_virt_beta) .AND. &
          PRESENT(homo_beta) .AND. PRESENT(e_fermi_beta) .AND. PRESENT(Eigenval_beta) .AND. PRESENT(Eigenval_last_beta) &
          .AND. PRESENT(Eigenval_scf_beta) .AND. PRESENT(vec_Sigma_x_gw_beta) .AND. &
          PRESENT(ic_corr_list_beta) .AND. PRESENT(fm_mo_coeff_occ_beta) .AND. PRESENT(fm_mo_coeff_virt_beta) .AND. &
          PRESENT(t_3c_overl_int_gw_RI_beta) .AND. PRESENT(t_3c_overl_int_gw_AO_beta)) THEN
         my_open_shell = .TRUE. ! todo: this should be refactored (and any similar occurence of my_open_shell)
      END IF

      ! postprocessing for cubic scaling GW calculation
      IF (do_im_time .AND. .NOT. do_kpoints_cubic_RPA) THEN
         num_points_corr = mp2_env%ri_g0w0%num_omega_points

         CALL compute_self_energy_cubic_gw(num_integ_points, num_integ_points_freq, nmo, tau_tj, tj, wj, &
                                           matrix_s, fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
                                           fm_mo_coeff_virt_scaled, fm_scaled_dm_occ_tau, &
                                           fm_scaled_dm_virt_tau, Eigenval, vec_omega_fit_gw, eps_filter, &
                                           e_fermi, fermi_level_offset, fm_mat_W, &
                                           gw_corr_lev_tot, gw_corr_lev_occ, gw_corr_lev_virt, homo, &
                                           count_ev_sc_GW, &
                                           t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                           mat_W, mat_SinvVSinv, mat_dm, stabilize_exp, &
                                           weights_cos_tf_t_to_w, weights_sin_tf_t_to_w, vec_Sigma_c_gw, &
                                           do_periodic, do_gw_im_time_mix_minimax_clenshaw, &
                                           num_points_corr, delta_corr, qs_env, para_env, para_env_RPA, &
                                           mp2_env, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                           first_cycle_periodic_correction, kpoints, num_fit_points, mo_coeff, &
                                           do_ri_Sigma_x, vec_Sigma_x_gw, unit_nr)

         IF (my_open_shell) THEN

            CALL compute_self_energy_cubic_gw(num_integ_points, num_integ_points_freq, nmo, tau_tj, tj, wj, &
                                              matrix_s, fm_mo_coeff_occ_beta, fm_mo_coeff_virt_beta, fm_mo_coeff_occ_scaled, &
                                              fm_mo_coeff_virt_scaled, fm_scaled_dm_occ_tau, &
                                              fm_scaled_dm_virt_tau, Eigenval_beta, vec_omega_fit_gw, eps_filter, &
                                              e_fermi_beta, fermi_level_offset, fm_mat_W, &
                                              gw_corr_lev_tot, gw_corr_lev_occ_beta, gw_corr_lev_virt_beta, homo_beta, &
                                              count_ev_sc_GW, &
                                              t_3c_overl_int_gw_RI_beta, t_3c_overl_int_gw_AO_beta, &
                                              mat_W, mat_SinvVSinv, mat_dm, stabilize_exp, &
                                              weights_cos_tf_t_to_w, weights_sin_tf_t_to_w, vec_Sigma_c_gw_beta, &
                                              do_periodic, do_gw_im_time_mix_minimax_clenshaw, &
                                              num_points_corr, delta_corr, qs_env, para_env, para_env_RPA, &
                                              mp2_env, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                              first_cycle_periodic_correction, kpoints, num_fit_points, mo_coeff, &
                                              do_ri_Sigma_x, vec_Sigma_x_gw_beta, unit_nr, do_beta=.TRUE.)

         END IF

      END IF

      IF (.NOT. do_im_time) THEN

         CALL mp_sum(vec_Sigma_c_gw, para_env%group)

      END IF

      IF (do_periodic .AND. mp2_env%ri_g0w0%do_average_deg_levels) THEN

         CALL average_degenerate_levels(vec_Sigma_c_gw, Eigenval(1 + homo - gw_corr_lev_occ:homo + gw_corr_lev_virt), &
                                        mp2_env%ri_g0w0%eps_eigenval)
         IF (my_open_shell) THEN
            CALL average_degenerate_levels(vec_Sigma_c_gw_beta, &
                                           Eigenval_beta(1 + homo_beta - gw_corr_lev_occ_beta: &
                                                         homo_beta + gw_corr_lev_virt_beta), &
                                           mp2_env%ri_g0w0%eps_eigenval)
         END IF
      END IF

      IF (my_open_shell .AND. .NOT. do_im_time) THEN
         CALL mp_sum(vec_Sigma_c_gw_beta, para_env%group)
      END IF

      CALL mp_sync(para_env%group)

      stop_crit = 1.0e-7
      num_poles = mp2_env%ri_g0w0%num_poles
      max_iter_fit = mp2_env%ri_g0w0%max_iter_fit
      check_fit = mp2_env%ri_g0w0%check_fit
      crossing_search = mp2_env%ri_g0w0%crossing_search

      ! arrays storing the correlation self-energy, stat. error and z-shot value
      ALLOCATE (vec_gw_energ(gw_corr_lev_tot))
      vec_gw_energ = 0.0_dp
      ALLOCATE (vec_gw_energ_error_fit(gw_corr_lev_tot))
      vec_gw_energ_error_fit = 0.0_dp
      ALLOCATE (z_value(gw_corr_lev_tot))
      z_value = 0.0_dp
      ALLOCATE (m_value(gw_corr_lev_tot))
      m_value = 0.0_dp

      ! the same for beta
      IF (my_open_shell) THEN
         ALLOCATE (vec_gw_energ_beta(gw_corr_lev_tot))
         vec_gw_energ_beta = 0.0_dp
         ALLOCATE (vec_gw_energ_error_fit_beta(gw_corr_lev_tot))
         vec_gw_energ_error_fit_beta = 0.0_dp
         ALLOCATE (z_value_beta(gw_corr_lev_tot))
         z_value_beta = 0.0_dp
         ALLOCATE (m_value_beta(gw_corr_lev_tot))
         m_value_beta = 0.0_dp
      END IF

      ! for the normal code for molecules or Gamma only: nkp = 1
      DO ikp = 1, nkp_self_energy

         IF (do_kpoints_cubic_RPA) THEN

            vec_gw_energ_error_fit = 0.0_dp
            vec_gw_energ = 0.0_dp
            z_value = 0.0_dp
            m_value = 0.0_dp

            CALL get_eigenval_for_conti(Eigenval, Eigenval_scf, Eigenval_kp, Eigenval_scf_kp, kpoints, &
                                        ikp, iter_ev_sc, my_open_shell)
         END IF

         ! fit the self-energy on imaginary frequency axis and evaluate the fit on the MO energy of the SCF
         DO n_level_gw = 1, gw_corr_lev_tot
            ! processes perform different fits
            IF (MODULO(n_level_gw, para_env%num_pe) /= para_env%mepos) CYCLE

            SELECT CASE (mp2_env%ri_g0w0%analytic_continuation)
            CASE (gw_two_pole_model)
               CALL fit_and_continuation_2pole(vec_gw_energ, vec_gw_energ_error_fit, vec_omega_fit_gw, &
                                               z_value, m_value, vec_Sigma_c_gw(:, :, ikp), &
                                               mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1, ikp), &
                                               Eigenval, Eigenval_scf, n_level_gw, gw_corr_lev_occ, num_poles, &
                                               num_fit_points, max_iter_fit, crossing_search, homo, check_fit, stop_crit, &
                                               fermi_level_offset, &
                                               do_im_time .AND. (.NOT. do_gw_im_time_mix_minimax_clenshaw))

            CASE (gw_pade_approx)
               CALL continuation_pade(vec_gw_energ, vec_omega_fit_gw, &
                                      z_value, m_value, vec_Sigma_c_gw(:, :, ikp), &
                                      mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1, ikp), &
                                      Eigenval, Eigenval_scf, n_level_gw, gw_corr_lev_occ, mp2_env%ri_g0w0%nparam_pade, &
                                      num_fit_points, crossing_search, homo, check_fit, &
                                      fermi_level_offset, &
                                      do_im_time .AND. (.NOT. do_gw_im_time_mix_minimax_clenshaw))

            CASE DEFAULT
               CPABORT("Only two-model and Pade approximation are implemented.")
            END SELECT

            IF (my_open_shell) THEN
               SELECT CASE (mp2_env%ri_g0w0%analytic_continuation)
               CASE (gw_two_pole_model)
                  CALL fit_and_continuation_2pole( &
                     vec_gw_energ_beta, vec_gw_energ_error_fit_beta, vec_omega_fit_gw, &
                     z_value_beta, m_value_beta, vec_Sigma_c_gw_beta(:, :, ikp), &
                     mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2, ikp), &
                     Eigenval_beta, Eigenval_scf_beta, n_level_gw, &
                     gw_corr_lev_occ_beta, num_poles, &
                     num_fit_points, max_iter_fit, crossing_search, homo_beta, check_fit, stop_crit, &
                     fermi_level_offset, do_im_time)
               CASE (gw_pade_approx)
                  CALL continuation_pade(vec_gw_energ_beta, vec_omega_fit_gw, &
                                         z_value_beta, m_value_beta, vec_Sigma_c_gw_beta(:, :, ikp), &
                                         mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2, ikp), &
                                         Eigenval_beta, Eigenval_scf_beta, n_level_gw, &
                                         gw_corr_lev_occ_beta, mp2_env%ri_g0w0%nparam_pade, &
                                         num_fit_points, crossing_search, homo_beta, check_fit, &
                                         fermi_level_offset, do_im_time)
               CASE DEFAULT
                  CPABORT("Only two-model and Pade approximation are implemented.")
               END SELECT

            END IF

         END DO ! n_level_gw

         CALL mp_sum(vec_gw_energ_error_fit, para_env%group)
         CALL mp_sum(vec_gw_energ, para_env%group)
         CALL mp_sum(z_value, para_env%group)
         CALL mp_sum(m_value, para_env%group)

         IF (my_open_shell) THEN
            CALL mp_sum(vec_gw_energ_error_fit_beta, para_env%group)
            CALL mp_sum(vec_gw_energ_beta, para_env%group)
            CALL mp_sum(z_value_beta, para_env%group)
            CALL mp_sum(m_value_beta, para_env%group)
         END IF

         remove_neg_virt_energies = mp2_env%ri_g0w0%remove_neg_virt_energies

         ! print the quasiparticle energies and update Eigenval in case you do eigenvalue self-consistent GW
         IF (my_open_shell) THEN

            CALL print_and_update_for_ev_sc( &
               vec_gw_energ, vec_gw_energ_error_fit, &
               z_value, m_value, mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1, ikp), Eigenval, &
               Eigenval_last, Eigenval_scf, gw_corr_lev_occ, gw_corr_lev_virt, gw_corr_lev_tot, &
               count_ev_sc_GW, crossing_search, homo, nmo, unit_nr, mp2_env%ri_g0w0%print_gw_details, &
               remove_neg_virt_energies, ikp, nkp_self_energy, kpoints, do_alpha=.TRUE.)

            CALL print_and_update_for_ev_sc( &
               vec_gw_energ_beta, vec_gw_energ_error_fit_beta, &
               z_value_beta, m_value_beta, mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2, ikp), Eigenval_beta, &
               Eigenval_last_beta, Eigenval_scf_beta, gw_corr_lev_occ_beta, gw_corr_lev_virt_beta, gw_corr_lev_tot, &
               count_ev_sc_GW, crossing_search, homo_beta, nmo, unit_nr, mp2_env%ri_g0w0%print_gw_details, &
               remove_neg_virt_energies, ikp, nkp_self_energy, kpoints, do_beta=.TRUE.)

            IF (do_apply_ic_corr_to_gw .AND. count_ev_sc_GW == 1) THEN

               CALL apply_ic_corr(Eigenval, Eigenval_scf, ic_corr_list, &
                                  gw_corr_lev_occ, gw_corr_lev_virt, gw_corr_lev_tot, &
                                  homo, nmo, unit_nr, do_alpha=.TRUE.)

               CALL apply_ic_corr(Eigenval_beta, Eigenval_scf_beta, ic_corr_list_beta, &
                                  gw_corr_lev_occ_beta, gw_corr_lev_virt_beta, gw_corr_lev_tot, &
                                  homo_beta, nmo, unit_nr, do_beta=.TRUE.)

            END IF

         ELSE

            CALL print_and_update_for_ev_sc( &
               vec_gw_energ, vec_gw_energ_error_fit, &
               z_value, m_value, mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1, ikp), Eigenval, &
               Eigenval_last, Eigenval_scf, gw_corr_lev_occ, gw_corr_lev_virt, gw_corr_lev_tot, &
               count_ev_sc_GW, crossing_search, homo, nmo, unit_nr, mp2_env%ri_g0w0%print_gw_details, &
               remove_neg_virt_energies, ikp, nkp_self_energy, kpoints)

            IF (do_apply_ic_corr_to_gw .AND. count_ev_sc_GW == 1) THEN

               CALL apply_ic_corr(Eigenval, Eigenval_scf, ic_corr_list, &
                                  gw_corr_lev_occ, gw_corr_lev_virt, gw_corr_lev_tot, &
                                  homo, nmo, unit_nr)

            END IF

         END IF

      END DO ! ikp

      DEALLOCATE (z_value)
      DEALLOCATE (m_value)
      DEALLOCATE (vec_gw_energ)
      DEALLOCATE (vec_gw_energ_error_fit)
      IF (my_open_shell) THEN
         DEALLOCATE (z_value_beta)
         DEALLOCATE (m_value_beta)
         DEALLOCATE (vec_gw_energ_beta)
         DEALLOCATE (vec_gw_energ_error_fit_beta)
      END IF

      CALL timestop(handle)

   END SUBROUTINE

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_periodic_correction', &
         routineP = moduleN//':'//routineN

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
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: homo
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: omega

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_eps_head_Berry', &
         routineP = moduleN//':'//routineN

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_berry_phase', &
         routineP = moduleN//':'//routineN

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

      END IF

      IF (do_aux_bas) THEN

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

         ! multiply eigenvectors with inverse sqrt of eigenvalues
         DO i_col_local = 1, ncol_local

            col_index = col_indices(i_col_local)

            fm_mat_scaled_eigv_P%local_data(:, i_col_local) = &
               fm_mat_scaled_eigv_P%local_data(:, i_col_local)*evals_P_sqrt_inv(col_index)

         END DO

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
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: eps_inv_head
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      INTEGER, INTENT(IN)                                :: homo, gw_corr_lev_occ, gw_corr_lev_virt
      TYPE(cp_para_env_type), OPTIONAL, POINTER          :: para_env_RPA
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
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: eps_head
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_eps_inv_head', &
         routineP = moduleN//':'//routineN

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
   SUBROUTINE average_degenerate_levels(vec_Sigma_c_gw, Eigenval_DFT, eps_eigenval)
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: vec_Sigma_c_gw
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
!> \param iter_ev_sc ...
!> \param my_open_shell ...
! **************************************************************************************************
   SUBROUTINE get_eigenval_for_conti(Eigenval, Eigenval_scf, Eigenval_kp, Eigenval_scf_kp, kpoints, ikp, iter_ev_sc, my_open_shell)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval, Eigenval_scf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: Eigenval_kp, Eigenval_scf_kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: ikp, iter_ev_sc
      LOGICAL, INTENT(IN)                                :: my_open_shell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_eigenval_for_conti', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, jkp, nkp, nmo, nspin
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      TYPE(mo_set_type), POINTER                         :: mo_set

      CALL timeset(routineN, handle)

      ! only implemented for non-evscGW
      CPASSERT(iter_ev_sc == 1)

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
!> \param vec_gw_energ_error_fit ...
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
!> \param max_iter_fit ...
!> \param crossing_search ...
!> \param homo ...
!> \param check_fit ...
!> \param stop_crit ...
!> \param fermi_level_offset ...
!> \param do_gw_im_time ...
! **************************************************************************************************
   SUBROUTINE fit_and_continuation_2pole(vec_gw_energ, vec_gw_energ_error_fit, vec_omega_fit_gw, &
                                         z_value, m_value, vec_Sigma_c_gw, vec_Sigma_x_minus_vxc_gw, &
                                         Eigenval, Eigenval_scf, n_level_gw, gw_corr_lev_occ, num_poles, &
                                         num_fit_points, max_iter_fit, crossing_search, homo, check_fit, stop_crit, &
                                         fermi_level_offset, do_gw_im_time)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: vec_gw_energ, vec_gw_energ_error_fit, &
                                                            vec_omega_fit_gw, z_value, m_value
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: vec_Sigma_c_gw
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vec_Sigma_x_minus_vxc_gw, Eigenval, &
                                                            Eigenval_scf
      INTEGER, INTENT(IN)                                :: n_level_gw, gw_corr_lev_occ, num_poles, &
                                                            num_fit_points, max_iter_fit, &
                                                            crossing_search, homo
      LOGICAL, INTENT(IN)                                :: check_fit
      REAL(KIND=dp), INTENT(IN)                          :: stop_crit, fermi_level_offset
      LOGICAL, INTENT(IN)                                :: do_gw_im_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fit_and_continuation_2pole', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: func_val, im_unit, one, re_unit, rho1, &
                                                            zero
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: dLambda, dLambda_2, Lambda, &
                                                            Lambda_without_offset, vec_b_gw, &
                                                            vec_b_gw_copy
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: mat_A_gw, mat_B_gw
      INTEGER                                            :: handle4, ierr, iii, iiter, info, &
                                                            integ_range, jjj, jquad, kkk, &
                                                            n_level_gw_ref, num_var, output_unit, &
                                                            xpos
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipiv
      LOGICAL                                            :: could_exit
      REAL(KIND=dp) :: chi2, chi2_old, delta, deriv_val_real, e_fermi, gw_energ, Ldown, &
         level_energ_GW, Lup, range_step, ScalParam, sign_occ_virt, stat_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Lambda_Im, Lambda_Re, stat_errors, &
                                                            vec_N_gw, vec_omega_fit_gw_sign
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_N_gw

      output_unit = cp_logger_get_default_io_unit()

      im_unit = (0.0_dp, 1.0_dp)
      re_unit = (1.0_dp, 0.0_dp)

      num_var = 2*num_poles + 1
      ALLOCATE (Lambda(num_var))
      Lambda = (0.0_dp, 0.0_dp)
      ALLOCATE (Lambda_without_offset(num_var))
      Lambda_without_offset = (0.0_dp, 0.0_dp)
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
         Lambda(iii) = Lambda_Re(iii) + im_unit*Lambda_Im(iii)
      END DO

      CALL calc_chi2(chi2_old, Lambda, vec_Sigma_c_gw, vec_omega_fit_gw_sign, num_poles, &
                     num_fit_points, n_level_gw)

      ALLOCATE (mat_A_gw(num_poles + 1, num_poles + 1))
      ALLOCATE (vec_b_gw(num_poles + 1))
      ALLOCATE (ipiv(num_poles + 1))
      mat_A_gw = (0.0_dp, 0.0_dp)
      vec_b_gw = 0.0_dp

      DO iii = 1, num_poles + 1
         mat_A_gw(iii, 1) = (1.0_dp, 0.0_dp)
      END DO
      integ_range = num_fit_points/num_poles
      DO kkk = 1, num_poles + 1
         xpos = (kkk - 1)*integ_range + 1
         xpos = MIN(xpos, num_fit_points)
         ! calculate coefficient at this point
         DO iii = 1, num_poles
            jjj = iii*2
            func_val = (1.0_dp, 0.0_dp)/(im_unit*vec_omega_fit_gw_sign(xpos) - &
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
            Lambda(iii) = Lambda_Re(iii) + im_unit*Lambda_Im(iii)
         END DO
         dLambda = (0.0_dp, 0.0_dp)

         DO kkk = 1, num_fit_points
            func_val = Lambda(1)
            DO iii = 1, num_poles
               jjj = iii*2
               func_val = func_val + Lambda(jjj)/(vec_omega_fit_gw_sign(kkk)*im_unit - Lambda(jjj + 1))
            END DO
            dLambda(kkk) = vec_Sigma_c_gw(n_level_gw, kkk) - func_val
         END DO
         rho1 = SUM(dLambda*dLambda)

         ! fill matrix
         mat_B_gw = (0.0_dp, 0.0_dp)
         DO iii = 1, num_fit_points
            mat_B_gw(iii, 1) = 1.0_dp
            mat_B_gw(iii, num_var + 1) = im_unit
         END DO
         DO iii = 1, num_poles
            jjj = iii*2
            DO kkk = 1, num_fit_points
               mat_B_gw(kkk, jjj) = 1.0_dp/(im_unit*vec_omega_fit_gw_sign(kkk) - Lambda(jjj + 1))
               mat_B_gw(kkk, jjj + num_var) = im_unit/(im_unit*vec_omega_fit_gw_sign(kkk) - Lambda(jjj + 1))
               mat_B_gw(kkk, jjj + 1) = Lambda(jjj)/(im_unit*vec_omega_fit_gw_sign(kkk) - Lambda(jjj + 1))**2
               mat_B_gw(kkk, jjj + 1 + num_var) = (-Lambda_Im(jjj) + im_unit*Lambda_Re(jjj))/ &
                                                  (im_unit*vec_omega_fit_gw_sign(kkk) - Lambda(jjj + 1))**2
            END DO
         END DO

         CALL timestop(handle4)

         CALL timeset(routineN//"_fit_matmul_1", handle4)

         one = (1.0_dp, 0.0_dp)
         zero = (0.0_dp, 0.0_dp)
         CALL zgemm('C', 'N', num_var*2, num_var*2, num_fit_points, one, mat_B_gw, num_fit_points, mat_B_gw, num_fit_points, &
                    zero, mat_A_gw, num_var*2)
         CALL timestop(handle4)

         CALL timeset(routineN//"_fit_zgemv_1", handle4)
         CALL zgemv('C', num_fit_points, num_var*2, one, mat_B_gw, num_fit_points, dLambda, 1, &
                    zero, vec_b_gw, 1)

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
            Lambda(iii) = Lambda_Re(iii) + im_unit*Lambda_Im(iii) + vec_b_gw(iii) + vec_b_gw(iii + num_var)
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

      ! print self_energy and fit on the imaginary frequency axis if required
      IF (check_fit) THEN

         IF (output_unit > 0) THEN

            WRITE (output_unit, *) ' '
            WRITE (output_unit, '(T3,A,I5)') 'Check the GW fit for molecular orbital', n_level_gw_ref
            WRITE (output_unit, '(T3,A)') '-------------------------------------------'
            WRITE (output_unit, *)
            WRITE (output_unit, '(T3,5A)') '  omega (i*eV)    ', 'Re(fit) (eV)    ', &
               'Im(fit) (eV)  ', 'Re(Sig_c) (eV)  ', &
               'Im(Sig_c) (eV)'

         END IF

         DO kkk = 1, num_fit_points
            func_val = Lambda(1)
            DO iii = 1, num_poles
               jjj = iii*2
               ! calculate value of the fit function
               func_val = func_val + Lambda(jjj)/(im_unit*vec_omega_fit_gw_sign(kkk) - Lambda(jjj + 1))
            END DO
            WRITE (output_unit, '(1F16.3,4F16.5)') vec_omega_fit_gw_sign(kkk)*evolt, REAL(func_val)*evolt, &
               AIMAG(func_val)*evolt, REAL(vec_Sigma_c_gw(n_level_gw, kkk))*evolt, &
               AIMAG(vec_Sigma_c_gw(n_level_gw, kkk))*evolt
         END DO

         WRITE (output_unit, *) ' '

      END IF

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

      ! either Z-shot or no crossing search for evaluating Sigma_c
      IF (crossing_search == ri_rpa_g0w0_crossing_none) THEN

         ! calculate func val on the real axis
         ! gw_energ = only correlation part of the self energy
         func_val = Lambda(1)
         DO iii = 1, num_poles
            jjj = iii*2
            func_val = func_val + Lambda(jjj)/(Eigenval(n_level_gw_ref) - e_fermi - Lambda(jjj + 1))
         END DO

         gw_energ = REAL(func_val)
         vec_gw_energ(n_level_gw) = gw_energ

      ELSE IF (crossing_search == ri_rpa_g0w0_crossing_z_shot .OR. &
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

      ! Compute error of Sigma_c on real axis according to error propagation

      vec_gw_energ_error_fit(n_level_gw) = 0.0_dp

      DO kkk = 1, num_poles
         vec_gw_energ_error_fit(n_level_gw) = vec_gw_energ_error_fit(n_level_gw) + &
                                              (stat_errors(4*kkk - 1) + stat_errors(4*kkk))* &
                                              ABS(1.0_dp/(Eigenval(n_level_gw_ref) - e_fermi - Lambda(2*kkk + 1)) - &
                                                  1.0_dp/(-Lambda(2*kkk + 1))) + &
                                              (stat_errors(4*kkk + 1) + stat_errors(4*kkk + 2))*ABS(Lambda(2*kkk))* &
                                              ABS(1.0_dp/(Eigenval(n_level_gw_ref) - e_fermi - Lambda(2*kkk + 1))**2 - &
                                                  1.0_dp/(-Lambda(2*kkk + 1))**2)
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
!> \param check_fit ...
!> \param fermi_level_offset ...
!> \param do_gw_im_time ...
! **************************************************************************************************
   SUBROUTINE continuation_pade(vec_gw_energ, vec_omega_fit_gw, &
                                z_value, m_value, vec_Sigma_c_gw, vec_Sigma_x_minus_vxc_gw, &
                                Eigenval, Eigenval_scf, n_level_gw, gw_corr_lev_occ, nparam_pade, &
                                num_fit_points, crossing_search, homo, check_fit, &
                                fermi_level_offset, do_gw_im_time)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: vec_gw_energ
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vec_omega_fit_gw
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: z_value, m_value
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: vec_Sigma_c_gw
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vec_Sigma_x_minus_vxc_gw, Eigenval, &
                                                            Eigenval_scf
      INTEGER, INTENT(IN)                                :: n_level_gw, gw_corr_lev_occ, &
                                                            nparam_pade, num_fit_points, &
                                                            crossing_search, homo
      LOGICAL, INTENT(IN)                                :: check_fit
      REAL(KIND=dp), INTENT(IN)                          :: fermi_level_offset
      LOGICAL, INTENT(IN)                                :: do_gw_im_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'continuation_pade', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: im_unit, re_unit, sigma_c_pade
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: coeff_pade, omega_points_pade, &
                                                            Sigma_c_gw_reorder
      INTEGER                                            :: handle, jquad, n_level_gw_ref, &
                                                            output_unit
      REAL(KIND=dp)                                      :: e_fermi, energy_val, level_energ_GW, &
                                                            sign_occ_virt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vec_omega_fit_gw_sign, &
                                                            vec_omega_fit_gw_sign_reorder

      CALL timeset(routineN, handle)

      output_unit = cp_logger_get_default_io_unit()

      im_unit = (0.0_dp, 1.0_dp)
      re_unit = (1.0_dp, 0.0_dp)

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
      IF (check_fit) THEN
         CALL check_fit_pade(vec_omega_fit_gw_sign, vec_Sigma_c_gw(n_level_gw, :), &
                             nparam_pade, omega_points_pade, coeff_pade, &
                             num_fit_points, n_level_gw_ref, output_unit)
      ENDIF

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
      CASE (ri_rpa_g0w0_crossing_none)
         energy_val = Eigenval(n_level_gw_ref) - e_fermi
         CALL evaluate_pade_function(energy_val, nparam_pade, omega_points_pade, &
                                     coeff_pade, sigma_c_pade)
         vec_gw_energ(n_level_gw) = REAL(sigma_c_pade)

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
                                         start_val=level_energ_GW)
         z_value(n_level_gw) = 1.0_dp
         m_value(n_level_gw) = 0.0_dp

      CASE (ri_rpa_g0w0_crossing_newton)
         CALL get_sigma_c_newton_pade(vec_gw_energ(n_level_gw), Eigenval_scf(n_level_gw_ref), &
                                      vec_Sigma_x_minus_vxc_gw(n_level_gw_ref), e_fermi, &
                                      nparam_pade, omega_points_pade, coeff_pade, &
                                      start_val=level_energ_GW)
         z_value(n_level_gw) = 1.0_dp
         m_value(n_level_gw) = 0.0_dp

      CASE DEFAULT
         CPABORT("Only NONE, Z_SHOT, NEWTON, and BISECTION crossing search implemented.")
      END SELECT

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
   SUBROUTINE get_pade_parameters(y, x, num_fit_points, nparam, xpoints, coeff)

      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: y
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: x
      INTEGER, INTENT(IN)                                :: num_fit_points, nparam
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(INOUT)      :: xpoints, coeff

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_pade_parameters', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: im_unit
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: ypoints
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: g_mat
      INTEGER                                            :: handle, idat, iparam, nstep

      CALL timeset(routineN, handle)

      im_unit = (0.0_dp, 1.0_dp)

      nstep = INT(num_fit_points/(nparam - 1))
      CPASSERT(LBOUND(x, 1) == 1)
      CPASSERT(LBOUND(y, 1) == 1)

      ALLOCATE (ypoints(nparam))
      !omega=i0 is in element x(1)
      idat = 1
      DO iparam = 1, nparam - 1
         xpoints(iparam) = im_unit*x(idat)
         ypoints(iparam) = y(idat)
         idat = idat + nstep
      ENDDO
      xpoints(nparam) = im_unit*x(num_fit_points)
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

      CALL timestop(handle)

   END SUBROUTINE get_pade_parameters

! **************************************************************************************************
!> \brief evalute pade function for a real value x_val
!> \param x_val real value
!> \param nparam number of pade parameters
!> \param xpoints selection of points of the original complex function, i.e. here of Sigma_c(iomega)
!> \param coeff pade coefficients
!> \param func_val function value
! **************************************************************************************************
   SUBROUTINE evaluate_pade_function(x_val, nparam, xpoints, coeff, func_val)

      REAL(KIND=dp), INTENT(IN)                          :: x_val
      INTEGER, INTENT(IN)                                :: nparam
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: xpoints, coeff
      COMPLEX(KIND=dp), INTENT(OUT)                      :: func_val

      CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_pade_function', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: im_unit, re_unit
      INTEGER                                            :: handle, iparam

      CALL timeset(routineN, handle)

      im_unit = (0.0_dp, 1.0_dp)
      re_unit = (1.0_dp, 0.0_dp)

      func_val = re_unit
      DO iparam = nparam, 2, -1
         func_val = re_unit + coeff(iparam)*(re_unit*x_val - xpoints(iparam - 1))/func_val
      ENDDO

      func_val = coeff(1)/func_val

      CALL timestop(handle)

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
   SUBROUTINE get_z_and_m_value_pade(x_val, nparam, xpoints, coeff, z_value, m_value)

      REAL(KIND=dp), INTENT(IN)                          :: x_val
      INTEGER, INTENT(IN)                                :: nparam
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: xpoints, coeff
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: z_value, m_value

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_z_and_m_value_pade', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: denominator, dev_denominator, &
                                                            dev_numerator, dev_val, func_val, &
                                                            im_unit, numerator, re_unit
      INTEGER                                            :: iparam

      im_unit = (0.0_dp, 1.0_dp)
      re_unit = (1.0_dp, 0.0_dp)

      func_val = re_unit
      dev_val = (0.0_dp, 0.0_dp)
      DO iparam = nparam, 2, -1
         numerator = coeff(iparam)*(re_unit*x_val - xpoints(iparam - 1))
         dev_numerator = coeff(iparam)*re_unit
         denominator = func_val
         dev_denominator = dev_val
         dev_val = dev_numerator/denominator - (numerator*dev_denominator)/(denominator**2)
         func_val = re_unit + coeff(iparam)*(re_unit*x_val - xpoints(iparam - 1))/func_val
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
!> \param start_val start value for the quasiparticle iteration
! **************************************************************************************************
   SUBROUTINE get_sigma_c_bisection_pade(gw_energ, Eigenval_scf, Sigma_x_minus_vxc_gw, e_fermi, &
                                         nparam_pade, omega_points_pade, coeff_pade, start_val)

      REAL(KIND=dp), INTENT(OUT)                         :: gw_energ
      REAL(KIND=dp), INTENT(IN)                          :: Eigenval_scf, Sigma_x_minus_vxc_gw, &
                                                            e_fermi
      INTEGER, INTENT(IN)                                :: nparam_pade
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: omega_points_pade, coeff_pade
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: start_val

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_sigma_c_bisection_pade', &
         routineP = moduleN//':'//routineN

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
            CPABORT("Self-consistent quasi-particle solution not found")
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
!> \param start_val start value for the quasiparticle iteration
! **************************************************************************************************
   SUBROUTINE get_sigma_c_newton_pade(gw_energ, Eigenval_scf, Sigma_x_minus_vxc_gw, e_fermi, &
                                      nparam_pade, omega_points_pade, coeff_pade, start_val)

      REAL(KIND=dp), INTENT(OUT)                         :: gw_energ
      REAL(KIND=dp), INTENT(IN)                          :: Eigenval_scf, Sigma_x_minus_vxc_gw, &
                                                            e_fermi
      INTEGER, INTENT(IN)                                :: nparam_pade
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: omega_points_pade, coeff_pade
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: start_val

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_sigma_c_newton_pade', &
         routineP = moduleN//':'//routineN

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
            CPABORT("Self-consistent quasi-particle solution not found")
            EXIT
         ENDIF
      ENDDO

      gw_energ = REAL(sigma_c)

      CALL timestop(handle)

   END SUBROUTINE get_sigma_c_newton_pade

! **************************************************************************************************
!> \brief check "fit" for analytic continuation with pade approximation
!> \param vec_omega_fit_gw_sign ...
!> \param Sigma_c_gw complex Sigma_c(iomega) for n_level_gw
!> \param nparam_pade number of pade parameters
!> \param omega_points_pade selection of frequency points of Sigma_c(iomega)
!> \param coeff_pade pade coefficients
!> \param num_fit_points total number of frequency points for the complex Sigma_c
!> \param n_level_gw_ref n_level_gw+homo-gw_corr_lev_occ
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE check_fit_pade(vec_omega_fit_gw_sign, Sigma_c_gw, &
                             nparam_pade, omega_points_pade, coeff_pade, &
                             num_fit_points, n_level_gw_ref, output_unit)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vec_omega_fit_gw_sign
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: Sigma_c_gw
      INTEGER, INTENT(IN)                                :: nparam_pade
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: omega_points_pade, coeff_pade
      INTEGER, INTENT(IN)                                :: num_fit_points, n_level_gw_ref, &
                                                            output_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_fit_pade', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: func_val, im_unit, re_unit
      INTEGER                                            :: iparam, kkk

      re_unit = (1.0_dp, 0.0_dp)
      im_unit = (0.0_dp, 1.0_dp)

      WRITE (output_unit, *) ' '
      WRITE (output_unit, '(T3,A,I5)') 'Check the GW fit for molecular orbital', n_level_gw_ref
      WRITE (output_unit, '(T3,A)') '-------------------------------------------'
      WRITE (output_unit, *)
      WRITE (output_unit, '(T3,5A)') '  omega (i*eV)    ', 'Re(fit) (eV)    ', &
         'Im(fit) (eV)  ', 'Re(Sig_c) (eV)  ', &
         'Im(Sig_c) (eV)'

      DO kkk = 1, num_fit_points
         func_val = re_unit
         DO iparam = nparam_pade, 2, -1
            func_val = re_unit + coeff_pade(iparam) &
                       *(im_unit*vec_omega_fit_gw_sign(kkk) - omega_points_pade(iparam - 1))/func_val
         ENDDO

         func_val = coeff_pade(1)/func_val

         WRITE (output_unit, '(1F16.3,4F16.5)') vec_omega_fit_gw_sign(kkk)*evolt, REAL(func_val)*evolt, &
            AIMAG(func_val)*evolt, REAL(Sigma_c_gw(kkk))*evolt, &
            AIMAG(Sigma_c_gw(kkk))*evolt
      END DO

      WRITE (output_unit, *) ' '

   END SUBROUTINE check_fit_pade

! **************************************************************************************************
!> \brief Prints the GW stuff to the output and optinally to an external file.
!>        Also updates the eigenvalues for eigenvalue-self-consistent GW
!> \param vec_gw_energ ...
!> \param vec_gw_energ_error_fit ...
!> \param z_value ...
!> \param m_value ...
!> \param vec_Sigma_x_minus_vxc_gw ...
!> \param Eigenval ...
!> \param Eigenval_last ...
!> \param Eigenval_scf ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param gw_corr_lev_tot ...
!> \param count_ev_sc_GW ...
!> \param crossing_search ...
!> \param homo ...
!> \param nmo ...
!> \param unit_nr ...
!> \param print_gw_details ...
!> \param remove_neg_virt_energies ...
!> \param ikp ...
!> \param nkp_self_energy ...
!> \param kpoints ...
!> \param do_alpha ...
!> \param do_beta ...
! **************************************************************************************************
   SUBROUTINE print_and_update_for_ev_sc(vec_gw_energ, vec_gw_energ_error_fit, &
                                         z_value, m_value, vec_Sigma_x_minus_vxc_gw, Eigenval, &
                                         Eigenval_last, Eigenval_scf, gw_corr_lev_occ, gw_corr_lev_virt, gw_corr_lev_tot, &
                                         count_ev_sc_GW, crossing_search, homo, nmo, unit_nr, print_gw_details, &
                                         remove_neg_virt_energies, ikp, nkp_self_energy, kpoints, do_alpha, do_beta)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: vec_gw_energ, vec_gw_energ_error_fit, &
                                                            z_value, m_value
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: vec_Sigma_x_minus_vxc_gw, Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: Eigenval_last, Eigenval_scf
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt, &
                                                            gw_corr_lev_tot, count_ev_sc_GW, &
                                                            crossing_search, homo, nmo, unit_nr
      LOGICAL, INTENT(IN)                                :: print_gw_details, &
                                                            remove_neg_virt_energies
      INTEGER, INTENT(IN)                                :: ikp, nkp_self_energy
      TYPE(kpoint_type), INTENT(IN), POINTER             :: kpoints
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_alpha, do_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_and_update_for_ev_sc', &
         routineP = moduleN//':'//routineN

      CHARACTER(4)                                       :: occ_virt
      INTEGER                                            :: handle, n_level_gw, n_level_gw_ref
      LOGICAL                                            :: do_closed_shell, do_kpoints, &
                                                            is_energy_okay, my_do_alpha, my_do_beta
      REAL(KIND=dp)                                      :: eigen_diff, new_energy

      CALL timeset(routineN, handle)

      IF (PRESENT(do_alpha)) THEN
         my_do_alpha = do_alpha
      ELSE
         my_do_alpha = .FALSE.
      END IF

      IF (PRESENT(do_beta)) THEN
         my_do_beta = do_beta
      ELSE
         my_do_beta = .FALSE.
      END IF

      do_closed_shell = .NOT. (my_do_alpha .OR. my_do_beta)
      do_kpoints = (nkp_self_energy > 1)

      Eigenval_last(:) = Eigenval(:)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, *) ' '

         IF (do_closed_shell) THEN
            WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies'
            WRITE (unit_nr, '(T3,A)') '-------------------------'
         ELSE IF (my_do_alpha) THEN
            WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies of alpha spins'
            WRITE (unit_nr, '(T3,A)') '----------------------------------------'
         ELSE IF (my_do_beta) THEN
            WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies of beta spins'
            WRITE (unit_nr, '(T3,A)') '---------------------------------------'
         END IF

         IF (do_kpoints) THEN
            WRITE (unit_nr, *) ' '
            WRITE (unit_nr, '(T3,A7,I3,A3,I3,A8,3F7.3,A12,3F7.3)') 'Kpoint ', ikp, '  /', nkp_self_energy, &
               '   xkp =', kpoints%xkp(1, ikp), kpoints%xkp(2, ikp), kpoints%xkp(3, ikp), &
               '  and  xkp =', -kpoints%xkp(1, ikp), -kpoints%xkp(2, ikp), -kpoints%xkp(3, ikp)
            WRITE (unit_nr, '(T3,A72)') '(Relative Brillouin zone size: [-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5])'
         END IF

      END IF

      IF (unit_nr > 0 .AND. (.NOT. print_gw_details)) THEN
         WRITE (unit_nr, *) ' '
         WRITE (unit_nr, '(T5,A)') 'Molecular orbital        MO energy after SCF (eV)        G0W0 QP energy (eV)'
      END IF

      IF (unit_nr > 0 .AND. print_gw_details) THEN
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A)') 'The GW quasiparticle energies are calculated according to: '
      END IF

      IF (crossing_search == ri_rpa_g0w0_crossing_none) THEN

         DO n_level_gw = 1, gw_corr_lev_tot
            n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
            new_energy = Eigenval_scf(n_level_gw_ref) + vec_gw_energ(n_level_gw) + &
                         vec_Sigma_x_minus_vxc_gw(n_level_gw_ref)

            is_energy_okay = .TRUE.

            IF (remove_neg_virt_energies .AND. (n_level_gw_ref > homo .AND. new_energy < Eigenval(homo))) THEN
               is_energy_okay = .FALSE.
            END IF

            IF (is_energy_okay) THEN
               Eigenval(n_level_gw_ref) = new_energy
            END IF
         END DO

         IF (unit_nr > 0 .AND. print_gw_details) THEN
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A)') 'E_GW = E_SCF +  Sigc(E_SCF) + Sigx - vxc'
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A)') 'The energy unit of the following table is eV. Sigc_fit is a very conservative'
            WRITE (unit_nr, '(T3,A)') 'estimate of the statistical error of the correlation self-energy caused by the'
            WRITE (unit_nr, '(T3,A)') 'fitting.'
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T14,2A)') 'MO        E_SCF         Sigc', &
               '     Sigc_fit     Sigx-vxc         E_GW'
         END IF

         DO n_level_gw = 1, gw_corr_lev_tot
            n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
            IF (n_level_gw <= gw_corr_lev_occ) THEN
               occ_virt = 'occ'
            ELSE
               occ_virt = 'vir'
            END IF

            IF (unit_nr > 0 .AND. (.NOT. print_gw_details)) THEN
               WRITE (unit_nr, '(T5,I9,3A,2F27.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ')     ', &
                  Eigenval_last(n_level_gw_ref)*evolt, &
                  Eigenval(n_level_gw_ref)*evolt
            END IF

            IF (unit_nr > 0 .AND. print_gw_details) THEN
               WRITE (unit_nr, '(T4,I4,3A,5F13.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ')', &
                  Eigenval_last(n_level_gw_ref)*evolt, &
                  vec_gw_energ(n_level_gw)*evolt, &
                  vec_gw_energ_error_fit(n_level_gw)*evolt, &
                  vec_Sigma_x_minus_vxc_gw(n_level_gw_ref)*evolt, &
                  Eigenval(n_level_gw_ref)*evolt
            END IF
         END DO

         ! z-shot
      ELSE

         DO n_level_gw = 1, gw_corr_lev_tot

            n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ

            new_energy = (Eigenval_scf(n_level_gw_ref) - &
                          m_value(n_level_gw)*Eigenval(n_level_gw_ref) + &
                          vec_gw_energ(n_level_gw) + &
                          vec_Sigma_x_minus_vxc_gw(n_level_gw_ref))* &
                         z_value(n_level_gw)

            is_energy_okay = .TRUE.

            IF (remove_neg_virt_energies .AND. (n_level_gw_ref > homo .AND. new_energy < Eigenval(homo))) THEN
               is_energy_okay = .FALSE.
            END IF

            IF (is_energy_okay) THEN
               Eigenval(n_level_gw_ref) = new_energy
            END IF

         END DO

         IF (unit_nr > 0 .AND. print_gw_details) THEN
            WRITE (unit_nr, '(T3,A)') 'E_GW = E_SCF + Z * ( Sigc(E_SCF) + Sigx - vxc )'
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A)') 'The energy unit of the following table is eV.  Sigc_fit is a very conservative'
            WRITE (unit_nr, '(T3,2A)') 'estimate of the statistical error of the fitting.'
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T13,2A)') 'MO      E_SCF       Sigc', &
               '   Sigc_fit   Sigx-vxc          Z       E_GW'
         END IF

         DO n_level_gw = 1, gw_corr_lev_tot
            n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
            IF (n_level_gw <= gw_corr_lev_occ) THEN
               occ_virt = 'occ'
            ELSE
               occ_virt = 'vir'
            END IF

            IF (unit_nr > 0 .AND. (.NOT. print_gw_details)) THEN
               WRITE (unit_nr, '(T5,I9,3A,2F27.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ')     ', &
                  Eigenval_last(n_level_gw_ref)*evolt, &
                  Eigenval(n_level_gw_ref)*evolt
            END IF

            IF (unit_nr > 0 .AND. print_gw_details) THEN
               WRITE (unit_nr, '(T3,I4,3A,6F11.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ')', &
                  Eigenval_last(n_level_gw_ref)*evolt, &
                  vec_gw_energ(n_level_gw)*evolt, &
                  vec_gw_energ_error_fit(n_level_gw)*evolt, &
                  vec_Sigma_x_minus_vxc_gw(n_level_gw_ref)*evolt, &
                  z_value(n_level_gw), &
                  Eigenval(n_level_gw_ref)*evolt
            END IF
         END DO

         IF (unit_nr > 0) THEN
            IF (do_closed_shell) THEN
               WRITE (unit_nr, '(T3,A)') ' '
               WRITE (unit_nr, '(T3,A,F57.2)') 'GW HOMO-LUMO gap (eV)', (Eigenval(homo + 1) - Eigenval(homo))*evolt
            ELSE IF (my_do_alpha) THEN
               WRITE (unit_nr, '(T3,A)') ' '
               WRITE (unit_nr, '(T3,A,F51.2)') 'Alpha GW HOMO-LUMO gap (eV)', (Eigenval(homo + 1) - Eigenval(homo))*evolt
            ELSE IF (my_do_beta) THEN
               WRITE (unit_nr, '(T3,A)') ' '
               WRITE (unit_nr, '(T3,A,F52.2)') 'Beta GW HOMO-LUMO gap (eV)', (Eigenval(homo + 1) - Eigenval(homo))*evolt
            END IF
         END IF

      END IF ! z-shot vs. no crossing

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) ' '
      END IF

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

      IF ((gw_corr_lev_occ == 0 .OR. gw_corr_lev_virt == 0) .AND. count_ev_sc_GW > 1) THEN
         CALL cp_warn(__LOCATION__, &
                      "Please increase for eigenvalue-self-consistent GW, the number of "// &
                      "corrected occupied and/or virtual orbitals above 0.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_and_update_for_ev_sc

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_mat_N', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: im_unit, re_unit
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: Lambda_tmp
      INTEGER                                            :: handle, num_var
      REAL(KIND=dp)                                      :: chi2, chi2_sum

      CALL timeset(routineN, handle)

      num_var = 2*num_poles + 1
      ALLOCATE (Lambda_tmp(num_var))
      Lambda_tmp = (0.0_dp, 0.0_dp)
      chi2_sum = 0.0_dp
      re_unit = (1.0_dp, 0.0_dp)
      im_unit = (0.0_dp, 1.0_dp)

      !test
      Lambda_tmp(:) = Lambda(:)
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)

      ! Fitting parameters with offset h
      Lambda_tmp(:) = Lambda(:)
      IF (MODULO(i, 2) == 0) THEN
         Lambda_tmp(i/2) = Lambda_tmp(i/2) + h*re_unit
      ELSE
         Lambda_tmp((i + 1)/2) = Lambda_tmp((i + 1)/2) + h*im_unit
      END IF
      IF (MODULO(j, 2) == 0) THEN
         Lambda_tmp(j/2) = Lambda_tmp(j/2) + h*re_unit
      ELSE
         Lambda_tmp((j + 1)/2) = Lambda_tmp((j + 1)/2) + h*im_unit
      END IF
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)
      chi2_sum = chi2_sum + chi2

      IF (MODULO(i, 2) == 0) THEN
         Lambda_tmp(i/2) = Lambda_tmp(i/2) - 2.0_dp*h*re_unit
      ELSE
         Lambda_tmp((i + 1)/2) = Lambda_tmp((i + 1)/2) - 2.0_dp*h*im_unit
      END IF
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)
      chi2_sum = chi2_sum - chi2

      IF (MODULO(j, 2) == 0) THEN
         Lambda_tmp(j/2) = Lambda_tmp(j/2) - 2.0_dp*h*re_unit
      ELSE
         Lambda_tmp((j + 1)/2) = Lambda_tmp((j + 1)/2) - 2.0_dp*h*im_unit
      END IF
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)
      chi2_sum = chi2_sum + chi2

      IF (MODULO(i, 2) == 0) THEN
         Lambda_tmp(i/2) = Lambda_tmp(i/2) + 2.0_dp*h*re_unit
      ELSE
         Lambda_tmp((i + 1)/2) = Lambda_tmp((i + 1)/2) + 2.0_dp*h*im_unit
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
   SUBROUTINE calc_chi2(chi2, Lambda, Sigma_c, vec_omega_fit_gw, num_poles, &
                        num_fit_points, n_level_gw)
      REAL(KIND=dp), INTENT(INOUT)                       :: chi2
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: Lambda
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: Sigma_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: vec_omega_fit_gw
      INTEGER, INTENT(IN)                                :: num_poles, num_fit_points, n_level_gw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_chi2', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: func_val, im_unit
      INTEGER                                            :: handle, iii, jjj, kkk

      CALL timeset(routineN, handle)

      im_unit = (0.0_dp, 1.0_dp)
      chi2 = 0.0_dp
      DO kkk = 1, num_fit_points
         func_val = Lambda(1)
         DO iii = 1, num_poles
            jjj = iii*2
            ! calculate value of the fit function
            func_val = func_val + Lambda(jjj)/(im_unit*vec_omega_fit_gw(kkk) - Lambda(jjj + 1))
         END DO
         chi2 = chi2 + (ABS(Sigma_c(n_level_gw, kkk) - func_val))**2
      END DO

      CALL timestop(handle)

   END SUBROUTINE calc_chi2

! **************************************************************************************************
!> \brief ...
!> \param Eigenval ...
!> \param Eigenval_scf ...
!> \param ic_corr_list ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param gw_corr_lev_tot ...
!> \param homo ...
!> \param nmo ...
!> \param unit_nr ...
!> \param do_alpha ...
!> \param do_beta ...
! **************************************************************************************************
   SUBROUTINE apply_ic_corr(Eigenval, Eigenval_scf, ic_corr_list, &
                            gw_corr_lev_occ, gw_corr_lev_virt, gw_corr_lev_tot, &
                            homo, nmo, unit_nr, do_alpha, do_beta)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval, Eigenval_scf
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: ic_corr_list
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt, &
                                                            gw_corr_lev_tot, homo, nmo, unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_alpha, do_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_ic_corr', routineP = moduleN//':'//routineN

      CHARACTER(4)                                       :: occ_virt
      INTEGER                                            :: handle, n_level_gw, n_level_gw_ref
      LOGICAL                                            :: do_closed_shell, my_do_alpha, my_do_beta
      REAL(KIND=dp)                                      :: eigen_diff

      CALL timeset(routineN, handle)

      IF (PRESENT(do_alpha)) THEN
         my_do_alpha = do_alpha
      ELSE
         my_do_alpha = .FALSE.
      END IF

      IF (PRESENT(do_beta)) THEN
         my_do_beta = do_beta
      ELSE
         my_do_beta = .FALSE.
      END IF

      do_closed_shell = .NOT. (my_do_alpha .OR. my_do_beta)

      ! check the number of input image charge corrected levels
      CPASSERT(SIZE(ic_corr_list) == gw_corr_lev_tot)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, *) ' '

         IF (do_closed_shell) THEN
            WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies with image charge (ic) correction'
            WRITE (unit_nr, '(T3,A)') '-----------------------------------------------------------'
         ELSE IF (my_do_alpha) THEN
            WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies of alpha spins with image charge (ic) correction'
            WRITE (unit_nr, '(T3,A)') '--------------------------------------------------------------------------'
         ELSE IF (my_do_beta) THEN
            WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies of beta spins with image charge (ic) correction'
            WRITE (unit_nr, '(T3,A)') '-------------------------------------------------------------------------'
         END IF

         WRITE (unit_nr, *) ' '

         DO n_level_gw = 1, gw_corr_lev_tot
            n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
            IF (n_level_gw <= gw_corr_lev_occ) THEN
               occ_virt = 'occ'
            ELSE
               occ_virt = 'vir'
            END IF

            WRITE (unit_nr, '(T4,I4,3A,3F21.3)') &
               n_level_gw_ref, ' ( ', occ_virt, ')  ', &
               Eigenval(n_level_gw_ref)*evolt, &
               ic_corr_list(n_level_gw)*evolt, &
               (Eigenval(n_level_gw_ref) + ic_corr_list(n_level_gw))*evolt

         END DO

         WRITE (unit_nr, *) ' '

      END IF

      Eigenval(homo - gw_corr_lev_occ + 1:homo + gw_corr_lev_virt) = Eigenval(homo - gw_corr_lev_occ + 1: &
                                                                              homo + gw_corr_lev_virt) &
                                                                     + ic_corr_list(1:gw_corr_lev_tot)

      Eigenval_scf(homo - gw_corr_lev_occ + 1:homo + gw_corr_lev_virt) = Eigenval_scf(homo - gw_corr_lev_occ + 1: &
                                                                                      homo + gw_corr_lev_virt) &
                                                                         + ic_corr_list(1:gw_corr_lev_tot)

      IF (unit_nr > 0) THEN

         IF (do_closed_shell) THEN
            WRITE (unit_nr, '(T3,A,F52.2)') 'G0W0 IC HOMO-LUMO gap (eV)', Eigenval(homo + 1) - Eigenval(homo)
         ELSE IF (my_do_alpha) THEN
            WRITE (unit_nr, '(T3,A,F46.2)') 'G0W0 Alpha IC HOMO-LUMO gap (eV)', Eigenval(homo + 1) - Eigenval(homo)
         ELSE IF (my_do_beta) THEN
            WRITE (unit_nr, '(T3,A,F47.2)') 'G0W0 Beta IC HOMO-LUMO gap (eV)', Eigenval(homo + 1) - Eigenval(homo)
         END IF

         WRITE (unit_nr, *) ' '

      END IF

      ! for eigenvalue self-consistent GW, all eigenvalues have to be corrected
      ! 1) the occupied; check if there are occupied MOs not being corrected by the IC model
      IF (gw_corr_lev_occ < homo .AND. gw_corr_lev_occ > 0) THEN

         ! calculate average IC contribution for occupied orbitals
         eigen_diff = 0.0_dp

         DO n_level_gw = 1, gw_corr_lev_occ
            eigen_diff = eigen_diff + ic_corr_list(n_level_gw)
         END DO
         eigen_diff = eigen_diff/gw_corr_lev_occ

         ! correct the eigenvalues of the occupied orbitals which have not been corrected by the IC model
         DO n_level_gw = 1, homo - gw_corr_lev_occ
            Eigenval(n_level_gw) = Eigenval(n_level_gw) + eigen_diff
            Eigenval_scf(n_level_gw) = Eigenval_scf(n_level_gw) + eigen_diff
         END DO

      END IF

      ! 2) the virtual: check if there are virtual orbitals not being corrected by the IC model
      IF (gw_corr_lev_virt < nmo - homo .AND. gw_corr_lev_virt > 0) THEN

         ! calculate average IC correction for virtual orbitals
         eigen_diff = 0.0_dp
         DO n_level_gw = gw_corr_lev_occ + 1, gw_corr_lev_tot
            eigen_diff = eigen_diff + ic_corr_list(n_level_gw)
         END DO
         eigen_diff = eigen_diff/gw_corr_lev_virt

         ! correct the eigenvalues of the virtual orbitals which have not been corrected by the IC model
         DO n_level_gw = homo + gw_corr_lev_virt + 1, nmo
            Eigenval(n_level_gw) = Eigenval(n_level_gw) + eigen_diff
            Eigenval_scf(n_level_gw) = Eigenval_scf(n_level_gw) + eigen_diff
         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE apply_ic_corr

! **************************************************************************************************
!> \brief ...
!> \param num_integ_points ...
!> \param num_integ_points_freq ...
!> \param nmo ...
!> \param tau_tj ...
!> \param tj ...
!> \param wj ...
!> \param matrix_s ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_mo_coeff_occ_scaled ...
!> \param fm_mo_coeff_virt_scaled ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param Eigenval ...
!> \param vec_omega_fit_gw ...
!> \param eps_filter ...
!> \param e_fermi ...
!> \param fermi_level_offset ...
!> \param fm_mat_W ...
!> \param gw_corr_lev_tot ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param count_ev_sc_GW ...
!> \param t_3c_overl_int_gw_RI ...
!> \param t_3c_overl_int_gw_AO ...
!> \param mat_W ...
!> \param mat_SinvVSinv ...
!> \param mat_dm ...
!> \param stabilize_exp ...
!> \param weights_cos_tf_t_to_w ...
!> \param weights_sin_tf_t_to_w ...
!> \param vec_Sigma_c_gw ...
!> \param do_periodic ...
!> \param do_gw_im_time_mix_minimax_clenshaw ...
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
!> \param do_beta ...
! **************************************************************************************************
   SUBROUTINE compute_self_energy_cubic_gw(num_integ_points, num_integ_points_freq, nmo, tau_tj, tj, wj, &
                                           matrix_s, fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
                                           fm_mo_coeff_virt_scaled, fm_scaled_dm_occ_tau, &
                                           fm_scaled_dm_virt_tau, Eigenval, vec_omega_fit_gw, eps_filter, &
                                           e_fermi, fermi_level_offset, fm_mat_W, &
                                           gw_corr_lev_tot, gw_corr_lev_occ, gw_corr_lev_virt, homo, count_ev_sc_GW, &
                                           t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                           mat_W, mat_SinvVSinv, mat_dm, stabilize_exp, &
                                           weights_cos_tf_t_to_w, weights_sin_tf_t_to_w, vec_Sigma_c_gw, &
                                           do_periodic, do_gw_im_time_mix_minimax_clenshaw, &
                                           num_points_corr, delta_corr, qs_env, para_env, para_env_RPA, &
                                           mp2_env, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                           first_cycle_periodic_correction, kpoints, num_fit_points, fm_mo_coeff, &
                                           do_ri_Sigma_x, vec_Sigma_x_gw, unit_nr, do_beta)
      INTEGER, INTENT(IN)                                :: num_integ_points, num_integ_points_freq, &
                                                            nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tj, wj
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(cp_fm_type), POINTER :: fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
         fm_mo_coeff_virt_scaled, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval, vec_omega_fit_gw
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter, e_fermi, fermi_level_offset
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fm_mat_W
      INTEGER, INTENT(IN)                                :: gw_corr_lev_tot, gw_corr_lev_occ, &
                                                            gw_corr_lev_virt, homo, count_ev_sc_GW
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_gw_RI, &
                                                            t_3c_overl_int_gw_AO
      TYPE(dbcsr_type), POINTER                          :: mat_W
      TYPE(dbcsr_p_type)                                 :: mat_SinvVSinv, mat_dm
      REAL(KIND=dp), INTENT(IN)                          :: stabilize_exp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: weights_cos_tf_t_to_w, &
                                                            weights_sin_tf_t_to_w
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: vec_Sigma_c_gw
      LOGICAL, INTENT(IN) :: do_periodic, do_gw_im_time_mix_minimax_clenshaw
      INTEGER, INTENT(IN)                                :: num_points_corr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: delta_corr
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      LOGICAL, INTENT(INOUT) :: first_cycle_periodic_correction
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: num_fit_points
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff
      LOGICAL, INTENT(IN)                                :: do_ri_Sigma_x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: vec_Sigma_x_gw
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_self_energy_cubic_gw', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: im_unit, re_im
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: delta_corr_omega
      INTEGER :: gw_lev_end, gw_lev_start, handle, handle3, i_occ_virt, iquad, jquad, mo_end, &
         mo_start, n, n_level_gw, n_level_gw_ref, real_imaginary, unit_nr_prv
      LOGICAL                                            :: do_occ, do_virt, memory_info, my_do_beta
      REAL(KIND=dp) :: ext_scaling, omega, omega_diff, omega_i, omega_sign, sign_occ_virt, &
         t_i_Clenshaw, tau, weight, weight_cos, weight_i, weight_sin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Sigma_c_gw_tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: vec_Sigma_c_gw_cos_omega, &
         vec_Sigma_c_gw_cos_tau, vec_Sigma_c_gw_neg_tau, vec_Sigma_c_gw_pos_tau, &
         vec_Sigma_c_gw_sin_omega, vec_Sigma_c_gw_sin_tau
      TYPE(dbcsr_t_type)                                 :: t_3c_ctr_RI
      TYPE(dbcsr_type), POINTER                          :: mat_greens_fct_occ, mat_greens_fct_virt

      CALL timeset(routineN, handle)

      memory_info = mp2_env%ri_rpa_im_time%memory_info
      IF (memory_info) THEN
         unit_nr_prv = unit_nr
      ELSE
         unit_nr_prv = 0
      ENDIF

      my_do_beta = .FALSE.
      IF (PRESENT(do_beta)) my_do_beta = do_beta

      im_unit = (0.0_dp, 1.0_dp)

      mo_start = homo - gw_corr_lev_occ + 1
      mo_end = homo + gw_corr_lev_virt
      CPASSERT(mo_end - mo_start + 1 == gw_corr_lev_tot)

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
      delta_corr_omega(:, :) = (0.0_dp, 0.0_dp)

      NULLIFY (mat_greens_fct_occ)
      CALL dbcsr_init_p(mat_greens_fct_occ)
      CALL dbcsr_create(matrix=mat_greens_fct_occ, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      NULLIFY (mat_greens_fct_virt)
      CALL dbcsr_init_p(mat_greens_fct_virt)
      CALL dbcsr_create(matrix=mat_greens_fct_virt, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      IF (do_gw_im_time_mix_minimax_clenshaw) THEN

         ! we perform an integration Sigma(iw) = int dw' G(i(w-w')) W(iw') here with jquad: omega' amd iquad: omega
         ! w' is stored in tj with entries > 0. Therefore, we split the integration in two parts for w' > 0 / w' = tj(jquad)
         ! w' < 0 / w' = -tj(jquad)
         DO jquad = 1, 2*num_integ_points_freq

            CALL dbcsr_set(mat_W, 0.0_dp)

            ! W(iw') = W(-iw')
            IF (jquad .LE. num_integ_points_freq) THEN
               CALL copy_fm_to_dbcsr(fm_mat_W(jquad)%matrix, mat_W, keep_sparsity=.FALSE.)
            ELSE
               CALL copy_fm_to_dbcsr(fm_mat_W(jquad - num_integ_points_freq)%matrix, mat_W, keep_sparsity=.FALSE.)
            END IF

            DO iquad = 1, num_fit_points

               IF (jquad .LE. num_integ_points_freq) THEN
                  omega_diff = vec_omega_fit_gw(iquad) - tj(jquad)
               ELSE
                  omega_diff = vec_omega_fit_gw(iquad) + tj(jquad - num_integ_points_freq)
               END IF

               DO real_imaginary = 1, 2

                  ! compute Green's function G(i*(w-/+w')) and store it to fm_scaled_dm_occ_tau and
                  ! fm_scaled_dm_virt_tau. In the occ-matrix, we store the Gf with e_f = e_HOMO + fermi_level_offset,
                  ! in the virt-matrix, we store the Gf with e_f = e_LUMO - fermi_level_offset.
                  ! The Gf is imaginary, so we do a loop real_imaginary = 1, 2
                  CALL compute_Greens_function_freq(mat_greens_fct_occ, mat_greens_fct_virt, fm_mo_coeff, &
                                                    fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                                    fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, Eigenval, fermi_level_offset, &
                                                    nmo, homo, eps_filter, omega_diff, real_imaginary)

                  DO i_occ_virt = 1, 2

                     do_occ = .FALSE.
                     do_virt = .FALSE.
                     IF (i_occ_virt == 1) do_occ = .TRUE.
                     IF (i_occ_virt == 2) do_virt = .TRUE.

                     IF (do_occ) THEN
                        mo_start = homo - gw_corr_lev_occ + 1
                        mo_end = homo
                        ALLOCATE (Sigma_c_gw_tmp(gw_corr_lev_occ))
                     ELSE IF (do_virt) THEN
                        mo_start = homo + 1
                        mo_end = homo + gw_corr_lev_virt
                        ALLOCATE (Sigma_c_gw_tmp(gw_corr_lev_virt))
                     END IF

                     IF (real_imaginary == 1) re_im = (1.0_dp, 0.0_dp)
                     IF (real_imaginary == 2) re_im = (0.0_dp, 1.0_dp)

                     Sigma_c_gw_tmp = 0.0_dp

                     ! 3-center * W(iw')
                     CALL compute_sigma_gw(t_3c_overl_int_gw_AO, t_3c_overl_int_gw_RI, &
                                           mat_greens_fct_occ, mat_W, [1.0_dp, 0.0_dp], &
                                           vec_Sigma_c_gw_neg_tau(:, 1), [mo_start, mo_end], para_env, unit_nr_prv, &
                                           t_3c_ctr_RI=t_3c_ctr_RI, t_3c_ctr_in=.FALSE.)

                     IF (do_occ) THEN

                        ! 3-center * G(i(w-w')) for an occupied level n (other Fermi level than virtuals)
                        CALL compute_sigma_gw(t_3c_overl_int_gw_AO, t_3c_overl_int_gw_RI, &
                                              mat_greens_fct_occ, mat_W, [1.0_dp, -1.0_dp], &
                                              Sigma_c_gw_tmp, [mo_start, mo_end], para_env, unit_nr_prv, &
                                              t_3c_ctr_RI=t_3c_ctr_RI, t_3c_ctr_in=.TRUE.)

                     END IF

                     IF (do_virt) THEN
                        ! 3-center * G(i(w-w')) for a virtual level n (other Fermi level than occupieds)
                        CALL compute_sigma_gw(t_3c_overl_int_gw_AO, t_3c_overl_int_gw_RI, &
                                              mat_greens_fct_virt, mat_W, [1.0_dp, 1.0_dp], &
                                              Sigma_c_gw_tmp, [mo_start, mo_end], para_env, unit_nr_prv, &
                                              t_3c_ctr_RI=t_3c_ctr_RI, t_3c_ctr_in=.TRUE.)
                     END IF

                     CALL dbcsr_t_destroy(t_3c_ctr_RI)

                     ! divide the integration weight by 2, because the integration
                     ! is from -infty to +infty and not just 0 to +infty as for RPA
                     IF (jquad .LE. num_integ_points_freq) THEN
                        weight = wj(jquad)/2.0_dp
                     ELSE
                        weight = wj(jquad - num_integ_points_freq)/2.0_dp
                     END IF

                     IF (do_occ) THEN
                        vec_Sigma_c_gw(1:gw_corr_lev_occ, iquad, 1) = &
                           vec_Sigma_c_gw(1:gw_corr_lev_occ, iquad, 1) - &
                           weight/twopi*re_im*Sigma_c_gw_tmp(:)
                     ELSE IF (do_virt) THEN
                        vec_Sigma_c_gw(gw_corr_lev_occ + 1:gw_corr_lev_tot, iquad, 1) = &
                           vec_Sigma_c_gw(gw_corr_lev_occ + 1:gw_corr_lev_tot, iquad, 1) - &
                           weight/twopi*re_im*Sigma_c_gw_tmp(:)
                     END IF

                     DEALLOCATE (Sigma_c_gw_tmp)

                  END DO

               END DO

            END DO

         END DO

         ! empirically found that the sign of real part of self-energy for occupied levels
         ! is wrong when comparing to N^4 code...
         n = gw_corr_lev_occ
         vec_Sigma_c_gw(1:n, :, :) = -REAL(vec_Sigma_c_gw(1:n, :, :))*(1.0_dp, 0.0_dp) &
                                     + AIMAG(vec_Sigma_c_gw(1:n, :, :))*(0.0_dp, 1.0_dp)

      ELSE

         DO jquad = 1, num_integ_points

            tau = tau_tj(jquad)

            CALL compute_Greens_function_time(mat_greens_fct_occ, mat_greens_fct_virt, &
                                              fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                              fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                              fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, Eigenval, &
                                              jquad, nmo, eps_filter, e_fermi, stabilize_exp, tau_tj)

            CALL dbcsr_set(mat_W, 0.0_dp)
            CALL copy_fm_to_dbcsr(fm_mat_W(jquad)%matrix, mat_W, keep_sparsity=.FALSE.)

            CALL compute_sigma_gw(t_3c_overl_int_gw_AO, t_3c_overl_int_gw_RI, &
                                  mat_greens_fct_occ, mat_W, [1.0_dp, -1.0_dp], &
                                  vec_Sigma_c_gw_neg_tau(:, jquad), [mo_start, mo_end], para_env, unit_nr_prv, &
                                  t_3c_ctr_RI=t_3c_ctr_RI, t_3c_ctr_in=.FALSE.)

            CALL compute_sigma_gw(t_3c_overl_int_gw_AO, t_3c_overl_int_gw_RI, &
                                  mat_greens_fct_virt, mat_W, [1.0_dp, 1.0_dp], &
                                  vec_Sigma_c_gw_pos_tau(:, jquad), [mo_start, mo_end], para_env, unit_nr_prv, &
                                  t_3c_ctr_RI=t_3c_ctr_RI, t_3c_ctr_in=.TRUE.)

            CALL dbcsr_t_destroy(t_3c_ctr_RI)

            vec_Sigma_c_gw_cos_tau(:, jquad) = 0.5_dp*(vec_Sigma_c_gw_pos_tau(:, jquad) + &
                                                       vec_Sigma_c_gw_neg_tau(:, jquad))

            vec_Sigma_c_gw_sin_tau(:, jquad) = 0.5_dp*(vec_Sigma_c_gw_pos_tau(:, jquad) - &
                                                       vec_Sigma_c_gw_neg_tau(:, jquad))

         END DO ! jquad (tau)

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
                                                  im_unit*vec_Sigma_c_gw_sin_omega(:, 1:num_fit_points)

      END IF

      CALL dbcsr_release_P(mat_greens_fct_occ)
      CALL dbcsr_release_P(mat_greens_fct_virt)

      IF (do_ri_Sigma_x .AND. count_ev_sc_GW == 1) THEN

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

         CALL compute_sigma_gw(t_3c_overl_int_gw_AO, t_3c_overl_int_gw_RI, &
                               mat_dm%matrix, mat_SinvVSinv%matrix, [1.0_dp, -1.0_dp], &
                               vec_Sigma_x_gw(mo_start:mo_end, 1), [mo_start, mo_end], para_env, unit_nr_prv)

         IF (my_do_beta) THEN

            mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2, 1) = &
               mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2, 1) + &
               vec_Sigma_x_gw(:, 1)

         ELSE

            mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1, 1) = &
               mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1, 1) + &
               vec_Sigma_x_gw(:, 1)

         END IF

      END IF

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
                     (1.0_dp/(im_unit*(omega_i + omega_sign) + e_fermi - Eigenval(n_level_gw_ref)) + &
                      1.0_dp/(im_unit*(-omega_i + omega_sign) + e_fermi - Eigenval(n_level_gw_ref)))

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
!> \param mat_AO ...
!> \param mat_RI ...
!> \param prefac ...
!> \param vec_Sigma ...
!> \param mo_bounds ...
!> \param para_env ...
!> \param unit_nr ...
!> \param t_3c_ctr_RI ...
!> \param t_3c_ctr_in ...
! **************************************************************************************************
   SUBROUTINE compute_sigma_gw(t_3c_overl_int_gw_AO, t_3c_overl_int_gw_RI, &
                               mat_AO, mat_RI, prefac, &
                               vec_Sigma, mo_bounds, para_env, unit_nr, &
                               t_3c_ctr_RI, t_3c_ctr_in)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t_3c_overl_int_gw_AO, &
                                                            t_3c_overl_int_gw_RI
      TYPE(dbcsr_type), INTENT(IN)                       :: mat_AO, mat_RI
      REAL(dp), DIMENSION(2), INTENT(IN)                 :: prefac
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: vec_Sigma
      INTEGER, DIMENSION(2), INTENT(IN)                  :: mo_bounds
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(dbcsr_t_type), INTENT(INOUT), OPTIONAL        :: t_3c_ctr_RI
      LOGICAL, INTENT(IN), OPTIONAL                      :: t_3c_ctr_in

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_sigma_gw', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: t_3c_ctr_in_prv
      TYPE(dbcsr_t_type)                                 :: t_3c_ctr_RI_prv, t_3c_overl_int_gw_copy, &
                                                            t_AO, t_RI

      CALL timeset(routineN, handle)

      IF (PRESENT(t_3c_ctr_in)) THEN
         t_3c_ctr_in_prv = t_3c_ctr_in
      ELSE
         t_3c_ctr_in_prv = .FALSE.
      ENDIF

      IF (t_3c_ctr_in_prv) THEN
         CPASSERT(PRESENT(t_3c_ctr_RI))
         t_3c_ctr_RI_prv = t_3c_ctr_RI
      ELSE

         CALL dbcsr_t_create(mat_RI, t_RI, name="(RI|RI)")
         CALL dbcsr_t_copy_matrix_to_tensor(mat_RI, t_RI)

         CALL dbcsr_t_create(t_3c_overl_int_gw_RI, t_3c_overl_int_gw_copy)

         CALL dbcsr_t_contract(dbcsr_scalar(prefac(1)), t_RI, t_3c_overl_int_gw_RI, dbcsr_scalar(0.0_dp), &
                               t_3c_overl_int_gw_copy, &
                               contract_1=[2], notcontract_1=[1], &
                               contract_2=[1], notcontract_2=[2, 3], &
                               map_1=[1], map_2=[2, 3], &
                               unit_nr=unit_nr)

         CALL dbcsr_t_create(t_3c_overl_int_gw_AO, t_3c_ctr_RI_prv)
         CALL dbcsr_t_copy(t_3c_overl_int_gw_copy, t_3c_ctr_RI_prv, order=[2, 1, 3])
         CALL dbcsr_t_destroy(t_3c_overl_int_gw_copy)
         CALL dbcsr_t_destroy(t_RI)

      ENDIF

      CALL dbcsr_t_create(mat_AO, t_AO, name="(AO|AO)")
      CALL dbcsr_t_copy_matrix_to_tensor(mat_AO, t_AO)

      CALL dbcsr_t_create(t_3c_overl_int_gw_AO, t_3c_overl_int_gw_copy)
      CALL dbcsr_t_contract(dbcsr_scalar(prefac(2)), t_AO, t_3c_overl_int_gw_AO, dbcsr_scalar(0.0_dp), &
                            t_3c_overl_int_gw_copy, &
                            contract_1=[2], notcontract_1=[1], &
                            contract_2=[1], notcontract_2=[2, 3], &
                            map_1=[1], map_2=[2, 3], &
                            unit_nr=unit_nr)

      CALL trace_sigma_gw(t_3c_ctr_RI_prv, t_3c_overl_int_gw_copy, vec_sigma, mo_bounds, para_env)
      CALL dbcsr_t_destroy(t_3c_overl_int_gw_copy)

      IF (PRESENT(t_3c_ctr_in)) THEN
         IF (.NOT. t_3c_ctr_in) THEN
            CPASSERT(PRESENT(t_3c_ctr_RI))
            t_3c_ctr_RI = t_3c_ctr_RI_prv
         ENDIF
      ELSE
         CALL dbcsr_t_destroy(t_3c_ctr_RI_prv)
      ENDIF

      CALL dbcsr_t_destroy(t_AO)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param t3c_1 ...
!> \param t3c_2 ...
!> \param vec_sigma ...
!> \param mo_bounds ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE trace_sigma_gw(t3c_1, t3c_2, vec_sigma, mo_bounds, para_env)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t3c_1, t3c_2
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: vec_Sigma
      INTEGER, DIMENSION(2), INTENT(IN)                  :: mo_bounds
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'trace_sigma_gw', routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, handle, n, n_end, n_end_block, &
                                                            n_start, n_start_block
      INTEGER, DIMENSION(1)                              :: trace_shape
      INTEGER, DIMENSION(3)                              :: boff, bsize, ind
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: block_1, block_2
      TYPE(dbcsr_t_iterator_type)                        :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_t_iterator_start(iter, t3c_1)
      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, ind, blk, blk_size=bsize, blk_offset=boff)
         CALL dbcsr_t_get_block(t3c_1, ind, block_1, found)
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
         vec_Sigma(n_start:n_end) = &
            vec_Sigma(n_start:n_end) + &
            (/(DOT_PRODUCT(RESHAPE(block_1(:, :, n), trace_shape), &
                           RESHAPE(block_2(:, :, n), trace_shape)), &
               n=n_start_block, n_end_block)/)
         DEALLOCATE (block_1, block_2)
      ENDDO
      CALL dbcsr_t_iterator_stop(iter)

      CALL mp_sum(vec_Sigma, para_env%group)

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
!> \param jquad ...
!> \param nmo ...
!> \param eps_filter ...
!> \param e_fermi ...
!> \param stabilize_exp ...
!> \param tau_tj ...
! **************************************************************************************************
   SUBROUTINE compute_Greens_function_time(mat_greens_fct_occ, mat_greens_fct_virt, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                           fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                           fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, Eigenval, jquad, nmo, &
                                           eps_filter, e_fermi, stabilize_exp, tau_tj)

      TYPE(dbcsr_type), POINTER                          :: mat_greens_fct_occ, mat_greens_fct_virt
      TYPE(cp_fm_type), POINTER :: fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
         fm_mo_coeff_virt_scaled, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: jquad, nmo
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter, e_fermi, stabilize_exp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Greens_function_time', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, jjB, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: tau

      CALL timeset(routineN, handle)

      tau = tau_tj(jquad)

      ! get info of fm_mo_coeff_occ
      CALL cp_fm_get_info(matrix=fm_mo_coeff_occ, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! Multiply the occupied and the virtual MO coefficients with the factor exp((-e_i-e_F)*tau/2).
      ! Then, we simply get the sum over all occ states and virt. states by a simple matrix-matrix
      ! multiplication.

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

! **************************************************************************************************
!> \brief ...
!> \param mat_greens_fct_occ ...
!> \param mat_greens_fct_virt ...
!> \param fm_mo_coeff ...
!> \param fm_mo_coeff_occ_scaled ...
!> \param fm_mo_coeff_virt_scaled ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param Eigenval ...
!> \param fermi_level_offset ...
!> \param nmo ...
!> \param homo ...
!> \param eps_filter ...
!> \param omega ...
!> \param real_imaginary ...
! **************************************************************************************************
   SUBROUTINE compute_Greens_function_freq(mat_greens_fct_occ, mat_greens_fct_virt, fm_mo_coeff, &
                                           fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, fm_scaled_dm_occ_tau, &
                                           fm_scaled_dm_virt_tau, Eigenval, fermi_level_offset, nmo, &
                                           homo, eps_filter, omega, real_imaginary)

      TYPE(dbcsr_type), POINTER                          :: mat_greens_fct_occ, mat_greens_fct_virt
      TYPE(cp_fm_type), POINTER :: fm_mo_coeff, fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
         fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: fermi_level_offset
      INTEGER, INTENT(IN)                                :: nmo, homo
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter, omega
      INTEGER, INTENT(IN)                                :: real_imaginary

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Greens_function_freq', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ncol_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices
      REAL(KIND=dp)                                      :: e_fermi

      CALL timeset(routineN, handle)

      ! get info of fm_mo_coeff
      CALL cp_fm_get_info(matrix=fm_mo_coeff, &
                          ncol_local=ncol_local, &
                          col_indices=col_indices)

      ! we have two different fermi levels for correcting occupied and virtual energy levels
      e_fermi = Eigenval(homo) + fermi_level_offset
      CALL scale_mo_coeff(fm_mo_coeff_occ_scaled, fm_mo_coeff, e_fermi, ncol_local, &
                          Eigenval, omega, real_imaginary, col_indices)

      e_fermi = Eigenval(homo + 1) - fermi_level_offset
      CALL scale_mo_coeff(fm_mo_coeff_virt_scaled, fm_mo_coeff, e_fermi, ncol_local, &
                          Eigenval, omega, real_imaginary, col_indices)

      CALL cp_gemm(transa="N", transb="T", m=nmo, n=nmo, k=nmo, alpha=1.0_dp, &
                   matrix_a=fm_mo_coeff, matrix_b=fm_mo_coeff_occ_scaled, beta=0.0_dp, &
                   matrix_c=fm_scaled_dm_occ_tau)

      CALL cp_gemm(transa="N", transb="T", m=nmo, n=nmo, k=nmo, alpha=1.0_dp, &
                   matrix_a=fm_mo_coeff, matrix_b=fm_mo_coeff_virt_scaled, beta=0.0_dp, &
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

   END SUBROUTINE compute_Greens_function_freq

! **************************************************************************************************
!> \brief ...
!> \param fm_mo_coeff_scaled ...
!> \param fm_mo_coeff ...
!> \param e_fermi ...
!> \param ncol_local ...
!> \param Eigenval ...
!> \param omega ...
!> \param real_imaginary ...
!> \param col_indices ...
! **************************************************************************************************
   SUBROUTINE scale_mo_coeff(fm_mo_coeff_scaled, fm_mo_coeff, e_fermi, ncol_local, &
                             Eigenval, omega, real_imaginary, col_indices)

      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff_scaled, fm_mo_coeff
      REAL(KIND=dp)                                      :: e_fermi
      INTEGER, INTENT(IN)                                :: ncol_local
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: omega
      INTEGER, INTENT(IN)                                :: real_imaginary
      INTEGER, DIMENSION(:), POINTER                     :: col_indices

      CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_mo_coeff', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB
      REAL(KIND=dp)                                      :: energy_diff

      CALL timeset(routineN, handle)

      DO iiB = 1, ncol_local

         i_global = col_indices(iiB)

         energy_diff = e_fermi - Eigenval(i_global)

         IF (real_imaginary == 1) THEN
            fm_mo_coeff_scaled%local_data(:, iiB) = fm_mo_coeff%local_data(:, iiB)* &
                                                    energy_diff/(omega**2 + energy_diff**2)
         ELSE IF (real_imaginary == 2) THEN
            fm_mo_coeff_scaled%local_data(:, iiB) = fm_mo_coeff%local_data(:, iiB)* &
                                                    (-omega/(omega**2 + energy_diff**2))
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE scale_mo_coeff

END MODULE rpa_gw

