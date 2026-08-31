!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief module that contains the definitions of the scf types
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE qs_scf_types
   USE cp_dbcsr_api,                    ONLY: dbcsr_deallocate_matrix,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE outer_scf_control_types,         ONLY: qs_outer_scf_type
   USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                              preconditioner_p_type
   USE qs_block_davidson_types,         ONLY: block_davidson_release,&
                                              davidson_type
   USE qs_density_mixing_types,         ONLY: direct_mixing_nr,&
                                              mixing_storage_create,&
                                              mixing_storage_release,&
                                              mixing_storage_type
   USE qs_diis_types,                   ONLY: qs_diis_b_release,&
                                              qs_diis_buffer_type
   USE qs_fb_env_types,                 ONLY: fb_env_has_data,&
                                              fb_env_nullify,&
                                              fb_env_obj,&
                                              fb_env_release
   USE qs_ot_types,                     ONLY: qs_ot_destroy,&
                                              qs_ot_restart_signature_nchar,&
                                              qs_ot_restart_signature_nint,&
                                              qs_ot_restart_signature_nlog,&
                                              qs_ot_restart_signature_nreal,&
                                              qs_ot_type
   USE qs_scf_subspace_types,           ONLY: qs_scf_subspace_buffer_release,&
                                              qs_scf_subspace_buffer_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_types'

   INTEGER, PARAMETER, PUBLIC :: general_diag_method_nr = 1, &
                                 special_diag_method_nr = 2, &
                                 ot_diag_method_nr = 3, &
                                 block_krylov_diag_method_nr = 4, &
                                 block_davidson_diag_method_nr = 5, &
                                 ot_method_nr = 10, &
                                 filter_matrix_diag_method_nr = 20, &
                                 smeagol_method_nr = 41

   PUBLIC :: qs_scf_env_type
   PUBLIC :: qs_ot_restart_dense_matrix_type
   PUBLIC :: scf_env_release, scf_env_create, scf_env_did_change
   PUBLIC :: qs_scf_ot_restart_clear
   PUBLIC :: krylov_space_type, krylov_space_create, subspace_env_type
   PUBLIC :: diag_subspace_env_create

! **************************************************************************************************
!> \brief wrapper for temporary and cached objects used in the scf iteration
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   TYPE krylov_space_type
      INTEGER :: nkrylov = -1, nblock = -1, nmo_conv = -1, nmo_nc = -1, max_iter = -1
      LOGICAL :: always_check_conv = .FALSE.
      REAL(dp) :: eps_std_diag = -1.0_dp, eps_conv = -1.0_dp, eps_adapt = -1.0_dp, max_res_norm = -1.0_dp, min_res_norm = -1.0_dp
      REAL(dp), DIMENSION(:), POINTER :: c_eval => NULL(), t_eval => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER :: v_mat => NULL(), mo_conv => NULL(), mo_refine => NULL()
      TYPE(cp_fm_type), POINTER ::  tmp_mat => NULL()
      !NOTE: the following matrices are small and could be used as standard array rather than distributed fm
      TYPE(cp_fm_type), POINTER :: block1_mat => NULL(), block2_mat => NULL(), block3_mat => NULL(), &
                                   block4_mat => NULL(), block5_mat => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER ::  c_vec => NULL(), chc_mat => NULL()
   END TYPE krylov_space_type

   TYPE subspace_env_type
      INTEGER :: max_iter = -1, mixing_method = -1
      REAL(dp) :: eps_diag_sub = -1.0_dp, eps_ene = -1.0_dp, eps_adapt = -1.0_dp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: p_matrix_store => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER :: p_matrix_mix => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER :: chc_mat => NULL(), c_vec => NULL(), c0 => NULL()
      TYPE(mixing_storage_type), POINTER :: mixing_store => NULL()
   END TYPE subspace_env_type

   TYPE floating_basis_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: gradient => NULL()
   END TYPE floating_basis_type

   TYPE qs_ot_restart_dense_matrix_type
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: data
   END TYPE qs_ot_restart_dense_matrix_type

   TYPE qs_ot_restart_channel_type
      LOGICAL :: valid = .FALSE.
      INTEGER :: spin_index = 0, kpoint_index = 0, local_kpoint_index = 0
      INTEGER :: nmo = 0, nmo_ortho = 0
      INTEGER :: diis_m = 0, history_m = 0, lbfgs_m = 0
      LOGICAL :: has_complex_kpoint_state = .FALSE., do_rotation = .FALSE., do_ener = .FALSE.
      LOGICAL :: has_matrix_x = .FALSE., has_matrix_x_im = .FALSE.
      LOGICAL :: has_rot_mat_x = .FALSE., has_rot_mat_x_im = .FALSE., has_ener_x = .FALSE.
      LOGICAL :: has_matrix_dx = .FALSE., has_matrix_dx_im = .FALSE.
      LOGICAL :: has_matrix_gx_old = .FALSE., has_matrix_gx_old_im = .FALSE.
      LOGICAL :: has_rot_mat_dx = .FALSE., has_rot_mat_dx_im = .FALSE.
      LOGICAL :: has_rot_mat_gx_old = .FALSE., has_rot_mat_gx_old_im = .FALSE.
      LOGICAL :: has_ener_dx = .FALSE., has_ener_gx_old = .FALSE.
      LOGICAL :: has_diis_scalars = .FALSE., has_matrix_history = .FALSE.
      LOGICAL :: has_matrix_history_im = .FALSE., has_rot_history = .FALSE.
      LOGICAL :: has_rot_history_im = .FALSE., has_ener_history = .FALSE.
      LOGICAL :: has_lbfgs_scalars = .FALSE.
      CHARACTER(LEN=20), DIMENSION(qs_ot_restart_signature_nchar) :: settings_char = ""
      INTEGER, DIMENSION(qs_ot_restart_signature_nint) :: settings_int = 0
      LOGICAL, DIMENSION(qs_ot_restart_signature_nlog) :: settings_log = .FALSE.
      REAL(KIND=dp), DIMENSION(qs_ot_restart_signature_nreal) :: settings_real = 0.0_dp
      REAL(KIND=dp) :: kpoint_weight = 0.0_dp
      REAL(KIND=dp) :: ot_pos(53) = 0.0_dp, ot_energy(53) = 0.0_dp, ot_grad(53) = 0.0_dp
      INTEGER :: line_search_left = -1, line_search_right = -1, line_search_mid = -1
      INTEGER :: line_search_count = -1, diis_iter = -1
      LOGICAL :: line_search_might_be_done = .FALSE., energy_only = .FALSE.
      REAL(KIND=dp) :: delta = 0.0_dp, gnorm = 0.0_dp, gnorm_old = 0.0_dp
      REAL(KIND=dp) :: etotal = 0.0_dp, gradient = 0.0_dp, ds_min = 0.0_dp
      REAL(KIND=dp) :: broyden_adaptive_sigma = 0.0_dp
      TYPE(qs_ot_restart_dense_matrix_type) :: matrix_x_dense, matrix_x_im_dense
      TYPE(qs_ot_restart_dense_matrix_type) :: rot_mat_x_dense, rot_mat_x_im_dense
      TYPE(qs_ot_restart_dense_matrix_type) :: matrix_dx_dense, matrix_dx_im_dense
      TYPE(qs_ot_restart_dense_matrix_type) :: matrix_gx_old_dense, matrix_gx_old_im_dense
      TYPE(qs_ot_restart_dense_matrix_type) :: rot_mat_dx_dense, rot_mat_dx_im_dense
      TYPE(qs_ot_restart_dense_matrix_type) :: rot_mat_gx_old_dense, rot_mat_gx_old_im_dense
      TYPE(qs_ot_restart_dense_matrix_type), DIMENSION(:), ALLOCATABLE :: matrix_h_e_dense, &
                                                                          matrix_h_x_dense
      TYPE(qs_ot_restart_dense_matrix_type), DIMENSION(:), ALLOCATABLE :: matrix_h_e_im_dense, &
                                                                          matrix_h_x_im_dense
      TYPE(qs_ot_restart_dense_matrix_type), DIMENSION(:), ALLOCATABLE :: rot_mat_h_e_dense, &
                                                                          rot_mat_h_x_dense
      TYPE(qs_ot_restart_dense_matrix_type), DIMENSION(:), ALLOCATABLE :: rot_mat_h_e_im_dense, &
                                                                          rot_mat_h_x_im_dense
      REAL(KIND=dp), DIMENSION(:), POINTER :: ener_x => NULL(), ener_dx => NULL(), ener_gx_old => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: c_diis => NULL(), c_broy => NULL(), energy_h => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: lbfgs_rho => NULL(), lbfgs_yy => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: lbfgs_sy_rotation => NULL(), lbfgs_yy_rotation => NULL()
      INTEGER, DIMENSION(:), POINTER :: ipivot => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: ls_diis => NULL(), lss_diis => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: ener_h_e => NULL(), ener_h_x => NULL()
   END TYPE qs_ot_restart_channel_type

   TYPE qs_ot_restart_cache_type
      LOGICAL :: valid = .FALSE.
      INTEGER :: format_version = 0, nchannel = 0
      TYPE(qs_ot_restart_channel_type), DIMENSION(:), POINTER :: channel => NULL()
   END TYPE qs_ot_restart_cache_type
! **************************************************************************************************
   TYPE qs_scf_env_type
      TYPE(qs_outer_scf_type) :: outer_scf = qs_outer_scf_type()
      INTEGER :: iter_count = -1
      INTEGER :: cholesky_method = -1, nelectron = -1, method = -1, mixing_method = -1, nskip_mixing = -1
      REAL(KIND=dp) :: iter_param = -1.0_dp, iter_delta = -1.0_dp, p_mix_alpha = -1.0_dp, sum_zeff_corr = -1.0_dp
      ! ADIIS diagnostics: accepted-step norm and raw-map delta.
      REAL(KIND=dp) :: step_norm = -1.0_dp, raw_map_delta = -1.0_dp
      ! ODA diagnostics for the most recent raw-step line search.
      REAL(KIND=dp) :: oda_lambda = 0.0_dp, oda_energy = 0.0_dp, oda_gradient0 = 0.0_dp, oda_gradient1 = 0.0_dp
      INTEGER :: oda_evaluations = 0
      CHARACTER(len=15) :: iter_method = "", oda_status = ""
      COMPLEX(KIND=dp), DIMENSION(:, :, :), POINTER :: cc_buffer => NULL()
      LOGICAL :: print_iter_line = .FALSE., skip_mixing = .FALSE., skip_diis = .FALSE., needs_ortho = .FALSE., &
                 adiis_check_next = .FALSE., adiis_validated = .FALSE., &
                 raw_map_delta_valid = .FALSE.
      TYPE(mixing_storage_type), POINTER :: mixing_store => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER :: scf_work1 => NULL(), scf_work1_red => NULL()
      TYPE(cp_fm_type), POINTER  :: scf_work2 => NULL(), ortho => NULL(), ortho_m1 => NULL(), &
                                    s_half => NULL(), s_minus_one => NULL(), &
                                    scf_work2_red => NULL(), ortho_red => NULL(), ortho_m1_red => NULL()
      TYPE(krylov_space_type), POINTER :: krylov_space => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER :: p_delta => NULL(), p_mix_new => NULL()
      TYPE(dbcsr_type), POINTER :: ortho_dbcsr => NULL(), buf1_dbcsr => NULL(), buf2_dbcsr => NULL()
      TYPE(preconditioner_p_type), DIMENSION(:), POINTER :: ot_preconditioner => NULL()
      TYPE(qs_ot_type), POINTER, DIMENSION(:)  :: qs_ot_env => NULL()
      TYPE(qs_ot_restart_cache_type) :: qs_ot_restart = qs_ot_restart_cache_type()
      TYPE(qs_diis_buffer_type), POINTER :: scf_diis_buffer => NULL()
      TYPE(qs_scf_subspace_buffer_type), POINTER :: scf_subspace_buffer => NULL()
      TYPE(subspace_env_type), POINTER :: subspace_env => NULL()
      TYPE(davidson_type), POINTER, DIMENSION(:) :: block_davidson_env => NULL()
      TYPE(fb_env_obj) :: filter_matrix_env = fb_env_obj()
      TYPE(floating_basis_type) :: floating_basis = floating_basis_type()
      !> reference molecular orbitals for the maximum overlap method
      TYPE(cp_fm_type), DIMENSION(:), POINTER          :: mom_ref_mo_coeff => NULL()
      !> MOM-related work matrices
      TYPE(cp_fm_type), DIMENSION(:), POINTER          :: mom_overlap => NULL(), mom_s_mo_coeff => NULL()
   END TYPE qs_scf_env_type

CONTAINS

! **************************************************************************************************
!> \brief allocates and initialize an scf_env
!> \param scf_env the scf env to initialize
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE scf_env_create(scf_env)
      TYPE(qs_scf_env_type), INTENT(OUT)                 :: scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'scf_env_create'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      scf_env%iter_count = 0
      scf_env%nelectron = 0
      scf_env%iter_param = 0.0_dp
      scf_env%iter_delta = 0.0_dp
      scf_env%step_norm = 0.0_dp
      scf_env%raw_map_delta = 0.0_dp
      scf_env%oda_lambda = 0.0_dp
      scf_env%oda_energy = 0.0_dp
      scf_env%oda_gradient0 = 0.0_dp
      scf_env%oda_gradient1 = 0.0_dp
      scf_env%oda_evaluations = 0
      scf_env%oda_status = ""
      scf_env%iter_method = ""
      scf_env%print_iter_line = .TRUE.
      scf_env%skip_mixing = .FALSE.
      scf_env%skip_diis = .FALSE.
      scf_env%needs_ortho = .FALSE.
      scf_env%adiis_check_next = .FALSE.
      scf_env%adiis_validated = .FALSE.
      scf_env%raw_map_delta_valid = .FALSE.
      scf_env%method = -1
      scf_env%mixing_method = -1
      scf_env%p_mix_alpha = 1.0_dp
      scf_env%cholesky_method = -1
      scf_env%outer_scf%iter_count = 0
      scf_env%sum_zeff_corr = 0.0_dp
      NULLIFY (scf_env%outer_scf%variables)
      NULLIFY (scf_env%outer_scf%gradient)
      NULLIFY (scf_env%outer_scf%energy)
      NULLIFY (scf_env%outer_scf%count)
      NULLIFY (scf_env%outer_scf%inv_jacobian)
      scf_env%outer_scf%deallocate_jacobian = .TRUE.
      NULLIFY (scf_env%scf_work1)
      NULLIFY (scf_env%scf_work2)
      NULLIFY (scf_env%scf_work1_red)
      NULLIFY (scf_env%scf_work2_red)
      NULLIFY (scf_env%ortho)
      NULLIFY (scf_env%ortho_red)
      NULLIFY (scf_env%ortho_dbcsr)
      NULLIFY (scf_env%ortho_m1)
      NULLIFY (scf_env%ortho_m1_red)
      NULLIFY (scf_env%p_mix_new)
      NULLIFY (scf_env%ot_preconditioner)
      NULLIFY (scf_env%qs_ot_env)
      scf_env%qs_ot_restart%valid = .FALSE.
      scf_env%qs_ot_restart%format_version = 0
      scf_env%qs_ot_restart%nchannel = 0
      NULLIFY (scf_env%qs_ot_restart%channel)
      NULLIFY (scf_env%scf_diis_buffer)
      NULLIFY (scf_env%scf_subspace_buffer)
      NULLIFY (scf_env%buf1_dbcsr)
      NULLIFY (scf_env%buf2_dbcsr)
      NULLIFY (scf_env%s_half)
      NULLIFY (scf_env%p_delta)
      NULLIFY (scf_env%s_minus_one)
      NULLIFY (scf_env%cc_buffer)
      NULLIFY (scf_env%mixing_store)
      NULLIFY (scf_env%krylov_space)
      NULLIFY (scf_env%subspace_env)
      NULLIFY (scf_env%block_davidson_env)
      NULLIFY (scf_env%floating_basis%gradient)
      CALL fb_env_nullify(scf_env%filter_matrix_env)
      NULLIFY (scf_env%mom_ref_mo_coeff)
      NULLIFY (scf_env%mom_overlap)
      NULLIFY (scf_env%mom_s_mo_coeff)

      CALL timestop(handle)

   END SUBROUTINE scf_env_create

! **************************************************************************************************
!> \brief function to be called to inform the scf_env about changes
!> \param scf_env the scf env to inform
!> \par History
!>      03.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE scf_env_did_change(scf_env)
      TYPE(qs_scf_env_type), INTENT(INOUT)               :: scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_did_change'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(scf_env%p_mix_new)) THEN
         CALL dbcsr_deallocate_matrix_set(scf_env%p_mix_new)
      END IF
      IF (ASSOCIATED(scf_env%p_delta)) THEN
         CALL dbcsr_deallocate_matrix_set(scf_env%p_delta)
      END IF
      IF (ASSOCIATED(scf_env%scf_subspace_buffer)) THEN
         CALL qs_scf_subspace_buffer_release(scf_env%scf_subspace_buffer)
         DEALLOCATE (scf_env%scf_subspace_buffer)
      END IF
      CALL cp_fm_release(scf_env%mom_ref_mo_coeff)

      CALL timestop(handle)

   END SUBROUTINE scf_env_did_change

! **************************************************************************************************
!> \brief releases an scf_env (see doc/ReferenceCounting.html)
!> \param scf_env the environment to release
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE scf_env_release(scf_env)

      TYPE(qs_scf_env_type), INTENT(INOUT)               :: scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'scf_env_release'

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      CALL cp_fm_release(scf_env%scf_work1)
      IF (ASSOCIATED(scf_env%scf_work1_red)) THEN
         CALL cp_fm_release(scf_env%scf_work1_red)
         DEALLOCATE (scf_env%scf_work1_red)
         NULLIFY (scf_env%scf_work1_red)
      END IF
      IF (ASSOCIATED(scf_env%scf_work2)) THEN
         CALL cp_fm_release(scf_env%scf_work2)
         DEALLOCATE (scf_env%scf_work2)
         NULLIFY (scf_env%scf_work2)
      END IF
      IF (ASSOCIATED(scf_env%scf_work2_red)) THEN
         CALL cp_fm_release(scf_env%scf_work2_red)
         DEALLOCATE (scf_env%scf_work2_red)
         NULLIFY (scf_env%scf_work2_red)
      END IF
      IF (ASSOCIATED(scf_env%ortho)) THEN
         CALL cp_fm_release(scf_env%ortho)
         DEALLOCATE (scf_env%ortho)
         NULLIFY (scf_env%ortho)
      END IF
      IF (ASSOCIATED(scf_env%ortho_red)) THEN
         CALL cp_fm_release(scf_env%ortho_red)
         DEALLOCATE (scf_env%ortho_red)
         NULLIFY (scf_env%ortho_red)
      END IF
      IF (ASSOCIATED(scf_env%ortho_m1)) THEN
         CALL cp_fm_release(scf_env%ortho_m1)
         DEALLOCATE (scf_env%ortho_m1)
         NULLIFY (scf_env%ortho_m1)
      END IF
      IF (ASSOCIATED(scf_env%ortho_m1_red)) THEN
         CALL cp_fm_release(scf_env%ortho_m1_red)
         DEALLOCATE (scf_env%ortho_m1_red)
         NULLIFY (scf_env%ortho_m1_red)
      END IF
      IF (ASSOCIATED(scf_env%ortho_dbcsr)) THEN
         ! we should not end up here, and give back using the pools
         CALL dbcsr_deallocate_matrix(scf_env%ortho_dbcsr)
      END IF
      IF (ASSOCIATED(scf_env%buf1_dbcsr)) THEN
         ! we should not end up here, and give back using the pools
         CALL dbcsr_deallocate_matrix(scf_env%buf1_dbcsr)
      END IF
      IF (ASSOCIATED(scf_env%buf2_dbcsr)) THEN
         ! we should not end up here, and give back using the pools
         CALL dbcsr_deallocate_matrix(scf_env%buf2_dbcsr)
      END IF
      IF (ASSOCIATED(scf_env%s_half)) THEN
         CALL cp_fm_release(scf_env%s_half)
         DEALLOCATE (scf_env%s_half)
      END IF
      IF (ASSOCIATED(scf_env%s_minus_one)) THEN
         CALL cp_fm_release(scf_env%s_minus_one)
         DEALLOCATE (scf_env%s_minus_one)
      END IF
      IF (ASSOCIATED(scf_env%p_mix_new)) THEN
         ! we should not end up here, and give back using the pools
         CALL dbcsr_deallocate_matrix_set(scf_env%p_mix_new)
      END IF
      IF (ASSOCIATED(scf_env%p_delta)) THEN
         ! we should not end up here, and give back using the pools
         CALL dbcsr_deallocate_matrix_set(scf_env%p_delta)
      END IF
      IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
         DO i = 1, SIZE(scf_env%ot_preconditioner)
            IF (ASSOCIATED(scf_env%ot_preconditioner(i)%preconditioner)) THEN
               CALL destroy_preconditioner(scf_env%ot_preconditioner(i)%preconditioner)
               DEALLOCATE (scf_env%ot_preconditioner(i)%preconditioner)
            END IF
         END DO
         DEALLOCATE (scf_env%ot_preconditioner)
      END IF
      IF (ASSOCIATED(scf_env%qs_ot_env)) THEN
         DO i = 1, SIZE(scf_env%qs_ot_env)
            CALL qs_ot_destroy(scf_env%qs_ot_env(i))
         END DO
         DEALLOCATE (scf_env%qs_ot_env)
      END IF
      CALL qs_scf_ot_restart_clear(scf_env)
      IF (ASSOCIATED(scf_env%scf_diis_buffer)) THEN
         CALL qs_diis_b_release(scf_env%scf_diis_buffer)
         DEALLOCATE (scf_env%scf_diis_buffer)
      END IF
      IF (ASSOCIATED(scf_env%scf_subspace_buffer)) THEN
         CALL qs_scf_subspace_buffer_release(scf_env%scf_subspace_buffer)
         DEALLOCATE (scf_env%scf_subspace_buffer)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%variables)) THEN
         DEALLOCATE (scf_env%outer_scf%variables)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%count)) THEN
         DEALLOCATE (scf_env%outer_scf%count)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%gradient)) THEN
         DEALLOCATE (scf_env%outer_scf%gradient)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
         DEALLOCATE (scf_env%outer_scf%inv_jacobian)
      END IF
      IF (ASSOCIATED(scf_env%outer_scf%energy)) THEN
         DEALLOCATE (scf_env%outer_scf%energy)
      END IF
      IF (ASSOCIATED(scf_env%cc_buffer)) THEN
         DEALLOCATE (scf_env%cc_buffer)
      END IF
      IF (ASSOCIATED(scf_env%mixing_store)) THEN
         CALL mixing_storage_release(scf_env%mixing_store)
         DEALLOCATE (scf_env%mixing_store)
      END IF
      IF (ASSOCIATED(scf_env%krylov_space)) THEN
         CALL krylov_space_release(scf_env%krylov_space)
      END IF
      IF (ASSOCIATED(scf_env%subspace_env)) THEN
         CALL diag_subspace_env_release(scf_env%subspace_env)
      END IF
      IF (ASSOCIATED(scf_env%block_davidson_env)) THEN
         CALL block_davidson_release(scf_env%block_davidson_env)
      END IF
      IF (fb_env_has_data(scf_env%filter_matrix_env)) THEN
         CALL fb_env_release(scf_env%filter_matrix_env)
      END IF
      IF (ASSOCIATED(scf_env%floating_basis%gradient)) THEN
         DEALLOCATE (scf_env%floating_basis%gradient)
      END IF
      CALL cp_fm_release(scf_env%mom_ref_mo_coeff)
      CALL cp_fm_release(scf_env%mom_overlap)
      CALL cp_fm_release(scf_env%mom_s_mo_coeff)

      CALL timestop(handle)

   END SUBROUTINE scf_env_release

! **************************************************************************************************
!> \brief clear cached k-point OT restart data
!> \param scf_env the environment containing the restart cache
! **************************************************************************************************
   SUBROUTINE qs_scf_ot_restart_clear(scf_env)
      TYPE(qs_scf_env_type), INTENT(INOUT)               :: scf_env

      INTEGER                                            :: ichannel, ihistory

      IF (ASSOCIATED(scf_env%qs_ot_restart%channel)) THEN
         DO ichannel = 1, SIZE(scf_env%qs_ot_restart%channel)
            CALL qs_ot_restart_dense_matrix_clear(scf_env%qs_ot_restart%channel(ichannel)%matrix_x_dense)
            CALL qs_ot_restart_dense_matrix_clear(scf_env%qs_ot_restart%channel(ichannel)%matrix_x_im_dense)
            CALL qs_ot_restart_dense_matrix_clear(scf_env%qs_ot_restart%channel(ichannel)%rot_mat_x_dense)
            CALL qs_ot_restart_dense_matrix_clear(scf_env%qs_ot_restart%channel(ichannel)%rot_mat_x_im_dense)
            CALL qs_ot_restart_dense_matrix_clear(scf_env%qs_ot_restart%channel(ichannel)%matrix_dx_dense)
            CALL qs_ot_restart_dense_matrix_clear(scf_env%qs_ot_restart%channel(ichannel)%matrix_dx_im_dense)
            CALL qs_ot_restart_dense_matrix_clear(scf_env%qs_ot_restart%channel(ichannel)%matrix_gx_old_dense)
            CALL qs_ot_restart_dense_matrix_clear(scf_env%qs_ot_restart%channel(ichannel)%matrix_gx_old_im_dense)
            CALL qs_ot_restart_dense_matrix_clear(scf_env%qs_ot_restart%channel(ichannel)%rot_mat_dx_dense)
            CALL qs_ot_restart_dense_matrix_clear(scf_env%qs_ot_restart%channel(ichannel)%rot_mat_dx_im_dense)
            CALL qs_ot_restart_dense_matrix_clear(scf_env%qs_ot_restart%channel(ichannel)%rot_mat_gx_old_dense)
            CALL qs_ot_restart_dense_matrix_clear(scf_env%qs_ot_restart%channel(ichannel)%rot_mat_gx_old_im_dense)
            IF (ALLOCATED(scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_dense)) THEN
               DO ihistory = 1, SIZE(scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_dense)
                  CALL qs_ot_restart_dense_matrix_clear( &
                     scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_dense(ihistory))
                  CALL qs_ot_restart_dense_matrix_clear( &
                     scf_env%qs_ot_restart%channel(ichannel)%matrix_h_x_dense(ihistory))
               END DO
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_dense, &
                           scf_env%qs_ot_restart%channel(ichannel)%matrix_h_x_dense)
            END IF
            IF (ALLOCATED(scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_im_dense)) THEN
               DO ihistory = 1, SIZE(scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_im_dense)
                  CALL qs_ot_restart_dense_matrix_clear( &
                     scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_im_dense(ihistory))
                  CALL qs_ot_restart_dense_matrix_clear( &
                     scf_env%qs_ot_restart%channel(ichannel)%matrix_h_x_im_dense(ihistory))
               END DO
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%matrix_h_e_im_dense, &
                           scf_env%qs_ot_restart%channel(ichannel)%matrix_h_x_im_dense)
            END IF
            IF (ALLOCATED(scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_dense)) THEN
               DO ihistory = 1, SIZE(scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_dense)
                  CALL qs_ot_restart_dense_matrix_clear( &
                     scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_dense(ihistory))
                  CALL qs_ot_restart_dense_matrix_clear( &
                     scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_x_dense(ihistory))
               END DO
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_dense, &
                           scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_x_dense)
            END IF
            IF (ALLOCATED(scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_im_dense)) THEN
               DO ihistory = 1, SIZE(scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_im_dense)
                  CALL qs_ot_restart_dense_matrix_clear( &
                     scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_im_dense(ihistory))
                  CALL qs_ot_restart_dense_matrix_clear( &
                     scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_x_im_dense(ihistory))
               END DO
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_e_im_dense, &
                           scf_env%qs_ot_restart%channel(ichannel)%rot_mat_h_x_im_dense)
            END IF
            IF (ASSOCIATED(scf_env%qs_ot_restart%channel(ichannel)%ener_x)) THEN
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ener_x)
            END IF
            IF (ASSOCIATED(scf_env%qs_ot_restart%channel(ichannel)%ener_dx)) THEN
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ener_dx)
            END IF
            IF (ASSOCIATED(scf_env%qs_ot_restart%channel(ichannel)%ener_gx_old)) THEN
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ener_gx_old)
            END IF
            IF (ASSOCIATED(scf_env%qs_ot_restart%channel(ichannel)%c_diis)) THEN
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%c_diis)
            END IF
            IF (ASSOCIATED(scf_env%qs_ot_restart%channel(ichannel)%c_broy)) THEN
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%c_broy)
            END IF
            IF (ASSOCIATED(scf_env%qs_ot_restart%channel(ichannel)%energy_h)) THEN
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%energy_h)
            END IF
            IF (ASSOCIATED(scf_env%qs_ot_restart%channel(ichannel)%lbfgs_rho)) THEN
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%lbfgs_rho, &
                           scf_env%qs_ot_restart%channel(ichannel)%lbfgs_yy, &
                           scf_env%qs_ot_restart%channel(ichannel)%lbfgs_sy_rotation, &
                           scf_env%qs_ot_restart%channel(ichannel)%lbfgs_yy_rotation)
            END IF
            IF (ASSOCIATED(scf_env%qs_ot_restart%channel(ichannel)%ipivot)) THEN
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ipivot)
            END IF
            IF (ASSOCIATED(scf_env%qs_ot_restart%channel(ichannel)%ls_diis)) THEN
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ls_diis)
            END IF
            IF (ASSOCIATED(scf_env%qs_ot_restart%channel(ichannel)%lss_diis)) THEN
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%lss_diis)
            END IF
            IF (ASSOCIATED(scf_env%qs_ot_restart%channel(ichannel)%ener_h_e)) THEN
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ener_h_e)
            END IF
            IF (ASSOCIATED(scf_env%qs_ot_restart%channel(ichannel)%ener_h_x)) THEN
               DEALLOCATE (scf_env%qs_ot_restart%channel(ichannel)%ener_h_x)
            END IF
         END DO
         DEALLOCATE (scf_env%qs_ot_restart%channel)
      END IF
      scf_env%qs_ot_restart%valid = .FALSE.
      scf_env%qs_ot_restart%format_version = 0
      scf_env%qs_ot_restart%nchannel = 0
      NULLIFY (scf_env%qs_ot_restart%channel)

   END SUBROUTINE qs_scf_ot_restart_clear

! **************************************************************************************************
!> \brief release storage owned by a dense OT restart matrix
!> \param matrix restart matrix to clear
! **************************************************************************************************
   SUBROUTINE qs_ot_restart_dense_matrix_clear(matrix)
      TYPE(qs_ot_restart_dense_matrix_type), &
         INTENT(INOUT)                                   :: matrix

      IF (ALLOCATED(matrix%data)) DEALLOCATE (matrix%data)

   END SUBROUTINE qs_ot_restart_dense_matrix_clear

! **************************************************************************************************
!> \brief  creates krylov space
!> \param krylov_space ...
!> \param scf_section ...
!> \par History
!>      05.2009 created [MI]
!> \author [MI]
! **************************************************************************************************
   SUBROUTINE krylov_space_create(krylov_space, scf_section)

      TYPE(krylov_space_type), POINTER                   :: krylov_space
      TYPE(section_vals_type), POINTER                   :: scf_section

      CPASSERT(.NOT. ASSOCIATED(krylov_space))
      ALLOCATE (krylov_space)

      NULLIFY (krylov_space%c_eval, krylov_space%t_eval)
      NULLIFY (krylov_space%v_mat)
      NULLIFY (krylov_space%mo_conv, krylov_space%mo_refine)
      NULLIFY (krylov_space%chc_mat, krylov_space%c_vec)
      NULLIFY (krylov_space%tmp_mat)
      NULLIFY (krylov_space%block1_mat, krylov_space%block2_mat)
      NULLIFY (krylov_space%block3_mat, krylov_space%block4_mat, krylov_space%block5_mat)

      CALL section_vals_val_get(scf_section, "DIAGONALIZATION%MAX_ITER", &
                                i_val=krylov_space%max_iter)
      CALL section_vals_val_get(scf_section, "DIAGONALIZATION%KRYLOV%NKRYLOV", &
                                i_val=krylov_space%nkrylov)
      CALL section_vals_val_get(scf_section, "DIAGONALIZATION%KRYLOV%NBLOCK", &
                                i_val=krylov_space%nblock)
      CALL section_vals_val_get(scf_section, "DIAGONALIZATION%KRYLOV%EPS_KRYLOV", &
                                r_val=krylov_space%eps_conv)
      CALL section_vals_val_get(scf_section, "DIAGONALIZATION%KRYLOV%EPS_STD_DIAG", &
                                r_val=krylov_space%eps_std_diag)
      CALL section_vals_val_get(scf_section, "DIAGONALIZATION%EPS_ADAPT", &
                                r_val=krylov_space%eps_adapt)
      CALL section_vals_val_get(scf_section, "DIAGONALIZATION%KRYLOV%CHECK_MOS_CONV", &
                                l_val=krylov_space%always_check_conv)

   END SUBROUTINE krylov_space_create

! **************************************************************************************************
!> \brief releases krylov space
!> \param krylov_space ...
!> \par History
!>      05.2009 created [MI]
!> \author [MI]
! **************************************************************************************************
   SUBROUTINE krylov_space_release(krylov_space)
      TYPE(krylov_space_type), POINTER                   :: krylov_space

      IF (ASSOCIATED(krylov_space)) THEN

         DEALLOCATE (krylov_space%c_eval)
         DEALLOCATE (krylov_space%t_eval)

         CALL cp_fm_release(krylov_space%v_mat)
         CALL cp_fm_release(krylov_space%mo_conv)
         CALL cp_fm_release(krylov_space%mo_refine)
         CALL cp_fm_release(krylov_space%chc_mat)
         CALL cp_fm_release(krylov_space%c_vec)

         IF (ASSOCIATED(krylov_space%tmp_mat)) THEN
            CALL cp_fm_release(krylov_space%tmp_mat)
            DEALLOCATE (krylov_space%tmp_mat)
         END IF
         IF (ASSOCIATED(krylov_space%block1_mat)) THEN
            CALL cp_fm_release(krylov_space%block1_mat)
            DEALLOCATE (krylov_space%block1_mat)
         END IF
         IF (ASSOCIATED(krylov_space%block2_mat)) THEN
            CALL cp_fm_release(krylov_space%block2_mat)
            DEALLOCATE (krylov_space%block2_mat)
         END IF
         IF (ASSOCIATED(krylov_space%block3_mat)) THEN
            CALL cp_fm_release(krylov_space%block3_mat)
            DEALLOCATE (krylov_space%block3_mat)
         END IF
         IF (ASSOCIATED(krylov_space%block4_mat)) THEN
            CALL cp_fm_release(krylov_space%block4_mat)
            DEALLOCATE (krylov_space%block4_mat)
         END IF
         IF (ASSOCIATED(krylov_space%block5_mat)) THEN
            CALL cp_fm_release(krylov_space%block5_mat)
            DEALLOCATE (krylov_space%block5_mat)
         END IF

         DEALLOCATE (krylov_space)

         NULLIFY (krylov_space)
      END IF

   END SUBROUTINE krylov_space_release

! **************************************************************************************************
!> \brief creates subspace-rotation environment
!> \param subspace_env ...
!> \param scf_section ...
!> \param ecut ...
!> \par History
!>      09.2009 created [MI]
!> \author [MI]
! **************************************************************************************************
   SUBROUTINE diag_subspace_env_create(subspace_env, scf_section, ecut)

      TYPE(subspace_env_type), POINTER                   :: subspace_env
      TYPE(section_vals_type), POINTER                   :: scf_section
      REAL(dp), INTENT(IN)                               :: ecut

      LOGICAL                                            :: do_mixing
      TYPE(section_vals_type), POINTER                   :: mixing_section

      CPASSERT(.NOT. ASSOCIATED(subspace_env))
      ALLOCATE (subspace_env)

      NULLIFY (subspace_env%p_matrix_store)
      NULLIFY (subspace_env%p_matrix_mix)
      NULLIFY (subspace_env%chc_mat)
      NULLIFY (subspace_env%c_vec)
      NULLIFY (subspace_env%c0)
      NULLIFY (subspace_env%mixing_store)
      NULLIFY (mixing_section)

      CALL section_vals_val_get(scf_section, "DIAGONALIZATION%DIAG_SUB_SCF%MAX_ITER", &
                                i_val=subspace_env%max_iter)
      CALL section_vals_val_get(scf_section, "DIAGONALIZATION%DIAG_SUB_SCF%EPS_ENE", &
                                r_val=subspace_env%eps_ene)
      CALL section_vals_val_get(scf_section, "DIAGONALIZATION%DIAG_SUB_SCF%EPS_SKIP_SUB_DIAG", &
                                r_val=subspace_env%eps_diag_sub)
      CALL section_vals_val_get(scf_section, "DIAGONALIZATION%DIAG_SUB_SCF%EPS_ADAPT_SCF", &
                                r_val=subspace_env%eps_adapt)
      subspace_env%mixing_method = 0
      do_mixing = .FALSE.
      mixing_section => section_vals_get_subs_vals(scf_section, "DIAGONALIZATION%DIAG_SUB_SCF%MIXING")
      CALL section_vals_val_get(mixing_section, "_SECTION_PARAMETERS_", &
                                l_val=do_mixing)
      IF (do_mixing) THEN
         CALL section_vals_val_get(mixing_section, "METHOD", &
                                   i_val=subspace_env%mixing_method)
         IF (subspace_env%mixing_method >= direct_mixing_nr) THEN
            ALLOCATE (subspace_env%mixing_store)
            CALL mixing_storage_create(subspace_env%mixing_store, mixing_section, &
                                       subspace_env%mixing_method, ecut=ecut)
         END IF
      END IF

   END SUBROUTINE diag_subspace_env_create

! **************************************************************************************************
!> \brief releases subspace-rotation environment
!> \param subspace_env ...
!> \par History
!>      09.2009 created [MI]
!> \author [MI]
! **************************************************************************************************
   SUBROUTINE diag_subspace_env_release(subspace_env)
      TYPE(subspace_env_type), POINTER                   :: subspace_env

      IF (ASSOCIATED(subspace_env)) THEN

         IF (ASSOCIATED(subspace_env%p_matrix_store)) THEN

            CPASSERT(.TRUE.)
            CALL dbcsr_deallocate_matrix_set(subspace_env%p_matrix_store)
         END IF
         CALL cp_fm_release(subspace_env%chc_mat)
         CALL cp_fm_release(subspace_env%c_vec)
         CALL cp_fm_release(subspace_env%c0)

         IF (ASSOCIATED(subspace_env%mixing_store)) THEN
            CALL mixing_storage_release(subspace_env%mixing_store)
            DEALLOCATE (subspace_env%mixing_store)
         END IF

         DEALLOCATE (subspace_env)
      END IF

   END SUBROUTINE diag_subspace_env_release

END MODULE qs_scf_types
