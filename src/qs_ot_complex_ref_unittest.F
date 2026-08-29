!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Regressions for complex k-point OT derivatives and preconditioner guards.
! **************************************************************************************************
PROGRAM qs_ot_complex_ref_unittest
   USE bibliography,                    ONLY: add_all_references
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_element,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_set_element,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_finalize_lib, dbcsr_get_block_p, dbcsr_init_lib, dbcsr_put_block, &
        dbcsr_release, dbcsr_reserve_blocks, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_logger_create,&
                                              cp_logger_release,&
                                              cp_logger_type,&
                                              cp_rm_default_logger
   USE input_constants,                 ONLY: &
        ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single_inverse, &
        ot_precond_none, ot_precond_s_inverse, ot_precond_solver_chebyshev, &
        ot_precond_solver_default, ot_precond_solver_direct, ot_precond_solver_inv_chol
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: default_output_unit
   USE mathlib,                         ONLY: diag_complex
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE preconditioner_makes,            ONLY: make_complex_full_kinetic,&
                                              make_complex_full_s_inverse,&
                                              make_complex_full_single_inverse
   USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                              init_preconditioner,&
                                              preconditioner_type
   USE qs_ot,                           ONLY: qs_ot_apply_complex_frechet_dbcsr,&
                                              qs_ot_complex_exp_frechet_kernel,&
                                              qs_ot_generate_rotation,&
                                              qs_ot_get_derivative_complex,&
                                              qs_ot_get_derivative_ref_complex,&
                                              qs_ot_get_orbitals_ref_complex,&
                                              qs_ot_get_p_complex,&
                                              qs_ot_rot_mat_derivative
   USE qs_ot_types,                     ONLY: qs_ot_kpoint_preconditioner_scale,&
                                              qs_ot_kpoint_preconditioner_solver_supported,&
                                              qs_ot_kpoint_preconditioner_supported,&
                                              qs_ot_type
   USE reference_manager,               ONLY: remove_all_references

   IMPLICIT NONE

   INTEGER, PARAMETER                                 :: k = 2, n = 3

   COMPLEX(KIND=dp), DIMENSION(n, k)                  :: direction, gradient, gradient_old, hc, q, x
   COMPLEX(KIND=dp), DIMENSION(n, n)                  :: hamiltonian
   COMPLEX(KIND=dp), DIMENSION(k, k)                  :: b, inverse_sqrt
   INTEGER                                            :: io_unit, mynode, nfail
   INTEGER, DIMENSION(:), POINTER                     :: col_dist, col_k, row_dist, row_k, row_n
   INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
   LOGICAL                                            :: found
   REAL(KIND=dp)                                      :: error, fd_slope, hstep, old_error, slope
   REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_p
   TYPE(cp_logger_type), POINTER                      :: logger
   TYPE(dbcsr_distribution_type)                     :: dist_kk, dist_nk
   TYPE(dbcsr_type), TARGET                           :: b_im, b_re, f_im, f_re, g_im, g_re, &
                                                         hc_im, hc_re, q_im, q_re, sx_im, sx_re, &
                                                         tmp_kk, x_im, x_re
   TYPE(dbcsr_type), POINTER                          :: hc_im_p, hc_re_p
   TYPE(mp_comm_type)                                 :: mp_comm
   TYPE(mp_para_env_type), POINTER                    :: para_env
   TYPE(qs_ot_type)                                   :: ot_env

   nfail = 0
   CALL test_preconditioner_support(nfail)
   NULLIFY (hc_im_p, hc_re_p, logger, para_env)
   CALL mp_world_init(mp_comm)
   mynode = mp_comm%mepos
   io_unit = -1
   IF (mynode == 0) io_unit = default_output_unit
   ALLOCATE (para_env)
   CALL para_env%from_dup(mp_comm)
   CALL cp_logger_create(logger, para_env=para_env, default_global_unit_nr=io_unit, &
                         close_global_unit_on_dealloc=.FALSE.)
   CALL cp_add_default_logger(logger)
   CALL add_all_references()
   CALL dbcsr_init_lib(mp_comm%get_handle(), io_unit)
   CALL test_real_rotation_frechet(para_env, nfail)
   CALL test_complex_rotation_frechet(para_env, nfail)
   CALL test_sparse_frechet_patterns(para_env, nfail)
   CALL test_complex_ref_rotation(para_env, nfail)
   CALL test_complex_preconditioner_gauge(para_env, nfail)

   x(:, 1) = [CMPLX(1.10_dp, 0.10_dp, KIND=dp), CMPLX(0.20_dp, -0.30_dp, KIND=dp), &
              CMPLX(-0.10_dp, 0.20_dp, KIND=dp)]
   x(:, 2) = [CMPLX(0.30_dp, 0.40_dp, KIND=dp), CMPLX(0.90_dp, -0.20_dp, KIND=dp), &
              CMPLX(0.25_dp, 0.10_dp, KIND=dp)]
   direction(:, 1) = [CMPLX(0.17_dp, -0.09_dp, KIND=dp), CMPLX(-0.21_dp, 0.13_dp, KIND=dp), &
                      CMPLX(0.08_dp, 0.19_dp, KIND=dp)]
   direction(:, 2) = [CMPLX(-0.12_dp, 0.23_dp, KIND=dp), CMPLX(0.16_dp, 0.07_dp, KIND=dp), &
                      CMPLX(-0.14_dp, -0.11_dp, KIND=dp)]

   hamiltonian(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
   hamiltonian(1, 1) = CMPLX(0.70_dp, 0.0_dp, KIND=dp)
   hamiltonian(2, 2) = CMPLX(-0.40_dp, 0.0_dp, KIND=dp)
   hamiltonian(3, 3) = CMPLX(1.10_dp, 0.0_dp, KIND=dp)
   hamiltonian(1, 2) = CMPLX(0.20_dp, 0.10_dp, KIND=dp)
   hamiltonian(2, 1) = CONJG(hamiltonian(1, 2))
   hamiltonian(1, 3) = CMPLX(-0.10_dp, 0.05_dp, KIND=dp)
   hamiltonian(3, 1) = CONJG(hamiltonian(1, 3))
   hamiltonian(2, 3) = CMPLX(0.15_dp, 0.20_dp, KIND=dp)
   hamiltonian(3, 2) = CONJG(hamiltonian(2, 3))

   CALL dense_inverse_sqrt(MATMUL(CONJG(TRANSPOSE(x)), x), inverse_sqrt)
   q = MATMUL(x, inverse_sqrt)
   hc = 2.0_dp*MATMUL(hamiltonian, q)

   ALLOCATE (pgrid(0:0, 0:0), row_dist(1), col_dist(1), row_n(1), row_k(1), col_k(1))
   pgrid(:, :) = 0
   row_dist(:) = 0
   col_dist(:) = 0
   row_n(:) = n
   row_k(:) = k
   col_k(:) = k
   CALL dbcsr_distribution_new(dist_nk, group=mp_comm_self%get_handle(), pgrid=pgrid, &
                               row_dist=row_dist, col_dist=col_dist)
   CALL dbcsr_distribution_new(dist_kk, group=mp_comm_self%get_handle(), pgrid=pgrid, &
                               row_dist=row_dist, col_dist=col_dist)
   CALL test_complex_strict_derivative(dist_nk, dist_kk, row_n, row_k, col_k, para_env, nfail)

   CALL create_one_block(x_re, "complex_ref_x_re", dist_nk, row_n, col_k)
   CALL create_one_block(x_im, "complex_ref_x_im", dist_nk, row_n, col_k)
   CALL create_one_block(sx_re, "complex_ref_sx_re", dist_nk, row_n, col_k)
   CALL create_one_block(sx_im, "complex_ref_sx_im", dist_nk, row_n, col_k)
   CALL create_one_block(hc_re, "complex_ref_hc_re", dist_nk, row_n, col_k)
   CALL create_one_block(hc_im, "complex_ref_hc_im", dist_nk, row_n, col_k)
   CALL create_one_block(g_re, "complex_ref_g_re", dist_nk, row_n, col_k)
   CALL create_one_block(g_im, "complex_ref_g_im", dist_nk, row_n, col_k)
   CALL create_one_block(q_re, "complex_ref_q_re", dist_nk, row_n, col_k)
   CALL create_one_block(q_im, "complex_ref_q_im", dist_nk, row_n, col_k)
   CALL create_one_block(f_re, "complex_ref_f_re", dist_kk, row_k, col_k)
   CALL create_one_block(f_im, "complex_ref_f_im", dist_kk, row_k, col_k)
   CALL create_one_block(b_re, "complex_ref_b_re", dist_kk, row_k, col_k)
   CALL create_one_block(b_im, "complex_ref_b_im", dist_kk, row_k, col_k)
   CALL create_one_block(tmp_kk, "complex_ref_tmp_kk", dist_kk, row_k, col_k)

   CALL put_complex_pair(x_re, x_im, x)
   CALL put_complex_pair(sx_re, sx_im, x)
   CALL put_complex_pair(hc_re, hc_im, hc)
   CALL put_complex_pair(f_re, f_im, inverse_sqrt)

   ot_env%has_complex_kpoint_state = .TRUE.
   ot_env%settings%eps_irac_filter_matrix = 0.0_dp
   ot_env%matrix_x => x_re
   ot_env%matrix_x_im => x_im
   ot_env%matrix_sx => sx_re
   ot_env%matrix_sx_im => sx_im
   ot_env%matrix_gx => g_re
   ot_env%matrix_gx_im => g_im
   ot_env%matrix_ref_inv_sqrt => f_re
   ot_env%matrix_ref_inv_sqrt_im => f_im
   ot_env%buf1_k_k_sym => b_re
   ot_env%buf2_k_k_sym => b_im
   ot_env%buf3_k_k_sym => tmp_kk
   ot_env%buf1_n_k => q_re
   ot_env%buf1_n_k_dp => q_im

   hc_re_p => hc_re
   hc_im_p => hc_im
   CALL qs_ot_get_derivative_ref_complex(hc_re_p, hc_im_p, ot_env)

   gradient(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
   CALL dbcsr_get_block_p(g_re, 1, 1, block_p, found)
   IF (found) gradient = gradient + CMPLX(block_p, 0.0_dp, KIND=dp)
   CALL dbcsr_get_block_p(g_im, 1, 1, block_p, found)
   IF (found) gradient = gradient + CMPLX(0.0_dp, block_p, KIND=dp)

   hstep = 1.0E-6_dp
   fd_slope = (polar_energy(x + hstep*direction, hamiltonian) - &
               polar_energy(x - hstep*direction, hamiltonian))/(2.0_dp*hstep)
   slope = REAL(SUM(CONJG(gradient)*direction), KIND=dp)
   error = ABS(slope - fd_slope)

   b = MATMUL(CONJG(TRANSPOSE(x)), hc)
   gradient_old = hc - MATMUL(x, b)
   old_error = ABS(REAL(SUM(CONJG(gradient_old)*direction), KIND=dp) - fd_slope)
   IF (error > 5.0E-8_dp .OR. old_error < 1.0E-3_dp) nfail = nfail + 1

   IF (mynode == 0) THEN
      WRITE (io_unit, '(A,3(1X,ES13.6))') "complex finite-REF pullback: fd/predicted/error", &
         fd_slope, slope, error
      WRITE (io_unit, '(A,1X,ES13.6)') "discarded projected-gradient error", old_error
   END IF

   CALL dbcsr_release(b_im)
   CALL dbcsr_release(b_re)
   CALL dbcsr_release(tmp_kk)
   CALL dbcsr_release(f_im)
   CALL dbcsr_release(f_re)
   CALL dbcsr_release(q_im)
   CALL dbcsr_release(q_re)
   CALL dbcsr_release(g_im)
   CALL dbcsr_release(g_re)
   CALL dbcsr_release(hc_im)
   CALL dbcsr_release(hc_re)
   CALL dbcsr_release(sx_im)
   CALL dbcsr_release(sx_re)
   CALL dbcsr_release(x_im)
   CALL dbcsr_release(x_re)
   CALL dbcsr_distribution_release(dist_kk)
   CALL dbcsr_distribution_release(dist_nk)
   DEALLOCATE (pgrid, row_dist, col_dist, row_n, row_k, col_k)

   CALL dbcsr_finalize_lib()
   CALL remove_all_references()
   CALL cp_rm_default_logger()
   CALL cp_logger_release(logger)
   CALL mp_para_env_release(para_env)
   CALL mp_world_finalize()

   IF (nfail > 0) ERROR STOP "qs_ot_complex_ref_unittest failed"

CONTAINS

! **************************************************************************************************
!> \brief Check the real finite rotation pullback used by Gamma-point REF OT.
!> \param para_env parallel environment
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_real_rotation_frechet(para_env, nfail)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: m = 3

      INTEGER                                            :: i
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, col_size, row_dist, row_size
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: error, fd_slope, hstep, predicted
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      REAL(KIND=dp), DIMENSION(m), TARGET                :: rotation_evals
      REAL(KIND=dp), DIMENSION(m, m)                     :: dedu, direction, generator, gradient, &
                                                            rotation_minus, rotation_plus
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_type), DIMENSION(7), TARGET             :: matrices
      TYPE(mp_para_env_type), POINTER                    :: local_para_env
      TYPE(qs_ot_type)                                   :: rotation_env

      NULLIFY (blacs_env, local_para_env)
      generator(:, :) = 0.0_dp
      generator(1, 2) = 0.23_dp
      generator(2, 1) = -generator(1, 2)
      generator(1, 3) = -0.11_dp
      generator(3, 1) = -generator(1, 3)
      generator(2, 3) = 0.17_dp
      generator(3, 2) = -generator(2, 3)

      direction(:, :) = 0.0_dp
      direction(1, 2) = -0.09_dp
      direction(2, 1) = -direction(1, 2)
      direction(1, 3) = 0.14_dp
      direction(3, 1) = -direction(1, 3)
      direction(2, 3) = 0.06_dp
      direction(3, 2) = -direction(2, 3)

      dedu = RESHAPE([0.31_dp, -0.27_dp, 0.08_dp, &
                      0.19_dp, 0.42_dp, -0.16_dp, &
                      -0.07_dp, 0.21_dp, 0.35_dp], [m, m])

      ALLOCATE (local_para_env)
      CALL local_para_env%from_dup(mp_comm_self)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=local_para_env)
      ALLOCATE (pgrid(0:0, 0:0), row_dist(1), col_dist(1), row_size(1), col_size(1))
      pgrid(:, :) = 0
      row_dist(:) = 0
      col_dist(:) = 0
      row_size(:) = m
      col_size(:) = m
      CALL dbcsr_distribution_new(dist, group=mp_comm_self%get_handle(), pgrid=pgrid, &
                                  row_dist=row_dist, col_dist=col_dist)
      DO i = 1, SIZE(matrices)
         CALL create_one_block(matrices(i), 'real_rotation', dist, row_size, col_size)
      END DO
      CALL dbcsr_put_block(matrices(2), 1, 1, generator)
      CALL dbcsr_put_block(matrices(3), 1, 1, dedu)

      rotation_env%para_env => local_para_env
      rotation_env%blacs_env => blacs_env
      rotation_env%rot_mat_u => matrices(1)
      rotation_env%rot_mat_x => matrices(2)
      rotation_env%rot_mat_dedu => matrices(3)
      rotation_env%rot_mat_gx => matrices(4)
      rotation_env%rot_mat_evec_re => matrices(5)
      rotation_env%rot_mat_evec_im => matrices(6)
      rotation_env%matrix_buf1 => matrices(7)
      rotation_env%rot_mat_evals => rotation_evals

      CALL qs_ot_generate_rotation(rotation_env)
      CALL qs_ot_rot_mat_derivative(rotation_env)
      gradient(:, :) = 0.0_dp
      CALL dbcsr_get_block_p(matrices(4), 1, 1, block, found)
      IF (found) gradient = block

      hstep = 1.0E-6_dp
      rotation_plus = REAL(dense_antihermitian_exp( &
                           CMPLX(generator + hstep*direction, 0.0_dp, KIND=dp)), KIND=dp)
      rotation_minus = REAL(dense_antihermitian_exp( &
                            CMPLX(generator - hstep*direction, 0.0_dp, KIND=dp)), KIND=dp)
      fd_slope = SUM(dedu*(rotation_plus - rotation_minus))/(2.0_dp*hstep)
      predicted = 0.5_dp*SUM(gradient*direction)
      error = ABS(predicted - fd_slope)
      IF (error > 5.0E-8_dp) nfail = nfail + 1
      IF (para_env%is_source()) THEN
         WRITE (io_unit, '(A,3(1X,ES13.6))') 'real finite rotation slope fd/predicted/error', &
            fd_slope, predicted, error
      END IF

      DO i = 1, SIZE(matrices)
         CALL dbcsr_release(matrices(i))
      END DO
      CALL dbcsr_distribution_release(dist)
      DEALLOCATE (pgrid, row_dist, col_dist, row_size, col_size)
      CALL cp_blacs_env_release(blacs_env)
      CALL mp_para_env_release(local_para_env)

   END SUBROUTINE test_real_rotation_frechet

! **************************************************************************************************
!> \brief Check the combined finite complex REF and rotation derivatives.
!> \param para_env parallel environment
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_complex_ref_rotation(para_env, nfail)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(INOUT)                             :: nfail

      COMPLEX(KIND=dp), DIMENSION(k, k)                  :: direction_a, generator, gradient_a, &
                                                            inverse_sqrt, rotation
      COMPLEX(KIND=dp), DIMENSION(n, k)                  :: c_dense, c_reference, direction_x, &
                                                            gradient_x, hc_dense, q_dense, x_dense
      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: hamiltonian, overlap
      INTEGER                                            :: i, j
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, col_k, col_n, row_dist, row_k, &
                                                            row_n
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      REAL(KIND=dp)                                      :: error_c, error_ref, error_rot, fd_ref, &
                                                            fd_rot, hstep, predicted_ref, &
                                                            predicted_rot
      REAL(KIND=dp), DIMENSION(k)                        :: occupation
      REAL(KIND=dp), DIMENSION(k), TARGET                :: rotation_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_distribution_type)                      :: dist_kk, dist_nk, dist_nn
      TYPE(dbcsr_type), DIMENSION(10), TARGET            :: rot_matrix
      TYPE(dbcsr_type), DIMENSION(2), TARGET             :: nk_buffer
      TYPE(dbcsr_type), DIMENSION(3), TARGET             :: nonsym_buffer
      TYPE(dbcsr_type), DIMENSION(4), TARGET             :: sym_buffer
      TYPE(dbcsr_type), POINTER                          :: c_im_p, c_re_p, hc_im_p, hc_re_p, &
                                                            s_im_p, s_re_p
      TYPE(dbcsr_type), TARGET                           :: c_im, c_re, f_im, f_re, g_im, g_re, &
                                                            hc_im, hc_re, s_im, s_re, sx_im, &
                                                            sx_re, x_im, x_re
      TYPE(mp_para_env_type), POINTER                    :: local_para_env
      TYPE(qs_ot_type)                                   :: rotation_env

      NULLIFY (blacs_env, c_im_p, c_re_p, hc_im_p, hc_re_p, local_para_env, &
               s_im_p, s_re_p)

      x_dense(:, 1) = [CMPLX(1.05_dp, 0.12_dp, KIND=dp), CMPLX(0.18_dp, -0.27_dp, KIND=dp), &
                       CMPLX(-0.14_dp, 0.16_dp, KIND=dp)]
      x_dense(:, 2) = [CMPLX(0.26_dp, 0.31_dp, KIND=dp), CMPLX(0.94_dp, -0.17_dp, KIND=dp), &
                       CMPLX(0.21_dp, 0.13_dp, KIND=dp)]
      direction_x(:, 1) = [CMPLX(0.11_dp, -0.07_dp, KIND=dp), CMPLX(-0.16_dp, 0.09_dp, KIND=dp), &
                           CMPLX(0.05_dp, 0.14_dp, KIND=dp)]
      direction_x(:, 2) = [CMPLX(-0.08_dp, 0.15_dp, KIND=dp), CMPLX(0.13_dp, 0.04_dp, KIND=dp), &
                           CMPLX(-0.10_dp, -0.06_dp, KIND=dp)]

      generator(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      generator(1, 1) = CMPLX(0.0_dp, 0.19_dp, KIND=dp)
      generator(2, 2) = CMPLX(0.0_dp, -0.08_dp, KIND=dp)
      generator(1, 2) = CMPLX(0.23_dp, -0.11_dp, KIND=dp)
      generator(2, 1) = -CONJG(generator(1, 2))
      direction_a(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      direction_a(1, 1) = CMPLX(0.0_dp, -0.07_dp, KIND=dp)
      direction_a(2, 2) = CMPLX(0.0_dp, 0.05_dp, KIND=dp)
      direction_a(1, 2) = CMPLX(-0.09_dp, 0.06_dp, KIND=dp)
      direction_a(2, 1) = -CONJG(direction_a(1, 2))

      hamiltonian(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      hamiltonian(1, 1) = CMPLX(0.65_dp, 0.0_dp, KIND=dp)
      hamiltonian(2, 2) = CMPLX(-0.35_dp, 0.0_dp, KIND=dp)
      hamiltonian(3, 3) = CMPLX(1.15_dp, 0.0_dp, KIND=dp)
      hamiltonian(1, 2) = CMPLX(0.18_dp, 0.12_dp, KIND=dp)
      hamiltonian(1, 3) = CMPLX(-0.09_dp, 0.07_dp, KIND=dp)
      hamiltonian(2, 3) = CMPLX(0.14_dp, -0.17_dp, KIND=dp)
      DO i = 1, n
         DO j = i + 1, n
            hamiltonian(j, i) = CONJG(hamiltonian(i, j))
         END DO
      END DO
      overlap(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, n
         overlap(i, i) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      END DO
      occupation = [1.0_dp, 0.37_dp]

      CALL dense_inverse_sqrt(MATMUL(CONJG(TRANSPOSE(x_dense)), x_dense), inverse_sqrt)
      q_dense = MATMUL(x_dense, inverse_sqrt)
      rotation = dense_antihermitian_exp(generator)
      c_reference = MATMUL(q_dense, rotation)

      ALLOCATE (local_para_env)
      CALL local_para_env%from_dup(mp_comm_self)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=local_para_env)

      ALLOCATE (pgrid(0:0, 0:0), row_dist(1), col_dist(1), row_n(1), row_k(1), &
                col_n(1), col_k(1))
      pgrid(:, :) = 0
      row_dist(:) = 0
      col_dist(:) = 0
      row_n(:) = n
      col_n(:) = n
      row_k(:) = k
      col_k(:) = k
      CALL dbcsr_distribution_new(dist_nk, group=mp_comm_self%get_handle(), pgrid=pgrid, &
                                  row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_distribution_new(dist_nn, group=mp_comm_self%get_handle(), pgrid=pgrid, &
                                  row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_distribution_new(dist_kk, group=mp_comm_self%get_handle(), pgrid=pgrid, &
                                  row_dist=row_dist, col_dist=col_dist)

      CALL create_one_block(x_re, 'rotation_ref_x_re', dist_nk, row_n, col_k)
      CALL create_one_block(x_im, 'rotation_ref_x_im', dist_nk, row_n, col_k)
      CALL create_one_block(sx_re, 'rotation_ref_sx_re', dist_nk, row_n, col_k)
      CALL create_one_block(sx_im, 'rotation_ref_sx_im', dist_nk, row_n, col_k)
      CALL create_one_block(c_re, 'rotation_ref_c_re', dist_nk, row_n, col_k)
      CALL create_one_block(c_im, 'rotation_ref_c_im', dist_nk, row_n, col_k)
      CALL create_one_block(hc_re, 'rotation_ref_hc_re', dist_nk, row_n, col_k)
      CALL create_one_block(hc_im, 'rotation_ref_hc_im', dist_nk, row_n, col_k)
      CALL create_one_block(g_re, 'rotation_ref_g_re', dist_nk, row_n, col_k)
      CALL create_one_block(g_im, 'rotation_ref_g_im', dist_nk, row_n, col_k)
      CALL create_one_block(s_re, 'rotation_ref_s_re', dist_nn, row_n, col_n)
      CALL create_one_block(s_im, 'rotation_ref_s_im', dist_nn, row_n, col_n)
      CALL create_one_block(f_re, 'rotation_ref_f_re', dist_kk, row_k, col_k)
      CALL create_one_block(f_im, 'rotation_ref_f_im', dist_kk, row_k, col_k)
      DO i = 1, SIZE(sym_buffer)
         CALL create_one_block(sym_buffer(i), 'rotation_ref_sym', dist_kk, row_k, col_k)
      END DO
      DO i = 1, SIZE(nonsym_buffer)
         CALL create_one_block(nonsym_buffer(i), 'rotation_ref_nonsym', dist_kk, row_k, col_k)
      END DO
      DO i = 1, SIZE(nk_buffer)
         CALL create_one_block(nk_buffer(i), 'rotation_ref_nk', dist_nk, row_n, col_k)
      END DO
      DO i = 1, SIZE(rot_matrix)
         CALL create_one_block(rot_matrix(i), 'rotation_ref_rot', dist_kk, row_k, col_k)
      END DO

      CALL put_complex_pair(x_re, x_im, x_dense)
      CALL put_complex_pair(s_re, s_im, overlap)
      CALL put_complex_pair(rot_matrix(3), rot_matrix(4), generator)

      rotation_env%has_complex_kpoint_state = .TRUE.
      rotation_env%settings%do_rotation = .TRUE.
      rotation_env%settings%eps_irac_filter_matrix = 0.0_dp
      rotation_env%settings%ot_method = 'LBFG'
      rotation_env%para_env => local_para_env
      rotation_env%blacs_env => blacs_env
      rotation_env%matrix_x => x_re
      rotation_env%matrix_x_im => x_im
      rotation_env%matrix_sx => sx_re
      rotation_env%matrix_sx_im => sx_im
      rotation_env%matrix_gx => g_re
      rotation_env%matrix_gx_im => g_im
      rotation_env%matrix_ref_inv_sqrt => f_re
      rotation_env%matrix_ref_inv_sqrt_im => f_im
      rotation_env%buf1_k_k_sym => sym_buffer(1)
      rotation_env%buf2_k_k_sym => sym_buffer(2)
      rotation_env%buf3_k_k_sym => sym_buffer(3)
      rotation_env%buf4_k_k_sym => sym_buffer(4)
      rotation_env%buf1_k_k_nosym => nonsym_buffer(1)
      rotation_env%buf2_k_k_nosym => nonsym_buffer(2)
      rotation_env%buf3_k_k_nosym => nonsym_buffer(3)
      rotation_env%buf1_n_k => nk_buffer(1)
      rotation_env%buf1_n_k_dp => nk_buffer(2)
      rotation_env%rot_mat_u => rot_matrix(1)
      rotation_env%rot_mat_u_im => rot_matrix(2)
      rotation_env%rot_mat_x => rot_matrix(3)
      rotation_env%rot_mat_x_im => rot_matrix(4)
      rotation_env%rot_mat_dedu => rot_matrix(5)
      rotation_env%rot_mat_dedu_im => rot_matrix(6)
      rotation_env%rot_mat_gx => rot_matrix(7)
      rotation_env%rot_mat_gx_im => rot_matrix(8)
      rotation_env%rot_mat_evec_re => rot_matrix(9)
      rotation_env%rot_mat_evec_im => rot_matrix(10)
      rotation_env%rot_mat_evals => rotation_evals

      c_re_p => c_re
      c_im_p => c_im
      s_re_p => s_re
      s_im_p => s_im
      CALL qs_ot_get_orbitals_ref_complex(c_re_p, c_im_p, s_re_p, s_im_p, rotation_env)
      CALL get_complex_pair(c_re, c_im, c_dense)
      error_c = MAXVAL(ABS(c_dense - c_reference))

      hc_dense = 2.0_dp*MATMUL(hamiltonian, c_dense)
      DO j = 1, k
         hc_dense(:, j) = occupation(j)*hc_dense(:, j)
      END DO
      CALL put_complex_pair(hc_re, hc_im, hc_dense)
      hc_re_p => hc_re
      hc_im_p => hc_im
      CALL qs_ot_get_derivative_ref_complex(hc_re_p, hc_im_p, rotation_env)
      CALL get_complex_pair(g_re, g_im, gradient_x)
      CALL get_complex_pair(rot_matrix(7), rot_matrix(8), gradient_a)

      hstep = 1.0E-6_dp
      fd_ref = (rotated_polar_energy(x_dense + hstep*direction_x, generator, &
                                     hamiltonian, occupation) - &
                rotated_polar_energy(x_dense - hstep*direction_x, generator, &
                                     hamiltonian, occupation))/(2.0_dp*hstep)
      predicted_ref = REAL(SUM(CONJG(gradient_x)*direction_x), KIND=dp)
      error_ref = ABS(predicted_ref - fd_ref)

      fd_rot = (rotated_polar_energy(x_dense, generator + hstep*direction_a, &
                                     hamiltonian, occupation) - &
                rotated_polar_energy(x_dense, generator - hstep*direction_a, &
                                     hamiltonian, occupation))/(2.0_dp*hstep)
      predicted_rot = 0.5_dp*(SUM(REAL(gradient_a, KIND=dp)*REAL(direction_a, KIND=dp)) + &
                              SUM(AIMAG(gradient_a)*AIMAG(direction_a)))
      error_rot = ABS(predicted_rot - fd_rot)

      IF (error_c > 5.0E-11_dp .OR. error_ref > 5.0E-8_dp .OR. &
          error_rot > 5.0E-8_dp) nfail = nfail + 1
      IF (para_env%is_source()) THEN
         WRITE (io_unit, '(A,1X,ES13.6)') 'complex finite REF rotation orbital error', error_c
         WRITE (io_unit, '(A,3(1X,ES13.6))') 'complex finite REF slope fd/predicted/error', &
            fd_ref, predicted_ref, error_ref
         WRITE (io_unit, '(A,3(1X,ES13.6))') 'complex finite rotation slope fd/predicted/error', &
            fd_rot, predicted_rot, error_rot
      END IF

      DO i = 1, SIZE(rot_matrix)
         CALL dbcsr_release(rot_matrix(i))
      END DO
      DO i = 1, SIZE(nk_buffer)
         CALL dbcsr_release(nk_buffer(i))
      END DO
      DO i = 1, SIZE(nonsym_buffer)
         CALL dbcsr_release(nonsym_buffer(i))
      END DO
      DO i = 1, SIZE(sym_buffer)
         CALL dbcsr_release(sym_buffer(i))
      END DO
      CALL dbcsr_release(f_im)
      CALL dbcsr_release(f_re)
      CALL dbcsr_release(s_im)
      CALL dbcsr_release(s_re)
      CALL dbcsr_release(g_im)
      CALL dbcsr_release(g_re)
      CALL dbcsr_release(hc_im)
      CALL dbcsr_release(hc_re)
      CALL dbcsr_release(c_im)
      CALL dbcsr_release(c_re)
      CALL dbcsr_release(sx_im)
      CALL dbcsr_release(sx_re)
      CALL dbcsr_release(x_im)
      CALL dbcsr_release(x_re)
      CALL dbcsr_distribution_release(dist_kk)
      CALL dbcsr_distribution_release(dist_nn)
      CALL dbcsr_distribution_release(dist_nk)
      DEALLOCATE (pgrid, row_dist, col_dist, row_n, row_k, col_n, col_k)
      CALL cp_blacs_env_release(blacs_env)
      CALL mp_para_env_release(local_para_env)

   END SUBROUTINE test_complex_ref_rotation

! **************************************************************************************************
!> \brief Check the finite complex rotation Frechet derivative against central differences.
!> \param para_env parallel environment
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_complex_rotation_frechet(para_env, nfail)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: m = 3

      COMPLEX(KIND=dp), DIMENSION(m, m)                  :: direction, fd, generator, predicted
      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: error, hstep
      REAL(KIND=dp), DIMENSION(m)                        :: evals

      evals = [-0.7_dp, 0.2000003_dp, 1.1_dp]
      generator(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, m
         generator(i, i) = CMPLX(0.0_dp, -evals(i), KIND=dp)
      END DO
      direction(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      direction(1, 2) = CMPLX(0.17_dp, -0.09_dp, KIND=dp)
      direction(2, 1) = -CONJG(direction(1, 2))
      direction(1, 3) = CMPLX(-0.06_dp, 0.13_dp, KIND=dp)
      direction(3, 1) = -CONJG(direction(1, 3))
      direction(2, 3) = CMPLX(0.11_dp, 0.04_dp, KIND=dp)
      direction(3, 2) = -CONJG(direction(2, 3))
      direction(1, 1) = CMPLX(0.0_dp, 0.08_dp, KIND=dp)
      direction(2, 2) = CMPLX(0.0_dp, -0.03_dp, KIND=dp)
      direction(3, 3) = CMPLX(0.0_dp, 0.05_dp, KIND=dp)

      DO i = 1, m
         DO j = 1, m
            predicted(i, j) = qs_ot_complex_exp_frechet_kernel(evals(i), evals(j))* &
                              direction(i, j)
         END DO
      END DO
      hstep = 1.0E-6_dp
      fd = (dense_antihermitian_exp(generator + hstep*direction) - &
            dense_antihermitian_exp(generator - hstep*direction))/(2.0_dp*hstep)
      error = MAXVAL(ABS(fd - predicted))
      IF (error > 5.0E-9_dp) nfail = nfail + 1
      IF (para_env%is_source()) THEN
         WRITE (io_unit, '(A,1X,ES13.6)') 'complex rotation Frechet finite-difference error', error
      END IF

   END SUBROUTINE test_complex_rotation_frechet

! **************************************************************************************************
!> \brief Check union handling for distinct real/imaginary DBCSR block patterns.
!> \param para_env parallel environment
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_sparse_frechet_patterns(para_env, nfail)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(INOUT)                             :: nfail

      COMPLEX(KIND=dp)                                   :: expected, value
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, col_size, row_dist, row_size
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: error
      REAL(KIND=dp), DIMENSION(2)                        :: evals
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_type)                                   :: inner_im, inner_re, outer_im, outer_re

      ALLOCATE (pgrid(0:0, 0:0), row_dist(2), col_dist(2), row_size(2), col_size(2))
      pgrid(:, :) = 0
      row_dist(:) = 0
      col_dist(:) = 0
      row_size(:) = 1
      col_size(:) = 1
      CALL dbcsr_distribution_new(dist, group=mp_comm_self%get_handle(), pgrid=pgrid, &
                                  row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_create(inner_re, 'frechet_sparse_re', dist, dbcsr_type_no_symmetry, row_size, col_size)
      CALL dbcsr_reserve_blocks(inner_re, rows=[1], cols=[2])
      CALL dbcsr_finalize(inner_re)
      CALL dbcsr_create(inner_im, 'frechet_sparse_im', dist, dbcsr_type_no_symmetry, row_size, col_size)
      CALL dbcsr_reserve_blocks(inner_im, rows=[2], cols=[1])
      CALL dbcsr_finalize(inner_im)
      CALL dbcsr_put_block(inner_re, 1, 2, RESHAPE([0.37_dp], [1, 1]))
      CALL dbcsr_put_block(inner_im, 2, 1, RESHAPE([-0.21_dp], [1, 1]))

      evals = [-0.4_dp, 0.9_dp]
      CALL qs_ot_apply_complex_frechet_dbcsr(evals, inner_re, inner_im, outer_re, outer_im)
      error = 0.0_dp
      expected = qs_ot_complex_exp_frechet_kernel(evals(1), evals(2))*CMPLX(0.37_dp, 0.0_dp, KIND=dp)
      value = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 1, 2, block, found)
      IF (found) value = value + CMPLX(block(1, 1), 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 1, 2, block, found)
      IF (found) value = value + CMPLX(0.0_dp, block(1, 1), KIND=dp)
      error = MAX(error, ABS(value - expected))

      expected = qs_ot_complex_exp_frechet_kernel(evals(2), evals(1))*CMPLX(0.0_dp, -0.21_dp, KIND=dp)
      value = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_re, 2, 1, block, found)
      IF (found) value = value + CMPLX(block(1, 1), 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(outer_im, 2, 1, block, found)
      IF (found) value = value + CMPLX(0.0_dp, block(1, 1), KIND=dp)
      error = MAX(error, ABS(value - expected))
      IF (error > 5.0E-14_dp) nfail = nfail + 1
      IF (para_env%is_source()) THEN
         WRITE (io_unit, '(A,1X,ES13.6)') 'sparse complex Frechet block-union error', error
      END IF

      CALL dbcsr_release(outer_im)
      CALL dbcsr_release(outer_re)
      CALL dbcsr_release(inner_im)
      CALL dbcsr_release(inner_re)
      CALL dbcsr_distribution_release(dist)
      DEALLOCATE (pgrid, row_dist, col_dist, row_size, col_size)

   END SUBROUTINE test_sparse_frechet_patterns
! **************************************************************************************************
!> \brief Check the finite complex STRICT chart and its analytic derivative.
!> \param dist_nk distribution for orbital-shaped matrices
!> \param dist_kk distribution for occupied-space matrices
!> \param row_n AO block sizes
!> \param row_k occupied-space row block sizes
!> \param col_k occupied-space column block sizes
!> \param para_env parallel environment
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_complex_strict_derivative(dist_nk, dist_kk, row_n, row_k, col_k, para_env, nfail)
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: dist_nk, dist_kk
      INTEGER, DIMENSION(:), INTENT(INOUT), POINTER      :: row_n, row_k, col_k
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(INOUT)                             :: nfail

      CHARACTER(LEN=32)                                  :: name
      COMPLEX(KIND=dp), DIMENSION(k, k)                  :: inverse_sqrt
      COMPLEX(KIND=dp), DIMENSION(n, k)                  :: c0, direction, gradient, hc, orbitals, &
                                                            sc0, sx, x
      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: hamiltonian, overlap
      INTEGER                                            :: i
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: error, fd_slope, hstep, ortho_error, &
                                                            slope
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_type), DIMENSION(15), TARGET            :: nk
      TYPE(dbcsr_type), DIMENSION(25), TARGET            :: kk
      TYPE(dbcsr_type), POINTER                          :: hc_im_p, hc_re_p, sx_im_p, sx_re_p, &
                                                            x_im_p, x_re_p
      TYPE(mp_para_env_type), POINTER                    :: strict_para_env
      TYPE(qs_ot_type)                                   :: strict_env

      NULLIFY (blacs_env, hc_im_p, hc_re_p, strict_para_env, x_im_p, x_re_p, sx_im_p, sx_re_p)
      ALLOCATE (strict_para_env)
      CALL strict_para_env%from_dup(mp_comm_self)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=strict_para_env)
      DO i = 1, SIZE(nk)
         WRITE (name, '(A,I0)') 'complex_strict_nk_', i
         CALL create_one_block(nk(i), TRIM(name), dist_nk, row_n, col_k)
      END DO
      DO i = 1, SIZE(kk)
         WRITE (name, '(A,I0)') 'complex_strict_kk_', i
         CALL create_one_block(kk(i), TRIM(name), dist_kk, row_k, col_k)
      END DO

      overlap(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      overlap(1, 1) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      overlap(2, 2) = CMPLX(1.2_dp, 0.0_dp, KIND=dp)
      overlap(3, 3) = CMPLX(1.4_dp, 0.0_dp, KIND=dp)
      overlap(1, 2) = CMPLX(0.06_dp, 0.03_dp, KIND=dp)
      overlap(2, 1) = CONJG(overlap(1, 2))
      overlap(1, 3) = CMPLX(-0.04_dp, 0.02_dp, KIND=dp)
      overlap(3, 1) = CONJG(overlap(1, 3))
      overlap(2, 3) = CMPLX(0.05_dp, -0.01_dp, KIND=dp)
      overlap(3, 2) = CONJG(overlap(2, 3))

      c0(:, 1) = [CMPLX(1.0_dp, 0.10_dp, KIND=dp), CMPLX(0.12_dp, -0.08_dp, KIND=dp), &
                  CMPLX(-0.07_dp, 0.04_dp, KIND=dp)]
      c0(:, 2) = [CMPLX(0.09_dp, 0.05_dp, KIND=dp), CMPLX(0.85_dp, -0.12_dp, KIND=dp), &
                  CMPLX(0.11_dp, 0.06_dp, KIND=dp)]
      CALL dense_inverse_sqrt(MATMUL(CONJG(TRANSPOSE(c0)), MATMUL(overlap, c0)), inverse_sqrt)
      c0 = MATMUL(c0, inverse_sqrt)

      x(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      x(3, 1) = CMPLX(0.23_dp, 0.17_dp, KIND=dp)
      x(3, 2) = CMPLX(-0.14_dp, 0.21_dp, KIND=dp)
      x = x - MATMUL(c0, MATMUL(CONJG(TRANSPOSE(c0)), MATMUL(overlap, x)))
      direction(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      direction(3, 1) = CMPLX(0.11_dp, -0.09_dp, KIND=dp)
      direction(3, 2) = CMPLX(0.07_dp, 0.13_dp, KIND=dp)
      direction = direction - MATMUL(c0, MATMUL(CONJG(TRANSPOSE(c0)), MATMUL(overlap, direction)))

      hamiltonian(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      hamiltonian(1, 1) = CMPLX(-0.70_dp, 0.0_dp, KIND=dp)
      hamiltonian(2, 2) = CMPLX(-0.20_dp, 0.0_dp, KIND=dp)
      hamiltonian(3, 3) = CMPLX(0.90_dp, 0.0_dp, KIND=dp)
      hamiltonian(1, 2) = CMPLX(0.08_dp, 0.04_dp, KIND=dp)
      hamiltonian(2, 1) = CONJG(hamiltonian(1, 2))
      hamiltonian(1, 3) = CMPLX(-0.12_dp, 0.07_dp, KIND=dp)
      hamiltonian(3, 1) = CONJG(hamiltonian(1, 3))
      hamiltonian(2, 3) = CMPLX(0.05_dp, -0.09_dp, KIND=dp)
      hamiltonian(3, 2) = CONJG(hamiltonian(2, 3))

      CALL strict_orbitals(c0, x, overlap, orbitals)
      hc = 2.0_dp*MATMUL(hamiltonian, orbitals)
      sx = MATMUL(overlap, x)
      sc0 = MATMUL(overlap, c0)
      CALL put_complex_pair(nk(1), nk(2), x)
      CALL put_complex_pair(nk(3), nk(4), sx)
      CALL put_complex_pair(nk(5), nk(6), hc)
      CALL put_complex_pair(nk(9), nk(10), c0)
      CALL put_complex_pair(nk(11), nk(12), sc0)

      strict_env%has_complex_kpoint_state = .TRUE.
      strict_env%settings%ot_algorithm = 'TOD'
      strict_env%settings%do_rotation = .FALSE.
      strict_env%para_env => strict_para_env
      strict_env%blacs_env => blacs_env
      strict_env%matrix_x => nk(1)
      strict_env%matrix_x_im => nk(2)
      strict_env%matrix_sx => nk(3)
      strict_env%matrix_sx_im => nk(4)
      strict_env%matrix_gx => nk(7)
      strict_env%matrix_gx_im => nk(8)
      strict_env%matrix_c0 => nk(9)
      strict_env%matrix_c0_im => nk(10)
      strict_env%matrix_sc0 => nk(11)
      strict_env%matrix_sc0_im => nk(12)
      strict_env%matrix_buf_nk => nk(13)
      strict_env%matrix_buf_nk_im => nk(14)
      strict_env%matrix_tmp_nk => nk(15)
      strict_env%matrix_p => kk(1)
      strict_env%matrix_p_im => kk(2)
      strict_env%matrix_r => kk(3)
      strict_env%matrix_r_im => kk(4)
      strict_env%matrix_sinp => kk(5)
      strict_env%matrix_sinp_im => kk(6)
      strict_env%matrix_cosp => kk(7)
      strict_env%matrix_cosp_im => kk(8)
      strict_env%matrix_sinp_b => kk(9)
      strict_env%matrix_cosp_b => kk(10)
      strict_env%matrix_buf1 => kk(11)
      strict_env%matrix_buf1_im => kk(12)
      strict_env%matrix_buf2 => kk(13)
      strict_env%matrix_buf2_im => kk(14)
      strict_env%matrix_buf3 => kk(15)
      strict_env%matrix_buf3_im => kk(16)
      strict_env%matrix_buf4 => kk(17)
      strict_env%matrix_buf4_im => kk(18)
      strict_env%matrix_os => kk(19)
      strict_env%matrix_os_im => kk(20)
      strict_env%matrix_buf1_ortho => kk(21)
      strict_env%matrix_buf1_ortho_im => kk(22)
      strict_env%matrix_buf2_ortho => kk(23)
      strict_env%matrix_buf2_ortho_im => kk(24)
      strict_env%matrix_tmp_ortho => kk(25)
      ALLOCATE (strict_env%evals(k), strict_env%dum(k))

      x_re_p => nk(1)
      x_im_p => nk(2)
      sx_re_p => nk(3)
      sx_im_p => nk(4)
      hc_re_p => nk(5)
      hc_im_p => nk(6)
      CALL qs_ot_get_p_complex(x_re_p, x_im_p, sx_re_p, sx_im_p, strict_env)
      CALL qs_ot_get_derivative_complex(hc_re_p, hc_im_p, strict_env)

      gradient(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(nk(7), 1, 1, block, found)
      IF (found) gradient = gradient + CMPLX(block, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(nk(8), 1, 1, block, found)
      IF (found) gradient = gradient + CMPLX(0.0_dp, block, KIND=dp)

      hstep = 1.0E-6_dp
      fd_slope = (strict_energy(c0, x + hstep*direction, overlap, hamiltonian) - &
                  strict_energy(c0, x - hstep*direction, overlap, hamiltonian))/(2.0_dp*hstep)
      slope = REAL(SUM(CONJG(gradient)*direction), KIND=dp)
      error = ABS(slope - fd_slope)
      ortho_error = MAXVAL(ABS(MATMUL(CONJG(TRANSPOSE(orbitals)), MATMUL(overlap, orbitals)) - &
                               RESHAPE([CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
                                        CMPLX(0.0_dp, 0.0_dp, KIND=dp), &
                                        CMPLX(0.0_dp, 0.0_dp, KIND=dp), &
                                        CMPLX(1.0_dp, 0.0_dp, KIND=dp)], [k, k])))
      IF (error > 5.0E-8_dp .OR. ortho_error > 5.0E-13_dp) nfail = nfail + 1
      IF (para_env%is_source()) THEN
         WRITE (io_unit, '(A,3(1X,ES13.6))') 'complex finite-STRICT derivative: fd/predicted/error', &
            fd_slope, slope, error
         WRITE (io_unit, '(A,1X,ES13.6)') 'complex finite-STRICT orthonormality error', ortho_error
      END IF

      DEALLOCATE (strict_env%evals, strict_env%dum)
      DO i = SIZE(kk), 1, -1
         CALL dbcsr_release(kk(i))
      END DO
      DO i = SIZE(nk), 1, -1
         CALL dbcsr_release(nk(i))
      END DO
      CALL cp_blacs_env_release(blacs_env)
      CALL mp_para_env_release(strict_para_env)

   END SUBROUTINE test_complex_strict_derivative
! **************************************************************************************************
!> \brief Check gauge invariance, Hermiticity, and positivity of the complex inverse operator.
!> \param para_env parallel environment
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_complex_preconditioner_gauge(para_env, nfail)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(INOUT)                             :: nfail

      COMPLEX(KIND=dp)                                   :: value
      COMPLEX(KIND=dp), DIMENSION(n, k)                  :: c0, c0_rot
      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: eigenvectors, h, inverse, inverse_rot, &
                                                            inverse_s, inverse_t, reference_s, &
                                                            reference_t, s, t
      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: gauge_error, hermitian_error, s_error, &
                                                            t_error
      REAL(KIND=dp), DIMENSION(n)                        :: eigenvalues
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: matrix_c0, matrix_c0_rot, matrix_h, &
                                                            matrix_s, matrix_t
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct_nk, fmstruct_nn
      TYPE(preconditioner_type)                          :: preconditioner, preconditioner_rot, &
                                                            preconditioner_s, preconditioner_t

      NULLIFY (blacs_env, fmstruct_nk, fmstruct_nn)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)
      CALL cp_fm_struct_create(fmstruct=fmstruct_nn, para_env=para_env, context=blacs_env, &
                               nrow_global=n, ncol_global=n)
      CALL cp_fm_struct_create(fmstruct=fmstruct_nk, para_env=para_env, context=blacs_env, &
                               nrow_global=n, ncol_global=k)
      CALL cp_cfm_create(matrix_h, fmstruct_nn, name='complex preconditioner H')
      CALL cp_cfm_create(matrix_s, fmstruct_nn, name='complex preconditioner S')
      CALL cp_cfm_create(matrix_t, fmstruct_nn, name='complex preconditioner T')
      CALL cp_cfm_create(matrix_c0, fmstruct_nk, name='complex preconditioner C')
      CALL cp_cfm_create(matrix_c0_rot, fmstruct_nk, name='complex preconditioner C U')
      CALL cp_cfm_set_all(matrix_h, CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      CALL cp_cfm_set_all(matrix_s, CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      CALL cp_cfm_set_all(matrix_t, CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      CALL cp_cfm_set_all(matrix_c0, CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      CALL cp_cfm_set_all(matrix_c0_rot, CMPLX(0.0_dp, 0.0_dp, KIND=dp))

      h(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      h(1, 1) = CMPLX(-0.35_dp, 0.0_dp, KIND=dp)
      h(2, 2) = CMPLX(0.15_dp, 0.0_dp, KIND=dp)
      h(3, 3) = CMPLX(0.90_dp, 0.0_dp, KIND=dp)
      h(1, 2) = CMPLX(0.08_dp, 0.04_dp, KIND=dp)
      h(2, 1) = CONJG(h(1, 2))
      h(1, 3) = CMPLX(-0.06_dp, 0.03_dp, KIND=dp)
      h(3, 1) = CONJG(h(1, 3))
      h(2, 3) = CMPLX(0.05_dp, -0.07_dp, KIND=dp)
      h(3, 2) = CONJG(h(2, 3))
      s(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      s(1, 1) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      s(2, 2) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      s(3, 3) = CMPLX(1.4_dp, 0.0_dp, KIND=dp)
      s(1, 3) = CMPLX(0.08_dp, 0.03_dp, KIND=dp)
      s(3, 1) = CONJG(s(1, 3))
      s(2, 3) = CMPLX(-0.04_dp, 0.02_dp, KIND=dp)
      s(3, 2) = CONJG(s(2, 3))
      t = h
      DO i = 1, n
         t(i, i) = t(i, i) + CMPLX(1.50_dp, 0.0_dp, KIND=dp)
      END DO
      c0(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      c0(1, 1) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      c0(2, 2) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      c0_rot(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      c0_rot(1, 1) = CMPLX(SQRT(0.5_dp), 0.0_dp, KIND=dp)
      c0_rot(1, 2) = CMPLX(0.0_dp, SQRT(0.5_dp), KIND=dp)
      c0_rot(2, 1) = CMPLX(0.0_dp, SQRT(0.5_dp), KIND=dp)
      c0_rot(2, 2) = CMPLX(SQRT(0.5_dp), 0.0_dp, KIND=dp)

      DO j = 1, n
         DO i = 1, n
            CALL cp_cfm_set_element(matrix_h, i, j, h(i, j))
            CALL cp_cfm_set_element(matrix_s, i, j, s(i, j))
            CALL cp_cfm_set_element(matrix_t, i, j, t(i, j))
         END DO
      END DO
      DO j = 1, k
         DO i = 1, n
            CALL cp_cfm_set_element(matrix_c0, i, j, c0(i, j))
            CALL cp_cfm_set_element(matrix_c0_rot, i, j, c0_rot(i, j))
         END DO
      END DO

      CALL init_preconditioner(preconditioner, para_env, blacs_env)
      CALL init_preconditioner(preconditioner_rot, para_env, blacs_env)
      CALL init_preconditioner(preconditioner_s, para_env, blacs_env)
      CALL init_preconditioner(preconditioner_t, para_env, blacs_env)
      CALL make_complex_full_single_inverse(preconditioner, matrix_c0, matrix_h, matrix_s, 0.20_dp)
      CALL make_complex_full_single_inverse(preconditioner_rot, matrix_c0_rot, matrix_h, matrix_s, 0.20_dp)
      CALL make_complex_full_s_inverse(preconditioner_s, matrix_s)
      CALL make_complex_full_kinetic(preconditioner_t, matrix_t, matrix_s, 0.20_dp)
      DO j = 1, n
         DO i = 1, n
            CALL cp_cfm_get_element(preconditioner%complex_fm, i, j, value)
            inverse(i, j) = value
            CALL cp_cfm_get_element(preconditioner_rot%complex_fm, i, j, value)
            inverse_rot(i, j) = value
            CALL cp_cfm_get_element(preconditioner_s%complex_fm, i, j, value)
            inverse_s(i, j) = value
            CALL cp_cfm_get_element(preconditioner_t%complex_fm, i, j, value)
            inverse_t(i, j) = value
         END DO
      END DO
      gauge_error = MAXVAL(ABS(inverse - inverse_rot))
      hermitian_error = MAXVAL(ABS(inverse - CONJG(TRANSPOSE(inverse))))
      CALL diag_complex(inverse, eigenvectors, eigenvalues)
      CALL dense_hermitian_inverse(s, reference_s)
      CALL dense_hermitian_inverse(t + 0.20_dp*s, reference_t)
      s_error = MAXVAL(ABS(inverse_s - reference_s))
      t_error = MAXVAL(ABS(inverse_t - reference_t))
      IF (gauge_error > 1.0E-11_dp .OR. hermitian_error > 1.0E-12_dp .OR. &
          MINVAL(eigenvalues) <= 0.0_dp) nfail = nfail + 1
      IF (s_error > 1.0E-11_dp .OR. t_error > 1.0E-11_dp) nfail = nfail + 1
      IF (para_env%is_source()) THEN
         WRITE (io_unit, '(A,3(1X,ES13.6))') 'complex FULL_SINGLE gauge/hermitian/min-eigenvalue', &
            gauge_error, hermitian_error, MINVAL(eigenvalues)
         WRITE (io_unit, '(A,2(1X,ES13.6))') 'complex S/KINETIC inverse errors', s_error, t_error
      END IF

      CALL destroy_preconditioner(preconditioner_t)
      CALL destroy_preconditioner(preconditioner_s)
      CALL destroy_preconditioner(preconditioner_rot)
      CALL destroy_preconditioner(preconditioner)
      CALL cp_cfm_release(matrix_c0_rot)
      CALL cp_cfm_release(matrix_c0)
      CALL cp_cfm_release(matrix_t)
      CALL cp_cfm_release(matrix_s)
      CALL cp_cfm_release(matrix_h)
      CALL cp_fm_struct_release(fmstruct_nk)
      CALL cp_fm_struct_release(fmstruct_nn)
      CALL cp_blacs_env_release(blacs_env)

   END SUBROUTINE test_complex_preconditioner_gauge

! **************************************************************************************************
!> \brief Check the intentionally narrow complex k-point preconditioner support matrix.
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_preconditioner_support(nfail)
      INTEGER, INTENT(INOUT)                             :: nfail

      IF (.NOT. qs_ot_kpoint_preconditioner_supported(ot_precond_none, .FALSE.)) nfail = nfail + 1
      IF (.NOT. qs_ot_kpoint_preconditioner_supported(ot_precond_full_single_inverse, .FALSE.)) THEN
         nfail = nfail + 1
      END IF
      IF (qs_ot_kpoint_preconditioner_supported(ot_precond_full_all, .FALSE.)) nfail = nfail + 1
      IF (.NOT. qs_ot_kpoint_preconditioner_supported(ot_precond_full_kinetic, .FALSE.)) nfail = nfail + 1
      IF (.NOT. qs_ot_kpoint_preconditioner_supported(ot_precond_s_inverse, .FALSE.)) nfail = nfail + 1
      IF (qs_ot_kpoint_preconditioner_supported(ot_precond_full_single_inverse, .TRUE.)) nfail = nfail + 1

      IF (.NOT. qs_ot_kpoint_preconditioner_solver_supported( &
          ot_precond_full_single_inverse, ot_precond_solver_default, .FALSE.)) nfail = nfail + 1
      IF (.NOT. qs_ot_kpoint_preconditioner_solver_supported( &
          ot_precond_full_single_inverse, ot_precond_solver_inv_chol, .FALSE.)) nfail = nfail + 1
      IF (qs_ot_kpoint_preconditioner_solver_supported( &
          ot_precond_full_single_inverse, ot_precond_solver_direct, .FALSE.)) nfail = nfail + 1
      IF (qs_ot_kpoint_preconditioner_solver_supported( &
          ot_precond_full_single_inverse, ot_precond_solver_chebyshev, .FALSE.)) nfail = nfail + 1
      IF (.NOT. qs_ot_kpoint_preconditioner_solver_supported( &
          ot_precond_full_kinetic, ot_precond_solver_default, .FALSE.)) nfail = nfail + 1
      IF (.NOT. qs_ot_kpoint_preconditioner_solver_supported( &
          ot_precond_full_kinetic, ot_precond_solver_inv_chol, .FALSE.)) nfail = nfail + 1
      IF (qs_ot_kpoint_preconditioner_solver_supported( &
          ot_precond_full_kinetic, ot_precond_solver_direct, .FALSE.)) nfail = nfail + 1
      IF (.NOT. qs_ot_kpoint_preconditioner_solver_supported( &
          ot_precond_s_inverse, ot_precond_solver_default, .FALSE.)) nfail = nfail + 1
      IF (.NOT. qs_ot_kpoint_preconditioner_solver_supported( &
          ot_precond_s_inverse, ot_precond_solver_inv_chol, .FALSE.)) nfail = nfail + 1
      IF (qs_ot_kpoint_preconditioner_solver_supported( &
          ot_precond_s_inverse, ot_precond_solver_chebyshev, .FALSE.)) nfail = nfail + 1
      IF (ABS(qs_ot_kpoint_preconditioner_scale(0.25_dp) - 4.0_dp) > EPSILON(1.0_dp)) THEN
         nfail = nfail + 1
      END IF

   END SUBROUTINE test_preconditioner_support

! **************************************************************************************************
!> \brief Create and reserve a one-block DBCSR matrix.
!> \param matrix ...
!> \param name ...
!> \param dist ...
!> \param row_size ...
!> \param col_size ...
! **************************************************************************************************
   SUBROUTINE create_one_block(matrix, name, dist, row_size, col_size)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: dist
      INTEGER, DIMENSION(:), INTENT(INOUT), POINTER      :: row_size, col_size

      CALL dbcsr_create(matrix, name, dist, dbcsr_type_no_symmetry, row_size, col_size)
      CALL dbcsr_reserve_blocks(matrix, rows=[1], cols=[1])
      CALL dbcsr_finalize(matrix)

   END SUBROUTINE create_one_block

! **************************************************************************************************
!> \brief Store a dense complex matrix in a real/imaginary DBCSR pair.
!> \param matrix_re ...
!> \param matrix_im ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE put_complex_pair(matrix_re, matrix_im, matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_re, matrix_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: matrix

      CALL dbcsr_put_block(matrix_re, 1, 1, REAL(matrix, KIND=dp))
      CALL dbcsr_put_block(matrix_im, 1, 1, AIMAG(matrix))

   END SUBROUTINE put_complex_pair

! **************************************************************************************************
!> \brief Retrieve a dense complex matrix from a one-block DBCSR pair.
!> \param matrix_re ...
!> \param matrix_im ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE get_complex_pair(matrix_re, matrix_im, matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_re, matrix_im
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: matrix

      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block

      matrix(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(matrix_re, 1, 1, block, found)
      IF (found) matrix = matrix + CMPLX(block, 0.0_dp, KIND=dp)
      CALL dbcsr_get_block_p(matrix_im, 1, 1, block, found)
      IF (found) matrix = matrix + CMPLX(0.0_dp, block, KIND=dp)

   END SUBROUTINE get_complex_pair

! **************************************************************************************************
!> \brief Dense fixed-occupation energy after finite REF orthogonalization and rotation.
!> \param x REF coordinate
!> \param generator anti-Hermitian rotation generator
!> \param hamiltonian Hermitian one-particle Hamiltonian
!> \param occupation orbital occupations
!> \return energy
! **************************************************************************************************
   FUNCTION rotated_polar_energy(x, generator, hamiltonian, occupation) RESULT(energy)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: x, generator, hamiltonian
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation
      REAL(KIND=dp)                                      :: energy

      COMPLEX(KIND=dp), &
         DIMENSION(SIZE(x, 1), SIZE(x, 2))               :: orbitals
      COMPLEX(KIND=dp), &
         DIMENSION(SIZE(x, 2), SIZE(x, 2))               :: inverse_sqrt, rotation
      INTEGER                                            :: i

      CALL dense_inverse_sqrt(MATMUL(CONJG(TRANSPOSE(x)), x), inverse_sqrt)
      rotation = dense_antihermitian_exp(generator)
      orbitals = MATMUL(MATMUL(x, inverse_sqrt), rotation)
      energy = 0.0_dp
      DO i = 1, SIZE(occupation)
         energy = energy + occupation(i)* &
                  REAL(DOT_PRODUCT(orbitals(:, i), MATMUL(hamiltonian, orbitals(:, i))), KIND=dp)
      END DO

   END FUNCTION rotated_polar_energy

! **************************************************************************************************
!> \brief Dense reference exponential for an anti-Hermitian matrix.
!> \param generator anti-Hermitian generator
!> \return its unitary exponential
! **************************************************************************************************
   FUNCTION dense_antihermitian_exp(generator) RESULT(rotation)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: generator
      COMPLEX(KIND=dp), DIMENSION(SIZE(generator, 1), &
         SIZE(generator, 2))                             :: rotation

      COMPLEX(KIND=dp), DIMENSION(SIZE(generator, 1), &
         SIZE(generator, 2))                             :: hermitian, vectors
      INTEGER                                            :: i
      REAL(KIND=dp), DIMENSION(SIZE(generator, 1))       :: eigenvalues

      hermitian = CMPLX(0.0_dp, 1.0_dp, KIND=dp)*generator
      CALL diag_complex(hermitian, vectors, eigenvalues)
      rotation(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, SIZE(eigenvalues)
         rotation = rotation + EXP(CMPLX(0.0_dp, -eigenvalues(i), KIND=dp))* &
                    MATMUL(RESHAPE(vectors(:, i), [SIZE(generator, 1), 1]), &
                           RESHAPE(CONJG(vectors(:, i)), [1, SIZE(generator, 1)]))
      END DO

   END FUNCTION dense_antihermitian_exp

! **************************************************************************************************
!> \brief Hermitian inverse for an independent dense reference.
!> \param matrix ...
!> \param inverse ...
! **************************************************************************************************
   SUBROUTINE dense_hermitian_inverse(matrix, inverse)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: matrix
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: inverse

      COMPLEX(KIND=dp), &
         DIMENSION(SIZE(matrix, 1), SIZE(matrix, 2))     :: vectors
      INTEGER                                            :: i
      REAL(KIND=dp), DIMENSION(SIZE(matrix, 1))          :: eigenvalues

      CALL diag_complex(matrix, vectors, eigenvalues)
      inverse(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, SIZE(eigenvalues)
         inverse = inverse + MATMUL(RESHAPE(vectors(:, i), [SIZE(matrix, 1), 1]), &
                                    RESHAPE(CONJG(vectors(:, i)), [1, SIZE(matrix, 1)]))/ &
                   eigenvalues(i)
      END DO

   END SUBROUTINE dense_hermitian_inverse

! **************************************************************************************************
!> \brief Hermitian inverse square root for the independent dense reference.
!> \param matrix ...
!> \param inverse ...
! **************************************************************************************************
   SUBROUTINE dense_inverse_sqrt(matrix, inverse)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: matrix
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: inverse

      COMPLEX(KIND=dp), &
         DIMENSION(SIZE(matrix, 1), SIZE(matrix, 2))     :: vectors
      INTEGER                                            :: i
      REAL(KIND=dp), DIMENSION(SIZE(matrix, 1))          :: eigenvalues

      CALL diag_complex(matrix, vectors, eigenvalues)
      inverse(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, SIZE(eigenvalues)
         inverse = inverse + MATMUL(RESHAPE(vectors(:, i), [SIZE(matrix, 1), 1]), &
                                    RESHAPE(CONJG(vectors(:, i)), [1, SIZE(matrix, 1)]))/ &
                   SQRT(eigenvalues(i))
      END DO

   END SUBROUTINE dense_inverse_sqrt

! **************************************************************************************************
!> \brief Energy of the polar-orthonormalized occupied subspace.
!> \param coordinate ...
!> \param hamiltonian ...
!> \return ...
! **************************************************************************************************
   FUNCTION polar_energy(coordinate, hamiltonian) RESULT(energy)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: coordinate, hamiltonian
      REAL(KIND=dp)                                      :: energy

      COMPLEX(KIND=dp), DIMENSION(SIZE(coordinate, 1), &
         SIZE(coordinate, 2))                            :: orbitals
      COMPLEX(KIND=dp), DIMENSION(SIZE(coordinate, 2), &
         SIZE(coordinate, 2))                            :: inverse

      CALL dense_inverse_sqrt(MATMUL(CONJG(TRANSPOSE(coordinate)), coordinate), inverse)
      orbitals = MATMUL(coordinate, inverse)
      energy = REAL(SUM(CONJG(orbitals)*MATMUL(hamiltonian, orbitals)), KIND=dp)

   END FUNCTION polar_energy

! **************************************************************************************************
!> \brief Independent dense finite STRICT transformation for a Hermitian overlap.
!> \param c0 reference orbitals
!> \param coordinate tangent coordinate
!> \param overlap Hermitian overlap matrix
!> \param orbitals transformed orbitals
! **************************************************************************************************
   SUBROUTINE strict_orbitals(c0, coordinate, overlap, orbitals)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: c0, coordinate, overlap
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: orbitals

      COMPLEX(KIND=dp), DIMENSION(SIZE(coordinate, 2), &
         SIZE(coordinate, 2))                            :: cosp, sinp, vectors
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: root
      REAL(KIND=dp), DIMENSION(SIZE(coordinate, 2))      :: eigenvalues

      CALL diag_complex(MATMUL(CONJG(TRANSPOSE(coordinate)), MATMUL(overlap, coordinate)), &
                        vectors, eigenvalues)
      cosp(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      sinp(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, SIZE(eigenvalues)
         root = SQRT(MAX(0.0_dp, eigenvalues(i)))
         cosp = cosp + COS(root)*MATMUL(RESHAPE(vectors(:, i), [SIZE(vectors, 1), 1]), &
                                        RESHAPE(CONJG(vectors(:, i)), [1, SIZE(vectors, 1)]))
         IF (root < 1.0E-12_dp) THEN
            sinp = sinp + MATMUL(RESHAPE(vectors(:, i), [SIZE(vectors, 1), 1]), &
                                 RESHAPE(CONJG(vectors(:, i)), [1, SIZE(vectors, 1)]))
         ELSE
            sinp = sinp + SIN(root)/root* &
                   MATMUL(RESHAPE(vectors(:, i), [SIZE(vectors, 1), 1]), &
                          RESHAPE(CONJG(vectors(:, i)), [1, SIZE(vectors, 1)]))
         END IF
      END DO
      orbitals = MATMUL(c0, cosp) + MATMUL(coordinate, sinp)

   END SUBROUTINE strict_orbitals

! **************************************************************************************************
!> \brief Energy of the independently transformed finite STRICT orbitals.
!> \param c0 reference orbitals
!> \param coordinate tangent coordinate
!> \param overlap Hermitian overlap matrix
!> \param hamiltonian Hermitian Hamiltonian
!> \return energy
! **************************************************************************************************
   FUNCTION strict_energy(c0, coordinate, overlap, hamiltonian) RESULT(energy)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: c0, coordinate, overlap, hamiltonian
      REAL(KIND=dp)                                      :: energy

      COMPLEX(KIND=dp), DIMENSION(SIZE(coordinate, 1), &
         SIZE(coordinate, 2))                            :: orbitals

      CALL strict_orbitals(c0, coordinate, overlap, orbitals)
      energy = REAL(SUM(CONJG(orbitals)*MATMUL(hamiltonian, orbitals)), KIND=dp)

   END FUNCTION strict_energy

END PROGRAM qs_ot_complex_ref_unittest
