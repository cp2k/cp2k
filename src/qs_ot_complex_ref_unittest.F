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
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
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
   USE kpoint_methods,                  ONLY: kpoint_ot_energy_weighted_density
   USE machine,                         ONLY: default_output_unit
   USE mathlib,                         ONLY: diag_complex
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE preconditioner,                  ONLY: apply_preconditioner
   USE preconditioner_makes,            ONLY: make_complex_full_all,&
                                              make_complex_full_kinetic,&
                                              make_complex_full_s_inverse,&
                                              make_complex_full_single_inverse
   USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                              init_preconditioner,&
                                              preconditioner_type
   USE qs_ot,                           ONLY: &
        qs_ot_antihermitian_spectral_norm, qs_ot_apply_complex_frechet_dbcsr, &
        qs_ot_complex_exp_frechet_kernel, qs_ot_density_secant_hessian, &
        qs_ot_density_secant_orbital_overlaps, qs_ot_density_secant_projected_hessian, &
        qs_ot_density_tangent, qs_ot_finite_rotation_response, qs_ot_fixed_n_energy_gradient, &
        qs_ot_fixed_n_energy_hessian, qs_ot_fixed_n_projector_frechet, qs_ot_fixed_n_schur_block, &
        qs_ot_generate_rotation, qs_ot_get_derivative_complex, qs_ot_get_derivative_ref_complex, &
        qs_ot_get_orbitals_ref_complex, qs_ot_get_p_complex, qs_ot_projected_response_update, &
        qs_ot_rot_mat_derivative, qs_ot_symmetric_abs_solve, qs_ot_symmetric_sr1_update
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
   CALL test_antihermitian_spectral_norm(mynode, nfail)
   CALL test_fixed_n_mermin_energy(mynode, nfail)
   CALL test_fixed_n_projector_frechet(mynode, nfail)
   CALL test_finite_rotation_response(mynode, nfail)
   CALL test_fixed_n_rotation_schur(mynode, nfail)
   CALL test_symmetric_abs_solve(mynode, nfail)
   CALL test_projected_response_update(mynode, nfail)
   CALL test_symmetric_sr1_update(mynode, nfail)
   CALL test_density_secant_hessian(mynode, nfail)
   CALL test_density_tangent(mynode, nfail)
   CALL test_density_secant_moving_subspace(mynode, nfail)
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
   CALL test_kpoint_ot_energy_weighted_density(para_env, nfail)

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
!> \brief Check the noncanonical complex OT energy-weighted density and its gauge covariance.
!> \param para_env parallel environment
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_kpoint_ot_energy_weighted_density(para_env, nfail)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nao = 4, nmo = 3

      COMPLEX(KIND=dp), DIMENSION(nao, nao)              :: actual, expected, hmat, rotated
      COMPLEX(KIND=dp), DIMENSION(nao, nmo)              :: coeff, coeff_rot, hc_dense, hc_rot
      COMPLEX(KIND=dp), DIMENSION(nmo, nmo)              :: hblock, lagrange, rotation
      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: error, gauge_error, hermitian_error, &
                                                            phase, pi
      REAL(KIND=dp), DIMENSION(nao, nao)                 :: actual_im, actual_re, rotated_im, &
                                                            rotated_re
      REAL(KIND=dp), DIMENSION(nmo)                      :: occupation, uniform_occupation
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: coeff_struct, hmat_struct
      TYPE(cp_fm_type)                                   :: coeff_im_fm, coeff_re_fm, hc_im_fm, &
                                                            hc_re_fm, wmat_im_fm, wmat_re_fm

      NULLIFY (blacs_env, coeff_struct, hmat_struct)
      pi = ACOS(-1.0_dp)
      hmat(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      hmat(1, 1) = CMPLX(-0.8_dp, 0.0_dp, KIND=dp)
      hmat(2, 2) = CMPLX(0.1_dp, 0.0_dp, KIND=dp)
      hmat(3, 3) = CMPLX(0.7_dp, 0.0_dp, KIND=dp)
      hmat(4, 4) = CMPLX(1.4_dp, 0.0_dp, KIND=dp)
      hmat(1, 2) = CMPLX(0.21_dp, 0.13_dp, KIND=dp)
      hmat(1, 3) = CMPLX(-0.17_dp, 0.09_dp, KIND=dp)
      hmat(1, 4) = CMPLX(0.08_dp, -0.14_dp, KIND=dp)
      hmat(2, 3) = CMPLX(0.31_dp, -0.11_dp, KIND=dp)
      hmat(2, 4) = CMPLX(-0.19_dp, -0.07_dp, KIND=dp)
      hmat(3, 4) = CMPLX(0.27_dp, 0.16_dp, KIND=dp)
      DO i = 1, nao
         DO j = i + 1, nao
            hmat(j, i) = CONJG(hmat(i, j))
         END DO
      END DO
      DO j = 1, nmo
         DO i = 1, nao
            phase = 2.0_dp*pi*REAL((i - 1)*(j - 1), KIND=dp)/REAL(nao, KIND=dp)
            coeff(i, j) = CMPLX(COS(phase), SIN(phase), KIND=dp)/SQRT(REAL(nao, KIND=dp))
         END DO
      END DO
      hc_dense = MATMUL(hmat, coeff)
      occupation = [1.0_dp, 0.63_dp, 0.11_dp]

      hblock = MATMUL(CONJG(TRANSPOSE(coeff)), hc_dense)
      DO j = 1, nmo
         lagrange(:, j) = hblock(:, j)*occupation(j)
      END DO
      lagrange = 0.5_dp*(lagrange + CONJG(TRANSPOSE(lagrange)))
      expected = MATMUL(MATMUL(coeff, lagrange), CONJG(TRANSPOSE(coeff)))

      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)
      CALL cp_fm_struct_create(hmat_struct, nrow_global=nao, ncol_global=nao, &
                               context=blacs_env, para_env=para_env)
      CALL cp_fm_struct_create(coeff_struct, nrow_global=nao, ncol_global=nmo, &
                               context=blacs_env, para_env=para_env)
      CALL cp_fm_create(coeff_re_fm, coeff_struct)
      CALL cp_fm_create(coeff_im_fm, coeff_struct)
      CALL cp_fm_create(hc_re_fm, coeff_struct)
      CALL cp_fm_create(hc_im_fm, coeff_struct)
      CALL cp_fm_create(wmat_re_fm, hmat_struct)
      CALL cp_fm_create(wmat_im_fm, hmat_struct)
      CALL cp_fm_set_submatrix(coeff_re_fm, REAL(coeff, KIND=dp))
      CALL cp_fm_set_submatrix(coeff_im_fm, AIMAG(coeff))
      CALL cp_fm_set_submatrix(hc_re_fm, REAL(hc_dense, KIND=dp))
      CALL cp_fm_set_submatrix(hc_im_fm, AIMAG(hc_dense))

      CALL kpoint_ot_energy_weighted_density( &
         coeff_re_fm, coeff_im_fm, hc_re_fm, hc_im_fm, occupation, &
         wmat_re_fm, wmat_im_fm)
      CALL cp_fm_get_submatrix(wmat_re_fm, actual_re)
      CALL cp_fm_get_submatrix(wmat_im_fm, actual_im)
      actual = CMPLX(actual_re, actual_im, KIND=dp)
      error = MAXVAL(ABS(actual - expected))
      hermitian_error = MAXVAL(ABS(actual - CONJG(TRANSPOSE(actual))))

      rotation(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      rotation(1, 1) = CMPLX(SQRT(0.5_dp), 0.0_dp, KIND=dp)
      rotation(1, 2) = CMPLX(0.0_dp, SQRT(0.5_dp), KIND=dp)
      rotation(2, 1) = CMPLX(0.0_dp, SQRT(0.5_dp), KIND=dp)
      rotation(2, 2) = CMPLX(SQRT(0.5_dp), 0.0_dp, KIND=dp)
      rotation(3, 3) = CMPLX(0.0_dp, 1.0_dp, KIND=dp)
      coeff_rot = MATMUL(coeff, rotation)
      hc_rot = MATMUL(hc_dense, rotation)
      uniform_occupation(:) = 1.0_dp
      CALL cp_fm_set_submatrix(coeff_re_fm, REAL(coeff_rot, KIND=dp))
      CALL cp_fm_set_submatrix(coeff_im_fm, AIMAG(coeff_rot))
      CALL cp_fm_set_submatrix(hc_re_fm, REAL(hc_rot, KIND=dp))
      CALL cp_fm_set_submatrix(hc_im_fm, AIMAG(hc_rot))
      CALL kpoint_ot_energy_weighted_density( &
         coeff_re_fm, coeff_im_fm, hc_re_fm, hc_im_fm, uniform_occupation, &
         wmat_re_fm, wmat_im_fm)
      CALL cp_fm_get_submatrix(wmat_re_fm, rotated_re)
      CALL cp_fm_get_submatrix(wmat_im_fm, rotated_im)
      rotated = CMPLX(rotated_re, rotated_im, KIND=dp)
      expected = MATMUL(MATMUL(coeff, hblock), CONJG(TRANSPOSE(coeff)))
      gauge_error = MAXVAL(ABS(rotated - expected))

      IF (error > 5.0E-12_dp .OR. hermitian_error > 5.0E-12_dp .OR. &
          gauge_error > 5.0E-12_dp) nfail = nfail + 1
      IF (para_env%is_source()) THEN
         WRITE (*, '(A,3(1X,ES13.6))') 'complex OT W error/hermitian/gauge', &
            error, hermitian_error, gauge_error
      END IF

      CALL cp_fm_release(wmat_im_fm)
      CALL cp_fm_release(wmat_re_fm)
      CALL cp_fm_release(hc_im_fm)
      CALL cp_fm_release(hc_re_fm)
      CALL cp_fm_release(coeff_im_fm)
      CALL cp_fm_release(coeff_re_fm)
      CALL cp_fm_struct_release(coeff_struct)
      CALL cp_fm_struct_release(hmat_struct)
      CALL cp_blacs_env_release(blacs_env)

   END SUBROUTINE test_kpoint_ot_energy_weighted_density
! **************************************************************************************************
!> \brief Check the complex rotation norm against an analytic pair and a unitary gauge change.
!> \param mynode MPI rank
!> \param nfail accumulated failures
! **************************************************************************************************
   SUBROUTINE test_antihermitian_spectral_norm(mynode, nfail)
      INTEGER, INTENT(IN)                                :: mynode
      INTEGER, INTENT(INOUT)                             :: nfail

      COMPLEX(KIND=dp), DIMENSION(3)                     :: phase
      COMPLEX(KIND=dp), DIMENSION(3, 3)                  :: generator, transformed
      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: error, norm, transformed_norm

      generator(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      generator(1, 2) = CMPLX(0.3_dp, 0.4_dp, KIND=dp)
      generator(2, 1) = -CONJG(generator(1, 2))
      phase = EXP(CMPLX(0.0_dp, 1.0_dp, KIND=dp)*[0.37_dp, -0.51_dp, 0.83_dp])
      DO j = 1, 3
         DO i = 1, 3
            transformed(i, j) = CONJG(phase(i))*generator(i, j)*phase(j)
         END DO
      END DO

      norm = qs_ot_antihermitian_spectral_norm(generator)
      transformed_norm = qs_ot_antihermitian_spectral_norm(transformed)
      error = MAX(ABS(norm - 0.5_dp), ABS(transformed_norm - norm))
      IF (mynode == 0) WRITE (default_output_unit, '(A,ES14.6)') &
         "complex rotation spectral-norm error: ", error
      IF (error > 1.0E-13_dp) nfail = nfail + 1

   END SUBROUTINE test_antihermitian_spectral_norm

! **************************************************************************************************
!> \brief Check the fixed-N Mermin gradient, dense gauge projector, and occupation Hessian.
!> \param mynode MPI rank
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_fixed_n_mermin_energy(mynode, nfail)
      INTEGER, INTENT(IN)                                :: mynode
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nbands = 4
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-6_dp, maxocc = 2.0_dp, &
                                                            target_electrons = 2.25_dp, &
                                                            temperature = 0.07_dp

      REAL(KIND=dp)                                      :: fd_slope, gradient_error, hessian_error, &
                                                            mu, predicted_slope, response_sum, &
                                                            weighted_residual
      REAL(KIND=dp), DIMENSION(nbands) :: direction, energies, gradient, gradient_minus, &
         gradient_plus, hessian_action, occupation, rayleigh, response, response_minus, &
         response_plus, weights
      REAL(KIND=dp), DIMENSION(nbands, nbands)           :: hessian

      energies(:) = [-0.24_dp, -0.05_dp, 0.08_dp, 0.31_dp]
      rayleigh(:) = [-0.19_dp, -0.02_dp, 0.04_dp, 0.27_dp]
      weights(:) = [0.35_dp, 0.35_dp, 0.65_dp, 0.65_dp]
      direction(:) = [0.17_dp, -0.11_dp, 0.08_dp, -0.05_dp]

      CALL fixed_n_fermi_occupations(energies, weights, target_electrons, temperature, &
                                     maxocc, occupation, mu)
      response(:) = weights(:)*occupation(:)*(maxocc - occupation(:))/(maxocc*temperature)
      response_sum = SUM(response)
      weighted_residual = DOT_PRODUCT(response, rayleigh - energies)
      CALL qs_ot_fixed_n_energy_gradient(rayleigh, energies, response, response_sum, &
                                         weighted_residual, gradient)

      fd_slope = (fixed_n_mermin_value(rayleigh, energies + fd_step*direction, weights, &
                                       target_electrons, temperature, maxocc) - &
                  fixed_n_mermin_value(rayleigh, energies - fd_step*direction, weights, &
                                       target_electrons, temperature, maxocc))/(2.0_dp*fd_step)
      predicted_slope = DOT_PRODUCT(gradient, direction)
      gradient_error = MAX(ABS(fd_slope - predicted_slope), ABS(SUM(gradient)))

      ! At the canonical point h_i=e_i, the exact local Hessian is the dense fixed-N projector.
      rayleigh(:) = energies(:)
      CALL qs_ot_fixed_n_energy_hessian(response, response_sum, hessian)
      hessian_action(:) = MATMUL(hessian, direction)

      CALL fixed_n_fermi_occupations(energies + fd_step*direction, weights, target_electrons, &
                                     temperature, maxocc, occupation, mu)
      response_plus(:) = weights(:)*occupation(:)*(maxocc - occupation(:))/(maxocc*temperature)
      weighted_residual = DOT_PRODUCT(response_plus, rayleigh - energies - fd_step*direction)
      CALL qs_ot_fixed_n_energy_gradient(rayleigh, energies + fd_step*direction, response_plus, &
                                         SUM(response_plus), weighted_residual, gradient_plus)

      CALL fixed_n_fermi_occupations(energies - fd_step*direction, weights, target_electrons, &
                                     temperature, maxocc, occupation, mu)
      response_minus(:) = weights(:)*occupation(:)*(maxocc - occupation(:))/(maxocc*temperature)
      weighted_residual = DOT_PRODUCT(response_minus, rayleigh - energies + fd_step*direction)
      CALL qs_ot_fixed_n_energy_gradient(rayleigh, energies - fd_step*direction, response_minus, &
                                         SUM(response_minus), weighted_residual, gradient_minus)

      hessian_error = MAXVAL(ABS((gradient_plus - gradient_minus)/(2.0_dp*fd_step) - &
                                 hessian_action))
      hessian_error = MAX(hessian_error, MAXVAL(ABS(hessian - TRANSPOSE(hessian))))
      hessian_error = MAX(hessian_error, MAXVAL(ABS(MATMUL(hessian, SPREAD(1.0_dp, 1, nbands)))))

      IF (gradient_error > 2.0E-8_dp .OR. hessian_error > 2.0E-7_dp) nfail = nfail + 1
      IF (mynode == 0) THEN
         WRITE (io_unit, '(A,2(1X,ES13.6))') "fixed-N Mermin gradient/Hessian errors:", &
            gradient_error, hessian_error
      END IF

   END SUBROUTINE test_fixed_n_mermin_energy

! **************************************************************************************************
!> \brief Check the complex fixed-N spectral occupation-projector derivative by finite differences.
!> \param mynode MPI rank
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_fixed_n_projector_frechet(mynode, nfail)
      INTEGER, INTENT(IN)                                :: mynode
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nbands = 3
      REAL(KIND=dp), PARAMETER :: fd_step = 1.0E-6_dp, kpoint_weight = 0.75_dp, maxocc = 2.0_dp, &
         target_electrons = 2.10_dp, temperature = 0.09_dp

      COMPLEX(KIND=dp), DIMENSION(nbands, nbands)        :: dchc, derivative, derivative_fd, &
                                                            hamiltonian, projector_minus, &
                                                            projector_plus, vectors
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: error, mu, trace_error
      REAL(KIND=dp), DIMENSION(nbands)                   :: eigenvalues, occupation, response, &
                                                            weights

      hamiltonian(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      hamiltonian(1, 1) = CMPLX(-0.22_dp, 0.0_dp, KIND=dp)
      hamiltonian(2, 2) = CMPLX(0.03_dp, 0.0_dp, KIND=dp)
      hamiltonian(3, 3) = CMPLX(0.28_dp, 0.0_dp, KIND=dp)
      hamiltonian(1, 2) = CMPLX(0.04_dp, -0.03_dp, KIND=dp)
      hamiltonian(2, 1) = CONJG(hamiltonian(1, 2))
      hamiltonian(2, 3) = CMPLX(-0.02_dp, 0.05_dp, KIND=dp)
      hamiltonian(3, 2) = CONJG(hamiltonian(2, 3))

      dchc(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      dchc(1, 1) = CMPLX(0.13_dp, 0.0_dp, KIND=dp)
      dchc(2, 2) = CMPLX(-0.08_dp, 0.0_dp, KIND=dp)
      dchc(3, 3) = CMPLX(0.05_dp, 0.0_dp, KIND=dp)
      dchc(1, 2) = CMPLX(0.17_dp, -0.11_dp, KIND=dp)
      dchc(2, 1) = CONJG(dchc(1, 2))
      dchc(2, 3) = CMPLX(0.12_dp, 0.09_dp, KIND=dp)
      dchc(3, 2) = CONJG(dchc(2, 3))

      CALL diag_complex(hamiltonian, vectors, eigenvalues)
      weights(:) = kpoint_weight
      CALL fixed_n_fermi_occupations(eigenvalues, weights, target_electrons, temperature, &
                                     maxocc, occupation, mu)
      response(:) = weights(:)*occupation(:)*(maxocc - occupation(:))/(maxocc*temperature)
      CALL qs_ot_fixed_n_projector_frechet(hamiltonian, dchc, occupation, kpoint_weight, &
                                           response, SUM(response), derivative)

      CALL fixed_n_spectral_projector(hamiltonian + fd_step*dchc, weights, target_electrons, &
                                      temperature, maxocc, projector_plus)
      CALL fixed_n_spectral_projector(hamiltonian - fd_step*dchc, weights, target_electrons, &
                                      temperature, maxocc, projector_minus)
      derivative_fd(:, :) = (projector_plus - projector_minus)/(2.0_dp*fd_step)
      error = MAXVAL(ABS(derivative - derivative_fd))
      trace_error = 0.0_dp
      DO i = 1, nbands
         trace_error = trace_error + REAL(derivative(i, i), KIND=dp)
      END DO
      error = MAX(error, ABS(trace_error))

      IF (error > 3.0E-7_dp) nfail = nfail + 1
      IF (mynode == 0) WRITE (io_unit, '(A,1X,ES13.6)') &
         "fixed-N complex occupation-projector error:", error

   END SUBROUTINE test_fixed_n_projector_frechet

! **************************************************************************************************
!> \brief Check the finite complex REF rotation Hessian against the scalar band energy.
!> \param mynode MPI rank
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_finite_rotation_response(mynode, nfail)
      INTEGER, INTENT(IN)                                :: mynode
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nbands = 3, nrotation = 6
      REAL(KIND=dp), PARAMETER                           :: fd_step = 2.0E-4_dp, &
                                                            kpoint_weight = 0.625_dp

      COMPLEX(KIND=dp), DIMENSION(nbands, nbands)        :: base_hamiltonian, chc, &
                                                            direction_generator, generator, &
                                                            rotation, rotation_minus, rotation_plus
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: energy, energy_minus, energy_plus, &
                                                            gradient_error, hessian_error, &
                                                            rayleigh_error, symmetry_error
      REAL(KIND=dp), DIMENSION(nbands)                   :: occupation, rayleigh_minus, rayleigh_plus
      REAL(KIND=dp), DIMENSION(nbands, nrotation)        :: rayleigh_response
      REAL(KIND=dp), DIMENSION(nrotation)                :: direction, gradient
      REAL(KIND=dp), DIMENSION(nrotation, nrotation)     :: hessian

      base_hamiltonian(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base_hamiltonian(1, 1) = CMPLX(-0.31_dp, 0.0_dp, KIND=dp)
      base_hamiltonian(2, 2) = CMPLX(0.07_dp, 0.0_dp, KIND=dp)
      base_hamiltonian(3, 3) = CMPLX(0.42_dp, 0.0_dp, KIND=dp)
      base_hamiltonian(1, 2) = CMPLX(0.09_dp, -0.06_dp, KIND=dp)
      base_hamiltonian(1, 3) = CMPLX(-0.04_dp, 0.08_dp, KIND=dp)
      base_hamiltonian(2, 3) = CMPLX(0.11_dp, 0.05_dp, KIND=dp)
      DO i = 1, nbands
         base_hamiltonian(i, 1:i - 1) = CONJG(base_hamiltonian(1:i - 1, i))
      END DO

      generator(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      generator(1, 2) = CMPLX(0.17_dp, -0.08_dp, KIND=dp)
      generator(1, 3) = CMPLX(-0.09_dp, 0.04_dp, KIND=dp)
      generator(2, 3) = CMPLX(0.12_dp, 0.07_dp, KIND=dp)
      DO i = 1, nbands
         generator(i, 1:i - 1) = -CONJG(generator(1:i - 1, i))
      END DO

      direction(:) = [0.14_dp, -0.07_dp, -0.11_dp, 0.09_dp, 0.05_dp, 0.13_dp]
      direction_generator(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      direction_generator(1, 2) = CMPLX(direction(1), direction(2), KIND=dp)
      direction_generator(2, 1) = -CONJG(direction_generator(1, 2))
      direction_generator(1, 3) = CMPLX(direction(3), direction(4), KIND=dp)
      direction_generator(3, 1) = -CONJG(direction_generator(1, 3))
      direction_generator(2, 3) = CMPLX(direction(5), direction(6), KIND=dp)
      direction_generator(3, 2) = -CONJG(direction_generator(2, 3))
      occupation(:) = [1.73_dp, 0.88_dp, 0.19_dp]

      rotation(:, :) = dense_antihermitian_exp(generator)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(rotation)), &
                         MATMUL(base_hamiltonian, rotation))
      CALL qs_ot_finite_rotation_response(chc, generator, occupation, kpoint_weight, &
                                          gradient, hessian, rayleigh_response)

      rotation_plus(:, :) = dense_antihermitian_exp(generator + fd_step*direction_generator)
      rotation_minus(:, :) = dense_antihermitian_exp(generator - fd_step*direction_generator)
      energy = finite_rotation_energy(rotation, base_hamiltonian, occupation, kpoint_weight)
      energy_plus = finite_rotation_energy(rotation_plus, base_hamiltonian, occupation, kpoint_weight)
      energy_minus = finite_rotation_energy(rotation_minus, base_hamiltonian, occupation, kpoint_weight)
      gradient_error = ABS((energy_plus - energy_minus)/(2.0_dp*fd_step) - &
                           DOT_PRODUCT(gradient, direction))
      hessian_error = ABS((energy_plus - 2.0_dp*energy + energy_minus)/fd_step**2 - &
                          DOT_PRODUCT(direction, MATMUL(hessian, direction)))
      symmetry_error = MAXVAL(ABS(hessian - TRANSPOSE(hessian)))

      DO i = 1, nbands
         rayleigh_plus(i) = REAL(DOT_PRODUCT(rotation_plus(:, i), &
                                             MATMUL(base_hamiltonian, rotation_plus(:, i))), KIND=dp)
         rayleigh_minus(i) = REAL(DOT_PRODUCT(rotation_minus(:, i), &
                                              MATMUL(base_hamiltonian, rotation_minus(:, i))), KIND=dp)
      END DO
      rayleigh_error = MAXVAL(ABS((rayleigh_plus - rayleigh_minus)/(2.0_dp*fd_step) - &
                                  MATMUL(rayleigh_response, direction)))

      IF (gradient_error > 5.0E-8_dp .OR. hessian_error > 2.0E-6_dp .OR. &
          rayleigh_error > 2.0E-7_dp .OR. symmetry_error > 1.0E-12_dp) nfail = nfail + 1
      IF (mynode == 0) THEN
         WRITE (io_unit, '(A,4(1X,ES13.6))') "finite complex rotation response errors:", &
            gradient_error, hessian_error, rayleigh_error, symmetry_error
      END IF

   END SUBROUTINE test_finite_rotation_response

! **************************************************************************************************
!> \brief Check the coupled fixed-N rotation/energy Hessian and its Schur decomposition.
!> \param mynode MPI rank
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_fixed_n_rotation_schur(mynode, nfail)
      INTEGER, INTENT(IN)                                :: mynode
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nbands = 3, nrotation = 6
      REAL(KIND=dp), PARAMETER                           :: fd_step = 1.0E-4_dp, &
                                                            kpoint_weight = 0.625_dp, &
                                                            maxocc = 2.0_dp, temperature = 0.08_dp

      COMPLEX(KIND=dp), DIMENSION(nbands, nbands)        :: base_hamiltonian, chc, &
                                                            direction_generator, generator, &
                                                            rotation, rotation_minus, rotation_plus
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: curvature, curvature_error, &
                                                            energy_minus, energy_plus, &
                                                            energy_zero, mu, response_sum, &
                                                            schur_error
      REAL(KIND=dp), DIMENSION(nbands) :: energy_coordinate, energy_direction, energy_gradient, &
         energy_minus_coordinate, energy_plus_coordinate, occupation, rayleigh, rayleigh_minus, &
         rayleigh_plus, response, weights
      REAL(KIND=dp), DIMENSION(nbands, nbands)           :: energy_hessian
      REAL(KIND=dp), DIMENSION(nbands, nrotation)        :: rayleigh_response
      REAL(KIND=dp), DIMENSION(nrotation)                :: coupling, rotation_direction, &
                                                            rotation_gradient, schur_rhs
      REAL(KIND=dp), DIMENSION(nrotation, nrotation)     :: rotation_hessian, schur, schur_block

      base_hamiltonian(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      base_hamiltonian(1, 1) = CMPLX(-0.31_dp, 0.0_dp, KIND=dp)
      base_hamiltonian(2, 2) = CMPLX(0.07_dp, 0.0_dp, KIND=dp)
      base_hamiltonian(3, 3) = CMPLX(0.42_dp, 0.0_dp, KIND=dp)
      base_hamiltonian(1, 2) = CMPLX(0.09_dp, -0.06_dp, KIND=dp)
      base_hamiltonian(1, 3) = CMPLX(-0.04_dp, 0.08_dp, KIND=dp)
      base_hamiltonian(2, 3) = CMPLX(0.11_dp, 0.05_dp, KIND=dp)
      DO i = 1, nbands
         base_hamiltonian(i, 1:i - 1) = CONJG(base_hamiltonian(1:i - 1, i))
      END DO

      generator(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      generator(1, 2) = CMPLX(0.17_dp, -0.08_dp, KIND=dp)
      generator(1, 3) = CMPLX(-0.09_dp, 0.04_dp, KIND=dp)
      generator(2, 3) = CMPLX(0.12_dp, 0.07_dp, KIND=dp)
      DO i = 1, nbands
         generator(i, 1:i - 1) = -CONJG(generator(1:i - 1, i))
      END DO
      rotation_direction(:) = [0.14_dp, -0.07_dp, -0.11_dp, 0.09_dp, 0.05_dp, 0.13_dp]
      energy_direction(:) = [0.08_dp, -0.11_dp, 0.03_dp]
      direction_generator(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      direction_generator(1, 2) = CMPLX(rotation_direction(1), rotation_direction(2), KIND=dp)
      direction_generator(2, 1) = -CONJG(direction_generator(1, 2))
      direction_generator(1, 3) = CMPLX(rotation_direction(3), rotation_direction(4), KIND=dp)
      direction_generator(3, 1) = -CONJG(direction_generator(1, 3))
      direction_generator(2, 3) = CMPLX(rotation_direction(5), rotation_direction(6), KIND=dp)
      direction_generator(3, 2) = -CONJG(direction_generator(2, 3))

      rotation(:, :) = dense_antihermitian_exp(generator)
      chc(:, :) = MATMUL(CONJG(TRANSPOSE(rotation)), MATMUL(base_hamiltonian, rotation))
      DO i = 1, nbands
         rayleigh(i) = REAL(chc(i, i), KIND=dp)
      END DO
      energy_coordinate(:) = rayleigh(:)
      weights(:) = kpoint_weight
      CALL fixed_n_fermi_occupations(energy_coordinate, weights, 1.7_dp, temperature, &
                                     maxocc, occupation, mu)
      response(:) = weights*occupation*(maxocc - occupation)/(maxocc*temperature)
      response_sum = SUM(response)
      CALL qs_ot_finite_rotation_response(chc, generator, occupation, kpoint_weight, &
                                          rotation_gradient, rotation_hessian, rayleigh_response)
      CALL qs_ot_fixed_n_energy_gradient(rayleigh, energy_coordinate, response, response_sum, &
                                         0.0_dp, energy_gradient)
      CALL qs_ot_fixed_n_energy_hessian(response, response_sum, energy_hessian)
      CALL qs_ot_fixed_n_schur_block( &
         rotation_hessian, rayleigh_response, response, rotation_gradient, energy_gradient, &
         schur_block, coupling, schur_rhs)
      schur(:, :) = schur_block + SPREAD(coupling, DIM=2, NCOPIES=nrotation)* &
                    SPREAD(coupling, DIM=1, NCOPIES=nrotation)/response_sum
      schur_error = MAXVAL(ABS(schur - rotation_hessian + &
                               MATMUL(TRANSPOSE(rayleigh_response), &
                                      MATMUL(energy_hessian, rayleigh_response))))

      rotation_plus(:, :) = dense_antihermitian_exp(generator + fd_step*direction_generator)
      rotation_minus(:, :) = dense_antihermitian_exp(generator - fd_step*direction_generator)
      DO i = 1, nbands
         rayleigh_plus(i) = REAL(DOT_PRODUCT(rotation_plus(:, i), &
                                             MATMUL(base_hamiltonian, rotation_plus(:, i))), KIND=dp)
         rayleigh_minus(i) = REAL(DOT_PRODUCT(rotation_minus(:, i), &
                                              MATMUL(base_hamiltonian, rotation_minus(:, i))), KIND=dp)
      END DO
      energy_plus_coordinate(:) = energy_coordinate + fd_step*energy_direction
      energy_minus_coordinate(:) = energy_coordinate - fd_step*energy_direction
      energy_zero = fixed_n_mermin_value(rayleigh, energy_coordinate, weights, 1.7_dp, &
                                         temperature, maxocc)
      energy_plus = fixed_n_mermin_value(rayleigh_plus, energy_plus_coordinate, weights, 1.7_dp, &
                                         temperature, maxocc)
      energy_minus = fixed_n_mermin_value(rayleigh_minus, energy_minus_coordinate, weights, 1.7_dp, &
                                          temperature, maxocc)
      curvature = DOT_PRODUCT(rotation_direction, MATMUL(rotation_hessian, rotation_direction)) - &
                  2.0_dp*DOT_PRODUCT(energy_direction, &
                                     MATMUL(energy_hessian, &
                                            MATMUL(rayleigh_response, rotation_direction))) + &
                  DOT_PRODUCT(energy_direction, MATMUL(energy_hessian, energy_direction))
      curvature_error = ABS((energy_plus - 2.0_dp*energy_zero + energy_minus)/fd_step**2 - curvature)

      IF (schur_error > 2.0E-11_dp .OR. curvature_error > 3.0E-6_dp) nfail = nfail + 1
      IF (mynode == 0) THEN
         WRITE (io_unit, '(A,2(1X,ES13.6))') "fixed-N rotation Schur errors:", &
            schur_error, curvature_error
      END IF

   END SUBROUTINE test_fixed_n_rotation_schur

! **************************************************************************************************
!> \brief Check that positive and negative response modes retain their spectral magnitude.
!> \param mynode MPI rank
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_symmetric_abs_solve(mynode, nfail)
      INTEGER, INTENT(IN)                                :: mynode
      INTEGER, INTENT(INOUT)                             :: nfail

      LOGICAL                                            :: valid
      REAL(KIND=dp)                                      :: error
      REAL(KIND=dp), DIMENSION(3)                        :: eigenvalues
      REAL(KIND=dp), DIMENSION(3, 2)                     :: expected, rhs, solution
      REAL(KIND=dp), DIMENSION(3, 3)                     :: absolute_matrix, matrix, vectors

      vectors(:, 1) = [0.8_dp, 0.6_dp, 0.0_dp]
      vectors(:, 2) = [-0.6_dp, 0.8_dp, 0.0_dp]
      vectors(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]
      eigenvalues(:) = [-0.4_dp, 1.7_dp, 0.25_dp]
      matrix(:, :) = MATMUL(vectors, MATMUL( &
                            RESHAPE([eigenvalues(1), 0.0_dp, 0.0_dp, 0.0_dp, eigenvalues(2), 0.0_dp, &
                                     0.0_dp, 0.0_dp, eigenvalues(3)], [3, 3]), TRANSPOSE(vectors)))
      absolute_matrix(:, :) = MATMUL(vectors, MATMUL( &
                                     RESHAPE([ABS(eigenvalues(1)), 0.0_dp, 0.0_dp, 0.0_dp, ABS(eigenvalues(2)), 0.0_dp, &
                                              0.0_dp, 0.0_dp, ABS(eigenvalues(3))], [3, 3]), TRANSPOSE(vectors)))
      rhs(:, 1) = [0.3_dp, -0.5_dp, 0.7_dp]
      rhs(:, 2) = [-0.2_dp, 0.4_dp, 0.1_dp]
      CALL qs_ot_symmetric_abs_solve(matrix, rhs, solution, valid)
      expected(:, :) = MATMUL(absolute_matrix, solution)
      error = MAXVAL(ABS(expected - rhs))
      IF (.NOT. valid .OR. error > 2.0E-13_dp) nfail = nfail + 1

      eigenvalues(:) = [-0.4_dp, 1.7_dp, 1.0E-12_dp]
      matrix(:, :) = MATMUL(vectors, MATMUL( &
                            RESHAPE([eigenvalues(1), 0.0_dp, 0.0_dp, 0.0_dp, eigenvalues(2), 0.0_dp, &
                                     0.0_dp, 0.0_dp, eigenvalues(3)], [3, 3]), TRANSPOSE(vectors)))
      rhs(:, 1) = [0.3_dp, -0.5_dp, 0.7_dp]
      rhs(:, 2) = [-0.2_dp, 0.4_dp, 0.1_dp]
      CALL qs_ot_symmetric_abs_solve(matrix, rhs, solution, valid)
      expected(:, :) = rhs(:, :)
      expected(3, :) = 0.0_dp
      absolute_matrix(:, :) = MATMUL(vectors(:, 1:2), MATMUL( &
                                     RESHAPE([ABS(eigenvalues(1)), 0.0_dp, 0.0_dp, ABS(eigenvalues(2))], [2, 2]), &
                                     TRANSPOSE(vectors(:, 1:2))))
      error = MAX(error, MAXVAL(ABS(MATMUL(absolute_matrix, solution) - expected)))
      IF (.NOT. valid .OR. error > 2.0E-13_dp .OR. MAXVAL(ABS(solution(3, :))) > 2.0E-13_dp) THEN
         nfail = nfail + 1
      END IF
      IF (mynode == 0) WRITE (io_unit, '(A,1X,ES13.6)') "symmetric absolute solve error:", error

   END SUBROUTINE test_symmetric_abs_solve

! **************************************************************************************************
!> \brief Check that the projected response update preserves the baseline and rejects indefiniteness.
!> \param mynode MPI rank
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_projected_response_update(mynode, nfail)
      INTEGER, INTENT(IN)                                :: mynode
      INTEGER, INTENT(INOUT)                             :: nfail

      LOGICAL                                            :: valid
      REAL(KIND=dp)                                      :: error, preservation_error
      REAL(KIND=dp), DIMENSION(2)                        :: coefficients, gradient, residual
      REAL(KIND=dp), DIMENSION(2, 2)                     :: correction, reference, total

      reference(:, :) = RESHAPE([4.0_dp, 0.7_dp, 0.7_dp, 2.5_dp], [2, 2])
      correction(:, :) = 0.0_dp
      CALL qs_ot_projected_response_update(reference, correction, coefficients, valid)
      preservation_error = MAXVAL(ABS(coefficients - [1.0_dp, 0.0_dp]))
      IF (.NOT. valid .OR. preservation_error > 1.0E-12_dp) nfail = nfail + 1

      correction(:, :) = RESHAPE([0.8_dp, -0.3_dp, -0.3_dp, 0.4_dp], [2, 2])
      gradient(:) = [1.2_dp, -0.4_dp]
      CALL qs_ot_projected_response_update( &
         reference, correction, coefficients, valid, projected_gradient=gradient)
      total(:, :) = reference + correction
      residual(:) = MATMUL(total, coefficients) - gradient
      error = MAXVAL(ABS(residual))
      IF (.NOT. valid .OR. error > 1.0E-12_dp) nfail = nfail + 1

      correction(:, :) = RESHAPE([-8.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2])
      CALL qs_ot_projected_response_update(reference, correction, coefficients, valid)
      IF (valid) nfail = nfail + 1
      IF (mynode == 0) WRITE (io_unit, '(A,2(1X,ES14.6))') &
         "projected response preservation/residual errors:", preservation_error, error

   END SUBROUTINE test_projected_response_update

! **************************************************************************************************
!> \brief Check exact signed SR1 secants and rejection of an unresolved orthogonal residual.
!> \param mynode MPI rank
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_symmetric_sr1_update(mynode, nfail)
      INTEGER, INTENT(IN)                                :: mynode
      INTEGER, INTENT(INOUT)                             :: nfail

      LOGICAL                                            :: valid
      REAL(KIND=dp)                                      :: rejected_error, secant_error, &
                                                            symmetry_error
      REAL(KIND=dp), DIMENSION(3)                        :: orthogonal_residual, response, step
      REAL(KIND=dp), DIMENSION(3, 3)                     :: matrix, updated

      matrix(:, :) = 0.0_dp
      matrix(1, 1) = 1.0_dp
      matrix(2, 2) = 2.0_dp
      matrix(3, 3) = 3.0_dp
      step(:) = [0.3_dp, -0.4_dp, 0.2_dp]
      response(:) = [0.0_dp, -0.3_dp, 0.8_dp]
      CALL qs_ot_symmetric_sr1_update(matrix, step, response, updated, valid)
      secant_error = MAXVAL(ABS(MATMUL(updated, step) - response))
      symmetry_error = MAXVAL(ABS(updated - TRANSPOSE(updated)))
      IF (.NOT. valid .OR. secant_error > 2.0E-14_dp .OR. symmetry_error > 2.0E-14_dp) THEN
         nfail = nfail + 1
      END IF

      orthogonal_residual(:) = [0.4_dp, 0.3_dp, 0.0_dp]
      response(:) = MATMUL(matrix, step) + orthogonal_residual
      CALL qs_ot_symmetric_sr1_update(matrix, step, response, updated, valid)
      rejected_error = MAXVAL(ABS(updated - matrix))
      IF (valid .OR. rejected_error > 2.0E-14_dp) nfail = nfail + 1
      IF (mynode == 0) WRITE (io_unit, '(A,3(1X,ES13.6))') &
         "symmetric SR1 secant/symmetry/rejection errors:", &
         secant_error, symmetry_error, rejected_error

   END SUBROUTINE test_symmetric_sr1_update

! **************************************************************************************************
!> \brief Check the complex gauge-invariant, sign-aware density-space response secant.
!> \param mynode MPI rank
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_density_secant_hessian(mynode, nfail)
      INTEGER, INTENT(IN)                                :: mynode
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: n = 3, nmode = 2
      REAL(KIND=dp), PARAMETER                           :: secant_position = 0.4_dp

      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: density_step, hamiltonian_step, unitary
      COMPLEX(KIND=dp), DIMENSION(n, n, nmode)           :: density_modes, transformed_modes
      LOGICAL                                            :: calibrated_valid, rejected_valid, valid
      REAL(KIND=dp) :: calibrated_error, density_norm, direct_error, gauge_error, overlap_error, &
         rejected_error, response_work, sign_error
      REAL(KIND=dp), DIMENSION(nmode)                    :: density_overlap, response_overlap
      REAL(KIND=dp), DIMENSION(nmode, nmode)             :: calibrated, correction, &
                                                            correction_gauge, correction_overlap, &
                                                            expected, rejected

      density_step(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      density_step(1, 1) = CMPLX(0.30_dp, 0.0_dp, KIND=dp)
      density_step(2, 2) = CMPLX(-0.18_dp, 0.0_dp, KIND=dp)
      density_step(3, 3) = CMPLX(-0.12_dp, 0.0_dp, KIND=dp)
      density_step(1, 2) = CMPLX(0.07_dp, -0.11_dp, KIND=dp)
      density_step(2, 1) = CONJG(density_step(1, 2))
      density_step(1, 3) = CMPLX(-0.04_dp, 0.05_dp, KIND=dp)
      density_step(3, 1) = CONJG(density_step(1, 3))

      hamiltonian_step(:, :) = -0.45_dp*density_step
      hamiltonian_step(2, 3) = CMPLX(0.06_dp, 0.03_dp, KIND=dp)
      hamiltonian_step(3, 2) = CONJG(hamiltonian_step(2, 3))
      density_modes(:, :, 1) = density_step
      density_modes(:, :, 2) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      density_modes(1, 1, 2) = CMPLX(-0.2_dp, 0.0_dp, KIND=dp)
      density_modes(2, 2, 2) = CMPLX(0.1_dp, 0.0_dp, KIND=dp)
      density_modes(3, 3, 2) = CMPLX(0.1_dp, 0.0_dp, KIND=dp)
      density_modes(1, 3, 2) = CMPLX(0.03_dp, 0.08_dp, KIND=dp)
      density_modes(3, 1, 2) = CONJG(density_modes(1, 3, 2))

      CALL qs_ot_density_secant_hessian( &
         density_step, hamiltonian_step, density_modes, correction, valid, &
         density_norm, response_work)
      density_overlap(1) = density_norm
      density_overlap(2) = SUM(REAL(CONJG(density_step)*density_modes(:, :, 2), KIND=dp))
      response_overlap(1) = response_work
      response_overlap(2) = SUM(REAL(CONJG(hamiltonian_step)*density_modes(:, :, 2), KIND=dp))
      CALL qs_ot_density_secant_projected_hessian( &
         density_norm, response_work, density_overlap, response_overlap, correction_overlap, valid)

      expected(1, 1) = response_work
      expected(1, 2) = response_overlap(2)
      expected(2, 1) = response_overlap(2)
      expected(2, 2) = (2.0_dp*density_overlap(2)*response_overlap(2) - &
                        response_work*density_overlap(2)**2/density_norm)/density_norm
      direct_error = MAXVAL(ABS(correction - expected))
      overlap_error = MAXVAL(ABS(correction_overlap - correction))
      sign_error = MAX(0.0_dp, correction(1, 1))

      CALL qs_ot_density_secant_projected_hessian( &
         density_norm, response_work, density_overlap, response_overlap, calibrated, &
         calibrated_valid, secant_mode=2, secant_position=secant_position)
      expected(1, 1) = response_work
      expected(1, 2) = response_work/secant_position
      expected(2, 1) = expected(1, 2)
      expected(2, 2) = response_work/secant_position**2
      calibrated_error = MAXVAL(ABS(calibrated - expected))
      CALL qs_ot_density_secant_projected_hessian( &
         density_norm, response_work, density_overlap, response_overlap, rejected, &
         rejected_valid, secant_mode=2, secant_position=0.0_dp)
      rejected_error = MAXVAL(ABS(rejected))

      unitary(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      unitary(1, 1) = CMPLX(SQRT(0.5_dp), 0.0_dp, KIND=dp)
      unitary(2, 1) = CMPLX(0.0_dp, SQRT(0.5_dp), KIND=dp)
      unitary(1, 2) = CMPLX(0.0_dp, SQRT(0.5_dp), KIND=dp)
      unitary(2, 2) = CMPLX(SQRT(0.5_dp), 0.0_dp, KIND=dp)
      unitary(3, 3) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      density_step = MATMUL(CONJG(TRANSPOSE(unitary)), MATMUL(density_step, unitary))
      hamiltonian_step = MATMUL(CONJG(TRANSPOSE(unitary)), MATMUL(hamiltonian_step, unitary))
      transformed_modes(:, :, 1) = &
         MATMUL(CONJG(TRANSPOSE(unitary)), MATMUL(density_modes(:, :, 1), unitary))
      transformed_modes(:, :, 2) = &
         MATMUL(CONJG(TRANSPOSE(unitary)), MATMUL(density_modes(:, :, 2), unitary))
      CALL qs_ot_density_secant_hessian( &
         density_step, hamiltonian_step, transformed_modes, correction_gauge, valid)
      gauge_error = MAXVAL(ABS(correction_gauge - correction))

      IF (.NOT. valid .OR. .NOT. calibrated_valid .OR. rejected_valid .OR. &
          direct_error > 2.0E-14_dp .OR. overlap_error > 2.0E-14_dp .OR. &
          gauge_error > 2.0E-14_dp .OR. sign_error > 2.0E-14_dp .OR. &
          calibrated_error > 2.0E-14_dp .OR. rejected_error > 2.0E-14_dp) nfail = nfail + 1
      IF (mynode == 0) WRITE (io_unit, '(A,6(1X,ES13.6))') &
         "density secant direct/overlap/gauge/sign/calibrated/rejected errors:", &
         direct_error, overlap_error, gauge_error, sign_error, calibrated_error, rejected_error

   END SUBROUTINE test_density_secant_hessian

! **************************************************************************************************
!> \brief Check the coupled finite-chart density tangent against a central finite difference.
!> \param mynode MPI rank
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_density_tangent(mynode, nfail)
      INTEGER, INTENT(IN)                                :: mynode
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nbands = 3, nrotation = 6
      REAL(KIND=dp), PARAMETER                           :: fd_step = 2.0E-6_dp, &
                                                            kpoint_weight = 0.625_dp

      COMPLEX(KIND=dp), DIMENSION(nbands, nbands) :: density_minus, density_plus, density_tangent, &
         direction_generator, generator, rotation, rotation_minus, rotation_plus, weighted_rotation
      INTEGER                                            :: i, j, r
      REAL(KIND=dp)                                      :: error, trace_error
      REAL(KIND=dp), DIMENSION(nbands)                   :: occupation, weighted_occupation_step
      REAL(KIND=dp), DIMENSION(nrotation)                :: rotation_step

      generator(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      generator(1, 2) = CMPLX(0.17_dp, -0.09_dp, KIND=dp)
      generator(2, 1) = -CONJG(generator(1, 2))
      generator(1, 3) = CMPLX(-0.06_dp, 0.11_dp, KIND=dp)
      generator(3, 1) = -CONJG(generator(1, 3))
      generator(2, 3) = CMPLX(0.08_dp, 0.04_dp, KIND=dp)
      generator(3, 2) = -CONJG(generator(2, 3))
      occupation(:) = [1.72_dp, 0.83_dp, 0.31_dp]
      rotation_step(:) = [0.21_dp, -0.14_dp, -0.08_dp, 0.17_dp, 0.11_dp, 0.06_dp]
      weighted_occupation_step(:) = [-0.037_dp, 0.052_dp, -0.015_dp]

      CALL qs_ot_density_tangent( &
         generator, occupation, kpoint_weight, rotation_step, weighted_occupation_step, &
         density_tangent, difference_step=fd_step)

      direction_generator(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      r = 0
      DO i = 1, nbands - 1
         DO j = i + 1, nbands
            r = r + 1
            direction_generator(i, j) = CMPLX(rotation_step(r), 0.0_dp, KIND=dp)
            direction_generator(j, i) = -direction_generator(i, j)
            r = r + 1
            direction_generator(i, j) = direction_generator(i, j) + &
                                        CMPLX(0.0_dp, rotation_step(r), KIND=dp)
            direction_generator(j, i) = direction_generator(j, i) + &
                                        CMPLX(0.0_dp, rotation_step(r), KIND=dp)
         END DO
      END DO
      rotation = dense_antihermitian_exp(generator)
      rotation_plus = dense_antihermitian_exp(generator + fd_step*direction_generator)
      rotation_minus = dense_antihermitian_exp(generator - fd_step*direction_generator)
      weighted_rotation(:, :) = rotation_plus
      DO j = 1, nbands
         weighted_rotation(:, j) = &
            (kpoint_weight*occupation(j) + fd_step*weighted_occupation_step(j))* &
            weighted_rotation(:, j)
      END DO
      density_plus(:, :) = MATMUL(weighted_rotation, CONJG(TRANSPOSE(rotation_plus)))
      weighted_rotation(:, :) = rotation_minus
      DO j = 1, nbands
         weighted_rotation(:, j) = &
            (kpoint_weight*occupation(j) - fd_step*weighted_occupation_step(j))* &
            weighted_rotation(:, j)
      END DO
      density_minus(:, :) = MATMUL(weighted_rotation, CONJG(TRANSPOSE(rotation_minus)))
      density_plus(:, :) = MATMUL(CONJG(TRANSPOSE(rotation)), &
                                  MATMUL((density_plus - density_minus)/(2.0_dp*fd_step), rotation))
      error = MAXVAL(ABS(density_tangent - density_plus))
      trace_error = ABS(REAL(SUM([(density_tangent(i, i), i=1, nbands)]), KIND=dp) - &
                        SUM(weighted_occupation_step))

      IF (error > 3.0E-9_dp .OR. trace_error > 3.0E-10_dp) nfail = nfail + 1
      IF (mynode == 0) WRITE (io_unit, '(A,2(1X,ES13.6))') &
         "finite complex density tangent errors:", error, trace_error

   END SUBROUTINE test_density_tangent

! **************************************************************************************************
!> \brief Check density-secants when accepted endpoint orbital subspaces differ.
!> \param mynode MPI rank
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_density_secant_moving_subspace(mynode, nfail)
      INTEGER, INTENT(IN)                                :: mynode
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nao = 4, nbands = 3, nmode = 2
      REAL(KIND=dp), PARAMETER                           :: weight = 0.625_dp

      COMPLEX(KIND=dp), DIMENSION(nao, nao)              :: density_step, hamiltonian_step
      COMPLEX(KIND=dp), DIMENSION(nao, nao, nmode)       :: density_modes_ao
      COMPLEX(KIND=dp), DIMENSION(nao, nbands)           :: c0, c1
      COMPLEX(KIND=dp), DIMENSION(nbands, nbands)        :: h0, h1, overlap
      COMPLEX(KIND=dp), DIMENSION(nbands, nbands, nmode) :: density_modes
      INTEGER                                            :: i
      LOGICAL                                            :: valid
      REAL(KIND=dp)                                      :: density_norm, density_norm_ref, error, &
                                                            response_work, response_work_ref
      REAL(KIND=dp), DIMENSION(nbands)                   :: occupation0, occupation1
      REAL(KIND=dp), DIMENSION(nmode)                    :: density_overlap, density_overlap_ref, &
                                                            response_overlap, response_overlap_ref

      c0(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      c0(1, 1) = 1.0_dp
      c0(2, 2) = 1.0_dp
      c0(3, 3) = 1.0_dp
      c1(:, :) = c0
      c1(1, 1) = COS(0.31_dp)
      c1(4, 1) = SIN(0.31_dp)
      c1(2, 2) = COS(0.23_dp)
      c1(3, 2) = CMPLX(0.0_dp, SIN(0.23_dp), KIND=dp)
      c1(2, 3) = CMPLX(0.0_dp, SIN(0.23_dp), KIND=dp)
      c1(3, 3) = COS(0.23_dp)
      occupation0(:) = [1.73_dp, 0.82_dp, 0.82_dp]
      occupation1(:) = [1.61_dp, 0.91_dp, 0.76_dp]

      hamiltonian_step(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      hamiltonian_step(1, 1) = 0.21_dp
      hamiltonian_step(2, 2) = -0.17_dp
      hamiltonian_step(3, 3) = 0.09_dp
      hamiltonian_step(4, 4) = -0.04_dp
      hamiltonian_step(1, 4) = CMPLX(0.07_dp, -0.03_dp, KIND=dp)
      hamiltonian_step(4, 1) = CONJG(hamiltonian_step(1, 4))
      hamiltonian_step(2, 3) = CMPLX(-0.05_dp, 0.08_dp, KIND=dp)
      hamiltonian_step(3, 2) = CONJG(hamiltonian_step(2, 3))
      density_modes(:, :, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      density_modes(1, 1, 1) = 0.13_dp
      density_modes(2, 2, 1) = -0.08_dp
      density_modes(3, 3, 1) = -0.05_dp
      density_modes(1, 2, 2) = CMPLX(0.04_dp, -0.06_dp, KIND=dp)
      density_modes(2, 1, 2) = CONJG(density_modes(1, 2, 2))

      overlap(:, :) = MATMUL(CONJG(TRANSPOSE(c0)), c1)
      h0(:, :) = MATMUL(CONJG(TRANSPOSE(c0)), MATMUL(hamiltonian_step, c0))
      h1(:, :) = MATMUL(CONJG(TRANSPOSE(c1)), MATMUL(hamiltonian_step, c1))
      CALL qs_ot_density_secant_orbital_overlaps( &
         overlap, occupation0, occupation1, h0, h1, density_modes, weight, density_norm, &
         response_work, density_overlap, response_overlap, valid)

      density_step(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i = 1, nbands
         density_step = density_step + weight*occupation1(i)* &
                        SPREAD(c1(:, i), DIM=2, NCOPIES=nao)* &
                        SPREAD(CONJG(c1(:, i)), DIM=1, NCOPIES=nao) - &
                        weight*occupation0(i)* &
                        SPREAD(c0(:, i), DIM=2, NCOPIES=nao)* &
                        SPREAD(CONJG(c0(:, i)), DIM=1, NCOPIES=nao)
      END DO
      DO i = 1, nmode
         density_modes_ao(:, :, i) = MATMUL(c1, &
                                            MATMUL(density_modes(:, :, i), &
                                                   CONJG(TRANSPOSE(c1))))
         density_overlap_ref(i) = REAL( &
                                  SUM(CONJG(density_step)*density_modes_ao(:, :, i)), KIND=dp)/weight
         response_overlap_ref(i) = &
            REAL(SUM(CONJG(hamiltonian_step)*density_modes_ao(:, :, i)), KIND=dp)
      END DO
      density_norm_ref = REAL(SUM(CONJG(density_step)*density_step), KIND=dp)/weight
      response_work_ref = REAL(SUM(CONJG(density_step)*hamiltonian_step), KIND=dp)
      error = MAX(ABS(density_norm - density_norm_ref), ABS(response_work - response_work_ref))
      error = MAX(error, MAXVAL(ABS(density_overlap - density_overlap_ref)))
      error = MAX(error, MAXVAL(ABS(response_overlap - response_overlap_ref)))

      IF (.NOT. valid .OR. error > 3.0E-13_dp) nfail = nfail + 1
      IF (mynode == 0) WRITE (io_unit, '(A,1X,ES13.6)') &
         "moving-subspace density secant error:", error

   END SUBROUTINE test_density_secant_moving_subspace

! **************************************************************************************************
!> \brief fixed-occupation band energy for a dense finite rotation
!> \param rotation unitary REF rotation
!> \param hamiltonian fixed Hermitian Hamiltonian in the REF basis
!> \param occupation band occupations
!> \param kpoint_weight irreducible-k-point weight
!> \return weighted band energy
! **************************************************************************************************
   FUNCTION finite_rotation_energy(rotation, hamiltonian, occupation, kpoint_weight) RESULT(energy)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: rotation, hamiltonian
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation
      REAL(KIND=dp), INTENT(IN)                          :: kpoint_weight
      REAL(KIND=dp)                                      :: energy

      INTEGER                                            :: i

      energy = 0.0_dp
      DO i = 1, SIZE(occupation)
         energy = energy + kpoint_weight*occupation(i)* &
                  REAL(DOT_PRODUCT(rotation(:, i), &
                                   MATMUL(hamiltonian, rotation(:, i))), KIND=dp)
      END DO

   END FUNCTION finite_rotation_energy

! **************************************************************************************************
!> \brief Solve fixed-N Fermi occupations by bisection.
!> \param energies band energies
!> \param weights integration weights
!> \param target_electrons requested electron number
!> \param temperature electronic temperature in energy units
!> \param maxocc maximum band occupation
!> \param occupation resulting occupations
!> \param mu resulting chemical potential
! **************************************************************************************************
   SUBROUTINE fixed_n_fermi_occupations(energies, weights, target_electrons, temperature, &
                                        maxocc, occupation, mu)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: energies, weights
      REAL(KIND=dp), INTENT(IN)                          :: target_electrons, temperature, maxocc
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: occupation
      REAL(KIND=dp), INTENT(OUT)                         :: mu

      INTEGER                                            :: i, iteration
      REAL(KIND=dp)                                      :: high, low, number

      low = MINVAL(energies) - 80.0_dp*temperature
      high = MAXVAL(energies) + 80.0_dp*temperature
      DO iteration = 1, 160
         mu = 0.5_dp*(low + high)
         number = 0.0_dp
         DO i = 1, SIZE(energies)
            occupation(i) = fermi_occupation(energies(i), mu, temperature, maxocc)
            number = number + weights(i)*occupation(i)
         END DO
         IF (number < target_electrons) THEN
            low = mu
         ELSE
            high = mu
         END IF
      END DO
      mu = 0.5_dp*(low + high)
      DO i = 1, SIZE(energies)
         occupation(i) = fermi_occupation(energies(i), mu, temperature, maxocc)
      END DO

   END SUBROUTINE fixed_n_fermi_occupations

! **************************************************************************************************
!> \brief Numerically stable Fermi occupation.
!> \param energy band energy
!> \param mu chemical potential
!> \param temperature electronic temperature in energy units
!> \param maxocc maximum band occupation
!> \return Fermi occupation
! **************************************************************************************************
   PURE FUNCTION fermi_occupation(energy, mu, temperature, maxocc) RESULT(occupation)
      REAL(KIND=dp), INTENT(IN)                          :: energy, mu, temperature, maxocc
      REAL(KIND=dp)                                      :: occupation

      REAL(KIND=dp)                                      :: x

      x = (energy - mu)/temperature
      IF (x > 40.0_dp) THEN
         occupation = maxocc*EXP(-x)
      ELSE IF (x < -40.0_dp) THEN
         occupation = maxocc*(1.0_dp - EXP(x))
      ELSE
         occupation = maxocc/(1.0_dp + EXP(x))
      END IF

   END FUNCTION fermi_occupation

! **************************************************************************************************
!> \brief Fixed-H fixed-N Mermin free energy used by the finite-difference test.
!> \param rayleigh Hamiltonian expectation values
!> \param energies auxiliary band energies
!> \param weights integration weights
!> \param target_electrons requested electron number
!> \param temperature electronic temperature in energy units
!> \param maxocc maximum band occupation
!> \return fixed-N Mermin free energy
! **************************************************************************************************
   FUNCTION fixed_n_mermin_value(rayleigh, energies, weights, target_electrons, &
                                 temperature, maxocc) RESULT(value)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rayleigh, energies, weights
      REAL(KIND=dp), INTENT(IN)                          :: target_electrons, temperature, maxocc
      REAL(KIND=dp)                                      :: value

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: mu, p
      REAL(KIND=dp), DIMENSION(SIZE(energies))           :: occupation

      CALL fixed_n_fermi_occupations(energies, weights, target_electrons, temperature, &
                                     maxocc, occupation, mu)
      value = 0.0_dp
      DO i = 1, SIZE(energies)
         p = MAX(EPSILON(1.0_dp), MIN(1.0_dp - EPSILON(1.0_dp), occupation(i)/maxocc))
         value = value + weights(i)*(occupation(i)*rayleigh(i) + &
                                     temperature*maxocc*(p*LOG(p) + (1.0_dp - p)*LOG(1.0_dp - p)))
      END DO

   END FUNCTION fixed_n_mermin_value

! **************************************************************************************************
!> \brief Canonical weighted occupation projector used by the finite-difference test.
!> \param hamiltonian Hermitian Hamiltonian
!> \param weights integration weights
!> \param target_electrons requested electron number
!> \param temperature electronic temperature in energy units
!> \param maxocc maximum band occupation
!> \param projector weighted spectral occupation projector
! **************************************************************************************************
   SUBROUTINE fixed_n_spectral_projector(hamiltonian, weights, target_electrons, temperature, &
                                         maxocc, projector)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: hamiltonian
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: weights
      REAL(KIND=dp), INTENT(IN)                          :: target_electrons, temperature, maxocc
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: projector

      COMPLEX(KIND=dp), DIMENSION(SIZE(hamiltonian, 1), &
         SIZE(hamiltonian, 1))                           :: vectors, weighted_vectors
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: mu
      REAL(KIND=dp), DIMENSION(SIZE(hamiltonian, 1))     :: eigenvalues, occupation

      CALL diag_complex(hamiltonian, vectors, eigenvalues)
      CALL fixed_n_fermi_occupations(eigenvalues, weights, target_electrons, temperature, &
                                     maxocc, occupation, mu)
      weighted_vectors(:, :) = vectors(:, :)
      DO i = 1, SIZE(eigenvalues)
         weighted_vectors(:, i) = weights(i)*occupation(i)*weighted_vectors(:, i)
      END DO
      projector(:, :) = MATMUL(weighted_vectors, CONJG(TRANSPOSE(vectors)))

   END SUBROUTINE fixed_n_spectral_projector

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
      COMPLEX(KIND=dp), DIMENSION(k, k)                  :: direction_a, generator, gradient_a, &
                                                            inverse_sqrt, rotation
      COMPLEX(KIND=dp), DIMENSION(n, k)                  :: c0, direction, gradient, hc, orbitals, &
                                                            orbitals_base, sc0, sx, x
      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: hamiltonian, overlap
      INTEGER                                            :: i
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: error, error_rot, fd_slope, &
                                                            fd_slope_rot, hstep, ortho_error, &
                                                            slope, slope_rot
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      REAL(KIND=dp), DIMENSION(k)                        :: occupation
      REAL(KIND=dp), DIMENSION(k), TARGET                :: rotation_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_type), DIMENSION(10), TARGET            :: rot_matrix
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
      DO i = 1, SIZE(rot_matrix)
         WRITE (name, '(A,I0)') 'complex_strict_rot_', i
         CALL create_one_block(rot_matrix(i), TRIM(name), dist_kk, row_k, col_k)
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

      generator(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      generator(1, 1) = CMPLX(0.0_dp, 0.16_dp, KIND=dp)
      generator(2, 2) = CMPLX(0.0_dp, -0.09_dp, KIND=dp)
      generator(1, 2) = CMPLX(0.18_dp, -0.07_dp, KIND=dp)
      generator(2, 1) = -CONJG(generator(1, 2))
      direction_a(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      direction_a(1, 1) = CMPLX(0.0_dp, -0.05_dp, KIND=dp)
      direction_a(2, 2) = CMPLX(0.0_dp, 0.04_dp, KIND=dp)
      direction_a(1, 2) = CMPLX(-0.08_dp, 0.06_dp, KIND=dp)
      direction_a(2, 1) = -CONJG(direction_a(1, 2))
      occupation = [1.0_dp, 0.34_dp]

      CALL strict_orbitals(c0, x, overlap, orbitals_base)
      rotation = dense_antihermitian_exp(generator)
      orbitals = MATMUL(orbitals_base, rotation)
      hc = 2.0_dp*MATMUL(hamiltonian, orbitals)
      DO i = 1, k
         hc(:, i) = occupation(i)*hc(:, i)
      END DO
      sx = MATMUL(overlap, x)
      sc0 = MATMUL(overlap, c0)
      CALL put_complex_pair(nk(1), nk(2), x)
      CALL put_complex_pair(nk(3), nk(4), sx)
      CALL put_complex_pair(nk(5), nk(6), hc)
      CALL put_complex_pair(nk(9), nk(10), c0)
      CALL put_complex_pair(nk(11), nk(12), sc0)
      CALL put_complex_pair(rot_matrix(3), rot_matrix(4), generator)

      strict_env%has_complex_kpoint_state = .TRUE.
      strict_env%settings%ot_algorithm = 'TOD'
      strict_env%settings%do_rotation = .TRUE.
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
      strict_env%rot_mat_u => rot_matrix(1)
      strict_env%rot_mat_u_im => rot_matrix(2)
      strict_env%rot_mat_x => rot_matrix(3)
      strict_env%rot_mat_x_im => rot_matrix(4)
      strict_env%rot_mat_dedu => rot_matrix(5)
      strict_env%rot_mat_dedu_im => rot_matrix(6)
      strict_env%rot_mat_gx => rot_matrix(7)
      strict_env%rot_mat_gx_im => rot_matrix(8)
      strict_env%rot_mat_evec_re => rot_matrix(9)
      strict_env%rot_mat_evec_im => rot_matrix(10)
      strict_env%rot_mat_evals => rotation_evals
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
      CALL get_complex_pair(rot_matrix(7), rot_matrix(8), gradient_a)

      hstep = 1.0E-6_dp
      fd_slope = (strict_rotated_energy(c0, x + hstep*direction, generator, overlap, &
                                        hamiltonian, occupation) - &
                  strict_rotated_energy(c0, x - hstep*direction, generator, overlap, &
                                        hamiltonian, occupation))/(2.0_dp*hstep)
      slope = REAL(SUM(CONJG(gradient)*direction), KIND=dp)
      error = ABS(slope - fd_slope)
      fd_slope_rot = (strict_rotated_energy(c0, x, generator + hstep*direction_a, overlap, &
                                            hamiltonian, occupation) - &
                      strict_rotated_energy(c0, x, generator - hstep*direction_a, overlap, &
                                            hamiltonian, occupation))/(2.0_dp*hstep)
      slope_rot = 0.5_dp*(SUM(REAL(gradient_a, KIND=dp)*REAL(direction_a, KIND=dp)) + &
                          SUM(AIMAG(gradient_a)*AIMAG(direction_a)))
      error_rot = ABS(slope_rot - fd_slope_rot)
      ortho_error = MAXVAL(ABS(MATMUL(CONJG(TRANSPOSE(orbitals)), MATMUL(overlap, orbitals)) - &
                               RESHAPE([CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
                                        CMPLX(0.0_dp, 0.0_dp, KIND=dp), &
                                        CMPLX(0.0_dp, 0.0_dp, KIND=dp), &
                                        CMPLX(1.0_dp, 0.0_dp, KIND=dp)], [k, k])))
      IF (error > 5.0E-8_dp .OR. error_rot > 5.0E-8_dp .OR. &
          ortho_error > 5.0E-13_dp) nfail = nfail + 1
      IF (para_env%is_source()) THEN
         WRITE (io_unit, '(A,3(1X,ES13.6))') 'complex finite-STRICT derivative: fd/predicted/error', &
            fd_slope, slope, error
         WRITE (io_unit, '(A,3(1X,ES13.6))') 'complex finite-STRICT rotation: fd/predicted/error', &
            fd_slope_rot, slope_rot, error_rot
         WRITE (io_unit, '(A,1X,ES13.6)') 'complex finite-STRICT orthonormality error', ortho_error
      END IF

      DEALLOCATE (strict_env%evals, strict_env%dum)
      DO i = SIZE(rot_matrix), 1, -1
         CALL dbcsr_release(rot_matrix(i))
      END DO
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
      COMPLEX(KIND=dp), DIMENSION(k)                     :: phases
      COMPLEX(KIND=dp), DIMENSION(n, k)                  :: c0, c0_phase, c0_rot, gradient, &
                                                            gradient_phase, output, output_phase
      COMPLEX(KIND=dp), DIMENSION(n, n) :: eigenvectors, full_all_vectors, h, inverse, &
         inverse_rot, inverse_s, inverse_t, reference_s, reference_t, s, spectral_metric, t
      INTEGER                                            :: i, j
      REAL(KIND=dp) :: full_all_descent, full_all_gauge_error, full_all_orthogonality_error, &
         gauge_error, hermitian_error, s_error, t_error
      REAL(KIND=dp), DIMENSION(k)                        :: c0_evals
      REAL(KIND=dp), DIMENSION(n)                        :: eigenvalues
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type) :: matrix_c0, matrix_c0_phase, matrix_c0_rot, matrix_gradient, &
         matrix_gradient_phase, matrix_h, matrix_output, matrix_output_phase, matrix_s, matrix_t
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct_nk, fmstruct_nn
      TYPE(preconditioner_type)                          :: full_all, full_all_phase, &
                                                            preconditioner, preconditioner_rot, &
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
      CALL cp_cfm_create(matrix_c0_phase, fmstruct_nk, name='complex preconditioner phased C')
      CALL cp_cfm_create(matrix_c0_rot, fmstruct_nk, name='complex preconditioner C U')
      CALL cp_cfm_create(matrix_gradient, fmstruct_nk, name='complex preconditioner gradient')
      CALL cp_cfm_create(matrix_gradient_phase, fmstruct_nk, &
                         name='complex preconditioner phased gradient')
      CALL cp_cfm_create(matrix_output, fmstruct_nk, name='complex preconditioner output')
      CALL cp_cfm_create(matrix_output_phase, fmstruct_nk, &
                         name='complex preconditioner phased output')
      CALL cp_cfm_set_all(matrix_h, CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      CALL cp_cfm_set_all(matrix_s, CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      CALL cp_cfm_set_all(matrix_t, CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      CALL cp_cfm_set_all(matrix_c0, CMPLX(0.0_dp, 0.0_dp, KIND=dp))
      CALL cp_cfm_set_all(matrix_c0_phase, CMPLX(0.0_dp, 0.0_dp, KIND=dp))
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
      phases = [EXP(CMPLX(0.0_dp, 0.43_dp, KIND=dp)), &
                EXP(CMPLX(0.0_dp, -0.67_dp, KIND=dp))]
      DO j = 1, k
         c0_phase(:, j) = phases(j)*c0(:, j)
      END DO
      gradient(:, 1) = [CMPLX(0.17_dp, 0.02_dp, KIND=dp), &
                        CMPLX(-0.08_dp, 0.05_dp, KIND=dp), &
                        CMPLX(0.11_dp, -0.04_dp, KIND=dp)]
      gradient(:, 2) = [CMPLX(-0.03_dp, 0.06_dp, KIND=dp), &
                        CMPLX(0.09_dp, -0.02_dp, KIND=dp), &
                        CMPLX(-0.14_dp, 0.05_dp, KIND=dp)]
      DO j = 1, k
         gradient_phase(:, j) = phases(j)*gradient(:, j)
      END DO
      c0_evals = [REAL(h(1, 1), KIND=dp), REAL(h(2, 2), KIND=dp)]

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
            CALL cp_cfm_set_element(matrix_c0_phase, i, j, c0_phase(i, j))
            CALL cp_cfm_set_element(matrix_c0_rot, i, j, c0_rot(i, j))
            CALL cp_cfm_set_element(matrix_gradient, i, j, gradient(i, j))
            CALL cp_cfm_set_element(matrix_gradient_phase, i, j, gradient_phase(i, j))
         END DO
      END DO

      CALL init_preconditioner(full_all, para_env, blacs_env)
      CALL init_preconditioner(full_all_phase, para_env, blacs_env)
      CALL init_preconditioner(preconditioner, para_env, blacs_env)
      CALL init_preconditioner(preconditioner_rot, para_env, blacs_env)
      CALL init_preconditioner(preconditioner_s, para_env, blacs_env)
      CALL init_preconditioner(preconditioner_t, para_env, blacs_env)
      CALL make_complex_full_single_inverse(preconditioner, matrix_c0, matrix_h, matrix_s, 0.20_dp)
      CALL make_complex_full_single_inverse(preconditioner_rot, matrix_c0_rot, matrix_h, matrix_s, 0.20_dp)
      CALL make_complex_full_s_inverse(preconditioner_s, matrix_s)
      CALL make_complex_full_kinetic(preconditioner_t, matrix_t, matrix_s, 0.20_dp)
      CALL make_complex_full_all(full_all, matrix_c0, matrix_h, matrix_s, c0_evals, 0.08_dp)
      CALL make_complex_full_all(full_all_phase, matrix_c0_phase, matrix_h, matrix_s, &
                                 c0_evals, 0.08_dp)
      CALL apply_preconditioner(full_all, matrix_gradient, matrix_output)
      CALL apply_preconditioner(full_all_phase, matrix_gradient_phase, matrix_output_phase)
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
      DO j = 1, k
         DO i = 1, n
            CALL cp_cfm_get_element(matrix_output, i, j, value)
            output(i, j) = value
            CALL cp_cfm_get_element(matrix_output_phase, i, j, value)
            output_phase(i, j) = value
         END DO
      END DO
      DO j = 1, n
         DO i = 1, n
            CALL cp_cfm_get_element(full_all%complex_fm, i, j, value)
            full_all_vectors(i, j) = value
         END DO
      END DO
      full_all_gauge_error = 0.0_dp
      DO j = 1, k
         full_all_gauge_error = MAX(full_all_gauge_error, &
                                    MAXVAL(ABS(output_phase(:, j) - phases(j)*output(:, j))))
      END DO
      spectral_metric = MATMUL(CONJG(TRANSPOSE(full_all_vectors)), &
                               MATMUL(s, full_all_vectors))
      full_all_orthogonality_error = MAXVAL(ABS(spectral_metric - &
                                                RESHAPE([CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
                                                         CMPLX(0.0_dp, 0.0_dp, KIND=dp), &
                                                         CMPLX(0.0_dp, 0.0_dp, KIND=dp), &
                                                         CMPLX(0.0_dp, 0.0_dp, KIND=dp), &
                                                         CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
                                                         CMPLX(0.0_dp, 0.0_dp, KIND=dp), &
                                                         CMPLX(0.0_dp, 0.0_dp, KIND=dp), &
                                                         CMPLX(0.0_dp, 0.0_dp, KIND=dp), &
                                                         CMPLX(1.0_dp, 0.0_dp, KIND=dp)], [n, n])))
      full_all_descent = REAL(SUM(CONJG(gradient)*output), KIND=dp)
      IF (gauge_error > 1.0E-11_dp .OR. hermitian_error > 1.0E-12_dp .OR. &
          MINVAL(eigenvalues) <= 0.0_dp) nfail = nfail + 1
      IF (s_error > 1.0E-11_dp .OR. t_error > 1.0E-11_dp) nfail = nfail + 1
      IF (full_all_gauge_error > 1.0E-11_dp .OR. &
          full_all_orthogonality_error > 1.0E-11_dp .OR. full_all_descent <= 0.0_dp) nfail = nfail + 1
      IF (para_env%is_source()) THEN
         WRITE (io_unit, '(A,3(1X,ES13.6))') 'complex FULL_SINGLE gauge/hermitian/min-eigenvalue', &
            gauge_error, hermitian_error, MINVAL(eigenvalues)
         WRITE (io_unit, '(A,2(1X,ES13.6))') 'complex S/KINETIC inverse errors', s_error, t_error
         WRITE (io_unit, '(A,3(1X,ES13.6))') 'complex FULL_ALL phase/orthogonality/descent', &
            full_all_gauge_error, full_all_orthogonality_error, full_all_descent
      END IF

      CALL destroy_preconditioner(full_all_phase)
      CALL destroy_preconditioner(full_all)
      CALL destroy_preconditioner(preconditioner_t)
      CALL destroy_preconditioner(preconditioner_s)
      CALL destroy_preconditioner(preconditioner_rot)
      CALL destroy_preconditioner(preconditioner)
      CALL cp_cfm_release(matrix_c0_rot)
      CALL cp_cfm_release(matrix_c0_phase)
      CALL cp_cfm_release(matrix_c0)
      CALL cp_cfm_release(matrix_output_phase)
      CALL cp_cfm_release(matrix_output)
      CALL cp_cfm_release(matrix_gradient_phase)
      CALL cp_cfm_release(matrix_gradient)
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
      IF (.NOT. qs_ot_kpoint_preconditioner_supported(ot_precond_full_all, .FALSE.)) nfail = nfail + 1
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
          ot_precond_full_all, ot_precond_solver_default, .FALSE.)) nfail = nfail + 1
      IF (qs_ot_kpoint_preconditioner_solver_supported( &
          ot_precond_full_all, ot_precond_solver_inv_chol, .FALSE.)) nfail = nfail + 1
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

! **************************************************************************************************
!> \brief Energy of finite STRICT orbitals followed by a complex occupied-space rotation.
!> \param c0 reference orbitals
!> \param coordinate tangent coordinate
!> \param generator anti-Hermitian rotation generator
!> \param overlap Hermitian overlap matrix
!> \param hamiltonian Hermitian Hamiltonian
!> \param occupation orbital occupations
!> \return energy
! **************************************************************************************************
   FUNCTION strict_rotated_energy(c0, coordinate, generator, overlap, hamiltonian, &
                                  occupation) RESULT(energy)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: c0, coordinate, generator, overlap, &
                                                            hamiltonian
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation
      REAL(KIND=dp)                                      :: energy

      COMPLEX(KIND=dp), DIMENSION(SIZE(coordinate, 1), &
         SIZE(coordinate, 2))                            :: orbitals
      INTEGER                                            :: i

      CALL strict_orbitals(c0, coordinate, overlap, orbitals)
      orbitals = MATMUL(orbitals, dense_antihermitian_exp(generator))
      energy = 0.0_dp
      DO i = 1, SIZE(occupation)
         energy = energy + occupation(i)* &
                  REAL(DOT_PRODUCT(orbitals(:, i), MATMUL(hamiltonian, orbitals(:, i))), KIND=dp)
      END DO

   END FUNCTION strict_rotated_energy

END PROGRAM qs_ot_complex_ref_unittest
