!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Focused regressions for complex k-point OT derivatives.
! **************************************************************************************************
PROGRAM qs_ot_complex_ref_unittest
   USE bibliography,                    ONLY: add_all_references
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_finalize_lib, dbcsr_get_block_p, dbcsr_init_lib, dbcsr_put_block, &
        dbcsr_release, dbcsr_reserve_blocks, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_logger_create,&
                                              cp_logger_release,&
                                              cp_logger_type,&
                                              cp_rm_default_logger
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: default_output_unit
   USE mathlib,                         ONLY: diag_complex
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE qs_ot,                           ONLY: qs_ot_apply_complex_frechet_dbcsr,&
                                              qs_ot_complex_exp_frechet_kernel,&
                                              qs_ot_density_tangent,&
                                              qs_ot_get_derivative_complex,&
                                              qs_ot_get_derivative_ref_complex,&
                                              qs_ot_get_orbitals_ref_complex,&
                                              qs_ot_get_p_complex
   USE qs_ot_types,                     ONLY: qs_ot_type
   USE reference_manager,               ONLY: remove_all_references

   IMPLICIT NONE

   INTEGER, PARAMETER                                 :: k = 2, n = 3

   INTEGER                                            :: io_unit, mynode, nfail
   INTEGER, DIMENSION(:), POINTER                     :: col_dist, col_k, row_dist, row_k, row_n
   INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
   TYPE(cp_logger_type), POINTER                      :: logger
   TYPE(dbcsr_distribution_type)                     :: dist
   TYPE(mp_comm_type)                                 :: mp_comm
   TYPE(mp_para_env_type), POINTER                    :: para_env

   nfail = 0
   NULLIFY (logger, para_env)
   CALL mp_world_init(mp_comm)
   mynode = mp_comm%mepos
   io_unit = -1
   IF (mynode == 0) io_unit = default_output_unit
   CALL test_density_tangent(mynode, nfail)
   ALLOCATE (para_env)
   CALL para_env%from_dup(mp_comm)
   CALL cp_logger_create(logger, para_env=para_env, default_global_unit_nr=io_unit, &
                         close_global_unit_on_dealloc=.FALSE.)
   CALL cp_add_default_logger(logger)
   CALL add_all_references()
   CALL dbcsr_init_lib(mp_comm%get_handle(), io_unit)
   ! Sparse matrix rotations and derivative pullbacks.
   CALL test_sparse_frechet_patterns(para_env, nfail)

   ALLOCATE (pgrid(0:0, 0:0), row_dist(1), col_dist(1), row_n(1), row_k(1), col_k(1))
   pgrid(:, :) = 0
   row_dist(:) = 0
   col_dist(:) = 0
   row_n(:) = n
   row_k(:) = k
   col_k(:) = k
   CALL dbcsr_distribution_new(dist, group=mp_comm_self%get_handle(), pgrid=pgrid, &
                               row_dist=row_dist, col_dist=col_dist)
   ! Both derivative fixtures use one block on MPI_COMM_SELF.
   CALL test_complex_ref_rotation(dist, row_n, row_k, col_k, para_env, nfail)
   CALL test_complex_strict_derivative(dist, row_n, row_k, col_k, para_env, nfail)
   CALL dbcsr_distribution_release(dist)
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
!> \brief Check the coupled finite-chart density tangent against a central finite difference.
!> \param mynode MPI rank
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_density_tangent(mynode, nfail)
      INTEGER, INTENT(IN)                                :: mynode
      INTEGER, INTENT(INOUT)                             :: nfail

      INTEGER, PARAMETER                                 :: nbands = 3, nrotation = 6
      REAL(KIND=dp), PARAMETER                           :: fd_step = 2.0E-5_dp, &
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
!> \brief Check the combined finite complex REF and rotation derivatives.
!> \param dist single-block matrix distribution
!> \param row_n AO block sizes
!> \param row_k occupied row block sizes
!> \param col_k occupied column block sizes
!> \param para_env parallel environment
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_complex_ref_rotation(dist, row_n, row_k, col_k, para_env, nfail)
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: dist
      INTEGER, DIMENSION(:), POINTER                     :: row_n, row_k, col_k
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(INOUT)                             :: nfail

      COMPLEX(KIND=dp), DIMENSION(k, k)                  :: direction_a, generator, gradient_a, &
                                                            inverse_sqrt, rotation
      COMPLEX(KIND=dp), DIMENSION(n, k)                  :: c_dense, c_reference, direction_x, &
                                                            gradient_x, hc_dense, q_dense, x_dense
      COMPLEX(KIND=dp), DIMENSION(n, n)                  :: hamiltonian, overlap
      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: error_c, error_ref, error_rot, fd_ref, &
                                                            fd_rot, hstep, predicted_ref, &
                                                            predicted_rot
      REAL(KIND=dp), DIMENSION(k)                        :: occupation
      REAL(KIND=dp), DIMENSION(k), TARGET                :: rotation_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
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

      CALL create_one_block(x_re, 'rotation_ref_x_re', dist, row_n, col_k)
      CALL create_one_block(x_im, 'rotation_ref_x_im', dist, row_n, col_k)
      CALL create_one_block(sx_re, 'rotation_ref_sx_re', dist, row_n, col_k)
      CALL create_one_block(sx_im, 'rotation_ref_sx_im', dist, row_n, col_k)
      CALL create_one_block(c_re, 'rotation_ref_c_re', dist, row_n, col_k)
      CALL create_one_block(c_im, 'rotation_ref_c_im', dist, row_n, col_k)
      CALL create_one_block(hc_re, 'rotation_ref_hc_re', dist, row_n, col_k)
      CALL create_one_block(hc_im, 'rotation_ref_hc_im', dist, row_n, col_k)
      CALL create_one_block(g_re, 'rotation_ref_g_re', dist, row_n, col_k)
      CALL create_one_block(g_im, 'rotation_ref_g_im', dist, row_n, col_k)
      CALL create_one_block(s_re, 'rotation_ref_s_re', dist, row_n, row_n)
      CALL create_one_block(s_im, 'rotation_ref_s_im', dist, row_n, row_n)
      CALL create_one_block(f_re, 'rotation_ref_f_re', dist, row_k, col_k)
      CALL create_one_block(f_im, 'rotation_ref_f_im', dist, row_k, col_k)
      DO i = 1, SIZE(sym_buffer)
         CALL create_one_block(sym_buffer(i), 'rotation_ref_sym', dist, row_k, col_k)
      END DO
      DO i = 1, SIZE(nonsym_buffer)
         CALL create_one_block(nonsym_buffer(i), 'rotation_ref_nonsym', dist, row_k, col_k)
      END DO
      DO i = 1, SIZE(nk_buffer)
         CALL create_one_block(nk_buffer(i), 'rotation_ref_nk', dist, row_n, col_k)
      END DO
      DO i = 1, SIZE(rot_matrix)
         CALL create_one_block(rot_matrix(i), 'rotation_ref_rot', dist, row_k, col_k)
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
      CALL cp_blacs_env_release(blacs_env)
      CALL mp_para_env_release(local_para_env)

   END SUBROUTINE test_complex_ref_rotation

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
      CALL check_error("sparse_frechet_patterns", error, 5.0E-14_dp, nfail)
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
!> \param dist single-block matrix distribution
!> \param row_n AO block sizes
!> \param row_k occupied-space row block sizes
!> \param col_k occupied-space column block sizes
!> \param para_env parallel environment
!> \param nfail accumulated number of failures
! **************************************************************************************************
   SUBROUTINE test_complex_strict_derivative(dist, row_n, row_k, col_k, para_env, nfail)
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: dist
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
      REAL(KIND=dp)                                      :: error, error_rot, fd_slope, &
                                                            fd_slope_rot, hstep, slope, slope_rot
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
         CALL create_one_block(nk(i), TRIM(name), dist, row_n, col_k)
      END DO
      DO i = 1, SIZE(kk)
         WRITE (name, '(A,I0)') 'complex_strict_kk_', i
         CALL create_one_block(kk(i), TRIM(name), dist, row_k, col_k)
      END DO
      DO i = 1, SIZE(rot_matrix)
         WRITE (name, '(A,I0)') 'complex_strict_rot_', i
         CALL create_one_block(rot_matrix(i), TRIM(name), dist, row_k, col_k)
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

      CALL get_complex_pair(nk(7), nk(8), gradient)
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
      IF (error > 5.0E-8_dp .OR. error_rot > 5.0E-8_dp) nfail = nfail + 1
      IF (para_env%is_source()) THEN
         WRITE (io_unit, '(A,3(1X,ES13.6))') 'complex finite-STRICT derivative: fd/predicted/error', &
            fd_slope, slope, error
         WRITE (io_unit, '(A,3(1X,ES13.6))') 'complex finite-STRICT rotation: fd/predicted/error', &
            fd_slope_rot, slope_rot, error_rot
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

! **************************************************************************************************
!> \brief Report a failed numerical check, including non-finite results.
!> \param label check name
!> \param error measured error
!> \param tolerance maximum error
!> \param nfail accumulated failures
! **************************************************************************************************
   SUBROUTINE check_error(label, error, tolerance, nfail)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: error, tolerance
      INTEGER, INTENT(INOUT)                             :: nfail

      IF (ieee_is_finite(error)) THEN
         IF (error <= tolerance) RETURN
      END IF
      nfail = nfail + 1
      IF (io_unit >= 0) WRITE (io_unit, '(A,2(1X,ES13.6))') TRIM(label)//': error/tolerance', error, tolerance
   END SUBROUTINE check_error

END PROGRAM qs_ot_complex_ref_unittest
