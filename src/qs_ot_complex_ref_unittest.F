!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Finite-difference regression for the complex finite-REF pullback.
! **************************************************************************************************
PROGRAM qs_ot_complex_ref_unittest
   USE bibliography,                    ONLY: add_all_references
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_finalize_lib, dbcsr_get_block_p, dbcsr_init_lib, dbcsr_put_block, &
        dbcsr_release, dbcsr_reserve_blocks, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_logger_create,&
                                              cp_logger_release,&
                                              cp_logger_type,&
                                              cp_rm_default_logger
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: default_output_unit
   USE mathlib,                         ONLY: diag_complex
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE qs_ot,                           ONLY: qs_ot_get_derivative_ref_complex
   USE qs_ot_types,                     ONLY: qs_ot_type
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
      WRITE (*, '(A,3(1X,ES13.6))') "complex finite-REF pullback: fd/predicted/error", &
         fd_slope, slope, error
      WRITE (*, '(A,1X,ES13.6)') "discarded projected-gradient error", old_error
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

END PROGRAM qs_ot_complex_ref_unittest
