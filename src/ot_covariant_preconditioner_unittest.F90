!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM ot_covariant_preconditioner_unittest
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: ot_low_rank_base_overlap,&
                                              ot_precond_fermi_low_rank,&
                                              ot_precond_full_all_covariant
   USE iso_fortran_env,                 ONLY: error_unit
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE preconditioner_apply,            ONLY: apply_preconditioner_fm
   USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                              init_preconditioner,&
                                              preconditioner_type

   IMPLICIT NONE

   INTEGER, PARAMETER                                  :: k = 3, n = 5
   REAL(KIND=dp), PARAMETER                            :: eps = 2.0E-11_dp

   INTEGER                                             :: i, j, rank_used
   REAL(KIND=dp), DIMENSION(n, k)                      :: gradient, gradient_rotated, hessian_x, &
                                                          output, output_reference, &
                                                          output_rotated, x_reference
   REAL(KIND=dp), DIMENSION(n)                         :: full_evals
   REAL(KIND=dp), DIMENSION(k)                         :: occ_evals
   REAL(KIND=dp), DIMENSION(k, k)                      :: occupied_h, q, rotation, rotation_2, &
                                                          rotation_gauge
   REAL(KIND=dp), DIMENSION(n, n)                      :: inverse, overlap_inverse, spectral_vectors
   TYPE(cp_blacs_env_type), POINTER                    :: blacs_env
   TYPE(cp_fm_struct_type), POINTER                    :: fm_struct
   TYPE(cp_fm_type)                                    :: matrix_in, matrix_out
   TYPE(mp_para_env_type), POINTER                     :: para_env
   TYPE(preconditioner_type)                           :: preconditioner_env

   NULLIFY (blacs_env, fm_struct, para_env)
   ALLOCATE (para_env)
   CALL mp_world_init(para_env)
   CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)
   CALL init_preconditioner(preconditioner_env, para_env, blacs_env)

   CALL cp_fm_struct_create(fm_struct, nrow_global=n, ncol_global=k, &
                            context=blacs_env, para_env=para_env)
   CALL cp_fm_create(matrix_in, fm_struct, name="covariant test input")
   CALL cp_fm_create(matrix_out, fm_struct, name="covariant test output")
   CALL cp_fm_struct_release(fm_struct)

   CALL cp_fm_struct_create(fm_struct, nrow_global=n, ncol_global=n, &
                            context=blacs_env, para_env=para_env)
   ALLOCATE (preconditioner_env%fm)
   CALL cp_fm_create(preconditioner_env%fm, fm_struct, name="spectral eigenvectors")
   CALL cp_fm_struct_release(fm_struct)

   CALL cp_fm_struct_create(fm_struct, nrow_global=k, ncol_global=k, &
                            context=blacs_env, para_env=para_env)
   ALLOCATE (preconditioner_env%occ_rotation)
   CALL cp_fm_create(preconditioner_env%occ_rotation, fm_struct, name="occupied rotation")
   CALL cp_fm_struct_release(fm_struct)

   preconditioner_env%in_use = ot_precond_full_all_covariant
   preconditioner_env%energy_gap = 0.10_dp
   ALLOCATE (preconditioner_env%full_evals(n), preconditioner_env%occ_evals(k))
   full_evals = [0.8_dp, 1.3_dp, 2.1_dp, 3.0_dp, 4.4_dp]
   occ_evals = [-0.9_dp, -0.2_dp, 0.35_dp]
   preconditioner_env%full_evals = full_evals
   preconditioner_env%occ_evals = occ_evals

   spectral_vectors = 0.0_dp
   DO i = 1, n
      spectral_vectors(i, i) = 1.0_dp
   END DO
   CALL cp_fm_set_submatrix(preconditioner_env%fm, spectral_vectors)
   CALL plane_rotation(0.43_dp, 1, 2, rotation)
   CALL plane_rotation(-0.31_dp, 2, 3, q)
   rotation = MATMUL(rotation, q)
   CALL cp_fm_set_submatrix(preconditioner_env%occ_rotation, rotation)

   DO j = 1, k
      DO i = 1, n
         gradient(i, j) = SIN(0.37_dp*REAL(2*i + j, dp)) + &
                          0.2_dp*COS(0.19_dp*REAL(i - 3*j, dp))
         x_reference(i, j) = COS(0.23_dp*REAL(i + 2*j, dp)) - &
                             0.1_dp*SIN(0.41_dp*REAL(3*i - j, dp))
      END DO
   END DO

   CALL apply_fm(preconditioner_env, matrix_in, matrix_out, gradient, output_reference)

   ! Invert the frozen-H Sylvester operator A X - X B exactly.
   occupied_h = 0.0_dp
   DO j = 1, k
      occupied_h = occupied_h + occ_evals(j)*outer_product(rotation(:, j), rotation(:, j))
   END DO
   hessian_x = -MATMUL(x_reference, occupied_h)
   DO i = 1, n
      hessian_x(i, :) = hessian_x(i, :) + full_evals(i)*x_reference(i, :)
   END DO
   CALL apply_fm(preconditioner_env, matrix_in, matrix_out, hessian_x, output)
   CALL assert_close(output, x_reference, eps, "Frozen-H Sylvester inverse")

   ! Rotate the occupied gauge and verify P_(C R)(G R) = P_C(G) R.
   CALL plane_rotation(-0.37_dp, 1, 2, rotation_gauge)
   CALL plane_rotation(0.28_dp, 2, 3, q)
   rotation_gauge = MATMUL(rotation_gauge, q)
   rotation_2 = MATMUL(TRANSPOSE(rotation_gauge), rotation)
   gradient_rotated = MATMUL(gradient, rotation_gauge)
   CALL cp_fm_set_submatrix(preconditioner_env%occ_rotation, rotation_2)
   CALL apply_fm(preconditioner_env, matrix_in, matrix_out, gradient_rotated, output_rotated)
   CALL assert_close(output_rotated, MATMUL(output_reference, rotation_gauge), eps, &
                     "Occupied-gauge covariance")

   ! The gap floor must leave the inverse positive definite.
   preconditioner_env%energy_gap = 0.75_dp
   preconditioner_env%occ_evals = [0.7_dp, 0.9_dp, 1.1_dp]
   CALL cp_fm_set_submatrix(preconditioner_env%occ_rotation, rotation)
   CALL apply_fm(preconditioner_env, matrix_in, matrix_out, gradient, output)
   IF (SUM(gradient*output) <= 0.0_dp) THEN
      ERROR STOP "Gap-floor inverse is not positive definite"
   END IF

   ! Exercise the production low-rank FM path with a nontrivial overlap base,
   ! including zero rank, truncation, and a gap-clipped complementary eigenvalue.
   preconditioner_env%in_use = ot_precond_fermi_low_rank
   preconditioner_env%low_rank_base = ot_low_rank_base_overlap
   preconditioner_env%energy_gap = 0.1_dp
   preconditioner_env%spectral_reference = full_evals(k + 1) - 0.04_dp
   preconditioner_env%spectral_base_scale = 1.25_dp
   CALL plane_rotation(0.29_dp, 2, 4, spectral_vectors)
   spectral_vectors(:, 2) = 1.2_dp*spectral_vectors(:, 2)
   spectral_vectors(:, 4) = 0.8_dp*spectral_vectors(:, 4)
   overlap_inverse = MATMUL(spectral_vectors, TRANSPOSE(spectral_vectors))
   CALL cp_fm_set_submatrix(preconditioner_env%fm, spectral_vectors)
   ALLOCATE (preconditioner_env%base_fm)
   CALL cp_fm_create(preconditioner_env%base_fm, preconditioner_env%fm%matrix_struct)
   CALL cp_fm_set_submatrix(preconditioner_env%base_fm, overlap_inverse)
   DO rank_used = 0, n - k
      preconditioner_env%spectral_rank = rank_used
      inverse = preconditioner_env%spectral_base_scale*overlap_inverse
      DO i = k + 1, k + rank_used
         inverse = inverse + (1.0_dp/MAX(preconditioner_env%energy_gap, &
                                         full_evals(i) - preconditioner_env%spectral_reference) - &
                              preconditioner_env%spectral_base_scale)* &
                   outer_product(spectral_vectors(:, i), spectral_vectors(:, i))
      END DO
      output_reference = MATMUL(inverse, gradient)
      CALL apply_fm(preconditioner_env, matrix_in, matrix_out, gradient, output)
      CALL assert_close(output, output_reference, eps, "Low-rank spectral inverse")
   END DO

   CALL cp_fm_release(matrix_in)
   CALL cp_fm_release(matrix_out)
   CALL destroy_preconditioner(preconditioner_env)
   CALL cp_blacs_env_release(blacs_env)
   CALL mp_world_finalize()
   DEALLOCATE (para_env)

CONTAINS

! **************************************************************************************************
!> \brief Apply the production FM path to a replicated test matrix.
!> \param preconditioner_env ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param input ...
!> \param RESULT ...
! **************************************************************************************************
   SUBROUTINE apply_fm(preconditioner_env, matrix_in, matrix_out, input, RESULT)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_in, matrix_out
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: input
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: result

      CALL cp_fm_set_submatrix(matrix_in, input)
      CALL apply_preconditioner_fm(preconditioner_env, matrix_in, matrix_out)
      CALL cp_fm_get_submatrix(matrix_out, RESULT)

   END SUBROUTINE apply_fm

! **************************************************************************************************
!> \brief Create an orthogonal plane rotation.
!> \param angle ...
!> \param axis_1 ...
!> \param axis_2 ...
!> \param rotation ...
! **************************************************************************************************
   SUBROUTINE plane_rotation(angle, axis_1, axis_2, rotation)

      REAL(KIND=dp), INTENT(IN)                          :: angle
      INTEGER, INTENT(IN)                                :: axis_1, axis_2
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: rotation

      INTEGER                                            :: i

      rotation = 0.0_dp
      DO i = 1, SIZE(rotation, 1)
         rotation(i, i) = 1.0_dp
      END DO
      rotation(axis_1, axis_1) = COS(angle)
      rotation(axis_2, axis_2) = COS(angle)
      rotation(axis_1, axis_2) = -SIN(angle)
      rotation(axis_2, axis_1) = SIN(angle)

   END SUBROUTINE plane_rotation

! **************************************************************************************************
!> \brief Form a real outer product.
!> \param left ...
!> \param right ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION outer_product(left, right) RESULT(product)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: left, right
      REAL(KIND=dp), DIMENSION(SIZE(left), SIZE(right))  :: product

      INTEGER                                            :: i, j

      DO j = 1, SIZE(right)
         DO i = 1, SIZE(left)
            PRODUCT(i, j) = left(i)*right(j)
         END DO
      END DO

   END FUNCTION outer_product

! **************************************************************************************************
!> \brief Abort when two matrices differ beyond the requested tolerance.
!> \param actual ...
!> \param reference ...
!> \param tolerance ...
!> \param label ...
! **************************************************************************************************
   SUBROUTINE assert_close(actual, reference, tolerance, label)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: actual, reference
      REAL(KIND=dp), INTENT(IN)                          :: tolerance
      CHARACTER(LEN=*), INTENT(IN)                       :: label

      INTEGER                                            :: io_unit

      IF (.NOT. ALL(ieee_is_finite(actual)) .OR. &
          .NOT. ALL(ieee_is_finite(reference)) .OR. MAXVAL(ABS(actual - reference)) > tolerance) THEN
         io_unit = error_unit
         WRITE (io_unit, '(A,2(1X,ES13.6))') TRIM(label)//": error/tolerance", &
            MAXVAL(ABS(actual - reference)), tolerance
         ERROR STOP "Matrix comparison failed"
      END IF

   END SUBROUTINE assert_close

END PROGRAM ot_covariant_preconditioner_unittest
