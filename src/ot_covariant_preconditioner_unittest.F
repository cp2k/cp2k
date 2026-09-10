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
   USE input_constants,                 ONLY: ot_precond_full_all_covariant
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

   INTEGER                                             :: i, j
   REAL(KIND=dp), DIMENSION(n, k)                      :: gradient, gradient_rotated, hessian_x, &
                                                          output, output_reference, &
                                                          output_rotated, x_reference
   REAL(KIND=dp), DIMENSION(n)                         :: full_evals
   REAL(KIND=dp), DIMENSION(k)                         :: occ_evals
   REAL(KIND=dp), DIMENSION(k, k)                      :: occupied_h, q, rotation, rotation_2, &
                                                          rotation_gauge
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

   CALL set_identity(preconditioner_env%fm, n)
   CALL make_rotation(0.43_dp, -0.31_dp, rotation)
   CALL cp_fm_set_submatrix(preconditioner_env%occ_rotation, rotation)

   DO j = 1, k
      DO i = 1, n
         gradient(i, j) = SIN(0.37_dp*REAL(2*i + j, dp)) + &
                          0.2_dp*COS(0.19_dp*REAL(i - 3*j, dp))
         x_reference(i, j) = COS(0.23_dp*REAL(i + 2*j, dp)) - &
                             0.1_dp*SIN(0.41_dp*REAL(3*i - j, dp))
      END DO
   END DO

   ! Check the production FM wrapper against its explicit spectral formula.
   CALL apply_model(gradient, rotation, full_evals, occ_evals, &
                    preconditioner_env%energy_gap, output_reference)
   CALL apply_fm(preconditioner_env, matrix_in, matrix_out, gradient, output)
   CALL assert_close(output, output_reference, eps, "Explicit covariant spectral inverse")

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
   CALL make_rotation(-0.37_dp, 0.28_dp, rotation_gauge)
   rotation_2 = MATMUL(TRANSPOSE(rotation_gauge), rotation)
   gradient_rotated = MATMUL(gradient, rotation_gauge)
   CALL cp_fm_set_submatrix(preconditioner_env%occ_rotation, rotation_2)
   CALL apply_fm(preconditioner_env, matrix_in, matrix_out, gradient_rotated, output_rotated)
   CALL assert_close(output_rotated, MATMUL(output_reference, rotation_gauge), eps, &
                     "Occupied-gauge covariance")

   ! Eigenvectors may rotate freely inside a degenerate occupied eigenspace.
   occ_evals = [-0.4_dp, -0.4_dp, 0.2_dp]
   preconditioner_env%occ_evals = occ_evals
   CALL cp_fm_set_submatrix(preconditioner_env%occ_rotation, rotation)
   CALL apply_fm(preconditioner_env, matrix_in, matrix_out, gradient, output)
   CALL plane_rotation(0.61_dp, 1, 2, q)
   rotation_2 = MATMUL(rotation, q)
   CALL cp_fm_set_submatrix(preconditioner_env%occ_rotation, rotation_2)
   CALL apply_fm(preconditioner_env, matrix_in, matrix_out, gradient, output_rotated)
   CALL assert_close(output_rotated, output, eps, "Degenerate occupied subspace")

   ! The gap floor must leave the inverse positive definite.
   preconditioner_env%energy_gap = 0.75_dp
   preconditioner_env%occ_evals = [0.7_dp, 0.9_dp, 1.1_dp]
   CALL cp_fm_set_submatrix(preconditioner_env%occ_rotation, rotation)
   CALL apply_fm(preconditioner_env, matrix_in, matrix_out, gradient, output)
   IF (SUM(gradient*output) <= 0.0_dp) THEN
      ERROR STOP "Gap-floor inverse is not positive definite"
   END IF

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
!> \brief Explicit reference for the covariant state-selective spectral inverse.
!> \param input ...
!> \param rotation ...
!> \param full_evals ...
!> \param occ_evals ...
!> \param gap ...
!> \param RESULT ...
! **************************************************************************************************
   SUBROUTINE apply_model(input, rotation, full_evals, occ_evals, gap, RESULT)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: input, rotation
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: full_evals, occ_evals
      REAL(KIND=dp), INTENT(IN)                          :: gap
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: result

      INTEGER                                            :: i, j
      REAL(KIND=dp), &
         DIMENSION(SIZE(input, 1), SIZE(input, 2))       :: canonical

      canonical = MATMUL(input, rotation)
      DO j = 1, SIZE(input, 2)
         DO i = 1, SIZE(input, 1)
            canonical(i, j) = canonical(i, j)/MAX(gap, full_evals(i) - occ_evals(j))
         END DO
      END DO
      RESULT = MATMUL(canonical, TRANSPOSE(rotation))

   END SUBROUTINE apply_model

! **************************************************************************************************
!> \brief Create a product of two plane rotations.
!> \param angle_12 ...
!> \param angle_23 ...
!> \param rotation ...
! **************************************************************************************************
   SUBROUTINE make_rotation(angle_12, angle_23, rotation)

      REAL(KIND=dp), INTENT(IN)                          :: angle_12, angle_23
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: rotation

      REAL(KIND=dp), &
         DIMENSION(SIZE(rotation, 1), SIZE(rotation, 2)) :: r12, r23

      CALL plane_rotation(angle_12, 1, 2, r12)
      CALL plane_rotation(angle_23, 2, 3, r23)
      rotation = MATMUL(r12, r23)

   END SUBROUTINE make_rotation

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
!> \brief Set a distributed full matrix to the identity.
!> \param matrix ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE set_identity(matrix, n)

      TYPE(cp_fm_type), INTENT(IN)                       :: matrix
      INTEGER, INTENT(IN)                                :: n

      INTEGER                                            :: i
      REAL(KIND=dp), DIMENSION(n, n)                     :: identity

      identity = 0.0_dp
      DO i = 1, n
         identity(i, i) = 1.0_dp
      END DO
      CALL cp_fm_set_submatrix(matrix, identity)

   END SUBROUTINE set_identity

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

      IF (MAXVAL(ABS(actual - reference)) > tolerance) THEN
         WRITE (*, '(A)') TRIM(label)
         ERROR STOP "Matrix comparison failed"
      END IF

   END SUBROUTINE assert_close

END PROGRAM ot_covariant_preconditioner_unittest
