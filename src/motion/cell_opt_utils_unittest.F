!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM cell_opt_utils_TEST
   USE cell_methods,                    ONLY: cell_create
   USE cell_opt_utils,                  ONLY: &
        apply_cell_constraints, cell_opt_constraint_label, comp_ax, comp_bx, comp_by, comp_cz, &
        keep_angles_conflict, keep_symmetry_conflict, mask_all, mask_b, mask_c, &
        parse_cell_constraint, parse_token, pres_crit_shear_pinned
   USE cell_types,                      ONLY: &
        cell_release, cell_sym_cubic, cell_sym_hexagonal_gamma_120, cell_sym_hexagonal_gamma_60, &
        cell_sym_monoclinic, cell_sym_monoclinic_gamma_ab, cell_sym_none, cell_sym_orthorhombic, &
        cell_sym_rhombohedral, cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, &
        cell_sym_tetragonal_bc, cell_sym_triclinic, cell_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE string_utilities,                ONLY: remove_word

   IMPLICIT NONE

   REAL(KIND=dp), PARAMETER :: initial(6) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
   REAL(KIND=dp), DIMENSION(:), POINTER :: gradient => NULL()
   TYPE(cell_type), POINTER :: cell => NULL()

   CALL cell_create(cell)
   ALLOCATE (gradient(6))
   CALL check_parser()
   CALL check_projections()
   CALL check_guards()
   DEALLOCATE (gradient)
   CALL cell_release(cell)

CONTAINS

! **************************************************************************************************
!> \brief Check token masks, label round trips, rejections and 2D defaults.
! **************************************************************************************************
   SUBROUTINE check_parser()
      CHARACTER(LEN=4), PARAMETER :: &
         legacy(7) = ["NONE", "   X", "   Y", "   Z", "  XY", "  XZ", "  YZ"], names(18) = [ &
         CHARACTER(LEN=4) :: "NONE", "X", "Y", "Z", "XY", "XZ", "YZ", "a", "AX", "B", "BX", "BY", &
         "C", "CX", "CY", "CZ", "none", " z "]
      CHARACTER(LEN=5), PARAMETER :: &
         invalid(5) = [CHARACTER(LEN=5) :: "AY", "AZ", "BZ", "ALPHA", "YX"]
      INTEGER, PARAMETER :: &
         masks(18) = [0, 11, 22, 56, 31, 59, 62, 1, 1, 6, 2, 4, 56, 8, 16, 32, 0, 56]

      CHARACTER(LEN=default_string_length)               :: label, reason, tokens(6)
      INTEGER                                            :: i, mask, n, wanted
      LOGICAL                                            :: accepted

! Bits 0..5 represent AX, BX, BY, CX, CY, CZ, respectively.

      cell%perd = 1
      DO i = 1, SIZE(names)
         CALL parse_cell_constraint([names(i)], cell, .TRUE., mask)
         IF (mask /= masks(i)) ERROR STOP "Incorrect constraint token mask"
         IF (i <= SIZE(legacy)) THEN
            IF (cell_opt_constraint_label(mask) /= legacy(i)) ERROR STOP "Legacy constraint label changed"
         END IF
      END DO
      CALL parse_cell_constraint(["AX", "BY", "CZ"], cell, .TRUE., mask)
      IF (mask /= 37) ERROR STOP "Constraint components were not combined"
      CALL parse_cell_constraint(["a ", "b ", "ax"], cell, .TRUE., mask)
      IF (mask /= 7) ERROR STOP "Redundant constraint changed the mask"
      CALL parse_cell_constraint(["X", "B"], cell, .TRUE., mask)
      IF (mask /= 15) ERROR STOP "Partially overlapping constraint lost components"

      DO wanted = 0, 63
         label = ADJUSTL(cell_opt_constraint_label(wanted))
         n = 0
         DO WHILE (LEN_TRIM(label) > 0)
            n = n + 1
            READ (label, *) tokens(n)
            CALL remove_word(label)
            label = ADJUSTL(label)
         END DO
         CALL parse_cell_constraint(tokens(:n), cell, .TRUE., mask)
         IF (mask /= wanted) ERROR STOP "Constraint label round trip failed"
      END DO
      DO i = 1, SIZE(invalid)
         CALL parse_token(invalid(i), accepted, mask, reason)
         IF (accepted .OR. LEN_TRIM(reason) == 0) ERROR STOP "Invalid constraint token accepted"
      END DO
      DO i = 1, 3
         cell%perd = 1
         cell%perd(i) = 0
         CALL parse_cell_constraint(["NONE"], cell, .FALSE., mask)
         IF (mask /= masks(i + 1)) ERROR STOP "Incorrect 2D default constraint"
         CALL parse_cell_constraint([names(i + 1)], cell, .TRUE., mask)
         IF (mask /= masks(i + 1)) ERROR STOP "Explicit 2D constraint changed"
         CALL parse_cell_constraint(["A", "B", "C"], cell, .TRUE., mask)
         IF (mask /= mask_all) ERROR STOP "Explicit 2D constraint superset changed"
      END DO
   END SUBROUTINE check_parser

! **************************************************************************************************
!> \brief Check each projection against a known gradient and apply it twice.
! **************************************************************************************************
   SUBROUTINE check_projections()
      REAL(KIND=dp)                                      :: a, normal(3), reference(6), root3

      root3 = SQRT(3.0_dp)
      cell%perd = 1
      cell%orthorhombic = .FALSE.
      cell%hmat = RESHAPE([2.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], [3, 3])
      CALL check_projection(cell_sym_none, initial, .FALSE., 0)
      CALL check_projection(cell_sym_triclinic, initial, .FALSE., 0)
      CALL check_projection(cell_sym_triclinic, &
                            [1.0_dp, 8.0_dp/5, 16.0_dp/5, 26.0_dp/9, 52.0_dp/9, 52.0_dp/9], .TRUE., 0)
      CALL check_projection(cell_sym_triclinic, [0.0_dp, 2.0_dp, 0.0_dp, 4.0_dp, 5.0_dp, 0.0_dp], &
                            .FALSE., comp_ax + comp_by + comp_cz)
      CALL check_projection(cell_sym_triclinic, [0.0_dp, 8.0_dp/5, 16.0_dp/5, 0.0_dp, 0.0_dp, 0.0_dp], &
                            .TRUE., comp_ax + mask_c)
      cell%hmat = RESHAPE([2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], [3, 3])
      cell%orthorhombic = .TRUE.
      CALL check_projection(cell_sym_cubic, [10.0_dp/3, 0.0_dp, 10.0_dp/3, 0.0_dp, 0.0_dp, 10.0_dp/3], .FALSE., 0)
      CALL check_projection(cell_sym_tetragonal_ab, [2.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 6.0_dp], .FALSE., 0)
      CALL check_projection(cell_sym_tetragonal_ac, [3.5_dp, 0.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 3.5_dp], .FALSE., 0)
      CALL check_projection(cell_sym_tetragonal_bc, [1.0_dp, 0.0_dp, 4.5_dp, 0.0_dp, 0.0_dp, 4.5_dp], .FALSE., 0)
      CALL check_projection(cell_sym_monoclinic, [1.0_dp, 0.0_dp, 3.0_dp, 4.0_dp, 0.0_dp, 6.0_dp], .FALSE., 0)
      cell%hmat = RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, root3/2, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], [3, 3])
      cell%orthorhombic = .FALSE.
      a = 1.0_dp + 0.75_dp*root3
      CALL check_projection(cell_sym_hexagonal_gamma_60, [a, a/2, root3*a/2, 0.0_dp, 0.0_dp, 6.0_dp], .FALSE., 0)
      a = 0.75_dp*root3
      cell%hmat(1, 2) = -0.5_dp
      CALL check_projection(cell_sym_hexagonal_gamma_120, [a, -a/2, root3*a/2, 0.0_dp, 0.0_dp, 6.0_dp], .FALSE., 0)

      cell%hmat = RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, root3/2, 0.0_dp, &
                           0.5_dp, root3/6, SQRT(2.0_dp/3)], [3, 3])
      ! Equal lengths and 60-degree angles leave one radial direction.
      reference = [1.0_dp, 0.5_dp, root3/2, 0.5_dp, root3/6, SQRT(2.0_dp/3)]
      reference = reference*DOT_PRODUCT(reference, initial)/3.0_dp
      CALL check_projection(cell_sym_rhombohedral, reference, .FALSE., 0)
      cell%hmat(:, 3) = [0.0_dp, 0.0_dp, 2.0_dp]
      ! Equal a/b lengths require -dAX + cos(gamma)*dBX + sin(gamma)*dBY = 0.
      normal = [-1.0_dp, 0.5_dp, root3/2]
      reference = initial
      reference(1:3) = initial(1:3) - normal*DOT_PRODUCT(normal, initial(1:3))/2.0_dp
      reference(4:5) = 0.0_dp
      CALL check_projection(cell_sym_monoclinic_gamma_ab, reference, .FALSE., 0)

      cell%hmat = 0.0_dp
      cell%hmat(1, 1) = 2.0_dp
      cell%hmat(2, 2) = 3.0_dp
      cell%hmat(3, 3) = 4.0_dp
      cell%orthorhombic = .TRUE.
      CALL check_projection(cell_sym_orthorhombic, [1.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 6.0_dp], .FALSE., 0)
      CALL check_projection(cell_sym_orthorhombic, [1.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 6.0_dp], .TRUE., 0)
      CALL check_projection(cell_sym_cubic, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], .TRUE., mask_all)
   END SUBROUTINE check_projections

! **************************************************************************************************
!> \brief Compare projected and expected gradients, including idempotence.
!> \param symmetry cell symmetry
!> \param expected projected gradient
!> \param keep_angles apply the angle projection instead of the symmetry projection
!> \param mask fixed components
! **************************************************************************************************
   SUBROUTINE check_projection(symmetry, expected, keep_angles, mask)
      INTEGER, INTENT(IN)                                :: symmetry
      REAL(KIND=dp), INTENT(IN)                          :: expected(6)
      LOGICAL, INTENT(IN)                                :: keep_angles
      INTEGER, INTENT(IN)                                :: mask

      INTEGER                                            :: repeat

      cell%symmetry_id = symmetry
      gradient = initial
      DO repeat = 1, 2
         CALL apply_cell_constraints(gradient, cell, keep_angles,.NOT. keep_angles, mask)
         IF (.NOT. ALL(ABS(gradient - expected) <= 1.0E-12_dp)) THEN
            PRINT *, "Incorrect cell projection: symmetry, mask, repetition:", symmetry, mask, repeat
            ERROR STOP 1
         END IF
      END DO
   END SUBROUTINE check_projection

! **************************************************************************************************
!> \brief Check representative legal/conflicting constraints and pressure warnings.
! **************************************************************************************************
   SUBROUTINE check_guards()
      INTEGER, PARAMETER :: &
         masks(8) = [0, comp_ax, mask_b, mask_c, comp_bx, comp_by, comp_cz, mask_all], &
         partial(12) = [1, 1, 1, 4, 2, 2, 2, 32, 1, 1, 1, 1], systems(12) = [cell_sym_cubic, &
         cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, cell_sym_tetragonal_bc, &
         cell_sym_hexagonal_gamma_60, cell_sym_hexagonal_gamma_120, cell_sym_monoclinic_gamma_ab, &
         cell_sym_rhombohedral, cell_sym_orthorhombic, cell_sym_monoclinic, cell_sym_triclinic, &
         cell_sym_none]
      LOGICAL, PARAMETER :: angle_conflict(8) = [.FALSE., .FALSE., .FALSE., .FALSE., .TRUE., .TRUE.&
         , .TRUE., .FALSE.], &
         pressure(8) = [.FALSE., .FALSE., .TRUE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE.]

      INTEGER                                            :: i

      DO i = 1, SIZE(systems)
         IF (keep_symmetry_conflict(0, systems(i))) ERROR STOP "Unconstrained symmetry rejected"
         IF (keep_symmetry_conflict(mask_all, systems(i))) ERROR STOP "Fully fixed symmetry rejected"
         IF (keep_symmetry_conflict(partial(i), systems(i)) .NEQV. (i <= 8)) ERROR STOP "Incorrect symmetry conflict"
      END DO
      IF (keep_symmetry_conflict(comp_ax + comp_by, cell_sym_tetragonal_ab)) ERROR STOP "Complete ab constraint rejected"
      DO i = 1, SIZE(masks)
         IF (keep_angles_conflict(masks(i), .FALSE.) .NEQV. angle_conflict(i)) ERROR STOP "Incorrect angle conflict"
         IF (keep_angles_conflict(masks(i), .TRUE.)) ERROR STOP "Orthorhombic angle constraint rejected"
         IF (pres_crit_shear_pinned(masks(i), .FALSE., .FALSE.) .NEQV. pressure(i)) ERROR STOP "Incorrect pressure warning"
         IF (pres_crit_shear_pinned(masks(i), .TRUE., .TRUE.)) ERROR STOP "Orthorhombic pressure warning"
      END DO
      IF (.NOT. pres_crit_shear_pinned(comp_ax, .FALSE., .TRUE.)) ERROR STOP "Missing pinned-angle pressure warning"
      IF (pres_crit_shear_pinned(37 + comp_bx, .FALSE., .TRUE.)) ERROR STOP "Pressure warning with all diagonals fixed"
   END SUBROUTINE check_guards

END PROGRAM cell_opt_utils_TEST
