!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM cell_opt_utils_TEST
   USE cell_manifold_utils,             ONLY: &
        cell_manifold_chart_eval, cell_manifold_geometry, cell_manifold_init, cell_manifold_type, &
        geom_a, geom_alpha, geom_b, geom_beta, geom_c, geom_gamma, info_ok, n_geom, n_slot, &
        slot_ax, slot_bx
   USE cell_methods,                    ONLY: cell_create
   USE cell_opt_utils,                  ONLY: cell_opt_constraint_label,&
                                              mask_all,&
                                              parse_cell_constraint,&
                                              parse_token,&
                                              spgr_point_group_symmetry_id,&
                                              symmetry_constraint_rows
   USE cell_types,                      ONLY: &
        cell_release, cell_sym_cubic, cell_sym_hexagonal_gamma_120, cell_sym_hexagonal_gamma_60, &
        cell_sym_monoclinic, cell_sym_monoclinic_gamma_ab, cell_sym_none, cell_sym_orthorhombic, &
        cell_sym_rhombohedral, cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, &
        cell_sym_tetragonal_bc, cell_sym_triclinic, cell_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi,&
                                              sqrt3
   USE space_groups_types,              ONLY: release_spgr_type,&
                                              spgr_type
   USE string_utilities,                ONLY: remove_word

   IMPLICIT NONE

   TYPE(cell_type), POINTER :: cell => NULL()
   TYPE(spgr_type), POINTER :: spgr => NULL()

   CALL cell_create(cell)
   ALLOCATE (spgr)
   CALL check_parser()
   CALL check_symmetry_manifolds()
   CALL check_spgr_point_group()
   CALL release_spgr_type(spgr)
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
!> \brief Check row counts, signs, 0-step flatness, and curved Newton steps
!>        across all crystal systems.
! **************************************************************************************************
   SUBROUTINE check_symmetry_manifolds()
      INTEGER, PARAMETER :: ndof_exp(12) = [3, 2, 2, 2, 1, 2, 2, 2, 4, 3, 6, 6], sid(12) = [ &
         cell_sym_orthorhombic, cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, &
         cell_sym_tetragonal_bc, cell_sym_cubic, cell_sym_hexagonal_gamma_60, &
         cell_sym_hexagonal_gamma_120, cell_sym_rhombohedral, cell_sym_monoclinic, &
         cell_sym_monoclinic_gamma_ab, cell_sym_triclinic, cell_sym_none]
      REAL(KIND=dp), PARAMETER :: g(n_slot) = [0.4_dp, -1.2_dp, 0.7_dp, 1.5_dp, -0.3_dp, 0.9_dp]

      INTEGER                                            :: eq_pair(2, n_geom), i, info, mask, n_eq, &
                                                            n_lin
      REAL(KIND=dp)                                      :: ct, cx, cy, cz, lin_w(n_slot, n_slot), &
                                                            q0(n_slot, 12), qout(n_slot), &
                                                            s(n_geom), st, u(n_slot), worst
      TYPE(cell_manifold_type)                           :: m

      ct = COS(70.0_dp*pi/180.0_dp)
      st = SIN(70.0_dp*pi/180.0_dp)
      cx = ct
      cy = (ct - ct*ct)/st
      cz = SQRT(MAX(1.0_dp - cx*cx - cy*cy, 0.0_dp))
      q0(:, 1) = [5.1_dp, 0.0_dp, 6.2_dp, 0.0_dp, 0.0_dp, 7.3_dp]
      q0(:, 2) = [5.1_dp, 0.0_dp, 5.1_dp, 0.0_dp, 0.0_dp, 6.3_dp]
      q0(:, 3) = [5.1_dp, 0.0_dp, 6.2_dp, 0.0_dp, 0.0_dp, 5.1_dp]
      q0(:, 4) = [6.2_dp, 0.0_dp, 5.1_dp, 0.0_dp, 0.0_dp, 5.1_dp]
      q0(:, 5) = [5.1_dp, 0.0_dp, 5.1_dp, 0.0_dp, 0.0_dp, 5.1_dp]
      q0(:, 6) = [5.1_dp, 2.55_dp, 0.5_dp*sqrt3*5.1_dp, 0.0_dp, 0.0_dp, 6.3_dp]
      q0(:, 7) = [5.1_dp, -2.55_dp, 0.5_dp*sqrt3*5.1_dp, 0.0_dp, 0.0_dp, 6.3_dp]
      q0(:, 8) = 5.1_dp*[1.0_dp, ct, st, cx, cy, cz]
      q0(:, 9) = [5.1_dp, 0.0_dp, 6.2_dp, 2.1_dp, 0.0_dp, 7.3_dp]
      q0(:, 10) = [5.1_dp, 5.1_dp*COS(100.0_dp*pi/180.0_dp), &
                   5.1_dp*SIN(100.0_dp*pi/180.0_dp), 0.0_dp, 0.0_dp, 7.3_dp]
      q0(:, 11) = [5.1_dp, 0.7_dp, 6.2_dp, 1.3_dp, 0.9_dp, 7.3_dp]
      q0(:, 12) = q0(:, 11)

      DO i = 1, SIZE(sid)
         CALL symmetry_constraint_rows(sid(i), mask, lin_w, n_lin, eq_pair, n_eq)
         IF (n_slot - (POPCNT(mask) + n_lin + n_eq) /= ndof_exp(i)) THEN
            ERROR STOP "Incorrect symmetry row count"
         END IF
         IF (sid(i) == cell_sym_hexagonal_gamma_120 .AND. (lin_w(slot_ax, 1) <= 0.0_dp .OR. &
                                                           lin_w(slot_bx, 1) /= 1.0_dp)) THEN
            ERROR STOP "Incorrect hexagonal 120 BX sign"
         END IF
         IF (sid(i) == cell_sym_hexagonal_gamma_60 .AND. lin_w(slot_ax, 1) >= 0.0_dp) THEN
            ERROR STOP "Incorrect hexagonal 60 BX sign"
         END IF

         CALL cell_manifold_init(m, q0(:, i), mask, info, lin_w=lin_w(:, 1:n_lin), &
                                 eq_pair=eq_pair(:, 1:n_eq))
         IF (info /= info_ok .OR. m%ndof /= ndof_exp(i)) THEN
            ERROR STOP "Symmetry constraint init failed"
         END IF
         u(1:m%ndof) = m%u0(1:m%ndof) + 0.05_dp*MATMUL(TRANSPOSE(m%Z(:, 1:m%ndof)), g)
         CALL cell_manifold_chart_eval(m, u(1:m%ndof), qout, info)
         IF (info /= info_ok) ERROR STOP "Symmetry chart evaluation failed"
         IF (sid(i) == cell_sym_rhombohedral) THEN
            CALL cell_manifold_geometry(qout, s, info)
            worst = MAX(ABS(s(geom_a) - s(geom_b)), ABS(s(geom_b) - s(geom_c)), &
                        ABS(s(geom_alpha) - s(geom_beta)), ABS(s(geom_beta) - s(geom_gamma)))
            IF (m%last_iter < 1 .OR. worst > 1.0E-11_dp*5.1_dp) THEN
               ERROR STOP "Rhombohedral chart evaluation failed"
            END IF
         ELSE IF (sid(i) == cell_sym_monoclinic_gamma_ab) THEN
            CALL cell_manifold_geometry(qout, s, info)
            IF (m%last_iter < 1 .OR. ABS(s(geom_a) - s(geom_b)) > 1.0E-11_dp*5.1_dp) THEN
               ERROR STOP "Monoclinic gamma_ab chart evaluation failed"
            END IF
         ELSE IF (m%last_iter /= 0) THEN
            ERROR STOP "Flat crystal system took nonzero Newton iterations"
         END IF
      END DO
   END SUBROUTINE check_symmetry_manifolds

! **************************************************************************************************
!> \brief Check space-group point-group classification across all 11 crystal systems and edge cases.
! **************************************************************************************************
   SUBROUTINE check_spgr_point_group()
      INTEGER, DIMENSION(3, 3), PARAMETER :: &
         r_c2x = RESHAPE([1, 0, 0, 0, -1, 0, 0, 0, -1], [3, 3]), &
         r_c2y = RESHAPE([-1, 0, 0, 0, 1, 0, 0, 0, -1], [3, 3]), &
         r_c2z = RESHAPE([-1, 0, 0, 0, -1, 0, 0, 0, 1], [3, 3]), &
         r_c3_111 = RESHAPE([0, 0, 1, 1, 0, 0, 0, 1, 0], [3, 3]), &
         r_c3_120 = RESHAPE([0, -1, 0, 1, -1, 0, 0, 0, 1], [3, 3]), &
         r_c4x = RESHAPE([1, 0, 0, 0, 0, -1, 0, 1, 0], [3, 3]), &
         r_c4y = RESHAPE([0, 0, 1, 0, 1, 0, -1, 0, 0], [3, 3]), &
         r_c4z = RESHAPE([0, -1, 0, 1, 0, 0, 0, 0, 1], [3, 3]), &
         r_c6_60 = RESHAPE([1, 1, 0, -1, 0, 0, 0, 0, 1], [3, 3]), &
         r_id = RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
      INTEGER, DIMENSION(3, 3), PARAMETER :: r_mz = RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, -1], [3, 3]), &
         r_swap_ab = RESHAPE([0, 1, 0, 1, 0, 0, 0, 0, 1], [3, 3])

      INTEGER                                            :: k
      INTEGER, DIMENSION(3, 3)                           :: rk
      INTEGER, DIMENSION(3, 3, 12)                       :: rots
      REAL(KIND=dp)                                      :: ct, cx, cy, cz, st

      cell%hmat = RESHAPE([5.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 5.0_dp], [3, 3])
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_triclinic) THEN
         ERROR STOP "Empty rotation set did not return cell_sym_triclinic"
      END IF
      cell%hmat = 0.0_dp
      CALL check_pg_cyclic(r_id, 1, cell_sym_triclinic)

      cell%hmat = RESHAPE([5.1_dp, 0.0_dp, 0.0_dp, 0.7_dp, 6.2_dp, 0.0_dp, &
                           1.3_dp, 0.9_dp, 7.3_dp], [3, 3])
      CALL check_pg_cyclic(r_id, 1, cell_sym_triclinic)
      CALL check_pg_cyclic(r_c2x, 2, cell_sym_triclinic)
      CALL check_pg_cyclic(r_c2y, 2, cell_sym_monoclinic)
      CALL check_pg_cyclic(r_c4z, 4, cell_sym_tetragonal_ab)
      CALL check_pg_cyclic(r_c4y, 4, cell_sym_tetragonal_ac)
      CALL check_pg_cyclic(r_c4x, 4, cell_sym_tetragonal_bc)

      rots(:, :, 1) = r_id
      rots(:, :, 2) = r_c2x
      rots(:, :, 3) = r_c2y
      rots(:, :, 4) = r_c2z
      CALL check_pg(rots(:, :, 1:4), cell_sym_orthorhombic)
      DO k = 1, 4
         rots(:, :, 4 + k) = MATMUL(rots(:, :, k), r_c3_111)
         rots(:, :, 8 + k) = MATMUL(rots(:, :, 4 + k), r_c3_111)
      END DO
      CALL check_pg(rots(:, :, 1:12), cell_sym_cubic)

      cell%hmat = RESHAPE([5.1_dp, 0.0_dp, 0.0_dp, 2.55_dp, 0.5_dp*sqrt3*5.1_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 6.3_dp], [3, 3])
      CALL check_pg_cyclic(r_c6_60, 6, cell_sym_hexagonal_gamma_60)

      cell%hmat = RESHAPE([5.1_dp, 0.0_dp, 0.0_dp, -2.55_dp, 0.5_dp*sqrt3*5.1_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 6.3_dp], [3, 3])
      rk = r_id
      DO k = 1, 3
         rots(:, :, k) = rk
         rots(:, :, 3 + k) = MATMUL(rk, r_mz)
         rk = MATMUL(rk, r_c3_120)
      END DO
      CALL check_pg(rots(:, :, 1:6), cell_sym_hexagonal_gamma_120)

      cell%hmat = RESHAPE([5.1_dp, 0.0_dp, 0.0_dp, &
                           5.1_dp*COS(100.0_dp*pi/180.0_dp), 5.1_dp*SIN(100.0_dp*pi/180.0_dp), &
                           0.0_dp, 0.0_dp, 0.0_dp, 7.3_dp], [3, 3])
      rots(:, :, 1) = r_id
      rots(:, :, 2) = r_swap_ab
      rots(:, :, 3) = r_mz
      rots(:, :, 4) = MATMUL(r_swap_ab, r_mz)
      CALL check_pg(rots(:, :, 1:4), cell_sym_monoclinic_gamma_ab)

      ct = COS(70.0_dp*pi/180.0_dp)
      st = SIN(70.0_dp*pi/180.0_dp)
      cx = ct
      cy = (ct - ct*ct)/st
      cz = SQRT(MAX(1.0_dp - cx*cx - cy*cy, 0.0_dp))
      cell%hmat = 5.1_dp*RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, ct, st, 0.0_dp, cx, cy, cz], [3, 3])
      CALL check_pg_cyclic(r_c3_111, 3, cell_sym_rhombohedral)
   END SUBROUTINE check_spgr_point_group

! **************************************************************************************************
!> \brief Check point-group classification for the cyclic group from gen.
!> \param gen integer rotation generator in lattice coordinates
!> \param order_n group order
!> \param expected expected cell_sym_* identifier
! **************************************************************************************************
   SUBROUTINE check_pg_cyclic(gen, order_n, expected)
      INTEGER, DIMENSION(3, 3), INTENT(IN)               :: gen
      INTEGER, INTENT(IN)                                :: order_n, expected

      INTEGER                                            :: k
      INTEGER, DIMENSION(3, 3)                           :: rk
      INTEGER, DIMENSION(3, 3, order_n)                  :: rots

      rk = RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
      DO k = 1, order_n
         rots(:, :, k) = rk
         rk = MATMUL(rk, gen)
      END DO
      CALL check_pg(rots, expected)
   END SUBROUTINE check_pg_cyclic

! **************************************************************************************************
!> \brief Populate spgr%rotations_subset and assert spgr_point_group_symmetry_id.
!> \param rots integer rotation matrices in lattice coordinates
!> \param expected expected cell_sym_* identifier
! **************************************************************************************************
   SUBROUTINE check_pg(rots, expected)
      INTEGER, DIMENSION(:, :, :), INTENT(IN)            :: rots
      INTEGER, INTENT(IN)                                :: expected

      IF (ALLOCATED(spgr%rotations_subset)) DEALLOCATE (spgr%rotations_subset)
      ALLOCATE (spgr%rotations_subset(3, 3, SIZE(rots, 3)))
      spgr%rotations_subset(:, :, :) = rots(:, :, :)
      spgr%n_operations_subset = SIZE(rots, 3)
      IF (spgr_point_group_symmetry_id(spgr, cell) /= expected) THEN
         ERROR STOP "Incorrect space-group point-group symmetry"
      END IF
   END SUBROUTINE check_pg

END PROGRAM cell_opt_utils_TEST
