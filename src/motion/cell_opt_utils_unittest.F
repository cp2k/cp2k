!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM cell_opt_utils_TEST
   USE cell_manifold_utils,             ONLY: &
        cell_manifold_active_jacobian, cell_manifold_chart_eval, cell_manifold_chart_inv, &
        cell_manifold_chart_pullback_grad, cell_manifold_geometry, cell_manifold_init, &
        cell_manifold_jacobian, cell_manifold_pack_hmat, cell_manifold_proj_pressure, &
        cell_manifold_residual, cell_manifold_status_label, cell_manifold_type, &
        cell_manifold_unpack_hmat, comp_ax, comp_bx, comp_cx, comp_cy, geom_a, geom_alpha, geom_b, &
        geom_beta, geom_c, geom_gamma, gfix_a, gfix_b, gfix_c, gfix_v, info_bad_input, &
        info_degenerate, info_licq, info_no_converge, info_ok, mask_a, mask_x, n_geom, n_slot, &
        slot_ax, slot_bx, slot_by, slot_cx, slot_cy
   USE cell_methods,                    ONLY: cell_create,&
                                              init_cell
   USE cell_opt_utils,                  ONLY: &
        cell_opt_constraint_label, get_dg_dh, gfix_alpha, gfix_beta, gfix_gamma, mask_all, &
        parse_cell_constraint, parse_token, spgr_point_group_symmetry_id, symmetry_constraint_rows
   USE cell_types,                      ONLY: &
        cell_release, cell_sym_cubic, cell_sym_hexagonal_gamma_120, cell_sym_hexagonal_gamma_60, &
        cell_sym_monoclinic, cell_sym_monoclinic_gamma_ab, cell_sym_none, cell_sym_orthorhombic, &
        cell_sym_rhombohedral, cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, &
        cell_sym_tetragonal_bc, cell_sym_triclinic, cell_type
   USE ieee_arithmetic,                 ONLY: ieee_positive_inf,&
                                              ieee_value
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi,&
                                              sqrt3
   USE space_groups_types,              ONLY: release_spgr_type,&
                                              spgr_type
   USE string_utilities,                ONLY: remove_word

   IMPLICIT NONE

   INTEGER, PARAMETER :: &
      eq_rhom(2, 4) = RESHAPE([geom_a, geom_b, geom_b, geom_c, &
                               geom_alpha, geom_beta, geom_beta, geom_gamma], [2, 4])
   REAL(KIND=dp), PARAMETER :: &
      pot_lin(n_slot) = [1.0_dp, -0.7_dp, 0.4_dp, 0.9_dp, -0.5_dp, 0.6_dp], &
      q_cubic(n_slot) = [3.5_dp, 0.0_dp, 3.5_dp, 0.0_dp, 0.0_dp, 3.5_dp], &
      q_tri(n_slot) = [5.0_dp, 1.2_dp, 4.8_dp, -0.7_dp, 1.9_dp, 5.5_dp], &
      step_dir(n_slot) = [0.13_dp, -0.27_dp, 0.11_dp, 0.31_dp, -0.07_dp, 0.19_dp]

   REAL(KIND=dp) :: ct, cy, q_rhom(n_slot), st
   TYPE(cell_type), POINTER :: cell => NULL()
   TYPE(spgr_type), POINTER :: spgr => NULL()

   ct = COS(70.0_dp*pi/180.0_dp)
   st = SIN(70.0_dp*pi/180.0_dp)
   cy = (ct - ct*ct)/st
   q_rhom = 5.1_dp*[1.0_dp, ct, st, ct, cy, SQRT(MAX(1.0_dp - ct*ct - cy*cy, 0.0_dp))]

   CALL cell_create(cell)
   ALLOCATE (spgr)
   CALL check_parser()
   CALL check_manifold_basis()
   CALL check_chart_and_pressure()
   CALL check_symmetry_and_spgr()
   CALL release_spgr_type(spgr)
   CALL cell_release(cell)

CONTAINS

! **************************************************************************************************
!> \brief Check token masks, label round trips, rejections and 2D defaults.
! **************************************************************************************************
   SUBROUTINE check_parser()
      CHARACTER(LEN=5), PARAMETER :: invalid(4) = [CHARACTER(LEN=5) :: "AY", "AZ", "BZ", "YX"], &
         legacy(7) = [CHARACTER(LEN=5) :: "NONE", "   X", "   Y", "   Z", "  XY", "  XZ", "  YZ"], &
         names(21) = [CHARACTER(LEN=5) :: "NONE", "X", "Y", "Z", "XY", "XZ", "YZ", "a", "AX", "B", &
         "BX", "BY", "C", "CX", "CY", "CZ", "ALPHA", "BETA", "GAMMA", "none", " z "]
      INTEGER, PARAMETER :: masks(21) = [0, 11, 22, 56, 31, 59, 62, 1, 1, 6, 2, 4, 56, 8, 16, 32, &
         gfix_alpha, gfix_beta, gfix_gamma, 0, 56]

      CHARACTER(LEN=default_string_length)               :: label, reason, tokens(9)
      INTEGER                                            :: i, mask, n, wanted
      LOGICAL                                            :: accepted

      ! Bits 0..5 represent AX, BX, BY, CX, CY, CZ, respectively.
      cell%perd = 1
      DO i = 1, SIZE(names)
         CALL parse_cell_constraint([names(i)], cell, .TRUE., mask)
         IF (mask /= masks(i)) ERROR STOP "Incorrect constraint token mask"
         IF (i <= SIZE(legacy)) THEN
            IF (cell_opt_constraint_label(mask) /= legacy(i) (1:4)) THEN
               ERROR STOP "Legacy constraint label changed"
            END IF
         END IF
      END DO
      CALL parse_cell_constraint( &
         [CHARACTER(LEN=5) :: "AX", "BY", "CZ", "ALPHA", "BETA", "GAMMA"], cell, .TRUE., mask)
      IF (mask /= 37 + gfix_alpha + gfix_beta + gfix_gamma .OR. &
          cell_opt_constraint_label(mask) /= "AX BY CZ ALPHA BETA GAMMA") THEN
         ERROR STOP "Angle constraint mask or label mismatch"
      END IF
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
      cell%perd = 1
   END SUBROUTINE check_parser

! **************************************************************************************************
!> \brief Check cell geometry, FD Jacobian, orthonormal [Y, Z] basis, and LICQ/degeneracy guards.
! **************************************************************************************************
   SUBROUTINE check_manifold_basis()
      INTEGER, PARAMETER :: &
         codes(6) = [info_ok, info_degenerate, info_licq, info_no_converge, info_bad_input, -99], &
         masks(7) = [gfix_gamma, gfix_alpha + gfix_beta + gfix_gamma, mask_x, mask_a + gfix_gamma, &
         gfix_a + gfix_b + gfix_c + gfix_v, mask_all, 0]

      INTEGER                                            :: c, i, info, j, pmax
      REAL(KIND=dp) :: A(n_slot, n_slot), eye(n_slot, n_slot), h, jac(n_geom, n_slot), qm(n_slot), &
         qp(n_slot), rel, s(n_geom), sm(n_geom), sp(n_geom), yz(n_slot, n_slot)
      TYPE(cell_manifold_type)                           :: m

      CALL cell_manifold_geometry(q_rhom, s, info)
      IF (info /= info_ok .OR. ANY(ABS(s(1:3) - 5.1_dp) > 1.0E-13_dp) .OR. &
          ANY(ABS(s(4:6) - 70.0_dp*pi/180.0_dp) > 1.0E-13_dp) .OR. &
          ABS(s(7) - q_rhom(1)*q_rhom(3)*q_rhom(6)) > 1.0E-12_dp) THEN
         ERROR STOP "Incorrect rhombohedral geometry"
      END IF
      CALL cell_manifold_geometry([0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], s, info)
      IF (info /= info_degenerate) ERROR STOP "Accepted zero-length cell"
      CALL cell_manifold_geometry([1.0_dp, 1.0E5_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], s, info)
      IF (info /= info_degenerate) ERROR STOP "Accepted nearly collinear cell vectors"
      CALL cell_manifold_geometry([1.0E+120_dp, 0.0_dp, 1.0E+120_dp, 0.0_dp, 0.0_dp, 1.0E+120_dp], &
                                  s, info)
      IF (info /= info_degenerate) ERROR STOP "Accepted overflow-prone cell"
      DO i = 1, SIZE(codes)
         IF (LEN_TRIM(cell_manifold_status_label(codes(i))) == 0) ERROR STOP "Empty status label"
      END DO

      CALL cell_manifold_jacobian(q_tri, jac, info)
      IF (info /= info_ok) ERROR STOP "Analytic Jacobian failed"
      rel = 0.0_dp
      DO j = 1, n_slot
         h = 1.0E-6_dp*MAX(1.0_dp, ABS(q_tri(j)))
         qp = q_tri
         qm = q_tri
         qp(j) = qp(j) + h
         qm(j) = qm(j) - h
         CALL cell_manifold_geometry(qp, sp, info)
         IF (info /= info_ok) ERROR STOP "FD +step geometry failed"
         CALL cell_manifold_geometry(qm, sm, info)
         IF (info /= info_ok) ERROR STOP "FD -step geometry failed"
         rel = MAX(rel, MAXVAL(ABS(jac(:, j) - (sp - sm)/(2.0_dp*h))/ &
                               MAX(1.0_dp, ABS((sp - sm)/(2.0_dp*h)))))
      END DO
      IF (rel > 1.0E-7_dp) ERROR STOP "Analytic Jacobian disagrees with finite differences"

      eye = 0.0_dp
      DO i = 1, n_slot
         eye(i, i) = 1.0_dp
      END DO
      DO c = 1, SIZE(masks)
         CALL cell_manifold_init(m, q_tri, masks(c), info)
         IF (info /= info_ok .OR. m%ndof /= n_slot - m%k) ERROR STOP "Chart split init failed"
         CALL cell_manifold_active_jacobian(m, q_tri, A(1:m%k, :), info)
         IF (info /= info_ok) ERROR STOP "Active Jacobian failed"
         yz(:, 1:m%k) = m%Y(:, 1:m%k)
         yz(:, m%k + 1:n_slot) = m%Z(:, 1:m%ndof)
         IF (MAXVAL(ABS(MATMUL(TRANSPOSE(yz), yz) - eye)) > 1.0E-12_dp) THEN
            ERROR STOP "[Y, Z] is not orthonormal"
         END IF
         IF (m%k > 0 .AND. m%ndof > 0) THEN
            IF (MAXVAL(ABS(MATMUL(A(1:m%k, :), m%Z(:, 1:m%ndof)))) > 1.0E-12_dp .OR. &
                MAXVAL(ABS(m%u0(1:m%ndof) - &
                           MATMUL(TRANSPOSE(m%Z(:, 1:m%ndof)), q_tri))) > 1.0E-14_dp) THEN
               ERROR STOP "Z is not in ker(A) or anchor u0 /= Z^T q0"
            END IF
         END IF
         IF (masks(c) == 0 .AND. (ANY(m%u0 /= q_tri) .OR. ANY(m%Z /= eye) .OR. &
                                  ANY(m%Y /= 0.0_dp))) THEN
            ERROR STOP "Unconstrained k = 0 basis is not canonical identity"
         END IF
         DO j = 1, n_slot
            pmax = MAXLOC(ABS(yz(:, j)), DIM=1)
            IF (yz(pmax, j) <= 0.0_dp) ERROR STOP "Basis column sign convention violated"
         END DO
      END DO

      CALL cell_manifold_init(m, q_tri, comp_ax + gfix_a, info)
      IF (info /= info_licq) ERROR STOP "Accepted redundant AX + gfix_a pair"
      CALL cell_manifold_init(m, q_tri, mask_all + gfix_v, info)
      IF (info /= info_licq) ERROR STOP "Accepted overconstrained >6-row manifold"
      CALL cell_manifold_init(m, q_cubic, gfix_a + gfix_b + gfix_c + gfix_v, info)
      IF (info /= info_licq) ERROR STOP "Missed rank collapse on cubic cell"
      CALL cell_manifold_init(m, [5.0_dp, 1.0E-6_dp, 5.0_dp, 0.0_dp, 0.0_dp, 5.0_dp], &
                              gfix_a + gfix_b + gfix_c + gfix_v, info)
      IF (info /= info_licq) ERROR STOP "Missed near rank collapse"
      CALL cell_manifold_init(m, [5.0_dp, 3.7_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp], 0, info)
      IF (info /= info_degenerate) ERROR STOP "Missed degenerate cell in init"
      CALL cell_manifold_init(m, q_tri, -1, info)
      IF (info /= info_bad_input) ERROR STOP "Accepted negative constraint mask"
      CALL cell_manifold_init(m, q_tri, 0, info, lin_w=RESHAPE([0.0_dp], [1, 1]))
      IF (info /= info_bad_input) ERROR STOP "Accepted wrong-shape lin_w"
      CALL cell_manifold_init(m, q_tri, 0, info, eq_pair=RESHAPE([geom_a, geom_a], [2, 1]))
      IF (info /= info_bad_input) ERROR STOP "Accepted self-equality pair"
   END SUBROUTINE check_manifold_basis

! **************************************************************************************************
!> \brief Check chart evaluation, curved FD pullback gradient, get_dg_dh, and projected pressure.
! **************************************************************************************************
   SUBROUTINE check_chart_and_pressure()
      REAL(KIND=dp), PARAMETER :: eye3(3, 3) = RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3]), h_fd = 1.0E-4_dp, zero3(3, 3) = 0.0_dp

      INTEGER                                            :: info, j, spec
      LOGICAL                                            :: projected
      REAL(KIND=dp) :: A4(4, n_slot), fd, g(n_slot), gu(n_slot), gz(0), lin(n_slot, 1), pdiff, &
         pint, ptens(3, 3), q(n_slot), qback(n_slot), qm(n_slot), qp(n_slot), r(4), ref, &
         u(n_slot), uback(n_slot), um(n_slot), up(n_slot), uv(n_slot), uz(0)
      TYPE(cell_manifold_type)                           :: m

      CALL cell_manifold_init(m, q_tri, 0, info)
      IF (info /= info_ok) ERROR STOP "k = 0 init failed"
      q = q_tri + 0.05_dp*step_dir
      CALL cell_manifold_chart_inv(m, q, u)
      CALL cell_manifold_chart_eval(m, u, qback, info)
      IF (info /= info_ok .OR. ANY(u /= q) .OR. ANY(qback /= q)) ERROR STOP "k = 0 chart failed"
      CALL cell_manifold_chart_pullback_grad(m, q, step_dir, gu, info)
      IF (info /= info_ok .OR. ANY(gu /= step_dir)) ERROR STOP "k = 0 pullback is not bit-identical"

      cell%hmat = -99.0_dp
      CALL cell_manifold_unpack_hmat(q_tri, cell%hmat)
      IF (cell%hmat(2, 1) /= 0.0_dp .OR. cell%hmat(3, 1) /= 0.0_dp .OR. &
          cell%hmat(3, 2) /= 0.0_dp .OR. ANY(cell_manifold_pack_hmat(cell%hmat) /= q_tri)) THEN
         ERROR STOP "hmat pack/unpack or lower-triangle zeroing failed"
      END IF
      CALL init_cell(cell)
      CALL get_dg_dh(g, 2.5_dp*eye3, 0.7_dp, cell, zero3, pint)
      CALL cell_manifold_proj_pressure(m, q_tri, g, pdiff, projected, info)
      IF (info /= info_ok .OR. .NOT. projected .OR. ABS(pint - 2.5_dp) > 1.0E-14_dp .OR. &
          MAXVAL(ABS(g([slot_bx, slot_cx, slot_cy]))) > 1.0E-13_dp .OR. &
          ABS(pdiff - 1.8_dp) > 1.0E-13_dp) THEN
         ERROR STOP "Unconstrained get_dg_dh / proj_pressure handshake failed"
      END IF
      uv = [q_tri(3)*q_tri(6), 0.0_dp, q_tri(1)*q_tri(6), 0.0_dp, 0.0_dp, q_tri(1)*q_tri(3)]
      CALL cell_manifold_proj_pressure(m, q_tri, step_dir, pdiff, projected, info)
      IF (info /= info_ok .OR. .NOT. projected .OR. &
          ABS(pdiff + DOT_PRODUCT(step_dir, uv)/DOT_PRODUCT(uv, uv)) > 1.0E-13_dp) THEN
         ERROR STOP "Unconstrained k = 0 projected pressure mismatch"
      END IF
      CALL cell_manifold_proj_pressure(m, [0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
                                       step_dir, pdiff, projected, info)
      IF (info /= info_degenerate .OR. projected) ERROR STOP "Missed degenerate cell in pressure"

      CALL cell_manifold_init(m, q_tri, mask_all, info)
      IF (info /= info_ok .OR. m%ndof /= 0) ERROR STOP "ndof = 0 init failed"
      qback = 0.0_dp
      CALL cell_manifold_chart_eval(m, uz, qback, info)
      IF (info /= info_ok .OR. ANY(qback /= m%q0)) ERROR STOP "ndof = 0 chart_eval failed"
      CALL cell_manifold_chart_inv(m, qback, uz)
      CALL cell_manifold_chart_pullback_grad(m, qback, step_dir, gz, info)
      IF (info /= info_ok) ERROR STOP "ndof = 0 pullback failed"
      CALL cell_manifold_proj_pressure(m, q_tri, step_dir, pdiff, projected, info)
      IF (info /= info_ok .OR. projected) ERROR STOP "Fully fixed cell pressure did not bypass"

      CALL cell_manifold_init(m, q_tri, gfix_beta, info)
      IF (info /= info_ok) ERROR STOP "gfix_beta init failed"
      u = 0.0_dp
      u(1) = ieee_value(1.0_dp, ieee_positive_inf)
      CALL cell_manifold_chart_eval(m, u(1:m%ndof), qback, info)
      IF (info /= info_no_converge) ERROR STOP "Accepted infinite chart coordinate"
      u(1) = 1.0E+120_dp
      CALL cell_manifold_chart_eval(m, u(1:m%ndof), qback, info)
      IF (info /= info_no_converge) ERROR STOP "Accepted 1e120 chart coordinate"
      m%max_iter = 0
      u(1:m%ndof) = m%u0(1:m%ndof) + 0.05_dp*step_dir(1:m%ndof)
      CALL cell_manifold_chart_eval(m, u(1:m%ndof), qback, info)
      IF (info /= info_no_converge) ERROR STOP "Missed Newton max_iter exhaustion"
      CALL cell_manifold_init(m, q_rhom, 0, info, eq_pair=eq_rhom)
      IF (info /= info_ok) ERROR STOP "Rhombohedral residual init failed"
      CALL cell_manifold_residual(m, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], r, info)
      IF (info /= info_degenerate) ERROR STOP "Accepted degenerate equality residual"

      lin = 0.0_dp
      lin(slot_ax, 1) = 1.0E9_dp
      lin(slot_by, 1) = -1.0E9_dp
      CALL cell_manifold_init(m, q_cubic, comp_bx + comp_cx + comp_cy, info, lin_w=lin)
      IF (info /= info_ok) ERROR STOP "Scaled linear coupling init failed"
      u(1:m%ndof) = m%u0(1:m%ndof) + 2.0_dp*step_dir(1:m%ndof)
      CALL cell_manifold_chart_eval(m, u(1:m%ndof), qback, info)
      IF (info /= info_ok .OR. m%last_iter /= 0 .OR. m%last_resid > 1.0E-12_dp) THEN
         ERROR STOP "Scaled linear coupling took nonzero Newton iterations"
      END IF

      CALL cell_manifold_init(m, q_tri, gfix_alpha + gfix_beta + gfix_gamma, info)
      IF (info /= info_ok) ERROR STOP "All-angles cone init failed"
      u(1:m%ndof) = m%u0(1:m%ndof) + 0.05_dp*step_dir(1:m%ndof)
      CALL cell_manifold_chart_eval(m, u(1:m%ndof), q, info)
      IF (info /= info_ok .OR. m%last_iter /= 0) ERROR STOP "All-angles cone took Newton steps"

      DO spec = 1, 2
         IF (spec == 1) THEN
            CALL cell_manifold_init(m, q_rhom, gfix_v, info, eq_pair=eq_rhom)
         ELSE
            CALL cell_manifold_init(m, q_tri, gfix_beta + gfix_gamma, info)
         END IF
         IF (info /= info_ok) ERROR STOP "Curved manifold init failed"
         m%tol = 1.0E-14_dp
         u(1:m%ndof) = m%u0(1:m%ndof) + 0.05_dp*step_dir(1:m%ndof)
         CALL cell_manifold_chart_eval(m, u(1:m%ndof), q, info)
         IF (info /= info_ok .OR. m%last_iter < 1 .OR. &
             m%last_resid > 1.0E-11_dp*MERGE(ABS(q(1)*q(3)*q(6)), 1.0_dp, spec == 1)) THEN
            ERROR STOP "Curved chart evaluation failed"
         END IF
         CALL cell_manifold_chart_inv(m, q, uback(1:m%ndof))
         IF (NORM2(uback(1:m%ndof) - u(1:m%ndof)) > 1.0E-13_dp) ERROR STOP "Chart round-trip failed"

         g = pot_lin + 0.4_dp*(q - m%q0)
         CALL cell_manifold_chart_pullback_grad(m, q, g, gu(1:m%ndof), info)
         IF (info /= info_ok) ERROR STOP "Reduced gradient failed"
         DO j = 1, m%ndof
            up(1:m%ndof) = u(1:m%ndof)
            um(1:m%ndof) = u(1:m%ndof)
            up(j) = up(j) + h_fd
            um(j) = um(j) - h_fd
            CALL cell_manifold_chart_eval(m, up(1:m%ndof), qp, info)
            IF (info /= info_ok) ERROR STOP "FD +u chart_eval failed"
            CALL cell_manifold_chart_eval(m, um(1:m%ndof), qm, info)
            IF (info /= info_ok) ERROR STOP "FD -u chart_eval failed"
            fd = (DOT_PRODUCT(pot_lin, qp - m%q0) + 0.2_dp*DOT_PRODUCT(qp - m%q0, qp - m%q0) - &
                  DOT_PRODUCT(pot_lin, qm - m%q0) - 0.2_dp*DOT_PRODUCT(qm - m%q0, qm - m%q0))/ &
                 (2.0_dp*h_fd)
            IF (ABS(fd - gu(j))/MAX(1.0_dp, ABS(fd)) > 1.0E-10_dp) THEN
               ERROR STOP "Reduced gradient disagrees with FD"
            END IF
         END DO
      END DO

      CALL cell_manifold_init(m, q_rhom, 0, info, eq_pair=eq_rhom)
      IF (info /= info_ok) ERROR STOP "Rhombohedral KKT init failed"
      u(1:m%ndof) = m%u0(1:m%ndof) + 0.05_dp*step_dir(1:m%ndof)
      CALL cell_manifold_chart_eval(m, u(1:m%ndof), q, info)
      IF (info /= info_ok) ERROR STOP "Rhombohedral KKT chart_eval failed"
      CALL cell_manifold_active_jacobian(m, q, A4, info)
      IF (info /= info_ok) ERROR STOP "Rhombohedral KKT active_jacobian failed"
      CALL cell_manifold_proj_pressure(m, q, &
                                       MATMUL(TRANSPOSE(A4), [0.4_dp, -1.2_dp, 0.9_dp, 0.5_dp]), &
                                       pdiff, projected, info)
      IF (info /= info_ok .OR. .NOT. projected .OR. ABS(pdiff) > 1.0E-13_dp) THEN
         ERROR STOP "Current-point KKT projected pressure is nonzero"
      END IF

      q = [5.0_dp, 0.0_dp, 6.0_dp, 0.0_dp, 0.0_dp, 7.0_dp]
      CALL cell_manifold_unpack_hmat(q, cell%hmat)
      CALL init_cell(cell)
      ptens = RESHAPE([99.0_dp, 1.3_dp, 0.0_dp, 1.3_dp, 2.7_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.7_dp], &
                      [3, 3])
      CALL get_dg_dh(g, ptens, 0.7_dp, cell, zero3, pint)
      CALL cell_manifold_init(m, q, comp_ax, info)
      IF (info /= info_ok) ERROR STOP "Pinned-AX init failed"
      ref = (2.0_dp*(q(1)*q(6))**2 + 5.0_dp*(q(1)*q(3))**2)/((q(1)*q(6))**2 + (q(1)*q(3))**2)
      CALL cell_manifold_proj_pressure(m, q, g, pdiff, projected, info)
      IF (info /= info_ok .OR. .NOT. projected .OR. &
          ABS(g(slot_bx) + 1.3_dp*q(1)*q(6)) > 1.0E-12_dp .OR. ABS(pdiff - ref) > 1.0E-12_dp*ref) &
         THEN
         ERROR STOP "Pinned-component get_dg_dh / projected pressure mismatch"
      END IF
      CALL cell_manifold_init(m, q_tri, gfix_v, info)
      IF (info /= info_ok) ERROR STOP "Constant-volume init failed"
      CALL cell_manifold_proj_pressure(m, q_tri, step_dir, pdiff, projected, info)
      IF (info /= info_ok .OR. projected) ERROR STOP "Constant-volume pressure did not bypass"
   END SUBROUTINE check_chart_and_pressure

! **************************************************************************************************
!> \brief Check crystal-system rows and space-group point-group classification.
! **************************************************************************************************
   SUBROUTINE check_symmetry_and_spgr()
      INTEGER, DIMENSION(3, 3), PARAMETER :: &
         r_c2x = RESHAPE([1, 0, 0, 0, -1, 0, 0, 0, -1], [3, 3]), &
         r_c2y = RESHAPE([-1, 0, 0, 0, 1, 0, 0, 0, -1], [3, 3]), &
         r_c3_111 = RESHAPE([0, 0, 1, 1, 0, 0, 0, 1, 0], [3, 3]), &
         r_c3_120 = RESHAPE([0, -1, 0, 1, -1, 0, 0, 0, 1], [3, 3]), &
         r_c4x = RESHAPE([1, 0, 0, 0, 0, -1, 0, 1, 0], [3, 3]), &
         r_c4y = RESHAPE([0, 0, 1, 0, 1, 0, -1, 0, 0], [3, 3]), &
         r_c4z = RESHAPE([0, -1, 0, 1, 0, 0, 0, 0, 1], [3, 3]), &
         r_c6_60 = RESHAPE([1, 1, 0, -1, 0, 0, 0, 0, 1], [3, 3]), &
         r_id = RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3]), &
         r_mz = RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, -1], [3, 3])
      INTEGER, DIMENSION(3, 3), PARAMETER :: &
         r_swap_ab = RESHAPE([0, 1, 0, 1, 0, 0, 0, 0, 1], [3, 3])
      INTEGER, PARAMETER :: ndof_exp(12) = [3, 2, 2, 2, 1, 2, 2, 2, 4, 3, 6, 6], pg_exp(12) = [ &
         cell_sym_triclinic, cell_sym_triclinic, cell_sym_monoclinic, cell_sym_orthorhombic, &
         cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, cell_sym_tetragonal_bc, cell_sym_cubic, &
         cell_sym_hexagonal_gamma_60, cell_sym_hexagonal_gamma_120, cell_sym_monoclinic_gamma_ab, &
         cell_sym_rhombohedral]
      INTEGER, PARAMETER :: sid(12) = [cell_sym_orthorhombic, cell_sym_tetragonal_ab, &
         cell_sym_tetragonal_ac, cell_sym_tetragonal_bc, cell_sym_cubic, &
         cell_sym_hexagonal_gamma_60, cell_sym_hexagonal_gamma_120, cell_sym_rhombohedral, &
         cell_sym_monoclinic, cell_sym_monoclinic_gamma_ab, cell_sym_triclinic, cell_sym_none]

      INTEGER                                            :: cand(3, 3), eq_pair(2, n_geom), &
                                                            gens(3, 3, 2, 12), i, ig, info, j, k, &
                                                            mask, n, n_eq, n_lin
      LOGICAL                                            :: curved
      REAL(KIND=dp)                                      :: lin_w(n_slot, n_slot), q0(n_slot), &
                                                            qout(n_slot), s(n_geom), u(n_slot)
      TYPE(cell_manifold_type)                           :: m

      DO i = 1, SIZE(sid)
         SELECT CASE (sid(i))
         CASE (cell_sym_hexagonal_gamma_60, cell_sym_hexagonal_gamma_120)
            q0 = [2.0_dp, MERGE(1.0_dp, -1.0_dp, sid(i) == cell_sym_hexagonal_gamma_60), &
                  sqrt3, 0.0_dp, 0.0_dp, 3.0_dp]
         CASE (cell_sym_rhombohedral)
            q0 = q_rhom
         CASE (cell_sym_monoclinic)
            q0 = [3.5_dp, 0.0_dp, 4.2_dp, 1.1_dp, 0.0_dp, 5.0_dp]
         CASE (cell_sym_monoclinic_gamma_ab)
            q0 = [5.0_dp, -1.0_dp, SQRT(24.0_dp), 0.0_dp, 0.0_dp, 6.0_dp]
         CASE DEFAULT
            q0 = q_cubic
         END SELECT

         CALL symmetry_constraint_rows(sid(i), mask, lin_w, n_lin, eq_pair, n_eq)
         CALL cell_manifold_init(m, q0, mask, info, lin_w=lin_w(:, 1:n_lin), &
                                 eq_pair=eq_pair(:, 1:n_eq))
         IF (info /= info_ok .OR. m%ndof /= ndof_exp(i)) ERROR STOP "Symmetry init failed"
         u(1:m%ndof) = m%u0(1:m%ndof) + 0.05_dp*MATMUL(TRANSPOSE(m%Z(:, 1:m%ndof)), step_dir)
         CALL cell_manifold_chart_eval(m, u(1:m%ndof), qout, info)
         IF (info /= info_ok) ERROR STOP "Symmetry chart evaluation failed"
         curved = (sid(i) == cell_sym_rhombohedral .OR. sid(i) == cell_sym_monoclinic_gamma_ab)
         IF ((curved .AND. m%last_iter < 1) .OR. (.NOT. curved .AND. m%last_iter /= 0)) THEN
            ERROR STOP "Unexpected Newton iteration count"
         END IF
         CALL cell_manifold_geometry(qout, s, info)
         IF (info /= info_ok) ERROR STOP "Symmetry chart left a degenerate cell"
         IF ((ANY(sid(i) == [cell_sym_tetragonal_ab, cell_sym_cubic, cell_sym_hexagonal_gamma_60, &
                             cell_sym_hexagonal_gamma_120, cell_sym_rhombohedral, &
                             cell_sym_monoclinic_gamma_ab]) .AND. &
              ABS(s(geom_a) - s(geom_b)) > 1.0E-11_dp*s(geom_a)) .OR. &
             (ANY(sid(i) == [cell_sym_tetragonal_ac, cell_sym_cubic, cell_sym_rhombohedral]) .AND. &
              ABS(s(geom_a) - s(geom_c)) > 1.0E-11_dp*s(geom_a)) .OR. &
             (ANY(sid(i) == [cell_sym_tetragonal_bc, cell_sym_cubic, cell_sym_rhombohedral]) .AND. &
              ABS(s(geom_b) - s(geom_c)) > 1.0E-11_dp*s(geom_b))) THEN
            ERROR STOP "Symmetry step broke equal cell lengths"
         END IF
         IF ((ANY(sid(i) == [cell_sym_orthorhombic, cell_sym_tetragonal_ab, &
                             cell_sym_tetragonal_ac, cell_sym_tetragonal_bc, cell_sym_cubic]) &
              .AND. MAXVAL(ABS(s(geom_alpha:geom_gamma) - 0.5_dp*pi)) > 1.0E-12_dp) .OR. &
             (sid(i) == cell_sym_monoclinic .AND. &
              MAX(ABS(s(geom_alpha) - 0.5_dp*pi), ABS(s(geom_gamma) - 0.5_dp*pi)) > 1.0E-12_dp) &
             .OR. (ANY(sid(i) == [cell_sym_monoclinic_gamma_ab, cell_sym_hexagonal_gamma_60, &
                                  cell_sym_hexagonal_gamma_120]) .AND. &
                   MAX(ABS(s(geom_alpha) - 0.5_dp*pi), ABS(s(geom_beta) - 0.5_dp*pi)) > &
                   1.0E-12_dp) .OR. &
             (sid(i) == cell_sym_hexagonal_gamma_60 .AND. ABS(s(geom_gamma) - pi/3.0_dp) > &
              1.0E-12_dp) .OR. &
             (sid(i) == cell_sym_hexagonal_gamma_120 .AND. ABS(s(geom_gamma) - 2.0_dp*pi/3.0_dp) > &
              1.0E-12_dp) .OR. &
             (sid(i) == cell_sym_rhombohedral .AND. &
              MAX(ABS(s(geom_alpha) - s(geom_beta)), ABS(s(geom_beta) - s(geom_gamma))) > &
              1.0E-11_dp)) THEN
            ERROR STOP "Symmetry step broke crystal angles"
         END IF
      END DO

      CALL cell_manifold_unpack_hmat(q_tri, cell%hmat)
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_triclinic) THEN
         ERROR STOP "Empty rotation set did not return cell_sym_triclinic"
      END IF
      gens = SPREAD(SPREAD(r_id, DIM=3, NCOPIES=2), DIM=4, NCOPIES=12)
      gens(:, :, 1, 2:12) = RESHAPE([r_c2x, r_c2y, r_c2x, r_c4z, r_c4y, r_c4x, r_c2x, &
                                     r_c6_60, r_c3_120, r_swap_ab, r_c3_111], [3, 3, 11])
      gens(:, :, 2, 4) = r_c2y
      gens(:, :, 2, 8) = r_c3_111
      gens(:, :, 2, 11) = r_mz

      ALLOCATE (spgr%rotations_subset(3, 3, 12))
      DO i = 1, SIZE(pg_exp)
         spgr%rotations_subset(:, :, 1) = r_id
         n = 1
         k = 1
         DO WHILE (k <= n)
            DO ig = 1, 2
               cand = MATMUL(spgr%rotations_subset(:, :, k), gens(:, :, ig, i))
               IF (.NOT. ANY([(ALL(spgr%rotations_subset(:, :, j) == cand), j=1, n)])) THEN
                  n = n + 1
                  spgr%rotations_subset(:, :, n) = cand
               END IF
            END DO
            k = k + 1
         END DO
         spgr%n_operations_subset = n
         IF (spgr_point_group_symmetry_id(spgr, cell) /= pg_exp(i)) THEN
            ERROR STOP "Incorrect space-group point-group symmetry"
         END IF
      END DO
      cell%hmat = 0.0_dp
      IF (spgr_point_group_symmetry_id(spgr, cell) /= cell_sym_triclinic) THEN
         ERROR STOP "Zero cell did not return cell_sym_triclinic"
      END IF
   END SUBROUTINE check_symmetry_and_spgr

END PROGRAM cell_opt_utils_TEST
