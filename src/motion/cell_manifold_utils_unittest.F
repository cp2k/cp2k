!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM cell_manifold_utils_TEST
   USE cell_manifold_utils,             ONLY: &
        cell_manifold_active_jacobian, cell_manifold_chart_eval, cell_manifold_chart_inv, &
        cell_manifold_chart_pullback_grad, cell_manifold_geometry, cell_manifold_init, &
        cell_manifold_jacobian, cell_manifold_pack_hmat, cell_manifold_proj_pressure, &
        cell_manifold_residual, cell_manifold_status_label, cell_manifold_type, &
        cell_manifold_unpack_hmat, comp_ax, comp_bx, comp_cx, comp_cy, geom_a, geom_alpha, geom_b, &
        geom_beta, geom_c, geom_gamma, gfix_a, gfix_alpha, gfix_b, gfix_beta, gfix_c, gfix_gamma, &
        gfix_v, info_bad_input, info_degenerate, info_licq, info_no_converge, info_ok, mask_a, &
        mask_all, mask_x, n_geom, n_slot, slot_ax, slot_by
   USE ieee_arithmetic,                 ONLY: ieee_positive_inf,&
                                              ieee_value
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi

   IMPLICIT NONE

   INTEGER, PARAMETER :: &
      eq_rhom(2, 4) = RESHAPE([geom_a, geom_b, geom_b, geom_c, &
                               geom_alpha, geom_beta, geom_beta, geom_gamma], [2, 4])
   REAL(KIND=dp), PARAMETER :: &
      pot_lin(n_slot) = [1.0_dp, -0.7_dp, 0.4_dp, 0.9_dp, -0.5_dp, 0.6_dp], &
      q_ortho(n_slot) = [3.5_dp, 0.0_dp, 3.5_dp, 0.0_dp, 0.0_dp, 3.5_dp], &
      q_tri(n_slot) = [5.0_dp, 1.2_dp, 4.8_dp, -0.7_dp, 1.9_dp, 5.5_dp], &
      step_dir(n_slot) = [0.13_dp, -0.27_dp, 0.11_dp, 0.31_dp, -0.07_dp, 0.19_dp]

   REAL(KIND=dp) :: ct, cx, cy, cz, q_mono(n_slot), q_rhom(n_slot), st

   q_mono = [5.1_dp, 0.0_dp, 4.9_dp, 6.2_dp*COS(105.0_dp*pi/180.0_dp), &
             0.0_dp, 6.2_dp*SIN(105.0_dp*pi/180.0_dp)]
   ct = COS(70.0_dp*pi/180.0_dp)
   st = SIN(70.0_dp*pi/180.0_dp)
   cx = ct
   cy = (ct - ct*ct)/st
   cz = SQRT(MAX(1.0_dp - cx*cx - cy*cy, 0.0_dp))
   q_rhom = 5.1_dp*[1.0_dp, ct, st, cx, cy, cz]

   CALL check_geometry_and_jacobian()
   CALL check_basis_and_licq()
   CALL check_chart_and_pullback()
   CALL check_proj_pressure()

CONTAINS

! **************************************************************************************************
!> \brief Check cell geometry, degeneracy/overflow guards, helpers, and FD Jacobians.
! **************************************************************************************************
   SUBROUTINE check_geometry_and_jacobian()
      INTEGER, PARAMETER :: &
         codes(5) = [info_ok, info_degenerate, info_licq, info_no_converge, info_bad_input]

      INTEGER                                            :: c, i, info, j
      REAL(KIND=dp) :: cells(n_slot, 2), h, hmat(3, 3), jac(n_geom, n_slot), qm(n_slot), &
         qp(n_slot), rel, s(n_geom), sm(n_geom), sp(n_geom)

      CALL cell_manifold_geometry(q_mono, s, info)
      IF (info /= info_ok) ERROR STOP "Monoclinic geometry failed"
      IF (ANY(ABS(s(1:3) - [5.1_dp, 4.9_dp, 6.2_dp]) > 1.0E-14_dp) .OR. &
          ABS(s(4) - 0.5_dp*pi) > 1.0E-12_dp .OR. ABS(s(5) - 105.0_dp*pi/180.0_dp) > 1.0E-12_dp &
          .OR. ABS(s(6) - 0.5_dp*pi) > 1.0E-12_dp .OR. ABS(s(7) - q_mono(1)*q_mono(3)*q_mono(6)) > &
          1.0E-12_dp) THEN
         ERROR STOP "Incorrect monoclinic lengths, angles, or volume"
      END IF
      CALL cell_manifold_geometry(q_rhom, s, info)
      IF (info /= info_ok .OR. ANY(ABS(s(1:3) - 5.1_dp) > 1.0E-13_dp) .OR. &
          ANY(ABS(s(4:6) - 70.0_dp*pi/180.0_dp) > 1.0E-13_dp)) THEN
         ERROR STOP "Incorrect rhombohedral lengths or skewed angles"
      END IF

      CALL cell_manifold_geometry([0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], s, info)
      IF (info /= info_degenerate) ERROR STOP "Accepted zero-length cell"
      CALL cell_manifold_geometry([1.0_dp, 1.0E5_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], s, info)
      IF (info /= info_degenerate) ERROR STOP "Accepted nearly collinear cell vectors"
      CALL cell_manifold_geometry([1.0E+120_dp, 0.0_dp, 1.0E+120_dp, 0.0_dp, 0.0_dp, 1.0E+120_dp], &
                                  s, info)
      IF (info /= info_degenerate) ERROR STOP "Accepted overflow-prone cell"

      hmat = -99.0_dp
      CALL cell_manifold_unpack_hmat(q_tri, hmat)
      IF (hmat(2, 1) /= 0.0_dp .OR. hmat(3, 1) /= 0.0_dp .OR. hmat(3, 2) /= 0.0_dp .OR. &
          ANY(cell_manifold_pack_hmat(hmat) /= q_tri)) THEN
         ERROR STOP "hmat pack/unpack round-trip failed"
      END IF
      DO i = 1, SIZE(codes)
         IF (LEN_TRIM(cell_manifold_status_label(codes(i))) == 0) ERROR STOP "Empty status label"
      END DO

      cells(:, 1) = q_tri
      cells(:, 2) = q_ortho
      rel = 0.0_dp
      DO c = 1, SIZE(cells, 2)
         CALL cell_manifold_jacobian(cells(:, c), jac, info)
         IF (info /= info_ok) ERROR STOP "Analytic Jacobian failed"
         DO j = 1, n_slot
            h = 1.0E-6_dp*MAX(1.0_dp, ABS(cells(j, c)))
            qp = cells(:, c)
            qm = cells(:, c)
            qp(j) = qp(j) + h
            qm(j) = qm(j) - h
            CALL cell_manifold_geometry(qp, sp, info)
            IF (info /= info_ok) ERROR STOP "FD +step geometry failed"
            CALL cell_manifold_geometry(qm, sm, info)
            IF (info /= info_ok) ERROR STOP "FD -step geometry failed"
            DO i = 1, n_geom
               rel = MAX(rel, ABS(jac(i, j) - (sp(i) - sm(i))/(2.0_dp*h))/ &
                         MAX(1.0_dp, ABS((sp(i) - sm(i))/(2.0_dp*h))))
            END DO
         END DO
      END DO
      IF (rel > 1.0E-7_dp) ERROR STOP "Analytic Jacobian disagrees with finite differences"
   END SUBROUTINE check_geometry_and_jacobian

! **************************************************************************************************
!> \brief Check orthonormal [Y, Z] split, sign convention, k = 0 identity, and rank checks.
! **************************************************************************************************
   SUBROUTINE check_basis_and_licq()
      INTEGER, PARAMETER :: masks(6) = [gfix_gamma, gfix_alpha + gfix_beta + gfix_gamma, mask_x, &
         mask_a + gfix_gamma, gfix_a + gfix_b + gfix_c + gfix_v, 0]

      INTEGER                                            :: i, info, j, mi, pmax
      REAL(KIND=dp)                                      :: A(n_slot, n_slot), eye(n_slot, n_slot)
      TYPE(cell_manifold_type)                           :: m

      eye = 0.0_dp
      DO i = 1, n_slot
         eye(i, i) = 1.0_dp
      END DO

      DO mi = 1, SIZE(masks)
         CALL cell_manifold_init(m, q_tri, masks(mi), info)
         IF (info /= info_ok .OR. m%ndof /= n_slot - m%k) ERROR STOP "Chart split init failed"
         CALL cell_manifold_active_jacobian(m, q_tri, A(1:m%k, :), info)
         IF (info /= info_ok) ERROR STOP "Active Jacobian failed"
         IF (m%k > 0 .AND. m%ndof > 0) THEN
            IF (MAXVAL(ABS(MATMUL(A(1:m%k, :), m%Z(:, 1:m%ndof)))) > 1.0E-12_dp) THEN
               ERROR STOP "Z is not in ker(A)"
            END IF
            IF (MAXVAL(ABS(MATMUL(TRANSPOSE(m%Y(:, 1:m%k)), m%Z(:, 1:m%ndof)))) > 1.0E-12_dp) THEN
               ERROR STOP "Y and Z are not orthogonal"
            END IF
         END IF
         IF (m%k > 0) THEN
            IF (MAXVAL(ABS(MATMUL(TRANSPOSE(m%Y(:, 1:m%k)), m%Y(:, 1:m%k)) - &
                           eye(1:m%k, 1:m%k))) > 1.0E-12_dp) THEN
               ERROR STOP "Y is not orthonormal"
            END IF
            IF (ANY(m%u0(1:m%ndof) /= MATMUL(TRANSPOSE(m%Z(:, 1:m%ndof)), q_tri))) THEN
               ERROR STOP "Anchor u0 /= Z^T q0"
            END IF
         END IF
         IF (m%ndof > 0) THEN
            IF (MAXVAL(ABS(MATMUL(TRANSPOSE(m%Z(:, 1:m%ndof)), m%Z(:, 1:m%ndof)) - &
                           eye(1:m%ndof, 1:m%ndof))) > 1.0E-12_dp) THEN
               ERROR STOP "Z is not orthonormal"
            END IF
         END IF
         DO j = 1, m%k
            pmax = MAXLOC(ABS(m%Y(:, j)), DIM=1)
            IF (m%Y(pmax, j) <= 0.0_dp) ERROR STOP "Y column sign convention violated"
         END DO
         DO j = 1, m%ndof
            pmax = MAXLOC(ABS(m%Z(:, j)), DIM=1)
            IF (m%Z(pmax, j) <= 0.0_dp) ERROR STOP "Z column sign convention violated"
         END DO
      END DO

      IF (ANY(m%u0 /= q_tri) .OR. ANY(m%Z /= eye) .OR. ANY(m%Y /= 0.0_dp)) THEN
         ERROR STOP "Unconstrained k = 0 basis is not canonical identity"
      END IF

      CALL cell_manifold_init(m, q_tri, comp_ax + gfix_a, info)
      IF (info /= info_licq) ERROR STOP "Accepted redundant AX + gfix_a pair"
      CALL cell_manifold_init(m, q_tri, -1, info)
      IF (info /= info_bad_input) ERROR STOP "Accepted negative constraint mask"
      CALL cell_manifold_init(m, q_tri, 0, info, lin_w=RESHAPE([0.0_dp], [1, 1]))
      IF (info /= info_bad_input) ERROR STOP "Accepted wrong-shape lin_w"
      CALL cell_manifold_init(m, q_tri, 0, info, eq_pair=RESHAPE([geom_a, geom_a], [2, 1]))
      IF (info /= info_bad_input) ERROR STOP "Accepted self-equality pair"

      CALL cell_manifold_init(m, [5.0_dp, 0.0_dp, 5.0_dp, 0.0_dp, 0.0_dp, 5.0_dp], &
                              gfix_a + gfix_b + gfix_c + gfix_v, info)
      IF (info /= info_licq) ERROR STOP "Missed exact rank collapse"
      CALL cell_manifold_init(m, [5.0_dp, 1.0E-6_dp, 5.0_dp, 0.0_dp, 0.0_dp, 5.0_dp], &
                              gfix_a + gfix_b + gfix_c + gfix_v, info)
      IF (info /= info_licq) ERROR STOP "Missed near rank collapse"
      CALL cell_manifold_init(m, [5.0_dp, 3.7_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp], &
                              gfix_a + gfix_b + gfix_c + gfix_v, info)
      IF (info /= info_degenerate) ERROR STOP "Missed degenerate cell in init"
   END SUBROUTINE check_basis_and_licq

! **************************************************************************************************
!> \brief Check chart edge cases, flat 0-step property, curved round-trip, and FD gradient.
! **************************************************************************************************
   SUBROUTINE check_chart_and_pullback()
      REAL(KIND=dp), PARAMETER                           :: h_fd = 1.0E-4_dp

      INTEGER                                            :: info, j, spec
      REAL(KIND=dp) :: derr, fd, g(n_slot), gback(n_slot), gu(n_slot), gz(0), lin(n_slot, n_slot), &
         q(n_slot), qback(n_slot), qm(n_slot), qp(n_slot), qpre(n_slot), r(n_slot), u(n_slot), &
         uback(n_slot), um(n_slot), up(n_slot), uz(0)
      TYPE(cell_manifold_type)                           :: m

      CALL cell_manifold_init(m, q_tri, 0, info)
      CALL cell_manifold_chart_inv(m, q_tri, u)
      CALL cell_manifold_chart_eval(m, u, qback, info)
      CALL cell_manifold_chart_pullback_grad(m, q_tri, step_dir, gback, info)
      IF (info /= info_ok .OR. ANY(u /= q_tri) .OR. ANY(qback /= q_tri) .OR. &
          ANY(gback /= step_dir)) ERROR STOP "k = 0 chart is not bit-identical"

      CALL cell_manifold_init(m, q_tri, mask_all, info)
      CALL cell_manifold_chart_eval(m, uz, qback, info)
      CALL cell_manifold_chart_inv(m, qback, uz)
      CALL cell_manifold_chart_pullback_grad(m, qback, step_dir, gz, info)
      IF (info /= info_ok .OR. m%ndof /= 0 .OR. ANY(qback /= m%q0)) THEN
         ERROR STOP "ndof = 0 chart boundary failed"
      END IF

      CALL cell_manifold_init(m, q_tri, gfix_beta, info)
      u = 0.0_dp
      u(1) = ieee_value(1.0_dp, ieee_positive_inf)
      CALL cell_manifold_chart_eval(m, u(1:m%ndof), qback, info)
      IF (info /= info_no_converge) ERROR STOP "Accepted infinite chart coordinate"
      u(1) = 1.0E+120_dp
      CALL cell_manifold_chart_eval(m, u(1:m%ndof), qback, info)
      IF (info /= info_no_converge) ERROR STOP "Accepted 1e120 chart coordinate"
      CALL cell_manifold_init(m, q_rhom, 0, info, eq_pair=eq_rhom)
      CALL cell_manifold_residual(m, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], r(1:4), info)
      IF (info /= info_degenerate) ERROR STOP "Accepted degenerate equality residual"

      ! 1e9 weight scaling tests unit-row normalization on a flat coupling.
      lin = 0.0_dp
      lin(slot_ax, 1) = 1.0E9_dp
      lin(slot_by, 1) = -1.0E9_dp
      CALL cell_manifold_init(m, [5.1_dp, 0.0_dp, 5.1_dp, 0.0_dp, 0.0_dp, 6.3_dp], &
                              comp_bx + comp_cx + comp_cy, info, lin_w=lin(:, 1:1))
      u(1:m%ndof) = m%u0(1:m%ndof) + 2.0_dp*step_dir(1:m%ndof)
      CALL cell_manifold_chart_eval(m, u(1:m%ndof), qback, info)
      IF (info /= info_ok .OR. m%last_iter /= 0 .OR. m%last_resid > 1.0E-12_dp) THEN
         ERROR STOP "Scaled linear coupling took nonzero Newton iterations"
      END IF

      DO spec = 1, 3
         SELECT CASE (spec)
         CASE (1)
            CALL cell_manifold_init(m, q_tri, gfix_alpha + gfix_beta + gfix_gamma, info)
         CASE (2)
            CALL cell_manifold_init(m, q_rhom, gfix_v, info, eq_pair=eq_rhom)
         CASE (3)
            CALL cell_manifold_init(m, q_tri, gfix_beta + gfix_gamma, info)
         CASE DEFAULT
            ERROR STOP "Invalid constraint test spec"
         END SELECT
         IF (info /= info_ok) ERROR STOP "Constraint init failed in chart loop"

         IF (spec == 3) THEN
            u(1:m%ndof) = m%u0(1:m%ndof) + MATMUL(TRANSPOSE(m%Z(:, 1:m%ndof)), step_dir)
            qpre = m%q0 + MATMUL(m%Z(:, 1:m%ndof), u(1:m%ndof) - m%u0(1:m%ndof))
            CALL cell_manifold_residual(m, qpre, r(1:m%k), info)
            IF (info /= info_ok .OR. NORM2(r(1:m%k)) <= 1.0E-8_dp) THEN
               ERROR STOP "Curved constraint showed no tangent drift"
            END IF
         ELSE
            u(1:m%ndof) = m%u0(1:m%ndof) + 0.05_dp*step_dir(1:m%ndof)
         END IF

         CALL cell_manifold_chart_eval(m, u(1:m%ndof), q, info)
         IF (info /= info_ok) ERROR STOP "Chart evaluation failed"
         IF (spec == 1 .AND. m%last_iter /= 0) THEN
            ERROR STOP "All-angles cone took nonzero Newton iterations"
         END IF
         IF (spec /= 1 .AND. m%last_iter < 1) THEN
            ERROR STOP "Curved constraint took 0 Newton iterations"
         END IF
         IF (m%last_resid > 1.0E-11_dp*MERGE(ABS(q(1)*q(3)*q(6)), 1.0_dp, spec == 2)) THEN
            ERROR STOP "Constraint residual exceeds tolerance"
         END IF

         CALL cell_manifold_chart_inv(m, q, uback(1:m%ndof))
         IF (NORM2(uback(1:m%ndof) - u(1:m%ndof)) > 1.0E-13_dp) ERROR STOP "Chart round-trip failed"
         CALL cell_manifold_chart_eval(m, m%u0(1:m%ndof), qback, info)
         CALL cell_manifold_chart_inv(m, qback, uback(1:m%ndof))
         IF (info /= info_ok .OR. ANY(uback(1:m%ndof) /= m%u0(1:m%ndof))) THEN
            ERROR STOP "Anchor round-trip is not exact"
         END IF

         m%tol = 1.0E-14_dp
         CALL cell_manifold_chart_eval(m, u(1:m%ndof), q, info)
         IF (info /= info_ok) ERROR STOP "Tight-tolerance chart evaluation failed"
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
            fd = (testpot(qp, m%q0) - testpot(qm, m%q0))/(2.0_dp*h_fd)
            derr = ABS(fd - gu(j))/MAX(1.0_dp, ABS(fd))
            IF (derr > 1.0E-10_dp) ERROR STOP "Reduced gradient disagrees with FD"
         END DO
      END DO
   END SUBROUTINE check_chart_and_pullback

! **************************************************************************************************
!> \brief Check KKT zeros, component/weighted formulas, bypasses, and errors.
! **************************************************************************************************
   SUBROUTINE check_proj_pressure()
      INTEGER                                            :: info
      LOGICAL                                            :: projected
      REAL(KIND=dp)                                      :: A4(4, n_slot), g(n_slot), pdiff, &
                                                            q(n_slot), ref, u(n_slot), uv(n_slot)
      TYPE(cell_manifold_type)                           :: m

      CALL cell_manifold_init(m, q_rhom, 0, info, eq_pair=eq_rhom)
      u(1:m%ndof) = m%u0(1:m%ndof) + 0.05_dp*step_dir(1:m%ndof)
      CALL cell_manifold_chart_eval(m, u(1:m%ndof), q, info)
      IF (info /= info_ok) ERROR STOP "Rhombohedral chart evaluation failed"
      CALL cell_manifold_active_jacobian(m, q, A4, info)
      IF (info /= info_ok) ERROR STOP "Rhombohedral active Jacobian failed"
      CALL cell_manifold_proj_pressure(m, q, &
                                       MATMUL(TRANSPOSE(A4), [0.4_dp, -1.2_dp, 0.9_dp, 0.5_dp]), &
                                       pdiff, projected, info)
      IF (info /= info_ok .OR. .NOT. projected .OR. ABS(pdiff) > 1.0E-13_dp) THEN
         ERROR STOP "Current-point KKT projected pressure is nonzero"
      END IF

      q = [5.0_dp, 0.0_dp, 6.0_dp, 0.0_dp, 0.0_dp, 7.0_dp]
      CALL cell_manifold_init(m, q, comp_ax, info)
      g = [99.0_dp, 0.0_dp, -2.0_dp*q(1)*q(6), 0.0_dp, 0.0_dp, -5.0_dp*q(1)*q(3)]
      ref = (2.0_dp*(q(1)*q(6))**2 + 5.0_dp*(q(1)*q(3))**2)/((q(1)*q(6))**2 + (q(1)*q(3))**2)
      CALL cell_manifold_proj_pressure(m, q, g, pdiff, projected, info)
      IF (info /= info_ok .OR. .NOT. projected .OR. ABS(pdiff - ref) > 1.0E-12_dp*ref) THEN
         ERROR STOP "Multi-component weighted projected pressure mismatch"
      END IF

      CALL cell_manifold_init(m, q_mono, 0, info)
      uv = [q_mono(3)*q_mono(6), 0.0_dp, q_mono(1)*q_mono(6), 0.0_dp, 0.0_dp, q_mono(1)*q_mono(3)]
      ref = -DOT_PRODUCT(step_dir, uv)/DOT_PRODUCT(uv, uv)
      CALL cell_manifold_proj_pressure(m, q_mono, step_dir, pdiff, projected, info)
      IF (info /= info_ok .OR. .NOT. projected .OR. ABS(pdiff - ref) > 1.0E-13_dp) THEN
         ERROR STOP "Unconstrained projected pressure mismatch"
      END IF

      CALL cell_manifold_init(m, q_mono, gfix_v, info)
      CALL cell_manifold_proj_pressure(m, q_mono, step_dir, pdiff, projected, info)
      IF (info /= info_ok .OR. projected) ERROR STOP "Constant-volume pressure did not bypass"
      CALL cell_manifold_init(m, q_mono, mask_all, info)
      CALL cell_manifold_proj_pressure(m, q_mono, step_dir, pdiff, projected, info)
      IF (info /= info_ok .OR. projected) ERROR STOP "Fully fixed cell pressure did not bypass"
      CALL cell_manifold_init(m, q_mono, 0, info)
      CALL cell_manifold_proj_pressure(m, [0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
                                       step_dir, pdiff, projected, info)
      IF (info /= info_degenerate .OR. projected) THEN
         ERROR STOP "Missed degenerate cell in proj_pressure"
      END IF
   END SUBROUTINE check_proj_pressure

! **************************************************************************************************
!> \brief Smooth test potential E(q) = c . (q - qref) + 0.2 * |q - qref|^2 for FD pullback checks.
!> \param q canonical cell coordinates
!> \param qref anchor cell coordinates
!> \return potential value
! **************************************************************************************************
   PURE FUNCTION testpot(q, qref) RESULT(e)
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot), qref(n_slot)
      REAL(KIND=dp)                                      :: e

      REAL(KIND=dp)                                      :: dq(n_slot)

      dq = q - qref
      e = DOT_PRODUCT(pot_lin, dq) + 0.2_dp*DOT_PRODUCT(dq, dq)
   END FUNCTION testpot

END PROGRAM cell_manifold_utils_TEST
