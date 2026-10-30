!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Reduced-coordinate parameterization of cell constraints for CELL_OPT.
!>
!> This pure-math module depends only on kinds, ieee_arithmetic, and LAPACK.
!> It has no force_env, cell_type, or input dependencies.
!>
!> Canonical cell vector:
!>   CP2K places lattice vector a along x and b in the xy plane:
!>     hmat(:, 1) = (AX,  0,  0)
!>     hmat(:, 2) = (BX, BY,  0)
!>     hmat(:, 3) = (CX, CY, CZ)
!>   The 6 upper-triangular entries form q = (AX, BX, BY, CX, CY, CZ).
!>   The 7 geometric quantities are s(q) = (a, b, c, alpha, beta, gamma, V),
!>   with angles in radians and volume V = AX * BY * CZ.
!>
!> Constraint equations and reduced coordinates:
!>   Active constraints define k equations r(q) = 0 with k x 6 Jacobian A(q) = dr/dq:
!>     1. Fixed cell components: q_i - q0_i = 0.
!>     2. Fixed geometric quantities: s_m(q) - s_m(q0) = 0.
!>     3. Equal lengths or angles: s_m(q) - s_n(q) = 0 (for example a = b).
!>     4. Linear couplings: w^T (q - q0) = 0 (for example AX = BY).
!>
!>   At the initial cell q0, QR factorization splits R^6 into two fixed
!>   orthonormal bases:
!>     - Y (6 x k) spans the constrained directions row(A(q0)).
!>     - Z (6 x ndof, ndof = 6 - k) spans the free directions ker(A(q0)).
!>
!>   The optimizer works only on the reduced vector u of length ndof:
!>     - Forward map (chart_eval): q(u) = q0 + Z (u - u0) + Y w(u), where
!>       Newton's method solves the k x k system r(q(u)) = 0 for w(u).
!>       Linear constraints satisfy r = 0 in zero Newton iterations.
!>     - Inverse map (chart_inv): u = u0 + Z^T (q - q0), because Z^T Y = 0.
!>     - Reduced gradient (chart_pullback_grad): g_u = (dq/du)^T g_q.
!>       Implicit differentiation gives g_u = Z^T (g_q - A(q)^T lam) with
!>       (A(q) Y)^T lam = Y^T g_q.
!>     - Unconstrained case (k = 0): Z = I and u0 = q0, so u = q and g_u = g_q
!>       with zero floating-point roundoff.
!>
!>   Sign convention: get_dg_dh returns g_q = +d(E + pres_ext*V)/dq, and
!>   cell_manifold_proj_pressure maps g_q to pdiff = P_constr - pres_ext.
! **************************************************************************************************
MODULE cell_manifold_utils

   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: default_string_length,&
                                              dp

   IMPLICIT NONE
   PRIVATE

   ! Indices of q = (AX, BX, BY, CX, CY, CZ) and s = (a, b, c, alpha, beta, gamma, V).
   INTEGER, PARAMETER, PUBLIC :: n_slot = 6, n_geom = 7
   INTEGER, PARAMETER, PUBLIC :: slot_ax = 1, &
                                 slot_bx = 2, &
                                 slot_by = 3, &
                                 slot_cx = 4, &
                                 slot_cy = 5, &
                                 slot_cz = 6
   INTEGER, PARAMETER, PUBLIC :: geom_a = 1, &
                                 geom_b = 2, &
                                 geom_c = 3, &
                                 geom_alpha = 4, &
                                 geom_beta = 5, &
                                 geom_gamma = 6, &
                                 geom_v = 7

   ! Bit flags for fixed canonical components and direction/vector shorthands.
   INTEGER, PARAMETER, PUBLIC :: comp_ax = 2**0, &
                                 comp_bx = 2**1, &
                                 comp_by = 2**2, &
                                 comp_cx = 2**3, &
                                 comp_cy = 2**4, &
                                 comp_cz = 2**5
   INTEGER, PARAMETER, PUBLIC :: mask_a = comp_ax, &
                                 mask_b = comp_bx + comp_by, &
                                 mask_c = comp_cx + comp_cy + comp_cz, &
                                 mask_x = comp_ax + comp_bx + comp_cx, &
                                 mask_y = comp_bx + comp_by + comp_cy, &
                                 mask_z = mask_c, &
                                 mask_xy = comp_ax + comp_bx + comp_by + comp_cx + comp_cy, &
                                 mask_xz = comp_ax + comp_bx + comp_cx + comp_cy + comp_cz, &
                                 mask_yz = comp_bx + comp_by + comp_cx + comp_cy + comp_cz, &
                                 mask_all = comp_ax + comp_bx + comp_by + comp_cx + &
                                 comp_cy + comp_cz

   ! Bit flags that fix s(geom_*) to its initial value s(q0).
   ! cell_opt_types sets these internal bits for KEEP_ANGLES and KEEP_VOLUME.
   ! The CONSTRAINT parser does not expose them.
   INTEGER, PARAMETER, PUBLIC :: gfix_a = 2**6, &
                                 gfix_b = 2**7, &
                                 gfix_c = 2**8, &
                                 gfix_alpha = 2**9, &
                                 gfix_beta = 2**10, &
                                 gfix_gamma = 2**11, &
                                 gfix_v = 2**12
   INTEGER, PARAMETER, PUBLIC :: gfix_all = gfix_a + gfix_b + gfix_c + gfix_alpha + &
                                 gfix_beta + gfix_gamma + gfix_v
   INTEGER, PARAMETER, PUBLIC :: n_constraint_bits = n_slot + n_geom

   ! max_cell_coord keeps AX*BY*CZ well below HUGE(1.0_dp) (~1.8E308) so volume
   ! evaluation cannot trap under -ffpe-trap=overflow.
   REAL(KIND=dp), PARAMETER, PUBLIC  :: min_cell_length = 1.0E-8_dp
   REAL(KIND=dp), PARAMETER, PRIVATE :: max_cell_coord = 1.0E90_dp
   REAL(KIND=dp), PARAMETER, PRIVATE :: collinear_cos_margin = 1.0E-8_dp
   REAL(KIND=dp), PARAMETER, PRIVATE :: licq_pivot_tol = 1.0E-6_dp
   REAL(KIND=dp), PARAMETER, PRIVATE :: default_chart_tol = 1.0E-12_dp
   INTEGER, PARAMETER, PRIVATE       :: default_chart_max_iter = 100
   REAL(KIND=dp), PARAMETER, PRIVATE :: proj_pressure_rel_tol = 1.0E-12_dp

   ! Status codes (negative values indicate an error).
   INTEGER, PARAMETER, PUBLIC :: info_ok = 0
   INTEGER, PARAMETER, PUBLIC :: info_degenerate = -1
   INTEGER, PARAMETER, PUBLIC :: info_licq = -2
   INTEGER, PARAMETER, PUBLIC :: info_no_converge = -3
   INTEGER, PARAMETER, PUBLIC :: info_bad_input = -4

   !> \brief Active constraint set and fixed coordinate bases anchored at q0.
   TYPE, PUBLIC :: cell_manifold_type
      INTEGER              :: mask = 0
      REAL(KIND=dp)        :: q0(n_slot) = 0.0_dp
      REAL(KIND=dp)        :: s0(n_geom) = 0.0_dp
      INTEGER              :: k = 0
      INTEGER              :: ndof = n_slot
      INTEGER              :: slot_rows(n_slot) = 0
      INTEGER              :: n_slot_rows = 0
      INTEGER              :: geom_rows(n_geom) = 0
      INTEGER              :: n_geom_rows = 0
      ! Column i holds the unit normal w_i of the linear constraint w_i^T (q - q0) = 0.
      REAL(KIND=dp)        :: lin_w(n_slot, n_slot) = 0.0_dp
      INTEGER              :: n_lin_rows = 0
      ! Column i holds a pair [m, n] in 1..6 for the equality s_m(q) - s_n(q) = 0.
      INTEGER              :: eq_pair(2, n_geom) = 0
      INTEGER              :: n_eq_rows = 0
      ! Fixed bases set once at init: Y(:, 1:k) spans row(A(q0)),
      ! Z(:, 1:ndof) spans ker(A(q0)), and u0(1:ndof) = Z(:, 1:ndof)^T q0.
      REAL(KIND=dp)        :: Y(n_slot, n_slot) = 0.0_dp
      REAL(KIND=dp)        :: Z(n_slot, n_slot) = 0.0_dp
      REAL(KIND=dp)        :: u0(n_slot) = 0.0_dp
      REAL(KIND=dp)        :: tol = default_chart_tol
      INTEGER              :: max_iter = default_chart_max_iter
      INTEGER              :: last_iter = 0
      REAL(KIND=dp)        :: last_resid = 0.0_dp
   END TYPE cell_manifold_type

   PUBLIC :: cell_manifold_init, &
             cell_manifold_geometry, &
             cell_manifold_jacobian, &
             cell_manifold_residual, &
             cell_manifold_active_jacobian, &
             cell_manifold_chart_eval, &
             cell_manifold_chart_inv, &
             cell_manifold_chart_pullback_grad, &
             cell_manifold_proj_pressure, &
             cell_manifold_pack_hmat, &
             cell_manifold_unpack_hmat, &
             cell_manifold_status_label

   EXTERNAL :: dgeqrf, dgetrf, dgetrs, dorgqr, dpotrf, dpotrs

CONTAINS

! **************************************************************************************************
!> \brief Record constraint targets at q0, check row independence, and build [Y, Z].
!> \param m initialized constraint state
!> \param q0 initial canonical cell components
!> \param mask bit union of comp_* and gfix_* flags
!> \param info info_ok, info_degenerate, info_licq, or info_bad_input
!> \param lin_w optional 6 x n_lin linear coupling weights, one constraint per column
!> \param eq_pair optional 2 x n_eq index pairs [m, n] in 1..6 for s_m(q) = s_n(q)
! **************************************************************************************************
   SUBROUTINE cell_manifold_init(m, q0, mask, info, lin_w, eq_pair)
      TYPE(cell_manifold_type), INTENT(OUT)              :: m
      REAL(KIND=dp), INTENT(IN)                          :: q0(n_slot)
      INTEGER, INTENT(IN)                                :: mask
      INTEGER, INTENT(OUT)                               :: info
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: lin_w(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: eq_pair(:, :)

      INTEGER, PARAMETER :: &
         comp_flags(n_slot) = [comp_ax, comp_bx, comp_by, comp_cx, comp_cy, comp_cz], &
         gfix_flags(n_geom) = [gfix_a, gfix_b, gfix_c, gfix_alpha, gfix_beta, gfix_gamma, gfix_v]

      INTEGER                                            :: i, n_eq, n_lin
      REAL(KIND=dp)                                      :: weight_norm

      m%mask = mask
      m%q0 = q0

      CALL cell_manifold_geometry(q0, m%s0, info)
      IF (info /= info_ok) RETURN

      IF (IAND(mask, NOT(mask_all + gfix_all)) /= 0) THEN
         info = info_bad_input
         RETURN
      END IF

      DO i = 1, n_slot
         IF (IAND(mask, comp_flags(i)) /= 0) THEN
            m%n_slot_rows = m%n_slot_rows + 1
            m%slot_rows(m%n_slot_rows) = i
         END IF
      END DO
      DO i = 1, n_geom
         IF (IAND(mask, gfix_flags(i)) /= 0) THEN
            m%n_geom_rows = m%n_geom_rows + 1
            m%geom_rows(m%n_geom_rows) = i
         END IF
      END DO

      IF (PRESENT(lin_w)) THEN
         n_lin = SIZE(lin_w, 2)
         IF (SIZE(lin_w, 1) /= n_slot .OR. n_lin > n_slot) THEN
            info = info_bad_input
            RETURN
         END IF
         DO i = 1, n_lin
            ! Normalize each weight vector so residual checks do not depend on caller scaling.
            weight_norm = NORM2(lin_w(:, i))
            IF (weight_norm <= TINY(1.0_dp)) THEN
               info = info_bad_input
               RETURN
            END IF
            m%lin_w(:, i) = lin_w(:, i)/weight_norm
         END DO
         m%n_lin_rows = n_lin
      END IF

      IF (PRESENT(eq_pair)) THEN
         n_eq = SIZE(eq_pair, 2)
         IF (SIZE(eq_pair, 1) /= 2 .OR. n_eq > n_geom - 1) THEN
            info = info_bad_input
            RETURN
         END IF
         DO i = 1, n_eq
            IF (ANY(eq_pair(:, i) < geom_a) .OR. ANY(eq_pair(:, i) > geom_gamma) .OR. &
                eq_pair(1, i) == eq_pair(2, i)) THEN
               info = info_bad_input
               RETURN
            END IF
            m%eq_pair(:, i) = eq_pair(:, i)
         END DO
         m%n_eq_rows = n_eq
      END IF

      m%k = m%n_slot_rows + m%n_geom_rows + m%n_lin_rows + m%n_eq_rows
      IF (m%k > n_slot) THEN
         info = info_licq
         RETURN
      END IF
      m%ndof = n_slot - m%k

      CALL build_chart_basis(m, info)

   END SUBROUTINE cell_manifold_init

! **************************************************************************************************
!> \brief Compute cell lengths, angles in radians, and volume s(q).
!> \param q canonical cell components (AX, BX, BY, CX, CY, CZ)
!> \param s geometric quantities (a, b, c, alpha, beta, gamma, V)
!> \param info info_ok or info_degenerate
! **************************************************************************************************
   SUBROUTINE cell_manifold_geometry(q, s, info)
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: s(n_geom)
      INTEGER, INTENT(OUT)                               :: info

      REAL(KIND=dp)                                      :: ax, bx, by, cos_ang(3), cx, cy, cz, &
                                                            vec_a(3), vec_b(3), vec_c(3)

      s = 0.0_dp
      IF (.NOT. ALL(ieee_is_finite(q))) THEN
         info = info_degenerate
         RETURN
      END IF
      IF (ANY(ABS(q) > max_cell_coord)) THEN
         info = info_degenerate
         RETURN
      END IF

      ax = q(slot_ax)
      bx = q(slot_bx)
      by = q(slot_by)
      cx = q(slot_cx)
      cy = q(slot_cy)
      cz = q(slot_cz)

      ! Positive diagonal entries ensure V = AX*BY*CZ > 0 and right-handed cell vectors.
      IF (MIN(ax, by, cz) < min_cell_length) THEN
         info = info_degenerate
         RETURN
      END IF

      vec_a = [ax, 0.0_dp, 0.0_dp]
      vec_b = [bx, by, 0.0_dp]
      vec_c = [cx, cy, cz]
      s(geom_a) = NORM2(vec_a)
      s(geom_b) = NORM2(vec_b)
      s(geom_c) = NORM2(vec_c)
      cos_ang(1) = DOT_PRODUCT(vec_b, vec_c)/(s(geom_b)*s(geom_c))
      cos_ang(2) = DOT_PRODUCT(vec_a, vec_c)/(s(geom_a)*s(geom_c))
      cos_ang(3) = DOT_PRODUCT(vec_a, vec_b)/(s(geom_a)*s(geom_b))
      IF (ANY(ABS(cos_ang) > 1.0_dp - collinear_cos_margin)) THEN
         info = info_degenerate
         s = 0.0_dp
         RETURN
      END IF
      s(geom_alpha:geom_gamma) = ACOS(MIN(1.0_dp, MAX(-1.0_dp, cos_ang)))
      s(geom_v) = ax*by*cz
      info = info_ok

   END SUBROUTINE cell_manifold_geometry

! **************************************************************************************************
!> \brief Compute the analytic 7 x 6 Jacobian J = ds/dq of the geometric quantities.
!> \param q canonical cell components (AX, BX, BY, CX, CY, CZ)
!> \param J Jacobian matrix ds/dq
!> \param info info_ok or info_degenerate
! **************************************************************************************************
   SUBROUTINE cell_manifold_jacobian(q, J, info)
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: J(n_geom, n_slot)
      INTEGER, INTENT(OUT)                               :: info

      REAL(KIND=dp) :: ax, bx, by, cos_ang(3), cx, cy, cz, db_dq(n_slot), dc_dq(n_slot), &
         ddot_dq(n_slot), len_a, len_b, len_c, s(n_geom), sin_ang(3)

      J = 0.0_dp
      CALL cell_manifold_geometry(q, s, info)
      IF (info /= info_ok) RETURN

      ax = q(slot_ax)
      bx = q(slot_bx)
      by = q(slot_by)
      cx = q(slot_cx)
      cy = q(slot_cy)
      cz = q(slot_cz)
      len_a = s(geom_a)
      len_b = s(geom_b)
      len_c = s(geom_c)

      J(geom_a, slot_ax) = SIGN(1.0_dp, ax)

      db_dq = [0.0_dp, bx/len_b, by/len_b, 0.0_dp, 0.0_dp, 0.0_dp]
      J(geom_b, :) = db_dq

      dc_dq = [0.0_dp, 0.0_dp, 0.0_dp, cx/len_c, cy/len_c, cz/len_c]
      J(geom_c, :) = dc_dq

      cos_ang = [DOT_PRODUCT([bx, by, 0.0_dp], [cx, cy, cz])/(len_b*len_c), &
                 (ax*cx)/(len_a*len_c), &
                 (ax*bx)/(len_a*len_b)]
      sin_ang = SIN(s(geom_alpha:geom_gamma))

      ! alpha = angle(b, c), with b.c = BX*CX + BY*CY.
      ddot_dq = [0.0_dp, cx, cy, bx, by, 0.0_dp]
      J(geom_alpha, :) = angle_grad_row(cos_ang(1), sin_ang(1), ddot_dq, len_b, len_c, db_dq, dc_dq)

      ! beta = angle(a, c), with a.c = AX*CX.
      ddot_dq = [cx, 0.0_dp, 0.0_dp, ax, 0.0_dp, 0.0_dp]
      J(geom_beta, :) = angle_grad_row(cos_ang(2), sin_ang(2), ddot_dq, len_a, len_c, &
                                       J(geom_a, :), dc_dq)

      ! gamma = angle(a, b), with a.b = AX*BX.
      ddot_dq = [bx, ax, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      J(geom_gamma, :) = angle_grad_row(cos_ang(3), sin_ang(3), ddot_dq, len_a, len_b, &
                                        J(geom_a, :), db_dq)

      J(geom_v, :) = [by*cz, 0.0_dp, ax*cz, 0.0_dp, 0.0_dp, ax*by]

   END SUBROUTINE cell_manifold_jacobian

! **************************************************************************************************
!> \brief Evaluate the k active constraint residuals r(q).
!> \param m constraint state
!> \param q canonical cell vector
!> \param r residual vector of length m%k (zero when q satisfies all constraints)
!> \param info info_ok or info_degenerate
! **************************************************************************************************
   SUBROUTINE cell_manifold_residual(m, q, r, info)
      TYPE(cell_manifold_type), INTENT(IN)               :: m
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: r(m%k)
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: i, row
      REAL(KIND=dp)                                      :: s(n_geom)

      r = 0.0_dp
      info = info_ok
      IF (.NOT. ALL(ieee_is_finite(q))) THEN
         info = info_degenerate
         RETURN
      END IF
      IF (m%k == 0) RETURN

      ! Check cell validity even when all active rows are linear so collapsed
      ! or inverted cells cannot pass unnoticed.
      CALL cell_manifold_geometry(q, s, info)
      IF (info /= info_ok) RETURN

      DO i = 1, m%n_slot_rows
         r(i) = q(m%slot_rows(i)) - m%q0(m%slot_rows(i))
      END DO
      DO i = 1, m%n_geom_rows
         row = m%n_slot_rows + i
         r(row) = s(m%geom_rows(i)) - m%s0(m%geom_rows(i))
      END DO
      DO i = 1, m%n_lin_rows
         row = m%n_slot_rows + m%n_geom_rows + i
         r(row) = DOT_PRODUCT(m%lin_w(:, i), q) - DOT_PRODUCT(m%lin_w(:, i), m%q0)
      END DO
      DO i = 1, m%n_eq_rows
         row = m%n_slot_rows + m%n_geom_rows + m%n_lin_rows + i
         r(row) = s(m%eq_pair(1, i)) - s(m%eq_pair(2, i))
      END DO

   END SUBROUTINE cell_manifold_residual

! **************************************************************************************************
!> \brief Assemble the k x 6 Jacobian A(q) = dr/dq of the active constraints.
!> \param m constraint state
!> \param q canonical cell vector
!> \param A active constraint Jacobian in the same row order as cell_manifold_residual
!> \param info info_ok or info_degenerate
! **************************************************************************************************
   SUBROUTINE cell_manifold_active_jacobian(m, q, A, info)
      TYPE(cell_manifold_type), INTENT(IN)               :: m
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: A(m%k, n_slot)
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: i, row
      REAL(KIND=dp)                                      :: J(n_geom, n_slot)

      A = 0.0_dp
      info = info_ok
      IF (m%k == 0) RETURN

      CALL cell_manifold_jacobian(q, J, info)
      IF (info /= info_ok) RETURN

      DO i = 1, m%n_slot_rows
         A(i, m%slot_rows(i)) = 1.0_dp
      END DO
      DO i = 1, m%n_geom_rows
         row = m%n_slot_rows + i
         A(row, :) = J(m%geom_rows(i), :)
      END DO
      DO i = 1, m%n_lin_rows
         row = m%n_slot_rows + m%n_geom_rows + i
         A(row, :) = m%lin_w(:, i)
      END DO
      DO i = 1, m%n_eq_rows
         row = m%n_slot_rows + m%n_geom_rows + m%n_lin_rows + i
         A(row, :) = J(m%eq_pair(1, i), :) - J(m%eq_pair(2, i), :)
      END DO

   END SUBROUTINE cell_manifold_active_jacobian

! **************************************************************************************************
!> \brief Map reduced coordinates u to a canonical cell q(u) that satisfies r(q) = 0.
!> \param m constraint state (only last_iter and last_resid are updated)
!> \param u reduced coordinates of length m%ndof
!> \param q canonical cell vector on the constraint surface
!> \param info info_ok, info_degenerate, info_licq, or info_no_converge
!> \note Starting from qt = q0 + Z (u - u0), Newton's method solves r(qt + Y w) = 0
!>       for the k-component correction w along Y. Linear constraints already
!>       satisfy r(qt) = 0 and exit in zero iterations. When k = 0, q is a direct
!>       copy of u with zero roundoff.
! **************************************************************************************************
   SUBROUTINE cell_manifold_chart_eval(m, u, q, info)
      TYPE(cell_manifold_type), INTENT(INOUT)            :: m
      REAL(KIND=dp), INTENT(IN)                          :: u(m%ndof)
      REAL(KIND=dp), INTENT(OUT)                         :: q(n_slot)
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: i, it, lapack_info
      REAL(KIND=dp)                                      :: qt(n_slot), resid_tol, rn, &
                                                            s_chk(n_geom), row_norms(m%k), r(m%k), &
                                                            dw(m%k), B(m%k, m%k), A(m%k, n_slot)
      INTEGER                                            :: ipiv(m%k)

      ! Separate checks prevent 0*Inf and NaN comparisons from trapping under -ffpe-trap.
      IF (.NOT. ALL(ieee_is_finite(u))) THEN
         q = m%q0
         m%last_iter = 0
         m%last_resid = max_cell_coord
         info = info_no_converge
         RETURN
      END IF
      IF (ANY(ABS(u) > max_cell_coord)) THEN
         q = m%q0
         m%last_iter = 0
         m%last_resid = max_cell_coord
         info = info_no_converge
         RETURN
      END IF

      m%last_iter = 0
      m%last_resid = 0.0_dp
      info = info_ok

      IF (m%ndof == 0) THEN
         q = m%q0
         RETURN
      END IF

      IF (m%k == 0) THEN
         q = u
         CALL cell_manifold_geometry(q, s_chk, info)
         RETURN
      END IF

      qt = m%q0 + MATMUL(m%Z(:, 1:m%ndof), u - m%u0(1:m%ndof))

      ! Scale the convergence threshold by the largest active geometric target so
      ! volume and length residuals match double-precision accuracy on large cells.
      resid_tol = m%tol
      IF (m%n_geom_rows > 0) THEN
         resid_tol = m%tol*MAX(1.0_dp, MAXVAL(ABS(m%s0(m%geom_rows(1:m%n_geom_rows)))))
      END IF
      DO i = 1, m%n_eq_rows
         resid_tol = MAX(resid_tol, m%tol*ABS(m%s0(m%eq_pair(1, i))), &
                         m%tol*ABS(m%s0(m%eq_pair(2, i))))
      END DO

      CALL cell_manifold_residual(m, qt, r, info)
      IF (info /= info_ok) THEN
         q = qt
         m%last_resid = max_cell_coord
         RETURN
      END IF
      rn = NORM2(r)
      m%last_resid = rn

      DO it = 1, m%max_iter
         IF (rn <= resid_tol) EXIT
         CALL normalized_active_jacobian(m, qt, A, row_norms, info)
         IF (info /= info_ok) EXIT
         CALL factor_normal_block(m%k, A, m%Y(:, 1:m%k), B, ipiv, info)
         IF (info /= info_ok) EXIT
         dw = r/row_norms
         CALL dgetrs('N', m%k, 1, B, m%k, ipiv, dw, m%k, lapack_info)
         IF (lapack_info /= 0) THEN
            info = info_licq
            EXIT
         END IF
         qt = qt - MATMUL(m%Y(:, 1:m%k), dw)
         m%last_iter = it
         CALL cell_manifold_residual(m, qt, r, info)
         IF (info /= info_ok) THEN
            m%last_resid = max_cell_coord
            EXIT
         END IF
         rn = NORM2(r)
         m%last_resid = rn
         IF (.NOT. ieee_is_finite(rn)) THEN
            info = info_no_converge
            EXIT
         END IF
      END DO
      q = qt
      IF (info == info_ok .AND. rn > resid_tol) info = info_no_converge

   END SUBROUTINE cell_manifold_chart_eval

! **************************************************************************************************
!> \brief Map a canonical cell q back to reduced coordinates u = u0 + Z^T (q - q0).
!> \param m constraint state
!> \param q canonical cell vector
!> \param u reduced coordinates of length m%ndof
!> \note Because Z^T Z = I and Z^T Y = 0, this affine projection is the exact
!>       inverse of cell_manifold_chart_eval. When k = 0, u is a direct copy of q.
! **************************************************************************************************
   SUBROUTINE cell_manifold_chart_inv(m, q, u)
      TYPE(cell_manifold_type), INTENT(IN)               :: m
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: u(m%ndof)

      IF (m%ndof == 0) RETURN
      IF (m%k == 0) THEN
         u = q
         RETURN
      END IF
      u = m%u0(1:m%ndof) + MATMUL(TRANSPOSE(m%Z(:, 1:m%ndof)), q - m%q0)

   END SUBROUTINE cell_manifold_chart_inv

! **************************************************************************************************
!> \brief Transform a 6-component cell gradient g_q into the reduced gradient g_u = (dq/du)^T g_q.
!> \param m constraint state
!> \param q current canonical cell vector
!> \param g_q enthalpy gradient +dH/dq in canonical cell components
!> \param g_u enthalpy gradient +dH/du in reduced coordinates
!> \param info info_ok, info_degenerate, or info_licq
!> \note Differentiating r(q0 + Z(u - u0) + Y w(u)) = 0 gives dq/du = Z - Y (A Y)^{-1} A Z.
!>       Therefore g_u = Z^T (g_q - A^T lam), where (A Y)^T lam = Y^T g_q is the
!>       transpose of the k x k matrix B = A Y factored in cell_manifold_chart_eval.
! **************************************************************************************************
   SUBROUTINE cell_manifold_chart_pullback_grad(m, q, g_q, g_u, info)
      TYPE(cell_manifold_type), INTENT(IN)               :: m
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot), g_q(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: g_u(m%ndof)
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: lapack_info
      REAL(KIND=dp)                                      :: row_norms(m%k), lam(m%k), B(m%k, m%k), &
                                                            A(m%k, n_slot)
      INTEGER                                            :: ipiv(m%k)

      g_u = 0.0_dp
      info = info_ok
      IF (m%ndof == 0) RETURN

      IF (m%k == 0) THEN
         g_u = g_q
         RETURN
      END IF

      CALL normalized_active_jacobian(m, q, A, row_norms, info)
      IF (info /= info_ok) RETURN

      CALL factor_normal_block(m%k, A, m%Y(:, 1:m%k), B, ipiv, info)
      IF (info /= info_ok) RETURN

      lam = MATMUL(TRANSPOSE(m%Y(:, 1:m%k)), g_q)
      CALL dgetrs('T', m%k, 1, B, m%k, ipiv, lam, m%k, lapack_info)
      IF (lapack_info /= 0) THEN
         info = info_licq
         RETURN
      END IF
      g_u = MATMUL(TRANSPOSE(m%Z(:, 1:m%ndof)), g_q - MATMUL(TRANSPOSE(A), lam))

   END SUBROUTINE cell_manifold_chart_pullback_grad

! **************************************************************************************************
!> \brief Compute the pressure residual P_constr - pres_ext along unconstrained cell directions.
!> \param m constraint state
!> \param q current canonical cell components
!> \param g_phys enthalpy gradient +d(E + pres_ext*V)/dq from get_dg_dh
!> \param pdiff residual pressure P_constr - pres_ext along the allowed volume direction
!> \param projected .TRUE. if free directions change volume, .FALSE. if volume is fixed
!> \param info info_ok, info_degenerate, or info_licq
!> \note On curved constraints the tangent space ker(A(q)) rotates away from q0, so
!>       this routine projects dV/dq onto ker(A(q)) at the current cell q:
!>       vt = dV/dq - A(q)^T (A(q) A(q)^T)^{-1} A(q) (dV/dq).
!>       At any constrained minimum g_phys lies in row(A(q)), so pdiff vanishes.
! **************************************************************************************************
   SUBROUTINE cell_manifold_proj_pressure(m, q, g_phys, pdiff, projected, info)
      TYPE(cell_manifold_type), INTENT(IN)               :: m
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot), g_phys(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: pdiff
      LOGICAL, INTENT(OUT)                               :: projected
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: lapack_info
      REAL(KIND=dp) :: ax, by, cz, denom, dv_dq(n_slot), dv_norm2, s_chk(n_geom), vt(n_slot), &
         row_norms(m%k), mu(m%k), chol(m%k, m%k), A(m%k, n_slot)

      pdiff = 0.0_dp
      projected = .FALSE.
      info = info_ok
      IF (m%ndof == 0) RETURN

      CALL cell_manifold_geometry(q, s_chk, info)
      IF (info /= info_ok) RETURN

      ax = q(slot_ax)
      by = q(slot_by)
      cz = q(slot_cz)
      dv_dq = [by*cz, 0.0_dp, ax*cz, 0.0_dp, 0.0_dp, ax*by]
      dv_norm2 = DOT_PRODUCT(dv_dq, dv_dq)
      vt = dv_dq

      IF (m%k > 0) THEN
         CALL normalized_active_jacobian(m, q, A, row_norms, info)
         IF (info /= info_ok) RETURN
         chol = MATMUL(A, TRANSPOSE(A))
         CALL dpotrf('U', m%k, chol, m%k, lapack_info)
         IF (lapack_info /= 0) THEN
            info = info_licq
            RETURN
         END IF
         mu = MATMUL(A, dv_dq)
         CALL dpotrs('U', m%k, 1, chol, m%k, mu, m%k, lapack_info)
         IF (lapack_info /= 0) THEN
            info = info_licq
            RETURN
         END IF
         vt = dv_dq - MATMUL(TRANSPOSE(A), mu)
      END IF

      denom = DOT_PRODUCT(vt, vt)
      IF (denom <= proj_pressure_rel_tol*dv_norm2) RETURN

      projected = .TRUE.
      ! For isotropic stress g_phys = -(P_int - pres_ext) * dV/dq, so the minus
      ! sign yields pdiff = P_constr - pres_ext.
      pdiff = -DOT_PRODUCT(g_phys, vt)/denom

   END SUBROUTINE cell_manifold_proj_pressure

! **************************************************************************************************
!> \brief Pack the upper-triangular canonical cell matrix into q(1:6).
!> \param hmat canonical 3x3 cell matrix (a along x, b in the xy plane)
!> \return q = (AX, BX, BY, CX, CY, CZ)
! **************************************************************************************************
   PURE FUNCTION cell_manifold_pack_hmat(hmat) RESULT(q)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat
      REAL(KIND=dp), DIMENSION(n_slot)                   :: q

      q = [hmat(1, 1), hmat(1, 2), hmat(2, 2), hmat(1, 3), hmat(2, 3), hmat(3, 3)]

   END FUNCTION cell_manifold_pack_hmat

! **************************************************************************************************
!> \brief Write q(1:6) into the upper-triangular canonical cell matrix.
!> \param q canonical cell vector (AX, BX, BY, CX, CY, CZ)
!> \param hmat canonical 3x3 cell matrix with zeroed lower triangle
! **************************************************************************************************
   PURE SUBROUTINE cell_manifold_unpack_hmat(q, hmat)
      REAL(KIND=dp), DIMENSION(n_slot), INTENT(IN)       :: q
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: hmat

      hmat = 0.0_dp
      hmat(1, 1) = q(slot_ax)
      hmat(1, 2) = q(slot_bx)
      hmat(2, 2) = q(slot_by)
      hmat(1, 3) = q(slot_cx)
      hmat(2, 3) = q(slot_cy)
      hmat(3, 3) = q(slot_cz)

   END SUBROUTINE cell_manifold_unpack_hmat

! **************************************************************************************************
!> \brief Return a short description of a cell_manifold_utils status code.
!> \param info status code from a cell_manifold_* routine
!> \return diagnostic text
! **************************************************************************************************
   PURE FUNCTION cell_manifold_status_label(info) RESULT(label)
      INTEGER, INTENT(IN)                                :: info
      CHARACTER(LEN=default_string_length)               :: label

      SELECT CASE (info)
      CASE (info_ok)
         label = "ok"
      CASE (info_degenerate)
         label = "degenerate or non-positive cell"
      CASE (info_licq)
         label = "cell constraints are linearly dependent or ill-conditioned"
      CASE (info_no_converge)
         label = "nonlinear cell constraint solver did not converge"
      CASE (info_bad_input)
         label = "invalid constraint mask or row specification"
      CASE DEFAULT
         label = "unknown cell-constraint status"
      END SELECT

   END FUNCTION cell_manifold_status_label

! **************************************************************************************************
!> \brief Build the fixed orthonormal split [Y, Z] and reduced anchor u0 at m%q0.
!> \param m constraint state
!> \param info info_ok, info_degenerate, or info_licq
!> \note Cholesky factorization of the row-normalized Gram matrix A A^T rejects both
!>       exact row dependence and near-collinear rows before QR builds [Y, Z].
!>       Orienting each QR column so its largest-magnitude entry is positive makes
!>       the basis independent of LAPACK Householder sign choices.
! **************************************************************************************************
   SUBROUTINE build_chart_basis(m, info)
      TYPE(cell_manifold_type), INTENT(INOUT)            :: m
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: i, info_lapack, j, lwork, pmax
      REAL(KIND=dp)                                      :: Qfull(n_slot, n_slot), work_query(1)
      REAL(KIND=dp), ALLOCATABLE                         :: work(:)
      REAL(KIND=dp)                                      :: tau(m%k), row_norms(m%k), &
                                                            gram(m%k, m%k), gdiag(m%k), &
                                                            AT(n_slot, m%k), A(m%k, n_slot)

      IF (m%k == 0) THEN
         m%ndof = n_slot
         m%Y = 0.0_dp
         m%Z = 0.0_dp
         DO i = 1, n_slot
            m%Z(i, i) = 1.0_dp
         END DO
         m%u0 = m%q0
         info = info_ok
         RETURN
      END IF

      CALL normalized_active_jacobian(m, m%q0, A, row_norms, info)
      IF (info /= info_ok) RETURN

      gram = MATMUL(A, TRANSPOSE(A))
      CALL dpotrf('U', m%k, gram, m%k, info_lapack)
      IF (info_lapack /= 0) THEN
         info = info_licq
         RETURN
      END IF
      gdiag = [(gram(i, i), i=1, m%k)]
      IF (MINVAL(gdiag) <= licq_pivot_tol*MAXVAL(gdiag)) THEN
         info = info_licq
         RETURN
      END IF

      AT = TRANSPOSE(A)
      CALL dgeqrf(n_slot, m%k, AT, n_slot, tau, work_query, -1, info_lapack)
      lwork = INT(work_query(1))
      ALLOCATE (work(MAX(1, lwork)))

      CALL dgeqrf(n_slot, m%k, AT, n_slot, tau, work, SIZE(work), info_lapack)
      IF (info_lapack /= 0) THEN
         info = info_licq
         DEALLOCATE (work)
         RETURN
      END IF
      Qfull = 0.0_dp
      Qfull(1:n_slot, 1:m%k) = AT(1:n_slot, 1:m%k)
      CALL dorgqr(n_slot, n_slot, m%k, Qfull, n_slot, tau, work, SIZE(work), info_lapack)
      DEALLOCATE (work)
      IF (info_lapack /= 0) THEN
         info = info_licq
         RETURN
      END IF

      DO j = 1, n_slot
         pmax = MAXLOC(ABS(Qfull(:, j)), DIM=1)
         IF (Qfull(pmax, j) < 0.0_dp) Qfull(:, j) = -Qfull(:, j)
      END DO

      m%ndof = n_slot - m%k
      m%Y = 0.0_dp
      m%Y(:, 1:m%k) = Qfull(:, 1:m%k)
      m%Z = 0.0_dp
      m%Z(:, 1:m%ndof) = Qfull(:, m%k + 1:n_slot)
      m%u0 = 0.0_dp
      m%u0(1:m%ndof) = MATMUL(TRANSPOSE(m%Z(:, 1:m%ndof)), m%q0)

   END SUBROUTINE build_chart_basis

! **************************************************************************************************
!> \brief Assemble A(q) and scale each row to unit 2-norm.
!> \param m constraint state (with m%k > 0)
!> \param q canonical cell vector
!> \param A row-normalized active Jacobian
!> \param row_norms original Euclidean norm of each row before normalization
!> \param info info_ok, info_degenerate, or info_licq
!> \note Normalizing rows keeps linear solves independent of unit/scale differences
!>       when large volume-gradient rows mix with unit component rows.
! **************************************************************************************************
   SUBROUTINE normalized_active_jacobian(m, q, A, row_norms, info)
      TYPE(cell_manifold_type), INTENT(IN)               :: m
      REAL(KIND=dp), INTENT(IN)                          :: q(n_slot)
      REAL(KIND=dp), INTENT(OUT)                         :: A(m%k, n_slot), row_norms(m%k)
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: i

      row_norms = 0.0_dp
      CALL cell_manifold_active_jacobian(m, q, A, info)
      IF (info /= info_ok) RETURN

      DO i = 1, m%k
         row_norms(i) = NORM2(A(i, :))
         IF (row_norms(i) <= TINY(1.0_dp)) THEN
            info = info_licq
            RETURN
         END IF
         A(i, :) = A(i, :)/row_norms(i)
      END DO

   END SUBROUTINE normalized_active_jacobian

! **************************************************************************************************
!> \brief Form and LU-factor the k x k block B = A(q) * Y.
!> \param k number of active constraints (k > 0)
!> \param A row-normalized k x 6 active Jacobian at q
!> \param Y fixed 6 x k constrained-direction basis from q0
!> \param B LU factors of A * Y written in-place by DGETRF
!> \param ipiv pivot indices from DGETRF
!> \param info info_ok or info_licq
!> \note After DGETRF, the diagonal of B holds the pivots of U. A tiny pivot ratio
!>       means ker(A(q)) now intersects span(Y), so q(u) is no longer single-valued
!>       along Y.
! **************************************************************************************************
   SUBROUTINE factor_normal_block(k, A, Y, B, ipiv, info)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: A(k, n_slot), Y(n_slot, k)
      REAL(KIND=dp), INTENT(OUT)                         :: B(k, k)
      INTEGER, INTENT(OUT)                               :: ipiv(k), info

      INTEGER                                            :: i, lapack_info
      REAL(KIND=dp)                                      :: u_pivots(k)

      B = MATMUL(A, Y)
      CALL dgetrf(k, k, B, k, ipiv, lapack_info)
      u_pivots = [(ABS(B(i, i)), i=1, k)]
      IF (lapack_info /= 0 .OR. MINVAL(u_pivots) <= licq_pivot_tol*MAXVAL(u_pivots)) THEN
         info = info_licq
      ELSE
         info = info_ok
      END IF

   END SUBROUTINE factor_normal_block

! **************************************************************************************************
!> \brief Differentiate angle(u, v) = acos((u.v) / (|u| |v|)) with respect to q.
!> \param cos_ang cos(angle)
!> \param sin_ang sin(angle)
!> \param ddot_dq derivative d(u.v)/dq of the dot product
!> \param len_u length |u|
!> \param len_v length |v|
!> \param dlen_u_dq derivative d|u|/dq
!> \param dlen_v_dq derivative d|v|/dq
!> \return row vector d(angle)/dq
! **************************************************************************************************
   PURE FUNCTION angle_grad_row(cos_ang, sin_ang, ddot_dq, len_u, len_v, dlen_u_dq, dlen_v_dq) &
      RESULT(d_ang_dq)
      REAL(KIND=dp), INTENT(IN)                          :: cos_ang, sin_ang, ddot_dq(n_slot), &
                                                            len_u, len_v, dlen_u_dq(n_slot), &
                                                            dlen_v_dq(n_slot)
      REAL(KIND=dp)                                      :: d_ang_dq(n_slot)

      REAL(KIND=dp)                                      :: dcos_dq(n_slot)

      dcos_dq = ddot_dq/(len_u*len_v) - cos_ang*(dlen_u_dq/len_u + dlen_v_dq/len_v)
      d_ang_dq = -dcos_dq/sin_ang

   END FUNCTION angle_grad_row

END MODULE cell_manifold_utils
