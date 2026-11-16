!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Local UPF matrix elements from the original radial potential.
!>        The Gaussian core electrostatics supplies -Z erf(sqrt(alpha)*r)/r.
!>        Its complement is integrated with analytic angles and radial quadrature.
!>        Gaussian products, center derivatives and RI integrals share this operator.
! **************************************************************************************************
MODULE upf_local_integrals
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE mathlib,                         ONLY: binomial
   USE orbital_pointers,                ONLY: coset,&
                                              indco,&
                                              ncoset
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_overlap,&
                                              upf_projector_type
   USE upf_types,                       ONLY: atom_upfpot_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'upf_local_integrals'

   TYPE upf_local_type
      TYPE(upf_projector_type) :: radial
      REAL(KIND=dp) :: eps = 0.0_dp
   END TYPE upf_local_type

   PUBLIC :: upf_local_type, upf_local_init, upf_local_matrix, upf_local_ri

CONTAINS

! **************************************************************************************************
!> \brief Prepare the short-range local potential outside parallel regions.
!> \param pot Original UPF data, with PP_LOCAL already converted to Hartree.
!> \param alpha Gaussian core-charge exponent, not its screening radius.
!> \param degree Largest Gaussian product degree, including derivatives or moments.
!> \param local Immutable integration data. Reinitialization retains its threshold.
!> \param eps Optional threshold for estimating the potential's neighbor-list radius.
!>        This is NOT an integration tolerance: EPS_PPL is normally 1e-2.
!>        All original radial samples, including the tail, enter the quadrature.
! **************************************************************************************************
   SUBROUTINE upf_local_init(pot, alpha, degree, local, eps)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      INTEGER, INTENT(IN)                                :: degree
      TYPE(upf_local_type), INTENT(INOUT)                :: local
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps

      INTEGER                                            :: i, last, n
      REAL(KIND=dp), ALLOCATABLE                         :: short_range(:)
      TYPE(atom_upfpot_type)                             :: scalar

      CPASSERT(alpha > 0.0_dp)
      CPASSERT(SIZE(pot%r) == SIZE(pot%vlocal))
      CPASSERT(SIZE(pot%r) >= 3)
      IF (PRESENT(eps)) local%eps = eps
      CPASSERT(local%eps >= 0.0_dp)
      ALLOCATE (short_range(SIZE(pot%vlocal)))
      short_range(:) = pot%vlocal
      DO i = 1, SIZE(pot%r)
         IF (pot%r(i) > 0.0_dp) THEN
            short_range(i) = short_range(i) + pot%zion*ERF(SQRT(alpha)*pot%r(i))/pot%r(i)
         ELSE
            short_range(i) = short_range(i) + 2.0_dp*pot%zion*SQRT(alpha/pi)
         END IF
      END DO
      last = 0
      DO i = SIZE(short_range), 1, -1
         IF (ABS(short_range(i)) > local%eps) THEN
            last = i
            EXIT
         END IF
      END DO
      n = SIZE(short_range)
      scalar%number_of_proj = 1
      ALLOCATE (scalar%r(n), scalar%rab(n), scalar%lbeta(1))
      scalar%r(:) = pot%r(:n)
      scalar%rab(:) = pot%rab(:n)
      scalar%lbeta = 0
      ALLOCATE (scalar%beta(n, 1), scalar%dion(1, 1))
      ! Y00 = 1/sqrt(4*pi); absorb it to integrate a scalar radial potential.
      scalar%beta(:, 1) = SQRT(fourpi)*short_range(:n)
      scalar%dion = 0.0_dp
      CALL upf_projector_init(scalar, degree, local%radial)
      ! Screening estimates the support in a pair-list bound. It must never
      ! truncate the function actually integrated within a retained pair.
      local%radial%radius = 0.0_dp
      IF (last > 0) local%radial%radius = pot%r(last)
   END SUBROUTINE upf_local_init

! **************************************************************************************************
!> \brief Evaluate primitive three-center local integrals by the Gaussian product theorem.
!> \param local Prepared radial potential.
!> \param la Largest Cartesian degree at A; all degrees 0..la are returned.
!> \param alpha Primitive exponents at A.
!> \param lb Largest Cartesian degree at B; all degrees 0..lb are returned.
!> \param beta Primitive exponents at B.
!> \param rab B minus A.
!> \param rac Potential center C minus A.
!> \param v Matrix elements in primitive-major Cartesian order.
!> \param da Optional derivatives with respect to A, at fixed B and C.
!> \param db Optional derivatives with respect to B, at fixed A and C.
!>        The derivative with respect to C is -(da+db).
! **************************************************************************************************
   SUBROUTINE upf_local_matrix(local, la, alpha, lb, beta, rab, rac, v, da, db)
      TYPE(upf_local_type), INTENT(IN)                   :: local
      INTEGER, INTENT(IN)                                :: la
      REAL(KIND=dp), INTENT(IN)                          :: alpha(:)
      INTEGER, INTENT(IN)                                :: lb
      REAL(KIND=dp), INTENT(IN)                          :: beta(:), rab(3), rac(3)
      REAL(KIND=dp), INTENT(OUT)                         :: v(:, :)
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: da(:, :, :), db(:, :, :)

      INTEGER                                            :: a, aa, axis, b, bb, da_max, db_max, &
                                                            degree, i, ia, ib, ix, iy, iz, j, k, &
                                                            pa(3), pb(3), shifted(3)
      REAL(KIND=dp)                                      :: center(3), p, pref, shift_a(3), &
                                                            shift_b(3), value
      REAL(KIND=dp), ALLOCATABLE                         :: moments(:, :, :), &
                                                            polynomial(:, :, :, :), work(:, :)

      CPASSERT(la >= 0 .AND. lb >= 0)
      CPASSERT(ALL(alpha > 0.0_dp) .AND. ALL(beta > 0.0_dp))
      CPASSERT(SIZE(v, 1) == SIZE(alpha)*ncoset(la))
      CPASSERT(SIZE(v, 2) == SIZE(beta)*ncoset(lb))
      degree = la + lb
      IF (PRESENT(da) .OR. PRESENT(db)) degree = degree + 1
      da_max = la
      db_max = lb
      IF (PRESENT(da)) THEN
         CPASSERT(ALL(SHAPE(da) == [SIZE(v, 1), SIZE(v, 2), 3]))
         da = 0.0_dp
         da_max = la + 1
      END IF
      IF (PRESENT(db)) THEN
         CPASSERT(ALL(SHAPE(db) == [SIZE(v, 1), SIZE(v, 2), 3]))
         db = 0.0_dp
         db_max = lb + 1
      END IF
      ALLOCATE (moments(ncoset(degree), 1, 1), polynomial(0:da_max, 0:db_max, 0:degree, 3))
      ALLOCATE (work(ncoset(da_max), ncoset(db_max)))
      v = 0.0_dp
      DO ib = 1, SIZE(beta)
         DO ia = 1, SIZE(alpha)
            p = alpha(ia) + beta(ib)
            shift_a = beta(ib)*rab/p
            shift_b = -alpha(ia)*rab/p
            center = rac - shift_a
            pref = EXP(-alpha(ia)*beta(ib)/p*SUM(rab**2))
            IF (pref == 0.0_dp) CYCLE
            CALL upf_projector_overlap(local%radial, degree, [p], center, 0, moments)
            polynomial = 0.0_dp
            DO axis = 1, 3
               DO a = 0, da_max
                  DO b = 0, MIN(db_max, degree - a)
                     DO i = 0, a
                        DO j = 0, b
                           polynomial(a, b, i + j, axis) = polynomial(a, b, i + j, axis) &
                                                       + binomial(a, i)*binomial(b, j)*shift_a(axis)**(a - i)*shift_b(axis)**(b - j)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            DO b = 1, ncoset(db_max)
               pb = indco(:, b)
               DO a = 1, ncoset(da_max)
                  pa = indco(:, a)
                  IF (SUM(pa + pb) > degree) CYCLE
                  value = 0.0_dp
                  DO iz = 0, pa(3) + pb(3)
                     DO iy = 0, pa(2) + pb(2)
                        DO ix = 0, pa(1) + pb(1)
                           k = coset(ix, iy, iz)
                           value = value + polynomial(pa(1), pb(1), ix, 1)*polynomial(pa(2), pb(2), iy, 2) &
                                   *polynomial(pa(3), pb(3), iz, 3)*moments(k, 1, 1)
                        END DO
                     END DO
                  END DO
                  work(a, b) = pref*value
               END DO
            END DO
            DO b = 1, ncoset(lb)
               pb = indco(:, b)
               bb = (ib - 1)*ncoset(lb) + b
               DO a = 1, ncoset(la)
                  pa = indco(:, a)
                  aa = (ia - 1)*ncoset(la) + a
                  v(aa, bb) = work(a, b)
                  DO axis = 1, 3
                     IF (PRESENT(da)) THEN
                        shifted = pa
                        shifted(axis) = shifted(axis) + 1
                        da(aa, bb, axis) = 2*alpha(ia)*work(coset(shifted(1), shifted(2), shifted(3)), b)
                        IF (pa(axis) > 0) THEN
                           shifted = pa
                           shifted(axis) = shifted(axis) - 1
                           da(aa, bb, axis) = da(aa, bb, axis) &
                                              - pa(axis)*work(coset(shifted(1), shifted(2), shifted(3)), b)
                        END IF
                     END IF
                     IF (PRESENT(db)) THEN
                        shifted = pb
                        shifted(axis) = shifted(axis) + 1
                        db(aa, bb, axis) = 2*beta(ib)*work(a, coset(shifted(1), shifted(2), shifted(3)))
                        IF (pb(axis) > 0) THEN
                           shifted = pb
                           shifted(axis) = shifted(axis) - 1
                           db(aa, bb, axis) = db(aa, bb, axis) &
                                              - pb(axis)*work(a, coset(shifted(1), shifted(2), shifted(3)))
                        END IF
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE upf_local_matrix

! **************************************************************************************************
!> \brief Local potential against a single Cartesian Gaussian for RI.
!> \param local Prepared radial potential.
!> \param la Largest Cartesian degree.
!> \param zeta Primitive exponents.
!> \param rac Potential center minus Gaussian center.
!> \param v Primitive-major Cartesian integrals.
!> \param dv Optional Gaussian-center derivatives, in x,y,z order.
! **************************************************************************************************
   SUBROUTINE upf_local_ri(local, la, zeta, rac, v, dv)
      TYPE(upf_local_type), INTENT(IN)                   :: local
      INTEGER, INTENT(IN)                                :: la
      REAL(KIND=dp), INTENT(IN)                          :: zeta(:), rac(3)
      REAL(KIND=dp), INTENT(OUT)                         :: v(:)
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: dv(:, :)

      INTEGER                                            :: n, nder
      REAL(KIND=dp), ALLOCATABLE                         :: integrals(:, :, :)

      nder = 0
      IF (PRESENT(dv)) nder = 1
      n = SIZE(zeta)*ncoset(la)
      CPASSERT(SIZE(v) == n)
      ALLOCATE (integrals(n, 1, ncoset(nder)))
      CALL upf_projector_overlap(local%radial, la, zeta, rac, nder, integrals)
      v = integrals(:, 1, 1)
      IF (PRESENT(dv)) THEN
         CPASSERT(ALL(SHAPE(dv) == [n, 3]))
         dv = integrals(:, 1, 2:4)
      END IF
   END SUBROUTINE upf_local_ri
END MODULE upf_local_integrals
