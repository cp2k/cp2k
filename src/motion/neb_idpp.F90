!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Initial paths from an image-dependent pair potential (IDPP).
! **************************************************************************************************
MODULE neb_idpp
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: bohr
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: neb_initialize_idpp

   REAL(KIND=dp), PARAMETER :: soft_radius = 0.01_dp ! Angstrom

CONTAINS

! **************************************************************************************************
!> \brief Relax generated images on the IDPP surface with weak elastic coupling.
!>        Pair distances are interpolated separately between consecutive anchors.
!> \param coords Cartesian coordinates in atomic units, indexed by component and image
!> \param anchor true for user-supplied images, which are left unchanged
!> \param cell fixed simulation cell
!> \param iw initialization output unit
! **************************************************************************************************
   SUBROUTINE neb_initialize_idpp(coords, anchor, cell, iw)
      REAL(KIND=dp), INTENT(INOUT)                       :: coords(:, :)
      LOGICAL, INTENT(IN)                                :: anchor(:)
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iw

      INTEGER                                            :: first, last
      REAL(KIND=dp), ALLOCATABLE                         :: path(:, :)

      IF (SIZE(coords, 1) <= 3) RETURN
      first = 1
      DO last = 2, SIZE(anchor)
         IF (.NOT. anchor(last)) CYCLE
         IF (last - first > 1) THEN
            ALLOCATE (path(SIZE(coords, 1), last - first + 1))
            path(:, :) = coords(:, first:last)/bohr
            CALL relax_segment(path, cell, iw)
            coords(:, first + 1:last - 1) = path(:, 2:last - first)*bohr
            DEALLOCATE (path)
         END IF
         first = last
      END DO
   END SUBROUTINE neb_initialize_idpp

! **************************************************************************************************
!> \brief Minimize the sum of regularized IDPP image energies and harmonic spring energies.
!> \param path path in Angstrom, including its fixed endpoints
!> \param cell fixed simulation cell
!> \param iw initialization output unit
! **************************************************************************************************
   SUBROUTINE relax_segment(path, cell, iw)
      REAL(KIND=dp), INTENT(INOUT)                       :: path(:, :)
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iw

      INTEGER, PARAMETER                                 :: max_iterations = 250, max_backtrack = 25
      REAL(KIND=dp), PARAMETER                            :: spring = 0.02_dp, tolerance = 1.0E-4_dp
      INTEGER                                            :: a, backtrack, b, dim_min, image, it, k, ni, ns
      REAL(KIND=dp)                                      :: alpha, e, e_trial, min_distance, min_initial, norm, slope
      REAL(KIND=dp)                                      :: delta(3), direction(3)
      REAL(KIND=dp), ALLOCATABLE                         :: anchor_distances(:, :), gradient(:, :), trial(:, :), trial_g(:, :)

      ns = SIZE(path, 1)
      ni = SIZE(path, 2)
      ALLOCATE (gradient(ns, ni), trial(ns, ni), trial_g(ns, ni), anchor_distances(ns/3*(ns/3 - 1)/2, 2))
      k = 0
      DO a = 1, ns - 3, 3
         DO b = a + 3, ns, 3
            k = k + 1
            anchor_distances(k, 1) = NORM2(pair_delta(path(a:a + 2, 1) - path(b:b + 2, 1)))
            anchor_distances(k, 2) = NORM2(pair_delta(path(a:a + 2, ni) - path(b:b + 2, ni)))
         END DO
      END DO
      CALL evaluate(path, e, gradient, min_initial)
      ! A perfectly overlapping pair has no defined separation direction. Break that
      ! symmetry deterministically, transverse to its orientation in the first anchor.
      DO image = 2, ni - 1
         DO a = 1, ns - 3, 3
            DO b = a + 3, ns, 3
               delta = pair_delta(path(a:a + 2, image) - path(b:b + 2, image))
               IF (NORM2(delta) >= soft_radius) CYCLE
               delta = pair_delta(path(a:a + 2, 1) - path(b:b + 2, 1))
               dim_min = MINLOC(ABS(delta), 1)
               direction = 0.0_dp
               direction(dim_min) = 1.0_dp
               IF (NORM2(delta) > soft_radius) THEN
                  delta = delta/NORM2(delta)
                  direction = direction - DOT_PRODUCT(direction, delta)*delta
               END IF
               direction = 3.0_dp*soft_radius*direction/NORM2(direction)
               path(a:a + 2, image) = path(a:a + 2, image) + direction
               path(b:b + 2, image) = path(b:b + 2, image) - direction
            END DO
         END DO
      END DO
      CALL evaluate(path, e, gradient, min_distance)
      alpha = 1.0_dp
      DO it = 1, max_iterations
         norm = MAXVAL(ABS(gradient))
         IF (norm < tolerance) EXIT
         alpha = MIN(alpha, 0.1_dp/norm)
         slope = SUM(gradient**2)
         DO backtrack = 1, max_backtrack
            trial(:, :) = path - alpha*gradient
            CALL evaluate(trial, e_trial, trial_g, min_distance)
            IF (e_trial <= e - 1.0E-4_dp*alpha*slope) EXIT
            alpha = 0.5_dp*alpha
         END DO
         IF (backtrack > max_backtrack) EXIT
         path(:, :) = trial
         gradient(:, :) = trial_g
         e = e_trial
         alpha = MIN(1.5_dp*alpha, 20.0_dp)
      END DO
      IF (iw > 0) THEN
         CALL evaluate(path, e, gradient, min_distance)
         WRITE (iw, '(T2,A,F14.8,A,F14.8)') &
            "NEB| IDPP minimum pair distance [angstrom]: ", min_initial, " -> ", min_distance
         IF (MAXVAL(ABS(gradient)) >= tolerance) THEN
            WRITE (iw, '(T2,A)') "NEB| IDPP reached its stopping limit. Using the last relaxed images."
         END IF
      END IF

   CONTAINS

! **************************************************************************************************
!> \brief Evaluate the auxiliary elastic-band objective and gradient.
!> \param coords path at which to evaluate the objective, in Angstrom
!> \param energy value of the objective
!> \param grad Cartesian gradient, zero at both anchors
!> \param shortest smallest unregularized pair distance in any generated image
! **************************************************************************************************
      SUBROUTINE evaluate(coords, energy, grad, shortest)
      REAL(KIND=dp), INTENT(IN)                          :: coords(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: energy, grad(:, :), shortest

      INTEGER                                            :: i, ia, ib, pair
      REAL(KIND=dp)                                      :: d, dg(3), dx(3), fraction, r, r2, weight

         energy = 0.5_dp*spring*SUM((coords(:, 2:ni) - coords(:, 1:ni - 1))**2)
         grad = 0.0_dp
         shortest = HUGE(1.0_dp)
         DO i = 2, ni - 1
            grad(:, i) = spring*(2.0_dp*coords(:, i) - coords(:, i - 1) - coords(:, i + 1))
            fraction = REAL(i - 1, dp)/REAL(ni - 1, dp)
            pair = 0
            DO ia = 1, ns - 3, 3
               DO ib = ia + 3, ns, 3
                  pair = pair + 1
                  dx = pair_delta(coords(ia:ia + 2, i) - coords(ib:ib + 2, i))
                  r2 = DOT_PRODUCT(dx, dx)
                  shortest = MIN(shortest, SQRT(r2))
                  r2 = r2 + soft_radius**2
                  r = SQRT(r2)
                  d = (1.0_dp - fraction)*anchor_distances(pair, 1) + fraction*anchor_distances(pair, 2)
                  d = SQRT(d*d + soft_radius**2)
                  weight = 1.0_dp/(r2*r2)
                  energy = energy + weight*(r - d)**2
                  dg = (2.0_dp*weight*(r - d)*(2.0_dp*d - r)/r2)*dx
                  grad(ia:ia + 2, i) = grad(ia:ia + 2, i) + dg
                  grad(ib:ib + 2, i) = grad(ib:ib + 2, i) - dg
               END DO
            END DO
         END DO
      END SUBROUTINE evaluate

! **************************************************************************************************
!> \brief Cartesian pair displacement in Angstrom, folded for periodic directions.
!> \param dx unreduced pair displacement
!> \return nearest-image pair displacement
! **************************************************************************************************
      FUNCTION pair_delta(dx) RESULT(dr)
      REAL(KIND=dp), INTENT(IN)                          :: dx(3)
      REAL(KIND=dp)                                      :: dr(3)

         dr = dx
         IF (ANY(cell%perd /= 0)) dr = pbc(dx*bohr, cell)/bohr
      END FUNCTION pair_delta
   END SUBROUTINE relax_segment
END MODULE neb_idpp
