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
!> \param sequential grow images from the anchors instead of using the linear guess
! **************************************************************************************************
   SUBROUTINE neb_initialize_idpp(coords, anchor, cell, iw, sequential)
      REAL(KIND=dp), INTENT(INOUT)                       :: coords(:, :)
      LOGICAL, INTENT(IN)                                :: anchor(:)
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iw
      LOGICAL, INTENT(IN)                                :: sequential

      INTEGER                                            :: first, last
      REAL(KIND=dp), ALLOCATABLE                         :: path(:, :)

      IF (SIZE(coords, 1) <= 3) RETURN
      first = 1
      DO last = 2, SIZE(anchor)
         IF (.NOT. anchor(last)) CYCLE
         IF (last - first > 1) THEN
            ALLOCATE (path(SIZE(coords, 1), last - first + 1))
            ! Both IDPP methods use angstrom internally, including the reported pair distances.
            path(:, :) = coords(:, first:last)/bohr
            CALL relax_segment(path, cell, iw, sequential)
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
!> \param sequential grow the segment using S-IDPP
! **************************************************************************************************
   SUBROUTINE relax_segment(path, cell, iw, sequential)
      REAL(KIND=dp), INTENT(INOUT)                       :: path(:, :)
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iw
      LOGICAL, INTENT(IN)                                :: sequential

      INTEGER, PARAMETER                                 :: max_iterations = 250, max_backtrack = 25
      REAL(KIND=dp), PARAMETER                            :: spring = 0.02_dp, tolerance = 1.0E-4_dp
      INTEGER                                            :: a, backtrack, b, image, it, k, ni, ns
      REAL(KIND=dp)                                      :: alpha, e, e_trial, min_distance, min_initial, norm, slope
      REAL(KIND=dp), ALLOCATABLE                         :: anchor_distances(:, :), gradient(:, :), trial(:, :), trial_g(:, :)

      ns = SIZE(path, 1)
      ni = SIZE(path, 2)
      ALLOCATE (anchor_distances(ns/3*(ns/3 - 1)/2, 2))
      k = 0
      DO a = 1, ns - 3, 3
         DO b = a + 3, ns, 3
            k = k + 1
            anchor_distances(k, 1) = NORM2(pair_delta(path(a:a + 2, 1) - path(b:b + 2, 1)))
            anchor_distances(k, 2) = NORM2(pair_delta(path(a:a + 2, ni) - path(b:b + 2, ni)))
         END DO
      END DO
      IF (sequential) THEN
         CALL grow_segment()
         RETURN
      END IF
      ALLOCATE (gradient(ns, ni), trial(ns, ni), trial_g(ns, ni))
      CALL evaluate(path, e, gradient, min_initial)
      DO image = 2, ni - 1
         CALL perturb_short_contacts(path(:, image))
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

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: fraction

         energy = 0.5_dp*spring*SUM((coords(:, 2:ni) - coords(:, 1:ni - 1))**2)
         grad = 0.0_dp
         shortest = HUGE(1.0_dp)
         DO i = 2, ni - 1
            grad(:, i) = spring*(2.0_dp*coords(:, i) - coords(:, i - 1) - coords(:, i + 1))
            fraction = REAL(i - 1, dp)/REAL(ni - 1, dp)
            CALL evaluate_image(coords(:, i), fraction, energy, grad(:, i), shortest)
         END DO
      END SUBROUTINE evaluate

! **************************************************************************************************
!> \brief Accumulate the regularized pair objective and its gradient for one image.
!> \param xyz image coordinates in Angstrom
!> \param fraction fractional position of the target distances between the anchors
!> \param energy objective accumulator
!> \param grad Cartesian gradient accumulator
!> \param shortest unregularized minimum pair distance accumulator
! **************************************************************************************************
      SUBROUTINE evaluate_image(xyz, fraction, energy, grad, shortest)
      REAL(KIND=dp), INTENT(IN)                          :: xyz(:), fraction
      REAL(KIND=dp), INTENT(INOUT)                       :: energy, grad(:), shortest

      INTEGER                                            :: ia, ib, pair
      REAL(KIND=dp)                                      :: d, dg(3), dx(3), r, r2, weight

         pair = 0
         DO ia = 1, ns - 3, 3
            DO ib = ia + 3, ns, 3
               pair = pair + 1
               dx = pair_delta(xyz(ia:ia + 2) - xyz(ib:ib + 2))
               r2 = DOT_PRODUCT(dx, dx)
               shortest = MIN(shortest, SQRT(r2))
               r2 = r2 + soft_radius**2
               r = SQRT(r2)
               d = (1.0_dp - fraction)*anchor_distances(pair, 1) + fraction*anchor_distances(pair, 2)
               d = SQRT(d*d + soft_radius**2)
               weight = 1.0_dp/(r2*r2)
               energy = energy + weight*(r - d)**2
               dg = (2.0_dp*weight*(r - d)*(2.0_dp*d - r)/r2)*dx
               grad(ia:ia + 2) = grad(ia:ia + 2) + dg
               grad(ib:ib + 2) = grad(ib:ib + 2) - dg
            END DO
         END DO
      END SUBROUTINE evaluate_image

! **************************************************************************************************
!> \brief Break the symmetry of overlapping or collinearly exchanging pairs.
!> \param xyz image coordinates in Angstrom
! **************************************************************************************************
      SUBROUTINE perturb_short_contacts(xyz)
      REAL(KIND=dp), INTENT(INOUT)                       :: xyz(:)

      INTEGER                                            :: dim_min, ia, ib, pair
      REAL(KIND=dp)                                      :: delta(3), direction(3), other(3), r

         pair = 0
         DO ia = 1, ns - 3, 3
            DO ib = ia + 3, ns, 3
               pair = pair + 1
               delta = pair_delta(xyz(ia:ia + 2) - xyz(ib:ib + 2))
               r = NORM2(delta)
               IF (r >= soft_radius) THEN
                  IF (.NOT. sequential) CYCLE
                  IF (r >= 0.5_dp*MINVAL(anchor_distances(pair, :))) CYCLE
                  ! A collinear exchange can put the collision between images.
                  ! The nudged force then removes the repulsion along the tangent.
                  delta = pair_delta(path(ia:ia + 2, 1) - path(ib:ib + 2, 1))
                  other = pair_delta(path(ia:ia + 2, ni) - path(ib:ib + 2, ni))
                  IF (DOT_PRODUCT(delta, other) > -0.999_dp*PRODUCT(anchor_distances(pair, :))) CYCLE
               END IF
               delta = pair_delta(path(ia:ia + 2, 1) - path(ib:ib + 2, 1))
               dim_min = MINLOC(ABS(delta), 1)
               direction = 0.0_dp
               direction(dim_min) = 1.0_dp
               IF (NORM2(delta) > soft_radius) THEN
                  delta = delta/NORM2(delta)
                  direction = direction - DOT_PRODUCT(direction, delta)*delta
               END IF
               direction = 3.0_dp*soft_radius*direction/NORM2(direction)
               xyz(ia:ia + 2) = xyz(ia:ia + 2) + direction
               xyz(ib:ib + 2) = xyz(ib:ib + 2) - direction
            END DO
         END DO
      END SUBROUTINE perturb_short_contacts

! **************************************************************************************************
!> \brief Grow an S-IDPP segment inward and relax it with uniform springs once full.
! **************************************************************************************************
      SUBROUTINE grow_segment()
      INTEGER                                            :: front, i, left, m, right
      LOGICAL                                            :: converged, hit_limit, incomplete
      REAL(KIND=dp)                                      :: gap, shortest, spacing, step, tangent(ns)
      REAL(KIND=dp), ALLOCATABLE                         :: band(:, :), energies(:), fractions(:), &
                                                            springs(:)

         ALLOCATE (band(ns, ni), energies(ni), fractions(ni), springs(ni - 1))
         IF (ni == 3) THEN
            left = 2
            right = 3
            path(:, 2) = 0.5_dp*(path(:, 1) + path(:, ni))
            CALL perturb_short_contacts(path(:, 2))
         ELSE
            left = 2
            right = ni - 1
            path(:, left) = path(:, 1) + (path(:, ni) - path(:, 1))/REAL(ni - 1, dp)
            path(:, right) = path(:, ni) - (path(:, ni) - path(:, 1))/REAL(ni - 1, dp)
            CALL perturb_short_contacts(path(:, left))
            CALL perturb_short_contacts(path(:, right))
         END IF
         hit_limit = .FALSE.
         DO
            m = left + ni - right + 1
            incomplete = m < ni
            band(:, 1:left) = path(:, 1:left)
            band(:, left + 1:m) = path(:, right:ni)
            DO i = 1, left
               fractions(i) = REAL(i - 1, dp)/REAL(ni - 1, dp)
            END DO
            DO i = left + 1, m
               fractions(i) = REAL(i - left + right - 2, dp)/REAL(ni - 1, dp)
            END DO
            spacing = SUM(NORM2(band(:, 2:m) - band(:, 1:m - 1), DIM=1))/REAL(ni - 1, dp)
            springs(1:m - 1) = 1.0_dp
            IF (incomplete) THEN
               gap = NORM2(band(:, left + 1) - band(:, left))
               ! Scale the spring across the unfilled interval to maintain the target tension.
               springs(left) = spacing/MAX(soft_radius, gap, spacing)
            END IF
            CALL relax_band(band(:, 1:m), fractions(1:m), springs(1:m - 1), left, incomplete, &
                            energies(1:m), front, shortest, converged)
            hit_limit = hit_limit .OR. .NOT. converged
            path(:, 1:left) = band(:, 1:left)
            path(:, right:ni) = band(:, left + 1:m)
            IF (.NOT. incomplete) EXIT

            gap = NORM2(band(:, left + 1) - band(:, left))
            spacing = SUM(NORM2(band(:, 2:m) - band(:, 1:m - 1), DIM=1))/REAL(ni - 1, dp)
            step = MIN(spacing, 0.5_dp*gap)
            CALL band_tangent(band(:, front - 1:front + 1), energies(front - 1:front + 1), &
                              tangent, weighted=.TRUE.)
            IF (front == left) THEN
               left = left + 1
               path(:, left) = band(:, front) + step*tangent
               CALL perturb_short_contacts(path(:, left))
            ELSE
               right = right - 1
               path(:, right) = band(:, front) - step*tangent
               CALL perturb_short_contacts(path(:, right))
            END IF
         END DO
         IF (iw > 0) THEN
            WRITE (iw, '(T2,A,I0,A,F14.8)') "NEB| S-IDPP ", ni, " images; minimum pair distance [angstrom]: ", shortest
            IF (hit_limit) WRITE (iw, '(T2,A)') "NEB| An S-IDPP stage reached its stopping limit. Using the last relaxed images."
         END IF
      END SUBROUTINE grow_segment

! **************************************************************************************************
!> \brief Relax an active band with velocity projection on the nudged IDPP force.
!> \param band active images in Angstrom
!> \param fractions fractional position of each active image in the final segment
!> \param springs spring constants of the active intervals
!> \param left index of the growing front on the left
!> \param incomplete whether the band is still growing
!> \param energies objective values at the final step
!> \param front active image to grow from when the band is incomplete
!> \param shortest minimum pair distance among the active interior images
!> \param converged whether the stage reached its force criterion
! **************************************************************************************************
      SUBROUTINE relax_band(band, fractions, springs, left, incomplete, energies, front, shortest, converged)
      REAL(KIND=dp), INTENT(INOUT)                       :: band(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: fractions(:), springs(:)
      INTEGER, INTENT(IN)                                :: left
      LOGICAL, INTENT(IN)                                :: incomplete
      REAL(KIND=dp), INTENT(OUT)                         :: energies(:)
      INTEGER, INTENT(OUT)                               :: front
      REAL(KIND=dp), INTENT(OUT)                         :: shortest
      LOGICAL, INTENT(OUT)                               :: converged

      INTEGER, PARAMETER                                 :: max_steps = 1200
      REAL(KIND=dp), PARAMETER                           :: max_step = 0.04_dp, max_tol = 0.01_dp, &
                                                            rms_tol = 0.005_dp, timestep = 0.10_dp

      INTEGER                                            :: i, it, m
      REAL(KIND=dp)                                      :: backward, forward, min_image, power, &
                                                            scale, tangent(ns)
      REAL(KIND=dp), ALLOCATABLE                         :: force(:, :), residual(:), velocity(:, :)

         m = SIZE(band, 2)
         ALLOCATE (force(ns, m), residual(m), velocity(ns, m))
         velocity = 0.0_dp
         DO it = 0, max_steps
            shortest = HUGE(1.0_dp)
            force = 0.0_dp
            energies = 0.0_dp
            DO i = 1, m
               min_image = HUGE(1.0_dp)
               CALL evaluate_image(band(:, i), fractions(i), energies(i), force(:, i), min_image)
               IF (i > 1 .AND. i < m) shortest = MIN(shortest, min_image)
            END DO
            residual = 0.0_dp
            DO i = 2, m - 1
               CALL band_tangent(band(:, i - 1:i + 1), energies(i - 1:i + 1), tangent, weighted=.FALSE.)
               backward = NORM2(band(:, i) - band(:, i - 1))
               forward = NORM2(band(:, i + 1) - band(:, i))
               force(:, i) = -force(:, i) + (DOT_PRODUCT(force(:, i), tangent) + &
                                             springs(i)*forward - springs(i - 1)*backward)*tangent
               residual(i) = MAX(MAXVAL(ABS(force(:, i)))/max_tol, NORM2(force(:, i))/(SQRT(REAL(ns, dp))*rms_tol))
            END DO
            front = left
            IF (incomplete) THEN
               IF (residual(left + 1) < residual(left)) front = left + 1
               converged = residual(front) <= 1.0_dp
            ELSE
               converged = MAXVAL(residual) <= 1.0_dp
            END IF
            IF (converged .OR. it == max_steps) EXIT
            DO i = 2, m - 1
               velocity(:, i) = velocity(:, i) + timestep*force(:, i)
               power = DOT_PRODUCT(velocity(:, i), force(:, i))
               IF (power > 0.0_dp) THEN
                  velocity(:, i) = (power/DOT_PRODUCT(force(:, i), force(:, i)))*force(:, i)
                  scale = MIN(timestep, max_step/MAX(MAXVAL(ABS(velocity(:, i))), TINY(1.0_dp)))
                  band(:, i) = band(:, i) + scale*velocity(:, i)
                  velocity(:, i) = (scale/timestep)*velocity(:, i)
               ELSE
                  velocity(:, i) = 0.0_dp
               END IF
            END DO
         END DO
      END SUBROUTINE relax_band

! **************************************************************************************************
!> \brief Estimate the tangent using normalized neighboring segments as in S-IDPP.
!> \param band three consecutive images
!> \param energies their IDPP objective values
!> \param tangent normalized tangent at the middle image
!> \param weighted always use both sides for seeding new images
! **************************************************************************************************
      SUBROUTINE band_tangent(band, energies, tangent, weighted)
      REAL(KIND=dp), INTENT(IN)                          :: band(:, :), energies(3)
      REAL(KIND=dp), INTENT(OUT)                         :: tangent(:)
      LOGICAL, INTENT(IN)                                :: weighted

      REAL(KIND=dp)                                      :: backward(ns), forward(ns), norm, wback, &
                                                            wfwd

         backward = band(:, 2) - band(:, 1)
         forward = band(:, 3) - band(:, 2)
         backward = backward/MAX(NORM2(backward), TINY(1.0_dp))
         forward = forward/MAX(NORM2(forward), TINY(1.0_dp))
         IF (.NOT. weighted .AND. energies(3) > energies(2) .AND. energies(2) > energies(1)) THEN
            tangent = forward
         ELSE IF (.NOT. weighted .AND. energies(3) < energies(2) .AND. energies(2) < energies(1)) THEN
            tangent = backward
         ELSE
            wfwd = ABS(energies(3) - energies(2))
            wback = ABS(energies(2) - energies(1))
            IF (energies(3) > energies(1)) THEN
               tangent = MAX(wfwd, wback)*forward + MIN(wfwd, wback)*backward
            ELSE IF (energies(3) < energies(1)) THEN
               tangent = MIN(wfwd, wback)*forward + MAX(wfwd, wback)*backward
            ELSE
               tangent = forward + backward
            END IF
         END IF
         norm = NORM2(tangent)
         IF (norm <= TINY(1.0_dp)) THEN
            tangent = forward
            IF (NORM2(tangent) <= TINY(1.0_dp)) tangent = backward
         ELSE
            tangent = tangent/norm
         END IF
      END SUBROUTINE band_tangent

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
