!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM skala_atom_grid_unittest
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: pi
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_comm_type,&
                                              mp_comm_world
   USE skala_atom_grid,                 ONLY: skala_atom_grid_build,&
                                              skala_atom_grid_direction,&
                                              skala_atom_grid_pullback,&
                                              skala_partition_init,&
                                              skala_partition_pullback,&
                                              skala_partition_type,&
                                              skala_partition_weights,&
                                              skala_radial_grid

   IMPLICIT NONE
   TYPE(mp_comm_type) :: comm, single
   TYPE(skala_partition_type) :: partition, translated, displaced
   REAL(dp) :: lattice(3, 3), inverse(3, 3), centres(3, 2), shifted(3, 2), x(3), offset(3), &
               p, d, p2, d2, total, scalar(4), serial(4), coarse_error, error, &
               radius(12), radial_weight(12), t(12), weight(12), moments(24)
   REAL(dp), ALLOCATABLE :: points(:, :), ew(:), dw(:), bc(:, :)
   INTEGER, ALLOCATABLE :: bs(:)
   INTEGER :: i, j, k, a, point, level, nr, na, blocks, empty, descriptor_only

   CALL cp2k_init()
   comm = mp_comm_world
   CALL comm%init()
   single = mp_comm_self
   CALL single%init()
   CALL skala_radial_grid(12, 1.7_dp, radius, radial_weight)
   t = radius/(radius + 1.7_dp)
   weight = radial_weight*(1.0_dp - t)**2/(radius**2*1.7_dp)
   DO i = 0, 23
      moments(i + 1) = SUM(weight*t**i) - 1.0_dp/REAL(i + 1, dp)
   END DO
   IF (MAXVAL(ABS(moments)) > 2.0e-14_dp) ERROR STOP "Gauss-Legendre moments or radial Jacobian are incorrect"
   CALL radial_grid_reuse()

   lattice = RESHAPE([8.0_dp, 0.0_dp, 0.0_dp, 1.3_dp, 7.0_dp, 0.0_dp, -0.7_dp, 0.8_dp, 9.0_dp], [3, 3])
   inverse = inv_3x3(lattice)
   centres(:, 1) = MATMUL(lattice, [0.13_dp, 0.31_dp, 0.21_dp])
   centres(:, 2) = MATMUL(lattice, [0.68_dp, 0.73_dp, 0.49_dp])
   CALL skala_partition_init(partition, lattice, centres, 1)
   shifted = centres
   shifted(:, 2) = centres(:, 2) + MATMUL(lattice, [2.0_dp, -1.0_dp, 1.0_dp])
   CALL skala_partition_init(translated, lattice, shifted, 1)
   shifted(:, 2) = centres(:, 2) + [0.2_dp, -0.3_dp, 0.1_dp]
   CALL skala_partition_init(displaced, lattice, shifted, 1)
   DO point = 1, 6
      x = MATMUL(lattice, [0.21_dp*point, 0.11_dp*point, 0.07_dp*point])
      IF (point == 5) x = centres(:, 1) + lattice(:, 1)
      IF (point == 6) x = centres(:, 1) + 2*lattice(:, 1) - [1.0e-7_dp, 0.0_dp, 0.0_dp]
      total = 0.0_dp
      DO a = 1, 2
         DO k = -3, 3
            DO j = -3, 3
               DO i = -3, 3
                  offset = x - centres(:, a) - MATMUL(lattice, REAL([i, j, k], dp))
                  CALL skala_partition_weights(partition, a, offset, p, d)
                  IF (p < 0.0_dp .OR. p > 1.0_dp .OR. d < 0.0_dp .OR. d > 1.0_dp) THEN
                     ERROR STOP "Invalid partition weight"
                  END IF
                  total = total + p
               END DO
            END DO
         END DO
      END DO
      IF (ABS(total - 1.0_dp) > 5.0e-13_dp) ERROR STOP "Periodic partition is not a partition of unity"
      offset = x - centres(:, 1)
      CALL skala_partition_weights(partition, 1, offset, p, d)
      CALL skala_partition_weights(translated, 1, offset, p2, d2)
      IF (ABS(p - p2) + ABS(d - d2) > 5.0e-13_dp) ERROR STOP "Rewrapping an atom changes the partition"
      CALL skala_partition_weights(displaced, 1, offset, p2, d2)
      IF (ABS(d - d2) > 5.0e-13_dp) ERROR STOP "Foreign atoms change the descriptor domain"
   END DO
   CALL skala_partition_weights(partition, 1, centres(:, 2) - centres(:, 1), p, d)
   IF (p /= 0.0_dp .OR. d <= 0.0_dp) ERROR STOP "A foreign nucleus must remain a descriptor-only target"
   CALL partition_derivatives(centres)
   CALL partition_derivatives(centres(:, :1))
   CALL quadrature_derivatives()
   CALL quadrature_direction()
   CALL cubic_quadrature_derivatives()

   coarse_error = HUGE(1.0_dp)
   DO level = 1, 2
      nr = 40
      na = 110
      IF (level == 2) THEN
         nr = 160
         na = 974
      END IF
      CALL skala_atom_grid_build(comm, lattice, centres, nr, na, 1, 1.0_dp, points, ew, dw, bs, bc)
      IF (SUM(bs) /= SIZE(ew) .OR. SIZE(dw) /= SIZE(ew)) ERROR STOP "Incomplete atom blocks"
      blocks = SIZE(bs)
      empty = 0
      IF (SIZE(ew) == 0) empty = 1
      descriptor_only = COUNT(ew == 0.0_dp .AND. dw > 0.0_dp)
      CALL comm%sum(blocks)
      CALL comm%sum(empty)
      CALL comm%sum(descriptor_only)
      IF (blocks /= 2 .OR. empty /= MAX(0, comm%num_pe - 2)) ERROR STOP "Atom-block ownership is not unique"
      scalar = integrals(points, ew)
      CALL comm%sum(scalar)
      ! The exact integrals are volume, 0, 0, and volume/2.
      serial = [ABS(det_3x3(lattice)), 0.0_dp, 0.0_dp, 0.5_dp*ABS(det_3x3(lattice))]
      error = MAXVAL(ABS(scalar - serial))/ABS(det_3x3(lattice))
      IF (comm%mepos == 0) WRITE (*, '(A,2I6,5ES16.7)') "Atom-grid integration (nr,na,integrals,error): ", nr, na, scalar, error
      IF (level == 1) coarse_error = error
      IF (level == 2 .AND. (error >= coarse_error .OR. error > 1.0e-6_dp)) THEN
         ERROR STOP "Periodic volume/Fourier quadrature did not converge"
      END IF
      IF (level == 1) THEN
         CALL skala_atom_grid_build(single, lattice, centres, nr, na, 1, 1.0_dp, points, ew, dw, bs, bc)
         serial = integrals(points, ew)
         IF (MAXVAL(ABS(scalar - serial)) > 1.0e-10_dp) ERROR STOP "Distributed atom-grid integrals differ"
      END IF
   END DO
   shifted = centres
   CALL skala_radial_grid(12, 1.0_dp, radius, radial_weight)
   shifted(:, 2) = centres(:, 1) + [radius(6), 0.0_dp, 0.0_dp]
   CALL skala_atom_grid_build(comm, lattice, shifted, 12, 6, 1, 1.0_dp, points, ew, dw, bs, bc)
   descriptor_only = COUNT(ew == 0.0_dp .AND. dw > 0.0_dp)
   CALL comm%sum(descriptor_only)
   IF (descriptor_only == 0) ERROR STOP "Descriptor-only quadrature rows have been discarded"
   IF (comm%mepos == 0) WRITE (*, '(A)') "Periodic Skala atom-grid checks passed"
   CALL cp2k_finalize()
CONTAINS
! **************************************************************************************************
!> \brief Check radial-rule keys, returned-array independence and threaded reuse.
! **************************************************************************************************
   SUBROUTINE radial_grid_reuse()
      INTEGER                                            :: trial
      REAL(dp) :: other_radius(13), other_weights(13), radius(12), reference_other_radius(13), &
         reference_other_weight(13), reference_radius(12), reference_weight(12), weights(12)

      CALL skala_radial_grid(12, 1.7_dp, reference_radius, reference_weight)
      CALL skala_radial_grid(13, 1.7_dp, reference_other_radius, reference_other_weight)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(reference_radius, reference_weight, reference_other_radius, reference_other_weight) &
!$OMP PRIVATE(trial, radius, weights, other_radius, other_weights)
      DO trial = 1, 8
         IF (MOD(trial, 2) == 0) THEN
            CALL skala_radial_grid(13, 1.7_dp, other_radius, other_weights)
            IF (ANY(other_radius /= reference_other_radius) .OR. ANY(other_weights /= reference_other_weight)) THEN
               ERROR STOP "Radial-rule point count is not respected"
            END IF
         ELSE
            CALL skala_radial_grid(12, 0.85_dp, radius, weights)
            IF (MAXVAL(ABS(2.0_dp*radius - reference_radius)) > 16*EPSILON(1.0_dp)*MAXVAL(reference_radius) .OR. &
                MAXVAL(ABS(8.0_dp*weights - reference_weight)) > 16*EPSILON(1.0_dp)*MAXVAL(reference_weight)) THEN
               ERROR STOP "Radial-rule scale is not respected"
            END IF
         END IF
         radius = -1.0_dp
         weights = -1.0_dp
         CALL skala_radial_grid(12, 1.7_dp, radius, weights)
         IF (ANY(radius /= reference_radius) .OR. ANY(weights /= reference_weight)) THEN
            ERROR STOP "Radial-rule reuse changes points or weights"
         END IF
      END DO
!$OMP END PARALLEL DO
   END SUBROUTINE radial_grid_reuse

! **************************************************************************************************
!> \brief Check periodic partition pullbacks against directional finite differences.
!> \param positions Cartesian positions, in bohr.
! **************************************************************************************************
   SUBROUTINE partition_derivatives(positions)
      REAL(dp), INTENT(IN)                               :: positions(:, :)

      INTEGER                                            :: atom, case_id, component, derivative, &
                                                            iteration, trial
      REAL(dp) :: ad, adjoint(2), dh(3, 3), dmw, dpw, dx(3), err, fd, gh(3, 3), gx(3), h, hi, &
         hm(3, 3), hp(3, 3), lo, maximum, mid, pm, pp, shifted_gh(3, 3), shifted_gx(3), xm(3), &
         xp(3), shifted_gr(3, SIZE(positions, 2)), rp(3, SIZE(positions, 2)), &
         rm(3, SIZE(positions, 2)), gr(3, SIZE(positions, 2)), dr(3, SIZE(positions, 2))
      TYPE(skala_partition_type)                         :: minus_partition, plus_partition, &
                                                            reference_partition

      dx = [0.13_dp, -0.17_dp, 0.11_dp]
      DO atom = 1, SIZE(positions, 2)
         dr(:, atom) = 0.2_dp*SIN(REAL(atom + [1, 2, 3], dp))
      END DO
      dh = RESHAPE([0.13_dp, 0.09_dp, -0.17_dp, 0.22_dp, -0.19_dp, 0.16_dp, 0.11_dp, -0.21_dp, 0.18_dp], [3, 3])
      CALL skala_partition_init(reference_partition, lattice, positions, 1)
      maximum = 0.0_dp
      DO atom = 1, SIZE(positions, 2)
         DO case_id = 1, 4
            SELECT CASE (case_id)
            CASE (1)
               x = MATMUL(lattice, [0.37_dp, 0.21_dp, 0.19_dp])
            CASE (2)
               x = positions(:, MOD(atom, SIZE(positions, 2)) + 1) - positions(:, atom)
            CASE (3)
               ! Exercise the narrow self-image descriptor turn-on rather than its flat region.
               lo = 0.0_dp
               hi = 1.0_dp
               DO iteration = 1, 55
                  mid = 0.5_dp*(lo + hi)
                  CALL skala_partition_weights(reference_partition, atom, mid*lattice(:, 1), pp, dpw)
                  IF (dpw > 0.5_dp) THEN
                     lo = mid
                  ELSE
                     hi = mid
                  END IF
               END DO
               x = mid*lattice(:, 1)
               IF (ABS(dpw - 0.5_dp) > 1.0e-3_dp) ERROR STOP "Descriptor transition fixture is not resolved"
            CASE (4)
               x = MATMUL(lattice, [2.1_dp, 0.0_dp, 0.0_dp])
            END SELECT
            DO derivative = 1, 2
               adjoint = 0.0_dp
               adjoint(derivative) = 1.0_dp
               CALL skala_partition_pullback(reference_partition, atom, x, adjoint(1), adjoint(2), gx, gr, gh)
               IF (MAXVAL(ABS(SUM(gr, DIM=2))) > 1.0e-10_dp) THEN
                  ERROR STOP "Partition centre derivatives violate translation invariance"
               END IF
               IF (derivative == 2 .AND. MAXVAL(ABS(gr)) > 1.0e-10_dp) THEN
                  ERROR STOP "Self-image descriptor depends on independent Cartesian atom centres"
               END IF
               IF (case_id == 1) THEN
                  rp = positions
                  rp(:, SIZE(positions, 2)) = rp(:, SIZE(positions, 2)) + MATMUL(lattice, [2.0_dp, -1.0_dp, 1.0_dp])
                  CALL skala_partition_init(plus_partition, lattice, rp, 1)
                  CALL skala_partition_pullback(plus_partition, atom, x, adjoint(1), adjoint(2), &
                                                shifted_gx, shifted_gr, shifted_gh)
                  shifted_gh = shifted_gh + SPREAD(gr(:, SIZE(positions, 2)), 2, 3)* &
                               SPREAD([2.0_dp, -1.0_dp, 1.0_dp], 1, 3)
                  IF (MAX(MAXVAL(ABS(gx - shifted_gx)), MAXVAL(ABS(gr - shifted_gr)), &
                          MAXVAL(ABS(gh - shifted_gh))) > 1.0e-10_dp) THEN
                     ERROR STOP "Rewrapping an atom changes partition geometry covectors"
                  END IF
               END IF
               DO component = 1, 3
                  ad = 0.0_dp
                  IF (component == 1) ad = SUM(gx*dx)
                  IF (component == 2) ad = SUM(gr*dr)
                  IF (component == 3) ad = SUM(gh*dh)
                  DO trial = 1, 2
                     h = 1.0e-5_dp/REAL(trial, dp)
                     xp = x
                     xm = x
                     rp = positions
                     rm = positions
                     hp = lattice
                     hm = lattice
                     SELECT CASE (component)
                     CASE (1)
                        xp = x + h*dx
                        xm = x - h*dx
                     CASE (2)
                        rp = positions + h*dr
                        rm = positions - h*dr
                     CASE (3)
                        hp = lattice + h*dh
                        hm = lattice - h*dh
                     END SELECT
                     CALL skala_partition_init(plus_partition, hp, rp, 1)
                     CALL skala_partition_init(minus_partition, hm, rm, 1)
                     CALL skala_partition_weights(plus_partition, atom, xp, pp, dpw)
                     CALL skala_partition_weights(minus_partition, atom, xm, pm, dmw)
                     fd = (adjoint(1)*(pp - pm) + adjoint(2)*(dpw - dmw))/(2.0_dp*h)
                     err = ABS(ad - fd)/MAX(1.0_dp, ABS(ad))
                     maximum = MAX(maximum, err)
                     IF (err > 2.0e-6_dp) THEN
                        WRITE (*, '(A,5I3,3ES18.8)') "Partition FD (atom,case,weight,component,step,AD,FD,error): ", &
                           atom, case_id, derivative, component, trial, ad, fd, err
                        ERROR STOP "Periodic partition derivative fails finite differences"
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
      IF (comm%mepos == 0) WRITE (*, '(A,ES12.4)') "Maximum partition derivative error: ", maximum
   END SUBROUTINE partition_derivatives

! **************************************************************************************************
!> \brief Check atom-grid coordinate and weight pullbacks.
! **************************************************************************************************
   SUBROUTINE quadrature_derivatives()
      INTEGER                                            :: component, trial
      REAL(dp)                                           :: ad, dh(3, 3), dr(3, 2), e, em, ep, err, &
                                                            es, gh(3, 3), gr(3, 2), h, maximum, &
                                                            sh(3, 3), sr(3, 2)

      dr = RESHAPE([0.21_dp, -0.11_dp, 0.09_dp, -0.17_dp, 0.23_dp, 0.14_dp], [3, 2])
      dh = RESHAPE([0.13_dp, 0.09_dp, -0.17_dp, 0.22_dp, -0.19_dp, 0.16_dp, 0.11_dp, -0.21_dp, 0.18_dp], [3, 3])
      e = grid_energy(comm, lattice, centres, gr, gh)
      es = grid_energy(single, lattice, centres, sr, sh)
      IF (MAX(ABS(e - es), MAXVAL(ABS(gr - sr)), MAXVAL(ABS(gh - sh))) > 1.0e-10_dp) THEN
         ERROR STOP "Distributed quadrature geometry covectors differ"
      END IF
      maximum = 0.0_dp
      DO component = 1, 2
         DO trial = 1, 2
            h = 1.0e-4_dp/REAL(trial, dp)
            IF (component == 1) THEN
               ad = SUM(gr*dr)
               ep = grid_energy(comm, lattice, centres + h*dr)
               em = grid_energy(comm, lattice, centres - h*dr)
            ELSE
               ad = SUM(gh*dh)
               ep = grid_energy(comm, lattice + h*dh, centres)
               em = grid_energy(comm, lattice - h*dh, centres)
            END IF
            err = ABS(ad - (ep - em)/(2.0_dp*h))/MAX(1.0_dp, ABS(ad))
            maximum = MAX(maximum, err)
            IF (err > 2.0e-6_dp) THEN
               WRITE (*, '(A,2I3,3ES18.8)') "Grid geometry FD (component,step,AD,FD,error): ", &
                  component, trial, ad, (ep - em)/(2.0_dp*h), err
               ERROR STOP "Atom-grid geometry derivative fails finite differences"
            END IF
         END DO
      END DO
      IF (comm%mepos == 0) WRITE (*, '(A,ES12.4)') "Maximum quadrature geometry derivative error: ", maximum
   END SUBROUTINE quadrature_derivatives

! **************************************************************************************************
!> \brief Check forward geometry tangents and their adjoint identity.
! **************************************************************************************************
   SUBROUTINE quadrature_direction()
      INTEGER                                            :: component, row, trial
      INTEGER, ALLOCATABLE                               :: minus_sizes(:), plus_sizes(:), sizes(:)
      REAL(dp)                                           :: ad, ah(3, 3), ar(3, 2), dh(3, 3), &
                                                            dr(3, 2), error, gh(3, 3), gr(3, 2), &
                                                            h, maximum, reverse
      REAL(dp), ALLOCATABLE :: ab(:, :), adesc(:), ae(:), ax(:, :), blocks(:, :), d(:), db(:, :), &
         dd(:), de(:), dx(:, :), e(:), mb(:, :), md(:), me(:), mx(:, :), pb(:, :), pd(:), pe(:), &
         px(:, :), xyz(:, :)

      CALL skala_atom_grid_build(comm, lattice, centres, 12, 14, 1, 1.0_dp, xyz, e, d, sizes, blocks)
      ALLOCATE (dx(3, SIZE(e)), de(SIZE(e)), dd(SIZE(e)), db(3, SIZE(sizes)))
      ALLOCATE (ax(3, SIZE(e)), ae(SIZE(e)), adesc(SIZE(e)), ab(3, SIZE(sizes)))
      ax = SIN(xyz)
      ae = COS(SUM(xyz, DIM=1))
      adesc = SIN(SUM(xyz, DIM=1))
      ab = COS(blocks)
      CALL skala_atom_grid_pullback(comm, lattice, centres, 12, 14, 1, 1.0_dp, sizes, &
                                    ax, ae, adesc, ab, gr, gh)
      dr = RESHAPE([0.21_dp, -0.11_dp, 0.09_dp, -0.17_dp, 0.23_dp, 0.14_dp], [3, 2])
      dh = RESHAPE([0.13_dp, 0.09_dp, -0.17_dp, 0.22_dp, -0.19_dp, 0.16_dp, 0.11_dp, -0.21_dp, 0.18_dp], [3, 3])
      maximum = 0.0_dp
      DO component = 1, 3
         ar = dr
         ah = dh
         IF (component == 1) ah = 0.0_dp
         IF (component == 2) ar = 0.0_dp
         CALL skala_atom_grid_direction(comm, lattice, centres, 12, 14, 1, 1.0_dp, sizes, &
                                        ar, ah, dx, de, dd, db)
         ad = SUM(ax*dx) + SUM(ae*de) + SUM(adesc*dd) + SUM(ab*db)
         CALL comm%sum(ad)
         reverse = SUM(gr*ar) + SUM(gh*ah)
         IF (ABS(ad - reverse)/MAX(1.0_dp, ABS(ad), ABS(reverse)) > 2.0e-12_dp) THEN
            ERROR STOP "Atom-grid forward and reverse geometry derivatives differ"
         END IF
         DO trial = 1, 2
            h = 1.0e-4_dp/REAL(trial, dp)
            CALL skala_atom_grid_build(comm, lattice + h*ah, centres + h*ar, 12, 14, 1, 1.0_dp, &
                                       px, pe, pd, plus_sizes, pb)
            CALL skala_atom_grid_build(comm, lattice - h*ah, centres - h*ar, 12, 14, 1, 1.0_dp, &
                                       mx, me, md, minus_sizes, mb)
            IF (SIZE(plus_sizes) /= SIZE(sizes) .OR. SIZE(minus_sizes) /= SIZE(sizes)) THEN
               ERROR STOP "Atom-grid direction changes block ownership"
            END IF
            IF (ANY(plus_sizes /= sizes) .OR. ANY(minus_sizes /= sizes)) THEN
               ERROR STOP "Atom-grid direction fixture changes row topology"
            END IF
            error = MAX(0.0_dp, MAXVAL(ABS(dx - (px - mx)/(2*h))), MAXVAL(ABS(db - (pb - mb)/(2*h))))
            DO row = 1, SIZE(e)
               error = MAX(error, ABS(de(row) - (pe(row) - me(row))/(2*h))/MAX(1.0_dp, ABS(de(row))), &
                           ABS(dd(row) - (pd(row) - md(row))/(2*h))/MAX(1.0_dp, ABS(dd(row))))
            END DO
            CALL comm%max(error)
            maximum = MAX(maximum, error)
            IF (error > 2.0e-6_dp) ERROR STOP "Atom-grid forward derivative fails finite differences"
         END DO
      END DO
      DO component = 1, 2
         ar(:, component) = [0.04_dp, -0.02_dp, 0.03_dp]
      END DO
      ah = 0.0_dp
      CALL skala_atom_grid_direction(comm, lattice, centres, 12, 14, 1, 1.0_dp, sizes, ar, ah, dx, de, dd, db)
      error = MAX(0.0_dp, MAXVAL(ABS(de)), MAXVAL(ABS(dd)))
      CALL comm%max(error)
      IF (error > 2.0e-12_dp) ERROR STOP "Rigid translation changes atom-grid weights"
      IF (comm%mepos == 0) WRITE (*, '(A,ES12.4)') "Maximum quadrature forward derivative error: ", maximum
   END SUBROUTINE quadrature_direction

! **************************************************************************************************
!> \brief Check cell and atom derivatives in a cubic cell.
! **************************************************************************************************
   SUBROUTINE cubic_quadrature_derivatives()
      INTEGER                                            :: component, trial
      REAL(dp)                                           :: ad, dh(3, 3), dr(3, 1), em, ep, err, &
                                                            gh(3, 3), gr(3, 1), h, hmat(3, 3), &
                                                            maximum, positions(3, 1), reference

      hmat = 0.0_dp
      DO component = 1, 3
         hmat(component, component) = 6.0_dp
      END DO
      positions(:, 1) = MATMUL(hmat, [0.3_dp, 0.4_dp, 0.3_dp])
      dr(:, 1) = [0.21_dp, -0.11_dp, 0.09_dp]
      dh = RESHAPE([0.13_dp, 0.09_dp, -0.17_dp, 0.22_dp, -0.19_dp, 0.16_dp, 0.11_dp, -0.21_dp, 0.18_dp], [3, 3])
      reference = grid_energy(comm, hmat, positions, gr, gh, 100, 110)
      maximum = 0.0_dp
      DO component = 1, 2
         DO trial = 1, 2
            h = 1.0e-4_dp/REAL(trial, dp)
            IF (component == 1) THEN
               ad = SUM(gr*dr)
               ep = grid_energy(comm, hmat, positions + h*dr, nradial=100, nangular=110)
               em = grid_energy(comm, hmat, positions - h*dr, nradial=100, nangular=110)
            ELSE
               ad = SUM(gh*dh)
               ep = grid_energy(comm, hmat + h*dh, positions, nradial=100, nangular=110)
               em = grid_energy(comm, hmat - h*dh, positions, nradial=100, nangular=110)
            END IF
            err = ABS(ad - (ep - em)/(2*h))/MAX(1.0_dp, ABS(ad))
            maximum = MAX(maximum, err)
            IF (err > 2.0e-6_dp) ERROR STOP "Cubic atom-grid geometry derivative fails finite differences"
         END DO
      END DO
      IF (comm%mepos == 0) THEN
         WRITE (*, '(A,ES12.4)') "Cubic quadrature reference energy: ", reference
         WRITE (*, '(A,ES12.4)') "Maximum cubic quadrature derivative error: ", maximum
      END IF
   END SUBROUTINE cubic_quadrature_derivatives

! **************************************************************************************************
!> \brief Evaluate a smooth test integral and its geometry derivatives.
!> \param group Context communicator.
!> \param hmat Cell vectors as columns, in bohr.
!> \param positions Cartesian positions, in bohr.
!> \param gr Atom-coordinate derivative.
!> \param gh Cell-vector derivative.
!> \param nradial Number of radial nodes.
!> \param nangular Lebedev point count.
!> \return Test integral over the context communicator.
! **************************************************************************************************
   FUNCTION grid_energy(group, hmat, positions, gr, gh, nradial, nangular) RESULT(e)
      TYPE(mp_comm_type), INTENT(IN)                     :: group
      REAL(dp), INTENT(IN)                               :: hmat(3, 3), positions(:, :)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: gr(:, :), gh(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: nradial, nangular
      REAL(dp)                                           :: e

      INTEGER                                            :: na, nr, row
      INTEGER, ALLOCATABLE                               :: sizes(:)
      REAL(dp), ALLOCATABLE                              :: blocks(:, :), descriptor_weights(:), &
                                                            energy_weights(:), pd(:), pe(:), &
                                                            px(:, :), xyz(:, :)

      nr = 12
      na = 14
      IF (PRESENT(nradial)) nr = nradial
      IF (PRESENT(nangular)) na = nangular
      CALL skala_atom_grid_build(group, hmat, positions, nr, na, 1, 1.0_dp, &
                                 xyz, energy_weights, descriptor_weights, sizes, blocks)
      ALLOCATE (pe(SIZE(energy_weights)), pd(SIZE(energy_weights)), px(3, SIZE(energy_weights)))
      DO row = 1, SIZE(energy_weights)
         pe(row) = EXP(-SUM(xyz(:, row)**2)/20.0_dp)
         pd(row) = 0.3_dp*EXP(-SUM(xyz(:, row)**2)/30.0_dp)
         px(:, row) = -xyz(:, row)*(energy_weights(row)*pe(row)/10.0_dp + descriptor_weights(row)*pd(row)/15.0_dp)
      END DO
      e = SUM(energy_weights*pe + descriptor_weights*pd) + 0.02_dp*SUM(blocks**2)
      CALL group%sum(e)
      IF (PRESENT(gr) .AND. PRESENT(gh)) THEN
         CALL skala_atom_grid_pullback(group, hmat, positions, nr, na, 1, 1.0_dp, sizes, &
                                       px, pe, pd, 0.04_dp*blocks, gr, gh)
      END IF
   END FUNCTION grid_energy

! **************************************************************************************************
!> \brief Integrate constant and Fourier test fields.
!> \param positions Cartesian positions, in bohr.
!> \param weights Quadrature integration weights.
!> \return Constant and Fourier integrals.
! **************************************************************************************************
   FUNCTION integrals(positions, weights) RESULT(values)
      REAL(dp), INTENT(IN)                               :: positions(:, :), weights(:)
      REAL(dp)                                           :: values(4)

      INTEGER                                            :: row
      REAL(dp)                                           :: f, s(3)

      values = 0.0_dp
      DO row = 1, SIZE(weights)
         s = MATMUL(inverse, positions(:, row))
         f = COS(2*pi*(s(1) - s(2)))
         values = values + weights(row)*[1.0_dp, COS(2*pi*s(1)), SIN(2*pi*s(2)), f*f]
      END DO
   END FUNCTION integrals
END PROGRAM skala_atom_grid_unittest
