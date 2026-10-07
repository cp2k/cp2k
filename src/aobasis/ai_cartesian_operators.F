!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Coordinate and basis-center derivative ladders for Cartesian Gaussian integrals.
!>        Derivatives act on a basis function with the integral operator held fixed.
!>        Electronic-coordinate derivatives have the opposite sign.
! **************************************************************************************************
MODULE ai_cartesian_operators
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: binomial_gen
   USE orbital_pointers,                ONLY: coset,&
                                              indco,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: cartesian_ladder_map, cartesian_ladder, apply_cartesian_operator, apply_cartesian_1d, &
             cartesian_polynomial_1d

CONTAINS

! **************************************************************************************************
!> \brief Construct base, raised, lowered and degree indices for a Cartesian ladder.
!> \param lmax highest output shell or one-dimensional degree
!> \param lmin lowest output shell (zero for one-dimensional tables)
!> \param axis Cartesian axis, or zero for a one-dimensional table
!> \param map rows contain base, raised, lowered and degree indices
!> \param source_lmin optional lowest input shell; defaults to lmin for Cartesian blocks
! **************************************************************************************************
   SUBROUTINE cartesian_ladder_map(lmax, lmin, axis, map, source_lmin)
      INTEGER, INTENT(IN)                                :: lmax, lmin, axis
      INTEGER, DIMENSION(:, :), INTENT(OUT)              :: map
      INTEGER, INTENT(IN), OPTIONAL                      :: source_lmin

      INTEGER                                            :: i, ico, lowered(3), offset, power(3), &
                                                            raised(3), source_offset

      IF (axis == 0) THEN
         DO i = 0, lmax
            map(:, i + 1) = [i + 1, i + 2, MAX(1, i), i]
         END DO
      ELSE
         offset = ncoset(lmin - 1)
         source_offset = offset
         IF (PRESENT(source_lmin)) source_offset = ncoset(source_lmin - 1)
         DO ico = offset + 1, ncoset(lmax)
            power = indco(:, ico)
            raised = power
            lowered = power
            raised(axis) = raised(axis) + 1
            lowered(axis) = MAX(0, lowered(axis) - 1)
            map(:, ico - offset) = [ico - source_offset, &
                                    coset(raised(1), raised(2), raised(3)) - source_offset, &
                                    MAX(1, coset(lowered(1), lowered(2), lowered(3)) - source_offset), power(axis)]
         END DO
      END IF
   END SUBROUTINE cartesian_ladder_map

! **************************************************************************************************
!> \brief Apply one coordinate or basis-center derivative ladder to indexed integrals.
!> \param source integral block including the raised functions
!> \param values transformed integral block
!> \param map base, raised, lowered and degree indices
!> \param derivative select a basis-center derivative instead of coordinate multiplication
!> \param exponent exponent of the differentiated Gaussian
!> \param origin displacement of the Gaussian center from the coordinate reference point
!> \param right apply the map to columns rather than rows
! **************************************************************************************************
   PURE SUBROUTINE cartesian_ladder(source, values, map, derivative, exponent, origin, right)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: source
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: values
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: map
      LOGICAL, INTENT(IN)                                :: derivative
      REAL(KIND=dp), INTENT(IN)                          :: exponent, origin
      LOGICAL, INTENT(IN)                                :: right

      INTEGER                                            :: i, j, second_index
      REAL(KIND=dp)                                      :: scale, second_scale(SIZE(map, 2))

      IF (derivative) THEN
         scale = 2.0_dp*exponent
         second_scale = -REAL(map(4, :), dp)
         second_index = 3
      ELSE
         scale = 1.0_dp
         second_scale = origin
         second_index = 1
      END IF

      ! Keep the first (contiguous) dimension innermost on both sides.
      IF (right) THEN
         DO j = 1, SIZE(map, 2)
            DO i = 1, SIZE(values, 1)
               values(i, j) = scale*source(i, map(2, j)) + &
                              second_scale(j)*source(i, map(second_index, j))
            END DO
         END DO
      ELSE
         DO j = 1, SIZE(values, 2)
            DO i = 1, SIZE(map, 2)
               values(i, j) = scale*source(map(2, i), j) + &
                              second_scale(i)*source(map(second_index, i), j)
            END DO
         END DO
      END IF
   END SUBROUTINE cartesian_ladder

! **************************************************************************************************
!> \brief Apply a Cartesian operator on either side of a primitive integral block.
!>        The input includes one extra shell on both sides. Coordinate multiplication
!>        uses the supplied shell-range layout; derivative inputs extend the lowest
!>        shell by one on the selected side for lowering. Screened primitive blocks
!>        must be zero in the input, including their raised and lowered functions.
!> \param la_max highest output shell on A
!> \param la_min lowest output shell on A
!> \param npgfa number of primitives on A
!> \param lb_max highest output shell on B
!> \param lb_min lowest output shell on B
!> \param npgfb number of primitives on B
!> \param ab input integral block with raised shells
!> \param values three Cartesian components in the output shell-range layout
!> \param right act on the ket rather than the bra
!> \param derivative select a basis-center derivative
!> \param exponents Gaussian exponents on the selected side, required for derivatives
!> \param origin displacement of the selected center from the coordinate reference point
! **************************************************************************************************
   SUBROUTINE apply_cartesian_operator(la_max, la_min, npgfa, lb_max, lb_min, npgfb, &
                                       ab, values, right, derivative, exponents, origin)
      INTEGER, INTENT(IN)                                :: la_max, la_min, npgfa, lb_max, lb_min, &
                                                            npgfb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: values
      LOGICAL, INTENT(IN)                                :: right, derivative
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: exponents
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: origin

      INTEGER                                            :: axis, ia, iap, ib, ibp, ipgf, jpgf, n, &
                                                            na, nap, nb, nbp, source_min_a, &
                                                            source_min_b
      INTEGER, ALLOCATABLE                               :: map(:, :, :)
      LOGICAL                                            :: valid_shape
      REAL(KIND=dp)                                      :: center(3), exponent

      na = ncoset(la_max) - ncoset(la_min - 1)
      nb = ncoset(lb_max) - ncoset(lb_min - 1)
      source_min_a = la_min
      source_min_b = lb_min
      IF (derivative) THEN
         IF (right) THEN
            source_min_b = MAX(0, lb_min - 1)
         ELSE
            source_min_a = MAX(0, la_min - 1)
         END IF
      END IF
      nap = ncoset(la_max + 1) - ncoset(source_min_a - 1)
      nbp = ncoset(lb_max + 1) - ncoset(source_min_b - 1)
      IF (derivative) THEN
         CPASSERT(PRESENT(exponents))
      END IF
      valid_shape = ALL(SHAPE(ab) >= [nap*npgfa, nbp*npgfb])
      CPASSERT(valid_shape)
      valid_shape = ALL(SHAPE(values) >= [na*npgfa, nb*npgfb, 3])
      CPASSERT(valid_shape)
      center = 0.0_dp
      IF (PRESENT(origin)) center = origin
      n = MERGE(nb, na, right)
      ALLOCATE (map(4, n, 3))
      DO axis = 1, 3
         IF (right) THEN
            CALL cartesian_ladder_map(lb_max, lb_min, axis, map(:, :, axis), source_min_b)
         ELSE
            CALL cartesian_ladder_map(la_max, la_min, axis, map(:, :, axis), source_min_a)
         END IF
      END DO

      values = 0.0_dp
      DO ipgf = 1, npgfa
         ia = (ipgf - 1)*na
         iap = (ipgf - 1)*nap
         DO jpgf = 1, npgfb
            ib = (jpgf - 1)*nb
            ibp = (jpgf - 1)*nbp
            exponent = 0.0_dp
            IF (derivative) exponent = exponents(MERGE(jpgf, ipgf, right))
            DO axis = 1, 3
               CALL cartesian_ladder(ab(iap + 1:iap + nap, ibp + 1:ibp + nbp), &
                                     values(ia + 1:ia + na, ib + 1:ib + nb, axis), map(:, :, axis), &
                                     derivative, exponent, center(axis), right)
            END DO
         END DO
      END DO
   END SUBROUTINE apply_cartesian_operator

! **************************************************************************************************
!> \brief Apply the common ladder to separable one-dimensional integral tables.
!> \param source input tables indexed by the degrees on A and B and the Cartesian axis
!> \param values transformed tables; entries outside the requested range are untouched
!> \param la_max highest output degree on A
!> \param lb_max highest output degree on B
!> \param right act on B rather than A
!> \param derivative select a basis-center derivative
!> \param exponent Gaussian exponent on the selected side
!> \param origin displacement of that center from the coordinate reference point
! **************************************************************************************************
   SUBROUTINE apply_cartesian_1d(source, values, la_max, lb_max, right, derivative, exponent, origin)
      REAL(KIND=dp), DIMENSION(0:, 0:, :), INTENT(IN)    :: source
      REAL(KIND=dp), DIMENSION(0:, 0:, :), INTENT(INOUT) :: values
      INTEGER, INTENT(IN)                                :: la_max, lb_max
      LOGICAL, INTENT(IN)                                :: right, derivative
      REAL(KIND=dp), INTENT(IN)                          :: exponent, origin(3)

      INTEGER :: axis, map(4, MERGE(lb_max, la_max, right) + 1)

      CALL cartesian_ladder_map(MERGE(lb_max, la_max, right), 0, 0, map)
      DO axis = 1, 3
         CALL cartesian_ladder(source(:, :, axis), values(0:la_max, 0:lb_max, axis), map, &
                               derivative, exponent, origin(axis), right)
      END DO
   END SUBROUTINE apply_cartesian_1d

! **************************************************************************************************
!> \brief Form requested coordinate powers from separable Cartesian integral tables.
!>        Only requested powers are stored. With the origin on A, each power is a
!>        direct shift of the bra degree; other origins use the binomial expansion.
!> \param source input tables with enough raised bra degrees for the largest power
!> \param values tables indexed by bra degree, ket degree, requested power and axis
!> \param la_max highest output degree on A
!> \param lb_max highest output degree on B
!> \param powers nonnegative coordinate powers to evaluate
!> \param origin displacement of A from the coordinate reference point
! **************************************************************************************************
   SUBROUTINE cartesian_polynomial_1d(source, values, la_max, lb_max, powers, origin)
      REAL(KIND=dp), DIMENSION(0:, 0:, :), INTENT(IN)    :: source
      REAL(KIND=dp), DIMENSION(0:, 0:, 0:, :), &
         INTENT(OUT)                                     :: values
      INTEGER, INTENT(IN)                                :: la_max, lb_max, powers(0:)
      REAL(KIND=dp), INTENT(IN)                          :: origin(3)

      INTEGER                                            :: axis, i, k, p
      REAL(KIND=dp)                                      :: coefficient

      DO axis = 1, 3
         DO i = 0, SIZE(powers) - 1
            p = powers(i)
            IF (ABS(origin(axis)) <= 0.0_dp) THEN
               values(:, :, i, axis) = source(p:p + la_max, 0:lb_max, axis)
            ELSE
               values(:, :, i, axis) = 0.0_dp
               DO k = 0, p
                  coefficient = binomial_gen(REAL(p, dp), k)*origin(axis)**(p - k)
                  values(:, :, i, axis) = values(:, :, i, axis) + &
                                          coefficient*source(k:k + la_max, 0:lb_max, axis)
               END DO
            END IF
         END DO
      END DO
   END SUBROUTINE cartesian_polynomial_1d

END MODULE ai_cartesian_operators
