!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM topology_snapshot_unittest
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE topology_snapshot,               ONLY: check_snapshot,&
                                              check_snapshot_seam,&
                                              snapshot_overlap,&
                                              snapshot_type

   IMPLICIT NONE
   REAL(KIND=dp), PARAMETER :: pi = 3.1415926535897932384626433832795_dp, tolerance = 1.e-12_dp
   TYPE(snapshot_type) :: a, b
   COMPLEX(KIND=dp) :: m(1, 1), reverse(1, 1), phase
   REAL(KIND=dp) :: metric, gap
   INTEGER :: i
   CALL init_orbital_pointers(1)
   a%nao = 1
   a%rank = 1
   a%nspin = 1
   a%channel = 1
   DO i = 1, 3
      a%cell(i, i) = 20.0_dp
      a%inverse(i, i) = 0.05_dp
   END DO
   ALLOCATE (a%atoms(1), a%bands(1), a%energies(2), a%coefficients(1, 1))
   a%bands(1) = 1
   a%energies(:) = [-1.0_dp, 1.0_dp]
   a%coefficients(:, :) = 1.0_dp
   a%atoms(1)%kind = 1
   a%atoms(1)%count = 1
   ALLOCATE (a%atoms(1)%shells(1))
   ASSOCIATE (s => a%atoms(1)%shells(1))
      s%first = 1
      s%count = 1
      s%radius = 5.0_dp
      ALLOCATE (s%exponent(1), s%radii(1), s%contraction(1, 1))
      s%exponent(1) = 1.0_dp
      s%radii(1) = 5.0_dp
      s%contraction(1, 1) = (2.0_dp/pi)**0.75_dp
   END ASSOCIATE
   CALL check_snapshot(a, tolerance, tolerance, metric, gap)
   IF (metric > tolerance .OR. ABS(gap - 2.0_dp) > tolerance) ERROR STOP 'Incorrect snapshot metric or gap'
   ! Derived-type assignment deliberately deep-copies allocatable frame components.
   b = a
   b%atoms(1)%center(1) = 0.7_dp
   b%atoms(1)%shells(1)%center(1) = 0.7_dp
   phase = EXP(CMPLX(0.0_dp, 0.3_dp, dp))
   b%coefficients(:, :) = phase
   CALL snapshot_overlap(a, b, a%k, b%k, m)
   IF (ABS(m(1, 1) - phase*EXP(-0.7_dp**2/2.0_dp)) > tolerance) ERROR STOP 'Moving Gaussian overlap failed'
   CALL snapshot_overlap(b, a, b%k, a%k, reverse)
   IF (ABS(m(1, 1) - CONJG(reverse(1, 1))) > tolerance) ERROR STOP 'Directed overlap adjoint failed'
   b%atoms(1)%center(:) = 0.0_dp
   b%atoms(1)%shells(1)%center(:) = 0.0_dp
   CALL check_snapshot_seam(a, b, [1], tolerance)
   CALL snapshot_overlap(a, b, a%k, b%k, m)
   IF (ABS(m(1, 1) - phase) > tolerance) ERROR STOP 'Endpoint gauge sewing failed'
END PROGRAM topology_snapshot_unittest
