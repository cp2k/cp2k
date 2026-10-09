!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Unit tests for the detection of degenerate orbital sets split by an active space
! **************************************************************************************************
PROGRAM qs_active_space_utils_unittest
   USE kinds,                           ONLY: dp
   USE qs_active_space_utils,           ONLY: find_degenerate_splits
#include "./base/base_uses.f90"

   IMPLICIT NONE

   REAL(KIND=dp), PARAMETER                           :: eps = 1.0E-5_dp

   INTEGER                                            :: first_active, first_other, npairs
   REAL(KIND=dp), DIMENSION(8)                        :: energies

   ! orbital energies: 1 | 2 3 4 (triplet) | 5 6 7 (triplet) | 8
   energies = [-0.44_dp, -0.17_dp, -0.17_dp, -0.17_dp, 0.15_dp, 0.15_dp, 0.15_dp, 0.41_dp]

   ! complete degenerate sets: no split
   CALL find_degenerate_splits(energies, [2, 3, 4, 5, 6, 7], eps, npairs, first_active, first_other)
   CALL check(npairs == 0 .AND. first_active == 0 .AND. first_other == 0, "complete sets")

   ! single orbitals that are not degenerate with anything
   CALL find_degenerate_splits(energies, [1, 8], eps, npairs, first_active, first_other)
   CALL check(npairs == 0, "non-degenerate orbitals")

   ! one orbital of the second triplet: two pairs (5-6, 5-7), first pair reported
   CALL find_degenerate_splits(energies, [4, 5], eps, npairs, first_active, first_other)
   CALL check(npairs == 4 .AND. first_active == 4 .AND. first_other == 2, &
              "two split triplets (4 with 2, 3 and 5 with 6, 7)")

   ! order of the active list does not matter for the count
   CALL find_degenerate_splits(energies, [5, 4], eps, npairs, first_active, first_other)
   CALL check(npairs == 4 .AND. first_active == 5 .AND. first_other == 6, "reversed active list")

   ! active indices beyond the available orbitals (virtuals computed later) are skipped
   CALL find_degenerate_splits(energies, [8, 9, 10], eps, npairs, first_active, first_other)
   CALL check(npairs == 0, "indices beyond the available orbitals")

   ! splitting below the tolerance counts as degenerate, above it does not
   energies(6) = energies(5) + 0.5_dp*eps
   energies(7) = energies(5) + 2.0_dp*eps
   CALL find_degenerate_splits(energies, [5], eps, npairs, first_active, first_other)
   CALL check(npairs == 1 .AND. first_other == 6, "tolerance")

   WRITE (*, '(A)') "qs_active_space_utils_unittest: all tests passed"

CONTAINS

! **************************************************************************************************
!> \brief Stops with an error message if a test fails
!> \param condition test result
!> \param label test name
! **************************************************************************************************
   SUBROUTINE check(condition, label)
      LOGICAL, INTENT(IN)                                :: condition
      CHARACTER(LEN=*), INTENT(IN)                       :: label

      IF (.NOT. condition) THEN
         WRITE (*, '(A,A)') "FAILED: ", label
         ERROR STOP 1
      END IF
      WRITE (*, '(A,A)') "OK: ", label

   END SUBROUTINE check

END PROGRAM qs_active_space_utils_unittest
