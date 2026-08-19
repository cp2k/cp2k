!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_trajectory_types
   USE trajana_kinds,                   ONLY: dp

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC :: cell_type
      REAL(dp) :: h(3, 3) = 0.0_dp
      LOGICAL :: periodic(3) = .TRUE.
      LOGICAL :: valid = .FALSE.
   END TYPE cell_type

   TYPE, PUBLIC :: frame_type
      INTEGER :: number = 0
      INTEGER :: natoms = 0
      CHARACTER(LEN=:), ALLOCATABLE :: comment
      CHARACTER(LEN=32), ALLOCATABLE :: label(:)
      REAL(dp), ALLOCATABLE :: value(:, :)
      TYPE(cell_type) :: cell = cell_type()
   END TYPE frame_type

END MODULE trajana_trajectory_types
