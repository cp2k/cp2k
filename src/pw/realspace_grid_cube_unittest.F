!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM realspace_grid_cube_unittest
   USE kinds,                           ONLY: dp
   USE realspace_grid_cube,             ONLY: cube_value_format,&
                                              cube_values_format

   IMPLICIT NONE

   CHARACTER(LEN=13)                                  :: value
   CHARACTER(LEN=78)                                  :: values
   REAL(KIND=dp), DIMENSION(6)                        :: buffer

   WRITE (value, cube_value_format) 0.33004E-101_dp
   IF (INDEX(value, "E-101") <= 0) &
      ERROR STOP "Cube value format must preserve three-digit negative exponents."

   buffer = [0.0_dp, -0.91458E-35_dp, 0.33004E-101_dp, &
             -0.32487E-103_dp, 1.0_dp, -1.0_dp]
   WRITE (values, cube_values_format) buffer
   IF (INDEX(values, "E-101") <= 0) &
      ERROR STOP "Cube line format must preserve E-101."
   IF (INDEX(values, "E-103") <= 0) &
      ERROR STOP "Cube line format must preserve E-103."

END PROGRAM realspace_grid_cube_unittest
