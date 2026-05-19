!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM realspace_grid_cube_unittest
   USE kinds,                           ONLY: dp
   USE realspace_grid_cube,             ONLY: cube_read_values,&
                                              cube_value_format,&
                                              cube_values_format

   IMPLICIT NONE

   CHARACTER(LEN=13)                                  :: value
   CHARACTER(LEN=92)                                  :: slice_values
   CHARACTER(LEN=79)                                  :: values
   REAL(KIND=dp), DIMENSION(6)                        :: buffer, reference
   REAL(KIND=dp), DIMENSION(7)                        :: slice_buffer, slice_reference

   WRITE (value, cube_value_format) 0.33004E-101_dp
   IF (INDEX(value, "E-101") <= 0) &
      ERROR STOP "Cube value format must preserve three-digit negative exponents."

   reference = [0.0_dp, -0.91458E-35_dp, 0.33004E-101_dp, &
                -0.32487E-103_dp, 1.0_dp, -1.0_dp]
   WRITE (values(1:78), cube_values_format) reference
   values(79:79) = NEW_LINE("C")
   IF (INDEX(values, "E-101") <= 0) &
      ERROR STOP "Cube line format must preserve E-101."
   IF (INDEX(values, "E-103") <= 0) &
      ERROR STOP "Cube line format must preserve E-103."
   CALL cube_read_values(values, buffer)
   IF (MAXVAL(ABS(buffer - reference)) > 1.0E-12_dp) &
      ERROR STOP "Cube reader must parse adjacent explicit-exponent values."

   slice_reference(1:6) = reference
   slice_reference(7) = 0.27183E-123_dp
   WRITE (slice_values(1:78), cube_values_format) slice_reference(1:6)
   slice_values(79:79) = NEW_LINE("C")
   WRITE (slice_values(80:92), cube_value_format) slice_reference(7)
   CALL cube_read_values(slice_values, slice_buffer)
   IF (MAXVAL(ABS(slice_buffer - slice_reference)) > 1.0E-12_dp) &
      ERROR STOP "Cube reader must parse multi-line z-slices."

END PROGRAM realspace_grid_cube_unittest
