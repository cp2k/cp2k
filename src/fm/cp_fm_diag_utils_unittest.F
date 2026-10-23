!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM cp_fm_diag_utils_unittest
   USE cp_fm_diag_utils,                ONLY: cp_fm_max_ncpu_non_zero_column
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE message_passing,                 ONLY: mp_para_env_type

   IMPLICIT NONE

   INTEGER, PARAMETER :: orders(2) = [509, 536]
   INTEGER :: i, itest, n
   INTEGER, EXTERNAL :: numroc
   TYPE(cp_fm_struct_type) :: matrix_struct
   TYPE(mp_para_env_type), TARGET :: para_env

   ! Layout metadata only: exercise large process grids even in two-rank CI.
   para_env%num_pe = 32
   matrix_struct%para_env => para_env
   matrix_struct%nrow_block = 64
   matrix_struct%ncol_block = 64
   ALLOCATE (matrix_struct%nrow_locals(0:3), matrix_struct%ncol_locals(0:7))

   DO itest = 1, SIZE(orders)
      n = orders(itest)
      matrix_struct%nrow_global = n
      matrix_struct%ncol_global = n
      DO i = 0, 3
         matrix_struct%nrow_locals(i) = numroc(n, 64, i, 0, 4)
      END DO
      DO i = 0, 7
         matrix_struct%ncol_locals(i) = numroc(n, 64, i, 0, 8)
      END DO

      IF (cp_fm_max_ncpu_non_zero_column(matrix_struct, 32, 64) /= 32) THEN
         ERROR STOP "Valid original redistribution layout changed"
      END IF
      ! 20 ranks use a 2x10 grid with empty columns; 18 ranks use a safe 3x6 grid.
      IF (cp_fm_max_ncpu_non_zero_column(matrix_struct, 20, 64) /= 18) THEN
         ERROR STOP "Unsafe forced redistribution layout accepted"
      END IF
      IF (cp_fm_max_ncpu_non_zero_column(matrix_struct, 12, 64) /= 12) THEN
         ERROR STOP "Safe redistribution heuristic target changed"
      END IF
      IF (cp_fm_max_ncpu_non_zero_column(matrix_struct, 8, 64) /= 8) THEN
         ERROR STOP "Safe smaller redistribution target changed"
      END IF
      IF (cp_fm_max_ncpu_non_zero_column(matrix_struct, 20, 32) /= 20) THEN
         ERROR STOP "Redistribution check ignored the final block size"
      END IF
   END DO

   ! Mandatory redistribution must still work when forced redistribution is off.
   matrix_struct%nrow_global = 1
   matrix_struct%ncol_global = 1
   matrix_struct%nrow_locals = 0
   matrix_struct%ncol_locals = 0
   matrix_struct%nrow_locals(0) = 1
   matrix_struct%ncol_locals(0) = 1
   IF (cp_fm_max_ncpu_non_zero_column(matrix_struct, 32, 64) /= 1) THEN
      ERROR STOP "Empty original redistribution layout was retained"
   END IF

   DEALLOCATE (matrix_struct%nrow_locals, matrix_struct%ncol_locals)

   ! Changing blocks also recreates the grid, even when the rank count is unchanged.
   para_env%num_pe = 20
   matrix_struct%nrow_global = 509
   matrix_struct%ncol_global = 509
   matrix_struct%nrow_block = 96
   matrix_struct%ncol_block = 96
   ALLOCATE (matrix_struct%nrow_locals(0:3), matrix_struct%ncol_locals(0:4))
   DO i = 0, 3
      matrix_struct%nrow_locals(i) = numroc(509, 96, i, 0, 4)
   END DO
   DO i = 0, 4
      matrix_struct%ncol_locals(i) = numroc(509, 96, i, 0, 5)
   END DO
   IF (cp_fm_max_ncpu_non_zero_column(matrix_struct, 20, 64) /= 18) THEN
      ERROR STOP "Block-size change skipped the target grid check"
   END IF
   DEALLOCATE (matrix_struct%nrow_locals, matrix_struct%ncol_locals)

   ! A safe default grid must not hide an invalid original custom grid.
   para_env%num_pe = 32
   matrix_struct%nrow_block = 64
   matrix_struct%ncol_block = 64
   ALLOCATE (matrix_struct%nrow_locals(0:1), matrix_struct%ncol_locals(0:15))
   DO i = 0, 1
      matrix_struct%nrow_locals(i) = numroc(509, 64, i, 0, 2)
   END DO
   DO i = 0, 15
      matrix_struct%ncol_locals(i) = numroc(509, 64, i, 0, 16)
   END DO
   IF (cp_fm_max_ncpu_non_zero_column(matrix_struct, 32, 64) /= 24) THEN
      ERROR STOP "Invalid original grid skipped mandatory redistribution"
   END IF
   DEALLOCATE (matrix_struct%nrow_locals, matrix_struct%ncol_locals)
END PROGRAM cp_fm_diag_utils_unittest
