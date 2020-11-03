!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM memory_utilities_TEST
   USE kinds,                           ONLY: dp
   USE memory_utilities,                ONLY: reallocate

   IMPLICIT NONE

   CALL check_real_rank1_allocated()
   CALL check_real_rank1_unallocated()

   CALL check_real_rank2_allocated()
   CALL check_real_rank2_unallocated()

CONTAINS
! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_real_rank1_allocated()
      INTEGER                                            :: idx
      REAL(KIND=dp), DIMENSION(:), POINTER               :: real_arr

      ALLOCATE (real_arr(10))
      real_arr = [(idx, idx=1, 10)]

      CALL reallocate(real_arr, 1, 20)

      IF (.NOT. ALL(real_arr(1:10) == [(idx, idx=1, 10)])) &
         ERROR STOP "check_real_rank1_allocated: reallocating changed the initial values"

      IF (.NOT. ALL(real_arr(11:20) == 0.)) &
         ERROR STOP "check_real_rank1_allocated: reallocation failed to initialise new values with 0."

      DEALLOCATE (real_arr)

      PRINT *, "check_real_rank1_allocated: OK"
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_real_rank1_unallocated()
      REAL(KIND=dp), DIMENSION(:), POINTER               :: real_arr => NULL()

      CALL reallocate(real_arr, 1, 20)

      IF (.NOT. ALL(real_arr(1:20) == 0.)) &
         ERROR STOP "check_real_rank1_unallocated: reallocation failed to initialise new values with 0."

      DEALLOCATE (real_arr)

      PRINT *, "check_real_rank1_unallocated: OK"
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_real_rank2_allocated()
      INTEGER                                            :: idx
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: real_arr

      ALLOCATE (real_arr(5, 2))
      real_arr = RESHAPE([(idx, idx=1, 10)], [5, 2])

      CALL reallocate(real_arr, 1, 10, 1, 5)

      IF (.NOT. (ALL(real_arr(1:5, 1) == [(idx, idx=1, 5)]) .AND. ALL(real_arr(1:5, 2) == [(idx, idx=6, 10)]))) &
         ERROR STOP "check_real_rank2_allocated: reallocating changed the initial values"

      IF (.NOT. (ALL(real_arr(6:10, 1:2) == 0.) .AND. ALL(real_arr(1:10, 3:5) == 0.))) &
         ERROR STOP "check_real_rank2_allocated: reallocation failed to initialise new values with 0."

      DEALLOCATE (real_arr)

      PRINT *, "check_real_rank1_allocated: OK"
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_real_rank2_unallocated()
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: real_arr => NULL()

      CALL reallocate(real_arr, 1, 10, 1, 5)

      IF (.NOT. ALL(real_arr(1:10, 1:5) == 0.)) &
         ERROR STOP "check_real_rank2_unallocated: reallocation failed to initialise new values with 0."

      DEALLOCATE (real_arr)

      PRINT *, "check_real_rank2_unallocated: OK"
   END SUBROUTINE
END PROGRAM
! vim: set ts=3 sw=3 tw=132 :
