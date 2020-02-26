!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE submatrix_methods

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: qsort, qsort_two

CONTAINS

! **************************************************************************************************
!> \brief sort integer array using quicksort
!> \param arr_a - input array
!> \param left - left boundary of region to be sorted
!> \param right - right boundary of region to be sorted
! **************************************************************************************************
   RECURSIVE PURE SUBROUTINE qsort(arr_a, left, right)

      INTEGER, DIMENSION(:), INTENT(inout)               :: arr_a
      INTEGER, INTENT(in)                                :: left, right

      INTEGER                                            :: i, j, pivot_a, tmp

      IF (right - left .LT. 1) RETURN

      i = left
      j = right - 1
      pivot_a = arr_a(right)

      DO
         DO WHILE (arr_a(i) .LT. pivot_a)
            i = i + 1
         ENDDO
         DO WHILE ((j .GT. i) .AND. (arr_a(j) .GE. pivot_a))
            j = j - 1
         ENDDO
         IF (i .LT. j) THEN
            tmp = arr_a(i)
            arr_a(i) = arr_a(j)
            arr_a(j) = tmp
         ELSE
            EXIT
         ENDIF
      ENDDO

      IF (arr_a(i) .GT. pivot_a) THEN
         tmp = arr_a(i)
         arr_a(i) = arr_a(right)
         arr_a(right) = tmp
      ENDIF

      IF (i - 1 .GT. left) CALL qsort(arr_a, left, i - 1)
      IF (i + 1 .LT. right) CALL qsort(arr_a, i + 1, right)

   END SUBROUTINE qsort

! **************************************************************************************************
!> \brief sort two integer arrays using quicksort, using the second list as second-level sorting criterion
!> \param arr_a - first input array
!> \param arr_b - second input array
!> \param left - left boundary of region to be sorted
!> \param right - right boundary of region to be sorted
! **************************************************************************************************
   RECURSIVE PURE SUBROUTINE qsort_two(arr_a, arr_b, left, right)

      INTEGER, DIMENSION(:), INTENT(inout)               :: arr_a, arr_b
      INTEGER, INTENT(in)                                :: left, right

      INTEGER                                            :: i, j, pivot_a, pivot_b, tmp

      IF (right - left .LT. 1) RETURN

      i = left
      j = right - 1
      pivot_a = arr_a(right)
      pivot_b = arr_b(right)

      DO
         DO WHILE ((arr_a(i) .LT. pivot_a) .OR. ((arr_a(i) .EQ. pivot_a) .AND. (arr_b(i) .LT. pivot_b)))
            i = i + 1
         ENDDO
         DO WHILE ((j .GT. i) .AND. ((arr_a(j) .GT. pivot_a) .OR. ((arr_a(j) .EQ. pivot_a) .AND. (arr_b(j) .GE. pivot_b))))
            j = j - 1
         ENDDO
         IF (i .LT. j) THEN
            tmp = arr_a(i)
            arr_a(i) = arr_a(j)
            arr_a(j) = tmp
            tmp = arr_b(i)
            arr_b(i) = arr_b(j)
            arr_b(j) = tmp
         ELSE
            EXIT
         ENDIF
      ENDDO

      IF ((arr_a(i) .GT. pivot_a) .OR. (arr_a(i) .EQ. pivot_a .AND. arr_b(i) .GT. pivot_b)) THEN
         tmp = arr_a(i)
         arr_a(i) = arr_a(right)
         arr_a(right) = tmp
         tmp = arr_b(i)
         arr_b(i) = arr_b(right)
         arr_b(right) = tmp
      ENDIF

      IF (i - 1 .GT. left) CALL qsort_two(arr_a, arr_b, left, i - 1)
      IF (i + 1 .LT. right) CALL qsort_two(arr_a, arr_b, i + 1, right)

   END SUBROUTINE qsort_two

END MODULE submatrix_methods
