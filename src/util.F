!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief All kind of helpfull little routines
!> \par History
!>      none
!> \author CJM & JGH
! *****************************************************************************
MODULE util
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_4,&
                                             int_8

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'util'
  PUBLIC :: sort,&
       get_limit,&
       locate,&
       find_boundary,&
       sort_unique

  INTERFACE sort
     MODULE PROCEDURE sort2, sort2i_int_4, sort2i_int_8
  END INTERFACE

  INTERFACE sort_unique
     MODULE PROCEDURE sort_unique1, sort_unique2
  END INTERFACE

  INTERFACE find_boundary
     MODULE PROCEDURE find_boundary1, find_boundary2, find_boundary3
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief Purpose: Given an array array(1:n), and given a value x, a value x_index
!>             is returned which is the index value of the array element equal
!>             to the value x: x = array(x_index)
!>             The array must be monotonic increasing.
!>             x_index = 0 is returned, if no array element equal to the value
!>             of x was found.
!> \par History
!>      Derived from the locate function described in
!>      Numerical Recipes in Fortran 90 (09.01.2004,MK)
! *****************************************************************************
  FUNCTION locate(array,x) RESULT(x_index)
    INTEGER, DIMENSION(:), INTENT(IN)        :: array
    INTEGER, INTENT(IN)                      :: x
    INTEGER                                  :: x_index

    CHARACTER(LEN=*), PARAMETER :: routineN = 'locate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: jl, jm, ju, n

    x_index = 0

    IF (x < array(1)) RETURN
    n = SIZE(array)
    IF (x > array(n)) RETURN
    jl = 0
    ju = n + 1
    DO WHILE (ju - jl > 1)
       jm = (ju + jl)/2
       IF (x >= array(jm)) THEN
          jl = jm
       ELSE
          ju = jm
       END IF
    END DO
    IF (x == array(jl)) x_index = jl
  END FUNCTION locate

! *****************************************************************************
!> \brief Sorts and returns a logical that checks if all elements are unique
!> \par History
!>      Teodoro Laino - Zurich University [tlaino] 04.2007
! *****************************************************************************
  SUBROUTINE sort_unique1( arr, unique)
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: arr
    LOGICAL, INTENT(OUT)                     :: unique

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sort_unique1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: wrk

    n = SIZE(arr)
    unique = .TRUE.
    ALLOCATE(wrk(n))
    CALL sort( arr, n, wrk)
    DO i = 2, n
       IF (arr(i)==arr(i-1)) THEN
          unique = .FALSE.
          EXIT
       END IF
    END DO
    DEALLOCATE(wrk)
  END SUBROUTINE sort_unique1

! *****************************************************************************
!> \brief Sorts and returns a logical that checks if all elements are unique
!> \par History
!>      Teodoro Laino - Zurich University [tlaino] 04.2007
! *****************************************************************************
  SUBROUTINE sort_unique2( arr, unique)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: arr
    LOGICAL, INTENT(OUT)                     :: unique

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sort_unique2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: wrk

    n = SIZE(arr)
    unique = .TRUE.
    ALLOCATE(wrk(n))
    CALL sort( arr, n, wrk)
    DO i = 2, n
       IF (arr(i)==arr(i-1)) THEN
          unique = .FALSE.
          EXIT
       END IF
    END DO
    DEALLOCATE(wrk)
  END SUBROUTINE sort_unique2

! *****************************************************************************
!> \brief Sorts an array of reals
! *****************************************************************************
  SUBROUTINE sort2 ( arr, n, index )
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(INOUT)             :: arr(1:n)
    INTEGER, INTENT(OUT)                     :: INDEX(1:n)

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sort2', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: m = 7, nstack = 50

    INTEGER                                  :: i, ib, ir, istack(1:nstack), &
                                                itemp, j, jstack, k, l
    REAL(KIND=dp)                            :: a, temp
#include "sort.f90"
  END SUBROUTINE sort2

! *****************************************************************************
!> \brief Sorts an array of integers
! *****************************************************************************
  SUBROUTINE sort2i_int_4 ( arr, n, index )
    INTEGER, INTENT(IN)                      :: n
    INTEGER(kind=int_4), INTENT(INOUT)                   :: arr(1:n)
    INTEGER, INTENT(OUT)                     :: INDEX(1:n)

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sort2i_int_4', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: m = 7, nstack = 50

    INTEGER(kind=int_4) :: a, temp
    INTEGER                                  :: i, ib, ir, &
                                                istack(1:nstack), itemp, j, &
                                                jstack, k, l
#include "sort.f90"
  END SUBROUTINE sort2i_int_4

! *****************************************************************************
!> \brief Sorts an array of int_8 integers
! *****************************************************************************
  SUBROUTINE sort2i_int_8 ( arr, n, index )
    INTEGER, INTENT(IN)                      :: n
    INTEGER(kind=int_8), INTENT(INOUT)       :: arr(1:n)
    INTEGER, INTENT(OUT)                     :: INDEX(1:n)

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sort2i_int_8', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: m = 7, nstack = 50

    INTEGER(kind=int_8)                      :: a,temp
    INTEGER                                  :: i, ib, ir, &
                                                istack(1:nstack), itemp, j, &
                                                jstack, k, l
#include "sort.f90"
  END SUBROUTINE sort2i_int_8

! *****************************************************************************
!> \brief divide m entries into n parts, return size of part me
! *****************************************************************************
  FUNCTION get_limit(m,n,me) RESULT (nlim)
    INTEGER, INTENT(IN)                      :: m, n, me
    INTEGER                                  :: nlim( 2 )

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_limit', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nl, nu
    REAL(KIND=dp)                            :: part

    part = float(m)/float(n)
    nl = NINT(float(me)*part) + 1
    nu = NINT(float(me+1)*part)
    nlim(1) = MAX(1,nl)
    nlim(2) = MIN(m,nu)

  END FUNCTION get_limit

! *****************************************************************************
!> \brief finds boundary where element search starts and ends in a 1D array
!>      array1:      XXXXXAAAAAAAAAXXDGFSFGWDDDDDDDAAAWE
!>                        |       |
!>                     start     end  (searching for A)
! *****************************************************************************
  SUBROUTINE find_boundary1(num_array,ntot,first,last,search)
    INTEGER, POINTER                         :: num_array(:)
    INTEGER                                  :: ntot, first, last, search

    CHARACTER(LEN=*), PARAMETER :: routineN = 'find_boundary1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: found

    found=.FALSE.
    first=0
    last=ntot

    DO i=1,ntot
       IF(num_array(i)==search) THEN
          IF(.NOT.found) THEN
             first = i
          END IF
          found=.TRUE.
       ELSE
          IF(found) THEN
             last = i-1
             EXIT
          END IF
          found=.FALSE.
       END IF
    END DO

  END SUBROUTINE find_boundary1

! *****************************************************************************
!> \brief finds boundary where element search1 starts and ends in array1 checking
!>      at the same time search2 in array2
!>      array1:      XXXXXAAAAAAAAAXXDGFSFGWDDDDDDDAAAWE
!>      array2:      XXXXASDEYYYYASDEFAAAARGASGASRGAWRRR
!>                           |  |
!>                       start  end  (searching for A and Y)
! *****************************************************************************
  SUBROUTINE find_boundary2(num_array1,num_array2,ntot,first,last,search1,search2)
    INTEGER, POINTER                         :: num_array1(:), num_array2(:)
    INTEGER                                  :: ntot, first, last, search1, &
                                                search2

    CHARACTER(LEN=*), PARAMETER :: routineN = 'find_boundary2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, tfirst, tlast
    LOGICAL                                  :: found

    found=.FALSE.
    first=0
    last=ntot

    CALL find_boundary(num_array1,ntot,tfirst,tlast,search1)
    last=tlast
    DO i=tfirst,tlast
       IF(num_array2(i)==search2) THEN
          IF(.NOT.found) THEN
             first = i
          END IF
          found=.TRUE.
       ELSE
          IF(found) THEN
             last = i-1
             EXIT
          END IF
          found=.FALSE.
       END IF
    END DO

  END SUBROUTINE find_boundary2

! *****************************************************************************
!> \brief finds boundary where element search1 starts and ends in array1 checking
!>      at the same time search2 in array2 and search3 in array3
!>      array1:      XXXXXAAAAAAAAAXXDGFSFGWDDDDDDDAAAWE
!>      array2:      XXXXASDEYYYYASDEFAAAARGASGASRGAWRRR
!>      array3:      ADFAWRESCXXSQWQWFQGEQGEQRQEHQEHQEHEQ
!>                            ||
!>                       start  end  (searching for A and Y and X)
! *****************************************************************************
  SUBROUTINE find_boundary3(num_array1,num_array2,num_array3,ntot,first,last,&
       search1,search2,search3)
    INTEGER, POINTER                         :: num_array1(:), num_array2(:), &
                                                num_array3(:)
    INTEGER                                  :: ntot, first, last, search1, &
                                                search2, search3

    CHARACTER(LEN=*), PARAMETER :: routineN = 'find_boundary3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, tfirst, tlast
    LOGICAL                                  :: found

    found=.FALSE.
    first=0
    last=ntot
    CALL find_boundary2(num_array1,num_array2,ntot,tfirst,tlast,search1,search2)
    last=tlast
    DO i=tfirst,tlast
       IF(num_array3(i)==search3) THEN
          IF(.NOT.found) THEN
             first = i
          END IF
          found=.TRUE.
       ELSE
          IF(found) THEN
             last = i-1
             EXIT
          END IF
          found=.FALSE.
       END IF
    END DO
  END SUBROUTINE find_boundary3

END MODULE util
