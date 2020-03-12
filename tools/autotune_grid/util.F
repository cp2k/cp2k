!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/util [1.0] *
!!
!!   NAME
!!     util
!!
!!   FUNCTION
!!     All kind of helpful little routines
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE util

  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sort
  
  INTERFACE sort
     MODULE PROCEDURE sort2, sort2i
  END INTERFACE

CONTAINS

SUBROUTINE sort2 ( arr, n, index )


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(INOUT)             :: arr(1:n)
    INTEGER, INTENT(OUT)                     :: INDEX(1:n)

    INTEGER, PARAMETER                       :: m = 7, nstack = 50

    INTEGER                                  :: i, ib, ir, istack(1:nstack), &
                                                itemp, j, jstack, k, l
    REAL(KIND=dp)                            :: a, temp

!------------------------------------------------------------------------------

  DO i = 1, n
     INDEX(i) = i
  END DO
  !
  !Temporary fix for INTEL compiler...
  IF (n==0) RETURN
  !Temporary fix for INTEL compiler...
  !
  IF (ALL(arr==arr(1))) RETURN ! Nothing to order..
  jstack = 0
  l = 1
  ir = n
1 IF (ir-l<m) THEN
     DO j = l + 1, ir
        a = arr(j)
        ib = INDEX(j)
        DO i = j - 1, 1, -1
           IF (arr(i)<=a) GO TO 2
           arr(i+1) = arr(i)
           INDEX(i+1) = INDEX(i)
        END DO
        i = 0
2       arr(i+1) = a
        INDEX(i+1) = ib
     END DO
     IF (jstack==0) RETURN
     ir = istack(jstack)
     l = istack(jstack-1)
     jstack = jstack - 2
  ELSE
     k = (l+ir)/2
     temp = arr(k)
     arr(k) = arr(l+1)
     arr(l+1) = temp
     itemp = INDEX(k)
     INDEX(k) = INDEX(l+1)
     INDEX(l+1) = itemp
     IF (arr(l+1)>arr(ir)) THEN
        temp = arr(l+1)
        arr(l+1) = arr(ir)
        arr(ir) = temp
        itemp = INDEX(l+1)
        INDEX(l+1) = INDEX(ir)
        INDEX(ir) = itemp
     END IF
     IF (arr(l)>arr(ir)) THEN
        temp = arr(l)
        arr(l) = arr(ir)
        arr(ir) = temp
        itemp = INDEX(l)
        INDEX(l) = INDEX(ir)
        INDEX(ir) = itemp
     END IF
     IF (arr(l+1)>arr(l)) THEN
        temp = arr(l+1)
        arr(l+1) = arr(l)
        arr(l) = temp
        itemp = INDEX(l+1)
        INDEX(l+1) = INDEX(l)
        INDEX(l) = itemp
     END IF
     i = l + 1
     j = ir
     a = arr(l)
     ib = INDEX(l)
3    CONTINUE
     i = i + 1
     IF (arr(i)<a) GO TO 3
4    CONTINUE
     j = j - 1
     IF (arr(j)>a) GO TO 4
     IF (j<i) GO TO 5
     temp = arr(i)
     arr(i) = arr(j)
     arr(j) = temp
     itemp = INDEX(i)
     INDEX(i) = INDEX(j)
     INDEX(j) = itemp
     GO TO 3
5    arr(l) = arr(j)
     arr(j) = a
     INDEX(l) = INDEX(j)
     INDEX(j) = ib
     jstack = jstack + 2
     IF (jstack>nstack) STOP ' Nstack too small in sortr'
     IF (ir-i+1>=j-l) THEN
        istack(jstack) = ir
        istack(jstack-1) = i
        ir = j - 1
     ELSE
        istack(jstack) = j - 1
        istack(jstack-1) = l
        l = i
     END IF
  END IF

  GO TO 1

END SUBROUTINE sort2

!******************************************************************************

SUBROUTINE sort2i ( iarr, n, index )


    INTEGER, INTENT(IN)                      :: n
    INTEGER, INTENT(INOUT)                   :: iarr(1:n)
    INTEGER, INTENT(OUT)                     :: INDEX(1:n)

    INTEGER, PARAMETER                       :: m = 7, nstack = 50

    INTEGER                                  :: a, i, ib, ir, &
                                                istack(1:nstack), itemp, j, &
                                                jstack, k, l, temp

!------------------------------------------------------------------------------

  DO i = 1, n
     INDEX(i) = i
  END DO
  !
  !Temporary fix for INTEL compiler...
  IF (n==0) RETURN
  !Temporary fix for INTEL compiler...
  !
  IF (ALL(iarr==iarr(1))) RETURN ! Nothing to order..
  jstack = 0
  l = 1
  ir = n
1 IF (ir-l<m) THEN
     DO j = l + 1, ir
        a = iarr(j)
        ib = INDEX(j)
        DO i = j - 1, 1, -1
           IF (iarr(i)<=a) GO TO 2
           iarr(i+1) = iarr(i)
           INDEX(i+1) = INDEX(i)
        END DO
        i = 0
2       iarr(i+1) = a
        INDEX(i+1) = ib
     END DO
     IF (jstack==0) RETURN
     ir = istack(jstack)
     l = istack(jstack-1)
     jstack = jstack - 2
  ELSE
     k = (l+ir)/2
     temp = iarr(k)
     iarr(k) = iarr(l+1)
     iarr(l+1) = temp
     itemp = INDEX(k)
     INDEX(k) = INDEX(l+1)
     INDEX(l+1) = itemp
     IF (iarr(l+1)>iarr(ir)) THEN
        temp = iarr(l+1)
        iarr(l+1) = iarr(ir)
        iarr(ir) = temp
        itemp = INDEX(l+1)
        INDEX(l+1) = INDEX(ir)
        INDEX(ir) = itemp
     END IF
     IF (iarr(l)>iarr(ir)) THEN
        temp = iarr(l)
        iarr(l) = iarr(ir)
        iarr(ir) = temp
        itemp = INDEX(l)
        INDEX(l) = INDEX(ir)
        INDEX(ir) = itemp
     END IF
     IF (iarr(l+1)>iarr(l)) THEN
        temp = iarr(l+1)
        iarr(l+1) = iarr(l)
        iarr(l) = temp
        itemp = INDEX(l+1)
        INDEX(l+1) = INDEX(l)
        INDEX(l) = itemp
     END IF
     i = l + 1
     j = ir
     a = iarr(l)
     ib = INDEX(l)
3    CONTINUE
     i = i + 1
     IF (iarr(i)<a) GO TO 3
4    CONTINUE
     j = j - 1
     IF (iarr(j)>a) GO TO 4
     IF (j<i) GO TO 5
     temp = iarr(i)
     iarr(i) = iarr(j)
     iarr(j) = temp
     itemp = INDEX(i)
     INDEX(i) = INDEX(j)
     INDEX(j) = itemp
     GO TO 3
5    iarr(l) = iarr(j)
     iarr(j) = a
     INDEX(l) = INDEX(j)
     INDEX(j) = ib
     jstack = jstack + 2
     IF (jstack>nstack) STOP ' Nstack too small in sortr'
     IF (ir-i+1>=j-l) THEN
        istack(jstack) = ir
        istack(jstack-1) = i
        ir = j - 1
     ELSE
        istack(jstack) = j - 1
        istack(jstack-1) = l
        l = i
     END IF
  END IF

  GO TO 1

END SUBROUTINE sort2i

END MODULE

