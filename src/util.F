!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/util [1.0] *
!!
!!   NAME
!!     util
!!
!!   FUNCTION
!!     All kind of helpfull little routines
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
  PUBLIC :: sort, ran2, get_unit, gasdev, close_unit, &
       get_share, get_limit, include_list, getinv3x3,  &
       matmul_3x3, dotprod_3d, crossprod_3d, matvec_3x3, &
       transpose_3d,locate, find_boundary

  INTERFACE sort
     MODULE PROCEDURE sort2, sort2i
  END INTERFACE

  INTERFACE find_boundary
     MODULE PROCEDURE find_boundary1, find_boundary2, find_boundary3
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  FUNCTION locate(array,x) RESULT(x_index)

!   Purpose: Given an array array(1:n), and given a value x, a value x_index
!            is returned which is the index value of the array element equal
!            to the value x: x = array(x_index)
!            The array must be monotonic increasing.
!            x_index = 0 is returned, if no array element equal to the value
!            of x was found.

!   History: - Derived from the locate function described in
!              Numerical Recipes in Fortran 90 (09.01.2004,MK)

!   ***************************************************************************

    INTEGER, DIMENSION(:), INTENT(IN)        :: array
    INTEGER, INTENT(IN)                      :: x
    INTEGER                                  :: x_index

    INTEGER                                  :: jl, jm, ju, n

!   ---------------------------------------------------------------------------

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

!******************************************************************************

!!< from second edition of numerical recipes
!!    call with negative idum to initialize
!!>   then do not alter idum from successive calls

FUNCTION ran2 ( idum, icontrol )


    INTEGER, INTENT(INOUT)                   :: idum
    INTEGER, INTENT(IN), OPTIONAL            :: icontrol
    REAL(KIND=dp)                            :: ran2

    INTEGER                                  :: ia1, ia2, im1, im2, imm1, &
                                                iq1, iq2, ir1, ir2, ndiv, ntab
    REAL(KIND=dp)                            :: am, eps, rnmx

  PARAMETER (im1=2147483563,im2=2147483399,am=1.0_dp/im1,imm1=im1-1, &
       ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab=32, &
       ndiv=1+imm1/ntab,eps=1.2E-7_dp,rnmx=1.0_dp-eps)
  INTEGER idum2, j, k, iv(ntab), iy, iounit
  SAVE iv, iy, idum2
  DATA idum2/123456789/, iv/ntab*0/, iy/0/, iounit/10/

!------------------------------------------------------------------------------

!
!       if icontrol=0, generate a random number normally
!       if icontrol=1, load old random number table, dont get a new #
!       if icontrol=2, save current random number table and dont get a
!                      new random number
!

  IF ( .NOT. PRESENT ( icontrol)) THEN
     IF (idum<=0) THEN
        idum = MAX(-idum,1)
        idum2 = idum
        DO j = ntab + 8, 1, -1
           k = idum/iq1
           idum = ia1*(idum-k*iq1) - k*ir1
           IF (idum<0) idum = idum + im1
           IF (j<=ntab) iv(j) = idum
        END DO
        iy = iv(1)
     END IF
     k = idum/iq1
     idum = ia1*(idum-k*iq1) - k*ir1
     IF (idum<0) idum = idum + im1
     k = idum2/iq2
     idum2 = ia2*(idum2-k*iq2) - k*ir2
     IF (idum2<0) idum2 = idum2 + im2
     j = 1 + iy/ndiv
     iy = iv(j) - idum2
     iv(j) = idum
     IF (iy<1) iy = iy + imm1
     ran2 = MIN(am*iy,rnmx)
  ELSE IF (icontrol==1) THEN
     OPEN (unit=iounit,file='random.table')
     DO j = 1, ntab
        READ (iounit,*) iv(j)
     END DO
     READ (iounit,*) iy
     READ (iounit,*) idum
     READ (iounit,*) idum2
     CLOSE (iounit)
  ELSE IF (icontrol==2) THEN
     OPEN (unit=iounit,file='random.table')
     DO j = 1, ntab
        WRITE (iounit,*) iv(j)
     END DO
     WRITE (iounit,*) iy
     WRITE (iounit,*) idum
     WRITE (iounit,*) idum2
     CLOSE (iounit)
  END IF

END FUNCTION ran2

!******************************************************************************

FUNCTION gasdev(idum)


! Return value
    INTEGER, INTENT(INOUT)                   :: idum
    REAL(KIND=dp)                            :: gasdev

    INTEGER                                  :: iset
    REAL(KIND=dp)                            :: fac, gset, r, v1, v2

  DATA iset/0/
  SAVE iset, gset

!------------------------------------------------------------------------------

  IF (iset==0) THEN
1    CONTINUE
     v1 = 2.0_dp*ran2(idum) - 1.0_dp
     v2 = 2.0_dp*ran2(idum) - 1.0_dp
     r = v1*v1 + v2*v2
     IF (r>=1.0_dp.OR.r==0._dp) GO TO 1
     fac = SQRT(-2.0_dp*LOG(r)/r)
     gset = v1*fac
     gasdev = v2*fac
     iset = 1
  ELSE
     gasdev = gset
     iset = 0
  END IF

END FUNCTION gasdev

!******************************************************************************

FUNCTION get_unit() RESULT (unit)


! Return value
    INTEGER                                  :: unit

    LOGICAL                                  :: exists

!------------------------------------------------------------------------------
!..find a not associated unit

  unit = 10
  DO
     exists = .FALSE.
     INQUIRE (unit=unit,opened=exists)
     IF ( .NOT. exists) EXIT
     unit = unit + 1
  END DO

END FUNCTION get_unit

!******************************************************************************

SUBROUTINE close_unit ( u1, u2 )


    INTEGER, INTENT(IN)                      :: u1, u2

    INTEGER                                  :: iunit
    LOGICAL                                  :: exists

!------------------------------------------------------------------------------
!..close all open units from u1 to u2

  DO iunit = u1, u2
     exists = .FALSE.
     INQUIRE (unit=iunit,opened=exists)
     IF (exists) CLOSE (iunit)
  END DO

END SUBROUTINE close_unit

!******************************************************************************

!..divide m entries into n parts, return size of part me

FUNCTION get_share ( m, n, me ) RESULT ( nme )


! Return value
    INTEGER, INTENT(IN)                      :: m, n, me
    INTEGER                                  :: nme

    INTEGER                                  :: nl, nu
    REAL(KIND=dp)                            :: part

!------------------------------------------------------------------------------

  part = float(m)/float(n)
  nl = NINT(float(me)*part) + 1
  nu = NINT(float(me+1)*part)
  nl = MAX(1,nl)
  nu = MIN(m,nu)
  nme = nu - nl + 1
  nme = MAX(nme,0)

END FUNCTION get_share

!******************************************************************************

!..divide m entries into n parts, return size of part me

FUNCTION get_limit(m,n,me) RESULT (nlim)


! Return value
    INTEGER, INTENT(IN)                      :: m, n, me
    INTEGER                                  :: nlim( 2 )

    INTEGER                                  :: nl, nu
    REAL(KIND=dp)                            :: part

!------------------------------------------------------------------------------

  part = float(m)/float(n)
  nl = NINT(float(me)*part) + 1
  nu = NINT(float(me+1)*part)
  nlim(1) = MAX(1,nl)
  nlim(2) = MIN(m,nu)

END FUNCTION get_limit

!******************************************************************************

FUNCTION include_list(i,j) RESULT (doit)


! Return value
    INTEGER, INTENT(IN)                      :: i, j
    LOGICAL                                  :: doit

!------------------------------------------------------------------------------

  IF ( i >= j ) THEN
     doit = MOD(i+j,2) == 0
  ELSE IF (i < j) THEN
     doit = MOD ( i + j, 2 ) /= 0
  ELSE
     doit = .FALSE.
  END IF

END FUNCTION include_list

!******************************************************************************

SUBROUTINE getinv3x3 ( hmat, hmati )


    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: hmat
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT)                            :: hmati

    REAL(KIND=dp)                            :: det, odet

!------------------------------------------------------------------------------

  det = hmat(1,1)*(hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2)) + &
       hmat(1,2)*(hmat(2,3)*hmat(3,1)-hmat(2,1)*hmat(3,3)) + &
       hmat(1,3)*(hmat(2,1)*hmat(3,2)-hmat(2,2)*hmat(3,1))
  odet = 1.0_dp / det
  hmati(1,1) = (hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2))*odet
  hmati(2,2) = (hmat(1,1)*hmat(3,3)-hmat(1,3)*hmat(3,1))*odet
  hmati(3,3) = (hmat(1,1)*hmat(2,2)-hmat(1,2)*hmat(2,1))*odet
  hmati(1,2) = (hmat(1,3)*hmat(3,2)-hmat(1,2)*hmat(3,3))*odet
  hmati(2,1) = (hmat(3,1)*hmat(2,3)-hmat(2,1)*hmat(3,3))*odet
  hmati(1,3) = (hmat(1,2)*hmat(2,3)-hmat(1,3)*hmat(2,2))*odet
  hmati(3,1) = (hmat(2,1)*hmat(3,2)-hmat(3,1)*hmat(2,2))*odet
  hmati(2,3) = (hmat(1,3)*hmat(2,1)-hmat(2,3)*hmat(1,1))*odet
  hmati(3,2) = (hmat(3,1)*hmat(1,2)-hmat(3,2)*hmat(1,1))*odet

END SUBROUTINE getinv3x3

!******************************************************************************

FUNCTION matvec_3x3 ( mat, vec )


! Return value
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: mat
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: vec
    REAL(KIND=dp), DIMENSION(3)              :: matvec_3x3

!------------------------------------------------------------------------------

  matvec_3x3(1)=mat(1,1)*vec(1)+mat(1,2)*vec(2)+mat(1,3)*vec(3)
  matvec_3x3(2)=mat(2,1)*vec(1)+mat(2,2)*vec(2)+mat(2,3)*vec(3)
  matvec_3x3(3)=mat(3,1)*vec(1)+mat(3,2)*vec(2)+mat(3,3)*vec(3)

END FUNCTION matvec_3x3

!******************************************************************************

FUNCTION dotprod_3d ( vec1, vec2 )


! Return value
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: vec1, vec2
    REAL(KIND=dp)                            :: dotprod_3d

!------------------------------------------------------------------------------

  dotprod_3d = &
       vec1 ( 1 ) * vec2 ( 1 ) &
       + vec1 ( 2 ) * vec2 ( 2 ) &
       + vec1 ( 3 ) * vec2 ( 3 )

END FUNCTION dotprod_3d

!******************************************************************************

FUNCTION crossprod_3d ( vec1, vec2 )


! Return value
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: vec1, vec2
    REAL(KIND=dp), DIMENSION(3)              :: crossprod_3d

!------------------------------------------------------------------------------

  crossprod_3d ( 1 ) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
  crossprod_3d ( 2 ) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
  crossprod_3d ( 3 ) = vec1(1) * vec2(2) - vec1(2) * vec2(1)

END FUNCTION crossprod_3d

!******************************************************************************

FUNCTION transpose_3d ( mat )


! Return value
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: mat
    REAL(KIND=dp), DIMENSION(3, 3)           :: transpose_3d

    INTEGER                                  :: i, j

!------------------------------------------------------------------------------

  DO i = 1, 3
     DO j = 1, 3
        transpose_3d ( j, i ) = mat ( i, j )
     END DO
  END DO

END FUNCTION transpose_3d

!******************************************************************************

FUNCTION matmul_3x3 ( mat1, mat2 )


! Return value
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: mat1, mat2
    REAL(KIND=dp), DIMENSION(3, 3)           :: matmul_3x3

!------------------------------------------------------------------------------

  matmul_3x3(1,1)=mat1(1,1)*mat2(1,1)+mat1(1,2)*mat2(2,1)+mat1(1,3)*mat2(3,1)
  matmul_3x3(1,2)=mat1(1,1)*mat2(1,2)+mat1(1,2)*mat2(2,2)+mat1(1,3)*mat2(3,2)
  matmul_3x3(1,3)=mat1(1,1)*mat2(1,3)+mat1(1,2)*mat2(2,3)+mat1(1,3)*mat2(3,3)
  matmul_3x3(2,1)=mat1(2,1)*mat2(1,1)+mat1(2,2)*mat2(2,1)+mat1(2,3)*mat2(3,1)
  matmul_3x3(2,2)=mat1(2,1)*mat2(1,2)+mat1(2,2)*mat2(2,2)+mat1(2,3)*mat2(3,2)
  matmul_3x3(2,3)=mat1(2,1)*mat2(1,3)+mat1(2,2)*mat2(2,3)+mat1(2,3)*mat2(3,3)
  matmul_3x3(3,1)=mat1(3,1)*mat2(1,1)+mat1(3,2)*mat2(2,1)+mat1(3,3)*mat2(3,1)
  matmul_3x3(3,2)=mat1(3,1)*mat2(1,2)+mat1(3,2)*mat2(2,2)+mat1(3,3)*mat2(3,2)
  matmul_3x3(3,3)=mat1(3,1)*mat2(1,3)+mat1(3,2)*mat2(2,3)+mat1(3,3)*mat2(3,3)

END FUNCTION matmul_3x3

!******************************************************************************
SUBROUTINE find_boundary1(num_array,ntot,first,last,search)


    INTEGER, POINTER                         :: num_array(:)
    INTEGER                                  :: ntot, first, last, search

    INTEGER                                  :: i
    LOGICAL                                  :: found

!------------------------------------------------------------------------------

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
        EXIT                  !Exit this do loop???
      END IF
      found=.FALSE.
    END IF
  END DO

END SUBROUTINE find_boundary1

!******************************************************************************
SUBROUTINE find_boundary2(num_array1,num_array2,ntot,first,last,search1,search2)


    INTEGER, POINTER                         :: num_array1(:), num_array2(:)
    INTEGER                                  :: ntot, first, last, search1, &
                                                search2

    INTEGER                                  :: i, tfirst, tlast
    LOGICAL                                  :: found

!------------------------------------------------------------------------------

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
        EXIT                  !Exit this do loop???
      END IF
      found=.FALSE.
    END IF
  END DO

END SUBROUTINE find_boundary2

!******************************************************************************
SUBROUTINE find_boundary3(num_array1,num_array2,num_array3,ntot,first,last,&
                          search1,search2,search3)


    INTEGER, POINTER                         :: num_array1(:), num_array2(:), &
                                                num_array3(:)
    INTEGER                                  :: ntot, first, last, search1, &
                                                search2, search3

    INTEGER                                  :: i, tfirst, tlast
    LOGICAL                                  :: found

!------------------------------------------------------------------------------

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
        EXIT                  !Exit this do loop???
      END IF
      found=.FALSE.
    END IF
  END DO

END SUBROUTINE find_boundary3


END MODULE util

!******************************************************************************
!!*****
