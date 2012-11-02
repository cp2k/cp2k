! --------------------------------------------------------------------------------------------------
!
! This file contains the compute intensive kernels for the Householder transformations.
! It should be compiled with the highest possible optimization level.
!
! On Intel use -O3 -xSSE4.2 (or the SSE level fitting to your CPU)
! 
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".
!
! --------------------------------------------------------------------------------------------------

SUBROUTINE single_hh_trafo_complex(q, hh, nb, nq, ldq)

   IMPLICIT NONE

   INTEGER, INTENT(in) :: nb, nq, ldq
   COMPLEX*16, INTENT(inout) :: q(ldq,*)
   COMPLEX*16, INTENT(in) :: hh(*)

   INTEGER i

   ! Safety only:

   IF(MOD(ldq,4) /= 0) STOP 'double_hh_trafo: ldq not divisible by 4!'

   ! Do the Householder transformations

   ! Always a multiple of 4 Q-rows is transformed, even if nq is smaller

   DO i=1,nq-8,12
      CALL hh_trafo_complex_kernel_12(q(i,1),hh, nb, ldq)
   ENDDO

   ! i > nq-8 now, i.e. at most 8 rows remain

   IF(nq-i+1 > 4) THEN
      CALL hh_trafo_complex_kernel_8(q(i,1),hh, nb, ldq)
   ELSE IF(nq-i+1 > 0) THEN
      CALL hh_trafo_complex_kernel_4(q(i,1),hh, nb, ldq)
   ENDIF

END

! --------------------------------------------------------------------------------------------------

SUBROUTINE hh_trafo_complex_kernel_12(q, hh, nb, ldq)

   IMPLICIT NONE

   INTEGER, INTENT(in) :: nb, ldq
   COMPLEX*16, INTENT(inout) :: q(ldq,*)
   COMPLEX*16, INTENT(in) :: hh(*)

   COMPLEX*16 x1, x2, x3, x4, x5, x6, x7, x8, x9, xa, xb, xc
   COMPLEX*16 h1, tau1
   INTEGER i


   x1 = q(1,1)
   x2 = q(2,1)
   x3 = q(3,1)
   x4 = q(4,1)
   x5 = q(5,1)
   x6 = q(6,1)
   x7 = q(7,1)
   x8 = q(8,1)
   x9 = q(9,1)
   xa = q(10,1)
   xb = q(11,1)
   xc = q(12,1)

!DEC$ VECTOR ALIGNED
   DO i=2,nb
      h1 = CONJG(hh(i))
      x1 = x1 + q(1,i)*h1
      x2 = x2 + q(2,i)*h1
      x3 = x3 + q(3,i)*h1
      x4 = x4 + q(4,i)*h1
      x5 = x5 + q(5,i)*h1
      x6 = x6 + q(6,i)*h1
      x7 = x7 + q(7,i)*h1
      x8 = x8 + q(8,i)*h1
      x9 = x9 + q(9,i)*h1
      xa = xa + q(10,i)*h1
      xb = xb + q(11,i)*h1
      xc = xc + q(12,i)*h1
   ENDDO

   tau1 = hh(1)

   h1 = -tau1
   x1 = x1*h1
   x2 = x2*h1
   x3 = x3*h1
   x4 = x4*h1
   x5 = x5*h1
   x6 = x6*h1
   x7 = x7*h1
   x8 = x8*h1
   x9 = x9*h1
   xa = xa*h1
   xb = xb*h1
   xc = xc*h1

   q(1,1) = q(1,1) + x1
   q(2,1) = q(2,1) + x2
   q(3,1) = q(3,1) + x3
   q(4,1) = q(4,1) + x4
   q(5,1) = q(5,1) + x5
   q(6,1) = q(6,1) + x6
   q(7,1) = q(7,1) + x7
   q(8,1) = q(8,1) + x8
   q(9,1) = q(9,1) + x9
   q(10,1) = q(10,1) + xa
   q(11,1) = q(11,1) + xb
   q(12,1) = q(12,1) + xc

!DEC$ VECTOR ALIGNED
   DO i=2,nb
      h1 = hh(i)
      q(1,i) = q(1,i) + x1*h1
      q(2,i) = q(2,i) + x2*h1
      q(3,i) = q(3,i) + x3*h1
      q(4,i) = q(4,i) + x4*h1
      q(5,i) = q(5,i) + x5*h1
      q(6,i) = q(6,i) + x6*h1
      q(7,i) = q(7,i) + x7*h1
      q(8,i) = q(8,i) + x8*h1
      q(9,i) = q(9,i) + x9*h1
      q(10,i) = q(10,i) + xa*h1
      q(11,i) = q(11,i) + xb*h1
      q(12,i) = q(12,i) + xc*h1
   ENDDO

END

! --------------------------------------------------------------------------------------------------

SUBROUTINE hh_trafo_complex_kernel_8(q, hh, nb, ldq)

   IMPLICIT NONE

   INTEGER, INTENT(in) :: nb, ldq
   COMPLEX*16, INTENT(inout) :: q(ldq,*)
   COMPLEX*16, INTENT(in) :: hh(*)

   COMPLEX*16 x1, x2, x3, x4, x5, x6, x7, x8
   COMPLEX*16 h1, tau1
   INTEGER i


   x1 = q(1,1)
   x2 = q(2,1)
   x3 = q(3,1)
   x4 = q(4,1)
   x5 = q(5,1)
   x6 = q(6,1)
   x7 = q(7,1)
   x8 = q(8,1)

!DEC$ VECTOR ALIGNED
   DO i=2,nb
      h1 = CONJG(hh(i))
      x1 = x1 + q(1,i)*h1
      x2 = x2 + q(2,i)*h1
      x3 = x3 + q(3,i)*h1
      x4 = x4 + q(4,i)*h1
      x5 = x5 + q(5,i)*h1
      x6 = x6 + q(6,i)*h1
      x7 = x7 + q(7,i)*h1
      x8 = x8 + q(8,i)*h1
   ENDDO

   tau1 = hh(1)

   h1 = -tau1
   x1 = x1*h1
   x2 = x2*h1
   x3 = x3*h1
   x4 = x4*h1
   x5 = x5*h1
   x6 = x6*h1
   x7 = x7*h1
   x8 = x8*h1

   q(1,1) = q(1,1) + x1
   q(2,1) = q(2,1) + x2
   q(3,1) = q(3,1) + x3
   q(4,1) = q(4,1) + x4
   q(5,1) = q(5,1) + x5
   q(6,1) = q(6,1) + x6
   q(7,1) = q(7,1) + x7
   q(8,1) = q(8,1) + x8

!DEC$ VECTOR ALIGNED
   DO i=2,nb
      h1 = hh(i)
      q(1,i) = q(1,i) + x1*h1
      q(2,i) = q(2,i) + x2*h1
      q(3,i) = q(3,i) + x3*h1
      q(4,i) = q(4,i) + x4*h1
      q(5,i) = q(5,i) + x5*h1
      q(6,i) = q(6,i) + x6*h1
      q(7,i) = q(7,i) + x7*h1
      q(8,i) = q(8,i) + x8*h1
   ENDDO

END

! --------------------------------------------------------------------------------------------------

SUBROUTINE hh_trafo_complex_kernel_4(q, hh, nb, ldq)

   IMPLICIT NONE

   INTEGER, INTENT(in) :: nb, ldq
   COMPLEX*16, INTENT(inout) :: q(ldq,*)
   COMPLEX*16, INTENT(in) :: hh(*)

   COMPLEX*16 x1, x2, x3, x4
   COMPLEX*16 h1, tau1
   INTEGER i


   x1 = q(1,1)
   x2 = q(2,1)
   x3 = q(3,1)
   x4 = q(4,1)

!DEC$ VECTOR ALIGNED
   DO i=2,nb
      h1 = CONJG(hh(i))
      x1 = x1 + q(1,i)*h1
      x2 = x2 + q(2,i)*h1
      x3 = x3 + q(3,i)*h1
      x4 = x4 + q(4,i)*h1
   ENDDO

   tau1 = hh(1)

   h1 = -tau1
   x1 = x1*h1
   x2 = x2*h1
   x3 = x3*h1
   x4 = x4*h1

   q(1,1) = q(1,1) + x1
   q(2,1) = q(2,1) + x2
   q(3,1) = q(3,1) + x3
   q(4,1) = q(4,1) + x4

!DEC$ VECTOR ALIGNED
   DO i=2,nb
      h1 = hh(i)
      q(1,i) = q(1,i) + x1*h1
      q(2,i) = q(2,i) + x2*h1
      q(3,i) = q(3,i) + x3*h1
      q(4,i) = q(4,i) + x4*h1
   ENDDO

END

! --------------------------------------------------------------------------------------------------
