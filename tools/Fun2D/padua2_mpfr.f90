!
!  These routines have been converted to use mpfr variables 
!  from an 'older' version of Padua2D.
!
!  In this way, arbitrary precision can be used to construct the 
!  coefficients and evaluate the interpolation
!
!  Joost VandeVondele, November 2008
!
!     ******************************************************************
!
!                                                           From Padua2D
!
!                                       Marco Caliari and Marco Vianello
!                             Department of Pure and Applied Mathematics
!                                            University of Padua (Italy)
!                                        {mcaliari,marcov}@math.unipd.it
!
!                                                      Stefano De Marchi
!                                         Department of Computer Science
!                                           University of Verona (Italy)
!                                              stefano.demarchi@univr.it
!
!                                                        June 19th, 2007
!
MODULE padua2_mpfr
      USE mpfr
      USE mpfr_ops
      IMPLICIT NONE
CONTAINS

!     ******************************************************************
!
!     PARAMS type                i=input, o=output, a=auxiliary
!
!     DEG    integer             i maximum degree
!     PT     mpfr                i point
!     TCHEB  mpfr(0:DEG)         o Normalized Chebyshev polynomials at
!                                  the point PT
!     S1     mpfr                i scale factor for TCHEB
!
!     **********************************************************************
      SUBROUTINE cheb(deg,pt,tcheb,s1)
         INTEGER :: deg
         TYPE(mpfr_type) :: pt,s1,tcheb(0:deg)

         TYPE(mpfr_type) :: sqrt2
         INTEGER         :: j

         sqrt2="2"
         sqrt2=sqrt(sqrt2)

         IF (deg .GE. 0) THEN
               tcheb(0) = s1
         END IF      
         IF (deg .GE. 1) THEN
            tcheb(1) = sqrt2 * pt * s1
         END IF
         IF (deg .GE. 2) THEN
            tcheb(2) = 2 * pt * tcheb(1) - sqrt2 * tcheb(0)
         END IF
         DO J = 3,DEG
            tcheb(j) = 2 * pt * tcheb(j - 1) - tcheb(j - 2)
         ENDDO

      END SUBROUTINE cheb

!     ******************************************************************
!
!     PARAMS type                i=input, o=output, a=auxiliary           
!     
!     DEG    integer             i degree of approximation
!     PD1    mpfr(*)             o array of the first coordinates of 
!                                  Padua points
!     PD2    mpfr(*)             o array of the second coordinates of 
!                                  Padua points
!     WPD    mpfr(*)             o array of the weights
!     NPD    integer             o number of Padua points = 
!                                  (DEG + 1) * (DEG + 2) / 2
!
!     ******************************************************************
      SUBROUTINE pdpts(deg,pd1,pd2,wpd,npd)
         INTEGER ::  deg,npd
         TYPE(mpfr_type) :: pd1(*),pd2(*),wpd(*)

         INTEGER j,k,deg1,itemp0,kold,retval
         TYPE(mpfr_type) :: pi,rtemp0

         pi=1.0D0
         pi=atan(pi)*4

         IF (deg .EQ. 0) THEN
            pd1(1) = -1
            pd2(1) = -1
            WPD(1) =  2
            npd = 1
            RETURN
         END IF   

         npd = 0
         deg1 = deg + 1
         itemp0 = deg * deg1
         rtemp0 = pi / itemp0
         k = 0
         npd = npd + 1

         !        
         ! The points are generated starting from bottom-left
         ! and following the generating curve                           
         !        
         pd1(npd) = -1
         pd2(npd) = -1
         wpd(npd) = .CONVERT."1" / (2*itemp0)
         DO k = 1,deg - 1
            npd = npd + 1
            pd1(npd) = -COS(deg1 * k * rtemp0)
            pd2(npd) = -COS(deg * k * rtemp0)
            wpd(npd) = .CONVERT."2" / itemp0
         ENDDO

         k = deg
         npd = npd + 1
         pd1(npd) = -COS(deg1 * k * rtemp0)
         pd2(npd) = -COS(deg * k * rtemp0)
         wpd(npd) = .CONVERT."1" / itemp0
         kold = deg1
         DO j = 0,deg - 2
            k = kold
            npd = npd + 1
            pd1(npd) = -COS(deg1 * k * rtemp0)
            pd2(npd) = -COS(deg * k * rtemp0)
            wpd(npd) = .CONVERT."1" / itemp0
            DO k = kold + 1,kold + deg - j - 2
               npd = npd + 1
               pd1(npd) = -COS(deg1 * k * rtemp0)
               pd2(npd) = -COS(deg * k * rtemp0)
               wpd(npd) = .CONVERT."2" / itemp0
            ENDDO
            k = kold + deg - j - 1
            npd = npd + 1
            pd1(npd) = -COS(deg1 * k * rtemp0)
            pd2(npd) = -COS(deg * k * rtemp0)
            wpd(npd) = .CONVERT."1" / itemp0
            kold = kold + deg1
         ENDDO
         k = itemp0
         npd = npd + 1
         pd1(npd) = -(2 * MOD(DEG,2) - 1)
         pd2(npd) = -pd1(npd)
         wpd(npd) = .CONVERT."1" / (2 * itemp0)

      END SUBROUTINE pdpts

!     ******************************************************************
!
!     PARAMS type                i=input, o=output, a=auxiliary
!
!     DEG    integer             i degree of approximation
!     DEGMAX integer             i maximum degree allowed. DEGMAX + 1 is
!                                  the leading dimension of T1, T2 and
!                                  C0
!     NPD    integer             i number of Padua points
!     PD1    mpfr(NPD)           i array of the first coordinates of
!                                    Padua points
!     PD2    mpfr(NPD)           i array of the second coordinates of
!                                    Padua points
!     WPD    mpfr(NPD)           i array of the weights
!     FPD    mpfr(NPD)           i array of the values of the function
!                                    at Padua points
!     T1     mpfr(0:DEG,NPD)     a Normalized Chebyshev polynomials at
!                                    the first coordinate of Padua points
!     T2     mpfr(0:DEG,NPD)     a Normalized Chebyshev polynomials at
!                                    the second coordinate of Padua points
!     C0     mpfr(0:DEG,0:DEG)   o matrix C0 of the coefficients
!     ESTERR mpfr                o estimated error
!
!     ******************************************************************
      SUBROUTINE padua2(deg,npd,pd1,pd2,wpd,fpd,t1,t2,c0,esterr)
         INTEGER  ::  deg,npd
         TYPE(mpfr_type) :: esterr
         TYPE(mpfr_type) :: pd1(npd),pd2(npd),wpd(npd),fpd(npd), &
              t1(0:deg,npd),t2(0:deg,npd),c0(0:deg,0:deg),rtmp

         INTEGER ::  i,j,k,retval

         rtmp=1
         DO i = 1,npd
            CALL cheb(deg,pd1(i),t1(0,i),wpd(i) * fpd(i))
            CALL cheb(deg,pd2(i),t2(0,i),rtmp)
         ENDDO

         ! B = T1 * T2'
         rtmp=0
         DO i=0,deg
         DO j=0,deg
            c0(i,j)=0
            DO k=1,npd
               ! reduce the impact of the memory leak problem, code explicitly 
               ! in mpfr operations the O(n^4) flops in this section
               ! c0(i,j)=c0(i,j)+t1(i,k)*t2(j,k)
               retval = mpfr_mul(rtmp,t1(i,k),t2(j,k),GMP_RNDN)
               retval = mpfr_add(c0(i,j),c0(i,j),rtmp,GMP_RNDN)
            ENDDO
         ENDDO
         ENDDO

         ! Halve the correct term in C0
         c0(deg,0) = c0(deg,0) / 2

         esterr = 0
         rtmp = 0
         DO j = 0,2
            DO i = 0,deg - j
               IF (mpfr_cmp(c0(i,deg - i - j),rtmp)>0) THEN
                  esterr = esterr + c0(i,deg - i - j)
               ELSE
                  esterr = esterr - c0(i,deg - i - j)
               ENDIF
            ENDDO
         ENDDO
         esterr = 2 * esterr

      END SUBROUTINE padua2

!     ******************************************************************
!
!     PARAMS type                i=input, o=output, a=auxiliary
!
!     DEG    integer             i degree of approximation
!     DEGMAX integer             i maximum degree allowed. DEGMAX + 1 is
!                                  the leading dimension of C0
!     C0     double(0:DEG,0:DEG) i matrix C0 of the coefficients
!     TG1    double              i first coordinate of the
!                                   target point
!     TG2    double              i second coordinate of the
!                                   target point
!     T1     double(0:DEG)       a Normalized Chebyshev polynomials at
!                                  the first coordinate of target point
!     T2     double(0:DEG)       a Normalized Chebyshev polynomials at
!                                  the second coordinate of target point
!
!     **********************************************************************
      FUNCTION pd2val(deg,c0,tg1,tg2,t1,t2) RESULT(res)
         INTEGER :: deg
         TYPE(mpfr_type) :: res,tg1,tg2,c0(0:deg,0:deg),t1(0:deg),t2(0:deg)

         INTEGER :: i,k
         TYPE(mpfr_type) :: dot,one
        
         one = 1 
         CALL cheb(deg,tg1,t1,one)
         CALL cheb(deg,tg2,t2,one)
 
         res = 0
         DO i = deg,0,-1
            dot = 0
            DO k=0,i
               dot = dot + t1(k)*c0(k,deg-i)
            ENDDO
            res = res + dot * t2(deg-i)
         ENDDO

      END FUNCTION pd2val

END MODULE padua2_mpfr

