!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/bessel_lib [1.0] *
!!
!!   NAME
!!     bessel_lib
!!
!!   FUNCTION
!!     Calculates Bessel functions
!!
!!   AUTHOR
!!     JGH (10-02-2001) 
!! 
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     Functions adapted from netlib 
!!
!!   SOURCE
!******************************************************************************

MODULE bessel_lib
  
  USE kinds,                           ONLY: dbl

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: bessi0, bessi1, bessj0, bessj1, bessk0, bessk1
  
!******************************************************************************

CONTAINS

!******************************************************************************

  FUNCTION bessj0 ( x )


    REAL(dbl), INTENT(IN)                    :: x
    REAL(dbl)                                :: bessj0

    REAL(dbl), PARAMETER :: p1 = 1._dbl, p2 = -0.1098628627e-2_dbl, &
      p3 = 0.2734510407e-4_dbl, p4 = -0.2073370639e-5_dbl, &
      p5 = 0.2093887211e-6_dbl, q1 = -0.1562499995e-1_dbl, &
      q2 = 0.1430488765e-3_dbl, q3 = -0.6911147651e-5_dbl, &
      q4 = 0.7621095161e-6_dbl, q5 = -0.934945152e-7_dbl, &
      r1 = 57568490574._dbl, r2 = -13362590354._dbl, r3 = 651619640.7_dbl, &
      r4 = -11214424.18_dbl, r5 = 77392.33017_dbl, r6 = -184.9052456_dbl, &
      s1 = 57568490411._dbl, s2 = 1029532985._dbl, s3 = 9494680.718_dbl, &
      s4 = 59272.64853_dbl, s5 = 267.8532712_dbl, s6 = 1._dbl

    REAL(dbl)                                :: ax, xx, y, z

    IF ( ABS ( x ) < 8._dbl ) THEN
      y = x * x
      bessj0 = ( r1 + y * ( r2 + y * ( r3 + y * ( r4 + y * &
        ( r5 + y * r6 ) ) ) ) ) / ( s1 + y * ( s2 + y * ( s3 + y * &
        ( s4 + y * ( s5 + y * s6 ) ) ) ) )
    ELSE
      ax = ABS ( x )
      z = 8._dbl / ax
      y = z * z
      xx = ax - 0.785398164_dbl
      bessj0 = SQRT ( 0.636619772_dbl / ax ) * ( COS ( xx ) * &
        ( p1 + y * ( p2 + y * ( p3 + y * ( p4 + y * p5 ) ) ) ) &
        - z * SIN ( xx ) * ( q1 + y * ( q2 + y * ( q3 + y * &
        ( q4 + y * q5 ) ) ) ) )
    END IF

  END FUNCTION bessj0

!******************************************************************************

  FUNCTION bessj1 ( x )


    REAL(dbl), INTENT(IN)                    :: x
    REAL(dbl)                                :: bessj1

    REAL(dbl), PARAMETER :: p1 = 1._dbl, p2 = -0.183105e-2_dbl, &
      p3 = -0.3516396496e-4_dbl, p4 = 0.2457520174e-5_dbl, &
      p5 = -0.240337019e-6_dbl, q1 = 0.04687499995_dbl, &
      q2 = -0.2002690873e-3_dbl, q3 = 0.8449199096e-5_dbl, &
      q4 = -0.88228987e-6_dbl, q5 = 0.105787412e-6_dbl, &
      r1 = 72362614232._dbl, r2 = -7895059235._dbl, r3 = 242396853.1_dbl, &
      r4 = -2972611.439_dbl, r5 = 15704.48260_dbl, r6 = -30.16036606_dbl, &
      s1 = 144725228442._dbl, s2 = 2300535178._dbl, s3 = 18583304.74_dbl, &
      s4 = 99447.43394_dbl, s5 = 376.9991397_dbl, s6 = 1._dbl

    REAL(dbl)                                :: ax, xx, y, z

      IF ( ABS ( x ) <  8._dbl ) THEN
        y = x * x
        bessj1 = x * ( r1 + y * ( r2 + y * ( r3 + y * ( r4 + y * &
          ( r5 + y * r6 ) ) ) ) ) / ( s1 + y * ( s2 + y * ( s3 + &
          y * ( s4 + y * ( s5 + y * s6 ) ) ) ) )
      ELSE
        ax = ABS ( x )
        z = 8._dbl / ax
        y = z * z
        xx = ax - 2.356194491_dbl
        bessj1 = SQRT ( 0.636619772_dbl / ax ) * ( COS ( xx ) * &
          ( p1 + y * ( p2 + y * ( p3 + y * ( p4 + y * p5 ) ) ) ) &
          - z * SIN ( xx ) * ( q1 + y * ( q2 + y * ( q3 + y * &
          ( q4 + y * q5 ) ) ) ) ) * SIGN ( 1._dbl, x )
      END IF
  END FUNCTION bessj1

!******************************************************************************

  FUNCTION bessk0 ( x )


    REAL(dbl), INTENT(IN)                    :: x
    REAL(dbl)                                :: bessk0

    REAL(dbl), PARAMETER :: p1 = -0.57721566_dbl, p2 = 0.42278420_dbl, &
      p3 = 0.23069756_dbl, p4 = 0.3488590e-1_dbl, p5 = 0.262698e-2_dbl, &
      p6 = 0.10750e-3_dbl, p7 = 0.74e-5_dbl, q1 = 1.25331414_dbl, &
      q2 = -0.7832358e-1_dbl, q3 = 0.2189568e-1_dbl, q4 = -0.1062446e-1_dbl, &
      q5 = 0.587872e-2_dbl, q6 = -0.251540e-2_dbl, q7 = 0.53208e-3_dbl

    REAL(dbl)                                :: y

    IF ( x < 2._dbl ) THEN
      y = x * x / 4._dbl
      bessk0 = ( -LOG ( x / 2._dbl ) * bessi0 ( x ) ) + ( p1 + y * &
        ( p2 + y * ( p3 + y * ( p4 + y * ( p5 + y * ( p6 + y * p7 ) ) ) ) ) )
    ELSE
        y = ( 2._dbl / x )
        bessk0 = ( EXP ( -x ) / SQRT ( x ) ) * ( q1 + y * ( q2 + y * &
          ( q3 + y * ( q4 + y * ( q5 + y * ( q6 + y * q7 ) ) ) ) ) )
    END IF

  END FUNCTION bessk0

!******************************************************************************

  FUNCTION bessk1 ( x )


    REAL(dbl), INTENT(IN)                    :: x
    REAL(dbl)                                :: bessk1

    REAL(dbl), PARAMETER :: p1 = 1._dbl, p2 = 0.15443144_dbl, &
      p3 = -0.67278579_dbl, p4 = -0.18156897_dbl, p5 = -0.1919402e-1_dbl, &
      p6 = -0.110404e-2_dbl, p7 = -0.4686e-4_dbl, q1 = 1.25331414_dbl, &
      q2 = 0.23498619_dbl, q3 = -0.3655620e-1_dbl, q4 = 0.1504268e-1_dbl, &
      q5 = -0.780353e-2_dbl, q6 = 0.325614e-2_dbl, q7 = -0.68245e-3_dbl

    REAL(dbl)                                :: y

    IF ( x < 2._dbl ) THEN
      y = x * x / 4._dbl
      bessk1 = ( LOG ( x / 2._dbl ) * bessi1 ( x ) ) + ( 1._dbl / x ) * &
        ( p1 + y * ( p2 + y * ( p3 + y * ( p4 + y * ( p5 + y * &
        ( p6 + y *p7 ) ) ) ) ) )
    ELSE
      y = 2._dbl / x
      bessk1 = ( EXP ( -x ) / SQRT ( x ) ) * ( q1 + y * ( q2 + y * &
        ( q3 + y * ( q4 + y * ( q5 + y * ( q6 + y * q7 ) ) ) ) ) )
    END IF

  END FUNCTION bessk1 

!******************************************************************************

  FUNCTION bessi0 ( x )


    REAL(dbl), INTENT(IN)                    :: x
    REAL(dbl)                                :: bessi0

    REAL(dbl), PARAMETER :: p1 = 1._dbl, p2 = 3.5156229_dbl, &
      p3 = 3.0899424_dbl, p4 = 1.2067492_dbl, p5 = 0.2659732_dbl, &
      p6 = 0.360768e-1_dbl, p7 = 0.45813e-2_dbl, q1 = 0.39894228_dbl, &
      q2 = 0.1328592e-1_dbl, q3 = 0.225319e-2_dbl, q4 = -0.157565e-2_dbl, &
      q5 = 0.916281e-2_dbl, q6 = -0.2057706e-1_dbl, q7 = 0.2635537e-1_dbl, &
      q8 = -0.1647633e-1_dbl, q9 = 0.392377e-2_dbl

    REAL(dbl)                                :: ax, y

    IF ( ABS ( x ) < 3.75_dbl ) THEN
      y = ( x / 3.75_dbl ) ** 2
      bessi0 = p1 + y * ( p2 + y * ( p3 + y * ( p4 + y * &
        ( p5 + y * ( p6 + y * p7 ) ) ) ) )
    ELSE
      ax = ABS ( x )
      y = 3.75_dbl / ax
      bessi0 = ( EXP ( ax ) / SQRT ( ax ) ) * ( q1 + y * ( q2 + y * &
        ( q3 + y * ( q4 + y * ( q5 + y * ( q6 + y * ( q7 + y * &
        ( q8 + y * q9 ) ) ) ) ) ) ) )
    END IF

  END FUNCTION bessi0 

!******************************************************************************

  FUNCTION bessi1 ( x )


    REAL(dbl), INTENT(IN)                    :: x
    REAL(dbl)                                :: bessi1

    REAL(dbl), PARAMETER :: p1 = 0.5_dbl, p2 = 0.87890594_dbl, &
      p3 = 0.51498869_dbl, p4 = 0.15084934e0_dbl, p5 = 0.2658733e-1_dbl, &
      p6 = 0.301532e-2_dbl, p7 = 0.32411e-3_dbl, q1 = 0.39894228_dbl, &
      q2 = -0.3988024e-1_dbl, q3 = -0.362018e-2_dbl, q4 = 0.163801e-2_dbl, &
      q5 = -0.1031555e-1_dbl, q6 = 0.2282967e-1_dbl, q7 = -0.2895312e-1_dbl, &
      q8 = 0.1787654e-1_dbl, q9 = -0.420059e-2_dbl

    REAL(dbl)                                :: ax, y

    IF ( ABS ( x ) < 3.75_dbl ) THEN
      y = ( x / 3.75_dbl ) ** 2
      bessi1 = p1 + y * ( p2 + y * ( p3 + y * ( p4 + y * &
        ( p5 + y * ( p6 + y * p7 ) ) ) ) )
    ELSE
      ax = ABS ( x )
      y = 3.75_dbl / ax
      bessi1 = ( EXP ( ax ) / SQRT ( ax ) ) * ( q1 + y * ( q2 + y * &
        ( q3 + y * ( q4 + y * ( q5 + y * ( q6 + y * ( q7 + y * &
        ( q8 + y * q9 ) ) ) ) ) ) ) )
      IF ( x < 0._dbl ) bessi1 = -bessi1
    END IF

  END FUNCTION bessi1

!!*****
!******************************************************************************

END MODULE bessel_lib

!******************************************************************************
