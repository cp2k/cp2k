!-----------------------------------------------------------------------------!
! Some F77 style routines to get full performance
!-----------------------------------------------------------------------------!
SUBROUTINE rankup ( n, za, cmat, zb, ex, ey, ez, scr )
!
! cmat(i,j,k) <- za * cmat(i,j,k) + ex(i) * ey(j) * ez(k)
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND ( 14, 200 )
  
  INTEGER, DIMENSION ( 3 ), INTENT ( IN ) :: n
  COMPLEX ( dbl ), INTENT ( IN ) :: za, zb
  COMPLEX ( dbl ), INTENT ( IN ) :: ex ( * )
  COMPLEX ( dbl ), INTENT ( IN ) :: ey ( * )
  COMPLEX ( dbl ), INTENT ( IN ) :: ez ( * )
  COMPLEX ( dbl ), INTENT ( INOUT ) :: cmat ( * )
  COMPLEX ( dbl ), INTENT ( INOUT ) :: scr ( * )

  INTEGER :: n2, n3, i
  COMPLEX ( dbl ) :: cone
  
  cone = CMPLX ( 1._dbl, KIND = dbl )
  n2 = n ( 1 ) * n ( 2 )
  n3 = n2 * n ( 3 )
  scr ( 1:n2 ) = CMPLX ( 0._dbl, KIND = dbl )
#if defined (__sp_lib) 
  CALL cgeru ( n ( 1 ), n ( 2 ), zb, ex, 1, ey, 1, scr, n ( 1 ) )
  CALL cscal ( n3, za, cmat, 1 )
  CALL cgeru ( n2, n ( 3 ), cone, scr, 1, ez, 1, cmat, n2 )
#else
  CALL zgeru ( n ( 1 ), n ( 2 ), zb, ex, 1, ey, 1, scr, n ( 1 ) )
  CALL zscal ( n3, za, cmat, 1 )
  CALL zgeru ( n2, n ( 3 ), cone, scr, 1, ez, 1, cmat, n2 )
#endif

END SUBROUTINE rankup 

!-----------------------------------------------------------------------------!

SUBROUTINE vr_x_vc ( n, rvec, cvec )
!
! cden(i) <- cden(i) * rden(i)
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND ( 14, 200 )

  INTEGER, INTENT ( IN ) :: n
  REAL ( dbl ), INTENT ( IN ) :: rvec (*)
  REAL ( dbl ), INTENT ( INOUT ) :: cvec (2,*)

  INTEGER :: i

  DO i = 1, n
    cvec ( 1, i ) = cvec ( 1, i ) * rvec ( i )
    cvec ( 2, i ) = cvec ( 2, i ) * rvec ( i )
  END DO

END SUBROUTINE vr_x_vc 

!-----------------------------------------------------------------------------!

SUBROUTINE copy_cri ( n, z, r1, r2 )
!
! r1 = real ( z )
! r2 = imag ( z )
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND ( 14, 200 )

  INTEGER, INTENT ( IN ) :: n
  REAL ( dbl ), INTENT ( OUT ) :: r1(*), r2(*)
  REAL ( dbl ), INTENT ( IN ) :: z(2,*)

  INTEGER :: i

  DO i = 1, n
    r1 ( i ) = z ( 1, i )
    r2 ( i ) = z ( 2, i )
  END DO

END SUBROUTINE copy_cri

!-----------------------------------------------------------------------------!

SUBROUTINE copy_cr ( n, z, r1 )
!
! r1 = real ( z )
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND ( 14, 200 )

  INTEGER, INTENT ( IN ) :: n
  REAL ( dbl ), INTENT ( OUT ) :: r1(*)
  REAL ( dbl ), INTENT ( IN ) :: z(2,*)

  INTEGER :: i

  DO i = 1, n
    r1 ( i ) = z ( 1, i )
  END DO

END SUBROUTINE copy_cr

!-----------------------------------------------------------------------------!

