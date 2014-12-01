!-----------------------------------------------------------------------------!
! Some F77 style routines to get full performance
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param za ...
!> \param cmat ...
!> \param zb ...
!> \param ex ...
!> \param ey ...
!> \param ez ...
!> \param scr ...
! *****************************************************************************
SUBROUTINE rankup ( n, za, cmat, zb, ex, ey, ez, scr )
!
! cmat(i,j,k) <- za * cmat(i,j,k) + ex(i) * ey(j) * ez(k)
!
  USE kinds,                        ONLY: dp
  IMPLICIT NONE

  INTEGER, DIMENSION ( 3 ), INTENT ( IN ) :: n
  COMPLEX (KIND=dp), INTENT ( IN ) :: za, zb
  COMPLEX (KIND=dp), INTENT ( IN ) :: ex ( * )
  COMPLEX (KIND=dp), INTENT ( IN ) :: ey ( * )
  COMPLEX (KIND=dp), INTENT ( IN ) :: ez ( * )
  COMPLEX (KIND=dp), INTENT ( INOUT ) :: cmat ( * )
  COMPLEX (KIND=dp), INTENT ( INOUT ) :: scr ( * )

  INTEGER :: n2, n3
  COMPLEX (KIND=dp) :: cone

  cone = CMPLX ( 1._dp, KIND = dp )
  n2 = n ( 1 ) * n ( 2 )
  n3 = n2 * n ( 3 )
  scr ( 1:n2 ) = CMPLX ( 0._dp, KIND = dp )
  CALL zgeru ( n ( 1 ), n ( 2 ), zb, ex, 1, ey, 1, scr, n ( 1 ) )
  CALL zscal ( n3, za, cmat, 1 )
  CALL zgeru ( n2, n ( 3 ), cone, scr, 1, ez, 1, cmat, n2 )

END SUBROUTINE rankup

!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param cvec2 ...
!> \param cvec ...
! *****************************************************************************
SUBROUTINE vc_x_vc ( n, cvec2, cvec )
!
! cden(i) <- cden(i) * cden2(i)
!
  USE kinds,                        ONLY: dp
  IMPLICIT NONE

  INTEGER, INTENT ( IN ) :: n
  REAL (KIND=dp), INTENT ( IN ) :: cvec2 (2,*)
  REAL (KIND=dp), INTENT ( INOUT ) :: cvec (2,*)

  INTEGER :: i
  REAL (KIND=dp) :: re_c, im_c, re_c2, im_c2

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(re_c,im_c,re_c2,im_c2) SHARED(cvec,cvec2,n)
  DO i = 1, n
    re_c = cvec ( 1, i )
    im_c = cvec ( 2, i )
    re_c2 = cvec2 ( 1, i )
    im_c2 = cvec2 ( 2, i )
    cvec ( 1, i ) = re_c * re_c2 - im_c * im_c2
    cvec ( 2, i ) = im_c * re_c2 + re_c * im_c2
  END DO

END SUBROUTINE vc_x_vc

!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param rvec ...
!> \param cvec ...
! *****************************************************************************
SUBROUTINE vr_x_vc ( n, rvec, cvec )
!
! cden(i) <- cden(i) * rden(i)
!
  USE kinds,                        ONLY: dp
  IMPLICIT NONE

  INTEGER, INTENT ( IN ) :: n
  REAL (KIND=dp), INTENT ( IN ) :: rvec (*)
  REAL (KIND=dp), INTENT ( INOUT ) :: cvec (2,*)

  INTEGER :: i

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(cvec,rvec,n) PRIVATE(i)
  DO i = 1, n
    cvec ( 1, i ) = cvec ( 1, i ) * rvec ( i )
    cvec ( 2, i ) = cvec ( 2, i ) * rvec ( i )
  END DO

END SUBROUTINE vr_x_vc

!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param z ...
!> \param r1 ...
!> \param r2 ...
! *****************************************************************************
SUBROUTINE copy_cri ( n, z, r1, r2 )
!
! r1 = real ( z )
! r2 = imag ( z )
!
  USE kinds,                        ONLY: dp
  IMPLICIT NONE

  INTEGER, INTENT ( IN ) :: n
  REAL (KIND=dp), INTENT ( OUT ) :: r1(*), r2(*)
  REAL (KIND=dp), INTENT ( IN ) :: z(2,*)

  INTEGER :: i

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(r1,r2,z,n) PRIVATE(i)
  DO i = 1, n
    r1 ( i ) = z ( 1, i )
    r2 ( i ) = z ( 2, i )
  END DO

END SUBROUTINE copy_cri

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param z ...
!> \param r1 ...
! *****************************************************************************
SUBROUTINE copy_cr ( n, z, r1 )
!
! r1 = real ( z )
!
  USE kinds,                        ONLY: dp
  IMPLICIT NONE

  INTEGER, INTENT ( IN ) :: n
  REAL (KIND=dp), INTENT ( OUT ) :: r1(n)
  REAL (KIND=dp), INTENT ( IN ) :: z(2,n)
  INTEGER :: i

!$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(r1,z,n)
  DO i=1,n
       r1(i)=z(1,i)
  ENDDO

END SUBROUTINE copy_cr

!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param r1 ...
!> \param z ...
! *****************************************************************************
SUBROUTINE copy_rc ( n, r1, z )
!
! z = r1
!
  USE kinds,                        ONLY: dp
  IMPLICIT NONE

  INTEGER, INTENT ( IN ) :: n
  REAL (KIND=dp), INTENT ( IN ) :: r1(n)
  REAL (KIND=dp), INTENT ( OUT ) :: z(2,n)
  INTEGER :: i
!$OMP PARALLEL DO DEFAULT(none), PRIVATE(i), SHARED(n,z,r1)
  DO i=1,n
       z(1,i)=r1(i)
       z(2,i)=0._dp
  ENDDO

END SUBROUTINE copy_rc

!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param z ...
! *****************************************************************************
SUBROUTINE zero_c ( n, z)
!
! z = ( 0.0_dp , 0.0_dp)
!

  USE kinds,                        ONLY: dp

  INTEGER, INTENT (IN)     :: n
  REAL (KIND=dp), INTENT(OUT) :: z(2,n)
  INTEGER                  :: i

!$OMP PARALLEL DO default(none), private(i), shared(z,n)
  DO i=1,n
     z(1,i) = 0.0_dp
     z(2,i) = 0.0_dp
  ENDDO

END SUBROUTINE zero_c

! *****************************************************************************
!-----------------------------------------------------------------------------!
