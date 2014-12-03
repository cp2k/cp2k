MODULE fast

  USE kinds,                           ONLY: dp

IMPLICIT NONE

PRIVATE

PUBLIC :: rankup, vc_x_vc, vr_x_vc, copy_cri, copy_cr, copy_rc, zero_c

  INTERFACE zero_c
     MODULE PROCEDURE zero_c2, zero_c3
  END INTERFACE

CONTAINS

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

    INTEGER, DIMENSION(3), INTENT(IN)        :: n
    COMPLEX(KIND=dp), INTENT(IN)             :: za
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: cmat
    COMPLEX(KIND=dp), INTENT(IN)             :: zb
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(IN)                             :: ex, ey, ez
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: scr

    COMPLEX(KIND=dp)                         :: cone
    INTEGER                                  :: n2, n3

    cone = CMPLX ( 1._dp, KIND = dp )
    n2 = n ( 1 ) * n ( 2 )
    n3 = n2 * n ( 3 )
    scr ( 1:n2 ) = CMPLX ( 0._dp, KIND = dp )
    CALL zgeru ( n ( 1 ), n ( 2 ), zb, ex, 1, ey, 1, scr, n ( 1 ) )
    CALL zscal ( n3, za, cmat, 1 )
    CALL zgeru ( n2, n ( 3 ), cone, scr, 1, ez, 1, cmat, n2 )

END SUBROUTINE rankup

! *****************************************************************************
!> \brief Multiply two complex 3D arrays element-wise
!> \param cvec2 the other complex array
!> \param cvec the complex array returned as the result 
! *****************************************************************************
SUBROUTINE vc_x_vc (cvec2, cvec )
!
! cvec(i) <- cvec(i) * cvec2(i)
!

    COMPLEX(KIND=dp), INTENT(IN)             :: cvec2(:,:,:)
    COMPLEX(KIND=dp), INTENT(INOUT)          :: cvec(:,:,:)

!$OMP PARALLEL WORKSHARE DEFAULT(none), SHARED(cvec,cvec2)
    cvec(:,:,:) = cvec(:,:,:) * cvec2(:,:,:)
!$OMP END PARALLEL WORKSHARE

END SUBROUTINE vc_x_vc

! *****************************************************************************
!> \brief Scale a complex array element-wise by a real array
!> \param rvec the real array
!> \param cvec the complex array
! *****************************************************************************
SUBROUTINE vr_x_vc ( rvec, cvec )
!
! cvec(i) <- cvec(i) * rvec(i)
!

    REAL(KIND=dp), INTENT(IN)                :: rvec(:,:,:)
    COMPLEX(KIND=dp), INTENT(INOUT)          :: cvec(:,:,:)

!$OMP PARALLEL WORKSHARE DEFAULT(none), SHARED(cvec,rvec)
    cvec(:,:,:) = cvec(:,:,:) * rvec(:,:,:)
!$OMP END PARALLEL WORKSHARE

END SUBROUTINE vr_x_vc

! *****************************************************************************
!> \brief Copy a the real and imag. parts of a complex 3D array into two real arrays
!> \param z the complex array
!> \param r1 the real array for the real part
!> \param r2 the real array for the imaginary part
! *****************************************************************************
SUBROUTINE copy_cri (z, r1, r2 )
!
! r1 = real ( z )
! r2 = imag ( z )
!

    COMPLEX(KIND=dp), INTENT(IN)             :: z(:,:,:)
    REAL(KIND=dp), INTENT(INOUT)             :: r1(:,:,:), r2(:,:,:)

!$OMP PARALLEL WORKSHARE DEFAULT(none), SHARED(r1,r2,z)
    r1(:,:,:) = REAL(z(:,:,:),KIND=dp)
    r2(:,:,:) = AIMAG(z(:,:,:))
!$OMP END PARALLEL WORKSHARE

END SUBROUTINE copy_cri

! *****************************************************************************
!> \brief Copy the real part of a complex 3D array to a real array
!> \param z the complex array
!> \param r1 the real array
! *****************************************************************************
SUBROUTINE copy_cr ( z, r1 )
!
! r1 = real ( z )
!

    COMPLEX(KIND=dp), INTENT(IN)             :: z(:,:,:)
    REAL(KIND=dp), INTENT(INOUT)             :: r1(:,:,:)

!$OMP PARALLEL WORKSHARE DEFAULT(none), SHARED(r1,z)
    r1(:,:,:) = REAL(z(:,:,:),KIND=dp)
!$OMP END PARALLEL WORKSHARE

END SUBROUTINE copy_cr

! *****************************************************************************
!> \brief Copy a real 3D array to complex 3D array
!> \param r1 the real array
!> \param z the complex array
! *****************************************************************************
SUBROUTINE copy_rc ( r1, z )
!
! z = r1
!

    REAL(KIND=dp), INTENT(IN)                :: r1(:,:,:)
    COMPLEX(KIND=dp), INTENT(INOUT)          :: z(:,:,:)

!$OMP PARALLEL WORKSHARE DEFAULT(none), SHARED(r1,z)
    z(:,:,:)=CMPLX(r1(:,:,:), 0.0_dp, KIND=dp)
!$OMP END PARALLEL WORKSHARE

END SUBROUTINE copy_rc

! *****************************************************************************
!> \brief Zero a complex 2D array (optionally with OpenMP)
!> \param z the array
! *****************************************************************************
SUBROUTINE zero_c2 (z)
!
! z = ( 0.0_dp , 0.0_dp)
!

    COMPLEX(KIND=dp), INTENT(INOUT)          :: z(:,:)

!$OMP PARALLEL WORKSHARE DEFAULT(none), SHARED(z)
    z(:,:) = 0.0_dp
!$OMP END PARALLEL WORKSHARE
END SUBROUTINE zero_c2

! *****************************************************************************
!> \brief Zero a complex 3D array (optionally with OpenMP)
!> \param z the array
! *****************************************************************************
SUBROUTINE zero_c3 (z)
!
! z = ( 0.0_dp , 0.0_dp)
!

    COMPLEX(KIND=dp), INTENT(INOUT)          :: z(:,:,:)

!$OMP PARALLEL WORKSHARE DEFAULT(none), SHARED(z)
    z(:,:,:) = 0.0_dp
!$OMP END PARALLEL WORKSHARE

END SUBROUTINE zero_c3

END MODULE fast
