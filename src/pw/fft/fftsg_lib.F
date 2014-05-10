!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

SUBROUTINE fftsg_do_cleanup()

 ! no cleanup needed

END SUBROUTINE


! *****************************************************************************
!> \par History
!>      Adapted to new interface structure
!> \author JGH
! *****************************************************************************
SUBROUTINE fftsg_get_lengths ( DATA, max_length )

  IMPLICIT NONE
  INTEGER, DIMENSION(*)                    :: DATA
  INTEGER, INTENT(INOUT)                   :: max_length

  INTEGER, PARAMETER                       :: rlen = 81
  INTEGER, DIMENSION(rlen), PARAMETER :: radix = (/ &
  &   2,   4,   6,   8,   9,  12,  15,  16,  18,  20,  24,  25,&
  &  27,  30,  32,  36,  40,  45,  48,  54,  60,  64,  72,  75,  80,  81,  90,&
  &  96, 100, 108, 120, 125, 128, 135, 144, 150, 160, 162, 180, 192, 200, 216,&
  & 225, 240, 243, 256, 270, 288, 300, 320, 324, 360, 375, 384, 400, 405, 432,&
  & 450, 480, 486, 500, 512, 540, 576, 600, 625, 640, 648, 675, 720, 729, 750,&
  & 768, 800, 810, 864, 900, 960, 972, 1000, 1024 /)

  INTEGER                                  :: ndata

!------------------------------------------------------------------------------

  ndata = MIN ( max_length, rlen )
  DATA ( 1:ndata ) = RADIX ( 1:ndata )
  max_length = ndata

END SUBROUTINE fftsg_get_lengths

! *****************************************************************************
SUBROUTINE fftsg3d ( fft_in_place, fsign, scale, n, zin, zout )

  USE fft_kinds, ONLY                                   : dp
  LOGICAL, INTENT(IN)                                  :: fft_in_place
  INTEGER, INTENT(INOUT)                               :: fsign
  REAL(KIND=dp), INTENT(IN)                            :: scale
  INTEGER, DIMENSION(*), INTENT(IN)                    :: n
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)        :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)        :: zout

  INTEGER                                              :: nx, ny, nz
  COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)          :: xf, yf

!------------------------------------------------------------------------------

#if defined ( __FFTSG )

  nx = n ( 1 )
  ny = n ( 2 )
  nz = n ( 3 )

  IF ( fft_in_place ) THEN

    ALLOCATE ( xf(nx*ny*nz), yf(nx*ny*nz) )

    CALL mltfftsg ( 'N', 'T', zin, nx, ny*nz, xf, ny*nz, nx, nx, &
                    ny*nz, fsign, 1.0_dp )
    CALL mltfftsg ( 'N', 'T', xf, ny, nx*nz, yf, nx*nz, ny, ny, &
                    nx*nz, fsign, 1.0_dp )
    CALL mltfftsg ( 'N', 'T', yf, nz, ny*nx, zin, ny*nx, nz, nz, &
                    ny*nx, fsign, scale)

    DEALLOCATE ( xf, yf )

  ELSE

    ALLOCATE ( xf(nx*ny*nz) )

    CALL mltfftsg ( 'N', 'T', zin, nx, ny*nz, zout, ny*nz, nx, nx, &
                    ny*nz, fsign, 1.0_dp )
    CALL mltfftsg ( 'N', 'T', zout, ny, nx*nz, xf, nx*nz, ny, ny, &
                    nx*nz, fsign, 1.0_dp )
    CALL mltfftsg ( 'N', 'T', xf, nz, ny*nx, zout, ny*nx, nz, nz, &
                    ny*nx, fsign, scale)

    DEALLOCATE ( xf )

  END IF

#else

  fsign = 0

#endif

END SUBROUTINE fftsg3d

! *****************************************************************************
SUBROUTINE fftsg1dm ( fsign, trans, n, m, zin, zout, scale )

  USE fft_kinds, ONLY                                   : dp
  INTEGER, INTENT ( INOUT )                          :: fsign
  LOGICAL, INTENT ( IN )                             :: trans
  INTEGER, INTENT ( IN )                             :: n, m
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( INOUT )  :: zin
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( OUT )    :: zout
  REAL (KIND=dp), INTENT ( IN )                      :: scale

!------------------------------------------------------------------------------

#if defined ( __FFTSG )

  IF ( trans ) THEN
    IF ( fsign > 0 ) THEN
      CALL mltfftsg ( "T", "N", zin, m, n, zout, n, m, n, m, fsign, scale )
    ELSE
      CALL mltfftsg ( "N", "T", zin, n, m, zout, m, n, n, m, fsign, scale )
    END IF
  ELSE
    CALL mltfftsg ( "N", "N", zin, n, m, zout, n, m, n, m, fsign, scale )
  END IF

#else

  fsign = 0

#endif

END SUBROUTINE fftsg1dm

