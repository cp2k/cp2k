!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
MODULE fftsg_lib
  USE fft_kinds,                       ONLY: dp
  USE mltfftsg_tools,                  ONLY: mltfftsg

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fftsg_lib'

  PUBLIC :: fftsg_do_init, fftsg_do_cleanup, fftsg_get_lengths, fftsg3d, fftsg1dm

CONTAINS

! *****************************************************************************
!> \brief ...
! *****************************************************************************
  SUBROUTINE fftsg_do_init()

   ! no init needed

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
! *****************************************************************************
  SUBROUTINE fftsg_do_cleanup()
  
   ! no cleanup needed
  
  END SUBROUTINE
  
  
! *****************************************************************************
!> \brief ...
!> \param DATA ...
!> \param max_length ...
!> \par History
!>      Adapted to new interface structure
!> \author JGH
! *****************************************************************************
  SUBROUTINE fftsg_get_lengths ( DATA, max_length )
  
    INTEGER, DIMENSION(*)                    :: DATA
    INTEGER, INTENT(INOUT)                   :: max_length

    INTEGER, PARAMETER                       :: rlen = 81
    INTEGER, DIMENSION(rlen), PARAMETER :: radix = (/    2,   4,   6,   8,   9&
      ,  12,  15,  16,  18,  20,  24,  25,  27,  30,  32,  36,  40,  45,  48, &
      54,  60,  64,  72,  75,  80,  81,  90,  96, 100, 108, 120, 125, 128, 135&
      , 144, 150, 160, 162, 180, 192, 200, 216, 225, 240, 243, 256, 270, 288, &
      300, 320, 324, 360, 375, 384, 400, 405, 432, 450, 480, 486, 500, 512, &
      540, 576, 600, 625, 640, 648, 675, 720, 729, 750, 768, 800, 810, 864, &
      900, 960, 972, 1000, 1024 /)

    INTEGER                                  :: ndata

!------------------------------------------------------------------------------

    ndata = MIN ( max_length, rlen )
    DATA ( 1:ndata ) = RADIX ( 1:ndata )
    max_length = ndata
  
  END SUBROUTINE fftsg_get_lengths
  
! *****************************************************************************
!> \brief ...
!> \param fft_in_place ...
!> \param fsign ...
!> \param scale ...
!> \param n ...
!> \param zin ...
!> \param zout ...
! *****************************************************************************
  SUBROUTINE fftsg3d ( fft_in_place, fsign, scale, n, zin, zout )
  
    LOGICAL, INTENT(IN)                      :: fft_in_place
    INTEGER, INTENT(INOUT)                   :: fsign
    REAL(KIND=dp), INTENT(IN)                :: scale
    INTEGER, DIMENSION(*), INTENT(IN)        :: n
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: zin, zout

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: xf, yf
    INTEGER                                  :: nx, ny, nz

!------------------------------------------------------------------------------

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
  
  END SUBROUTINE fftsg3d
  
! *****************************************************************************
!> \brief ...
!> \param fsign ...
!> \param trans ...
!> \param n ...
!> \param m ...
!> \param zin ...
!> \param zout ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE fftsg1dm ( fsign, trans, n, m, zin, zout, scale )
  
    INTEGER, INTENT(INOUT)                   :: fsign
    LOGICAL, INTENT(IN)                      :: trans
    INTEGER, INTENT(IN)                      :: n, m
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: zin
    COMPLEX(KIND=dp), DIMENSION(*), &
      INTENT(OUT)                            :: zout
    REAL(KIND=dp), INTENT(IN)                :: scale

!------------------------------------------------------------------------------

    IF ( trans ) THEN
      IF ( fsign > 0 ) THEN
        CALL mltfftsg ( "T", "N", zin, m, n, zout, n, m, n, m, fsign, scale )
      ELSE
        CALL mltfftsg ( "N", "T", zin, n, m, zout, m, n, n, m, fsign, scale )
      END IF
    ELSE
      CALL mltfftsg ( "N", "N", zin, n, m, zout, n, m, n, m, fsign, scale )
    END IF
  
  END SUBROUTINE fftsg1dm

END MODULE

