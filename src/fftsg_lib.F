!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fftsg_lib [1.0] *
!!
!!   NAME
!!     fftsg_lib
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE fftsg_lib
  
  USE kinds, ONLY: dbl
  USE stop_program, ONLY : stop_memory, stop_prg
  
  IMPLICIT NONE 
  
  PRIVATE
  PUBLIC :: fft3d
  PUBLIC :: fft_get_lengths
  
!!*****
!******************************************************************************

CONTAINS
  
!******************************************************************************
!!****** fftsg_lib/fft_get_lengths [1.0] *
!!
!!   NAME
!!     fft_get_lengths
!!
!!   FUNCTION
!!     Give the allowed lengths of FFT's   '''
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE fft_get_lengths ( data, max_length )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: max_length
  INTEGER, DIMENSION ( : ), POINTER :: data
  
! Locals
  INTEGER, PARAMETER :: rlen = 82 
  INTEGER, DIMENSION ( rlen ), PARAMETER :: radix = &
   (/ 2, 4, 5, 6, 8, 9, 12, 15, 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, &
     45, 48, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128, &
     135, 144, 150, 160, 162, 180, 192, 200, 216, 225, 240, 243, 256, 270, &
     288, 300, 320, 324, 360, 375, 384, 400, 405, 432, 450, 480, 486, 500, &
     512, 540, 576, 600, 625, 640, 648, 675, 720, 729, 750, 768, 800, 810, &
     864, 900, 960, 972, 1000, 1024 /)
  INTEGER :: i, allocstat, ndata
  
!------------------------------------------------------------------------------
  
  ndata = 0
  DO i = 1, rlen
    IF ( radix ( i ) > max_length ) EXIT
    ndata = ndata + 1
  END DO
  
  ALLOCATE ( data ( ndata ), STAT = allocstat )
  IF ( allocstat /= 0 ) THEN
     CALL stop_memory ( "fft_get_lengths", "data", ndata )
  END IF
  
  data ( 1:ndata ) = radix ( 1:ndata )
  
END SUBROUTINE fft_get_lengths

!!*****
!******************************************************************************
!!****** fftsg_lib/fft_3d [1.0] *
!!
!!   NAME
!!     fft_3d
!!
!!   FUNCTION
!!     Routine with wrapper for all 3dfft call:
!!     Does transform with exp(+ig.r*sign):
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE fft3d ( fsign, scale, n, zg, zg_out )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( INOUT ) :: fsign
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: scale
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: n
  COMPLEX ( dbl ), DIMENSION(:,:,:), INTENT ( INOUT ) :: zg
  COMPLEX ( dbl ), DIMENSION(:,:,:), INTENT ( INOUT ), OPTIONAL :: zg_out
  
! Locals
  INTEGER :: sign_fft, ldx, ldy, ldz, ldox, ldoy, ldoz, ierr
  INTEGER :: nx, ny, nz
  COMPLEX ( dbl ), DIMENSION(:), ALLOCATABLE :: xf, yf
  LOGICAL :: fft_in_place
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( zg_out ) ) THEN
     fft_in_place = .false.
  ELSE
     fft_in_place = .true.
  END IF
  
  sign_fft = fsign

  nx = n ( 1 )
  ny = n ( 2 )
  nz = n ( 3 )

  ldx = SIZE ( zg (:,1,1) )
  ldy = SIZE ( zg (1,:,1) )
  ldz = SIZE ( zg (1,1,:) )
  
#if defined ( __FFTSG )

  IF ( fft_in_place ) THEN

    ALLOCATE ( xf ( ldx*ldy*ldz ), STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "fft3d", "xf", ldx*ldy*ldz )
    ALLOCATE ( yf ( ldx*ldy*ldz ), STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "fft3d", "yf", ldx*ldy*ldz )

    CALL mltfftsg ( 'N', 'T', zg, ldx, ldy*ldz, xf, ldy*ldz, ldx, nx, &
                    ldy*ldz, sign_fft, 1._dbl )
    CALL mltfftsg ( 'N', 'T', xf, ldy, ldx*ldz, yf, ldx*ldz, ldy, ny, &
                    ldx*ldz, sign_fft, 1._dbl )
    CALL mltfftsg ( 'N', 'T', yf, ldz, ldy*ldx, zg, ldy*ldx, ldz, nz, &
                    ldy*ldx, sign_fft, scale)

    DEALLOCATE ( xf, STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "fft3d", "xf" )
    DEALLOCATE ( yf, STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "fft3d", "yf" )

  ELSE

    ldox = SIZE ( zg_out (:,1,1) )
    ldoy = SIZE ( zg_out (1,:,1) )
    ldoz = SIZE ( zg_out (1,1,:) )

    ALLOCATE ( xf ( ldx*ldy*ldz ), STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "fft3d", "xf", ldx*ldy*ldz )

    CALL mltfftsg ( 'N', 'T', zg, ldx, ldy*ldz, zg_out, ldy*ldz, ldx, nx, &
                    ldy*ldz, sign_fft, 1._dbl )
    CALL mltfftsg ( 'N', 'T', zg_out, ldy, ldx*ldz, xf, ldx*ldz, ldy, ny, &
                    ldx*ldz, sign_fft, 1._dbl )
    CALL mltfftsg ( 'N', 'T', xf, ldz, ldy*ldx, zg_out, ldy*ldx, ldz, nz, &
                    ldy*ldx, sign_fft, scale)

    DEALLOCATE ( xf, STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "fft3d", "xf" )

  END IF

#else

  fsign = 0

#endif

END SUBROUTINE fft3d

!!*****
!******************************************************************************

END MODULE fftsg_lib

!******************************************************************************
