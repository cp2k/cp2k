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

  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fft3d, mfft2d, mltfft
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


    INTEGER, DIMENSION(:), POINTER           :: data
    INTEGER, INTENT(IN)                      :: max_length

    INTEGER, PARAMETER                       :: rlen = 81
    INTEGER, DIMENSION(rlen), PARAMETER :: radix = (/ &
         &   2,   4,   6,   8,   9,  12,  15,  16,  18,  20,  24,  25,&
         &  27,  30,  32,  36,  40,  45,  48,  54,  60,  64,  72,  75,  80,  81,  90,&
         &  96, 100, 108, 120, 125, 128, 135, 144, 150, 160, 162, 180, 192, 200, 216,&
         & 225, 240, 243, 256, 270, 288, 300, 320, 324, 360, 375, 384, 400, 405, 432,&
         & 450, 480, 486, 500, 512, 540, 576, 600, 625, 640, 648, 675, 720, 729, 750,&
         & 768, 800, 810, 864, 900, 960, 972, 1000, 1024 /)

    INTEGER                                  :: allocstat, i, ndata

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
!!****** fftsg_lib/fft3d [1.0] *
!!
!!   NAME
!!     fft3d
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


    INTEGER, INTENT(INOUT)                   :: fsign
    REAL(KIND=dp), INTENT(IN), OPTIONAL          :: scale
    INTEGER, DIMENSION(:), INTENT(IN)        :: n
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: zg
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL                :: zg_out

    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: xf, yf
    INTEGER                                  :: ierr, ldox, ldoy, ldoz, ldx, &
                                                ldy, ldz, nx, ny, nz, sign_fft
    LOGICAL                                  :: fft_in_place

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

  ldx = SIZE ( zg, 1 )
  ldy = SIZE ( zg, 2 )
  ldz = SIZE ( zg, 3 )
  
#if defined ( __FFTSG )

  IF ( fft_in_place ) THEN

    ALLOCATE ( xf ( ldx*ldy*ldz ), STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "fft3d", "xf", ldx*ldy*ldz )
    ALLOCATE ( yf ( ldx*ldy*ldz ), STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "fft3d", "yf", ldx*ldy*ldz )

    CALL mltfftsg ( 'N', 'T', zg, ldx, ldy*ldz, xf, ldy*ldz, ldx, nx, &
                    ldy*ldz, sign_fft, 1.0_dp )
    CALL mltfftsg ( 'N', 'T', xf, ldy, ldx*ldz, yf, ldx*ldz, ldy, ny, &
                    ldx*ldz, sign_fft, 1.0_dp )
    IF (PRESENT ( scale ) ) THEN
      CALL mltfftsg ( 'N', 'T', yf, ldz, ldy*ldx, zg, ldy*ldx, ldz, nz, &
                      ldy*ldx, sign_fft, scale)
    ELSE
      CALL mltfftsg ( 'N', 'T', yf, ldz, ldy*ldx, zg, ldy*ldx, ldz, nz, &
                      ldy*ldx, sign_fft, 1.0_dp)
    END IF

    DEALLOCATE ( xf, STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "fft3d", "xf" )
    DEALLOCATE ( yf, STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "fft3d", "yf" )

  ELSE

    ldox = SIZE ( zg_out, 1 )
    ldoy = SIZE ( zg_out, 2 )
    ldoz = SIZE ( zg_out, 3 )

    ALLOCATE ( xf ( ldx*ldy*ldz ), STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "fft3d", "xf", ldx*ldy*ldz )

    CALL mltfftsg ( 'N', 'T', zg, ldx, ldy*ldz, zg_out, ldy*ldz, ldx, nx, &
                    ldy*ldz, sign_fft, 1.0_dp )
    CALL mltfftsg ( 'N', 'T', zg_out, ldy, ldx*ldz, xf, ldx*ldz, ldy, ny, &
                    ldx*ldz, sign_fft, 1.0_dp )
    IF (PRESENT ( scale ) ) THEN
      CALL mltfftsg ( 'N', 'T', xf, ldz, ldy*ldx, zg_out, ldy*ldx, ldz, nz, &
                      ldy*ldx, sign_fft, scale)
    ELSE
      CALL mltfftsg ( 'N', 'T', xf, ldz, ldy*ldx, zg_out, ldy*ldx, ldz, nz, &
                      ldy*ldx, sign_fft, 1.0_dp)
    END IF

    DEALLOCATE ( xf, STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "fft3d", "xf" )

  END IF

#else

  fsign = 0

#endif

END SUBROUTINE fft3d

!!*****
!******************************************************************************
!!****** fftsg_lib/mfft2d [1.0] *
!!
!!   NAME
!!     mfft2d
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (15-Jan-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE mfft2d ( rin, rout, fsign, scale, n1, n2, nm, zin, zout )

  IMPLICIT NONE

! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: rin, rout
  INTEGER, INTENT ( INOUT ) :: fsign
  REAL (KIND=dp), INTENT ( IN ) :: scale
  INTEGER, INTENT ( IN ) :: n1, n2, nm
  COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT ( IN ) :: zin
  COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT ( OUT ) :: zout
  
! Locals
  COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: xf
  INTEGER :: li1, li2, li3, lo1, lo2, lo3, ll, ierr, sign_fft
  
!------------------------------------------------------------------------------
  
#if defined ( __FFTSG )

  sign_fft = fsign

  li1 = SIZE ( zin, 1 )
  li2 = SIZE ( zin, 2 )
  li3 = SIZE ( zin, 3 )
  lo1 = SIZE ( zout, 1 )
  lo2 = SIZE ( zout, 2 )
  lo3 = SIZE ( zout, 3 )

  ll = MAX ( li1*li2*li3, lo1*lo2*lo3 )
  ALLOCATE ( xf ( ll ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "mfft2d", "xf", ll )

  IF ( rin(1:1) == "T" .OR.  rin(1:1) == "t" ) THEN
    IF ( rout(1:1) == "T" .OR.  rout(1:1) == "t" ) THEN
!..T-T
      CALL stop_program ( "mfft2d@fftsg_lib","T-T case not implemented")
    ELSE
!..T-N
      CALL mltfftsg ( "T", "N", zin, li1*li2, li3, xf, li3, li1*li2, &
                      n2, li1*li2, sign_fft, 1.0_dp )
      CALL mltfftsg ( "T", "N", xf, li1*li3, li2, zout, lo1, lo2*lo3, &
                      n1, li1*li3, sign_fft, scale )
    END IF
  ELSE
    IF ( rout(1:1) == "T" .OR.  rout(1:1) == "t" ) THEN
!..N-T
      CALL mltfftsg ( "N", "T", zin, li1, li2*li3, xf, li2*li3, li1, &
                      n1, li2*li3, sign_fft, 1.0_dp )
      CALL mltfftsg ( "N", "T", xf, li2, li1*li3, zout, lo1*lo2, lo3, &
                      n2, li1*li3, sign_fft, scale )
    ELSE
!..N-N
      CALL stop_program ( "mfft2d@fftsg_lib","N-N case not implemented")
    END IF
  END IF

  DEALLOCATE ( xf, STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "mfft2d", "xf" )

#else

  fsign = 0

#endif

END SUBROUTINE mfft2d

!!*****
!******************************************************************************
!!****** fftsg_lib/mltfft [1.0] *
!!
!!   NAME
!!     mltfft
!!
!!   FUNCTION
!!     Calls multiple 1d FFT from the SG library
!!
!!   AUTHOR
!!     JGH (8-Jan-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE mltfft ( tin, tout, fsign, scale, n, m, zin, zout )
  
  IMPLICIT NONE

! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: tin, tout
  INTEGER, INTENT ( INOUT ) :: fsign
  REAL (KIND=dp), INTENT ( IN ) :: scale
  INTEGER, INTENT ( IN ) :: n, m
  COMPLEX (KIND=dp), DIMENSION(:,:), INTENT ( IN ) :: zin
  COMPLEX (KIND=dp), DIMENSION(:,:), INTENT ( OUT ) :: zout
  
! Locals
  INTEGER :: ldi, lmi, ldo, lmo
  
!------------------------------------------------------------------------------
  
#if defined ( __FFTSG )

  ldi = SIZE ( zin, 1 )
  lmi = SIZE ( zin, 2 )
  ldo = SIZE ( zout, 1 )
  lmo = SIZE ( zout, 2 )

  CALL mltfftsg ( tin, tout, zin, ldi, lmi, zout, ldo, lmo, &
                  n, m, fsign, scale )

#else

  fsign = 0

#endif

END SUBROUTINE mltfft

!!*****
!******************************************************************************

END MODULE fftsg_lib

!******************************************************************************
