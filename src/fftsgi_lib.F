!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fftsgi_lib [1.0] *
!!
!!   NAME
!!     fftsgi_lib
!!
!!   FUNCTION
!!     Interface to FFT routines in the complib.sgimath Library of SGI
!!
!!   AUTHOR
!!     JGH 5-Jan-2001
!!
!!   MODIFICATION HISTORY
!!     JGH 20-Feb-2001 : mltfft, mfft2d routines added
!!
!!   SOURCE
!******************************************************************************

MODULE fftsgi_lib

  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE util,                            ONLY: sort

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fft3d, mltfft, mfft2d
  PUBLIC :: fft_get_lengths

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** fftsgi_lib/fft_get_lengths [1.0] *
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

    INTEGER :: allocstat, h, i, j, k, m, maxn, maxn_elevens, maxn_fives, &
      maxn_sevens, maxn_thirteens, maxn_threes, maxn_twos, ndata, nmax, number
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: idx

!------------------------------------------------------------------------------
! compute ndata
! FFTSGI can do arbitrary(?) lenghts, maybe you want to limit them to some
!    powers of small prime numbers though...
! For the moment we keep it the same as for the ESSL

  maxn_twos = 15
  maxn_threes = 2
  maxn_fives = 1
  maxn_sevens = 1
  maxn_elevens = 1
  maxn_thirteens = 0
  maxn = MIN ( max_length, 37748736 )

  ndata = 0
  DO h = 0, maxn_twos
     nmax = HUGE(0) / 2**h
     DO i = 0, maxn_threes
        DO j = 0, maxn_fives
           DO k = 0, maxn_sevens
              DO m = 0, maxn_elevens
                 number = (3**i) * (5**j) * (7**k) * (11**m)

                 IF ( number > nmax ) CYCLE

                 number = number * 2 ** h
                 IF ( number >= maxn ) CYCLE

                 ndata = ndata + 1
              END DO
           END DO
        END DO
     END DO
  END DO

  ALLOCATE ( data ( ndata ), idx ( ndata ), STAT = allocstat )
  IF ( allocstat /= 0 ) THEN
     CALL stop_memory ( "fft_get_lengths", "data, idx", 2*ndata )
  END IF

  ndata = 0
  data ( : ) = 0
  DO h = 0, maxn_twos
     nmax = HUGE(0) / 2**h
     DO i = 0, maxn_threes
        DO j = 0, maxn_fives
           DO k = 0, maxn_sevens
              DO m = 0, maxn_elevens
                 number = (3**i) * (5**j) * (7**k) * (11**m)

                 IF ( number > nmax ) CYCLE

                 number = number * 2 ** h
                 IF ( number >= maxn ) CYCLE

                 ndata = ndata + 1
                 data ( ndata ) = number
              END DO
           END DO
        END DO
     END DO
  END DO

  CALL sort ( data, ndata, idx )

  DEALLOCATE ( idx, STAT = allocstat )
  IF ( allocstat /= 0 ) THEN
     CALL stop_memory ( "fft_get_lengths", "idx" )
  END IF

END SUBROUTINE fft_get_lengths

!!*****
!******************************************************************************
!!****** fftsgi_lib/fft_3d [1.0] *
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


    INTEGER, INTENT(INOUT)                   :: fsign
    REAL(dbl), INTENT(IN), OPTIONAL          :: scale
    INTEGER, DIMENSION(:), INTENT(IN)        :: n
    COMPLEX(dbl), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: zg
    COMPLEX(dbl), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL                :: zg_out

    COMPLEX(dbl), ALLOCATABLE, DIMENSION(:)  :: coef
    INTEGER                                  :: isos, ldox, ldoy, ldoz, ldx, &
                                                ldy, ldz, ncoef, nx, ny, nz, &
                                                sign_fft
    LOGICAL                                  :: fft_in_place

!------------------------------------------------------------------------------

  IF ( PRESENT ( zg_out ) ) THEN
     fft_in_place = .false.
  ELSE
     fft_in_place = .true.
  END IF

  sign_fft = -fsign

  nx = n ( 1 )
  ny = n ( 2 )
  nz = n ( 3 )

  ldx = SIZE ( zg, 1 )
  ldy = SIZE ( zg, 2 )
  ldz = SIZE ( zg, 3 )
  
  ncoef = ( nx + 15 ) + ( ny + 15 ) + ( nz + 15 )

#if defined ( __FFTSGI )

  ALLOCATE ( coef ( ncoef ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fft3d','coef',ncoef )

  IF ( fft_in_place ) THEN

    IF ( scale /= 1._dbl ) call zscal3d ( nx, ny, nz, scale, zg, ldx, ldy )
    CALL zfft3di ( nx, ny, nz, coef )
    CALL zfft3d ( sign_fft, nx, ny, nz, zg, ldx, ldy, coef )

  ELSE

    ldox = SIZE ( zg_out, 1 )
    ldoy = SIZE ( zg_out, 2 )
    ldoz = SIZE ( zg_out, 3 )

    IF ( scale /= 1._dbl ) THEN
      zg_out ( 1:nx, 1:ny, 1:nz ) = zg ( 1:nx, 1:ny, 1:nz )
      call zscal3d ( nx, ny, nz, scale, zg_out, ldox, ldoy )
    ELSE
      zg_out ( 1:nx, 1:ny, 1:nz ) = zg ( 1:nx, 1:ny, 1:nz )
    END IF

    CALL zfft3di ( nx, ny, nz, coef )
    CALL zfft3d ( sign_fft, nx, ny, nz, zg_out, ldox, ldoy, coef )

  END IF

  DEALLOCATE ( coef, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fft3d','coef' )

#else

  fsign = 0

#endif

END SUBROUTINE fft3d

!!*****
!******************************************************************************
!!****** fftsgi_lib/mfft2d [1.0] *
!!
!!   NAME
!!     mfft2d
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (11-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE mfft2d ( rin, rout, fsign, scale, n1, n2, nm, zin, zout )


    CHARACTER(LEN=*), INTENT(IN)             :: rin, rout
    INTEGER, INTENT(INOUT)                   :: fsign
    REAL(dbl), INTENT(IN)                    :: scale
    INTEGER, INTENT(IN)                      :: n1, n2, nm
    COMPLEX(dbl), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zin
    COMPLEX(dbl), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: zout

    COMPLEX(dbl), ALLOCATABLE, DIMENSION(:)  :: coef
    COMPLEX(dbl), ALLOCATABLE, &
      DIMENSION(:, :)                        :: xf
    INTEGER                                  :: i, ierr, ncoef, sign_fft

!------------------------------------------------------------------------------

  sign_fft = -fsign

  ncoef = ( n1 + 15 ) + ( n2 + 15 ) 

#if defined ( __FFTSGI )

  ALLOCATE ( coef ( ncoef ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'mfft2d','coef',ncoef )
  ALLOCATE ( xf ( n1, n2 ), STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "mfft2d", "xf", n1*n2 )

  IF ( rin(1:1) == "T" .OR.  rin(1:1) == "t" ) THEN
    IF ( rout(1:1) == "T" .OR.  rout(1:1) == "t" ) THEN
!..T-T
      CALL stop_program ( "mfft2d@fftsgi_lib","T-T case not implemented")
    ELSE
!..T-N
      DO i = 1, nm
         IF ( scale /= 1._dbl ) THEN
            xf ( 1:n1, 1:n2 ) = scale * zin ( i, 1:n1, 1:n2 )
         ELSE
            xf ( 1:n1, 1:n2 ) = zin ( i, 1:n1, 1:n2 )
         END IF
         CALL zfft2di ( n1, n2, coef )
         CALL zfft2d ( sign_fft, n1, n2, xf, n1, coef )
         zout ( 1:n1, 1:n2, i ) = xf ( 1:n1, 1:n2 )
      END DO
    END IF
  ELSE
    IF ( rout(1:1) == "T" .OR.  rout(1:1) == "t" ) THEN
!..N-T
      DO i = 1, nm
         IF ( scale /= 1._dbl ) THEN
            xf ( 1:n1, 1:n2 ) = scale * zin ( 1:n1, 1:n2, i )
         ELSE
            xf ( 1:n1, 1:n2 ) = zin ( 1:n1, 1:n2, i )
         END IF
         CALL zfft2di ( n1, n2, coef )
         CALL zfft2d ( sign_fft, n1, n2, xf, n1, coef )
         zout ( i, 1:n1, 1:n2 ) = xf ( 1:n1, 1:n2 )
      END DO
    ELSE
!..N-N
      CALL stop_program ( "mfft2d@fftsgi_lib","N-N case not implemented")
    END IF
  END IF

  DEALLOCATE ( xf, STAT = ierr )
  IF ( ierr /= 0 ) call stop_memory ( "mfft2d", "xf" )
  DEALLOCATE ( coef, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'mfft2d','coef' )

#else

  fsign = 0

#endif

END SUBROUTINE mfft2d

!!*****
!******************************************************************************
!!****** fftsgi_lib/mltfft [1.0] *
!!
!!   NAME
!!     mltfft
!!
!!   FUNCTION
!!     Calls multiple 1d FFT from the SGI library
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
  

    CHARACTER(LEN=*), INTENT(IN)             :: tin, tout
    INTEGER, INTENT(INOUT)                   :: fsign
    REAL(dbl), INTENT(IN)                    :: scale
    INTEGER, INTENT(IN)                      :: n, m
    COMPLEX(dbl), DIMENSION(:, :), &
      INTENT(IN)                             :: zin
    COMPLEX(dbl), DIMENSION(:, :), &
      INTENT(OUT)                            :: zout

    COMPLEX(dbl), ALLOCATABLE, DIMENSION(:)  :: coef
    COMPLEX(dbl), ALLOCATABLE, &
      DIMENSION(:, :)                        :: zz
    INTEGER                                  :: isos, ldi, ldo, lmi, lmo, &
                                                ncoef

!------------------------------------------------------------------------------

  ldi = SIZE ( zin, 1 )
  lmi = SIZE ( zin, 2 )
  ldo = SIZE ( zout, 1 )
  lmo = SIZE ( zout, 2 )
  
#if defined ( __FFTSGI )

  ncoef = n + 15
  ALLOCATE ( coef ( ncoef ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fft2dm','coef',ncoef )

  CALL zfftm1di ( n, coef )

  IF ( tin(1:1) == 'N' .OR. tin(1:1) == 'n' ) THEN

    IF ( tout(1:1) == 'N' .OR. tout(1:1) == 'n' ) THEN
      
      IF ( scale /= 1._dbl ) THEN
        zout(1:n,1:m) = scale * zin(1:n,1:m)
      ELSE
        zout(1:n,1:m) = zin(1:n,1:m)
      END IF
      CALL zfftm1d ( -fsign, n, m, zout, 1, ldi, coef )

    ELSE IF ( tout(1:1) == 'T' .OR. tout(1:1) == 't' ) THEN

      ALLOCATE ( zz ( ldi, lmi ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'fft2dm','zz',ldi*lmi )
      IF ( scale /= 1._dbl ) THEN
        zz(1:n,1:m) = scale * zin(1:n,1:m)
      ELSE
        zz(1:n,1:m) = zin(1:n,1:m)
      END IF
      CALL zfftm1d ( -fsign, n, m, zz, 1, ldi, coef )
      zout(1:m,1:n) = TRANSPOSE ( zz(1:n,1:m) )
      DEALLOCATE ( zz, STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'fft2dm','zz' )

    ELSE

      CALL stop_program ( " fft1dm ", "illegal tout" )

    END IF

  ELSE IF ( tin(1:1) == 'T' .OR. tin(1:1) == 't' ) THEN

    IF ( tout(1:1) == 'N' .OR. tout(1:1) == 'n' ) THEN

      IF ( scale /= 1._dbl ) THEN
        zout(1:n,1:m) = scale * TRANSPOSE ( zin(1:m,1:n) )
      ELSE
        zout(1:n,1:m) = TRANSPOSE ( zin(1:m,1:n) )
      END IF
      CALL zfftm1d ( -fsign, n, m, zout, 1, ldo, coef )

    ELSE IF ( tout(1:1) == 'T' .OR. tout(1:1) == 't' ) THEN

      ALLOCATE ( zz ( lmi, ldi ), STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'fft2dm','zz',ldi*lmi )
      IF ( scale /= 1._dbl ) THEN
        zz(1:n,1:m) = scale * TRANSPOSE ( zin(1:m,1:n) )
      ELSE
        zz(1:n,1:m) = TRANSPOSE ( zin(1:m,1:n) )
      END IF
      CALL zfftm1d ( -fsign, n, m, zz, 1, lmi, coef )
      zout(1:m,1:n) = TRANSPOSE ( zz(1:n,1:m) )
      DEALLOCATE ( zz, STAT = isos )
      IF ( isos /= 0 ) CALL stop_memory ( 'fft2dm','zz' )

    ELSE

      CALL stop_program ( " fft1dm ", "illegal tout" )

    END IF

  ELSE

    CALL stop_program ( " fft1dm ", "illegal tin" )

  END IF

  DEALLOCATE ( coef, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fft1dm','coef' )

#else

  fsign = 0

#endif

END SUBROUTINE mltfft

!!*****
!******************************************************************************

END MODULE fftsgi_lib

!******************************************************************************
