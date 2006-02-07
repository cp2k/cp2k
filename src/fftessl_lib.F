!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fftessl_lib [1.0] *
!!
!!   NAME
!!     fftessl_lib
!!
!!   FUNCTION
!!     Interface to FFT routines in ESSL library
!!
!!   AUTHOR
!!     JGH 30-Nov-2000
!!
!!   MODIFICATION HISTORY
!!     JGH 20-Feb-2001 : mltfft, mfft2d routines added
!!
!!   SOURCE
!******************************************************************************

MODULE fftessl_lib

  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory
  USE util,                            ONLY: sort

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fft3d, mltfft, mfft2d
  PUBLIC :: fft_get_lengths

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** fftessl_lib/fft_get_lengths [1.0] *
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

SUBROUTINE fft_get_lengths ( DATA, max_length )


    INTEGER, DIMENSION(:), POINTER           :: DATA
    INTEGER, INTENT(IN)                      :: max_length

    INTEGER :: allocstat, h, i, j, k, m, maxn, maxn_elevens, maxn_fives, &
      maxn_sevens, maxn_thirteens, maxn_threes, maxn_twos, ndata, nmax, number
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: idx

!------------------------------------------------------------------------------
! compute ndata

  maxn_twos = 15
  maxn_threes = 2
  maxn_fives = 1
  maxn_sevens = 1
  maxn_elevens = 1
  maxn_thirteens = 0
  maxn = MIN ( max_length, 37748736 )

  ndata = 0
  DO h = 1, maxn_twos
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

  ALLOCATE ( DATA ( ndata ), idx ( ndata ), STAT = allocstat )
  IF ( allocstat /= 0 ) THEN
     CALL stop_memory ( "fft_get_lengths", "data, idx", 2*ndata )
  END IF

  ndata = 0
  DATA ( : ) = 0
  DO h = 1, maxn_twos
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
                 DATA ( ndata ) = number
              END DO
           END DO
        END DO
     END DO
  END DO

  CALL sort ( DATA, ndata, idx )

  DEALLOCATE ( idx, STAT = allocstat )
  IF ( allocstat /= 0 ) THEN
     CALL stop_memory ( "fft_get_lengths", "idx" )
  END IF

END SUBROUTINE fft_get_lengths

!!*****
!******************************************************************************
!!****** fftessl_lib/fft3d [1.0] *
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
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: scale
    INTEGER, DIMENSION(:), INTENT(IN)        :: n
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: zg
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL                :: zg_out

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: aux
    INTEGER                                  :: isos, ldox, ldoy, ldoz, ldx, &
                                                ldy, ldz, na1, na2, naux, nx, &
                                                ny, nz, sign_fft
    LOGICAL                                  :: fft_in_place

!------------------------------------------------------------------------------

  IF ( PRESENT ( zg_out ) ) THEN
     fft_in_place = .FALSE.
  ELSE
     fft_in_place = .TRUE.
  END IF

  sign_fft = fsign

  nx = n ( 1 )
  ny = n ( 2 )
  nz = n ( 3 )

  ldx = SIZE ( zg, 1 )
  ldy = SIZE ( zg, 2 )
  ldz = SIZE ( zg, 3 )
  
  IF( MAX ( ny, nz ) < 252 ) THEN
     IF( nx <= 2048 ) THEN
        naux = 60000
     ELSE
        naux = 60000 + NINT( 4.56_dp * nx )
     END IF
  ELSE
     IF( nx <= 2048 ) THEN
        na1 = 60000 + ( 2 * ny + 256 ) * ( MIN ( 64, nx ) + 4.56 )
        na2 = 60000 + ( 2 * nz + 256 ) * ( MIN ( 64, nx * ny ) + 4.56 )
     ELSE
        na1 = 60000 + 4.56 * nx + ( 2 * ny + 256 ) * ( MIN ( 64, nx ) + 4.56 )
        na2 = 60000 + 4.56 * nx + ( 2 * nz + 256 ) * ( MIN ( 64, nx * ny ) + 4.56 )
     END IF
     IF( ny >= 252 .AND. nz < 252 ) THEN
        naux=na1
     ELSE IF( ny < 252 .AND. nz >= 252 ) THEN
        naux = na2
     ELSE
        naux = MAX ( na1, na2 )
     END IF
  END IF

#if defined ( __FFTESSL )

  ALLOCATE ( aux ( naux ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fft3d','aux',naux )

  IF ( fft_in_place ) THEN

    CALL dcft3(zg,ldx,ldx*ldy,zg,ldx,ldx*ldy,nx,ny,nz,sign_fft,scale,aux,naux)

  ELSE

    ldox = SIZE ( zg_out, 1 )
    ldoy = SIZE ( zg_out, 2 )
    ldoz = SIZE ( zg_out, 3 )

    CALL dcft3(zg,ldx,ldx*ldy,zg_out,ldox,ldox*ldoy,nx,ny,nz,sign_fft,scale,aux,naux)

  END IF

  DEALLOCATE ( aux, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fft3d','aux' )

#else

  fsign = 0

#endif

END SUBROUTINE fft3d

!!*****
!******************************************************************************
!!****** fftessl_lib/mfft2d [1.0] *
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
    REAL(KIND=dp), INTENT(IN)                :: scale
    INTEGER, INTENT(IN)                      :: n1, n2, nm
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: zin
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: zout

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: aux1, aux2
    INTEGER :: i, i1, i2, i3, ierr, inc1x, inc1y, inc2x, inc2y, init, li1, &
      li2, li3, lo1, lo2, lo3, naux1, naux2, o1, o2, o3, sign_fft
    REAL(KIND=dp)                            :: r, s

!------------------------------------------------------------------------------

  sign_fft = fsign

  li1 = SIZE ( zin, 1 )
  li2 = SIZE ( zin, 2 )
  li3 = SIZE ( zin, 3 )

  lo1 = SIZE ( zout, 1 )
  lo2 = SIZE ( zout, 2 )
  lo3 = SIZE ( zout, 3 )

  IF( MAX ( n1, n2 ) < 2048 ) THEN
     naux1 = 40000
  ELSE
     naux1 = 40000 + NINT ( 2.28_dp * REAL ( n1 + n2,KIND=dp) )
  END IF
  IF( MAX ( n1, n2 ) < 252 ) THEN
     naux2 = 20000
  ELSE
     r = REAL ( MAX ( n1, n2 ),KIND=dp)
     s = REAL ( MIN ( 64, n1, n2 ),KIND=dp)
     naux2 = 20000 + NINT ( ( 2.0_dp * r + 256.0_dp ) * ( s + 2.28_dp ) )
  END IF

  IF ( rin(1:1) == "T" .OR. rin(1:1) == "t" ) THEN
    inc1x = li1
    inc2x = li1*li2
  ELSE
    inc1x = 1
    inc2x = li1
  END IF
  IF ( rout(1:1) == "T" .OR. rout(1:1) == "t" ) THEN
    inc1y = lo1
    inc2y = lo1*lo2
  ELSE
    inc1y = 1
    inc2y = lo1
  END IF

#if defined ( __FFTESSL )

  ALLOCATE ( aux1 ( naux1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'mfft2d','aux1',naux1 )
  ALLOCATE ( aux2 ( naux2 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'mfft2d','aux1',naux2 )

  init = 1
  CALL dcft2(init,zin,inc1x,inc2x,zout,inc1y,inc2y, &
             n1,n2,sign_fft,scale,aux1,naux1,aux2,naux2)
  init = 0
  DO i = 1, nm
    IF ( rin(1:1) == "T" .OR. rin(1:1) == "t" ) THEN
      i1 = i
      i2 = 1
      i3 = 1
    ELSE
      i1 = 1
      i2 = 1
      i3 = i
    END IF
    IF ( rout(1:1) == "T" .OR. rout(1:1) == "t" ) THEN
      o1 = i
      o2 = 1
      o3 = 1
    ELSE
      o1 = 1
      o2 = 1
      o3 = i
    END IF
    CALL dcft2(init,zin(i1,i2,i3),inc1x,inc2x,zout(o1,o2,o3),inc1y,inc2y, &
               n1,n2,sign_fft,scale,aux1,naux1,aux2,naux2)
  END DO

  DEALLOCATE ( aux1, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'mfft2d','aux1' )
  DEALLOCATE ( aux2, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'mfft2d','aux2' )

#else

  fsign = 0

#endif

END SUBROUTINE mfft2d

!!*****
!******************************************************************************
!!****** fftessl_lib/mltfft [1.0] *
!!
!!   NAME
!!     mltfft
!!
!!   FUNCTION
!!     Calls multiple 1d FFT from the ESSL library
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
    REAL(KIND=dp), INTENT(IN)                :: scale
    INTEGER, INTENT(IN)                      :: n, m
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: zin
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: zout

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: aux1, aux2
    INTEGER                                  :: ierr, inc1x, inc1y, inc2x, &
                                                inc2y, init, isign, ldi, ldo, &
                                                lmi, lmo, naux1, naux2

!------------------------------------------------------------------------------

  ldi = SIZE ( zin, 1 )
  lmi = SIZE ( zin, 2 )
  ldo = SIZE ( zout, 1 )
  lmo = SIZE ( zout, 2 )
  
  IF( n < 2048 ) THEN
     naux1 = 20000
  ELSE
     naux1 = 20000 + NINT ( 2.28_dp * REAL ( n,KIND=dp) )
  END IF
  IF ( tin(1:1) == "T" .OR. tin(1:1) == "t" .OR. &
       tout(1:1) == "T" .OR. tout(1:1) == "t" ) THEN
     naux1 = naux1 + ( 2 * n + 256 ) * MIN ( 64, n )
  END IF
  naux2 = naux1

  IF ( tin(1:1) == "T" .OR. tin(1:1) == "t" ) THEN
    inc1x = ldi
    inc2x = 1
  ELSE
    inc1x = 1
    inc2x = ldi
  END IF
  IF ( tout(1:1) == "T" .OR. tout(1:1) == "t" ) THEN
    inc1y = ldo
    inc2y = 1
  ELSE
    inc1y = 1
    inc2y = ldo
  END IF

  isign = fsign

#if defined ( __FFTESSL )

  ALLOCATE ( aux1 ( naux1 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'mltfft','aux1',naux1 )
  ALLOCATE ( aux2 ( naux2 ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'mltfft','aux1',naux2 )

  init = 1
  CALL DCFT ( init, zin, inc1x, inc2x, zout, inc1y, inc2y, &
              n, m, isign, scale, aux1, naux1, aux2, naux2 )
  init = 0
  CALL DCFT ( init, zin, inc1x, inc2x, zout, inc1y, inc2y, &
              n, m, isign, scale, aux1, naux1, aux2, naux2 )

  DEALLOCATE ( aux1, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'mltfft','aux1' )
  DEALLOCATE ( aux2, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( 'mltfft','aux2' )

#else

  fsign = 0

#endif

END SUBROUTINE mltfft

!!*****
!******************************************************************************

END MODULE fftessl_lib

!******************************************************************************
