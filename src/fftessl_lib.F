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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE fftessl_lib
  
  USE kinds, ONLY: dbl
  USE stop_program, ONLY : stop_memory, stop_prg
  USE util, ONLY : sort
  
  IMPLICIT NONE 
  
  PRIVATE
  PUBLIC :: fft3d
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

SUBROUTINE fft_get_lengths ( data, max_length )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: max_length
  INTEGER, DIMENSION ( : ), POINTER :: data

! Locals
  INTEGER :: iloc
  INTEGER, DIMENSION ( : ), ALLOCATABLE :: idx
  INTEGER :: h, i, j, k, m, number, ndata, nmax, allocstat, maxn
  INTEGER :: maxn_twos, maxn_threes, maxn_fives
  INTEGER :: maxn_sevens, maxn_elevens, maxn_thirteens

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
!!****** fftessl_lib/fft_3d [1.0] *
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
  INTEGER :: nx, ny, nz, na1, na2, isos, naux
  COMPLEX ( dbl ), DIMENSION(:), ALLOCATABLE :: aux
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
  
  IF( MAX ( ny, nz ) < 252 ) THEN
     IF( nx <= 2048 ) THEN
        naux = 60000
     ELSE
        naux = 60000 + NINT( 4.56_dbl * nx )
     END IF
  ELSE
     IF( nx <= 2048 ) THEN
        na1 = 60000 + ( 2 * ny + 256 ) * ( MIN ( 64, nx ) + 4.56 )
        na2 = 60000 + ( 2 * nz + 256 ) * ( MIN ( 64, nx * ny ) + 4.56 )
     ELSE
        na1 = 60000 + 4.56 * nx + ( 2 * ny + 256 ) * ( MIN ( 64, nx ) + 4.56 )
        naux = 60000 + 4.56 * nx + ( 2 * nz + 256 ) * ( MIN ( 64, nx * ny ) + 4.56 )
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

!    sign_fft = fsign

  IF ( fft_in_place ) THEN

    CALL dcft3(zg,ldx,ldx*ldy,zg,ldx,ldx*ldy,nx,ny,nz,sign_fft,scale,aux,naux)

  ELSE

    ldox = SIZE ( zg_out (:,1,1) )
    ldoy = SIZE ( zg_out (1,:,1) )
    ldoz = SIZE ( zg_out (1,1,:) )

    CALL dcft3(zg,ldx,ldx*ldy,zg_out,ldox,ldox*ldoz,nx,ny,nz,sign_fft,scale,aux,naux)

  END IF

  DEALLOCATE ( aux, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'fft3d','aux' )

#else

  fsign = 0

#endif

END SUBROUTINE fft3d

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

SUBROUTINE mltfft ( fsign, trans, scale, n, m, zin, zout )
  
  IMPLICIT NONE
  
! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: trans
  INTEGER, INTENT ( INOUT ) :: fsign
  REAL ( dbl ), INTENT ( IN ) :: scale
  INTEGER, INTENT ( IN ) :: n, m
  COMPLEX ( dbl ), DIMENSION(:,:), INTENT ( IN ) :: zin
  COMPLEX ( dbl ), DIMENSION(:,:), INTENT ( OUT ) :: zout
  
! Locals
  INTEGER :: ldi, lmi, ldo, lmo
  
!------------------------------------------------------------------------------
  
  ldi = SIZE ( zin (:,1) )
  lmi = SIZE ( zin (1,:) )
  ldo = SIZE ( zin (:,1) )
  lmo = SIZE ( zin (1,:) )
  
#if defined ( __FFTESSL )

  stop " fftessl_lib "

#else

  fsign = 0

#endif

END SUBROUTINE mltfft

!!*****
!******************************************************************************

END MODULE fftessl_lib

!******************************************************************************
