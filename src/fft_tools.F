!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE fft_tools
  
  USE kinds, ONLY: dbl, sgl
!TMPTMPTMP  USE timings, ONLY: timeset, timestop
  USE stop_program, ONLY : stop_prg
  USE util, ONLY:  sort
  
  IMPLICIT NONE 
  
  PRIVATE
  PUBLIC :: fft_wrap, FWFFT, BWFFT
  PUBLIC :: FFT_RADIX_CLOSEST, FFT_RADIX_NEXT, FFT_RADIX_ALLOWED
  PUBLIC :: FFT_RADIX_DISALLOWED
  PUBLIC :: fft_radix_operations
  
  INTEGER, PARAMETER :: FWFFT = +1, BWFFT = -1
  INTEGER, PARAMETER :: FFT_RADIX_CLOSEST = 493, FFT_RADIX_NEXT = 494
  INTEGER, PARAMETER :: FFT_RADIX_ALLOWED = 495, FFT_RADIX_DISALLOWED = 496
  
  INTERFACE fft_wrap
     MODULE PROCEDURE fft_wrap_dbl, fft_wrap_sgl
  END INTERFACE
  
#if defined ( __FFTW )
!*apsi 220500 This is actually the file fftw_f77.i ...
!     This file contains PARAMETER statements for various constants
!     that can be passed to FFTW routines.  You should include
!     this file in any FORTRAN program that calls the fftw_f77
!     routines (either directly or with an #include statement
!     if you use the C preprocessor).

      integer FFTW_FORWARD,FFTW_BACKWARD
      parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

      integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
      parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

      integer FFTW_ESTIMATE,FFTW_MEASURE
      parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)

      integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
      parameter (FFTW_OUT_OF_PLACE=0)
      parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

      integer FFTW_THREADSAFE
      parameter (FFTW_THREADSAFE=128)

!     Constants for the MPI wrappers:
      integer FFTW_TRANSPOSED_ORDER, FFTW_NORMAL_ORDER
      integer FFTW_SCRAMBLED_INPUT, FFTW_SCRAMBLED_OUTPUT
      parameter(FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0)
      parameter(FFTW_SCRAMBLED_INPUT=8192)
      parameter(FFTW_SCRAMBLED_OUTPUT=16384)
#endif

CONTAINS
  
!******************************************************************************

!! Give the allowed lengths of FFT's   '''

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
!! FFTW can do arbitrary(?) lenghts, maybe you want to limit them to some
!!    powers of small prime numbers though...
  
#if defined ( __FFTW )
  maxn_twos = 15
  maxn_threes = 3
  maxn_fives = 2
  maxn_sevens = 1
  maxn_elevens = 1
  maxn_thirteens = 0
  maxn = MIN ( max_length, 37748736 )
  
#elif defined ( __AIX )
  maxn_twos = 15
  maxn_threes = 2
  maxn_fives = 1
  maxn_sevens = 1
  maxn_elevens = 1
  maxn_thirteens = 0
  maxn = MIN ( max_length, 37748736 )
  
#else
  CALL stop_prg ( "fft_get_lengths", "no fft basis defined" )
#endif
  
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
     CALL stop_prg ( "fft_get_lengths", "error allocating data, idx" )
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
     CALL stop_prg ( "fft_get_lengths", "error deallocating idx" )
  END IF
  
END SUBROUTINE fft_get_lengths

!******************************************************************************

SUBROUTINE fft_radix_operations ( radix_in, radix_out, operation )
!! Determine the allowed lengths of FFT's   '''
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: radix_in, operation
  INTEGER, INTENT ( OUT ) :: radix_out
  
! Locals
  INTEGER :: i, iloc
  INTEGER, DIMENSION ( : ), POINTER :: data
  
!------------------------------------------------------------------------------
  
  CALL fft_get_lengths ( data, max_length = 1024 )
  
  iloc = 0
  DO i = 1, SIZE ( data )
     IF ( data ( i ) == radix_in ) THEN
        iloc = i
        EXIT
     ELSE
        IF ( OPERATION == FFT_RADIX_ALLOWED ) THEN
           CYCLE
        ELSE IF ( data ( i ) > radix_in ) THEN
           iloc = i
           EXIT
        END IF
     END IF
  END DO
  
  IF ( iloc == 0 ) THEN
     CALL stop_prg ( "fft_radix_operations", "index to radix array not found" )
  END IF
  
  IF ( OPERATION == FFT_RADIX_ALLOWED ) THEN
     IF ( data ( iloc ) == radix_in ) THEN
        radix_out = FFT_RADIX_ALLOWED
     ELSE
        radix_out = FFT_RADIX_DISALLOWED
     END IF
     
  ELSE IF ( OPERATION == FFT_RADIX_CLOSEST ) THEN
     IF ( data ( iloc ) == radix_in ) THEN
        radix_out = data ( iloc )
     ELSE
        IF ( ABS ( data ( iloc - 1 ) - radix_in ) <= &
             ABS ( data ( iloc ) - radix_in ) ) THEN
           radix_out = data ( iloc - 1 )
        ELSE
           radix_out = data ( iloc )
        END IF
     END IF
     
  ELSE IF ( OPERATION == FFT_RADIX_NEXT ) THEN
     radix_out = data ( iloc )
     
  ELSE
     CALL stop_prg ( "fft_radix_operations", "disallowed radix operation" )
  END IF
  
END SUBROUTINE fft_radix_operations

!******************************************************************************

#if defined ( __FFTSG )

SUBROUTINE fft_wrap_sgl ( fsign, n, zg, zgout, scale )
! routine with wrapper for all fft calls:
! Does transform with exp(+ig.r*sign):
  IMPLICIT NONE
  INTEGER, INTENT ( IN ):: fsign
  INTEGER, INTENT ( IN ), DIMENSION ( : ) :: n
  COMPLEX(sgl), DIMENSION(0:,0:,0:), INTENT ( INOUT ):: zg
  COMPLEX(sgl), DIMENSION(0:,0:,0:), OPTIONAL, INTENT ( INOUT ):: zgout
  REAL(sgl), OPTIONAL, INTENT ( IN ):: scale
! locals

  stop "SG sFFT not implemented"

END SUBROUTINE fft_wrap_sgl

!******************************************************************************

SUBROUTINE fft_wrap_dbl ( fsign, n, zg, zgout, scale )
! routine with wrapper for all fft calls:
! Does transform with exp(+ig.r*sign):
  IMPLICIT NONE
  INTEGER, INTENT ( IN ):: fsign
  INTEGER, INTENT ( IN ), DIMENSION ( : ) :: n
  COMPLEX ( dbl ), DIMENSION(0:,0:,0:), INTENT ( INOUT ):: zg
  COMPLEX ( dbl ), DIMENSION(0:,0:,0:), OPTIONAL, INTENT ( INOUT ):: zgout
  REAL ( dbl ), OPTIONAL, INTENT ( IN ):: scale
! locals
  INTEGER :: ldx, ldy, ldz, nx, ny, nz, isign, norm
  COMPLEX ( dbl ), DIMENSION(SIZE(zg)) :: xf, yf

  ldx = n(1)
  ldy = n(2)
  ldz = n(3)
  nx = n(1)
  ny = n(1)
  nz = n(1)
  norm = 1._dbl

!  CALL mltfft( 'N','T',zg,ldx,ldy*ldz,xf,ldy*ldz,ldx,nx, &
!       ldy*ldz,isign,1.d0)
!  CALL mltfft( 'N','T',xf,ldy,ldx*ldz,yf,ldx*ldz,ldy,ny, &
!       ldx*ldz,isign,1.d0)
!  IF ( PRESENT ( zgout ) ) THEN
!     CALL mltfft( 'N','T',yf,ldz,ldy*ldx,zgout,ldy*ldx,ldz,nz, &
!          ldy*ldx,isign,scale)
!  ELSE
!     CALL mltfft( 'N','T',yf,ldz,ldy*ldx,zg,ldy*ldx,ldz,nz, &
!          ldy*ldx,isign,scale)
!  END IF

END SUBROUTINE fft_wrap_dbl

!******************************************************************************

#elif defined ( __FFTW ) && defined ( __AIX ) && defined ( __SPRECFFT )

!*** Single precision from FFTW ...

SUBROUTINE fft_wrap_sgl ( fsign, n, zg, zgout, scale )
! routine with wrapper for all fft calls:
! Does transform with exp(+ig.r*sign):
  IMPLICIT NONE
  INTEGER, INTENT ( IN ):: fsign
  INTEGER, INTENT ( IN ), DIMENSION ( : ) :: n
  COMPLEX(sgl), DIMENSION(0:,0:,0:), INTENT ( INOUT ):: zg
  COMPLEX(sgl), DIMENSION(0:,0:,0:), OPTIONAL, INTENT ( INOUT ):: zgout
  REAL(sgl), OPTIONAL, INTENT ( IN ):: scale
! locals
  INTEGER, SAVE :: n1a_save = -1, n2a_save = -1, n3a_save = -1
  INTEGER, SAVE :: n1b_save = -1, n2b_save = -1, n3b_save = -1
  LOGICAL, SAVE :: ffta_in_place = .true., fftb_in_place = .true.
  LOGICAL :: fft_in_place
  INTEGER :: sign_fft, n1, n2, n3
  INTEGER ( KIND = 4 ), SAVE :: plan_a_fw, plan_a_bw
  INTEGER ( KIND = 4 ), SAVE :: plan_b_fw, plan_b_bw
  REAL(sgl) :: norm

  IF( PRESENT ( scale) ) THEN
     norm=scale
  ELSE
     norm=1._dbl
  END IF

  n1 = n(1)
  n2 = n(2)
  n3 = n(3)

  IF ( PRESENT ( zgout ) ) THEN
     fft_in_place = .false.
  ELSE
     fft_in_place = .true.
  END IF

  sign_fft = fsign

  IF ( n1a_save == n1 .AND. n2a_save == n2 .AND. n3a_save == n3 &
       .AND. ( fft_in_place .eqv. ffta_in_place ) ) THEN
     IF ( sign_fft == +1 ) THEN
        CALL fftwnd_f77_one ( plan_a_fw, zg, zgout )
     ELSE
        CALL fftwnd_f77_one ( plan_a_bw, zg, zgout )
     END IF
  ELSE IF ( n1b_save == n1 .AND. n2b_save == n2 .AND. n3b_save == n3 &
       .AND. ( fft_in_place .eqv. fftb_in_place ) ) THEN
     IF ( sign_fft == +1 ) THEN
        CALL fftwnd_f77_one ( plan_b_fw, zg, zgout )
     ELSE
        CALL fftwnd_f77_one ( plan_b_bw, zg, zgout )
     END IF
  ELSE IF ( n1a_save == -1 .OR. &
       ( n1a_save == n1 .AND. n2a_save == n2 .AND. n3a_save == n3 &
       .AND. ( fft_in_place .neqv. ffta_in_place ) .AND. &
       n1b_save /= -1 ) ) THEN ! Initialise 'a'
     write(6,*) "Init a"
     IF ( fft_in_place ) THEN
        CALL fftw3d_f77_create_plan ( plan_a_fw, n1, n2, n3, FFTW_FORWARD, &
             FFTW_ESTIMATE + FFTW_IN_PLACE )
        CALL fftw3d_f77_create_plan ( plan_a_bw, n1, n2, n3, FFTW_BACKWARD, &
             FFTW_ESTIMATE + FFTW_IN_PLACE )
     ELSE
        CALL fftw3d_f77_create_plan ( plan_a_fw, n1, n2, n3, FFTW_FORWARD, &
             FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
        CALL fftw3d_f77_create_plan ( plan_a_bw, n1, n2, n3, FFTW_BACKWARD, &
             FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
     END IF
     n1a_save = n1
     n2a_save = n2
     n3a_save = n3
     ffta_in_place = fft_in_place

     IF ( sign_fft == +1 ) THEN
        CALL fftwnd_f77_one ( plan_a_fw, zg, zgout )
     ELSE
        CALL fftwnd_f77_one ( plan_a_bw, zg, zgout )
     END IF
  ELSE ! Initialise 'b'
     write(6,*) "Init b"
     IF ( fft_in_place ) THEN
        CALL fftw3d_f77_create_plan ( plan_b_fw, n1, n2, n3, FFTW_FORWARD, &
             FFTW_ESTIMATE + FFTW_IN_PLACE )
        CALL fftw3d_f77_create_plan ( plan_b_bw, n1, n2, n3, FFTW_BACKWARD, &
             FFTW_ESTIMATE + FFTW_IN_PLACE )
     ELSE
        CALL fftw3d_f77_create_plan ( plan_b_fw, n1, n2, n3, FFTW_FORWARD, &
             FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
        CALL fftw3d_f77_create_plan ( plan_b_bw, n1, n2, n3, FFTW_BACKWARD, &
             FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
     END IF
     n1b_save = n1
     n2b_save = n2
     n3b_save = n3
     fftb_in_place = fft_in_place

     IF ( sign_fft == +1 ) THEN
        CALL fftwnd_f77_one ( plan_b_fw, zg, zgout )
     ELSE
        CALL fftwnd_f77_one ( plan_b_bw, zg, zgout )
     END IF
  END IF

  IF ( norm /= 1._dbl ) THEN
     IF ( fft_in_place ) THEN
        zg = zg * norm
     ELSE
        zgout = zgout * norm
     END IF
  END IF

END SUBROUTINE fft_wrap_sgl

!*** ... and double precision from ESSL

SUBROUTINE fft_wrap_dbl(fsign,n,zg,zgout,scale)
! routine with wrapper for all fft calls:
! Does transform with exp(+ig.r*sign):
  IMPLICIT NONE
  INTEGER, INTENT ( IN ):: fsign
  INTEGER, INTENT ( IN ), DIMENSION ( : ) :: n
  COMPLEX ( dbl ), DIMENSION(0:,0:,0:), INTENT ( INOUT ):: zg
  COMPLEX ( dbl ), DIMENSION(0:,0:,0:), OPTIONAL, INTENT ( INOUT ):: zgout
  REAL ( dbl ), OPTIONAL, INTENT ( IN ):: scale
! locals
!    REAL(kind=dbl), DIMENSION ( : ), ALLOCATABLE:: aux
  REAL ( dbl ):: norm
  INTEGER:: la1,la2,la3,sign_fft,na1,na2,naux,isos
  INTEGER :: n1,n2,n3

  IF( PRESENT ( scale) ) THEN
     norm=scale
  ELSE
     norm=1._dbl
  END IF
! IBM FFT CALL
  n1 = n(1)
  n2 = n(2)
  n3 = n(3)
  la1 = n1
  la2 = n1*n2
  IF( MAX(n2,n3)<252 ) THEN
     IF( n1<=2048 ) THEN
        naux=60000
     ELSE
        naux=60000+4.56*n1
     END IF
  ELSE
     IF( n1<=2048 ) THEN
        na1=60000+(2*n2+256)*(MIN(64,n1)+4.56)
        na2=60000+(2*n3+256)*(MIN(64,n1*n2)+4.56)
     ELSE
        na1=60000+4.56*n1+(2*n2+256)*(MIN(64,n1)+4.56)
        naux=60000+4.56*n1+(2*n3+256)*(MIN(64,n1*n2)+4.56) 
     END IF
     IF( n2>=252 .AND. n3<252 ) THEN
        naux=na1
     ELSE IF( n2<252 .AND. n3>=252 ) THEN
        naux=na2
     ELSE
        naux=MAX(na1,na2)
     END IF
  END IF

!    ALLOCATE(aux(naux),STAT=isos)
!    IF ( isos /= 0 ) CALL stop_prg( 'fft_wrap','failed to allocate aux')

!    sign_fft = fsign
!    CALL dcft3(zg,la1,la2,zgout,la1,la2,n1,n2,n3,sign_fft,norm,aux,naux)

  IF ( PRESENT ( zgout ) ) THEN
     CALL do_ibm_fft_ab ( naux )
  ELSE
     CALL do_ibm_fft_aa ( naux )
  END IF

!    DEALLOCATE(aux,STAT=isos)
!    IF ( isos /= 0 ) CALL stop_prg( 'fft_wrap','failed to deallocate aux')

CONTAINS

  SUBROUTINE do_ibm_fft_ab ( naux )

    INTEGER, INTENT ( IN ) :: naux

    REAL(kind=dbl), DIMENSION ( naux ) :: aux

    sign_fft = fsign
    CALL dcft3(zg,la1,la2,zgout,la1,la2,n1,n2,n3,sign_fft,norm,aux,naux)

  END SUBROUTINE do_ibm_fft_ab

  SUBROUTINE do_ibm_fft_aa ( naux )

    INTEGER, INTENT ( IN ) :: naux

    REAL(kind=dbl), DIMENSION ( naux ) :: aux

    sign_fft = fsign
    CALL dcft3(zg,la1,la2,zg,la1,la2,n1,n2,n3,sign_fft,norm,aux,naux)

  END SUBROUTINE do_ibm_fft_aa

END SUBROUTINE fft_wrap_dbl

!******************************************************************************

#elif defined ( __FFTW ) && ! defined ( __SPRECFFT )

SUBROUTINE fft_wrap_sgl ( fsign, n, zg, zgout, scale )
! routine with wrapper for all fft calls:
! Does transform with exp(+ig.r*sign):
  IMPLICIT NONE
  INTEGER, INTENT ( IN ):: fsign
  INTEGER, INTENT ( IN ), DIMENSION ( : ) :: n
  COMPLEX(sgl), DIMENSION(0:,0:,0:), INTENT ( INOUT ):: zg
  COMPLEX(sgl), DIMENSION(0:,0:,0:), OPTIONAL, INTENT ( INOUT ):: zgout
  REAL(sgl), OPTIONAL, INTENT ( IN ):: scale

  stop "fftw cannot be used in both single and double precision..."
END SUBROUTINE fft_wrap_sgl

!******************************************************************************

! routine with wrapper for all fft calls:
! Does transform with exp(+ig.r*sign):

SUBROUTINE fft_wrap_dbl ( fsign, n, zg, zg_out, scale )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: fsign
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: n
  COMPLEX ( dbl ), DIMENSION(:,:,:), INTENT ( INOUT ) :: zg
  COMPLEX ( dbl ), DIMENSION(:,:,:), INTENT ( INOUT ), OPTIONAL, TARGET :: zg_out
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: scale
  
! Locals
  INTEGER, SAVE :: n1a_save = -1, n2a_save = -1, n3a_save = -1
  INTEGER, SAVE :: n1b_save = -1, n2b_save = -1, n3b_save = -1
  INTEGER, SAVE :: n1c_save = -1, n2c_save = -1, n3c_save = -1
  LOGICAL, SAVE :: ffta_in_place = .TRUE., fftb_in_place = .TRUE.
  LOGICAL, SAVE :: fftc_in_place = .TRUE.
  LOGICAL :: fft_in_place
  INTEGER :: sign_fft, n1, n2, n3
  INTEGER ( KIND = 8 ), SAVE :: plan_a_fw, plan_a_bw
  INTEGER ( KIND = 8 ), SAVE :: plan_b_fw, plan_b_bw
  INTEGER ( KIND = 8 ), SAVE :: plan_c_fw, plan_c_bw
  REAL ( dbl ) :: norm
  
  ! Just due to DEC... apsi
  COMPLEX ( dbl ), DIMENSION(:,:,:), POINTER :: zgout
  COMPLEX ( dbl ), DIMENSION(1,1,1), TARGET :: zgout_tmp
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( scale ) ) THEN
     norm = scale
  ELSE
     norm = 1._dbl
  END IF
  
  n1 = n(1)
  n2 = n(2)
  n3 = n(3)
  
  IF ( PRESENT ( zg_out ) ) THEN
     fft_in_place = .false.
     zgout => zg_out
  ELSE
     fft_in_place = .true.
     zgout => zgout_tmp
  END IF
  
  sign_fft = fsign
  
  IF ( n1a_save == n1 .AND. n2a_save == n2 .AND. n3a_save == n3 &
       .AND. ( fft_in_place .eqv. ffta_in_place ) ) THEN
     IF ( sign_fft == +1 ) THEN
        CALL fftwnd_f77_one ( plan_a_fw, zg, zgout )
     ELSE
        CALL fftwnd_f77_one ( plan_a_bw, zg, zgout )
     END IF
  ELSE IF ( n1b_save == n1 .AND. n2b_save == n2 .AND. n3b_save == n3 &
       .AND. ( fft_in_place .eqv. fftb_in_place ) ) THEN
     IF ( sign_fft == +1 ) THEN
        CALL fftwnd_f77_one ( plan_b_fw, zg, zgout )
     ELSE
        CALL fftwnd_f77_one ( plan_b_bw, zg, zgout )
     END IF
  ELSE IF ( n1c_save == n1 .AND. n2c_save == n2 .AND. n3c_save == n3 &
       .AND. ( fft_in_place .eqv. fftc_in_place ) ) THEN
     IF ( sign_fft == +1 ) THEN
        CALL fftwnd_f77_one ( plan_c_fw, zg, zgout )
     ELSE
        CALL fftwnd_f77_one ( plan_c_bw, zg, zgout )
     END IF
  ELSE IF ( n1a_save == -1 .OR. &
       ( n1a_save == n1 .AND. n2a_save == n2 .AND. n3a_save == n3 &
       .AND. ( fft_in_place .neqv. ffta_in_place ) .AND. &
       ( n1b_save /= -1 .AND. n1c_save /= -1 ) ) ) THEN ! Initialise 'a'
     write(6,*) "Init a"
     IF ( fft_in_place ) THEN
        CALL fftw3d_f77_create_plan ( plan_a_fw, n1, n2, n3, FFTW_FORWARD, &
             FFTW_ESTIMATE + FFTW_IN_PLACE )
        CALL fftw3d_f77_create_plan ( plan_a_bw, n1, n2, n3, FFTW_BACKWARD, &
             FFTW_ESTIMATE + FFTW_IN_PLACE )
     ELSE
        CALL fftw3d_f77_create_plan ( plan_a_fw, n1, n2, n3, FFTW_FORWARD, &
             FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
        CALL fftw3d_f77_create_plan ( plan_a_bw, n1, n2, n3, FFTW_BACKWARD, &
             FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
     END IF
     n1a_save = n1
     n2a_save = n2
     n3a_save = n3
     ffta_in_place = fft_in_place

     IF ( sign_fft == +1 ) THEN
        CALL fftwnd_f77_one ( plan_a_fw, zg, zgout )
     ELSE
        CALL fftwnd_f77_one ( plan_a_bw, zg, zgout )
     END IF
     
  ELSE IF ( n1b_save == -1 .OR. &
       ( n1b_save == n1 .AND. n2b_save == n2 .AND. n3b_save == n3 &
       .AND. ( fft_in_place .neqv. fftb_in_place ) .AND. &
       n1c_save /= -1 ) ) THEN ! Initialise 'a'
     write(6,*) "Init b"
     IF ( fft_in_place ) THEN
        CALL fftw3d_f77_create_plan ( plan_b_fw, n1, n2, n3, FFTW_FORWARD, &
             FFTW_ESTIMATE + FFTW_IN_PLACE )
        CALL fftw3d_f77_create_plan ( plan_b_bw, n1, n2, n3, FFTW_BACKWARD, &
             FFTW_ESTIMATE + FFTW_IN_PLACE )
     ELSE
        CALL fftw3d_f77_create_plan ( plan_b_fw, n1, n2, n3, FFTW_FORWARD, &
             FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
        CALL fftw3d_f77_create_plan ( plan_b_bw, n1, n2, n3, FFTW_BACKWARD, &
             FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
     END IF
     n1b_save = n1
     n2b_save = n2
     n3b_save = n3
     fftb_in_place = fft_in_place

     IF ( sign_fft == +1 ) THEN
        CALL fftwnd_f77_one ( plan_b_fw, zg, zgout )
     ELSE
        CALL fftwnd_f77_one ( plan_b_bw, zg, zgout )
     END IF

  ELSE ! Initialise 'c'
     write(6,*) "Init c"
     IF ( fft_in_place ) THEN
        CALL fftw3d_f77_create_plan ( plan_c_fw, n1, n2, n3, FFTW_FORWARD, &
             FFTW_ESTIMATE + FFTW_IN_PLACE )
        CALL fftw3d_f77_create_plan ( plan_c_bw, n1, n2, n3, FFTW_BACKWARD, &
             FFTW_ESTIMATE + FFTW_IN_PLACE )
     ELSE
        CALL fftw3d_f77_create_plan ( plan_c_fw, n1, n2, n3, FFTW_FORWARD, &
             FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
        CALL fftw3d_f77_create_plan ( plan_c_bw, n1, n2, n3, FFTW_BACKWARD, &
             FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
     END IF
     n1c_save = n1
     n2c_save = n2
     n3c_save = n3
     fftc_in_place = fft_in_place

     IF ( sign_fft == +1 ) THEN
        CALL fftwnd_f77_one ( plan_c_fw, zg, zgout )
     ELSE
        CALL fftwnd_f77_one ( plan_c_bw, zg, zgout )
     END IF
  END IF

  IF ( norm /= 1._dbl ) THEN
     IF ( fft_in_place ) THEN
        zg = zg * norm
     ELSE
        zgout = zgout * norm
     END IF
  END IF

END SUBROUTINE fft_wrap_dbl

!******************************************************************************

#elif defined ( __AIX )

SUBROUTINE fft_wrap_dbl(fsign,n,zg,zgout,scale)
! routine with wrapper for all fft calls:
! Does transform with exp(+ig.r*sign):
  IMPLICIT NONE
  INTEGER, INTENT ( IN ):: fsign
  INTEGER, INTENT ( IN ), DIMENSION ( : ) :: n
  COMPLEX ( dbl ), DIMENSION(0:,0:,0:), INTENT ( INOUT ):: zg
  COMPLEX ( dbl ), DIMENSION(0:,0:,0:), OPTIONAL, INTENT ( INOUT ):: zgout
  REAL ( dbl ), OPTIONAL, INTENT ( IN ):: scale
! locals
!    REAL(kind=dbl), DIMENSION ( : ), ALLOCATABLE:: aux
  REAL ( dbl ):: norm
  INTEGER:: la1,la2,la3,sign_fft,na1,na2,naux,isos
  INTEGER :: n1,n2,n3

  IF( PRESENT ( scale) ) THEN
     norm=scale
  ELSE
     norm=1._dbl
  END IF
! IBM FFT CALL
  n1 = n(1)
  n2 = n(2)
  n3 = n(3)
  la1 = n1
  la2 = n1*n2
  IF( MAX(n2,n3)<252 ) THEN
     IF( n1<=2048 ) THEN
        naux=60000
     ELSE
        naux=60000+4.56*n1
     END IF
  ELSE
     IF( n1<=2048 ) THEN
        na1=60000+(2*n2+256)*(MIN(64,n1)+4.56)
        na2=60000+(2*n3+256)*(MIN(64,n1*n2)+4.56)
     ELSE
        na1=60000+4.56*n1+(2*n2+256)*(MIN(64,n1)+4.56)
        naux=60000+4.56*n1+(2*n3+256)*(MIN(64,n1*n2)+4.56) 
     END IF
     IF( n2>=252 .AND. n3<252 ) THEN
        naux=na1
     ELSE IF( n2<252 .AND. n3>=252 ) THEN
        naux=na2
     ELSE
        naux=MAX(na1,na2)
     END IF
  END IF

!    ALLOCATE(aux(naux),STAT=isos)
!    IF ( isos /= 0 ) CALL stop_prg( 'fft_wrap','failed to allocate aux')

!    sign_fft = fsign
!    CALL dcft3(zg,la1,la2,zgout,la1,la2,n1,n2,n3,sign_fft,norm,aux,naux)

  IF ( PRESENT ( zgout ) ) THEN
     CALL do_ibm_fft_ab ( naux )
  ELSE
     CALL do_ibm_fft_aa ( naux )
  END IF

!    DEALLOCATE(aux,STAT=isos)
!    IF ( isos /= 0 ) CALL stop_prg( 'fft_wrap','failed to deallocate aux')

CONTAINS

  SUBROUTINE do_ibm_fft_ab ( naux )

    INTEGER, INTENT ( IN ) :: naux

    REAL(kind=dbl), DIMENSION ( naux ) :: aux

    sign_fft = fsign
    CALL dcft3(zg,la1,la2,zgout,la1,la2,n1,n2,n3,sign_fft,norm,aux,naux)

  END SUBROUTINE do_ibm_fft_ab

  SUBROUTINE do_ibm_fft_aa ( naux )

    INTEGER, INTENT ( IN ) :: naux

    REAL(kind=dbl), DIMENSION ( naux ) :: aux

    sign_fft = fsign
    CALL dcft3(zg,la1,la2,zg,la1,la2,n1,n2,n3,sign_fft,norm,aux,naux)

  END SUBROUTINE do_ibm_fft_aa

END SUBROUTINE fft_wrap_dbl

!******************************************************************************

SUBROUTINE fft_wrap_sgl(fsign,n,zg,zgout,scale)
! routine with wrapper for all fft calls:
! Does transform with exp(+ig.r*sign):
  IMPLICIT NONE
  INTEGER, INTENT ( IN ):: fsign
  INTEGER, INTENT ( IN ), DIMENSION ( : ) :: n
  COMPLEX(sgl), DIMENSION(0:,0:,0:), INTENT ( INOUT ):: zg
  COMPLEX(sgl), DIMENSION(0:,0:,0:), OPTIONAL, INTENT ( INOUT ):: zgout
  REAL(sgl), OPTIONAL, INTENT ( IN ):: scale
! locals
!    REAL(kind=dbl), DIMENSION ( : ), ALLOCATABLE:: aux
  REAL(sgl):: norm
  INTEGER:: la1,la2,la3,sign_fft,na1,na2,naux,isos
  INTEGER :: n1,n2,n3

  IF( PRESENT ( scale) ) THEN
     norm=scale
  ELSE
     norm=1._dbl
  END IF
! IBM FFT CALL
  n1 = n(1)
  n2 = n(2)
  n3 = n(3)
  la1 = n1
  la2 = n1*n2
  IF( MAX(n2,n3)<252 ) THEN
     IF( n1<=2048 ) THEN
        naux=60000
     ELSE
        naux=60000+4.56*n1
     END IF
  ELSE
     IF( n1<=2048 ) THEN
        na1=60000+(2*n2+256)*(MIN(64,n1)+4.56)
        na2=60000+(2*n3+256)*(MIN(64,n1*n2)+4.56)
     ELSE
        na1=60000+4.56*n1+(2*n2+256)*(MIN(64,n1)+4.56)
        naux=60000+4.56*n1+(2*n3+256)*(MIN(64,n1*n2)+4.56) 
     END IF
     IF( n2>=252 .AND. n3<252 ) THEN
        naux=na1
     ELSE IF( n2<252 .AND. n3>=252 ) THEN
        naux=na2
     ELSE
        naux=MAX(na1,na2)
     END IF
  END IF

!    ALLOCATE(aux(naux),STAT=isos)
!    IF ( isos /= 0 ) CALL stop_prg( 'fft_wrap','failed to allocate aux')

!    sign_fft = fsign
!    CALL dcft3(zg,la1,la2,zgout,la1,la2,n1,n2,n3,sign_fft,norm,aux,naux)

  IF ( PRESENT ( zgout ) ) THEN
     CALL do_ibm_fft_ab ( naux )
  ELSE
     CALL do_ibm_fft_aa ( naux )
  END IF

!    DEALLOCATE(aux,STAT=isos)
!    IF ( isos /= 0 ) CALL stop_prg( 'fft_wrap','failed to deallocate aux')

CONTAINS

  SUBROUTINE do_ibm_fft_ab ( naux )

    INTEGER, INTENT ( IN ) :: naux

    REAL(kind=dbl), DIMENSION ( naux ) :: aux

    sign_fft = fsign
    CALL scft3(zg,la1,la2,zgout,la1,la2,n1,n2,n3,sign_fft,norm,aux,naux)

  END SUBROUTINE do_ibm_fft_ab

  SUBROUTINE do_ibm_fft_aa ( naux )

    INTEGER, INTENT ( IN ) :: naux

    REAL(kind=dbl), DIMENSION ( naux ) :: aux

    sign_fft = fsign
    CALL scft3(zg,la1,la2,zg,la1,la2,n1,n2,n3,sign_fft,norm,aux,naux)

  END SUBROUTINE do_ibm_fft_aa

END SUBROUTINE fft_wrap_sgl

!******************************************************************************

#elif defined ( __T3E )

SUBROUTINE fft_wrap_ab ( fsign, n, zg, zgout, scale )
! routine with wrapper for all fft calls:
! Does transform with exp(+ig.r*sign):
  IMPLICIT NONE
  INTEGER, INTENT ( IN ):: fsign
  INTEGER, INTENT ( IN ), DIMENSION ( : ) :: n
  COMPLEX ( dbl ), DIMENSION(0:,0:,0:), INTENT ( IN ):: zg
  COMPLEX ( dbl ), DIMENSION(0:,0:,0:), INTENT ( INOUT ):: zgout
  REAL ( dbl ), OPTIONAL, INTENT ( IN ):: scale
! locals
  INTEGER, SAVE :: n1save = -1, n2save = -1, n3save = -1
  INTEGER:: la1,la2,la3,sign_fft,na1,na2,isos, isys ( 4 )
  INTEGER :: n1,n2,n3
  REAL ( dbl ):: norm
  REAL(kind=dbl), DIMENSION ( : ), ALLOCATABLE :: work
  REAL(kind=dbl), DIMENSION ( : ), ALLOCATABLE, SAVE :: table

  IF( PRESENT ( scale) ) THEN
     norm=scale
  ELSE
     norm=1._dbl
  END IF

  n1 = n(1)
  n2 = n(2)
  n3 = n(3)
  la1 = n1
  la2 = n2
  isys ( 1 ) = 3
  isys ( 2 ) = 0
  isys ( 3 ) = 0
  isys ( 4 ) = 0

  ALLOCATE ( work ( 2 * n1 * n2 * n3 ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_prg( 'fft_wrap','failed to allocate work')
  sign_fft = fsign

  IF ( n1save /= n1 .OR. n2save /= n2 .OR. n3save /= n3 ) THEN
     IF ( ALLOCATED ( table ) ) DEALLOCATE ( table )
     ALLOCATE ( table ( 2 * ( n1 + n2 + n3 ) ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_prg( 'fft_wrap','failed to allocate table')
     CALL ccfft3d ( 0, n1, n2, n3, norm, zg, la1, la2, &
          zgout, la1, la2, table, work, isys )
     n1save = n1
     n2save = n2
     n3save = n3
  END IF

  CALL ccfft3d ( sign_fft, n1, n2, n3, norm, zg, la1, la2, &
       zgout, la1, la2, table, work, isys )

  DEALLOCATE ( work, STAT = isos )
  IF ( isos /= 0 ) CALL stop_prg( 'fft_wrap','failed to deallocate work')

END SUBROUTINE fft_wrap_ab

!******************************************************************************

SUBROUTINE fft_wrap_aa ( fsign, n, zg, scale )
  
! routine with wrapper for all fft calls:
! Does transform with exp(+ig.r*sign):
  
  IMPLICIT NONE
  INTEGER, INTENT ( IN ):: fsign
  INTEGER, INTENT ( IN ), DIMENSION ( : ) :: n
  COMPLEX ( dbl ), DIMENSION(0:,0:,0:), INTENT ( INOUT ):: zg
  REAL ( dbl ), OPTIONAL, INTENT ( IN ):: scale
! locals
  INTEGER, SAVE :: n1save = -1, n2save = -1, n3save = -1
  INTEGER:: la1,la2,la3,sign_fft,na1,na2,isos, isys ( 4 )
  INTEGER :: n1,n2,n3
  REAL ( dbl ):: norm
  REAL(kind=dbl), DIMENSION ( : ), ALLOCATABLE :: work
  REAL(kind=dbl), DIMENSION ( : ), ALLOCATABLE, SAVE :: table

  IF( PRESENT ( scale) ) THEN
     norm=scale
  ELSE
     norm=1._dbl
  END IF

  n1 = n(1)
  n2 = n(2)
  n3 = n(3)
  la1 = n1
  la2 = n2
  isys ( 1 ) = 3
  isys ( 2 ) = 0
  isys ( 3 ) = 0
  isys ( 4 ) = 0

  ALLOCATE ( work ( 2 * n1 * n2 * n3 ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_prg( 'fft_wrap','failed to allocate work')
  sign_fft = fsign

  IF ( n1save /= n1 .OR. n2save /= n2 .OR. n3save /= n3 ) THEN
     IF ( ALLOCATED ( table ) ) DEALLOCATE ( table )
     ALLOCATE ( table ( 2 * ( n1 + n2 + n3 ) ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_prg( 'fft_wrap','failed to allocate table')
     CALL ccfft3d ( 0, n1, n2, n3, norm, zg, la1, la2, &
          zg, la1, la2, table, work, isys )
     n1save = n1
     n2save = n2
     n3save = n3
  END IF

  CALL ccfft3d ( sign_fft, n1, n2, n3, norm, zg, la1, la2, &
       zg, la1, la2, table, work, isys )

  DEALLOCATE ( work, STAT = isos )
  IF ( isos /= 0 ) CALL stop_prg( 'fft_wrap','failed to deallocate work')

END SUBROUTINE fft_wrap_aa

#endif

END MODULE fft_tools
