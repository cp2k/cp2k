!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fftw_lib [1.0] *
!!
!!   NAME
!!     fftw_lib
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

MODULE fftw_lib
  
  USE kinds, ONLY: dbl
  USE stop_program, ONLY : stop_memory, stop_prg
  USE util, ONLY:  sort
  
  IMPLICIT NONE 
  
  PRIVATE
  PUBLIC :: fft3d
  PUBLIC :: fft_get_lengths
  
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

!!*****
!******************************************************************************

CONTAINS
  
!******************************************************************************
!!****** fftw_lib/fft_get_lengths [1.0] *
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
!! FFTW can do arbitrary(?) lenghts, maybe you want to limit them to some
!!    powers of small prime numbers though...
  
  maxn_twos = 15
  maxn_threes = 3
  maxn_fives = 2
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
!!****** fftw_lib/fft3d [1.0] *
!!
!!   NAME
!!     fft3d
!!
!!   FUNCTION
!!     Routine with wrapper for 3dfft calls:
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
  REAL ( dbl ), INTENT ( IN ) :: scale
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: n
  COMPLEX ( dbl ), DIMENSION(:,:,:), INTENT ( INOUT ) :: zg
  COMPLEX ( dbl ), DIMENSION(:,:,:), INTENT ( INOUT ), OPTIONAL, TARGET :: zg_out
  
! Locals
  INTEGER, SAVE :: n1a_save = -1, n2a_save = -1, n3a_save = -1
  INTEGER, SAVE :: n1b_save = -1, n2b_save = -1, n3b_save = -1
  INTEGER, SAVE :: n1c_save = -1, n2c_save = -1, n3c_save = -1
  LOGICAL, SAVE :: ffta_in_place = .TRUE., fftb_in_place = .TRUE.
  LOGICAL, SAVE :: fftc_in_place = .TRUE.
  LOGICAL :: fft_in_place
  INTEGER :: sign_fft, n1, n2, n3
!
! KIND=8 does not work on all machines !!
!
  INTEGER ( KIND = 8 ), SAVE :: plan_a_fw, plan_a_bw
  INTEGER ( KIND = 8 ), SAVE :: plan_b_fw, plan_b_bw
  INTEGER ( KIND = 8 ), SAVE :: plan_c_fw, plan_c_bw
  REAL ( dbl ) :: norm
  
  ! Just due to DEC... apsi
  COMPLEX ( dbl ), DIMENSION(:,:,:), POINTER :: zgout
  COMPLEX ( dbl ), DIMENSION(1,1,1), TARGET :: zgout_tmp
  
!------------------------------------------------------------------------------
  
  norm = scale
  
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

#if defined ( __FFTW )  

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

#else

  fsign = 0

#endif

  IF ( norm /= 1._dbl ) THEN
     IF ( fft_in_place ) THEN
        zg = zg * norm
     ELSE
        zgout = zgout * norm
     END IF
  END IF

END SUBROUTINE fft3d

!!*****
!******************************************************************************
!!****** fftw_lib/mltfft [1.0] *
!!
!!   NAME
!!     mltfft
!!
!!   FUNCTION
!!     Calls multiple 1d FFT from the FFTw library
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
  
#if defined ( __FFTW )

  stop " fftw_lib "

#else

  fsign = 0

#endif

END SUBROUTINE mltfft

!!*****
!******************************************************************************

END MODULE fftw_lib

!******************************************************************************
