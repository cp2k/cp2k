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
!!     JGH 20-Feb-2001 : mltfft, mfft2d routines added
!!                       integer_pointer_kind
!!
!!   SOURCE
!******************************************************************************

MODULE fftw_lib

  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory
  USE util,                            ONLY: sort

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fft3d, mltfft, mfft2d
  PUBLIC :: fft_get_lengths
  
! integer_pointer_kind: kind parameter for pointer variables (addresses)
!                       this is usually an integer with 32 or 64 bit size
! needs improving, this will cause problems
#if defined (__ABSOFT)
  INTEGER, PARAMETER :: integer_pointer_kind = 4
#elif defined (__NAG)
  INTEGER*4          :: dummy_32_bit                  
  INTEGER*8          :: dummy_64_bit                  
  INTEGER, PARAMETER :: integer_pointer_kind = KIND(dummy_64_bit)
#else
  INTEGER, PARAMETER :: integer_pointer_kind = 8
#endif

!*apsi 220500 This is actually the file fftw_f77.i ...
!     This file contains PARAMETER statements for various constants
!     that can be passed to FFTW routines.  You should include
!     this file in any FORTRAN program that calls the fftw_f77
!     routines (either directly or with an #include statement
!     if you use the C preprocessor).

      INTEGER FFTW_FORWARD,FFTW_BACKWARD
      PARAMETER (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

      INTEGER FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
      PARAMETER (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

      INTEGER FFTW_ESTIMATE,FFTW_MEASURE
      PARAMETER (FFTW_ESTIMATE=0,FFTW_MEASURE=1)

      INTEGER FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
      PARAMETER (FFTW_OUT_OF_PLACE=0)
      PARAMETER (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

      INTEGER FFTW_THREADSAFE
      PARAMETER (FFTW_THREADSAFE=128)

!     Constants for the MPI wrappers:
      INTEGER FFTW_TRANSPOSED_ORDER, FFTW_NORMAL_ORDER
      INTEGER FFTW_SCRAMBLED_INPUT, FFTW_SCRAMBLED_OUTPUT
      PARAMETER(FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0)
      PARAMETER(FFTW_SCRAMBLED_INPUT=8192)
      PARAMETER(FFTW_SCRAMBLED_OUTPUT=16384)

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

SUBROUTINE fft_get_lengths ( DATA, max_length )


    INTEGER, DIMENSION(:), POINTER           :: DATA
    INTEGER, INTENT(IN)                      :: max_length

    INTEGER :: allocstat, h, i, j, k, m, maxn, maxn_elevens, maxn_fives, &
      maxn_sevens, maxn_thirteens, maxn_threes, maxn_twos, ndata, nmax, number
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: idx

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

  ALLOCATE ( DATA ( ndata ), idx ( ndata ), STAT = allocstat )
  IF ( allocstat /= 0 ) THEN
     CALL stop_memory ( "fft_get_lengths", "data, idx", 2*ndata )
  END IF

  ndata = 0
  DATA ( : ) = 0
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


    INTEGER, INTENT(INOUT)                   :: fsign
    REAL(KIND=dp), INTENT(IN)                :: scale
    INTEGER, DIMENSION(:), INTENT(IN)        :: n
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: zg
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT), OPTIONAL                :: zg_out

    COMPLEX(KIND=dp), DIMENSION(1, 1, 1)     :: zgout_tmp
    INTEGER                                  :: n1, n2, n3, sign_fft
    INTEGER(KIND=integer_pointer_kind), SAVE :: plan_a_bw = 0, plan_a_fw = 0, &
      plan_b_bw = 0, plan_b_fw = 0, plan_c_bw = 0, plan_c_fw = 0
    INTEGER, SAVE :: n1a_save = -1, n1b_save = -1, n1c_save = -1, &
      n2a_save = -1, n2b_save = -1, n2c_save = -1, n3a_save = -1, &
      n3b_save = -1, n3c_save = -1
    LOGICAL                                  :: fft_in_place
    LOGICAL, SAVE                            :: ffta_in_place = .TRUE., &
                                                fftb_in_place = .TRUE., &
                                                fftc_in_place = .TRUE.
    REAL(KIND=dp)                            :: norm

!
! Just due to DEC... apsi
!------------------------------------------------------------------------------

  norm = scale

  n1 = n(1)
  n2 = n(2)
  n3 = n(3)
  IF ( PRESENT ( zg_out ) ) THEN
     fft_in_place = .FALSE.
  ELSE
     fft_in_place = .TRUE.
  END IF

  sign_fft = fsign

#if defined ( __FFTW )

  IF ( n1a_save == n1 .AND. n2a_save == n2 .AND. n3a_save == n3 &
       .AND. ( fft_in_place .EQV. ffta_in_place ) ) THEN
     IF ( sign_fft == +1 ) THEN
        IF (fft_in_place) THEN
           CALL fftwnd_f77_one ( plan_a_fw, zg, zgout_tmp )
        ELSE
           CALL fftwnd_f77_one ( plan_a_fw, zg, zg_out )
        ENDIF
     ELSE
        IF (fft_in_place) THEN
            CALL fftwnd_f77_one ( plan_a_bw, zg, zgout_tmp )
        ELSE
            CALL fftwnd_f77_one ( plan_a_bw, zg, zg_out )
        ENDIF
     END IF
  ELSE IF ( n1b_save == n1 .AND. n2b_save == n2 .AND. n3b_save == n3 &
       .AND. ( fft_in_place .EQV. fftb_in_place ) ) THEN
     IF ( sign_fft == +1 ) THEN
        IF (fft_in_place) THEN
            CALL fftwnd_f77_one ( plan_b_fw, zg, zgout_tmp )
        ELSE
            CALL fftwnd_f77_one ( plan_b_fw, zg, zg_out )
        ENDIF
     ELSE
        IF (fft_in_place) THEN
            CALL fftwnd_f77_one ( plan_b_bw, zg, zgout_tmp )
        ELSE
            CALL fftwnd_f77_one ( plan_b_bw, zg, zg_out )
        ENDIF
     END IF
  ELSE IF ( n1c_save == n1 .AND. n2c_save == n2 .AND. n3c_save == n3 &
       .AND. ( fft_in_place .EQV. fftc_in_place ) ) THEN
     IF ( sign_fft == +1 ) THEN
        IF (fft_in_place) THEN
            CALL fftwnd_f77_one ( plan_c_fw, zg, zgout_tmp )
        ELSE
            CALL fftwnd_f77_one ( plan_c_fw, zg, zg_out )
        ENDIF
     ELSE
        IF (fft_in_place) THEN
            CALL fftwnd_f77_one ( plan_c_bw, zg, zgout_tmp )
        ELSE
            CALL fftwnd_f77_one ( plan_c_bw, zg, zg_out )
        ENDIF
     END IF
  ELSE IF ( n1a_save == -1 .OR. &
       ( n1a_save == n1 .AND. n2a_save == n2 .AND. n3a_save == n3 &
       .AND. ( fft_in_place .NEQV. ffta_in_place ) .AND. &
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
         IF (fft_in_place) THEN
             CALL fftwnd_f77_one ( plan_a_fw, zg, zgout_tmp )
         ELSE
             CALL fftwnd_f77_one ( plan_a_fw, zg, zg_out )
         ENDIF
     ELSE
         IF (fft_in_place) THEN
             CALL fftwnd_f77_one ( plan_a_bw, zg, zgout_tmp )
         ELSE
             CALL fftwnd_f77_one ( plan_a_bw, zg, zg_out )
         ENDIF
     END IF

  ELSE IF ( n1b_save == -1 .OR. &
       ( n1b_save == n1 .AND. n2b_save == n2 .AND. n3b_save == n3 &
       .AND. ( fft_in_place .NEQV. fftb_in_place ) .AND. &
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
        IF (fft_in_place) THEN
            CALL fftwnd_f77_one ( plan_b_fw, zg, zgout_tmp )
        ELSE
            CALL fftwnd_f77_one ( plan_b_fw, zg, zg_out )
        ENDIF
     ELSE
        IF (fft_in_place) THEN
            CALL fftwnd_f77_one ( plan_b_bw, zg, zgout_tmp )
        ELSE
            CALL fftwnd_f77_one ( plan_b_bw, zg, zg_out )
        ENDIF
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
        IF (fft_in_place) THEN
            CALL fftwnd_f77_one ( plan_c_fw, zg, zgout_tmp )
        ELSE
            CALL fftwnd_f77_one ( plan_c_fw, zg, zg_out )
        ENDIF
     ELSE
        IF (fft_in_place) THEN
            CALL fftwnd_f77_one ( plan_c_bw, zg, zgout_tmp )
        ELSE
            CALL fftwnd_f77_one ( plan_c_bw, zg, zg_out )
        ENDIF
     END IF
  END IF

#else

  fsign = 0

#endif

  IF ( norm /= 1.0_dp ) THEN
     IF ( fft_in_place ) THEN
        zg = zg * norm
     ELSE
        zg_out = zg_out * norm
     END IF
  END IF

END SUBROUTINE fft3d

!!*****
!******************************************************************************
!!****** fftw_lib/mfft2d [1.0] *
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

  IMPLICIT NONE

! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: rin, rout
  INTEGER, INTENT ( INOUT ) :: fsign
  REAL (KIND=dp), INTENT ( IN ) :: scale
  INTEGER, INTENT ( IN ) :: n1, n2, nm
  COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT ( IN ) :: zin
  COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT ( OUT ) :: zout

! Locals
  INTEGER :: li1, li2, lo1, lo2, ierr, li, lo, i
  COMPLEX (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: zz
  INTEGER (KIND=integer_pointer_kind) :: plan

!------------------------------------------------------------------------------

#if defined ( __FFTW )

  li1 = SIZE ( zin, 1 )
  li2 = SIZE ( zin, 2 )
  lo1 = SIZE ( zout, 1 )
  lo2 = SIZE ( zout, 2 )

  IF ( fsign == +1 ) THEN
    CALL fftw2d_f77_create_plan ( plan, n1, n2, FFTW_FORWARD, &
       FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
  ELSE
    CALL fftw2d_f77_create_plan ( plan, n1, n2, FFTW_BACKWARD, &
       FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
  END IF

  IF ( rin(1:1) == "T" .OR.  rin(1:1) == "t" ) THEN
    IF ( rout(1:1) == "T" .OR.  rout(1:1) == "t" ) THEN
!..T-T
      li = li1 * li2
      lo = lo1 * lo2
      CALL fftwnd_f77 ( plan, nm, zin, li1, 1, zout, lo1, 1 )
      IF ( scale /= 1.0_dp ) THEN
        zout ( 1:nm ,1:n1, 1:n2 ) = scale * zout ( 1:nm, 1:n1, 1:n2 )
      END IF

    ELSE
!..T-N
      li = li1 * li2
      lo = lo1 * lo2
      CALL fftwnd_f77 ( plan, nm, zin, li1, 1, zout, 1, lo )
      IF ( scale /= 1.0_dp ) THEN
        zout ( 1:n1, 1:n2, 1:nm ) = scale * zout ( 1:n1, 1:n2, 1:nm )
      END IF

    END IF
  ELSE
    IF ( rout(1:1) == "T" .OR.  rout(1:1) == "t" ) THEN
!..N-T
      li = li1 * li2
      lo = lo1 * lo2
!     ALLOCATE ( zz ( li1, li2, nm ), STAT = ierr )
!     IF ( ierr /= 0 ) call stop_memory ( "fft2dm", "zz", li*nm )
!     CALL fftwnd_f77 ( plan, nm, zin, 1, li, zz, 1, li )
!     IF ( scale /= 1.0_dp ) THEN
!       DO i = 1, nm
!         zout ( i, 1:n1, 1:n2 ) = scale * zz ( 1:n1, 1:n2, i )
!       END DO
!     ELSE
!       DO i = 1, nm
!         zout ( i, 1:n1, 1:n2 ) = zz ( 1:n1, 1:n2, i )
!       END DO
!     END IF
!     DEALLOCATE ( zz, STAT = ierr )
!     IF ( ierr /= 0 ) call stop_memory ( "fft2dm", "zz" )
      li = li1 * li2
      lo = lo1 * lo2
      CALL fftwnd_f77 ( plan, nm, zin, 1, li, zout, lo1, 1 )
      IF ( scale /= 1.0_dp ) THEN
        zout ( 1:nm, 1:n1, 1:n2 ) = scale * zout ( 1:nm, 1:n1, 1:n2 )
      END IF

    ELSE
!..N-N
      li = li1 * li2
      lo = lo1 * lo2
      CALL fftwnd_f77 ( plan, nm, zin, 1, li, zout, 1, lo )
      IF ( scale /= 1.0_dp ) THEN
        zout ( 1:n1, 1:n2, 1:nm ) = scale * zout ( 1:n1, 1:n2, 1:nm )
      END IF

    END IF
  END IF

  CALL fftwnd_f77_destroy_plan(plan)

#else

  fsign = 0

#endif

END SUBROUTINE mfft2d

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

SUBROUTINE mltfft ( tin, tout, fsign, scale, n, m, zin, zout )
  

    CHARACTER(LEN=*), INTENT(IN)             :: tin, tout
    INTEGER, INTENT(INOUT)                   :: fsign
    REAL(KIND=dp), INTENT(IN)                :: scale
    INTEGER, INTENT(IN)                      :: n, m
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: zin
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: zout

    INTEGER                                  :: di, DO, ii, io, ldi, ldo
    INTEGER(KIND=integer_pointer_kind)       :: plan

!------------------------------------------------------------------------------

  ldi = SIZE ( zin, 1 )
  ldo = SIZE ( zout, 1 )
  
#if defined ( __FFTW )

  IF ( fsign == +1 ) THEN
    CALL fftw_f77_create_plan(plan,n,FFTW_FORWARD,FFTW_ESTIMATE)
  ELSE
    CALL fftw_f77_create_plan(plan,n,FFTW_BACKWARD,FFTW_ESTIMATE)
  END IF
  IF ( tin(1:1) == 'N' .OR. tin(1:1) == 'N' ) THEN
    ii = 1
    di = ldi
  ELSE
    ii = ldi
    di = 1
  ENDIF
  IF ( tout(1:1) == 'N' .OR. tout(1:1) == 'N' ) THEN
    io = 1
    DO = ldo
  ELSE
    io = ldo
    DO = 1
  ENDIF
  CALL fftw_f77(plan,m,zin,ii,di,zout,io,DO)
  CALL fftw_f77_destroy_plan(plan)
  IF ( scale /= 1.0_dp ) zout = scale * zout

#else

  fsign = 0

#endif

END SUBROUTINE mltfft

!!*****
!******************************************************************************

END MODULE fftw_lib

!******************************************************************************
