!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****** cp2k/fftw2_lib [1.0] *
!!
!!   NAME
!!     fftw2_lib
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH 20-Feb-2001 : mltfft, mfft2d routines added
!!                       integer_pointer_kind
!!     Adapted to new interface style
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE fftw2_get_lengths ( data, max_length )

  IMPLICIT NONE 

  INTEGER, DIMENSION(*)                    :: data
  INTEGER, INTENT(INOUT)                   :: max_length

  INTEGER, DIMENSION(:), ALLOCATABLE       :: dlocal, idx
  INTEGER                                  :: h, i, j, k, m, maxn, &
                                              maxn_elevens, maxn_fives, &
                                              maxn_sevens, maxn_thirteens, &
                                              maxn_threes, maxn_twos, &
                                              ndata, nmax, number

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

  ALLOCATE ( dlocal ( ndata ), idx ( ndata ) )

  ndata = 0
  dlocal ( : ) = 0
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
                 dlocal ( ndata ) = number
              END DO
           END DO
        END DO
     END DO
  END DO

  CALL sortint ( dlocal, ndata, idx )
  ndata = MIN ( ndata, max_length )
  data(1:ndata) = dlocal(1:ndata)
  max_length = ndata

  DEALLOCATE ( dlocal, idx )

END SUBROUTINE fftw2_get_lengths

!******************************************************************************

SUBROUTINE fftw23d ( fft_in_place, fsign, scale, n, zin, zout )

  IMPLICIT NONE
#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER :: dp = KIND(0.0)
#else
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
#endif
  LOGICAL, INTENT(IN)                                  :: fft_in_place
  INTEGER, INTENT(INOUT)                               :: fsign
  REAL(KIND=dp), INTENT(IN)                            :: scale
  INTEGER, DIMENSION(*), INTENT(IN)                    :: n
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)        :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)        :: zout

! this is a verbatim include of the fortran include file, as this is
! not installed with fftw2 :
! fftw-2.1.5/fortran> cat fftw_f77.i
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
! end include fftw_f77.i



! integer_pointer_kind: kind parameter for pointer variables (addresses)
!                       this is usually an integer with 32 or 64 bit size
!                       needs improving, this will cause problems
#if defined (__NAG)
  INTEGER*8          :: dummy_64_bit
  INTEGER, PARAMETER :: integer_pointer_kind = KIND(dummy_64_bit)
#else
  INTEGER, PARAMETER :: integer_pointer_kind = 8
#endif

  INTEGER                                      :: n1, n2, n3, sign_fft
  INTEGER(KIND=integer_pointer_kind)           :: plan

!------------------------------------------------------------------------------

  n1 = n(1)
  n2 = n(2)
  n3 = n(3)

  sign_fft = fsign

#if defined ( __FFTW2 )

  IF ( fft_in_place ) THEN
    IF ( sign_fft == +1 ) THEN
      CALL fftw3d_f77_create_plan ( plan, n1, n2, n3, FFTW_FORWARD, &
                                    FFTW_ESTIMATE + FFTW_IN_PLACE )
    ELSE
      CALL fftw3d_f77_create_plan ( plan, n1, n2, n3, FFTW_BACKWARD, &
                                    FFTW_ESTIMATE + FFTW_IN_PLACE )
    END IF
  ELSE
    IF ( sign_fft == +1 ) THEN
      CALL fftw3d_f77_create_plan ( plan, n1, n2, n3, FFTW_FORWARD, &
                                    FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
    ELSE
      CALL fftw3d_f77_create_plan ( plan, n1, n2, n3, FFTW_BACKWARD, &
                                    FFTW_ESTIMATE + FFTW_OUT_OF_PLACE )
    END IF
  END IF

  CALL fftwnd_f77_one ( plan, zin, zout )
  CALL fftw_f77_destroy_plan(plan)

  IF ( scale /= 1.0_dp ) THEN
#if defined(__SGL) || defined(__FFTSGL)
    IF ( fft_in_place ) THEN
      CALL csscal(n1*n2*n3,scale,zin,1)
    ELSE
      CALL csscal(n1*n2*n3,scale,zout,1)
    END IF
#else
    IF ( fft_in_place ) THEN
      CALL zdscal(n1*n2*n3,scale,zin,1)
    ELSE
      CALL zdscal(n1*n2*n3,scale,zout,1)
    END IF
#endif
  END IF

#else

  fsign = 0

#endif

END SUBROUTINE fftw23d

!******************************************************************************

SUBROUTINE fftw21dm ( fsign, trans, n, m, zin, zout, scale )

  IMPLICIT NONE
#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER :: dp = KIND(0.0)
#else
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
#endif

! this is a verbatim include of the fortran include file, as this is
! not installed with fftw2 :
! fftw-2.1.5/fortran> cat fftw_f77.i
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
! end include fftw_f77.i



! integer_pointer_kind: kind parameter for pointer variables (addresses)
!                       this is usually an integer with 32 or 64 bit size
!                       needs improving, this will cause problems
#if defined (__NAG)
  INTEGER*8          :: dummy_64_bit
  INTEGER, PARAMETER :: integer_pointer_kind = KIND(dummy_64_bit)
#else
  INTEGER, PARAMETER :: integer_pointer_kind = 8
#endif

  INTEGER, INTENT ( INOUT )                          :: fsign
  LOGICAL, INTENT ( IN )                             :: trans
  INTEGER, INTENT ( IN )                             :: n, m
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( INOUT )  :: zin
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( OUT )    :: zout
  REAL (KIND=dp), INTENT ( IN )                      :: scale

  INTEGER                                  :: di, do, ii, io
  INTEGER(KIND=integer_pointer_kind)       :: plan

!------------------------------------------------------------------------------

#if defined ( __FFTW2 )

  ii = 1
  di = n
  io = 1
  do = n
  IF ( fsign == +1 .AND. trans ) THEN
    ii = m
    di = 1
  ELSEIF ( fsign == -1 .AND. trans ) THEN
    io = m
    do = 1
  END IF
  IF ( fsign == +1 ) THEN
    CALL fftw_f77_create_plan(plan,n,FFTW_FORWARD,FFTW_ESTIMATE)
  ELSE
    CALL fftw_f77_create_plan(plan,n,FFTW_BACKWARD,FFTW_ESTIMATE)
  END IF
  CALL fftw_f77(plan,m,zin,ii,di,zout,io,do)
  CALL fftw_f77_destroy_plan(plan)
#if defined(__SGL) || defined(__FFTSGL)
  IF ( scale /= 1.0_dp ) CALL csscal(n*m,scale,zout,1)
#else
  IF ( scale /= 1.0_dp ) CALL zdscal(n*m,scale,zout,1)
#endif

#else

  fsign = 0

#endif

END SUBROUTINE fftw21dm
