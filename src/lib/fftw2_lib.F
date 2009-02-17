!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH 20-Feb-2001 : mltfft, mfft2d routines added
!>                        integer_pointer_kind
!>      Adapted to new interface style
!>      IAB 09-Jan-2009 : Modified to cache plans in fft_plan_type
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!> \author JGH
! *****************************************************************************
SUBROUTINE fftw2_get_lengths ( DATA, max_length )

  IMPLICIT NONE 

  INTEGER, DIMENSION(*)                    :: DATA
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
  maxn = 37748736

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
  DATA(1:ndata) = dlocal(1:ndata)
  max_length = ndata

  DEALLOCATE ( dlocal, idx )

END SUBROUTINE fftw2_get_lengths

! *****************************************************************************
SUBROUTINE fftw2_create_plan_3d ( plan, plan_style )

  USE fft_plan,                           ONLY: fft_plan_type, dp

  IMPLICIT NONE

! this is a verbatim include of the fortran include file, as this is
! not installed with fftw2 :
! fftw-2.1.5/fortran> cat fftw_f77.i
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
! end include fftw_f77.i

  TYPE(fft_plan_type), INTENT ( INOUT )              :: plan
  INTEGER, INTENT (IN)                               :: plan_style
  INTEGER                                            :: n1, n2, n3, sign_fft

  INTEGER :: fftw_plan_type
  SELECT CASE(plan_style)
  CASE(1)
         fftw_plan_type = FFTW_ESTIMATE
  CASE(2)
         fftw_plan_type = FFTW_MEASURE
  CASE(3) ! patient does not exist, use measure instead
         fftw_plan_type = FFTW_MEASURE
  CASE(4) ! exhaustive does not exist, use measure instead
         fftw_plan_type = FFTW_MEASURE
  CASE DEFAULT
   STOP "fftw2_create_plan_3d"
  END SELECT

  n1 = plan%n_3d(1)
  n2 = plan%n_3d(2)
  n3 = plan%n_3d(3)

  sign_fft = plan%fsign

#if defined ( __FFTW2 )

  IF ( plan%fft_in_place ) THEN
    IF ( sign_fft == +1 ) THEN
      CALL fftw3d_f77_create_plan ( plan%fftw_plan, n1, n2, n3, FFTW_FORWARD, &
                                    fftw_plan_type + FFTW_IN_PLACE )
    ELSE
      CALL fftw3d_f77_create_plan ( plan%fftw_plan, n1, n2, n3, FFTW_BACKWARD, &
                                    fftw_plan_type + FFTW_IN_PLACE )
    END IF
  ELSE
    IF ( sign_fft == +1 ) THEN
      CALL fftw3d_f77_create_plan ( plan%fftw_plan, n1, n2, n3, FFTW_FORWARD, &
                                    fftw_plan_type )
    ELSE
      CALL fftw3d_f77_create_plan ( plan%fftw_plan, n1, n2, n3, FFTW_BACKWARD, &
                                    fftw_plan_type )
    END IF
  END IF

#endif

END SUBROUTINE fftw2_create_plan_3d

SUBROUTINE fftw23d ( plan, scale, zin, zout )

  USE fft_plan,                           ONLY: fft_plan_type, dp

  IMPLICIT NONE

  TYPE(fft_plan_type), INTENT ( INOUT )              :: plan
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)         :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)      :: zout
  REAL (KIND=dp), INTENT ( IN )                      :: scale
  INTEGER                                            :: n1,n2,n3

#if defined(__FFTW2)

  n1 = plan%n_3d(1)
  n2 = plan%n_3d(2)
  n3 = plan%n_3d(3)

  CALL fftwnd_f77_one ( plan%fftw_plan, zin, zout )

  IF ( scale /= 1.0_dp ) THEN
#if defined(__SGL) || defined(__FFTSGL)
    IF ( plan%fft_in_place ) THEN
      CALL csscal(n1*n2*n3,scale,zin,1)
    ELSE
      CALL csscal(n1*n2*n3,scale,zout,1)
    END IF
#else
    IF ( plan%fft_in_place ) THEN
      CALL zdscal(n1*n2*n3,scale,zin,1)
    ELSE
      CALL zdscal(n1*n2*n3,scale,zout,1)
    END IF
#endif
  END IF

#else

  plan%fsign = 0

#endif

END SUBROUTINE fftw23d

! *****************************************************************************
SUBROUTINE fftw2_create_plan_1dm ( plan, plan_style )

  USE fft_plan,                           ONLY: fft_plan_type, dp

  IMPLICIT NONE

! this is a verbatim include of the fortran include file, as this is
! not installed with fftw2 :
! fftw-2.1.5/fortran> cat fftw_f77.i
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
! end include fftw_f77.i

  TYPE(fft_plan_type), INTENT ( INOUT )              :: plan
  INTEGER :: plan_style
  
  INTEGER :: fftw_plan_type
  SELECT CASE(plan_style)
  CASE(1)
         fftw_plan_type = FFTW_ESTIMATE
  CASE(2)
         fftw_plan_type = FFTW_MEASURE
  CASE(3) ! patient does not exist, use measure instead
         fftw_plan_type = FFTW_MEASURE
  CASE(4) ! exhaustive does not exist, use measure instead
         fftw_plan_type = FFTW_MEASURE
  CASE DEFAULT
   STOP "fftw2_create_plan_1dm"
  END SELECT

#if defined ( __FFTW2 )
  IF ( plan%fsign == +1 ) THEN
    CALL fftw_f77_create_plan(plan%fftw_plan,plan%n,FFTW_FORWARD,fftw_plan_type)
  ELSE
    CALL fftw_f77_create_plan(plan%fftw_plan,plan%n,FFTW_BACKWARD,fftw_plan_type)
  END IF
#endif

END SUBROUTINE fftw2_create_plan_1dm

SUBROUTINE fftw2_destroy_plan ( plan )

  USE fft_plan,                           ONLY: fft_plan_type

  IMPLICIT NONE

  TYPE(fft_plan_type), INTENT (INOUT)   :: plan

#if defined ( __FFTW2 )
  CALL fftw_f77_destroy_plan(plan%fftw_plan)
#endif

END SUBROUTINE fftw2_destroy_plan

SUBROUTINE fftw21dm ( plan, zin, zout, scale )

  USE fft_plan,                           ONLY: fft_plan_type, dp

  IMPLICIT NONE

  TYPE(fft_plan_type), INTENT ( INOUT )              :: plan  
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN)         :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)      :: zout
  REAL (KIND=dp), INTENT ( IN )                      :: scale

  INTEGER                                  :: di, DO, ii, io

!------------------------------------------------------------------------------

#if defined ( __FFTW2 )

  ii = 1
  di = plan%n
  io = 1
  DO = plan%n
  IF ( plan%fsign == +1 .AND. plan%trans ) THEN
    ii = plan%m
    di = 1
  ELSEIF ( plan%fsign == -1 .AND. plan%trans ) THEN
    io = plan%m
    DO = 1
  END IF

  CALL fftw_f77(plan%fftw_plan,plan%m,zin,ii,di,zout,io,DO)

#if defined(__SGL) || defined(__FFTSGL)
  IF ( scale /= 1.0_dp ) CALL csscal(plan%n*plan%m,scale,zout,1)
#else
  IF ( scale /= 1.0_dp ) CALL zdscal(plan%n*plan%m,scale,zout,1)
#endif

#else

  plan%fsign = 0

#endif

END SUBROUTINE fftw21dm
