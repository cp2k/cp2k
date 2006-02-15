!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fftw3_lib [1.0] *
!!
!!   NAME
!!     fftw3_lib
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH 23-Jan-2006 : initial version
!!     Adapted for new interface
!!
!!*****
!******************************************************************************

SUBROUTINE fftw3_get_lengths ( data, max_length )

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

END SUBROUTINE fftw3_get_lengths

!******************************************************************************

SUBROUTINE fftw33d ( fft_in_place, fsign, scale, n, zin, zout )

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

#if defined ( __FFTW3 )
#include "fftw3.f"
#endif

#if defined (__NAG)
  INTEGER*8          :: dummy_64_bit
  INTEGER, PARAMETER :: integer8_kind = KIND(dummy_64_bit)
#else
  INTEGER, PARAMETER :: integer8_kind = 8
#endif

  INTEGER                                  :: n1, n2, n3
  INTEGER(KIND=integer8_kind)              :: plan

!------------------------------------------------------------------------------

  n1 = n(1)
  n2 = n(2)
  n3 = n(3)

#if defined ( __FFTW3 )

#if defined(__SGL) || defined(__FFTSGL)
  IF ( fft_in_place .AND. fsign == +1 ) THEN
    CALL sfftw_plan_dft_3d(plan,n1,n2,n3,zin,zin,FFTW_FORWARD,FFTW_ESTIMATE)
  ELSE IF ( fft_in_place .AND. fsign == -1 ) THEN
    CALL sfftw_plan_dft_3d(plan,n1,n2,n3,zin,zin,FFTW_BACKWARD,FFTW_ESTIMATE)
  ELSE IF ( .NOT. fft_in_place .AND. fsign == +1 ) THEN
    CALL sfftw_plan_dft_3d(plan,n1,n2,n3,zin,zout,FFTW_FORWARD,FFTW_ESTIMATE)
  ELSE IF ( .NOT. fft_in_place .AND. fsign == -1 ) THEN
    CALL sfftw_plan_dft_3d(plan,n1,n2,n3,zin,zout,FFTW_BACKWARD,FFTW_ESTIMATE)
  END IF
  CALL sfftw_execute(plan)
  CALL sfftw_destroy_plan(plan)
#else
  IF ( fft_in_place .AND. fsign == +1 ) THEN
    CALL dfftw_plan_dft_3d(plan,n1,n2,n3,zin,zin,FFTW_FORWARD,FFTW_ESTIMATE)
  ELSE IF ( fft_in_place .AND. fsign == -1 ) THEN
    CALL dfftw_plan_dft_3d(plan,n1,n2,n3,zin,zin,FFTW_BACKWARD,FFTW_ESTIMATE)
  ELSE IF ( .NOT. fft_in_place .AND. fsign == +1 ) THEN
    CALL dfftw_plan_dft_3d(plan,n1,n2,n3,zin,zout,FFTW_FORWARD,FFTW_ESTIMATE)
  ELSE IF ( .NOT. fft_in_place .AND. fsign == -1 ) THEN
    CALL dfftw_plan_dft_3d(plan,n1,n2,n3,zin,zout,FFTW_BACKWARD,FFTW_ESTIMATE)
  END IF
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
#endif

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

END SUBROUTINE fftw33d

!******************************************************************************

SUBROUTINE fftw31dm ( fsign, trans, n, m, zin, zout, scale )

  IMPLICIT NONE
#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER :: dp = KIND(0.0)
#else
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
#endif

#if defined ( __FFTW3 )
#include "fftw3.f"
#endif

  INTEGER, INTENT ( INOUT )                          :: fsign
  LOGICAL, INTENT ( IN )                             :: trans
  INTEGER, INTENT ( IN )                             :: n, m
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( INOUT )  :: zin
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( OUT )    :: zout
  REAL (KIND=dp), INTENT ( IN )                      :: scale

#if defined (__NAG)
  INTEGER*8          :: dummy_64_bit
  INTEGER, PARAMETER :: integer8_kind = KIND(dummy_64_bit)
#else
  INTEGER, PARAMETER :: integer8_kind = 8
#endif
  INTEGER                                            :: ii,di,io,do
  INTEGER(KIND=integer8_kind)                        :: plan

!------------------------------------------------------------------------------
  
#if defined ( __FFTW3 )

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

#if defined(__SGL) || defined(__FFTSGL)
  IF ( fsign == +1 ) THEN
    CALL sfftw_plan_many_dft(plan,1,n,m,zin,0,ii,di,&
              zout,0,io,do,FFTW_FORWARD,FFTW_ESTIMATE)
  ELSE
    CALL sfftw_plan_many_dft(plan,1,n,m,zin,0,ii,di,&
              zout,0,io,do,FFTW_BACKWARD,FFTW_ESTIMATE)
  END IF
  CALL sfftw_execute(plan)
  CALL sfftw_destroy_plan(plan)
  IF ( scale /= 1.0_dp ) CALL csscal(n*m,scale,zout,1)
#else
  IF ( fsign == +1 ) THEN
    CALL dfftw_plan_many_dft(plan,1,n,m,zin,0,ii,di,&
              zout,0,io,do,FFTW_FORWARD,FFTW_ESTIMATE)
  ELSE
    CALL dfftw_plan_many_dft(plan,1,n,m,zin,0,ii,di,&
              zout,0,io,do,FFTW_BACKWARD,FFTW_ESTIMATE)
  END IF
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  IF ( scale /= 1.0_dp ) CALL zdscal(n*m,scale,zout,1)
#endif

#else

  fsign = 0

#endif

END SUBROUTINE fftw31dm

!******************************************************************************
