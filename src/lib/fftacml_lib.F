!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fftacml_lib [1.0] *
!!
!!   NAME
!!     fftacml_lib
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH 3-Feb-2006 : initial version
!!     Adapted for new interface 
!!
!!*****
!******************************************************************************

SUBROUTINE fftacml_get_lengths ( data, max_length )

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

END SUBROUTINE fftacml_get_lengths

!******************************************************************************

SUBROUTINE fftacml3d ( fft_in_place, fsign, scale, n, zin, zout )

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

  INTEGER                                              :: n1, n2, n3, ll, info
  COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE          :: comm

!------------------------------------------------------------------------------

  n1 = n(1)
  n2 = n(2)
  n3 = n(3)

#if defined ( __FFTACML )

  ll = n1*n2*n3 + 5*(n1+n2+n3) + 300
  ALLOCATE ( comm ( ll ) )

#if defined(__SGL) || defined(__FFTSGL)
  CALL cfft3dx(0,scale,.TRUE.,fft_in_place,n1,n2,n3,zin,zout,comm,info)
  CALL cfft3dx(-fsign,scale,.TRUE.,fft_in_place,n1,n2,n3,zin,zout,comm,info)
#else
  CALL zfft3dx(0,scale,.TRUE.,fft_in_place,n1,n2,n3,zin,zout,comm,info)
  CALL zfft3dx(-fsign,scale,.TRUE.,fft_in_place,n1,n2,n3,zin,zout,comm,info)
#endif

  IF ( info /= 0 ) fsign = 0

  DEALLOCATE ( comm )
  
#else

  fsign = 0

#endif

END SUBROUTINE fftacml3d

!******************************************************************************

SUBROUTINE fftacml1dm ( fsign, trans, n, m, zin, zout, scale )

  IMPLICIT NONE
#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER :: dp = KIND(0.0)
#else
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
#endif

! Arguments
  INTEGER, INTENT ( INOUT )                          :: fsign
  LOGICAL, INTENT ( IN )                             :: trans
  INTEGER, INTENT ( IN )                             :: n, m
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( INOUT )  :: zin
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( OUT )    :: zout
  REAL (KIND=dp), INTENT ( IN )                      :: scale

  INTEGER                                     :: info, ll
  INTEGER                                     :: incx1, incx2, incy1, incy2
  COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE :: comm

!------------------------------------------------------------------------------

  incx1 = 1
  incx2 = n
  incy1 = 1
  incy2 = n
  IF ( fsign == +1 .AND. trans ) THEN
    incx1 = m
    incx2 = 1
  ELSEIF ( fsign == -1 .AND. trans ) THEN
    incy1 = m
    incy2 = 1
  END IF

#if defined ( __FFTACML )

  ll = 3*n + 100
  ALLOCATE ( comm ( ll ) )

#if defined(__SGL) || defined(__FFTSGL)
  CALL cfft1mx (0,scale,.FALSE.,m,n,zin,incx1,incx2,zout,incy1,incy2,comm,info)
  CALL cfft1mx (-fsign,scale,.FALSE.,m,n,zin,incx1,incx2,zout,incy1,incy2,&
                comm,info)
#else
  CALL zfft1mx (0,scale,.FALSE.,m,n,zin,incx1,incx2,zout,incy1,incy2,comm,info)
  CALL zfft1mx (-fsign,scale,.FALSE.,m,n,zin,incx1,incx2,zout,incy1,incy2,&
                comm,info)
#endif

  IF ( info /= 0 ) fsign=0
  
  DEALLOCATE ( comm )

#else

  fsign = 0

#endif

END SUBROUTINE fftacml1dm

!******************************************************************************
