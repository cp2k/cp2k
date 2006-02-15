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
!!     Adapted for new interface
!!
!******************************************************************************

SUBROUTINE fftessl_get_lengths ( data, max_length )

  IMPLICIT NONE 

  INTEGER, DIMENSION(*)                    :: data
  INTEGER, INTENT(INOUT)                   :: max_length

  INTEGER, PARAMETER                       :: rlen =135
  INTEGER, DIMENSION(rlen), PARAMETER :: radix = (/ &
    2,   4,   6,   8,  10,  12,  14,  16,  18,  20,  22,  24,  28,  30,  32,&
   36,  40,  42,  44,  48,  56,  60,  64,  66,  70,  72,  80,  84,  88,  90,&
   96, 110, 112, 120, 126, 128, 132, 140, 144, 154, 160, 168, 176, 180, 192,&
  198, 210, 220, 224, 240, 252, 256, 264, 280, 288, 308, 320, 330, 336, 352,&
  360, 384, 396, 420, 440, 448, 462, 480, 504, 512, 528, 560, 576, 616, 630,& 
  640, 660, 672, 704, 720, 768, 770, 792, 840, 880, 896, 924, 960, 990,1008,&
 1024,1056,1120,1152,1232,1260,1280,1320,1344,1386,1408,1440,1536,1540,1584,&
 1680,1760,1792,1848,1920,1980,2016,2048,2112,2240,2304,2310,2464,2520,2560,&
 2640,2688,2772,2816,2880,3072,3080,3168,3360,3520,3584,3696,3840,3960,4032 /)

 INTEGER                                   :: ndata

!------------------------------------------------------------------------------

  ndata = MIN ( rlen, max_length )
  data(1:ndata) = radix(1:ndata)
  max_length = ndata

END SUBROUTINE fftessl_get_lengths

!******************************************************************************

SUBROUTINE fftessl3d ( fft_in_place, fsign, scale, n, zin, zout )

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

  COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)          :: aux
  INTEGER                                              :: naux, nx, ny, nz
  INTEGER                                              :: na1, na2

!------------------------------------------------------------------------------

  nx = n ( 1 )
  ny = n ( 2 )
  nz = n ( 3 )

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

  ALLOCATE ( aux ( naux ) )

  IF ( fft_in_place ) THEN

#if defined ( __SGL ) || defined ( __FFTSGL )
    CALL scft3(zin,nx,nx*ny,zin,nx,nx*ny,nx,ny,nz,fsign,scale,aux,naux)
#else
    CALL dcft3(zin,nx,nx*ny,zin,nx,nx*ny,nx,ny,nz,fsign,scale,aux,naux)
#endif

  ELSE

#if defined ( __SGL ) || defined ( __FFTSGL )
    CALL scft3(zin,nx,nx*ny,zout,nx,nx*ny,nx,ny,nz,fsign,scale,aux,naux)
#else
    CALL dcft3(zin,nx,nx*ny,zout,nx,nx*ny,nx,ny,nz,fsign,scale,aux,naux)
#endif

  END IF

  DEALLOCATE ( aux )

#else

  fsign = 0

#endif

END SUBROUTINE fftessl3d

!******************************************************************************

SUBROUTINE fftessl1dm ( fsign, trans, n, m, zin, zout, scale )

  IMPLICIT NONE
#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER :: dp = KIND(0.0)
#else
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
#endif

  INTEGER, INTENT ( INOUT )                          :: fsign
  LOGICAL, INTENT ( IN )                             :: trans
  INTEGER, INTENT ( IN )                             :: n, m
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( INOUT )  :: zin
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( OUT )    :: zout
  REAL (KIND=dp), INTENT ( IN )                      :: scale

  COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: aux1, aux2
  INTEGER                                            :: naux1, naux2, init
  INTEGER                                            :: inc1x, inc2x, inc1y, inc2y

!------------------------------------------------------------------------------

  IF( n < 2048 ) THEN
     naux1 = 20000
  ELSE
     naux1 = 20000 + NINT ( 2.28_dp * REAL ( n,KIND=dp) )
  END IF
  IF ( trans ) THEN
     naux1 = naux1 + ( 2 * n + 256 ) * MIN ( 64, n )
  END IF
  naux2 = naux1

  inc1x = 1
  inc2x = n
  inc1y = 1
  inc2y = n
  IF ( fsign == +1 .AND. trans ) THEN
    inc1x = m
    inc2x = 1
  ELSEIF ( fsign == -1 .AND. trans ) THEN
    inc1y = m
    inc2y = 1
  END IF

#if defined ( __FFTESSL )

  ALLOCATE ( aux1 ( naux1 ), aux2 ( naux2 ) )

#if defined ( __SGL ) || defined ( __FFTSGL )
  init = 1
  CALL SCFT ( init, zin, inc1x, inc2x, zout, inc1y, inc2y, &
              n, m, fsign, scale, aux1, naux1, aux2, naux2 )
  init = 0
  CALL SCFT ( init, zin, inc1x, inc2x, zout, inc1y, inc2y, &
              n, m, fsign, scale, aux1, naux1, aux2, naux2 )
#else
  init = 1
  CALL DCFT ( init, zin, inc1x, inc2x, zout, inc1y, inc2y, &
              n, m, fsign, scale, aux1, naux1, aux2, naux2 )
  init = 0
  CALL DCFT ( init, zin, inc1x, inc2x, zout, inc1y, inc2y, &
              n, m, fsign, scale, aux1, naux1, aux2, naux2 )
#endif

  DEALLOCATE ( aux1, aux2 )

#else

  fsign = 0

#endif

END SUBROUTINE fftessl1dm

!******************************************************************************
