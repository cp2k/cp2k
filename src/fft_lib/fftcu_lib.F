!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH 10-Feb-2006 : initial version
!> \par
!>  This file is part of the mkl distribution and can be found in the
!>  corresponding include directory
!>  add this to the path to the -I flags of the CPP
!> \author JGH
! *****************************************************************************

SUBROUTINE fftcu_get_lengths ( DATA, max_length )

  IMPLICIT NONE

  INTEGER, DIMENSION(*)                    :: DATA
  INTEGER, INTENT(INOUT)                   :: max_length

  INTEGER, DIMENSION(:), ALLOCATABLE       :: dlocal, idx
  INTEGER                                  :: i, maxn, &
                                              maxn_fac, &
                                              ndata, number

  maxn_fac = 15
  maxn = 37748736

  ndata = 0
  DO i=0,maxn_fac
    number = 2**i
    IF ( number >= maxn ) CYCLE
    ndata = ndata + 1
  ENDDO

  DO i=0,maxn_fac
    number = 3**i
    IF ( number >= maxn ) CYCLE
    ndata = ndata + 1
  ENDDO

  DO i=0,maxn_fac
    number = 5**i
    IF ( number >= maxn ) CYCLE
    ndata = ndata + 1
  ENDDO

  DO i=0,maxn_fac
    number = 7**i
    IF ( number >= maxn ) CYCLE
    ndata = ndata + 1
  ENDDO

  ALLOCATE ( dlocal ( ndata ), idx ( ndata ) )

  ndata = 0
  dlocal ( : ) = 0
  DO i=0,maxn_fac
    number = 2**i
    IF ( number >= maxn ) CYCLE
    ndata = ndata + 1
    dlocal ( ndata ) = number
  ENDDO

  DO i=0,maxn_fac
    number = 3**i
    IF ( number >= maxn ) CYCLE
    ndata = ndata + 1
    dlocal ( ndata ) = number
  ENDDO

  DO i=0,maxn_fac
    number = 5**i
    IF ( number >= maxn ) CYCLE
    ndata = ndata + 1
    dlocal ( ndata ) = number
  ENDDO

  DO i=0,maxn_fac
    number = 7**i
    IF ( number >= maxn ) CYCLE
    ndata = ndata + 1
    dlocal ( ndata ) = number
  ENDDO

  CALL sortint ( dlocal, ndata, idx )
  ndata = MIN ( ndata, max_length )
  DATA(1:ndata) = dlocal(1:ndata)
  max_length = ndata

  DEALLOCATE ( dlocal, idx )

END SUBROUTINE fftcu_get_lengths

! *****************************************************************************
SUBROUTINE fftcu3d ( fft_in_place, fsign, scale, n, zin, zout )

  USE fft_kinds, ONLY                                   : dp
  LOGICAL, INTENT(IN)                                  :: fft_in_place
  INTEGER, INTENT(INOUT)                               :: fsign
  REAL(KIND=dp), INTENT(IN)                            :: scale
  INTEGER, DIMENSION(*), INTENT(IN)                    :: n
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)        :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)        :: zout

#if defined ( __FFTCU )  || defined ( __CUDA )

  WRITE(*,*) 'FFTCU only works in single precision'
  CALL flush(6)
  STOP

#else
  fsign = 0
#endif

END SUBROUTINE fftcu3d

! *****************************************************************************
SUBROUTINE fftcu1dm ( fsign, trans, n, m, zin, zout, scale )

  USE fft_kinds, ONLY                                   : dp
  INTEGER, INTENT ( INOUT )                          :: fsign
  LOGICAL, INTENT ( IN )                             :: trans
  INTEGER, INTENT ( IN )                             :: n, m
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( INOUT )  :: zin
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( OUT )    :: zout
  REAL (KIND=dp), INTENT ( IN )                      :: scale


#if defined ( __FFTCU ) || defined ( __CUDA )

  WRITE(*,*) 'FFTCU only works in single precision'
  CALL flush(6)
  STOP

#else

  fsign = 0

#endif

END SUBROUTINE fftcu1dm

