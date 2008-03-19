!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006 - 2008 CP2K developers group                           !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH 10-Feb-2006 : initial version
!> \par
!>  This file is part of the mkl distribution and can be found in the 
!>  corresponding include directory
!>  add this to the path to the -I flags of the CPP
!> \author JGH
!#if defined ( __FFTCU ) || defined ( __CUDA )
!#include "mkl_dfti.f90"

!#endif

! *****************************************************************************
SUBROUTINE fftcu_get_lengths ( DATA, max_length )

  IMPLICIT NONE 

  INTEGER, DIMENSION(*)                    :: DATA
  INTEGER, INTENT(INOUT)                   :: max_length

  INTEGER, DIMENSION(:), ALLOCATABLE       :: dlocal, idx
  INTEGER                                  :: h, i, j, k, m, maxn, &
                                              maxn_fac, &
                                              ndata, nmax, number

!------------------------------------------------------------------------------

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

  INTEGER                                              :: ifft_in_place

  REAL(KIND=dp), POINTER                               :: ptr
  INTEGER                                              :: status, lmem
!------------------------------------------------------------------------------

#if defined ( __FFTCU )  || defined ( __CUDA )
#if defined(__SGL) || defined(__FFTSGL)

  ifft_in_place=0
  if (fft_in_place) ifft_in_place=1

  CALL fftcu3d_cu(ifft_in_place, fsign, scale, n, zin, zout)

#else

  write(*,*) 'FFTCU only works in single precision'
  call flush(6)
  stop

#endif
#else
  fsign = 0
#endif

END SUBROUTINE fftcu3d

! *****************************************************************************
SUBROUTINE fftcu1dm ( fsign, trans, n, m, zin, zout, scale )

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
  INTEGER                                              :: status, itrans
  INTEGER                                            :: ii(3),di,io(3),DO

!------------------------------------------------------------------------------

#if defined ( __FFTCU ) || defined ( __CUDA )
#if defined(__SGL) || defined(__FFTSGL)

  itrans=0
!  if (trans) itrans=1
  if (trans) then
    write(*,*) 'FFTCU 1D tranpose FFT has not been implimented'
    call flush(6)
    stop
  endif

  write(*,*) 'WARNING: FFTCU 1D FFT has not been tested!'
  CALL fftcu1dm_cu(fsign, itrans, n, m, zin, zout, scale)

#else

  write(*,*) 'FFTCU only works in single precision'
  call flush(6)
  stop

#endif
#else

  fsign = 0

#endif

END SUBROUTINE fftcu1dm

