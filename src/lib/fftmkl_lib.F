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
#if defined ( __FFTMKL )
#include "mkl_dfti.f90"
#endif

! *****************************************************************************
SUBROUTINE fftmkl_get_lengths ( DATA, max_length )

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

END SUBROUTINE fftmkl_get_lengths

! *****************************************************************************
SUBROUTINE fftmkl3d ( fft_in_place, fsign, scale, n, zin, zout )

#if defined ( __FFTMKL )
  USE MKL_DFTI
#endif
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

#if defined ( __FFTMKL )
  TYPE(DFTI_DESCRIPTOR), POINTER                       :: desc
#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER :: dft_prec = DFTi_SINGLE
#else
  INTEGER, PARAMETER :: dft_prec = DFTi_DOUBLE
#endif
#endif
  INTEGER                                              :: status
!------------------------------------------------------------------------------

#if defined ( __FFTMKL )

  status = DftiCreateDescriptor( desc, dft_prec, DFTI_COMPLEX, 3, n )
  IF (fsign == +1) THEN
    IF ( scale /= 1.0_dp ) status = DftiSetValue( desc, DFTI_FORWARD_SCALE, scale )
    IF ( fft_in_place ) THEN
      status = DftiSetValue( desc, DFTI_PLACEMENT, DFTI_INPLACE )
      status = DftiCommitDescriptor( desc )
      status = DftiComputeForward( desc, zin )
    ELSE
      status = DftiSetValue( desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
      status = DftiCommitDescriptor( desc )
      status = DftiComputeForward( desc, zin, zout )
    ENDIF
  ELSE
    IF ( scale /= 1.0_dp ) status = DftiSetValue( desc, DFTI_BACKWARD_SCALE, scale )
    IF ( fft_in_place ) THEN
      status = DftiSetValue( desc, DFTI_PLACEMENT, DFTI_INPLACE )
      status = DftiCommitDescriptor( desc )
      status = DftiComputeBackward( desc, zin )
    ELSE
      status = DftiSetValue( desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
      status = DftiCommitDescriptor( desc )
      status = DftiComputeBackward( desc, zin, zout )
    ENDIF
  END IF
  status = DftiFreeDescriptor( desc )
  
#else

  fsign = 0

#endif

END SUBROUTINE fftmkl3d

! *****************************************************************************
SUBROUTINE fftmkl1dm ( fsign, trans, n, m, zin, zout, scale )

#if defined ( __FFTMKL )
  USE MKL_DFTI
#endif
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
#if defined ( __FFTMKL )
  TYPE(DFTI_DESCRIPTOR), POINTER                       :: desc
#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER :: dft_prec = DFTi_SINGLE
#else
  INTEGER, PARAMETER :: dft_prec = DFTi_DOUBLE
#endif
#endif
  INTEGER                                              :: status
  INTEGER                                            :: ii(3),di,io(3),DO

!------------------------------------------------------------------------------

#if defined ( __FFTMKL )

  ii(1) = 0
  ii(2) = 1
  ii(3) = 0
  di = n
  io(1) = 0
  io(2) = 1
  io(3) = 0
  DO = n
  IF ( fsign == +1 .AND. trans ) THEN
    ii(2) = m
    di = 1
  ELSEIF ( fsign == -1 .AND. trans ) THEN
    io(2) = m
    DO = 1
  END IF

  status = DftiCreateDescriptor( desc, dft_prec, DFTI_COMPLEX, 1, n )
  status = DftiSetValue( desc, DFTI_INPUT_STRIDES, ii )
  status = DftiSetValue( desc, DFTI_OUTPUT_STRIDES, io )
  status = DftiSetValue( desc, DFTI_NUMBER_OF_TRANSFORMS, m )
  status = DftiSetValue( desc, DFTI_INPUT_DISTANCE, di )
  status = DftiSetValue( desc, DFTI_OUTPUT_DISTANCE, DO )
  status = DftiSetValue( desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
  IF (fsign == +1) THEN
    IF ( scale /= 1.0_dp ) status = DftiSetValue( desc, DFTI_FORWARD_SCALE, scale )
    status = DftiCommitDescriptor( desc )
    status = DftiComputeForward( desc, zin, zout )
  ELSE
    IF ( scale /= 1.0_dp ) status = DftiSetValue( desc, DFTI_BACKWARD_SCALE, scale )
    status = DftiCommitDescriptor( desc )
    status = DftiComputeBackward( desc, zin, zout )
  END IF
  status = DftiFreeDescriptor( desc )

#else

  fsign = 0

#endif

END SUBROUTINE fftmkl1dm

