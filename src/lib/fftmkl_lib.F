!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fftmkl_lib [1.0] *
!!
!!   NAME
!!     fftmkl_lib
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH 10-Feb-2006 : initial version
!!
!!*****
!******************************************************************************

SUBROUTINE fftmkl_get_lengths ( data, max_length )

  IMPLICIT NONE 

  INTEGER, DIMENSION(*)                    :: data
  INTEGER, INTENT(INOUT)                   :: max_length

END SUBROUTINE fftmkl_get_lengths

!******************************************************************************

SUBROUTINE fftmkl3d ( fft_in_place, fsign, scale, n, zin, zout )

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

!------------------------------------------------------------------------------

#if defined ( __FFTMKL )

  fsign = 0
  
#else

  fsign = 0

#endif

END SUBROUTINE fftmkl3d

!******************************************************************************

SUBROUTINE fftmkl1dm ( fsign, trans, n, m, zin, zout, scale )

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

!------------------------------------------------------------------------------

#if defined ( __FFTMKL )

  fsign = 0

#else

  fsign = 0

#endif

END SUBROUTINE fftmkl1dm

!******************************************************************************
