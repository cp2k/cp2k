!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH 10-Feb-2006 : initial version
!> \author JGH
! *****************************************************************************
SUBROUTINE fftsci_get_lengths ( DATA, max_length )

  IMPLICIT NONE 

  INTEGER, DIMENSION(*)                    :: DATA
  INTEGER, INTENT(INOUT)                   :: max_length

!------------------------------------------------------------------------------

END SUBROUTINE fftsci_get_lengths

! *****************************************************************************
SUBROUTINE fftsci3d ( fft_in_place, fsign, scale, n, zin, zout )

  USE fft_kinds, ONLY                                   : dp
  LOGICAL, INTENT(IN)                                  :: fft_in_place
  INTEGER, INTENT(INOUT)                               :: fsign
  REAL(KIND=dp), INTENT(IN)                            :: scale
  INTEGER, DIMENSION(*), INTENT(IN)                    :: n
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)        :: zin
  COMPLEX(KIND=dp), DIMENSION(*), INTENT(INOUT)        :: zout

!------------------------------------------------------------------------------

#if defined ( __FFTSCI )

  fsign = 0
  
#else

  fsign = 0

#endif

END SUBROUTINE fftsci3d

! *****************************************************************************
SUBROUTINE fftsci1dm ( fsign, trans, n, m, zin, zout, scale )

  USE fft_kinds, ONLY                                   : dp
  INTEGER, INTENT ( INOUT )                          :: fsign
  LOGICAL, INTENT ( IN )                             :: trans
  INTEGER, INTENT ( IN )                             :: n, m
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( INOUT )  :: zin
  COMPLEX (KIND=dp), DIMENSION(*), INTENT ( OUT )    :: zout
  REAL (KIND=dp), INTENT ( IN )                      :: scale

!------------------------------------------------------------------------------

#if defined ( __FFTSCI )

  fsign = 0

#else

  fsign = 0

#endif

END SUBROUTINE fftsci1dm

