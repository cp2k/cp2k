!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/pw_grid_types [1.0] *
!!
!!   NAME
!!     pw_grid_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE pw_grid_types
  
  USE kinds, ONLY : dbl
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: pw_grid_type, map_pn
  PUBLIC :: HALFSPACE, FULLSPACE
  
  !! maps to positive and negative g-vectors
  TYPE map_pn
     INTEGER, DIMENSION ( : ), POINTER :: pos, neg
  END TYPE map_pn
  
  !! all you always wanted to know about grids, but were...
  TYPE pw_grid_type
     INTEGER :: ngpts                             ! # grid points
     INTEGER :: ngpts_cut                         ! # grid points within cutoff
     INTEGER :: ndimension                        ! # of dimensions
     INTEGER, DIMENSION ( 2, 3 ) :: bounds        ! lower and upper bounds
     INTEGER, DIMENSION ( 3 ) :: npts             ! # point in all directions
     REAL ( dbl ), DIMENSION ( 3 ) :: dr          ! grid spacing
     REAL ( dbl ) :: dvol, vol                    ! volume element, volume
     TYPE ( map_pn ) :: mapl, mapm, mapn          ! mapping 1D => 3D
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: g   ! grid point vectors
     REAL ( dbl ), DIMENSION ( : ), POINTER :: gsq   ! squared vector lengths
     INTEGER, DIMENSION ( :, : ), POINTER :: g_hat    ! grid point indices (Miller)
!    INTEGER :: grid_kind                         ! type of grid: equidistance,
     !   logarithmic, ...
     INTEGER :: grid_span                         ! type of equidistance grid:
                         !   is f^*(g) = f(-g)
     
     LOGICAL :: have_g0                           ! whether I have G = [0,0,0]
     INTEGER :: first_gne0                        ! first g index /= 0 [1/2]
     
     INTEGER :: nglengths                         ! # different lengths/shells
     REAL ( dbl ), DIMENSION ( : ), POINTER :: glsq ! squared vector lengths
     INTEGER, DIMENSION ( : ), POINTER :: glidx   ! index to the full arrays
     
  END TYPE pw_grid_type
  
!! (only for reciprocal grid:) fill in half or full space
  INTEGER, PARAMETER :: HALFSPACE = 211, FULLSPACE = 212
  
END MODULE pw_grid_types

!******************************************************************************
!!*****
