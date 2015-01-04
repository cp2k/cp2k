!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH (20-12-2000) : Parallel data layout
!> \author APSI
! *****************************************************************************
MODULE pw_grid_types
  
  USE kinds,                           ONLY: dp,&
                                             int_8
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pw_grid_type, map_pn

  ! (only for reciprocal grid:) fill in half or full space
  INTEGER, PARAMETER, PUBLIC :: HALFSPACE = 211, FULLSPACE = 212
  INTEGER, PARAMETER, PUBLIC :: PW_MODE_LOCAL = 0, PW_MODE_DISTRIBUTED=1


  ! maps to positive and negative g-vectors
! *****************************************************************************
  TYPE map_pn
     INTEGER, DIMENSION ( : ), POINTER :: pos, neg
  END TYPE map_pn

! info on parallelisation
  ! contains only significant information if mode == PW_MODE_DISTRIBUTED
! *****************************************************************************
  TYPE pw_para_type
     INTEGER :: mode                 ! 0 = local = PW_MODE_LOCAL ; 1 = distributed = PW_MODE_DISTRIBUTED
     LOGICAL :: ray_distribution     ! block or pencil distribution
     LOGICAL :: blocked              ! block or pencil distribution
     INTEGER :: group                ! MPI group id
     INTEGER :: my_pos               ! Position within group
     INTEGER :: group_size           ! # of Processors in group
     LOGICAL :: group_head           ! Master process within group
     INTEGER :: group_head_id        ! Id of group_head
     INTEGER, DIMENSION (:,:,:), POINTER :: yzp  ! g-space rays (xy,k,pe)
     INTEGER, DIMENSION (:,:), POINTER :: yzq    ! local inverse pointer of yzp
     INTEGER, DIMENSION (:), POINTER :: nyzray   ! number of g-space rays (pe)
     INTEGER :: rs_group             ! real space group (2-dim cart)
     INTEGER :: rs_mpo               ! real space group position
     INTEGER, DIMENSION (2) :: rs_dims    ! real space group dimensions
     INTEGER, DIMENSION (2) :: rs_pos     ! real space group positions in grid
     INTEGER, DIMENSION (:,:,:,:), POINTER :: bo ! list of axis distribution
     INTEGER, DIMENSION(:), POINTER :: pos_of_x ! what my_pos holds a given x plane....should go since this hard-codes to plane distributed
  END TYPE pw_para_type

  ! all you always wanted to know about grids, but were...
! *****************************************************************************
  TYPE pw_grid_type
     INTEGER(int_8) :: ngpts                             ! # grid points
     INTEGER(int_8) :: ngpts_cut                         ! # grid points within cutoff
     INTEGER, DIMENSION ( 2, 3 ) :: bounds        ! lower and upper bounds
     INTEGER, DIMENSION ( 3 ) :: npts             ! # point in all directions
     INTEGER :: ngpts_local                       ! # grid points
     INTEGER :: ngpts_cut_local                   ! # grid points within cutoff
     INTEGER, DIMENSION ( 2, 3 ) :: bounds_local  ! bounds on local process
     INTEGER, DIMENSION ( 3 ) :: npts_local       ! local version of npts
     REAL (KIND=dp), DIMENSION ( 3 ) :: dr        ! grid spacing
     REAL (KIND=dp), DIMENSION ( 3, 3 ) :: dh     ! incremental cell matrix
     REAL (KIND=dp), DIMENSION ( 3, 3 ) :: dh_inv ! inverse incremental cell matrix
     LOGICAL :: orthorhombic                      ! cell symmetry
     REAL (KIND=dp) :: dvol, vol                  ! volume element, volume
     REAL (KIND=dp) :: cutoff                     ! cutoff in a.u.
     TYPE ( map_pn ) :: mapl, mapm, mapn          ! mapping 1D => 3D
     TYPE ( pw_para_type ) :: para                ! information on parallelisation
     REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: g   ! grid point vectors
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: gsq    ! squared vector lengths
     INTEGER, DIMENSION ( :, : ), POINTER :: g_hat    ! grid point indices (Miller)
     INTEGER, DIMENSION ( :, : ), POINTER :: g_hatmap ! mapped grid point indices (Miller) [CUDA]
     INTEGER :: grid_span                         ! type HALFSPACE/FULLSPACE
     LOGICAL :: have_g0                           ! whether I have G = [0,0,0]
     INTEGER :: first_gne0                        ! first g index /= 0 [1/2]
     INTEGER :: id_nr                             ! tag of this grid
     INTEGER :: reference                         ! reference grid identifier
     INTEGER, DIMENSION ( : ), POINTER :: gidx    ! ref grid index
     INTEGER :: ref_count                         ! reference count
     LOGICAL :: spherical                         ! spherical cutoff?
     COMPLEX (KIND=dp), DIMENSION ( :, : ), POINTER :: grays ! used by parallel 3D FFT routine
  END TYPE pw_grid_type

END MODULE pw_grid_types

