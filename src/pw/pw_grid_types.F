!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      JGH (20-12-2000) : Parallel data layout
!> \author APSI
! **************************************************************************************************
MODULE pw_grid_types

   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_cart_type,&
                                              mp_comm_null,&
                                              mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: pw_grid_type, map_pn

   ! (only for reciprocal grid:) fill in half or full space
   INTEGER, PARAMETER, PUBLIC :: HALFSPACE = 211, FULLSPACE = 212
   INTEGER, PARAMETER, PUBLIC :: PW_MODE_LOCAL = 0, PW_MODE_DISTRIBUTED = 1

   ! maps to positive and negative g-vectors
! **************************************************************************************************
   TYPE map_pn
      INTEGER, DIMENSION(:), ALLOCATABLE :: pos, neg
   END TYPE map_pn

! info on parallelisation
   ! contains only significant information if mode == PW_MODE_DISTRIBUTED
! **************************************************************************************************
   TYPE pw_para_type
      INTEGER :: mode = PW_MODE_LOCAL ! 0 = local = PW_MODE_LOCAL ; 1 = distributed = PW_MODE_DISTRIBUTED
      LOGICAL :: ray_distribution = .FALSE. ! block or pencil distribution
      LOGICAL :: blocked = .FALSE. ! block or pencil distribution
      TYPE(mp_comm_type) :: group = mp_comm_null ! MPI group
      INTEGER :: my_pos = -1! Position within group
      INTEGER :: group_size = -1! # of Processors in group
      LOGICAL :: group_head = .FALSE. ! Master process within group
      INTEGER :: group_head_id = 0 ! Id of group_head
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: yzp ! g-space rays (xy,k,pe)
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: yzq ! local inverse pointer of yzp
      INTEGER, DIMENSION(:), ALLOCATABLE :: nyzray ! number of g-space rays (pe)
      TYPE(mp_cart_type) :: rs_group = mp_cart_type() ! real space group (2-dim cart)
      INTEGER :: rs_mpo = -1 ! real space group position
      INTEGER, DIMENSION(2) :: rs_dims = 0 ! real space group dimensions
      INTEGER, DIMENSION(2) :: rs_pos = -1 ! real space group positions in grid
      INTEGER, DIMENSION(:, :, :, :), ALLOCATABLE :: bo ! list of axis distribution
      INTEGER, DIMENSION(:), ALLOCATABLE :: pos_of_x ! what my_pos holds a given x plane....should go: hard-codes to plane distributed
   END TYPE pw_para_type

   ! all you always wanted to know about grids, but were...
! **************************************************************************************************
   TYPE pw_grid_type
      INTEGER(int_8) :: ngpts = 0_int_8 ! # grid points
      INTEGER(int_8) :: ngpts_cut = 0_int_8 ! # grid points within cutoff
      INTEGER, DIMENSION(2, 3) :: bounds = 0 ! lower and upper bounds
      INTEGER, DIMENSION(3) :: npts = 0 ! # point in all directions
      INTEGER :: ngpts_local = 0 ! # grid points
      INTEGER :: ngpts_cut_local = 0 ! # grid points within cutoff
      INTEGER, DIMENSION(2, 3) :: bounds_local = 0 ! bounds on local process
      INTEGER, DIMENSION(3) :: npts_local = 0 ! local version of npts
      REAL(KIND=dp), DIMENSION(3) :: dr = 0.0_dp ! grid spacing
      REAL(KIND=dp), DIMENSION(3, 3) :: dh = 0.0_dp ! incremental cell matrix
      REAL(KIND=dp), DIMENSION(3, 3) :: dh_inv = 0.0_dp ! inverse incremental cell matrix
      LOGICAL :: orthorhombic = .TRUE. ! cell symmetry
      REAL(KIND=dp) :: dvol = 0.0_dp, vol = 0.0_dp ! volume element, volume
      REAL(KIND=dp) :: cutoff = 0.0_dp ! cutoff in a.u.
      TYPE(map_pn) :: mapl = map_pn(), mapm = map_pn(), mapn = map_pn() ! mapping 1D => 3D
      TYPE(pw_para_type) :: para = pw_para_type() ! information on parallelisation
      REAL(KIND=dp), DIMENSION(:, :), POINTER, CONTIGUOUS :: g => NULL() ! grid point vectors
      REAL(KIND=dp), DIMENSION(:), POINTER, CONTIGUOUS :: gsq => NULL() ! squared vector lengths
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: g_hat ! grid point indices (Miller)
      INTEGER, DIMENSION(:, :), POINTER, CONTIGUOUS :: g_hatmap => NULL() ! mapped grid point indices (Miller) [CUDA]
      INTEGER :: grid_span = FULLSPACE ! type HALFSPACE/FULLSPACE
      LOGICAL :: have_g0 = .TRUE. ! whether I have G = [0,0,0]
      INTEGER :: first_gne0 = 0 ! first g index /= 0 [1/2]
      INTEGER :: id_nr = -1 ! tag of this grid
      INTEGER :: reference = 0 ! reference grid identifier
      INTEGER, DIMENSION(:), POINTER :: gidx => NULL() ! ref grid index
      INTEGER :: ref_count = 0 ! reference count
      LOGICAL :: spherical = .FALSE. ! spherical cutoff?
      COMPLEX(KIND=dp), DIMENSION(:, :), CONTIGUOUS, POINTER :: grays => NULL() ! used by parallel 3D FFT routine
   END TYPE pw_grid_type

END MODULE pw_grid_types

