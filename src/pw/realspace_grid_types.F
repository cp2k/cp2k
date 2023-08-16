!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \note
!>      Basic type for real space grid methods
!> \par History
!>      JGH (22-May-2002) : New routine rs_grid_zero
!>      JGH (12-Jun-2002) : Bug fix for mpi groups
!>      JGH (19-Jun-2003) : Added routine for task distribution
!>      JGH (23-Nov-2003) : Added routine for task loop separation
!> \author JGH (18-Mar-2001)
! **************************************************************************************************
MODULE realspace_grid_types
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE fast,                            ONLY: copy_cr,&
                                              copy_rc
   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE machine,                         ONLY: m_memory
   USE mathlib,                         ONLY: det_3x3
   USE message_passing,                 ONLY: mp_comm_null,&
                                              mp_comm_type,&
                                              mp_request_null,&
                                              mp_request_type,&
                                              mp_waitall,&
                                              mp_waitany
   USE offload_api,                     ONLY: offload_buffer_type,&
                                              offload_create_buffer,&
                                              offload_free_buffer
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL,&
                                              pw_grid_type
   USE pw_grids,                        ONLY: pw_grid_release,&
                                              pw_grid_retain
   USE pw_methods,                      ONLY: pw_integrate_function
   USE pw_types,                        ONLY: COMPLEXDATA3D,&
                                              REALDATA3D,&
                                              pw_type
   USE util,                            ONLY: get_limit

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: realspace_grid_type, &
             realspace_grid_desc_type, &
             realspace_grid_p_type, &
             realspace_grid_desc_p_type, &
             realspace_grid_input_type

   PUBLIC :: transfer_rs2pw, &
             transfer_pw2rs, &
             rs_grid_zero, &
             rs_grid_set_box, &
             rs_grid_create, &
             rs_grid_create_descriptor, &
             rs_grid_retain_descriptor, &
             rs_grid_release, &
             rs_grid_release_descriptor, &
             rs_grid_reorder_ranks, &
             rs_grid_print, &
             rs_grid_locate_rank, &
             rs_grid_max_ngpts, &
             rs_grid_mult_and_add, &
             map_gaussian_here

   INTEGER, PARAMETER, PUBLIC               :: rsgrid_distributed = 0, &
                                               rsgrid_replicated = 1, &
                                               rsgrid_automatic = 2

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'realspace_grid_types'

! **************************************************************************************************
   TYPE realspace_grid_input_type
      INTEGER       :: distribution_type = rsgrid_replicated
      INTEGER       :: distribution_layout(3) = -1
      REAL(KIND=dp) :: memory_factor = 0.0_dp
      LOGICAL       :: lock_distribution = .FALSE.
      INTEGER       :: nsmax = -1
      REAL(KIND=dp) :: halo_reduction_factor = 1.0_dp
   END TYPE realspace_grid_input_type

! **************************************************************************************************
   TYPE realspace_grid_desc_type
      TYPE(pw_grid_type), POINTER   :: pw => NULL() ! the pw grid

      INTEGER :: ref_count = 0 ! reference count

      INTEGER(int_8) :: ngpts = 0_int_8 ! # grid points
      INTEGER, DIMENSION(3) :: npts = 0 ! # grid points per dimension
      INTEGER, DIMENSION(3) :: lb = 0 ! lower bounds
      INTEGER, DIMENSION(3) :: ub = 0 ! upper bounds

      INTEGER :: border = 0 ! border points

      INTEGER, DIMENSION(3) :: perd = -1 ! periodicity enforced
      REAL(KIND=dp), DIMENSION(3, 3) :: dh = 0.0_dp ! incremental grid matrix
      REAL(KIND=dp), DIMENSION(3, 3) :: dh_inv = 0.0_dp ! inverse incremental grid matrix
      LOGICAL :: orthorhombic = .TRUE. ! grid symmetry

      LOGICAL :: parallel = .TRUE. ! whether the corresponding pw grid is distributed
      LOGICAL :: distributed = .TRUE. ! whether the rs grid is distributed
      ! these MPI related quantities are only meaningful depending on how the grid has been laid out
      ! they are most useful for fully distributed grids, where they reflect the topology of the grid
      TYPE(mp_comm_type) :: group = mp_comm_null
      INTEGER :: my_pos = -1
      LOGICAL :: group_head = .FALSE.
      INTEGER :: group_size = 0
      INTEGER, DIMENSION(3) :: group_dim = -1
      INTEGER, DIMENSION(3) :: group_coor = -1
      INTEGER, DIMENSION(3) :: neighbours = -1
      ! only meaningful on distributed grids
      ! a list of bounds for each CPU
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: lb_global
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: ub_global
      ! a mapping from linear rank to 3d coord
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: rank2coord
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: coord2rank
      ! a mapping from index to rank (which allows to figure out easily on which rank a given point of the grid is)
      INTEGER, DIMENSION(:), ALLOCATABLE :: x2coord
      INTEGER, DIMENSION(:), ALLOCATABLE :: y2coord
      INTEGER, DIMENSION(:), ALLOCATABLE :: z2coord

      INTEGER                :: my_virtual_pos = -1
      INTEGER, DIMENSION(3) :: virtual_group_coor = -1

      INTEGER, DIMENSION(:), ALLOCATABLE :: virtual2real, real2virtual

   END TYPE realspace_grid_desc_type

   TYPE realspace_grid_type

      TYPE(realspace_grid_desc_type), POINTER :: desc => NULL()

      INTEGER :: ngpts_local = -1 ! local dimensions
      INTEGER, DIMENSION(3) :: npts_local = -1
      INTEGER, DIMENSION(3) :: lb_local = -1
      INTEGER, DIMENSION(3) :: ub_local = -1
      INTEGER, DIMENSION(3) :: lb_real = -1 ! lower bounds of the real local data
      INTEGER, DIMENSION(3) :: ub_real = -1 ! upper bounds of the real local data

      INTEGER, DIMENSION(:), ALLOCATABLE         :: px, py, pz ! index translators
      TYPE(offload_buffer_type)                  :: buffer = offload_buffer_type() ! owner of the grid's memory
      REAL(KIND=dp), DIMENSION(:, :, :), CONTIGUOUS, POINTER :: r => NULL() ! the grid (pointer to buffer%host_buffer)

   END TYPE realspace_grid_type

! **************************************************************************************************
   TYPE realspace_grid_p_type
      TYPE(realspace_grid_type), POINTER :: rs_grid => NULL()
   END TYPE realspace_grid_p_type

   TYPE realspace_grid_desc_p_type
      TYPE(realspace_grid_desc_type), POINTER :: rs_desc => NULL()
   END TYPE realspace_grid_desc_p_type

CONTAINS

! **************************************************************************************************
!> \brief returns the 1D rank of the task which is a cartesian shift away from 1D rank rank_in
!>        only possible if rs_grid is a distributed grid
!> \param rs_desc ...
!> \param rank_in ...
!> \param shift ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION rs_grid_locate_rank(rs_desc, rank_in, shift) RESULT(rank_out)
      TYPE(realspace_grid_desc_type), INTENT(IN)         :: rs_desc
      INTEGER, INTENT(IN)                                :: rank_in
      INTEGER, DIMENSION(3), INTENT(IN)                  :: shift
      INTEGER                                            :: rank_out

      INTEGER                                            :: coord(3)

      coord = MODULO(rs_desc%rank2coord(:, rank_in) + shift, rs_desc%group_dim)
      rank_out = rs_desc%coord2rank(coord(1), coord(2), coord(3))
   END FUNCTION rs_grid_locate_rank

! **************************************************************************************************
!> \brief Determine the setup of real space grids - this is divided up into the
!>        creation of a descriptor and the actual grid itself (see rs_grid_create)
!> \param desc ...
!> \param pw_grid ...
!> \param input_settings ...
!> \param border_points ...
!> \par History
!>      JGH (08-Jun-2003) : nsmax <= 0 indicates fully replicated grid
!>      Iain Bethune (05-Sep-2008) : modified cut heuristic
!>      (c) The Numerical Algorithms Group (NAG) Ltd, 2008 on behalf of the HECToR project
!>      - Create a descriptor for realspace grids with a number of border
!>        points as exactly given by the optional argument border_points.
!>        These grids are always distributed.
!>        (27.11.2013, Matthias Krack)
!> \author JGH (18-Mar-2001)
! **************************************************************************************************
   SUBROUTINE rs_grid_create_descriptor(desc, pw_grid, input_settings, border_points)
      TYPE(realspace_grid_desc_type), POINTER            :: desc
      TYPE(pw_grid_type), INTENT(INOUT), TARGET          :: pw_grid
      TYPE(realspace_grid_input_type), INTENT(IN)        :: input_settings
      INTEGER, INTENT(IN), OPTIONAL                      :: border_points

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rs_grid_create_descriptor'

      INTEGER                                            :: border_size, dir, handle, i, j, k, l, &
                                                            lb(2), min_npts_real, n_slices(3), &
                                                            n_slices_tmp(3), nmin
      LOGICAL                                            :: overlap
      REAL(KIND=dp)                                      :: ratio, ratio_best, volume, volume_dist

      CALL timeset(routineN, handle)

      IF (PRESENT(border_points)) THEN
         border_size = border_points
      ELSE
         border_size = 0
      END IF

      ALLOCATE (desc)

      CALL pw_grid%para%group%sync()

      desc%pw => pw_grid
      CALL pw_grid_retain(desc%pw)

      desc%dh = pw_grid%dh
      desc%dh_inv = pw_grid%dh_inv
      desc%orthorhombic = pw_grid%orthorhombic
      desc%ref_count = 1

      IF (pw_grid%para%mode == PW_MODE_LOCAL) THEN
         ! The corresponding group has dimension 1
         ! All operations will be done locally
         desc%npts = pw_grid%npts
         desc%ngpts = PRODUCT(INT(desc%npts, KIND=int_8))
         desc%lb = pw_grid%bounds(1, :)
         desc%ub = pw_grid%bounds(2, :)
         desc%border = border_size
         IF (border_size == 0) THEN
            desc%perd = 1
         ELSE
            desc%perd = 0
         END IF
         desc%parallel = .FALSE.
         desc%distributed = .FALSE.
         desc%group = mp_comm_null
         desc%group_size = 1
         desc%group_head = .TRUE.
         desc%group_dim = 1
         desc%group_coor = 0
         desc%my_pos = 0
      ELSE
         ! group size of desc grid
         ! global grid dimensions are still the same
         desc%group_size = pw_grid%para%group_size
         desc%npts = pw_grid%npts
         desc%ngpts = PRODUCT(INT(desc%npts, KIND=int_8))
         desc%lb = pw_grid%bounds(1, :)
         desc%ub = pw_grid%bounds(2, :)

         ! this is the eventual border size
         IF (border_size == 0) THEN
            nmin = (input_settings%nsmax + 1)/2
            nmin = MAX(0, NINT(nmin*input_settings%halo_reduction_factor))
         ELSE
            ! Set explicitly the requested border size
            nmin = border_size
         END IF

         IF (input_settings%distribution_type == rsgrid_replicated) THEN

            n_slices = 1
            IF (border_size > 0) THEN
               CALL cp_abort(__LOCATION__, &
                             "An explicit border size > 0 is not yet working for "// &
                             "replicated realspace grids. Request DISTRIBUTION_TYPE "// &
                             "distributed for RS_GRID explicitly.")
            END IF

         ELSE

            n_slices = 1
            ratio_best = -HUGE(ratio_best)

            ! don't allow distributions with more processors than real grid points
            DO k = 1, MIN(desc%npts(3), desc%group_size)
            DO j = 1, MIN(desc%npts(2), desc%group_size)
               i = MIN(desc%npts(1), desc%group_size/(j*k))
               n_slices_tmp = (/i, j, k/)

               ! we don't match the actual number of CPUs
               IF (PRODUCT(n_slices_tmp) .NE. desc%group_size) CYCLE

               ! we see if there has been a input constraint
               ! i.e. if the layout is not -1 we need to fullfil it
               IF (.NOT. ALL(PACK(n_slices_tmp == input_settings%distribution_layout, &
                                  (/-1, -1, -1/) /= input_settings%distribution_layout) &
                             )) CYCLE

               ! We can not work with a grid that has more local than global grid points.
               ! This can happen when a halo region wraps around and overlaps with the other halo.
               overlap = .FALSE.
               DO dir = 1, 3
                  IF (n_slices_tmp(dir) > 1) THEN
                     DO l = 0, n_slices_tmp(dir) - 1
                        lb = get_limit(desc%npts(dir), n_slices_tmp(dir), l)
                        IF (lb(2) - lb(1) + 1 + 2*nmin > desc%npts(dir)) overlap = .TRUE.
                     END DO
                  END IF
               END DO
               IF (overlap) CYCLE

               ! a heuristic optimisation to reduce the memory usage
               ! we go for the smallest local to real volume
               ! volume of the box without the wings / volume of the box with the wings
               ! with prefactodesc to promote less cuts in Z dimension
               ratio = PRODUCT(REAL(desc%npts, KIND=dp)/n_slices_tmp)/ &
                       PRODUCT(REAL(desc%npts, KIND=dp)/n_slices_tmp + &
                               MERGE((/0.0, 0.0, 0.0/), 2*(/1.06*nmin, 1.05*nmin, 1.03*nmin/), n_slices_tmp == (/1, 1, 1/)))
               IF (ratio > ratio_best) THEN
                  ratio_best = ratio
                  n_slices = n_slices_tmp
               END IF

            END DO
            END DO

            ! if automatic we can still decide this is a replicated grid
            ! if the memory gain (or the gain is messages) is too small.
            IF (input_settings%distribution_type == rsgrid_automatic) THEN
               volume = PRODUCT(REAL(desc%npts, KIND=dp))
               volume_dist = PRODUCT(REAL(desc%npts, KIND=dp)/n_slices + &
                                     MERGE((/0, 0, 0/), 2*(/nmin, nmin, nmin/), n_slices == (/1, 1, 1/)))
               IF (volume < volume_dist*input_settings%memory_factor) THEN
                  n_slices = 1
               END IF
            END IF

         END IF

         desc%group_dim(:) = n_slices(:)
         CALL desc%group%from_dup(pw_grid%para%group)
         desc%group_size = desc%group%num_pe
         desc%my_pos = desc%group%mepos

         IF (ALL(n_slices == 1)) THEN
            ! CASE 1 : only one slice: we do not need overlapping regions and special
            !          recombination of the total density
            desc%border = border_size
            IF (border_size == 0) THEN
               desc%perd = 1
            ELSE
               desc%perd = 0
            END IF
            desc%distributed = .FALSE.
            desc%parallel = .TRUE.
            desc%group_head = pw_grid%para%group_head
            desc%group_coor(:) = 0
            desc%my_virtual_pos = 0

            ALLOCATE (desc%virtual2real(0:desc%group_size - 1))
            ALLOCATE (desc%real2virtual(0:desc%group_size - 1))
            ! Start with no reordering
            DO i = 0, desc%group_size - 1
               desc%virtual2real(i) = i
               desc%real2virtual(i) = i
            END DO
         ELSE
            ! CASE 2 : general case
            ! periodicity is no longer enforced arbritary directions
            IF (border_size == 0) THEN
               desc%perd = 1
               DO dir = 1, 3
                  IF (n_slices(dir) > 1) desc%perd(dir) = 0
               END DO
            ELSE
               desc%perd(:) = 0
            END IF
            ! we keep a border of nmin points
            desc%border = nmin
            ! we are going parallel on the real space grid
            desc%parallel = .TRUE.
            desc%distributed = .TRUE.
            desc%group_head = (desc%my_pos == 0)

            ! set up global info about the distribution
            ALLOCATE (desc%rank2coord(3, 0:desc%group_size - 1))
            ALLOCATE (desc%coord2rank(0:desc%group_dim(1) - 1, 0:desc%group_dim(2) - 1, 0:desc%group_dim(3) - 1))
            ALLOCATE (desc%lb_global(3, 0:desc%group_size - 1))
            ALLOCATE (desc%ub_global(3, 0:desc%group_size - 1))
            ALLOCATE (desc%x2coord(desc%lb(1):desc%ub(1)))
            ALLOCATE (desc%y2coord(desc%lb(2):desc%ub(2)))
            ALLOCATE (desc%z2coord(desc%lb(3):desc%ub(3)))

            DO i = 0, desc%group_size - 1
               ! Calculate coordinates in a row-major order (to be SMP-friendly)
               desc%rank2coord(1, i) = i/(desc%group_dim(2)*desc%group_dim(3))
               desc%rank2coord(2, i) = MODULO(i, desc%group_dim(2)*desc%group_dim(3)) &
                                       /desc%group_dim(3)
               desc%rank2coord(3, i) = MODULO(i, desc%group_dim(3))

               IF (i == desc%my_pos) THEN
                  desc%group_coor = desc%rank2coord(:, i)
               END IF

               desc%coord2rank(desc%rank2coord(1, i), desc%rank2coord(2, i), desc%rank2coord(3, i)) = i
               ! the lb_global and ub_global correspond to lb_real and ub_real of each task
               desc%lb_global(:, i) = desc%lb
               desc%ub_global(:, i) = desc%ub
               DO dir = 1, 3
                  IF (desc%group_dim(dir) .GT. 1) THEN
                     lb = get_limit(desc%npts(dir), desc%group_dim(dir), desc%rank2coord(dir, i))
                     desc%lb_global(dir, i) = lb(1) + desc%lb(dir) - 1
                     desc%ub_global(dir, i) = lb(2) + desc%lb(dir) - 1
                  END IF
               END DO
            END DO

            ! map a grid point to a CPU coord
            DO dir = 1, 3
               DO l = 0, desc%group_dim(dir) - 1
                  IF (desc%group_dim(dir) .GT. 1) THEN
                     lb = get_limit(desc%npts(dir), desc%group_dim(dir), l)
                     lb = lb + desc%lb(dir) - 1
                  ELSE
                     lb(1) = desc%lb(dir)
                     lb(2) = desc%ub(dir)
                  END IF
                  SELECT CASE (dir)
                  CASE (1)
                     desc%x2coord(lb(1):lb(2)) = l
                  CASE (2)
                     desc%y2coord(lb(1):lb(2)) = l
                  CASE (3)
                     desc%z2coord(lb(1):lb(2)) = l
                  END SELECT
               END DO
            END DO

            ! an upper bound for the number of neighbours the border is overlapping with
            DO dir = 1, 3
               desc%neighbours(dir) = 0
               IF ((n_slices(dir) > 1) .OR. (border_size > 0)) THEN
                  min_npts_real = HUGE(0)
                  DO l = 0, n_slices(dir) - 1
                     lb = get_limit(desc%npts(dir), n_slices(dir), l)
                     min_npts_real = MIN(lb(2) - lb(1) + 1, min_npts_real)
                  END DO
                  desc%neighbours(dir) = (desc%border + min_npts_real - 1)/min_npts_real
               END IF
            END DO

            ALLOCATE (desc%virtual2real(0:desc%group_size - 1))
            ALLOCATE (desc%real2virtual(0:desc%group_size - 1))
            ! Start with no reordering
            DO i = 0, desc%group_size - 1
               desc%virtual2real(i) = i
               desc%real2virtual(i) = i
            END DO

            desc%my_virtual_pos = desc%real2virtual(desc%my_pos)
            desc%virtual_group_coor(:) = desc%rank2coord(:, desc%my_virtual_pos)

         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE rs_grid_create_descriptor

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \param desc ...
! **************************************************************************************************
   SUBROUTINE rs_grid_create(rs, desc)
      TYPE(realspace_grid_type), INTENT(OUT)             :: rs
      TYPE(realspace_grid_desc_type), INTENT(INOUT), &
         TARGET                                          :: desc

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rs_grid_create'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      rs%desc => desc
      CALL rs_grid_retain_descriptor(rs%desc)

      IF (desc%pw%para%mode == PW_MODE_LOCAL) THEN
         ! The corresponding group has dimension 1
         ! All operations will be done locally
         rs%lb_real = desc%lb
         rs%ub_real = desc%ub
         rs%lb_local = rs%lb_real - desc%border*(1 - desc%perd)
         rs%ub_local = rs%ub_real + desc%border*(1 - desc%perd)
         rs%npts_local = rs%ub_local - rs%lb_local + 1
         rs%ngpts_local = PRODUCT(rs%npts_local)
      END IF

      IF (ALL(rs%desc%group_dim == 1)) THEN
         ! CASE 1 : only one slice: we do not need overlapping regions and special
         !          recombination of the total density
         rs%lb_real = desc%lb
         rs%ub_real = desc%ub
         rs%lb_local = rs%lb_real - desc%border*(1 - desc%perd)
         rs%ub_local = rs%ub_real + desc%border*(1 - desc%perd)
         rs%npts_local = rs%ub_local - rs%lb_local + 1
         rs%ngpts_local = PRODUCT(rs%npts_local)
      ELSE
         ! CASE 2 : general case
         ! extract some more derived quantities about the local grid
         rs%lb_real = desc%lb_global(:, desc%my_virtual_pos)
         rs%ub_real = desc%ub_global(:, desc%my_virtual_pos)
         rs%lb_local = rs%lb_real - desc%border*(1 - desc%perd)
         rs%ub_local = rs%ub_real + desc%border*(1 - desc%perd)
         rs%npts_local = rs%ub_local - rs%lb_local + 1
         rs%ngpts_local = PRODUCT(rs%npts_local)
      END IF

      CALL offload_create_buffer(rs%ngpts_local, rs%buffer)
      rs%r(rs%lb_local(1):rs%ub_local(1), &
           rs%lb_local(2):rs%ub_local(2), &
           rs%lb_local(3):rs%ub_local(3)) => rs%buffer%host_buffer

      ALLOCATE (rs%px(desc%npts(1)))
      ALLOCATE (rs%py(desc%npts(2)))
      ALLOCATE (rs%pz(desc%npts(3)))

      CALL timestop(handle)

   END SUBROUTINE rs_grid_create

! **************************************************************************************************
!> \brief Defines a new ordering of ranks on this realspace grid, recalculating
!>        the data bounds and reallocating the grid.  As a result, each MPI process
!>        now has a real rank (i.e., its rank in the MPI communicator from the pw grid)
!>        and a virtual rank (the rank of the process where the data now owned by this
!>        process would reside in an ordinary cartesian distribution).
!>        NB. Since the grid size required may change, the caller should be sure to release
!>        and recreate the corresponding rs_grids
!>        The desc%real2virtual and desc%virtual2real arrays can be used to map
!>        a physical rank to the 'rank' of data owned by that process and vice versa
!> \param desc ...
!> \param real2virtual ...
!> \par History
!>        04-2009 created [Iain Bethune]
!>          (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
! **************************************************************************************************
   PURE SUBROUTINE rs_grid_reorder_ranks(desc, real2virtual)

      TYPE(realspace_grid_desc_type), INTENT(INOUT)      :: desc
      INTEGER, DIMENSION(:), INTENT(IN)                  :: real2virtual

      INTEGER                                            :: i

      desc%real2virtual(:) = real2virtual

      DO i = 0, desc%group_size - 1
         desc%virtual2real(desc%real2virtual(i)) = i
      END DO

      desc%my_virtual_pos = desc%real2virtual(desc%my_pos)

      IF (.NOT. ALL(desc%group_dim == 1)) THEN
         desc%virtual_group_coor(:) = desc%rank2coord(:, desc%my_virtual_pos)
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief Print information on grids to output
!> \param rs ...
!> \param iounit ...
!> \author JGH (17-May-2007)
! **************************************************************************************************
   SUBROUTINE rs_grid_print(rs, iounit)
      TYPE(realspace_grid_type), INTENT(IN)              :: rs
      INTEGER, INTENT(in)                                :: iounit

      INTEGER                                            :: dir, i, nn
      REAL(KIND=dp)                                      :: pp(3)

      IF (rs%desc%parallel) THEN
         IF (iounit > 0) THEN
            WRITE (iounit, '(/,A,T71,I10)') &
               " RS_GRID| Information for grid number ", rs%desc%pw%id_nr
            DO i = 1, 3
               WRITE (iounit, '(A,I3,T30,2I8,T62,A,T71,I10)') " RS_GRID|   Bounds ", &
                  i, rs%desc%lb(i), rs%desc%ub(i), "Points:", rs%desc%npts(i)
            END DO
            IF (.NOT. rs%desc%distributed) THEN
               WRITE (iounit, '(A)') " RS_GRID| Real space fully replicated"
               WRITE (iounit, '(A,T71,I10)') &
                  " RS_GRID| Group size ", rs%desc%group_dim(2)
            ELSE
               DO dir = 1, 3
                  IF (rs%desc%perd(dir) /= 1) THEN
                     WRITE (iounit, '(A,T71,I3,A)') &
                        " RS_GRID| Real space distribution over ", rs%desc%group_dim(dir), " groups"
                     WRITE (iounit, '(A,T71,I10)') &
                        " RS_GRID| Real space distribution along direction ", dir
                     WRITE (iounit, '(A,T71,I10)') &
                        " RS_GRID| Border size ", rs%desc%border
                  END IF
               END DO
            END IF
         END IF
         IF (rs%desc%distributed) THEN
            DO dir = 1, 3
               IF (rs%desc%perd(dir) /= 1) THEN
                  nn = rs%npts_local(dir)
                  CALL rs%desc%group%sum(nn)
                  pp(1) = REAL(nn, KIND=dp)/REAL(PRODUCT(rs%desc%group_dim), KIND=dp)
                  nn = rs%npts_local(dir)
                  CALL rs%desc%group%max(nn)
                  pp(2) = REAL(nn, KIND=dp)
                  nn = rs%npts_local(dir)
                  CALL rs%desc%group%min(nn)
                  pp(3) = REAL(nn, KIND=dp)
                  IF (iounit > 0) THEN
                     WRITE (iounit, '(A,T48,A)') " RS_GRID|   Distribution", &
                        "  Average         Max         Min"
                     WRITE (iounit, '(A,T45,F12.1,2I12)') " RS_GRID|   Planes   ", &
                        pp(1), NINT(pp(2)), NINT(pp(3))
                  END IF
               END IF
            END DO
!          WRITE ( iounit, '(/)' )
         END IF
      ELSE
         IF (iounit > 0) THEN
            WRITE (iounit, '(/,A,T71,I10)') &
               " RS_GRID| Information for grid number ", rs%desc%pw%id_nr
            DO i = 1, 3
               WRITE (iounit, '(A,I3,T30,2I8,T62,A,T71,I10)') " RS_GRID|   Bounds ", &
                  i, rs%desc%lb(i), rs%desc%ub(i), "Points:", rs%desc%npts(i)
            END DO
!         WRITE ( iounit, '(/)' )
         END IF
      END IF

   END SUBROUTINE rs_grid_print

! **************************************************************************************************
!> \brief Copy a function to a PW grid type from a real space type
!> \param rs ...
!> \param pw ...
! **************************************************************************************************
   SUBROUTINE transfer_rs2pw(rs, pw)

      TYPE(realspace_grid_type), INTENT(IN)              :: rs
      TYPE(pw_type), INTENT(IN)                          :: pw

      CHARACTER(len=*), PARAMETER                        :: routineN = 'transfer_rs2pw'

      INTEGER                                            :: handle, handle2, i

      CALL timeset(routineN, handle2)
      CALL timeset(routineN//"_"//TRIM(ADJUSTL(cp_to_string(CEILING(pw%pw_grid%cutoff/10)*10))), handle)

      IF (.NOT. ASSOCIATED(rs%desc%pw, pw%pw_grid)) &
         CPABORT("Different rs and pw indentifiers")
      IF ((pw%in_use .NE. REALDATA3D) .AND. (pw%in_use .NE. COMPLEXDATA3D)) &
         CPABORT("Need REALDATA3D or COMPLEXDATA3D")

      IF (rs%desc%distributed) THEN
         CALL transfer_rs2pw_distributed(rs, pw)
      ELSE IF (rs%desc%parallel) THEN
         CALL transfer_rs2pw_replicated(rs, pw)
      ELSE ! treat simple serial case locally
         IF (pw%in_use == REALDATA3D) THEN
            IF (rs%desc%border == 0) THEN
               CALL dcopy(SIZE(rs%r), rs%r, 1, pw%cr3d, 1)
            ELSE
               CPASSERT(LBOUND(pw%cr3d, 3) .EQ. rs%lb_real(3))
!$OMP          PARALLEL DO DEFAULT(NONE) SHARED(pw,rs)
               DO i = rs%lb_real(3), rs%ub_real(3)
                  pw%cr3d(:, :, i) = rs%r(rs%lb_real(1):rs%ub_real(1), &
                                          rs%lb_real(2):rs%ub_real(2), i)
               END DO
!$OMP          END PARALLEL DO
            END IF
         ELSEIF (pw%in_use == COMPLEXDATA3D) THEN
            IF (rs%desc%border == 0) THEN
               CALL copy_rc(rs%r, pw%cc3d)
            ELSE
               CPASSERT(LBOUND(pw%cr3d, 3) .EQ. rs%lb_real(3))
!$OMP          PARALLEL DO DEFAULT(NONE) SHARED(pw,rs)
               DO i = rs%lb_real(3), rs%ub_real(3)
                  pw%cc3d(:, :, i) = CMPLX(rs%r(rs%lb_real(1):rs%ub_real(1), &
                                                rs%lb_real(2):rs%ub_real(2), i), &
                                           0.0_dp, KIND=dp)
               END DO
!$OMP          END PARALLEL DO
            END IF
         ELSE
            CPABORT("PW type not compatible")
         END IF
      END IF

      CALL timestop(handle)
      CALL timestop(handle2)

   END SUBROUTINE transfer_rs2pw

! **************************************************************************************************
!> \brief Copy a function from/to a PW grid type to/from a real space type
!> \param rs ...
!> \param pw ...
! **************************************************************************************************
   SUBROUTINE transfer_pw2rs(rs, pw)

      TYPE(realspace_grid_type), INTENT(IN)              :: rs
      TYPE(pw_type), INTENT(IN)                          :: pw

      CHARACTER(len=*), PARAMETER                        :: routineN = 'transfer_pw2rs'

      INTEGER                                            :: handle, handle2, i, im, j, jm, k, km

      CALL timeset(routineN, handle2)
      CALL timeset(routineN//"_"//TRIM(ADJUSTL(cp_to_string(CEILING(pw%pw_grid%cutoff/10)*10))), handle)

      IF (.NOT. ASSOCIATED(rs%desc%pw, pw%pw_grid)) &
         CPABORT("Different rs and pw indentifiers")
      IF ((pw%in_use .NE. REALDATA3D) .AND. (pw%in_use .NE. COMPLEXDATA3D)) &
         CPABORT("Need REALDATA3D or COMPLEXDATA3D")

      IF (rs%desc%distributed) THEN
         CALL transfer_pw2rs_distributed(rs, pw)
      ELSE IF (rs%desc%parallel) THEN
         CALL transfer_pw2rs_replicated(rs, pw)
      ELSE ! treat simple serial case locally
         IF (pw%in_use == REALDATA3D) THEN
            IF (rs%desc%border == 0) THEN
               CALL dcopy(SIZE(rs%r), pw%cr3d, 1, rs%r, 1)
            ELSE
!$OMP          PARALLEL DO DEFAULT(NONE) &
!$OMP                      PRIVATE(i,im,j,jm,k,km) &
!$OMP                      SHARED(pw,rs)
               DO k = rs%lb_local(3), rs%ub_local(3)
                  IF (k < rs%lb_real(3)) THEN
                     km = k + rs%desc%npts(3)
                  ELSE IF (k > rs%ub_real(3)) THEN
                     km = k - rs%desc%npts(3)
                  ELSE
                     km = k
                  END IF
                  DO j = rs%lb_local(2), rs%ub_local(2)
                     IF (j < rs%lb_real(2)) THEN
                        jm = j + rs%desc%npts(2)
                     ELSE IF (j > rs%ub_real(2)) THEN
                        jm = j - rs%desc%npts(2)
                     ELSE
                        jm = j
                     END IF
                     DO i = rs%lb_local(1), rs%ub_local(1)
                        IF (i < rs%lb_real(1)) THEN
                           im = i + rs%desc%npts(1)
                        ELSE IF (i > rs%ub_real(1)) THEN
                           im = i - rs%desc%npts(1)
                        ELSE
                           im = i
                        END IF
                        rs%r(i, j, k) = pw%cr3d(im, jm, km)
                     END DO
                  END DO
               END DO
!$OMP          END PARALLEL DO
            END IF
         ELSEIF (pw%in_use == COMPLEXDATA3D) THEN
            IF (rs%desc%border == 0) THEN
               CALL copy_cr(pw%cc3d, rs%r)
            ELSE
!$OMP          PARALLEL DO DEFAULT(NONE) &
!$OMP                      PRIVATE(i,im,j,jm,k,km) &
!$OMP                      SHARED(pw,rs)
               DO k = rs%lb_local(3), rs%ub_local(3)
                  IF (k < rs%lb_real(3)) THEN
                     km = k + rs%desc%npts(3)
                  ELSE IF (k > rs%ub_real(3)) THEN
                     km = k - rs%desc%npts(3)
                  ELSE
                     km = k
                  END IF
                  DO j = rs%lb_local(2), rs%ub_local(2)
                     IF (j < rs%lb_real(2)) THEN
                        jm = j + rs%desc%npts(2)
                     ELSE IF (j > rs%ub_real(2)) THEN
                        jm = j - rs%desc%npts(2)
                     ELSE
                        jm = j
                     END IF
                     DO i = rs%lb_local(1), rs%ub_local(1)
                        IF (i < rs%lb_real(1)) THEN
                           im = i + rs%desc%npts(1)
                        ELSE IF (i > rs%ub_real(1)) THEN
                           im = i - rs%desc%npts(1)
                        ELSE
                           im = i
                        END IF
                        rs%r(i, j, k) = REAL(pw%cc3d(im, jm, km), KIND=dp)
                     END DO
                  END DO
               END DO
!$OMP          END PARALLEL DO
            END IF
         ELSE
            CPABORT("PW type not compatible")
         END IF
      END IF

      CALL timestop(handle)
      CALL timestop(handle2)

   END SUBROUTINE transfer_pw2rs

! **************************************************************************************************
!> \brief transfer from a realspace grid to a planewave grid
!> \param rs ...
!> \param pw ...
! **************************************************************************************************
   SUBROUTINE transfer_rs2pw_replicated(rs, pw)
      TYPE(realspace_grid_type), INTENT(IN)              :: rs
      TYPE(pw_type), INTENT(IN)                          :: pw

      INTEGER                                            :: dest, ii, ip, ix, iy, iz, nma, nn, s(3), &
                                                            source
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rcount
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: recvbuf, sendbuf, swaparray

      ASSOCIATE (np => pw%pw_grid%para%group_size, bo => pw%pw_grid%para%bo(1:2, 1:3, 0:pw%pw_grid%para%group_size - 1, 1), &
                 pbo => pw%pw_grid%bounds, group => pw%pw_grid%para%rs_group, mepos => pw%pw_grid%para%rs_mpo, &
                 grid => rs%r)
         ALLOCATE (rcount(0:np - 1))
         DO ip = 1, np
            rcount(ip - 1) = PRODUCT(bo(2, :, ip) - bo(1, :, ip) + 1)
         END DO
         nma = MAXVAL(rcount(0:np - 1))
         ALLOCATE (sendbuf(nma), recvbuf(nma))
         sendbuf = 1.0E99_dp; recvbuf = 1.0E99_dp ! init mpi'ed buffers to silence warnings under valgrind

         !sample peak memory
         CALL m_memory()

         dest = MODULO(mepos + 1, np)
         source = MODULO(mepos - 1, np)
         sendbuf = 0.0_dp

         DO ip = 1, np

            lb = pbo(1, :) + bo(1, :, MODULO(mepos - ip, np) + 1) - 1
            ub = pbo(1, :) + bo(2, :, MODULO(mepos - ip, np) + 1) - 1
            ! this loop takes about the same time as the message passing call
            ! notice that the range of ix is only a small fraction of the first index of grid
            ! therefore it seems faster to have the second index as the innermost loop
            ! if this runs on many cpus
            ! tested on itanium, pentium4, opteron, ultrasparc...
            s = ub - lb + 1
            DO iz = lb(3), ub(3)
               DO ix = lb(1), ub(1)
                  ii = (iz - lb(3))*s(1)*s(2) + (ix - lb(1)) + 1
                  DO iy = lb(2), ub(2)
                     sendbuf(ii) = sendbuf(ii) + grid(ix, iy, iz)
                     ii = ii + s(1)
                  END DO
               END DO
            END DO
            IF (ip .EQ. np) EXIT
            CALL group%sendrecv(sendbuf, dest, recvbuf, source, 13)
            CALL MOVE_ALLOC(sendbuf, swaparray)
            CALL MOVE_ALLOC(recvbuf, sendbuf)
            CALL MOVE_ALLOC(swaparray, recvbuf)
         END DO
         nn = rcount(mepos)
      END ASSOCIATE

      IF (pw%in_use == REALDATA3D) THEN
         CALL dcopy(nn, sendbuf, 1, pw%cr3d, 1)
      ELSEIF (pw%in_use == COMPLEXDATA3D) THEN
         CALL copy_rc(RESHAPE(sendbuf, SHAPE(pw%cc3d)), pw%cc3d)
      ELSE
         CPABORT("PW type not compatible")
      END IF

      DEALLOCATE (rcount)
      DEALLOCATE (sendbuf)
      DEALLOCATE (recvbuf)

   END SUBROUTINE transfer_rs2pw_replicated

! **************************************************************************************************
!> \brief transfer from a planewave grid to a realspace grid
!> \param rs ...
!> \param pw ...
! **************************************************************************************************
   SUBROUTINE transfer_pw2rs_replicated(rs, pw)
      TYPE(realspace_grid_type), INTENT(IN)              :: rs
      TYPE(pw_type), INTENT(IN)                          :: pw

      INTEGER                                            :: dest, i, ii, im, ip, ix, iy, iz, j, jm, &
                                                            k, km, nma, nn, source
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rcount
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: recvbuf, sendbuf, swaparray
      TYPE(mp_request_type), DIMENSION(2)                :: req

      ASSOCIATE (np => pw%pw_grid%para%group_size, bo => pw%pw_grid%para%bo(1:2, 1:3, 0:pw%pw_grid%para%group_size - 1, 1), &
                 pbo => pw%pw_grid%bounds, group => pw%pw_grid%para%rs_group, mepos => pw%pw_grid%para%rs_mpo, &
                 grid => rs%r)
         ALLOCATE (rcount(0:np - 1))
         DO ip = 1, np
            rcount(ip - 1) = PRODUCT(bo(2, :, ip) - bo(1, :, ip) + 1)
         END DO
         nma = MAXVAL(rcount(0:np - 1))
         ALLOCATE (sendbuf(nma), recvbuf(nma))
         sendbuf = 1.0E99_dp; recvbuf = 1.0E99_dp ! init mpi'ed buffers to silence warnings under valgrind

         !sample peak memory
         CALL m_memory()

         nn = rcount(mepos)
         IF (pw%in_use == REALDATA3D) THEN
            CALL dcopy(nn, pw%cr3d, 1, sendbuf, 1)
         ELSEIF (pw%in_use == COMPLEXDATA3D) THEN
            CALL dcopy(nn, TRANSFER(pw%cc3d, sendbuf), 2, sendbuf, 1)
         ELSE
            CPABORT("PW type not compatible")
         END IF

         dest = MODULO(mepos + 1, np)
         source = MODULO(mepos - 1, np)

         DO ip = 0, np - 1
            ! we must shift the buffer only np-1 times around
            IF (ip .NE. np - 1) THEN
               CALL group%isendrecv(sendbuf, dest, recvbuf, source, &
                                    req(1), req(2), 13)
            END IF
            lb = pbo(1, :) + bo(1, :, MODULO(mepos - ip, np) + 1) - 1
            ub = pbo(1, :) + bo(2, :, MODULO(mepos - ip, np) + 1) - 1
            ii = 0
            ! this loop takes about the same time as the message passing call
            ! If I read the code correctly then:
            DO iz = lb(3), ub(3)
               DO iy = lb(2), ub(2)
                  DO ix = lb(1), ub(1)
                     ii = ii + 1
                     grid(ix, iy, iz) = sendbuf(ii)
                  END DO
               END DO
            END DO
            IF (ip .NE. np - 1) THEN
               CALL mp_waitall(req)
            END IF
            CALL MOVE_ALLOC(sendbuf, swaparray)
            CALL MOVE_ALLOC(recvbuf, sendbuf)
            CALL MOVE_ALLOC(swaparray, recvbuf)
         END DO
         IF (rs%desc%border > 0) THEN
!$OMP       PARALLEL DO DEFAULT(NONE) &
!$OMP                   PRIVATE(i,im,j,jm,k,km) &
!$OMP                   SHARED(rs)
            DO k = rs%lb_local(3), rs%ub_local(3)
               IF (k < rs%lb_real(3)) THEN
                  km = k + rs%desc%npts(3)
               ELSE IF (k > rs%ub_real(3)) THEN
                  km = k - rs%desc%npts(3)
               ELSE
                  km = k
               END IF
               DO j = rs%lb_local(2), rs%ub_local(2)
                  IF (j < rs%lb_real(2)) THEN
                     jm = j + rs%desc%npts(2)
                  ELSE IF (j > rs%ub_real(2)) THEN
                     jm = j - rs%desc%npts(2)
                  ELSE
                     jm = j
                  END IF
                  DO i = rs%lb_local(1), rs%ub_local(1)
                     IF (i < rs%lb_real(1)) THEN
                        im = i + rs%desc%npts(1)
                     ELSE IF (i > rs%ub_real(1)) THEN
                        im = i - rs%desc%npts(1)
                     ELSE
                        im = i
                     END IF
                     rs%r(i, j, k) = rs%r(im, jm, km)
                  END DO
               END DO
            END DO
!$OMP       END PARALLEL DO
         END IF
      END ASSOCIATE

      DEALLOCATE (rcount)
      DEALLOCATE (sendbuf)
      DEALLOCATE (recvbuf)

   END SUBROUTINE transfer_pw2rs_replicated

! **************************************************************************************************
!> \brief does the rs2pw transfer in the case where the rs grid is
!>       distributed (3D domain decomposition)
!> \param rs ...
!> \param pw ...
!> \par History
!>      12.2007 created [Matt Watkins]
!>      9.2008 reduced amount of halo data sent [Iain Bethune]
!>      10.2008 added non-blocking communication [Iain Bethune]
!>      4.2009 added support for rank-reordering on the grid [Iain Bethune]
!>      12.2009 added OMP and sparse alltoall [Iain Bethune]
!>              (c) The Numerical Algorithms Group (NAG) Ltd, 2008-2009 on behalf of the HECToR project
!> \note
!>       the transfer is a two step procedure. For example, for the rs2pw transfer:
!>
!>       1) Halo-exchange in 3D so that the local part of the rs_grid contains the full data
!>       2) an alltoall communication to redistribute the local rs_grid to the local pw_grid
!>
!>       the halo exchange is most expensive on a large number of CPUs. Particular in this halo
!>       exchange is that the border region is rather large (e.g. 20 points) and that it might overlap
!>       with the central domain of several CPUs (i.e. next nearest neighbors)
! **************************************************************************************************
   SUBROUTINE transfer_rs2pw_distributed(rs, pw)
      TYPE(realspace_grid_type), INTENT(IN)              :: rs
      TYPE(pw_type), INTENT(IN)                          :: pw

      CHARACTER(LEN=200)                                 :: error_string
      INTEGER :: completed, dest_down, dest_up, i, idir, j, k, lb, my_id, my_pw_rank, my_rs_rank, &
         n_shifts, nn, num_threads, position, source_down, source_up, ub, x, y, z
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dshifts, recv_disps, recv_sizes, &
                                                            send_disps, send_sizes, ushifts
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: bounds, recv_tasks, send_tasks
      INTEGER, DIMENSION(2)                              :: neighbours, pos
      INTEGER, DIMENSION(3) :: coords, lb_recv, lb_recv_down, lb_recv_up, lb_send, lb_send_down, &
         lb_send_up, ub_recv, ub_recv_down, ub_recv_up, ub_send, ub_send_down, ub_send_up
      LOGICAL, DIMENSION(3)                              :: halo_swapped
      REAL(KIND=dp)                                      :: pw_sum, rs_sum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: recv_buf_3d_down, recv_buf_3d_up, &
                                                            send_buf_3d_down, send_buf_3d_up
      TYPE(cp_1d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: recv_bufs, send_bufs
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: recv_reqs, send_reqs
      TYPE(mp_request_type), DIMENSION(4)                :: req

      num_threads = 1
      my_id = 0

      ! safety check, to be removed once we're absolute sure the routine is correct
      IF (debug_this_module) THEN
         rs_sum = accurate_sum(rs%r)*ABS(det_3x3(rs%desc%dh))
         CALL rs%desc%group%sum(rs_sum)
      END IF

      halo_swapped = .FALSE.
      ! We don't need to send the 'edges' of the halos that have already been sent
      ! Halos are contiguous in memory in z-direction only, so swap these first,
      ! and send less data in the y and x directions which are more expensive

      DO idir = 3, 1, -1

         IF (rs%desc%perd(idir) .NE. 1) THEN

            ALLOCATE (dshifts(0:rs%desc%neighbours(idir)))
            ALLOCATE (ushifts(0:rs%desc%neighbours(idir)))

            ushifts = 0
            dshifts = 0

            ! check that we don't try to send data to ourself
            DO n_shifts = 1, MIN(rs%desc%neighbours(idir), rs%desc%group_dim(idir) - 1)

               ! need to take into account the possible varying widths of neighbouring cells
               ! offset_up and offset_down hold the real size of the neighbouring cells
               position = MODULO(rs%desc%virtual_group_coor(idir) - n_shifts, rs%desc%group_dim(idir))
               neighbours = get_limit(rs%desc%npts(idir), rs%desc%group_dim(idir), position)
               dshifts(n_shifts) = dshifts(n_shifts - 1) + (neighbours(2) - neighbours(1) + 1)

               position = MODULO(rs%desc%virtual_group_coor(idir) + n_shifts, rs%desc%group_dim(idir))
               neighbours = get_limit(rs%desc%npts(idir), rs%desc%group_dim(idir), position)
               ushifts(n_shifts) = ushifts(n_shifts - 1) + (neighbours(2) - neighbours(1) + 1)

               ! The border data has to be send/received from the neighbours
               ! First we calculate the source and destination processes for the shift
               ! We do both shifts at once to allow for more overlap of communication and buffer packing/unpacking

               CALL cart_shift(rs, idir, -1*n_shifts, source_down, dest_down)

               lb_send_down(:) = rs%lb_local(:)
               lb_recv_down(:) = rs%lb_local(:)
               ub_recv_down(:) = rs%ub_local(:)
               ub_send_down(:) = rs%ub_local(:)

               IF (dshifts(n_shifts - 1) .LE. rs%desc%border) THEN
                  ub_send_down(idir) = lb_send_down(idir) + rs%desc%border - 1 - dshifts(n_shifts - 1)
                  lb_send_down(idir) = MAX(lb_send_down(idir), &
                                           lb_send_down(idir) + rs%desc%border - dshifts(n_shifts))

                  ub_recv_down(idir) = ub_recv_down(idir) - rs%desc%border
                  lb_recv_down(idir) = MAX(lb_recv_down(idir) + rs%desc%border, &
                                           ub_recv_down(idir) - rs%desc%border + 1 + ushifts(n_shifts - 1))
               ELSE
                  lb_send_down(idir) = 0
                  ub_send_down(idir) = -1
                  lb_recv_down(idir) = 0
                  ub_recv_down(idir) = -1
               END IF

               DO i = 1, 3
                  IF (halo_swapped(i)) THEN
                     lb_send_down(i) = rs%lb_real(i)
                     ub_send_down(i) = rs%ub_real(i)
                     lb_recv_down(i) = rs%lb_real(i)
                     ub_recv_down(i) = rs%ub_real(i)
                  END IF
               END DO

               ! post the receive
               ALLOCATE (recv_buf_3d_down(lb_recv_down(1):ub_recv_down(1), &
                                          lb_recv_down(2):ub_recv_down(2), lb_recv_down(3):ub_recv_down(3)))
               CALL rs%desc%group%irecv(recv_buf_3d_down, source_down, req(1))

               ! now allocate, pack and send the send buffer
               nn = PRODUCT(ub_send_down - lb_send_down + 1)
               ALLOCATE (send_buf_3d_down(lb_send_down(1):ub_send_down(1), &
                                          lb_send_down(2):ub_send_down(2), lb_send_down(3):ub_send_down(3)))

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP          PRIVATE(lb,ub,my_id,NUM_THREADS), &
!$OMP          SHARED(send_buf_3d_down,rs,lb_send_down,ub_send_down)
!$             num_threads = MIN(omp_get_max_threads(), ub_send_down(3) - lb_send_down(3) + 1)
!$             my_id = omp_get_thread_num()
               IF (my_id < num_threads) THEN
                  lb = lb_send_down(3) + ((ub_send_down(3) - lb_send_down(3) + 1)*my_id)/num_threads
                  ub = lb_send_down(3) + ((ub_send_down(3) - lb_send_down(3) + 1)*(my_id + 1))/num_threads - 1

                  send_buf_3d_down(lb_send_down(1):ub_send_down(1), lb_send_down(2):ub_send_down(2), &
                                   lb:ub) = rs%r(lb_send_down(1):ub_send_down(1), &
                                                 lb_send_down(2):ub_send_down(2), lb:ub)
               END IF
!$OMP END PARALLEL

               CALL rs%desc%group%isend(send_buf_3d_down, dest_down, req(3))

               ! Now for the other direction
               CALL cart_shift(rs, idir, n_shifts, source_up, dest_up)

               lb_send_up(:) = rs%lb_local(:)
               lb_recv_up(:) = rs%lb_local(:)
               ub_recv_up(:) = rs%ub_local(:)
               ub_send_up(:) = rs%ub_local(:)

               IF (ushifts(n_shifts - 1) .LE. rs%desc%border) THEN

                  lb_send_up(idir) = ub_send_up(idir) - rs%desc%border + 1 + ushifts(n_shifts - 1)
                  ub_send_up(idir) = MIN(ub_send_up(idir), &
                                         ub_send_up(idir) - rs%desc%border + ushifts(n_shifts))

                  lb_recv_up(idir) = lb_recv_up(idir) + rs%desc%border
                  ub_recv_up(idir) = MIN(ub_recv_up(idir) - rs%desc%border, &
                                         lb_recv_up(idir) + rs%desc%border - 1 - dshifts(n_shifts - 1))
               ELSE
                  lb_send_up(idir) = 0
                  ub_send_up(idir) = -1
                  lb_recv_up(idir) = 0
                  ub_recv_up(idir) = -1
               END IF

               DO i = 1, 3
                  IF (halo_swapped(i)) THEN
                     lb_send_up(i) = rs%lb_real(i)
                     ub_send_up(i) = rs%ub_real(i)
                     lb_recv_up(i) = rs%lb_real(i)
                     ub_recv_up(i) = rs%ub_real(i)
                  END IF
               END DO

               ! post the receive
               ALLOCATE (recv_buf_3d_up(lb_recv_up(1):ub_recv_up(1), &
                                        lb_recv_up(2):ub_recv_up(2), lb_recv_up(3):ub_recv_up(3)))
               CALL rs%desc%group%irecv(recv_buf_3d_up, source_up, req(2))

               ! now allocate,pack and send the send buffer
               nn = PRODUCT(ub_send_up - lb_send_up + 1)
               ALLOCATE (send_buf_3d_up(lb_send_up(1):ub_send_up(1), &
                                        lb_send_up(2):ub_send_up(2), lb_send_up(3):ub_send_up(3)))

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP          PRIVATE(lb,ub,my_id,NUM_THREADS), &
!$OMP          SHARED(send_buf_3d_up,rs,lb_send_up,ub_send_up)
!$             num_threads = MIN(omp_get_max_threads(), ub_send_up(3) - lb_send_up(3) + 1)
!$             my_id = omp_get_thread_num()
               IF (my_id < num_threads) THEN
                  lb = lb_send_up(3) + ((ub_send_up(3) - lb_send_up(3) + 1)*my_id)/num_threads
                  ub = lb_send_up(3) + ((ub_send_up(3) - lb_send_up(3) + 1)*(my_id + 1))/num_threads - 1

                  send_buf_3d_up(lb_send_up(1):ub_send_up(1), lb_send_up(2):ub_send_up(2), &
                                 lb:ub) = rs%r(lb_send_up(1):ub_send_up(1), &
                                               lb_send_up(2):ub_send_up(2), lb:ub)
               END IF
!$OMP END PARALLEL

               CALL rs%desc%group%isend(send_buf_3d_up, dest_up, req(4))

               ! wait for a recv to complete, then we can unpack

               DO i = 1, 2

                  CALL mp_waitany(req(1:2), completed)

                  IF (completed .EQ. 1) THEN

                     ! only some procs may need later shifts
                     IF (ub_recv_down(idir) .GE. lb_recv_down(idir)) THEN
                        ! Sum the data in the RS Grid
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP          PRIVATE(lb,ub,my_id,NUM_THREADS), &
!$OMP          SHARED(recv_buf_3d_down,rs,lb_recv_down,ub_recv_down)
!$                      num_threads = MIN(omp_get_max_threads(), ub_recv_down(3) - lb_recv_down(3) + 1)
!$                      my_id = omp_get_thread_num()
                        IF (my_id < num_threads) THEN
                           lb = lb_recv_down(3) + ((ub_recv_down(3) - lb_recv_down(3) + 1)*my_id)/num_threads
                           ub = lb_recv_down(3) + ((ub_recv_down(3) - lb_recv_down(3) + 1)*(my_id + 1))/num_threads - 1

                           rs%r(lb_recv_down(1):ub_recv_down(1), &
                                lb_recv_down(2):ub_recv_down(2), lb:ub) = &
                              rs%r(lb_recv_down(1):ub_recv_down(1), &
                                   lb_recv_down(2):ub_recv_down(2), lb:ub) + &
                              recv_buf_3d_down(:, :, lb:ub)
                        END IF
!$OMP END PARALLEL
                     END IF
                     DEALLOCATE (recv_buf_3d_down)
                  ELSE

                     ! only some procs may need later shifts
                     IF (ub_recv_up(idir) .GE. lb_recv_up(idir)) THEN
                        ! Sum the data in the RS Grid
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP          PRIVATE(lb,ub,my_id,NUM_THREADS), &
!$OMP          SHARED(recv_buf_3d_up,rs,lb_recv_up,ub_recv_up)
!$                      num_threads = MIN(omp_get_max_threads(), ub_recv_up(3) - lb_recv_up(3) + 1)
!$                      my_id = omp_get_thread_num()
                        IF (my_id < num_threads) THEN
                           lb = lb_recv_up(3) + ((ub_recv_up(3) - lb_recv_up(3) + 1)*my_id)/num_threads
                           ub = lb_recv_up(3) + ((ub_recv_up(3) - lb_recv_up(3) + 1)*(my_id + 1))/num_threads - 1

                           rs%r(lb_recv_up(1):ub_recv_up(1), &
                                lb_recv_up(2):ub_recv_up(2), lb:ub) = &
                              rs%r(lb_recv_up(1):ub_recv_up(1), &
                                   lb_recv_up(2):ub_recv_up(2), lb:ub) + &
                              recv_buf_3d_up(:, :, lb:ub)
                        END IF
!$OMP END PARALLEL
                     END IF
                     DEALLOCATE (recv_buf_3d_up)
                  END IF

               END DO

               ! make sure the sends have completed before we deallocate

               CALL mp_waitall(req(3:4))

               DEALLOCATE (send_buf_3d_down)
               DEALLOCATE (send_buf_3d_up)
            END DO

            DEALLOCATE (dshifts)
            DEALLOCATE (ushifts)

         END IF

         halo_swapped(idir) = .TRUE.

      END DO

      ! This is the real redistribution
      ALLOCATE (bounds(0:pw%pw_grid%para%group_size - 1, 1:4))

      ! work out the pw grid points each proc holds
      DO i = 0, pw%pw_grid%para%group_size - 1
         bounds(i, 1:2) = pw%pw_grid%para%bo(1:2, 1, i, 1)
         bounds(i, 3:4) = pw%pw_grid%para%bo(1:2, 2, i, 1)
         bounds(i, 1:2) = bounds(i, 1:2) - pw%pw_grid%npts(1)/2 - 1
         bounds(i, 3:4) = bounds(i, 3:4) - pw%pw_grid%npts(2)/2 - 1
      END DO

      ALLOCATE (send_tasks(0:pw%pw_grid%para%group_size - 1, 1:6))
      ALLOCATE (send_sizes(0:pw%pw_grid%para%group_size - 1))
      ALLOCATE (send_disps(0:pw%pw_grid%para%group_size - 1))
      ALLOCATE (recv_tasks(0:pw%pw_grid%para%group_size - 1, 1:6))
      ALLOCATE (recv_sizes(0:pw%pw_grid%para%group_size - 1))
      ALLOCATE (recv_disps(0:pw%pw_grid%para%group_size - 1))
      send_tasks(:, 1) = 1
      send_tasks(:, 2) = 0
      send_tasks(:, 3) = 1
      send_tasks(:, 4) = 0
      send_tasks(:, 5) = 1
      send_tasks(:, 6) = 0
      send_sizes = 0
      recv_sizes = 0

      my_rs_rank = rs%desc%my_pos
      my_pw_rank = pw%pw_grid%para%rs_mpo

      ! find the processors that should hold our data
      ! should be part of the rs grid type
      ! this is a loop over real ranks (i.e. the in-order cartesian ranks)
      ! do the recv and send tasks in two separate loops which will
      ! load balance better for OpenMP with large numbers of MPI tasks

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(coords,idir,pos,lb_send,ub_send), &
!$OMP             SHARED(rs,bounds,my_rs_rank,recv_tasks,recv_sizes)
      DO i = 0, rs%desc%group_size - 1

         coords(:) = rs%desc%rank2coord(:, rs%desc%real2virtual(i))
         !calculate the rs grid points on each processor
         !coords is the part of the grid that rank i actually holds
         DO idir = 1, 3
            pos(:) = get_limit(rs%desc%npts(idir), rs%desc%group_dim(idir), coords(idir))
            pos(:) = pos(:) - rs%desc%npts(idir)/2 - 1
            lb_send(idir) = pos(1)
            ub_send(idir) = pos(2)
         END DO

         IF (lb_send(1) .GT. bounds(my_rs_rank, 2)) CYCLE
         IF (ub_send(1) .LT. bounds(my_rs_rank, 1)) CYCLE
         IF (lb_send(2) .GT. bounds(my_rs_rank, 4)) CYCLE
         IF (ub_send(2) .LT. bounds(my_rs_rank, 3)) CYCLE

         recv_tasks(i, 1) = MAX(lb_send(1), bounds(my_rs_rank, 1))
         recv_tasks(i, 2) = MIN(ub_send(1), bounds(my_rs_rank, 2))
         recv_tasks(i, 3) = MAX(lb_send(2), bounds(my_rs_rank, 3))
         recv_tasks(i, 4) = MIN(ub_send(2), bounds(my_rs_rank, 4))
         recv_tasks(i, 5) = lb_send(3)
         recv_tasks(i, 6) = ub_send(3)
         recv_sizes(i) = (recv_tasks(i, 2) - recv_tasks(i, 1) + 1)* &
                         (recv_tasks(i, 4) - recv_tasks(i, 3) + 1)*(recv_tasks(i, 6) - recv_tasks(i, 5) + 1)

      END DO
!$OMP END PARALLEL DO

      coords(:) = rs%desc%rank2coord(:, rs%desc%real2virtual(my_rs_rank))
      DO idir = 1, 3
         pos(:) = get_limit(rs%desc%npts(idir), rs%desc%group_dim(idir), coords(idir))
         pos(:) = pos(:) - rs%desc%npts(idir)/2 - 1
         lb_send(idir) = pos(1)
         ub_send(idir) = pos(2)
      END DO

      lb_recv(:) = lb_send(:)
      ub_recv(:) = ub_send(:)
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             SHARED(pw,lb_send,ub_send,bounds,send_tasks,send_sizes)
      DO j = 0, pw%pw_grid%para%group_size - 1

         IF (lb_send(1) .GT. bounds(j, 2)) CYCLE
         IF (ub_send(1) .LT. bounds(j, 1)) CYCLE
         IF (lb_send(2) .GT. bounds(j, 4)) CYCLE
         IF (ub_send(2) .LT. bounds(j, 3)) CYCLE

         send_tasks(j, 1) = MAX(lb_send(1), bounds(j, 1))
         send_tasks(j, 2) = MIN(ub_send(1), bounds(j, 2))
         send_tasks(j, 3) = MAX(lb_send(2), bounds(j, 3))
         send_tasks(j, 4) = MIN(ub_send(2), bounds(j, 4))
         send_tasks(j, 5) = lb_send(3)
         send_tasks(j, 6) = ub_send(3)
         send_sizes(j) = (send_tasks(j, 2) - send_tasks(j, 1) + 1)* &
                         (send_tasks(j, 4) - send_tasks(j, 3) + 1)*(send_tasks(j, 6) - send_tasks(j, 5) + 1)

      END DO
!$OMP END PARALLEL DO

      send_disps(0) = 0
      recv_disps(0) = 0
      DO i = 1, pw%pw_grid%para%group_size - 1
         send_disps(i) = send_disps(i - 1) + send_sizes(i - 1)
         recv_disps(i) = recv_disps(i - 1) + recv_sizes(i - 1)
      END DO

      CPASSERT(SUM(send_sizes) == PRODUCT(ub_recv - lb_recv + 1))

      ALLOCATE (send_bufs(0:rs%desc%group_size - 1))
      ALLOCATE (recv_bufs(0:rs%desc%group_size - 1))

      DO i = 0, rs%desc%group_size - 1
         IF (send_sizes(i) .NE. 0) THEN
            ALLOCATE (send_bufs(i)%array(send_sizes(i)))
         ELSE
            NULLIFY (send_bufs(i)%array)
         END IF
         IF (recv_sizes(i) .NE. 0) THEN
            ALLOCATE (recv_bufs(i)%array(recv_sizes(i)))
         ELSE
            NULLIFY (recv_bufs(i)%array)
         END IF
      END DO

      ALLOCATE (recv_reqs(0:rs%desc%group_size - 1))
      recv_reqs = mp_request_null

      DO i = 0, rs%desc%group_size - 1
         IF (recv_sizes(i) .NE. 0) THEN
            CALL rs%desc%group%irecv(recv_bufs(i)%array, i, recv_reqs(i))
         END IF
      END DO

      ! do packing
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(k,z,y,x), &
!$OMP             SHARED(rs,send_tasks,send_bufs,send_disps)
      DO i = 0, rs%desc%group_size - 1
         k = 0
         DO z = send_tasks(i, 5), send_tasks(i, 6)
            DO y = send_tasks(i, 3), send_tasks(i, 4)
               DO x = send_tasks(i, 1), send_tasks(i, 2)
                  k = k + 1
                  send_bufs(i)%array(k) = rs%r(x, y, z)
               END DO
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      ALLOCATE (send_reqs(0:rs%desc%group_size - 1))
      send_reqs = mp_request_null

      DO i = 0, rs%desc%group_size - 1
         IF (send_sizes(i) .NE. 0) THEN
            CALL rs%desc%group%isend(send_bufs(i)%array, i, send_reqs(i))
         END IF
      END DO

      ! do unpacking
      ! no OMP here so we can unpack each message as it arrives
      DO i = 0, rs%desc%group_size - 1
         IF (recv_sizes(i) .EQ. 0) CYCLE

         CALL mp_waitany(recv_reqs, completed)
         k = 0
         DO z = recv_tasks(completed - 1, 5), recv_tasks(completed - 1, 6)
            DO y = recv_tasks(completed - 1, 3), recv_tasks(completed - 1, 4)
               DO x = recv_tasks(completed - 1, 1), recv_tasks(completed - 1, 2)
                  k = k + 1
                  pw%cr3d(x, y, z) = recv_bufs(completed - 1)%array(k)
               END DO
            END DO
         END DO
      END DO

      CALL mp_waitall(send_reqs)

      DEALLOCATE (recv_reqs)
      DEALLOCATE (send_reqs)

      DO i = 0, rs%desc%group_size - 1
         IF (ASSOCIATED(send_bufs(i)%array)) THEN
            DEALLOCATE (send_bufs(i)%array)
         END IF
         IF (ASSOCIATED(recv_bufs(i)%array)) THEN
            DEALLOCATE (recv_bufs(i)%array)
         END IF
      END DO

      DEALLOCATE (send_bufs)
      DEALLOCATE (recv_bufs)
      DEALLOCATE (send_tasks)
      DEALLOCATE (send_sizes)
      DEALLOCATE (send_disps)
      DEALLOCATE (recv_tasks)
      DEALLOCATE (recv_sizes)
      DEALLOCATE (recv_disps)

      IF (debug_this_module) THEN
         ! safety check, to be removed once we're absolute sure the routine is correct
         pw_sum = pw_integrate_function(pw)
         IF (ABS(pw_sum - rs_sum)/MAX(1.0_dp, ABS(pw_sum), ABS(rs_sum)) > EPSILON(rs_sum)*1000) THEN
            WRITE (error_string, '(A,6(1X,I4.4),3F25.16)') "rs_pw_transfer_distributed", &
               rs%desc%npts, rs%desc%group_dim, pw_sum, rs_sum, ABS(pw_sum - rs_sum)
            CALL cp_abort(__LOCATION__, &
                          error_string//" Please report this bug ... quick workaround: use "// &
                          "DISTRIBUTION_TYPE REPLICATED")
         END IF
      END IF

   END SUBROUTINE transfer_rs2pw_distributed

! **************************************************************************************************
!> \brief does the pw2rs transfer in the case where the rs grid is
!>       distributed (3D domain decomposition)
!> \param rs ...
!> \param pw ...
!> \par History
!>      12.2007 created [Matt Watkins]
!>      9.2008 reduced amount of halo data sent [Iain Bethune]
!>      10.2008 added non-blocking communication [Iain Bethune]
!>      4.2009 added support for rank-reordering on the grid [Iain Bethune]
!>      12.2009 added OMP and sparse alltoall [Iain Bethune]
!>              (c) The Numerical Algorithms Group (NAG) Ltd, 2008-2009 on behalf of the HECToR project
!> \note
!>       the transfer is a two step procedure. For example, for the rs2pw transfer:
!>
!>       1) Halo-exchange in 3D so that the local part of the rs_grid contains the full data
!>       2) an alltoall communication to redistribute the local rs_grid to the local pw_grid
!>
!>       the halo exchange is most expensive on a large number of CPUs. Particular in this halo
!>       exchange is that the border region is rather large (e.g. 20 points) and that it might overlap
!>       with the central domain of several CPUs (i.e. next nearest neighbors)
! **************************************************************************************************
   SUBROUTINE transfer_pw2rs_distributed(rs, pw)
      TYPE(realspace_grid_type), INTENT(IN)              :: rs
      TYPE(pw_type), INTENT(IN)                          :: pw

      INTEGER :: completed, dest_down, dest_up, i, idir, j, k, lb, my_id, my_pw_rank, my_rs_rank, &
         n_shifts, nn, num_threads, position, source_down, source_up, ub, x, y, z
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dshifts, recv_disps, recv_sizes, &
                                                            send_disps, send_sizes, ushifts
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: bounds, recv_tasks, send_tasks
      INTEGER, DIMENSION(2)                              :: neighbours, pos
      INTEGER, DIMENSION(3) :: coords, lb_recv, lb_recv_down, lb_recv_up, lb_send, lb_send_down, &
         lb_send_up, ub_recv, ub_recv_down, ub_recv_up, ub_send, ub_send_down, ub_send_up
      LOGICAL, DIMENSION(3)                              :: halo_swapped
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: recv_buf_3d_down, recv_buf_3d_up, &
                                                            send_buf_3d_down, send_buf_3d_up
      TYPE(cp_1d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: recv_bufs, send_bufs
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: recv_reqs, send_reqs
      TYPE(mp_request_type), DIMENSION(4)                :: req

      num_threads = 1
      my_id = 0

      CALL rs_grid_zero(rs)

      ! This is the real redistribution

      ALLOCATE (bounds(0:pw%pw_grid%para%group_size - 1, 1:4))

      DO i = 0, pw%pw_grid%para%group_size - 1
         bounds(i, 1:2) = pw%pw_grid%para%bo(1:2, 1, i, 1)
         bounds(i, 3:4) = pw%pw_grid%para%bo(1:2, 2, i, 1)
         bounds(i, 1:2) = bounds(i, 1:2) - pw%pw_grid%npts(1)/2 - 1
         bounds(i, 3:4) = bounds(i, 3:4) - pw%pw_grid%npts(2)/2 - 1
      END DO

      ALLOCATE (send_tasks(0:pw%pw_grid%para%group_size - 1, 1:6))
      ALLOCATE (send_sizes(0:pw%pw_grid%para%group_size - 1))
      ALLOCATE (send_disps(0:pw%pw_grid%para%group_size - 1))
      ALLOCATE (recv_tasks(0:pw%pw_grid%para%group_size - 1, 1:6))
      ALLOCATE (recv_sizes(0:pw%pw_grid%para%group_size - 1))
      ALLOCATE (recv_disps(0:pw%pw_grid%para%group_size - 1))

      send_tasks = 0
      send_tasks(:, 1) = 1
      send_tasks(:, 2) = 0
      send_tasks(:, 3) = 1
      send_tasks(:, 4) = 0
      send_tasks(:, 5) = 1
      send_tasks(:, 6) = 0
      send_sizes = 0

      recv_tasks = 0
      recv_tasks(:, 1) = 1
      recv_tasks(:, 2) = 0
      send_tasks(:, 3) = 1
      send_tasks(:, 4) = 0
      send_tasks(:, 5) = 1
      send_tasks(:, 6) = 0
      recv_sizes = 0

      my_rs_rank = rs%desc%my_pos
      my_pw_rank = pw%pw_grid%para%rs_mpo

      ! find the processors that should hold our data
      ! should be part of the rs grid type
      ! this is a loop over real ranks (i.e. the in-order cartesian ranks)
      ! do the recv and send tasks in two separate loops which will
      ! load balance better for OpenMP with large numbers of MPI tasks

      ! this is the reverse of rs2pw: what were the sends are now the recvs

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(coords,idir,pos,lb_send,ub_send), &
!$OMP             SHARED(rs,bounds,my_rs_rank,send_tasks,send_sizes,pw)
      DO i = 0, pw%pw_grid%para%group_size - 1

         coords(:) = rs%desc%rank2coord(:, rs%desc%real2virtual(i))
         !calculate the real rs grid points on each processor
         !coords is the part of the grid that rank i actually holds
         DO idir = 1, 3
            pos(:) = get_limit(rs%desc%npts(idir), rs%desc%group_dim(idir), coords(idir))
            pos(:) = pos(:) - rs%desc%npts(idir)/2 - 1
            lb_send(idir) = pos(1)
            ub_send(idir) = pos(2)
         END DO

         IF (ub_send(1) .LT. bounds(my_rs_rank, 1)) CYCLE
         IF (lb_send(1) .GT. bounds(my_rs_rank, 2)) CYCLE
         IF (ub_send(2) .LT. bounds(my_rs_rank, 3)) CYCLE
         IF (lb_send(2) .GT. bounds(my_rs_rank, 4)) CYCLE

         send_tasks(i, 1) = MAX(lb_send(1), bounds(my_rs_rank, 1))
         send_tasks(i, 2) = MIN(ub_send(1), bounds(my_rs_rank, 2))
         send_tasks(i, 3) = MAX(lb_send(2), bounds(my_rs_rank, 3))
         send_tasks(i, 4) = MIN(ub_send(2), bounds(my_rs_rank, 4))
         send_tasks(i, 5) = lb_send(3)
         send_tasks(i, 6) = ub_send(3)
         send_sizes(i) = (send_tasks(i, 2) - send_tasks(i, 1) + 1)* &
                         (send_tasks(i, 4) - send_tasks(i, 3) + 1)*(send_tasks(i, 6) - send_tasks(i, 5) + 1)

      END DO
!$OMP END PARALLEL DO

      coords(:) = rs%desc%rank2coord(:, rs%desc%real2virtual(my_rs_rank))
      DO idir = 1, 3
         pos(:) = get_limit(rs%desc%npts(idir), rs%desc%group_dim(idir), coords(idir))
         pos(:) = pos(:) - rs%desc%npts(idir)/2 - 1
         lb_send(idir) = pos(1)
         ub_send(idir) = pos(2)
      END DO

      lb_recv(:) = lb_send(:)
      ub_recv(:) = ub_send(:)

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             SHARED(pw,lb_send,ub_send,bounds,recv_tasks,recv_sizes)
      DO j = 0, pw%pw_grid%para%group_size - 1

         IF (ub_send(1) .LT. bounds(j, 1)) CYCLE
         IF (lb_send(1) .GT. bounds(j, 2)) CYCLE
         IF (ub_send(2) .LT. bounds(j, 3)) CYCLE
         IF (lb_send(2) .GT. bounds(j, 4)) CYCLE

         recv_tasks(j, 1) = MAX(lb_send(1), bounds(j, 1))
         recv_tasks(j, 2) = MIN(ub_send(1), bounds(j, 2))
         recv_tasks(j, 3) = MAX(lb_send(2), bounds(j, 3))
         recv_tasks(j, 4) = MIN(ub_send(2), bounds(j, 4))
         recv_tasks(j, 5) = lb_send(3)
         recv_tasks(j, 6) = ub_send(3)
         recv_sizes(j) = (recv_tasks(j, 2) - recv_tasks(j, 1) + 1)* &
                         (recv_tasks(j, 4) - recv_tasks(j, 3) + 1)*(recv_tasks(j, 6) - recv_tasks(j, 5) + 1)

      END DO
!$OMP END PARALLEL DO

      send_disps(0) = 0
      recv_disps(0) = 0
      DO i = 1, pw%pw_grid%para%group_size - 1
         send_disps(i) = send_disps(i - 1) + send_sizes(i - 1)
         recv_disps(i) = recv_disps(i - 1) + recv_sizes(i - 1)
      END DO

      CPASSERT(SUM(recv_sizes) == PRODUCT(ub_recv - lb_recv + 1))

      ALLOCATE (send_bufs(0:rs%desc%group_size - 1))
      ALLOCATE (recv_bufs(0:rs%desc%group_size - 1))

      DO i = 0, rs%desc%group_size - 1
         IF (send_sizes(i) .NE. 0) THEN
            ALLOCATE (send_bufs(i)%array(send_sizes(i)))
         ELSE
            NULLIFY (send_bufs(i)%array)
         END IF
         IF (recv_sizes(i) .NE. 0) THEN
            ALLOCATE (recv_bufs(i)%array(recv_sizes(i)))
         ELSE
            NULLIFY (recv_bufs(i)%array)
         END IF
      END DO

      ALLOCATE (recv_reqs(0:rs%desc%group_size - 1))
      recv_reqs = mp_request_null

      DO i = 0, rs%desc%group_size - 1
         IF (recv_sizes(i) .NE. 0) THEN
            CALL rs%desc%group%irecv(recv_bufs(i)%array, i, recv_reqs(i))
         END IF
      END DO

      ! do packing
!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(k,z,y,x), &
!$OMP             SHARED(pw,rs,send_tasks,send_bufs,send_disps)
      DO i = 0, rs%desc%group_size - 1
         k = 0
         DO z = send_tasks(i, 5), send_tasks(i, 6)
            DO y = send_tasks(i, 3), send_tasks(i, 4)
               DO x = send_tasks(i, 1), send_tasks(i, 2)
                  k = k + 1
                  send_bufs(i)%array(k) = pw%cr3d(x, y, z)
               END DO
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      ALLOCATE (send_reqs(0:rs%desc%group_size - 1))
      send_reqs = mp_request_null

      DO i = 0, rs%desc%group_size - 1
         IF (send_sizes(i) .NE. 0) THEN
            CALL rs%desc%group%isend(send_bufs(i)%array, i, send_reqs(i))
         END IF
      END DO

      ! do unpacking
      ! no OMP here so we can unpack each message as it arrives

      DO i = 0, rs%desc%group_size - 1
         IF (recv_sizes(i) .EQ. 0) CYCLE

         CALL mp_waitany(recv_reqs, completed)
         k = 0
         DO z = recv_tasks(completed - 1, 5), recv_tasks(completed - 1, 6)
            DO y = recv_tasks(completed - 1, 3), recv_tasks(completed - 1, 4)
               DO x = recv_tasks(completed - 1, 1), recv_tasks(completed - 1, 2)
                  k = k + 1
                  rs%r(x, y, z) = recv_bufs(completed - 1)%array(k)
               END DO
            END DO
         END DO
      END DO

      CALL mp_waitall(send_reqs)

      DEALLOCATE (recv_reqs)
      DEALLOCATE (send_reqs)

      DO i = 0, rs%desc%group_size - 1
         IF (ASSOCIATED(send_bufs(i)%array)) THEN
            DEALLOCATE (send_bufs(i)%array)
         END IF
         IF (ASSOCIATED(recv_bufs(i)%array)) THEN
            DEALLOCATE (recv_bufs(i)%array)
         END IF
      END DO

      DEALLOCATE (send_bufs)
      DEALLOCATE (recv_bufs)
      DEALLOCATE (send_tasks)
      DEALLOCATE (send_sizes)
      DEALLOCATE (send_disps)
      DEALLOCATE (recv_tasks)
      DEALLOCATE (recv_sizes)
      DEALLOCATE (recv_disps)

      ! now pass wings around
      halo_swapped = .FALSE.

      DO idir = 1, 3

         IF (rs%desc%perd(idir) /= 1) THEN

            ALLOCATE (dshifts(0:rs%desc%neighbours(idir)))
            ALLOCATE (ushifts(0:rs%desc%neighbours(idir)))
            ushifts = 0
            dshifts = 0

            DO n_shifts = 1, rs%desc%neighbours(idir)

               ! need to take into account the possible varying widths of neighbouring cells
               ! ushifts and dshifts hold the real size of the neighbouring cells

               position = MODULO(rs%desc%virtual_group_coor(idir) - n_shifts, rs%desc%group_dim(idir))
               neighbours = get_limit(rs%desc%npts(idir), rs%desc%group_dim(idir), position)
               dshifts(n_shifts) = dshifts(n_shifts - 1) + (neighbours(2) - neighbours(1) + 1)

               position = MODULO(rs%desc%virtual_group_coor(idir) + n_shifts, rs%desc%group_dim(idir))
               neighbours = get_limit(rs%desc%npts(idir), rs%desc%group_dim(idir), position)
               ushifts(n_shifts) = ushifts(n_shifts - 1) + (neighbours(2) - neighbours(1) + 1)

               ! The border data has to be send/received from the neighbors
               ! First we calculate the source and destination processes for the shift
               ! The first shift is "downwards"

               CALL cart_shift(rs, idir, -1*n_shifts, source_down, dest_down)

               lb_send_down(:) = rs%lb_local(:)
               ub_send_down(:) = rs%ub_local(:)
               lb_recv_down(:) = rs%lb_local(:)
               ub_recv_down(:) = rs%ub_local(:)

               IF (dshifts(n_shifts - 1) .LE. rs%desc%border) THEN
                  lb_send_down(idir) = lb_send_down(idir) + rs%desc%border
                  ub_send_down(idir) = MIN(ub_send_down(idir) - rs%desc%border, &
                                           lb_send_down(idir) + rs%desc%border - 1 - dshifts(n_shifts - 1))

                  lb_recv_down(idir) = ub_recv_down(idir) - rs%desc%border + 1 + ushifts(n_shifts - 1)
                  ub_recv_down(idir) = MIN(ub_recv_down(idir), &
                                           ub_recv_down(idir) - rs%desc%border + ushifts(n_shifts))
               ELSE
                  lb_send_down(idir) = 0
                  ub_send_down(idir) = -1
                  lb_recv_down(idir) = 0
                  ub_recv_down(idir) = -1
               END IF

               DO i = 1, 3
                  IF (.NOT. (halo_swapped(i) .OR. i .EQ. idir)) THEN
                     lb_send_down(i) = rs%lb_real(i)
                     ub_send_down(i) = rs%ub_real(i)
                     lb_recv_down(i) = rs%lb_real(i)
                     ub_recv_down(i) = rs%ub_real(i)
                  END IF
               END DO

               ! allocate the recv buffer
               nn = PRODUCT(ub_recv_down - lb_recv_down + 1)
               ALLOCATE (recv_buf_3d_down(lb_recv_down(1):ub_recv_down(1), &
                                          lb_recv_down(2):ub_recv_down(2), lb_recv_down(3):ub_recv_down(3)))

               ! recv buffer is now ready, so post the receive
               CALL rs%desc%group%irecv(recv_buf_3d_down, source_down, req(1))

               ! now allocate,pack and send the send buffer
               nn = PRODUCT(ub_send_down - lb_send_down + 1)
               ALLOCATE (send_buf_3d_down(lb_send_down(1):ub_send_down(1), &
                                          lb_send_down(2):ub_send_down(2), lb_send_down(3):ub_send_down(3)))

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP          PRIVATE(lb,ub,my_id,NUM_THREADS), &
!$OMP          SHARED(send_buf_3d_down,rs,lb_send_down,ub_send_down)
!$             num_threads = MIN(omp_get_max_threads(), ub_send_down(3) - lb_send_down(3) + 1)
!$             my_id = omp_get_thread_num()
               IF (my_id < num_threads) THEN
                  lb = lb_send_down(3) + ((ub_send_down(3) - lb_send_down(3) + 1)*my_id)/num_threads
                  ub = lb_send_down(3) + ((ub_send_down(3) - lb_send_down(3) + 1)*(my_id + 1))/num_threads - 1

                  send_buf_3d_down(lb_send_down(1):ub_send_down(1), lb_send_down(2):ub_send_down(2), &
                                   lb:ub) = rs%r(lb_send_down(1):ub_send_down(1), &
                                                 lb_send_down(2):ub_send_down(2), lb:ub)
               END IF
!$OMP END PARALLEL

               CALL rs%desc%group%isend(send_buf_3d_down, dest_down, req(3))

               ! Now for the other direction

               CALL cart_shift(rs, idir, n_shifts, source_up, dest_up)

               lb_send_up(:) = rs%lb_local(:)
               ub_send_up(:) = rs%ub_local(:)
               lb_recv_up(:) = rs%lb_local(:)
               ub_recv_up(:) = rs%ub_local(:)

               IF (ushifts(n_shifts - 1) .LE. rs%desc%border) THEN
                  ub_send_up(idir) = ub_send_up(idir) - rs%desc%border
                  lb_send_up(idir) = MAX(lb_send_up(idir) + rs%desc%border, &
                                         ub_send_up(idir) - rs%desc%border + 1 + ushifts(n_shifts - 1))

                  ub_recv_up(idir) = lb_recv_up(idir) + rs%desc%border - 1 - dshifts(n_shifts - 1)
                  lb_recv_up(idir) = MAX(lb_recv_up(idir), &
                                         lb_recv_up(idir) + rs%desc%border - dshifts(n_shifts))
               ELSE
                  lb_send_up(idir) = 0
                  ub_send_up(idir) = -1
                  lb_recv_up(idir) = 0
                  ub_recv_up(idir) = -1
               END IF

               DO i = 1, 3
                  IF (.NOT. (halo_swapped(i) .OR. i .EQ. idir)) THEN
                     lb_send_up(i) = rs%lb_real(i)
                     ub_send_up(i) = rs%ub_real(i)
                     lb_recv_up(i) = rs%lb_real(i)
                     ub_recv_up(i) = rs%ub_real(i)
                  END IF
               END DO

               ! allocate the recv buffer
               nn = PRODUCT(ub_recv_up - lb_recv_up + 1)
               ALLOCATE (recv_buf_3d_up(lb_recv_up(1):ub_recv_up(1), &
                                        lb_recv_up(2):ub_recv_up(2), lb_recv_up(3):ub_recv_up(3)))

               ! recv buffer is now ready, so post the receive

               CALL rs%desc%group%irecv(recv_buf_3d_up, source_up, req(2))

               ! now allocate,pack and send the send buffer
               nn = PRODUCT(ub_send_up - lb_send_up + 1)
               ALLOCATE (send_buf_3d_up(lb_send_up(1):ub_send_up(1), &
                                        lb_send_up(2):ub_send_up(2), lb_send_up(3):ub_send_up(3)))

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP          PRIVATE(lb,ub,my_id,NUM_THREADS), &
!$OMP          SHARED(send_buf_3d_up,rs,lb_send_up,ub_send_up)
!$             num_threads = MIN(omp_get_max_threads(), ub_send_up(3) - lb_send_up(3) + 1)
!$             my_id = omp_get_thread_num()
               IF (my_id < num_threads) THEN
                  lb = lb_send_up(3) + ((ub_send_up(3) - lb_send_up(3) + 1)*my_id)/num_threads
                  ub = lb_send_up(3) + ((ub_send_up(3) - lb_send_up(3) + 1)*(my_id + 1))/num_threads - 1

                  send_buf_3d_up(lb_send_up(1):ub_send_up(1), lb_send_up(2):ub_send_up(2), &
                                 lb:ub) = rs%r(lb_send_up(1):ub_send_up(1), &
                                               lb_send_up(2):ub_send_up(2), lb:ub)
               END IF
!$OMP END PARALLEL

               CALL rs%desc%group%isend(send_buf_3d_up, dest_up, req(4))

               ! wait for a recv to complete, then we can unpack

               DO i = 1, 2

                  CALL mp_waitany(req(1:2), completed)

                  IF (completed .EQ. 1) THEN

                     ! only some procs may need later shifts
                     IF (ub_recv_down(idir) .GE. lb_recv_down(idir)) THEN

                        ! Add the data to the RS Grid
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP          PRIVATE(lb,ub,my_id,NUM_THREADS), &
!$OMP          SHARED(recv_buf_3d_down,rs,lb_recv_down,ub_recv_down)
!$                      num_threads = MIN(omp_get_max_threads(), ub_recv_down(3) - lb_recv_down(3) + 1)
!$                      my_id = omp_get_thread_num()
                        IF (my_id < num_threads) THEN
                           lb = lb_recv_down(3) + ((ub_recv_down(3) - lb_recv_down(3) + 1)*my_id)/num_threads
                           ub = lb_recv_down(3) + ((ub_recv_down(3) - lb_recv_down(3) + 1)*(my_id + 1))/num_threads - 1

                           rs%r(lb_recv_down(1):ub_recv_down(1), lb_recv_down(2):ub_recv_down(2), &
                                lb:ub) = recv_buf_3d_down(:, :, lb:ub)
                        END IF
!$OMP END PARALLEL
                     END IF

                     DEALLOCATE (recv_buf_3d_down)
                  ELSE

                     ! only some procs may need later shifts
                     IF (ub_recv_up(idir) .GE. lb_recv_up(idir)) THEN

                        ! Add the data to the RS Grid
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP          PRIVATE(lb,ub,my_id,NUM_THREADS), &
!$OMP          SHARED(recv_buf_3d_up,rs,lb_recv_up,ub_recv_up)
!$                      num_threads = MIN(omp_get_max_threads(), ub_recv_up(3) - lb_recv_up(3) + 1)
!$                      my_id = omp_get_thread_num()
                        IF (my_id < num_threads) THEN
                           lb = lb_recv_up(3) + ((ub_recv_up(3) - lb_recv_up(3) + 1)*my_id)/num_threads
                           ub = lb_recv_up(3) + ((ub_recv_up(3) - lb_recv_up(3) + 1)*(my_id + 1))/num_threads - 1

                           rs%r(lb_recv_up(1):ub_recv_up(1), lb_recv_up(2):ub_recv_up(2), &
                                lb:ub) = recv_buf_3d_up(:, :, lb:ub)
                        END IF
!$OMP END PARALLEL
                     END IF

                     DEALLOCATE (recv_buf_3d_up)
                  END IF
               END DO

               CALL mp_waitall(req(3:4))

               DEALLOCATE (send_buf_3d_down)
               DEALLOCATE (send_buf_3d_up)
            END DO

            DEALLOCATE (ushifts)
            DEALLOCATE (dshifts)
         END IF

         halo_swapped(idir) = .TRUE.

      END DO

   END SUBROUTINE transfer_pw2rs_distributed

! **************************************************************************************************
!> \brief Initialize grid to zero
!> \param rs ...
!> \par History
!>      none
!> \author JGH (23-Mar-2002)
! **************************************************************************************************
   SUBROUTINE rs_grid_zero(rs)

      TYPE(realspace_grid_type), INTENT(IN)              :: rs

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rs_grid_zero'

      INTEGER                                            :: handle, i, j, k, l(3), u(3)

      CALL timeset(routineN, handle)
      l(1) = LBOUND(rs%r, 1); l(2) = LBOUND(rs%r, 2); l(3) = LBOUND(rs%r, 3)
      u(1) = UBOUND(rs%r, 1); u(2) = UBOUND(rs%r, 2); u(3) = UBOUND(rs%r, 3)
!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(3) &
!$OMP             PRIVATE(i,j,k) &
!$OMP             SHARED(rs,l,u)
      DO k = l(3), u(3)
      DO j = l(2), u(2)
      DO i = l(1), u(1)
         rs%r(i, j, k) = 0.0_dp
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
      CALL timestop(handle)

   END SUBROUTINE rs_grid_zero

! **************************************************************************************************
!> \brief rs1(i) = rs1(i) + rs2(i)*rs3(i)
!> \param rs1 ...
!> \param rs2 ...
!> \param rs3 ...
!> \param scalar ...
!> \par History
!>      none
!> \author
! **************************************************************************************************
   SUBROUTINE rs_grid_mult_and_add(rs1, rs2, rs3, scalar)

      TYPE(realspace_grid_type), INTENT(IN)              :: rs1, rs2, rs3
      REAL(dp), INTENT(IN)                               :: scalar

      CHARACTER(len=*), PARAMETER :: routineN = 'rs_grid_mult_and_add'

      INTEGER                                            :: handle, i, j, k, l(3), u(3)

!-----------------------------------------------------------------------------!

      CALL timeset(routineN, handle)
      IF (scalar .NE. 0.0_dp) THEN
         l(1) = LBOUND(rs1%r, 1); l(2) = LBOUND(rs1%r, 2); l(3) = LBOUND(rs1%r, 3)
         u(1) = UBOUND(rs1%r, 1); u(2) = UBOUND(rs1%r, 2); u(3) = UBOUND(rs1%r, 3)
!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(3) &
!$OMP             PRIVATE(i,j,k) &
!$OMP             SHARED(rs1,rs2,rs3,scalar,l,u)
         DO k = l(3), u(3)
         DO j = l(2), u(2)
         DO i = l(1), u(1)
            rs1%r(i, j, k) = rs1%r(i, j, k) + scalar*rs2%r(i, j, k)*rs3%r(i, j, k)
         END DO
         END DO
         END DO
!$OMP END PARALLEL DO
      END IF
      CALL timestop(handle)
   END SUBROUTINE rs_grid_mult_and_add

! **************************************************************************************************
!> \brief Set box matrix info for real space grid
!>      This is needed for variable cell simulations
!> \param pw_grid ...
!> \param rs ...
!> \par History
!>      none
!> \author JGH (15-May-2007)
! **************************************************************************************************
   SUBROUTINE rs_grid_set_box(pw_grid, rs)

      TYPE(pw_grid_type), INTENT(IN), TARGET             :: pw_grid
      TYPE(realspace_grid_type), INTENT(IN)              :: rs

      CPASSERT(ASSOCIATED(rs%desc%pw, pw_grid))
      rs%desc%dh = pw_grid%dh
      rs%desc%dh_inv = pw_grid%dh_inv

   END SUBROUTINE rs_grid_set_box

! **************************************************************************************************
!> \brief retains the given rs grid descriptor (see doc/ReferenceCounting.html)
!> \param rs_desc the grid descriptor to retain
!> \par History
!>      04.2009 created [Iain Bethune]
!>        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
! **************************************************************************************************
   SUBROUTINE rs_grid_retain_descriptor(rs_desc)
      TYPE(realspace_grid_desc_type), INTENT(INOUT)      :: rs_desc

      CPASSERT(rs_desc%ref_count > 0)
      rs_desc%ref_count = rs_desc%ref_count + 1
   END SUBROUTINE rs_grid_retain_descriptor

! **************************************************************************************************
!> \brief releases the given rs grid (see doc/ReferenceCounting.html)
!> \param rs_grid the rs grid to release
!> \par History
!>      03.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE rs_grid_release(rs_grid)
      TYPE(realspace_grid_type), INTENT(INOUT)           :: rs_grid

      CALL rs_grid_release_descriptor(rs_grid%desc)

      CALL offload_free_buffer(rs_grid%buffer)
      NULLIFY (rs_grid%r)

      IF (ALLOCATED(rs_grid%px)) DEALLOCATE (rs_grid%px)
      IF (ALLOCATED(rs_grid%py)) DEALLOCATE (rs_grid%py)
      IF (ALLOCATED(rs_grid%pz)) DEALLOCATE (rs_grid%pz)
   END SUBROUTINE rs_grid_release

! **************************************************************************************************
!> \brief releases the given rs grid descriptor (see doc/ReferenceCounting.html)
!> \param rs_desc the rs grid descriptor to release
!> \par History
!>      04.2009 created [Iain Bethune]
!>        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
! **************************************************************************************************
   SUBROUTINE rs_grid_release_descriptor(rs_desc)
      TYPE(realspace_grid_desc_type), POINTER            :: rs_desc

      IF (ASSOCIATED(rs_desc)) THEN
         CPASSERT(rs_desc%ref_count > 0)
         rs_desc%ref_count = rs_desc%ref_count - 1
         IF (rs_desc%ref_count == 0) THEN

            CALL pw_grid_release(rs_desc%pw)

            IF (rs_desc%parallel) THEN
               ! release the group communicator
               CALL rs_desc%group%free()

               DEALLOCATE (rs_desc%virtual2real)
               DEALLOCATE (rs_desc%real2virtual)
            END IF

            IF (rs_desc%distributed) THEN
               DEALLOCATE (rs_desc%rank2coord)
               DEALLOCATE (rs_desc%coord2rank)
               DEALLOCATE (rs_desc%lb_global)
               DEALLOCATE (rs_desc%ub_global)
               DEALLOCATE (rs_desc%x2coord)
               DEALLOCATE (rs_desc%y2coord)
               DEALLOCATE (rs_desc%z2coord)
            END IF

            DEALLOCATE (rs_desc)
         END IF
      END IF
      NULLIFY (rs_desc)
   END SUBROUTINE rs_grid_release_descriptor

! **************************************************************************************************
!> \brief emulates the function of an MPI_cart_shift operation, but the shift is
!>        done in virtual coordinates, and the corresponding real ranks are returned
!> \param rs_grid ...
!> \param dir ...
!> \param disp ...
!> \param source ...
!> \param dest ...
!> \par History
!>      04.2009 created [Iain Bethune]
!>        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
! **************************************************************************************************
   PURE SUBROUTINE cart_shift(rs_grid, dir, disp, source, dest)

      TYPE(realspace_grid_type), INTENT(IN)              :: rs_grid
      INTEGER, INTENT(IN)                                :: dir, disp
      INTEGER, INTENT(OUT)                               :: source, dest

      INTEGER, DIMENSION(3)                              :: shift_coords

      shift_coords = rs_grid%desc%virtual_group_coor
      shift_coords(dir) = MODULO(shift_coords(dir) + disp, rs_grid%desc%group_dim(dir))
      dest = rs_grid%desc%virtual2real(rs_grid%desc%coord2rank(shift_coords(1), shift_coords(2), shift_coords(3)))
      shift_coords = rs_grid%desc%virtual_group_coor
      shift_coords(dir) = MODULO(shift_coords(dir) - disp, rs_grid%desc%group_dim(dir))
      source = rs_grid%desc%virtual2real(rs_grid%desc%coord2rank(shift_coords(1), shift_coords(2), shift_coords(3)))

   END SUBROUTINE

! **************************************************************************************************
!> \brief returns the maximum number of points in the local grid of any process
!>        to account for the case where the grid may later be reordered
!> \param desc ...
!> \return ...
!> \par History
!>      10.2011 created [Iain Bethune]
! **************************************************************************************************
   FUNCTION rs_grid_max_ngpts(desc) RESULT(max_ngpts)
      TYPE(realspace_grid_desc_type), INTENT(IN)         :: desc
      INTEGER                                            :: max_ngpts

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rs_grid_max_ngpts'

      INTEGER                                            :: handle, i
      INTEGER, DIMENSION(3)                              :: lb, ub

      CALL timeset(routineN, handle)

      max_ngpts = 0
      IF ((desc%pw%para%mode == PW_MODE_LOCAL) .OR. &
          (ALL(desc%group_dim == 1))) THEN
         CPASSERT(PRODUCT(INT(desc%npts, KIND=int_8)) < HUGE(1))
         max_ngpts = PRODUCT(desc%npts)
      ELSE
         DO i = 0, desc%group_size - 1
            lb = desc%lb_global(:, i)
            ub = desc%ub_global(:, i)
            lb = lb - desc%border*(1 - desc%perd)
            ub = ub + desc%border*(1 - desc%perd)
            CPASSERT(PRODUCT(INT(ub - lb + 1, KIND=int_8)) < HUGE(1))
            max_ngpts = MAX(max_ngpts, PRODUCT(ub - lb + 1))
         END DO
      END IF

      CALL timestop(handle)

   END FUNCTION rs_grid_max_ngpts

! **************************************************************************************************
!> \brief ...
!> \param rs_grid ...
!> \param h_inv ...
!> \param ra ...
!> \param offset ...
!> \param group_size ...
!> \param my_pos ...
!> \return ...
! **************************************************************************************************
   PURE LOGICAL FUNCTION map_gaussian_here(rs_grid, h_inv, ra, offset, group_size, my_pos) RESULT(res)
      TYPE(realspace_grid_type), INTENT(IN)              :: rs_grid
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: h_inv
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra
      INTEGER, INTENT(IN), OPTIONAL                      :: offset, group_size, my_pos

      INTEGER                                            :: dir, lb(3), location(3), tp(3), ub(3)

      res = .FALSE.

      IF (.NOT. ALL(rs_grid%desc%perd == 1)) THEN
         DO dir = 1, 3
            ! bounds of local grid (i.e. removing the 'wings'), if periodic
            tp(dir) = FLOOR(DOT_PRODUCT(h_inv(dir, :), ra)*rs_grid%desc%npts(dir))
            tp(dir) = MODULO(tp(dir), rs_grid%desc%npts(dir))
            IF (rs_grid%desc%perd(dir) .NE. 1) THEN
               lb(dir) = rs_grid%lb_local(dir) + rs_grid%desc%border
               ub(dir) = rs_grid%ub_local(dir) - rs_grid%desc%border
            ELSE
               lb(dir) = rs_grid%lb_local(dir)
               ub(dir) = rs_grid%ub_local(dir)
            END IF
            ! distributed grid, only map if it is local to the grid
            location(dir) = tp(dir) + rs_grid%desc%lb(dir)
         END DO
         IF (ALL(lb(:) <= location(:)) .AND. ALL(location(:) <= ub(:))) THEN
            res = .TRUE.
         END IF
      ELSE
         IF (PRESENT(offset) .AND. PRESENT(group_size) .AND. PRESENT(my_pos)) THEN
            ! not distributed, just a round-robin distribution over the full set of CPUs
            IF (MODULO(offset, group_size) == my_pos) res = .TRUE.
         END IF
      END IF

   END FUNCTION map_gaussian_here

END MODULE realspace_grid_types
