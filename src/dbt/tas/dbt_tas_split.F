!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief methods to split tall-and-skinny matrices along longest dimension.
!>        Basically, we are splitting process grid and each subgrid holds its own DBM matrix.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_tas_split
   USE dbt_tas_global,                  ONLY: dbt_tas_distribution
   USE dbt_tas_types,                   ONLY: dbt_tas_distribution_type,&
                                              dbt_tas_split_info
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_cart_create,&
                                              mp_comm_dup,&
                                              mp_comm_free,&
                                              mp_comm_split_direct,&
                                              mp_dims_create,&
                                              mp_environ
   USE util,                            ONLY: sort
#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: &
      dbt_index_global_to_local, &
      dbt_index_local_to_global, &
      colsplit, &
      dbt_tas_get_split_info, &
      dbt_tas_info_hold, &
      dbt_tas_mp_comm, &
      dbt_tas_mp_dims, &
      dbt_tas_release_info, &
      dbt_tas_create_split, &
      dbt_tas_create_split_rows_or_cols, &
      dbt_tas_set_strict_split, &
      group_to_mrowcol, &
      group_to_world_proc_map, &
      rowsplit, &
      world_to_group_proc_map, &
      accept_pgrid_dims, &
      default_nsplit_accept_ratio, &
      default_pdims_accept_ratio

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_tas_split'

   INTEGER, PARAMETER :: rowsplit = 1, colsplit = 2
   REAL(dp), PARAMETER :: default_pdims_accept_ratio = 1.2_dp
   REAL(dp), PARAMETER :: default_nsplit_accept_ratio = 3.0_dp

   INTERFACE dbt_tas_mp_comm
      MODULE PROCEDURE dbt_tas_mp_comm
      MODULE PROCEDURE dbt_tas_mp_comm_from_matrix_sizes
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief split mpi grid by rows or columns
!> \param split_info ...
!> \param mp_comm global mpi communicator with a 2d cartesian grid
!> \param ngroup number of groups
!> \param igroup my group ID
!> \param split_rowcol split rows or columns
!> \param own_comm Whether split_info should own communicator
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_create_split_rows_or_cols(split_info, mp_comm, ngroup, igroup, split_rowcol, own_comm)
      TYPE(dbt_tas_split_info), INTENT(OUT)              :: split_info
      INTEGER, INTENT(IN)                                :: mp_comm
      INTEGER, INTENT(INOUT)                             :: ngroup
      INTEGER, INTENT(IN)                                :: igroup, split_rowcol
      LOGICAL, INTENT(IN), OPTIONAL                      :: own_comm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_tas_create_split_rows_or_cols'

      INTEGER :: handle, igroup_check, iproc, iproc_group, iproc_group_check, mp_comm_group, &
         numproc, numproc_group, numproc_group_check
      INTEGER, DIMENSION(2)                              :: pcoord, pcoord_group, pdims, pdims_group
      LOGICAL                                            :: own_comm_prv, to_assert

      CALL timeset(routineN, handle)

      IF (PRESENT(own_comm)) THEN
         own_comm_prv = own_comm
      ELSE
         own_comm_prv = .FALSE.
      END IF

      IF (own_comm_prv) THEN
         split_info%mp_comm = mp_comm
      ELSE
         CALL mp_comm_dup(mp_comm, split_info%mp_comm)
      END IF

      split_info%igroup = igroup
      split_info%split_rowcol = split_rowcol

      CALL mp_comm_split_direct(mp_comm, mp_comm_group, igroup)

      CALL mp_environ(numproc, iproc, mp_comm)
      CALL mp_environ(numproc, pdims, pcoord, mp_comm)
      split_info%pdims = pdims

      CALL mp_environ(numproc_group, iproc_group, mp_comm_group)

      IF (iproc == 0) THEN
         to_assert = MOD(numproc_group, pdims(MOD(split_rowcol, 2) + 1)) == 0
         CPASSERT(to_assert)
         split_info%pgrid_split_size = numproc_group/pdims(MOD(split_rowcol, 2) + 1)
      END IF
      CALL mp_bcast(split_info%pgrid_split_size, 0, split_info%mp_comm)

      ngroup = (pdims(split_rowcol) + split_info%pgrid_split_size - 1)/split_info%pgrid_split_size
      split_info%ngroup = ngroup
      split_info%group_size = split_info%pgrid_split_size*pdims(MOD(split_rowcol, 2) + 1)

      CALL world_to_group_proc_map(iproc, pdims, split_rowcol, split_info%pgrid_split_size, igroup_check, pdims_group, iproc_group)

      IF (igroup_check .NE. split_info%igroup) THEN
         CPABORT('inconsistent subgroups')
      END IF

      CALL mp_cart_create(mp_comm_group, 2, pdims_group, pcoord_group, split_info%mp_comm_group)

      CALL mp_environ(numproc_group_check, iproc_group_check, split_info%mp_comm_group)

      CPASSERT(iproc_group_check .EQ. iproc_group)

      CALL mp_comm_free(mp_comm_group)

      ALLOCATE (split_info%refcount)
      split_info%refcount = 1

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Create default cartesian process grid that is consistent with default split heuristic
!>        of dbt_tas_create_split
!> \param mp_comm ...
!> \param split_rowcol ...
!> \param nsplit ...
!> \return new communicator
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_mp_comm(mp_comm, split_rowcol, nsplit)
      INTEGER, INTENT(IN)                                :: mp_comm, split_rowcol, nsplit
      INTEGER                                            :: dbt_tas_mp_comm

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbt_tas_mp_comm'

      INTEGER                                            :: handle, iproc, numproc
      INTEGER, DIMENSION(2)                              :: myploc, npdims

      CALL timeset(routineN, handle)

      CALL mp_environ(numproc, iproc, mp_comm)

      npdims = dbt_tas_mp_dims(numproc, split_rowcol, nsplit)

      CALL mp_cart_create(mp_comm, 2, npdims, myploc, dbt_tas_mp_comm)

      CALL timestop(handle)
   END FUNCTION

! **************************************************************************************************
!> \brief Get optimal process grid dimensions consistent with dbt_tas_create_split
!> \param numproc ...
!> \param split_rowcol ...
!> \param nsplit ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_mp_dims(numproc, split_rowcol, nsplit)
      INTEGER, INTENT(IN)                                :: numproc, split_rowcol, nsplit
      INTEGER, DIMENSION(2)                              :: dbt_tas_mp_dims

      INTEGER                                            :: group_size, nsplit_opt
      INTEGER, DIMENSION(2)                              :: group_dims

      nsplit_opt = get_opt_nsplit(numproc, nsplit, split_pgrid=.FALSE.)

      group_size = numproc/nsplit_opt
      group_dims(:) = 0

      CALL mp_dims_create(group_size, group_dims)

      ! here we choose order of group dims s.t. a split factor < nsplit_opt is favoured w.r.t.
      ! optimal subgrid dimensions
      SELECT CASE (split_rowcol)
      CASE (rowsplit)
         group_dims = [MINVAL(group_dims), MAXVAL(group_dims)]
      CASE (colsplit)
         group_dims = [MAXVAL(group_dims), MINVAL(group_dims)]
      END SELECT

      SELECT CASE (split_rowcol)
      CASE (rowsplit)
         dbt_tas_mp_dims(:) = [group_dims(1)*nsplit_opt, group_dims(2)]
      CASE (colsplit)
         dbt_tas_mp_dims(:) = [group_dims(1), group_dims(2)*nsplit_opt]
      END SELECT

   END FUNCTION

! **************************************************************************************************
!> \brief Heuristic to get good split factor for a given process grid OR a given number of processes
!> \param numproc total number of processes or (if split_pgrid) process grid dimension to split
!> \param nsplit Desired split factor
!> \param split_pgrid whether to split process grid
!> \param pdim_nonsplit if split_pgrid: other process grid dimension
!> \return split factor consistent with process grid or number of processes
!> \param
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION get_opt_nsplit(numproc, nsplit, split_pgrid, pdim_nonsplit)
      INTEGER, INTENT(IN)                                :: numproc, nsplit
      LOGICAL, INTENT(IN)                                :: split_pgrid
      INTEGER, INTENT(IN), OPTIONAL                      :: pdim_nonsplit
      INTEGER                                            :: get_opt_nsplit

      INTEGER                                            :: count, count_accept, count_square, lb, &
                                                            minpos, split, ub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nsplit_list, nsplit_list_accept, &
                                                            nsplit_list_square
      INTEGER, DIMENSION(2)                              :: dims_sub

      CPASSERT(nsplit > 0)

      IF (split_pgrid) THEN
         CPASSERT(PRESENT(pdim_nonsplit))
      END IF

      lb = CEILING(REAL(nsplit, dp)/default_nsplit_accept_ratio)
      ub = FLOOR(REAL(nsplit, dp)*default_nsplit_accept_ratio)

      IF (ub < lb) ub = lb

      ALLOCATE (nsplit_list(1:ub - lb + 1), nsplit_list_square(1:ub - lb + 1), nsplit_list_accept(1:ub - lb + 1))
      count = 0
      count_square = 0
      count_accept = 0
      DO split = lb, ub
         IF (MOD(numproc, split) == 0) THEN
            count = count + 1
            nsplit_list(count) = split

            dims_sub = 0
            IF (.NOT. split_pgrid) THEN
               CALL mp_dims_create(numproc/split, dims_sub)
            ELSE
               dims_sub = [numproc/split, pdim_nonsplit]
            END IF

            IF (dims_sub(1) == dims_sub(2)) THEN
               count_square = count_square + 1
               nsplit_list_square(count_square) = split
               count_accept = count_accept + 1
               nsplit_list_accept(count_accept) = split
            ELSEIF (accept_pgrid_dims(dims_sub, relative=.FALSE.)) THEN
               count_accept = count_accept + 1
               nsplit_list_accept(count_accept) = split
            END IF

         END IF
      END DO

      IF (count_square > 0) THEN
         minpos = MINLOC(ABS(nsplit_list_square(1:count_square) - nsplit), DIM=1)
         get_opt_nsplit = nsplit_list_square(minpos)
      ELSEIF (count_accept > 0) THEN
         minpos = MINLOC(ABS(nsplit_list_accept(1:count_accept) - nsplit), DIM=1)
         get_opt_nsplit = nsplit_list_accept(minpos)
      ELSEIF (count > 0) THEN
         minpos = MINLOC(ABS(nsplit_list(1:count) - nsplit), DIM=1)
         get_opt_nsplit = nsplit_list(minpos)
      ELSE
         get_opt_nsplit = nsplit
         DO WHILE (MOD(numproc, get_opt_nsplit) .NE. 0)
            get_opt_nsplit = get_opt_nsplit - 1
         END DO
      END IF

   END FUNCTION

! **************************************************************************************************
!> \brief Derive optimal cartesian process grid from matrix sizes. This ensures optimality for
!>        dense matrices only
!> \param mp_comm ...
!> \param nblkrows total number of block rows
!> \param nblkcols total number of block columns
!> \return MPI communicator
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_mp_comm_from_matrix_sizes(mp_comm, nblkrows, nblkcols) RESULT(mp_comm_new)
      INTEGER, INTENT(IN)                                :: mp_comm
      INTEGER(KIND=int_8), INTENT(IN)                    :: nblkrows, nblkcols
      INTEGER                                            :: mp_comm_new

      INTEGER                                            :: nsplit, split_rowcol

      IF (nblkrows >= nblkcols) THEN
         split_rowcol = rowsplit
         nsplit = INT((nblkrows - 1)/nblkcols + 1)
      ELSE
         split_rowcol = colsplit
         nsplit = INT((nblkcols - 1)/nblkrows + 1)
      END IF

      mp_comm_new = dbt_tas_mp_comm(mp_comm, split_rowcol, nsplit)
   END FUNCTION

! **************************************************************************************************
!> \brief Split Cartesian process grid using a default split heuristic.
!> \param split_info object storing all data corresponding to split, submatrices and parallelization
!> \param mp_comm MPI communicator with associated cartesian grid
!> \param split_rowcol split rows or columns
!> \param nsplit desired split factor, set to 0 if split factor of exactly 1 is required
!> \param own_comm whether split_info should own communicator
!> \param opt_nsplit whether nsplit should be optimized to process grid
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_create_split(split_info, mp_comm, split_rowcol, nsplit, own_comm, opt_nsplit)
      TYPE(dbt_tas_split_info), INTENT(OUT)              :: split_info
      INTEGER, INTENT(IN)                                :: mp_comm, split_rowcol, nsplit
      LOGICAL, INTENT(IN), OPTIONAL                      :: own_comm, opt_nsplit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_tas_create_split'

      INTEGER                                            :: handle, igroup, iproc, nsplit_opt, &
                                                            numproc, pdim_nonsplit, pdim_split
      INTEGER, DIMENSION(2)                              :: pcoord, pdims, pdims_group
      LOGICAL                                            :: opt_nsplit_prv

      CALL timeset(routineN, handle)

      IF (PRESENT(opt_nsplit)) THEN
         opt_nsplit_prv = opt_nsplit
      ELSE
         opt_nsplit_prv = .TRUE.
      END IF

      CPASSERT(nsplit > 0)

      CALL mp_environ(numproc, iproc, mp_comm)
      CALL mp_environ(numproc, pdims, pcoord, mp_comm)

      SELECT CASE (split_rowcol)
      CASE (rowsplit)
         pdim_split = pdims(1)
         pdim_nonsplit = pdims(2)
      CASE (colsplit)
         pdim_split = pdims(2)
         pdim_nonsplit = pdims(1)
      END SELECT

      IF (opt_nsplit_prv) THEN
         nsplit_opt = get_opt_nsplit(pdim_split, nsplit, split_pgrid=.TRUE., pdim_nonsplit=pdim_nonsplit)
      ELSE
         IF (MOD(pdims(split_rowcol), nsplit) .NE. 0) THEN
            CPABORT("Split factor does not divide process grid dimension")
         END IF
         nsplit_opt = nsplit
      END IF

      pdims_group = pdims
      pdims_group(split_rowcol) = pdims_group(split_rowcol)/nsplit_opt

      igroup = pcoord(split_rowcol)/pdims_group(split_rowcol)

      CALL dbt_tas_create_split_rows_or_cols(split_info, mp_comm, nsplit_opt, igroup, split_rowcol, own_comm=own_comm)

      IF (nsplit > 0) THEN
         ALLOCATE (split_info%ngroup_opt, SOURCE=nsplit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Whether to accept proposed process grid dimensions (based on ratio of dimensions)
!> \param dims ...
!> \param relative ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION accept_pgrid_dims(dims, relative)
      INTEGER, DIMENSION(2), INTENT(IN)                  :: dims
      LOGICAL, INTENT(IN)                                :: relative
      LOGICAL                                            :: accept_pgrid_dims

      INTEGER, DIMENSION(2)                              :: dims_opt

      IF (relative) THEN
         dims_opt = 0
         CALL mp_dims_create(PRODUCT(dims), dims_opt)
         accept_pgrid_dims = (MAXVAL(REAL(dims, dp))/MAXVAL(dims_opt) .LT. default_pdims_accept_ratio)
      ELSE
         accept_pgrid_dims = (MAXVAL(REAL(dims, dp))/MINVAL(dims) .LT. default_pdims_accept_ratio**2)
      END IF
   END FUNCTION

! **************************************************************************************************
!> \brief Get info on split
!> \param info ...
!> \param mp_comm communicator (global process grid)
!> \param nsplit split factor
!> \param igroup which group do I belong to
!> \param mp_comm_group subgroup communicator (group-local process grid)
!> \param split_rowcol split rows or columns
!> \param pgrid_offset group-local offset in process grid
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_get_split_info(info, mp_comm, nsplit, igroup, mp_comm_group, split_rowcol, pgrid_offset)
      TYPE(dbt_tas_split_info), INTENT(IN)               :: info
      INTEGER, INTENT(OUT), OPTIONAL                     :: mp_comm, nsplit, igroup, mp_comm_group, &
                                                            split_rowcol
      INTEGER, DIMENSION(2), INTENT(OUT), OPTIONAL       :: pgrid_offset

      IF (PRESENT(mp_comm)) mp_comm = info%mp_comm
      IF (PRESENT(mp_comm_group)) mp_comm_group = info%mp_comm_group
      IF (PRESENT(split_rowcol)) split_rowcol = info%split_rowcol
      IF (PRESENT(igroup)) igroup = info%igroup
      IF (PRESENT(nsplit)) nsplit = info%ngroup

      IF (PRESENT(pgrid_offset)) THEN
         SELECT CASE (info%split_rowcol)
         CASE (rowsplit)
            pgrid_offset(:) = [info%igroup*info%pgrid_split_size, 0]
         CASE (colsplit)
            pgrid_offset(:) = [0, info%igroup*info%pgrid_split_size]
         END SELECT
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param split_info ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_release_info(split_info)
      TYPE(dbt_tas_split_info), INTENT(INOUT)            :: split_info

      LOGICAL                                            :: abort

      abort = .FALSE.

      IF (.NOT. ASSOCIATED(split_info%refcount)) THEN
         abort = .TRUE.
      ELSEIF (split_info%refcount < 1) THEN
         abort = .TRUE.
      END IF

      IF (abort) THEN
         CPABORT("can not destroy non-existing split_info")
      END IF

      split_info%refcount = split_info%refcount - 1

      IF (split_info%refcount == 0) THEN
         CALL mp_comm_free(split_info%mp_comm_group)
         CALL mp_comm_free(split_info%mp_comm)
         DEALLOCATE (split_info%refcount)
      END IF

      split_info%pdims = 0

      IF (ALLOCATED(split_info%ngroup_opt)) DEALLOCATE (split_info%ngroup_opt)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param split_info ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_info_hold(split_info)
      TYPE(dbt_tas_split_info), INTENT(IN)               :: split_info

      INTEGER, POINTER                                   :: ref

      IF (split_info%refcount < 1) THEN
         CPABORT("can not hold non-existing split_info")
      END IF
      ref => split_info%refcount
      ref = ref + 1
   END SUBROUTINE

! **************************************************************************************************
!> \brief map global process info to group
!> \param iproc global process ID
!> \param pdims global process dimensions
!> \param split_rowcol split rows or column
!> \param pgrid_split_size how many process rows/cols per group
!> \param igroup group ID
!> \param pdims_group local process grid dimensions
!> \param iproc_group group local process ID
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE world_to_group_proc_map(iproc, pdims, split_rowcol, pgrid_split_size, igroup, &
                                      pdims_group, iproc_group)
      INTEGER, INTENT(IN)                                :: iproc
      INTEGER, DIMENSION(2), INTENT(IN)                  :: pdims
      INTEGER, INTENT(IN)                                :: split_rowcol, pgrid_split_size
      INTEGER, INTENT(OUT)                               :: igroup
      INTEGER, DIMENSION(2), INTENT(OUT), OPTIONAL       :: pdims_group
      INTEGER, INTENT(OUT), OPTIONAL                     :: iproc_group

      INTEGER, DIMENSION(2)                              :: pcoord, pcoord_group

      IF (PRESENT(iproc_group)) THEN
         CPASSERT(PRESENT(pdims_group))
      END IF

      pcoord = [iproc/pdims(2), MOD(iproc, pdims(2))]

      igroup = pcoord(split_rowcol)/pgrid_split_size

      SELECT CASE (split_rowcol)
      CASE (rowsplit)
         IF (PRESENT(pdims_group)) pdims_group = [pgrid_split_size, pdims(2)]
         IF (PRESENT(iproc_group)) pcoord_group = [MOD(pcoord(1), pgrid_split_size), pcoord(2)]
      CASE (colsplit)
         IF (PRESENT(pdims_group)) pdims_group = [pdims(1), pgrid_split_size]
         IF (PRESENT(iproc_group)) pcoord_group = [pcoord(1), MOD(pcoord(2), pgrid_split_size)]
      END SELECT
      IF (PRESENT(iproc_group)) iproc_group = pcoord_group(1)*pdims_group(2) + pcoord_group(2)
   END SUBROUTINE

! **************************************************************************************************
!> \brief map local process info to global info
!> \param iproc global process id
!> \param pdims global process grid dimensions
!> \param split_rowcol split rows or colum
!> \param pgrid_split_size how many process rows/cols per group
!> \param igroup group ID
!> \param iproc_group local process ID
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE group_to_world_proc_map(iproc, pdims, split_rowcol, pgrid_split_size, &
                                      igroup, iproc_group)
      INTEGER, INTENT(OUT)                               :: iproc
      INTEGER, DIMENSION(2), INTENT(IN)                  :: pdims
      INTEGER, INTENT(IN)                                :: split_rowcol, pgrid_split_size, igroup, &
                                                            iproc_group

      INTEGER, DIMENSION(2)                              :: pcoord, pcoord_group, pdims_group

      SELECT CASE (split_rowcol)
      CASE (rowsplit)
         pdims_group = [pgrid_split_size, pdims(2)]
      CASE (colsplit)
         pdims_group = [pdims(1), pgrid_split_size]
      END SELECT

      pcoord_group = [iproc_group/pdims_group(2), MOD(iproc_group, pdims_group(2))]

      SELECT CASE (split_rowcol)
      CASE (rowsplit)
         pcoord = [igroup*pgrid_split_size + pcoord_group(1), pcoord_group(2)]
      CASE (colsplit)
         pcoord = [pcoord_group(1), igroup*pgrid_split_size + pcoord_group(2)]
      END SELECT
      iproc = pcoord(1)*pdims(2) + pcoord(2)
   END SUBROUTINE

! **************************************************************************************************
!> \brief map group local block index to global matrix index
!> \param info ...
!> \param dist ...
!> \param row_group group local row block index
!> \param column_group group local column block index
!> \param row global block row
!> \param column global block column
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_index_local_to_global(info, dist, row_group, column_group, row, column)
      TYPE(dbt_tas_split_info), INTENT(IN)               :: info
      TYPE(dbt_tas_distribution_type), INTENT(IN)        :: dist
      INTEGER, INTENT(IN), OPTIONAL                      :: row_group, column_group
      INTEGER(KIND=int_8), INTENT(OUT), OPTIONAL         :: row, column

      SELECT CASE (info%split_rowcol)
      CASE (rowsplit)
         ASSOCIATE (rows => dist%local_rowcols)
            IF (PRESENT(row)) row = rows(row_group)
            IF (PRESENT(column)) column = column_group
         END ASSOCIATE
      CASE (colsplit)
         ASSOCIATE (cols => dist%local_rowcols)
            IF (PRESENT(row)) row = row_group
            IF (PRESENT(column)) column = cols(column_group)
         END ASSOCIATE
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief map global block index to group local index
!> \param info ...
!> \param dist ...
!> \param row ...
!> \param column ...
!> \param row_group ...
!> \param column_group ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_index_global_to_local(info, dist, row, column, row_group, column_group)
      TYPE(dbt_tas_split_info), INTENT(IN)               :: info
      TYPE(dbt_tas_distribution_type), INTENT(IN)        :: dist
      INTEGER(KIND=int_8), INTENT(IN), OPTIONAL          :: row, column
      INTEGER, INTENT(OUT), OPTIONAL                     :: row_group, column_group

      SELECT CASE (info%split_rowcol)
      CASE (rowsplit)
         IF (PRESENT(row_group)) row_group = i8_bsearch(dist%local_rowcols, row)
         IF (PRESENT(column_group)) column_group = INT(column)
      CASE (colsplit)
         IF (PRESENT(row_group)) row_group = INT(row)
         IF (PRESENT(column_group)) column_group = i8_bsearch(dist%local_rowcols, column)
      END SELECT

   END SUBROUTINE

! **************************************************************************************************
!> \brief binary search for 8-byte integers
!> \param array ...
!> \param el ...
!> \param l_index ...
!> \param u_index ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION i8_bsearch(array, el, l_index, u_index) RESULT(res)
      INTEGER(KIND=int_8), INTENT(in)                    :: array(:), el
      INTEGER, INTENT(in), OPTIONAL                      :: l_index, u_index
      INTEGER                                            :: res

      INTEGER                                            :: aindex, lindex, uindex

      lindex = 1
      uindex = SIZE(array)
      IF (PRESENT(l_index)) lindex = l_index
      IF (PRESENT(u_index)) uindex = u_index
      DO WHILE (lindex <= uindex)
         aindex = (lindex + uindex)/2
         IF (array(aindex) < el) THEN
            lindex = aindex + 1
         ELSE
            uindex = aindex - 1
         END IF
      END DO
      res = lindex
   END FUNCTION

! **************************************************************************************************
!> \brief maps a process subgroup to matrix rows/columns
!> \param info ...
!> \param rowcol_dist ...
!> \param igroup group ID
!> \param rowcols rows/ columns on this group
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE group_to_mrowcol(info, rowcol_dist, igroup, rowcols)
      TYPE(dbt_tas_split_info), INTENT(IN)               :: info

      CLASS(dbt_tas_distribution), INTENT(IN)                   :: rowcol_dist
      INTEGER, INTENT(IN)                                         :: igroup
      INTEGER(KIND=int_8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: rowcols
      INTEGER, DIMENSION(0:info%pgrid_split_size - 1)             :: nrowcols_group
      INTEGER                                                     :: pcoord, nrowcols, count, pcoord_group
      INTEGER, DIMENSION(:), ALLOCATABLE                          :: sort_indices

      nrowcols_group(:) = 0
      DO pcoord = igroup*info%pgrid_split_size, (igroup + 1)*info%pgrid_split_size - 1
         pcoord_group = pcoord - igroup*info%pgrid_split_size
         nrowcols_group(pcoord_group) = SIZE(rowcol_dist%rowcols(pcoord))
      END DO
      nrowcols = SUM(nrowcols_group)

      ALLOCATE (rowcols(nrowcols))

      count = 0
      DO pcoord = igroup*info%pgrid_split_size, (igroup + 1)*info%pgrid_split_size - 1
         pcoord_group = pcoord - igroup*info%pgrid_split_size
         rowcols(count + 1:count + nrowcols_group(pcoord_group)) = rowcol_dist%rowcols(pcoord)
         count = count + nrowcols_group(pcoord_group)
      END DO

      ALLOCATE (sort_indices(nrowcols))
      CALL sort(rowcols, nrowcols, sort_indices)
   END SUBROUTINE

! **************************************************************************************************
!> \brief freeze current split factor such that it is never changed during multiplication
!> \param info ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_set_strict_split(info)
      TYPE(dbt_tas_split_info), INTENT(INOUT)            :: info

      info%strict_split = [.TRUE., .TRUE.]
   END SUBROUTINE

END MODULE
