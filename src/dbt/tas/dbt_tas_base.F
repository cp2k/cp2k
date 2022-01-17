!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Tall-and-skinny matrices: base routines similar to DBM API,
!>        mostly wrappers around existing DBM routines.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_tas_base
   USE dbm_api,                         ONLY: &
        dbm_clear, dbm_create, dbm_create_from_template, dbm_distribution_col_dist, &
        dbm_distribution_hold, dbm_distribution_new, dbm_distribution_obj, &
        dbm_distribution_release, dbm_distribution_row_dist, dbm_filter, dbm_finalize, &
        dbm_get_block_p, dbm_get_col_block_sizes, dbm_get_distribution, dbm_get_local_cols, &
        dbm_get_local_rows, dbm_get_name, dbm_get_num_blocks, dbm_get_nze, &
        dbm_get_row_block_sizes, dbm_iterator, dbm_iterator_blocks_left, dbm_iterator_next_block, &
        dbm_iterator_start, dbm_iterator_stop, dbm_put_block, dbm_release, dbm_reserve_blocks, &
        dbm_type
   USE dbt_tas_global,                  ONLY: dbt_tas_blk_size_arb,&
                                              dbt_tas_dist_arb,&
                                              dbt_tas_distribution,&
                                              dbt_tas_rowcol_data
   USE dbt_tas_split,                   ONLY: colsplit,&
                                              dbt_index_global_to_local,&
                                              dbt_index_local_to_global,&
                                              dbt_tas_create_split,&
                                              dbt_tas_info_hold,&
                                              dbt_tas_release_info,&
                                              group_to_mrowcol,&
                                              rowsplit
   USE dbt_tas_types,                   ONLY: dbt_tas_distribution_type,&
                                              dbt_tas_iterator,&
                                              dbt_tas_split_info,&
                                              dbt_tas_type
   USE dbt_tas_util,                    ONLY: index_unique
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_cart_rank,&
                                              mp_environ,&
                                              mp_sum
#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_tas_base'

   ! DBM wrappers / interface routines
   PUBLIC :: &
      dbt_tas_blk_sizes, &
      dbt_tas_clear, &
      dbt_tas_copy, &
      dbt_tas_create, &
      dbt_tas_destroy, &
      dbt_tas_distribution_destroy, &
      dbt_tas_distribution_new, &
      dbt_tas_filter, &
      dbt_tas_finalize, &
      dbt_tas_get_block_p, &
      dbt_tas_get_info, &
      dbt_tas_get_num_blocks, &
      dbt_tas_get_nze, &
      dbt_tas_get_nze_total, &
      dbt_tas_get_num_blocks_total, &
      dbt_tas_get_stored_coordinates, &
      dbt_tas_info, &
      dbt_tas_iterator_blocks_left, &
      dbt_tas_iterator_next_block, &
      dbt_tas_iterator_start, &
      dbt_tas_iterator_stop, &
      dbt_tas_nblkcols_local, &
      dbt_tas_nblkcols_total, &
      dbt_tas_nblkrows_local, &
      dbt_tas_nblkrows_total, &
      dbt_tas_nfullrows_total, &
      dbt_tas_nfullcols_total, &
      dbt_tas_put_block, &
      dbt_tas_reserve_blocks, &
      dbt_repl_get_stored_coordinates

   ! conversion routines
   PUBLIC :: &
      dbt_tas_convert_to_dbm, &
      dbt_tas_convert_to_tas

   INTERFACE dbt_tas_create
      MODULE PROCEDURE dbt_tas_create_new
      MODULE PROCEDURE dbt_tas_create_template
   END INTERFACE

   INTERFACE dbt_tas_reserve_blocks
      MODULE PROCEDURE dbt_tas_reserve_blocks_template
      MODULE PROCEDURE dbt_tas_reserve_blocks_index
   END INTERFACE

   INTERFACE dbt_tas_iterator_next_block
      MODULE PROCEDURE dbt_tas_iterator_next_block_d
      MODULE PROCEDURE dbt_tas_iterator_next_block_index
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Create new tall-and-skinny matrix.
!>        Exactly like dbt_create_new but with custom types for row_blk_size and col_blk_size
!>        instead of arrays.
!> \param matrix ...
!> \param name ...
!> \param dist ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param own_dist whether matrix should own distribution
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_create_new(matrix, name, dist, row_blk_size, col_blk_size, own_dist)
      TYPE(dbt_tas_type), INTENT(OUT)                    :: matrix
      CHARACTER(len=*), INTENT(IN)                       :: name
      TYPE(dbt_tas_distribution_type), INTENT(INOUT)     :: dist

      CLASS(dbt_tas_rowcol_data), INTENT(IN)       :: row_blk_size, col_blk_size
      LOGICAL, INTENT(IN), OPTIONAL                  :: own_dist

      TYPE(dbt_tas_split_info)                     :: info

      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS     :: row_blk_size_vec, col_blk_size_vec
      INTEGER                                        :: nrows, ncols, irow, col, icol, row
      CHARACTER(LEN=*), PARAMETER                    :: routineN = 'dbt_tas_create_new'
      INTEGER                                        :: handle

      CALL timeset(routineN, handle)

      CALL dbt_tas_copy_distribution(dist, matrix%dist, own_dist)
      matrix%nblkrows = row_blk_size%nmrowcol
      matrix%nblkcols = col_blk_size%nmrowcol

      CPASSERT(matrix%nblkrows == dist%row_dist%nmrowcol)
      CPASSERT(matrix%nblkcols == dist%col_dist%nmrowcol)

      matrix%nfullrows = row_blk_size%nfullrowcol
      matrix%nfullcols = col_blk_size%nfullrowcol

      ALLOCATE (matrix%row_blk_size, source=row_blk_size)
      ALLOCATE (matrix%col_blk_size, source=col_blk_size)

      info = dbt_tas_info(matrix)

      SELECT CASE (info%split_rowcol)
      CASE (rowsplit)
         matrix%nblkrowscols_split = matrix%nblkrows

         ASSOCIATE (rows => dist%local_rowcols)
            nrows = SIZE(rows)
            ncols = INT(dist%col_dist%nmrowcol)
            ALLOCATE (row_blk_size_vec(nrows))
            ALLOCATE (col_blk_size_vec(ncols))
            DO irow = 1, nrows
               row_blk_size_vec(irow) = row_blk_size%data(rows(irow))
            END DO
            DO col = 1, ncols
               col_blk_size_vec(col) = col_blk_size%data(INT(col, KIND=int_8))
            END DO
         END ASSOCIATE
      CASE (colsplit)
         matrix%nblkrowscols_split = matrix%nblkcols

         ASSOCIATE (cols => dist%local_rowcols)
            ncols = SIZE(cols)
            nrows = INT(dist%row_dist%nmrowcol)
            ALLOCATE (row_blk_size_vec(nrows))
            ALLOCATE (col_blk_size_vec(ncols))
            DO icol = 1, ncols
               col_blk_size_vec(icol) = col_blk_size%data(cols(icol))
            END DO
            DO row = 1, nrows
               row_blk_size_vec(row) = row_blk_size%data(INT(row, KIND=int_8))
            END DO
         END ASSOCIATE
      END SELECT

      CALL dbm_create(matrix=matrix%matrix, &
                      name=name, &
                      dist=dist%dbm_dist, &
                      row_block_sizes=row_blk_size_vec, &
                      col_block_sizes=col_blk_size_vec)

      DEALLOCATE (row_blk_size_vec, col_blk_size_vec)
      matrix%valid = .TRUE.
      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Create matrix from template
!> \param matrix_in ...
!> \param matrix ...
!> \param name ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_create_template(matrix_in, matrix, name)
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix_in
      TYPE(dbt_tas_type), INTENT(OUT)                    :: matrix
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: name

      CALL dbm_create_from_template(matrix%matrix, name=name, template=matrix_in%matrix)
      CALL dbm_finalize(matrix%matrix)

      CALL dbt_tas_copy_distribution(matrix_in%dist, matrix%dist)
      ALLOCATE (matrix%row_blk_size, source=matrix_in%row_blk_size)
      ALLOCATE (matrix%col_blk_size, source=matrix_in%col_blk_size)
      matrix%nblkrows = matrix_in%nblkrows
      matrix%nblkcols = matrix_in%nblkcols
      matrix%nblkrowscols_split = matrix_in%nblkrowscols_split
      matrix%nfullrows = matrix_in%nfullrows
      matrix%nfullcols = matrix_in%nfullcols
      matrix%valid = .TRUE.

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_destroy(matrix)
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix

      CALL dbm_release(matrix%matrix)
      CALL dbt_tas_distribution_destroy(matrix%dist)
      DEALLOCATE (matrix%row_blk_size)
      DEALLOCATE (matrix%col_blk_size)
      matrix%valid = .FALSE.
   END SUBROUTINE

! **************************************************************************************************
!> \brief Copy matrix_a to matrix_b
!> \param matrix_b ...
!> \param matrix_a ...
!> \param summation Whether to sum matrices b = a + b
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_copy(matrix_b, matrix_a, summation)
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix_b
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix_a
      LOGICAL, INTENT(IN), OPTIONAL                      :: summation

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbt_tas_copy'

      INTEGER                                            :: handle
      INTEGER(KIND=int_8)                                :: column, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbt_tas_iterator)                             :: iter

      CALL timeset(routineN, handle)
      CPASSERT(matrix_b%valid)

      IF (PRESENT(summation)) THEN
         IF (.NOT. summation) CALL dbt_tas_clear(matrix_b)
      ELSE
         CALL dbt_tas_clear(matrix_b)
      END IF

      CALL dbt_tas_reserve_blocks(matrix_a, matrix_b)

      CALL dbt_tas_iterator_start(iter, matrix_a)
      DO WHILE (dbt_tas_iterator_blocks_left(iter))
         CALL dbt_tas_iterator_next_block(iter, row, column, block)
         CALL dbt_tas_put_block(matrix_b, row, column, block, summation=summation)
      END DO
      CALL dbt_tas_iterator_stop(iter)
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Make sure that matrix_out has same blocks reserved as matrix_in.
!>         This assumes that both matrices have same number of block rows and block columns.
!> \param matrix_in ...
!> \param matrix_out ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_reserve_blocks_template(matrix_in, matrix_out)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix_in
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix_out

      INTEGER                                            :: iblk, nblk
      INTEGER(KIND=int_8)                                :: column, row
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: col_res, row_res
      TYPE(dbt_tas_iterator)                             :: iter

      nblk = dbt_tas_get_num_blocks(matrix_in)
      ALLOCATE (row_res(nblk), col_res(nblk))
      CALL dbt_tas_iterator_start(iter, matrix_in)
      DO iblk = 1, nblk
         CALL dbt_tas_iterator_next_block(iter, row, column)
         row_res(iblk) = row
         col_res(iblk) = column
      END DO
      CPASSERT(.NOT. dbt_tas_iterator_blocks_left(iter))
      CALL dbt_tas_iterator_stop(iter)

      CALL dbt_tas_reserve_blocks(matrix_out, row_res, col_res)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_finalize(matrix)
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix

      CALL dbm_finalize(matrix%matrix)
   END SUBROUTINE

! **************************************************************************************************
!> \brief create new distribution.
!>        Exactly like dbm_distribution_new but with custom types for row_dist and col_dist
!>        instead of arrays.
!> \param dist ...
!> \param mp_comm ...
!> \param row_dist ...
!> \param col_dist ...
!> \param split_info Strategy of how to split process grid (optional).
!>        If not present a default split heuristic is applied.
!> \param nosplit if .TRUE. don't split process grid (optional)
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_distribution_new(dist, mp_comm, row_dist, col_dist, split_info, nosplit)
      TYPE(dbt_tas_distribution_type), INTENT(OUT)       :: dist
      INTEGER, INTENT(IN)                                :: mp_comm

      CLASS(dbt_tas_distribution), INTENT(IN)        :: row_dist, col_dist
      TYPE(dbt_tas_split_info), INTENT(IN), OPTIONAL :: split_info
      !!
      LOGICAL, INTENT(IN), OPTIONAL                    :: nosplit
      !LOGICAL, INTENT(IN), OPTIONAL                    :: strict_split

      TYPE(dbt_tas_split_info)                       :: split_info_prv

      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS       :: row_dist_vec
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS       :: col_dist_vec
      INTEGER                                          :: nrows, ncols, irow, col, icol, row, &
                                                          split_rowcol, nsplit, handle
      LOGICAL                                          :: opt_nsplit
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_tas_distribution_new'

      CALL timeset(routineN, handle)
      IF (PRESENT(split_info)) THEN
         CALL dbt_tas_info_hold(split_info)
         split_info_prv = split_info
      ELSE
         ! default split heuristic: split into submatrices that have roughly same block dimensions
         IF (row_dist%nmrowcol >= col_dist%nmrowcol) THEN
            split_rowcol = rowsplit
            nsplit = INT((row_dist%nmrowcol - 1)/col_dist%nmrowcol + 1)
         ELSE
            split_rowcol = colsplit
            nsplit = INT((col_dist%nmrowcol - 1)/row_dist%nmrowcol + 1)
         END IF
         opt_nsplit = .TRUE.
         IF (PRESENT(nosplit)) THEN
            IF (nosplit) THEN
               nsplit = 1
               opt_nsplit = .FALSE.
            END IF
         END IF
         CALL dbt_tas_create_split(split_info_prv, mp_comm, split_rowcol, nsplit=nsplit, opt_nsplit=opt_nsplit)
      END IF

      SELECT CASE (split_info_prv%split_rowcol)
      CASE (rowsplit)
         CALL group_to_mrowcol(split_info_prv, row_dist, split_info_prv%igroup, dist%local_rowcols)
         nrows = SIZE(dist%local_rowcols)
         ncols = INT(col_dist%nmrowcol)
         ALLOCATE (row_dist_vec(nrows))
         ALLOCATE (col_dist_vec(ncols))
         DO irow = 1, nrows
            row_dist_vec(irow) = row_dist%dist(dist%local_rowcols(irow)) - split_info_prv%pgrid_split_size*split_info_prv%igroup
         END DO
         DO col = 1, ncols
            col_dist_vec(col) = col_dist%dist(INT(col, KIND=int_8))
         END DO
      CASE (colsplit)
         CALL group_to_mrowcol(split_info_prv, col_dist, split_info_prv%igroup, dist%local_rowcols)
         ncols = SIZE(dist%local_rowcols)
         nrows = INT(row_dist%nmrowcol)
         ALLOCATE (col_dist_vec(ncols))
         ALLOCATE (row_dist_vec(nrows))
         DO icol = 1, ncols
            col_dist_vec(icol) = col_dist%dist(dist%local_rowcols(icol)) - split_info_prv%pgrid_split_size*split_info_prv%igroup
         END DO
         DO row = 1, nrows
            row_dist_vec(row) = row_dist%dist(INT(row, KIND=int_8))
         END DO
      END SELECT

      dist%info = split_info_prv

      CALL dbm_distribution_new(dist%dbm_dist, split_info_prv%mp_comm_group, &
                                row_dist_vec, col_dist_vec)
      DEALLOCATE (row_dist_vec, col_dist_vec)
      ALLOCATE (dist%row_dist, source=row_dist)
      ALLOCATE (dist%col_dist, source=col_dist)

      !IF(PRESENT(strict_split)) dist%strict_split = strict_split

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param dist ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_distribution_destroy(dist)
      TYPE(dbt_tas_distribution_type), INTENT(INOUT)     :: dist

      ! Note: Issue with Cray CCE compiler
      ! commented out the following deallocate statements on polymorphic variables,
      ! these cause segfaults with CCE compiler at a later point

      !IF (ALLOCATED(dist%row_dist)) THEN
      !   DEALLOCATE (dist%row_dist)
      !ENDIF
      !IF (ALLOCATED(dist%col_dist)) THEN
      !   DEALLOCATE (dist%col_dist)
      !ENDIF

      IF (ALLOCATED(dist%local_rowcols)) THEN
         DEALLOCATE (dist%local_rowcols)
      END IF
      CALL dbt_tas_release_info(dist%info)
      CALL dbm_distribution_release(dist%dbm_dist)
   END SUBROUTINE

! **************************************************************************************************
!> \brief As dbt_get_stored_coordinates
!> \param matrix ...
!> \param row global matrix blocked row
!> \param column global matrix blocked column
!> \param processor process ID
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_get_stored_coordinates(matrix, row, column, processor)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8), INTENT(IN)                    :: row, column
      INTEGER, INTENT(OUT)                               :: processor

      INTEGER, DIMENSION(2)                              :: pcoord
      TYPE(dbt_tas_split_info)                           :: info

      pcoord(1) = matrix%dist%row_dist%dist(row)
      pcoord(2) = matrix%dist%col_dist%dist(column)
      info = dbt_tas_info(matrix)

      ! workaround for inefficient mpi_cart_rank
      processor = pcoord(1)*info%pdims(2) + pcoord(2)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Get all processors for a given row/col combination if matrix is replicated on each process
!>        subgroup.
!> \param matrix tall-and-skinny matrix whose DBM submatrices are replicated matrices
!> \param row row of a submatrix
!> \param column column of a submatrix
!> \param processors ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_repl_get_stored_coordinates(matrix, row, column, processors)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER, INTENT(IN)                                :: row, column
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: processors

      INTEGER                                            :: igroup
      INTEGER(KIND=int_8)                                :: col_s, row_s
      INTEGER, DIMENSION(2)                              :: pcoord
      TYPE(dbt_tas_split_info)                           :: info

      row_s = INT(row, KIND=int_8); col_s = INT(column, KIND=int_8)

      info = dbt_tas_info(matrix)
      pcoord(1) = matrix%dist%row_dist%dist(row_s)
      pcoord(2) = matrix%dist%col_dist%dist(col_s)

      DO igroup = 0, info%ngroup - 1
         CALL mp_cart_rank(info%mp_comm, pcoord, processors(igroup + 1))
         SELECT CASE (info%split_rowcol)
         CASE (rowsplit)
            row_s = row_s + dbt_tas_nblkrows_local(matrix)
            pcoord(1) = matrix%dist%row_dist%dist(row_s)
         CASE (colsplit)
            col_s = col_s + dbt_tas_nblkcols_local(matrix)
            pcoord(2) = matrix%dist%col_dist%dist(col_s)
         END SELECT
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief Convert a tall-and-skinny matrix into a normal DBM matrix.
!>        This is not recommended for matrices with a very large dimension.
!> \param matrix_rect ...
!> \param matrix_dbm ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_convert_to_dbm(matrix_rect, matrix_dbm)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix_rect
      TYPE(dbm_type), INTENT(OUT)                        :: matrix_dbm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_tas_convert_to_dbm'

      INTEGER                                            :: handle, nblks_local, rb_count
      INTEGER(KIND=int_8)                                :: col, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nz_cols, nz_rows
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: col_dist_vec, col_size_vec, &
                                                            row_dist_vec, row_size_vec
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbm_distribution_obj)                         :: dist
      TYPE(dbt_tas_iterator)                             :: iter
      TYPE(dbt_tas_split_info)                           :: info

      CALL timeset(routineN, handle)

      info = dbt_tas_info(matrix_rect)

      ALLOCATE (row_dist_vec(matrix_rect%nblkrows))
      ALLOCATE (row_size_vec(matrix_rect%nblkrows))
      ALLOCATE (col_dist_vec(matrix_rect%nblkcols))
      ALLOCATE (col_size_vec(matrix_rect%nblkcols))

      DO row = 1, matrix_rect%nblkrows
         row_dist_vec(row) = matrix_rect%dist%row_dist%dist(row)
         row_size_vec(row) = matrix_rect%row_blk_size%data(row)
      END DO

      DO col = 1, matrix_rect%nblkcols
         col_dist_vec(col) = matrix_rect%dist%col_dist%dist(col)
         col_size_vec(col) = matrix_rect%col_blk_size%data(col)
      END DO

      CALL dbm_distribution_new(dist, info%mp_comm, row_dist_vec, col_dist_vec)
      DEALLOCATE (row_dist_vec, col_dist_vec)

      CALL dbm_create(matrix=matrix_dbm, &
                      name=TRIM(dbm_get_name(matrix_rect%matrix)), &
                      dist=dist, &
                      row_block_sizes=row_size_vec, &
                      col_block_sizes=col_size_vec)

      CALL dbm_distribution_release(dist)

      DEALLOCATE (row_size_vec, col_size_vec)
      nblks_local = dbt_tas_get_num_blocks(matrix_rect)
      CALL dbt_tas_iterator_start(iter, matrix_rect)
      ALLOCATE (nz_rows(nblks_local), nz_cols(nblks_local))
      rb_count = 0
      DO WHILE (dbt_tas_iterator_blocks_left(iter))
         CALL dbt_tas_iterator_next_block(iter, row, col)
         rb_count = rb_count + 1
         nz_rows(rb_count) = INT(row)
         nz_cols(rb_count) = INT(col)
      END DO

      CALL dbm_reserve_blocks(matrix_dbm, nz_rows, nz_cols)

      CALL dbt_tas_iterator_stop(iter)
      CALL dbt_tas_iterator_start(iter, matrix_rect)
      DO WHILE (dbt_tas_iterator_blocks_left(iter))
         CALL dbt_tas_iterator_next_block(iter, row, col, block)
         CALL dbm_put_block(matrix_dbm, INT(row), INT(col), block)
      END DO

      CALL dbt_tas_iterator_stop(iter)

      CALL dbm_finalize(matrix_dbm)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Converts a DBM matrix into the tall-and-skinny matrix type.
!> \param info Strategy of how to split process grid
!> \param matrix_rect ...
!> \param matrix_dbm ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_convert_to_tas(info, matrix_rect, matrix_dbm)
      TYPE(dbt_tas_split_info), INTENT(IN)               :: info
      TYPE(dbt_tas_type), INTENT(OUT)                    :: matrix_rect
      TYPE(dbm_type), INTENT(IN)                         :: matrix_dbm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_tas_convert_to_tas'

      CHARACTER(len=default_string_length)               :: name
      INTEGER                                            :: col, handle, numnodes, row
      INTEGER(KIND=int_8)                                :: nbcols, nbrows
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: col_blk_size, row_blk_size
      INTEGER, DIMENSION(2)                              :: pcoord, pdims
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbm_distribution_obj)                         :: dbm_dist
      TYPE(dbm_iterator)                                 :: iter
      TYPE(dbt_tas_blk_size_arb)                         :: col_blk_size_obj, row_blk_size_obj
      TYPE(dbt_tas_dist_arb)                             :: col_dist_obj, row_dist_obj
      TYPE(dbt_tas_distribution_type)                    :: dist

      NULLIFY (col_blk_size, row_blk_size)
      CALL timeset(routineN, handle)
      CALL mp_environ(numnodes, pdims, pcoord, info%mp_comm)

      name = dbm_get_name(matrix_dbm)
      row_blk_size => dbm_get_row_block_sizes(matrix_dbm)
      col_blk_size => dbm_get_col_block_sizes(matrix_dbm)

      nbrows = SIZE(row_blk_size)
      nbcols = SIZE(col_blk_size)

      dbm_dist = dbm_get_distribution(matrix_dbm)
      row_dist_obj = dbt_tas_dist_arb(dbm_distribution_row_dist(dbm_dist), pdims(1), nbrows)
      col_dist_obj = dbt_tas_dist_arb(dbm_distribution_col_dist(dbm_dist), pdims(2), nbcols)

      row_blk_size_obj = dbt_tas_blk_size_arb(row_blk_size)
      col_blk_size_obj = dbt_tas_blk_size_arb(col_blk_size)

      CALL dbt_tas_distribution_new(dist, info%mp_comm, row_dist_obj, col_dist_obj)

      CALL dbt_tas_create(matrix_rect, TRIM(name)//"_compressed", &
                          dist, row_blk_size_obj, col_blk_size_obj)

      CALL dbm_iterator_start(iter, matrix_dbm)
      DO WHILE (dbm_iterator_blocks_left(iter))
         CALL dbm_iterator_next_block(iter, row, col, block)
         CALL dbt_tas_put_block(matrix_rect, INT(row, KIND=int_8), INT(col, KIND=int_8), block)
      END DO
      CALL dbm_iterator_stop(iter)

      CALL dbt_tas_finalize(matrix_rect)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief As dbm_iterator_start
!> \param iter ...
!> \param matrix_in ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_iterator_start(iter, matrix_in)
      TYPE(dbt_tas_iterator), INTENT(INOUT)              :: iter
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix_in

      CALL dbm_iterator_start(iter%iter, matrix_in%matrix)

      iter%dist = matrix_in%dist
   END SUBROUTINE

! **************************************************************************************************
!> \brief As dbm_iterator_blocks_left
!> \param iter ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_iterator_blocks_left(iter)
      TYPE(dbt_tas_iterator), INTENT(IN)                 :: iter
      LOGICAL                                            :: dbt_tas_iterator_blocks_left

      dbt_tas_iterator_blocks_left = dbm_iterator_blocks_left(iter%iter)
   END FUNCTION

! **************************************************************************************************
!> \brief As dbm_iterator_stop
!> \param iter ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_iterator_stop(iter)
      TYPE(dbt_tas_iterator), INTENT(INOUT)              :: iter

      CALL dbm_iterator_stop(iter%iter)
   END SUBROUTINE

! **************************************************************************************************
!> \brief As dbm_iterator_next_block
!> \param iterator ...
!> \param row global block row
!> \param column global block column
!> \param row_size ...
!> \param col_size ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_iterator_next_block_index(iterator, row, column, row_size, col_size)
      TYPE(dbt_tas_iterator), INTENT(INOUT)              :: iterator
      INTEGER(KIND=int_8), INTENT(OUT)                   :: row, column
      INTEGER, INTENT(OUT), OPTIONAL                     :: row_size, col_size

      INTEGER                                            :: column_group, row_group

      CALL dbm_iterator_next_block(iterator%iter, row=row_group, column=column_group, &
                                   row_size=row_size, col_size=col_size)

      CALL dbt_index_local_to_global(iterator%dist%info, iterator%dist, row_group=row_group, column_group=column_group, &
                                     row=row, column=column)

   END SUBROUTINE

! **************************************************************************************************
!> \brief As dbm_reserve_blocks
!> \param matrix ...
!> \param rows ...
!> \param columns ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_reserve_blocks_index(matrix, rows, columns)
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix
      INTEGER(KIND=int_8), DIMENSION(:), INTENT(IN)      :: rows, columns

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_tas_reserve_blocks_index'

      INTEGER                                            :: handle, icol, irow
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: rowcols_group_u
      INTEGER, DIMENSION(SIZE(rows), 2)                  :: rowcols_group

      CALL timeset(routineN, handle)

      DO irow = 1, SIZE(rows)
         CALL dbt_index_global_to_local(dbt_tas_info(matrix), matrix%dist, row=rows(irow), &
                                        row_group=rowcols_group(irow, 1))
      END DO

      DO icol = 1, SIZE(columns)
         CALL dbt_index_global_to_local(dbt_tas_info(matrix), matrix%dist, column=columns(icol), &
                                        column_group=rowcols_group(icol, 2))
      END DO

      ! TODO: While this was needed for DBCSR, it is no longer necessary for DBM.
      CALL index_unique(rowcols_group, rowcols_group_u) ! make sure that index is unique, not sure
      ! if this is really needed or whether DBM
      ! takes care of duplicate indices

      CALL dbm_reserve_blocks(matrix%matrix, rowcols_group_u(:, 1), rowcols_group_u(:, 2))

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Copy a distribution
!> \param dist_in ...
!> \param dist_out ...
!> \param own_dist Whether distribution should be owned by dist_out
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_copy_distribution(dist_in, dist_out, own_dist)
      TYPE(dbt_tas_distribution_type), INTENT(INOUT)     :: dist_in
      TYPE(dbt_tas_distribution_type), INTENT(OUT)       :: dist_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: own_dist

      LOGICAL                                            :: own_dist_prv

      IF (PRESENT(own_dist)) THEN
         own_dist_prv = own_dist
      ELSE
         own_dist_prv = .FALSE.
      END IF

      IF (.NOT. own_dist_prv) THEN
         CALL dbm_distribution_hold(dist_in%dbm_dist)
         CALL dbt_tas_info_hold(dist_in%info)
      END IF

      dist_out = dist_in
   END SUBROUTINE

! **************************************************************************************************
!> \brief Get block size for a given row & column
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param row_size ...
!> \param col_size ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_blk_sizes(matrix, row, col, row_size, col_size)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8), INTENT(IN)                    :: row, col
      INTEGER, INTENT(OUT)                               :: row_size, col_size

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbt_tas_blk_sizes'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      row_size = matrix%row_blk_size%data(row)
      col_size = matrix%col_blk_size%data(col)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief get info on mpi grid splitting
!> \param matrix ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_info(matrix)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      TYPE(dbt_tas_split_info)                           :: dbt_tas_info

      dbt_tas_info = matrix%dist%info
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_nblkrows_total(matrix) RESULT(nblkrows_total)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8)                                :: nblkrows_total

      nblkrows_total = matrix%nblkrows
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_nfullrows_total(matrix) RESULT(nfullrows_total)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8)                                :: nfullrows_total

      nfullrows_total = matrix%nfullrows
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_nblkcols_total(matrix) RESULT(nblkcols_total)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8)                                :: nblkcols_total

      nblkcols_total = matrix%nblkcols
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_nfullcols_total(matrix) RESULT(nfullcols_total)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8)                                :: nfullcols_total

      nfullcols_total = matrix%nfullcols
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_nblkcols_local(matrix) RESULT(nblkcols_local)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER                                            :: nblkcols_local

      nblkcols_local = SIZE(dbm_get_col_block_sizes(matrix%matrix))
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_nblkrows_local(matrix) RESULT(nblkrows_local)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER                                            :: nblkrows_local

      nblkrows_local = SIZE(dbm_get_row_block_sizes(matrix%matrix))
   END FUNCTION

! **************************************************************************************************
!> \brief As dbt_get_num_blocks: get number of local blocks
!> \param matrix ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION dbt_tas_get_num_blocks(matrix) RESULT(num_blocks)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER                                            :: num_blocks

      num_blocks = dbm_get_num_blocks(matrix%matrix)
   END FUNCTION

! **************************************************************************************************
!> \brief get total number of blocks
!> \param matrix ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_get_num_blocks_total(matrix) RESULT(num_blocks)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8)                                :: num_blocks

      TYPE(dbt_tas_split_info)                           :: info

      info = dbt_tas_info(matrix)
      num_blocks = dbt_tas_get_num_blocks(matrix)
      CALL mp_sum(num_blocks, info%mp_comm)

   END FUNCTION

! **************************************************************************************************
!> \brief As dbt_get_nze: get number of local non-zero elements
!> \param matrix ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION dbt_tas_get_nze(matrix)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER                                            :: dbt_tas_get_nze

      dbt_tas_get_nze = dbm_get_nze(matrix%matrix)

   END FUNCTION

! **************************************************************************************************
!> \brief Get total number of non-zero elements
!> \param matrix ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_tas_get_nze_total(matrix)
      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8)                                :: dbt_tas_get_nze_total

      TYPE(dbt_tas_split_info)                           :: info

      dbt_tas_get_nze_total = dbt_tas_get_nze(matrix)
      info = dbt_tas_info(matrix)
      CALL mp_sum(dbt_tas_get_nze_total, info%mp_comm)
   END FUNCTION

! **************************************************************************************************
!> \brief Clear matrix (erase all data)
!> \param matrix ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_clear(matrix)
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix

      CALL dbm_clear(matrix%matrix)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param nblkrows_total ...
!> \param nblkcols_total ...
!> \param local_rows ...
!> \param local_cols ...
!> \param proc_row_dist ...
!> \param proc_col_dist ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param distribution ...
!> \param name ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_get_info(matrix, &
                               nblkrows_total, nblkcols_total, &
                               local_rows, local_cols, &
                               proc_row_dist, proc_col_dist, &
                               row_blk_size, col_blk_size, distribution, name)

      TYPE(dbt_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8), INTENT(OUT), OPTIONAL         :: nblkrows_total, nblkcols_total
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:), &
         OPTIONAL                                        :: local_rows, local_cols

      CLASS(dbt_tas_distribution), ALLOCATABLE, OPTIONAL, &
         INTENT(OUT)                                                  :: proc_row_dist, proc_col_dist
      CLASS(dbt_tas_rowcol_data), ALLOCATABLE, OPTIONAL, &
         INTENT(OUT)                                                  :: row_blk_size, col_blk_size
      TYPE(dbt_tas_distribution_type), OPTIONAL                     :: distribution
      CHARACTER(len=*), INTENT(OUT), OPTIONAL                         :: name

      TYPE(dbt_tas_split_info)                                      :: info
      INTEGER                                                       :: irow, icol
      INTEGER, ALLOCATABLE, DIMENSION(:)                            :: local_rows_local, local_cols_local

      info = dbt_tas_info(matrix)

      IF (PRESENT(local_rows)) THEN
         CALL dbm_get_local_rows(matrix%matrix, local_rows_local)
         ALLOCATE (local_rows(SIZE(local_rows_local)))
         DO irow = 1, SIZE(local_rows_local)
            CALL dbt_index_local_to_global(info, matrix%dist, row_group=local_rows_local(irow), row=local_rows(irow))
         END DO
      END IF

      IF (PRESENT(local_cols)) THEN
         CALL dbm_get_local_cols(matrix%matrix, local_cols_local)
         ALLOCATE (local_cols(SIZE(local_cols_local)))
         DO icol = 1, SIZE(local_cols_local)
            CALL dbt_index_local_to_global(info, matrix%dist, column_group=local_cols_local(icol), column=local_cols(icol))
         END DO
      END IF

      IF (PRESENT(name)) name = dbm_get_name(matrix%matrix)
      IF (PRESENT(nblkrows_total)) nblkrows_total = dbt_tas_nblkrows_total(matrix)
      IF (PRESENT(nblkcols_total)) nblkcols_total = dbt_tas_nblkcols_total(matrix)
      IF (PRESENT(proc_row_dist)) ALLOCATE (proc_row_dist, SOURCE=matrix%dist%row_dist)
      IF (PRESENT(proc_col_dist)) ALLOCATE (proc_col_dist, SOURCE=matrix%dist%col_dist)
      IF (PRESENT(row_blk_size)) ALLOCATE (row_blk_size, SOURCE=matrix%row_blk_size)
      IF (PRESENT(col_blk_size)) ALLOCATE (col_blk_size, SOURCE=matrix%col_blk_size)
      IF (PRESENT(distribution)) distribution = matrix%dist

   END SUBROUTINE

! **************************************************************************************************
!> \brief As dbm_iterator_next_block
!> \param iterator ...
!> \param row ...
!> \param column ...
!> \param block ...
!> \param row_size ...
!> \param col_size ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_iterator_next_block_d(iterator, row, column, block, row_size, col_size)
      TYPE(dbt_tas_iterator), INTENT(INOUT)              :: iterator
      INTEGER(KIND=int_8), INTENT(OUT)                   :: row, column
      REAL(dp), DIMENSION(:, :), POINTER                 :: block
      INTEGER, INTENT(OUT), OPTIONAL                     :: row_size, col_size

      INTEGER                                            :: column_group, row_group

      CALL dbm_iterator_next_block(iterator%iter, row_group, column_group, block, &
                                   row_size=row_size, col_size=col_size)

      CALL dbt_index_local_to_global(iterator%dist%info, iterator%dist, row_group=row_group, column_group=column_group, &
                                     row=row, column=column)

   END SUBROUTINE

! **************************************************************************************************
!> \brief As dbm_put_block
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param summation ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_put_block(matrix, row, col, block, summation)
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix
      INTEGER(KIND=int_8), INTENT(IN)                    :: row, col
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: block
      LOGICAL, INTENT(IN), OPTIONAL                      :: summation

      INTEGER                                            :: col_group, row_group

      CALL dbt_index_global_to_local(matrix%dist%info, matrix%dist, row=row, column=col, &
                                     row_group=row_group, column_group=col_group)

      CALL dbm_put_block(matrix%matrix, row_group, col_group, block, summation=summation)

   END SUBROUTINE

! **************************************************************************************************
!> \brief As dbm_get_block_p
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param row_size ...
!> \param col_size ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_get_block_p(matrix, row, col, block, row_size, col_size)
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix
      INTEGER(KIND=int_8), INTENT(IN)                    :: row, col
      REAL(dp), DIMENSION(:, :), POINTER                 :: block
      INTEGER, INTENT(OUT), OPTIONAL                     :: row_size, col_size

      INTEGER                                            :: col_group, row_group

      CALL dbt_index_global_to_local(matrix%dist%info, matrix%dist, row=row, column=col, &
                                     row_group=row_group, column_group=col_group)

      CALL dbm_get_block_p(matrix%matrix, row_group, col_group, block, &
                           row_size=row_size, col_size=col_size)

   END SUBROUTINE

! **************************************************************************************************
!> \brief As dbm_filter
!> \param matrix ...
!> \param eps ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_filter(matrix, eps)
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix
      REAL(dp), INTENT(IN)                               :: eps

      CALL dbm_filter(matrix%matrix, eps)

   END SUBROUTINE

END MODULE
