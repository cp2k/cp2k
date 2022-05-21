!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: BSD-3-Clause                                                          !
!--------------------------------------------------------------------------------------------------!

MODULE dbm_api
   USE ISO_C_BINDING, ONLY: &
      C_ASSOCIATED, C_BOOL, C_CHAR, C_DOUBLE, C_F_POINTER, C_INT, C_INT64_T, C_NULL_CHAR, C_NULL_PTR, C_PTR
   USE kinds, ONLY: default_string_length, &
                    dp, &
                    int_8
   USE message_passing, ONLY: mp_cart_rank, &
                              mp_environ

! Uncomment the following line to enable validation.
!#define DBM_VALIDATE_AGAINST_DBCSR
#define DBM_VALIDATE_NBLOCKS_MATCH .TRUE.
#define DBM_VALIDATE_THRESHOLD 5e-10_dp

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
   USE dbcsr_block_access, ONLY: dbcsr_get_block_p, &
                                 dbcsr_put_block, &
                                 dbcsr_reserve_blocks
   USE dbcsr_dist_methods, ONLY: dbcsr_distribution_col_dist, &
                                 dbcsr_distribution_hold, &
                                 dbcsr_distribution_new, &
                                 dbcsr_distribution_release, &
                                 dbcsr_distribution_row_dist
   USE dbcsr_dist_operations, ONLY: dbcsr_get_stored_coordinates
   USE dbcsr_dist_util, ONLY: dbcsr_checksum
   USE dbcsr_iterator_operations, ONLY: dbcsr_iterator_blocks_left, &
                                        dbcsr_iterator_next_block, &
                                        dbcsr_iterator_start, &
                                        dbcsr_iterator_stop
   USE dbcsr_methods, ONLY: dbcsr_col_block_sizes, &
                            dbcsr_get_num_blocks, &
                            dbcsr_get_nze, &
                            dbcsr_mp_release, &
                            dbcsr_release, &
                            dbcsr_row_block_sizes
   USE dbcsr_mp_methods, ONLY: dbcsr_mp_new
   USE dbcsr_multiply_api, ONLY: dbcsr_multiply
   USE dbcsr_operations, ONLY: dbcsr_add, &
                               dbcsr_clear, &
                               dbcsr_copy, &
                               dbcsr_filter, &
                               dbcsr_get_info, &
                               dbcsr_maxabs, &
                               dbcsr_scale, &
                               dbcsr_zero
   USE dbcsr_transformations, ONLY: dbcsr_redistribute
   USE dbcsr_types, ONLY: dbcsr_distribution_obj, &
                          dbcsr_iterator, &
                          dbcsr_mp_obj, &
                          dbcsr_no_transpose, &
                          dbcsr_transpose, &
                          dbcsr_type, &
                          dbcsr_type_no_symmetry, &
                          dbcsr_type_real_8
   USE dbcsr_work_operations, ONLY: dbcsr_create, &
                                    dbcsr_finalize
   USE dbcsr_data_methods, ONLY: dbcsr_scalar

#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbm_api'

   PUBLIC :: dbm_distribution_obj
   PUBLIC :: dbm_distribution_new
   PUBLIC :: dbm_distribution_hold
   PUBLIC :: dbm_distribution_release
   PUBLIC :: dbm_distribution_col_dist
   PUBLIC :: dbm_distribution_row_dist

   PUBLIC :: dbm_iterator
   PUBLIC :: dbm_iterator_start
   PUBLIC :: dbm_iterator_stop
   PUBLIC :: dbm_iterator_num_blocks
   PUBLIC :: dbm_iterator_blocks_left
   PUBLIC :: dbm_iterator_next_block

   PUBLIC :: dbm_type
   PUBLIC :: dbm_release
   PUBLIC :: dbm_create
   PUBLIC :: dbm_create_from_template
   PUBLIC :: dbm_clear
   PUBLIC :: dbm_scale
   PUBLIC :: dbm_get_block_p
   PUBLIC :: dbm_put_block
   PUBLIC :: dbm_reserve_blocks
   PUBLIC :: dbm_filter
   PUBLIC :: dbm_finalize
   PUBLIC :: dbm_multiply
   PUBLIC :: dbm_redistribute
   PUBLIC :: dbm_copy
   PUBLIC :: dbm_add
   PUBLIC :: dbm_maxabs
   PUBLIC :: dbm_zero
   PUBLIC :: dbm_checksum
   PUBLIC :: dbm_get_name
   PUBLIC :: dbm_get_distribution
   PUBLIC :: dbm_get_num_blocks
   PUBLIC :: dbm_get_nze
   PUBLIC :: dbm_get_stored_coordinates
   PUBLIC :: dbm_get_row_block_sizes
   PUBLIC :: dbm_get_col_block_sizes
   PUBLIC :: dbm_get_local_rows
   PUBLIC :: dbm_get_local_cols

   PUBLIC :: dbm_library_init
   PUBLIC :: dbm_library_finalize

   TYPE dbm_distribution_obj
      PRIVATE
      TYPE(C_PTR)                          :: c_ptr = C_NULL_PTR
#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      TYPE(dbcsr_distribution_obj)         :: dbcsr
#endif
   END TYPE dbm_distribution_obj

   TYPE dbm_type
      PRIVATE
      TYPE(C_PTR)                          :: c_ptr = C_NULL_PTR
#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      TYPE(dbcsr_type)                     :: dbcsr
#endif
   END TYPE dbm_type

   TYPE dbm_iterator
      PRIVATE
      TYPE(C_PTR)                          :: c_ptr = C_NULL_PTR
   END TYPE dbm_iterator

CONTAINS

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
! **************************************************************************************************
!> \brief Compates the given DBM matrix against its shadow DBCSR matrics.
!> \param matrix ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE validate(matrix)
      TYPE(dbm_type), INTENT(IN)                         :: matrix

      INTEGER                                            :: col, col_size, col_size_dbcsr, i, j, &
                                                            num_blocks, num_blocks_dbcsr, &
                                                            num_blocks_diff, row, row_size, &
                                                            row_size_dbcsr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: local_cols, local_rows
      LOGICAL                                            :: transposed
      REAL(dp)                                           :: norm2, rel_diff
      REAL(dp), DIMENSION(:, :), POINTER                 :: block, block_dbcsr
      TYPE(C_PTR)                                        :: block_c
      TYPE(dbcsr_iterator)                               :: iter
      INTERFACE
         SUBROUTINE dbm_get_block_p_c(matrix, row, col, block, row_size, col_size) &
            BIND(C, name="dbm_get_block_p")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                        :: matrix
            INTEGER(kind=C_INT), VALUE                :: row
            INTEGER(kind=C_INT), VALUE                :: col
            TYPE(C_PTR)                               :: block
            INTEGER(kind=C_INT)                       :: row_size
            INTEGER(kind=C_INT)                       :: col_size
         END SUBROUTINE dbm_get_block_p_c
      END INTERFACE

      ! Call some getters to run their validation code.
      CALL dbm_get_local_rows(matrix, local_rows)
      CALL dbm_get_local_cols(matrix, local_cols)

      num_blocks_dbcsr = dbcsr_get_num_blocks(matrix%dbcsr)
      num_blocks = dbm_get_num_blocks(matrix)
      num_blocks_diff = ABS(num_blocks - num_blocks_dbcsr)
      IF (num_blocks_diff /= 0) THEN
         WRITE (*, *) "num_blocks mismatch dbcsr:", num_blocks_dbcsr, "new:", num_blocks
         IF (DBM_VALIDATE_NBLOCKS_MATCH) &
            CPABORT("num_blocks mismatch")
      END IF

      IF (DBM_VALIDATE_NBLOCKS_MATCH) THEN
         CPASSERT(dbm_get_nze(matrix) == dbcsr_get_nze(matrix%dbcsr))
      END IF

      ! check all dbcsr blocks
      norm2 = 0.0_dp
      CALL dbcsr_iterator_start(iter, matrix%dbcsr)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row=row, column=col, block=block_dbcsr, &
                                        transposed=transposed, &
                                        row_size=row_size_dbcsr, col_size=col_size_dbcsr)
         CPASSERT(.NOT. transposed)
         CALL dbm_get_block_p_c(matrix=matrix%c_ptr, row=row - 1, col=col - 1, &
                                block=block_c, row_size=row_size, col_size=col_size)

         CPASSERT(row_size == row_size_dbcsr .AND. col_size == col_size_dbcsr)
         IF (SIZE(block_dbcsr) == 0) THEN
            CYCLE
         END IF
         IF (.NOT. C_ASSOCIATED(block_c)) THEN
            CPASSERT(MAXVAL(ABS(block_dbcsr)) < DBM_VALIDATE_THRESHOLD)
            CYCLE
         END IF

         CALL C_F_POINTER(block_c, block, shape=(/row_size, col_size/))
         DO i = 1, row_size
            DO j = 1, col_size
               rel_diff = ABS(block(i, j) - block_dbcsr(i, j))/MAX(1.0_dp, ABS(block_dbcsr(i, j)))
               IF (rel_diff > DBM_VALIDATE_THRESHOLD) THEN
                  WRITE (*, *) "row:", row, "col:", col, "i:", i, "j:", j, "rel_diff:", rel_diff
                  WRITE (*, *) "values dbcsr:", block_dbcsr(i, j), "new:", block(i, j)
                  CPABORT("block value mismatch")
               END IF
            END DO
         END DO
         norm2 = norm2 + SUM(block**2)
         block_dbcsr(:, :) = block(:, :) ! quench numerical noise
      END DO
      CALL dbcsr_iterator_stop(iter)

      ! Can not call dbcsr_get_block_p because it's INTENT(INOUT) :-(

      !! At least check that the norm (=checksum) of excesive blocks is small.
      !TODO: sum norm2 across all mpi ranks.
      !TODO: re-add INTERFACE to dbm_checksum_c, which got removed by prettify.
      !rel_diff = ABS(dbm_checksum_c(matrix%c_ptr) - norm2)/MAX(1.0_dp, norm2)
      !IF (rel_diff > DBM_VALIDATE_THRESHOLD) THEN
      !   WRITE (*, *) "num_blocks dbcsr:", num_blocks_dbcsr, "new:", num_blocks
      !   WRITE (*, *) "norm2: ", norm2
      !   WRITE (*, *) "relative residual norm diff: ", rel_diff
      !   CPABORT("residual norm diff")
      !END IF
   END SUBROUTINE validate

#else

! **************************************************************************************************
!> \brief Dummy for when DBM_VALIDATE_AGAINST_DBCSR is not defined.
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE validate(matrix)
      TYPE(dbm_type), INTENT(IN)                         :: matrix

      MARK_USED(matrix)
   END SUBROUTINE validate
#endif

! **************************************************************************************************
!> \brief Creates a new matrix from given template, reusing dist and row/col_block_sizes.
!> \param matrix ...
!> \param name ...
!> \param template ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_create_from_template(matrix, name, template)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix
      CHARACTER(len=*), INTENT(IN)                       :: name
      TYPE(dbm_type), INTENT(IN)                         :: template

      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: col_block_sizes, row_block_sizes

      ! Store pointers in intermediate variables to workaround a CCE error.
      row_block_sizes => dbm_get_row_block_sizes(template)
      col_block_sizes => dbm_get_col_block_sizes(template)

      CALL dbm_create(matrix, &
                      name=name, &
                      dist=dbm_get_distribution(template), &
                      row_block_sizes=row_block_sizes, &
                      col_block_sizes=col_block_sizes)

   END SUBROUTINE dbm_create_from_template

! **************************************************************************************************
!> \brief Creates a new matrix.
!> \param matrix ...
!> \param name ...
!> \param dist ...
!> \param row_block_sizes ...
!> \param col_block_sizes ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_create(matrix, name, dist, row_block_sizes, col_block_sizes)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix
      CHARACTER(len=*), INTENT(IN)                       :: name
      TYPE(dbm_distribution_obj), INTENT(IN)             :: dist
      INTEGER, CONTIGUOUS, DIMENSION(:), INTENT(IN), &
         POINTER                                         :: row_block_sizes, col_block_sizes

      INTERFACE
         SUBROUTINE dbm_create_c(matrix, dist, name, nrows, ncols, row_sizes, col_sizes) &
            BIND(C, name="dbm_create")
            IMPORT :: C_PTR, C_CHAR, C_INT
            TYPE(C_PTR)                               :: matrix
            TYPE(C_PTR), VALUE                        :: dist
            CHARACTER(kind=C_CHAR), DIMENSION(*)      :: name
            INTEGER(kind=C_INT), VALUE                :: nrows
            INTEGER(kind=C_INT), VALUE                :: ncols
            INTEGER(kind=C_INT), DIMENSION(*)         :: row_sizes
            INTEGER(kind=C_INT), DIMENSION(*)         :: col_sizes
         END SUBROUTINE dbm_create_c
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(matrix%c_ptr))
      CALL dbm_create_c(matrix=matrix%c_ptr, &
                        dist=dist%c_ptr, &
                        name=TRIM(name)//C_NULL_CHAR, &
                        nrows=SIZE(row_block_sizes), &
                        ncols=SIZE(col_block_sizes), &
                        row_sizes=row_block_sizes, &
                        col_sizes=col_block_sizes)
      CPASSERT(C_ASSOCIATED(matrix%c_ptr))

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_create(matrix%dbcsr, name=name, dist=dist%dbcsr, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=row_block_sizes, col_blk_size=col_block_sizes, &
                        data_type=dbcsr_type_real_8)

      CALL validate(matrix)
#endif
   END SUBROUTINE dbm_create

! **************************************************************************************************
!> \brief Needed to be called for DBCSR after blocks where inserted. For DBM it's a no-opt.
!> \param matrix ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_finalize(matrix)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix

      MARK_USED(matrix) ! New implementation does not need finalize.

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_finalize(matrix%dbcsr)
#endif
   END SUBROUTINE dbm_finalize

! **************************************************************************************************
!> \brief Releases a matrix and all its ressources.
!> \param matrix ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_release(matrix)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix

      INTERFACE
         SUBROUTINE dbm_release_c(matrix) &
            BIND(C, name="dbm_release")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                               :: matrix
         END SUBROUTINE dbm_release_c
      END INTERFACE

      CALL dbm_release_c(matrix=matrix%c_ptr)
      matrix%c_ptr = C_NULL_PTR

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_release(matrix%dbcsr)
#endif
   END SUBROUTINE dbm_release

! **************************************************************************************************
!> \brief Copies content of matrix_b into matrix_a.
!>        Matrices must have the same row/col block sizes and distribution.
!> \param matrix_a ...
!> \param matrix_b ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_copy(matrix_a, matrix_b)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix_a
      TYPE(dbm_type), INTENT(IN)                         :: matrix_b

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbm_copy'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE dbm_copy_c(matrix_a, matrix_b) &
            BIND(C, name="dbm_copy")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                               :: matrix_a
            TYPE(C_PTR), VALUE                               :: matrix_b
         END SUBROUTINE dbm_copy_c
      END INTERFACE

      CALL timeset(routineN, handle)
      CALL dbm_copy_c(matrix_a=matrix_a%c_ptr, matrix_b=matrix_b%c_ptr)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_copy(matrix_a%dbcsr, matrix_b%dbcsr)
      CALL validate(matrix_a)
#endif
      CALL timestop(handle)
   END SUBROUTINE dbm_copy

! **************************************************************************************************
!> \brief Copies content of matrix_b into matrix_a. Matrices may have different distributions.
!> \param matrix ...
!> \param redist ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_redistribute(matrix, redist)
      TYPE(dbm_type), INTENT(IN)                         :: matrix
      TYPE(dbm_type), INTENT(INOUT)                      :: redist

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbm_redistribute'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE dbm_redistribute_c(matrix, redist) &
            BIND(C, name="dbm_redistribute")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                               :: matrix
            TYPE(C_PTR), VALUE                               :: redist
         END SUBROUTINE dbm_redistribute_c
      END INTERFACE

      CALL timeset(routineN, handle)
      CALL dbm_redistribute_c(matrix=matrix%c_ptr, redist=redist%c_ptr)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_redistribute(matrix%dbcsr, redist%dbcsr)
      CALL validate(redist)
#endif
      CALL timestop(handle)
   END SUBROUTINE dbm_redistribute

! **************************************************************************************************
!> \brief Looks up a block from given matrics. This routine is thread-safe.
!>        If the block is not found then a null pointer is returned.
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param row_size ...
!> \param col_size ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_get_block_p(matrix, row, col, block, row_size, col_size)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix
      INTEGER, INTENT(IN)                                :: row, col
      REAL(dp), DIMENSION(:, :), INTENT(OUT), POINTER    :: block
      INTEGER, INTENT(OUT), OPTIONAL                     :: row_size, col_size

      INTEGER                                            :: my_col_size, my_row_size
      TYPE(C_PTR)                                        :: block_c
      INTERFACE
         SUBROUTINE dbm_get_block_p_c(matrix, row, col, block, row_size, col_size) &
            BIND(C, name="dbm_get_block_p")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                        :: matrix
            INTEGER(kind=C_INT), VALUE                :: row
            INTEGER(kind=C_INT), VALUE                :: col
            TYPE(C_PTR)                               :: block
            INTEGER(kind=C_INT)                       :: row_size
            INTEGER(kind=C_INT)                       :: col_size
         END SUBROUTINE dbm_get_block_p_c
      END INTERFACE

      CALL dbm_get_block_p_c(matrix=matrix%c_ptr, row=row - 1, col=col - 1, &
                             block=block_c, row_size=my_row_size, col_size=my_col_size)
      IF (C_ASSOCIATED(block_c)) THEN
         CALL C_F_POINTER(block_c, block, shape=(/my_row_size, my_col_size/))
      ELSE
         NULLIFY (block)  ! block not found
      END IF
      IF (PRESENT(row_size)) row_size = my_row_size
      IF (PRESENT(col_size)) col_size = my_col_size
   END SUBROUTINE dbm_get_block_p

! **************************************************************************************************
!> \brief Adds a block to given matrix. This routine is thread-safe.
!>        If block already exist then it gets overwritten (or summed).
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param summation ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_put_block(matrix, row, col, block, summation)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix
      INTEGER, INTENT(IN)                                :: row, col
      REAL(dp), CONTIGUOUS, DIMENSION(:, :), INTENT(IN)  :: block
      LOGICAL, INTENT(IN), OPTIONAL                      :: summation

      LOGICAL                                            :: my_summation
      INTERFACE
         SUBROUTINE dbm_put_block_c(matrix, row, col, summation, block) &
            BIND(C, name="dbm_put_block")
            IMPORT :: C_PTR, C_INT, C_BOOL, C_DOUBLE
            TYPE(C_PTR), VALUE                        :: matrix
            INTEGER(kind=C_INT), VALUE                :: row
            INTEGER(kind=C_INT), VALUE                :: col
            LOGICAL(kind=C_BOOL), VALUE               :: summation
            REAL(kind=C_DOUBLE), DIMENSION(*)         :: block
         END SUBROUTINE dbm_put_block_c
      END INTERFACE

      my_summation = .FALSE.
      IF (PRESENT(summation)) my_summation = summation

      CALL dbm_put_block_c(matrix=matrix%c_ptr, &
                           row=row - 1, col=col - 1, &
                           summation=LOGICAL(my_summation, C_BOOL), &
                           block=block)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_put_block(matrix%dbcsr, row, col, block, summation=summation)
      ! Can not call validate(matrix) because the dbcsr matrix needs to be finalized first.
#endif
   END SUBROUTINE dbm_put_block

! **************************************************************************************************
!> \brief Remove all blocks from given matrix, but does not release the underlying memory.
!> \param matrix ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_clear(matrix)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix

      INTERFACE
         SUBROUTINE dbm_clear_c(matrix) &
            BIND(C, name="dbm_clear")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                               :: matrix
         END SUBROUTINE dbm_clear_c
      END INTERFACE

      CALL dbm_clear_c(matrix=matrix%c_ptr)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_clear(matrix%dbcsr)
      CALL validate(matrix)
#endif
   END SUBROUTINE dbm_clear

! **************************************************************************************************
!> \brief Removes all blocks from the given matrix whose block norm is below the given threshold.
!>        Blocks of size zero are always kept.
!> \param matrix ...
!> \param eps ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_filter(matrix, eps)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix
      REAL(dp), INTENT(IN)                               :: eps

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbm_filter'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE dbm_filter_c(matrix, eps) &
            BIND(C, name="dbm_filter")
            IMPORT :: C_PTR, C_DOUBLE
            TYPE(C_PTR), VALUE                        :: matrix
            REAL(kind=C_DOUBLE), VALUE                :: eps
         END SUBROUTINE dbm_filter_c
      END INTERFACE

      CALL timeset(routineN, handle)
      CALL validate(matrix)
      CALL dbm_filter_c(matrix=matrix%c_ptr, eps=eps)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_filter(matrix%dbcsr, eps)
      CALL validate(matrix)
#endif
      CALL timestop(handle)
   END SUBROUTINE dbm_filter

! **************************************************************************************************
!> \brief Adds given list of blocks efficiently. The blocks will be filled with zeros.
!> \param matrix ...
!> \param rows ...
!> \param cols ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_reserve_blocks(matrix, rows, cols)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix
      INTEGER, DIMENSION(:), INTENT(IN)                  :: rows, cols

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbm_reserve_blocks'

      INTEGER                                            :: handle
      INTEGER(kind=C_INT), DIMENSION(SIZE(rows))         :: cols_c, rows_c
      INTERFACE
         SUBROUTINE dbm_reserve_blocks_c(matrix, nblocks, rows, cols) &
            BIND(C, name="dbm_reserve_blocks")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                        :: matrix
            INTEGER(kind=C_INT), VALUE                :: nblocks
            INTEGER(kind=C_INT), DIMENSION(*)         :: rows
            INTEGER(kind=C_INT), DIMENSION(*)         :: cols
         END SUBROUTINE dbm_reserve_blocks_c
      END INTERFACE

      CALL timeset(routineN, handle)
      CPASSERT(SIZE(rows) == SIZE(cols))
      rows_c = rows - 1
      cols_c = cols - 1

      CALL dbm_reserve_blocks_c(matrix=matrix%c_ptr, &
                                nblocks=SIZE(rows), &
                                rows=rows_c, &
                                cols=cols_c)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_reserve_blocks(matrix%dbcsr, rows, cols)
      CALL validate(matrix)
#endif
      CALL timestop(handle)
   END SUBROUTINE dbm_reserve_blocks

! **************************************************************************************************
!> \brief Multiplies all entries in the given matrix by the given factor alpha.
!> \param matrix ...
!> \param alpha ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_scale(matrix, alpha)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix
      REAL(dp), INTENT(IN)                               :: alpha

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbm_scale'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE dbm_scale_c(matrix, alpha) &
            BIND(C, name="dbm_scale")
            IMPORT :: C_PTR, C_DOUBLE
            TYPE(C_PTR), VALUE                              :: matrix
            REAL(kind=C_DOUBLE), VALUE                      :: alpha
         END SUBROUTINE dbm_scale_c
      END INTERFACE

      CALL timeset(routineN, handle)
      CALL dbm_scale_c(matrix=matrix%c_ptr, alpha=alpha)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_scale(matrix%dbcsr, alpha)
      CALL validate(matrix)
#endif
      CALL timestop(handle)
   END SUBROUTINE dbm_scale

! **************************************************************************************************
!> \brief Sets all blocks in the given matrix to zero.
!> \param matrix ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_zero(matrix)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbm_zero'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE dbm_zero_c(matrix) &
            BIND(C, name="dbm_zero")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                               :: matrix
         END SUBROUTINE dbm_zero_c
      END INTERFACE

      CALL timeset(routineN, handle)
      CALL dbm_zero_c(matrix=matrix%c_ptr)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_zero(matrix%dbcsr)
      CALL validate(matrix)
#endif
      CALL timestop(handle)
   END SUBROUTINE dbm_zero

! **************************************************************************************************
!> \brief Adds matrix_b to matrix_a.
!> \param matrix_a ...
!> \param matrix_b ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_add(matrix_a, matrix_b)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix_a
      TYPE(dbm_type), INTENT(IN)                         :: matrix_b

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbm_add'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE dbm_add_c(matrix_a, matrix_b) &
            BIND(C, name="dbm_add")
            IMPORT :: C_PTR, C_DOUBLE
            TYPE(C_PTR), VALUE                               :: matrix_a
            TYPE(C_PTR), VALUE                               :: matrix_b
         END SUBROUTINE dbm_add_c
      END INTERFACE

      CALL timeset(routineN, handle)
      CALL validate(matrix_a)
      CALL validate(matrix_b)
      CALL dbm_add_c(matrix_a=matrix_a%c_ptr, matrix_b=matrix_b%c_ptr)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_add(matrix_a%dbcsr, matrix_b%dbcsr)
      CALL validate(matrix_a)
#endif
      CALL timestop(handle)
   END SUBROUTINE dbm_add

! **************************************************************************************************
!> \brief Computes matrix product: matrix_c = alpha * matrix_a * matrix_b + beta * matrix_c.
!> \param transa ...
!> \param transb ...
!> \param alpha ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param beta ...
!> \param matrix_c ...
!> \param retain_sparsity ...
!> \param filter_eps ...
!> \param flop ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_multiply(transa, transb, &
                           alpha, matrix_a, matrix_b, beta, matrix_c, &
                           retain_sparsity, filter_eps, flop)
      LOGICAL, INTENT(IN)                                :: transa, transb
      REAL(kind=dp), INTENT(IN)                          :: alpha
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix_a, matrix_b
      REAL(kind=dp), INTENT(IN)                          :: beta
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix_c
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_sparsity
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: filter_eps
      INTEGER(int_8), INTENT(OUT), OPTIONAL              :: flop

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbm_multiply'

      CHARACTER(LEN=1)                                   :: transa_char, transb_char
      INTEGER                                            :: handle
      INTEGER(int_8)                                     :: flop_dbcsr, my_flop
      LOGICAL                                            :: my_retain_sparsity
      REAL(kind=dp)                                      :: my_filter_eps
      INTERFACE
         SUBROUTINE dbm_multiply_c(transa, transb, alpha, &
                                   matrix_a, matrix_b, &
                                   beta, matrix_c, &
                                   retain_sparsity, filter_eps, flop) &
            BIND(C, name="dbm_multiply")
            IMPORT :: C_PTR, C_DOUBLE, C_BOOL, C_INT64_T
            LOGICAL(kind=C_BOOL), VALUE                      :: transa
            LOGICAL(kind=C_BOOL), VALUE                      :: transb
            REAL(kind=C_DOUBLE), VALUE                       :: alpha
            TYPE(C_PTR), VALUE                               :: matrix_a
            TYPE(C_PTR), VALUE                               :: matrix_b
            REAL(kind=C_DOUBLE), VALUE                       :: beta
            TYPE(C_PTR), VALUE                               :: matrix_c
            LOGICAL(kind=C_BOOL), VALUE                      :: retain_sparsity
            REAL(kind=C_DOUBLE), VALUE                       :: filter_eps
            INTEGER(kind=C_INT64_T)                          :: flop
         END SUBROUTINE dbm_multiply_c
      END INTERFACE

      CALL timeset(routineN, handle)

      IF (PRESENT(retain_sparsity)) THEN
         my_retain_sparsity = retain_sparsity
      ELSE
         my_retain_sparsity = .FALSE.
      END IF

      IF (PRESENT(filter_eps)) THEN
         my_filter_eps = filter_eps
      ELSE
         my_filter_eps = 0.0_dp
      END IF

      CALL validate(matrix_a)
      CALL validate(matrix_b)
      CALL validate(matrix_c)
      CALL dbm_multiply_c(transa=LOGICAL(transa, C_BOOL), &
                          transb=LOGICAL(transb, C_BOOL), &
                          alpha=alpha, &
                          matrix_a=matrix_a%c_ptr, &
                          matrix_b=matrix_b%c_ptr, &
                          beta=beta, &
                          matrix_c=matrix_c%c_ptr, &
                          retain_sparsity=LOGICAL(my_retain_sparsity, C_BOOL), &
                          filter_eps=my_filter_eps, &
                          flop=my_flop)

      IF (PRESENT(flop)) THEN
         flop = my_flop
      END IF

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      IF (transa) THEN
         transa_char = dbcsr_transpose
      ELSE
         transa_char = dbcsr_no_transpose
      END IF
      IF (transb) THEN
         transb_char = dbcsr_transpose
      ELSE
         transb_char = dbcsr_no_transpose
      END IF
      CALL dbcsr_multiply(transa=transa_char, transb=transb_char, &
                          alpha=alpha, matrix_a=matrix_a%dbcsr, &
                          matrix_b=matrix_b%dbcsr, beta=beta, matrix_c=matrix_c%dbcsr, &
                          retain_sparsity=retain_sparsity, filter_eps=filter_eps, flop=flop_dbcsr)
      CPASSERT(my_flop == flop_dbcsr)
      CALL validate(matrix_c)
#else
      ! Can not use preprocessor's ifdefs before INTERFACE because it confuses prettify.
      MARK_USED(transa_char)
      MARK_USED(transb_char)
      MARK_USED(flop_dbcsr)
#endif
      CALL timestop(handle)
   END SUBROUTINE dbm_multiply

! **************************************************************************************************
!> \brief Creates an iterator for the blocks of the given matrix. The iteration order is not stable.
!> \param iterator ...
!> \param matrix ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_iterator_start(iterator, matrix)
      TYPE(dbm_iterator), INTENT(OUT)                    :: iterator
      TYPE(dbm_type), INTENT(IN)                         :: matrix

      INTERFACE
         SUBROUTINE dbm_iterator_start_c(iterator, matrix) &
            BIND(C, name="dbm_iterator_start")
            IMPORT :: C_PTR
            TYPE(C_PTR)                               :: iterator
            TYPE(C_PTR), VALUE                        :: matrix
         END SUBROUTINE dbm_iterator_start_c
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(iterator%c_ptr))
      CALL dbm_iterator_start_c(iterator=iterator%c_ptr, matrix=matrix%c_ptr)
      CPASSERT(C_ASSOCIATED(iterator%c_ptr))
      CALL validate(matrix)
   END SUBROUTINE dbm_iterator_start

! **************************************************************************************************
!> \brief Returns number of blocks the iterator will provide to calling thread.
!> \param iterator ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION dbm_iterator_num_blocks(iterator) RESULT(num_blocks)
      TYPE(dbm_iterator), INTENT(IN)                     :: iterator
      INTEGER                                            :: num_blocks

      INTERFACE
         FUNCTION dbm_iterator_num_blocks_c(iterator) &
            BIND(C, name="dbm_iterator_num_blocks")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                        :: iterator
            INTEGER(kind=C_INT)                       :: dbm_iterator_num_blocks_c
         END FUNCTION dbm_iterator_num_blocks_c
      END INTERFACE

      num_blocks = dbm_iterator_num_blocks_c(iterator%c_ptr)
   END FUNCTION dbm_iterator_num_blocks

! **************************************************************************************************
!> \brief Tests whether the given iterator has any block left.
!> \param iterator ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION dbm_iterator_blocks_left(iterator) RESULT(blocks_left)
      TYPE(dbm_iterator), INTENT(IN)                     :: iterator
      LOGICAL                                            :: blocks_left

      INTERFACE
         FUNCTION dbm_iterator_blocks_left_c(iterator) &
            BIND(C, name="dbm_iterator_blocks_left")
            IMPORT :: C_PTR, C_BOOL
            TYPE(C_PTR), VALUE                        :: iterator
            LOGICAL(C_BOOL)                           :: dbm_iterator_blocks_left_c
         END FUNCTION dbm_iterator_blocks_left_c
      END INTERFACE

      blocks_left = dbm_iterator_blocks_left_c(iterator%c_ptr)
   END FUNCTION dbm_iterator_blocks_left

! **************************************************************************************************
!> \brief Returns the next block from the given iterator.
!> \param iterator ...
!> \param row ...
!> \param column ...
!> \param block ...
!> \param row_size ...
!> \param col_size ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_iterator_next_block(iterator, row, column, block, row_size, col_size)
      TYPE(dbm_iterator), INTENT(INOUT)                  :: iterator
      INTEGER, INTENT(OUT)                               :: row, column
      REAL(dp), DIMENSION(:, :), INTENT(OUT), OPTIONAL, &
         POINTER                                         :: block
      INTEGER, INTENT(OUT), OPTIONAL                     :: row_size, col_size

      INTEGER                                            :: col0, my_col_size, my_row_size, row0
      TYPE(C_PTR)                                        :: block_c
      INTERFACE
         SUBROUTINE dbm_iterator_next_block_c(iterator, row, col, block, row_size, col_size) &
            BIND(C, name="dbm_iterator_next_block")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                        :: iterator
            INTEGER(kind=C_INT)                       :: row
            INTEGER(kind=C_INT)                       :: col
            TYPE(C_PTR)                               :: block
            INTEGER(kind=C_INT)                       :: row_size
            INTEGER(kind=C_INT)                       :: col_size
         END SUBROUTINE dbm_iterator_next_block_c
      END INTERFACE

      CALL dbm_iterator_next_block_c(iterator%c_ptr, row=row0, col=col0, block=block_c, &
                                     row_size=my_row_size, col_size=my_col_size)

      CPASSERT(C_ASSOCIATED(block_c))
      IF (PRESENT(block)) CALL C_F_POINTER(block_c, block, shape=(/my_row_size, my_col_size/))
      row = row0 + 1
      column = col0 + 1
      IF (PRESENT(row_size)) row_size = my_row_size
      IF (PRESENT(col_size)) col_size = my_col_size
   END SUBROUTINE dbm_iterator_next_block

! **************************************************************************************************
!> \brief Releases the given iterator.
!> \param iterator ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_iterator_stop(iterator)
      TYPE(dbm_iterator), INTENT(INOUT)                  :: iterator

      INTERFACE
         SUBROUTINE dbm_iterator_stop_c(iterator) &
            BIND(C, name="dbm_iterator_stop")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: iterator
         END SUBROUTINE dbm_iterator_stop_c
      END INTERFACE

      CALL dbm_iterator_stop_c(iterator%c_ptr)
      iterator%c_ptr = C_NULL_PTR
   END SUBROUTINE dbm_iterator_stop

! **************************************************************************************************
!> \brief Computes a checksum of the given matrix.
!> \param matrix ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION dbm_checksum(matrix) RESULT(res)
      TYPE(dbm_type), INTENT(IN)                         :: matrix
      REAL(KIND=dp)                                      :: res

      INTERFACE
         FUNCTION dbm_checksum_c(matrix) &
            BIND(C, name="dbm_checksum")
            IMPORT :: C_PTR, C_DOUBLE
            TYPE(C_PTR), VALUE                        :: matrix
            REAL(C_DOUBLE)                            :: dbm_checksum_c
         END FUNCTION dbm_checksum_c
      END INTERFACE

      CALL validate(matrix)
      res = dbm_checksum_c(matrix%c_ptr)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CPASSERT(ABS(res - dbcsr_checksum(matrix%dbcsr))/MAX(1.0_dp, ABS(res)) < DBM_VALIDATE_THRESHOLD)
#endif
   END FUNCTION dbm_checksum

! **************************************************************************************************
!> \brief Returns the absolute value of the larges element of the entire given matrix.
!> \param matrix ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION dbm_maxabs(matrix) RESULT(res)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix
      REAL(KIND=dp)                                      :: res

      INTERFACE
         FUNCTION dbm_maxabs_c(matrix) &
            BIND(C, name="dbm_maxabs")
            IMPORT :: C_PTR, C_DOUBLE
            TYPE(C_PTR), VALUE                        :: matrix
            REAL(C_DOUBLE)                            :: dbm_maxabs_c
         END FUNCTION dbm_maxabs_c
      END INTERFACE

      CALL validate(matrix)
      res = dbm_maxabs_c(matrix%c_ptr)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CPASSERT(ABS(res - dbcsr_maxabs(matrix%dbcsr))/MAX(1.0_dp, ABS(res)) < DBM_VALIDATE_THRESHOLD)
#endif
   END FUNCTION dbm_maxabs

! **************************************************************************************************
!> \brief Returns the name of the matrix of the given matrix.
!> \param matrix ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION dbm_get_name(matrix) RESULT(res)
      TYPE(dbm_type), INTENT(IN)                         :: matrix
      CHARACTER(len=default_string_length)               :: res

      CHARACTER(LEN=1, KIND=C_CHAR), DIMENSION(:), &
         POINTER                                         :: name_f
      INTEGER                                            :: i
      TYPE(C_PTR)                                        :: name_c
      INTERFACE
         FUNCTION dbm_get_name_c(matrix) BIND(C, name="dbm_get_name")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: matrix
            TYPE(C_PTR)                               :: dbm_get_name_c
         END FUNCTION dbm_get_name_c
      END INTERFACE

      name_c = dbm_get_name_c(matrix%c_ptr)

      CALL C_F_POINTER(name_c, name_f, shape=(/default_string_length/))

      res = ""
      DO i = 1, default_string_length
         IF (name_f(i) == C_NULL_CHAR) EXIT
         res(i:i) = name_f(i)
      END DO

   END FUNCTION dbm_get_name

! **************************************************************************************************
!> \brief Returns the number of local Non-Zero Elements of the given matrix.
!> \param matrix ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   PURE FUNCTION dbm_get_nze(matrix) RESULT(res)
      TYPE(dbm_type), INTENT(IN)                         :: matrix
      INTEGER                                            :: res

      INTERFACE
         PURE FUNCTION dbm_get_nze_c(matrix) &
            BIND(C, name="dbm_get_nze")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE, INTENT(IN)            :: matrix
            INTEGER(C_INT)                            :: dbm_get_nze_c
         END FUNCTION dbm_get_nze_c
      END INTERFACE

      res = dbm_get_nze_c(matrix%c_ptr)

   END FUNCTION dbm_get_nze

! **************************************************************************************************
!> \brief Returns the number of local blocks of the given matrix.
!> \param matrix ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   PURE FUNCTION dbm_get_num_blocks(matrix) RESULT(res)
      TYPE(dbm_type), INTENT(IN)                         :: matrix
      INTEGER                                            :: res

      INTERFACE
         PURE FUNCTION dbm_get_num_blocks_c(matrix) &
            BIND(C, name="dbm_get_num_blocks")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE, INTENT(IN)            :: matrix
            INTEGER(C_INT)                            :: dbm_get_num_blocks_c
         END FUNCTION dbm_get_num_blocks_c
      END INTERFACE

      res = dbm_get_num_blocks_c(matrix%c_ptr)

   END FUNCTION dbm_get_num_blocks

! **************************************************************************************************
!> \brief Returns the row block sizes of the given matrix.
!> \param matrix ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION dbm_get_row_block_sizes(matrix) RESULT(res)
      TYPE(dbm_type), INTENT(IN)                         :: matrix
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: res

      INTEGER                                            :: nrows
      TYPE(C_PTR)                                        :: row_sizes
      INTERFACE
         SUBROUTINE dbm_get_row_sizes_c(matrix, nrows, row_sizes) &
            BIND(C, name="dbm_get_row_sizes")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                        :: matrix
            INTEGER(C_INT)                            :: nrows
            TYPE(C_PTR)                               :: row_sizes
         END SUBROUTINE dbm_get_row_sizes_c
      END INTERFACE

      CALL dbm_get_row_sizes_c(matrix%c_ptr, nrows, row_sizes)
      CALL C_F_POINTER(row_sizes, res, shape=(/nrows/))
      ! TODO: maybe return an ALLOCATABLE
   END FUNCTION dbm_get_row_block_sizes

! **************************************************************************************************
!> \brief Returns the column block sizes of the given matrix.
!> \param matrix ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION dbm_get_col_block_sizes(matrix) RESULT(res)
      TYPE(dbm_type), INTENT(IN)                         :: matrix
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: res

      INTEGER                                            :: ncols
      TYPE(C_PTR)                                        :: col_sizes
      INTERFACE
         SUBROUTINE dbm_get_col_sizes_c(matrix, ncols, col_sizes) &
            BIND(C, name="dbm_get_col_sizes")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                        :: matrix
            INTEGER(C_INT)                            :: ncols
            TYPE(C_PTR)                               :: col_sizes
         END SUBROUTINE dbm_get_col_sizes_c
      END INTERFACE

      CALL dbm_get_col_sizes_c(matrix%c_ptr, ncols, col_sizes)
      CALL C_F_POINTER(col_sizes, res, shape=(/ncols/))
      ! TODO: maybe return an ALLOCATABLE
   END FUNCTION dbm_get_col_block_sizes

! **************************************************************************************************
!> \brief Returns the local row block sizes of the given matrix.
!> \param matrix ...
!> \param local_rows ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_get_local_rows(matrix, local_rows)
      TYPE(dbm_type), INTENT(IN)                         :: matrix
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: local_rows

      INTEGER                                            :: nlocal_rows
      INTEGER, DIMENSION(:), POINTER                     :: local_rows_dbcsr, local_rows_ptr
      TYPE(C_PTR)                                        :: local_rows_c
      INTERFACE
         SUBROUTINE dbm_get_local_rows_c(matrix, nlocal_rows, local_rows) &
            BIND(C, name="dbm_get_local_rows")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                        :: matrix
            INTEGER(C_INT)                            :: nlocal_rows
            TYPE(C_PTR)                               :: local_rows
         END SUBROUTINE dbm_get_local_rows_c
      END INTERFACE

      CALL dbm_get_local_rows_c(matrix%c_ptr, nlocal_rows, local_rows_c)
      CALL C_F_POINTER(local_rows_c, local_rows_ptr, shape=(/nlocal_rows/))
      ALLOCATE (local_rows(nlocal_rows))
      local_rows(:) = local_rows_ptr(:) + 1

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_get_info(matrix%dbcsr, local_rows=local_rows_dbcsr)
      CPASSERT(ALL(local_rows == local_rows_dbcsr))
#else
      MARK_USED(local_rows_dbcsr)
#endif
   END SUBROUTINE dbm_get_local_rows

! **************************************************************************************************
!> \brief Returns the local column block sizes of the given matrix.
!> \param matrix ...
!> \param local_cols ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_get_local_cols(matrix, local_cols)
      TYPE(dbm_type), INTENT(IN)                         :: matrix
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: local_cols

      INTEGER                                            :: nlocal_cols
      INTEGER, DIMENSION(:), POINTER                     :: local_cols_dbcsr, local_cols_ptr
      TYPE(C_PTR)                                        :: local_cols_c
      INTERFACE
         SUBROUTINE dbm_get_local_cols_c(matrix, nlocal_cols, local_cols) &
            BIND(C, name="dbm_get_local_cols")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                        :: matrix
            INTEGER(C_INT)                            :: nlocal_cols
            TYPE(C_PTR)                               :: local_cols
         END SUBROUTINE dbm_get_local_cols_c
      END INTERFACE

      CALL dbm_get_local_cols_c(matrix%c_ptr, nlocal_cols, local_cols_c)
      CALL C_F_POINTER(local_cols_c, local_cols_ptr, shape=(/nlocal_cols/))
      ALLOCATE (local_cols(nlocal_cols))
      local_cols(:) = local_cols_ptr(:) + 1

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_get_info(matrix%dbcsr, local_cols=local_cols_dbcsr)
      CPASSERT(ALL(local_cols == local_cols_dbcsr))
#else
      MARK_USED(local_cols_dbcsr)
#endif
   END SUBROUTINE dbm_get_local_cols

! **************************************************************************************************
!> \brief Returns the MPI rank on which the given block should be stored.
!> \param matrix ...
!> \param row ...
!> \param column ...
!> \param processor ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_get_stored_coordinates(matrix, row, column, processor)
      TYPE(dbm_type), INTENT(IN)                         :: matrix
      INTEGER, INTENT(IN)                                :: row, column
      INTEGER, INTENT(OUT)                               :: processor

      INTEGER                                            :: processor_dbcsr
      INTERFACE
         PURE FUNCTION dbm_get_stored_coordinates_c(matrix, row, col) &
            BIND(C, name="dbm_get_stored_coordinates")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE, INTENT(IN)            :: matrix
            INTEGER(C_INT), VALUE, INTENT(IN)         :: row
            INTEGER(C_INT), VALUE, INTENT(IN)         :: col
            INTEGER(C_INT)                            :: dbm_get_stored_coordinates_c
         END FUNCTION dbm_get_stored_coordinates_c
      END INTERFACE

      processor = dbm_get_stored_coordinates_c(matrix%c_ptr, row=row - 1, col=column - 1)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_get_stored_coordinates(matrix%dbcsr, row, column, processor_dbcsr)
      CPASSERT(processor == processor_dbcsr)
#else
      MARK_USED(processor_dbcsr)
#endif
   END SUBROUTINE dbm_get_stored_coordinates

! **************************************************************************************************
!> \brief Returns the distribution of the given matrix.
!> \param matrix ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION dbm_get_distribution(matrix) RESULT(res)
      TYPE(dbm_type), INTENT(IN)                         :: matrix
      TYPE(dbm_distribution_obj)                         :: res

      INTERFACE
         FUNCTION dbm_get_distribution_c(matrix) BIND(C, name="dbm_get_distribution")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: matrix
            TYPE(C_PTR)                               :: dbm_get_distribution_c
         END FUNCTION dbm_get_distribution_c
      END INTERFACE

      res%c_ptr = dbm_get_distribution_c(matrix%c_ptr)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_get_info(matrix%dbcsr, distribution=res%dbcsr)
#endif

   END FUNCTION dbm_get_distribution

! **************************************************************************************************
!> \brief Creates a new two dimensional distribution.
!> \param dist ...
!> \param mp_comm ...
!> \param row_dist_block ...
!> \param col_dist_block ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_distribution_new(dist, mp_comm, row_dist_block, col_dist_block)
      TYPE(dbm_distribution_obj), INTENT(OUT)            :: dist
      INTEGER, INTENT(IN)                                :: mp_comm
      INTEGER, CONTIGUOUS, DIMENSION(:), INTENT(IN), &
         POINTER                                         :: row_dist_block, col_dist_block

      INTERFACE
         SUBROUTINE dbm_distribution_new_c(dist, fortran_comm, nrows, ncols, row_dist, col_dist) &
            BIND(C, name="dbm_distribution_new")
            IMPORT :: C_PTR, C_CHAR, C_INT
            TYPE(C_PTR)                               :: dist
            INTEGER(kind=C_INT), VALUE                :: fortran_comm
            INTEGER(kind=C_INT), VALUE                :: nrows
            INTEGER(kind=C_INT), VALUE                :: ncols
            INTEGER(kind=C_INT), DIMENSION(*)         :: row_dist
            INTEGER(kind=C_INT), DIMENSION(*)         :: col_dist
         END SUBROUTINE dbm_distribution_new_c
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(dist%c_ptr))
      CALL dbm_distribution_new_c(dist=dist%c_ptr, &
                                  fortran_comm=mp_comm, &
                                  nrows=SIZE(row_dist_block), &
                                  ncols=SIZE(col_dist_block), &
                                  row_dist=row_dist_block, &
                                  col_dist=col_dist_block)
      CPASSERT(C_ASSOCIATED(dist%c_ptr))

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_distribution_new_wrapper(dist, mp_comm, row_dist_block, col_dist_block)
#endif
   END SUBROUTINE dbm_distribution_new

! **************************************************************************************************
!> \brief Helper for creating a new DBCSR distribution. Only needed for DBM_VALIDATE_AGAINST_DBCSR.
!> \param dist ...
!> \param mp_comm ...
!> \param row_dist_block ...
!> \param col_dist_block ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbcsr_distribution_new_wrapper(dist, mp_comm, row_dist_block, col_dist_block)
      TYPE(dbm_distribution_obj), INTENT(INOUT)          :: dist
      INTEGER, INTENT(IN)                                :: mp_comm
      INTEGER, CONTIGUOUS, DIMENSION(:), INTENT(IN), &
         POINTER                                         :: row_dist_block, col_dist_block

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      INTEGER                                            :: mynode, numnodes, pcol, prow
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: pgrid
      INTEGER, DIMENSION(2)                              :: coord, mycoord, npdims
      TYPE(dbcsr_mp_obj)                                 :: mp_env

      ! Create a dbcsr mp environment from communicator
      CALL mp_environ(numnodes, npdims, mycoord, mp_comm)
      CALL mp_environ(numnodes, mynode, mp_comm)
      ALLOCATE (pgrid(0:npdims(1) - 1, 0:npdims(2) - 1))
      DO prow = 0, npdims(1) - 1
         DO pcol = 0, npdims(2) - 1
            coord = (/prow, pcol/)
            CALL mp_cart_rank(mp_comm, coord, pgrid(prow, pcol))
         END DO
      END DO
      CPASSERT(mynode == pgrid(mycoord(1), mycoord(2)))

      CALL dbcsr_mp_new(mp_env, mp_comm, pgrid, mynode, numnodes, mycoord(1), mycoord(2))
      CALL dbcsr_distribution_new(dist=dist%dbcsr, mp_env=mp_env, &
                                  row_dist_block=row_dist_block, col_dist_block=col_dist_block)
      CALL dbcsr_mp_release(mp_env)
#else
      MARK_USED(dist)
      MARK_USED(mp_comm)
      MARK_USED(row_dist_block)
      MARK_USED(col_dist_block)
#endif
   END SUBROUTINE dbcsr_distribution_new_wrapper

! **************************************************************************************************
!> \brief Increases the reference counter of the given distribution.
!> \param dist ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_distribution_hold(dist)
      TYPE(dbm_distribution_obj)                         :: dist

      INTERFACE
         SUBROUTINE dbm_distribution_hold_c(dist) &
            BIND(C, name="dbm_distribution_hold")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: dist
         END SUBROUTINE dbm_distribution_hold_c
      END INTERFACE

      CALL dbm_distribution_hold_c(dist%c_ptr)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_distribution_hold(dist%dbcsr)
#endif
   END SUBROUTINE dbm_distribution_hold

! **************************************************************************************************
!> \brief Decreases the reference counter of the given distribution.
!> \param dist ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_distribution_release(dist)
      TYPE(dbm_distribution_obj)                         :: dist

      INTERFACE
         SUBROUTINE dbm_distribution_release_c(dist) &
            BIND(C, name="dbm_distribution_release")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: dist
         END SUBROUTINE dbm_distribution_release_c
      END INTERFACE

      CALL dbm_distribution_release_c(dist%c_ptr)

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CALL dbcsr_distribution_release(dist%dbcsr)
#endif
   END SUBROUTINE dbm_distribution_release

! **************************************************************************************************
!> \brief Returns the rows of the given distribution.
!> \param dist ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION dbm_distribution_row_dist(dist) RESULT(res)
      TYPE(dbm_distribution_obj), INTENT(IN)             :: dist
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: res

      INTEGER                                            :: nrows
      TYPE(C_PTR)                                        :: row_dist
      INTERFACE
         SUBROUTINE dbm_distribution_row_dist_c(dist, nrows, row_dist) &
            BIND(C, name="dbm_distribution_row_dist")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                        :: dist
            INTEGER(C_INT)                            :: nrows
            TYPE(C_PTR)                               :: row_dist
         END SUBROUTINE dbm_distribution_row_dist_c
      END INTERFACE

      CALL dbm_distribution_row_dist_c(dist%c_ptr, nrows, row_dist)
      CALL C_F_POINTER(row_dist, res, shape=(/nrows/))

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CPASSERT(ALL(res == dbcsr_distribution_row_dist(dist%dbcsr)))
#endif
   END FUNCTION dbm_distribution_row_dist

! **************************************************************************************************
!> \brief Returns the columns of the given distribution.
!> \param dist ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION dbm_distribution_col_dist(dist) RESULT(res)
      TYPE(dbm_distribution_obj), INTENT(IN)             :: dist
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: res

      INTEGER                                            :: ncols
      TYPE(C_PTR)                                        :: col_dist
      INTERFACE
         SUBROUTINE dbm_distribution_col_dist_c(dist, ncols, col_dist) &
            BIND(C, name="dbm_distribution_col_dist")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                        :: dist
            INTEGER(C_INT)                            :: ncols
            TYPE(C_PTR)                               :: col_dist
         END SUBROUTINE dbm_distribution_col_dist_c
      END INTERFACE

      CALL dbm_distribution_col_dist_c(dist%c_ptr, ncols, col_dist)
      CALL C_F_POINTER(col_dist, res, shape=(/ncols/))

#if defined(DBM_VALIDATE_AGAINST_DBCSR)
      CPASSERT(ALL(res == dbcsr_distribution_col_dist(dist%dbcsr)))
#endif
   END FUNCTION dbm_distribution_col_dist

! **************************************************************************************************
!> \brief Initialize DBM library
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_library_init()
      INTERFACE
         SUBROUTINE dbm_library_init_c() BIND(C, name="dbm_library_init")
         END SUBROUTINE dbm_library_init_c
      END INTERFACE

      CALL dbm_library_init_c()

   END SUBROUTINE dbm_library_init

! **************************************************************************************************
!> \brief Finalize DBM library
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbm_library_finalize()
      INTERFACE
         SUBROUTINE dbm_library_finalize_c() BIND(C, name="dbm_library_finalize")
         END SUBROUTINE dbm_library_finalize_c
      END INTERFACE

      CALL dbm_library_finalize_c()

   END SUBROUTINE dbm_library_finalize

END MODULE dbm_api
