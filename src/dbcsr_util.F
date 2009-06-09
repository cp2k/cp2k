!!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR sparse matrix utility routines
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - Created
! *****************************************************************************
MODULE dbcsr_util

  USE array_types,                     ONLY: array_data,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release,&
                                             array_size
  USE cp_array_i_utils,                ONLY: cp_1d_i_p_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_type,&
                                             cp_blacs_env_write,&
                                             get_blacs_info
  USE cp_fm_struct,                    ONLY: cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE dbcsr_methods,                   ONLY: &
       dbcsr_distribution_col_dist, dbcsr_distribution_hold, &
       dbcsr_distribution_init, dbcsr_distribution_local_cols, &
       dbcsr_distribution_local_rows, dbcsr_distribution_mp, &
       dbcsr_distribution_ncols, dbcsr_distribution_new, &
       dbcsr_distribution_nlocal_cols, dbcsr_distribution_nlocal_rows, &
       dbcsr_distribution_nrows, dbcsr_distribution_processor, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_mp_group, dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_new, &
       dbcsr_mp_new_transposed, dbcsr_mp_npcols, dbcsr_mp_nprows, &
       dbcsr_mp_numnodes, dbcsr_mp_pgrid, dbcsr_mp_release
  USE dbcsr_types,                     ONLY: &
       dbcsr_array_type, dbcsr_distribution_obj, &
       dbcsr_imagedistribution_type, dbcsr_meta_size, dbcsr_mp_obj, &
       dbcsr_num_slots, dbcsr_slot_blk_p, dbcsr_slot_col_i, dbcsr_slot_nblks, &
       dbcsr_slot_nze, dbcsr_slot_row_p, dbcsr_slot_size, dbcsr_type, &
       dbcsr_work_type
  USE distribution_2d_types,           ONLY: distribution_2d_create,&
                                             distribution_2d_get,&
                                             distribution_2d_type,&
                                             distribution_2d_write
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_4,&
                                             int_8
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_allocate,&
                                             mp_alltoall,&
                                             mp_deallocate,&
                                             mp_sum
  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_util'

  INTERFACE addto_array
     MODULE PROCEDURE addto_array_i, addto_array_d
  END INTERFACE

  INTERFACE swap
     MODULE PROCEDURE iswap
  END INTERFACE

  INTERFACE ensure_array_size
     MODULE PROCEDURE ensure_array_size_i, ensure_array_size_d, ensure_array_size_i_2d
  END INTERFACE

  INTERFACE dbcsr_block_transpose
     MODULE PROCEDURE dbcsr_block_transpose_inplace, dbcsr_block_transpose_copy
  END INTERFACE


  PUBLIC :: dbcsr_create, dbcsr_work_create, dbcsr_finalize,&
       dbcsr_create_dist2d_r, dbcsr_create_dist2d_r_unrot,&
       dbcsr_create_dist_r_rot, dbcsr_create_dist_r_unrot,&
       dbcsr_create_dist_l_rot, dbcsr_create_dist_l_unrot,&
       dbcsr_from_fm, dbcsr_from_sm, dbcsr_get_info,&
       dbcsr_checksum, dbcsr_verify_matrix,&
       dbcsr_printmat, dbcsr_print2dmat, dbcsr_repoint_index,&
       dbcsr_pack_meta, dbcsr_unpack_i8_2i4, dbcsr_init,&
       dbcsr_create_dist2d_transposed, dbcsr_new_transposed,&
       dbcsr_transpose_distribution
  PUBLIC :: dbcsr_addto_index_array,&
       addto_array, swap, ensure_array_size, dbcsr_ensure_data_size,&
       dbcsr_set_debug,&
       dbcsr_destroy, dbcsr_work_destroy,&
       checker_tr,&
       dbcsr_get_lookup, dbcsr_get_block_info,&
       dbcsr_print, dbcsr_normalize, dbcsr_order, add_work_coordinate,&
       dbcsr_destroy_array, dbcsr_copy_block_sizes, dbcsr_block_transpose
  PUBLIC :: dbcsr_dist2d_to_dist, dbcsr_create_image_dist_left,&
       dbcsr_create_image_dist, dbcsr_destroy_image_dist
  PUBLIC :: dbcsr_set_diagonal, dbcsr_extract_column



#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL dbcsr_set_debug(dbg, info=info)

  LOGICAL, PARAMETER :: bcsr_debug =   .TRUE.
  LOGICAL, PARAMETER :: bcsr_info =    .FALSE.
  LOGICAL, PARAMETER :: bcsr_verbose = .FALSE.

CONTAINS

! *****************************************************************************
!> \brief Determines whether a transpose must be applied
!> \param[in] row   The absolute matrix row.
!> \param[in] column          The absolute matrix column.
! *****************************************************************************
  ELEMENTAL FUNCTION checker_tr(row, column) RESULT(transpose)
    INTEGER, INTENT(IN)                      :: row, column
    LOGICAL                                  :: transpose

!transpose = MOD(column+row,2).NE.0 .EQV. column.GE.row

    transpose = BTEST(column+row, 0) .EQV. column.GE.row
  END FUNCTION checker_tr


! *****************************************************************************
!> \brief Sets the correct source matrix, row, column and possible data
!>        transposition for a given matrix and matrix logical row and
!>        column.
!> \param[in] matrix          DBCSR matrix
!> \param[in,out] row         input is logical row, output is lookup row
!> \param[in,out] column      input is logical column, output is lookup
!>                            column
!> \param[out] transpose      whether the data must be transposed
!> \param[out] processor      (optional) returns the processor on which this
!>                            block resides
! *****************************************************************************
  PURE SUBROUTINE dbcsr_get_lookup(matrix, row, column, transpose, processor)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(INOUT)                   :: row, column
    LOGICAL, INTENT(OUT)                     :: transpose
    INTEGER, INTENT(OUT), OPTIONAL           :: processor

    INTEGER                                  :: tmp
    LOGICAL                                  :: straight

!   ---------------------------------------------------------------------------

    straight = matrix%matrix_type .NE. 'N'
    straight = (.NOT.&
         (matrix%matrix_type.EQ.'S'.OR.matrix%matrix_type.EQ.'A'))&
         .OR.&
         (straight .EQV. .NOT.checker_tr(row, column))
    transpose = .NOT. straight
    IF (.NOT. straight) THEN
       tmp = row
       row = column
       column = tmp
    ENDIF
    IF (PRESENT (processor)) THEN
       processor = dbcsr_distribution_processor (matrix%dist, row, column)
    END IF
  END SUBROUTINE dbcsr_get_lookup


! *****************************************************************************
!> \brief Returns the index to a queried block, given a real (stored) row and
!>        column
!> \param[in] matrix          bcsr matrix
!> \param[in] row             input is logical row, output is lookup row
!> \param[in] column          input is logical column, output is lookup
!>                            column
!> \param[out] block_number   returns the block number of the row
!>                            and column; 0 if nonexistent
!> \param[out] data_offset    (optional) data offset for the block; 0 if
!>                            nonexistent
! *****************************************************************************
  PURE SUBROUTINE dbcsr_get_block_info(matrix, row, column,&
       block_number, data_offset)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(IN)                      :: row, column
    INTEGER, INTENT(OUT)                     :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: data_offset

    INTEGER                                  :: blk
    LOGICAL                                  :: found

!   ---------------------------------------------------------------------------
!CALL timeset("dbcsr_get_block_info", timing_handle)

    found = .FALSE.
    blk = matrix%row_p(row)+1
    DO WHILE (.NOT. found)
       IF (blk .GT. matrix%row_p(row+1)) THEN
          blk = 0
          EXIT
       ENDIF
       IF (matrix%col_i(blk) .EQ. column) EXIT
       IF (matrix%col_i(blk) .GT. column) THEN
          blk = 0
          EXIT
       ENDIF
       blk = blk + 1
    ENDDO
    block_number = blk
    IF (PRESENT (data_offset)) THEN
       IF (blk .GT. 0) THEN
          data_offset = matrix%blk_p(blk)
       ELSE
          data_offset = 0
       ENDIF
    ENDIF
    !CALL timestop (timing_handle)
  END SUBROUTINE dbcsr_get_block_info


! *****************************************************************************
!> \brief Initializes a DBCSR matrix but does not create it.
!> \param[out] matrix         uninitialized matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_init(matrix, error)
    TYPE(dbcsr_type), INTENT(OUT)            :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    matrix%initialized = .TRUE.
    matrix%valid = .FALSE.
    ! Nullifies all pointers.
    NULLIFY (matrix%DATA, matrix%index, matrix%row_p, matrix%col_i,&
         matrix%blk_p)
    CALL dbcsr_distribution_init (matrix%dist)
    NULLIFY (matrix%w%DATA, matrix%w%row_i, matrix%w%col_i, matrix%w%blk_p)
  END SUBROUTINE dbcsr_init

! *****************************************************************************
!> \brief Creates a matrix, allocating the essentials.
!> \par The matrix itself is allocated, as well as the essential parts of
!>      the index. When passed the nze argument, the data is also allocated
!>      to that size.
!> \param[in,out] matrix      new matrix
!> \param[in] dist            distribution_2d distribution
!> \param[in] matrix_type     'N' for normal, 'T' for transposed, 'S' for
!>                            symmetric, and 'A' for antisymmetric
!> \param[in] nblks           number of blocks
!> \param[in] nze             number of elements
!> \param[in] special         (optional) allocate indices and data using
!>                            special memory
!> \param[in] reuse           (optional) reuses an existing matrix, default
!>                            is to create a fresh one
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create(matrix, name, dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, special, reuse, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    CHARACTER(len=*), INTENT(IN)             :: name
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    CHARACTER, INTENT(IN)                    :: matrix_type
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN)                      :: nblks, nze
    LOGICAL, INTENT(IN), OPTIONAL            :: special, reuse
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: mt
    INTEGER                                  :: stat
    INTEGER, DIMENSION(dbcsr_meta_size)      :: new_meta
    LOGICAL                                  :: hijack

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    DBGI 'Allocating matrix', name
    hijack = .FALSE.
    IF (PRESENT (reuse)) hijack = reuse
    ! Why would you want to not reuse an existing matrix?
    CALL cp_assert (hijack .EQV. matrix%valid, cp_warning_level,&
         cp_caller_error, routineP, "Reuse and matrix validity do not match.",&
         error=error)
    matrix%valid = .TRUE.
    matrix%name = name
    IF (.NOT. hijack) THEN
       NULLIFY(matrix%data)
       NULLIFY(matrix%index)
    ENDIF
    NULLIFY(matrix%row_p, matrix%col_i, matrix%blk_p)
    matrix%row_blk_size = row_blk_size
    CALL array_hold (matrix%row_blk_size)
    matrix%col_blk_size = col_blk_size
    CALL array_hold (matrix%col_blk_size)
    matrix%dist = dist
    CALL dbcsr_distribution_hold (matrix%dist)
    CALL meta_from_dist (new_meta, dist, row_blk_size, col_blk_size,&
         error=error)
!    IF (PRESENT (meta)) THEN
!       DO stat = 3, 10
!318       FORMAT ("@",I2,": ",I7," vs. ",I7)
!          IF (new_meta(stat) .NE. meta(stat)) &
!               WRITE(*,318)stat, new_meta(stat), meta(stat)
!          CALL cp_assert (new_meta(stat) .EQ. meta(stat), cp_warning_level,&
!               cp_wrong_args_error, routineP, "Meta does not match dist2.",&
!               error=error)
!       ENDDO
!    ENDIF
    matrix%nblkrows_total  = new_meta(3)
    matrix%nblkcols_total  = new_meta(4)
    matrix%nfullrows_total = new_meta(5)
    matrix%nfullcols_total = new_meta(6)
    matrix%nblkrows_local  = new_meta(7)
    matrix%nblkcols_local  = new_meta(8)
    matrix%nfullrows_local = new_meta(9)
    matrix%nfullcols_local = new_meta(10)
    matrix%nblks = nblks
    matrix%nze = nze
    matrix%special_memory = .FALSE.
    IF (PRESENT (special)) matrix%special_memory = special
    !
    DBG 'special allocation',matrix%special_memory
    IF (.NOT. hijack) THEN
       IF (matrix%special_memory) THEN
          CALL mp_allocate(matrix%data, nze, stat=stat)
          IF (stat /= 0) CALL stop_memory (routineP, 'matrix%data', nze)
          CALL mp_allocate(matrix%index, dbcsr_num_slots, stat=stat)
          IF (stat /= 0) CALL stop_memory (routineP, 'matrix%data', nze)
       ELSE
          ALLOCATE(matrix%data(nze), stat=stat)
          IF (stat /= 0) CALL stop_memory (routineP, 'matrix%data', nze)
          ALLOCATE(matrix%index(dbcsr_num_slots), stat=stat)
          IF (stat /= 0) CALL stop_memory (routineP, 'matrix%data', nze)
       ENDIF
    ENDIF
    matrix%index(:) = 0
    matrix%index(dbcsr_slot_size) = dbcsr_num_slots
    !
    SELECT CASE (matrix_type)
    CASE ('N', 'n')
       mt = 'N'
    CASE ('T', 't')
       mt = 'T'
    CASE ('S', 's')
       mt = 'S'
    CASE DEFAULT
       CALL cp_assert(.FALSE., cp_failure_level,&
            cp_wrong_args_error, routineP, "Invalid matrix type.", error)
    END SELECT
    matrix%matrix_type = mt
    matrix%w%lastblk = 0
    matrix%w%datasize = 0
    DBGV ' data allocated', ASSOCIATED(matrix%data)
    DBGV ' index allocated', ASSOCIATED(matrix%index)
  END SUBROUTINE dbcsr_create


! *****************************************************************************
!> \brief Gets information about a matrix
!> \param[in] matrix          matrix to query
!> \param[out] nblkrows_total
!> \param[out] nblkcols_total
!> \param[out] nfullrows_total
!> \param[out] nfullcols_total
!> \param[out] nblkrows_local
!> \param[out] nblkrows_local
!> \param[out] nfullrows_local
!> \param[out] nfullrows_local
!> \param[out] my_prow, my_pcol
!> \param[out] local_rows, local_cols
!> \param[out] proc_row_dist, proc_col_dist
!> \param[out] row_blk_size
!> \param[out] col_blk_size
!> \param[out] meta           the metadata on the matrix as an array
!> \param[out] distribution   the data distribution of the matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_get_info(matrix, nblkrows_total, nblkcols_total,&
       nfullrows_total, nfullcols_total,&
       nblkrows_local, nblkcols_local,&
       nfullrows_local, nfullcols_local,&
       my_prow, my_pcol,&
       local_rows, local_cols, proc_row_dist, proc_col_dist,&
       row_blk_size, col_blk_size, meta, distribution,&
       error)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(OUT), OPTIONAL :: nblkrows_total, nblkcols_total, &
      nfullrows_total, nfullcols_total, nblkrows_local, nblkcols_local, &
      nfullrows_local, nfullcols_local, my_prow, my_pcol
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: local_rows, local_cols, &
                                                proc_row_dist, proc_col_dist
    TYPE(array_i1d_obj), INTENT(OUT), &
      OPTIONAL                               :: row_blk_size, col_blk_size
    INTEGER, DIMENSION(dbcsr_meta_size), &
      INTENT(OUT), OPTIONAL                  :: meta
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT), OPTIONAL                  :: distribution
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_info', &
      routineP = moduleN//':'//routineN

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL cp_assert (matrix%valid, cp_warning_level, cp_caller_error,&
         routineP,"Invalid matrix", error=error)
    IF (PRESENT (nblkrows_total)) nblkrows_total = matrix%nblkrows_total
    IF (PRESENT (nblkcols_total)) nblkcols_total = matrix%nblkcols_total
    IF (PRESENT (nfullrows_total)) nfullrows_total = matrix%nfullrows_total
    IF (PRESENT (nfullcols_total)) nfullcols_total = matrix%nfullcols_total
    IF (PRESENT (nblkrows_local)) nblkrows_local = matrix%nblkrows_local
    IF (PRESENT (nblkcols_local)) nblkcols_local = matrix%nblkcols_local
    IF (PRESENT (nfullrows_local)) nfullrows_local = matrix%nfullrows_local
    IF (PRESENT (nfullcols_local)) nfullcols_local = matrix%nfullcols_local
    IF (PRESENT (row_blk_size)) row_blk_size = matrix%row_blk_size
    IF (PRESENT (col_blk_size)) col_blk_size = matrix%col_blk_size
    IF (PRESENT (meta)) CALL dbcsr_pack_meta (matrix, meta, error)
    IF (PRESENT (distribution)) distribution = matrix%dist
    !IF (PRESENT (dist2d)) dist2d => matrix%dist2d
!    CALL distribution_2d_get(dist2d, row_distribution=prow_dist,&
!         col_distribution=pcol_dist, error=error, blacs_env=blacs_env,&
!         n_flat_local_rows=nblkrows_local, n_flat_local_cols=nblkcols_local,&
!         n_row_distribution=nblkrows_total, n_col_distribution=nblkcols_total,&
!         flat_local_rows=blkrows_local, flat_local_cols=blkcols_local)
    IF (PRESENT (local_rows)) &
         !CALL distribution_2d_get(matrix%dist2d, flat_local_rows=local_rows,&
         !error=error)
         local_rows => array_data (dbcsr_distribution_local_rows (matrix%dist))
    IF (PRESENT (local_cols)) &
         !CALL distribution_2d_get(matrix%dist2d, flat_local_cols=local_cols,&
         !error=error)
         local_cols => array_data (dbcsr_distribution_local_cols (matrix%dist))
    IF (PRESENT (proc_row_dist)) &
         !CALL distribution_2d_get(matrix%dist2d,&
         !row_distribution = proc_row_dist, error=error)
         proc_row_dist => array_data (dbcsr_distribution_row_dist (matrix%dist))
    IF (PRESENT (proc_col_dist)) &
         !CALL distribution_2d_get(matrix%dist2d,&
         !col_distribution = proc_col_dist, error=error)
         proc_col_dist => array_data (dbcsr_distribution_col_dist (matrix%dist))
    !IF (PRESENT (my_pcol) .OR. PRESENT (my_prow)) THEN
    !   CALL distribution_2d_get(matrix%dist2d, blacs_env=blacs_env, error=error)
    !   !CALL get_blacs_info(blacs_env, para_env=para_env, blacs2mpi=blacs2mpi,&
    !   !  mpi2blacs=mpi2blacs, my_process_row=myprow, my_process_column=mypcol,&
    !   !  number_of_process_rows=nprows, number_of_process_columns=npcols)
    !   IF (PRESENT (my_prow)) &
    !        CALL get_blacs_info (blacs_env, my_process_row = my_prow)
    !   IF (PRESENT (my_pcol)) &
    !        CALL get_blacs_info (blacs_env, my_process_column = my_pcol)
    !ENDIF
    IF (PRESENT (my_prow)) &
       my_prow = dbcsr_mp_myprow (dbcsr_distribution_mp (matrix%dist))
    IF (PRESENT (my_pcol)) &
       my_pcol = dbcsr_mp_mypcol (dbcsr_distribution_mp (matrix%dist))
  END SUBROUTINE dbcsr_get_info


! *****************************************************************************
!> \brief Creates a working matrix.
!> \param[out] matrix         new matrix
!> \param[in] nblks_guess     (optional) estimated number of blocks
!> \param[in] sizedata_guess  (optional) estimated size of data
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_work_create(matrix, nblks_guess, sizedata_guess, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN), OPTIONAL            :: nblks_guess, sizedata_guess
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_work_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nblks, nze, stat

!   ---------------------------------------------------------------------------

    matrix%w%lastblk = 0
    matrix%w%datasize = 0
    IF(PRESENT(nblks_guess)) THEN
       nblks = nblks_guess
    ELSE
       nblks = CEILING (REAL (matrix%nblkrows_local * matrix%nblkcols_local)&
            / REAL (dbcsr_mp_numnodes (dbcsr_distribution_mp (matrix%dist))))
    ENDIF
    IF(PRESENT(sizedata_guess)) THEN
       nze = sizedata_guess
    ELSE
       nze = CEILING(REAL(matrix%nfullrows_local * matrix%nfullcols_local)&
            / REAL(dbcsr_mp_numnodes (dbcsr_distribution_mp (&
            matrix%dist))))
    ENDIF
    NULLIFY(matrix%w%data)
    ALLOCATE(matrix%w%data(nze), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'matrix%w%data', nze)
    !
    ALLOCATE(matrix%w%row_i(nblks), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'matrix%w%row_i', nblks)
    !
    ALLOCATE(matrix%w%col_i(nblks), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'matrix%w%col_i', nblks)
    !
    ALLOCATE(matrix%w%blk_p(nblks), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'matrix%w%blk_p', nblks)
  END SUBROUTINE dbcsr_work_create


! *****************************************************************************
!> \brief Creates the final dbcsr_type matrix from the working matrix.
!> \par The final matrix is present as the element work_matrix%m.
!> \note The columns and rows must be ordered! We don't support appending
!>       data to an existing matrix.
!> \param[in,out] matrix      final matrix
!> \param[in] resort          whether the indices should be sorted, default
!>                            is true
!> \param[in] reshuffle       whether the data should be reshuffled,
!>                            default is false
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_finalize(matrix, resort, reshuffle, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: resort, reshuffle
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p_p, col_i_p, &
                                                prev_row, row, row_p_p, stat, &
                                                timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sizes
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: shfl, srt
    TYPE(dbcsr_work_type)                    :: wm

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    IF (dbg) WRITE(*,*)routineN//' finalizing matrix', matrix%name
    wm = matrix%w
    srt = .TRUE.
    IF (PRESENT (resort)) srt = resort
    IF (srt) THEN
       CALL dbcsr_sort_indices(wm%lastblk, wm%row_i, wm%col_i, wm%blk_p,&
            error=error)
       !CALL dbcsr_wm_merge_blocks(wm, error)
    ENDIF
    IF(.NOT.ASSOCIATED(matrix%data)) THEN
       IF (matrix%special_memory) THEN
          CALL mp_allocate(matrix%data, wm%datasize, stat=stat)
       ELSE
          ALLOCATE(matrix%data(wm%datasize), stat=stat)
       ENDIF
       IF (stat /= 0) CALL stop_memory(routineP, 'matrix%data', wm%datasize)
    ENDIF
    IF (SIZE(matrix%data) .LT. wm%datasize) THEN
       IF (matrix%special_memory) THEN
          CALL mp_deallocate(matrix%data, stat=stat)
          CALL mp_allocate(matrix%data, wm%datasize, stat=stat)
       ELSE
          DEALLOCATE(matrix%data, stat=stat)
          ALLOCATE(matrix%data(wm%datasize), stat=stat)
       ENDIF
       IF (stat /= 0) CALL stop_memory(routineP, 'm%data', wm%datasize)
    ENDIF
    IF (dbg) THEN
       !DO i = LBOUND(wm%data,1), UBOUND(wm%data,1)
       !   write(*,*)'array data, index',i,'=',wm%data(i)
       !ENDDO
       WRITE(*,*)routineP//' non-finalized DP(',LBOUND(wm%data,1),'-',&
            UBOUND(wm%data,1),') size=',wm%datasize,'=>',&
            DOT_PRODUCT(wm%data(1:wm%datasize), wm%data(1:wm%datasize))
    ENDIF
    IF (dbg) WRITE(*,*)routineP//' copyfinalized DP=',&
         DOT_PRODUCT(matrix%data(1:wm%datasize), matrix%data(1:wm%datasize))
    matrix%nblks = wm%lastblk
    matrix%nze = wm%datasize
    matrix%index(dbcsr_slot_nblks) = matrix%nblks
    matrix%index(dbcsr_slot_nze) = matrix%nze
    DBG routineP//' index',matrix%index(:dbcsr_num_slots)
    DBG routineP//' # blocks=',matrix%nblks
    !> We first figure out how large the indices should be.
    CALL dbcsr_addto_index_array(matrix, dbcsr_slot_row_p,&
         reservation=matrix%nblkrows_total+1, extra=wm%lastblk*2, error=error)
    matrix%row_p(1) = 0
    CALL dbcsr_addto_index_array(matrix, dbcsr_slot_col_i,&
         reservation=wm%lastblk, error=error)
    CALL dbcsr_addto_index_array(matrix, dbcsr_slot_blk_p,&
         reservation=wm%lastblk, error=error)
    matrix%row_p(1) = 0
    CALL dbcsr_repoint_index(matrix, error=error)
    !
    !> Now go through the indices and copy the data.
    row_p_p = 1
    col_i_p = 1
    blk_p_p = 1
    matrix%row_p(1) = 0
    prev_row = 1
    DO blk = 1, wm%lastblk
       !> If the current block's row is different (i.e., higher), then
       !> update the row_p array.
       row = wm%row_i(blk)
       IF (dbg.AND.bcsr_verbose) WRITE(*,*)routineP//' blk',blk,'row,col,blk',&
            wm%row_i(blk), wm%col_i(blk), wm%blk_p(blk)
       CALL cp_assert(row.GE.prev_row, cp_failure_level,&
            cp_unimplemented_error_nr, routineP, 'Unordered rows', error)
       IF (row .NE. prev_row) THEN
          matrix%row_p(prev_row+1:row) = blk-1
          prev_row = row
       ENDIF
    ENDDO !blk
    matrix%row_p(prev_row+1:matrix%nblkrows_total+1) = wm%lastblk
    IF(dbg.AND.bcsr_verbose)&
         WRITE(*,*)routineP//' row_p=',matrix%row_p,'size=',SIZE(matrix%row_p)
    IF(dbg.AND.bcsr_verbose)&
         WRITE(*,*)routineP//' col_i=',wm%col_i
    matrix%col_i(:) = wm%col_i(1:matrix%nblks)
    shfl = .FALSE.
    IF (PRESENT (reshuffle)) shfl = reshuffle
    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)
    IF (.NOT.shfl) THEN
       matrix%blk_p(:) = wm%blk_p(1:matrix%nblks)
       matrix%data(1:wm%datasize) = wm%data(1:wm%datasize)
    ELSE
       ALLOCATE (sizes(matrix%nblks), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineP, "sizes", matrix%nblks)
       CALL calc_blk_sizes (sizes, matrix%row_p, matrix%col_i,&
            array_data(matrix%row_blk_size), array_data(matrix%col_blk_size))
       CALL dbcsr_sort_data(matrix%blk_p, wm%blk_p(1:matrix%nblks),&
            sizes, matrix%data, wm%data, error)
       DEALLOCATE(sizes)
    ENDIF

    wm%datasize = 0
    wm%lastblk = 0
    matrix%w = wm
    CALL dbcsr_work_destroy(matrix%w, keepfinal=.TRUE., keepfinaldata=.TRUE.,&
         error=error)
    matrix%valid = .TRUE.
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_finalize


! *****************************************************************************
!> \brief Deallocates and destroys a matrix.
!> \param[in,out] matrix      matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_destroy(matrix, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_destroy', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF (matrix%special_memory) THEN
       CALL MP_DEALLOCATE(matrix%index)
    ELSE
       DEALLOCATE(matrix%index)
    ENDIF
    IF (matrix%special_memory) THEN
       CALL MP_DEALLOCATE(matrix%data)
    ELSE
       DEALLOCATE(matrix%data)
    ENDIF
    CALL array_release (matrix%row_blk_size)
    CALL array_release (matrix%col_blk_size)
    CALL dbcsr_distribution_release(matrix%dist)
    matrix%valid = .FALSE.
  END SUBROUTINE dbcsr_destroy


! *****************************************************************************
!> \brief Deallocates and destroys a work matrix.
!> \param[in,out] wm          work matrix
!> \param[in,out] error       cp2k error
!> \param keepdata    do not deallocate data
!> \param keepfinal   do not destroy the final, non-work matrix
!> \param keepfinaldata       do not destroy the data in the final,
!>                            non-work matrix
! *****************************************************************************
  SUBROUTINE dbcsr_work_destroy(wm, error, keepdata, keepfinal,&
       keepfinaldata)
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: wm
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: keepdata, keepfinal, &
                                                keepfinaldata

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_work_destroy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: datasharing, dbg, lkeepdata, &
                                                lkeepfinal, lkeepfinaldata

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    lkeepdata = .FALSE.
    IF (PRESENT (keepdata)) lkeepdata = keepdata
    lkeepfinal = .TRUE.
    IF (PRESENT (keepfinal)) lkeepfinal = keepfinal
    lkeepfinaldata = .FALSE.
    IF (PRESENT (keepfinaldata)) lkeepfinaldata = keepfinaldata
    datasharing = .FALSE.
    !datasharing = ASSOCIATED(wm%data, wm%m%data) .AND. lkeepfinaldata
    IF (dbg) WRITE(*,*)routineP//' keepdata, keepfinal, keepfinaldata, share',&
         lkeepdata, lkeepfinal, lkeepfinaldata, datasharing
    !IF (.NOT. lkeepfinal .AND. ASSOCIATED (wm%m)) &
    !     CALL dbcsr_destroy(wm%m, error, keepdata=lkeepfinaldata)
    IF (ASSOCIATED (wm%data).AND..NOT.datasharing.AND..NOT.lkeepdata) THEN
       DEALLOCATE(wm%data)
       NULLIFY (wm%data)
    ENDIF
    IF (ASSOCIATED (wm%row_i)) THEN
       DEALLOCATE(wm%row_i)
       NULLIFY (wm%row_i)
    ENDIF
    IF (ASSOCIATED (wm%col_i)) THEN
       DEALLOCATE(wm%col_i)
       NULLIFY (wm%row_i)
    ENDIF
    IF (ASSOCIATED (wm%blk_p)) THEN
       DEALLOCATE(wm%blk_p)
       NULLIFY (wm%blk_p)
    ENDIF
  END SUBROUTINE dbcsr_work_destroy


! *****************************************************************************
!> \brief Converts a CP full matrix into a BCSR matrix.
!> \param[out] matrix         the created BCSR matrix
!> \param[in] fm    the cpk full matrix
!> \param[in] threshold       the threshold for determining sparsity
!> \param[in] row_blk_begin,row_blk_end Starting and ending full rows for
!>                                      blocked rows
!> \param[in] col_blk_begin,col_blk_end Starting and ending full columns for
!>                                      blocked columns
!> \param[in] dist2d          the distribution_2d to use for the new matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_from_fm(matrix, fm, threshold, dist2d, row_blk_size,&
       col_blk_size, error)
    TYPE(dbcsr_type), INTENT(OUT)            :: matrix
    TYPE(cp_fm_type), POINTER                :: fm
    REAL(KIND=dp), INTENT(IN)                :: threshold
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: dist2d
    INTEGER, DIMENSION(:), INTENT(IN)        :: row_blk_size, col_blk_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_from_fm', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=80)                        :: fm_name
    INTEGER :: blkcol, blkcol_l, blkrow, blkrow_l, fullrow, i, j, m, mpicomm, &
      mynod, mypcol, myprow, n, nblkcols_local, nblkcols_total, &
      nblkrows_local, nblkrows_total, nfullcols_local, nfullcols_total, &
      nfullrows_local, nfullrows_total, npcols, nprows, numnod, prev_blkrow, &
      prow, stat, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blk_ps, col_blk_begin, &
                                                row_blk_begin
    INTEGER, DIMENSION(:), POINTER           :: blkcols_local, blkrows_local, &
                                                fullcols_local, &
                                                fullrows_local, pcol_dist, &
                                                prow_dist
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi, mpi2blacs
    LOGICAL                                  :: dbg, failure = .FALSE.
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: thefmrow
    TYPE(array_i1d_obj)                      :: col_blk_size_array, &
                                                row_blk_size_array
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dbcsr_work_type)                    :: wm

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    CALL timeset(routineN, timing_handle)
    CALL cp_fm_get_info(fm, name=fm_name, matrix_struct=fm_struct,&
         para_env=para_env, context=blacs_env,&
         row_indices=fullrows_local, col_indices=fullcols_local,&
         nrow_global=nfullrows_total, ncol_global=nfullcols_total,&
         error=error)
    CALL distribution_2d_get(dist2d, row_distribution=prow_dist,&
         col_distribution=pcol_dist, error=error, blacs_env=blacs_env,&
         n_flat_local_rows=nblkrows_local, n_flat_local_cols=nblkcols_local,&
         n_row_distribution=nblkrows_total, n_col_distribution=nblkcols_total,&
         flat_local_rows=blkrows_local, flat_local_cols=blkcols_local)
    CALL get_blacs_info(blacs_env, para_env=para_env, blacs2mpi=blacs2mpi,&
         mpi2blacs=mpi2blacs, my_process_row=myprow, my_process_column=mypcol,&
         number_of_process_rows=nprows, number_of_process_columns=npcols)
    mynod = blacs_env%para_env%mepos
    numnod = blacs_env%para_env%num_pe
    mpicomm = blacs_env%para_env%group
    nfullrows_local = 0
    nfullcols_local = 0
    IF (dbg) WRITE(*,*)routineP//' myprow, mypcol',myprow,'/',nprows,';',&
         mypcol,'/',npcols
    IF (dbg) CALL cp_blacs_env_write(blacs_env, 6, error)
    IF (dbg) WRITE(*,*)routineP//' row_blk_size=',row_blk_size
    IF (dbg) WRITE(*,*)routineP//' blkrows_local=',blkrows_local
    IF (dbg) WRITE(*,*)routineP//' blkcols_local=',blkcols_local
    ALLOCATE(row_blk_begin(nblkrows_total))
    row_blk_begin(1) = 1
    DO blkrow = 2, nblkrows_total
       row_blk_begin(blkrow) = row_blk_begin(blkrow-1)+row_blk_size(blkrow-1)
    ENDDO
    ALLOCATE(col_blk_begin(nblkcols_total))
    col_blk_begin(1) = 1
    DO blkcol = 2, nblkcols_total
       col_blk_begin(blkcol) = col_blk_begin(blkcol-1)+col_blk_size(blkcol-1)
    ENDDO
    DO blkrow_l = 1, nblkrows_local
       blkrow = blkrows_local(blkrow_l)
       nfullrows_local = nfullrows_local +&
            row_blk_size(blkrow)
    ENDDO
    DO blkcol_l = 1, nblkcols_local
       blkcol = blkcols_local(blkcol_l)
       nfullcols_local = nfullcols_local +&
            col_blk_size(blkcol)
    ENDDO
    ! Create local work matrix.
    CALL array_new (row_blk_size_array, row_blk_size, lb=1)
    CALL array_new (col_blk_size_array, col_blk_size, lb=1)
    CALL dbcsr_distribution_init (dist)
    CALL dbcsr_dist2d_to_dist (dist2d, dist, error)
    CALL dbcsr_init (matrix, error)
    CALL dbcsr_create(matrix, fm%name, dist, 'N',&
         row_blk_size_array, col_blk_size_array,&
         0, 0, error=error)
    CALL array_release (row_blk_size_array)
    CALL array_release (col_blk_size_array)
    CALL dbcsr_distribution_release (dist)
    ! Somehow assign row & column block sizes.
    CALL dbcsr_work_create(matrix, error=error)
    wm = matrix%w
    ALLOCATE(blk_ps(nblkcols_local), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'blk_ps', nblkcols_local+1)
    ! Figure out the maximum row & column block sizes
    m = MAXVAL(row_blk_size)
    n = MAXVAL(col_blk_size)
    IF(dbg) WRITE(*,*)routineP//' max size',m,'x',n
    ALLOCATE(thefmrow(m,nfullcols_total), stat=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    !> This is what happens here: I go through the rows and obtain the FM block.
    !> How do I update the indexes of the wm? I just add to the array in
    !> the wm matrix and add some data to the data array.
    prev_blkrow = 0
    fullrow = 1
    wm%lastblk = 0
    wm%datasize = 0
    wm%data(:) = 999
    IF(dbg) WRITE(*,*)routineP//' my local blocks',blkrows_local
    DO blkrow = 1, nblkrows_total
       prow = prow_dist(blkrow)
       ! This processor must participate in forming the thefmrow, but I don't
       ! need it until it's one of my local rows.
       IF(dbg) WRITE(*,*)routineP//' blocked row', blkrow
       m = row_blk_size(blkrow)
       IF(dbg) WRITE(*,*)routineP//' Submatrix',row_blk_size(blkrow),' x 1:',&
            nfullcols_total,'size',m,'x',&
            nfullcols_total
       DEALLOCATE(thefmrow)
       ALLOCATE(thefmrow(m,nfullcols_total), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP, 'thefmrow', m*n)
       thefmrow(:,:) = 888
       CALL cp_fm_get_submatrix(fm, thefmrow,&
            start_row=row_blk_begin(blkrow), n_rows=m,&
            start_col=1, n_cols=nfullcols_total, error=error)
       !IF (dbg) WRITE(*,*)thefmrow
       IF(dbg) WRITE(*,*)routineP//' Chksum=',&
            DOT_PRODUCT(RESHAPE(thefmrow,(/m*nfullcols_total/)),&
            &           RESHAPE(thefmrow,(/m*nfullcols_total/)))
       IF (myprow .EQ. prow) THEN
          ! The array-adding should be done in one pass
          IF(dbg) WRITE(*,*)routineP//' adding row',blkrow,'sized',m,&
               'for local columns', blkcols_local,'fullrows',fullrow
          CALL addto_array(wm%row_i, (/ (blkrow, i=1,nblkcols_local) /),&
               wm%lastblk+1, error=error)
          CALL addto_array(wm%col_i, blkcols_local, wm%lastblk+1, error=error)
          DO blkcol_l = 1, nblkcols_local
             blkcol = blkcols_local(blkcol_l)
             i = col_blk_begin(blkcol)
             j = i + col_blk_size(blkcol)-1
             n = col_blk_size(blkcol)
             blk_ps(blkcol_l) = wm%datasize + 1
             IF(dbg) WRITE(*,*)routineP//' adding block at column',blkcol,&
                  'sized',n,'(',i,':',j,')'
             IF(dbg) WRITE(*,*)routineP//' with chksum',&
                  DOT_PRODUCT(RESHAPE(thefmrow(1:m,i:j),(/m*n/)),&
                  &           RESHAPE(thefmrow(1:m,i:j),(/m*n/)))
             CALL addto_array(wm%data, RESHAPE(thefmrow(1:m,i:j),(/m*n/)),&
                  offset=blk_ps(blkcol_l), error=error)
             wm%datasize = wm%datasize + m*n
          ENDDO
          IF(dbg) WRITE(*,*)routineP//' with blk_ps',blk_ps
          CALL addto_array(wm%blk_p, blk_ps(1:nblkcols_local), wm%lastblk+1,&
               error=error)
          wm%lastblk = wm%lastblk + nblkcols_local
       ENDIF ! myprow = prow
    ENDDO !blkrow
    DEALLOCATE(blk_ps, thefmrow)
    ! Construct the final matrix.
    matrix%w = wm
    CALL dbcsr_finalize(matrix, error=error)
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_from_fm


! *****************************************************************************
!> \brief Create a DBCSR matrix from a real_matrix_type
!> \param[out] matrix         bcsr matrix
!> \param[in] sm    cp2k sparse matrix
!> \param[in,out] error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_from_sm(matrix, sm, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    TYPE(real_matrix_type), POINTER          :: sm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_from_sm', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: data_exchange_info_cnt = 3, &
                                                metalen = 2

    INTEGER :: absrow, blk, blk_l, blk_p, block_col, col, dst_p, mp_group, &
      mynod, mypcol, myprow, nblkcols_local, nblkcols_total, nblkrows_local, &
      nblkrows_total, nblks, newpcol, newprow, nfullcols_local, &
      nfullcols_total, nfullrows_local, nfullrows_total, npcols, nprows, &
      numproc, nze, p, row, row_l, src_p, stat, stored_col, stored_row, &
      timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                sd_disp, sdp, send_meta, &
                                                sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count
    INTEGER, DIMENSION(:), POINTER :: blkcols_local, blkrows_local, &
      col_blk_begin, col_blk_end, col_blk_size, pcol_dist, prow_dist, &
      row_blk_begin, row_blk_end, row_blk_size
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: symmetric, tr
    REAL(KIND=dp)                            :: tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: send_data
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: recv_data
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blkdata_2d
    TYPE(array_i1d_obj)                      :: cbs_array, rbs_array
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(distribution_2d_type), POINTER      :: dist2d
    TYPE(real_block_node_type), POINTER      :: sm_block_node

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    tstart = m_walltime ()
    dist2d => sm%distribution_2d
    CALL dbcsr_distribution_init (dist)
    CALL dbcsr_dist2d_to_dist (dist2d, dist, error)
    symmetric = (sm%symmetry .EQ. "symmetric")
    CALL cp_assert (symmetric, cp_warning_level, cp_caller_error,&
         routineN, "Matrix not symmetric.", error=error)
    matrix%matrix_type = 'S'
    matrix%dist2d => dist2d
    nblkrows_total = sm%nblock_row
    nblkcols_total = sm%nblock_col
    row_blk_begin => sm%first_row
    row_blk_end => sm%last_row
    col_blk_begin => sm%first_col
    col_blk_end => sm%last_col
    CALL distribution_2d_get(dist2d, row_distribution=prow_dist,&
         col_distribution=pcol_dist, error=error, blacs_env=blacs_env,&
         n_flat_local_rows=nblkrows_local, n_flat_local_cols=nblkcols_local,&
         n_row_distribution=nblkrows_total, n_col_distribution=nblkcols_total,&
         flat_local_rows=blkrows_local, flat_local_cols=blkcols_local)
    nfullrows_total = MAXVAL(row_blk_end)
    nfullcols_total = MAXVAL(col_blk_end)
    nfullrows_local = 0
    nfullcols_local = 0
    DO row_l = 1, nblkrows_local
       absrow = blkrows_local(row_l)
       nfullrows_local = nfullrows_local +&
            row_blk_end(absrow) - row_blk_begin(absrow) + 1
    ENDDO
    DO row_l = 1, nblkcols_local
       absrow = blkcols_local(row_l)
       nfullcols_local = nfullcols_local +&
            col_blk_end(absrow) - col_blk_begin(absrow) + 1
    ENDDO
    CALL get_blacs_info(blacs_env, para_env=para_env, blacs2mpi=blacs2mpi,&
         my_process_row=myprow, my_process_column=mypcol,&
         number_of_process_rows=nprows, number_of_process_columns=npcols)
    IF (dbg) CALL cp_blacs_env_write(blacs_env, 6, error)
    IF (dbg) WRITE(*,*)routineP//' myprow, mypcol',myprow,'/',nprows,';',&
         mypcol,'/',npcols
    mynod = blacs_env%para_env%mepos
    numproc = blacs_env%para_env%num_pe
    mp_group = para_env%group
    CALL cp_assert(SIZE(prow_dist) .EQ. SIZE(pcol_dist), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Unequal row and column distributions for symmetric matrix.', error)
    CALL cp_assert(MAXVAL(prow_dist).LE.UBOUND(blacs2mpi,1), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows', error)
    CALL cp_assert(MAXVAL(prow_dist).EQ.UBOUND(blacs2mpi,1), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of row distribution not equal to processor rows', error)
    CALL cp_assert(MAXVAL(pcol_dist).LE.UBOUND(blacs2mpi,2), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols', error)
    CALL cp_assert(MAXVAL(pcol_dist).EQ.UBOUND(blacs2mpi,2), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of col distribution not equal to processor cols', error)
    ALLOCATE (row_blk_size (nblkrows_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, 'row_blk_size', nblkrows_total)
    ALLOCATE (col_blk_size (nblkcols_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, 'col_blk_size', nblkcols_total)
    row_blk_size = row_blk_end(:)-row_blk_begin(:)+1
    col_blk_size = col_blk_end(:)-col_blk_begin(:)+1
    ALLOCATE (send_count(2, 0:numproc-1))
    ALLOCATE (recv_count(2, 0:numproc-1))
    ALLOCATE (sdp(0:numproc-1))
    ALLOCATE (sd_disp(0:numproc-1))
    ALLOCATE (smp(0:numproc-1))
    ALLOCATE (sm_disp(0:numproc-1))
    ALLOCATE (rd_disp(0:numproc-1))
    ALLOCATE (rm_disp(0:numproc-1))
    ! Count what to send to others.
    send_count(:,:) = 0
    DO absrow = 1, nblkrows_total
       newprow = prow_dist(absrow)
       sm_block_node => first_block_node(sm, absrow)
       DO WHILE (ASSOCIATED(sm_block_node))
          CALL get_block_node(block_node=sm_block_node, block_col=block_col,&
               block=blkdata_2d)
          stored_row = absrow
          stored_col = block_col
          !tr = checker_tr (stored_row, stored_col)
          !IF (tr) CALL swap (stored_row, stored_col)
          CALL dbcsr_get_lookup(matrix, stored_row, stored_col, tr)
          newprow = prow_dist(stored_row)
          newpcol = pcol_dist(stored_col)
          dst_p = blacs2mpi(newprow, newpcol)
          nze = row_blk_size(absrow) * col_blk_size(block_col)
          send_count(1, dst_p) = send_count(1, dst_p) + 1
          send_count(2, dst_p) = send_count(2, dst_p) + nze
          IF (dbg) WRITE(*,*)routineP//' To send block from',mynod,'to',&
               dst_p,'block at (',absrow,block_col,') with',nze,&
               'elements, position','?','as (',stored_row,',',stored_col,')',tr
          sm_block_node => next_block_node(sm_block_node)
       END DO ! associated
    ENDDO ! absrow
    DBG 'B'
    ! Exchange data on what to send
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
    ! Allocate data structures needed for data exchange.
    ALLOCATE (recv_data(SUM(recv_count(2, :))))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :))))
    ALLOCATE (send_data(SUM(send_count(2, :))))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :))))
    DBG 'send data len',SIZE (send_data)
    DBG 'send meta len',SIZE (send_meta)
    DBG 'recv data len',SIZE (recv_data)
    DBG 'recv meta len',SIZE (recv_meta)
    ! Fill in the meta data structures and copy the data.
    sd_disp = -1 ; sm_disp = -1
    rd_disp = -1 ; rm_disp = -1
    sd_disp(0) = 1 ; sm_disp(0) = 1
    rd_disp(0) = 1 ; rm_disp(0) = 1
    DO dst_p = 1, numproc-1
       sm_disp(dst_p) = sm_disp(dst_p-1) + metalen*send_count(1, dst_p-1)
       sd_disp(dst_p) = sd_disp(dst_p-1) + send_count(2, dst_p-1)
       rm_disp(dst_p) = rm_disp(dst_p-1) + metalen*recv_count(1, dst_p-1)
       rd_disp(dst_p) = rd_disp(dst_p-1) + recv_count(2, dst_p-1)
    ENDDO
    sdp = sd_disp ; smp = sm_disp
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    DBG 'D'
    !> Now exchange information on the upcoming transfers.
    ! Perform exchanges
    DO absrow = 1, nblkrows_total
       newprow = prow_dist(absrow)
       sm_block_node => first_block_node(sm, absrow)
       DO WHILE (ASSOCIATED(sm_block_node))
          CALL get_block_node(block_node=sm_block_node, block_col=block_col,&
               block=blkdata_2d)
          stored_row = absrow
          stored_col = block_col
          !tr = checker_tr (stored_row, stored_col)
          !IF (tr) CALL swap (stored_row, stored_col)
          CALL dbcsr_get_lookup(matrix, stored_row, stored_col, tr)
          newprow = prow_dist(stored_row)
          newpcol = pcol_dist(stored_col)
          dst_p = blacs2mpi(newprow, newpcol)
          nze = row_blk_size(stored_row) * col_blk_size(stored_col)
          send_meta(smp(dst_p)) = stored_row
          send_meta(smp(dst_p)+1) = stored_col
          IF (.NOT. tr) THEN
             send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                  RESHAPE (blkdata_2d, (/nze/))
                  !sm%data(blk_p:blk_p+nze-1)
          ELSE
             !send_meta(smp(dst_p)) = -stored_row
             send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                  RESHAPE (TRANSPOSE(blkdata_2d), (/nze/))
                  !sm%data(blk_p:blk_p+nze-1)
          ENDIF
          smp(dst_p) = smp(dst_p) + metalen
          sdp(dst_p) = sdp(dst_p) + nze
          !
          sm_block_node => next_block_node(sm_block_node)
       END DO ! associated
    ENDDO ! absrow
    !
    ! Exchange the data and metadata structures.
    DBG 'E'
    CALL mp_alltoall(send_data(:), send_count(2,:), sd_disp(:)-1,&
         recv_data(:), recv_count(2,:), rd_disp(:)-1, mp_group)
    DBG 'F'
    DBG 'send_meta',send_meta
    CALL mp_alltoall(send_meta(:), metalen*send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*recv_count(1,:), rm_disp(:)-1, mp_group)
    DBG 'recv_meta',recv_meta
    DBG 'G'
    ! Now fill in the data.
    blk = SUM(recv_count(1,:))
    nze = SUM(recv_count(2,:))
    CALL array_new (rbs_array, row_blk_size, gift=.TRUE.)
    CALL array_new (cbs_array, col_blk_size, gift=.TRUE.)
    CALL dbcsr_init (matrix, error)
    CALL dbcsr_create(matrix, sm%name, dist, sm%symmetry,&
         rbs_array, cbs_array, nblks=blk, nze=nze, error=error)
    CALL array_release (rbs_array)
    CALL array_release (cbs_array)
    CALL dbcsr_work_create(matrix,&
            nblks_guess=blk,&
            sizedata_guess=nze, error=error)
    !DBG 'send disps data',sd_disp
    !DBG 'send disps meta',sm_disp
    !DBG 'recv disps data',rd_disp
    !DBG 'recv disps meta',rm_disp
    blk_p = 1
    blk = 1
    DO src_p = 0, numproc-1
       nze = recv_count(2, src_p)
       matrix%w%data(blk_p:blk_p+nze-1) =&
            recv_data(rd_disp(src_p):rd_disp(src_p)+nze-1)
       DBG 'processing from',src_p,'counts',send_count(1, src_p)
       DO blk_l = 1, recv_count(1, src_p)
          stored_row = recv_meta(rm_disp(src_p)+metalen*(blk_l-1))
          stored_col = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+1)
          !DBG 'blk: p',blk,src_p,'; row, col',stored_row,stored_col,&
          !     '; blk_p',blk_p,'; row_rep, rptdist',row_rep,&
          !     row_rep_dist (ABS (stored_row))
          DBG 'blk: p',blk,src_p,'; row, col',&
               stored_row,stored_col,&
               '; blk_p',blk_p
          matrix%w%row_i(blk) = stored_row ! It's always positive here
          !matrix%w%row_i(blk) = ABS(stored_row)
          matrix%w%col_i(blk) = stored_col
          matrix%w%blk_p(blk) = SIGN(blk_p, stored_row)
          !matrix%w%blk_p(blk) = SIGN(blk_p, stored_row)
          nze = row_blk_size(ABS(stored_row))&
               * col_blk_size(stored_col)
          blk_p = blk_p + nze
          blk = blk + 1
       ENDDO ! blk_l
    ENDDO ! src_p
    matrix%w%lastblk = blk - 1
    matrix%w%datasize = blk_p - 1
    !
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_data)
    DEALLOCATE(send_meta)
    CALL dbcsr_finalize (matrix, resort=.TRUE., reshuffle=.TRUE.,&
         error=error)
    CALL dbcsr_distribution_release (dist)
    tstop = m_walltime ()
    DBG 'time:',tstop-tstart
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_from_sm


! *****************************************************************************
!> \brief Transposes a DBCSR matrix.
!> \param[out] transposed     transposed DBCSR matrix
!> \param[in] normal          input DBCSR matrix
!> \param[in] reuse           (optional) reuse an existing matrix for the
!>                            transpose, default is False
!> \param[in] reshuffle       (optional) reshuffle data, default is True
!> \param[in] transpose_data  (optional) transpose data blocks, default is True
!> \param[in] transpose_distribution    (optional) transpose the distribution
!>                                      from the input matrix, default is True
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_new_transposed (transposed, normal, reshuffle,&
       transpose_data, transpose_distribution, reuse, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: transposed
    TYPE(dbcsr_type), INTENT(IN)             :: normal
    LOGICAL, INTENT(IN), OPTIONAL            :: reshuffle, transpose_data, &
                                                transpose_distribution, reuse
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_new_transposed', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: new_type
    INTEGER                                  :: stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blk_p
    LOGICAL                                  :: hijack, shuffle, tr_blocks, &
                                                tr_dist
    TYPE(dbcsr_distribution_obj)             :: new_dist

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    ! Internalize options
    hijack = .FALSE.
    IF (PRESENT (reuse)) hijack = reuse
    shuffle = .TRUE.
    IF (PRESENT (reshuffle)) shuffle = reshuffle
    tr_blocks = .TRUE.
    IF (PRESENT (transpose_data)) tr_blocks = transpose_data
    tr_dist = .TRUE.
    IF (PRESENT (transpose_distribution)) tr_dist = transpose_distribution
    ! Prepare the distribution for the transposed matrix
    CALL cp_assert (tr_dist, cp_failure_level, cp_unimplemented_error_nr,&
         routineP, "Distribution must be transposed.", error)
    !CALL cp_assert (tr_dist.OR.(normal%nblkrows_total.EQ.nblkcols_total),&
    !     cp_failure_level, cp_caller_error, routineP,&
    !     "Distribution of non-square matrix should be transposed", error)
    IF (tr_dist) THEN
       !write(*,*)' transposing dist...'
       CALL dbcsr_transpose_distribution (new_dist, normal%dist)
    ELSE
       new_dist = normal%dist
    ENDIF
    ! Create the transposed matrix
    IF (.NOT. hijack) CALL dbcsr_init (transposed, error)
    !CALL dbcsr_transposed_type (new_type, normal%matrix_type)
    new_type = normal%matrix_type
    CALL dbcsr_create (transposed, "transposed "//normal%name, new_dist,&
         new_type, normal%col_blk_size, normal%row_blk_size, normal%nblks,&
         normal%nze, normal%special_memory, reuse=hijack, error=error)
    ! Reserve the space for the new indices.
    CALL dbcsr_addto_index_array (transposed, dbcsr_slot_row_p,&
         reservation=transposed%nblkrows_total+1, extra=transposed%nblks*2,&
         error=error)
    CALL dbcsr_addto_index_array (transposed, dbcsr_slot_col_i,&
         reservation=normal%nblks, error=error)
    CALL dbcsr_addto_index_array (transposed, dbcsr_slot_blk_p,&
         reservation=normal%nblks, error=error)
    CALL dbcsr_repoint_index (transposed, error=error)
    CALL dbcsr_ensure_data_size (transposed, normal%nze, nocopy=.TRUE.,&
         error=error)
    !
    transposed%nblks = normal%nblks
    transposed%nze = normal%nze
    transposed%index(dbcsr_slot_nblks) = normal%nblks
    transposed%index(dbcsr_slot_nze) = normal%nze
    ! Transpose the local index.
    DBG 'meta', transposed%nblkrows_total+1, normal%nblks
    DBG 'size', SIZE (transposed%row_p), SIZE (transposed%col_i)
    ALLOCATE (blk_p (normal%nblks), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, "blk_p", normal%nblks)
    CALL transpose_index_local (transposed%row_p, transposed%col_i,&
         normal%row_p, normal%col_i, blk_p, normal%blk_p)
    ! Copy the data
    CALL dbcsr_copy_sort_data (transposed%blk_p, blk_p, transposed%row_p,&
         transposed%col_i, array_data (transposed%row_blk_size),&
         array_data (transposed%col_blk_size),&
         transposed%data, normal%data, transpose_blocks=.TRUE., error=error)
    !
1315 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
    IF (dbg) THEN
       WRITE(*,*)'new FINAL index'
       WRITE(*,1315)transposed%row_p
       WRITE(*,1315)transposed%col_i
       WRITE(*,1315)transposed%blk_p
    ENDIF
    !
    DEALLOCATE (blk_p)
    IF (tr_dist) CALL dbcsr_distribution_release (new_dist)
  END SUBROUTINE dbcsr_new_transposed


! *****************************************************************************
!> \brief Determines the correct transposed type of a DBCSR matrix.
!> \param[out] new_type       new matrix type
!> \param[in] old_type        current matrix_type
! *****************************************************************************
  ELEMENTAL SUBROUTINE dbcsr_transposed_type (new_type, old_type)
    CHARACTER, INTENT(OUT)                   :: new_type
    CHARACTER, INTENT(IN)                    :: old_type

    SELECT CASE (old_type)
    CASE ('N')
       new_type = 'T'
    CASE ('T')
       new_type = 'N'
    CASE DEFAULT
       new_type = old_type
    END SELECT
  END SUBROUTINE dbcsr_transposed_type

! *****************************************************************************
!> \brief Re-indexes row_p and blk_i according to columns.
!> \par The re-indexing is equivalent to a local-only transpose.
!> \param[out] new_col_p      new column pointer
!> \param[out] new_row_i      new row index
!> \param[in] old_row_p       old row pointer
!> \param[in] old_col_i       old column index
!> \param[out] new_blk_p      (optional) new block pointer
!> \param[in] old_blk_p       (optional) old block pointer
! *****************************************************************************
  SUBROUTINE transpose_index_local (new_col_p, new_row_i, old_row_p,&
       old_col_i, new_blk_p, old_blk_p)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: new_col_p, new_row_i
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_row_p, old_col_i
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: new_blk_p
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: old_blk_p

    CHARACTER(len=*), PARAMETER :: routineN = 'transpose_index_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: curcol, mincol, nblks, &
                                                new_ncols, old_max_col, &
                                                old_nrows, old_row, prev_col
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: row_p
    LOGICAL                                  :: column_le

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY

1323 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
    new_ncols = SIZE (new_col_p)-1
    old_nrows = SIZE (old_row_p)-1
    ALLOCATE (row_p (SIZE (old_row_p)))
    row_p(:) = old_row_p(:)+1
    IF (SIZE (old_col_i) .GT. 0) THEN
       old_max_col = MAXVAL (old_col_i)
    ELSE
       old_max_col = 0
    ENDIF
    mincol = old_max_col
    !new_col_p(:) = 1
    IF (dbg) THEN
       WRITE(*,*)'old index'
       WRITE(*,1323)old_row_p
       WRITE(*,1323)old_col_i
       IF (PRESENT (old_blk_p)) THEN
          WRITE(*,*)'old blkp'
          WRITE(*,1323)old_blk_p
       ENDIF
    ENDIF
    ! Scan through the columns of the given matrix.
    ! row_p is an index, for each row, of the current column scanning position
    ! appropriate to the column being scanned.
    ! mincol is use to jump ahead for very sparse matrices
    nblks = 0
    new_col_p(1) = 0
    new_col_p(:) = 0
    prev_col = 0
    curcol = 1
    DO WHILE (curcol .LE. new_ncols)
       DBG 'curcol', curcol
       IF (dbg) THEN
          WRITE(*,*)'rowp'
          WRITE(*,1323)row_p
       ENDIF
       DO old_row = 1, old_nrows
          ! Condition 1:
          ! row_p(old_row) .LE. old_row_p(old_row+1)
          !           => there are still blocks in the row
          ! Condition 2:
          ! old_col_i(row_p(old_row)) .LE. curcol
          !           => the examined column in this row is <= to the
          !              curcol current column
          column_le = row_p(old_row) .LE. old_row_p(old_row+1)
          IF (column_le) column_le = old_col_i(row_p(old_row)) .LE. curcol
          ! In every row, try to find the column appropriate to the current
          ! column.
          DO WHILE (column_le)
             IF (old_col_i(row_p(old_row)) .EQ. curcol) THEN
                nblks = nblks + 1
                DBG 'Adding block',nblks,'at old row/new col',old_row
                !new_col_p(curcol+1) = i = new_col_p(curcol+1) + 1
                new_row_i(nblks) = old_row
                IF (PRESENT (new_blk_p) .AND. PRESENT (old_blk_p)) THEN
                   new_blk_p (nblks) = old_blk_p(row_p(old_row))
                ENDIF
             ENDIF
             row_p(old_row) = row_p(old_row) + 1
             ! We've come to the end of this row
             column_le = row_p(old_row) .LE. old_row_p(old_row+1)
             IF (column_le) column_le = old_col_i(row_p(old_row)) .LE. curcol
          ENDDO
          IF (row_p(old_row) .LE. old_row_p(old_row+1)) THEN
             mincol = MIN (mincol, old_col_i(row_p(old_row)))
             !DBG 'for row',old_row,'leftoff column is',old_col_i(row_p(old_row))
          ELSE
             !DBG 'for row',old_row,'forced column is ',old_col_i(MIN(row_p(old_row), old_max_col))
          ENDIF
       ENDDO
       !! Fill in missing indices.
       DBG 'prev_col, curcol, mincol, nblks', prev_col, curcol, mincol, nblks
       IF (prev_col+1 .LE. mincol) THEN
          DBG 'Skipped!'
       ENDIF
       !new_col_p(prev_col+2:mincol+1) = nblks
       new_col_p(prev_col+2:MIN(mincol+1,new_ncols+1)) = nblks
       prev_col = curcol
       curcol = MAX (mincol, curcol+1)
       !curcol = mincol
       mincol = mincol+1
       IF (.NOT. (prev_col .LT. curcol)) THEN
          DBG 'BAD new state: prev_col, curcol, mincol',&
               prev_col, curcol, mincol
          STOP
       ENDIF
       !IF(curcol+2 .LE. mincol .AND. curcol+2 .LE. old_max_col+1) THEN
       !   new_col_p(curcol+2:mincol) = new_col_p(curcol+1)
       !END IF
       !curcol = MAX(mincol,curcol+1)
       !mincol = old_max_col
    ENDDO
    IF (dbg) THEN
       WRITE(*,*)'new index'
       WRITE(*,1323)new_col_p
       WRITE(*,1323)new_row_i
       IF (PRESENT (new_blk_p)) THEN
          WRITE(*,*)'new blkp'
          WRITE(*,1323)new_blk_p
       ENDIF
    ENDIF
  END SUBROUTINE transpose_index_local


! *****************************************************************************
!> \brief Fills meta information from a given distribution_2d
!> \param[out] meta           meta information array to fill
!> \param dist                processor distribution
!> \param[in] row_blk_size    row block sizes
!> \param[in] col_blk_size    column block sizes
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE meta_from_dist (meta, dist, row_blk_size, col_blk_size, error)
    INTEGER, DIMENSION(dbcsr_meta_size), &
      INTENT(OUT)                            :: meta
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_from_dist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nfullcols_local, &
                                                nfullcols_total, &
                                                nfullrows_local, &
                                                nfullrows_total
    INTEGER, DIMENSION(:), POINTER           :: blkcols_local, blkrows_local, &
                                                cbs, rbs

!INTEGER                             :: row_l, absrow, stat

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    !CALL distribution_2d_get(dist2d, error=error,&
    !     n_flat_local_rows=nblkrows_local, n_flat_local_cols=nblkcols_local,&
    !     n_row_distribution=nblkrows_total, n_col_distribution=nblkcols_total,&
    !     flat_local_rows=blkrows_local, flat_local_cols=blkcols_local)
    blkrows_local => array_data (dbcsr_distribution_local_rows (dist))
    blkcols_local => array_data (dbcsr_distribution_local_cols (dist))
    rbs => array_data (row_blk_size)
    cbs => array_data (col_blk_size)
    nfullrows_total = SUM (rbs)
    nfullcols_total = SUM (cbs)
    nfullrows_local = 0
    nfullcols_local = 0
    DO i = 1, dbcsr_distribution_nlocal_rows (dist)
       nfullrows_local = nfullrows_local + rbs(blkrows_local(i))
    ENDDO
    DO i = 1, dbcsr_distribution_nlocal_cols (dist)
       nfullcols_local = nfullcols_local + cbs(blkcols_local(i))
    ENDDO
    meta(:) = 0
    meta(3)  = dbcsr_distribution_nrows (dist)
    meta(4)  = dbcsr_distribution_ncols (dist)
    meta(5)  = nfullrows_total
    meta(6)  = nfullcols_total
    meta(7)  = dbcsr_distribution_nlocal_rows (dist)
    meta(8)  = dbcsr_distribution_nlocal_cols (dist)
    meta(9)  = nfullrows_local
    meta(10) = nfullcols_local
  END SUBROUTINE meta_from_dist

  SUBROUTINE test_create_bl_distribution(error)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER, PARAMETER                       :: rpt = 10

    INTEGER                                  :: bins, elements, i
    REAL                                     :: rn
    TYPE(array_i1d_obj)                      :: block_distribution, block_size

    DO i = 1, rpt
       CALL RANDOM_NUMBER (rn)
       elements = INT(rn * 500)
       CALL RANDOM_NUMBER (rn)
       bins = MAX(1, INT(rn*8))
       CALL create_bl_distribution (block_distribution,&
            block_size, elements, bins, error)
       CALL array_release (block_distribution)
       CALL array_release (block_size)
    ENDDO
  END SUBROUTINE test_create_bl_distribution

! *****************************************************************************
!> \brief Distributes elements into blocks and into bins
!> \param[out] block_distribution       block distribution to bins
!> \param[out] block_size       sizes of blocks
!> \param[in] elements          number of elements to bin
!> \param[in] nbins             number of bins
! *****************************************************************************
  SUBROUTINE create_bl_distribution (block_distribution,&
       block_size, nelements, nbins, error)
    TYPE(array_i1d_obj), INTENT(OUT)         :: block_distribution, block_size
    INTEGER, INTENT(IN)                      :: nelements, nbins
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_bl_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER :: bin, blk_layer, element_stack, els, estimated_blocks, &
      max_block_size = 32, max_blocks_per_bin, max_elements_per_bin, &
      max_elements_per_block, nblks, stat
    INTEGER, DIMENSION(:), POINTER           :: blk_dist, blk_sizes

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL array_nullify (block_distribution)
    CALL array_nullify (block_size)
    ! Define the sizes on which we build the distribution.
    max_elements_per_bin = CEILING (REAL(nelements)/REAL(nbins))
    max_elements_per_block = MIN (max_elements_per_bin, max_block_size)
    max_blocks_per_bin = CEILING (REAL(max_elements_per_bin)/REAL(max_elements_per_block))
    IF (dbg) bin = max_elements_per_block
    max_elements_per_block = CEILING (REAL(max_elements_per_bin)/REAL(max_blocks_per_bin))
    IF (max_elements_per_block*max_blocks_per_bin*(nbins-1) .GE. nelements) THEN
       max_elements_per_block = MIN(1, max_elements_per_block-1)
       max_blocks_per_bin = max_blocks_per_bin+1
       DBG " small number FIX"
    ENDIF
    IF (dbg) THEN
       WRITE(*,'(1X,A,1X,A,I7,A,I7,A)')routineP,"For",nelements,&
            " elements and",nbins," bins"
       WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
            max_elements_per_bin," max elements/bin"
       WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
            bin," true max elements/block"
       WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
            max_blocks_per_bin," max blocks/bin"
       WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
            max_elements_per_block," max elements/block"
    ENDIF
    ! Prepare block distribution and sizes array, which may be
    ! oversized at this point.
    estimated_blocks = max_blocks_per_bin * nbins
    ALLOCATE (blk_dist (estimated_blocks), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, "blk_dist", estimated_blocks)
    ALLOCATE (blk_sizes (estimated_blocks), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, "blk_sizes", estimated_blocks)
    !
    element_stack = 0
    nblks = 0
    DO blk_layer = 1, max_blocks_per_bin
       DO bin = 0, nbins-1
          els = MIN (max_elements_per_block, nelements-element_stack)
          IF (els .GT. 0) THEN
             element_stack = element_stack + els
             nblks = nblks + 1
             blk_dist(nblks) = bin
             blk_sizes(nblks) = els
             IF (dbg) WRITE(*,'(1X,A,I5,A,I5,A,I5)')routineP//" Assigning",&
                  els," elements as block",nblks," to bin",bin
          ENDIF
       ENDDO
    ENDDO
    ! Create the output arrays.
    DBG "Wound up with",nblks,'blocks instead of',estimated_blocks
    IF (nblks .EQ. estimated_blocks) THEN
       CALL array_new (block_distribution, blk_dist, gift=.TRUE.)
       CALL array_new (block_size, blk_sizes, gift=.TRUE.)
    ELSE
       CALL array_new (block_distribution, blk_dist(1:nblks), lb=1)
       CALL array_new (block_size, blk_sizes(1:nblks), lb=1)
    ENDIF
1579 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
    IF (dbg) THEN
       WRITE(*,'(1X,A,A)')routineP//" Distribution"
       WRITE(*,1579)array_data (block_distribution)
       WRITE(*,'(1X,A,A)')routineP//" Sizes"
       WRITE(*,1579)array_data (block_size)
    ENDIF
  END SUBROUTINE create_bl_distribution


! *****************************************************************************
!> \brief Create the distribution_2d for the right matrix in a matrix-matrix
!>        multiplication, given the distribution_2d of the left matrix.
!> \par The new row distribution for the right matrix mirrors the column
!>      distribution of the left matrix, while the column distribution is
!>      created so that it is appropriate to the parallel environment. The
!>      columns are blocked together and each block is assigned to its own
!>      atom type.
!> \param[out] dist2d         new distribution_2d for the right matrix
!> \param[in] dist2d_left     the distribution_2d of the left matrix
!> \param[in] ncolumns        number of columns in right matrix
!> \param row_blk_sizes, col_blk_sizes  (optional) sizes of blocks for the
!>                                      row or column. Calculated if
!>                                      specified but not associated.
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist2d_r(dist2d, dist2d_left, ncolumns,&
       col_blk_size, error)
    TYPE(distribution_2d_type), POINTER      :: dist2d, dist2d_left
    INTEGER, INTENT(IN)                      :: ncolumns
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: col_blk_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist2d_r', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: blksize_default = 3

    INTEGER :: at, at_high, at_low, column_stack, i, k, mynod, mypcol, &
      mypcol_left, myprow, myprow_left, nblkcols, nblkcols_left, nblkrows, &
      nblks, ncols, npcols, npcols_left, nprows, nprows_left, pcol, pcsize, &
      stat
    INTEGER, DIMENSION(:), POINTER :: col_dist, col_dist_left, &
      n_at_local_cols, n_at_local_cols_left, n_at_local_rows, row_dist, &
      row_dist_left
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi, mpi2blacs
    LOGICAL                                  :: dbg
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      POINTER                                :: at_local_cols, &
                                                at_local_cols_left, &
                                                at_local_rows
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env, blacs_env_right
    TYPE(cp_para_env_type), POINTER          :: para_env

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    IF (dbg) CALL distribution_2d_write(dist2d_left, 6,&
         long_description=.TRUE., error=error, local=.TRUE.)
    IF (dbg) &
         CALL distribution_2d_write(dist2d_left, 6, long_description=.TRUE.,&
         error=error, local=.TRUE.)
    CALL distribution_2d_get(dist2d_left, error=error,&
         row_distribution=row_dist_left, col_distribution=col_dist_left,&
         local_cols=at_local_cols_left, n_local_cols=n_at_local_cols_left,&
         n_col_distribution=nblkcols_left, blacs_env=blacs_env)
    CALL get_blacs_info(blacs_env, my_process_row=myprow_left,&
         my_process_column=mypcol_left, number_of_process_rows=nprows_left,&
         number_of_process_columns=npcols_left, para_env=para_env,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    !  myprow = mypcol_left
    !  mypcol = myprow_left
    !  nprows = npcols_left
    !  npcols = nprows_left

    !  myprow = myprow_left
    !  mypcol = mypcol_left
    !  nprows = nprows_left
    !  npcols = npcols_left

    !
    !> Create the transposed blacs_env
    IF (dbg) WRITE(*,*)routineP//' LEFT prows, pcols:', nprows_left, npcols_left
    CALL cp_transpose_blacs_env(blacs_env, blacs_env_right, error)
    CALL get_blacs_info(blacs_env_right, my_process_row=myprow,&
         my_process_column=mypcol, number_of_process_rows=nprows,&
         number_of_process_columns=npcols,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    IF (dbg) WRITE(*,*)routineP//' RIGHT prows, pcols:', nprows, npcols
    IF (dbg) WRITE(*,*)routineP//' myprow, mypcol',myprow,'/',nprows,';',&
         mypcol,'/',npcols
    IF (dbg) CALL cp_blacs_env_write(blacs_env_right, 6, error)
    !
    mynod = para_env%mepos
    !> Handle row distribution. This is based on the column distribution of the
    !> left matrix.
    ALLOCATE(at_local_rows(LBOUND(at_local_cols_left,1):&
         UBOUND(at_local_cols_left,1)), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows',&
         SIZE(at_local_cols_left,1))
    at_low = LBOUND(n_at_local_cols_left,1) 
    at_high = UBOUND(n_at_local_cols_left,1)
    ALLOCATE(n_at_local_rows(at_low:at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_rows',&
         SIZE(n_at_local_cols_left,1))
    nblkrows = nblkcols_left
    ALLOCATE(row_dist(nblkrows), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'row_dist', nblkrows)
    n_at_local_rows(:) = n_at_local_cols_left(:)
    !> The rows and columns of the right matrix are the transposes of the left,
    !> so we can safely use the right-matrix's column data for our row data.
    !  IF(nprows.EQ.npcols) THEN
    !at_local_rows%array(:) = at_local_cols_left%array(:)
    DO i = LBOUND(at_local_cols_left,1), UBOUND(at_local_cols_left,1)
       ALLOCATE(at_local_rows(i)%array(&
            LBOUND(at_local_cols_left(i)%array,1):&
            UBOUND(at_local_cols_left(i)%array,1)), stat=stat)
       IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows(?)%array',&
            SIZE(at_local_cols_left(i)%array))
       at_local_rows(i)%array(:) = at_local_cols_left(i)%array(:)
    ENDDO
    n_at_local_rows(:) = n_at_local_cols_left(:)
    row_dist(:) = col_dist_left(:)
    ! Handle column distribution.
    ! Estimate the number of full columns per processor column, but it
    ! must be at least 1.
    pcsize = CEILING((1.0*ncolumns)/(1.0*npcols)) ! Full columns per processor column.
    column_stack = ncolumns
    nblkcols = 0
    DO pcol = 0, npcols-1
       at_low = nblkcols+1
       at_high = nblkcols
       IF (column_stack .GT. 0) THEN
          ncols = MIN(column_stack, pcsize)
          column_stack = column_stack - ncols
          nblks = CEILING ((1.0*ncols)/(1.0*blksize_default))
          nblkcols = nblkcols + nblks
          at_high = nblkcols
       ENDIF
       IF (dbg) WRITE(*,*)routineP//' pcol',pcol,'has blkcols',at_low,':',&
            at_high,'col stack is',column_stack
    ENDDO

    IF (dbg) WRITE(*,*)routineP//' There are',nblkcols,'blocked columns in',&
         npcols,'pcols for',ncolumns,'full columns.'

    at_low = 1
    at_high = nblkcols
    ALLOCATE(col_dist(nblkcols), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'col_dist', nblkcols)
    col_dist(:) = -999
    ALLOCATE(n_at_local_cols(at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_cols', at_high)
    ALLOCATE(at_local_cols(at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols', at_high)
    n_at_local_cols(:) = 0

    IF (PRESENT(col_blk_size)) THEN
       ALLOCATE(col_blk_size(nblkcols), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP,&
            'col_blk_begin, col_blk_end', nblkcols)
    ENDIF
    column_stack = ncolumns
    nblkcols = 0
    i = 1
    DO pcol = 0, npcols-1
       at_low = nblkcols+1
       at_high = nblkcols
       nblks = 0
       IF (column_stack .GT. 0) THEN
          ncols = MIN(column_stack, pcsize)
          column_stack = column_stack - ncols
          nblks = CEILING ((1.0*ncols)/(1.0*blksize_default))
          nblkcols = nblkcols + nblks
          at_high = nblkcols
       ENDIF
       IF (dbg) WRITE(*,*)routineP//' pcol',pcol,'has blkcols',at_low,':',&
            at_high,'col stack is',column_stack
       col_dist(at_low:at_high) = pcol
       DO at = at_low, at_high
          IF (pcol.EQ.mypcol) THEN
             ALLOCATE(at_local_cols(at)%array(1), stat=stat)
             IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols(?)', 1)
             at_local_cols(at)%array(1) = at
             n_at_local_cols(at) = 1
          ELSE
             ALLOCATE(at_local_cols(at)%array(0), stat=stat)
             IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols(?)', 1)
          ENDIF
          IF (PRESENT(col_blk_size)) THEN
             k = MIN(ncols,blksize_default)
             col_blk_size(at) = k
             i = i + k
             ncols = ncols - k
             IF (dbg) WRITE(*,*)routineP//' col',at,'has',k-1,&
                  'full colums; ncols, i',ncols, i,&
                  'nblks',nblks
          ENDIF
       ENDDO
    ENDDO

    IF(dbg) WRITE(*,*)routineP//' Col. distribution',col_dist
    IF (PRESENT (col_blk_size)) THEN
       IF(dbg) WRITE(*,*)routineP//' Col. block sizes',col_blk_size
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Col. at_local_cols=',n_at_local_cols
    IF(dbg) WRITE(*,*)routineP//' Col. n_at_local_cols=',n_at_local_cols

    NULLIFY(dist2d)
    CALL distribution_2d_create(dist2d, blacs_env_right, error=error,&
         local_rows_ptr=at_local_rows, n_local_rows=n_at_local_rows,&
         local_cols_ptr=at_local_cols, n_local_cols=n_at_local_cols,&
         row_distribution_ptr=row_dist, col_distribution_ptr=col_dist)
    DEALLOCATE(n_at_local_rows)
    DEALLOCATE(n_at_local_cols)
    !CALL cp_blacs_env_release(blacs_env_right, error)
    IF (dbg) CALL distribution_2d_write(dist2d, 6, long_description=.TRUE.,&
        error=error, local=.TRUE.)
    IF (dbg) &
         CALL distribution_2d_write(dist2d, 6, long_description=.TRUE.,&
         error=error, local=.TRUE.)
  END SUBROUTINE dbcsr_create_dist2d_r

! *****************************************************************************
!> \brief Create a transposed distribution_2d
!> \param[out] dist2d         transposed distribution_2d
!> \param[in] dist2d_normal   the base distribution_2d
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist2d_transposed(dist2d, dist2d_normal,  error)
    TYPE(distribution_2d_type), POINTER      :: dist2d, dist2d_normal
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_create_dist2d_transposed', &
      routineP = moduleN//':'//routineN

    INTEGER :: at_high, at_low, i, mynod, mypcol_normal, mypcol_tr, &
      myprow_normal, myprow_tr, nblkcols_normal, nblkcols_tr, &
      nblkrows_normal, nblkrows_tr, npcols_normal, npcols_tr, nprows_normal, &
      nprows_tr, stat
    INTEGER, DIMENSION(:), POINTER :: col_dist_normal, col_dist_tr, &
      n_at_local_cols_normal, n_at_local_cols_tr, n_at_local_rows_normal, &
      n_at_local_rows_tr, row_dist_normal, row_dist_tr
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi, mpi2blacs
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      POINTER                                :: at_local_cols_normal, &
                                                at_local_cols_tr, &
                                                at_local_rows_normal, &
                                                at_local_rows_tr
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_normal, blacs_env_tr
    TYPE(cp_para_env_type), POINTER          :: para_env

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    IF (dbg) CALL distribution_2d_write(dist2d_normal, 6,&
         long_description=.TRUE.,&
         error=error, local=.TRUE.)
    CALL dbcsr_set_debug(dbg)
    CALL distribution_2d_get(dist2d_normal, error=error,&
         row_distribution=row_dist_normal, col_distribution=col_dist_normal,&
         local_rows=at_local_rows_normal, n_local_rows=n_at_local_rows_normal,&
         local_cols=at_local_cols_normal, n_local_cols=n_at_local_cols_normal,&
         n_row_distribution=nblkrows_normal,&
         n_col_distribution=nblkcols_normal,&
         blacs_env=blacs_env_normal)
    CALL get_blacs_info(blacs_env_normal, my_process_row=myprow_normal,&
         my_process_column=mypcol_normal, number_of_process_rows=nprows_normal,&
         number_of_process_columns=npcols_normal, para_env=para_env,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    ! Create the transposed blacs_env
    CALL cp_transpose_blacs_env(blacs_env_normal, blacs_env_tr, error)
    CALL get_blacs_info(blacs_env_tr, my_process_row=myprow_tr,&
         my_process_column=mypcol_tr, number_of_process_rows=nprows_tr,&
         number_of_process_columns=npcols_tr,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    mynod = para_env%mepos
    !> Handle row distribution. This is based on the column distribution of the
    !> given matrix.
    at_low = LBOUND(n_at_local_cols_normal,1) 
    at_high = UBOUND(n_at_local_cols_normal,1)
    ALLOCATE(at_local_rows_tr(at_low:at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows_tr',&
         at_high-at_low+1)
    ALLOCATE(n_at_local_rows_tr(at_low:at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_rows_tr',&
         at_high-at_low+1)
    nblkrows_tr = nblkcols_normal
    ALLOCATE(row_dist_tr(nblkrows_tr), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'row_dist_tr', nblkrows_tr)
    n_at_local_rows_tr(:) = n_at_local_cols_normal(:)
    DO i = at_low, at_high
       ALLOCATE(at_local_rows_tr(i)%array(&
            LBOUND(at_local_cols_normal(i)%array,1):&
            UBOUND(at_local_cols_normal(i)%array,1)), stat=stat)
       IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows_tr(?)%array',&
            SIZE(at_local_cols_normal(i)%array))
       at_local_rows_tr(i)%array(:) = at_local_cols_normal(i)%array(:)
    END DO
    n_at_local_rows_tr(:) = n_at_local_cols_normal(:)
    row_dist_tr(:) = col_dist_normal(:)
    !
    ! Handle column distribution.
    at_low = LBOUND(n_at_local_rows_normal,1) 
    at_high = UBOUND(n_at_local_rows_normal,1)
    ALLOCATE(at_local_cols_tr(at_low:at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols_tr',&
         at_high-at_low+1)
    ALLOCATE(n_at_local_cols_tr(at_low:at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_cols_tr',&
         at_high-at_low+1)
    nblkcols_tr = nblkrows_normal
    ALLOCATE(col_dist_tr(nblkcols_tr), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'col_dist_tr', nblkcols_tr)
    n_at_local_cols_tr(:) = n_at_local_rows_normal(:)
    DO i = at_low, at_high
       ALLOCATE(at_local_cols_tr(i)%array(&
            LBOUND(at_local_rows_normal(i)%array,1):&
            UBOUND(at_local_rows_normal(i)%array,1)), stat=stat)
       IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols_tr(?)%array',&
            SIZE(at_local_rows_normal(i)%array))
       at_local_cols_tr(i)%array(:) = at_local_rows_normal(i)%array(:)
    END DO
    n_at_local_cols_tr(:) = n_at_local_rows_normal(:)
    col_dist_tr(:) = row_dist_normal(:)
    !
    NULLIFY(dist2d)
    CALL distribution_2d_create(dist2d, blacs_env_tr, error=error,&
         local_rows_ptr=at_local_rows_tr, n_local_rows=n_at_local_rows_tr,&
         local_cols_ptr=at_local_cols_tr, n_local_cols=n_at_local_cols_tr,&
         row_distribution_ptr=row_dist_tr, col_distribution_ptr=col_dist_tr)
    DEALLOCATE(n_at_local_rows_tr)
    DEALLOCATE(n_at_local_cols_tr)
    IF (dbg) CALL distribution_2d_write(dist2d, 6, long_description=.TRUE.,&
        error=error, local=.TRUE.)
  END SUBROUTINE dbcsr_create_dist2d_transposed



! *****************************************************************************
!> \brief Create the distribution_2d for the right matrix in a matrix-matrix
!>        multiplication, given the distribution_2d of the left matrix.
!> \par The new row distribution for the right matrix is the same as the row
!>      distribution of the left matrix, while the column distribution is
!>      created so that it is appropriate to the parallel environment. The
!>      columns are blocked together and each block is assigned to its own
!>      atom type.
!> \param[out] dist2d         new distribution_2d for the right matrix
!> \param[in] dist2d_left     the distribution_2d of the left matrix
!> \param[in] ncolumns        number of columns in right matrix
!> \param row_blk_sizes, col_blk_sizes  (optional) sizes of blocks for the
!>                                      row or column. Calculated if
!>                                      specified but not associated.
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist2d_r_unrot(dist2d, dist2d_left, ncolumns,&
       col_blk_size, error)
    TYPE(distribution_2d_type), POINTER      :: dist2d, dist2d_left
    INTEGER, INTENT(IN)                      :: ncolumns
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: col_blk_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist2d_r_unrot', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: blksize_default = 3

    INTEGER :: at, at_high, at_low, column_stack, i, k, mynod, mypcol, &
      mypcol_left, myprow, myprow_left, nblkcols, nblkrows, nblkrows_left, &
      nblks, ncols, npcols, npcols_left, nprows, nprows_left, pcol, pcsize, &
      stat
    INTEGER, DIMENSION(:), POINTER :: col_dist, col_dist_left, &
      n_at_local_cols, n_at_local_rows, n_at_local_rows_left, row_dist, &
      row_dist_left
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi, mpi2blacs
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      POINTER                                :: at_local_cols, at_local_rows, &
                                                at_local_rows_left
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env, blacs_env_right
    TYPE(cp_para_env_type), POINTER          :: para_env

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    IF (dbg) CALL distribution_2d_write(dist2d_left, 6,&
         long_description=.TRUE., error=error, local=.TRUE.)
    IF (dbg) &
         CALL distribution_2d_write(dist2d_left, 6, long_description=.TRUE.,&
         error=error, local=.TRUE.)
    CALL distribution_2d_get(dist2d_left, error=error,&
         row_distribution=row_dist_left, col_distribution=col_dist_left,&
         local_rows=at_local_rows_left, n_local_rows=n_at_local_rows_left,&
         n_row_distribution=nblkrows_left, blacs_env=blacs_env)
    CALL get_blacs_info(blacs_env, my_process_row=myprow_left,&
         my_process_column=mypcol_left, number_of_process_rows=nprows_left,&
         number_of_process_columns=npcols_left, para_env=para_env,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    !
    IF (dbg) WRITE(*,*)routineP//' LEFT prows, pcols:', nprows_left, npcols_left
    blacs_env_right => blacs_env
    !CALL cp_blacs_env_retain(blacs_env, error)
    CALL get_blacs_info(blacs_env_right, my_process_row=myprow,&
         my_process_column=mypcol, number_of_process_rows=nprows,&
         number_of_process_columns=npcols,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    IF (dbg) WRITE(*,*)routineP//' RIGHT prows, pcols:', nprows, npcols
    IF (dbg) WRITE(*,*)routineP//' myprow, mypcol',myprow,'/',nprows,';',&
         mypcol,'/',npcols
    IF (dbg) CALL cp_blacs_env_write(blacs_env_right, 6, error)
    !
    mynod = para_env%mepos
    !> Handle row distribution. This is based on the column distribution of the
    !> left matrix.
    ALLOCATE(at_local_rows(LBOUND(at_local_rows_left,1):&
         UBOUND(at_local_rows_left,1)), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows',&
         SIZE(at_local_rows_left,1))
    at_low = LBOUND(n_at_local_rows_left,1) 
    at_high = UBOUND(n_at_local_rows_left,1)
    ALLOCATE(n_at_local_rows(at_low:at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_rows',&
         SIZE(n_at_local_rows_left,1))
    nblkrows = nblkrows_left
    ALLOCATE(row_dist(nblkrows), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'row_dist', nblkrows)
    n_at_local_rows(:) = n_at_local_rows_left(:)
    DO i = LBOUND(at_local_rows_left,1), UBOUND(at_local_rows_left,1)
       ALLOCATE(at_local_rows(i)%array(&
            LBOUND(at_local_rows_left(i)%array,1):&
            UBOUND(at_local_rows_left(i)%array,1)), stat=stat)
       IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows(?)%array',&
            SIZE(at_local_rows_left(i)%array))
       at_local_rows(i)%array(:) = at_local_rows_left(i)%array(:)
    ENDDO
    n_at_local_rows(:) = n_at_local_rows_left(:)
    row_dist(:) = row_dist_left(:)
    ! Handle column distribution.
    ! Estimate the number of full columns per processor column, but it
    ! must be at least 1.
    pcsize = CEILING((1.0*ncolumns)/(1.0*npcols)) ! Full columns per processor column.
    column_stack = ncolumns
    nblkcols = 0
    DO pcol = 0, npcols-1
       at_low = nblkcols+1
       at_high = nblkcols
       IF (column_stack .GT. 0) THEN
          ncols = MIN(column_stack, pcsize)
          column_stack = column_stack - ncols
          nblks = CEILING ((1.0*ncols)/(1.0*blksize_default))
          nblkcols = nblkcols + nblks
          at_high = nblkcols
       ENDIF
       IF (dbg) WRITE(*,*)routineP//' pcol',pcol,'has blkcols',at_low,':',&
            at_high,'col stack is',column_stack
    ENDDO

    IF (dbg) WRITE(*,*)routineP//' There are',nblkcols,'blocked columns in',&
         npcols,'pcols for',ncolumns,'full columns.'

    at_low = 1
    at_high = nblkcols
    ALLOCATE(col_dist(nblkcols), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'col_dist', nblkcols)
    col_dist(:) = -999
    ALLOCATE(n_at_local_cols(at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_cols', at_high)
    ALLOCATE(at_local_cols(at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols', at_high)
    n_at_local_cols(:) = 0

    IF (PRESENT(col_blk_size)) THEN
       ALLOCATE(col_blk_size(nblkcols), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP,&
            'col_blk_begin, col_blk_end', nblkcols)
    ENDIF
    column_stack = ncolumns
    nblkcols = 0
    i = 1
    DO pcol = 0, npcols-1
       at_low = nblkcols+1
       at_high = nblkcols
       nblks = 0
       IF (column_stack .GT. 0) THEN
          ncols = MIN(column_stack, pcsize)
          column_stack = column_stack - ncols
          nblks = CEILING ((1.0*ncols)/(1.0*blksize_default))
          nblkcols = nblkcols + nblks
          at_high = nblkcols
       ENDIF
       IF (dbg) WRITE(*,*)routineP//' pcol',pcol,'has blkcols',at_low,':',&
            at_high,'col stack is',column_stack
       col_dist(at_low:at_high) = pcol
       DO at = at_low, at_high
          IF (pcol.EQ.mypcol) THEN
             ALLOCATE(at_local_cols(at)%array(1), stat=stat)
             IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols(?)', 1)
             at_local_cols(at)%array(1) = at
             n_at_local_cols(at) = 1
          ELSE
             ALLOCATE(at_local_cols(at)%array(0), stat=stat)
             IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols(?)', 1)
          ENDIF
          IF (PRESENT(col_blk_size)) THEN
             k = MIN(ncols,blksize_default)
             col_blk_size(at) = k
             i = i + k
             ncols = ncols - k
             IF (dbg) WRITE(*,*)routineP//' col',at,'has',k-1,&
                  'full colums; ncols, i',ncols, i,&
                  'nblks',nblks
          ENDIF
       ENDDO
    ENDDO

    IF(dbg) WRITE(*,*)routineP//' Col. distribution',col_dist
    IF (PRESENT (col_blk_size)) THEN
       IF(dbg) WRITE(*,*)routineP//' Col. block sizes',col_blk_size
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Col. at_local_cols=',n_at_local_cols
    IF(dbg) WRITE(*,*)routineP//' Col. n_at_local_cols=',n_at_local_cols

    NULLIFY(dist2d)
    CALL distribution_2d_create(dist2d, blacs_env_right, error=error,&
         local_rows_ptr=at_local_rows, n_local_rows=n_at_local_rows,&
         local_cols_ptr=at_local_cols, n_local_cols=n_at_local_cols,&
         row_distribution_ptr=row_dist, col_distribution_ptr=col_dist)
    DEALLOCATE(n_at_local_rows)
    DEALLOCATE(n_at_local_cols)
    !CALL cp_blacs_env_release(blacs_env_right, error)
    IF (dbg) CALL distribution_2d_write(dist2d, 6, long_description=.TRUE.,&
        error=error, local=.TRUE.)
    IF (dbg) &
         CALL distribution_2d_write(dist2d, 6, long_description=.TRUE.,&
         error=error, local=.TRUE.)
  END SUBROUTINE dbcsr_create_dist2d_r_unrot


!! *****************************************************************************
!!> \brief Create a distribution_2d for a row as a selection of an existing
!!>        distribution_2d.
!!> \par The column distribution is the same as the original distribution_2d
!!> \param[out] dist2d         new distribution_2d for a row matrix
!!> \param[in] dist2d_left     the distribution_2d of the left matrix
!!> \param error     cp2k error
!! *****************************************************************************
!  SUBROUTINE dbcsr_create_row_dist2d_from_dist2d(dist2d, dist2d_main, error)
!    TYPE(distribution_2d_type), POINTER      :: dist2d, dist2d_main
!    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: row_blk_size
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_row_dist2d_from_dist2d', &
!      routineP = moduleN//':'//routineN
!    INTEGER, PARAMETER                       :: blksize_default = 3
!
!    INTEGER :: at, at_high, at_low, column_stack, i, k, mynod, mypcol, &
!      mypcol_main, myprow, myprow_main, nblkcols, nblkrows, nblkrows_main, &
!      nblks, ncols, npcols, npcols_main, nprows, nprows_main, pcol, pcsize
!      stat
!    INTEGER                                  :: stat
!    INTEGER                                  :: nrows
!    INTEGER, DIMENSION(:), POINTER :: col_dist, col_dist_main, &
!      n_at_local_cols, n_at_local_rows, n_at_local_rows_main, row_dist, &
!      row_dist_main
!    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi, mpi2blacs
!    TYPE(cp_1d_i_p_type), DIMENSION(:), &
!      POINTER                                :: at_local_cols, at_local_rows, &
!                                                at_local_rows_main
!    TYPE(cp_blacs_env_type), POINTER         :: blacs_env, blacs_env_right
!    TYPE(cp_para_env_type), POINTER          :: para_env
!
!    DEBUG_HEADER
!!   ---------------------------------------------------------------------------
!    DEBUG_BODYY
!    IF (dbg) &
!         CALL distribution_2d_write(dist2d_main, 6, long_description=.TRUE.,&
!         error=error, local=.TRUE.)
!    CALL distribution_2d_get(dist2d_main, error=error,&
!         row_distribution=row_dist_main, col_distribution=col_dist_main,&
!         local_rows=at_local_rows_main, n_local_rows=n_at_local_rows_main,&
!         n_row_distribution=nblkrows_main, blacs_env=blacs_env)
!    CALL get_blacs_info(blacs_env, my_process_row=myprow_main,&
!         my_process_column=mypcol_main, number_of_process_rows=nprows_main,&
!         number_of_process_columns=npcols_main, para_env=para_env,&
!         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
!    !
!    IF (dbg) WRITE(*,*)routineP//' MAIN prows, pcols:', nprows_left, npcols_left
!    blacs_env_right => blacs_env
!    CALL cp_blacs_env_retain(blacs_env, error)
!    CALL get_blacs_info(blacs_env_right, my_process_row=myprow,&
!         my_process_column=mypcol, number_of_process_rows=nprows,&
!         number_of_process_columns=npcols,&
!         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
!    IF (dbg) WRITE(*,*)routineP//' RIGHT prows, pcols:', nprows, npcols
!    IF (dbg) WRITE(*,*)routineP//' myprow, mypcol',myprow,'/',nprows,';',&
!         mypcol,'/',npcols
!    IF (dbg) CALL cp_blacs_env_write(blacs_env_right, 6, error)
!    !
!    mynod = para_env%mepos
!    !> Handle row distribution. This is based on the column distribution of the
!    !> left matrix.
!    ALLOCATE(at_local_rows(LBOUND(at_local_rows_left,1):&
!         UBOUND(at_local_rows_left,1)), stat=stat)
!    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows',&
!         SIZE(at_local_rows_left,1))
!    at_low = LBOUND(n_at_local_rows_left,1) 
!    at_high = UBOUND(n_at_local_rows_left,1)
!    ALLOCATE(n_at_local_rows(at_low:at_high), stat=stat)
!    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_rows',&
!         SIZE(n_at_local_rows_left,1))
!    nblkrows = nblkrows_left
!    ALLOCATE(row_dist(nblkrows), stat=stat)
!    IF(stat /= 0) CALL stop_memory(routineP, 'row_dist', nblkrows)
!    n_at_local_rows(:) = n_at_local_rows_left(:)
!    DO i = LBOUND(at_local_rows_left,1), UBOUND(at_local_rows_left,1)
!       ALLOCATE(at_local_rows(i)%array(&
!            LBOUND(at_local_rows_left(i)%array,1):&
!            UBOUND(at_local_rows_left(i)%array,1)), stat=stat)
!       IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows(?)%array',&
!            SIZE(at_local_rows_left(i)%array))
!       at_local_rows(i)%array(:) = at_local_rows_left(i)%array(:)
!    ENDDO
!    n_at_local_rows(:) = n_at_local_rows_left(:)
!    row_dist(:) = row_dist_left(:)
!    ! Handle column distribution.
!    ! Estimate the number of full columns per processor column, but it
!    ! must be at least 1.
!    pcsize = CEILING((1.0*ncolumns)/(1.0*npcols)) ! Full columns per processor column.
!    column_stack = ncolumns
!    nblkcols = 0
!    DO pcol = 0, npcols-1
!       at_low = nblkcols+1
!       at_high = nblkcols
!       IF (column_stack .GT. 0) THEN
!          ncols = MIN(column_stack, pcsize)
!          column_stack = column_stack - ncols
!          nblks = CEILING ((1.0*ncols)/(1.0*blksize_default))
!          nblkcols = nblkcols + nblks
!          at_high = nblkcols
!       ENDIF
!       IF (dbg) WRITE(*,*)routineP//' pcol',pcol,'has blkcols',at_low,':',&
!            at_high,'col stack is',column_stack
!    ENDDO
!
!    IF (dbg) WRITE(*,*)routineP//' There are',nblkcols,'blocked columns in',&
!         npcols,'pcols for',ncolumns,'full columns.'
!
!    at_low = 1
!    at_high = nblkcols
!    ALLOCATE(col_dist(nblkcols), stat=stat)
!    IF(stat /= 0) CALL stop_memory(routineP, 'col_dist', nblkcols)
!    col_dist(:) = -999
!    ALLOCATE(n_at_local_cols(at_high), stat=stat)
!    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_cols', at_high)
!    ALLOCATE(at_local_cols(at_high), stat=stat)
!    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols', at_high)
!    n_at_local_cols(:) = 0
!
!    IF (PRESENT(col_blk_size)) THEN
!       ALLOCATE(col_blk_size(nblkcols), stat=stat)
!       IF (stat /= 0) CALL stop_memory(routineP,&
!            'col_blk_begin, col_blk_end', nblkcols)
!    ENDIF
!    column_stack = ncolumns
!    nblkcols = 0
!    i = 1
!    DO pcol = 0, npcols-1
!       at_low = nblkcols+1
!       at_high = nblkcols
!       nblks = 0
!       IF (column_stack .GT. 0) THEN
!          ncols = MIN(column_stack, pcsize)
!          column_stack = column_stack - ncols
!          nblks = CEILING ((1.0*ncols)/(1.0*blksize_default))
!          nblkcols = nblkcols + nblks
!          at_high = nblkcols
!       ENDIF
!       IF (dbg) WRITE(*,*)routineP//' pcol',pcol,'has blkcols',at_low,':',&
!            at_high,'col stack is',column_stack
!       col_dist(at_low:at_high) = pcol
!       DO at = at_low, at_high
!          IF (pcol.EQ.mypcol) THEN
!             ALLOCATE(at_local_cols(at)%array(1), stat=stat)
!             IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols(?)', 1)
!             at_local_cols(at)%array(1) = at
!             n_at_local_cols(at) = 1
!          ELSE
!             ALLOCATE(at_local_cols(at)%array(0), stat=stat)
!             IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols(?)', 1)
!          ENDIF
!          IF (PRESENT(col_blk_size)) THEN
!             k = MIN(ncols,blksize_default)
!             col_blk_size(at) = k
!             i = i + k
!             ncols = ncols - k
!             IF (dbg) WRITE(*,*)routineP//' col',at,'has',k-1,&
!                  'full colums; ncols, i',ncols, i,&
!                  'nblks',nblks
!          ENDIF
!       ENDDO
!    ENDDO
!
!    IF(dbg) WRITE(*,*)routineP//' Col. distribution',col_dist
!    IF (PRESENT (col_blk_size)) THEN
!       IF(dbg) WRITE(*,*)routineP//' Col. block sizes',col_blk_size
!    ENDIF
!    IF(dbg) WRITE(*,*)routineP//' Col. at_local_cols=',n_at_local_cols
!    IF(dbg) WRITE(*,*)routineP//' Col. n_at_local_cols=',n_at_local_cols
!
!    NULLIFY(dist2d)
!    CALL distribution_2d_create(dist2d, blacs_env_right, error=error,&
!         local_rows_ptr=at_local_rows, n_local_rows=n_at_local_rows,&
!         local_cols_ptr=at_local_cols, n_local_cols=n_at_local_cols,&
!         row_distribution_ptr=row_dist, col_distribution_ptr=col_dist)
!    DEALLOCATE(n_at_local_rows)
!    DEALLOCATE(n_at_local_cols)
!    !CALL cp_blacs_env_release(blacs_env_right, error)
!    IF (dbg) CALL distribution_2d_write(dist2d, 6, long_description=.TRUE.,&
!        error=error, local=.TRUE.)
!    IF (dbg) &
!         CALL distribution_2d_write(dist2d, 6, long_description=.TRUE.,&
!         error=error, local=.TRUE.)
!  END SUBROUTINE dbcsr_create_row_dist2d_from_dist2d


! *****************************************************************************
!> \brief Creates a new distribution for the right matrix in a matrix
!>        multiplication with unrotated grid.
!> \par The new row distribution for the right matrix is the same as the row
!>      distribution of the left matrix, while the column distribution is
!>      created so that it is appropriate to the parallel environment.
!> \param[out] dist_right     new distribution for the right matrix
!> \param[in] dist_left       the distribution of the left matrix
!> \param[in] ncolumns        number of columns in right matrix
!> \param[out] right_col_blk_sizes      sizes of blocks in the created column
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_r_unrot(dist_right, dist_left, ncolumns,&
       right_col_blk_sizes, error)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_right
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_left
    INTEGER, INTENT(IN)                      :: ncolumns
    TYPE(array_i1d_obj), INTENT(OUT)         :: right_col_blk_sizes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_r_unrot', &
      routineP = moduleN//':'//routineN

    TYPE(array_i1d_obj)                      :: right_col_dist
    TYPE(dbcsr_mp_obj)                       :: right_mp_env

!   ---------------------------------------------------------------------------

    right_mp_env = dbcsr_distribution_mp (dist_left)
    CALL create_bl_distribution (right_col_dist,&
         right_col_blk_sizes, ncolumns, dbcsr_mp_npcols (right_mp_env), error)
    CALL dbcsr_distribution_new (dist_right, right_mp_env,&
         dbcsr_distribution_row_dist (dist_left), right_col_dist,&
         local_rows=dbcsr_distribution_local_rows (dist_left))
    CALL array_release (right_col_dist)
  END SUBROUTINE dbcsr_create_dist_r_unrot

! *****************************************************************************
!> \brief Creates a new distribution for the right matrix in a matrix
!>        multiplication with rotated grid.
!> \par The new row distribution for the right matrix is the same as the
!>      column distribution of the left matrix, while the column distribution
!>      is created so that it is appropriate to the parallel environment.
!> \param[out] dist_right     new distribution for the right matrix
!> \param[in] dist_left       the distribution of the left matrix
!> \param[in] ncolumns        number of columns in right matrix
!> \param[out] right_col_blk_sizes      sizes of blocks in the created column
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_r_rot(dist_right, dist_left, ncolumns,&
       right_col_blk_sizes, error)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_right
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_left
    INTEGER, INTENT(IN)                      :: ncolumns
    TYPE(array_i1d_obj), INTENT(OUT)         :: right_col_blk_sizes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_r_rot', &
      routineP = moduleN//':'//routineN

    TYPE(array_i1d_obj)                      :: right_col_dist
    TYPE(dbcsr_mp_obj)                       :: right_mp_env

!   ---------------------------------------------------------------------------
!CALL test_create_bl_distribution (error)

    CALL dbcsr_mp_new_transposed (right_mp_env,&
         dbcsr_distribution_mp (dist_left))
    CALL create_bl_distribution (right_col_dist,&
         right_col_blk_sizes, ncolumns, dbcsr_mp_npcols (right_mp_env), error)
    CALL dbcsr_distribution_new (dist_right, right_mp_env,&
         dbcsr_distribution_col_dist (dist_left), right_col_dist,&
         local_rows=dbcsr_distribution_local_cols (dist_left))
    CALL dbcsr_mp_release (right_mp_env)
    CALL array_release (right_col_dist)
  END SUBROUTINE dbcsr_create_dist_r_rot

! *****************************************************************************
!> \brief Creates a new distribution for the left matrix in a matrix
!>        multiplication with unrotated grid.
!> \par The new column distribution for the left matrix is the same as the
!>      column distribution of the right matrix, while the row distribution is
!>      created so that it is appropriate to the parallel environment.
!> \param[out] dist_left      new distribution for the left matrix
!> \param[in] dist_right      old distribution of the right matrix
!> \param[in] nrows           number of rows in the left matrix
!> \param[out] left_row_blk_sizes       sizes of blocks in the created rows
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_l_unrot(dist_left, dist_right, nrows,&
       left_row_blk_sizes, error)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_left
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_right
    INTEGER, INTENT(IN)                      :: nrows
    TYPE(array_i1d_obj), INTENT(OUT)         :: left_row_blk_sizes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_l_unrot', &
      routineP = moduleN//':'//routineN

    TYPE(array_i1d_obj)                      :: left_row_dist
    TYPE(dbcsr_mp_obj)                       :: left_mp_env

!   ---------------------------------------------------------------------------

    left_mp_env = dbcsr_distribution_mp (dist_left)
    CALL create_bl_distribution (left_row_dist,&
         left_row_blk_sizes, nrows, dbcsr_mp_nprows (left_mp_env), error)
    CALL dbcsr_distribution_new (dist_left, left_mp_env,&
         left_row_dist, dbcsr_distribution_col_dist (dist_right),&
         local_cols=dbcsr_distribution_local_cols (dist_right))
    CALL array_release (left_row_dist)
  END SUBROUTINE dbcsr_create_dist_l_unrot

! *****************************************************************************
!> \brief Creates a new distribution for the left matrix in a matrix
!>        multiplication with rotated grid.
!> \par The new column distribution for the left matrix is the same as the
!>      row distribution of the right matrix, while the row distribution is
!>      created so that it is appropriate to the parallel environment.
!> \param[out] dist_left      new distribution for the left matrix
!> \param[in] dist_right      the distribution of the right matrix
!> \param[in] nrows           number of rows in left matrix
!> \param[out] left_row_blk_sizes      sizes of blocks in the created rows
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_l_rot(dist_right, dist_left, nrows,&
       left_row_blk_sizes, error)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_right
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_left
    INTEGER, INTENT(IN)                      :: nrows
    TYPE(array_i1d_obj), INTENT(OUT)         :: left_row_blk_sizes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_l_rot', &
      routineP = moduleN//':'//routineN

    TYPE(array_i1d_obj)                      :: left_row_dist
    TYPE(dbcsr_mp_obj)                       :: left_mp_env

!   ---------------------------------------------------------------------------

    CALL dbcsr_mp_new_transposed (left_mp_env,&
         dbcsr_distribution_mp (dist_right))
    CALL create_bl_distribution (left_row_dist,&
         left_row_blk_sizes, nrows, dbcsr_mp_nprows (left_mp_env), error)
    CALL dbcsr_distribution_new (dist_left, left_mp_env,&
         left_row_dist, dbcsr_distribution_row_dist (dist_right),&
         local_cols=dbcsr_distribution_local_rows (dist_right))
    CALL dbcsr_mp_release (left_mp_env)
    CALL array_release (left_row_dist)
  END SUBROUTINE dbcsr_create_dist_l_rot


! *****************************************************************************
!> \brief Creates a transposed BLACS environment.
!> \par This routine assumes that all BLACS environments are created equal.
!>      It just fiddles with the generated data.
!> \param[in] blacs_env       original BLACS environment.
!> \param[out] blacs_env_t    transposed BLACS environment.
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE cp_transpose_blacs_env(blacs_env, blacs_env_t, error)
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env, blacs_env_t
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_transpose_blacs_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi, mpi2blacs
    LOGICAL                                  :: dbg
    TYPE(cp_para_env_type), POINTER          :: para_env

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    CALL get_blacs_info(blacs_env, para_env=para_env,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    !
    !> Create a default
    NULLIFY(blacs_env_t)
    CALL cp_blacs_env_create(blacs_env_t, para_env, error)
    DEALLOCATE(blacs_env_t%blacs2mpi, blacs_env_t%mpi2blacs)
    blacs_env_t%num_pe(1:2) = blacs_env_t%num_pe(2:1:-1)
    blacs_env_t%mepos(1:2) = blacs_env_t%mepos(2:1:-1)
    !
    ALLOCATE(blacs_env_t%blacs2mpi(LBOUND(blacs2mpi,2):UBOUND(blacs2mpi,2),&
         &                         LBOUND(blacs2mpi,1):UBOUND(blacs2mpi,1)),&
         stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'blacs2mpi_t', SIZE(blacs2mpi))
    blacs_env_t%blacs2mpi = TRANSPOSE(blacs2mpi)
    !
    IF (dbg) WRITE(*,*)'blacs2mpi new',SIZE(blacs_env_t%blacs2mpi,1),&
         SIZE(blacs_env_t%blacs2mpi,2),&
         ';',blacs_env_t%blacs2mpi
    !
    ALLOCATE(blacs_env_t%mpi2blacs(LBOUND(mpi2blacs,2):UBOUND(mpi2blacs,2),&
         &                         LBOUND(mpi2blacs,1):UBOUND(mpi2blacs,1)),&
         stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'mpi2blacs_t', SIZE(mpi2blacs))
    blacs_env_t%mpi2blacs = TRANSPOSE(mpi2blacs)
  END SUBROUTINE cp_transpose_blacs_env


! *****************************************************************************
!> \brief Adds data to the index. Incleases the index size when neccessary.
!> \param[in,out] matrix      bcsr matrix
!> \param[in] slot  which index array to add (e.g., dbcsr_slot_row_blk_sizes)
!> \param[in] data  (optional) array holding the index data to add to the index
!>                  array
!> \param[in] reservation     (optional) only reserve space for subsequent
!>                            array
!> \param[in] extra (optional) reserve extra space for later additions
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_addto_index_array(matrix, slot, DATA, reservation, extra,&
     error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN)                      :: slot
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: DATA
    INTEGER, INTENT(IN), OPTIONAL            :: reservation, extra
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_addto_index_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: deplus, space, timing_handle, &
                                                ub, ub_new
    LOGICAL                                  :: dbg, failure = .FALSE.

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    CALL timeset(routineN, timing_handle)
    CPPrecondition(ASSOCIATED(matrix%index), cp_fatal_level, routineP, error, failure)
    CPPrecondition(UBOUND(matrix%index,1).GE.dbcsr_num_slots, cp_failure_level, routineP, error, failure)
    CALL cp_assert(PRESENT(DATA) .OR. PRESENT(reservation), cp_fatal_level,&
         cp_caller_error, routineP,&
         'Either an array or its size must be specified.', error=error)
    IF(dbg) WRITE(*,*)routineP//' index',&
         matrix%index(:dbcsr_num_slots)
    IF (PRESENT (reservation)) THEN
       space = reservation
    ELSE
       space = SIZE(DATA)
    ENDIF
    IF (PRESENT (extra)) THEN
       deplus = extra
    ELSE
       deplus = 0
    ENDIF
    ub = UBOUND(matrix%index,1)
    !> The data area was not defined or the new area is greater than the old.
    IF (matrix%index(slot).EQ.0 .OR.&
         space.GT.matrix%index(slot+1)-matrix%index(slot)+1) THEN
       IF(dbg) WRITE(*,*)routineP//' Slot',slot,'not filled, adding at',&
            matrix%index(dbcsr_slot_size)+1,'sized',space
       matrix%index(slot) = matrix%index(dbcsr_slot_size)+1
       matrix%index(slot+1) = matrix%index(slot) + space - 1
       matrix%index(dbcsr_slot_size) = matrix%index(slot+1)
    ENDIF
    ub_new = matrix%index(slot+1) + deplus
    IF(dbg) WRITE(*,*)routineP//' need',space,'at',matrix%index(slot),&
         'to',matrix%index(slot+1),'(',ub_new,')','have',ub
    IF (ub_new .GT. ub) THEN
       IF(dbg) WRITE(*,*)routineP//' Reallocating index to ubound', ub_new
       !CALL reallocate(matrix%index, 1, ub_new)
       CALL ensure_array_size(matrix%index, 1, ub_new, factor=1.0,&
            nocopy=.FALSE., special=matrix%special_memory, error=error)
       CALL dbcsr_repoint_index(matrix, error=error)
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Adding slot',slot,'at',&
         matrix%index(slot),'sized',space
    CALL dbcsr_repoint_index(matrix, slot, error)
    IF (PRESENT(DATA)) &
         matrix%index(matrix%index(slot):matrix%index(slot+1)) = DATA(:)
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_addto_index_array


! *****************************************************************************
!> \brief Updates the index pointers of a bcsr matrix
!> \param[in,out] m           matrix for which index pointers are updated
!> \param[in] slot            (optional) only repoint this index
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_repoint_index(m, slot, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: m
    INTEGER, INTENT(IN), OPTIONAL            :: slot
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_repoint_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: s
    LOGICAL                                  :: all

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    DBG 'Repointing index for matrix ',m%name,'with #blocks',&
         m%nblks,m%index(dbcsr_slot_nblks)
    IF (m%nblks .NE. m%index(dbcsr_slot_nblks)) THEN
       m%nblks = m%index(dbcsr_slot_nblks)
       m%nze = m%index(dbcsr_slot_nze)
    ENDIF
    all = .TRUE.
    IF (PRESENT (slot)) THEN
       all = .FALSE.
       s = slot
       DBG 'only slot',slot, m%index(slot)
    ELSE
       s = 0
    ENDIF
    DBGV 'index',m%index(1:dbcsr_num_slots)
    DBG 'all?',all,'s=',s

    IF (m%index(dbcsr_slot_row_p).GT.0.AND.all .OR.&
         s.EQ.dbcsr_slot_row_p) THEN
       DBG ' setting pointer row_p'
       m%row_p => m%index(m%index(dbcsr_slot_row_p):&
            &                  m%index(dbcsr_slot_row_p+1))
       DBG ' set pointer row_p, size',&
            m%index(dbcsr_slot_row_p+1)-m%index(dbcsr_slot_row_p)+1
    ENDIF
    IF (m%index(dbcsr_slot_col_i).GT.0.AND.all .OR.&
         s.EQ.dbcsr_slot_col_i) THEN
       DBG ' setting pointer col_i'
       m%col_i => m%index(m%index(dbcsr_slot_col_i):&
            &                  m%index(dbcsr_slot_col_i+1))
       DBG ' set pointer col_i, size',&
            m%index(dbcsr_slot_col_i+1)-m%index(dbcsr_slot_col_i)+1
    ENDIF
    IF (m%index(dbcsr_slot_blk_p).GT.0.AND.all .OR.&
         s.EQ.dbcsr_slot_blk_p) THEN
       DBG ' setting pointer blk_p'
       m%blk_p => m%index(m%index(dbcsr_slot_blk_p):&
            &                  m%index(dbcsr_slot_blk_p+1))
       DBG ' set pointer blk_p, size',&
            m%index(dbcsr_slot_blk_p+1)-m%index(dbcsr_slot_blk_p)+1
    ENDIF
    IF (all) THEN
       m%index(dbcsr_slot_nblks) = m%nblks
       m%index(dbcsr_slot_nze) = m%nze
    ENDIF
  END SUBROUTINE dbcsr_repoint_index


! *****************************************************************************
!> \brief Calculates the checksum of a DBCSR matrix.
!> \param[in] matrix          matrix
!> \param[out] chksum         calculated checksum
!> \param[in,out] error       cp2k error
! *****************************************************************************
  FUNCTION dbcsr_checksum(matrix, error) RESULT(checksum)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    TYPE(cp_error_type), INTENT(IN)          :: error
    REAL(KIND=dp)                            :: checksum

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_checksum', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bc, blk, blk_p, br, m, mn, n, &
                                                timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    REAL(KIND=dp)                            :: blk_cs, local_cs
    TYPE(cp_error_type)                      :: e

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    e = error
    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)
    local_cs = 0.0_dp
    DO br = 1, matrix%nblkrows_total
       m = row_blk_size(br)
       DO blk = matrix%row_p(br)+1, matrix%row_p(br+1)
          bc = matrix%col_i(blk)
          n = col_blk_size(bc)
          mn = m*n
          blk_p = ABS(matrix%blk_p(blk))
          ! Calculate DDOT
          !blk_cs = bl_dot(mn, matrix%data(blk_p:blk_p+mn-1), 1,&
          !     &              matrix%data(blk_p:blk_p+mn-1), 1)
          !write(*,*)matrix%data(blk_p:blk_p+mn-1)
          blk_cs = DOT_PRODUCT(matrix%data(blk_p:blk_p+mn-1),&
               &               matrix%data(blk_p:blk_p+mn-1))
          !IF (dbg .AND. bcsr_verbose)&
          !     CALL dbcsr_printmat(matrix%data(blk_p:blk_p+mn-1),m,n)
          DBGV ' chksum for (',br,',',bc,') at', blk_p,'l',mn,'= ', blk_cs
          local_cs = local_cs + blk_cs
       ENDDO
    ENDDO
    DBG 'local chksum=',local_cs
    CALL mp_sum(local_cs, dbcsr_mp_group (dbcsr_distribution_mp (matrix%dist)))
    checksum = local_cs
    DBGI 'local chksum=',local_cs,'global=',checksum
    CALL timestop(timing_handle)
  END FUNCTION dbcsr_checksum

! *****************************************************************************
!> \brief Prints a BCSR matrix (block-style, not full)
!> \param[in] matrix          matrix
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_print(matrix, error)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_print', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bc, blk, blk_p, br, m, mn, n, &
                                                timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: tr
    REAL(KIND=dp)                            :: blk_cs

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    WRITE(*,*)routineP//' Contents of matrix named ', matrix%name
    DBG 'index', matrix%index
    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)
    IF (matrix%nblks .GT. 0) THEN
       DO br = 1, matrix%nblkrows_total
          m = row_blk_size(br)
          DO blk = matrix%row_p(br)+1, matrix%row_p(br+1)
             bc = matrix%col_i(blk)
             n = col_blk_size(bc)
             mn = m*n
             blk_p = ABS(matrix%blk_p(blk))
             tr = matrix%blk_p(blk) .LT. 0
             blk_cs = DOT_PRODUCT(matrix%data(blk_p:blk_p+mn-1),&
                  &               matrix%data(blk_p:blk_p+mn-1))
             WRITE(*,*)routineP//' chksum for (',br,',',bc,') at',&
                  blk_p,'l',mn,'= ', blk_cs
             CALL &
                  dbcsr_printmat(matrix%data(blk_p:blk_p+mn-1),m,n, tr=tr)
          ENDDO
       ENDDO
    ENDIF
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_print



! *****************************************************************************
!> \brief Verify the correctness of a BCSR matrix.
!> \param[in] m     bcsr matrix
!> \param[in,out] error     cp2k error
!> \param[in] verbosity     how detailed errors are; 0=nothing; 1=summary at
!>                          end if matrix not consistent; 2=also individual
!>                          errors; 3=always print info about matrix; >3=even
!>                          more info
! *****************************************************************************
  SUBROUTINE dbcsr_verify_matrix(m, error, verbosity)
    TYPE(dbcsr_type), INTENT(IN)             :: m
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER, INTENT(IN), OPTIONAL            :: verbosity

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_verify_matrix', &
      r = moduleN//':'//routineN

    INTEGER :: bc, blk, blk_p, br, dbg, i, mb, mn, n, n_full_blocks_total, &
      n_have_blocks_local, n_have_blocks_total, prev_br, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    REAL(KIND=dp)                            :: sparsity_total

!   ---------------------------------------------------------------------------

    dbg = 2
    IF (PRESENT (verbosity)) dbg = verbosity
    CALL timeset(routineN, timing_handle)
    IF (dbg.GE.3) WRITE(*,'(1X,A,A,A,A)')r//'Matrix name: ', m%name,&
         " of type ",m%matrix_type
    IF (dbg.GE.3) THEN
       WRITE(*,'(1X,A,I5,"x",I5,A,I5,"x",I5)')r//' Size blocked',&
            m%nblkrows_total, m%nblkcols_total, ", full ",&
            m%nfullrows_total, m%nfullcols_total
    ENDIF
    row_blk_size => array_data (m%row_blk_size)
    col_blk_size => array_data (m%col_blk_size)
    !
    IF(m%matrix_type.EQ.'N' .OR. m%matrix_type.EQ.'T') THEN
       n_full_blocks_total = m%nblkrows_total * m%nblkcols_total
    ELSE
       CALL bcsr_assert(dbg, m%nblkrows_total.EQ.m%nblkcols_total, r,&
            'Symmetric matrix is not square')
       n_full_blocks_total = m%nblkrows_total * (m%nblkrows_total+1)/2
    ENDIF
    n_have_blocks_local = m%nblks
2045 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
2047 FORMAT (I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7)
    IF (dbg.GE.4) THEN
       WRITE(*,'(1X,A)')r//' index='
       WRITE(*,2045)m%index(:dbcsr_num_slots)
    ENDIF
    CALL bcsr_assert(dbg, m%index(1).GT.0, r, 'Index size 0')
    DO i = dbcsr_slot_row_p, dbcsr_num_slots
       CALL bcsr_assert(dbg, m%index(i) .GT. 0, r, 'Index member is 0')
       IF (i .LE. dbcsr_num_slots-2) &
            CALL bcsr_assert(dbg,m%index(i) .LE. m%index(1), r,&
            'Index member is greater than size')
    ENDDO
    !
    IF(dbg.GE.4) WRITE(*,*)r//' row_p extents',m%index(dbcsr_slot_row_p+1),&
         m%index(dbcsr_slot_row_p), SIZE(m%row_p)
    CALL bcsr_assert(dbg,m%index(dbcsr_slot_row_p+1)-m%index(dbcsr_slot_row_p)+1&
         .EQ.m%nblkrows_total+1, r,&
         'Size of row_p index inconsistent with number of rows')
    CALL bcsr_assert(dbg,SIZE(m%row_p).EQ.m%nblkrows_total+1, r,&
         'Size of row_p inconsistent with number of rows')
    !
    IF(dbg.GE.4) WRITE(*,*)r//' col_i extents',m%index(dbcsr_slot_col_i+1),&
         m%index(dbcsr_slot_col_i), SIZE(m%col_i)
    CALL bcsr_assert(dbg,m%index(dbcsr_slot_col_i+1)-m%index(dbcsr_slot_col_i)+1&
         .EQ.m%nblks, r,&
         'Size of col_i index inconsistent with number of blocks')
    CALL bcsr_assert(dbg,SIZE(m%col_i).EQ.m%nblks, r,&
         'Size of col inconsistent with number of blocks')
    !
    IF(dbg.GE.4) WRITE(*,*)r//' blk_p extents',m%index(dbcsr_slot_blk_p+1),&
         m%index(dbcsr_slot_blk_p), SIZE(m%blk_p)
    CALL bcsr_assert(dbg,m%index(dbcsr_slot_blk_p+1)-m%index(dbcsr_slot_blk_p)+1&
         .EQ.m%nblks, r,&
         'Size of blk_p index inconsistent with number of blocks')
    CALL bcsr_assert(dbg,SIZE(m%col_i).EQ.m%nblks, r,&
         'Size of blk_p inconsistent with number of blocks')
    !
!    IF(dbg) WRITE(*,*)r//' row_blk_begin',m%index(dbcsr_slot_row_blk_begin+1),&
!         m%index(dbcsr_slot_row_blk_begin), SIZE(m%row_blk_begin)
!    CALL bcsr_assert(dbg,m%index(dbcsr_slot_row_blk_begin+1)-&
!         m%index(dbcsr_slot_row_blk_begin)+1&
!         .EQ.m%nblkrows_total, r,&
!         'Size of row_blk_begin index inconsistent with number of rows')
!    CALL bcsr_assert(dbg,SIZE(m%row_blk_begin).EQ.m%nblkrows_total, r,&
!         'Size of row_blk_begin inconsistent with number of rows')
!    !
!    IF(dbg) WRITE(*,*)r//' row_blk_end',m%index(dbcsr_slot_row_blk_end+1),&
!         m%index(dbcsr_slot_row_blk_end), SIZE(m%row_blk_end)
!    CALL bcsr_assert(dbg,m%index(dbcsr_slot_row_blk_end+1)-&
!         m%index(dbcsr_slot_row_blk_end)+1&
!         .EQ.m%nblkrows_total, r,&
!         'Size of row_blk_end index inconsistent with number of rows')
!    CALL bcsr_assert(dbg,SIZE(m%row_blk_end).EQ.m%nblkrows_total, r,&
!         'Size of row_blk_end inconsistent with number of rows')
!    !
!    IF(dbg) WRITE(*,*)r//' col_blk_begin',m%index(dbcsr_slot_col_blk_begin+1),&
!         m%index(dbcsr_slot_col_blk_begin), SIZE(m%col_blk_begin)
!    CALL bcsr_assert(dbg,m%index(dbcsr_slot_col_blk_begin+1)-&
!         m%index(dbcsr_slot_col_blk_begin)+1&
!         .EQ.m%nblkcols_total, r,&
!         'Size of col_blk_begin index inconsistent with number of rows')
!    CALL bcsr_assert(dbg,SIZE(m%col_blk_begin).EQ.m%nblkcols_total, r,&
!         'Size of col_blk_begin inconsistent with number of columns')
!    !
!    IF(dbg) WRITE(*,*)r//' col_blk_end',m%index(dbcsr_slot_col_blk_end+1),&
!         m%index(dbcsr_slot_col_blk_end), SIZE(m%col_blk_end)
!    CALL bcsr_assert(dbg,m%index(dbcsr_slot_col_blk_end+1)-&
!         m%index(dbcsr_slot_col_blk_end)+1&
!         .EQ.m%nblkcols_total, r,&
!         'Size of col_blk_end index inconsistent with number of columns')
!    CALL bcsr_assert(dbg,SIZE(m%col_blk_end).EQ.m%nblkcols_total, r,&
!         'Size of col_blk_end inconsistent with number of columns')
    CALL bcsr_assert(dbg,SIZE(row_blk_size).EQ.m%nblkrows_total, r,&
         'Row block size array inconsistent with number of blocked rows')
    CALL bcsr_assert(dbg,SIZE(col_blk_size).EQ.m%nblkcols_total, r,&
         'Column block size array inconsistent with number of blocked columns')
    !
    IF(dbg.GE.4) THEN
       IF (ASSOCIATED (m%data))&
            WRITE(*,'(1X,A,I7,A,I7)')r//' nze=',m%nze,'data size',SIZE(m%data)
    ENDIF
    IF (ASSOCIATED (m%data))&
         CALL bcsr_assert(dbg,SIZE(m%data).GE.m%nze, r, 'Data storage too small')
    prev_br = 0
    DO br = 1, m%nblkrows_total
       CALL bcsr_assert(dbg,m%row_p(br).GE.0, r, 'row_p less than zero')
       IF(br.GT.1) CALL bcsr_assert(dbg,m%row_p(br).GE.m%row_p(prev_br),&
            r, 'row_p decreases')
       mb = row_blk_size(br)
       CALL bcsr_assert(dbg,mb.GT.0, r, 'Row blocked size is zero')
       DO blk = m%row_p(br)+1, m%row_p(br+1)
          CALL bcsr_assert(dbg,blk.GT.0, r, 'Block number is zero')
          CALL bcsr_assert(dbg,blk.LE.m%nblks, r, 'Block number too high')
          bc = m%col_i(blk)
          IF (dbg.GE.5) THEN
             WRITE(*,'(1X,A,I7,"(",I5,",",I5,")")')r//' block',blk,br,bc
          ENDIF
          CALL bcsr_assert(dbg,bc.GT.0, r, 'col_i is zero')
          CALL bcsr_assert(dbg,bc.LE.m%nblkcols_total, r, 'col_i too high')
          n = col_blk_size(bc)
          CALL bcsr_assert(dbg,n.GT.0, r, 'Column blocked size is zero')
          blk_p = m%blk_p(blk)
          mn = mb*n
          !CALL bcsr_assert(dbg,blk_p.GT.0, r, 'Block pointer is negative')
          CALL bcsr_assert(dbg,blk_p.LE.m%nze, r, 'Block pointer too large')
          CALL bcsr_assert(dbg,blk_p+mn-1.LE.m%nze, r, 'Block extends too far')
       ENDDO
       prev_br = br
    ENDDO
    IF (dbg.GE.5) THEN
       WRITE (*,'(1X,A,I7,A)') r//' size=',SIZE(m%row_p),' row_p='
       WRITE (*,2047) m%row_p(1:m%nblkrows_total)
       WRITE (*,'(1X,A)') r//' col_i='
       WRITE (*,2047) m%col_i(1:m%nblks)
       WRITE (*,'(1X,A)') r//' blk_p='
       WRITE (*,2047) m%blk_p(1:m%nblks)
    ENDIF
    IF (dbg.GE.3) THEN
       CALL mp_sum(n_have_blocks_local, dbcsr_mp_group (dbcsr_distribution_mp (&
            m%dist)))
       n_have_blocks_total = n_have_blocks_local
       sparsity_total = REAL(n_have_blocks_total)&
            / REAL(n_full_blocks_total)*100.0
       !WRITE(*,FMT='(30A,F5.1,A)')r//' Sparsity: ', sparsity_total,'%'
       WRITE(*,FMT='(1X,A,F5.1,A)')r//' Non-sparsity: ',&
            sparsity_total,'%'
    ENDIF
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_verify_matrix


! *****************************************************************************
!> \brief Glue for the real assert code
!> \param[in] l     Assertion
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
! *****************************************************************************
  SUBROUTINE bcsr_assert(level, l, routine, msg)
    INTEGER, INTENT(IN)                      :: level
    LOGICAL, INTENT(IN)                      :: l
    CHARACTER(len=*), INTENT(IN)             :: routine, msg

!   ---------------------------------------------------------------------------

    IF (level .GE. 2) THEN
       CALL cp_assert(l, cp_warning_level, cp_internal_error, routine, msg)
    ENDIF
  END SUBROUTINE bcsr_assert


! *****************************************************************************
!> \brief Makes a matrix canonical: sorts the indices and makes the data
!>        sequential.
!> \param[in,out] matrix      matrix in which to sort the data blocks
!> \param[in,out] error       cp2k error
!> \param[in] onlydata        (optional) only sort data
! *****************************************************************************
  SUBROUTINE dbcsr_order(matrix, error, onlydata)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: onlydata

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_order', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col_size, row, row_size, &
                                                stat, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sizes
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    REAL(KIND=dp), DIMENSION(:), POINTER     :: new_data

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    ! We should also check that the indices actually exist.
    CALL cp_assert(ASSOCIATED(matrix%data), cp_failure_level, cp_caller_error,&
         routineN, "Matrix has no data.", error)
    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)
    !
    IF (PRESENT (onlydata)) THEN
       IF (onlydata) THEN
          CALL dbcsr_sort_indices(matrix%nblks, matrix%row_p, matrix%col_i,&
               matrix%blk_p, error)
       ENDIF
    ENDIF
    ALLOCATE(sizes(matrix%nblks), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineN, "sizes", matrix%nblks)
!$OMP PARALLEL DO PRIVATE(row, blk, row_size, col_size)
    DO row = 1, matrix%nblkrows_total
       row_size = row_blk_size(row)
       DO blk = matrix%row_p(row)+1, matrix%row_p(row+1)
          col_size = col_blk_size(matrix%col_i(blk))
          sizes(blk) = row_size*col_size
       ENDDO
    ENDDO
!$OMP END PARALLEL DO
    IF (matrix%special_memory) THEN
       CALL mp_allocate(new_data, SIZE(matrix%data), stat=stat)
    ELSE
       ALLOCATE(new_data(SIZE(matrix%data)), stat=stat)
    ENDIF
    IF (stat /= 0) CALL stop_memory(routineN, "new_data", SIZE(matrix%data))
    CALL dbcsr_sort_data(matrix%blk_p, matrix%blk_p, sizes,&
         new_data, matrix%data, error)
    IF (matrix%special_memory) THEN
       CALL mp_deallocate(matrix%data)
    ELSE
       DEALLOCATE(matrix%data)
    ENDIF
    matrix%data => new_data
    DEALLOCATE(sizes)
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_order


! *****************************************************************************
!> \brief Sorts the rows & columns of a work matrix
!> \par Sorts the row and column indices so that the rows monotonically
!>      increase and the columns monotonically increase within each row.
!>      Passing the blk_p array rearranges the block pointers accordingly.
!> \param[in] n     number of blocks (elements) to sort
!> \param[in,out] row_i       row indices
!> \param[in,out] col_i       column indices
!> \param[in,out] blk_p       (optional) block pointers
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_sort_indices(n, row_i, col_i, blk_p, error)
    INTEGER, INTENT(IN)                      :: n
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: row_i, col_i
    INTEGER, DIMENSION(:), INTENT(INOUT), &
      OPTIONAL                               :: blk_p
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sort_indices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: buf

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    CALL cp_assert(SIZE(row_i).GE.n, cp_failure_level,&
         cp_caller_error, routineP, 'row_i too small', error=error)
    CALL cp_assert(SIZE(col_i).GE.n, cp_failure_level,&
         cp_caller_error, routineP, 'col_i too small', error=error)
    IF (PRESENT (blk_p)) CALL cp_assert(SIZE(blk_p).GE.n, cp_failure_level,&
         cp_caller_error, routineP, 'blk_p too small', error=error)
    IF (PRESENT (blk_p)) THEN
       ALLOCATE(buf(n), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP, 'buf', n)
       buf(1:n) = blk_p(1:n)
    ENDIF
    !> Create an ordering for both rows and columns while. If the blk_p must be
    !> rearranged, then the col_i array will be used as a permutation vector.
    row_i(:) = ISHFT(row_i(:), 16) + col_i
    IF (PRESENT (blk_p)) col_i(1:n) = (/ (i, i=1,n) /)
    !> Now do a nice quicksort.
    DBG ' sort, row_i=',row_i
    !col_i = (/ (i, i=1, n) /)
    CALL sort(row_i, n, col_i)
    DO i = 1, n
       IF (PRESENT (blk_p)) THEN
          blk_p(i) = buf(col_i(i))
       ENDIF
       col_i(i) = IAND(row_i(i), (ISHFT(1,16)-1))
       row_i(i) = ISHFT(row_i(i), -16)
    ENDDO
    DBG ' sort, row_i=',row_i
    DBG ' sort, col_i=',col_i
    IF(dbg.AND.bcsr_verbose.AND.PRESENT(blk_p))&
         WRITE(*,*)routineP//' sort, blk_p =',blk_p
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_sort_indices


  SUBROUTINE calc_blk_sizes(sizes, row_p, col_i, rbs, cbs)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: sizes
    INTEGER, DIMENSION(:), INTENT(IN)        :: row_p, col_i, rbs, cbs

    INTEGER                                  :: blk, col_size, nrows, row, &
                                                row_size

    nrows = SIZE (row_p)-1
!$OMP PARALLEL DO PRIVATE (row, blk, row_size, col_size)
    DO row = 1, nrows
       row_size = rbs(row)
       FORALL (blk = row_p(row)+1 : row_p(row+1))
          sizes(blk) = row_size * cbs(col_i(blk))
       END FORALL
    ENDDO
!$OMP END PARALLEL DO
  END SUBROUTINE calc_blk_sizes


! *****************************************************************************
!> \brief Sorts the data in a matrix so that the data blocks follow
!>        sequentially.
!> \param[out] blk_p          re-arranged block pointers reflecting the new
!>                            data order
!> \param[in] old_blk_p       current block pointers
!> \param[in] sizes           sizes of the data blocks
!> \param[out] dst            sorted data
!> \param[in] src             existing unordered data
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_sort_data(blk_p, old_blk_p, sizes, dst, src, error)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: blk_p
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_blk_p, sizes
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: dst
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: src
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sort_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nblks, timing_handle

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    nblks = SIZE(old_blk_p)
    DBG 'blk sizes', SIZE(blk_p), nblks
    DBG 'data sizes', SIZE(dst), SIZE(src)
    CALL cp_assert(SIZE(blk_p).GE.nblks, cp_failure_level,&
         cp_wrong_args_error, routineN, 'Destination blk_p too small.', error)
    CALL cp_assert(SIZE(dst).GE.SIZE(src), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Destination data space may be too small.',&
         error)
    DBGV 'old_blk_p',old_blk_p
    DBGV 'sizes', sizes
    IF (SIZE (old_blk_p) .GT. 0) THEN
       blk_p(1) = SIGN(1, old_blk_p(1))
       DO i = 2, nblks
          blk_p(i) = SIGN(ABS(blk_p(i-1)) + sizes(i-1), old_blk_p(i))
       ENDDO
       !$OMP PARALLEL DO PRIVATE(i)
       DO i = 1, nblks
          dst(ABS(blk_p(i)):ABS(blk_p(i))+sizes(i)-1) =&
               src(ABS(old_blk_p(i)):ABS(old_blk_p(i))+sizes(i)-1)
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF
    DBGV 'new blk_p =',blk_p
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_sort_data


! *****************************************************************************
!> \brief Determines the sign for a block pointer.
!> \param[in] n     block offset
!> \param[in] old   sign for the new block offset
!> \param[in] val   reverse the sign of the new block offset
! *****************************************************************************
  ELEMENTAL FUNCTION sgn (n, old, x) RESULT (val)
    INTEGER, INTENT(IN)                      :: n, old
    LOGICAL, INTENT(IN)                      :: x
    INTEGER                                  :: val

    IF (.NOT.x) THEN
       val = SIGN (n, old)
    ELSE
       val = -SIGN (n, old)
    ENDIF
  END FUNCTION sgn


! *****************************************************************************
!> \brief Sorts the data in a matrix so that the data blocks follow
!>        sequentially and does various transposing options.
!> \param[out] blk_p          re-arranged block pointers reflecting the new
!>                            data order
!> \param[in] old_blk_p       current block pointers
!> \param[in] rbs, cbs        sizes of the blocked rows and columns
!> \param[in] row_p, col_i    index
!> \param[out] dst            sorted data
!> \param[in] src             existing unordered data
!> \param[in] mark_transposed           mark data as transposed by negating the
!>                                      blk_p index entries
!> \param[in] transpose_blocks          transpose data blocks
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_copy_sort_data(blk_p, old_blk_p, row_p, col_i, rbs, cbs,&
       dst, src, mark_transposed, transpose_blocks, error)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: blk_p
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_blk_p, row_p, col_i, rbs, &
                                                cbs
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: dst
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: src
    LOGICAL, INTENT(IN), OPTIONAL            :: mark_transposed, &
                                                transpose_blocks
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_sort_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col_size, nblks, nrows, &
                                                nze, nze_prev, row, row_size
    LOGICAL                                  :: mark, trb

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    ! Analyze parameters
    mark = .FALSE.
    IF (PRESENT (mark_transposed)) mark = mark_transposed
    trb = .FALSE.
    IF (PRESENT (transpose_blocks)) trb = transpose_blocks
    !
    nblks = SIZE (old_blk_p)
    nrows = SIZE (row_p)-1
    CALL cp_assert(SIZE(blk_p).GE.nblks, cp_failure_level,&
         cp_wrong_args_error, routineN, 'Destination blk_p too small.', error)
    CALL cp_assert(SIZE(dst).GE.SIZE(src), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Destination data space may be too small.',&
         error)
    IF (nblks .GE. 1) &
         blk_p(1) = SGN(1, old_blk_p(1), mark)
    nze_prev = 0
    DO row = 1, nrows
       row_size = rbs (row)
       DO blk = row_p(row)+1, row_p(row+1)
          col_size = cbs(col_i(blk))
          nze = row_size * col_size
          IF (blk .GT. 1) THEN
             blk_p(blk) = SGN(ABS(blk_p(blk-1)) + nze_prev, old_blk_p(blk),&
                  mark)
          ENDIF
          IF (.NOT. trb) THEN
             dst(ABS(blk_p(blk)):ABS(blk_p(blk))+nze-1) =&
                  src(ABS(old_blk_p(blk)):ABS(old_blk_p(blk))+nze-1)
          ELSE
             CALL dbcsr_block_transpose (&
                  dst(ABS(blk_p(blk)):ABS(blk_p(blk))+nze-1),&
                  src(ABS(old_blk_p(blk)):ABS(old_blk_p(blk))+nze-1),&
                  row_size, col_size)
          ENDIF
          nze_prev = nze
       ENDDO
    ENDDO
  END SUBROUTINE dbcsr_copy_sort_data



! *****************************************************************************
!> \brief Prepares a matrix to be multipiled. Chooses the correct normalization
!>        routine based on the type of matrix and what normalization must take
!>        place
!> \param[in] source          input matrix
!> \param[in,out] normalized  image array of the normalized matrix
!> \param[in] target_image_dist          normalize to this image distribution
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_normalize(source, normalized, target_image_dist, error)
    TYPE(dbcsr_type), INTENT(IN)             :: source
    TYPE(dbcsr_array_type), INTENT(OUT)      :: normalized
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(IN)                             :: target_image_dist
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_normalize', &
      routineP = moduleN//':'//routineN

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    DBG 'Normalizing',source%matrix_type,source%name
    normalized%image_dist = target_image_dist
    SELECT CASE (source%matrix_type)
    CASE ('N')
       CALL dbcsr_normalize_normal_images(source, normalized, target_image_dist, error)
    CASE ('S')
       CALL dbcsr_normalize_symmetric(source, normalized, target_image_dist,&
            error)
    CASE ('T', 'A')
       CALL cp_assert(.FALSE., cp_failure_level,&
            cp_unimplemented_error_nr, routineP,&
            "Normalization of transposed or antisymmetric matrices not "//&
            "supported yet.", error)
    CASE DEFAULT
       CALL cp_assert(.FALSE., cp_failure_level,&
            cp_wrong_args_error, routineP, "Invalid matrix type.", error)
    END SELECT
    normalized%image_dist = target_image_dist
    DBG 'parent ums(1) name', normalized%mats(1)%name
  END SUBROUTINE dbcsr_normalize


! *****************************************************************************
!> \brief Prepares a matrix to be multiplied
!> \param[in] sm    input symmetric matrix
!> \param[in,out] ums         normalied matrices
!> \param[in] row_dist, col_dist        normalize to these distributions
!> \param[in] blacs2mpi       normalize to this processor grid
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_normalize_symmetric(sm, ums, target_imgdist, error)
    TYPE(dbcsr_type), INTENT(IN)             :: sm
    TYPE(dbcsr_array_type), INTENT(OUT), &
      TARGET                                 :: ums
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(IN)                             :: target_imgdist
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_normalize_symmetric', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: metalen = 2

    INTEGER :: blk, blk_l, blk_p, col, dst_p, i, mp_group, numproc, nze, &
      pcol, prow, row, row_rep, row_repeats, src_p, stored_col, stored_row, &
      timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blk_ps, blks, recv_meta, &
                                                send_meta
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: rd_disp, rm_disp, sd_disp, &
                                                sdp, sm_disp, smp, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: recv_count, send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist, &
                                                row_rep_dist
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    INTEGER, DIMENSION(dbcsr_meta_size)      :: new_meta
    LOGICAL                                  :: tr
    REAL(KIND=dp)                            :: tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: send_data
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: recv_data
    TYPE(dbcsr_distribution_obj)             :: target_dist
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_work_type), ALLOCATABLE, &
      DIMENSION(:)                           :: wms

!TYPE(dbcsr_work_type)                    :: wm
!TYPE(dbcsr_type), POINTER                :: um

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    tstart = m_walltime ()
    row_blk_size => array_data (sm%row_blk_size)
    col_blk_size => array_data (sm%col_blk_size)
    target_dist = target_imgdist%main
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    row_rep_dist => array_data (target_imgdist%row_image)
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numproc = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    DBGV 'row_dist', row_dist
    DBGV 'col_dist', col_dist
    DBG 'row, col maxval:', MAXVAL(row_dist), MAXVAL(col_dist)
    DBG 'blacs2mpi',blacs2mpi
    DBG 'blacs2mpi LB',LBOUND(blacs2mpi)
    DBG 'blacs2mpi UB',UBOUND(blacs2mpi)
    row_repeats = target_imgdist%nrow_images
    DBG 'row_repeats',row_repeats
    ALLOCATE (ums%mats(row_repeats))
    ALLOCATE (wms(row_repeats))
    !um => ums%mats(1)
    CALL cp_assert(SIZE(row_dist) .EQ. SIZE(col_dist), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Unequal row and column distributions for symmetric matrix.', error)
    CALL cp_assert(MAXVAL(row_dist).LE.UBOUND(blacs2mpi,1), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows', error)
    CALL cp_assert(MAXVAL(row_dist).EQ.UBOUND(blacs2mpi,1), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of row distribution not equal to processor rows', error)
    CALL cp_assert(MAXVAL(col_dist).LE.UBOUND(blacs2mpi,2), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols', error)
    CALL cp_assert(MAXVAL(col_dist).EQ.UBOUND(blacs2mpi,2), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of col distribution not equal to processor cols', error)
    ALLOCATE (send_count(2, row_repeats, 0:numproc-1))
    ALLOCATE (recv_count(2, row_repeats, 0:numproc-1))
    ALLOCATE (total_send_count(2, 0:numproc-1))
    ALLOCATE (total_recv_count(2, 0:numproc-1))
    ALLOCATE (sdp(row_repeats, 0:numproc-1))
    ALLOCATE (sd_disp(row_repeats, 0:numproc-1))
    ALLOCATE (smp(row_repeats, 0:numproc-1))
    ALLOCATE (sm_disp(row_repeats, 0:numproc-1))
    ALLOCATE (rd_disp(row_repeats, 0:numproc-1))
    ALLOCATE (rm_disp(row_repeats, 0:numproc-1))
    ! We know the upper bound of how big the new matrix will be, so we can
    ! easily allocate the new one.
    CALL dbcsr_pack_meta(sm, new_meta, error)
    new_meta(1:2) = 2 * new_meta(1:2)
    DBG 'A0'
    DO row_rep = 1, row_repeats
       CALL dbcsr_init (ums%mats(row_rep), error)
       CALL dbcsr_create(ums%mats(row_rep), 'norm '//sm%name, sm%dist,&
            'N', sm%row_blk_size, sm%col_blk_size,&
            0, 0, error=error, special=.TRUE.)
    ENDDO
    ! Count initial sizes for sending.
    DBG 'A'
    send_count(:,:,:) = 0
    DO row = 1, sm%nblkrows_total
       DO blk = sm%row_p(row)+1, sm%row_p(row+1)
          col = sm%col_i(blk)
          DO i = 1, 2
             IF (i .EQ. 1) THEN
                stored_row = row ; stored_col = col
             ELSE
                IF (row .EQ. col) CYCLE
                stored_row = col ; stored_col = row
             ENDIF
             ! Where do we send this block?
             prow = row_dist(stored_row)
             pcol = col_dist(stored_col)
             dst_p = blacs2mpi(prow, pcol)
             row_rep = row_rep_dist (stored_row)
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             send_count(1, row_rep, dst_p) = send_count(1, row_rep, dst_p) + 1
             send_count(2, row_rep, dst_p) = send_count(2, row_rep, dst_p) + nze
          ENDDO !i
       ENDDO ! col_i
    ENDDO ! row
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2*row_repeats, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
    ! Allocate data structures needed for data exchange.
    ALLOCATE (recv_data(SUM(recv_count(2, :, :))))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :, :))))
    ALLOCATE (send_data(SUM(send_count(2, :, :))))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :, :))))
    DBG 'send data len',SIZE (send_data)
    DBG 'send meta len',SIZE (send_meta)
    DBG 'recv data len',SIZE (recv_data)
    DBG 'recv meta len',SIZE (recv_meta)
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numproc-1
       total_send_count(1, dst_p) = SUM (send_count (1, :, dst_p))
       total_send_count(2, dst_p) = SUM (send_count (2, :, dst_p))
       total_recv_count(1, dst_p) = SUM (recv_count (1, :, dst_p))
       total_recv_count(2, dst_p) = SUM (recv_count (2, :, dst_p))
    ENDDO
    sd_disp = -1 ; sm_disp = -1
    rd_disp = -1 ; rm_disp = -1
    sd_disp(1,0) = 1 ; sm_disp(1,0) = 1
    rd_disp(1,0) = 1 ; rm_disp(1,0) = 1
    DO dst_p = 1, numproc-1
       sm_disp(1,dst_p) = sm_disp(1,dst_p-1)&
                        + metalen*total_send_count(1, dst_p-1)
       sd_disp(1,dst_p) = sd_disp(1,dst_p-1)&
                        + total_send_count(2, dst_p-1)
       rm_disp(1,dst_p) = rm_disp(1,dst_p-1)&
                        + metalen*total_recv_count(1, dst_p-1)
       rd_disp(1,dst_p) = rd_disp(1,dst_p-1)&
                        + total_recv_count(2, dst_p-1)
    ENDDO
    DO dst_p = 0, numproc-1
       DO row_rep = 2, row_repeats
          sm_disp (row_rep, dst_p) = sm_disp (row_rep-1, dst_p)&
                                   + metalen*send_count(1, row_rep-1, dst_p)
          sd_disp (row_rep, dst_p) = sd_disp (row_rep-1, dst_p)&
                                   + send_count(2, row_rep-1, dst_p)
          rm_disp (row_rep, dst_p) = rm_disp (row_rep-1, dst_p)&
                                   + metalen*recv_count(1, row_rep-1, dst_p)
          rd_disp (row_rep, dst_p) = rd_disp (row_rep-1, dst_p)&
                                   + recv_count(2, row_rep-1, dst_p)
       ENDDO
    ENDDO
    sdp = sd_disp ; smp = sm_disp
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    DBG 'D'
    DO row = 1, sm%nblkrows_total
       DO blk = sm%row_p(row)+1, sm%row_p(row+1)
          col = sm%col_i(blk)
          blk_p = sm%blk_p(blk)
          DO i = 1, 2
             IF (i .EQ. 1) THEN
                stored_row = row ; stored_col = col; tr = .FALSE.
             ELSE                                                
                IF (row .EQ. col) CYCLE
                stored_row = col ; stored_col = row; tr = .TRUE.
             ENDIF
             ! Where do we send this block?
             prow = row_dist(stored_row)
             pcol = col_dist(stored_col)
             dst_p = blacs2mpi(prow, pcol)
             row_rep = row_rep_dist(stored_row)
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             send_meta(smp(row_rep, dst_p)) = stored_row
             send_meta(smp(row_rep, dst_p)+1) = stored_col
             DBG stored_row, stored_col, blk_p, tr, nze
             IF (.NOT. tr) THEN
                send_data(sdp(row_rep, dst_p):sdp(row_rep, dst_p)+nze-1) =&
                     sm%data(blk_p:blk_p+nze-1)
             ELSE
                send_meta(smp(row_rep, dst_p)) = -stored_row
                send_data(sdp(row_rep, dst_p):sdp(row_rep, dst_p)+nze-1) =&
                     sm%data(blk_p:blk_p+nze-1)
             ENDIF
             smp(row_rep, dst_p) = smp(row_rep, dst_p) + metalen
             sdp(row_rep, dst_p) = sdp(row_rep, dst_p) + nze
          ENDDO ! i
       ENDDO ! col_i
    ENDDO ! row
    ! Exchange the data and metadata structures.
    DBG 'E'
    CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(1,:)-1,&
         recv_data(:), total_recv_count(2,:), rd_disp(1,:)-1, mp_group)
    DBG 'F'
    DBG 'send_meta',send_meta
    CALL mp_alltoall(send_meta(:), metalen*total_send_count(1,:), sm_disp(1,:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(1,:)-1, mp_group)
    DBG 'recv_meta',recv_meta
    DBG 'G'
    ! Now fill in the data.
    DO row_rep = 1, row_repeats
       CALL dbcsr_work_create(ums%mats(row_rep),&
            SUM(recv_count(1,row_rep,:)),&
            SUM(recv_count(2,row_rep,:)), error=error)
       wms(row_rep) = ums%mats(row_rep)%w
    ENDDO
    !DBG 'send disps data',sd_disp
    !DBG 'send disps meta',sm_disp
    !DBG 'recv disps data',rd_disp
    !DBG 'recv disps meta',rm_disp
    ALLOCATE (blk_ps(row_repeats))
    ALLOCATE (blks (row_repeats))
    !blk_p = 1
    !blk = 1
    blk_ps(:) = 1
    blks(:) = 1
    DO src_p = 0, numproc-1
       !wm%data(blk_p : blk_p+recv_count(2,src_p)-1) =&
       !     recv_data(rd_disp(src_p) :&
       !               rd_disp(src_p)+recv_count(2,src_p)-1)
       DO row_rep = 1, row_repeats
          nze = recv_count(2, row_rep, src_p)
          wms(row_rep)%data(blk_ps(row_rep):blk_ps(row_rep)+nze-1) =&
               recv_data(rd_disp(row_rep, src_p):rd_disp(row_rep, src_p)+nze-1)
          DBG 'processing from',src_p,'counts',send_count(1, :, src_p)
          DO blk_l = 1, recv_count(1, row_rep, src_p)
             stored_row = recv_meta(rm_disp(row_rep,src_p)+metalen*(blk_l-1))
             stored_col = recv_meta(rm_disp(row_rep,src_p)+metalen*(blk_l-1)+1)
             !DBG 'blk: p',blk,src_p,'; row, col',stored_row,stored_col,&
             !     '; blk_p',blk_p,'; row_rep, imgdist',row_rep,&
             !     row_rep_dist (ABS (stored_row))
             DBG 'blk: p',blks(row_rep),src_p,'; row, col',&
                  stored_row,stored_col,&
                  '; blk_ps',blk_ps,'; row_rep, imgdist',row_rep,&
                  row_rep_dist (ABS (stored_row))
             !row_rep = row_rep_dist (ABS(stored_row))
             wms(row_rep)%row_i(blks(row_rep)) = ABS(stored_row)
             wms(row_rep)%col_i(blks(row_rep)) = stored_col
             wms(row_rep)%blk_p(blks(row_rep)) =&
                  SIGN(blk_ps(row_rep), stored_row)
             nze = row_blk_size(ABS(stored_row))&
                  * col_blk_size(stored_col)
             !blk_p = blk_p + nze
             !blk = blk + 1
             blk_ps(row_rep) = blk_ps(row_rep) + nze
             blks(row_rep) = blks(row_rep) + 1
          ENDDO
       ENDDO
    ENDDO
    !
    DO row_rep = 1, row_repeats
       wms(row_rep)%lastblk = blks(row_rep) - 1
       wms(row_rep)%datasize = blk_ps(row_rep) - 1
       DBG 'Finalizing normalization'
       ums%mats(row_rep)%w = wms(row_rep)
       CALL dbcsr_finalize(ums%mats(row_rep),&
            resort=.TRUE., reshuffle=.TRUE.,&
            error=error)
    ENDDO
    tstop = m_walltime ()
    DBG 'time:',tstop-tstart
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_data)
    DEALLOCATE(send_meta)
    DEALLOCATE (wms)
    DBG 'ums(1) name',ums%mats(1)%name
    DBG 'Done desymmetrizing!'
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_normalize_symmetric


! *****************************************************************************
!> \brief Prepares a matrix to be multiplied
!> \param[in] sm    input matrix
!> \param[in,out] um          normalied matrix
!> \param[in] row_dist, col_dist        normalize to these distributions
!> \param[in] blacs2mpi       normalize to this processor grid
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_normalize_normal(sm, ums, target_imgdist, error)
    TYPE(dbcsr_type), INTENT(IN)             :: sm
    TYPE(dbcsr_array_type), INTENT(OUT), &
      TARGET                                 :: ums
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(IN)                             :: target_imgdist
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_normalize_normal', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, blk_l, blk_p, col, dst_p, mp_group, numproc, nze, pcol, &
      prow, row, src_p, stored_col, stored_row, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                sd_disp, sdp, send_meta, &
                                                sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    INTEGER, DIMENSION(dbcsr_meta_size)      :: new_meta
    REAL(KIND=dp)                            :: tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: send_data
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: recv_data
    TYPE(dbcsr_distribution_obj)             :: target_dist
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_type), POINTER                :: um
    TYPE(dbcsr_work_type)                    :: wm

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    tstart = m_walltime ()
    row_blk_size => array_data (sm%row_blk_size)
    col_blk_size => array_data (sm%col_blk_size)
    target_dist = target_imgdist%main
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numproc = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    DBG 'blacs2mpi',blacs2mpi
    ALLOCATE (send_count(2, 0:numproc-1))
    ALLOCATE (recv_count(2, 0:numproc-1))
    ALLOCATE (sdp(0:numproc-1))
    ALLOCATE (sd_disp(0:numproc-1))
    ALLOCATE (smp(0:numproc-1))
    ALLOCATE (sm_disp(0:numproc-1))
    ALLOCATE (rd_disp(0:numproc-1))
    ALLOCATE (rm_disp(0:numproc-1))
    ALLOCATE (ums%mats(1))
    um => ums%mats(1)
    ! We know the upper bound of how big the new matrix will be, so we can
    ! easily allocate the new one.
    CALL dbcsr_pack_meta(sm, new_meta, error)
    new_meta(1:2) = 2 * new_meta(1:2)
    DBG 'A0'
    CALL dbcsr_init (ums%mats(1), error)
    CALL dbcsr_create(ums%mats(1), 'norm '//sm%name, sm%dist, 'N',&
         sm%row_blk_size, sm%col_blk_size, 0, 0, error=error)
    um => ums%mats(1)
    ! Count initial sizes for sending.
    DBG 'A'
    send_count(:,:) = 0
    DO row = 1, sm%nblkrows_total
       DO blk = sm%row_p(row)+1, sm%row_p(row+1)
          col = sm%col_i(blk)
          stored_row = row ; stored_col = col
          ! Where do we send this block?
          prow = row_dist(stored_row)
          pcol = col_dist(stored_col)
          dst_p = blacs2mpi(prow, pcol)
          nze = row_blk_size(stored_row) * col_blk_size(stored_col)
          send_count(1, dst_p) = send_count(1, dst_p) + 1
          send_count(2, dst_p) = send_count(2, dst_p) + nze
       ENDDO ! col_i
    ENDDO ! row
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
    ! Allocate data structures needed for data exchange.
    ALLOCATE (recv_data(SUM(recv_count(2, :))))
    ALLOCATE (recv_meta(2*SUM(recv_count(1, :))))
    ALLOCATE (send_data(SUM(send_count(2, :))))
    ALLOCATE (send_meta(2*SUM(send_count(1, :))))
    ! Fill in the meta data structures and copy the data.
    sd_disp(0) = 1 ; sm_disp(0) = 1
    rd_disp(0) = 1 ; rm_disp(0) = 1
    DO dst_p = 1, numproc-1
       sm_disp(dst_p) = sm_disp(dst_p-1) + 2*send_count(1, dst_p-1)
       sd_disp(dst_p) = sd_disp(dst_p-1) + send_count(2, dst_p-1)
       rm_disp(dst_p) = rm_disp(dst_p-1) + 2*recv_count(1, dst_p-1)
       rd_disp(dst_p) = rd_disp(dst_p-1) + recv_count(2, dst_p-1)
    ENDDO
    sdp = sd_disp ; smp = sm_disp
    DBG 'D'
    DO row = 1, sm%nblkrows_total
       DO blk = sm%row_p(row)+1, sm%row_p(row+1)
          col = sm%col_i(blk)
          blk_p = sm%blk_p(blk)
          stored_row = row ; stored_col = col
          ! Where do we send this block?
          prow = row_dist(stored_row)
          pcol = col_dist(stored_col)
          dst_p = blacs2mpi(prow, pcol)
          nze = row_blk_size(stored_row) * col_blk_size(stored_col)
          send_meta(smp(dst_p)) = stored_row
          send_meta(smp(dst_p)+1) = stored_col
          smp(dst_p) = smp(dst_p) + 2
          send_data(sdp(dst_p):sdp(dst_p)+nze-1) = sm%data(blk_p:blk_p+nze-1)
          sdp(dst_p) = sdp(dst_p) + nze
       ENDDO ! col_i
    ENDDO ! row
    ! Exchange the data and metadata structures.
    DBG 'E'
    CALL mp_alltoall(send_data(:), send_count(2,:), sd_disp(:)-1,&
         recv_data(:), recv_count(2,:), rd_disp(:)-1, mp_group)
    DBG 'F'
    DBG 'send_meta',send_meta
    CALL mp_alltoall(send_meta(:), 2*send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), 2*recv_count(1,:), rm_disp(:)-1, mp_group)
    DBG 'recv_meta',recv_meta
    DBG 'G'
    ! Now fill in the data.
    CALL dbcsr_work_create(um, SUM(recv_count(1,:)),&
         SUM(recv_count(2,:)), error=error)
    wm = um%w
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    blk_p = 1
    blk = 1
    DO src_p = 0, numproc-1
       wm%data(blk_p : blk_p+recv_count(2,src_p)-1) =&
            recv_data(rd_disp(src_p) :&
                      rd_disp(src_p)+recv_count(2,src_p)-1)
       DBG 'processing from',src_p,'count',send_count(1, src_p)
       DO blk_l = 1, recv_count(1, src_p)
          stored_row = recv_meta(rm_disp(src_p)+2*(blk_l-1))
          stored_col = recv_meta(rm_disp(src_p)+2*(blk_l-1)+1)
          DBG 'blk: p',blk,src_p,'; row, col',stored_row,stored_col,&
               '; blk_p',blk_p
          wm%row_i(blk) = stored_row
          wm%col_i(blk) = stored_col
          wm%blk_p(blk) = blk_p
          nze = row_blk_size(stored_row) * col_blk_size(stored_col)
          blk_p = blk_p + nze
          blk = blk + 1
       ENDDO
    ENDDO
    !
    wm%lastblk = blk-1
    wm%datasize = blk_p - 1
    DBG 'Finalizing normalization'
    um%w = wm
    CALL dbcsr_finalize(um, resort=.TRUE., reshuffle=.TRUE., error=error)
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_data)
    DEALLOCATE(send_meta)
    tstop = m_walltime ()
    DBG 'time:',tstop-tstart
    DBG 'Done normalizing.'
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_normalize_normal

! *****************************************************************************
!> \brief Prepares a matrix to be multiplied
!> \param[in] sm    input matrix
!> \param[in,out] um          normalied matrix
!> \param[in] row_dist, col_dist        normalize to these distributions
!> \param[in] blacs2mpi       normalize to this processor grid
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_normalize_normal_images(sm, ums, target_imgdist, error)
    TYPE(dbcsr_type), INTENT(IN)             :: sm
    TYPE(dbcsr_array_type), INTENT(OUT), &
      TARGET                                 :: ums
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(IN)                             :: target_imgdist
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_normalize_normal_images', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: metalen = 2

    INTEGER :: blk, blk_l, blk_p, col, dst_p, i, mp_group, numproc, nze, &
      pcol, prow, row, row_rep, row_repeats, src_p, stored_col, stored_row, &
      timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blk_ps, blks, recv_meta, &
                                                send_meta
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: rd_disp, rm_disp, sd_disp, &
                                                sdp, sm_disp, smp, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: recv_count, send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist, &
                                                row_rep_dist
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    INTEGER, DIMENSION(dbcsr_meta_size)      :: new_meta
    LOGICAL                                  :: tr
    REAL(KIND=dp)                            :: tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: send_data
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: recv_data
    TYPE(dbcsr_distribution_obj)             :: target_dist
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_work_type), ALLOCATABLE, &
      DIMENSION(:)                           :: wms

!TYPE(dbcsr_work_type)                    :: wm
!TYPE(dbcsr_type), POINTER                :: um

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    tstart = m_walltime ()
    row_blk_size => array_data (sm%row_blk_size)
    col_blk_size => array_data (sm%col_blk_size)
    target_dist = target_imgdist%main
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    row_rep_dist => array_data (target_imgdist%row_image)
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numproc = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    DBGV 'row_dist', row_dist
    DBGV 'col_dist', col_dist
    DBG 'row, col maxval:', MAXVAL(row_dist), MAXVAL(col_dist)
    DBG 'blacs2mpi',blacs2mpi
    DBG 'blacs2mpi LB',LBOUND(blacs2mpi)
    DBG 'blacs2mpi UB',UBOUND(blacs2mpi)
    row_repeats = target_imgdist%nrow_images
    DBG 'row_repeats',row_repeats
    ALLOCATE (ums%mats(row_repeats))
    ALLOCATE (wms(row_repeats))
    !um => ums%mats(1)
    CALL cp_assert(MAXVAL(row_dist).LE.UBOUND(blacs2mpi,1), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows', error)
    CALL cp_assert(MAXVAL(row_dist).EQ.UBOUND(blacs2mpi,1), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of row distribution not equal to processor rows', error)
    CALL cp_assert(MAXVAL(col_dist).LE.UBOUND(blacs2mpi,2), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols', error)
    CALL cp_assert(MAXVAL(col_dist).EQ.UBOUND(blacs2mpi,2), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of col distribution not equal to processor cols', error)
    ALLOCATE (send_count(2, row_repeats, 0:numproc-1))
    ALLOCATE (recv_count(2, row_repeats, 0:numproc-1))
    ALLOCATE (total_send_count(2, 0:numproc-1))
    ALLOCATE (total_recv_count(2, 0:numproc-1))
    ALLOCATE (sdp(row_repeats, 0:numproc-1))
    ALLOCATE (sd_disp(row_repeats, 0:numproc-1))
    ALLOCATE (smp(row_repeats, 0:numproc-1))
    ALLOCATE (sm_disp(row_repeats, 0:numproc-1))
    ALLOCATE (rd_disp(row_repeats, 0:numproc-1))
    ALLOCATE (rm_disp(row_repeats, 0:numproc-1))
    ! We know the upper bound of how big the new matrix will be, so we can
    ! easily allocate the new one.
    CALL dbcsr_pack_meta(sm, new_meta, error)
    new_meta(1:2) = 2 * new_meta(1:2)
    DBG 'A0'
    DO row_rep = 1, row_repeats
       CALL dbcsr_init (ums%mats(row_rep), error)
       CALL dbcsr_create(ums%mats(row_rep), 'norm '//sm%name, sm%dist,&
            'N', sm%row_blk_size, sm%col_blk_size,&
            0, 0, error=error)
    ENDDO
    ! Count initial sizes for sending.
    DBG 'A'
    send_count(:,:,:) = 0
    DO row = 1, sm%nblkrows_total
       DO blk = sm%row_p(row)+1, sm%row_p(row+1)
          col = sm%col_i(blk)
          DO i = 1, 1
             IF (i .EQ. 1) THEN
                stored_row = row ; stored_col = col
             ELSE
                IF (row .EQ. col) CYCLE
                stored_row = col ; stored_col = row
             ENDIF
             ! Where do we send this block?
             prow = row_dist(stored_row)
             pcol = col_dist(stored_col)
             dst_p = blacs2mpi(prow, pcol)
             row_rep = row_rep_dist (stored_row)
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             send_count(1, row_rep, dst_p) = send_count(1, row_rep, dst_p) + 1
             send_count(2, row_rep, dst_p) = send_count(2, row_rep, dst_p) + nze
          ENDDO !i
       ENDDO ! col_i
    ENDDO ! row
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2*row_repeats, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
    ! Allocate data structures needed for data exchange.
    ALLOCATE (recv_data(SUM(recv_count(2, :, :))))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :, :))))
    ALLOCATE (send_data(SUM(send_count(2, :, :))))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :, :))))
    DBG 'send data len',SIZE (send_data)
    DBG 'send meta len',SIZE (send_meta)
    DBG 'recv data len',SIZE (recv_data)
    DBG 'recv meta len',SIZE (recv_meta)
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numproc-1
       total_send_count(1, dst_p) = SUM (send_count (1, :, dst_p))
       total_send_count(2, dst_p) = SUM (send_count (2, :, dst_p))
       total_recv_count(1, dst_p) = SUM (recv_count (1, :, dst_p))
       total_recv_count(2, dst_p) = SUM (recv_count (2, :, dst_p))
    ENDDO
    sd_disp = -1 ; sm_disp = -1
    rd_disp = -1 ; rm_disp = -1
    sd_disp(1,0) = 1 ; sm_disp(1,0) = 1
    rd_disp(1,0) = 1 ; rm_disp(1,0) = 1
    DO dst_p = 1, numproc-1
       sm_disp(1,dst_p) = sm_disp(1,dst_p-1)&
                        + metalen*total_send_count(1, dst_p-1)
       sd_disp(1,dst_p) = sd_disp(1,dst_p-1)&
                        + total_send_count(2, dst_p-1)
       rm_disp(1,dst_p) = rm_disp(1,dst_p-1)&
                        + metalen*total_recv_count(1, dst_p-1)
       rd_disp(1,dst_p) = rd_disp(1,dst_p-1)&
                        + total_recv_count(2, dst_p-1)
    ENDDO
    DO dst_p = 0, numproc-1
       DO row_rep = 2, row_repeats
          sm_disp (row_rep, dst_p) = sm_disp (row_rep-1, dst_p)&
                                   + metalen*send_count(1, row_rep-1, dst_p)
          sd_disp (row_rep, dst_p) = sd_disp (row_rep-1, dst_p)&
                                   + send_count(2, row_rep-1, dst_p)
          rm_disp (row_rep, dst_p) = rm_disp (row_rep-1, dst_p)&
                                   + metalen*recv_count(1, row_rep-1, dst_p)
          rd_disp (row_rep, dst_p) = rd_disp (row_rep-1, dst_p)&
                                   + recv_count(2, row_rep-1, dst_p)
       ENDDO
    ENDDO
    sdp = sd_disp ; smp = sm_disp
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    DBG 'D'
    DO row = 1, sm%nblkrows_total
       DO blk = sm%row_p(row)+1, sm%row_p(row+1)
          col = sm%col_i(blk)
          blk_p = sm%blk_p(blk)
          DO i = 1, 1
             IF (i .EQ. 1) THEN
                stored_row = row ; stored_col = col; tr = .FALSE.
             ELSE                                                
                IF (row .EQ. col) CYCLE
                stored_row = col ; stored_col = row; tr = .TRUE.
             ENDIF
             ! Where do we send this block?
             prow = row_dist(stored_row)
             pcol = col_dist(stored_col)
             dst_p = blacs2mpi(prow, pcol)
             row_rep = row_rep_dist(stored_row)
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             send_meta(smp(row_rep, dst_p)) = stored_row
             send_meta(smp(row_rep, dst_p)+1) = stored_col
             DBG stored_row, stored_col, blk_p, tr, nze
             IF (.NOT. tr) THEN
                send_data(sdp(row_rep, dst_p):sdp(row_rep, dst_p)+nze-1) =&
                     sm%data(blk_p:blk_p+nze-1)
             ELSE
                send_meta(smp(row_rep, dst_p)) = -stored_row
                send_data(sdp(row_rep, dst_p):sdp(row_rep, dst_p)+nze-1) =&
                     sm%data(blk_p:blk_p+nze-1)
             ENDIF
             smp(row_rep, dst_p) = smp(row_rep, dst_p) + metalen
             sdp(row_rep, dst_p) = sdp(row_rep, dst_p) + nze
          ENDDO ! i
       ENDDO ! col_i
    ENDDO ! row
    ! Exchange the data and metadata structures.
    DBG 'E'
    CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(1,:)-1,&
         recv_data(:), total_recv_count(2,:), rd_disp(1,:)-1, mp_group)
    DBG 'F'
    DBG 'send_meta',send_meta
    CALL mp_alltoall(send_meta(:), metalen*total_send_count(1,:), sm_disp(1,:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(1,:)-1, mp_group)
    DBG 'recv_meta',recv_meta
    DBG 'G'
    ! Now fill in the data.
    DO row_rep = 1, row_repeats
       CALL dbcsr_work_create(ums%mats(row_rep),&
            SUM(recv_count(1,row_rep,:)),&
            SUM(recv_count(2,row_rep,:)), error=error)
       wms(row_rep) = ums%mats(row_rep)%w
    ENDDO
    !DBG 'send disps data',sd_disp
    !DBG 'send disps meta',sm_disp
    !DBG 'recv disps data',rd_disp
    !DBG 'recv disps meta',rm_disp
    ALLOCATE (blk_ps(row_repeats))
    ALLOCATE (blks (row_repeats))
    !blk_p = 1
    !blk = 1
    blk_ps(:) = 1
    blks(:) = 1
    DO src_p = 0, numproc-1
       !wm%data(blk_p : blk_p+recv_count(2,src_p)-1) =&
       !     recv_data(rd_disp(src_p) :&
       !               rd_disp(src_p)+recv_count(2,src_p)-1)
       DO row_rep = 1, row_repeats
          nze = recv_count(2, row_rep, src_p)
          wms(row_rep)%data(blk_ps(row_rep):blk_ps(row_rep)+nze-1) =&
               recv_data(rd_disp(row_rep, src_p):rd_disp(row_rep, src_p)+nze-1)
          DBG 'processing from',src_p,'counts',send_count(1, :, src_p)
          DO blk_l = 1, recv_count(1, row_rep, src_p)
             stored_row = recv_meta(rm_disp(row_rep,src_p)+metalen*(blk_l-1))
             stored_col = recv_meta(rm_disp(row_rep,src_p)+metalen*(blk_l-1)+1)
             !DBG 'blk: p',blk,src_p,'; row, col',stored_row,stored_col,&
             !     '; blk_p',blk_p,'; row_rep, imgdist',row_rep,&
             !     row_rep_dist (ABS (stored_row))
             DBG 'blk: p',blks(row_rep),src_p,'; row, col',&
                  stored_row,stored_col,&
                  '; blk_ps',blk_ps,'; row_rep, imgdist',row_rep,&
                  row_rep_dist (ABS (stored_row))
             !row_rep = row_rep_dist (ABS(stored_row))
             wms(row_rep)%row_i(blks(row_rep)) = ABS(stored_row)
             wms(row_rep)%col_i(blks(row_rep)) = stored_col
             wms(row_rep)%blk_p(blks(row_rep)) =&
                  SIGN(blk_ps(row_rep), stored_row)
             nze = row_blk_size(ABS(stored_row))&
                  * col_blk_size(stored_col)
             !blk_p = blk_p + nze
             !blk = blk + 1
             blk_ps(row_rep) = blk_ps(row_rep) + nze
             blks(row_rep) = blks(row_rep) + 1
          ENDDO
       ENDDO
    ENDDO
    !
    DO row_rep = 1, row_repeats
       wms(row_rep)%lastblk = blks(row_rep) - 1
       wms(row_rep)%datasize = blk_ps(row_rep) - 1
       DBG 'Finalizing normalization'
       ums%mats(row_rep)%w = wms(row_rep)
       CALL dbcsr_finalize(ums%mats(row_rep),&
            resort=.TRUE., reshuffle=.TRUE.,&
            error=error)
    ENDDO
    tstop = m_walltime ()
    DBG 'time:',tstop-tstart
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_data)
    DEALLOCATE(send_meta)
    DEALLOCATE (wms)
    DBG 'ums(1) name',ums%mats(1)%name
    DBG 'Done desymmetrizing!'
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_normalize_normal_images


! *****************************************************************************
!> \brief Copies an array into another array, ensuring that the target array
!>        is large enough.
!> \par Makes use of FORTRAN's array bounds for array sections.
!> \param to        array to copy to (target array)
!> \param[in] from  array to copy from (source array)
!> \param[in] offset          starting position in the target array
!> \param[in] factor          factor for size increase
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE addto_array_i(to, from, offset, factor, special, error)
    INTEGER, DIMENSION(:), POINTER           :: to
    INTEGER, DIMENSION(:), INTENT(IN)        :: from
    INTEGER, INTENT(IN), OPTIONAL            :: offset
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: special
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'addto_array_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: last, lb_t, s_f, s_t, start, &
                                                timing_handle, ub_t
    LOGICAL                                  :: spec
    REAL                                     :: f

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    s_f = SIZE(from, 1)
    ub_t = UBOUND(to, 1)
    lb_t = LBOUND(to, 1)
    s_t = SIZE(to)
    start = lb_t
    IF (PRESENT (offset)) start = offset
    ! assert offset > lb_t
    last = start+s_f-1
    IF(dbg) WRITE(*,*)routineP//' target array bounds',lb_t,':',ub_t,&
         'starting at',start,'last',last,'source size', s_f
    IF (ub_t .LT. last) THEN
       IF (PRESENT (factor)) THEN
          f = factor
       ELSE
          f = 1.4142136
       ENDIF
       ub_t = MAX(start + INT(s_t * f), last)
       IF(dbg) WRITE(*,*)routineP//' expanding array from :',UBOUND(to),&
            'to',ub_t
       !CALL reallocate(to, lb_t, ub_t)
       CALL ensure_array_size(to, lb_t, ub_t, factor=1.0, nocopy=.FALSE.,&
            special=spec, error=error)
    ENDIF
    to(start:last) = from(:)
    CALL timestop(timing_handle)
  END SUBROUTINE addto_array_i


! *****************************************************************************
!> \brief Copies an array into another array, ensuring that the target array
!>        large enough.
!> \par Makes use of FORTRAN's array bounds for array sections.
!> \param to        array to copy to
!> \param[in] from  array to copy from
!> \param[in] offset          starting position in the target array
!> \param[in] factor          factor for size increase
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE addto_array_d(to, from, offset, factor, special, error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: to
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: from
    INTEGER, INTENT(IN), OPTIONAL            :: offset
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: special
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'addto_array_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: last, lb_t, s_f, s_t, start, &
                                                timing_handle, ub_t
    LOGICAL                                  :: dbg, spec
    REAL                                     :: f

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    CALL timeset(routineN, timing_handle)
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    s_f = SIZE(from, 1)
    IF (dbg) THEN
       WRITE(*,*)routineP//' from size',SIZE(from),';', SIZE(from,1)
       WRITE(*,*)routineP//' from lbound',LBOUND(from),';', LBOUND(from,1)
       WRITE(*,*)routineP//' from ubound',UBOUND(from),';', UBOUND(from,1)
    ENDIF
    ub_t = UBOUND(to, 1)
    lb_t = LBOUND(to, 1)
    s_t = SIZE(to)
    start = lb_t
    IF (PRESENT (offset)) start = offset
    ! assert offset > lb_t
    last = start+s_f-1
    IF(dbg) WRITE(*,*)routineP//' target array bounds',lb_t,':',ub_t,&
         'starting at',start,'last',last,'source size', s_f
    IF (ub_t .LT. last) THEN
       IF (PRESENT (factor)) THEN
          f = factor
       ELSE
          f = 1.4142136
       ENDIF
       !ub_t = MAX(lb_t + INT(s_t * f), lb_t+s_f-1)
       ub_t = MAX(start + INT(s_t * f), last)
       IF(dbg) WRITE(*,*)routineP//' expanding array from :',UBOUND(to),&
            'to',ub_t
       !CALL reallocate(to, lb_t, ub_t)
       CALL ensure_array_size(to, lb_t, ub_t, factor=1.0, nocopy=.FALSE.,&
            special=spec, error=error)
    ENDIF
    to(start:last) = from(:)
    CALL timestop(timing_handle)
  END SUBROUTINE addto_array_d


! *****************************************************************************
!> \brief Ensures that an array is appropriately large.
!> \param[in,out] array       array to verify and possibly resize
!> \param[in] lb    (optional) desired array lower bound
!> \param[in] ub    desired array upper bound
!> \param[in] factor          (optional) factor by which to exagerrate
!>                            enlargements
!> \param[in] nocopy          (optional) copy array on enlargement; default
!>                            is to copy
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE ensure_array_size_i(array, lb, ub, factor, nocopy, special, error)
    INTEGER, DIMENSION(:), POINTER           :: array
    INTEGER, INTENT(IN), OPTIONAL            :: lb
    INTEGER, INTENT(IN)                      :: ub
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, special
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ensure_array_size_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lb_new, lb_orig, stat, &
                                                timing_handle, ub_new, ub_orig
    INTEGER, DIMENSION(:), POINTER           :: newarray
    LOGICAL                                  :: dbg, docopy = .TRUE., &
                                                failure = .TRUE., info, spec

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    info = .FALSE.
    CALL dbcsr_set_debug(dbg, info=info)
    CALL timeset(routineN, timing_handle)
    docopy = .TRUE.
    IF (PRESENT (nocopy)) docopy = .NOT. nocopy
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    lb_new = 1
    IF (PRESENT (lb)) lb_new = lb
    !> Creates a new array if it doesn't yet exist.
    IF (.NOT.ASSOCIATED(array)) THEN
       failure = .FALSE.
       !CPPrecondition(.NOT.docopy, cp_warning_level, routineP, error, failure)
       failure = .TRUE.
       IF (spec) THEN
          CPPostcondition(lb_new.EQ.1, cp_warning_level, routineP, error, failure)
          CALL mp_allocate(array, ub-lb_new+1, stat=stat)
       ELSE
          ALLOCATE(array(lb_new:ub), stat=stat)
       ENDIF
       IF (stat /= 0) CALL stop_memory(routineP, 'array', ub-lb_new+1)
       CALL timestop(timing_handle)
       RETURN
    ENDIF
    lb_orig = LBOUND(array,1)
    ub_orig = UBOUND(array,1)
    IF(dbg) WRITE(*,*)routineP//' Current bounds are',lb_orig,':',ub_orig
    IF (lb_orig.LE.lb_new .AND. ub_orig.GE.ub) THEN
       CALL timestop(timing_handle)
       RETURN
    ENDIF
    IF (lb_orig.GT.lb_new) THEN
       IF (PRESENT(factor)) THEN
          lb_new = MIN(lb_orig, lb_new - NINT(REAL(lb_orig - lb_new)*factor))
       ENDIF
    ENDIF
    IF (ub_orig.LT.ub) THEN
       IF (PRESENT(factor)) THEN
          ub_new = MAX(ub_orig, ub + NINT(REAL(ub-ub_orig)*factor))
       ELSE
          ub_new = ub
       ENDIF
    ELSE
       ub_new = ub
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Resizing to bounds',lb_new,':',ub_new
    IF(.NOT.docopy) THEN
       IF (spec) THEN
          CALL mp_deallocate(array)
       ELSE
          DEALLOCATE(array)
       ENDIF
    ENDIF
    IF (spec) THEN
       CPPostcondition(lb.EQ.1, cp_warning_level, routineP, error, failure)
       CALL mp_allocate(newarray, ub_new-lb_new+1, stat=stat)
    ELSE
       ALLOCATE(newarray(lb_new:ub_new), stat=stat)
    ENDIF
    IF(stat /= 0) CALL stop_memory(routineP, 'new', ub_new-lb_new+1)
    IF(docopy) THEN
       IF(dbg) CALL cp_assert(lb_new.LE.lb_orig .AND. ub_new.GE.ub_orig,&
            cp_failure_level, cp_internal_error, routineP,&
            'Old extent exceeds the new one.', error, failure)
       newarray(lb_orig:ub_orig) = array(lb_orig:ub_orig)
       IF (spec) THEN
          CALL mp_deallocate(array)
       ELSE
          DEALLOCATE(array)
       ENDIF
    END IF
    array => newarray
    CALL timestop(timing_handle)
  END SUBROUTINE ensure_array_size_i


! *****************************************************************************
!> \brief Ensures that an array is appropriately large.
!> \param[in,out] array       array to verify and possibly resize
!> \param[in] lb    (optional) desired array lower bound
!> \param[in] ub    desired array upper bound
!> \param[in] factor          (optional) factor by which to exagerrate
!>                            enlargements
!> \param[in] nocopy          (optional) copy array on enlargement; default
!>                            is to copy
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE ensure_array_size_i_2d(array, dim, lb, ub, factor, nocopy,&
       special, error)
    INTEGER, DIMENSION(:, :), POINTER        :: array
    INTEGER, INTENT(IN)                      :: dim
    INTEGER, INTENT(IN), OPTIONAL            :: lb
    INTEGER, INTENT(IN)                      :: ub
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, special
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ensure_array_size_i_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lb_new, lb_orig, osize, stat, &
                                                timing_handle, ub_new, ub_orig
    INTEGER, DIMENSION(:, :), POINTER        :: newarray
    LOGICAL                                  :: dbg, docopy = .TRUE., &
                                                failure = .TRUE., info

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    info = .FALSE.
    CALL dbcsr_set_debug(dbg, info=info)
    CALL timeset(routineN, timing_handle)
    docopy = .TRUE.
    IF (PRESENT (nocopy)) docopy = .NOT. nocopy
    lb_new = 1
    IF (PRESENT (lb)) lb_new = lb
    !> Creates a new array if it doesn't yet exist.
    CALL cp_assert(ASSOCIATED(array), cp_failure_level, cp_internal_error, routineP,&
         'cannot resize empty multi-dimensional array', error=error, failure=failure)
    lb_orig = LBOUND(array,dim)
    ub_orig = UBOUND(array,dim)
    IF(dbg) WRITE(*,*)routineP//' Current bounds are',lb_orig,':',ub_orig
    IF (lb_orig.LE.lb_new .AND. ub_orig.GE.ub) THEN
       CALL timestop(timing_handle)
       RETURN
    ENDIF
    IF (lb_orig.GT.lb_new) THEN
       IF (PRESENT(factor)) THEN
          lb_new = MIN(lb_orig, lb_new - NINT(REAL(lb_orig - lb_new)*factor))
       ENDIF
    ENDIF
    IF (ub_orig.LT.ub) THEN
       IF (PRESENT(factor)) THEN
          ub_new = MAX(ub_orig, ub + NINT(REAL(ub-ub_orig)*factor))
       ELSE
          ub_new = ub
       ENDIF
    ELSE
       ub_new = ub
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Resizing to bounds',lb_new,':',ub_new
    IF(.NOT.docopy) DEALLOCATE(array)
    CALL cp_assert(dim.EQ.1, cp_failure_level, cp_unimplemented_error_nr,&
         routineP, 'Only dimension 1 can be specified.', error=error,&
         failure=failure)
    osize = SIZE(array,2)
    ALLOCATE(newarray(lb_new:ub_new,osize), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'new', ub_new-lb_new+1)
    IF(docopy) THEN
       IF(dbg) CALL cp_assert(lb_new.LE.lb_orig .AND. ub_new.GE.ub_orig,&
            cp_failure_level, cp_internal_error, routineP,&
            'Old extent exceeds the new one.', error, failure)
       newarray(lb_orig:ub_orig,:) = array(lb_orig:ub_orig,:)
       DEALLOCATE(array)
    END IF
    array => newarray
    CALL timestop(timing_handle)
  END SUBROUTINE ensure_array_size_i_2d


! *****************************************************************************
!> \brief Adds a coordinate (or other data) into a work matrix's row_i and
!>        col_i arrays and returns its position.
!> \note  Uses the matrix%lastblk to keep track of the current position.
!> \param[in,out] matrix      work matrix
!> \param[in] row,col         row, col data to add
!> \param[in] blk   (optional) block pointer to add
!> \param[out] index          (optional) saved position
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE add_work_coordinate(matrix, row, col, blk, index, error)
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    INTEGER, INTENT(IN), OPTIONAL            :: blk
    INTEGER, INTENT(OUT), OPTIONAL           :: index
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_work_coordinate', &
      routineP = moduleN//':'//routineN

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    dbg = .FALSE.
    !CALL timeset(routineN, timing_handle)
    matrix%lastblk = matrix%lastblk+1
    CALL ensure_array_size(matrix%row_i, ub=matrix%lastblk, error=error)
    CALL ensure_array_size(matrix%col_i, ub=matrix%lastblk, error=error)
    matrix%row_i(matrix%lastblk) = row
    matrix%col_i(matrix%lastblk) = col
    IF (PRESENT(blk)) THEN
       CALL ensure_array_size(matrix%blk_p, ub=matrix%lastblk, error=error)
       matrix%blk_p(matrix%lastblk) = blk
    ENDIF
    IF(dbg.AND.PRESENT(blk))&
         WRITE(*,*)routineP//' Adding',row,col,blk,'at',matrix%lastblk
    IF (dbg.AND.bcsr_verbose) THEN
       WRITE(*,*)routineP//' row_i=',matrix%row_i(1:matrix%lastblk)
       WRITE(*,*)routineP//' col_i=',matrix%col_i(1:matrix%lastblk)
       WRITE(*,*)routineP//' blk_p=',matrix%blk_p(1:matrix%lastblk)
    ENDIF
    IF (PRESENT (index)) index = matrix%lastblk
    !CALL timestop(timing_handle)
  END SUBROUTINE add_work_coordinate


! *****************************************************************************
!> \brief Ensures that the data array is appropriately large.
!> \param[in,out] matrix      ensure data size of this matrix
!> \param[in] size            desired data size
!> \param[in] factor          (optional) factor by which to exagerrate
!>                            enlargements, default is no exagerration
!> \param[in] nocopy          (optional) copy array on enlargement; default
!>                            is to copy
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_ensure_data_size (matrix, size, factor, nocopy, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN)                      :: size
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ensure_data_size', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL ensure_array_size (matrix%data, ub=size, factor=factor, nocopy=nocopy,&
         special=matrix%special_memory, error=error)
  END SUBROUTINE dbcsr_ensure_data_size


! *****************************************************************************
!> \brief Ensures that an array is appropriately large.
!> \param[in,out] array       array to verify and possibly resize
!> \param[in] lb    (optional) desired array lower bound
!> \param[in] ub    desired array upper bound
!> \param[in] factor          (optional) factor by which to exagerrate
!>                            enlargements
!> \param[in] nocopy          (optional) copy array on enlargement; default
!>                            is to copy
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE ensure_array_size_d(array, lb, ub, factor, nocopy, special, error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: array
    INTEGER, INTENT(IN), OPTIONAL            :: lb
    INTEGER, INTENT(IN)                      :: ub
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, special
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ensure_array_size_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lb_new, lb_orig, stat, &
                                                timing_handle, ub_new, ub_orig
    LOGICAL                                  :: dbg, docopy = .TRUE., &
                                                failure = .TRUE., info, spec
    REAL(KIND=dp), DIMENSION(:), POINTER     :: newarray

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    info = .FALSE.
    CALL dbcsr_set_debug(dbg, info=info)
    CALL timeset(routineN,timing_handle)

    docopy = .TRUE.
    IF (PRESENT (nocopy)) docopy = .NOT. nocopy
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    lb_new = 1
    IF (PRESENT (lb)) lb_new = lb
    !> Creates a new array if it doesn't yet exist.
    IF (.NOT.ASSOCIATED(array)) THEN
       failure = .FALSE.
       !CPPrecondition(.NOT.docopy, cp_warning_level, routineP, error, failure)
       !failure = .TRUE.
       IF (spec) THEN
          CPPostcondition(lb_new.EQ.1, cp_warning_level, routineP, error, failure)
          CALL mp_allocate(array, ub-lb_new+1, stat=stat)
       ELSE
          ALLOCATE(array(lb_new:ub), stat=stat)
       ENDIF
       IF (stat /= 0) CALL stop_memory(routineP, 'array', ub-lb+1)
       CALL timestop(timing_handle)
       RETURN
    ENDIF
    lb_orig = LBOUND(array,1)
    ub_orig = UBOUND(array,1)
    IF(dbg) WRITE(*,*)routineP//' Current bounds are',lb_orig,':',ub_orig,&
         '; special?',spec
    IF (lb_orig.LE.lb_new .AND. ub_orig.GE.ub) THEN
       CALL timestop(timing_handle)
       RETURN
    ENDIF
    IF (lb_orig.GT.lb_new) THEN
       IF (PRESENT(factor)) THEN
          lb_new = MIN(lb_orig, lb_new - NINT(REAL(lb_orig - lb_new)*factor))
       ENDIF
    ENDIF
    IF (ub_orig.LT.ub) THEN
       IF (PRESENT(factor)) THEN
          ub_new = MAX(ub_orig, ub + NINT(REAL(ub-ub_orig)*factor))
       ELSE
          ub_new = ub
       ENDIF
    ELSE
       ub_new = ub
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Resizing to bounds',lb_new,':',ub_new
    IF(.NOT.docopy) THEN
       IF (spec) THEN
          CALL mp_deallocate(array)
       ELSE
          DEALLOCATE(array)
       ENDIF
    ENDIF
    IF (spec) THEN
       CPPostcondition(lb.EQ.1, cp_warning_level, routineP, error, failure)
       CALL mp_allocate(newarray, ub_new-lb_new+1, stat=stat)
    ELSE
       ALLOCATE(newarray(lb_new:ub_new), stat=stat)
    ENDIF
    IF(stat /= 0) CALL stop_memory(routineP, 'newarray', ub_new-lb_new+1)
    IF(docopy) THEN
       IF(dbg) CALL cp_assert(lb_new.LE.lb_orig .AND. ub_new.GE.ub_orig,&
            cp_failure_level, cp_internal_error, routineP,&
            'Old extent exceeds the new one.', error, failure)
       !IF(dbg) WRITE(*,*)routineP//' Old array',array(lb_orig:ub_orig)
       newarray(lb_orig:ub_orig) = array(lb_orig:ub_orig)
       IF (spec) THEN
          CALL mp_deallocate(array)
       ELSE
          DEALLOCATE(array)
       ENDIF
    END IF
    array => newarray
    IF(dbg.AND.bcsr_verbose)&
         WRITE(*,*)routineP//' New array',array(lb_orig:ub_orig)
    IF (dbg) WRITE(*,*)routineP//' New array size', SIZE(array)
    CALL timestop(timing_handle)
  END SUBROUTINE ensure_array_size_d


! *****************************************************************************
!> \brief Prints the elements of a matrix.
!> \param[in] rows,columns    is the logical (possibly detransposed) matrix
!>                            size, not the stored size.
!> \param[in] tr    (optional) specifies whether the elements are stored
!>                  transposed
! *****************************************************************************
  SUBROUTINE dbcsr_printmat(matrix, rows, cols, title, tr)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

    CHARACTER(30)                            :: f
    INTEGER                                  :: m, n, r
    LOGICAL                                  :: t
    REAL(KIND=dp)                            :: bit_bucket

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    !IF (.NOT.dbg.OR..NOT.bcsr_verbose) RETURN
    m = rows
    n = cols
    t = .FALSE.
    IF(PRESENT(title)) WRITE(*,*)title
    IF(PRESENT(tr)) THEN
       IF(tr) THEN
          t = .TRUE.
          m = cols
          n = rows
       ENDIF
    ENDIF
    IF (dbg) WRITE(*,*)'Print Matrix 1D:',rows,cols,t,'range',&
         LBOUND(matrix),':',UBOUND(matrix)
    DO r = LBOUND(matrix,1), UBOUND(matrix,1)
       bit_bucket = matrix(r)
    ENDDO
    IF (dbg) WRITE(*,*)'Range OK'
    bit_bucket = 0.0_dp
    DO r = LBOUND(matrix,1), UBOUND(matrix,1)
       !IF (dbg.AND.bcsr_verbose) WRITE(*,*)':@',r,'=',matrix(r)
       bit_bucket = bit_bucket + matrix(r)
    ENDDO
    IF(m.GT.10000) m=0
    IF(n.GT.10000) n=0
    IF(m*n.LT.1 .OR. m*n.GT.SIZE(matrix)) RETURN
    !WRITE(f,FMT="((I4),(F9.4)))")n
    WRITE(f,FMT="('(',I4,'(F9.4))')")cols
    !matrix(:) = (/ (real(r), r=1,m*n) /)
    DO r = 1, rows
       !WRITE(*,FMT=f)matrix(r,1:(n-1)*m+r:n)
       IF(.NOT.t) THEN
          WRITE(*,FMT=f)matrix(r:r+(cols-1)*rows:rows)
       ELSE
          WRITE(*,FMT=f)matrix((r-1)*cols+1:r*cols)
       ENDIF
    ENDDO
  END SUBROUTINE dbcsr_printmat


! *****************************************************************************
!> \brief Prints the elements of a matrix.
!> rows and columns is the matrix size. tr specifies whether
!> the elements are stored transposed
! *****************************************************************************
  SUBROUTINE dbcsr_print2dmat(matrix, title, tr)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: matrix
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

    CHARACTER(30)                            :: f
    INTEGER                                  :: cols, m, n, r, rows
    LOGICAL                                  :: dbg, t

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg, on=.FALSE.)
    !IF (.NOT.dbg.OR..NOT.bcsr_verbose) RETURN
    IF (.NOT.dbg) RETURN
    rows = SIZE(matrix,1)
    cols = SIZE(matrix,2)
    m = rows
    n = cols
    t = .FALSE.
    IF(PRESENT(title)) WRITE(*,*)title
    IF(PRESENT(tr)) THEN
       IF(tr) THEN
          t = .TRUE.
          m = cols
          n = rows
       ENDIF
    ENDIF
    WRITE(*,*)'Print Matrix 2D:',rows,cols,t,'range',&
         LBOUND(matrix),':',UBOUND(matrix)
    IF(m.GT.10000) m=0
    IF(n.GT.10000) n=0
    IF(m*n.LT.1 .OR. m*n.GT.SIZE(matrix)) RETURN
    WRITE(f,FMT="('(',I4,'(F9.4))')")cols
    DO r = 1, rows
       !WRITE(*,FMT=f)matrix(r,1:(n-1)*m+r:n)
       IF(.NOT.t) THEN
          WRITE(*,FMT=f)matrix(r,:)
       ELSE
          WRITE(*,FMT=f)matrix(:,r)
       ENDIF
    ENDDO
  END SUBROUTINE dbcsr_print2dmat


! *****************************************************************************
!> \brief Swaps two integers
!> \param[in,out] a,b         Integers to swap
! *****************************************************************************
  SUBROUTINE iswap(a, b)
    INTEGER, INTENT(INOUT)                   :: a, b

    INTEGER                                  :: tmp

!   ---------------------------------------------------------------------------

    tmp = a
    a = b
    b = tmp
  END SUBROUTINE iswap


! *****************************************************************************
!> \brief In-place block transpose.
!> \param[in,out] extent      Matrix in the form of a 1-d array
!> \param[in] rows, columns   Input matrix size
!> \param error     cp2k error
! *****************************************************************************
  PURE SUBROUTINE dbcsr_block_transpose_inplace(extent, rows, columns)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: extent
    INTEGER, INTENT(IN)                      :: rows, columns

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_block_transpose_inplace', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    extent = RESHAPE(TRANSPOSE(&
         RESHAPE(extent, (/rows, columns/))), (/rows*columns/))
  END SUBROUTINE dbcsr_block_transpose_inplace
! *****************************************************************************
!> \brief Copy and transpose block
!> \param[out] extent_out     output matrix in the form of a 1-d array
!> \param[in] extent_in       input matrix in the form of a 1-d array
!> \param[in] rows, columns   input matrix size
!> \param error     cp2k error
! *****************************************************************************
  PURE SUBROUTINE dbcsr_block_transpose_copy(extent_out, extent_in,&
       rows, columns)
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: extent_out
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: extent_in
    INTEGER, INTENT(IN)                      :: rows, columns

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_block_transpose_copy', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    extent_out = RESHAPE(TRANSPOSE(&
         RESHAPE(extent_in, (/rows, columns/))), (/rows*columns/))
  END SUBROUTINE dbcsr_block_transpose_copy



! *****************************************************************************
!> \brief Copies metadata into an array.
!> \param[in] matrix      Matrix
!> \param[out] meta       Metadata elements
!> \param[in,out] error   cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_pack_meta(matrix, meta, error)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, DIMENSION(dbcsr_meta_size), &
      INTENT(OUT)                            :: meta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_pack_meta', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    meta(1)  = matrix%nblks
    meta(2)  = matrix%nze
    meta(3)  = matrix%nblkrows_total
    meta(4)  = matrix%nblkcols_total
    meta(5)  = matrix%nfullrows_total
    meta(6)  = matrix%nfullcols_total
    meta(7)  = matrix%nblkrows_local
    meta(8)  = matrix%nblkcols_local
    meta(9)  = matrix%nfullrows_local
    meta(10) = matrix%nfullcols_local
    SELECT CASE (matrix%matrix_type)
    CASE ('N')
       meta(11) = 0
    CASE ('T')
       meta(11) = 1
    CASE ('S')
       meta(11) = 2
    CASE ('A')
       meta(11) = 3
    CASE DEFAULT
       meta(11) = -1
       CALL cp_assert(.FALSE., cp_failure_level,&
            cp_wrong_args_error, routineP, "Invalid matrix type.", error)
    END SELECT
  END SUBROUTINE dbcsr_pack_meta


! *****************************************************************************
!> \brief Copies row and column block sizes from another matrix.
!> \param[in,out] matrix      target matrix
!> \param[in] meta            source matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_copy_block_sizes(dst, src)
    TYPE(dbcsr_type), INTENT(INOUT)          :: dst
    TYPE(dbcsr_type), INTENT(IN)             :: src

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_block_sizes', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    dst%row_blk_size = src%row_blk_size
    CALL array_hold (dst%row_blk_size)
    dst%col_blk_size = src%col_blk_size
    CALL array_hold (dst%col_blk_size)
  END SUBROUTINE dbcsr_copy_block_sizes


! *****************************************************************************
!> \brief Sets debugging for an individual subroutine depending on the
!>        value of the bcsr_debug module variable, default initialization in
!>        the declared subroutine, or 
!> \param[in,out] dbg         debug flag
!> \param[in] on    (optional) turn debugging on or off in this subroutine
!> \param[in,out] info        print informational messages (less than debug)
! *****************************************************************************
  PURE SUBROUTINE dbcsr_set_debug(dbg, on, info)
    LOGICAL, INTENT(INOUT)                   :: dbg
    LOGICAL, INTENT(IN), OPTIONAL            :: on
    LOGICAL, INTENT(INOUT), OPTIONAL         :: info

!   ---------------------------------------------------------------------------

    IF (PRESENT (info)) info = info .AND. bcsr_info
    IF (PRESENT (on)) THEN
       dbg = on
       RETURN
    ENDIF
    IF (.NOT. bcsr_debug) THEN
       dbg = .FALSE.
       RETURN
    ENDIF
  END SUBROUTINE dbcsr_set_debug


! *****************************************************************************
!> \brief Creates a DBCSR distribution from a distribution_2d
!> \param[in] dist2d          distribution_2d
!> \param[out] dist           DBCSR distribution
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_dist2d_to_dist(dist2d, dist, error)
    TYPE(distribution_2d_type), INTENT(IN), &
      TARGET                                 :: dist2d
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: mypcol, myproc, myprow, &
                                                numproc
    INTEGER, DIMENSION(:), POINTER           :: col_dist_data, row_dist_data
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    TYPE(array_i1d_obj)                      :: cd, rd
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_mp_obj)                       :: mp_env
    TYPE(distribution_2d_type), POINTER      :: dist2d_p

!
!   ---------------------------------------------------------------------------

    dist2d_p => dist2d
    CALL distribution_2d_get(dist2d_p, error=error,&
         row_distribution=row_dist_data, col_distribution=col_dist_data,&
         blacs_env=blacs_env)
    CALL get_blacs_info(blacs_env, para_env=para_env,&
         my_process_row=myprow, my_process_column=mypcol,&
         blacs2mpi=pgrid)
    myproc = para_env%mepos
    numproc = para_env%num_pe
    CALL dbcsr_mp_new(mp_env, pgrid, para_env%group, myproc, numproc,&
         myprow, mypcol)
    CALL array_nullify (rd)
    CALL array_nullify (cd)
    CALL array_new(rd, row_dist_data)
    CALL array_new(cd, col_dist_data)
    CALL dbcsr_distribution_new(dist, mp_env, rd, cd)
    CALL dbcsr_mp_release (mp_env)
    CALL array_release (rd)
    CALL array_release (cd)
  END SUBROUTINE dbcsr_dist2d_to_dist

! *****************************************************************************
!> \brief Transposes a distribution
!> \param[out] dist_tr        transposed distribution
!> \param[in] dist2d          current distribution
! *****************************************************************************
  SUBROUTINE dbcsr_transpose_distribution(dist_tr, dist_normal)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_tr
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_normal

    TYPE(array_i1d_obj)                      :: cd, lc, lr, rd
    TYPE(dbcsr_mp_obj)                       :: mp_env_tr

!   ---------------------------------------------------------------------------

    rd = dbcsr_distribution_col_dist (dist_normal)
    cd = dbcsr_distribution_row_dist (dist_normal)
    lr = dbcsr_distribution_local_cols (dist_normal)
    lc = dbcsr_distribution_local_rows (dist_normal)
    CALL dbcsr_mp_new_transposed (mp_env_tr, dbcsr_distribution_mp (&
         dist_normal))
    CALL dbcsr_distribution_new(dist_tr, mp_env_tr, rd, cd, lr, lc)
    CALL dbcsr_mp_release (mp_env_tr)
  END SUBROUTINE dbcsr_transpose_distribution



! *****************************************************************************
!> \brief Creates a left-matrix DBCSR repeated distribution for a matrix
!>        multiplication based on the right matrix
!> \par The processor grid is a transpose of the right matrix and the column
!>      distribution is is right matrix's row distribution. The row
!>      is a cycle of the *new* processor rows (i.e., same as processor
!>      columns of the right while the old row distribution becomes the new
!>      image distribution. In other words, the number of row images are the
!>      same as the previous number of processor rows.
!> \param[out] imgdist_product          product distribution repetition
!> \param[in] dist_left, dist_right     distributions of the left and right
!>                                      matrices
! *****************************************************************************
  SUBROUTINE dbcsr_create_image_dist_left(imgdist_left,&
       dist_left, dist_right)
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(OUT)                            :: imgdist_left
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_left, dist_right

    INTEGER                                  :: i, nprows, nrows, row
    INTEGER, DIMENSION(:), POINTER           :: rd, row_dist_data, &
                                                row_rep_data
    INTEGER, DIMENSION(:, :), POINTER        :: right_pgrid
    TYPE(array_i1d_obj)                      :: col_dist, row_dist, &
                                                row_repetition
    TYPE(dbcsr_distribution_obj)             :: new_dist_left
    TYPE(dbcsr_mp_obj)                       :: mp_env, old_mp_env

!
!   ---------------------------------------------------------------------------
! Set up the basic data: the row repetition is the left matrix row
! distribution.
!row_repetition = dbcsr_distribution_row_dist (dist_left)

    CALL dbcsr_distribution_init (new_dist_left)
    CALL dbcsr_distribution_init (imgdist_left%main)
    CALL array_nullify (imgdist_left%col_image)
    CALL array_nullify (imgdist_left%row_image)
    rd => array_data (dbcsr_distribution_row_dist (dist_left))
    ALLOCATE (row_rep_data (SIZE (rd)))
    row_rep_data(:) = rd(:) + 1
    CALL array_new (row_repetition, row_rep_data, gift=.TRUE.)
    nrows = array_size (row_repetition)
    right_pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (dist_right))
    imgdist_left%nrow_images = SIZE (right_pgrid, 1)
    imgdist_left%ncol_images = 1
    ! Calculate the row distribution. It would be nicer to calculate
    ! this in a distributed manner and more load-balanced. For now we must
    ! just be sure all processors calculate the same thing.
    ALLOCATE (row_dist_data (nrows))
    nprows = SIZE (right_pgrid, 2)
    i = 0
    DO row = 1, nrows
       row_dist_data (row) = i
       i = i+1
       IF (i .GE. nprows) i = 0
    ENDDO
    CALL array_new (row_dist, row_dist_data, gift=.TRUE.)
    !
    col_dist = dbcsr_distribution_row_dist (dist_right)
    old_mp_env = dbcsr_distribution_mp (dist_right)
    CALL dbcsr_mp_new_transposed (mp_env, old_mp_env)
    CALL dbcsr_distribution_new(new_dist_left,&
         mp_env,&
         row_dist, col_dist)
    WRITE(*,*)'left pgrid shape',&
         SHAPE (dbcsr_mp_pgrid (dbcsr_distribution_mp (dist_left)))
    ! Now finish the distribution repetition.
    imgdist_left%main = new_dist_left
    !CALL array_hold (row_repetition)
    imgdist_left%row_image = row_repetition
    CALL dbcsr_mp_release (mp_env)
    CALL array_release (row_dist)
  END SUBROUTINE dbcsr_create_image_dist_left

! *****************************************************************************
!> \brief Creates a DBCSR distribution for a matrix multiplication based on
!>        the right matrix
!> \param[out] imgdist_product          product distribution repetition
!> \param[in] dist                      standard distribution on which the
!>                                      images are based
! *****************************************************************************
  SUBROUTINE dbcsr_create_image_dist(imgdist, dist)
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(OUT)                            :: imgdist
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    imgdist%nrow_images = 1
    imgdist%ncol_images = 1
    !CALL array_nullify(imgdist%row_image)
    !CALL array_nullify(imgdist%col_image)
    CALL array_new (imgdist%row_image,&
         (/ (1, i=1, dbcsr_distribution_nrows (dist)) /), lb=1)
    CALL array_new (imgdist%col_image,&
         (/ (1, i=1, dbcsr_distribution_ncols (dist)) /), lb=1)
    imgdist%main = dist
    CALL dbcsr_distribution_hold (imgdist%main)
  END SUBROUTINE dbcsr_create_image_dist

! *****************************************************************************
!> \brief Destroys a DBCSR distribution for a matrix multiplication based on
!>        the right matrix
!> \param[inout] imgdist_product        product distribution repetition
! *****************************************************************************
  SUBROUTINE dbcsr_destroy_image_dist(imgdist)
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(INOUT)                          :: imgdist

!
!   ---------------------------------------------------------------------------

    CALL array_release (imgdist%row_image)
    CALL array_release (imgdist%col_image)
    CALL dbcsr_distribution_release (imgdist%main)
  END SUBROUTINE dbcsr_destroy_image_dist

! *****************************************************************************
!> \brief Splits an array of int8 values into two int4 arrays.
!> \param[in] merged          array of merged values
!> \param[out] array_upper    array to fill with the upper bytes of the merged
!>                            values
!> \param[out] array_lower    array to fill with the lower bytes of the merged
!>                            values
! *****************************************************************************
  SUBROUTINE dbcsr_unpack_i8_2i4 (merged, array_upper, array_lower)
    INTEGER(KIND=int_8), DIMENSION(:), &
      INTENT(IN)                             :: merged
    INTEGER(KIND=int_4), DIMENSION(:), &
      INTENT(OUT)                            :: array_upper, array_lower

    INTEGER(KIND=int_8), PARAMETER           :: lmask8 = 4294967295_int_8

    INTEGER                                  :: i

!
!   ---------------------------------------------------------------------------
! Lmask is used to filter in the lower 4 bytes and so its lower 32 bits are
! set to 1: lmask8 = 2^32-1.
! Umask is used to filter in the higher 4 bytes and so its higher 32 bits
! are set to 1: umask8 = 2^32-1 << 32
!lmask8 = 4294967295 ! 2^32-1
!umask8 = 18446744069414584320 ! (2^32-1) * 2^32 = (2^64-1)-(2^32-1)

!$OMP DO SCHEDULE (STATIC)
    DO i = 1, SIZE (merged)
       array_upper(i) = ISHFT (merged(i), -32)
       array_lower(i) = IAND (merged(i), lmask8)
    ENDDO
!$OMP END DO
  END SUBROUTINE dbcsr_unpack_i8_2i4


! *****************************************************************************
!> \brief Releases all matrices in an arrray.
!> \param[in] source          input matrix
!> \param[in,out] marray      matrix array
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_destroy_array(marray, error)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: marray
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_destroy_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    DO i = LBOUND(marray%mats,1), UBOUND (marray%mats,1)
       CALL dbcsr_destroy (marray%mats(i), error)
    ENDDO
    DEALLOCATE (marray%mats)
    !CALL dbcsr_destroy_image_dist(marray%image_dist)
  END SUBROUTINE dbcsr_destroy_array


! *****************************************************************************
!> \brief Extracts a column from a matrix into a new matrix, local data only.
!> \param[inout] vector_m     matrix into which the column is written
!> \param[in] matrix          matrix in which the column resides
!> \param[in] column          column to extract
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_extract_column (vector_m, matrix, column, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: vector_m
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(IN)                      :: column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_extract_column', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: block_nze, col, nblks, &
                                                next_data_p, row, src_blk, &
                                                src_blk_p, src_nblk_rows_total
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    INTEGER, DIMENSION(dbcsr_meta_size)      :: meta
    LOGICAL                                  :: symmetric_src, transposed_src

!   ---------------------------------------------------------------------------

    symmetric_src = (matrix%matrix_type.EQ.'S'.OR.matrix%matrix_type.EQ.'A')
    transposed_src = (matrix%matrix_type.EQ.'T')
    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)
    !> We first ensure there is enough space in the vector matrix,
    !> then we build it up.
    CALL dbcsr_get_info (matrix, meta=meta, error=error,&
         nblkrows_total = src_nblk_rows_total)
    !> Now go through the rows and find the appropriate column. If it
    !> exists, place the data into the vector matrix and update the
    !> index.
    vector_m%nze = 0
    nblks = 0
    next_data_p = 1
    CALL cp_assert(dbcsr_valid_index (vector_m), cp_warning_level,&
         cp_wrong_args_error, routineP, "Vector matrix must be valid.",&
         error=error)
    vector_m%row_p(1) = 0
    DO row = 1, src_nblk_rows_total
       vector_m%row_p(row+1) = nblks
       DO src_blk = matrix%row_p(row)+1, matrix%row_p(row+1)
          col = matrix%col_i(src_blk)
          !WRITE(*,*)col,'vs.',column,transposed_src,symmetric_src
          IF ((.NOT. transposed_src .AND. col .EQ. column) .OR.&
               (symmetric_src .AND. row .EQ. column) .OR.&
               (transposed_src .AND. row .EQ. column)) THEN
             nblks = nblks + 1
             vector_m%row_p(row+1) = nblks
             vector_m%blk_p(nblks) = next_data_p
             vector_m%col_i(nblks) = 1
             block_nze = row_blk_size(row) * col_blk_size(col)
             src_blk_p = matrix%blk_p(src_blk)
             vector_m%data(next_data_p : next_data_p+block_nze-1) =&
                  matrix%data(src_blk_p : src_blk_p+block_nze-1)
             next_data_p = next_data_p + block_nze
             vector_m%nze = vector_m%nze + block_nze
             WRITE(*,*)'added block',nblks,'in row',row
          ENDIF ! coordinates match
       ENDDO ! src_blk
    ENDDO ! blk_row
    !WRITE(*,*)nblks,'blocks',lbound(vector_m%row_p),ubound(vector_m%row_p)
    !WRITE(*,*)'rowp=',vector_m%row_p
    !WRITE(*,*)'col_i=',vector_m%col_i(1:nblks)
    !WRITE(*,*)'blk_p=',vector_m%blk_p(1:nblks)
    !vector_m%row_p(src_nblk_rows_total+1) = nblks
    vector_m%nblks = nblks
    vector_m%index(dbcsr_slot_nblks) = nblks
    vector_m%index(dbcsr_slot_nze) = vector_m%nze
    ! Are there any other elements we have to update, especially
    ! regarding the number of blocks and the data size?
  END SUBROUTINE dbcsr_extract_column


! *****************************************************************************
!> \brief Returns whether the index structure of the matrix is valid.
!> \param[in] matrix          verify index validity of this matrix
!> \param[out] valid_index    index validity
! *****************************************************************************
  PURE FUNCTION dbcsr_valid_index (matrix) RESULT (valid_index)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    LOGICAL                                  :: valid_index

!   ---------------------------------------------------------------------------

    valid_index = .FALSE.
    IF (ASSOCIATED (matrix%row_p)) THEN
       valid_index = SIZE (matrix%row_p) .GT. 0
    ENDIF
  END FUNCTION dbcsr_valid_index

! *****************************************************************************
!> \brief Sets the diagonal of a DBCSR matrix.
!> \param[in] matrix          matrix in which to set the diagonal.
!> \param[in] value           value of the diagonal elements
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_set_diagonal (matrix, value, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    REAL(KIND=dp), INTENT(IN)                :: value
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_diagonal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p, col, local_row, &
                                                my_pcol, nze, row, row_size
    INTEGER, DIMENSION(:), POINTER           :: local_rows, pcol_dist, &
                                                row_blk_size
    LOGICAL                                  :: in_this_row, valid_index

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL cp_assert(matrix%nblkrows_total .EQ. matrix%nblkcols_total,&
         cp_failure_level, cp_caller_error, routineP,&
         "Can not set diagonal of non-square matrix.", error=error)
    row_blk_size => array_data (matrix%row_blk_size)
    CALL dbcsr_get_info (matrix, local_rows = local_rows,&
         proc_col_dist = pcol_dist,&
         my_pcol = my_pcol, error=error)
    valid_index = dbcsr_valid_index (matrix)
    ! Try to find the column element in each row. If it's not found,
    ! then add it.
    DO local_row = 1, matrix%nblkrows_local
       row = local_rows(local_row)
       ! Only diagonal elements come into play.
       IF (pcol_dist(row) .NE. my_pcol) CYCLE
       in_this_row = .FALSE.
       row_size = row_blk_size(row)
       nze = row_size * row_size
       IF (valid_index) THEN
          blk_loop: DO blk = matrix%row_p(row)+1, matrix%row_p(row+1)
             col = matrix%col_i(blk)
             IF (col .LT. row) CYCLE blk_loop
             IF (col .GT. row) EXIT
             IF (col .EQ. row) THEN
                ! Sets the diagonal in the existing block.
                blk_p = matrix%blk_p(blk)
                CALL set_block_diagonal (matrix%data(blk_p:blk_p+nze-1),&
                     value, row_size)
                in_this_row = .TRUE.
                EXIT blk_loop ! Don't bother with the remaining columns
             ENDIF
          ENDDO blk_loop
       ENDIF ! valid_index
       ! Add the column if it doesn't already exist.
       IF (.NOT. in_this_row) THEN
          DBG 'Adding at row',row,'dimen',row_size,'size',nze,'pos',&
               matrix%w%datasize+1
          CALL add_work_coordinate (matrix%w, row, row, matrix%w%datasize+1,&
               error=error)
          ! add_work_coordinate increments the block counter
          CALL ensure_array_size (matrix%w%data, ub=matrix%w%datasize+nze,&
               error=error)
          CALL set_block_diagonal (&
               matrix%w%data(matrix%w%datasize+1 : matrix%w%datasize+nze),&
               value, row_size)
          matrix%w%datasize = matrix%w%datasize + nze
       ENDIF
    ENDDO ! row
  END SUBROUTINE dbcsr_set_diagonal

! *****************************************************************************
!> \brief Sets the diagonal of a square data block represented as a 1-D array.
!>        Non-diagonal elements are set to 0.
!> \param[out] block_data     sets diagonal in this data block
!> \param[in] value           value of the diagonal elements
!> \param[in] d               dimension of the square data block
! *****************************************************************************
  PURE SUBROUTINE set_block_diagonal (block_data, value, d)
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: block_data
    REAL(KIND=dp), INTENT(IN)                :: value
    INTEGER, INTENT(IN)                      :: d

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    block_data(:) = 0
    FORALL (i = 0 : d-1)
       block_data(1 + i*d + i) = value
    END FORALL
  END SUBROUTINE set_block_diagonal

END MODULE dbcsr_util


