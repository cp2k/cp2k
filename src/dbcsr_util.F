!-----------------------------------------------------------------------------!
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
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release,&
                                             array_size
  USE cp_array_i_utils,                ONLY: cp_1d_i_p_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_retain,&
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
       dbcsr_distribution_mp, dbcsr_distribution_new, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_mp_group, dbcsr_mp_new, dbcsr_mp_numproc, dbcsr_mp_pgrid, &
       dbcsr_mp_release, dbcsr_mp_transposed
  USE dbcsr_types,                     ONLY: &
       dbcsr_array_type, dbcsr_distribution_obj, dbcsr_meta_size, &
       dbcsr_mp_obj, dbcsr_num_slots, dbcsr_rptdistribution_type, &
       dbcsr_slot_blk_p, dbcsr_slot_col_i, dbcsr_slot_nblks, dbcsr_slot_nze, &
       dbcsr_slot_row_p, dbcsr_slot_size, dbcsr_type, dbcsr_work_type
  USE distribution_2d_types,           ONLY: distribution_2d_create,&
                                             distribution_2d_get,&
                                             distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type,&
                                             distribution_2d_write
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_4,&
                                             int_8
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: &
       mp_allocate, mp_alltoall, mp_deallocate, mp_rma_put, mp_sum, &
       mp_window_create, mp_window_destroy, mp_window_fence, mp_window_type
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


  PUBLIC :: dbcsr_create, dbcsr_work_create, dbcsr_finalize,&
       dbcsr_create_dist2d_r, dbcsr_create_dist2d_r_unrot,&
       dbcsr_from_fm, dbcsr_from_sm,&
       dbcsr_checksum, dbcsr_verify_matrix,&
       dbcsr_printmat, dbcsr_print2dmat, dbcsr_repoint_index,&
       dbcsr_pack_meta, dbcsr_unpack_i8_2i4
  PUBLIC :: dbcsr_addto_index_array,&
       addto_array, swap, ensure_array_size,&
       dbcsr_set_debug,&
       dbcsr_destroy, dbcsr_work_destroy,&
       checker_tr,&
       dbcsr_get_lookup, dbcsr_get_block_info,&
       dbcsr_print, dbcsr_normalize, dbcsr_order, add_work_coordinate,&
       dbcsr_destroy_array
  PUBLIC :: dbcsr_dist2d_to_dist, dbcsr_create_dist_left,&
       dbcsr_create_rdist, dbcsr_destroy_rdist



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
  FUNCTION checker_tr(row, column) RESULT(transpose)
    INTEGER, INTENT(IN)                      :: row, column
    LOGICAL                                  :: transpose

!transpose = MOD(column+row,2).NE.0 .EQV. column.GE.row

    transpose = BTEST(column+row, 0) .EQV. column.GE.row
  END FUNCTION checker_tr


! *****************************************************************************
!> \brief Sets the correct source matrix, row, column and possible data
!>        transposition for a given matrix and matrix logical row and
!>        column.
!> \param[in] matrix          bcsr matrix
!> \param[in,out] row         input is logical row, output is lookup row
!> \param[in,out] column      input is logical column, output is lookup
!>                            column
!> \param[out] transpose      whether the data must be transposed
!> \param[out] processor      (optional) returns the processor on which this
!>                            block resides
! *****************************************************************************
  SUBROUTINE dbcsr_get_lookup(matrix, row, column, transpose, processor)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(INOUT)                   :: row, column
    LOGICAL, INTENT(OUT)                     :: transpose
    INTEGER, INTENT(OUT), OPTIONAL           :: processor

    INTEGER                                  :: timing_handle, tmp
    LOGICAL                                  :: straight

!   ---------------------------------------------------------------------------

    CALL timeset("dbcsr_get_lookup", timing_handle)
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
       processor = matrix%dist2d%blacs_env%blacs2mpi(&
            matrix%dist2d%row_distribution(row),&
            matrix%dist2d%col_distribution(column))
    ENDIF
    CALL timestop (timing_handle)
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
  SUBROUTINE dbcsr_get_block_info(matrix, row, column,&
       block_number, data_offset)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(IN)                      :: row, column
    INTEGER, INTENT(OUT)                     :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: data_offset

    INTEGER                                  :: blk, timing_handle
    LOGICAL                                  :: found

!   ---------------------------------------------------------------------------

    CALL timeset("dbcsr_get_block_info", timing_handle)
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
    CALL timestop (timing_handle)
  END SUBROUTINE dbcsr_get_block_info


! *****************************************************************************
!> \brief Creates a matrix, allocating the essentials.
!> \par The matrix itself is allocated, as well as the essential parts of
!>      the index. When passed the nze argument, the data is also allocated
!>      to that size.
!> \param[out] matrix         new matrix
!> \param[in] dist2d          distribution_2d distribution
!> \param[in] matrix_type     'N' for normal, 'T' for transposed, 'S' for
!>                            symmetric, and 'A' for antisymmetric
!> \param[in] meta            sizes for the matrix (equivalent to the meta
!>                            fields of the dbcsr type) but the nze, nblks,
!>                            and type fields are ignored
!> \param[in] nblks           number of blocks
!> \param[in] nze             number of elements
!> \param[in] special         (optional) allocate indices and data using
!>                            special memory
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create(matrix, name, dist2d, matrix_type,&
       meta, nblks, nze, special, error)
    TYPE(dbcsr_type), INTENT(OUT)            :: matrix
    CHARACTER(len=*), INTENT(IN)             :: name
    TYPE(distribution_2d_type), POINTER      :: dist2d
    CHARACTER, INTENT(IN)                    :: matrix_type
    INTEGER, DIMENSION(:), INTENT(IN)        :: meta
    INTEGER, INTENT(IN)                      :: nblks, nze
    LOGICAL, INTENT(IN), OPTIONAL            :: special
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: mt
    INTEGER                                  :: stat

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    DBGI 'Allocating matrix', name
    matrix%valid = .FALSE.
    matrix%name = name
    NULLIFY(matrix%data)
    NULLIFY(matrix%index)
    NULLIFY(matrix%row_p, matrix%col_i, matrix%blk_p,&
         matrix%row_blk_size, matrix%col_blk_size)
    matrix%nblkrows_total = meta(3)
    matrix%nblkcols_total = meta(4)
    matrix%nfullrows_total = meta(5)
    matrix%nfullcols_total = meta(6)
    matrix%nblkrows_local = meta(7)
    matrix%nblkcols_local = meta(8)
    matrix%nfullrows_local = meta(9)
    matrix%nfullcols_local = meta(10)
    matrix%nblks = nblks
    matrix%nze = nze
    matrix%special_memory = .FALSE.
    IF (PRESENT (special)) matrix%special_memory = special
    !
    DBG 'special allocation',matrix%special_memory
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
    matrix%index(:) = 0
    matrix%index(dbcsr_slot_size) = dbcsr_num_slots
    !
    CALL distribution_2d_retain(dist2d, error)
    matrix%dist2d => dist2d
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
    DBGV ' data allocated', ASSOCIATED(matrix%data)
    DBGV ' index allocated', ASSOCIATED(matrix%index)
  END SUBROUTINE dbcsr_create


! *****************************************************************************
!> \brief Creates a working matrix.
!> \param[out] matrix         new matrix
!> \param[in] final           final matrix
!> \param[in] nblks_guess     (optional) estimated number of blocks
!> \param[in] sizedata_guess  (optional) estimated size of data
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_work_create(matrix, final, nblks_guess,&
     sizedata_guess, error)
    TYPE(dbcsr_work_type), INTENT(OUT)       :: matrix
    TYPE(dbcsr_type), INTENT(IN)             :: final
    INTEGER, INTENT(IN), OPTIONAL            :: nblks_guess, sizedata_guess
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_work_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nblks, nze, stat, &
                                                timing_handle

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    matrix%lastblk = 0
    matrix%datasize = 0
    IF(PRESENT(nblks_guess)) THEN
       nblks = nblks_guess
    ELSE
       nblks = final%nblkrows_local * final%nblkcols_local /&
            final%dist2d%blacs_env%para_env%num_pe
    ENDIF
    IF(PRESENT(sizedata_guess)) THEN
       nze = sizedata_guess
    ELSE
       nze = CEILING(REAL(final%nfullrows_local * final%nfullcols_local)/&
            REAL(final%dist2d%blacs_env%para_env%num_pe))
    ENDIF
    NULLIFY(matrix%data)
    ALLOCATE(matrix%data(nze), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'matrix%data', nze)
    !
    ALLOCATE(matrix%row_i(nblks), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'matrix%row_i', nblks)
    !
    ALLOCATE(matrix%col_i(nblks), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'matrix%col_i', nblks)
    !
    ALLOCATE(matrix%blk_p(nblks), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'matrix%blk_p', nblks)
    !
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_work_create


! *****************************************************************************
!> \brief Creates the final dbcsr_type matrix from the working matrix.
!> \par The final matrix is present as the element work_matrix%m.
!> \note The columns and rows must be ordered! We don't support appending
!>       data to an existing matrix.
!> \param[in] wm    work matrix
!> \param[in,out] matrix      final matrix
!> \param[in] resort          whether the indices should be sorted, default
!>                            is true
!> \param[in] reshuffle       whether the data should be reshuffled,
!>                            default is false
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_finalize(wm, matrix, resort, reshuffle, error)
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: wm
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: resort, reshuffle
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p_p, col_i_p, &
                                                col_size, prev_row, row, &
                                                row_p_p, row_size, stat, &
                                                timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sizes
    LOGICAL                                  :: shfl, srt

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    IF (dbg) WRITE(*,*)routineN//' finalizing matrix', matrix%name
    srt = .TRUE.
    IF (PRESENT (resort)) srt = resort
    IF (srt) THEN
       CALL dbcsr_sort_indices(wm%lastblk, wm%row_i, wm%col_i, wm%blk_p,&
            error=error)
       !CALL dbcsr_wm_merge_blocks(wm, error)
    ENDIF
    !CALL cp_assert(matrix%special_memory, cp_warning_level, cp_internal_error,&
    !     routineP,&
    !     'Final matrix data and index may not have been specially allocated',&
    !     error)
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
    IF (.NOT.shfl) THEN
       matrix%blk_p(:) = wm%blk_p(1:matrix%nblks)
       matrix%data(1:wm%datasize) = wm%data(1:wm%datasize)
    ELSE
       ALLOCATE(sizes(matrix%nblks), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineN, "sizes", matrix%nblks)
!$OMP PARALLEL DO PRIVATE(row, blk, row_size, col_size)
       DO row = 1, matrix%nblkrows_total
          row_size = matrix%row_blk_size(row)
          DO blk = matrix%row_p(row)+1, matrix%row_p(row+1)
             col_size = matrix%col_blk_size(matrix%col_i(blk))
             sizes(blk) = row_size*col_size
          ENDDO
       ENDDO
!$OMP END PARALLEL DO
       CALL dbcsr_sort_data(matrix%blk_p, wm%blk_p(1:matrix%nblks),&
            sizes, matrix%data, wm%data, error)
       DEALLOCATE(sizes)
    ENDIF

    wm%datasize = 0
    wm%lastblk = 0
    CALL dbcsr_work_destroy(wm, keepfinal=.TRUE., keepfinaldata=.TRUE.,&
         error=error)
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_finalize


! *****************************************************************************
!> \brief Deallocates and destroys a matrix.
!> \param[in,out] matrix      matrix
!> \param[in,out] error       cp2k error
!> \param[in] no_sizes        do not deallocate row and column block sizes
! *****************************************************************************
  SUBROUTINE dbcsr_destroy(matrix, error, no_sizes)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: no_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_destroy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: nosd

!   ---------------------------------------------------------------------------

    nosd = .FALSE.
    IF (PRESENT (no_sizes)) nosd = no_sizes
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
    IF (.NOT. nosd) THEN
       IF (ASSOCIATED (matrix%row_blk_size)) &
            DEALLOCATE (matrix%row_blk_size)
       IF (ASSOCIATED (matrix%col_blk_size)) &
            DEALLOCATE (matrix%col_blk_size)
    ENDIF
    CALL distribution_2d_release(matrix%dist2d, error)
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
    IF (ASSOCIATED (wm%data).AND..NOT.datasharing.AND..NOT.lkeepdata) &
         DEALLOCATE(wm%data)
    IF (ASSOCIATED (wm%row_i)) DEALLOCATE(wm%row_i)
    IF (ASSOCIATED (wm%col_i)) DEALLOCATE(wm%col_i)
    IF (ASSOCIATED (wm%blk_p)) DEALLOCATE(wm%blk_p)
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
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_para_env_type), POINTER          :: para_env
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
    CALL dbcsr_create(matrix, fm%name, dist2d, 'N', (/0, 0,&
         nblkrows_total, nblkcols_total,&
         nfullrows_total, nfullcols_total,&
         nblkrows_local, nblkcols_local,&
         nfullrows_local, nfullcols_local/),&
         0, 0, error=error)
    ! Somehow assign row & column block sizes.
    ALLOCATE(matrix%row_blk_size(SIZE(row_blk_size)),&
         matrix%col_blk_size(SIZE(col_blk_size)))
    matrix%row_blk_size(:) = row_blk_size(:)
    matrix%col_blk_size(:) = col_blk_size(:)
    CALL dbcsr_work_create(wm, matrix, error=error)
    ALLOCATE(blk_ps(nblkcols_local), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'blk_ps', nblkcols_local+1)
    ! Figure out the maximum row & column block sizes
    m = MAXVAL(matrix%row_blk_size)
    n = MAXVAL(matrix%col_blk_size)
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
             n = matrix%col_blk_size(blkcol)
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
    CALL dbcsr_finalize(wm, matrix, error=error)
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_from_fm



! *****************************************************************************
!> \brief Create a bcsr sparse matrix from an old-style cp2k sparse matrix.
!> \param[out] matrix         bcsr matrix
!> \param[in] sm    cp2k sparse matrix
!> \param[in,out] error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_from_sm(matrix, sm, error)
    TYPE(dbcsr_type), INTENT(OUT)            :: matrix
    TYPE(real_matrix_type), POINTER          :: sm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_from_sm', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: data_exchange_info_cnt = 3

    INTEGER :: absrow, blk, block_col, block_row, datap, i, j, mpicomm, &
      mynod, mypcol, myprow, nblkcols_local, nblkcols_total, nblkrows_local, &
      nblkrows_total, nblks, newp, newpcol, newprow, nfullcols_local, &
      nfullcols_total, nfullrows_local, nfullrows_total, nnon0, npcols, &
      nprows, numnod, nze, p, row_l, send_blkdata_p, sndlen, stat, store_col, &
      store_row, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: data_torecv_cnt, &
                                                data_tosend_cnt
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: data_torecv_info, &
                                                data_tosend_info
    INTEGER, DIMENSION(:), POINTER :: blkcols_local, blkrows_local, &
      col_blk_begin, col_blk_end, pcol_dist, prow_dist, row_blk_begin, &
      row_blk_end
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: dbg, empty, &
                                                failure = .FALSE., symmetric, &
                                                tr
    REAL(KIND=dp)                            :: blk_cs
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: send_blkdata
    REAL(KIND=dp), DIMENSION(:), POINTER     :: blkdata_1d
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blkdata_2d
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_work_type)                    :: wm
    TYPE(distribution_2d_type), POINTER      :: dist2d
    TYPE(mp_window_type)                     :: window
    TYPE(real_block_node_type), POINTER      :: sm_block_node

!   ---------------------------------------------------------------------------
!> The distribution_2d remains the same

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    CALL timeset(routineN, timing_handle)
    dist2d => sm%distribution_2d
    symmetric = (sm%symmetry .EQ. "symmetric")
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
    !CALL get_blacs_info(blacs_env, para_env=para_env, blacs2mpi=blacs2mpi)
    CALL get_blacs_info(blacs_env, para_env=para_env, blacs2mpi=blacs2mpi,&
         my_process_row=myprow, my_process_column=mypcol,&
         number_of_process_rows=nprows, number_of_process_columns=npcols)
    IF (dbg) CALL cp_blacs_env_write(blacs_env, 6, error)
    IF (dbg) WRITE(*,*)routineP//' myprow, mypcol',myprow,'/',nprows,';',&
         mypcol,'/',npcols
    mynod = blacs_env%para_env%mepos
    numnod = blacs_env%para_env%num_pe
    mpicomm = para_env%group
    IF (dbg) WRITE(*,*)'left0=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    CALL dbcsr_create(matrix, sm%name, dist2d, sm%symmetry, (/0, 0,&
         0,0,0,0, 0,0,0,0/), 0, 0, error=error)
    IF (dbg) WRITE(*,*)'left1=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    !> Determine what to send, and build up information about sends.
    ALLOCATE(data_tosend_cnt(0:numnod-1), data_torecv_cnt(0:numnod-1), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'data_to{send,recv}_cnt|pcp',&
         numnod)
    sndlen = MAX(2, (nblkrows_total*nblkcols_total*4)/numnod)
    ALLOCATE(data_tosend_info(data_exchange_info_cnt, sndlen, 0:numnod-1),&
         stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP,'data_tosend_info',&
         data_exchange_info_cnt*sndlen*numnod)
    ALLOCATE(data_torecv_info(data_exchange_info_cnt, sndlen, 0:numnod-1),&
         stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP,'data_torecv_info',&
         data_exchange_info_cnt*sndlen*numnod)
    data_tosend_cnt(:) = 0
    data_tosend_info(1,1,:) = 0
    nze = 0
    DO absrow = 1, nblkrows_total
       newprow = prow_dist(absrow)
       sm_block_node => first_block_node(sm, absrow)
       DO WHILE (ASSOCIATED(sm_block_node))
          CALL get_block_node(block_node=sm_block_node, block_col=block_col,&
               block=blkdata_2d)
          store_row = absrow
          store_col = block_col
          !IF (MOD(store_row + store_col,3).NE.0) THEN
          CALL dbcsr_get_lookup(matrix, store_row, store_col, tr)
          newprow = prow_dist(store_row)
          newpcol = pcol_dist(store_col)
          newp = blacs2mpi(newprow, newpcol)
          nze = (row_blk_end(absrow)-row_blk_begin(absrow)+1) *&
               (col_blk_end(block_col)-col_blk_begin(block_col)+1)
          data_tosend_cnt(newp) = data_tosend_cnt(newp) + nze
          j = data_tosend_info(1, 1, newp) + 1
          data_tosend_info(1, 1, newp) = j
          IF (dbg) WRITE(*,*)routineP//' To send block from',mynod,'to',&
               newp,'block at (',absrow,block_col,') with',nze,&
               'elements, position',j+1,'as (',store_row,',',store_col,')',tr
          CALL cp_assert(j+1.LE.sndlen, cp_failure_level, cp_internal_error,&
               routineP, 'Send buffer not large enough', error=error,&
               failure=failure)
          IF (tr) nze = -nze
          data_tosend_info(:, j+1, newp) = (/ store_row, store_col, nze /)
          !ENDIF
          sm_block_node => next_block_node(sm_block_node)
       END DO ! associated
    ENDDO ! absrow
    !> Now exchange information on the upcoming transfers.
    CALL mp_alltoall(data_tosend_cnt, data_torecv_cnt, 1, para_env%group)
    nnon0 = data_exchange_info_cnt*sndlen
    CALL mp_alltoall(data_tosend_info(:,:,0:numnod-1),&
         &           data_torecv_info(:,:,0:numnod-1),&
         nnon0, para_env%group)
    !> Tally our space requirements.
    DBGV 'data_tosend_cnt',data_tosend_cnt
    DBGV 'data_torecv_cnt',data_torecv_cnt
    nnon0 = SUM(ABS(data_torecv_cnt))
    nblks = SUM(data_torecv_info(1,1,:))
    empty = nblks .EQ. 0
    IF (dbg) WRITE(*,*)'left2=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    CALL dbcsr_destroy(matrix, error)
    IF (dbg) WRITE(*,*)'left3=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    CALL dbcsr_create(matrix, sm%name, dist2d, sm%symmetry, &
         (/nblks, nnon0,&
         nblkrows_total, nblkcols_total, nfullrows_total, nfullcols_total,&
         nblkrows_local, nblkcols_local, nfullrows_local, nfullcols_local/),&
         nblks, nnon0, error=error)
    IF (dbg) WRITE(*,*)'left4=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    CALL dbcsr_work_create(wm, matrix, nblks, nnon0,&
         error=error)
    wm%lastblk = nblks
    wm%datasize = nnon0
    !> Setup local indices
    CALL cp_assert(SIZE(wm%row_i).GE.nblks .AND. SIZE(wm%col_i).GE.nblks .AND.&
         SIZE(wm%blk_p).GE.nblks, cp_failure_level, cp_internal_error, routineP,&
         'Work matrix indices not large enough', error=error, failure=failure)
    blk = 1
    datap = 1
    DO p = 0, numnod-1
       j = data_torecv_info(1,1,p)
       DO i = 1, j
          block_row = data_torecv_info(1, i+1, p)
          block_col = data_torecv_info(2, i+1, p)
          nze = ABS(data_torecv_info(3, i+1, p))
          tr = data_torecv_info(3, i+1, p) .LT. 0
          wm%row_i(blk) = block_row
          wm%col_i(blk) = block_col
          wm%blk_p(blk) = datap
          IF (dbg) WRITE(*,*)routineP//' To recv block from',p,'to',mynod,&
               'block at (',block_row,block_col,') with',nze,&
               'elements, position', i+1,'ptr',datap,tr
          CALL cp_assert(block_row.LE.SIZE(row_blk_end) .AND.&
               block_row.LE.SIZE(row_blk_begin).AND.&
               block_col.LE.SIZE(col_blk_end).AND.&
               block_col.LE.SIZE(col_blk_begin), cp_warning_level,&
               cp_internal_error, routineP,&
               'Row or Column extent arrays too small.', error=error,&
               failure=failure)
          nze = (row_blk_end(block_row)-row_blk_begin(block_row)+1) *&
               (col_blk_end(block_col)-col_blk_begin(block_col)+1)
          CALL cp_assert(ABS(data_torecv_info(3, i+1, p)).EQ.nze,&
               cp_warning_level, cp_internal_error, routineP,&
               'Prediction of data size is invalid.', error=error,&
               failure=failure)
          IF(.NOT.tr) THEN
             data_torecv_info(3, i+1, p) = datap
          ELSE
             data_torecv_info(3, i+1, p) = -datap
          ENDIF
          datap = datap + nze
          blk = blk + 1
       ENDDO
    ENDDO
    CALL dbcsr_sort_indices(nblks, wm%row_i, wm%col_i, wm%blk_p, error=error)
    IF (dbg) WRITE(*,*)routineP//' Data length counted=', datap-1,&
         'specified=', nnon0
    IF (dbg.AND.bcsr_verbose) WRITE(*,*)routineP//' row_i=',wm%row_i
    IF (dbg.AND.bcsr_verbose) WRITE(*,*)routineP//' col_i=',wm%col_i
    CALL cp_assert(datap-1.EQ.nnon0, cp_warning_level,&
         cp_internal_error, routineP, 'Data lengths do not match.',&
         error=error, failure=failure)
    IF (dbg) WRITE(*,*)routineP//' Data length allocated=',&
         SIZE(wm%data),'specified=', datap-1
    CALL cp_assert(SIZE(wm%data).GE.datap-1, cp_warning_level,&
         cp_internal_error, routineP, 'Data storage not large enough.',&
         error=error, failure=failure)
    !> Exchange destination pointers
    CALL mp_alltoall(data_torecv_info(:,:,0:numnod-1),&
         &           data_tosend_info(:,:,0:numnod-1),&
         data_exchange_info_cnt*sndlen, para_env%group)
    !> Setup MPI Windows
    !wm%data(:) = 999
    !matrix%data(:) = 999
    ALLOCATE (send_blkdata(SUM (ABS (data_tosend_cnt))))
    send_blkdata_p = 1
    nnon0 = MAX(nnon0,1)
    CALL mp_window_create(window, wm%data, nnon0, para_env%group)
    IF(dbg) WRITE(*,*)routineP//' Creating window for FP data, sized',&
         SIZE(wm%data),' of kind=',dp,'declaring size',nnon0,'bs=',i
    CALL mp_window_fence(window, nothing_stored=.TRUE., no_pre_access=.TRUE.)
    ! Perform exchanges
    data_tosend_info(1,1,:) = 0
    DO absrow = 1, nblkrows_total
       newprow = prow_dist(absrow)
       sm_block_node => first_block_node(sm, absrow)
       DO WHILE (ASSOCIATED(sm_block_node))
          CALL get_block_node(block_node=sm_block_node, block_col=block_col,&
               block=blkdata_2d)
          store_row = absrow
          store_col = block_col
          !IF (MOD(store_row + store_col,3).NE.0) THEN
          CALL dbcsr_get_lookup(matrix, store_row, store_col, tr)
          newprow = prow_dist(store_row)
          newpcol = pcol_dist(store_col)
          newp = blacs2mpi(newprow, newpcol)
          nze = (row_blk_end(store_row)-row_blk_begin(store_row)+1) *&
               (col_blk_end(store_col)-col_blk_begin(store_col)+1)
          !nze = (row_blk_end(absrow)-row_blk_begin(absrow)+1) *&
          !     (col_blk_end(block_col)-col_blk_begin(block_col)+1)
          data_tosend_cnt(newp) = data_tosend_cnt(newp) + nze
          j = data_tosend_info(1, 1, newp) + 1
          CALL cp_assert(j+1.LE.sndlen, cp_failure_level, cp_internal_error,&
               routineP, 'Send buffer not large enough', error=error,&
               failure=failure)
          CALL cp_assert(data_tosend_info(1, j+1, newp).EQ.store_row .AND.&
               data_tosend_info(2, j+1, newp).EQ.store_col, cp_warning_level,&
               cp_internal_error, routineP, 'Block coordinates do not match.',&
               error=error, failure=failure)
          datap = data_tosend_info(3, j+1, newp)
          CALL cp_assert(datap.LT.0.EQV.tr, cp_warning_level, cp_internal_error,&
               routineP, 'Transpose not consistent', error=error,&
               failure=failure)
          datap = ABS(datap)
          !blk_cs = bl_dot(nze, RESHAPE(blkdata_2d, (/1/)), 1, &
          !     &               RESHAPE(blkdata_2d, (/1/)), 1)
          blk_cs = DOT_PRODUCT(RESHAPE(blkdata_2d, (/nze/)),&
               RESHAPE(blkdata_2d, (/nze/)))
          IF (dbg) WRITE(*,*)routineP//' Putting block from',mynod,'to',newp,&
               'block at (',absrow,block_col,'),(',store_row,store_col,&
               ') with',nze,'els, to', datap, 'sum=',blk_cs,tr
          CALL cp_assert(SIZE(blkdata_2d).EQ.nze, cp_warning_level,&
               cp_internal_error, routineP, 'Block size inequality', error=error,&
               failure=failure)
          IF (.NOT. tr) THEN
#if defined(__parallel)
             !CALL mp_rma_put(window, RESHAPE(blkdata_2d, (/nze/)), nze,&
             !     newp, datap-1)
             blkdata_1d => send_blkdata (send_blkdata_p:send_blkdata_p+nze-1)
             blkdata_1d = RESHAPE(blkdata_2d, (/nze/))
             !IF (newp .NE. para_env%mepos) &
             CALL mp_rma_put(window, blkdata_1d, nze,&
                  newp, datap-1)
             send_blkdata_p = send_blkdata_p + nze
#else
             wm%data(datap:datap+nze-1) = RESHAPE(blkdata_2d, (/nze/))
#endif
             !CALL dbcsr_print2dmat(blkdata_2d)
          ELSE
#if defined(__parallel)
             !CALL mp_rma_put(window, RESHAPE(TRANSPOSE(blkdata_2d), (/nze/)),&
             !     nze, newp, datap-1)
             blkdata_1d => send_blkdata (send_blkdata_p:send_blkdata_p+nze-1)
             blkdata_1d = RESHAPE(TRANSPOSE(blkdata_2d), (/nze/))
             !IF (newp .NE. para_env%mepos) &
             CALL mp_rma_put(window, blkdata_1d,&
                  nze, newp, datap-1)
             send_blkdata_p = send_blkdata_p + nze
#else
             wm%data(datap:datap+nze-1) =&
                  RESHAPE(TRANSPOSE(blkdata_2d), (/nze/))
#endif
             !CALL dbcsr_print2dmat(blkdata_2d, tr=.TRUE.)
          ENDIF
          data_tosend_info(1, 1, newp) = j
          !ENDIF
          sm_block_node => next_block_node(sm_block_node)
       END DO ! associated
    ENDDO
    CALL mp_window_fence(window, no_put_follows=.TRUE., no_post_access=.TRUE.)
    CALL mp_window_destroy(window)
    DEALLOCATE(data_tosend_cnt, data_torecv_cnt, data_tosend_info,&
         data_torecv_info)
    !
    ALLOCATE(matrix%row_blk_size(SIZE(row_blk_begin)),&
         matrix%col_blk_size(SIZE(col_blk_begin)))
    matrix%row_blk_size(:) = row_blk_end(:) - row_blk_begin(:) +1
    matrix%col_blk_size(:) = col_blk_end(:) - col_blk_begin(:) +1
    DBGV 'blocks, datasize',wm%lastblk, wm%datasize,'vs',nblks,nnon0
    !
    CALL dbcsr_finalize(wm, matrix, resort=.FALSE., error=error)
    IF (dbg) WRITE(*,*)'left5=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_from_sm


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
    CALL cp_blacs_env_retain(blacs_env, error)
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
    DBG 'Repointing index for matrix',m%name,'with #blocks',&
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
    !IF (m%index(dbcsr_slot_row_blk_begin).GT.0.AND.all .OR.&
    !     s.EQ.dbcsr_slot_row_blk_begin) THEN
    !   IF(dbg) WRITE(*,*)routineP//' setting pointer rbb'
    !   m%row_blk_begin => m%index(m%index(dbcsr_slot_row_blk_begin):&
    !        &                          m%index(dbcsr_slot_row_blk_begin+1))
    !ENDIF
    !IF (m%index(dbcsr_slot_row_blk_end).GT.0.AND.all .OR.&
    !     s.EQ.dbcsr_slot_row_blk_end) THEN
    !   IF(dbg) WRITE(*,*)routineP//' setting pointer rbe'
    !   m%row_blk_end => m%index(m%index(dbcsr_slot_row_blk_end):&
    !        &                        m%index(dbcsr_slot_row_blk_end+1))
    !ENDIF
    !IF (m%index(dbcsr_slot_col_blk_begin).GT.0.AND.all .OR.&
    !     s.EQ.dbcsr_slot_col_blk_begin) THEN
    !   IF(dbg) WRITE(*,*)routineP//' setting pointer cbb'
    !   m%col_blk_begin => m%index(m%index(dbcsr_slot_col_blk_begin):&
    !        &                           m%index(dbcsr_slot_col_blk_begin+1))
    !ENDIF
    !IF (m%index(dbcsr_slot_col_blk_end).GT.0.AND.all .OR.&
    !     s.EQ.dbcsr_slot_col_blk_end) THEN
    !   IF(dbg) WRITE(*,*)routineP//' setting pointer cbe'
    !   m%col_blk_end => m%index(m%index(dbcsr_slot_col_blk_end):&
    !        &                        m%index(dbcsr_slot_col_blk_end+1))
    !ENDIF

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
    REAL(KIND=dp)                            :: blk_cs, local_cs
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type)                      :: e
    TYPE(cp_para_env_type), POINTER          :: para_env

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    e = error
    CALL distribution_2d_get(matrix%dist2d, blacs_env=blacs_env, error=e)
    CALL get_blacs_info(blacs_env, para_env=para_env)
    local_cs = 0.0d0
    DO br = 1, matrix%nblkrows_total
       m = matrix%row_blk_size(br)
       DO blk = matrix%row_p(br)+1, matrix%row_p(br+1)
          bc = matrix%col_i(blk)
          n = matrix%col_blk_size(bc)
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
    CALL mp_sum(local_cs, para_env%group)
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
    LOGICAL                                  :: tr
    REAL(KIND=dp)                            :: blk_cs
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    IF (.NOT.dbg) RETURN
    CALL timeset(routineN, timing_handle)
    DBG 'Contents of matrix named ', matrix%name
    DBG 'index', matrix%index
    CALL distribution_2d_get(matrix%dist2d, blacs_env=blacs_env, error=error)
    CALL get_blacs_info(blacs_env, para_env=para_env)
    IF (matrix%nblks .GT. 0) THEN
       DO br = 1, matrix%nblkrows_total
          m = matrix%row_blk_size(br)
          DO blk = matrix%row_p(br)+1, matrix%row_p(br+1)
             bc = matrix%col_i(blk)
             n = matrix%col_blk_size(bc)
             mn = m*n
             blk_p = ABS(matrix%blk_p(blk))
             tr = matrix%blk_p(blk) .LT. 0
             blk_cs = DOT_PRODUCT(matrix%data(blk_p:blk_p+mn-1),&
                  &               matrix%data(blk_p:blk_p+mn-1))
             IF(dbg) WRITE(*,*)routineP//' chksum for (',br,',',bc,') at',&
                  blk_p,'l',mn,'= ', blk_cs
             IF(dbg) CALL &
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
! *****************************************************************************
  SUBROUTINE dbcsr_verify_matrix(m, error)
    TYPE(dbcsr_type), INTENT(IN)             :: m
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_verify_matrix', &
      r = moduleN//':'//routineN

    INTEGER :: bc, blk, blk_p, br, i, mb, mn, n, n_full_blocks_total, &
      n_have_blocks_local, n_have_blocks_total, prev_br, timing_handle
    REAL(KIND=dp)                            :: sparsity_total
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    DBG 'Matrix name: ', m%name
    DBGV 'Matrix type: ', m%matrix_type
    CALL distribution_2d_get(m%dist2d, blacs_env=blacs_env, error=error)
    CALL get_blacs_info(blacs_env, para_env=para_env)
    !
    IF(m%matrix_type.EQ.'N' .OR. m%matrix_type.EQ.'T') THEN
       n_full_blocks_total = m%nblkrows_total * m%nblkcols_total
    ELSE
       CALL bcsr_assert(m%nblkrows_total.EQ.m%nblkcols_total, r,&
            'Symmetric matrix is not square')
       n_full_blocks_total = m%nblkrows_total * (m%nblkrows_total+1)/2
    ENDIF
    n_have_blocks_local = m%nblks
    IF(dbg) WRITE(*,*)r//' index=',m%index(:dbcsr_num_slots)
    CALL bcsr_assert(m%index(1).GT.0, r, 'Index size 0')
    DO i = dbcsr_slot_row_p, dbcsr_num_slots
       CALL bcsr_assert(m%index(i) .GT. 0, r, 'Index member is 0')
       IF (i .LE. dbcsr_num_slots-2) &
            CALL bcsr_assert(m%index(i) .LE. m%index(1), r,&
            'Index member is greater than size')
    ENDDO
    !
    IF(dbg) WRITE(*,*)r//' row_p extents',m%index(dbcsr_slot_row_p+1),&
         m%index(dbcsr_slot_row_p), SIZE(m%row_p)
    CALL bcsr_assert(m%index(dbcsr_slot_row_p+1)-m%index(dbcsr_slot_row_p)+1&
         .EQ.m%nblkrows_total+1, r,&
         'Size of row_p index inconsistent with number of rows')
    CALL bcsr_assert(SIZE(m%row_p).EQ.m%nblkrows_total+1, r,&
         'Size of row_p inconsistent with number of rows')
    !
    IF(dbg) WRITE(*,*)r//' col_i extents',m%index(dbcsr_slot_col_i+1),&
         m%index(dbcsr_slot_col_i), SIZE(m%col_i)
    CALL bcsr_assert(m%index(dbcsr_slot_col_i+1)-m%index(dbcsr_slot_col_i)+1&
         .EQ.m%nblks, r,&
         'Size of col_i index inconsistent with number of blocks')
    CALL bcsr_assert(SIZE(m%col_i).EQ.m%nblks, r,&
         'Size of col inconsistent with number of blocks')
    !
    IF(dbg) WRITE(*,*)r//' blk_p extents',m%index(dbcsr_slot_blk_p+1),&
         m%index(dbcsr_slot_blk_p), SIZE(m%blk_p)
    CALL bcsr_assert(m%index(dbcsr_slot_blk_p+1)-m%index(dbcsr_slot_blk_p)+1&
         .EQ.m%nblks, r,&
         'Size of blk_p index inconsistent with number of blocks')
    CALL bcsr_assert(SIZE(m%col_i).EQ.m%nblks, r,&
         'Size of blk_p inconsistent with number of blocks')
    !
!    IF(dbg) WRITE(*,*)r//' row_blk_begin',m%index(dbcsr_slot_row_blk_begin+1),&
!         m%index(dbcsr_slot_row_blk_begin), SIZE(m%row_blk_begin)
!    CALL bcsr_assert(m%index(dbcsr_slot_row_blk_begin+1)-&
!         m%index(dbcsr_slot_row_blk_begin)+1&
!         .EQ.m%nblkrows_total, r,&
!         'Size of row_blk_begin index inconsistent with number of rows')
!    CALL bcsr_assert(SIZE(m%row_blk_begin).EQ.m%nblkrows_total, r,&
!         'Size of row_blk_begin inconsistent with number of rows')
!    !
!    IF(dbg) WRITE(*,*)r//' row_blk_end',m%index(dbcsr_slot_row_blk_end+1),&
!         m%index(dbcsr_slot_row_blk_end), SIZE(m%row_blk_end)
!    CALL bcsr_assert(m%index(dbcsr_slot_row_blk_end+1)-&
!         m%index(dbcsr_slot_row_blk_end)+1&
!         .EQ.m%nblkrows_total, r,&
!         'Size of row_blk_end index inconsistent with number of rows')
!    CALL bcsr_assert(SIZE(m%row_blk_end).EQ.m%nblkrows_total, r,&
!         'Size of row_blk_end inconsistent with number of rows')
!    !
!    IF(dbg) WRITE(*,*)r//' col_blk_begin',m%index(dbcsr_slot_col_blk_begin+1),&
!         m%index(dbcsr_slot_col_blk_begin), SIZE(m%col_blk_begin)
!    CALL bcsr_assert(m%index(dbcsr_slot_col_blk_begin+1)-&
!         m%index(dbcsr_slot_col_blk_begin)+1&
!         .EQ.m%nblkcols_total, r,&
!         'Size of col_blk_begin index inconsistent with number of rows')
!    CALL bcsr_assert(SIZE(m%col_blk_begin).EQ.m%nblkcols_total, r,&
!         'Size of col_blk_begin inconsistent with number of columns')
!    !
!    IF(dbg) WRITE(*,*)r//' col_blk_end',m%index(dbcsr_slot_col_blk_end+1),&
!         m%index(dbcsr_slot_col_blk_end), SIZE(m%col_blk_end)
!    CALL bcsr_assert(m%index(dbcsr_slot_col_blk_end+1)-&
!         m%index(dbcsr_slot_col_blk_end)+1&
!         .EQ.m%nblkcols_total, r,&
!         'Size of col_blk_end index inconsistent with number of columns')
!    CALL bcsr_assert(SIZE(m%col_blk_end).EQ.m%nblkcols_total, r,&
!         'Size of col_blk_end inconsistent with number of columns')
    !
    IF(dbg) THEN
       IF (ASSOCIATED (m%data))&
            WRITE(*,*)r//' nze=',m%nze,'data size',SIZE(m%data)
    ENDIF
    IF (ASSOCIATED (m%data))&
         CALL bcsr_assert(SIZE(m%data).GE.m%nze, r, 'Data storage too small')
    prev_br = 0
    DO br = 1, m%nblkrows_total
       CALL bcsr_assert(m%row_p(br).GE.0, r, 'row_p less than zero')
       IF(br.GT.1) CALL bcsr_assert(m%row_p(br).GE.m%row_p(prev_br),&
            r, 'row_p decreases')
       mb = m%row_blk_size(br)
       CALL bcsr_assert(mb.GT.0, r, 'Row blocked size is zero')
       DO blk = m%row_p(br)+1, m%row_p(br+1)
          CALL bcsr_assert(blk.GT.0, r, 'Block number is zero')
          CALL bcsr_assert(blk.LE.m%nblks, r, 'Block number too high')
          bc = m%col_i(blk)
          CALL bcsr_assert(bc.GT.0, r, 'col_i is zero')
          CALL bcsr_assert(bc.LE.m%nblkcols_total, r, 'col_i too high')
          n = m%col_blk_size(bc)
          CALL bcsr_assert(n.GT.0, r, 'Column blocked size is zero')
          blk_p = m%blk_p(blk)
          mn = mb*n
          CALL bcsr_assert(blk_p.LE.m%nze, r, 'Block pointer too large')
          CALL bcsr_assert(blk_p+mn-1.LE.m%nze, r, 'Block extends too far')
       ENDDO
       prev_br = br
    ENDDO
    DBGV ' row_p=',m%row_p,'size=',SIZE(m%row_p)
    DBGV ' col_i=',m%col_i
    DBGV ' blk_p=',m%blk_p
    CALL mp_sum(n_have_blocks_local, para_env%group)
    n_have_blocks_total = n_have_blocks_local
    sparsity_total = REAL(n_have_blocks_total)/REAL(n_full_blocks_total)*100.0
    !WRITE(*,FMT='(30A,F5.1,A)')r//' Sparsity: ', sparsity_total,'%'
    IF (dbg) WRITE(*,FMT='(A,F5.1,A)')' '//r//' Non-sparsity: ',&
         sparsity_total,'%'
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_verify_matrix


! *****************************************************************************
!> \brief Glue for the real assert code
!> \param[in] l     Assertion
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
! *****************************************************************************
  SUBROUTINE bcsr_assert(l, routine, msg)
    LOGICAL, INTENT(IN)                      :: l
    CHARACTER(len=*), INTENT(IN)             :: routine, msg

!   ---------------------------------------------------------------------------

    CALL cp_assert(l, cp_warning_level, cp_internal_error, routine, msg)
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
    REAL(KIND=dp), DIMENSION(:), POINTER     :: new_data

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    ! We should also check that the indices actually exist.
    CALL cp_assert(ASSOCIATED(matrix%data), cp_failure_level, cp_caller_error,&
         routineN, "Matrix has no data.", error)
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
       row_size = matrix%row_blk_size(row)
       DO blk = matrix%row_p(row)+1, matrix%row_p(row+1)
          col_size = matrix%col_blk_size(matrix%col_i(blk))
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
    CALL cp_assert(SIZE(dst).GE.SIZE(src), cp_failure_level,&
         cp_wrong_args_error, routineN, 'Destination data space too small.',&
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
!> \brief Prepares a matrix to be multipiled. Chooses the correct normalization
!>        routine based on the type of matrix and what normalization must take
!>        place
!> \param[in] source          input matrix
!> \param[in,out] normalized  normalied matrix
!> \param[in] target_rptdist     normalize to which distribution repetition
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_normalize(source, normalized, target_rptdist, error)
    TYPE(dbcsr_type), INTENT(IN)             :: source
    TYPE(dbcsr_array_type), INTENT(OUT)      :: normalized
    TYPE(dbcsr_rptdistribution_type), &
      INTENT(IN)                             :: target_rptdist
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_normalize', &
      routineP = moduleN//':'//routineN

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    DBG 'Normalizing',source%matrix_type,source%name
    SELECT CASE (source%matrix_type)
    CASE ('N')
       CALL dbcsr_normalize_normal(source, normalized, target_rptdist, error)
    CASE ('S')
       CALL dbcsr_normalize_symmetric(source, normalized, target_rptdist,&
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
  SUBROUTINE dbcsr_normalize_symmetric(sm, ums, target_rptdist, error)
    TYPE(dbcsr_type), INTENT(IN)             :: sm
    TYPE(dbcsr_array_type), INTENT(OUT), &
      TARGET                                 :: ums
    TYPE(dbcsr_rptdistribution_type), &
      INTENT(IN)                             :: target_rptdist
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
    INTEGER, DIMENSION(:), POINTER           :: col_dist, row_dist, &
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
    target_dist = target_rptdist%main
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    row_rep_dist => array_data (target_rptdist%row_repetition)
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numproc = dbcsr_mp_numproc (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    DBGV 'row_dist', row_dist
    DBGV 'col_dist', col_dist
    DBG 'row, col maxval:', MAXVAL(row_dist), MAXVAL(col_dist)
    DBG 'blacs2mpi',blacs2mpi
    DBG 'blacs2mpi LB',LBOUND(blacs2mpi)
    DBG 'blacs2mpi UB',UBOUND(blacs2mpi)
    row_repeats = target_rptdist%nrow_repeats
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
       CALL dbcsr_create(ums%mats(row_rep), 'norm '//sm%name, sm%dist2d,&
            'N', new_meta,&
            0, 0, error=error)
       CALL dbcsr_copy_block_sizes(ums%mats(row_rep), sm)
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
             nze = sm%row_blk_size(stored_row) * sm%col_blk_size(stored_col)
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
             nze = sm%row_blk_size(stored_row) * sm%col_blk_size(stored_col)
             send_meta(smp(row_rep, dst_p)) = stored_row
             send_meta(smp(row_rep, dst_p)+1) = stored_col
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
       CALL dbcsr_work_create(wms(row_rep), ums%mats(row_rep),&
            SUM(recv_count(1,row_rep,:)),&
            SUM(recv_count(2,row_rep,:)), error=error)
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
             !     '; blk_p',blk_p,'; row_rep, rptdist',row_rep,&
             !     row_rep_dist (ABS (stored_row))
             DBG 'blk: p',blks(row_rep),src_p,'; row, col',&
                  stored_row,stored_col,&
                  '; blk_ps',blk_ps,'; row_rep, rptdist',row_rep,&
                  row_rep_dist (ABS (stored_row))
             !row_rep = row_rep_dist (ABS(stored_row))
             wms(row_rep)%row_i(blks(row_rep)) = ABS(stored_row)
             wms(row_rep)%col_i(blks(row_rep)) = stored_col
             wms(row_rep)%blk_p(blks(row_rep)) =&
                  SIGN(blk_ps(row_rep), stored_row)
             nze = sm%row_blk_size(ABS(stored_row))&
                  * sm%col_blk_size(stored_col)
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
       CALL dbcsr_finalize(wms(row_rep), ums%mats(row_rep),&
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
  SUBROUTINE dbcsr_normalize_normal(sm, ums, target_rptdist, error)
    TYPE(dbcsr_type), INTENT(IN)             :: sm
    TYPE(dbcsr_array_type), INTENT(OUT), &
      TARGET                                 :: ums
    TYPE(dbcsr_rptdistribution_type), &
      INTENT(IN)                             :: target_rptdist
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_normalize_normal', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, blk_l, blk_p, col, dst_p, mp_group, numproc, nze, pcol, &
      prow, row, src_p, stored_col, stored_row, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                sd_disp, sdp, send_meta, &
                                                sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count
    INTEGER, DIMENSION(:), POINTER           :: col_dist, row_dist
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
    target_dist = target_rptdist%main
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numproc = dbcsr_mp_numproc (mp_obj)
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
    CALL dbcsr_create(ums%mats(1), 'norm '//sm%name, sm%dist2d, 'N', new_meta,&
         0, 0, error=error)
    um => ums%mats(1)
    CALL dbcsr_copy_block_sizes(um, sm)
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
          nze = sm%row_blk_size(stored_row) * sm%col_blk_size(stored_col)
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
          nze = sm%row_blk_size(stored_row) * sm%col_blk_size(stored_col)
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
    CALL dbcsr_work_create(wm, um, SUM(recv_count(1,:)),&
         SUM(recv_count(2,:)), error=error)
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
          nze = sm%row_blk_size(stored_row) * sm%col_blk_size(stored_col)
          blk_p = blk_p + nze
          blk = blk + 1
       ENDDO
    ENDDO
    !
    wm%lastblk = blk-1
    wm%datasize = blk_p - 1
    DBG 'Finalizing normalization'
    CALL dbcsr_finalize(wm, um, resort=.TRUE., reshuffle=.TRUE., error=error)
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
       CPPrecondition(.NOT.docopy, cp_warning_level, routineP, error, failure)
       failure = .TRUE.
       IF (spec) THEN
          CPPostcondition(lb.EQ.1, cp_warning_level, routineP, error, failure)
          CALL mp_allocate(array, ub-lb+1, stat=stat)
       ELSE
          ALLOCATE(array(lb:ub), stat=stat)
       ENDIF
       IF (stat /= 0) CALL stop_memory(routineP, 'array', ub-lb+1)
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

    INTEGER                                  :: timing_handle

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    dbg = .FALSE.
    CALL timeset(routineN, timing_handle)
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
    CALL timestop(timing_handle)
  END SUBROUTINE add_work_coordinate


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
       CPPrecondition(.NOT.docopy, cp_warning_level, routineP, error, failure)
       !failure = .TRUE.
       IF (spec) THEN
          CPPostcondition(lb.EQ.1, cp_warning_level, routineP, error, failure)
          CALL mp_allocate(array, ub-lb+1, stat=stat)
       ELSE
          ALLOCATE(array(lb:ub), stat=stat)
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
    IF (.NOT.dbg) RETURN
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
    bit_bucket = 0.0d0
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
  SUBROUTINE dbcsr_block_transpose(extent, rows, columns, error)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: extent
    INTEGER, INTENT(IN)                      :: rows, columns
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_block_transpose', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: dbg

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    IF (PRESENT (error)) THEN
       IF (dbg) WRITE(*,*)routineP//' Array size=',SIZE(extent),&
            'declared size=',rows*columns,'=',rows,'x',columns
       CALL cp_assert(SIZE(extent).GE.rows*columns,&
            cp_failure_level, cp_precondition_failed, routineP,&
            'Array to transpose is smaller than declared size!')
    ENDIF
    extent = RESHAPE(TRANSPOSE(&
         RESHAPE(extent, (/rows, columns/))), (/rows*columns/))
  END SUBROUTINE dbcsr_block_transpose


! *****************************************************************************
!> \brief Copies metadata into an array.
!> \param[in] matrix      Matrix
!> \param[out] meta       Metadata elements
!> \param[in,out] error   cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_pack_meta(matrix, meta, error)
    TYPE(dbcsr_type)                         :: matrix
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

    INTEGER                                  :: stat

!   ---------------------------------------------------------------------------

    ALLOCATE(dst%row_blk_size(SIZE(src%row_blk_size)), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineN, "dst%row_blk_sizes",&
         SIZE(src%row_blk_size))
    ALLOCATE(dst%col_blk_size(SIZE(src%col_blk_size)), stat=stat)
    IF (stat /= 0) CALL stop_memory(routineN, "dst%col_blk_sizes",&
         SIZE(src%col_blk_size))
    dst%row_blk_size(:) = src%row_blk_size(:)
    dst%col_blk_size(:) = src%col_blk_size(:)
  END SUBROUTINE dbcsr_copy_block_sizes


! *****************************************************************************
!> \brief Sets debugging for an individual subroutine depending on the
!>        value of the bcsr_debug module variable, default initialization in
!>        the declared subroutine, or 
!> \param[in,out] dbg         debug flag
!> \param[in] on    (optional) turn debugging on or off in this subroutine
!> \param[in,out] info        print informational messages (less than debug)
! *****************************************************************************
  SUBROUTINE dbcsr_set_debug(dbg, on, info)
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
    CALL array_new(rd, row_dist_data)
    CALL array_new(cd, col_dist_data)
    CALL dbcsr_distribution_new(dist, mp_env, rd, cd)
    CALL dbcsr_mp_release (mp_env)
    CALL array_release (rd)
    CALL array_release (cd)
  END SUBROUTINE dbcsr_dist2d_to_dist


! *****************************************************************************
!> \brief Creates a left-matrix DBCSR repeated distribution for a matrix
!>        multiplication based on the right matrix
!> \param[out] rptdist_product          product distribution repetition
!> \param[in] dist_left, dist_right     distributions of the left and right
!>                                      matrices
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_left(rptdist_left,&
       dist_left, dist_right)
    TYPE(dbcsr_rptdistribution_type), &
      INTENT(OUT)                            :: rptdist_left
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

    rd => array_data (dbcsr_distribution_row_dist (dist_left))
    ALLOCATE (row_rep_data (SIZE (rd)))
    row_rep_data(:) = rd(:) + 1
    CALL array_new (row_repetition, row_rep_data, gift=.TRUE.)
    nrows = array_size (row_repetition)
    right_pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (dist_right))
    rptdist_left%nrow_repeats = SIZE (right_pgrid, 1)
    rptdist_left%ncol_repeats = 1
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
    mp_env = dbcsr_mp_transposed (old_mp_env)
    CALL dbcsr_distribution_new(new_dist_left,&
         mp_env,&
         row_dist, col_dist)
    WRITE(*,*)'left pgrid shape',&
         SHAPE (dbcsr_mp_pgrid (dbcsr_distribution_mp (dist_left)))
    ! Now finish the distribution repetition.
    rptdist_left%main = new_dist_left
    !CALL array_hold (row_repetition)
    rptdist_left%row_repetition = row_repetition
    CALL array_nullify (rptdist_left%col_repetition)
    CALL dbcsr_mp_release (mp_env)
    CALL array_release (row_dist)
  END SUBROUTINE dbcsr_create_dist_left

! *****************************************************************************
!> \brief Creates a DBCSR distribution for a matrix multiplication based on
!>        the right matrix
!> \param[out] rptdist_product          product distribution repetition
!> \param[in] dist                      standard distribution on which the
!>                                      repetition is based
! *****************************************************************************
  SUBROUTINE dbcsr_create_rdist(rptdist, dist)
    TYPE(dbcsr_rptdistribution_type), &
      INTENT(OUT)                            :: rptdist
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

!
!TYPE(array_i1d_obj)                           :: row_repetition,&
!                                                 row_dist
!   ---------------------------------------------------------------------------

    rptdist%nrow_repeats = 1
    rptdist%ncol_repeats = 1
    CALL array_nullify(rptdist%row_repetition)
    CALL array_nullify(rptdist%col_repetition)
    CALL dbcsr_distribution_hold (dist)
    rptdist%main = dist
  END SUBROUTINE dbcsr_create_rdist

! *****************************************************************************
!> \brief Destroys a DBCSR distribution for a matrix multiplication based on
!>        the right matrix
!> \param[inout] rptdist_product        product distribution repetition
! *****************************************************************************
  SUBROUTINE dbcsr_destroy_rdist(rptdist)
    TYPE(dbcsr_rptdistribution_type), &
      INTENT(INOUT)                          :: rptdist

!
!   ---------------------------------------------------------------------------

    CALL array_release (rptdist%row_repetition)
    CALL array_release (rptdist%col_repetition)
    CALL dbcsr_distribution_release (rptdist%main)
  END SUBROUTINE dbcsr_destroy_rdist

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
!> \param[in] no_sizes        do not deallocate row and column block sizes
! *****************************************************************************
  SUBROUTINE dbcsr_destroy_array(marray, error, no_sizes)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: marray
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: no_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_destroy_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: nosd

!   ---------------------------------------------------------------------------

    nosd = .FALSE.
    IF (PRESENT (no_sizes)) nosd = no_sizes
    DO i = LBOUND(marray%mats,1), UBOUND (marray%mats,1)
       CALL dbcsr_destroy (marray%mats(i), error, no_sizes = nosd)
    ENDDO
    DEALLOCATE (marray%mats)
  END SUBROUTINE dbcsr_destroy_array


END MODULE dbcsr_util


