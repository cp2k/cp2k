!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR sparse matrix utility routines
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - Created 2009-05-12
! *****************************************************************************
MODULE dbcsr_util

  USE array_types,                     ONLY: array_data,&
                                             array_hold,&
                                             array_i1d_obj
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: dp,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE dbcsr_machine,                   ONLY: m_walltime
  USE dbcsr_message_passing,           ONLY: mp_sum
  USE dbcsr_methods,                   ONLY: &
       dbcsr_distribution_local_cols, dbcsr_distribution_local_rows, &
       dbcsr_distribution_mp, dbcsr_distribution_ncols, &
       dbcsr_distribution_nlocal_cols, dbcsr_distribution_nlocal_rows, &
       dbcsr_distribution_nrows, dbcsr_get_data, dbcsr_get_data_size, &
       dbcsr_mp_group, dbcsr_valid_index
  USE dbcsr_types,                     ONLY: &
       dbcsr_distribution_obj, dbcsr_magic_number, dbcsr_meta_size, &
       dbcsr_num_slots, dbcsr_obj, dbcsr_slot_blk_p, dbcsr_slot_col_i, &
       dbcsr_slot_dense, dbcsr_slot_nblkcols_local, &
       dbcsr_slot_nblkcols_total, dbcsr_slot_nblkrows_local, &
       dbcsr_slot_nblkrows_total, dbcsr_slot_nblks, &
       dbcsr_slot_nfullcols_local, dbcsr_slot_nfullcols_total, &
       dbcsr_slot_nfullrows_local, dbcsr_slot_nfullrows_total, &
       dbcsr_slot_nze, dbcsr_slot_row_p, dbcsr_slot_type, dbcsr_type, &
       dbcsr_type_complex_4, dbcsr_type_complex_8, dbcsr_type_real_4, &
       dbcsr_type_real_8
  USE f77_blas
  USE f77_blas_generic

  !$ USE OMP_LIB
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_util'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  ! Main
  PUBLIC :: dbcsr_checksum, dbcsr_verify_matrix,&
            dbcsr_pack_meta, dbcsr_unpack_meta, meta_from_dist
  PUBLIC :: dbcsr_set_debug, dbcsr_sleep
  ! Block sizes and arrays
  PUBLIC :: dbcsr_copy_block_sizes,&
            convert_sizes_to_offsets, convert_offsets_to_sizes,&
            global_offsets_to_local,&
            dbcsr_unpack_i8_2i4, make_coordinate_tuple,&
            nfull_elements
  PUBLIC :: swap, sgn
  ! math routines
  PUBLIC :: gcd, lcm

  INTERFACE swap
     MODULE PROCEDURE iswap
  END INTERFACE

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
!> \brief Finds block to which a full element belongs.
!> \par Assumptions
!>        It is assumed that block_start and block_end are sorted and
!>        that hint is in the range [0, nblocks].
!> \param[in] full            full element
!> \param[out] block          block to which full belongs
!> \param[in] block_start     starting full elements of blocks
!> \param[out] block_end      ending full elements of blocks
!> \param[in] hint            where to start looking; ignored if 0
! *****************************************************************************
  SUBROUTINE find_block_of_element(full, block, nblocks,&
       block_start, block_end, hint, error)
    INTEGER, INTENT(in)                      :: full
    INTEGER, INTENT(out)                     :: block
    INTEGER, INTENT(in)                      :: nblocks
    INTEGER, DIMENSION(1:nblocks), &
      INTENT(in)                             :: block_start, block_end
    INTEGER, INTENT(in)                      :: hint
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    INTEGER                                  :: count

    IF (hint .NE. 0) THEN
       block = hint
    ELSE
       block = MAX(1,(nblocks+1)/2)
    ENDIF
    count = 0
    DO WHILE (block_start(block) .GT. full .OR. block_end(block) .LT. full)
       IF (block_start(block) .GT. full) THEN
          block = block-1
       ELSEIF (block_end(block) .LT. full) THEN
          block = block+1
       ENDIF
       count = count+1
       IF (count .GT. nblocks .OR. block .LT. 1 .OR. block .GT. nblocks) THEN
          WRITE(*,'(1X,A,I9,A,I9,A)')"Want to find block",&
               block," of",nblocks," blocks"
          !WRITE(*,'(10(1X,I7))')block_start
          !WRITE(*,'(10(1X,I7))')block_end
          !WRITE(*,*)"full, nblocks, count, hint", full, nblocks, count, hint
          CALL dbcsr_assert (count .LE. nblocks, dbcsr_failure_level,&
               dbcsr_internal_error, "find_block_of_element",&
               "Too many searches",__LINE__,error)
       ENDIF
    ENDDO
  END SUBROUTINE find_block_of_element



! *****************************************************************************
!> \brief The sum of a subset of rows/columns
!> \param[in] all_offsets     ordered offsets of all the elements
!> \param[in] local_elements  enumerated local elements
!> \result nfull_elements     sum of sizes of local elemetns
!> \note Used for making matrices dense/undense
! *****************************************************************************
  PURE FUNCTION nfull_elements (all_offsets, local_elements)
    INTEGER, DIMENSION(:), INTENT(IN)        :: all_offsets, local_elements
    INTEGER                                  :: nfull_elements

    INTEGER                                  :: el, lel

    nfull_elements = 0
    DO lel = 1, SIZE(local_elements)
       el = local_elements(lel)
       nfull_elements = nfull_elements + all_offsets(el+1) - all_offsets(el)
    ENDDO
  END FUNCTION nfull_elements



! *****************************************************************************
!> \brief Converts sizes to offsets
!>
!> \param[in] sizes           array with sizes
!> \param[out] offsets_start  offsets of starts
!> \param[out] offsets_stop   (optional) offsets of ends
! *****************************************************************************

  PURE SUBROUTINE convert_sizes_to_offsets (sizes,&
       offsets_start, offsets_stop)
    INTEGER, DIMENSION(:), INTENT(IN)        :: sizes
    INTEGER, DIMENSION(:), INTENT(OUT)       :: offsets_start
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: offsets_stop

    CHARACTER(len=*), PARAMETER :: routineN = 'convert_sizes_to_offsets', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n

!   ---------------------------------------------------------------------------

    n = SIZE (sizes)
    !IF (bcsr_debug) THEN
    !   CALL dbcsr_assert (SIZE (offsets_start) .GE. n, dbcsr_failure_level,&
    !        dbcsr_caller_error, routineN, "Start offset array too small")
    !   CALL dbcsr_assert (SIZE (offsets_start) .EQ. n, dbcsr_warning_level,&
    !        dbcsr_caller_error, routineN, "Start offset array too big")
    !   IF (PRESENT (offsests_stop)) THEN
    !      CALL dbcsr_assert (SIZE (offsets_stop) .GE. n, dbcsr_failure_level,&
    !           dbcsr_caller_error, routineN, "Stop offset array too small")
    !      CALL dbcsr_assert (SIZE (offsets_stop) .EQ. n, dbcsr_warning_level,&
    !           dbcsr_caller_error, routineN, "Stop offset array too big")
    !   ENDIF
    !ENDIF
    IF (n .GT. 0) THEN
       offsets_start(1) = 1
       IF (PRESENT (offsets_stop)) offsets_stop(1) = sizes(1)
       IF (.NOT. PRESENT (offsets_stop)) THEN
          DO i = 2, n
             offsets_start(i) = offsets_start(i-1) + sizes(i-1)
          ENDDO
          IF(SIZE(offsets_start).GT.n) &
               offsets_start(n+1) = offsets_start(n) + sizes(n)
       ELSE
          DO i = 2, n
             offsets_start(i) = offsets_start(i-1) + sizes(i-1)
             offsets_stop(i) = offsets_stop(i-1) + sizes(i)
          ENDDO
          IF(SIZE(offsets_start).GT.n) &
               offsets_start(n+1) = offsets_start(n) + sizes(n)
       ENDIF
    ENDIF
  END SUBROUTINE convert_sizes_to_offsets

! *****************************************************************************
!> \brief Converts offsets to sizes
!>
!> If the offsets of ends are not given, then the array of sizes is assumed
!> to be one greater than the desired sizes.
!>
!> \param[in] offsets_start  offsets of starts
!> \param[out] sizes         array with sizes
!> \param[in] offsets_stop   (optional) offsets of ends
! *****************************************************************************

  PURE SUBROUTINE convert_offsets_to_sizes (offsets_start, sizes, offsets_stop)
    INTEGER, DIMENSION(:), INTENT(IN)        :: offsets_start
    INTEGER, DIMENSION(:), INTENT(OUT)       :: sizes
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: offsets_stop

    CHARACTER(len=*), PARAMETER :: routineN = 'convert_offsets_to_sizes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n

!   ---------------------------------------------------------------------------

    n = SIZE (offsets_start)
    IF (PRESENT (offsets_stop)) THEN
       sizes(:) = offsets_stop(:) - offsets_start(:) + 1
    ELSE
       IF (n .GT. 1) THEN
          DO i = 1, n-1
             sizes(i) = sizes(i+1) - sizes(i)
          ENDDO
       ENDIF
    ENDIF
  END SUBROUTINE convert_offsets_to_sizes


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
!> \brief Converts global offsets to local
!> \par Global vs. Local Indexing
!>      local_offsets may be sized according to the
!>      local index (|local_elements+|1) or the
!>      global index (|global_offsets|).
!> \param[in] global_offsets   Offsets of elements in the global grid
!> \param[in] local_elements   Which elements are local
!> \param[out] local_offsets   Offsets of local elements.
! *****************************************************************************
  SUBROUTINE global_offsets_to_local (global_offsets,&
       local_elements, local_offsets)
    INTEGER, DIMENSION(:), INTENT(IN)        :: global_offsets, local_elements
    INTEGER, DIMENSION(:), INTENT(OUT)       :: local_offsets

    CHARACTER(len=*), PARAMETER :: routineN = 'global_offsets_to_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acc, el, lel, nglobal, nlo, &
                                                nlocal, prev_el, sz
    LOGICAL                                  :: local
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    nglobal = SIZE (global_offsets)-1
    nlocal = SIZE (local_elements)
    nlo = SIZE (local_offsets)-1
    local = .NOT. (nglobal .EQ. nlo)
    IF (local) THEN
       CALL dbcsr_assert (nlocal, "EQ", nlo,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid size for local offsets",__LINE__,error)
    ENDIF
    IF (local) THEN
       acc = 1
       DO lel = 1, nlocal
          local_offsets(lel) = acc
          el = local_elements(lel)
          sz = global_offsets(el+1)-global_offsets(el)
          acc = acc + sz
       ENDDO
       local_offsets(nlocal+1) = acc
    ELSE
       acc = 1
       prev_el=0
       DO lel = 1, nlocal
          el = local_elements(lel)
          local_offsets(prev_el+1:el) = acc
          sz = global_offsets(el+1)-global_offsets(el)
          acc = acc + sz
          prev_el = el
       ENDDO
       local_offsets(prev_el+1:nglobal+1) = acc
    ENDIF
  END SUBROUTINE global_offsets_to_local


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
!> \brief Calculates explicit sizes for all data blocks.
!> \param[out] sizes          sizes of all data blocks
!> \param[in] row_p, col_i    index structure
!> \param[in] rbs, cbs        row block sizes and column block sizes
! *****************************************************************************
  SUBROUTINE calc_blk_sizes(sizes, row_p, col_i, rbs, cbs)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: sizes
    INTEGER, DIMENSION(:), INTENT(IN)        :: row_p, col_i, rbs, cbs

    INTEGER                                  :: blk, nrows, row, row_size

    nrows = SIZE (row_p)-1
    !$OMP PARALLEL DO PRIVATE (row, blk, row_size)
    DO row = 1, nrows
       row_size = rbs(row)
       FORALL (blk = row_p(row)+1 : row_p(row+1))
          sizes(blk) = row_size * cbs(col_i(blk))
       END FORALL
    ENDDO
    !$OMP END PARALLEL DO
  END SUBROUTINE calc_blk_sizes



! *****************************************************************************
!> \brief Calculates the intersections of blocks
!> \par nints output format
!>      The ints array should be up to twice as large as the number of
!>      intersecting blocks. Each entry is comprised of the target
!>      block and the common length.
!> \par n_src_dsts format
!>      This arrays stored the number of intersecting blocks in common
!>      (position 2) and the offset of the first common intersecting
!>      block (position 1).
!> \note This routine is used in the counting and sending loops in
!>       dbcsr_complete_redistribute
!> \param[out] ints           intersections of blocks
!> \param[in] num_ints        maximum number of expected intersections
!> \param[out] n_src_dsts     offset and number intersections belonging
!>                            to source blocks
!> \param[in] src_sizes       sizes of source blocks
!> \param[in] dst_sizes       sizes of target blocks
! *****************************************************************************
  SUBROUTINE dbcsr_reblocking_targets (ints, numints, n_src_dsts,&
       src_sizes, dst_sizes)
    INTEGER, INTENT(INOUT)                   :: numints
    INTEGER, DIMENSION(2, numints), &
      INTENT(OUT)                            :: ints
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: n_src_dsts
    INTEGER, DIMENSION(:), INTENT(IN)        :: src_sizes, dst_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reblocking_targets', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: common_extent, current_dst, &
                                                current_int, current_src, &
                                                n_dst, n_src, s_dst, s_src

!   ---------------------------------------------------------------------------

    n_src = SIZE(src_sizes)
    n_dst = SIZE(dst_sizes)
    current_int = 1
    current_src = 1
    current_dst = 1
    IF(n_src.GT.0) THEN
       s_src = src_sizes (current_src)
       n_src_dsts(:,1) = (/ 1, 0 /)
    ELSE
       s_src = HUGE(0)
    ENDIF
    IF(n_dst.GT.0) THEN
       s_dst = dst_sizes (current_dst)
    ELSE
       s_dst = HUGE(0)
    ENDIF
    DO WHILE (current_src .LE. n_src .AND. current_dst .LE. n_dst)
       ! Calculate many elements do the current blocks have in common
       ! and record these as going to the current target block.
       common_extent = MIN (s_src, s_dst)
       ints(1, current_int) = current_dst ! target block
       ints(2, current_int) = common_extent
       ! We've used up the common extents.
       s_src = s_src - common_extent
       s_dst = s_dst - common_extent
       n_src_dsts(2,current_src) = n_src_dsts(2,current_src) + 1
       IF (s_src .EQ. 0) THEN
          current_src = current_src + 1
          IF (current_src .LE. n_src) THEN
             s_src = src_sizes (current_src)
             n_src_dsts(:, current_src) = (/ current_int+1, 0 /)
          ENDIF
       ENDIF
       IF (s_dst .EQ. 0) THEN
          current_dst = current_dst + 1
          IF (current_dst .LE. n_dst) s_dst = dst_sizes (current_dst)
       ENDIF
       current_int = current_int + 1
    ENDDO
    numints = current_int - 1
  END SUBROUTINE dbcsr_reblocking_targets




! *****************************************************************************
!> \brief Determines the sign for a block pointer.
!> \param[in] n          block offset
!> \param[in] oldsign    sign for the new block offset
!> \param[in] x          reverse the sign of the new block offset
! *****************************************************************************
  ELEMENTAL FUNCTION sgn (n, oldsign, x) RESULT (val)
    INTEGER, INTENT(IN)                      :: n, oldsign
    LOGICAL, INTENT(IN)                      :: x
    INTEGER                                  :: val

    IF (.NOT.x) THEN
       val = SIGN (n, oldsign)
    ELSE
       val = -SIGN (n, oldsign)
    ENDIF
  END FUNCTION sgn





! *****************************************************************************
!> \brief Fills meta information from a given distribution_2d
!> \param[out] meta           meta information array to fill
!> \param dist                processor distribution
!> \param[in] row_blk_size    row block sizes
!> \param[in] col_blk_size    column block sizes
! *****************************************************************************
  SUBROUTINE meta_from_dist (meta, dist, row_blk_size, col_blk_size)
    INTEGER, DIMENSION(dbcsr_meta_size), &
      INTENT(OUT)                            :: meta
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_from_dist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nfullcols_local, &
                                                nfullcols_total, &
                                                nfullrows_local, &
                                                nfullrows_total
    INTEGER, DIMENSION(:), POINTER           :: blkcols_local, blkrows_local, &
                                                cbs, rbs

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
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
    meta(5)  = dbcsr_distribution_nrows (dist)
    meta(6)  = dbcsr_distribution_ncols (dist)
    meta(7)  = nfullrows_total
    meta(8)  = nfullcols_total
    meta(9)  = dbcsr_distribution_nlocal_rows (dist)
    meta(10) = dbcsr_distribution_nlocal_cols (dist)
    meta(11) = nfullrows_local
    meta(12) = nfullcols_local
  END SUBROUTINE meta_from_dist


! *****************************************************************************
!> \brief Copies metadata into an array.
!> \param[in] matrix      Matrix
!> \param[out] meta       Metadata elements
! *****************************************************************************
  SUBROUTINE dbcsr_pack_meta(matrix, meta)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, DIMENSION(dbcsr_meta_size), &
      INTENT(OUT)                            :: meta

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_pack_meta', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    meta(dbcsr_slot_nblks)  = matrix%nblks
    meta(dbcsr_slot_nze)    = matrix%nze
    meta(dbcsr_slot_nblkrows_total ) = matrix%nblkrows_total
    meta(dbcsr_slot_nblkcols_total ) = matrix%nblkcols_total
    meta(dbcsr_slot_nfullrows_total) = matrix%nfullrows_total
    meta(dbcsr_slot_nfullcols_total) = matrix%nfullcols_total
    meta(dbcsr_slot_nblkrows_local ) = matrix%nblkrows_local
    meta(dbcsr_slot_nblkcols_local ) = matrix%nblkcols_local
    meta(dbcsr_slot_nfullrows_local) = matrix%nfullrows_local
    meta(dbcsr_slot_nfullcols_local) = matrix%nfullcols_local
    meta(dbcsr_slot_dense) = 0
    SELECT CASE (matrix%matrix_type)
    CASE ('N')
       meta(dbcsr_slot_type) = 0
    CASE ('T')
       meta(dbcsr_slot_type) = 1
    CASE ('S')
       meta(dbcsr_slot_type) = 2
    CASE ('A')
       meta(dbcsr_slot_type) = 3
    CASE DEFAULT
       meta(dbcsr_slot_type) = -1
       CALL dbcsr_assert(.FALSE.,  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, &
            "Invalid matrix type.",__LINE__,error)
    END SELECT
  END SUBROUTINE dbcsr_pack_meta

! *****************************************************************************
!> \brief Sets metadata form an array.
!> \param[in,out] matrix  Matrix
!> \param[in] meta        Metadata elements
! *****************************************************************************
  SUBROUTINE dbcsr_unpack_meta(matrix, meta)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, DIMENSION(dbcsr_meta_size), &
      INTENT(IN)                             :: meta

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_unpack_meta', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    matrix%nblks = meta(dbcsr_slot_nblks)
    matrix%nze   = meta(dbcsr_slot_nze)
    matrix%nblkrows_total  = meta(dbcsr_slot_nblkrows_total )
    matrix%nblkcols_total  = meta(dbcsr_slot_nblkcols_total )
    matrix%nfullrows_total = meta(dbcsr_slot_nfullrows_total)
    matrix%nfullcols_total = meta(dbcsr_slot_nfullcols_total)
    matrix%nblkrows_local  = meta(dbcsr_slot_nblkrows_local )
    matrix%nblkcols_local  = meta(dbcsr_slot_nblkcols_local )
    matrix%nfullrows_local = meta(dbcsr_slot_nfullrows_local)
    matrix%nfullcols_local = meta(dbcsr_slot_nfullcols_local)
    matrix%index(dbcsr_slot_dense) = 0
    SELECT CASE (meta(dbcsr_slot_type))
    CASE (0)
       matrix%matrix_type = 'N'
    CASE (1)
       matrix%matrix_type = 'T'
    CASE (2)
       matrix%matrix_type = 'S'
    CASE (3)
       matrix%matrix_type = 'A'
    CASE DEFAULT
       matrix%matrix_type = '0'
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid matrix type.",__LINE__,error)
    END SELECT
  END SUBROUTINE dbcsr_unpack_meta



  ELEMENTAL FUNCTION make_coordinate_tuple (most, least) RESULT (tuple)
    INTEGER, INTENT(IN)                      :: most, least
    INTEGER(KIND=int_8)                      :: tuple

!tuple = IOR (ISHFT (most, 32), least)

    tuple = most
    tuple = IOR (ISHFT(tuple, 32), INT(least, int_8))
  END FUNCTION make_coordinate_tuple



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
!> \brief Busy-loop sleep.
!> \param[in] duration        time to sleep/s.
! *****************************************************************************
  SUBROUTINE dbcsr_sleep(duration)
    REAL(KIND=dp), INTENT(IN)                :: duration

    REAL(KIND=dp)                            :: t_now, t_start

!   ---------------------------------------------------------------------------

    t_start = m_walltime ()
    t_now = t_start
    DO WHILE (t_now .LT. t_start + duration)
       t_now = m_walltime ()
    ENDDO
  END SUBROUTINE dbcsr_sleep


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






  ELEMENTAL FUNCTION gcd(a,b)
    INTEGER, INTENT(IN)                      :: a, b
    INTEGER                                  :: gcd

    INTEGER                                  :: aa, ab, l, rem, s

    aa=ABS(a)
    ab=ABS(b)
    IF (aa<ab) THEN
       s=aa
       l=ab
    ELSE
       s=ab
       l=aa
    ENDIF
    IF (s.NE.0) THEN
       DO
          rem=MOD(l,s)
          IF (rem==0) EXIT
          l=s
          s=rem
       ENDDO
       GCD=s
    ELSE
       GCD=l
    ENDIF
  END FUNCTION gcd

  ELEMENTAL FUNCTION lcm(a,b)
    INTEGER, INTENT(IN)                      :: a, b
    INTEGER                                  :: lcm

    INTEGER                                  :: tmp

    tmp = gcd(a,b)
    IF (tmp==0) THEN
       lcm = 0
    ELSE
       ! could still overflow if the true lcm is larger than maxint
       lcm = ABS((a/tmp)*b)
    END IF
  END FUNCTION lcm


! *****************************************************************************
!> \brief Calculates the checksum of a DBCSR matrix.
!> \param[in] matrix          matrix
!> \param[out] chksum         calculated checksum
!> \param[in,out] error       cp2k error
!> \param[in] local           (optional) no global communication
! *****************************************************************************
  FUNCTION dbcsr_checksum(matrix, local, error) RESULT(checksum)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: local
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    REAL(KIND=dp)                            :: checksum

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_checksum', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: c_sp
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: c_dp
    INTEGER                                  :: bc, blk, blk_p, br, &
                                                error_handler, m, mn, n
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: nocomm
    REAL(KIND=dp)                            :: blk_cs, local_cs
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r_sp
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r_dp

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index(matrix),&
         dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "Invalid matrix.",__LINE__,error)
    nocomm = .FALSE.
    IF (PRESENT (local)) nocomm = local
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    local_cs = 0.0_dp
    SELECT CASE (matrix%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix%m%data_area, r_sp)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix%m%data_area, c_dp)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix%m%data_area, c_sp)
    END SELECT
    DO br = 1, matrix%m%nblkrows_total
       m = row_blk_size(br)
       DO blk = matrix%m%row_p(br)+1, matrix%m%row_p(br+1)
          bc = matrix%m%col_i(blk)
          n = col_blk_size(bc)
          mn = m*n
          blk_p = ABS(matrix%m%blk_p(blk))
          IF (blk_p .NE. 0) THEN
             ! Calculate DDOT
             SELECT CASE (matrix%m%data_type)
             CASE (dbcsr_type_real_8)
                blk_cs = DOT_PRODUCT(r_dp(blk_p:blk_p+mn-1),&
                     &               r_dp(blk_p:blk_p+mn-1))
             CASE (dbcsr_type_real_4)
                blk_cs = DOT_PRODUCT(r_sp(blk_p:blk_p+mn-1),&
                     &               r_sp(blk_p:blk_p+mn-1))
             CASE (dbcsr_type_complex_8)
                blk_cs = DOT_PRODUCT(c_dp(blk_p:blk_p+mn-1),&
                     &               c_dp(blk_p:blk_p+mn-1))
             CASE (dbcsr_type_complex_4)
                blk_cs = DOT_PRODUCT(c_sp(blk_p:blk_p+mn-1),&
                     &               c_sp(blk_p:blk_p+mn-1))
             END SELECT
             DBGV ' chksum for (',br,',',bc,') at', blk_p,'l',mn,'= ', blk_cs
             local_cs = local_cs + blk_cs
          ENDIF
       ENDDO
    ENDDO
    DBG 'local chksum=',local_cs
    checksum = local_cs
    IF (.NOT. nocomm) THEN
       CALL mp_sum(local_cs, dbcsr_mp_group (dbcsr_distribution_mp (&
            matrix%m%dist)))
       checksum = local_cs
       DBGI 'local chksum=',local_cs,'global=',checksum
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END FUNCTION dbcsr_checksum

! *****************************************************************************
!> \brief Verify the correctness of a BCSR matrix.
!> \param[in] m     bcsr matrix
!> \param[in,out] error     cp2k error
!> \param[in] verbosity     how detailed errors are; 0=nothing; 1=summary at
!>                          end if matrix not consistent; 2=also individual
!>                          errors; 3=always print info about matrix; >3=even
!>                          more info
!> \param[in] local         (optional) no global communication
! *****************************************************************************
  SUBROUTINE dbcsr_verify_matrix(m, verbosity, local, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: m
    INTEGER, INTENT(IN), OPTIONAL            :: verbosity
    LOGICAL, INTENT(IN), OPTIONAL            :: local
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_verify_matrix', &
      r = moduleN//':'//routineN

    INTEGER :: bc, blk, blk_p, br, dbg, error_handler, i, mb, mn, n, &
      n_full_blocks_total, n_have_blocks_local, n_have_blocks_total, prev_br
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: nocomm
    REAL(KIND=dp)                            :: sparsity_total

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (m%m%initialized, 'EQ', dbcsr_magic_number,&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not verify non-existing matrix object.",__LINE__,error)
    dbg = 2
    nocomm = .FALSE.
    IF (PRESENT (local)) nocomm = local
    IF (PRESENT (verbosity)) dbg = verbosity
    IF (dbg.GE.3) WRITE(*,'(1X,A,A,A,A)')r//'Matrix name: ', m%m%name,&
         " of type ",m%m%matrix_type
    IF (dbg.GE.3) THEN
       WRITE(*,'(1X,A,I5,"x",I5,A,I5,"x",I5)')r//' Size blocked',&
            m%m%nblkrows_total, m%m%nblkcols_total, ", full ",&
            m%m%nfullrows_total, m%m%nfullcols_total
    ENDIF
    row_blk_size => array_data (m%m%row_blk_size)
    col_blk_size => array_data (m%m%col_blk_size)
    !
    IF(m%m%matrix_type.EQ.'N' .OR. m%m%matrix_type.EQ.'T') THEN
       n_full_blocks_total = m%m%nblkrows_total * m%m%nblkcols_total
    ELSE
       CALL dbcsr_assert(m%m%nblkrows_total, "EQ", m%m%nblkcols_total,&
            dbcsr_warning_level, dbcsr_internal_error, r,&
            'Symmetric matrix is not square',__LINE__,error)
       n_full_blocks_total = m%m%nblkrows_total * (m%m%nblkrows_total+1)/2
    ENDIF
    n_have_blocks_local = m%m%nblks
2045 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
2047 FORMAT (I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7)
    IF (dbg.GE.4) THEN
       WRITE(*,'(1X,A)')r//' index='
       WRITE(*,2045)m%m%index(:dbcsr_num_slots)
    ENDIF
    CALL dbcsr_assert(m%m%index(1), "GT", 0,&
         dbcsr_warning_level, dbcsr_internal_error, r,&
         'Index size 0',__LINE__,error)
    DO i = dbcsr_slot_row_p, dbcsr_num_slots
       CALL dbcsr_assert(m%m%index(i), "GT", 0,&
            dbcsr_warning_level, dbcsr_internal_error, r,&
            'Index member is 0',__LINE__,error)
       IF (.NOT. (i.EQ.dbcsr_slot_col_i.OR.i.EQ.dbcsr_slot_blk_p)) THEN
            CALL dbcsr_assert(m%m%index(i), "LE", m%m%index(1),&
                 dbcsr_warning_level, dbcsr_internal_error, r,&
                 'Index member is greater than size',__LINE__,error)
         ENDIF
    ENDDO
    !
    IF(dbg.GE.4) WRITE(*,*)r//' row_p extents',m%m%index(dbcsr_slot_row_p+1),&
         m%m%index(dbcsr_slot_row_p), SIZE(m%m%row_p)
    CALL dbcsr_assert(m%m%index(dbcsr_slot_row_p+1)-m%m%index(dbcsr_slot_row_p)+1,&
         "EQ", m%m%nblkrows_total+1,&
         dbcsr_warning_level, dbcsr_internal_error, r,&
         'Size of row_p index inconsistent with number of rows',__LINE__,error)
    CALL dbcsr_assert(SIZE(m%m%row_p), "EQ", m%m%nblkrows_total+1,&
         dbcsr_warning_level, dbcsr_internal_error, r,&
         'Size of row_p inconsistent with number of rows',__LINE__,error)
    !
    IF(dbg.GE.4) WRITE(*,*)r//' col_i extents',m%m%index(dbcsr_slot_col_i+1),&
         m%m%index(dbcsr_slot_col_i), SIZE(m%m%col_i)
    CALL dbcsr_assert(m%m%index(dbcsr_slot_col_i+1)-m%m%index(dbcsr_slot_col_i)+1,&
         "EQ", m%m%nblks,&
         dbcsr_warning_level, dbcsr_internal_error, r,&
         'Size of col_i index inconsistent with number of blocks',__LINE__,error)
    CALL dbcsr_assert(SIZE(m%m%col_i), "EQ", m%m%nblks,&
         dbcsr_warning_level, dbcsr_internal_error, r,&
         'Size of col inconsistent with number of blocks',__LINE__,error)
    !
    IF(dbg.GE.4) WRITE(*,*)r//' blk_p extents',m%m%index(dbcsr_slot_blk_p+1),&
         m%m%index(dbcsr_slot_blk_p), SIZE(m%m%blk_p)
    CALL dbcsr_assert(m%m%index(dbcsr_slot_blk_p+1)-m%m%index(dbcsr_slot_blk_p)+1,&
         "EQ", m%m%nblks,&
         dbcsr_warning_level, dbcsr_internal_error, r,&
         'Size of blk_p index inconsistent with number of blocks',__LINE__,error)
    CALL dbcsr_assert(SIZE(m%m%col_i), "EQ", m%m%nblks,&
         dbcsr_warning_level, dbcsr_internal_error, r,&
         'Size of blk_p inconsistent with number of blocks',__LINE__,error)
    !
    CALL dbcsr_assert(SIZE(row_blk_size), "EQ", m%m%nblkrows_total,&
         dbcsr_warning_level, dbcsr_internal_error, r,&
         'Row block size array inconsistent with number of blocked rows',__LINE__,error)
    CALL dbcsr_assert(SIZE(col_blk_size), "EQ", m%m%nblkcols_total,&
         dbcsr_warning_level, dbcsr_internal_error, r,&
         'Column block size array inconsistent with number of blocked columns',__LINE__,error)
    !
    IF(dbg.GE.4) THEN
       WRITE(*,'(1X,A,I7,A,I7)')r//' nze=',m%m%nze,'data size',&
            dbcsr_get_data_size(m%m%data_area)
    ENDIF
    CALL dbcsr_assert(dbcsr_get_data_size(m%m%data_area), "GE", m%m%nze,&
         dbcsr_warning_level, dbcsr_internal_error, r,&
         'Data storage may be too small.',__LINE__,error)
    IF (dbg.GE.5) THEN
       WRITE (*,'(1X,A,I7,A)') r//' size=',SIZE(m%m%row_p),' row_p='
       WRITE (*,2047) m%m%row_p(1:m%m%nblkrows_total+1)
       WRITE (*,'(1X,A)') r//' col_i='
       WRITE (*,2047) m%m%col_i(1:m%m%nblks)
       WRITE (*,'(1X,A)') r//' blk_p='
       WRITE (*,2047) m%m%blk_p(1:m%m%nblks)
    ENDIF
    prev_br = 0
    DO br = 1, m%m%nblkrows_total
       CALL dbcsr_assert(m%m%row_p(br), "GE", 0,&
            dbcsr_warning_level, dbcsr_internal_error, r,&
            'row_p less than zero',__LINE__,error)
       IF(br.GT.1) CALL dbcsr_assert(m%m%row_p(br), "GE", m%m%row_p(prev_br),&
            dbcsr_warning_level, dbcsr_internal_error, r,&
            'row_p decreases',__LINE__,error)
       mb = row_blk_size(br)
       CALL dbcsr_assert(mb, "GE", 0,&
            dbcsr_warning_level, dbcsr_internal_error, r,&
            'Row blocked size is negative',__LINE__,error)
       DO blk = m%m%row_p(br)+1, m%m%row_p(br+1)
          CALL dbcsr_assert(blk, "GT", 0,&
               dbcsr_warning_level, dbcsr_internal_error, r,&
               'Block number is zero',__LINE__,error)
          CALL dbcsr_assert(blk, "LE", m%m%nblks,&
               dbcsr_warning_level, dbcsr_internal_error, r,&
               'Block number too high',__LINE__,error)
          bc = m%m%col_i(blk)
          IF (dbg.GE.5) THEN
             WRITE(*,'(1X,A,I7,"(",I5,",",I5,")")')r//' block',blk,br,bc
          ENDIF
          CALL dbcsr_assert(bc, "GT", 0,&
               dbcsr_warning_level, dbcsr_internal_error, r,&
               'col_i is zero',__LINE__,error)
          CALL dbcsr_assert(bc, "LE", m%m%nblkcols_total,&
               dbcsr_warning_level, dbcsr_internal_error, r,&
               'col_i too high',__LINE__,error)
          n = col_blk_size(bc)
          CALL dbcsr_assert(n, "GE", 0,&
               dbcsr_warning_level, dbcsr_internal_error, r,&
               'Column blocked size is negative',__LINE__,error)
          blk_p = m%m%blk_p(blk)
          mn = mb*n
          !CALL dbcsr_assert(dbg,blk_p.GT.0, r, 'Block pointer is negative')
          CALL dbcsr_assert(blk_p, "LE", m%m%nze,&
               dbcsr_warning_level, dbcsr_internal_error, r,&
               'Block pointer too large',__LINE__,error)
          CALL dbcsr_assert(blk_p+mn-1, "LE", m%m%nze,&
               dbcsr_warning_level, dbcsr_internal_error, r,&
               'Block extends too far',__LINE__,error)
       ENDDO
       prev_br = br
    ENDDO
    IF (dbg.GE.3.AND..NOT.nocomm) THEN
       CALL mp_sum(n_have_blocks_local, dbcsr_mp_group (dbcsr_distribution_mp (&
            m%m%dist)))
       n_have_blocks_total = n_have_blocks_local
       sparsity_total = REAL(n_have_blocks_total)&
            / REAL(n_full_blocks_total)*100.0
       !WRITE(*,FMT='(30A,F5.1,A)')r//' Sparsity: ', sparsity_total,'%'
       WRITE(*,FMT='(1X,A,F5.1,A)')r//' Non-sparsity: ',&
            sparsity_total,'%'
    ENDIF

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_verify_matrix

END MODULE dbcsr_util
