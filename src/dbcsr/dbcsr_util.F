!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE dbcsr_data_methods,              ONLY: dbcsr_data_get_size,&
                                             dbcsr_get_data
  USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_local_cols,&
                                             dbcsr_distribution_local_rows,&
                                             dbcsr_distribution_mp,&
                                             dbcsr_distribution_ncols,&
                                             dbcsr_distribution_nlocal_cols,&
                                             dbcsr_distribution_nlocal_rows,&
                                             dbcsr_distribution_nrows
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_failure_level, dbcsr_fatal_level, &
       dbcsr_internal_error, dbcsr_wrong_args_error
  USE dbcsr_methods,                   ONLY: dbcsr_blk_col_offset,&
                                             dbcsr_blk_row_offset,&
                                             dbcsr_get_data_size_referenced,&
                                             dbcsr_has_symmetry,&
                                             dbcsr_is_initialized,&
                                             dbcsr_valid_index
  USE dbcsr_mp_methods,                ONLY: dbcsr_mp_group,&
                                             dbcsr_mp_mypcol,&
                                             dbcsr_mp_myprow
  USE dbcsr_toollib,                   ONLY: sort,&
                                             swap
  USE dbcsr_types,                     ONLY: &
       dbcsr_distribution_obj, dbcsr_meta_size, dbcsr_num_slots, dbcsr_obj, &
       dbcsr_slot_blk_p, dbcsr_slot_col_i, dbcsr_slot_dense, &
       dbcsr_slot_home_coli, dbcsr_slot_home_pcol, dbcsr_slot_home_prow, &
       dbcsr_slot_home_rowi, dbcsr_slot_home_vpcol, dbcsr_slot_home_vprow, &
       dbcsr_slot_nblkcols_local, dbcsr_slot_nblkcols_total, &
       dbcsr_slot_nblkrows_local, dbcsr_slot_nblkrows_total, &
       dbcsr_slot_nblks, dbcsr_slot_nfullcols_local, &
       dbcsr_slot_nfullcols_total, dbcsr_slot_nfullrows_local, &
       dbcsr_slot_nfullrows_total, dbcsr_slot_nze, dbcsr_slot_row_p, &
       dbcsr_slot_type, dbcsr_type, dbcsr_type_complex_4, &
       dbcsr_type_complex_8, dbcsr_type_real_4, dbcsr_type_real_8
  USE kinds,                           ONLY: dp,&
                                             int_4,&
                                             real_4,&
                                             real_8
  USE message_passing,                 ONLY: mp_sum

  !$ USE OMP_LIB
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_util'

  ! Main
  PUBLIC :: dbcsr_checksum, dbcsr_verify_matrix,&
            dbcsr_pack_meta, dbcsr_unpack_meta, meta_from_dist
  ! Block sizes and arrays
  PUBLIC :: convert_sizes_to_offsets, convert_offsets_to_sizes,&
            global_offsets_to_local,&
            nfull_elements,&
            dbcsr_calc_block_sizes,&
            find_block_of_element,&
            get_internal_offsets,&
            map_most_common
  ! utility routines
  PUBLIC :: count_bins, sgn


  LOGICAL, PARAMETER :: bcsr_info =    .FALSE.
  LOGICAL, PARAMETER :: bcsr_verbose = .FALSE.

CONTAINS


! *****************************************************************************
!> \brief Finds block to which a full element belongs.
!> \param[in] full            full element
!> \param[out] block          block to which full belongs
!> \param nblocks ...
!> \param[in] block_offsets   starting full elements of blocks
!> \param[in] hint            where to start looking; ignored if 0
!> \param error ...
!> \par Assumptions
!>        It is assumed that block_start and block_end are sorted and
!>        that hint is in the range [0, nblocks].
! *****************************************************************************
  SUBROUTINE find_block_of_element(full, block, nblocks,&
       block_offsets, hint, error)
    INTEGER, INTENT(in)                      :: full
    INTEGER, INTENT(out)                     :: block
    INTEGER, INTENT(in)                      :: nblocks
    INTEGER, DIMENSION(1:nblocks+1), &
      INTENT(in)                             :: block_offsets
    INTEGER, INTENT(in)                      :: hint
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: count

    IF (hint .NE. 0) THEN
       block = hint
    ELSE
       block = MAX(1,(nblocks+1)/2)
    ENDIF
    count = 0
    DO WHILE (block_offsets(block) .GT. full .OR. block_offsets(block+1)-1 .LT. full)
       IF (block_offsets(block) .GT. full) THEN
          block = block-1
       ELSEIF (block_offsets(block+1)-1 .LT. full) THEN
          block = block+1
       ENDIF
       count = count+1
       IF (dbg) THEN
          IF (count .GT. nblocks .OR. block .LT. 1 .OR. block .GT. nblocks) THEN
             WRITE(*,'(1X,A,I9,A,I9,A)')"Want to find block",&
                  block," of",nblocks," blocks"
             CALL dbcsr_assert (count .LE. nblocks, dbcsr_failure_level,&
                  dbcsr_internal_error, "find_block_of_element",&
                  "Too many searches",__LINE__,error)
          ENDIF
       ENDIF
    ENDDO
  END SUBROUTINE find_block_of_element



! *****************************************************************************
!> \brief The sum of a subset of rows/columns
!> \param[in] all_offsets     ordered offsets of all the elements
!> \param[in] local_elements  enumerated local elements
!> \retval nfull_elements     sum of sizes of local elemetns 
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
    ELSE
       IF (.NOT. PRESENT (offsets_stop)) THEN
          offsets_start(1) = 0
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
!> \brief Converts global offsets to local
!> \param[in] global_offsets   Offsets of elements in the global grid
!> \param[in] local_elements   Which elements are local
!> \param[out] local_offsets   Offsets of local elements.
!> \par Global vs. Local Indexing
!>      local_offsets may be sized according to the
!>      local index (|local_elements+|1) or the
!>      global index (|global_offsets|).
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
!> \brief Finds internal offsets
!>
!> For all local blocks in blk_local_els, it calculates its offset in
!> the dense block to which it belongs.
!> \param blk_local_els ...
!> \param el_map ...
!> \param blk_el_offsets ...
!> \param dense_el_offsets ...
!> \param internal_offsets ...
! *****************************************************************************
  SUBROUTINE get_internal_offsets(blk_local_els, el_map, blk_el_offsets,&
       dense_el_offsets, internal_offsets)
    INTEGER, DIMENSION(:), INTENT(IN)        :: blk_local_els, el_map, &
                                                blk_el_offsets, &
                                                dense_el_offsets
    INTEGER, DIMENSION(:), INTENT(OUT)       :: internal_offsets

    CHARACTER(len=*), PARAMETER :: routineN = 'get_internal_offsets', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_el, d_el, i, ndense, nlblk
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: off_acc

!   ---------------------------------------------------------------------------

    nlblk = SIZE(blk_local_els)
    ndense = SIZE(dense_el_offsets)
    ALLOCATE (off_acc (ndense))
    off_acc(:) = 0
    internal_offsets(:) = 0
    DO i = 1, nlblk
       blk_el = blk_local_els(i)
       d_el = el_map(blk_el)
       internal_offsets(blk_el) = off_acc(d_el)
       off_acc(d_el) = off_acc(d_el) + blk_el_offsets(blk_el+1) - blk_el_offsets(blk_el)
    ENDDO
    DEALLOCATE (off_acc)
  END SUBROUTINE get_internal_offsets


! *****************************************************************************
!> \brief Copies row and column block sizes from another matrix.
!> \param[in,out] dst      target matrix 
!> \param[in] src          source matrix 
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
!> \param[in] row_p index structure
!> \param[in] col_i index structure
!> \param[in] rbs row block sizes
!> \param[in] cbs column block sizes
! *****************************************************************************
  SUBROUTINE dbcsr_calc_block_sizes(sizes, row_p, col_i, rbs, cbs)
    INTEGER, DIMENSION(*), INTENT(OUT)       :: sizes
    INTEGER, DIMENSION(:), INTENT(IN)        :: row_p
    INTEGER, DIMENSION(*), INTENT(IN)        :: col_i, rbs, cbs

    INTEGER                                  :: blk, nrows, row, row_size

    nrows = SIZE(row_p)-1
    !$OMP DO
    DO row = 1, nrows
       row_size = rbs(row)
       DO blk = row_p(row)+1, row_p(row+1)
          sizes(blk) = row_size * cbs(col_i(blk))
       END DO
    ENDDO
    !$OMP END DO
  END SUBROUTINE dbcsr_calc_block_sizes


! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param oldsign ...
!> \param x ...
!> \retval val ...
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

!   ---------------------------------------------------------------------------

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
    meta(dbcsr_slot_home_prow) = dbcsr_mp_myprow (dbcsr_distribution_mp (dist))
    meta(dbcsr_slot_home_rowi) = 1
    meta(dbcsr_slot_home_pcol) = dbcsr_mp_mypcol (dbcsr_distribution_mp (dist))
    meta(dbcsr_slot_home_coli) = 1
    meta(dbcsr_slot_home_vprow) = -1
    meta(dbcsr_slot_home_vpcol) = -1
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
    meta(dbcsr_slot_type) = 0
    !IF (matrix%transpose)&
    !     meta(dbcsr_slot_type) = IBSET (meta(dbcsr_slot_type), 0)
    IF (matrix%symmetry)&
         meta(dbcsr_slot_type) = IBSET (meta(dbcsr_slot_type), 1)
    IF (matrix%negate_real)&
         meta(dbcsr_slot_type) = IBSET (meta(dbcsr_slot_type), 2)
    IF (matrix%negate_imaginary)&
         meta(dbcsr_slot_type) = IBSET (meta(dbcsr_slot_type), 3)
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
    !matrix%transpose = BTEST (meta(dbcsr_slot_type), 0)
    matrix%symmetry = BTEST (meta(dbcsr_slot_type), 1)
    matrix%negate_real = BTEST (meta(dbcsr_slot_type), 2)
    matrix%negate_imaginary = BTEST (meta(dbcsr_slot_type), 3)
  END SUBROUTINE dbcsr_unpack_meta

! *****************************************************************************
!> \brief Calculates the checksum of a DBCSR matrix.
!> \param[in] matrix          matrix
!> \param[in] local           (optional) no global communication
!> \param[in] pos             (optional) position-dependent checksum
!> \param[in,out] error       cp2k error
!> \retval checksum calculated checksum
! *****************************************************************************
  FUNCTION dbcsr_checksum(matrix, local, pos, error) RESULT(checksum)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: local, pos
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    REAL(KIND=dp)                            :: checksum

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_checksum', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: c_sp
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: c_dp
    INTEGER                                  :: bc, blk, blk_p, br, co, &
                                                error_handler, m, mn, n, ro
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: nocomm, pd, tr
    REAL(KIND=dp)                            :: blk_cs, local_cs, local_cs_row
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r_sp
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r_dp

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index(matrix),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Invalid matrix.",__LINE__,error)
    nocomm = .FALSE.
    IF (PRESENT (local)) nocomm = local
    IF (PRESENT (pos)) THEN
       pd = pos
    ELSE
       pd = .FALSE.
    ENDIF
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
       ro = dbcsr_blk_row_offset (matrix, br)
       local_cs_row=0
       !$OMP PARALLEL DO DEFAULT(NONE) &
       !$OMP          PRIVATE(bc,m,n,mn,blk_p,blk_cs,tr,co) &
       !$OMP          SHARED(pd,br,matrix,ro,row_blk_size,col_blk_size,r_dp, r_sp, c_dp,c_sp) &
       !$OMP           REDUCTION(+:local_cs_row)
       DO blk = matrix%m%row_p(br)+1, matrix%m%row_p(br+1)
          bc = matrix%m%col_i(blk)
          m = row_blk_size(br)
          n = col_blk_size(bc)
          mn = m*n
          blk_p = ABS(matrix%m%blk_p(blk))
          tr = matrix%m%blk_p(blk) .LT. 0
          IF (blk_p .NE. 0) THEN
             IF (mn .GT. 0) THEN
                IF (tr) CALL swap (m, n)
                co = dbcsr_blk_col_offset (matrix, bc)
                ! Calculate DDOT
                SELECT CASE (matrix%m%data_type)
                CASE (dbcsr_type_real_8)
                   IF (pd) THEN
                      blk_cs = pd_blk_cs (m, n, r_dp(blk_p:blk_p+mn-1),&
                           tr, ro, co)
                   ELSE
                      blk_cs = REAL(DOT_PRODUCT(r_dp(blk_p:blk_p+mn-1),&
                                                r_dp(blk_p:blk_p+mn-1)), KIND=dp)
                   ENDIF
                CASE (dbcsr_type_real_4)
                   IF (pd) THEN
                      blk_cs = pd_blk_cs (m,n, REAL(r_sp(blk_p:blk_p+mn-1), KIND=dp),&
                           tr, ro, co)
                   ELSE
                      blk_cs = REAL(DOT_PRODUCT(r_sp(blk_p:blk_p+mn-1),&
                                                r_sp(blk_p:blk_p+mn-1)), KIND=dp)
                   ENDIF
                CASE (dbcsr_type_complex_8)
                   IF (pd) THEN
                      blk_cs = pd_blk_cs (m, n, REAL(c_dp(blk_p:blk_p+mn-1), KIND=dp),&
                           tr, ro, co)
                   ELSE
                      blk_cs = REAL(DOT_PRODUCT(c_dp(blk_p:blk_p+mn-1),&
                                                c_dp(blk_p:blk_p+mn-1)), KIND=dp)
                   ENDIF
                CASE (dbcsr_type_complex_4)
                   IF (pd) THEN
                      blk_cs = pd_blk_cs (m, n, REAL(c_sp(blk_p:blk_p+mn-1), KIND=dp),&
                           tr, ro, co)
                   ELSE
                      blk_cs = REAL(DOT_PRODUCT(c_sp(blk_p:blk_p+mn-1),&
                                                c_sp(blk_p:blk_p+mn-1)), KIND=dp)
                   ENDIF
                CASE default
                   blk_cs = 0.0_dp
                END SELECT
             ELSE
                blk_cs = 0.0_dp
             ENDIF
             local_cs_row = local_cs_row + blk_cs
          ENDIF
       ENDDO
       local_cs=local_cs+local_cs_row
    ENDDO
    checksum = local_cs
    IF (.NOT. nocomm) THEN
       CALL mp_sum(local_cs, dbcsr_mp_group (dbcsr_distribution_mp (&
            matrix%m%dist)))
       checksum = local_cs
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END FUNCTION dbcsr_checksum

! *****************************************************************************
!> \brief ...
!> \param ld ...
!> \param od ...
!> \param DATA ...
!> \param tr ...
!> \param ro ...
!> \param co ...
!> \retval pd_cs ...
! *****************************************************************************
  PURE FUNCTION pd_blk_cs (ld, od, DATA, tr, ro, co) RESULT (pd_cs)
    INTEGER, INTENT(IN)                      :: ld, od
    REAL(KIND=dp), DIMENSION(ld, od), &
      INTENT(IN)                             :: DATA
    LOGICAL, INTENT(IN)                      :: tr
    INTEGER, INTENT(IN)                      :: ro, co
    REAL(KIND=dp)                            :: pd_cs

    INTEGER                                  :: c, cs, r, rs

    pd_cs = 0.0_dp
    rs = ld ; cs = od
    IF (tr) THEN
       CALL swap (rs, cs)
       DO r = 1, rs
          DO c = 1, cs
             pd_cs = pd_cs + DATA(c, r)*LOG(ABS(REAL((ro+r-1),KIND=dp)*REAL((co+c-1), KIND=dp)))
          ENDDO
       ENDDO
    ELSE
       DO c = 1, cs
          DO r = 1, rs
             pd_cs = pd_cs + DATA(r, c)*LOG(ABS(REAL((ro+r-1),KIND=dp)*REAL((co+c-1), KIND=dp)))
          ENDDO
       ENDDO
    ENDIF
  END FUNCTION pd_blk_cs


! *****************************************************************************
!> \brief Verify the correctness of a BCSR matrix.
!> \param[in] m     bcsr matrix
!> \param[in] verbosity     how detailed errors are; 0=nothing; 1=summary at
!>                          end if matrix not consistent; 2=also individual
!>                          errors; 3=always print info about matrix; >3=even
!>                          more info
!> \param[in] local         (optional) no global communication
!> \param[in,out] error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_verify_matrix(m, verbosity, local, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: m
    INTEGER, INTENT(IN), OPTIONAL            :: verbosity
    LOGICAL, INTENT(IN), OPTIONAL            :: local
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_verify_matrix', &
      r = moduleN//':'//routineN

    INTEGER :: bc, blk, blk_p, br, data_size_referenced, dbg, error_handler, &
      i, mb, mn, n, n_full_blocks_total, n_have_blocks_local, &
      n_have_blocks_total, prev_br
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: nocomm
    REAL(KIND=dp)                            :: sparsity_total

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_is_initialized(m),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not verify non-existing matrix object.",__LINE__,error)
    dbg = 2
    nocomm = .FALSE.
    IF (PRESENT (local)) nocomm = local
    IF (PRESENT (verbosity)) dbg = verbosity
    IF (dbg.GE.3) WRITE(*,'(1X,A,A,A,3(L1))')r//'Matrix name: ', m%m%name,&
         " of types ", m%m%symmetry, m%m%negate_real,&
         m%m%negate_imaginary
    IF (dbg.GE.3) THEN
       WRITE(*,'(1X,A,I5,"x",I5,A,I5,"x",I5)')r//' Size blocked',&
            m%m%nblkrows_total, m%m%nblkcols_total, ", full ",&
            m%m%nfullrows_total, m%m%nfullcols_total
    ENDIF
    row_blk_size => array_data (m%m%row_blk_size)
    col_blk_size => array_data (m%m%col_blk_size)
    !
    IF(.NOT.dbcsr_has_symmetry (m)) THEN
       n_full_blocks_total = m%m%nblkrows_total * m%m%nblkcols_total
    ELSE
       CALL dbcsr_assert(m%m%nblkrows_total, "EQ", m%m%nblkcols_total,&
            dbcsr_fatal_level, dbcsr_internal_error, r,&
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
         dbcsr_fatal_level, dbcsr_internal_error, r,&
         'Index size 0',__LINE__,error)
    DO i = dbcsr_slot_row_p, dbcsr_num_slots
       !CALL dbcsr_assert(m%m%index(i), "GT", 0,&
       !     dbcsr_fatal_level, dbcsr_internal_error, r,&
       !     'Index member is 0',__LINE__,error)
       IF (.NOT. (i.EQ.dbcsr_slot_col_i.OR.i.EQ.dbcsr_slot_blk_p)) THEN
            CALL dbcsr_assert(m%m%index(i), "LE", m%m%index(1),&
                 dbcsr_fatal_level, dbcsr_internal_error, r,&
                 'Index member is greater than size',__LINE__,error)
         ENDIF
    ENDDO
    !
    IF(dbg.GE.4) WRITE(*,*)r//' row_p extents',m%m%index(dbcsr_slot_row_p+1),&
         m%m%index(dbcsr_slot_row_p), SIZE(m%m%row_p)
    CALL dbcsr_assert(m%m%index(dbcsr_slot_row_p+1)-m%m%index(dbcsr_slot_row_p)+1,&
         "EQ", m%m%nblkrows_total+1,&
         dbcsr_fatal_level, dbcsr_internal_error, r,&
         'Size of row_p index inconsistent with number of rows',__LINE__,error)
    CALL dbcsr_assert(SIZE(m%m%row_p), "EQ", m%m%nblkrows_total+1,&
         dbcsr_fatal_level, dbcsr_internal_error, r,&
         'Size of row_p inconsistent with number of rows',__LINE__,error)
    !
    IF(dbg.GE.4) WRITE(*,*)r//' col_i extents',m%m%index(dbcsr_slot_col_i+1),&
         m%m%index(dbcsr_slot_col_i), SIZE(m%m%col_i)
    CALL dbcsr_assert(m%m%index(dbcsr_slot_col_i+1)-m%m%index(dbcsr_slot_col_i)+1,&
         "EQ", m%m%nblks,&
         dbcsr_fatal_level, dbcsr_internal_error, r,&
         'Size of col_i index inconsistent with number of blocks',__LINE__,error)
    CALL dbcsr_assert(SIZE(m%m%col_i), "EQ", m%m%nblks,&
         dbcsr_fatal_level, dbcsr_internal_error, r,&
         'Size of col inconsistent with number of blocks',__LINE__,error)
    !
    IF(dbg.GE.4) WRITE(*,*)r//' blk_p extents',m%m%index(dbcsr_slot_blk_p+1),&
         m%m%index(dbcsr_slot_blk_p), SIZE(m%m%blk_p)
    CALL dbcsr_assert(m%m%index(dbcsr_slot_blk_p+1)-m%m%index(dbcsr_slot_blk_p)+1,&
         "EQ", m%m%nblks,&
         dbcsr_fatal_level, dbcsr_internal_error, r,&
         'Size of blk_p index inconsistent with number of blocks',__LINE__,error)
    CALL dbcsr_assert(SIZE(m%m%col_i), "EQ", m%m%nblks,&
         dbcsr_fatal_level, dbcsr_internal_error, r,&
         'Size of blk_p inconsistent with number of blocks',__LINE__,error)
    !
    CALL dbcsr_assert(SIZE(row_blk_size), "EQ", m%m%nblkrows_total,&
         dbcsr_fatal_level, dbcsr_internal_error, r,&
         'Row block size array inconsistent with number of blocked rows',__LINE__,error)
    CALL dbcsr_assert(SIZE(col_blk_size), "EQ", m%m%nblkcols_total,&
         dbcsr_fatal_level, dbcsr_internal_error, r,&
         'Column block size array inconsistent with number of blocked columns',__LINE__,error)
    !
    IF(dbg.GE.4) THEN
       WRITE(*,'(1X,A,I7,A,I7)')r//' nze=',m%m%nze,'data size',&
            dbcsr_data_get_size(m%m%data_area)
    ENDIF
    data_size_referenced = dbcsr_get_data_size_referenced (m)
    !This tends to be too verbose and usually untrue for symmetric
    !matrices.
    !CALL dbcsr_assert(dbcsr_get_data_size(m%m%data_area), "GE", m%m%nze,&
    !     dbcsr_fatal_level, dbcsr_internal_error, r,&
    !     'Data storage may be too small.',__LINE__,error)
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
            dbcsr_fatal_level, dbcsr_internal_error, r,&
            'row_p less than zero',__LINE__,error)
       IF(br.GT.1) CALL dbcsr_assert(m%m%row_p(br), "GE", m%m%row_p(prev_br),&
            dbcsr_fatal_level, dbcsr_internal_error, r,&
            'row_p decreases',__LINE__,error)
       mb = row_blk_size(br)
       CALL dbcsr_assert(mb, "GE", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, r,&
            'Row blocked size is negative',__LINE__,error)
       DO blk = m%m%row_p(br)+1, m%m%row_p(br+1)
          CALL dbcsr_assert(blk, "GT", 0,&
               dbcsr_fatal_level, dbcsr_internal_error, r,&
               'Block number is zero',__LINE__,error)
          CALL dbcsr_assert(blk, "LE", m%m%nblks,&
               dbcsr_fatal_level, dbcsr_internal_error, r,&
               'Block number too high',__LINE__,error)
          bc = m%m%col_i(blk)
          IF (dbg.GE.5) THEN
             WRITE(*,'(1X,A,I7,"(",I5,",",I5,")")')r//' block',blk,br,bc
          ENDIF
          CALL dbcsr_assert(bc, "GT", 0,&
               dbcsr_fatal_level, dbcsr_internal_error, r,&
               'col_i is zero',__LINE__,error)
          CALL dbcsr_assert(bc, "LE", m%m%nblkcols_total,&
               dbcsr_fatal_level, dbcsr_internal_error, r,&
               'col_i too high',__LINE__,error)
          n = col_blk_size(bc)
          CALL dbcsr_assert(n, "GE", 0,&
               dbcsr_fatal_level, dbcsr_internal_error, r,&
               'Column blocked size is negative',__LINE__,error)
          blk_p = m%m%blk_p(blk)
          mn = mb*n
          !CALL dbcsr_assert(dbg,blk_p.GT.0, r, 'Block pointer is negative')
          !CALL dbcsr_assert(blk_p, "LE", m%m%nze,&
          !     dbcsr_fatal_level, dbcsr_internal_error, r,&
          !     'Block pointer too large',__LINE__,error)
          !CALL dbcsr_assert(blk_p+mn-1, "LE", m%m%nze,&
          !     dbcsr_fatal_level, dbcsr_internal_error, r,&
          !     'Block extends too far',__LINE__,error)
          IF (mn .GT. 0) THEN
             CALL dbcsr_assert (ABS(blk_p), "LE", data_size_referenced,&
                  dbcsr_fatal_level, dbcsr_internal_error, r,&
                  "Block pointer pointso outside of declared referenced area",&
                  __LINE__, error)
          ENDIF
          CALL dbcsr_assert (ABS(blk_p)+mn-1, "LE", data_size_referenced,&
               dbcsr_fatal_level, dbcsr_internal_error, r,&
               "Block extends outside of declared referenced area",&
               __LINE__, error)
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


! *****************************************************************************
!> \brief ...
!> \param nelements ...
!> \param bins ...
!> \param nbins ...
!> \param bin_counts ...
! *****************************************************************************
  PURE SUBROUTINE count_bins(nelements, bins, nbins, bin_counts)
    INTEGER, INTENT(IN)                      :: nelements
    INTEGER, DIMENSION(1:nelements), &
      INTENT(IN)                             :: bins
    INTEGER, INTENT(IN)                      :: nbins
    INTEGER, DIMENSION(1:nbins), INTENT(OUT) :: bin_counts

    INTEGER                                  :: el

    bin_counts(:) = 0
    DO el = 1, nelements
       bin_counts(bins(el)) = bin_counts(bins(el)) + 1
    ENDDO
  END SUBROUTINE count_bins

! *****************************************************************************
!> \brief Makes a lookup table from the most common elements.
!> \param[in] array                  Array for which to find the most common
!>                                   elements.
!> \param[out] most_common_map       Ranking of the most common elements
!>                                   in array
!> \param[in] nmost_common           The number of most common elements
!> \param[out] most_common_elements  The most common elements in array
!> \param[in] size_limit             Limit maximum size to this value
!> \param max_val ...
!> \par Lookup table
!>      The lookup table is indexed by the most common array values
!>      (i.e., block sizes).  Its values are the order of their frequency.
! *****************************************************************************
  SUBROUTINE map_most_common(array, most_common_map, nmost_common,&
       most_common_elements, size_limit, max_val)
    INTEGER, DIMENSION(:), INTENT(IN)        :: array
    INTEGER(KIND=int_4), DIMENSION(:), &
      POINTER                                :: most_common_map
    INTEGER, INTENT(IN)                      :: nmost_common
    INTEGER, DIMENSION(:), INTENT(OUT)       :: most_common_elements
    INTEGER, INTENT(IN)                      :: size_limit
    INTEGER, INTENT(OUT)                     :: max_val

    CHARACTER(LEN=*), PARAMETER :: routineN = 'map_most_common', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, max_val_l, nmc
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: permutation, size_counts

!   ---------------------------------------------------------------------------

    IF (SIZE (array) .GT. 0) THEN
       max_val = MAXVAL (array)
       max_val_l = MIN(MIN (size_limit, max_val), INT(HUGE(most_common_map)))
    ELSE
       max_val = 0
       max_val_l = 0
    ENDIF
    ! Count the frequency of all block sizes up to max_val_l.
    ALLOCATE (size_counts (0:max_val_l))
    ALLOCATE (permutation (0:max_val_l))
    size_counts = 0
    permutation = 0
    DO i = 1, SIZE(array)
       ! Counts are decreased to easily get a reverse sort order.
       IF (array(i) .LE. max_val_l) &
            size_counts(array(i)) = size_counts(array(i)) - 1
    END DO
    IF (SIZE(array) .GT. 0) THEN
       CALL sort (size_counts, max_val_l+1, permutation)
    ENDIF
    ! Limiting nmc to max_val_l prevents out-of-bounds.
    nmc = MIN(nmost_common, max_val_l)
    ! Determine the biggest block size and allocate the map.
    ALLOCATE (most_common_map(0:max_val_l))
    ! Create the mapping from block size to order.
    most_common_map = nmost_common+1
    DO i = 1, nmc
       most_common_map(permutation(i-1)-1) = i
    END DO
    ! Copy the most common elements
    most_common_elements(:) = 0
    most_common_elements(1:nmc) = permutation(0:nmc-1)-1
    DEALLOCATE (size_counts)
    DEALLOCATE (permutation)
  END SUBROUTINE map_most_common


END MODULE dbcsr_util
