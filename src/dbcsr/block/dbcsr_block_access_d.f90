!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Gets a 2-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-2 array)
!> \param[in] tr      whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_2d_block_d(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:,:), INTENT(OUT)     :: block
    LOGICAL, INTENT(IN)                      :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    TYPE(dbcsr_error_type)                   :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_2d_block_d', &
      routineP = moduleN//':'//routineN

    REAL(kind=real_8), DIMENSION(:), POINTER           :: block_1d
    INTEGER                                  :: rsize, csize,&
                                                blk, nze, offset,&
                                                stored_row,&
                                                stored_col, iw, nwms
    INTEGER                                  :: error_handle
    TYPE(btree_2d_data_d)          :: data_block
    LOGICAL                                  :: stored_tr

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set (routineN, error_handle, error=error)
    IF (debug_mod) THEN
       CALL dbcsr_assert (matrix%m%data_type, "EQ", dbcsr_type_real_8,&
            dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Data type mismatch for requested block.",__LINE__,error)
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    IF(found) THEN
       nze = rsize*csize
       IF (nze .GT. 0) THEN
          !
          ! let's copy the block
          block_1d => pointer_view (dbcsr_get_data_p (&
               matrix%m%data_area, 0.0_real_8), offset, offset+nze-1)
          IF (tr .EQV. stored_tr) THEN
             CALL dbcsr_block_copy (block(:,:), block_1d, rsize, csize)
          ELSE
             IF (stored_tr) CALL swap (rsize, csize)
             CALL dbcsr_block_transpose (block,&
                  block_1d,&
                  rsize, csize)
          ENDIF
       ENDIF
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",__LINE__,error)
       nwms = SIZE(matrix%m%wms)
       iw = 1
!$     CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$        dbcsr_fatal_level, dbcsr_internal_error,&
!$        routineN, "Number of work matrices not equal to number of threads", &
!$        __LINE__, error=error)
!$     iw = omp_get_thread_num () + 1
       IF (dbcsr_use_mutable (matrix%m)) THEN
          IF (.NOT. dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%m%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          CALL btree_get_d (&
               matrix%m%wms(iw)%mutable%m%btree_d,&
               make_coordinate_tuple(stored_row, stored_col),&
               data_block, found)
          IF (found) THEN
             IF (tr .EQV. data_block%tr) THEN
                CALL dcopy(nze,&
                     data_block%p(1,1), 1, block(1,1), 1)
             ELSE
                IF (data_block%tr) CALL swap (rsize, csize)
                block(:,:) = TRANSPOSE (data_block%p)
             ENDIF
          ENDIF
       ENDIF
    ENDIF
    IF (careful_mod) CALL dbcsr_error_stop (error_handle, error=error)
  END SUBROUTINE dbcsr_get_2d_block_d

! *****************************************************************************
!> \brief Gets a 2-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-2 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_2d_block_p_d(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    TYPE(dbcsr_error_type)                   :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_2d_block_p_d', &
      routineP = moduleN//':'//routineN

    REAL(kind=real_8), DIMENSION(:), POINTER           :: block_1d
    INTEGER                                  :: rsize, csize,&
                                                blk, nze, offset,&
                                                stored_row,&
                                                stored_col, iw, nwms
    INTEGER                                  :: error_handle
    TYPE(btree_2d_data_d)          :: data_block
    LOGICAL                                  :: stored_tr
    REAL(kind=real_8), DIMENSION(1,1), TARGET, SAVE    :: block0
!   ---------------------------------------------------------------------------
    IF (careful_mod) CALL dbcsr_error_set (routineN, error_handle, error=error)
    IF (debug_mod) THEN
       CALL dbcsr_assert (matrix%m%data_type, "EQ", dbcsr_type_real_8,&
            dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Data type mismatch for requested block.",__LINE__,error)
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)
    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    NULLIFY (block)
    IF(found) THEN
       nze = rsize*csize
       IF(nze.eq.0) THEN
          found = .TRUE.
          block => block0(1:0, 1:0)
       ELSE
          block_1d => pointer_view (dbcsr_get_data_p (&
               matrix%m%data_area, 0.0_real_8), offset, offset+nze-1)
          CALL dbcsr_set_block_pointer (matrix, block, stored_row, stored_col,&
               rsize, csize, stored_tr, offset, buffer_tr=stored_tr,&
               contiguous_pointers=.TRUE., error=error)
       ENDIF
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       nwms = SIZE(matrix%m%wms)
       iw = 1
!$     CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$        dbcsr_fatal_level, dbcsr_internal_error,&
!$        routineN, "Number of work matrices not equal to number of threads", &
!$        __LINE__, error=error)
!$     iw = omp_get_thread_num () + 1
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",__LINE__,error)
       IF (dbcsr_use_mutable (matrix%m)) THEN
          IF (.NOT. dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%m%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          CALL btree_get_d (&
               matrix%m%wms(iw)%mutable%m%btree_d,&
               make_coordinate_tuple(stored_row, stored_col),&
               data_block, found)
          IF (found) THEN
             block => data_block%p
          ENDIF
       ENDIF
    ENDIF
    IF (careful_mod) CALL dbcsr_error_stop (error_handle, error=error)
  END SUBROUTINE dbcsr_get_2d_block_p_d


! *****************************************************************************
!> \brief Gets a 1-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-1 array)
!> \param[in] tr      whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_d(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:), INTENT(OUT)       :: block
    LOGICAL, INTENT(IN)                      :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    TYPE(dbcsr_error_type)                   :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, csize, &
                                                nze, offset, &
                                                rsize, stored_row,&
                                                stored_col, iw, nwms
    LOGICAL                                  :: stored_tr
    TYPE(btree_2d_data_d)          :: data_block
    REAL(kind=real_8), DIMENSION(:), POINTER           :: block_1d

!   ---------------------------------------------------------------------------

    IF (debug_mod) THEN
       CALL dbcsr_assert (matrix%m%data_type, "EQ", dbcsr_type_real_8,&
            dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Data type mismatch for requested block.",__LINE__,error)
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    IF(found) THEN
       nze = rsize*csize
       IF (nze .GT. 0) THEN
          !
          ! let's copy the block
          block_1d => pointer_view (dbcsr_get_data_p (&
               matrix%m%data_area, 0.0_real_8), offset, offset+nze-1)
          IF (tr .EQV. stored_tr) THEN
             block(1:nze) = block_1d(1:nze)
          ELSE
             IF (stored_tr) CALL swap (rsize, csize)
             CALL dbcsr_block_transpose (block,&
                  block_1d,&
                  rsize, csize)
          ENDIF
       ENDIF
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       nwms = SIZE(matrix%m%wms)
       iw = 1
!$     IF (debug_mod) THEN
!$        CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$            dbcsr_fatal_level, dbcsr_internal_error,&
!$            routineN, "Number of work matrices not equal to number of threads", &
!$            __LINE__, error=error)
!$     ENDIF
!$     iw = omp_get_thread_num () + 1
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",__LINE__,error)
       IF (dbcsr_use_mutable (matrix%m)) THEN
          CALL btree_get_d (&
               matrix%m%wms(iw)%mutable%m%btree_d,&
               make_coordinate_tuple(stored_row, stored_col),&
               data_block, found)
          IF (found) THEN
             IF (tr .EQV. data_block%tr) THEN
                !CALL dcopy(nze,&
                !     data_block%p(1,1) 1, block(iw), 1)
                block(:) = RESHAPE (data_block%p(:,:), (/SIZE(data_block%p)/))
             ELSE
                IF (data_block%tr) CALL swap (rsize, csize)
                CALL dbcsr_block_transpose (block,&
                     data_block%p,&
                     rsize, csize)
             ENDIF
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_get_block_d

! *****************************************************************************
!> \brief Gets a 1-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-1 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_p_d(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:), POINTER           :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size
    TYPE(dbcsr_error_type)                   :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_p_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, csize, &
                                                nze, offset, &
                                                rsize, stored_row,&
                                                stored_col
    LOGICAL                                  :: stored_tr

!   ---------------------------------------------------------------------------

    IF (debug_mod) THEN
       CALL dbcsr_assert (matrix%m%data_type, "EQ", dbcsr_type_real_8,&
            dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Data type mismatch for requested block.",__LINE__,error)
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)
    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    NULLIFY (block)
    IF(found) THEN
       nze = rsize*csize
       !
       block => pointer_view (&
            dbcsr_get_data_p (matrix%m%data_area, 0.0_real_8), offset, offset+nze-1&
            )
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",__LINE__,error)
       CALL dbcsr_assert ("NOT", dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve rank-1 block pointers from mutable work matrices.",__LINE__,error)
    ENDIF
  END SUBROUTINE dbcsr_get_block_p_d


! *****************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix using the btree
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to reserve; added if not NULL
!> \param[in] transposed      the block holds transposed data
!> \param[out] existed        (optional) block already existed
! *****************************************************************************
  SUBROUTINE dbcsr_reserve_block2d_d(matrix, row, col, block,&
       transposed, existed)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_block2d_d', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_d)          :: data_block, data_block2
    INTEGER                                  :: col_size, row_size, &
                                                stored_row, stored_col, &
                                                iw, nwms
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found, gift, tr, sym_tr
    REAL(kind=real_8), DIMENSION(:,:), POINTER         :: original_block
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    gift = ASSOCIATED (block)
    IF (gift) THEN
       original_block => block
    ELSE
       NULLIFY (original_block)
    ENDIF
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    row_size = row_blk_size(row)
    col_size = col_blk_size(col)

    stored_row = row ; stored_col = col
    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    sym_tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix, stored_row, stored_col, sym_tr)
    !call cp_assert (associated (matrix%m%wms), cp_fatal_level,&
    !     cp_caller_error, routineN, "Work matrices not prepared")
    IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_create (matrix, work_mutable=.TRUE., error=error)
       !$OMP MASTER
       matrix%m%valid = .FALSE.
       !$OMP END MASTER
       !$OMP BARRIER
    ENDIF

    NULLIFY (data_block%p)
    IF (.NOT. gift) THEN
       ALLOCATE (data_block%p (row_size, col_size))
       block => data_block%p
    ELSE
       data_block%p => block
    ENDIF
    data_block%tr = tr

    nwms = SIZE(matrix%m%wms)
    iw = 1
!$  CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$     dbcsr_fatal_level, dbcsr_internal_error,&
!$     routineN, "Number of work matrices not equal to number of threads", &
!$     __LINE__, error=error)
!$  iw = omp_get_thread_num () + 1
    CALL btree_add_d (matrix%m%wms(iw)%mutable%m%btree_d,&
         make_coordinate_tuple(stored_row, stored_col),&
         data_block, found, data_block2)

    IF (.NOT. found) THEN
!$OMP CRITICAL (critical_reserve_block2d)
       matrix%m%valid = .FALSE.
!$OMP END CRITICAL (critical_reserve_block2d)
       matrix%m%wms(iw)%lastblk = matrix%m%wms(iw)%lastblk + 1
       matrix%m%wms(iw)%datasize = matrix%m%wms(iw)%datasize + row_size*col_size
    ELSE
       IF (.NOT. gift) THEN
          DEALLOCATE (data_block%p)
       ELSE
          DEALLOCATE (original_block)
       ENDIF
       !CALL cp_assert (ASSOCIATED (data_block%p), cp_warning_level,&
       !     cp_internal_error, routineN,&
       !     "Existing block has no associated pointer.")
       block => data_block2%p
    ENDIF
    IF (PRESENT (existed)) existed = found
  END SUBROUTINE dbcsr_reserve_block2d_d

! *****************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to put
!> \param[in]  transposed     the block is transposed
!> \param[in]  summation      (optional) if block exists, then sum the new
!>                            block to the old one instead of replacing it
!> \param[in]  scale          (optional) scale the block being added
! *****************************************************************************
  SUBROUTINE dbcsr_put_block2d_d(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:,:), INTENT(IN)      :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    REAL(kind=real_8), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block2d_d', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: tr, do_sum

    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    IF (PRESENT (summation)) THEN
       do_sum = summation
    ELSE
       do_sum = .FALSE.
    ENDIF
    IF (PRESENT (scale)) THEN
       CALL dbcsr_put_block (matrix, row, col,&
            RESHAPE (block, (/SIZE(block)/)), tr, do_sum, scale)
    ELSE
       CALL dbcsr_put_block (matrix, row, col,&
            RESHAPE (block, (/SIZE(block)/)), tr, do_sum)
    ENDIF
  END SUBROUTINE dbcsr_put_block2d_d

! *****************************************************************************
!> \brief Inserts a block in a dbcsr matrix.
!>
!> If the block exists, the current data is overwritten.
!> \param[in]  matrix         DBCSR matrix
!> \param[in]  row            the logical row
!> \param[in]  col            the logical column
!> \param[in]  block          the block to put
!> \param[in]  transposed     (optional) the block is transposed
!> \param[in]  summation      (optional) if block exists, then sum the new
!>                            block to the old one instead of replacing it
!> \param[in]  scale          (optional) scale the block being added
! *****************************************************************************
  SUBROUTINE dbcsr_put_block_d(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:), INTENT(IN)        :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    REAL(kind=real_8), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block_d', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_d)          :: data_block, data_block2
    INTEGER                                  :: blk, col_size, &
                                                nze, offset, &
                                                row_size, blk_p,&
                                                stored_row, stored_col,&
                                                iw, nwms
    LOGICAL                                  :: found, tr, do_sum, tr_diff,&
                                                sym_tr
    REAL(kind=real_8), DIMENSION(:), POINTER           :: block_1d
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------
    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    IF (PRESENT (summation)) THEN
       do_sum = summation
    ELSE
       do_sum = .FALSE.
    ENDIF
    row_size = dbcsr_blk_row_size(matrix, row)
    col_size = dbcsr_blk_column_size(matrix, col)
    IF (tr) CALL swap (row_size, col_size)

    stored_row = row ; stored_col = col; sym_tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix%m, stored_row, stored_col, sym_tr)
    nze = row_size*col_size
    !
    IF (debug_mod) THEN
       CALL dbcsr_assert (SIZE(block), "GE", nze, dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "Invalid block dimensions",__LINE__,error)
    ENDIF
    CALL dbcsr_get_stored_block_info (matrix%m, stored_row, stored_col,&
         found, blk, offset)
    IF(found) THEN
       ! let's copy the block
       offset = ABS (offset)
       ! Fix the index if the new block's transpose flag is different
       ! from the old one.
       tr_diff = .FALSE.
       IF (matrix%m%blk_p(blk).LT.0 .NEQV. tr) THEN
          tr_diff = .TRUE.
          matrix%m%blk_p(blk) = -matrix%m%blk_p(blk)
       ENDIF
       block_1d => pointer_view (dbcsr_get_data_p (&
            matrix%m%data_area, 0.0_real_8), offset, offset+nze-1)
       IF (nze .GT. 0) THEN
          IF (do_sum) THEN
             IF(tr_diff) &
                  block_1d = RESHAPE(TRANSPOSE(RESHAPE(block_1d,(/col_size,row_size/))),(/nze/))
             IF (PRESENT (scale)) THEN
                CALL daxpy (nze, scale, block(1:nze), 1,&
                     block_1d, 1)
             ELSE
                CALL daxpy (nze, 1.0_real_8, block(1:nze), 1,&
                     block_1d, 1)
             ENDIF
          ELSE
             IF (PRESENT (scale)) THEN
                CALL dcopy (nze, scale*block(1:nze), 1,&
                     block_1d, 1)
             ELSE
                CALL dcopy (nze, block(1:nze), 1,&
                     block_1d, 1)
             ENDIF
          ENDIF
       ENDIF
    ELSE
       !!@@@
       !call cp_assert (associated (matrix%m%wms), cp_fatal_level,&
       !     cp_caller_error, routineN, "Work matrices not prepared")
       IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
          CALL dbcsr_work_create (matrix, nblks_guess=1,&
               sizedata_guess=SIZE(block),error=error)
       ENDIF
       nwms = SIZE(matrix%m%wms)
       iw = 1
!$     IF (debug_mod) THEN
!$     CALL dbcsr_assert (nwms, "GE", omp_get_num_threads(),&
!$        dbcsr_fatal_level, dbcsr_internal_error,&
!$        routineN, "Number of work matrices not equal to number of threads", &
!$        __LINE__, error=error)
!$     ENDIF
!$     iw = omp_get_thread_num () + 1
       blk_p = matrix%m%wms(iw)%datasize + 1
       IF (.NOT.dbcsr_wm_use_mutable (matrix%m%wms(iw))) THEN
          IF (tr) blk_p = -blk_p
          CALL add_work_coordinate (matrix%m%wms(iw), row, col, blk_p, error=error)
          CALL dbcsr_data_ensure_size (matrix%m%wms(iw)%data_area,&
               matrix%m%wms(iw)%datasize+SIZE(block),&
               factor=default_resize_factor, error=error)
          IF (PRESENT (scale)) THEN
             CALL dbcsr_data_set (matrix%m%wms(iw)%data_area, ABS(blk_p),&
                  data_size=SIZE(block), src=scale*block, source_lb=1)
          ELSE
             CALL dbcsr_data_set (matrix%m%wms(iw)%data_area, ABS(blk_p),&
                  data_size=SIZE(block), src=block, source_lb=1)
          ENDIF
       ELSE
          ALLOCATE (data_block%p (row_size, col_size))
          IF (PRESENT (scale)) THEN
             data_block%p(:,:) = scale*RESHAPE (block, (/row_size, col_size/))
          ELSE
             data_block%p(:,:) = RESHAPE (block, (/row_size, col_size/))
          ENDIF
          data_block%tr = tr
          IF (.NOT. dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%m%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          IF (.NOT. do_sum) THEN
             CALL btree_add_d (&
                  matrix%m%wms(iw)%mutable%m%btree_d,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.TRUE.)
             IF (found) THEN
                CALL dbcsr_assert (ASSOCIATED (data_block2%p), dbcsr_warning_level,&
                     dbcsr_internal_error, routineN,&
                     "Data was not present in block",__LINE__,error)
                IF (ASSOCIATED (data_block2%p)) DEALLOCATE (data_block2%p)
             ENDIF
          ELSE
             CALL btree_add_d (&
                  matrix%m%wms(iw)%mutable%m%btree_d,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.FALSE.)
             IF (found .AND. nze .GT. 0) THEN
                CALL daxpy (nze, 1.0_real_8, block(1), 1,&
                     data_block2%p(1,1), 1)
                CALL dbcsr_assert (ASSOCIATED (data_block%p), dbcsr_warning_level,&
                     dbcsr_internal_error, routineN,&
                     "Data was not present in block",__LINE__,error)
                IF (ASSOCIATED (data_block%p)) DEALLOCATE (data_block%p)
             ENDIF
          ENDIF
          IF (.NOT. found) THEN
             matrix%m%wms(iw)%lastblk = matrix%m%wms(iw)%lastblk + 1
          ENDIF
       ENDIF
       IF (.NOT. found) THEN
          matrix%m%wms(iw)%datasize = matrix%m%wms(iw)%datasize + SIZE (block)
       ELSE
       ENDIF
!$OMP CRITICAL (dbcsr_put_block_critical)
       matrix%m%valid = .FALSE.
!$OMP END CRITICAL (dbcsr_put_block_critical)
    ENDIF
  END SUBROUTINE dbcsr_put_block_d


! *****************************************************************************
!> \brief Sets a pointer, possibly using the buffers.
!> \param[in] matrix           Matrix to use
!> \param pointer_any The pointer to set
!> \param row Row of block to point to
!> \param col Column of block to point to
!> \param rsize Row size of block to point to
!> \param csize Column size of block to point to
!> \param[in] main_tr          Whether block is transposed in the matrix
!> \param[in] base_offset      The block pointer
!> \param[in] buffer_tr        Whether buffer should be transposed
!> \param[in] contiguous_pointers  (optional) Whether pointers should be made
!>                                 contiguous
!> \param[in] read_only        (optional) User promise not to change data
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE dbcsr_set_block_pointer_2d_d (&
       matrix, pointer_any, row, col,&
       rsize, csize, main_tr, base_offset, buffer_tr, contiguous_pointers,&
       read_only, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    REAL(kind=real_8), DIMENSION(:,:), POINTER         :: pointer_any
    INTEGER, INTENT(IN)                      :: row, col, rsize, csize
    LOGICAL, INTENT(IN)                      :: main_tr
    INTEGER, INTENT(IN)                      :: base_offset
    LOGICAL, INTENT(IN)                      :: buffer_tr
    LOGICAL, INTENT(IN), OPTIONAL            :: contiguous_pointers, read_only
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_set_block_pointer_2d_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    REAL(kind=real_8), DIMENSION(:), POINTER           :: lin_blk_p

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set (routineN, error_handler, error)
    CALL dbcsr_get_data (matrix%m%data_area, lin_blk_p,&
         lb=base_offset, ub=base_offset+rsize*csize-1)
    CALL pointer_d_rank_remap2 (pointer_any, rsize, csize,&
         lin_blk_p)
    IF (careful_mod) CALL dbcsr_error_stop (error_handler, error)
  END SUBROUTINE dbcsr_set_block_pointer_2d_d
