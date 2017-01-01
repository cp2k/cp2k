!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Gets a 2-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-2 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! **************************************************************************************************
  SUBROUTINE dbcsr_get_2d_block_p_c(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_2d_block_p_c', &
      routineP = moduleN//':'//routineN

    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: block_1d
    INTEGER                                  :: rsize, csize,&
                                                blk, nze, offset,&
                                                stored_row,&
                                                stored_col, iw, nwms
    INTEGER                                  :: error_handle
    TYPE(btree_2d_data_c)          :: data_block
    LOGICAL                                  :: stored_tr
    COMPLEX(kind=real_4), DIMENSION(1,1), TARGET, SAVE    :: block0
!   ---------------------------------------------------------------------------
    IF (careful_mod) CALL timeset (routineN, error_handle)
    IF (debug_mod) THEN
       IF(matrix%data_type /= dbcsr_type_complex_4) &
          CPABORT("Data type mismatch for requested block.")
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)
    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix, stored_row)
    csize = dbcsr_blk_column_size (matrix, stored_col)
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
               matrix%data_area, CMPLX(0.0, 0.0, real_4)), offset, offset+nze-1)
          CALL dbcsr_set_block_pointer (matrix, block, rsize, csize, offset)
       ENDIF
    ELSEIF (ASSOCIATED (matrix%wms)) THEN
       nwms = SIZE(matrix%wms)
       iw = 1
!$     IF(nwms < omp_get_num_threads()) &
!$        CPABORT("Number of work matrices not equal to number of threads")
!$     iw = omp_get_thread_num () + 1
       IF(.NOT.dbcsr_use_mutable (matrix))&
          CPABORT("Can not retrieve blocks from non-mutable work matrices.")
       IF (dbcsr_use_mutable (matrix)) THEN
          IF (.NOT. dbcsr_mutable_instantiated(matrix%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          CALL btree_get_c (&
               matrix%wms(iw)%mutable%m%btree_c,&
               make_coordinate_tuple(stored_row, stored_col),&
               data_block, found)
          IF (found) THEN
             block => data_block%p
          ENDIF
       ENDIF
    ENDIF
    IF (careful_mod) CALL timestop (error_handle)
  END SUBROUTINE dbcsr_get_2d_block_p_c


! **************************************************************************************************
!> \brief Gets a 1-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-1 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! **************************************************************************************************
  SUBROUTINE dbcsr_get_block_p_c(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_type), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_p_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, csize, &
                                                nze, offset, &
                                                rsize, stored_row,&
                                                stored_col
    LOGICAL                                  :: stored_tr

!   ---------------------------------------------------------------------------

    IF (debug_mod) THEN
       IF(matrix%data_type /= dbcsr_type_complex_4) &
          CPABORT("Data type mismatch for requested block.")
    ENDIF

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)
    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix, stored_row)
    csize = dbcsr_blk_column_size (matrix, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    NULLIFY (block)
    IF(found) THEN
       nze = rsize*csize
       !
       block => pointer_view (&
            dbcsr_get_data_p (matrix%data_area, CMPLX(0.0, 0.0, real_4)), offset, offset+nze-1&
            )
    ELSEIF (ASSOCIATED (matrix%wms)) THEN
       IF(.NOT.dbcsr_use_mutable (matrix))&
          CPABORT("Can not retrieve blocks from non-mutable work matrices.")
       IF(dbcsr_use_mutable (matrix))&
          CPABORT("Can not retrieve rank-1 block pointers from mutable work matrices.")
    ENDIF
  END SUBROUTINE dbcsr_get_block_p_c


! **************************************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix using the btree
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to reserve; added if not NULL
!> \param[in] transposed      the block holds transposed data
!> \param[out] existed        (optional) block already existed
! **************************************************************************************************
  SUBROUTINE dbcsr_reserve_block2d_c(matrix, row, col, block,&
       transposed, existed)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_block2d_c', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_c)          :: data_block, data_block2
    INTEGER                                  :: col_size, row_size, &
                                                stored_row, stored_col, &
                                                iw, nwms
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found, gift, tr, sym_tr
    COMPLEX(kind=real_4), DIMENSION(:,:), POINTER         :: original_block

!   ---------------------------------------------------------------------------

    gift = ASSOCIATED (block)
    IF (gift) THEN
       original_block => block
    ELSE
       NULLIFY (original_block)
    ENDIF
    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)
    row_size = row_blk_size(row)
    col_size = col_blk_size(col)

    stored_row = row ; stored_col = col
    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    sym_tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix, stored_row, stored_col)
    IF (.NOT.ASSOCIATED (matrix%wms)) THEN
       CALL dbcsr_work_create (matrix, work_mutable=.TRUE.)
       !$OMP MASTER
       matrix%valid = .FALSE.
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

    nwms = SIZE(matrix%wms)
    iw = 1
!$  IF(nwms < omp_get_num_threads()) &
!$     CPABORT("Number of work matrices not equal to number of threads")
!$  iw = omp_get_thread_num () + 1
    CALL btree_add_c (matrix%wms(iw)%mutable%m%btree_c,&
         make_coordinate_tuple(stored_row, stored_col),&
         data_block, found, data_block2)

    IF (.NOT. found) THEN
!$OMP CRITICAL (critical_reserve_block2d)
       matrix%valid = .FALSE.
!$OMP END CRITICAL (critical_reserve_block2d)
       matrix%wms(iw)%lastblk = matrix%wms(iw)%lastblk + 1
       matrix%wms(iw)%datasize = matrix%wms(iw)%datasize + row_size*col_size
    ELSE
       IF (.NOT. gift) THEN
          DEALLOCATE (data_block%p)
       ELSE
          DEALLOCATE (original_block)
       ENDIF
       block => data_block2%p
    ENDIF
    IF (PRESENT (existed)) existed = found
  END SUBROUTINE dbcsr_reserve_block2d_c

! **************************************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to put
!> \param[in]  transposed     the block is transposed
!> \param[in]  summation      (optional) if block exists, then sum the new
!>                            block to the old one instead of replacing it
!> \param[in]  scale          (optional) scale the block being added
! **************************************************************************************************
  SUBROUTINE dbcsr_put_block2d_c(matrix, row, col, block, lb_row_col, transposed,&
       summation, flop, scale)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:,:), INTENT(IN)      :: block
    INTEGER, DIMENSION(2), OPTIONAL, INTENT(INOUT) :: lb_row_col
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    INTEGER(KIND=int_8), INTENT(INOUT), OPTIONAL :: flop
    COMPLEX(kind=real_4), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block2d_c', &
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
            RESHAPE (block, (/SIZE(block)/)), lb_row_col, tr, do_sum, flop, scale)
    ELSE
       CALL dbcsr_put_block (matrix, row, col,&
            RESHAPE (block, (/SIZE(block)/)), lb_row_col, tr, do_sum, flop)
    ENDIF
  END SUBROUTINE dbcsr_put_block2d_c

! **************************************************************************************************
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
!> \param[in]  scale          (optional) scale the OBblock being added
! **************************************************************************************************
  SUBROUTINE dbcsr_put_block_c(matrix, row, col, block, lb_row_col, transposed,&
       summation, flop, scale)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:), INTENT(IN)        :: block
    INTEGER, DIMENSION(2), OPTIONAL, INTENT(INOUT) :: lb_row_col
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    INTEGER(KIND=int_8), INTENT(INOUT), OPTIONAL :: flop
    COMPLEX(kind=real_4), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block_c', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_c)          :: data_block, data_block2
    INTEGER                                  :: blk, col_size, &
                                                nze, offset, &
                                                row_size, blk_p,&
                                                stored_row, stored_col,&
                                                iw, nwms
    LOGICAL                                  :: found, tr, do_sum, tr_diff,&
                                                sym_tr
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: block_1d
    INTEGER(KIND=int_8)                      :: my_flop

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
    my_flop = 0
    row_size = dbcsr_blk_row_size(matrix, row)
    col_size = dbcsr_blk_column_size(matrix, col)
    IF (tr) CALL swap (row_size, col_size)

    stored_row = row ; stored_col = col; sym_tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix, stored_row, stored_col)
    nze = row_size*col_size
    !
    IF (debug_mod .AND. SIZE(block) < nze) &
       CPABORT("Invalid block dimensions")
    CALL dbcsr_get_stored_block_info (matrix, stored_row, stored_col,&
         found, blk, lb_row_col, offset)
    IF(found) THEN
       ! let's copy the block
       offset = ABS (offset)
       ! Fix the index if the new block's transpose flag is different
       ! from the old one.
       tr_diff = .FALSE.
       IF (matrix%blk_p(blk).LT.0 .NEQV. tr) THEN
          tr_diff = .TRUE.
          matrix%blk_p(blk) = -matrix%blk_p(blk)
       ENDIF
       block_1d => pointer_view (dbcsr_get_data_p (&
            matrix%data_area, CMPLX(0.0, 0.0, real_4)), offset, offset+nze-1)
       IF (nze .GT. 0) THEN
          IF (do_sum) THEN
             IF(tr_diff) &
                  block_1d = RESHAPE(TRANSPOSE(RESHAPE(block_1d,(/col_size,row_size/))),(/nze/))
             IF (PRESENT (scale)) THEN
                CALL caxpy (nze, scale, block(1:nze), 1,&
                     block_1d, 1)
             ELSE
                CALL caxpy (nze, CMPLX(1.0, 0.0, real_4), block(1:nze), 1,&
                     block_1d, 1)
             ENDIF
             my_flop = my_flop + nze * 2
          ELSE
             IF (PRESENT (scale)) THEN
                CALL ccopy (nze, scale*block(1:nze), 1,&
                     block_1d, 1)
             ELSE
                CALL ccopy (nze, block(1:nze), 1,&
                     block_1d, 1)
             ENDIF
          ENDIF
       ENDIF
    ELSE
       !!@@@
       !call cp_assert (associated (matrix%wms), cp_fatal_level,&
       !     cp_caller_error, routineN, "Work matrices not prepared")
       IF (.NOT.ASSOCIATED (matrix%wms)) THEN
          CALL dbcsr_work_create (matrix, nblks_guess=1,&
               sizedata_guess=SIZE(block))
       ENDIF
       nwms = SIZE(matrix%wms)
       iw = 1
!$     IF(debug_mod .AND. nwms < omp_get_num_threads()) &
!$        CPABORT("Number of work matrices not equal to number of threads")
!$     iw = omp_get_thread_num () + 1
       blk_p = matrix%wms(iw)%datasize + 1
       IF (.NOT.dbcsr_wm_use_mutable (matrix%wms(iw))) THEN
          IF (tr) blk_p = -blk_p
          CALL add_work_coordinate (matrix%wms(iw), row, col, blk_p)
          CALL dbcsr_data_ensure_size (matrix%wms(iw)%data_area,&
               matrix%wms(iw)%datasize+SIZE(block),&
               factor=default_resize_factor)
          IF (PRESENT (scale)) THEN
             CALL dbcsr_data_set (matrix%wms(iw)%data_area, ABS(blk_p),&
                  data_size=SIZE(block), src=scale*block, source_lb=1)
          ELSE
             CALL dbcsr_data_set (matrix%wms(iw)%data_area, ABS(blk_p),&
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
          IF (.NOT. dbcsr_mutable_instantiated(matrix%wms(iw)%mutable)) THEN
             CALL dbcsr_mutable_new(matrix%wms(iw)%mutable,&
                  dbcsr_get_data_type(matrix))
          ENDIF
          IF (.NOT. do_sum) THEN
             CALL btree_add_c (&
                  matrix%wms(iw)%mutable%m%btree_c,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.TRUE.)
             IF (found) THEN
                IF(.NOT.ASSOCIATED(data_block2%p))&
                   CPWARN("Data was not present in block")
                IF (ASSOCIATED (data_block2%p)) DEALLOCATE (data_block2%p)
             ENDIF
          ELSE
             CALL btree_add_c (&
                  matrix%wms(iw)%mutable%m%btree_c,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.FALSE.)
             IF (found) THEN
                IF(nze > 0) &
                   CALL caxpy (nze, CMPLX(1.0, 0.0, real_4), block(1), 1,&
                        data_block2%p(1,1), 1)
                IF(.NOT.ASSOCIATED(data_block%p))&
                   CPWARN("Data was not present in block")
                IF (ASSOCIATED (data_block%p)) DEALLOCATE (data_block%p)
             ENDIF
          ENDIF
          IF (.NOT. found) THEN
             matrix%wms(iw)%lastblk = matrix%wms(iw)%lastblk + 1
          ENDIF
       ENDIF
       IF (.NOT. found) THEN
          matrix%wms(iw)%datasize = matrix%wms(iw)%datasize + SIZE (block)
       ELSE
       ENDIF
!$OMP CRITICAL (dbcsr_put_block_critical)
       matrix%valid = .FALSE.
!$OMP END CRITICAL (dbcsr_put_block_critical)
    ENDIF
    IF (PRESENT(flop)) flop = flop + my_flop
  END SUBROUTINE dbcsr_put_block_c


! **************************************************************************************************
!> \brief Sets a pointer, possibly using the buffers.
!> \param[in] matrix           Matrix to use
!> \param pointer_any The pointer to set
!> \param rsize Row size of block to point to
!> \param csize Column size of block to point to
!> \param[in] base_offset      The block pointer
! **************************************************************************************************
  SUBROUTINE dbcsr_set_block_pointer_2d_c (&
       matrix, pointer_any, rsize, csize, base_offset)
    TYPE(dbcsr_type), INTENT(IN)              :: matrix
    COMPLEX(kind=real_4), DIMENSION(:,:), POINTER         :: pointer_any
    INTEGER, INTENT(IN)                      :: rsize, csize
    INTEGER, INTENT(IN)                      :: base_offset

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_set_block_pointer_2d_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: lin_blk_p

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL timeset (routineN, error_handler)
    CALL dbcsr_get_data (matrix%data_area, lin_blk_p,&
         lb=base_offset, ub=base_offset+rsize*csize-1)
    CALL pointer_rank_remap2 (pointer_any, rsize, csize,&
         lin_blk_p)
    IF (careful_mod) CALL timestop (error_handler)
  END SUBROUTINE dbcsr_set_block_pointer_2d_c
