! *****************************************************************************
!> \brief traces a DBCSR matrix
!> \param[in] matrix_a       DBCSR matrix
!> \param[out] trace         the trace of the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_trace_a_d(matrix_a, trace, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(kind=real_8), INTENT(INOUT)                   :: trace
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_a_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_row, a_row_size, i, &
                                                mynode, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size,&
                                                row_dist, col_dist
    REAL(kind=real_8), DIMENSION(:), POINTER           :: a_data, data_p
    INTEGER, DIMENSION(:,:), POINTER         :: pgrid
    TYPE(dbcsr_distribution_obj)             :: dist

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (matrix_a%m%row_blk_size)
    col_blk_size => array_data (matrix_a%m%col_blk_size)
    CALL dbcsr_get_data (matrix_a%m%data_area, data_p)
    dist = dbcsr_distribution (matrix_a)
    mynode = dbcsr_mp_mynode (dbcsr_distribution_mp (dist))
    pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (dist))
    row_dist => array_data (dbcsr_distribution_row_dist (dist))
    col_dist => array_data (dbcsr_distribution_col_dist (dist))
    !
    ! let's go
    trace = 0.0_dp
    DO a_row = 1, matrix_a%m%nblkrows_total
       a_row_size = row_blk_size(a_row)
       DO a_blk = matrix_a%m%row_p(a_row)+1,matrix_a%m%row_p(a_row+1)
          a_col = matrix_a%m%col_i(a_blk)
          IF(a_col.ne.a_row) CYCLE
          ! We must skip non-local blocks in a replicated matrix.
          IF(matrix_a%m%replication_type .NE. dbcsr_repl_full) THEN
             IF (mynode .NE. checker_square_proc (a_row, a_col, pgrid,&
                  row_dist, col_dist)) CYCLE
          ENDIF
          a_col_size = col_blk_size(a_col)
          IF(a_row_size.NE.a_col_size) CALL stop_program(routineN,"is that a square matrix?")
          a_nze = a_row_size**2
          a_data => pointer_view (data_p, ABS(matrix_a%m%blk_p(a_blk)),&
               ABS(matrix_a%m%blk_p(a_blk))+a_nze-1)
          !data_a => matrix_a%m%data(ABS(matrix_a%m%blk_p(a_blk)):ABS(matrix_a%m%blk_p(a_blk))+a_nze-1)
          !
          ! let's trace the block
          DO i = 1,a_row_size
             trace = trace + a_data((i-1)*a_row_size+i)
          ENDDO
       ENDDO ! a_col
    ENDDO ! a_row
    !
    ! summe
    CALL mp_sum(trace,dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%m%dist)))

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_trace_a_d

! *****************************************************************************
!> \brief traces a product of DBCSR matrices
!> \param[in] matrix_a,matrix_b  DBCSR matrices
!> \param[out] trace             the trace of the product of the matrices
!> \param[in] trans_a            (optional) is matrix_a transposed or not?
!> \param[in] trans_b            (optional) is matrix_b transposed or not?
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_trace_ab_d(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a, matrix_b
    REAL(kind=real_8), INTENT(INOUT)                   :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    LOGICAL, INTENT(IN), OPTIONAL            :: local_sum
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_ab_d', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: my_trans_a, my_trans_b
    INTEGER :: a_blk, a_col, a_col_size, a_row_size, b_blk, b_col_size, &
      b_frst_blk, b_last_blk, b_row_size, nze, row, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, &
                                                a_row_blk_size, &
                                                b_col_blk_size, b_row_blk_size
    REAL(kind=real_8)                                  :: sym_fac, fac
    LOGICAL                                  :: found, my_local_sum
    REAL(kind=real_8), EXTERNAL                        :: DDOT
    REAL(kind=real_8), DIMENSION(:), POINTER           :: a_data_p, b_data_p

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    my_local_sum = .FALSE.
    IF(PRESENT(local_sum)) my_local_sum=local_sum
    IF(.NOT.my_local_sum) THEN!fixme
       CALL cp_assert (matrix_a%m%replication_type .EQ. dbcsr_repl_none&
            .AND. matrix_b%m%replication_type .EQ. dbcsr_repl_none,&
            cp_failure_level, cp_unimplemented_error_nr, routineN,&
            "Trace of product of replicated matrices not yet possible.")
    ENDIF

    sym_fac = 1.0_dp
    IF((dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric.OR.&
        dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric).AND.&
       (dbcsr_get_matrix_type(matrix_b).EQ.dbcsr_type_symmetric.OR.&
        dbcsr_get_matrix_type(matrix_b).EQ.dbcsr_type_antisymmetric)) sym_fac = 2.0_dp

    a_row_blk_size => array_data (matrix_a%m%row_blk_size)
    a_col_blk_size => array_data (matrix_a%m%col_blk_size)
    b_row_blk_size => array_data (matrix_b%m%row_blk_size)
    b_col_blk_size => array_data (matrix_b%m%col_blk_size)
    CALL dbcsr_get_data (matrix_a%m%data_area, a_data_p)
    CALL dbcsr_get_data (matrix_b%m%data_area, b_data_p)
    
    my_trans_a = 'T'
    IF(PRESENT(trans_a)) my_trans_a = trans_a
    my_trans_b = 'N'
    IF(PRESENT(trans_b)) my_trans_b = trans_b
    IF(my_trans_a.NE.'T'.OR.my_trans_b.NE.'N') CALL stop_program(routineN,"NYI")
    !
    ! let's go
    trace = 0.0_dp
    IF(matrix_a%m%nblkrows_total.NE.matrix_b%m%nblkrows_total) THEN
       CALL stop_program(routineN,"matrices not consistent")
    ENDIF
    DO row = 1, matrix_a%m%nblkrows_total
       a_row_size = a_row_blk_size(row)
       b_row_size = b_row_blk_size(row)
       IF(a_row_size.NE.b_row_size) THEN
          WRITE(*,*) 'a_row_size',a_row_size
          WRITE(*,*) 'b_row_size',b_row_size
          CALL stop_program(routineN,"matrices not consistent")
       ENDIF
       b_blk = matrix_b%m%row_p(row)+1
       b_frst_blk = matrix_b%m%row_p(row)+1
       b_last_blk = matrix_b%m%row_p(row+1)
       DO a_blk = matrix_a%m%row_p(row)+1,matrix_a%m%row_p(row+1)
          IF (matrix_a%m%blk_p(a_blk) .EQ. 0) CYCLE ! Deleted block
          a_col = matrix_a%m%col_i(a_blk)
          a_col_size = a_col_blk_size(a_col)
          !
          ! find the b_blk we assume here that the colums are ordered !
          CALL dbcsr_find_column(a_col,b_frst_blk,b_last_blk,matrix_b%m%col_i,&
               matrix_b%m%blk_p,b_blk,found)
          IF(found) THEN
             b_col_size = b_col_blk_size(a_col)
             IF(a_col_size.NE.b_col_size)  THEN
                WRITE(*,*) 'a_col_size',a_col_size
                WRITE(*,*) 'b_col_size',b_col_size
                CALL stop_program(routineN,"matrices not consistent")
             ENDIF
             !
             nze = a_row_size*a_col_size
             !
             IF(nze.GT.0) THEN
                !
                ! let's trace the blocks
                fac = 1.0_dp
                IF(row.NE.a_col) fac = sym_fac
                trace = trace + fac * DDOT (nze,&
                        a_data_p(ABS(matrix_a%m%blk_p(a_blk))),1,&
                        b_data_p(ABS(matrix_b%m%blk_p(b_blk))),1)
             ENDIF
          ENDIF
       ENDDO ! a_col
    ENDDO ! a_row
    !
    ! summe
    IF(.NOT.my_local_sum) &
         & CALL mp_sum(trace,dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%m%dist)))

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_trace_ab_d


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

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_2d_block_d', &
      routineP = moduleN//':'//routineN
    
    REAL(kind=real_8), DIMENSION(:), POINTER           :: block_1d
    INTEGER                                  :: rsize, csize,&
                                                blk, last_blk, nze, offset,&
                                                frst_blk, stored_row,&
                                                stored_col, iw
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    TYPE(btree_2d_data_d)          :: data_block
    LOGICAL                                  :: stored_tr
!   ---------------------------------------------------------------------------
    CALL cp_assert (matrix%m%data_type .EQ. dbcsr_type_real_8,&
         cp_fatal_level, cp_caller_error,&
         routineN, "Data type mismatch for requested block.")

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
       CALL cp_assert (dbcsr_use_mutable (matrix%m), cp_failure_level,&
            cp_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.")
       iw = 1
!$     iw = omp_get_thread_num ()+1
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

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_2d_block_p_d', &
      routineP = moduleN//':'//routineN
    
    REAL(kind=real_8), DIMENSION(:), POINTER           :: block_1d
    INTEGER                                  :: rsize, csize,&
                                                blk, last_blk, nze, offset,&
                                                frst_blk, stored_row,&
                                                stored_col, iw
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    TYPE(btree_2d_data_d)          :: data_block
    LOGICAL                                  :: stored_tr
!   ---------------------------------------------------------------------------
    CALL cp_assert (matrix%m%data_type .EQ. dbcsr_type_real_8,&
         cp_fatal_level, cp_caller_error,&
         routineN, "Data type mismatch for requested block.")

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
          found = .FALSE.
          NULLIFY (block)
       ELSE
          block_1d => pointer_view (dbcsr_get_data_p (&
               matrix%m%data_area, 0.0_real_8), offset, offset+nze-1)
#if defined(__PTR_RANK_REMAP)
          block(1:rsize, 1:csize) => block_1d
#elif !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
          CALL pointer_d_rank_remap2 (block, rsize, csize, block_1d)
#endif
       ENDIF
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       iw = 1
!$     iw = omp_get_thread_num()+1
       CALL cp_assert (dbcsr_use_mutable (matrix%m), cp_failure_level,&
            cp_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.")
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
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:), INTENT(OUT)       :: block
    LOGICAL, INTENT(IN)                      :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size


    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, csize, frst_blk, &
                                                last_blk, nze, offset, &
                                                rsize, stored_row,&
                                                stored_col, iw
    LOGICAL                                  :: stored_tr
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    TYPE(btree_2d_data_d)          :: data_block
    REAL(kind=real_8), DIMENSION(:), POINTER           :: block_1d
!   ---------------------------------------------------------------------------
    CALL cp_assert (matrix%m%data_type .EQ. dbcsr_type_real_8,&
         cp_fatal_level, cp_caller_error,&
         routineN, "Data type mismatch for requested block.")

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
       iw = 1
!$     iw = omp_get_thread_num()+1
       CALL cp_assert (dbcsr_use_mutable (matrix%m), cp_failure_level,&
            cp_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.")
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
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:), POINTER           :: block
    LOGICAL, INTENT(OUT)                     :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_p_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, csize, frst_blk, &
                                                last_blk, nze, offset, &
                                                rsize, stored_row,&
                                                stored_col
    LOGICAL                                  :: stored_tr
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    TYPE(btree_2d_data_d)          :: data_block
!   ---------------------------------------------------------------------------
    CALL cp_assert (matrix%m%data_type .EQ. dbcsr_type_real_8,&
         cp_fatal_level, cp_caller_error,&
         routineN, "Data type mismatch for requested block.")

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
            dbcsr_get_data_p (matrix%m%data_area, 0.0_real_8), offset, offset&
            )
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       CALL cp_assert (dbcsr_use_mutable (matrix%m), cp_failure_level,&
            cp_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.")
       CALL cp_assert (.NOT.dbcsr_use_mutable (matrix%m), cp_failure_level,&
            cp_caller_error, routineN,&
            "Can not retrieve rank-1 block pointers from mutable work matrices.")
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
    INTEGER                                  :: blk, col_size, frst_blk, &
                                                last_blk, nze, offset, &
                                                row_size, blk_p,&
                                                stored_row, stored_col, i,&
                                                iw
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found, gift, tr, tr2
    REAL(kind=real_8), DIMENSION(:,:), POINTER         :: original_block

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
    CALL dbcsr_get_stored_coordinates (matrix, stored_row, stored_col, tr)
    !@@@
    !call cp_assert (associated (matrix%m%wms), cp_fatal_level,&
    !     cp_caller_error, routineN, "Work matrices not prepared")
    IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_create (matrix, work_mutable=.TRUE.)
       matrix%m%valid = .FALSE.
    ENDIF

    NULLIFY (data_block%p)
    IF (.NOT. gift) THEN
       ALLOCATE (data_block%p (row_size, col_size))
       block => data_block%p
    ELSE
       data_block%p => block
    ENDIF
    data_block%tr = tr

    iw = 1
!$  iw = omp_get_thread_num()+1
    CALL btree_add_d (matrix%m%wms(iw)%mutable%m%btree_d,&
         make_coordinate_tuple(stored_row, stored_col),&
         data_block, found, data_block2)

    IF (.NOT. found) THEN
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
!> \brief Get a 2D block from the working matrix
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          pointer to data block
!> \param[out] tr             block is transposed
!> \param[out] found          block was found
! *****************************************************************************
  SUBROUTINE dbcsr_get_tmp_block2d_d(matrix, row, col, block, tr,&
       found)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(OUT)                     :: tr, found

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_tmp_block2d_d', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_d)          :: data_block
    INTEGER                                  :: stored_row, stored_col, iw

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
       found = .FALSE.
       NULLIFY (block)
       RETURN
    ENDIF

    stored_row = row ; stored_col = col
    tr = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix, stored_row, stored_col, tr)
    iw = 1
!$  iw = omp_get_thread_num()+1
    CALL btree_get_d (matrix%m%wms(iw)%mutable%m%btree_d,&
         make_coordinate_tuple(stored_row, stored_col),&
         data_block, found)

    IF (.NOT. found) THEN
       !WRITE(*,*)routineN//" Block did not exist"
       NULLIFY (block)
       !block => data_block%p
    ELSE
       !WRITE(*,*)routineN//" Found block", ASSOCIATED (data_block%p)
       block => data_block%p
       tr = data_block%tr
    ENDIF
  END SUBROUTINE dbcsr_get_tmp_block2d_d



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
    INTEGER                                  :: blk, col_size, frst_blk, &
                                                last_blk, nze, offset, &
                                                row_size, blk_p,&
                                                stored_row, stored_col,&
                                                iw
    LOGICAL                                  :: found, tr, do_sum
    REAL(kind=real_8), DIMENSION(:), POINTER           :: block_1d

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

    stored_row = row ; stored_col = col
    CALL dbcsr_get_stored_coordinates (matrix%m, stored_row, stored_col, tr)
    nze = row_size*col_size
    !
    ! check for 0 block
    IF(nze.LE.0) RETURN

    CALL cp_assert (SIZE(block) .GE. nze, cp_fatal_level,&
         cp_caller_error, routineN, "Invalid block dimensions")
    CALL dbcsr_get_stored_block_info (matrix%m, stored_row, stored_col,&
         found, blk, offset)
    IF(found) THEN
       !
       ! let's copy the block
       offset = ABS (offset)
       ! Fix the index if the new block's transpose flag is different
       ! from the old one.
       IF (matrix%m%blk_p(blk).LT.0 .NEQV. tr) THEN
          matrix%m%blk_p(blk) = -matrix%m%blk_p(blk)
       ENDIF
       block_1d => pointer_view (dbcsr_get_data_p (&
            matrix%m%data_area, 0.0_real_8), offset, offset+nze-1)
       IF (do_sum) THEN
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
    ELSE
       !!@@@
       !call cp_assert (associated (matrix%m%wms), cp_fatal_level,&
       !     cp_caller_error, routineN, "Work matrices not prepared")
       IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
          CALL dbcsr_work_create (matrix, nblks_guess=1,&
               sizedata_guess=SIZE(block))
       ENDIF
       iw = 1
!$     iw = omp_get_thread_num()+1
       blk_p = matrix%m%wms(iw)%datasize + 1
       IF (.NOT.dbcsr_wm_use_mutable (matrix%m%wms(iw))) THEN
          IF (tr) blk_p = -blk_p
          CALL add_work_coordinate (matrix%m%wms(iw), row, col, blk_p)
          CALL dbcsr_data_ensure_size (matrix%m%wms(iw)%data_area,&
               matrix%m%wms(iw)%datasize+SIZE(block))
          CALL addto_array (matrix%m%wms(iw)%data_area, block, ABS(blk_p))
       ELSE
          ALLOCATE (data_block%p (row_size, col_size))
          data_block%p(:,:) = RESHAPE (block, (/row_size, col_size/))
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
                CALL cp_assert (ASSOCIATED (data_block2%p), cp_warning_level,&
                     cp_internal_error, routineN,&
                     "Data was not present in block")
                IF (ASSOCIATED (data_block2%p)) DEALLOCATE (data_block2%p)
             ENDIF
          ELSE
             CALL btree_add_d (&
                  matrix%m%wms(iw)%mutable%m%btree_d,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block, found, data_block2, replace=.FALSE.)
             IF (found) THEN
                CALL daxpy (nze, 1.0_real_8, block(1), 1,&
                     data_block2%p(1,1), 1)
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
       matrix%m%valid = .FALSE.
    ENDIF
  END SUBROUTINE dbcsr_put_block_d


  SUBROUTINE bcast_block_d(blk, source, scope, mp_obj, error)
    REAL(kind=real_8), DIMENSION(:), INTENT(INOUT)     :: blk
    INTEGER, INTENT(IN)                      :: source
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'bcast_block_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol, irow, mp_group, mynode, &
                                                mypcol, myprow, npcols, &
                                                nprows, numnodes, src
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    REAL(kind=real_8), ALLOCATABLE, DIMENSION(:)       :: buff

!   ---------------------------------------------------------------------------

    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)

    ALLOCATE(buff(SIZE(blk)))
    SELECT CASE(scope)
    CASE('rowise')
       !
       ! simple hack
       DO icol = 0,npcols-1
          buff = blk
          src = blacs2mpi(source,icol)
          CALL mp_bcast(buff,src,mp_group)
          IF(mypcol.EQ.icol) THEN
             blk = buff
          ENDIF
       ENDDO
    CASE('columnwise')
       !
       ! simple hack
       DO irow = 0,nprows-1
          buff = blk
          src = blacs2mpi(irow,source)
          CALL mp_bcast(buff,src,mp_group)
          IF(myprow.EQ.irow) THEN
             blk = buff
          ENDIF
       ENDDO
    CASE('all')
       CALL stop_program(routineN,"more work here")
    CASE DEFAULT
       CALL stop_program(routineN,"how do you wanna bcast")
    END SELECT
    DEALLOCATE(buff)

  END SUBROUTINE bcast_block_d

  SUBROUTINE block_add_on_diag_d(m, blk, alpha)
    INTEGER, INTENT(in)                      :: m
    REAL(kind=real_8), INTENT(inout)                   :: blk(m,m)
    REAL(kind=real_8), INTENT(in)                      :: alpha

    CHARACTER(len=*), PARAMETER :: routineN = 'block_add_on_diag_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    FORALL (i = 1:m)
       blk(i,i) = blk(i,i) + alpha
    END FORALL
  END SUBROUTINE block_add_on_diag_d

  SUBROUTINE block_set_d (m, n, blk, alpha, beta)
    INTEGER                                  :: m, n
    REAL(kind=real_8)                                  :: blk(m,n)
    REAL(kind=real_8), OPTIONAL                        :: alpha, beta

    CHARACTER(len=*), PARAMETER :: routineN = 'block_set_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(dp)                                 :: my_alpha, my_beta

!   ---------------------------------------------------------------------------

    my_beta = 0.0_dp
    my_alpha = 0.0_dp
    IF(PRESENT(alpha)) my_alpha = alpha
    IF(PRESENT(beta)) my_beta = beta
    blk(:,:) = my_beta
    IF(m.EQ.n) THEN
       FORALL (i = 1:m)
          blk(i,i) = my_alpha
       END FORALL
    ENDIF
  END SUBROUTINE block_set_d

  SUBROUTINE block_chol_inv_d(m, blk)
    INTEGER                                  :: m
    REAL(kind=real_8)                                  :: blk(m,m)

    CHARACTER(len=*), PARAMETER :: routineN = 'block_chol_inv_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, info, j

!   ---------------------------------------------------------------------------

    CALL dpotrf( 'U', m, blk, m, info )
    IF(info.NE.0)CALL stop_program(routineN,"error in dpotrf")
    CALL dpotri( 'U', m, blk, m, info )
    IF(info.NE.0)CALL stop_program(routineN,"error in dpotri")
    !
    ! symmetrize
    DO i=1,m
       DO j=i,m
          blk(j,i) = blk(i,j)
       ENDDO
    ENDDO
  END SUBROUTINE block_chol_inv_d


#if !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data using ugly hacks.
! *****************************************************************************
  SUBROUTINE pointer_d_rank_remap2 (r2p, d1, d2, r1p)
    REAL(kind=real_8), DIMENSION(:,:), POINTER :: r2p
    INTEGER, INTENT(IN)              :: d1, d2
    REAL(kind=real_8), DIMENSION(d1,*), TARGET :: r1p
    r2p => r1p(1:d1, 1:d2)
  END SUBROUTINE pointer_d_rank_remap2
#endif

