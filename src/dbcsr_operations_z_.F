! *****************************************************************************
!> \brief traces a DBCSR matrix
!> \param[in] matrix_a       DBCSR matrix
!> \param[out] trace         the trace of the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_trace_a_z(matrix_a, trace, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: trace
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_a_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_row, a_row_size, i, &
                                                timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: a_data, data_p

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (matrix_a%m%row_blk_size)
    col_blk_size => array_data (matrix_a%m%col_blk_size)
    CALL dbcsr_get_data (matrix_a%m%data_area, data_p)
    !
    ! let's go
    trace = 0.0_dp
    DO a_row = 1, matrix_a%m%nblkrows_total
       a_row_size = row_blk_size(a_row)
       DO a_blk = matrix_a%m%row_p(a_row)+1,matrix_a%m%row_p(a_row+1)
          a_col = matrix_a%m%col_i(a_blk)
          IF(a_col.ne.a_row) CYCLE
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

  END SUBROUTINE dbcsr_trace_a_z

! *****************************************************************************
!> \brief traces a product of DBCSR matrices
!> \param[in] matrix_a,matrix_b  DBCSR matrices
!> \param[out] trace             the trace of the product of the matrices
!> \param[in] trans_a            (optional) is matrix_a transposed or not?
!> \param[in] trans_b            (optional) is matrix_b transposed or not?
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_trace_ab_z(matrix_a, matrix_b, trace, trans_a, trans_b, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a, matrix_b
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_ab_z', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: my_trans_a, my_trans_b
    INTEGER :: a_blk, a_col, a_col_size, a_row_size, b_blk, b_col_size, &
      b_frst_blk, b_last_blk, b_row_size, nze, row, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, &
                                                a_row_blk_size, &
                                                b_col_blk_size, b_row_blk_size
    LOGICAL                                  :: found
    COMPLEX(kind=real_8), EXTERNAL                        :: ZDOTU
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: a_data_p, b_data_p



!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
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
       IF(a_row_size.NE.b_row_size) CALL stop_program(routineN,"matrices not consistent")
       b_blk = matrix_b%m%row_p(row)+1
       b_frst_blk = matrix_b%m%row_p(row)+1
       b_last_blk = matrix_b%m%row_p(row+1)
       DO a_blk = matrix_a%m%row_p(row)+1,matrix_a%m%row_p(row+1)
          a_col = matrix_a%m%col_i(a_blk)
          a_col_size = a_col_blk_size(a_col)
          !
          ! find the b_blk we assume here that the colums are ordered !
          CALL dbcsr_find_column(a_col,b_frst_blk,b_last_blk,matrix_b%m%col_i,b_blk,found)
          IF(found) THEN
             b_col_size = b_col_blk_size(a_col)
             IF(a_col_size.NE.b_col_size) CALL stop_program(routineN,"matrices not consistent")
             !
             nze = a_row_size*a_col_size
             !
             ! let's trace the blocks
             !trace = trace + DOT_PRODUCT(a_data,b_data)
                trace = trace + ZDOTU (nze,&
                     a_data_p(ABS(matrix_a%m%blk_p(a_blk))),1,&
                     b_data_p(ABS(matrix_b%m%blk_p(b_blk))),1)
             !trace = trace + ddot(nze,matrix_a%m%data(ABS(matrix_a%m%blk_p(a_blk))),1,&
             !     matrix_b%m%data(ABS(matrix_b%m%blk_p(b_blk))),1)
          ENDIF
       ENDDO ! a_col
    ENDDO ! a_row
    !
    ! summe
    CALL mp_sum(trace,dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%m%dist)))

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_trace_ab_z


! *****************************************************************************
!> \brief Gets a 2-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-2 array)
!> \param[out] found  whether the block exists in the matrix
! *****************************************************************************
  SUBROUTINE dbcsr_get_2d_block_z(matrix,row,col,block,found)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:,:), INTENT(OUT)     :: block
    LOGICAL, INTENT(OUT)                     :: found
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_z', &
      routineP = moduleN//':'//routineN
    
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: block_1d
    INTEGER                                  :: row_size, col_size,&
                                                blk, last_blk, nze, offset,&
                                                frst_blk
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size

!   ---------------------------------------------------------------------------
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    CALL cp_assert (matrix%m%data_type .EQ. dbcsr_type_complex_8,&
         cp_fatal_level, cp_caller_error,&
         routineN, "Data type mismatch for requested block.")
    row_size = row_blk_size(row)
    col_size = col_blk_size(col)
    blk = matrix%m%row_p(row)+1
    frst_blk = matrix%m%row_p(row)+1
    last_blk = matrix%m%row_p(row+1)
    CALL dbcsr_find_column(col,frst_blk,last_blk,matrix%m%col_i,blk,found)
    IF(found) THEN
       nze = row_size*col_size
       !
       ! let's copy the block
       offset = ABS(matrix%m%blk_p(blk))
       block_1d = pointer_view ( dbcsr_get_data_p (&
                     matrix%m%data_area, CMPLX(0.0, 0.0, real_8)), offset, offset)
       block(1:row_size, 1:col_size)&
            = RESHAPE (block_1d, (/row_size, col_size/))
    ENDIF
  END SUBROUTINE dbcsr_get_2d_block_z

! *****************************************************************************
!> \brief Gets a 2-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-1 array)
!> \param[out] found  whether the block exists in the matrix
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_z(matrix,row,col,block,found)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(OUT)       :: block
    LOGICAL, INTENT(OUT)                     :: found

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col_size, frst_blk, &
                                                last_blk, nze, offset, &
                                                row_size
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
!   ---------------------------------------------------------------------------
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    CALL cp_assert (matrix%m%data_type .EQ. dbcsr_type_complex_8,&
         cp_fatal_level, cp_caller_error,&
         routineN, "Data type mismatch for requested block.")

    row_size = row_blk_size(row)
    col_size = col_blk_size(col)
    blk = matrix%m%row_p(row)+1
    frst_blk = matrix%m%row_p(row)+1
    last_blk = matrix%m%row_p(row+1)
    CALL dbcsr_find_column(col,frst_blk,last_blk,matrix%m%col_i,blk,found)
    IF(found) THEN
       nze = row_size*col_size
       !
       ! let's copy the block
       offset = ABS(matrix%m%blk_p(blk))
       CALL zcopy(nze,&
            pointer_view (&
               dbcsr_get_data_p (matrix%m%data_area, CMPLX(0.0, 0.0, real_8)), offset, offset&
            ), 1, block(1), 1)
       !CALL dcopy(nze,matrix%m%data(offset),1,block(1),1)
    ENDIF
  END SUBROUTINE dbcsr_get_block_z


! *****************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix using the btree
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to reserve; added if not NULL
!> \param[in] transposed      the block holds transposed data
!> \param[out] existed        (optional) block already existed
! *****************************************************************************
  SUBROUTINE dbcsr_reserve_block2d_z(matrix, row, col, block,&
       transposed, existed)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_block2d_z', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_z)          :: data_block, data_block2
    INTEGER                                  :: blk, col_size, frst_blk, &
                                                last_blk, nze, offset, &
                                                row_size, blk_p,&
                                                stored_row, stored_col
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found, gift, tr, tr2
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER         :: original_block

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
    IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_create (matrix)
       matrix%m%valid = .FALSE.
       CALL dbcsr_mutable_new (matrix%m%wms(1)%mutable,&
            dbcsr_get_data_type (matrix))
    ENDIF

    NULLIFY (data_block%p)
    IF (.NOT. gift) THEN
       ALLOCATE (data_block%p (row_size, col_size))
       block => data_block%p
    ELSE
       data_block%p => block
    ENDIF
    data_block%tr = tr

    CALL btree_add_z (matrix%m%wms(1)%mutable%m%btree_z,&
         make_coordinate_tuple(stored_row, stored_col),&
         data_block, found, data_block2)

    IF (.NOT. found) THEN
       matrix%m%wms(1)%lastblk = matrix%m%wms(1)%lastblk + 1
       matrix%m%wms(1)%datasize = matrix%m%wms(1)%datasize + row_size*col_size
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
  END SUBROUTINE dbcsr_reserve_block2d_z

! *****************************************************************************
!> \brief Get a 2D block from the working matrix
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          pointer to data block
!> \param[out] existed        block was found
! *****************************************************************************
  SUBROUTINE dbcsr_get_tmp_block2d_z(matrix, row, col, block, found)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER         :: block
    LOGICAL, INTENT(OUT)                     :: found

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_tmp_block2d_z', &
      routineP = moduleN//':'//routineN

    TYPE(btree_2d_data_z)          :: data_block

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
       found = .FALSE.
       NULLIFY (block)
    ENDIF

    CALL btree_get_z (matrix%m%wms(1)%mutable%m%btree_z,&
         make_coordinate_tuple(row, col),&
         data_block, found)

    IF (.NOT. found) THEN
       !WRITE(*,*)routineN//" Block did not exist"
       NULLIFY (block)
       !block => data_block%p
    ELSE
       !WRITE(*,*)routineN//" Found block", ASSOCIATED (data_block%p)
       block => data_block%p
    ENDIF
  END SUBROUTINE dbcsr_get_tmp_block2d_z



! *****************************************************************************
!> \brief Put a 2-D block in a DBCSR matrix
!> \param[in.out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to put
!> \param[in]  transposed     the block is transposed
! *****************************************************************************
  SUBROUTINE dbcsr_put_block2d_z(matrix, row, col, block, transposed)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:,:), INTENT(IN)      :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block2d_z', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: tr
    
    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    CALL dbcsr_put_block (matrix, row, col, RESHAPE (block, (/SIZE(block)/)),&
      tr)
  END SUBROUTINE dbcsr_put_block2d_z

! *****************************************************************************
!> \brief put a block in a dbcsr matrix
!> \param[in]  matrix         DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to put
!> \param[in]  transposed     the block is transposed
! *****************************************************************************
  SUBROUTINE dbcsr_put_block_z(matrix, row, col, block, transposed)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN)        :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col_size, frst_blk, &
                                                last_blk, nze, offset, &
                                                row_size, blk_p
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found, tr

! IN PARALLEL WE SHOULD ENSURE THAT THE BLOCK GOES TO THE RIGHT PLACE !

    IF (PRESENT (transposed)) THEN
       tr = transposed
    ELSE
       tr = .FALSE.
    ENDIF
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)

    row_size = row_blk_size(row)
    col_size = col_blk_size(col)
    blk = matrix%m%row_p(row)+1
    frst_blk = matrix%m%row_p(row)+1
    last_blk = matrix%m%row_p(row+1)
    CALL dbcsr_find_column(col,frst_blk,last_blk,matrix%m%col_i,blk,found)
    IF(found) THEN
       nze = row_size*col_size
       !
       ! let's copy the block
       offset = ABS(matrix%m%blk_p(blk))
       IF (matrix%m%blk_p(blk).LT.0 .NEQV. tr) THEN
          matrix%m%blk_p(blk) = -matrix%m%blk_p(blk)
       ENDIF
       CALL zcopy (nze, block(1), 1,&
            pointer_view (&
            dbcsr_get_data_p (matrix%m%data_area, CMPLX(0.0, 0.0, real_8)), offset, offset&
            ), 1)
       !CALL dcopy(nze,block(1),1,matrix%m%data(offset),1)
    ELSE
       !WRITE(*,*)routineN//" Adding block to matrix"
       IF (.NOT.ASSOCIATED (matrix%m%wms)) THEN
          !WRITE(*,*)routineN//" Creating new work matrix to hold new data"
          CALL dbcsr_work_create (matrix, nblks_guess=1,&
               sizedata_guess=SIZE(block))
       ENDIF
       blk_p = matrix%m%wms(1)%datasize + 1
       IF (tr) blk_p = -blk_p
       CALL add_work_coordinate (matrix%m%wms(1), row, col, blk_p)
       CALL addto_array (matrix%m%wms(1)%data_area, block, blk_p)
       matrix%m%wms(1)%datasize = matrix%m%wms(1)%datasize + SIZE (block)
       !matrix%m%wms(1)%datasize = dbcsr_get_data_size (matrix%m%wms(1)%data_area)
       !WRITE(*,*)routineN//" Volume of new data is",matrix%m%wms(1)%datasize,&
       !     dbcsr_get_data_size (matrix%m%wms(1)%data_area)
       matrix%m%valid = .FALSE.
    ENDIF
  END SUBROUTINE dbcsr_put_block_z


  SUBROUTINE bcast_block_z(blk, source, scope, mp_obj, error)
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(INOUT)     :: blk
    INTEGER, INTENT(IN)                      :: source
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'bcast_block_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol, irow, mp_group, mynode, &
                                                mypcol, myprow, npcols, &
                                                nprows, numnodes, src
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    COMPLEX(kind=real_8), ALLOCATABLE, DIMENSION(:)       :: buff

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

  END SUBROUTINE bcast_block_z

  SUBROUTINE block_add_on_diag_z(m, blk, alpha)
    INTEGER, INTENT(in)                      :: m
    COMPLEX(kind=real_8), INTENT(inout)                   :: blk(m,m)
    COMPLEX(kind=real_8), INTENT(in)                      :: alpha

    CHARACTER(len=*), PARAMETER :: routineN = 'block_add_on_diag_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    FORALL (i = 1:m)
       blk(i,i) = blk(i,i) + alpha
    END FORALL
  END SUBROUTINE block_add_on_diag_z

  SUBROUTINE block_set_z (m, n, blk, alpha, beta)
    INTEGER                                  :: m, n
    COMPLEX(kind=real_8)                                  :: blk(m,n)
    COMPLEX(kind=real_8), OPTIONAL                        :: alpha, beta

    CHARACTER(len=*), PARAMETER :: routineN = 'block_set_z', &
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
  END SUBROUTINE block_set_z

  SUBROUTINE block_chol_inv_z(m, blk)
    INTEGER                                  :: m
    COMPLEX(kind=real_8)                                  :: blk(m,m)

    CHARACTER(len=*), PARAMETER :: routineN = 'block_chol_inv_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, info, j

!   ---------------------------------------------------------------------------

    CALL zpotrf( 'U', m, blk, m, info )
    IF(info.NE.0)CALL stop_program(routineN,"error in dpotrf")
    CALL zpotri( 'U', m, blk, m, info )
    IF(info.NE.0)CALL stop_program(routineN,"error in dpotri")
    !
    ! symmetrize
    DO i=1,m
       DO j=i,m
          blk(j,i) = blk(i,j)
       ENDDO
    ENDDO
  END SUBROUTINE block_chol_inv_z
