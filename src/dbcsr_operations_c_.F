! *****************************************************************************
!> \brief traces a DBCSR matrix
!> \param[in] matrix_a       DBCSR matrix
!> \param[out] trace         the trace of the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_trace_a_c(matrix_a, trace, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: trace
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_a_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_row, a_row_size, i, &
                                                timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: a_data, data_p

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (matrix_a%row_blk_size)
    col_blk_size => array_data (matrix_a%col_blk_size)
    CALL dbcsr_get_data (matrix_a%data_area, data_p)
    !
    ! let's go
    trace = 0.0_sp
    DO a_row = 1, matrix_a%nblkrows_total
       a_row_size = row_blk_size(a_row)
       DO a_blk = matrix_a%row_p(a_row)+1,matrix_a%row_p(a_row+1)
          a_col = matrix_a%col_i(a_blk)
          IF(a_col.ne.a_row) CYCLE
          a_col_size = col_blk_size(a_col)
          IF(a_row_size.NE.a_col_size) CALL stop_program(routineN,"is that a square matrix?")
          a_nze = a_row_size**2
          a_data => pointer_view (data_p, ABS(matrix_a%blk_p(a_blk)),&
               ABS(matrix_a%blk_p(a_blk))+a_nze-1)
          !data_a => matrix_a%data(ABS(matrix_a%blk_p(a_blk)):ABS(matrix_a%blk_p(a_blk))+a_nze-1)
          !
          ! let's trace the block
          DO i = 1,a_row_size
             trace = trace + a_data((i-1)*a_row_size+i)
          ENDDO
       ENDDO ! a_col
    ENDDO ! a_row
    !
    ! summe
    CALL mp_sum(trace,dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_trace_a_c

! *****************************************************************************
!> \brief traces a product of DBCSR matrices
!> \param[in] matrix_a,matrix_b  DBCSR matrices
!> \param[out] trace             the trace of the product of the matrices
!> \param[in] trans_a            (optional) is matrix_a transposed or not?
!> \param[in] trans_b            (optional) is matrix_b transposed or not?
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_trace_ab_c(matrix_a, matrix_b, trace, trans_a, trans_b, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a, matrix_b
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_ab_c', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: my_trans_a, my_trans_b
    INTEGER :: a_blk, a_col, a_col_size, a_row_size, b_blk, b_col_size, &
      b_frst_blk, b_last_blk, b_row_size, nze, row, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, &
                                                a_row_blk_size, &
                                                b_col_blk_size, b_row_blk_size
    LOGICAL                                  :: found
    COMPLEX(kind=real_4), EXTERNAL                        :: CDOTU
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: a_data_p, b_data_p



!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    a_row_blk_size => array_data (matrix_a%row_blk_size)
    a_col_blk_size => array_data (matrix_a%col_blk_size)
    b_row_blk_size => array_data (matrix_b%row_blk_size)
    b_col_blk_size => array_data (matrix_b%col_blk_size)
    CALL dbcsr_get_data (matrix_a%data_area, a_data_p)
    CALL dbcsr_get_data (matrix_b%data_area, b_data_p)
    
    my_trans_a = 'T'
    IF(PRESENT(trans_a)) my_trans_a = trans_a
    my_trans_b = 'N'
    IF(PRESENT(trans_b)) my_trans_b = trans_b
    IF(my_trans_a.NE.'T'.OR.my_trans_b.NE.'N') CALL stop_program(routineN,"NYI")
    !
    ! let's go
    trace = 0.0_sp
    IF(matrix_a%nblkrows_total.NE.matrix_b%nblkrows_total) THEN
       CALL stop_program(routineN,"matrices not consistent")
    ENDIF
    DO row = 1, matrix_a%nblkrows_total
       a_row_size = a_row_blk_size(row)
       b_row_size = b_row_blk_size(row)
       IF(a_row_size.NE.b_row_size) CALL stop_program(routineN,"matrices not consistent")
       b_blk = matrix_b%row_p(row)+1
       b_frst_blk = matrix_b%row_p(row)+1
       b_last_blk = matrix_b%row_p(row+1)
       DO a_blk = matrix_a%row_p(row)+1,matrix_a%row_p(row+1)
          a_col = matrix_a%col_i(a_blk)
          a_col_size = a_col_blk_size(a_col)
          !
          ! find the b_blk we assume here that the colums are ordered !
          CALL dbcsr_find_column(a_col,b_frst_blk,b_last_blk,matrix_b%col_i,b_blk,found)
          IF(found) THEN
             b_col_size = b_col_blk_size(a_col)
             IF(a_col_size.NE.b_col_size) CALL stop_program(routineN,"matrices not consistent")
             !
             nze = a_row_size*a_col_size
             !
             ! let's trace the blocks
             !trace = trace + DOT_PRODUCT(a_data,b_data)
                trace = trace + CDOTU (nze,&
                     a_data_p(ABS(matrix_a%blk_p(a_blk))),1,&
                     b_data_p(ABS(matrix_b%blk_p(b_blk))),1)
             !trace = trace + ddot(nze,matrix_a%data(ABS(matrix_a%blk_p(a_blk))),1,&
             !     matrix_b%data(ABS(matrix_b%blk_p(b_blk))),1)
          ENDIF
       ENDDO ! a_col
    ENDDO ! a_row
    !
    ! summe
    CALL mp_sum(trace,dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_trace_ab_c


! *****************************************************************************
!> \brief get a block in a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  block  the block to get
!> \param[in]  row    the row
!> \param[in]  col    the column
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_c(matrix,block,row,col,found,error)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    COMPLEX(kind=real_4), DIMENSION(:), INTENT(INOUT) :: block
    INTEGER, INTENT(IN)                      :: row, col
    LOGICAL, INTENT(OUT)                     :: found
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col_size, frst_blk, &
                                                last_blk, nze, offset, &
                                                row_size
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size

    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)
    CALL cp_assert (matrix%data_type .EQ. 'c',&
         cp_fatal_level, cp_caller_error,&
         routineN, "Data type mismatch for requested block.")

    row_size = row_blk_size(row)
    col_size = col_blk_size(col)
    blk = matrix%row_p(row)+1
    frst_blk = matrix%row_p(row)+1
    last_blk = matrix%row_p(row+1)
    CALL dbcsr_find_column(col,frst_blk,last_blk,matrix%col_i,blk,found)
    IF(found) THEN
       nze = row_size*col_size
       !
       ! let's copy the block
       offset = ABS(matrix%blk_p(blk))
       CALL ccopy(nze,&
            pointer_view (&
               dbcsr_get_data_p (matrix%data_area, CMPLX(0.0, 0.0, real_4)), offset, offset&
            ), 1, block(1), 1)
       !CALL dcopy(nze,matrix%data(offset),1,block(1),1)
    ENDIF

  END SUBROUTINE dbcsr_get_block_c

! *****************************************************************************
!> \brief put a block in a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  block  the block to put
!> \param[in]  row    the row
!> \param[in]  col    the column
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_put_block_c(matrix, block, row, col, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    COMPLEX(kind=real_4), DIMENSION(:), INTENT(IN)        :: block
    INTEGER, INTENT(IN)                      :: row, col
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col_size, frst_blk, &
                                                last_blk, nze, offset, &
                                                row_size
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found

! IN PARALLEL WE SHOULD ENSURE THAT THE BLOCK GOES TO THE RIGHT PLACE !

    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)

    row_size = row_blk_size(row)
    col_size = col_blk_size(col)
    blk = matrix%row_p(row)+1
    frst_blk = matrix%row_p(row)+1
    last_blk = matrix%row_p(row+1)
    CALL dbcsr_find_column(col,frst_blk,last_blk,matrix%col_i,blk,found)
    IF(found) THEN
       nze = row_size*col_size
       !
       ! let's copy the block
       offset = ABS(matrix%blk_p(blk))
       CALL dcopy (nze, block(1), 1,&
            pointer_view (&
               dbcsr_get_data_p (matrix%data_area, CMPLX(0.0, 0.0, real_4)), offset, offset&
            ), 1)
       !CALL dcopy(nze,block(1),1,matrix%data(offset),1)
    ELSE
       CALL stop_program(routineN,"need ot create the block, NYI")
    ENDIF

  END SUBROUTINE dbcsr_put_block_c


  SUBROUTINE bcast_block_c(blk, source, scope, mp_obj, error)
    COMPLEX(kind=real_4), DIMENSION(:), INTENT(INOUT)     :: blk
    INTEGER, INTENT(IN)                      :: source
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'bcast_block_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol, irow, mp_group, mynode, &
                                                mypcol, myprow, npcols, &
                                                nprows, numnodes, src
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    COMPLEX(kind=real_4), ALLOCATABLE, DIMENSION(:)       :: buff

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

  END SUBROUTINE bcast_block_c

  SUBROUTINE block_add_on_diag_c(m, blk, alpha)
    INTEGER, INTENT(in)                      :: m
    COMPLEX(kind=real_4), INTENT(inout)                   :: blk(m,m)
    COMPLEX(kind=real_4), INTENT(in)                      :: alpha

    CHARACTER(len=*), PARAMETER :: routineN = 'block_add_on_diag_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    FORALL (i = 1:m)
       blk(i,i) = blk(i,i) + alpha
    END FORALL
  END SUBROUTINE block_add_on_diag_c

  SUBROUTINE block_set_c (m, n, blk, alpha, beta)
    INTEGER                                  :: m, n
    COMPLEX(kind=real_4)                                  :: blk(m,n)
    COMPLEX(kind=real_4), OPTIONAL                        :: alpha, beta

    CHARACTER(len=*), PARAMETER :: routineN = 'block_set_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(dp)                                 :: my_alpha, my_beta

!   ---------------------------------------------------------------------------

    my_beta = 0.0_sp
    my_alpha = 0.0_sp
    IF(PRESENT(alpha)) my_alpha = alpha
    IF(PRESENT(beta)) my_beta = beta
    blk(:,:) = my_beta
    IF(m.EQ.n) THEN
       FORALL (i = 1:m)
          blk(i,i) = my_alpha
       END FORALL
    ENDIF
  END SUBROUTINE block_set_c

  SUBROUTINE block_chol_inv_c(m, blk)
    INTEGER                                  :: m
    COMPLEX(kind=real_4)                                  :: blk(m,m)

    CHARACTER(len=*), PARAMETER :: routineN = 'block_chol_inv_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, info, j

!   ---------------------------------------------------------------------------

    CALL cpotrf( 'U', m, blk, m, info )
    IF(info.NE.0)CALL stop_program(routineN,"error in dpotrf")
    CALL cpotri( 'U', m, blk, m, info )
    IF(info.NE.0)CALL stop_program(routineN,"error in dpotri")
    !
    ! symmetrize
    FORALL (i = 1:m, j = 1:m)
       !DO i=1,m
       !DO j=i,m
       blk(j,i) = blk(i,j)
       !ENDDO
       !ENDDO
    END FORALL
  END SUBROUTINE block_chol_inv_c

