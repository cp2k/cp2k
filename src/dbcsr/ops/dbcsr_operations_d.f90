!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief traces a DBCSR matrix
!> \param[in] matrix_a       DBCSR matrix
!> \param[out] trace         the trace of the matrix
!>
! **************************************************************************************************
  SUBROUTINE dbcsr_trace_a_d(matrix_a, trace)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix_a
    REAL(kind=real_8), INTENT(INOUT)                   :: trace

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_a_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_row, a_row_size, i, &
                                                mynode, error_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size,&
                                                row_dist, col_dist
    REAL(kind=real_8), DIMENSION(:), POINTER           :: a_data, data_p
    INTEGER, DIMENSION(:,:), POINTER         :: pgrid
    TYPE(dbcsr_distribution_obj)             :: dist

!   ---------------------------------------------------------------------------
    CALL timeset(routineN, error_handle)

    row_blk_size => array_data (matrix_a%row_blk_size)
    col_blk_size => array_data (matrix_a%col_blk_size)
    IF(dbcsr_get_data_type (matrix_a) /=  dbcsr_type_real_8) &
       CPABORT("Incompatible data types")
    CALL dbcsr_get_data (matrix_a%data_area, data_p)
    dist = dbcsr_distribution (matrix_a)
    mynode = dbcsr_mp_mynode (dbcsr_distribution_mp (dist))
    pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (dist))
    row_dist => dbcsr_distribution_row_dist (dist)
    col_dist => dbcsr_distribution_col_dist (dist)
    !
    ! let's go
    trace = REAL(0.0,real_8)
    DO a_row = 1, matrix_a%nblkrows_total
       a_row_size = row_blk_size(a_row)
       DO a_blk = matrix_a%row_p(a_row)+1,matrix_a%row_p(a_row+1)
          IF (a_blk .EQ. 0) CYCLE
          a_col = matrix_a%col_i(a_blk)
          IF(a_col.ne.a_row) CYCLE
          ! We must skip non-local blocks in a replicated matrix.
          IF(matrix_a%replication_type .NE. dbcsr_repl_full) THEN
             IF (mynode .NE. checker_square_proc (a_row, a_col, pgrid, row_dist, col_dist)) &
                CYCLE
          ENDIF
          a_col_size = col_blk_size(a_col)
          IF(a_row_size.NE.a_col_size)&
             CPABORT("is that a square matrix?")
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

    CALL timestop(error_handle)
  END SUBROUTINE dbcsr_trace_a_d

! **************************************************************************************************
!> \brief traces a product of DBCSR matrices
!> \param[in] matrix_a DBCSR matrices
!> \param[in] matrix_b DBCSR matrices
!> \param[out] trace             the trace of the product of the matrices
! **************************************************************************************************
  SUBROUTINE dbcsr_trace_ab_d(matrix_a, matrix_b, trace)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix_a, matrix_b
    REAL(kind=real_8), INTENT(INOUT)                   :: trace

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_ab_d', &
      routineP = moduleN//':'//routineN

    INTEGER :: a_blk, a_col, a_col_size, a_row_size, b_blk, b_col_size, &
      b_frst_blk, b_last_blk, b_row_size, nze, row, a_beg, a_end, b_beg, b_end
    CHARACTER                                :: matrix_a_type, matrix_b_type
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, &
                                                a_row_blk_size, &
                                                b_col_blk_size, b_row_blk_size
    REAL(kind=real_8)                                  :: sym_fac, fac
    LOGICAL                                  :: found, matrix_a_symm, matrix_b_symm
    REAL(real_4), EXTERNAL                   :: SDOT
    REAL(real_8), EXTERNAL                   :: DDOT
    COMPLEX(real_4), EXTERNAL                :: CDOTU
    COMPLEX(real_8), EXTERNAL                :: ZDOTU
    REAL(kind=real_8), DIMENSION(:), POINTER           :: a_data, b_data

!   ---------------------------------------------------------------------------

    IF(matrix_a%replication_type .NE. dbcsr_repl_none&
         .OR. matrix_b%replication_type .NE. dbcsr_repl_none)&
         CPABORT("Trace of product of replicated matrices not yet possible.")

    sym_fac = REAL(1.0,real_8)
    matrix_a_type = dbcsr_get_matrix_type(matrix_a)
    matrix_b_type = dbcsr_get_matrix_type(matrix_b)
    matrix_a_symm = matrix_a_type==dbcsr_type_symmetric.OR.matrix_a_type==dbcsr_type_antisymmetric
    matrix_b_symm = matrix_b_type==dbcsr_type_symmetric.OR.matrix_b_type==dbcsr_type_antisymmetric

    IF(matrix_a_symm .AND. matrix_b_symm) sym_fac = REAL(2.0,real_8)

    ! tracing a symmetric with a general matrix is not implemented, as it would require communication of blocks
    IF(matrix_a_symm .NEQV. matrix_b_symm) &
       CPABORT("Tracing general with symmetric matrix NYI")

    a_row_blk_size => array_data (matrix_a%row_blk_size)
    a_col_blk_size => array_data (matrix_a%col_blk_size)
    b_row_blk_size => array_data (matrix_b%row_blk_size)
    b_col_blk_size => array_data (matrix_b%col_blk_size)

    CALL dbcsr_get_data (matrix_a%data_area, a_data)
    CALL dbcsr_get_data (matrix_b%data_area, b_data)

    ! let's go
    trace = REAL(0.0,real_8)
    IF(matrix_a%nblkrows_total.NE.matrix_b%nblkrows_total)&
       CPABORT("this combination of transpose is NYI")
    DO row = 1, matrix_a%nblkrows_total
       a_row_size = a_row_blk_size(row)
       b_row_size = b_row_blk_size(row)
       IF(a_row_size.NE.b_row_size) CPABORT("matrices not consistent")
       b_blk = matrix_b%row_p(row)+1
       b_frst_blk = matrix_b%row_p(row)+1
       b_last_blk = matrix_b%row_p(row+1)
       DO a_blk = matrix_a%row_p(row)+1,matrix_a%row_p(row+1)
          IF (matrix_a%blk_p(a_blk) .EQ. 0) CYCLE ! Deleted block
          a_col = matrix_a%col_i(a_blk)
          a_col_size = a_col_blk_size(a_col)
          !
          ! find the b_blk we assume here that the colums are ordered !
          CALL dbcsr_find_column(a_col,b_frst_blk,b_last_blk,matrix_b%col_i,&
               matrix_b%blk_p,b_blk,found)
          IF(found) THEN
             b_col_size = b_col_blk_size(a_col)
             IF(a_col_size.NE.b_col_size) CPABORT("matrices not consistent")
             !
             nze = a_row_size*a_col_size
             !
             IF(nze.GT.0) THEN
                !
                ! let's trace the blocks
                a_beg = ABS(matrix_a%blk_p(a_blk))
                a_end = a_beg + nze - 1
                b_beg = ABS(matrix_b%blk_p(b_blk))
                b_end = b_beg + nze - 1
                fac = REAL(1.0,real_8)
                IF(row.NE.a_col) fac = sym_fac

                trace = trace + fac * SUM ( a_data(a_beg:a_end) * b_data(b_beg:b_end) )

             ENDIF
          ENDIF
       ENDDO ! a_col
    ENDDO ! a_row
    !
    ! sum
    CALL mp_sum(trace,dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))

  END SUBROUTINE dbcsr_trace_ab_d


! **************************************************************************************************
!> \brief Interface for matrix scaling by a scalar
!> \param matrix_a ...
!> \param alpha_scalar ...
!> \param last_column ...
! **************************************************************************************************
  SUBROUTINE dbcsr_scale_d(matrix_a, alpha_scalar, last_column)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix_a
    REAL(kind=real_8), INTENT(IN)                      :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    TYPE(dbcsr_scalar_type)                  :: sc

    sc = dbcsr_scalar (alpha_scalar)
    CALL dbcsr_scalar_fill_all (sc)
    sc%data_type = dbcsr_get_data_type (matrix_a)
    CALL timeset(routineN, error_handler)
    IF (PRESENT (last_column)) THEN
       CALL dbcsr_scale_anytype(matrix_a,&
            alpha_scalar=sc,&
            limits=(/0,0,0,last_column/))
    ELSE
       CALL dbcsr_scale_anytype(matrix_a, alpha_scalar=sc)
    ENDIF
    CALL timestop(error_handler)
  END SUBROUTINE dbcsr_scale_d

! **************************************************************************************************
!> \brief Interface for matrix scaling by a vector
!> \param matrix_a ...
!> \param alpha ...
!> \param side ...
! **************************************************************************************************
  SUBROUTINE dbcsr_scale_by_vector_d(matrix_a, alpha, side)
    TYPE(dbcsr_type), INTENT(INOUT)            :: matrix_a
    REAL(kind=real_8), DIMENSION(:), INTENT(IN), TARGET :: alpha
    CHARACTER(LEN=*), INTENT(IN)              :: side
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_by_vector_d', &
      routineP = moduleN//':'//routineN
    REAL(kind=real_8), DIMENSION(:), POINTER            :: tmp_p
    TYPE(dbcsr_data_obj)                      :: enc_alpha_vec

    CALL dbcsr_data_init (enc_alpha_vec)
    CALL dbcsr_data_new (enc_alpha_vec, dbcsr_type_real_8)
    tmp_p => alpha
    CALL dbcsr_data_set_pointer (enc_alpha_vec, tmp_p)
    CALL dbcsr_scale_by_vector_anytype(matrix_a, enc_alpha_vec, side)
    CALL dbcsr_data_clear_pointer (enc_alpha_vec)
    CALL dbcsr_data_release (enc_alpha_vec)
  END SUBROUTINE dbcsr_scale_by_vector_d


! **************************************************************************************************
!> \brief Interface for dbcsr_set
!> \param matrix ...
!> \param alpha ...
! **************************************************************************************************
  SUBROUTINE dbcsr_set_d(matrix, alpha)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix
    REAL(kind=real_8), INTENT(IN)                      :: alpha

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set'

      INTEGER                                            :: col, handle, row
      TYPE(dbcsr_iterator)                               :: iter
      REAL(kind=real_8), DIMENSION(:,:), POINTER                   :: block
      LOGICAL                                            :: tr

    CALL timeset(routineN, handle)

    IF (alpha==0.0_real_8) THEN
      CALL dbcsr_zero(matrix)
    ELSE
       IF(dbcsr_get_data_type (matrix) /=  dbcsr_type_real_8) &
         CPABORT("Incompatible data types")

      !TODO: could be speedup by direct assigment to data_area, similar to dbcsr_zero()
      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, tr)
         block(:,:) = alpha
      ENDDO
      CALL dbcsr_iterator_stop(iter)
    ENDIF

    CALL timestop(handle)
  END SUBROUTINE dbcsr_set_d

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eps ...
!> \param method ...
!> \param use_absolute ...
!> \param filter_diag ...
!> \param quick ...
! **************************************************************************************************
  SUBROUTINE dbcsr_filter_d (matrix, eps, method, use_absolute, &
       filter_diag, quick)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix
    REAL(kind=real_8), INTENT(IN)                      :: eps
    INTEGER, INTENT(IN), OPTIONAL            :: method
    LOGICAL, INTENT(in), OPTIONAL            :: use_absolute, filter_diag, quick
    CALL dbcsr_filter_anytype (matrix, dbcsr_scalar(eps), method, &
         use_absolute, filter_diag, quick)
  END SUBROUTINE dbcsr_filter_d

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
  SUBROUTINE dbcsr_set_diag_d(matrix, diag)
    TYPE(dbcsr_type), INTENT(INOUT)            :: matrix
    REAL(kind=real_8), DIMENSION(:), INTENT(IN)          :: diag

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_diag'

    INTEGER                                            :: icol, irow, row_offset, handle, i
    LOGICAL                                            :: tr
    TYPE(dbcsr_iterator)                               :: iter
    REAL(kind=real_8), DIMENSION(:,:), POINTER                   :: block


    CALL timeset(routineN, handle)

    IF(dbcsr_get_data_type (matrix) /=  dbcsr_type_real_8) &
         CPABORT("Incompatible data types")

    IF (dbcsr_nfullrows_total(matrix) /= SIZE(diag)) &
         CPABORT("Diagonal has wrong size")

    IF (.NOT. array_equality(dbcsr_row_block_offsets(matrix), dbcsr_row_block_offsets(matrix))) &
        CPABORT("matrix not quadratic")

    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, irow, icol, block, tr, row_offset=row_offset)
       IF (irow /= icol) CYCLE

       IF(sIZE(block, 1) /= sIZE(block, 2)) &
          CPABORT("Diagonal block non-squared")

       DO i = 1 , sIZE(block, 1)
          block(i,i) = diag(row_offset+i-1)
       END DO
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    CALL timestop(handle)
  END SUBROUTINE dbcsr_set_diag_d

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
  SUBROUTINE dbcsr_get_diag_d(matrix, diag)
    TYPE(dbcsr_type), INTENT(IN)               :: matrix
    REAL(kind=real_8), DIMENSION(:), INTENT(OUT)         :: diag

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_diag'

    INTEGER                                            :: icol, irow, row_offset, handle, i
    LOGICAL                                            :: tr
    TYPE(dbcsr_iterator)                               :: iter
    REAL(kind=real_8), DIMENSION(:,:), POINTER                   :: block


    CALL timeset(routineN, handle)

    IF(dbcsr_get_data_type (matrix) /=  dbcsr_type_real_8) &
         CPABORT("Incompatible data types")

    IF (dbcsr_nfullrows_total(matrix) /= SIZE(diag)) &
         CPABORT("Diagonal has wrong size")

    IF (.NOT. array_equality(dbcsr_row_block_offsets(matrix), dbcsr_row_block_offsets(matrix))) &
        CPABORT("matrix not quadratic")

    diag(:) = 0.0_real_8

    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, irow, icol, block, tr, row_offset=row_offset)
       IF (irow /= icol) CYCLE

       IF(sIZE(block, 1) /= sIZE(block, 2)) &
          CPABORT("Diagonal block non-squared")

       DO i = 1 , sIZE(block, 1)
          diag(row_offset+i-1) = block(i,i)
       END DO
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    CALL timestop(handle)
  END SUBROUTINE dbcsr_get_diag_d

! **************************************************************************************************
!> \brief add a constant to the diagonal of a matrix
!> \param[inout] matrix       DBCSR matrix
!> \param[in]    alpha scalar
! **************************************************************************************************
   SUBROUTINE dbcsr_add_on_diag_d(matrix, alpha)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(kind=real_8), INTENT(IN)                                :: alpha


      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_on_diag'

      INTEGER                                            :: handle, mynode, node, irow, i, row_size
      LOGICAL                                            :: found, tr
      REAL(kind=real_8), DIMENSION(:,:), POINTER                   :: block

      CALL timeset(routineN, handle)

      IF(dbcsr_get_data_type (matrix) /=  dbcsr_type_real_8) &
         CPABORT("Incompatible data types")

      IF (.NOT. array_equality(dbcsr_row_block_offsets(matrix), dbcsr_row_block_offsets(matrix))) &
         CPABORT("matrix not quadratic")

      mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(dbcsr_distribution(matrix)))

      CALL dbcsr_work_create(matrix, work_mutable=.TRUE.)

      DO irow = 1, dbcsr_nblkrows_total(matrix)
         CALL dbcsr_get_stored_coordinates(matrix, irow, irow, node)
         IF (node /= mynode) CYCLE

         CALL dbcsr_get_block_p(matrix, irow, irow, block, tr, found, row_size=row_size)
         IF (.NOT.found) THEN
            ALLOCATE(block(row_size,row_size))
            block(:,:) = 0.0_real_8
         ENDIF

         DO i = 1, row_size
             block(i,i) = block(i,i) + alpha
         END DO

         IF (.NOT.found) THEN
            CALL dbcsr_put_block(matrix, irow, irow, block)
            DEALLOCATE(block)
         ENDIF
      ENDDO

      CALL dbcsr_finalize(matrix)
      CALL timestop(handle)
   END SUBROUTINE dbcsr_add_on_diag_d

