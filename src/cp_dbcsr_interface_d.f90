!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Encapsulates a given scalar value and makes it conformant to the
!>        type of the matrix.
! *****************************************************************************
  FUNCTION make_conformant_scalar_d (scalar, matrix, error) RESULT (encapsulated)
    REAL(kind=real_8), INTENT(IN)                      :: scalar
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_conformant_scalar_d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_scalar_type)                  :: encapsulated
    INTEGER                                  :: data_type, scalar_data_type

    encapsulated = dbcsr_scalar (scalar)
    CALL dbcsr_scalar_fill_all (encapsulated)
    data_type = dbcsr_get_data_type (matrix%matrix)
    scalar_data_type = dbcsr_scalar_get_type(encapsulated)
    IF (scalar_data_type .EQ. dbcsr_type_complex_4 .OR.&
        scalar_data_type .EQ. dbcsr_type_complex_8) THEN
       CALL cp_assert (data_type .EQ. dbcsr_type_complex_4 .OR.&
            data_type .EQ. dbcsr_type_complex_8,&
            cp_fatal_level, cp_wrong_args_error, routineN,&
            "Can not conform a complex to a real number", error=error)
    END IF
    CALL dbcsr_scalar_set_type (encapsulated,data_type)
  END FUNCTION make_conformant_scalar_d


  SUBROUTINE cp_dbcsr_reserve_block2d_d (matrix, row, col, block,&
       transposed, existed) 
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:, :), POINTER        :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_reserve_block2d_d', &
      routineP = moduleN//':'//routineN

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL dbcsr_reserve_block2d(matrix%matrix, row, col, block,&
         transposed, existed)

  END SUBROUTINE cp_dbcsr_reserve_block2d_d


  SUBROUTINE cp_iterator_next_2d_block_d (iterator, row, column,&
       block,&
       block_number, row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    REAL(kind=real_8), DIMENSION(:, :), POINTER        :: block
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number, row_size, &
                                                col_size, row_offset, &
                                                col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterator_next_2d_block_d', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: transposed

    CALL dbcsr_iterator_next_block (iterator, row, column,&
       block, transposed,&
       block_number, row_size, col_size, row_offset, col_offset)
    CALL cp_assert(.NOT. transposed, cp_fatal_level, cp_internal_error,&
         routineN, "CP2K does not handle transposed blocks.")

  END SUBROUTINE cp_iterator_next_2d_block_d


  SUBROUTINE cp_iterator_next_1d_block_d (iterator, row, column, block,&
       block_number, row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)    :: iterator
    INTEGER, INTENT(OUT)                      :: row, column
    REAL(kind=real_8), DIMENSION(:), POINTER            :: block
    INTEGER, INTENT(OUT), OPTIONAL            :: block_number, row_size, &
                                                 col_size, row_offset, &
                                                 col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterator_next_1d_block_d', &
      routineP = moduleN//':'//routineN

    LOGICAL                                   :: transposed

    CALL dbcsr_iterator_next_block (iterator, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)
    CALL cp_assert(.NOT. transposed, cp_fatal_level, cp_internal_error,&
         routineN, "CP2K does not handle transposed blocks.")

  END SUBROUTINE cp_iterator_next_1d_block_d


  SUBROUTINE cp_dbcsr_put_block2d_d (matrix, row, col, block,&
       summation, scale)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:, :), INTENT(IN)     :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: summation
    REAL(kind=real_8), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_put_block2d_d', &
      routineP = moduleN//':'//routineN

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL dbcsr_put_block(matrix%matrix, row, col, block,&
       summation=summation, scale=scale)

  END SUBROUTINE cp_dbcsr_put_block2d_d


  SUBROUTINE cp_dbcsr_put_block_d (matrix, row, col, block,&
       summation, scale)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:), INTENT(IN)        :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: summation
    REAL(kind=real_8), INTENT(IN), OPTIONAL            :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_put_block_d', &
      routineP = moduleN//':'//routineN

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL dbcsr_put_block(matrix%matrix, row, col, block,&
       summation=summation, scale=scale)

  END SUBROUTINE cp_dbcsr_put_block_d


  SUBROUTINE cp_dbcsr_get_block_d (matrix,row,col,block,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:), INTENT(OUT)      :: block
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_d', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: tr
    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    tr=.FALSE.
    CALL dbcsr_get_block(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_block_d


  SUBROUTINE cp_dbcsr_get_2d_block_d (matrix,row,col,block,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:, :), INTENT(OUT)    :: block
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_2d_block_d', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: tr
    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    tr=.FALSE.
    CALL dbcsr_get_block(matrix%matrix,row,col,block,tr,found,&
         row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_2d_block_d


  SUBROUTINE cp_dbcsr_get_2d_block_p_d (matrix,row,col,block,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:, :), POINTER        :: block
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_2d_block_p_d', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: tr

    CALL dbcsr_get_block_p(matrix%matrix,row,col,block,tr,found,&
         row_size, col_size)
    CALL cp_assert(.NOT. tr, cp_fatal_level, cp_internal_error,&
         routineN, "CP2K does not handle transposed blocks.")
  END SUBROUTINE cp_dbcsr_get_2d_block_p_d


  SUBROUTINE cp_dbcsr_get_block_p_d (matrix,row,col,block,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(IN)           :: matrix
    INTEGER, INTENT(IN)                       :: row, col
    REAL(kind=real_8), DIMENSION(:), POINTER            :: block
    LOGICAL, INTENT(OUT)                      :: found
    INTEGER, INTENT(OUT), OPTIONAL            :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_p_d', &
      routineP = moduleN//':'//routineN

    LOGICAL                                   :: tr
    TYPE(cp_error_type)                       :: error

    CALL cp_error_init (error)
    CALL dbcsr_get_block_p(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)
    CALL cp_assert(.NOT. tr, cp_fatal_level, cp_internal_error,&
         routineN, "CP2K does not handle transposed blocks.")

  END SUBROUTINE cp_dbcsr_get_block_p_d


  SUBROUTINE cp_dbcsr_trace_a_d (matrix_a, trace, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(kind=real_8), INTENT(OUT)                     :: trace
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_trace_a_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_scalar_type)                  :: trace_scalar

    CALL timeset(routineN, timing_handle)
    trace_scalar = dbcsr_scalar_zero (cp_dbcsr_get_data_type(matrix_a))
    CALL dbcsr_trace(matrix_a%matrix, trace_scalar, dbcsr_error)
    CALL dbcsr_scalar_fill_all (trace_scalar)
    CALL dbcsr_scalar_get_value (trace_scalar, trace)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_trace_a_d


  SUBROUTINE cp_dbcsr_trace_ab_d (matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a, matrix_b
    REAL(kind=real_8), INTENT(INOUT)                   :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    LOGICAL, INTENT(IN), OPTIONAL            :: local_sum
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_trace_ab_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_trace(matrix_a%matrix, matrix_b%matrix, trace, trans_a, trans_b, local_sum, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_trace_ab_d


  SUBROUTINE cp_dbcsr_multiply_d (transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, &
       filter_eps,&
       error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(kind=real_8), INTENT(IN)                      :: alpha
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a, matrix_b
    REAL(kind=real_8), INTENT(IN)                      :: beta
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(kind=real_8), INTENT(IN), OPTIONAL :: filter_eps
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(int_8), INTENT(OUT), OPTIONAL    :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_d', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: prnt = .FALSE., &
                                                verify = .FALSE.

    CHARACTER(LEN=1)                         :: shape_a, shape_b, trans_a, &
                                                trans_b
    INTEGER                                  :: timing_handle, &
                                                timing_handle_detail
    LOGICAL                                  :: new_a_is_new, new_b_is_new
    REAL(kind=real_8)                       :: cs_b, cs_c
    TYPE(cp_dbcsr_type)                      :: new_a, new_b
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    !
    trans_a = transa
    trans_b = transb
    CALL uppercase(trans_a)
    CALL uppercase(trans_b)
    shape_a='R'
    IF(cp_dbcsr_nfullcols_total(matrix_a).EQ.cp_dbcsr_nfullrows_total(matrix_a)) shape_a='S'
    shape_b='R'
    IF(cp_dbcsr_nfullcols_total(matrix_b).EQ.cp_dbcsr_nfullrows_total(matrix_b)) shape_b='S'
    CALL matrix_match_sizes (matrix_c, matrix_a, transa, matrix_b, transb,&
         new_a, new_b, new_a_is_new, new_b_is_new, error)
    CALL timeset('cp_dbcsr_mult_'//trans_a//shape_a//'_'&
         //trans_b//shape_b, timing_handle_detail)
    CALL dbcsr_multiply(transa, transb,&
         alpha, new_a%matrix, new_b%matrix, beta, matrix_c%matrix,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, &
         filter_eps=filter_eps,&
         error=dbcsr_error, flop=flop)
    IF (new_a_is_new) THEN
       CALL cp_dbcsr_release (new_a, error=error)
    ENDIF
    IF (new_b_is_new) THEN
       CALL cp_dbcsr_release (new_b, error=error)
    ENDIF
    IF (prnt) THEN
       CALL cp_dbcsr_print (matrix_c, matlab_format=.TRUE.,&
            variable_name="mpo", error=error)
    ENDIF
    IF (verify) cs_b = cp_dbcsr_checksum (matrix_c, error=error)
    CALL timestop(timing_handle_detail)

    CALL timestop(timing_handle)
    IF (verify) THEN
       WRITE(*,'(A,4(1X,E9.3))')routineN//" checksums", cs_c, cs_b,&
            cs_c-cs_b, ABS(cs_c-cs_b)/cs_b
       WRITE(*,*)routineN//" multiply type",&
            trans_a//shape_a//'_'&
            //trans_b//shape_b

       IF (ABS(cs_c-cs_b) .GT. 0.00001) STOP "Bad multiply"
    ENDIF
  END SUBROUTINE cp_dbcsr_multiply_d


  SUBROUTINE cp_dbcsr_scale_by_vector_d (matrix_a, alpha, side, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)        :: matrix_a
    REAL(kind=real_8), DIMENSION(:), INTENT(IN), TARGET :: alpha
    CHARACTER(LEN=*), INTENT(IN)              :: side
    TYPE(cp_error_type), INTENT(INOUT)        :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_by_vector_d ', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: timing_handle
    TYPE(dbcsr_error_type)                    :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale_by_vector(matrix_a%matrix, alpha, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_by_vector_d


  SUBROUTINE cp_dbcsr_scale_d (matrix_a, alpha_scalar, last_column, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(kind=real_8), INTENT(IN)                      :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a%matrix, alpha_scalar, last_column, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_d


  SUBROUTINE cp_dbcsr_scale_d_m (matrix_a, alpha_matrix, side, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)        :: matrix_a
    REAL(kind=real_8), DIMENSION(:), INTENT(IN), TARGET :: alpha_matrix
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL    :: side
    TYPE(dbcsr_error_type), INTENT(INOUT)     :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_d_m', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: timing_handle
    TYPE(cp_error_type)                       :: cp_error
    TYPE(dbcsr_error_type)                    :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL cp_error_init (cp_error)
    dbcsr_error = error
    CALL dbcsr_scale_mat(matrix_a%matrix, alpha_matrix, side, error=dbcsr_error)
    error = dbcsr_error
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_d_m


  SUBROUTINE cp_dbcsr_set_d (matrix, alpha, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(kind=real_8), INTENT(IN)                      :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_set_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_set(matrix%matrix, cp_dbcsr_conform_scalar (alpha, matrix, error), dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_set_d


  SUBROUTINE cp_dbcsr_add_d (matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_b
    REAL(kind=real_8), INTENT(IN)                      :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_add(matrix_a%matrix, matrix_b%matrix, alpha_scalar, beta_scalar, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_d
