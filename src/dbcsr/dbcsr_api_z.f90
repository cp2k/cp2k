!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Encapsulates a given scalar value and makes it conformant to the
!>        type of the matrix.
!> \param scalar ...
!> \param matrix ...
!> \retval encapsulated ...
! **************************************************************************************************
  FUNCTION make_conformant_scalar_z (scalar, matrix) RESULT (encapsulated)
    COMPLEX(kind=real_8), INTENT(IN)                      :: scalar
    TYPE(dbcsr_type), INTENT(IN)             :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'make_conformant_scalar_z', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_scalar_type)                  :: encapsulated
    INTEGER                                  :: data_type, scalar_data_type

    encapsulated = dbcsr_scalar (scalar)
    CALL dbcsr_scalar_fill_all (encapsulated)
    data_type = dbcsr_get_data_type_prv(matrix%prv)
    scalar_data_type = dbcsr_scalar_get_type(encapsulated)
    IF (scalar_data_type .EQ. dbcsr_type_complex_4 .OR.&
        scalar_data_type .EQ. dbcsr_type_complex_8) THEN
       IF(data_type .NE. dbcsr_type_complex_4 .AND. data_type .NE. dbcsr_type_complex_8)&
          CPABORT("Can not conform a complex to a real number")
    END IF
    CALL dbcsr_scalar_set_type (encapsulated,data_type)
  END FUNCTION make_conformant_scalar_z


! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param transposed ...
!> \param existed ...
! **************************************************************************************************
  SUBROUTINE dbcsr_reserve_block2d_z (matrix, row, col, block, transposed, existed)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:, :), POINTER        :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_block2d_z', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_reserve_block2d_prv(matrix%prv, row, col, block,&
         transposed, existed)

  END SUBROUTINE dbcsr_reserve_block2d_z


! **************************************************************************************************
!> \brief ...
!> \param iterator ...
!> \param row ...
!> \param column ...
!> \param block ...
!> \param block_number ...
!> \param row_size ...
!> \param col_size ...
!> \param row_offset ...
!> \param col_offset ...
! **************************************************************************************************
  SUBROUTINE dbcsr_iterator_next_2d_block_z (iterator, row, column, block,&
       block_number, row_size, col_size, row_offset, col_offset)
    TYPE(dbcsr_iterator_type), INTENT(INOUT) :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    COMPLEX(kind=real_8), DIMENSION(:, :), POINTER        :: block
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number, row_size, &
                                                col_size, row_offset, &
                                                col_offset

    LOGICAL                                  :: transposed

    CALL dbcsr_iterator_next_block_prv (iterator%prv, row, column, block, transposed,&
       block_number, row_size, col_size, row_offset, col_offset)
    IF(transposed)&
       CPABORT("CP2K does not handle transposed blocks.")

  END SUBROUTINE dbcsr_iterator_next_2d_block_z


! **************************************************************************************************
!> \brief ...
!> \param iterator ...
!> \param row ...
!> \param column ...
!> \param block ...
!> \param block_number ...
!> \param row_size ...
!> \param col_size ...
!> \param row_offset ...
!> \param col_offset ...
! *************************************************************************************************
  SUBROUTINE dbcsr_iterator_next_1d_block_z (iterator, row, column, block,&
       block_number, row_size, col_size, row_offset, col_offset)
    TYPE(dbcsr_iterator_type), INTENT(INOUT)  :: iterator
    INTEGER, INTENT(OUT)                      :: row, column
    COMPLEX(kind=real_8), DIMENSION(:), POINTER            :: block
    INTEGER, INTENT(OUT), OPTIONAL            :: block_number, row_size, &
                                                 col_size, row_offset, &
                                                 col_offset

    LOGICAL                                   :: transposed

    CALL dbcsr_iterator_next_block_prv(iterator%prv, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)
    IF(transposed)&
       CPABORT("CP2K does not handle transposed blocks.")

  END SUBROUTINE dbcsr_iterator_next_1d_block_z


! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param summation ...
!> \param scale ...
! **************************************************************************************************
  SUBROUTINE dbcsr_put_block2d_z (matrix, row, col, block,&
       summation, scale)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:, :), INTENT(IN)     :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: summation
    COMPLEX(kind=real_8), INTENT(IN), OPTIONAL            :: scale

    CALL dbcsr_put_block_prv(matrix%prv, row, col, block, summation=summation, scale=scale)
  END SUBROUTINE dbcsr_put_block2d_z


! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param summation ...
!> \param scale ...
! **************************************************************************************************
  SUBROUTINE dbcsr_put_block_z (matrix, row, col, block,&
       summation, scale)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN)        :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: summation
    COMPLEX(kind=real_8), INTENT(IN), OPTIONAL            :: scale

    CALL dbcsr_put_block_prv(matrix%prv, row, col, block, summation=summation, scale=scale)
  END SUBROUTINE dbcsr_put_block_z


! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param found ...
!> \param row_size ...
!> \param col_size ...
! **************************************************************************************************
  SUBROUTINE dbcsr_get_2d_block_p_z (matrix,row,col,block,found, row_size, col_size)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:, :), POINTER        :: block
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    LOGICAL                                  :: tr

    CALL dbcsr_get_block_p_prv(matrix%prv,row,col,block,tr,found, row_size, col_size)
    IF(tr)&
       CPABORT("CP2K does not handle transposed blocks.")
  END SUBROUTINE dbcsr_get_2d_block_p_z


! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param found ...
!> \param row_size ...
!> \param col_size ...
! **************************************************************************************************
  SUBROUTINE dbcsr_get_block_p_z (matrix,row,col,block,found, row_size, col_size)
    TYPE(dbcsr_type), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                       :: row, col
    COMPLEX(kind=real_8), DIMENSION(:), POINTER            :: block
    LOGICAL, INTENT(OUT)                      :: found
    INTEGER, INTENT(OUT), OPTIONAL            :: row_size, col_size

    LOGICAL                                   :: tr

    CALL dbcsr_get_block_p_prv(matrix%prv,row,col,block,tr,found, row_size, col_size)
    IF(tr)&
       CPABORT("CP2K does not handle transposed blocks.")

  END SUBROUTINE dbcsr_get_block_p_z


! **************************************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param trace ...
! **************************************************************************************************
  SUBROUTINE dbcsr_trace_a_z (matrix_a, trace)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a
    COMPLEX(kind=real_8), INTENT(OUT)                     :: trace

    TYPE(dbcsr_scalar_type)                  :: trace_scalar

    trace_scalar = dbcsr_scalar_zero (dbcsr_get_data_type(matrix_a))
    CALL dbcsr_trace_prv(matrix_a%prv, trace_scalar)
    CALL dbcsr_scalar_fill_all(trace_scalar)
    CALL dbcsr_scalar_get_value(trace_scalar, trace)
  END SUBROUTINE dbcsr_trace_a_z


! **************************************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param trace ...
! **************************************************************************************************
  SUBROUTINE dbcsr_trace_ab_z (matrix_a, matrix_b, trace)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a, matrix_b
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: trace

    CALL dbcsr_trace_prv(matrix_a%prv, matrix_b%prv, trace)
  END SUBROUTINE dbcsr_trace_ab_z


! **************************************************************************************************
!> \brief ...
!> \param transa ...
!> \param transb ...
!> \param alpha ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param beta ...
!> \param matrix_c ...
!> \param first_row ...
!> \param last_row ...
!> \param first_column ...
!> \param last_column ...
!> \param first_k ...
!> \param last_k ...
!> \param retain_sparsity ...
!> \param match_matrix_sizes Enables BLAS XGEMM-style multiplication
!>        of matrices with incompatible dimensions. By default it's disabled.
!> \param filter_eps ...
!> \param flop ...
! **************************************************************************************************
  SUBROUTINE dbcsr_multiply_z (transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, filter_eps, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    COMPLEX(kind=real_8), INTENT(IN)                      :: alpha
    TYPE(dbcsr_type), INTENT(IN)             :: matrix_a, matrix_b
    COMPLEX(kind=real_8), INTENT(IN)                      :: beta
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(kind=real_8), INTENT(IN), OPTIONAL :: filter_eps
    INTEGER(int_8), INTENT(OUT), OPTIONAL    :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_z', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: trans_a, trans_b

    trans_a = transa
    trans_b = transb
    CALL uppercase(trans_a)
    CALL uppercase(trans_b)

    CALL dbcsr_multiply_prv(transa, transb,&
         alpha, matrix_a%prv, matrix_b%prv, beta, matrix_c%prv,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, &
         filter_eps=filter_eps,&
         flop=flop)

  END SUBROUTINE dbcsr_multiply_z


! **************************************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param alpha ...
!> \param side ...
! **************************************************************************************************
  SUBROUTINE dbcsr_scale_by_vector_z (matrix_a, alpha, side)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix_a
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN), TARGET :: alpha
    CHARACTER(LEN=*), INTENT(IN)              :: side

    CALL dbcsr_scale_by_vector_prv(matrix_a%prv, alpha, side)
  END SUBROUTINE dbcsr_scale_by_vector_z


! **************************************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param alpha_scalar ...
!> \param last_column ...
! **************************************************************************************************
  SUBROUTINE dbcsr_scale_z (matrix_a, alpha_scalar, last_column)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a
    COMPLEX(kind=real_8), INTENT(IN)                      :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column

    CALL dbcsr_scale_prv(matrix_a%prv, alpha_scalar, last_column)
  END SUBROUTINE dbcsr_scale_z


! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param alpha ...
! **************************************************************************************************
  SUBROUTINE dbcsr_set_z (matrix, alpha)
    TYPE(dbcsr_type), INTENT(INOUT)       :: matrix
    COMPLEX(kind=real_8), INTENT(IN)                      :: alpha

    CALL dbcsr_set_prv(matrix%prv, alpha)
  END SUBROUTINE dbcsr_set_z


! **************************************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param alpha_scalar ...
!> \param beta_scalar ...
! **************************************************************************************************
  SUBROUTINE dbcsr_add_z (matrix_a, matrix_b, alpha_scalar, beta_scalar)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a
    TYPE(dbcsr_type), INTENT(IN)             :: matrix_b
    COMPLEX(kind=real_8), INTENT(IN)                      :: alpha_scalar, beta_scalar

    CALL dbcsr_add_prv(matrix_a%prv, matrix_b%prv, alpha_scalar, beta_scalar)
  END SUBROUTINE dbcsr_add_z

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param alpha_scalar ...
!> \param first_row ...
!> \param last_row ...
! **************************************************************************************************
   SUBROUTINE dbcsr_add_on_diag_z(matrix, alpha_scalar)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      COMPLEX(kind=real_8), INTENT(IN)                                :: alpha_scalar

      CALL dbcsr_add_on_diag_prv(matrix%prv, alpha_scalar)
   END SUBROUTINE dbcsr_add_on_diag_z

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
   SUBROUTINE dbcsr_set_diag_z(matrix, diag)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN)                  :: diag

      CALL dbcsr_set_diag_prv(matrix%prv, diag)
   END SUBROUTINE dbcsr_set_diag_z

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_diag_z(matrix, diag)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      COMPLEX(kind=real_8), DIMENSION(:), INTENT(OUT)                 :: diag

      CALL dbcsr_get_diag_prv(matrix%prv, diag)
   END SUBROUTINE dbcsr_get_diag_z

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param index_matrix ...
!> \param lb ...
!> \param ub ...
!> \retval DATA ...
! **************************************************************************************************
  FUNCTION dbcsr_get_wms_data_z (matrix, index_matrix, select_data_type, lb, ub) RESULT (DATA)
    TYPE(dbcsr_type), INTENT(IN)     :: matrix
    INTEGER, INTENT(IN)              :: index_matrix
    COMPLEX(kind=real_8), INTENT(IN)              :: select_data_type
    COMPLEX(kind=real_8), DIMENSION(:), POINTER   :: DATA
    INTEGER, INTENT(IN), OPTIONAL    :: lb, ub

    DATA => dbcsr_get_data_p_prv(matrix%prv%wms(index_matrix)%data_area,select_data_type,lb,ub)

  END FUNCTION dbcsr_get_wms_data_z

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param index_matrix ...
!> \param lb ...
!> \param ub ...
!> \retval DATA ...
! **************************************************************************************************
  FUNCTION dbcsr_get_data_z (matrix, select_data_type, lb, ub) RESULT (DATA)
    TYPE(dbcsr_type), INTENT(IN)     :: matrix
    COMPLEX(kind=real_8), INTENT(IN)              :: select_data_type
    COMPLEX(kind=real_8), DIMENSION(:), POINTER   :: DATA
    INTEGER, INTENT(IN), OPTIONAL    :: lb, ub

    DATA => dbcsr_get_data_p_prv(matrix%prv%data_area,select_data_type,lb,ub)

  END FUNCTION dbcsr_get_data_z

