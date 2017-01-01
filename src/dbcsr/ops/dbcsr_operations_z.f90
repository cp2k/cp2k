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
  SUBROUTINE dbcsr_trace_a_z(matrix_a, trace)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix_a
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: trace

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_a_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_row, a_row_size, i, &
                                                mynode, error_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size,&
                                                row_dist, col_dist
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: a_data, data_p
    INTEGER, DIMENSION(:,:), POINTER         :: pgrid
    TYPE(dbcsr_distribution_obj)             :: dist

!   ---------------------------------------------------------------------------
    CALL timeset(routineN, error_handle)

    row_blk_size => array_data (matrix_a%row_blk_size)
    col_blk_size => array_data (matrix_a%col_blk_size)
    IF(dbcsr_get_data_type (matrix_a) /=  dbcsr_type_complex_8) &
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
             IF (mynode .NE. checker_square_proc (a_row, a_col, pgrid,&
                  row_dist, col_dist)) CYCLE
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
  END SUBROUTINE dbcsr_trace_a_z

! **************************************************************************************************
!> \brief traces a product of DBCSR matrices
!> \param[in] matrix_a DBCSR matrices
!> \param[in] matrix_b DBCSR matrices
!> \param[out] trace             the trace of the product of the matrices
! **************************************************************************************************
  SUBROUTINE dbcsr_trace_ab_z(matrix_a, matrix_b, trace)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix_a, matrix_b
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: trace

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_ab_z', &
      routineP = moduleN//':'//routineN

    INTEGER :: a_blk, a_col, a_col_size, a_row_size, b_blk, b_col_size, &
      b_frst_blk, b_last_blk, b_row_size, nze, row, a_beg, a_end, b_beg, b_end
    CHARACTER                                :: matrix_a_type, matrix_b_type
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, &
                                                a_row_blk_size, &
                                                b_col_blk_size, b_row_blk_size
    COMPLEX(kind=real_8)                                  :: sym_fac, fac
    LOGICAL                                  :: found, matrix_a_symm, matrix_b_symm
    REAL(real_4), EXTERNAL                   :: SDOT
    REAL(real_8), EXTERNAL                   :: DDOT
    COMPLEX(real_4), EXTERNAL                :: CDOTU
    COMPLEX(real_8), EXTERNAL                :: ZDOTU
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: a_data, b_data

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

  END SUBROUTINE dbcsr_trace_ab_z


! **************************************************************************************************
!> \brief Interface for matrix scaling by a scalar
!> \param matrix_a ...
!> \param alpha_scalar ...
!> \param last_column ...
! **************************************************************************************************
  SUBROUTINE dbcsr_scale_z(matrix_a, alpha_scalar, last_column)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix_a
    COMPLEX(kind=real_8), INTENT(IN)                      :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_z', &
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
  END SUBROUTINE dbcsr_scale_z

! **************************************************************************************************
!> \brief Interface for matrix scaling by a vector
!> \param matrix_a ...
!> \param alpha ...
!> \param side ...
! **************************************************************************************************
  SUBROUTINE dbcsr_scale_by_vector_z(matrix_a, alpha, side)
    TYPE(dbcsr_type), INTENT(INOUT)            :: matrix_a
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN), TARGET :: alpha
    CHARACTER(LEN=*), INTENT(IN)              :: side
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_by_vector_z', &
      routineP = moduleN//':'//routineN
    COMPLEX(kind=real_8), DIMENSION(:), POINTER            :: tmp_p
    TYPE(dbcsr_data_obj)                      :: enc_alpha_vec

    CALL dbcsr_data_init (enc_alpha_vec)
    CALL dbcsr_data_new (enc_alpha_vec, dbcsr_type_complex_8)
    tmp_p => alpha
    CALL dbcsr_data_set_pointer (enc_alpha_vec, tmp_p)
    CALL dbcsr_scale_by_vector_anytype(matrix_a, enc_alpha_vec, side)
    CALL dbcsr_data_clear_pointer (enc_alpha_vec)
    CALL dbcsr_data_release (enc_alpha_vec)
  END SUBROUTINE dbcsr_scale_by_vector_z


! **************************************************************************************************
!> \brief Interface for dbcsr_set
!> \param matrix ...
!> \param alpha ...
! **************************************************************************************************
  SUBROUTINE dbcsr_set_z(matrix, alpha)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix
    COMPLEX(kind=real_8), INTENT(IN)                      :: alpha
    IF (alpha==CMPLX(0.0, 0.0, real_8)) THEN
      CALL dbcsr_zero(matrix)
    ELSE
      CALL dbcsr_set_anytype(matrix, dbcsr_scalar(alpha))
    ENDIF
  END SUBROUTINE dbcsr_set_z

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eps ...
!> \param method ...
!> \param use_absolute ...
!> \param filter_diag ...
!> \param quick ...
! **************************************************************************************************
  SUBROUTINE dbcsr_filter_z (matrix, eps, method, use_absolute, &
       filter_diag, quick)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix
    COMPLEX(kind=real_8), INTENT(IN)                      :: eps
    INTEGER, INTENT(IN), OPTIONAL            :: method
    LOGICAL, INTENT(in), OPTIONAL            :: use_absolute, filter_diag, quick
    CALL dbcsr_filter_anytype (matrix, dbcsr_scalar(eps), method, &
         use_absolute, filter_diag, quick)
  END SUBROUTINE dbcsr_filter_z

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
  SUBROUTINE dbcsr_set_diag_z(matrix, diag)
    TYPE(dbcsr_type), INTENT(INOUT)            :: matrix
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN), TARGET :: diag

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_diag_z', &
      routineP = moduleN//':'//routineN

    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: diag_p
    TYPE(dbcsr_data_obj)                     :: diag_a

    diag_p => diag
    CALL dbcsr_data_init (diag_a)
    CALL dbcsr_data_new (diag_a, dbcsr_get_data_type(matrix))
    CALL dbcsr_data_set_pointer (diag_a, diag_p)
    CALL dbcsr_set_diag(matrix, diag_a)
    CALL dbcsr_data_clear_pointer (diag_a)
    CALL dbcsr_data_release (diag_a)
  END SUBROUTINE dbcsr_set_diag_z

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
  SUBROUTINE dbcsr_get_diag_z(matrix, diag)

    TYPE(dbcsr_type), INTENT(IN)                    :: matrix
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(INOUT), TARGET   :: diag

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_diag_z', &
      routineP = moduleN//':'//routineN

    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: diag_p
    TYPE(dbcsr_data_obj)                     :: diag_a

    diag_p => diag
    CALL dbcsr_data_init (diag_a)
    CALL dbcsr_data_new (diag_a, dbcsr_get_data_type(matrix))
    CALL dbcsr_data_set_pointer (diag_a, diag_p)
    CALL dbcsr_get_diag(matrix, diag_a)
    CALL dbcsr_data_clear_pointer (diag_a)
    CALL dbcsr_data_release (diag_a)
  END SUBROUTINE dbcsr_get_diag_z


! **************************************************************************************************
!> \brief add a constant to the diagonal of a matrix
!> \param[inout] matrix       DBCSR matrix
!> \param[in]    alpha_scalar scalar
!> \param first_row ...
!> \param last_row ...
! **************************************************************************************************
  SUBROUTINE dbcsr_add_on_diag_z(matrix, alpha_scalar, first_row, last_row)
    TYPE(dbcsr_type), INTENT(INOUT)           :: matrix
    COMPLEX(kind=real_8), INTENT(IN)                      :: alpha_scalar
    integer, intent(in), optional            :: first_row, last_row

    CALL dbcsr_add_on_diag(matrix, dbcsr_scalar(alpha_scalar), first_row, last_row)
  END SUBROUTINE dbcsr_add_on_diag_z
