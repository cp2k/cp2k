!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief traces a DBCSR matrix
!> \param[in] matrix_a       DBCSR matrix
!> \param[out] trace         the trace of the matrix
!>
! **************************************************************************************************
  SUBROUTINE dbcsr_trace_a_s(matrix_a, trace)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(kind=real_4), INTENT(INOUT)                   :: trace

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_a_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_row, a_row_size, i, &
                                                mynode, error_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size,&
                                                row_dist, col_dist
    REAL(kind=real_4), DIMENSION(:), POINTER           :: a_data, data_p
    INTEGER, DIMENSION(:,:), POINTER         :: pgrid
    TYPE(dbcsr_distribution_obj)             :: dist

!   ---------------------------------------------------------------------------
    CALL timeset(routineN, error_handle)

    row_blk_size => array_data (matrix_a%m%row_blk_size)
    col_blk_size => array_data (matrix_a%m%col_blk_size)
    CALL dbcsr_assert (dbcsr_get_data_type (matrix_a), "EQ", dbcsr_type_real_4,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Incompatible data types", __LINE__)
    CALL dbcsr_get_data (matrix_a%m%data_area, data_p)
    dist = dbcsr_distribution (matrix_a)
    mynode = dbcsr_mp_mynode (dbcsr_distribution_mp (dist))
    pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (dist))
    row_dist => dbcsr_distribution_row_dist (dist)
    col_dist => dbcsr_distribution_col_dist (dist)
    !
    ! let's go
    trace = REAL(0.0,real_4)
    DO a_row = 1, matrix_a%m%nblkrows_total
       a_row_size = row_blk_size(a_row)
       DO a_blk = matrix_a%m%row_p(a_row)+1,matrix_a%m%row_p(a_row+1)
          IF (a_blk .EQ. 0) CYCLE
          a_col = matrix_a%m%col_i(a_blk)
          IF(a_col.ne.a_row) CYCLE
          ! We must skip non-local blocks in a replicated matrix.
          IF(matrix_a%m%replication_type .NE. dbcsr_repl_full) THEN
             IF (mynode .NE. checker_square_proc (a_row, a_col, pgrid,&
                  row_dist, col_dist)) CYCLE
          ENDIF
          a_col_size = col_blk_size(a_col)
          IF(a_row_size.NE.a_col_size)&
             CPABORT("is that a square matrix?")
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

    CALL timestop(error_handle)
  END SUBROUTINE dbcsr_trace_a_s

! **************************************************************************************************
!> \brief traces a product of DBCSR matrices
!> \param[in] matrix_a DBCSR matrices
!> \param[in] matrix_b DBCSR matrices
!> \param[out] trace             the trace of the product of the matrices
!> \param[in] trans_a            (optional) is matrix_a transposed or not?
!> \param[in] trans_b            (optional) is matrix_b transposed or not?
!> \param local_sum ...
! **************************************************************************************************
  SUBROUTINE dbcsr_trace_ab_s(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a, matrix_b
    REAL(kind=real_4), INTENT(INOUT)                   :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    LOGICAL, INTENT(IN), OPTIONAL            :: local_sum

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_ab_s', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: my_trans_a, my_trans_b
    INTEGER :: a_blk, a_col, a_col_size, a_row_size, b_blk, b_col_size, &
      b_frst_blk, b_last_blk, b_row_size, nze, row, a_beg, a_end, b_beg, b_end
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, &
                                                a_row_blk_size, &
                                                b_col_blk_size, b_row_blk_size
    REAL(kind=real_4)                                  :: sym_fac, fac
    LOGICAL                                  :: found, my_local_sum
    REAL(real_4), EXTERNAL                   :: SDOT
    REAL(real_8), EXTERNAL                   :: DDOT
    COMPLEX(real_4), EXTERNAL                :: CDOTU
    COMPLEX(real_8), EXTERNAL                :: ZDOTU
    REAL(real_4), DIMENSION(:), POINTER      :: a_data_r, b_data_r
    REAL(real_8), DIMENSION(:), POINTER      :: a_data_d, b_data_d
    COMPLEX(real_4), DIMENSION(:), POINTER   :: a_data_c, b_data_c
    COMPLEX(real_8), DIMENSION(:), POINTER   :: a_data_z, b_data_z

!   ---------------------------------------------------------------------------

    my_local_sum = .FALSE.
    IF(PRESENT(local_sum)) my_local_sum=local_sum
    IF(.NOT.my_local_sum) THEN!fixme
       IF(matrix_a%m%replication_type .NE. dbcsr_repl_none&
            .OR. matrix_b%m%replication_type .NE. dbcsr_repl_none)&
            CPABORT("Trace of product of replicated matrices not yet possible.")
    ENDIF

    sym_fac = REAL(1.0,real_4)
    IF((dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric.OR.&
        dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric).AND.&
       (dbcsr_get_matrix_type(matrix_b).EQ.dbcsr_type_symmetric.OR.&
        dbcsr_get_matrix_type(matrix_b).EQ.dbcsr_type_antisymmetric)) sym_fac = REAL(2.0,real_4)

    ! tracing a symmetric with a general matrix is not implemented, as it would require communication of blocks
    IF((dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric.OR.&
        dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric).NEQV. &
       (dbcsr_get_matrix_type(matrix_b).EQ.dbcsr_type_symmetric.OR.&
        dbcsr_get_matrix_type(matrix_b).EQ.dbcsr_type_antisymmetric)) THEN
       CPABORT("Tracing general with symmetric matrix NYI")
    ENDIF

    a_row_blk_size => array_data (matrix_a%m%row_blk_size)
    a_col_blk_size => array_data (matrix_a%m%col_blk_size)
    b_row_blk_size => array_data (matrix_b%m%row_blk_size)
    b_col_blk_size => array_data (matrix_b%m%col_blk_size)
    CALL dbcsr_get_data (matrix_a%m%data_area, a_data_r)
    CALL dbcsr_get_data (matrix_b%m%data_area, b_data_r)
    CALL dbcsr_get_data (matrix_a%m%data_area, a_data_d)
    CALL dbcsr_get_data (matrix_b%m%data_area, b_data_d)
    CALL dbcsr_get_data (matrix_a%m%data_area, a_data_c)
    CALL dbcsr_get_data (matrix_b%m%data_area, b_data_c)
    CALL dbcsr_get_data (matrix_a%m%data_area, a_data_z)
    CALL dbcsr_get_data (matrix_b%m%data_area, b_data_z)

    my_trans_a = 'T'
    IF(PRESENT(trans_a)) my_trans_a = trans_a
    my_trans_b = 'N'
    IF(PRESENT(trans_b)) my_trans_b = trans_b
    IF(my_trans_a.NE.'T'.OR.my_trans_b.NE.'N')&
       CPABORT("this combination of transpose is NYI")
    !
    ! let's go
    trace = REAL(0.0,real_4)
    IF(matrix_a%m%nblkrows_total.NE.matrix_b%m%nblkrows_total)&
       CPABORT("this combination of transpose is NYI")
    DO row = 1, matrix_a%m%nblkrows_total
       a_row_size = a_row_blk_size(row)
       b_row_size = b_row_blk_size(row)
       IF(a_row_size.NE.b_row_size) THEN
          WRITE(*,*) 'a_row_size',a_row_size
          WRITE(*,*) 'b_row_size',b_row_size
          CPABORT("matrices not consistent")
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
                CPABORT("matrices not consistent")
             ENDIF
             !
             nze = a_row_size*a_col_size
             !
             IF(nze.GT.0) THEN
                !
                ! let's trace the blocks
                a_beg = ABS(matrix_a%m%blk_p(a_blk))
                a_end = a_beg + nze - 1
                b_beg = ABS(matrix_b%m%blk_p(b_blk))
                b_end = b_beg + nze - 1
                fac = REAL(1.0,real_4)
                IF(row.NE.a_col) fac = sym_fac
                ! this is a mess, no need to support that..
                IF(    matrix_a%m%data_type.EQ.dbcsr_type_real_4.AND.&
                       matrix_b%m%data_type.EQ.dbcsr_type_real_4) THEN
                   trace = trace +REAL(fac * SDOT (nze,&
                        a_data_r(ABS(matrix_a%m%blk_p(a_blk))),1,&
                        b_data_r(ABS(matrix_b%m%blk_p(b_blk))),1),kind=real_4)
                ELSEIF(matrix_a%m%data_type.EQ.dbcsr_type_real_4.AND.&
                     matrix_b%m%data_type.EQ.dbcsr_type_real_8) THEN
                   trace = trace + &
                        REAL(fac * SUM ( a_data_r(a_beg:a_end) * b_data_d(b_beg:b_end) ),kind=real_4)
                ELSEIF(matrix_a%m%data_type.EQ.dbcsr_type_real_8.AND.&
                       matrix_b%m%data_type.EQ.dbcsr_type_real_4) THEN
                   trace = trace + &
                        REAL(fac * SUM ( a_data_d(a_beg:a_end) * b_data_r(b_beg:b_end) ),kind=real_4)
                ELSEIF(matrix_a%m%data_type.EQ.dbcsr_type_real_8.AND.&
                       matrix_b%m%data_type.EQ.dbcsr_type_real_8) THEN
                   trace = trace + REAL(fac * DDOT (nze,&
                        a_data_d(ABS(matrix_a%m%blk_p(a_blk))),1,&
                        b_data_d(ABS(matrix_b%m%blk_p(b_blk))),1),kind=real_4)
                ELSEIF(matrix_a%m%data_type.EQ.dbcsr_type_complex_4.AND.&
                       matrix_b%m%data_type.EQ.dbcsr_type_complex_4) THEN
                   trace = trace + REAL(fac * CDOTU (nze,&
                        a_data_c(ABS(matrix_a%m%blk_p(a_blk))),1,&
                        b_data_c(ABS(matrix_b%m%blk_p(b_blk))),1),kind=real_4)
                ELSEIF(matrix_a%m%data_type.EQ.dbcsr_type_complex_8.AND.&
                       matrix_b%m%data_type.EQ.dbcsr_type_complex_8) THEN
                   trace = trace + REAL(fac * ZDOTU (nze,&
                        a_data_z(ABS(matrix_a%m%blk_p(a_blk))),1,&
                        b_data_z(ABS(matrix_b%m%blk_p(b_blk))),1),kind=real_4)
                ELSE
                   CPABORT("combination of types NYI")
                ENDIF
             ENDIF
          ENDIF
       ENDDO ! a_col
    ENDDO ! a_row
    !
    ! summe
    IF(.NOT.my_local_sum) &
           CALL mp_sum(trace,dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%m%dist)))

  END SUBROUTINE dbcsr_trace_ab_s


! **************************************************************************************************
!> \brief Interface for matrix scaling by a scalar
!> \param matrix_a ...
!> \param alpha_scalar ...
!> \param last_column ...
! **************************************************************************************************
  SUBROUTINE dbcsr_scale_s(matrix_a, alpha_scalar, last_column)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(kind=real_4), INTENT(IN)                      :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_s', &
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
  END SUBROUTINE dbcsr_scale_s

! **************************************************************************************************
!> \brief Interface for matrix scaling by a vector
!> \param matrix_a ...
!> \param alpha ...
!> \param side ...
! **************************************************************************************************
  SUBROUTINE dbcsr_scale_by_vector_s(matrix_a, alpha, side)
    TYPE(dbcsr_obj), INTENT(INOUT)            :: matrix_a
    REAL(kind=real_4), DIMENSION(:), INTENT(IN), TARGET :: alpha
    CHARACTER(LEN=*), INTENT(IN)              :: side
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_by_vector_s', &
      routineP = moduleN//':'//routineN
    REAL(kind=real_4), DIMENSION(:), POINTER            :: tmp_p
    TYPE(dbcsr_data_obj)                      :: enc_alpha_vec

    CALL dbcsr_data_init (enc_alpha_vec)
    CALL dbcsr_data_new (enc_alpha_vec, dbcsr_type_real_4)
    tmp_p => alpha
    CALL dbcsr_data_set_pointer (enc_alpha_vec, tmp_p)
    CALL dbcsr_scale_by_vector_anytype(matrix_a, enc_alpha_vec, side)
    CALL dbcsr_data_clear_pointer (enc_alpha_vec)
    CALL dbcsr_data_release (enc_alpha_vec)
  END SUBROUTINE dbcsr_scale_by_vector_s


! **************************************************************************************************
!> \brief Interface for dbcsr_set
!> \param matrix ...
!> \param alpha ...
! **************************************************************************************************
  SUBROUTINE dbcsr_set_s(matrix, alpha)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(kind=real_4), INTENT(IN)                      :: alpha
    IF (alpha==0.0_real_4) THEN
      CALL dbcsr_zero(matrix)
    ELSE
      CALL dbcsr_set_anytype(matrix, dbcsr_scalar(alpha))
    ENDIF
  END SUBROUTINE dbcsr_set_s

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eps ...
!> \param method ...
!> \param use_absolute ...
!> \param filter_diag ...
!> \param quick ...
! **************************************************************************************************
  SUBROUTINE dbcsr_filter_s (matrix, eps, method, use_absolute, &
       filter_diag, quick)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(kind=real_4), INTENT(IN)                      :: eps
    INTEGER, INTENT(IN), OPTIONAL            :: method
    LOGICAL, INTENT(in), OPTIONAL            :: use_absolute, filter_diag, quick
    CALL dbcsr_filter_anytype (matrix, dbcsr_scalar(eps), method, &
         use_absolute, filter_diag, quick)
  END SUBROUTINE dbcsr_filter_s

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
  SUBROUTINE dbcsr_set_diag_s(matrix, diag)
    TYPE(dbcsr_obj), INTENT(INOUT)            :: matrix
    REAL(kind=real_4), DIMENSION(:), INTENT(IN), TARGET :: diag

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_diag_s', &
      routineP = moduleN//':'//routineN

    REAL(kind=real_4), DIMENSION(:), POINTER           :: diag_p
    TYPE(dbcsr_data_obj)                     :: diag_a

    diag_p => diag
    CALL dbcsr_data_init (diag_a)
    CALL dbcsr_data_new (diag_a, dbcsr_get_data_type(matrix))
    CALL dbcsr_data_set_pointer (diag_a, diag_p)
    CALL dbcsr_set_diag(matrix, diag_a)
    CALL dbcsr_data_clear_pointer (diag_a)
    CALL dbcsr_data_release (diag_a)
  END SUBROUTINE dbcsr_set_diag_s

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
  SUBROUTINE dbcsr_get_diag_s(matrix, diag)

    TYPE(dbcsr_obj), INTENT(IN)                    :: matrix
    REAL(kind=real_4), DIMENSION(:), INTENT(INOUT), TARGET   :: diag

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_diag_s', &
      routineP = moduleN//':'//routineN

    REAL(kind=real_4), DIMENSION(:), POINTER           :: diag_p
    TYPE(dbcsr_data_obj)                     :: diag_a

    diag_p => diag
    CALL dbcsr_data_init (diag_a)
    CALL dbcsr_data_new (diag_a, dbcsr_get_data_type(matrix))
    CALL dbcsr_data_set_pointer (diag_a, diag_p)
    CALL dbcsr_get_diag(matrix, diag_a)
    CALL dbcsr_data_clear_pointer (diag_a)
    CALL dbcsr_data_release (diag_a)
  END SUBROUTINE dbcsr_get_diag_s


! **************************************************************************************************
!> \brief add a constant to the diagonal of a matrix
!> \param[inout] matrix       DBCSR matrix
!> \param[in]    alpha_scalar scalar
!> \param first_row ...
!> \param last_row ...
! **************************************************************************************************
  SUBROUTINE dbcsr_add_on_diag_s(matrix, alpha_scalar, first_row, last_row)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(kind=real_4), INTENT(IN)                      :: alpha_scalar
    integer, intent(in), optional            :: first_row, last_row

    CALL dbcsr_add_on_diag(matrix, dbcsr_scalar(alpha_scalar), first_row, last_row)
  END SUBROUTINE dbcsr_add_on_diag_s
