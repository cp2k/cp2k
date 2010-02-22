! *****************************************************************************
!> \brief traces a DBCSR matrix
!> \param[in] matrix_a       DBCSR matrix
!> \param[out] trace         the trace of the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_trace_a_d(matrix_a, trace, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(kind=real_8), INTENT(INOUT)                   :: trace
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_a_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_row, a_row_size, i, &
                                                mynode
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size,&
                                                row_dist, col_dist
    REAL(kind=real_8), DIMENSION(:), POINTER           :: a_data, data_p
    INTEGER, DIMENSION(:,:), POINTER         :: pgrid
    TYPE(dbcsr_distribution_obj)             :: dist

!   ---------------------------------------------------------------------------

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
    trace = REAL(0.0,real_8)
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
          CALL dbcsr_assert (a_row_size.EQ.a_col_size, dbcsr_fatal_level,&
               dbcsr_internal_error, routineN, "is that a square matrix?",__LINE__,error)
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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_ab_d', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: my_trans_a, my_trans_b
    INTEGER :: a_blk, a_col, a_col_size, a_row_size, b_blk, b_col_size, &
      b_frst_blk, b_last_blk, b_row_size, nze, row, a_beg, a_end, b_beg, b_end
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, &
                                                a_row_blk_size, &
                                                b_col_blk_size, b_row_blk_size
    REAL(kind=real_8)                                  :: sym_fac, fac
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
       CALL dbcsr_assert (matrix_a%m%replication_type .EQ. dbcsr_repl_none&
            .AND. matrix_b%m%replication_type .EQ. dbcsr_repl_none,&
            dbcsr_failure_level, dbcsr_unimplemented_error_nr, routineN,&
            "Trace of product of replicated matrices not yet possible.",__LINE__,error)
    ENDIF

    sym_fac = REAL(1.0,real_8)
    IF((dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric.OR.&
        dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric).AND.&
       (dbcsr_get_matrix_type(matrix_b).EQ.dbcsr_type_symmetric.OR.&
        dbcsr_get_matrix_type(matrix_b).EQ.dbcsr_type_antisymmetric)) sym_fac = REAL(2.0,real_8)

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
    CALL dbcsr_assert (my_trans_a.EQ.'T'.AND.my_trans_b.EQ.'N', dbcsr_fatal_level,&
         dbcsr_unimplemented_error_nr, routineN, "this combination of transpose is NYI",&
         __LINE__,error)
    !
    ! let's go
    trace = REAL(0.0,real_8)
    CALL dbcsr_assert (matrix_a%m%nblkrows_total.EQ.matrix_b%m%nblkrows_total,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN, "this combination of transpose is NYI",&
         __LINE__,error)
    DO row = 1, matrix_a%m%nblkrows_total
       a_row_size = a_row_blk_size(row)
       b_row_size = b_row_blk_size(row)
       IF(a_row_size.NE.b_row_size) THEN
          WRITE(*,*) 'a_row_size',a_row_size
          WRITE(*,*) 'b_row_size',b_row_size
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error, &
               routineN, "matrices not consistent",__LINE__,error)
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
                CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error, &
                     routineN, "matrices not consistent",__LINE__,error)
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
                fac = REAL(1.0,real_8)
                IF(row.NE.a_col) fac = sym_fac
                IF(    matrix_a%m%data_type.EQ.dbcsr_type_real_4.AND.&
                       matrix_b%m%data_type.EQ.dbcsr_type_real_4) THEN
                   trace = trace + fac * SDOT (nze,&
                        a_data_r(ABS(matrix_a%m%blk_p(a_blk))),1,&
                        b_data_r(ABS(matrix_b%m%blk_p(b_blk))),1)
                ELSEIF(matrix_a%m%data_type.EQ.dbcsr_type_real_4.AND.&
                     matrix_b%m%data_type.EQ.dbcsr_type_real_8) THEN
                   trace = trace + &
                        fac * SUM ( a_data_r(a_beg:a_end) * b_data_d(b_beg:b_end) )
                ELSEIF(matrix_a%m%data_type.EQ.dbcsr_type_real_8.AND.&
                       matrix_b%m%data_type.EQ.dbcsr_type_real_4) THEN
                   trace = trace + &
                        fac * SUM ( a_data_d(a_beg:a_end) * b_data_r(b_beg:b_end) )
                ELSEIF(matrix_a%m%data_type.EQ.dbcsr_type_real_8.AND.&
                       matrix_b%m%data_type.EQ.dbcsr_type_real_8) THEN
                   trace = trace + fac * DDOT (nze,&
                        a_data_d(ABS(matrix_a%m%blk_p(a_blk))),1,&
                        b_data_d(ABS(matrix_b%m%blk_p(b_blk))),1)
                ELSEIF(matrix_a%m%data_type.EQ.dbcsr_type_complex_4.AND.&
                       matrix_b%m%data_type.EQ.dbcsr_type_complex_4) THEN
                   trace = trace + fac * CDOTU (nze,&
                        a_data_c(ABS(matrix_a%m%blk_p(a_blk))),1,&
                        b_data_c(ABS(matrix_b%m%blk_p(b_blk))),1)
                ELSEIF(matrix_a%m%data_type.EQ.dbcsr_type_complex_8.AND.&
                       matrix_b%m%data_type.EQ.dbcsr_type_complex_8) THEN
                   trace = trace + fac * ZDOTU (nze,&
                        a_data_z(ABS(matrix_a%m%blk_p(a_blk))),1,&
                        b_data_z(ABS(matrix_b%m%blk_p(b_blk))),1)
                ELSE
                   CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                        routineN, "combination of types NYI",__LINE__,error)
                ENDIF
             ENDIF
          ENDIF
       ENDDO ! a_col
    ENDDO ! a_row
    !
    ! summe
    IF(.NOT.my_local_sum) &
         & CALL mp_sum(trace,dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%m%dist)))

  END SUBROUTINE dbcsr_trace_ab_d


  SUBROUTINE bcast_block_d(blk, source, scope, mp_obj, error)
    REAL(kind=real_8), DIMENSION(:), INTENT(INOUT)     :: blk
    INTEGER, INTENT(IN)                      :: source
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_error_type), INTENT(INOUT)       :: error

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
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
            routineN, "combination of types NYI",__LINE__,error)
    CASE DEFAULT
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "how do you wanna bcast",__LINE__,error)
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

    my_beta = REAL(0.0,real_8)
    my_alpha = REAL(0.0,real_8)
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
    TYPE(dbcsr_error_type)                   :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'block_chol_inv_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, info, j

!   ---------------------------------------------------------------------------

    CALL dpotrf( 'U', m, blk, m, info )
    CALL dbcsr_assert (info.EQ.0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "error in dpotrf",__LINE__,error)
    CALL dpotri( 'U', m, blk, m, info )
    CALL dbcsr_assert (info.EQ.0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "error in dpotri",__LINE__,error)
    !
    ! symmetrize
    DO i=1,m
       DO j=i,m
          blk(j,i) = blk(i,j)
       ENDDO
    ENDDO
  END SUBROUTINE block_chol_inv_d

  !> \brief Interface for matrix scaling by a vector
  SUBROUTINE dbcsr_scale_by_vector_d(matrix_a, alpha, side, error)
    TYPE(dbcsr_obj), INTENT(INOUT)            :: matrix_a
    REAL(kind=real_8), DIMENSION(:), INTENT(IN), TARGET :: alpha
    CHARACTER(LEN=*), INTENT(IN)              :: side
    TYPE(dbcsr_error_type), INTENT(INOUT)     :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_by_vector_d', &
      routineP = moduleN//':'//routineN
    REAL(kind=real_8), DIMENSION(:), POINTER            :: tmp_p
    CHARACTER                                 :: scale_type
    TYPE(dbcsr_data_obj)                      :: enc_alpha_vec

    scale_type = 'V'
    CALL dbcsr_data_init (enc_alpha_vec)
    CALL dbcsr_data_new (enc_alpha_vec, dbcsr_type_real_8)
    tmp_p => alpha
    CALL dbcsr_data_set_pointer (enc_alpha_vec, tmp_p)
    CALL dbcsr_scale_by_vector_anytype(matrix_a, enc_alpha_vec, side, error)
    CALL dbcsr_data_clear_pointer (enc_alpha_vec)
    CALL dbcsr_data_release (enc_alpha_vec)
  END SUBROUTINE dbcsr_scale_by_vector_d

  !> \brief Interface for matrix scaling by a matrix
  SUBROUTINE dbcsr_scale_d_m(matrix_a, alpha_matrix,&
       last_column, side, error)
    TYPE(dbcsr_obj), INTENT(INOUT)             :: matrix_a
    REAL(kind=real_8), DIMENSION(:), INTENT(IN), TARGET  :: alpha_matrix
    INTEGER, INTENT(IN), OPTIONAL              :: last_column
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL     :: side
    TYPE(dbcsr_error_type), INTENT(INOUT)      :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_d_m', &
      routineP = moduleN//':'//routineN
    REAL(kind=real_8), DIMENSION(:), POINTER             :: tmp_p
    CHARACTER                                  :: scale_type
    TYPE(dbcsr_data_obj)                       :: enc_alpha_mat

    scale_type = 'M'
    CALL dbcsr_data_init (enc_alpha_mat)
    CALL dbcsr_data_new (enc_alpha_mat, dbcsr_type_real_8)
    tmp_p => alpha_matrix
    CALL dbcsr_data_set_pointer (enc_alpha_mat, tmp_p)
    CALL dbcsr_scale_anytype(matrix_a, last_column, error, scale_type,&
         side=side, alpha_matrix = enc_alpha_mat)
    CALL dbcsr_data_clear_pointer (enc_alpha_mat)
    CALL dbcsr_data_release (enc_alpha_mat)
  END SUBROUTINE dbcsr_scale_d_m

  !> \brief Interface for dbcsr_set
  SUBROUTINE dbcsr_set_d(matrix, alpha, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(kind=real_8), INTENT(IN)                      :: alpha
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    CALL dbcsr_set_anytype(matrix, dbcsr_scalar(alpha), error)
  END SUBROUTINE dbcsr_set_d

  !> \brief Interface for dbcsr_add THIS IS NOT NEEDED ANYMORE (alpha and beta ALWAYS REQUESTED)
  SUBROUTINE dbcsr_add_alpha_d(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    REAL(kind=real_8), INTENT(IN)                      :: alpha_scalar
    REAL(kind=real_8), INTENT(IN), OPTIONAL            :: beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    IF (PRESENT(beta_scalar)) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSE
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar), error=error)
    ENDIF
  END SUBROUTINE dbcsr_add_alpha_d
  !> \brief Interface for dbcsr_add THIS IS NOT NEEDED ANYMORE (alpha and beta ALWAYS REQUESTED)
  SUBROUTINE dbcsr_add_beta_d(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    REAL(kind=real_8), INTENT(IN), OPTIONAL            :: alpha_scalar
    REAL(kind=real_8), INTENT(IN)                      :: beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    IF (PRESENT(alpha_scalar)) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSE
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ENDIF
  END SUBROUTINE dbcsr_add_beta_d


! *****************************************************************************
!> \brief Issues actual DGEMM calls.
!>
!> \param[in] params           Stack of DGEMM parameters
!> \param[in] n                Number of parameters
!> \param[in] left_data_a      Left-matrix data
!> \param[in] right_data_a     Right-matrix data
!> \param[in,out] product_data_area  Data for results
!> \param[out] lflop           (optional) Number of FLOPs used by DGEMM
! *****************************************************************************
  SUBROUTINE process_dgemm_stack_d(params,&
       params_a, params_b, params_c,&
       stack_size, stack_size_a, stack_size_b, stack_size_c,&
       left_data_a, right_data_a, product_data_area, lflop, ltime, error)
    INTEGER, INTENT(IN)                       :: stack_size, stack_size_a,&
                                                 stack_size_b, stack_size_c
    TYPE(dgemm_join), DIMENSION(1:stack_size), &
      INTENT(IN)                              :: params
    TYPE(block_parameters), DIMENSION(1:stack_size_a), INTENT(IN) :: params_a
    TYPE(block_parameters), DIMENSION(1:stack_size_b), INTENT(IN) :: params_b
    TYPE(block_parameters), DIMENSION(1:stack_size_c), INTENT(IN) :: params_c
    REAL(kind=real_8), DIMENSION(:), INTENT(IN)         :: left_data_a, &
                                                 right_data_a
    TYPE(dbcsr_data_obj), INTENT(INOUT)       :: product_data_area
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                                :: lflop
    REAL(KIND=real_8), INTENT(OUT), OPTIONAL  :: ltime
    TYPE(dbcsr_error_type), INTENT(inout)     :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_dgemm_stack_d', &
      routineP = moduleN//':'//routineN
    REAL, PARAMETER                           :: resize_factor = 1.618034

    INTEGER                                   :: c, lda, ldb, ldc, maxs, r, sp,&
                                                 poff
    INTEGER                                   :: sp_c, sp_a, sp_b
    LOGICAL                                   :: do_resize

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY

    IF (PRESENT (lflop)) lflop = INT(0, int_8)
    IF (PRESENT (ltime)) ltime = 0.0_real_8
    ! Increase product data area size if necessary.
    do_resize = .FALSE.
    maxs = dbcsr_get_data_size(product_data_area)
    DO sp = 1, stack_size
       sp_c = params(sp)%p_c
       IF (params_c(sp_c)%offset+params_c(sp_c)%nze-1 .GT. maxs) THEN
          maxs = params_c(sp_c)%offset+params_c(sp_c)%nze-1
          do_resize = .TRUE.
       ENDIF
    ENDDO
    IF (do_resize) THEN
       CALL dbcsr_data_ensure_size (product_data_area,&
            maxs, factor=resize_factor, error=error)
    ENDIF
    ! Setup encapsulated data area
    DO sp = 1, stack_size
       sp_c = params(sp)%p_c
       sp_a = params(sp)%p_a
       sp_b = params(sp)%p_b
       IF (params(sp)%last_k .EQ. 0) CYCLE
       IF (params_c(sp_c)%nze .EQ. 0&
            .or.params_a(sp_a)%nze .EQ. 0&
            .or.params_b(sp_b)%nze .EQ. 0) THEN
          CYCLE
       ENDIF
       !
       poff = params_c(sp_c)%offset
       !
       ! Avoid nagging compilers
       IF (params(sp)%beta%r_dp .EQ. 0.0_real_8) THEN
          IF (.NOT. params_c(sp_c)%tr &
               .OR. params(sp)%last_n .EQ. params_c(sp_c)%logical_cols) THEN
             CALL dbcsr_data_clear (product_data_area,&
                  lb=poff,&
                  ub=poff-1+params_c(sp_c)%logical_rows*params(sp)%last_n)
          ELSE
             ! It must be transposed.
             !### optimize out the inner loop (turn into a range)
             DO r = 1, params_c(sp_c)%logical_rows
                DO c = 1, params(sp)%last_n
                   CALL dbcsr_data_clear(product_data_area,&
                        lb=poff-1+(r-1)*params_c(sp_c)%logical_cols+c,&
                        ub=poff-1+(r-1)*params_c(sp_c)%logical_cols+c)
                ENDDO
             ENDDO
             !FORALL (r = 1:params(sp)%lrows_c, c = 1:params(sp)%last_n)
             !   product_data_d((r-1)*params(sp)%lcols_c+c) = 0.0_dp
             !END FORALL
          ENDIF
       ENDIF
       ! Now prepare and issue the multiplies.
       lda = select_n_or_t (params_a(sp_a)%tr,&
            params_a(sp_a)%logical_rows, params_a(sp_a)%logical_cols)
       ldb = select_n_or_t (params_b(sp_b)%tr,&
            params_b(sp_b)%logical_rows, params_b(sp_b)%logical_cols)
       ldc = select_n_or_t (params_c(sp_c)%tr,&
            params_c(sp_c)%logical_rows, params_c(sp_c)%logical_cols)
       !write(*,*)"LDx", lda, ldb, ldc
       IF (.NOT. params_c(sp_c)%tr) THEN
          IF (PRESENT (ltime)) ltime = ltime + m_walltime()
          CALL DGEMM(&
               blas_mat_type(params_a(sp_a)%tr),&
               blas_mat_type(params_b(sp_b)%tr),&
               params_c(sp_c)%logical_rows, params(sp)%last_n,& !m, n
               params(sp)%last_k,& ! k
               params(sp)%alpha%r_dp,&
               left_data_a(params_a(sp_a)%offset), lda,&
               right_data_a(params_b(sp_b)%offset), ldb,&
               params(sp)%beta%r_dp,&
               product_data_area%d%r_dp(poff), ldc)
          IF (PRESENT (ltime)) ltime = m_walltime() - ltime
          IF (PRESENT (lflop)) &
               lflop = lflop + INT(2, int_8) * &
               INT(params_c(sp_c)%logical_rows &
               * params(sp)%last_n * params(sp)%last_k, int_8)
       ELSE
          IF (PRESENT (ltime)) ltime = ltime + m_walltime()
          CALL DGEMM(&
               blas_mat_type(.NOT.params_b(sp_b)%tr),&
               blas_mat_type(.NOT.params_a(sp_a)%tr),&
               params_c(sp_c)%logical_cols, params_c(sp_c)%logical_rows,& !m, n (switched)
               params(sp)%last_k,& ! k
               params(sp)%alpha%r_dp,&
               right_data_a(params_b(sp_b)%offset), ldb,&
               left_data_a(params_a(sp_a)%offset), lda,&
               params(sp)%beta%r_dp,&
               product_data_area%d%r_dp(poff), ldc)
          IF (PRESENT (ltime)) ltime = m_walltime() - ltime
          IF (PRESENT (lflop)) &
               lflop = lflop + INT(2, int_8) * &
               INT(params(sp)%last_n &
               * params_c(sp_c)%logical_cols * params(sp)%last_k, int_8)
       ENDIF
    ENDDO
  END SUBROUTINE process_dgemm_stack_d


