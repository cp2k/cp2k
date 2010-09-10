!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Tests for DBCSR multiply
!> \author  VW
!> \date    2010
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2010
! *****************************************************************************
MODULE dbcsr_test_multiply
  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_release
  USE dbcsr_block_access
  USE dbcsr_block_operations,          ONLY: dbcsr_block_partial_copy,&
                                             dbcsr_data_clear
  USE dbcsr_data_methods
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_caller_error,&
                                             dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_internal_error,&
                                             dbcsr_wrong_args_error
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_kinds,                     ONLY: real_4,&
                                             real_8
  USE dbcsr_message_passing,           ONLY: mp_cart_create,&
                                             mp_cart_rank,&
                                             mp_comm_null,&
                                             mp_dims_create,&
                                             mp_environ,&
                                             mp_sync
  USE dbcsr_methods
  USE dbcsr_operations
  USE dbcsr_ptr_util
  USE dbcsr_test_methods,              ONLY: atoi,&
                                             atol,&
                                             ator,&
                                             dbcsr_make_random_block_sizes,&
                                             dbcsr_make_random_matrix,&
                                             dbcsr_random_dist
  USE dbcsr_transformations
  USE dbcsr_types
  USE dbcsr_util
  USE dbcsr_work_operations

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dbcsr_test_multiplies

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_test_multiply'

CONTAINS

  SUBROUTINE dbcsr_test_multiplies( group, ionode, nproc, narg, args, error)

    INTEGER                                  :: group, ionode
    INTEGER, DIMENSION(:), POINTER           :: nproc
    INTEGER                                  :: narg
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(IN)                             :: args
    TYPE(dbcsr_error_type)                   :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_test_multiplies', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: symmetries(3)
    INTEGER                                  :: limits(6), matrix_sizes(3)
    INTEGER, DIMENSION(:), POINTER           :: bs_k, bs_m, bs_n
    LOGICAL                                  :: retain_sparsity
    REAL(real_8)                             :: alpha(2), beta(2), &
                                                sparsities(3)

!
! parsing

    CALL dbcsr_assert( narg.GE.21 , dbcsr_fatal_level, dbcsr_wrong_args_error, &
         routineN, "narg not correct", __LINE__, error )

     matrix_sizes(1) = atoi(args( 2))
     matrix_sizes(2) = atoi(args( 3))
     matrix_sizes(3) = atoi(args( 4))
     sparsities(1)   = ator(args( 5))
     sparsities(2)   = ator(args( 6))
     sparsities(3)   = ator(args( 7))
     alpha(1)        = ator(args( 8))
     alpha(2)        = ator(args( 9))
     beta(1)         = ator(args(10))
     beta(2)         = ator(args(11))
     limits(1)       = atoi(args(12))
     limits(2)       = atoi(args(13))
     limits(3)       = atoi(args(14))
     limits(4)       = atoi(args(15))
     limits(5)       = atoi(args(16))
     limits(6)       = atoi(args(17))
     symmetries(1)   = args(18)
     symmetries(2)   = args(19)
     symmetries(3)   = args(20)
     retain_sparsity= atol(args(21))

     ! lets go !
     NULLIFY(bs_m, bs_n, bs_k,nproc)
     CALL dbcsr_test_mult( group, ionode, nproc, matrix_sizes, &
          bs_m, bs_n, bs_k, sparsities, symmetries, &
          alpha, beta, limits, retain_sparsity, error=error)

   END SUBROUTINE dbcsr_test_multiplies



! *****************************************************************************
!> \brief Performs a variety of matrix multiplies of same matrices on different
!>        processor grids
!> \param[in] mp_group          MPI communicator
!> \param[in] io_unit           which unit to write to, if not negative
!> \param[in] nproc             number of processors to test on
!> \param[in] matrix_sizes      size of matrices to test
!> \param[in] matrix_types      types of matrices to create
!> \param[in] trs               transposes of the two matrices
!> \param[in] bs_m, bs_n, bs_k  block sizes of the 3 dimensions
!> \param[in] sparsities        sparsities of matrices to create
!> \param[in] alpha, beta       alpha and beta values to use in multiply
!> \param[in,out] error         cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_test_mult (mp_group, io_unit, nproc,&
       matrix_sizes, bs_m, bs_n, bs_k, sparsities, symmetries, &
       alpha_in, beta_in, limits, retain_sparsity, error)
    INTEGER, INTENT(IN)                      :: mp_group, io_unit
    INTEGER, DIMENSION(:), POINTER           :: nproc
    INTEGER, DIMENSION(:), INTENT(in)        :: matrix_sizes
    INTEGER, DIMENSION(:), POINTER           :: bs_m, bs_n, bs_k
    REAL(real_8), DIMENSION(3), INTENT(in)   :: sparsities
    CHARACTER, DIMENSION(3), INTENT(in)      :: symmetries
    REAL(real_8), DIMENSION(2), INTENT(in)   :: alpha_in, beta_in
    INTEGER, DIMENSION(6), INTENT(in)        :: limits
    LOGICAL, INTENT(in)                      :: retain_sparsity
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: &
      fmt_desc = '(A,3(1X,I6),1X,A,2(1X,I5),1X,A,2(1X,L1))', &
      routineN = 'dbcsr_test_mult', routineP = moduleN//':'//routineN
    CHARACTER, DIMENSION(3), PARAMETER :: trans = (/dbcsr_type_normal, &
      dbcsr_type_transposed, dbcsr_type_conjugate/)
    CHARACTER, DIMENSION(4), PARAMETER :: types = (/dbcsr_type_real_4, &
      dbcsr_type_real_8, dbcsr_type_complex_4, dbcsr_type_complex_8/)

    CHARACTER                                :: transa, transb, TYPE
    INTEGER                                  :: a_c, a_r, a_tr, b_c, b_r, &
                                                b_tr, bmax, bmin, c_c, c_r, &
                                                error_handler, istat, itype, &
                                                mynode, n_tr, numnodes
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: group_sizes
    LOGICAL                                  :: do_complex, pgiven
    TYPE(array_i1d_obj)                      :: sizes_k, sizes_m, sizes_n
    TYPE(dbcsr_data_obj)                     :: data_a, data_b, data_c, &
                                                data_c_dbcsr
    TYPE(dbcsr_obj)                          :: matrix_a, matrix_b, matrix_c
    TYPE(dbcsr_scalar_type)                  :: alpha, beta

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    !
    ! loop over types
    DO itype = 1,SIZE(types)
       TYPE = types(itype)

       do_complex = TYPE.EQ.dbcsr_type_complex_4.OR.TYPE.EQ.dbcsr_type_complex_8

       SELECT CASE (TYPE)
       CASE (dbcsr_type_real_4)
          alpha = dbcsr_scalar (REAL(alpha_in(1), real_4))
          beta  = dbcsr_scalar (REAL( beta_in(1), real_4))
       CASE (dbcsr_type_real_8)
          alpha = dbcsr_scalar (REAL(alpha_in(1), real_8))
          beta  = dbcsr_scalar (REAL( beta_in(1), real_8))
       CASE (dbcsr_type_complex_4)
          alpha = dbcsr_scalar (CMPLX(alpha_in(1), alpha_in(2), real_4))
          beta  = dbcsr_scalar (CMPLX( beta_in(1),  beta_in(2), real_4))
       CASE (dbcsr_type_complex_8)
          alpha = dbcsr_scalar (CMPLX(alpha_in(1), alpha_in(2), real_8))
          beta  = dbcsr_scalar (CMPLX( beta_in(1),  beta_in(2), real_8))
       END SELECT
       !
       ! loop over transpositions
       n_tr = 3
       !IF(do_complex) n_tr = 3
       DO a_tr = 1,n_tr
          DO b_tr = 1,n_tr
             transa = trans(a_tr)
             transb = trans(b_tr)

             ! Create the row/column block sizes.
             IF (ASSOCIATED (bs_m)) THEN
                bmin = MINVAL (bs_m(2::2))
                bmax = MAXVAL (bs_m(2::2))
                CALL dbcsr_make_random_block_sizes (sizes_m, matrix_sizes(1), bs_m)
             ELSE
                CALL dbcsr_make_random_block_sizes (sizes_m, matrix_sizes(1), (/ 1, 13, 2, 5 /))
                bmin = 5 ; bmax = 13
             ENDIF
             IF (ASSOCIATED (bs_n)) THEN
                bmin = MIN (bmin, MINVAL (bs_n(2::2)))
                bmax = MAX (bmax, MAXVAL (bs_n(2::2)))
                CALL dbcsr_make_random_block_sizes (sizes_n, matrix_sizes(2), bs_n)
             ELSE
                CALL dbcsr_make_random_block_sizes (sizes_n, matrix_sizes(2), (/ 1, 13, 2, 5 /))
                bmin = MIN (bmin, 5) ; bmax = MAX (bmax, 13)
             ENDIF
             IF (ASSOCIATED (bs_k)) THEN
                bmin = MIN (bmin, MINVAL (bs_k(2::2)))
                bmax = MAX (bmax, MAXVAL (bs_k(2::2)))
                CALL dbcsr_make_random_block_sizes (sizes_k, matrix_sizes(3), bs_k)
             ELSE
                CALL dbcsr_make_random_block_sizes (sizes_k, matrix_sizes(3), (/ 1, 13, 2, 5 /))
                bmin = MIN (bmin, 5) ; bmax = MAX (bmax, 13)
             ENDIF
             ! Create the undistributed matrices.


             CALL dbcsr_make_random_matrix (matrix_c, sizes_m, sizes_n, "Matrix C",&
                  sparsities(3),&
                  mp_group, data_type=TYPE, symmetry=symmetries(3), error=error)

             IF (transa.EQ.'T') THEN
                CALL dbcsr_make_random_matrix (matrix_a, sizes_k, sizes_m, "Matrix A",&
                     sparsities(1),&
                     mp_group, data_type=TYPE, symmetry=symmetries(1), error=error)
             ELSE
                CALL dbcsr_make_random_matrix (matrix_a, sizes_m, sizes_k, "Matrix A",&
                     sparsities(1),&
                     mp_group, data_type=TYPE, symmetry=symmetries(1), error=error)
             ENDIF
             IF (transb.EQ.'T') THEN
                CALL dbcsr_make_random_matrix (matrix_b, sizes_n, sizes_k, "Matrix B",&
                     sparsities(2),&
                     mp_group, data_type=TYPE, symmetry=symmetries(2), error=error)
             ELSE
                CALL dbcsr_make_random_matrix (matrix_b, sizes_k, sizes_n, "Matrix B",&
                     sparsities(2),&
                     mp_group, data_type=TYPE, symmetry=symmetries(2), error=error)
             ENDIF

             CALL array_release (sizes_m)
             CALL array_release (sizes_n)
             CALL array_release (sizes_k)

             !
             ! convert the dbcsr matrices to denses
             a_r = dbcsr_nfullrows_total(matrix_a); a_c = dbcsr_nfullcols_total(matrix_a)
             b_r = dbcsr_nfullrows_total(matrix_b); b_c = dbcsr_nfullcols_total(matrix_b)
             c_r = dbcsr_nfullrows_total(matrix_c); c_c = dbcsr_nfullcols_total(matrix_c)

             CALL dbcsr_data_init ( data_a )
             CALL dbcsr_data_init ( data_b )
             CALL dbcsr_data_init ( data_c )
             CALL dbcsr_data_init ( data_c_dbcsr )
             CALL dbcsr_data_new ( data_a, dbcsr_type_1d_to_2d(TYPE), data_size=a_r, data_size2=a_c )
             CALL dbcsr_data_new ( data_b, dbcsr_type_1d_to_2d(TYPE), data_size=b_r, data_size2=b_c )
             CALL dbcsr_data_new ( data_c, dbcsr_type_1d_to_2d(TYPE), data_size=c_r, data_size2=c_c )
             CALL dbcsr_data_new ( data_c_dbcsr, dbcsr_type_1d_to_2d(TYPE), data_size=c_r, data_size2=c_c )
             CALL dbcsr_to_dense_local( matrix_a, data_a, limits(1:2), limits(5:6), error=error )
             CALL dbcsr_to_dense_local( matrix_b, data_b, limits(5:6), limits(3:4), error=error )
             CALL dbcsr_to_dense_local( matrix_c, data_c, limits(1:2), limits(3:4), error=error )
             !
             ! Prepare test parameters
             CALL mp_environ (numnodes, mynode, mp_group)
             pgiven = ASSOCIATED (nproc)
             IF (pgiven) pgiven = nproc(1) .NE. 0
             IF (pgiven) THEN
                ALLOCATE (group_sizes (SIZE (nproc), 2))
                group_sizes(:,1) = nproc(:)
                group_sizes(:,2) = 0
             ELSE
                !ALLOCATE (group_sizes (numnodes, 2))
                !DO test = numnodes, 1, -1
                !   group_sizes(1+numnodes-test, 1:2) = (/ test, 0 /)
                !ENDDO
                ALLOCATE (group_sizes (1, 2))
                group_sizes(1, 1:2) = (/ numnodes, 0 /)
             ENDIF
             CALL test_multiplies_multiproc (mp_group, group_sizes,&
                  matrix_a, matrix_b, matrix_c, &
                  data_a, data_b, data_c, data_c_dbcsr, &
                  transa, transb,&
                  alpha, beta,&
                  limits, retain_sparsity, &
                  io_unit = io_unit, error=error)
             !
             ! cleanup
             CALL dbcsr_release (matrix_a)
             CALL dbcsr_release (matrix_b)
             CALL dbcsr_release (matrix_c)

             DEALLOCATE(group_sizes, STAT=istat)
             CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
                  routineN, "deallocattion problem",__LINE__,error)
             CALL dbcsr_data_release( data_a )
             CALL dbcsr_data_release( data_b )
             CALL dbcsr_data_release( data_c )
             CALL dbcsr_data_release( data_c_dbcsr )

          ENDDO
       ENDDO

    ENDDO ! itype

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_test_mult

! *****************************************************************************
!> \brief Performs a variety of matrix multiplies of same matrices on different
!>        processor grids
!> \param[in] mp_group          MPI communicator
!> \param[in] group_sizes       array of (sub) communicator
!>                              sizes to test (2-D)
!> \param[in] matrix_a, matrix_b, matrix_c    matrices to multiply
!> \param[in] io_unit           which unit to write to, if not negative
! *****************************************************************************
  SUBROUTINE test_multiplies_multiproc (mp_group, group_sizes,&
       matrix_a, matrix_b, matrix_c,&
       data_a, data_b, data_c, data_c_dbcsr, &
       transa, transb, alpha, beta, limits, retain_sparsity,&
       io_unit, error)
    INTEGER, INTENT(IN)                      :: mp_group
    INTEGER, DIMENSION(:, :)                 :: group_sizes
    TYPE(dbcsr_obj), INTENT(in)              :: matrix_a, matrix_b, matrix_c
    TYPE(dbcsr_data_obj)                     :: data_a, data_b, data_c, &
                                                data_c_dbcsr
    CHARACTER, INTENT(in)                    :: transa, transb
    TYPE(dbcsr_scalar_type), INTENT(in)      :: alpha, beta
    INTEGER, DIMENSION(6), INTENT(in)        :: limits
    LOGICAL, INTENT(in)                      :: retain_sparsity
    INTEGER, INTENT(IN)                      :: io_unit
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: fmt_per_run_local = &
      '(A,1X,I5,1X,A,1X,F9.3,1X,"s,",1X,EN12.3,1X,"FLOP/s")', &
      fmt_per_run_total = &
      '(A,1X,I5,1X,A,1X,F9.3,1X,"s,",1X,EN12.3,1X,"FLOP/s",1X,E13.5)', &
      routineN = 'test_multiplies_multiproc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, mynode, numnodes, &
                                                pcol, prow, test
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: pgrid
    INTEGER, DIMENSION(2)                    :: myploc, npdims
    LOGICAL                                  :: i_am_alive, success
    TYPE(array_i1d_obj)                      :: col_dist_a, col_dist_b, &
                                                col_dist_c, row_dist_a, &
                                                row_dist_b, row_dist_c
    TYPE(dbcsr_distribution_obj)             :: dist_a, dist_b, dist_c
    TYPE(dbcsr_mp_obj)                       :: mp_env
    TYPE(dbcsr_obj)                          :: m_a, m_b, m_c

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (SIZE(group_sizes, 2), "EQ", 2, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN,&
         "second dimension of group_sizes must be 2",__LINE__,error)
    CALL dbcsr_init (m_a)
    CALL dbcsr_init (m_b)
    CALL dbcsr_init (m_c)
    DO test = 1, SIZE(group_sizes, 1)
       npdims(1:2) = group_sizes(test, 1:2)
       numnodes = npdims(1) * npdims(2)
       CALL dbcsr_assert (numnodes, "GE", 0, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Cartesian sides must be greater or equal to 0",__LINE__,error)
       ! Let MPI decide our process grid
       IF (numnodes .EQ. 0) THEN
          numnodes = MAXVAL (npdims)
          npdims(1:2) = 0
          CALL mp_dims_create (numnodes, npdims)
       ENDIF
       ! Create communicator with given cartesian grid
       CALL mp_cart_create (mp_group, 2, npdims, myploc, group)
       i_am_alive = group .NE. MP_COMM_NULL
       IF (.NOT. i_am_alive) CYCLE
       CALL mp_environ (numnodes, mynode, group)
       ! Create process grid
       IF (ALLOCATED (pgrid)) DEALLOCATE (pgrid)
       ALLOCATE (pgrid(0:npdims(1)-1, 0:npdims(2)-1))
       DO prow = 0, npdims(1)-1
          DO pcol = 0, npdims(2)-1
             CALL mp_cart_rank (group, (/ prow, pcol /), pgrid(prow, pcol))
          ENDDO
       ENDDO
       ! Create the dbcsr_mp_obj
       CALL dbcsr_mp_new (mp_env, pgrid, group, mynode, numnodes,&
            myprow=myploc(1), mypcol=myploc(2))
       ! Row & column distributions
       CALL dbcsr_random_dist (row_dist_a, dbcsr_nblkrows_total (matrix_a), npdims(1))
       CALL dbcsr_random_dist (col_dist_a, dbcsr_nblkcols_total (matrix_a), npdims(2))
       CALL dbcsr_random_dist (row_dist_b, dbcsr_nblkrows_total (matrix_b), npdims(1))
       CALL dbcsr_random_dist (col_dist_b, dbcsr_nblkcols_total (matrix_b), npdims(2))
       CALL dbcsr_random_dist (row_dist_c, dbcsr_nblkrows_total (matrix_c), npdims(1))
       CALL dbcsr_random_dist (col_dist_c, dbcsr_nblkcols_total (matrix_c), npdims(2))
       CALL dbcsr_distribution_new (dist_a, mp_env, row_dist_a, col_dist_a)
       CALL dbcsr_distribution_new (dist_b, mp_env, row_dist_b, col_dist_b)
       CALL dbcsr_distribution_new (dist_c, mp_env, row_dist_c, col_dist_c)
       CALL array_release (row_dist_a)
       CALL array_release (col_dist_a)
       CALL array_release (row_dist_b)
       CALL array_release (col_dist_b)
       CALL array_release (row_dist_c)
       CALL array_release (col_dist_c)
       ! Redistribute the matrices
       ! A
       CALL dbcsr_create (m_a, "Test for "//TRIM(dbcsr_name (matrix_a)),&
            dist_a, dbcsr_get_matrix_type(matrix_a),&
            dbcsr_row_block_sizes (matrix_a),&
            dbcsr_col_block_sizes (matrix_a),&
            data_type=dbcsr_get_data_type (matrix_a),&
            error=error)
       CALL dbcsr_distribution_release (dist_a)
       CALL dbcsr_redistribute (matrix_a, m_a, error=error)
       ! B
       CALL dbcsr_create (m_b, "Test for "//TRIM(dbcsr_name (matrix_b)),&
            dist_b, dbcsr_get_matrix_type(matrix_b),&
            dbcsr_row_block_sizes (matrix_b),&
            dbcsr_col_block_sizes (matrix_b),&
            data_type=dbcsr_get_data_type (matrix_b),&
            error=error)
       CALL dbcsr_distribution_release (dist_b)
       CALL dbcsr_redistribute (matrix_b, m_b, error=error)
       ! C
       CALL dbcsr_create (m_c, "Test for "//TRIM(dbcsr_name (matrix_c)),&
            dist_c, dbcsr_get_matrix_type(matrix_c),&
            dbcsr_row_block_sizes (matrix_c),&
            dbcsr_col_block_sizes (matrix_c),&
            data_type=dbcsr_get_data_type (matrix_c),&
            error=error)
       CALL dbcsr_distribution_release (dist_c)
       CALL dbcsr_redistribute (matrix_c, m_c, error=error)
       !
       ! Perform multiply
       CALL mp_sync (group)
       IF (ALL (limits==0)) THEN
          CALL dbcsr_multiply (transa, transb, alpha,&
               m_a, m_b, beta, m_c,&
               retain_sparsity=retain_sparsity, error=error)
       ELSE
          CALL dbcsr_multiply (transa, transb, alpha,&
               m_a, m_b, beta, m_c,&
               first_row = limits(1),&
               last_row = limits(2),&
               first_column = limits(3),&
               last_column = limits(4),&
               first_k = limits(5),&
               last_k = limits(6),&
               retain_sparsity=retain_sparsity, error=error)
       ENDIF
       CALL dbcsr_replicate_all( m_c, error=error )
       CALL dbcsr_to_dense_local( m_c, data_c_dbcsr, limits(1:2), limits(3:4), error=error)
       CALL dbcsr_check_multiply( data_c_dbcsr, data_a, data_b, data_c, &
            transa, transb, alpha, beta, limits, retain_sparsity, io_unit, &
            success, error=error )
       ! Release
       IF (io_unit .GT. 0) THEN
          IF (success) THEN
             WRITE(io_unit, *) REPEAT("*",70)
             WRITE(io_unit, *) " -- TESTING dbcsr_multiply (",transa,", ",transb,&
                  ", ",dbcsr_get_data_type( m_a ),&
                  ", ",dbcsr_get_matrix_type( m_a ),&
                  ", ",dbcsr_get_matrix_type( m_b ),&
                  ", ",dbcsr_get_matrix_type( m_c ),&
                  ") ............... PASSED !"
             WRITE(io_unit, *) REPEAT("*",70)
          ELSE
             WRITE(io_unit, *) REPEAT("*",70)
             WRITE(io_unit, *) " -- TESTING dbcsr_multiply (",transa,", ",transb,&
                  ", ",dbcsr_get_data_type( m_a ),&
                  ", ",dbcsr_get_matrix_type( m_a ),&
                  ", ",dbcsr_get_matrix_type( m_b ),&
                  ", ",dbcsr_get_matrix_type( m_c ),&
             ") ... FAILED !"
             WRITE(io_unit, *) REPEAT("*",70)
          ENDIF
       ENDIF
       CALL dbcsr_mp_release (mp_env)
       CALL dbcsr_release (m_a)
       CALL dbcsr_release (m_b)
       CALL dbcsr_release (m_c)
    ENDDO
  END SUBROUTINE test_multiplies_multiproc


! *****************************************************************************
!> \brief Performs a check of matrix multiplies 
!> \param[in] dense_c_dbcsr             dense result of the dbcsr_multiply
!> \param[in] dense_a, dense_b, dense_c input dense matrices
!> \param[in] transa, transb            transposition status
!> \param[in] alpha, beta               coefficients for the gemm
!> \param[in] limits                    limits for the gemm
!> \param[in] io_unit                   io unit for printing
!> \param[out] success                  if passed the check success=T
!> \param[inout] error                  dbcsr error 
!>
!> !!!! missing retain_sparsity !!!!
! *****************************************************************************
  SUBROUTINE dbcsr_check_multiply(dense_c_dbcsr, dense_a, dense_b, dense_c,&
       transa, transb, alpha, beta, limits, retain_sparsity, io_unit, &
       success, error)

    TYPE(dbcsr_data_obj), INTENT(inout)      :: dense_c_dbcsr, dense_a, &
                                                dense_b, dense_c
    CHARACTER, INTENT(in)                    :: transa, transb
    TYPE(dbcsr_scalar_type), INTENT(in)      :: alpha, beta
    INTEGER, DIMENSION(6), INTENT(in)        :: limits
    LOGICAL, INTENT(in)                      :: retain_sparsity
    INTEGER, INTENT(IN)                      :: io_unit
    LOGICAL, INTENT(out)                     :: success
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_check_multiply', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: a_col, a_m, a_n, a_row, &
                                                b_col, b_m, b_n, b_row, &
                                                c_col, c_row, i, istat, k, &
                                                lda, ldb, ldc, m, n
    CHARACTER, PARAMETER                     :: norm = 'I'

    LOGICAL                                  :: valid
    REAL(real_4), ALLOCATABLE, DIMENSION(:)  :: work_sp
    REAL(real_4), EXTERNAL                   :: clange, slamch, slange
    REAL(real_8)                             :: a_norm, b_norm, c_norm_dbcsr, &
                                                c_norm_in, c_norm_out, eps, &
                                                residual
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: work
    REAL(real_8), EXTERNAL                   :: dlamch, dlange, zlange

    CALL dbcsr_data_get_sizes (dense_a, lda, i, valid, error=error)
    CALL dbcsr_data_get_sizes (dense_b, ldb, i, valid, error=error)
    CALL dbcsr_data_get_sizes (dense_c, ldc, i, valid, error=error)

    IF(ANY(limits==0)) THEN
       CALL dbcsr_data_get_sizes (dense_a, m, k, valid, error=error)
       CALL dbcsr_data_get_sizes (dense_b, n, i, valid, error=error)
       a_row = 1; a_col = 1
       b_row = 1; b_col = 1
       c_row = 1; c_col = 1
    ELSE
       m = limits(2) - limits(1) + 1
       n = limits(4) - limits(3) + 1
       k = limits(6) - limits(5) + 1
       a_row = limits(1); a_col = limits(5)
       b_row = limits(5); b_col = limits(3)
       c_row = limits(1); c_col = limits(3)
    ENDIF
    !
    !
    IF(transA == 'N') THEN
       a_m = m
       a_n = k
    ELSE
       a_m = k
       a_n = m
    ENDIF
    IF(transB == 'N') THEN
       b_m = k
       b_n = n
    ELSE
       b_m = n
       b_n = k
    ENDIF
    !
    !
    SELECT CASE (dense_a%d%data_type)
    CASE (dbcsr_type_real_8_2d)
       ALLOCATE(work(MAX(K,MAX(M, N))), STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "allocattion problem",__LINE__,error)
       eps          = dlamch('eps')
       a_norm       = dlange(norm, a_m, a_n, dense_a%d%r2_dp(a_row,a_col), lda, work)
       b_norm       = dlange(norm, b_m, b_n, dense_b%d%r2_dp(b_row,b_col), ldb, work)
       c_norm_in    = dlange(norm,   m,   n, dense_c%d%r2_dp(c_row,c_col), ldc, work)
       c_norm_dbcsr = dlange(norm,   m,   n, dense_c_dbcsr%d%r2_dp(c_row,c_col), ldc, work)
       !
       CALL dgemm(transa, transb, m, n, k, alpha%r_dp, dense_a%d%r2_dp(a_row,a_col), lda, &
            dense_b%d%r2_dp(b_row,b_col), ldb, beta%r_dp, dense_c%d%r2_dp(c_row,c_col), ldc)
       !
       c_norm_out   = dlange(norm, m, n, dense_c%d%r2_dp(c_row,c_col), ldc, work)
       !
       ! take the difference dense/sparse
       dense_c%d%r2_dp(c_row:c_row+m-1,c_col:c_col+n-1) = dense_c%d%r2_dp(c_row:c_row+m-1,c_col:c_col+n-1) &
            - dense_c_dbcsr%d%r2_dp(c_row:c_row+m-1,c_col:c_col+n-1)
       !
       ! compute the residual
       residual     = dlange(norm, m, n, dense_c%d%r2_dp(c_row,c_col), ldc, work)
       DEALLOCATE(work, STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "deallocattion problem",__LINE__,error)
    CASE (dbcsr_type_real_4_2d)
       ALLOCATE(work_sp(MAX(K,MAX(M, N))), STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "allocattion problem",__LINE__,error)
       eps          = slamch('eps')
       a_norm       = slange(norm, a_m, a_n, dense_a%d%r2_sp(a_row,a_col), lda, work_sp)
       b_norm       = slange(norm, b_m, b_n, dense_b%d%r2_sp(b_row,b_col), ldb, work_sp)
       c_norm_in    = slange(norm,   m,   n, dense_c%d%r2_sp(c_row,c_col), ldc, work_sp)
       c_norm_dbcsr = slange(norm,   m,   n, dense_c_dbcsr%d%r2_sp(c_row,c_col), ldc, work_sp)
       !
       CALL sgemm(transa, transb, m, n, k, alpha%r_sp, dense_a%d%r2_sp(a_row,a_col), lda, &
            dense_b%d%r2_sp(b_row,b_col), ldb, beta%r_sp, dense_c%d%r2_sp(c_row,c_col), ldc)
       !
       c_norm_out   = slange(norm, m, n, dense_c%d%r2_sp(c_row,c_col), ldc, work_sp)
       !
       ! take the difference dense/sparse
       dense_c%d%r2_sp(c_row:c_row+m-1,c_col:c_col+n-1) = dense_c%d%r2_sp(c_row:c_row+m-1,c_col:c_col+n-1) &
            - dense_c_dbcsr%d%r2_sp(c_row:c_row+m-1,c_col:c_col+n-1)
       !
       ! compute the residual
       residual     = slange(norm, m, n, dense_c%d%r2_sp(c_row,c_col), ldc, work_sp)
       DEALLOCATE(work_sp, STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "deallocattion problem",__LINE__,error)

    CASE (dbcsr_type_complex_8_2d)
       ALLOCATE(work(MAX(K,MAX(M, N))), STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "allocattion problem",__LINE__,error)
       eps          = dlamch('eps')
       a_norm       = zlange(norm, a_m, a_n, dense_a%d%c2_dp(a_row,a_col), lda, work)
       b_norm       = zlange(norm, b_m, b_n, dense_b%d%c2_dp(b_row,b_col), ldb, work)
       c_norm_in    = zlange(norm,   m,   n, dense_c%d%c2_dp(c_row,c_col), ldc, work)
       c_norm_dbcsr = zlange(norm,   m,   n, dense_c_dbcsr%d%c2_dp(c_row,c_col), ldc, work)
       !
       CALL zgemm(transa, transb, m, n, k, alpha%c_dp, dense_a%d%c2_dp(a_row,a_col), lda, &
            dense_b%d%c2_dp(b_row,b_col), ldb, beta%c_dp, dense_c%d%c2_dp(c_row,c_col), ldc)
       !
       c_norm_out   = zlange(norm, m, n, dense_c%d%c2_dp(c_row,c_col), ldc, work)
       !
       ! take the difference dense/sparse
       dense_c%d%c2_dp(c_row:c_row+m-1,c_col:c_col+n-1) = dense_c%d%c2_dp(c_row:c_row+m-1,c_col:c_col+n-1) &
            - dense_c_dbcsr%d%c2_dp(c_row:c_row+m-1,c_col:c_col+n-1)
       !
       ! compute the residual
       residual     = zlange(norm, m, n, dense_c%d%c2_dp(c_row,c_col), ldc, work)
       DEALLOCATE(work, STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "deallocattion problem",__LINE__,error)
    CASE (dbcsr_type_complex_4_2d)
       ALLOCATE(work_sp(MAX(K,MAX(M, N))), STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "allocattion problem",__LINE__,error)
       eps          = slamch('eps')
       a_norm       = clange('f', a_m, a_n, dense_a%d%c2_sp(a_row,a_col), lda, work_sp)
       b_norm       = clange(norm, b_m, b_n, dense_b%d%c2_sp(b_row,b_col), ldb, work_sp)
       c_norm_in    = clange(norm,   m,   n, dense_c%d%c2_sp(c_row,c_col), ldc, work_sp)
       c_norm_dbcsr = clange(norm,   m,   n, dense_c_dbcsr%d%c2_sp(c_row,c_col), ldc, work_sp)
       !
       CALL cgemm(transa, transb, m, n, k, alpha%c_sp, dense_a%d%c2_sp(a_row,a_col), lda, &
            dense_b%d%c2_sp(b_row,b_col), ldb, beta%c_sp, dense_c%d%c2_sp(c_row,c_col), ldc)
       !
       c_norm_out   = clange(norm, m, n, dense_c%d%c2_sp(c_row,c_col), ldc, work_sp)
       !
       ! take the difference dense/sparse
       dense_c%d%c2_sp(c_row:c_row+m-1,c_col:c_col+n-1) = dense_c%d%c2_sp(c_row:c_row+m-1,c_col:c_col+n-1) &
            - dense_c_dbcsr%d%c2_sp(c_row:c_row+m-1,c_col:c_col+n-1)
       !
       ! compute the residual
       residual     = clange(norm, m, n, dense_c%d%c2_sp(c_row,c_col), ldc, work_sp)
       DEALLOCATE(work_sp, STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "deallocattion problem",__LINE__,error)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Incorrect or 1-D data type", __LINE__, error )
    END SELECT
    !
    ! printing
    IF(io_unit .GT. 0) THEN
       !WRITE(io_unit,*) 'residual ',residual,', a_norm ',a_norm,', b_norm ',b_norm, &
       !     ' c_norm_in ',c_norm_in,' c_norm_out ',c_norm_out,' c_norm_dbcsr ',c_norm_dbcsr
       WRITE(io_unit,*) 'Checking the norm of the difference against reference GEMM '
       WRITE(io_unit,*) '-- ||C_dbcsr - C_dense||_oo/((||A||_oo+||B||_oo+||C||_oo).N.eps) =', &
            residual / ((a_norm + b_norm + c_norm_in) * n * eps)
       !
       ! check for nan or inf here
       IF(residual / ((a_norm + b_norm + c_norm_in) * n * eps).GT.10.0_real_8) THEN
          WRITE(io_unit,*) ' The solution is suspicious !'
          success = .FALSE.
       ELSE
          WRITE(io_unit,*) ' The solution is CORRECT !'
          success = .TRUE.
       ENDIF

    ENDIF

  END SUBROUTINE dbcsr_check_multiply

  SUBROUTINE dbcsr_to_dense_local (sparse, dense, row_limits, col_limits, error)
    TYPE(dbcsr_obj), INTENT(in)              :: sparse
    TYPE(dbcsr_data_obj), INTENT(inout)      :: dense
    INTEGER, DIMENSION(2), INTENT(in)        :: row_limits, col_limits
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_to_dense_local', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: symm
    INTEGER                                  :: col, col_offset, col_size, &
                                                dense_col_size, &
                                                dense_row_size, row, &
                                                row_offset, row_size
    LOGICAL                                  :: tr, valid
    TYPE(dbcsr_data_obj)                     :: block
    TYPE(dbcsr_iterator)                     :: iter

    CALL dbcsr_assert (ALL(row_limits==0).AND.ALL(col_limits==0), dbcsr_fatal_level, &
         dbcsr_wrong_args_error, routineN, "limits NYI",__LINE__,error)

    CALL dbcsr_data_get_sizes (dense, dense_row_size, dense_col_size, valid, error)
    CALL dbcsr_assert (valid, dbcsr_fatal_level, dbcsr_internal_error, &
         routineN, "dense matrix not valid",__LINE__,error)

    symm = dbcsr_get_matrix_type(sparse)
    CALL dbcsr_assert (symm==dbcsr_type_normal, dbcsr_fatal_level, dbcsr_internal_error, &
         routineN, "antisymmetric NYI",__LINE__,error)

    CALL dbcsr_data_clear(dense)
    CALL dbcsr_data_init (block)
    CALL dbcsr_data_new (block, dbcsr_type_1d_to_2d(dbcsr_data_get_type(dense)))
    CALL dbcsr_iterator_start( iter, sparse )
    DO WHILE( dbcsr_iterator_blocks_left(iter) )
       CALL dbcsr_iterator_next_block( iter, row, col, block, tr, &
            row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset )
       CALL dbcsr_block_partial_copy(dst=dense,&
            dst_rs=dense_row_size, dst_cs=dense_col_size, dst_tr=.FALSE.,&
            dst_r_lb = row_offset, dst_c_lb = col_offset,&
            src=block, src_rs=row_size, src_cs=col_size, src_tr=tr,&
            src_r_lb=1, src_c_lb=1, nrow=row_size, ncol=col_size)
    END DO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_data_clear_pointer (block)
    CALL dbcsr_data_release (block)

  END SUBROUTINE dbcsr_to_dense_local

END MODULE dbcsr_test_multiply
