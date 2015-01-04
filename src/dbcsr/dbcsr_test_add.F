!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Tests for DBCSR add
!> \author  VW
!> \date    2010
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2010
! *****************************************************************************
MODULE dbcsr_test_add
  USE array_types,                     ONLY: array_data,&
                                             array_i1d_obj,&
                                             array_release
  USE dbcsr_data_methods,              ONLY: dbcsr_data_get_sizes,&
                                             dbcsr_data_init,&
                                             dbcsr_data_new,&
                                             dbcsr_data_release,&
                                             dbcsr_scalar,&
                                             dbcsr_type_1d_to_2d
  USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_new,&
                                             dbcsr_distribution_release
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_caller_error,&
                                             dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_internal_error,&
                                             dbcsr_wrong_args_error
  USE dbcsr_io,                        ONLY: dbcsr_print
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_sizes, dbcsr_get_data_type, dbcsr_get_matrix_type, &
       dbcsr_init, dbcsr_name, dbcsr_nblkcols_total, dbcsr_nblkrows_total, &
       dbcsr_nfullcols_total, dbcsr_nfullrows_total, dbcsr_release, &
       dbcsr_row_block_sizes
  USE dbcsr_operations,                ONLY: dbcsr_add,&
                                             dbcsr_get_occupation
  USE dbcsr_test_methods,              ONLY: atoi,&
                                             atol,&
                                             ator,&
                                             dbcsr_impose_sparsity,&
                                             dbcsr_make_random_block_sizes,&
                                             dbcsr_make_random_matrix,&
                                             dbcsr_random_dist,&
                                             dbcsr_to_dense_local
  USE dbcsr_transformations,           ONLY: dbcsr_redistribute,&
                                             dbcsr_replicate_all
  USE dbcsr_types,                     ONLY: &
       dbcsr_conjugate_transpose, dbcsr_data_obj, dbcsr_distribution_obj, &
       dbcsr_mp_obj, dbcsr_no_transpose, dbcsr_obj, dbcsr_scalar_type, &
       dbcsr_transpose, dbcsr_type_antihermitian, dbcsr_type_antisymmetric, &
       dbcsr_type_complex_4, dbcsr_type_complex_4_2d, dbcsr_type_complex_8, &
       dbcsr_type_complex_8_2d, dbcsr_type_hermitian, dbcsr_type_no_symmetry, &
       dbcsr_type_real_4, dbcsr_type_real_4_2d, dbcsr_type_real_8, &
       dbcsr_type_real_8_2d, dbcsr_type_symmetric
  USE dbcsr_work_operations,           ONLY: dbcsr_create
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_environ

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dbcsr_test_adds

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_test_add'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param group ...
!> \param mp_env ...
!> \param npdims ...
!> \param ionode ...
!> \param narg ...
!> \param args ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_test_adds( group, mp_env, npdims, ionode, narg, args, error)

    INTEGER                                  :: group
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, DIMENSION(2), INTENT(in)        :: npdims
    INTEGER                                  :: ionode, narg
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(IN)                             :: args
    TYPE(dbcsr_error_type)                   :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_test_adds', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iblk, istat, limits(4), &
                                                matrix_sizes(2), &
                                                mblk_to_read, nblk_to_read
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bs_m, bs_n
    LOGICAL                                  :: retain_sparsity
    REAL(real_8)                             :: alpha(2), beta(2), &
                                                sparsities(2)

!
! parsing

    CALL dbcsr_assert( narg.GE.16 , dbcsr_fatal_level, dbcsr_wrong_args_error, &
         routineN, "narg not correct", __LINE__, error )

     matrix_sizes(1) = atoi(args( 2))
     matrix_sizes(2) = atoi(args( 3))
     sparsities(1)   = ator(args( 4))
     sparsities(2)   = ator(args( 5))
     alpha(1)        = ator(args( 6))
     alpha(2)        = ator(args( 7))
     beta(1)         = ator(args( 8))
     beta(2)         = ator(args( 9))
     limits(1)       = atoi(args(10))
     limits(2)       = atoi(args(11))
     limits(3)       = atoi(args(12))
     limits(4)       = atoi(args(13))
     retain_sparsity = atol(args(14))
     mblk_to_read    = atoi(args(15))
     nblk_to_read    = atoi(args(16))

     CALL dbcsr_assert( narg.GE.16+2*(mblk_to_read+nblk_to_read) , &
          dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, "narg not correct", &
          __LINE__, error )

     ALLOCATE(bs_m(2*mblk_to_read), bs_n(2*nblk_to_read), STAT=istat)
     CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
          routineN, "allocattion problem",__LINE__,error)

     i = 16
     DO iblk = 1,mblk_to_read
        i = i + 1
        bs_m(2*(iblk-1)+1) = atoi(args(i))
        i = i + 1
        bs_m(2*(iblk-1)+2) = atoi(args(i))
     ENDDO
     DO iblk = 1,nblk_to_read
        i = i + 1
        bs_n(2*(iblk-1)+1) = atoi(args(i))
        i = i + 1
        bs_n(2*(iblk-1)+2) = atoi(args(i))
     ENDDO

     !
     ! do checks here

     !
     ! if the limits are not specified (i.e 0), we set them here
     IF(limits(1).EQ.0) limits(1) = 1
     IF(limits(2).EQ.0) limits(2) = matrix_sizes(1)
     IF(limits(3).EQ.0) limits(3) = 1
     IF(limits(4).EQ.0) limits(4) = matrix_sizes(2)

     !
     ! lets go !
     CALL dbcsr_test_add_low( group, mp_env, npdims, ionode, matrix_sizes, &
          bs_m, bs_n, sparsities, &
          alpha, beta, limits, retain_sparsity, error=error)

     DEALLOCATE(bs_m, bs_n, STAT=istat)
     CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
          routineN, "deallocattion problem",__LINE__,error)

   END SUBROUTINE dbcsr_test_adds


! *****************************************************************************
!> \brief Performs a variety of matrix add
!> \param[in] mp_group          MPI communicator
!> \param mp_env ...
!> \param npdims ...
!> \param[in] io_unit           which unit to write to, if not negative
!> \param[in] matrix_sizes      size of matrices to test
!> \param[in] bs_m block sizes of the 3 dimensions
!> \param[in] bs_n block sizes of the 3 dimensions
!> \param[in] sparsities        sparsities of matrices to create
!> \param[in] alpha_in alpha value to use in multiply
!> \param[in] beta_in alpha value to use in multiply
!> \param limits ...
!> \param retain_sparsity ...
!> \param[in,out] error         cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_test_add_low (mp_group, mp_env, npdims, io_unit, &
       matrix_sizes, bs_m, bs_n, sparsities, &
       alpha_in, beta_in, limits, retain_sparsity, error)
    INTEGER, INTENT(IN)                      :: mp_group
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, DIMENSION(2), INTENT(in)        :: npdims
    INTEGER, INTENT(IN)                      :: io_unit
    INTEGER, DIMENSION(:), INTENT(in)        :: matrix_sizes, bs_m, bs_n
    REAL(real_8), DIMENSION(2), INTENT(in)   :: sparsities, alpha_in, beta_in
    INTEGER, DIMENSION(4), INTENT(in)        :: limits
    LOGICAL, INTENT(in)                      :: retain_sparsity
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: &
      fmt_desc = '(A,3(1X,I6),1X,A,2(1X,I5),1X,A,2(1X,L1))', &
      routineN = 'dbcsr_test_add_low', routineP = moduleN//':'//routineN
    CHARACTER, DIMENSION(2, 5), PARAMETER :: symmetries = RESHAPE((/&
      dbcsr_type_no_symmetry  , dbcsr_type_no_symmetry  ,dbcsr_type_symmetric &
      , dbcsr_type_symmetric    ,dbcsr_type_antisymmetric, &
      dbcsr_type_antisymmetric, dbcsr_type_hermitian    , dbcsr_type_hermitian&
      , dbcsr_type_antihermitian    , dbcsr_type_antihermitian /), (/2,5/) )
    CHARACTER, DIMENSION(3), PARAMETER :: trans = (/dbcsr_no_transpose, &
      dbcsr_transpose, dbcsr_conjugate_transpose/)
    INTEGER, DIMENSION(4), PARAMETER :: types = (/dbcsr_type_real_4, &
      dbcsr_type_real_8, dbcsr_type_complex_4, dbcsr_type_complex_8/)

    CHARACTER                                :: a_symm, b_symm
    INTEGER                                  :: a_c, a_r, b_c, b_r, &
                                                error_handler, isymm, itype, &
                                                mynode, numnodes, TYPE
    INTEGER, DIMENSION(:), POINTER           :: blk_sizes
    LOGICAL                                  :: do_complex
    TYPE(array_i1d_obj)                      :: my_sizes_m, my_sizes_n, &
                                                sizes_m, sizes_n
    TYPE(dbcsr_data_obj)                     :: data_a, data_a_dbcsr, data_b
    TYPE(dbcsr_obj)                          :: matrix_a, matrix_b
    TYPE(dbcsr_scalar_type)                  :: alpha, beta

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    !
    ! print
    CALL mp_environ (numnodes, mynode, mp_group)
    IF (io_unit .GT. 0) THEN
       WRITE(io_unit,*) 'numnodes',numnodes
       WRITE(io_unit,*) 'matrix_sizes',matrix_sizes
       WRITE(io_unit,*) 'sparsities',sparsities
       WRITE(io_unit,*) 'alpha_in',alpha_in
       WRITE(io_unit,*) 'beta_in',beta_in
       WRITE(io_unit,*) 'limits',limits
       WRITE(io_unit,*) 'retain_sparsity',retain_sparsity
       WRITE(io_unit,*) 'bs_m',bs_m
       WRITE(io_unit,*) 'bs_n',bs_n
    ENDIF
    !
    !
    ! loop over symmetry
    DO isymm = 1,SIZE(symmetries,2)
       a_symm = symmetries(1, isymm)
       b_symm = symmetries(2, isymm)

       IF( a_symm.NE.dbcsr_type_no_symmetry .AND. matrix_sizes(1).NE.matrix_sizes(2) ) CYCLE
       IF( b_symm.NE.dbcsr_type_no_symmetry .AND. matrix_sizes(1).NE.matrix_sizes(2) ) CYCLE

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
             IF(a_symm.EQ.dbcsr_type_hermitian.OR.a_symm.EQ.dbcsr_type_antihermitian) &
                  alpha = dbcsr_scalar (CMPLX(alpha_in(1), 0.0_real_8, real_4))
             IF(b_symm.EQ.dbcsr_type_hermitian.OR.b_symm.EQ.dbcsr_type_antihermitian) beta = &
                  dbcsr_scalar (CMPLX(beta_in(1), 0.0_real_8, real_4))
          CASE (dbcsr_type_complex_8)
             alpha = dbcsr_scalar (CMPLX(alpha_in(1), alpha_in(2), real_8))
             beta  = dbcsr_scalar (CMPLX( beta_in(1),  beta_in(2), real_8))
             IF(a_symm.EQ.dbcsr_type_hermitian.OR.a_symm.EQ.dbcsr_type_antihermitian) &
                  alpha = dbcsr_scalar (CMPLX(alpha_in(1), 0.0_real_8, real_8))
             IF(b_symm.EQ.dbcsr_type_hermitian.OR.b_symm.EQ.dbcsr_type_antihermitian) &
                  beta = dbcsr_scalar (CMPLX(beta_in(1), 0.0_real_8, real_8))
          END SELECT

          !
          ! Create the row/column block sizes.
          CALL dbcsr_make_random_block_sizes (sizes_m, matrix_sizes(1), bs_m)
          CALL dbcsr_make_random_block_sizes (sizes_n, matrix_sizes(2), bs_n)


          !
          ! If A/B has symmetry we need the same row/col blocking
          my_sizes_m = sizes_m
          my_sizes_n = sizes_n
          IF(a_symm.NE.dbcsr_type_no_symmetry) THEN
             my_sizes_n = sizes_m
          ENDIF
          IF(b_symm.NE.dbcsr_type_no_symmetry) THEN
             my_sizes_n = sizes_m
          ENDIF

          IF(.FALSE.) THEN
             blk_sizes => array_data(my_sizes_m)
             WRITE(*,*) 'sizes_m',blk_sizes
             WRITE(*,*) 'sum(sizes_m)',SUM(blk_sizes),' matrix_sizes(1)',matrix_sizes(1)
             blk_sizes => array_data(my_sizes_n)
             WRITE(*,*) 'sizes_n',blk_sizes
             WRITE(*,*) 'sum(sizes_n)',SUM(blk_sizes),' matrix_sizes(2)',matrix_sizes(2)
          ENDIF

          !
          ! Create the undistributed matrices.
          CALL dbcsr_make_random_matrix (matrix_a, my_sizes_m, my_sizes_n, "Matrix A",&
               sparsities(1),&
               mp_group, data_type=TYPE, symmetry=a_symm, error=error)

          CALL dbcsr_make_random_matrix (matrix_b, my_sizes_m, my_sizes_n, "Matrix B",&
               sparsities(2),&
               mp_group, data_type=TYPE, symmetry=b_symm, error=error)

          CALL array_release (sizes_m)
          CALL array_release (sizes_n)

          !
          ! convert the dbcsr matrices to denses
          a_r = dbcsr_nfullrows_total(matrix_a); a_c = dbcsr_nfullcols_total(matrix_a)
          b_r = dbcsr_nfullrows_total(matrix_b); b_c = dbcsr_nfullcols_total(matrix_b)
          CALL dbcsr_data_init ( data_a )
          CALL dbcsr_data_init ( data_b )
          CALL dbcsr_data_init ( data_a_dbcsr )
          CALL dbcsr_data_new ( data_a, dbcsr_type_1d_to_2d(TYPE), data_size=a_r, data_size2=a_c )
          CALL dbcsr_data_new ( data_b, dbcsr_type_1d_to_2d(TYPE), data_size=b_r, data_size2=b_c )
          CALL dbcsr_data_new ( data_a_dbcsr, dbcsr_type_1d_to_2d(TYPE), data_size=a_r, data_size2=a_c )
          CALL dbcsr_to_dense_local ( matrix_a, data_a, error=error )
          CALL dbcsr_to_dense_local ( matrix_b, data_b, error=error )

          IF(.FALSE.) THEN
             CALL dbcsr_print(matrix_a, matlab_format=.TRUE., variable_name='a0', error=error)
             CALL dbcsr_print(matrix_b, matlab_format=.TRUE., variable_name='b0', error=error)
          ENDIF

          !
          ! Prepare test parameters
          CALL test_add (mp_group, mp_env, npdims, io_unit,&
               matrix_a, matrix_b, &
               data_a, data_b, data_a_dbcsr, &
               alpha, beta,&
               limits, retain_sparsity, &
               error=error)
          !
          ! cleanup
          CALL dbcsr_release (matrix_a)
          CALL dbcsr_release (matrix_b)
          CALL dbcsr_data_release( data_a )
          CALL dbcsr_data_release( data_b )
          CALL dbcsr_data_release( data_a_dbcsr )

       ENDDO ! itype

    ENDDO !isymm

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_test_add_low

! *****************************************************************************
!> \brief Performs a variety of matrix add
!> \param[in] mp_group          MPI communicator
!> \param mp_env ...
!> \param npdims ...
!> \param[in] io_unit           which unit to write to, if not negative
!> \param[in] matrix_a matrices to add
!> \param[in] matrix_b matrices to add
!> \param data_a ...
!> \param data_b ...
!> \param data_a_dbcsr ...
!> \param alpha ...
!> \param beta ...
!> \param limits ...
!> \param retain_sparsity ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE test_add (mp_group, mp_env, npdims, io_unit,&
       matrix_a, matrix_b, &
       data_a, data_b, data_a_dbcsr, &
       alpha, beta, limits, retain_sparsity,&
       error)
    INTEGER, INTENT(IN)                      :: mp_group
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, DIMENSION(2), INTENT(in)        :: npdims
    INTEGER, INTENT(IN)                      :: io_unit
    TYPE(dbcsr_obj), INTENT(in)              :: matrix_a, matrix_b
    TYPE(dbcsr_data_obj)                     :: data_a, data_b, data_a_dbcsr
    TYPE(dbcsr_scalar_type), INTENT(in)      :: alpha, beta
    INTEGER, DIMENSION(4), INTENT(in)        :: limits
    LOGICAL, INTENT(in)                      :: retain_sparsity
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'test_add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: c, error_handler, r
    LOGICAL                                  :: success
    REAL(real_8)                             :: occ_a_in, occ_a_out, occ_b
    TYPE(array_i1d_obj)                      :: col_dist_a, col_dist_b, &
                                                row_dist_a, row_dist_b
    TYPE(dbcsr_distribution_obj)             :: dist_a, dist_b
    TYPE(dbcsr_obj)                          :: m_a, m_b

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_init (m_a)
    CALL dbcsr_init (m_b)

    ! Row & column distributions
    CALL dbcsr_random_dist (row_dist_a, dbcsr_nblkrows_total (matrix_a), npdims(1))
    CALL dbcsr_random_dist (col_dist_a, dbcsr_nblkcols_total (matrix_a), npdims(2))
    CALL dbcsr_random_dist (row_dist_b, dbcsr_nblkrows_total (matrix_b), npdims(1))
    CALL dbcsr_random_dist (col_dist_b, dbcsr_nblkcols_total (matrix_b), npdims(2))
    CALL dbcsr_distribution_new (dist_a, mp_env, row_dist_a, col_dist_a)
    CALL dbcsr_distribution_new (dist_b, mp_env, row_dist_b, col_dist_b)
    CALL array_release (row_dist_a)
    CALL array_release (col_dist_a)
    CALL array_release (row_dist_b)
    CALL array_release (col_dist_b)

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

    IF(.FALSE.) THEN
       CALL dbcsr_print(m_a, matlab_format=.FALSE., variable_name='a_in_', error=error)
       CALL dbcsr_print(m_b, matlab_format=.FALSE., variable_name='b_', error=error)
    ENDIF

    occ_a_in = dbcsr_get_occupation(m_a)
    occ_b    = dbcsr_get_occupation(m_b)

    !
    ! Perform add
    IF (ALL (limits==0)) THEN
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, " limits shouldnt be 0 ",__LINE__,error)
    ELSE
       CALL  dbcsr_add(m_a, m_b, alpha, beta, error=error)
    ENDIF

    occ_a_out = dbcsr_get_occupation(m_a)

    IF(.FALSE.) THEN
       PRINT*,'retain_sparsity',retain_sparsity,occ_a_in,occ_b,occ_a_out
       CALL dbcsr_print(m_a, matlab_format=.TRUE., variable_name='a_out_', error=error)
    ENDIF

    CALL dbcsr_replicate_all( m_a, error=error )
    CALL dbcsr_to_dense_local( m_a, data_a_dbcsr, error=error)
    CALL dbcsr_check_add( m_a, data_a_dbcsr, data_a, data_b, &
         alpha, beta, limits, retain_sparsity, io_unit, mp_group, &
         success, error=error )

    r = dbcsr_nfullrows_total( m_a )
    c = dbcsr_nfullcols_total( m_a )

    IF (io_unit .GT. 0) THEN
       IF (success) THEN
          WRITE(io_unit, *) REPEAT("*",70)
          WRITE(io_unit, *) " -- TESTING dbcsr_add (",&
                    dbcsr_get_data_type( m_a ),&
               ", ",dbcsr_get_matrix_type( m_a ),&
               ", ",dbcsr_get_matrix_type( m_b ),&
               ") ............................. PASSED !"
          WRITE(io_unit, *) REPEAT("*",70)
       ELSE
          WRITE(io_unit, *) REPEAT("*",70)
          WRITE(io_unit, *) " -- TESTING dbcsr_add (",&
                     dbcsr_get_data_type( m_a ),&
                ", ",dbcsr_get_matrix_type( m_a ),&
                ", ",dbcsr_get_matrix_type( m_b ),&
                ") ................. FAILED !"
          WRITE(io_unit, *) REPEAT("*",70)
          !CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error, &
          !     routineN, 'for the moment, we stop', __LINE__, error)
       ENDIF
    ENDIF

    CALL dbcsr_release (m_a)
    CALL dbcsr_release (m_b)

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE test_add


! *****************************************************************************
!> \brief Performs a check of matrix adds
!> \param matrix_a ...
!> \param dense_a_dbcsr ...
!> \param dense_a input dense matrices
!> \param dense_b input dense matrices
!> \param[in] alpha coefficients for the add
!> \param[in] beta coefficients for the add
!> \param[in] limits                    limits for the add
!> \param retain_sparsity ...
!> \param[in] io_unit                   io unit for printing
!> \param mp_group ...
!> \param[out] success                  if passed the check success=T
!> \param[inout] error                  dbcsr error
!>
! *****************************************************************************
  SUBROUTINE dbcsr_check_add(matrix_a, dense_a_dbcsr, dense_a, dense_b,&
       alpha, beta, limits, retain_sparsity, io_unit, mp_group, &
       success, error)

    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    TYPE(dbcsr_data_obj), INTENT(inout)      :: dense_a_dbcsr, dense_a, &
                                                dense_b
    TYPE(dbcsr_scalar_type), INTENT(in)      :: alpha, beta
    INTEGER, DIMENSION(4), INTENT(in)        :: limits
    LOGICAL, INTENT(in)                      :: retain_sparsity
    INTEGER, INTENT(IN)                      :: io_unit, mp_group
    LOGICAL, INTENT(out)                     :: success
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_check_add', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: col, col_size, error_handler, &
                                                i, istat, j, ld, lwork, m, &
                                                mynode, n, numnodes, row, &
                                                row_size
    CHARACTER, PARAMETER                     :: norm = 'I'

    LOGICAL                                  :: valid
    REAL(real_4), ALLOCATABLE, DIMENSION(:)  :: work_sp
    REAL(real_4), EXTERNAL                   :: clange, slamch, slange
    REAL(real_8)                             :: a_norm_dbcsr, a_norm_in, &
                                                a_norm_out, b_norm, eps, &
                                                residual
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: work
    REAL(real_8), EXTERNAL                   :: dlamch, dlange, zlange

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL mp_environ (numnodes, mynode, mp_group)

    CALL dbcsr_data_get_sizes (dense_a, row_size, col_size, valid, error)
    CALL dbcsr_assert (valid, dbcsr_fatal_level, dbcsr_internal_error, &
         routineN, "dense matrix not valid",__LINE__,error)
    !
    !
    m = limits(2) - limits(1) + 1
    n = limits(4) - limits(3) + 1
    row = limits(1); col = limits(3)

    !
    ! set the size of the work array
    lwork = row_size
    ld = row_size
    !
    !
    SELECT CASE (dense_a%d%data_type)
    CASE (dbcsr_type_real_8_2d)
       ALLOCATE(work(lwork), STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "allocattion problem",__LINE__,error)
       eps          = dlamch('eps')
       a_norm_in    = dlange(norm, row_size, col_size, dense_a%d%r2_dp(1,1), ld, work)
       b_norm       = dlange(norm, row_size, col_size, dense_b%d%r2_dp(1,1), ld, work)
       a_norm_dbcsr = dlange(norm, row_size, col_size, dense_a_dbcsr%d%r2_dp(1,1), ld, work)
       !
       dense_a%d%r2_dp(row:row+m-1,col:col+n-1) = alpha%r_dp * dense_a%d%r2_dp(row:row+m-1,col:col+n-1) + &
                                                   beta%r_dp * dense_b%d%r2_dp(row:row+m-1,col:col+n-1)
       !
       ! impose the sparsity if needed
       IF(retain_sparsity) CALL dbcsr_impose_sparsity(matrix_a, dense_a, error=error)
       !
       a_norm_out   = dlange(norm, row_size, col_size, dense_a%d%r2_dp(1,1), ld, work)
       !
       ! take the difference dense/sparse
       dense_a%d%r2_dp = dense_a%d%r2_dp - dense_a_dbcsr%d%r2_dp
       !
       ! compute the residual
       residual     = dlange(norm, row_size, col_size, dense_a%d%r2_dp(1,1), ld, work)
       DEALLOCATE(work, STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "deallocattion problem",__LINE__,error)
    CASE (dbcsr_type_real_4_2d)
       ALLOCATE(work_sp(lwork), STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "allocattion problem",__LINE__,error)
       eps          = slamch('eps')
       a_norm_in    = slange(norm, row_size, col_size, dense_a%d%r2_sp(1,1), ld, work_sp)
       b_norm       = slange(norm, row_size, col_size, dense_b%d%r2_sp(1,1), ld, work_sp)
       a_norm_dbcsr = slange(norm, row_size, col_size, dense_a_dbcsr%d%r2_sp(1,1), ld, work_sp)
       !
       dense_a%d%r2_sp(row:row+m-1,col:col+n-1) = alpha%r_sp * dense_a%d%r2_sp(row:row+m-1,col:col+n-1) + &
                                                   beta%r_sp * dense_b%d%r2_sp(row:row+m-1,col:col+n-1)
       !
       ! impose the sparsity if needed
       IF(retain_sparsity) CALL dbcsr_impose_sparsity(matrix_a, dense_a, error=error)
       !
       a_norm_out   = slange(norm, row_size, col_size, dense_a%d%r2_sp(1,1), ld, work_sp)
       !
       ! take the difference dense/sparse
       dense_a%d%r2_sp = dense_a%d%r2_sp - dense_a_dbcsr%d%r2_sp
       !
       ! compute the residual
       residual     = REAL(slange(norm, row_size, col_size, dense_a%d%r2_sp(1,1), ld, work_sp), real_8)
       DEALLOCATE(work_sp, STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "deallocattion problem",__LINE__,error)
    CASE (dbcsr_type_complex_8_2d)
       ALLOCATE(work(lwork), STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "allocattion problem",__LINE__,error)
       eps          = dlamch('eps')
       a_norm_in    = zlange(norm, row_size, col_size, dense_a%d%c2_dp(1,1), ld, work)
       b_norm       = zlange(norm, row_size, col_size, dense_b%d%c2_dp(1,1), ld, work)
       a_norm_dbcsr = zlange(norm, row_size, col_size, dense_a_dbcsr%d%c2_dp(1,1), ld, work)
       !
       dense_a%d%c2_dp(row:row+m-1,col:col+n-1) = alpha%c_dp * dense_a%d%c2_dp(row:row+m-1,col:col+n-1) + &
                                                   beta%c_dp * dense_b%d%c2_dp(row:row+m-1,col:col+n-1)
       !
       ! impose the sparsity if needed
       IF(retain_sparsity) CALL dbcsr_impose_sparsity(matrix_a, dense_a, error=error)
       !
       a_norm_out   = zlange(norm, row_size, col_size, dense_a%d%c2_dp(1,1), ld, work)
       !
       ! take the difference dense/sparse
       dense_a%d%c2_dp = dense_a%d%c2_dp - dense_a_dbcsr%d%c2_dp
       !
       ! compute the residual
       residual     = zlange(norm, row_size, col_size, dense_a%d%c2_dp(1,1), ld, work)
       DEALLOCATE(work, STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "deallocattion problem",__LINE__,error)
    CASE (dbcsr_type_complex_4_2d)
       ALLOCATE(work_sp(lwork), STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "allocattion problem",__LINE__,error)
       eps          = REAL(slamch('eps'),real_8)
       a_norm_in    = clange(norm, row_size, col_size, dense_a%d%c2_sp(1,1), ld, work_sp)
       b_norm       = clange(norm, row_size, col_size, dense_b%d%c2_sp(1,1), ld, work_sp)
       a_norm_dbcsr = clange(norm, row_size, col_size, dense_a_dbcsr%d%c2_sp(1,1), ld, work_sp)
       !
       IF(.FALSE.) THEN
       !IF(io_unit .GT. 0) THEN
          DO j=1,SIZE(dense_a%d%c2_sp,2)
          DO i=1,SIZE(dense_a%d%c2_sp,1)
             WRITE(*,'(A,I3,A,I3,A,E15.7,A,E15.7,A)')'a_in(',i,',',j,')=',REAL(dense_a%d%c2_sp(i,j)),'+',&
                  AIMAG(dense_a%d%c2_sp(i,j)),'i;'
          ENDDO
          ENDDO
          DO j=1,SIZE(dense_b%d%c2_sp,2)
          DO i=1,SIZE(dense_b%d%c2_sp,1)
             WRITE(*,'(A,I3,A,I3,A,E15.7,A,E15.7,A)')'b(',i,',',j,')=',REAL(dense_b%d%c2_sp(i,j)),'+',&
                  AIMAG(dense_b%d%c2_sp(i,j)),'i;'
          ENDDO
          ENDDO
       ENDIF

       dense_a%d%c2_sp(row:row+m-1,col:col+n-1) = alpha%c_sp * dense_a%d%c2_sp(row:row+m-1,col:col+n-1) + &
                                                   beta%c_sp * dense_b%d%c2_sp(row:row+m-1,col:col+n-1)
       !
       ! impose the sparsity if needed
       IF(retain_sparsity) CALL dbcsr_impose_sparsity(matrix_a, dense_a, error=error)
       !

       IF(.FALSE.) THEN
       !IF(io_unit .GT. 0) THEN
          DO j=1,SIZE(dense_a%d%c2_sp,2)
          DO i=1,SIZE(dense_a%d%c2_sp,1)
             WRITE(*,'(A,I3,A,I3,A,E15.7,A,E15.7,A)')'a_out(',i,',',j,')=',REAL(dense_a%d%c2_sp(i,j)),'+',&
                  AIMAG(dense_a%d%c2_sp(i,j)),'i;'
          ENDDO
          ENDDO
          DO j=1,SIZE(dense_a_dbcsr%d%c2_sp,2)
          DO i=1,SIZE(dense_a_dbcsr%d%c2_sp,1)
             WRITE(*,'(A,I3,A,I3,A,E15.7,A,E15.7,A)')'a_dbcsr(',i,',',j,')=',REAL(dense_a_dbcsr%d%c2_sp(i,j)),'+',&
                  AIMAG(dense_a_dbcsr%d%c2_sp(i,j)),'i;'
          ENDDO
          ENDDO
       ENDIF

       a_norm_out   = clange(norm, row_size, col_size, dense_a%d%c2_sp(1,1), ld, work_sp)
       !
       ! take the difference dense/sparse
       dense_a%d%c2_sp = dense_a%d%c2_sp - dense_a_dbcsr%d%c2_sp
       !
       ! compute the residual
       residual     = REAL(clange(norm, row_size, col_size, dense_a%d%c2_sp(1,1), ld, work_sp), real_8)
       DEALLOCATE(work_sp, STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "deallocattion problem",__LINE__,error)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Incorrect or 1-D data type", __LINE__, error )
    END SELECT

    IF(mynode.EQ.0) THEN
       IF(residual / ((a_norm_in + b_norm ) * REAL(row_size,real_8) * eps).GT.10.0_real_8) THEN
          success = .FALSE.
       ELSE
          success = .TRUE.
       ENDIF
    ENDIF
    !
    ! syncronize the result...
    CALL mp_bcast(success, 0, mp_group)
    !
    ! printing
    IF(io_unit .GT. 0) THEN
       WRITE(io_unit,'(2(A,E12.5))') ' residual ',residual,', b_norm ',b_norm
       WRITE(io_unit,'(3(A,E12.5))') ' a_norm_in ',a_norm_in,', a_norm_out ',a_norm_out,&
            ', a_norm_dbcsr ',a_norm_dbcsr
       WRITE(io_unit,'(A)') ' Checking the norm of the difference against reference ADD '
       WRITE(io_unit,'(A,E12.5)') ' -- ||A_dbcsr-A_dense||_oo/((||A||_oo+||B||_oo).N.eps)=', &
            residual / ((a_norm_in + b_norm) * n * eps)
       !
       ! check for nan or inf here
       IF(success) THEN
          WRITE(io_unit,'(A)') ' The solution is CORRECT !'
       ELSE
          WRITE(io_unit,'(A)') ' The solution is suspicious !'
       ENDIF

    ENDIF

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_check_add

END MODULE dbcsr_test_add
