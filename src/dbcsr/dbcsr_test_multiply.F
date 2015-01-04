!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE array_types,                     ONLY: array_data,&
                                             array_i1d_obj,&
                                             array_release
  USE dbcsr_data_methods,              ONLY: dbcsr_data_get_sizes,&
                                             dbcsr_data_init,&
                                             dbcsr_data_new,&
                                             dbcsr_data_release,&
                                             dbcsr_scalar,&
                                             dbcsr_scalar_negative,&
                                             dbcsr_scalar_one,&
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
       dbcsr_col_block_offsets, dbcsr_col_block_sizes, dbcsr_get_data_type, &
       dbcsr_get_matrix_type, dbcsr_init, dbcsr_name, dbcsr_nblkcols_total, &
       dbcsr_nblkrows_total, dbcsr_nfullcols_total, dbcsr_nfullrows_total, &
       dbcsr_release, dbcsr_row_block_offsets, dbcsr_row_block_sizes
  USE dbcsr_multiply_api,              ONLY: dbcsr_multiply
  USE dbcsr_operations,                ONLY: dbcsr_copy,&
                                             dbcsr_get_occupation,&
                                             dbcsr_scale
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
       dbcsr_transpose, dbcsr_type_antisymmetric, dbcsr_type_complex_4, &
       dbcsr_type_complex_4_2d, dbcsr_type_complex_8, &
       dbcsr_type_complex_8_2d, dbcsr_type_no_symmetry, dbcsr_type_real_4, &
       dbcsr_type_real_4_2d, dbcsr_type_real_8, dbcsr_type_real_8_2d, &
       dbcsr_type_symmetric
  USE dbcsr_work_operations,           ONLY: dbcsr_create
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_environ

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dbcsr_test_multiplies

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_test_multiply'

  LOGICAL, PARAMETER :: debug_mod = .FALSE.

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param group ...
!> \param mp_env ...
!> \param npdims ...
!> \param io_unit ...
!> \param narg ...
!> \param args ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_test_multiplies (group, mp_env, npdims, io_unit, narg, args, error)

    INTEGER                                  :: group
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, DIMENSION(2), INTENT(in)        :: npdims
    INTEGER                                  :: io_unit, narg
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(IN)                             :: args
    TYPE(dbcsr_error_type)                   :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_test_multiplies', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iblk, istat, kblk_to_read, &
                                                limits(6), matrix_sizes(3), &
                                                mblk_to_read, nblk_to_read
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bs_k, bs_m, bs_n
    LOGICAL                                  :: retain_sparsity
    REAL(real_8)                             :: alpha(2), beta(2), &
                                                sparsities(3)

!
! parsing

    CALL dbcsr_assert (narg.GE.21 , dbcsr_fatal_level, dbcsr_wrong_args_error, &
         routineN, "narg not correct", __LINE__, error)

     matrix_sizes(1) = atoi(args(2))
     matrix_sizes(2) = atoi(args(3))
     matrix_sizes(3) = atoi(args(4))
     sparsities(1)   = ator(args(5))
     sparsities(2)   = ator(args(6))
     sparsities(3)   = ator(args(7))
     alpha(1)        = ator(args(8))
     alpha(2)        = ator(args(9))
     beta(1)         = ator(args(10))
     beta(2)         = ator(args(11))
     limits(1)       = atoi(args(12))
     limits(2)       = atoi(args(13))
     limits(3)       = atoi(args(14))
     limits(4)       = atoi(args(15))
     limits(5)       = atoi(args(16))
     limits(6)       = atoi(args(17))
     retain_sparsity = atol(args(18))
     mblk_to_read    = atoi(args(19))
     nblk_to_read    = atoi(args(20))
     kblk_to_read    = atoi(args(21))

     CALL dbcsr_assert (narg.GE.21+2*(mblk_to_read+nblk_to_read+kblk_to_read) , &
          dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, "narg not correct", &
          __LINE__, error)

     ALLOCATE(bs_m(2*mblk_to_read), bs_n(2*nblk_to_read), bs_k(2*kblk_to_read), STAT=istat)
     CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
          routineN, "allocattion problem", __LINE__, error)

     i = 21
     DO iblk = 1,mblk_to_read
        i = i + 1
        bs_m(2*(iblk-1)+1) = atoi(args(i))
        i = i + 1
        bs_m(2*(iblk-1)+2) = atoi(args(i))
     END DO
     DO iblk = 1,nblk_to_read
        i = i + 1
        bs_n(2*(iblk-1)+1) = atoi(args(i))
        i = i + 1
        bs_n(2*(iblk-1)+2) = atoi(args(i))
     END DO
     DO iblk = 1,kblk_to_read
        i = i + 1
        bs_k(2*(iblk-1)+1) = atoi(args(i))
        i = i + 1
        bs_k(2*(iblk-1)+2) = atoi(args(i))
     END DO

     !
     ! do checks here

     !
     ! if the limits are not specified (i.e 0), we set them here
     IF (limits(1).EQ.0) limits(1) = 1
     IF (limits(2).EQ.0) limits(2) = matrix_sizes(1)
     IF (limits(3).EQ.0) limits(3) = 1
     IF (limits(4).EQ.0) limits(4) = matrix_sizes(2)
     IF (limits(5).EQ.0) limits(5) = 1
     IF (limits(6).EQ.0) limits(6) = matrix_sizes(3)

     !
     ! lets go !
     CALL dbcsr_test_multiply_low (group, mp_env, npdims, io_unit, matrix_sizes, &
          bs_m, bs_n, bs_k, sparsities, &
          alpha, beta, limits, retain_sparsity, error=error)

     DEALLOCATE(bs_m, bs_n, bs_k, STAT=istat)
     CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
          routineN, "deallocattion problem",__LINE__,error)

   END SUBROUTINE dbcsr_test_multiplies


! *****************************************************************************
!> \brief Performs a variety of matrix multiplies of same matrices on different
!>        processor grids
!> \param[in] mp_group          MPI communicator
!> \param mp_env ...
!> \param npdims ...
!> \param[in] io_unit           which unit to write to, if not negative
!> \param[in] matrix_sizes      size of matrices to test
!> \param bs_m block sizes of the 3 dimension
!> \param bs_n block sizes of the 3 dimension
!> \param bs_k block sizes of the 3 dimension
!> \param[in] sparsities        sparsities of matrices to create
!> \param alpha_in alpha value to use in multiply
!> \param beta_in beta value to use in multiply
!> \param limits ...
!> \param retain_sparsity ...
!> \param[in,out] error         cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_test_multiply_low (mp_group, mp_env, npdims, io_unit, &
       matrix_sizes, bs_m, bs_n, bs_k, sparsities, &
       alpha_in, beta_in, limits, retain_sparsity, error)
    INTEGER, INTENT(IN)                      :: mp_group
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, DIMENSION(2), INTENT(in)        :: npdims
    INTEGER, INTENT(IN)                      :: io_unit
    INTEGER, DIMENSION(:), INTENT(in)        :: matrix_sizes, bs_m, bs_n, bs_k
    REAL(real_8), DIMENSION(3), INTENT(in)   :: sparsities
    REAL(real_8), DIMENSION(2), INTENT(in)   :: alpha_in, beta_in
    INTEGER, DIMENSION(6), INTENT(in)        :: limits
    LOGICAL, INTENT(in)                      :: retain_sparsity
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: &
      fmt_desc = '(A,3(1X,I6),1X,A,2(1X,I5),1X,A,2(1X,L1))', &
      routineN = 'dbcsr_test_multiply_low', routineP = moduleN//':'//routineN
    CHARACTER, DIMENSION(3), PARAMETER :: trans = (/dbcsr_no_transpose, &
      dbcsr_transpose, dbcsr_conjugate_transpose/)
    CHARACTER, DIMENSION(3, 12), PARAMETER :: symmetries = RESHAPE((/&
      dbcsr_type_no_symmetry,   dbcsr_type_no_symmetry,   &
      dbcsr_type_no_symmetry, dbcsr_type_symmetric,     dbcsr_type_no_symmetry&
      ,   dbcsr_type_no_symmetry, dbcsr_type_antisymmetric, &
      dbcsr_type_no_symmetry,   dbcsr_type_no_symmetry, dbcsr_type_no_symmetry&
      ,   dbcsr_type_symmetric,     dbcsr_type_no_symmetry, &
      dbcsr_type_symmetric,     dbcsr_type_symmetric,     &
      dbcsr_type_no_symmetry, dbcsr_type_antisymmetric, dbcsr_type_symmetric, &
      dbcsr_type_no_symmetry, dbcsr_type_no_symmetry,   &
      dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, dbcsr_type_symmetric, &
      dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, &
      dbcsr_type_antisymmetric, dbcsr_type_antisymmetric, &
      dbcsr_type_no_symmetry, dbcsr_type_no_symmetry,   dbcsr_type_no_symmetry&
      ,   dbcsr_type_symmetric, dbcsr_type_symmetric,     dbcsr_type_symmetric&
      ,     dbcsr_type_symmetric, dbcsr_type_antisymmetric, &
      dbcsr_type_antisymmetric, dbcsr_type_symmetric /), (/3,12/))
    INTEGER, DIMENSION(4), PARAMETER :: types = (/dbcsr_type_real_4, &
      dbcsr_type_real_8, dbcsr_type_complex_4, dbcsr_type_complex_8/)

    CHARACTER                                :: a_symm, b_symm, c_symm, &
                                                transa, transb
    INTEGER                                  :: a_c, a_r, a_tr, b_c, b_r, &
                                                b_tr, c_c, c_r, &
                                                error_handler, isymm, itype, &
                                                mynode, numnodes, TYPE
    INTEGER, DIMENSION(:), POINTER           :: blk_sizes
    LOGICAL                                  :: do_complex
    TYPE(array_i1d_obj)                      :: my_sizes_k, my_sizes_m, &
                                                my_sizes_n, sizes_k, sizes_m, &
                                                sizes_n
    TYPE(dbcsr_data_obj)                     :: data_a, data_b, data_c, &
                                                data_c_dbcsr
    TYPE(dbcsr_obj)                          :: matrix_a, matrix_b, matrix_c
    TYPE(dbcsr_scalar_type)                  :: alpha, beta

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)

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
       WRITE(io_unit,*) 'bs_k',bs_k
    END IF
    !
    !
    ! loop over symmetry
    DO isymm = 1,SIZE(symmetries,2)
       a_symm = symmetries(1, isymm)
       b_symm = symmetries(2, isymm)
       c_symm = symmetries(3, isymm)

       IF (a_symm.NE.dbcsr_type_no_symmetry .AND. matrix_sizes(1).NE.matrix_sizes(3)) CYCLE
       IF (b_symm.NE.dbcsr_type_no_symmetry .AND. matrix_sizes(2).NE.matrix_sizes(3)) CYCLE
       IF (c_symm.NE.dbcsr_type_no_symmetry .AND. matrix_sizes(1).NE.matrix_sizes(2)) CYCLE

       !
       ! loop over types
       DO itype = 1,SIZE(types)
          TYPE = types(itype)

          do_complex = TYPE.EQ.dbcsr_type_complex_4.OR.TYPE.EQ.dbcsr_type_complex_8

          SELECT CASE (TYPE)
          CASE (dbcsr_type_real_4)
             alpha = dbcsr_scalar(REAL(alpha_in(1), real_4))
             beta  = dbcsr_scalar(REAL(beta_in(1), real_4))
          CASE (dbcsr_type_real_8)
             alpha = dbcsr_scalar(REAL(alpha_in(1), real_8))
             beta  = dbcsr_scalar(REAL(beta_in(1), real_8))
          CASE (dbcsr_type_complex_4)
             alpha = dbcsr_scalar(CMPLX(alpha_in(1), alpha_in(2), real_4))
             beta  = dbcsr_scalar(CMPLX(beta_in(1),  beta_in(2), real_4))
          CASE (dbcsr_type_complex_8)
             alpha = dbcsr_scalar(CMPLX(alpha_in(1), alpha_in(2), real_8))
             beta  = dbcsr_scalar(CMPLX(beta_in(1),  beta_in(2), real_8))
          END SELECT

          !
          ! loop over transpositions
          DO a_tr = 1,SIZE(trans)
          DO b_tr = 1,SIZE(trans)
             transa = trans(a_tr)
             transb = trans(b_tr)

             !
             ! if C has a symmetry, we need special transpositions
             IF (c_symm.NE.dbcsr_type_no_symmetry) THEN
                IF (.NOT.(transa.EQ.dbcsr_no_transpose.AND.transb.EQ.dbcsr_transpose.OR.&
                          transa.EQ.dbcsr_transpose.AND.transb.EQ.dbcsr_no_transpose.OR.&
                          transa.EQ.dbcsr_no_transpose.AND.transb.EQ.dbcsr_conjugate_transpose.AND.&
                          .NOT.do_complex.OR.&
                          transa.EQ.dbcsr_conjugate_transpose.AND.transb.EQ.dbcsr_no_transpose.AND.&
                          .NOT.do_complex)) CYCLE
             END IF
             !
             ! if C has symmetry and special limits
             IF (c_symm.NE.dbcsr_type_no_symmetry) THEN
                IF (limits(1).NE.1.OR.limits(2).NE.matrix_sizes(1).OR.&
                   limits(3).NE.1.OR.limits(4).NE.matrix_sizes(2)) CYCLE
             END IF

             !
             ! Create the row/column block sizes.
             CALL dbcsr_make_random_block_sizes (sizes_m, matrix_sizes(1), bs_m)
             CALL dbcsr_make_random_block_sizes (sizes_n, matrix_sizes(2), bs_n)
             CALL dbcsr_make_random_block_sizes (sizes_k, matrix_sizes(3), bs_k)

             !
             ! if we have symmetry the row and column block sizes hae to match
             IF (c_symm.NE.dbcsr_type_no_symmetry.AND.a_symm.NE.dbcsr_type_no_symmetry.AND.&
                 b_symm.NE.dbcsr_type_no_symmetry) THEN
                my_sizes_m = sizes_m
                my_sizes_n = sizes_m
                my_sizes_k = sizes_m
             ELSE IF ((c_symm.EQ.dbcsr_type_no_symmetry.AND.a_symm.NE.dbcsr_type_no_symmetry.AND.&
                       b_symm.NE.dbcsr_type_no_symmetry).OR.&
                      (c_symm.NE.dbcsr_type_no_symmetry.AND.a_symm.EQ.dbcsr_type_no_symmetry.AND.&
                       b_symm.NE.dbcsr_type_no_symmetry).OR.&
                      (c_symm.NE.dbcsr_type_no_symmetry.AND.a_symm.NE.dbcsr_type_no_symmetry.AND.&
                       b_symm.EQ.dbcsr_type_no_symmetry)) THEN
                my_sizes_m = sizes_m
                my_sizes_n = sizes_m
                my_sizes_k = sizes_m
             ELSE IF (c_symm.EQ.dbcsr_type_no_symmetry.AND.a_symm.EQ.dbcsr_type_no_symmetry.AND.&
                      b_symm.NE.dbcsr_type_no_symmetry) THEN
                my_sizes_m = sizes_m
                my_sizes_n = sizes_n
                my_sizes_k = sizes_n
             ELSE IF (c_symm.EQ.dbcsr_type_no_symmetry.AND.a_symm.NE.dbcsr_type_no_symmetry.AND.&
                      b_symm.EQ.dbcsr_type_no_symmetry) THEN
                my_sizes_m = sizes_m
                my_sizes_n = sizes_n
                my_sizes_k = sizes_m
             ELSE IF (c_symm.NE.dbcsr_type_no_symmetry.AND.a_symm.EQ.dbcsr_type_no_symmetry.AND.&
                      b_symm.EQ.dbcsr_type_no_symmetry) THEN
                my_sizes_m = sizes_m
                my_sizes_n = sizes_m
                my_sizes_k = sizes_k
             ELSE IF (c_symm.EQ.dbcsr_type_no_symmetry.AND.a_symm.EQ.dbcsr_type_no_symmetry.AND.&
                      b_symm.EQ.dbcsr_type_no_symmetry)THEN
                my_sizes_m = sizes_m
                my_sizes_n = sizes_n
                my_sizes_k = sizes_k
             ELSE
                CALL dbcsr_assert(.FALSE., dbcsr_fatal_level, dbcsr_internal_error, &
                     routineN, " something wrong here... ", __LINE__, error)
             END IF

             IF (.FALSE.) THEN
                blk_sizes => array_data(my_sizes_m)
                WRITE(*,*) 'sizes_m',blk_sizes
                WRITE(*,*) 'sum(sizes_m)',SUM(blk_sizes),' matrix_sizes(1)',matrix_sizes(1)
                blk_sizes => array_data(my_sizes_n)
                WRITE(*,*) 'sizes_n',blk_sizes
                WRITE(*,*) 'sum(sizes_n)',SUM(blk_sizes),' matrix_sizes(2)',matrix_sizes(2)
                blk_sizes => array_data(my_sizes_k)
                WRITE(*,*) 'sizes_k',blk_sizes
                WRITE(*,*) 'sum(sizes_k)',SUM(blk_sizes),' matrix_sizes(3)',matrix_sizes(3)
             END IF

             !
             ! Create the undistributed matrices.
             CALL dbcsr_make_random_matrix (matrix_c, my_sizes_m, my_sizes_n, "Matrix C",&
                  sparsities(3),&
                  mp_group, data_type=TYPE, symmetry=c_symm, error=error)

             IF (transa.NE.dbcsr_no_transpose) THEN
                CALL dbcsr_make_random_matrix (matrix_a, my_sizes_k, my_sizes_m, "Matrix A",&
                     sparsities(1),&
                     mp_group, data_type=TYPE, symmetry=a_symm, error=error)
             ELSE
                CALL dbcsr_make_random_matrix (matrix_a, my_sizes_m, my_sizes_k, "Matrix A",&
                     sparsities(1),&
                     mp_group, data_type=TYPE, symmetry=a_symm, error=error)
             END IF
             IF (transb.NE.dbcsr_no_transpose) THEN
                CALL dbcsr_make_random_matrix (matrix_b, my_sizes_n, my_sizes_k, "Matrix B",&
                     sparsities(2),&
                     mp_group, data_type=TYPE, symmetry=b_symm, error=error)
             ELSE
                CALL dbcsr_make_random_matrix (matrix_b, my_sizes_k, my_sizes_n, "Matrix B",&
                     sparsities(2),&
                     mp_group, data_type=TYPE, symmetry=b_symm, error=error)
             END IF

             CALL array_release (sizes_m)
             CALL array_release (sizes_n)
             CALL array_release (sizes_k)

             !
             ! if C has a symmetry, we build it accordingly, i.e. C=A*A and C=A*(-A)
             IF (c_symm.NE.dbcsr_type_no_symmetry) THEN
                CALL dbcsr_copy (matrix_b, matrix_a, error=error)
                !print*, a_symm,b_symm,dbcsr_get_matrix_type(matrix_a),dbcsr_get_matrix_type(matrix_b)
                IF (c_symm.EQ.dbcsr_type_antisymmetric) THEN
                   CALL dbcsr_scale (matrix_b, &
                        alpha_scalar=dbcsr_scalar_negative(&
                                dbcsr_scalar_one(TYPE)), error=error)
                END IF
             END IF

             !
             ! convert the dbcsr matrices to denses
             a_r = dbcsr_nfullrows_total(matrix_a)
             a_c = dbcsr_nfullcols_total(matrix_a)
             b_r = dbcsr_nfullrows_total(matrix_b)
             b_c = dbcsr_nfullcols_total(matrix_b)
             c_r = dbcsr_nfullrows_total(matrix_c)
             c_c = dbcsr_nfullcols_total(matrix_c)
             CALL dbcsr_data_init (data_a)
             CALL dbcsr_data_init (data_b)
             CALL dbcsr_data_init (data_c)
             CALL dbcsr_data_init (data_c_dbcsr)
             CALL dbcsr_data_new (data_a, dbcsr_type_1d_to_2d(TYPE), data_size=a_r, data_size2=a_c)
             CALL dbcsr_data_new (data_b, dbcsr_type_1d_to_2d(TYPE), data_size=b_r, data_size2=b_c)
             CALL dbcsr_data_new (data_c, dbcsr_type_1d_to_2d(TYPE), data_size=c_r, data_size2=c_c)
             CALL dbcsr_data_new (data_c_dbcsr, dbcsr_type_1d_to_2d(TYPE), data_size=c_r, data_size2=c_c)
             CALL dbcsr_to_dense_local (matrix_a, data_a, error=error)
             CALL dbcsr_to_dense_local (matrix_b, data_b, error=error)
             CALL dbcsr_to_dense_local (matrix_c, data_c, error=error)

             !
             ! Prepare test parameters
             CALL test_multiply (mp_group, mp_env, npdims, io_unit,&
                  matrix_a, matrix_b, matrix_c, &
                  data_a, data_b, data_c, data_c_dbcsr, &
                  transa, transb,&
                  alpha, beta,&
                  limits, retain_sparsity, &
                  error=error)
             !
             ! cleanup
             CALL dbcsr_release (matrix_a)
             CALL dbcsr_release (matrix_b)
             CALL dbcsr_release (matrix_c)
             CALL dbcsr_data_release (data_a)
             CALL dbcsr_data_release (data_b)
             CALL dbcsr_data_release (data_c)
             CALL dbcsr_data_release (data_c_dbcsr)

          END DO
          END DO

       END DO ! itype

    END DO !isymm

    CALL dbcsr_error_stop (error_handler, error)

  END SUBROUTINE dbcsr_test_multiply_low

! *****************************************************************************
!> \brief Performs a variety of matrix multiplies of same matrices on different
!>        processor grids
!> \param[in] mp_group          MPI communicator
!> \param mp_env ...
!> \param npdims ...
!> \param[in] io_unit           which unit to write to, if not negative
!> \param[in] matrix_a matrices to multiply
!> \param[in] matrix_b matrices to multiply
!> \param[in] matrix_c matrices to multiply
!> \param data_a ...
!> \param data_b ...
!> \param data_c ...
!> \param data_c_dbcsr ...
!> \param transa ...
!> \param transb ...
!> \param alpha ...
!> \param beta ...
!> \param limits ...
!> \param retain_sparsity ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE test_multiply (mp_group, mp_env, npdims, io_unit,&
       matrix_a, matrix_b, matrix_c,&
       data_a, data_b, data_c, data_c_dbcsr, &
       transa, transb, alpha, beta, limits, retain_sparsity,&
       error)
    INTEGER, INTENT(IN)                      :: mp_group
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, DIMENSION(2), INTENT(in)        :: npdims
    INTEGER, INTENT(IN)                      :: io_unit
    TYPE(dbcsr_obj), INTENT(in)              :: matrix_a, matrix_b, matrix_c
    TYPE(dbcsr_data_obj)                     :: data_a, data_b, data_c, &
                                                data_c_dbcsr
    CHARACTER, INTENT(in)                    :: transa, transb
    TYPE(dbcsr_scalar_type), INTENT(in)      :: alpha, beta
    INTEGER, DIMENSION(6), INTENT(in)        :: limits
    LOGICAL, INTENT(in)                      :: retain_sparsity
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: fmt_per_run_local = &
      '(A,1X,I5,1X,A,1X,F9.3,1X,"s,",1X,EN12.3,1X,"FLOP/s")', &
      fmt_per_run_total = &
      '(A,1X,I5,1X,A,1X,F9.3,1X,"s,",1X,EN12.3,1X,"FLOP/s",1X,E13.5)', &
      routineN = 'test_multiply', routineP = moduleN//':'//routineN

    INTEGER                                  :: c_a, c_b, c_c, error_handler, &
                                                r_a, r_b, r_c
    INTEGER, DIMENSION(:), POINTER           :: blk_offsets
    LOGICAL                                  :: success
    REAL(real_8)                             :: occ_a, occ_b, occ_c_in, &
                                                occ_c_out
    TYPE(array_i1d_obj)                      :: col_dist_a, col_dist_b, &
                                                col_dist_c, row_dist_a, &
                                                row_dist_b, row_dist_c
    TYPE(dbcsr_distribution_obj)             :: dist_a, dist_b, dist_c
    TYPE(dbcsr_obj)                          :: m_a, m_b, m_c

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)

    IF (debug_mod .AND. io_unit .GT. 0) THEN
       WRITE(io_unit, *) REPEAT("*",70)
       WRITE(io_unit, *) " -- TESTING dbcsr_multiply (",transa,", ",transb,&
            ", ",dbcsr_get_data_type(m_a),&
            ", ",dbcsr_get_matrix_type(m_a),&
            ", ",dbcsr_get_matrix_type(m_b),&
            ", ",dbcsr_get_matrix_type(m_c),&
            ") ............... !"
       WRITE(io_unit, *) REPEAT("*",70)
    END IF

    CALL dbcsr_init (m_a)
    CALL dbcsr_init (m_b)
    CALL dbcsr_init (m_c)

    ! Row & column distributions
    CALL dbcsr_random_dist (row_dist_a, dbcsr_nblkrows_total(matrix_a), npdims(1))
    CALL dbcsr_random_dist (col_dist_a, dbcsr_nblkcols_total(matrix_a), npdims(2))
    CALL dbcsr_random_dist (row_dist_b, dbcsr_nblkrows_total(matrix_b), npdims(1))
    CALL dbcsr_random_dist (col_dist_b, dbcsr_nblkcols_total(matrix_b), npdims(2))
    CALL dbcsr_random_dist (row_dist_c, dbcsr_nblkrows_total(matrix_c), npdims(1))
    CALL dbcsr_random_dist (col_dist_c, dbcsr_nblkcols_total(matrix_c), npdims(2))
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
    CALL dbcsr_create (m_a, "Test for "//TRIM(dbcsr_name(matrix_a)),&
         dist_a, dbcsr_get_matrix_type(matrix_a),&
         dbcsr_row_block_sizes(matrix_a),&
         dbcsr_col_block_sizes(matrix_a),&
         data_type=dbcsr_get_data_type(matrix_a),&
         error=error)
    CALL dbcsr_distribution_release (dist_a)
    CALL dbcsr_redistribute (matrix_a, m_a, error=error)
    ! B
    CALL dbcsr_create (m_b, "Test for "//TRIM(dbcsr_name(matrix_b)),&
         dist_b, dbcsr_get_matrix_type(matrix_b),&
         dbcsr_row_block_sizes(matrix_b),&
         dbcsr_col_block_sizes(matrix_b),&
         data_type=dbcsr_get_data_type(matrix_b),&
         error=error)
    CALL dbcsr_distribution_release (dist_b)
    CALL dbcsr_redistribute (matrix_b, m_b, error=error)
    ! C
    CALL dbcsr_create (m_c, "Test for "//TRIM(dbcsr_name(matrix_c)),&
         dist_c, dbcsr_get_matrix_type(matrix_c),&
         dbcsr_row_block_sizes(matrix_c),&
         dbcsr_col_block_sizes(matrix_c),&
         data_type=dbcsr_get_data_type(matrix_c),&
         error=error)
    CALL dbcsr_distribution_release (dist_c)
    CALL dbcsr_redistribute (matrix_c, m_c, error=error)

    IF (.FALSE.) THEN
       blk_offsets => array_data(dbcsr_row_block_offsets(matrix_c))
       WRITE(*,*) 'row_block_offsets(matrix_c)',blk_offsets
       blk_offsets => array_data(dbcsr_col_block_offsets(matrix_c))
       WRITE(*,*) 'col_block_offsets(matrix_c)',blk_offsets
    END IF

    IF (.FALSE.) THEN
       CALL dbcsr_print (m_c, matlab_format=.FALSE., variable_name='c_in_', error=error)
       CALL dbcsr_print (m_a, matlab_format=.FALSE., variable_name='a_', error=error)
       CALL dbcsr_print (m_b, matlab_format=.FALSE., variable_name='b_', error=error)
       CALL dbcsr_print (m_c, matlab_format=.FALSE., variable_name='c_out_', error=error)
    END IF

    occ_a    = dbcsr_get_occupation(m_a)
    occ_b    = dbcsr_get_occupation(m_b)
    occ_c_in = dbcsr_get_occupation(m_c)

    !
    ! Perform multiply
    IF (ALL(limits==0)) THEN
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, " limits shouldnt be 0 ",__LINE__,error)
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
    END IF

    occ_c_out = dbcsr_get_occupation(m_c)

    IF (.FALSE.) THEN
       PRINT*,'retain_sparsity',retain_sparsity,occ_a,occ_b,occ_c_in,occ_c_out
       CALL dbcsr_print (m_a, matlab_format=.TRUE., variable_name='a_', error=error)
       CALL dbcsr_print (m_b, matlab_format=.TRUE., variable_name='b_', error=error)
       CALL dbcsr_print (m_c, matlab_format=.FALSE., variable_name='c_out_', error=error)
    END IF

    CALL dbcsr_replicate_all (m_c, error=error)
    CALL dbcsr_to_dense_local (m_c, data_c_dbcsr, error=error)
    CALL dbcsr_check_multiply (m_c, data_c_dbcsr, data_a, data_b, data_c, &
         transa, transb, alpha, beta, limits, retain_sparsity, io_unit, mp_group, &
         success, error=error)

    r_a = dbcsr_nfullrows_total(m_a)
    c_a = dbcsr_nfullcols_total(m_a)
    r_b = dbcsr_nfullrows_total(m_b)
    c_b = dbcsr_nfullcols_total(m_b)
    r_c = dbcsr_nfullrows_total(m_c)
    c_c = dbcsr_nfullcols_total(m_c)
    IF (io_unit .GT. 0) THEN
       IF (success) THEN
          WRITE(io_unit, *) REPEAT("*",70)
          WRITE(io_unit, *) " -- TESTING dbcsr_multiply (",transa,", ",transb,&
               ", ",dbcsr_get_data_type(m_a),&
               ", ",dbcsr_get_matrix_type(m_a),&
               ", ",dbcsr_get_matrix_type(m_b),&
               ", ",dbcsr_get_matrix_type(m_c),&
               ") ............... PASSED !"
          WRITE(io_unit, *) REPEAT("*",70)
       ELSE
          WRITE(io_unit, *) REPEAT("*",70)
          WRITE(io_unit, *) " -- TESTING dbcsr_multiply (",transa,", ",transb,&
               ", ",dbcsr_get_data_type(m_a),&
               ", ",dbcsr_get_matrix_type(m_a),&
               ", ",dbcsr_get_matrix_type(m_b),&
               ", ",dbcsr_get_matrix_type(m_c),&
               ") ... FAILED !"
          WRITE(io_unit, *) REPEAT("*",70)
          !CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error, &
          !     routineN, 'for the moment, we stop', __LINE__, error)
       END IF
    END IF

    CALL dbcsr_release (m_a)
    CALL dbcsr_release (m_b)
    CALL dbcsr_release (m_c)

    CALL dbcsr_error_stop (error_handler, error)

  END SUBROUTINE test_multiply


! *****************************************************************************
!> \brief Performs a check of matrix multiplies
!> \param matrix_c ...
!> \param[in] dense_c_dbcsr             dense result of the dbcsr_multiply
!> \param dense_a input dense matrices
!> \param dense_b input dense matrices
!> \param dense_c input dense matrices
!> \param transa transposition status
!> \param transb transposition status
!> \param alpha coefficients for the gemm
!> \param beta coefficients for the gemm
!> \param[in] limits                    limits for the gemm
!> \param retain_sparsity ...
!> \param[in] io_unit                   io unit for printing
!> \param mp_group ...
!> \param[out] success                  if passed the check success=T
!> \param[inout] error                  dbcsr error
!>
! *****************************************************************************
  SUBROUTINE dbcsr_check_multiply (matrix_c, dense_c_dbcsr, dense_a, dense_b, dense_c,&
       transa, transb, alpha, beta, limits, retain_sparsity, io_unit, mp_group, &
       success, error)

    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_c
    TYPE(dbcsr_data_obj), INTENT(inout)      :: dense_c_dbcsr, dense_a, &
                                                dense_b, dense_c
    CHARACTER, INTENT(in)                    :: transa, transb
    TYPE(dbcsr_scalar_type), INTENT(in)      :: alpha, beta
    INTEGER, DIMENSION(6), INTENT(in)        :: limits
    LOGICAL, INTENT(in)                      :: retain_sparsity
    INTEGER, INTENT(IN)                      :: io_unit, mp_group
    LOGICAL, INTENT(out)                     :: success
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_check_multiply', &
      routineP = moduleN//':'//routineN
    INTEGER :: a_col, a_m, a_n, a_row, b_col, b_m, b_n, b_row, c_col, &
      c_col_size, c_row, c_row_size, error_handler, i, istat, j, k, lda, ldb, &
      ldc, lwork, m, mynode, n, numnodes
    CHARACTER, PARAMETER                     :: norm = 'I'

    LOGICAL                                  :: valid
    REAL(real_4), ALLOCATABLE, DIMENSION(:)  :: work_sp
    REAL(real_4), EXTERNAL                   :: clange, slamch, slange
    REAL(real_8)                             :: a_norm, b_norm, c_norm_dbcsr, &
                                                c_norm_in, c_norm_out, eps, &
                                                residual
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: work
    REAL(real_8), EXTERNAL                   :: dlamch, dlange, zlange

    CALL dbcsr_error_set (routineN, error_handler, error)

    CALL mp_environ (numnodes, mynode, mp_group)

    CALL dbcsr_data_get_sizes (dense_c, c_row_size, c_col_size, valid, error)
    CALL dbcsr_assert (valid, dbcsr_fatal_level, dbcsr_internal_error, &
         routineN, "dense matrix not valid",__LINE__,error)
    CALL dbcsr_data_get_sizes (dense_c, ldc, i, valid, error=error)
    CALL dbcsr_assert (valid, dbcsr_fatal_level, dbcsr_internal_error, &
         routineN, "dense matrix not valid",__LINE__,error)
    CALL dbcsr_data_get_sizes (dense_a, lda, i, valid, error=error)
    CALL dbcsr_assert (valid, dbcsr_fatal_level, dbcsr_internal_error, &
         routineN, "dense matrix not valid",__LINE__,error)
    CALL dbcsr_data_get_sizes (dense_b, ldb, i, valid, error=error)
    CALL dbcsr_assert (valid, dbcsr_fatal_level, dbcsr_internal_error, &
         routineN, "dense matrix not valid",__LINE__,error)
    !
    !
    m = limits(2) - limits(1) + 1
    n = limits(4) - limits(3) + 1
    k = limits(6) - limits(5) + 1
    a_row = limits(1); a_col = limits(5)
    b_row = limits(5); b_col = limits(3)
    c_row = limits(1); c_col = limits(3)
    !
    !
    IF (transA == dbcsr_no_transpose) THEN
       a_m = m
       a_n = k
    ELSE
       a_m = k
       a_n = m
       i = a_row
       a_row = a_col
       a_col = i
    END IF
    IF (transB == dbcsr_no_transpose) THEN
       b_m = k
       b_n = n
    ELSE
       b_m = n
       b_n = k
       i = b_row
       b_row = b_col
       b_col = i
    END IF
    !
    ! set the size of the work array
    lwork = MAXVAL((/lda, ldb, ldc/))
    !
    !
    SELECT CASE (dense_a%d%data_type)
    CASE (dbcsr_type_real_8_2d)
       ALLOCATE(work(lwork), STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "allocation problem",__LINE__,error)
       eps          = dlamch('eps')
       a_norm       = dlange(norm, a_m, a_n, dense_a%d%r2_dp(a_row,a_col), lda, work)
       b_norm       = dlange(norm, b_m, b_n, dense_b%d%r2_dp(b_row,b_col), ldb, work)
       c_norm_in    = dlange(norm, c_row_size, c_col_size, dense_c%d%r2_dp(1,1), ldc, work)
       c_norm_dbcsr = dlange(norm, c_row_size, c_col_size, dense_c_dbcsr%d%r2_dp(1,1), ldc, work)
       !
       CALL dgemm (transa, transb, m, n, k, alpha%r_dp, dense_a%d%r2_dp(a_row,a_col), lda, &
            dense_b%d%r2_dp(b_row,b_col), ldb, beta%r_dp, dense_c%d%r2_dp(c_row,c_col), ldc)
       !
       ! impose the sparsity if needed
       IF (retain_sparsity) CALL dbcsr_impose_sparsity (matrix_c, dense_c, error=error)
       !
       c_norm_out   = dlange(norm, m, n, dense_c%d%r2_dp(c_row,c_col), ldc, work)
       !
       ! take the difference dense/sparse
       dense_c%d%r2_dp = dense_c%d%r2_dp - dense_c_dbcsr%d%r2_dp
       !
       ! compute the residual
       residual     = dlange(norm, c_row_size, c_col_size, dense_c%d%r2_dp(1,1), ldc, work)
       DEALLOCATE(work, STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "deallocation problem",__LINE__,error)
    CASE (dbcsr_type_real_4_2d)
       ALLOCATE(work_sp(lwork), STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "allocation problem",__LINE__,error)
       eps          = REAL(slamch('eps'),real_8)
       a_norm       = slange(norm, a_m, a_n, dense_a%d%r2_sp(a_row,a_col), lda, work_sp)
       b_norm       = slange(norm, b_m, b_n, dense_b%d%r2_sp(b_row,b_col), ldb, work_sp)
       c_norm_in    = slange(norm, c_row_size, c_col_size, dense_c%d%r2_sp(1,1), ldc, work_sp)
       c_norm_dbcsr = slange(norm, c_row_size, c_col_size, dense_c_dbcsr%d%r2_sp(1,1), ldc, work_sp)
       !

       IF (.FALSE.) THEN
       !IF (io_unit .GT. 0) THEN
          DO j=1,SIZE(dense_a%d%r2_sp,2)
             DO i=1,SIZE(dense_a%d%r2_sp,1)
                WRITE(*,'(A,I3,A,I3,A,E15.7,A)')'a(',i,',',j,')=',dense_a%d%r2_sp(i,j),';'
             END DO
          END DO
          DO j=1,SIZE(dense_b%d%r2_sp,2)
             DO i=1,SIZE(dense_b%d%r2_sp,1)
                WRITE(*,'(A,I3,A,I3,A,E15.7,A)')'b(',i,',',j,')=',dense_b%d%r2_sp(i,j),';'
             END DO
          END DO
          DO j=1,SIZE(dense_c%d%r2_sp,2)
             DO i=1,SIZE(dense_c%d%r2_sp,1)
                WRITE(*,'(A,I3,A,I3,A,E15.7,A)')'c_in(',i,',',j,')=',dense_c%d%r2_sp(i,j),';'
             END DO
          END DO
       END IF

       CALL sgemm (transa, transb, m, n, k, alpha%r_sp, dense_a%d%r2_sp(a_row,a_col), lda, &
            dense_b%d%r2_sp(b_row,b_col), ldb, beta%r_sp, dense_c%d%r2_sp(c_row,c_col), ldc)
       !
       ! impose the sparsity if needed
       IF (retain_sparsity) CALL dbcsr_impose_sparsity (matrix_c, dense_c, error=error)

       IF (.FALSE.) THEN
       !IF (io_unit .GT. 0) THEN
          DO j=1,SIZE(dense_c%d%r2_sp,2)
             DO i=1,SIZE(dense_c%d%r2_sp,1)
                WRITE(*,'(A,I3,A,I3,A,E15.7,A)')'c_out(',i,',',j,')=',dense_c%d%r2_sp(i,j),';'
             END DO
          END DO
          DO j=1,SIZE(dense_c_dbcsr%d%r2_sp,2)
             DO i=1,SIZE(dense_c_dbcsr%d%r2_sp,1)
                WRITE(*,'(A,I3,A,I3,A,E15.7,A)')'c_dbcsr(',i,',',j,')=',dense_c_dbcsr%d%r2_sp(i,j),';'
             END DO
          END DO
       END IF
       !
       c_norm_out   = slange(norm, m, n, dense_c%d%r2_sp(c_row,c_col), ldc, work_sp)
       !
       ! take the difference dense/sparse
       dense_c%d%r2_sp = dense_c%d%r2_sp - dense_c_dbcsr%d%r2_sp
       !
       ! compute the residual
       residual     = REAL(slange(norm, c_row_size, c_col_size, dense_c%d%r2_sp(1,1), ldc, work_sp),real_8)
       DEALLOCATE(work_sp, STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "deallocation problem",__LINE__,error)
    CASE (dbcsr_type_complex_8_2d)
       ALLOCATE(work(lwork), STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "allocation problem",__LINE__,error)
       eps          = dlamch('eps')
       a_norm       = zlange(norm, a_m, a_n, dense_a%d%c2_dp(a_row,a_col), lda, work)
       b_norm       = zlange(norm, b_m, b_n, dense_b%d%c2_dp(b_row,b_col), ldb, work)
       c_norm_in    = zlange(norm, c_row_size, c_col_size, dense_c%d%c2_dp(1,1), ldc, work)
       c_norm_dbcsr = zlange(norm, c_row_size, c_col_size, dense_c_dbcsr%d%c2_dp(1,1), ldc, work)
       !
       CALL zgemm (transa, transb, m, n, k, alpha%c_dp, dense_a%d%c2_dp(a_row,a_col), lda, &
            dense_b%d%c2_dp(b_row,b_col), ldb, beta%c_dp, dense_c%d%c2_dp(c_row,c_col), ldc)
       !
       ! impose the sparsity if needed
       IF (retain_sparsity) CALL dbcsr_impose_sparsity (matrix_c, dense_c, error=error)
       !
       c_norm_out   = zlange(norm, m, n, dense_c%d%c2_dp(c_row,c_col), ldc, work)
       !
       ! take the difference dense/sparse
       dense_c%d%c2_dp = dense_c%d%c2_dp - dense_c_dbcsr%d%c2_dp
       !
       ! compute the residual
       residual     = zlange(norm, c_row_size, c_col_size, dense_c%d%c2_dp(1,1), ldc, work)
       DEALLOCATE(work, STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "deallocation problem",__LINE__,error)
    CASE (dbcsr_type_complex_4_2d)
       ALLOCATE(work_sp(lwork), STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "allocation problem",__LINE__,error)
       eps          = REAL(slamch('eps'),real_8)
       a_norm       = clange(norm, a_m, a_n, dense_a%d%c2_sp(a_row,a_col), lda, work_sp)
       b_norm       = clange(norm, b_m, b_n, dense_b%d%c2_sp(b_row,b_col), ldb, work_sp)
       c_norm_in    = clange(norm, c_row_size, c_col_size, dense_c%d%c2_sp(1,1), ldc, work_sp)
       c_norm_dbcsr = clange(norm, c_row_size, c_col_size, dense_c_dbcsr%d%c2_sp(1,1), ldc, work_sp)
       !
       CALL cgemm (transa, transb, m, n, k, alpha%c_sp, dense_a%d%c2_sp(a_row,a_col), lda, &
            dense_b%d%c2_sp(b_row,b_col), ldb, beta%c_sp, dense_c%d%c2_sp(c_row,c_col), ldc)
       !
       ! impose the sparsity if needed
       IF (retain_sparsity) CALL dbcsr_impose_sparsity (matrix_c, dense_c, error=error)
       !
       c_norm_out   = clange(norm, m, n, dense_c%d%c2_sp(c_row,c_col), ldc, work_sp)
       !
       ! take the difference dense/sparse
       dense_c%d%c2_sp = dense_c%d%c2_sp - dense_c_dbcsr%d%c2_sp
       !
       ! compute the residual
       residual     = clange(norm, c_row_size, c_col_size, dense_c%d%c2_sp(1,1), ldc, work_sp)
       DEALLOCATE(work_sp, STAT=istat)
       CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "deallocation problem",__LINE__,error)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Incorrect or 1-D data type", __LINE__, error)
    END SELECT

    IF (mynode.EQ.0) THEN
       IF (residual / ((a_norm + b_norm + c_norm_in) * REAL(n,real_8) * eps).GT.10.0_real_8) THEN
          success = .FALSE.
       ELSE
          success = .TRUE.
       END IF
    END IF
    !
    ! syncronize the result...
    CALL mp_bcast (success, 0, mp_group)
    !
    ! printing
    IF (io_unit .GT. 0) THEN
       WRITE(io_unit,'(3(A,E12.5))') ' residual ',residual,', a_norm ',a_norm,', b_norm ',b_norm
       WRITE(io_unit,'(3(A,E12.5))') ' c_norm_in ',c_norm_in,', c_norm_out ',c_norm_out,&
            ', c_norm_dbcsr ',c_norm_dbcsr
       WRITE(io_unit,'(A)') ' Checking the norm of the difference against reference GEMM '
       WRITE(io_unit,'(A,E12.5)') ' -- ||C_dbcsr-C_dense||_oo/((||A||_oo+||B||_oo+||C||_oo).N.eps)=', &
            residual / ((a_norm + b_norm + c_norm_in) * n * eps)
       !
       ! check for nan or inf here
       IF (success) THEN
          WRITE(io_unit,'(A)') ' The solution is CORRECT !'
       ELSE
          WRITE(io_unit,'(A)') ' The solution is suspicious !'
       END IF

    END IF

    CALL dbcsr_error_stop (error_handler, error)

  END SUBROUTINE dbcsr_check_multiply

END MODULE dbcsr_test_multiply
