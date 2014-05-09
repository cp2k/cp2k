!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Tests for CP2K DBCSR operations
!> \author  Urban Borstnik
!> \date    2010-02-08
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2010-02-08
! *****************************************************************************
MODULE dbcsr_test_methods
  USE array_types,                     ONLY: array_data,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release,&
                                             array_size
  USE dbcsr_blas_operations,           ONLY: dbcsr_lapack_larnv
  USE dbcsr_block_access,              ONLY: dbcsr_put_block
  USE dbcsr_block_operations,          ONLY: dbcsr_block_conjg,&
                                             dbcsr_block_partial_copy,&
                                             dbcsr_block_scale,&
                                             dbcsr_block_transpose,&
                                             dbcsr_data_clear,&
                                             dbcsr_data_set
  USE dbcsr_data_methods,              ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_get_sizes, dbcsr_data_get_type, &
       dbcsr_data_init, dbcsr_data_new, dbcsr_data_release, &
       dbcsr_scalar_negative, dbcsr_scalar_one, dbcsr_type_1d_to_2d, &
       dbcsr_type_2d_to_1d
  USE dbcsr_dist_operations,           ONLY: dbcsr_get_stored_coordinates
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
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
  USE dbcsr_methods,                   ONLY: &
       dbcsr_distribution_hold, dbcsr_distribution_mp, &
       dbcsr_distribution_new, dbcsr_distribution_release, &
       dbcsr_get_matrix_type, dbcsr_init, dbcsr_max_col_size, &
       dbcsr_max_row_size, dbcsr_mp_mynode, dbcsr_mp_new, dbcsr_mp_numnodes, &
       dbcsr_mp_release, dbcsr_nblkcols_total, dbcsr_nblkrows_total, &
       dbcsr_nfullcols_total, dbcsr_nfullrows_total
!  USE dbcsr_transformations, only:dbcsr_put_block
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE dbcsr_types,                     ONLY: &
       dbcsr_data_obj, dbcsr_distribution_obj, dbcsr_iterator, dbcsr_mp_obj, &
       dbcsr_obj, dbcsr_type_antihermitian, dbcsr_type_antisymmetric, &
       dbcsr_type_hermitian, dbcsr_type_no_symmetry, dbcsr_type_real_default, &
       dbcsr_type_symmetric
  USE dbcsr_util,                      ONLY: dbcsr_verify_matrix
  USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_work_create
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             real_8
  USE message_passing,                 ONLY: mp_comm_null,&
                                             mp_environ

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dbcsr_to_dense_local, dbcsr_impose_sparsity
  PUBLIC :: dbcsr_random_dist, dbcsr_make_random_matrix, dbcsr_make_null_mp, &
       dbcsr_make_random_block_sizes, dbcsr_make_null_dist, dbcsr_test_read_args, &
       atoi, atol, ator

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_test_methods'

CONTAINS

  SUBROUTINE dbcsr_test_read_args(narg, args)
    INTEGER, INTENT(out)                     :: narg
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(out)                            :: args

    CHARACTER(len=default_string_length)     :: line
    INTEGER                                  :: istat

    narg = 0
    DO
       READ(5,*,IOSTAT=istat) line
       IF(istat.NE.0) EXIT
       IF(line(1:1).EQ.'#') CYCLE
       narg = narg + 1
       args(narg) = line
       !WRITE(*,*) 'we read <'//TRIM(args(narg))//'>'
    ENDDO

  END SUBROUTINE dbcsr_test_read_args

! *****************************************************************************
!> \brief Impose sparsity on a dense matrix based on a dbcsr
!> \param[in] sparse          sparse matrix
!> \param[inout] dense        dense matrix
!> \param[inout] error        dbcsr error
!>
!> Take into account the symmetry of the sparse matrix.
!> The dense matrix need to be valid.
!> The operation is done localy.
! *****************************************************************************
  SUBROUTINE dbcsr_impose_sparsity(sparse, dense, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: sparse
    TYPE(dbcsr_data_obj), INTENT(inout)      :: dense
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_impose_sparsity', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: symm
    INTEGER :: blk, col, col_offset, col_size, data_type, dense_col_size, &
      dense_row_size, error_handler, row, row_offset, row_size
    LOGICAL                                  :: valid
    TYPE(dbcsr_data_obj)                     :: tmp
    TYPE(dbcsr_iterator)                     :: iter

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_data_get_sizes ( dense, dense_row_size, dense_col_size, valid, error )
    CALL dbcsr_assert (valid, dbcsr_fatal_level, dbcsr_internal_error, &
         routineN, "dense matrix not valid",__LINE__,error)
    data_type = dbcsr_data_get_type( dense )
    symm = dbcsr_get_matrix_type(sparse)

    CALL dbcsr_data_init ( tmp )
    CALL dbcsr_data_new ( tmp, dbcsr_type_1d_to_2d(data_type), data_size=dense_row_size, &
         data_size2=dense_col_size )
    CALL dbcsr_data_set ( dst=tmp, lb=1, data_size=dense_row_size, src=dense, source_lb=1,&
         lb2=1, data_size2=dense_col_size, source_lb2=1 )
    CALL dbcsr_data_clear(dense)

    CALL dbcsr_iterator_start( iter, sparse )
    DO WHILE( dbcsr_iterator_blocks_left(iter) )
       CALL dbcsr_iterator_next_block( iter, row, col, blk, &
            row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset )
       CALL dbcsr_block_partial_copy(&
            dst=dense,&
            dst_rs=dense_row_size, dst_cs=dense_col_size, dst_tr=.FALSE.,&
            dst_r_lb=row_offset, dst_c_lb=col_offset,&
            src=tmp,&
            src_rs=dense_row_size, src_cs=dense_col_size, src_tr=.FALSE.,&
            src_r_lb=row_offset, src_c_lb=col_offset,&
            nrow=row_size, ncol=col_size)
       IF(symm.NE.dbcsr_type_no_symmetry) THEN
          SELECT CASE(symm)
          CASE(dbcsr_type_symmetric)
             CALL dbcsr_block_partial_copy(&
                  dst=dense,&
                  dst_rs=dense_row_size, dst_cs=dense_col_size, dst_tr=.TRUE.,&
                  dst_r_lb=row_offset, dst_c_lb=col_offset,&
                  src=tmp,&
                  src_rs=dense_row_size, src_cs=dense_col_size, src_tr=.FALSE.,&
                  src_r_lb=row_offset, src_c_lb=col_offset,&
                  nrow=row_size, ncol=col_size)
          CASE(dbcsr_type_antisymmetric)
             CALL dbcsr_block_partial_copy(&
                  dst=dense,&
                  dst_rs=dense_row_size, dst_cs=dense_col_size, dst_tr=.TRUE.,&
                  dst_r_lb=row_offset, dst_c_lb=col_offset,&
                  src=tmp,&
                  src_rs=dense_row_size, src_cs=dense_col_size, src_tr=.FALSE.,&
                  src_r_lb=row_offset, src_c_lb=col_offset,&
                  nrow=row_size, ncol=col_size)
             CALL dbcsr_block_scale (dense, dbcsr_scalar_negative(dbcsr_scalar_one(&
                  dbcsr_type_2d_to_1d(data_type))),&
                  row_size=col_size, col_size=row_size, &
                  lb=col_offset, lb2=row_offset, error=error)
          CASE(dbcsr_type_hermitian)
             CALL dbcsr_block_partial_copy(&
                  dst=dense,&
                  dst_rs=dense_row_size, dst_cs=dense_col_size, dst_tr=.TRUE.,&
                  dst_r_lb=row_offset, dst_c_lb=col_offset,&
                  src=tmp,&
                  src_rs=dense_row_size, src_cs=dense_col_size, src_tr=.FALSE.,&
                  src_r_lb=row_offset, src_c_lb=col_offset,&
                  nrow=row_size, ncol=col_size)
             CALL dbcsr_block_conjg (dense, row_size=col_size, col_size=row_size, &
                  lb=col_offset, lb2=row_offset, error=error)
          CASE(dbcsr_type_antihermitian)
             CALL dbcsr_block_partial_copy(&
                  dst=dense,&
                  dst_rs=dense_row_size, dst_cs=dense_col_size, dst_tr=.TRUE.,&
                  dst_r_lb=row_offset, dst_c_lb=col_offset,&
                  src=tmp,&
                  src_rs=dense_row_size, src_cs=dense_col_size, src_tr=.FALSE.,&
                  src_r_lb=row_offset, src_c_lb=col_offset,&
                  nrow=row_size, ncol=col_size)
             CALL dbcsr_block_scale (dense, dbcsr_scalar_negative(dbcsr_scalar_one(&
                  dbcsr_type_2d_to_1d(data_type))),&
                  row_size=col_size, col_size=row_size, &
                  lb=col_offset, lb2=row_offset, error=error)
             CALL dbcsr_block_conjg (dense, row_size=col_size, col_size=row_size, &
                  lb=col_offset, lb2=row_offset, error=error)
          CASE DEFAULT
             CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, &
                  dbcsr_wrong_args_error, routineN, "wrong matrix symmetry",__LINE__,error)
          END SELECT
       ENDIF
    ENDDO
    CALL dbcsr_iterator_stop ( iter )

    CALL dbcsr_data_release ( tmp )

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_impose_sparsity

! *****************************************************************************
!> \brief Convert a sparse matrix to a dense matrix
!> \param[in] sparse          sparse matrix
!> \param[inout] dense        dense matrix
!> \param[inout] error        dbcsr error
!>
!> Take into account the symmetry of the sparse matrix.
!> The dense matrix need to be valid.
!> The operation is done localy.
! *****************************************************************************
  SUBROUTINE dbcsr_to_dense_local (sparse, dense, error)
    TYPE(dbcsr_obj), INTENT(in)              :: sparse
    TYPE(dbcsr_data_obj), INTENT(inout)      :: dense
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_to_dense_local', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: symm
    INTEGER :: col, col_offset, col_size, data_type, dense_col_size, &
      dense_row_size, error_handler, row, row_offset, row_size
    LOGICAL                                  :: tr, valid
    TYPE(dbcsr_data_obj)                     :: block
    TYPE(dbcsr_iterator)                     :: iter

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_data_get_sizes (dense, dense_row_size, dense_col_size, valid, error)
    CALL dbcsr_assert (valid, dbcsr_fatal_level, dbcsr_internal_error, &
         routineN, "dense matrix not valid",__LINE__,error)

    symm = dbcsr_get_matrix_type(sparse)
    data_type = dbcsr_data_get_type(dense)

    CALL dbcsr_data_clear(dense)
    CALL dbcsr_data_init (block)
    CALL dbcsr_data_new (block, dbcsr_type_1d_to_2d(data_type))
    CALL dbcsr_iterator_start( iter, sparse )
    DO WHILE( dbcsr_iterator_blocks_left(iter) )
       CALL dbcsr_iterator_next_block( iter, row, col, block, tr, &
            row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset )
          CALL dbcsr_block_partial_copy(dst=dense,&
               dst_rs=dense_row_size, dst_cs=dense_col_size, dst_tr=.FALSE.,&
               dst_r_lb=row_offset, dst_c_lb=col_offset,&
               src=block, src_rs=row_size, src_cs=col_size, src_tr=tr,&
               src_r_lb=1, src_c_lb=1, nrow=row_size, ncol=col_size)
       IF(symm.NE.dbcsr_type_no_symmetry.AND.row.NE.col) THEN
          SELECT CASE(symm)
          CASE(dbcsr_type_symmetric)
             CALL dbcsr_block_partial_copy(dst=dense,&
                  dst_rs=dense_row_size, dst_cs=dense_col_size, dst_tr=.TRUE.,&
                  dst_r_lb=row_offset, dst_c_lb=col_offset,&
                  src=block, src_rs=row_size, src_cs=col_size, src_tr=tr,&
                  src_r_lb=1, src_c_lb=1, nrow=row_size, ncol=col_size)
          CASE(dbcsr_type_antisymmetric)
             CALL dbcsr_block_partial_copy(dst=dense,&
                  dst_rs=dense_row_size, dst_cs=dense_col_size, dst_tr=.TRUE.,&
                  dst_r_lb=row_offset, dst_c_lb=col_offset,&
                  src=block, src_rs=row_size, src_cs=col_size, src_tr=tr,&
                     src_r_lb=1, src_c_lb=1, nrow=row_size, ncol=col_size)
             CALL dbcsr_block_scale (dense, dbcsr_scalar_negative(dbcsr_scalar_one(&
                  dbcsr_type_2d_to_1d(data_type))),&
                  row_size=col_size, col_size=row_size, &
                  lb=col_offset, lb2=row_offset, error=error)
          CASE(dbcsr_type_hermitian)
             CALL dbcsr_block_partial_copy(dst=dense,&
                  dst_rs=dense_row_size, dst_cs=dense_col_size, dst_tr=.TRUE.,&
                  dst_r_lb=row_offset, dst_c_lb=col_offset,&
                  src=block, src_rs=row_size, src_cs=col_size, src_tr=tr,&
                     src_r_lb=1, src_c_lb=1, nrow=row_size, ncol=col_size)
             CALL dbcsr_block_conjg (dense, row_size=col_size, col_size=row_size, &
                  lb=col_offset, lb2=row_offset, error=error)
          CASE(dbcsr_type_antihermitian)
             CALL dbcsr_block_partial_copy(dst=dense,&
                  dst_rs=dense_row_size, dst_cs=dense_col_size, dst_tr=.TRUE.,&
                  dst_r_lb=row_offset, dst_c_lb=col_offset,&
                  src=block, src_rs=row_size, src_cs=col_size, src_tr=tr,&
                  src_r_lb=1, src_c_lb=1, nrow=row_size, ncol=col_size)
             CALL dbcsr_block_scale (dense, dbcsr_scalar_negative(dbcsr_scalar_one(&
                  dbcsr_type_2d_to_1d(data_type))),&
                  row_size=col_size, col_size=row_size, &
                  lb=col_offset, lb2=row_offset, error=error)
             CALL dbcsr_block_conjg (dense, row_size=col_size, col_size=row_size, &
                  lb=col_offset, lb2=row_offset, error=error)
          CASE DEFAULT
             CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, &
                  dbcsr_wrong_args_error, routineN, "wrong matrix symmetry",__LINE__,error)
          END SELECT
       ENDIF
    END DO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_data_clear_pointer (block)
    CALL dbcsr_data_release (block)

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_to_dense_local

  SUBROUTINE dbcsr_random_dist (dist_array, dist_size, nbins)
    TYPE(array_i1d_obj), INTENT(out)         :: dist_array
    INTEGER, INTENT(in)                      :: dist_size, nbins

    INTEGER                                  :: i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: grid_dist

!REAL, ALLOCATABLE, DIMENSION(:)          :: grid_dist
!

    ALLOCATE (grid_dist(dist_size))
    CALL array_nullify (dist_array)
    !CALL RANDOM_NUMBER (grid_dist)
    FORALL (i = 1 : dist_size)
       grid_dist(i) = MODULO (nbins-i, nbins)
    END FORALL
    !CALL array_new (dist_array, INT(grid_dist*REAL(nbins)), lb=1)
    CALL array_new (dist_array, grid_dist, lb=1)
  END SUBROUTINE dbcsr_random_dist


!> \brief Creates a random matrix.
!> \param[in] nondet     (optional) Specifies that deterministic behavior is
!>                       unneeded.
  SUBROUTINE dbcsr_make_random_matrix (matrix, row_blk_sizes, col_blk_sizes,&
       name, sparsity, mp_group, data_type, symmetry, dist, nondet, error)
    TYPE(dbcsr_obj), INTENT(out)             :: matrix
    TYPE(array_i1d_obj), INTENT(in)          :: row_blk_sizes, col_blk_sizes
    CHARACTER(len=*), INTENT(in)             :: name
    REAL(kind=real_8), INTENT(in)            :: sparsity
    INTEGER, INTENT(in)                      :: mp_group
    INTEGER, INTENT(in), OPTIONAL            :: data_type
    CHARACTER, INTENT(in), OPTIONAL          :: symmetry
    TYPE(dbcsr_distribution_obj), &
      INTENT(IN), OPTIONAL                   :: dist
    LOGICAL, INTENT(in), OPTIONAL            :: nondet
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_random_matrix', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: my_symmetry
    INTEGER                                  :: col, error_handle, max_nze, &
                                                my_data_type, my_proc, &
                                                numproc, nze, p, row, s_col, &
                                                s_row
    INTEGER, DIMENSION(4), SAVE              :: iseed = (/4,3,2,1/)
    INTEGER, DIMENSION(:), POINTER           :: cbs, rbs
    LOGICAL                                  :: quick, tr
    REAL(kind=real_8)                        :: my_sparsity
    REAL(kind=real_8), DIMENSION(1)          :: value
    TYPE(dbcsr_data_obj)                     :: data_values, data_values_tr
    TYPE(dbcsr_distribution_obj)             :: new_dist

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    !
    ! Create the matrix
    IF (PRESENT (dist)) THEN
       new_dist = dist
       CALL dbcsr_distribution_hold (new_dist)
    ELSE
       CALL dbcsr_make_null_dist (new_dist, array_size (row_blk_sizes),&
            array_size (col_blk_sizes), group=mp_group)
    ENDIF
    IF (PRESENT (nondet)) THEN
       quick = nondet
    ELSE
       quick = .FALSE.
    ENDIF
    CALL dbcsr_init (matrix)
    my_data_type = dbcsr_type_real_default
    IF(PRESENT(data_type)) my_data_type = data_type
    my_symmetry = dbcsr_type_no_symmetry
    IF(PRESENT(symmetry)) my_symmetry = symmetry
    CALL dbcsr_create (matrix, name,&
         new_dist, my_symmetry,&
         row_blk_sizes,&
         col_blk_sizes,&
         data_type=my_data_type,&
         error=error)
    rbs => array_data (row_blk_sizes)
    cbs => array_data (col_blk_sizes)
    numproc = dbcsr_mp_numnodes (dbcsr_distribution_mp (new_dist))
    !
    IF (sparsity .GT. 1) THEN
       my_sparsity = sparsity / 100.0
    ELSE
       my_sparsity = sparsity
    ENDIF
    max_nze = dbcsr_max_row_size (matrix) * dbcsr_max_col_size (matrix)
    CALL dbcsr_work_create (matrix,&
         nblks_guess=INT(REAL(dbcsr_nblkrows_total(matrix),KIND=dp)&
                        *REAL(dbcsr_nblkcols_total(matrix),KIND=dp)&
                        *(1.0_dp-sparsity)*1.1_dp/numproc),&
         sizedata_guess=INT(REAL(dbcsr_nfullrows_total(matrix),KIND=dp)&
                           *REAL(dbcsr_nfullcols_total(matrix),KIND=dp)&
                           *(1.0_dp-sparsity)*1.1_dp/numproc),&
         work_mutable = .TRUE.,error=error)

    CALL dbcsr_data_init (data_values)
    CALL dbcsr_data_new (data_values, my_data_type, data_size=max_nze)
    CALL dbcsr_data_init (data_values_tr)
    CALL dbcsr_data_new (data_values_tr, my_data_type, data_size=max_nze)
    my_proc = dbcsr_mp_mynode (dbcsr_distribution_mp (new_dist))
    DO row = 1, dbcsr_nblkrows_total (matrix)
       DO col = 1, dbcsr_nblkcols_total (matrix)
          ! We have to extract random number generation outside so
          ! that results are the same regardless of the process count.
          IF (.NOT. quick) THEN
             CALL dlarnv( 1, iseed, 1, value )
             nze = rbs(row) * cbs(col)
             CALL dbcsr_lapack_larnv( 1, iseed, nze, data_values, error )
          ENDIF
          !
          ! build the upper matrix if some symmetry.
          s_row=row ; s_col=col
          IF (PRESENT (dist)) THEN
             tr = .FALSE.
             CALL dbcsr_get_stored_coordinates(matrix, s_row, s_col, tr, p)
             IF(my_symmetry.NE.dbcsr_type_no_symmetry.AND.s_col.LT.s_row) CYCLE
             IF (p .NE. my_proc) CYCLE
          ELSE
             IF(my_symmetry.NE.dbcsr_type_no_symmetry.AND.s_col.LT.s_row) CYCLE
          ENDIF
          IF (quick) CALL dlarnv( 1, iseed, 1, value )
          IF (.NOT. PRESENT (dist) .AND. my_proc .NE. 0) CYCLE
          IF (value(1) .LT. my_sparsity) CYCLE
          IF (quick) THEN
             nze = rbs(s_row) * cbs(s_col)
             CALL dbcsr_lapack_larnv( 1, iseed, nze, data_values, error )
          ENDIF
          CALL dbcsr_put_block (matrix, s_row, s_col, data_values)
          IF(my_symmetry.NE.dbcsr_type_no_symmetry.AND.s_col.EQ.s_row) THEN
             SELECT CASE(my_symmetry)
             CASE(dbcsr_type_symmetric)
                CALL dbcsr_block_transpose (data_values_tr, data_values, &
                     row_size=rbs(s_row), col_size=cbs(s_col), lb=1, source_lb=1, &
                     error=error)
             CASE(dbcsr_type_antisymmetric)
                CALL dbcsr_block_transpose (data_values_tr, data_values, &
                     row_size=rbs(s_row), col_size=cbs(s_col), lb=1, source_lb=1, &
                     scale=dbcsr_scalar_negative(dbcsr_scalar_one(my_data_type)), &
                     error=error)
             CASE(dbcsr_type_hermitian)
                CALL dbcsr_block_transpose (data_values_tr, data_values, &
                     row_size=rbs(s_row), col_size=cbs(s_col), lb=1, source_lb=1, &
                     error=error)
                CALL dbcsr_block_conjg (data_values_tr, row_size=cbs(s_col), col_size=rbs(s_row), &
                     lb=1, error=error)
             CASE(dbcsr_type_antihermitian)
                CALL dbcsr_block_transpose (data_values_tr, data_values, &
                     row_size=rbs(s_row), col_size=cbs(s_col), lb=1, source_lb=1, &
                     scale=dbcsr_scalar_negative(dbcsr_scalar_one(my_data_type)), &
                     error=error)
                CALL dbcsr_block_conjg (data_values_tr, row_size=cbs(s_col), col_size=rbs(s_row), &
                     lb=1, error=error)
             CASE DEFAULT
                CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, &
                     dbcsr_wrong_args_error, routineN, "wrong matrix symmetry",&
                     __LINE__,error)
             END SELECT
             CALL dbcsr_put_block (matrix, s_row, s_col, data_values_tr, summation=.TRUE.)
          ENDIF
       ENDDO
    ENDDO

    CALL dbcsr_data_release (data_values)
    CALL dbcsr_data_release (data_values_tr)

    CALL dbcsr_distribution_release (new_dist)
    CALL dbcsr_finalize (matrix, error=error)
    CALL dbcsr_verify_matrix (matrix, error=error)
    !
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_make_random_matrix

  SUBROUTINE dbcsr_make_random_block_sizes(block_sizes, size_sum, size_mix)
    TYPE(array_i1d_obj), INTENT(out)         :: block_sizes
    INTEGER, INTENT(in)                      :: size_sum
    INTEGER, DIMENSION(:), INTENT(in)        :: size_mix

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_make_random_block_sizes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: block_size, current_sum, &
                                                istat, nblocks, nsize_mix, &
                                                selector
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: mixer
    INTEGER, DIMENSION(:), POINTER           :: sizes
    TYPE(dbcsr_error_type)                   :: error

!

    NULLIFY (sizes)
    nsize_mix = SIZE (size_mix) / 2
    ALLOCATE (mixer (3, nsize_mix))
    mixer(1, :) = size_mix(1:nsize_mix*2-1:2)
    mixer(2, :) = size_mix(2:nsize_mix*2:2)
    mixer(3, :) = 1
    nblocks = 0
    current_sum = 0
    CALL ensure_array_size (sizes, lb=1, ub=1, error=error)

    selector = 1
    !
    DO WHILE (current_sum .LT. size_sum)
       nblocks = nblocks+1
       !CALL RANDOM_NUMBER(value)
       !block_size = MIN (INT (value(1) * size_max),&
       !                  size_sum - current_sum)
       block_size = MIN (mixer(2, selector),&
            size_sum - current_sum)
       sizes(nblocks) = block_size
       current_sum = current_sum + block_size
       CALL ensure_array_size (sizes, ub=nblocks+1, factor=2.0_dp, error=error)
       mixer(3, selector) = mixer(3, selector) + 1
       IF (mixer(3, selector) .GT. mixer(1, selector)) THEN
          mixer(3, selector) = 1
          selector = MOD (selector, nsize_mix)+1
       ENDIF
    ENDDO
    CALL array_new (block_sizes, sizes(1:nblocks), lb=1)
    current_sum = SUM (array_data (block_sizes))
    CALL dbcsr_assert (current_sum, "EQ", size_sum, dbcsr_fatal_level,&
         dbcsr_internal_error, routineN, "Incorrect block sizes",__LINE__,error)
    DEALLOCATE(mixer, sizes, STAT=istat)
    CALL dbcsr_assert (istat==0, dbcsr_fatal_level, dbcsr_internal_error, &
         routineN, "deallocation error",__LINE__,error)

  END SUBROUTINE dbcsr_make_random_block_sizes


  SUBROUTINE dbcsr_make_null_mp (mp_env, group)
    TYPE(dbcsr_mp_obj), INTENT(out)          :: mp_env
    INTEGER, INTENT(in), OPTIONAL            :: group

    INTEGER                                  :: mynode, numnodes

    IF (PRESENT (group)) THEN
       CALL mp_environ (numnodes, mynode, group)
       CALL dbcsr_mp_new (mp_env, RESHAPE( (/ 1 /), (/1,1/)),&
            group, mynode, numnodes,&
            myprow=0, mypcol=0)
    ELSE
       CALL dbcsr_mp_new (mp_env, RESHAPE( (/ 1 /), (/1,1/)),&
            MP_COMM_NULL, 0, 1,&
            myprow=0, mypcol=0)
    ENDIF
  END SUBROUTINE dbcsr_make_null_mp
  !
  SUBROUTINE dbcsr_make_null_dist (distribution, nblkrows, nblkcols, group)
    TYPE(dbcsr_distribution_obj), &
      INTENT(out)                            :: distribution
    INTEGER, INTENT(in)                      :: nblkrows, nblkcols
    INTEGER, INTENT(in), OPTIONAL            :: group

    INTEGER                                  :: i
    TYPE(array_i1d_obj)                      :: col_dist, row_dist
    TYPE(dbcsr_mp_obj)                       :: mp_env

    CALL dbcsr_make_null_mp (mp_env, group=group)
    CALL array_new (row_dist, (/ (0, i = 1, nblkrows) /), lb=1)
    CALL array_new (col_dist, (/ (0, i = 1, nblkcols) /), lb=1)
    CALL dbcsr_distribution_new (distribution, mp_env,&
         row_dist, col_dist)
    CALL array_release (row_dist)
    CALL array_release (col_dist)
    CALL dbcsr_mp_release (mp_env)
  END SUBROUTINE dbcsr_make_null_dist

  FUNCTION atoi(a)
    CHARACTER(len=*), INTENT(in)             :: a
    INTEGER                                  :: atoi

    READ(a,'(I6)') atoi
  END FUNCTION atoi

  FUNCTION atol(a)
    CHARACTER(len=*), INTENT(in)             :: a
    LOGICAL                                  :: atol

    READ(a,'(L1)') atol
  END FUNCTION atol

  FUNCTION ator(a)
    CHARACTER(len=*), INTENT(in)             :: a
    REAL(real_8)                             :: ator

    READ(a,'(E26.15)') ator
  END FUNCTION ator

END MODULE dbcsr_test_methods
