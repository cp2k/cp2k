!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
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
MODULE dbcsr_tests
  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_release
  USE dbcsr_block_access
  USE dbcsr_data_methods
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_wrong_args_error
  USE dbcsr_message_passing,           ONLY: dmp_max,&
                                             mp_cart_create,&
                                             mp_cart_rank,&
                                             mp_comm_null,&
                                             mp_dims_create,&
                                             mp_environ,&
                                             mp_sum,&
                                             mp_sync
  USE dbcsr_methods
  USE dbcsr_operations
  USE dbcsr_ptr_util
  USE dbcsr_test_methods
  USE dbcsr_transformations
  USE dbcsr_types
  USE dbcsr_util
  USE dbcsr_work_operations
  USE kinds,                           ONLY: dp,&
                                             int_8,&
                                             real_8
  USE machine,                         ONLY: m_walltime

  !$ USE OMP_LIB
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cp_test_multiplies

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tests'

CONTAINS

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
  SUBROUTINE cp_test_multiplies (mp_group, io_unit, nproc,&
       matrix_sizes, matrix_types, trs, &
       bs_m, bs_n, bs_k, sparsities, alpha, beta, error)
    INTEGER, INTENT(IN)                      :: mp_group, io_unit
    INTEGER, DIMENSION(:), POINTER           :: nproc
    INTEGER, DIMENSION(:), INTENT(in)        :: matrix_sizes
    CHARACTER, DIMENSION(3), INTENT(in)      :: matrix_types
    LOGICAL, DIMENSION(2), INTENT(in)        :: trs
    INTEGER, DIMENSION(:), POINTER           :: bs_m, bs_n, bs_k
    REAL(kind=dp), DIMENSION(3), INTENT(in)  :: sparsities
    REAL(kind=dp), INTENT(in)                :: alpha, beta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      fmt_desc = '(A,3(1X,I6),1X,A,2(1X,I5),1X,A,2(1X,L1))', &
      routineN = 'cp_test_multiplies', routineP = moduleN//':'//routineN

    CHARACTER                                :: t_a, t_b
    INTEGER                                  :: bmax, bmin, mynode, numnodes
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: group_sizes
    LOGICAL                                  :: pgiven
    TYPE(array_i1d_obj)                      :: sizes_k, sizes_m, sizes_n
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_obj)                          :: matrix_a, matrix_b, matrix_c

!   ---------------------------------------------------------------------------
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
         REAL(sparsities(3), real_8),&
         mp_group, error=dbcsr_error)
    IF (trs(1)) THEN
       CALL dbcsr_make_random_matrix (matrix_a, sizes_k, sizes_m, "Matrix A",&
            REAL(sparsities(1), real_8),&
            mp_group, error=dbcsr_error)
    ELSE
       CALL dbcsr_make_random_matrix (matrix_a, sizes_m, sizes_k, "Matrix A",&
            REAL(sparsities(1), real_8),&
            mp_group, error=dbcsr_error)
    ENDIF
    IF (trs(2)) THEN
       CALL dbcsr_make_random_matrix (matrix_b, sizes_n, sizes_k, "Matrix B",&
            REAL(sparsities(2), real_8),&
            mp_group, error=dbcsr_error)
    ELSE
       CALL dbcsr_make_random_matrix (matrix_b, sizes_k, sizes_n, "Matrix B",&
            REAL(sparsities(2), real_8),&
            mp_group, error=dbcsr_error)
    ENDIF
    CALL array_release (sizes_m)
    CALL array_release (sizes_n)
    CALL array_release (sizes_k)
    ! Prepare test parameters
    IF (io_unit .GT. 0) THEN
       WRITE(io_unit, fmt_desc)"Multiplication with sizes",matrix_sizes(1:3),&
            "min/max block sizes", bmin, bmax, "transposed?", trs(1:2)
    ENDIF
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
    t_a = 'N' ; IF (trs(1)) t_a = 'T'
    t_b = 'N' ; IF (trs(2)) t_b = 'T'
    CALL test_multiplies_multiproc (mp_group, group_sizes,&
         matrix_a, matrix_b, matrix_c, t_a, t_b,&
         dbcsr_scalar (REAL(alpha, real_8)), dbcsr_scalar (REAL(beta, real_8)),&
         io_unit = io_unit, error=error)
    CALL dbcsr_release (matrix_a)
    CALL dbcsr_release (matrix_b)
    CALL dbcsr_release (matrix_c)
  END SUBROUTINE cp_test_multiplies

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
       transa, transb, alpha, beta, limits, retain_sparsity,&
       io_unit, error)
    INTEGER, INTENT(IN)                      :: mp_group
    INTEGER, DIMENSION(:, :)                 :: group_sizes
    TYPE(dbcsr_obj), INTENT(in)              :: matrix_a, matrix_b, matrix_c
    CHARACTER, INTENT(in)                    :: transa, transb
    TYPE(dbcsr_scalar_type), INTENT(in)      :: alpha, beta
    INTEGER, DIMENSION(6), INTENT(in), &
      OPTIONAL                               :: limits
    LOGICAL, INTENT(in), OPTIONAL            :: retain_sparsity
    INTEGER, INTENT(IN)                      :: io_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'test_multiplies_multiproc', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: group, mynode, numnodes, &
                                                pcol, prow, test
    INTEGER(kind=int_8)                      :: flop, flop_sum
    CHARACTER(len=*), PARAMETER :: fmt_per_run_total = &
      '(A,1X,I5,1X,A,1X,F9.3,1X,"s,",1X,EN12.3,1X,"FLOP/s",1X,E13.5)', &
      fmt_per_run_local = &
      '(A,1X,I5,1X,A,1X,F9.3,1X,"s,",1X,EN12.3,1X,"FLOP/s")'

    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: pgrid
    INTEGER, DIMENSION(2)                    :: myploc, npdims
    LOGICAL                                  :: i_am_alive
    REAL(kind=real_8)                        :: cs, flops, flops_all, t1, t2, &
                                                t_max
    TYPE(array_i1d_obj)                      :: col_dist_a, col_dist_b, &
                                                col_dist_c, row_dist_a, &
                                                row_dist_b, row_dist_c
    TYPE(dbcsr_distribution_obj)             :: dist_a, dist_b, dist_c
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_mp_obj)                       :: mp_env
    TYPE(dbcsr_obj)                          :: m_a, m_b, m_c

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (SIZE(group_sizes, 2), "EQ", 2, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN,&
         "second dimension of group_sizes must be 2",__LINE__,dbcsr_error)
    CALL dbcsr_init (m_a)
    CALL dbcsr_init (m_b)
    CALL dbcsr_init (m_c)
    DO test = 1, SIZE(group_sizes, 1)
       npdims(1:2) = group_sizes(test, 1:2)
       numnodes = npdims(1) * npdims(2)
       CALL dbcsr_assert (numnodes, "GE", 0, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Cartesian sides must be greater or equal to 0",__LINE__,dbcsr_error)
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
            dist_a, dbcsr_type_normal,&
            dbcsr_row_block_sizes (matrix_a),&
            dbcsr_col_block_sizes (matrix_a),&
            data_type=dbcsr_get_data_type (matrix_a),&
            error=dbcsr_error)
       CALL dbcsr_distribution_release (dist_a)
       CALL dbcsr_redistribute (matrix_a, m_a, error=dbcsr_error)
       ! B
       CALL dbcsr_create (m_b, "Test for "//TRIM(dbcsr_name (matrix_b)),&
            dist_b, dbcsr_type_normal,&
            dbcsr_row_block_sizes (matrix_b),&
            dbcsr_col_block_sizes (matrix_b),&
            data_type=dbcsr_get_data_type (matrix_b),&
            error=dbcsr_error)
       CALL dbcsr_distribution_release (dist_b)
       CALL dbcsr_redistribute (matrix_b, m_b, error=dbcsr_error)
       ! C
       CALL dbcsr_create (m_c, "Test for "//TRIM(dbcsr_name (matrix_c)),&
            dist_c, dbcsr_type_normal,&
            dbcsr_row_block_sizes (matrix_c),&
            dbcsr_col_block_sizes (matrix_c),&
            data_type=dbcsr_get_data_type (matrix_c),&
            error=dbcsr_error)
       CALL dbcsr_distribution_release (dist_c)
       CALL dbcsr_redistribute (matrix_c, m_c, error=dbcsr_error)
       ! Perform multiply
       CALL mp_sync (group)
       t1 = m_walltime()
       IF (PRESENT (limits)) THEN
          CALL dbcsr_multiply (transa, transb, alpha,&
               m_a, m_b, beta, m_c,&
               first_row = limits(1),&
               last_row = limits(2),&
               first_column = limits(3),&
               last_column = limits(4),&
               first_k = limits(5),&
               last_k = limits(6),&
               retain_sparsity=retain_sparsity, flop=flop, error=dbcsr_error)
       ELSE
          CALL dbcsr_multiply (transa, transb, alpha,&
               m_a, m_b, beta, m_c,&
               retain_sparsity=retain_sparsity, flop=flop, error=dbcsr_error)
       ENDIF
       t2 = m_walltime()
       t_max = t2 - t1
       CALL dmp_max (t_max, group)
       flop_sum = flop
       CALL mp_sum (flop_sum, group)
       t_max = MAX (t_max, 0.001_real_8)
       flops_all = REAL(flop_sum, KIND=real_8) / t_max
       cs = dbcsr_checksum (m_c, error=dbcsr_error)
       ! Release
       IF (io_unit .GT. 0) THEN
          flops = REAL(flop, KIND=real_8) / (t2-t1)
          !WRITE(io_unit,fmt_per_run_local)&
          !     "Local Multiplication",&
          !     numnodes,"processors:",&
          !     t2-t1, flops
          WRITE(io_unit,fmt_per_run_total)&
               "Total Multiplication",&
               numnodes,"processors:",&
               t_max, flops_all, cs
       ENDIF
       CALL dbcsr_mp_release (mp_env)
       CALL dbcsr_release (m_a)
       CALL dbcsr_release (m_b)
       CALL dbcsr_release (m_c)
    ENDDO
  END SUBROUTINE test_multiplies_multiproc

END MODULE dbcsr_tests
