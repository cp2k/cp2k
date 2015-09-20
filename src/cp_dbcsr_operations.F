!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR operations in CP2K
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.8
!>
!> <b>Modification history:</b>
!> - Created 2009-05-12
! *****************************************************************************
MODULE cp_dbcsr_operations
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type,&
                                             get_blacs_info
  USE cp_cfm_types,                    ONLY: cp_cfm_type
  USE cp_dbcsr_interface,              ONLY: &
       convert_sizes_to_offsets, cp_add_work_coordinate, &
       cp_create_bl_distribution, cp_dbcsr_add, cp_dbcsr_col_block_sizes, &
       cp_dbcsr_complete_redistribute, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_distribution, cp_dbcsr_distribution_new, cp_dbcsr_finalize, &
       cp_dbcsr_get_data_p, cp_dbcsr_get_data_size, cp_dbcsr_get_data_type, &
       cp_dbcsr_get_info, cp_dbcsr_get_matrix_type, cp_dbcsr_init, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_mp_new, cp_dbcsr_mp_release, &
       cp_dbcsr_multiply, cp_dbcsr_name, cp_dbcsr_norm, cp_dbcsr_release, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_scale, cp_dbcsr_type, &
       cp_dbcsr_valid_index, cp_dbcsr_verify_matrix, cp_dbcsr_work_create, &
       cp_set_work_size, dbcsr_create_dist_block_cyclic, &
       dbcsr_create_dist_r_unrot, dbcsr_distribution_col_clusters, &
       dbcsr_distribution_col_dist, dbcsr_distribution_init, &
       dbcsr_distribution_local_cols, dbcsr_distribution_local_rows, &
       dbcsr_distribution_mp, dbcsr_distribution_ncols, &
       dbcsr_distribution_nlocal_cols, dbcsr_distribution_nlocal_rows, &
       dbcsr_distribution_nrows, dbcsr_distribution_obj, &
       dbcsr_distribution_release, dbcsr_distribution_row_clusters, &
       dbcsr_distribution_row_dist, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_mp_hold, dbcsr_mp_npcols, dbcsr_mp_nprows, &
       dbcsr_mp_obj, dbcsr_norm_frobenius, dbcsr_type_antisymmetric, &
       dbcsr_type_complex_8, dbcsr_type_no_symmetry, dbcsr_type_real_4, &
       dbcsr_type_real_8, dbcsr_type_symmetric
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_get,&
                                             distribution_2d_type
  USE kinds,                           ONLY: dp,&
                                             real_4,&
                                             sp
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
#include "common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_operations'


  PUBLIC :: cp_dbcsr_multiply_local

  INTERFACE cp_dbcsr_multiply_local
     MODULE PROCEDURE cp_dbcsr_multiply_local_d,&
                      cp_dbcsr_multiply_local_s
  END INTERFACE

  ! CP2K API emulation
  PUBLIC :: copy_fm_to_dbcsr, copy_dbcsr_to_fm,&
            copy_dbcsr_to_cfm, copy_cfm_to_dbcsr,&
            cp_dbcsr_sm_fm_multiply, cp_dbcsr_plus_fm_fm_t,&
            copy_dbcsr_to_fm_bc, copy_fm_to_dbcsr_bc, cp_fm_to_dbcsr_row_template,&
            cp_dbcsr_m_by_n_from_template, cp_dbcsr_m_by_n_from_row_template

  ! distribution_2d_type compatibility
  PUBLIC :: cp_dbcsr_dist2d_to_dist

  PUBLIC :: cp_dbcsr_copy_columns_hack

  INTERFACE cp_dbcsr_plus_fm_fm_t
     MODULE PROCEDURE cp_dbcsr_plus_fm_fm_t_native
  END INTERFACE

  PRIVATE

CONTAINS

! *****************************************************************************
!> \brief hack for dbcsr_copy_columns
!> \param matrix_b ...
!> \param matrix_a ...
!> \param ncol ...
!> \param source_start ...
!> \param target_start ...
!> \param para_env ...
!> \param blacs_env ...
!> \author vw
! *****************************************************************************
  SUBROUTINE cp_dbcsr_copy_columns_hack(matrix_b, matrix_a,&
       ncol, source_start, target_start, para_env, blacs_env)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_b
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    INTEGER, INTENT(IN)                      :: ncol, source_start, &
                                                target_start
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_copy_columns_hack', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nfullcols_total, &
                                                nfullrows_total
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_matrix_a, fm_matrix_b

    NULLIFY(fm_matrix_a, fm_matrix_b, fm_struct)
    CALL cp_dbcsr_get_info(matrix_a,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env)
    CALL cp_fm_create(fm_matrix_a,fm_struct,name="fm_matrix_a")
    CALL cp_fm_struct_release(fm_struct)

    CALL cp_dbcsr_get_info(matrix_b,nfullrows_total=nfullrows_total,nfullcols_total=nfullcols_total)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=nfullrows_total,&
         ncol_global=nfullcols_total,para_env=para_env)
    CALL cp_fm_create(fm_matrix_b,fm_struct,name="fm_matrix_b")
    CALL cp_fm_struct_release(fm_struct)

    CALL copy_dbcsr_to_fm(matrix_a, fm_matrix_a)
    CALL copy_dbcsr_to_fm(matrix_b, fm_matrix_b)

    CALL cp_fm_to_fm(fm_matrix_a, fm_matrix_b, ncol, source_start, target_start)

    CALL copy_fm_to_dbcsr(fm_matrix_b, matrix_b)

    CALL cp_fm_release(fm_matrix_a)
    CALL cp_fm_release(fm_matrix_b)

  END SUBROUTINE cp_dbcsr_copy_columns_hack


! *****************************************************************************
!> \brief Creates a DBCSR distribution from a distribution_2d
!> \param[in] dist2d          distribution_2d
!> \param[out] dist           DBCSR distribution
!> \param unit_nr ...
!> \param mp_obj ...
!> \par History
!>    move form dbcsr_operation 01.2010
! *****************************************************************************
  SUBROUTINE cp_dbcsr_dist2d_to_dist(dist2d, dist, unit_nr,mp_obj)
    TYPE(distribution_2d_type), INTENT(IN), &
      TARGET                                 :: dist2d
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    INTEGER, INTENT(IN)                      :: unit_nr
    TYPE(dbcsr_mp_obj), INTENT(IN), OPTIONAL :: mp_obj

    INTEGER                                  :: mypcol, myproc, myprow, &
                                                numproc
    INTEGER, DIMENSION(:, :), POINTER        :: col_dist_data, pgrid, &
                                                row_dist_data
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_mp_obj)                       :: mp_env
    TYPE(distribution_2d_type), POINTER      :: dist2d_p

!
!   ---------------------------------------------------------------------------

    dist2d_p => dist2d
    CALL distribution_2d_get(dist2d_p,&
         row_distribution=row_dist_data, col_distribution=col_dist_data,&
         blacs_env=blacs_env)
    CALL get_blacs_info(blacs_env, para_env=para_env,&
         my_process_row=myprow, my_process_column=mypcol,&
         blacs2mpi=pgrid)
    myproc = para_env%mepos
    numproc = para_env%num_pe
    IF (PRESENT (mp_obj)) THEN
       mp_env = mp_obj
       CALL dbcsr_mp_hold (mp_env)
    ELSE
       CALL cp_dbcsr_mp_new(mp_env, pgrid, para_env%group, myproc, numproc,&
            myprow, mypcol)
    ENDIF
    CALL cp_dbcsr_distribution_new(dist, mp_env, row_dist_data, col_dist_data, unit_nr=unit_nr)
    CALL cp_dbcsr_mp_release (mp_env)
  END SUBROUTINE cp_dbcsr_dist2d_to_dist

! *****************************************************************************
!> \brief multiply a dbcsr with a replicated array
!>        c = alpha_scalar * A (dbscr) * b + c
!> \param[in] matrix_a DBSCR matrxx
!> \param[in]  vec_b        vectors b
!> \param[inout] vec_c      vectors c
!> \param[in]  ncol         nbr of columns
!> \param[in]  alpha        alpha
!>
! *****************************************************************************
  SUBROUTINE cp_dbcsr_multiply_local_d(matrix_a, vec_b, vec_c, ncol, alpha)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: vec_b
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: vec_c
    INTEGER, INTENT(in), OPTIONAL            :: ncol
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_local_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, coloff, my_ncol, &
                                                row, rowoff, timing_handle
    LOGICAL                                  :: has_symm
    REAL(dp)                                 :: my_alpha, my_alpha2
    REAL(dp), DIMENSION(:, :), POINTER       :: data_d
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN, timing_handle)


    my_alpha = 1.0_dp
    IF (PRESENT(alpha)) my_alpha = alpha

    my_ncol = SIZE(vec_b,2)
    IF(PRESENT(ncol)) my_ncol = ncol

    my_alpha2 = 0.0_dp
    IF(cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric) my_alpha2 = my_alpha
    IF(cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric) my_alpha2 = -my_alpha

    has_symm=(cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric.OR.&
          cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric)

    !$OMP PARALLEL DEFAULT(NONE) SHARED(matrix_a,vec_b,vec_c,ncol,my_alpha2,my_alpha,my_ncol,has_symm) &
    !$OMP          PRIVATE(iter,row,col,data_d,blk,rowoff,coloff)
    CALL cp_dbcsr_iterator_start(iter, matrix_a, read_only=.TRUE., dynamic = .TRUE., dynamic_byrows=.TRUE.)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, data_d, blk, row_offset=rowoff, col_offset=coloff)
       IF (my_ncol.ne.1) THEN
          CALL dgemm('N','N',&
               SIZE(data_d,1),my_ncol,SIZE(data_d,2),&
               my_alpha, data_d(1,1),          SIZE(data_d,1),&
               vec_b(coloff,1), SIZE(vec_b,1), &
               1.0_dp,   vec_c(rowoff,1), SIZE(vec_c,1))
       ELSE
          CALL dgemv('N',SIZE(data_d,1),SIZE(data_d,2),&
               my_alpha, data_d(1,1),          SIZE(data_d,1),&
               vec_b(coloff,1),1,&
               1.0_dp,   vec_c(rowoff,1),1)
       ENDIF
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)
    !$OMP END PARALLEL

    ! FIXME ... in the symmetric case, the writes to vec_c depend on the column, not the row. This makes OMP-ing more difficult
    ! needs e.g. a buffer for vec_c and a reduction of that buffer.
    IF (has_symm) THEN
       CALL cp_dbcsr_iterator_start(iter, matrix_a)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, row, col, data_d, blk, row_offset=rowoff, col_offset=coloff)
          IF(row.NE.col) THEN
             IF (my_ncol.ne.1) THEN
                CALL dgemm('T','N', &
                     SIZE(data_d,2), my_ncol, SIZE(data_d,1),&
                     my_alpha2, data_d(1,1),          SIZE(data_d,1), &
                     vec_b(rowoff,1), SIZE(vec_b,1),  &
                     1.0_dp,    vec_c(coloff,1), SIZE(vec_c,1))
             ELSE
                CALL dgemv('T',SIZE(data_d,1),SIZE(data_d,2),&
                     my_alpha2, data_d(1,1),          SIZE(data_d,1), &
                     vec_b(rowoff,1),1,  &
                     1.0_dp,    vec_c(coloff,1), 1)
             ENDIF
          ENDIF
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
    ENDIF

    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_local_d

! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param vec_b ...
!> \param vec_c ...
!> \param ncol ...
!> \param alpha ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_multiply_local_s(matrix_a, vec_b, vec_c, ncol, alpha)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    REAL(KIND=real_4), DIMENSION(:, :), &
      INTENT(IN)                             :: vec_b
    REAL(KIND=real_4), DIMENSION(:, :), &
      INTENT(INOUT)                          :: vec_c
    INTEGER, INTENT(in), OPTIONAL            :: ncol
    REAL(KIND=real_4), INTENT(IN), OPTIONAL  :: alpha

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_local_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, coloff, my_ncol, &
                                                row, rowoff, timing_handle
    REAL(KIND=real_4)                        :: my_alpha, my_alpha2
    REAL(KIND=real_4), DIMENSION(:, :), &
      POINTER                                :: data_d
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN, timing_handle)

    my_alpha = 1.0
    IF (PRESENT(alpha)) my_alpha = alpha

    my_ncol = SIZE(vec_b,2)
    IF(PRESENT(ncol)) my_ncol = ncol

    my_alpha2 = 0.0
    IF(cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric) my_alpha2 = my_alpha
    IF(cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric) my_alpha2 = -my_alpha

    CALL cp_dbcsr_iterator_start(iter, matrix_a)

    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

       CALL cp_dbcsr_iterator_next_block(iter, row, col, data_d, blk, row_offset=rowoff, col_offset=coloff)

       CALL sgemm('N','N',&
                  SIZE(data_d,1),my_ncol,SIZE(data_d,2),&
                  my_alpha, data_d(1,1),          SIZE(data_d,1),&
                            vec_b(coloff,1), SIZE(vec_b,1), &
                  1.0,      vec_c(rowoff,1), SIZE(vec_c,1))

       IF((cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric.OR.&
          cp_dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric)) THEN
          IF(row.NE.col) THEN
             CALL sgemm('T','N', &
                        SIZE(data_d,2), my_ncol, SIZE(data_d,1),&
                        my_alpha2, data_d(1,1),          SIZE(data_d,1), &
                                   vec_b(rowoff,1), SIZE(vec_b,1),  &
                        1.0,       vec_c(coloff,1), SIZE(vec_c,1))
          ENDIF
       ENDIF
    ENDDO

    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_local_s



! *****************************************************************************
!> \brief multiply a dbcsr with a fm matrix
!>
!> For backwards compatibility with BLAS XGEMM, this routine supports
!> the multiplication of matrices with incompatible dimensions.
!>
!> \param[in]  matrix         DBCSR matrix
!> \param fm_in full matrix
!> \param fm_out full matrix
!> \param[in]  ncol           nbr of columns
!> \param[in]  alpha          alpha
!> \param[in]  beta           beta
!>
! *****************************************************************************
  SUBROUTINE cp_dbcsr_sm_fm_multiply(matrix, fm_in, fm_out, ncol, alpha, beta)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_fm_type), POINTER                :: fm_in, fm_out
    INTEGER, INTENT(IN)                      :: ncol
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha, beta

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_sm_fm_multiply', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k_in, k_out, timing_handle, &
                                                timing_handle_mult
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size_right_in, &
                                                col_blk_size_right_out
    TYPE(cp_dbcsr_type)                      :: in, out
    REAL(dp)                                 :: my_alpha, my_beta
    TYPE(dbcsr_distribution_obj)             :: dist_right_in, product_dist

    CALL timeset(routineN, timing_handle)

    my_alpha=1.0_dp
    my_beta=0.0_dp
    IF (PRESENT(alpha)) my_alpha=alpha
    IF (PRESENT(beta)) my_beta=beta

    CALL cp_fm_get_info(fm_out, ncol_global=k_out)

    CALL cp_fm_get_info(fm_in, ncol_global=k_in)
    !write(*,*)routineN//" -----------------------------------"
    !IF (k_in .NE. k_out) &
    !   WRITE(*,'(3(A,I5,1X),2(A,F5.2,1X))')&
    !   routineN//" ncol", ncol,'k_in',k_in,'k_out',k_out,&
    !   'alpha',my_alpha,'beta',my_beta

    IF (ncol.GT.0.AND.k_out.GT.0.AND.k_in.GT.0) THEN
       CALL dbcsr_create_dist_r_unrot (dist_right_in, cp_dbcsr_distribution(matrix), k_in, &
            col_blk_size_right_in)

       CALL cp_dbcsr_init(in)
       CALL cp_dbcsr_create(in, "D", dist_right_in, dbcsr_type_no_symmetry, &
            cp_dbcsr_row_block_sizes(matrix), col_blk_size_right_in,&
            nze=0)

       CALL cp_dbcsr_init(out)
       CALL cp_dbcsr_distribution_new (product_dist,&
            dbcsr_distribution_mp (cp_dbcsr_distribution(matrix)),&
            dbcsr_distribution_row_dist (cp_dbcsr_distribution(matrix)),&
            dbcsr_distribution_col_dist (dist_right_in),&
            dbcsr_distribution_row_clusters (cp_dbcsr_distribution(matrix)),&
            dbcsr_distribution_col_clusters (dist_right_in))
       ALLOCATE(col_blk_size_right_out(SIZE(col_blk_size_right_in)))
       col_blk_size_right_out = col_blk_size_right_in
       CALL match_col_sizes (col_blk_size_right_out, col_blk_size_right_in, k_out)

       !if (k_in .ne. k_out) then
       !   write(*,*)routineN//" in cs", col_blk_size_right_in
       !   write(*,*)routineN//" out cs", col_blk_size_right_out
       !endif

       CALL cp_dbcsr_create(out, "D", product_dist, dbcsr_type_no_symmetry, &
            cp_dbcsr_row_block_sizes(matrix), col_blk_size_right_out,&
            nze=0)

       CALL copy_fm_to_dbcsr(fm_in, in)
       IF(ncol.NE.k_out.OR.my_beta.NE.0.0_dp) &
            CALL copy_fm_to_dbcsr(fm_out, out)

       CALL timeset(routineN//'_core', timing_handle_mult)
       CALL cp_dbcsr_multiply("N", "N", my_alpha, matrix, in, my_beta, out,&
            last_column=ncol, match_matrix_sizes=.TRUE.)
       CALL timestop(timing_handle_mult)

       CALL copy_dbcsr_to_fm(out, fm_out)

       CALL cp_dbcsr_release(in)
       CALL cp_dbcsr_release(out)
       DEALLOCATE(col_blk_size_right_in,col_blk_size_right_out)
       CALL dbcsr_distribution_release(dist_right_in)
       CALL dbcsr_distribution_release(product_dist)

    ENDIF

    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_sm_fm_multiply

! *****************************************************************************
!> \brief ...
!> \param sizes1 ...
!> \param sizes2 ...
!> \param full_num ...
! *****************************************************************************
  SUBROUTINE match_col_sizes (sizes1, sizes2, full_num)
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: sizes1
    INTEGER, DIMENSION(:), INTENT(IN)        :: sizes2
    INTEGER, INTENT(IN)                      :: full_num

    INTEGER                                  :: left, n1, n2, p, rm, used

    n1 = SIZE(sizes1)
    n2 = SIZE(sizes2)
    CALL cp_assert (n1 .EQ. n2, cp_fatal_level, cp_caller_error,&
         "match_col_sizes", "distributions must be equal!")
    sizes1(1:n1) = sizes2(1:n1)
    used = SUM (sizes1(1:n1))
    ! If sizes1 does not cover everything, then we increase the
    ! size of the last block; otherwise we reduce the blocks
    ! (from the end) until it is small enough.
    IF (used .LT. full_num) THEN
       sizes1(n1) = sizes1(n1) + full_num-used
    ELSE
       left = used - full_num
       p = n1
       DO WHILE (left .GT. 0 .AND. p .GT. 0)
          rm = MIN(left, sizes1(p))
          sizes1(p) = sizes1(p) - rm
          left = left - rm
          p = p-1
       ENDDO
    ENDIF
  END SUBROUTINE match_col_sizes
! *****************************************************************************
!> \brief performs the multiplication sparse_matrix+dense_mat*dens_mat^T
!>        if matrix_g is not explicitly given, matrix_v^T will be used
!>        this can be important to save the necessary redistribute for a 
!>        different matrix_g and increase performance.
!> \param sparse_matrix ...
!> \param matrix_v ...
!> \param matrix_g ...
!> \param ncol ...
!> \param alpha ...
!> \param keep_sparsity Determines if the sparsity of sparse_matrix is retained
!>        by default it is TRUE
! *****************************************************************************
  SUBROUTINE cp_dbcsr_plus_fm_fm_t_native(sparse_matrix,matrix_v,matrix_g,ncol,&
       alpha,keep_sparsity)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: sparse_matrix
    TYPE(cp_fm_type), POINTER                :: matrix_v
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: matrix_g
    INTEGER, INTENT(IN)                      :: ncol
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: alpha
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_plus_fm_fm_t_native', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k, nao, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size_left, &
                                                col_dist_cluster, &
                                                col_dist_left
    LOGICAL                                  :: check_product, &
                                                my_keep_sparsity
    REAL(KIND=dp)                            :: my_alpha, norm
    TYPE(cp_dbcsr_type)                      :: mat_g, mat_v, sparse_matrix2, &
                                                sparse_matrix3
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: fm_matrix
    TYPE(dbcsr_distribution_obj)             :: dist_left
    TYPE(dbcsr_mp_obj)                       :: mp

    check_product = .FALSE.

    CALL timeset(routineN, timing_handle)

    my_keep_sparsity=.TRUE.
    IF(PRESENT(keep_sparsity)) my_keep_sparsity=keep_sparsity
    NULLIFY(col_dist_left)

    IF (ncol .GT. 0) THEN
       CALL cp_assert (cp_dbcsr_valid_index (sparse_matrix), cp_fatal_level,&
            cp_caller_error, routineN, "sparse_matrix must pre-exist")
       !
       ! Setup matrix_v
       CALL cp_fm_get_info(matrix_v, ncol_global=k)
       !WRITE(*,*)routineN//'truncated mult k, ncol',k,ncol,' PRESENT (matrix_g)',PRESENT (matrix_g)
       mp = dbcsr_distribution_mp (cp_dbcsr_distribution(sparse_matrix))
       CALL cp_create_bl_distribution (col_dist_left, col_blk_size_left,&
            k, dbcsr_mp_npcols (mp))
       NULLIFY(col_dist_cluster)
       CALL cp_dbcsr_distribution_new (dist_left, mp,&
            dbcsr_distribution_row_dist (cp_dbcsr_distribution(sparse_matrix)),&
            col_dist_left,&
            dbcsr_distribution_row_clusters(cp_dbcsr_distribution(sparse_matrix)),&
            col_dist_cluster)
       DEALLOCATE(col_dist_left)
       CALL cp_dbcsr_init (mat_v)
       CALL cp_dbcsr_create(mat_v, "DBCSR matrix_v", dist_left, dbcsr_type_no_symmetry,&
            cp_dbcsr_row_block_sizes (sparse_matrix), col_blk_size_left, nze=0,&
            data_type=cp_dbcsr_get_data_type (sparse_matrix))
       CALL copy_fm_to_dbcsr(matrix_v, mat_v)
       CALL cp_dbcsr_verify_matrix(mat_v)
       !
       ! Setup matrix_g
       IF(PRESENT (matrix_g)) THEN
          CALL cp_dbcsr_init(mat_g)
          CALL cp_dbcsr_create(mat_g, "DBCSR matrix_g", dist_left,&
               dbcsr_type_no_symmetry,&
               cp_dbcsr_row_block_sizes (sparse_matrix),&
               cp_dbcsr_col_block_sizes (mat_v),&
               data_type=cp_dbcsr_get_data_type (sparse_matrix))
          CALL copy_fm_to_dbcsr(matrix_g, mat_g)
       ENDIF
       !
       DEALLOCATE (col_blk_size_left)
       CALL dbcsr_distribution_release (dist_left)
       !
       !
       IF(check_product) THEN
          NULLIFY(fm_matrix)
          CALL cp_fm_get_info(matrix_v,nrow_global=nao)
          CALL cp_fm_struct_create(fm_struct_tmp,context=matrix_v%matrix_struct%context,nrow_global=nao,&
               ncol_global=nao,para_env=matrix_v%matrix_struct%para_env)
          CALL cp_fm_create(fm_matrix,fm_struct_tmp,name="fm matrix")
          CALL cp_fm_struct_release(fm_struct_tmp)
          CALL copy_dbcsr_to_fm(sparse_matrix,fm_matrix)
          CALL cp_dbcsr_init(sparse_matrix3)
          CALL cp_dbcsr_copy(sparse_matrix3,sparse_matrix)
       ENDIF
       !
       my_alpha = 1.0_dp
       IF(PRESENT (alpha)) my_alpha = alpha
       IF(PRESENT (matrix_g)) THEN
          CALL cp_dbcsr_multiply("N", "T", my_alpha, mat_v, mat_g,&
               1.0_dp, sparse_matrix,&
               retain_sparsity=my_keep_sparsity,&
               last_k = ncol)
       ELSE
          CALL cp_dbcsr_multiply("N", "T", my_alpha, mat_v, mat_v,&
               1.0_dp, sparse_matrix,&
               retain_sparsity=my_keep_sparsity,&
               last_k = ncol)
       ENDIF

       IF(check_product) THEN
          IF(PRESENT (matrix_g)) THEN
             CALL cp_fm_gemm("N","T",nao,nao,ncol,my_alpha,matrix_v,matrix_g,&
                  1.0_dp,fm_matrix)
          ELSE
             CALL cp_fm_gemm("N","T",nao,nao,ncol,my_alpha,matrix_v,matrix_v,&
                  1.0_dp,fm_matrix)
          ENDIF

          CALL cp_dbcsr_init(sparse_matrix2)
          CALL cp_dbcsr_copy(sparse_matrix2,sparse_matrix)
          CALL cp_dbcsr_scale(sparse_matrix2,alpha_scalar=0.0_dp)
          CALL copy_fm_to_dbcsr(fm_matrix,sparse_matrix2,keep_sparsity=my_keep_sparsity)
          CALL cp_dbcsr_add(sparse_matrix2,sparse_matrix,alpha_scalar=1.0_dp,&
               beta_scalar=-1.0_dp)
          CALL cp_dbcsr_norm(sparse_matrix2,which_norm=dbcsr_norm_frobenius,&
               norm_scalar=norm)
          WRITE(*,*) 'nao=',nao,' k=',k,' ncol=',ncol,' my_alpha=',my_alpha
          WRITE(*,*) 'PRESENT (matrix_g)',PRESENT (matrix_g)
          WRITE(*,*) 'matrix_type=',cp_dbcsr_get_matrix_type(sparse_matrix)
          WRITE(*,*) 'norm(sm+alpha*v*g^t - fm+alpha*v*g^t)/n=',norm/REAL(nao,dp)
          IF(norm/REAL(nao,dp).GT.1e-12_dp) THEN
             !WRITE(*,*) 'fm_matrix'
             !DO j=1,SIZE(fm_matrix%local_data,2)
             !   DO i=1,SIZE(fm_matrix%local_data,1)
             !      WRITE(*,'(A,I3,A,I3,A,E26.16,A)') 'a(',i,',',j,')=',fm_matrix%local_data(i,j),';'
             !   ENDDO
             !ENDDO
             !WRITE(*,*) 'mat_v'
             !CALL cp_dbcsr_print(mat_v,matlab_format=.TRUE.)
             !WRITE(*,*) 'mat_g'
             !CALL cp_dbcsr_print(mat_g,matlab_format=.TRUE.)
             !WRITE(*,*) 'sparse_matrix'
             !CALL cp_dbcsr_print(sparse_matrix,matlab_format=.TRUE.)
             !WRITE(*,*) 'sparse_matrix2 (-sm + sparse(fm))'
             !CALL cp_dbcsr_print(sparse_matrix2,matlab_format=.TRUE.)
             !WRITE(*,*) 'sparse_matrix3 (copy of sm input)'
             !CALL cp_dbcsr_print(sparse_matrix3,matlab_format=.TRUE.)
             !stop
          ENDIF
          CALL cp_dbcsr_release(sparse_matrix2)
          CALL cp_dbcsr_release(sparse_matrix3)
          CALL cp_fm_release(fm_matrix)
       ENDIF
       CALL cp_dbcsr_release (mat_v)
       IF(PRESENT (matrix_g)) CALL cp_dbcsr_release (mat_g)
    ENDIF
    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_plus_fm_fm_t_native

! *****************************************************************************
!> \brief   Copy a BLACS matrix to a dbcsr matrix.
!>
!>          real_matrix=beta*real_matrix+alpha*fm
!>          beta defaults to 0, alpha to 1
!> \param[in] fm              full matrix
!> \param[out] matrix         DBCSR matrix
!> \param[in] alpha           (optional) scaling of FM
!> \param[in] beta            (optional) scaling of existing SM
!> \param[in] keep_sparsity   (optional) retains the sparsity of the input
!>                            matrix
!> \date    2009-10-13
!> \par History
!>          2009-10-13 rewritten based on copy_dbcsr_to_fm
!> \author  Urban Borstnik
!> \version 2.0
! *****************************************************************************
  SUBROUTINE copy_fm_to_dbcsr(fm,matrix,alpha,beta,keep_sparsity)
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(kind=dp), INTENT(IN), OPTIONAL      :: alpha, beta
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_fm_to_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk_p, col, col_l, col_size, error_handler, group, handle, &
      nblkcols_local, nblkcols_total, nblkrows_local, nblkrows_total, &
      ncol_block, ncol_global, nfullcols_local, nfullcols_total, &
      nfullrows_local, nfullrows_total, nrow_block, nrow_global, nze, row, &
      row_l, row_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_col, first_row, &
                                                last_col, last_row, &
                                                local_col_sizes, &
                                                local_row_sizes
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, local_cols, &
                                                local_rows, row_blk_size
    REAL(kind=dp)                            :: my_beta
    REAL(KIND=dp), DIMENSION(:), POINTER     :: blk_1d_dp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fm_block
    REAL(KIND=sp), DIMENSION(:), POINTER     :: blk_1d_sp
    REAL(kind=sp), DIMENSION(:, :), POINTER  :: fm_block_sp
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_dbcsr_type)                      :: bc_mat
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, dbcsr_error)
    CALL timeset(routineN,handle)

    NULLIFY(col_blk_size,row_blk_size)

    my_beta=0._dp
    IF (PRESENT(beta)) THEN
       CALL cp_assert (beta .EQ. my_beta, cp_fatal_level,&
            cp_unimplemented_error_nr, routineN,&
            "beta not supported, use matrix addition instead")
       my_beta=beta
    ENDIF
    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global
    ncol_global = fm%matrix_struct%ncol_global

    CALL cp_dbcsr_get_info(matrix,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total)

    CALL cp_assert (nrow_global.eq.nfullrows_total,&
         cp_fatal_level, cp_caller_error, routineN,&
         "FM and DBCSR matrix sizes do not match in rows")
    CALL cp_assert (ncol_global.eq.nfullcols_total,&
         cp_fatal_level, cp_caller_error, routineN,&
         "FM and DBCSR matrix sizes do not match in columns")

    ! Create a block-cyclic distribution compatible with the FM matrix.
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,& ! Actual full matrix size
         nrow_block, ncol_block,&           ! BLACS parameters
         dbcsr_distribution_mp (cp_dbcsr_distribution (matrix)),&
         row_blk_size, col_blk_size)        ! block-cyclic row/col sizes

    ! Create the block-cyclic DBCSR matrix
    CALL cp_dbcsr_init (bc_mat)
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic "//cp_dbcsr_name(matrix), bc_dist,&
         cp_dbcsr_get_matrix_type(matrix), row_blk_size, col_blk_size, nze=0,&
         reuse_arrays=.TRUE., data_type=cp_dbcsr_get_data_type(matrix))

    !call dbcsr_finalize (bc_mat)
    CALL dbcsr_distribution_release (bc_dist)
    NULLIFY(col_blk_size,row_blk_size)

    CALL cp_dbcsr_get_info(bc_mat,&
         nblkrows_total=nblkrows_total,&
         nblkcols_total=nblkcols_total,&
         nblkrows_local=nblkrows_local,&
         nblkcols_local=nblkcols_local,&
         nfullrows_local=nfullrows_local,&
         nfullcols_local=nfullcols_local,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total,&
         local_rows=local_rows,&
         local_cols=local_cols,&
         row_blk_size=row_blk_size,&
         col_blk_size=col_blk_size)

    ALLOCATE (local_row_sizes (nblkrows_total))
    ALLOCATE (local_col_sizes (nblkcols_total))
    local_row_sizes(:) = 0
    IF (nblkrows_local .GE. 1) THEN
       DO row = 1 , nblkrows_local
          local_row_sizes(local_rows(row)) = row_blk_size(local_rows(row))
       END DO
    ENDIF
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       DO col = 1 , nblkcols_local
          local_col_sizes(local_cols(col)) = col_blk_size(local_cols(col))
       END DO
    ENDIF

    ALLOCATE (first_row(nblkrows_total),last_row(nblkrows_total))
    ALLOCATE (first_col(nblkcols_total),last_col(nblkcols_total))
    CALL convert_sizes_to_offsets (local_row_sizes, first_row, last_row)
    CALL convert_sizes_to_offsets (local_col_sizes, first_col, last_col)

    ! Copy the FM data to the block-cyclic DBCSR matrix.  This step
    ! could be skipped with appropriate DBCSR index manipulation.
    fm_block => fm%local_data
    fm_block_sp => fm%local_data_sp

    CALL cp_dbcsr_work_create (bc_mat, nblks_guess=nblkrows_local*nblkcols_local,&
         sizedata_guess=nfullrows_local*nfullcols_local, work_mutable=.FALSE.,&
         n=1)
    blk_p = 1
    bc_rows: DO row_l = 1, nblkrows_local
       row = local_rows (row_l)
       row_size = row_blk_size(row)
       bc_cols: DO col_l = 1, nblkcols_local
          col = local_cols (col_l)
          col_size = col_blk_size(col)
          nze = row_size*col_size
          !WRITE(*,*)routineN//" Adding block",row,col,"size",nze
          CALL cp_add_work_coordinate(bc_mat, 1, row, col, blk_p, error=dbcsr_error)
          IF (fm%use_sp) THEN
             !blk_1d_sp => bc_mat%m%wms(1)%data_area%d%r_sp(blk_p:blk_p+nze-1)
             blk_1d_sp => cp_dbcsr_get_data_p (bc_mat, 1,&
                  select_data_type=REAL(0.0, KIND=sp), lb=blk_p, ub=blk_p+nze-1)
          ELSE
             !blk_1d_dp => bc_mat%m%wms(1)%data_area%d%r_dp(blk_p:blk_p+nze-1)
             blk_1d_dp => cp_dbcsr_get_data_p (bc_mat, 1,&
                  select_data_type=REAL(0.0, KIND=dp), lb=blk_p, ub=blk_p+nze-1)
          ENDIF
          CALL cp_assert (nze .EQ. (last_row(row)-first_row(row)+1)*(last_col(col)-first_col(col)+1),&
               cp_fatal_level, cp_internal_error, routineN,&
               "Block size does not match block row/col sizes")
          IF (fm%use_sp) THEN
             blk_1d_sp(1:nze) = RESHAPE(&
                  fm_block_sp(&
                     first_row(row):last_row(row),first_col(col):last_col(col)&
                  ), (/ nze /))
          ELSE
             blk_1d_dp(1:nze) = RESHAPE(&
                  fm_block(&
                     first_row(row):last_row(row),first_col(col):last_col(col)&
                  ), (/ nze /))
          ENDIF
          blk_p = blk_p + nze
       ENDDO bc_cols
    ENDDO bc_rows
    CALL cp_set_work_size(bc_mat, 1, blk_p - 1)
    CALL cp_dbcsr_finalize (bc_mat, reshuffle=.FALSE.)

    ! Now convert to the desired matrix distribution
    IF (PRESENT (alpha)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,'no more alpha... clean me')
    ELSE
       CALL cp_dbcsr_complete_redistribute (bc_mat, matrix,&
            keep_sparsity=keep_sparsity)
    ENDIF
    CALL cp_dbcsr_release (bc_mat)

    CALL timestop(handle)
    CALL dbcsr_error_stop(error_handler, dbcsr_error)
  END SUBROUTINE copy_fm_to_dbcsr

! *****************************************************************************
!> \brief ...
!> \param fm ...
!> \param bc_mat ...
!> \param beta ...
! *****************************************************************************
  SUBROUTINE copy_fm_to_dbcsr_bc(fm,bc_mat,beta)
    TYPE(cp_fm_type), POINTER                :: fm
    TYPE(cp_dbcsr_type)                      :: bc_mat
    REAL(kind=dp), INTENT(IN), OPTIONAL      :: beta

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_fm_to_dbcsr_bc', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk_p, col, col_l, col_size, data_type, error_handler, group, &
      handle, my_blk, mypcol, myproc, myprow, nblkcols_local, nblkcols_total, &
      nblkrows_local, nblkrows_total, ncol_block, ncol_global, &
      nfullcols_local, nfullcols_total, nfullrows_local, nfullrows_total, &
      nrow_block, nrow_global, numproc, nze, row, row_l, row_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_col, first_row, &
                                                last_col, last_row, &
                                                local_col_sizes, &
                                                local_row_sizes
    INTEGER, DIMENSION(:), POINTER           :: blk_map, col_blk_size, &
                                                local_cols, local_rows, &
                                                row_blk_size
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    REAL(kind=dp)                            :: my_beta
    REAL(KIND=dp), DIMENSION(:), POINTER     :: blk_1d_dp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fm_block
    REAL(KIND=sp), DIMENSION(:), POINTER     :: blk_1d_sp
    REAL(kind=sp), DIMENSION(:, :), POINTER  :: fm_block_sp
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_mp_obj)                       :: mp_env

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, dbcsr_error)
    CALL timeset(routineN,handle)

    NULLIFY(col_blk_size,row_blk_size)

    my_beta=0._dp
    IF (PRESENT(beta)) THEN
       CALL cp_assert (beta .EQ. my_beta, cp_fatal_level,&
            cp_unimplemented_error_nr, routineN,&
            "beta not supported, use matrix addition instead")
       my_beta=beta
    ENDIF
    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global
    ncol_global = fm%matrix_struct%ncol_global
    myproc = fm%matrix_struct%para_env%mepos
    numproc = fm%matrix_struct%para_env%num_pe
    myprow = context%mepos(1)
    mypcol = context%mepos(2)
    pgrid => context%blacs2mpi

    nfullrows_total=nrow_global
    nfullcols_total=ncol_global
    CALL cp_dbcsr_mp_new(mp_env, pgrid, group, myproc, numproc,&
                      myprow, mypcol)

    ! Create a block-cyclic distribution compatible with the FM matrix.
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,& ! Actual full matrix size
         nrow_block, ncol_block,&           ! BLACS parameters
         mp_env,&
         row_blk_size, col_blk_size)        ! block-cyclic row/col sizes

    CALL cp_dbcsr_mp_release (mp_env)

    ! Create the block-cyclic DBCSR matrix
    data_type=dbcsr_type_real_8
    IF(fm%use_sp)data_type=dbcsr_type_real_4
    CALL cp_dbcsr_init (bc_mat)
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic ", bc_dist,&
         dbcsr_type_no_symmetry, row_blk_size, col_blk_size, nze=0,&
         reuse_arrays=.TRUE., data_type=data_type)

    !call dbcsr_finalize (bc_mat)
    CALL dbcsr_distribution_release (bc_dist)
    NULLIFY(col_blk_size,row_blk_size)

    CALL cp_dbcsr_get_info(bc_mat,&
         nblkrows_total=nblkrows_total,&
         nblkcols_total=nblkcols_total,&
         nblkrows_local=nblkrows_local,&
         nblkcols_local=nblkcols_local,&
         nfullrows_local=nfullrows_local,&
         nfullcols_local=nfullcols_local,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total,&
         local_rows=local_rows,&
         local_cols=local_cols,&
         row_blk_size=row_blk_size,&
         col_blk_size=col_blk_size)

    ALLOCATE (local_row_sizes (nblkrows_total))
    ALLOCATE (local_col_sizes (nblkcols_total))
    local_row_sizes(:) = 0
    IF (nblkrows_local .GE. 1) THEN
       DO row = 1 , nblkrows_local
          local_row_sizes(local_rows(row)) = row_blk_size(local_rows(row))
       END DO
    ENDIF
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       DO col = 1 , nblkcols_local
          local_col_sizes(local_cols(col)) = col_blk_size(local_cols(col))
       END DO
    ENDIF

    ALLOCATE (first_row(nblkrows_total),last_row(nblkrows_total))
    ALLOCATE (first_col(nblkcols_total),last_col(nblkcols_total))
    ALLOCATE (blk_map(0:nblkrows_total*nblkcols_total))
    blk_map(0)=1
    CALL convert_sizes_to_offsets (local_row_sizes, first_row, last_row)
    CALL convert_sizes_to_offsets (local_col_sizes, first_col, last_col)

    ! Copy the FM data to the block-cyclic DBCSR matrix.  This step
    ! could be skipped with appropriate DBCSR index manipulation.
    fm_block => fm%local_data
    fm_block_sp => fm%local_data_sp

    CALL cp_dbcsr_work_create (bc_mat, nblks_guess=nblkrows_local*nblkcols_local,&
         sizedata_guess=nfullrows_local*nfullcols_local, work_mutable=.FALSE.,&
         n=1)
    blk_p = 1
    bc_rows: DO row_l = 1, nblkrows_local
       row = local_rows (row_l)
       row_size = row_blk_size(row)
       bc_cols: DO col_l = 1, nblkcols_local
          col = local_cols (col_l)
          col_size = col_blk_size(col)
          my_blk=(row_l-1)*nblkcols_local+col_l-1
          nze = row_size*col_size
          CALL cp_add_work_coordinate(bc_mat, 1, row, col, blk_map(my_blk), error=dbcsr_error)
          blk_map(my_blk+1)=blk_map(my_blk)+nze
       ENDDO bc_cols
    ENDDO bc_rows

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(my_blk,nze,row_l,col_l,row,col,blk_1d_dp,blk_1d_sp) &
!$OMP SHARED(local_rows, nblkrows_local, nblkcols_local, local_cols, blk_map, fm, bc_mat,&
!$OMP        last_row, first_row, last_col, first_col, fm_block, fm_block_sp)
    bc_rows1: DO row_l = 1, nblkrows_local
       row = local_rows (row_l)
       bc_cols1: DO col_l = 1, nblkcols_local
          col = local_cols (col_l)
          my_blk=(row_l-1)*nblkcols_local+col_l-1
          nze=blk_map(my_blk+1)-blk_map(my_blk)
          IF (fm%use_sp) THEN
             !blk_1d_sp => bc_mat%m%wms(1)%data_area%d%r_sp(blk_p:blk_p+nze-1)
             blk_1d_sp => cp_dbcsr_get_data_p (bc_mat, 1,&
                  select_data_type=REAL(0.0, KIND=sp), lb=blk_map(my_blk), ub=blk_map(my_blk)+nze-1)
          ELSE
             !blk_1d_dp => bc_mat%m%wms(1)%data_area%d%r_dp(blk_p:blk_p+nze-1)
             blk_1d_dp => cp_dbcsr_get_data_p (bc_mat, 1,&
                  select_data_type=REAL(0.0, KIND=dp), lb=blk_map(my_blk), ub=blk_map(my_blk)+nze-1)
          ENDIF
          CALL cp_assert (nze .EQ. (last_row(row)-first_row(row)+1)*(last_col(col)-first_col(col)+1),&
               cp_fatal_level, cp_internal_error, routineN,&
               "Block size does not match block row/col sizes")
          IF (fm%use_sp) THEN
             blk_1d_sp(1:nze) = RESHAPE(&
                  fm_block_sp(&
                     first_row(row):last_row(row),first_col(col):last_col(col)&
                  ), (/ nze /))
          ELSE
             blk_1d_dp(1:nze) = RESHAPE(&
                  fm_block(&
                     first_row(row):last_row(row),first_col(col):last_col(col)&
                  ), (/ nze /))
          ENDIF
       ENDDO bc_cols1
    ENDDO bc_rows1
!$OMP END PARALLEL DO
    CALL cp_set_work_size(bc_mat, 1, blk_map(nblkrows_local*nblkcols_local) - 1)
    CALL cp_dbcsr_finalize (bc_mat, reshuffle=.FALSE.)
    DEALLOCATE(blk_map)

    CALL timestop(handle)
    CALL dbcsr_error_stop(error_handler, dbcsr_error)
  END SUBROUTINE copy_fm_to_dbcsr_bc
! *****************************************************************************
!> \brief   Copy a BLACS matrix to a dbcsr matrix.
!>
!>          real_matrix=beta*real_matrix+alpha*fm
!>          beta defaults to 0, alpha to 1
!> \param[in] fm              full matrix
!> \param[out] matrix         DBCSR matrix
!> \param[in] keep_sparsity   (optional) retains the sparsity of the input
!>                            matrix
!> \date    2010
!> \par History
!>          2010  copied from copy_dbcsr_to_fm
!> \author  VW
!> \version 2.0
! *****************************************************************************
  SUBROUTINE copy_cfm_to_dbcsr(fm,matrix,keep_sparsity)
    TYPE(cp_cfm_type), POINTER               :: fm
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_cfm_to_dbcsr', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: blk_1d
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: fm_block
    INTEGER :: blk_p, col, col_l, col_size, error_handler, group, handle, &
      nblkcols_local, nblkcols_total, nblkrows_local, nblkrows_total, &
      ncol_block, ncol_global, nfullcols_local, nfullcols_total, &
      nfullrows_local, nfullrows_total, nrow_block, nrow_global, nze, row, &
      row_l, row_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_col, first_row, &
                                                last_col, last_row, &
                                                local_col_sizes, &
                                                local_row_sizes
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, local_cols, &
                                                local_rows, row_blk_size
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_dbcsr_type)                      :: bc_mat
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, dbcsr_error)
    CALL timeset(routineN,handle)

    NULLIFY(col_blk_size,row_blk_size)

    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global
    ncol_global = fm%matrix_struct%ncol_global

    CALL cp_dbcsr_get_info(matrix,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total)

    CALL cp_assert (nrow_global.eq.nfullrows_total,&
         cp_fatal_level, cp_caller_error, routineN,&
         "FM and DBCSR matrix sizes do not match in rows")
    CALL cp_assert (ncol_global.eq.nfullcols_total,&
         cp_fatal_level, cp_caller_error, routineN,&
         "FM and DBCSR matrix sizes do not match in columns")

    ! Create a block-cyclic distribution compatible with the FM matrix.
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,& ! Actual full matrix size
         nrow_block, ncol_block,&           ! BLACS parameters
         dbcsr_distribution_mp (cp_dbcsr_distribution (matrix)),&
         row_blk_size, col_blk_size)        ! block-cyclic row/col sizes

    ! Create the block-cyclic DBCSR matrix
    CALL cp_dbcsr_init (bc_mat)
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic "//cp_dbcsr_name(matrix), bc_dist,&
         cp_dbcsr_get_matrix_type(matrix), row_blk_size, col_blk_size, nze=0,&
         data_type=dbcsr_type_complex_8,reuse_arrays=.TRUE.) ! type hard coded !
    !call dbcsr_finalize (bc_mat)
    CALL dbcsr_distribution_release (bc_dist)
    NULLIFY(col_blk_size,row_blk_size)

    CALL cp_dbcsr_get_info(bc_mat,&
         nblkrows_total=nblkrows_total,&
         nblkcols_total=nblkcols_total,&
         nblkrows_local=nblkrows_local,&
         nblkcols_local=nblkcols_local,&
         nfullrows_local=nfullrows_local,&
         nfullcols_local=nfullcols_local,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total,&
         local_rows=local_rows,&
         local_cols=local_cols,&
         row_blk_size=row_blk_size,&
         col_blk_size=col_blk_size)

    ALLOCATE (local_row_sizes (nblkrows_total))
    ALLOCATE (local_col_sizes (nblkcols_total))
    local_row_sizes(:) = 0
    IF (nblkrows_local .GE. 1) THEN
       DO row = 1 , nblkrows_local
          local_row_sizes(local_rows(row)) = row_blk_size(local_rows(row))
       END DO
    ENDIF
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       DO col = 1 , nblkcols_local
          local_col_sizes(local_cols(col)) = col_blk_size(local_cols(col))
       END DO
    ENDIF

    ALLOCATE (first_row(nblkrows_total),last_row(nblkrows_total))
    ALLOCATE (first_col(nblkcols_total),last_col(nblkcols_total))
    CALL convert_sizes_to_offsets (local_row_sizes, first_row, last_row)
    CALL convert_sizes_to_offsets (local_col_sizes, first_col, last_col)

    ! Copy the FM data to the block-cyclic DBCSR matrix.  This step
    ! could be skipped with appropriate DBCSR index manipulation.
    fm_block => fm%local_data
    CALL cp_dbcsr_work_create (bc_mat, nblks_guess=nblkrows_local*nblkcols_local,&
         sizedata_guess=nfullrows_local*nfullcols_local, work_mutable=.FALSE.,&
         n=1)
    blk_p = 1
    bc_rows: DO row_l = 1, nblkrows_local
       row = local_rows (row_l)
       row_size = row_blk_size(row)
       bc_cols: DO col_l = 1, nblkcols_local
          col = local_cols (col_l)
          col_size = col_blk_size(col)
          nze = row_size*col_size
          !WRITE(*,*)routineN//" Adding block",row,col,"size",nze
          CALL cp_add_work_coordinate(bc_mat, 1, row, col, blk_p, error=dbcsr_error)
          !blk_1d => bc_mat%m%wms(1)%data_area%d%c_dp(blk_p:blk_p+nze-1)
          blk_1d => cp_dbcsr_get_data_p (bc_mat, 1,&
               select_data_type=CMPLX(0.0, KIND=dp), lb=blk_p, ub=blk_p+nze-1)
          CALL cp_assert (nze .EQ. (last_row(row)-first_row(row)+1)*(last_col(col)-first_col(col)+1),&
               cp_fatal_level, cp_internal_error, routineN,&
               "Block size does not match block row/col sizes")
          blk_1d(1:nze) = RESHAPE(&
               fm_block(&
               first_row(row):last_row(row),first_col(col):last_col(col)&
               ), (/ nze /))
          blk_p = blk_p + nze
       ENDDO bc_cols
    ENDDO bc_rows
    CALL cp_set_work_size(bc_mat, 1, blk_p - 1)
    CALL cp_dbcsr_finalize (bc_mat, reshuffle=.FALSE.)

    ! Now convert to the desired matrix distribution
    CALL cp_dbcsr_complete_redistribute (bc_mat, matrix, keep_sparsity=keep_sparsity)
    CALL cp_dbcsr_release (bc_mat)

    CALL timestop(handle)
    CALL dbcsr_error_stop(error_handler, dbcsr_error)

  END SUBROUTINE copy_cfm_to_dbcsr
! *****************************************************************************
!> \brief   Copy a DBCSR matrix to a BLACS matrix
!> \param[in] matrix          DBCSR matrix
!> \param[out] fm             full matrix
! *****************************************************************************
  SUBROUTINE copy_dbcsr_to_fm(matrix, fm)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_fm_type), POINTER                :: fm

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_dbcsr_to_fm', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER :: col, error_handle, group, handle, mypcol, mype, myprow, &
      nblkcols_local, nblkcols_total, nblkrows_local, nblkrows_total, &
      ncol_block, ncol_global, nfullcols_total, nfullrows_total, npcol, npe, &
      nprow, nrow_block, nrow_global, row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_col, first_row, &
                                                last_col, last_row, &
                                                local_col_sizes, &
                                                local_row_sizes
    INTEGER, DIMENSION(:), POINTER :: col_blk_size, col_blk_size_tmp, &
      local_cols, local_rows, row_blk_size, row_blk_size_tmp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blk_2d, fm_block
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: fm_block_sp
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type)                      :: bc_mat
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,handle)
    CALL dbcsr_error_set (routineN, error_handle, error=dbcsr_error)

    NULLIFY(col_blk_size,row_blk_size)

    ! info about the full matrix
    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global
    ncol_global = fm%matrix_struct%ncol_global

    CALL cp_dbcsr_get_info(matrix,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total)

    ! Convert DBCSR to a block-cyclic one
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,&
         nrow_block, ncol_block,&
         dbcsr_distribution_mp (cp_dbcsr_distribution(matrix)),&
         row_blk_size, col_blk_size)

    CALL cp_dbcsr_init (bc_mat)
    row_blk_size_tmp => row_blk_size
    col_blk_size_tmp => col_blk_size
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic"//cp_dbcsr_name(matrix), bc_dist,&
         dbcsr_type_no_symmetry, row_blk_size_tmp, col_blk_size_tmp, nze=0, &
         reuse_arrays=.TRUE.)
    CALL dbcsr_distribution_release (bc_dist)
    CALL cp_dbcsr_complete_redistribute (matrix, bc_mat)

    ! Find the local extents of the local blocked rows so that index lookups
    ! into the FM matrix work correctly.
    local_rows => dbcsr_distribution_local_rows (bc_dist)
    local_cols => dbcsr_distribution_local_cols (bc_dist)
    ALLOCATE (local_row_sizes (dbcsr_distribution_nrows (bc_dist)))
    ALLOCATE (local_col_sizes (dbcsr_distribution_ncols (bc_dist)))
    nblkrows_local = dbcsr_distribution_nlocal_rows (bc_dist)
    nblkcols_local = dbcsr_distribution_nlocal_cols (bc_dist)
    local_row_sizes(:) = 0
    IF (nblkrows_local .GE. 1) THEN
       DO row = 1 , nblkrows_local
          local_row_sizes(local_rows(row)) = row_blk_size(local_rows(row))
       END DO
    ENDIF
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       DO col = 1 , nblkcols_local
          local_col_sizes(local_cols(col)) = col_blk_size(local_cols(col))
       END DO
    ENDIF
    nblkrows_total = dbcsr_distribution_nrows (bc_dist)
    nblkcols_total = dbcsr_distribution_ncols (bc_dist)
    ALLOCATE (first_row(nblkrows_total),last_row(nblkrows_total))
    ALLOCATE (first_col(nblkcols_total),last_col(nblkcols_total))
    CALL convert_sizes_to_offsets (local_row_sizes, first_row, last_row)
    CALL convert_sizes_to_offsets (local_col_sizes, first_col, last_col)
    !
    ! Now copy data to the FM matrix
    fm_block => fm%local_data
    fm_block_sp => fm%local_data_sp
    IF(fm%use_sp) THEN
       fm_block_sp=0.0_sp
    ELSE
       fm_block=0.0_dp
    ENDIF

    IF(dbg) THEN
       WRITE(*,*)routineN//" FM data size is", UBOUND(fm_block)
       WRITE(*,*)routineN//" dbcsr data size is", cp_dbcsr_get_data_size(bc_mat)
       WRITE(*,*)routineN//" FM block sizes are",nrow_block,'/',nfullrows_total
       WRITE(*,*)routineN//" FM block sizes are",ncol_block,'/',nfullcols_total
       WRITE(*,*)routineN//" dbcsr row sizes are",cp_dbcsr_row_block_sizes(bc_mat)
       WRITE(*,*)routineN//" dbcsr col sizes are",cp_dbcsr_col_block_sizes(bc_mat)
    ENDIF
    !
    CALL cp_dbcsr_iterator_start(iter, bc_mat)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block (iter, row, col, blk_2d)
       ! Convert absolute coordinates to FM-local coordinates
       IF(fm%use_sp) THEN
          fm_block_sp(first_row(row):last_row(row),first_col(col):last_col(col))&
               = REAL(blk_2d(:,:),sp)
       ELSE
          IF (dbg) THEN
             WRITE(*,*)routineN//" blk2d size",UBOUND(blk_2d)
             WRITE(*,*)routineN//" want to set coor.",row,col
             WRITE(*,*)routineN//" local extents",&
                  first_row(row),last_row(row),first_col(col),last_col(col)
          ENDIF
          fm_block(first_row(row):last_row(row),first_col(col):last_col(col))&
               = blk_2d(:,:)
       ENDIF
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL cp_dbcsr_release (bc_mat)
    CALL dbcsr_error_stop (error_handle, error=dbcsr_error)

    CALL timestop(handle)
  END SUBROUTINE copy_dbcsr_to_fm

! *****************************************************************************
!> \brief   Copy a DBCSR matrix to a BLACS matrix
!> \param[in] matrix          DBCSR matrix
!> \param[out] fm             full matrix
! *****************************************************************************
  SUBROUTINE copy_dbcsr_to_cfm(matrix, fm)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_cfm_type), POINTER               :: fm

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_dbcsr_to_cfm', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: blk_2d, fm_block
    INTEGER :: col, error_handle, group, handle, mypcol, mype, myprow, &
      nblkcols_local, nblkcols_total, nblkrows_local, nblkrows_total, &
      ncol_block, ncol_global, nfullcols_total, nfullrows_total, npcol, npe, &
      nprow, nrow_block, nrow_global, row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_col, first_row, &
                                                last_col, last_row, &
                                                local_col_sizes, &
                                                local_row_sizes
    INTEGER, DIMENSION(:), POINTER :: col_blk_size, col_blk_size_tmp, &
      local_cols, local_rows, row_blk_size, row_blk_size_tmp
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type)                      :: bc_mat
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,handle)
    CALL dbcsr_error_set (routineN, error_handle, error=dbcsr_error)

    NULLIFY(col_blk_size,row_blk_size)

    ! info about the full matrix
    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global
    ncol_global = fm%matrix_struct%ncol_global

    CALL cp_dbcsr_get_info(matrix,&
         nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total)

    ! Convert DBCSR to a block-cyclic one
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,&
         nrow_block, ncol_block,&
         dbcsr_distribution_mp (cp_dbcsr_distribution(matrix)),&
         row_blk_size, col_blk_size)

    CALL cp_dbcsr_init (bc_mat)
    row_blk_size_tmp => row_blk_size
    col_blk_size_tmp => col_blk_size
    CALL cp_dbcsr_create (bc_mat, "Block-cyclic"//cp_dbcsr_name(matrix), bc_dist,&
         dbcsr_type_no_symmetry, row_blk_size_tmp, col_blk_size_tmp, nze=0, data_type=cp_dbcsr_get_data_type(matrix),&
         reuse_arrays=.TRUE.)
    CALL dbcsr_distribution_release (bc_dist)
    CALL cp_dbcsr_complete_redistribute (matrix, bc_mat)

    ! Find the local extents of the local blocked rows so that index lookups
    ! into the FM matrix work correctly.
    local_rows => dbcsr_distribution_local_rows (bc_dist)
    local_cols => dbcsr_distribution_local_cols (bc_dist)
    ALLOCATE (local_row_sizes (dbcsr_distribution_nrows (bc_dist)))
    ALLOCATE (local_col_sizes (dbcsr_distribution_ncols (bc_dist)))
    nblkrows_local = dbcsr_distribution_nlocal_rows (bc_dist)
    nblkcols_local = dbcsr_distribution_nlocal_cols (bc_dist)
    local_row_sizes(:) = 0
    IF (nblkrows_local .GE. 1) THEN
       DO row = 1 , nblkrows_local
          local_row_sizes(local_rows(row)) = row_blk_size(local_rows(row))
       END DO
    ENDIF
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       DO col = 1 , nblkcols_local
          local_col_sizes(local_cols(col)) = col_blk_size(local_cols(col))
       END DO
    ENDIF
    nblkrows_total = dbcsr_distribution_nrows (bc_dist)
    nblkcols_total = dbcsr_distribution_ncols (bc_dist)
    ALLOCATE (first_row(nblkrows_total),last_row(nblkrows_total))
    ALLOCATE (first_col(nblkcols_total),last_col(nblkcols_total))
    CALL convert_sizes_to_offsets (local_row_sizes, first_row, last_row)
    CALL convert_sizes_to_offsets (local_col_sizes, first_col, last_col)
    !
    ! Now copy data to the FM matrix
    fm_block => fm%local_data
    fm_block=(0.0_dp,0.0_dp)

    IF(dbg) THEN
       WRITE(*,*)routineN//" FM data size is", UBOUND(fm_block)
       WRITE(*,*)routineN//" dbcsr data size is", cp_dbcsr_get_data_size(bc_mat)
       WRITE(*,*)routineN//" FM block sizes are",nrow_block,'/',nfullrows_total
       WRITE(*,*)routineN//" FM block sizes are",ncol_block,'/',nfullcols_total
       WRITE(*,*)routineN//" dbcsr row sizes are",cp_dbcsr_row_block_sizes(bc_mat)
       WRITE(*,*)routineN//" dbcsr col sizes are",cp_dbcsr_col_block_sizes(bc_mat)
    ENDIF
    !
    CALL cp_dbcsr_iterator_start(iter, bc_mat)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block (iter, row, col, blk_2d)
       ! Convert absolute coordinates to FM-local coordinates
       IF (dbg) THEN
          WRITE(*,*)routineN//" blk2d size",UBOUND(blk_2d)
          WRITE(*,*)routineN//" want to set coor.",row,col
          WRITE(*,*)routineN//" local extents",&
               first_row(row),last_row(row),first_col(col),last_col(col)
       ENDIF
       fm_block(first_row(row):last_row(row),first_col(col):last_col(col))&
            = blk_2d(:,:)
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL cp_dbcsr_release (bc_mat)
    CALL dbcsr_error_stop (error_handle, error=dbcsr_error)

    CALL timestop(handle)
  END SUBROUTINE copy_dbcsr_to_cfm
! *****************************************************************************
!> \brief   Copy a DBCSR_BLACS matrix to a BLACS matrix 
!> \param bc_mat DBCSR matrix
!> \param[out] fm             full matrix
! *****************************************************************************
  SUBROUTINE copy_dbcsr_to_fm_bc(bc_mat, fm)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: bc_mat
    TYPE(cp_fm_type), POINTER                :: fm

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_dbcsr_to_fm_bc', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER :: col, error_handle, group, handle, mypcol, mype, myprow, &
      nblkcols_local, nblkcols_total, nblkrows_local, nblkrows_total, &
      ncol_block, ncol_global, nfullcols_total, nfullrows_total, npcol, npe, &
      nprow, nrow_block, nrow_global, row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_col, first_row, &
                                                last_col, last_row, &
                                                local_col_sizes, &
                                                local_row_sizes
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, local_cols, &
                                                local_rows, row_blk_size
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blk_2d, fm_block
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: blk_2d_sp, fm_block_sp
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(dbcsr_distribution_obj)             :: bc_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,handle)
    CALL dbcsr_error_set (routineN, error_handle, error=dbcsr_error)

    NULLIFY(col_blk_size,row_blk_size)

    ! info about the full matrix
    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context
    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1); mypcol=context%mepos(2)
    nprow=context%num_pe(1); npcol=context%num_pe(2)
    nrow_block = fm%matrix_struct%nrow_block; ncol_block = fm%matrix_struct%ncol_block
    nrow_global = fm%matrix_struct%nrow_global; ncol_global = fm%matrix_struct%ncol_global
    nfullrows_total = nrow_global; nfullcols_total=ncol_global

    ! Convert DBCSR to a block-cyclic one
    CALL dbcsr_distribution_init (bc_dist)
    CALL dbcsr_create_dist_block_cyclic (bc_dist,&
         nfullrows_total, nfullcols_total,&
         nrow_block, ncol_block,&
         dbcsr_distribution_mp (cp_dbcsr_distribution(bc_mat)),&
         row_blk_size, col_blk_size)

    ! Find the local extents of the local blocked rows so that index lookups
    ! into the FM matrix work correctly.
    local_rows => dbcsr_distribution_local_rows (bc_dist)
    local_cols => dbcsr_distribution_local_cols (bc_dist)
    ALLOCATE (local_row_sizes (dbcsr_distribution_nrows (bc_dist)))
    ALLOCATE (local_col_sizes (dbcsr_distribution_ncols (bc_dist)))
    nblkrows_local = dbcsr_distribution_nlocal_rows (bc_dist)
    nblkcols_local = dbcsr_distribution_nlocal_cols (bc_dist)
    local_row_sizes(:) = 0
    IF (nblkrows_local .GE. 1) THEN
       DO row = 1 , nblkrows_local
          local_row_sizes(local_rows(row)) = row_blk_size(local_rows(row))
       END DO
    ENDIF
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       DO col = 1 , nblkcols_local
          local_col_sizes(local_cols(col)) = col_blk_size(local_cols(col))
       END DO
    ENDIF
    nblkrows_total = dbcsr_distribution_nrows (bc_dist)
    nblkcols_total = dbcsr_distribution_ncols (bc_dist)
    ALLOCATE (first_row(nblkrows_total),last_row(nblkrows_total))
    ALLOCATE (first_col(nblkcols_total),last_col(nblkcols_total))
    CALL convert_sizes_to_offsets (local_row_sizes, first_row, last_row)
    CALL convert_sizes_to_offsets (local_col_sizes, first_col, last_col)
    !
    ! Now copy data to the FM matrix
    fm_block => fm%local_data
    fm_block_sp => fm%local_data_sp
    IF(fm%use_sp) THEN
       fm_block_sp=0.0_sp
    ELSE
       fm_block=0.0_dp
    ENDIF

    IF(dbg) THEN
       WRITE(*,*)routineN//" FM data size is", UBOUND(fm_block_sp), cp_dbcsr_get_data_type(bc_mat)
       WRITE(*,*)routineN//" dbcsr data size is", cp_dbcsr_get_data_size(bc_mat)
       WRITE(*,*)routineN//" FM block sizes are",nrow_block,'/',nfullrows_total
       WRITE(*,*)routineN//" FM block sizes are",ncol_block,'/',nfullcols_total
       WRITE(*,*)routineN//" dbcsr row sizes are",cp_dbcsr_row_block_sizes(bc_mat)
       WRITE(*,*)routineN//" dbcsr col sizes are",cp_dbcsr_col_block_sizes(bc_mat)
    ENDIF
    !
    CALL cp_dbcsr_iterator_start(iter, bc_mat)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       IF(cp_dbcsr_get_data_type (bc_mat)== dbcsr_type_real_4)THEN
          CALL cp_dbcsr_iterator_next_block (iter, row, col, blk_2d_sp)
         ! Convert absolute coordinates to FM-local coordinates
          IF(fm%use_sp) THEN
            fm_block_sp(first_row(row):last_row(row),first_col(col):last_col(col))&
                 = blk_2d_sp(:,:)
          ELSE
            fm_block(first_row(row):last_row(row),first_col(col):last_col(col))&
                 =REAL(blk_2d_sp(:,:),dp)
          END IF
       ELSE
          CALL cp_dbcsr_iterator_next_block (iter, row, col, blk_2d)
          IF(fm%use_sp) THEN
            fm_block_sp(first_row(row):last_row(row),first_col(col):last_col(col))&
                 = REAL(blk_2d_sp(:,:),sp)
          ELSE
            fm_block(first_row(row):last_row(row),first_col(col):last_col(col))&
                 = blk_2d(:,:)
          END IF
       ENDIF
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)
    DEALLOCATE(row_blk_size, col_blk_size)

    CALL dbcsr_distribution_release (bc_dist)
    CALL dbcsr_error_stop (error_handle, error=dbcsr_error)

    CALL timestop(handle)
  END SUBROUTINE copy_dbcsr_to_fm_bc

! *****************************************************************************
!> \brief Utility function to copy a specially shaped fm to dbcsr_matrix
!>        The result matrix will be the matrix in dbcsr format
!>        with the row blocks sizes according to the block_sizes of the template
!>        and the col blocks sizes evenly blocked with the internal dbcsr conversion
!>        size (32 is the current default)
!> \param matrix ...
!> \param fm_in ...
!> \param template ...
! *****************************************************************************
  SUBROUTINE cp_fm_to_dbcsr_row_template(matrix, fm_in, template)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_fm_type), POINTER                :: fm_in
    TYPE(cp_dbcsr_type), INTENT(IN)          :: template

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_to_dbcsr_row_template', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k_in
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size_right_in
    TYPE(dbcsr_distribution_obj)             :: dist_right_in

    CALL cp_fm_get_info(fm_in, ncol_global=k_in)

    CALL dbcsr_create_dist_r_unrot (dist_right_in, cp_dbcsr_distribution(template), k_in, &
         col_blk_size_right_in)
    CALL cp_dbcsr_init(matrix)
    CALL cp_dbcsr_create(matrix, "D", dist_right_in, dbcsr_type_no_symmetry, &
         cp_dbcsr_row_block_sizes(template), col_blk_size_right_in,&
         nze=0, data_type=cp_dbcsr_get_data_type(template))

    CALL copy_fm_to_dbcsr(fm_in, matrix)
    DEALLOCATE(col_blk_size_right_in)
    CALL dbcsr_distribution_release(dist_right_in)

  END SUBROUTINE cp_fm_to_dbcsr_row_template

! *****************************************************************************
!> \brief Utility function to create an arbitrary shaped dbcsr matrix
!>        with the same processor grid as the template matrix
!>        both row sizes and col sizes are evenly blocked with the internal
!>        dbcsr_conversion size (32 is the current default)
!> \param matrix dbcsr matrix to be created
!> \param template template dbcsr matrix giving its mp_env
!> \param m global row size of output matrix
!> \param n global col size of output matrix
!> \param sym ...
!> \param data_type ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_m_by_n_from_template(matrix,template,m,n,sym,data_type)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix, template
    INTEGER                                  :: m, n
    CHARACTER, OPTIONAL                      :: sym
    INTEGER, OPTIONAL                        :: data_type

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'cp_dbcsr_m_by_n_from_template', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: mysym
    INTEGER                                  :: my_data_type
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, &
                                                col_cluster_dist, col_dist, &
                                                row_blk_size, &
                                                row_cluster_dist, row_dist
    TYPE(dbcsr_distribution_obj)             :: dist_m_n

    mysym=cp_dbcsr_get_matrix_type(template)
    IF(PRESENT(sym))mysym=sym
    my_data_type=cp_dbcsr_get_data_type(template)
    IF(PRESENT(data_type))my_data_type=data_type

    NULLIFY(row_dist,col_dist)
    NULLIFY(row_blk_size,col_blk_size)
    NULLIFY(row_cluster_dist, col_cluster_dist)

    CALL cp_create_bl_distribution (row_dist, row_blk_size, m, &
         dbcsr_mp_nprows(dbcsr_distribution_mp(cp_dbcsr_distribution(template))))
    CALL cp_create_bl_distribution (col_dist, col_blk_size, n, &
         dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(template))))
    CALL cp_dbcsr_distribution_new (dist_m_n, dbcsr_distribution_mp (cp_dbcsr_distribution(template)),&
         row_dist, col_dist, &
         row_cluster_dist, col_cluster_dist, &
         reuse_arrays=.TRUE.)

    CALL cp_dbcsr_create(matrix, "m_n_template", dist_m_n, mysym,&
            row_blk_size, col_blk_size, nze=0, data_type=my_data_type,&
            reuse_arrays=.TRUE.)
    CALL dbcsr_distribution_release(dist_m_n)

  END SUBROUTINE cp_dbcsr_m_by_n_from_template

! *****************************************************************************
!> \brief Utility function to create dbcsr matrix, m x n matrix (n arbitrary)
!>        with the same processor grid and row distribution  as the template matrix
!>        col sizes are evenly blocked with the internal
!>        dbcsr_conversion size (32 is the current default)
!> \param matrix dbcsr matrix to be created
!> \param template template dbcsr matrix giving its mp_env
!> \param n global col size of output matrix
!> \param sym ...
!> \param data_type ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_m_by_n_from_row_template(matrix,template,n,sym,data_type)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix, template
    INTEGER                                  :: n
    CHARACTER, OPTIONAL                      :: sym
    INTEGER, OPTIONAL                        :: data_type

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'cp_dbcsr_m_by_n_from_row_template', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: mysym
    INTEGER                                  :: my_data_type
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, &
                                                col_cluster_dist, col_dist
    TYPE(dbcsr_distribution_obj)             :: dist_m_n

    mysym=cp_dbcsr_get_matrix_type(template)
    IF(PRESENT(sym))mysym=sym
    my_data_type=cp_dbcsr_get_data_type(template)
    IF(PRESENT(data_type))my_data_type=data_type

    NULLIFY(col_dist,col_blk_size)

    CALL cp_create_bl_distribution (col_dist, col_blk_size, n, &
         dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(template))))
    NULLIFY(col_cluster_dist)
    CALL cp_dbcsr_distribution_new (dist_m_n, dbcsr_distribution_mp (cp_dbcsr_distribution(template)),&
         dbcsr_distribution_row_dist(cp_dbcsr_distribution(template)), col_dist, &
         dbcsr_distribution_row_clusters(cp_dbcsr_distribution(template)), &
         col_cluster_dist)

    CALL cp_dbcsr_create(matrix, "m_n_template", dist_m_n, mysym,&
            cp_dbcsr_row_block_sizes(template), col_blk_size, nze=0, data_type=my_data_type)

    DEALLOCATE(col_dist,col_blk_size)
    CALL dbcsr_distribution_release(dist_m_n)

  END SUBROUTINE cp_dbcsr_m_by_n_from_row_template  
    


END MODULE cp_dbcsr_operations
