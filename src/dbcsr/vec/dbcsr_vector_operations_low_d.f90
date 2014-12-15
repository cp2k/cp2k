! *****************************************************************************
!> \brief the real driver routine for the multiply, not all symmetries implemented yet
!> \param matrix ...
!> \param vec_in ...
!> \param vec_out ...
!> \param alpha ...
!> \param beta ...
!> \param work_row ...
!> \param work_col ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_matrix_colvec_multiply_low_d(matrix, vec_in, vec_out, alpha, beta, work_row, work_col, error)
    TYPE(dbcsr_obj)                          :: matrix, vec_in, vec_out
    REAL(kind=real_8)                          :: alpha, beta
    TYPE(dbcsr_obj)                          :: work_row, work_col
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_matrix_colvec_multiply_low', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: matrix_type

    matrix_type=dbcsr_get_matrix_type(matrix)
    SELECT CASE(matrix_type)
    CASE(dbcsr_type_no_symmetry)
       CALL dbcsr_matrix_vector_mult_d(matrix, vec_in, vec_out, alpha, beta, work_row, work_col, error)
    CASE(dbcsr_type_symmetric)
       CALL dbcsr_sym_matrix_vector_mult_d(matrix, vec_in, vec_out, alpha, beta, work_row, work_col, error)
    CASE(dbcsr_type_antisymmetric)
        ! Not yet implemented, should mainly be some prefactor magic, but who knows how antisymmetric matrices are stored???
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error, &
            routineN, "NYI, antisymmetric matrix not permitted", __LINE__, error)
    CASE DEFAULT
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error, &
            routineN, "Unknown matrix type, ...", __LINE__, error)
    END SELECT

  END SUBROUTINE dbcsr_matrix_colvec_multiply_low_d

! *****************************************************************************
!> \brief low level routines for matrix vector multiplies
!> \param matrix ...
!> \param vec_in ...
!> \param vec_out ...
!> \param alpha ...
!> \param beta ...
!> \param work_row ...
!> \param work_col ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_matrix_vector_mult_d(matrix, vec_in, vec_out, alpha, beta, work_row, work_col, error)
    TYPE(dbcsr_obj)                          :: matrix, vec_in, vec_out
    REAL(kind=real_8)                          :: alpha, beta
    TYPE(dbcsr_obj)                          :: work_row, work_col
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_matrix_vector_mult', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, col_size, mypcol, &
                                                myprow, nblk_col, nblk_row, &
                                                ncols, pcol_group, nrows, &
                                                prow_group, row, row_size, &
                                                handle, handle1, ithread 
    LOGICAL                                  :: transposed
    REAL(kind=real_8), DIMENSION(:), POINTER          :: data_vec
    REAL(kind=real_8), DIMENSION(:, :), POINTER       :: data_d, vec_res
    TYPE(block_ptr_d), ALLOCATABLE, &
      DIMENSION(:)                           :: blk_map_col, blk_map_row
    TYPE(dbcsr_distribution_obj)             :: distri
    TYPE(dbcsr_iterator)                     :: iter

    CALL dbcsr_error_set(routineN, handle, error)
    ithread=0

! Collect some data about the parallel environment. We will use them later to move the vector around
    CALL dbcsr_get_info(matrix=matrix, distribution=distri)
    prow_group=distri%d%mp_env%mp%prow_group 
    pcol_group=distri%d%mp_env%mp%pcol_group 
    mypcol=distri%d%mp_env%mp%mypcol
    myprow=distri%d%mp_env%mp%myprow

! Create pointers to the content of the vectors, this is simply done for efficiency in the later steps
! However, it restricts the sparsity of the vector as it assumes the sparsity pattern does not change
    CALL dbcsr_get_info(matrix=work_row, nblkcols_total=nblk_col)
    CALL dbcsr_get_info(matrix=work_col, nblkrows_total=nblk_row, nfullcols_local=ncols)
    ALLOCATE(blk_map_row(nblk_col));  ALLOCATE(blk_map_col(nblk_row))
    CALL assign_row_vec_block_ptr_d(work_row, blk_map_row, error)
    CALL assign_col_vec_block_ptr_d(work_col, blk_map_col, error)

! Transfer the correct parts of the input vector to the correct locations so we can do a local multiply
    CALL dbcsr_col_vec_to_rep_row_d(vec_in, work_col, work_row, blk_map_col, error)

! Set the work vector for the results to 0
    CALL dbcsr_set(work_col, 0.0_real_8, error=error)

! Perform the local multiply. Here we exploit, that we have the blocks replicated on the mpi processes
! It is important to note, that the input and result vector are sitributed differently (row wise, col wise respectively)
    CALL dbcsr_error_set(routineN//"_local_mm", handle1, error)

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(row,col,iter,data_d,row_size,col_size,transposed,ithread) &
!$OMP          SHARED(matrix,blk_map_row,blk_map_col,ncols)
    !$ ithread = omp_get_thread_num ()
    CALL dbcsr_iterator_start(iter, matrix, shared=.FALSE.)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_d, transposed, row_size=row_size, col_size=col_size)
       IF(blk_map_col(row)%assigned_thread .NE. ithread ) CYCLE
       blk_map_col(row)%ptr=blk_map_col(row)%ptr+MATMUL(data_d,TRANSPOSE(blk_map_row(col)%ptr))
    END DO
    CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

    CALL dbcsr_error_stop(handle1, error)

! sum all the data onto the first processor col where the original vector is stored
    data_vec => dbcsr_get_data_p (work_col%m%data_area, coersion=0.0_real_8)
    CALL dbcsr_get_info(matrix=work_col, nfullrows_local=nrows, nfullcols_local=ncols)
    CALL mp_sum(data_vec(1:nrows*ncols), prow_group)

! Local copy on the first mpi col (as this is the localtion of the vec_res blocks) of the result vector
! from the replicated to the original vector. Let's play it safe and use the iterator
    CALL dbcsr_iterator_start(iter, vec_out)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_res, transposed, row_size=row_size)
       IF(ASSOCIATED(blk_map_col(row)%ptr))THEN
          vec_res(:, :)= beta*vec_res(:, :)+alpha*blk_map_col(row)%ptr(:, :)     
       ELSE
          vec_res(:, :)= beta*vec_res(:, :)
       END IF
    END DO
    CALL dbcsr_iterator_stop(iter)
    DEALLOCATE(blk_map_col, blk_map_row)

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE dbcsr_matrix_vector_mult_d

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param vec_in ...
!> \param vec_out ...
!> \param alpha ...
!> \param beta ...
!> \param work_row ...
!> \param work_col ...
!> \param skip_diag ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_matrixT_vector_mult_d(matrix, vec_in, vec_out, alpha, beta, work_row, work_col, skip_diag, error)
    TYPE(dbcsr_obj)                          :: matrix, vec_in, vec_out
    REAL(kind=real_8)                          :: alpha, beta
    TYPE(dbcsr_obj)                          :: work_row, work_col
    LOGICAL                                  :: skip_diag
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_matrixT_vector_mult', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, col_size, mypcol, &
                                                myprow, nblk_col, nblk_row, &
                                                ncols, pcol_group, nrows, &
                                                prow_group, row, row_size, &
                                                handle, handle1, ithread
    LOGICAL                                  :: transposed
    REAL(kind=real_8), DIMENSION(:), POINTER          :: data_vec
    REAL(kind=real_8), DIMENSION(:, :), POINTER       :: data_d, vec_bl, vec_res
    TYPE(block_ptr_d), ALLOCATABLE, &
      DIMENSION(:)                           :: blk_map_col, blk_map_row
    TYPE(dbcsr_distribution_obj)             :: distri
    TYPE(dbcsr_iterator)                     :: iter

    CALL dbcsr_error_set(routineN, handle, error)
    ithread=0

! Collect some data about the parallel environment. We will use them later to move the vector around
    CALL dbcsr_get_info(matrix=matrix, distribution=distri)
    prow_group=distri%d%mp_env%mp%prow_group; pcol_group=distri%d%mp_env%mp%pcol_group
    mypcol=distri%d%mp_env%mp%mypcol; myprow=distri%d%mp_env%mp%myprow

! Perform the local multiply. Here we exploit, that we have the blocks replicated on the mpi processes
! It is important to note, that the input and result vector are sitributed differently (row wise, col wise respectively)
    CALL dbcsr_get_info(matrix=work_row, nblkcols_total=nblk_col)
    CALL dbcsr_get_info(matrix=work_col, nblkrows_total=nblk_row)
    ALLOCATE(blk_map_row(nblk_col));  ALLOCATE(blk_map_col(nblk_row))
    CALL assign_row_vec_block_ptr_d(work_row, blk_map_row, error)
    CALL assign_col_vec_block_ptr_d(work_col, blk_map_col, error)

! Set the work vector for the results to 0
    CALL dbcsr_set(work_row, 0.0_real_8, error=error)

! Transfer the correct parts of the input vector to the replicated vector on proc_col 0
    CALL dbcsr_iterator_start(iter, vec_in)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_bl, transposed, row_size=row_size, col_size=col_size)
       blk_map_col(row)%ptr(1:row_size, 1:col_size)= vec_bl(1:row_size, 1:col_size)
    END DO
    CALL dbcsr_iterator_stop(iter)
! Replicate the data on all processore in the row
    data_vec => dbcsr_get_data_p (work_col%m%data_area, coersion=0.0_real_8)
    CALL mp_bcast(data_vec, 0, prow_group)

! Perform the local multiply. Here it is obvious why the vectors are replicated on the mpi rows and cols
    CALL dbcsr_error_set(routineN//"local_mm", handle1, error)
    CALL dbcsr_get_info(matrix=work_col, nfullcols_local=ncols)
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(row,col,iter,data_d,row_size,col_size,transposed,ithread) &
!$OMP          SHARED(matrix,blk_map_row,blk_map_col,skip_diag,ncols)
    !$ ithread = omp_get_thread_num ()
    CALL dbcsr_iterator_start(iter, matrix, shared=.FALSE.)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_d, transposed, row_size=row_size, col_size=col_size)
       IF(skip_diag.AND.col==row)CYCLE
       IF(ASSOCIATED(blk_map_row(col)%ptr).AND.ASSOCIATED(blk_map_col(row)%ptr))THEN
          IF(blk_map_row(col)%assigned_thread .NE. ithread ) CYCLE
          blk_map_row(col)%ptr=blk_map_row(col)%ptr+MATMUL(TRANSPOSE(blk_map_col(row)%ptr),data_d)
       ELSE
          IF(blk_map_row(row)%assigned_thread .NE. ithread ) CYCLE
          blk_map_row(row)%ptr=blk_map_row(row)%ptr+MATMUL(TRANSPOSE(blk_map_col(col)%ptr),TRANSPOSE(data_d))
       END IF
    END DO
    CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

    CALL dbcsr_error_stop(handle1, error)

! sum all the data within a processor column to obtain the replicated result
    data_vec => dbcsr_get_data_p (work_row%m%data_area, coersion=0.0_real_8)
! we use the replicated vector but the final answer is only summed to proc_col 0 for effiniency
    CALL dbcsr_get_info(matrix=work_row, nfullrows_local=nrows, nfullcols_local=ncols)
    CALL mp_sum(data_vec(1:nrows*ncols), pcol_group)

! Convert the result to a column wise distribution
    CALL dbcsr_rep_row_to_rep_col_vec_d(work_col, work_row, blk_map_row, error=error)
    
! Create_the final vector by summing it to the result vector which lives on proc_col 0
    CALL dbcsr_iterator_start(iter, vec_out)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_res, transposed, row_size=row_size)
       IF(ASSOCIATED(blk_map_col(row)%ptr))THEN
          vec_res(:, :)= beta*vec_res(:, :)+alpha*blk_map_col(row)%ptr(:, :)  
       ELSE
          vec_res(:, :)= beta*vec_res(:, :)
       END IF
    END DO
    CALL dbcsr_iterator_stop(iter)
    DEALLOCATE(blk_map_col, blk_map_row)

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE dbcsr_matrixT_vector_mult_d 

! *****************************************************************************
!> \brief ...
!> \param vec_in ...
!> \param rep_col_vec ...
!> \param rep_row_vec ...
!> \param blk_map_col ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_col_vec_to_rep_row_d(vec_in, rep_col_vec, rep_row_vec, blk_map_col, error)
    TYPE(dbcsr_obj)                          :: vec_in, rep_col_vec, &
                                                rep_row_vec
    TYPE(block_ptr_d), ALLOCATABLE, &
      DIMENSION(:)                           :: blk_map_col
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_col_vec_to_rep_row', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, mypcol, myprow, ncols, &
                                                nrows, pcol_group, &
                                                prow_group, row, handle
    INTEGER, DIMENSION(:), POINTER           :: local_cols, row_dist
    LOGICAL                                  :: transposed
    REAL(kind=real_8), DIMENSION(:), POINTER          :: data_vec, data_vec_rep
    REAL(kind=real_8), DIMENSION(:, :), POINTER       :: vec_row
    TYPE(dbcsr_distribution_obj)             :: distri
    TYPE(dbcsr_iterator)                     :: iter

    CALL dbcsr_error_set(routineN, handle, error)

! get information about the parallel environment
    CALL dbcsr_get_info(matrix=vec_in, distribution=distri)
    prow_group=distri%d%mp_env%mp%prow_group
    pcol_group=distri%d%mp_env%mp%pcol_group
    mypcol=distri%d%mp_env%mp%mypcol
    myprow=distri%d%mp_env%mp%myprow

! Get the vector which tells us which blocks are local to which processor row in the col vec
    row_dist=>array_data(dbcsr_distribution_row_dist (dbcsr_distribution(rep_col_vec)))

! Copy the local vector to the replicated on the first processor column (this is where vec_in lives)
    CALL dbcsr_get_info(matrix=rep_col_vec, nfullrows_local=nrows, nfullcols_local=ncols)
    data_vec_rep => dbcsr_get_data_p (rep_col_vec%m%data_area, coersion=0.0_real_8)
    data_vec => dbcsr_get_data_p (vec_in%m%data_area, coersion=0.0_real_8)
    IF(mypcol==0)data_vec_rep(1:nrows*ncols)=data_vec(1:nrows*ncols)
! Replicate the data along the row
    CALL mp_bcast(data_vec_rep(1:nrows*ncols), 0, prow_group)

! Here it gets a bit tricky as we are dealing with two different parallel layouts:
! The rep_col_vec contains all blocks local to the row distribution of the vector. 
! The rep_row_vec only needs the fraction which is local to the col distribution.
! However in most cases this won't the complete set of block which can be obtained from col_vector p_row i
! Anyway, as the blocks don't repeat in the col_vec, a different fraction of the row vec will be available
! on every replica in the processor column, by summing along the column we end up with the complete vector everywhere
! Hope this clarifies the idea
    CALL dbcsr_set(rep_row_vec, 0.0_real_8, error=error)
    CALL dbcsr_get_info(matrix=rep_row_vec, nfullrows_local=nrows, local_cols=local_cols, nfullcols_local=ncols)
    CALL dbcsr_iterator_start(iter, rep_row_vec)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_row, transposed)
       IF(row_dist(col)==myprow)THEN
          vec_row=TRANSPOSE(blk_map_col(col)%ptr)
       END IF
    END DO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_get_info(matrix=rep_row_vec, nfullrows_local=nrows, nfullcols_local=ncols)
    data_vec_rep => dbcsr_get_data_p (rep_row_vec%m%data_area, coersion=0.0_real_8)
    CALL mp_sum(data_vec_rep(1:ncols*nrows), pcol_group)

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE dbcsr_col_vec_to_rep_row_d    
       
! *****************************************************************************
!> \brief ...
!> \param rep_col_vec ...
!> \param rep_row_vec ...
!> \param blk_map_row ...
!> \param blk_map_col_add ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_rep_row_to_rep_col_vec_d(rep_col_vec, rep_row_vec, blk_map_row, blk_map_col_add, error)
    TYPE(dbcsr_obj)                          :: rep_col_vec, rep_row_vec
    TYPE(block_ptr_d),DIMENSION(:), OPTIONAL     :: blk_map_col_add
    TYPE(block_ptr_d), ALLOCATABLE, &
      DIMENSION(:)                           :: blk_map_row
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_rep_row_to_rep_col_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, mypcol, myprow, ncols, &
                                                nrows, pcol_group, &
                                                prow_group, row, handle
    INTEGER, DIMENSION(:), POINTER           :: col_dist
    LOGICAL                                  :: transposed
    REAL(kind=real_8), DIMENSION(:), POINTER          :: data_vec_rep
    REAL(kind=real_8), DIMENSION(:, :), POINTER       :: vec_col
    TYPE(dbcsr_distribution_obj)             :: distri
    TYPE(dbcsr_iterator)                     :: iter

    CALL dbcsr_error_set(routineN, handle, error)

! get information about the parallel environment
    CALL dbcsr_get_info(matrix=rep_col_vec, distribution=distri)
    prow_group=distri%d%mp_env%mp%prow_group
    pcol_group=distri%d%mp_env%mp%pcol_group
    mypcol=distri%d%mp_env%mp%mypcol
    myprow=distri%d%mp_env%mp%myprow
! Get the vector which tells us which blocks are local to which processor col in the row vec
    col_dist=>array_data(dbcsr_distribution_col_dist (dbcsr_distribution(rep_row_vec)))

! The same trick as described above with opposite direction
    CALL dbcsr_set(rep_col_vec, 0.0_real_8, error=error)
    CALL dbcsr_iterator_start(iter, rep_col_vec)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_col, transposed)
       IF(col_dist(row)==mypcol)THEN
          vec_col=TRANSPOSE(blk_map_row(row)%ptr)
       END IF
       ! this one is special and allows to add the elements of a not yet summed replicated 
       ! column vector as it appears in M*V(row_rep) as result. Save an mp_sum in the symmetric case
       IF(PRESENT(blk_map_col_add))vec_col=vec_col+blk_map_col_add(row)%ptr(:,:)
    END DO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_get_info(matrix=rep_col_vec, nfullrows_local=nrows, nfullcols_local=ncols)
    data_vec_rep => dbcsr_get_data_p (rep_col_vec%m%data_area, coersion=0.0_real_8)
    CALL mp_sum(data_vec_rep(1:nrows*ncols), prow_group)

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE dbcsr_rep_row_to_rep_col_vec_d

! *****************************************************************************
!> \brief ...
!> \param row_vec ...
!> \param row_blk_ptr ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE assign_row_vec_block_ptr_d(row_vec, row_blk_ptr, error)  
    TYPE(dbcsr_obj)                          :: row_vec
    TYPE(block_ptr_d), DIMENSION(:)            :: row_blk_ptr
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'assign_row_vec_block_ptr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, row, handle, iblock, nthreads
    LOGICAL                                  :: transposed
    REAL(kind=real_8), DIMENSION(:, :), POINTER       :: vec_bl
    TYPE(dbcsr_iterator)                     :: iter

    nthreads = 1
!$OMP PARALLEL DEFAULT(NONE) SHARED(nthreads)
!$OMP MASTER
    !$ nthreads = OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

    CALL dbcsr_error_set(routineN, handle, error)

    iblock=0
    CALL dbcsr_iterator_start(iter, row_vec)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_bl, transposed)
       iblock=iblock+1
       row_blk_ptr(col)%ptr=>vec_bl
       row_blk_ptr(col)%assigned_thread=MOD(iblock,nthreads)
    END DO
    CALL dbcsr_iterator_stop(iter)

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE assign_row_vec_block_ptr_d

! *****************************************************************************
!> \brief ...
!> \param col_vec ...
!> \param col_blk_ptr ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE assign_col_vec_block_ptr_d(col_vec, col_blk_ptr, error)
    TYPE(dbcsr_obj)                          :: col_vec
    TYPE(block_ptr_d), DIMENSION(:)            :: col_blk_ptr
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'assign_col_vec_block_ptr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, row, handle, iblock, nthreads
    LOGICAL                                  :: transposed
    REAL(kind=real_8), DIMENSION(:, :), POINTER       :: vec_bl
    TYPE(dbcsr_iterator)                     :: iter

    CALL dbcsr_error_set(routineN, handle, error)

    nthreads = 1
!$OMP PARALLEL DEFAULT(NONE) SHARED(nthreads)
!$OMP MASTER
    !$ nthreads = OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

    iblock=0
    CALL dbcsr_iterator_start(iter, col_vec)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_bl, transposed)
       iblock=iblock+1
       col_blk_ptr(row)%ptr=>vec_bl
       col_blk_ptr(row)%assigned_thread=MOD(iblock,nthreads)
    END DO
    CALL dbcsr_iterator_stop(iter)

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE assign_col_vec_block_ptr_d

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param vec_in ...
!> \param vec_out ...
!> \param alpha ...
!> \param beta ...
!> \param work_row ...
!> \param work_col ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_sym_matrix_vector_mult_d(matrix, vec_in, vec_out, alpha, beta, work_row, work_col, error)
    TYPE(dbcsr_obj)                          :: matrix, vec_in, vec_out
    REAL(kind=real_8)                          :: alpha, beta
    TYPE(dbcsr_obj)                          :: work_row, work_col
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_sym_m_v_mult', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, col_size, mypcol, &
                                                myprow, nblk_col, nblk_row, &
                                                pcol_group, nrows, ncols,&
                                                prow_group, row, row_size, &
                                                handle, handle1, ithread, vec_dim
    LOGICAL                                  :: transposed
    REAL(kind=real_8), DIMENSION(:), POINTER          :: data_vec
    REAL(kind=real_8), DIMENSION(:, :), POINTER       :: data_d, vec_res
    TYPE(block_ptr_d), ALLOCATABLE, &
      DIMENSION(:)                           :: blk_map_col, blk_map_row, res_blk_map_col, res_blk_map_row
    TYPE(dbcsr_distribution_obj)             :: distri
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_obj)                          :: result_row, result_col

    CALL dbcsr_error_set(routineN, handle, error)
    ithread=0
! We need some work matrices as we try to exploit operations on the replicated vectors which are duplicated otherwise
    CALL dbcsr_get_info(matrix=vec_in,nfullcols_total=vec_dim)
! This is a performance hack as the new creation of a replicated vector is a fair bit more expensive
    CALL dbcsr_init(result_col)
    CALL dbcsr_set(work_col, 0.0_real_8, error=error)
    CALL dbcsr_copy(result_col, work_col, error=error)
    CALL dbcsr_init(result_row)
    CALL dbcsr_set(work_row, 0.0_real_8, error=error)
    CALL dbcsr_copy(result_row, work_row, error=error)

! Collect some data about the parallel environment. We will use them later to move the vector around
    CALL dbcsr_get_info(matrix=matrix, distribution=distri)
    prow_group=distri%d%mp_env%mp%prow_group; pcol_group=distri%d%mp_env%mp%pcol_group
    mypcol=distri%d%mp_env%mp%mypcol; myprow=distri%d%mp_env%mp%myprow

! Perform the local multiply. Here we exploit, that we have the blocks replicated on the mpi processes
! It is important to note, that the input and result vector are sitributed differently (row wise, col wise respectively)
    CALL dbcsr_get_info(matrix=work_row, nblkcols_total=nblk_col)
    CALL dbcsr_get_info(matrix=work_col, nblkrows_total=nblk_row)
    ALLOCATE(blk_map_row(nblk_col));  ALLOCATE(blk_map_col(nblk_row))
    ALLOCATE(res_blk_map_row(nblk_col));  ALLOCATE(res_blk_map_col(nblk_row))
    CALL assign_row_vec_block_ptr_d(work_row, blk_map_row, error)
    CALL assign_col_vec_block_ptr_d(work_col, blk_map_col, error)
    CALL assign_row_vec_block_ptr_d(result_row, res_blk_map_row, error)
    CALL assign_col_vec_block_ptr_d(result_col, res_blk_map_col, error)

! Transfer the correct parts of the input vector to the correct locations so we can do a local multiply
    CALL dbcsr_col_vec_to_rep_row_d(vec_in, work_col, work_row, blk_map_col, error)

! Probably I should rename the routine above as it delivers both the replicated row and column vector

! Perform the local multiply. Here we exploit, that we have the blocks replicated on the mpi processes
! It is important to note, that the input and result vector are distributed differently (row wise, col wise respectively)
    CALL dbcsr_error_set(routineN//"_local_mm", handle1, error)

!------ perform the multiplication, we have to take car to take the correct blocks ----------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(row,col,iter,data_d,row_size,col_size,transposed,ithread) &
!$OMP          SHARED(matrix,blk_map_row,res_blk_map_col,res_blk_map_row,blk_map_col)
    !$ ithread = omp_get_thread_num ()
    CALL dbcsr_iterator_start(iter, matrix, shared=.FALSE.)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_d, transposed, row_size=row_size, col_size=col_size)
       IF(ASSOCIATED(blk_map_row(col)%ptr).AND.ASSOCIATED(res_blk_map_col(row)%ptr))THEN
          IF(res_blk_map_col(row)%assigned_thread .EQ. ithread ) THEN
             res_blk_map_col(row)%ptr=res_blk_map_col(row)%ptr+MATMUL(data_d,TRANSPOSE(blk_map_row(col)%ptr))
          END IF
          IF(res_blk_map_row(col)%assigned_thread .EQ. ithread .AND. row .NE. col) THEN
             res_blk_map_row(col)%ptr=res_blk_map_row(col)%ptr+MATMUL(TRANSPOSE(blk_map_col(row)%ptr),data_d)
          END IF
       ELSE
          IF(res_blk_map_col(col)%assigned_thread .EQ. ithread ) THEN
             res_blk_map_col(col)%ptr=res_blk_map_col(col)%ptr+TRANSPOSE(MATMUL(blk_map_row(row)%ptr,data_d))
          END IF
          IF(res_blk_map_row(row)%assigned_thread .EQ. ithread  .AND. row .NE. col ) THEN
             res_blk_map_row(row)%ptr=res_blk_map_row(row)%ptr+TRANSPOSE(MATMUL(data_d,blk_map_col(col)%ptr))
          END IF
       END IF 
    END DO
    CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

    CALL dbcsr_error_stop(handle1, error)

    ! sum all the data within a processor column to obtain the replicated result from lower
    data_vec => dbcsr_get_data_p (result_row%m%data_area, coersion=0.0_real_8)
    CALL dbcsr_get_info(matrix=result_row, nfullrows_local=nrows, nfullcols_local=ncols)

    CALL mp_sum(data_vec(1:nrows*ncols), pcol_group)
!
!! Convert the results to a column wise distribution, this is a bit involved as the result_row is fully replicated
!! While the result_col still has the partial results in parallel. The routine below takes care of that and saves an
!! mp_sum. Of the res_row vectors are created only taking the approriate element (0 otherwise) while the res_col
!! parallel bits are locally added. The mp_sum magically creates the correct vector
    CALL dbcsr_rep_row_to_rep_col_vec_d(work_col, result_row, res_blk_map_row, res_blk_map_col, error)

!    ! Create_the final vector by summing it to the result vector which lives on proc_col 0 lower
    CALL dbcsr_iterator_start(iter, vec_out)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_res, transposed, row_size=row_size)
       IF(ASSOCIATED(blk_map_col(row)%ptr))THEN
          vec_res(:, :)= beta*vec_res(:, :)+alpha*(blk_map_col(row)%ptr(:, :)) 
       ELSE
          vec_res(:, :)= beta*vec_res(:, :)
       END IF
    END DO
    CALL dbcsr_iterator_stop(iter)

    DEALLOCATE(blk_map_col, blk_map_row, res_blk_map_row, res_blk_map_col)

    CALL dbcsr_release(result_row); CALL dbcsr_release(result_col)

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE dbcsr_sym_matrix_vector_mult_d

