! **************************************************************************************************
!> \brief Calculates norms with minimal overhead.
!> \param buffer ...
!> \param norms ...
! **************************************************************************************************
  SUBROUTINE calc_norms_z(matrix, norms)
  TYPE(dbcsr_type), INTENT(IN)               :: matrix
  REAL(kind=sp), DIMENSION(:), INTENT(INOUT) :: norms

  INTEGER, DIMENSION(:), POINTER    :: row, col, bps, rbs, cbs, &
                                       local_rows, local_cols
  COMPLEX(kind=real_8), DIMENSION(:), POINTER    :: data
  INTEGER                           :: blk, bpe

  row => matrix%coo_l(1::3)
  col => matrix%coo_l(2::3)
  bps => matrix%coo_l(3::3)
  rbs => array_data(matrix%row_blk_size)
  cbs => array_data(matrix%col_blk_size)
  local_rows => array_data(matrix%local_rows)
  local_cols => array_data(matrix%local_cols)
  data => dbcsr_get_data_p_z (matrix%data_area)
  !$omp parallel default(none) &
  !$omp shared(row,col,bps,data,&
  !$omp        rbs,cbs,local_rows,local_cols, &
  !$omp        norms,matrix) &
  !$omp private(blk,bpe)
  !$omp do
  DO blk = 1, matrix%nblks
     bpe = bps(blk) + rbs(local_rows(row(blk))) * cbs(local_cols(col(blk))) - 1
     norms(blk) = SQRT (REAL (SUM(ABS(data(bps(blk):bpe))**2), KIND=sp))
  ENDDO
  !$omp end do
  !$omp end parallel

END SUBROUTINE calc_norms_z

