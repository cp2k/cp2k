!> \brief Calculates norms of the entire matrix with minimal overhead.
!> \param norms ...
!> \param nrows ...
!> \param row_p ...
!> \param col_i ...
!> \param blk_p ...
!> \param rbs ...
!> \param cbs ...
!> \param DATA ...
!> \param local ...
!> \param local2global ...
!> \param max_val ...
! *****************************************************************************
  SUBROUTINE calc_norms_z(norms, nrows,&
       row_p, col_i, blk_p, rbs, cbs, DATA, local, local2global,&
       max_val)
    REAL(kind=sp), DIMENSION(:), INTENT(OUT) :: norms
    INTEGER, INTENT(IN)                      :: nrows
    INTEGER, DIMENSION(1:nrows+1), &
      INTENT(IN)                             :: row_p
    INTEGER, DIMENSION(*), INTENT(IN)        :: col_i, blk_p, rbs, cbs
    COMPLEX(kind=real_8), DIMENSION(*), &
      INTENT(IN)                             :: DATA
    LOGICAL, INTENT(IN)                      :: local
    INTEGER, DIMENSION(*), INTENT(IN)        :: local2global
    REAL(kind=sp), INTENT(OUT)               :: max_val

    INTEGER                                  :: blk, bp, bpe, row, row_i, &
                                                row_size
    LOGICAL                                  :: is_in_parallel

!   ---------------------------------------------------------------------------

    max_val = 0
    is_in_parallel = .FALSE.
!$  is_in_parallel = omp_in_parallel()

    !$omp parallel default(none) &
    !$omp          private (row_i, row, row_size, blk, bp, bpe) &
    !$omp          shared (nrows, local) &
    !$omp          shared (local2global, rbs, cbs, row_p, col_i, blk_p, &
    !$omp                  data, norms) &
    !$omp          reduction (max:max_val) if (.not.is_in_parallel)
    !$omp do
    DO row_i = 1, nrows
       IF (local) THEN
          row = local2global(row_i)
       ELSE
          row = row_i
       ENDIF
       row_size = rbs(row)
       DO blk = row_p(row_i)+1, row_p(row_i+1)
          IF (blk_p(blk) .NE. 0) THEN
             bp = ABS(blk_p(blk))
             bpe = bp + row_size * cbs(col_i(blk)) - 1
             norms(blk) = SQRT (REAL (SUM(ABS(DATA(bp:bpe))**2), KIND=sp))
          ELSE
             norms(blk) = 0.0_sp
          ENDIF
          max_val = MAX(max_val,norms(blk))
       ENDDO
    ENDDO
    !$omp end do
    !$omp end parallel
  END SUBROUTINE calc_norms_z

!> \brief Calculates norms of the entire matrix with minimal overhead.
!> \param norms ...
!> \param nblks ...
!> \param blki ...
!> \param rbs ...
!> \param cbs ...
!> \param DATA ...
!> \param local ...
!> \param local2global_rows ...
!> \param local2global_cols ...
!> \param max_val ...
! *****************************************************************************
  SUBROUTINE calc_norms_list_z(norms, nblks,&
       blki, rbs, cbs, DATA, local, local2global_rows, local2global_cols,&
       max_val)
    REAL(kind=sp), DIMENSION(:), INTENT(OUT) :: norms
    INTEGER, INTENT(IN)                      :: nblks
    INTEGER, DIMENSION(3,nblks), INTENT(IN)  :: blki
    INTEGER, DIMENSION(:), INTENT(IN)        :: rbs, cbs
    COMPLEX(kind=real_8), DIMENSION(:), &
      INTENT(IN)                             :: DATA
    LOGICAL, INTENT(IN)                      :: local
    INTEGER, DIMENSION(:), INTENT(IN)        :: local2global_rows
    INTEGER, DIMENSION(:), INTENT(IN)        :: local2global_cols
    REAL(kind=sp), INTENT(OUT)               :: max_val

    INTEGER                                  :: blk, bp, bpe, row, col
    LOGICAL                                  :: is_in_parallel

!   ---------------------------------------------------------------------------

    max_val = 0
    is_in_parallel = .FALSE.
!$  is_in_parallel = omp_in_parallel()

    !$omp parallel default(none) &
    !$omp          private (row, col, blk, bp, bpe) &
    !$omp          shared (local, nblks) &
    !$omp          shared (rbs, cbs, blki, &
    !$omp                  data, norms, local2global_rows, local2global_cols) &
    !$omp          reduction (max:max_val) if (.not.is_in_parallel)
    !$omp do
    DO blk = 1, nblks
       IF (blki(3,blk) .NE. 0) THEN
          bp = blki(3,blk)
          IF (local) THEN
             row = local2global_rows(blki(1,blk))
             col = local2global_cols(blki(2,blk))
          ELSE
             row = blki(1,blk)
             col = blki(2,blk)
          ENDIF
          bpe = bp + rbs(row) * cbs(col) - 1
          norms(blk) = SQRT (REAL (SUM(ABS(DATA(bp:bpe))**2), KIND=sp))
       ELSE
          norms(blk) = 0.0_sp
       ENDIF
       max_val = MAX(max_val,norms(blk))
    ENDDO
    !$omp end do
    !$omp end parallel
  END SUBROUTINE calc_norms_list_z


