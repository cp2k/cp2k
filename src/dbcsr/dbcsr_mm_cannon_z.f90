!> \brief Calculates norms of the entire matrix with minimal overhead.
  SUBROUTINE calc_norms_z(norms, nrows,&
       row_p, col_i, blk_p, rbs, cbs, DATA, local, local2global)
    REAL(kind=sp), DIMENSION(:), INTENT(OUT) :: norms
    INTEGER, INTENT(IN)                      :: nrows
    INTEGER, DIMENSION(1:nrows+1), &
      INTENT(IN)                             :: row_p
    INTEGER, DIMENSION(*), INTENT(IN)        :: col_i, blk_p, rbs, cbs
    COMPLEX(kind=real_8), DIMENSION(*), &
      INTENT(IN)                             :: DATA
    LOGICAL, INTENT(IN)                      :: local
    INTEGER, DIMENSION(*), INTENT(IN)        :: local2global

    INTEGER                                  :: blk, bp, bpe, row, row_i, &
                                                row_size

!   ---------------------------------------------------------------------------

    !$omp parallel default(none) &
    !$omp          private (row_i, row, row_size, blk, bp, bpe) &
    !$omp          shared (nrows, local) &
    !$omp          shared (local2global, rbs, cbs, row_p, col_i, blk_p, &
    !$omp                  data, norms)
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
       ENDDO
    ENDDO
    !$omp end do
    !$omp end parallel
  END SUBROUTINE calc_norms_z

!> \brief Calculates norms of the entire matrix with minimal overhead.
  SUBROUTINE calc_norms_list_z(norms, nblks,&
       blki, rbs, cbs, DATA, local, local2global_rows, local2global_cols)
    REAL(kind=sp), DIMENSION(:), INTENT(OUT) :: norms
    INTEGER, INTENT(IN)                      :: nblks
    INTEGER, DIMENSION(3,nblks), INTENT(IN)  :: blki
    INTEGER, DIMENSION(:), INTENT(IN)        :: rbs, cbs
    COMPLEX(kind=real_8), DIMENSION(:), &
      INTENT(IN)                             :: DATA
    LOGICAL, INTENT(IN)                      :: local
    INTEGER, DIMENSION(:), INTENT(IN)        :: local2global_rows
    INTEGER, DIMENSION(:), INTENT(IN)        :: local2global_cols

    INTEGER                                  :: blk, bp, bpe, row, col

!   ---------------------------------------------------------------------------

    !$omp parallel default(none) &
    !$omp          private (row, col, blk, bp, bpe) &
    !$omp          shared (local, nblks) &
    !$omp          shared (rbs, cbs, blki, &
    !$omp                  data, norms, local2global_rows, local2global_cols)
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
    ENDDO
    !$omp end do
    !$omp end parallel
  END SUBROUTINE calc_norms_list_z


