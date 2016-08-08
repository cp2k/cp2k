! **************************************************************************************************
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
! **************************************************************************************************
  SUBROUTINE calc_norms_list_s(norms, nblks,&
       blki, rbs, cbs, DATA, local2global_rows, local2global_cols)
    REAL(kind=sp), DIMENSION(:), INTENT(OUT) :: norms
    INTEGER, INTENT(IN)                      :: nblks
    INTEGER, DIMENSION(3,nblks), INTENT(IN)  :: blki
    INTEGER, DIMENSION(:), INTENT(IN)        :: rbs, cbs
    REAL(kind=real_4), DIMENSION(:), &
      INTENT(IN)                             :: DATA
    INTEGER, DIMENSION(:), INTENT(IN)        :: local2global_rows
    INTEGER, DIMENSION(:), INTENT(IN)        :: local2global_cols

    INTEGER                                  :: blk, bp, bpe, row, col
    REAL(kind=sp)                            :: val

!   ---------------------------------------------------------------------------

    !$omp parallel default(none) &
    !$omp          private (row, col, blk, bp, bpe, val) &
    !$omp          shared (nblks) &
    !$omp          shared (rbs, cbs, blki, &
    !$omp                  data, norms, local2global_rows, local2global_cols)
    !$omp do
    DO blk = 1, nblks
       IF (blki(3,blk) .NE. 0) THEN
          bp = blki(3,blk)
          row = local2global_rows(blki(1,blk))
          col = local2global_cols(blki(2,blk))
          bpe = bp + rbs(row) * cbs(col) - 1
          val = SQRT (REAL (SUM(ABS(DATA(bp:bpe))**2), KIND=sp))
       ELSE
          val = 0.0_sp
       ENDIF
       norms(blk) = val
    ENDDO
    !$omp end do
    !$omp end parallel
  END SUBROUTINE calc_norms_list_s

