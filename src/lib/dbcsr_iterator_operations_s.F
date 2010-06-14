! *****************************************************************************
!> \brief Gets the next data block, single/double precision real/complex
!> \param[in,out] iterator   the iterator
!> \param[out] row           row of the data block
!> \param[out] column        column of the data block
!> \param[out] block         pointer to the data block
!> \param[out] transposed    whether the block data is transposed
!> \param[out] block_number  (optional) block number
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE iterator_next_1d_block_s (iterator, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    REAL(kind=real_4), DIMENSION(:), POINTER :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size,&
         row_offset, col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'iterator_next_1d_block_real(kind=real_4)', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_p, bp, csize, ithread, &
                                                nze, rsize, pos
    REAL(kind=real_4), DIMENSION(:), POINTER :: rbp

!   ---------------------------------------------------------------------------
! If we're pointing to a valid block, return that block.


    IF (iterator%pos .LE. iterator%nblks&
         .AND. iterator%pos .NE. 0) THEN
       row = iterator%row
       column = iterator%col_i(iterator%pos)
       IF (iterator%transpose) CALL swap (row, column)
       blk_p = iterator%blk_p(iterator%pos)
       transposed = blk_p .LT. 0
       bp = ABS (blk_p)
       rsize = iterator%row_size
       csize = iterator%cbs(column)
       nze = rsize * csize
       IF (PRESENT (row_size)) row_size = rsize
       IF (PRESENT (col_size)) col_size = csize
       IF (PRESENT (row_offset)) row_offset = iterator%row_offset
       IF (PRESENT (col_offset)) col_offset = iterator%coff(column)
       rbp => dbcsr_get_data_p (iterator%data_area, 0.0_real_4)
       block => rbp(bp:bp+nze-1)
       IF (PRESENT (block_number)) block_number = iterator%pos
       ! Move to the next non-deleted position.
       CALL iterator_advance (iterator)
       CALL update_row_info (iterator)
    ELSE
       row = 0
       column = 0
       NULLIFY (block)
       IF (PRESENT (block_number)) block_number = 0
    ENDIF
  END SUBROUTINE iterator_next_1d_block_s



! *****************************************************************************
!> \brief Gets the next data block, single/double precision real/complex
!> \param[in,out] iterator   the iterator
!> \param[out] row           row of the data block
!> \param[out] column        column of the data block
!> \param[out] block         pointer to the data block
!> \param[out] transposed    whether the block data is transposed
!> \param[out] block_number  (optional) block number
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE iterator_next_2d_block_s (iterator, row, column,&
       block, transposed,&
       block_number, row_size, col_size, row_offset, col_offset)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    REAL(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size, row_offset, col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'iterator_next_block_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_p, bp, csize, ithread, &
                                                nze, rsize, pos
    REAL(kind=real_4), DIMENSION(:), POINTER :: lin_blk_p, rbp

!   ---------------------------------------------------------------------------
! If we're pointing to a valid block, return that block.

    IF (iterator%pos .LE. iterator%nblks&
         .AND. iterator%pos .NE. 0) THEN
       row = iterator%row
       column = iterator%col_i(iterator%pos)
       IF (iterator%transpose) CALL swap (row, column)
       blk_p = iterator%blk_p(iterator%pos)
       transposed = blk_p .LT. 0
       bp = ABS (blk_p)
       rsize = iterator%row_size
       csize = iterator%cbs(column)
       IF (PRESENT (row_size)) row_size = rsize
       IF (PRESENT (col_size)) col_size = csize
       IF (PRESENT (row_offset)) row_offset = iterator%row_offset
       IF (PRESENT (col_offset)) col_offset = iterator%coff(column)
       nze = rsize * csize
       rbp => dbcsr_get_data_p (iterator%data_area, 0.0_real_4)
       IF (transposed) CALL swap (rsize, csize)
#if defined(__PTR_RANK_REMAP)
       block(1:rsize, 1:csize) => rbp(bp:bp+nze-1)
#elif !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
       lin_blk_p => rbp(bp:bp+nze-1)
       CALL pointer_s_rank_remap2 (block, rsize, csize, lin_blk_p)
#else
       IF (iterator%contiguous_pointers) THEN
          IF (SIZE (iterator%buffer_2d%d%r2_sp, 1) .NE. rsize) THEN
             DEALLOCATE (iterator%buffer_2d%d%r2_sp)
             ALLOCATE (iterator%buffer_2d%d%r2_sp(&
                  rsize, MAX(iterator%max_cbs, iterator%max_rbs)))
          ENDIF
       ENDIF
       iterator%buffer_2d%d%r2_sp =&
            RESHAPE (rbp(bp:bp+nze-1), (/rsize, csize/))
       block => iterator%buffer_2d%d%r2_sp(1:rsize, 1:csize)
#endif /* __PTR_RANK_REMAP || !__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION */
       IF (PRESENT (block_number)) block_number = iterator%pos
       ! Move to the next non-deleted position.
       CALL iterator_advance (iterator)
       CALL update_row_info (iterator)
    ELSE
       row = 0
       column = 0
       NULLIFY (block)
       IF (PRESENT (block_number)) block_number = 0
    ENDIF
  END SUBROUTINE iterator_next_2d_block_s



#if !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
#if defined(__PTR_RANK_REMAP)
! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data using Fortran 2003 pointer
!>        rank remapping.
!>
!> Merge this with dbcsr_util* upon refactoring.
! *****************************************************************************
  SUBROUTINE pointer_s_rank_remap2 (r2p, d1, d2, r1p)
    INTEGER, INTENT(IN)                      :: d1, d2
    REAL(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: r2p
    REAL(kind=real_4), DIMENSION(d1*d2), &
      TARGET                                 :: r1p

    r2p(1:d1,1:d2) => r1p(1:d1*d2)
  END SUBROUTINE pointer_s_rank_remap2
#else
! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data using ugly hacks.
!>
!> Merge this with the same routine in dbcsr_util* upon refactoring.
! *****************************************************************************
  SUBROUTINE pointer_s_rank_remap2 (r2p, d1, d2, r1p)
    REAL(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: r2p
    INTEGER, INTENT(IN)                      :: d1, d2
    REAL(kind=real_4), DIMENSION(d1, *), &
      TARGET                                 :: r1p

    r2p => r1p(1:d1, 1:d2)
  END SUBROUTINE pointer_s_rank_remap2
#endif
#endif

