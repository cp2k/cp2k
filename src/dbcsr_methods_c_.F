! *****************************************************************************
!> \brief Gets the next data block, single/double precision real/complex
!> \param[in,out] iterator   the iterator
!> \param[out] row           row of the data block
!> \param[out] column        column of the data block
!> \param[out] block         pointer to the data block
!> \param[out] transposed    whether the block data is transposed
!> \param[out] block_number  (optional) block number
! *****************************************************************************
  SUBROUTINE iterator_next_1d_block_c (iterator, row, column, block,&
       transposed, block_number)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    COMPLEX(kind=real_4), DIMENSION(:), POINTER :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number

    CHARACTER(len=*), PARAMETER :: routineN = 'iterator_next_1d_block_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_p, bp, csize, nze, rsize
    COMPLEX(kind=real_4), DIMENSION(:), POINTER :: rbp

!   ---------------------------------------------------------------------------
! If we're pointing to a valid block, return that block.


    IF (iterator%pos .LE. iterator%nblks&
         .AND. iterator%pos .NE. 0) THEN
       row = iterator%row
       column = iterator%col_i(iterator%pos)
       blk_p = iterator%blk_p(iterator%pos)
       transposed = blk_p .LT. 0
       bp = ABS (blk_p)
       rsize = iterator%row_size
       csize = iterator%cbs(column)
       nze = rsize * csize
       rbp => dbcsr_get_data_p (iterator%data_area, CMPLX(0.0, 0.0, real_4))
       block => rbp(bp:bp+nze-1)
       IF (PRESENT (block_number)) block_number = iterator%pos
       ! Move to the next non-deleted position.
       iterator%pos = iterator%pos+1
       CALL find_first_valid_block (iterator%pos, iterator%nblks,&
            iterator%blk_p)
       CALL find_proper_row (iterator%pos, iterator%row,&
            iterator%nblkrows_total, iterator%row_p)
       IF (iterator%row .GT. 0) iterator%row_size = iterator%rbs(iterator%row)
    ELSE
       row = 0
       column = 0
       NULLIFY (block)
       IF (PRESENT (block_number)) block_number = 0
    ENDIF
  END SUBROUTINE iterator_next_1d_block_c

! *****************************************************************************
!> \brief Sets a data pointer.
!>
!> Assumes that no memory will be lost when repointing the pointer in the data
!> area and that the area is initialized.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
! *****************************************************************************
  SUBROUTINE set_data_p_c (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    COMPLEX(kind=real_4), DIMENSION(:), POINTER :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_complex(kind=real_4)', &
      routineP = moduleN//':'//routineN

    area%d%c_sp => p
  END SUBROUTINE set_data_p_c

! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_c_c (area, coersion) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(kind=real_4), INTENT(IN)            :: coersion
    COMPLEX(kind=real_4), DIMENSION(:), POINTER :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%c_sp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_c_c

! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_c (area) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(kind=real_4), DIMENSION(:), POINTER :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%c_sp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_c

! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_m_c (matrix, DATA)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    COMPLEX(kind=real_4), DIMENSION(:), POINTER :: DATA

    IF (ASSOCIATED (matrix%m%data_area%d)) THEN
       DATA => matrix%m%data_area%d%c_sp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_m_c

! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_c (area, DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(kind=real_4), DIMENSION(:), POINTER :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%c_sp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_c

! *****************************************************************************
!> \brief Gets the next data block, single/double precision real/complex
!> \param[in,out] iterator   the iterator
!> \param[out] row           row of the data block
!> \param[out] column        column of the data block
!> \param[out] block         pointer to the data block
!> \param[out] transposed    whether the block data is transposed
!> \param[out] block_number  (optional) block number
! *****************************************************************************
  SUBROUTINE iterator_next_block_c (iterator, row, column, block, transposed,&
       block_number)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    COMPLEX(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number

    CHARACTER(len=*), PARAMETER :: routineN = 'iterator_next_block_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_p, bp, csize, nze, rsize
    COMPLEX(kind=real_4), DIMENSION(:), POINTER :: lin_blk_p, rbp

!   ---------------------------------------------------------------------------
! If we're pointing to a valid block, return that block.

    IF (iterator%pos .LE. iterator%nblks&
         .AND. iterator%pos .NE. 0) THEN
       row = iterator%row
       column = iterator%col_i(iterator%pos)
       blk_p = iterator%blk_p(iterator%pos)
       transposed = blk_p .LT. 0
       bp = ABS (blk_p)
       rsize = iterator%row_size
       csize = iterator%cbs(column)
       nze = rsize * csize
       rbp => dbcsr_get_data_p (iterator%data_area, CMPLX(0.0, 0.0, real_4))
#if defined(__PTR_RANK_REMAP)
       block(1:rsize, 1:csize) => rbp(bp:bp+nze-1)
#elif !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
       lin_blk_p => rbp(bp:bp+nze-1)
       CALL pointer_c_rank_remap2 (block, rsize, csize, lin_blk_p)
#else
       IF (iterator%contiguous_pointers) THEN
          IF (SIZE (iterator%buffer_c_sp, 1) .NE. rsize) THEN
             DEALLOCATE (iterator%buffer_c_sp)
             ALLOCATE (iterator%buffer_c_sp(&
                  rsize, MAX(iterator%max_cbs, iterator%max_rbs)))
          ENDIF
       ENDIF
       iterator%buffer_c_sp =&
            RESHAPE (rbp(bp:bp+nze-1), (/rsize, csize/))
       block => iterator%buffer_c_sp(1:rsize, 1:csize)
#endif /* __PTR_RANK_REMAP || !__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION */
       IF (PRESENT (block_number)) block_number = iterator%pos
       ! Move to the next non-deleted position.
       iterator%pos = iterator%pos+1
       CALL find_first_valid_block (iterator%pos, iterator%nblks,&
            iterator%blk_p)
       CALL find_proper_row (iterator%pos, iterator%row,&
            iterator%nblkrows_total, iterator%row_p)
       IF (iterator%row .GT. 0) iterator%row_size = iterator%rbs(iterator%row)
    ELSE
       row = 0
       column = 0
       NULLIFY (block)
       IF (PRESENT (block_number)) block_number = 0
    ENDIF
  END SUBROUTINE iterator_next_block_c

#if !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data using ugly hacks.
!>
!> Merge this with the same routine in dbcsr_util* upon refactoring.
! *****************************************************************************
  SUBROUTINE pointer_c_rank_remap2 (r2p, d1, d2, r1p)
    COMPLEX(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: r2p
    INTEGER, INTENT(IN)                      :: d1, d2
    COMPLEX(kind=real_4), DIMENSION(d1, *), &
      TARGET                                 :: r1p

    r2p => r1p(1:d1, 1:d2)
  END SUBROUTINE pointer_c_rank_remap2
#endif

