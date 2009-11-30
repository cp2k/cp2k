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
  SUBROUTINE iterator_next_1d_block_z (iterator, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size,&
         row_offset, col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'iterator_next_1d_block_complex(kind=real_8)', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_p, bp, csize, ithread, &
                                                nze, rsize
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: rbp

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
       IF (PRESENT (row_size)) row_size = rsize
       IF (PRESENT (col_size)) col_size = csize
       IF (PRESENT (row_offset)) row_offset = iterator%row_offset
       IF (PRESENT (col_offset)) col_offset = iterator%coff(column)
       rbp => dbcsr_get_data_p (iterator%data_area, CMPLX(0.0, 0.0, real_8))
       block => rbp(bp:bp+nze-1)
       IF (PRESENT (block_number)) block_number = iterator%pos
       ! Move to the next non-deleted position.
       iterator%pos = iterator%pos+1
       IF (.NOT. iterator%shared) THEN
          CALL find_proper_position(iterator%pos,&
               iterator%row, iterator%nblks, iterator%nblkrows_total,&
               iterator%blk_p, iterator%row_p)
       ELSE
          ithread = 0
!$        ithread = OMP_GET_THREAD_NUM()
          CALL find_proper_position(iterator%pos,&
               iterator%row, iterator%nblks, iterator%nblkrows_total,&
               iterator%blk_p, iterator%row_p,&
               tdist=iterator%tdist, tid=ithread)
       ENDIF
       IF (iterator%row .GT. 0) iterator%row_size = iterator%rbs(iterator%row)
       IF (iterator%row .GT. 0) iterator%row_offset = iterator%roff(iterator%row)
    ELSE
       row = 0
       column = 0
       NULLIFY (block)
       IF (PRESENT (block_number)) block_number = 0
    ENDIF
  END SUBROUTINE iterator_next_1d_block_z

! *****************************************************************************
!> \brief Sets a data pointer.
!> \par Assumptions
!>      Assumes that no memory will be lost when repointing the
!>      pointer in the data area and that the area is initialized.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
! *****************************************************************************
  SUBROUTINE set_data_p_z (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_complex(kind=real_8)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    area%d%c_dp => p
  END SUBROUTINE set_data_p_z

! *****************************************************************************
!> \brief Sets a data pointer.
!> \par Assumptions
!>      Assumes that no memory will be lost when repointing the
!>      pointer in the data area and that the area is initialized.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
! *****************************************************************************
  SUBROUTINE set_data_p_2d_z (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER         :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_2d_complex(kind=real_8)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    area%d%c2_dp => p
  END SUBROUTINE set_data_p_2d_z


! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \par Calling
!>      This routine is hidden behind the dbcsr_get_data interface, hence the
!>      need for the coersion argument.
!> \sa dbcsr_get_data_p_z
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \retval data          pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_c_z (area, coersion, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(kind=real_8), INTENT(IN)            :: coersion
    INTEGER, INTENT(IN), OPTIONAL  :: lb, ub
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: DATA

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_c_z', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
    LOGICAL, PARAMETER             :: dbg = .TRUE.

!   ---------------------------------------------------------------------------

    ! The coersion argument is needed to make this function unique
    ! enough to use in the interface. The following line is added to
    ! remove compiler warnings about unused arguments.
    l = coersion
    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c_dp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c_dp, 1)
          IF (PRESENT (ub)) u = ub
          IF (dbg) THEN
             CALL cp_assert (l .GE. LBOUND (area%d%c_dp, 1),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
             CALL cp_assert (u .LE. UBOUND (area%d%c_dp, 1),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
          ENDIF
          DATA => area%d%c_dp(l:u)
       ELSE
          DATA => area%d%c_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_c_z

! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \par Calling
!>      This routine can be called explicitly.
!> \brief dbcsr_get_data_c_z
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_z (area, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: DATA
    INTEGER, INTENT(IN), OPTIONAL  :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_p_z', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
    LOGICAL, PARAMETER             :: dbg = .TRUE.

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c_dp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c_dp, 1)
          IF (PRESENT (ub)) u = ub
          IF (dbg) THEN
             CALL cp_assert (l .GE. LBOUND (area%d%c_dp, 1),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
             CALL cp_assert (u .LE. UBOUND (area%d%c_dp, 1),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
          ENDIF
          DATA => area%d%c_dp(l:u)
       ELSE
          DATA => area%d%c_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_z

! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \par Calling
!>      This routine can be called explicitly.
!> \brief dbcsr_get_data_c_z
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_2d_z (area, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)            :: area
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER            :: DATA
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_p_2d_z', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(2)          :: l, u
    LOGICAL, PARAMETER             :: dbg = .TRUE.

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c2_dp)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c2_dp)
          IF (PRESENT (ub)) u = ub
          IF (dbg) THEN
             CALL cp_assert (l(1) .GE. LBOUND (area%d%c2_dp, 1),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
             CALL cp_assert (l(2) .GE. LBOUND (area%d%c2_dp, 2),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
             CALL cp_assert (u(1) .LE. UBOUND (area%d%c2_dp, 1),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
             CALL cp_assert (u(2) .LE. UBOUND (area%d%c2_dp, 2),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
          ENDIF
          DATA => area%d%c2_dp(l(1):u(1), l(2):u(2))
       ELSE
          DATA => area%d%c2_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_2d_z



! *****************************************************************************
!> \brief Returns the entire data for a matrix.
!> \par Warning
!>      This routine should only be used within DBCSR code.
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_m_z (matrix, DATA)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: DATA

    CALL get_data_z (matrix%m%data_area, DATA)
  END SUBROUTINE get_data_m_z


! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[out] data      pointer to data
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
! *****************************************************************************
  SUBROUTINE get_data_z (area, DATA, lb, ub)
    TYPE(dbcsr_data_obj), INTENT(IN)  :: area
    COMPLEX(kind=real_8), DIMENSION(:), POINTER    :: DATA
    INTEGER, INTENT(IN), OPTIONAL     :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'get_data_z', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
    LOGICAL, PARAMETER             :: dbg = .TRUE.

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c_dp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c_dp, 1)
          IF (PRESENT (ub)) u = ub
          IF (dbg) THEN
             CALL cp_assert (l .GE. LBOUND (area%d%c_dp, 1),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
             CALL cp_assert (u .LE. UBOUND (area%d%c_dp, 1),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
          ENDIF
          DATA => area%d%c_dp(l:u)
       ELSE
          DATA => area%d%c_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_z


! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[out] data      pointer to data
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
! *****************************************************************************
  SUBROUTINE get_data_2d_z (area, DATA, lb, ub)
    TYPE(dbcsr_data_obj), INTENT(IN)            :: area
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER            :: DATA
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'get_data_2d_z', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(2)          :: l, u
    LOGICAL, PARAMETER             :: dbg = .TRUE.

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c2_dp)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c2_dp)
          IF (PRESENT (ub)) u = ub
          IF (dbg) THEN
             CALL cp_assert (l(1) .GE. LBOUND (area%d%c2_dp, 1),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
             CALL cp_assert (l(2) .GE. LBOUND (area%d%c2_dp, 2),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
             CALL cp_assert (u(1) .LE. UBOUND (area%d%c2_dp, 1),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
             CALL cp_assert (u(2) .LE. UBOUND (area%d%c2_dp, 2),&
                  cp_fatal_level, cp_wrong_args_error, routineN,&
                  "Out of bounds")
          ENDIF
          DATA => area%d%c2_dp(l(1):u(1), l(2):u(2))
       ELSE
          DATA => area%d%c2_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_2d_z

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
  SUBROUTINE iterator_next_2d_block_z (iterator, row, column,&
       block, transposed,&
       block_number, row_size, col_size, row_offset, col_offset)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    COMPLEX(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size, row_offset, col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'iterator_next_block_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_p, bp, csize, ithread, &
                                                nze, rsize
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: lin_blk_p, rbp

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
       IF (PRESENT (row_size)) row_size = rsize
       IF (PRESENT (col_size)) col_size = csize
       IF (PRESENT (row_offset)) row_offset = iterator%row_offset
       IF (PRESENT (col_offset)) col_offset = iterator%coff(column)
       nze = rsize * csize
       rbp => dbcsr_get_data_p (iterator%data_area, CMPLX(0.0, 0.0, real_8))
#if defined(__PTR_RANK_REMAP)
       block(1:rsize, 1:csize) => rbp(bp:bp+nze-1)
#elif !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
       lin_blk_p => rbp(bp:bp+nze-1)
       CALL pointer_z_rank_remap2 (block, rsize, csize, lin_blk_p)
#else
       IF (iterator%contiguous_pointers) THEN
          IF (SIZE (iterator%buffer_2d%d%c2_dp, 1) .NE. rsize) THEN
             DEALLOCATE (iterator%buffer_2d%d%c2_dp)
             ALLOCATE (iterator%buffer_2d%d%c2_dp(&
                  rsize, MAX(iterator%max_cbs, iterator%max_rbs)))
          ENDIF
       ENDIF
       iterator%buffer_2d%d%c2_dp =&
            RESHAPE (rbp(bp:bp+nze-1), (/rsize, csize/))
       block => iterator%buffer_2d%d%c2_dp(1:rsize, 1:csize)
#endif /* __PTR_RANK_REMAP || !__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION */
       IF (PRESENT (block_number)) block_number = iterator%pos
       ! Move to the next non-deleted position.
       iterator%pos = iterator%pos+1
       IF (.NOT.iterator%shared) THEN
          CALL find_proper_position(iterator%pos,&
               iterator%row, iterator%nblks, iterator%nblkrows_total,&
               iterator%blk_p, iterator%row_p)
       ELSE
          ithread = 0
!$        ithread = OMP_GET_THREAD_NUM()
          CALL find_proper_position(iterator%pos,&
               iterator%row, iterator%nblks, iterator%nblkrows_total,&
               iterator%blk_p, iterator%row_p,&
               tdist=iterator%tdist, tid=ithread)
       ENDIF
       IF (iterator%row .GT. 0) iterator%row_size = iterator%rbs(iterator%row)
       IF (iterator%row .GT. 0) iterator%row_offset = iterator%roff(iterator%row)
    ELSE
       row = 0
       column = 0
       NULLIFY (block)
       IF (PRESENT (block_number)) block_number = 0
    ENDIF
  END SUBROUTINE iterator_next_2d_block_z

#if !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
#if defined(__PTR_RANK_REMAP)
! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data using Fortran 2003 pointer
!>        rank remapping.
!>
!> Merge this with dbcsr_util* upon refactoring.
! *****************************************************************************
  SUBROUTINE pointer_z_rank_remap2 (r2p, d1, d2, r1p)
    INTEGER, INTENT(IN)                      :: d1, d2
    COMPLEX(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: r2p
    COMPLEX(kind=real_8), DIMENSION(d1*d2), &
      TARGET                                 :: r1p

    r2p(1:d1,1:d2) => r1p(1:d1*d2)
  END SUBROUTINE pointer_z_rank_remap2
#else
! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data using ugly hacks.
!>
!> Merge this with the same routine in dbcsr_util* upon refactoring.
! *****************************************************************************
  SUBROUTINE pointer_z_rank_remap2 (r2p, d1, d2, r1p)
    COMPLEX(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: r2p
    INTEGER, INTENT(IN)                      :: d1, d2
    COMPLEX(kind=real_8), DIMENSION(d1, *), &
      TARGET                                 :: r1p

    r2p => r1p(1:d1, 1:d2)
  END SUBROUTINE pointer_z_rank_remap2
#endif
#endif

