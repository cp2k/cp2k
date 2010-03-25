! *****************************************************************************
!> \brief Converts mutable data to linear (array) type.
!>
!> \param[in,out] wm      work matrix to convert
! *****************************************************************************
  SUBROUTINE tree_to_linear_d(wm, error)
    USE btree_I8_k_dp2d_v,&
        ONLY: btree_2d_data_d => dp2d,&
              btree_destroy_d => btree_delete,&
              btree_size_d => btree_get_entries
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: wm
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tree_to_linear_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p, treesize, &
                                                error_handler
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: keys
    REAL(kind=real_8), DIMENSION(:, :), POINTER        :: block_2d
    TYPE(btree_2d_data_d), ALLOCATABLE, &
      DIMENSION(:)                           :: values

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    ! srt = .TRUE. ! Not needed because of the copy
    treesize = btree_size_d(wm%mutable%m%btree_d)
    IF (dbg) WRITE(*,*)routineN//" Declared size",wm%lastblk,&
         "present blocks",treesize
    CALL dbcsr_assert(wm%lastblk .EQ. treesize,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Mismatch in number of blocks",__LINE__,error)
    ALLOCATE (keys(treesize), values(treesize))
    CALL btree_destroy_d (wm%mutable%m%btree_d, keys, values)
    CALL ensure_array_size (wm%row_i, ub=treesize, error=error)
    CALL ensure_array_size (wm%col_i, ub=treesize, error=error)
    CALL dbcsr_unpack_i8_2i4 (keys, wm%row_i,&
         wm%col_i)
    ! For now we also fill the data, sloooowly, but this should
    ! be avoided and the data should be copied directly from the
    ! source in the subroutine's main loop.
    CALL ensure_array_size (wm%blk_p, ub=treesize, error=error)
    CALL dbcsr_data_ensure_size (wm%data_area,&
         wm%datasize, error=error)
    blk_p = 1
    DO blk = 1, treesize
       block_2d => values(blk)%p
       !write(*,*)"block", blk
       !write(*,'(10(1X,F5.1))')block_2d

       IF (.NOT. values(blk)%tr) THEN
          wm%blk_p(blk) = blk_p
       ELSE
          wm%blk_p(blk) = -blk_p
       ENDIF

       CALL dbcsr_data_set (wm%data_area, blk_p,&
            SIZE(block_2d),&
            RESHAPE (block_2d, (/SIZE(block_2d)/)), source_lb=1)
       blk_p = blk_p + SIZE(block_2d)
       DEALLOCATE (block_2d)
    ENDDO
    DEALLOCATE (keys, values)
    CALL dbcsr_mutable_release (wm%mutable)
    !write(*,*)"final data",dbcsr_get_data_size (wm%data_area),&
    !     treesize
    !write(*,'(10(1X,F5.1))')&
    !dbcsr_get_data_p (wm%data_area, 0.0_dp)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE tree_to_linear_d

