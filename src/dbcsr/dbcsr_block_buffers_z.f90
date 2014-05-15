!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

  SUBROUTINE dbcsr_buffers_set_p_2d_z (pointer_2d,&
       row, col, rsize, csize, main_tr, base_offset, buffers, buff_tr, error)
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER      :: pointer_2d
    INTEGER, INTENT(IN)                   :: row, col, rsize, csize, base_offset
    LOGICAL, INTENT(IN)                   :: main_tr, buff_tr
    TYPE(dbcsr_block_buffer_obj), INTENT(INOUT) :: buffers
    TYPE(dbcsr_error_type), INTENT(INOUT) :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_set_p_2d_z', &
      routineP = moduleN//':'//routineN
    INTEGER                               :: error_handler, ibuff, rs, cs,&
                                             buff_rs, buff_cs, transposes,&
                                             main_rs, main_cs
    TYPE(dbcsr_data_obj)                  :: buff_area
    LOGICAL :: valid

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    ibuff = dbcsr_buffers_which_is_my (buffers)
    buff_area = dbcsr_buffers_get_area (buffers, ibuff, error=error)
    !
    !CALL dbcsr_assert (main_tr, "EQV", buff_tr, dbcsr_fatal_level, &
    !     dbcsr_wrong_args_error, routineN,&
    !     "Source and buffer should have same transposed status.",&
    !     __LINE__, error=error)
    IF (careful_mod) THEN
       CALL dbcsr_assert (ASSOCIATED (buff_area%d%c2_dp),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Buffer not associated.",__LINE__,error)
    ENDIF
    CALL dbcsr_data_get_sizes (buff_area, rs, cs, valid, error=error)
    CALL dbcsr_assert (valid, dbcsr_fatal_level, dbcsr_internal_error,&
         routineN, "Buffer not setup", __LINE__, error=error)
    IF (main_tr) THEN
       main_rs = csize
       main_cs = rsize
    ELSE
       main_rs = rsize
       main_cs = csize
    ENDIF
    IF (buff_tr) THEN
       buff_rs = csize
       buff_cs = rsize
    ELSE
       buff_rs = rsize
       buff_cs = csize
    ENDIF
    CALL dbcsr_block_partial_copy (dst=buff_area,&
         dst_rs=rsize, dst_cs=csize, dst_tr=buff_tr,&
         src=buffers%b%main, src_offset=ABS(base_offset)-1,&
         src_rs=rsize, src_cs=csize, src_tr=main_tr,&
         dst_r_lb=1, dst_c_lb=1, src_r_lb=1, src_c_lb=1,&
         nrow = rsize, ncol = csize)
    pointer_2d => buff_area%d%c2_dp(1:buff_rs, 1:buff_cs)
    transposes = 0
    IF (buff_tr) transposes = IBSET (transposes, 0)
    IF (main_tr) transposes = IBSET (transposes, 1)
    buffers%b%rcb(1, ibuff) = row
    buffers%b%rcb(2, ibuff) = col
    buffers%b%rcb(3, ibuff) = base_offset
    buffers%b%rcb(4, ibuff) = rsize
    buffers%b%rcb(5, ibuff) = csize
    buffers%b%rcb(6, ibuff) = transposes
    buffers%b%dirty(ibuff) = .TRUE.
    CALL dbcsr_error_stop (error_handler, error)
  END SUBROUTINE dbcsr_buffers_set_p_2d_z
