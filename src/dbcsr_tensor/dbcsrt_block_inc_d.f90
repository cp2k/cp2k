! **************************************************************************************************
!> \brief Create block from array, array can be n-dimensional.
!> \param block ...
!> \param sizes ...
!> \param array ...
! **************************************************************************************************
SUBROUTINE create_block_data_d(block, sizes, array)
   TYPE(block_nd), INTENT(OUT)                        :: block
   INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes
   REAL(kind=real_8), DIMENSION(PRODUCT(sizes)),  INTENT(IN) :: array

   ASSOCIATE(blk => block%r_dp)
   block%data_type = dbcsr_type_real_8
   ALLOCATE(blk%sizes(SIZE(sizes)))
   blk%sizes(:) = sizes(:)
   ALLOCATE(blk%blk(PRODUCT(sizes)))
   blk%blk(:) = array(:)
   END ASSOCIATE
END SUBROUTINE

! **************************************************************************************************
!> \brief Create and allocate block, but no data.
!> \param block ...
!> \param sizes ...
! **************************************************************************************************
SUBROUTINE create_block_nodata_d(block, sizes)
   INTEGER, INTENT(IN), DIMENSION(:)                  :: sizes
   TYPE(block_nd_d), INTENT(OUT)            :: block
   ALLOCATE(block%sizes(SIZE(sizes)))
   block%sizes(:) = sizes(:)
   ALLOCATE(block%blk(PRODUCT(sizes)))
END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param block ...
! **************************************************************************************************
SUBROUTINE destroy_block_d (block)
   TYPE(block_nd_d), INTENT(INOUT) :: block
   DEALLOCATE(block%blk)
   DEALLOCATE(block%sizes)
END SUBROUTINE

! **************************************************************************************************
!> \brief add block to buffer.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
! **************************************************************************************************
SUBROUTINE block_buffer_add_block_d(buffer, ndata, index, block)
   TYPE(block_buffer_type), INTENT(INOUT)             :: buffer
   INTEGER, INTENT(IN)                                :: ndata
   REAL(kind=real_8), DIMENSION(ndata), INTENT(IN)              :: block
   INTEGER, DIMENSION(ndims_buffer(buffer)), INTENT(IN) :: index
   INTEGER                                            :: p, ndims, p_data
   CPASSERT(buffer%data_type .EQ. dbcsr_type_real_8)
   ndims = ndims_buffer(buffer)
   p = buffer%endpos
   IF (p .EQ. 0) THEN
      p_data = 0
   ELSE
      p_data = buffer%indx(p, ndims+1)
   ENDIF

   buffer%msg_r_dp(p_data+1:p_data+ndata) = block(:)
   buffer%indx(p+1, 1:ndims) = index(:)
   IF (p > 0) THEN
      buffer%indx(p+1,ndims+1) = buffer%indx(p,ndims+1)+ndata
   ELSE
      buffer%indx(p+1, ndims+1) = ndata
   ENDIF
   buffer%endpos = buffer%endpos+1
END SUBROUTINE

! **************************************************************************************************
!> \brief get next block from buffer. Iterator is advanced only if block is retrieved or advance_iter.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
! **************************************************************************************************
SUBROUTINE block_buffer_get_next_block_d(buffer, ndata, index, block, advance_iter)
   TYPE(block_buffer_type), INTENT(INOUT)             :: buffer
   INTEGER, INTENT(OUT)                                :: ndata
   REAL(kind=real_8), DIMENSION(:), ALLOCATABLE, OPTIONAL, INTENT(OUT) :: block
   INTEGER, DIMENSION(ndims_buffer(buffer)), INTENT(OUT)     :: index
   INTEGER :: p, ndims, p_data
   LOGICAL, INTENT(IN), OPTIONAL                             :: advance_iter
   LOGICAL                                                   :: do_advance

   do_advance = .FALSE.
   IF (PRESENT(advance_iter)) THEN
      do_advance = advance_iter
   ELSE IF (PRESENT(block)) THEN
      do_advance = .TRUE.
   ENDIF

   CPASSERT(buffer%data_type .EQ. dbcsr_type_real_8)
   ndims = ndims_buffer(buffer)
   p = buffer%endpos
   IF (p .EQ. 0) THEN
      p_data = 0
   ELSE
      p_data = buffer%indx(p, ndims+1)
   ENDIF
   IF (p > 0) THEN
      ndata = buffer%indx(p+1, ndims+1)-buffer%indx(p, ndims+1)
   ELSE
      ndata = buffer%indx(p+1, ndims+1)
   ENDIF
   index(:) = buffer%indx(p+1,1:ndims)
   IF (PRESENT(block)) THEN
      ALLOCATE (block(ndata))
      block(:) = buffer%msg_r_dp (p_data+1:p_data+ndata)
   ENDIF

   IF(do_advance) buffer%endpos = buffer%endpos+1
END SUBROUTINE
