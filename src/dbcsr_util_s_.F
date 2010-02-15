! *****************************************************************************
!> \brief Sets the diagonal of a DBCSR matrix.
!> \param[in] matrix          matrix in which to set the diagonal.
!> \param[in] value           value of the diagonal elements
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_set_diagonal_s (matrix, value, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(kind=real_4), INTENT(IN)                      :: value
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_diagonal_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p, col, local_row, &
                                                my_pcol, nze, row, row_size,&
                                                error_handler
    INTEGER, DIMENSION(:), POINTER           :: local_rows, pcol_dist, &
                                                row_blk_size
    LOGICAL                                  :: in_this_row, valid_index
    REAL(kind=real_4), DIMENSION(:), POINTER           :: data_p

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index (matrix),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Matrix not initialized.",__LINE__,error)
    CALL dbcsr_assert(matrix%m%nblkrows_total, 'EQ', matrix%m%nblkcols_total,&
         dbcsr_failure_level, dbcsr_caller_error, routineP,&
         "Can not set diagonal of non-square matrix.",__LINE__,error)
    CALL dbcsr_assert (dbcsr_type_real_4, 'EQ', matrix%m%data_type, dbcsr_failure_level,&
         dbcsr_caller_error, routineN, "Datatype mismatch.",__LINE__,error)
    IF (.NOT. ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_create (matrix, n=1, error=error)
    ENDIF
    row_blk_size => array_data (matrix%m%row_blk_size)
    CALL dbcsr_get_info (matrix, local_rows = local_rows,&
         proc_col_dist = pcol_dist,&
         my_pcol = my_pcol)
    valid_index = dbcsr_valid_index (matrix)
    ! Try to find the column element in each row. If it's not found,
    ! then add it.
    DO local_row = 1, matrix%m%nblkrows_local
       row = local_rows(local_row)
       ! Only diagonal elements come into play.
       IF (pcol_dist(row) .NE. my_pcol) CYCLE
       in_this_row = .FALSE.
       row_size = row_blk_size(row)
       nze = row_size * row_size
       IF (valid_index) THEN
          blk_loop: DO blk = matrix%m%row_p(row)+1, matrix%m%row_p(row+1)
             col = matrix%m%col_i(blk)
             IF (col .LT. row) CYCLE blk_loop
             IF (col .GT. row) EXIT
             IF (col .EQ. row) THEN
                ! Sets the diagonal in the existing block.
                blk_p = matrix%m%blk_p(blk)
                data_p => pointer_view (&
                        dbcsr_get_data_p (matrix%m%data_area, 0.0_real_4),&
                        blk_p, blk_p+nze-1&
                     )
                CALL set_block_diagonal (data_p, value, row_size)
                !CALL set_block_diagonal (matrix%m%data(blk_p:blk_p+nze-1),&
                !     value, row_size)
                in_this_row = .TRUE.
                EXIT blk_loop ! Don't bother with the remaining columns
             ENDIF
          ENDDO blk_loop
       ENDIF ! valid_index
       ! Add the column if it doesn't already exist.
       IF (.NOT. in_this_row) THEN
          DBG 'Adding at row',row,'dimen',row_size,'size',nze,'pos',&
               matrix%m%wms(1)%datasize+1
          CALL add_work_coordinate (matrix%m%wms(1), row, row,&
               matrix%m%wms(1)%datasize+1,&
               error=error)
          ! add_work_coordinate increments the block counter
          CALL dbcsr_data_ensure_size (matrix%m%wms(1)%data_area,&
               matrix%m%wms(1)%datasize + nze, error=error)
          !CALL ensure_array_size (matrix%m%w%data, ub=matrix%m%w%datasize+nze,&
          !     error=error)
          data_p => pointer_view (&
               dbcsr_get_data_p (matrix%m%wms(1)%data_area, 0.0_real_4),&
               matrix%m%wms(1)%datasize+1, matrix%m%wms(1)%datasize+nze)
          CALL set_block_diagonal(data_p, value, row_size)
          !CALL set_block_diagonal (&
          !     matrix%m%w%data(matrix%m%w%datasize+1 : matrix%m%w%datasize+nze),&
          !     value, row_size)
          matrix%m%wms(1)%datasize = matrix%m%wms(1)%datasize + nze
       ENDIF
    ENDDO ! row
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_set_diagonal_s

! *****************************************************************************
!> \brief Sets the diagonal of a square data block represented as a 1-D array.
!>        Non-diagonal elements are set to 0.
!> \param[out] block_data     sets diagonal in this data block
!> \param[in] value           value of the diagonal elements
!> \param[in] d               dimension of the square data block
! *****************************************************************************
  PURE SUBROUTINE set_block_diagonal_s (block_data, value, d)
    REAL(kind=real_4), DIMENSION(:), INTENT(OUT) :: block_data
    REAL(kind=real_4), INTENT(IN)                :: value
    INTEGER, INTENT(IN)                :: d

    INTEGER                            :: i

!   ---------------------------------------------------------------------------

    block_data(:) = 0
    FORALL (i = 0 : d-1)
       block_data(1 + i*d + i) = value
    END FORALL
  END SUBROUTINE set_block_diagonal_s


! *****************************************************************************
!> \brief Copies a block subset
!> \see block_partial_copy_a
! *****************************************************************************
  SUBROUTINE block_partial_copy_s(dst, dst_rs, dst_cs, dst_tr,&
       src, src_rs, src_cs, src_tr,&
       dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
       dst_offset, src_offset)
    REAL(kind=real_4), DIMENSION(:), &
      INTENT(INOUT)                          :: dst
    INTEGER, INTENT(IN)                      :: dst_rs, dst_cs
    INTEGER, INTENT(IN)                      :: src_offset, dst_offset
    LOGICAL                                  :: dst_tr
    REAL(kind=real_4), DIMENSION(:), &
      INTENT(IN)                             :: src
    INTEGER, INTENT(IN)                      :: src_rs, src_cs
    LOGICAL                                  :: src_tr
    INTEGER, INTENT(IN)                      :: dst_r_lb, dst_c_lb, src_r_lb, &
                                                src_c_lb, nrow, ncol

    CHARACTER(len=*), PARAMETER :: routineN = 'block_partial_copy_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, row

!   ---------------------------------------------------------------------------
! Factors out the 4 combinations to remove branches from the inner loop. rs is the logical row size so it always remains the leading dimension.
    IF (.NOT. dst_tr .AND. .NOT. src_tr) THEN
       FORALL (col = 0:ncol-1, row=0:nrow-1)
          dst(dst_offset+dst_r_lb+row + (dst_c_lb+col-1)*dst_rs) &
                = src(src_offset+src_r_lb+row+(src_c_lb+col-1)*src_rs)
       END FORALL
       !DO col = 1, ncol
       !   DO row = 1, nrow
       !      dst(dst_offset+dst_r_lb+row + (dst_c_lb+col-1)*dst_rs) &
       !    = src(src_offset+src_r_lb+row + (src_c_lb+col-1)*src_rs)
       !   ENDDO
       !ENDDO
    ELSEIF (dst_tr .AND. .NOT. src_tr) THEN
       FORALL (col = 0:ncol-1, row=0:nrow-1)
          dst(dst_offset+dst_c_lb+col + (dst_r_lb+row-1)*dst_rs) &
              = src(src_offset+src_r_lb+row+(src_c_lb+col-1)*src_rs)
       END FORALL
    ELSEIF (.NOT. dst_tr .AND. src_tr) THEN
       FORALL (col = 0:ncol-1, row=0:nrow-1)
          dst(dst_offset+dst_r_lb+row + (dst_c_lb+col-1)*dst_rs) &
             = src(src_offset+src_c_lb+col+(src_r_lb+row-1)*src_rs)
       END FORALL
    ELSEIF (dst_tr .AND. src_tr) THEN
       FORALL (col = 0:ncol-1, row=0:nrow-1)
          dst(dst_offset+dst_c_lb+col + (dst_r_lb+row-1)*dst_rs)&
             = src(src_offset + src_c_lb+col+(src_r_lb+row-1)*src_rs)
       END FORALL
    ENDIF
  END SUBROUTINE block_partial_copy_s



! *****************************************************************************
!> \brief Returns a pointer with different bounds.
!> \param[in] original   original data pointer
!> \param[in] lb, ub     lower and upper bound for the new pointer view
!> \param[out] view      new pointer
! *****************************************************************************
  FUNCTION pointer_view_s (original, lb, ub) RESULT (view)
    REAL(kind=real_4), DIMENSION(:), POINTER :: original, view
    INTEGER, INTENT(IN)                  :: lb, ub
    view => original(lb:ub)
  END FUNCTION pointer_view_s


#if !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data using ugly hacks.
! *****************************************************************************
  SUBROUTINE pointer_s_rank_remap2 (r2p, d1, d2, r1p)
    REAL(kind=real_4), DIMENSION(:,:), POINTER :: r2p
    INTEGER, INTENT(IN)              :: d1, d2
    REAL(kind=real_4), DIMENSION(d1,*), TARGET :: r1p
    r2p => r1p(1:d1, 1:d2)
  END SUBROUTINE pointer_s_rank_remap2
#endif

! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data.
!> \param[out] r2p       pointer with rank 2
!> \param[in] d1, d2     dimensions for rank 2 pointer
!> \param[out] r1p       existing pointer with rank1
! *****************************************************************************
  SUBROUTINE pointer_set_r2_from_r1_s (r2p, shp, r1p, r2buff)
    REAL(kind=real_4), DIMENSION(:,:), POINTER    :: r2p
    INTEGER, DIMENSION(2), INTENT(IN)   :: shp
    REAL(kind=real_4), DIMENSION(:), POINTER      :: r1p
    REAL(kind=real_4), DIMENSION(:,:), POINTER    :: r2buff
    INTEGER                             :: m, n

    m = shp(1)
    n = shp(2)
#if defined(__PTR_RANK_REMAP)
    r2p(1:m, 1:n) => r1p(1:m*n)
#elif !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
    CALL pointer_s_rank_remap2 (r2p, m, n, r1p)
#else
    r2buff(1:m, 1:n) = RESHAPE (r1p, shp)
    r2p => r2buff(1:m, 1:n)
#endif
  END SUBROUTINE pointer_set_r2_from_r1_s

! *****************************************************************************
!> \brief Copy and transpose block.
!> \param[out] extent_out     output matrix in the form of a 1-d array
!> \param[in] extent_in       input matrix in the form of a 1-d array
!> \param[in] rows, columns   input matrix size
! *****************************************************************************
  PURE SUBROUTINE block_transpose_copy_s(extent_out, extent_in,&
       rows, columns)
    REAL(kind=real_4), DIMENSION(:), INTENT(OUT) :: extent_out
    REAL(kind=real_4), DIMENSION(:), INTENT(IN)  :: extent_in
    INTEGER, INTENT(IN)                :: rows, columns

    CHARACTER(len=*), PARAMETER :: routineN = 'block_transpose_copy_s', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    extent_out = RESHAPE(TRANSPOSE(&
         RESHAPE(extent_in, (/rows, columns/))), (/rows*columns/))
  END SUBROUTINE block_transpose_copy_s

! *****************************************************************************
!> \brief Copy a block
!> \param[out] extent_out     output matrix in the form of a 2-d array
!> \param[in] extent_in       input matrix in the form of a 1-d array
!> \param[in] rows, columns   input matrix size
! *****************************************************************************
  PURE SUBROUTINE block_copy_2d1d_s(extent_out, extent_in,&
       rows, columns)
    INTEGER, INTENT(IN)                           :: rows, columns
    REAL(kind=real_4), DIMENSION(rows,columns), INTENT(OUT) :: extent_out
    REAL(kind=real_4), DIMENSION(:), INTENT(IN)             :: extent_in

    CHARACTER(len=*), PARAMETER :: routineN = 'block_transpose_copy_2d1d_s', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    extent_out = RESHAPE(extent_in, (/rows, columns/))
  END SUBROUTINE block_copy_2d1d_s


! *****************************************************************************
!> \brief Copy and transpose block.
!> \param[out] extent_out     output matrix in the form of a 2-d array
!> \param[in] extent_in       input matrix in the form of a 1-d array
!> \param[in] rows, columns   input matrix size
! *****************************************************************************
  PURE SUBROUTINE block_transpose_copy_2d1d_s(extent_out, extent_in,&
       rows, columns)
    INTEGER, INTENT(IN)                           :: rows, columns
    REAL(kind=real_4), DIMENSION(columns,rows), INTENT(OUT) :: extent_out
    REAL(kind=real_4), DIMENSION(:), INTENT(IN)             :: extent_in

    CHARACTER(len=*), PARAMETER :: routineN = 'block_transpose_copy_2d1d_s', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    extent_out = TRANSPOSE(RESHAPE(extent_in, (/rows, columns/)))
  END SUBROUTINE block_transpose_copy_2d1d_s

! *****************************************************************************
!> \brief Copy and transpose block.
!> \param[out] extent_out     output matrix in the form of a 1-d array
!> \param[in] extent_in       input matrix in the form of a 2-d array
!> \param[in] rows, columns   input matrix size
! *****************************************************************************
  PURE SUBROUTINE block_copy_1d2d_s(extent_out, extent_in,&
       rows, columns)
    REAL(kind=real_4), DIMENSION(:), INTENT(OUT)            :: extent_out
    INTEGER, INTENT(IN)                           :: rows, columns
    REAL(kind=real_4), DIMENSION(rows,columns), INTENT(IN)  :: extent_in

    CHARACTER(len=*), PARAMETER :: routineN = 'block_copy_1d2d_s', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    extent_out = RESHAPE(extent_in, (/rows*columns/))
  END SUBROUTINE block_copy_1d2d_s


! *****************************************************************************
!> \brief Copy and transpose block.
!> \param[out] extent_out     output matrix in the form of a 1-d array
!> \param[in] extent_in       input matrix in the form of a 2-d array
!> \param[in] rows, columns   input matrix size
! *****************************************************************************
  PURE SUBROUTINE block_transpose_copy_1d2d_s(extent_out, extent_in,&
       rows, columns)
    REAL(kind=real_4), DIMENSION(:), INTENT(OUT)            :: extent_out
    INTEGER, INTENT(IN)                           :: rows, columns
    REAL(kind=real_4), DIMENSION(rows,columns), INTENT(IN)  :: extent_in

    CHARACTER(len=*), PARAMETER :: routineN = 'block_transpose_copy_1d2d_s', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    extent_out = RESHAPE(TRANSPOSE(extent_in), (/rows*columns/))
  END SUBROUTINE block_transpose_copy_1d2d_s


! *****************************************************************************
!> \brief In-place block transpose.
!> \param[in,out] extent      Matrix in the form of a 1-d array
!> \param[in] rows, columns   Input matrix size
!> \param error     cp2k error
! *****************************************************************************
  PURE SUBROUTINE block_transpose_inplace_s(extent, rows, columns)
    REAL(kind=real_4), DIMENSION(:), &
      INTENT(INOUT)                          :: extent
    INTEGER, INTENT(IN)                      :: rows, columns

    CHARACTER(len=*), PARAMETER :: routineN = 'block_transpose_inplace_s', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    extent = RESHAPE(TRANSPOSE(&
         RESHAPE(extent, (/rows, columns/))), (/rows*columns/))
  END SUBROUTINE block_transpose_inplace_s


! *****************************************************************************
!> \brief Copy data from a double real array to a data area
!>
!> There are no checks done for correctness!
!> \param[in] dst        destination data area
!> \param[in] lb         lower bound for destination (and source if
!>                       not given explicity)
!> \param[in] data_size  number of elements to copy
!> \param[in] src        source data array
!> \param[in] source_lb  (optional) lower bound of source
! *****************************************************************************
  SUBROUTINE dbcsr_data_set_as (dst, lb, data_size, src, source_lb)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    INTEGER, INTENT(IN)                      :: lb, data_size
    REAL(kind=real_4), DIMENSION(:), INTENT(IN)        :: src
    INTEGER, INTENT(IN), OPTIONAL            :: source_lb
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_set_as', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: lb_s, ub, ub_s
    TYPE(dbcsr_error_type)                   :: error
!   ---------------------------------------------------------------------------
    IF (bcsr_debug) THEN
       CALL dbcsr_assert (ASSOCIATED (dst%d),&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Target data area must be setup.",__LINE__,error)
       CALL dbcsr_assert (SIZE(src) .GE. data_size,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Not enough source data.",__LINE__,error)
       CALL dbcsr_assert (dst%d%data_type .EQ. dbcsr_type_real_4, dbcsr_failure_level,&
            dbcsr_caller_error, routineN, "Data type mismatch.",__LINE__,error)
    ENDIF
    ub = lb + data_size - 1
    IF (PRESENT (source_lb)) THEN
       lb_s = source_lb
       ub_s = source_lb + data_size-1
    ELSE
       lb_s = lb
       ub_s = ub
    ENDIF
    dst%d%r_sp(lb:ub) = src(lb_s:ub_s)
  END SUBROUTINE dbcsr_data_set_as


! *****************************************************************************
!> \brief Ensures that an array is appropriately large.
!> \param[in,out] array       array to verify and possibly resize
!> \param[in] lb    (optional) desired array lower bound
!> \param[in] ub    desired array upper bound
!> \param[in] factor          (optional) factor by which to exagerrate
!>                            enlargements
!> \param[in] nocopy          (optional) copy array on enlargement; default
!>                            is to copy
!> \param[in] special         (optional) use MPI-allocated memory; default
!>                            is no
!> \param[in] zero_pad        (optional) zero new allocations; default is to
!>                            write nothing
! *****************************************************************************
  SUBROUTINE ensure_array_size_s(array, lb, ub, factor,&
       nocopy, special, zero_pad, error)
    REAL(kind=real_4), DIMENSION(:), POINTER           :: array
    INTEGER, INTENT(IN), OPTIONAL            :: lb
    INTEGER, INTENT(IN)                      :: ub
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, special, zero_pad
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ensure_array_size_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lb_new, lb_orig, stat, &
                                                ub_new, ub_orig, old_size,&
                                                size_increase, error_handler
    LOGICAL                                  :: dbg, docopy = .TRUE., &
                                                failure = .TRUE., info, spec,&
                                                pad
    REAL(kind=real_4), DIMENSION(:), POINTER           :: newarray

!   ---------------------------------------------------------------------------
    CALL dbcsr_error_set(routineN, error_handler, error)
    dbg = .FALSE.
    info = .FALSE.
    CALL dbcsr_set_debug(dbg, info=info)

    docopy = .TRUE.
    IF (PRESENT (nocopy)) docopy = .NOT. nocopy
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    lb_new = 1
    IF (PRESENT (lb)) lb_new = lb
    pad = .FALSE.
    IF (PRESENT (zero_pad)) pad = zero_pad
    !> Creates a new array if it doesn't yet exist.
    IF (.NOT.ASSOCIATED(array)) THEN
       failure = .FALSE.
       IF (spec) THEN
          CALL dbcsr_assert (lb.EQ.1, dbcsr_warning_level, dbcsr_internal_error,&
               routineN, "Special memory can only start at 1.",__LINE__,error)
          CALL mp_allocate(array, ub-lb_new+1, stat=stat)
       ELSE
          ALLOCATE(array(lb_new:ub), stat=stat)
       ENDIF
       IF (stat /= 0) CALL stop_memory(routineP, 'array', ub-lb+1)
       IF (pad) array(:) = 0.0_real_4
       CALL dbcsr_error_stop(error_handler, error)
       RETURN
    ENDIF
    lb_orig = LBOUND(array,1)
    ub_orig = UBOUND(array,1)
    old_size = ub_orig - lb_orig + 1
    ! The existing array is big enough.
    IF (lb_orig.LE.lb_new .AND. ub_orig.GE.ub) THEN
       CALL dbcsr_error_stop(error_handler, error)
       RETURN
    ENDIF
    ! A reallocation must be performed
    IF(dbg) WRITE(*,*)routineP//' Current bounds are',lb_orig,':',ub_orig,&
         '; special?',spec
    !CALL timeset(routineN,timing_handle)
    IF (lb_orig.GT.lb_new) THEN
       IF (PRESENT(factor)) THEN
          size_increase = lb_orig - lb_new
          size_increase = MAX (NINT(REAL(size_increase)*factor),&
                               NINT(REAL(old_size)*factor))
          lb_new = MIN (lb_orig, lb_new - size_increase)
       ELSE
          lb_new = lb_orig
       ENDIF
    ENDIF
    IF (ub_orig.LT.ub) THEN
       IF (PRESENT(factor)) THEN
          size_increase = ub - ub_orig
          size_increase = MAX (NINT(REAL(size_increase)*factor),&
                               NINT(REAL(old_size)*factor))
          ub_new = MAX (ub_orig, ub + size_increase)
       ELSE
          ub_new = ub
       ENDIF
    ELSE
       ub_new = ub
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Resizing to bounds',lb_new,':',ub_new
    IF(.NOT.docopy) THEN
       IF (spec) THEN
          CALL mp_deallocate(array)
       ELSE
          DEALLOCATE(array)
       ENDIF
    ENDIF
    IF (spec) THEN
       CALL dbcsr_assert (lb_new.EQ.1, dbcsr_warning_level, dbcsr_internal_error,&
            routineN, "Special memory can only start at 1.",__LINE__,error)
       CALL mp_allocate(newarray, ub_new-lb_new+1, stat=stat)
    ELSE
       ALLOCATE(newarray(lb_new:ub_new), stat=stat)
    ENDIF
    IF(stat /= 0) CALL stop_memory(routineP, 'newarray', ub_new-lb_new+1)
    IF(docopy) THEN
       IF(dbg) CALL dbcsr_assert(lb_new.LE.lb_orig .AND. ub_new.GE.ub_orig,&
            dbcsr_failure_level, dbcsr_internal_error, routineP,&
            "Old extent exceeds the new one.",__LINE__,error)
       newarray(lb_orig:ub_orig) = array(lb_orig:ub_orig)
       IF (spec) THEN
          CALL mp_deallocate(array)
       ELSE
          DEALLOCATE(array)
       ENDIF
       IF (pad) THEN
          array(ub_orig+1:ub_new) = 0.0_real_4
          array(lb_new:lb_orig-1) = 0.0_real_4
       ENDIF
    ELSEIF (pad) THEN
       array(:) = 0.0_real_4
    ENDIF
    array => newarray
    IF(dbg.AND.bcsr_verbose)&
         WRITE(*,*)routineP//' New array',array(lb_orig:ub_orig)
    IF (dbg) WRITE(*,*)routineP//' New array size', SIZE(array)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE ensure_array_size_s


! *****************************************************************************
!> \brief Copies an array into another array, ensuring that the target array
!>        large enough.
!> \par Makes use of FORTRAN's array bounds for array sections.
!> \param to        array to copy to
!> \param[in] from  array to copy from
!> \param[in] offset          starting position in the target array
!> \param[in] factor          factor for size increase
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE addto_array_as(to, from, offset, factor, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: to
    REAL(kind=real_4), DIMENSION(:), INTENT(IN)        :: from
    INTEGER, INTENT(IN), OPTIONAL            :: offset
    REAL, INTENT(IN), OPTIONAL               :: factor
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'addto_array_as', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: last, lb_t, s_f, s_t, start, &
                                                timing_handle, ub_t, error_handler
    LOGICAL                                  :: dbg, spec
    REAL                                     :: f

!   ---------------------------------------------------------------------------
    CALL dbcsr_error_set(routineN, error_handler, error)
    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    spec = to%d%special
    s_f = SIZE(from, 1)
    IF (dbg) THEN
       WRITE(*,*)routineP//' from size',SIZE(from),';', SIZE(from,1)
       WRITE(*,*)routineP//' from lbound',LBOUND(from),';', LBOUND(from,1)
       WRITE(*,*)routineP//' from ubound',UBOUND(from),';', UBOUND(from,1)
    ENDIF
    IF (.NOT. ASSOCIATED (to%d%r_sp)) THEN
       IF (PRESENT (offset)) THEN
          start = offset
          last = offset + SIZE(from)-1
          CALL ensure_array_size(to%d%r_sp, lb=start,&
               ub=last,&
               factor=1.0, nocopy=.TRUE., special=spec, error=error)
       ELSE
          start = 1
          last = SIZE(from)
          CALL ensure_array_size(to%d%r_sp, ub=last,&
               factor=1.0, nocopy=.TRUE., special=spec, error=error)
       ENDIF
    ELSE
       ub_t = UBOUND(to%d%r_sp, 1)
       lb_t = LBOUND(to%d%r_sp, 1)
       s_t = SIZE(to%d%r_sp)
       start = lb_t
       IF (PRESENT (offset)) start = offset
       ! assert offset > lb_t
       last = start+s_f-1
       IF(dbg) WRITE(*,*)routineP//' target array bounds',lb_t,':',ub_t,&
            'starting at',start,'last',last,'source size', s_f
       IF (ub_t .LT. last) THEN
          IF (PRESENT (factor)) THEN
             f = factor
          ELSE
             f = 1.4142136
          ENDIF
          !ub_t = MAX(lb_t + INT(s_t * f), lb_t+s_f-1)
          ub_t = MAX(start + INT(s_t * f), last)
          IF(dbg) WRITE(*,*)routineP//' expanding array from :',&
               UBOUND(to%d%r_sp),&
               'to',ub_t
          !CALL reallocate(to, lb_t, ub_t)
          CALL ensure_array_size(to%d%r_sp, lb_t, ub_t,&
               factor=1.0, nocopy=.FALSE.,&
               special=spec, error=error)
       ENDIF
    ENDIF
    to%d%r_sp(start:last) = from(:)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE addto_array_as

! *****************************************************************************
!> \brief Copies an array into another array, ensuring that the target array
!>        is large enough.
!> \par Makes use of FORTRAN's array bounds for array sections.
!> \param to        array to copy to
!> \param[in] from  array to copy from
!> \param[in] offset          starting position in the target array
!> \param[in] factor          factor for size increase
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE addto_array_s(to, from, offset, factor, special, error)
    REAL(kind=real_4), DIMENSION(:), POINTER           :: to
    REAL(kind=real_4), DIMENSION(:), INTENT(IN)        :: from
    INTEGER, INTENT(IN), OPTIONAL            :: offset
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: special
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'addto_array_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: last, lb_t, s_f, s_t, start, &
                                                ub_t, error_handler
    LOGICAL                                  :: dbg, spec
    REAL                                     :: f

!   ---------------------------------------------------------------------------
    CALL dbcsr_error_set(routineN, error_handler, error)
    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    s_f = SIZE(from, 1)
    IF (dbg) THEN
       WRITE(*,*)routineP//' from size',SIZE(from),';', SIZE(from,1)
       WRITE(*,*)routineP//' from lbound',LBOUND(from),';', LBOUND(from,1)
       WRITE(*,*)routineP//' from ubound',UBOUND(from),';', UBOUND(from,1)
    ENDIF
    ub_t = UBOUND(to, 1)
    lb_t = LBOUND(to, 1)
    s_t = SIZE(to)
    start = lb_t
    IF (PRESENT (offset)) start = offset
    ! assert offset > lb_t
    last = start+s_f-1
    IF(dbg) WRITE(*,*)routineP//' target array bounds',lb_t,':',ub_t,&
         'starting at',start,'last',last,'source size', s_f
    IF (ub_t .LT. last) THEN
       IF (PRESENT (factor)) THEN
          f = factor
       ELSE
          f = 1.4142136
       ENDIF
       !ub_t = MAX(lb_t + INT(s_t * f), lb_t+s_f-1)
       ub_t = MAX(start + INT(s_t * f), last)
       IF(dbg) WRITE(*,*)routineP//' expanding array from :',UBOUND(to),&
            'to',ub_t
       !CALL reallocate(to, lb_t, ub_t)
       CALL ensure_array_size(to, lb_t, ub_t, factor=1.0, nocopy=.FALSE.,&
            special=spec, error=error)
    ENDIF
    to(start:last) = from(:)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE addto_array_s


! *****************************************************************************
!> \brief Converts mutable data to linear (array) type.
!>
!> \param[in,out] wm      work matrix to convert
! *****************************************************************************
  SUBROUTINE tree_to_linear_s(wm, error)
    USE btree_I8_k_sp2d_v,&
        ONLY: btree_2d_data_s => sp2d,&
              btree_destroy_s => btree_delete,&
              btree_size_s => btree_get_entries
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: wm
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tree_to_linear_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p, i, treesize, &
         error_handler
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: keys
    REAL(kind=real_4), DIMENSION(:, :), POINTER        :: block_2d
    TYPE(btree_2d_data_s), ALLOCATABLE, &
      DIMENSION(:)                           :: values

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    IF (dbg) WRITE(*,'(1X,A,I5,A)')routineN//" Work matrix",i,&
         " has mutable data"
    ! srt = .TRUE. ! Not needed because of the copy
    treesize = btree_size_s(wm%mutable%m%btree_s)
    IF (dbg) WRITE(*,*)routineN//" Declared size",wm%lastblk,&
         "present blocks",treesize
    CALL dbcsr_assert(wm%lastblk .EQ. treesize,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Mismatch in number of blocks",__LINE__,error)
    ALLOCATE (keys(treesize), values(treesize))
    CALL btree_destroy_s (wm%mutable%m%btree_s, keys, values)
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
  END SUBROUTINE tree_to_linear_s

! *****************************************************************************
!> \brief Row/column and global all-to-all
!> \par Communicator selection
!>      Uses row and column communicators for row/column
!>      sends. Remaining sends are performed using the global
!>      communicator.  Point-to-point isend/irecv are used if ptp is
!>      set, otherwise a alltoall collective call is issued.
!> \param[in] mp_env         MP Environment
!> \param[in] most_ptp       (optional) Use point-to-point for row/column;
!>                           default is no
!> \param[in] remainder_ptp  (optional) Use point-to-point for remaining;
!>                           default is no
!> \param[in] no_hybrid      (optional) Use regular global collective; default
!>                           is no
!> \see mp_alltoall
! *****************************************************************************
  SUBROUTINE hybrid_alltoall_s1 (sb, nsb, scount, sdispl,&
       rb, nrb, rcount, rdispl, mp_env, most_ptp, remainder_ptp, no_hybrid)
    INTEGER, INTENT(in)                      :: nsb
    REAL(kind=real_4), DIMENSION(1:nsb), &
      INTENT(in), TARGET                     :: sb
    INTEGER, DIMENSION(:), INTENT(IN)        :: scount, sdispl
    INTEGER, INTENT(in)                      :: nrb
    REAL(kind=real_4), DIMENSION(1:nrb), &
      INTENT(INOUT), TARGET                  :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    LOGICAL, INTENT(in), OPTIONAL            :: most_ptp, remainder_ptp,&
                                                no_hybrid

    CHARACTER(len=*), PARAMETER :: routineN = 'hybrid_alltoall_s', &
      routineP = moduleN//':'//routineN

    INTEGER :: all_group, mynode, mypcol, myprow, nall_rr, nall_sr, ncol_rr, &
      ncol_sr, npcols, nprows, nrow_rr, nrow_sr, numnodes, dst, src,&
      prow, pcol, send_cnt, recv_cnt, tag, grp, i
    INTEGER, ALLOCATABLE, DIMENSION(:) :: all_rr, all_sr, col_rr, col_sr, &
      new_rcount, new_rdispl, new_scount, new_sdispl, row_rr, row_sr
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    LOGICAL                                  :: most_collective, &
                                                remainder_collective, no_h
    REAL(kind=real_4), DIMENSION(:), POINTER           :: send_data_p, recv_data_p
    TYPE(dbcsr_mp_obj)                       :: mpe

    !CALL dbcsr_assert (mp_env%mp%subgroups_defined, dbcsr_warning_level,&
    !     dbcsr_caller_error, routineN, "Row/col communicators undefined.")
    IF (.NOT. dbcsr_mp_has_subgroups (mp_env)) THEN
       mpe = mp_env
       CALL dbcsr_mp_grid_setup (mpe)
    ENDIF
    most_collective = .TRUE.
    remainder_collective = .TRUE.
    no_h = .FALSE.
    IF (PRESENT (most_ptp)) most_collective = .NOT. most_ptp
    IF (PRESENT (remainder_ptp)) remainder_collective = .NOT. remainder_ptp
    IF (PRESENT (no_hybrid)) no_h = no_hybrid
    all_group = dbcsr_mp_group (mp_env)
    subgrouped: IF (mp_env%mp%subgroups_defined .AND. .NOT. no_h) THEN
       mynode = dbcsr_mp_mynode (mp_env)
       numnodes = dbcsr_mp_numnodes (mp_env)
       nprows = dbcsr_mp_nprows (mp_env)
       npcols = dbcsr_mp_npcols (mp_env)
       myprow = dbcsr_mp_myprow (mp_env)
       mypcol = dbcsr_mp_mypcol (mp_env)
       pgrid => dbcsr_mp_pgrid (mp_env)
       ALLOCATE (row_sr(0:npcols-1)) ; nrow_sr = 0
       ALLOCATE (row_rr(0:npcols-1)) ; nrow_rr = 0
       ALLOCATE (col_sr(0:nprows-1)) ; ncol_sr = 0
       ALLOCATE (col_rr(0:nprows-1)) ; ncol_rr = 0
       ALLOCATE (all_sr(0:numnodes-1)) ; nall_sr = 0
       ALLOCATE (all_rr(0:numnodes-1)) ; nall_rr = 0
       ALLOCATE (new_scount(numnodes), new_rcount(numnodes))
       ALLOCATE (new_sdispl(numnodes), new_rdispl(numnodes))
       IF (.NOT.remainder_collective) THEN
          CALL remainder_point_to_point ()
       ENDIF
       IF (.NOT.most_collective) THEN
          CALL most_point_to_point ()
       ELSE
          CALL most_alltoall ()
       ENDIF
       IF (remainder_collective) THEN
          CALL remainder_alltoall ()
       ENDIF
       ! Wait for all issued sends and receives.
       IF (.NOT.most_collective) THEN
          CALL mp_waitall (row_sr(0:nrow_sr-1))
          CALL mp_waitall (col_sr(0:ncol_sr-1))
          CALL mp_waitall (row_rr(0:nrow_rr-1))
          CALL mp_waitall (col_rr(0:ncol_rr-1))
       END IF
       IF (.NOT.remainder_collective) THEN
          CALL mp_waitall (all_sr(1:nall_sr))
          CALL mp_waitall (all_rr(1:nall_rr))
       ENDIF
    ELSE
       CALL mp_alltoall (sb, scount, sdispl,&
            rb, rcount, rdispl,&
            all_group)
    ENDIF subgrouped
  CONTAINS
    SUBROUTINE most_alltoall()
      FORALL (pcol = 0 : npcols-1)
         new_scount(1+pcol) = scount(1+pgrid(myprow, pcol))
         new_rcount(1+pcol) = rcount(1+pgrid(myprow, pcol))
         new_sdispl(1+pcol) = sdispl(1+pgrid(myprow, pcol))
         new_rdispl(1+pcol) = rdispl(1+pgrid(myprow, pcol))
      END FORALL
      CALL mp_alltoall (sb, new_scount(1:npcols), new_sdispl(1:npcols),&
           rb, new_rcount(1:npcols), new_rdispl(1:npcols),&
           dbcsr_mp_my_row_group (mp_env))
      FORALL (prow = 0 : nprows-1)
         new_scount(1+prow) = scount(1+pgrid(prow, mypcol))
         new_rcount(1+prow) = rcount(1+pgrid(prow, mypcol))
         new_sdispl(1+prow) = sdispl(1+pgrid(prow, mypcol))
         new_rdispl(1+prow) = rdispl(1+pgrid(prow, mypcol))
      END FORALL
      CALL mp_alltoall (sb, new_scount(1:nprows), new_sdispl(1:nprows),&
           rb, new_rcount(1:nprows), new_rdispl(1:nprows),&
           dbcsr_mp_my_col_group (mp_env))
    END SUBROUTINE most_alltoall
    SUBROUTINE most_point_to_point ()
      ! Go through my prow and exchange.
      DO i = 0, npcols - 1
         pcol = MOD (mypcol+i, npcols)
         grp = dbcsr_mp_my_row_group (mp_env)
         !
         dst = dbcsr_mp_get_process (mp_env, myprow, pcol)
         send_cnt = scount(dst+1)
         send_data_p => sb( 1+sdispl(dst+1) : 1+sdispl(dst+1)+send_cnt-1 )
         tag = 4*mypcol
         IF (send_cnt .GT. 0) THEN
            CALL mp_isend (send_data_p, pcol, grp, row_sr(nrow_sr), tag)
            nrow_sr = nrow_sr+1
         ENDIF
         !
         pcol = MODULO (mypcol-i, npcols)
         src = dbcsr_mp_get_process (mp_env, myprow, pcol)
         recv_cnt = rcount(src+1)
         recv_data_p => rb( 1+rdispl(src+1) : 1+rdispl(src+1)+recv_cnt-1 )
         tag = 4*pcol
         IF (recv_cnt .GT. 0) THEN
            CALL mp_irecv (recv_data_p, pcol, grp, row_rr(nrow_rr), tag)
            nrow_rr = nrow_rr+1
         ENDIF
      ENDDO
      ! go through my pcol and exchange
      DO i = 0, nprows - 1
         prow = MOD (myprow+i, nprows)
         grp = dbcsr_mp_my_col_group (mp_env)
         !
         dst = dbcsr_mp_get_process (mp_env, prow, mypcol)
         send_cnt = scount(dst+1)
         send_data_p => sb( 1+sdispl(dst+1) : 1+sdispl(dst+1)+send_cnt-1 )
         tag = 4*myprow+1
         IF (send_cnt .GT. 0) THEN
            CALL mp_isend (send_data_p, prow, grp, col_sr(ncol_sr), tag)
            ncol_sr = ncol_sr + 1
         ENDIF
         !
         prow = MODULO (myprow-i, nprows)
         src = dbcsr_mp_get_process (mp_env, prow, mypcol)
         recv_cnt = rcount(src+1)
         recv_data_p => rb( 1+rdispl(src+1) : 1+rdispl(src+1)+recv_cnt-1 )
         tag = 4*prow+1
         IF (recv_cnt .GT. 0) THEN
            CALL mp_irecv (recv_data_p, prow, grp, col_rr(ncol_rr), tag)
            ncol_rr = ncol_rr + 1
         ENDIF
      ENDDO
    END SUBROUTINE most_point_to_point
    SUBROUTINE remainder_alltoall ()
      new_scount(:) = scount(:)
      new_rcount(:) = rcount(:)
      FORALL (prow = 0:nprows-1)
         new_scount(1+pgrid(prow, mypcol)) = 0
         new_rcount(1+pgrid(prow, mypcol)) = 0
      END FORALL
      FORALL (pcol = 0:npcols-1)
         new_scount(1+pgrid(myprow, pcol)) = 0
         new_rcount(1+pgrid(myprow, pcol)) = 0
      END FORALL
      CALL mp_alltoall (sb, new_scount, sdispl,&
           rb, new_rcount, rdispl, all_group)
    END SUBROUTINE remainder_alltoall
    SUBROUTINE remainder_point_to_point()
    INTEGER                                  :: col, row

      DO row = 0, nprows-1
         prow = MOD(row+myprow, nprows)
         IF (prow .EQ. myprow) CYCLE
         DO col = 0, npcols-1
            pcol = MOD (col+mypcol, npcols)
            IF (pcol .EQ. mypcol) CYCLE
            dst = dbcsr_mp_get_process (mp_env, prow, pcol)
            send_cnt = scount(dst+1)
            send_data_p => sb( 1+sdispl(dst+1) : 1+sdispl(dst+1)+send_cnt-1 )
            tag = 4*mynode+2
            IF (send_cnt .GT. 0) THEN
               CALL mp_isend (send_data_p, dst, all_group, all_sr(nall_sr+1), tag)
               nall_sr = nall_sr + 1
            ENDIF
            !
            src = dbcsr_mp_get_process (mp_env, prow, pcol)
            recv_cnt = rcount(src+1)
            recv_data_p => rb( 1+rdispl(src+1) : 1+rdispl(src+1)+recv_cnt-1 )
            tag = 4*src+2
            IF (recv_cnt .GT. 0) THEN
               CALL mp_irecv (recv_data_p, src, all_group, all_rr(nall_rr+1), tag)
               nall_rr = nall_rr+1
            ENDIF
         ENDDO
      ENDDO
    END SUBROUTINE remainder_point_to_point
  END SUBROUTINE hybrid_alltoall_s1

