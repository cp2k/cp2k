! *****************************************************************************
!> \brief Sets the diagonal of a DBCSR matrix.
!> \param[in] matrix          matrix in which to set the diagonal.
!> \param[in] value           value of the diagonal elements
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_set_diagonal_z (matrix, value, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    COMPLEX(kind=real_8), INTENT(IN)                      :: value
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_diagonal_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p, col, local_row, &
                                                my_pcol, nze, row, row_size
    INTEGER, DIMENSION(:), POINTER           :: local_rows, pcol_dist, &
                                                row_blk_size
    LOGICAL                                  :: in_this_row, valid_index
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: data_p

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL cp_assert(matrix%nblkrows_total .EQ. matrix%nblkcols_total,&
         cp_failure_level, cp_caller_error, routineP,&
         "Can not set diagonal of non-square matrix.", error=error)
    CALL cp_assert ('C' .EQ. matrix%data_type, cp_failure_level,&
         cp_caller_error, routineN, "Datatype mismatch.")
    IF (.NOT. ASSOCIATED (matrix%wms)) THEN
       CALL dbcsr_work_create (matrix, error=error)
    ENDIF
    row_blk_size => array_data (matrix%row_blk_size)
    CALL dbcsr_get_info (matrix, local_rows = local_rows,&
         proc_col_dist = pcol_dist,&
         my_pcol = my_pcol, error=error)
    valid_index = dbcsr_valid_index (matrix)
    ! Try to find the column element in each row. If it's not found,
    ! then add it.
    DO local_row = 1, matrix%nblkrows_local
       row = local_rows(local_row)
       ! Only diagonal elements come into play.
       IF (pcol_dist(row) .NE. my_pcol) CYCLE
       in_this_row = .FALSE.
       row_size = row_blk_size(row)
       nze = row_size * row_size
       IF (valid_index) THEN
          blk_loop: DO blk = matrix%row_p(row)+1, matrix%row_p(row+1)
             col = matrix%col_i(blk)
             IF (col .LT. row) CYCLE blk_loop
             IF (col .GT. row) EXIT
             IF (col .EQ. row) THEN
                ! Sets the diagonal in the existing block.
                blk_p = matrix%blk_p(blk)
                data_p => pointer_view (&
                        dbcsr_get_data_p (matrix%data_area, CMPLX(0.0, 0.0, real_8)),&
                        blk_p, blk_p+nze-1&
                     )
                CALL set_block_diagonal (data_p, value, row_size)
                !CALL set_block_diagonal (matrix%data(blk_p:blk_p+nze-1),&
                !     value, row_size)
                in_this_row = .TRUE.
                EXIT blk_loop ! Don't bother with the remaining columns
             ENDIF
          ENDDO blk_loop
       ENDIF ! valid_index
       ! Add the column if it doesn't already exist.
       IF (.NOT. in_this_row) THEN
          DBG 'Adding at row',row,'dimen',row_size,'size',nze,'pos',&
               matrix%wms(1)%datasize+1
          CALL add_work_coordinate (matrix%wms(1), row, row,&
               matrix%wms(1)%datasize+1,&
               error=error)
          ! add_work_coordinate increments the block counter
          CALL dbcsr_data_ensure_size (matrix%wms(1)%data_area,&
               matrix%wms(1)%datasize + nze)
          !CALL ensure_array_size (matrix%w%data, ub=matrix%w%datasize+nze,&
          !     error=error)
          data_p => pointer_view (&
               dbcsr_get_data_p (matrix%wms(1)%data_area, CMPLX(0.0, 0.0, real_8)),&
               matrix%wms(1)%datasize+1, matrix%wms(1)%datasize+nze)
          CALL set_block_diagonal(data_p, value, row_size)
          !CALL set_block_diagonal (&
          !     matrix%w%data(matrix%w%datasize+1 : matrix%w%datasize+nze),&
          !     value, row_size)
          matrix%wms(1)%datasize = matrix%wms(1)%datasize + nze
       ENDIF
    ENDDO ! row
  END SUBROUTINE dbcsr_set_diagonal_z

! *****************************************************************************
!> \brief Sets the diagonal of a square data block represented as a 1-D array.
!>        Non-diagonal elements are set to 0.
!> \param[out] block_data     sets diagonal in this data block
!> \param[in] value           value of the diagonal elements
!> \param[in] d               dimension of the square data block
! *****************************************************************************
  PURE SUBROUTINE set_block_diagonal_z (block_data, value, d)
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(OUT) :: block_data
    COMPLEX(kind=real_8), INTENT(IN)                :: value
    INTEGER, INTENT(IN)                :: d

    INTEGER                            :: i

!   ---------------------------------------------------------------------------

    block_data(:) = 0
    FORALL (i = 0 : d-1)
       block_data(1 + i*d + i) = value
    END FORALL
  END SUBROUTINE set_block_diagonal_z


! *****************************************************************************
!> \brief Returns a pointer with different bounds.
!> \param[in] original   original data pointer
!> \param[in] lb, ub     lower and upper bound for the new pointer view
!> \param[out] view      new pointer
! *****************************************************************************
  FUNCTION pointer_view_z (original, lb, ub) RESULT (view)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: original, view
    INTEGER, INTENT(IN)                  :: lb, ub
    view => original(lb:ub)
  END FUNCTION pointer_view_z


! *****************************************************************************
!> \brief Copy and transpose block.
!> \param[out] extent_out     output matrix in the form of a 1-d array
!> \param[in] extent_in       input matrix in the form of a 1-d array
!> \param[in] rows, columns   input matrix size
!> \param error     cp2k error
! *****************************************************************************
  PURE SUBROUTINE block_transpose_copy_z(extent_out, extent_in,&
       rows, columns)
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(OUT) :: extent_out
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN)  :: extent_in
    INTEGER, INTENT(IN)                :: rows, columns

    CHARACTER(len=*), PARAMETER :: routineN = 'block_transpose_copy_z', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    extent_out = RESHAPE(TRANSPOSE(&
         RESHAPE(extent_in, (/rows, columns/))), (/rows*columns/))
  END SUBROUTINE block_transpose_copy_z


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
  SUBROUTINE dbcsr_data_set_az (dst, lb, data_size, src, source_lb)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    INTEGER, INTENT(IN)                      :: lb, data_size
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN)        :: src
    INTEGER, INTENT(IN), OPTIONAL            :: source_lb
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_set_az', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: lb_s, ub, ub_s
!   ---------------------------------------------------------------------------
    IF (bcsr_debug) THEN
       CALL cp_assert (ASSOCIATED (dst%d),&
            cp_fatal_level, cp_caller_error, routineN,&
            "Target data area must be setup.")
       CALL cp_assert (SIZE(src) .GE. data_size,&
            cp_fatal_level, cp_caller_error, routineN,&
            "Not enough source data.")
       CALL cp_assert (dst%d%data_type .EQ. 'C', cp_failure_level,&
            cp_caller_error, routineN, "Data type mismatch.")
    ENDIF
    ub = lb + data_size - 1
    IF (PRESENT (source_lb)) THEN
       lb_s = source_lb
       ub_s = source_lb + data_size-1
    ELSE
       lb_s = lb
       ub_s = ub
    ENDIF
    dst%d%c_dp(lb:ub) = src(lb_s:ub_s)
  END SUBROUTINE dbcsr_data_set_az


! *****************************************************************************
!> \brief Ensures that an array is appropriately large.
!> \param[in,out] array       array to verify and possibly resize
!> \param[in] lb    (optional) desired array lower bound
!> \param[in] ub    desired array upper bound
!> \param[in] factor          (optional) factor by which to exagerrate
!>                            enlargements
!> \param[in] nocopy          (optional) copy array on enlargement; default
!>                            is to copy
! *****************************************************************************
  SUBROUTINE ensure_array_size_z(array, lb, ub, factor, nocopy, special)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: array
    INTEGER, INTENT(IN), OPTIONAL            :: lb
    INTEGER, INTENT(IN)                      :: ub
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, special

    CHARACTER(len=*), PARAMETER :: routineN = 'ensure_array_size_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lb_new, lb_orig, stat, &
                                                ub_new, ub_orig
    LOGICAL                                  :: dbg, docopy = .TRUE., &
                                                failure = .TRUE., info, spec
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: newarray

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    info = .FALSE.
    CALL dbcsr_set_debug(dbg, info=info)

    docopy = .TRUE.
    IF (PRESENT (nocopy)) docopy = .NOT. nocopy
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    lb_new = 1
    IF (PRESENT (lb)) lb_new = lb
    !> Creates a new array if it doesn't yet exist.
    IF (.NOT.ASSOCIATED(array)) THEN
       failure = .FALSE.
       IF (spec) THEN
          CALL cp_assert (lb.EQ.1, cp_warning_level, cp_internal_error,&
               routineN, "Special memory can only start at 1.")
          CALL mp_allocate(array, ub-lb_new+1, stat=stat)
       ELSE
          ALLOCATE(array(lb_new:ub), stat=stat)
       ENDIF
       IF (stat /= 0) CALL stop_memory(routineP, 'array', ub-lb+1)
       RETURN
    ENDIF
    lb_orig = LBOUND(array,1)
    ub_orig = UBOUND(array,1)
    IF(dbg) WRITE(*,*)routineP//' Current bounds are',lb_orig,':',ub_orig,&
         '; special?',spec
    ! The existing array is big enough.
    IF (lb_orig.LE.lb_new .AND. ub_orig.GE.ub) THEN
       RETURN
    ENDIF
    ! A reallocation must be performed
    !CALL timeset(routineN,timing_handle)
    IF (lb_orig.GT.lb_new) THEN
       IF (PRESENT(factor)) THEN
          lb_new = MIN(lb_orig, lb_new - NINT(REAL(lb_orig - lb_new)*factor))
       ENDIF
    ENDIF
    IF (ub_orig.LT.ub) THEN
       IF (PRESENT(factor)) THEN
          ub_new = MAX(ub_orig, ub + NINT(REAL(ub-ub_orig)*factor))
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
       CALL cp_assert (lb_new.EQ.1, cp_warning_level, cp_internal_error,&
            routineN, "Special memory can only start at 1.")
       CALL mp_allocate(newarray, ub_new-lb_new+1, stat=stat)
    ELSE
       ALLOCATE(newarray(lb_new:ub_new), stat=stat)
    ENDIF
    IF(stat /= 0) CALL stop_memory(routineP, 'newarray', ub_new-lb_new+1)
    IF(docopy) THEN
       IF(dbg) CALL cp_assert(lb_new.LE.lb_orig .AND. ub_new.GE.ub_orig,&
            cp_failure_level, cp_internal_error, routineP,&
            "Old extent exceeds the new one.")
       newarray(lb_orig:ub_orig) = array(lb_orig:ub_orig)
       IF (spec) THEN
          CALL mp_deallocate(array)
       ELSE
          DEALLOCATE(array)
       ENDIF
    END IF
    array => newarray
    IF(dbg.AND.bcsr_verbose)&
         WRITE(*,*)routineP//' New array',array(lb_orig:ub_orig)
    IF (dbg) WRITE(*,*)routineP//' New array size', SIZE(array)
  END SUBROUTINE ensure_array_size_z


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
  SUBROUTINE addto_array_az(to, from, offset, factor)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: to
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN)        :: from
    INTEGER, INTENT(IN), OPTIONAL            :: offset
    REAL, INTENT(IN), OPTIONAL               :: factor

    CHARACTER(len=*), PARAMETER :: routineN = 'addto_array_az', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: last, lb_t, s_f, s_t, start, &
                                                timing_handle, ub_t
    LOGICAL                                  :: dbg, spec
    REAL                                     :: f

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    spec = to%d%special
    s_f = SIZE(from, 1)
    IF (dbg) THEN
       WRITE(*,*)routineP//' from size',SIZE(from),';', SIZE(from,1)
       WRITE(*,*)routineP//' from lbound',LBOUND(from),';', LBOUND(from,1)
       WRITE(*,*)routineP//' from ubound',UBOUND(from),';', UBOUND(from,1)
    ENDIF
    IF (.NOT. ASSOCIATED (to%d%c_dp)) THEN
       IF (PRESENT (offset)) THEN
          start = offset
          last = offset + SIZE(from)-1
          CALL ensure_array_size(to%d%c_dp, lb=start,&
               ub=last,&
               factor=1.0, nocopy=.TRUE., special=spec)
       ELSE
          start = 1
          last = SIZE(from)
          CALL ensure_array_size(to%d%c_dp, ub=last,&
               factor=1.0, nocopy=.TRUE., special=spec)          
       ENDIF
    ELSE
       ub_t = UBOUND(to%d%c_dp, 1)
       lb_t = LBOUND(to%d%c_dp, 1)
       s_t = SIZE(to%d%c_dp)
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
               UBOUND(to%d%c_dp),&
               'to',ub_t
          !CALL reallocate(to, lb_t, ub_t)
          CALL ensure_array_size(to%d%c_dp, lb_t, ub_t,&
               factor=1.0, nocopy=.FALSE.,&
               special=spec)
       ENDIF
    ENDIF
    to%d%c_dp(start:last) = from(:)
  END SUBROUTINE addto_array_az

