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
