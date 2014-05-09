!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR data operations
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - 2010-02-18 Moved from dbcsr_util
! *****************************************************************************
MODULE dbcsr_data_operations

  USE dbcsr_block_operations,          ONLY: dbcsr_block_transpose,&
                                             dbcsr_data_copy,&
                                             dbcsr_data_set
  USE dbcsr_data_methods,              ONLY: dbcsr_data_get_size,&
                                             dbcsr_data_get_size_referenced,&
                                             dbcsr_data_hold,&
                                             dbcsr_data_set_size_referenced,&
                                             dbcsr_get_data
  USE dbcsr_error_handling
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE dbcsr_util,                      ONLY: sgn
  USE kinds,                           ONLY: real_4,&
                                             real_8

  !$ USE OMP_LIB

  IMPLICIT NONE


  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_data_operations'


  PUBLIC :: dbcsr_data_copyall, dbcsr_data_convert,&
            dbcsr_copy_sort_data,&
            dbcsr_sort_data


CONTAINS



! *****************************************************************************
!> \brief Copies a data area, deep by default.
!> \param[in,out] target_area     target data area
!> \param[in]     source_area     source data area
!> \param[in]     shallow         shallow copy (default is deep)
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_data_copyall (target_area, source_area, shallow, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: target_area
    TYPE(dbcsr_data_obj), INTENT(IN)         :: source_area
    LOGICAL, INTENT(IN), OPTIONAL            :: shallow
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_copyall', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, i, n
    LOGICAL                                  :: shallow_copy

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    CALL dbcsr_assert (ASSOCIATED (source_area%d), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Attempt to copy unassigned data",__LINE__,error)
    CALL dbcsr_assert (source_area%d%refcount .GT. 0, dbcsr_warning_level,&
         dbcsr_caller_error, routineN, "Attempt to copy unheld data",__LINE__,error)
    shallow_copy = .FALSE.
    IF (PRESENT (shallow)) shallow_copy = shallow
    IF (shallow_copy) THEN
       target_area = source_area
       CALL dbcsr_data_hold (target_area)
    ELSE
       CALL dbcsr_assert (ASSOCIATED (target_area%d), dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "Target area does not exist.",&
            __LINE__, error=error)
       CALL dbcsr_data_set_size_referenced(target_area,&
            dbcsr_data_get_size_referenced(source_area))
       n = dbcsr_data_get_size_referenced (source_area)
       SELECT CASE (target_area%d%data_type)
       CASE (dbcsr_type_real_4)
          !$OMP PARALLEL DO
          DO i=1,n
             target_area%d%r_sp(i)=source_area%d%r_sp(i)
          ENDDO
       CASE (dbcsr_type_real_8)
          !$OMP PARALLEL DO
          DO i=1,n
             target_area%d%r_dp(i)=source_area%d%r_dp(i)
          ENDDO
       CASE (dbcsr_type_complex_4)
          !$OMP PARALLEL DO
          DO i=1,n
             target_area%d%c_sp(i)=source_area%d%c_sp(i)
          ENDDO
       CASE (dbcsr_type_complex_8)
          !$OMP PARALLEL DO
          DO i=1,n
             target_area%d%c_dp(i)=source_area%d%c_dp(i)
          ENDDO
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_data_copyall

! *****************************************************************************
!> \brief Copies a data area, converting data type
!> \param[inout] target_area     target data area
!> \param[in]    source_area     source data area
!> \param[in]    drop_real       (optional) drops real part of complex
!>                               numbers instead of the imaginary part; default
!>                               is false
!> \param[in]    multiply_by_i   (optional) converts real to complex by placing
!>                               into imaginary instead of real part
! *****************************************************************************
  SUBROUTINE dbcsr_data_convert (target_area, source_area, drop_real,&
       multiply_by_i)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: target_area
    TYPE(dbcsr_data_obj), INTENT(IN)         :: source_area
    LOGICAL, INTENT(IN), OPTIONAL            :: drop_real, multiply_by_i

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_convert', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: s_data_c, t_data_c
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: s_data_z, t_data_z
    INTEGER                                  :: n, ns, nt
    LOGICAL                                  :: keep_real, noimult
    REAL(KIND=real_4), DIMENSION(:), POINTER :: s_data_r, t_data_r
    REAL(KIND=real_8), DIMENSION(:), POINTER :: s_data_d, t_data_d
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (source_area%d), dbcsr_warning_level,&
         dbcsr_caller_error, routineN, "Attempt to copy unassigned data",__LINE__,error)
    CALL dbcsr_assert (source_area%d%refcount .GT. 0, dbcsr_warning_level,&
         dbcsr_caller_error, routineN, "Attempt to copy unheld data",__LINE__,error)
    IF (.NOT. ASSOCIATED (source_area%d)) THEN
       RETURN
    ENDIF
    keep_real = .TRUE.
    IF (PRESENT (drop_real)) keep_real = .NOT. drop_real
    noimult = .TRUE. 
    IF (PRESENT (multiply_by_i)) noimult = .NOT. multiply_by_i
    ns = dbcsr_data_get_size (source_area)
    nt = dbcsr_data_get_size (source_area)
    n = MIN (ns, nt)
    IF (n .GT. 0) THEN
       SELECT CASE (source_area%d%data_type)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_get_data (source_area, s_data_d)
          SELECT CASE (target_area%d%data_type)
          CASE(dbcsr_type_real_8)
             CALL dbcsr_get_data (target_area, t_data_d)
             t_data_d(1:n) = s_data_d(1:n)
          CASE(dbcsr_type_real_4)
             CALL dbcsr_get_data (target_area, t_data_r)
             t_data_r(1:n) = REAL(s_data_d(1:n), KIND=real_4)
          CASE(dbcsr_type_complex_8)
             CALL dbcsr_get_data (target_area, t_data_z)
             IF (noimult) THEN
                t_data_z(1:n) = CMPLX(s_data_d(1:n), KIND=real_8)
             ELSE
                t_data_z(1:n) = CMPLX(0.0, s_data_d(1:n), KIND=real_8)
             ENDIF
          CASE(dbcsr_type_complex_4)
             CALL dbcsr_get_data (target_area, t_data_c)
             IF (noimult) THEN
                t_data_c(1:n) = CMPLX(s_data_d(1:n), KIND=real_4)
             ELSE
                t_data_c(1:n) = CMPLX(0.0, s_data_d(1:n), KIND=real_4)
             ENDIF
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error, &
                  routineN, "Invalid data type",__LINE__,error)
          END SELECT
       CASE (dbcsr_type_real_4)
          CALL dbcsr_get_data (source_area, s_data_r)
          SELECT CASE (target_area%d%data_type)
          CASE(dbcsr_type_real_8)
             CALL dbcsr_get_data (target_area, t_data_d)
             t_data_d(1:n) = REAL(s_data_r(1:n), KIND=real_8)
          CASE(dbcsr_type_real_4)
             CALL dbcsr_get_data (target_area, t_data_r)
             t_data_r(1:n) = s_data_r(1:n)
          CASE(dbcsr_type_complex_8)
             CALL dbcsr_get_data (target_area, t_data_z)
             IF (noimult) THEN
                t_data_z(1:n) = CMPLX(s_data_r(1:n), KIND=real_8)
             ELSE
                t_data_z(1:n) = CMPLX(0.0, s_data_r(1:n), KIND=real_8)
             ENDIF
          CASE(dbcsr_type_complex_4)
             CALL dbcsr_get_data (target_area, t_data_c)
             IF (noimult) THEN
                t_data_c(1:n) = CMPLX(s_data_r(1:n), KIND=real_4)
             ELSE
                t_data_c(1:n) = CMPLX(0.0, s_data_r(1:n), KIND=real_4)
             ENDIF
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error, &
                  routineN, "Invalid data type",__LINE__,error)
          END SELECT
       CASE (dbcsr_type_complex_8)
          CALL dbcsr_get_data (source_area, s_data_z)
          SELECT CASE (target_area%d%data_type)
          CASE(dbcsr_type_real_8)
             CALL dbcsr_get_data (target_area, t_data_d)
             IF (keep_real) THEN
                t_data_d(1:n) = REAL(s_data_z(1:n), KIND=real_8)
             ELSE
                t_data_d(1:n) = AIMAG(s_data_z(1:n))
             ENDIF
          CASE(dbcsr_type_real_4)
             CALL dbcsr_get_data (target_area, t_data_r)
             IF (keep_real) THEN
                t_data_r(1:n) = REAL(s_data_z(1:n), KIND=real_4)
             ELSE
                t_data_r(1:n) = REAL(AIMAG (s_data_z(1:n)), KIND=real_4)
             ENDIF
          CASE(dbcsr_type_complex_8)
             CALL dbcsr_get_data (target_area, t_data_z)
             t_data_z(1:n) = s_data_z(1:n)
          CASE(dbcsr_type_complex_4)
             CALL dbcsr_get_data (target_area, t_data_c)
             t_data_c(1:n) = CMPLX(s_data_z(1:n), KIND=real_4)
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error, &
                  routineN, "Invalid data type",__LINE__,error)
          END SELECT
       CASE (dbcsr_type_complex_4)
          CALL dbcsr_get_data (source_area, s_data_c)
          SELECT CASE (target_area%d%data_type)
          CASE(dbcsr_type_real_8)
             CALL dbcsr_get_data (target_area, t_data_d)
             IF (keep_real) THEN
                t_data_d(1:n) = REAL(s_data_c(1:n), KIND=real_8)
             ELSE
                t_data_d(1:n) = REAL(AIMAG(s_data_c(1:n)), KIND=real_8)
             ENDIF
          CASE(dbcsr_type_real_4)
             CALL dbcsr_get_data (target_area, t_data_r)
             IF (keep_real) THEN
                t_data_r(1:n) = REAL(s_data_c(1:n), KIND=real_4)
             ELSE
                t_data_r(1:n) = AIMAG(s_data_c(1:n))
             ENDIF
          CASE(dbcsr_type_complex_8)
             CALL dbcsr_get_data (target_area, t_data_z)
             t_data_z(1:n) = CMPLX(s_data_c(1:n), KIND=real_8)
          CASE(dbcsr_type_complex_4)
             CALL dbcsr_get_data (target_area, t_data_c)
             t_data_c(1:n) = s_data_c(1:n)
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error, &
                  routineN, "Invalid data type",__LINE__,error)
          END SELECT
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error, &
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    ENDIF
  END SUBROUTINE dbcsr_data_convert



! *****************************************************************************
!> \brief Sorts the data in a matrix so that the data blocks follow
!>        sequentially and does various transposing options.
!>
!> As opposed to dbcsr_sort_data, this routine calculates block sizes
!> \param[out] blk_p          re-arranged block pointers reflecting the new
!>                            data order
!> \param[in] old_blk_p       current block pointers
!> \param[in] rbs, cbs        sizes of the blocked rows and columns
!> \param[in] row_p, col_i    index
!> \param[out] dst            sorted data
!> \param[in] src             existing unordered data
!> \param[in] mark_transposed           mark data as transposed by negating the
!>                                      blk_p index entries
!> \param[in] transpose_blocks          transpose data blocks
! *****************************************************************************
  SUBROUTINE dbcsr_copy_sort_data(blk_p, old_blk_p, row_p, col_i, rbs, cbs,&
       dst, src, mark_transposed, transpose_blocks)
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: blk_p
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_blk_p, row_p, col_i, rbs, &
                                                cbs
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    TYPE(dbcsr_data_obj), INTENT(IN)         :: src
    LOGICAL, INTENT(IN), OPTIONAL            :: mark_transposed, &
                                                transpose_blocks

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_sort_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col_size, nblks, nrows, &
                                                nze, nze_prev, row, row_size
    LOGICAL                                  :: mark, trb
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------
! Analyze parameters

    mark = .FALSE.
    IF (PRESENT (mark_transposed)) mark = mark_transposed
    trb = .FALSE.
    IF (PRESENT (transpose_blocks)) trb = transpose_blocks
    !
    nblks = SIZE (old_blk_p)
    nrows = SIZE (row_p)-1
    CALL dbcsr_assert(SIZE(blk_p),'GE',nblks, dbcsr_failure_level,&
         dbcsr_wrong_args_error, routineN, 'Destination blk_p too small.',__LINE__,error)
    IF (nblks .GE. 1) &
         blk_p(1) = SGN(1, old_blk_p(1), mark)
    nze_prev = 0
    DO row = 1, nrows
       row_size = rbs (row)
       DO blk = row_p(row)+1, row_p(row+1)
          IF (old_blk_p(blk) .NE. 0) THEN
             col_size = cbs(col_i(blk))
             nze = row_size * col_size
             IF (blk .GT. 1) THEN
                blk_p(blk) = SGN(ABS(blk_p(blk-1)) + nze_prev, old_blk_p(blk),&
                     mark)
             ENDIF
             CALL dbcsr_assert(&
                  ABS(blk_p(blk))+nze-1,'LE', dbcsr_data_get_size(dst),&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  'Destination data space is too small.',__LINE__,error)
             IF (.NOT. trb) THEN
                CALL dbcsr_data_copy (dst=dst, dst_lb=(/ABS(blk_p(blk))/),&
                     dst_sizes = (/ nze /),&
                     src=src, src_lb=(/ ABS(old_blk_p(blk)) /),&
                     src_sizes=(/ nze /), error=error)
                !CALL dbcsr_data_set (dst, ABS(blk_p(blk)), nze,&
                !     src, source_lb=ABS(old_blk_p(blk)))
             ELSE
                CALL dbcsr_block_transpose (dst, src,&
                     col_size, row_size,&
                     lb=ABS(blk_p(blk)), source_lb=ABS(old_blk_p(blk)), &
                     error=error)
             ENDIF
             nze_prev = nze
          ENDIF ! blk exists
       ENDDO ! blk
    ENDDO ! row
  END SUBROUTINE dbcsr_copy_sort_data



! *****************************************************************************
!> \brief Sorts the data in a matrix so that the data blocks follow
!>        sequentially.
!> \param[out] blk_p          re-arranged block pointers reflecting the new
!>                            data order
!> \param[in] old_blk_p       current block pointers
!> \param[in] sizes           sizes of the data blocks
!> \param[out] dst            sorted data
!> \param[in] src             existing unordered data
!> \param[in] srcs            (optional) multiple source areas
!> \param[in] blk_d           (optional) mapping from blocks to source areas
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_sort_data(blk_p, old_blk_p, sizes, dsts, src,&
       srcs, old_blk_d, error)
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: blk_p
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_blk_p, sizes
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dsts
    TYPE(dbcsr_data_obj), INTENT(IN)         :: src
    TYPE(dbcsr_data_obj), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: srcs
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: old_blk_d
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sort_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, i, nblks
    LOGICAL                                  :: multidata

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    multidata = PRESENT (srcs) .AND. PRESENT (old_blk_d)
    nblks = SIZE(old_blk_p)
    IF (nblks .GT. 0) THEN
       !$OMP BARRIER
       !$OMP MASTER
       blk_p(1) = SIGN(1, old_blk_p(1))
       DO i = 2, nblks
          blk_p(i) = SIGN(ABS(blk_p(i-1)) + sizes(i-1), old_blk_p(i))
       ENDDO
       CALL dbcsr_data_set_size_referenced (dsts,&
            ABS(blk_p(nblks)) + sizes(nblks)-1)
       !$OMP END MASTER
       !$OMP BARRIER
       !$OMP DO
       DO i = 1, nblks
          IF (old_blk_p(i) .NE. 0) THEN
             IF (.NOT. multidata) THEN
                CALL dbcsr_data_set (dsts,&
                     ABS(blk_p(i)), sizes(i),&
                     src, source_lb=ABS(old_blk_p(i)))
                !dst(ABS(blk_p(i)):ABS(blk_p(i))+sizes(i)-1) =&
                !     src(ABS(old_blk_p(i)):ABS(old_blk_p(i))+sizes(i)-1)
             ELSE
                CALL dbcsr_data_set (dsts,&
                     ABS(blk_p(i)), sizes(i),&
                     srcs(old_blk_d(i)), source_lb=ABS(old_blk_p(i)))
                !dst(ABS(blk_p(i)):ABS(blk_p(i))+sizes(i)-1) =&
                !     srcs(old_blk_d(i))%d&
                !     %r_dp(ABS(old_blk_p(i)):ABS(old_blk_p(i))+sizes(i)-1)
             ENDIF
          ENDIF
       ENDDO
       !$OMP END DO NOWAIT
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_sort_data




END MODULE dbcsr_data_operations
