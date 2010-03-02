!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
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
                                             dbcsr_data_set
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: real_4,&
                                             real_8
  USE dbcsr_methods,                   ONLY: dbcsr_data_hold,&
                                             dbcsr_get_data,&
                                             dbcsr_get_data_size
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE dbcsr_util,                      ONLY: dbcsr_set_debug,&
                                             sgn
  USE f77_blas
  USE f77_blas_generic

  !$ USE OMP_LIB

  IMPLICIT NONE


  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_data_operations'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  PUBLIC :: dbcsr_data_copyall, dbcsr_data_convert,&
            dbcsr_copy_sort_data,&
            dbcsr_sort_data

#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL dbcsr_set_debug(dbg, info=info)

  LOGICAL, PARAMETER :: bcsr_debug =   .TRUE.
  LOGICAL, PARAMETER :: bcsr_info =    .FALSE.
  LOGICAL, PARAMETER :: bcsr_verbose = .FALSE.

CONTAINS



! *****************************************************************************
!> \brief Copies a data area, deep by default.
!> \param[inout] target_area     target data area
!> \param[in]    source_area     source data area
!> \param[in]    shallow      shallow copy (default is deep)
! *****************************************************************************
  SUBROUTINE dbcsr_data_copyall (target_area, source_area, shallow)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: target_area
    TYPE(dbcsr_data_obj), INTENT(IN)         :: source_area
    LOGICAL, INTENT(IN), OPTIONAL            :: shallow

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_copyall', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: s_data_c, t_data_c
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: s_data_z, t_data_z
    INTEGER                                  :: n
    LOGICAL                                  :: shallow_copy
    REAL(KIND=real_4), DIMENSION(:), POINTER :: s_data_r, t_data_r
    REAL(KIND=real_8), DIMENSION(:), POINTER :: s_data_d, t_data_d
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (source_area%d), dbcsr_warning_level,&
         dbcsr_caller_error, routineN, "Attempt to copy unassigned data",__LINE__,error)
    CALL dbcsr_assert (source_area%d%refcount .GT. 0, dbcsr_warning_level,&
         dbcsr_caller_error, routineN, "Attempt to copy unheld data",__LINE__,error)
    shallow_copy = .FALSE.
    IF (PRESENT (shallow)) shallow_copy = shallow
    IF (shallow_copy) THEN
       target_area = source_area
       CALL dbcsr_data_hold (target_area)
    ELSE
       IF (.NOT. ASSOCIATED (source_area%d)) THEN
          RETURN
       ENDIF
       SELECT CASE (source_area%d%data_type)
          CASE (dbcsr_type_real_8)
             n = dbcsr_get_data_size (source_area)
             CALL dbcsr_get_data (target_area, t_data_d)
             CALL dbcsr_get_data (source_area, s_data_d)
             IF(n.GT.0)t_data_d(1:n) = s_data_d(1:n)
          CASE (dbcsr_type_real_4)
             n = dbcsr_get_data_size (source_area)
             CALL dbcsr_get_data (target_area, t_data_r)
             CALL dbcsr_get_data (source_area, s_data_r)
             IF(n.GT.0)t_data_r(1:n) = s_data_r(1:n)
          CASE (dbcsr_type_complex_8)
             n = dbcsr_get_data_size (source_area)
             CALL dbcsr_get_data (target_area, t_data_z)
             CALL dbcsr_get_data (source_area, s_data_z)
             IF(n.GT.0)t_data_z(1:n) = s_data_z(1:n)
          CASE (dbcsr_type_complex_4)
             n = dbcsr_get_data_size (source_area)
             CALL dbcsr_get_data (target_area, t_data_c)
             CALL dbcsr_get_data (source_area, s_data_c)
             IF(n.GT.0)t_data_c(1:n) = s_data_c(1:n)
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error, &
                  routineN, "Invalid data type",__LINE__,error)
          END SELECT
    ENDIF
    ! If we're releasing the last reference, then free the memory.
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
    noimult = .FALSE.
    IF (PRESENT (multiply_by_i)) noimult = .NOT. multiply_by_i
    ns = dbcsr_get_data_size (source_area)
    nt = dbcsr_get_data_size (source_area)
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
                t_data_c(1:n) = CMPLX(s_data_d(1:n), KIND=real_8)
             ELSE
                t_data_c(1:n) = CMPLX(0.0, s_data_d(1:n), KIND=real_8)
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
                t_data_c(1:n) = CMPLX(s_data_r(1:n), KIND=real_8)
             ELSE
                t_data_c(1:n) = CMPLX(0.0, s_data_r(1:n), KIND=real_8)
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
             t_data_c(1:n) = CMPLX(s_data_z(1:n), KIND=real_8)
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
    INTEGER, DIMENSION(:), INTENT(OUT)       :: blk_p
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

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
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
                  ABS(blk_p(blk))+nze-1,'LE', dbcsr_get_data_size(dst),&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  'Destination data space is too small.',__LINE__,error)
             IF (.NOT. trb) THEN
                CALL dbcsr_data_set (dst, ABS(blk_p(blk)), nze,&
                     src, source_lb=ABS(old_blk_p(blk)))
             ELSE
                CALL dbcsr_block_transpose (dst, src,&
                     col_size, row_size,&
                     lb=ABS(blk_p(blk)), source_lb=ABS(old_blk_p(blk)))
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
    INTEGER, DIMENSION(:), INTENT(OUT)       :: blk_p
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

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    multidata = PRESENT (srcs) .AND. PRESENT (old_blk_d)
    nblks = SIZE(old_blk_p)
    DBG 'blk sizes', SIZE(blk_p), nblks
    DBG 'dst data sizes', dbcsr_get_data_size(dsts)
    CALL dbcsr_assert(SIZE(blk_p),'GE',nblks, dbcsr_failure_level,&
         dbcsr_wrong_args_error, routineN, 'Destination blk_p too small.',__LINE__,error)
    ! If the source is deliberately bigger, we don't want any warnings.
    !CALL dbcsr_assert(SIZE(dst).GE.SIZE(src), dbcsr_warning_level,&
    !     dbcsr_wrong_args_error, routineN,&
    !     'Destination data space may be too small.',&
    !     error)
    DBGV 'old_blk_p',old_blk_p
    DBGV 'sizes', sizes
    IF (SIZE (old_blk_p) .GT. 0) THEN
       blk_p(1) = SIGN(1, old_blk_p(1))
       DO i = 2, nblks
          blk_p(i) = SIGN(ABS(blk_p(i-1)) + sizes(i-1), old_blk_p(i))
       ENDDO
       !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i)
       DO i = 1, nblks
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
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF
    DBGV 'new blk_p =',blk_p
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_sort_data


END MODULE dbcsr_data_operations
