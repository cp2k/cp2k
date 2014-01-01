!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Operations with block buffers
!> \author  Urban Borstnik
!> \date    2010-04-22
!> \version 0.95
!>
!> <b>Modification history:</b>
!> - Created
! *****************************************************************************
MODULE dbcsr_block_buffers
  USE dbcsr_block_operations,          ONLY: dbcsr_block_partial_copy
  USE dbcsr_data_methods,              ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_get_sizes, dbcsr_data_get_type, &
       dbcsr_data_hold, dbcsr_data_init, dbcsr_data_new, dbcsr_data_release, &
       dbcsr_data_resize, dbcsr_type_1d_to_2d, dbcsr_type_is_2d
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: real_4,&
                                             real_8
  USE dbcsr_ptr_util,                  ONLY: dbcsr_ptr_remapping
  USE dbcsr_types,                     ONLY: dbcsr_block_buffer_obj,&
                                             dbcsr_data_obj,&
                                             dbcsr_type_complex_4_2d,&
                                             dbcsr_type_complex_8_2d,&
                                             dbcsr_type_real_4_2d,&
                                             dbcsr_type_real_8_2d

  !$ USE OMP_LIB
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_block_buffers'

  LOGICAL, PARAMETER, PRIVATE :: debug_mod = .FALSE., careful_mod = .FALSE.

  LOGICAL, PARAMETER :: dbcsr_buffers_2d_needed = .NOT. dbcsr_ptr_remapping

  PUBLIC :: dbcsr_buffers_init, dbcsr_buffers_valid, dbcsr_buffers_new,&
            dbcsr_buffers_hold, dbcsr_buffers_release

  PUBLIC :: dbcsr_buffers_are_dirty, dbcsr_buffers_is_dirty, &
            dbcsr_buffers_flush, dbcsr_buffers_mark_dirty

  PUBLIC :: dbcsr_buffers_resize
  PUBLIC :: dbcsr_buffers_set_pointer_2d

  PUBLIC :: dbcsr_buffers_2d_needed

  ! For dealing with OpenMP
  PUBLIC :: dbcsr_buffers_get_num, dbcsr_buffers_which_is_my
  PUBLIC :: dbcsr_buffers_get_area


  INTERFACE dbcsr_buffers_set_pointer_2d
     MODULE PROCEDURE dbcsr_buffers_set_p_2d_any
     MODULE PROCEDURE dbcsr_buffers_set_p_2d_s,&
                      dbcsr_buffers_set_p_2d_d,&
                      dbcsr_buffers_set_p_2d_c,&
                      dbcsr_buffers_set_p_2d_z
  END INTERFACE
  !INTERFACE dbcsr_buffer_make_size
  !   MODULE PROCEDURE dbcsr_destroy_2d_array, dbcsr_destroy_1d_array
  !END INTERFACE

CONTAINS


! *****************************************************************************
!> \brief Initializes a buffer object
!> \param[out] buffers            Buffers to create
! *****************************************************************************
  SUBROUTINE dbcsr_buffers_init (buffers)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(OUT)                            :: buffers

    NULLIFY (buffers%b)
  END SUBROUTINE dbcsr_buffers_init


! *****************************************************************************
!> \brief Verifies whether buffers exist
!> \param[in] buffers            Buffers to verify
!> \result exist                 Buffers exist
! *****************************************************************************
  FUNCTION dbcsr_buffers_valid (buffers) RESULT (valid)
    TYPE(dbcsr_block_buffer_obj), INTENT(IN) :: buffers
    LOGICAL                                  :: valid

    valid = ASSOCIATED (buffers%b)
  END FUNCTION dbcsr_buffers_valid

! *****************************************************************************
!> \brief Creates a buffer object.
!>
!>        Allocates the buffer array for the maximum number of threads.
!> \param[out] buffers            Buffers to create
!> \param[in] data_backing        Data backed by the buffers
!> \param[in,out] error           error
! *****************************************************************************
  SUBROUTINE dbcsr_buffers_new (buffers, data_backing, error)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    TYPE(dbcsr_data_obj), INTENT(IN)         :: data_backing
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, i, n, stat

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
!$OMP BARRIER
!$OMP MASTER
    ALLOCATE (buffers%b)
    n = 1
!$  n = omp_get_max_threads()
    buffers%b%data_type = dbcsr_data_get_type (data_backing)
    buffers%b%main = data_backing
    CALL dbcsr_data_hold (buffers%b%main)
    buffers%b%refcount = 1
    ALLOCATE (buffers%b%dirty(n), stat=stat)
    buffers%b%dirty(:) = .FALSE.
    ALLOCATE (buffers%b%rcb(10,n))
    buffers%b%rcb(:,:) = 0
    ALLOCATE (buffers%b%buffers(n))
    ALLOCATE (buffers%b%backing(n))
    DO i = 1, n
       CALL dbcsr_data_init (buffers%b%buffers(i))
       CALL dbcsr_data_new (buffers%b%buffers(i),&
            dbcsr_type_1d_to_2d(dbcsr_data_get_type (data_backing)))
       CALL dbcsr_data_init (buffers%b%backing(i))
       CALL dbcsr_data_new (buffers%b%backing(i),&
            dbcsr_type_1d_to_2d(dbcsr_data_get_type (data_backing)))
    ENDDO
!$OMP END MASTER
!$OMP BARRIER
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_buffers_new


! *****************************************************************************
!> \brief Registers another reference to a buffer object
!> \param[in,out] buffers          Buffers for which to register another
!>                                 reference
!> \param[in,out] error            error
! *****************************************************************************
  SUBROUTINE dbcsr_buffers_hold (buffers, error)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_hold', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    IF (careful_mod) &
         CALL dbcsr_assert (ASSOCIATED (buffers%b), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Buffers not valid", __LINE__,&
         error=error)
    buffers%b%refcount = buffers%b%refcount + 1
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_buffers_hold


! *****************************************************************************
!> \brief Releases a reference to a buffer object and cleans up if it is the
!>        final release.
!> \param[in,out] buffers         Buffers to release
!> \param[in,out] error            error
! *****************************************************************************
  SUBROUTINE dbcsr_buffers_release (buffers, error)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, i, n

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    !$OMP BARRIER
    !$OMP MASTER
    CALL dbcsr_assert (ASSOCIATED (buffers%b), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Buffers not valid", __LINE__,&
         error=error)
    CALL dbcsr_assert ("NOT", ANY (buffers%b%dirty(:)),&
         dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Should not release dirty buffers.",&
         __LINE__, error=error)
    buffers%b%refcount = buffers%b%refcount - 1
    n = 1
    !$ n = omp_get_max_threads()
    IF (buffers%b%refcount .LE. 0) THEN
       DO i = 1, n
          CALL dbcsr_data_release (buffers%b%buffers(i))
       ENDDO
       CALL dbcsr_data_release (buffers%b%main)
       DEALLOCATE (buffers%b%dirty)
       DEALLOCATE (buffers%b%buffers)
       DEALLOCATE (buffers%b%rcb)
       NULLIFY (buffers%b%dirty)
       NULLIFY (buffers%b%buffers)
       NULLIFY (buffers%b%rcb)
       IF (ASSOCIATED (buffers%b%backing)) THEN
          DO i = 1, n
             CALL dbcsr_data_clear_pointer (buffers%b%backing(i))
             CALL dbcsr_data_release (buffers%b%backing(i))
          ENDDO
          DEALLOCATE (buffers%b%backing)
          NULLIFY (buffers%b%backing)
       ENDIF
       DEALLOCATE (buffers%b)
       NULLIFY (buffers%b)
    ENDIF
    !$OMP END MASTER
    !$OMP BARRIER
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_buffers_release


! *****************************************************************************
!> \brief Checks whether any buffer is dirty
!> \param[in] buffers         Buffers to check
!> \retval dirty              Whether any buffer is dirty
! *****************************************************************************
  FUNCTION dbcsr_buffers_are_dirty (buffers, error) RESULT (dirty)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    LOGICAL                                  :: dirty

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_are_dirty', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    IF (careful_mod) &
         CALL dbcsr_assert (ASSOCIATED (buffers%b), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Buffers not valid", __LINE__,&
         error=error)
    dirty = ANY (buffers%b%dirty(:))
    CALL dbcsr_error_stop(error_handler, error)
  END FUNCTION dbcsr_buffers_are_dirty


  FUNCTION dbcsr_buffers_is_dirty (buffers, which_buffer, error) RESULT (dirty)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    INTEGER, INTENT(IN)                      :: which_buffer
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    LOGICAL                                  :: dirty

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_is_dirty', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

!   ---------------------------------------------------------------------------

    IF (careful_mod) &
         CALL dbcsr_error_set (routineN, error_handler, error)
    IF (debug_mod) &
         CALL dbcsr_assert (ASSOCIATED (buffers%b), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Buffers not valid", __LINE__,&
         error=error)
    dirty = buffers%b%dirty(which_buffer)
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handler, error)
  END FUNCTION dbcsr_buffers_is_dirty

! *****************************************************************************
!> \brief Sets buffer dirty flag
!> \param[in,out] buffers     The buffers to mark
!> \param[in] which_buffer    (optional) Which buffer to mark
!> \param[in] dirty           (optional) Whether to mark buffer as dirty
!>                            (TRUE, default) or clean (FALSE).
!> \param[in,out]             error
! *****************************************************************************
  SUBROUTINE dbcsr_buffers_mark_dirty (buffers, which_buffer, dirty, error)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    INTEGER, INTENT(IN), OPTIONAL            :: which_buffer
    LOGICAL, INTENT(IN), OPTIONAL            :: dirty
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_mark_dirty', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, which

!   ---------------------------------------------------------------------------

    IF (careful_mod) &
         CALL dbcsr_error_set (routineN, error_handler, error)
    IF (debug_mod) &
         CALL dbcsr_assert (ASSOCIATED (buffers%b), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Buffers not valid", __LINE__,&
         error=error)
    IF (PRESENT (which_buffer)) THEN
       which = which_buffer
    ELSE
       which = dbcsr_buffers_which_is_my (buffers)
    ENDIF
    IF (PRESENT (dirty)) THEN
       buffers%b%dirty(which) = dirty
    ELSE
       buffers%b%dirty(which) = .TRUE.
    ENDIF
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_buffers_mark_dirty

! *****************************************************************************
!> \brief Returns the number of buffers
!> \param[in] buffers
!> \result nbuffers
! *****************************************************************************
  FUNCTION dbcsr_buffers_get_num (buffers) RESULT (nbuffers)
    TYPE(dbcsr_block_buffer_obj), INTENT(IN) :: buffers
    INTEGER                                  :: nbuffers

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_get_num', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (debug_mod) &
         CALL dbcsr_assert (ASSOCIATED (buffers%b), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Buffers not valid", __LINE__,&
         error=error)
    nbuffers = SIZE (buffers%b%buffers)
  END FUNCTION dbcsr_buffers_get_num


! *****************************************************************************
!> \brief Returns which buffer belongs to the calling thread
!> \param[in] buffers
!> \result which
! *****************************************************************************
  FUNCTION dbcsr_buffers_which_is_my (buffers) RESULT (which)
    TYPE(dbcsr_block_buffer_obj), INTENT(IN) :: buffers
    INTEGER                                  :: which

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_which_is_my', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (debug_mod) &
         CALL dbcsr_assert (ASSOCIATED (buffers%b), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Buffers not valid", __LINE__,&
         error=error)
!$  which = omp_get_thread_num() + 1
!$  RETURN
    which = 1
  END FUNCTION dbcsr_buffers_which_is_my


! *****************************************************************************
!> \brief Returns the data area of a certain buffer
!> \param[in] buffers          buffers
!> \param[in] which_buffer     return area for this buffer
!> \param[in,out] error        error
!> \result area                Area belonging to buffer
! *****************************************************************************
  FUNCTION dbcsr_buffers_get_area (buffers, which_buffer, error) RESULT (area)
    TYPE(dbcsr_block_buffer_obj), INTENT(IN) :: buffers
    INTEGER, INTENT(IN)                      :: which_buffer
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    TYPE(dbcsr_data_obj)                     :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_get_area', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    CALL dbcsr_assert (ASSOCIATED (buffers%b), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Buffers not valid", __LINE__,&
         error=error)
    IF (debug_mod) THEN
       CALL dbcsr_assert (which_buffer, "GE", LBOUND (buffers%b%buffers,1),&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Invalid buffer number", __LINE__, error=error)
       CALL dbcsr_assert (which_buffer, "LE", UBOUND (buffers%b%buffers,1),&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Invalid buffer number", __LINE__, error=error)
    ENDIF
    area = buffers%b%buffers(which_buffer)
    CALL dbcsr_error_stop(error_handler, error)
  END FUNCTION dbcsr_buffers_get_area

! *****************************************************************************
!> \brief Resizes the buffer to the correct size.
!> \param[in,out] buffers     Buffers to resize
!> \param[in] row_size        (optional) Row size
!> \param[in] col_size        (optional) Column size
!> \param[in] minimum         (optional) If set, the sizes are set to at least
!>                            to those specified
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE dbcsr_buffers_resize (buffers, row_size, col_size, minimum, error)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    INTEGER, INTENT(IN), OPTIONAL            :: row_size, col_size
    LOGICAL, INTENT(IN), OPTIONAL            :: minimum
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_resize', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler, ibuff, new_cs, &
                                                new_rs, old_cs, old_rs
    LOGICAL                                  :: exact, existing_ok, valid
    TYPE(dbcsr_data_obj)                     :: area

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    IF (careful_mod) &
         CALL dbcsr_assert (ASSOCIATED (buffers%b), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Buffers not valid", __LINE__,&
         error=error)
    ibuff = dbcsr_buffers_which_is_my (buffers)
    CALL dbcsr_assert ("NOT", dbcsr_buffers_is_dirty(buffers, ibuff, error=error),&
         dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Can not resize dirty buffer!",&
         __LINE__,&
         error=error)
    IF (PRESENT (minimum)) THEN
       exact = .NOT. minimum
    ELSE
       exact = .TRUE.
    ENDIF
    area = dbcsr_buffers_get_area(buffers, ibuff, error=error)
    ! Get existing size.
    CALL dbcsr_data_get_sizes (area, old_rs, old_cs, valid, error=error)
    existing_ok = valid
    new_rs = old_rs
    new_cs = old_cs
    IF (PRESENT (row_size)) THEN
       new_rs = row_size
       IF (exact .AND. new_rs .NE. old_rs) existing_ok = .FALSE.
       IF (.NOT. exact .AND. new_rs .GT. old_rs) existing_ok = .FALSE.
    ENDIF
    IF (PRESENT (col_size)) THEN
       new_cs = col_size
       IF (exact .AND. new_cs .NE. old_cs) existing_ok = .FALSE.
       IF (.NOT. exact .AND. new_cs .GT. old_cs) existing_ok = .FALSE.
    ENDIF
    IF (.NOT. existing_ok) THEN
       CALL dbcsr_data_resize (area, (/ new_rs, new_cs /), error=error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_buffers_resize


! *****************************************************************************
!> \brief Flushes buffers
!> \param[in,out] buffers     Buffers to risez
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE dbcsr_buffers_flush (buffers, error)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_flush', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, ibuff, nthreads
    LOGICAL                                  :: reap, valid

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    IF (careful_mod) THEN
       CALL dbcsr_assert (ASSOCIATED (buffers%b), dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "Buffers not valid", __LINE__,&
            error=error)
    ENDIF
    reap = .FALSE.
    nthreads = 1
!$  nthreads = omp_get_num_threads()
    reap = nthreads .LT. dbcsr_buffers_get_num(buffers)
!$  reap = reap .OR. .NOT. omp_in_parallel()
    IF (reap) THEN
!$omp master
       DO ibuff = 1, dbcsr_buffers_get_num(buffers)
          IF (buffers%b%dirty(ibuff)) THEN
             CALL flush_single_buffer (buffers, ibuff, error=error)
          ENDIF
       ENDDO
!$omp end master
    ELSE
       ibuff = dbcsr_buffers_which_is_my (buffers)
       IF (dbcsr_buffers_is_dirty (buffers, ibuff, error=error)) THEN
          CALL flush_single_buffer (buffers, ibuff, error=error)
       ENDIF
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_buffers_flush


! *****************************************************************************
!> \brief Flushes one buffers
!> \param[in,out] buffers     Buffers
!> \param[in] which_buffer    Which buffer to flush
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE flush_single_buffer (buffers, which_buffer, error)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    INTEGER, INTENT(IN)                      :: which_buffer
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'flush_single_buffer', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: col, cs, dst_cs, dst_rs, &
                                                error_handler, nze, offset, &
                                                row, rs
    INTEGER, DIMENSION(2)                    :: bs
    LOGICAL                                  :: buff_tr, main_tr, valid
    TYPE(dbcsr_data_obj)                     :: area

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    IF (debug_mod) THEN ! This is already checked
       CALL dbcsr_assert (ASSOCIATED (buffers%b), dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "Buffers not valid", __LINE__,&
            error=error)
    ENDIF
    IF (careful_mod) THEN
       CALL dbcsr_assert (buffers%b%dirty(which_buffer), dbcsr_warning_level,&
            dbcsr_caller_error, routineN, "Should not flush clean buffer.",&
            __LINE__, error=error)
    ENDIF
    IF (buffers%b%dirty(which_buffer)) THEN
       area    = buffers%b%buffers(which_buffer)
       row     = buffers%b%rcb(1, which_buffer)
       col     = buffers%b%rcb(2, which_buffer)
       offset  = buffers%b%rcb(3, which_buffer)
       rs      = buffers%b%rcb(4, which_buffer)
       cs      = buffers%b%rcb(5, which_buffer)
       main_tr = BTEST (buffers%b%rcb(6, which_buffer), 1)
       buff_tr = BTEST (buffers%b%rcb(6, which_buffer), 0)
       nze = rs * cs
       IF (offset .GE. 0) THEN
          dst_rs = rs
          dst_cs = cs
       ELSE
          dst_rs = cs
          dst_cs = rs
       ENDIF
       CALL dbcsr_data_get_sizes (area, bs, valid, error=error)
       CALL dbcsr_block_partial_copy (dst=buffers%b%main,&
            dst_offset=ABS(offset)-1,&
            dst_rs = rs, dst_cs = cs, dst_tr = main_tr,&
            src = area, src_rs = rs, src_cs = cs, src_tr = buff_tr,&
            dst_r_lb=1, dst_c_lb=1, src_r_lb=1, src_c_lb=1,&
            nrow = rs, ncol = cs)
       buffers%b%dirty(which_buffer) = .FALSE.
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE flush_single_buffer


  SUBROUTINE dbcsr_buffers_set_p_2d_any (pointer_2d, row, col,&
       rsize, csize, main_tr, base_offset, buffers, buff_tr, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: pointer_2d
    INTEGER, INTENT(IN)                      :: row, col, rsize, csize
    LOGICAL, INTENT(IN)                      :: main_tr
    INTEGER, INTENT(IN)                      :: base_offset
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    LOGICAL, INTENT(IN)                      :: buff_tr
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_buffers_set_p_2d_any', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, ibuff
    TYPE(dbcsr_data_obj)                     :: buff_area

!    TYPE(dbcsr_data_obj), INTENT(IN)         :: base_1d_pointer
!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    ibuff = dbcsr_buffers_which_is_my (buffers)
    CALL dbcsr_buffers_flush (buffers, error=error)
    IF (dbcsr_buffers_is_dirty (buffers, ibuff, error=error)) THEN
       CALL flush_single_buffer (buffers, ibuff, error=error)
    ENDIF
    buff_area = dbcsr_buffers_get_area (buffers, ibuff, error=error)
    !
    CALL dbcsr_assert (dbcsr_data_get_type (pointer_2d), "EQ",&
         dbcsr_data_get_type (buff_area),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Incompatible data_types.",__LINE__,error)
    !CALL dbcsr_assert (dbcsr_data_get_type (base_1d_pointer), "EQ",&
    !     dbcsr_type_2d_to_1d (dbcsr_data_get_type (buff_area)),&
    !     dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
    !     "Incompatible data_types.",__LINE__,error)
    CALL dbcsr_assert (dbcsr_type_is_2d (dbcsr_data_get_type (pointer_2d)),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Pointer must be 2-D.",__LINE__,error)
    !CALL dbcsr_assert ("NOT", dbcsr_type_is_2d (dbcsr_data_get_type(base_1d_pointer)),&
    !     dbcsr_fatal_level, dbcsr_internal_error, routineN,&
    !     "Base pointer must be 1-D.",__LINE__,error)
    !
    SELECT CASE (dbcsr_data_get_type (pointer_2d))
    CASE (dbcsr_type_real_4_2d)
       CALL dbcsr_buffers_set_p_2d_s (pointer_2d%d%r2_sp, row, col,&
            rsize, csize, main_tr, base_offset,&
            buffers, buff_tr, error=error)
    CASE (dbcsr_type_real_8_2d)
       CALL dbcsr_buffers_set_p_2d_d (pointer_2d%d%r2_dp, row, col,&
            rsize, csize, main_tr, base_offset,&
            buffers, buff_tr, error=error)
    CASE (dbcsr_type_complex_4_2d)
       CALL dbcsr_buffers_set_p_2d_c (pointer_2d%d%c2_sp, row, col,&
            rsize, csize, main_tr, base_offset,&
            buffers, buff_tr, error=error)
    CASE (dbcsr_type_complex_8_2d)
       CALL dbcsr_buffers_set_p_2d_z (pointer_2d%d%c2_dp, row, col,&
            rsize, csize, main_tr, base_offset,&
            buffers, buff_tr, error=error)
    CASE default
       CALL dbcsr_assert (.FALSE.,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid pointer_2d.",__LINE__,error)
    END SELECT
    CALL dbcsr_error_stop (error_handler, error)
  END SUBROUTINE dbcsr_buffers_set_p_2d_any


#include "dbcsr_block_buffers_d.F"
#include "dbcsr_block_buffers_z.F"
#include "dbcsr_block_buffers_s.F"
#include "dbcsr_block_buffers_c.F"

END MODULE dbcsr_block_buffers
