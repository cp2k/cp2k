!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR data methods
!> \author  Urban Borstnik
!> \date    2010-06-15
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - 2010-02-18 Moved from dbcsr_methods
! *****************************************************************************
MODULE dbcsr_data_methods

  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_failure_level, dbcsr_fatal_level, &
       dbcsr_internal_error, dbcsr_postcondition_failed, &
       dbcsr_unimplemented_error_nr, dbcsr_warning_level, &
       dbcsr_wrong_args_error
  USE dbcsr_kinds,                     ONLY: real_4,&
                                             real_8
  USE dbcsr_message_passing,           ONLY: mp_allocate,&
                                             mp_deallocate
  USE dbcsr_ptr_util,                  ONLY: dbcsr_ptr_remapping,&
                                             ensure_array_size,&
                                             pointer_c_rank_remap2,&
                                             pointer_d_rank_remap2,&
                                             pointer_s_rank_remap2,&
                                             pointer_z_rank_remap2
  USE dbcsr_types,                     ONLY: &
       dbcsr_data_area_type, dbcsr_data_obj, dbcsr_obj, dbcsr_scalar_type, &
       dbcsr_type_complex_4, dbcsr_type_complex_4_2d, dbcsr_type_complex_8, &
       dbcsr_type_complex_8_2d, dbcsr_type_real_4, dbcsr_type_real_4_2d, &
       dbcsr_type_real_8, dbcsr_type_real_8_2d

  !$ USE OMP_LIB

  IMPLICIT NONE


  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_data_methods'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  PUBLIC :: dbcsr_type_is_2d, dbcsr_type_2d_to_1d, dbcsr_type_1d_to_2d
  PUBLIC :: dbcsr_scalar, dbcsr_scalar_one, dbcsr_scalar_zero,&
            dbcsr_scalar_are_equal, dbcsr_scalar_negative
  PUBLIC :: dbcsr_data_init, dbcsr_data_new, dbcsr_data_hold,&
            dbcsr_data_release, dbcsr_data_get_size, dbcsr_data_get_type,&
            dbcsr_data_reset_type, dbcsr_data_query_type
  PUBLIC :: dbcsr_data_resize
  PUBLIC :: dbcsr_get_data, &
            dbcsr_data_set_pointer,&
            dbcsr_data_clear_pointer, dbcsr_data_set_2d_pointer,&
            dbcsr_data_clear_2d_pointer, dbcsr_data_ensure_size,&
            dbcsr_data_get_sizes, dbcsr_data_verify_bounds,&
            dbcsr_data_exists, dbcsr_data_valid
  PUBLIC :: dbcsr_data_set_size_referenced, dbcsr_data_get_size_referenced
  PUBLIC :: dbcsr_get_data_p, dbcsr_get_data_p_s, dbcsr_get_data_p_c,&
            dbcsr_get_data_p_d, dbcsr_get_data_p_z,&
            dbcsr_get_data_p_2d_s, dbcsr_get_data_p_2d_d,&
            dbcsr_get_data_p_2d_c, dbcsr_get_data_p_2d_z
  !> \brief Encapsulates a scalar.
  INTERFACE dbcsr_scalar
     MODULE PROCEDURE dbcsr_scalar_s, dbcsr_scalar_d,&
                      dbcsr_scalar_c, dbcsr_scalar_z
  END INTERFACE

  INTERFACE dbcsr_data_set_pointer
     MODULE PROCEDURE set_data_p_s, set_data_p_d, set_data_p_c, set_data_p_z
     MODULE PROCEDURE set_data_p_2d_s, set_data_p_2d_d,&
                      set_data_p_2d_c, set_data_p_2d_z
     MODULE PROCEDURE set_data_area_area
  END INTERFACE

  INTERFACE dbcsr_get_data
     MODULE PROCEDURE get_data_s, get_data_d, get_data_c, get_data_z
     MODULE PROCEDURE get_data_m_s, get_data_m_d, get_data_m_c, get_data_m_z
     MODULE PROCEDURE get_data_2d_s, get_data_2d_d, get_data_2d_c, get_data_2d_z
  END INTERFACE


  INTERFACE dbcsr_get_data_p
     MODULE PROCEDURE dbcsr_get_data_c_s, dbcsr_get_data_c_c,&
                      dbcsr_get_data_c_d, dbcsr_get_data_c_z
  END INTERFACE

  INTERFACE dbcsr_data_get_sizes
     MODULE PROCEDURE dbcsr_data_get_sizes_any
     MODULE PROCEDURE dbcsr_data_get_sizes_1, dbcsr_data_get_sizes_2
  END INTERFACE

  INTERFACE dbcsr_data_query_type
     MODULE PROCEDURE query_type_s_1d, query_type_d_1d,&
                      query_type_c_1d, query_type_z_1d
     MODULE PROCEDURE query_type_s_2d, query_type_d_2d,&
                      query_type_c_2d, query_type_z_2d
  END INTERFACE

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
!> \brief Returns data type of a data area
!> \param[in] area         data area
!> \result data_type       data type of the data area
! *****************************************************************************
  PURE FUNCTION dbcsr_data_get_type (area) RESULT (data_type)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    CHARACTER                                :: data_type

    data_type = area%d%data_type
  END FUNCTION dbcsr_data_get_type

! Data type transformations
  FUNCTION data_type_2d_from_1d (type_1d) RESULT (type_2d)
    CHARACTER, INTENT(in)                    :: type_1d
    CHARACTER                                :: type_2d

    CHARACTER(len=*), PARAMETER :: routineN = 'data_type_2d_from_1d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

    SELECT CASE (type_1d)
    CASE (dbcsr_type_real_4)
       type_2d = dbcsr_type_real_4_2d
    CASE (dbcsr_type_real_8)
       type_2d = dbcsr_type_real_8_2d
    CASE (dbcsr_type_complex_4)
       type_2d = dbcsr_type_complex_4_2d
    CASE (dbcsr_type_complex_8)
       type_2d = dbcsr_type_complex_8_2d
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
            routineN, "Invalid data type.",__LINE__,error)
    END SELECT
  END FUNCTION data_type_2d_from_1d

! Data type transformations
  FUNCTION data_type_1d_from_2d (type_2d) RESULT (type_1d)
    CHARACTER, INTENT(IN)                    :: type_2d
    CHARACTER                                :: type_1d

    CHARACTER(len=*), PARAMETER :: routineN = 'data_type_1d_from_2d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

    SELECT CASE (type_2d)
    CASE (dbcsr_type_real_4_2d)
       type_1d = dbcsr_type_real_4
    CASE (dbcsr_type_real_8_2d)
       type_1d = dbcsr_type_real_8
    CASE (dbcsr_type_complex_4_2d)
       type_1d = dbcsr_type_complex_4
    CASE (dbcsr_type_complex_8_2d)
       type_1d = dbcsr_type_complex_8
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
            routineN, "Invalid data type.",__LINE__,error)
    END SELECT
  END FUNCTION data_type_1d_from_2d


! *****************************************************************************
!> \brief Initializes a data area
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_init (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    NULLIFY (area%d)
  END SUBROUTINE dbcsr_data_init

! *****************************************************************************
!> \brief Allocates pointers in the data type
!> \param[in,out] area        internal structure holding array pointers
!> \param[in] data_type       selects which array to allocate
!> \param[in] sizes           sizes to allocate to
!> \param[in] special_memory  whether to use special memory
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE internal_data_allocate (area, data_type, sizes,&
       special_memory, error)
    TYPE(dbcsr_data_area_type), &
      INTENT(INOUT)                          :: area
    CHARACTER, INTENT(IN)                    :: data_type
    INTEGER, DIMENSION(:), INTENT(IN)        :: sizes
    LOGICAL, INTENT(IN)                      :: special_memory
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_data_allocate', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: d, error_handler, stat

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    IF (dbg) &
         WRITE(*,*)routineN//" Setting to sizes", sizes
    IF (dbcsr_type_is_2d (data_type)) THEN
       CALL dbcsr_assert (SIZE(sizes), "EQ", 2, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Sizes must have 2 elements for 2-D data", __LINE__, error=error)
    ELSE
       CALL dbcsr_assert (SIZE(sizes), "EQ", 1, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Sizes must have 1 elements for 1-D data", __LINE__, error=error)
    ENDIF
    !
    stat = 0
    IF (special_memory) THEN
       SELECT CASE (data_type)
       CASE (dbcsr_type_real_4)
          CALL mp_allocate (area%r_sp, sizes(1), stat=stat)
       CASE (dbcsr_type_real_8)
          CALL mp_allocate (area%r_dp, sizes(1), stat=stat)
       CASE (dbcsr_type_complex_4)
          CALL mp_allocate (area%c_sp, sizes(1), stat=stat)
       CASE (dbcsr_type_complex_8)
          CALL mp_allocate (area%c_dp, sizes(1), stat=stat)
       CASE (dbcsr_type_real_8_2d, dbcsr_type_real_4_2d,&
            dbcsr_type_complex_8_2d, dbcsr_type_complex_4_2d)
          CALL dbcsr_assert ("NOT", special_memory,&
               dbcsr_fatal_level, dbcsr_caller_error, routineN,&
               "Can not use special memory with 2D data areas.",&
               __LINE__, error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type.",__LINE__,error)
       END SELECT
    ELSE
       SELECT CASE (data_type)
       CASE (dbcsr_type_real_8)
          ALLOCATE (area%r_dp(sizes(1)), stat=stat)
       CASE (dbcsr_type_real_4)
          ALLOCATE (area%r_sp(sizes(1)), stat=stat)
       CASE (dbcsr_type_complex_8)
          ALLOCATE (area%c_dp(sizes(1)), stat=stat)
       CASE (dbcsr_type_complex_4)
          ALLOCATE (area%c_sp(sizes(1)), stat=stat)
       CASE (dbcsr_type_real_8_2d)
          ALLOCATE (area%r2_dp(sizes(1), sizes(2)), stat=stat)
       CASE (dbcsr_type_real_4_2d)
          ALLOCATE (area%r2_sp(sizes(1), sizes(2)), stat=stat)
       CASE (dbcsr_type_complex_8_2d)
          ALLOCATE (area%c2_dp(sizes(1), sizes(2)), stat=stat)
       CASE (dbcsr_type_complex_4_2d)
          ALLOCATE (area%c2_sp(sizes(1), sizes(2)), stat=stat)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type.", __LINE__, error=error)
       END SELECT
    ENDIF
    CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level,&
         dbcsr_postcondition_failed, routineN,&
         "Error allocating memory", __LINE__, error=error)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE internal_data_allocate

! *****************************************************************************
!> \brief Allocates pointers in the data type
!> \param[in,out] area        internal structure holding array pointers
!> \param[in] data_type       selects which array to allocate
!> \param[in] sizes           sizes to allocate to
!> \param[in] special_memory  whether to use special memory
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE internal_data_deallocate (area, data_type,&
       special_memory, error)
    TYPE(dbcsr_data_area_type), &
      INTENT(INOUT)                          :: area
    CHARACTER, INTENT(IN)                    :: data_type
    LOGICAL, INTENT(IN)                      :: special_memory
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_data_deallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, stat

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    stat = 0
    IF (special_memory) THEN
       SELECT CASE (data_type)
       CASE (dbcsr_type_real_4)
          CALL mp_deallocate (area%r_sp, stat=stat)
          NULLIFY (area%r_sp)
       CASE (dbcsr_type_real_8)
          CALL mp_deallocate (area%r_dp, stat=stat)
          NULLIFY (area%r_dp)
       CASE (dbcsr_type_complex_4)
          CALL mp_deallocate (area%c_sp, stat=stat)
          NULLIFY (area%c_sp)
       CASE (dbcsr_type_complex_8)
          CALL mp_deallocate (area%c_dp, stat=stat)
          NULLIFY (area%c_dp)
       CASE (dbcsr_type_real_8_2d, dbcsr_type_real_4_2d,&
            dbcsr_type_complex_8_2d, dbcsr_type_complex_4_2d)
          CALL dbcsr_assert ("NOT", special_memory,&
               dbcsr_fatal_level, dbcsr_caller_error, routineN,&
               "Can not use special memory with 2D data areas.",&
               __LINE__, error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type.",__LINE__,error)
       END SELECT
    ELSE
       SELECT CASE (data_type)
       CASE (dbcsr_type_real_8)
          DEALLOCATE (area%r_dp, stat=stat)
          NULLIFY (area%r_dp)
       CASE (dbcsr_type_real_4)
          DEALLOCATE (area%r_sp, stat=stat)
          NULLIFY (area%r_sp)
       CASE (dbcsr_type_complex_8)
          DEALLOCATE (area%c_dp, stat=stat)
          NULLIFY (area%c_dp)
       CASE (dbcsr_type_complex_4)
          DEALLOCATE (area%c_sp, stat=stat)
          NULLIFY (area%c_sp)
       CASE (dbcsr_type_real_8_2d)
          DEALLOCATE (area%r2_dp, stat=stat)
          NULLIFY (area%r2_dp)
       CASE (dbcsr_type_real_4_2d)
          DEALLOCATE (area%r2_sp, stat=stat)
          NULLIFY (area%r2_sp)
       CASE (dbcsr_type_complex_8_2d)
          DEALLOCATE (area%c2_dp, stat=stat)
          NULLIFY (area%c2_dp)
       CASE (dbcsr_type_complex_4_2d)
          DEALLOCATE (area%c2_sp, stat=stat)
          NULLIFY (area%c2_sp)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type.", __LINE__, error=error)
       END SELECT
    ENDIF
    CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level,&
         dbcsr_postcondition_failed, routineN,&
         "Error deallocating memory", __LINE__, error=error)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE internal_data_deallocate


! *****************************************************************************
!> \brief Initializes a data area and all the actual data pointers
!> \param[inout] area         data area
!> \param[in] data_type       select data type to use
!> \param[in] data_size       (optional) allocate this much data
!> \param[in] data_size2      (optional) second dimension data size
!> \param[in] special_memory  (optional) whether to use specially-allocated
!>                            memory, default is False
! *****************************************************************************
  SUBROUTINE dbcsr_data_new (area, data_type, data_size, data_size2,&
       special_memory)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    CHARACTER, INTENT(IN)                    :: data_type
    INTEGER, INTENT(IN), OPTIONAL            :: data_size, data_size2
    LOGICAL, INTENT(IN), OPTIONAL            :: special_memory

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: d, total_size
    INTEGER, DIMENSION(2)                    :: sizes
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (.NOT. ASSOCIATED (area%d)) THEN
       ALLOCATE (area%d)
    ENDIF
    NULLIFY (area%d%r_sp)
    NULLIFY (area%d%r_dp)
    NULLIFY (area%d%c_sp)
    NULLIFY (area%d%c_dp)
    NULLIFY (area%d%r2_sp)
    NULLIFY (area%d%r2_dp)
    NULLIFY (area%d%c2_sp)
    NULLIFY (area%d%c2_dp)
    area%d%refcount = 1
    IF (PRESENT (special_memory)) THEN
       area%d%special = special_memory
    ELSE
       area%d%special = .FALSE.
    ENDIF
    IF (PRESENT (data_size)) THEN
       IF (dbcsr_type_is_2d (data_type)) THEN
          CALL dbcsr_assert (PRESENT (data_size2), dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "Must specify 2 sizes for 2-D data", __LINE__, error=error)
          d = 2
          sizes(1) = data_size
          sizes(2) = data_size2
          total_size = data_size * data_size2
       ELSE
          d = 1
          sizes(1) = data_size
          total_size = data_size
       ENDIF
       CALL internal_data_allocate (area%d, data_type, sizes(1:d),&
            special_memory = area%d%special, error=error)
       CALL dbcsr_data_set_size_referenced (area, total_size)
    ELSE
       CALL dbcsr_data_set_size_referenced (area, 0)
    ENDIF
    area%d%data_type = data_type
  END SUBROUTINE dbcsr_data_new

! *****************************************************************************
!> \brief Removes a reference and/or clears the data area.
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_release (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_release', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (.NOT. ASSOCIATED (area%d)) THEN
       RETURN
    ENDIF
    CALL dbcsr_assert (area%d%refcount .GT. 0, dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "Data seems to be unreferenced.",__LINE__,error)
    IF (area%d%refcount .GT. 0) THEN
       !$OMP ATOMIC
       area%d%refcount = area%d%refcount - 1
    ELSE
       RETURN
    ENDIF
    ! If we're releasing the last reference, then free the memory.
    IF (area%d%refcount .EQ. 0) THEN
       IF (area%d%special) THEN
          SELECT CASE (area%d%data_type)
          CASE (dbcsr_type_real_4)
             CALL mp_deallocate (area%d%r_sp)
             NULLIFY (area%d%r_sp)
          CASE (dbcsr_type_real_8)
             CALL mp_deallocate (area%d%r_dp)
             NULLIFY (area%d%r_dp)
          CASE (dbcsr_type_complex_4)
             CALL mp_deallocate (area%d%c_sp)
             NULLIFY (area%d%c_sp)
          CASE (dbcsr_type_complex_8)
             CALL mp_deallocate (area%d%c_dp)
             NULLIFY (area%d%c_dp)
          CASE (dbcsr_type_real_8_2d, dbcsr_type_real_4_2d,&
               dbcsr_type_complex_8_2d, dbcsr_type_complex_4_2d)
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_internal_error,&
                  routineN, "How could 2D data areas use special memory?",__LINE__,error)
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
                  routineN, "Invalid data type.",__LINE__,error)
          END SELECT
       ELSE
          SELECT CASE (area%d%data_type)
          CASE (dbcsr_type_real_8)
             IF (ASSOCIATED (area%d%r_dp)) &
                  DEALLOCATE (area%d%r_dp)
             NULLIFY (area%d%r_dp)
          CASE (dbcsr_type_real_4)
             IF (ASSOCIATED (area%d%r_sp)) &
                  DEALLOCATE (area%d%r_sp)
             NULLIFY (area%d%r_sp)
          CASE (dbcsr_type_complex_8)
             IF (ASSOCIATED (area%d%c_dp)) &
                  DEALLOCATE (area%d%c_dp)
             NULLIFY (area%d%c_dp)
          CASE (dbcsr_type_complex_4)
             IF (ASSOCIATED (area%d%c_sp)) &
                  DEALLOCATE (area%d%c_sp)
             NULLIFY (area%d%c_sp)
          CASE (dbcsr_type_real_8_2d)
             IF (ASSOCIATED (area%d%r2_dp)) &
                  DEALLOCATE (area%d%r2_dp)
             NULLIFY (area%d%r2_dp)
          CASE (dbcsr_type_real_4_2d)
             IF (ASSOCIATED (area%d%r2_sp)) &
                  DEALLOCATE (area%d%r2_sp)
             NULLIFY (area%d%r2_sp)
          CASE (dbcsr_type_complex_8_2d)
             IF (ASSOCIATED (area%d%c2_dp)) &
                  DEALLOCATE (area%d%c2_dp)
             NULLIFY (area%d%c2_dp)
          CASE (dbcsr_type_complex_4_2d)
             IF (ASSOCIATED (area%d%c2_sp)) &
                  DEALLOCATE (area%d%c2_sp)
             NULLIFY (area%d%c2_sp)
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
                  routineN, "Invalid data type.",__LINE__,error)
          END SELECT
       ENDIF
       DEALLOCATE (area%d)
       NULLIFY (area%d)
    ENDIF
  END SUBROUTINE dbcsr_data_release

! *****************************************************************************
!> \brief Clears pointers from the data area.
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_clear_pointer (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_clear_pointer', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (.NOT. ASSOCIATED (area%d)) THEN
       RETURN
    ENDIF
    CALL dbcsr_assert (area%d%refcount .GT. 0, dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "Data seems to be unreferenced.",__LINE__,error)
    SELECT CASE (area%d%data_type)
    CASE (dbcsr_type_real_4)
       NULLIFY (area%d%r_sp)
    CASE (dbcsr_type_real_8)
       NULLIFY (area%d%r_dp)
    CASE (dbcsr_type_complex_4)
       NULLIFY (area%d%c_sp)
    CASE (dbcsr_type_complex_8)
       NULLIFY (area%d%c_dp)
    CASE (dbcsr_type_real_8_2d)
       NULLIFY (area%d%r2_dp)
    CASE (dbcsr_type_real_4_2d)
       NULLIFY (area%d%r2_sp)
    CASE (dbcsr_type_complex_8_2d)
       NULLIFY (area%d%c2_dp)
    CASE (dbcsr_type_complex_4_2d)
       NULLIFY (area%d%c2_sp)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
            routineN, "Invalid data type.",__LINE__,error)
    END SELECT
  END SUBROUTINE dbcsr_data_clear_pointer

! *****************************************************************************
!> \brief Clears pointers from the data area.
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_reset_type (area, new_type)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    CHARACTER, INTENT(IN)                    :: new_type

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_reset_type', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (.NOT. ASSOCIATED (area%d)) THEN
       CALL dbcsr_assert (ASSOCIATED (area%d), dbcsr_fatal_level,&
            dbcsr_caller_error,&
            routineN, "Data not initialized.",__LINE__,error)
    ENDIF
    CALL dbcsr_data_clear_pointer (area)
    area%d%data_type = new_type
  END SUBROUTINE dbcsr_data_reset_type

! *****************************************************************************
!> \brief Checks whether a data area is valid
!> \param[in] area         data area
!> \result valid           whether the data area is valid
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_data_valid (area) RESULT (valid)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    LOGICAL                                  :: valid

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_valid', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    valid = ASSOCIATED (area%d)
  END FUNCTION dbcsr_data_valid

! *****************************************************************************
!> \brief Checks whether a data pointer exists
!> \param[in] area         data area
!> \param[in,out] error    error
!> \result valid           whether the data pointer exists
! *****************************************************************************
  FUNCTION dbcsr_data_exists (area, error) RESULT (valid)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    LOGICAL                                  :: valid

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_exists', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    valid = dbcsr_data_valid (area)
    CALL dbcsr_assert (valid, dbcsr_fatal_level, dbcsr_wrong_args_error,&
         routineN, "Data area is invalid.", __LINE__, error)
    SELECT CASE (area%d%data_type)
    CASE (dbcsr_type_real_8)
       valid = ASSOCIATED (area%d%r_dp)
    CASE (dbcsr_type_real_4)
       valid = ASSOCIATED (area%d%r_sp)
    CASE (dbcsr_type_complex_8)
       valid = ASSOCIATED (area%d%c_dp)
    CASE (dbcsr_type_complex_4)
       valid = ASSOCIATED (area%d%c_sp)
    CASE (dbcsr_type_real_8_2d)
       valid = ASSOCIATED (area%d%r2_dp)
    CASE (dbcsr_type_real_4_2d)
       valid = ASSOCIATED (area%d%r2_sp)
    CASE (dbcsr_type_complex_8_2d)
       valid = ASSOCIATED (area%d%c2_dp)
    CASE (dbcsr_type_complex_4_2d)
       valid = ASSOCIATED (area%d%c2_sp)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Incorrect data type", __LINE__, error)
    END SELECT
  END FUNCTION dbcsr_data_exists

! *****************************************************************************
!> \brief Registers another use of the data area
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_hold (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_hold', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (area%d), dbcsr_caller_error, dbcsr_warning_level,&
         routineN, "Can't hold an empty data area.",__LINE__,error)
    CALL dbcsr_assert (area%d%refcount .GT. 0, dbcsr_caller_error, dbcsr_warning_level,&
         routineN, "Should not hold an area with zero references.",__LINE__,error)
    IF (.NOT. ASSOCIATED (area%d)) THEN
       RETURN
    ENDIF
    !$OMP ATOMIC
    area%d%refcount = area%d%refcount + 1
  END SUBROUTINE dbcsr_data_hold



! *****************************************************************************
!> \brief Points data area data pointers to another data area
!>
!> Assumes that no memory will be lost when repointing the pointer in the data
!> area and that the area is initialized.
!> \param[in,out] area         data area to repoint
!> \param[in] rsize, csize     size of data area to point to
!> \param[in] pointee          data area to point to
!> \param[in] source_lb        (optional) point to this offset in pointee
! *****************************************************************************
  SUBROUTINE set_data_area_area (area, rsize, csize, pointee, source_lb)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, INTENT(IN)                      :: rsize, csize
    TYPE(dbcsr_data_obj), INTENT(IN)         :: pointee
    INTEGER, INTENT(IN), OPTIONAL            :: source_lb

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_area_area', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: dt1, dt2
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: c_sp
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: c_dp
    INTEGER                                  :: bp, nze
    LOGICAL                                  :: compatible, pointee_is_2d, rmp
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r_sp
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r_dp
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    bp = 1 ; IF (PRESENT (source_lb)) bp = source_lb
    nze = rsize*csize
    dt1 = area%d%data_type
    dt2 = pointee%d%data_type
    compatible = dt1 .EQ. dt2 .OR. dt1 .EQ. dbcsr_type_1d_to_2d (dt2)
    CALL dbcsr_assert (compatible, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN,&
         "Can not point 1-d pointer to 2-d data",__LINE__,error)
    pointee_is_2d = dbcsr_type_is_2d (dt2)
    CALL dbcsr_assert (.NOT. PRESENT (source_lb) .OR. .NOT. pointee_is_2d, &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Lower bound specification not possible with 2-d data",__LINE__,error)
    ! Check if size is OK.
    CALL dbcsr_assert (bp, "GE", 1,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Attempt to point out of bounds",__LINE__,error)
    CALL dbcsr_assert (bp + nze - 1, "LE", dbcsr_data_get_size (pointee),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Attempt to point out of bounds",__LINE__,error)
    ! There's a remap if the ranks are compatible but not equal.
    rmp = dt1 .NE. dt2
    IF (.NOT. dbcsr_ptr_remapping) &
         CALL dbcsr_assert (.NOT. rmp, dbcsr_fatal_level, dbcsr_internal_error,&
         routineN, "Compiler does not support pointer rank remapping.",&
         __LINE__, error=error)
    SELECT CASE (dt2)
    CASE (dbcsr_type_real_4_2d)
       area%d%r2_sp => pointee%d%r2_sp(1:rsize,1:csize)
    CASE (dbcsr_type_real_4)
       IF (rmp) THEN
          r_sp => dbcsr_get_data_p_s (pointee, bp, bp+nze-1)
          CALL pointer_s_rank_remap2(area%d%r2_sp, rsize, csize,&
               r_sp)
       ELSE
          area%d%r_sp => dbcsr_get_data_p_s (pointee, bp, bp+nze-1)
       ENDIF
    CASE (dbcsr_type_real_8_2d)
       area%d%r2_dp => pointee%d%r2_dp(1:rsize,1:csize)
    CASE (dbcsr_type_real_8)
       IF (rmp) THEN
          r_dp => dbcsr_get_data_p_d (pointee, bp, bp+nze-1)
          CALL pointer_d_rank_remap2(area%d%r2_dp, rsize, csize,&
               r_dp)
       ELSE
          area%d%r_dp => dbcsr_get_data_p_d (pointee, bp, bp+nze-1)
       ENDIF
    CASE (dbcsr_type_complex_4_2d)
       area%d%c2_sp => pointee%d%c2_sp(1:rsize,1:csize)
    CASE (dbcsr_type_complex_4)
       IF (rmp) THEN
          c_sp => dbcsr_get_data_p_c (pointee, bp, bp+nze-1)
          CALL pointer_c_rank_remap2(area%d%c2_sp, rsize, csize,&
               c_sp)
       ELSE
          area%d%c_sp => dbcsr_get_data_p_c (pointee, bp, bp+nze-1)
       ENDIF
    CASE (dbcsr_type_complex_8_2d)
       area%d%c2_dp => pointee%d%c2_dp(1:rsize,1:csize)
    CASE (dbcsr_type_complex_8)
       IF (rmp) THEN
          c_dp => dbcsr_get_data_p_z (pointee, bp, bp+nze-1)
          CALL pointer_z_rank_remap2(area%d%c2_dp, rsize, csize,&
               c_dp)
       ELSE
          area%d%c_dp => dbcsr_get_data_p_z (pointee, bp, bp+nze-1)
       ENDIF
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
    CALL dbcsr_data_set_size_referenced (area, rsize*csize)
    CALL dbcsr_assert (dbcsr_data_get_size_referenced (area), "EQ",&
         dbcsr_data_get_size (area), dbcsr_fatal_level, dbcsr_internal_error,&
         routineN, "Size mismatch", __LINE__, error)
  END SUBROUTINE set_data_area_area


! *****************************************************************************
!> \brief Returns the allocated data size
!> \param[in] area       data area
!> \retval data_size      size of data
! *****************************************************************************
  FUNCTION dbcsr_data_get_size (area) RESULT (data_size)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER                                  :: data_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_get_size', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

    data_size = 0
    IF (ASSOCIATED (area%d)) THEN
       SELECT CASE (area%d%data_type)
          CASE (dbcsr_type_real_8)
             IF (ASSOCIATED (area%d%r_dp))&
                  data_size = SIZE (area%d%r_dp)
          CASE (dbcsr_type_real_4)
             IF (ASSOCIATED (area%d%r_sp))&
                  data_size = SIZE (area%d%r_sp)
          CASE (dbcsr_type_complex_8)
             IF (ASSOCIATED (area%d%c_dp))&
                  data_size = SIZE (area%d%c_dp)
          CASE (dbcsr_type_complex_4)
             IF (ASSOCIATED (area%d%c_sp))&
                  data_size = SIZE (area%d%c_sp)
          CASE (dbcsr_type_real_8_2d)
             IF (ASSOCIATED (area%d%r2_dp))&
                  data_size = SIZE (area%d%r2_dp)
          CASE (dbcsr_type_real_4_2d)
             IF (ASSOCIATED (area%d%r2_sp))&
                  data_size = SIZE (area%d%r2_sp)
          CASE (dbcsr_type_complex_8_2d)
             IF (ASSOCIATED (area%d%c2_dp))&
                  data_size = SIZE (area%d%c2_dp)
          CASE (dbcsr_type_complex_4_2d)
             IF (ASSOCIATED (area%d%c2_sp))&
                  data_size = SIZE (area%d%c2_sp)
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
                  routineN, "Incorrect data type",__LINE__,error)
          END SELECT
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_warning_level, dbcsr_caller_error, routineN,&
            "Uninitialized data area",__LINE__,error)
       data_size = 0
    ENDIF
  END FUNCTION dbcsr_data_get_size


! *****************************************************************************
!> \brief Ensures a minimum size of a previously-setup data area.
!>
!> The data area must have been previously setup with dbcsr_data_new.
!> \param[inout] area         data area
!> \param[in] data_size       allocate this much data
!> \param[in] nocopy          (optional) do not keep potentially existing data,
!>                            default is to keep it
!> \param[in] zero_pad        (optional) pad new data with zeros
!> \param[in] factor          (optional) increase size by this factor
! *****************************************************************************
  SUBROUTINE dbcsr_data_ensure_size (area, data_size, nocopy, zero_pad, factor, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, INTENT(IN)                      :: data_size
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, zero_pad
    REAL, INTENT(IN), OPTIONAL               :: factor
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_ensure_size', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: current_size, error_handler
    LOGICAL                                  :: nocp, pad

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert(ASSOCIATED (area%d), dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Data area must be setup.",__LINE__,error)
    current_size = dbcsr_data_get_size (area)
    CALL dbcsr_data_set_size_referenced (area, data_size)
    IF (current_size .GT. 1 .AND. current_size .GE. data_size) THEN
       CALL dbcsr_error_stop(error_handler, error)
       RETURN
    ENDIF
    !
    nocp = .FALSE.
    IF (PRESENT (nocopy)) nocp = nocopy
    pad = .FALSE.
    IF (PRESENT (zero_pad)) pad = zero_pad
    !
    SELECT CASE (area%d%data_type)
    CASE (dbcsr_type_real_8)
       IF (.NOT. ASSOCIATED (area%d%r_dp)) THEN
          IF (area%d%special) THEN
             CALL mp_allocate (area%d%r_dp, data_size)
          ELSE
             ALLOCATE (area%d%r_dp(data_size))
          ENDIF
          IF (pad) area%d%r_dp(:) = 0.0_real_8
       ELSE
          CALL ensure_array_size (area%d%r_dp, ub=data_size,&
               special=area%d%special, nocopy=nocp, zero_pad=zero_pad,&
               factor=factor,error=error)
       ENDIF
    CASE (dbcsr_type_real_4)
       IF (.NOT. ASSOCIATED (area%d%r_sp)) THEN
          IF (area%d%special) THEN
             CALL mp_allocate (area%d%r_sp, data_size)
          ELSE
             ALLOCATE (area%d%r_sp(data_size))
          ENDIF
          IF (pad) area%d%r_sp(:) = 0.0_real_4
       ELSE
          CALL ensure_array_size (area%d%r_sp, ub=data_size,&
               special=area%d%special, nocopy=nocp, zero_pad=zero_pad,&
               factor=factor,error=error)
       ENDIF
    CASE (dbcsr_type_complex_8)
       IF (.NOT. ASSOCIATED (area%d%c_dp)) THEN
          IF (area%d%special) THEN
             CALL mp_allocate (area%d%c_dp, data_size)
          ELSE
             ALLOCATE (area%d%c_dp(data_size))
          ENDIF
          IF (pad) area%d%c_dp(:) = 0.0_real_8
       ELSE
          CALL ensure_array_size (area%d%c_dp, ub=data_size,&
               special=area%d%special, nocopy=nocp, zero_pad=zero_pad,&
               factor=factor,error=error)
       ENDIF
    CASE (dbcsr_type_complex_4)
       IF (.NOT. ASSOCIATED (area%d%c_sp)) THEN
          IF (area%d%special) THEN
             CALL mp_allocate (area%d%c_sp, data_size)
          ELSE
             ALLOCATE (area%d%c_sp(data_size))
          ENDIF
          IF (pad) area%d%c_sp(:) = 0.0_real_4
       ELSE
          CALL ensure_array_size (area%d%c_sp, ub=data_size,&
               special=area%d%special, nocopy=nocp, zero_pad=zero_pad,&
               factor=factor,error=error)
       ENDIF
    CASE default
       CALL dbcsr_assert(.FALSE., dbcsr_failure_level,&
            dbcsr_unimplemented_error_nr, routineN,&
            "Invalid data type are supported",__LINE__,error)
    END SELECT
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_ensure_size


! *****************************************************************************
!> \brief Verifies bounds of a data area
!> \param[in] area             Data area
!> \param[in] lb               lower bounds
!> \param[in] ub               upper bounds
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE dbcsr_data_verify_bounds (area, lb, ub, error)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER, DIMENSION(:), INTENT(IN)        :: lb, ub
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_verify_bounds', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = .TRUE., dbg = .TRUE.

    CHARACTER                                :: data_type
    INTEGER                                  :: error_handler

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    data_type = dbcsr_data_get_type (area)
    IF (dbcsr_type_is_2d (data_type)) THEN
       CALL dbcsr_assert (SIZE (lb), "EQ", 2, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "size must be 2 for 2-d lb", __LINE__, error=error)
       CALL dbcsr_assert (SIZE (ub), "EQ", 2, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "size must be 2 for 2-d ub", __LINE__, error=error)
    ELSE
       CALL dbcsr_assert (SIZE (lb), "EQ", 1, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "size must be 1 for 1-d lb", __LINE__, error=error)
       CALL dbcsr_assert (SIZE (ub), "EQ", 1, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "size must be 1 for 1-d ub", __LINE__, error=error)
    ENDIF
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_assert(lb(1), "GE", LBOUND(area%d%r_sp,1),&
            dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "lb r_sp", __LINE__, error)
       CALL dbcsr_assert(ub(1), "LE", UBOUND(area%d%r_sp,1),&
            dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "ub r_sp",__LINE__,error)
    CASE (dbcsr_type_real_4_2d)
       CALL dbcsr_assert(lb(1), "GE", LBOUND(area%d%r2_sp,1),&
            dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "lb r_sp 2d", __LINE__, error)
       CALL dbcsr_assert(ub(1), "LE", UBOUND(area%d%r2_sp,1),&
            dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "ub r_sp 2d",__LINE__,error)
       CALL dbcsr_assert(lb(2), "GE", LBOUND(area%d%r2_sp,2),&
            dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "lb r_sp 2d", __LINE__, error)
       CALL dbcsr_assert(ub(2), "LE", UBOUND(area%d%r2_sp,2),&
            dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "ub r_sp 2d",__LINE__,error)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_assert(lb(1), "GE", LBOUND(area%d%r_dp,1),&
            dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "lb r_dp", __LINE__, error)
       CALL dbcsr_assert(ub(1), "LE", UBOUND(area%d%r_dp,1),&
            dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "ub r_dp",__LINE__,error)
    CASE (dbcsr_type_real_8_2d)
       CALL dbcsr_assert(lb(1), "GE", LBOUND(area%d%r2_dp,1),&
            dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "lb r_dp 2d", __LINE__, error)
       CALL dbcsr_assert(ub(1), "LE", UBOUND(area%d%r2_dp,1),&
            dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "ub r_dp 2d",__LINE__,error)
       CALL dbcsr_assert(lb(2), "GE", LBOUND(area%d%r2_dp,2),&
            dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "lb r_dp 2d", __LINE__, error)
       CALL dbcsr_assert(ub(2), "LE", UBOUND(area%d%r2_dp,2),&
            dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "ub r_dp 2d",__LINE__,error)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_assert(lb(1), "GE", LBOUND(area%d%c_sp,1),&
            dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "lb c_sp", __LINE__, error)
       CALL dbcsr_assert(ub(1), "LE", UBOUND(area%d%c_sp,1),&
            dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "ub c_sp",__LINE__,error)
    CASE (dbcsr_type_complex_4_2d)
       CALL dbcsr_assert(lb(1), "GE", LBOUND(area%d%c2_sp,1),&
            dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "lb c_sp 2d", __LINE__, error)
       CALL dbcsr_assert(ub(1), "LE", UBOUND(area%d%c2_sp,1),&
            dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "ub c_sp 2d",__LINE__,error)
       CALL dbcsr_assert(lb(2), "GE", LBOUND(area%d%c2_sp,2),&
            dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "lb c_sp 2d", __LINE__, error)
       CALL dbcsr_assert(ub(2), "LE", UBOUND(area%d%c2_sp,2),&
            dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "ub c_sp 2d",__LINE__,error)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_assert(lb(1), "GE", LBOUND(area%d%c_dp,1),&
            dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "lb c_dp", __LINE__, error)
       CALL dbcsr_assert(ub(1), "LE", UBOUND(area%d%c_dp,1),&
            dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "ub c_dp",__LINE__,error)
    CASE (dbcsr_type_complex_8_2d)
       CALL dbcsr_assert(lb(1), "GE", LBOUND(area%d%c2_dp,1),&
            dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "lb c_dp 2d", __LINE__, error)
       CALL dbcsr_assert(ub(1), "LE", UBOUND(area%d%c2_dp,1),&
            dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "ub c_dp 2d",__LINE__,error)
       CALL dbcsr_assert(lb(2), "GE", LBOUND(area%d%c2_dp,2),&
            dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "lb c_dp 2d", __LINE__, error)
       CALL dbcsr_assert(ub(2), "LE", UBOUND(area%d%c2_dp,2),&
            dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "ub c_dp 2d",__LINE__,error)
    CASE default
       CALL dbcsr_assert(.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type", __LINE__, error=error)
    END SELECT
    CALL dbcsr_error_stop (error_handler, error)
  END SUBROUTINE dbcsr_data_verify_bounds


! *****************************************************************************
!> \brief Points a 2-d pointer to current 1-d data
!>
!> \param[in,out] area         data area to work on
!> \param[in] rowsize          row size of 2-d data
!> \param[in] colsize          (optional) column size of 2-d data
!> \param[in] offset           (optional) offset
! *****************************************************************************
  SUBROUTINE dbcsr_data_set_2d_pointer (area, rowsize, colsize, offset)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, INTENT(IN)                      :: rowsize
    INTEGER, INTENT(IN), OPTIONAL            :: colsize, offset

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_set_2d_pointer', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: dt1
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: c_sp
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: c_dp
    INTEGER                                  :: bp, csize, nze
    LOGICAL                                  :: compatible
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r_sp
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r_dp
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    dt1 = area%d%data_type
    compatible = dt1 .EQ. dbcsr_type_real_4&
            .OR. dt1 .EQ. dbcsr_type_real_8&
            .OR. dt1 .EQ. dbcsr_type_complex_4&
            .OR. dt1 .EQ. dbcsr_type_complex_8
    CALL dbcsr_assert (compatible, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Must target 1-D data",__LINE__,error)
    IF (PRESENT (offset)) THEN
       bp = offset
       nze = dbcsr_data_get_size (area) - offset + 1
    ELSE
       bp = 1
       nze = dbcsr_data_get_size (area)
    ENDIF
    IF (PRESENT (colsize)) THEN
       csize = colsize
    ELSE
       csize = nze / rowsize
    ENDIF
    SELECT CASE (dt1)
    CASE (dbcsr_type_real_4)
       r_sp => dbcsr_get_data_p_s (area, bp, bp+nze-1)
       CALL pointer_s_rank_remap2(area%d%r2_sp, rowsize, csize,&
            r_sp)
    CASE (dbcsr_type_real_8)
       r_dp => dbcsr_get_data_p_d (area, bp, bp+nze-1)
       CALL pointer_d_rank_remap2(area%d%r2_dp, rowsize, csize,&
            r_dp)
    CASE (dbcsr_type_complex_4)
       c_sp => dbcsr_get_data_p_c (area, bp, bp+nze-1)
       CALL pointer_c_rank_remap2(area%d%c2_sp, rowsize, csize,&
            c_sp)
    CASE (dbcsr_type_complex_8)
       c_dp => dbcsr_get_data_p_z (area, bp, bp+nze-1)
       CALL pointer_z_rank_remap2(area%d%c2_dp, rowsize, csize,&
            c_dp)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
  END SUBROUTINE dbcsr_data_set_2d_pointer

! *****************************************************************************
!> \brief Nullifies a 2-d pointer to current 1-d data
!>
!> \param[in,out] area         data area to work on
! *****************************************************************************
  SUBROUTINE dbcsr_data_clear_2d_pointer (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_clear_2d_pointer', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: dt1
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    dt1 = area%d%data_type
    SELECT CASE (dt1)
    CASE (dbcsr_type_real_4)
       NULLIFY (area%d%r2_sp)
    CASE (dbcsr_type_real_8)
       NULLIFY (area%d%r2_dp)
    CASE (dbcsr_type_complex_4)
       NULLIFY (area%d%c2_sp)
    CASE (dbcsr_type_complex_8)
       NULLIFY (area%d%c2_dp)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
  END SUBROUTINE dbcsr_data_clear_2d_pointer
 

! *****************************************************************************
!> \brief Returns the allocated data size
!> \param[in] area       data area to query for size
!> \param[out] sizes     array with the data sizes
!> \param[out] valid     whether the data is actually allocated
!> \param[in,out] error  error
! *****************************************************************************
  SUBROUTINE dbcsr_data_get_sizes_any (area, sizes, valid, error)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER, DIMENSION(:), INTENT(OUT)       :: sizes
    LOGICAL, INTENT(OUT)                     :: valid
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_get_sizes_any', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: d, error_handler

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)

    valid = .FALSE.
    sizes(:) = 0
    IF (ASSOCIATED (area%d)) THEN
       IF (dbcsr_type_is_2d (area%d%data_type)) THEN
          CALL dbcsr_assert (SIZE(sizes), "EQ", 2, dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "Sizes must have 2 elements for 2-D data", __LINE__, error=error)
       ELSE
          CALL dbcsr_assert (SIZE(sizes), "EQ", 1, dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "Sizes must have 1 elements for 1-D data", __LINE__, error=error)
       ENDIF
       SELECT CASE (area%d%data_type)
       CASE (dbcsr_type_real_8)
          IF (ASSOCIATED (area%d%r_dp)) THEN
             sizes(1) = SIZE (area%d%r_dp)
             valid = .TRUE.
          ENDIF
       CASE (dbcsr_type_real_4)
          IF (ASSOCIATED (area%d%r_sp)) THEN
             sizes(1) = SIZE (area%d%r_sp)
             valid = .TRUE.
          ENDIF
       CASE (dbcsr_type_complex_8)
          IF (ASSOCIATED (area%d%c_dp)) THEN
             sizes(1) = SIZE (area%d%c_dp)
             valid = .TRUE.
          ENDIF
       CASE (dbcsr_type_complex_4)
          IF (ASSOCIATED (area%d%c_sp)) THEN
             sizes(1) = SIZE (area%d%c_sp)
             valid = .TRUE.
          ENDIF
       CASE (dbcsr_type_real_8_2d)
          IF (ASSOCIATED (area%d%r2_dp)) THEN
             sizes(1) = SIZE (area%d%r2_dp,1)
             sizes(2) = SIZE (area%d%r2_dp,2)
             valid = .TRUE.
          ENDIF
       CASE (dbcsr_type_real_4_2d)
          IF (ASSOCIATED (area%d%r2_sp)) THEN
             sizes(1) = SIZE (area%d%r2_sp,1)
             sizes(2) = SIZE (area%d%r2_sp,2)
             valid = .TRUE.
          ENDIF
       CASE (dbcsr_type_complex_8_2d)
          IF (ASSOCIATED (area%d%c2_dp)) THEN
             sizes(1) = SIZE (area%d%c2_dp,1)
             sizes(2) = SIZE (area%d%c2_dp,2)
             valid = .TRUE.
          ENDIF
       CASE (dbcsr_type_complex_4_2d)
          IF (ASSOCIATED (area%d%c2_sp)) THEN
             sizes(1) = SIZE (area%d%c2_sp,1)
             sizes(2) = SIZE (area%d%c2_sp,2)
             valid = .TRUE.
          ENDIF
       CASE default
          CALL dbcsr_assert (.FALSE.,&
               dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Incorrect data type", __LINE__, error=error)
       END SELECT
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_get_sizes_any

! *****************************************************************************
!> \brief Returns the allocated data size
!> \param[in] area       data area to query for size, should be 2-D
!> \param[out] row_size  row size
!> \param[out] row_size  column size
!> \param[out] valid     whether the data is actually allocated
!> \param[in,out] error  error
! *****************************************************************************
  SUBROUTINE dbcsr_data_get_sizes_2 (area, row_size, col_size, valid, error)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER, INTENT(OUT)                     :: row_size, col_size
    LOGICAL, INTENT(OUT)                     :: valid
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_get_sizes_2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: d, error_handler

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)

    valid = .FALSE.
    row_size = 0
    col_size = 0
    IF (ASSOCIATED (area%d)) THEN
       CALL dbcsr_assert (dbcsr_type_is_2d (area%d%data_type),&
            dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "1-D data can not have column size", __LINE__, error=error)
       SELECT CASE (area%d%data_type)
       CASE (dbcsr_type_real_8_2d)
          IF (ASSOCIATED (area%d%r2_dp)) THEN
             row_size = SIZE (area%d%r2_dp,1)
             col_size = SIZE (area%d%r2_dp,2)
             valid = .TRUE.
          ENDIF
       CASE (dbcsr_type_real_4_2d)
          IF (ASSOCIATED (area%d%r2_sp)) THEN
             row_size = SIZE (area%d%r2_sp,1)
             col_size = SIZE (area%d%r2_sp,2)
             valid = .TRUE.
          ENDIF
       CASE (dbcsr_type_complex_8_2d)
          IF (ASSOCIATED (area%d%c2_dp)) THEN
             row_size = SIZE (area%d%c2_dp,1)
             col_size = SIZE (area%d%c2_dp,2)
              valid = .TRUE.
         ENDIF
       CASE (dbcsr_type_complex_4_2d)
          IF (ASSOCIATED (area%d%c2_sp)) THEN
             row_size = SIZE (area%d%c2_sp,1)
             col_size = SIZE (area%d%c2_sp,2)
             valid = .TRUE.
          ENDIF
       CASE default
          CALL dbcsr_assert (.FALSE.,&
               dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Incorrect or 1-D data type", __LINE__, error=error)
       END SELECT
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_get_sizes_2

! *****************************************************************************
!> \brief Returns the allocated data size
!> \param[in] area         data area to query for size
!> \param[out] total_size  size of array
!> \param[out] valid       whether the data is actually allocated
!> \param[in,out] error    error
! *****************************************************************************
  SUBROUTINE dbcsr_data_get_sizes_1 (area, total_size, valid, error)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER, INTENT(OUT)                     :: total_size
    LOGICAL, INTENT(OUT)                     :: valid
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_get_sizes_1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: d, error_handler

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)

    valid = .FALSE.
    total_size = 0
    IF (ASSOCIATED (area%d)) THEN
       CALL dbcsr_assert ("NOT", dbcsr_type_is_2d (area%d%data_type),&
            dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Should not use 2-D data", __LINE__, error=error)
       SELECT CASE (area%d%data_type)
       CASE (dbcsr_type_real_8)
          IF (ASSOCIATED (area%d%r_dp)) THEN
             total_size = SIZE (area%d%r_dp)
             valid = .TRUE.
          ENDIF
       CASE (dbcsr_type_real_4)
          IF (ASSOCIATED (area%d%r_sp)) THEN
             total_size = SIZE (area%d%r_sp)
             valid = .TRUE.
          ENDIF
       CASE (dbcsr_type_complex_8)
          IF (ASSOCIATED (area%d%c_dp)) THEN
             total_size = SIZE (area%d%c_dp)
             valid = .TRUE.
          ENDIF
       CASE (dbcsr_type_complex_4)
          IF (ASSOCIATED (area%d%c_sp)) THEN
             total_size = SIZE (area%d%c_sp)
             valid = .TRUE.
          ENDIF
       CASE default
          CALL dbcsr_assert (.FALSE.,&
               dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Incorrect data type", __LINE__, error=error)
       END SELECT
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_get_sizes_1


! *****************************************************************************
!> \brief Resizes the data stored in an object destructively.
!> \param[in,out] area   data area to resize
!> \param[in] sizes      array with the data sizes
!> \param[in,out] error  error
! *****************************************************************************
  SUBROUTINE dbcsr_data_resize (area, sizes, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, DIMENSION(:), INTENT(IN)        :: sizes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_resize', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: d, error_handler
    INTEGER, DIMENSION(2)                    :: old_sizes
    LOGICAL                                  :: old_valid

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    !
    CALL dbcsr_assert (ASSOCIATED (area%d), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN,&
         "Invalid data area", __LINE__, error=error)
    IF (dbcsr_type_is_2d (area%d%data_type)) THEN
       d = 2
       CALL dbcsr_assert (SIZE(sizes), "EQ", 2, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Sizes must have 2 elements for 2-D data", __LINE__, error=error)
    ELSE
       d = 1
       CALL dbcsr_assert (SIZE(sizes), "EQ", 1, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Sizes must have 1 elements for 1-D data", __LINE__, error=error)
    ENDIF
    CALL dbcsr_data_get_sizes (area, old_sizes(1:d),&
         valid=old_valid, error=error)
    IF (.NOT. old_valid .OR. ANY(old_sizes(1:d) .NE. sizes(1:d))) THEN
       IF (old_valid) THEN
          CALL internal_data_deallocate (area%d, data_type=area%d%data_type,&
               special_memory=area%d%special, error=error)
       ENDIF
       CALL internal_data_allocate (area%d, data_type=area%d%data_type,&
            sizes=sizes,&
            special_memory=area%d%special, error=error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_resize


! *****************************************************************************
!> \brief Returns an encapsulated scalar "1"
!> \param[in] data_type       use the data type
!> \result one                enpsulated value of one
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_one (data_type) RESULT (one)
    CHARACTER, INTENT(IN)                    :: data_type
    TYPE(dbcsr_scalar_type)                  :: one

    one%data_type = data_type
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       one%r_sp = 1.0_real_4
    CASE (dbcsr_type_real_8)
       one%r_dp = 1.0_real_8
    CASE (dbcsr_type_complex_4)
       one%c_sp = CMPLX(1.0, 0.0, real_4)
    CASE (dbcsr_type_complex_8)
       one%c_dp = CMPLX(1.0, 0.0, real_8)
    END SELECT
  END FUNCTION dbcsr_scalar_one

! *****************************************************************************
!> \brief Returns an encapsulated scalar "0"
!> \param[in] data_type       use the data type
!> \result zero               enpsulated value of zero
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_zero (data_type) RESULT (zero)
    CHARACTER, INTENT(IN)                    :: data_type
    TYPE(dbcsr_scalar_type)                  :: zero

    zero%data_type = data_type
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       zero%r_sp = 0.0_real_4
    CASE (dbcsr_type_real_8)
       zero%r_dp = 0.0_real_8
    CASE (dbcsr_type_complex_4)
       zero%c_sp = CMPLX(0.0, 0.0, real_4)
    CASE (dbcsr_type_complex_8)
       zero%c_dp = CMPLX(0.0, 0.0, real_8)
    END SELECT
  END FUNCTION dbcsr_scalar_zero


! *****************************************************************************
!> \brief Returns whether an encapsulated scalar is equal to another value
!> \param[in] s1                    one value
!> \param[in] s2                    another value
!> \result are_equal                whether values are equal
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_are_equal (s1, s2) RESULT (are_equal)
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: s1, s2
    LOGICAL                                  :: are_equal

    IF (s1%data_type .NE. s2%data_type) THEN
       are_equal = .FALSE.
    ELSE
       SELECT CASE (s1%data_type)
       CASE (dbcsr_type_real_4)
          are_equal = s1%r_sp .EQ. s2%r_sp
       CASE (dbcsr_type_real_8)
          are_equal = s1%r_dp .EQ. s2%r_dp
       CASE (dbcsr_type_complex_4)
          are_equal = s1%c_sp .EQ. s2%c_sp
       CASE (dbcsr_type_complex_8)
          are_equal = s1%c_dp .EQ. s2%c_dp
       CASE default
          are_equal = .FALSE.
       END SELECT
    ENDIF
  END FUNCTION dbcsr_scalar_are_equal

! *****************************************************************************
!> \brief Returns an encapsulated scalar as a negation of the given
!> \param[in] v                    given value
!> \result negated                 negated value
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_negative (s) RESULT (negated)
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: s
    TYPE(dbcsr_scalar_type)                  :: negated

    negated%data_type = s%data_type
    SELECT CASE (s%data_type)
    CASE (dbcsr_type_real_4)
       negated%r_sp = -s%r_sp
    CASE (dbcsr_type_real_8)
       negated%r_dp = -s%r_dp
    CASE (dbcsr_type_complex_4)
       negated%c_sp = -s%c_sp
    CASE (dbcsr_type_complex_8)
       negated%c_dp = -s%c_dp
    CASE default
       negated = dbcsr_scalar_zero (s%data_type)
    END SELECT
  END FUNCTION dbcsr_scalar_negative


! *****************************************************************************
!> \brief Checks whether the data type is 2-D.
!> \result Data type is 2-D.
! *****************************************************************************
  PURE FUNCTION dbcsr_type_is_2d (data_type)
    CHARACTER, INTENT(IN)                    :: data_type
    LOGICAL                                  :: dbcsr_type_is_2d

    dbcsr_type_is_2d = data_type .EQ. dbcsr_type_real_4_2d .OR.&
         data_type .EQ. dbcsr_type_real_8_2d .OR.&
         data_type .EQ. dbcsr_type_complex_4_2d .OR.&
         data_type .EQ. dbcsr_type_complex_8_2d
  END FUNCTION dbcsr_type_is_2d

! *****************************************************************************
!> \brief Returns 1-d data type corresponding to the given 2-D one.
!> \result 1-D data type
! *****************************************************************************
  PURE FUNCTION dbcsr_type_2d_to_1d (data_type)
    CHARACTER, INTENT(IN)                    :: data_type
    CHARACTER                                :: dbcsr_type_2d_to_1d

    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4_2d)
       dbcsr_type_2d_to_1d = dbcsr_type_real_4
    CASE (dbcsr_type_real_8_2d)
       dbcsr_type_2d_to_1d = dbcsr_type_real_8
    CASE (dbcsr_type_complex_4_2d)
       dbcsr_type_2d_to_1d = dbcsr_type_complex_4
    CASE (dbcsr_type_complex_8_2d)
       dbcsr_type_2d_to_1d = dbcsr_type_complex_8
    CASE (dbcsr_type_real_4)
       dbcsr_type_2d_to_1d = dbcsr_type_real_4
    CASE (dbcsr_type_real_8)
       dbcsr_type_2d_to_1d = dbcsr_type_real_8
    CASE (dbcsr_type_complex_4)
       dbcsr_type_2d_to_1d = dbcsr_type_complex_4
    CASE (dbcsr_type_complex_8)
       dbcsr_type_2d_to_1d = dbcsr_type_complex_8
    CASE default
       dbcsr_type_2d_to_1d = '0'
    END SELECT
  END FUNCTION dbcsr_type_2d_to_1d

! *****************************************************************************
!> \brief Returns 2-D data type corresponding to the given 1-D one.
!> \result 2-D data type
! *****************************************************************************
  PURE FUNCTION dbcsr_type_1d_to_2d (data_type)
    CHARACTER, INTENT(IN)                    :: data_type
    CHARACTER                                :: dbcsr_type_1d_to_2d

    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       dbcsr_type_1d_to_2d = dbcsr_type_real_4_2d
    CASE (dbcsr_type_real_8)
       dbcsr_type_1d_to_2d = dbcsr_type_real_8_2d
    CASE (dbcsr_type_complex_4)
       dbcsr_type_1d_to_2d = dbcsr_type_complex_4_2d
    CASE (dbcsr_type_complex_8)
       dbcsr_type_1d_to_2d = dbcsr_type_complex_8_2d
    CASE (dbcsr_type_real_4_2d)
       dbcsr_type_1d_to_2d = dbcsr_type_real_4_2d
    CASE (dbcsr_type_real_8_2d)
       dbcsr_type_1d_to_2d = dbcsr_type_real_8_2d
    CASE (dbcsr_type_complex_4_2d)
       dbcsr_type_1d_to_2d = dbcsr_type_complex_4_2d
    CASE (dbcsr_type_complex_8_2d)
       dbcsr_type_1d_to_2d = dbcsr_type_complex_8_2d
    CASE default
       dbcsr_type_1d_to_2d = '0'
    END SELECT
  END FUNCTION dbcsr_type_1d_to_2d



! *****************************************************************************
!> \brief Get actual data storage used for matrix
!> \param[in] area            Count data of this matrix
!> \result data_size          Data size used by matrix
! *****************************************************************************
  PURE FUNCTION dbcsr_data_get_size_referenced (area) RESULT (data_size_referenced)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER                                  :: data_size_referenced

    IF (ASSOCIATED (area%d)) THEN
       data_size_referenced = area%d%ref_size
    ELSE
       data_size_referenced = 0
    ENDIF
  END FUNCTION dbcsr_data_get_size_referenced

! *****************************************************************************
!> \brief Sets the referenced size of the data area
!> \param[in,out] data_area  area for which to set referenced data size
!> \param[in]                set referenced data size to this value
! *****************************************************************************
  PURE SUBROUTINE dbcsr_data_set_size_referenced (data_area, referenced_size)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: data_area
    INTEGER, INTENT(IN)                      :: referenced_size

    data_area%d%ref_size = referenced_size
  END SUBROUTINE dbcsr_data_set_size_referenced






#include "dbcsr_data_methods_d.F"
#include "dbcsr_data_methods_z.F"
#include "dbcsr_data_methods_s.F"
#include "dbcsr_data_methods_c.F"
END MODULE dbcsr_data_methods
