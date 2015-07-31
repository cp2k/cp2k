!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
MODULE dbcsr_data_methods_low

  USE acc_devmem,                      ONLY: acc_devmem_allocate_bytes,&
                                             acc_devmem_allocated,&
                                             acc_devmem_deallocate
  USE acc_event,                       ONLY: acc_event_create,&
                                             acc_event_destroy
  USE dbcsr_data_types,                ONLY: &
       dbcsr_data_area_type, dbcsr_data_obj, dbcsr_datatype_sizeof, &
       dbcsr_memtype_type, dbcsr_scalar_type, dbcsr_type_complex_4, &
       dbcsr_type_complex_4_2d, dbcsr_type_complex_8, &
       dbcsr_type_complex_8_2d, dbcsr_type_int_4, dbcsr_type_int_8, &
       dbcsr_type_real_4, dbcsr_type_real_4_2d, dbcsr_type_real_8, &
       dbcsr_type_real_8_2d
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_failure_level, dbcsr_fatal_level, &
       dbcsr_internal_error, dbcsr_warning_level, dbcsr_wrong_args_error
  USE dbcsr_ptr_util,                  ONLY: memory_allocate,&
                                             memory_deallocate,&
                                             memory_zero,&
                                             pointer_rank_remap2
  USE kinds,                           ONLY: dp,&
                                             real_4,&
                                             real_8

  !$ USE OMP_LIB

  IMPLICIT NONE


  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_data_methods_low'


  PUBLIC :: dbcsr_type_is_2d, dbcsr_type_2d_to_1d, dbcsr_type_1d_to_2d
  PUBLIC :: dbcsr_scalar, dbcsr_scalar_one, dbcsr_scalar_zero,&
            dbcsr_scalar_are_equal, dbcsr_scalar_negative,&
            dbcsr_scalar_get_type, dbcsr_scalar_set_type,&
            dbcsr_scalar_fill_all, dbcsr_scalar_get_value,&
            dbcsr_data_valid
  PUBLIC :: dbcsr_data_init, dbcsr_data_hold,&
            dbcsr_data_get_size, dbcsr_data_get_type
  PUBLIC :: dbcsr_get_data,&
            dbcsr_data_set_pointer,&
            dbcsr_data_clear_pointer,&
            dbcsr_data_get_sizes, dbcsr_data_verify_bounds,&
            dbcsr_data_exists, dbcsr_data_get_memory_type
  PUBLIC :: dbcsr_data_zero
  PUBLIC :: dbcsr_data_set_size_referenced, dbcsr_data_get_size_referenced
  PUBLIC :: dbcsr_get_data_p, dbcsr_get_data_p_s, dbcsr_get_data_p_c,&
            dbcsr_get_data_p_d, dbcsr_get_data_p_z,&
            dbcsr_get_data_p_2d_s, dbcsr_get_data_p_2d_d,&
            dbcsr_get_data_p_2d_c, dbcsr_get_data_p_2d_z

  PUBLIC :: internal_data_allocate, internal_data_deallocate

!> \brief Encapsulates a scalar.
  INTERFACE dbcsr_scalar
     MODULE PROCEDURE dbcsr_scalar_s, dbcsr_scalar_d,&
                      dbcsr_scalar_c, dbcsr_scalar_z
  END INTERFACE

!> \brief Encapsulates a scalar.
  INTERFACE dbcsr_scalar_get_value
     MODULE PROCEDURE dbcsr_scalar_get_value_s, dbcsr_scalar_get_value_d,&
                      dbcsr_scalar_get_value_c, dbcsr_scalar_get_value_z
  END INTERFACE

  INTERFACE dbcsr_data_set_pointer
     MODULE PROCEDURE set_data_p_s, set_data_p_d, set_data_p_c, set_data_p_z
     MODULE PROCEDURE set_data_p_2d_s, set_data_p_2d_d,&
                      set_data_p_2d_c, set_data_p_2d_z
     MODULE PROCEDURE set_data_area_area
  END INTERFACE

  INTERFACE dbcsr_get_data
     MODULE PROCEDURE get_data_s, get_data_d, get_data_c, get_data_z
     MODULE PROCEDURE get_data_2d_s, get_data_2d_d, get_data_2d_c, get_data_2d_z
  END INTERFACE


  INTERFACE dbcsr_get_data_p
     MODULE PROCEDURE dbcsr_get_data_c_s, dbcsr_get_data_c_c,&
                      dbcsr_get_data_c_d, dbcsr_get_data_c_z
  END INTERFACE

  INTERFACE dbcsr_get_data_cptr
     MODULE PROCEDURE dbcsr_get_data_c_s, dbcsr_get_data_c_c,&
                      dbcsr_get_data_c_d, dbcsr_get_data_c_z
  END INTERFACE

  INTERFACE dbcsr_data_get_sizes
     MODULE PROCEDURE dbcsr_data_get_sizes_any
     MODULE PROCEDURE dbcsr_data_get_sizes_1, dbcsr_data_get_sizes_2
  END INTERFACE

  LOGICAL, PARAMETER :: careful_mod = .FALSE.
  LOGICAL, PARAMETER :: debug_mod = .FALSE.

  CONTAINS

! *****************************************************************************
!> \brief Returns data type of a data area
!> \param[in] area         data area
!> \retval data_type       data type of the data area
! *****************************************************************************
  PURE FUNCTION dbcsr_data_get_type (area) RESULT (data_type)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER                                  :: data_type

    data_type = area%d%data_type
  END FUNCTION dbcsr_data_get_type


! *****************************************************************************
!> \brief ...
!> \param area ...
!> \retval memory_type ...
! *****************************************************************************
  FUNCTION dbcsr_data_get_memory_type (area) RESULT (memory_type)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    TYPE(dbcsr_memtype_type)                 :: memory_type

    memory_type = area%d%memory_type
  END FUNCTION dbcsr_data_get_memory_type

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
!> \param[in] sizes           sizes to allocate to
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE internal_data_allocate (area, sizes, error)
    TYPE(dbcsr_data_area_type), &
      INTENT(INOUT)                          :: area
    INTEGER, DIMENSION(:), INTENT(IN)        :: sizes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_data_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    IF (debug_mod) &
         WRITE(*,*)routineN//" Setting to sizes", sizes
    IF (dbcsr_type_is_2d (area%data_type)) THEN
       CALL dbcsr_assert (SIZE(sizes), "EQ", 2, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Sizes must have 2 elements for 2-D data", __LINE__, error=error)
    ELSE
       CALL dbcsr_assert (SIZE(sizes), "EQ", 1, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Sizes must have 1 elements for 1-D data", __LINE__, error=error)
    ENDIF

    SELECT CASE (area%data_type)
       CASE (dbcsr_type_int_4)
          CALL memory_allocate(area%i4, n=sizes(1), mem_type=area%memory_type, error=error)
       CASE (dbcsr_type_int_8)
          CALL memory_allocate(area%i8, n=sizes(1), mem_type=area%memory_type, error=error)
       CASE (dbcsr_type_real_4)
          CALL memory_allocate(area%r_sp, n=sizes(1), mem_type=area%memory_type, error=error)
       CASE (dbcsr_type_real_8)
          CALL memory_allocate(area%r_dp, n=sizes(1), mem_type=area%memory_type, error=error)
       CASE (dbcsr_type_complex_4)
          CALL memory_allocate(area%c_sp, n=sizes(1), mem_type=area%memory_type, error=error)
       CASE (dbcsr_type_complex_8)
          CALL memory_allocate(area%c_dp, n=sizes(1), mem_type=area%memory_type, error=error)

       CASE (dbcsr_type_real_4_2d)
          CALL memory_allocate(area%r2_sp, sizes=sizes, mem_type=area%memory_type, error=error)
       CASE (dbcsr_type_real_8_2d)
          CALL memory_allocate(area%r2_dp, sizes=sizes, mem_type=area%memory_type, error=error)
       CASE (dbcsr_type_complex_4_2d)
          CALL memory_allocate(area%c2_sp, sizes=sizes, mem_type=area%memory_type, error=error)
       CASE (dbcsr_type_complex_8_2d)
          CALL memory_allocate(area%c2_dp, sizes=sizes, mem_type=area%memory_type, error=error)

       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type.",__LINE__,error)
       END SELECT



    IF(area%memory_type%acc_devalloc) THEN
       IF(sizes(1) > 1) &
          CALL acc_devmem_allocate_bytes(area%acc_devmem, dbcsr_datatype_sizeof(area%data_type)*sizes(1))
       CALL acc_event_create(area%acc_ready)
    ENDIF

    CALL dbcsr_error_stop(error_handle, error)

  END SUBROUTINE internal_data_allocate


! *****************************************************************************
!> \brief Allocates pointers in the data type
!> \param[in,out] area        internal structure holding array pointers
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE internal_data_deallocate (area, error)
    TYPE(dbcsr_data_area_type), &
      INTENT(INOUT)                          :: area
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_data_deallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)

    SELECT CASE (area%data_type)
       CASE (dbcsr_type_int_4)
          CALL memory_deallocate(area%i4, mem_type=area%memory_type, error=error)
          NULLIFY (area%i4)
       CASE (dbcsr_type_int_8)
          CALL memory_deallocate(area%i8, mem_type=area%memory_type, error=error)
          NULLIFY (area%i8)
       CASE (dbcsr_type_real_4)
          CALL memory_deallocate(area%r_sp, mem_type=area%memory_type, error=error)
          NULLIFY (area%r_sp)
       CASE (dbcsr_type_real_8)
          CALL memory_deallocate(area%r_dp, mem_type=area%memory_type, error=error)
          NULLIFY (area%r_dp)
       CASE (dbcsr_type_complex_4)
          CALL memory_deallocate(area%c_sp, mem_type=area%memory_type, error=error)
          NULLIFY (area%c_sp)
       CASE (dbcsr_type_complex_8)
          CALL memory_deallocate(area%c_dp, mem_type=area%memory_type, error=error)
          NULLIFY (area%c_dp)

       CASE (dbcsr_type_real_4_2d)
          CALL memory_deallocate(area%r2_sp, mem_type=area%memory_type, error=error)
          NULLIFY (area%r2_sp)
       CASE (dbcsr_type_real_8_2d)
          CALL memory_deallocate(area%r2_dp, mem_type=area%memory_type, error=error)
          NULLIFY (area%r2_dp)
       CASE (dbcsr_type_complex_4_2d)
          CALL memory_deallocate(area%c2_sp, mem_type=area%memory_type, error=error)
          NULLIFY (area%c2_sp)
       CASE (dbcsr_type_complex_8_2d)
          CALL memory_deallocate(area%c2_dp, mem_type=area%memory_type, error=error)
          NULLIFY (area%c2_dp)

       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type.",__LINE__,error)
    END SELECT

    IF(area%memory_type%acc_devalloc) THEN
       IF(acc_devmem_allocated(area%acc_devmem)) &
          CALL acc_devmem_deallocate(area%acc_devmem)
       CALL acc_event_destroy(area%acc_ready)
    ENDIF

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE internal_data_deallocate



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
    CASE (dbcsr_type_int_4)
       NULLIFY (area%d%i4)
    CASE (dbcsr_type_int_8)
       NULLIFY (area%d%i8)
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
!> \brief Checks whether a data area is valid
!> \param[in] area         data area
!> \retval valid           whether the data area is valid
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
!> \retval valid           whether the data pointer exists 
! *****************************************************************************
  FUNCTION dbcsr_data_exists (area, error) RESULT (valid)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    LOGICAL                                  :: valid

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_exists', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    IF (careful_mod) THEN
       CALL dbcsr_error_set (routineN, error_handle, error)
    ENDIF
    !
    valid = dbcsr_data_valid (area)
    CALL dbcsr_assert (valid, dbcsr_fatal_level, dbcsr_wrong_args_error,&
         routineN, "Data area is invalid.", __LINE__, error)

    SELECT CASE (area%d%data_type)
       CASE (dbcsr_type_int_4)
          valid = ASSOCIATED (area%d%i4)
       CASE (dbcsr_type_int_8)
          valid = ASSOCIATED (area%d%i8)
       CASE (dbcsr_type_real_4)
          valid = ASSOCIATED (area%d%r_sp)
       CASE (dbcsr_type_real_8)
          valid = ASSOCIATED (area%d%r_dp)
       CASE (dbcsr_type_complex_4)
          valid = ASSOCIATED (area%d%c_sp)
       CASE (dbcsr_type_complex_8)
          valid = ASSOCIATED (area%d%c_dp)
       CASE (dbcsr_type_real_4_2d)
          valid = ASSOCIATED (area%d%r2_sp)
       CASE (dbcsr_type_real_8_2d)
          valid = ASSOCIATED (area%d%r2_dp)
       CASE (dbcsr_type_complex_4_2d)
          valid = ASSOCIATED (area%d%c2_sp)
       CASE (dbcsr_type_complex_8_2d)
          valid = ASSOCIATED (area%d%c2_dp)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type.",__LINE__,error)
    END SELECT

    IF (careful_mod) THEN
       CALL dbcsr_error_stop(error_handle, error)
    ENDIF
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

    IF (careful_mod) THEN
       CALL dbcsr_assert (ASSOCIATED (area%d),&
            dbcsr_caller_error, dbcsr_fatal_level,&
            routineN, "Can't hold an empty data area.",__LINE__,error)
       CALL dbcsr_assert (area%d%refcount, "GT", 0,&
            dbcsr_caller_error, dbcsr_fatal_level,&
            routineN, "Should not hold an area with zero references.",&
            __LINE__,error)
    ENDIF
    IF (.NOT. ASSOCIATED (area%d)) THEN
       RETURN
    ENDIF
    area%d%refcount = area%d%refcount + 1
  END SUBROUTINE dbcsr_data_hold



! *****************************************************************************
!> \brief Points data area data pointers to another data area
!>
!> Assumes that no memory will be lost when repointing the pointer in the data
!> area and that the area is initialized.
!> \param[in,out] area         data area to repoint
!> \param[in] rsize size of data area to point to
!> \param[in] csize size of data area to point to
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

    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: c_sp
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: c_dp
    INTEGER                                  :: bp, dt1, dt2, nze
    LOGICAL                                  :: compatible, pointee_is_2d, rmp
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r_sp
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r_dp
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    bp = 1 ; IF (PRESENT (source_lb)) bp = source_lb
    nze = rsize*csize
    dt1 = area%d%data_type
    dt2 = pointee%d%data_type
    IF (careful_mod) THEN
        compatible = dt1 .EQ. dt2 .OR. dt1 .EQ. dbcsr_type_1d_to_2d (dt2)
        CALL dbcsr_assert (compatible, dbcsr_fatal_level,&
             dbcsr_wrong_args_error, routineN,&
             "Can not point 1-d pointer to 2-d data",__LINE__,error)

    ENDIF
    pointee_is_2d = dbcsr_type_is_2d (dt2)
    IF (careful_mod) THEN
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
    ENDIF
    ! There's a remap if the ranks are compatible but not equal.
    rmp = dt1 .NE. dt2
    SELECT CASE (dt2)
    CASE (dbcsr_type_real_4_2d)
       area%d%r2_sp => pointee%d%r2_sp(1:rsize,1:csize)
    CASE (dbcsr_type_real_4)
       IF (rmp) THEN
          r_sp => dbcsr_get_data_p_s (pointee, bp, bp+nze-1)
          CALL pointer_rank_remap2(area%d%r2_sp, rsize, csize,&
               r_sp)
       ELSE
          area%d%r_sp => dbcsr_get_data_p_s (pointee, bp, bp+nze-1)
       ENDIF
    CASE (dbcsr_type_real_8_2d)
       area%d%r2_dp => pointee%d%r2_dp(1:rsize,1:csize)
    CASE (dbcsr_type_real_8)
       IF (rmp) THEN
          r_dp => dbcsr_get_data_p_d (pointee, bp, bp+nze-1)
          CALL pointer_rank_remap2(area%d%r2_dp, rsize, csize,&
               r_dp)
       ELSE
          area%d%r_dp => dbcsr_get_data_p_d (pointee, bp, bp+nze-1)
       ENDIF
    CASE (dbcsr_type_complex_4_2d)
       area%d%c2_sp => pointee%d%c2_sp(1:rsize,1:csize)
    CASE (dbcsr_type_complex_4)
       IF (rmp) THEN
          c_sp => dbcsr_get_data_p_c (pointee, bp, bp+nze-1)
          CALL pointer_rank_remap2(area%d%c2_sp, rsize, csize,&
               c_sp)
       ELSE
          area%d%c_sp => dbcsr_get_data_p_c (pointee, bp, bp+nze-1)
       ENDIF
    CASE (dbcsr_type_complex_8_2d)
       area%d%c2_dp => pointee%d%c2_dp(1:rsize,1:csize)
    CASE (dbcsr_type_complex_8)
       IF (rmp) THEN
          c_dp => dbcsr_get_data_p_z (pointee, bp, bp+nze-1)
          CALL pointer_rank_remap2(area%d%c2_dp, rsize, csize,&
               c_dp)
       ELSE
          area%d%c_dp => dbcsr_get_data_p_z (pointee, bp, bp+nze-1)
       ENDIF
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
    CALL dbcsr_data_set_size_referenced (area, rsize*csize)
    IF (debug_mod) THEN
       CALL dbcsr_assert (dbcsr_data_get_size_referenced (area), "EQ",&
            dbcsr_data_get_size (area), dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Size mismatch", __LINE__, error)
    ENDIF
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
          CASE (dbcsr_type_int_4)
             IF (ASSOCIATED (area%d%i4))&
                  data_size = SIZE (area%d%i4)
          CASE (dbcsr_type_int_8)
             IF (ASSOCIATED (area%d%i8))&
                  data_size = SIZE (area%d%i8)
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

    INTEGER                                  :: data_type, error_handler

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
!> \brief Clears a part of the data area
!> \param[in] area         data area
!> \param lb ...
!> \param ub ...
!> \param[in,out] error    error
!> \note Optimized for clearing big 1-D data areas from all data types.
! *****************************************************************************
  SUBROUTINE dbcsr_data_zero (area, lb, ub, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, DIMENSION(:), INTENT(in)        :: lb, ub
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_zero', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    IF (careful_mod) THEN
       CALL dbcsr_error_set (routineN, error_handle, error)
    ENDIF

    SELECT CASE (area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL memory_zero (area%d%r_sp(lb(1):), ub(1))
       CASE (dbcsr_type_real_8)
          CALL memory_zero (area%d%r_dp(lb(1):), ub(1))
       CASE (dbcsr_type_complex_4)
          CALL memory_zero (area%d%c_sp(lb(1):), ub(1))
       CASE (dbcsr_type_complex_8)
          CALL memory_zero (area%d%c_dp(lb(1):), ub(1))
       CASE (dbcsr_type_real_4_2d)
          area%d%r2_sp(lb(1):ub(1), lb(2):ub(2)) = 0.0_real_4
       CASE (dbcsr_type_real_8_2d)
          area%d%r2_dp(lb(1):ub(1), lb(2):ub(2)) = 0.0_real_8
       CASE (dbcsr_type_complex_4_2d)
          area%d%c2_sp(lb(1):ub(1), lb(2):ub(2)) = 0.0_real_4
       CASE (dbcsr_type_complex_8_2d)
          area%d%c2_dp(lb(1):ub(1), lb(2):ub(2)) = 0.0_real_8
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type.",__LINE__,error)
    END SELECT

    IF(area%d%memory_type%acc_devalloc) THEN
       WRITE(*,*) "dbcsr_data_zero: not yet supported for acc devmem"
       STOP 1
    ENDIF

    IF (careful_mod) THEN
       CALL dbcsr_error_stop(error_handle, error)
    ENDIF
  END SUBROUTINE dbcsr_data_zero

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

    INTEGER                                  :: error_handler

!   ---------------------------------------------------------------------------

    IF (careful_mod) &
         CALL dbcsr_error_set (routineN, error_handler, error)

    valid = .FALSE.
    sizes(:) = 0
    IF (ASSOCIATED (area%d)) THEN
       IF (careful_mod) THEN
          IF (dbcsr_type_is_2d (area%d%data_type)) THEN
             CALL dbcsr_assert (SIZE(sizes), "EQ", 2, dbcsr_fatal_level,&
                  dbcsr_wrong_args_error, routineN,&
                  "Sizes must have 2 elements for 2-D data", __LINE__, error=error)
          ELSE
             CALL dbcsr_assert (SIZE(sizes), "EQ", 1, dbcsr_fatal_level,&
                  dbcsr_wrong_args_error, routineN,&
                  "Sizes must have 1 elements for 1-D data", __LINE__, error=error)
          ENDIF
       ENDIF
       valid = dbcsr_data_exists (area, error=error)
       IF (valid) THEN
          SELECT CASE (area%d%data_type)
             CASE (dbcsr_type_real_8)
                   sizes(1) = SIZE (area%d%r_dp)
             CASE (dbcsr_type_real_4)
                   sizes(1) = SIZE (area%d%r_sp)
             CASE (dbcsr_type_complex_8)
                   sizes(1) = SIZE (area%d%c_dp)
             CASE (dbcsr_type_complex_4)
                   sizes(1) = SIZE (area%d%c_sp)
             CASE (dbcsr_type_real_8_2d)
                   sizes(1) = SIZE (area%d%r2_dp,1)
                   sizes(2) = SIZE (area%d%r2_dp,2)
             CASE (dbcsr_type_real_4_2d)
                   sizes(1) = SIZE (area%d%r2_sp,1)
                   sizes(2) = SIZE (area%d%r2_sp,2)
             CASE (dbcsr_type_complex_8_2d)
                   sizes(1) = SIZE (area%d%c2_dp,1)
                   sizes(2) = SIZE (area%d%c2_dp,2)
             CASE (dbcsr_type_complex_4_2d)
                   sizes(1) = SIZE (area%d%c2_sp,1)
                   sizes(2) = SIZE (area%d%c2_sp,2)
             CASE default
                CALL dbcsr_assert (.FALSE.,&
                     dbcsr_fatal_level, dbcsr_caller_error,&
                     routineN, "Incorrect data type", __LINE__, error=error)
          END SELECT
       ENDIF
    ENDIF
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_get_sizes_any

! *****************************************************************************
!> \brief Returns the allocated data size
!> \param[in] area       data area to query for size, should be 2-D
!> \param[out] row_size  row size
!> \param[out] col_size  column size 
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

    INTEGER                                  :: error_handler
    INTEGER, DIMENSION(2)                    :: s

!   ---------------------------------------------------------------------------

    IF (careful_mod) &
       CALL dbcsr_error_set (routineN, error_handler, error)
    IF (ASSOCIATED (area%d)) THEN
       IF (careful_mod) &
            CALL dbcsr_assert (dbcsr_type_is_2d (area%d%data_type),&
            dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "1-D data can not have column size", __LINE__, error=error)
       CALL dbcsr_data_get_sizes_any (area, s, valid, error=error)
       row_size = s(1)
       col_size = s(2)
    ELSE
       valid = .FALSE.
       row_size = 0
       col_size = 0
    ENDIF
    IF (careful_mod) &
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

    INTEGER                                  :: error_handler
    INTEGER, DIMENSION(1)                    :: s

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)

    IF (ASSOCIATED (area%d)) THEN
       IF (careful_mod) &
            CALL dbcsr_assert ("NOT", dbcsr_type_is_2d (area%d%data_type),&
            dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Should not use 2-D data", __LINE__, error=error)
       CALL dbcsr_data_get_sizes_any (area, s, valid, error=error)
       total_size = s(1)
    ELSE
       valid = .FALSE.
       total_size = 0
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_get_sizes_1


! *****************************************************************************
!> \brief Returns an encapsulated scalar "1"
!> \param[in] data_type       use the data type
!> \retval one                encapsulated value of one 
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_one (data_type) RESULT (one)
    INTEGER, INTENT(IN)                      :: data_type
    TYPE(dbcsr_scalar_type)                  :: one

    one = dbcsr_scalar_zero (data_type)
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
!> \brief Returns an encapsulated scalar "sqrt(-1)"
!> \param[in] data_type       use the data type
!> \retval i                 encapsulated value of i
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_i (data_type) RESULT (i)
    INTEGER, INTENT(IN)                      :: data_type
    TYPE(dbcsr_scalar_type)                  :: i

    i = dbcsr_scalar_zero (data_type)
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       i%r_sp = 0.0_real_4
    CASE (dbcsr_type_real_8)
       i%r_dp = 0.0_real_8
    CASE (dbcsr_type_complex_4)
       i%c_sp = CMPLX(0.0, 1.0, real_4)
    CASE (dbcsr_type_complex_8)
       i%c_dp = CMPLX(0.0, 1.0, real_8)
    END SELECT
  END FUNCTION dbcsr_scalar_i

! *****************************************************************************
!> \brief Returns an encapsulated scalar "0"
!> \param[in] data_type       use the data type
!> \retval zero               encapsulated value of zero
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_zero (data_type) RESULT (zero)
    INTEGER, INTENT(IN)                      :: data_type
    TYPE(dbcsr_scalar_type)                  :: zero

    zero%data_type = data_type
    zero%r_sp = 0.0_real_4
    zero%r_dp = 0.0_real_8
    zero%c_sp = CMPLX(0.0, 0.0, real_4)
    zero%c_dp = CMPLX(0.0, 0.0, real_8)
  END FUNCTION dbcsr_scalar_zero


! *****************************************************************************
!> \brief Returns whether an encapsulated scalar is equal to another value
!> \param[in] s1                    one value
!> \param[in] s2                    another value
!> \retval are_equal                whether values are equal 
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
!> \param[in] s                    given value
!> \retval negated                 negated value 
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_negative (s) RESULT (negated)
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: s
    TYPE(dbcsr_scalar_type)                  :: negated

    negated = dbcsr_scalar_zero (s%data_type)
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
!> \brief Returns an encapsulated scalar as a negation of the given
!> \param s1 ...
!> \param s2 ...
!> \retval s_sum ...
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_add (s1, s2) RESULT (s_sum)
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: s1, s2
    TYPE(dbcsr_scalar_type)                  :: s_sum

    s_sum = dbcsr_scalar_zero (s1%data_type)
    SELECT CASE (s1%data_type)
    CASE (dbcsr_type_real_4)
       s_sum%r_sp = s1%r_sp + s2%r_sp
    CASE (dbcsr_type_real_8)
       s_sum%r_dp = s1%r_dp + s2%r_dp
    CASE (dbcsr_type_complex_4)
       s_sum%c_sp = s1%c_sp + s2%c_sp
    CASE (dbcsr_type_complex_8)
       s_sum%c_dp = s1%c_dp + s2%c_dp
    CASE default
       s_sum = dbcsr_scalar_zero (s1%data_type)
    END SELECT
  END FUNCTION dbcsr_scalar_add

! *****************************************************************************
!> \brief ...
!> \param s1 ...
!> \param s2 ...
!> \retval s_product ...
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_multiply (s1, s2) RESULT (s_product)
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: s1, s2
    TYPE(dbcsr_scalar_type)                  :: s_product

    s_product = dbcsr_scalar_zero (s1%data_type)
    SELECT CASE (s1%data_type)
    CASE (dbcsr_type_real_4)
       s_product%r_sp = s1%r_sp * s2%r_sp
    CASE (dbcsr_type_real_8)
       s_product%r_dp = s1%r_dp * s2%r_dp
    CASE (dbcsr_type_complex_4)
       s_product%c_sp = s1%c_sp * s2%c_sp
    CASE (dbcsr_type_complex_8)
       s_product%c_dp = s1%c_dp * s2%c_dp
    CASE default
       s_product = dbcsr_scalar_zero (s1%data_type)
    END SELECT
  END FUNCTION dbcsr_scalar_multiply


! *****************************************************************************
!> \brief Returns data type of a scalar
!> \param[in] scalar       scalar
!> \retval data_type       data type of the scalar
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_get_type (scalar) RESULT (data_type)
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: scalar
    INTEGER                                  :: data_type

    data_type = scalar%data_type
  END FUNCTION dbcsr_scalar_get_type

! *****************************************************************************
!> \brief Sets data type of a scalar
!> \param[in] scalar       scalar
!> \param data_type ...
! *****************************************************************************
  ELEMENTAL SUBROUTINE dbcsr_scalar_set_type (scalar, data_type)
    TYPE(dbcsr_scalar_type), INTENT(INOUT)   :: scalar
    INTEGER, INTENT(IN)                      :: data_type

    scalar%data_type = data_type
  END SUBROUTINE dbcsr_scalar_set_type


! *****************************************************************************
!> \brief Fills all data and precision types from the set one
!> \param[in,out] scalar         data area
! *****************************************************************************
  ELEMENTAL SUBROUTINE dbcsr_scalar_fill_all (scalar)
    TYPE(dbcsr_scalar_type), INTENT(INOUT)   :: scalar

    SELECT CASE(scalar%data_type)
       CASE (dbcsr_type_real_4)
          !scalar%r_sp = 0
          scalar%r_dp = REAL(scalar%r_sp, KIND=real_8)
          scalar%c_sp = CMPLX(scalar%r_sp, 0, KIND=real_4)
          scalar%c_dp = CMPLX(scalar%r_sp, 0, KIND=real_8)
       CASE (dbcsr_type_real_8)
          scalar%r_sp = REAL(scalar%r_dp, KIND=real_4)
          !scalar%r_dp = REAL(scalar%r_dp, KIND=real_8)
          scalar%c_sp = CMPLX(scalar%r_dp, 0, KIND=real_4)
          scalar%c_dp = CMPLX(scalar%r_dp, 0, KIND=real_8)
       CASE (dbcsr_type_complex_4)
          scalar%r_sp = REAL(scalar%c_sp, KIND=real_4)
          scalar%r_dp = REAL(scalar%c_sp, KIND=real_8)
          !scalar%c_sp = CMPLX(scalar%c_sp, KIND=real_4)
          scalar%c_dp = CMPLX(scalar%c_sp, KIND=real_8)
       CASE (dbcsr_type_complex_8)
          scalar%r_sp = REAL(scalar%c_dp, KIND=real_4)
          scalar%r_dp = REAL(scalar%c_dp, KIND=real_8)
          scalar%c_sp = CMPLX(scalar%c_dp, KIND=real_4)
          !scalar%c_dp = CMPLX(scalar%r_dp, KIND=real_8)
    END SELECT
  END SUBROUTINE dbcsr_scalar_fill_all

! *****************************************************************************
!> \brief Checks whether the data type is 2-D.
!> \param data_type ...
!> \retval dbcsr_type_is_2d Data type is 2-D.
! *****************************************************************************
  PURE FUNCTION dbcsr_type_is_2d (data_type)
    INTEGER, INTENT(IN)                      :: data_type
    LOGICAL                                  :: dbcsr_type_is_2d

    dbcsr_type_is_2d = data_type .EQ. dbcsr_type_real_4_2d .OR.&
         data_type .EQ. dbcsr_type_real_8_2d .OR.&
         data_type .EQ. dbcsr_type_complex_4_2d .OR.&
         data_type .EQ. dbcsr_type_complex_8_2d
  END FUNCTION dbcsr_type_is_2d

! *****************************************************************************
!> \brief Returns 1-d data type corresponding to the given 2-D one.
!> \param data_type ...
!> \retval dbcsr_type_2d_to_1d 1-D data type
! *****************************************************************************
  PURE FUNCTION dbcsr_type_2d_to_1d (data_type)
    INTEGER, INTENT(IN)                      :: data_type
    INTEGER                                  :: dbcsr_type_2d_to_1d

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
       dbcsr_type_2d_to_1d = -1
    END SELECT
  END FUNCTION dbcsr_type_2d_to_1d

! *****************************************************************************
!> \brief Returns 2-D data type corresponding to the given 1-D one.
!> \param data_type ...
!> \retval dbcsr_type_1d_to_2d 2-D data type 
! *****************************************************************************
  PURE FUNCTION dbcsr_type_1d_to_2d (data_type)
    INTEGER, INTENT(IN)                      :: data_type
    INTEGER                                  :: dbcsr_type_1d_to_2d

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
       dbcsr_type_1d_to_2d = -1
    END SELECT
  END FUNCTION dbcsr_type_1d_to_2d



! *****************************************************************************
!> \brief Get actual data storage used for matrix
!> \param[in] area            Count data of this matrix
!> \retval data_size_referenced Data size used by matrix
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
!> \param[in] referenced_size     set referenced data size to this value  
! *****************************************************************************
  PURE SUBROUTINE dbcsr_data_set_size_referenced (data_area, referenced_size)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: data_area
    INTEGER, INTENT(IN)                      :: referenced_size

    data_area%d%ref_size = referenced_size
  END SUBROUTINE dbcsr_data_set_size_referenced


! *****************************************************************************
#include "dbcsr_data_methods_low_d.f90"
#include "dbcsr_data_methods_low_z.f90"
#include "dbcsr_data_methods_low_s.f90"
#include "dbcsr_data_methods_low_c.f90"
END MODULE dbcsr_data_methods_low
