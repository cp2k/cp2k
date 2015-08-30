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
MODULE dbcsr_data_methods
  USE acc_devmem,                      ONLY: acc_devmem_allocate_bytes,&
                                             acc_devmem_allocated,&
                                             acc_devmem_dev2host,&
                                             acc_devmem_ensure_size_bytes,&
                                             acc_devmem_host2dev,&
                                             acc_devmem_setzero_bytes,&
                                             acc_devmem_size_in_bytes
  USE acc_event,                       ONLY: acc_event_record
  USE dbcsr_data_methods_low,          ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_exists, &
       dbcsr_data_get_memory_type, dbcsr_data_get_size, &
       dbcsr_data_get_size_referenced, dbcsr_data_get_sizes, &
       dbcsr_data_get_type, dbcsr_data_hold, dbcsr_data_init, &
       dbcsr_data_set_pointer, dbcsr_data_set_size_referenced, &
       dbcsr_data_valid, dbcsr_data_zero, dbcsr_get_data, dbcsr_get_data_p, &
       dbcsr_get_data_p_c, dbcsr_get_data_p_d, dbcsr_get_data_p_s, &
       dbcsr_get_data_p_z, dbcsr_scalar, dbcsr_scalar_are_equal, &
       dbcsr_scalar_fill_all, dbcsr_scalar_get_type, dbcsr_scalar_get_value, &
       dbcsr_scalar_negative, dbcsr_scalar_one, dbcsr_scalar_set_type, &
       dbcsr_scalar_zero, dbcsr_type_1d_to_2d, dbcsr_type_2d_to_1d, &
       dbcsr_type_is_2d, internal_data_allocate, internal_data_deallocate
  USE dbcsr_data_types,                ONLY: &
       dbcsr_data_obj, dbcsr_datatype_sizeof, dbcsr_memtype_default, &
       dbcsr_memtype_type, dbcsr_type_complex_4, dbcsr_type_complex_8, &
       dbcsr_type_int_4, dbcsr_type_int_8, dbcsr_type_real_4, &
       dbcsr_type_real_8
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_failure_level, dbcsr_fatal_level, &
       dbcsr_unimplemented_error_nr, dbcsr_warning_level, &
       dbcsr_wrong_args_error
  USE dbcsr_mem_methods,               ONLY: dbcsr_mempool_add,&
                                             dbcsr_mempool_get
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE kinds,                           ONLY: dp,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8

  IMPLICIT NONE


  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_data_methods'
  LOGICAL, PARAMETER :: careful_mod = .FALSE.

  INTEGER, SAVE                        :: id = 0

  PUBLIC :: dbcsr_type_2d_to_1d, dbcsr_type_1d_to_2d
  PUBLIC :: dbcsr_scalar, dbcsr_scalar_one, dbcsr_scalar_zero,&
            dbcsr_scalar_are_equal, dbcsr_scalar_negative,&
            dbcsr_scalar_get_type, dbcsr_scalar_set_type,&
            dbcsr_scalar_fill_all, dbcsr_scalar_get_value, &
            dbcsr_data_valid
  PUBLIC :: dbcsr_data_init, dbcsr_data_new, dbcsr_data_copy, dbcsr_data_hold,&
            dbcsr_data_release, dbcsr_data_get_size, dbcsr_data_get_type
  PUBLIC :: dbcsr_get_data,&
            dbcsr_data_set_pointer,&
            dbcsr_data_clear_pointer,&
            dbcsr_data_ensure_size,&
            dbcsr_data_get_sizes,&
            dbcsr_data_get_memory_type
  PUBLIC :: dbcsr_data_set_size_referenced, dbcsr_data_get_size_referenced
  PUBLIC :: dbcsr_get_data_p, dbcsr_get_data_p_s, dbcsr_get_data_p_c,&
            dbcsr_get_data_p_d, dbcsr_get_data_p_z
  PUBLIC :: dbcsr_data_host2dev, dbcsr_data_dev2host

  INTERFACE dbcsr_data_copy
     MODULE PROCEDURE dbcsr_data_copy_frombuffer_obj
     MODULE PROCEDURE dbcsr_data_copy_frombuffer_i
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief Transfers data from host- to device-buffer, asynchronously.
!> \param area ...
!> \param error ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_data_host2dev(area, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_host2dev', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: c_sp
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: c_dp
    INTEGER(KIND=int_4), DIMENSION(:), &
      POINTER                                :: i4
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: i8
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r_sp
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r_dp

    IF(.NOT. acc_devmem_allocated(area%d%acc_devmem)) RETURN !nothing to do
    IF(area%d%ref_size==0) RETURN !nothing to do

    SELECT CASE (area%d%data_type)
       CASE (dbcsr_type_int_4)
         i4 => area%d%i4(:area%d%ref_size)
         CALL acc_devmem_host2dev(area%d%acc_devmem, hostmem=i4,   stream=area%d%memory_type%acc_stream)
       CASE (dbcsr_type_int_8)
         i8 => area%d%i8(:area%d%ref_size)
         CALL acc_devmem_host2dev(area%d%acc_devmem, hostmem=i8,   stream=area%d%memory_type%acc_stream)
       CASE (dbcsr_type_real_4)
         r_sp => area%d%r_sp(:area%d%ref_size)
         CALL acc_devmem_host2dev(area%d%acc_devmem, hostmem=r_sp, stream=area%d%memory_type%acc_stream)
       CASE (dbcsr_type_real_8)
         r_dp => area%d%r_dp(:area%d%ref_size)
         CALL acc_devmem_host2dev(area%d%acc_devmem, hostmem=r_dp, stream=area%d%memory_type%acc_stream)
       CASE (dbcsr_type_complex_4)
         c_sp => area%d%c_sp(:area%d%ref_size)
         CALL acc_devmem_host2dev(area%d%acc_devmem, hostmem=c_sp, stream=area%d%memory_type%acc_stream)
       CASE (dbcsr_type_complex_8)
         c_dp => area%d%c_dp(:area%d%ref_size)
         CALL acc_devmem_host2dev(area%d%acc_devmem, hostmem=c_dp, stream=area%d%memory_type%acc_stream)
       CASE default
         CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type.",__LINE__,error)
    END SELECT

    CALL acc_event_record(area%d%acc_ready, area%d%memory_type%acc_stream)
  END SUBROUTINE dbcsr_data_host2dev


! *****************************************************************************
!> \brief Transfers data from device- to host-buffer, asynchronously.
!> \param area ...
!> \param error ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_data_dev2host(area, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_dev2host', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: c_sp
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: c_dp
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r_sp
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r_dp

    IF(area%d%ref_size==0) RETURN !nothing to do

    SELECT CASE (area%d%data_type)
       CASE (dbcsr_type_real_4)
         r_sp => area%d%r_sp(:area%d%ref_size)
         CALL acc_devmem_dev2host(area%d%acc_devmem, hostmem=r_sp, stream=area%d%memory_type%acc_stream)
       CASE (dbcsr_type_real_8)
         r_dp => area%d%r_dp(:area%d%ref_size)
         CALL acc_devmem_dev2host(area%d%acc_devmem, hostmem=r_dp, stream=area%d%memory_type%acc_stream)
       CASE (dbcsr_type_complex_4)
         c_sp => area%d%c_sp(:area%d%ref_size)
         CALL acc_devmem_dev2host(area%d%acc_devmem, hostmem=c_sp, stream=area%d%memory_type%acc_stream)
       CASE (dbcsr_type_complex_8)
         c_dp => area%d%c_dp(:area%d%ref_size)
         CALL acc_devmem_dev2host(area%d%acc_devmem, hostmem=c_dp, stream=area%d%memory_type%acc_stream)
       CASE default
         CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type.",__LINE__,error)
    END SELECT

  END SUBROUTINE dbcsr_data_dev2host


! *****************************************************************************
!> \brief Initializes a data area and all the actual data pointers
!> \param[inout] area         data area
!> \param[in] data_type       select data type to use
!> \param[in] data_size       (optional) allocate this much data
!> \param[in] data_size2      (optional) second dimension data size
!> \param[in] memory_type     (optional) type of memory to use
! *****************************************************************************
  SUBROUTINE dbcsr_data_new (area, data_type, data_size, data_size2,&
       memory_type)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, INTENT(IN)                      :: data_type
    INTEGER, INTENT(IN), OPTIONAL            :: data_size, data_size2
    TYPE(dbcsr_memtype_type), INTENT(IN), &
      OPTIONAL                               :: memory_type

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: d, error_handler, &
                                                total_size_oversized, &
                                                total_size_requested
    INTEGER, DIMENSION(2)                    :: sizes_oversized, &
                                                sizes_requested
    TYPE(dbcsr_error_type)                   :: error
    TYPE(dbcsr_memtype_type)                 :: my_memory_type

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert(.NOT.ASSOCIATED(area%d), dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "area already associcated", __LINE__, error=error)

    my_memory_type = dbcsr_memtype_default
    IF (PRESENT (memory_type)) my_memory_type = memory_type

    sizes_requested(:)=0; d=1
    IF (PRESENT (data_size)) sizes_requested(1)=data_size

    IF (dbcsr_type_is_2d(data_type)) THEN
       d=2
       IF (PRESENT(data_size2)) sizes_requested(2)=data_size2

       CALL dbcsr_assert(PRESENT(data_size).EQV.PRESENT(data_size2), &
               dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
               "Must specify 2 sizes for 2-D data", __LINE__, error=error)
    ENDIF

    sizes_oversized = INT(sizes_requested * my_memory_type%oversize_factor)
    total_size_requested = PRODUCT(sizes_requested(1:d))
    total_size_oversized = PRODUCT(sizes_oversized(1:d))

    IF(total_size_requested>1 .AND. ASSOCIATED(my_memory_type%pool)) THEN
        area = dbcsr_mempool_get(my_memory_type, data_type, total_size_requested, error)
    ENDIF

    IF(.NOT. ASSOCIATED(area%d)) THEN
       ALLOCATE(area%d)
       !$OMP CRITICAL (crit_area_id)
       id = id + 1
       area%d%id = id
       !$OMP END CRITICAL (crit_area_id)
       area%d%refcount = 1
       area%d%memory_type = my_memory_type
       area%d%data_type = data_type
       area%d%ref_size = 0
       IF(PRESENT(data_size)) THEN
          area%d%ref_size = total_size_oversized
          CALL internal_data_allocate (area%d, sizes_oversized(1:d), error=error)
       END IF
    ENDIF

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_new

! *****************************************************************************
!> \brief Copy data from an input template obj
!> \param[inout] area         data area
!> \param[in] template_area   template data area (dbcsr_data_obj)
!> \param[in] size_data       data area size
!> \param[in] displ           initial position in the data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_copy_frombuffer_obj (area, template_area, size_data, displ)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    TYPE(dbcsr_data_obj), INTENT(IN)         :: template_area
    INTEGER, INTENT(IN)                      :: size_data
    INTEGER(kind=int_8), INTENT(IN)          :: displ

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_data_copy_frombuffer_obj', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    INTEGER(kind=int_8)                      :: displ_end
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert(ASSOCIATED(area%d), dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN,&
         "area not associated", __LINE__, error=error)

    CALL dbcsr_assert(ASSOCIATED(template_area%d), dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN,&
         "template area not associated", __LINE__, error=error)

    CALL dbcsr_assert(area%d%data_type, "EQ", template_area%d%data_type, &
         dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN,&
         "area and template area type mismatch", __LINE__, error=error)

    displ_end = displ+size_data-1
    CALL dbcsr_assert (INT(displ_end), "LE", dbcsr_data_get_size_referenced (area),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Attempt to copy out of bounds",__LINE__,error)

    IF (size_data.gt.0) THEN
       SELECT CASE (area%d%data_type)
       CASE (dbcsr_type_int_8)
          area%d%i8(displ:displ_end) = template_area%d%i8(1:size_data)
       CASE (dbcsr_type_int_4)
          area%d%i4(displ:displ_end) = template_area%d%i4(1:size_data)
       CASE (dbcsr_type_real_8)
          area%d%r_dp(displ:displ_end) = template_area%d%r_dp(1:size_data)
       CASE (dbcsr_type_real_4)
          area%d%r_sp(displ:displ_end) = template_area%d%r_sp(1:size_data)
       CASE (dbcsr_type_complex_8)
          area%d%c_dp(displ:displ_end) = template_area%d%c_dp(1:size_data)
       CASE (dbcsr_type_complex_4)
          area%d%c_sp(displ:displ_end) = template_area%d%c_sp(1:size_data)
       CASE default
          CALL dbcsr_assert(.FALSE., dbcsr_failure_level,&
               dbcsr_unimplemented_error_nr, routineN,&
               "Invalid data type are supported",__LINE__,error)
       END SELECT
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_copy_frombuffer_obj

! *****************************************************************************
!> \brief Copy data from an input template int array
!> \param[inout] area         data area
!> \param[in] template_area   template data area (dbcsr_data_obj)
!> \param[in] size_data       data area size
!> \param[in] displ           initial position in the data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_copy_frombuffer_i (area, template_area, size_data, displ)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, DIMENSION(:), INTENT(IN)        :: template_area
    INTEGER, INTENT(IN)                      :: size_data
    INTEGER(kind=int_8), INTENT(IN)          :: displ

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_copy_frombuffer_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    INTEGER(kind=int_8)                      :: displ_end
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert(ASSOCIATED(area%d), dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN,&
         "area not associated", __LINE__, error=error)

    CALL dbcsr_assert(area%d%data_type, "EQ", dbcsr_type_int_4, &
         dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN,&
         "area and template area type mismatch", __LINE__, error=error)

    displ_end = displ+size_data-1
    CALL dbcsr_assert (INT(displ_end), "LE", dbcsr_data_get_size_referenced (area),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Attempt to copy out of bounds",__LINE__,error)

    area%d%i4(displ:displ_end) = template_area(1:size_data)

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_copy_frombuffer_i

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
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_data_ensure_size (area, data_size, nocopy, zero_pad, factor, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, INTENT(IN)                      :: data_size
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, zero_pad
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: factor
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_ensure_size', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: current_size, error_handler, &
                                                wanted_size
    LOGICAL                                  :: nocp, pad
    TYPE(dbcsr_data_obj)                     :: area_tmp

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert(ASSOCIATED (area%d), dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Data area must be setup.",__LINE__,error)
    current_size = dbcsr_data_get_size (area)

    ! allocate some more as padding for libsmm kernels which read over the end.
    IF (data_size.gt.0)  THEN
       wanted_size = data_size + 10
    ELSE
       wanted_size = data_size
    ENDIF

    !IF(area%d%memory_type%acc_devalloc) THEN
    !    IF(current_size==acc_devmem_size(area%d%acc_devmem)) &
    !      WRITE (*,*) "dbcsr_data_ensure_size: Host and device buffer differ in size."
    !END IF
    !   CALL dbcsr_assert(current_size==acc_devmem_size(area%d%acc_devmem),&
    !          dbcsr_fatal_level, dbcsr_caller_error,&
    !          routineN, "Host and device buffer differ in size.",__LINE__,error)

    CALL dbcsr_data_set_size_referenced (area, data_size)
    IF (current_size .GT. 1 .AND. current_size .GE. wanted_size) THEN
       IF (careful_mod) CALL dbcsr_error_stop(error_handler, error)
       RETURN
    ENDIF
    !
    nocp = .FALSE.
    IF (PRESENT (nocopy)) nocp = nocopy
    pad = .FALSE.
    IF (PRESENT (zero_pad)) pad = zero_pad

    IF(dbcsr_data_exists(area, error=error)) THEN
      IF(nocp .AND. dbcsr_data_get_size(area) <= 1) &
         CALL internal_data_deallocate(area%d, error)
    END IF

    IF (.NOT. dbcsr_data_exists (area, error=error)) THEN
        IF(ASSOCIATED(area%d%memory_type%pool)) THEN
           area_tmp = dbcsr_mempool_get(area%d%memory_type, area%d%data_type, wanted_size, error)
           IF(ASSOCIATED(area_tmp%d)) THEN
              area_tmp%d%ref_size = wanted_size
              area_tmp%d%refcount = area%d%refcount
              DEALLOCATE(area%d)
              area = area_tmp
           END IF
        END IF

        IF (.NOT. dbcsr_data_exists (area, error=error)) &
           CALL internal_data_allocate (area%d, (/ wanted_size /), error=error)

        IF(pad) CALL dbcsr_data_zero (area, (/ 1 /), (/ wanted_size /), error=error)
    ELSE
       SELECT CASE (area%d%data_type)
          CASE (dbcsr_type_int_8)
             CALL ensure_array_size(area%d%i8, ub=wanted_size,&
                  memory_type=area%d%memory_type,&
                  nocopy=nocp, zero_pad=zero_pad,&
                  factor=factor,error=error)
          CASE (dbcsr_type_int_4)
             CALL ensure_array_size(area%d%i4, ub=wanted_size,&
                  memory_type=area%d%memory_type,&
                  nocopy=nocp, zero_pad=zero_pad,&
                  factor=factor,error=error)
          CASE (dbcsr_type_real_8)
             CALL ensure_array_size(area%d%r_dp, ub=wanted_size,&
                  memory_type=area%d%memory_type,&
                  nocopy=nocp, zero_pad=zero_pad,&
                  factor=factor,error=error)
          CASE (dbcsr_type_real_4)
             CALL ensure_array_size (area%d%r_sp, ub=wanted_size,&
                  memory_type=area%d%memory_type,&
                  nocopy=nocp, zero_pad=zero_pad,&
                  factor=factor,error=error)
          CASE (dbcsr_type_complex_8)
             CALL ensure_array_size (area%d%c_dp, ub=wanted_size,&
                  memory_type=area%d%memory_type,&
                  nocopy=nocp, zero_pad=zero_pad,&
                  factor=factor,error=error)
          CASE (dbcsr_type_complex_4)
             CALL ensure_array_size (area%d%c_sp, ub=wanted_size,&
                  memory_type=area%d%memory_type,&
                  nocopy=nocp, zero_pad=zero_pad,&
                  factor=factor,error=error)
          CASE default
             CALL dbcsr_assert(.FALSE., dbcsr_failure_level,&
                  dbcsr_unimplemented_error_nr, routineN,&
                  "Invalid data type are supported",__LINE__,error)
       END SELECT


       IF(area%d%memory_type%acc_devalloc) THEN
          IF(.NOT.acc_devmem_allocated(area%d%acc_devmem)) THEN
             CALL acc_devmem_allocate_bytes(area%d%acc_devmem, &
                dbcsr_datatype_sizeof(area%d%data_type)*dbcsr_data_get_size(area))
             IF(pad) CALL acc_devmem_setzero_bytes(area%d%acc_devmem, stream=area%d%memory_type%acc_stream)
          ELSE
             CALL acc_devmem_ensure_size_bytes(area%d%acc_devmem, &
                   area%d%memory_type%acc_stream, &
                   dbcsr_datatype_sizeof(area%d%data_type)*dbcsr_data_get_size(area), &
                   nocopy, zero_pad)
          END IF
          CALL acc_event_record(area%d%acc_ready, area%d%memory_type%acc_stream)
          CALL dbcsr_assert(&
                 dbcsr_datatype_sizeof(area%d%data_type)*dbcsr_data_get_size(area)&
                  == acc_devmem_size_in_bytes(area%d%acc_devmem),&
                 dbcsr_fatal_level, dbcsr_caller_error,&
                 routineN, "Host and device buffer differ in size.",__LINE__,error)
       ENDIF

    ENDIF
    IF (careful_mod) CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_ensure_size



! *****************************************************************************
!> \brief Removes a reference and/or clears the data area.
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_release (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (ASSOCIATED (area%d), &
         dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "Data seems to be unreferenced.",__LINE__,error)
    IF (ASSOCIATED (area%d)) THEN
       !
       IF (careful_mod) &
            CALL dbcsr_assert (area%d%refcount, "GT", 0,&
            dbcsr_warning_level, dbcsr_caller_error,&
            routineN, "Data seems to be unreferenced.",__LINE__,error)
       !
       area%d%refcount = area%d%refcount - 1
       ! If we're releasing the last reference, then free the memory.
       IF (area%d%refcount .EQ. 0) THEN
          IF(.NOT.dbcsr_data_exists(area,error))THEN
              DEALLOCATE (area%d)
          ELSE IF(dbcsr_data_get_size(area)>1 .AND. ASSOCIATED(area%d%memory_type%pool)) THEN
              area%d%ref_size = 0
              CALL dbcsr_mempool_add(area, error)
          ELSE
              CALL internal_data_deallocate(area%d, error)
              DEALLOCATE (area%d)
          ENDIF
          NULLIFY (area%d)
       ENDIF
    ENDIF

    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_data_release

END MODULE dbcsr_data_methods
