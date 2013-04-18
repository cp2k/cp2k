!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR pointer and unmanaged array utilities
!> \author  Urban Borstndik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - 2010-02-18 Moved from dbcsr_util
!> - 2010-06-18 Moved all pointer methods into here.
! *****************************************************************************
MODULE dbcsr_ptr_util
#if !defined (__HAS_NO_ISO_C_BINDING)
  USE ISO_C_BINDING
#endif

  USE dbcsr_cuda_memory,               ONLY: dbcsr_cuda_host_mem_alloc,&
                                             dbcsr_cuda_host_mem_dealloc
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: real_4,&
                                             real_8
  USE dbcsr_message_passing,           ONLY: mp_allocate,&
                                             mp_deallocate
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_memory_CUDA_host_pinned,&
                                             dbcsr_memory_MPI,&
                                             dbcsr_memory_default,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_ptr_util'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034

  LOGICAL, PARAMETER :: careful_mod = .FALSE.

#if defined(__PTR_RANK_REMAP) || !defined(__HAS_NO_ISO_C_BINDING)
  ! True pointer rank remapping or safe pointer rank remapping using
  ! ISO_C_BINDING can be used.
  LOGICAL, PARAMETER :: dbcsr_ptr_remapping = .TRUE.
#elif defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
  ! Use buffers
  LOGICAL, PARAMETER :: dbcsr_ptr_remapping = .FALSE.
  ERROR "CP2K must have some way of mapping rank-2 pointers onto rank-1 pointers."
#else
  ! Use crazy Fortran hacks.
  ! This can be very unsafe!
  LOGICAL, PARAMETER :: dbcsr_ptr_remapping = .TRUE.
#endif


  PUBLIC :: ensure_array_size, default_resize_factor
  PUBLIC :: memory_allocate, memory_deallocate
  PUBLIC :: memory_copy, memory_zero
  PUBLIC :: pointer_replace, pointer_view
  PUBLIC :: pointer_rank_remap2, dbcsr_ptr_remapping,&
            pointer_s_rank_remap2, pointer_d_rank_remap2,&
            pointer_c_rank_remap2, pointer_z_rank_remap2


  INTERFACE ensure_array_size
     MODULE PROCEDURE ensure_array_size_i
     MODULE PROCEDURE ensure_array_size_s, ensure_array_size_d,&
                      ensure_array_size_c, ensure_array_size_z
  END INTERFACE

  ! Ugly fortran hack
  INTERFACE pointer_view
     MODULE PROCEDURE pointer_view_s, pointer_view_d,&
                      pointer_view_c, pointer_view_z
     MODULE PROCEDURE pointer_view_i
     MODULE PROCEDURE pointer_view_a
  END INTERFACE

  INTERFACE pointer_replace
     MODULE PROCEDURE pointer_replace_i
  END INTERFACE

  INTERFACE pointer_rank_remap2
     MODULE PROCEDURE pointer_s_rank_remap2, pointer_d_rank_remap2,&
                      pointer_c_rank_remap2, pointer_z_rank_remap2
  END INTERFACE

  INTERFACE memory_copy
     MODULE PROCEDURE mem_copy_i
     MODULE PROCEDURE mem_copy_s, mem_copy_d, mem_copy_c, mem_copy_z
  END INTERFACE

  INTERFACE memory_zero
     MODULE PROCEDURE mem_zero_i
     MODULE PROCEDURE mem_zero_s, mem_zero_d, mem_zero_c, mem_zero_z
  END INTERFACE

  INTERFACE memory_allocate
     MODULE PROCEDURE mem_alloc_i
     MODULE PROCEDURE mem_alloc_s, mem_alloc_d, mem_alloc_c, mem_alloc_z
  END INTERFACE

  INTERFACE memory_deallocate
     MODULE PROCEDURE mem_dealloc_i
     MODULE PROCEDURE mem_dealloc_s, mem_dealloc_d, mem_dealloc_c, mem_dealloc_z
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief Ensures that an array is appropriately large.
!> \param[in,out] array       array to verify and possibly resize
!> \param[in] lb    (optional) desired array lower bound
!> \param[in] ub    desired array upper bound
!> \param[in] factor          (optional) factor by which to exagerrate
!>                            enlargements
!> \param[in] nocopy          (optional) copy array on enlargement; default
!>                            is to copy
!> \param[in] zero_pad        (optional) zero new allocations; default is to
!>                            write nothing
!> \param[in] memory_type     (optional) type of memory
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE ensure_array_size_i(array, lb, ub, factor, nocopy, memory_type,&
       zero_pad, error)
    INTEGER, DIMENSION(:), POINTER           :: array
    INTEGER, INTENT(IN), OPTIONAL            :: lb
    INTEGER, INTENT(IN)                      :: ub
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy
    INTEGER, INTENT(IN), OPTIONAL            :: memory_type
    LOGICAL, INTENT(IN), OPTIONAL            :: zero_pad
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ensure_array_size_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, lb_new, &
                                                lb_orig, mem_type, old_size, &
                                                size_increase, ub_new, ub_orig
    INTEGER, DIMENSION(:), POINTER           :: newarray
    LOGICAL                                  :: dbg, docopy, pad

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handler, error)
    dbg = .FALSE.

    IF (PRESENT (nocopy)) THEN
       docopy = .NOT. nocopy
    ELSE
       docopy = .TRUE.
    ENDIF
    IF (PRESENT (memory_type)) THEN
       mem_type = memory_type
    ELSE
       mem_type = dbcsr_memory_default
    ENDIF
    lb_new = 1
    IF (PRESENT (lb)) lb_new = lb
    pad = .FALSE.
    IF (PRESENT (zero_pad)) pad = zero_pad
    !> Creates a new array if it doesn't yet exist.
    IF (.NOT.ASSOCIATED(array)) THEN
       CALL dbcsr_assert (lb_new, "EQ", 1, &
            dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
            "Arrays must start at 1", __LINE__, error=error)
       CALL mem_alloc_i (array, ub, mem_type=mem_type, error=error)
       IF (pad .AND. ub .GT. 0) CALL mem_zero_i (array, ub)
       IF (careful_mod) CALL dbcsr_error_stop(error_handler, error)
       RETURN
    ENDIF
    lb_orig = LBOUND(array,1)
    ub_orig = UBOUND(array,1)
    old_size = ub_orig - lb_orig + 1
    ! The existing array is big enough.
    IF (lb_orig.LE.lb_new .AND. ub_orig.GE.ub) THEN
       IF (careful_mod) CALL dbcsr_error_stop(error_handler, error)
       RETURN
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Current bounds are',lb_orig,':',ub_orig
    ! A reallocation must be performed.
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
    IF(dbg) WRITE(*,*)routineP//' Resizing to bounds',lb_new,':',ub_new,'v',ub
    !
    ! Deallocates the old array if it's not needed to copy the old data.
    IF(.NOT.docopy) THEN
       CALL mem_dealloc_i (array, mem_type=mem_type,error=error)
    ENDIF
    !
    ! Allocates the new array
    CALL dbcsr_assert (lb_new, "EQ", 1, &
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Arrays must start at 1", __LINE__, error=error)
    CALL mem_alloc_i (newarray, ub_new-lb_new+1, mem_type, error=error)
    !
    ! Now copy and/or zero pad.
    IF(docopy) THEN
       IF(dbg) CALL dbcsr_assert(lb_new.LE.lb_orig .AND. ub_new.GE.ub_orig,&
            dbcsr_failure_level, dbcsr_internal_error, routineP,&
            "Old extent exceeds the new one.",__LINE__,error)
       IF (ub_orig-lb_orig+1 .GT. 0) THEN
          !newarray(lb_orig:ub_orig) = array(lb_orig:ub_orig)
          CALL mem_copy_i (newarray(lb_orig:ub_orig),&
               array(lb_orig:ub_orig), ub_orig-lb_orig+1)
       ENDIF
       IF (pad) THEN
          !newarray(lb_new:lb_orig-1) = 0
          CALL mem_zero_i (newarray(lb_new:lb_orig-1), (lb_orig-1)-lb_new+1)
          !newarray(ub_orig+1:ub_new) = 0
          CALL mem_zero_i (newarray(ub_orig+1:ub_new), ub_new-(ub_orig+1)+1)
       ENDIF
       CALL mem_dealloc_i (array, mem_type, error=error)
    ELSEIF (pad) THEN
       !newarray(:) = 0
       CALL mem_zero_i (newarray, SIZE(newarray))
    END IF
    array => newarray
    IF (careful_mod) CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE ensure_array_size_i


! *****************************************************************************
!> \brief Copies memory area
!> \param[out] dst   destination memory
!> \param[in] src    source memory
!> \param[in] n      length of copy
! *****************************************************************************
  SUBROUTINE mem_copy_i (dst, src, n)
    INTEGER, INTENT(IN)                      :: n
    INTEGER, DIMENSION(1:n), INTENT(IN)      :: src
    INTEGER, DIMENSION(1:n), INTENT(OUT)     :: dst

    dst(:) = src(:)
  END SUBROUTINE mem_copy_i

! *****************************************************************************
!> \brief Zeros memory area
!> \param[out] dst   destination memory
!> \param[in] n      length of elements to zero
! *****************************************************************************
  SUBROUTINE mem_zero_i (dst, n)
    INTEGER, INTENT(IN)                      :: n
    INTEGER, DIMENSION(1:n), INTENT(OUT)     :: dst

    dst(:) = 0
  END SUBROUTINE mem_zero_i

! *****************************************************************************
!> \brief Allocates memory
!> \param[out] mem        memory to allocate
!> \param[in] n           length of elements to allocate
!> \param[in] mem_type    memory type
!> \param[in,out] error   error
! *****************************************************************************
  SUBROUTINE mem_alloc_i (mem, n, mem_type, error)
    INTEGER, DIMENSION(:), POINTER           :: mem
    INTEGER, INTENT(IN)                      :: n, mem_type
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mem_alloc_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, stat

!   ---------------------------------------------------------------------------

    IF (careful_mod) &
       CALL dbcsr_error_set (routineN, error_handle, error=error)
    !
    SELECT CASE (mem_type)
    CASE (dbcsr_memory_default)
       ALLOCATE(mem(n), stat=stat)
    CASE (dbcsr_memory_MPI)
       CALL mp_allocate(mem, n, stat=stat)
    CASE (dbcsr_memory_CUDA_host_pinned)
       CALL dbcsr_cuda_host_mem_alloc(mem, n, error=error)
    CASE default
       CALL dbcsr_assert (.FALSE.,&
            dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
            "Unsupported memory type",__LINE__,error)
    END SELECT
    !
    CALL dbcsr_assert (stat == 0, dbcsr_warning_level, dbcsr_internal_error,&
         routineN, "memory",__LINE__,error)
    !
    IF (careful_mod) &
       CALL dbcsr_error_stop (error_handle, error=error)
  END SUBROUTINE mem_alloc_i

! *****************************************************************************
!> \brief Deallocates memory
!> \param[out] mem        memory to allocate
!> \param[in] n           length of elements to allocate
!> \param[in] mem_type    memory type
!> \param[in,out] error   error
! *****************************************************************************
  SUBROUTINE mem_dealloc_i (mem, mem_type, error)
    INTEGER, DIMENSION(:), POINTER           :: mem
    INTEGER, INTENT(IN)                      :: mem_type
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mem_dealloc_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, stat

!   ---------------------------------------------------------------------------

    IF (careful_mod) &
       CALL dbcsr_error_set (routineN, error_handle, error=error)
    !
    SELECT CASE (mem_type)
    CASE (dbcsr_memory_default)
       DEALLOCATE(mem, stat=stat)
    CASE (dbcsr_memory_MPI)
       CALL mp_deallocate(mem, stat=stat)
    CASE (dbcsr_memory_CUDA_host_pinned)
       CALL dbcsr_cuda_host_mem_dealloc(mem, error=error)
    CASE default
       CALL dbcsr_assert (.FALSE.,&
               dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
               "Unsupported memory type",__LINE__,error)
    END SELECT
    !
    CALL dbcsr_assert (stat == 0, dbcsr_warning_level, dbcsr_internal_error,&
         routineN, "memory",__LINE__,error)
    !
    IF (careful_mod) &
       CALL dbcsr_error_stop (error_handle, error=error)
  END SUBROUTINE mem_dealloc_i

! *****************************************************************************
!> \brief Replaces an existing pointer with a new one, freeing memory as
!>        required.
!> \param[in,out] original_p  original pointer, to be replaced
!> \param[in] new_p           replacement pointer
! *****************************************************************************
  SUBROUTINE pointer_replace_i (original_p, new_p)
    INTEGER, DIMENSION(:), POINTER           :: original_p, new_p

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (original_p)) DEALLOCATE (original_p)
    original_p => new_p
  END SUBROUTINE pointer_replace_i


! *****************************************************************************
!> \brief Returns a pointer with different bounds.
!> \param[in] original   original data pointer
!> \param[in] lb, ub     lower and upper bound for the new pointer view
!> \param[out] view      new pointer
! *****************************************************************************
  FUNCTION pointer_view_i (original, lb, ub) RESULT (view)
    INTEGER, DIMENSION(:), POINTER           :: original
    INTEGER, INTENT(IN)                      :: lb, ub
    INTEGER, DIMENSION(:), POINTER           :: view

    view => original(lb:ub)
  END FUNCTION pointer_view_i


! *****************************************************************************
!> \brief Repoints a pointer into a part of a data area
!> \param[in,out] new_area    repoints this encapsulated pointer
!> \param[in] area            area to point into
!> \param[in] offset          point to this offset in area
!> \param[in] len             (optional) length of data area to point to
!> \result narea2             copy of new_area
! *****************************************************************************
  FUNCTION pointer_view_a (new_area, area, offset, len) RESULT (narea2)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: new_area
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER, INTENT(IN)                      :: offset
    INTEGER, INTENT(IN), OPTIONAL            :: len
    TYPE(dbcsr_data_obj)                     :: narea2

    CHARACTER(len=*), PARAMETER :: routineN = 'pointer_view_a', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

    CALL dbcsr_assert (area%d%data_type, "EQ", new_area%d%data_type,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Incompatible data types.",__LINE__,error)
    IF (PRESENT (len)) THEN
       SELECT CASE (area%d%data_type)
       CASE (dbcsr_type_real_4)
          new_area%d%r_sp => area%d%r_sp(offset:offset+len-1)
       CASE (dbcsr_type_real_8)
          new_area%d%r_dp => area%d%r_dp(offset:offset+len-1)
       CASE (dbcsr_type_complex_4)
          new_area%d%c_sp => area%d%c_sp(offset:offset+len-1)
       CASE (dbcsr_type_complex_8)
          new_area%d%c_dp => area%d%c_dp(offset:offset+len-1)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
               routineN, "Invalid data type.",__LINE__,error)
       END SELECT
    ELSE
       SELECT CASE (area%d%data_type)
       CASE (dbcsr_type_real_4)
          new_area%d%r_sp => area%d%r_sp(offset:)
       CASE (dbcsr_type_real_8)
          new_area%d%r_dp => area%d%r_dp(offset:)
       CASE (dbcsr_type_complex_4)
          new_area%d%c_sp => area%d%c_sp(offset:)
       CASE (dbcsr_type_complex_8)
          new_area%d%c_dp => area%d%c_dp(offset:)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
               routineN, "Invalid data type.",__LINE__,error)
       END SELECT
    ENDIF
    narea2 = new_area
  END FUNCTION pointer_view_a

#include "dbcsr_ptr_util_d.F"
#include "dbcsr_ptr_util_z.F"
#include "dbcsr_ptr_util_s.F"
#include "dbcsr_ptr_util_c.F"

END MODULE dbcsr_ptr_util
