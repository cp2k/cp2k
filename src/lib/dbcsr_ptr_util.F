!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
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

  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: real_4,&
                                             real_8
  USE dbcsr_message_passing,           ONLY: mp_allocate,&
                                             mp_deallocate
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_ptr_util'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034

#if defined(__PTR_RANK_REMAP)
  ! True pointer rank remapping can be used.
  LOGICAL, PARAMETER :: dbcsr_ptr_remapping = .TRUE.
#else
#if defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
  ! Use buffers
  LOGICAL, PARAMETER :: dbcsr_ptr_remapping = .FALSE.
#else
  ! Use crazy Fortran hacks.
  ! This can be very unsafe!
  LOGICAL, PARAMETER :: dbcsr_ptr_remapping = .TRUE.
#endif
#endif


  PUBLIC :: addto_array, ensure_array_size
  PUBLIC :: pointer_replace, pointer_view
  PUBLIC :: pointer_rank_remap2, dbcsr_ptr_remapping,&
            pointer_s_rank_remap2, pointer_d_rank_remap2,&
            pointer_c_rank_remap2, pointer_z_rank_remap2


  INTERFACE addto_array
     MODULE PROCEDURE addto_array_i
     MODULE PROCEDURE addto_array_s, addto_array_d,&
                      addto_array_c, addto_array_z
     MODULE PROCEDURE addto_array_as, addto_array_ad,&
                      addto_array_ac, addto_array_az
  END INTERFACE

  INTERFACE ensure_array_size
     MODULE PROCEDURE ensure_array_size_i,&! ensure_array_size_i_2d,&
          ensure_array_size_s, ensure_array_size_d, ensure_array_size_c,&
          ensure_array_size_z
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

CONTAINS

! *****************************************************************************
!> \brief Copies an array into another array, ensuring that the target array
!>        is large enough.
!> \par Makes use of FORTRAN's array bounds for array sections.
!> \param to        array to copy to (target array)
!> \param[in] from  array to copy from (source array)
!> \param[in] offset          starting position in the target array
!> \param[in] factor          factor for size increase
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE addto_array_i(to, from, offset, factor, special, error)
    INTEGER, DIMENSION(:), POINTER           :: to
    INTEGER, DIMENSION(:), INTENT(IN)        :: from
    INTEGER, INTENT(IN), OPTIONAL            :: offset
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: special
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'addto_array_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, last, lb_t, &
                                                s_f, s_t, start, ub_t
    LOGICAL                                  :: dbg, spec
    REAL                                     :: f

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    dbg = .FALSE.
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    s_f = SIZE(from, 1)
    ub_t = UBOUND(to, 1)
    lb_t = LBOUND(to, 1)
    s_t = SIZE(to)
    start = lb_t
    IF (PRESENT (offset)) start = offset
    ! assert offset > lb_t
    last = start+s_f-1
    IF (ub_t .LT. last) THEN
       IF (PRESENT (factor)) THEN
          f = factor
       ELSE
          f = 1.4142136
       ENDIF
       ub_t = MAX(start + INT(s_t * f), last)
       IF(dbg) WRITE(*,*)routineP//' expanding array from :',UBOUND(to),&
            'to',ub_t
       CALL ensure_array_size(to, lb_t, ub_t, factor=1.0, nocopy=.FALSE.,&
            special=spec, error=error)
    ENDIF
    to(start:last) = from(:)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE addto_array_i


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
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE ensure_array_size_i(array, lb, ub, factor, nocopy, special,&
       zero_pad, error)
    INTEGER, DIMENSION(:), POINTER           :: array
    INTEGER, INTENT(IN), OPTIONAL            :: lb
    INTEGER, INTENT(IN)                      :: ub
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, special, zero_pad
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ensure_array_size_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, lb_new, &
                                                lb_orig, old_size, &
                                                size_increase, stat, ub_new, &
                                                ub_orig
    INTEGER, DIMENSION(:), POINTER           :: newarray
    LOGICAL                                  :: dbg, docopy = .TRUE., pad, &
                                                spec

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    dbg = .FALSE.

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
       !CPPrecondition(.NOT.docopy, dbcsr_warning_level, routineP, error, failure)
       IF (spec) THEN
          CALL dbcsr_assert (lb_new,'EQ',1, dbcsr_warning_level,&
               dbcsr_internal_error,&
               routineN, "Special memory can only start at 1.",__LINE__,error)
          CALL mp_allocate(array, ub-lb_new+1, stat=stat)
       ELSE
          ALLOCATE(array(lb_new:ub), stat=stat)
       ENDIF
       CALL dbcsr_assert (stat == 0, dbcsr_warning_level, dbcsr_internal_error,&
            routineN, "array",__LINE__,error)
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
    IF(.NOT.docopy) THEN
       IF (spec) THEN
          CALL mp_deallocate(array)
       ELSE
          DEALLOCATE(array)
       ENDIF
    ENDIF
    IF (spec) THEN
       CALL dbcsr_assert (lb_new,'EQ',1, dbcsr_warning_level, dbcsr_internal_error,&
            routineN, "Special memory can only start at 1.",__LINE__,error)
       CALL mp_allocate(newarray, ub_new-lb_new+1, stat=stat)
    ELSE
       ALLOCATE(newarray(lb_new:ub_new), stat=stat)
    ENDIF
    CALL dbcsr_assert (stat == 0, dbcsr_warning_level, dbcsr_internal_error,&
         routineN, "newarray",__LINE__,error)
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
          array(ub_orig+1:ub_new) = 0
          array(lb_new:lb_orig-1) = 0
       ENDIF
    ELSEIF (pad) THEN
       newarray(:) = 0
    END IF
    array => newarray
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE ensure_array_size_i


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
