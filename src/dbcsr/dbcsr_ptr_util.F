!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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

  USE acc_hostmem,                     ONLY: acc_hostmem_allocate,&
                                             acc_hostmem_deallocate
  USE dbcsr_error_handling
  USE kinds,                     ONLY: real_4,int_4, int_8, real_8, dp
  USE dbcsr_message_passing,           ONLY: mp_allocate,&
                                             mp_deallocate
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_memtype_type,&
                                             dbcsr_memtype_default,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE dbcsr_config,                    ONLY: default_resize_factor

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_ptr_util'

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


  PUBLIC :: ensure_array_size
  PUBLIC :: memory_allocate, memory_deallocate
  PUBLIC :: memory_copy, memory_zero
  PUBLIC :: pointer_replace, pointer_view
  PUBLIC :: pointer_rank_remap2, dbcsr_ptr_remapping,&
            pointer_s_rank_remap2, pointer_d_rank_remap2,&
            pointer_c_rank_remap2, pointer_z_rank_remap2


  INTERFACE ensure_array_size
     MODULE PROCEDURE ensure_array_size_i, ensure_array_size_l
     MODULE PROCEDURE ensure_array_size_s, ensure_array_size_d,&
                      ensure_array_size_c, ensure_array_size_z
  END INTERFACE

  ! Ugly fortran hack
  INTERFACE pointer_view
     MODULE PROCEDURE pointer_view_s, pointer_view_d,&
                      pointer_view_c, pointer_view_z
     MODULE PROCEDURE pointer_view_i, pointer_view_l
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
     MODULE PROCEDURE mem_copy_i, mem_copy_l
     MODULE PROCEDURE mem_copy_s, mem_copy_d, mem_copy_c, mem_copy_z
  END INTERFACE

  INTERFACE memory_zero
     MODULE PROCEDURE mem_zero_i, mem_zero_l
     MODULE PROCEDURE mem_zero_s, mem_zero_d, mem_zero_c, mem_zero_z
  END INTERFACE

  INTERFACE memory_allocate
     MODULE PROCEDURE mem_alloc_i,mem_alloc_l, mem_alloc_s,    mem_alloc_d,    mem_alloc_c,    mem_alloc_z
     MODULE PROCEDURE mem_alloc_i_2d, mem_alloc_s_2d, mem_alloc_d_2d, mem_alloc_c_2d, mem_alloc_z_2d
  END INTERFACE

  INTERFACE memory_deallocate
     MODULE PROCEDURE mem_dealloc_i, mem_dealloc_l, mem_dealloc_s, mem_dealloc_d, mem_dealloc_c, mem_dealloc_z
     MODULE PROCEDURE mem_dealloc_i_2d, mem_dealloc_s_2d, mem_dealloc_d_2d, mem_dealloc_c_2d, mem_dealloc_z_2d
  END INTERFACE

 CONTAINS

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

#include "dbcsr_ptr_util_i.f90"
#include "dbcsr_ptr_util_l.f90"
#include "dbcsr_ptr_util_d.f90"
#include "dbcsr_ptr_util_z.f90"
#include "dbcsr_ptr_util_s.f90"
#include "dbcsr_ptr_util_c.f90"

END MODULE dbcsr_ptr_util
