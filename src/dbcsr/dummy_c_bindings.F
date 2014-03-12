!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Dummy C bindings
!> \author  Urban Borstnik
!> \date    2011-04-14
!> \version 1.0
!>
!> \par Justification
!>      This module defines alternative types to those found in the
!>      ISO_C_BINDING module.  They have no meaning but exist solely
!>      to enable compilation without peppering the source with
!>      preprocessor commands.
!>
!> <b>Modification history:</b>
!> - Created 2011-04-14
! *****************************************************************************
MODULE dummy_c_bindings

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dummy_c_bindings'

#if defined (__HAS_NO_ISO_C_BINDING)
#if defined (__DBCSR_CUDA)
  ! Compile-time error.  These two flags are incompatible.
  STOP "CUDA needs ISO_C_BINDING.  Compile without __HAS_NO_ISO_C_BINDING"
#endif
  PUBLIC :: c_ptr
  PUBLIC :: c_int
  PUBLIC :: c_size_t

  TYPE c_ptr
     INTEGER :: dummy
  END TYPE c_ptr

  INTEGER, PARAMETER :: c_int    = KIND(0)
  INTEGER, PARAMETER :: c_size_t = KIND(0)
#endif

END MODULE dummy_c_bindings
