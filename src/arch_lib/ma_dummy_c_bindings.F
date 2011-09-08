!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  Urban Borstnik and the CP2K developers group   !
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
MODULE ma_dummy_c_bindings

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ma_dummy_c_bindings'

#if defined (__HAS_NO_ISO_C_BINDING)
#if defined (__HWLOC) || defined(__LIBNUMA)
  ! Compile-time error.  These two flags are incompatible.
  STOP "Arch Lib needs ISO_C_BINDING.  Compile without __HAS_NO_ISO_C_BINDING"
#endif
#endif

#if defined (__HWLOC) && defined(__LIBNUMA)
  ! Compile-time error.  These two flags are incompatible.
  STOP "HWLOC and LIBNUMA libraries are incompatible! Try without one of them."
#endif

#if defined (__HAS_NO_ISO_C_BINDING)
#if !defined (__HWLOC) || !defined (__LIBNUMA)
  PUBLIC :: c_int 
  PUBLIC :: c_char

  TYPE c_char
   CHARACTER(len=1) c
  END TYPE c_char

  INTEGER, PARAMETER :: c_int    = KIND(0)
#endif
#endif

END MODULE ma_dummy_c_bindings
