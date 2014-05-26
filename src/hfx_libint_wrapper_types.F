!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Provides types that are used to interface with the libint library
!> \par History
!>      10.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE hfx_libint_wrapper_types

  USE, INTRINSIC :: ISO_C_BINDING
  USE kinds,                           ONLY: dp, int_8

  IMPLICIT NONE
  PRIVATE
  PUBLIC lib_int,&
         lib_deriv,&
         prim_data,&
         prim_data_f_size,&
         build_eri_size,&
         build_deriv1_eri_size

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_libint_wrapper_types'

! *** Set default libint parameters if not specified via arch file

#if !defined (__LIBINT_MAX_AM)
#define __LIBINT_MAX_AM 5
#endif
#if !defined (__LIBDERIV_MAX_AM1)
#define __LIBDERIV_MAX_AM1 4
#endif

  INTEGER, PARAMETER :: libint_max_am = __LIBINT_MAX_AM
  INTEGER, PARAMETER :: libderiv_max_am1 = __LIBDERIV_MAX_AM1
  INTEGER, PARAMETER :: prim_data_f_size = 4 * (libint_max_am - 1) + 1
  INTEGER, PARAMETER :: libint_vrr_classes_size = 2 * (libint_max_am - 1) + 1
  INTEGER, PARAMETER :: libint_dvrr_classes_size = 2 * (libderiv_max_am1 - 1) + 1
  INTEGER, PARAMETER :: build_eri_size = libint_max_am - 1
  INTEGER, PARAMETER :: build_deriv1_eri_size = libderiv_max_am1 - 1

  TYPE, BIND(C) :: prim_data
    REAL(C_DOUBLE) :: F(prim_data_f_size)
    REAL(C_DOUBLE) :: U(3,6)
    REAL(C_DOUBLE) :: twozeta_a
    REAL(C_DOUBLE) :: twozeta_b
    REAL(C_DOUBLE) :: twozeta_c
    REAL(C_DOUBLE) :: twozeta_d
    REAL(C_DOUBLE) :: oo2z
    REAL(C_DOUBLE) :: oo2n
    REAL(C_DOUBLE) :: oo2zn
    REAL(C_DOUBLE) :: poz
    REAL(C_DOUBLE) :: pon
    REAL(C_DOUBLE) :: oo2p
    REAL(C_DOUBLE) :: ss_r12_ss
  END TYPE prim_data

  TYPE, BIND(C) :: lib_int
    TYPE(C_PTR)     :: int_stack
    TYPE(C_PTR)     :: PrimQuartet
    REAL(C_DOUBLE)  :: AB(3)
    REAL(C_DOUBLE)  :: CD(3)
    TYPE(C_PTR)     :: vrr_classes(libint_vrr_classes_size,libint_vrr_classes_size)
    TYPE(C_PTR)     :: vrr_stack
  END TYPE lib_int

  TYPE, BIND(C) :: lib_deriv
    TYPE(C_PTR)     :: int_stack
    TYPE(C_PTR)     :: PrimQuartet
    TYPE(C_PTR)     :: zero_stack
    TYPE(C_PTR)     :: ABCD(156)
    REAL(C_DOUBLE)  :: AB(3)
    REAL(C_DOUBLE)  :: CD(3)
    TYPE(C_PTR)     :: deriv_classes(12,libint_dvrr_classes_size,libint_dvrr_classes_size)
    TYPE(C_PTR)     :: deriv2_classes(144,libint_dvrr_classes_size,libint_dvrr_classes_size)
    TYPE(C_PTR)     :: dvrr_classes(libint_dvrr_classes_size,libint_dvrr_classes_size)
    TYPE(C_PTR)     :: dvtt_stack
  END TYPE lib_deriv

END MODULE hfx_libint_wrapper_types
