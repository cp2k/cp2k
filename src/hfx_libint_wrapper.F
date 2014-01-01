!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to the Libint-Library or a c++ wrapper.
!> \par History
!>      11.2007 created [Manuel Guidon]
!>      10.2009 refactored [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE hfx_libint_wrapper

#if defined (__LIBINT)

#if !defined (__HAS_NO_ISO_C_BINDING)
#if !defined (__HAS_ISO_C_BINDING)
#define __HAS_ISO_C_BINDING
#endif
#endif

#if defined (__HAS_ISO_C_BINDING)
  USE, INTRINSIC :: ISO_C_BINDING
#endif
#endif

  USE kinds,                           ONLY: dp
  USE mathconstants
  USE hfx_libint_wrapper_types
  USE orbital_pointers

  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC initialize_libint, terminate_libint,&
         initialize_libderiv, &
         get_eris, get_derivs, terminate_libderiv

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_libint_wrapper'

#if defined (__LIBINT)

!****************************************************************************!
!****************************************************************************!
!***                                                                      ***!
!***  WHAT FOLLOWS IS CODE FOR COMPILERS THAT EITHER FULLY SUPPORT THE    ***!
!***  ISO_C_BINDING FORTRAN 2003 STANDARD OR NOT. IF NOT, ALL CALLS TO    ***!
!***  LIBINT ARE REDIRECTED TO A C++ WRAPPER                              ***!
!***                                                                      ***!
!****************************************************************************!
!****************************************************************************!

#if defined (__HAS_ISO_C_BINDING)
  TYPE(C_FUNPTR), DIMENSION(0:build_eri_size,0:build_eri_size,0:build_eri_size,0:build_eri_size), BIND(C) :: build_eri
  TYPE(C_FUNPTR), DIMENSION(0:build_deriv1_eri_size,0:build_deriv1_eri_size,&
                            0:build_deriv1_eri_size,0:build_deriv1_eri_size), BIND(C) :: build_deriv1_eri

  INTERFACE
    FUNCTION build(lib, np) BIND(C)
      USE hfx_libint_wrapper_types
      USE, INTRINSIC                 :: ISO_C_BINDING
      TYPE(C_PTR)                    :: build
      TYPE(lib_int)                  :: lib
      INTEGER(KIND=C_INT), VALUE     :: np
    END FUNCTION build

    FUNCTION init_lib(lib, max_am, np) BIND(C, name="init_libint")
      USE hfx_libint_wrapper_types
      USE, INTRINSIC                 :: ISO_C_BINDING
      INTEGER(KIND=C_INT)            :: init_lib
      TYPE(lib_int)                  :: lib
      INTEGER(KIND=C_INT), VALUE     :: max_am
      INTEGER(KIND=C_INT), VALUE     :: np
    END FUNCTION init_lib

    SUBROUTINE init_base() BIND(C, name="init_libint_base")
    END SUBROUTINE init_base

    SUBROUTINE free_lib_int(lib) BIND(C, name="free_libint")
      USE hfx_libint_wrapper_types
      USE, INTRINSIC                 :: ISO_C_BINDING
      TYPE(lib_int)                  :: lib
    END SUBROUTINE free_lib_int

    SUBROUTINE init_deriv_base() BIND(C, name="init_libderiv_base")
    END SUBROUTINE init_deriv_base

    FUNCTION init_deriv1(deriv, max_am, np, ccs) BIND(C, name="init_libderiv1")
      USE hfx_libint_wrapper_types
      USE, INTRINSIC                 :: ISO_C_BINDING
      INTEGER(KIND=C_INT)            :: init_deriv1
      TYPE(lib_deriv)                :: deriv
      INTEGER(KIND=C_INT), VALUE     :: max_am
      INTEGER(KIND=C_INT), VALUE     :: np
      INTEGER(KIND=C_INT), VALUE     :: ccs
    END FUNCTION init_deriv1

    SUBROUTINE build_deriv1(deriv, np) BIND(C)
      USE hfx_libint_wrapper_types
      USE, INTRINSIC                 :: ISO_C_BINDING
      TYPE(lib_deriv)                :: deriv
      INTEGER(KIND=C_INT),VALUE      :: np
    END SUBROUTINE build_deriv1

    SUBROUTINE free_lib_deriv(deriv) BIND(C, name="free_libderiv")
      USE hfx_libint_wrapper_types
      USE, INTRINSIC                 :: ISO_C_BINDING
      TYPE(lib_deriv)                  :: deriv
    END SUBROUTINE free_lib_deriv
  END INTERFACE
#endif

  CONTAINS

  SUBROUTINE initialize_libint(lib,max_am,error)
    TYPE(lib_int)                            :: lib
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libint', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
#if defined (__HAS_ISO_C_BINDING)
    INTEGER(KIND=C_INT)                      :: lib_storage, max_am_local, &
                                                max_prim
#else
    INTEGER                                  :: lib_storage, max_am_local, &
                                                max_prim
    EXTERNAL wrapper_init_lib
#endif

    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
#if defined (__HAS_ISO_C_BINDING)
    CALL init_base()
    lib_storage = init_lib(lib, max_am_local, max_prim)
#else
    CALL wrapper_init_lib(lib, max_am_local, max_prim, lib_storage)
#endif

    IF (lib_storage<0) THEN
      CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                     " the angular momentum needed exceeds the value assumed when configuring libint ", &
                     error,failure)
    ENDIF
  END SUBROUTINE initialize_libint

  SUBROUTINE initialize_libderiv(deriv,max_am,error)
    TYPE(lib_deriv)                          :: deriv
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libderiv', &
      routineP = moduleN//':'//routineN
    LOGICAL                                  :: failure
#if defined (__HAS_ISO_C_BINDING)
    INTEGER(KIND=C_INT)                      :: lib_deriv_storage, &
                                                max_am_local, max_classes, &
                                                max_prim
#else
    INTEGER                                  :: lib_deriv_storage, &
                                                max_am_local, max_classes, &
                                                max_prim
    EXTERNAL wrapper_init_deriv
#endif

    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    max_classes = nco(max_am)**4

#if defined (__HAS_ISO_C_BINDING)
    CALL init_deriv_base()
    lib_deriv_storage = init_deriv1(deriv, max_am_local, max_prim, max_classes)
#else
    CALL wrapper_init_deriv(deriv, max_am_local, max_prim, max_classes, lib_deriv_storage)
#endif
    IF (lib_deriv_storage<0) THEN
      CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                     " the angular momentum needed exceeds the value assumed when configuring libderiv ", &
                     error,failure)
    ENDIF
  END SUBROUTINE initialize_libderiv

  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                            :: lib
#if defined (__HAS_ISO_C_BINDING)
      CALL free_lib_int(lib)
#else
    EXTERNAL wrapper_free_libint
      CALL wrapper_free_libint(lib)
#endif
  END SUBROUTINE terminate_libint

  SUBROUTINE terminate_libderiv(deriv)
    TYPE(lib_deriv)                          :: deriv

#if defined (__HAS_ISO_C_BINDING)
      CALL free_lib_deriv(deriv)
#else
    EXTERNAL wrapper_free_libderiv
      CALL wrapper_free_libderiv(deriv)
#endif
  END SUBROUTINE terminate_libderiv

  SUBROUTINE get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_int)                            :: lib
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(:), POINTER          :: p_work
    INTEGER                                  :: a_mysize(1)

#if defined (__HAS_ISO_C_BINDING)
    PROCEDURE(build), POINTER               :: pbuild
    TYPE(C_PTR)                             :: pc_result
    REAL(dp), DIMENSION(:), POINTER         :: p_tmp
#else
    EXTERNAL wrapper_build_eri
#endif

#if defined (__HAS_ISO_C_BINDING)
    lib%PrimQuartet = C_LOC(prim)
    CALL C_F_PROCPOINTER(build_eri(n_d,n_c,n_b,n_a),pbuild)
    pc_result = pbuild(lib,1)
    CALL C_F_POINTER(pc_result, p_tmp, a_mysize)
    p_work => p_tmp
#else
    CALL wrapper_build_eri(n_a, n_b, n_c, n_d, lib, a_mysize(1), p_work(1), prim)
#endif
  END SUBROUTINE get_eris

  SUBROUTINE get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_deriv)                          :: deriv
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*&
      nco(n_c)*nco(n_d), 12)                 :: work_forces
    INTEGER                                  :: a_mysize(1)

#if defined (__HAS_ISO_C_BINDING)
    PROCEDURE(build_deriv1), POINTER         :: pbuild_deriv1
    TYPE(C_PTR)                              :: pc_result
    REAL(C_DOUBLE), DIMENSION(:), POINTER    :: tmp_data
#else
    EXTERNAL wrapper_build_deriv1_eri
#endif
    INTEGER                                  :: i, k

#if defined (__HAS_ISO_C_BINDING)
    deriv%PrimQuartet = C_LOC(prim)
    CALL C_F_PROCPOINTER(build_deriv1_eri(n_d,n_c,n_b,n_a),pbuild_deriv1)
    CALL pbuild_deriv1(deriv,1)

    DO k=1,12
      IF(k==4 .OR. k==5 .OR. k==6) CYCLE
      pc_result = deriv%ABCD(k)
      CALL C_F_POINTER(pc_result, tmp_data , a_mysize)
      DO i=1,a_mysize(1)
        work_forces(i,k) = tmp_data(i)
      ENDDO
    END DO
#else
    CALL wrapper_build_deriv1_eri(n_a, n_b, n_c, n_d, deriv, a_mysize(1), work_forces(1,1),prim)
#endif
  END SUBROUTINE get_derivs

#else

!****************************************************************************!
!****************************************************************************!
!***                                                                      ***!
!***  WHAT FOLLOWS IS CODE THAT USES BOGUS SUBROUTINES AND TYPES IN       ***!
!***  ORDER TO ALLOW ALL COMPILERS TO COMPILE CP2K                        ***!
!***                                                                      ***!
!****************************************************************************!
!****************************************************************************!

  CONTAINS

  SUBROUTINE initialize_libint(lib,max_am,error)
    TYPE(lib_int)                            :: lib
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libint', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                   " This CP2K executable has not been linked against the library libint, required for HFX.", &
                   error,failure)

  END SUBROUTINE initialize_libint

  SUBROUTINE initialize_libderiv(deriv,max_am,error)
    TYPE(lib_deriv)                          :: deriv
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libderiv', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                   " This CP2K executable has not been linked against the library libint, required for HFX.", &
                   error,failure)

  END SUBROUTINE initialize_libderiv

  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                            :: lib

  END SUBROUTINE terminate_libint

  SUBROUTINE terminate_libderiv(deriv)
    TYPE(lib_deriv)                          :: deriv

  END SUBROUTINE terminate_libderiv

  SUBROUTINE get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_int)                            :: lib
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(:), POINTER          :: p_work
    INTEGER                                  :: a_mysize(1)

  END SUBROUTINE get_eris

  SUBROUTINE get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_deriv)                          :: deriv
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*&
      nco(n_c)*nco(n_d), 12)                 :: work_forces
    INTEGER                                  :: a_mysize(1)

  END SUBROUTINE get_derivs
#endif

END MODULE hfx_libint_wrapper
