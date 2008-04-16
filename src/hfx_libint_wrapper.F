!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to the Libint-Library or a c++ wrapper. 
!> \note
!>      IMPORTANT NOTE : this file currently is for a libint configured for
!>                       LIBINT_MAX_AM 5 and LIBINT_MAX_AM1 4
!> \par History
!>      11.2007 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE hfx_libint_wrapper

#if defined (__LIBINT)
#if defined (__G95)
  USE, INTRINSIC :: ISO_C_BINDING
#elif defined (__INTEL) || defined (__GFORTRAN)
  USE, INTRINSIC :: ISO_C_BINDING
#elif defined (__NAG)
  USE, INTRINSIC :: ISO_C_BINDING
#else
  !! We do not need it
!COMPILERS
#endif
!LIBINT
#endif 
  USE kinds,                           ONLY: dp, int_8
  USE mathconstants
  USE orbital_pointers

  USE f77_blas
#include "cp_common_uses.h"
 
  IMPLICIT NONE
  PRIVATE
  PUBLIC initialize_libint, lib_int, terminate_libint,&
         initialize_libderiv, lib_deriv, prim_data,&
         get_eris, get_derivs

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_libint_wrapper'

  INTEGER, PARAMETER :: libint_max_am = 5
  INTEGER, PARAMETER :: build_eri_size = libint_max_am-1
  INTEGER, PARAMETER :: prim_data_f_size = 17
  INTEGER, PARAMETER :: lib_int_vrr_classes_size = 9
  INTEGER, PARAMETER :: libderiv_max_am1 = 4
  INTEGER, PARAMETER :: build_deriv1_eri_size = libderiv_max_am1 - 1 

#if defined (__LIBINT)
!****************************************************************************!
!****************************************************************************!
!***                                                                      ***!
!***  WHAT FOLLOWS IS CODE FOR G95 COMPILER THAT FULLY SUPPORTS THE       ***!
!***  ISO_C_BINDING FORTRAN 2003 STANDARD                                 ***!
!***                                                                      ***!
!****************************************************************************!
!****************************************************************************!

#if defined (__G95)
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
    TYPE(C_PTR)     :: vrr_classes(lib_int_vrr_classes_size,lib_int_vrr_classes_size)
    TYPE(C_PTR)     :: vrr_stack
  END TYPE lib_int

  TYPE, BIND(C) :: lib_deriv
    TYPE(C_PTR)     :: int_stack
    TYPE(C_PTR)     :: PrimQuartet
    TYPE(C_PTR)     :: zero_stack
    TYPE(C_PTR)     :: ABCD(156)
    REAL(C_DOUBLE)  :: AB(3)
    REAL(C_DOUBLE)  :: CD(3)
    TYPE(C_PTR)     :: deriv_classes(12,7,7)
    TYPE(C_PTR)     :: deriv2_classes(144,7,7)
    TYPE(C_PTR)     :: dvrr_classes(7,7)
    TYPE(C_PTR)     :: dvtt_stack
  END TYPE lib_deriv

  !global vars
  TYPE(C_FUNPTR), DIMENSION(0:build_eri_size,0:build_eri_size,0:build_eri_size,0:build_eri_size), BIND(C) :: build_eri
  TYPE(C_FUNPTR), DIMENSION(0:build_deriv1_eri_size,0:build_deriv1_eri_size,&
                            0:build_deriv1_eri_size,0:build_deriv1_eri_size), BIND(C) :: build_deriv1_eri

  INTERFACE
! *****************************************************************************
    FUNCTION build(lib, np) BIND(C)
      IMPORT
      TYPE(C_PTR)                    :: build
      TYPE(lib_int)                  :: lib
      INTEGER(KIND=C_INT), VALUE     :: np
    END FUNCTION build
  
! *****************************************************************************
    FUNCTION init_lib(lib, max_am, np) BIND(C, name="init_libint")
      IMPORT
      INTEGER(KIND=C_INT)            :: init_lib
      TYPE(lib_int)                  :: lib
      INTEGER(KIND=C_INT), VALUE     :: max_am
      INTEGER(KIND=C_INT), VALUE     :: np
    END FUNCTION init_lib
    
! *****************************************************************************
    SUBROUTINE init_base() BIND(C, name="init_libint_base")
    END SUBROUTINE init_base

! *****************************************************************************
    SUBROUTINE free_lib_int(lib) BIND(C, name="free_libint")
      IMPORT
      TYPE(lib_int)                  :: lib
    END SUBROUTINE free_lib_int

! *****************************************************************************
    SUBROUTINE init_deriv_base() BIND(C, name="init_libderiv_base")
    END SUBROUTINE init_deriv_base

! *****************************************************************************
    FUNCTION init_deriv1(deriv, max_am, np, ccs) BIND(C, name="init_libderiv1")
      IMPORT
      INTEGER(KIND=C_INT)            :: init_deriv1
      TYPE(lib_deriv)                :: deriv
      INTEGER(KIND=C_INT), VALUE     :: max_am
      INTEGER(KIND=C_INT), VALUE     :: np
      INTEGER(KIND=C_INT), VALUE     :: ccs
    END FUNCTION init_deriv1

! *****************************************************************************
    SUBROUTINE build_deriv1(deriv, np) BIND(C)
      IMPORT
      TYPE(lib_deriv)                :: deriv
      INTEGER(KIND=C_INT),VALUE      :: np
    END SUBROUTINE build_deriv1

! *****************************************************************************
    SUBROUTINE free_lib_deriv(deriv) BIND(C, name="free_libderiv")
      IMPORT
      TYPE(lib_deriv)                  :: deriv
    END SUBROUTINE free_lib_deriv

  END INTERFACE

!***

  CONTAINS

! *****************************************************************************
  SUBROUTINE initialize_libint(lib,max_am,error)
    TYPE(lib_int)                            :: lib
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libint', &
      routineP = moduleN//':'//routineN

    INTEGER(KIND=C_INT)                      :: lib_storage, max_am_local, &
                                                max_prim
    LOGICAL                                  :: failure

    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    CALL init_base()
    lib_storage = init_lib(lib, max_am_local, max_prim)
    
    IF (lib_storage<0) THEN
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " the angular momentum needed exceeds the value assumed when configuring libint ", &
                        error,failure)
    ENDIF
  END SUBROUTINE initialize_libint

! *****************************************************************************
  SUBROUTINE initialize_libderiv(deriv,max_am,error)
    TYPE(lib_deriv)                          :: deriv
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libderiv', &
      routineP = moduleN//':'//routineN

    INTEGER(KIND=C_INT)                      :: lib_deriv_storage, &
                                                max_am_local, max_classes, &
                                                max_prim
    LOGICAL                                  :: failure

    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    max_classes = 10000!max_am**4

    CALL init_deriv_base()
    lib_deriv_storage = init_deriv1(deriv, max_am_local, max_prim, max_classes)

    IF (lib_deriv_storage<0) THEN
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " the angular momentum needed exceeds the value assumed when configuring libderiv ", &
                        error,failure)
    ENDIF
  END SUBROUTINE initialize_libderiv

! *****************************************************************************
  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                            :: lib

      CALL free_lib_int(lib)
  END SUBROUTINE terminate_libint

! *****************************************************************************
  SUBROUTINE terminate_libderiv(deriv)
    TYPE(lib_deriv)                          :: deriv

      CALL free_lib_deriv(deriv)
  END SUBROUTINE terminate_libderiv
 
! *****************************************************************************
  SUBROUTINE get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_int)                            :: lib
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(:), POINTER          :: p_work
    INTEGER                                  :: a_mysize(1)

    PROCEDURE(build), POINTER               :: pbuild
    TYPE(C_PTR)                             :: pc_result 
    REAL(dp), DIMENSION(:), POINTER         :: p_tmp
  
    lib%PrimQuartet = C_LOC(prim) 
    CALL C_F_PROCPOINTER(build_eri(n_d,n_c,n_b,n_a),pbuild)
    pc_result = pbuild(lib,1)
    CALL C_F_POINTER(pc_result, p_tmp, a_mysize)
    p_work => p_tmp
  END SUBROUTINE get_eris 

! *****************************************************************************
  SUBROUTINE get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_deriv)                          :: deriv
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*&
      nco(n_c)*nco(n_d), 12)                 :: work_forces
    INTEGER                                  :: a_mysize(1)

    PROCEDURE(build_deriv1), POINTER        :: pbuild_deriv1
    TYPE(C_PTR)                             :: pc_result 
    INTEGER                                 :: i, k
    REAL(C_DOUBLE), DIMENSION(:), POINTER   :: tmp_data

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
  END SUBROUTINE get_derivs 

!****************************************************************************!
!****************************************************************************!
!***                                                                      ***!
!***  WHAT FOLLOWS IS CODE FOR INTEL COMPILERS THAT PARTIALLY SUPPORT THE ***!
!***  ISO_C_BINDING FORTRAN 2003 STANDARD                                 ***!
!***                                                                      ***!
!***  NOTES: - Only C Procedure Pointers are not implemented yet. All     ***!
!***           other functionality is provided                            ***!
!***         - Compiler Version testet l_fc_c_10.0.025                    ***!
!***                                                                      ***!
!****************************************************************************!
!****************************************************************************!

#elif defined (__INTEL) || defined(__GFORTRAN)
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
    TYPE(C_PTR)     :: vrr_classes(lib_int_vrr_classes_size,lib_int_vrr_classes_size)
    TYPE(C_PTR)     :: vrr_stack
  END TYPE lib_int

  TYPE, BIND(C) :: lib_deriv
    TYPE(C_PTR)     :: int_stack
    TYPE(C_PTR)     :: PrimQuartet
    TYPE(C_PTR)     :: zero_stack
    TYPE(C_PTR)     :: ABCD(156)
    REAL(C_DOUBLE)  :: AB(3)
    REAL(C_DOUBLE)  :: CD(3)
    TYPE(C_PTR)     :: deriv_classes(12,7,7)
    TYPE(C_PTR)     :: deriv2_classes(144,7,7)
    TYPE(C_PTR)     :: dvrr_classes(7,7)
    TYPE(C_PTR)     :: dvtt_stack
  END TYPE lib_deriv

  !global vars
  TYPE(C_FUNPTR), DIMENSION(0:build_eri_size,0:build_eri_size,0:build_eri_size,0:build_eri_size), BIND(C) :: build_eri
  TYPE(C_FUNPTR), DIMENSION(0:build_deriv1_eri_size,0:build_deriv1_eri_size,&
                            0:build_deriv1_eri_size,0:build_deriv1_eri_size), BIND(C) :: build_deriv1_eri
  
  INTERFACE
! *****************************************************************************
    FUNCTION wrapper_build_eri(na, nb, nc, nd, lib) BIND(C, name="wrap_to_build_eri")
     IMPORT
     TYPE(C_PTR)                    :: wrapper_build_eri
     INTEGER(KIND=C_INT)            :: na, nb, nc, nd
     TYPE(lib_int)                  :: lib
   END FUNCTION wrapper_build_eri
  
! *****************************************************************************
    FUNCTION init_lib(lib, max_am, np) BIND(C, name="init_libint")
      IMPORT
      INTEGER(KIND=C_INT)            :: init_lib
      TYPE(lib_int)                  :: lib
      INTEGER(KIND=C_INT), VALUE     :: max_am
      INTEGER(KIND=C_INT), VALUE     :: np
    END FUNCTION init_lib
    
! *****************************************************************************
    SUBROUTINE init_base() BIND(C, name="init_libint_base")
    END SUBROUTINE init_base

! *****************************************************************************
    SUBROUTINE free_lib_int(lib) BIND(C, name="free_libint")
      IMPORT
      TYPE(lib_int)                  :: lib
    END SUBROUTINE free_lib_int

! *****************************************************************************
    SUBROUTINE init_deriv_base() BIND(C, name="init_libderiv_base")
    END SUBROUTINE init_deriv_base

! *****************************************************************************
    FUNCTION init_deriv1(deriv, max_am, np, ccs) BIND(C, name="init_libderiv1")
      IMPORT
      INTEGER(KIND=C_INT)            :: init_deriv1
      TYPE(lib_deriv)                :: deriv
      INTEGER(KIND=C_INT), VALUE     :: max_am
      INTEGER(KIND=C_INT), VALUE     :: np
      INTEGER(KIND=C_INT), VALUE     :: ccs
    END FUNCTION init_deriv1

! *****************************************************************************
    SUBROUTINE wrapper_build_deriv1_eri(na, nb, nc, nd, deriv) BIND(C, name="wrap_to_build_deriv1_eri")
      IMPORT
      INTEGER(KIND=C_INT)            :: na, nb, nc, nd
      TYPE(lib_deriv)                :: deriv
    END SUBROUTINE wrapper_build_deriv1_eri

! *****************************************************************************
    SUBROUTINE free_lib_deriv(deriv) BIND(C, name="free_libderiv")
      IMPORT
      TYPE(lib_deriv)                  :: deriv
    END SUBROUTINE free_lib_deriv

  END INTERFACE

!***

  CONTAINS

! *****************************************************************************
  SUBROUTINE initialize_libint(lib,max_am,error)
    TYPE(lib_int)                            :: lib
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libint', &
      routineP = moduleN//':'//routineN

    INTEGER(KIND=C_INT)                      :: lib_storage, max_am_local, &
                                                max_prim
    LOGICAL                                  :: failure

    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    CALL init_base()
    lib_storage = init_lib(lib, max_am_local, max_prim)
    IF (lib_storage<0) THEN
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " the angular momentum needed exceeds the value assumed when configuring libint ", &
                        error,failure)
    ENDIF
  END SUBROUTINE initialize_libint

! *****************************************************************************
  SUBROUTINE initialize_libderiv(deriv,max_am,error)
    TYPE(lib_deriv)                          :: deriv
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libderiv', &
      routineP = moduleN//':'//routineN

    INTEGER(KIND=C_INT)                      :: lib_deriv_storage, &
                                                max_am_local, max_classes, &
                                                max_prim
    LOGICAL                                  :: failure

    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    max_classes = 10000!max_am**4

    CALL init_deriv_base()
    lib_deriv_storage = init_deriv1(deriv, max_am_local, max_prim, max_classes)

    IF (lib_deriv_storage<0) THEN
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " the angular momentum needed exceeds the value assumed when configuring libderiv ", &
                        error,failure)
    ENDIF
  END SUBROUTINE initialize_libderiv

! *****************************************************************************
  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                            :: lib

      CALL free_lib_int(lib)
  END SUBROUTINE terminate_libint

! *****************************************************************************
  SUBROUTINE terminate_libderiv(deriv)
    TYPE(lib_deriv)                          :: deriv

      CALL free_lib_deriv(deriv)
  END SUBROUTINE terminate_libderiv
 
! *****************************************************************************
  SUBROUTINE get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_int)                            :: lib
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(:), POINTER          :: p_work
    INTEGER                                  :: a_mysize(1)

    REAL(dp), DIMENSION(:), POINTER          :: p_tmp
    TYPE(C_PTR)                              :: pc_result

    lib%PrimQuartet = C_LOC(prim)
    pc_result= wrapper_build_eri(n_a, n_b, n_c, n_d, lib)
    CALL C_F_POINTER(pc_result, p_tmp, a_mysize)
    p_work => p_tmp 
  END SUBROUTINE get_eris 

! *****************************************************************************
  SUBROUTINE get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_deriv)                          :: deriv
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*&
      nco(n_c)*nco(n_d), 12)                 :: work_forces
    INTEGER                                  :: a_mysize(1)

    INTEGER                                  :: i, k
    REAL(C_DOUBLE), DIMENSION(:), POINTER    :: tmp_data
    TYPE(C_PTR)                              :: pc_result

    deriv%PrimQuartet = C_LOC(prim)
    CALL wrapper_build_deriv1_eri(n_a, n_b, n_c, n_d, deriv)
    DO k=1,12
      IF(k==4 .OR. k==5 .OR. k==6) CYCLE
      pc_result = deriv%ABCD(k)
      CALL C_F_POINTER(pc_result, tmp_data , a_mysize)
      DO i=1,a_mysize(1)
        work_forces(i,k) = tmp_data(i)
      ENDDO
    END DO
  END SUBROUTINE get_derivs 

!****************************************************************************!
!****************************************************************************!
!***                                                                      ***!
!***  WHAT FOLLOWS IS CODE FOR NAG COMPILERS THAT PARTIALLY SUPPORT THE   ***!
!***  ISO_C_BINDING FORTRAN 2003 STANDARD                                 ***!
!***                                                                      ***!
!***  NOTES: - Only C Procedure Pointers are not implemented yet. All     ***!
!***           other functionality is provided                            ***!
!***         - Compiler Version testet 5.1(327)                           ***!
!***                                                                      ***!
!****************************************************************************!
!****************************************************************************!

#elif defined (__NAG)
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
    TYPE(C_PTR)     :: vrr_classes(lib_int_vrr_classes_size,lib_int_vrr_classes_size)
    TYPE(C_PTR)     :: vrr_stack
  END TYPE lib_int

  TYPE, BIND(C) :: lib_deriv
    TYPE(C_PTR)     :: int_stack
    TYPE(C_PTR)     :: PrimQuartet
    TYPE(C_PTR)     :: zero_stack
    TYPE(C_PTR)     :: ABCD(156)
    REAL(C_DOUBLE)  :: AB(3)
    REAL(C_DOUBLE)  :: CD(3)
    TYPE(C_PTR)     :: deriv_classes(12,7,7)
    TYPE(C_PTR)     :: deriv2_classes(144,7,7)
    TYPE(C_PTR)     :: dvrr_classes(7,7)
    TYPE(C_PTR)     :: dvtt_stack
  END TYPE lib_deriv

  !global vars
  TYPE(C_FUNPTR), DIMENSION(0:build_eri_size,0:build_eri_size,0:build_eri_size,0:build_eri_size), BIND(C) :: build_eri
  TYPE(C_FUNPTR), DIMENSION(0:build_deriv1_eri_size,0:build_deriv1_eri_size,&
                            0:build_deriv1_eri_size,0:build_deriv1_eri_size), BIND(C) :: build_deriv1_eri
  
  INTERFACE
! *****************************************************************************
    FUNCTION wrapper_build_eri(na, nb, nc, nd, lib) BIND(C, name="wrap_to_build_eri")
     IMPORT
     TYPE(C_PTR)                    :: wrapper_build_eri
     INTEGER(KIND=C_INT)            :: na, nb, nc, nd
     TYPE(lib_int)                  :: lib
   END FUNCTION wrapper_build_eri
  
! *****************************************************************************
    FUNCTION init_lib(lib, max_am, np) BIND(C, name="init_libint")
      IMPORT
      INTEGER(KIND=C_INT)            :: init_lib
      TYPE(lib_int)                  :: lib
      INTEGER(KIND=C_INT), VALUE     :: max_am
      INTEGER(KIND=C_INT), VALUE     :: np
    END FUNCTION init_lib
    
! *****************************************************************************
    SUBROUTINE init_base() BIND(C, name="init_libint_base")
    END SUBROUTINE init_base

! *****************************************************************************
    SUBROUTINE free_lib_int(lib) BIND(C, name="free_libint")
      IMPORT
      TYPE(lib_int)                  :: lib
    END SUBROUTINE free_lib_int

! *****************************************************************************
    SUBROUTINE init_deriv_base() BIND(C, name="init_libderiv_base")
    END SUBROUTINE init_deriv_base

! *****************************************************************************
    FUNCTION init_deriv1(deriv, max_am, np, ccs) BIND(C, name="init_libderiv1")
      IMPORT
      INTEGER(KIND=C_INT)            :: init_deriv1
      TYPE(lib_deriv)                :: deriv
      INTEGER(KIND=C_INT), VALUE     :: max_am
      INTEGER(KIND=C_INT), VALUE     :: np
      INTEGER(KIND=C_INT), VALUE     :: ccs
    END FUNCTION init_deriv1

! *****************************************************************************
    SUBROUTINE wrapper_build_deriv1_eri(na, nb, nc, nd, deriv) BIND(C, name="wrap_to_build_deriv1_eri")
      IMPORT
      INTEGER(KIND=C_INT)            :: na, nb, nc, nd
      TYPE(lib_deriv)                :: deriv
    END SUBROUTINE wrapper_build_deriv1_eri

! *****************************************************************************
    SUBROUTINE free_lib_deriv(deriv) BIND(C, name="free_libderiv")
      IMPORT
      TYPE(lib_deriv)                  :: deriv
    END SUBROUTINE free_lib_deriv

  END INTERFACE

!***

  CONTAINS

! *****************************************************************************
  SUBROUTINE initialize_libint(lib,max_am,error)
    TYPE(lib_int)                            :: lib
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libint', &
      routineP = moduleN//':'//routineN

    INTEGER(KIND=C_INT)                      :: lib_storage, max_am_local, &
                                                max_prim
    LOGICAL                                  :: failure

    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    CALL init_base()
    lib_storage = init_lib(lib, max_am_local, max_prim)
    IF (lib_storage<0) THEN
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " the angular momentum needed exceeds the value assumed when configuring libint ", &
                        error,failure)
    ENDIF
  END SUBROUTINE initialize_libint

! *****************************************************************************
  SUBROUTINE initialize_libderiv(deriv,max_am,error)
    TYPE(lib_deriv)                          :: deriv
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libderiv', &
      routineP = moduleN//':'//routineN

    INTEGER(KIND=C_INT)                      :: lib_deriv_storage, &
                                                max_am_local, max_classes, &
                                                max_prim
    LOGICAL                                  :: failure

    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    max_classes = 10000!max_am**4

    CALL init_deriv_base()
    lib_deriv_storage = init_deriv1(deriv, max_am_local, max_prim, max_classes)

    IF (lib_deriv_storage<0) THEN
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " the angular momentum needed exceeds the value assumed when configuring libderiv ", &
                        error,failure)
    ENDIF
  END SUBROUTINE initialize_libderiv

! *****************************************************************************
  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                            :: lib

      CALL free_lib_int(lib)
  END SUBROUTINE terminate_libint

! *****************************************************************************
  SUBROUTINE terminate_libderiv(deriv)
    TYPE(lib_deriv)                          :: deriv

      CALL free_lib_deriv(deriv)
  END SUBROUTINE terminate_libderiv
 
! *****************************************************************************
  SUBROUTINE get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_int)                            :: lib
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(:), POINTER          :: p_work
    INTEGER                                  :: a_mysize(1)

    REAL(dp), DIMENSION(:), POINTER          :: p_tmp
    TYPE(C_PTR)                              :: pc_result

    lib%PrimQuartet = C_LOC(prim)
    pc_result= wrapper_build_eri(n_a, n_b, n_c, n_d, lib)
    CALL C_F_POINTER(pc_result, p_tmp, a_mysize)
    p_work => p_tmp 
  END SUBROUTINE get_eris 

! *****************************************************************************
  SUBROUTINE get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_deriv)                          :: deriv
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*&
      nco(n_c)*nco(n_d), 12)                 :: work_forces
    INTEGER                                  :: a_mysize(1)

    INTEGER                                  :: i, k
    REAL(C_DOUBLE), DIMENSION(:), POINTER    :: tmp_data
    TYPE(C_PTR)                              :: pc_result

    deriv%PrimQuartet = C_LOC(prim)
    CALL wrapper_build_deriv1_eri(n_a, n_b, n_c, n_d, deriv)
    DO k=1,12
      IF(k==4 .OR. k==5 .OR. k==6) CYCLE
      pc_result = deriv%ABCD(k)
      CALL C_F_POINTER(pc_result, tmp_data , a_mysize)
      DO i=1,a_mysize(1)
        work_forces(i,k) = tmp_data(i)
      ENDDO
    END DO
  END SUBROUTINE get_derivs 

!****************************************************************************!
!****************************************************************************!
!***                                                                      ***!
!***  WHAT FOLLOWS IS CODE FOR PGI COMPILERS THAT DO NOT SUPPORT THE      ***!
!***  ISO_C_BINDING FORTRAN 2003 STANDARD AT ALL                          ***!
!***                                                                      ***!
!***  NOTES: - C_PTR has been replaced by INTEGER(int_8)                  ***!
!***         - The member prim_data of lib_int/lib_deriv is assigned      ***!
!***           in the C++ wrapping function                               ***!
!***         - Compiler Version testet:                                   ***!
!***             pgf90 6.1-4 64-bit target on catamount target            ***!
!***         - UGLY_BUT_FAST uses CRAY POINTERS. Its faster, because      ***!
!***           no arrays have to be allocated twice in memory             ***!
!***           If not defined, the array p_work has to be allocated       ***!
!***           in the Module hfx_libint_interface::evaluate_eri and       ***!
!***           hfx_libint_interface::evaluate_eri_screen                  ***!
!***                                                                      ***!
!****************************************************************************!
!****************************************************************************!

#elif defined (__XT3)
! *****************************************************************************
  TYPE prim_data
    REAL(dp)        :: F(prim_data_f_size)
    REAL(dp)        :: U(3,6)
    REAL(dp)        :: twozeta_a
    REAL(dp)        :: twozeta_b
    REAL(dp)        :: twozeta_c
    REAL(dp)        :: twozeta_d
    REAL(dp)        :: oo2z
    REAL(dp)        :: oo2n
    REAL(dp)        :: oo2zn
    REAL(dp)        :: poz
    REAL(dp)        :: pon
    REAL(dp)        :: oo2p
    REAL(dp)        :: ss_r12_ss
  END TYPE prim_data
        
! *****************************************************************************
  TYPE lib_int
    INTEGER(int_8)  :: int_stack
    INTEGER(int_8)  :: PrimQuartet
    REAL(dp)        :: AB(3)
    REAL(dp)        :: CD(3)
    INTEGER(int_8)  :: vrr_classes(lib_int_vrr_classes_size,lib_int_vrr_classes_size)
    INTEGER(int_8)  :: vrr_stack
  END TYPE lib_int

! *****************************************************************************
  TYPE lib_deriv
    INTEGER(int_8)  :: int_stack
    INTEGER(int_8)  :: PrimQuartet
    INTEGER(int_8)  :: zero_stack
    REAL(dp)        :: ABCD(156)
    REAL(dp)        :: AB(3)
    REAL(dp)        :: CD(3)
    INTEGER(int_8)  :: deriv_classes(12,7,7)
    INTEGER(int_8)  :: deriv2_classes(144,7,7)
    INTEGER(int_8)  :: dvrr_classes(7,7)
    INTEGER(int_8)  :: dvtt_stack
  END TYPE lib_deriv
!***
  
  CONTAINS
  
! *****************************************************************************
  SUBROUTINE initialize_libint(lib,max_am,error)
    TYPE(lib_int)                            :: lib
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lib_storage, max_am_local, &
                                                max_prim
    LOGICAL                                  :: failure

    EXTERNAL wrapper_init_lib
  
    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    CALL wrapper_init_lib(lib, max_am_local, max_prim, lib_storage)
    IF (lib_storage<0) THEN
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " the angular momentum needed exceeds the value assumed when configuring libint ", &
                        error,failure)
    ENDIF

  END SUBROUTINE initialize_libint

! *****************************************************************************
  SUBROUTINE initialize_libderiv(deriv,max_am,error)
    TYPE(lib_deriv)                          :: deriv
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libderiv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lib_deriv_storage, &
                                                max_am_local, max_classes, &
                                                max_prim
    LOGICAL                                  :: failure

    EXTERNAL wrapper_init_deriv

    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    max_classes = 10000!max_am**4

    CALL wrapper_init_deriv(deriv, max_am_local, max_prim, max_classes, lib_deriv_storage)

    IF (lib_deriv_storage<0) THEN
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " the angular momentum needed exceeds the value assumed when configuring libderiv ", &
                        error,failure)
    ENDIF
  END SUBROUTINE initialize_libderiv

! *****************************************************************************
  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                            :: lib

    EXTERNAL wrapper_free_libint

    CALL wrapper_free_libint(lib)
  END SUBROUTINE terminate_libint

! *****************************************************************************
  SUBROUTINE terminate_libderiv(deriv)
    TYPE(lib_deriv)                          :: deriv

    EXTERNAL wrapper_free_libderiv

    CALL wrapper_free_libderiv(deriv)
  END SUBROUTINE terminate_libderiv

! *****************************************************************************
  SUBROUTINE get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
    INTEGER , INTENT(IN)                    :: n_d, n_c, n_b, n_a
    TYPE(lib_int)                           :: lib
    TYPE(prim_data), TARGET                 :: prim
    REAL(dp), DIMENSION(:), POINTER         :: p_work
    INTEGER                                 :: a_mysize(1)

#if defined (__UGLY_BUT_FAST)   
    REAL(dp), TARGET                        :: tmp(*)
    POINTER(POINTER,tmp)
    EXTERNAL wrapper_build_eri_cray

    CALL wrapper_build_eri_cray(n_a, n_b, n_c, n_d, lib, prim, POINTER)
    p_work => tmp(1:a_mysize(1))
#else
    EXTERNAL wrapper_build_eri

    CALL wrapper_build_eri(n_a, n_b, n_c, n_d, lib, a_mysize(1), p_work(1), prim)
#endif
    
  END SUBROUTINE get_eris

! *****************************************************************************
  SUBROUTINE get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_deriv)                          :: deriv
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*&
      nco(n_c)*nco(n_d), 12)                 :: work_forces
    INTEGER                                  :: a_mysize(1)

    EXTERNAL wrapper_build_deriv1_eri
    
    CALL wrapper_build_deriv1_eri(n_a, n_b, n_c, n_d, deriv, a_mysize(1), work_forces(1,1),prim)

  END SUBROUTINE get_derivs

#elif defined (__PGI) || defined (__AIX)
! *****************************************************************************
  TYPE prim_data
    REAL(dp)        :: F(prim_data_f_size)
    REAL(dp)        :: U(3,6)
    REAL(dp)        :: twozeta_a
    REAL(dp)        :: twozeta_b
    REAL(dp)        :: twozeta_c
    REAL(dp)        :: twozeta_d
    REAL(dp)        :: oo2z
    REAL(dp)        :: oo2n
    REAL(dp)        :: oo2zn
    REAL(dp)        :: poz
    REAL(dp)        :: pon
    REAL(dp)        :: oo2p
    REAL(dp)        :: ss_r12_ss
  END TYPE prim_data
        
! *****************************************************************************
  TYPE lib_int
    INTEGER(int_8)  :: int_stack
    INTEGER(int_8)  :: PrimQuartet
    REAL(dp)        :: AB(3)
    REAL(dp)        :: CD(3)
    INTEGER(int_8)  :: vrr_classes(lib_int_vrr_classes_size,lib_int_vrr_classes_size)
    INTEGER(int_8)  :: vrr_stack
  END TYPE lib_int

! *****************************************************************************
  TYPE lib_deriv
    INTEGER(int_8)  :: int_stack
    INTEGER(int_8)  :: PrimQuartet
    INTEGER(int_8)  :: zero_stack
    REAL(dp)        :: ABCD(156)
    REAL(dp)        :: AB(3)
    REAL(dp)        :: CD(3)
    INTEGER(int_8)  :: deriv_classes(12,7,7)
    INTEGER(int_8)  :: deriv2_classes(144,7,7)
    INTEGER(int_8)  :: dvrr_classes(7,7)
    INTEGER(int_8)  :: dvtt_stack
  END TYPE lib_deriv
!***
  
  CONTAINS
  
! *****************************************************************************
  SUBROUTINE initialize_libint(lib,max_am,error)
    TYPE(lib_int)                            :: lib
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lib_storage, max_am_local, &
                                                max_prim
    LOGICAL                                  :: failure

    EXTERNAL wrapper_init_lib
  
    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    CALL wrapper_init_lib(lib, max_am_local, max_prim, lib_storage)
    IF (lib_storage<0) THEN
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " the angular momentum needed exceeds the value assumed when configuring libint ", &
                        error,failure)
    ENDIF

  END SUBROUTINE initialize_libint

! *****************************************************************************
  SUBROUTINE initialize_libderiv(deriv,max_am,error)
    TYPE(lib_deriv)                          :: deriv
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libderiv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lib_deriv_storage, &
                                                max_am_local, max_classes, &
                                                max_prim
    LOGICAL                                  :: failure

    EXTERNAL wrapper_init_deriv

    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    max_classes = 10000!max_am**4

    CALL wrapper_init_deriv(deriv, max_am_local, max_prim, max_classes, lib_deriv_storage)

    IF (lib_deriv_storage<0) THEN
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " the angular momentum needed exceeds the value assumed when configuring libderiv ", &
                        error,failure)
    ENDIF
  END SUBROUTINE initialize_libderiv

! *****************************************************************************
  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                            :: lib

    EXTERNAL wrapper_free_libint

    CALL wrapper_free_libint(lib)
  END SUBROUTINE terminate_libint

! *****************************************************************************
  SUBROUTINE terminate_libderiv(deriv)
    TYPE(lib_deriv)                          :: deriv

    EXTERNAL wrapper_free_libderiv

    CALL wrapper_free_libderiv(deriv)
  END SUBROUTINE terminate_libderiv

! *****************************************************************************
  SUBROUTINE get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
    INTEGER , INTENT(IN)                    :: n_d, n_c, n_b, n_a
    TYPE(lib_int)                           :: lib
    TYPE(prim_data), TARGET                 :: prim
    REAL(dp), DIMENSION(:), POINTER         :: p_work
    INTEGER                                 :: a_mysize(1)

#if defined (__UGLY_BUT_FAST)   
    REAL(dp), TARGET                        :: tmp(*)
    POINTER(POINTER,tmp)
    EXTERNAL wrapper_build_eri_cray

    CALL wrapper_build_eri_cray(n_a, n_b, n_c, n_d, lib, prim, POINTER)
    p_work => tmp(1:a_mysize(1))
#else
    EXTERNAL wrapper_build_eri
    CALL wrapper_build_eri(n_a, n_b, n_c, n_d, lib, a_mysize(1), p_work(1), prim)
#endif
    
  END SUBROUTINE get_eris

! *****************************************************************************
  SUBROUTINE get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_deriv)                          :: deriv
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*&
      nco(n_c)*nco(n_d), 12)                 :: work_forces
    INTEGER                                  :: a_mysize(1)

    EXTERNAL wrapper_build_deriv1_eri
    
    CALL wrapper_build_deriv1_eri(n_a, n_b, n_c, n_d, deriv, a_mysize(1), work_forces(1,1),prim)

  END SUBROUTINE get_derivs

!COMPILERS
#endif 

!****************************************************************************!
!****************************************************************************!
!***                                                                      ***!
!***  WHAT FOLLOWS IS CODE THAT USES BOGUS SUBROUTINES AND TYPES IN       ***!
!***  ORDER TO ALLOW ALL COMPILERS TO COMPILE CP2K                        ***!
!***                                                                      ***!
!****************************************************************************!
!****************************************************************************!
!LIBINT
#else
  ! some bogus types to be able to compile without libint being present
! *****************************************************************************
  TYPE prim_data
    REAL(dp)        :: F(prim_data_f_size)
    REAL(dp)        :: U(3,6)
    REAL(dp)        :: twozeta_a
    REAL(dp)        :: twozeta_b
    REAL(dp)        :: twozeta_c
    REAL(dp)        :: twozeta_d
    REAL(dp)        :: oo2z
    REAL(dp)        :: oo2n
    REAL(dp)        :: oo2zn
    REAL(dp)        :: poz
    REAL(dp)        :: pon
    REAL(dp)        :: oo2p
    REAL(dp)        :: ss_r12_ss
  END TYPE prim_data
        
! *****************************************************************************
  TYPE lib_int
    REAL(dp)        :: int_stack
    TYPE(prim_data) :: PrimQuartet
    REAL(dp)        :: AB(3)
    REAL(dp)        :: CD(3)
    REAL(dp)        :: vrr_classes(lib_int_vrr_classes_size,lib_int_vrr_classes_size)
    REAL(dp)        :: vrr_stack
  END TYPE lib_int

! *****************************************************************************
  TYPE lib_deriv
    REAL(dp)        :: int_stack
    TYPE(prim_data) :: PrimQuartet
    REAL(dp)        :: zero_stack
    REAL(dp)        :: ABCD(156)
    REAL(dp)        :: AB(3)
    REAL(dp)        :: CD(3)
    REAL(dp)        :: deriv_classes(12,7,7)
    REAL(dp)        :: deriv2_classes(144,7,7)
    REAL(dp)        :: dvrr_classes(7,7)
    REAL(dp)        :: dvtt_stack
  END TYPE lib_deriv
!***

  CONTAINS

! *****************************************************************************
  SUBROUTINE initialize_libint(lib,max_am,error)
    TYPE(lib_int)                            :: lib
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libint', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                           " This compiler does not support interfacing to libint." ,error,failure)

  END SUBROUTINE initialize_libint

! *****************************************************************************
  SUBROUTINE initialize_libderiv(deriv,max_am,error)
    TYPE(lib_deriv)                          :: deriv
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libderiv', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                           " This compiler does not support interfacing to libint." ,error,failure)

  END SUBROUTINE initialize_libderiv

! *****************************************************************************
  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                            :: lib

  END SUBROUTINE terminate_libint

! *****************************************************************************
  SUBROUTINE terminate_libderiv(deriv)
    TYPE(lib_deriv)                          :: deriv

  END SUBROUTINE terminate_libderiv

! *****************************************************************************
  SUBROUTINE get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_int)                            :: lib
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(:), POINTER          :: p_work
    INTEGER                                  :: a_mysize(1)

  END SUBROUTINE get_eris

! *****************************************************************************
  SUBROUTINE get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_deriv)                          :: deriv
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*&
      nco(n_c)*nco(n_d), 12)                 :: work_forces
    INTEGER                                  :: a_mysize(1)

  END SUBROUTINE get_derivs
!LIBINT
#endif 

END MODULE hfx_libint_wrapper
