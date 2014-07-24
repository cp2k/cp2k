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

  USE, INTRINSIC :: ISO_C_BINDING
  USE kinds,                           ONLY: dp
  USE hfx_libint_wrapper_types,        ONLY: lib_int,&
                                             lib_deriv,&
                                             prim_data,&
                                             build_eri_size,&
                                             build_deriv1_eri_size
  USE orbital_pointers,                ONLY: nco

  USE f77_blas
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_libint, terminate_libint,&
         initialize_libderiv, get_eris, get_derivs, terminate_libderiv

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_libint_wrapper'

#if(__LIBINT)

  TYPE(C_FUNPTR), DIMENSION(0:build_eri_size,0:build_eri_size,0:build_eri_size,0:build_eri_size), BIND(C) :: build_eri
  TYPE(C_FUNPTR), DIMENSION(0:build_deriv1_eri_size,0:build_deriv1_eri_size,&
                            0:build_deriv1_eri_size,0:build_deriv1_eri_size), BIND(C) :: build_deriv1_eri

  INTERFACE
    FUNCTION build(lib, np) BIND(C) 
      IMPORT
      TYPE(C_PTR)                    :: build
      TYPE(lib_int)                  :: lib
      INTEGER(KIND=C_INT), VALUE     :: np
    END FUNCTION build

    FUNCTION init_lib(lib, max_am, np) BIND(C, name="init_libint")
      IMPORT
      INTEGER(KIND=C_INT)            :: init_lib
      TYPE(lib_int)                  :: lib
      INTEGER(KIND=C_INT), VALUE     :: max_am
      INTEGER(KIND=C_INT), VALUE     :: np
    END FUNCTION init_lib

    SUBROUTINE init_base() BIND(C, name="init_libint_base")
    END SUBROUTINE init_base

    SUBROUTINE free_lib_int(lib) BIND(C, name="free_libint")
      IMPORT
      TYPE(lib_int)                  :: lib
    END SUBROUTINE free_lib_int

    SUBROUTINE init_deriv_base() BIND(C, name="init_libderiv_base")
    END SUBROUTINE init_deriv_base

    FUNCTION init_deriv1(deriv, max_am, np, ccs) BIND(C, name="init_libderiv1")
      IMPORT
      INTEGER(KIND=C_INT)            :: init_deriv1
      TYPE(lib_deriv)                :: deriv
      INTEGER(KIND=C_INT), VALUE     :: max_am
      INTEGER(KIND=C_INT), VALUE     :: np
      INTEGER(KIND=C_INT), VALUE     :: ccs
    END FUNCTION init_deriv1

    SUBROUTINE build_deriv1(deriv, np) BIND(C)
      IMPORT
      TYPE(lib_deriv)                :: deriv
      INTEGER(KIND=C_INT),VALUE      :: np
    END SUBROUTINE build_deriv1

    SUBROUTINE free_lib_deriv(deriv) BIND(C, name="free_libderiv")
      IMPORT
      TYPE(lib_deriv)                  :: deriv
    END SUBROUTINE free_lib_deriv
  END INTERFACE

  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param lib ...
!> \param max_am ...
!> \param error ...
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
!> \brief ...
!> \param deriv ...
!> \param max_am ...
!> \param error ...
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
    max_classes = nco(max_am)**4

    CALL init_deriv_base()
    lib_deriv_storage = init_deriv1(deriv, max_am_local, max_prim, max_classes)
    IF (lib_deriv_storage<0) THEN
      CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                     " the angular momentum needed exceeds the value assumed when configuring libderiv ", &
                     error,failure)
    ENDIF
  END SUBROUTINE initialize_libderiv

! *****************************************************************************
!> \brief ...
!> \param lib ...
! *****************************************************************************
  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                            :: lib

    CALL free_lib_int(lib)
  END SUBROUTINE terminate_libint

! *****************************************************************************
!> \brief ...
!> \param deriv ...
! *****************************************************************************
  SUBROUTINE terminate_libderiv(deriv)
    TYPE(lib_deriv)                          :: deriv

    CALL free_lib_deriv(deriv)
  END SUBROUTINE terminate_libderiv

! *****************************************************************************
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param prim ...
!> \param p_work ...
!> \param a_mysize ...
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
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param deriv ...
!> \param prim ...
!> \param work_forces ...
!> \param a_mysize ...
! *****************************************************************************
  SUBROUTINE get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_deriv)                          :: deriv
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*&
      nco(n_c)*nco(n_d), 12)                 :: work_forces
    INTEGER                                  :: a_mysize(1)

    PROCEDURE(build_deriv1), POINTER         :: pbuild_deriv1
    TYPE(C_PTR)                              :: pc_result
    REAL(C_DOUBLE), DIMENSION(:), POINTER    :: tmp_data
    INTEGER                                  :: i, k

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

#else

!****************************************************************************!
!****************************************************************************!
!***                                                                      ***!
!***  WHAT FOLLOWS IS CODE TO PROVIDE STUB ROUTINES IN ABSENCE OF __LIBINT **!
!***                                                                      ***!
!****************************************************************************!
!****************************************************************************!

  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param lib ...
!> \param max_am ...
!> \param error ...
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
                   " This CP2K executable has not been linked against the library libint, required for HFX.", &
                   error,failure)

  END SUBROUTINE initialize_libint

! *****************************************************************************
!> \brief ...
!> \param deriv ...
!> \param max_am ...
!> \param error ...
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
                   " This CP2K executable has not been linked against the library libint, required for HFX.", &
                   error,failure)

  END SUBROUTINE initialize_libderiv

! *****************************************************************************
!> \brief ...
!> \param lib ...
! *****************************************************************************
  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                            :: lib

  END SUBROUTINE terminate_libint

! *****************************************************************************
!> \brief ...
!> \param deriv ...
! *****************************************************************************
  SUBROUTINE terminate_libderiv(deriv)
    TYPE(lib_deriv)                          :: deriv

  END SUBROUTINE terminate_libderiv

! *****************************************************************************
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param prim ...
!> \param p_work ...
!> \param a_mysize ...
! *****************************************************************************
  SUBROUTINE get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
    INTEGER, INTENT(IN)                      :: n_d, n_c, n_b, n_a
    TYPE(lib_int)                            :: lib
    TYPE(prim_data), TARGET                  :: prim
    REAL(dp), DIMENSION(:), POINTER          :: p_work
    INTEGER                                  :: a_mysize(1)

  END SUBROUTINE get_eris

! *****************************************************************************
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param deriv ...
!> \param prim ...
!> \param work_forces ...
!> \param a_mysize ...
! *****************************************************************************
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
