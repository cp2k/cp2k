!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the Libint-Library or a c++ wrapper.
!> \par History
!>      11.2007 created [Manuel Guidon]
!>      10.2009 refactored [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
MODULE hfx_libint_wrapper

   USE ISO_C_BINDING,                   ONLY: C_F_POINTER,&
                                              C_F_PROCPOINTER
   USE hfx_libint_wrapper_types,        ONLY: libint_t
   USE kinds,                           ONLY: dp
   USE libint_f,                        ONLY: &
        libint2_build, libint2_build_eri, libint2_build_eri1, libint2_cleanup_eri, &
        libint2_cleanup_eri1, libint2_init_eri, libint2_init_eri1, libint2_static_cleanup, &
        libint2_static_init, libint_t
   USE orbital_pointers,                ONLY: nco
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: get_eris, get_derivs, libint2_init_eri, libint2_init_eri1, libint2_cleanup_eri, libint2_cleanup_eri1, &
             libint2_static_init, libint2_static_cleanup

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_libint_wrapper'

#if(__LIBINT)

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE get_eris(n_d, n_c, n_b, n_a, lib, p_work, a_mysize)
      INTEGER, INTENT(IN)                                :: n_d, n_c, n_b, n_a
      TYPE(libint_t), DIMENSION(1)                       :: lib
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      INTEGER                                            :: a_mysize(1)

      PROCEDURE(libint2_build), POINTER               :: pbuild

      CALL C_F_PROCPOINTER(libint2_build_eri(n_d, n_c, n_b, n_a), pbuild)
      CALL pbuild(lib)

      CALL C_F_POINTER(lib(1)%targets(1), p_work, SHAPE=a_mysize)

   END SUBROUTINE get_eris

! **************************************************************************************************
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param work_forces ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE get_derivs(n_d, n_c, n_b, n_a, lib, work_forces, a_mysize)
      INTEGER, INTENT(IN)                                :: n_d, n_c, n_b, n_a
      TYPE(libint_t), DIMENSION(1)                       :: lib
      REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*nco(n_c)*nco&
         (n_d), 12)                                      :: work_forces
      INTEGER                                            :: a_mysize(1)

      REAL(dp), DIMENSION(:), POINTER                    :: p_work

      PROCEDURE(libint2_build), POINTER               :: pbuild
      INTEGER :: i, k

      CALL C_F_PROCPOINTER(libint2_build_eri1(n_d, n_c, n_b, n_a), pbuild)
      CALL pbuild(lib)

      DO k = 1, 12
         IF (k == 4 .OR. k == 5 .OR. k == 6) CYCLE
         CALL C_F_POINTER(lib(1)%targets(k), p_work, SHAPE=a_mysize)
         DO i = 1, a_mysize(1)
            work_forces(i, k) = p_work(i)
         ENDDO
      ENDDO

   END SUBROUTINE get_derivs

#else

! TODO: recover fall back routine in absence of __LIBINT
!!****************************************************************************!
!!****************************************************************************!
!!***                                                                      ***!
!!***  WHAT FOLLOWS IS CODE TO PROVIDE STUB ROUTINES IN ABSENCE OF __LIBINT **!
!!***                                                                      ***!
!!****************************************************************************!
!!****************************************************************************!
!
!CONTAINS
!
!! **************************************************************************************************
!!> \brief ...
!!> \param lib ...
!!> \param max_am ...
!! **************************************************************************************************
!   !SUBROUTINE initialize_libint(lib, max_am)
!   !   TYPE(lib_int)                                      :: lib
!   !   INTEGER                                            :: max_am
!
!   !   CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libint', &
!   !      routineP = moduleN//':'//routineN
!
!   !   MARK_USED(lib)
!   !   MARK_USED(max_am)
!   !   CPABORT("This CP2K executable has not been linked against the required library libint (HFX).")
!
!   !END SUBROUTINE initialize_libint
!
!! **************************************************************************************************
!!> \brief ...
!!> \param deriv ...
!!> \param max_am ...
!! **************************************************************************************************
!   !SUBROUTINE initialize_libderiv(deriv, max_am)
!   !   TYPE(lib_deriv)                                    :: deriv
!   !   INTEGER                                            :: max_am
!
!   !   CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libderiv', &
!   !      routineP = moduleN//':'//routineN
!
!   !   MARK_USED(deriv)
!   !   MARK_USED(max_am)
!   !   CPABORT("This CP2K executable has not been linked against the required library libint (HFX).")
!
!   !END SUBROUTINE initialize_libderiv
!
!! **************************************************************************************************
!!> \brief ...
!!> \param lib ...
!! **************************************************************************************************
!   !SUBROUTINE terminate_libint(lib)
!   !   TYPE(lib_int)                                      :: lib
!
!   !   MARK_USED(lib)
!   !   CPABORT("This CP2K executable has not been linked against the required library libint (HFX).")
!
!   !END SUBROUTINE terminate_libint
!
!! **************************************************************************************************
!!> \brief ...
!!> \param deriv ...
!! **************************************************************************************************
!   !SUBROUTINE terminate_libderiv(deriv)
!   !   TYPE(lib_deriv)                                    :: deriv
!
!   !   MARK_USED(deriv)
!   !   CPABORT("This CP2K executable has not been linked against the required library libint (HFX).")
!
!   !END SUBROUTINE terminate_libderiv
!
!! **************************************************************************************************
!!> \brief ...
!!> \param n_d ...
!!> \param n_c ...
!!> \param n_b ...
!!> \param n_a ...
!!> \param lib ...
!!> \param prim ...
!!> \param p_work ...
!!> \param a_mysize ...
!! **************************************************************************************************
!   SUBROUTINE get_eris(n_d, n_c, n_b, n_a, lib, p_work, a_mysize)
!      INTEGER, INTENT(IN)                                :: n_d, n_c, n_b, n_a
!      TYPE(libint_t), DIMENSION(1)                        :: lib
!      REAL(dp), DIMENSION(:), POINTER                    :: p_work
!      INTEGER                                            :: a_mysize(1)
!
!      PROCEDURE(libint2_build), POINTER               :: pbuild
!
!
!      MARK_USED(n_a)
!      MARK_USED(n_b)
!      MARK_USED(n_c)
!      MARK_USED(n_d)
!      MARK_USED(lib)
!      MARK_USED(p_work)
!      MARK_USED(a_mysize)
!      CPABORT("This CP2K executable has not been linked against the required library libint (HFX).")
!
!   END SUBROUTINE get_eris
!
!! **************************************************************************************************
!!> \brief ...
!!> \param n_d ...
!!> \param n_c ...
!!> \param n_b ...
!!> \param n_a ...
!!> \param deriv ...
!!> \param prim ...
!!> \param work_forces ...
!!> \param a_mysize ...
!! **************************************************************************************************
!   SUBROUTINE get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
!      INTEGER, INTENT(IN)                                :: n_d, n_c, n_b, n_a
!      TYPE(lib_deriv)                                    :: deriv
!      TYPE(prim_data), TARGET                            :: prim
!      REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*nco(n_c)*nco&
!         (n_d), 12)                                      :: work_forces
!      INTEGER                                            :: a_mysize(1)
!
!      MARK_USED(n_a)
!      MARK_USED(n_b)
!      MARK_USED(n_c)
!      MARK_USED(n_d)
!      MARK_USED(deriv)
!      MARK_USED(prim)
!      MARK_USED(work_forces)
!      MARK_USED(a_mysize)
!      CPABORT("This CP2K executable has not been linked against the required library libint (HFX).")
!
!   END SUBROUTINE get_derivs
#endif

END MODULE hfx_libint_wrapper
