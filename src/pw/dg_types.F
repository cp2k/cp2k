!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      none
! **************************************************************************************************
MODULE dg_types

   USE dg_rho0_types,                   ONLY: dg_rho0_create,&
                                              dg_rho0_release,&
                                              dg_rho0_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

! Information on the assignment function for Ewald
! **************************************************************************************************
   TYPE dg_type
      PRIVATE
      TYPE(dg_rho0_type), POINTER :: dg_rho0 => NULL()
   END TYPE dg_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dg_types'

   PRIVATE
   PUBLIC :: dg_type, dg_get, &
             dg_set, dg_release, &
             dg_create

CONTAINS

! **************************************************************************************************
!> \brief   Get the dg_type
!> \param dg ...
!> \param dg_rho0 ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE dg_get(dg, dg_rho0)
      TYPE(dg_type), POINTER                             :: dg
      TYPE(dg_rho0_type), OPTIONAL, POINTER              :: dg_rho0

      IF (PRESENT(dg_rho0)) dg_rho0 => dg%dg_rho0

   END SUBROUTINE dg_get

! **************************************************************************************************
!> \brief   create the dg structure
!> \param dg ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE dg_create(dg)
      TYPE(dg_type), INTENT(OUT)                         :: dg

      CALL dg_rho0_create(dg%dg_rho0)

   END SUBROUTINE dg_create

! **************************************************************************************************
!> \brief releases the given dg_type
!> \param dg the dg_type to release
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
!> \note
!>      see doc/ReferenceCounting.html
! **************************************************************************************************
   SUBROUTINE dg_release(dg)
      TYPE(dg_type), INTENT(INOUT)                       :: dg

      CALL dg_rho0_release(dg%dg_rho0)
   END SUBROUTINE dg_release

! **************************************************************************************************
!> \brief   Set the double grid environment
!> \param dg ...
!> \param dg_rho0 ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE dg_set(dg, dg_rho0)
      TYPE(dg_type), POINTER                             :: dg
      TYPE(dg_rho0_type), OPTIONAL, POINTER              :: dg_rho0

      IF (PRESENT(dg_rho0)) THEN
         CALL dg_rho0_release(dg%dg_rho0)
         dg%dg_rho0 => dg_rho0
      END IF
   END SUBROUTINE dg_set

END MODULE dg_types
