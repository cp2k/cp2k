!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Holds information on atomic properties
!> \par History
!>      07.2011 created
!> \author JHU
! **************************************************************************************************
MODULE atprop_types

   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atprop_types'

   PUBLIC :: atprop_type, atprop_create, atprop_release, atprop_init
   PUBLIC :: atprop_array_init, atprop_array_add

! **************************************************************************************************
!> \brief type for the atomic properties
!> \author fawzi
! **************************************************************************************************
   TYPE atprop_type
      LOGICAL                                   :: energy = .FALSE.
      LOGICAL                                   :: stress = .FALSE.
      REAL(KIND=dp), DIMENSION(:), POINTER      :: atener => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER      :: ateb => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER      :: atexc => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER      :: ateself => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER      :: atevdw => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER      :: ategcp => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER      :: atecc => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER      :: ate1c => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER      :: atecoul => NULL()
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER  :: atstress => NULL()
   END TYPE atprop_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param atprop_env ...
! **************************************************************************************************
   SUBROUTINE atprop_create(atprop_env)
      TYPE(atprop_type), POINTER                         :: atprop_env

      CALL atprop_release(atprop_env)
      ALLOCATE (atprop_env)

   END SUBROUTINE atprop_create

! **************************************************************************************************
!> \brief ...
!> \param atprop_env ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE atprop_init(atprop_env, natom)
      TYPE(atprop_type), INTENT(INOUT)                   :: atprop_env
      INTEGER, INTENT(IN)                                :: natom

      IF (atprop_env%energy) THEN
         CALL atprop_array_init(atprop_env%atener, natom)
         CALL atprop_array_release(atprop_env%ateb)
         CALL atprop_array_release(atprop_env%atevdw)
         CALL atprop_array_release(atprop_env%ategcp)
         CALL atprop_array_release(atprop_env%atecc)
         CALL atprop_array_release(atprop_env%atecoul)
         CALL atprop_array_release(atprop_env%ateself)
         CALL atprop_array_release(atprop_env%atexc)
         CALL atprop_array_release(atprop_env%ate1c)
      END IF

      IF (atprop_env%stress) THEN
         IF (ASSOCIATED(atprop_env%atstress)) THEN
            CPASSERT(SIZE(atprop_env%atstress, 3) == natom)
         ELSE
            ALLOCATE (atprop_env%atstress(3, 3, natom))
         END IF
         atprop_env%atstress = 0._dp
      END IF

   END SUBROUTINE atprop_init

! **************************************************************************************************
!> \brief ...
!> \param atarray ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE atprop_array_init(atarray, natom)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atarray
      INTEGER, INTENT(IN)                                :: natom

      IF (ASSOCIATED(atarray)) THEN
         CPASSERT(SIZE(atarray) == natom)
      ELSE
         ALLOCATE (atarray(natom))
      END IF
      atarray = 0._dp

   END SUBROUTINE atprop_array_init

! **************************************************************************************************
!> \brief ...
!> \param atarray ...
! **************************************************************************************************
   SUBROUTINE atprop_array_release(atarray)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atarray

      IF (ASSOCIATED(atarray)) THEN
         DEALLOCATE (atarray)
      END IF

   END SUBROUTINE atprop_array_release

! **************************************************************************************************
!> \brief ...
!> \param array_a ...
!> \param array_b ...
! **************************************************************************************************
   SUBROUTINE atprop_array_add(array_a, array_b)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: array_a
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), POINTER   :: array_b

      IF (ASSOCIATED(array_b)) THEN
         array_a = array_a + array_b
      END IF

   END SUBROUTINE atprop_array_add

! **************************************************************************************************
!> \brief releases the atprop
!> \param atprop_env the object to release
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE atprop_release(atprop_env)
      TYPE(atprop_type), POINTER                         :: atprop_env

      IF (ASSOCIATED(atprop_env)) THEN
         ! energy
         CALL atprop_array_release(atprop_env%atener)
         CALL atprop_array_release(atprop_env%ateb)
         CALL atprop_array_release(atprop_env%ateself)
         CALL atprop_array_release(atprop_env%atexc)
         CALL atprop_array_release(atprop_env%atevdw)
         CALL atprop_array_release(atprop_env%ategcp)
         CALL atprop_array_release(atprop_env%atecc)
         CALL atprop_array_release(atprop_env%ate1c)
         CALL atprop_array_release(atprop_env%atecoul)
         ! stress
         IF (ASSOCIATED(atprop_env%atstress)) THEN
            DEALLOCATE (atprop_env%atstress)
         END IF
         ! atprop type
         DEALLOCATE (atprop_env)
      END IF
      NULLIFY (atprop_env)
   END SUBROUTINE atprop_release

END MODULE atprop_types
