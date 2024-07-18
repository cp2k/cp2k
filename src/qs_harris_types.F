!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for a for a Harris model calculation
!> \par History
!>       2024.07 created
!> \author JGH
! **************************************************************************************************
MODULE qs_harris_types
#include "./base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_harris_types'

   PUBLIC :: harris_type, harris_env_release

! *****************************************************************************
!> \brief Contains information on the Harris method
!> \par History
!>       07.2024 created
!> \author JGH
! *****************************************************************************
   TYPE harris_type
      INTEGER                                          :: energy_functional = 0
      INTEGER                                          :: density_source = 0
      INTEGER                                          :: orbital_basis = 0
      !
      LOGICAL                                          :: debug_forces = .FALSE.
      LOGICAL                                          :: debug_stress = .FALSE.
   END TYPE harris_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param harris_env ...
! **************************************************************************************************
   SUBROUTINE harris_env_release(harris_env)
      TYPE(harris_type), POINTER                         :: harris_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'harris_env_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(harris_env)) THEN

      END IF

      NULLIFY (harris_env)

      CALL timestop(handle)

   END SUBROUTINE harris_env_release

END MODULE qs_harris_types
