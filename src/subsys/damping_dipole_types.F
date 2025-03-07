!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \author Rodolphe Vuilleumier (29.12.2009)
! **************************************************************************************************
MODULE damping_dipole_types

   USE kinds,                           ONLY: default_string_length,&
                                              dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters (only in this module)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'damping_dipole_types'

! *** Global public parameters

   INTEGER, PUBLIC, PARAMETER :: no_damping = -1, &
                                 tang_toennies = 1

! *** Define the damping types ***
! **************************************************************************************************
   TYPE damping_info_type
      CHARACTER(LEN=default_string_length)   :: atm_name1 = "", atm_name2 = ""
      CHARACTER(LEN=default_string_length)   :: dtype = ""
      INTEGER                                 :: order = -1
      REAL(KIND=dp)                           :: bij = 0.0_dp, cij = 0.0_dp
   END TYPE damping_info_type
! **************************************************************************************************
   TYPE damping_type
      INTEGER :: itype = no_damping
      INTEGER :: order = 1
      REAL(KIND=dp) :: bij = HUGE(0.0_dp), cij = 0.0_dp
   END TYPE damping_type

   TYPE damping_p_type
      TYPE(damping_type), DIMENSION(:), POINTER :: damp => NULL()
   END TYPE

! **************************************************************************************************

! *** Public data types ***

   PUBLIC :: damping_info_type, damping_type

! *** Public subroutines ***

   PUBLIC :: damping_p_type, damping_p_create, damping_p_release

CONTAINS

! **************************************************************************************************
!> \brief Creates Data-structure that contains damping information
!> \param damping ...
!> \param nkinds ...
!> \author Rodolphe Vuilleumier
! **************************************************************************************************
   SUBROUTINE damping_p_create(damping, nkinds)
      TYPE(damping_p_type), POINTER                      :: damping
      INTEGER, INTENT(IN)                                :: nkinds

      CPASSERT(.NOT. ASSOCIATED(damping))
      ALLOCATE (damping)
      ALLOCATE (damping%damp(nkinds))

   END SUBROUTINE damping_p_create
! **************************************************************************************************
!> \brief Release Data-structure that contains damping information
!> \param damping ...
!> \author Rodolphe Vuilleumier [RV]
! **************************************************************************************************
   SUBROUTINE damping_p_release(damping)
      TYPE(damping_p_type), POINTER                      :: damping

      IF (ASSOCIATED(damping)) THEN
         IF (ASSOCIATED(damping%damp)) THEN
            DEALLOCATE (damping%damp)
         END IF
         DEALLOCATE (damping)
      END IF
      NULLIFY (damping)

   END SUBROUTINE damping_p_release

! **************************************************************************************************
END MODULE damping_dipole_types
