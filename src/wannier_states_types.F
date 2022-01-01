!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief defines the type needed for computing wannier states expectations
!> \author alin m elena
! **************************************************************************************************
MODULE wannier_states_types
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: wannier_centres_type

   TYPE wannier_centres_type
      REAL(kind=dp), POINTER, DIMENSION(:)        :: WannierHamDiag => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:, :)      :: centres => NULL()
   END TYPE wannier_centres_type

END MODULE wannier_states_types
