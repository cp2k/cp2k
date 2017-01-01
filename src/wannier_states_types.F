!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
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
