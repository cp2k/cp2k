!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Feature vectors for describing chemical environments in a rotationally invariant fashion.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_ml_descriptor
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
   USE pao_potentials,                  ONLY: pao_build_V_neighbors
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_ml_descriptor'

   PUBLIC :: pao_ml_calc_descriptor

CONTAINS

! **************************************************************************************************
!> \brief Calulates a descriptor based on the eigenvalues of V_neighbors
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param descriptor ...
! **************************************************************************************************
   SUBROUTINE pao_ml_calc_descriptor(particle_set, qs_kind_set, cell, iatom, descriptor)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: descriptor

      INTEGER                                            :: N
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: V_neighbors

      CPASSERT(.NOT. ALLOCATED(descriptor))

      CALL pao_build_V_neighbors(particle_set, qs_kind_set, cell, iatom, V_neighbors)
      N = SIZE(V_neighbors, 1)
      ALLOCATE (descriptor(N))
      CALL diamat_all(V_neighbors, descriptor)
      DEALLOCATE (V_neighbors)

   END SUBROUTINE pao_ml_calc_descriptor

END MODULE pao_ml_descriptor
