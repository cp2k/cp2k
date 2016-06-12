!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Feature vectors for describing chemical environments in a rotationally invariant fashion.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_ml_descriptor
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
   USE pao_potentials,                  ONLY: calc_gaussian
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
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

      INTEGER                                            :: ikind, jatom, N, natoms
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: V_neighbors
      REAL(dp), DIMENSION(3)                             :: Ra, Rab, Rb
      TYPE(gto_basis_set_type), POINTER                  :: basis_set

      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
      N = basis_set%nsgf
      natoms = SIZE(particle_set)

      ALLOCATE (V_neighbors(N, N))
      V_neighbors = 0.0_dp
      DO jatom = 1, natoms
         IF (jatom == iatom) CYCLE
         Ra = particle_set(iatom)%r
         Rb = particle_set(jatom)%r
         Rab = pbc(ra, rb, cell)
         ! TODO read beta and weight from input
         CALL calc_gaussian(basis_set, V_neighbors, Rab, lpot=0, beta=1.0_dp, weight=1.0_dp)
      ENDDO

      ALLOCATE (descriptor(N))
      CALL diamat_all(V_neighbors, descriptor)
      DEALLOCATE (V_neighbors)

   END SUBROUTINE pao_ml_calc_descriptor

END MODULE pao_ml_descriptor
