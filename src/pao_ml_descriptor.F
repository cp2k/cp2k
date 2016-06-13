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
                                              pao_descriptor_type,&
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

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_ml_calc_descriptor'

      INTEGER                                            :: handle, idesc, ikind, jatom, jkind, N, &
                                                            natoms, ndesc
      REAL(dp)                                           :: beta, weight
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: V_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: V_block
      REAL(dp), DIMENSION(3)                             :: Ra, Rab, Rb
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(pao_descriptor_type), DIMENSION(:), POINTER   :: pao_descriptors

      CALL timeset(routineN, handle)

      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set, pao_descriptors=pao_descriptors)

      N = basis_set%nsgf
      natoms = SIZE(particle_set)
      ndesc = SIZE(pao_descriptors)
      IF (ndesc == 0) CPABORT("No PAO_DESCRIPTOR section found")

      ALLOCATE (descriptor(N*ndesc), V_block(N, N), V_evals(N))

      DO idesc = 1, ndesc

         ! construct matrix V_block from neighboring atoms
         V_block = 0.0_dp
         DO jatom = 1, natoms
            IF (jatom == iatom) CYCLE
            Ra = particle_set(iatom)%r
            Rb = particle_set(jatom)%r
            Rab = pbc(ra, rb, cell)
            CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=jkind)
            CALL get_qs_kind(qs_kind_set(jkind), pao_descriptors=pao_descriptors)
            IF (SIZE(pao_descriptors) /= ndesc) &
               CPABORT("Not all KINDs have the same number of PAO_DESCRIPTOR sections")
            weight = pao_descriptors(idesc)%weight
            beta = pao_descriptors(idesc)%beta
            CALL calc_gaussian(basis_set, V_block, Rab, lpot=0, beta=beta, weight=weight)
         ENDDO

         ! diagonalize V_block
         CALL diamat_all(V_block, V_evals)

         ! use eigenvalues of V_block as descriptor
         descriptor((idesc-1)*N+1:idesc*N) = V_evals(:)
      ENDDO

      DEALLOCATE (V_block, V_evals)

      CALL timestop(handle)
   END SUBROUTINE pao_ml_calc_descriptor

END MODULE pao_ml_descriptor
