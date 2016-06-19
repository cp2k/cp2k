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
   USE pao_input,                       ONLY: pao_ml_desc_pot,&
                                              pao_ml_desc_r12
   USE pao_potentials,                  ONLY: pao_calc_gaussian
   USE pao_types,                       ONLY: pao_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              pao_descriptor_type,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_ml_descriptor'

   PUBLIC :: pao_ml_calc_descriptor, pao_ml_descriptor_forces

CONTAINS

! **************************************************************************************************
!> \brief Calulates a descriptor for chemical environment of given atom
!> \param pao ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param descriptor ...
! **************************************************************************************************
   SUBROUTINE pao_ml_calc_descriptor(pao, particle_set, qs_kind_set, cell, iatom, descriptor)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: descriptor

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_ml_calc_descriptor'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%ml_descriptor)
      CASE (pao_ml_desc_pot)
         CALL calc_descriptor_pot(particle_set, qs_kind_set, cell, iatom, descriptor)
      CASE (pao_ml_desc_r12)
         CALL calc_descriptor_r12(particle_set, qs_kind_set, cell, iatom, descriptor)
      CASE DEFAULT
         CPABORT("PAO: unkown descriptor")
      END SELECT

      CALL timestop(handle)
   END SUBROUTINE pao_ml_calc_descriptor

! **************************************************************************************************
!> \brief Calulates a descriptor for chemical environment of given atom
!> \param pao ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param descriptor_grad ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE pao_ml_descriptor_forces(pao, particle_set, qs_kind_set, cell, iatom, descriptor_grad, forces)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: descriptor_grad
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: forces

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_ml_descriptor_forces'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%ml_descriptor)
      CASE (pao_ml_desc_pot)
         CALL forces_descriptor_pot(particle_set, qs_kind_set, cell, iatom, descriptor_grad, forces)
      CASE (pao_ml_desc_r12)
         CALL forces_descriptor_r12(particle_set, qs_kind_set, cell, iatom, descriptor_grad, forces)
      CASE DEFAULT
         CPABORT("PAO: unkown descriptor")
      END SELECT

      CALL timestop(handle)
   END SUBROUTINE pao_ml_descriptor_forces

! **************************************************************************************************
!> \brief Calulates a descriptor based on the eigenvalues of V_neighbors
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param descriptor ...
! **************************************************************************************************
   SUBROUTINE calc_descriptor_pot(particle_set, qs_kind_set, cell, iatom, descriptor)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: descriptor

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_descriptor_pot'

      INTEGER                                            :: handle, idesc, ikind, jatom, jkind, N, &
                                                            natoms, ndesc
      REAL(dp)                                           :: beta, weight
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: V_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: block_V
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

      ALLOCATE (descriptor(N*ndesc), block_V(N, N), V_evals(N))

      DO idesc = 1, ndesc

         ! construct matrix V_block from neighboring atoms
         block_V = 0.0_dp
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
            CALL pao_calc_gaussian(basis_set, block_V=block_V, Rab=Rab, lpot=0, beta=beta, weight=weight)
         ENDDO

         ! diagonalize block_V
         CALL diamat_all(block_V, V_evals)

         ! use eigenvalues of V_block as descriptor
         descriptor((idesc-1)*N+1:idesc*N) = V_evals(:)
      ENDDO

      DEALLOCATE (block_V, V_evals)

      CALL timestop(handle)
   END SUBROUTINE calc_descriptor_pot

! **************************************************************************************************
!> \brief Add forces contributed by V_neighbors descriptor
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param descriptor_grad ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE forces_descriptor_pot(particle_set, qs_kind_set, cell, iatom, descriptor_grad, forces)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: descriptor_grad
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: forces

      CHARACTER(len=*), PARAMETER :: routineN = 'forces_descriptor_pot'

      INTEGER                                            :: handle, i, idesc, ikind, jatom, jkind, &
                                                            k, N, natoms, ndesc
      REAL(dp)                                           :: beta, w, weight
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: V_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: block_M, block_V, V_evecs
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: block_D
      REAL(dp), DIMENSION(3)                             :: Ra, Rab, Rb
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(pao_descriptor_type), DIMENSION(:), POINTER   :: pao_descriptors

      CALL timeset(routineN, handle)

      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set, pao_descriptors=pao_descriptors)

      N = basis_set%nsgf
      natoms = SIZE(particle_set)
      ndesc = SIZE(pao_descriptors)
      CPASSERT(ndesc > 0)
      CPASSERT(SIZE(descriptor_grad) == ndesc*N)

      ALLOCATE (block_V(N, N), V_evecs(N, N), V_evals(N), block_D(N, N, 3), block_M(N, N))

      DO idesc = 1, ndesc
         ! construct block_V from neighboring atoms
         block_V = 0.0_dp
         DO jatom = 1, natoms
            IF (jatom == iatom) CYCLE
            Ra = particle_set(iatom)%r
            Rb = particle_set(jatom)%r
            Rab = pbc(ra, rb, cell)
            CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=jkind)
            CALL get_qs_kind(qs_kind_set(jkind), pao_descriptors=pao_descriptors)
            CPASSERT(SIZE(pao_descriptors) == ndesc)
            weight = pao_descriptors(idesc)%weight
            beta = pao_descriptors(idesc)%beta
            CALL pao_calc_gaussian(basis_set, block_V=block_V, Rab=Rab, lpot=0, beta=beta, weight=weight)
         ENDDO

         ! diagonalize block_V
         V_evecs(:, :) = block_V(:, :)
         CALL diamat_all(V_evecs, V_evals)

         ! TURNING POINT ---------------------------------------------------------------------------

         block_M = 0.0_dp
         DO k = 1, N
            w = descriptor_grad((idesc-1)*N+k)
            block_M(:, :) = block_M(:, :)+w*MATMUL(V_evecs(:, k:k), TRANSPOSE(V_evecs(:, k:k)))
         ENDDO

         DO jatom = 1, natoms
            IF (jatom == iatom) CYCLE
            Ra = particle_set(iatom)%r
            Rb = particle_set(jatom)%r
            Rab = pbc(ra, rb, cell)
            CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=jkind)
            CALL get_qs_kind(qs_kind_set(jkind), pao_descriptors=pao_descriptors)
            CPASSERT(SIZE(pao_descriptors) == ndesc)
            weight = pao_descriptors(idesc)%weight
            beta = pao_descriptors(idesc)%beta
            block_D = 0.0_dp
            CALL pao_calc_gaussian(basis_set, block_D=block_D, Rab=Rab, lpot=0, beta=beta, weight=weight)
            DO i = 1, 3
               forces(iatom, i) = forces(iatom, i)-SUM(block_M*block_D(:, :, i))
               forces(jatom, i) = forces(jatom, i)+SUM(block_M*block_D(:, :, i))
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE (block_V, V_evecs, V_evals, block_D, block_M)

      CALL timestop(handle)
   END SUBROUTINE forces_descriptor_pot

! **************************************************************************************************
!> \brief Calulates a descriptor based on distance between two atoms
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param descriptor ...
! **************************************************************************************************
   SUBROUTINE calc_descriptor_r12(particle_set, qs_kind_set, cell, iatom, descriptor)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: descriptor

      REAL(dp), DIMENSION(3)                             :: R1, R12, R2

      CPASSERT(SIZE(particle_set) == 2)

      MARK_USED(qs_kind_set)
      MARK_USED(iatom)
      MARK_USED(cell)

      R1 = particle_set(1)%r
      R2 = particle_set(2)%r
      R12 = pbc(R1, R2, cell)

      ALLOCATE (descriptor(1))
      descriptor(1) = SQRT(SUM(R12**2))

   END SUBROUTINE calc_descriptor_r12

! **************************************************************************************************
!> \brief Add forces contributed by atom-distance descriptor
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param descriptor_grad ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE forces_descriptor_r12(particle_set, qs_kind_set, cell, iatom, descriptor_grad, forces)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: descriptor_grad
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: forces

      REAL(dp), DIMENSION(3)                             :: G, R1, R12, R2

      MARK_USED(qs_kind_set)
      MARK_USED(iatom)

      R1 = particle_set(1)%r
      R2 = particle_set(2)%r
      R12 = pbc(R1, R2, cell)

      G = R12/SQRT(SUM(R12**2))*descriptor_grad(1)
      forces(1, :) = forces(1, :)+G
      forces(2, :) = forces(2, :)-G
   END SUBROUTINE forces_descriptor_r12

END MODULE pao_ml_descriptor
