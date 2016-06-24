!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Feature vectors for describing chemical environments in a rotationally invariant fashion.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_ml_descriptor
   USE ao_util,                         ONLY: exp_radius
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              rootpi
   USE mathlib,                         ONLY: diamat_all
   USE pao_input,                       ONLY: pao_ml_desc_overlap,&
                                              pao_ml_desc_pot,&
                                              pao_ml_desc_r12
   USE pao_potentials,                  ONLY: pao_calc_gaussian
   USE pao_types,                       ONLY: pao_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              pao_descriptor_type,&
                                              qs_kind_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_ml_descriptor'

   PUBLIC :: pao_ml_calc_descriptor

CONTAINS

! **************************************************************************************************
!> \brief Calulates a descriptor for chemical environment of given atom
!> \param pao ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param descriptor ...
!> \param descr_grad ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE pao_ml_calc_descriptor(pao, particle_set, qs_kind_set, cell, iatom, descriptor, descr_grad, forces)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), ALLOCATABLE, DIMENSION(:), OPTIONAL      :: descriptor
      REAL(dp), DIMENSION(:), INTENT(IN), OPTIONAL       :: descr_grad
      REAL(dp), DIMENSION(:, :), INTENT(INOUT), OPTIONAL :: forces

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_ml_calc_descriptor'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(PRESENT(forces) .EQV. PRESENT(descr_grad))

      SELECT CASE (pao%ml_descriptor)
      CASE (pao_ml_desc_pot)
         CALL calc_descriptor_pot(particle_set, qs_kind_set, cell, iatom, descriptor, descr_grad, forces)
      CASE (pao_ml_desc_overlap)
         CALL calc_descriptor_overlap(particle_set, qs_kind_set, cell, iatom, descriptor, descr_grad, forces)
      CASE (pao_ml_desc_r12)
         CALL calc_descriptor_r12(particle_set, qs_kind_set, cell, iatom, descriptor, descr_grad, forces)
      CASE DEFAULT
         CPABORT("PAO: unkown descriptor")
      END SELECT

      CALL timestop(handle)
   END SUBROUTINE pao_ml_calc_descriptor

! **************************************************************************************************
!> \brief Calulates a descriptor based on the eigenvalues of V_neighbors
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param descriptor ...
!> \param descr_grad ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE calc_descriptor_pot(particle_set, qs_kind_set, cell, iatom, descriptor, descr_grad, forces)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), ALLOCATABLE, DIMENSION(:), OPTIONAL      :: descriptor
      REAL(dp), DIMENSION(:), INTENT(IN), OPTIONAL       :: descr_grad
      REAL(dp), DIMENSION(:, :), INTENT(INOUT), OPTIONAL :: forces

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_descriptor_pot'

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
      IF (ndesc == 0) CPABORT("No PAO_DESCRIPTOR section found")

      ALLOCATE (block_V(N, N), V_evecs(N, N), V_evals(N))
      IF (PRESENT(descriptor)) ALLOCATE (descriptor(N*ndesc))
      IF (PRESENT(forces)) ALLOCATE (block_D(N, N, 3), block_M(N, N))

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
         V_evecs(:, :) = block_V(:, :)
         CALL diamat_all(V_evecs, V_evals)

         ! use eigenvalues of V_block as descriptor
         IF (PRESENT(descriptor)) &
            descriptor((idesc-1)*N+1:idesc*N) = V_evals(:)

         ! FORCES ----------------------------------------------------------------------------------
         IF (PRESENT(forces)) THEN
            CPASSERT(PRESENT(descr_grad))
            block_M = 0.0_dp
            DO k = 1, N
               w = descr_grad((idesc-1)*N+k)
               block_M(:, :) = block_M(:, :)+w*MATMUL(V_evecs(:, k:k), TRANSPOSE(V_evecs(:, k:k)))
            ENDDO
            DO jatom = 1, natoms
               IF (jatom == iatom) CYCLE
               Ra = particle_set(iatom)%r
               Rb = particle_set(jatom)%r
               Rab = pbc(ra, rb, cell)
               CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=jkind)
               CALL get_qs_kind(qs_kind_set(jkind), pao_descriptors=pao_descriptors)
               weight = pao_descriptors(idesc)%weight
               beta = pao_descriptors(idesc)%beta
               block_D = 0.0_dp
               CALL pao_calc_gaussian(basis_set, block_D=block_D, Rab=Rab, lpot=0, beta=beta, weight=weight)
               DO i = 1, 3
                  forces(iatom, i) = forces(iatom, i)-SUM(block_M*block_D(:, :, i))
                  forces(jatom, i) = forces(jatom, i)+SUM(block_M*block_D(:, :, i))
               ENDDO
            ENDDO
         ENDIF

      ENDDO

      CALL timestop(handle)
   END SUBROUTINE calc_descriptor_pot

! **************************************************************************************************
!> \brief Calulates a descriptor based on the eigenvalues of local overlap matrix
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param descriptor ...
!> \param descr_grad ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE calc_descriptor_overlap(particle_set, qs_kind_set, cell, iatom, descriptor, descr_grad, forces)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), ALLOCATABLE, DIMENSION(:), OPTIONAL      :: descriptor
      REAL(dp), DIMENSION(:), INTENT(IN), OPTIONAL       :: descr_grad
      REAL(dp), DIMENSION(:, :), INTENT(INOUT), OPTIONAL :: forces

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_descriptor_overlap'

      INTEGER                                            :: handle, idesc, ikind, j, jatom, jkind, &
                                                            k, katom, kkind, N, natoms, ndesc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: neighbor_order
      REAL(dp) :: beta_sum, deriv, exponent, integral, jbeta, jweight, kbeta, kweight, &
         normalization, Rij2, Rik2, Rjk2, sbeta, sbeta_min, screening_radius, screening_volume, w
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: S_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: block_M, block_S, S_evecs
      REAL(dp), DIMENSION(3)                             :: Ri, Rij, Rik, Rj, Rjk, Rk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: neighbor_dist
      TYPE(pao_descriptor_type), DIMENSION(:), POINTER   :: ipao_descriptors, jpao_descriptors, &
                                                            kpao_descriptors

      CALL timeset(routineN, handle)

      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), pao_descriptors=ipao_descriptors)

      natoms = SIZE(particle_set)
      ndesc = SIZE(ipao_descriptors)
      IF (ndesc == 0) CPABORT("No PAO_DESCRIPTOR section found")

      ! determine smallest screening beta
      sbeta_min = HUGE(1.0_dp)
      DO idesc = 1, ndesc
         sbeta_min = MIN(sbeta_min, ipao_descriptors(idesc)%screening)
      ENDDO

      ! estimate maximum number of neighbors within screening
      screening_radius = exp_radius(0, sbeta_min, 1.0E-12_dp, 1.0_dp) ! TODO get the EPS parameter from somewhere / precompute this elsewhere
      screening_volume = fourpi/3.0_dp*screening_radius**3
      N = INT(screening_volume/45.0_dp) ! rule of thumb

      ALLOCATE (block_S(N, N), S_evals(N), S_evecs(N, N))
      IF (PRESENT(descriptor)) ALLOCATE (descriptor(N*ndesc))
      IF (PRESENT(forces)) ALLOCATE (block_M(N, N))

      !find neighbors
      !TODO: this is a quadratic algorithm, use a neighbor-list instead
      ALLOCATE (neighbor_dist(natoms), neighbor_order(natoms))
      Ri = particle_set(iatom)%r
      DO jatom = 1, natoms
         Rj = particle_set(jatom)%r
         Rij = pbc(Ri, Rj, cell)
         neighbor_dist(jatom) = SQRT(SUM(Rij**2))
      ENDDO
      CALL sort(neighbor_dist, natoms, neighbor_order)
      CPASSERT(neighbor_order(1) == iatom) !central atom should be closesd to itself

      ! check if N was choosen large enough
      IF (natoms > N) THEN
         IF (neighbor_dist(N+1) < screening_radius) &
            CPABORT("PAO heuristic for descriptor size broke down")
      ENDIF

      DO idesc = 1, ndesc
         sbeta = ipao_descriptors(idesc)%screening

         ! construct matrix S_block from neighboring atoms
         block_S = 0.0_dp
         DO j = 1, MIN(natoms, N)
         DO k = 1, MIN(natoms, N)
            jatom = neighbor_order(j)
            katom = neighbor_order(k)

            ! get weigths and betas
            CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=jkind)
            CALL get_qs_kind(qs_kind_set(jkind), pao_descriptors=jpao_descriptors)
            CALL get_atomic_kind(particle_set(katom)%atomic_kind, kind_number=kkind)
            CALL get_qs_kind(qs_kind_set(kkind), pao_descriptors=kpao_descriptors)
            IF (SIZE(jpao_descriptors) /= ndesc .OR. SIZE(kpao_descriptors) /= ndesc) &
               CPABORT("Not all KINDs have the same number of PAO_DESCRIPTOR sections")
            jweight = jpao_descriptors(idesc)%weight
            jbeta = jpao_descriptors(idesc)%beta
            kweight = kpao_descriptors(idesc)%weight
            kbeta = kpao_descriptors(idesc)%beta
            beta_sum = sbeta+jbeta+kbeta

            ! get distances
            Rj = particle_set(jatom)%r
            Rk = particle_set(katom)%r
            Rij = pbc(Ri, Rj, cell)
            Rik = pbc(Ri, Rk, cell)
            Rjk = pbc(Rj, Rk, cell)
            Rij2 = SUM(Rij**2)
            Rik2 = SUM(Rik**2)
            Rjk2 = SUM(Rjk**2)

            ! calculate integral over three Gaussians
            exponent = -(sbeta*jbeta*Rij2+sbeta*kbeta*Rik2+jbeta*kbeta*Rjk2)/beta_sum
            integral = EXP(exponent)*rootpi/SQRT(beta_sum)
            normalization = SQRT(jbeta*kbeta)/rootpi**2
            block_S(j, k) = jweight*kweight*normalization*integral
         ENDDO
         ENDDO

         ! diagonalize V_block
         S_evecs(:, :) = block_S(:, :)
         CALL diamat_all(S_evecs, S_evals)

         ! use eigenvalues of S_block as descriptor
         IF (PRESENT(descriptor)) &
            descriptor((idesc-1)*N+1:idesc*N) = S_evals(:)

         ! FORCES ----------------------------------------------------------------------------------
         IF (PRESENT(forces)) THEN
            CPASSERT(PRESENT(descr_grad))
            block_M = 0.0_dp
            DO k = 1, N
               w = descr_grad((idesc-1)*N+k)
               block_M(:, :) = block_M(:, :)+w*MATMUL(S_evecs(:, k:k), TRANSPOSE(S_evecs(:, k:k)))
            ENDDO

            DO j = 1, MIN(natoms, N)
            DO k = 1, MIN(natoms, N)
               jatom = neighbor_order(j)
               katom = neighbor_order(k)

               ! get weigths and betas
               CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=jkind)
               CALL get_qs_kind(qs_kind_set(jkind), pao_descriptors=jpao_descriptors)
               CALL get_atomic_kind(particle_set(katom)%atomic_kind, kind_number=kkind)
               CALL get_qs_kind(qs_kind_set(kkind), pao_descriptors=kpao_descriptors)
               jweight = jpao_descriptors(idesc)%weight
               jbeta = jpao_descriptors(idesc)%beta
               kweight = kpao_descriptors(idesc)%weight
               kbeta = kpao_descriptors(idesc)%beta
               beta_sum = sbeta+jbeta+kbeta

               ! get distances
               Rj = particle_set(jatom)%r
               Rk = particle_set(katom)%r
               Rij = pbc(Ri, Rj, cell)
               Rik = pbc(Ri, Rk, cell)
               Rjk = pbc(Rj, Rk, cell)
               Rij2 = SUM(Rij**2)
               Rik2 = SUM(Rik**2)
               Rjk2 = SUM(Rjk**2)

               ! calculate integral over three Gaussians
               exponent = -(sbeta*jbeta*Rij2+sbeta*kbeta*Rik2+jbeta*kbeta*Rjk2)/beta_sum
               integral = EXP(exponent)*rootpi/SQRT(beta_sum)
               normalization = SQRT(jbeta*kbeta)/rootpi**2
               deriv = 2.0_dp/beta_sum*block_M(j, k)
               w = jweight*kweight*normalization*integral*deriv
               forces(iatom, :) = forces(iatom, :)-sbeta*jbeta*Rij*w
               forces(jatom, :) = forces(jatom, :)+sbeta*jbeta*Rij*w
               forces(iatom, :) = forces(iatom, :)-sbeta*kbeta*Rik*w
               forces(katom, :) = forces(katom, :)+sbeta*kbeta*Rik*w
               forces(jatom, :) = forces(jatom, :)-jbeta*kbeta*Rjk*w
               forces(katom, :) = forces(katom, :)+jbeta*kbeta*Rjk*w
            ENDDO
            ENDDO
         ENDIF
      ENDDO

      CALL timestop(handle)
   END SUBROUTINE calc_descriptor_overlap

! **************************************************************************************************
!> \brief Calulates a descriptor based on distance between two atoms
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param descriptor ...
!> \param descr_grad ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE calc_descriptor_r12(particle_set, qs_kind_set, cell, iatom, descriptor, descr_grad, forces)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), ALLOCATABLE, DIMENSION(:), OPTIONAL      :: descriptor
      REAL(dp), DIMENSION(:), INTENT(IN), OPTIONAL       :: descr_grad
      REAL(dp), DIMENSION(:, :), INTENT(INOUT), OPTIONAL :: forces

      REAL(dp), DIMENSION(3)                             :: G, R1, R12, R2

      CPASSERT(SIZE(particle_set) == 2)

      MARK_USED(qs_kind_set)
      MARK_USED(iatom)
      MARK_USED(cell)

      R1 = particle_set(1)%r
      R2 = particle_set(2)%r
      R12 = pbc(R1, R2, cell)

      IF (PRESENT(descriptor)) THEN
         ALLOCATE (descriptor(1))
         descriptor(1) = SQRT(SUM(R12**2))
      ENDIF

      IF (PRESENT(forces)) THEN
         CPASSERT(PRESENT(descr_grad))
         G = R12/SQRT(SUM(R12**2))*descr_grad(1)
         forces(1, :) = forces(1, :)+G
         forces(2, :) = forces(2, :)-G
      ENDIF
   END SUBROUTINE calc_descriptor_r12

END MODULE pao_ml_descriptor
