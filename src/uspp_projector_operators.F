!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Distributed real scalar-USPP operators in the Gaussian and projector bases.
!>        Each column block of B=<chi|beta> belongs to one augmentation center.
!>        Reciprocal vectors are distributed without replication over B's communicator.
! **************************************************************************************************
MODULE uspp_projector_operators
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_block_p, &
        dbcsr_get_info, dbcsr_get_stored_coordinates, dbcsr_multiply, dbcsr_release, &
        dbcsr_reserve_blocks, dbcsr_set, dbcsr_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_collocate,&
                                              uspp_augmentation_integrate,&
                                              uspp_augmentation_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'uspp_projector_operators'

   PUBLIC :: uspp_projector_density, uspp_projector_add_operator, uspp_projector_metric, &
             uspp_projector_collocate, uspp_projector_integrate, uspp_projector_pullback, uspp_projector_ao_pattern
   PUBLIC :: uspp_projector_create_onsite

CONTAINS

! **************************************************************************************************
!> \brief Union of Gaussian blocks and all AO block pairs sharing an augmentation center.
!>        A zero integrated q_ij does not imply a zero Q_ij(G) or screened operator.
!>        Therefore the pattern must depend on projector connectivity, not metric values.
!> \param overlaps B, with one block column per augmentation center.
!> \param gaussian_overlap Gaussian overlap defining the ordinary AO blocks.
!> \param pattern New AO matrix with unit entries on the required block pattern.
! **************************************************************************************************
   SUBROUTINE uspp_projector_ao_pattern(overlaps, gaussian_overlap, pattern)
      TYPE(dbcsr_type), INTENT(IN)                       :: overlaps, gaussian_overlap
      TYPE(dbcsr_type), INTENT(INOUT)                    :: pattern

      TYPE(dbcsr_type)                                   :: connectivity

      CALL dbcsr_copy(connectivity, overlaps)
      CALL dbcsr_set(connectivity, 1.0_dp)
      CALL dbcsr_copy(pattern, gaussian_overlap)
      CALL dbcsr_set(pattern, 1.0_dp)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, connectivity, connectivity, 1.0_dp, pattern)
      CALL dbcsr_set(pattern, 1.0_dp)
      CALL dbcsr_release(connectivity)
   END SUBROUTINE uspp_projector_ao_pattern

! **************************************************************************************************
!> \brief Create an empty block-diagonal projector matrix on B's processor grid.
!> \param overlaps Rectangular Gaussian/projector matrix B.
!> \param onsite New matrix with exactly one diagonal block per center.
! **************************************************************************************************
   SUBROUTINE uspp_projector_create_onsite(overlaps, onsite)
      TYPE(dbcsr_type), INTENT(IN)                       :: overlaps
      TYPE(dbcsr_type), INTENT(INOUT)                    :: onsite

      INTEGER                                            :: i, nprows, owner
      INTEGER, POINTER                                   :: col_dist(:), row_dist(:), sizes(:)
      TYPE(dbcsr_distribution_type)                      :: dist, onsite_dist
      TYPE(mp_comm_type)                                 :: group

      CALL dbcsr_get_info(overlaps, distribution=dist, col_blk_size=sizes, group=group)
      CALL dbcsr_distribution_get(dist, nprows=nprows, col_dist=col_dist)
      ALLOCATE (row_dist(SIZE(sizes)))
      DO i = 1, SIZE(sizes)
         row_dist(i) = MODULO(i - 1, nprows)
      END DO
      CALL dbcsr_distribution_new(onsite_dist, template=dist, row_dist=row_dist, col_dist=col_dist)
      ! N stores full on-site blocks and also supports rectangular MPI processor grids.
      CALL dbcsr_create(onsite, 'USPP onsite', onsite_dist, 'N', sizes, sizes)
      CALL dbcsr_distribution_release(onsite_dist)
      DEALLOCATE (row_dist)
      DO i = 1, SIZE(sizes)
         CALL dbcsr_get_stored_coordinates(onsite, i, i, owner)
         IF (owner == group%mepos) CALL dbcsr_reserve_blocks(onsite, rows=[i], cols=[i])
      END DO
      CALL dbcsr_finalize(onsite)
      CALL dbcsr_set(onsite, 0.0_dp)
   END SUBROUTINE uspp_projector_create_onsite

! **************************************************************************************************
!> \brief Compute only the atomic diagonal blocks of B^T P B.
!> \param overlaps B, with one column block per augmentation center.
!> \param density Gaussian AO density matrix P, including spin occupations.
!> \param projected New block-diagonal projector density matrix.
! **************************************************************************************************
   SUBROUTINE uspp_projector_density(overlaps, density, projected)
      TYPE(dbcsr_type), INTENT(IN)                       :: overlaps, density
      TYPE(dbcsr_type), INTENT(INOUT)                    :: projected

      TYPE(dbcsr_type)                                   :: work

      CALL uspp_projector_create_onsite(overlaps, projected)
      CALL dbcsr_create(work, template=overlaps)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, density, overlaps, 0.0_dp, work)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, overlaps, work, 0.0_dp, projected, retain_sparsity=.TRUE.)
      CALL dbcsr_release(work)
   END SUBROUTINE uspp_projector_density

! **************************************************************************************************
!> \brief Add B A B^T to an AO operator, allowing new AO blocks.
!> \param overlaps B.
!> \param onsite Block-diagonal A (overlap q or screened nonlocal potential).
!> \param ao_operator Existing AO matrix, updated without imposing its old sparsity pattern.
! **************************************************************************************************
   SUBROUTINE uspp_projector_add_operator(overlaps, onsite, ao_operator)
      TYPE(dbcsr_type), INTENT(IN)                       :: overlaps, onsite
      TYPE(dbcsr_type), INTENT(INOUT)                    :: ao_operator

      TYPE(dbcsr_type)                                   :: work

      CALL dbcsr_create(work, template=overlaps)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, overlaps, onsite, 0.0_dp, work)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, work, overlaps, 1.0_dp, ao_operator)
      CALL dbcsr_release(work)
   END SUBROUTINE uspp_projector_add_operator

! **************************************************************************************************
!> \brief Construct S_G+B q B^T while preserving the Gaussian overlap S_G.
!> \param overlaps B.
!> \param augmentation Prepared radial/angular data per kind.
!> \param kind Kind index for each column block of B.
!> \param gaussian_overlap Ordinary Gaussian overlap, unchanged.
!> \param metric New generalized overlap matrix.
! **************************************************************************************************
   SUBROUTINE uspp_projector_metric(overlaps, augmentation, kind, gaussian_overlap, metric)
      TYPE(dbcsr_type), INTENT(IN)                       :: overlaps
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation(:)
      INTEGER, INTENT(IN)                                :: kind(:)
      TYPE(dbcsr_type), INTENT(IN)                       :: gaussian_overlap
      TYPE(dbcsr_type), INTENT(INOUT)                    :: metric

      INTEGER                                            :: i
      LOGICAL                                            :: found
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_type)                                   :: onsite

      CALL check_centers(overlaps, augmentation, kind)
      CALL uspp_projector_create_onsite(overlaps, onsite)
      DO i = 1, SIZE(kind)
         CALL dbcsr_get_block_p(onsite, i, i, block, found)
         IF (found) block(:, :) = augmentation(KIND(i))%q
      END DO
      CALL dbcsr_copy(metric, gaussian_overlap)
      CALL uspp_projector_add_operator(overlaps, onsite, metric)
      CALL dbcsr_release(onsite)
   END SUBROUTINE uspp_projector_metric

! **************************************************************************************************
!> \brief Collocate the distributed AO density's augmentation on a local reciprocal list.
!> \param overlaps B.
!> \param density P for one spin, including occupations.
!> \param augmentation Prepared data per kind.
!> \param kind Kind index for each projector center.
!> \param position Center positions (3,ncenter), in bohr.
!> \param gvecs This rank's reciprocal vectors (3,nglocal), in inverse bohr.
!> \param volume Cell volume.
!> \param rho Local Fourier coefficients, integral/volume, accumulated in place.
!> \param charge Total augmentation charge, replicated on every rank.
! **************************************************************************************************
   SUBROUTINE uspp_projector_collocate(overlaps, density, augmentation, kind, position, gvecs, volume, rho, charge)
      TYPE(dbcsr_type), INTENT(IN)                       :: overlaps, density
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation(:)
      INTEGER, INTENT(IN)                                :: kind(:)
      REAL(KIND=dp), INTENT(IN)                          :: position(:, :), gvecs(:, :), volume
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: rho(:)
      REAL(KIND=dp), INTENT(OUT)                         :: charge

      INTEGER                                            :: i, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE                         :: matrix(:, :)
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_type)                                   :: projected
      TYPE(mp_comm_type)                                 :: group

      CALL check_centers(overlaps, augmentation, kind)
      CPASSERT(ALL(SHAPE(position) == [3, SIZE(kind)]))
      CALL dbcsr_get_info(overlaps, group=group)
      CALL uspp_projector_density(overlaps, density, projected)
      charge = 0.0_dp
      DO i = 1, SIZE(kind)
         n = augmentation(KIND(i))%nproj
         ALLOCATE (matrix(n, n))
         matrix = 0.0_dp
         CALL dbcsr_get_block_p(projected, i, i, block, found)
         IF (found) matrix(:, :) = block
         ! Only the owner contributes; each grid rank needs this atom's density.
         CALL group%sum(matrix)
         CALL uspp_augmentation_collocate(augmentation(KIND(i)), gvecs, position(:, i), volume, matrix, rho)
         charge = charge + SUM(matrix*augmentation(KIND(i))%q)
         DEALLOCATE (matrix)
      END DO
      CALL dbcsr_release(projected)
   END SUBROUTINE uspp_projector_collocate

! **************************************************************************************************
!> \brief Integrate a distributed local potential and add B A[V] B^T to the Hamiltonian.
!> \param overlaps B.
!> \param augmentation Prepared data per kind.
!> \param kind Kind index for each projector center.
!> \param position Center positions (3,ncenter).
!> \param gvecs Local reciprocal vectors; every retained G occurs on exactly one rank.
!> \param volume Cell volume.
!> \param vg Local Fourier potential coefficients, integral/volume (no grid-volume weight).
!> \param multiplicity Conjugate-pair weights, one for self-conjugate coefficients.
!> \param hamiltonian AO matrix, incremented by the screened augmentation operator.
!> \param density Optional AO density for the derivatives; requires both gradient outputs.
!> \param gradient Optional new B-shaped derivative matrix, distributed like overlaps.
!> \param center_gradient Optional local center derivatives (3,ncenter), not MPI-reduced.
!> \param cell_gradient Optional local density strain contribution, incremented without MPI reduction.
! **************************************************************************************************
   SUBROUTINE uspp_projector_integrate(overlaps, augmentation, kind, position, gvecs, volume, vg, multiplicity, hamiltonian, &
                                       density, gradient, center_gradient, cell_gradient)
      TYPE(dbcsr_type), INTENT(IN)                       :: overlaps
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation(:)
      INTEGER, INTENT(IN)                                :: kind(:)
      REAL(KIND=dp), INTENT(IN)                          :: position(:, :), gvecs(:, :), volume
      COMPLEX(KIND=dp), INTENT(IN)                       :: vg(:)
      REAL(KIND=dp), INTENT(IN)                          :: multiplicity(:)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: hamiltonian
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: density
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: gradient
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: center_gradient(:, :)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: cell_gradient(3, 3)

      INTEGER                                            :: i, idir, jdir, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE                         :: derivative(:, :, :), dm(:, :), &
                                                            matrix(:, :), strain(:, :, :, :)
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_type)                                   :: onsite, projected
      TYPE(mp_comm_type)                                 :: group

      CALL check_centers(overlaps, augmentation, kind)
      CPASSERT(ALL(SHAPE(position) == [3, SIZE(kind)]))
      CALL dbcsr_get_info(overlaps, group=group)
      CALL uspp_projector_create_onsite(overlaps, onsite)
      CPASSERT(PRESENT(density) .EQV. PRESENT(gradient))
      CPASSERT(PRESENT(density) .EQV. PRESENT(center_gradient))
      IF (PRESENT(cell_gradient)) THEN
         CPASSERT(PRESENT(density))
      END IF
      IF (PRESENT(density)) THEN
         CPASSERT(ALL(SHAPE(center_gradient) == [3, SIZE(kind)]))
         CALL uspp_projector_density(overlaps, density, projected)
         center_gradient = 0.0_dp
      END IF
      DO i = 1, SIZE(kind)
         n = augmentation(KIND(i))%nproj
         ALLOCATE (matrix(n, n))
         IF (PRESENT(density)) THEN
            ALLOCATE (dm(n, n), derivative(n, n, 3))
            IF (PRESENT(cell_gradient)) ALLOCATE (strain(n, n, 3, 3))
            dm = 0.0_dp
            CALL dbcsr_get_block_p(projected, i, i, block, found)
            IF (found) dm(:, :) = block
            CALL group%sum(dm)
            CALL uspp_augmentation_integrate(augmentation(KIND(i)), gvecs, position(:, i), volume, vg, multiplicity, &
                                             matrix, derivative, strain)
            DO idir = 1, 3
               center_gradient(idir, i) = SUM(dm*derivative(:, :, idir))
            END DO
            IF (PRESENT(cell_gradient)) THEN
               DO jdir = 1, 3
                  DO idir = 1, 3
                     cell_gradient(idir, jdir) = cell_gradient(idir, jdir) + SUM(dm*strain(:, :, idir, jdir))
                  END DO
               END DO
               DEALLOCATE (strain)
            END IF
            DEALLOCATE (dm, derivative)
         ELSE
            CALL uspp_augmentation_integrate(augmentation(KIND(i)), gvecs, position(:, i), volume, vg, multiplicity, matrix)
         END IF
         CALL group%sum(matrix)
         CALL dbcsr_get_block_p(onsite, i, i, block, found)
         IF (found) block(:, :) = matrix
         DEALLOCATE (matrix)
      END DO
      CALL uspp_projector_add_operator(overlaps, onsite, hamiltonian)
      IF (PRESENT(density)) THEN
         CALL uspp_projector_pullback(overlaps, density, onsite, 1.0_dp, gradient)
         CALL dbcsr_release(projected)
      END IF
      CALL dbcsr_release(onsite)
   END SUBROUTINE uspp_projector_integrate

! **************************************************************************************************
!> \brief Derivative of alpha*Tr(P B A B^T) with respect to B, at fixed P and A.
!>        Applying this to P,A=D and W,A=q with opposite signs gives the projector
!>        part of the variational force, including the generalized-metric Pulay term.
!> \param overlaps B.
!> \param density Symmetric P (or energy-weighted W).
!> \param onsite Symmetric on-site A.
!> \param alpha Prefactor, e.g. minus one for the metric term.
!> \param gradient New B-shaped matrix containing 2*alpha*P B A.
! **************************************************************************************************
   SUBROUTINE uspp_projector_pullback(overlaps, density, onsite, alpha, gradient)
      TYPE(dbcsr_type), INTENT(IN)                       :: overlaps, density, onsite
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(dbcsr_type), INTENT(INOUT)                    :: gradient

      TYPE(dbcsr_type)                                   :: work

      CALL dbcsr_create(work, template=overlaps)
      CALL dbcsr_create(gradient, template=overlaps)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, density, overlaps, 0.0_dp, work)
      CALL dbcsr_multiply('N', 'N', 2.0_dp*alpha, work, onsite, 0.0_dp, gradient)
      CALL dbcsr_release(work)
   END SUBROUTINE uspp_projector_pullback

! **************************************************************************************************
!> \brief Validate the mapping from projector column blocks to augmentation kinds.
!> \param overlaps B.
!> \param augmentation Prepared data per kind.
!> \param kind Kind index for each column block of B.
! **************************************************************************************************
   SUBROUTINE check_centers(overlaps, augmentation, kind)
      TYPE(dbcsr_type), INTENT(IN)                       :: overlaps
      TYPE(uspp_augmentation_type), INTENT(IN)           :: augmentation(:)
      INTEGER, INTENT(IN)                                :: kind(:)

      INTEGER                                            :: i
      INTEGER, POINTER                                   :: sizes(:)

      CALL dbcsr_get_info(overlaps, col_blk_size=sizes)
      CPASSERT(SIZE(sizes) == SIZE(kind))
      CPASSERT(ALL(kind >= 1) .AND. ALL(kind <= SIZE(augmentation)))
      DO i = 1, SIZE(kind)
         CPASSERT(sizes(i) == augmentation(KIND(i))%nproj)
      END DO
   END SUBROUTINE check_centers

END MODULE uspp_projector_operators
