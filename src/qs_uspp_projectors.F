!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief USPP projector assembly and image-resolved contractions from Quickstep overlap lists.
! **************************************************************************************************
MODULE qs_uspp_projectors
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_get_readonly_block_p, dbcsr_get_stored_coordinates, dbcsr_p_type, &
        dbcsr_reserve_blocks, dbcsr_set, dbcsr_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE sap_kind_types,                  ONLY: alist_type,&
                                              get_alist,&
                                              sap_int_type
   USE uspp_projector_operators,        ONLY: uspp_projector_create_onsite
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_uspp_projectors'
   PUBLIC :: qs_uspp_projector_matrix, qs_uspp_projector_gradient
   PUBLIC :: qs_uspp_projector_pair, qs_uspp_projector_images

CONTAINS

! **************************************************************************************************
!> \brief Contract one AO pair and one augmentation center without summing away lattice images.
!>        For the AO pair (a,0),(b,R), only projector images T,U with T-U=R contribute:
!>        H_ab(R) = sum_T B_ac(T) A_c B_bc(T-R)^T. The optional reverse contraction
!>        accumulates the real symmetric on-site density, the adjoint for scalar USPP.
!>        Atom blocks use their natural a,b orientation, even when a>b; the caller handles
!>        triangular storage, neighbor-list multiplicities and MPI reduction. No k weight
!>        is applied here: real-space density blocks must already contain the weighted k sum.
!> \param alist_a Uncompressed overlaps and AO-center derivatives for atom a.
!> \param alist_b Uncompressed overlaps and AO-center derivatives for atom b.
!> \param center_atom Global atom index of the selected non-ghost augmentation center.
!> \param cell_b Lattice translation R of atom b relative to atom a.
!> \param onsite Symmetric on-site operator A_c, required for operator and energy derivatives.
!> \param ao_operator AO block H_ab(R), incremented in natural atom order.
!> \param density Real-space AO density block in the same natural order.
!> \param projected Real symmetric on-site density, incremented as the adjoint of ao_operator.
!> \param weight AO-pair multiplicity for density and energy derivatives (default one).
!> \param atom_gradient Local atomic energy gradients (3,natom), incremented without reduction.
!> \param cell_gradient Local strain energy gradient at fixed reduced Bloch vectors.
! **************************************************************************************************
   SUBROUTINE qs_uspp_projector_pair(alist_a, alist_b, center_atom, cell_b, onsite, ao_operator, density, &
                                     projected, weight, atom_gradient, cell_gradient)
      TYPE(alist_type), INTENT(IN)                       :: alist_a, alist_b
      INTEGER, INTENT(IN)                                :: center_atom, cell_b(3)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: onsite(:, :)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: ao_operator(:, :)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: density(:, :)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: projected(:, :)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: weight
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: atom_gradient(:, :), cell_gradient(3, 3)

      INTEGER                                            :: a, b, idir, na, nb, np
      LOGICAL                                            :: derivatives
      REAL(KIND=dp)                                      :: factor, ga, gb
      REAL(KIND=dp), ALLOCATABLE                         :: local_density(:, :), weighted_a(:, :)

      CPASSERT(center_atom > 0)
      derivatives = PRESENT(atom_gradient) .OR. PRESENT(cell_gradient)
      IF (PRESENT(ao_operator) .OR. derivatives) THEN
         CPASSERT(PRESENT(onsite))
      END IF
      IF (PRESENT(projected) .OR. derivatives) THEN
         CPASSERT(PRESENT(density))
      END IF
      IF (.NOT. ASSOCIATED(alist_a%clist) .OR. .NOT. ASSOCIATED(alist_b%clist)) RETURN
      factor = 1.0_dp
      IF (PRESENT(weight)) factor = weight
      DO a = 1, SIZE(alist_a%clist)
         IF (alist_a%clist(a)%catom /= center_atom) CYCLE
         DO b = 1, SIZE(alist_b%clist)
            IF (alist_b%clist(b)%catom /= center_atom) CYCLE
            IF (ANY(cell_b + alist_b%clist(b)%cell - alist_a%clist(a)%cell /= 0)) CYCLE
            CPASSERT(ASSOCIATED(alist_a%clist(a)%acint))
            CPASSERT(ASSOCIATED(alist_b%clist(b)%acint))
            ASSOCIATE (ba => alist_a%clist(a)%acint, bb => alist_b%clist(b)%acint)
               na = SIZE(ba, 1)
               nb = SIZE(bb, 1)
               np = SIZE(ba, 2)
               CPASSERT(SIZE(bb, 2) == np .AND. SIZE(ba, 3) >= 1 .AND. SIZE(bb, 3) >= 1)
               IF (PRESENT(density)) THEN
                  CPASSERT(ALL(SHAPE(density) == [na, nb]))
               END IF
               IF (PRESENT(projected)) THEN
                  CPASSERT(ALL(SHAPE(projected) == [np, np]))
                  ALLOCATE (local_density(np, np))
                  local_density(:, :) = MATMUL(TRANSPOSE(ba(:, :, 1)), MATMUL(density, bb(:, :, 1)))
                  projected = projected + 0.5_dp*factor*(local_density + TRANSPOSE(local_density))
                  DEALLOCATE (local_density)
               END IF
               IF (PRESENT(ao_operator) .OR. derivatives) THEN
                  CPASSERT(ALL(SHAPE(onsite) == [np, np]))
                  ALLOCATE (weighted_a(na, np))
                  weighted_a(:, :) = MATMUL(ba(:, :, 1), onsite)
                  IF (PRESENT(ao_operator)) THEN
                     CPASSERT(ALL(SHAPE(ao_operator) == [na, nb]))
                     ao_operator = ao_operator + MATMUL(weighted_a, TRANSPOSE(bb(:, :, 1)))
                  END IF
                  IF (derivatives) THEN
                     CPASSERT(SIZE(ba, 3) >= 4 .AND. SIZE(bb, 3) >= 4)
                     DO idir = 1, 3
                        ga = factor*SUM(density*MATMUL(MATMUL(ba(:, :, idir + 1), onsite), TRANSPOSE(bb(:, :, 1))))
                        gb = factor*SUM(density*MATMUL(weighted_a, TRANSPOSE(bb(:, :, idir + 1))))
                        IF (PRESENT(atom_gradient)) THEN
                           CPASSERT(SIZE(atom_gradient, 1) == 3)
                           CPASSERT(MIN(alist_a%aatom, alist_b%aatom) > 0)
                           CPASSERT(alist_a%aatom <= SIZE(atom_gradient, 2))
                           CPASSERT(alist_b%aatom <= SIZE(atom_gradient, 2))
                           CPASSERT(center_atom <= SIZE(atom_gradient, 2))
                           atom_gradient(idir, alist_a%aatom) = atom_gradient(idir, alist_a%aatom) + ga
                           atom_gradient(idir, alist_b%aatom) = atom_gradient(idir, alist_b%aatom) + gb
                           atom_gradient(idir, center_atom) = atom_gradient(idir, center_atom) - ga - gb
                        END IF
                        IF (PRESENT(cell_gradient)) THEN
                           ! At fixed reduced k, the image phases do not change with strain.
                           cell_gradient(idir, :) = cell_gradient(idir, :) &
                                                    - ga*alist_a%clist(a)%rac - gb*alist_b%clist(b)%rac
                        END IF
                     END DO
                  END IF
                  DEALLOCATE (weighted_a)
               END IF
            END ASSOCIATE
         END DO
      END DO
   END SUBROUTINE qs_uspp_projector_pair

! **************************************************************************************************
!> \brief Distributed image-pair operator, adjoint density and energy gradients.
!>        Pair columns contain a,b,kind(a),kind(b),image,R(1:3). The list and AO matrices
!>        use the same distribution. Each stored AO block is contracted on its owner only.
!> \param template Rectangular projector matrix defining center blocks and communicator.
!> \param sap Sorted image-resolved overlaps, with derivatives when gradients are requested.
!> \param pairs Local real-space AO neighbor pairs, including the cell index.
!> \param nkind Number of kinds used to index sap.
!> \param center_atoms Global atoms corresponding to projector columns.
!> \param center_kinds Kinds corresponding to projector columns.
!> \param onsite Optional distributed real on-site operator.
!> \param ao_operator Optional AO image matrices, incremented in place.
!> \param density Optional AO image densities, including spin occupations and k weights.
!> \param projected Optional new block-diagonal projected density, globally reduced once.
!> \param atom_gradient Optional local atomic energy gradients, not MPI-reduced.
!> \param cell_gradient Optional local strain energy gradient, not MPI-reduced.
!> \param alpha Optional prefactor for operator and its gradients, e.g. -1 for the metric.
! **************************************************************************************************
   SUBROUTINE qs_uspp_projector_images(template, sap, pairs, nkind, center_atoms, center_kinds, &
                                       onsite, ao_operator, density, projected, atom_gradient, cell_gradient, alpha)
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(sap_int_type), INTENT(IN)                     :: sap(:)
      INTEGER, INTENT(IN)                                :: pairs(:, :), nkind, center_atoms(:), &
                                                            center_kinds(:)
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: onsite
      TYPE(dbcsr_p_type), INTENT(INOUT), OPTIONAL        :: ao_operator(:)
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: density(:)
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: projected
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: atom_gradient(:, :), cell_gradient(3, 3)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha

      INTEGER                                            :: center, col, ENTRY, ia, iac, ib, ibc, &
                                                            image, n, owner, row
      INTEGER, POINTER                                   :: row_sizes(:), sizes(:)
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: factor, scale
      REAL(KIND=dp), ALLOCATABLE                         :: dm(:, :), matrix(:, :), pd(:, :), &
                                                            ph(:, :)
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(alist_type), POINTER                          :: alist_a, alist_b
      TYPE(mp_comm_type)                                 :: group

      CPASSERT(SIZE(pairs, 1) == 8)
      CPASSERT(SIZE(sap) == nkind*nkind)
      CALL dbcsr_get_info(template, group=group, row_blk_size=row_sizes, col_blk_size=sizes)
      CPASSERT(SIZE(center_atoms) == SIZE(sizes))
      CPASSERT(SIZE(center_kinds) == SIZE(sizes))
      scale = 1.0_dp
      IF (PRESENT(alpha)) scale = alpha
      IF (PRESENT(projected)) CALL uspp_projector_create_onsite(template, projected)
      DO center = 1, SIZE(center_atoms)
         n = sizes(center)
         IF (PRESENT(onsite)) THEN
            ALLOCATE (matrix(n, n))
            matrix = 0.0_dp
            CALL dbcsr_get_readonly_block_p(onsite, center, center, block=block, found=found)
            IF (found) matrix(:, :) = scale*block
            CALL group%sum(matrix)
         END IF
         IF (PRESENT(projected)) THEN
            ALLOCATE (dm(n, n))
            dm = 0.0_dp
         END IF
         DO ENTRY = 1, SIZE(pairs, 2)
            ia = pairs(1, ENTRY)
            ib = pairs(2, ENTRY)
            image = pairs(5, ENTRY)
            row = MIN(ia, ib)
            col = MAX(ia, ib)
            IF (PRESENT(ao_operator)) THEN
               CALL dbcsr_get_stored_coordinates(ao_operator(image)%matrix, row, col, owner)
            ELSE
               CPASSERT(PRESENT(density))
               CALL dbcsr_get_stored_coordinates(density(image)%matrix, row, col, owner)
            END IF
            IF (owner /= group%mepos) CYCLE
            iac = pairs(3, ENTRY) + nkind*(center_kinds(center) - 1)
            ibc = pairs(4, ENTRY) + nkind*(center_kinds(center) - 1)
            IF (.NOT. ASSOCIATED(sap(iac)%alist) .OR. .NOT. ASSOCIATED(sap(ibc)%alist)) CYCLE
            CALL get_alist(sap(iac), alist_a, ia)
            CALL get_alist(sap(ibc), alist_b, ib)
            IF (.NOT. ASSOCIATED(alist_a) .OR. .NOT. ASSOCIATED(alist_b)) CYCLE
            IF (PRESENT(ao_operator)) THEN
               ALLOCATE (ph(row_sizes(ia), row_sizes(ib)))
               ph = 0.0_dp
            END IF
            IF (PRESENT(density)) THEN
               ALLOCATE (pd(row_sizes(ia), row_sizes(ib)))
               pd = 0.0_dp
               CALL dbcsr_get_block_p(density(image)%matrix, row, col, block, found)
               IF (found) THEN
                  IF (ia <= ib) THEN
                     pd(:, :) = block
                  ELSE
                     pd(:, :) = TRANSPOSE(block)
                  END IF
               END IF
            END IF
            factor = MERGE(1.0_dp, 2.0_dp, ia == ib)
            CALL qs_uspp_projector_pair(alist_a, alist_b, center_atoms(center), pairs(6:8, ENTRY), &
                                        onsite=matrix, ao_operator=ph, density=pd, projected=dm, weight=factor, &
                                        atom_gradient=atom_gradient, cell_gradient=cell_gradient)
            IF (PRESENT(ao_operator)) THEN
               CALL dbcsr_get_block_p(ao_operator(image)%matrix, row, col, block, found)
               IF (found) THEN
                  IF (ia <= ib) THEN
                     block(:, :) = block + ph
                  ELSE
                     block(:, :) = block + TRANSPOSE(ph)
                  END IF
               ELSE
                  ! A nonzero projector contribution must never be dropped by AO screening.
                  CPASSERT(ALL(ph == 0.0_dp))
               END IF
               DEALLOCATE (ph)
            END IF
            IF (PRESENT(density)) DEALLOCATE (pd)
         END DO
         IF (PRESENT(projected)) THEN
            CALL group%sum(dm)
            CALL dbcsr_get_block_p(projected, center, center, block, found)
            IF (found) block(:, :) = dm
            DEALLOCATE (dm)
         END IF
         IF (PRESENT(onsite)) DEALLOCATE (matrix)
      END DO
   END SUBROUTINE qs_uspp_projector_images

! **************************************************************************************************
!> \brief Sum all periodic images of <chi_mu|beta_i> into a rectangular DBCSR matrix.
!>        The row distribution must be that used to build the ABBA neighbor list.
!>        ABBA supplies every projector neighbor on each rank owning the AO's row or
!>        column. Selecting B's unique block owner therefore removes replicated work
!>        without discarding images. No filter or Gaussian-overlap sparsity is imposed.
!> \param template AO matrix whose atom blocks and row distribution generated sap_int.
!> \param sap_int Uncompressed overlaps from build_sap_ints (not H-weighted achint).
!> \param center_atoms Global atom indices for the USPP centers, in column-block order.
!> \param projector_sizes Magnetic-projector count per selected center.
!> \param overlaps New Gamma-point B matrix. Non-Gamma phases require a separate path.
! **************************************************************************************************
   SUBROUTINE qs_uspp_projector_matrix(template, sap_int, center_atoms, projector_sizes, overlaps)
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(sap_int_type), INTENT(IN)                     :: sap_int(:)
      INTEGER, INTENT(IN)                                :: center_atoms(:), projector_sizes(:)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: overlaps

      INTEGER                                            :: atom, center, i, j, k, npcols, owner, &
                                                            pass
      INTEGER, ALLOCATABLE                               :: atom_to_center(:)
      INTEGER, POINTER                                   :: col_dist(:), col_sizes(:), row_dist(:), &
                                                            row_sizes(:)
      LOGICAL                                            :: found
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_distribution_type)                      :: dist, projector_dist
      TYPE(mp_comm_type)                                 :: group

      CALL dbcsr_get_info(template, distribution=dist, row_blk_size=row_sizes, group=group)
      CALL dbcsr_distribution_get(dist, row_dist=row_dist, npcols=npcols)
      CPASSERT(SIZE(center_atoms) == SIZE(projector_sizes))
      CPASSERT(SIZE(center_atoms) > 0 .AND. ALL(projector_sizes > 0))
      CPASSERT(ALL(center_atoms >= 1) .AND. ALL(center_atoms <= SIZE(row_sizes)))
      ALLOCATE (atom_to_center(SIZE(row_sizes)), col_sizes(SIZE(center_atoms)), col_dist(SIZE(center_atoms)))
      atom_to_center = 0
      DO i = 1, SIZE(center_atoms)
         CPASSERT(atom_to_center(center_atoms(i)) == 0)
         atom_to_center(center_atoms(i)) = i
         col_dist(i) = MODULO(i - 1, npcols)
      END DO
      col_sizes(:) = projector_sizes
      CALL dbcsr_distribution_new(projector_dist, template=dist, row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_create(overlaps, 'USPP Gaussian projector overlaps', projector_dist, 'N', row_sizes, col_sizes)
      CALL dbcsr_distribution_release(projector_dist)
      DEALLOCATE (col_sizes, col_dist)
      DO pass = 1, 2
      DO i = 1, SIZE(sap_int)
         IF (.NOT. ASSOCIATED(sap_int(i)%alist)) CYCLE
         DO j = 1, SIZE(sap_int(i)%alist)
            IF (.NOT. ASSOCIATED(sap_int(i)%alist(j)%clist)) CYCLE
            atom = sap_int(i)%alist(j)%aatom
            CPASSERT(atom >= 1 .AND. atom <= SIZE(row_sizes))
            DO k = 1, SIZE(sap_int(i)%alist(j)%clist)
               ! Unused entries have catom=0 and no integral storage.
               center = sap_int(i)%alist(j)%clist(k)%catom
               IF (center == 0) CYCLE
               CPASSERT(center >= 1 .AND. center <= SIZE(atom_to_center))
               center = atom_to_center(center)
               IF (center == 0) CYCLE
               CALL dbcsr_get_stored_coordinates(overlaps, atom, center, owner)
               IF (owner /= group%mepos) CYCLE
               CPASSERT(ASSOCIATED(sap_int(i)%alist(j)%clist(k)%acint))
               ASSOCIATE (integrals => sap_int(i)%alist(j)%clist(k)%acint)
                  CPASSERT(SIZE(integrals, 1) == row_sizes(atom))
                  CPASSERT(SIZE(integrals, 2) == projector_sizes(center))
                  CPASSERT(SIZE(integrals, 3) >= 1)
                  IF (pass == 1) THEN
                     CALL dbcsr_reserve_blocks(overlaps, rows=[atom], cols=[center])
                  ELSE
                     CALL dbcsr_get_block_p(overlaps, atom, center, block, found)
                     CPASSERT(found)
                     block(:, :) = block + integrals(:, :, 1)
                  END IF
               END ASSOCIATE
            END DO
         END DO
      END DO
      IF (pass == 1) THEN
         CALL dbcsr_finalize(overlaps)
         CALL dbcsr_set(overlaps, 0.0_dp)
      END IF
      END DO
   END SUBROUTINE qs_uspp_projector_matrix

! **************************************************************************************************
!> \brief Contract dE/dB with all image-resolved AO-center projector derivatives.
!>        Returns energy gradients, with equal and opposite AO/projector-center contributions.
!>        The unique owner of each rectangular block supplies its contribution; no MPI sum here.
!> \param gradient Distributed dE/dB, with the same row/column distribution as B.
!> \param sap_int Uncompressed nder=1 projector overlaps from the ABBA neighbor list.
!> \param center_atoms Global atom indices for the projector column blocks.
!> \param atom_gradient Local energy gradients (3,natom), incremented in place.
!> \param cell_gradient Optional strain energy gradient, including each periodic image.
! **************************************************************************************************
   SUBROUTINE qs_uspp_projector_gradient(gradient, sap_int, center_atoms, atom_gradient, cell_gradient)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: gradient
      TYPE(sap_int_type), INTENT(IN)                     :: sap_int(:)
      INTEGER, INTENT(IN)                                :: center_atoms(:)
      REAL(KIND=dp), INTENT(INOUT)                       :: atom_gradient(:, :)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: cell_gradient(3, 3)

      INTEGER                                            :: atom, catom, center, i, idir, j, k, owner
      INTEGER, ALLOCATABLE                               :: atom_to_center(:)
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: value
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(mp_comm_type)                                 :: group

      CPASSERT(SIZE(atom_gradient, 1) == 3)
      ALLOCATE (atom_to_center(SIZE(atom_gradient, 2)))
      atom_to_center = 0
      DO center = 1, SIZE(center_atoms)
         atom_to_center(center_atoms(center)) = center
      END DO
      CALL dbcsr_get_info(gradient, group=group)
      DO i = 1, SIZE(sap_int)
         IF (.NOT. ASSOCIATED(sap_int(i)%alist)) CYCLE
         DO j = 1, SIZE(sap_int(i)%alist)
            IF (.NOT. ASSOCIATED(sap_int(i)%alist(j)%clist)) CYCLE
            atom = sap_int(i)%alist(j)%aatom
            DO k = 1, SIZE(sap_int(i)%alist(j)%clist)
               catom = sap_int(i)%alist(j)%clist(k)%catom
               IF (catom == 0) CYCLE
               CPASSERT(catom >= 1 .AND. catom <= SIZE(atom_to_center))
               center = atom_to_center(catom)
               IF (center == 0) CYCLE
               CALL dbcsr_get_stored_coordinates(gradient, atom, center, owner)
               IF (owner /= group%mepos) CYCLE
               CALL dbcsr_get_block_p(gradient, atom, center, block, found)
               IF (.NOT. found) CYCLE
               CPASSERT(ASSOCIATED(sap_int(i)%alist(j)%clist(k)%acint))
               ASSOCIATE (integrals => sap_int(i)%alist(j)%clist(k)%acint)
                  CPASSERT(SIZE(integrals, 3) >= 4)
                  DO idir = 1, 3
                     value = SUM(block*integrals(:, :, idir + 1))
                     atom_gradient(idir, atom) = atom_gradient(idir, atom) + value
                     atom_gradient(idir, catom) = atom_gradient(idir, catom) - value
                     IF (PRESENT(cell_gradient)) THEN
                        ! rac points from the AO center to this image of the projector center.
                        cell_gradient(idir, :) = cell_gradient(idir, :) - value*sap_int(i)%alist(j)%clist(k)%rac
                     END IF
                  END DO
               END ASSOCIATE
            END DO
         END DO
      END DO
   END SUBROUTINE qs_uspp_projector_gradient

END MODULE qs_uspp_projectors
