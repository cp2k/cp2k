!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Geometry-dependent scalar USPP state, separate from the Gaussian overlap.
! **************************************************************************************************
MODULE qs_uspp_types
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_release_p,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE kpoint_transitional,             ONLY: kpoint_transitional_release,&
                                              kpoint_transitional_type,&
                                              set_1d_pointer,&
                                              set_2d_pointer
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_uspp_projectors,              ONLY: qs_uspp_projector_images,&
                                              qs_uspp_projector_matrix
   USE sap_kind_types,                  ONLY: release_sap_int,&
                                              sap_int_type,&
                                              sap_sort
   USE upf_projector_integrals,         ONLY: upf_projector_init
   USE uspp_augmentation,               ONLY: uspp_augmentation_init,&
                                              uspp_augmentation_type
   USE uspp_projector_operators,        ONLY: uspp_projector_ao_pattern,&
                                              uspp_projector_create_onsite,&
                                              uspp_projector_metric
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_uspp_types'

   TYPE qs_uspp_type
      INTEGER, ALLOCATABLE :: center_atoms(:), center_kinds(:), projector_sizes(:)
      REAL(KIND=dp), ALLOCATABLE :: position(:, :)
      TYPE(uspp_augmentation_type), ALLOCATABLE :: augmentation(:)
      TYPE(dbcsr_type), POINTER :: overlaps => NULL(), ao_pattern => NULL()
      TYPE(kpoint_transitional_type) :: metric
      TYPE(sap_int_type), POINTER :: image_integrals(:) => NULL()
      INTEGER, ALLOCATABLE :: image_pairs(:, :)
   END TYPE qs_uspp_type

   PUBLIC :: qs_uspp_type, qs_uspp_init, qs_uspp_update, qs_uspp_update_images, qs_uspp_release

CONTAINS

! **************************************************************************************************
!> \brief Select actual USPP centers and prepare per-kind augmentation in native beta order.
!>        No state is allocated for a system without USPP centers. Ghost atoms never augment.
!> \param uspp Owned state, replaced if already present.
!> \param qs_kind_set Kinds with original scalar USPP data and native projector caches.
!> \param particle_set Global particles, indexed in the same order as AO matrix blocks.
! **************************************************************************************************
   SUBROUTINE qs_uspp_init(uspp, qs_kind_set, particle_set)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(qs_kind_type), INTENT(INOUT)                  :: qs_kind_set(:)
      TYPE(particle_type), INTENT(IN)                    :: particle_set(:)

      INTEGER                                            :: atom, center, degree, ikind, ncenter
      LOGICAL                                            :: active(SIZE(qs_kind_set))

      CALL qs_uspp_release(uspp)
      active = .FALSE.
      DO ikind = 1, SIZE(qs_kind_set)
         IF (qs_kind_set(ikind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(qs_kind_set(ikind)%sgp_potential)) CYCLE
         ASSOCIATE (potential => qs_kind_set(ikind)%sgp_potential)
            IF (.NOT. ALLOCATED(potential%upf_data)) CYCLE
            active(ikind) = potential%upf_data%is_ultrasoft
         END ASSOCIATE
      END DO
      IF (.NOT. ANY(active)) RETURN
      ncenter = 0
      DO atom = 1, SIZE(particle_set)
         CPASSERT(ASSOCIATED(particle_set(atom)%atomic_kind))
         ikind = particle_set(atom)%atomic_kind%kind_number
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(qs_kind_set))
         IF (active(ikind)) ncenter = ncenter + 1
      END DO
      IF (ncenter == 0) RETURN
      ALLOCATE (uspp)
      ALLOCATE (uspp%center_atoms(ncenter), uspp%center_kinds(ncenter), uspp%projector_sizes(ncenter))
      ALLOCATE (uspp%position(3, ncenter), uspp%augmentation(SIZE(qs_kind_set)))
      center = 0
      DO atom = 1, SIZE(particle_set)
         ikind = particle_set(atom)%atomic_kind%kind_number
         IF (.NOT. active(ikind)) CYCLE
         ASSOCIATE (potential => qs_kind_set(ikind)%sgp_potential)
            IF (uspp%augmentation(ikind)%nproj == 0) THEN
               CPASSERT(potential%upf_data%number_of_proj > 0)
               degree = MAXVAL(potential%upf_data%lbeta)
               IF (potential%upf_projectors%degree < degree) THEN
                  CALL upf_projector_init(potential%upf_data, degree, potential%upf_projectors)
               END IF
               CALL uspp_augmentation_init(potential%upf_data, potential%upf_projectors, uspp%augmentation(ikind))
            END IF
         END ASSOCIATE
         center = center + 1
         uspp%center_atoms(center) = atom
         uspp%center_kinds(center) = ikind
         uspp%projector_sizes(center) = uspp%augmentation(ikind)%nproj
         uspp%position(:, center) = particle_set(atom)%r
      END DO
   END SUBROUTINE qs_uspp_init

! **************************************************************************************************
!> \brief Rebuild B and S from current neighbor integrals, allowing projector-mediated fill-in.
!>        The Gaussian overlap and its derivatives remain owned by the ordinary KS environment.
!> \param uspp Initialized USPP state.
!> \param gaussian_overlap Current Gaussian overlap and AO block distribution.
!> \param sap_int Current uncompressed AO/projector integrals including periodic images.
!> \param particle_set Current particle positions; the kind mapping must be unchanged.
! **************************************************************************************************
   SUBROUTINE qs_uspp_update(uspp, gaussian_overlap, sap_int, particle_set)
      TYPE(qs_uspp_type), INTENT(INOUT)                  :: uspp
      TYPE(dbcsr_type), INTENT(IN)                       :: gaussian_overlap
      TYPE(sap_int_type), INTENT(IN)                     :: sap_int(:)
      TYPE(particle_type), INTENT(IN)                    :: particle_set(:)

      INTEGER                                            :: atom, center
      TYPE(dbcsr_p_type), POINTER                        :: metric(:)

      CALL dbcsr_release_p(uspp%overlaps)
      CALL dbcsr_release_p(uspp%ao_pattern)
      CALL kpoint_transitional_release(uspp%metric)
      IF (ASSOCIATED(uspp%image_integrals)) CALL release_sap_int(uspp%image_integrals)
      IF (ALLOCATED(uspp%image_pairs)) DEALLOCATE (uspp%image_pairs)
      DO center = 1, SIZE(uspp%center_atoms)
         atom = uspp%center_atoms(center)
         CPASSERT(atom >= 1 .AND. atom <= SIZE(particle_set))
         CPASSERT(ASSOCIATED(particle_set(atom)%atomic_kind))
         CPASSERT(particle_set(atom)%atomic_kind%kind_number == uspp%center_kinds(center))
         uspp%position(:, center) = particle_set(atom)%r
      END DO
      ALLOCATE (uspp%overlaps)
      CALL qs_uspp_projector_matrix(gaussian_overlap, sap_int, uspp%center_atoms, uspp%projector_sizes, uspp%overlaps)
      ALLOCATE (uspp%ao_pattern)
      CALL uspp_projector_ao_pattern(uspp%overlaps, gaussian_overlap, uspp%ao_pattern)
      ALLOCATE (metric(1))
      ALLOCATE (metric(1)%matrix)
      CALL uspp_projector_metric(uspp%overlaps, uspp%augmentation, uspp%center_kinds, gaussian_overlap, metric(1)%matrix)
      CALL set_1d_pointer(uspp%metric, metric)
   END SUBROUTINE qs_uspp_update

! **************************************************************************************************
!> \brief Retain image-resolved overlaps and build the generalized real-space metric.
!> \param uspp Initialized USPP state.
!> \param gaussian_overlap Gaussian overlap matrices for every lattice image.
!> \param sap_int Overlaps whose ownership transfers to uspp; null on return.
!> \param particle_set Current atoms.
!> \param sab Complete AO pair list used by the lattice Fourier transformation.
!> \param cell_to_index Current lattice image map.
! **************************************************************************************************
   SUBROUTINE qs_uspp_update_images(uspp, gaussian_overlap, sap_int, particle_set, sab, cell_to_index)
      TYPE(qs_uspp_type), INTENT(INOUT)                  :: uspp
      TYPE(dbcsr_p_type), INTENT(IN)                     :: gaussian_overlap(:)
      TYPE(sap_int_type), POINTER                        :: sap_int(:)
      TYPE(particle_type), INTENT(IN)                    :: particle_set(:)
      TYPE(neighbor_list_set_p_type), INTENT(IN)         :: sab(:)
      INTEGER, POINTER                                   :: cell_to_index(:, :, :)

      INTEGER                                            :: cell(3), center, image, slot
      LOGICAL                                            :: found
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_p_type), POINTER                        :: metric(:, :)
      TYPE(dbcsr_type)                                   :: onsite

      CALL qs_uspp_update(uspp, gaussian_overlap(1)%matrix, sap_int, particle_set)
      uspp%image_integrals => sap_int
      NULLIFY (sap_int)
      CALL sap_sort(uspp%image_integrals)
      ALLOCATE (uspp%image_pairs(8, sab(1)%nl_size))
      DO slot = 1, sab(1)%nl_size
         cell = sab(1)%nlist_task(slot)%cell
         image = cell_to_index(cell(1), cell(2), cell(3))
         CPASSERT(image >= 1 .AND. image <= SIZE(gaussian_overlap))
         uspp%image_pairs(:, slot) = [sab(1)%nlist_task(slot)%iatom, sab(1)%nlist_task(slot)%jatom, &
                                      sab(1)%nlist_task(slot)%ikind, sab(1)%nlist_task(slot)%jkind, image, cell]
      END DO
      CALL kpoint_transitional_release(uspp%metric)
      ALLOCATE (metric(1, SIZE(gaussian_overlap)))
      DO image = 1, SIZE(gaussian_overlap)
         ALLOCATE (metric(1, image)%matrix)
         CALL dbcsr_copy(metric(1, image)%matrix, uspp%ao_pattern)
         CALL dbcsr_set(metric(1, image)%matrix, 0.0_dp)
         CALL dbcsr_add(metric(1, image)%matrix, gaussian_overlap(image)%matrix, 1.0_dp, 1.0_dp)
      END DO
      CALL uspp_projector_create_onsite(uspp%overlaps, onsite)
      DO center = 1, SIZE(uspp%center_atoms)
         CALL dbcsr_get_block_p(onsite, center, center, block, found)
         IF (found) block(:, :) = uspp%augmentation(uspp%center_kinds(center))%q
      END DO
      CALL qs_uspp_projector_images(uspp%overlaps, uspp%image_integrals, uspp%image_pairs, SIZE(uspp%augmentation), &
                                    uspp%center_atoms, uspp%center_kinds, onsite=onsite, ao_operator=metric(1, :))
      CALL dbcsr_release(onsite)
      CALL set_2d_pointer(uspp%metric, metric)
   END SUBROUTINE qs_uspp_update_images

! **************************************************************************************************
!> \brief Release only USPP-owned matrices and data; safe also for an absent state.
!> \param uspp Owned state, null on return.
! **************************************************************************************************
   SUBROUTINE qs_uspp_release(uspp)
      TYPE(qs_uspp_type), POINTER                        :: uspp

      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL dbcsr_release_p(uspp%overlaps)
      CALL dbcsr_release_p(uspp%ao_pattern)
      CALL kpoint_transitional_release(uspp%metric)
      IF (ASSOCIATED(uspp%image_integrals)) CALL release_sap_int(uspp%image_integrals)
      DEALLOCATE (uspp)
   END SUBROUTINE qs_uspp_release

END MODULE qs_uspp_types
