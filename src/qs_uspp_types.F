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
   USE cp_dbcsr_api,                    ONLY: dbcsr_deallocate_matrix,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE upf_projector_integrals,         ONLY: upf_projector_init
   USE uspp_augmentation,               ONLY: uspp_augmentation_init,&
                                              uspp_augmentation_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_uspp_types'

   TYPE qs_uspp_type
      INTEGER, ALLOCATABLE :: center_atoms(:), center_kinds(:), projector_sizes(:)
      REAL(KIND=dp), ALLOCATABLE :: position(:, :)
      TYPE(uspp_augmentation_type), ALLOCATABLE :: augmentation(:)
      TYPE(dbcsr_type), POINTER :: ao_pattern => NULL()
   END TYPE qs_uspp_type

   PUBLIC :: qs_uspp_type, qs_uspp_init, qs_uspp_release

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param uspp ...
!> \param qs_kind_set ...
!> \param particle_set ...
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
            active(ikind) = potential%upf_data%is_ultrasoft .OR. potential%upf_data%is_paw
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
      ALLOCATE (uspp%ao_pattern)
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
!> \brief Release owned augmentation and the inverse AO allocation pattern.
!> \param uspp Owned state, nullified on return.
! **************************************************************************************************
   SUBROUTINE qs_uspp_release(uspp)
      TYPE(qs_uspp_type), POINTER                        :: uspp

      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL dbcsr_deallocate_matrix(uspp%ao_pattern)
      DEALLOCATE (uspp)
      NULLIFY (uspp)
   END SUBROUTINE qs_uspp_release
END MODULE qs_uspp_types
