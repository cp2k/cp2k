!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to deal with vectors in 3-D real space.
! **************************************************************************************************

MODULE negf_vectors
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_vectors'

   PUBLIC :: contact_direction_vector, projection_on_direction_vector

CONTAINS

! **************************************************************************************************
!> \brief compute direction vector of the given contact
!> \param origin                origin
!> \param direction_vector      direction vector
!> \param origin_bias           origin which will be used to apply the external bias
!>                              (in contrast with 'origin' it does not include screening region)
!> \param direction_vector_bias direction vector which will be used to apply the external bias
!>                              (together with 'origin_bias' it defines a contact region where
!>                              the external potential is kept constant)
!> \param atomlist_screening    atoms belonging to the contact's screening region
!> \param atomlist_bulk         atoms belonging to the contact's bulk region
!> \param subsys                QuickStep subsystem
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE contact_direction_vector(origin, direction_vector, origin_bias, direction_vector_bias, &
                                       atomlist_screening, atomlist_bulk, subsys)
      REAL(kind=dp), DIMENSION(3), INTENT(out)           :: origin, direction_vector, origin_bias, &
                                                            direction_vector_bias
      INTEGER, DIMENSION(:), INTENT(in)                  :: atomlist_screening, atomlist_bulk
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contact_direction_vector'

      INTEGER                                            :: handle, iatom, natoms_bulk, &
                                                            natoms_screening, nparticles
      REAL(kind=dp)                                      :: proj, proj_max, proj_min, proj_min_bias
      REAL(kind=dp), DIMENSION(3)                        :: vector
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      CALL qs_subsys_get(subsys, particle_set=particle_set)

      natoms_screening = SIZE(atomlist_screening)
      natoms_bulk = SIZE(atomlist_bulk)
      nparticles = SIZE(particle_set)

      CPASSERT(natoms_screening > 0)
      CPASSERT(natoms_bulk > 0)
      CPASSERT(nparticles > 0)

      ! geometrical centre of the first contact unit cell
      origin = particle_set(atomlist_screening(1))%r
      DO iatom = 2, natoms_screening
         origin = origin + particle_set(atomlist_screening(iatom))%r
      END DO
      origin = origin/REAL(natoms_screening, kind=dp)

      ! geometrical centre of the second contact unit cell
      direction_vector = particle_set(atomlist_bulk(1))%r
      DO iatom = 2, natoms_bulk
         direction_vector = direction_vector + particle_set(atomlist_bulk(iatom))%r
      END DO
      direction_vector = direction_vector/REAL(natoms_bulk, kind=dp)

      ! vector between the geometrical centers of the first and the second contact unit cells
      direction_vector = direction_vector - origin

      ! the point 'center_of_coords0' belongs to the first unit cell, so the lowest projection of any point
      ! from the first unit cell on the direction vector 'center_of_coords1 - center_of_coords0' should be <= 0
      proj_min = 0.0_dp
      proj_min_bias = 0.0_dp
      DO iatom = 1, natoms_screening
         vector = particle_set(atomlist_screening(iatom))%r - origin
         proj = projection_on_direction_vector(vector, direction_vector)

         IF (proj < proj_min) proj_min = proj
         IF (proj > proj_min_bias) proj_min_bias = proj
      END DO

      ! the point 'center_of_coords1' belongs to the given contact, so the highest projection should be >= 1
      proj_max = 1.0_dp
      DO iatom = 1, nparticles
         vector = particle_set(iatom)%r - origin
         proj = projection_on_direction_vector(vector, direction_vector)

         IF (proj > proj_max) proj_max = proj
      END DO

      ! adjust the origin, so it lies on a plane between the given contact and the scattering region
      origin_bias = origin + proj_min_bias*direction_vector
      origin = origin + proj_min*direction_vector

      ! rescale the vector, so the last atom of the given contact and the point 'origin + direction_vector' lie on
      ! the same plane parallel to the 'origin' plane -- which separates the contact from the scattering region.
      direction_vector_bias = (proj_max - proj_min_bias)*direction_vector
      direction_vector = (proj_max - proj_min)*direction_vector

      CALL timestop(handle)
   END SUBROUTINE contact_direction_vector

! **************************************************************************************************
!> \brief project the 'vector' onto the direction 'vector0'. Both vectors should have the same origin.
!> \param vector   vector to project
!> \param vector0  direction vector
!> \return projection (in terms of the direction vector's length)
!> \author Sergey Chulkov
!> \note \parblock
!>          & vector
!>         /
!>        /
!>       /
!>      *---|------> vector0
!>       l  = proj * | vector0 |
!>\endparblock
! **************************************************************************************************
   PURE FUNCTION projection_on_direction_vector(vector, vector0) RESULT(proj)
      REAL(kind=dp), DIMENSION(3), INTENT(in)            :: vector, vector0
      REAL(kind=dp)                                      :: proj

      REAL(kind=dp)                                      :: len2, len2_v0, len2_v1
      REAL(kind=dp), DIMENSION(3)                        :: vector1

      vector1 = vector - vector0
      len2 = DOT_PRODUCT(vector, vector)
      len2_v0 = DOT_PRODUCT(vector0, vector0)
      len2_v1 = DOT_PRODUCT(vector1, vector1)

      proj = 0.5_dp*((len2 - len2_v1)/len2_v0 + 1.0_dp)
   END FUNCTION projection_on_direction_vector
END MODULE negf_vectors
