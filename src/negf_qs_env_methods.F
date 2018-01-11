!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Derived types for NEGF based quantum transport calculations
! **************************************************************************************************
MODULE negf_qs_env_methods
   USE dbcsr_api,                       ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE negf_qs_env_types,               ONLY: negf_qs_env_contact_type,&
                                              negf_qs_env_type
   USE negf_vectors,                    ONLY: projection_on_direction_vector
   USE pw_types,                        ONLY: pw_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_qs_env_methods'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   PUBLIC :: negf_update_ks_matrix, negf_modify_hartree_pot

CONTAINS

! **************************************************************************************************
!> \brief Ensure that Kohn-Sham matrix blocks which corresponds to a bulk contact region are
!>        identical to the actual Kohn-Sham matrix blocks from bulk calculation.
!> \param matrix_ks_kp  Kohn-Sham matrix to update
!> \param negf_qs_env   NEGF_QS environment
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_update_ks_matrix(matrix_ks_kp, negf_qs_env)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp
      TYPE(negf_qs_env_type), POINTER                    :: negf_qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_update_ks_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER :: block_index, handle, iatom_col, iatom_row, iblock, icontact, image, ispin, &
         ispin_contact, nblocks, ncontacts, nimages, nspins, nspins_contacts
      LOGICAL                                            :: found, is_transpose
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: sm_block
      TYPE(negf_qs_env_contact_type), POINTER            :: contact

      CALL timeset(routineN, handle)

      IF (ALLOCATED(negf_qs_env%contacts)) THEN
         nspins = SIZE(matrix_ks_kp, 1)
         nimages = SIZE(matrix_ks_kp, 2)

         ncontacts = SIZE(negf_qs_env%contacts)

         DO icontact = 1, ncontacts
            contact => negf_qs_env%contacts(icontact)
            nblocks = SIZE(contact%ks_block_map)
            nspins_contacts = SIZE(contact%ks_block_storage, 2)

            ispin_contact = 0
            DO ispin = 1, nspins
               IF (ispin_contact < nspins_contacts) &
                  ispin_contact = ispin_contact+1

               DO iblock = 1, nblocks
                  iatom_row = contact%ks_block_map(iblock)%iatom_row
                  iatom_col = contact%ks_block_map(iblock)%iatom_col
                  image = contact%ks_block_map(iblock)%image

                  CPASSERT(image <= nimages)

                  block_index = contact%ks_block_map(iblock)%block_index
                  is_transpose = block_index < 0
                  block_index = ABS(block_index)

                  CALL dbcsr_get_block_p(matrix=matrix_ks_kp(ispin, image)%matrix, &
                                         row=iatom_row, col=iatom_col, block=sm_block, found=found)
                  CPASSERT(found)

                  IF (block_index > 0) THEN
                     IF (is_transpose) THEN
                        sm_block(:, :) = &
                           TRANSPOSE(contact%ks_block_storage(block_index, ispin_contact)%matrix)
                     ELSE
                        sm_block(:, :) = &
                           contact%ks_block_storage(block_index, ispin_contact)%matrix(:, :)
                     END IF
                  ELSE
                     sm_block(:, :) = 0.0_dp
                  END IF
               END DO
            END DO
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_update_ks_matrix

! **************************************************************************************************
!> \brief Modify the Hartree potential by applying an external potential.
!> \param v_hartree     Hartree potential to update
!> \param negf_qs_env   NEGF_QS environment
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_modify_hartree_pot(v_hartree, negf_qs_env)
      TYPE(pw_p_type), INTENT(inout)                     :: v_hartree
      TYPE(negf_qs_env_type), POINTER                    :: negf_qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_modify_hartree_pot', &
         routineP = moduleN//':'//routineN
      REAL(kind=dp), PARAMETER :: threshold = 16.0_dp*EPSILON(0.0_dp)

      INTEGER                                            :: dx, dy, dz, handle, icontact, ix, iy, &
                                                            iz, lx, ly, lz, ncontacts, ux, uy, uz
      REAL(kind=dp)                                      :: pot, proj, v1, v2
      REAL(kind=dp), DIMENSION(3)                        :: dirvector_bias, point_coord, &
                                                            point_indices, vector

      CALL timeset(routineN, handle)

      IF (ALLOCATED(negf_qs_env%contacts)) THEN
         ncontacts = SIZE(negf_qs_env%contacts)
         CPASSERT(ncontacts == 2)

         dirvector_bias = negf_qs_env%contacts(2)%origin_scattering-negf_qs_env%contacts(1)%origin_scattering
         v1 = negf_qs_env%contacts(1)%v_external
         v2 = negf_qs_env%contacts(2)%v_external

         lx = v_hartree%pw%pw_grid%bounds_local(1, 1)
         ux = v_hartree%pw%pw_grid%bounds_local(2, 1)
         ly = v_hartree%pw%pw_grid%bounds_local(1, 2)
         uy = v_hartree%pw%pw_grid%bounds_local(2, 2)
         lz = v_hartree%pw%pw_grid%bounds_local(1, 3)
         uz = v_hartree%pw%pw_grid%bounds_local(2, 3)

         dx = v_hartree%pw%pw_grid%npts(1)/2
         dy = v_hartree%pw%pw_grid%npts(2)/2
         dz = v_hartree%pw%pw_grid%npts(3)/2

         DO iz = lz, uz
            point_indices(3) = REAL(iz+dz, kind=dp)
            DO iy = ly, uy
               point_indices(2) = REAL(iy+dy, kind=dp)

               DO ix = lx, ux
                  point_indices(1) = REAL(ix+dx, kind=dp)
                  point_coord(:) = MATMUL(v_hartree%pw%pw_grid%dh, point_indices)

                  vector = point_coord-negf_qs_env%contacts(1)%origin_scattering
                  proj = projection_on_direction_vector(vector, dirvector_bias)
                  IF (proj+threshold >= 0.0_dp .AND. proj-threshold <= 1.0_dp) THEN
                     ! scattering region
                     ! proj == 0   we are at the first contact boundary
                     ! proj == 1   we are at the second contact boundary
                     IF (proj < 0.0_dp) THEN
                        proj = 0.0_dp
                     ELSE IF (proj > 1.0_dp) THEN
                        proj = 1.0_dp
                     END IF
                     pot = v1+(v2-v1)*proj
                  ELSE
                     pot = 0.0_dp
                     DO icontact = 1, ncontacts
                        vector = point_coord-negf_qs_env%contacts(icontact)%origin_scattering
                        proj = projection_on_direction_vector(vector, negf_qs_env%contacts(icontact)%direction_vector)

                        IF (proj+threshold >= 0.0_dp .AND. proj-threshold <= 1.0_dp) THEN
                           pot = negf_qs_env%contacts(icontact)%v_external
                           EXIT
                        END IF
                     END DO
                  END IF

                  v_hartree%pw%cr3d(ix, iy, iz) = v_hartree%pw%cr3d(ix, iy, iz)+pot*v_hartree%pw%pw_grid%dvol
               END DO
            END DO
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_modify_hartree_pot
END MODULE negf_qs_env_methods
