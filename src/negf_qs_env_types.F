!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Part of the NEGF environment that needs to be stored within qs_environment.
!>        It is used to update the Hartree potential and the Kohn-Sham matrix within qs_scf().
! **************************************************************************************************
MODULE negf_qs_env_types
   USE kinds,                           ONLY: dp
   USE negf_alloc_types,                ONLY: negf_allocatable_rmatrix
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_qs_env_types'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

! **************************************************************************************************
!> \brief ...
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE negf_matrix_block_map_type
      !> column and row atomic indices
      INTEGER                                            :: iatom_col, iatom_row
      !> cell replica
      INTEGER                                            :: image
      !> absolute value is the index of the block within the block data storage;
      !>  0  all matrix elements of the given block are equal to zero,
      !> >0  block is stored simply as block_storage(block_index),
      !> <0  block is stored as TRANSPOSE(block_storage(-block_index))
      INTEGER                                            :: block_index
   END TYPE negf_matrix_block_map_type

! **************************************************************************************************
!> \brief Part of the NEGF environment stored within the QuickStep environment.
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE negf_qs_env_contact_type
      REAL(kind=dp), DIMENSION(3)                        :: direction_vector, origin_scattering
      REAL(kind=dp)                                      :: v_external
      !> atom_row, atom_col, image, and the index of the actual data block in 'block_storage'
      !> all indices are identical for both spins
      TYPE(negf_matrix_block_map_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: ks_block_map
      !> collection of DBSCR Kohn-Sham matrix blocks ([index, spin])
      TYPE(negf_allocatable_rmatrix), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: ks_block_storage
   END TYPE negf_qs_env_contact_type

! **************************************************************************************************
!> \brief data to update the Hartree potential and the Kohn-Sham matrix within qs_scf().
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE negf_qs_env_type
      LOGICAL                                            :: apply_bias
      LOGICAL                                            :: update_ks_matrix
      TYPE(negf_qs_env_contact_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: contacts
      !> reference Hartree potential (with no bias)
      TYPE(pw_type), POINTER                             :: v_hartree_ref
   END TYPE negf_qs_env_type

   PUBLIC :: negf_qs_env_type, negf_qs_env_contact_type, negf_matrix_block_map_type
   PUBLIC :: negf_qs_env_create, negf_qs_env_release

CONTAINS

! **************************************************************************************************
!> \brief Allocate control options for Non-equilibrium Green's Function calculation.
!> \param negf_qs_env an object to create
!> \par History
!>    * 10.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_qs_env_create(negf_qs_env)
      TYPE(negf_qs_env_type), POINTER                    :: negf_qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_qs_env_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CPASSERT(.NOT. ASSOCIATED(negf_qs_env))
      CALL timeset(routineN, handle)

      ALLOCATE (negf_qs_env)
      negf_qs_env%apply_bias = .FALSE.
      negf_qs_env%update_ks_matrix = .FALSE.
      NULLIFY (negf_qs_env%v_hartree_ref)

      CALL timestop(handle)
   END SUBROUTINE negf_qs_env_create

! **************************************************************************************************
!> \brief Release memory allocated for NEGF control options.
!> \param negf_qs_env an object to release
!> \param pw_env      plane-wave environment used to create the reference Hartree potential
!> \par History
!>    * 10.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_qs_env_release(negf_qs_env, pw_env)
      TYPE(negf_qs_env_type), POINTER                    :: negf_qs_env
      TYPE(pw_env_type), POINTER                         :: pw_env

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_qs_env_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iblock, icontact, ispin
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(negf_qs_env)) THEN
         IF (ALLOCATED(negf_qs_env%contacts)) THEN
            DO icontact = SIZE(negf_qs_env%contacts), 1, -1
               IF (ALLOCATED(negf_qs_env%contacts(icontact)%ks_block_map)) &
                  DEALLOCATE (negf_qs_env%contacts(icontact)%ks_block_map)

               IF (ALLOCATED(negf_qs_env%contacts(icontact)%ks_block_storage)) THEN
                  DO ispin = SIZE(negf_qs_env%contacts(icontact)%ks_block_storage, 2), 1, -1
                     DO iblock = SIZE(negf_qs_env%contacts(icontact)%ks_block_storage, 1), 1, -1
                        IF (ALLOCATED(negf_qs_env%contacts(icontact)%ks_block_storage(iblock, ispin)%matrix)) &
                           DEALLOCATE (negf_qs_env%contacts(icontact)%ks_block_storage(iblock, ispin)%matrix)
                     END DO
                  END DO

                  DEALLOCATE (negf_qs_env%contacts(icontact)%ks_block_storage)
               END IF
            END DO

            DEALLOCATE (negf_qs_env%contacts)
         END IF

         IF (ASSOCIATED(negf_qs_env%v_hartree_ref)) THEN
            CALL pw_env_get(pw_env, auxbas_pw_pool=pw_pool)
            CALL pw_pool_give_back_pw(pw_pool, negf_qs_env%v_hartree_ref)
         END IF

         DEALLOCATE (negf_qs_env)
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_qs_env_release
END MODULE negf_qs_env_types
