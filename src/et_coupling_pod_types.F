!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief data types and related routines used in POD method
!> \author Z. Futera (03.2019)
! **************************************************************************************************
MODULE et_coupling_pod_types

   USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              mo_set_p_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_pod_types'

   PUBLIC :: et_cpl, &
             et_cpl_block, & 
             et_cpl_atom

   ! Electronic-coupling calculation data structure
   !
   ! n_atoms  - number of atoms in the blocks
   ! n_blocks - number of atomic blocks (donor,acceptor,bridge,...)
   ! fermi    - system Fermi level (alpha/beta spin component)
   ! tmatf    - transformation matrix for basis-set orthogonalization (S^{-1/2})
   ! tmatr    - inversion transformation matrix
   ! block    - atomic data blocks
   TYPE et_cpl
      INTEGER                                            :: n_atoms
      INTEGER                                            :: n_blocks
      REAL(KIND=dp), DIMENSION(:), POINTER               :: fermi
      TYPE(cp_fm_type), POINTER                          :: tmatf
      TYPE(cp_fm_type), POINTER                          :: tmatr
      TYPE(et_cpl_block), DIMENSION(:), POINTER          :: block
   END TYPE et_cpl

   ! Electronic-coupling data block
   !
   ! n_atoms     - number of atoms
   ! n_electrons - number of electrons
   ! n_ao        - number of AO basis functions
   ! atom        - list of atoms
   ! mo          - electronic states
   ! hab         - electronic-coupling elements
   TYPE et_cpl_block
      INTEGER                                            :: n_atoms
      INTEGER                                            :: n_electrons
      INTEGER                                            :: n_ao
      INTEGER                                            :: n_hab_blocks
      TYPE(et_cpl_atom), DIMENSION(:), POINTER           :: atom
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: hab
   END TYPE et_cpl_block

   ! Electronic-coupling block-atom data
   ! id     - atom ID
   ! n_ao   - number of AO basis functions
   ! ao_pos - position of atom in array of AO functions
   TYPE et_cpl_atom
      INTEGER                                            :: id
      INTEGER                                            :: n_ao
      INTEGER                                            :: ao_pos
   END TYPE et_cpl_atom

   PUBLIC :: pod_data_release

CONTAINS

! **************************************************************************************************
!> \brief Release memory allocate for electronic coupling data structures
!> \param ec electronic coupling data structure
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_data_release(ec)

      ! Routine arguments
      TYPE(et_cpl), POINTER                              :: ec

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_data_release', &
                                     routineP = moduleN//':'//routineN
      ! Local variables
      INTEGER                                            :: i, j, k

      IF (ASSOCIATED(ec)) THEN

         IF (ASSOCIATED(ec%fermi)) &
            DEALLOCATE (ec%fermi)
         IF (ASSOCIATED(ec%tmatf)) &
            CALL cp_fm_release(matrix=ec%tmatf)
         IF (ASSOCIATED(ec%tmatr)) &
            CALL cp_fm_release(matrix=ec%tmatr)

         IF (ASSOCIATED(ec%block)) THEN

            DO i = 1, SIZE(ec%block)
               IF (ASSOCIATED(ec%block(i)%atom)) &
                  DEALLOCATE (ec%block(i)%atom)
               IF (ASSOCIATED(ec%block(i)%mo)) THEN
                  DO j = 1, SIZE(ec%block(i)%mo)
                     IF (ASSOCIATED(ec%block(i)%mo(j)%mo_set)) &
                        CALL deallocate_mo_set(ec%block(i)%mo(j)%mo_set)
                  END DO
                  DEALLOCATE (ec%block(i)%mo)
               END IF
               IF (ASSOCIATED(ec%block(i)%hab)) THEN
                  DO j = 1, SIZE(ec%block(i)%hab, 1)
                     DO k = 1, SIZE(ec%block(i)%hab, 2)
                        IF (ASSOCIATED(ec%block(i)%hab(j, k)%matrix)) &
                           CALL cp_fm_release(matrix=ec%block(i)%hab(j, k)%matrix)
                     END DO
                  END DO
                  DEALLOCATE (ec%block(i)%hab)
               END IF
            END DO

            DEALLOCATE (ec%block)

         END IF

         DEALLOCATE (ec)

      END IF

   END SUBROUTINE pod_data_release

END MODULE et_coupling_pod_types
