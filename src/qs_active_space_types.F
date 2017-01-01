!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The types needed for the calculation of active space Hamiltonians
!> \par History
!>      04.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE qs_active_space_types

   USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                              cp_fm_release
   USE dbcsr_api,                       ONLY: csr_destroy,&
                                              csr_p_type,&
                                              dbcsr_deallocate_matrix_set,&
                                              dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              mo_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_active_space_types'

   PUBLIC :: active_space_type, eri_type
   PUBLIC :: create_active_space_type, release_active_space_type

! **************************************************************************************************
!> \brief Quantities needed for AS determination
!> \author JGH
! **************************************************************************************************
   TYPE eri_gpw_type
      LOGICAL                       :: redo_poisson
      REAL(KIND=dp)                 :: cutoff
      REAL(KIND=dp)                 :: rel_cutoff
      REAL(KIND=dp)                 :: eps_grid
      INTEGER                       :: print_level
      LOGICAL                       :: store_wfn
   END TYPE eri_gpw_type

   TYPE eri_type
      INTEGER                       :: method
      INTEGER                       :: OPERATOR
      REAL(KIND=dp)                 :: operator_parameter
      INTEGER, DIMENSION(3)         :: periodicity
      REAL(KIND=dp)                 :: cutoff_radius
      REAL(KIND=dp)                 :: eps_integral
      TYPE(eri_gpw_type)            :: eri_gpw
      TYPE(csr_p_type), &
         DIMENSION(:), POINTER      :: eri => NULL()
   END TYPE eri_type

   TYPE active_space_type
      INTEGER                       :: nactive
      INTEGER                       :: ninactive
      INTEGER, DIMENSION(2)         :: ninspin
      INTEGER, DIMENSION(2)         :: nelectrons
      INTEGER                       :: multiplicity
      INTEGER                       :: nspins
      LOGICAL                       :: molecule
      INTEGER                       :: model
      REAL(KIND=dp)                 :: energy_total
      REAL(KIND=dp)                 :: energy_ref
      REAL(KIND=dp)                 :: energy_inactive
      REAL(KIND=dp)                 :: energy_active
      TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos_active
      TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos_inactive
      TYPE(eri_type)                :: eri
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER   :: p_ref
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER   :: ks_sub
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER   :: vxc_sub
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER   :: h_sub
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER   :: fock_sub
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER:: pmat_inactive
   END TYPE active_space_type

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param active_space_env ...
! **************************************************************************************************
   SUBROUTINE create_active_space_type(active_space_env)
      TYPE(active_space_type), POINTER                   :: active_space_env

      CHARACTER(len=*), PARAMETER :: routineN = 'create_active_space_type', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(active_space_env)) THEN
         CALL release_active_space_type(active_space_env)
      END IF

      ALLOCATE (active_space_env)
      NULLIFY (active_space_env%mos_active, active_space_env%mos_inactive)
      NULLIFY (active_space_env%ks_sub, active_space_env%p_ref)
      NULLIFY (active_space_env%vxc_sub, active_space_env%h_sub)
      NULLIFY (active_space_env%fock_sub, active_space_env%pmat_inactive)

   END SUBROUTINE create_active_space_type

! **************************************************************************************************
!> \brief ...
!> \param active_space_env ...
! **************************************************************************************************
   SUBROUTINE release_active_space_type(active_space_env)
      TYPE(active_space_type), POINTER                   :: active_space_env

      CHARACTER(len=*), PARAMETER :: routineN = 'release_active_space_type', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: imo, isp

      IF (ASSOCIATED(active_space_env)) THEN

         IF (ASSOCIATED(active_space_env%mos_active)) THEN
            DO imo = 1, SIZE(active_space_env%mos_active)
               CALL deallocate_mo_set(active_space_env%mos_active(imo)%mo_set)
            END DO
            DEALLOCATE (active_space_env%mos_active)
         END IF

         IF (ASSOCIATED(active_space_env%mos_inactive)) THEN
            DO imo = 1, SIZE(active_space_env%mos_inactive)
               CALL deallocate_mo_set(active_space_env%mos_inactive(imo)%mo_set)
            END DO
            DEALLOCATE (active_space_env%mos_inactive)
         END IF

         CALL release_eri_type(active_space_env%eri)

         IF (ASSOCIATED(active_space_env%p_ref)) THEN
            DO isp = 1, SIZE(active_space_env%p_ref)
               CALL cp_fm_release(active_space_env%p_ref(isp)%matrix)
            END DO
            DEALLOCATE (active_space_env%p_ref)
         END IF

         IF (ASSOCIATED(active_space_env%ks_sub)) THEN
            DO isp = 1, SIZE(active_space_env%ks_sub)
               CALL cp_fm_release(active_space_env%ks_sub(isp)%matrix)
            END DO
            DEALLOCATE (active_space_env%ks_sub)
         END IF

         IF (ASSOCIATED(active_space_env%vxc_sub)) THEN
            DO isp = 1, SIZE(active_space_env%vxc_sub)
               CALL cp_fm_release(active_space_env%vxc_sub(isp)%matrix)
            END DO
            DEALLOCATE (active_space_env%vxc_sub)
         END IF

         IF (ASSOCIATED(active_space_env%h_sub)) THEN
            DO isp = 1, SIZE(active_space_env%h_sub)
               CALL cp_fm_release(active_space_env%h_sub(isp)%matrix)
            END DO
            DEALLOCATE (active_space_env%h_sub)
         END IF

         IF (ASSOCIATED(active_space_env%fock_sub)) THEN
            DO isp = 1, SIZE(active_space_env%fock_sub)
               CALL cp_fm_release(active_space_env%fock_sub(isp)%matrix)
            END DO
            DEALLOCATE (active_space_env%fock_sub)
         END IF

         IF (ASSOCIATED(active_space_env%pmat_inactive)) &
            CALL dbcsr_deallocate_matrix_set(active_space_env%pmat_inactive)

         DEALLOCATE (active_space_env)
      END IF

   END SUBROUTINE release_active_space_type

! **************************************************************************************************
!> \brief ...
!> \param eri_env ...
! **************************************************************************************************
   SUBROUTINE release_eri_type(eri_env)
      TYPE(eri_type)                                     :: eri_env

      CHARACTER(len=*), PARAMETER :: routineN = 'release_eri_type', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      IF (ASSOCIATED(eri_env%eri)) THEN

         DO i = 1, SIZE(eri_env%eri)
            CALL csr_destroy(eri_env%eri(i)%csr_mat)
            DEALLOCATE (eri_env%eri(i)%csr_mat)
         END DO
         DEALLOCATE (eri_env%eri)

      END IF

   END SUBROUTINE release_eri_type

! **************************************************************************************************

END MODULE qs_active_space_types
