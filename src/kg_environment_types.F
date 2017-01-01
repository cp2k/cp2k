!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for a for a Kim-Gordon-like partitioning into molecular
!>        subunits
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! **************************************************************************************************
MODULE kg_environment_types
   USE dbcsr_api,                       ONLY: dbcsr_deallocate_matrix_set,&
                                              dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_release,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE molecule_types_new,              ONLY: molecule_type
   USE qs_dispersion_types,             ONLY: qs_dispersion_release,&
                                              qs_dispersion_type
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              neighbor_list_set_p_type
   USE task_list_types,                 ONLY: deallocate_task_list,&
                                              task_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_environment_types'

   PUBLIC :: kg_environment_type, kg_env_release, energy_correction_type

   TYPE subset_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
      TYPE(task_list_type), POINTER :: task_list
   END TYPE subset_type

! *****************************************************************************
!> \brief Contains information on the energy correction functional for KG
!> \par History
!>       03.2014 created
!> \author JGH
! *****************************************************************************
   TYPE energy_correction_type
      CHARACTER(len=20)                                :: ec_name
      INTEGER                                          :: energy_functional
      INTEGER                                          :: ks_solver
      INTEGER                                          :: factorization
      REAL(KIND=dp)                                    :: eps_default
      ! basis set
      CHARACTER(len=20)                                :: basis
      LOGICAL                                          :: mao
      INTEGER                                          :: mao_max_iter
      REAL(KIND=dp)                                    :: mao_eps_grad
      ! energy components
      REAL(KIND=dp)                                    :: etotal
      REAL(KIND=dp)                                    :: eband, exc, ehartree, vhxc
      REAL(KIND=dp)                                    :: edispersion
      ! full neighbor lists and corresponding task list
      TYPE(neighbor_list_set_p_type), &
         DIMENSION(:), POINTER                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(task_list_type), POINTER                    :: task_list
      ! the XC function to be used for the correction, dispersion info
      TYPE(section_vals_type), POINTER                 :: xc_section
      TYPE(qs_dispersion_type), POINTER                :: dispersion_env
      ! matrices in complete basis
      ! KS: Kohn-Sham; H: Core; S: overlap; T: kinetic energy;
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_ks
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_h
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_t
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_p
      ! reduce basis
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER     :: mao_coef
   END TYPE energy_correction_type

! **************************************************************************************************
!> \brief Contains all the info needed for KG runs...
!> \param xc_section_kg:    XC section with only the KE functional
!> \param molecule_set:     set of molecular entities as in qs_env
!> \param sab_orb_full:     full neighborlist (build with molecular=.FALSE.)
!>                          needed for the coloring
!> \param subset_of_mol:    ith entry contains the index of the subset, the ith
!>                          molecule belongs to
!> \param subset:   task list and neighbor list of each subset of molecules
!> \param nsubsets: number of subsets
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! **************************************************************************************************
   TYPE kg_environment_type
      INTEGER                                         :: nspins
      INTEGER                                         :: natom
      TYPE(section_vals_type), POINTER                :: xc_section_kg
      INTEGER, ALLOCATABLE, DIMENSION(:)              :: atom_to_molecule
      TYPE(molecule_type), DIMENSION(:), POINTER      :: molecule_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER    :: tnadd_mat
      TYPE(neighbor_list_set_p_type), &
         DIMENSION(:), POINTER                        :: sab_orb_full, sac_kin
      INTEGER, DIMENSION(:), POINTER                  :: subset_of_mol
      TYPE(subset_type), DIMENSION(:), POINTER        :: subset
      INTEGER                                         :: nsubsets
      INTEGER                                         :: maxdegree
      INTEGER                                         :: coloring_method
      INTEGER                                         :: tnadd_method
      LOGICAL                                         :: energy_correction
      TYPE(energy_correction_type)                    :: ec_env
   END TYPE kg_environment_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param kg_env ...
! **************************************************************************************************
   SUBROUTINE kg_env_release(kg_env)
      TYPE(kg_environment_type), POINTER                 :: kg_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_env_release', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iab, isub

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(kg_env))

      CALL section_vals_release(kg_env%xc_section_kg)

      IF (ASSOCIATED(kg_env%sab_orb_full)) THEN
         DO iab = 1, SIZE(kg_env%sab_orb_full)
            CALL deallocate_neighbor_list_set(kg_env%sab_orb_full(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (kg_env%sab_orb_full)
      ENDIF
      IF (ASSOCIATED(kg_env%sac_kin)) THEN
         DO iab = 1, SIZE(kg_env%sac_kin)
            CALL deallocate_neighbor_list_set(kg_env%sac_kin(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (kg_env%sac_kin)
      END IF

      IF (ASSOCIATED(kg_env%tnadd_mat)) THEN
         CALL dbcsr_deallocate_matrix_set(kg_env%tnadd_mat)
      ENDIF

      DO isub = 1, kg_env%nsubsets
         DO iab = 1, SIZE(kg_env%subset(isub)%sab_orb)
            CALL deallocate_neighbor_list_set(kg_env%subset(isub)%sab_orb(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (kg_env%subset(isub)%sab_orb)
         CALL deallocate_task_list(kg_env%subset(isub)%task_list)
      END DO

      IF (ASSOCIATED(kg_env%subset_of_mol)) DEALLOCATE (kg_env%subset_of_mol)
      IF (ASSOCIATED(kg_env%subset)) DEALLOCATE (kg_env%subset)

      IF (ALLOCATED(kg_env%atom_to_molecule)) DEALLOCATE (kg_env%atom_to_molecule)

      ! energy correction
      IF (kg_env%energy_correction) THEN
         ! neighbor lists
         IF (ASSOCIATED(kg_env%ec_env%sab_orb)) THEN
            DO iab = 1, SIZE(kg_env%ec_env%sab_orb)
               CALL deallocate_neighbor_list_set(kg_env%ec_env%sab_orb(iab)%neighbor_list_set)
            END DO
            DEALLOCATE (kg_env%ec_env%sab_orb)
         END IF
         IF (ASSOCIATED(kg_env%ec_env%sac_ppl)) THEN
            DO iab = 1, SIZE(kg_env%ec_env%sac_ppl)
               CALL deallocate_neighbor_list_set(kg_env%ec_env%sac_ppl(iab)%neighbor_list_set)
            END DO
            DEALLOCATE (kg_env%ec_env%sac_ppl)
         END IF
         IF (ASSOCIATED(kg_env%ec_env%sap_ppnl)) THEN
            DO iab = 1, SIZE(kg_env%ec_env%sap_ppnl)
               CALL deallocate_neighbor_list_set(kg_env%ec_env%sap_ppnl(iab)%neighbor_list_set)
            END DO
            DEALLOCATE (kg_env%ec_env%sap_ppnl)
         END IF
         ! operator matrices
         IF (ASSOCIATED(kg_env%ec_env%matrix_ks)) CALL dbcsr_deallocate_matrix_set(kg_env%ec_env%matrix_ks)
         IF (ASSOCIATED(kg_env%ec_env%matrix_h)) CALL dbcsr_deallocate_matrix_set(kg_env%ec_env%matrix_h)
         IF (ASSOCIATED(kg_env%ec_env%matrix_s)) CALL dbcsr_deallocate_matrix_set(kg_env%ec_env%matrix_s)
         IF (ASSOCIATED(kg_env%ec_env%matrix_t)) CALL dbcsr_deallocate_matrix_set(kg_env%ec_env%matrix_t)
         IF (ASSOCIATED(kg_env%ec_env%matrix_p)) CALL dbcsr_deallocate_matrix_set(kg_env%ec_env%matrix_p)
         IF (ASSOCIATED(kg_env%ec_env%task_list)) THEN
            CALL deallocate_task_list(kg_env%ec_env%task_list)
         END IF
         ! reduced basis
         IF (ASSOCIATED(kg_env%ec_env%mao_coef)) CALL dbcsr_deallocate_matrix_set(kg_env%ec_env%mao_coef)
         ! dispersion environment
         IF (ASSOCIATED(kg_env%ec_env%dispersion_env)) THEN
            CALL qs_dispersion_release(kg_env%ec_env%dispersion_env)
         END IF
      END IF

      DEALLOCATE (kg_env)

      CALL timestop(handle)

   END SUBROUTINE kg_env_release

END MODULE kg_environment_types
