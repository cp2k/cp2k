!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for a for a Kim-Gordon-like partitioning into molecular
!>        subunits
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! **************************************************************************************************
MODULE kg_environment_types
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_type
   USE integration_grid_types,          ONLY: deallocate_intgrid,&
                                              integration_grid_type
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_density_release,&
                                              lri_density_type,&
                                              lri_env_release,&
                                              lri_environment_type
   USE molecule_types,                  ONLY: molecule_type
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_grid_atom,                    ONLY: atom_integration_grid_type,&
                                              deallocate_atom_int_grid
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
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
      INTEGER                                         :: tnadd_method
      TYPE(neighbor_list_set_p_type), &
         DIMENSION(:), POINTER                        :: sab_orb_full, sac_kin
      !
      INTEGER, DIMENSION(:), POINTER                  :: subset_of_mol
      TYPE(subset_type), DIMENSION(:), POINTER        :: subset
      INTEGER                                         :: nsubsets
      INTEGER                                         :: maxdegree
      INTEGER                                         :: coloring_method
      !
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: tnadd_mat
      ! LRI
      TYPE(lri_environment_type), POINTER             :: lri_env
      TYPE(lri_density_type), POINTER                 :: lri_density
      ! atomic grid
      TYPE(atom_integration_grid_type), POINTER       :: int_grid_atom
      TYPE(integration_grid_type), POINTER            :: int_grid_molecules
      TYPE(integration_grid_type), POINTER            :: int_grid_full
   END TYPE kg_environment_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param kg_env ...
! **************************************************************************************************
   SUBROUTINE kg_env_release(kg_env)
      TYPE(kg_environment_type), POINTER                 :: kg_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'kg_env_release'

      INTEGER                                            :: handle, isub

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(kg_env))

      CALL release_neighbor_list_sets(kg_env%sab_orb_full)
      CALL release_neighbor_list_sets(kg_env%sac_kin)

      IF (ASSOCIATED(kg_env%tnadd_mat)) THEN
         CALL dbcsr_deallocate_matrix_set(kg_env%tnadd_mat)
      ENDIF

      DO isub = 1, kg_env%nsubsets
         CALL release_neighbor_list_sets(kg_env%subset(isub)%sab_orb)
         CALL deallocate_task_list(kg_env%subset(isub)%task_list)
      END DO

      IF (ASSOCIATED(kg_env%subset_of_mol)) DEALLOCATE (kg_env%subset_of_mol)
      IF (ASSOCIATED(kg_env%subset)) DEALLOCATE (kg_env%subset)

      IF (ALLOCATED(kg_env%atom_to_molecule)) DEALLOCATE (kg_env%atom_to_molecule)

      ! LRI
      IF (ASSOCIATED(kg_env%lri_env)) THEN
         CALL lri_env_release(kg_env%lri_env)
      END IF
      IF (ASSOCIATED(kg_env%lri_density)) THEN
         CALL lri_density_release(kg_env%lri_density)
      END IF
      ! atom grids
      IF (ASSOCIATED(kg_env%int_grid_atom)) THEN
         CALL deallocate_atom_int_grid(kg_env%int_grid_atom)
      END IF
      IF (ASSOCIATED(kg_env%int_grid_molecules)) THEN
         CALL deallocate_intgrid(kg_env%int_grid_molecules)
      END IF
      IF (ASSOCIATED(kg_env%int_grid_full)) THEN
         CALL deallocate_intgrid(kg_env%int_grid_full)
      END IF

      DEALLOCATE (kg_env)

      CALL timestop(handle)

   END SUBROUTINE kg_env_release

END MODULE kg_environment_types
