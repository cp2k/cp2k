!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types needed for a for a Kim-Gordon-like partitioning into molecular 
!>        subunits 
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! *****************************************************************************
MODULE kg_environment_types
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_type
  USE molecule_types_new,              ONLY: molecule_type
  USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                             neighbor_list_set_p_type
  USE task_list_types,                 ONLY: deallocate_task_list,&
                                             task_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_environment_types'

  PUBLIC :: kg_environment_type, kg_env_release, block_list_type, subset_type

  ! all blocks belonging to one molecule
  TYPE block_list_type
    INTEGER, DIMENSION(:), ALLOCATABLE :: rows, cols
  END TYPE block_list_type

  TYPE subset_type
    INTEGER, DIMENSION(:), POINTER :: molecules
    TYPE (neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb 
    TYPE (task_list_type), POINTER :: task_list
  END TYPE subset_type

! *****************************************************************************
!> \brief Contains all the info needed for KG runs... 
!> \param xc_section_kg:    XC section with only the KE functional
!> \param molecule_set:     set of molecular entities as in qs_env
!> \param molecular_blocks: ith entry contains all sparse matrix blocks that
!>                          belong to the ith molecule
!> \param sab_orb_full:     the full neighborlist (build with molecular=.FALSE.)
!> \param subset_of_mol:    ith entry contains the index of the subset, the ith
!>                          molecule belongs to
!> \param subset: the list of molecules, the task list and the neighbor list of
!>                each subset of molecules
!> \param nsubsets: number of subsets
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! *****************************************************************************  
  TYPE kg_environment_type
    INTEGER                                       :: nspins
    INTEGER                                       :: natom
    TYPE(section_vals_type), POINTER              :: xc_section_kg
    INTEGER, ALLOCATABLE, DIMENSION(:)            :: atom_to_molecule 
    TYPE(molecule_type), DIMENSION(:), POINTER    :: molecule_set
    TYPE(block_list_type), DIMENSION(:), POINTER  :: molecular_blocks
    TYPE(neighbor_list_set_p_type),&
      DIMENSION(:), POINTER                       :: sab_orb_full
    INTEGER, DIMENSION(:), POINTER                :: subset_of_mol
    TYPE(subset_type), DIMENSION(:), POINTER      :: subset
    INTEGER                                       :: nsubsets
    INTEGER                                       :: maxdegree
  END TYPE kg_environment_type

CONTAINS

  SUBROUTINE kg_env_release(kg_env, error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iab, imol, isub, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN, handle)

    CPPostcondition(ASSOCIATED(kg_env),cp_failure_level,routineP,error,failure)

    DO imol=1,SIZE(kg_env%molecular_blocks)

      DEALLOCATE(kg_env%molecular_blocks(imol)%rows)
      DEALLOCATE(kg_env%molecular_blocks(imol)%cols)

    END DO

    DEALLOCATE(kg_env%molecular_blocks)

    CALL section_vals_release(kg_env%xc_section_kg, error=error)

    DO iab=1,SIZE(kg_env%sab_orb_full)
      CALL deallocate_neighbor_list_set(kg_env%sab_orb_full(iab)%neighbor_list_set)
    END DO

    DEALLOCATE(kg_env%sab_orb_full)

    DO isub=1,kg_env%nsubsets

      DO iab=1,SIZE(kg_env%subset(isub)%sab_orb)

        CALL deallocate_neighbor_list_set(kg_env%subset(isub)%sab_orb(iab)%neighbor_list_set)
      END DO

      DEALLOCATE(kg_env%subset(isub)%sab_orb)

      CALL deallocate_task_list(kg_env%subset(isub)%task_list, error)

      DEALLOCATE(kg_env%subset(isub)%molecules)

    END DO

    DEALLOCATE(kg_env%subset_of_mol)
    DEALLOCATE(kg_env%subset)

    DEALLOCATE(kg_env%atom_to_molecule) 

    DEALLOCATE(kg_env,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE kg_env_release

END MODULE kg_environment_types
