!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief types that represent a quickstep subsys
!> \author Ole Schuett
! *****************************************************************************
MODULE qs_subsys_types 
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atprop_types,                    ONLY: atprop_type
  USE cell_types,                      ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type
  USE colvar_types,                    ONLY: colvar_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_retain,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE multipole_types,                 ONLY: multipole_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: deallocate_qs_energy,&
                                             qs_energy_type
  USE qs_force_types,                  ONLY: deallocate_qs_force,&
                                             qs_force_type
  USE qs_kind_types,                   ONLY: deallocate_qs_kind_set,&
                                             qs_kind_type
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_subsys_types'

  PUBLIC :: qs_subsys_type

  PUBLIC :: qs_subsys_retain,&
            qs_subsys_release,&
            qs_subsys_get,&
            qs_subsys_set


  TYPE qs_subsys_type
     PRIVATE
     INTEGER                                               :: ref_count = 1
     INTEGER                                               :: nelectron_total = -1
     INTEGER                                               :: nelectron_spin(2) = -1
     TYPE(cp_subsys_type), POINTER                         :: cp_subsys => Null()
     TYPE(qs_kind_type), DIMENSION(:), POINTER             :: qs_kind_set => Null()
     TYPE(cell_type), POINTER                              :: cell_ref => Null()
     LOGICAL                                               :: use_ref_cell = .FALSE.
     TYPE(qs_energy_type), POINTER                         :: energy => Null()
     TYPE(qs_force_type), DIMENSION(:), POINTER            :: force => Null()
  END TYPE qs_subsys_type


CONTAINS

! *****************************************************************************
!> \brief retains a subsys (see doc/ReferenceCounting.html)
!> \param subsys the subsys to retain
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qs_subsys_retain(subsys)
    TYPE(qs_subsys_type), POINTER            :: subsys

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_subsys_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(ASSOCIATED(subsys))
    CPASSERT(subsys%ref_count>0)
    subsys%ref_count=subsys%ref_count+1
  END SUBROUTINE qs_subsys_retain

! *****************************************************************************
!> \brief releases a subsys (see doc/ReferenceCounting.html)
!> \param subsys the subsys to release
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qs_subsys_release(subsys)
    TYPE(qs_subsys_type), POINTER            :: subsys

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_subsys_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(subsys)) THEN
       CPASSERT(subsys%ref_count>0)
       subsys%ref_count=subsys%ref_count-1
       IF (subsys%ref_count==0) THEN
          CALL cp_subsys_release(subsys%cp_subsys)
          CALL cell_release(subsys%cell_ref)
          IF (ASSOCIATED(subsys%qs_kind_set))&
             CALL deallocate_qs_kind_set(subsys%qs_kind_set)
          IF (ASSOCIATED(subsys%energy))&
             CALL deallocate_qs_energy(subsys%energy)
          IF (ASSOCIATED(subsys%force))&
             CALL deallocate_qs_force(subsys%force)
          DEALLOCATE(subsys)
       END IF
    END IF
  END SUBROUTINE qs_subsys_release

! *****************************************************************************
!> \brief ...
!> \param subsys ...
!> \param atomic_kinds ...
!> \param atomic_kind_set ...
!> \param particles ...
!> \param particle_set ...
!> \param local_particles ...
!> \param molecules_new ...
!> \param molecule_set ...
!> \param molecule_kinds_new ...
!> \param molecule_kind_set ...
!> \param local_molecules_new ...
!> \param para_env ...
!> \param colvar_p ...
!> \param shell_particles ...
!> \param core_particles ...
!> \param gci ...
!> \param multipoles ...
!> \param natom ...
!> \param nparticle ...
!> \param ncore ...
!> \param nshell ...
!> \param nkind ...
!> \param atprop ...
!> \param virial ...
!> \param results ...
!> \param cell ...
!> \param cell_ref ...
!> \param use_ref_cell ...
!> \param energy ...
!> \param force ...
!> \param qs_kind_set ...
!> \param cp_subsys ...
!> \param nelectron_total ...
!> \param nelectron_spin ...
! *****************************************************************************
  SUBROUTINE qs_subsys_get(subsys, atomic_kinds, atomic_kind_set, particles,particle_set,&
                           local_particles, molecules_new, molecule_set,&
                           molecule_kinds_new,molecule_kind_set,&
                           local_molecules_new, para_env, colvar_p,&
                           shell_particles, core_particles, gci, multipoles,&
                           natom, nparticle, ncore, nshell, nkind, atprop, virial, &
                           results, cell, cell_ref, use_ref_cell, energy, force,&
                           qs_kind_set,cp_subsys,nelectron_total,nelectron_spin)
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(atomic_kind_list_type), OPTIONAL, &
      POINTER                                :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles
    TYPE(mol_new_list_type), OPTIONAL, &
      POINTER                                :: molecules_new
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_set
    TYPE(mol_kind_new_list_type), OPTIONAL, &
      POINTER                                :: molecule_kinds_new
    TYPE(molecule_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_kind_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_molecules_new
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(colvar_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: colvar_p
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: shell_particles, &
                                                core_particles
    TYPE(global_constraint_type), OPTIONAL, &
      POINTER                                :: gci
    TYPE(multipole_type), OPTIONAL, POINTER  :: multipoles
    INTEGER, INTENT(out), OPTIONAL           :: natom, nparticle, ncore, &
                                                nshell, nkind
    TYPE(atprop_type), OPTIONAL, POINTER     :: atprop
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(cp_result_type), OPTIONAL, POINTER  :: results
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    LOGICAL, OPTIONAL                        :: use_ref_cell
    TYPE(qs_energy_type), OPTIONAL, POINTER  :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: qs_kind_set
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: cp_subsys
    INTEGER, OPTIONAL                        :: nelectron_total
    INTEGER, DIMENSION(2), OPTIONAL          :: nelectron_spin

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_subsys_get', &
      routineP = moduleN//':'//routineN

    CALL cp_subsys_get(subsys%cp_subsys,&
                       atomic_kinds=atomic_kinds,&
                       atomic_kind_set=atomic_kind_set,&
                       particles=particles,&
                       particle_set=particle_set,&
                       local_particles=local_particles,&
                       molecules_new=molecules_new,&
                       molecule_set=molecule_set,&
                       molecule_kinds_new=molecule_kinds_new,&
                       molecule_kind_set=molecule_kind_set,&
                       local_molecules_new=local_molecules_new,&
                       para_env=para_env,&
                       colvar_p=colvar_p,&
                       shell_particles=shell_particles,&
                       core_particles=core_particles,&
                       gci=gci,&
                       multipoles=multipoles,&
                       natom=natom,&
                       nkind=nkind,&
                       nparticle=nparticle,&
                       ncore=ncore,&
                       nshell=nshell,&
                       atprop=atprop,&
                       virial=virial,&
                       results=results,&
                       cell=cell)

     IF (PRESENT(cell_ref)) cell_ref => subsys%cell_ref
     IF (PRESENT(use_ref_cell)) use_ref_cell = subsys%use_ref_cell
     IF (PRESENT(energy)) energy => subsys%energy
     IF (PRESENT(force)) force => subsys%force
     IF (PRESENT(qs_kind_set)) qs_kind_set => subsys%qs_kind_set
     IF (PRESENT(cp_subsys)) cp_subsys => subsys%cp_subsys
     IF (PRESENT(nelectron_total)) nelectron_total=subsys%nelectron_total
     IF (PRESENT(nelectron_spin))  nelectron_spin=subsys%nelectron_spin
  END SUBROUTINE qs_subsys_get


! *****************************************************************************
!> \brief ...
!> \param subsys ...
!> \param cp_subsys ...
!> \param local_particles ...
!> \param local_molecules_new ...
!> \param cell ...
!> \param cell_ref ...
!> \param use_ref_cell ...
!> \param energy ...
!> \param force ...
!> \param qs_kind_set ...
!> \param nelectron_total ...
!> \param nelectron_spin ...
! *****************************************************************************
SUBROUTINE qs_subsys_set(subsys,cp_subsys,&
                         local_particles,local_molecules_new,cell,&
                         cell_ref, use_ref_cell, energy, force,&
                         qs_kind_set, nelectron_total, nelectron_spin)
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: cp_subsys
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles, &
                                                local_molecules_new
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    LOGICAL, OPTIONAL                        :: use_ref_cell
    TYPE(qs_energy_type), OPTIONAL, POINTER  :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: qs_kind_set
    INTEGER, OPTIONAL                        :: nelectron_total
    INTEGER, DIMENSION(2), OPTIONAL          :: nelectron_spin

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_subsys_set', &
      routineP = moduleN//':'//routineN

    IF (PRESENT(cp_subsys)) THEN
       CALL cp_subsys_retain(cp_subsys)
       CALL cp_subsys_release(subsys%cp_subsys)
       subsys%cp_subsys => cp_subsys
    END IF

    CALL cp_subsys_set(subsys%cp_subsys,&
                       local_particles=local_particles,&
                       local_molecules_new=local_molecules_new,&
                       cell=cell)

    IF (PRESENT(cell_ref)) THEN
       CALL cell_retain(cell_ref)
       CALL cell_release(subsys%cell_ref)
       subsys%cell_ref => cell_ref
    END IF

    IF (PRESENT(use_ref_cell)) subsys%use_ref_cell = use_ref_cell
    IF (PRESENT(energy)) subsys%energy => energy
    ! if intels checking (-C) complains here, you have rediscovered a bug in the intel
    ! compiler (present in at least 10.0.025). A testcase has been submitted to intel.
    IF (PRESENT(force)) subsys%force => force
    IF (PRESENT(qs_kind_set)) subsys%qs_kind_set => qs_kind_set
    IF (PRESENT(nelectron_total)) subsys%nelectron_total = nelectron_total
    IF (PRESENT(nelectron_spin)) subsys%nelectron_spin = nelectron_spin
END SUBROUTINE qs_subsys_set

END MODULE qs_subsys_types
