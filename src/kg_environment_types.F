!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004 CP2K developers group                           !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/kg_environment_types [1.0] *
!!
!!   NAME
!!     kg_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gt  SEPT-24-02 
!!
!!   MODIFICATION HISTORY
!!      gt 13-11-03 rewritten           
!!      JGH 20-05-04 use QS as polarization method
!!
!!   SOURCE
!******************************************************************************

MODULE kg_environment_types

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_create
  USE cp_control_types,                ONLY: dft_control_release,&
                                             dft_control_retain,&
                                             dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_retain,&
                                             cp_subsys_set,&
                                             cp_subsystem_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_release,&
                                             dyn_coeff_set_retain,&
                                             dyn_coeff_set_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_retain,&
                                             section_vals_release
  USE kg_energy_types,                 ONLY: deallocate_kg_energy,&
                                             kg_energy_type
  USE kg_gpw_fm_mol_types,             ONLY: kg_fm_mol_set_release,&
                                             kg_fm_mol_set_retain,&
                                             kg_fm_p_type 
  USE kg_force_types,                  ONLY: deallocate_kg_force,&
                                             kg_force_type
  USE kg_rho_types,                    ONLY: kg_rho_release,&
                                             kg_rho_retain,&
                                             kg_rho_type
  USE kg_rspw_types,                   ONLY: kg_rspw_release,&
                                             kg_rspw_retain,&
                                             kg_rspw_type
  USE kg_scf_types,                    ONLY: kg_scf_env_type,&
                                             kg_scf_env_retain,&
                                             kg_scf_env_release
  USE kinds,                           ONLY: dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: qs_environment_type,&
                                             qs_env_retain,&
                                             qs_env_release
  USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                             neighbor_list_set_p_type
  USE scf_control_types,               ONLY: scf_c_release,&
                                             scf_c_retain,&
                                             scf_control_type
  USE simulation_cell,                 ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type
  USE virial_types,                    ONLY: virial_type,&
                                             virial_retain,&
                                             virial_release,&
					     virial_create
  PRIVATE

  character(len=*), private, parameter :: moduleN='kg_environment_types'
  INTEGER, PRIVATE, SAVE :: last_kg_env_id_nr=0

  TYPE kg_environment_type
    INTEGER :: id_nr, ref_count
    LOGICAL :: use_ref_cell
    TYPE(cell_type), POINTER                              :: cell
    TYPE(cell_type), POINTER                              :: cell_ref
    TYPE(virial_type), POINTER                            :: virial
    TYPE(cp_blacs_env_type), POINTER                      :: blacs_env
    TYPE(cp_para_env_type), POINTER                       :: para_env
    TYPE(dft_control_type), POINTER                       :: dft_control
    TYPE(distribution_2d_type), POINTER                   :: distribution_2d
    TYPE(dyn_coeff_set_type), POINTER                     :: dyn_coeff_set
    TYPE(kg_energy_type), POINTER                         :: energy
    TYPE(kg_fm_p_type), POINTER                           :: kg_fm_set
    TYPE(kg_force_type), DIMENSION(:), POINTER            :: force
    TYPE(kg_rho_type),POINTER                             :: rho 
    TYPE(kg_rspw_type),POINTER                            :: rspw 
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sac_ppl 
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sac_ppl_aux 
    TYPE(kg_scf_env_type), POINTER                        :: scf_env
    TYPE(scf_control_type),POINTER                        :: scf_control
    TYPE(cp_subsystem_type), POINTER                      :: subsys
    TYPE(qs_environment_type), POINTER                    :: sub_qs_env
    TYPE(section_vals_type), POINTER                      :: input
  END TYPE kg_environment_type

! *** Public data types ***

  PUBLIC :: kg_environment_type 

! *** Public subroutines ***

  PUBLIC :: get_kg_env,&
            init_kg_env,&
            set_kg_env, &
            kg_env_create, &
            kg_env_retain, &
            kg_env_release
!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_kg_env( kg_env,atomic_kind_set,cell,cell_ref, &
                         use_ref_cell,sub_qs_env,&
                         dft_control,distribution_2d,&
                         dyn_coeff_set,energy,force,input,&
                         kg_fm_set,nao_global,&
                         nelectron_global,maxocc_global,&
                         local_particles,local_molecules,&
                         molecule_kind_set,molecule_set,id_nr,para_env,&
                         blacs_env,particle_set,rho,rspw,scf_env,&
                         sac_ppl,sac_ppl_aux,scf_control,subsys,&
                         virial, error)

!   Purpose: Get the KG environment.

!   ***************************************************************************

    TYPE(kg_environment_type), POINTER         :: kg_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: atomic_kind_set
    TYPE(cell_type), OPTIONAL, POINTER         :: cell, cell_ref
    LOGICAL, OPTIONAL                          :: use_ref_cell
    TYPE(dft_control_type), OPTIONAL, &
      POINTER                                  :: dft_control
    TYPE(qs_environment_type), &
      OPTIONAL, POINTER                        :: sub_qs_env
    TYPE(dyn_coeff_set_type), OPTIONAL, &
      POINTER                                  :: dyn_coeff_set
    TYPE(kg_energy_type), OPTIONAL, POINTER    :: energy
    TYPE(kg_force_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: force
    TYPE(kg_fm_p_type), OPTIONAL, POINTER      :: kg_fm_set
    TYPE(section_vals_type), POINTER, optional :: input
    INTEGER, INTENT(OUT), OPTIONAL             :: nao_global
    INTEGER, INTENT(OUT), OPTIONAL             :: nelectron_global
    REAL(dp), DIMENSION(2), INTENT(OUT), &
      OPTIONAL                                 :: maxocc_global
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                  :: local_particles, &
                                                  local_molecules
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                  :: distribution_2d
    TYPE(molecule_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: molecule_set
    INTEGER, INTENT(out), OPTIONAL             :: id_nr
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                  :: para_env
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                  :: blacs_env
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: particle_set
    TYPE(kg_rho_type), OPTIONAL, POINTER       :: rho
    TYPE(kg_rspw_type), OPTIONAL, POINTER      :: rspw
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER          :: sac_ppl, sac_ppl_aux
    TYPE(scf_control_type), OPTIONAL, &
      POINTER                                  :: scf_control
    TYPE(kg_scf_env_type), OPTIONAL, &
      POINTER                                  :: scf_env
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                  :: subsys
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_kg_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                    :: failure
    TYPE(atomic_kind_list_type), POINTER       :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER      :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER           :: molecules_new
    TYPE(particle_list_type), POINTER          :: particles

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(atomic_kinds, particles,molecule_kinds_new, molecules_new) 
!    CPPrecondition(ASSOCIATED(kg_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(kg_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(kg_env%subsys),cp_failure_level,routineP,error,failure)

    IF (PRESENT(virial)) virial => kg_env%virial
    IF (PRESENT(blacs_env)) blacs_env => kg_env%blacs_env
    IF (PRESENT(cell)) cell => kg_env%cell
    IF (PRESENT(cell_ref)) cell_ref => kg_env%cell_ref
    IF (PRESENT(use_ref_cell)) use_ref_cell = kg_env%use_ref_cell
    IF (PRESENT(dft_control)) dft_control => kg_env % dft_control
    IF (PRESENT(distribution_2d)) distribution_2d => kg_env%distribution_2d
    IF (PRESENT(sub_qs_env)) sub_qs_env => kg_env % sub_qs_env
    IF (PRESENT(dyn_coeff_set)) dyn_coeff_set => kg_env % dyn_coeff_set
    IF (PRESENT(energy)) energy => kg_env%energy
    IF (PRESENT(force)) force => kg_env%force
    IF (PRESENT(id_nr)) id_nr=kg_env%id_nr
    IF (PRESENT(kg_fm_set))  kg_fm_set => kg_env%kg_fm_set
    IF (PRESENT(nao_global)) THEN
        CPPrecondition(ASSOCIATED(kg_env%kg_fm_set),cp_failure_level,routineP,error,failure)
        nao_global = kg_env%kg_fm_set%nao_global
    END IF
    IF (PRESENT(nelectron_global)) THEN
        CPPrecondition(ASSOCIATED(kg_env%kg_fm_set),cp_failure_level,routineP,error,failure)
        nelectron_global = kg_env%kg_fm_set%nelectron_global
    END IF
    IF (PRESENT(maxocc_global)) THEN
        CPPrecondition(ASSOCIATED(kg_env%kg_fm_set),cp_failure_level,routineP,error,failure)
        maxocc_global = kg_env%kg_fm_set%maxocc_global
    END IF
    IF (PRESENT(subsys)) subsys => kg_env%subsys
    IF (PRESENT(para_env)) para_env => kg_env%para_env
    IF (PRESENT(scf_env)) scf_env => kg_env%scf_env
    IF (PRESENT(rho)) rho => kg_env%rho
    IF (PRESENT(rspw)) rspw => kg_env%rspw
    IF (PRESENT(sac_ppl)) sac_ppl => kg_env%sac_ppl
    IF (PRESENT(sac_ppl_aux)) sac_ppl_aux => kg_env%sac_ppl_aux
    IF (PRESENT(scf_control)) scf_control => kg_env%scf_control
    CALL cp_subsys_get(kg_env%subsys,&
                      atomic_kinds=atomic_kinds,&
                      local_molecules_new=local_molecules,&
                      local_particles=local_particles,&
                      particles=particles,&
                      molecule_kinds_new=molecule_kinds_new,&
                      molecules_new=molecules_new,&
                      error=error)
    IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
    IF (PRESENT(particle_set)) particle_set => particles%els
    IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds_new%els
    IF (PRESENT(molecule_set)) molecule_set => molecules_new%els
    IF (PRESENT(input)) input => kg_env%input

  END SUBROUTINE get_kg_env

! *****************************************************************************

  SUBROUTINE init_kg_env ( kg_env, para_env, error )

!   Purpose: Initialise the KG environment.

!   ***************************************************************************

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

!   ---------------------------------------------------------------------------

    NULLIFY (kg_env%blacs_env)
    NULLIFY (kg_env%cell)
    NULLIFY (kg_env%cell_ref)
    kg_env%use_ref_cell = .FALSE.
    NULLIFY (kg_env%dft_control)
    NULLIFY (kg_env%distribution_2d)
    NULLIFY (kg_env%dyn_coeff_set)
    NULLIFY (kg_env%energy)
    NULLIFY (kg_env%input)
    NULLIFY (kg_env%force)
    NULLIFY (kg_env%kg_fm_set)
    NULLIFY (kg_env%rho)
    NULLIFY (kg_env%rspw)
    NULLIFY (kg_env%sac_ppl)
    NULLIFY (kg_env%sac_ppl_aux)
    NULLIFY (kg_env%scf_control)
    NULLIFY (kg_env%scf_env)
    NULLIFY (kg_env%subsys)
    NULLIFY (kg_env%sub_qs_env)
    NULLIFY (kg_env%virial)

    CALL cp_subsys_create(kg_env%subsys, para_env=para_env, error=error)
    CALL cp_para_env_retain(para_env,error=error)
    kg_env%para_env => para_env
    CALL cp_blacs_env_create(kg_env%blacs_env,kg_env%para_env,error=error)

    kg_env%ref_count=1
    last_kg_env_id_nr=last_kg_env_id_nr+1
    kg_env%id_nr=last_kg_env_id_nr

  END SUBROUTINE init_kg_env

! *****************************************************************************

  SUBROUTINE set_kg_env( kg_env,atomic_kind_set,cell,cell_ref,use_ref_cell,  &
                         dft_control,distribution_2d,dyn_coeff_set,energy,input,force, &
                         kg_fm_set,local_particles,local_molecules,&
                         molecule_kind_set,molecule_set,&
                         sub_qs_env,particle_set,rho,rspw,scf_env,&
                         sac_ppl,sac_ppl_aux,scf_control,subsys,virial,&
                         error)

!   Purpose: set the KG environment.

!   ***************************************************************************

    TYPE(kg_environment_type), &
      POINTER                                  :: kg_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: atomic_kind_set
    TYPE(cell_type), OPTIONAL, POINTER         :: cell
    TYPE(cell_type), OPTIONAL, POINTER         :: cell_ref
    LOGICAL, OPTIONAL                          :: use_ref_cell
    TYPE(dft_control_type), OPTIONAL, &
      POINTER                                  :: dft_control
    TYPE(dyn_coeff_set_type), OPTIONAL, &
      POINTER                                  :: dyn_coeff_set
    TYPE(kg_energy_type), OPTIONAL, POINTER    :: energy
    TYPE(kg_force_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: force
    TYPE(kg_fm_p_type), OPTIONAL, POINTER      :: kg_fm_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                  :: local_particles, &
                                                  local_molecules
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                  :: distribution_2d
    TYPE(molecule_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: molecule_set
    TYPE(qs_environment_type), &
      OPTIONAL, POINTER                        :: sub_qs_env
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: particle_set
    TYPE(kg_rho_type), OPTIONAL, POINTER       :: rho
    TYPE(kg_rspw_type), OPTIONAL, POINTER      :: rspw
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER          :: sac_ppl, sac_ppl_aux
    TYPE(scf_control_type), OPTIONAL, &
      POINTER                                  :: scf_control
    TYPE(kg_scf_env_type), OPTIONAL, &
      POINTER                                  :: scf_env
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                  :: subsys
    TYPE(section_vals_type), POINTER, optional :: input    
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routineN = "set_kg_env", &
      routineP = moduleN//':'//routineN

    LOGICAL                                    :: failure
    TYPE(atomic_kind_list_type), POINTER       :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER      :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER           :: molecules_new
    TYPE(particle_list_type), POINTER          :: particles

!   ---------------------------------------------------------------------------
!   ---------------------------------------------------------------------------

    failure=.FALSE.
!    CPPrecondition(ASSOCIATED(kg_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(kg_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(virial)) THEN
       IF (ASSOCIATED(virial)) THEN
          CALL virial_retain(virial,error=error)
       END IF
       CALL virial_release(kg_env%virial,error=error)
       kg_env%virial => virial
    END IF
    IF (PRESENT(cell)) THEN
       CALL cell_retain(cell, error=error)
       CALL cell_release(kg_env%cell,error=error)
       kg_env%cell => cell
    END IF
    IF (PRESENT(cell_ref)) THEN
       CALL cell_retain(cell_ref, error=error)
       CALL cell_release(kg_env%cell_ref,error=error)
       kg_env%cell_ref => cell_ref
    END IF
    IF (PRESENT(use_ref_cell)) kg_env%use_ref_cell=use_ref_cell
    IF (PRESENT(dft_control)) THEN
       CALL dft_control_retain(dft_control, error=error)
       CALL dft_control_release(kg_env%dft_control, error=error)
       kg_env%dft_control => dft_control
    END IF
    IF (PRESENT(scf_control)) THEN 
       CALL scf_c_retain(scf_control,error=error)
       CALL scf_c_release(kg_env%scf_control,error=error)
       kg_env%scf_control => scf_control
    END IF
    IF (PRESENT(dyn_coeff_set)) THEN
       CALL dyn_coeff_set_retain(dyn_coeff_set, error=error)
       CALL dyn_coeff_set_release(kg_env%dyn_coeff_set, error=error)
       kg_env%dyn_coeff_set => dyn_coeff_set
    END IF
    IF (PRESENT(energy)) kg_env%energy => energy
    IF (PRESENT(force)) kg_env%force => force
    IF (PRESENT(kg_fm_set)) THEN
       CALL kg_fm_mol_set_retain(kg_fm_set, error=error)
       CALL kg_fm_mol_set_release(kg_env%kg_fm_set, error=error)
       kg_env%kg_fm_set => kg_fm_set
    END IF
    IF (PRESENT(rho)) THEN ! accepts also null pointers !
      IF (ASSOCIATED(rho)) CALL kg_rho_retain(rho,error=error)
      CALL kg_rho_release(kg_env%rho,error=error)
      kg_env%rho => rho
    END IF
    IF (PRESENT(rspw)) THEN ! accepts also null pointers !
      IF (ASSOCIATED(rspw)) CALL kg_rspw_retain(rspw,error=error)
      CALL kg_rspw_release(kg_env%rspw,error=error)
      kg_env%rspw => rspw
    END IF
    IF (PRESENT(scf_env)) THEN 
       CALL kg_scf_env_retain(scf_env,error=error)
       CALL kg_scf_env_release(kg_env%scf_env, error=error)
       kg_env%scf_env => scf_env
    END IF
    IF ( PRESENT ( sac_ppl ) ) kg_env % sac_ppl => sac_ppl
    IF ( PRESENT ( sac_ppl_aux ) ) kg_env % sac_ppl_aux => sac_ppl_aux
    IF (PRESENT(distribution_2d)) THEN
       CALL distribution_2d_retain(distribution_2d)
       CALL distribution_2d_release(kg_env%distribution_2d)
       kg_env%distribution_2d => distribution_2d
    ENDIF
    IF (PRESENT(subsys)) THEN
       CALL cp_subsys_retain(subsys,error=error)
       CALL cp_subsys_release(kg_env%subsys,error=error)
       kg_env%subsys => subsys
    END IF
    IF (PRESENT(atomic_kind_set)) THEN
      CALL atomic_kind_list_create(atomic_kinds,&
                                   els_ptr=atomic_kind_set,&
                                   error=error)
      CALL cp_subsys_set(kg_env%subsys,&
                        atomic_kinds=atomic_kinds,&
                        error=error)
      CALL atomic_kind_list_release(atomic_kinds,error=error)
    END IF
    IF (PRESENT(sub_qs_env)) THEN
       CALL qs_env_retain(sub_qs_env,error=error)
       kg_env%sub_qs_env => sub_qs_env
       CALL qs_env_release(sub_qs_env,error=error)
    END IF
    IF (PRESENT(particle_set)) THEN
      CALL particle_list_create(particles,&
                                els_ptr=particle_set,&
                                error=error)
      CALL cp_subsys_set(kg_env%subsys,&
                        particles=particles,&
                        error=error)
      CALL particle_list_release(particles,error=error)
    END IF
    IF (PRESENT(local_particles)) THEN
      CALL cp_subsys_set(kg_env%subsys,&
                        local_particles=local_particles,&
                        error=error)
    END IF
    IF (PRESENT(local_molecules)) THEN
      CALL cp_subsys_set(kg_env%subsys,&
                        local_molecules_new=local_molecules,&
                        error=error)
    END IF

    IF (PRESENT(molecule_kind_set)) THEN
      CALL mol_kind_new_list_create(molecule_kinds_new,&
                                    els_ptr=molecule_kind_set,&
                                    error=error)
      CALL cp_subsys_set(kg_env%subsys,&
                        molecule_kinds_new=molecule_kinds_new,&
                        error=error)
      CALL mol_kind_new_list_release(molecule_kinds_new,error=error)
    END IF
    IF (PRESENT(molecule_set)) THEN
      CALL mol_new_list_create(molecules_new,&
                               els_ptr=molecule_set,&
                               error=error)
      CALL cp_subsys_set(kg_env%subsys,&
                        molecules_new=molecules_new,&
                        error=error)
      CALL mol_new_list_release(molecules_new,error=error)
    END IF
    IF (PRESENT(input)) THEN
       CALL section_vals_retain(input,error=error)
       CALL section_vals_release(kg_env%input,error=error)
       kg_env%input => input
    END IF

  END SUBROUTINE set_kg_env

!***************************************************************************
! *****************************************************************************

!!****f* kg_environment_types/kg_env_create [1.0] *
!!
!!   NAME
!!     kg_env_create
!!
!!   SYNOPSIS
!!     Subroutine kg_env_create(kg_env, para_env,  error)
!!       Type(kg_environment_type), Pointer:: kg_env
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine kg_env_create
!!
!!   FUNCTION
!!     allocates and intitializes a kg_env
!!
!!   ARGUMENTS
!!     - kg_env: the object to create
!!     - para_env: the parallel environement for the kg_env
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!*** **********************************************************************
SUBROUTINE kg_env_create(kg_env,para_env,error)
  TYPE(kg_environment_type), POINTER             :: kg_env
  TYPE(cp_para_env_type), POINTER                :: para_env
  TYPE(cp_error_type), INTENT(inout), OPTIONAL   :: error

  CHARACTER(len=*), PARAMETER :: routineN = 'kg_env_create', &
      routineP = moduleN//':'//routineN

  INTEGER                                        :: stat
  LOGICAL                                        :: failure

  failure=.FALSE.

  ALLOCATE(kg_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL init_kg_env(kg_env,para_env=para_env,error=error)         
  END IF
END SUBROUTINE kg_env_create
!***************************************************************************
!!****f* kg_environment_types/kg_env_retain [1.0] *
!!
!!   NAME
!!     kg_env_retain
!!
!!   SYNOPSIS
!!     Subroutine kg_env_retain(kg_env, error)
!!       Type(kg_environment_type), Pointer:: kg_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine kg_env_retain
!!
!!   FUNCTION
!!     retains the given kg_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - kg_env: the object to retain
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!***************************************************************************
SUBROUTINE kg_env_retain(kg_env,error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(kg_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(kg_env%ref_count>0,cp_failure_level,routineP,error,failure)
     kg_env%ref_count=kg_env%ref_count+1
  END IF
END SUBROUTINE kg_env_retain
!***************************************************************************
!!****f* kg_environment_types/kg_env_release [1.0] *
!!
!!   NAME
!!     kg_env_release
!!
!!   SYNOPSIS
!!     Subroutine kg_env_release(kg_env, error)
!!       Type(kg_environment_type), Pointer:: kg_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine kg_env_release
!!
!!   FUNCTION
!!     releases the given kg_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - kg_env: the object to release
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!*** **********************************************************************
SUBROUTINE kg_env_release(kg_env,error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(kg_env)) THEN
     CPPrecondition(kg_env%ref_count>0,cp_failure_level,routineP,error,failure)
     kg_env%ref_count=kg_env%ref_count-1
     IF (kg_env%ref_count<1) THEN
        CALL cell_release(kg_env%cell)
        CALL cell_release(kg_env%cell_ref)
        CALL dft_control_release(kg_env%dft_control, error=error)
        IF (ASSOCIATED(kg_env%sac_ppl)) THEN
           DO i=1,SIZE(kg_env%sac_ppl)
              CALL deallocate_neighbor_list_set(&
                   kg_env%sac_ppl(i)%neighbor_list_set)
           END DO
           DEALLOCATE(kg_env%sac_ppl,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(kg_env%sac_ppl_aux)) THEN
           DO i=1,SIZE(kg_env%sac_ppl_aux)
              CALL deallocate_neighbor_list_set(&
                   kg_env%sac_ppl_aux(i)%neighbor_list_set)
           END DO
           DEALLOCATE(kg_env%sac_ppl_aux,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(kg_env%energy)) &
             CALL deallocate_kg_energy(kg_env%energy)
        IF (ASSOCIATED(kg_env%force))&
             CALL deallocate_kg_force(kg_env%force)
        IF (ASSOCIATED(kg_env%distribution_2d)) THEN
           CALL distribution_2d_release(kg_env%distribution_2d)
        ENDIF
        IF (ASSOCIATED(kg_env%kg_fm_set)) THEN 
           CALL kg_fm_mol_set_release(kg_env%kg_fm_set)
        END IF
        IF (ASSOCIATED(kg_env%scf_env))  THEN
           CALL kg_scf_env_release(kg_env%scf_env, error=error)
        END IF
        IF(ASSOCIATED(kg_env%virial)) THEN
          CALL virial_release(kg_env%virial,error=error)
        END IF
 
        CALL kg_rho_release(kg_env%rho, error=error)
        CALL kg_rspw_release(kg_env%rspw, error=error)
        CALL dyn_coeff_set_release(kg_env%dyn_coeff_set,error=error)
        CALL scf_c_release(kg_env%scf_control, error=error)
        CALL cp_para_env_release(kg_env%para_env,error=error)
        CALL cp_blacs_env_release(kg_env%blacs_env, error=error)
        CALL cp_subsys_release(kg_env%subsys,error=error)
        CALL section_vals_release(kg_env%input,error=error)
        IF (ASSOCIATED(kg_env%sub_qs_env)) THEN
          CALL qs_env_release(kg_env%sub_qs_env,error=error)
        END IF

        DEALLOCATE(kg_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(kg_env)
END SUBROUTINE kg_env_release

!***************************************************************************

END MODULE kg_environment_types
