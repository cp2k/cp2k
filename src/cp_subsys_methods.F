!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Initialize a small environment for a particular calculation
!> \par History
!>      5.2004 created [fawzi]
!>      9.2007 cleaned [tlaino] - University of Zurich
!> \author Teodoro Laino
! *****************************************************************************
MODULE cp_subsys_methods
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atprop_types,                    ONLY: atprop_create
  USE cell_types,                      ONLY: cell_retain,&
                                             cell_type
  USE colvar_methods,                  ONLY: colvar_read
  USE cp_para_env,                     ONLY: cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_create
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE exclusion_types,                 ONLY: exclusion_type
  USE input_constants,                 ONLY: do_conn_off,&
                                             do_stress_analytical,&
                                             do_stress_diagonal_anal,&
                                             do_stress_diagonal_numer,&
                                             do_stress_none,&
                                             do_stress_numerical
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
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
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE topology,                        ONLY: connectivity_control,&
                                             topology_control
  USE topology_connectivity_util,      ONLY: topology_connectivity_pack
  USE topology_coordinate_util,        ONLY: topology_coordinate_pack
  USE topology_types,                  ONLY: deallocate_topology,&
                                             init_topology,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: check_subsys_element
  USE virial_types,                    ONLY: virial_create,&
                                             virial_set
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_subsys_methods'

  PUBLIC :: create_small_subsys, cp_subsys_create

CONTAINS

! *****************************************************************************
!> \brief Creates allocates and fills subsys from given input.
!> \param subsys ...
!> \param para_env ...
!> \param root_section ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \param qmmm ...
!> \param qmmm_env ...
!> \param exclusions ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cp_subsys_create(subsys, para_env, &
                              root_section, force_env_section, subsys_section,&
                              use_motion_section, qmmm, qmmm_env, exclusions, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN), OPTIONAL            :: use_motion_section
    LOGICAL, OPTIONAL                        :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(exclusion_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: exclusions
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stress_tensor
    LOGICAL :: atomic_energy, atomic_stress, failure, my_use_motion_section, &
      pv_availability, pv_diagonal, pv_numerical
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(mol_kind_new_list_type), POINTER    :: mol_kinds
    TYPE(mol_new_list_type), POINTER         :: mols
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: colvar_section, &
                                                my_force_env_section, &
                                                my_subsys_section

    CPPostcondition(.NOT.ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
    ALLOCATE(subsys)

    CALL cp_para_env_retain(para_env,error=error)
    subsys%para_env => para_env

    my_use_motion_section = .FALSE.
    IF(PRESENT(use_motion_section)) &
      my_use_motion_section = use_motion_section

    my_force_env_section => section_vals_get_subs_vals(root_section,"FORCE_EVAL",error=error)
    IF(PRESENT(force_env_section)) &
      my_force_env_section => force_env_section

    my_subsys_section => section_vals_get_subs_vals(my_force_env_section,"SUBSYS",error=error)
    IF(PRESENT(subsys_section)) &
      my_subsys_section => subsys_section

    colvar_section => section_vals_get_subs_vals(my_subsys_section,"COLVAR",error=error)

    CALL cp_subsys_read_colvar(subsys, colvar_section, error=error)

    !   *** Read the particle coordinates and allocate the atomic kind, ***
    !   *** the molecule kind, and the molecule data structures         ***
    CALL topology_control(atomic_kind_set, particle_set, molecule_kind_set, molecule_set,&
                       subsys%colvar_p, subsys%gci, root_section, para_env,&
                       force_env_section=my_force_env_section,&
                       subsys_section=my_subsys_section,use_motion_section=my_use_motion_section,&
                       qmmm=qmmm, qmmm_env=qmmm_env, exclusions=exclusions,&
                       error=error)

    CALL particle_list_create(particles,els_ptr=particle_set,error=error)
    CALL atomic_kind_list_create(atomic_kinds,els_ptr=atomic_kind_set,error=error)
    CALL mol_new_list_create(mols,els_ptr=molecule_set,error=error)
    CALL mol_kind_new_list_create(mol_kinds,els_ptr=molecule_kind_set,error=error)

    CALL cp_subsys_set(subsys,particles=particles,atomic_kinds=atomic_kinds,&
         molecules_new=mols,molecule_kinds_new=mol_kinds,error=error)

    CALL particle_list_release(particles,error=error)
    CALL atomic_kind_list_release(atomic_kinds,error=error)
    CALL mol_new_list_release(mols,error=error)
    CALL mol_kind_new_list_release(mol_kinds,error=error)

    ! Should we compute the virial?
    CALL section_vals_val_get(my_force_env_section,"STRESS_TENSOR",i_val=stress_tensor,error=error)
    SELECT CASE(stress_tensor)
    CASE(do_stress_none)
       pv_availability=.FALSE.
       pv_numerical=.FALSE.
       pv_diagonal=.FALSE.
    CASE(do_stress_analytical)
       pv_availability=.TRUE.
       pv_numerical=.FALSE.
       pv_diagonal=.FALSE.
    CASE(do_stress_numerical)
       pv_availability=.TRUE.
       pv_numerical=.TRUE.
       pv_diagonal=.FALSE.
    CASE(do_stress_diagonal_anal)
       pv_availability=.TRUE.
       pv_numerical=.FALSE.
       pv_diagonal=.TRUE.
    CASE(do_stress_diagonal_numer)
       pv_availability=.TRUE.
       pv_numerical=.TRUE.
       pv_diagonal=.TRUE.
    END SELECT

    CALL virial_create(subsys%virial, error=error)
    CALL virial_set(virial=subsys%virial,&
                   pv_availability=pv_availability,&
                   pv_numer=pv_numerical,&
                   pv_diagonal=pv_diagonal)

    ! Should we compute atomic properties?
    CALL atprop_create(subsys%atprop,error)
    CALL section_vals_val_get(my_force_env_section,"PROPERTIES%ATOMIC%ENERGY",l_val=atomic_energy,error=error)
    subsys%atprop%energy = atomic_energy
    CALL section_vals_val_get(my_force_env_section,"PROPERTIES%ATOMIC%PRESSURE",l_val=atomic_stress,error=error)
    IF (atomic_stress) THEN
       CPPrecondition(pv_availability,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.pv_numerical,cp_failure_level,routineP,error,failure)
    END IF
    subsys%atprop%stress = atomic_stress

    CALL cp_result_create(subsys%results,error=error)
  END SUBROUTINE cp_subsys_create

! *****************************************************************************
!> \brief reads the colvar section of the colvar
!> \param subsys ...
!> \param colvar_section ...
!> \param error ...
!> \par History
!>      2006.01 Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_subsys_read_colvar(subsys,colvar_section,error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: colvar_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_read_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ig, ncol, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CALL section_vals_get(colvar_section,n_repetition=ncol,error=error)
    ALLOCATE(subsys%colvar_p(ncol),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ig= 1, ncol
       NULLIFY(subsys%colvar_p(ig)%colvar)
       CALL colvar_read(subsys%colvar_p(ig)%colvar,ig,colvar_section, subsys%para_env, error=error)
    ENDDO
  END SUBROUTINE cp_subsys_read_colvar

! *****************************************************************************
!> \brief updates the molecule information of the given subsys
!> \param small_subsys the subsys to create
!> \param big_subsys the superset of small_subsys
!> \param small_cell ...
!> \param small_para_env the parallel environment for the new (small)
!>        subsys
!> \param sub_atom_index indexes of the atoms that should be in small_subsys
!> \param sub_atom_kind_name ...
!> \param para_env ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param ignore_outside_box ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2004 created [fawzi]
!> \author Fawzi Mohamed, Teodoro Laino
!> \note
!>      not really ready to be used with different para_envs for the small
!>      and big part
! *****************************************************************************
  SUBROUTINE create_small_subsys(small_subsys, big_subsys, small_cell,&
       small_para_env, sub_atom_index, sub_atom_kind_name, &
       para_env, force_env_section, subsys_section, ignore_outside_box, error)

    TYPE(cp_subsys_type), POINTER            :: small_subsys, big_subsys
    TYPE(cell_type), POINTER                 :: small_cell
    TYPE(cp_para_env_type), POINTER          :: small_para_env
    INTEGER, DIMENSION(:), INTENT(in)        :: sub_atom_index
    CHARACTER(len=default_string_length), &
      DIMENSION(:), INTENT(in)               :: sub_atom_kind_name
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(in), OPTIONAL            :: ignore_outside_box
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_small_subsys', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: my_element, strtmp1
    INTEGER                                  :: iat, id_, nat, stat
    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(mol_kind_new_list_type), POINTER    :: mol_kinds
    TYPE(mol_new_list_type), POINTER         :: mols
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(topology_parameters_type)           :: topology

    failure=.FALSE.
    NULLIFY(mol_kinds,mols,particles,atomic_kinds, atomic_kind_set, particle_set, &
         molecule_kind_set,molecule_set,particles,atomic_kinds)

    CPPrecondition(.NOT.ASSOCIATED(small_subsys),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(big_subsys),cp_failure_level,routineP,error,failure)
    CALL cp_assert(big_subsys%para_env%group==small_para_env%group, &
         cp_failure_level, cp_assertion_failed, routineP, &
         "big_subsys%para_env%group==small_para_env%group",error,failure)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Initialize the topology structure type
    !-----------------------------------------------------------------------------
    CALL init_topology(topology,error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Get the cell info
    !-----------------------------------------------------------------------------
    topology%cell=> small_cell
    CALL cell_retain(small_cell,error=error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Initialize atom coords from the bigger system
    !-----------------------------------------------------------------------------
    nat=SIZE(sub_atom_index)
    topology%natoms=nat
    CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%r),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%id_atmname),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%id_molname),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%id_resname),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%atm_mass),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%atm_charge),cp_failure_level,routineP,error,failure)
    ALLOCATE(topology%atom_info%r(3,nat),topology%atom_info%id_atmname(nat),&
         topology%atom_info%id_molname(nat),topology%atom_info%id_resname(nat),&
         topology%atom_info%id_element(nat), topology%atom_info%atm_mass(nat),&
         topology%atom_info%atm_charge(nat),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL cp_subsys_get(big_subsys, particles=particles, error=error)
    DO iat=1,nat
       topology%atom_info%r(:,iat)=particles%els(sub_atom_index(iat))%r
       topology%atom_info%id_atmname(iat)=str2id(s2s(sub_atom_kind_name(iat)))
       topology%atom_info%id_molname(iat)=topology%atom_info%id_atmname(iat)
       topology%atom_info%id_resname(iat)=topology%atom_info%id_atmname(iat)
       !
       ! Defining element
       !
       id_ = INDEX(id2str(topology%atom_info%id_atmname(iat)),"_")-1
       IF ( id_ == -1 ) id_ = LEN_TRIM(id2str(topology%atom_info%id_atmname(iat)))
       strtmp1=id2str(topology%atom_info%id_atmname(iat))
       strtmp1=strtmp1(1:id_)
       CALL check_subsys_element(strtmp1,strtmp1, my_element,&
            subsys_section,use_mm_map_first=.FALSE.,error=error)
       topology%atom_info%id_element(iat)=str2id(s2s(my_element))
       topology%atom_info%atm_mass(iat)=0._dp
       topology%atom_info%atm_charge(iat)=0._dp
    END DO
    topology%conn_type=do_conn_off

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Read in or generate the molecular connectivity
    !-----------------------------------------------------------------------------
    CALL connectivity_control(topology,para_env,subsys_section=subsys_section,&
         force_env_section=force_env_section,error=error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Pack everything into the molecular types
    !-----------------------------------------------------------------------------
    CALL topology_connectivity_pack(molecule_kind_set,molecule_set,&
         topology,subsys_section=subsys_section,error=error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Pack everything into the atomic types
    !-----------------------------------------------------------------------------
    CALL topology_coordinate_pack(particle_set,atomic_kind_set,&
         molecule_kind_set,molecule_set,topology,subsys_section=subsys_section,&
         force_env_section=force_env_section,ignore_outside_box=ignore_outside_box,&
         error=error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. Cleanup the topology structure type
    !-----------------------------------------------------------------------------
    CALL deallocate_topology(topology,error=error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 8. Allocate new subsys
    !-----------------------------------------------------------------------------
    ALLOCATE(small_subsys)
    CALL cp_para_env_retain(para_env,error=error)
    small_subsys%para_env => para_env
    CALL particle_list_create(particles,els_ptr=particle_set,error=error)
    CALL atomic_kind_list_create(atomic_kinds,els_ptr=atomic_kind_set,error=error)
    CALL mol_new_list_create(mols,els_ptr=molecule_set,error=error)
    CALL mol_kind_new_list_create(mol_kinds,els_ptr=molecule_kind_set,error=error)
    CALL cp_subsys_set(small_subsys,particles=particles,atomic_kinds=atomic_kinds,&
         molecules_new=mols,molecule_kinds_new=mol_kinds,error=error)
    CALL particle_list_release(particles,error=error)
    CALL atomic_kind_list_release(atomic_kinds,error=error)
    CALL mol_new_list_release(mols,error=error)
    CALL mol_kind_new_list_release(mol_kinds,error=error)

    CALL virial_create(small_subsys%virial, error=error)
    CALL atprop_create(small_subsys%atprop,error)
    CALL cp_result_create(small_subsys%results,error=error)
  END SUBROUTINE create_small_subsys

END MODULE cp_subsys_methods
