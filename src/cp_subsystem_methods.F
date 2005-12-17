!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_subsystem_methods [1.0] *
!!
!!   NAME
!!     cp_subsystem_methods
!!
!!   FUNCTION
!!     Initialize a small environment for a particular calculation
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     5.2004 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_subsystem_methods
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_get,&
                                             cp_subsys_set,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_type
  USE global_types,                    ONLY: global_environment_type
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
  USE simulation_cell,                 ONLY: cell_retain,&
                                             cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: connectivity_control
  USE topology_types,                  ONLY: deallocate_topology,&
                                             init_topology,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: topology_connectivity_pack,&
                                             topology_coordinate_pack
  USE periodic_table,                  ONLY: ptable
  USE string_utilities,                ONLY: uppercase
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_get_subs_vals,&
                                             section_vals_get,&
                                             section_vals_val_get,&
                                             section_vals_release

  IMPLICIT NONE 
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_subsystem_methods'

  PUBLIC :: create_small_subsystem

CONTAINS

!!****f* cp_subsystem_methods/create_small_subsystem [1.0] *
!!
!!   NAME
!!     create_small_subsystem
!!
!!   FUNCTION
!!     updates the molecule information of the given subsys
!!
!!   NOTES
!!     not really ready to be used with different para_envs for the small
!!     and big part
!!
!!   INPUTS
!!     - small_subsys: the subsystem to create
!!     - big_subsys: the superset of small_subsys
!!     - small_para_env: the parallel environment for the new (small)
!!       subsys
!!     - sub_atom_index: indexes of the atoms that should be in small_subsys
!!     - sub_subsys_index: indexes of the subsys of big_subsys from which
!!       the corresponding atom should be taken
!!     - globenv: place where to find the input file
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed, Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE create_small_subsystem(small_subsystem,big_subsys,small_cell,&
       small_para_env,sub_atom_index,sub_subsys_index,sub_atom_kind_name,&
       globenv,error)
    TYPE(cp_subsystem_type), POINTER         :: small_subsystem
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: big_subsys
    TYPE(cell_type), POINTER                 :: small_cell
    TYPE(cp_para_env_type), POINTER          :: small_para_env
    INTEGER, DIMENSION(:), INTENT(in)        :: sub_atom_index, &
                                                sub_subsys_index
    CHARACTER(len=default_string_length), &
      DIMENSION(:), INTENT(in)               :: sub_atom_kind_name
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_small_subsystem', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, id_, isubsys, nat, stat
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
    CHARACTER(len=default_string_length)     :: my_element

    failure=.FALSE.
    NULLIFY(mol_kinds,mols,particles,atomic_kinds, atomic_kind_set, particle_set, &
         molecule_kind_set,molecule_set,particles,atomic_kinds)

    CPPrecondition(.NOT.ASSOCIATED(small_subsystem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(big_subsys),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
    CALL cp_assert(big_subsys(1)%subsys%para_env%group==small_para_env%group, &
                   cp_failure_level,cp_assertion_failed,routineP, &
                   "big_subsys(1)%subsys%para_env%group==small_para_env%group", &
                   error,failure)

    IF (.NOT. failure) THEN
       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 1. Initialize the topology structure type
       !-----------------------------------------------------------------------------
       CALL init_topology(topology)


       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 2. Get the cell info
       !-----------------------------------------------------------------------------
       topology%cell=> small_cell
       CALL cell_retain(small_cell,error=error)

       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 3. Read in the topology section in the input file if any
       !-----------------------------------------------------------------------------

       ! initialize atom coords from the bigger system
       nat=SIZE(sub_atom_index)
       topology%natoms=nat
       CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%r),cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%label_atmname),cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%label_molname),cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%label_resname),cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%atm_mass),cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%atm_charge),cp_failure_level,routineP,error,failure)
       ALLOCATE(topology%atom_info%r(3,nat),topology%atom_info%label_atmname(nat),&
            topology%atom_info%label_molname(nat),topology%atom_info%label_resname(nat),&
            topology%atom_info%element(nat),&
            topology%atom_info%atm_mass(nat), topology%atom_info%atm_charge(nat),stat=stat)
       isubsys=-1
       DO iat=1,nat
          IF (isubsys/=sub_subsys_index(iat)) THEN
             isubsys=sub_subsys_index(iat)
             CALL cp_subsys_get(big_subsys(isubsys)%subsys,&
                  particles=particles, error=error)
          END IF
          topology%atom_info%r(:,iat)=particles%els(sub_atom_index(iat))%r
          topology%atom_info%label_atmname(iat)=sub_atom_kind_name(iat)
          topology%atom_info%label_molname(iat)=topology%atom_info%label_atmname(iat)
          topology%atom_info%label_resname(iat)=topology%atom_info%label_atmname(iat)
          !
          ! Defining element
          !
          id_ = INDEX(topology%atom_info%label_atmname(iat),"_")-1
          IF ( id_ == -1 ) id_ = LEN_TRIM(topology%atom_info%label_atmname(iat))
          CALL check_subsys_element(topology%atom_info%label_atmname(iat)(1:id_),my_element,globenv,error)
          topology%atom_info%element(iat)=my_element 
          topology%atom_info%atm_mass(iat)=0._dp
          topology%atom_info%atm_charge(iat)=0._dp
       END DO
       topology%bondparm_factor=0._dp


       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 5. Read in or generate the molecular connectivity
       !-----------------------------------------------------------------------------
       CALL connectivity_control(topology,globenv)


       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 7. Pack everything into the molecular types
       !-----------------------------------------------------------------------------
       CALL topology_connectivity_pack(molecule_kind_set,molecule_set,&
            topology,globenv)


       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 8. Pack everything into the atomic types
       !-----------------------------------------------------------------------------
       CALL topology_coordinate_pack(particle_set,atomic_kind_set,&
            molecule_kind_set,molecule_set,topology,globenv)

       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 12. Cleanup the topology structure type
       !-----------------------------------------------------------------------------
       CALL deallocate_topology(topology)


       CALL cp_subsys_create(small_subsystem,para_env=small_para_env,error=error)
       CALL particle_list_create(particles,els_ptr=particle_set,error=error)
       CALL atomic_kind_list_create(atomic_kinds,els_ptr=atomic_kind_set,error=error)
       CALL mol_new_list_create(mols,els_ptr=molecule_set,error=error)
       CALL mol_kind_new_list_create(mol_kinds,els_ptr=molecule_kind_set,error=error)
       CALL cp_subsys_set(small_subsystem,particles=particles,atomic_kinds=atomic_kinds,&
            molecules_new=mols,molecule_kinds_new=mol_kinds,error=error)
       CALL particle_list_release(particles,error=error)
       CALL atomic_kind_list_release(atomic_kinds,error=error)
       CALL mol_new_list_release(mols,error=error)
       CALL mol_kind_new_list_release(mol_kinds,error=error)
    END IF
  END SUBROUTINE create_small_subsystem

!!****f* cp_subsystem_methods/check_subsys_element [1.0] *
!!
!!   NAME
!!     check_subsys_element
!!
!!   FUNCTION
!!     Check and returns the ELEMENT label
!!
!!   NOTES
!!
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     12.2005 created [teo]
!!
!!*** **********************************************************************
  SUBROUTINE check_subsys_element(element_in, element_out, globenv, error)
    IMPLICIT NONE
    ! Arguments
    CHARACTER(len=*), INTENT(IN)                          :: element_in
    CHARACTER(len=default_string_length), INTENT(OUT)     :: element_out
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'check_subsys_element', &
      routineP = moduleN//':'//routineN    
    CHARACTER(len=default_string_length)     :: element_symbol, current_symbol, keyword
    INTEGER :: i, n_rep, i_rep
    LOGICAL :: found, failure
    TYPE(section_vals_type), pointer :: kind_section
    
    found = .FALSE.
    failure = .FALSE.
    element_symbol = element_in
    CALL uppercase(element_symbol)
    DO i=LBOUND(ptable,1),UBOUND(ptable,1)
       current_symbol = ptable(i)%symbol
       CALL uppercase(current_symbol)
       IF (TRIM(current_symbol) == TRIM(element_symbol)) THEN
          element_out = TRIM(element_symbol)
          found = .TRUE.
          EXIT
       END IF
    END DO
    IF (.NOT.found) THEN
       kind_section => section_vals_get_subs_vals(globenv%input_file,&
            "FORCE_EVAL%SUBSYS%KIND",error=error)
       CALL section_vals_get(kind_section,n_repetition=n_rep,error=error)
       DO i_rep=1,n_rep
          CALL section_vals_val_get(kind_section,"_SECTION_PARAMETERS_",&
               c_val=keyword,i_rep_section=i_rep,error=error)
          CALL uppercase(keyword)
          IF (TRIM(keyword)==TRIM(element_symbol)) THEN
             CALL section_vals_val_get(kind_section,i_rep_section=i_rep,&
                  keyword_name="ELEMENT",n_rep_val=i,error=error)
             IF (i>0) THEN
                CALL section_vals_val_get(kind_section,i_rep_section=i_rep,&
                     keyword_name="ELEMENT",c_val=element_symbol,error=error)
                element_out = TRIM(element_symbol)
                found = .TRUE.                
                EXIT
             END IF
          END IF
       END DO
    END IF
    IF (.NOT.found) THEN
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)       
    END IF
  END SUBROUTINE check_subsys_element

END MODULE cp_subsystem_methods
