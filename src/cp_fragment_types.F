!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_fragment_types [1.0] *
!!
!!   NAME
!!     cp_fragment_types
!!
!!   FUNCTION
!!     types that represent a fragment, i.e. a part of the system
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_fragment_types
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_release,&
                                             atomic_kind_list_retain,&
                                             atomic_kind_list_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_retain,&
                                             distribution_1d_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE mol_kind_list_types,             ONLY: mol_kind_list_release,&
                                             mol_kind_list_retain,&
                                             mol_kind_list_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_release,&
                                             mol_kind_new_list_retain,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_release,&
                                             mol_new_list_retain,&
                                             mol_new_list_type
  USE mol_struct_list_types,           ONLY: mol_struct_list_release,&
                                             mol_struct_list_retain,&
                                             mol_struct_list_type
  USE molecule_types,                  ONLY: particle_node_type
  USE particle_list_types,             ONLY: particle_list_release,&
                                             particle_list_retain,&
                                             particle_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fragment_types'
  INTEGER, PRIVATE, SAVE :: last_fragment_id=0

  PUBLIC :: cp_fragment_type, cp_fragment_p_type
  PUBLIC :: fragment_create, fragment_retain, fragment_release,&
       fragment_get, fragment_set

!***
!****************************************************************************

!!****s* force_control/cp_fragment_type [1.0] *
!!
!!   NAME
!!     cp_fragment_type
!!
!!   FUNCTION
!!     represents a part of a system: atoms, molecules, their pos,vel,...
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - atomic_kinds: list with all the kinds in the actual fragment
!!     - particles: list with the particles of the actual fragment
!!     - local_particles: the particles that are local to the actual processor
!!     - molecule_kinds: list with the molecule kinds
!!     - local_molecules: the molecule structures of the actual fragment
!!       that are local to this processor
!!     - para_env: the parallel environment of the actual fragment
!!     - pnode: *to be removed* !
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_fragment_type
!     PRIVATE !FM removed for convenience, but to as if it where there
     INTEGER :: ref_count, id_nr
     TYPE (atomic_kind_list_type), POINTER :: atomic_kinds
     TYPE (particle_list_type), POINTER :: particles
     TYPE (distribution_1d_type), POINTER :: local_particles
     TYPE (cp_para_env_type), POINTER :: para_env

!old
     TYPE (mol_kind_list_type), POINTER :: molecule_kinds
     TYPE (mol_struct_list_type), POINTER :: local_molecules

!new
     TYPE (mol_new_list_type), POINTER :: molecules_new
     TYPE (mol_kind_new_list_type), POINTER :: molecule_kinds_new
     TYPE (distribution_1d_type), POINTER :: local_molecules_new

!to remove
     TYPE (particle_node_type), DIMENSION(:), POINTER :: pnode
  END TYPE cp_fragment_type
!!***
  !****************************************************************************

!!****s* cp_fragment_types/cp_fragment_p_type [1.0] *
!!
!!   NAME
!!     cp_fragment_p_type
!!
!!   FUNCTION
!!     represent a pointer to a fragment, to be able to create arrays
!!     of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - fragment: the pointer to the fragment
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_fragment_p_type
     TYPE(cp_fragment_type), POINTER :: fragment
  END TYPE cp_fragment_p_type
!!***
  !****************************************************************************

CONTAINS

!!****f* cp_fragment_types/fragment_create [1.0] *
!!
!!   NAME
!!     fragment_create
!!
!!   SYNOPSIS
!!     Subroutine fragment_create(fragment, para_env, error)
!!       Type(cp_fragment_type), Pointer:: fragment
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fragment_create
!!
!!   FUNCTION
!!     creates a fragment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - fragment: the fragment to allocate and initialize
!!     - para_env: the parallel environment of this fragment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE fragment_create(fragment, para_env, error)
  TYPE(cp_fragment_type), POINTER :: fragment
  TYPE(cp_para_env_type), POINTER :: para_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fragment_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.
  
  ALLOCATE(fragment, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     last_fragment_id=last_fragment_id+1
     fragment%id_nr=last_fragment_id
     fragment%ref_count=1
     CALL cp_para_env_retain(para_env,error=error)
     fragment%para_env => para_env
     NULLIFY(fragment%atomic_kinds, fragment%particles, &
          fragment%local_particles, fragment%molecule_kinds,&
          fragment%local_molecules, fragment%molecules_new,&
          fragment%molecule_kinds_new, fragment%local_molecules_new,&
          fragment%pnode)
  END IF
END SUBROUTINE fragment_create
!***************************************************************************

!!****f* cp_fragment_types/fragment_retain [1.0] *
!!
!!   NAME
!!     fragment_retain
!!
!!   SYNOPSIS
!!     Subroutine fragment_retain(fragment, error)
!!       Type(cp_fragment_type), Pointer:: fragment
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fragment_retain
!!
!!   FUNCTION
!!     retains a fragment (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - fragment: the fragment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE fragment_retain(fragment, error)
  TYPE(cp_fragment_type), POINTER :: fragment
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fragment_retain',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(fragment),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(fragment%ref_count>0,cp_failure_level,routineP,error)
     fragment%ref_count=fragment%ref_count+1
  END IF
END SUBROUTINE fragment_retain
!***************************************************************************

!!****f* cp_fragment_types/fragment_release [1.0] *
!!
!!   NAME
!!     fragment_release
!!
!!   SYNOPSIS
!!     Subroutine fragment_release(fragment, error)
!!       Type(cp_fragment_type), Pointer:: fragment
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fragment_release
!!
!!   FUNCTION
!!     releases a fragment (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - fragment: the fragment to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE fragment_release(fragment, error)
  TYPE(cp_fragment_type), POINTER :: fragment
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fragment_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.
  
  IF (ASSOCIATED(fragment)) THEN
     CPPreconditionNoFail(fragment%ref_count>0,cp_failure_level,routineP,error)
     fragment%ref_count=fragment%ref_count-1
     IF (fragment%ref_count==0) THEN
        CALL atomic_kind_list_release(fragment%atomic_kinds,error=error)
        CALL particle_list_release(fragment%particles, error=error)
        CALL distribution_1d_release(fragment%local_particles, error=error)
        CALL mol_kind_list_release(fragment%molecule_kinds, error=error)
        CALL mol_struct_list_release(fragment%local_molecules, error=error)
        CALL mol_kind_new_list_release(fragment%molecule_kinds_new, error=error)
        CALL mol_new_list_release(fragment%molecules_new, error=error)
        CALL distribution_1d_release(fragment%local_molecules_new,error=error)
        CALL cp_para_env_release(fragment%para_env, error=error)
        NULLIFY(fragment%pnode)
        DEALLOCATE(fragment, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(fragment)
END SUBROUTINE fragment_release
!***************************************************************************

!!****f* cp_fragment_types/fragment_set [1.0] *
!!
!!   NAME
!!     fragment_set
!!
!!   SYNOPSIS
!!     Subroutine fragment_set(fragment, atomic_kinds, particles,&
!!         local_particles, molecule_kinds, local_molecules, molecules_new,&
!!         molecule_kinds_new, local_molecules_new, pnode, para_env, error)
!!       Type(cp_fragment_type), Pointer:: fragment
!!       Type(atomic_kind_list_type), Pointer, Optional:: atomic_kinds
!!       Type(particle_list_type), Pointer, Optional:: particles
!!       Type(distribution_1d_type), Pointer, Optional:: local_particles
!!       Type(mol_kind_list_type), Pointer, Optional:: molecule_kinds
!!       Type(mol_struct_list_type), Pointer, Optional:: local_molecules
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Type(mol_new_list_type), Pointer, Optional:: molecules_new
!!       Type(mol_kind_new_list_type), Pointer, Optional::&
!!         molecule_kinds_new
!!       Type(distribution_1d_type), Pointer, Optional::&
!!         local_molecules_new
!!       Type(particle_node_type), Dimension(:), Pointer, Optional:: pnode
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fragment_set
!!
!!   FUNCTION
!!     sets various propreties of the fragment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - fragment: the fragment you want to modify
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!     For the description of the other arguments see cp_fragment_type
!!     attributes.
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE fragment_set(fragment, atomic_kinds, particles, local_particles,&
     molecule_kinds, local_molecules,molecules_new, molecule_kinds_new,&
     local_molecules_new,pnode,para_env, error)
  TYPE (cp_fragment_type), POINTER :: fragment
  TYPE (atomic_kind_list_type), POINTER, OPTIONAL :: atomic_kinds
  TYPE (particle_list_type), POINTER, OPTIONAL :: particles
  TYPE (distribution_1d_type), POINTER, OPTIONAL :: local_particles
  TYPE (mol_kind_list_type), POINTER, OPTIONAL :: molecule_kinds
  TYPE (mol_struct_list_type), POINTER, OPTIONAL :: local_molecules
  TYPE (cp_para_env_type), POINTER, OPTIONAL :: para_env
  TYPE (mol_new_list_type), POINTER, OPTIONAL :: molecules_new
  TYPE (mol_kind_new_list_type), POINTER, OPTIONAL :: molecule_kinds_new
  TYPE (distribution_1d_type), POINTER, OPTIONAL :: local_molecules_new
  TYPE (particle_node_type), DIMENSION(:), POINTER, OPTIONAL :: pnode
  TYPE (cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fragment_set',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(fragment),cp_failure_level,routineP,error,failure)
  CPPrecondition(fragment%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(atomic_kinds)) THEN
        CALL atomic_kind_list_retain(atomic_kinds,error=error)
        CALL atomic_kind_list_release(fragment%atomic_kinds, error=error)
        fragment%atomic_kinds => atomic_kinds
     END IF
     IF (PRESENT(particles)) THEN
        CALL particle_list_retain(particles, error=error)
        CALL particle_list_release(fragment%particles, error=error)
        fragment%particles => particles
     END IF
     IF (PRESENT(local_particles)) THEN
        CALL distribution_1d_retain(local_particles,error=error)
        CALL distribution_1d_release(fragment%local_particles,error=error)
        fragment%local_particles => local_particles
     END IF
     IF (PRESENT(molecule_kinds)) THEN
        CALL mol_kind_list_retain(molecule_kinds, error=error)
        CALL mol_kind_list_release(fragment%molecule_kinds, error=error)
        fragment%molecule_kinds => molecule_kinds
     END IF
     IF (PRESENT(local_molecules)) THEN
        CALL mol_struct_list_retain(local_molecules, error=error)
        CALL mol_struct_list_release(fragment%local_molecules, error=error)
        fragment%local_molecules => local_molecules
     END IF
     IF (PRESENT(local_molecules_new)) THEN
        CALL distribution_1d_retain(local_molecules_new,error=error)
        CALL distribution_1d_release(fragment%local_molecules_new,error=error)
        fragment%local_molecules_new => local_molecules_new
     END IF
     IF (PRESENT(molecule_kinds_new)) THEN
        CALL mol_kind_new_list_retain(molecule_kinds_new, error=error)
        CALL mol_kind_new_list_release(fragment%molecule_kinds_new, error=error)
        fragment%molecule_kinds_new => molecule_kinds_new
     END IF
     IF (PRESENT(molecules_new)) THEN
        CALL mol_new_list_retain(molecules_new, error=error)
        CALL mol_new_list_release(fragment%molecules_new, error=error)
        fragment%molecules_new => molecules_new
     END IF
     IF (PRESENT(para_env)) THEN
        CALL cp_para_env_retain(para_env, error=error)
        CALL cp_para_env_release(fragment%para_env, error=error)
        fragment%para_env => para_env
     END IF
     IF (PRESENT(pnode)) fragment%pnode => pnode
  END IF
END SUBROUTINE fragment_set
!***************************************************************************

!!****f* cp_fragment_types/fragment_get [1.0] *
!!
!!   NAME
!!     fragment_get
!!
!!   SYNOPSIS
!!     Subroutine fragment_get(fragment, id_nr, ref_count, atomic_kinds,&
!!         particles, local_particles, molecule_kinds, local_molecules,&
!!         molecules_new, molecule_kinds_new, local_molecules_new, para_env,&
!!         pnode, error)
!!       Type(cp_fragment_type), Pointer:: fragment
!!       Integer, Intent (OUT), Optional:: id_nr, ref_count
!!       Type(atomic_kind_list_type), Pointer, Optional:: atomic_kinds
!!       Type(particle_list_type), Pointer, Optional:: particles
!!       Type(distribution_1d_type), Pointer, Optional:: local_particles
!!       Type(mol_kind_list_type), Pointer, Optional:: molecule_kinds
!!       Type(mol_struct_list_type), Pointer, Optional:: local_molecules
!!       Type(mol_new_list_type), Pointer, Optional:: molecules_new
!!       Type(mol_kind_new_list_type), Pointer, Optional::&
!!         molecule_kinds_new
!!       Type(distribution_1d_type), Pointer, Optional::&
!!         local_molecules_new
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Type(particle_node_type), Dimension(:), Pointer, Optional:: pnode
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine fragment_get
!!
!!   FUNCTION
!!     returns information about various attributes of the given fragment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - fragment: the fragment you want info about
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!     For the description of the other arguments see cp_fragment_type
!!     attributes.
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE fragment_get(fragment, id_nr, ref_count, atomic_kinds, particles,&
     local_particles, molecule_kinds, local_molecules,molecules_new,&
     molecule_kinds_new, local_molecules_new, para_env,pnode, error)
  TYPE (cp_fragment_type), POINTER :: fragment
  INTEGER, INTENT(out), OPTIONAL :: id_nr, ref_count
  TYPE (atomic_kind_list_type), POINTER, OPTIONAL :: atomic_kinds
  TYPE (particle_list_type), POINTER, OPTIONAL :: particles
  TYPE (distribution_1d_type), POINTER, OPTIONAL :: local_particles
  TYPE (mol_kind_list_type), POINTER, OPTIONAL :: molecule_kinds
  TYPE (mol_struct_list_type), POINTER, OPTIONAL :: local_molecules
  TYPE (mol_new_list_type), POINTER, OPTIONAL :: molecules_new
  TYPE (mol_kind_new_list_type), POINTER, OPTIONAL :: molecule_kinds_new
  TYPE (distribution_1d_type), POINTER, OPTIONAL :: local_molecules_new
  TYPE (cp_para_env_type), POINTER, OPTIONAL :: para_env
  TYPE (particle_node_type), DIMENSION(:), POINTER, OPTIONAL :: pnode
  TYPE (cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fragment_get',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(fragment),cp_failure_level,routineP,error,failure)
  CPPrecondition(fragment%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(id_nr)) id_nr = fragment%id_nr
     IF (PRESENT(ref_count)) ref_count = fragment%ref_count
     IF (PRESENT(atomic_kinds)) atomic_kinds => fragment%atomic_kinds
     IF (PRESENT(particles)) particles => fragment%particles
     IF (PRESENT(local_particles)) local_particles => fragment%local_particles
     IF (PRESENT(molecule_kinds)) molecule_kinds => fragment%molecule_kinds
     IF (PRESENT(local_molecules)) local_molecules => fragment%local_molecules
     IF (PRESENT(molecules_new)) molecules_new => fragment%molecules_new
     IF (PRESENT(molecule_kinds_new)) molecule_kinds_new => fragment%molecule_kinds_new
     IF (PRESENT(local_molecules_new)) local_molecules_new => fragment%local_molecules_new
     IF (PRESENT(para_env)) para_env => fragment%para_env
     IF (PRESENT(pnode)) pnode => fragment%pnode
  END IF
END SUBROUTINE fragment_get
!***************************************************************************

! to do: read from input

END MODULE cp_fragment_types
