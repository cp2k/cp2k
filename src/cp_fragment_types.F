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
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_logger_type
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_a_l, cp_error_check
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  USE mol_struct_list_types, ONLY: mol_struct_list_type, mol_struct_list_retain,&
       mol_struct_list_release
  USE mol_kind_list_types, ONLY: mol_kind_list_type, mol_kind_list_retain,&
       mol_kind_list_release
  USE particle_list_types, ONLY: particle_list_type, particle_list_retain,&
       particle_list_release
  USE atomic_kind_list_types, ONLY: atomic_kind_list_type, atomic_kind_list_retain,&
       atomic_kind_list_release
  USE cp_para_types, ONLY: cp_para_env_type
  USE cp_para_env, ONLY: cp_para_env_retain, cp_para_env_release
  USE cp_local_lists_types, ONLY: cp_local_lists_type, cp_ll_release, cp_ll_retain
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_fragment_types'
  INTEGER, PRIVATE, SAVE :: last_fragment_id=0

  PUBLIC :: cp_fragment_type, cp_fragment_p_type
  PUBLIC :: fragment_create, fragment_retain, fragment_release

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
  !!     -
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
     INTEGER :: ref_count, id_nr
     TYPE (atomic_kind_list_type), POINTER :: atomic_kinds
     TYPE (particle_list_type), POINTER :: particles
     TYPE (cp_local_lists_type), POINTER :: local_particles
     TYPE (mol_kind_list_type), POINTER :: molecule_kinds
     TYPE (mol_struct_list_type), POINTER :: local_molecules
     TYPE (cp_para_env_type), POINTER :: para_env
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
     TYPE(cp_fragment_type), pointer :: fragment
  END TYPE cp_fragment_p_type
  !!***
  !****************************************************************************

contains

!!****f* cp_fragment_types/fragment_create [1.0] *
!!
!!   NAME
!!     fragment_create
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
  TYPE(cp_fragment_type), pointer :: fragment
  TYPE(cp_para_env_type), pointer :: para_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='fragment_create',&
        routineP=moduleN//':'//routineN
  integer :: stat

  failure=.false.
  
  ALLOCATE(fragment, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     last_fragment_id=last_fragment_id+1
     fragment%id_nr=last_fragment_id
     fragment%ref_count=1
     fragment%para_env => para_env
     CALL cp_para_env_retain(para_env,error=error)
     NULLIFY(fragment%atomic_kinds, fragment%particles, &
          fragment%local_particles, fragment%molecule_kinds,&
          fragment%local_molecules)
  end if
END SUBROUTINE fragment_create
!***************************************************************************

!!****f* cp_fragment_types/fragment_retain [1.0] *
!!
!!   NAME
!!     fragment_retain
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
  TYPE(cp_fragment_type), pointer :: fragment
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='fragment_retain',&
        routineP=moduleN//':'//routineN

  failure=.false.
  
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
  TYPE(cp_fragment_type), pointer :: fragment
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='fragment_release',&
        routineP=moduleN//':'//routineN
  integer :: stat

  failure=.false.
  
  IF (ASSOCIATED(fragment)) THEN
     CPPreconditionNoFail(fragment%ref_count>0,cp_failure_level,routineP,error)
     fragment%ref_count=fragment%ref_count-1
     IF (fragment%ref_count==0) THEN
        CALL atomic_kind_list_release(fragment%atomic_kinds,error=error)
        CALL particle_list_release(fragment%particles, error=error)
        CALL cp_ll_release(fragment%local_particles, error=error)
        CALL mol_kind_list_release(fragment%molecule_kinds, error=error)
        CALL mol_struct_list_release(fragment%local_molecules, error=error)
        CALL cp_para_env_release(fragment%para_env, error=error)
        DEALLOCATE(fragment, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  nullify(fragment)
END SUBROUTINE fragment_release
!***************************************************************************

! to do: read from input, set

END MODULE cp_fragment_types
