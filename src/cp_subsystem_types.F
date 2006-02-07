!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/cp_subsystem_types [1.0] *
!!
!!   NAME
!!     cp_subsystem_types
!!
!!   FUNCTION
!!     types that represent a subsys, i.e. a part of the system
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
MODULE cp_subsystem_types
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_release,&
                                             atomic_kind_list_retain,&
                                             atomic_kind_list_type
  USE colvar_types,                    ONLY: colvar_dealloc_ref,&
                                             colvar_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_retain,&
                                             distribution_1d_type
  USE kinds,                           ONLY: dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_release,&
                                             mol_kind_new_list_retain,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_release,&
                                             mol_new_list_retain,&
                                             mol_new_list_type
  USE particle_list_types,             ONLY: particle_list_release,&
                                             particle_list_retain,&
                                             particle_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_subsystem_types'
  INTEGER, PRIVATE, SAVE :: last_fragment_id=0

  PUBLIC :: cp_subsystem_type, cp_subsystem_p_type,cp_subsystem_pp_type
  PUBLIC :: cp_subsys_create, cp_subsys_retain, cp_subsys_release,&
       cp_subsys_get, cp_subsys_set

!***
!****************************************************************************

!!****s* force_env_types/cp_subsystem_type [1.0] *
!!
!!   NAME
!!     cp_subsystem_type
!!
!!   FUNCTION
!!     represents a part of a system: atoms, molecules, their pos,vel,...
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - atomic_kinds: list with all the kinds in the actual subsys
!!     - particles: list with the particles of the actual subsys
!!     - local_particles: the particles that are local to the actual processor
!!     - molecule_kinds: list with the molecule kinds
!!     - local_molecules: the molecule structures of the actual subsys
!!       that are local to this processor
!!     - para_env: the parallel environment of the actual subsys
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_subsystem_type
!     PRIVATE !FM removed for convenience, but to as if it where there
     INTEGER :: ref_count, id_nr
     TYPE (atomic_kind_list_type), POINTER :: atomic_kinds
     TYPE (particle_list_type), POINTER :: particles
     TYPE (distribution_1d_type), POINTER :: local_particles
     TYPE (cp_para_env_type), POINTER :: para_env
!new
     TYPE (mol_new_list_type), POINTER :: molecules_new
     TYPE (mol_kind_new_list_type), POINTER :: molecule_kinds_new
     TYPE (distribution_1d_type), POINTER :: local_molecules_new
! definitions of the collective variables
     TYPE (colvar_type), DIMENSION(:), POINTER :: colvar
  END TYPE cp_subsystem_type
!!***
  !****************************************************************************

!!****s* cp_subsystem_types/cp_subsystem_p_type [1.0] *
!!
!!   NAME
!!     cp_subsystem_p_type
!!
!!   FUNCTION
!!     represent a pointer to a subsys, to be able to create arrays
!!     of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - subsys: the pointer to the subsys
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE cp_subsystem_p_type
     TYPE(cp_subsystem_type), POINTER :: subsys
  END TYPE cp_subsystem_p_type

!!****s* cp_subsystem_types/cp_subsystem_pp_type [1.0] *
!!
!!   NAME
!!     cp_subsystem_pp_type
!!
!!   FUNCTION
!!     added to make an array of subsys arrays...used for getting the subsys
!!     of multiple force env
!!
!!   ATTRIBUTES
!!     - subsys: the pointer to the array of subsys
!!
!!   AUTHOR
!!     MJM
!!
  !***************************************************************************
  TYPE cp_subsystem_pp_type
     TYPE(cp_subsystem_p_type), DIMENSION(:),POINTER :: subsys
  END TYPE cp_subsystem_pp_type

!!***
  !****************************************************************************

CONTAINS

!!****f* cp_subsystem_types/cp_subsys_create [1.0] *
!!
!!   NAME
!!     cp_subsys_create
!!
!!   SYNOPSIS
!!     Subroutine cp_subsys_create(subsys, para_env, error)
!!       Type(cp_subsystem_type), Pointer:: subsys
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_subsys_create
!!
!!   FUNCTION
!!     creates a subsys
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - subsys: the subsys to allocate and initialize
!!     - para_env: the parallel environment of this subsys
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
SUBROUTINE cp_subsys_create(subsys, para_env, error)
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  ALLOCATE(subsys, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     last_fragment_id=last_fragment_id+1
     subsys%id_nr=last_fragment_id
     subsys%ref_count=1
     CALL cp_para_env_retain(para_env,error=error)
     subsys%para_env => para_env
     NULLIFY(subsys%atomic_kinds, subsys%particles, &
          subsys%local_particles, subsys%molecules_new,&
          subsys%molecule_kinds_new, subsys%local_molecules_new)
     NULLIFY(subsys%colvar)
  END IF
END SUBROUTINE cp_subsys_create
!***************************************************************************

!!****f* cp_subsystem_types/cp_subsys_retain [1.0] *
!!
!!   NAME
!!     cp_subsys_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_subsys_retain(subsys, error)
!!       Type(cp_subsystem_type), Pointer:: subsys
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_subsys_retain
!!
!!   FUNCTION
!!     retains a subsys (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - subsys: the subsys to retain
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
SUBROUTINE cp_subsys_retain(subsys, error)
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(subsys%ref_count>0,cp_failure_level,routineP,error)
     subsys%ref_count=subsys%ref_count+1
  END IF
END SUBROUTINE cp_subsys_retain
!***************************************************************************

!!****f* cp_subsystem_types/cp_subsys_release [1.0] *
!!
!!   NAME
!!     cp_subsys_release
!!
!!   SYNOPSIS
!!     Subroutine cp_subsys_release(subsys, error)
!!       Type(cp_subsystem_type), Pointer:: subsys
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_subsys_release
!!
!!   FUNCTION
!!     releases a subsys (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - subsys: the subsys to release
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
SUBROUTINE cp_subsys_release(subsys, error)
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(subsys)) THEN
     CPPreconditionNoFail(subsys%ref_count>0,cp_failure_level,routineP,error)
     subsys%ref_count=subsys%ref_count-1
     IF (subsys%ref_count==0) THEN
        CALL atomic_kind_list_release(subsys%atomic_kinds,error=error)
        CALL particle_list_release(subsys%particles, error=error)
        CALL distribution_1d_release(subsys%local_particles, error=error)
        CALL mol_kind_new_list_release(subsys%molecule_kinds_new, error=error)
        CALL mol_new_list_release(subsys%molecules_new, error=error)
        CALL distribution_1d_release(subsys%local_molecules_new,error=error)
        CALL cp_para_env_release(subsys%para_env, error=error)
        IF(ASSOCIATED(subsys%colvar)) THEN
           DO i=1,SIZE(subsys%colvar)
              CALL colvar_dealloc_ref(subsys%colvar(i),error)
           ENDDO
           DEALLOCATE(subsys%colvar,STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
        DEALLOCATE(subsys, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(subsys)

END SUBROUTINE cp_subsys_release
!***************************************************************************

!!****f* cp_subsystem_types/cp_subsys_set [1.0] *
!!
!!   NAME
!!     cp_subsys_set
!!
!!   SYNOPSIS
!!     Subroutine cp_subsys_set(subsys, atomic_kinds, particles,&
!!         local_particles, molecule_kinds, local_molecules, molecules_new,&
!!         molecule_kinds_new, local_molecules_new, pnode, para_env, error)
!!       Type(cp_subsystem_type), Pointer:: subsys
!!       Type(atomic_kind_list_type), Pointer, Optional:: atomic_kinds
!!       Type(particle_list_type), Pointer, Optional:: particles
!!       Type(distribution_1d_type), Pointer, Optional:: local_particles
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Type(mol_new_list_type), Pointer, Optional:: molecules_new
!!       Type(mol_kind_new_list_type), Pointer, Optional::&
!!         molecule_kinds_new
!!       Type(distribution_1d_type), Pointer, Optional::&
!!         local_molecules_new
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_subsys_set
!!
!!   FUNCTION
!!     sets various propreties of the subsys
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - subsys: the subsys you want to modify
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!     For the description of the other arguments see cp_subsystem_type
!!     attributes.
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_subsys_set(subsys, atomic_kinds, particles, local_particles,&
                         molecules_new, molecule_kinds_new, local_molecules_new, &
                         para_env, colvar, error)
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(atomic_kind_list_type), OPTIONAL, &
      POINTER                                :: atomic_kinds
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles
    TYPE(mol_new_list_type), OPTIONAL, &
      POINTER                                :: molecules_new
    TYPE(mol_kind_new_list_type), OPTIONAL, &
      POINTER                                :: molecule_kinds_new
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_molecules_new
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(colvar_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: colvar
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
  CPPrecondition(subsys%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(atomic_kinds)) THEN
        CALL atomic_kind_list_retain(atomic_kinds,error=error)
        CALL atomic_kind_list_release(subsys%atomic_kinds, error=error)
        subsys%atomic_kinds => atomic_kinds
     END IF
     IF (PRESENT(particles)) THEN
        CALL particle_list_retain(particles, error=error)
        CALL particle_list_release(subsys%particles, error=error)
        subsys%particles => particles
     END IF
     IF (PRESENT(local_particles)) THEN
        CALL distribution_1d_retain(local_particles,error=error)
        CALL distribution_1d_release(subsys%local_particles,error=error)
        subsys%local_particles => local_particles
     END IF
     IF (PRESENT(local_molecules_new)) THEN
        CALL distribution_1d_retain(local_molecules_new,error=error)
        CALL distribution_1d_release(subsys%local_molecules_new,error=error)
        subsys%local_molecules_new => local_molecules_new
     END IF
     IF (PRESENT(molecule_kinds_new)) THEN
        CALL mol_kind_new_list_retain(molecule_kinds_new, error=error)
        CALL mol_kind_new_list_release(subsys%molecule_kinds_new, error=error)
        subsys%molecule_kinds_new => molecule_kinds_new
     END IF
     IF (PRESENT(molecules_new)) THEN
        CALL mol_new_list_retain(molecules_new, error=error)
        CALL mol_new_list_release(subsys%molecules_new, error=error)
        subsys%molecules_new => molecules_new
     END IF
     IF (PRESENT(para_env)) THEN
        CALL cp_para_env_retain(para_env, error=error)
        CALL cp_para_env_release(subsys%para_env, error=error)
        subsys%para_env => para_env
     END IF
     IF (PRESENT(colvar)) THEN
        CPPrecondition(.NOT.ASSOCIATED(subsys%colvar),cp_failure_level,routineP,error,failure)
        subsys%colvar=>colvar
     ENDIF
  END IF
END SUBROUTINE cp_subsys_set
!***************************************************************************

!!****f* cp_subsystem_types/cp_subsys_get [1.0] *
!!
!!   NAME
!!     cp_subsys_get
!!
!!   SYNOPSIS
!!     Subroutine cp_subsys_get(subsys, id_nr, ref_count, atomic_kinds,&
!!         particles, local_particles, &
!!         molecules_new, molecule_kinds_new, local_molecules_new, para_env,&
!!         error)
!!       Type(cp_subsystem_type), Pointer:: subsys
!!       Integer, Intent (OUT), Optional:: id_nr, ref_count
!!       Type(atomic_kind_list_type), Pointer, Optional:: atomic_kinds
!!       Type(particle_list_type), Pointer, Optional:: particles
!!       Type(distribution_1d_type), Pointer, Optional:: local_particles
!!       Type(mol_new_list_type), Pointer, Optional:: molecules_new
!!       Type(mol_kind_new_list_type), Pointer, Optional::&
!!         molecule_kinds_new
!!       Type(distribution_1d_type), Pointer, Optional::&
!!         local_molecules_new
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_subsys_get
!!
!!   FUNCTION
!!     returns information about various attributes of the given subsys
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - subsys: the subsys you want info about
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!     For the description of the other arguments see cp_subsystem_type
!!     attributes.
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_subsys_get(subsys, id_nr, ref_count, atomic_kinds, particles,&
     local_particles, molecules_new, molecule_kinds_new, local_molecules_new, &
     para_env, colvar, error)
    TYPE(cp_subsystem_type), POINTER         :: subsys
    INTEGER, INTENT(out), OPTIONAL           :: id_nr, ref_count
    TYPE(atomic_kind_list_type), OPTIONAL, &
      POINTER                                :: atomic_kinds
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles
    TYPE(mol_new_list_type), OPTIONAL, &
      POINTER                                :: molecules_new
    TYPE(mol_kind_new_list_type), OPTIONAL, &
      POINTER                                :: molecule_kinds_new
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_molecules_new
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(colvar_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: colvar
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
  CPPrecondition(subsys%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(id_nr)) id_nr = subsys%id_nr
     IF (PRESENT(ref_count)) ref_count = subsys%ref_count
     IF (PRESENT(atomic_kinds)) atomic_kinds => subsys%atomic_kinds
     IF (PRESENT(particles)) particles => subsys%particles
     IF (PRESENT(local_particles)) local_particles => subsys%local_particles
     IF (PRESENT(molecules_new)) molecules_new => subsys%molecules_new
     IF (PRESENT(molecule_kinds_new)) molecule_kinds_new => subsys%molecule_kinds_new
     IF (PRESENT(local_molecules_new)) local_molecules_new => subsys%local_molecules_new
     IF (PRESENT(para_env)) para_env => subsys%para_env
     IF (PRESENT(colvar)) colvar => subsys%colvar
  END IF
END SUBROUTINE cp_subsys_get
!***************************************************************************

! to do: read from input

END MODULE cp_subsystem_types
