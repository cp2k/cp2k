!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/particle_list_types [1.0] *
!!
!!   NAME
!!     particle_list_types
!!
!!   FUNCTION
!!     represent a simple array based list of the given type
!!
!!   NOTES
!!     This is a template
!!
!!     **** DO NOT MODIFY THE .F FILES ****
!!     modify the .template and .instantition
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE particle_list_types
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
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: deallocate_particle_set,&
                                             particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='particle_list_types'
  INTEGER, PRIVATE, SAVE :: last_particle_list_id=0

  PUBLIC :: particle_list_type, particle_list_p_type
  PUBLIC :: particle_list_create, particle_list_retain,&
       particle_list_release
       
!***
!****************************************************************************

!!****s* particle_list_types/particle_list_type [1.0] *
!!
!!   NAME
!!     particle_list_type
!!
!!   FUNCTION
!!     represent a list of objects
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: identification number of this list
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - n_el: the number of elements in the list
!!     - owns_list: if the els are owned by this structure, and
!!       should be deallocated by it
!!     - list: the array of object, might be oversized,
!!       only the fist n_el have some meaning
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE particle_list_type
     INTEGER :: id_nr, ref_count, n_els
     LOGICAL :: owns_els
     TYPE(particle_type), DIMENSION(:), POINTER :: els
  END TYPE particle_list_type
!!***
!****************************************************************************

!!****s* particle_list_types/particle_list_p_type [1.0] *
!!
!!   NAME
!!     particle_list_p_type
!!
!!   FUNCTION
!!     represents a pointer to a list
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - list: the pointer to the list
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE particle_list_p_type
     TYPE(particle_list_type), POINTER :: list
  END TYPE particle_list_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* particle_list_types/particle_list_create [1.0] *
!!
!!   NAME
!!     particle_list_create
!!
!!   SYNOPSIS
!!     Subroutine particle_list_create(list, els_ptr, owns_els, n_els,&
!!         error)
!!       Type(particle_list_type), Pointer, Optional:: list
!!       Type(particle_type), Dimension(:), Pointer, Optional:: els_ptr
!!       Logical, Intent (IN), Optional:: owns_els
!!       Integer, Intent (IN), Optional:: n_els
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine particle_list_create
!!
!!   FUNCTION
!!     creates a list
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - list: the list to allocate and initialize
!!     - els_ptr: the elements to store in the list (the array is only,
!!       referenced, not copied!)
!!     - owns_els: if the list takes the ownership of els_ptr and
!!       will deallocate it (defaults to true)
!!     - n_els: number of elements in the list (at least one els_ptr or
!!       n_els should be given)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE particle_list_create(list, els_ptr, &
     owns_els, n_els, error)
  TYPE(particle_list_type), POINTER, OPTIONAL :: list
  TYPE(particle_type), DIMENSION(:), POINTER, OPTIONAL :: els_ptr
  LOGICAL, INTENT(in), OPTIONAL :: owns_els
  INTEGER, INTENT(in), OPTIONAL :: n_els
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='particle_list_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.

  CPPrecondition(PRESENT(els_ptr).OR.PRESENT(n_els),cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     ALLOCATE(list, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     last_particle_list_id=last_particle_list_id+1
     list%id_nr=last_particle_list_id
     list%ref_count=1
     list%owns_els=.TRUE.
     list%n_els=0
     IF (PRESENT(owns_els)) list%owns_els=owns_els
     NULLIFY(list%els)
     IF (PRESENT(els_ptr)) THEN
        list%els => els_ptr
        IF (ASSOCIATED(els_ptr)) THEN
           list%n_els=SIZE(els_ptr)
        END IF
     END IF
     IF (PRESENT(n_els)) list%n_els=n_els
     IF (.NOT.ASSOCIATED(list%els)) THEN
        ALLOCATE(list%els(list%n_els),stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        CPPreconditionNoFail(list%owns_els,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE particle_list_create
!***************************************************************************

!!****f* particle_list_types/particle_list_retain [1.0] *
!!
!!   NAME
!!     particle_list_retain
!!
!!   SYNOPSIS
!!     Subroutine particle_list_retain(list, error)
!!       Type(particle_list_type), Pointer:: list
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine particle_list_retain
!!
!!   FUNCTION
!!     retains a list (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - list: the list to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE particle_list_retain(list, error)
  TYPE(particle_list_type), POINTER :: list
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='particle_list_retain',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(list),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(list%ref_count>0,cp_failure_level,routineP,error,failure)
     list%ref_count=list%ref_count+1
  END IF
END SUBROUTINE particle_list_retain
!***************************************************************************

!!****f* particle_list_types/particle_list_release [1.0] *
!!
!!   NAME
!!     particle_list_release
!!
!!   SYNOPSIS
!!     Subroutine particle_list_release(list, error)
!!       Type(particle_list_type), Pointer:: list
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine particle_list_release
!!
!!   FUNCTION
!!     releases a  list (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - list: the list to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE particle_list_release(list, error)
  TYPE(particle_list_type), POINTER :: list
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='particle_list_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.
  
  IF (ASSOCIATED(list)) THEN
     CPPrecondition(list%ref_count>0,cp_failure_level,routineP,error,failure)
     list%ref_count=list%ref_count-1
     IF (list%ref_count==0) THEN
        IF (list%owns_els) THEN
           IF (ASSOCIATED(list%els)) THEN
              CALL deallocate_particle_set(list%els)
           END IF
        END IF
        NULLIFY(list%els)
        DEALLOCATE(list,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(list)
END SUBROUTINE particle_list_release
!***************************************************************************

! template def put here so that line numbers in template and derived 
! files are almost the same (multi-line use change it a bit)
! [template(el_typename,el_type,USE,deallocate_els_code)]
! ARGS:
!  USE = "use particle_types, only: particle_type, deallocate_particle_set"
!  deallocate_els_code = "call deallocate_particle_set(list%els)"
!  el_type = "type(particle_type)"
!  el_typename = "particle"


END MODULE particle_list_types
