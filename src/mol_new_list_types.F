!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/mol_new_list_types [1.0] *
!!
!!   NAME
!!     mol_new_list_types
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
MODULE mol_new_list_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE molecule_types_new,              ONLY: deallocate_molecule_set,&
                                             molecule_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='mol_new_list_types'
  INTEGER, PRIVATE, SAVE :: last_mol_new_list_id=0

  PUBLIC :: mol_new_list_type, mol_new_list_p_type
  PUBLIC :: mol_new_list_create, mol_new_list_retain,&
       mol_new_list_release
       
!***
!****************************************************************************

!!****s* mol_new_list_types/mol_new_list_type [1.0] *
!!
!!   NAME
!!     mol_new_list_type
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
  TYPE mol_new_list_type
     INTEGER :: id_nr, ref_count, n_els
     LOGICAL :: owns_els
     TYPE(molecule_type), DIMENSION(:), POINTER :: els
  END TYPE mol_new_list_type
!!***
!****************************************************************************

!!****s* mol_new_list_types/mol_new_list_p_type [1.0] *
!!
!!   NAME
!!     mol_new_list_p_type
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
  TYPE mol_new_list_p_type
     TYPE(mol_new_list_type), POINTER :: list
  END TYPE mol_new_list_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* mol_new_list_types/mol_new_list_create [1.0] *
!!
!!   NAME
!!     mol_new_list_create
!!
!!   SYNOPSIS
!!     Subroutine mol_new_list_create(list, els_ptr, owns_els, n_els,&
!!         error)
!!       Type(mol_new_list_type), Optional, Pointer:: list
!!       Type(molecule_type), Dimension(:), Optional, Pointer:: els_ptr
!!       Logical, Intent (IN), Optional:: owns_els
!!       Integer, Intent (IN), Optional:: n_els
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine mol_new_list_create
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
SUBROUTINE mol_new_list_create(list, els_ptr, &
     owns_els, n_els, error)
    TYPE(mol_new_list_type), OPTIONAL, &
      POINTER                                :: list
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: els_ptr
    LOGICAL, INTENT(in), OPTIONAL            :: owns_els
    INTEGER, INTENT(in), OPTIONAL            :: n_els
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mol_new_list_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(PRESENT(els_ptr).OR.PRESENT(n_els),cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     ALLOCATE(list, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     last_mol_new_list_id=last_mol_new_list_id+1
     list%id_nr=last_mol_new_list_id
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
END SUBROUTINE mol_new_list_create
!***************************************************************************

!!****f* mol_new_list_types/mol_new_list_retain [1.0] *
!!
!!   NAME
!!     mol_new_list_retain
!!
!!   SYNOPSIS
!!     Subroutine mol_new_list_retain(list, error)
!!       Type(mol_new_list_type), Pointer:: list
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine mol_new_list_retain
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
SUBROUTINE mol_new_list_retain(list, error)
    TYPE(mol_new_list_type), POINTER         :: list
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mol_new_list_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(list),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(list%ref_count>0,cp_failure_level,routineP,error,failure)
     list%ref_count=list%ref_count+1
  END IF
END SUBROUTINE mol_new_list_retain
!***************************************************************************

!!****f* mol_new_list_types/mol_new_list_release [1.0] *
!!
!!   NAME
!!     mol_new_list_release
!!
!!   SYNOPSIS
!!     Subroutine mol_new_list_release(list, error)
!!       Type(mol_new_list_type), Pointer:: list
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine mol_new_list_release
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
SUBROUTINE mol_new_list_release(list, error)
    TYPE(mol_new_list_type), POINTER         :: list
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mol_new_list_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(list)) THEN
     CPPrecondition(list%ref_count>0,cp_failure_level,routineP,error,failure)
     list%ref_count=list%ref_count-1
     IF (list%ref_count==0) THEN
        IF (list%owns_els) THEN
           IF (ASSOCIATED(list%els)) THEN
              CALL deallocate_molecule_set(list%els)
           END IF
        END IF
        NULLIFY(list%els)
        DEALLOCATE(list,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(list)
END SUBROUTINE mol_new_list_release
!***************************************************************************

! template def put here so that line numbers in template and derived 
! files are almost the same (multi-line use change it a bit)
! [template(el_typename,el_type,USE,deallocate_els_code)]
! ARGS:
!  USE = "use molecule_types_new, only: molecule_type, deallocate_molecule_set"
!  deallocate_els_code = "call deallocate_molecule_set(list%els)"
!  el_type = "type(molecule_type)"
!  el_typename = "mol_new"


END MODULE mol_new_list_types
