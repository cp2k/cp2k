!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief represent a simple array based list of the given type
!> \note
!>     ____              _ _     __  __           _ _  __         _____ _     _       _____ _ _      _
!>    |  _ \  ___  _ __ ( ) |_  |  \/  | ___   __| (_)/ _|_   _  |_   _| |__ (_)___  |  ___(_) | ___| |
!>    | | | |/ _ \| '_ \|/| __| | |\/| |/ _ \ / _` | | |_| | | |   | | | '_ \| / __| | |_  | | |/ _ \ |
!>    | |_| | (_) | | | | | |_  | |  | | (_) | (_| | |  _| |_| |   | | | | | | \__ \ |  _| | | |  __/_|
!>    |____/ \___/|_| |_|  \__| |_|  |_|\___/ \__,_|_|_|  \__, |   |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!>                                                        |___/
!>      ____ _                  ___                              _ _       _       _
!>     / ___| | ___  ___  ___  |_ _|_ __ ___  _ __ ___   ___  __| (_) __ _| |_ ___| |_   _
!>    | |   | |/ _ \/ __|/ _ \  | || '_ ` _ \| '_ ` _ \ / _ \/ _` | |/ _` | __/ _ \ | | | |
!>    | |___| | (_) \__ \  __/  | || | | | | | | | | | |  __/ (_| | | (_| | ||  __/ | |_| |
!>     \____|_|\___/|___/\___| |___|_| |_| |_|_| |_| |_|\___|\__,_|_|\__,_|\__\___|_|\__, |
!>                                                                                   |___/
!>     _____ _     _       _____ _ _      _
!>    |_   _| |__ (_)___  |  ___(_) | ___| |
!>      | | | '_ \| / __| | |_  | | |/ _ \ |
!>      | | | | | | \__ \ |  _| | | |  __/_|
!>      |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!>
!>      This is a template
!>
!>      **** DO NOT MODIFY THE .F FILES ****
!>      modify the .template and .instantition instead
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE mol_new_list_types
  USE molecule_types_new,              ONLY: deallocate_molecule_set,&
                                             molecule_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mol_new_list_types'
  INTEGER, PRIVATE, SAVE :: last_mol_new_list_id=0

  !API
  PUBLIC :: mol_new_list_type, mol_new_list_p_type,&
            mol_new_list_create, mol_new_list_retain,&
            mol_new_list_release

!***

! *****************************************************************************
!> \brief represent a list of objects
!> \param id_nr identification number of this list
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param n_el the number of elements in the list
!> \param owns_list if the els are owned by this structure, and
!>        should be deallocated by it
!> \param list the array of object, might be oversized,
!>        only the fist n_el have some meaning
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE mol_new_list_type
     INTEGER :: id_nr, ref_count, n_els
     LOGICAL :: owns_els
     TYPE(molecule_type), DIMENSION(:), POINTER :: els
  END TYPE mol_new_list_type

! *****************************************************************************
!> \brief represents a pointer to a list
!> \param list the pointer to the list
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE mol_new_list_p_type
     TYPE(mol_new_list_type), POINTER :: list
  END TYPE mol_new_list_p_type

CONTAINS

! *****************************************************************************
!> \brief creates a list
!> \param list the list to allocate and initialize
!> \param els_ptr the elements to store in the list (the array is only,
!>        referenced, not copied!)
!> \param owns_els if the list takes the ownership of els_ptr and
!>        will deallocate it (defaults to true)
!> \param n_els number of elements in the list (at least one els_ptr or
!>        n_els should be given)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE mol_new_list_create(list, els_ptr, &
     owns_els, n_els, error)
    TYPE(mol_new_list_type), OPTIONAL, &
      POINTER                                :: list
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: els_ptr
    LOGICAL, INTENT(in), OPTIONAL            :: owns_els
    INTEGER, INTENT(in), OPTIONAL            :: n_els
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief retains a list (see doc/ReferenceCounting.html)
!> \param list the list to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE mol_new_list_retain(list, error)
    TYPE(mol_new_list_type), POINTER         :: list
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief releases a  list (see doc/ReferenceCounting.html)
!> \param list the list to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE mol_new_list_release(list, error)
    TYPE(mol_new_list_type), POINTER         :: list
    TYPE(cp_error_type), INTENT(inout)       :: error

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
              CALL deallocate_molecule_set(list%els,error=error)
           END IF
        END IF
        NULLIFY(list%els)
        DEALLOCATE(list,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(list)
END SUBROUTINE mol_new_list_release

! template def put here so that line numbers in template and derived
! files are almost the same (multi-line use change it a bit)
! [template(el_typename,el_type,USE,deallocate_els_code)]
! ARGS:
!  USE = "use molecule_types_new, only: molecule_type, deallocate_molecule_set"
!  deallocate_els_code = "call deallocate_molecule_set(list%els,error=error)"
!  el_type = "type(molecule_type)"
!  el_typename = "mol_new"


END MODULE mol_new_list_types
