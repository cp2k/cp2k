!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/input_file_types [1.0] *
!!
!!   NAME
!!     input_file_types
!!
!!   FUNCTION
!!     objects that reperesnt the structure and content of an input file
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE input_file_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE input_section_types,             ONLY: section_create,&
                                             section_describe,&
                                             section_release,&
                                             section_type,&
                                             section_vals_create,&
                                             section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_write
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='input_file_types'
  INTEGER, SAVE, PRIVATE :: last_structure_id=0, last_input_id=0

  PUBLIC :: input_structure_type, input_file_type
  PUBLIC :: input_structure_create, input_structure_retain, &
       input_structure_release, input_structure_describe
  PUBLIC :: input_file_create, input_file_retain, input_file_release,&
       input_file_write
!!***
!****************************************************************************

!!****s* input_file_types/input_structure_type *
!!
!!   NAME
!!     input_structure_type
!!
!!   FUNCTION
!!     represent the structure of an input file
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: identification number (unique)
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - root_section: pseudo section that represents the whole file
!!
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
  !***************************************************************************
  TYPE input_structure_type
     INTEGER :: id_nr, ref_count
     TYPE(section_type), POINTER :: root_section
  END TYPE input_structure_type
!!***
  !****************************************************************************

!!****s* input_file_types/input_file_type [1.0] *
!!
!!   NAME
!!     input_file_type
!!
!!   FUNCTION
!!     represents a parsed input file
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
  !***************************************************************************
  TYPE input_file_type
     INTEGER :: id_nr, ref_count
     TYPE(input_structure_type), POINTER :: structure
     CHARACTER(len=default_path_length) :: path
     TYPE(section_vals_type), POINTER :: root_vals
  END TYPE input_file_type
!!***
  !****************************************************************************

CONTAINS

!!****f* input_file_types/input_structure_create *
!!
!!   NAME
!!     input_structure_create
!!
!!   FUNCTION
!!     creates an objects that represents the structure of the input
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - structure: the input structure to create
!!     - description: a description of the structure
!!     - n_sections: a hint about how many sections will be added to this
!!       file structure
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE input_structure_create(structure,description,n_sections,error)
    TYPE(input_structure_type), POINTER      :: structure
    CHARACTER(len=*), INTENT(in)             :: description
    INTEGER, INTENT(in), OPTIONAL            :: n_sections
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'input_structure_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_n_sections, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  my_n_sections=10
  IF (PRESENT(n_sections)) my_n_sections=n_sections

  CPPrecondition(.NOT.ASSOCIATED(structure),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(structure,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     
     structure%ref_count=1
     last_structure_id=last_structure_id+1
     structure%id_nr=last_structure_id
     NULLIFY(structure%root_section)
     CALL section_create(structure%root_section,name="_ROOT_",&
          description=description,n_keywords=0,n_subsections=n_sections,&
          repeats=.FALSE.,error=error)
  END IF
END SUBROUTINE input_structure_create
!***************************************************************************

!!****f* input_file_types/input_structure_retain *
!!
!!   NAME
!!     input_structure_retain
!!
!!   FUNCTION
!!     retains the given input structure (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - structure: the input structure to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE input_structure_retain(structure,error)
    TYPE(input_structure_type), POINTER      :: structure
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'input_structure_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(structure),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(structure%ref_count>0,cp_failure_level,routineP,error)
     structure%ref_count=structure%ref_count+1
  END IF
END SUBROUTINE input_structure_retain
!***************************************************************************

!!****f* input_file_types/input_structure_release *
!!
!!   NAME
!!     input_structure_release
!!
!!   FUNCTION
!!     releases the given input structure (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - structure: the input structure to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE input_structure_release(structure,error)
    TYPE(input_structure_type), POINTER      :: structure
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'input_structure_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(structure)) THEN
     CPPreconditionNoFail(structure%ref_count>0,cp_failure_level,routineP,error)
     structure%ref_count=structure%ref_count-1
     IF (structure%ref_count==0) THEN
        CALL section_release(structure%root_section,error=error)
        DEALLOCATE(structure,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(structure)
END SUBROUTINE input_structure_release
!***************************************************************************

!!****f* input_file_types/input_structure_describe *
!!
!!   NAME
!!     input_structure_describe
!!
!!   FUNCTION
!!     describes the structure of the input
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - structure: the structure to describe
!!     - unit_nr: the unit where to write
!!     - level: the level of output (0: only section names, 1: + keywords
!!       names,... see keyword_describe)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE input_structure_describe(structure, unit_nr, level, error)
    TYPE(input_structure_type), POINTER      :: structure
    INTEGER, INTENT(in)                      :: unit_nr, level
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'input_structure_describe', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(structure),cp_failure_level,routineP,error,failure)
  CPPrecondition(structure%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     WRITE (unit=unit_nr, fmt=&
          "('***************************************************************')")
     CALL section_describe(structure%root_section,unit_nr,level,&
             hide_root=.TRUE.,error=error)
     WRITE (unit=unit_nr, fmt=&
          "('***************************************************************')")
  END IF
END SUBROUTINE input_structure_describe
!***************************************************************************

!!****f* input_file_types/input_file_create *
!!
!!   NAME
!!     input_file_create
!!
!!   FUNCTION
!!     creates a place to store the values of an input file
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - input_file: the object to be created
!!     - structure: the structure of the input file
!!     - path: the path of the input file (empty if there is no real file)
!!     - root_vals: the values stored in this object (if not given or not
!!       associated an empty container is created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE input_file_create(input_file,structure,path,root_vals,error)
    TYPE(input_file_type), POINTER           :: input_file
    TYPE(input_structure_type), POINTER      :: structure
    CHARACTER(len=*), INTENT(in)             :: path
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: root_vals
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'input_file_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(.NOT.ASSOCIATED(input_file),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(input_file,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     last_input_id=last_input_id+1
     input_file%id_nr=last_input_id
     input_file%ref_count=1
     CPPrecondition(LEN(input_file%path)>=LEN_TRIM(path),cp_failure_level,routineP,error,failure)
     input_file%path=path
     input_file%structure => structure
     CALL input_structure_retain(structure,error=error)
     NULLIFY(input_file%root_vals)
     IF (PRESENT(root_vals)) THEN
        input_file%root_vals => root_vals
     END IF
     IF (ASSOCIATED(input_file%root_vals)) THEN
        CALL section_vals_retain(input_file%root_vals,error=error)        
     ELSE
        CALL section_vals_create(input_file%root_vals,&
             structure%root_section,error=error)
     END IF
  END IF
END SUBROUTINE input_file_create
!***************************************************************************

!!****f* input_file_types/input_file_retain *
!!
!!   NAME
!!     input_file_retain
!!
!!   FUNCTION
!!     retains the given input file
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - input_file: the object to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE input_file_retain(input_file,error)
    TYPE(input_file_type), POINTER           :: input_file
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'input_file_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(input_file),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(input_file%ref_count>0,cp_failure_level,routineP,error)
     input_file%ref_count=input_file%ref_count+1
  END IF
END SUBROUTINE input_file_retain
!***************************************************************************

!!****f* input_file_types/input_file_release *
!!
!!   NAME
!!     input_file_release
!!
!!   FUNCTION
!!     releases the given input file
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - input_file: the object to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE input_file_release(input_file,error)
    TYPE(input_file_type), POINTER           :: input_file
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'input_file_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(input_file)) THEN
     CPPreconditionNoFail(input_file%ref_count>0,cp_failure_level,routineP,error)
     input_file%ref_count=input_file%ref_count+1
     IF (input_file%ref_count==0) THEN
        CALL input_structure_release(input_file%structure, error=error)
        CALL section_vals_release(input_file%root_vals,error=error)
        DEALLOCATE(input_file,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(input_file)
END SUBROUTINE input_file_release
!***************************************************************************

!!****f* input_file_types/input_file_write [1.0] *
!!
!!   NAME
!!     input_file_write
!!
!!   FUNCTION
!!     writes out the actual input file in a way that can be parsed back in,
!!     and making all defaults explicit-
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - input_file: the input file to write out
!!     - unit_nr: unit where to write the input file
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE input_file_write(input_file, unit_nr, error)
    TYPE(input_file_type), POINTER           :: input_file
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'input_file_write', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(input_file),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(input_file%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     CALL section_vals_write(input_file%root_vals,unit_nr,hide_root=.TRUE.,&
          error=error)
  END IF
END SUBROUTINE input_file_write
!***************************************************************************

END MODULE input_file_types
