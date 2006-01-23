!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****h* cp2k/global_types [1.0] *
!!
!!   NAME
!!     global_types
!!
!!   FUNCTION
!!     represent the global information of a run: the input file, 
!!     parallel environment (and partially output).
!!     Use sparingly, try not to send it too deep in your structures.
!!
!!   AUTHOR
!!     JGH,MK,fawzi
!!
!!   MODIFICATION HISTORY
!!     - print keys, basis_set_file name and potential_file_name added to the
!!       global type (27.02.2001, MK)
!!     - JGH (28.11.2001) : added pp_library_path to type
!!     - Merged with MODULE print_keys (17.01.2002, MK)
!!     - input_file, reference counting, create (08.2004, fawzi)
!!
!!   SOURCE
!******************************************************************************

MODULE global_types

  USE cp_files,                        ONLY: close_file
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: m_walltime

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Public parameters ***

! *** Print level parameters ***
  INTEGER, PARAMETER, PUBLIC :: SILENT = 0,&
                        LOW = 1,&
                        MEDIUM = 2,&
                        HIGH = 3,&
                        DEBUG = 4

  PUBLIC :: BLACS_GRID_SQUARE, BLACS_GRID_ROW, BLACS_GRID_COL

! *** Public data types

  PUBLIC :: global_environment_type

! *** Public subroutines ***

  PUBLIC :: cp_para_env_from_globenv, globenv_create, globenv_retain, &
       globenv_release

  INTEGER, PARAMETER :: BLACS_GRID_SQUARE=1, &
                        BLACS_GRID_ROW=2, &
                        BLACS_GRID_COL=3

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="global_types"
  INTEGER, SAVE, PRIVATE :: last_globenv_id=0

!!*** ************************************************************************

!!****s* global_types/global_environment_type [1.0] *
!!
!!   NAME
!!     global_environment_type
!!
!!   FUNCTION
!!     contains the initally parsed file and the initial parallel environment
!!     
!!   NOTES
!!     This is not but really should be passed as pointer and use reference
!!     counting. Use it accordingly wherever possible.
!!
!!   ATTRIBUTES
!!     - id_nr: identification number (unique)
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - para_env: parallel environment
!!     - input_file: the input file (what you need is probably 
!!       input_file%root_vals, or a subsection, try to use that not
!!       the whole input_file
!!
!!     Personally I think that all the other attributes should go away
!!     (and maybe add  a logger)[fawzi]
!!
!!   SOURCE
!!
  TYPE global_environment_type
     INTEGER :: id_nr, ref_count
     TYPE(cp_para_env_type), POINTER         :: para_env
     CHARACTER(LEN=default_path_length)      :: input_file_name
     CHARACTER(LEN=default_string_length)    :: default_fft_library,&
                                                project_name,&
                                                coord_file_name
     INTEGER :: group        !! group handle for all processors
     INTEGER :: idum         !! random number seed
     INTEGER :: mepos        !! position number for local processor
     INTEGER :: num_pe       !! total number of processors
     INTEGER :: print_level  !! how much information is written
     INTEGER :: prog_name_id !! index to define the type of program
     INTEGER :: run_type_id  !! index to define the run_tupe
     INTEGER :: scr          !! output unit
     INTEGER :: source       !! processor id for global output node
     LOGICAL :: ionode       !! this processor performs global output
     INTEGER :: blacs_grid_layout !! will store the user preference for the blacs grid
     REAL(KIND=dp) :: cp2k_start_time, cp2k_target_time
     TYPE(section_vals_type), POINTER :: input_file
     LOGICAL :: close_output_unit
  END TYPE global_environment_type
!!***
! *****************************************************************************

CONTAINS

!!****f* global_types/globenv_create [1.0] *
!!
!!   NAME
!!     globenv_create
!!
!!   FUNCTION
!!     creates a globenv
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - globenv: the globenv to create
!!     - para_env: the parallel environment of the globenv
!!     - input_file_name: the input file name (remove?)
!!     - output_unit: the output unit
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE globenv_create(globenv, para_env,  output_unit, &
       input_file_name, close_output_unit, error)
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(in)                      :: output_unit
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: input_file_name
    LOGICAL, INTENT(in), OPTIONAL            :: close_output_unit
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'globenv_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
    ALLOCATE(globenv,stat=stat)
    IF (.NOT. failure) THEN
       last_globenv_id=last_globenv_id+1
       globenv%id_nr=last_globenv_id
       globenv%ref_count=1
       globenv%para_env=> para_env
       globenv%close_output_unit=.FALSE.
       IF(PRESENT(close_output_unit))globenv%close_output_unit=close_output_unit
       CALL cp_para_env_retain(para_env,error=error)
       globenv%input_file_name=""
       IF (PRESENT(input_file_name)) globenv%input_file_name=input_file_name
       globenv%run_type_id=0
       globenv%default_fft_library="FFTSG"
       globenv%prog_name_id=0
       globenv%project_name=""
       globenv%group=para_env%group
       globenv%idum=0                 !! random number seed
       globenv%coord_file_name=""
       globenv%mepos=para_env%mepos
       globenv%num_pe=para_env%num_pe
       globenv%print_level=1
       globenv%scr=output_unit
       globenv%source=para_env%source
       globenv%ionode=(para_env%source==para_env%mepos)
       globenv%blacs_grid_layout=BLACS_GRID_SQUARE
       globenv%cp2k_start_time=m_walltime()
       NULLIFY(globenv%input_file)
    END IF
  END SUBROUTINE globenv_create
!***************************************************************************


!!****f* global_types/globenv_retain [1.0] *
!!
!!   NAME
!!     globenv_retain
!!
!!   FUNCTION
!!     retains the global environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - globenv: the global environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE globenv_retain(globenv, error)
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'globenv_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(globenv%ref_count>0,cp_failure_level,routineP,error,failure)
     globenv%ref_count=globenv%ref_count+1
  END IF
END SUBROUTINE globenv_retain
!***************************************************************************

!!****f* global_types/globenv_release [1.0] *
!!
!!   NAME
!!     globenv_release
!!
!!   FUNCTION
!!     releases the global environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - globenv: the global environment to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE globenv_release(globenv, error)
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'globenv_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(globenv)) THEN
     CPPrecondition(globenv%ref_count>0,cp_failure_level,routineP,error,failure)
     globenv%ref_count=globenv%ref_count-1
     IF (globenv%ref_count==0) THEN
        CALL section_vals_release(globenv%input_file,error=error)
        CALL cp_para_env_release(globenv%para_env,error=error)
        IF (globenv%close_output_unit) THEN
           IF (globenv%ionode) CALL close_file(globenv%scr)
        END IF
        DEALLOCATE(globenv,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(globenv)
END SUBROUTINE globenv_release
!***************************************************************************

!!****f* global_types/cp_para_env_from_globenv [1.0] *
!!
!!   NAME
!!     cp_para_env_from_globenv
!!
!!   FUNCTION
!!     returns the parallel environment stored in the global environment
!!     gives a warning if not initialized
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - globenv: the global environment that contains the para_env
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
FUNCTION cp_para_env_from_globenv(globenv, error) RESULT(res)
  TYPE(global_environment_type),POINTER      :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cp_para_env_type), POINTER          :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_para_env_from_globenv', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(globenv%para_env),cp_failure_level,routineP,error,failure)
  res => globenv%para_env
END FUNCTION cp_para_env_from_globenv
!***************************************************************************

END MODULE global_types
