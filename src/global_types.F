!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****h* cp2k/global_types [1.0] *
!!
!!   NAME
!!     global_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     - print keys, basis_set_file name and potential_file_name added to the
!!       global type (27.02.2001, MK)
!!     - JGH (28.11.2001) : added pp_library_path to type
!!     - Merged with MODULE print_keys (17.01.2002, MK)
!!
!!   SOURCE
!******************************************************************************

MODULE global_types

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
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: default_string_length
  USE print_key_types,                 ONLY: DEBUG,&
                                             HIGH,&
                                             LOW,&
                                             MEDIUM,&
                                             SILENT,&
                                             WANNIER_CRAZY,&
                                             WANNIER_DIRECT,&
                                             WANNIER_JACOBI,&
                                             WANNIER_NONE,&
                                             init_print_keys,&
                                             print_key_type

  IMPLICIT NONE

  PRIVATE

! *** Public parameters ***

  PUBLIC :: DEBUG,HIGH,LOW,MEDIUM,SILENT
  PUBLIC :: BLACS_GRID_SQUARE, BLACS_GRID_ROW, BLACS_GRID_COL
  PUBLIC :: WANNIER_JACOBI, WANNIER_CRAZY, WANNIER_DIRECT, WANNIER_NONE

! *** Public data types

  PUBLIC :: global_environment_type,&
            print_key_type

! *** Public subroutines ***

  PUBLIC :: init_print_keys
  PUBLIC :: cp_para_env_from_globenv

  INTEGER, PARAMETER :: BLACS_GRID_SQUARE=1, &
                        BLACS_GRID_ROW=2, &
                        BLACS_GRID_COL=3

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="global_types"

!!*****
!!****s* global_types/global_environment_type [1.0] *
!!
!!   NAME
!!     global_environment_type
!!
!!   FUNCTION
!!     
!!   NOTES
!!     -
!!
!!   SOURCE
!!
  TYPE global_environment_type
     TYPE(cp_para_env_type), POINTER         :: para_env
     TYPE(print_key_type)                    :: PRINT
     CHARACTER(LEN=default_string_length)    :: basis_set_file_name,&
                                                coord_file_name,&
                                                default_fft_library,&
                                                error_file_name,&
                                                input_file_name,&
                                                output_file_name,&
                                                potential_file_name,&
                                                mm_potential_file_name,&
                                                pp_library_path,&
                                                program_name,&
                                                project_name,&
                                                restart_file_name
     CHARACTER(LEN=default_string_length)    :: run_type
     INTEGER :: group        !! group handle for all processors
     INTEGER :: idum         !! random number seed
     INTEGER :: mepos        !! position number for local processor
     INTEGER :: num_pe       !! total number of processors
     INTEGER :: print_level  !! how much information is written
     INTEGER :: scr          !! output unit
     INTEGER :: source       !! processor id for global output node
     LOGICAL :: ionode       !! this processor performs global output
     INTEGER :: blacs_grid_layout !! will store the user preference for the blacs grid
  END TYPE global_environment_type

! *****************************************************************************

CONTAINS

!!*****

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
!!     should try to generate a para env??
!!
!!   INPUTS
!!     - global_env: the global environment that contains the para_env
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!*** **********************************************************************
FUNCTION cp_para_env_from_globenv(global_env, error) RESULT(res)
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: global_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cp_para_env_type), POINTER          :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_para_env_from_globenv', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(global_env%para_env),cp_failure_level,routineP,error,failure)
  res => global_env%para_env
END FUNCTION cp_para_env_from_globenv
!***************************************************************************

END MODULE global_types
!!*****
