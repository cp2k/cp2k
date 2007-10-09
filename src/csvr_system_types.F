!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!***** cp2k/csvr_system_types [1.0] *
!!
!!   NAME
!!     csvr_system_types
!!
!!   FUNCTION
!!     Type for the canonical sampling through velocity rescaling     
!!
!!   AUTHOR
!!     Teodoro Laino - 09.2007 University of Zurich [tlaino]
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE csvr_system_types
  USE extended_system_types,           ONLY: create_map_info_type,&
                                             map_info_type,&
                                             release_map_info_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: csvr_system_type,&
            csvr_info_type,&
            csvr_init,&
            csvr_dealloc

  TYPE csvr_info_type
     REAL (KIND=dp)                          :: nkt
  END TYPE csvr_info_type

  TYPE csvr_system_type
     INTEGER                                 :: region, num_csvr
     REAL(KIND=dp)                           :: tau_csvr
     REAL(KIND=dp)                           :: thermostat_energy
     TYPE(csvr_info_type), POINTER           :: nvt(:)
     TYPE(map_info_type), POINTER            :: map_info
  END TYPE csvr_system_type

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'csvr_system_types'

CONTAINS

! *****************************************************************************
!!****f* extended_system_types/csvr_init [1.0] *
!!
!!   NAME
!!    csvr_init
!!
!!   FUNCTION
!!     Initialize type for Canonical Sampling through Velocity Rescaling (CSVR)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] 10.2007- University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!
!!*** **********************************************************************
  SUBROUTINE csvr_init(csvr, section, error)
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_init', &
      routineP = moduleN//':'//routineN

    NULLIFY(csvr%nvt)
    NULLIFY(csvr%map_info)
    csvr%num_csvr=0
    CALL section_vals_val_get(section,"TIMECON",r_val=csvr%tau_csvr,error=error)
    CALL create_map_info_type(csvr%map_info, error)

  END SUBROUTINE csvr_init

! *****************************************************************************
!!****f* extended_system_types/csvr_dealloc [1.0] *
!!
!!   NAME
!!     csvr_dealloc
!!
!!   FUNCTION
!!     Deallocate type for CSVR thermostat
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] 10.2007- University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!
!!*** **********************************************************************
  SUBROUTINE csvr_dealloc ( csvr, error )
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(csvr)) THEN
       IF (ASSOCIATED (csvr%nvt)) THEN
          DEALLOCATE (csvr%nvt, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       ENDIF
       CALL release_map_info_type(csvr%map_info, error)
       DEALLOCATE (csvr, STAT=stat)
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    ENDIF

  END SUBROUTINE csvr_dealloc

END MODULE csvr_system_types

