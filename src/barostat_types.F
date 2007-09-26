!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!***** cp2k/barostat_types [1.0] *
!!
!!   NAME
!!     barostat_types
!!
!!   FUNCTION
!!     Barostat structure: module containing barostat available for MD
!!
!!   AUTHOR
!!     teo [tlaino] - University of Zurich - 09.2007
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE barostat_types
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE extended_system_types,           ONLY: npt_info_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type
  USE global_types,                    ONLY: global_environment_type
  USE init_extended_system_types,      ONLY: initialize_npt
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: barostat_type,&
            create_barostat_type,&
            release_barostat_type,&
            retain_barostat_type

  TYPE barostat_type
     INTEGER :: id_nr, ref_count, type_of_barostat
     TYPE ( npt_info_type ), POINTER :: npt (:,:)
  END TYPE barostat_type

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'barostat_types'
  INTEGER, PRIVATE, SAVE               :: last_barostat_id_nr=0

CONTAINS

!!****f* barostat_types/create_barostat_type [1.0] *
!!
!!   NAME
!!     create_barostat_type
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE create_barostat_type( barostat, motion_section, force_env_p,&
       simpar, para_env, globenv, error )
    TYPE(barostat_type), POINTER             :: barostat
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_barostat_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: check, explicit, failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: barostat_section

    failure =.FALSE.
    check = .NOT.ASSOCIATED(barostat)
    CPPrecondition(check,cp_failure_level,routineP,error,failure)
    barostat_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT",error=error)
    CALL section_vals_get(barostat_section, explicit=explicit, error=error)
    IF (explicit) THEN
       ALLOCATE(barostat, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       last_barostat_id_nr = last_barostat_id_nr + 1
       barostat%id_nr      = last_barostat_id_nr
       barostat%ref_count  =  1
       NULLIFY(barostat%npt)
       CPPrecondition(SIZE(force_env_p)==1,cp_failure_level,routineP,error,failure)
       CALL force_env_get( force_env_p(1)%force_env, cell=cell, error=error )
       ! Initialize or possibly restart Barostat
       CALL initialize_npt (simpar, para_env, globenv, barostat%npt,&
            cell, work_section=barostat_section, error=error)
    END IF
    
  END SUBROUTINE create_barostat_type

!!****f* barostat_types/retain_barostat_type [1.0] *
!!
!!   NAME
!!     retain_barostat_type
!!
!!   FUNCTION
!!     retains the given barostat
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     -
!!
!!   AUTHOR
!!    Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [tlaino]
!! 
!!*** **********************************************************************
  SUBROUTINE retain_barostat_type(barostat,error)
    TYPE(barostat_type), POINTER             :: barostat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_barostat_type', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    
    IF (.NOT. failure) THEN
       IF (ASSOCIATED(barostat)) THEN
          CPPrecondition(barostat%ref_count>0,cp_failure_level,routineP,error,failure)
          barostat%ref_count=barostat%ref_count+1
       END IF
    END IF
  END SUBROUTINE retain_barostat_type

!!****f* barostat_types/release_barostat_type [1.0] *
!!
!!   NAME
!!     release_barostat_type
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE release_barostat_type(barostat, error)
    TYPE(barostat_type), POINTER             :: barostat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_barostat_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: check, failure

    failure =.FALSE.
    IF (ASSOCIATED(barostat)) THEN
       check = barostat%ref_count>0
       CPPrecondition(check,cp_failure_level,routineP,error,failure)
       barostat%ref_count=barostat%ref_count-1
       IF (barostat%ref_count<1) THEN
          IF ( ASSOCIATED ( barostat%npt ) ) THEN
             DEALLOCATE (barostat%npt , stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          DEALLOCATE(barostat, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF

  END SUBROUTINE release_barostat_type
!!*****
!******************************************************************************

END MODULE barostat_types
!******************************************************************************
