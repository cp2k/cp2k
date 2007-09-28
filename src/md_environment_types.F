!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****s* cp2k/md_environment_types [1.0] *
!!
!!   NAME
!!     md_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM SEPT-12-02
!!
!!   MODIFICATION HISTORY
!!     give the md_env its own para_env Joost VandeVondele 07.2003
!!     Teodoro Laino - 09.2007 - University of Zurich - generalizing thermostats
!!                               and barostats
!!
!!   SOURCE
!******************************************************************************

MODULE md_environment_types
  USE barostat_types,                  ONLY: barostat_type,&
                                             release_barostat_type,&
                                             retain_barostat_type
  USE cell_types,                      ONLY: cell_type
  USE cp_files,                        ONLY: close_file
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE force_env_types,                 ONLY: force_env_release,&
                                             force_env_retain,&
                                             force_env_type
  USE free_energy_types,               ONLY: fe_env_release,&
                                             free_energy_type
  USE input_constants,                 ONLY: do_thermo_nose
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE thermostat_types,                ONLY: release_thermostat_type,&
                                             retain_thermostat_type,&
                                             thermostat_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  TYPE md_environment_type
     ! para_env is the parallel enviroment of the MD,  i.e. the systems 
     ! that are dealt with by the integrator e.g in the PIMD this could 
     ! be parent of every bead.
     PRIVATE
     INTEGER :: id_nr, ref_count
     LOGICAL :: init, first_time
     INTEGER, POINTER :: itimes, md_approach_nr, accepted_md
     REAL (KIND=dp), POINTER            :: used_time, t
     REAL (KIND=dp), POINTER            :: constant
     TYPE (cp_para_env_type), POINTER   :: para_env
     TYPE (cell_type), POINTER          :: cell
     TYPE (force_env_type), POINTER     :: force_env
     TYPE (thermostat_type), POINTER    :: thermostat
     TYPE (barostat_type), POINTER      :: barostat
     TYPE (reftraj_type), POINTER       :: reftraj_info
     TYPE (free_energy_type) , POINTER  :: fe_env
     TYPE (simulation_parameters_type), POINTER :: simpar
  END TYPE md_environment_type

  TYPE reftraj_type
     INTEGER :: unit
  END TYPE reftraj_type

! *** Public subroutines and data types ***
  PUBLIC :: md_environment_type,&
            reftraj_type,&
            set_md_env,&
            get_md_env,&
            md_env_create,&
            md_env_release,&
            md_env_retain

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_environment_types'
  INTEGER, SAVE, PRIVATE :: last_md_env_id=0

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****f* md_environment_types/md_env_create [1.0] *
!!
!!   NAME
!!     md_env_create
!!
!!   FUNCTION
!!     Creates MD environment
!!     Purpose: Initialise the integrator environment.
!!     retain the para_env for this environment (should be used for parallel
!!     communications)
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - md_env: the force environment to retain
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE md_env_create ( md_env, para_env, error )
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'md_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    ALLOCATE ( md_env, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    last_md_env_id=last_md_env_id+1
    md_env%id_nr=last_md_env_id
    md_env%ref_count=1

    NULLIFY ( md_env % itimes )
    NULLIFY ( md_env % constant )
    NULLIFY ( md_env % t )
    NULLIFY ( md_env % md_approach_nr )
    NULLIFY ( md_env % accepted_md )
    NULLIFY ( md_env % cell )
    NULLIFY ( md_env % simpar )
    NULLIFY ( md_env % thermostat )
    NULLIFY ( md_env % barostat )
    NULLIFY ( md_env % force_env )
    NULLIFY ( md_env % fe_env )
    md_env % para_env => para_env
    CALL cp_para_env_retain ( md_env % para_env ,error=error)
    ALLOCATE( md_env % itimes , stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    ALLOCATE( md_env % md_approach_nr , stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    ALLOCATE( md_env % accepted_md , stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    ALLOCATE( md_env % constant , stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    ALLOCATE( md_env % used_time , stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    ALLOCATE( md_env % t , stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    md_env % itimes    = -1
    md_env % md_approach_nr = -1
    md_env % accepted_md = 0
    md_env % constant  = 0.0_dp
    md_env % used_time = 0.0_dp
    md_env % t         = 0.0_dp
    ALLOCATE ( md_env % reftraj_info )
    md_env % reftraj_info % unit =-1
    md_env%init = .TRUE.
    md_env%first_time = .TRUE.

  END SUBROUTINE md_env_create
!***************************************************************************

!!****f* md_environment_types/md_env_retain [1.0] *
!!
!!   NAME
!!     md_env_retain
!!
!!   FUNCTION
!!     retains the given md env
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - md_env: the force environment to retain
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE md_env_retain(md_env, error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'md_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(md_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(md_env%ref_count>0,cp_failure_level,routineP,error)
     md_env%ref_count=md_env%ref_count+1
  END IF
END SUBROUTINE md_env_retain
!***************************************************************************
!!****f* md_environment_types/md_env_release [1.0] *
!!
!!   NAME
!!     md_env_release
!!
!!   FUNCTION
!!     releases the given md env
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - md_env: the md environment to release
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE md_env_release(md_env, error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'md_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

  IF (ASSOCIATED(md_env)) THEN
     CPPreconditionNoFail(md_env%ref_count>0,cp_failure_level,routineP,error)
     md_env%ref_count=md_env%ref_count-1
     IF (md_env%ref_count==0) THEN
        md_env%ref_count=1
        CALL fe_env_release( md_env%fe_env, error=error)
        CALL cp_para_env_release ( md_env % para_env, error = error )
        DEALLOCATE( md_env % itimes , stat=istat )
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        DEALLOCATE( md_env % md_approach_nr , stat=istat )
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        DEALLOCATE( md_env % accepted_md , stat=istat )
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        DEALLOCATE( md_env % constant , stat=istat )
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        DEALLOCATE( md_env % used_time , stat=istat )
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        DEALLOCATE( md_env % t , stat=istat )
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        IF ( md_env%reftraj_info%unit /= -1 ) THEN 
           CALL close_file(unit_number=md_env%reftraj_info%unit)
        END IF
        DEALLOCATE( md_env % reftraj_info , stat=istat )
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        NULLIFY ( md_env % cell )
        NULLIFY ( md_env % simpar )
        CALL release_barostat_type(md_env%barostat, error)
        CALL release_thermostat_type(md_env%thermostat, error)
        CALL force_env_release(md_env%force_env, error=error)
        md_env%ref_count=0
        DEALLOCATE(md_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(md_env)
END SUBROUTINE md_env_release

! *****************************************************************************
!!****f* md_environment_types/get_md_env [1.0] *
!!
!!   NAME
!!     get_md_env
!!
!!   FUNCTION
!!     get components of MD environment type
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - md_env: the force environment to retain
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE get_md_env ( md_env, itimes, constant, used_time, cell,  &
                          simpar, nhc_part, nhc_coef, nhc_baro, nhc_shell, &
                          npt, force_env, para_env, reftraj_info, t,&
                          md_approach_nr, accepted_md, init, first_time,&
                          fe_env, thermostat, barostat, error)

    TYPE(md_environment_type), POINTER       :: md_env
    INTEGER, OPTIONAL, POINTER               :: itimes
    REAL(KIND=dp), OPTIONAL, POINTER         :: constant, used_time
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(simulation_parameters_type), &
      OPTIONAL, POINTER                      :: simpar
    TYPE(lnhc_parameters_type), OPTIONAL, &
      POINTER                                :: nhc_part, nhc_coef, nhc_baro, &
                                                nhc_shell
    TYPE(npt_info_type), OPTIONAL, POINTER   :: npt( :, : )
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(reftraj_type), OPTIONAL, POINTER    :: reftraj_info
    REAL(KIND=dp), OPTIONAL, POINTER         :: t
    INTEGER, OPTIONAL, POINTER               :: md_approach_nr, accepted_md
    LOGICAL, OPTIONAL                        :: init, first_time
    TYPE(free_energy_type), OPTIONAL, &
      POINTER                                :: fe_env
    TYPE(thermostat_type), OPTIONAL, POINTER :: thermostat
    TYPE(barostat_type), OPTIONAL, POINTER   :: barostat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_md_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check, failure

    failure = .FALSE.
    check = ASSOCIATED(md_env)
    CPPrecondition(check,cp_failure_level,routineP,error,failure)
    IF (PRESENT ( itimes ) ) itimes => md_env % itimes
    IF (PRESENT ( fe_env ) ) fe_env => md_env % fe_env
    IF (PRESENT ( md_approach_nr ) ) md_approach_nr => md_env % md_approach_nr
    IF (PRESENT ( accepted_md ) ) accepted_md => md_env % accepted_md
    IF (PRESENT ( constant ) ) constant => md_env % constant
    IF (PRESENT ( used_time ) ) used_time => md_env % used_time
    IF (PRESENT ( t ) ) t => md_env % t
    IF (PRESENT ( cell ) ) cell => md_env % cell
    IF (PRESENT ( simpar ) ) simpar => md_env % simpar
    IF (PRESENT ( thermostat ) ) thermostat => md_env % thermostat
    IF (PRESENT ( barostat ) ) barostat => md_env % barostat
    IF (PRESENT(nhc_part).OR.PRESENT(nhc_coef).OR.PRESENT(nhc_baro).OR.&
        PRESENT(nhc_shell)) THEN
       IF (ASSOCIATED(md_env%thermostat)) THEN
          IF (md_env%thermostat%type_of_thermostat==do_thermo_nose) THEN
             IF (PRESENT ( nhc_part ) )  nhc_part => md_env%thermostat%nhc%nhc_part
             IF (PRESENT ( nhc_coef ) )  nhc_coef => md_env%thermostat%nhc%nhc_coef
             IF (PRESENT ( nhc_baro ) )  nhc_baro => md_env%thermostat%nhc%nhc_baro
             IF (PRESENT ( nhc_shell ) ) nhc_shell => md_env%thermostat%nhc%nhc_shell
          END IF
       END IF
    END IF
    IF (PRESENT(npt)) THEN
       IF (ASSOCIATED(md_env%barostat)) THEN
          npt => md_env%barostat%npt
       END IF
    END IF
    IF (PRESENT ( force_env) ) force_env => md_env % force_env
    IF (PRESENT ( para_env ) ) para_env => md_env % para_env
    IF (PRESENT ( reftraj_info ) ) reftraj_info => md_env % reftraj_info
    IF (PRESENT (init)) init=md_env%init 
    IF (PRESENT (first_time)) first_time=md_env%first_time 

  END SUBROUTINE get_md_env

! *****************************************************************************
!!****f* md_environment_types/set_md_env [1.0] *
!!
!!   NAME
!!     set_md_env
!!
!!   FUNCTION
!!     Set the integrator environment to the correct program.
!!
!!   NOTES
!!     
!!
!!   ARGUMENTS
!!     - md_env: the force environment to retain
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE set_md_env ( md_env, itimes, constant, cell, simpar, fe_env, force_env,&
                          para_env, accepted_md, md_approach_nr , init, first_time,&
                          thermostat, barostat, error)

    TYPE(md_environment_type), POINTER       :: md_env
    INTEGER, OPTIONAL, POINTER               :: itimes
    REAL(KIND=dp), OPTIONAL, POINTER         :: constant
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(simulation_parameters_type), &
      OPTIONAL, POINTER                      :: simpar
    TYPE(free_energy_type), OPTIONAL, &
      POINTER                                :: fe_env
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    INTEGER, OPTIONAL, POINTER               :: accepted_md, md_approach_nr
    LOGICAL, OPTIONAL                        :: init, first_time
    TYPE(thermostat_type), OPTIONAL, POINTER :: thermostat
    TYPE(barostat_type), OPTIONAL, POINTER   :: barostat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_md_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check, failure

    failure = .FALSE.
    check = ASSOCIATED(md_env)
    CPPrecondition(check,cp_failure_level,routineP,error,failure)
    IF ( PRESENT (init)) md_env%init = init
    IF ( PRESENT (first_time)) md_env%first_time = first_time
    IF ( PRESENT ( cell ) ) md_env % cell => cell
    IF ( PRESENT ( barostat ) ) THEN 
       CALL release_barostat_type(md_env % barostat, error=error)
       CALL retain_barostat_type(barostat, error=error)
       md_env % barostat => barostat
    END IF
    IF ( PRESENT ( thermostat ) ) THEN
       CALL release_thermostat_type(md_env % thermostat, error=error)
       CALL retain_thermostat_type(thermostat, error=error)
       md_env % thermostat => thermostat
    END IF
    IF ( PRESENT ( simpar ) ) md_env % simpar => simpar
    IF ( PRESENT ( itimes ) ) md_env % itimes => itimes
    IF ( PRESENT ( fe_env ) ) md_env % fe_env => fe_env
    IF ( PRESENT ( accepted_md ) ) md_env % accepted_md => accepted_md
    IF ( PRESENT ( md_approach_nr ) ) md_env % md_approach_nr => md_approach_nr
    IF ( PRESENT ( constant ) ) md_env % constant => constant
    IF ( PRESENT ( force_env ) ) THEN
       IF (ASSOCIATED(force_env)) THEN
          CALL force_env_retain (force_env, error=error)
       END IF
       IF (ASSOCIATED(md_env%force_env)) THEN
          CALL force_env_release (md_env%force_env, error=error)
       END IF
       md_env%force_env => force_env
    END IF

  END SUBROUTINE set_md_env

! *****************************************************************************

END MODULE md_environment_types
!******************************************************************************
