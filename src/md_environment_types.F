!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      give the md_env its own para_env Joost VandeVondele 07.2003
!>      Teodoro Laino - 09.2007 - University of Zurich - generalizing thermostats
!>                                and barostats
!> \author CJM SEPT-12-02
! *****************************************************************************
MODULE md_environment_types
  USE barostat_types,                  ONLY: barostat_type,&
                                             release_barostat_type,&
                                             retain_barostat_type
  USE cell_types,                      ONLY: cell_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE extended_system_types,           ONLY: npt_info_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_release,&
                                             force_env_retain,&
                                             force_env_type
  USE free_energy_types,               ONLY: fe_env_release,&
                                             free_energy_type
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE reftraj_types,                   ONLY: reftraj_type,&
                                             release_reftraj,&
                                             retain_reftraj
  USE thermostat_types,                ONLY: release_thermostats,&
                                             retain_thermostats,&
                                             thermostat_type,&
                                             thermostats_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
  TYPE md_ener_type
     INTEGER       :: ref_count
     INTEGER       :: nfree, nfree_qm, nfree_shell
     REAL(KIND=dp) :: constant
     REAL(KIND=dp) :: delta_cons, delta_epot
     REAL(KIND=dp) :: epot
     REAL(KIND=dp) :: ekin, ekin_qm
     REAL(KIND=dp) :: temp_part, temp_qm
     REAL(KIND=dp) :: temp_baro
     REAL(KIND=dp) :: ekin_coefs
     REAL(KIND=dp) :: temp_coefs
     REAL(KIND=dp) :: ekin_shell, temp_shell
     REAL(KIND=dp) :: thermostat_part_kin, thermostat_part_pot
     REAL(KIND=dp) :: thermostat_baro_kin, thermostat_baro_pot
     REAL(KIND=dp) :: thermostat_coef_kin, thermostat_coef_pot
     REAL(KIND=dp) :: thermostat_shell_kin, thermostat_shell_pot
     REAL(KIND=dp) :: baro_kin, baro_pot
     REAL(KIND=dp) :: vcom(3), total_mass
     REAL(KIND=dp), DIMENSION(:), POINTER :: ekin_kind
     REAL(KIND=dp), DIMENSION(:), POINTER :: temp_kind
     INTEGER, DIMENSION(:), POINTER :: nfree_kind
     REAL(KIND=dp), DIMENSION(:), POINTER :: ekin_shell_kind
     REAL(KIND=dp), DIMENSION(:), POINTER :: temp_shell_kind
     INTEGER, DIMENSION(:), POINTER :: nfree_shell_kind
  END TYPE md_ener_type

! *****************************************************************************
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
     TYPE (md_ener_type), POINTER       :: md_ener
     TYPE (thermostats_type), POINTER   :: thermostats
     TYPE (barostat_type), POINTER      :: barostat
     TYPE (reftraj_type), POINTER       :: reftraj
     TYPE (free_energy_type) , POINTER  :: fe_env
     TYPE (simulation_parameters_type), POINTER :: simpar
  END TYPE md_environment_type

! *** Public subroutines and data types ***
  PUBLIC :: create_md_ener,&
            md_ener_type,&
            md_environment_type,&
            set_md_env,&
            get_md_env,&
            md_env_create,&
            md_env_release,&
            md_env_retain,&
            release_md_ener,&
            retain_md_ener,&
            zero_md_ener

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_environment_types'
  INTEGER, SAVE, PRIVATE :: last_md_env_id=0

CONTAINS

! *****************************************************************************
!> \brief Creates MD environment
!>      Purpose: Initialise the integrator environment.
!>      retain the para_env for this environment (should be used for parallel
!>      communications)
!> \param md_env the force environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
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
    NULLIFY ( md_env % thermostats )
    NULLIFY ( md_env % barostat )
    NULLIFY ( md_env % force_env )
    NULLIFY ( md_env % fe_env )
    NULLIFY ( md_env % md_ener )
    NULLIFY ( md_env % reftraj )
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
    md_env%init = .TRUE.
    md_env%first_time = .TRUE.

  END SUBROUTINE md_env_create

! *****************************************************************************
!> \brief retains the given md env
!> \param md_env the force environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      see doc/ReferenceCounting.html
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
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

! *****************************************************************************
!> \brief releases the given md env
!> \param md_env the md environment to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      see doc/ReferenceCounting.html
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
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
        NULLIFY ( md_env % cell )
        NULLIFY ( md_env % simpar )
        CALL release_barostat_type(md_env%barostat, error)
        CALL release_thermostats(md_env%thermostats, error)
        CALL release_reftraj(md_env%reftraj, error)
        CALL release_md_ener(md_env%md_ener, error)
        CALL force_env_release(md_env%force_env, error=error)
        md_env%ref_count=0
        DEALLOCATE(md_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(md_env)
END SUBROUTINE md_env_release

! *****************************************************************************
!> \brief get components of MD environment type
!> \param md_env the force environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE get_md_env ( md_env, itimes, constant, used_time, cell,  &
                          simpar, npt, force_env, para_env, reftraj, t,&
                          md_approach_nr, accepted_md, init, first_time,&
                          fe_env, thermostats, barostat, thermostat_coeff, &
                          thermostat_part, thermostat_shell, thermostat_baro,&
                          md_ener, error)

    TYPE(md_environment_type), POINTER       :: md_env
    INTEGER, OPTIONAL, POINTER               :: itimes
    REAL(KIND=dp), OPTIONAL, POINTER         :: constant, used_time
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(simulation_parameters_type), &
      OPTIONAL, POINTER                      :: simpar
    TYPE(npt_info_type), OPTIONAL, POINTER   :: npt( :, : )
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(reftraj_type), OPTIONAL, POINTER    :: reftraj
    REAL(KIND=dp), OPTIONAL, POINTER         :: t
    INTEGER, OPTIONAL, POINTER               :: md_approach_nr, accepted_md
    LOGICAL, OPTIONAL                        :: init, first_time
    TYPE(free_energy_type), OPTIONAL, &
      POINTER                                :: fe_env
    TYPE(thermostats_type), OPTIONAL, &
      POINTER                                :: thermostats
    TYPE(barostat_type), OPTIONAL, POINTER   :: barostat
    TYPE(thermostat_type), OPTIONAL, POINTER :: thermostat_coeff, &
                                                thermostat_part, &
                                                thermostat_shell, &
                                                thermostat_baro
    TYPE(md_ener_type), OPTIONAL, POINTER    :: md_ener
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
    IF (PRESENT ( thermostats ) ) thermostats => md_env % thermostats
    IF (PRESENT ( barostat ) ) barostat => md_env % barostat
    IF  (PRESENT(thermostat_part).OR.PRESENT(thermostat_coeff).OR.&
         PRESENT(thermostat_baro).OR.PRESENT(thermostat_shell)) THEN
       IF (ASSOCIATED(md_env%thermostats)) THEN
          IF (PRESENT ( thermostat_part ) ) THEN
             thermostat_part => md_env%thermostats%thermostat_part
          END IF
          IF (PRESENT ( thermostat_coeff ) ) THEN
             thermostat_coeff => md_env%thermostats%thermostat_coef
          END IF
          IF (PRESENT ( thermostat_shell ) ) THEN
             thermostat_shell => md_env%thermostats%thermostat_shell
          END IF
          IF (PRESENT ( thermostat_baro ) ) THEN
             thermostat_baro  => md_env%thermostats%thermostat_baro
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
    IF (PRESENT ( reftraj ) ) reftraj => md_env % reftraj
    IF (PRESENT ( md_ener ) ) md_ener => md_env % md_ener
    IF (PRESENT (init)) init=md_env%init 
    IF (PRESENT (first_time)) first_time=md_env%first_time 

  END SUBROUTINE get_md_env

! *****************************************************************************
!> \brief Set the integrator environment to the correct program.
!> \param md_env the force environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE set_md_env ( md_env, itimes, constant, cell, simpar, fe_env, force_env,&
                          para_env, accepted_md, md_approach_nr , init, first_time,&
                          thermostats, barostat, reftraj, md_ener, error)

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
    TYPE(thermostats_type), OPTIONAL, &
      POINTER                                :: thermostats
    TYPE(barostat_type), OPTIONAL, POINTER   :: barostat
    TYPE(reftraj_type), OPTIONAL, POINTER    :: reftraj
    TYPE(md_ener_type), OPTIONAL, POINTER    :: md_ener
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
    IF ( PRESENT ( thermostats ) ) THEN
       CALL release_thermostats(md_env%thermostats, error=error)
       CALL retain_thermostats(thermostats, error=error)
       md_env%thermostats => thermostats
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
    IF ( PRESENT ( reftraj ) ) THEN
       CALL release_reftraj(md_env%reftraj, error=error)
       CALL retain_reftraj(reftraj, error=error)
       md_env%reftraj => reftraj
    END IF
    IF ( PRESENT ( md_ener ) ) THEN
       CALL release_md_ener(md_env%md_ener, error=error)
       CALL retain_md_ener(md_ener, error=error)
       md_env%md_ener => md_ener
    END IF

  END SUBROUTINE set_md_env

! *****************************************************************************
!> \brief retains the given md_ener structure
!> \param md_ener :
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [MI]
!> \author MI
! *****************************************************************************
SUBROUTINE create_md_ener(md_ener, error)
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_md_ener', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(.NOT.ASSOCIATED(md_ener),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(md_ener,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    md_ener%ref_count = 1

    NULLIFY(md_ener%ekin_kind)
    NULLIFY(md_ener%temp_kind)
    NULLIFY(md_ener%nfree_kind)
    NULLIFY(md_ener%ekin_shell_kind)
    NULLIFY(md_ener%temp_shell_kind)
    NULLIFY(md_ener%nfree_shell_kind)

  END IF
END SUBROUTINE create_md_ener
  
! *****************************************************************************
!> \brief retains the given md_ener structure
!> \param md_ener :
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [MI]
!> \author MI
! *****************************************************************************
SUBROUTINE retain_md_ener(md_ener, error)
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_md_ener', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(md_ener),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(md_ener%ref_count>0,cp_failure_level,routineP,error)
     md_ener%ref_count=md_ener%ref_count+1
  END IF
END SUBROUTINE retain_md_ener

! *****************************************************************************
!> \brief releases the given md_ener structure
!> \param md_ener :
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [MI]
!> \author MI
! *****************************************************************************
SUBROUTINE release_md_ener(md_ener, error)
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_md_ener', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(md_ener)) THEN
     CPPreconditionNoFail(md_ener%ref_count>0,cp_failure_level,routineP,error)
     md_ener%ref_count=md_ener%ref_count-1
     IF (md_ener%ref_count==0) THEN
        IF(ASSOCIATED(md_ener%temp_kind))THEN
          DEALLOCATE(md_ener%temp_kind,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF(ASSOCIATED(md_ener%ekin_kind))THEN
          DEALLOCATE(md_ener%ekin_kind,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF(ASSOCIATED(md_ener%nfree_kind))THEN
          DEALLOCATE(md_ener%nfree_kind,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF(ASSOCIATED(md_ener%temp_shell_kind))THEN
          DEALLOCATE(md_ener%temp_shell_kind,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF(ASSOCIATED(md_ener%ekin_shell_kind))THEN
          DEALLOCATE(md_ener%ekin_shell_kind,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF(ASSOCIATED(md_ener%nfree_shell_kind))THEN
          DEALLOCATE(md_ener%nfree_shell_kind,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
  
        DEALLOCATE(md_ener,stat=istat)
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE release_md_ener

! *****************************************************************************
!> \brief initialize to zero energies and temperatures  
!> \param md_ener :
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [MI]
!> \author MI
! *****************************************************************************
SUBROUTINE zero_md_ener(md_ener, tkind, tshell, error)
    TYPE(md_ener_type), POINTER              :: md_ener
    LOGICAL, INTENT(IN)                      :: tkind, tshell
    TYPE(cp_error_type), INTENT(inout)       :: error

    md_ener%ekin        = 0.0_dp
    md_ener%temp_part   = 0.0_dp
    md_ener%temp_baro   = 0.0_dp
    md_ener%ekin_coefs  = 0.0_dp
    md_ener%temp_coefs  = 0.0_dp
    md_ener%ekin_qm     = 0.0_dp
    md_ener%temp_qm     = 0.0_dp
    md_ener%ekin_shell  = 0.0_dp
    md_ener%temp_shell  = 0.0_dp    
    md_ener%constant    = 0.0_dp
    md_ener%delta_cons  = 0.0_dp
    md_ener%delta_epot  = 0.0_dp
    md_ener%thermostat_part_kin  = 0.0_dp
    md_ener%thermostat_part_pot  = 0.0_dp
    md_ener%thermostat_coef_kin  = 0.0_dp
    md_ener%thermostat_coef_pot  = 0.0_dp
    md_ener%thermostat_baro_kin  = 0.0_dp
    md_ener%thermostat_baro_pot  = 0.0_dp
    md_ener%thermostat_shell_kin  = 0.0_dp
    md_ener%thermostat_shell_pot  = 0.0_dp
    md_ener%baro_kin = 0.0_dp
    md_ener%baro_pot = 0.0_dp
    IF(tkind) THEN
       md_ener%temp_kind  = 0.0_dp
       md_ener%ekin_kind  = 0.0_dp

       IF(tshell) THEN
          md_ener%temp_shell_kind  = 0.0_dp
          md_ener%ekin_shell_kind  = 0.0_dp 
       END IF
    END IF        
    md_ener%vcom(:) = 0.0_dp
    md_ener%total_mass = 0.0_dp
    
END SUBROUTINE zero_md_ener

END MODULE md_environment_types
