!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
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
!!
!!   SOURCE
!******************************************************************************

MODULE md_environment_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE force_env_methods,               ONLY: force_env_release,&
                                             force_env_retain,&
                                             force_env_type
  USE kinds,                           ONLY: dbl
  USE md,                              ONLY: simulation_parameters_type,&
                                             virial_type
  USE simulation_cell,                 ONLY: cell_type
!USE termination,                     ONLY: stop_memory


  IMPLICIT NONE

  PRIVATE

  ! para_env is the parallel enviroment of the MD,
  ! i.e. the systems that are dealt with by the integrator
  ! e.g in the PIMD this could be parent of every bead

  TYPE md_environment_type
     PRIVATE
     INTEGER :: id_nr, ref_count, in_use
     TYPE (cp_para_env_type), POINTER :: para_env 
     INTEGER, POINTER :: itimes
     REAL ( dbl ), POINTER  :: used_time
     REAL ( dbl ), POINTER :: constant
     TYPE ( cell_type ), POINTER :: cell
     TYPE ( force_env_type ), POINTER :: force_env  !? should this also be an array for PIMD
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_part ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_baro ( : )
     TYPE ( npt_info_type ), POINTER :: npt ( :, : )
   !  TYPE ( pimd_environment_type ), POINTER :: pimd_env
     TYPE ( simulation_parameters_type ), POINTER :: simpar
     TYPE ( virial_type ), POINTER :: virial
  END TYPE md_environment_type

! *** Public subroutines and data types ***
  PUBLIC :: md_environment_type, set_md_env,   &
            virial_type, get_md_env, md_env_create,  &
            md_env_release, md_env_retain, zero_virial

! *** Global parameters ***

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='md_environment_types'
  INTEGER, SAVE, PRIVATE :: last_md_env_id=0

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

!  SUBROUTINE get_md_env ( md_env, itimes, constant, cell, simpar, pimd_env, &
!                          virial, &
!                          nhc_part, nhc_coef, nhc_baro, npt, &
!                          force_env, para_env )
  SUBROUTINE get_md_env ( md_env, itimes, constant, used_time, cell,  &
                          simpar, virial, &
                          nhc_part, nhc_coef, nhc_baro, npt, force_env, &
                          para_env )

!   Purpose: Get the integrator environment

!   ***************************************************************************

    TYPE(md_environment_type), POINTER       :: md_env
    INTEGER, OPTIONAL, POINTER               :: itimes
    REAL(dbl), OPTIONAL, POINTER             :: constant, used_time
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(simulation_parameters_type), &
      OPTIONAL, POINTER                      :: simpar
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(lnhc_parameters_type), OPTIONAL, &
      POINTER                                :: nhc_part( : ), nhc_coef( : ), &
                                                nhc_baro( : )
    TYPE(npt_info_type), OPTIONAL, POINTER   :: npt( :, : )
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_md_env", &
      routineP = moduleN//':'//routine_name

!    TYPE ( pimd_environment_type ), POINTER :: pimd_env
!   ---------------------------------------------------------------------------

    IF ( PRESENT ( itimes ) ) itimes => md_env % itimes
    IF ( PRESENT ( constant ) ) constant => md_env % constant
    IF ( PRESENT ( used_time ) ) used_time => md_env % used_time
    IF ( PRESENT ( cell ) ) cell => md_env % cell 
    IF ( PRESENT ( simpar ) ) simpar => md_env % simpar
    IF ( PRESENT ( nhc_part ) ) nhc_part => md_env % nhc_part 
    IF ( PRESENT ( nhc_coef ) ) nhc_coef => md_env % nhc_coef
    IF ( PRESENT ( nhc_baro ) ) nhc_baro => md_env % nhc_baro
!    IF ( PRESENT ( pimd_env ) ) pimd_env => md_env % pimd_env
    IF ( PRESENT ( npt ) ) npt => md_env % npt
    IF ( PRESENT ( virial ) ) virial => md_env % virial
    IF ( PRESENT ( force_env ) ) force_env => md_env % force_env
    IF ( PRESENT ( para_env ) ) para_env => md_env % para_env

  END SUBROUTINE get_md_env

! *****************************************************************************

  SUBROUTINE md_env_create ( md_env, para_env, error )

!   Purpose: Initialise the integrator environment.
!   retain the para_env for this environment (should be used for parallel
!   communications)

!   ***************************************************************************

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine_name = "md_env_create", &
      routineP = moduleN//':'//routine_name

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    ALLOCATE ( md_env, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    last_md_env_id=last_md_env_id+1
    md_env%id_nr=last_md_env_id
    md_env%ref_count=1
    md_env%in_use=0

    NULLIFY ( md_env % itimes )
    NULLIFY ( md_env % constant )
    NULLIFY ( md_env % cell )
    NULLIFY ( md_env % simpar )
    NULLIFY ( md_env % nhc_part )
    NULLIFY ( md_env % nhc_coef )
    NULLIFY ( md_env % nhc_baro )
    NULLIFY ( md_env % npt )
    NULLIFY ( md_env % virial )
    NULLIFY ( md_env % force_env )
!    NULLIFY ( md_env % pimd_env )
    md_env % para_env => para_env 
    CALL cp_para_env_retain ( md_env % para_env )
    ALLOCATE( md_env % itimes , stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    ALLOCATE( md_env % constant , stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    ALLOCATE( md_env % used_time , stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    md_env % itimes    = 0
    md_env % constant  = 0
    md_env % used_time = 0

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
        CALL cp_para_env_release ( md_env % para_env, error = error )
        DEALLOCATE( md_env % itimes , stat=istat )
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        DEALLOCATE( md_env % constant , stat=istat )
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        DEALLOCATE( md_env % used_time , stat=istat )
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        NULLIFY ( md_env % cell )
        NULLIFY ( md_env % simpar )
        NULLIFY ( md_env % nhc_part )
        NULLIFY ( md_env % nhc_coef )
        NULLIFY ( md_env % nhc_baro )
        NULLIFY ( md_env % npt )
        NULLIFY ( md_env % virial )
        CALL force_env_release( md_env % force_env) 
        md_env%ref_count=0
        DEALLOCATE(md_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(md_env)
END SUBROUTINE md_env_release

! *****************************************************************************
!  SUBROUTINE set_md_env ( md_env, itimes, constant, cell, simpar, pimd_env, &
!                          coef_pos, coef_vel, coef_force, virial, &
!                          nhc_part, nhc_coef, nhc_baro, npt, &
!                          force_env, para_env )
  SUBROUTINE set_md_env ( md_env, itimes, constant, cell, simpar, &
                          nhc_part, nhc_coef, nhc_baro, npt, virial, &
                          force_env, para_env )

!   Purpose: Set the integrator environment to the correct program.

!   ***************************************************************************

    TYPE(md_environment_type), POINTER       :: md_env
    INTEGER, OPTIONAL, POINTER               :: itimes
    REAL(dbl), OPTIONAL, POINTER             :: constant
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(simulation_parameters_type), &
      OPTIONAL, POINTER                      :: simpar
    TYPE(lnhc_parameters_type), OPTIONAL, &
      POINTER                                :: nhc_part( : ), nhc_coef( : ), &
                                                nhc_baro( : )
    TYPE(npt_info_type), OPTIONAL, POINTER   :: npt( :, : )
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_md_env", &
      routineP = moduleN//':'//routine_name

!    TYPE ( pimd_environment_type ), POINTER :: pimd_env

    IF ( PRESENT ( cell ) ) md_env % cell => cell
!    IF ( PRESENT ( pimd_env ) ) md_env % pimd_env => pimd_env
    IF ( PRESENT ( simpar ) ) md_env % simpar => simpar
    IF ( PRESENT ( nhc_part ) ) md_env % nhc_part => nhc_part
    IF ( PRESENT ( nhc_coef ) ) md_env % nhc_coef => nhc_coef
    IF ( PRESENT ( nhc_baro ) ) md_env % nhc_baro => nhc_baro
    IF ( PRESENT ( npt ) ) md_env % npt => npt
    IF ( PRESENT ( virial ) ) md_env % virial => virial
    IF ( PRESENT ( itimes ) ) md_env % itimes => itimes
    IF ( PRESENT ( constant ) ) md_env % constant => constant
    IF ( PRESENT ( force_env ) ) THEN
       CALL force_env_retain ( force_env ) ! accept null pointers?
       CALL force_env_release ( md_env % force_env )
       md_env % force_env => force_env
    END IF
    
         
!   ---------------------------------------------------------------------------

  END SUBROUTINE set_md_env

! *****************************************************************************

  SUBROUTINE zero_virial ( virial )

!   Purpose:  zero the virials

!   ***************************************************************************

    TYPE(virial_type), INTENT(OUT)           :: virial

    CHARACTER(LEN=*), PARAMETER :: routine_name = "zero_virial", &
      routineP = moduleN//':'//routine_name

    virial % pv_total = 0._dbl 
    virial % pv_kinetic = 0._dbl 
    virial % pv_virial = 0._dbl 
    virial % pv_constraint = 0._dbl 
    virial % pv_availability = .TRUE.

! *****************************************************************************
  END SUBROUTINE  zero_virial
! *****************************************************************************

END MODULE md_environment_types
!******************************************************************************
