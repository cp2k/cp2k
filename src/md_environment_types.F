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
  USE coefficient_types,               ONLY: coeff_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE force_env_types,                 ONLY: force_env_release,&
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
     TYPE ( coeff_type ), POINTER  :: coef_pos ( : )        
     TYPE ( coeff_type ), POINTER  :: coef_vel ( : )        
     TYPE ( coeff_type ), POINTER  :: coef_force ( : )      
     TYPE ( force_env_type ), POINTER :: force_env  !? should this also be an array for PIMD
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_part ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef 
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
  INTEGER, SAVE, PRIVATE :: last_md_env_id

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

!  SUBROUTINE get_md_env ( md_env, itimes, constant, cell, simpar, pimd_env, &
!                          coef_pos, coef_vel, coef_force, virial, &
!                          nhc_part, nhc_coef, nhc_baro, npt, &
!                          force_env, para_env )
  SUBROUTINE get_md_env ( md_env, itimes, constant, used_time, cell,  &
                          simpar,coef_pos, coef_vel, coef_force, virial, &
                          nhc_part, nhc_coef, nhc_baro, npt, force_env, &
                          para_env )

    IMPLICIT NONE
!   Purpose: Get the integrator environment

!   ***************************************************************************

    TYPE ( md_environment_type ), POINTER :: md_env
    TYPE ( virial_type ), POINTER, OPTIONAL :: virial
    TYPE ( cell_type ), POINTER, OPTIONAL :: cell
    TYPE ( simulation_parameters_type ), POINTER, OPTIONAL :: simpar
    TYPE ( coeff_type ), OPTIONAL, POINTER :: coef_pos ( : )
    TYPE ( coeff_type ), OPTIONAL, POINTER :: coef_vel ( : )
    TYPE ( coeff_type ), OPTIONAL, POINTER :: coef_force ( : )
!    TYPE ( pimd_environment_type ), POINTER :: pimd_env
    TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: nhc_part ( : )
    TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: nhc_coef 
    TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: nhc_baro ( : )
    TYPE ( npt_info_type ), POINTER, OPTIONAL :: npt ( :, : )
    TYPE ( force_env_type ), POINTER, OPTIONAL :: force_env
    INTEGER, POINTER, OPTIONAL :: itimes
    REAL ( dbl ), POINTER, OPTIONAL :: constant
    REAL ( dbl ), POINTER, OPTIONAL :: used_time
    TYPE(cp_para_env_type), POINTER, OPTIONAL :: para_env 
!   ---------------------------------------------------------------------------
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_md_env",&
         routineP=moduleN//':'//routine_name

    IF ( PRESENT ( itimes ) ) itimes => md_env % itimes
    IF ( PRESENT ( constant ) ) constant => md_env % constant
    IF ( PRESENT ( used_time ) ) used_time => md_env % used_time
    IF ( PRESENT ( cell ) ) cell => md_env % cell 
    IF ( PRESENT ( coef_pos ) ) coef_pos => md_env % coef_pos
    IF ( PRESENT ( coef_vel ) ) coef_vel => md_env % coef_vel
    IF ( PRESENT ( coef_force ) ) coef_force => md_env % coef_force
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

    TYPE( md_environment_type ), POINTER :: md_env
    TYPE ( cp_para_env_type ), POINTER :: para_env 
    TYPE ( cp_error_type ), INTENT ( INOUT ), OPTIONAL :: error
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "md_env_create",&
         routineP=moduleN//':'//routine_name
    INTEGER :: istat
    INTEGER, TARGET :: itimes_local = 0
    REAL ( dbl ), TARGET :: constant_local=0._dbl, used_time_local=0._dbl

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
    NULLIFY ( md_env % coef_pos )
    NULLIFY ( md_env % coef_vel )
    NULLIFY ( md_env % coef_force )
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
    constant_local = 0._dbl
    used_time_local = 0._dbl
    md_env % itimes => itimes_local
    md_env % constant => constant_local
    md_env % used_time => used_time_local

  END SUBROUTINE md_env_create
!***************************************************************************

!!****f* md_environment_types/md_env_retain [1.0] *
!!
!!   NAME
!!     md_env_retain
!!
!!   SYNOPSIS
!!     Subroutine md_env_retain(md_env, error)
!!       Type(md_env_type), Pointer:: md_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine md_env_retain
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
  TYPE(md_environment_type), POINTER :: md_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='md_env_retain',&
        routineP=moduleN//':'//routineN

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
!!   SYNOPSIS
!!     Subroutine md_env_release(md_env, error)
!!       Type(md_env_type), Pointer:: md_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine md_env_release
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
  TYPE(md_environment_type), POINTER :: md_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='md_env_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i

  failure=.FALSE.
  
  IF (ASSOCIATED(md_env)) THEN
     CPPreconditionNoFail(md_env%ref_count>0,cp_failure_level,routineP,error)
     md_env%ref_count=md_env%ref_count-1
     IF (md_env%ref_count==0) THEN
        md_env%ref_count=1
        CALL cp_para_env_release ( md_env % para_env, error = error )
        NULLIFY ( md_env % itimes )
        NULLIFY ( md_env % constant )
        NULLIFY ( md_env % cell )
        NULLIFY ( md_env % simpar )
        NULLIFY ( md_env % coef_pos )
        NULLIFY ( md_env % coef_vel )
        NULLIFY ( md_env % coef_force )
        NULLIFY ( md_env % nhc_part )
        NULLIFY ( md_env % nhc_coef )
        NULLIFY ( md_env % nhc_baro )
        NULLIFY ( md_env % npt )
        NULLIFY ( md_env % virial )
        NULLIFY ( md_env % force_env )
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
                          coef_pos, coef_vel, coef_force, virial, &
                          nhc_part, nhc_coef, nhc_baro, npt, &
                          force_env, para_env )

!   Purpose: Set the integrator environment to the correct program.

!   ***************************************************************************

    TYPE ( md_environment_type ), POINTER :: md_env
    TYPE ( virial_type ), POINTER, OPTIONAL :: virial
    TYPE ( cell_type ), POINTER, OPTIONAL :: cell
    TYPE ( simulation_parameters_type ), POINTER, OPTIONAL :: simpar
    TYPE ( coeff_type ), POINTER, OPTIONAL :: coef_pos ( : )
    TYPE ( coeff_type ), POINTER, OPTIONAL :: coef_vel ( : )
    TYPE ( coeff_type ), POINTER, OPTIONAL :: coef_force ( : )
!    TYPE ( pimd_environment_type ), POINTER :: pimd_env
    TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: nhc_part ( : )
    TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: nhc_coef 
    TYPE ( lnhc_parameters_type ), POINTER, OPTIONAL :: nhc_baro ( : )
    TYPE ( npt_info_type ), POINTER, OPTIONAL :: npt ( :, : )
    TYPE ( force_env_type ), POINTER, OPTIONAL :: force_env
    INTEGER, POINTER, OPTIONAL :: itimes
    REAL ( dbl ), POINTER, OPTIONAL :: constant
    TYPE(cp_para_env_type), POINTER, OPTIONAL :: para_env 

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_md_env",&
         routineP=moduleN//':'//routine_name

    INTEGER, TARGET :: itimes_local
    REAL ( dbl ), TARGET :: constant_local, used_time_local

    IF ( PRESENT ( cell ) ) md_env % cell => cell
    IF ( PRESENT ( coef_pos ) ) md_env % coef_pos => coef_pos
    IF ( PRESENT ( coef_vel ) ) md_env % coef_vel => coef_vel
    IF ( PRESENT ( coef_force ) ) md_env % coef_force => coef_force
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

    TYPE (virial_type), INTENT(OUT) :: virial
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "zero_virial",&
         routineP=moduleN//':'//routine_name

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
