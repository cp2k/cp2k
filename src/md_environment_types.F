!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
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
!!
!!   SOURCE
!******************************************************************************

MODULE md_environment_types
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE force_env_types,                 ONLY: force_env_release,&
                                             force_env_retain,&
                                             force_env_p_type
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE cell_types,                 ONLY: cell_type
#include "cp_common_uses.h"


  IMPLICIT NONE

  PRIVATE

  ! para_env is the parallel enviroment of the MD,
  ! i.e. the systems that are dealt with by the integrator
  ! e.g in the PIMD this could be parent of every bead

  TYPE md_environment_type
     PRIVATE
     INTEGER :: id_nr, ref_count
     TYPE (cp_para_env_type), POINTER :: para_env 
     INTEGER, POINTER :: itimes, md_approach_nr, accepted_md
     REAL (KIND=dp), POINTER  :: used_time, t
     REAL (KIND=dp), POINTER :: constant
     TYPE ( cell_type ), POINTER :: cell
     TYPE ( force_env_p_type ), POINTER :: force_env_p ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_part ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_coef ( : )
     TYPE ( lnhc_parameters_type ), POINTER :: nhc_baro ( : )
     TYPE ( npt_info_type ), POINTER :: npt ( :, : )
     TYPE ( simulation_parameters_type ), POINTER :: simpar
     TYPE ( reftraj_type ), POINTER :: reftraj_info
  END TYPE md_environment_type

  TYPE reftraj_type
     INTEGER :: unit
  END TYPE

! *** Public subroutines and data types ***
  PUBLIC :: md_environment_type, reftraj_type, set_md_env,   &
            get_md_env, md_env_create,  &
            md_env_release, md_env_retain

! *** Global parameters ***

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='md_environment_types'
  INTEGER, SAVE, PRIVATE :: last_md_env_id=0

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_md_env ( md_env, itimes, constant, used_time, cell,  &
                          simpar, nhc_part, nhc_coef, nhc_baro, &
                          npt, force_env_p, para_env, reftraj_info, t,&
                          md_approach_nr, accepted_md)

!   Purpose: Get the integrator environment

!   ***************************************************************************

    TYPE(md_environment_type), POINTER       :: md_env
    INTEGER, OPTIONAL, POINTER               :: itimes, md_approach_nr, accepted_md
    REAL(KIND=dp), OPTIONAL, POINTER             :: constant, used_time, t
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(simulation_parameters_type), &
      OPTIONAL, POINTER                      :: simpar
    TYPE(lnhc_parameters_type), OPTIONAL, &
      POINTER                                :: nhc_part( : ), nhc_coef( : ), &
                                                nhc_baro( : )
    TYPE(npt_info_type), OPTIONAL, POINTER   :: npt( :, : )
    TYPE(force_env_p_type), OPTIONAL, POINTER  :: force_env_p ( : )
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(reftraj_type), OPTIONAL, POINTER    :: reftraj_info

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_md_env", &
      routineP = moduleN//':'//routine_name

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( itimes ) ) itimes => md_env % itimes
    IF ( PRESENT ( md_approach_nr ) ) md_approach_nr => md_env % md_approach_nr
    IF ( PRESENT ( accepted_md ) ) accepted_md => md_env % accepted_md
    IF ( PRESENT ( constant ) ) constant => md_env % constant
    IF ( PRESENT ( used_time ) ) used_time => md_env % used_time
    IF ( PRESENT ( t ) ) t => md_env % t
    IF ( PRESENT ( cell ) ) cell => md_env % cell 
    IF ( PRESENT ( simpar ) ) simpar => md_env % simpar
    IF ( PRESENT ( nhc_part ) ) nhc_part => md_env % nhc_part 
    IF ( PRESENT ( nhc_coef ) ) nhc_coef => md_env % nhc_coef
    IF ( PRESENT ( nhc_baro ) ) nhc_baro => md_env % nhc_baro
    IF ( PRESENT ( npt ) ) npt => md_env % npt
    IF ( PRESENT ( force_env_p ) ) force_env_p => md_env % force_env_p
    IF ( PRESENT ( para_env ) ) para_env => md_env % para_env
    IF ( PRESENT ( reftraj_info ) ) reftraj_info => md_env % reftraj_info

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

    NULLIFY ( md_env % itimes )
    NULLIFY ( md_env % constant )
    nullify ( md_env % t )
    NULLIFY ( md_env % md_approach_nr )
    NULLIFY ( md_env % accepted_md )
    NULLIFY ( md_env % cell )
    NULLIFY ( md_env % simpar )
    NULLIFY ( md_env % nhc_part )
    NULLIFY ( md_env % nhc_coef )
    NULLIFY ( md_env % nhc_baro )
    NULLIFY ( md_env % npt )
    NULLIFY ( md_env % force_env_p )
    md_env % para_env => para_env 
    CALL cp_para_env_retain ( md_env % para_env )
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
    md_env % itimes    = 0
    md_env % md_approach_nr = 0
    md_env % accepted_md = 0
    md_env % constant  = 0.0_dp
    md_env % used_time = 0.0_dp
    md_env % t         = 0.0_dp
    ALLOCATE ( md_env % reftraj_info )
    md_env % reftraj_info % unit =-1

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

    INTEGER                                  :: istat, stat, i
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
        DEALLOCATE( md_env % reftraj_info , stat=istat )
        CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
        NULLIFY ( md_env % cell )
        NULLIFY ( md_env % simpar )
        NULLIFY ( md_env % nhc_part )
        NULLIFY ( md_env % nhc_coef )
        NULLIFY ( md_env % nhc_baro )
        NULLIFY ( md_env % npt )
        IF (ASSOCIATED(md_env%force_env_p)) THEN
           DO i=1,SIZE(md_env%force_env_p)
              CALL force_env_release(md_env%force_env_p(i)%force_env, error=error)
           END DO
           DEALLOCATE(md_env%force_env_p, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        md_env%ref_count=0
        DEALLOCATE(md_env,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(md_env)
END SUBROUTINE md_env_release

! *****************************************************************************
  SUBROUTINE set_md_env ( md_env, itimes, constant, cell, simpar, &
                          nhc_part, nhc_coef, nhc_baro, &
                          npt, force_env_p, para_env, accepted_md, md_approach_nr )

!   Purpose: Set the integrator environment to the correct program.

!   ***************************************************************************

    TYPE(md_environment_type), POINTER       :: md_env
    INTEGER, OPTIONAL, POINTER               :: itimes, accepted_md, md_approach_nr
    REAL(KIND=dp), OPTIONAL, POINTER             :: constant
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(simulation_parameters_type), &
      OPTIONAL, POINTER                      :: simpar
    TYPE(lnhc_parameters_type), OPTIONAL, &
      POINTER                                :: nhc_part( : ), nhc_coef( : ), &
                                                nhc_baro( : )
    TYPE(npt_info_type), OPTIONAL, POINTER   :: npt( :, : )
    TYPE(force_env_p_type), OPTIONAL,     &
      POINTER                                :: force_env_p ( : )
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env

    INTEGER :: i

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_md_env", &
      routineP = moduleN//':'//routine_name

    IF ( PRESENT ( cell ) ) md_env % cell => cell
    IF ( PRESENT ( simpar ) ) md_env % simpar => simpar
    IF ( PRESENT ( nhc_part ) ) md_env % nhc_part => nhc_part
    IF ( PRESENT ( nhc_coef ) ) md_env % nhc_coef => nhc_coef
    IF ( PRESENT ( nhc_baro ) ) md_env % nhc_baro => nhc_baro
    IF ( PRESENT ( npt ) ) md_env % npt => npt
    IF ( PRESENT ( itimes ) ) md_env % itimes => itimes
    IF ( PRESENT ( accepted_md ) ) md_env % accepted_md => accepted_md
    IF ( PRESENT ( md_approach_nr ) ) md_env % md_approach_nr => md_approach_nr
    IF ( PRESENT ( constant ) ) md_env % constant => constant
    IF ( PRESENT ( force_env_p ) ) THEN
! accept null pointers?
       IF (ASSOCIATED(force_env_p)) THEN
          DO i=1,SIZE(force_env_p)
             CALL force_env_retain ( force_env_p(i) % force_env )
          END DO
       END IF
       IF (ASSOCIATED(md_env%force_env_p)) THEN
          DO i=1,SIZE(md_env%force_env_p)
             CALL force_env_release ( md_env % force_env_p(i) % force_env )
          END DO
       END IF
       md_env % force_env_p => force_env_p
    END IF
    
         
!   ---------------------------------------------------------------------------

  END SUBROUTINE set_md_env

! *****************************************************************************

END MODULE md_environment_types
!******************************************************************************
