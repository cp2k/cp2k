!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/dump [1.0] *
!!
!!   NAME
!!     dump
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     CJM 2/01:  Now dumps out all extended_system variables
!!     and the box_ref
!!     - Simple dump for QS without extended variables (MK,15.09.2003)
!!
!!   SOURCE
!******************************************************************************

MODULE dump
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_create,&
                                             cp_logger_release,&
                                             cp_logger_type,&
                                             cp_to_string
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE force_env_types,                   ONLY: force_env_get,&
                                             force_env_type,&
                                             force_env_write_variables
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE md_environment_types,            ONLY: md_environment_type, &
                                             get_md_env
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_unit
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="dump"

  PUBLIC :: dump_variables

!!*****
!******************************************************************************

CONTAINS

! *****************************************************************************

!!****f* dump/dump_variables [1.0] *
!!
!!   NAME
!!     dump_force_env_variables
!!
!!   FUNCTION
!!     dumps the variables 
!!
!!   NOTES
!!     this call call the right dump_*_variables, but I find it more clean
!!     to call the force_env routine
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************

  SUBROUTINE dump_variables(md_env,dump_file_name,error)

    TYPE(md_environment_type), POINTER        :: md_env
    CHARACTER(LEN=*), INTENT(IN)                 :: dump_file_name
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL :: error

    CHARACTER(LEN=*), PARAMETER :: routineN="dump_force_env_variables",&
                                   routineP=moduleN//":"//routineN

    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(cp_logger_type), POINTER   :: logger,new_logger

    TYPE(cp_error_type) :: new_error
    INTEGER             :: dump_unit,handle,i,j,natom,size1,size2

    TYPE(lnhc_parameters_type), DIMENSION(:), POINTER :: nhc_part
    TYPE(force_env_type), POINTER :: force_env

!   ---------------------------------------------------------------------------

    CALL get_md_env ( md_env, force_env = force_env )
    CALL timeset(routineN,"I","",handle)

    NULLIFY (nhc_part,para_env,new_logger)

    CALL force_env_get(force_env,para_env=para_env)

    logger => cp_error_get_logger(error)
    CALL cp_logger_create(new_logger,para_env=para_env,template_logger=logger)
    CALL cp_error_init(new_error,logger=new_logger,template_error=error)

    IF (para_env%mepos==para_env%source) THEN
      CALL open_file(file_name=dump_file_name,&
                     file_action="WRITE",&
                     file_form="FORMATTED",&
                     file_status="REPLACE",&
                     unit_number=dump_unit)
    END IF

    CALL force_env_write_variables(force_env,unit_nr=dump_unit,error=new_error)

    IF (para_env%mepos==para_env%source) THEN

!     *** Particle thermostat positions, velocities, masses ***

!MK      nhc_part => md_env%nhc_part
!MK
!MK      IF (ASSOCIATED(nhc_part)) THEN
!MK
!MK        size1 = SIZE(nhc_part(1)%nvt,1)
!MK        size2 = SIZE(nhc_part(1)%nvt,2)
!MK
!MK        WRITE (UNIT=dump_unit,FMT=*) size1,size2
!MK
!MK        DO i=1,size1
!MK          DO j=1,size2
!MK            WRITE (UNIT=dump_unit,FMT=*) nhc_part(1)%nvt(i,j)%eta
!MK          END DO
!MK        END DO
!MK
!MK        DO i=1,size1
!MK          DO j=1,size2
!MK            WRITE (UNIT=dump_unit,FMT=*) nhc_part(1)%nvt(i,j)%v
!MK          END DO
!MK        END DO
!MK
!MK        DO i=1,size1
!MK          DO j=1,size2
!MK            WRITE (UNIT=dump_unit,FMT=*) nhc_part(1)%nvt(i,j)%mass
!MK          END DO
!MK        END DO
!MK
!MK      END IF
!MK
!MK!     *** Close dump unit and keep it ***

      CALL close_file (unit_number=dump_unit)

    END IF

    CALL cp_error_dealloc_ref(new_error)
    CALL cp_logger_release(new_logger)

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE dump_variables

!***************************************************************************

END MODULE dump
