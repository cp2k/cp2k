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
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_write_variables
  USE kinds,                           ONLY: dbl
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE timings,                         ONLY: timeset,&
                                             timestop

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

    TYPE(md_environment_type), POINTER       :: md_env
    CHARACTER(LEN=*), INTENT(IN)             :: dump_file_name
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "dump_force_env_variables", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: dump_unit, handle, i, j, size1, size2
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger, new_logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(force_env_p_type), DIMENSION(:), POINTER  :: force_env_p
    TYPE(lnhc_parameters_type), POINTER      :: nhc_part( : ), nhc_coef( : ), &
                                                nhc_baro( : )
    TYPE(npt_info_type), POINTER   :: npt( :, : )

!   ---------------------------------------------------------------------------

    CALL get_md_env ( md_env, force_env_p = force_env_p, nhc_part = nhc_part, &
                      nhc_baro = nhc_baro, npt = npt )
    CALL timeset(routineN,"I","",handle)

    NULLIFY (para_env,new_logger)

    CALL get_md_env ( md_env, para_env = para_env )

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

    IF (ASSOCIATED(force_env_p)) THEN
       DO i = 1, SIZE(force_env_p)
          CALL force_env_write_variables(force_env_p(i)%force_env, &
                                           unit_nr=dump_unit,error=new_error)
       END DO
    END IF

    IF (para_env%mepos==para_env%source) THEN

!     *** Particle thermostat positions, velocities, masses ***



       size1 = SIZE(nhc_part(1)%nvt,1)
       size2 = SIZE(nhc_part(1)%nvt,2)

       WRITE (UNIT=dump_unit,FMT=*) size1,size2

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_part(1)%nvt(i,j)%eta
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_part(1)%nvt(i,j)%v
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_part(1)%nvt(i,j)%f
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_part(1)%nvt(i,j)%mass
         END DO
       END DO


       size1 = SIZE(nhc_baro(1)%nvt,1)
       size2 = SIZE(nhc_baro(1)%nvt,2)

       WRITE (UNIT=dump_unit,FMT=*) size1,size2

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_baro(1)%nvt(i,j)%eta
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_baro(1)%nvt(i,j)%v
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_baro(1)%nvt(i,j)%f
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) nhc_baro(1)%nvt(i,j)%mass
         END DO
       END DO

!     *** barostat velocities, masses ***
       size1 = SIZE(npt,1)
       size2 = SIZE(npt,2)

       WRITE (UNIT=dump_unit,FMT=*) size1,size2

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) npt(i,j)%v
         END DO
       END DO

       DO i=1,size1
         DO j=1,size2
           WRITE (UNIT=dump_unit,FMT=*) npt(i,j)%mass
         END DO
       END DO

!MK!     *** Close dump unit and keep it ***

      CALL close_file (unit_number=dump_unit)

    END IF

    CALL cp_error_dealloc_ref(new_error)
    CALL cp_logger_release(new_logger)

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE dump_variables

!***************************************************************************

END MODULE dump
