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
  USE extended_system_dynamics,        ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE fist_environment_types,          ONLY: replica_environment_type
  USE force_env_types,                   ONLY: force_env_get,&
                                             force_env_type,&
                                             force_env_write_variables
  USE global_types,                    ONLY: global_environment_type
  USE kg_environment_types,            ONLY: kg_environment_type
  USE kinds,                           ONLY: dbl
  USE md_environment_types,            ONLY: md_environment_type
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE simulation_cell,                 ONLY: cell_type
  USE string_utilities,                ONLY: integer_to_string,&
                                             xstring
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_unit
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="dump"

  PUBLIC :: dump_variables

  INTERFACE dump_variables
    MODULE PROCEDURE dump_fist_variables,&
                     dump_kg_variables,&
                     dump_qs_variables,&
                     dump_force_env_variables
  END INTERFACE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** dump/dump_fist_variables [1.0] *
!!
!!   NAME
!!     dump_fist_variables
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************

SUBROUTINE dump_fist_variables ( int_env, rep_env, dump_file_name, globenv, is_pimd )

! routine to dump all simulation variables periodically and at the last step.

  IMPLICIT NONE

! Arguments
  TYPE ( md_environment_type ), INTENT ( IN ) :: int_env
  TYPE ( replica_environment_type ), INTENT ( IN ) :: rep_env ( : )
  CHARACTER ( LEN = * ), INTENT ( IN ) :: dump_file_name
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  LOGICAL, INTENT ( IN ) :: is_pimd

! Locals
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part ( : )
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_baro ( : )
  TYPE ( npt_info_type ), POINTER  :: npt_info ( :, : )
  INTEGER :: i, j, size1, size2, natoms, dmp, ibead, beads, ia, ib
  CHARACTER ( LEN = 80 ) :: file_name
  CHARACTER ( LEN = 5 ) :: tag

!------------------------------------------------------------------------------
! Assigning local pointers
  nhc_part => int_env % nhc_part
  nhc_baro => int_env % nhc_baro
  npt_info => int_env % npt

  beads = SIZE ( rep_env )

  DO ibead = 1, beads
    IF ( globenv % ionode ) THEN
       dmp = get_unit()
       IF ( is_pimd ) THEN
         CALL xstring ( dump_file_name, ia, ib )
         CALL integer_to_string ( ibead, tag )
         file_name = dump_file_name(ia:ib) // '_' // ADJUSTL(tag)
         OPEN ( UNIT = dmp, FILE = file_name, STATUS = 'unknown', &
              POSITION = 'rewind' )
       ELSE
         OPEN ( UNIT = dmp, FILE = dump_file_name, STATUS = 'unknown', &
              POSITION = 'rewind' )
       END IF
    END IF

    IF (globenv % ionode) THEN
       natoms = SIZE(rep_env(ibead) % part)
       WRITE (dmp,*) natoms
! position
       DO i = 1, natoms
          WRITE (dmp,*) rep_env(ibead) % part(i) %r
       END DO
! hmat
       DO i = 1, 3
          WRITE (dmp,*) rep_env(ibead) % box % hmat(i,:)
       END DO
! hmat_ref
       DO i = 1, 3
          WRITE (dmp,*) rep_env(ibead) % box_ref % hmat(i,:)
       END DO
! velocities
       DO i = 1, natoms
          WRITE (dmp,*) rep_env(ibead) % part(i) % v
       END DO
! particle thermostat postions, velocities, masses
       size1 = SIZE ( nhc_part(ibead) % nvt, 1 )
       size2 = SIZE ( nhc_part(ibead) % nvt, 2 )
       WRITE ( dmp, * ) size1, size2
       DO i = 1, size1
         DO j = 1, size2
           WRITE ( dmp, * ) nhc_part(ibead) % nvt ( i, j ) % eta
         ENDDO
       ENDDO
       DO i = 1, size1
         DO j = 1, size2
           WRITE ( dmp, * ) nhc_part(ibead) % nvt ( i, j ) % v
         ENDDO
       ENDDO
       DO i = 1, size1
         DO j = 1, size2
           WRITE ( dmp, * ) nhc_part(ibead) % nvt ( i, j ) % mass
         ENDDO
       ENDDO
! barostat thermostat postions, velocities, masses
       size1 = SIZE ( nhc_baro(ibead) % nvt, 1 )
       size2 = SIZE ( nhc_baro(ibead) % nvt, 2 )
       WRITE ( dmp, * ) size1, size2
       DO i = 1, size1
         DO j = 1, size2
           WRITE ( dmp, * ) nhc_baro(ibead) % nvt ( i, j ) % eta
         ENDDO
       ENDDO
       DO i = 1, size1
         DO j = 1, size2
           WRITE ( dmp, * ) nhc_baro(ibead) % nvt ( i, j ) % v
         ENDDO
       ENDDO
       DO i = 1, size1
         DO j = 1, size2
           WRITE ( dmp, * ) nhc_baro(ibead) % nvt ( i, j ) % mass
         ENDDO
       ENDDO
! barostat velocities, masses
       size1 = SIZE ( npt_info, 1 )
       size2 = SIZE ( npt_info, 2 )
       WRITE ( dmp, * ) size1, size2
       DO i = 1, size1
         DO j = 1, size2
           WRITE ( dmp, * ) npt_info ( i, j ) % v
         ENDDO
       ENDDO
       DO i = 1, size1
         DO j = 1, size2
           WRITE ( dmp, * ) npt_info ( i, j ) % mass
         ENDDO
       ENDDO

       CLOSE (dmp)
    END IF
  END DO

END SUBROUTINE dump_fist_variables

!******************************************************************************
!******************************************************************************
!!****** dump/dump_kg_variables [1.0] *
!!
!!   NAME
!!     dump_kg_variables
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************
SUBROUTINE dump_kg_variables ( int_env, kg_env, dump_file_name, globenv )

! routine to dump all simulation variables periodically and at the last step.

  IMPLICIT NONE

! Arguments
  TYPE ( md_environment_type ), INTENT ( IN ) :: int_env
  TYPE ( kg_environment_type ), INTENT ( IN ) :: kg_env
  CHARACTER ( LEN = * ), INTENT ( IN ) :: dump_file_name
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part ( : )
  TYPE ( npt_info_type ), POINTER  :: npt_info ( :, : )
  INTEGER :: i, j, size1, size2, natoms, dmp, ia, ib
  CHARACTER ( LEN = 80 ) :: file_name
  CHARACTER ( LEN = 5 ) :: tag

!------------------------------------------------------------------------------
! Assigning local pointers
  nhc_part => int_env % nhc_part
  npt_info => int_env % npt

  IF ( globenv % ionode ) THEN
     dmp = get_unit()
     OPEN ( UNIT = dmp, FILE = dump_file_name, STATUS = 'unknown', &
           POSITION = 'rewind' )
  END IF

  IF (globenv % ionode) THEN
     natoms = SIZE(kg_env % part)
     WRITE (dmp,*) natoms
! position
     DO i = 1, natoms
        WRITE (dmp,*) kg_env % part(i) %r
     END DO
! hmat
     DO i = 1, 3
        WRITE (dmp,*) kg_env % box % hmat(i,:)
     END DO
! hmat_ref
     DO i = 1, 3
        WRITE (dmp,*) kg_env % box_ref % hmat(i,:)
     END DO
! velocities
     DO i = 1, natoms
        WRITE (dmp,*) kg_env % part(i) % v
     END DO
! particle thermostat postions, velocities, masses
     size1 = SIZE ( nhc_part(1) % nvt, 1 )
     size2 = SIZE ( nhc_part(1) % nvt, 2 )
     WRITE ( dmp, * ) size1, size2
     DO i = 1, size1
       DO j = 1, size2
         WRITE ( dmp, * ) nhc_part(1) % nvt ( i, j ) % eta
       ENDDO
     ENDDO
     DO i = 1, size1
       DO j = 1, size2
         WRITE ( dmp, * ) nhc_part(1) % nvt ( i, j ) % v
       ENDDO
     ENDDO
     DO i = 1, size1
       DO j = 1, size2
         WRITE ( dmp, * ) nhc_part(1) % nvt ( i, j ) % mass
       ENDDO
     ENDDO
     CLOSE (dmp)
  END IF

END SUBROUTINE dump_kg_variables

!******************************************************************************
!!****** dump/dump_qs_variables [1.0] *
!!
!!   NAME
!!     dump_qs_variables
!!
!!   FUNCTION
!!     Dump all simulation variables periodically and at the last step for QS.
!!
!!   AUTHOR
!!     MK (15.09.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*****************************************************************************

  SUBROUTINE dump_qs_variables (md_env,qs_env,dump_file_name,globenv)

    TYPE(md_environment_type), INTENT(IN)     :: md_env
    TYPE(qs_environment_type), INTENT(IN)     :: qs_env
    CHARACTER(LEN=*), INTENT(IN)              :: dump_file_name
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    TYPE(cell_type), POINTER :: cell

    INTEGER :: dump_unit,i,j,natom,size1,size2

    TYPE(lnhc_parameters_type), DIMENSION(:), POINTER :: nhc_part
    TYPE(particle_type), DIMENSION(:), POINTER        :: particle_set

!   ---------------------------------------------------------------------------

    IF (globenv%ionode) THEN

      CALL open_file(file_name=dump_file_name,&
                     file_action="WRITE",&
                     file_form="FORMATTED",&
                     file_status="REPLACE",&
                     unit_number=dump_unit)

      CALL get_qs_env(qs_env=qs_env,&
                      cell=cell,&
                      particle_set=particle_set)

      natom = SIZE(particle_set)

!MK   WRITE (UNIT=dump_unit,FMT=*) natom

!     *** Coordinates (a.u.) ***

      DO i=1,natom
        WRITE (UNIT=dump_unit,FMT=*) particle_set(i)%r(1:3)
      END DO

!     *** Cell vectors (a.u.) ***

      DO i=1,3
        WRITE (UNIT=dump_unit,FMT=*) cell%hmat(i,1:3)
      END DO

!     *** Reference cell vectors (NYI for QS) ***

      DO i=1,3
        WRITE (UNIT=dump_unit,FMT=*) cell%hmat(i,1:3)
      END DO

!     *** Velocities ***

      DO i=1,natom
        WRITE (UNIT=dump_unit,FMT=*) particle_set(i)%v(1:3)
      END DO

!     *** Particle thermostat positions, velocities, masses ***

      nhc_part => md_env%nhc_part

      IF (ASSOCIATED(nhc_part(1)%nvt)) THEN

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
        ENDDO
      ENDDO

      DO i=1,size1
        DO j=1,size2
          WRITE (UNIT=dump_unit,FMT=*) nhc_part(1)%nvt(i,j)%mass
        END DO
      END DO

      END IF

!     *** Close dump unit and keep it ***

      CALL close_file(unit_number=dump_unit)

    END IF

  END SUBROUTINE dump_qs_variables

! *****************************************************************************

!!****f* dump/dump_force_env_variables [1.0] *
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

  SUBROUTINE dump_force_env_variables(md_env,force_env,dump_file_name,error)

    TYPE(md_environment_type), INTENT(IN)        :: md_env
    TYPE(force_env_type), POINTER                :: force_env
    CHARACTER(LEN=*), INTENT(IN)                 :: dump_file_name
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL :: error

    CHARACTER(LEN=*), PARAMETER :: routineN="dump_force_env_variables",&
                                   routineP=moduleN//":"//routineN

    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(cp_logger_type), POINTER   :: logger,new_logger

    TYPE(cp_error_type) :: new_error
    INTEGER             :: dump_unit,handle,i,j,natom,size1,size2

    TYPE(lnhc_parameters_type), DIMENSION(:), POINTER :: nhc_part

!   ---------------------------------------------------------------------------

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

  END SUBROUTINE dump_force_env_variables

!***************************************************************************

END MODULE dump
