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
!!
!!   SOURCE
!******************************************************************************

MODULE dump

  USE fist_environment_types, ONLY : replica_environment_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE md_environment_types, ONLY : md_environment_type
  USE nose, ONLY : lnhc_parameters_type, npt_info_type
  USE util, ONLY : get_unit
  USE string_utilities, ONLY : xstring, integer_to_string

  PRIVATE
  PUBLIC :: dump_variables

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** dump/dump_variables [1.0] *
!!
!!   NAME
!!     dump_variables
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE dump_variables ( int_env, rep_env, dump_file_name, globenv, is_pimd )

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
       natoms = size(rep_env(ibead) % part)
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

END SUBROUTINE dump_variables

!!*****
!******************************************************************************

END MODULE dump
