!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE dump
  
  USE kinds, ONLY : dbl
  USE global_types, ONLY : global_environment_type
  USE stop_program, ONLY : stop_prg
  USE structure_types, ONLY : structure_type
  USE util, ONLY : get_unit
  
  PRIVATE
  PUBLIC :: dump_variables
  
CONTAINS

!******************************************************************************

SUBROUTINE dump_variables ( struc, dump_file_name, globenv )
  
! routine to dump all simulation variables periodically and at the last step.
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( structure_type ), INTENT ( IN ) :: struc
  CHARACTER ( LEN = * ), INTENT ( IN ) :: dump_file_name
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  
! Locals
  INTEGER :: i, natoms, dmp
  
!------------------------------------------------------------------------------
  
  IF ( globenv % ionode ) THEN
     dmp = get_unit()
     OPEN ( UNIT = dmp, FILE = dump_file_name, STATUS = 'unknown', &
          POSITION = 'rewind' )
  END IF
  
  IF (globenv % ionode) THEN
     natoms = size(struc % part)
     WRITE (dmp,*) natoms
     DO i = 1, natoms
        WRITE (dmp,*) struc % part(i) %r
     END DO
     DO i = 1, 3
        WRITE (dmp,*) struc % box % hmat(i,:)
     END DO
     DO i = 1, natoms
        WRITE (dmp,*) struc % part(i) % v
     END DO
     CLOSE (dmp)
  END IF
  
END SUBROUTINE dump_variables

!******************************************************************************

END MODULE dump
