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

  USE kinds, ONLY : dbl
  USE global_types, ONLY : global_environment_type
  USE nose, ONLY : lnhc_parameters_type, npt_info_type
  USE structure_types, ONLY : structure_type
  USE util, ONLY : get_unit

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

SUBROUTINE dump_variables ( struc, nhc_part, nhc_baro, npt_info,  &
                            dump_file_name, globenv )

! routine to dump all simulation variables periodically and at the last step.

  IMPLICIT NONE

! Arguments
  TYPE ( structure_type ), INTENT ( IN ) :: struc
  CHARACTER ( LEN = * ), INTENT ( IN ) :: dump_file_name
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE ( lnhc_parameters_type ), INTENT ( IN ) :: nhc_part, nhc_baro
  TYPE ( npt_info_type ), DIMENSION ( :, : ), INTENT ( IN ) :: npt_info

! Locals
  INTEGER :: i, j, size1, size2, natoms, dmp

!------------------------------------------------------------------------------

  IF ( globenv % ionode ) THEN
     dmp = get_unit()
     OPEN ( UNIT = dmp, FILE = dump_file_name, STATUS = 'unknown', &
          POSITION = 'rewind' )
  END IF

  IF (globenv % ionode) THEN
     natoms = size(struc % part)
     WRITE (dmp,*) natoms
! position
     DO i = 1, natoms
        WRITE (dmp,*) struc % part(i) %r
     END DO
! hmat
     DO i = 1, 3
        WRITE (dmp,*) struc % box % hmat(i,:)
     END DO
! hmat_ref
     DO i = 1, 3
        WRITE (dmp,*) struc % box_ref % hmat(i,:)
     END DO
! velocities
     DO i = 1, natoms
        WRITE (dmp,*) struc % part(i) % v
     END DO
! particle thermostat postions, velocities, masses
     size1 = SIZE ( nhc_part % nvt, 1 )
     size2 = SIZE ( nhc_part % nvt, 2 )
     WRITE ( dmp, * ) size1, size2
     DO i = 1, size1
       DO j = 1, size2
         WRITE ( dmp, * ) nhc_part % nvt ( i, j ) % eta
       ENDDO
     ENDDO
     DO i = 1, size1
       DO j = 1, size2
         WRITE ( dmp, * ) nhc_part % nvt ( i, j ) % v
       ENDDO
     ENDDO
     DO i = 1, size1
       DO j = 1, size2
         WRITE ( dmp, * ) nhc_part % nvt ( i, j ) % mass
       ENDDO
     ENDDO
! barostat thermostat postions, velocities, masses
     size1 = SIZE ( nhc_baro % nvt, 1 )
     size2 = SIZE ( nhc_baro % nvt, 2 )
     WRITE ( dmp, * ) size1, size2
     DO i = 1, size1
       DO j = 1, size2
         WRITE ( dmp, * ) nhc_baro % nvt ( i, j ) % eta
       ENDDO
     ENDDO
     DO i = 1, size1
       DO j = 1, size2
         WRITE ( dmp, * ) nhc_baro % nvt ( i, j ) % v
       ENDDO
     ENDDO
     DO i = 1, size1
       DO j = 1, size2
         WRITE ( dmp, * ) nhc_baro % nvt ( i, j ) % mass
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

END SUBROUTINE dump_variables

!!*****
!******************************************************************************

END MODULE dump
