!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE atoms_input
  
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE mp, ONLY : mp_bcast
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       cfield, p_error, get_real, get_int
  USE stop_program, ONLY : stop_prg, stop_memory
  USE string_utilities, ONLY : xstring, uppercase
  USE util, ONLY : get_unit
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: read_coord_vel, system_type
  
  TYPE system_type
     REAL ( dbl ) :: box(3,3)
     CHARACTER ( LEN = 3 ) :: ptype
     CHARACTER ( LEN = 20 ) :: rtype
     INTEGER :: n
     REAL ( dbl ), POINTER :: c ( :, : )
     REAL ( dbl ), POINTER :: v ( :, : )
  END TYPE system_type
  
CONTAINS

!!>---------------------------------------------------------------------------!
!!                                                                            !
!!  OPTIONS: either read the section &atoms in the input file or              !
!!           read the coordinates from the file project_name.dat              !
!!  SECTION: &atoms ... &end                                                  !
!!      file:   file_name                                                     !
!!      cell:   b11 b12 b13 &                                                 !
!!              b21 b22 b23 &                                                 !
!!              b31 b32 b33                                                   !
!!                                                                            !
!!                                                                            !
!!                                                                            !
!!<---------------------------------------------------------------------------!

SUBROUTINE read_coord_vel ( atype, filen, globenv )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( system_type ), INTENT ( INOUT ) :: atype
  CHARACTER ( LEN = * ), INTENT ( INOUT ) :: filen
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  
! Locals
  INTEGER :: ierror, ilen, iw, source, allgrp, ia, ie, i, ios
  CHARACTER ( LEN = 6 ) :: string
  CHARACTER ( LEN = 5 ) :: label
  
!------------------------------------------------------------------------------
  
  iw = globenv % scr
     
!..defaults
!..parse the input section
  label = '&ATOMS'
  CALL parser_init ( globenv % input_file_name, label, ierror, globenv )
  IF ( ierror /= 0 ) THEN

     IF( globenv % ionode) THEN
       WRITE ( iw, '( A )' ) ' ATOM| No input section &ATOMS found '
       IF ( filen == ' ' ) THEN
          CALL xstring ( globenv % project_name, ia, ie )
          filen = globenv % project_name ( ia:ie ) // '.dat'
       END IF
       ia = MIN ( LEN ( filen ), 20 )
       WRITE ( iw, '( A, T61, A )' ) ' ATOM| Try to read default file ', &
            ADJUSTR ( filen ( 1:ia ) )
       CALL read_file ( filen, atype )
     END IF
!..broadcast the input data to all nodes
#if defined(__parallel)
     source = globenv % source
     allgrp = globenv % group
     CALL mp_bcast ( atype % box, source, allgrp )
     CALL mp_bcast ( atype % ptype, source, allgrp )
     CALL mp_bcast ( atype % rtype, source, allgrp )
     CALL mp_bcast ( atype % n, source, allgrp )
     IF ( .NOT. globenv % ionode ) THEN
        ALLOCATE ( atype % c(1:3,1:atype % n ), STAT = ios )
        IF ( ios /= 0 ) CALL stop_memory ( 'ATOM', 'atype%c', 3 * atype % n )
     END IF
     CALL mp_bcast ( atype % c, source, allgrp )
     IF ( atype % rtype == 'POSVEL' ) THEN
        IF ( .NOT. globenv % ionode ) THEN
           ALLOCATE (atype % v(1:3,1:atype % n),STAT=ios)
           IF ( ios /= 0 ) CALL stop_memory ( 'ATOM', 'atype%v', 3 * atype % n )
        END IF
        CALL mp_bcast ( atype % v, source, allgrp )
     END IF
#endif
        
     ELSE

        CALL stop_prg ( 'atom','this part of the code not yet written' )
        
     END IF
     CALL parser_end
!..end of parsing the input section
     
!..write some information to output
     IF (globenv%ionode .and. globenv%print_level>0) THEN
        WRITE ( iw, '( A )' ) ' ATOM| Box parameters [Angstrom]'
        WRITE ( iw, '( A, T36, 3F15.5 )' ) &
             ' ATOM| ', ( atype % box ( 1, i ), i = 1, 3 )
        WRITE ( iw, '( A, T36, 3F15.5 )' ) &
             ' ATOM| ', ( atype % box ( 2, i ), i = 1, 3 )
        WRITE ( iw, '( A, T36, 3F15.5 )' ) &
             ' ATOM| ', ( atype % box ( 3, i ), i = 1, 3 )
        WRITE ( iw, '( A, T71, I10 )' ) &
             ' ATOM| Number of atoms read ', atype % n
        IF ( globenv % print_level > 4 ) THEN
           IF ( atype % rtype == 'POS' ) THEN
              CALL print_c ( iw, atype % c )
           ELSE IF ( atype % rtype == 'POSVEL' ) THEN
              CALL print_cv ( iw, atype % c, atype % v )
           END IF
        END IF
        WRITE ( iw, '()' )
     END IF
  
END SUBROUTINE read_coord_vel

!******************************************************************************

SUBROUTINE read_file ( filen, atype )
  
  IMPLICIT NONE
  
! Arguments
  CHARACTER ( LEN = * ) :: filen
  TYPE ( system_type ), INTENT ( INOUT ) :: atype
  
! Locals
  INTEGER :: iunit, i, j, ios
  LOGICAL :: exists
  
!------------------------------------------------------------------------------
  
  INQUIRE ( FILE = filen, EXIST = exists )
  IF ( exists ) THEN
     iunit = get_unit()
     OPEN ( iunit, file = filen )
     READ ( iunit, * ) atype % n
        
     ALLOCATE ( atype % c ( 1:3, 1:atype % n ), STAT = ios )
     IF ( ios /= 0 ) CALL stop_memory ( 'ATOM', 'atype%c', 3 * atype % n )
     NULLIFY ( atype % v )
     
     IF ( atype % rtype == 'POS' ) THEN
        DO i = 1, atype % n
           READ ( iunit, * ) atype % c ( 1:3, i )
        END DO
        READ ( iunit, * ) atype % box ( 1, 1:3 )
        READ ( iunit, * ) atype % box ( 2, 1:3 )
        READ ( iunit, * ) atype % box ( 3, 1:3 )
        
     ELSE IF ( atype % rtype == 'POSVEL' ) THEN
        ALLOCATE ( atype % v ( 1:3, 1:atype % n ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%v', 3 * atype % n )
        
        DO i = 1, atype % n
           READ ( iunit, * ) atype % c ( 1:3, i )
        END DO
        READ ( iunit, * ) atype % box ( 1, 1:3 )
        READ ( iunit, * ) atype % box ( 2, 1:3 )
        READ ( iunit, * ) atype % box ( 3, 1:3 )
        DO i = 1, atype % n
           READ ( iunit, * ) atype % v ( 1:3, i )
        END DO
        
     ELSE
        CALL stop_prg ( 'ATOM', 'this rtype not programmed' )
     END IF

  ELSE
     
     CALL stop_prg ( 'ATOM', 'No information on atoms found ' )

  END IF
  
END SUBROUTINE read_file

!******************************************************************************

SUBROUTINE print_c ( iw, c )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: iw
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: c
  
! Locals
  INTEGER :: i, n
  
!------------------------------------------------------------------------------
  
  n = SIZE ( c, 2 )
  WRITE ( iw, '( A )' ) ' ATOM| Atom coordinates [Angstrom]'
  DO i = 1, n
     WRITE ( iw, '( A, T26, I10, 3F15.5 )' ) ' ATOM| ', i, c ( :, i )
  END DO
  
END SUBROUTINE print_c

!******************************************************************************

SUBROUTINE print_cv ( iw, c, v )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: iw
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: c
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: v
  
! Locals
  INTEGER :: i, n
  
!------------------------------------------------------------------------------
  
  n = SIZE ( c, 2 )
  WRITE ( iw, '( A )' ) ' ATOM| Atom coordinates [Angstrom]'
  DO i = 1, n
     WRITE ( iw, '( A, T8, I8, 3F10.4, 5X, 3F10.4 )' ) &
          ' ATOM| ', i, c ( :, i ), v ( :, i )
  END DO
  
END SUBROUTINE print_cv

!******************************************************************************

END MODULE atoms_input
