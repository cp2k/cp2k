!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/atoms_input [1.0] *
!!
!!   NAME
!!     atoms_input
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!   cjm, Feb-20-2001 : added all the extended variables to
!!   system_type
!!
!!   SOURCE
!******************************************************************************

MODULE atoms_input

  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE message_passing, ONLY : mp_bcast
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       p_error, get_next
  USE string_utilities, ONLY : xstring, uppercase
  USE termination, ONLY : stop_memory, stop_program
  USE util, ONLY : get_unit

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_coord_vel, system_type

  TYPE system_type
     REAL ( dbl ) :: box ( 3, 3 ), box_ref ( 3, 3 )
     CHARACTER ( LEN = 3 ) :: ptype
     CHARACTER ( LEN = 20 ) :: rtype
     INTEGER :: n, size1_nvt_part, size2_nvt_part
     INTEGER :: size1_nvt_baro, size2_nvt_baro
     INTEGER :: size1_npt, size2_npt
     REAL ( dbl ), POINTER :: c ( :, : )
     REAL ( dbl ), POINTER :: v ( :, : )
     REAL ( dbl ), POINTER :: veps ( :, : )
     REAL ( dbl ), POINTER :: mass_eps ( :, : )
     REAL ( dbl ), POINTER :: eta_part ( :, : )
     REAL ( dbl ), POINTER :: veta_part ( :, : )
     REAL ( dbl ), POINTER :: mass_eta_part ( :, : )
     REAL ( dbl ), POINTER :: eta_baro ( :, : )
     REAL ( dbl ), POINTER :: veta_baro ( :, : )
     REAL ( dbl ), POINTER :: mass_eta_baro ( :, : )
  END TYPE system_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** atoms_input/read_coord_vel [1.0] *
!!
!!   NAME
!!     read_coord_vel
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!   MODIFICATION HISTORY
!!   cjm, Feb-20-2001 : broadcasts extended variables in
!!   system_type
!!
!!  OPTIONS
!!    either read the section &atoms in the input file or
!!    read the coordinates from the file project_name.dat
!!
!!  INPUTS
!!I---------------------------------------------------------------------------I
!!I SECTION: &atoms ... &end                                                  I
!!I     file:   file_name                                                     I
!!I     cell:   b11 b12 b13 &                                                 I
!!I             b21 b22 b23 &                                                 I
!!I             b31 b32 b33                                                   I
!!I---------------------------------------------------------------------------I
!!   SOURCE
!******************************************************************************

SUBROUTINE read_coord_vel ( atype, filen, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( system_type ), INTENT ( INOUT ) :: atype
  CHARACTER ( LEN = * ), INTENT ( INOUT ) :: filen
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: ierror, ilen, iw, source, allgrp, ia, ie, i, ios
  INTEGER :: size1, size2
  CHARACTER ( LEN = 6 ) :: string
  CHARACTER ( LEN = 5 ) :: label

!------------------------------------------------------------------------------

  iw = globenv % scr

!..defaults
  NULLIFY ( atype % c )
  NULLIFY ( atype % v )
  NULLIFY ( atype % veps )
  NULLIFY ( atype % mass_eps )
  NULLIFY ( atype % eta_part )
  NULLIFY ( atype % veta_part )
  NULLIFY ( atype % mass_eta_part )
  NULLIFY ( atype % eta_baro )
  NULLIFY ( atype % veta_baro )
  NULLIFY ( atype % mass_eta_baro )

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
     CALL mp_bcast ( atype % box_ref, source, allgrp )
     CALL mp_bcast ( atype % ptype, source, allgrp )
     CALL mp_bcast ( atype % rtype, source, allgrp )
     CALL mp_bcast ( atype % n, source, allgrp )
     IF ( .NOT. ASSOCIATED ( atype % c ) ) THEN
        ALLOCATE ( atype % c(1:3,1:atype % n ), STAT = ios )
        IF ( ios /= 0 ) CALL stop_memory ( 'ATOM', 'atype%c', 3 * atype % n )
     END IF
     CALL mp_bcast ( atype % c, source, allgrp )
     IF ( atype % rtype == 'POSVEL' ) THEN
        IF ( .NOT. ASSOCIATED ( atype % v ) ) THEN
           ALLOCATE (atype % v(1:3,1:atype % n),STAT=ios)
           IF ( ios /= 0 ) CALL stop_memory ( 'ATOM', 'atype%v', 3 * atype % n )
        END IF
        CALL mp_bcast ( atype % v, source, allgrp )
     END IF
     IF ( atype % rtype == 'ALL' ) THEN
        IF ( .NOT. globenv % ionode ) THEN
           ALLOCATE (atype % v(1:3,1:atype % n),STAT=ios)
           IF ( ios /= 0 ) CALL stop_memory ( 'ATOM', 'atype%v', 3 * atype % n )
           size1 = atype % size1_nvt_part
           size2 = atype % size2_nvt_part
           ALLOCATE (atype % veta_part (size1, size2 ),STAT=ios)
           IF ( ios /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%veta_part', size1 * size2 )
           ALLOCATE (atype % eta_part (size1, size2 ),STAT=ios)
           IF ( ios /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%eta_part', size1 * size2 )
           ALLOCATE (atype % mass_eta_part (size1, size2 ),STAT=ios)
           IF ( ios /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%mass_eta_part', size1 * size2 )
           ALLOCATE (atype % veta_baro (size1, size2 ),STAT=ios)
           size1 = atype % size1_nvt_baro
           size2 = atype % size2_nvt_baro
           IF ( ios /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%veta_baro', size1 * size2 )
           ALLOCATE (atype % eta_baro (size1, size2 ),STAT=ios)
           IF ( ios /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%eta_baro', size1 * size2 )
           ALLOCATE (atype % mass_eta_baro (size1, size2 ),STAT=ios)
           IF ( ios /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%mass_eta_baro', size1 * size2 )
           size1 = atype % size1_npt
           size2 = atype % size2_npt
           ALLOCATE (atype % veps (size1, size2 ),STAT=ios)
           IF ( ios /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%veps', size1 * size2 )
           ALLOCATE (atype % mass_eps (size1, size2 ),STAT=ios)
           IF ( ios /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%mass_eps', size1 * size2 )
        END IF
        CALL mp_bcast ( atype % v, source, allgrp )
        CALL mp_bcast ( atype % veta_part, source, allgrp )
        CALL mp_bcast ( atype % eta_part, source, allgrp )
        CALL mp_bcast ( atype % veta_baro, source, allgrp )
        CALL mp_bcast ( atype % eta_baro, source, allgrp )
        CALL mp_bcast ( atype % veps, source, allgrp )
        CALL mp_bcast ( atype % mass_eps, source, allgrp )
        CALL mp_bcast ( atype % mass_eta_part, source, allgrp )
        CALL mp_bcast ( atype % mass_eta_baro, source, allgrp )
     END IF
#endif

     ELSE

        CALL stop_program ( 'atom','this part of the code not yet written' )

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

!!*****
!******************************************************************************

SUBROUTINE read_file ( filen, atype )

  IMPLICIT NONE

! Arguments
  CHARACTER ( LEN = * ) :: filen
  TYPE ( system_type ), INTENT ( INOUT ) :: atype

! Locals
  INTEGER :: iunit, i, j, ios, size1, size2
  LOGICAL :: exists

!------------------------------------------------------------------------------

  INQUIRE ( FILE = filen, EXIST = exists )
  IF ( exists ) THEN
     iunit = get_unit()
     OPEN ( iunit, file = filen )
     READ ( iunit, * ) atype % n

     IF (.NOT. ASSOCIATED (atype % c) ) THEN
        ALLOCATE ( atype % c ( 1:3, 1:atype % n ), STAT = ios )
        IF ( ios /= 0 ) CALL stop_memory ( 'ATOM', 'atype%c', 3 * atype % n )
        NULLIFY ( atype % v )
     ELSE IF (atype % n /= SIZE( atype % c, 2 ) ) THEN
        CALL stop_program ( 'atoms_input','Inconsistent number of atoms' )
     END IF

     IF ( atype % rtype == 'INIT' ) THEN
        DO i = 1, atype % n
           READ ( iunit, * ) atype % c ( 1:3, i )
        END DO
        READ ( iunit, * ) atype % box ( 1, 1:3 )
        READ ( iunit, * ) atype % box ( 2, 1:3 )
        READ ( iunit, * ) atype % box ( 3, 1:3 )

     ELSE IF ( atype % rtype == 'POS' )THEN
        DO i = 1, atype % n
           READ ( iunit, * ) atype % c ( 1:3, i )
        END DO
        READ ( iunit, * ) atype % box ( 1, 1:3 )
        READ ( iunit, * ) atype % box ( 2, 1:3 )
        READ ( iunit, * ) atype % box ( 3, 1:3 )

        READ ( iunit, * ) atype % box_ref ( 1, 1:3 )
        READ ( iunit, * ) atype % box_ref ( 2, 1:3 )
        READ ( iunit, * ) atype % box_ref ( 3, 1:3 )
     ELSE IF ( atype % rtype == 'POSVEL' ) THEN
        IF (.NOT. ASSOCIATED (atype % v) ) THEN
           ALLOCATE ( atype % v ( 1:3, 1:atype % n ), STAT = ios )
           IF ( ios /= 0 ) &
                CALL stop_memory ( 'ATOM', 'atype%v', 3 * atype % n )
        END IF

        DO i = 1, atype % n
           READ ( iunit, * ) atype % c ( 1:3, i )
        END DO
        READ ( iunit, * ) atype % box ( 1, 1:3 )
        READ ( iunit, * ) atype % box ( 2, 1:3 )
        READ ( iunit, * ) atype % box ( 3, 1:3 )

        READ ( iunit, * ) atype % box_ref ( 1, 1:3 )
        READ ( iunit, * ) atype % box_ref ( 2, 1:3 )
        READ ( iunit, * ) atype % box_ref ( 3, 1:3 )
        DO i = 1, atype % n
           READ ( iunit, * ) atype % v ( 1:3, i )
        END DO
     ELSE IF ( atype % rtype == 'ALL' ) THEN
        IF (.NOT. ASSOCIATED (atype % v) ) THEN
           ALLOCATE ( atype % v ( 1:3, 1:atype % n ), STAT = ios )
           IF ( ios /= 0 ) &
                CALL stop_memory ( 'ATOM', 'atype%v', 3 * atype % n )
        END IF

        DO i = 1, atype % n
           READ ( iunit, * ) atype % c ( 1:3, i )
        END DO
        READ ( iunit, * ) atype % box ( 1, 1:3 )
        READ ( iunit, * ) atype % box ( 2, 1:3 )
        READ ( iunit, * ) atype % box ( 3, 1:3 )

        READ ( iunit, * ) atype % box_ref ( 1, 1:3 )
        READ ( iunit, * ) atype % box_ref ( 2, 1:3 )
        READ ( iunit, * ) atype % box_ref ( 3, 1:3 )
        DO i = 1, atype % n
           READ ( iunit, * ) atype % v ( 1:3, i )
        END DO
! particle thermostat positions and velocities
        READ ( iunit, * ) atype % size1_nvt_part, atype % size2_nvt_part
        size1 = atype % size1_nvt_part
        size2 = atype % size2_nvt_part
        ALLOCATE ( atype % veta_part ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%veta_part', size1 * size2 )
        ALLOCATE ( atype % eta_part ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%eta_part', size1 * size2 )
        ALLOCATE ( atype % mass_eta_part ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%mass_eta_part', size1 * size2 )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % eta_part ( i, j )
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % veta_part ( i, j )
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % mass_eta_part ( i, j )
          ENDDO
        ENDDO
! barostat thermostat positions and velocities
        READ ( iunit, * ) atype % size1_nvt_baro, atype % size2_nvt_baro
        size1 = atype % size1_nvt_baro
        size2 = atype % size2_nvt_baro
        ALLOCATE ( atype % veta_baro ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%veta_baro', size1 * size2 )
        ALLOCATE ( atype % eta_baro ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%eta_baro', size1 * size2 )
        ALLOCATE ( atype % mass_eta_baro ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%mass_eta_baro', size1 * size2 )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % eta_baro ( i, j )
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % veta_baro ( i, j )
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % mass_eta_baro ( i, j )
          ENDDO
        ENDDO
! barostat velocities
        READ ( iunit, * ) atype % size1_npt, atype % size2_npt
        size1 = atype % size1_npt
        size2 = atype % size2_npt
        ALLOCATE ( atype % veps ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%veps', size1 * size2 )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % veps ( i, j )
          ENDDO
        ENDDO
        ALLOCATE ( atype % mass_eps ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%mass_eps', size1 * size2 )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % mass_eps ( i, j )
          ENDDO
        ENDDO
     ELSE
        CALL stop_program ( 'ATOM', 'this rtype not programmed' )
     END IF

     CLOSE ( iunit )

  ELSE

     CALL stop_program ( 'ATOM', 'No information on atoms found ' )

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
