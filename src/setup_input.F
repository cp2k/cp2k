!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/setup_input [1.0] *
!!
!!   NAME
!!     setup_input
!!
!!   FUNCTION
!!     Reads the input sections "setup" and "molecule"
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) Added read_topology_section
!!
!!   SOURCE
!******************************************************************************

MODULE setup_input
  USE global_types,                    ONLY: global_environment_type
  USE input_types,                     ONLY: setup_parameters_type
  USE kinds,                           ONLY: dbl
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE string_utilities,                ONLY: uppercase

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_setup_section

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** setup_input/read_setup_section [1.0] *
!!
!!   NAME
!!     read_setup_section
!!
!!   SYNOPSIS
!!     Subroutine read_setup_section(mol_setup, setup, globenv)
!!       Implicit None
!!       Type(molecule_type), Dimension(:), Pointer:: mol_setup
!!       Type(setup_parameters_type):: setup
!!       Type(global_environment_type), Intent (IN):: globenv
!!     End Subroutine read_setup_section
!!
!!   FUNCTION
!!     reads the input section setup
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (30.11.2001) : first search on input file, then on set file
!!                        additional input (charge, multiplicity etc)
!!     JGH (01.12.2001) : new keyword to set default for thermostat option
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &setup ... &end                                                  I
!!I                                                                           I
!!I charge           charge                                                   I
!!I multiplicity     multiplicity                                             I
!!I spin_moment      spinmoment                                               I
!!I states           states                                                   I
!!I spins            nspins                                                   I
!!I LSD              nspins=2                                                 I
!!I etemp            etemp                                                    I
!!I broadening       bmethod                                                  I
!!I thermostats      default_thermostat                                       I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_setup_section ( setup, globenv )


    TYPE(setup_parameters_type)              :: setup
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=3), PARAMETER              :: yn( 0:1 ) = (/ ' NO', 'YES' /)

    CHARACTER(LEN=12)                        :: default_thermostat
    CHARACTER(LEN=20)                        :: string
    CHARACTER(LEN=6)                         :: label
    INTEGER                                  :: icount, ierror, ilen, iw
    LOGICAL                                  :: section_found

  default_thermostat = "GLOBAL"
  iw = globenv % scr
  icount = 0
!
  setup % states = -1
  setup % nspins = 1
  setup % multiplicity = 1
  setup % spinmoment = 0._dbl
  setup % etemp = 0._dbl
  setup % charge = 0._dbl
  setup % broadening_method = ""

!..parse the input section
  label = '&SETUP'
! try to find the input section on the main input file
  IF (globenv%ionode) WRITE(iw,*) "PARSING ",setup%set_file_name
  CALL parser_init ( setup % set_file_name, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
     IF ( globenv % ionode .AND. globenv % print_level >= 2 ) THEN
        WRITE ( iw, '( A )' ) ' No input section &SETUP found on files '
        WRITE ( iw, '( T2, A )' ) setup % set_file_name
     END IF
     section_found = .FALSE.
  ELSE
     section_found = .TRUE.
  END IF
  CALL parser_end

  IF ( section_found ) THEN
! reinitiate the parser at the beginning of the section
     CALL parser_init ( setup % set_file_name, globenv )
     CALL search_label ( label, ierror, ignore_case=.TRUE. )
     IF ( ierror /= 0 ) CALL stop_parser ('read_setup_section','reinit failed')
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 8
        CALL get_next ( string, ilen )
        CALL uppercase ( string )
        SELECT CASE (string)
        CASE DEFAULT
           CALL p_error()
           CALL stop_parser ( 'read_setup_section','unknown option')
        CASE ( 'CHARGE')
           CALL get_next ( setup % charge )
        CASE ( 'MULTIPLI')
           CALL get_next ( setup % multiplicity )
        CASE ( 'SPIN_MOM')
           CALL get_next ( setup % spinmoment )
        CASE ( 'STATES')
           CALL get_next ( setup % states )
        CASE ( 'SPINS')
           CALL get_next ( setup % nspins )
        CASE ( 'LSD')
           setup % nspins = 2
        CASE ( 'ETEMP')
           CALL get_next ( setup % etemp )
        CASE ( 'BROADENI')
           ilen = 0
           CALL get_next(setup % broadening_method,ilen)
        CASE ( 'THERMOST')
           ilen = 0
           CALL get_next(default_thermostat,ilen)
        END SELECT
        CALL read_line
     END DO

     CALL parser_end
  END IF

! write some information to output
  IF ( globenv % ionode ) THEN
     IF (globenv%print_level>=0) THEN
     END IF
  END IF

END SUBROUTINE read_setup_section


END MODULE setup_input

!******************************************************************************
