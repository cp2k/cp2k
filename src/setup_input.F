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
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE global_types,                    ONLY: global_environment_type
  USE input_types,                     ONLY: setup_parameters_type
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: pi,&
                                             zero
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE periodic_table,                  ONLY: ptable
  USE string_utilities,                ONLY: make_tuple,&
                                             str_comp,&
                                             str_search,&
                                             uppercase,&
                                             xstring
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
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
!!I cell                                                                      I
!!I end cell                                                                  I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_setup_section ( setup, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: ilen, iw, source, group, icount, nmol_type, ios
  INTEGER :: ierror, i, j, nmol_topo
  LOGICAL :: sim, section_found
  CHARACTER ( LEN = 20 ) :: string, string2
  CHARACTER ( LEN = 6 ) :: label
  CHARACTER ( LEN = 12 ) :: default_thermostat
  CHARACTER ( LEN = 3 ), PARAMETER :: yn ( 0:1 ) = (/ ' NO', 'YES' /)

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
  setup % reference_cell % hmat = 0._dbl
  setup % reference_cell % poisson_solver = "STANDARD"
  setup % reference_cell % perd = 1
  setup % simulation_cell % hmat = 0._dbl
  setup % simulation_cell % poisson_solver = "STANDARD"
  setup % simulation_cell % perd = 1

!..parse the input section
  label = '&SETUP'
! try to find the input section on the main input file
  WRITE(iw,*) "PARSING ",setup%set_file_name
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
        CASE ( 'CELL')
           sim = .TRUE.
           DO
              CALL read_line
              ilen = 12
              CALL get_next ( string2, ilen )
              CALL uppercase ( string2 )
              SELECT CASE ( string2 )
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_setup_section','unknown option')
              CASE ( 'REFERENCE')
                 sim = .FALSE.
              CASE ( 'SIMULATION')
                 sim = .TRUE.
              CASE ( 'PARAMETER')
                 DO i = 1, 3
                   DO j = 1, 3
                     IF ( sim ) THEN
                       CALL get_next(setup % simulation_cell % hmat(i,j))
                     ELSE
                       CALL get_next(setup % reference_cell % hmat(i,j))
                     END IF
                   END DO
                 END DO
              CASE ( 'POISSON_SOLV')
                 ilen=0
                 IF ( sim ) THEN
                   CALL get_next(setup % simulation_cell % poisson_solver,ilen)
                 ELSE
                   CALL get_next(setup % reference_cell % poisson_solver,ilen)
                 END IF
                 CALL uppercase ( setup % reference_cell % poisson_solver )
              CASE ( 'PERIODIC')
                 IF ( sim ) THEN
                   CALL get_next(setup % simulation_cell % perd(1))
                   CALL get_next(setup % simulation_cell % perd(2))
                   CALL get_next(setup % simulation_cell % perd(3))
                 ELSE
                   CALL get_next(setup % reference_cell % perd(1))
                   CALL get_next(setup % reference_cell % perd(2))
                   CALL get_next(setup % reference_cell % perd(3))
                 END IF
              CASE ( 'END')
                 ilen = 4
                 CALL get_next(string2,ilen)
                 EXIT
              END SELECT
           END DO
        END SELECT
        CALL read_line
     END DO

     CALL parser_end
  END IF

! cell setup
!  IF ( SUM ( ABS (setup % reference_cell % hmat) ) == 0._dbl ) THEN
!    setup % reference_cell = setup % simulation_cell
!  ELSE IF ( SUM ( ABS (setup % reference_cell % hmat) ) == 0._dbl ) THEN
!    setup % simulation_cell = setup % reference_cell
!  END IF
!  IF ( SUM ( ABS (setup % reference_cell % hmat) ) == 0._dbl ) THEN
!    CALL stop_program ( "read_setup_section", "No cell information available" )
!  END IF

! write some information to output
  IF ( globenv % ionode ) THEN
     IF (globenv%print_level>=0) THEN
        WRITE ( iw, '()' )
        WRITE ( iw, '( A )' ) ' CELL| Reference Cell '
        WRITE ( iw, '( A, T78, A )' ) ' CELL| Periodic in X direction ', &
             yn(setup % reference_cell % perd(1))
        WRITE ( iw, '( A, T78, A )' ) ' CELL| Periodic in Y direction ', &
             yn(setup % reference_cell % perd(2))
        WRITE ( iw, '( A, T78, A )' ) ' CELL| Periodic in Z direction ', &
             yn(setup % reference_cell % perd(3))
        WRITE ( iw, '( A,T61,A20 )' ) ' CELL| Poisson Solver ', &
            ADJUSTR(setup % reference_cell % poisson_solver)
        WRITE ( iw, '( A,T51,3F10.4 )' ) ' CELL| Parameters  ', &
            setup % reference_cell % hmat(1,1:3)
        WRITE ( iw, '( A,T51,3F10.4 )' ) ' CELL| ', &
            setup % reference_cell % hmat(2,1:3)
        WRITE ( iw, '( A,T51,3F10.4 )' ) ' CELL| ', &
            setup % reference_cell % hmat(3,1:3)
        WRITE ( iw, '()' )
     END IF
  END IF

END SUBROUTINE read_setup_section


END MODULE setup_input

!******************************************************************************
