!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/molecule_input [1.0] *
!!
!!   NAME
!!     molecule_input
!!
!!   FUNCTION
!!     Reads the input sections "setup" and "molecule"
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE molecule_input

  USE global_types, ONLY : global_environment_type
  USE input_types, ONLY : setup_parameters_type
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi, zero
  USE molecule_types, ONLY : molecule_type
  USE message_passing, ONLY : mp_bcast
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       p_error, get_next, stop_parser, search_label
  USE periodic_table, ONLY : ptable
  USE string_utilities, ONLY : uppercase, xstring, str_search, &
       str_comp, make_tuple
  USE termination, ONLY : stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_setup_section, read_molecule_section

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** molecule_input/read_setup_section [1.0] *
!!
!!   NAME
!!     read_setup_section
!!
!!   FUNCTION
!!     reads the input section setup
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (30.11.2001) : automatic detection of nmol_type
!!                        first search on input file, then on set file
!!                        added short form for molecule and atom
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
!!I atom             name numat {nhcopt}                                      I
!!I molecule         name numol {nhcopt}                                      I
!!I molecule                                                                  I
!!I    name          molname                                                  I
!!I    nummol        nmol                                                     I
!!I    nhcopt        [molecule,global,massive,...]                            I
!!I end molecule                                                              I
!!I cell                                                                      I
!!I end cell                                                                  I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_setup_section ( mol_setup, setup, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: ilen, iw, source, group, icount, nmol_type, ios
  INTEGER :: ierror, ierr2, i, j
  LOGICAL :: sim
  CHARACTER ( LEN = 20 ) :: string, string2
  CHARACTER ( LEN = 6 ) :: label
  CHARACTER ( LEN = 12 ) :: default_thermostat
  CHARACTER ( LEN = 3 ), PARAMETER :: yn ( 0:1 ) = (/ ' NO', 'YES' /)

!------------------------------------------------------------------------------

  default_thermostat = "GLOBAL"
  iw = globenv % scr
  icount = 0
  nmol_type = -1
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
  CALL parser_init ( setup % set_file_name, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
     IF ( globenv % ionode .AND. globenv % print_level >= 2 ) THEN
        WRITE ( iw, '( A )' ) ' No input section &SETUP found on files '
        WRITE ( iw, '( T2, A )' ) setup % set_file_name
     END IF
  ELSE
!    count number of molecule entries
     nmol_type = 0
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 8
        CALL get_next ( string, ilen )
        CALL uppercase ( string )
        SELECT CASE (string)
        CASE ( 'MOLECULE')
           nmol_type = nmol_type + 1
           IF ( test_next() /= 'C' ) THEN
              DO
                 CALL read_line
                 ilen = 6
                 CALL get_next ( string2, ilen )
                 CALL uppercase ( string2 )
                 IF ( string2 == "END" ) THEN
                    ilen = 8
                    CALL get_next(string2,ilen)
                    EXIT
                 END IF
              END DO
           END IF
        CASE ( 'ATOM')
           nmol_type = nmol_type + 1
        END SELECT
        CALL read_line
     END DO
!    allocate space for the molecules
     ALLOCATE (mol_setup(nmol_type ),STAT=ios)
     IF ( ios /= 0 ) CALL stop_memory ( 'read_setup_section', &
        'mol_setup', nmol_type )
!    reinitiate the parser at the beginning of the section
     CALL search_label ( label, ierror, ignore_case=.TRUE. )
     IF ( ierr2 /= 0 ) CALL stop_parser ('read_setup_section','reinit failed')
     IF ( .NOT. ASSOCIATED ( mol_setup ) ) &
          CALL stop_program ( 'read_setup_section', &
          'number of molecule types not yet specified')
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
        CASE ( 'ATOM')
           icount = icount + 1
           IF (icount>nmol_type ) CALL stop_program ( 'read_setup_section', &
                'too many molecule types specified')
           ilen = 0
           CALL get_next ( mol_setup(icount) % molpar % name, ilen )
           CALL get_next ( mol_setup(icount) % num_mol )
!          set default for thermostats
           mol_setup(icount) % nhcopt = default_thermostat
           IF ( test_next() == 'C' ) THEN
              ilen = 0
              CALL get_next(mol_setup(icount) % nhcopt,ilen)
              CALL uppercase(mol_setup(icount) % nhcopt)
           ENDIF
        CASE ( 'MOLECULE')
           icount = icount + 1
           IF (icount>nmol_type ) CALL stop_program ( 'read_setup_section', &
                'too many molecule types specified')
!          set default for thermostats
           mol_setup(icount) % nhcopt = default_thermostat
           IF ( test_next() == 'C' ) THEN
              CALL get_next ( mol_setup(icount) % molpar % name, ilen )
              CALL get_next ( mol_setup(icount) % num_mol )
              IF ( test_next() == 'C' ) THEN
                 ilen = 0
                 CALL get_next(mol_setup(icount) % nhcopt,ilen)
                 CALL uppercase(mol_setup(icount) % nhcopt)
              ENDIF
           ELSE
              DO
                 CALL read_line
                 ilen = 6
                 CALL get_next ( string2, ilen )
                 CALL uppercase ( string2 )
                 SELECT CASE ( string2 )
                 CASE DEFAULT
                    CALL p_error()
                    CALL stop_parser ( 'read_setup_section','unknown option')
                 CASE ( 'NAME')
                    ilen = 0
                    CALL get_next(mol_setup(icount) % molpar % name,ilen)
                 CASE ( 'NUMMOL')
                    CALL get_next ( mol_setup(icount) % num_mol )
                 CASE ( 'NHCOPT')
                    ilen = 0
                    CALL get_next(mol_setup(icount) % nhcopt,ilen)
                    CALL uppercase(mol_setup(icount) % nhcopt)
                 CASE ( 'END')
                    ilen = 8
                    CALL get_next(string2,ilen)
                    EXIT
                 END SELECT
              END DO
           END IF
        END SELECT
        CALL read_line
     END DO
  END IF
  CALL parser_end

!..cell setup
  IF ( SUM ( ABS (setup % reference_cell % hmat) ) == 0._dbl ) THEN
    setup % reference_cell = setup % simulation_cell
  ELSE IF ( SUM ( ABS (setup % reference_cell % hmat) ) == 0._dbl ) THEN
    setup % simulation_cell = setup % reference_cell
  END IF
  IF ( SUM ( ABS (setup % reference_cell % hmat) ) == 0._dbl ) THEN
    CALL stop_program ( "read_setup_section", "No cell information available" )
  END IF
!..write some information to output
  IF ( globenv % ionode ) THEN
     IF (globenv%print_level>=0) THEN
        IF ( nmol_type > 0 ) THEN
          WRITE ( iw, '()' )
          WRITE ( iw, '( A,T71,I10 )' ) ' MOLECULE| Number of molecule types ', &
               nmol_type
          WRITE ( iw, '( A,T25,A,T50,A,T69,A )' ) ' MOLECULE| Molecule type', &
               ' Molecule name ', 'Molecules', ' Nose option'
          DO i = 1, nmol_type
             WRITE ( iw, '( A,I6,T19,A,T50,I6,T61,A )' ) ' MOLECULE| ', i, &
                  ADJUSTR ( mol_setup(i) % molpar % name), &
                  mol_setup(i) % num_mol, &
                  ADJUSTR ( mol_setup(i) % nhcopt)
          END DO
        END IF
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
     END IF
  END IF

END SUBROUTINE read_setup_section

!!*****
!******************************************************************************
!!****** cp2k/molecule_input [1.0] *
!!
!!   NAME
!!     molecule_input
!!
!!   FUNCTION
!!     Read input section molecule
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (30.11.2001) : search on input file and set file for molecule section
!!                        short form input
!!                        default settings for non-specified molecules (mass<0)
!!                        determine setup variables (nmol, natoms, natom_type)
!!                        and atom_names (this was in force_fields before)
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &molecule ... &end                                               I
!!I                                                                           I
!!I moleculename name {mass charge}                                           I
!!I                                                                           I
!!I moleculename                                                              I
!!I    numat         n                                                        I
!!I    atom          name label mass charge                                   I
!!I    bonds         nbonds                                                   I
!!I    bends         nbends                                                   I
!!I    torsion       ntorsion                                                 I
!!I    distcon       ndcon                                                    I
!!I    3x3con        n3x3con                                                  I
!!I    basis         name library nmax epsilon                                I
!!I    drho          name library nmax epsilon                                I
!!I end moleculename                                                          I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_molecule_section ( mol_setup, setup, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: ierror, ilen, iw, source, group, icount, nmol_type, ios
  INTEGER :: i, j, imol, na, nat, j1, j2, nf, nt, ie, ia, i1, i2, i3, i4
  INTEGER :: isos, nbasis, ndrho, allgrp, npol, nrho0, anum, n
  CHARACTER ( LEN = 20 ) :: libstring, string, string2, at1, at2, at3, at4
  CHARACTER ( LEN = 6 ) :: label
  CHARACTER ( LEN = 20 ), DIMENSION ( : ), ALLOCATABLE :: molname
  CHARACTER ( LEN = 20 ), DIMENSION ( : ), ALLOCATABLE :: tuple
  CHARACTER ( LEN = 80 ) :: bl = ' '
  CHARACTER ( LEN = 2 ) :: pt,la

!------------------------------------------------------------------------------

  IF ( .NOT. ASSOCIATED ( mol_setup ) ) RETURN

  nmol_type = SIZE ( mol_setup )

  iw = globenv % scr
  ALLOCATE ( molname ( nmol_type ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
       'molname', nmol_type )
  DO i = 1, nmol_type
     molname(i) = mol_setup(i) %molpar%name
  END DO

!..parse the input section
  label = '&MOLEC'
  CALL parser_init ( setup % set_file_name, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
     IF ( globenv % ionode .AND. globenv % print_level >= 2 ) THEN
        WRITE ( iw, '( A )' ) ' No input section &MOLECULE found on file '
        WRITE ( iw, '( T2, A )' ) setup % set_file_name
     END IF
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 8
        CALL get_next(string,ilen)
        imol = str_search(molname,nmol_type,string)
        nat = -1
        IF (imol>0) THEN
!          we have found a valid molecule description
           mol_setup(imol) %molpar%natom_type = 0
           mol_setup(imol) %molpar%nbonds = 0
           mol_setup(imol) %molpar%nbends = 0
           mol_setup(imol) %molpar%ntorsions = 0
           mol_setup(imol) %molpar%ndcon = 0
           mol_setup(imol) %molpar%n3x3con = 0
           mol_setup(imol) %molpar%nbasis_pol = 0
           mol_setup(imol) %molpar%nbasis_rho0 = 0

           IF (test_next()=="C") THEN
!             short form, single atom molecules
              mol_setup(imol) %molpar%natom_type = 1
              mol_setup(imol) %molpar%natom = 1
              ALLOCATE (mol_setup(imol)%molpar%aname(1),STAT=ios)
              IF ( ios /= 0 ) CALL stop_memory &
                   ( 'read_molecule_section', 'aname', 1 )
              ALLOCATE (mol_setup(imol)%molpar%alabel(1),STAT=ios)
              IF ( ios /= 0 ) CALL stop_memory &
                   ( 'read_molecule_section', 'alabel', 1 )
              ALLOCATE (mol_setup(imol)%molpar%aweight(1),STAT=ios)
              IF ( ios /= 0 ) CALL stop_memory &
                   ( 'read_molecule_section', 'aweight', 1 )
              ALLOCATE (mol_setup(imol)%molpar%acharge(1),STAT=ios)
              IF ( ios /= 0 ) CALL stop_memory &
                   ( 'read_molecule_section', 'acharge', 1 )
              ilen = 0
              CALL get_next(mol_setup(imol)%molpar%aname(1),ilen)
              la = mol_setup(imol)%molpar%aname(1)
              CALL uppercase ( la )
              DO i = 0, SIZE ( ptable )
                pt = ptable ( i ) % symbol
                CALL uppercase ( pt )
                IF ( pt == la ) THEN
                  mol_setup(imol)%molpar%alabel(1) = ptable ( i ) % symbol
                  anum = i
                  EXIT
                ENDIF
              END DO
              mol_setup(imol) %molpar%acharge(1) = 0._dbl
              IF (test_next()=="N") THEN
                 CALL get_next ( mol_setup(imol) %molpar%aweight(1) )
                 IF (test_next()=="N") THEN
                    CALL get_next ( mol_setup(imol) %molpar%acharge(1) )
                 END IF
              ELSE
                 mol_setup(imol) % molpar%aweight(1) = ptable ( anum ) % amass
              END IF
           ELSE
!             extended form
              DO
                 CALL read_line
                 ilen = 8
                 CALL get_next(string2,ilen)
                 CALL uppercase ( string2 )
   
                 SELECT CASE (string2)
                 CASE DEFAULT
                    CALL p_error()
                    CALL stop_parser ( 'read_molecule_section','unknown option')
   
                 CASE ( 'NUMAT')
                    CALL get_next ( na )
                    mol_setup(imol) %molpar%natom = na
                    ALLOCATE (mol_setup(imol)%molpar%aname(na),STAT=ios)
                    IF ( ios /= 0 ) CALL stop_memory &
                         ( 'read_molecule_section', 'aname', na )
                    ALLOCATE (mol_setup(imol)%molpar%alabel(na),STAT=ios)
                    IF ( ios /= 0 ) CALL stop_memory &
                         ( 'read_molecule_section', 'alabel', na )
                    ALLOCATE (mol_setup(imol)%molpar%aweight(na),STAT=ios)
                    IF ( ios /= 0 ) CALL stop_memory &
                         ( 'read_molecule_section', 'aweight', na )
                    ALLOCATE (mol_setup(imol)%molpar%acharge(na),STAT=ios)
                    IF ( ios /= 0 ) CALL stop_memory &
                         ( 'read_molecule_section', 'acharge', na )
                    nat = 0
   
                 CASE ( 'ATOM' )
                    IF ( nat < 0 ) CALL stop_parser ( 'read_molecule_section', &
                      'aweight' ) !TMPTMPTMP
                    nat = nat + 1
                    ilen = 0
                    CALL get_next(mol_setup(imol)%molpar%aname(nat),ilen)
                    ilen = 0
                    CALL get_next(mol_setup(imol)%molpar%alabel(nat),ilen)
                    CALL get_next ( mol_setup(imol) %molpar%aweight(nat) )
                    CALL get_next ( mol_setup(imol) %molpar%acharge(nat) )
   
                 CASE ( 'BONDS')
                    CALL get_next ( mol_setup(imol) %molpar%nbonds )
                    na = mol_setup(imol) %molpar%nbonds
                    IF (na>=0) THEN
                       ALLOCATE (mol_setup(imol)%molpar%bonds(2,na),STAT=ios)
                       IF (ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'bonds', 2 * na )
                       DO i = 1, na
                          CALL get_next ( mol_setup(imol) %molpar%bonds(1,i) )
                          CALL get_next ( mol_setup(imol) %molpar%bonds(2,i) )
                       END DO
                    END IF
   
                 CASE ( 'BENDS')
                    CALL get_next ( mol_setup(imol) %molpar%nbends )
                    na = mol_setup(imol) %molpar%nbends
                    IF (na>=0) THEN
                       ALLOCATE (mol_setup(imol)%molpar%bends(3,na),STAT=ios)
                       IF (ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'bends', 3 * na )
                       DO i = 1, na
                          CALL get_next ( mol_setup(imol) %molpar%bends(1,i) )
                          CALL get_next ( mol_setup(imol) %molpar%bends(2,i) )
                          CALL get_next ( mol_setup(imol) %molpar%bends(3,i) )
                       END DO
                    END IF
   
                 CASE ( 'TORSIO')
                    CALL get_next ( mol_setup(imol) %molpar%ntorsions )
                    na = mol_setup(imol) %molpar%ntorsions
                    IF (na>=0) THEN
                       ALLOCATE ( &
                            mol_setup ( imol ) % molpar % torsions ( 4, na ), &
                            STAT = ios )
                       IF ( ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'torsions', 4 * na )
                       DO i = 1, na
                          CALL get_next ( mol_setup(imol) %molpar%torsions(1,i) )
                          CALL get_next ( mol_setup(imol) %molpar%torsions(2,i) )
                          CALL get_next ( mol_setup(imol) %molpar%torsions(3,i) )
                          CALL get_next ( mol_setup(imol) %molpar%torsions(4,i) )
                       END DO
                    END IF
   
                 CASE ( 'DISTCO')
                    CALL get_next ( mol_setup(imol) %molpar%ndcon )
                    na = mol_setup(imol) %molpar%ndcon
                    IF (na>=0) THEN
                       ALLOCATE (mol_setup(imol)%molpar%ndc(2,na),STAT=ios)
                       IF ( ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'ndc', 2 * na )
                       DO i = 1, na
                          CALL get_next ( mol_setup(imol) %molpar%ndc(1,i) )
                          CALL get_next ( mol_setup(imol) %molpar%ndc(2,i) )
                       END DO
                    END IF
   
                 CASE ( '3X3CON')
                    CALL get_next ( mol_setup(imol) %molpar%n3x3con )
                    na = mol_setup(imol) %molpar%n3x3con
                    IF (na>=0) THEN
                       ALLOCATE (mol_setup(imol)%molpar%n3x3c(3,na),STAT=ios)
                       IF (ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'n3x3c', 3 * na )
                       DO i = 1, na
                          CALL get_next ( mol_setup(imol) %molpar%n3x3c(1,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n3x3c(2,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n3x3c(3,i) )
                       END DO
                    END IF
    
                 CASE ('RHO0_SET')
        
                    na = mol_setup(imol) %molpar%natom  
                    IF (.NOT. ASSOCIATED (mol_setup(imol)%molpar%rho0_basis)) THEN
                      ALLOCATE (mol_setup(imol)%molpar%rho0_basis(na),STAT=ios)
                      IF ( ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'rho0_basis',  na )
                    END IF
                    nrho0 = mol_setup(imol)%molpar%nbasis_rho0  + 1     
                    mol_setup(imol)%molpar%nbasis_rho0 = nrho0
                    ilen=0
                    CALL get_next(mol_setup(imol)%molpar%rho0_basis(nrho0)%aname,ilen)
                    ilen=0
                    CALL get_next(mol_setup(imol)%molpar%rho0_basis(nrho0)%libname,ilen)
                    CALL get_next ( mol_setup(imol)%molpar%rho0_basis(nrho0)%nmax )
                    CALL get_next ( mol_setup(imol)%molpar%rho0_basis(nrho0)%eps )
   
                 CASE ('POL_SET')
                    na = mol_setup(imol) %molpar%natom  
                    IF (.NOT. ASSOCIATED (mol_setup(imol)%molpar%pol_basis)) THEN
                       ALLOCATE (mol_setup(imol)%molpar%pol_basis(na),STAT=ios)
                       IF ( ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'pol_basis',  na )
                    END IF
                    npol = mol_setup(imol)%molpar%nbasis_pol + 1     
                    mol_setup(imol)%molpar%nbasis_pol = npol
                    ilen=0
                    CALL get_next(mol_setup(imol)%molpar%pol_basis(npol)%aname,ilen)
                    ilen=0
                    CALL get_next(mol_setup(imol)%molpar%pol_basis(npol)%libname,ilen)
                    CALL get_next ( mol_setup(imol)%molpar%pol_basis(npol)%nmax )
                    CALL get_next ( mol_setup(imol)%molpar%pol_basis(npol)%eps )
   
                 CASE ( 'END')
                    ilen = 0
                    CALL get_next(string2,ilen)
                    EXIT
                 END SELECT
   
              END DO

!..count number of atom types
              mol_setup(imol) %molpar%natom_type = 1
              DO i = 2, mol_setup(imol) %molpar%natom
                 j = str_search ( &
                      mol_setup ( imol ) % molpar % aname ( 1:i-1 ), i-1, &
                      mol_setup ( imol ) % molpar % aname ( i ) )
                 IF ( j == 0 ) mol_setup ( imol ) % molpar % natom_type &
                      = mol_setup ( imol ) % molpar % natom_type + 1
              END DO
           END IF
        ELSE
           WRITE ( iw, '( A )' ) ' molecule type not found '
           WRITE ( iw, '( T2, A )' ) string
           CALL stop_parser ( 'read_molecule_section', ' ' )
        END IF
        CALL read_line
     END DO
  END IF
  CALL parser_end

! add some defaults to the non specified molecules
  DO imol = 1, nmol_type
     IF ( mol_setup(imol) %molpar%natom == 0 ) THEN
        mol_setup(imol) %molpar%natom_type = 1
        mol_setup(imol) %molpar%natom = 1
        ALLOCATE (mol_setup(imol)%molpar%aname(1),STAT=ios)
        IF ( ios /= 0 ) CALL stop_memory &
             ( 'read_molecule_section', 'aname', 1 )
        ALLOCATE (mol_setup(imol)%molpar%alabel(1),STAT=ios)
        IF ( ios /= 0 ) CALL stop_memory &
             ( 'read_molecule_section', 'alabel', 1 )
        ALLOCATE (mol_setup(imol)%molpar%aweight(1),STAT=ios)
        IF ( ios /= 0 ) CALL stop_memory &
             ( 'read_molecule_section', 'aweight', 1 )
        ALLOCATE (mol_setup(imol)%molpar%acharge(1),STAT=ios)
        IF ( ios /= 0 ) CALL stop_memory &
             ( 'read_molecule_section', 'acharge', 1 )
        mol_setup(imol) %molpar%aname(1) = molname(imol)
        mol_setup(imol) %molpar%alabel(1) = ""
        mol_setup(imol) %molpar%aweight(1) = -1._dbl
        mol_setup(imol) %molpar%acharge(1) = 0._dbl
        mol_setup(imol) %molpar%nbonds = 0
        mol_setup(imol) %molpar%nbends = 0
        mol_setup(imol) %molpar%ntorsions = 0
        mol_setup(imol) %molpar%ndcon = 0
        mol_setup(imol) %molpar%n3x3con = 0
        mol_setup(imol) %molpar%nbasis_pol = 0
        mol_setup(imol) %molpar%nbasis_rho0 = 0
     END IF
  END DO

  DEALLOCATE (molname,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', 'molname')

!..count total number of molecules and atoms
  setup % nmol = 0
  setup % natoms = 0
  setup % natom_type = 0
  DO i = 1, nmol_type
     setup % nmol = setup % nmol + mol_setup(i) % num_mol
     setup % natoms = setup % natoms + mol_setup(i) % num_mol *  &
          mol_setup(i) % molpar % natom
     setup % natom_type = setup % natom_type +  &
          mol_setup(i) % molpar % natom_type
  END DO

!..generate unique atom list
  ALLOCATE ( setup % atom_names ( setup % natom_type ), STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'read_molecule_section', &
       'setup % atom_names', setup % natom_type )

  n = 0
  DO i = 1, nmol_type
     DO j = 1, mol_setup(i) % molpar % natom
        IF (str_search(setup % atom_names,n,mol_setup(i) % molpar%aname(j))==0) THEN
           n = n + 1
           setup % atom_names(n) = mol_setup(i) % molpar%aname(j)
        END IF
     END DO
  END DO

  IF ( n /= setup % natom_type ) CALL stop_program ( 'read_molecule_section', &
       'natom_type mismatch' )

! write some information to output
  IF (globenv%ionode) THEN
     IF (globenv%print_level>=0) THEN
        DO imol = 1, nmol_type
           IF ( mol_setup(imol) %molpar%aweight(1) < 0._dbl ) CYCLE
           WRITE ( iw,'( )' )
           WRITE ( iw, '( A,T61,A )' ) ' MOLECULE| Molecule name ', &
                ADJUSTR ( mol_setup(imol)%molpar%name)
           WRITE ( iw, '( A,T71,I10 )' ) ' MOLECULE| Number of atoms ', &
                mol_setup(imol) %molpar%natom
           WRITE ( iw, '( A,T71,I10 )' ) ' MOLECULE| Number of atom types ', &
                mol_setup(imol) %molpar%natom_type
           WRITE ( iw, '( A,T40,A,T55,A )' ) ' MOLECULE|', 'Atom', &
                '       Weight       Charge'
           DO i = 1, mol_setup(imol) %molpar%natom
              CALL xstring(mol_setup(imol)%molpar%aname(i),ia,ie)
              WRITE ( iw, '( A,T40,A,T55,F13.5,F13.5 )' ) ' MOLECULE| ', &
                   mol_setup(imol) %molpar%aname(i) (ia:ie), &
                   mol_setup(imol) %molpar%aweight(i), &
                   mol_setup(imol) %molpar%acharge(i)
           END DO

           na = mol_setup(imol) %molpar%nbonds
           IF (na>0) THEN
              WRITE ( iw, '( A,T71,I10 )' ) ' MOLECULE| Number of bonds ', &
                   mol_setup(imol) %molpar%nbonds
              ALLOCATE ( tuple ( na ), STAT = isos )
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple', na )
              CALL make_tuple(mol_setup(imol)%molpar%bonds,2,na,tuple)
              CALL xstring(tuple(1),ia,ie)
              ilen = ie - ia + 1
              nt = 69/ilen
              DO i = 1, na, nt
                 j1 = i
                 j2 = min(j1+nt-1,na)
                 nf = 69 - ilen*(j2-j1+1)
                 WRITE ( iw, '( 80(A) )' ) ' MOLECULE| ', bl(1:nf), &
                      (tuple(j)(ia:ie),j=j1,j2)
              END DO
              DEALLOCATE (tuple,STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple' )
           END IF

           na = mol_setup(imol) %molpar%nbends
           IF (na>0) THEN
              WRITE ( iw, '( A,T71,I10 )' ) ' MOLECULE| Number of bends ', &
                   mol_setup(imol) %molpar%nbends
              ALLOCATE (tuple(na),STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple', na )
              CALL make_tuple(mol_setup(imol)%molpar%bends,3,na,tuple)
              CALL xstring(tuple(1),ia,ie)
              ilen = ie - ia + 1
              nt = 69/ilen
              DO i = 1, na, nt
                 j1 = i
                 j2 = min(j1+nt-1,na)
                 nf = 69 - ilen*(j2-j1+1)
                 WRITE ( iw, '( 80(A) )' ) ' MOLECULE| ', bl(1:nf), &
                      (tuple(j)(ia:ie),j=j1,j2)
              END DO
              DEALLOCATE (tuple,STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple')
           END IF

           na = mol_setup(imol) %molpar%ntorsions
           IF (na>0) THEN
              WRITE ( iw, '( A,T71,I10 )' ) ' MOLECULE| Number of torsions ', &
                   mol_setup(imol) %molpar%ntorsions
              ALLOCATE (tuple(na),STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple', na )
              CALL make_tuple(mol_setup(imol)%molpar%torsions,4,na,tuple)
              CALL xstring(tuple(1),ia,ie)
              ilen = ie - ia + 1
              nt = 69/ilen
              DO i = 1, na, nt
                 j1 = i
                 j2 = min(j1+nt-1,na)
                 nf = 69 - ilen*(j2-j1+1)
                 WRITE ( iw, '( 80(A) )' ) ' MOLECULE| ', bl(1:nf), &
                      (tuple(j)(ia:ie),j=j1,j2)
              END DO
              DEALLOCATE (tuple,STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple')
           END IF

           na = mol_setup(imol) %molpar%ndcon
           IF (na>0) THEN
              WRITE ( iw, '( A,T71,I10 )' ) &
                   ' MOLECULE| Number of distance constraints', &
                   mol_setup(imol) %molpar%ndcon
              ALLOCATE (tuple(na),STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple', na )
              CALL make_tuple(mol_setup(imol)%molpar%ndc,2,na,tuple)
              CALL xstring(tuple(1),ia,ie)
              ilen = ie - ia + 1
              nt = 69/ilen
              DO i = 1, na, nt
                 j1 = i
                 j2 = min(j1+nt-1,na)
                 nf = 69 - ilen*(j2-j1+1)
                 WRITE ( iw, '( 80(A) )' ) ' MOLECULE| ', bl(1:nf), &
                      (tuple(j)(ia:ie),j=j1,j2)
              END DO
              DEALLOCATE (tuple,STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple')
           END IF

           na = mol_setup(imol) %molpar%n3x3con
           IF (na>0) THEN
              WRITE ( iw, '( A,T71,I10 )' ) &
                   ' MOLECULE| Number of 3x3 constraints', &
                   mol_setup(imol) %molpar%n3x3con
              ALLOCATE (tuple(na),STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple', na )
              CALL make_tuple(mol_setup(imol)%molpar%n3x3c,3,na,tuple)
              CALL xstring(tuple(1),ia,ie)
              ilen = ie - ia + 1
              nt = 69/ilen
              DO i = 1, na, nt
                 j1 = i
                 j2 = min(j1+nt-1,na)
                 nf = 69 - ilen*(j2-j1+1)
                 WRITE ( iw, '( 80(A) )' ) ' MOLECULE| ', bl(1:nf), &
                      (tuple(j)(ia:ie),j=j1,j2)
              END DO
              DEALLOCATE (tuple,STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple' )
           END IF
   
           na = mol_setup(imol) %molpar%nbasis_rho0
           IF (na>0) THEN
              WRITE ( iw, '( A,T71,I10 )' ) &
                   ' MOLECULE| Number of effective densities', &
                   mol_setup(imol) %molpar%nbasis_rho0
               DO i = 1, mol_setup(imol)%molpar%nbasis_rho0
                  WRITE (iw,'(A,T61,A)') ' MOLECULE| name of rho0 basis set ', &
                   mol_setup(imol)%molpar%rho0_basis(i)%aname
                  WRITE (iw,'(A,T61,A)') ' MOLECULE| library for basis set ', &
                   mol_setup(imol)%molpar%rho0_basis(i)%libname   
                  WRITE (iw,'(A,T71,I10)') ' MOLECULE| max grid points ', &
                   mol_setup(imol)%molpar%rho0_basis(i)%nmax  
                  WRITE (iw,'(A,T68,F12.6)') ' MOLECULE| epsilon ', &
                   mol_setup(imol)%molpar%rho0_basis(i)%eps   
               END DO
           END IF

           na = mol_setup(imol) %molpar%nbasis_pol
           IF (na>0) THEN
              WRITE ( iw, '( A,T71,I10 )' ) &
                   ' MOLECULE| Number of polarization densities', &
                   mol_setup(imol) %molpar%nbasis_pol
              DO i = 1, mol_setup(imol)%molpar%nbasis_pol
                  WRITE (iw,'(A,T61,A)') ' MOLECULE| name of polarizable basis set ', &
                   mol_setup(imol)%molpar%pol_basis(i)%aname
                  WRITE (iw,'(A,T61,A)') ' MOLECULE| library for basis set ', &
                   mol_setup(imol)%molpar%pol_basis(i)%libname   
                  WRITE (iw,'(A,T71,I10)') ' MOLECULE| max grid points ', &
                   mol_setup(imol)%molpar%pol_basis(i)%nmax  
                  WRITE (iw,'(A,T68,F12.6)') ' MOLECULE| epsilon ', &
                   mol_setup(imol)%molpar%pol_basis(i)%eps   
              END DO
           END IF

        END DO
     END IF
  END IF

END SUBROUTINE read_molecule_section

!!*****
!******************************************************************************

END MODULE molecule_input

!******************************************************************************
