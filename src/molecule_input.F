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
!!     JGH (26-01-2002) Added read_topology_section
!!
!!   SOURCE
!******************************************************************************

MODULE molecule_input
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE global_types,                    ONLY: global_environment_type
  USE input_types,                     ONLY: setup_parameters_type
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: pi,&
                                             zero
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_bcast
  USE molecule_types,                  ONLY: init_molecule_type,&
                                             molecule_type,&
                                             topology_type
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
  USE topology_control,                ONLY: read_topology_control
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_setup_section, read_molecule_section, read_topology_section

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** molecule_input/read_setup_section [1.0] *
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
  INTEGER :: ierror, i, j, nmol_topo
  LOGICAL :: sim, section_found
  CHARACTER ( LEN = 20 ) :: string, string2
  CHARACTER ( LEN = 6 ) :: label
  CHARACTER ( LEN = 12 ) :: default_thermostat
  CHARACTER ( LEN = 3 ), PARAMETER :: yn ( 0:1 ) = (/ ' NO', 'YES' /)

  default_thermostat = "GLOBAL"
  iw = globenv % scr
  icount = 0
  nmol_type = -1
  nmol_topo = 0
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

  IF(setup%topology) THEN
    CALL read_topology_control(mol_setup,setup,globenv)
    nmol_topo = SIZE(mol_setup)
  END IF

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
!    molecule types defined on topology file
     IF ( nmol_topo > 0 ) THEN
        IF ( nmol_topo < nmol_type ) &
          CALL stop_program ( 'read_setup_section', &
                ' inconsistent number of molecules specified ')
        nmol_type = nmol_topo
     END IF
!    allocate space for the molecules
     IF ( nmol_topo == 0 ) THEN
        ALLOCATE (mol_setup(nmol_type ),STAT=ios)
        IF ( ios /= 0 ) CALL stop_memory ( 'read_setup_section', &
           'mol_setup', nmol_type )
        DO i=1,SIZE(mol_setup)
           CALL init_molecule_type(mol_setup(i))
        ENDDO
     END IF
  END IF
  CALL parser_end

!!-------------------------------------------------------------------------
!!Commented out....
  !mol_setup (:) % molpar % name = "#"
!!-------------------------------------------------------------------------

! process the information from the topology file
  !CALL set_molecule_type_topology ( mol_setup, default_thermostat )

  IF ( section_found ) THEN
! reinitiate the parser at the beginning of the section
     CALL parser_init ( setup % set_file_name, globenv )
     CALL search_label ( label, ierror, ignore_case=.TRUE. )
     IF ( ierror /= 0 ) CALL stop_parser ('read_setup_section','reinit failed')
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

     CALL parser_end
  END IF

! cell setup
  IF ( SUM ( ABS (setup % reference_cell % hmat) ) == 0._dbl ) THEN
    setup % reference_cell = setup % simulation_cell
  ELSE IF ( SUM ( ABS (setup % reference_cell % hmat) ) == 0._dbl ) THEN
    setup % simulation_cell = setup % reference_cell
  END IF
  IF ( SUM ( ABS (setup % reference_cell % hmat) ) == 0._dbl ) THEN
    CALL stop_program ( "read_setup_section", "No cell information available" )
  END IF

! write some information to output
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
        WRITE ( iw, '()' )
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
!!I    4x6con        n4x6con                                                  I
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
  INTEGER :: isos, allgrp, npol, nrho0, anum, n
  CHARACTER ( LEN = 20 ) :: libstring, string, string2, at1, at2, at3, at4
  CHARACTER ( LEN = 6 ) :: label
  CHARACTER ( LEN = 20 ), DIMENSION ( : ), ALLOCATABLE :: molname
  CHARACTER ( LEN = 20 ), DIMENSION ( : ), ALLOCATABLE :: tuple
  CHARACTER ( LEN = 80 ) :: bl = ' '
  CHARACTER ( LEN = 2 ) :: pt,la
  CHARACTER ( LEN = 20 ), DIMENSION(:), ALLOCATABLE :: work

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

  IF(setup%topology) THEN
    ! get molecule structures from topology file
    !CALL get_molecule_structure_topology ( molname, mol_setup )
  ELSE

  ! parse the input section
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
           ! we have found a valid molecule description
           ! it would be better to first deallocate all structures already
           ! allocated at this point
           mol_setup(imol) %molpar%natom_type = 0
           mol_setup(imol) %molpar%nbonds = 0
           mol_setup(imol) %molpar%nbends = 0
           mol_setup(imol) %molpar%ntorsions = 0
           mol_setup(imol) %molpar%ndcon = 0
           mol_setup(imol) %molpar%n3x3con = 0
           mol_setup(imol) %molpar%n4x6con = 0

           IF (test_next()=="C") THEN
              ! short form, single atom molecules
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
                       ALLOCATE (mol_setup(imol)%molpar%bonds(3,na),STAT=ios)
                       IF (ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'bonds', 3 * na )
                       DO i = 1, na
                          CALL get_next ( mol_setup(imol) %molpar%bonds(1,i) )
                          CALL get_next ( mol_setup(imol) %molpar%bonds(2,i) )
                          CALL get_next ( mol_setup(imol) %molpar%bonds(3,i) )
                       END DO
                    END IF
   
                 CASE ( 'BENDS')
                    CALL get_next ( mol_setup(imol) %molpar%nbends )
                    na = mol_setup(imol) %molpar%nbends
                    IF (na>=0) THEN
                       ALLOCATE (mol_setup(imol)%molpar%bends(4,na),STAT=ios)
                       IF (ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'bends', 4 * na )
                       DO i = 1, na
                          CALL get_next ( mol_setup(imol) %molpar%bends(1,i) )
                          CALL get_next ( mol_setup(imol) %molpar%bends(2,i) )
                          CALL get_next ( mol_setup(imol) %molpar%bends(3,i) )
                          CALL get_next ( mol_setup(imol) %molpar%bends(4,i) )
                       END DO
                    END IF
   
                 CASE ( 'TORSIO')
                    CALL get_next ( mol_setup(imol) %molpar%ntorsions )
                    na = mol_setup(imol) %molpar%ntorsions
                    IF (na>=0) THEN
                       ALLOCATE ( &
                            mol_setup ( imol ) % molpar % torsions ( 5, na ), &
                            STAT = ios )
                       IF ( ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'torsions', 5 * na )
                       DO i = 1, na
                          CALL get_next ( mol_setup(imol) %molpar%torsions(1,i) )
                          CALL get_next ( mol_setup(imol) %molpar%torsions(2,i) )
                          CALL get_next ( mol_setup(imol) %molpar%torsions(3,i) )
                          CALL get_next ( mol_setup(imol) %molpar%torsions(4,i) )
                          CALL get_next ( mol_setup(imol) %molpar%torsions(5,i) )
                       END DO
                    END IF
   
                 CASE ( 'DISTCO')
                    CALL get_next ( mol_setup(imol) %molpar%ndcon )
                    na = mol_setup(imol) %molpar%ndcon
                    IF (na>=0) THEN
                       ALLOCATE (mol_setup(imol)%molpar%ndc(3,na),STAT=ios)
                       IF ( ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'ndc', 3 * na )
                       DO i = 1, na
                          CALL get_next ( mol_setup(imol) %molpar%ndc(1,i) )
                          CALL get_next ( mol_setup(imol) %molpar%ndc(2,i) )
                          CALL get_next ( mol_setup(imol) %molpar%ndc(3,i) )
                       END DO
                    END IF
   
                 CASE ( '3X3CON')
                    CALL get_next ( mol_setup(imol) %molpar%n3x3con )
                    na = mol_setup(imol) %molpar%n3x3con
                    IF (na>=0) THEN
                       ALLOCATE (mol_setup(imol)%molpar%n3x3c(6,na),STAT=ios)
                       IF (ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'n3x3c', 6 * na )
                       DO i = 1, na
                          CALL get_next ( mol_setup(imol) %molpar%n3x3c(1,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n3x3c(2,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n3x3c(3,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n3x3c(4,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n3x3c(5,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n3x3c(6,i) )
                       END DO
                    END IF

                 CASE ( '4X6CON')
                    CALL get_next ( mol_setup(imol) %molpar%n4x6con )
                    na = mol_setup(imol) %molpar%n4x6con
                    IF (na>=0) THEN
                       ALLOCATE (mol_setup(imol)%molpar%n4x6c(10,na),STAT=ios)
                       IF (ios /= 0 ) CALL stop_memory &
                            ( 'read_molecule_section', 'n4x6c', 10 * na )
                       DO i = 1, na
                          CALL get_next ( mol_setup(imol) %molpar%n4x6c(1,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n4x6c(2,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n4x6c(3,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n4x6c(4,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n4x6c(5,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n4x6c(6,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n4x6c(7,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n4x6c(8,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n4x6c(9,i) )
                          CALL get_next ( mol_setup(imol) %molpar%n4x6c(10,i) )
                       END DO
                    END IF
    
                 CASE ( 'END')
                    ilen = 0
                    CALL get_next(string2,ilen)
                    EXIT
                 END SELECT
   
              END DO

!             count number of atom types
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
  END IF

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
        mol_setup(imol) %molpar%n4x6con = 0
     END IF
  END DO

  DEALLOCATE (molname,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', 'molname')

! count total number of molecules and atoms
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

! generate unique atom list
  ALLOCATE ( work( setup % natom_type ), STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'read_molecule_section', &
       'work', setup % natom_type )

  n = 0
  DO i = 1, nmol_type
     DO j = 1, mol_setup(i) % molpar % natom
        IF(str_search(work,n,mol_setup(i)%molpar%aname(j))==0)THEN
           n = n + 1
           work(n) = mol_setup(i) % molpar%aname(j)
        END IF
     END DO
  END DO

  ALLOCATE ( setup%atom_names(n), STAT=ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'read_molecule_section', &
       'setup%atom_names', n )
  setup % natom_type = n
  setup % atom_names(1:n) = work(1:n)
  DEALLOCATE( work, STAT=ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'read_molecule_section', &
       'work' )
  
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
                 j2 = MIN(j1+nt-1,na)
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
                 j2 = MIN(j1+nt-1,na)
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
                 j2 = MIN(j1+nt-1,na)
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
                 j2 = MIN(j1+nt-1,na)
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
                 j2 = MIN(j1+nt-1,na)
                 nf = 69 - ilen*(j2-j1+1)
                 WRITE ( iw, '( 80(A) )' ) ' MOLECULE| ', bl(1:nf), &
                      (tuple(j)(ia:ie),j=j1,j2)
              END DO
              DEALLOCATE (tuple,STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple' )
           END IF

           na = mol_setup(imol) %molpar%n4x6con
           IF (na>0) THEN
              WRITE ( iw, '( A,T71,I10 )' ) &
                   ' MOLECULE| Number of 4x6 constraints', &
                   mol_setup(imol) %molpar%n4x6con
              ALLOCATE (tuple(na),STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple', na )
              CALL make_tuple(mol_setup(imol)%molpar%n4x6c,4,na,tuple)
              CALL xstring(tuple(1),ia,ie)
              ilen = ie - ia + 1
              nt = 69/ilen
              DO i = 1, na, nt
                 j1 = i
                 j2 = MIN(j1+nt-1,na)
                 nf = 69 - ilen*(j2-j1+1)
                 WRITE ( iw, '( 80(A) )' ) ' MOLECULE| ', bl(1:nf), &
                      (tuple(j)(ia:ie),j=j1,j2)
              END DO
              DEALLOCATE (tuple,STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', &
                   'tuple' )
           END IF
   
        END DO
     END IF
  END IF
! convert amu -> au
  DO imol = 1, nmol_type
    DO i = 1, SIZE ( mol_setup ( imol ) % molpar % aweight )
       CALL convert_to_cp2k_units ( mass = mol_setup ( imol ) % molpar % aweight ( i ) )
    END DO
  END DO 

END SUBROUTINE read_molecule_section

!!*****
!******************************************************************************
!!****** molecule_input/read_topology_section [1.0] *
!!
!!   NAME
!!     read_topology_section
!!
!!   SYNOPSIS
!!     Subroutine read_topology_section(topo, setup, globenv)
!!       Implicit None
!!       Type(topology_type):: topo
!!       Type(setup_parameters_type):: setup
!!       Type(global_environment_type), Intent (IN):: globenv
!!     End Subroutine read_topology_section
!!
!!   FUNCTION
!!     reads the input section topology
!!
!!   AUTHOR
!!     JGH (26-01-2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &topology ... &end                                               I
!!I                                                                           I
!!I    bonds         nbonds                                                   I
!!I     iat jat nb                                                            I
!!I    bends         nbends                                                   I
!!I     iat jat kat nb                                                        I
!!I    torsion       ntors                                                    I
!!I     iat jat kat lat nt                                                    I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!******************************************************************************
SUBROUTINE read_topology_section ( topo, setup, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( topology_type ) :: topo
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: ierror,i,n,ilen,iw
  CHARACTER ( LEN = 20 ) :: string
  CHARACTER ( LEN = 6 ) :: label

!------------------------------------------------------------------------------

  iw = globenv % scr

  topo%nbond = 0
  topo%nbend = 0
  topo%ntors = 0

  IF(setup%topology) THEN
    !CALL get_topology_from_file ( )
  ELSE

  ! parse the input section
  label = '&TOPOL'
  CALL parser_init ( setup % set_file_name, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
     IF ( globenv % ionode .AND. globenv % print_level >= 2 ) THEN
        WRITE ( iw, '( A )' ) ' No input section &TOPOLOGY found on file '
        WRITE ( iw, '( T2, A )' ) setup % set_file_name
     END IF
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 0
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        SELECT CASE ( string )
        CASE ( 'BOND')
          CALL get_next ( n )
          topo%nbond = n
          ALLOCATE ( topo%bond(3,n), STAT=ierror )
          IF ( ierror /= 0 ) CALL stop_memory ( "read_topology_section", &
            "topo%bond", 3*n )
          CALL read_line
          DO i = 1, n
            CALL get_next ( topo%bond(1,i) )
            CALL get_next ( topo%bond(2,i) )
            CALL get_next ( topo%bond(3,i) )
          END DO
        CASE ( 'BEND')
          CALL get_next ( n )
          topo%nbend = n
          ALLOCATE ( topo%bend(4,n), STAT=ierror )
          IF ( ierror /= 0 ) CALL stop_memory ( "read_topology_section", &
            "topo%bend", 4*n )
          CALL read_line
          DO i = 1, n
            CALL get_next ( topo%bend(1,i) )
            CALL get_next ( topo%bend(2,i) )
            CALL get_next ( topo%bend(3,i) )
            CALL get_next ( topo%bend(4,i) )
          END DO
        CASE ( 'TORSION')
          CALL get_next ( n )
          topo%ntors = n
          ALLOCATE ( topo%torsion(5,n), STAT=ierror )
          IF ( ierror /= 0 ) CALL stop_memory ( "read_topology_section", &
            "topo%torsion", 5*n )
          CALL read_line
          DO i = 1, n
            CALL get_next ( topo%torsion(1,i) )
            CALL get_next ( topo%torsion(2,i) )
            CALL get_next ( topo%torsion(3,i) )
            CALL get_next ( topo%torsion(4,i) )
            CALL get_next ( topo%torsion(5,i) )
          END DO
        CASE DEFAULT
           CALL stop_program ( "read_topology_section", &
               "Unknown option "//string )
        END SELECT
     END DO
  END IF

  CALL parser_end
  END IF

! write some information to output
  IF (globenv%ionode) THEN
     IF (globenv%print_level>=0) THEN
       IF (topo%nbond > 0) THEN
         WRITE ( iw,'( )' )
         DO i=1,topo%nbond
           WRITE ( iw,'( A,T40,A,I7,T55,A,I7,T70,A,I3 )' ) &
                 "TOPOLOGY| Bond","atom",topo%bond(1,i),"atom",&
                 topo%bond(2,i),"bondtype",topo%bond(3,i)
         END DO
       ENDIF
       IF (topo%nbend > 0) THEN
         WRITE ( iw,'( )' )
         DO i=1,topo%nbond
           WRITE ( iw,'(A,T25,A,I7,T40,A,I7,T55,A,I7,T70,A,I3 )' ) &
                 "TOPOLOGY| Bend","atom",topo%bend(1,i),"atom",&
                 topo%bend(2,i),"atom",topo%bend(3,i),"bendtype",topo%bend(4,i)
         END DO
       ENDIF
       IF (topo%ntors > 0) THEN
         WRITE ( iw,'( )' )
         DO i=1,topo%ntors
           WRITE ( iw,'(A,T20,A,I6,T33,A,I6,T46,A,I6,T59,A,i6,T72,A,I5 )' ) &
             "TOPOLOGY| Torsion","atom",topo%torsion(1,i),"atom",&
             topo%torsion(2,i),"atom",topo%torsion(3,i),"atom",&
             topo%torsion(4,i),"type",topo%torsion(5,i)
         END DO
       ENDIF
     ENDIF
  ENDIF

END SUBROUTINE read_topology_section

!!*****
!******************************************************************************

END MODULE molecule_input

!******************************************************************************
