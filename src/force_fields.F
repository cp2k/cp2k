!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/force_fields [1.0] *
!!
!!   NAME
!!     force_fields
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     Subroutine input_torsions changed (DG) 05-Dec-2000
!!     Output formats changed (DG) 05-Dec-2000
!!     JGH (26-01-2002) : force field parameters stored in tables, not in
!!       matrices. Input changed to have parameters labeled by the position
!!       and not atom pairs (triples etc)
!!
!!   SOURCE
!******************************************************************************

MODULE force_fields

  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: pi
  USE molecule_types,                  ONLY: bend_parameters_type,&
                                             bond_parameters_type,&
                                             intra_parameters_type,&
                                             torsion_parameters_type
  USE pair_potential_types,            ONLY: pair_potential_type
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE string_utilities,                ONLY: str_search,&
                                             uppercase
  USE termination,                     ONLY: stop_memory
  USE topology_types,                  ONLY: topology_parameters_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_force_field_section

  INTEGER :: scr
  TYPE ( global_environment_type ) :: globenv

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** force_fields/read_force_field_section [1.0] *
!!
!!   NAME
!!     read_force_field_section
!!
!!   FUNCTION
!!     Reads the force_field input section
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (30.11.2001) : moved determination of setup variables to 
!!                        molecule_input 
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_force_field_section ( topology, set_fn, intra_param, &
                                      potparm, globenv_input )


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    CHARACTER(LEN=*), INTENT(IN)             :: set_fn
    TYPE(intra_parameters_type), INTENT(OUT) :: intra_param
    TYPE(pair_potential_type), &
      DIMENSION(:, :), POINTER               :: potparm
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv_input

    CHARACTER(LEN=12)                        :: label
    CHARACTER(LEN=20)                        :: string
    CHARACTER(LEN=9)                         :: units
    INTEGER                                  :: ierror, ilen, ios, iw, &
                                                natom_type
    LOGICAL                                  :: units_set = .FALSE.

!------------------------------------------------------------------------------

  globenv = globenv_input
  scr = globenv % scr
  natom_type = topology % natom_type

  iw = scr

!..start the parser
  CALL parser_init(set_fn,globenv_input)

!..parse the input section for units
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror /= 0 ) THEN
     IF (globenv%ionode .AND. globenv%print_level >= 2) THEN
        WRITE ( iw, '( A )' ) ' No input section &FORCE_FIELD found on file '
        WRITE ( iw, '( T2, A )' ) set_fn
     END IF
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 20
        CALL get_next ( string, ilen )
        CALL uppercase ( string )
        SELECT CASE ( string )
        CASE ( 'UNITS' )
           ilen = 20
           CALL get_next ( units, ilen )
           CALL uppercase ( units )
           units_set = .TRUE.
           EXIT
        END SELECT
        CALL read_line
     END DO
  END IF
  CALL parser_end ()

  IF ( units_set ) THEN 
    CALL parser_init(set_fn,globenv_input)

!..read in non-bonded interaction information
    ALLOCATE (potparm(natom_type,natom_type),STAT=ios)
    IF ( ios /= 0 ) CALL stop_memory ( 'fist_input', &
       'potparm', natom_type ** 2 )
    CALL input_nonbonded(units,potparm,topology % atom_names,set_fn)

!..read the bond parameters
    NULLIFY (intra_param % bond_param)
    CALL input_bonds(units,intra_param % bond_param,set_fn)

!..read the bend parameters
    NULLIFY (intra_param % bend_param)
    CALL input_bends(units,intra_param % bend_param,set_fn)

!..read the torsion parameters
    NULLIFY (intra_param % torsion_param)
    CALL input_torsions(units,intra_param % torsion_param,set_fn)

!..read the distance constraint parameters
    NULLIFY (intra_param % constraint_distance)
    CALL input_distance_constraints (units, intra_param % constraint_distance, set_fn )

!..end the parser session
    CALL parser_end()
! 
    CALL convert_to_cp2k_units ( unit_char = units, potparm = potparm, intra_param = intra_param )
  ELSE
    CALL p_error()
    CALL stop_parser( 'read_force_field_section','UNITS not specified')
  ENDIF

END SUBROUTINE read_force_field_section

!!*****
!******************************************************************************
!!****** force_fields/input_nonbonded [1.0] *
!!
!!   NAME
!!     input_nonbonded
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I nonbonded                                                                 I
!!I    lennard-jones atom1 atom2 epsilon sigma rcut                           I
!!I    williams      atom1 atom2 a b c rcut                                   I
!!I    goodwin       atom1 atom2 vr0 d dc m mc rcut                           I
!!I end nonbonded                                                             I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_nonbonded (energy_unit, potparm, atom_names, set_fn )


    CHARACTER(LEN=*), INTENT(IN)             :: energy_unit
    TYPE(pair_potential_type), INTENT(INOUT) :: potparm( :, : )
    CHARACTER(LEN=*), INTENT(IN)             :: atom_names( : ), set_fn

    CHARACTER(LEN=10)                        :: at1, at2
    CHARACTER(LEN=12)                        :: label
    CHARACTER(LEN=40)                        :: string, string2
    CHARACTER(LEN=8)                         :: length_unit
    INTEGER                                  :: i1, i2, ierror, ilen, iw, &
                                                natom_types, np
    REAL(dbl)                                :: a, b, c, d, dc, epsilon, m, &
                                                mc, rcut, sigma, vr0

!------------------------------------------------------------------------------

  IF ( energy_unit == 'KCAL' ) length_unit = 'ANGSTROM'
  IF ( energy_unit == 'HARTREE' ) length_unit = 'BOHR'
  IF ( energy_unit == 'AU' ) length_unit = 'BOHR'
  IF ( energy_unit == 'ATOMIC' ) length_unit = 'BOHR'

  natom_types = SIZE ( atom_names )
  np = 0
  DO i1 = 1, natom_types
     DO i2 = 1, natom_types
        potparm(i1,i2) %type = 'NULL'
        !CALL uppercase(potparm(i1,i2) %type )
     END DO
  END DO

! initialize parameters in potparm
  potparm ( :, : ) % lj % epsilon = 0._dbl
  potparm ( :, : ) % lj % sigma6 = 0._dbl
  potparm ( :, : ) % lj % sigma12 = 0._dbl
  potparm ( :, : ) % willis % a = 0._dbl
  potparm ( :, : ) % willis % b = 0._dbl
  potparm ( :, : ) % willis % c = 0._dbl
  potparm ( :, : ) % goodwin % vr0 = 0._dbl
  potparm ( :, : ) % goodwin % d = 0._dbl
  potparm ( :, : ) % goodwin % dc = 0._dbl
  potparm ( :, : ) % goodwin % m = 0._dbl
  potparm ( :, : ) % goodwin % mc = 0._dbl
  iw = scr
!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror /= 0 ) THEN
     IF (globenv%ionode .AND. globenv%print_level >= 2) THEN
        WRITE ( iw, '( A )' ) ' No input section &FORCE_FIELD found on file '
        WRITE ( iw, '( T2, A )' ) set_fn
     END IF
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 9
        CALL get_next(string,ilen)
        CALL uppercase(string)
        IF (index(string,'NONBONDED') /= 0 ) THEN
           DO
              CALL read_line
              ilen = 0
              CALL get_next(string2,ilen)
              CALL uppercase(string2)
              SELECT CASE (string2)
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser( 'input_nonbonded','unknown option')
              CASE ( 'LENNARD-JONES')
                 np = np + 1
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( epsilon )
                 CALL get_next ( sigma )
                 CALL get_next ( rcut )
                 i1 = str_search(atom_names,natom_types,at1)
                 IF (i1==0) CALL stop_parser( 'input_nonbonded', &
                      'atom not found')
                 i2 = str_search(atom_names,natom_types,at2)
                 IF (i2==0) CALL stop_parser( 'input_nonbonded', &
                      'atom not found')
                 potparm(i1,i2) %type = string2
                 potparm(i1,i2) %lj%epsilon = epsilon
                 potparm(i1,i2) %lj%sigma6 = sigma**6
                 potparm(i1,i2) %lj%sigma12 = sigma**12
                 potparm(i1,i2) %rcutsq = rcut*rcut
                 potparm(i2,i1) = potparm(i1,i2)
              CASE ( 'WILLIAMS')
                 np = np + 1
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( a )
                 CALL get_next ( b )
                 CALL get_next ( c )
                 CALL get_next ( rcut )
                 i1 = str_search(atom_names,natom_types,at1)
                 IF (i1==0) CALL stop_parser( 'input_nonbonded', &
                      'atom not found')
                 i2 = str_search(atom_names,natom_types,at2)
                 IF (i2==0) CALL stop_parser( 'input_nonbonded', &
                      'atom not found')
                 potparm(i1,i2) %type = string2
                 potparm(i1,i2) %willis%a = a
                 potparm(i1,i2) %willis%b = b
                 potparm(i1,i2) %willis%c = c
                 potparm(i1,i2) %rcutsq = rcut*rcut
                 potparm(i2,i1) = potparm(i1,i2)
              CASE ( 'GOODWIN')
                 np = np + 1
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( vr0 )
                 CALL get_next ( d )
                 CALL get_next ( dc )
                 CALL get_next ( m )
                 CALL get_next ( mc )
                 CALL get_next ( rcut )
                 i1 = str_search(atom_names,natom_types,at1)
                 IF (i1==0) CALL stop_parser( 'input_nonbonded', &
                      'atom not found')
                 i2 = str_search(atom_names,natom_types,at2)
                 IF (i2==0) CALL stop_parser( 'input_nonbonded', &
                      'atom not found')
                 potparm(i1,i2) %type = string2
                 potparm(i1,i2) %goodwin%vr0 = vr0
                 potparm(i1,i2) %goodwin%d = d
                 potparm(i1,i2) %goodwin%dc = dc
                 potparm(i1,i2) %goodwin%m = m
                 potparm(i1,i2) %goodwin%mc = mc
                 potparm(i1,i2) %rcutsq = rcut*rcut
                 potparm(i2,i1) = potparm(i1,i2)
              CASE ( 'END')
                 ilen = 0
                 CALL get_next(string2,ilen)
                 EXIT
              END SELECT
           END DO
        END IF
        CALL read_line
     END DO
  END IF

!..write some information to output
  IF (globenv%ionode) THEN
     IF (globenv % print_level>=0 .AND. np>0) THEN
        WRITE ( iw, '( A,A )' ) ' FORCE FIELD| Non-bonded interaction parameters in ', energy_unit
        DO i1 = 1, natom_types
           DO i2 = i1, natom_types
              IF (potparm(i1,i2)%type=='LENNARD-JONES') THEN
                 WRITE ( iw, '( A,A,T30,A,T37,A,T58,A,T71,F10.4 )' ) &
                      ' FORCE FIELD| ', 'LJ (12-6)', &
                      adjustl(atom_names(i1)), adjustl(atom_names(i2)), &
                      ' epsilon ', potparm(i1,i2) %lj%epsilon
                 WRITE ( iw, '( A,T58,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
                      '   sigma ', (potparm(i1,i2)%lj%sigma6)**(1./6.)
                 WRITE ( iw, '( A,T58,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
                      '    rcut ', sqrt(potparm(i1,i2)%rcutsq)
              ELSE IF (potparm(i1,i2)%type=='WILLIAMS') THEN
                 WRITE ( iw, '( A,A,T30,A,T37,A,T55,A,T68,E13.4 )' ) &
                      ' FORCE FIELD| ', 'Williams', adjustl(atom_names(i1)), &
                      adjustl(atom_names(i2)), 'A', potparm(i1,i2) %willis%a
                 WRITE ( iw, '( A,T47,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
                      '    B', potparm(i1,i2) %willis%b
                 WRITE ( iw, '( A,T47,A,T68,E13.4 )' ) ' FORCE FIELD| ', &
                      '   C', potparm(i1,i2) %willis%c
                 WRITE ( iw, '( A,T47,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
                      '    rcut ', sqrt(potparm(i1,i2)%rcutsq)
              ELSE IF (potparm(i1,i2)%type=='GOODWIN') THEN
                 WRITE ( iw, '( A,A,T30,A,T37,A,T55,A,T68,E13.4 )' ) &
                      ' FORCE FIELD| ', 'Goodwin', adjustl(atom_names(i1)), &
                      adjustl(atom_names(i2)), 'Ve0', &
                      potparm(i1,i2) %goodwin%vr0
                 WRITE ( iw, '( A,T44,A,T61,2F10.4 )' ) ' FORCE FIELD| ', &
                      '  d,dc ', potparm(i1,i2) %goodwin%d, &
                      potparm(i1,i2) %goodwin%dc
                 WRITE ( iw, '( A,T44,A,T61,2F10.4 )' ) ' FORCE FIELD| ', &
                      '  m,mc   ', potparm(i1,i2) %goodwin%m, &
                      potparm(i1,i2) %goodwin%mc
                 WRITE ( iw, '( A,T44,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
                      '  rcut', sqrt(potparm(i1,i2)%rcutsq)
              END IF
           END DO
        END DO
        WRITE ( iw, '( )' )
     END IF
  END IF

END SUBROUTINE input_nonbonded

!!*****
!******************************************************************************
!!****** force_fields/input_bonds [1.0] *
!!
!!   NAME
!!     input_bonds
!!
!!   FUNCTION
!!     Read all the bond parameters. Put them in the bond_parm table.
!!     While we read the bond_pot_type,
!!     we (so far) will do nothing with this information. It appears that
!!     all major parameters sets, such as CHARMM or AMBER, all have only
!!     one type of bonded interaction (harmonic).
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) Bond parameters are now stored in tables. The position
!!        within the table is used as handle for the topology
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I bonds nbonds                                                              I
!!I    harmonic bp1 bp2                                                       I
!!I end bonds                                                                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_bonds(energy_unit,bond_parm,set_fn)

    CHARACTER(LEN=*), INTENT(IN)             :: energy_unit
    TYPE(bond_parameters_type), POINTER      :: bond_parm( : )
    CHARACTER(LEN=*), INTENT(IN)             :: set_fn

    CHARACTER(LEN=12)                        :: label
    CHARACTER(LEN=40)                        :: string, string2
    INTEGER                                  :: i, ierror, ilen, ios, iw, nb
    REAL(dbl)                                :: bp1, bp2

!------------------------------------------------------------------------------

  nb = 0
  iw = scr
!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror == 0 ) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 5
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        IF (index(string,'BONDS') /= 0 ) THEN
           CALL get_next(nb)
           ALLOCATE (bond_parm(nb),STAT=ios)
           IF ( ios /= 0 ) CALL stop_memory ( 'input_bonds', 'bond_parm', nb )
           bond_parm ( : ) %type = 'null'
           bond_parm ( : ) %r0 = 0._dbl
           bond_parm ( : ) %k = 0._dbl
           i = 0
           DO
              CALL read_line
              ilen = 0
              CALL get_next(string2,ilen)
              CALL uppercase(string2)
              WRITE(*,*) "string2 ",string2
              SELECT CASE (string2)
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser( 'input_bond','unknown option')
              CASE ( 'HARMONIC')
                 i = i + 1
                 CALL get_next ( bp1 )
                 CALL get_next ( bp2 )
                 bond_parm(i) %type = string2
                 bond_parm(i) %r0 = bp1
                 bond_parm(i) %k = bp2
              CASE ( 'END')
                 ilen = 0
                 CALL get_next(string2,ilen)
                 EXIT
              END SELECT
           END DO
        END IF
        CALL read_line
     END DO
  END IF

!..write some information to output
  IF ( globenv % ionode ) THEN
     IF (globenv % print_level>=0 .AND. nb>0) THEN
        WRITE ( iw, '( A,A )' ) ' FORCE FIELD| Bond interaction parameters in ', energy_unit
        DO i = 1, nb
          IF (index(bond_parm(i)%type,'HARMONIC') /= 0 ) THEN
                 WRITE ( iw, '( A, T30, A, I5, T58, A, T71, F10.2 )' ) &
                      ' FORCE FIELD| harmonic', ' bond nr. ',i,&
                      '      r0 ', bond_parm ( i ) % r0
                 WRITE ( iw, '( A, T58, A, T71, F10.2 )' ) ' FORCE FIELD| ', &
                      ' k ', bond_parm ( i ) % k
          END IF
        END DO
        WRITE ( iw, '( )' )
     END IF
  END IF

END SUBROUTINE input_bonds

!!*****
!******************************************************************************
!!****** force_fields/input_bends [1.0] *
!!
!!   NAME
!!     input_bends
!!
!!   FUNCTION
!!     Read all the bend parameters. Put them in the bend_parm table.
!!     While we read the bend_pot_type,
!!     we (so far) will do nothing with this information. It appears that
!!     all major parameters sets, such as CHARMM or AMBER, all have only
!!     only one type of bend interaction (harmonic).
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) Bend parameters are now stored in tables. The position
!!        within the table is used as handle for the topology
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I bends nbends                                                              I
!!I    harmonic bp1 bp2                                                       I
!!I end bends                                                                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_bends ( energy_unit, bend_parm, set_fn )


    CHARACTER(LEN=*), INTENT(IN)             :: energy_unit
    TYPE(bend_parameters_type), POINTER      :: bend_parm(:)
    CHARACTER(LEN=*), INTENT(IN)             :: set_fn

    CHARACTER(LEN=12)                        :: label
    CHARACTER(LEN=40)                        :: string, string2
    INTEGER                                  :: i, ierror, ilen, ios, iw, nb
    REAL(dbl)                                :: bp1, bp2, todeg, torad

!------------------------------------------------------------------------------

  torad = pi / 180.0_dbl
  todeg = 180.0_dbl/pi

  nb = 0
  iw = scr
!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror == 0 ) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 5
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        IF (index(string,'BENDS') /= 0 ) THEN
           CALL get_next ( nb )
           ALLOCATE (bend_parm(nb),STAT=ios)
           IF ( ios /= 0 ) CALL stop_memory ( 'input_bends', 'bend_parm', nb )
           bend_parm(:) %type = 'null'
           bend_parm(:) %theta0 = 0._dbl
           bend_parm(:) %k = 0._dbl
           i=0
           DO
              CALL read_line
              ilen = 0
              CALL get_next(string2,ilen)
              CALL uppercase(string2)
              SELECT CASE (string2)
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser( 'input_bond','unknown option')
              CASE ( 'HARMONIC')
                 i = i + 1
                 CALL get_next ( bp1 )
                 CALL get_next ( bp2 )
                 bend_parm(i) %type = string2
                 bend_parm(i) %theta0 = bp1*torad
                 bend_parm(i) %k = bp2
              CASE ( 'END')
                 ilen = 0
                 CALL get_next(string2,ilen)
                 EXIT
              END SELECT
           END DO
        END IF
        CALL read_line
     END DO
  END IF

!..write some information to output
  IF ( globenv % ionode ) THEN
     IF (globenv % print_level>=0 .AND. nb>0) THEN
        WRITE ( iw, '( A, A )' ) ' FORCE FIELD| Bend interaction parameters in ', energy_unit
        DO i = 1, nb
           IF (index(bend_parm(i)%type,'HARMONIC') /= 0 ) THEN
              WRITE ( iw, '( A,T30,A,I5,T58,A,T71,F10.2 )' ) &
                   ' FORCE FIELD| harmonic', ' bend nr ',i, &
                   'theta0 ', bend_parm(i) %theta0*todeg
              WRITE ( iw, '( A,T58,A,T71,F10.2 )' ) ' FORCE FIELD| ', &
                   '       k ', bend_parm(i) %k
           END IF
        END DO
        WRITE ( iw, '( )' )
     END IF
  END IF

END SUBROUTINE input_bends

!!*****
!******************************************************************************
!!****** force_fields/input_torsions [1.0] *
!!
!!   NAME
!!     input_torsions
!!
!!   FUNCTION
!!     Read all the torsion parameters. Put them in the torsion_parm table.
!!     We have defined the
!!     torsion_parameter type with 7 coefficients for 0 to 6. Hopefully
!!     this will be generic enough to handle all torsion potentials.
!!     Since different types will be used, we need to store the type of torsion
!!     as well.
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) Torsion parameters are now stored in tables. The position
!!        within the table is used as handle for the topology
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I torsions  ntors                                                           I
!!I    jorgensen tp[1-4]                                                      I
!!I end torsions                                                              I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_torsions ( energy_unit,  torsion_parm, set_fn )


    CHARACTER(LEN=*), INTENT(IN)             :: energy_unit
    TYPE(torsion_parameters_type), POINTER   :: torsion_parm( : )
    CHARACTER(LEN=*), INTENT(IN)             :: set_fn

    CHARACTER(LEN=12)                        :: label
    CHARACTER(LEN=20)                        :: type
    CHARACTER(LEN=40)                        :: string, string2
    INTEGER                                  :: i, ierror, ilen, ios, iw, j, &
                                                nb
    REAL(dbl), DIMENSION(0:3)                :: tp

!------------------------------------------------------------------------------

  nb = 0
  iw = scr
!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror == 0 ) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 8
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        IF (index(string,'TORSIONS') /= 0 ) THEN
           CALL get_next ( nb )
           ALLOCATE ( torsion_parm ( nb ), STAT = ios )
           IF ( ios /= 0 ) CALL stop_memory ( 'input_torsions', &
              'torsion_parm', nb )
           torsion_parm ( : ) % type = 'null'
           DO i = 0, SIZE ( torsion_parm(1) % a ) - 1
              torsion_parm(:) %a(i) = 0._dbl
           END DO
           i = 0
           DO
              CALL read_line
              ilen = 0
              CALL get_next(string2,ilen)
              CALL uppercase(string2)
              SELECT CASE (string2)
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser( 'input_torsion','unknown option')
              CASE ( 'JORGENSEN')
                 i = i + 1
                 torsion_parm(i) %type = string2
                 DO j = 0, SIZE(tp) - 1
                    CALL get_next ( torsion_parm(i) %a(j) )
                 END DO
              CASE ( 'END')
                 ilen = 0
                 CALL get_next(string2,ilen)
                 EXIT
              END SELECT
           END DO
        END IF
        CALL read_line()
     END DO
  END IF

!..write some information to output
  IF ( globenv % ionode ) THEN
     IF ( globenv % print_level >= 0 .AND. nb > 0 ) THEN
        WRITE ( iw, '( A, A )' ) &
             ' FORCE FIELD| Torsion interaction parameters in ', energy_unit
        DO i = 1, nb
           IF ( index ( torsion_parm(i) % type, 'JORGENSEN' ) /= 0 ) THEN
              WRITE ( iw, '( A, T30, A, i5, T58, A, T71, F10.2 )' ) &
                ' FORCE FIELD| Jorgensen', ' torsion nr ',i, &
                            '      a0 ', torsion_parm ( i ) % a(0)
                       WRITE ( iw, '( A, T58, A, T71, F10.2 )' ) ' FORCE FIELD| ', &
                            '      a1 ', torsion_parm ( i ) % a(1)
                       WRITE ( iw, '( A, T58, A, T71, F10.2 )' ) ' FORCE FIELD| ', &
                            '      a2 ', torsion_parm ( i ) % a(2)
                       WRITE ( iw, '( A, T58, A, T71, F10.2 )' ) ' FORCE FIELD| ', &
                            '      a3 ', torsion_parm ( i ) % a(3)
           END IF
        END DO
        WRITE ( iw, '( )' )
     END IF
  END IF

END SUBROUTINE input_torsions

!!*****
!******************************************************************************
!!****** force_fields/input_distance_constraints [1.0] *
!!
!!   NAME
!!     input_distance_constraints
!!
!!   FUNCTION
!!     Read all the distance parameters. Put them in the
!!     constraint_distance array.
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) Distance parameters are now stored in tables. The position
!!        within the table is used as handle for the topology
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I distcon  ndist                                                            I
!!I     dp                                                                    I
!!I end distcon                                                               I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_distance_constraints ( energy_unit, constraint_distance, set_fn )


    CHARACTER(LEN=*), INTENT(IN)             :: energy_unit
    REAL(dbl), POINTER                       :: constraint_distance( : )
    CHARACTER(LEN=*), INTENT(IN)             :: set_fn

    CHARACTER(LEN=12)                        :: label
    CHARACTER(LEN=40)                        :: string
    INTEGER                                  :: i, ierror, ilen, ios, iw, nb

!------------------------------------------------------------------------------

  nb = 0
  iw = scr
!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror == 0 ) THEN
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 7
        CALL get_next ( string, ilen )
        CALL uppercase ( string )

        IF ( INDEX ( string, 'DISTCON' ) /= 0 ) THEN
           CALL get_next ( nb )
           ALLOCATE (constraint_distance(nb),STAT=ios)
           IF ( ios /= 0 ) CALL stop_memory ( 'input_distance_constraints', &
                'constraint_distance', nb )
           constraint_distance ( : ) = 0.0_dbl
           CALL read_line()
           DO i = 1, nb
              CALL get_next ( constraint_distance ( i ) )
           END DO
        END IF
        CALL read_line
     END DO
  END IF

!..write some information to output
  IF ( globenv % ionode ) THEN
     IF ( globenv % print_level >= 0 .AND. nb > 0 ) THEN
        WRITE ( iw, '( A, A )' ) ' FORCE FIELD| Distance constraints in ', energy_unit
        DO i = 1, nb
          IF ( constraint_distance ( i ) /= 0.0_dbl ) THEN
            WRITE ( iw, '( A, T30, A, T37, I6, T71, F10.4 )' ) &
              ' FORCE FIELD| Constraint: ', ' distcon nr. ',i, &
                      constraint_distance ( i )
          END IF
        END DO
        WRITE ( iw, '()' )
     END IF
  END IF

END SUBROUTINE input_distance_constraints

!!*****
!******************************************************************************

END MODULE force_fields

!******************************************************************************
