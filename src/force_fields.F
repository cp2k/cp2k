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
  USE kinds,                           ONLY: dbl,&
                                             default_string_length
  USE mathconstants,                   ONLY: pi
  USE memory_utilities,                ONLY: reallocate
  USE pair_potential_types,            ONLY: reallocate_potparm
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
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
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

SUBROUTINE read_force_field_section ( topology , globenv )


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=default_string_length)                        :: label
    CHARACTER(LEN=default_string_length)                        :: string
    CHARACTER(LEN=default_string_length)                         :: units
    INTEGER                                  :: ierror, ilen, ios, iw, &
                                                natom_type
    LOGICAL                                  :: units_set = .FALSE.

!------------------------------------------------------------------------------

  scr = globenv % scr
  natom_type = topology % natom_type

  iw = scr

!!..start the parser
!  CALL parser_init(globenv%input_file_name,globenv)
!
!!..parse the input section for units
!  label = '&FORCE_FIELD'
!  CALL search_label ( label, ierror, ignore_case=.TRUE. )
!  IF (ierror /= 0 ) THEN
!     IF (globenv%ionode .AND. globenv%print_level >= 2) THEN
!        WRITE ( iw, '( A )' ) ' No input section &FORCE_FIELD found on file '
!        WRITE ( iw, '( T2, A )' ) globenv%input_file_name
!     END IF
!  ELSE
!     CALL read_line
!     DO WHILE ( test_next() /= 'X' )
!        ilen = 20
!        CALL get_next ( string, ilen )
!        CALL uppercase ( string )
!        SELECT CASE ( string )
!        CASE ( 'UNITS' )
!           ilen = 20
!           CALL get_next ( units, ilen )
!           CALL uppercase ( units )
!           units_set = .TRUE.
!           EXIT
!        END SELECT
!        CALL read_line
!     END DO
!  END IF
!  CALL parser_end ()
!
!  IF ( units_set ) THEN 
    CALL parser_init(globenv%input_file_name,globenv)

!!..read in non-bonded interaction information
    CALL input_nonbonded(topology,globenv)
!
!..read the bond parameters
    CALL input_bonds(topology,globenv)

!..read the bend parameters
    CALL input_bends(topology,globenv)

!..read the torsion parameters
    CALL input_torsions(topology,globenv)

!..read the distance constraint parameters
    CALL input_constraints(topology,globenv)

!..end the parser session
    CALL parser_end()
 
!  ELSE
!    CALL p_error()
!    CALL stop_parser( 'read_force_field_section','UNITS not specified')
!  ENDIF

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

SUBROUTINE input_nonbonded(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=default_string_length)     :: at1, at2
    CHARACTER(LEN=default_string_length)     :: label
    CHARACTER(LEN=default_string_length)     :: string, string2
    INTEGER                                  :: i, ierror, ilen, iw, &
                                                natom_types, np
    REAL(dbl)                                :: a, b, c, d, dc, epsilon, m, &
                                                mc, rcut, sigma, vr0

!------------------------------------------------------------------------------

  np = 0
  iw = scr
!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror /= 0 ) THEN
     IF (globenv%ionode .AND. globenv%print_level >= 2) THEN
        WRITE ( iw, '( A )' ) ' No input section &FORCE_FIELD found on file '
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
                 CALL reallocate_potparm(topology%ff_nonbonded,1,np)
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( epsilon )
                 CALL get_next ( sigma )
                 CALL get_next ( rcut )
                 topology%ff_nonbonded(np) %type = string2
                 topology%ff_nonbonded(np) %at1 = at1
                 topology%ff_nonbonded(np) %at2 = at2
                 topology%ff_nonbonded(np) %lj%epsilon = epsilon
                 topology%ff_nonbonded(np) %lj%sigma6 = sigma**6
                 topology%ff_nonbonded(np) %lj%sigma12 = sigma**12
                 topology%ff_nonbonded(np) %rcutsq = rcut*rcut
              CASE ( 'WILLIAMS')
                 np = np + 1
                 CALL reallocate_potparm(topology%ff_nonbonded,1,np)
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( a )
                 CALL get_next ( b )
                 CALL get_next ( c )
                 CALL get_next ( rcut )
                 topology%ff_nonbonded(np) %type = string2
                 topology%ff_nonbonded(np) %at1 = at1
                 topology%ff_nonbonded(np) %at2 = at2
                 topology%ff_nonbonded(np) %willis%a = a
                 topology%ff_nonbonded(np) %willis%b = b
                 topology%ff_nonbonded(np) %willis%c = c
                 topology%ff_nonbonded(np) %rcutsq = rcut*rcut
              CASE ( 'GOODWIN')
                 np = np + 1
                 CALL reallocate_potparm(topology%ff_nonbonded,1,np)
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
                 topology%ff_nonbonded(np) %type = string2
                 topology%ff_nonbonded(np) %at1 = at1
                 topology%ff_nonbonded(np) %at2 = at2
                 topology%ff_nonbonded(np) %goodwin%vr0 = vr0
                 topology%ff_nonbonded(np) %goodwin%d = d
                 topology%ff_nonbonded(np) %goodwin%dc = dc
                 topology%ff_nonbonded(np) %goodwin%m = m
                 topology%ff_nonbonded(np) %goodwin%mc = mc
                 topology%ff_nonbonded(np) %rcutsq = rcut*rcut
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

!!..write some information to output
!  IF (globenv%ionode) THEN
!     IF (globenv % print_level>=0 .AND. np>0) THEN
!        WRITE ( iw, '( A,A )' ) ' FORCE FIELD| Non-bonded interaction parameters in ', energy_unit
!        DO i1 = 1, natom_types
!           DO i2 = i1, natom_types
!              IF (potparm(i1,i2)%type=='LENNARD-JONES') THEN
!                 WRITE ( iw, '( A,A,T30,A,T37,A,T58,A,T71,F10.4 )' ) &
!                      ' FORCE FIELD| ', 'LJ (12-6)', &
!                      adjustl(atom_names(i1)), adjustl(atom_names(i2)), &
!                      ' epsilon ', potparm(i1,i2) %lj%epsilon
!                 WRITE ( iw, '( A,T58,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
!                      '   sigma ', (potparm(i1,i2)%lj%sigma6)**(1./6.)
!                 WRITE ( iw, '( A,T58,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
!                      '    rcut ', sqrt(potparm(i1,i2)%rcutsq)
!              ELSE IF (potparm(i1,i2)%type=='WILLIAMS') THEN
!                 WRITE ( iw, '( A,A,T30,A,T37,A,T55,A,T68,E13.4 )' ) &
!                      ' FORCE FIELD| ', 'Williams', adjustl(atom_names(i1)), &
!                      adjustl(atom_names(i2)), 'A', potparm(i1,i2) %willis%a
!                 WRITE ( iw, '( A,T47,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
!                      '    B', potparm(i1,i2) %willis%b
!                 WRITE ( iw, '( A,T47,A,T68,E13.4 )' ) ' FORCE FIELD| ', &
!                      '   C', potparm(i1,i2) %willis%c
!                 WRITE ( iw, '( A,T47,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
!                      '    rcut ', sqrt(potparm(i1,i2)%rcutsq)
!              ELSE IF (potparm(i1,i2)%type=='GOODWIN') THEN
!                 WRITE ( iw, '( A,A,T30,A,T37,A,T55,A,T68,E13.4 )' ) &
!                      ' FORCE FIELD| ', 'Goodwin', adjustl(atom_names(i1)), &
!                      adjustl(atom_names(i2)), 'Ve0', &
!                      potparm(i1,i2) %goodwin%vr0
!                 WRITE ( iw, '( A,T44,A,T61,2F10.4 )' ) ' FORCE FIELD| ', &
!                      '  d,dc ', potparm(i1,i2) %goodwin%d, &
!                      potparm(i1,i2) %goodwin%dc
!                 WRITE ( iw, '( A,T44,A,T61,2F10.4 )' ) ' FORCE FIELD| ', &
!                      '  m,mc   ', potparm(i1,i2) %goodwin%m, &
!                      potparm(i1,i2) %goodwin%mc
!                 WRITE ( iw, '( A,T44,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
!                      '  rcut', sqrt(potparm(i1,i2)%rcutsq)
!              END IF
!           END DO
!        END DO
!        WRITE ( iw, '( )' )
!     END IF
!  END IF

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
!!I bonds
!!I    harmonic ATMNAE1 ATMNAME2 k r0
!!I end bonds                                                                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_bonds(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=default_string_length)              :: label
    CHARACTER(LEN=default_string_length)              :: string, string2
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
                 topology%ff = .TRUE.
                 i = i + 1
                 CALL reallocate(topology%ff_bond_a,1,i)
                 CALL reallocate(topology%ff_bond_b,1,i)
                 CALL reallocate(topology%ff_bond_k,1,i)
                 CALL reallocate(topology%ff_bond_r0,1,i)

                 ilen = default_string_length
                 CALL get_next ( topology%ff_bond_a(i), ilen )
                 CALL get_next ( topology%ff_bond_b(i), ilen )
                 CALL get_next ( topology%ff_bond_k(i) )
                 CALL get_next ( topology%ff_bond_r0(i) )
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

!!..write some information to output
!  IF ( globenv % ionode ) THEN
!     IF (globenv % print_level>=0 .AND. nb>0) THEN
!        WRITE ( iw, '( A,A )' ) ' FORCE FIELD| Bond interaction parameters in ', energy_unit
!        DO i = 1, nb
!          IF (index(bond_parm(i)%type,'HARMONIC') /= 0 ) THEN
!                 WRITE ( iw, '( A, T30, A, I5, T58, A, T71, F10.2 )' ) &
!                      ' FORCE FIELD| harmonic', ' bond nr. ',i,&
!                      '      r0 ', bond_parm ( i ) % r0
!                 WRITE ( iw, '( A, T58, A, T71, F10.2 )' ) ' FORCE FIELD| ', &
!                      ' k ', bond_parm ( i ) % k
!          END IF
!        END DO
!        WRITE ( iw, '( )' )
!     END IF
!  END IF

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
!!I bends
!!I    harmonic bp1 bp2                                                       I
!!I end bends                                                                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_bends(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

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
                 topology%ff = .TRUE.
                 i = i + 1
                 CALL reallocate(topology%ff_bend_a,1,i)
                 CALL reallocate(topology%ff_bend_b,1,i)
                 CALL reallocate(topology%ff_bend_c,1,i)
                 CALL reallocate(topology%ff_bend_k,1,i)
                 CALL reallocate(topology%ff_bend_theta0,1,i)

                 ilen = default_string_length
                 CALL get_next ( topology%ff_bend_a(i), ilen )
                 CALL get_next ( topology%ff_bend_b(i), ilen )
                 CALL get_next ( topology%ff_bend_c(i), ilen )
                 CALL get_next ( topology%ff_bend_k(i) )
                 CALL get_next ( topology%ff_bend_theta0(i) )
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

!!..write some information to output
!  IF ( globenv % ionode ) THEN
!     IF (globenv % print_level>=0 .AND. nb>0) THEN
!        WRITE ( iw, '( A, A )' ) ' FORCE FIELD| Bend interaction parameters in ', energy_unit
!        DO i = 1, nb
!           IF (index(bend_parm(i)%type,'HARMONIC') /= 0 ) THEN
!              WRITE ( iw, '( A,T30,A,I5,T58,A,T71,F10.2 )' ) &
!                   ' FORCE FIELD| harmonic', ' bend nr ',i, &
!                   'theta0 ', bend_parm(i) %theta0*todeg
!              WRITE ( iw, '( A,T58,A,T71,F10.2 )' ) ' FORCE FIELD| ', &
!                   '       k ', bend_parm(i) %k
!           END IF
!        END DO
!        WRITE ( iw, '( )' )
!     END IF
!  END IF

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
!!I torsions
!!I    jorgensen tp[1-4]                                                      I
!!I end torsions                                                              I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_torsions(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=default_string_length)  :: label
    CHARACTER(LEN=default_string_length)  :: type
    CHARACTER(LEN=default_string_length)  :: string, string2
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
                 topology%ff = .TRUE.
                 i = i + 1
                 CALL reallocate(topology%ff_torsion_a,1,i)
                 CALL reallocate(topology%ff_torsion_b,1,i)
                 CALL reallocate(topology%ff_torsion_c,1,i)
                 CALL reallocate(topology%ff_torsion_d,1,i)
                 CALL reallocate(topology%ff_torsion_k,1,i)
                 CALL reallocate(topology%ff_torsion_phi0,1,i)

                 ilen = default_string_length
                 CALL get_next ( topology%ff_torsion_a(i), ilen )
                 CALL get_next ( topology%ff_torsion_b(i), ilen )
                 CALL get_next ( topology%ff_torsion_c(i), ilen )
                 CALL get_next ( topology%ff_torsion_d(i), ilen )
                 CALL get_next ( topology%ff_torsion_k(i) )
                 CALL get_next ( topology%ff_torsion_phi0(i) )
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

!!..write some information to output
!  IF ( globenv % ionode ) THEN
!     IF ( globenv % print_level >= 0 .AND. nb > 0 ) THEN
!        WRITE ( iw, '( A, A )' ) &
!             ' FORCE FIELD| Torsion interaction parameters in ', energy_unit
!        DO i = 1, nb
!           IF ( index ( torsion_parm(i) % type, 'JORGENSEN' ) /= 0 ) THEN
!              WRITE ( iw, '( A, T30, A, i5, T58, A, T71, F10.2 )' ) &
!                ' FORCE FIELD| Jorgensen', ' torsion nr ',i, &
!                            '      a0 ', torsion_parm ( i ) % a(0)
!                       WRITE ( iw, '( A, T58, A, T71, F10.2 )' ) ' FORCE FIELD| ', &
!                            '      a1 ', torsion_parm ( i ) % a(1)
!                       WRITE ( iw, '( A, T58, A, T71, F10.2 )' ) ' FORCE FIELD| ', &
!                            '      a2 ', torsion_parm ( i ) % a(2)
!                       WRITE ( iw, '( A, T58, A, T71, F10.2 )' ) ' FORCE FIELD| ', &
!                            '      a3 ', torsion_parm ( i ) % a(3)
!           END IF
!        END DO
!        WRITE ( iw, '( )' )
!     END IF
!  END IF

END SUBROUTINE input_torsions

!!*****
!******************************************************************************
!!****** force_fields/input_constraints [1.0] *
!!
!!   NAME
!!     input_constraints
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
!!I constraints                                                               I
!!I   HYDROGEN                                                                I
!!I   ATOM ??????                                                             I
!!I   DIST MOLKIND ATM_A ATM_B D_AB                                           I
!!I   G3x3 MOLKIND ATM_A ATM_B ATM_C D_AB D_AC D_BC                           I
!!I                                                                           I
!!I end constraints                                                           I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_constraints(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

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

        IF ( INDEX ( string, 'ATOM' ) /= 0 ) THEN
           topology%constraint=.TRUE.
           topology%const_atom=.TRUE.
           CALL stop_program ( "input_constraints", &
                               "CONATOM not done yet " )
        ELSE IF ( INDEX ( string, 'CONHYDR' ) /= 0 ) THEN
           topology%constraint=.TRUE.
           topology%const_hydr=.TRUE.
           CALL stop_program ( "input_constraints", &
               "CONHYDR will not work until force field is completed " )
        ELSE IF ( INDEX ( string, 'DIST' ) /= 0 ) THEN
           CALL read_line()
           topology%constraint=.TRUE.
           topology%const_dist=.TRUE.
           ilen = 40
           CALL get_next(string,ilen)
           CALL uppercase(string)
           i=topology%nconst_dist + 1
           CALL reallocate(topology%const_dist_mol,1,i)
           CALL reallocate(topology%const_dist_a,1,i)
           CALL reallocate(topology%const_dist_b,1,i)
           CALL reallocate(topology%const_dist_dab,1,i)
           CALL get_next(topology%const_dist_mol(i))
           CALL get_next(topology%const_dist_a(i))
           CALL get_next(topology%const_dist_b(i))
           CALL get_next(topology%const_dist_dab(i))
           topology%nconst_dist = i
        ELSE IF ( INDEX ( string, 'G3X3' ) /= 0 ) THEN
           topology%constraint=.TRUE.
           topology%const_33=.TRUE.
           ilen = 40
           CALL get_next(string,ilen)
           CALL uppercase(string)
           i=topology%nconst_g33 + 1
           CALL reallocate(topology%const_g33_mol,1,i)
           CALL reallocate(topology%const_g33_a,1,i)
           CALL reallocate(topology%const_g33_b,1,i)
           CALL reallocate(topology%const_g33_c,1,i)
           CALL reallocate(topology%const_g33_dab,1,i)
           CALL reallocate(topology%const_g33_dac,1,i)
           CALL reallocate(topology%const_g33_dbc,1,i)
           CALL get_next(topology%const_g33_mol(i))
           CALL get_next(topology%const_g33_a(i))
           CALL get_next(topology%const_g33_b(i))
           CALL get_next(topology%const_g33_c(i))
           CALL get_next(topology%const_g33_dab(i))
           CALL get_next(topology%const_g33_dac(i))
           CALL get_next(topology%const_g33_dbc(i))
        END IF
        CALL read_line
     END DO
  END IF

!!..write some information to output
!  IF ( globenv % ionode ) THEN
!     IF ( globenv % print_level >= 0 .AND. nb > 0 ) THEN
!        WRITE ( iw, '( A, A )' ) ' FORCE FIELD| Distance constraints in ', energy_unit
!        DO i = 1, nb
!          IF ( constraint_distance ( i ) /= 0.0_dbl ) THEN
!            WRITE ( iw, '( A, T30, A, T37, I6, T71, F10.4 )' ) &
!              ' FORCE FIELD| Constraint: ', ' distcon nr. ',i, &
!                      constraint_distance ( i )
!          END IF
!        END DO
!        WRITE ( iw, '()' )
!     END IF
!  END IF
!
END SUBROUTINE input_constraints

!!*****
!******************************************************************************

END MODULE force_fields

!******************************************************************************
