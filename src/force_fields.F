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
!!
!!   SOURCE
!******************************************************************************

MODULE force_fields

  USE kinds, ONLY : dbl
  USE global_types, ONLY : global_environment_type
  USE input_types, ONLY : setup_parameters_type
  USE particle_types, ONLY : particle_prop_type
  USE pair_potential, ONLY : potentialparm_type
  USE mathconstants, ONLY : pi
  USE molecule_types, ONLY : bond_parameters_type, bend_parameters_type, &
       torsion_parameters_type, intra_parameters_type, &
       molecule_type
  USE message_passing, ONLY : mp_bcast
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       p_error, get_next, stop_parser
  USE string_utilities, ONLY : uppercase, xstring, str_search,  &
       str_comp, make_tuple
  USE termination, ONLY : stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_force_field_section

!*apsi* Because of a feature(?) in DEC compiler
  PUBLIC :: ATOMNAMESLENGTH
  INTEGER, PARAMETER :: ATOMNAMESLENGTH = 20

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
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_force_field_section ( setup, mol_setup, set_fn, intra_param, &
     potparm, atom_names, pstat, globenv_input )

  IMPLICIT NONE

! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: set_fn
  CHARACTER ( LEN = ATOMNAMESLENGTH ), DIMENSION ( : ), POINTER :: atom_names
  TYPE (intra_parameters_type ), INTENT ( OUT ) :: intra_param
  TYPE (potentialparm_type ), DIMENSION ( :, : ), POINTER :: potparm
  TYPE (molecule_type ), DIMENSION ( : ), INTENT ( IN ) :: mol_setup
  TYPE (particle_prop_type ), DIMENSION ( : ), POINTER :: pstat
  TYPE (setup_parameters_type ), INTENT ( INOUT ) :: setup
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv_input

! Locals
  INTEGER :: i, j, n, natom_types, ios, nmol_type, iw

!------------------------------------------------------------------------------

  globenv = globenv_input
  scr = globenv % scr

  iw = scr

!..count total number of molecules and atoms
  setup % nmol = 0
  setup % natoms = 0
  setup % natom_type = 0
  nmol_type = size(mol_setup)
  DO i = 1, nmol_type
     setup % nmol = setup % nmol + mol_setup(i) % num_mol
     setup % natoms = setup % natoms + mol_setup(i) % num_mol *  &
          mol_setup(i) % molpar % natom
     setup % natom_type = setup % natom_type +  &
          mol_setup(i) % molpar % natom_type
  END DO

!..generate unique atom list
  ALLOCATE ( atom_names ( setup % natom_type ), STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'input_fist', &
       'atom_names', setup % natom_type )

  n = 0
  DO i = 1, nmol_type
     DO j = 1, mol_setup(i) % molpar % natom
        IF (str_search(atom_names,n,mol_setup(i) % molpar%aname(j))==0) THEN
           n = n + 1
           atom_names(n) = mol_setup(i) % molpar%aname(j)
        END IF
     END DO
  END DO

  IF ( n /= setup % natom_type ) CALL stop_program ( 'input_fist', &
       'natom_type mismatch' )

!..define pstat to hold all information (static) on different particle types
  ALLOCATE ( pstat ( setup % natom_type ), STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'input_fist', &
       'pstat', setup % natom_type )
  natom_types = SIZE ( atom_names )

  DO i = 1, nmol_type
     DO j = 1, mol_setup(i) % molpar % natom
        n = str_search(atom_names,natom_types,mol_setup(i) % molpar % aname(j))
        pstat(n) % charge = mol_setup(i) % molpar % acharge(j)
        pstat(n) % mass = mol_setup(i) % molpar % aweight(j)
        pstat(n) % ptype = n
     END DO
  END DO

!..read in non-bonded interaction information
  natom_types = size(atom_names)
  ALLOCATE (potparm(natom_types,natom_types),STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory ( 'fist_input', &
       'potparm', natom_types ** 2 )
  CALL input_nonbonded(potparm,atom_names,set_fn)

!..read the bond parameters
  CALL input_bonds(intra_param % bond_param,atom_names,set_fn)

!..read the bend parameters
  CALL input_bends(intra_param % bend_param,atom_names,set_fn)

!..read the torsion parameters
  CALL input_torsions(intra_param % torsion_param,atom_names,set_fn)

!..read the distance constraint parameters
  CALL input_distance_constraints ( intra_param % constraint_distance, &
       atom_names, set_fn )

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

SUBROUTINE input_nonbonded ( potparm, atom_names, set_fn )

  IMPLICIT NONE

! Arguments
  TYPE (potentialparm_type ), INTENT ( INOUT ) :: potparm ( :, : )
  CHARACTER ( LEN = * ), INTENT ( IN ) :: atom_names ( : )
  CHARACTER ( LEN = * ), INTENT ( IN ) :: set_fn

! Locals
  INTEGER :: ios, i, i1, i2, msglen, np, iw, source, allgrp
  INTEGER :: ierror, ilen, natom_types
  LOGICAL :: exists
  CHARACTER ( LEN = 10 ) :: at1, at2
  CHARACTER ( LEN = 40 ) :: string, string2
  CHARACTER ( LEN = 12 ) :: label
  REAL ( dbl ) :: epsilon, sigma, rcut, a, b, c, c0, c1, c2, d
  REAL ( dbl ) :: vr0, dc, m, mc

!------------------------------------------------------------------------------

  natom_types = SIZE ( atom_names )
  np = 0
  DO i1 = 1, natom_types
     DO i2 = 1, natom_types
        potparm(i1,i2) %type = 'NULL'
        CALL uppercase(potparm(i1,i2) %type )
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
  CALL parser_init(set_fn,label,ierror,globenv)
  IF (ierror /= 0 ) THEN
     IF (globenv%ionode) THEN
        WRITE ( iw, '( a )' ) ' No input section &FORCE_FIELD found on file '
        WRITE ( iw, '( a )' ) set_fn
     END IF
     CALL stop_parser( 'input_nonbond','FORCE_FIELD')
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
  CALL parser_end

!..write some information to output
  IF (globenv%ionode) THEN
     IF (globenv % print_level>=0 .AND. np>0) THEN
        WRITE ( iw, '( A )' ) ' FORCE FIELD| Non-bonded interaction parameters'
        DO i1 = 1, natom_types
           DO i2 = i1, natom_types
              IF (potparm(i1,i2)%type=='LENNARD-JONES') THEN
                 WRITE ( iw, '( A,A,T30,A,T37,A,T58,A,T71,F10.4 )' ) &
                      ' FORCE FIELD| ', 'LJ (12-6)', &
                      adjustl(atom_names(i1)), adjustl(atom_names(i2)), &
                      ' epsilon [K]', potparm(i1,i2) %lj%epsilon
                 WRITE ( iw, '( A,T58,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
                      '   sigma [A]', (potparm(i1,i2)%lj%sigma6)**(1./6.)
                 WRITE ( iw, '( A,T58,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
                      '    rcut [A]', sqrt(potparm(i1,i2)%rcutsq)
              ELSE IF (potparm(i1,i2)%type=='WILLIAMS') THEN
                 WRITE ( iw, '( A,A,T30,A,T37,A,T55,A,T68,E13.4 )' ) &
                      ' FORCE FIELD| ', 'Williams', adjustl(atom_names(i1)), &
                      adjustl(atom_names(i2)), 'A[K]', potparm(i1,i2) %willis%a
                 WRITE ( iw, '( A,T47,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
                      '    B [A^-1]', potparm(i1,i2) %willis%b
                 WRITE ( iw, '( A,T47,A,T68,E13.4 )' ) ' FORCE FIELD| ', &
                      '   C [K A^6]', potparm(i1,i2) %willis%c
                 WRITE ( iw, '( A,T47,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
                      '    rcut [A]', sqrt(potparm(i1,i2)%rcutsq)
              ELSE IF (potparm(i1,i2)%type=='GOODWIN') THEN
                 WRITE ( iw, '( A,A,T30,A,T37,A,T55,A,T68,E13.4 )' ) &
                      ' FORCE FIELD| ', 'Goodwin', adjustl(atom_names(i1)), &
                      adjustl(atom_names(i2)), 'Ve0[eV]', &
                      potparm(i1,i2) %goodwin%vr0
                 WRITE ( iw, '( A,T44,A,T61,2F10.4 )' ) ' FORCE FIELD| ', &
                      '  d,dc [A]', potparm(i1,i2) %goodwin%d, &
                      potparm(i1,i2) %goodwin%dc
                 WRITE ( iw, '( A,T44,A,T61,2F10.4 )' ) ' FORCE FIELD| ', &
                      '  m,mc   ', potparm(i1,i2) %goodwin%m, &
                      potparm(i1,i2) %goodwin%mc
                 WRITE ( iw, '( A,T44,A,T71,F10.4 )' ) ' FORCE FIELD| ', &
                      '  rcut [A]', sqrt(potparm(i1,i2)%rcutsq)
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
!!     Read all the bond parameters. Put them in the
!!     bond_parm matrix. Symmetrize this matrix since a bond between atoms
!!     A and B is the same as between B and A. While we read the bond_pot_type,
!!     we (so far) will do nothing with this information. It appears that
!!     all major parameters sets, such as CHARMM or AMBER, all have only
!!     one type of bonded interaction (harmonic).
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
!!I bonds                                                                     I
!!I    harmonic atom1 atom2 bp1 bp2                                           I
!!I end bonds                                                                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_bonds(bond_parm,atom_names,set_fn)
  IMPLICIT NONE

! Arguments
  TYPE (bond_parameters_type ), POINTER :: bond_parm ( :, : )
  CHARACTER ( LEN = * ), INTENT ( IN ) :: atom_names ( : )
  CHARACTER ( LEN = * ), INTENT ( IN ) :: set_fn

! Locals
  INTEGER :: i1, i2, ierror, ilen, ios
  INTEGER :: nb, iw, allgrp, source, natom_types
  LOGICAL :: exists
  CHARACTER ( LEN = 20 ) :: type, at1, at2
  CHARACTER ( LEN = 40 ) :: string, string2
  CHARACTER ( LEN = 12 ) :: label
  REAL ( dbl ) :: bp1, bp2

!------------------------------------------------------------------------------

  natom_types=SIZE(atom_names)
  nb = 0
  ALLOCATE (bond_parm(natom_types,natom_types),STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory ( 'input_bonds', &
       'bond_parm', natom_types ** 2 )

  bond_parm ( :, : ) %type = 'null'
  bond_parm ( :, : ) %r0 = 0._dbl
  bond_parm ( :, : ) %k = 0._dbl
  iw = scr
!..parse the input section
  label = '&FORCE_FIELD'
  CALL parser_init(set_fn,label,ierror,globenv)
  IF (ierror /= 0 ) THEN
     IF ( globenv % ionode ) THEN
        WRITE ( iw, '( a )' ) ' No input section &FORCE_FIELD found on file '
        WRITE ( iw, '( a )' ) set_fn
     END IF
     CALL stop_parser( 'input_bond', '&FORCE_FIELD' )
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 5
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        IF (index(string,'BONDS') /= 0 ) THEN
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
                 nb = nb + 1
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( bp1 )
                 CALL get_next ( bp2 )
                 i1 = str_search(atom_names,natom_types,at1)
                 IF (i1==0) CALL stop_parser( 'input_bonds','atom not found')
                 i2 = str_search(atom_names,natom_types,at2)
                 IF (i2==0) CALL stop_parser( 'input_bonds','atom not found')
                 bond_parm(i1,i2) %type = string2
                 bond_parm(i1,i2) %r0 = bp1
                 bond_parm(i1,i2) %k = bp2
                 bond_parm(i2,i1) = bond_parm(i1,i2)
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
  CALL parser_end

!..write some information to output
  IF ( globenv % ionode ) THEN
     IF (globenv % print_level>=0 .AND. nb>0) THEN
        WRITE ( iw, '( A )' ) ' FORCE FIELD| Bond interaction parameters'
        DO i1 = 1, natom_types
           DO i2 = i1, natom_types
              IF (index(bond_parm(i1,i2)%type,'HARMONIC') /= 0 ) THEN
                 WRITE ( iw, '( A, T30, A, T37, A, T58, A, T71, F10.2 )' ) &
                      ' FORCE FIELD| harmonic', &
                      ADJUSTL ( atom_names ( i1 ) ), ADJUSTL ( atom_names ( i2 ) ), &
                      '      r0 [A]', bond_parm ( i1, i2 ) % r0
                 WRITE ( iw, '( A, T58, A, T71, F10.2 )' ) ' FORCE FIELD| ', &
                      ' k [K*A^-2)]', bond_parm ( i1, i2 ) % k
              END IF
           END DO
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
!!     Read all the bend parameters. Put them in the
!!     bend_parm matrix. Fill this matrix since a bend between atoms
!!     A-B-C is the same as between C-B-A. While we read the bend_pot_type,
!!     we (so far) will do nothing with this information. It appears that
!!     all major parameters sets, such as CHARMM or AMBER, all have only
!!     only one type of bend interaction (harmonic).
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
!!I bends                                                                     I
!!I    harmonic atom1 atom2 atom3 bp1 bp2                                     I
!!I end bends                                                                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_bends ( bend_parm, atom_names, set_fn )

  IMPLICIT NONE

! Arguments
  TYPE (bend_parameters_type ), POINTER  :: bend_parm(:,:,:)
  CHARACTER ( LEN = * ), INTENT ( IN ) :: atom_names ( : )
  CHARACTER ( LEN = * ), INTENT ( IN ) :: set_fn

! Locals
  INTEGER :: i1, i2, i3, ierror, ilen, ios
  INTEGER ::  nb, iw, allgrp, source, natom_types
  CHARACTER ( LEN = 20 ) :: type, at1, at2, at3
  CHARACTER ( LEN = 40 ) :: string, string2
  CHARACTER ( LEN = 12 ) :: label
  REAL ( dbl ) :: bp1, bp2, torad, todeg

!------------------------------------------------------------------------------

  natom_types=SIZE(atom_names)
  nb = 0
  torad = pi / 180.0_dbl
  todeg = 180.0_dbl/pi
  ALLOCATE (bend_parm(natom_types,natom_types,natom_types),STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory ( 'input_bends', &
       'bend_parm', natom_types ** 3 )
  bend_parm(:,:,:) %type = 'null'
  bend_parm(:,:,:) %theta0 = 0._dbl
  bend_parm(:,:,:) %k = 0._dbl

  iw = scr
!..parse the input section
  label = '&FORCE_FIELD'
  CALL parser_init(set_fn,label,ierror,globenv)
  IF (ierror /= 0 ) THEN
     IF ( globenv % ionode ) THEN
        WRITE ( iw, '( a )' ) ' No input section &FORCE_FIELD found on file '
        WRITE ( iw, '( a )' ) set_fn
     END IF
     CALL stop_parser( 'input_bend', '&FORCE_FIELD' )
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 5
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        IF (index(string,'BENDS') /= 0 ) THEN
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
                 nb = nb + 1
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 ilen = 0
                 CALL get_next(at3,ilen)
                 CALL get_next ( bp1 )
                 CALL get_next ( bp2 )
                 i1 = str_search(atom_names,natom_types,at1)
                 IF (i1==0) CALL stop_parser( 'input_bends','atom not found')
                 i2 = str_search(atom_names,natom_types,at2)
                 IF (i2==0) CALL stop_parser( 'input_bends','atom not found')
                 i3 = str_search(atom_names,natom_types,at3)
                 IF (i3==0) CALL stop_parser( 'input_bends','atom not found')
                 bend_parm(i1,i2,i3) %type = string2
                 bend_parm(i1,i2,i3) %theta0 = bp1*torad
                 bend_parm(i1,i2,i3) %k = bp2
                 bend_parm(i3,i2,i1) = bend_parm(i1,i2,i3)
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
  CALL parser_end

!..write some information to output
  IF ( globenv % ionode ) THEN
     IF (globenv % print_level>=0 .AND. nb>0) THEN
        WRITE ( iw, '( A )' ) ' FORCE FIELD| Bend interaction parameters'
        DO i1 = 1, natom_types
           DO i2 = 1, natom_types
              DO i3 = i1, natom_types
                 IF (index(bend_parm(i1,i2,i3)%type,'HARMONIC') /= 0 ) THEN
                    WRITE ( iw, '( A,T30,A,T37,A,T44,A,T58,A,T71,F10.2 )' ) &
                         ' FORCE FIELD| harmonic', &
                         adjustl(atom_names(i1)), adjustl(atom_names(i2)), adjustl(atom_names(i3)), &
                         'theta0 [deg]', bend_parm(i1,i2,i3) %theta0*todeg
                    WRITE ( iw, '( A,T58,A,T71,F10.2 )' ) ' FORCE FIELD| ', &
                         '       k [K]', bend_parm(i1,i2,i3) %k
                 END IF
              END DO
           END DO
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
!!     Read all the torsion parameters. Put them in the
!!     torsion_parm matrix. Fill this matrix since a torsion between atoms
!!     A-B-C-D is the same as between D-C-B-A. We have defined the
!!     torsion_parameter type with 7 coefficients for 0 to 6. Hopefully
!!     this will be generic enough to handle all torsion potentials.
!!     Since different types will be used, we need to store the type of torsion
!!     as well.
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
!!I torsions                                                                  I
!!I    jorgensen atom1 atom2 atom3 atom4 tp[1-4]                              I
!!I end torsions                                                              I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_torsions ( torsion_parm, atom_names, set_fn )

  IMPLICIT NONE

! Arguments
  TYPE ( torsion_parameters_type ), POINTER :: torsion_parm ( :, :, :, : )
  CHARACTER ( LEN = * ), INTENT ( IN ) :: atom_names ( : )
  CHARACTER ( LEN = * ), INTENT ( IN ) :: set_fn

! Locals
  INTEGER :: i, i1, i2, i3, i4, ierror, ilen, ios
  INTEGER :: nb, iw, allgrp, source, natom_types
  CHARACTER ( LEN = 20 ) :: type, at1, at2, at3, at4
  CHARACTER ( LEN = 40 ) :: string, string2
  CHARACTER ( LEN = 12 ) :: label
  REAL ( dbl ), DIMENSION (0:3) :: tp

!------------------------------------------------------------------------------

  natom_types=SIZE(atom_names)
  nb = 0
  ALLOCATE ( &
       torsion_parm ( natom_types, natom_types, natom_types, natom_types ), &
       STAT = ios )
  IF ( ios /= 0 ) CALL stop_memory ( 'input_torsions', &
       'torsion_parm', natom_types ** 4 )
  torsion_parm ( :, :, :, : ) % type = 'null'
  DO i = 0, SIZE ( torsion_parm(1,1,1,1) % a ) - 1
     torsion_parm(:,:,:,:) %a(i) = 0._dbl
  END DO

  iw = scr

!..parse the input section
  label = '&FORCE_FIELD'
  CALL parser_init(set_fn,label,ierror,globenv)
  IF ( ierror /= 0 ) THEN
     IF ( globenv % ionode ) THEN
        WRITE ( iw, '( A )' ) &
             ' No input section &FORCE_FIELD found on file '
        WRITE ( iw, '( A )' ) set_fn
     END IF
     CALL stop_parser( 'input_torsion', '&FORCE_FIELD' )

  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 8
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        IF (index(string,'TORSIONS') /= 0 ) THEN
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
                 nb = nb + 1
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 ilen = 0
                 CALL get_next(at3,ilen)
                 ilen = 0
                 CALL get_next(at4,ilen)
                 DO i = 0, SIZE(tp) - 1
                    CALL get_next ( tp(i) )
                 END DO
                 i1 = str_search(atom_names,natom_types,at1)
                 IF (i1==0) CALL stop_parser( 'input_torsions', &
                      'atom not found')
                 i2 = str_search(atom_names,natom_types,at2)
                 IF (i2==0) CALL stop_parser( 'input_torsions', &
                      'atom not found')
                 i3 = str_search(atom_names,natom_types,at3)
                 IF (i3==0) CALL stop_parser( 'input_torsions', &
                      'atom not found')
                 i4 = str_search(atom_names,natom_types,at4)
                 IF (i4==0) CALL stop_parser( 'input_torsions', &
                      'atom not found')
                 torsion_parm(i1,i2,i3,i4) %type = string2
                 DO i = 1, SIZE(tp) - 1
                 torsion_parm(i1,i2,i3,i4) %a(i) = tp(i)
                 END DO
                 torsion_parm(i4,i3,i2,i1) = torsion_parm(i1,i2,i3,i4)
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
  CALL parser_end()

!..write some information to output
  IF ( globenv % ionode ) THEN
     IF ( globenv % print_level >= 0 .AND. nb > 0 ) THEN
        WRITE ( iw, '( A )' ) &
             ' FORCE FIELD| Torsion interaction parameters'
        DO i1 = 1, natom_types
           DO i2 = 1, natom_types
              DO i3 = i2, natom_types
                 DO i4 = i1, natom_types
                    IF ( index ( torsion_parm(i1,i2,i3,i4) % type, &
                         'JORGENSEN' ) /= 0 ) THEN
                       WRITE ( iw, '( A, T30, A, T37, A, T44, A, T51, A, T58, A, T71, F10.2 )' ) &
                            ' FORCE FIELD| Jorgensen', &
                            ADJUSTL ( atom_names(i1)), &
                            ADJUSTL ( atom_names(i2)), &
                            ADJUSTL ( atom_names(i3)), &
                            ADJUSTL ( atom_names(i4)), &
                            '      a0 [K]', torsion_parm ( i1, i2, i3, i4 ) % a(0)
                       WRITE ( iw, '( A, T58, A, T71, F10.2 )' ) ' FORCE FIELD| ', &
                            '      a1 [K]', torsion_parm ( i1, i2, i3, i4 ) % a(1)
                       WRITE ( iw, '( A, T58, A, T71, F10.2 )' ) ' FORCE FIELD| ', &
                            '      a2 [K]', torsion_parm ( i1, i2, i3, i4 ) % a(2)
                       WRITE ( iw, '( A, T58, A, T71, F10.2 )' ) ' FORCE FIELD| ', &
                            '      a3 [K]', torsion_parm ( i1, i2, i3, i4 ) % a(3)
                    END IF
                 END DO
              END DO
           END DO
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
!!     constraint_distance array. Symmetrize this matrix since a bond between
!!     atoms A and B is the same as between B and A.
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
!!I distcon                                                                   I
!!I     atom1 atom2 dp                                                        I
!!I end distcon                                                               I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_distance_constraints ( constraint_distance, atom_names, &
     set_fn )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), POINTER :: constraint_distance ( :, : )
  CHARACTER ( LEN = * ), INTENT ( IN ) :: atom_names ( : )
  CHARACTER ( LEN = * ), INTENT ( IN ) :: set_fn

! Locals
  INTEGER :: i1, i2, ierror, ilen, ios
  INTEGER :: nb, iw, allgrp, source, natom_types
  CHARACTER ( LEN = 20 ) :: at1, at2
  CHARACTER ( LEN = 40 ) :: string, string2, string2_tmp
  CHARACTER ( LEN = 12 ) :: label
  REAL ( dbl ) :: dp

!------------------------------------------------------------------------------

  natom_types=SIZE(atom_names)
  nb = 0
  ALLOCATE (constraint_distance(natom_types,natom_types),STAT=ios)
  IF ( ios /= 0 ) CALL stop_memory ( 'input_distance_constraints', &
       'constraint_distance', natom_types ** 2 )
  constraint_distance ( :, : ) = 0.0_dbl

  iw = scr

!..parse the input section
  label = '&FORCE_FIELD'
  CALL parser_init(set_fn,label,ierror,globenv)
  IF (ierror /= 0 ) THEN
     IF ( globenv % ionode ) THEN
        WRITE ( iw, '( A )' ) ' No input section &FORCE_FIELD found on file '
        WRITE ( iw, '( A )' ) set_fn
     END IF
     CALL stop_parser ( 'input_distance_constraints', '&FORCE_FIELD' )
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 7
        CALL get_next ( string, ilen )
        CALL uppercase ( string )

        IF ( INDEX ( string, 'DISTCON' ) /= 0 ) THEN
           DO
              CALL read_line()

              ilen = 0
              CALL get_next ( string2, ilen )
              string2_tmp = string2
              CALL uppercase ( string2_tmp )
              IF ( INDEX ( string2_tmp, 'END' ) /= 0 ) EXIT

              at1 = string2
              nb = nb + 1
              ilen = 0
              CALL get_next ( at2, ilen )
              i1 = str_search ( atom_names, natom_types, at1 )
              IF ( i1 == 0 ) CALL stop_parser ( 'input_distance_constraints', &
                   'atom not found' )
              i2 = str_search ( atom_names, natom_types, at2 )
              IF ( i2 == 0 ) CALL stop_parser ( 'input_distance_constraints', &
                   'atom not found' )
              CALL get_next ( dp )
              constraint_distance ( i1, i2 ) = dp
              constraint_distance ( i2, i1 ) = dp
           END DO
        END IF
        CALL read_line
     END DO
  END IF
  CALL parser_end

!..write some information to output
  IF ( globenv % ionode ) THEN
     IF ( globenv % print_level >= 0 .AND. nb > 0 ) THEN
        WRITE ( iw, '( A )' ) ' FORCE FIELD| Distance constraints '
        DO i1 = 1, natom_types
           DO i2 = i1, natom_types
              IF ( constraint_distance ( i1, i2 ) /= 0.0_dbl ) THEN
                 WRITE ( iw, '( A, T30, A, T37, A, T71, F10.4 )' ) &
                      ' FORCE FIELD| Constraint: ', &
                      ADJUSTL ( atom_names ( i1 ) ), &
                      ADJUSTL ( atom_names ( i2 ) ), &
                      constraint_distance ( i1, i2 )
              END IF
           END DO
        END DO
        WRITE ( iw, '()' )
     END IF
  END IF

END SUBROUTINE input_distance_constraints

!!*****
!******************************************************************************

END MODULE force_fields

!******************************************************************************
