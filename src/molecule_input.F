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
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE molecule_input

  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi, zero
  USE molecule_types, ONLY : molecule_type
  USE message_passing, ONLY : mp_bcast
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       cfield, p_error, get_real, get_int, stop_parser
  USE string_utilities, ONLY : uppercase, xstring, str_search, &
       str_comp, make_tuple
  USE termination, ONLY : stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_setup_section, read_molecule_section, charge, multiplicity

! total charge of the system

  REAL ( dbl ) :: charge

! multiplicity of the system

  INTEGER :: multiplicity

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
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &setup ... &end                                                  I
!!I                                                                           I
!!I charge           charge                                                   I
!!I multiplicity     multiplicity                                             I
!!I moltype          nmol_types                                               I
!!I molecule                                                                  I
!!I    name          molname                                                  I
!!I    nummol        nmol                                                     I
!!I    nhcopt        [molecule,global,massive,...]                            I
!!I end molecule                                                              I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_setup_section ( mol_setup, set_fn, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
  CHARACTER ( LEN = * ), INTENT ( IN ) :: set_fn
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: ierror, ilen, iw, source, group, icount, nmol_type, ios, i
  CHARACTER ( LEN = 20 ) :: string, string2
  CHARACTER ( LEN = 6 ) :: label

!------------------------------------------------------------------------------

  charge = 0.0_dbl
  multiplicity = 0
  iw = globenv % scr
  icount = 0

!..parse the input section
  label = '&SETUP'
  CALL parser_init ( set_fn, label, ierror, globenv )
  IF ( ierror /= 0 ) THEN
     IF ( globenv % ionode ) THEN
        WRITE ( iw, '( A )' ) ' No input section &SETUP found on file '
        WRITE ( iw, '( T2, A )' ) set_fn
     END IF
     CALL stop_parser ( 'read_setup_section', '&SETUP' )
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 8
        CALL cfield ( string, ilen )
        CALL uppercase ( string )
        SELECT CASE (string)
        CASE DEFAULT
           CALL p_error()
           CALL stop_parser ( 'read_setup_section','unknown option')
        CASE ( 'CHARGE')
           charge = get_real()
        CASE ( 'MULTIPLI')
           multiplicity = get_int()
        CASE ( 'MOLTYPE')
           nmol_type = get_int()
           ALLOCATE (mol_setup(nmol_type ),STAT=ios)
           IF ( ios /= 0 ) CALL stop_memory ( 'read_setup_section', &
                'mol_setup', nmol_type )
        CASE ( 'MOLECULE')
           IF ( .NOT. ASSOCIATED ( mol_setup ) ) &
                CALL stop_program ( 'read_setup_section', &
                'number of molecule types not yet specified')
           icount = icount + 1
           IF (icount>nmol_type ) CALL stop_program ( 'read_setup_section', &
                'too many molecule types specified')
           DO
              CALL read_line
              ilen = 6
              CALL cfield ( string2, ilen )
              CALL uppercase ( string2 )
              SELECT CASE ( string2 )
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_setup_section','unknown option')
              CASE ( 'NAME')
                 ilen = 0
                 CALL cfield(mol_setup(icount) % molpar % name,ilen)
              CASE ( 'NUMMOL')
                 mol_setup(icount) % num_mol = get_int()
              CASE ( 'NHCOPT')
                 ilen = 0
                 CALL cfield(mol_setup(icount) % nhcopt,ilen)
                 CALL uppercase(mol_setup(icount) % nhcopt)
              CASE ( 'END')
                 ilen = 8
                 CALL cfield(string2,ilen)
                 EXIT
              END SELECT
           END DO
        END SELECT
        CALL read_line
     END DO
  END IF
  CALL parser_end

!..write some information to output
  IF ( globenv % ionode ) THEN
     IF (globenv%print_level>=0) THEN
        WRITE ( iw, '( A,T71,F10.4 )' ) ' MOLECULE| Total charge of system ', &
             charge
        WRITE ( iw, '( A,T71,I10 )' ) ' MOLECULE| Multiplicity of system ', &
             multiplicity
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
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &molecule ... &end                                               I
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

SUBROUTINE read_molecule_section ( mol_setup, set_fn, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( INOUT ) :: mol_setup
  CHARACTER ( LEN = * ), INTENT ( IN ) :: set_fn
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: ierror, ilen, iw, source, group, icount, nmol_type, ios
  INTEGER :: i, j, imol, na, nat, j1, j2, nf, nt, ie, ia, i1, i2, i3, i4
  INTEGER :: isos, nbasis, ndrho, allgrp, npol, nrho0
  CHARACTER ( LEN = 20 ) :: libstring, string, string2, at1, at2, at3, at4
  CHARACTER ( LEN = 6 ) :: label
  CHARACTER ( LEN = 20 ), DIMENSION ( : ), ALLOCATABLE :: molname
  CHARACTER ( LEN = 20 ), DIMENSION ( : ), ALLOCATABLE :: tuple
  CHARACTER ( LEN = 80 ) :: bl = ' '

!------------------------------------------------------------------------------

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
  CALL parser_init(set_fn,label,ierror,globenv)
  IF (ierror /= 0 ) THEN
     IF ( globenv % ionode ) THEN
        WRITE ( iw, '( A )' ) ' No input section &MOLECULE found on file '
        WRITE ( iw, '( T2, A )' ) set_fn
     END IF
     CALL stop_parser ( 'read_molecule_section', '&MOLEC' )
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 8
        CALL cfield(string,ilen)
        imol = str_search(molname,nmol_type,string)
        nat = -1
        IF (imol>0) THEN
           mol_setup(imol) %molpar%natom_type = 0
           mol_setup(imol) %molpar%nbonds = 0
           mol_setup(imol) %molpar%nbends = 0
           mol_setup(imol) %molpar%ntorsions = 0
           mol_setup(imol) %molpar%ndcon = 0
           mol_setup(imol) %molpar%n3x3con = 0
           mol_setup(imol) %molpar%nbasis_pol = 0
           mol_setup(imol) %molpar%nbasis_rho0 = 0

           DO
              CALL read_line
              ilen = 8
              CALL cfield(string2,ilen)
              CALL uppercase ( string2 )

              SELECT CASE (string2)
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_molecule_section','unknown option')

              CASE ( 'NUMAT')
                 na = get_int()
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
                 CALL cfield(mol_setup(imol)%molpar%aname(nat),ilen)
                 ilen = 0
                 CALL cfield(mol_setup(imol)%molpar%alabel(nat),ilen)
                 mol_setup(imol) %molpar%aweight(nat) = get_real()
                 mol_setup(imol) %molpar%acharge(nat) = get_real()

              CASE ( 'BONDS')
                 mol_setup(imol) %molpar%nbonds = get_int()
                 na = mol_setup(imol) %molpar%nbonds
                 IF (na>=0) THEN
                    ALLOCATE (mol_setup(imol)%molpar%bonds(2,na),STAT=ios)
                    IF (ios /= 0 ) CALL stop_memory &
                         ( 'read_molecule_section', 'bonds', 2 * na )
                    DO i = 1, na
                       mol_setup(imol) %molpar%bonds(1,i) = get_int()
                       mol_setup(imol) %molpar%bonds(2,i) = get_int()
                    END DO
                 END IF

              CASE ( 'BENDS')
                 mol_setup(imol) %molpar%nbends = get_int()
                 na = mol_setup(imol) %molpar%nbends
                 IF (na>=0) THEN
                    ALLOCATE (mol_setup(imol)%molpar%bends(3,na),STAT=ios)
                    IF (ios /= 0 ) CALL stop_memory &
                         ( 'read_molecule_section', 'bends', 3 * na )
                    DO i = 1, na
                       mol_setup(imol) %molpar%bends(1,i) = get_int()
                       mol_setup(imol) %molpar%bends(2,i) = get_int()
                       mol_setup(imol) %molpar%bends(3,i) = get_int()
                    END DO
                 END IF

              CASE ( 'TORSIO')
                 mol_setup(imol) %molpar%ntorsions = get_int()
                 na = mol_setup(imol) %molpar%ntorsions
                 IF (na>=0) THEN
                    ALLOCATE ( &
                         mol_setup ( imol ) % molpar % torsions ( 4, na ), &
                         STAT = ios )
                    IF ( ios /= 0 ) CALL stop_memory &
                         ( 'read_molecule_section', 'torsions', 4 * na )
                    DO i = 1, na
                       mol_setup(imol) %molpar%torsions(1,i) = get_int()
                       mol_setup(imol) %molpar%torsions(2,i) = get_int()
                       mol_setup(imol) %molpar%torsions(3,i) = get_int()
                       mol_setup(imol) %molpar%torsions(4,i) = get_int()
                    END DO
                 END IF

              CASE ( 'DISTCO')
                 mol_setup(imol) %molpar%ndcon = get_int()
                 na = mol_setup(imol) %molpar%ndcon
                 IF (na>=0) THEN
                    ALLOCATE (mol_setup(imol)%molpar%ndc(2,na),STAT=ios)
                    IF ( ios /= 0 ) CALL stop_memory &
                         ( 'read_molecule_section', 'ndc', 2 * na )
                    DO i = 1, na
                       mol_setup(imol) %molpar%ndc(1,i) = get_int()
                       mol_setup(imol) %molpar%ndc(2,i) = get_int()
                    END DO
                 END IF

              CASE ( '3X3CON')
                 mol_setup(imol) %molpar%n3x3con = get_int()
                 na = mol_setup(imol) %molpar%n3x3con
                 IF (na>=0) THEN
                    ALLOCATE (mol_setup(imol)%molpar%n3x3c(3,na),STAT=ios)
                    IF (ios /= 0 ) CALL stop_memory &
                         ( 'read_molecule_section', 'n3x3c', 3 * na )
                    DO i = 1, na
                       mol_setup(imol) %molpar%n3x3c(1,i) = get_int()
                       mol_setup(imol) %molpar%n3x3c(2,i) = get_int()
                       mol_setup(imol) %molpar%n3x3c(3,i) = get_int()
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
                 CALL cfield(mol_setup(imol)%molpar%rho0_basis(nrho0)%aname,ilen)
                 ilen=0
                 CALL cfield(mol_setup(imol)%molpar%rho0_basis(nrho0)%libname,ilen)
                 mol_setup(imol)%molpar%rho0_basis(nrho0)%nmax=get_int()
                 mol_setup(imol)%molpar%rho0_basis(nrho0)%eps=get_real()

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
                 CALL cfield(mol_setup(imol)%molpar%pol_basis(npol)%aname,ilen)
                 ilen=0
                 CALL cfield(mol_setup(imol)%molpar%pol_basis(npol)%libname,ilen)
                 mol_setup(imol)%molpar%pol_basis(npol)%nmax=get_int()
                 mol_setup(imol)%molpar%pol_basis(npol)%eps=get_real()

              CASE ( 'END')
                 ilen = 0
                 CALL cfield(string2,ilen)
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
        ELSE
           WRITE ( iw, '( A )' ) ' molecule type not found '
           WRITE ( iw, '( T2, A )' ) string
           CALL stop_parser ( 'read_molecule_section', ' ' )
        END IF
        CALL read_line
     END DO
  END IF
  CALL parser_end

  DEALLOCATE (molname,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'read_molecule_section', 'molname')

! ..write some information to output
  IF (globenv%ionode) THEN
     IF (globenv%print_level>=0) THEN
        DO imol = 1, nmol_type
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

           WRITE ( iw,'( )' )
        END DO
     END IF
  END IF

END SUBROUTINE read_molecule_section

!!*****
!******************************************************************************

END MODULE molecule_input

!******************************************************************************
