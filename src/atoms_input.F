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
!!   gt 23-09-2002 : major changes. Pointer part is allocated/deallocated
!!                   and initialized here. Atomic coordinates can now be 
!!                   read also from &COORD section in the input file.
!!                   If &COORD is not found, .dat file is read.  
!!                   If & coord is found and .NOT. 'INIT', parsing of the .dat
!!                   is performed to get the proper coords/vel/eta variables
!!    CJM 31-7-03  : Major rewrite.  No more atype
!!
!!   SOURCE
!******************************************************************************

MODULE atoms_input

  USE atomic_kind_types,        ONLY : allocate_atomic_kind_set,&
                                       atomic_kind_type,&
                                       get_atomic_kind,&
                                       init_atomic_kind_set,&
                                       set_atomic_kind
  USE convert_units,            ONLY : convert_to_cp2k_units
  USE external_potential_types, ONLY : epc_potential_type,&
                                       set_potential
  USE global_types,             ONLY : global_environment_type
  USE input_types,              ONLY : setup_parameters_type
  USE kinds,                    ONLY : dbl
  USE memory_utilities,         ONLY : reallocate
  USE molecule_types,           ONLY : molecule_type
  USE extended_system_dynamics, ONLY : lnhc_parameters_type, &
                                       npt_info_type
  USE particle_types,           ONLY : particle_type, &
                                       allocate_particle_set
  USE parser,                   ONLY : parser_init, &
                                       parser_end, &
                                       read_line, &
                                       test_next, &
                                       p_error, &
                                       get_next, &
                                       search_label
  USE periodic_table,           ONLY : ptable
  USE simulation_cell,          ONLY : cell_type
  USE string_utilities,         ONLY : xstring,&
                                       integer_to_string, &
                                       string_to_integer, &
                                       str_search
  USE termination,              ONLY : stop_memory, &
                                       stop_program
  USE topology_control,         ONLY : read_coordinate_control
  USE util,                     ONLY : get_unit

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_atoms_input, read_md_restart, read_pimd_atoms_input

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** atoms_input/read_atoms_input [1.0] *
!!
!!   NAME
!!     read_atoms_input
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    CJM
!!
!!   MODIFICATION HISTORY
!!
!!  OPTIONS
!!
!!  INPUTS
!!*** *************************************************************************

SUBROUTINE read_atoms_input ( unit_type, particle_set, atomic_kind_set, &
                              mol_setup, setup, globenv )

  IMPLICIT NONE

! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: unit_type
  TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE ( molecule_type ),DIMENSION(:), INTENT(IN) :: mol_setup
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE ( setup_parameters_type ) :: setup



! Locals
  INTEGER :: ierror, ilen, iw,  ia, ie, i
  INTEGER, PARAMETER :: max_name_length = 60
  INTEGER, PARAMETER :: nblock_atom = 1000,&
                        nblock_kind = 10

  TYPE(atomic_kind_type), POINTER :: atomic_kind
  TYPE(epc_potential_type), POINTER :: epc_potential
  CHARACTER(LEN=max_name_length)  :: current_label,label,string
  CHARACTER(LEN=200)              :: message
  CHARACTER(LEN=8)                :: unit_name
  CHARACTER(LEN=3)                :: test_result
  INTEGER                         :: iatom,ikind,istat,natom,nkind,z,n,iat,&
                                     natom_types,nmol_type,k,j
  LOGICAL                         :: new_kind,coords,&
                                       init

  INTEGER, DIMENSION(:), POINTER :: kind_of,natom_of_kind,natom_of_kind_check
  INTEGER, DIMENSION(:,:), POINTER  :: atom_list,kind_label
  REAL(dbl), DIMENSION (:), POINTER :: charge, mass
  REAL(dbl), DIMENSION(:,:), POINTER :: r
  CHARACTER ( len = 20 ), DIMENSION ( : ), POINTER :: atom_names

!   ---------------------------------------------------------------------------

  nkind = 0
  natom = 0
  iat = 0

  iw = globenv % scr

  atom_names => setup % atom_names
!..defaults
  NULLIFY (atom_list,kind_label,kind_of,natom_of_kind,charge,mass,r)

  CALL reallocate(atom_list,1,nblock_kind,1,nblock_atom)
  CALL reallocate(kind_label,1,max_name_length,1,nblock_atom)
  CALL reallocate(kind_of,1,nblock_atom)
  CALL reallocate(natom_of_kind,1,nblock_kind)
  CALL reallocate(mass,1,nblock_kind)
  CALL reallocate(charge,1,nblock_kind)
  CALL reallocate(r,1,3,1,nblock_atom)

!..parse the input section
  label = '&COORD'
  CALL parser_init ( globenv % input_file_name, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
     WRITE(iw,*) "Calling coorindate_control"
     IF (setup%coordinate) THEN
       CALL read_coordinate_control(particle_set,atomic_kind_set,&
                                    mol_setup,setup,globenv)
     ELSE IF( globenv % ionode) THEN 
       CALL stop_program ( 'read_atoms_input','No &COORD section found')
     END IF
  ELSE
!   ******************************************************************
!   **** Reads cordinates from &COORD section ***  
!   ******************************************************************
     CALL read_line
     DO WHILE (test_next()/='X')
       ilen = 8
       IF (test_next()=="C") THEN
          CALL get_next ( current_label, ilen )
       ELSE IF (test_next()=="N") THEN
          CALL get_next ( z )
          current_label = ptable ( z ) % symbol
       END IF

       natom = natom + 1

       IF (natom > SIZE(atom_list,2)) THEN
        CALL reallocate(atom_list,1,nkind,1,natom+nblock_atom)
        CALL reallocate(kind_of,1,natom+nblock_atom)
        CALL reallocate(r,1,3,1,natom+nblock_atom)
       END IF

!     *** Read the atomic coordinates ***

       CALL get_next(r(1,natom))
       CALL get_next(r(2,natom))
       CALL get_next(r(3,natom))
       label = current_label

       new_kind = .TRUE.

       DO ikind=1,nkind
         CALL integer_to_string(kind_label(:,ikind),string)
         IF (string == label) THEN
           new_kind = .FALSE.
           EXIT
         END IF
       END DO

       IF (new_kind) THEN
         nkind = nkind + 1
         IF (nkind > SIZE(atom_list,1)) THEN
           CALL reallocate(atom_list,1,nkind+nblock_kind,1,natom)
           CALL reallocate(kind_label,1,max_name_length,1,nkind+nblock_kind)
           CALL reallocate(natom_of_kind,1,nkind+nblock_kind)
           CALL reallocate(charge,1,nkind+nblock_kind)
           CALL reallocate(mass,1,nkind+nblock_kind)
         END IF
         ikind = nkind
         CALL string_to_integer(current_label,kind_label(:,ikind))
       END IF

       kind_of(natom) = ikind
       natom_of_kind(ikind) = natom_of_kind(ikind) + 1
       atom_list(ikind,natom_of_kind(ikind)) = natom

     call read_line
     END DO

    CALL parser_end
  END IF
!..end of parsing the input section

  IF (.NOT. setup%coordinate) THEN
    CALL reallocate(atom_list,1,nkind,1,natom)
    CALL reallocate(kind_label,1,max_name_length,1,nkind)
    CALL reallocate(kind_of,1,natom)
    CALL reallocate(natom_of_kind,1,nkind)
    CALL reallocate(charge,1,nkind)
    CALL reallocate(mass,1,nkind)
    CALL reallocate(r,1,3,1,natom)

!   *** Allocate the data structure for the atomic kind information ***

    CALL allocate_atomic_kind_set(atomic_kind_set,nkind)

!   *** Initialize the data structure for atomic kind information ***

    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL integer_to_string(kind_label(:,ikind),string)
      CALL set_atomic_kind(atomic_kind=atomic_kind,&
                           kind_number=ikind,&
                           name=string,&
                           atom_list=atom_list(ikind,1:natom_of_kind(ikind)))
    END DO
!   *** Allocate the data structure for the atomic information ***

    CALL allocate_particle_set(particle_set,natom)

!   *** Initialize the data structure for the atomic information ***
    DO iatom=1,natom
      ikind = kind_of(iatom)
      particle_set(iatom) % abase = 0
      particle_set(iatom) % grid_id = 0
      particle_set(iatom)%atomic_kind => atomic_kind_set(ikind)
      particle_set(iatom)%iatom = iatom
      particle_set(iatom)%f(:) = 0.0_dbl
      particle_set(iatom)%r(:) = r(:,iatom)
      particle_set(iatom)%v(:) = 0.0_dbl
      NULLIFY ( particle_set ( iatom ) % coef_list )
    END DO

!   *** check if the labels in &COORD match with the &MOLECULE setup
!   *** and store the charge and mass read in the &MOLECULE section

    ALLOCATE (natom_of_kind_check(nkind), STAT=istat) 
    natom_of_kind_check = 0
    nmol_type = SIZE (mol_setup) 
    natom_types = SIZE ( atom_names )
    iat = 0
    DO i = 1, nmol_type
      DO k = 1, mol_setup(i) % num_mol
        DO j = 1, mol_setup(i) % molpar % natom
          iat = iat + 1
          n = str_search ( atom_names, natom_types, &
             mol_setup ( i ) % molpar % aname ( j ) )
          charge(n)= mol_setup(i) % molpar % acharge(j)
! Convert mass from amu -> atomic units with massunit
          mass(n)= mol_setup(i) % molpar % aweight(j)
          CALL get_atomic_kind(particle_set(iat)%atomic_kind,&
                               kind_number=ikind)
          IF (ikind/=n) CALL stop_program ( 'read_atoms_input', &
             'kind label ordering in &COORD and &MOL dont match ')
          CALL integer_to_string(kind_label(:,ikind),string)
          IF (string/=mol_setup(i)%molpar%aname(j)) CALL stop_program ( 'read_atoms_input', &
              'inconsistent kind labels in &COORD and &MOL - check input')
          natom_of_kind_check(ikind) = natom_of_kind_check(ikind) + 1
          iatom =  atom_list(ikind,natom_of_kind_check(ikind)) 
          IF (iatom/=iat) CALL stop_program ( 'read_atoms_input', &
             'atom lists in &COORD and &MOL dont match ')
        END DO
      END DO
    END DO
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL set_atomic_kind(atomic_kind=atomic_kind,&
                            mass = mass ( ikind ) ) 
    END DO
 
    IF ( globenv%program_name == 'FIST' )  THEN
       CALL  init_atomic_kind_set(atomic_kind_set,globenv)
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                          epc_potential=epc_potential) 
          CALL set_potential(potential=epc_potential,&
                          qeff=charge(ikind))
       END DO
     END IF
  
     CALL convert_to_cp2k_units ( unit_type, part = particle_set )

!   *** Release work storage ***

     DEALLOCATE (atom_list,kind_label,kind_of,natom_of_kind,natom_of_kind_check,r,STAT=istat)
     DEALLOCATE (charge,mass,STAT=istat)
     IF (istat /= 0) THEN
       CALL stop_memory('read_atoms_input',&
       "atom_list,charge,kind_label,kind_of,mass,natom_of_kind_check,natom_of_kind,r")
     END IF
  END IF
END SUBROUTINE read_atoms_input

!******************************************************************************
!!****** atoms_input/read_md_restart [1.0] *
!!
!!   NAME
!!     read_md_restart
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    CJM
!!
!!*** *************************************************************************
SUBROUTINE read_md_restart ( rtype, globenv, particle_set, cell, cell_ref,  &
                             nhc_part, nhc_baro, npt_info ) 
   IMPLICIT NONE
! Arguments
  CHARACTER ( len = * ), INTENT ( IN ) :: rtype
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE ( particle_type ), DIMENSION ( : ), POINTER    :: particle_set
  TYPE ( cell_type ), POINTER :: cell, cell_ref
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part, nhc_baro
  TYPE ( npt_info_type ), POINTER :: npt_info ( :, : )
! locals
  INTEGER :: ierror, ilen, iw,  ia, ie, i, j, size1, size2, iunit
  LOGICAL :: exists
  CHARACTER ( len = 20 ) :: filen

  iw = globenv % scr

  CALL xstring ( globenv % project_name, ia, ie )
  filen = globenv % project_name ( ia:ie ) // '.restart'
  ia = MIN ( LEN ( filen ), 20 )
  WRITE ( iw, '( A, T61, A )' ) ' READ PIMD| Read pos/vel/eta from file ', &
           ADJUSTR ( filen ( 1:ia ) )

  INQUIRE ( FILE = filen, EXIST = exists )
  IF ( exists ) THEN
     iunit = get_unit()
     OPEN ( iunit, file = filen )

     IF ( rtype == 'INIT' ) THEN
        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % r ( : )
           particle_set ( i ) % f = 0._dbl
        END DO
        READ ( iunit, * ) cell % hmat ( 1, 1:3 )
        READ ( iunit, * ) cell % hmat ( 2, 1:3 )
        READ ( iunit, * ) cell % hmat ( 3, 1:3 )

     ELSE IF ( rtype == 'POS' )THEN
        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % r ( : )
           particle_set ( i ) % f = 0._dbl
        END DO
        READ ( iunit, * ) cell % hmat( 1, 1:3 )
        READ ( iunit, * ) cell % hmat( 2, 1:3 )
        READ ( iunit, * ) cell % hmat( 3, 1:3 )

        READ ( iunit, * ) cell_ref % hmat( 1, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 2, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 3, 1:3 )
     ELSE IF ( rtype == 'POSVEL' ) THEN

        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % r ( : )
        END DO
        READ ( iunit, * ) cell % hmat( 1, 1:3 )
        READ ( iunit, * ) cell % hmat( 2, 1:3 )
        READ ( iunit, * ) cell % hmat( 3, 1:3 )

        READ ( iunit, * ) cell_ref % hmat( 1, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 2, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 3, 1:3 )

        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % v ( : )
        END DO

     ELSE IF ( rtype == 'ALL' ) THEN

        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % r ( : )
           particle_set ( i ) % f = 0._dbl
        END DO
        READ ( iunit, * ) cell % hmat( 1, 1:3 )
        READ ( iunit, * ) cell % hmat( 2, 1:3 )
        READ ( iunit, * ) cell % hmat( 3, 1:3 )

        READ ( iunit, * ) cell_ref % hmat( 1, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 2, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 3, 1:3 )

        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % v ( : )
        END DO

! particle thermostat positions and velocities
        READ ( iunit, * ) size1, size2
        IF ( size1 /= SIZE ( nhc_part % nvt, 1 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for nhc_part' )
        IF ( size2 /= SIZE ( nhc_part % nvt, 2 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for nhc_part' )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_part % nvt ( i, j ) % eta 
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_part % nvt ( i, j ) % v 
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_part % nvt ( i, j ) % mass
          ENDDO
        ENDDO
! barostat thermostat positions and velocities
        READ ( iunit, * ) size1, size2
        IF ( size1 /= SIZE ( nhc_baro % nvt, 1 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for nhc_baro' )
        IF ( size2 /= SIZE ( nhc_baro % nvt, 2 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for nhc_baro' )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_baro % nvt ( i, j ) % eta 
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_baro % nvt ( i, j ) % v
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_baro % nvt ( i, j ) % mass 
          ENDDO
        ENDDO
! barostat velocities
        READ ( iunit, * ) size1, size2
        IF ( size1 /= SIZE ( npt_info, 1 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for veps' )
        IF ( size2 /= SIZE ( npt_info, 2 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for veps' )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) npt_info ( i, j ) % v
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * )  npt_info ( i, j ) % mass
          ENDDO
        ENDDO
     ELSE
        CALL stop_program ( 'read_md_restart', 'this RESTART type not programmed' )
     END IF

     CLOSE ( iunit )

  ELSE

     CALL stop_program ( 'read_md_restart', 'No MD RESTART file found ' )

  END IF


!..write some information to output
 IF ( globenv%ionode .AND. globenv%print_level>0 ) THEN
    WRITE ( iw, '( A )' ) ' RESTART MD | Box parameters [Angstrom]'
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' RESTART MD| ', ( cell % hmat ( 1, i ), i = 1, 3 )
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' RESTART MD| ', ( cell % hmat ( 2, i ), i = 1, 3 )
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' RESTART MD| ', ( cell % hmat ( 3, i ), i = 1, 3 )
    WRITE ( iw, '( A, T71, I10 )' ) &
        ' RESTART MD| Number of atoms read ', SIZE ( particle_set )
    IF ( globenv % print_level > 4 ) THEN
         IF ( rtype == 'POS' ) THEN
           CALL print_c ( iw, particle_set )
         ELSE IF ( rtype == 'POSVEL' ) THEN
            CALL print_cv ( iw, particle_set )
         END IF
    END IF
    WRITE ( iw, '()' )
  END IF

END SUBROUTINE read_md_restart

!******************************************************************************
!!****** atoms_input/read_pimd_atoms_input [1.0] *
!!
!!   NAME
!!     read_pimd_atoms_input
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    CJM
!!
!!*** *************************************************************************
SUBROUTINE read_pimd_atoms_input ( rtype, globenv, particle_set, cell,  &
                                  cell_ref, nhc_part, nhc_baro, npt_info ) 
   IMPLICIT NONE

! Arguments
  CHARACTER ( len = * ), INTENT ( IN ) :: rtype
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE ( particle_type ), DIMENSION ( : ), POINTER    :: particle_set
  TYPE ( cell_type ), POINTER :: cell, cell_ref
  TYPE ( lnhc_parameters_type ), POINTER :: nhc_part, nhc_baro
  TYPE ( npt_info_type ), POINTER :: npt_info ( :, : )
! locals
  INTEGER :: ierror, ilen, iw,  ia, ie, i, j, size1, size2, iunit
  LOGICAL :: exists
  CHARACTER ( len = 20 ) :: filen

  iw = globenv % scr

  CALL xstring ( globenv % project_name, ia, ie )
  filen = globenv % project_name ( ia:ie ) // '.dat'
  ia = MIN ( LEN ( filen ), 20 )
  WRITE ( iw, '( A, T61, A )' ) ' READ_PIMD| Read pos/vel/eta from file ', &
           ADJUSTR ( filen ( 1:ia ) )

  INQUIRE ( FILE = filen, EXIST = exists )
  IF ( exists ) THEN
     iunit = get_unit()
     OPEN ( iunit, file = filen )

     IF ( rtype == 'INIT' ) THEN
        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % r ( : )
           particle_set ( i ) % f = 0._dbl
        END DO
        READ ( iunit, * ) cell % hmat( 1, 1:3 )
        READ ( iunit, * ) cell % hmat( 2, 1:3 )
        READ ( iunit, * ) cell % hmat( 3, 1:3 )

     ELSE IF ( rtype == 'POS' )THEN
        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % r ( : )
           particle_set ( i ) % f = 0._dbl
        END DO
        READ ( iunit, * ) cell % hmat( 1, 1:3 )
        READ ( iunit, * ) cell % hmat( 2, 1:3 )
        READ ( iunit, * ) cell % hmat( 3, 1:3 )

        READ ( iunit, * ) cell_ref % hmat( 1, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 2, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 3, 1:3 )
     ELSE IF ( rtype == 'POSVEL' ) THEN

        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % r ( : )
        END DO
        READ ( iunit, * ) cell % hmat( 1, 1:3 )
        READ ( iunit, * ) cell % hmat( 2, 1:3 )
        READ ( iunit, * ) cell % hmat( 3, 1:3 )

        READ ( iunit, * ) cell_ref % hmat( 1, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 2, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 3, 1:3 )

        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % v ( : )
        END DO

     ELSE IF ( rtype == 'ALL' ) THEN

        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % r ( : )
           particle_set ( i ) % f = 0._dbl
        END DO
        READ ( iunit, * ) cell % hmat( 1, 1:3 )
        READ ( iunit, * ) cell % hmat( 2, 1:3 )
        READ ( iunit, * ) cell % hmat( 3, 1:3 )

        READ ( iunit, * ) cell_ref % hmat( 1, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 2, 1:3 )
        READ ( iunit, * ) cell_ref % hmat( 3, 1:3 )

        DO i = 1, SIZE ( particle_set )
           READ ( iunit, * ) particle_set ( i ) % v ( : )
        END DO

! particle thermostat positions and velocities
        READ ( iunit, * ) size1, size2
        IF ( size1 /= SIZE ( nhc_part % nvt, 1 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for nhc_part' )
        IF ( size2 /= SIZE ( nhc_part % nvt, 2 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for nhc_part' )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_part % nvt ( i, j ) % eta 
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_part % nvt ( i, j ) % v 
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_part % nvt ( i, j ) % mass 
          ENDDO
        ENDDO
! barostat thermostat positions and velocities
        READ ( iunit, * ) size1, size2
        IF ( size1 /= SIZE ( nhc_baro % nvt, 1 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for nhc_baro' )
        IF ( size2 /= SIZE ( nhc_baro % nvt, 2 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for nhc_baro' )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_baro % nvt ( i, j ) % eta 
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_baro % nvt ( i, j ) % v 
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) nhc_baro % nvt ( i, j ) % mass 
          ENDDO
        ENDDO
! barostat velocities
        READ ( iunit, * ) size1, size2
        IF ( size1 /= SIZE ( npt_info, 1 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for veps' )
        IF ( size2 /= SIZE ( npt_info, 2 ) ) &
        CALL stop_program ( 'read_md_restart', 'memory mismatch for veps' )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) npt_info ( i, j ) % v
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * )  npt_info ( i, j ) % mass 
          ENDDO
        ENDDO
     ELSE
        CALL stop_program ( 'read_md_restart', 'this RESTART type not programmed' )
     END IF

     CLOSE ( iunit )

  ELSE

     CALL stop_program ( 'read_md_restart', 'No MD RESTART file found ' )

  END IF

!..write some information to output
 IF ( globenv%ionode .AND. globenv%print_level>4 ) THEN
    WRITE ( iw, '( A )' ) ' READ PIMD| Box parameters [Angstrom]'
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' READ PIMD| ', ( cell % hmat ( 1, i ), i = 1, 3 )
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' RESTART MD| ', ( cell % hmat ( 2, i ), i = 1, 3 )
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' READ PIMD| ', ( cell % hmat ( 3, i ), i = 1, 3 )
    WRITE ( iw, '( A, T71, I10 )' ) &
        ' READ PIMD| Number of atoms read ', SIZE ( particle_set )
    IF ( globenv % print_level > 4 ) THEN
         IF ( rtype == 'POS' ) THEN
           CALL print_c ( iw, particle_set )
         ELSE IF ( rtype == 'POSVEL' ) THEN
            CALL print_cv ( iw, particle_set )
         END IF
    END IF
    WRITE ( iw, '()' )
  END IF

END SUBROUTINE read_pimd_atoms_input

!******************************************************************************

SUBROUTINE print_c ( iw, particle_set )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_type ), INTENT ( IN ) :: particle_set ( : )
  INTEGER, INTENT ( IN ) :: iw

! Locals
  INTEGER :: i, n

!------------------------------------------------------------------------------

  n = SIZE ( particle_set )
  WRITE ( iw, '( A )' ) ' RESTART MD | Atom coordinates [Angstrom]'
  DO i = 1, n
     WRITE ( iw, '( A, T26, I10, 3F15.5 )' ) ' RESTART MD| ', i, particle_set ( i ) % r ( : )
  END DO

END SUBROUTINE print_c

!******************************************************************************

SUBROUTINE print_cv ( iw, particle_set )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_type ), INTENT ( IN ) :: particle_set ( : )
  INTEGER, INTENT ( IN ) :: iw

! Locals
  INTEGER :: i, n

!------------------------------------------------------------------------------

  n = SIZE ( particle_set )
  WRITE ( iw, '( A )' ) ' RESTART MD | Atom coordinates [Angstrom]'
  DO i = 1, n
     WRITE ( iw, '( A, T8, I8, 3F10.4, 5X, 3F10.4 )' ) &
          ' RESTART MD| ', i, particle_set ( i ) % r ( : ), &
                        particle_set ( i ) % v ( : )
  END DO

END SUBROUTINE print_cv

!******************************************************************************

END MODULE atoms_input
