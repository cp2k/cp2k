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
  USE global_types,             ONLY : global_environment_type
  USE input_types,              ONLY : setup_parameters_type
  USE kinds,                    ONLY : dbl,&
                                       default_string_length
  USE memory_utilities,         ONLY : reallocate
  USE molecule_types,           ONLY : molecule_type
  USE extended_system_types,   ONLY : lnhc_parameters_type, &
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
                                       str_search
  USE termination,              ONLY : stop_memory, &
                                       stop_program
  USE topology_types,           ONLY : topology_parameters_type
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

SUBROUTINE read_atoms_input ( topology, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( topology_parameters_type ) :: topology
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: ierror, ilen, iw
  INTEGER, PARAMETER :: nblock_atom = 1000

  CHARACTER(LEN=6)  :: label
  INTEGER                         :: iatom,istat,natom

  INTEGER :: atm_equal_mol

!   ---------------------------------------------------------------------------

  natom = 0

  iw = globenv % scr

  NULLIFY(topology%label_resname)
  NULLIFY(topology%label_atmname)
  NULLIFY(topology%r)
  NULLIFY(topology%atm_mass)
  NULLIFY(topology%atm_charge)
!..defaults

  CALL reallocate(topology%label_resname,1,nblock_atom)
  CALL reallocate(topology%label_atmname,1,nblock_atom)
  CALL reallocate(topology%r,1,3,1,nblock_atom)
  CALL reallocate(topology%atm_mass,1,nblock_atom)
  CALL reallocate(topology%atm_charge,1,nblock_atom)

!..parse the input section
  label = '&COORD'
  CALL parser_init ( globenv % input_file_name, globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
     CALL stop_program ( 'read_atoms_input','No &COORD section found')
  ELSE
!   ******************************************************************
!   **** Reads coordinates from &COORD section ***  
!   ******************************************************************
     CALL read_line
     DO WHILE (test_next()/='X')

       natom = natom + 1

       IF (natom > SIZE(topology%r)) THEN
        CALL reallocate(topology%label_resname,1,natom+nblock_atom)
        CALL reallocate(topology%label_atmname,1,natom+nblock_atom)
        CALL reallocate(topology%r,1,3,1,natom+nblock_atom)
        CALL reallocate(topology%atm_mass,1,natom+nblock_atom)
        CALL reallocate(topology%atm_charge,1,natom+nblock_atom)
       END IF

       ilen = default_string_length
       topology%label_resname(natom) = ""
       topology%label_atmname(natom) = ""

       ! For default, set atom name to molecule name
       CALL get_next(topology%label_atmname(natom),ilen)
       topology%label_resname(natom) = topology%label_atmname(natom)

!     *** Read the atomic coordinates ***
       CALL get_next(topology%r(1,natom))
       CALL get_next(topology%r(2,natom))
       CALL get_next(topology%r(3,natom))

       !Future user optional inputs
       IF (test_next()=='C') THEN
         CALL get_next(topology%label_resname(natom),ilen)
       END IF
       topology%atm_mass(natom) = 0.0_dbl 
       topology%atm_charge(natom) = 0.0_dbl 

     call read_line
     END DO

    CALL parser_end
  END IF
!..end of parsing the input section

  topology%natoms = natom
  CALL reallocate(topology%label_resname,1,natom)
  CALL reallocate(topology%label_atmname,1,natom)
  CALL reallocate(topology%r,1,3,1,natom)
  CALL reallocate(topology%atm_mass,1,natom)
  CALL reallocate(topology%atm_charge,1,natom)

  DO iatom=1,natom
    CALL convert_to_cp2k_units(topology%cell%unit_of_length_name,&
                               length=topology%r(1,iatom))
    CALL convert_to_cp2k_units(topology%cell%unit_of_length_name,&
                               length=topology%r(2,iatom))
    CALL convert_to_cp2k_units(topology%cell%unit_of_length_name,&
                               length=topology%r(3,iatom))
  END DO

  atm_equal_mol = 0
  DO iatom=1,natom
    IF(topology%label_resname(iatom)==topology%label_atmname(iatom))THEN
      atm_equal_mol = atm_equal_mol + 1
      topology%bondparm_factor=0.0_dbl
    END IF
  END DO

  IF ((atm_equal_mol /= 0) .AND. (atm_equal_mol /= natom)) THEN
    CALL stop_program ('read_atoms_input','check coord optional molname')
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
