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
!!
!!   SOURCE
!******************************************************************************

MODULE atoms_input

  USE atomic_kind_types,        ONLY : allocate_atomic_kind_set,&
                                       atomic_kind_type,&
                                       get_atomic_kind,&
                                       init_atomic_kind_set,&
                                       set_atomic_kind
  USE external_potential_types, ONLY : epc_potential_type,&
                                       set_potential
  USE global_types,             ONLY : global_environment_type
  USE input_types,              ONLY : setup_parameters_type
  USE kinds,                    ONLY : dbl
  USE memory_utilities,         ONLY : reallocate
  USE message_passing,          ONLY : mp_bcast
  USE molecule_types,           ONLY : molecule_type
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
  USE physcon,                   ONLY : massunit
  USE string_utilities,         ONLY : xstring,&
                                       uppercase, &
                                       integer_to_string, &
                                       string_to_integer, &
                                       str_search
  USE termination,              ONLY : stop_memory, &
                                       stop_program
  USE util,                     ONLY : get_unit

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_coord_vel, system_type

  TYPE system_type
     REAL ( dbl ) :: box ( 3, 3 ), box_ref ( 3, 3 )
     CHARACTER ( LEN = 3 ) :: ptype
     CHARACTER ( LEN = 20 ) :: rtype
     INTEGER :: n, size1_nvt_part, size2_nvt_part
     INTEGER :: size1_nvt_baro, size2_nvt_baro
     INTEGER :: size1_npt, size2_npt
     REAL ( dbl ), POINTER :: c ( :, : )
     REAL ( dbl ), POINTER :: v ( :, : )
     REAL ( dbl ), POINTER :: veps ( :, : )
     REAL ( dbl ), POINTER :: mass_eps ( :, : )
     REAL ( dbl ), POINTER :: eta_part ( :, : )
     REAL ( dbl ), POINTER :: veta_part ( :, : )
     REAL ( dbl ), POINTER :: mass_eta_part ( :, : )
     REAL ( dbl ), POINTER :: eta_baro ( :, : )
     REAL ( dbl ), POINTER :: veta_baro ( :, : )
     REAL ( dbl ), POINTER :: mass_eta_baro ( :, : )
  END TYPE system_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** atoms_input/set_particle_properties [1.0] *
!!
!!   NAME
!!     set_particle_properties
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   CJM SEPT-26-2002 
!!
!!*** *************************************************************************

!!SUBROUTINE set_particle_properties ( setup, mol_setup, pstat )
!!  IMPLICIT NONE 
!!  TYPE ( setup_parameters_type ), INTENT ( IN )  :: setup
!!  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( IN )  :: mol_setup
!!  TYPE ( particle_prop_type ), DIMENSION ( : ), POINTER :: pstat 

!locals
!!  INTEGER :: nmol_type, natom_type, ios, i, j, n  

!!  nmol_type = size(mol_setup)
!!  natom_type = setup % natom_type

!..define pstat to hold all information (static) on different particle types
!!  ALLOCATE ( pstat ( natom_type ), STAT = ios )
!!   IF ( ios /= 0 ) CALL stop_memory ( 'atoms_input', 'pstat', natom_type )

!..this defines the order of the atoms!
!!  DO i = 1, nmol_type
!!     DO j = 1, mol_setup(i) % molpar % natom
!!        n = str_search ( setup % atom_names,natom_type, &
!!                         mol_setup(i) % molpar % aname(j) )
!!        pstat(n) % charge = mol_setup(i) % molpar % acharge(j)
!!        pstat(n) % mass = mol_setup(i) % molpar % aweight(j)
!!        pstat(n) % ptype = n
!!     END DO
!!  END DO

!!END SUBROUTINE set_particle_properties 

!******************************************************************************
!!****** atoms_input/read_coord_vel [1.0] *
!!
!!   NAME
!!     read_coord_vel
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!   MODIFICATION HISTORY
!!   cjm, Feb-20-2001 : broadcasts extended variables in
!!   system_type
!!
!!  OPTIONS
!!    either read the section &atoms in the input file or
!!    read the coordinates from the file project_name.dat
!!
!!  INPUTS
!!I---------------------------------------------------------------------------I
!!I SECTION: &atoms ... &end                                                  I
!!I     file:   file_name                                                     I
!!I     cell:   b11 b12 b13 &                                                 I
!!I             b21 b22 b23 &                                                 I
!!I             b31 b32 b33                                                   I
!!I---------------------------------------------------------------------------I
!!*** *************************************************************************

SUBROUTINE read_coord_vel (particle_set, atomic_kind_set, atype, filen, &
                           mol_setup, atom_names, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( system_type ), INTENT ( INOUT ) :: atype
  CHARACTER ( LEN = * ), INTENT ( INOUT ) :: filen
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( IN ) :: mol_setup
  CHARACTER ( len = * ), DIMENSION ( : ), INTENT ( IN ) :: atom_names


! Locals
  INTEGER :: ierror, ilen, iw, source, allgrp, ia, ie, i
  INTEGER :: size1, size2
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

!   ---------------------------------------------------------------------------

  nkind = 0
  natom = 0
  iat = 0
  COORDS = .FALSE.
  INIT = .TRUE.

  iw = globenv % scr

!..defaults
  NULLIFY ( atype % c )
  NULLIFY ( atype % v )
  NULLIFY ( atype % veps )
  NULLIFY ( atype % mass_eps )
  NULLIFY ( atype % eta_part )
  NULLIFY ( atype % veta_part )
  NULLIFY ( atype % mass_eta_part )
  NULLIFY ( atype % eta_baro )
  NULLIFY ( atype % veta_baro )
  NULLIFY ( atype % mass_eta_baro )
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

     IF( globenv % ionode) THEN
       WRITE ( iw, '( A )' ) ' ATOM| No input section &COORD found '
       IF ( filen == ' ' ) THEN
          CALL xstring ( globenv % project_name, ia, ie )
          filen = globenv % project_name ( ia:ie ) // '.dat'
       END IF
       ia = MIN ( LEN ( filen ), 20 )
       WRITE ( iw, '( A, T61, A )' ) ' ATOM| Try to read default file ', &
            ADJUSTR ( filen ( 1:ia ) )
       CALL read_file ( filen, atype )
     END IF
!..broadcast the input data to all nodes
#if defined(__parallel)
     source = globenv % source
     allgrp = globenv % group
     CALL mp_bcast ( atype % box, source, allgrp )
     CALL mp_bcast ( atype % box_ref, source, allgrp )
     CALL mp_bcast ( atype % ptype, source, allgrp )
     CALL mp_bcast ( atype % rtype, source, allgrp )
     CALL mp_bcast ( atype % n, source, allgrp )
     IF ( .NOT. ASSOCIATED ( atype % c ) ) THEN
        ALLOCATE ( atype % c(1:3,1:atype % n ), STAT = istat )
        IF ( istat /= 0 ) CALL stop_memory ( 'ATOM', 'atype%c', 3 * atype % n )
     END IF
     CALL mp_bcast ( atype % c, source, allgrp )
     IF ( atype % rtype == 'POSVEL' ) THEN
        IF ( .NOT. ASSOCIATED ( atype % v ) ) THEN
           ALLOCATE (atype % v(1:3,1:atype % n),STAT=istat)
           IF ( istat /= 0 ) CALL stop_memory ( 'ATOM', 'atype%v', 3 * atype % n )
        END IF
        CALL mp_bcast ( atype % v, source, allgrp )
     END IF
     IF ( atype % rtype == 'ALL' ) THEN
        IF ( .NOT. globenv % ionode ) THEN
           ALLOCATE (atype % v(1:3,1:atype % n),STAT=istat)
           IF ( istat /= 0 ) CALL stop_memory ( 'ATOM', 'atype%v', 3 * atype % n )
           size1 = atype % size1_nvt_part
           size2 = atype % size2_nvt_part
           ALLOCATE (atype % veta_part (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%veta_part', size1 * size2 )
           ALLOCATE (atype % eta_part (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%eta_part', size1 * size2 )
           ALLOCATE (atype % mass_eta_part (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%mass_eta_part', size1 * size2 )
           ALLOCATE (atype % veta_baro (size1, size2 ),STAT=istat)
           size1 = atype % size1_nvt_baro
           size2 = atype % size2_nvt_baro
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%veta_baro', size1 * size2 )
           ALLOCATE (atype % eta_baro (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%eta_baro', size1 * size2 )
           ALLOCATE (atype % mass_eta_baro (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%mass_eta_baro', size1 * size2 )
           size1 = atype % size1_npt
           size2 = atype % size2_npt
           ALLOCATE (atype % veps (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%veps', size1 * size2 )
           ALLOCATE (atype % mass_eps (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%mass_eps', size1 * size2 )
        END IF
        CALL mp_bcast ( atype % v, source, allgrp )
        CALL mp_bcast ( atype % veta_part, source, allgrp )
        CALL mp_bcast ( atype % eta_part, source, allgrp )
        CALL mp_bcast ( atype % veta_baro, source, allgrp )
        CALL mp_bcast ( atype % eta_baro, source, allgrp )
        CALL mp_bcast ( atype % veps, source, allgrp )
        CALL mp_bcast ( atype % mass_eps, source, allgrp )
        CALL mp_bcast ( atype % mass_eta_part, source, allgrp )
        CALL mp_bcast ( atype % mass_eta_baro, source, allgrp )
     END IF
#endif
! allocate particle set

     CALL allocate_particle_set(particle_set, atype % n)
! fill position and velocities
     CALL put_posvel_in_particle_set(particle_set, atype )
     nkind = SIZE(atom_names)
     natom = atype % n
     CALL reallocate(atom_list,1,nkind,1,natom)
     CALL reallocate(kind_label,1,max_name_length,1,nkind)
     CALL reallocate(kind_of,1,natom)
     CALL reallocate(natom_of_kind,1,nkind)
     CALL reallocate(charge,1,nkind)
     CALL reallocate(mass,1,nkind)

!   *** Allocate the data structure for the atomic kind information ***

    CALL allocate_atomic_kind_set(atomic_kind_set,nkind)

!   *** setup atomic_kind_based atom lists
 
     nmol_type = SIZE(mol_setup)
     natom_types = SIZE(atom_names)
     DO i = 1, nmol_type
       DO k = 1, mol_setup(i) % num_mol
         DO j = 1, mol_setup(i) % molpar % natom
           iat = iat + 1
           n = str_search ( atom_names, natom_types, &
                mol_setup ( i ) % molpar % aname ( j ) )
           charge(n)= mol_setup(i) % molpar % acharge(j)
!dbg gt this should be done using general unit conversion routine
! au:
!            mass(n)= massunit*mol_setup(i) % molpar % aweight(j)
! kelvin:
            mass(n)= mol_setup(i) % molpar % aweight(j)
!dbg
           natom_of_kind(n) = natom_of_kind(n) + 1
           atom_list(n,natom_of_kind(n)) = iat
           kind_of(iat) = n
        END DO
      END DO
    END DO

!   *** Initialize the data structure for atomic kind information ***

     DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       string = atom_names (ikind) 
       CALL set_atomic_kind(atomic_kind=atomic_kind,&
                            kind_number=ikind,&
                            mass=mass(ikind),&
                            name=string,&
                            atom_list=atom_list(ikind,1:natom_of_kind(ikind)))
     END DO

     IF (globenv%program_name == "FIST") THEN
       CALL  init_atomic_kind_set(atomic_kind_set,globenv)
       DO ikind=1,nkind
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           epc_potential=epc_potential) 
         CALL set_potential(potential=epc_potential,&
                           qeff=charge(ikind))
       END DO
     END IF

      DO iatom=1,natom
        ikind = kind_of(iatom)
        particle_set(iatom) % abase = 0
        particle_set(iatom) % kind = ikind
        particle_set(iatom) % grid_id = 0
        particle_set(iatom)%atomic_kind => atomic_kind_set(ikind)
        particle_set(iatom)%iatom = iatom
        particle_set(iatom)%first_cgf = 0
        particle_set(iatom)%first_sgf = 0
        NULLIFY ( particle_set ( iatom ) % coef_list )
     END DO

   ELSE

!   ******************************************************************
!   **** Reads cordinates from the Input file in QS-friendly style ***  
!   ******************************************************************

     COORDS = .TRUE.
     CALL read_line
     DO WHILE (test_next()/='X')
       ilen = 8
       IF (test_next()=="C") THEN
          CALL get_next(current_label,ilen)
       ELSE IF (test_next()=="N") THEN
          CALL get_next (z)
          current_label = ptable(z)%symbol
       END IF

       CALL uppercase(current_label)
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
       CALL uppercase(label)

       new_kind = .TRUE.

       DO ikind=1,nkind
         CALL integer_to_string(kind_label(:,ikind),string)
         CALL uppercase(string)
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

   END IF
   CALL parser_end
!..end of parsing the input section

  IF (COORDS) THEN
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
       particle_set(iatom) % kind = ikind
       particle_set(iatom) % grid_id = 0
       particle_set(iatom)%atomic_kind => atomic_kind_set(ikind)
       particle_set(iatom)%iatom = iatom
       particle_set(iatom)%first_cgf = 0
       particle_set(iatom)%first_sgf = 0
       particle_set(iatom)%f(:) = 0.0_dbl
       particle_set(iatom)%r(:) = r(:,iatom)
       particle_set(iatom)%v(:) = 0.0_dbl
       NULLIFY ( particle_set ( iatom ) % coef_list )
    END DO

!   *** check if the lables in &COORD match with the &MOLECULE setup
!   *** and store the charge and mass read in the &MOLECULE section

     ALLOCATE (natom_of_kind_check(nkind), STAT=istat) 
     nmol_type = SIZE (mol_setup) 
     natom_types = SIZE(atom_names)
     DO i = 1, nmol_type
       DO k = 1, mol_setup(i) % num_mol
         DO j = 1, mol_setup(i) % molpar % natom
           iat = iat + 1
           n = str_search ( atom_names, natom_types, &
                mol_setup ( i ) % molpar % aname ( j ) )
           charge(n)= mol_setup(i) % molpar % acharge(j)
!dbg this should be done using general unit conversion routine
!AU
!            mass(n)= massunit* mol_setup(i) % molpar % aweight(j)
! kelvin
             mass(n)= mol_setup(i) % molpar % aweight(j)
! gt
           ikind = particle_set ( iat) % kind        
           IF (ikind/=n) CALL stop_program ( 'read_coord_vel', &
                'kind label ordering in &COORD and &MOL dont match ')
           CALL integer_to_string(kind_label(:,ikind),string)
           IF (string/=mol_setup(i)%molpar%aname(j)) CALL stop_program ( 'read_coord_vel', &
                'inconsistent kind labels in &COORD and &MOL- check input')
           natom_of_kind_check(ikind) = natom_of_kind_check(ikind) + 1
           iatom =  atom_list(ikind,natom_of_kind_check(ikind)) 
           IF (iatom/=iat) CALL stop_program ( 'read_coord_vel', &
                'atom lists in &COORD and &MOL dont match ')
        END DO
      END DO
    END DO
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL set_atomic_kind(atomic_kind=atomic_kind,&
                            mass=mass(ikind)) 
    END DO
    IF (globenv%program_name == "FIST") THEN
       CALL  init_atomic_kind_set(atomic_kind_set,globenv)
       DO ikind=1,nkind
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           epc_potential=epc_potential) 
         CALL set_potential(potential=epc_potential,&
                           qeff=charge(ikind))
       END DO
    END IF

!   *** Release work storage ***

    DEALLOCATE (atom_list,kind_label,kind_of,natom_of_kind,natom_of_kind_check,r,STAT=istat)
    DEALLOCATE (charge,mass,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory('ATOM',&
                 "atom_list,charge,kind_label,kind_of,mass,natom_of_kind_check,natom_of_kind,r")
    END IF

 END IF

! now that particle set has been initialized using the input cooridnates,
! check if we have to restart from the position/velocities in the .dat file 
 
 IF (atype%rtype /= 'INIT') THEN 
   INIT = .FALSE.
 END IF

  IF ((COORDS).AND.(.NOT.(INIT))) THEN
   CALL xstring ( globenv % project_name, ia, ie )
   filen = globenv % project_name ( ia:ie ) // '.dat'
   ia = MIN ( LEN ( filen ), 20 )
   WRITE ( iw, '( A, T61, A )' ) ' ATOM| Read pos/vel/eta from file ', &
            ADJUSTR ( filen ( 1:ia ) )
   CALL read_file ( filen, atype )

!..broadcast the input data to all nodes
#if defined(__parallel)
     source = globenv % source
     allgrp = globenv % group
     CALL mp_bcast ( atype % box, source, allgrp )
     CALL mp_bcast ( atype % box_ref, source, allgrp )
     CALL mp_bcast ( atype % ptype, source, allgrp )
     CALL mp_bcast ( atype % rtype, source, allgrp )
     CALL mp_bcast ( atype % n, source, allgrp )
     IF ( .NOT. ASSOCIATED ( atype % c ) ) THEN
        ALLOCATE ( atype % c(1:3,1:atype % n ), STAT = istat )
        IF ( istat /= 0 ) CALL stop_memory ( 'ATOM', 'atype%c', 3 * atype % n )
     END IF
     CALL mp_bcast ( atype % c, source, allgrp )
     IF ( atype % rtype == 'POSVEL' ) THEN
        IF ( .NOT. ASSOCIATED ( atype % v ) ) THEN
           ALLOCATE (atype % v(1:3,1:atype % n),STAT=istat)
           IF ( istat /= 0 ) CALL stop_memory ( 'ATOM', 'atype%v', 3 * atype % n )
        END IF
        CALL mp_bcast ( atype % v, source, allgrp )
     END IF
     IF ( atype % rtype == 'ALL' ) THEN
        IF ( .NOT. globenv % ionode ) THEN
           ALLOCATE (atype % v(1:3,1:atype % n),STAT=istat)
           IF ( istat /= 0 ) CALL stop_memory ( 'ATOM', 'atype%v', 3 * atype % n )
           size1 = atype % size1_nvt_part
           size2 = atype % size2_nvt_part
           ALLOCATE (atype % veta_part (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%veta_part', size1 * size2 )
           ALLOCATE (atype % eta_part (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%eta_part', size1 * size2 )
           ALLOCATE (atype % mass_eta_part (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%mass_eta_part', size1 * size2 )
           ALLOCATE (atype % veta_baro (size1, size2 ),STAT=istat)
           size1 = atype % size1_nvt_baro
           size2 = atype % size2_nvt_baro
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%veta_baro', size1 * size2 )
           ALLOCATE (atype % eta_baro (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%eta_baro', size1 * size2 )
           ALLOCATE (atype % mass_eta_baro (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%mass_eta_baro', size1 * size2 )
           size1 = atype % size1_npt
           size2 = atype % size2_npt
           ALLOCATE (atype % veps (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%veps', size1 * size2 )
           ALLOCATE (atype % mass_eps (size1, size2 ),STAT=istat)
           IF ( istat /= 0 )  &
            CALL stop_memory ( 'ATOM', 'atype%mass_eps', size1 * size2 )
        END IF
        CALL mp_bcast ( atype % v, source, allgrp )
        CALL mp_bcast ( atype % veta_part, source, allgrp )
        CALL mp_bcast ( atype % eta_part, source, allgrp )
        CALL mp_bcast ( atype % veta_baro, source, allgrp )
        CALL mp_bcast ( atype % eta_baro, source, allgrp )
        CALL mp_bcast ( atype % veps, source, allgrp )
        CALL mp_bcast ( atype % mass_eps, source, allgrp )
        CALL mp_bcast ( atype % mass_eta_part, source, allgrp )
        CALL mp_bcast ( atype % mass_eta_baro, source, allgrp )
     END IF
#endif
! Check if the informations in the .dat files are consistent
! with the ones in the &COORD section
   IF (atype % n /= natom)  CALL stop_program ( 'read_coord_vel', &
                ' inconsistent number of atoms specified ')
   
! if yes, fill the position and/or velocities in particle_set
   CALL put_posvel_in_particle_set(particle_set, atype )

 END IF

!..write some information to output
 IF (globenv%ionode .and. globenv%print_level>0 .and.(.not.INIT)) THEN
    WRITE ( iw, '( A )' ) ' ATOM| Box parameters [Angstrom]'
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' ATOM| ', ( atype % box ( 1, i ), i = 1, 3 )
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' ATOM| ', ( atype % box ( 2, i ), i = 1, 3 )
    WRITE ( iw, '( A, T36, 3F15.5 )' ) &
         ' ATOM| ', ( atype % box ( 3, i ), i = 1, 3 )
    WRITE ( iw, '( A, T71, I10 )' ) &
        ' ATOM| Number of atoms read ', atype % n
    IF ( globenv % print_level > 4 ) THEN
         IF ( atype % rtype == 'POS' ) THEN
           CALL print_c ( iw, atype % c )
         ELSE IF ( atype % rtype == 'POSVEL' ) THEN
            CALL print_cv ( iw, atype % c, atype % v )
         END IF
    END IF
    WRITE ( iw, '()' )
  END IF

END SUBROUTINE read_coord_vel

!******************************************************************************

SUBROUTINE read_file ( filen, atype )

  IMPLICIT NONE

! Arguments
  CHARACTER ( LEN = * ) :: filen
  TYPE ( system_type ), INTENT ( INOUT ) :: atype

! Locals
  INTEGER :: iunit, i, j, ios, size1, size2
  LOGICAL :: exists

!------------------------------------------------------------------------------

  INQUIRE ( FILE = filen, EXIST = exists )
  IF ( exists ) THEN
     iunit = get_unit()
     OPEN ( iunit, file = filen )
     READ ( iunit, * ) atype % n

     IF (.NOT. ASSOCIATED (atype % c) ) THEN
        ALLOCATE ( atype % c ( 1:3, 1:atype % n ), STAT = ios )
        IF ( ios /= 0 ) CALL stop_memory ( 'ATOM', 'atype%c', 3 * atype % n )
        NULLIFY ( atype % v )
     ELSE IF (atype % n /= SIZE( atype % c, 2 ) ) THEN
        CALL stop_program ( 'atoms_input','Inconsistent number of atoms' )
     END IF

     IF ( atype % rtype == 'INIT' ) THEN
        DO i = 1, atype % n
           READ ( iunit, * ) atype % c ( 1:3, i )
        END DO
        READ ( iunit, * ) atype % box ( 1, 1:3 )
        READ ( iunit, * ) atype % box ( 2, 1:3 )
        READ ( iunit, * ) atype % box ( 3, 1:3 )

     ELSE IF ( atype % rtype == 'POS' )THEN
        DO i = 1, atype % n
           READ ( iunit, * ) atype % c ( 1:3, i )
        END DO
        READ ( iunit, * ) atype % box ( 1, 1:3 )
        READ ( iunit, * ) atype % box ( 2, 1:3 )
        READ ( iunit, * ) atype % box ( 3, 1:3 )

        READ ( iunit, * ) atype % box_ref ( 1, 1:3 )
        READ ( iunit, * ) atype % box_ref ( 2, 1:3 )
        READ ( iunit, * ) atype % box_ref ( 3, 1:3 )
     ELSE IF ( atype % rtype == 'POSVEL' ) THEN
        IF (.NOT. ASSOCIATED (atype % v) ) THEN
           ALLOCATE ( atype % v ( 1:3, 1:atype % n ), STAT = ios )
           IF ( ios /= 0 ) &
                CALL stop_memory ( 'ATOM', 'atype%v', 3 * atype % n )
        END IF

        DO i = 1, atype % n
           READ ( iunit, * ) atype % c ( 1:3, i )
        END DO
        READ ( iunit, * ) atype % box ( 1, 1:3 )
        READ ( iunit, * ) atype % box ( 2, 1:3 )
        READ ( iunit, * ) atype % box ( 3, 1:3 )

        READ ( iunit, * ) atype % box_ref ( 1, 1:3 )
        READ ( iunit, * ) atype % box_ref ( 2, 1:3 )
        READ ( iunit, * ) atype % box_ref ( 3, 1:3 )
        DO i = 1, atype % n
           READ ( iunit, * ) atype % v ( 1:3, i )
        END DO
     ELSE IF ( atype % rtype == 'ALL' ) THEN
        IF (.NOT. ASSOCIATED (atype % v) ) THEN
           ALLOCATE ( atype % v ( 1:3, 1:atype % n ), STAT = ios )
           IF ( ios /= 0 ) &
                CALL stop_memory ( 'ATOM', 'atype%v', 3 * atype % n )
        END IF

        DO i = 1, atype % n
           READ ( iunit, * ) atype % c ( 1:3, i )
        END DO
        READ ( iunit, * ) atype % box ( 1, 1:3 )
        READ ( iunit, * ) atype % box ( 2, 1:3 )
        READ ( iunit, * ) atype % box ( 3, 1:3 )

        READ ( iunit, * ) atype % box_ref ( 1, 1:3 )
        READ ( iunit, * ) atype % box_ref ( 2, 1:3 )
        READ ( iunit, * ) atype % box_ref ( 3, 1:3 )
        DO i = 1, atype % n
           READ ( iunit, * ) atype % v ( 1:3, i )
        END DO
! particle thermostat positions and velocities
        READ ( iunit, * ) atype % size1_nvt_part, atype % size2_nvt_part
        size1 = atype % size1_nvt_part
        size2 = atype % size2_nvt_part
        ALLOCATE ( atype % veta_part ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%veta_part', size1 * size2 )
        ALLOCATE ( atype % eta_part ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%eta_part', size1 * size2 )
        ALLOCATE ( atype % mass_eta_part ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%mass_eta_part', size1 * size2 )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % eta_part ( i, j )
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % veta_part ( i, j )
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % mass_eta_part ( i, j )
          ENDDO
        ENDDO
! barostat thermostat positions and velocities
        READ ( iunit, * ) atype % size1_nvt_baro, atype % size2_nvt_baro
        size1 = atype % size1_nvt_baro
        size2 = atype % size2_nvt_baro
        ALLOCATE ( atype % veta_baro ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%veta_baro', size1 * size2 )
        ALLOCATE ( atype % eta_baro ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%eta_baro', size1 * size2 )
        ALLOCATE ( atype % mass_eta_baro ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%mass_eta_baro', size1 * size2 )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % eta_baro ( i, j )
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % veta_baro ( i, j )
          ENDDO
        ENDDO
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % mass_eta_baro ( i, j )
          ENDDO
        ENDDO
! barostat velocities
        READ ( iunit, * ) atype % size1_npt, atype % size2_npt
        size1 = atype % size1_npt
        size2 = atype % size2_npt
        ALLOCATE ( atype % veps ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%veps', size1 * size2 )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % veps ( i, j )
          ENDDO
        ENDDO
        ALLOCATE ( atype % mass_eps ( size1, size2 ), STAT = ios )
        IF ( ios /= 0 ) &
             CALL stop_memory ( 'ATOM', 'atype%mass_eps', size1 * size2 )
        DO i = 1, size1
          DO j = 1, size2
            READ ( iunit, * ) atype % mass_eps ( i, j )
          ENDDO
        ENDDO
     ELSE
        CALL stop_program ( 'ATOM', 'this rtype not programmed' )
     END IF

     CLOSE ( iunit )

  ELSE

     CALL stop_program ( 'ATOM', 'No information on atoms found ' )

  END IF

END SUBROUTINE read_file

!******************************************************************************
SUBROUTINE put_posvel_in_particle_set ( part, atype )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
  TYPE ( system_type ), INTENT ( INOUT ) :: atype

!locals
  INTEGER :: natoms,iat

!------------------------------------------------------------------------------

  natoms = SIZE ( atype % c, 2 )
  DO iat = 1, natoms
     part(iat) % iatom = iat
     part(iat) % r(1) = atype % c(1,iat)
     part(iat) % r(2) = atype % c(2,iat)
     part(iat) % r(3) = atype % c(3,iat)
  END DO

! set initial velocities to 0
  part ( : ) % v ( 1 ) =  0._dbl
  part ( : ) % v ( 2 ) =  0._dbl
  part ( : ) % v ( 3 ) =  0._dbl

! set forces to 0
  part ( : ) % f ( 1 ) =  0._dbl
  part ( : ) % f ( 2 ) =  0._dbl
  part ( : ) % f ( 3 ) =  0._dbl

  IF ( ( atype%rtype=='POSVEL' ) .OR. ( atype%rtype=='ALL' ) ) THEN
     DO iat = 1, natoms
        part ( iat ) % v ( 1 ) = atype % v ( 1, iat )
        part ( iat ) % v ( 2 ) = atype % v ( 2, iat )
        part ( iat ) % v ( 3 ) = atype % v ( 3, iat )
     END DO
  END IF

 END SUBROUTINE put_posvel_in_particle_set
!******************************************************************************

SUBROUTINE print_c ( iw, c )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: iw
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: c

! Locals
  INTEGER :: i, n

!------------------------------------------------------------------------------

  n = SIZE ( c, 2 )
  WRITE ( iw, '( A )' ) ' ATOM| Atom coordinates [Angstrom]'
  DO i = 1, n
     WRITE ( iw, '( A, T26, I10, 3F15.5 )' ) ' ATOM| ', i, c ( :, i )
  END DO

END SUBROUTINE print_c

!******************************************************************************

SUBROUTINE print_cv ( iw, c, v )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: iw
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: c
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: v

! Locals
  INTEGER :: i, n

!------------------------------------------------------------------------------

  n = SIZE ( c, 2 )
  WRITE ( iw, '( A )' ) ' ATOM| Atom coordinates [Angstrom]'
  DO i = 1, n
     WRITE ( iw, '( A, T8, I8, 3F10.4, 5X, 3F10.4 )' ) &
          ' ATOM| ', i, c ( :, i ), v ( :, i )
  END DO

END SUBROUTINE print_cv

!******************************************************************************

END MODULE atoms_input
