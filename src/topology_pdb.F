!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/topology_pdb [1.0] *
!!
!!   NAME
!!     topology_pdb
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
MODULE topology_pdb

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
  USE physcon,                  ONLY : massunit
  USE simulation_cell,          ONLY : cell_type
  USE string_utilities,         ONLY : xstring,&
                                       integer_to_string, &
                                       string_to_integer, &
                                       str_search, &
                                       uppercase
  USE termination,              ONLY : stop_memory, &
                                       stop_program
  USE topology_util,            ONLY : azero
  USE util,                     ONLY : get_unit

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_coordinate_pdb

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** topology_pdb/read_coordinate_pdb [1.0] *
!!
!!   NAME
!!     read_coordinate_pdb
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!  OPTIONS
!!
!!  INPUTS
!!*** *************************************************************************

SUBROUTINE read_coordinate_pdb ( particle_set, atomic_kind_set, &
                              mol_setup, atom_names, setup, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( IN ) :: mol_setup
  CHARACTER ( len = * ), DIMENSION ( : ), INTENT ( IN ) :: atom_names
  TYPE ( setup_parameters_type ) :: setup



! Locals
  INTEGER :: ierror, ilen, iw,  ia, ie, i
  INTEGER, PARAMETER :: max_name_length = 10
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

  CHARACTER (LEN=60), POINTER :: label1(:),label2(:)
  INTEGER, POINTER :: num1(:),num2(:)
  REAL (dbl), POINTER :: pos1(:),pos2(:),pos3(:)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  iw = globenv % scr

  nkind = 0
  natom = 0
  iat = 0
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Change to REALLOCATE later
  ALLOCATE(num1(100000))
  ALLOCATE(num2(100000))
  ALLOCATE(label1(100000))
  ALLOCATE(label2(100000))
  ALLOCATE(pos1(100000))
  ALLOCATE(pos2(100000))
  ALLOCATE(pos3(100000))
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Parse the file right now and put all info into temporary storage
  WRITE(iw,*) "READING IN PDB FILE ",setup%coord_file_name
  natom = 0
  label = 'REMARK PDB'
  CALL parser_init(setup%coord_file_name,globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
       CALL stop_program ( 'read_coordinate_pdb', 'REMARK PDB'  )
  ELSE
    CALL read_line
    DO WHILE (test_next() /= 'X')
      ilen = 10
      CALL get_next(string,ilen)
      CALL uppercase(string)
      SELECT CASE (string)
      CASE ('ATOM')
        !WRITE (iw,*) "Case ATOM"
        natom = natom + 1
        CALL get_next(num1(natom))
        CALL get_next(label1(natom),ilen)
        CALL get_next(label2(natom),ilen)
        CALL get_next(num2(natom))
        CALL get_next(pos1(natom))
        CALL get_next(pos2(natom))
        CALL get_next(pos3(natom))
        !WRITE(iw,*) num1(natom),label1(natom),label2(natom),num2(natom),&
        !            pos1(natom),pos2(natom),pos3(natom)
        !WRITE(iw,*) pos1(natom),pos2(natom),pos3(natom)
        !WRITE(iw,*) label1(natom)
      CASE ('HEATOM')
        !WRITE (iw,*) "Case HEATOM"
        CALL stop_program ('read_coordinate_pdb','HEATOM not implemented yet')
      CASE ('END')
        !WRITE (iw,*) "Case END"
        EXIT
      CASE DEFAULT
        !WRITE(iw,*) "PDB| ",string
      END SELECT
      CALL read_line
    END DO
  END IF
  CALL parser_end

  WRITE(iw,*) "NATOM = ",natom

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Put all the info from temporary storage info compatible FIST storage

  !Determine the number of kinds (nkind)
  nkind = 0
  CALL reallocate(kind_label,1,max_name_length,1,nblock_atom)
  CALL reallocate(kind_of,1,natom)
  CALL reallocate(natom_of_kind,1,natom)
  CALL reallocate(r,1,3,1,natom)
  CALL azero(kind_of,natom)
  CALL azero(natom_of_kind,natom)
  DO iatom=1,natom
    new_kind = .TRUE.
    DO ikind=1,nkind
      CALL integer_to_string(kind_label(:,ikind),string)
      IF(string==label1(iatom)) THEN
        new_kind = .FALSE.
        kind_of(iatom)=ikind
        natom_of_kind(ikind) = natom_of_kind(ikind)+1
        EXIT
      END IF
    END DO
    IF(new_kind) THEN
      nkind = nkind + 1
      kind_of(iatom)=nkind
      CALL string_to_integer(label1(iatom),kind_label(:,nkind))
      natom_of_kind(ikind) = natom_of_kind(ikind)+1
    END IF
    r(1,iatom) = pos1(iatom)
    r(2,iatom) = pos2(iatom)
    r(3,iatom) = pos3(iatom)
  END DO
  !WRITE(iw,*) "NKIND = ",nkind
  CALL reallocate(kind_label,1,max_name_length,1,nkind)
  CALL reallocate(natom_of_kind,1,nkind)
  CALL reallocate(charge,1,nkind)
  CALL reallocate(mass,1,nkind)

  !DO ikind=1,nkind
  !  WRITE(iw,*) "iatom = ",iatom," nkind = ",kind_of(iatom)
  !  WRITE(iw,*) "ikind = ",ikind," number = ",natom_of_kind(ikind)
  !END DO

  !Stuff the atom_list(1:nkind,1:natom_of_kind(ikind))
  !WRITE(iw,*) "MAXVAL(natom_of_kind)",MAXVAL(natom_of_kind)
  CALL reallocate(atom_list,1,nkind,1,MAXVAL(natom_of_kind))
  DO ikind=1,nkind
    iat = 0
    DO iatom = 1,natom
      IF (kind_of(iatom)==ikind) THEN
        iat = iat + 1
        atom_list(ikind,iat) = iatom
      END IF
    END DO
  END DO

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

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

!   *** check if the labels in &COORD match with the &MOLECULE setup
!   *** and store the charge and mass read in the &MOLECULE section

  ALLOCATE (natom_of_kind_check(nkind), STAT=istat) 
  natom_of_kind_check = 0
  nmol_type = SIZE (mol_setup) 
  natom_types = SIZE(atom_names)
  iat = 0
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
        IF (ikind/=n) CALL stop_program ( 'read_coordinate_pdb', &
             'kind label ordering in coord and mol dont match ')
        CALL integer_to_string(kind_label(:,ikind),string)
        IF (string/=mol_setup(i)%molpar%aname(j)) CALL stop_program ( &
             'read_coordinate_pdb', &
             'inconsistent kind labels in coord and mol')
        natom_of_kind_check(ikind) = natom_of_kind_check(ikind) + 1
        iatom =  atom_list(ikind,natom_of_kind_check(ikind)) 
        IF (iatom/=iat) CALL stop_program ( 'read_coordinate_pdb', &
             'atom lists in coord and mol dont match ')
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
  DEALLOCATE (atom_list,kind_label,kind_of,natom_of_kind,&
              natom_of_kind_check,r,STAT=istat)
  DEALLOCATE (charge,mass,STAT=istat)
  IF (istat /= 0) THEN
    CALL stop_memory('read_coordinate_pdb',&
              'release work storage')
  END IF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  DEALLOCATE(num1)
  DEALLOCATE(num2)
  DEALLOCATE(label1)
  DEALLOCATE(label2)
  DEALLOCATE(pos1)
  DEALLOCATE(pos2)
  DEALLOCATE(pos3)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  !CALL stop_program ('read_coordinate_pdb','checkpoint')

END SUBROUTINE read_coordinate_pdb

END MODULE topology_pdb
