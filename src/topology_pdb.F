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
!PDB Format Description Version 2.2 from http://www.rcsb.org
!COLUMNS       DATA TYPE       FIELD         DEFINITION
!-------------------------------------------------------------------------------
! 1 -  6       Record name     "ATOM  "
! 7 - 11       Integer         serial        Atom serial number.
!13 - 16       Atom            name          Atom name.
!17            Character       altLoc        Alternate location indicator.
!18 - 20       Residue name    resName       Residue name.
!22            Character       chainID       Chain identifier.
!23 - 26       Integer         resSeq        Residue sequence number.
!27            AChar           iCode         Code for insertion of residues.
!31 - 38       Real(8.3)       x             Orthogonal coordinates for X in
!                                            Angstroms.
!39 - 46       Real(8.3)       y             Orthogonal coordinates for Y in
!                                            Angstroms.
!47 - 54       Real(8.3)       z             Orthogonal coordinates for Z in
!                                            Angstroms.
!55 - 60       Real(6.2)       occupancy     Occupancy.
!61 - 66       Real(6.2)       tempFactor    Temperature factor.
!73 - 76       LString(4)      segID         Segment identifier, left-justified.
!77 - 78       LString(2)      element       Element symbol, right-justified.
!79 - 80       LString(2)      charge        Charge on the atom.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
MODULE topology_pdb

  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE memory_utilities,                ONLY: reallocate
  USE physcon,                         ONLY: bohr
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: topology_parameters_type,&
                                             atom_info_type

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology_pdb"

  PRIVATE
  PUBLIC :: read_coordinate_pdb,write_coordinate_pdb

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
!!     TLAINO 05.2004 - Added the TER option to use different non-bonded molecules
!!                         
!!
!!   OPTIONS
!!
!!   INPUTS
!!*** *************************************************************************
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE read_coordinate_pdb (topology,globenv)


    TYPE(topology_parameters_type)           :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    TYPE (atom_info_type),POINTER            :: atom_info
    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_pdb', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock_atom = 1000

    CHARACTER(LEN=default_string_length)     :: line, record, root_mol_name
    INTEGER                                  :: handle, istat,&
                                                iw, natom,i, inum_mol
    INTEGER                                  :: file_unit

    LOGICAL                                  :: PM1, PM2, PM3, PM4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "  Entering read_coordinate_pdb"

  CALL write_checkpoint_information("entering "//routineN,globenv)
  CALL timeset(routineN,'I','',handle)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  atom_info => topology%atom_info

  CALL reallocate(atom_info%label_molname,1,nblock_atom)
  CALL reallocate(atom_info%label_resname,1,nblock_atom)
  CALL reallocate(atom_info%label_resid,1,nblock_atom)
  CALL reallocate(atom_info%label_atmname,1,nblock_atom)
  CALL reallocate(atom_info%r,1,3,1,nblock_atom)
  CALL reallocate(atom_info%atm_mass,1,nblock_atom)
  CALL reallocate(atom_info%atm_charge,1,nblock_atom)
  CALL reallocate(atom_info%occup,1,nblock_atom)
  CALL reallocate(atom_info%beta,1,nblock_atom)
  CALL reallocate(atom_info%element,1,nblock_atom)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  IF(PM1) WRITE(iw,*) "    Reading in PDB file ",TRIM(topology%coord_file_name)
  CALL open_file(file_name=topology%coord_file_name,&
                 file_action="READ",&
                 unit_number=file_unit)

  natom         = 0
  inum_mol      = 1
  WRITE(root_mol_name,'(A3,I7)')"MOL",inum_mol
  DO i =  1, 10
     if (root_mol_name(i:i) == " ") root_mol_name(i:i) = "0" 
  END DO

  DO
    line = ''
    READ(file_unit,'(A80)',IOSTAT=istat) line 
    IF(istat /= 0) EXIT

    record = line(1:6)
    record = TRIM(record)

    IF(record == "ATOM" .OR. record == "HETATM") THEN
      natom = natom + 1
      topology%natoms = natom
      IF(natom > SIZE(atom_info%label_atmname)) THEN
        CALL reallocate(atom_info%label_molname,1,natom+nblock_atom)
        CALL reallocate(atom_info%label_resname,1,natom+nblock_atom)
        CALL reallocate(atom_info%label_resid,1,natom+nblock_atom)
        CALL reallocate(atom_info%label_atmname,1,natom+nblock_atom)
        CALL reallocate(atom_info%r,1,3,1,natom+nblock_atom)
        CALL reallocate(atom_info%atm_mass,1,natom+nblock_atom)
        CALL reallocate(atom_info%atm_charge,1,natom+nblock_atom)
        CALL reallocate(atom_info%occup,1,natom+nblock_atom)
        CALL reallocate(atom_info%beta,1,natom+nblock_atom)
        CALL reallocate(atom_info%element,1,natom+nblock_atom)
      END IF
    END IF

    SELECT CASE (record)
    CASE ("ATOM","HETATM") 
      atom_info%label_molname(natom) = TRIM(root_mol_name)
      READ(line(13:16),*) atom_info%label_atmname(natom)
      CALL uppercase(atom_info%label_atmname(natom))
      READ(line(18:20),*,IOSTAT=istat) atom_info%label_resname(natom)
      READ(line(23:26),*,IOSTAT=istat) atom_info%label_resid(natom)
      READ(line(31:38),*,IOSTAT=istat) atom_info%r(1,natom)
      READ(line(39:46),*,IOSTAT=istat) atom_info%r(2,natom)
      READ(line(47:54),*,IOSTAT=istat) atom_info%r(3,natom)
      READ(line(55:60),*,IOSTAT=istat) atom_info%occup(natom)
      READ(line(61:66),*,IOSTAT=istat) atom_info%beta(natom)
      READ(line(73:76),*,IOSTAT=istat) atom_info%label_molname(natom)
      READ(line(77:78),*,IOSTAT=istat) atom_info%element(natom)
      IF(LEN_TRIM(atom_info%element(natom)) == 0 )  atom_info%element(natom) = atom_info%label_atmname(natom)
      atom_info%atm_mass(natom)  =0.0_dp
      atom_info%atm_charge(natom)=-10000.0_dp
      IF(topology%charge_occup) & 
        atom_info%atm_charge(natom) = atom_info%occup(natom)
      IF(topology%charge_beta) &
        atom_info%atm_charge(natom) = atom_info%beta(natom)
      IF(PM2) WRITE(iw,*) natom,line(13:16),line(77:78)
      IF(PM2) WRITE(iw,'(5X,A4,5X,A2,3F8.3,A4,I4)') &
                   TRIM(atom_info%label_atmname(natom)),&
                   TRIM(atom_info%element(natom)),&
                   atom_info%r(1,natom),&
                   atom_info%r(2,natom),&
                   atom_info%r(3,natom),&
                   TRIM(atom_info%label_resname(natom)),&
                   atom_info%label_resid(natom)
      CALL convert_to_cp2k_units("ANGSTROM",length=atom_info%r(1,natom))
      CALL convert_to_cp2k_units("ANGSTROM",length=atom_info%r(2,natom))
      CALL convert_to_cp2k_units("ANGSTROM",length=atom_info%r(3,natom))
    CASE ("TER")
      inum_mol = inum_mol + 1
      WRITE(root_mol_name,'(A3,I7)')"MOL",inum_mol
      DO i =  1, 10
         if (root_mol_name(i:i) == " ") root_mol_name(i:i) = "0" 
      END DO
    CASE ("REMARK")
      IF(PM1) WRITE(iw,*) TRIM(line)
    CASE ("END")
      EXIT
    CASE DEFAULT
    END SELECT
  END DO
  CALL close_file(unit_number=file_unit)

  CALL reallocate(atom_info%label_molname,1,natom)
  CALL reallocate(atom_info%label_resname,1,natom)
  CALL reallocate(atom_info%label_resid,1,natom)
  CALL reallocate(atom_info%label_atmname,1,natom)
  CALL reallocate(atom_info%r,1,3,1,natom)
  CALL reallocate(atom_info%atm_mass,1,natom)
  CALL reallocate(atom_info%atm_charge,1,natom)
  CALL reallocate(atom_info%occup,1,natom)
  CALL reallocate(atom_info%beta,1,natom)
  CALL reallocate(atom_info%element,1,natom)

  IF(.NOT.topology%para_res) atom_info%label_resid(:) = 1

  topology%natoms = natom


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting read_coordinate_pdb"

  CALL timestop(0.0_dp,handle)
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("exiting "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END SUBROUTINE read_coordinate_pdb


!******************************************************************************
!!****** topology_pdb/write_coordinate_pdb [1.0] *
!!
!!   NAME
!!     write_coordinate_pdb
!!
!!   FUNCTION
!!
!!*** *************************************************************************
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE write_coordinate_pdb (topology,globenv)


    TYPE(topology_parameters_type)           :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    TYPE (atom_info_type),POINTER            :: atom_info
    CHARACTER(len=*), PARAMETER :: routineN = 'write_coordinate_pdb', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: line, record
    INTEGER                                  :: handle, istat,&
                                                iw, natom,i
    INTEGER                                  :: file_unit,imol,id1,id2

    LOGICAL                                  :: PM1, PM2, PM3, PM4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "  Entering write_coordinate_pdb"

  CALL write_checkpoint_information("entering "//routineN,globenv)
  CALL timeset(routineN,'I','',handle)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  atom_info => topology%atom_info

  record = "topology_dump.pdb"

  IF(PM1) WRITE(iw,*) "    Writing out PDB file ",TRIM(record)

  CALL open_file(file_name=record,&
                 file_status="REPLACE",&
                 file_action="WRITE",&
                 unit_number=file_unit)

  !line = ''
  !line(1:6) = "TITLE  "
  !line(11:70) = "CP2K generated DUMP of coordinates"
  !WRITE(file_unit,*) TRIM(line)

  natom = topology%natoms
  imol = 0
  id1 = 0 
  id2 = 0

  DO i=1,natom
    IF((id1/=atom_info%map_mol_num(i)).OR.(id2/=atom_info%map_mol_typ(i))) THEN
      imol = imol + 1
      id1 = atom_info%map_mol_num(i)
      id2 = atom_info%map_mol_typ(i)
    END IF
    
    WRITE(file_unit,"(A6,I5,1X,A4,1X,A3,2X,I4,4X,3F8.3,2F6.2)") &
                       "ATOM  ",i,&
                       TRIM(atom_info%label_atmname(i)),&
                       TRIM(atom_info%label_resname(i)),&
                       imol,&
                       (atom_info%r(1,i)/bohr),&
                       (atom_info%r(2,i)/bohr),&
                       (atom_info%r(3,i)/bohr),&
                       1.0,1.0
  END DO
  CALL close_file(unit_number=file_unit)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  IF(PM1) WRITE(iw,*) "  Exiting write_coordinate_pdb"

  CALL timestop(0.0_dp,handle)
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("exiting "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END SUBROUTINE write_coordinate_pdb

END MODULE topology_pdb
