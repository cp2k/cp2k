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
!!     IKUO   08.16.03
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
MODULE topology_pdb

  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             default_string_length
  USE memory_utilities,                ONLY: reallocate
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: topology_parameters_type

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology_pdb"

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
!!     IKUO   08.16.03
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*** *************************************************************************
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
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE read_coordinate_pdb (topology,globenv)


    TYPE(topology_parameters_type)           :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_pdb', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock_atom = 1000

    CHARACTER(LEN=default_string_length)     :: line, record
    INTEGER                                  :: handle, istat,&
                                                iw, natom
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

  NULLIFY(topology%label_resname)
  NULLIFY(topology%label_atmname)
  NULLIFY(topology%r)
  NULLIFY(topology%atm_mass)
  NULLIFY(topology%atm_charge)

  CALL reallocate(topology%label_resname,1,nblock_atom)
  CALL reallocate(topology%label_atmname,1,nblock_atom)
  CALL reallocate(topology%r,1,3,1,nblock_atom)
  CALL reallocate(topology%atm_mass,1,nblock_atom)
  CALL reallocate(topology%atm_charge,1,nblock_atom)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  IF(PM1) WRITE(iw,*) "    Reading in PDB file ",TRIM(topology%coord_file_name)
  CALL open_file(file_name=topology%coord_file_name,&
                 file_action="READ",&
                 unit_number=file_unit)

  natom =  0

  DO
    READ(file_unit,'(A80)',IOSTAT=istat) line 
    IF(istat /= 0) EXIT

    record = line(1:6)
    record = TRIM(record)

    IF(record == "ATOM") THEN
      natom = natom + 1
      topology%natoms = natom
      IF(natom > SIZE(topology%label_resname)) THEN
        CALL reallocate(topology%label_resname,1,natom+nblock_atom)
        CALL reallocate(topology%label_atmname,1,natom+nblock_atom)
        CALL reallocate(topology%r,1,3,1,natom+nblock_atom)
        CALL reallocate(topology%atm_mass,1,natom+nblock_atom)
        CALL reallocate(topology%atm_charge,1,natom+nblock_atom)
      END IF
    END IF

    SELECT CASE (record)
    CASE ("ATOM") 
      READ(line(13:16),*) topology%label_atmname(natom)
      READ(line(18:20),*) topology%label_resname(natom)
      READ(line(31:38),*) topology%r(1,natom)
      READ(line(39:46),*) topology%r(2,natom)
      READ(line(47:54),*) topology%r(3,natom)
      topology%atm_mass(natom)  =0.0_dbl
      topology%atm_charge(natom)=0.0_dbl
      IF(PM2) WRITE(iw,'(5X,A4,3F8.3,A4)') &
                   TRIM(topology%label_atmname(natom)),&
                   topology%r(1,natom),&
                   topology%r(2,natom),&
                   topology%r(3,natom),&
                   TRIM(topology%label_resname(natom))
      CALL convert_to_cp2k_units("ANGSTROM",length=topology%r(1,natom))
      CALL convert_to_cp2k_units("ANGSTROM",length=topology%r(2,natom))
      CALL convert_to_cp2k_units("ANGSTROM",length=topology%r(3,natom))
    CASE ("END")
      EXIT
    CASE DEFAULT
    END SELECT
  END DO
  CALL close_file(unit_number=file_unit)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  IF(PM1) WRITE(iw,*) "  Exiting read_coordinate_pdb"

  CALL timestop(0.0_dbl,handle)
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("exiting "//routineN,globenv)
  END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END SUBROUTINE read_coordinate_pdb

END MODULE topology_pdb
