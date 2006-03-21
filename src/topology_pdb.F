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

  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_parser,                       ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_get_next_line,&
                                             parser_release
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             topology_parameters_type
#include "cp_common_uses.h"  

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_pdb'

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

  SUBROUTINE read_coordinate_pdb (topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type)           :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_pdb', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock_atom = 1000

    CHARACTER(LEN=default_string_length)     :: line, record, root_mol_name
    INTEGER                                  :: handle, i, inum_mol, istat, &
                                                iw, natom
    LOGICAL                                  :: my_end
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    NULLIFY(parser, logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/PDB_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)

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

    IF (iw>0)  WRITE(iw,*) "    Reading in PDB file ",TRIM(topology%coord_file_name)
    CALL parser_create(parser,topology%coord_file_name,&
         para_env=globenv%para_env,error=error)  

    natom         = 0
    inum_mol      = 1
    WRITE(root_mol_name,'(A3,I7)')"MOL",inum_mol
    DO i =  1, 10
       IF (root_mol_name(i:i) == " ") root_mol_name(i:i) = "0" 
    END DO
    
    DO
       line = ''
       CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
       IF(my_end) EXIT
       line   = parser%input_line(1:default_string_length)
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
          atom_info%atm_mass(natom)  = 0.0_dp
          atom_info%atm_charge(natom)= 0.0_dp
          IF(topology%charge_occup) & 
               atom_info%atm_charge(natom) = atom_info%occup(natom)
          IF(topology%charge_beta) &
               atom_info%atm_charge(natom) = atom_info%beta(natom)
          IF(iw>0) WRITE(iw,*) natom,line(13:16),line(77:78)
          IF(iw>0) WRITE(iw,'(5X,A4,5X,A2,3F8.3,A4,I4)') &
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
             IF (root_mol_name(i:i) == " ") root_mol_name(i:i) = "0" 
          END DO
       CASE ("REMARK")
          IF(iw>0) WRITE(iw,*) TRIM(line)
       CASE ("END")
          EXIT
       CASE DEFAULT
       END SELECT
    END DO
    CALL parser_release(parser,error=error)
    
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
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/PDB_INFO",error=error)
    CALL timestop(0.0_dp,handle)

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

  SUBROUTINE write_coordinate_pdb (topology,globenv,subsys_section,error)
  
    TYPE(topology_parameters_type)           :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_coordinate_pdb', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: line, my_tag1, my_tag2, &
                                                my_tag3, my_tag4, record
    INTEGER                                  :: file_unit, handle, i, id1, &
                                                id2, idres, iw, natom
    LOGICAL                                  :: ldum
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/PDB_INFO",&
         extension=".subsysLog",error=error)    
    
    CALL timeset(routineN,'I','',handle)

    atom_info => topology%atom_info
    record = "topology_dump.pdb"

    IF(iw>0) WRITE(iw,*) "    Writing out PDB file ",TRIM(record)

    CALL open_file(file_name=record,&
         file_status="REPLACE",&
         file_action="WRITE",&
         unit_number=file_unit)
    
    natom = topology%natoms
    idres = 0
    id1   = 0
    id2   = 0
    DO i=1,natom

       IF (topology%para_res) THEN
          idres = atom_info%label_resid(i)
       ELSE
          IF((id1/=atom_info%map_mol_num(i)).OR.(id2/=atom_info%map_mol_typ(i))) THEN
             idres = idres + 1
             id1 = atom_info%map_mol_num(i)
             id2 = atom_info%map_mol_typ(i)
          END IF
       END IF

       line = ""
       my_tag1 = atom_info%label_atmname(i); ldum=qmmm_ff_precond_only_qm(my_tag1)
       my_tag2 = atom_info%label_resname(i); ldum=qmmm_ff_precond_only_qm(my_tag2)
       my_tag3 = atom_info%label_molname(i); ldum=qmmm_ff_precond_only_qm(my_tag3)
       my_tag4 = atom_info%element(i)      ; ldum=qmmm_ff_precond_only_qm(my_tag4)
       WRITE(line(1 : 6),'(A6)')"ATOM  "
       WRITE(line(7 :11),'(I5)')i
       WRITE(line(13:16),'(A4)')TRIM(my_tag1)
       WRITE(line(18:20),'(A3)')TRIM(my_tag2)
       WRITE(line(23:26),'(I4)')idres
       WRITE(line(31:38),'(F8.3)')atom_info%r(1,i)
       WRITE(line(39:46),'(F8.3)')atom_info%r(2,i)
       WRITE(line(47:54),'(F8.3)')atom_info%r(3,i)
       WRITE(line(55:60),'(F6.2)')atom_info%occup(i)
       WRITE(line(61:66),'(F6.2)')atom_info%beta(i)
       WRITE(line(73:76),'(A4)')TRIM(my_tag3)
       WRITE(line(77:78),'(A2)')TRIM(my_tag4)

       WRITE(file_unit,"(A)")line
    END DO
    CALL close_file(unit_number=file_unit)

    IF(iw>0) WRITE(iw,*) "  Exiting write_coordinate_pdb"

    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/PDB_INFO",error=error)

  END SUBROUTINE write_coordinate_pdb

END MODULE topology_pdb
