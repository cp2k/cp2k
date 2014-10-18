!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Handles PDB files
!>
!> PDB Format Description Version 2.2 from http://www.rcsb.org
!> COLUMNS       DATA TYPE       FIELD         DEFINITION
!>
!>  1 -  6       Record name     "ATOM  "
!>  7 - 11       Integer         serial        Atom serial number.
!> 13 - 16       Atom            name          Atom name.
!> 17            Character       altLoc        Alternate location indicator.
!> 18 - 20       Residue name    resName       Residue name.
!> 22            Character       chainID       Chain identifier.
!> 23 - 26       Integer         resSeq        Residue sequence number.
!> 27            AChar           iCode         Code for insertion of residues.
!> 31 - 38       Real(8.3)       x             Orthogonal coordinates for X in
!>                                             Angstroms.
!> 39 - 46       Real(8.3)       y             Orthogonal coordinates for Y in
!>                                             Angstroms.
!> 47 - 54       Real(8.3)       z             Orthogonal coordinates for Z in
!>                                             Angstroms.
!> 55 - 60       Real(6.2)       occupancy     Occupancy.
!> 61 - 66       Real(6.2)       tempFactor    Temperature factor.
!> 73 - 76       LString(4)      segID         Segment identifier, left-justified.
!> 77 - 78       LString(2)      element       Element symbol, right-justified.
!> 79 - 80       LString(2)      charge        Charge on the atom.
!>
!> 81 -          Real(*)         Charge Ext.   This last field is an extenstion to
!>                                             standard PDB to provide a full charge
!>                                             without limitation of digits.
!>
!>  1 -  6       Record name    "CRYST1"
!>  7 - 15       Real(9.3)      a (Angstroms)
!> 16 - 24       Real(9.3)      b (Angstroms)
!> 25 - 33       Real(9.3)      c (Angstroms)
!> 34 - 40       Real(7.2)      alpha (degrees)
!> 41 - 47       Real(7.2)      beta (degrees)
!> 48 - 54       Real(7.2)      gamma (degrees)
!> 56 - 66       LString        Space group
!> 67 - 70       Integer        Z value
! *****************************************************************************
MODULE topology_pdb
  USE cell_types,                      ONLY: get_cell
  USE cp2k_info,                       ONLY: compile_revision,&
                                             cp2k_version,&
                                             r_datx,&
                                             r_host_name,&
                                             r_user_name
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_constants,                 ONLY: do_conn_user
  USE input_section_types,             ONLY: section_get_rval,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE physcon,                         ONLY: angstrom
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             topology_parameters_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_pdb'

  PRIVATE
  PUBLIC :: read_coordinate_pdb,write_coordinate_pdb

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param topology ...
!> \param para_env ...
!> \param subsys_section ...
!> \param error ...
!> \par History
!>      TLAINO 05.2004 - Added the TER option to use different non-bonded molecules
! *****************************************************************************
  SUBROUTINE read_coordinate_pdb (topology,para_env,subsys_section,error)
    TYPE(topology_parameters_type)           :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_pdb', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000

    CHARACTER(LEN=default_path_length)       :: line
    CHARACTER(LEN=default_string_length)     :: record, root_mol_name, strtmp
    INTEGER                                  :: handle, id0, inum_mol, istat, &
                                                iw, natom, newsize
    LOGICAL                                  :: my_end
    REAL(KIND=dp)                            :: pfactor
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    NULLIFY(parser, logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/PDB_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)

    pfactor = section_get_rval(subsys_section,"TOPOLOGY%MEMORY_PROGRESSION_FACTOR",error)
    atom_info => topology%atom_info
    CALL reallocate(atom_info%id_molname,1,nblock)
    CALL reallocate(atom_info%id_resname,1,nblock)
    CALL reallocate(atom_info%resid,1,nblock)
    CALL reallocate(atom_info%id_atmname,1,nblock)
    CALL reallocate(atom_info%r,1,3,1,nblock)
    CALL reallocate(atom_info%atm_mass,1,nblock)
    CALL reallocate(atom_info%atm_charge,1,nblock)
    CALL reallocate(atom_info%occup,1,nblock)
    CALL reallocate(atom_info%beta,1,nblock)
    CALL reallocate(atom_info%id_element,1,nblock)

    IF (iw > 0) THEN
       WRITE (UNIT=iw,FMT="(T2,A)")&
        "BEGIN of PDB data read from file "//TRIM(topology%coord_file_name)
    END IF

    id0 = str2id(s2s(""))
    topology%molname_generated = .FALSE.

    CALL parser_create(parser,topology%coord_file_name,para_env=para_env,error=error)

    natom         = 0
    inum_mol      = 1
    WRITE (UNIT=root_mol_name,FMT='(A3,I0)') "MOL",inum_mol
    DO
       line = ""
       CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
       IF(my_end) EXIT
       line   = parser%input_line(1:default_path_length)
       record = line(1:6)
       record = TRIM(record)

       IF ((record == "ATOM").OR.(record == "HETATM")) THEN
          natom = natom + 1
          topology%natoms = natom
          IF (natom > SIZE(atom_info%id_atmname)) THEN
             newsize = pfactor*natom
             CALL reallocate(atom_info%id_molname,1,newsize)
             CALL reallocate(atom_info%id_resname,1,newsize)
             CALL reallocate(atom_info%resid,1,newsize)
             CALL reallocate(atom_info%id_atmname,1,newsize)
             CALL reallocate(atom_info%r,1,3,1,newsize)
             CALL reallocate(atom_info%atm_mass,1,newsize)
             CALL reallocate(atom_info%atm_charge,1,newsize)
             CALL reallocate(atom_info%occup,1,newsize)
             CALL reallocate(atom_info%beta,1,newsize)
             CALL reallocate(atom_info%id_element,1,newsize)
          END IF
       END IF

       SELECT CASE (record)
       CASE ("ATOM","HETATM")
          READ (UNIT=line(13:16),FMT=*) strtmp
          atom_info%id_atmname(natom) = str2id(s2s(strtmp))
          READ (UNIT=line(18:20),FMT=*,IOSTAT=istat) strtmp
          IF (istat == 0) THEN
             atom_info%id_resname(natom) = str2id(s2s(strtmp))
          ELSE
             atom_info%id_resname(natom) = id0
          END IF
          READ (UNIT=line(23:26),FMT=*,IOSTAT=istat) atom_info%resid(natom)
          READ (UNIT=line(31:38),FMT=*,IOSTAT=istat) atom_info%r(1,natom)
          READ (UNIT=line(39:46),FMT=*,IOSTAT=istat) atom_info%r(2,natom)
          READ (UNIT=line(47:54),FMT=*,IOSTAT=istat) atom_info%r(3,natom)
          READ (UNIT=line(55:60),FMT=*,IOSTAT=istat) atom_info%occup(natom)
          READ (UNIT=line(61:66),FMT=*,IOSTAT=istat) atom_info%beta(natom)
          READ (UNIT=line(73:76),FMT=*,IOSTAT=istat) strtmp
          IF (istat == 0) THEN
             atom_info%id_molname(natom) = str2id(s2s(strtmp))
          ELSE
             atom_info%id_molname(natom) = str2id(s2s(root_mol_name))
             topology%molname_generated = .TRUE.
          END IF
          READ (UNIT=line(77:78),FMT=*,IOSTAT=istat) strtmp
          IF (istat == 0) THEN
             atom_info%id_element(natom) = str2id(s2s(strtmp))
          ELSE
             atom_info%id_element(natom) = id0
          END IF
          atom_info%atm_mass(natom)  = 0.0_dp
          atom_info%atm_charge(natom)= -HUGE(0.0_dp)
          IF (topology%charge_occup) atom_info%atm_charge(natom) = atom_info%occup(natom)
          IF (topology%charge_beta)  atom_info%atm_charge(natom) = atom_info%beta(natom)
          IF (topology%charge_extended) THEN
             READ (UNIT=line(81:),FMT=*,IOSTAT=istat) atom_info%atm_charge(natom)
          END IF

          IF (atom_info%id_element(natom) == id0) THEN
             ! Element is assigned on the basis of the atm_name
             topology%aa_element = .TRUE.
             atom_info%id_element(natom) = atom_info%id_atmname(natom)
          END IF

          IF (iw > 0) THEN
             WRITE (UNIT=iw,FMT="(A6,I5,T13,A4,T18,A3,T23,I4,T31,3F8.3,T73,A4,T77,A2)")&
               record,natom,&
               TRIM(id2str(atom_info%id_atmname(natom))),&
               TRIM(id2str(atom_info%id_resname(natom))),&
               atom_info%resid(natom),&
               atom_info%r(1,natom),&
               atom_info%r(2,natom),&
               atom_info%r(3,natom),&
               ADJUSTL(TRIM(id2str(atom_info%id_molname(natom)))),&
               ADJUSTR(TRIM(id2str(atom_info%id_element(natom))))
          END IF
          atom_info%r(1,natom) = cp_unit_to_cp2k(atom_info%r(1,natom),"angstrom",error=error)
          atom_info%r(2,natom) = cp_unit_to_cp2k(atom_info%r(2,natom),"angstrom",error=error)
          atom_info%r(3,natom) = cp_unit_to_cp2k(atom_info%r(3,natom),"angstrom",error=error)
       CASE ("TER")
          inum_mol = inum_mol + 1
          WRITE (UNIT=root_mol_name,FMT='(A3,I0)') "MOL",inum_mol
       CASE ("REMARK")
          IF (iw > 0) WRITE (UNIT=iw,FMT=*) TRIM(line)
       CASE ("END")
          EXIT
       CASE DEFAULT
       END SELECT
    END DO
    CALL parser_release(parser,error=error)

    CALL reallocate(atom_info%id_molname,1,natom)
    CALL reallocate(atom_info%id_resname,1,natom)
    CALL reallocate(atom_info%resid,1,natom)
    CALL reallocate(atom_info%id_atmname,1,natom)
    CALL reallocate(atom_info%r,1,3,1,natom)
    CALL reallocate(atom_info%atm_mass,1,natom)
    CALL reallocate(atom_info%atm_charge,1,natom)
    CALL reallocate(atom_info%occup,1,natom)
    CALL reallocate(atom_info%beta,1,natom)
    CALL reallocate(atom_info%id_element,1,natom)

    IF (topology%conn_type /= do_conn_user) THEN
       IF (.NOT.topology%para_res) atom_info%resid(:) = 1
    END IF

    IF (iw > 0) THEN
       WRITE (UNIT=iw,FMT="(T2,A)")&
         "END of PDB data read from file "//TRIM(topology%coord_file_name)
    END IF

    topology%natoms = natom
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/PDB_INFO",error=error)
    CALL timestop(handle)

  END SUBROUTINE read_coordinate_pdb

! *****************************************************************************
!> \brief ...
!> \param file_unit ...
!> \param topology ...
!> \param subsys_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE write_coordinate_pdb (file_unit,topology,subsys_section,error)

    INTEGER, INTENT(IN)                      :: file_unit
    TYPE(topology_parameters_type)           :: topology
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_coordinate_pdb', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=120)                       :: line
    CHARACTER(LEN=default_string_length)     :: my_tag1, my_tag2, my_tag3, &
                                                my_tag4, record
    INTEGER                                  :: handle, i, id1, id2, idres, &
                                                iw, natom
    LOGICAL                                  :: charge_beta, charge_extended, &
                                                charge_occup, ldum
    REAL(KIND=dp)                            :: angle_alpha, angle_beta, &
                                                angle_gamma
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    NULLIFY (logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/PDB_INFO",&
         extension=".subsysLog",error=error)
    print_key => section_vals_get_subs_vals(subsys_section,"TOPOLOGY%DUMP_PDB",error=error)
    CALL timeset(routineN,handle)

    CALL section_vals_val_get(print_key,"CHARGE_OCCUP",l_val=charge_occup,error=error)
    CALL section_vals_val_get(print_key,"CHARGE_BETA",l_val=charge_beta,error=error)
    CALL section_vals_val_get(print_key,"CHARGE_EXTENDED",l_val=charge_extended,error=error)
    i = COUNT((/charge_occup,charge_beta,charge_extended/))
    CALL cp_assert((i <= 1),cp_failure_level,cp_assertion_failed,routineP,&
                   "Either only CHARGE_OCCUP, CHARGE_BETA, or CHARGE_EXTENDED can be selected, "//&
CPSourceFileRef,&
                   only_ionode=.TRUE.)

    atom_info => topology%atom_info
    record = cp_print_key_generate_filename(logger,print_key,&
                                            extension=".pdb",&
                                            my_local=.FALSE.,&
                                            error=error)

    IF (iw > 0) WRITE (UNIT=iw,FMT=*) "    Writing out PDB file ",TRIM(record)

    ! Write file header
    WRITE (UNIT=file_unit,FMT="(A6,T11,A)")&
     "TITLE ","PDB file created by "//TRIM(cp2k_version)//" (revision "//TRIM(compile_revision)//")",&
     "AUTHOR",TRIM(r_user_name)//"@"//TRIM(r_host_name)//" "//r_datx(1:19)
    ! Write cell information
    CALL get_cell(cell=topology%cell,alpha=angle_alpha,beta=angle_beta,gamma=angle_gamma,abc=abc)
    WRITE (UNIT=file_unit,FMT="(A6,3F9.3,3F7.2)")&
     "CRYST1",abc(1:3)*angstrom,angle_alpha,angle_beta,angle_gamma

    natom = topology%natoms
    idres = 0
    id1   = 0
    id2   = 0

    DO i=1,natom

       IF (topology%para_res) THEN
          idres = atom_info%resid(i)
       ELSE
          IF((id1 /= atom_info%map_mol_num(i)).OR.(id2 /= atom_info%map_mol_typ(i))) THEN
             idres = idres + 1
             id1 = atom_info%map_mol_num(i)
             id2 = atom_info%map_mol_typ(i)
          END IF
       END IF

       line = ""
       my_tag1 = id2str(atom_info%id_atmname(i)); ldum=qmmm_ff_precond_only_qm(my_tag1)
       my_tag2 = id2str(atom_info%id_resname(i)); ldum=qmmm_ff_precond_only_qm(my_tag2)
       my_tag3 = id2str(atom_info%id_molname(i)); ldum=qmmm_ff_precond_only_qm(my_tag3)
       my_tag4 = id2str(atom_info%id_element(i)); ldum=qmmm_ff_precond_only_qm(my_tag4)

       WRITE (UNIT=line( 1: 6),FMT="(A6)") "ATOM  "
       WRITE (UNIT=line( 7:11),FMT="(I5)") MODULO(i,100000)
       WRITE (UNIT=line(13:16),FMT="(A4)") ADJUSTL(my_tag1(1:4))
       WRITE (UNIT=line(18:20),FMT="(A3)") TRIM(my_tag2)
       WRITE (UNIT=line(23:26),FMT="(I4)") MODULO(idres,10000)
       WRITE (UNIT=line(31:54),FMT="(3F8.3)") atom_info%r(1:3,i)*angstrom
       IF (ASSOCIATED(atom_info%occup)) THEN
          WRITE (UNIT=line(55:60),FMT="(F6.2)") atom_info%occup(i)
       ELSE
          WRITE (UNIT=line(55:60),FMT="(F6.2)") 0.0_dp
       END IF
       IF (ASSOCIATED(atom_info%beta)) THEN
           WRITE (UNIT=line(61:66),FMT="(F6.2)") atom_info%beta(i)
       ELSE
           WRITE (UNIT=line(61:66),FMT="(F6.2)") 0.0_dp
       END IF
       IF (ASSOCIATED(atom_info%atm_charge)) THEN
          CALL cp_assert(.NOT.(ANY((/charge_occup,charge_beta,charge_extended/)).AND.&
                         (atom_info%atm_charge(i) == -HUGE(0.0_dp))),&
                         cp_failure_level,cp_assertion_failed,routineP,&
                         "No atomic charges found yet (after the topology setup), "//&
CPSourceFileRef,&
                         only_ionode=.TRUE.)
          IF (charge_occup) THEN
             WRITE (UNIT=line(55:60),FMT="(F6.2)") atom_info%atm_charge(i)
          ELSE IF (charge_beta) THEN
             WRITE (UNIT=line(61:66),FMT="(F6.2)") atom_info%atm_charge(i)
          ELSE IF (charge_extended) THEN
             WRITE (UNIT=line(81:),FMT="(F20.16)") atom_info%atm_charge(i)
          ELSE
             ! Wrtie no atomic charge
          END IF
       END IF
       WRITE (UNIT=line(73:76),FMT="(A4)") ADJUSTL(my_tag3)
       WRITE (UNIT=line(77:78),FMT="(A2)") TRIM(my_tag4)
       WRITE (UNIT=file_unit,FMT="(A)") TRIM(line)
    END DO
    WRITE (UNIT=file_unit,FMT="(A3)") "END"

    IF (iw > 0) WRITE (UNIT=iw,FMT=*) "  Exiting "//routineN

    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/PDB_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE write_coordinate_pdb

END MODULE topology_pdb
