!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****** cp2k/topology_gromos [1.0] *
!!
!!   NAME
!!     topology_gromos
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE topology_gromos
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE cp_files,                        ONLY: open_file, close_file
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_parser,                       ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_get_next_line,&
                                             parser_get_object,&
                                             parser_release,&
                                             parser_search_string
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
                                             connectivity_info_type,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: reorder_structure,&
                                             atom_bond_list_type
  USE input_constants,                 ONLY: do_skip_14

#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_gromos'

  PRIVATE
  PUBLIC :: pre_read_topology_gromos, read_topology_gromos, read_coordinate_g96, &
            write_topology_gromos

!!*****
!******************************************************************************

CONTAINS

!!*****
!******************************************************************************
!!****** topology/pre_read_topology_gromos [1.0] *
!!
!!   NAME
!!     pre_read_topology_gromos
!!
!!   FUNCTION
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE pre_read_topology_gromos (topology,subsys_section, error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pre_read_topology_gromos', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, stat
    LOGICAL                                  :: failure
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GTOP_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)
    
    atom_info => topology%atom_info

    IF (ASSOCIATED(atom_info%label_molname)) &
      DEALLOCATE(atom_info%label_molname,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(atom_info%label_molname)
  
    IF (ASSOCIATED(atom_info%label_resid)) &
      DEALLOCATE(atom_info%label_resid,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(atom_info%label_resid)
  
    IF (ASSOCIATED(atom_info%label_resname)) &
      DEALLOCATE(atom_info%label_resname,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(atom_info%label_resname)
    
    IF (ASSOCIATED(atom_info%label_atmname)) &
      DEALLOCATE(atom_info%label_atmname,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(atom_info%label_atmname)
  
    IF (ASSOCIATED(atom_info%atm_charge)) &
      DEALLOCATE(atom_info%atm_charge,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(atom_info%atm_charge)
  
    IF (ASSOCIATED(atom_info%atm_mass)) &
      DEALLOCATE(atom_info%atm_mass,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(atom_info%atm_mass)

    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/GTOP_INFO",error=error)

  END SUBROUTINE pre_read_topology_gromos

!!*****
!******************************************************************************
!!****** topology/read_topology_gromos [1.0] *
!!
!!   NAME
!!     read_topology_gromos
!!
!!   FUNCTION
!!     Read GROMOS topology file
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE read_topology_gromos (file_name,topology,globenv,subsys_section, error)
    CHARACTER(LEN=*), INTENT(IN)             :: file_name
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_topology_gromos', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: dummy_field, label,string
    INTEGER :: handle, i, iatom, ibond, index_now, iphi, itheta, iw, jatom, &
      natom, natom_prev, nbond, nbond_prev, nphi, nphi_prev, ntheta, &
      ntheta_prev, nonfo, ionfo, nonfo_prev, stat, N

    INTEGER :: itype,ntype

    LOGICAL                                  :: failure, found
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(atom_bond_list_type), DIMENSION(:), POINTER   :: ex_bond_list,&
                                                          ex_bend_list
    LOGICAL                                  :: do_onfo

    CHARACTER(LEN=default_string_length), &
      DIMENSION(20)                          :: avail_section
    CHARACTER(LEN=default_string_length), &
      POINTER                                :: namearray1(:),namearray2(:)
    INTEGER:: irestyp,iatyp,iresid,itemp,offset

    NULLIFY(parser, logger)
    NULLIFY(namearray1,namearray2)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GTOP_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)

    avail_section( 1) = "TITLE" 
    avail_section( 2) = "TOPPHYSCON"
    avail_section( 3) = "TOPVERSION"
    avail_section( 4) = "ATOMTYPENAME"
    avail_section( 5) = "RESNAME"
    avail_section( 6) = "SOLUTEATOM"
    avail_section( 7) = "BONDTYPE"
    avail_section( 8) = "BONDH"
    avail_section( 9) = "BOND"
    avail_section(10) = "BONDANGLETYPE"
    avail_section(11) = "BONDANGLEH" 
    avail_section(12) = "BONDANGLE"
    avail_section(13) = "IMPDIHEDRALTYPE"
    avail_section(14) = "IMPDIHEDRALH"
    avail_section(15) = "IMPDIHEDRAL"
    avail_section(16) = "DIHEDRALTYPE"
    avail_section(17) = "DIHEDRALH"
    avail_section(18) = "DIHEDRAL"
    avail_section(19) = "LJPARAMETERS"
    avail_section(20) = "SOLVENTATOM"

    atom_info => topology%atom_info
    conn_info => topology%conn_info
    do_onfo = (topology%exclude_vdw/=do_skip_14).OR.(topology%exclude_ei/=do_skip_14)

    natom_prev = 0
    IF(ASSOCIATED(atom_info%label_molname)) natom_prev = SIZE(atom_info%label_molname)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TITLE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the TITLE section'    
    CALL parser_create(parser,file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(1))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      DO 
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,string,string_length=80)
        IF(string == TRIM("END")) EXIT
        IF(iw>0) WRITE(iw,*) "GTOP_INFO| ",TRIM(string)
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TOPPHYSCON SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the TOPPHYSCON section'    
    CALL parser_create(parser,file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(2))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      DO 
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,string)
        IF(string == TRIM("END")) EXIT
        IF(iw>0) WRITE(iw,*) "GTOP_INFO|  ",TRIM(string)
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TOPVERSION SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the TOPVERSION section'    
    CALL parser_create(parser,file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(3))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      DO 
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,string)
        IF(string == TRIM("END")) EXIT
        IF(iw>0) WRITE(iw,*) "GTOP_INFO|  ",TRIM(string)
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ATOMTYPENAME SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the ATOMTYPENAME section'    
    CALL parser_create(parser,file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(4))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,ntype)
      CALL reallocate(namearray1,1,ntype)
      !ALLOCATE(namearray1(ntype),stat=stat)
      !CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO itype=1,ntype
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,namearray1(itype))
        IF(iw>0) WRITE(iw,*) "GTOP_INFO|  ",TRIM(namearray1(itype))
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! RESNAME SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the RESNAME section'    
    CALL parser_create(parser,file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(5))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,ntype)
      CALL reallocate(namearray2,1,ntype)
      !ALLOCATE(namearray2(ntype),stat=stat)
      !CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO itype=1,ntype
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,namearray2(itype))
        IF(iw>0) WRITE(iw,*) "GTOP_INFO|  ",TRIM(namearray2(itype))
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SOLUTEATOM SECTION
    iresid=1
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the SOLUTEATOM section'    
    CALL parser_create(parser,file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(6))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,natom)
      CALL reallocate(atom_info%label_molname,1,natom_prev+natom)
      CALL reallocate(atom_info%label_resid,1,natom_prev+natom)
      CALL reallocate(atom_info%label_resname,1,natom_prev+natom)
      CALL reallocate(atom_info%label_atmname,1,natom_prev+natom)
      CALL reallocate(atom_info%atm_charge,1,natom_prev+natom)
      CALL reallocate(atom_info%atm_mass,1,natom_prev+natom)
      DO iatom=1,natom
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,itemp)
        CALL parser_get_object(parser,itemp)
        atom_info%label_molname(iatom) = namearray2(itemp)
        atom_info%label_resname(iatom) = namearray2(itemp)
        CALL parser_get_object(parser,string)
        CALL parser_get_object(parser,itemp)
        atom_info%label_atmname(iatom) = namearray1(itemp)

        CALL parser_get_object(parser,atom_info%atm_mass(iatom))
        CALL parser_get_object(parser,atom_info%atm_charge(iatom))
        CALL parser_get_object(parser,itemp)
        atom_info%label_resid(iatom)=iresid
        IF(itemp==1) iresid=iresid+1
        IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT SOLUTEATOM INFO HERE!!!!"
        CALL parser_get_next_line(parser,1,error=error)
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SOLVENTATOM SECTION
    iresid=1
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the SOLVENTATOM section'    
    CALL parser_create(parser,file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(20))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,natom)
      offset = 0
      IF(ASSOCIATED(atom_info%label_molname)) offset = SIZE(atom_info%label_molname)
      CALL reallocate(atom_info%label_molname,1,offset+natom)
      CALL reallocate(atom_info%label_resid,1,offset+natom)
      CALL reallocate(atom_info%label_resname,1,offset+natom)
      CALL reallocate(atom_info%label_atmname,1,offset+natom)
      CALL reallocate(atom_info%atm_charge,1,offset+natom)
      CALL reallocate(atom_info%atm_mass,1,offset+natom)
      DO iatom=1,natom
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,itemp)
        CALL parser_get_object(parser,string)
        CALL parser_get_object(parser,itemp)
        atom_info%label_atmname(iatom) = namearray1(itemp)
        atom_info%label_resname(iatom) = "SOL"
        atom_info%label_molname(iatom) = "SOL"
        CALL parser_get_object(parser,atom_info%atm_mass(iatom))
        CALL parser_get_object(parser,atom_info%atm_charge(iatom))
        IF(itemp==1) iresid=iresid+1
        atom_info%label_resid(iatom)=iresid
        IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT SOLUTEATOM INFO HERE!!!!"
        CALL parser_get_next_line(parser,1,error=error)
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BONDH SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the BONDH section'    
    CALL parser_create(parser,file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(8))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,ntype)
      offset = 0
      IF(ASSOCIATED(conn_info%bond_a)) offset = SIZE(conn_info%bond_a)
      CALL reallocate(conn_info%bond_a,1,offset+ntype)
      CALL reallocate(conn_info%bond_b,1,offset+ntype)
      DO itype=1,ntype
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,conn_info%bond_a(offset+itype))
        CALL parser_get_object(parser,conn_info%bond_b(offset+itype))
        IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT BONDH INFO HERE!!!!"
      END DO
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BOND SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the BOND section'
    label = TRIM(avail_section(9))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,ntype)
      offset = 0
      IF(ASSOCIATED(conn_info%bond_a)) offset = SIZE(conn_info%bond_a)
      CALL reallocate(conn_info%bond_a,1,offset+ntype)
      CALL reallocate(conn_info%bond_b,1,offset+ntype)
      DO itype=1,ntype
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,conn_info%bond_a(offset+itype))
        CALL parser_get_object(parser,conn_info%bond_b(offset+itype))
        IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT BOND INFO HERE!!!!"
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BONDANGLEH SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the BONDANGLEH section'    
    CALL parser_create(parser,file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(11))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,ntype)
      offset = 0
      IF(ASSOCIATED(conn_info%theta_a)) offset = SIZE(conn_info%theta_a)
      CALL reallocate(conn_info%theta_a,1,offset+ntype)
      CALL reallocate(conn_info%theta_b,1,offset+ntype)
      CALL reallocate(conn_info%theta_c,1,offset+ntype)
      DO itype=1,ntype
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,conn_info%theta_a(offset+itype))
        CALL parser_get_object(parser,conn_info%theta_b(offset+itype))
        CALL parser_get_object(parser,conn_info%theta_c(offset+itype))
        IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT BONDANGLEH INFO HERE!!!!"
      END DO
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BONDANGLE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the BONDANGLE section'
    label = TRIM(avail_section(12))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,ntype)
      offset = SIZE(conn_info%theta_a)
      CALL reallocate(conn_info%theta_a,1,offset+ntype)
      CALL reallocate(conn_info%theta_b,1,offset+ntype)
      CALL reallocate(conn_info%theta_c,1,offset+ntype)
      DO itype=1,ntype
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,conn_info%theta_a(offset+itype))
        CALL parser_get_object(parser,conn_info%theta_b(offset+itype))
        CALL parser_get_object(parser,conn_info%theta_c(offset+itype))
        IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT BONDANGLE INFO HERE!!!!"
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! IMPDIHEDRALH SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the IMPDIHEDRALH section'    
    CALL parser_create(parser,file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(14))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,ntype)
      offset = 0
      IF(ASSOCIATED(conn_info%impr_a)) offset = SIZE(conn_info%impr_a)
      CALL reallocate(conn_info%impr_a,1,offset+ntype)
      CALL reallocate(conn_info%impr_b,1,offset+ntype)
      CALL reallocate(conn_info%impr_c,1,offset+ntype)
      CALL reallocate(conn_info%impr_d,1,offset+ntype)
      DO itype=1,ntype
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,conn_info%impr_a(offset+itype))
        CALL parser_get_object(parser,conn_info%impr_b(offset+itype))
        CALL parser_get_object(parser,conn_info%impr_c(offset+itype))
        CALL parser_get_object(parser,conn_info%impr_d(offset+itype))
        IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT IMPDIHEDRALH INFO HERE!!!!"
      END DO
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! IMPDIHEDRAL SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the IMPDIHEDRAL section'
    label = TRIM(avail_section(12))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,ntype)
      offset = 0
      IF(ASSOCIATED(conn_info%impr_a)) offset = SIZE(conn_info%impr_a)
      CALL reallocate(conn_info%impr_a,1,offset+ntype)
      CALL reallocate(conn_info%impr_b,1,offset+ntype)
      CALL reallocate(conn_info%impr_c,1,offset+ntype)
      CALL reallocate(conn_info%impr_d,1,offset+ntype)
      DO itype=1,ntype
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,conn_info%impr_a(offset+itype))
        CALL parser_get_object(parser,conn_info%impr_b(offset+itype))
        CALL parser_get_object(parser,conn_info%impr_c(offset+itype))
        CALL parser_get_object(parser,conn_info%impr_d(offset+itype))
        IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT IMPDIHEDRAL INFO HERE!!!!"
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DIHEDRALH SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the DIHEDRALH section'    
    CALL parser_create(parser,file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(14))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,ntype)
      offset = 0
      IF(ASSOCIATED(conn_info%phi_a)) offset = SIZE(conn_info%phi_a)
      CALL reallocate(conn_info%phi_a,1,offset+ntype)
      CALL reallocate(conn_info%phi_b,1,offset+ntype)
      CALL reallocate(conn_info%phi_c,1,offset+ntype)
      CALL reallocate(conn_info%phi_d,1,offset+ntype)
      DO itype=1,ntype
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,conn_info%phi_a(offset+itype))
        CALL parser_get_object(parser,conn_info%phi_b(offset+itype))
        CALL parser_get_object(parser,conn_info%phi_c(offset+itype))
        CALL parser_get_object(parser,conn_info%phi_d(offset+itype))
        IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT DIHEDRALH INFO HERE!!!!"
      END DO
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DIHEDRAL SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the DIHEDRAL section'
    label = TRIM(avail_section(12))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,ntype)
      offset = 0
      IF(ASSOCIATED(conn_info%phi_a)) offset = SIZE(conn_info%phi_a)
      CALL reallocate(conn_info%phi_a,1,offset+ntype)
      CALL reallocate(conn_info%phi_b,1,offset+ntype)
      CALL reallocate(conn_info%phi_c,1,offset+ntype)
      CALL reallocate(conn_info%phi_d,1,offset+ntype)
      DO itype=1,ntype
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,conn_info%phi_a(offset+itype))
        CALL parser_get_object(parser,conn_info%phi_b(offset+itype))
        CALL parser_get_object(parser,conn_info%phi_c(offset+itype))
        CALL parser_get_object(parser,conn_info%phi_d(offset+itype))
        IF(iw>0) WRITE(iw,*) "GTOP_INFO| PUT DIHEDRAL INFO HERE!!!!"
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! ikuo - Need to be debugged so that PARA_RES works...
    ! PARA_RES structure
    i = 0
    nbond_prev = 0
    IF(ASSOCIATED(conn_info%c_bond_a)) i = SIZE(conn_info%c_bond_a)
    nbond = SIZE(conn_info%bond_a)
    DO ibond=1+nbond_prev,nbond+nbond_prev
       iatom = conn_info%bond_a(ibond)
       jatom = conn_info%bond_b(ibond)
       IF(topology%para_res) THEN
          IF((atom_info%label_molname(iatom)/=atom_info%label_molname(jatom)).OR.&
               (atom_info%label_resid(iatom)/=atom_info%label_resid(jatom)).OR.&
               (atom_info%label_resname(iatom)/=atom_info%label_resname(jatom)))THEN
             IF(iw>0) WRITE(iw,'(T2,A)') "GTOP_INFO| PARA_RES, bond between molecules atom ",&
                  iatom,jatom
             i = i + 1
             CALL reallocate(conn_info%c_bond_a,1,i+1)
             CALL reallocate(conn_info%c_bond_b,1,i+1)
             conn_info%c_bond_a(i) = iatom
             conn_info%c_bond_b(i) = jatom
          END IF
       ELSE
          IF(atom_info%label_molname(iatom)/=atom_info%label_molname(jatom)) THEN
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
       END IF
    END DO


    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/GTOP_INFO",error=error)

  END SUBROUTINE read_topology_gromos

!!*****
!******************************************************************************
!!****** topology/write_topology_gromos [1.0] *
!!
!!   NAME
!!     write_topology_gromos
!!
!!   FUNCTION
!!     Write GROMOS topology file
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE write_topology_gromos (topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_topology_gromos', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_tag1, my_tag2, my_tag3, &
                                                record
    INTEGER                                  :: file_unit, handle, i, iw, j, &
                                                natom, nbond, nphi, ntheta
    LOGICAL                                  :: failure, ldum
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(connectivity_info_type), POINTER    :: conn_info
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/GTOP_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)

    atom_info => topology%atom_info
    conn_info => topology%conn_info
    
    record = "topology_dump.psf"

    IF(iw>0) WRITE(iw,*) "    Writing out PSF file ",TRIM(record)

    CALL open_file(file_name=record,&
         file_status="REPLACE",&
         file_action="WRITE",&
         unit_number=file_unit)

    WRITE(file_unit,*) "PSF"
    WRITE(file_unit,*) ""
    WRITE(file_unit,*) " 1 !NTITLE"
    WRITE(file_unit,*) "   CP2K generated DUMP of connectivity"
    WRITE(file_unit,*) ""
    
    WRITE(file_unit,*) " ",topology%natoms," !NATOM"
    DO i=1,topology%natoms
       my_tag1 = atom_info%label_molname(i)
       my_tag2 = atom_info%label_resname(i)
       my_tag3 = atom_info%label_atmname(i)
       ldum = qmmm_ff_precond_only_qm(my_tag1)
       ldum = qmmm_ff_precond_only_qm(my_tag2)
       ldum = qmmm_ff_precond_only_qm(my_tag3)
       WRITE(file_unit,'(I7,1X,A5,I7,1X,A5,A5,A5,F10.5,F10.5,I7)')&
            i,&
            TRIM(my_tag1),&
            atom_info%map_mol_num(i),&
            TRIM(my_tag2),&
            TRIM(my_tag3),&
            TRIM(my_tag3),&
            atom_info%atm_charge(i),&
            atom_info%atm_mass(i),0
    END DO
    WRITE(file_unit,*) ""
    
    WRITE(file_unit,*) " ",SIZE(conn_info%bond_a)," !NBOND"
    DO i=1,SIZE(conn_info%bond_a),4
       j=0
       DO WHILE ((j<4).AND.((i+j)<=SIZE(conn_info%bond_a)))
          WRITE(file_unit,'(I8,I8)',ADVANCE="NO") &
               conn_info%bond_a(i+j),conn_info%bond_b(i+j)
          j=j+1
       END DO
       WRITE(file_unit,*) ""
    END DO
    WRITE(file_unit,*) ""
    
    WRITE(file_unit,*) " ",SIZE(conn_info%theta_a)," !NTHETA"
    DO i=1,SIZE(conn_info%theta_a),3
       j=0
       DO WHILE ((j<3).AND.((i+j)<=SIZE(conn_info%theta_a)))
          WRITE(file_unit,'(I8,I8,I8)',ADVANCE="NO") &
               conn_info%theta_a(i+j),conn_info%theta_b(i+j),&
               conn_info%theta_c(i+j)
          j=j+1
       END DO
       WRITE(file_unit,*) ""
    END DO
    WRITE(file_unit,*) ""
    
    WRITE(file_unit,*) " ",SIZE(conn_info%phi_a)," !NPHI"
    DO i=1,SIZE(conn_info%phi_a),2
       j=0
       DO WHILE ((j<2).AND.((i+j)<=SIZE(conn_info%phi_a)))
          WRITE(file_unit,'(I8,I8,I8,I8)',ADVANCE="NO") &
               conn_info%phi_a(i+j),conn_info%phi_b(i+j),&
               conn_info%phi_c(i+j),conn_info%phi_d(i+j)
          j=j+1
       END DO
       WRITE(file_unit,*) ""
    END DO
    WRITE(file_unit,*) ""
    
    WRITE(file_unit,*) " ",SIZE(conn_info%impr_a)," !NIMPHI"
    DO i=1,SIZE(conn_info%impr_a),2
       j=0
       DO WHILE ((j<2).AND.((i+j)<=SIZE(conn_info%impr_a)))
          WRITE(file_unit,'(I8,I8,I8,I8)',ADVANCE="NO") &
               conn_info%impr_a(i+j),conn_info%impr_b(i+j),&
               conn_info%impr_c(i+j),conn_info%impr_d(i+j)
          j=j+1
       END DO
       WRITE(file_unit,*) ""
    END DO
    WRITE(file_unit,*) ""
    
    WRITE(file_unit,*) " 0 !NDON"
    WRITE(file_unit,*) ""
    WRITE(file_unit,*) " 0 !NACC"
    WRITE(file_unit,*) ""
    WRITE(file_unit,*) " 0 !NNB"
    WRITE(file_unit,*) ""

    CALL close_file(file_unit)

    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/GTOP_INFO",error=error)  
    CALL timestop(0.0_dp,handle)
    
  END SUBROUTINE write_topology_gromos















!******************************************************************************
!!****** topology_gromos/read_coordinate_g96 [1.0] *
!!
!!   NAME
!!     read_coordinate_g96
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *************************************************************************
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE read_coordinate_g96 (topology,globenv,subsys_section,error)
    TYPE(topology_parameters_type)           :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_g96', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock_atom = 10000

    INTEGER                                  :: handle, &
                                                iw, iatom, natom
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    CHARACTER(LEN=default_string_length), &
      DIMENSION(20)                          :: avail_section

    CHARACTER(LEN=default_string_length)     :: dummy_field, label,string
    LOGICAL                                  :: found

    NULLIFY(parser, logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/G96_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)

    atom_info => topology%atom_info


    IF (iw>0)  WRITE(iw,*) "    Reading in G96 file ",TRIM(topology%coord_file_name)
    CALL parser_create(parser,topology%coord_file_name,&
         para_env=globenv%para_env,error=error)  


    avail_section( 1) = "TITLE"
    avail_section( 2) = "TIMESTEP"
    avail_section( 3) = "POSITION"
    avail_section( 5) = "VELOCITY"
    avail_section( 6) = "BOX"




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TITLE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'G96_INFO| Parsing the TITLE section'    
    CALL parser_create(parser,topology%coord_file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(1))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      DO 
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,string,string_length=80)
        IF(string == TRIM("END")) EXIT
        IF(iw>0) WRITE(iw,*) "G96_INFO| ",TRIM(string)
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! POSITION SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'G96_INFO| Parsing the POSITION section'    
    CALL parser_create(parser,topology%coord_file_name,para_env=globenv%para_env)
    label = TRIM(avail_section(3))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.)
    IF(found) THEN
      CALL parser_get_next_line(parser,1,error=error)
      CALL parser_get_object(parser,natom)
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
      DO iatom=1,natom
        CALL parser_get_next_line(parser,1,error=error)
        CALL parser_get_object(parser,atom_info%label_resid(iatom))
        CALL parser_get_object(parser,atom_info%label_resname(iatom))
        CALL parser_get_object(parser,atom_info%label_atmname(iatom))
        CALL parser_get_object(parser,atom_info%r(1,iatom))
        CALL parser_get_object(parser,atom_info%r(2,iatom))
        CALL parser_get_object(parser,atom_info%r(3,iatom))
        atom_info%label_molname(iatom)=atom_info%label_resname(iatom)
        IF(iw>0) WRITE(iw,*) "G96_INFO| PUT POSITION INFO HERE!!!!"
      END DO
    END IF
    CALL parser_release(parser,error=error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    IF(.NOT.topology%para_res) atom_info%label_resid(:) = 1
    
    topology%natoms = natom
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/G96_INFO",error=error)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE read_coordinate_g96


END MODULE topology_gromos

!******************************************************************************
